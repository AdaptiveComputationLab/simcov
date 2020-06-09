#pragma once

#include <stdarg.h>
#include <math.h>
#include <iostream>
#include <map>
#include <fstream>
#include <numeric>
#include <algorithm>
#include <upcxx/upcxx.hpp>

#include "upcxx_utils/log.hpp"
#include "upcxx_utils/progress_bar.hpp"
#include "upcxx_utils/flat_aggr_store.hpp"
#include "upcxx_utils/timers.hpp"

#include "utils.hpp"


using upcxx::rank_me;
using upcxx::rank_n;

enum class TCellStatus { Activated, Circulating, Emigrating, Chemotaxing, Dead };

struct TCell {
  int64_t id;
  TCellStatus status = TCellStatus::Activated;
};

enum class EpiCellStatus { Healthy, Incubating, Secreting, Apoptotic, Dead };

struct EpiCell {
  int64_t id;
  EpiCellStatus status = EpiCellStatus::Healthy;

  string str() {
    return std::to_string(id);
  }
};

struct GridCoords {
  int64_t x, y, z;

  GridCoords() {}
  
  GridCoords(int64_t x, int64_t y, int64_t z) : x(x), y(y), z(z) {}

  // create a grid point from 1d
  GridCoords(int64_t i, const GridCoords &grid_size) {
    z = i / (grid_size.x * grid_size.y);
    i = i % (grid_size.x * grid_size.y);
    y = i / grid_size.x;
    x = i % grid_size.x;
  }

  // create a random grid point
  GridCoords(Random rnd_gen, const GridCoords &grid_size) {
    x = rnd_gen.get(0, grid_size.x);
    y = rnd_gen.get(0, grid_size.y);
    z = rnd_gen.get(0, grid_size.z);
  }
  
  int64_t to_1d(const GridCoords &grid_size) const {
    return x + y * grid_size.x + z * grid_size.x * grid_size.y;
  }

  string str() {
    return "(" + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) + ")";
  }
};

struct GridPoint {
  int64_t id;
  GridCoords coords;
  // empty space is nullptr
  EpiCell *epicell = nullptr;
  vector<TCell> tcells = {};
  double cytokines = 0.0;
  double virus = 0.0;

  string str() {
    ostringstream oss;
    oss << id << " " << coords.str() << " " << epicell->str() << " " << cytokines << " " << virus;
    return oss.str();
  }
};

class Tissue {
 private:
  // these are static so they don't have to be passed in RPCs
  static int block_size;
  static GridCoords grid_size;
  
  using grid_points_t = upcxx::dist_object<vector<GridPoint>>;
  grid_points_t grid_points;
  int64_t num_infected = 0;
  int64_t num_tcells = 0;

  // TODO: maintain a list of active cells (i.e. those that need to be updated, e.g because of cytokines, viruses, T-cells).
  // This will reduce computation by avoiding having to iterate through the whole space.
  
  intrank_t get_rank_for_grid_point(const GridCoords &coords) {
    int64_t id = coords.to_1d(Tissue::grid_size);
    int64_t block_i = id / Tissue::block_size;
    return block_i % rank_n();
  }
  
  static GridPoint *get_local_grid_point(grid_points_t &grid_points, int64_t id, const GridCoords &coords) {
    int64_t block_i = id / Tissue::block_size / rank_n();
    int64_t i = id % Tissue::block_size + block_i * Tissue::block_size;
    assert(i < grid_points->size());
    GridPoint *grid_point = &(*grid_points)[i];
    assert(grid_point->id == id);
    assert(grid_point->coords.x == coords.x);
    assert(grid_point->coords.y == coords.y);
    assert(grid_point->coords.z == coords.z);
    return grid_point;
  }    
    
  void set_infection(const GridCoords &coords, double concentration) {
    int64_t id = coords.to_1d(Tissue::grid_size);
    num_infected++;
    upcxx::rpc(get_rank_for_grid_point(coords),
               [](grid_points_t &grid_points, int64_t id, GridCoords coords, double concentration) {
                 auto grid_point = Tissue::get_local_grid_point(grid_points, id, coords);
                 grid_point->virus = concentration;
               }, grid_points, id, coords, concentration).wait();
  }

  void add_tcell(const GridCoords &coords, int64_t tcell_id) {
    // FIXME: this should be an aggregating store
    int64_t id = coords.to_1d(Tissue::grid_size);
    num_tcells++;
    upcxx::rpc(get_rank_for_grid_point(coords),
               [](grid_points_t &grid_points, int64_t id, GridCoords coords, int64_t tcell_id) {
                 auto grid_point = Tissue::get_local_grid_point(grid_points, id, coords);
                 grid_point->tcells.push_back({.id = tcell_id, .status = TCellStatus::Activated});
               }, grid_points, id, coords, tcell_id).wait();
  }
  
 public:
  Tissue() : grid_points({}) {}

  ~Tissue() {}

  int64_t get_num_grid_points() {
    return Tissue::grid_size.x * Tissue::grid_size.y * Tissue::grid_size.z;
  }
  
  void construct(GridCoords grid_size) {
    BarrierTimer timer(__FILEFUNC__, false, true);
    Tissue::grid_size = grid_size;
    int64_t num_grid_points = get_num_grid_points();
    int64_t max_block_size = num_grid_points / rank_n();
    // we want the blocks to be cubes, so find a size that divides all three dimensions
    auto block_dim = gcd(gcd(Tissue::grid_size.x, Tissue::grid_size.y), Tissue::grid_size.z);
    Tissue::block_size = block_dim * block_dim * block_dim;
    // reduce block size until we can have at least one per process
    while (Tissue::block_size > max_block_size) {
      block_dim /= 2;
      Tissue::block_size = block_dim * block_dim * block_dim;
    }
    if (block_dim == 1)
      WARN("Using a block size of 1: this will result in a lot of communication. You should change the dimensions.");
    int64_t x_blocks = Tissue::grid_size.x / block_dim;
    int64_t y_blocks = Tissue::grid_size.y / block_dim;
    int64_t z_blocks = Tissue::grid_size.z / block_dim;
    int64_t num_blocks = num_grid_points / Tissue::block_size;
    SLOG("Dividing ", num_grid_points, " grid points into ", num_blocks, " blocks of size ", Tissue::block_size,
         " (", block_dim, "^3)\n");
    int64_t blocks_per_rank = ceil((double)num_blocks / rank_n());
    grid_points->reserve(blocks_per_rank * Tissue::block_size);
    auto mem_reqd = sizeof(GridPoint) * blocks_per_rank * Tissue::block_size;
    SLOG("Total initial memory required per process is ", get_size_str(mem_reqd), "\n");
    // FIXME: it may be more efficient (less communication) to have contiguous blocks
    // this is the quick & dirty approach
    for (int64_t i = 0; i < blocks_per_rank; i++) {
      int64_t start_id = (i * rank_n() + rank_me()) * Tissue::block_size;
      if (start_id >= num_grid_points) break;
      for (auto id = start_id; id < start_id + Tissue::block_size; id++) {
        GridCoords coords(id, Tissue::grid_size);
        // FIXME: this is an epicall on every grid point.
        // They should be placed according to the underlying lung structure (gaps, etc)
        grid_points->push_back({.id = id, .coords = coords, .epicell = new EpiCell({.id = id, .status = EpiCellStatus::Healthy}),
                                .tcells = {}, .cytokines = 0, .virus = 0});
      }
    }
#ifdef DEBUG
    for (auto &grid_point : (*grid_points)) {
      DBG("grid point ", grid_point.str(), "\n");
    }
#endif
    barrier();
  }
  
  void infect(int num_infections) {
    BarrierTimer timer(__FILEFUNC__, false, true);
    Random rnd_gen;
    // each rank generates a block of infected cells
    // for large dimensions, the chance of repeat sampling for a few points is very small
    int local_num_infections = ceil((double)num_infections / rank_n());
    int64_t min_i = rank_me() * local_num_infections;
    int64_t max_i = (rank_me() + 1) * local_num_infections;
    if (max_i > num_infections) max_i = num_infections;
    if (min_i > num_infections) min_i = num_infections;
    barrier();
    for (int i = min_i; i < max_i; i++) {
      GridCoords coords(rnd_gen, Tissue::grid_size);
      DBG("infection: ", coords.str() + "\n");
      set_infection(coords, 1.0);
      upcxx::progress();
    }
    barrier();
    SLOG("Starting with ", get_num_infected(), " infected cells\n");
  }

  int64_t get_num_infected() {
    return upcxx::reduce_one(num_infected, upcxx::op_fast_add, 0).wait();
  }

  void generate_tcells(int num_tcells) {
    BarrierTimer timer(__FILEFUNC__, false, true);
    Random rnd_gen;
    // each rank generates a block of tcells
    int local_num_tcells = ceil((double)num_tcells / rank_n());
    int64_t min_tcell_id = rank_me() * local_num_tcells;
    int64_t max_tcell_id = (rank_me() + 1) * local_num_tcells;
    if (max_tcell_id > num_tcells) max_tcell_id = num_tcells;
    if (min_tcell_id > num_tcells) min_tcell_id = num_tcells;
    for (int64_t tcell_id = min_tcell_id; tcell_id < max_tcell_id; tcell_id++) {
      GridCoords coords(rnd_gen, Tissue::grid_size);
      DBG("tcell: ", coords.str(), "\n");
      add_tcell(coords, tcell_id);
      upcxx::progress();
    }
    barrier();
    SLOG("Starting with ", get_num_tcells(), " infected cells\n");
  }
  
  int64_t get_num_tcells() {
    return upcxx::reduce_one(num_tcells, upcxx::op_fast_add, 0).wait();
  }

};

inline GridCoords Tissue::grid_size;
inline int Tissue::block_size;

