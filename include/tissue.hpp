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


struct TCell {
  int64_t id;
  //int64_t x, y, z;
};

struct EpiCell {
  int64_t id;
  int64_t x, y, z;
  vector<TCell> tcells;
  double cytokines;
  double virus;

  string to_string() {
    ostringstream oss;
    oss << id << " (" << x << ", " << y << ", " << z << ") " << cytokines << " " << virus;
    return oss.str();
  }
};

class Tissue {
 private:
  static int block_size;
  static int64_t x_dim, y_dim, z_dim;
  
  using epi_cells_t = upcxx::dist_object<vector<EpiCell>>;
  epi_cells_t epi_cells;
  int64_t num_infected = 0;
  int64_t num_tcells = 0;

  // TODO: maintain a list of active cells (i.e. those that need to be updated, e.g because of cytokines, viruses, T-cells).
  // This will reduce computation by avoiding having to iterate through the whole space.
  
  int64_t map_3d_to_1d(int64_t x, int64_t y, int64_t z) {
    return x + y * Tissue::x_dim + z * Tissue::x_dim * Tissue::y_dim;
  }

  std::tuple<int64_t, int64_t, int64_t> map_1d_to_3d(int64_t i) {
    int64_t z = i / (Tissue::x_dim * Tissue::y_dim);
    i = i % (Tissue::x_dim * Tissue::y_dim);
    int64_t y = i / Tissue::x_dim;
    int64_t x = i % Tissue::x_dim;
    return {x, y, z};
  }

  intrank_t get_rank_for_cell(int64_t x, int64_t y, int64_t z) {
    int64_t id = map_3d_to_1d(x, y, z);
    int64_t block_i = id / Tissue::block_size;
    return block_i % rank_n();
  }
  
  int64_t gcd(int64_t x, int64_t y) {
    return (!x ? y : gcd(y % x, x));
  }

  static EpiCell *get_local_cell(epi_cells_t &epi_cells, int64_t id, int64_t x, int64_t y, int64_t z) {
    int64_t block_i = id / Tissue::block_size / rank_n();
    int64_t i = id % Tissue::block_size + block_i * Tissue::block_size;
    assert(i < epi_cells->size());
    EpiCell *epi_cell = &(*epi_cells)[i];
    assert(epi_cell->id == id);
    assert(epi_cell->x == x);
    assert(epi_cell->y == y);
    assert(epi_cell->z == z);
    return epi_cell;
  }    
    
  void set_infection(int64_t x, int64_t y, int64_t z, double concentration) {
    int64_t id = map_3d_to_1d(x, y, z);
    num_infected++;
    upcxx::rpc(get_rank_for_cell(x, y, z),
               [](epi_cells_t &epi_cells, int64_t id, int64_t x, int64_t y, int64_t z, double concentration) {
                 auto epi_cell = Tissue::get_local_cell(epi_cells, id, x, y, z);
                 epi_cell->virus = concentration;
               }, epi_cells, id, x, y, z, concentration).wait();
  }

  void add_tcell(int64_t x, int64_t y, int64_t z, int64_t tcell_id) {
    // FIXME: this should be an aggregating store
    int64_t id = map_3d_to_1d(x, y, z);
    num_tcells++;
    upcxx::rpc(get_rank_for_cell(x, y, z),
               [](epi_cells_t &epi_cells, int64_t id, int64_t x, int64_t y, int64_t z, int64_t tcell_id) {
                 auto epi_cell = Tissue::get_local_cell(epi_cells, id, x, y, z);
                 epi_cell->tcells.push_back({.id = tcell_id});
               }, epi_cells, id, x, y, z, tcell_id).wait();
  }
  
 public:
  Tissue() : epi_cells({}) {}

  ~Tissue() {}

  int64_t get_num_cells() {
    return Tissue::x_dim * Tissue::y_dim * Tissue::z_dim;
  }
  
  void construct(int64_t x_dim, int64_t y_dim, int64_t z_dim) {
    BarrierTimer timer(__FILEFUNC__, false, true);
    Tissue::x_dim = x_dim;
    Tissue::y_dim = y_dim;
    Tissue::z_dim = z_dim;
    int64_t num_cells = Tissue::x_dim * Tissue::y_dim * Tissue::z_dim;
    int64_t max_block_size = num_cells / rank_n();
    // we want the blocks to be cubes, so find a size that divides all three dimensions
    auto block_dim = gcd(gcd(Tissue::x_dim, Tissue::y_dim), Tissue::z_dim);
    Tissue::block_size = block_dim * block_dim * block_dim;
    // reduce block size until we can have at least one per process
    while (Tissue::block_size > max_block_size) {
      block_dim /= 2;
      Tissue::block_size = block_dim * block_dim * block_dim;
    }
    if (block_dim == 1)
      WARN("Using a block size of 1: this will result in a lot of communication. You should change the dimensions.");
    int64_t x_blocks = Tissue::x_dim / block_dim;
    int64_t y_blocks = Tissue::y_dim / block_dim;
    int64_t z_blocks = Tissue::z_dim / block_dim;
    int64_t num_blocks = num_cells / Tissue::block_size;
    SLOG("Dividing ", num_cells, " cells into ", num_blocks, " blocks of size ", Tissue::block_size, " (", block_dim, "^3)\n");
    int64_t blocks_per_rank = ceil((double)num_blocks / rank_n());
    epi_cells->reserve(blocks_per_rank * Tissue::block_size);
    auto mem_reqd = sizeof(EpiCell) * blocks_per_rank * Tissue::block_size;
    SLOG("Total initial memory required per process is ", get_size_str(mem_reqd), "\n");
    // FIXME: it may be more efficient (less communication) to have contiguous blocks
    // this is the quick & dirty approach
    for (int64_t i = 0; i < blocks_per_rank; i++) {
      int64_t start_id = (i * rank_n() + rank_me()) * Tissue::block_size;
      if (start_id >= num_cells) break;
      for (auto id = start_id; id < start_id + Tissue::block_size; id++) {
        auto [x, y, z] = map_1d_to_3d(id);
        epi_cells->push_back({.id = id, .x = x, .y = y, .z = z, .tcells = {}, .cytokines = 0, .virus = 0});
      }
    }
#ifdef DEBUG
    for (auto &epi_cell : (*epi_cells)) {
      DBG("cell ", epi_cell.to_string(), "\n");
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
      auto x = rnd_gen.get(0, Tissue::x_dim);
      auto y = rnd_gen.get(0, Tissue::y_dim);
      auto z = rnd_gen.get(0, Tissue::z_dim);
      DBG("infection: ", x, ", ", y, ", ", z, "\n");
      set_infection(x, y, z, 1.0);
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
      auto x = rnd_gen.get(0, Tissue::x_dim);
      auto y = rnd_gen.get(0, Tissue::y_dim);
      auto z = rnd_gen.get(0, Tissue::z_dim);
      DBG("tcell: ", x, ", ", y, ", ", z, "\n");
      add_tcell(x, y, z, tcell_id);
      upcxx::progress();
    }
    barrier();
    SLOG("Starting with ", get_num_tcells(), " infected cells\n");
  }
  
  int64_t get_num_tcells() {
    return upcxx::reduce_one(num_tcells, upcxx::op_fast_add, 0).wait();
  }

};

inline int64_t Tissue::x_dim, Tissue::y_dim, Tissue::z_dim;
inline int Tissue::block_size;

