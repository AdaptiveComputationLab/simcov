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


#define SANITY_CHECK_INDEXES


using upcxx::rank_me;
using upcxx::rank_n;


struct TCell {
  int64_t id;
  int64_t x, y, z;
};

struct EpiCell {
  int64_t id;
  int64_t x, y, z;
  vector<TCell> t_cells;
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
  using epi_cells_t = upcxx::dist_object<vector<EpiCell>>;
  epi_cells_t epi_cells;
  int block_dim;
  int block_size;
  int64_t x_dim, y_dim, z_dim;

  int64_t map_3d_to_1d(int64_t x, int64_t y, int64_t z) {
    return x + y * x_dim + z * x_dim * y_dim;
  }

  std::tuple<int64_t, int64_t, int64_t> map_1d_to_3d(int64_t i) {
    int64_t z = i / (x_dim * y_dim);
    i = i % (x_dim * y_dim);
    int64_t y = i / x_dim;
    int64_t x = i % x_dim;
    return {x, y, z};
  }

  intrank_t get_rank_for_cell(int64_t x, int64_t y, int64_t z) {
    int64_t id = map_3d_to_1d(x, y, z);
    int64_t block_i = id / block_size;
    return block_i % rank_n();
  }
  
  int64_t gcd(int64_t x, int64_t y) {
    return (!x ? y : gcd(y % x, x));
  }

#ifdef SANITY_CHECK_INDEXES
  void sanity_check_cell(int64_t id, int64_t x, int64_t y, int64_t z) {
    upcxx::rpc(get_rank_for_cell(x, y, z),
               [](epi_cells_t &epi_cells, int64_t id, int64_t x, int64_t y, int64_t z, int block_size) {
                 // find the entry in our local epi_cells array and compare id, x, y and z
                 int64_t block_i = id / block_size / rank_n();
                 int64_t i = id % block_size + block_i * block_size;
                 EpiCell *epi_cell = &(*epi_cells)[i];
                 if (epi_cell->id != id) WARN("ID not equal ", id, " != ", epi_cell->id, "\n");
                 if (epi_cell->x != x) WARN("X not equal ", x, " != ", epi_cell->x, "\n");
                 if (epi_cell->y != y) WARN("Y not equal ", y, " != ", epi_cell->y, "\n");
                 if (epi_cell->z != z) WARN("Z not equal ", z, " != ", epi_cell->z, "\n");
               }, epi_cells, id, x, y, z, block_size).wait();
  }
#endif
  
 public:
  Tissue() : epi_cells({}) {}

  ~Tissue() {}

  void construct(int64_t x_dim, int64_t y_dim, int64_t z_dim) {
    BarrierTimer timer(__FILEFUNC__, false, true);
    this->x_dim = x_dim;
    this->y_dim = y_dim;
    this->z_dim = z_dim;
    int64_t num_cells = x_dim * y_dim * z_dim;
    int64_t max_block_size = num_cells / rank_n();
    // we want the blocks to be cubes, so find a size that divides all three dimensions
    block_dim = gcd(gcd(x_dim, y_dim), z_dim);
    block_size = block_dim * block_dim * block_dim;
    // reduce block size until we can have at least one per process
    while (block_size > max_block_size) {
      block_dim /= 2;
      block_size = block_dim * block_dim * block_dim;
    }
    int64_t x_blocks = x_dim / block_dim;
    int64_t y_blocks = y_dim / block_dim;
    int64_t z_blocks = z_dim / block_dim;
    int64_t num_blocks = num_cells / block_size;
    SLOG("Dividing ", num_cells, " cells into ", num_blocks, " blocks of size ", block_size, " (", block_dim, "^3)\n");
    int64_t blocks_per_rank = ceil((double)num_blocks / rank_n());
    epi_cells->reserve(blocks_per_rank * block_size);
    auto mem_reqd = sizeof(EpiCell) * blocks_per_rank * block_size;
    SLOG("Total initial memory required per process is ", get_size_str(mem_reqd), "\n");
    // FIXME: it may be more efficient (less communication) to have contiguous blocks
    // this is the quick & dirty approach
    for (int64_t i = 0; i < blocks_per_rank; i++) {
      int64_t start_id = (i * rank_n() + rank_me()) * block_size;
      if (start_id >= num_cells) break;
      for (auto id = start_id; id < start_id + block_size; id++) {
        auto [x, y, z] = map_1d_to_3d(id);
        epi_cells->push_back({.id = id, .x = x, .y = y, .z = z, .t_cells = {}, .cytokines = 0, .virus = 0});
      }
    }
#ifdef DEBUG
    for (auto &epi_cell : (*epi_cells)) {
      DBG("cell ", rank_me(), " ", epi_cell.to_string(), "\n");
    }
#endif
#ifdef SANITY_CHECK_INDEXES
    if (!rank_me()) {
      for (int xi = 0; xi < x_dim; xi++) {
        for (int yi = 0; yi < y_dim; yi++) {
          for (int zi = 0; zi < z_dim; zi++) {
            int64_t id = map_3d_to_1d(xi, yi, zi);
            sanity_check_cell(id, xi, yi, zi);
          }
        }
      }
    }
#endif
    barrier();
  }
  
  void infect(int num_infections) {
    BarrierTimer timer(__FILEFUNC__, false, true);
  }
};
