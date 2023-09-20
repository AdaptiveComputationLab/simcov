#pragma once

#include <math.h>
#include <stdarg.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <upcxx/upcxx.hpp>

#include "upcxx_utils/flat_aggr_store.hpp"
#include "upcxx_utils/log.hpp"
#include "upcxx_utils/progress_bar.hpp"
#include "upcxx_utils/timers.hpp"
#include "utils.hpp"
#include "options.hpp"

extern shared_ptr<Options> _options;

using upcxx::rank_me;
using upcxx::rank_n;

using std::array;
using std::list;
using std::pair;
using std::shared_ptr;
using std::to_string;
using std::vector;

enum class ViewObject { VIRUS, TCELL_TISSUE, EPICELL, CHEMOKINE };

inline string view_object_str(ViewObject view_object) {
  switch (view_object) {
    case ViewObject::TCELL_TISSUE: return "tcelltissue";
    case ViewObject::VIRUS: return "virus";
    case ViewObject::EPICELL: return "epicell";
    case ViewObject::CHEMOKINE: return "chemokine";
    default: DIE("Unknown view object");
  }
  return "";
}

struct GridBlocks {
  int64_t block_size, num_x, num_y, num_z, size_x, size_y, size_z;
};

inline GridBlocks _grid_blocks;

struct GridCoords;
inline shared_ptr<GridCoords> _grid_size = nullptr;

struct GridCoords {
  int x, y, z;

  GridCoords() {}

  GridCoords(int64_t x, int64_t y, int64_t z)
      : x(x)
      , y(y)
      , z(z) {}

  // create a grid point from 1d
  GridCoords(int64_t i);

  // create a random grid point
  GridCoords(shared_ptr<Random> rnd_gen);

  void set_rnd(shared_ptr<Random> rnd_gen);

  bool operator==(const GridCoords &coords) {
    return x == coords.x && y == coords.y && z == coords.z;
  }

  bool operator!=(const GridCoords &coords) {
    return x != coords.x || y != coords.y || z != coords.z;
  }

  int64_t to_1d() const;

  static int64_t to_1d(int x, int y, int z);

  // convert linear coord system to block - needed to use with external lung model data
  static int64_t linear_to_block(int64_t i);

  string str() const {
    return "(" + to_string(x) + ", " + to_string(y) + ", " + to_string(z) + ")";
  }
};

struct TCell {
  string id;
  int binding_period = -1;
  int tissue_time_steps = -1;
  bool moved = true;
  int bindTarget = -1;
  int moveTarget = -1;
  int tie_break = -1;

  UPCXX_SERIALIZED_FIELDS(id, binding_period, tissue_time_steps, moved, moveTarget, bindTarget, tie_break);

  TCell(const string &id);

  TCell();
};

enum class EpiCellStatus { HEALTHY = 0, INCUBATING = 1, EXPRESSING = 2, APOPTOTIC = 3, DEAD = 4 };
const string EpiCellStatusStr[] = {"HEALTHY", "INCUBATING", "EXPRESSING", "APOPTOTIC", "DEAD"};
enum class EpiCellType { NONE, AIRWAY, ALVEOLI };

class EpiCell {
  int id;
  int incubation_time_steps = -1;
  int expressing_time_steps = -1;
  int apoptotic_time_steps = -1;

 public:
  EpiCellStatus status = EpiCellStatus::HEALTHY;
  EpiCellType type = EpiCellType::AIRWAY;
  bool infectable = true;
  bool boundTo = false;

  EpiCell(int id);

  string str();

  void infect();
  bool transition_to_expressing();
  bool apoptosis_death();
  bool infection_death();
  bool is_active();
  double get_binding_prob();
  bool was_expressing();
};

// 3D size 4*3+26*8+16+8+16=372
struct GridPoint {
  GridCoords coords;
  // empty space is nullptr
  EpiCell *epicell = nullptr;
  TCell *tcell = nullptr;
  // starts off empty and if calculated because this grid point becomes active, it is saved
  vector<int64_t> *neighbors = nullptr;
  float chemokine = 0, nb_chemokine = 0;
  float virions = 0, nb_virions = 0;

  string str() const;

  int tie_break = -1;

  bool is_active();
};

struct SampleData {
  double tcells = 0;
  bool has_epicell = false;
  EpiCellStatus epicell_status = EpiCellStatus::HEALTHY;
  float virions = 0;
  float chemokine = 0;
};

inline int64_t get_num_grid_points() {
  return (int64_t)_grid_size->x * (int64_t)_grid_size->y * (int64_t)_grid_size->z;
}

class Tissue {
 private:
  using grid_points_t = upcxx::dist_object<vector<GridPoint>>;
  grid_points_t grid_points;
  vector<GridPoint>::iterator grid_point_iter;

  // keeps track of all grid points that need to be updated
  using new_active_grid_points_t = upcxx::dist_object<HASH_TABLE<GridPoint *, bool>>;
  new_active_grid_points_t new_active_grid_points;

  HASH_TABLE<GridPoint *, bool> active_grid_points;
  HASH_TABLE<GridPoint *, bool>::iterator active_grid_point_iter;

  int64_t num_circulating_tcells;
  upcxx::dist_object<int64_t> tcells_generated;
  std::vector<EpiCellType> lung_cells;

  // this is static for ease of use in rpcs
  static GridPoint *get_local_grid_point(grid_points_t &grid_points, int64_t grid_i);

  int load_data_file(const string &fname, int num_grid_points, EpiCellType epicell_type);
  vector<int> get_model_dims(const string &fname);

 public:
  Tissue();

  ~Tissue() {}

  int64_t get_num_local_grid_points();

  intrank_t get_rank_for_grid_point(int64_t grid_i);

  vector<int64_t> *get_neighbors(GridCoords c);

  bool set_initial_infection(int64_t grid_i);

  void accumulate_chemokines(HASH_TABLE<int64_t, float> &chemokines_to_update,
                             IntermittentTimer &timer);

  void accumulate_virions(HASH_TABLE<int64_t, float> &virions_to_update, IntermittentTimer &timer);

  float get_chemokine(int64_t grid_i);

  bool tcells_in_neighborhood(GridPoint *grid_point);

  int64_t get_num_circulating_tcells();

  void change_num_circulating_tcells(int num);

  bool try_add_new_tissue_tcell(int64_t grid_i);

  bool try_add_tissue_tcell(int64_t grid_i, TCell &tcell);
  void prepare_tcell_move(int64_t grid_i, TCell &tcell);
  bool try_move_tissue_tcell(int64_t grid_i, TCell &tcell);

  bool try_bind_tcell(int64_t grid_i);

  GridPoint *get_first_local_grid_point();
  GridPoint *get_next_local_grid_point();

  GridPoint *get_first_active_grid_point();
  GridPoint *get_next_active_grid_point();

  void set_active(GridPoint *grid_point);
  void erase_active(GridPoint *grid_point);

  void add_new_actives(IntermittentTimer &timer);

  size_t get_num_actives();

  SampleData get_grid_point_sample_data(int64_t grid_i);

  // int64_t get_random_airway_epicell_location();

#ifdef DEBUG
  void check_actives(int time_step);
#endif
};
