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

using std::vector;
using std::array;
using std::list;
using std::pair;
using std::shared_ptr;
using std::to_string;

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
  int64_t x, y, z;

  GridCoords() {}

  GridCoords(int64_t x, int64_t y, int64_t z) : x(x), y(y), z(z) {}

  // create a grid point from 1d
  GridCoords(int64_t i);

  // create a random grid point
  GridCoords(shared_ptr<Random> rnd_gen);

  void from_1d_linear(int64_t i);

  void set_rnd(shared_ptr<Random> rnd_gen);

  bool operator==(const GridCoords &coords) {
    return x == coords.x && y == coords.y && z == coords.z;
  }

  bool operator!=(const GridCoords &coords) {
    return x != coords.x || y != coords.y || z != coords.z;
  }

  int64_t to_1d() const;
  int64_t to_1d_linear() const;

  string str() const {
    return "(" + to_string(x) + ", " + to_string(y) + ", " + to_string(z) + ")";
  }
};

struct TCell {
  string id;
  int binding_period = -1;
  int vascular_time_steps = -1;
  int tissue_time_steps = -1;
  bool moved = true;

  UPCXX_SERIALIZED_FIELDS(id, binding_period, moved);

  TCell(const string &id);

  TCell();
};

enum class EpiCellStatus { HEALTHY, INCUBATING, EXPRESSING, APOPTOTIC, DEAD };
const string EpiCellStatusStr[] = {"HEALTHY", "INCUBATING", "EXPRESSING", "APOPTOTIC", "DEAD" };

class EpiCell {
  int id;
  int infection_time_step = -1;
  int incubation_time_steps = -1;
  int expressing_time_steps = -1;
  int apoptotic_time_steps = -1;
  bool is_expressing = false;

 public:
  EpiCellStatus status = EpiCellStatus::HEALTHY;
  bool infectable = true;

  EpiCell(int id);

  string str();

  void infect(int time_stamp);
  bool transition_to_expressing();
  bool apoptosis_death();
  bool infection_death();
  bool is_active();
  double get_binding_prob(int time_step);
  bool was_expressing();
};

struct GridPoint {
  int64_t id;
  GridCoords coords;
  // vector for connectivity
  vector<GridCoords> neighbors;
  // empty space is nullptr
  EpiCell *epicell = nullptr;
  TCell *tcell = nullptr;
  int virions = 0, nb_virions = 0;
  double chemokine = 0, nb_chemokine = 0;

  string str() const;

  bool is_active();
};

struct SampleData {
  int64_t id = -1;
  bool has_tcell = false;
  bool has_epicell = false;
  EpiCellStatus epicell_status = EpiCellStatus::HEALTHY;
  int virions = 0;
  double chemokine = 0;
};

inline int64_t get_num_grid_points() { return _grid_size->x * _grid_size->y * _grid_size->z; }

class Tissue {
 private:

  using grid_points_t = upcxx::dist_object<vector<GridPoint>>;
  grid_points_t grid_points;
  vector<GridPoint>::iterator grid_point_iter;

  // keeps track of all grid points that need to be updated
  using new_active_grid_points_t = upcxx::dist_object<HASH_TABLE<GridPoint*, bool>>;
  new_active_grid_points_t new_active_grid_points;

  HASH_TABLE<GridPoint*, bool> active_grid_points;
  HASH_TABLE<GridPoint*, bool>::iterator active_grid_point_iter;

  upcxx::dist_object<list<TCell>> circulating_tcells;

  // this is static for ease of use in rpcs
  static GridPoint *get_local_grid_point(grid_points_t &grid_points, const GridCoords &coords);

  SampleData get_grid_point_sample_data(const GridCoords &coords);

  vector<GridCoords> get_neighbors(GridCoords c);

 public:
  int64_t tcells_generated = 0;

  Tissue();

  ~Tissue() {}

  int64_t get_num_local_grid_points();

  intrank_t get_rank_for_grid_point(const GridCoords &coords);

  bool set_initial_infection(GridCoords coords);

  void accumulate_chemokines(HASH_TABLE<int64_t, double> &chemokines_to_update,
                             IntermittentTimer &timer);

  void accumulate_virions(HASH_TABLE<int64_t, int> &virions_to_update, IntermittentTimer &timer);

  double get_chemokine(GridCoords coords);

  list<TCell> *get_circulating_tcells();

  bool tcells_in_neighborhood(GridPoint *grid_point);

  void add_circulating_tcell(GridCoords coords, TCell tcell);

  bool try_add_tissue_tcell(GridCoords coords, TCell tcell, bool extravasate);

  EpiCellStatus try_bind_tcell(GridCoords coords, int time_step);

  GridPoint *get_first_local_grid_point();
  GridPoint *get_next_local_grid_point();

  GridPoint *get_first_active_grid_point();
  GridPoint *get_next_active_grid_point();

  void set_active(GridPoint *grid_point);
  void erase_active(GridPoint *grid_point);

  void add_new_actives(IntermittentTimer &timer);

  size_t get_num_actives();

  void get_samples(vector<SampleData> &samples, int64_t &start_id);

#ifdef DEBUG
  void check_actives(int time_step);
#endif

};

