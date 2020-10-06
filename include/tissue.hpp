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

using std::unordered_map;
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

struct GridCoords {
  int64_t x, y, z;

  GridCoords() {}

  GridCoords(int64_t x, int64_t y, int64_t z) : x(x), y(y), z(z) {}

  // create a grid point from 1d
  GridCoords(int64_t i, const GridCoords &grid_size);

  // create a random grid point
  GridCoords(shared_ptr<Random> rnd_gen, const GridCoords &grid_size);

  void set_rnd(shared_ptr<Random> rnd_gen, const GridCoords &grid_size);

  bool operator==(const GridCoords &coords) {
    return x == coords.x && y == coords.y && z == coords.z;
  }

  bool operator!=(const GridCoords &coords) {
    return x != coords.x || y != coords.y || z != coords.z;
  }

  int64_t to_1d(const GridCoords &grid_size) const {
    if (x >= grid_size.x || y >= grid_size.y || z >= grid_size.z)
      DIE("Grid point is out of range: ", str(), " max size ", grid_size.str());
    return x + y * grid_size.x + z * grid_size.x * grid_size.y;
  }

  string str() const {
    return "(" + to_string(x) + ", " + to_string(y) + ", " + to_string(z) + ")";
  }
};

struct TCell {
  string id;
  int vascular_period = -1;
  int tissue_period = -1;
  int binding_period = -1;
  bool moved = true;

  UPCXX_SERIALIZED_FIELDS(id, vascular_period, tissue_period, binding_period, moved);

  TCell(const string &id, GridCoords &coords);

  TCell() {}
};

enum class EpiCellStatus { HEALTHY, INCUBATING, EXPRESSING, APOPTOTIC, DEAD };
const string EpiCellStatusStr[] = {"HEALTHY", "INCUBATING", "EXPRESSING", "APOPTOTIC", "DEAD" };

class EpiCell {
  int64_t id;
  int initial_incubation_period = -1;
  int incubation_period = -1;
  int apoptosis_period = -1;
  int infection_period = -1;

 public:
  EpiCellStatus status = EpiCellStatus::HEALTHY;
  bool infectable = false;

  EpiCell(int id) : id(id) {};

  string str();

  void infect();
  bool transition_to_expressing();
  bool apoptosis_death();
  bool infection_death();
  bool is_active();
  double get_binding_prob();
};

struct GridPoint {
  int64_t id;
  GridCoords coords;
  // vector for connectivity
  vector<GridCoords> neighbors;
  // empty space is nullptr
  EpiCell *epicell = nullptr;
  TCell *tcell = nullptr;
  double virus = 0, nb_virus = 0;
  double chemokine = 0, nb_chemokine = 0;

  string str() const;

  bool is_active();
};


using grid_to_conc_map_t = unordered_map<int64_t, array<double, 2>>;

class Tissue {
 private:
  // these are static so they don't have to be passed in RPCs
  static int block_size;

  using grid_points_t = upcxx::dist_object<vector<GridPoint>>;
  grid_points_t grid_points;
  vector<GridPoint>::iterator grid_point_iter;

  // keeps track of all grid points that need to be updated
  using new_active_grid_points_t = upcxx::dist_object<unordered_map<GridPoint*, bool>>;
  new_active_grid_points_t new_active_grid_points;

  unordered_map<GridPoint*, bool> active_grid_points;
  unordered_map<GridPoint*, bool>::iterator active_grid_point_iter;

  upcxx::dist_object<list<TCell>> circulating_tcells;

  static GridPoint *get_local_grid_point(grid_points_t &grid_points, const GridCoords &coords);

  vector<GridCoords> get_neighbors(GridCoords c, GridCoords grid_size);

 public:
  static GridCoords grid_size;

  int64_t tcells_generated = 0;

  Tissue() : grid_points({}), new_active_grid_points({}), circulating_tcells({}) {};

  ~Tissue() {}

  static int64_t get_num_grid_points();

  int64_t get_num_local_grid_points();

  intrank_t get_rank_for_grid_point(const GridCoords &coords);

  void set_initial_infection(GridCoords coords);

  void accumulate_concentrations(grid_to_conc_map_t &concs_to_update, IntermittentTimer &timer);

  double get_chemokine(GridCoords coords);

  list<TCell> *get_circulating_tcells();

  bool tcells_in_neighborhood(GridPoint *grid_point);

  void add_circulating_tcell(GridCoords coords, TCell tcell);

  bool try_add_tissue_tcell(GridCoords coords, TCell tcell, bool extravasate);

  void construct(GridCoords grid_size);

  pair<size_t, size_t> dump_blocks(const string &fname, const string &header_str,
                                   ViewObject view_object);

  GridPoint *get_first_local_grid_point();
  GridPoint *get_next_local_grid_point();

  GridPoint *get_first_active_grid_point();
  GridPoint *get_next_active_grid_point();

  void set_active(GridPoint *grid_point);
  void erase_active(GridPoint *grid_point);

  void add_new_actives(IntermittentTimer &timer);

  size_t get_num_actives();

#ifdef DEBUG
  void check_actives(int time_step);
#endif

};

inline GridCoords Tissue::grid_size;
inline int Tissue::block_size;
