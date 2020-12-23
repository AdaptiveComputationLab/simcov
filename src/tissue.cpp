#include "tissue.hpp"

using namespace std;
using upcxx::make_view;
using upcxx::progress;
using upcxx::rpc;
using upcxx::view;

GridCoords::GridCoords(int64_t i) {
#ifdef BLOCK_PARTITION
  int64_t blocknum = i / _grid_blocks.block_size;
  int64_t block_z = blocknum / (_grid_blocks.num_x * _grid_blocks.num_y);
  blocknum -= block_z * _grid_blocks.num_x * _grid_blocks.num_y;
  int64_t block_y = blocknum / _grid_blocks.num_x;
  int64_t block_x = blocknum % _grid_blocks.num_x;
  int64_t in_block_id = i % _grid_blocks.block_size;
  block_x *= _grid_blocks.size_x;
  block_y *= _grid_blocks.size_y;
  block_z *= _grid_blocks.size_z;
  int64_t dz = in_block_id / (_grid_blocks.size_x * _grid_blocks.size_y);
  in_block_id -= (dz * _grid_blocks.size_x * _grid_blocks.size_y);
  int64_t dy = in_block_id / _grid_blocks.size_x;
  int64_t dx = in_block_id % _grid_blocks.size_x;
  x = block_x + dx;
  y = block_y + dy;
  z = block_z + dz;
#else
  from_1d_linear(i);
#endif
}

GridCoords::GridCoords(shared_ptr<Random> rnd_gen) {
  x = rnd_gen->get(0, _grid_size->x);
  y = rnd_gen->get(0, _grid_size->y);
  z = rnd_gen->get(0, _grid_size->z);
}

void GridCoords::from_1d_linear(int64_t i) {
  z = i / (_grid_size->x * _grid_size->y);
  i = i % (_grid_size->x * _grid_size->y);
  y = i / _grid_size->x;
  x = i % _grid_size->x;
}

int64_t GridCoords::to_1d() const {
  if (x >= _grid_size->x || y >= _grid_size->y || z >= _grid_size->z)
    DIE("Grid point is out of range: ", str(), " max size ", _grid_size->str());
#ifdef BLOCK_PARTITION
  int64_t block_x = x / _grid_blocks.size_x;
  int64_t block_y = y / _grid_blocks.size_y;
  int64_t block_z = z / _grid_blocks.size_z;
  int64_t block_id =
      block_x + block_y * _grid_blocks.num_x + block_z * _grid_blocks.num_x * _grid_blocks.num_y;
  int64_t in_block_x = x % _grid_blocks.size_x;
  int64_t in_block_y = y % _grid_blocks.size_y;
  int64_t in_block_z = z % _grid_blocks.size_z;
  int64_t in_block_id = in_block_x + in_block_y * _grid_blocks.size_x +
                        in_block_z * _grid_blocks.size_x * _grid_blocks.size_y;
  return in_block_id + block_id * _grid_blocks.block_size;
#else
  return to_1d_linear();
#endif
}

int64_t GridCoords::to_1d_linear() const {
  if (x >= _grid_size->x || y >= _grid_size->y || z >= _grid_size->z)
    DIE("Grid point is out of range: ", str(), " max size ", _grid_size->str());
  return x + y * _grid_size->x + z * _grid_size->x * _grid_size->y;
}

void GridCoords::set_rnd(shared_ptr<Random> rnd_gen) {
  x = rnd_gen->get(0, _grid_size->x);
  y = rnd_gen->get(0, _grid_size->y);
  z = rnd_gen->get(0, _grid_size->z);
}

TCell::TCell(const string &id)
    : id(id) {
  tissue_time_steps = _rnd_gen->get_poisson(_options->tcell_tissue_period);
  DBG("init tcell ", id, " ", tissue_time_steps, "\n");
}

TCell::TCell() { tissue_time_steps = _rnd_gen->get_poisson(_options->tcell_tissue_period); }

EpiCell::EpiCell(int id)
    : id(id) {
  incubation_time_steps = _rnd_gen->get_poisson(_options->incubation_period);
  expressing_time_steps = _rnd_gen->get_poisson(_options->expressing_period);
  apoptotic_time_steps = _rnd_gen->get_poisson(_options->apoptosis_period);
  DBG("init epicell ", str(), "\n");
}

string EpiCell::str() {
  ostringstream oss;
  oss << id << " " << EpiCellStatusStr[(int)status] << " " << incubation_time_steps << " "
      << expressing_time_steps << " " << apoptotic_time_steps;
  return oss.str();
}

void EpiCell::infect() {
  assert(status == EpiCellStatus::HEALTHY);
  assert(infectable);
  status = EpiCellStatus::INCUBATING;
}

bool EpiCell::transition_to_expressing() {
  assert(status == EpiCellStatus::INCUBATING);
  incubation_time_steps--;
  if (incubation_time_steps > 0) return false;
  status = EpiCellStatus::EXPRESSING;
  return true;
}

bool EpiCell::was_expressing() {
  // this is used to determine if the epicell was expressing before apoptosis was induced
  assert(status == EpiCellStatus::APOPTOTIC);
  return (incubation_time_steps == 0);
}

bool EpiCell::apoptosis_death() {
  assert(status == EpiCellStatus::APOPTOTIC);
  apoptotic_time_steps--;
  if (apoptotic_time_steps > 0) return false;
  status = EpiCellStatus::DEAD;
  return true;
}

bool EpiCell::infection_death() {
  expressing_time_steps--;
  if (expressing_time_steps > 0) return false;
  status = EpiCellStatus::DEAD;
  return true;
}

bool EpiCell::is_active() {
  return (status != EpiCellStatus::HEALTHY && status != EpiCellStatus::DEAD);
}

double EpiCell::get_binding_prob() {
  // binding prob is linearly scaled from 0 to 1 for incubating cells over the course of the
  // incubation period, but is always 1 for expressing cells
  if (status == EpiCellStatus::EXPRESSING || status == EpiCellStatus::APOPTOTIC)
    return _options->max_binding_prob;
  double scaling = 1.0 - (double)incubation_time_steps / _options->incubation_period;
  if (scaling < 0) scaling = 0;
  double prob = _options->max_binding_prob * scaling;
  return min(prob, _options->max_binding_prob);
}

string GridPoint::str() const {
  ostringstream oss;
  oss << "id " << id << ", xyz " << coords.str() << ", epi " << (epicell ? epicell->str() : "none")
      << ", v " << virions << ", c " << chemokine;
  return oss.str();
}

bool GridPoint::is_active() {
  // it could be incubating but without anything else set
  return ((epicell && epicell->is_active()) || virions || chemokine > 0 || tcell);
}

static int get_cube_block_dim(int64_t num_grid_points) {
  int block_dim = 1;
  int min_dim = min(min(_grid_size->x, _grid_size->y), _grid_size->z);
  for (int i = 1; i < min_dim; i++) {
    // only allow dims that divide each main dimension perfectly
    if (remainder(_grid_size->x, i) || remainder(_grid_size->y, i) || remainder(_grid_size->z, i)) {
      DBG("dim ", i, " does not divide all main dimensions cleanly\n");
      continue;
    }
    size_t cube = (size_t)pow((double)i, 3.0);
    size_t num_cubes = num_grid_points / cube;
    DBG("cube size ", cube, " num cubes ", num_cubes, "\n");
    if (num_cubes < rank_n() * _options->min_blocks_per_proc) {
      DBG("not enough cubes ", num_cubes, " < ", rank_n() * _options->min_blocks_per_proc, "\n");
      break;
    }
    // there is a remainder - this is not a perfect division
    if (remainder(num_grid_points, cube)) {
      DBG("there is a remainder - don't use\n");
      continue;
    } else {
      DBG("selected dim ", i, "\n");
      block_dim = i;
    }
  }
  return block_dim;
}

static int get_square_block_dim(int64_t num_grid_points) {
  int block_dim = 1;
  int min_dim = min(_grid_size->x, _grid_size->y);
  for (int i = 1; i < min_dim; i++) {
    // only allow dims that divide each main dimension perfectly
    if (remainder(_grid_size->x, i) || remainder(_grid_size->y, i)) {
      DBG("dim ", i, " does not divide all main dimensions cleanly\n");
      continue;
    }
    size_t square = (size_t)pow((double)i, 2.0);
    size_t num_squares = num_grid_points / square;
    DBG("square size ", square, " num squares ", num_squares, "\n");
    if (num_squares < rank_n() * _options->min_blocks_per_proc) {
      DBG("not enough squares ", num_squares, " < ", rank_n() * _options->min_blocks_per_proc,
          "\n");
      break;
    }
    // there is a remainder - this is not a perfect division
    if (remainder(num_grid_points, square)) {
      DBG("there is a remainder - don't use\n");
      continue;
    } else {
      DBG("selected dim ", i, "\n");
      block_dim = i;
    }
  }
  return block_dim;
}

Tissue::Tissue()
    : grid_points({})
    , new_active_grid_points({})
    , num_circulating_tcells(0)
    , tcells_generated({0}) {
  auto remainder = [](int64_t numerator, int64_t denominator) -> bool {
    return ((double)numerator / denominator - (numerator / denominator) != 0);
  };
  BarrierTimer timer(__FILEFUNC__, false, true);
  _grid_size = make_shared<GridCoords>(
      GridCoords(_options->dimensions[0], _options->dimensions[1], _options->dimensions[2]));
  int64_t num_grid_points = get_num_grid_points();
  // find the biggest cube that perfectly divides the grid and gives enough
  // data for at least two cubes per rank (for load
  // balance)
  // This is a trade-off: the more data is blocked, the better the locality,
  // but load balance could be a problem if not all
  // ranks get the same number of cubes. Also, having bigger cubes could lead
  // to load imbalance if all of the computation is
  // happening within a cube.
  int block_dim = (_grid_size->z > 1 ? get_cube_block_dim(num_grid_points) :
                                       get_square_block_dim(num_grid_points));
  if (block_dim == 1)
    SWARN("Using a block size of 1: this will result in a lot of "
          "communication. You should change the dimensions.");

  _grid_blocks.block_size =
      (_grid_size->z > 1 ? block_dim * block_dim * block_dim : block_dim * block_dim);
  _grid_blocks.num_x = _grid_size->x / block_dim;
  _grid_blocks.num_y = _grid_size->y / block_dim;
  _grid_blocks.num_z = (_grid_size->z > 1 ? _grid_size->z / block_dim : 1);
  _grid_blocks.size_x = block_dim;
  _grid_blocks.size_y = block_dim;
  _grid_blocks.size_z = (_grid_size->z > 1 ? block_dim : 1);

  int64_t num_blocks = num_grid_points / _grid_blocks.block_size;

  int64_t blocks_per_rank = ceil((double)num_blocks / rank_n());

  bool threeD = _grid_size->z > 1;
  SLOG("Dividing ", num_grid_points, " grid points into ", num_blocks,
       (threeD ? " blocks" : " squares"), " of size ", _grid_blocks.block_size, " (", block_dim,
       "^", (threeD ? 3 : 2), "), with ", blocks_per_rank, " per process\n");

  grid_points->reserve(blocks_per_rank * _grid_blocks.block_size);
  auto mem_reqd = sizeof(GridPoint) * blocks_per_rank * _grid_blocks.block_size;
  SLOG("Total initial memory required per process is a max of ", get_size_str(mem_reqd), "\n");
  // FIXME: it may be more efficient (less communication) to have contiguous blocks
  // this is the quick & dirty approach
  for (int64_t i = 0; i < blocks_per_rank; i++) {
    int64_t start_id = (i * rank_n() + rank_me()) * _grid_blocks.block_size;
    if (start_id >= num_grid_points) break;
    for (auto id = start_id; id < start_id + _grid_blocks.block_size; id++) {
      assert(id < num_grid_points);
      GridCoords coords(id);
      auto neighbors = get_neighbors(coords);
      // infectable epicells should be placed according to the underlying lung structure
      // (gaps, etc)
      EpiCell *epicell = new EpiCell(id);
      if ((coords.x + coords.y + coords.z) % _options->infectable_spacing != 0)
        epicell->infectable = false;
      grid_points->emplace_back(GridPoint({id, coords, neighbors, epicell}));
#ifdef DEBUG
      DBG("adding grid point ", id, " at ", coords.str(), "\n");
      auto id_1d = coords.to_1d();
      if (id_1d != id) DIE("id ", id, " is not same as returned by to_1d ", id_1d);
      ostringstream oss;
      for (auto nb_coords : neighbors) {
        oss << nb_coords.str() << " ";
      }
      DBG("nbs: ", oss.str(), "\n");
#endif
    }
  }
  barrier();
}

intrank_t Tissue::get_rank_for_grid_point(const GridCoords &coords) {
  int64_t id = coords.to_1d();
  int64_t block_i = id / _grid_blocks.block_size;
  return block_i % rank_n();
}

GridPoint *Tissue::get_local_grid_point(grid_points_t &grid_points, const GridCoords &coords) {
  auto id = coords.to_1d();
  int64_t block_i = id / _grid_blocks.block_size / rank_n();
  int64_t i = id % _grid_blocks.block_size + block_i * _grid_blocks.block_size;
  assert(i < grid_points->size());
  GridPoint *grid_point = &(*grid_points)[i];
  assert(grid_point->id == id);
  assert(grid_point->coords.x == coords.x);
  assert(grid_point->coords.y == coords.y);
  assert(grid_point->coords.z == coords.z);
  return grid_point;
}

SampleData Tissue::get_grid_point_sample_data(const GridCoords &coords) {
  return rpc(get_rank_for_grid_point(coords),
             [](grid_points_t &grid_points, GridCoords coords) {
               GridPoint *grid_point = Tissue::get_local_grid_point(grid_points, coords);
               SampleData sample;
               sample.id = grid_point->id;
               if (grid_point->tcell) sample.has_tcell = true;
               if (grid_point->epicell) {
                 sample.has_epicell = true;
                 sample.epicell_status = grid_point->epicell->status;
               }
               sample.virions = grid_point->virions;
               sample.chemokine = grid_point->chemokine;
               return sample;
             },
             grid_points, coords)
      .wait();
}

vector<GridCoords> Tissue::get_neighbors(GridCoords c) {
  vector<GridCoords> n = {};
  int64_t newx, newy, newz;
  for (int64_t i = -1; i <= 1; i++) {
    for (int64_t j = -1; j <= 1; j++) {
      for (int64_t k = -1; k <= 1; k++) {
        newx = c.x + i;
        newy = c.y + j;
        newz = c.z + k;
        if ((newx >= 0 && newx < _grid_size->x) && (newy >= 0 && newy < _grid_size->y) &&
            (newz >= 0 && newz < _grid_size->z)) {
          if (newx != c.x || newy != c.y || newz != c.z) n.push_back(GridCoords(newx, newy, newz));
        }
      }
    }
  }
  return n;
}

int64_t Tissue::get_num_local_grid_points() { return grid_points->size(); }

bool Tissue::set_initial_infection(GridCoords coords) {
  return rpc(get_rank_for_grid_point(coords),
             [](grid_points_t &grid_points, new_active_grid_points_t &new_active_grid_points,
                GridCoords coords) {
               GridPoint *grid_point = Tissue::get_local_grid_point(grid_points, coords);
               DBG("set infected for grid point ", grid_point, " ", grid_point->str(), "\n");
               // ensure we have at least on infected cell to start
               /*
               if (grid_point->epicell->status != EpiCellStatus::HEALTHY) return false;
               grid_point->epicell->infectable = true;
               grid_point->epicell->infect();
               */
               grid_point->virions = _options->initial_infection;
               new_active_grid_points->insert({grid_point, true});
               return true;
             },
             grid_points, new_active_grid_points, coords)
      .wait();
}

void Tissue::accumulate_chemokines(HASH_TABLE<int64_t, double> &chemokines_to_update,
                                   IntermittentTimer &timer) {
  timer.start();
  // accumulate updates for each target rank
  HASH_TABLE<intrank_t, vector<pair<GridCoords, double>>> target_rank_updates;
  for (auto &[coords_1d, chemokines] : chemokines_to_update) {
    progress();
    GridCoords coords(coords_1d);
    target_rank_updates[get_rank_for_grid_point(coords)].push_back({coords, chemokines});
  }
  future<> fut_chain = make_future<>();
  // dispatch all updates to each target rank in turn
  for (auto &[target_rank, update_vector] : target_rank_updates) {
    progress();
    auto fut = rpc(target_rank,
                   [](grid_points_t &grid_points, new_active_grid_points_t &new_active_grid_points,
                      view<pair<GridCoords, double>> update_vector) {
                     for (auto &[coords, chemokine] : update_vector) {
                       GridPoint *grid_point = Tissue::get_local_grid_point(grid_points, coords);
                       new_active_grid_points->insert({grid_point, true});
                       // just accumulate the concentrations. We will adjust them to be the average
                       // of all neighbors later
                       grid_point->nb_chemokine += chemokine;
                     }
                   },
                   grid_points, new_active_grid_points, make_view(update_vector));
    fut_chain = when_all(fut_chain, fut);
  }
  fut_chain.wait();
  timer.stop();
}

void Tissue::accumulate_virions(HASH_TABLE<int64_t, int> &virions_to_update,
                                IntermittentTimer &timer) {
  timer.start();
  // accumulate updates for each target rank
  HASH_TABLE<intrank_t, vector<pair<GridCoords, int>>> target_rank_updates;
  for (auto &[coords_1d, virions] : virions_to_update) {
    progress();
    GridCoords coords(coords_1d);
    target_rank_updates[get_rank_for_grid_point(coords)].push_back({coords, virions});
  }
  future<> fut_chain = make_future<>();
  // dispatch all updates to each target rank in turn
  for (auto &[target_rank, update_vector] : target_rank_updates) {
    progress();
    auto fut = rpc(target_rank,
                   [](grid_points_t &grid_points, new_active_grid_points_t &new_active_grid_points,
                      view<pair<GridCoords, int>> update_vector) {
                     for (auto &[coords, virions] : update_vector) {
                       GridPoint *grid_point = Tissue::get_local_grid_point(grid_points, coords);
                       new_active_grid_points->insert({grid_point, true});
                       grid_point->nb_virions += virions;
                     }
                   },
                   grid_points, new_active_grid_points, make_view(update_vector));
    fut_chain = when_all(fut_chain, fut);
  }
  fut_chain.wait();
  timer.stop();
}

double Tissue::get_chemokine(GridCoords coords) {
  return rpc(get_rank_for_grid_point(coords),
             [](grid_points_t &grid_points, GridCoords coords) {
               GridPoint *grid_point = Tissue::get_local_grid_point(grid_points, coords);
               return grid_point->chemokine;
             },
             grid_points, coords)
      .wait();
}

int64_t Tissue::get_num_circulating_tcells() { return num_circulating_tcells; }

void Tissue::change_num_circulating_tcells(int num) {
  num_circulating_tcells += num;
  if (num_circulating_tcells < 0) num_circulating_tcells = 0;
}

bool Tissue::try_add_new_tissue_tcell(GridCoords coords) {
  auto res = rpc(get_rank_for_grid_point(coords),
                 [](grid_points_t &grid_points, new_active_grid_points_t &new_active_grid_points,
                    GridCoords coords, dist_object<int64_t> &tcells_generated) {
                   GridPoint *grid_point = Tissue::get_local_grid_point(grid_points, coords);
                   // grid point is already occupied by a tcell, don't add
                   if (grid_point->tcell) return false;
                   if (grid_point->chemokine < _options->min_chemokine) return false;
                   new_active_grid_points->insert({grid_point, true});
                   string tcell_id = to_string(rank_me()) + "-" + to_string(*tcells_generated);
                   (*tcells_generated)++;
                   grid_point->tcell = new TCell(tcell_id);
                   grid_point->tcell->moved = true;
                   return true;
                 },
                 grid_points, new_active_grid_points, coords, tcells_generated)
                 .wait();
  if (res) num_circulating_tcells--;
  assert(num_circulating_tcells >= 0);
  return res;
}

bool Tissue::try_add_tissue_tcell(GridCoords coords, TCell &tcell) {
  return rpc(get_rank_for_grid_point(coords),
             [](grid_points_t &grid_points, new_active_grid_points_t &new_active_grid_points,
                GridCoords coords, TCell tcell) {
               GridPoint *grid_point = Tissue::get_local_grid_point(grid_points, coords);
               // grid point is already occupied by a tcell, don't add
               if (grid_point->tcell) return false;
               new_active_grid_points->insert({grid_point, true});
               tcell.moved = true;
               grid_point->tcell = new TCell(tcell);
               return true;
             },
             grid_points, new_active_grid_points, coords, tcell)
      .wait();
}

EpiCellStatus Tissue::try_bind_tcell(GridCoords coords) {
  return rpc(get_rank_for_grid_point(coords),
             [](grid_points_t &grid_points, new_active_grid_points_t &new_active_grid_points_t,
                GridCoords coords) {
               GridPoint *grid_point = Tissue::get_local_grid_point(grid_points, coords);
               if (!grid_point->epicell) return EpiCellStatus::DEAD;
               if (grid_point->epicell->status == EpiCellStatus::HEALTHY ||
                   grid_point->epicell->status == EpiCellStatus::DEAD)
                 return grid_point->epicell->status;

               // if (grid_point->epicell->status == EpiCellStatus::DEAD) return
               // EpiCellStatus::DEAD;

               double binding_prob = grid_point->epicell->get_binding_prob();
               if (_rnd_gen->trial_success(binding_prob)) {
                 auto prev_status = grid_point->epicell->status;
                 grid_point->epicell->status = EpiCellStatus::APOPTOTIC;
                 return prev_status;
               }
               return EpiCellStatus::DEAD;
             },
             grid_points, new_active_grid_points, coords)
      .wait();
}

void Tissue::get_samples(vector<SampleData> &samples, int64_t &start_id) {
  int64_t num_points = get_num_grid_points();
  int64_t num_points_per_rank = ceil((double)num_points / rank_n());
  start_id = rank_me() * num_points_per_rank;
  int64_t end_id = min((rank_me() + 1) * num_points_per_rank, num_points);
  size_t buf_size = end_id - start_id;
  samples.clear();
  samples.reserve(buf_size);
  // FIXME: this should be aggregated fetches per target rank
  for (int64_t id = start_id; id < end_id; id++) {
    GridCoords coords;
    coords.from_1d_linear(id);
    samples.emplace_back(get_grid_point_sample_data(coords));
  }
}

GridPoint *Tissue::get_first_local_grid_point() {
  grid_point_iter = grid_points->begin();
  if (grid_point_iter == grid_points->end()) return nullptr;
  auto grid_point = &(*grid_point_iter);
  ++grid_point_iter;
  return grid_point;
}

GridPoint *Tissue::get_next_local_grid_point() {
  if (grid_point_iter == grid_points->end()) return nullptr;
  auto grid_point = &(*grid_point_iter);
  ++grid_point_iter;
  return grid_point;
}

GridPoint *Tissue::get_first_active_grid_point() {
  active_grid_point_iter = active_grid_points.begin();
  if (active_grid_point_iter == active_grid_points.end()) return nullptr;
  auto grid_point = active_grid_point_iter->first;
  ++active_grid_point_iter;
  return grid_point;
}

GridPoint *Tissue::get_next_active_grid_point() {
  if (active_grid_point_iter == active_grid_points.end()) return nullptr;
  auto grid_point = active_grid_point_iter->first;
  ++active_grid_point_iter;
  return grid_point;
}

void Tissue::set_active(GridPoint *grid_point) {
  new_active_grid_points->insert({grid_point, true});
}

void Tissue::erase_active(GridPoint *grid_point) { active_grid_points.erase(grid_point); }

void Tissue::add_new_actives(IntermittentTimer &timer) {
  timer.start();
  DBG("add ", new_active_grid_points->size(), " new active grid points\n");
  for (auto elem : *new_active_grid_points) {
    // DBG("inserting from new active ", elem.first, " ", elem.first->str(), "\n");
    active_grid_points.insert(elem);
  }
  new_active_grid_points->clear();
  timer.stop();
}

size_t Tissue::get_num_actives() { return active_grid_points.size(); }

#ifdef DEBUG
void Tissue::check_actives(int time_step) {
  for (int64_t i = 0; i < grid_points->size(); i++) {
    GridPoint *grid_point = &(*grid_points)[i];
    if (time_step == 0)
      DBG("coords ", grid_point->coords.str(), " to 1d ", grid_point->coords.to_1d(), "\n");
    bool found_active = (active_grid_points.find(grid_point) != active_grid_points.end());
    if (grid_point->is_active() && !found_active)
      DIE(time_step, ": active grid point ", grid_point->str(), " not found in active list");
    if (!grid_point->is_active() && found_active)
      DIE(time_step, ": inactive grid point ", grid_point->str(), " found in active list");
  }
}
#endif
