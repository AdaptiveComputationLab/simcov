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
  z = i / (_grid_size->x * _grid_size->y);
  i = i % (_grid_size->x * _grid_size->y);
  y = i / _grid_size->x;
  x = i % _grid_size->x;
#endif
}

GridCoords::GridCoords(shared_ptr<Random> rnd_gen) {
  x = rnd_gen->get(0, _grid_size->x);
  y = rnd_gen->get(0, _grid_size->y);
  z = rnd_gen->get(0, _grid_size->z);
}

int64_t GridCoords::to_1d(int x, int y, int z) {
  if (x >= _grid_size->x || y >= _grid_size->y || z >= _grid_size->z)
    DIE("Grid point is out of range: ", x, " ", y, " ", z, " max size ", _grid_size->str());
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
  return (int64_t)x + (int64_t)y * _grid_size->x + (int64_t)z * _grid_size->x * _grid_size->y;
#endif
}

int64_t GridCoords::linear_to_block(int64_t i) {
  int z = i / (_grid_size->x * _grid_size->y);
  i = i % (_grid_size->x * _grid_size->y);
  int y = i / _grid_size->x;
  int x = i % _grid_size->x;
  return GridCoords::to_1d(x, y, z);
}

int64_t GridCoords::to_1d() const { return GridCoords::to_1d(x, y, z); }

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
  oss << "xyz " << coords.str() << ", epi " << (epicell ? epicell->str() : "none") << ", v "
      << virions << ", c " << chemokine;
  return oss.str();
}

bool GridPoint::is_active() {
  // it could be incubating but without anything else set
  return ((epicell && epicell->is_active()) || virions > 0 || chemokine > 0 || tcell);
}

static int get_cube_block_dim(int64_t num_grid_points) {
  int block_dim = 1;
  for (int d = 1; d < _options->max_block_dim; d++) {
    // only allow dims that divide each main dimension perfectly
    if (remainder(_grid_size->x, d) || remainder(_grid_size->y, d) || remainder(_grid_size->z, d)) {
      DBG("dim ", d, " does not divide all main dimensions cleanly\n");
      continue;
    }
    size_t cube = (size_t)pow((double)d, 3.0);
    size_t num_cubes = num_grid_points / cube;
    DBG("cube size ", cube, " num cubes ", num_cubes, "\n");
    if (num_cubes < rank_n() * MIN_BLOCKS_PER_PROC) {
      DBG("not enough cubes ", num_cubes, " < ", rank_n() * MIN_BLOCKS_PER_PROC, "\n");
      break;
    }
    // there is a remainder - this is not a perfect division
    if (remainder(num_grid_points, cube)) {
      DBG("there is a remainder - don't use\n");
      continue;
    }
    // skip sizes that distribute the blocks in columns
    if (d > 1 && (_grid_size->x % (d * rank_n()) == 0 || _grid_size->y % (d * rank_n()) == 0)) {
      DBG("dim ", d, " gives perfect division of all blocks into x axis - skip\n");
      continue;
    }
    DBG("selected dim ", d, "\n");
    block_dim = d;
  }
  return block_dim;
}

static int get_square_block_dim(int64_t num_grid_points) {
  int block_dim = 1;
  for (int d = 1; d < _options->max_block_dim; d++) {
    // only allow dims that divide each main dimension perfectly
    if (remainder(_grid_size->x, d) || remainder(_grid_size->y, d)) {
      DBG("dim ", d, " does not divide all main dimensions cleanly\n");
      continue;
    }
    size_t square = (size_t)pow((double)d, 2.0);
    size_t num_squares = num_grid_points / square;
    DBG("square size ", square, " num squares ", num_squares, "\n");
    if (num_squares < rank_n() * MIN_BLOCKS_PER_PROC) {
      DBG("not enough squares ", num_squares, " < ", rank_n() * MIN_BLOCKS_PER_PROC, "\n");
      break;
    }
    // there is a remainder - this is not a perfect division
    if (remainder(num_grid_points, square)) {
      DBG("there is a remainder - don't use\n");
      continue;
    }
    // skip sizes that distribute the blocks in columns
    if (d > 1 && _grid_size->x % (d * rank_n()) == 0) {
      DBG("dim ", d, " gives perfect division of all blocks into x axis - skip\n");
      continue;
    }
    DBG("selected dim ", d, "\n");
    block_dim = d;
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
  int64_t block_dim = (_grid_size->z > 1 ? get_cube_block_dim(num_grid_points) :
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
  double sz_grid_point = sizeof(GridPoint) + (double)sizeof(EpiCell);
  auto mem_reqd = sz_grid_point * blocks_per_rank * _grid_blocks.block_size;
  SLOG("Total initial memory required per process is at least ", get_size_str(mem_reqd),
       " with each grid point requiring on average ", sz_grid_point, " bytes\n");
  int64_t num_lung_cells = 0;
  if (!_options->lung_model_dir.empty()) {
    lung_cells.resize(num_grid_points, EpiCellType::NONE);
    Timer t_load_lung_model("load lung model");
    t_load_lung_model.start();
    // Read alveolus epithileal cells
    num_lung_cells += load_data_file(_options->lung_model_dir + "/alveolus.dat", num_grid_points,
                                     EpiCellType::ALVEOLI);
    // Read bronchiole epithileal cells
    num_lung_cells += load_data_file(_options->lung_model_dir + "/bronchiole.dat", num_grid_points,
                                     EpiCellType::AIRWAY);

    t_load_lung_model.stop();
    SLOG("Lung model loaded ", num_lung_cells, " epithileal cells in ", fixed, setprecision(2),
         t_load_lung_model.get_elapsed(), " s\n");
  }

  // FIXME: these blocks need to be stride distributed to better load balance
  grid_points->reserve(blocks_per_rank * _grid_blocks.block_size);
  for (int64_t i = 0; i < blocks_per_rank; i++) {
    int64_t start_id = (i * rank_n() + rank_me()) * _grid_blocks.block_size;
    if (start_id >= num_grid_points) break;
    for (auto id = start_id; id < start_id + _grid_blocks.block_size; id++) {
      assert(id < num_grid_points);
      GridCoords coords(id);
      if (num_lung_cells) {
        if (lung_cells[id] != EpiCellType::NONE) {
          EpiCell *epicell = new EpiCell(id);
          epicell->type = lung_cells[id];
          epicell->infectable = true;
          grid_points->emplace_back(GridPoint({coords, epicell}));
        } else {  // Add empty space == air
          grid_points->emplace_back(GridPoint({coords, nullptr}));
        }
      } else {
        EpiCell *epicell = new EpiCell(id);
        epicell->type = EpiCellType::ALVEOLI;
        // epicell->status = static_cast<EpiCellStatus>(rank_me() % 4);
        epicell->infectable = true;
        grid_points->emplace_back(GridPoint({coords, epicell}));
      }
#ifdef DEBUG
      DBG("adding grid point ", id, " at ", coords.str(), "\n");
      auto id_1d = coords.to_1d();
      if (id_1d != id) DIE("id ", id, " is not same as returned by to_1d ", id_1d);
      auto nbs = get_neighbors(coords);
      ostringstream oss;
      for (auto nb_grid_i : *nbs) {
        oss << GridCoords(nb_grid_i).str() << " ";
      }
      DBG("nbs: ", oss.str(), "\n");
#endif
    }
  }
  barrier();
}

int Tissue::load_data_file(const string &fname, int num_grid_points, EpiCellType epicell_type) {
  ifstream f(fname, ios::in | ios::binary);
  if (!f) SDIE("Couldn't open file ", fname);
  f.seekg(0, ios::end);
  auto fsize = f.tellg();
  auto num_ids = fsize / sizeof(int);
  if (num_ids > num_grid_points + 3) DIE("Too many ids in ", fname, " max is ", num_grid_points);
  f.clear();
  f.seekg(0, ios::beg);
  vector<int> id_buf(num_ids);
  if (!f.read(reinterpret_cast<char *>(&(id_buf[0])), fsize))
    DIE("Couldn't read all bytes in ", fname);
  int num_lung_cells = 0;
  // skip first three wwhich are dimensions
  for (int i = 3; i < id_buf.size(); i++) {
    auto id = id_buf[i];
#ifdef BLOCK_PARTITION
    id = GridCoords::linear_to_block(id);
#endif
    lung_cells[id] = epicell_type;
    num_lung_cells++;
  }
  f.close();
  return num_lung_cells;
}

intrank_t Tissue::get_rank_for_grid_point(int64_t grid_i) {
  int64_t block_i = grid_i / _grid_blocks.block_size;
  return block_i % rank_n();
}

GridPoint *Tissue::get_local_grid_point(grid_points_t &grid_points, int64_t grid_i) {
  int64_t block_i = grid_i / _grid_blocks.block_size / rank_n();
  int64_t i = grid_i % _grid_blocks.block_size + block_i * _grid_blocks.block_size;
  assert(i < grid_points->size());
  GridPoint *grid_point = &(*grid_points)[i];
  if (grid_point->coords.to_1d() != grid_i)
    DIE("mismatched coords to grid i ", grid_point->coords.to_1d(), " != ", grid_i);
  return grid_point;
}

SampleData Tissue::get_grid_point_sample_data(int64_t grid_i) {
  return rpc(
             get_rank_for_grid_point(grid_i),
             [](grid_points_t &grid_points, int64_t grid_i) {
               GridPoint *grid_point = Tissue::get_local_grid_point(grid_points, grid_i);
               SampleData sample;
               if (grid_point->tcell) sample.tcells = 1;
               if (grid_point->epicell) {
                 sample.has_epicell = true;
                 sample.epicell_status = grid_point->epicell->status;
               }
               sample.virions = grid_point->virions;
               sample.chemokine = grid_point->chemokine;
               return sample;
             },
             grid_points, grid_i)
      .wait();
}

vector<int64_t> *Tissue::get_neighbors(GridCoords c) {
  GridPoint *grid_point = Tissue::get_local_grid_point(grid_points, c.to_1d());
  if (!grid_point->neighbors) {
    grid_point->neighbors = new vector<int64_t>;
    int newx, newy, newz;
    for (int i = -1; i <= 1; i++) {
      for (int j = -1; j <= 1; j++) {
        for (int k = -1; k <= 1; k++) {
          newx = c.x + i;
          newy = c.y + j;
          newz = c.z + k;
          if ((newx >= 0 && newx < _grid_size->x) && (newy >= 0 && newy < _grid_size->y) &&
              (newz >= 0 && newz < _grid_size->z)) {
            if (newx != c.x || newy != c.y || newz != c.z) {
              grid_point->neighbors->push_back(GridCoords::to_1d(newx, newy, newz));
            }
          }
        }
      }
    }
  }
  return grid_point->neighbors;
}

int64_t Tissue::get_num_local_grid_points() { return grid_points->size(); }

/*
int64_t Tissue::get_random_airway_epicell_location() {
  std::set<int>::iterator it = airway.begin();
  std::advance(it, airway.size() / 2);
  return *it;
}
*/

bool Tissue::set_initial_infection(int64_t grid_i) {
  return rpc(
             get_rank_for_grid_point(grid_i),
             [](grid_points_t &grid_points, new_active_grid_points_t &new_active_grid_points,
                int64_t grid_i) {
               GridPoint *grid_point = Tissue::get_local_grid_point(grid_points, grid_i);
               DBG("set infected for grid point ", grid_point, " ", grid_point->str(), "\n");
               if (!grid_point->epicell) return false;
               grid_point->virions = _options->initial_infection;
               new_active_grid_points->insert({grid_point, true});
               return true;
             },
             grid_points, new_active_grid_points, grid_i)
      .wait();
}

void Tissue::accumulate_chemokines(HASH_TABLE<int64_t, float> &chemokines_to_update,
                                   IntermittentTimer &timer) {
  timer.start();
  // accumulate updates for each target rank
  HASH_TABLE<intrank_t, vector<pair<int64_t, float>>> target_rank_updates;
  for (auto &[coords_1d, chemokines] : chemokines_to_update) {
    progress();
    target_rank_updates[get_rank_for_grid_point(coords_1d)].push_back({coords_1d, chemokines});
  }
  future<> fut_chain = make_future<>();
  // dispatch all updates to each target rank in turn
  for (auto &[target_rank, update_vector] : target_rank_updates) {
    progress();
    auto fut = rpc(
        target_rank,
        [](grid_points_t &grid_points, new_active_grid_points_t &new_active_grid_points,
           view<pair<int64_t, float>> update_vector) {
          for (auto &[grid_i, chemokine] : update_vector) {
            GridPoint *grid_point = Tissue::get_local_grid_point(grid_points, grid_i);
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

void Tissue::accumulate_virions(HASH_TABLE<int64_t, float> &virions_to_update,
                                IntermittentTimer &timer) {
  timer.start();
  // accumulate updates for each target rank
  HASH_TABLE<intrank_t, vector<pair<int64_t, float>>> target_rank_updates;
  for (auto &[coords_1d, virions] : virions_to_update) {
    progress();
    target_rank_updates[get_rank_for_grid_point(coords_1d)].push_back({coords_1d, virions});
  }
  future<> fut_chain = make_future<>();
  // dispatch all updates to each target rank in turn
  for (auto &[target_rank, update_vector] : target_rank_updates) {
    progress();
    auto fut = rpc(
        target_rank,
        [](grid_points_t &grid_points, new_active_grid_points_t &new_active_grid_points,
           view<pair<int64_t, float>> update_vector) {
          for (auto &[grid_i, virions] : update_vector) {
            GridPoint *grid_point = Tissue::get_local_grid_point(grid_points, grid_i);
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

float Tissue::get_chemokine(int64_t grid_i) {
  return rpc(
             get_rank_for_grid_point(grid_i),
             [](grid_points_t &grid_points, int64_t grid_i) {
               GridPoint *grid_point = Tissue::get_local_grid_point(grid_points, grid_i);
               return grid_point->chemokine;
             },
             grid_points, grid_i)
      .wait();
}

int64_t Tissue::get_num_circulating_tcells() { return num_circulating_tcells; }

void Tissue::change_num_circulating_tcells(int num) {
  num_circulating_tcells += num;
  if (num_circulating_tcells < 0) num_circulating_tcells = 0;
}

bool Tissue::try_add_new_tissue_tcell(int64_t grid_i) {
  auto res = rpc(
                 get_rank_for_grid_point(grid_i),
                 [](grid_points_t &grid_points, new_active_grid_points_t &new_active_grid_points,
                    int64_t grid_i, dist_object<int64_t> &tcells_generated) {
                   GridPoint *grid_point = Tissue::get_local_grid_point(grid_points, grid_i);
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
                 grid_points, new_active_grid_points, grid_i, tcells_generated)
                 .wait();
  if (res) num_circulating_tcells--;
  assert(num_circulating_tcells >= 0);
  return res;
}

bool Tissue::try_add_tissue_tcell(int64_t grid_i, TCell &tcell) {
  return rpc(
             get_rank_for_grid_point(grid_i),
             [](grid_points_t &grid_points, new_active_grid_points_t &new_active_grid_points,
                int64_t grid_i, TCell tcell) {
               GridPoint *grid_point = Tissue::get_local_grid_point(grid_points, grid_i);
               // grid point is already occupied by a tcell, don't add
               if (grid_point->tcell) return false;
               new_active_grid_points->insert({grid_point, true});
               tcell.moved = true;
               grid_point->tcell = new TCell(tcell);
               return true;
             },
             grid_points, new_active_grid_points, grid_i, tcell)
      .wait();
}

EpiCellStatus Tissue::try_bind_tcell(int64_t grid_i) {
  return rpc(
             get_rank_for_grid_point(grid_i),
             [](grid_points_t &grid_points, new_active_grid_points_t &new_active_grid_points_t,
                int64_t grid_i) {
               GridPoint *grid_point = Tissue::get_local_grid_point(grid_points, grid_i);
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
             grid_points, new_active_grid_points, grid_i)
      .wait();
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
