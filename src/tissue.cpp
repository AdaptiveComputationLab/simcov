#include "tissue.hpp"

GridCoords::GridCoords(int64_t i, const GridCoords &grid_size) {
  z = i / (grid_size.x * grid_size.y);
  i = i % (grid_size.x * grid_size.y);
  y = i / grid_size.x;
  x = i % grid_size.x;
}

GridCoords::GridCoords(std::shared_ptr<Random> rnd_gen, const GridCoords &grid_size) {
  x = rnd_gen->get(0, grid_size.x);
  y = rnd_gen->get(0, grid_size.y);
  z = rnd_gen->get(0, grid_size.z);
}

string GridPoint::str() const {
  ostringstream oss;
  oss << id << " " << coords.str() << " " << epicell->str() << " " << virus;
  return oss.str();
}

void GridPoint::switch_tcells_vector() {
  auto switch_vectors = [](vector<TCell> *tcells_backing_new,
                           vector<TCell> *tcells_backing_old) -> vector<TCell> * {
    tcells_backing_old->clear();
    return tcells_backing_new;
  };
  if (tcells == &tcells_backing_1)
    tcells = switch_vectors(&tcells_backing_2, &tcells_backing_1);
  else
    tcells = switch_vectors(&tcells_backing_1, &tcells_backing_2);
}

void GridPoint::add_tcell(TCell tcell) {
  auto next_tcells = (tcells == &tcells_backing_1 ? &tcells_backing_2 : &tcells_backing_1);
  next_tcells->push_back(tcell);
}

intrank_t Tissue::get_rank_for_grid_point(const GridCoords &coords) {
  int64_t id = coords.to_1d(Tissue::grid_size);
  int64_t block_i = id / Tissue::block_size;
  return block_i % rank_n();
}

GridPoint &Tissue::get_local_grid_point(grid_points_t &grid_points, int64_t id,
                                        const GridCoords &coords) {
  int64_t block_i = id / Tissue::block_size / rank_n();
  int64_t i = id % Tissue::block_size + block_i * Tissue::block_size;
  assert(i < grid_points->size());
  GridPoint &grid_point = (*grid_points)[i];
  assert(grid_point.id == id);
  assert(grid_point.coords.x == coords.x);
  assert(grid_point.coords.y == coords.y);
  assert(grid_point.coords.z == coords.z);
  return grid_point;
}

vector<GridCoords> Tissue::get_neighbors(GridCoords c, GridCoords grid_size) {
  vector<GridCoords> n = {};
  int64_t newx, newy, newz;
  for (int64_t i = -1; i <= 1; i++) {
    for (int64_t j = -1; j <= 1; j++) {
      for (int64_t k = -1; k <= 1; k++) {
        newx = c.x + i;
        newy = c.y + j;
        newz = c.z + k;
        if ((newx >= 0 && newx < grid_size.x) && (newy >= 0 && newy < grid_size.y) &&
            (newz >= 0 && newz < grid_size.z)) {
          n.push_back(GridCoords(newx, newy, newz));
        }
      }
    }
  }
  return n;
}

int64_t Tissue::get_num_grid_points() {
  return Tissue::grid_size.x * Tissue::grid_size.y * Tissue::grid_size.z;
}

void Tissue::inc_incoming_virus(GridCoords coords, double virus) {
  upcxx::rpc(
      get_rank_for_grid_point(coords),
      [](grid_points_t &grid_points, int64_t id, GridCoords coords, double virus) {
        GridPoint &grid_point = Tissue::get_local_grid_point(grid_points, id, coords);
        if (grid_point.virus == 0 && grid_point.epicell->status == EpiCellStatus::HEALTHY) {
          grid_point.incoming_virus += virus;
          if (grid_point.incoming_virus > 1) grid_point.incoming_virus = 1;
        }
      },
      grid_points, coords.to_1d(Tissue::grid_size), coords, virus)
      .wait();
}

void Tissue::inc_incoming_chemokines(GridCoords coords, double chemokine) {
  upcxx::rpc(
      get_rank_for_grid_point(coords),
      [](grid_points_t &grid_points, int64_t id, GridCoords coords, double chemokine) {
        GridPoint &grid_point = Tissue::get_local_grid_point(grid_points, id, coords);
        if (grid_point.epicell->status != EpiCellStatus::EXPRESSING) {
          grid_point.incoming_chemokine += chemokine;
          //grid_point.incoming_chemokine = std::max(chemokine, grid_point.incoming_chemokine);
          if (grid_point.incoming_chemokine > 1) grid_point.incoming_chemokine = 1;
        }
      },
      grid_points, coords.to_1d(Tissue::grid_size), coords, chemokine)
      .wait();
}

void Tissue::inc_incoming_icytokines(GridCoords coords, double icytokine) {
  upcxx::rpc(
      get_rank_for_grid_point(coords),
      [](grid_points_t &grid_points, int64_t id, GridCoords coords, double icytokine) {
        GridPoint &grid_point = Tissue::get_local_grid_point(grid_points, id, coords);
        if (grid_point.epicell->status != EpiCellStatus::EXPRESSING) {
          grid_point.incoming_icytokine += icytokine;
          //grid_point.incoming_icytokine = std::max(icytokine, grid_point.incoming_icytokine);
          if (grid_point.incoming_icytokine > 1) grid_point.incoming_icytokine = 1;
        }
      },
      grid_points, coords.to_1d(Tissue::grid_size), coords, icytokine)
      .wait();
}

double Tissue::get_chemokine(GridCoords coords) {
  return upcxx::rpc(
             get_rank_for_grid_point(coords),
             [](grid_points_t &grid_points, int64_t id, GridCoords coords) {
               GridPoint &grid_point = Tissue::get_local_grid_point(grid_points, id, coords);
               if (grid_point.epicell->status != EpiCellStatus::DEAD) return 0.0;
               return grid_point.chemokine;
             },
             grid_points, coords.to_1d(Tissue::grid_size), coords)
      .wait();
}

void Tissue::add_tcell(GridCoords coords, TCell tcell) {
  // FIXME: this should be an aggregating store
  upcxx::rpc(
      get_rank_for_grid_point(coords),
      [](grid_points_t &grid_points, int64_t id, GridCoords coords, TCell tcell) {
        GridPoint &grid_point = Tissue::get_local_grid_point(grid_points, id, coords);
        assert(grid_point.tcells != nullptr);
        // add the tcell to the future tcells vector, i.e. not
        // the one pointed to by tcells
        grid_point.add_tcell(tcell);
      },
      grid_points, coords.to_1d(Tissue::grid_size), coords, tcell)
      .wait();
}

void Tissue::construct(GridCoords grid_size) {
  auto remainder = [](int64_t numerator, int64_t denominator) -> bool {
    return ((double)numerator / denominator - (numerator / denominator) != 0);
  };
  BarrierTimer timer(__FILEFUNC__, false, true);
  Tissue::grid_size = grid_size;
  int64_t num_grid_points = get_num_grid_points();
  // find the biggest cube that perfectly divides the grid and gives enough
  // data for at least two cubes per rank (for load
  // balance)
  // This is a trade-off: the more data is blocked, the better the locality,
  // but load balance could be a problem if not all
  // ranks get the same number of cubes. Also, having bigger cubes could lead
  // to load imbalance if all of the computation is
  // happening within a cube.
  int min_cubes_per_rank = 2;
  int block_dim = 1;
  int min_dim = std::min(std::min(Tissue::grid_size.x, Tissue::grid_size.y), Tissue::grid_size.z);
  for (int i = 1; i < min_dim; i++) {
    // only allow dims that divide each main dimension perfectly
    if (remainder(Tissue::grid_size.x, i) || remainder(Tissue::grid_size.y, i) ||
        remainder(Tissue::grid_size.z, i)) {
      DBG("dim ", i, " does not divide all main dimensions cleanly\n");
      continue;
    }
    size_t cube = (size_t)pow((double)i, 3.0);
    size_t num_cubes = num_grid_points / cube;
    DBG("cube size ", cube, " num cubes ", num_cubes, "\n");
    if (num_cubes < rank_n() * min_cubes_per_rank) {
      DBG("not enough cubes ", num_cubes, " < ", rank_n() * min_cubes_per_rank, "\n");
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
  Tissue::block_size = block_dim * block_dim * block_dim;
  if (block_dim == 1)
    SWARN("Using a block size of 1: this will result in a lot of "
          "communication. You should change the dimensions.");
  int64_t x_blocks = Tissue::grid_size.x / block_dim;
  int64_t y_blocks = Tissue::grid_size.y / block_dim;
  int64_t z_blocks = Tissue::grid_size.z / block_dim;
  int64_t num_blocks = num_grid_points / Tissue::block_size;
  SLOG("Dividing ", num_grid_points, " grid points into ", num_blocks, " blocks of size ",
       Tissue::block_size, " (", block_dim, "^3)\n");
  int64_t blocks_per_rank = ceil((double)num_blocks / rank_n());
  SLOG_VERBOSE("Each process has ", blocks_per_rank, " blocks\n");
  grid_points->reserve(blocks_per_rank * Tissue::block_size);
  auto mem_reqd = sizeof(GridPoint) * blocks_per_rank * Tissue::block_size;
  SLOG("Total initial memory required per process is a max of ", get_size_str(mem_reqd), "\n");
  // FIXME: it may be more efficient (less communication) to have contiguous
  // blocks
  // this is the quick & dirty approach
  for (int64_t i = 0; i < blocks_per_rank; i++) {
    int64_t start_id = (i * rank_n() + rank_me()) * Tissue::block_size;
    if (start_id >= num_grid_points) break;
    for (auto id = start_id; id < start_id + Tissue::block_size; id++) {
      assert(id < num_grid_points);
      GridCoords coords(id, Tissue::grid_size);
      auto neighbors = get_neighbors(coords, grid_size);
      // FIXME: this is an epicall on every grid point.
      // They should be placed according to the underlying lung structure
      // (gaps, etc)
      GridPoint grid_point;
      grid_points->push_back(std::move(grid_point));
      grid_points->back().init(
          id, coords, neighbors,
          new EpiCell({.id = id, .status = EpiCellStatus::HEALTHY, .num_steps_infected = 0}));
      DBG("adding grid point ", id, " at ", coords.str(), "\n");
    }
  }
  /*
  #ifdef DEBUG
      for (auto &grid_point : (*grid_points)) {
        DBG("grid point ", grid_point.str(), "\n");
      }
  #endif
  */
  barrier();
}

std::pair<size_t, size_t> Tissue::dump_blocks(const string &fname, const string &header_str,
                                              ViewObject view_object) {
  auto fileno = open(fname.c_str(), O_WRONLY, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
  if (fileno == -1) DIE("Cannot open file ", fname, ": ", strerror(errno), "\n");
  if (!upcxx::rank_me()) {
    auto bytes_written = pwrite(fileno, header_str.c_str(), header_str.length(), 0);
    if (bytes_written != header_str.length())
      DIE("Could not write all ", header_str.length(), " bytes: only wrote ", bytes_written);
  }
  DBG("Writing samples to ", fname, "\n");
  DBG("header size is ", header_str.length(), "\n");
  size_t tot_bytes_written = 0;
  size_t grid_points_written = 0;
  int64_t num_grid_points = get_num_grid_points();
  int64_t num_blocks = num_grid_points / Tissue::block_size;
  int64_t blocks_per_rank = ceil((double)num_blocks / rank_n());
  unsigned char *buf = new unsigned char[Tissue::block_size + 1];
  for (int64_t i = 0; i < blocks_per_rank; i++) {
    int64_t start_id = (i * rank_n() + rank_me()) * Tissue::block_size;
    if (start_id >= num_grid_points) break;
    for (auto id = start_id; id < start_id + Tissue::block_size; id++) {
      assert(id < num_grid_points);
      GridCoords coords(id, Tissue::grid_size);
      GridPoint &grid_point = Tissue::get_local_grid_point(grid_points, id, coords);
      // a dead cell is 0 value (blue), otherwise a normal cell is 127, which
      // is grey. The view objects add to the 127 color.
      int val = 0;
      switch (view_object) {
        case ViewObject::TCELL_VAS:
          for (auto tcell : *grid_point.tcells) {
            if (tcell.in_vasculature) {
              val++;
              if (val == 255) break;
            }
          }
          break;
        case ViewObject::TCELL_TISSUE:
          for (auto tcell : *grid_point.tcells) {
            if (!tcell.in_vasculature) {
              val++;
              if (val == 255) break;
            }
          }
          break;
        case ViewObject::VIRUS:
          if (grid_point.virus < 0)
            val = 0;
          else
            val = 255 * grid_point.virus;
          break;
        case ViewObject::EPICELL:
          switch (grid_point.epicell->status) {
            case EpiCellStatus::INCUBATING: val = 1; break;
            case EpiCellStatus::EXPRESSING: val = 2; break;
            case EpiCellStatus::APOPTOTIC: val = 3; break;
            case EpiCellStatus::DEAD: val = 4; break;
            default: val = 0;
          }
          break;
        case ViewObject::CHEMOKINE:
          if (grid_point.chemokine < 0)
            val = 0;
          else
            val = 255 * grid_point.chemokine;
          break;
        case ViewObject::ICYTOKINE:
          if (grid_point.icytokine < 0)
            val = 0;
          else
            val = 255 * grid_point.icytokine;
          break;
      }
      buf[id - start_id] = (unsigned char)val;
      grid_points_written++;
    }
    buf[Tissue::block_size] = 0;
    size_t fpos = start_id + header_str.length();
    auto bytes_written = pwrite(fileno, buf, Tissue::block_size, fpos);
    if (bytes_written != Tissue::block_size)
      DIE("Could not write all ", Tissue::block_size, " bytes; only wrote ", bytes_written, "\n");
    tot_bytes_written += bytes_written;
    DBG("wrote block ", i, ", ", bytes_written, " bytes at position ", fpos, "\n");
  }
  delete[] buf;
  close(fileno);
  return {tot_bytes_written, grid_points_written};
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
