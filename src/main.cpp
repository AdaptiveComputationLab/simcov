// SimCov
//
// Steven Hofmeyr, LBNL May 2020

#include <iostream>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <string>
#include <math.h>
#include <stdarg.h>
#include <unistd.h>
#include <fcntl.h>
#include <upcxx/upcxx.hpp>

using namespace std;

#include "upcxx_utils.hpp"
#include "options.hpp"
#include "utils.hpp"
#include "tissue.hpp"


using namespace upcxx;
using namespace upcxx_utils;

ofstream _logstream;
bool _verbose = false;
shared_ptr<Options> _options;
shared_ptr<Random> _rnd_gen;

IntermittentTimer _update_tcell_timer(__FILENAME__ + string(":") + "update_tcell");
IntermittentTimer _update_virus_timer(__FILENAME__ + string(":") + "update_virus");
IntermittentTimer _sample_timer(__FILENAME__ + string(":") + "sample");


int64_t initial_infection(int num_infections, Tissue &tissue) {
  BarrierTimer timer(__FILEFUNC__, false, true);
  // each rank generates a block of infected cells
  // for large dimensions, the chance of repeat sampling for a few points is very small
  int local_num_infections = ceil((double)num_infections / rank_n());
  int64_t min_i = rank_me() * local_num_infections;
  int64_t max_i = (rank_me() + 1) * local_num_infections;
  if (max_i > num_infections) max_i = num_infections;
  if (min_i > num_infections) min_i = num_infections;
  int64_t num_infected = 0;
  barrier();
  DBG("Infecting ", (max_i - min_i), " cells from ", min_i, " to ", max_i, "\n");
  future<> fchain = make_future<>();
  ProgressBar progbar(max_i - min_i, "Setting initial infections");
  for (int i = min_i; i < max_i; i++) {
    progbar.update();
    GridCoords coords(_rnd_gen, Tissue::grid_size);
    DBG("infection: ", coords.str() + "\n");
    auto f = tissue.infect_epicell(coords).then([&num_infected](InfectionResult result) {
      if (result == InfectionResult::Success) num_infected++;
    });
    fchain = when_all(fchain, f);
    upcxx::progress();
  }
  fchain.wait();
  progbar.done();
  barrier();
  SLOG("Starting with ", reduce_one(num_infected, op_fast_add, 0).wait(), " infected cells\n");
  return num_infected;
}

int64_t generate_tcells(int num_tcells, Tissue &tissue) {
  BarrierTimer timer(__FILEFUNC__, false, true);
  // each rank generates a block of tcells
  int local_num_tcells = ceil((double)num_tcells / rank_n());
  int64_t min_tcell_id = rank_me() * local_num_tcells;
  int64_t max_tcell_id = (rank_me() + 1) * local_num_tcells;
  if (max_tcell_id > num_tcells) max_tcell_id = num_tcells;
  if (min_tcell_id > num_tcells) min_tcell_id = num_tcells;
  int64_t num_tcells_added = 0;
  barrier();
  DBG("Generating tcells from ", min_tcell_id, " to ", max_tcell_id, "\n");
  future<> fchain = upcxx::make_future<>();
  ProgressBar progbar(max_tcell_id - min_tcell_id, "Generating t-cells");
  for (int64_t tcell_id = min_tcell_id; tcell_id < max_tcell_id; tcell_id++) {
    progbar.update();
    GridCoords coords(_rnd_gen, Tissue::grid_size);
    DBG("tcell ", tcell_id, " starting at ", coords.str(), "\n");
    auto f = tissue.add_tcell(coords, tcell_id);
    fchain = upcxx::when_all(fchain, f);
    num_tcells_added++;
    upcxx::progress();
  }
  fchain.wait();
  progbar.done();
  barrier();
  // now set all tcells to be current
  tissue.update_tcells(/*check_switch=*/false);
  barrier();
  SLOG("Generated ", reduce_one(num_tcells_added, op_fast_add, 0).wait(), " t-cells\n");
  return num_tcells_added;
}

int64_t get_rnd_coord(int64_t x, int64_t max_x) {
  int64_t new_x = x + _rnd_gen->get(0, 3) - 1;
  if (new_x < 0) new_x = 0;
  if (new_x >= max_x) new_x = max_x - 1;
  return new_x;
}

void update_tcell(int time_step, Tissue &tissue, GridPoint *grid_point, TCell &tcell, int64_t &num_viral_kills, future<> &fchain) {
  _update_tcell_timer.start();
  if (grid_point->virus) {
    // kill the virus
    DBG(time_step, ": tcell ", tcell.id, " at ", grid_point->coords.str(), " killed virus\n");
    grid_point->virus = false;
    num_viral_kills++;
    // need to add here to ensure it gets to the next time step in the vector swap
    auto f = tissue.add_tcell(grid_point->coords, tcell.id);
    fchain = when_all(fchain, f);
  } else {
    // move to a random neighbor (or don't move)
    // GridCoords coords(get_rnd_coord(grid_point->coords.x, Tissue::grid_size.x),
    //                   get_rnd_coord(grid_point->coords.y, Tissue::grid_size.y),
    //                   get_rnd_coord(grid_point->coords.z, Tissue::grid_size.z));
    GridCoords coords = grid_point->neighbors[_rnd_gen->get(0, int64_t(grid_point->neighbors.size()))];
    DBG(time_step, ": tcell ", tcell.id, " at ", grid_point->coords.str(), " moving to ", coords.str(), "\n");
    tissue.add_tcell(coords, tcell.id);
  }
  _update_tcell_timer.stop();
}

void update_virus(int time_step, Tissue &tissue, GridPoint *grid_point, int64_t &num_infected, int64_t &num_dead_epicells,
                  future<> &fchain) {
  _update_virus_timer.start();
  assert(grid_point->epicell->status == EpiCellStatus::Incubating);
  grid_point->epicell->num_steps_infected++;
  if (grid_point->epicell->num_steps_infected > _options->incubation_period) {
    grid_point->epicell->status = EpiCellStatus::Dead;
    grid_point->virus = false;
    num_dead_epicells++;
    // spread virus to healthy neighbors with probability _options->spread_prob
    for (int64_t i = 0; i < grid_point->neighbors.size(); i++) {
      GridCoords coords = grid_point->neighbors[i];
      if (coords == grid_point->coords) continue;
      if (_rnd_gen->get_prob() >= _options->spread_prob) continue;
      auto grid_point_coords = grid_point->coords;
      auto f = tissue.infect_epicell(coords).then([&num_infected, time_step, coords,
                                                   grid_point_coords](InfectionResult result) {
        if (result == InfectionResult::Success) {
          DBG(time_step, ": virus at ", grid_point_coords.str(), " spread to ", coords.str(), "\n");
          num_infected++;
        } else {
          DBG(time_step, ": virus at ", grid_point_coords.str(), " could not spread to ", coords.str(), " (",
              (result == InfectionResult::NotHealthy ? " not healthy " : " already infected"), ")\n");
        }
      });
      fchain = when_all(fchain, f);
    }
    // for (int64_t x = 0; x < 3; x++) {
    //   auto vx = grid_point->coords.x + x - 1;
    //   if (vx < 0 || vx >= Tissue::grid_size.x) continue;
    //   for (int64_t y = 0; y < 3; y++) {
    //     auto vy = grid_point->coords.y + y - 1;
    //     if (vy < 0 || vy >= Tissue::grid_size.y) continue;
    //     for (int64_t z = 0; z < 3; z++) {
    //       auto vz = grid_point->coords.z + z - 1;
    //       if (vz < 0 || vz >= Tissue::grid_size.z) continue;
    //       GridCoords coords(vx, vy, vz);
    //       if (coords == grid_point->coords) continue;
    //       if (_rnd_gen->get_prob() >= _options->spread_prob) continue;
    //       auto grid_point_coords = grid_point->coords;
    //       auto f = tissue.infect_epicell(coords).then([&num_infected, time_step, coords,
    //                                                    grid_point_coords](InfectionResult result) {
    //         if (result == InfectionResult::Success) {
    //           DBG(time_step, ": virus at ", grid_point_coords.str(), " spread to ", coords.str(), "\n");
    //           num_infected++;
    //         } else {
    //           DBG(time_step, ": virus at ", grid_point_coords.str(), " could not spread to ", coords.str(), " (",
    //               (result == InfectionResult::NotHealthy ? " not healthy " : " already infected"), ")\n");
    //         }
    //       });
    //       fchain = when_all(fchain, f);
    //     }
    //   }
    // }
  }
  _update_virus_timer.stop();
}

void sample(int time_step, Tissue &tissue) {
  // each rank writes its own blocks into the file at the appropriate locations. We compute the location knowing that
  // each position takes exactly 4 characters (up to 3 numbers and a space), so we can work out how much space a block
  // takes up, in order to position the data correctly.
  // Each point is represented by a scalar char. From 1-127 are used to indicate viral load, with a color
  // of red, going from faint to strong. From -127 to -1 are used to indicate tcell count, with a color of blue,
  // going from faint to strong. The last color (0) is used to indicate a dead epicell (how do we make this a separate color?)
  _sample_timer.start();
  // each grid point takes up a single char
  size_t tot_sz = tissue.get_num_grid_points() * 1;
  string fname = "samples/sample_" + to_string(time_step) + ".vtk";
  int x_dim = _options->dimensions[0];
  int y_dim = _options->dimensions[1];
  int z_dim = _options->dimensions[2];
  ostringstream header_oss;
  header_oss << "# vtk DataFile Version 4.2\n"
             << "SimCov sample " << basename(_options->output_dir.c_str()) << time_step << "\n"
             //<< "ASCII\n"
             << "BINARY\n"
             << "DATASET STRUCTURED_POINTS\n"
             // we add one in each dimension because these are for drawing the visualization points, and our visualization
             // entities are cells
             << "DIMENSIONS " << (x_dim + 1) << " " << (y_dim + 1) << " " << (z_dim + 1) << "\n"
             << "SPACING 1 1 1\n"
             << "ORIGIN 0 0 0\n"
             << "CELL_DATA " << (x_dim * y_dim * z_dim) << "\n"
             << "SCALARS virus_or_tcell unsigned_char 1\n"
             << "LOOKUP_TABLE default\n";
  if (!rank_me()) {
    // all ranks have the same size since we zero pad ASCII chars to ensure we always write out 3 chars per point
    tot_sz += header_oss.str().size();
    // rank 0 creates the file and truncates it to the correct length
    auto fileno = open(fname.c_str(), O_WRONLY | O_CREAT, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
    if (fileno == -1) SDIE("Cannot create file ", fname, ": ", strerror(errno), "\n");
    if (ftruncate(fileno, tot_sz) == -1) DIE("Could not truncate ", fname, " to ", tot_sz, " bytes\n");
    close(fileno);
    DBG("Truncated sample file ", fname, " to ", tot_sz, " bytes\n");
  }
  upcxx::barrier();
  // wait until rank 0 has finished setting up the file
  auto [bytes_written, grid_points_written] = tissue.dump_blocks(fname, header_oss.str());
  upcxx::barrier();
  auto tot_bytes_written = upcxx::reduce_one(bytes_written, upcxx::op_fast_add, 0).wait();
  auto tot_grid_points_written = upcxx::reduce_one(grid_points_written, upcxx::op_fast_add, 0).wait();
  if (!rank_me()) assert(tot_grid_points_written == tissue.get_num_grid_points());
  SLOG_VERBOSE("Successfully wrote ", tot_grid_points_written, " grid points, with size ", get_size_str(tot_bytes_written), " to ",
               fname, "\n");
  _sample_timer.stop();
}

void run_sim(Tissue &tissue, int64_t &num_tcells, int64_t &tot_num_infected) {
  BarrierTimer timer(__FILEFUNC__, false, true);
  // this is a trivial test case:
  // when a virus first arrives in an epicell, it spends a period of time there
  // after the period has elapsed, the epicell dies and the virus spreads in all directions
  // the virus cannot spread to an already dead cell
  // a t-cell moves to a new point at random
  // when a t-cell arrives at a point where there is a virus, the virus goes to 0 and the epicell dies
  int64_t num_dead_epicells = 0;
  int64_t tot_num_viral_kills = 0;
  double time_step_ticks = (double)_options->num_iters / 20;
  auto start_t = chrono::high_resolution_clock::now();
  auto curr_t = start_t;
  for (int time_step = 0; time_step < _options->num_iters; time_step++) {
    int64_t num_infected = 0;
    int64_t num_viral_kills = 0;
    future<> fchain = make_future<>();
    // iterate through all local grid points
    // FIXME: this should be iteration through the local active grid points only
    for (auto grid_point = tissue.get_first_local_grid_point(); grid_point; grid_point = tissue.get_next_local_grid_point()) {
      upcxx::progress();
      // the tcells are moved (added to the new list, but only cleared out at the end of all updates)
      if (grid_point->tcells) {
        for (auto &tcell : *grid_point->tcells) {
          update_tcell(time_step, tissue, grid_point, tcell, num_viral_kills, fchain);
        }
      }
      if (grid_point->virus) {
        if (!grid_point->epicell->num_steps_infected) {
          // we will process this next round
          grid_point->epicell->num_steps_infected = 1;
          continue;
        }
        update_virus(time_step, tissue, grid_point, num_infected, num_dead_epicells, fchain);
      }
    }
    fchain.wait();
    barrier();
    num_tcells = tissue.update_tcells();
    barrier();
    tot_num_infected += num_infected;
    tot_num_viral_kills += num_viral_kills;
    // print every 5% of the iterations
    if (time_step >= time_step_ticks || time_step == _options->num_iters - 1) {
      chrono::duration<double> t_elapsed = chrono::high_resolution_clock::now() - curr_t;
      curr_t = chrono::high_resolution_clock::now();
      // memory doesn't really change so don't report it every iteration
      SLOG(setw(5), left, time_step,
           " [", get_current_time(true), " ", setprecision(2), fixed, setw(5), right, t_elapsed.count(), "s",
           /*setw(9), left, get_size_str(get_free_mem()),*/ "]: ",
           " infections ", reduce_one(tot_num_infected, op_fast_add, 0).wait(), "+",
           reduce_one(num_infected, op_fast_add, 0).wait(),
           ", dead epicells ", perc_str(reduce_one(num_dead_epicells, op_fast_add, 0).wait(), tissue.get_num_grid_points()),
           ", viral kills ", reduce_one(tot_num_viral_kills, op_fast_add, 0).wait(), "+",
           reduce_one(num_viral_kills, op_fast_add, 0).wait(),
           ", t-cells: ", reduce_one(num_tcells, op_fast_add, 0).wait(), "\n");
      time_step_ticks += (double)_options->num_iters / 20;
    }
    // sample
    if (time_step % _options->sample_period == 0) {
      sample(time_step, tissue);
    }
  }
  _update_tcell_timer.done_all();
  _update_virus_timer.done_all();
  _sample_timer.done_all();
  chrono::duration<double> t_elapsed = chrono::high_resolution_clock::now() - start_t;
  SLOG("Finished ", _options->num_iters, " time steps in ", setprecision(4), fixed, t_elapsed.count(), " s (",
       (double)t_elapsed.count() / _options->num_iters, " s per step)\n");
}

int main(int argc, char **argv) {
  upcxx::init();
  auto start_t = chrono::high_resolution_clock::now();
  _options = make_shared<Options>();
  if (!_options->load(argc, argv)) return 0;
  ProgressBar::SHOW_PROGRESS = _options->show_progress;

  _rnd_gen = make_shared<Random>(_options->rnd_seed + rank_me());

  if (pin_thread(getpid(), local_team().rank_me()) == -1) WARN("Could not pin process ", getpid(), " to core ", rank_me());
  else SLOG_VERBOSE("Pinned processes, with process 0 (pid ", getpid(), ") pinned to core ", local_team().rank_me(), "\n");

  MemoryTrackerThread memory_tracker;
  memory_tracker.start();
  auto start_free_mem = get_free_mem();
  SLOG(KBLUE, "Starting with ", get_size_str(start_free_mem), " free on node 0", KNORM, "\n");
  Tissue tissue;
  tissue.construct({_options->dimensions[0], _options->dimensions[1], _options->dimensions[2]});
  int64_t num_infected = initial_infection(_options->num_infections, tissue);
  int64_t num_tcells = generate_tcells(_options->num_tcells, tissue);
  SLOG(KBLUE, "Memory used on node 0 after initialization is  ", get_size_str(start_free_mem - get_free_mem()), KNORM, "\n");

  run_sim(tissue, num_tcells, num_infected);

  memory_tracker.stop();
  chrono::duration<double> t_elapsed = chrono::high_resolution_clock::now() - start_t;
  SLOG("Finished in ", setprecision(2), fixed, t_elapsed.count(), " s at ", get_current_time(),
       " for SimCov version ", SIMCOV_VERSION, "\n");
  barrier();
  upcxx::finalize();
  return 0;
}
