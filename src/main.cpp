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

int64_t initial_infection(int num_infections, Tissue &tissue) {
  BarrierTimer timer(__FILEFUNC__, false, true);
  Random rnd_gen;
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
  ProgressBar progbar(max_i - min_i, "Setting initial infections");
  for (int i = min_i; i < max_i; i++) {
    progbar.update();
    GridCoords coords(rnd_gen, Tissue::grid_size);
    DBG("infection: ", coords.str() + "\n");
    if (tissue.infect_epicell(coords) == InfectionResult::Success) num_infected++;
    upcxx::progress();
  }
  progbar.done();
  barrier();
  SLOG("Starting with ", reduce_one(num_infected, op_fast_add, 0).wait(), " infected cells\n");
  return num_infected;
}

int64_t generate_tcells(int num_tcells, Tissue &tissue) {
  BarrierTimer timer(__FILEFUNC__, false, true);
  Random rnd_gen;
  // each rank generates a block of tcells
  int local_num_tcells = ceil((double)num_tcells / rank_n());
  int64_t min_tcell_id = rank_me() * local_num_tcells;
  int64_t max_tcell_id = (rank_me() + 1) * local_num_tcells;
  if (max_tcell_id > num_tcells) max_tcell_id = num_tcells;
  if (min_tcell_id > num_tcells) min_tcell_id = num_tcells;
  int64_t num_tcells_added = 0;
  barrier();
  DBG("Generating tcells from ", min_tcell_id, " to ", max_tcell_id, "\n");
  ProgressBar progbar(max_tcell_id - min_tcell_id, "Generating t-cells");
  for (int64_t tcell_id = min_tcell_id; tcell_id < max_tcell_id; tcell_id++) {
    progbar.update();
    GridCoords coords(rnd_gen, Tissue::grid_size);
    DBG("tcell: ", coords.str(), "\n");
    tissue.add_tcell(coords, tcell_id, 0);
    num_tcells_added++;
    upcxx::progress();
  }
  progbar.done();
  barrier();
  SLOG("Generated ", reduce_one(num_tcells_added, op_fast_add, 0).wait(), " t-cells\n");
  return num_tcells_added;
}

int64_t get_rnd_coord(Random &rnd_gen, int64_t x, int64_t max_x) {
  int64_t new_x = x + rnd_gen.get(0, 3) - 1;
  if (new_x < 0) new_x = 0;
  if (new_x >= max_x) new_x = max_x - 1;
  return new_x;
}

void run_sim(Tissue &tissue, int64_t &num_tcells, int64_t &num_infected, shared_ptr<Options> options) {
  BarrierTimer timer(__FILEFUNC__, false, true);
  // this is a trivial test case:
  // when a virus first arrives in an epicell, it spends a period of time there
  // after the period has elapsed, the epicell dies and the virus spreads in all directions
  // the virus cannot spread to an already dead cell
  // a t-cell moves to a new point at random
  // when a t-cell arrives at a point where there is a virus, the virus goes to 0 and the epicell dies
  Random rnd_gen;
  int64_t num_viral_kills = 0;
  int64_t num_dead_epicells = 0;
  double step_ticks = (double)options->num_iters / 20;
  for (int step = 0; step < options->num_iters; step++) {
    num_infected = 0;
    // iterate through all local grid points
    // FIXME: this should be iteration through the local active grid points only
    for (auto grid_point = tissue.get_first_local_grid_point(); grid_point; grid_point = tissue.get_next_local_grid_point()) {
      // first iterate through all the tcells, killing viruses and moving the tcells
      // note that we move the tcell to the destination, but don't remove from the target because we don't want to modify the
      // tcell vectors while they are being iterated
      // And because t-cells can be added during this loop, we iterate explicitly through the cell indexes, so we don't process
      // the added cells
      int num_tcells = grid_point->tcells.size();
      for (int i = 0; i < num_tcells; i++) {
        TCell &tcell = grid_point->tcells[i];
        // tcell is not processed yet
        if (tcell.step != step) {
          tcell.step = step;
          if (grid_point->virus) {
            // kill the virus
            DBG(step, ": tcell at ", grid_point->coords.str(), " killed virus\n");
            grid_point->virus = false;
            num_viral_kills++;
          } else {
            // move to a random neighbor
            GridCoords coords(get_rnd_coord(rnd_gen, grid_point->coords.x, Tissue::grid_size.x),
                              get_rnd_coord(rnd_gen, grid_point->coords.y, Tissue::grid_size.y),
                              get_rnd_coord(rnd_gen, grid_point->coords.z, Tissue::grid_size.z));
            if (coords != grid_point->coords) {
              DBG(step, ": tcell at ", grid_point->coords.str(), " moving to ", coords.str(), "\n");
              tcell.coords = coords;
              tissue.add_tcell(coords, tcell.id, step);
            }
          }
          upcxx::progress();
        }
      }
      if (grid_point->virus) {
        assert(grid_point->epicell->status == EpiCellStatus::Incubating);
        grid_point->epicell->infection_steps++;
        if (grid_point->epicell->infection_steps > options->incubation_period) {
          grid_point->epicell->status == EpiCellStatus::Dead;
          grid_point->virus = 0;
          num_dead_epicells++;
          // spread virus to healthy neighbors
          for (int64_t x = 0; x < 3; x++) {
            auto vx = grid_point->coords.x + x - 1;
            if (vx < 0 || vx >= Tissue::grid_size.x) continue;
            for (int64_t y = 0; y < 3; y++) {
              auto vy = grid_point->coords.y + y - 1;
              if (vy < 0 || vy >= Tissue::grid_size.y) continue;
              for (int64_t z = 0; z < 3; z++) {
                auto vz = grid_point->coords.z + z - 1;
                if (vz < 0 || vz >= Tissue::grid_size.z) continue;
                GridCoords coords(vx, vy, vz);
                if (coords == grid_point->coords) continue;
                auto res = tissue.infect_epicell(coords);
                if (res == InfectionResult::Success) {
                  DBG(step, ": virus at ", grid_point->coords.str(), " spread to ", coords.str(), "\n");
                  num_infected++;
                } else {
                  DBG(step, ": virus at ", grid_point->coords.str(), " could not spread to ", coords.str(),
                      " (", (res == InfectionResult::NotHealthy ? " not healthy " : " already infected"), ")\n");
                }
              }
            }
          }
        }
      }
      upcxx::progress();
    }
    barrier();
    auto all_num_tcells = reduce_one(num_tcells, op_fast_add, 0).wait();
    int64_t num_found_tcells = 0, num_deleted_tcells = 0;
    // when every rank has finished processing, remove tcells from old grid points
    for (auto grid_point = tissue.get_first_local_grid_point(); grid_point; grid_point = tissue.get_next_local_grid_point()) {
      for (auto it = grid_point->tcells.begin(); it != grid_point->tcells.end(); ) {
        auto tcell = make_shared<TCell>(*it);
        if (tcell->coords != grid_point->coords) {
          it = grid_point->tcells.erase(it);
          num_deleted_tcells++;
        } else {
          ++it;
          num_found_tcells++;
        }
      }
    }
    DBG("Deleted ", num_deleted_tcells, " tcells that have moved\n");
    barrier();
#ifdef DEBUG
    auto all_num_found_tcells = reduce_one(num_found_tcells, op_fast_add, 0).wait();
    if (!rank_me() && all_num_found_tcells != all_num_tcells)
      SDIE("Found ", all_num_found_tcells, " but expected ", all_num_tcells);
#endif
    // print every 5% of the iterations
    if (step >= step_ticks) {
      SLOG(step, ": ", "infections ", perc_str(reduce_one(num_infected, op_fast_add, 0).wait(), tissue.get_num_grid_points()),
           ", dead epicells ", perc_str(reduce_one(num_dead_epicells, op_fast_add, 0).wait(), tissue.get_num_grid_points()),
           ", viral kills ", reduce_one(num_viral_kills, op_fast_add, 0).wait(),
           ", t-cells: ", all_num_tcells, "\n");
      step_ticks += (double)options->num_iters / 20;
    }
  }
}

int main(int argc, char **argv) {
  upcxx::init();
  auto start_t = chrono::high_resolution_clock::now();
  auto options = make_shared<Options>();
  if (!options->load(argc, argv)) return 0;
  ProgressBar::SHOW_PROGRESS = options->show_progress;

  if (pin_thread(getpid(), local_team().rank_me()) == -1) WARN("Could not pin process ", getpid(), " to core ", rank_me());
  else SLOG_VERBOSE("Pinned processes, with process 0 (pid ", getpid(), ") pinned to core ", local_team().rank_me(), "\n");

  MemoryTrackerThread memory_tracker;
  memory_tracker.start();
  SLOG(KBLUE, "Starting with ", get_size_str(get_free_mem()), " free on node 0", KNORM, "\n");

  Tissue tissue;
  tissue.construct({options->dimensions[0], options->dimensions[1], options->dimensions[2]});
  int64_t num_infected = initial_infection(options->num_infections, tissue);
  int64_t num_tcells = generate_tcells(options->num_tcells, tissue);

  run_sim(tissue, num_tcells, num_infected, options);

  memory_tracker.stop();
  chrono::duration<double> t_elapsed = chrono::high_resolution_clock::now() - start_t;
  SLOG("Finished in ", setprecision(2), fixed, t_elapsed.count(), " s at ", get_current_time(),
       " for SimCov version ", SIMCOV_VERSION, "\n");
  barrier();
  upcxx::finalize();
  return 0;
}
