// SimCov
//
// Steven Hofmeyr, LBNL May 2020

#include <fcntl.h>
#include <math.h>
#include <stdarg.h>
#include <unistd.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <chrono>
#include <upcxx/upcxx.hpp>

using namespace std;

#include "options.hpp"
#include "tissue.hpp"
#include "upcxx_utils.hpp"
#include "utils.hpp"

using namespace upcxx;
using namespace upcxx_utils;

#define NOW chrono::high_resolution_clock::now

struct SimStats {
  int64_t num_infected = 0;
  int64_t num_dead_epicells = 0;
  int64_t num_tcells_in_vasculature = 0;
  int64_t num_tcells_in_tissue = 0;
  int64_t tot_num_viral_kills = 0;

  void clear() {
    num_infected = 0;
    num_dead_epicells = 0;
    num_tcells_in_vasculature = 0;
    num_tcells_in_tissue = 0;
  }

  string to_str(int64_t num_grid_points) {
    ostringstream oss;
    oss << " infections " << upcxx::reduce_one(num_infected, upcxx::op_fast_add, 0).wait()
        << " dead epicells "
        << perc_str(upcxx::reduce_one(num_dead_epicells, upcxx::op_fast_add, 0).wait(),
                    num_grid_points)
        << " viral kills " << upcxx::reduce_one(tot_num_viral_kills, upcxx::op_fast_add, 0).wait()
        << " t-cells in vasculature "
        << upcxx::reduce_one(num_tcells_in_vasculature, upcxx::op_fast_add, 0).wait()
        << " t-cells in tissue "
        << upcxx::reduce_one(num_tcells_in_tissue, upcxx::op_fast_add, 0).wait();
    return oss.str();
  }
};

ofstream _logstream;
bool _verbose = false;
SimStats _sim_stats;
shared_ptr<Options> _options;

IntermittentTimer _generate_tcell_timer(__FILENAME__ + string(":") + "generate_tcells");
IntermittentTimer _update_tcell_timer(__FILENAME__ + string(":") + "update_tcell");
IntermittentTimer _update_epicell_timer(__FILENAME__ + string(":") + "update_epicell");
IntermittentTimer _update_concentration_timer(__FILENAME__ + string(":") + "update_concentration");
IntermittentTimer _sample_timer(__FILENAME__ + string(":") + "sample");
IntermittentTimer _finish_round_timer(__FILENAME__ + string(":") + "finish_round");

void initial_infection(Tissue &tissue) {
  BarrierTimer timer(__FILEFUNC__);
  // each rank generates a block of infected cells
  // for large dimensions, the chance of repeat sampling for a few points is very small
  int local_num_infections = _options->num_infections / rank_n();
  int remaining_infections = _options->num_infections - local_num_infections * rank_n();
  if (rank_me() < remaining_infections) local_num_infections++;
  DBG("Infecting ", local_num_infections, " epicells\n");
  ProgressBar progbar(local_num_infections, "Setting initial infections");
  for (int i = 0; i < local_num_infections; i++) {
    progbar.update();
    GridCoords coords(_rnd_gen, Tissue::grid_size);
    DBG("infection: ", coords.str() + "\n");
    tissue.inc_incoming_virus(coords, 1.0);
    upcxx::progress();
  }
  progbar.done();
  barrier();
  SLOG("Initially infected ", reduce_one(local_num_infections, op_fast_add, 0).wait(),
       " epicells\n");
  int num_infections_found = 0;
  for (auto grid_point = tissue.get_first_active_grid_point(); grid_point;
       grid_point = tissue.get_next_active_grid_point()) {
    if (grid_point->incoming_virus > 0) {
      grid_point->virus = 1.0;
      grid_point->incoming_virus = 0;
      grid_point->epicell->infect();
      num_infections_found++;
    }
  }
  barrier();
  tissue.clear_active();
  barrier();
  int tot_num_infections_found = reduce_one(num_infections_found, op_fast_add, 0).wait();
  if (!rank_me() && tot_num_infections_found != _options->num_infections)
    WARN("Generated fewer initial infections that expected, ", tot_num_infections_found, " < ",
         _options->num_infections);
}

void generate_tcells(Tissue &tissue) {
  _generate_tcell_timer.start();
  // each rank could generate some t-cells
  int local_num_tcells = _options->tcell_generation_rate / rank_n();
  int remaining_tcells = _options->tcell_generation_rate - local_num_tcells * rank_n();
  if (rank_me() < remaining_tcells) local_num_tcells++;
  if (local_num_tcells) {
    for (int i = 0; i < local_num_tcells; i++) {
      GridCoords coords(_rnd_gen, Tissue::grid_size);
      string tcell_id = to_string(rank_me()) + "-" + to_string(tissue.tcells_generated);
      tissue.tcells_generated++;
      // the initial position is actually meaningless
      tissue.add_tcell(coords, TCell(tcell_id, coords));
      upcxx::progress();
    }
  }
  _generate_tcell_timer.stop();
  barrier();
}

int64_t get_rnd_coord(int64_t x, int64_t max_x) {
  int64_t new_x = x + _rnd_gen->get(0, 3) - 1;
  if (new_x < 0) new_x = 0;
  if (new_x >= max_x) new_x = max_x - 1;
  return new_x;
}

void update_tcell(int time_step, Tissue &tissue, GridPoint *grid_point, TCell &tcell) {
  _update_tcell_timer.start();
  if (grid_point->icytokine > 0) {
    tcell.in_vasculature = false;
    tcell.prev_coords = grid_point->coords;
  }
  if (tcell.in_vasculature) {
    tcell.vascular_period--;
    if (tcell.vascular_period > 0) {
      // tcell is still alive - moves to any random location in the grid
      tissue.add_tcell({_rnd_gen, Tissue::grid_size}, tcell);
    }
  } else {
    tcell.tissue_period--;
    // tcell is in the tissue
    if (tcell.tissue_period > 0) {
      // still alive
      if (grid_point->epicell->status == EpiCellStatus::EXPRESSING) {
        // induce apoptosis
        DBG(time_step, ": tcell ", tcell.id, " at ", grid_point->coords.str(),
            " induced apoptosis\n");
        grid_point->epicell->induce_apoptosis();
        // add here to ensure it gets to the next time step in the vector swap, in the same spot
        tcell.prev_coords = grid_point->coords;
        tissue.add_tcell(grid_point->coords, tcell);
      } else {
        // no apoptosis to induce, will move, either following gradient or at random
        GridCoords selected_coords;
        double highest_chemokine = 0;
        for (auto &nb_coords : grid_point->neighbors) {
          if (nb_coords == grid_point->coords) continue;
          double chemokine = tissue.get_chemokine(nb_coords);
          if (chemokine > highest_chemokine) {
            highest_chemokine = chemokine;
            selected_coords = nb_coords;
          }
        }
        if (highest_chemokine == 0) {
          // no chemokine found, move at random, except back to previous coords
          do {
            auto rnd_nb_i = _rnd_gen->get(0, (int64_t)grid_point->neighbors.size());
            selected_coords = grid_point->neighbors[rnd_nb_i];
          } while (selected_coords == tcell.prev_coords);
        } else {
          DBG(time_step, ": highest nb chemokine for tcell ", tcell.id, " at ",
              grid_point->coords.str(), " is at ", selected_coords.str(), " with ",
              highest_chemokine, "\n");
        }
        DBG(time_step, ": tcell ", tcell.id, " at ", grid_point->coords.str(), " moving to ",
            selected_coords.str(), "\n");
        tcell.prev_coords = grid_point->coords;
        tissue.add_tcell(selected_coords, tcell);
      }
    }
  }
  _update_tcell_timer.stop();
}

void update_epicell(int time_step, Tissue &tissue, GridPoint *grid_point) {
  if (grid_point->epicell->status == EpiCellStatus::DEAD) return;

  _update_epicell_timer.start();
  bool spreading = false;
  switch (grid_point->epicell->status) {
    case EpiCellStatus::INCUBATING:
      grid_point->epicell->transition_to_expressing();
      break;
    case EpiCellStatus::APOPTOTIC:
      if (grid_point->epicell->apoptosis_death()) {
        grid_point->virus = 0;
        _sim_stats.tot_num_viral_kills++;
        break;
      }
      spreading = true;
      break;
    case EpiCellStatus::EXPRESSING:
      if (grid_point->epicell->infection_death()) {
        grid_point->virus = 0;
        break;
      }
      spreading = true;
      break;
    case EpiCellStatus::HEALTHY:
      if (grid_point->virus > 0) {
        double infection_prob = grid_point->virus * _options->virus_infection_prob;
        if (_rnd_gen->trial_success(infection_prob)) grid_point->epicell->infect();
      }
      break;
    default: break;
  }
  if (spreading) {
    // always gives out max concentrations every time step
    grid_point->chemokine = 1.0;
    grid_point->incoming_chemokine = 1.0;
    grid_point->icytokine = 1.0;
    grid_point->incoming_icytokine = 1.0;
    grid_point->virus = 1.0;
    grid_point->incoming_virus = 1.0;
  }
  _update_epicell_timer.stop();
}

void update_concentration(int time_step, GridPoint *grid_point,
                          double &concentration, double decay_rate, double diffusion_coef,
                          Tissue &tissue, void (Tissue::*inc_incoming)(GridCoords, double)) {
  _update_concentration_timer.start();
  if (concentration > 0) {
    if (grid_point->epicell->status != EpiCellStatus::EXPRESSING) {
      concentration -= decay_rate;
      if (concentration < 0) concentration = 0;
    }
    if (concentration > 0) {
      // diffuses, reducing concentration at source
      double diffusion_amount = concentration * diffusion_coef;
      concentration -= diffusion_amount;
      // the amount diffusing spreads uniformly to all neighbors
      double amount_to_nb = diffusion_amount / grid_point->neighbors.size();
      for (auto &nb_coords : grid_point->neighbors) {
        if (nb_coords == grid_point->coords) continue;
        (tissue.*inc_incoming)(nb_coords, amount_to_nb);
      }
    }
  }
  _update_concentration_timer.stop();
}

void finish_round(Tissue &tissue, int time_step) {
  _finish_round_timer.start();
  barrier();
  _sim_stats.clear();
  // FIXME: go through active list
  for (auto grid_point = tissue.get_first_active_grid_point(); grid_point;
       grid_point = tissue.get_next_active_grid_point()) {
    if (grid_point->tcells) {
      grid_point->switch_tcells_vector();
      for (auto &tcell : *grid_point->tcells) {
        if (tcell.in_vasculature) _sim_stats.num_tcells_in_vasculature++;
        else _sim_stats.num_tcells_in_tissue++;
      }
    }
    if (grid_point->incoming_virus) {
      grid_point->virus += grid_point->incoming_virus;
      if (grid_point->virus > 1) grid_point->virus = 1;
      grid_point->incoming_virus = 0;
    }
    if (grid_point->incoming_chemokine) {
      grid_point->chemokine += grid_point->incoming_chemokine;
      if (grid_point->chemokine > 1) grid_point->chemokine = 1;
      grid_point->incoming_chemokine = 0;
    }
    if (grid_point->incoming_icytokine) {
      grid_point->icytokine += grid_point->incoming_icytokine;
      if (grid_point->icytokine > 1) grid_point->icytokine = 1;
      grid_point->incoming_icytokine = 0;
    }
    switch (grid_point->epicell->status) {
      case EpiCellStatus::INCUBATING:
      case EpiCellStatus::EXPRESSING:
      case EpiCellStatus::APOPTOTIC: _sim_stats.num_infected++; break;
      case EpiCellStatus::DEAD: _sim_stats.num_dead_epicells++; break;
      case EpiCellStatus::HEALTHY: break;
    }
  }
  barrier();
  tissue.clear_active();
  barrier();
  _finish_round_timer.stop();
}

void write_paraview_state() {
  /*
  if (rank_me() > 0) return;
  ofstream outf("sample_state.py")
  outf << "from paraview.simple import *\n"
  */
}

void sample(int time_step, Tissue &tissue, ViewObject view_object) {
  // each rank writes its own blocks into the file at the appropriate locations. We compute the
  // location knowing that each position takes exactly 4 characters (up to 3 numbers and a space),
  // so we can work out how much space a block takes up, in order to position the data correctly.
  // Each point is represented by a scalar char. From 1-127 are used to indicate viral load, with a
  // color of red, going from faint to strong. From -127 to -1 are used to indicate tcell count,
  // with a color of blue, going from faint to strong. The last color (0) is used to indicate a dead
  // epicell (how do we make this a separate color?)
  _sample_timer.start();
  // each grid point takes up a single char
  size_t tot_sz = tissue.get_num_grid_points();
  string fname = "samples/sample_" + view_object_str(view_object) + "_" + to_string(time_step) +
                 ".vtk";
  int x_dim = _options->dimensions[0];
  int y_dim = _options->dimensions[1];
  int z_dim = _options->dimensions[2];
  ostringstream header_oss;
  header_oss << "# vtk DataFile Version 4.2\n"
             << "SimCov sample " << basename(_options->output_dir.c_str()) << time_step
             << "\n"
             //<< "ASCII\n"
             << "BINARY\n"
             << "DATASET STRUCTURED_POINTS\n"
             // we add one in each dimension because these are for drawing the
             // visualization points, and our visualization entities are cells
             << "DIMENSIONS " << (x_dim + 1) << " " << (y_dim + 1) << " " << (z_dim + 1)
             << "\n"
             // each cell is 5 microns
             << "SPACING 5 5 5\n"
             << "ORIGIN 0 0 0\n"
             << "CELL_DATA " << (x_dim * y_dim * z_dim) << "\n"
             << "SCALARS ";
  switch (view_object) {
    case ViewObject::VIRUS: header_oss << "virus"; break;
    case ViewObject::TCELL_TISSUE: header_oss << "t-cell-tissue"; break;
    case ViewObject::EPICELL: header_oss << "epicell"; break;
    case ViewObject::ICYTOKINE: header_oss << "icytokine"; break;
    case ViewObject::CHEMOKINE: header_oss << "chemokine"; break;
    default: SDIE("unknown view object");
  }
  header_oss << " unsigned_char\n"
             << "LOOKUP_TABLE default\n";
  if (!rank_me()) {
    tot_sz += header_oss.str().size();
    // rank 0 creates the file and truncates it to the correct length
    auto fileno = open(fname.c_str(), O_WRONLY | O_CREAT, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
    if (fileno == -1) SDIE("Cannot create file ", fname, ": ", strerror(errno), "\n");
    if (ftruncate(fileno, tot_sz) == -1)
      DIE("Could not truncate ", fname, " to ", tot_sz, " bytes\n");
    close(fileno);
    DBG("Truncated sample file ", fname, " to ", tot_sz, " bytes\n");
  }
  upcxx::barrier();
  // wait until rank 0 has finished setting up the file
  auto [bytes_written, grid_points_written] = tissue.dump_blocks(fname, header_oss.str(),
                                                                 view_object);
  assert(bytes_written == tot_sz);
  upcxx::barrier();
  auto tot_bytes_written = reduce_one(bytes_written, op_fast_add, 0).wait();
  auto tot_grid_points_written = reduce_one(grid_points_written, op_fast_add, 0).wait();
  if (!rank_me()) assert(tot_grid_points_written == tissue.get_num_grid_points());
  _sample_timer.stop();
}

void run_sim(Tissue &tissue) {
  BarrierTimer timer(__FILEFUNC__);
  auto start_t = NOW();
  auto curr_t = start_t;
  auto five_perc = _options->num_timesteps / 20;
  int tick = 0;
  for (int time_step = 0; time_step < _options->num_timesteps; time_step++) {
    if (time_step > _options->tcell_initial_delay) generate_tcells(tissue);
    // iterate through all local grid points
    // FIXME: this should be iteration through the local active grid points only
    for (auto grid_point = tissue.get_first_local_grid_point(); grid_point;
         grid_point = tissue.get_next_local_grid_point()) {
      upcxx::progress();
      // the tcells are moved (added to the new list, but only cleared out at the end of all
      // updates)
      if (grid_point->tcells && grid_point->tcells->size()) {
        for (auto &tcell : *grid_point->tcells) {
          update_tcell(time_step, tissue, grid_point, tcell);
        }
      }
      update_epicell(time_step, tissue, grid_point);
      update_concentration(time_step, grid_point, grid_point->chemokine,
                           _options->chemokine_decay_rate, _options->chemokine_diffusion_coef,
                           tissue, &Tissue::inc_incoming_chemokines);
      update_concentration(time_step, grid_point, grid_point->icytokine,
                           _options->icytokine_decay_rate, _options->icytokine_diffusion_coef,
                           tissue, &Tissue::inc_incoming_icytokines);
      update_concentration(time_step, grid_point, grid_point->virus,
                           _options->virus_decay_rate, _options->virus_diffusion_coef,
                           tissue, &Tissue::inc_incoming_virus);
    }
    finish_round(tissue, time_step);
    if (time_step % _options->sample_period == 0 || time_step == _options->num_timesteps - 1) {
      chrono::duration<double> t_elapsed = NOW() - curr_t;
      if (_options->verbose) curr_t = NOW();
      SLOG_VERBOSE(setw(5), left, time_step, " [", get_current_time(true), " ", setprecision(2),
                   fixed, setw(5), right, t_elapsed.count(),
                   " s]: ", _sim_stats.to_str(tissue.get_num_grid_points()), "\n");
    }
    if (!_options->verbose && tick % five_perc == 0) {
      chrono::duration<double> t_elapsed = NOW() - curr_t;
      curr_t = NOW();
      SLOG(setw(5), left, time_step, " [", get_current_time(true), " ", setprecision(2), fixed,
           setw(5), right, t_elapsed.count(),
           " s]: ", _sim_stats.to_str(tissue.get_num_grid_points()), "\n");
    }
    tick++;
    // sample
    if (time_step % _options->sample_period == 0) {
      sample(time_step, tissue, ViewObject::TCELL_TISSUE);
      sample(time_step, tissue, ViewObject::VIRUS);
      sample(time_step, tissue, ViewObject::EPICELL);
      sample(time_step, tissue, ViewObject::ICYTOKINE);
      sample(time_step, tissue, ViewObject::CHEMOKINE);
    }
  }
  _generate_tcell_timer.done_all();
  _update_tcell_timer.done_all();
  _update_epicell_timer.done_all();
  _update_concentration_timer.done_all();
  _sample_timer.done_all();
  _finish_round_timer.done_all();

  chrono::duration<double> t_elapsed = NOW() - start_t;
  SLOG("Finished ", _options->num_timesteps, " time steps in ", setprecision(4), fixed,
       t_elapsed.count(), " s (", (double)t_elapsed.count() / _options->num_timesteps,
       " s per step)\n");
  // write out paraview state file for easy preformatted viewing
  write_paraview_state();
}

int main(int argc, char **argv) {
  upcxx::init();
  auto start_t = NOW();
  _options = make_shared<Options>();
  if (!_options->load(argc, argv)) return 0;
  ProgressBar::SHOW_PROGRESS = _options->show_progress;

  _rnd_gen = make_shared<Random>(_options->rnd_seed + rank_me());

  if (pin_thread(getpid(), local_team().rank_me()) == -1)
    WARN("Could not pin process ", getpid(), " to core ", rank_me());
  else
    SLOG_VERBOSE("Pinned processes, with process 0 (pid ", getpid(), ") pinned to core ",
                 local_team().rank_me(), "\n");

  MemoryTrackerThread memory_tracker;
  memory_tracker.start();
  auto start_free_mem = get_free_mem();
  SLOG(KBLUE, "Starting with ", get_size_str(start_free_mem), " free on node 0", KNORM, "\n");
  Tissue tissue;
  tissue.construct({_options->dimensions[0], _options->dimensions[1], _options->dimensions[2]});
  initial_infection(tissue);
  SLOG(KBLUE, "Memory used on node 0 after initialization is  ",
       get_size_str(start_free_mem - get_free_mem()), KNORM, "\n");

  run_sim(tissue);

  memory_tracker.stop();
  chrono::duration<double> t_elapsed = NOW() - start_t;
  SLOG("Finished in ", setprecision(2), fixed, t_elapsed.count(), " s at ", get_current_time(),
       " for SimCov version ", SIMCOV_VERSION, "\n");
  barrier();
  upcxx::finalize();
  return 0;
}
