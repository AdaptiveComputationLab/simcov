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

class SimStats {

 private:
  ofstream log_file;

 public:
  int64_t incubating = 0;
  int64_t expressing = 0;
  int64_t apoptotic = 0;
  int64_t dead = 0;
  int64_t tcells_vasculature = 0;
  int64_t tcells_tissue = 0;
  double chemokines = 0;
  double icytokines = 0;
  double virus = 0;

  void init() {
    if (!rank_me()) {
      log_file.open("simcov.stats");
      log_file << "# " << header() << endl;
    }
  }

  string header() {
    ostringstream oss;
    oss << "incb\t"
        << "expr\t"
        << "apop\t"
        << "dead\t"
        << "tvas\t"
        << "ttis\t"
        << "chem\t"
        << "icyt\t"
        << "virs";
    return oss.str();
  }

  string to_str() {
    auto tot_incubating = reduce_one(incubating, op_fast_add, 0).wait();
    auto tot_expressing = reduce_one(expressing, op_fast_add, 0).wait();
    auto tot_apoptotic = reduce_one(apoptotic, op_fast_add, 0).wait();
    auto tot_dead = reduce_one(dead, op_fast_add, 0).wait();
    auto tot_tcells_vasculature = reduce_one(tcells_vasculature, op_fast_add, 0).wait();
    auto tot_tcells_tissue = reduce_one(tcells_tissue, op_fast_add, 0).wait();
    int64_t num_grid_points = Tissue::grid_size.x * Tissue::grid_size.y * Tissue::grid_size.z;
    auto tot_chemokines = reduce_one(chemokines, op_fast_add, 0).wait() / num_grid_points;
    auto tot_icytokines = reduce_one(icytokines, op_fast_add, 0).wait() / num_grid_points;
    auto tot_virus = reduce_one(virus, op_fast_add, 0).wait() / num_grid_points;
    ostringstream oss;
    oss << left << tot_incubating << "\t" << tot_expressing << "\t" << tot_apoptotic << "\t"
        << tot_dead << "\t" << tot_tcells_vasculature << "\t" << tot_tcells_tissue << "\t" << fixed
        << setprecision(4) << tot_chemokines << "\t" << tot_icytokines << "\t" << tot_virus;
    return oss.str();
  }

  void log() {
    string s = to_str();
    if (!rank_me()) log_file << s << endl;
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
IntermittentTimer _switch_round_timer(__FILENAME__ + string(":") + "switch_round");

void initial_infection(Tissue &tissue) {
  BarrierTimer timer(__FILEFUNC__);
  int local_num_infections = 0;
  if (_options->infection_coords[0] != -1) {
    if (!rank_me()) {
      local_num_infections = 1;
      tissue.inc_incoming_virus({_options->infection_coords[0], _options->infection_coords[1],
                                 _options->infection_coords[2]},
                                1.0);
    } else {
      local_num_infections = 0;
    }
  } else {
    // each rank generates a block of infected cells
    // for large dimensions, the chance of repeat sampling for a few points is very small
    local_num_infections = _options->num_infections / rank_n();
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
  }
  barrier();
  SLOG("Initially infected ", reduce_one(local_num_infections, op_fast_add, 0).wait(),
       " epicells\n");
  tissue.add_new_actives();
  barrier();
  int num_infections_found = 0;
  for (auto grid_point = tissue.get_first_active_grid_point(); grid_point;
       grid_point = tissue.get_next_active_grid_point()) {
    if (grid_point->incoming_virus > 0) {
      grid_point->virus = 1.0;
      grid_point->incoming_virus = 0;
      grid_point->epicell->infect();
      num_infections_found++;
      _sim_stats.incubating++;
    }
  }
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
      _sim_stats.tcells_vasculature++;
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
  if (grid_point->icytokine > 0 && tcell.in_vasculature) {
    tcell.in_vasculature = false;
    tcell.prev_coords = grid_point->coords;
    _sim_stats.tcells_vasculature--;
    _sim_stats.tcells_tissue++;
    DBG(time_step, ": tcell ", tcell.id, " extravasates at ", grid_point->coords.str(),"\n");
  }
  if (tcell.in_vasculature) {
    tcell.vascular_period--;
    if (tcell.vascular_period > 0) {
      // tcell is still alive - moves to any random location in the grid
      tissue.add_tcell({_rnd_gen, Tissue::grid_size}, tcell);
    } else {
      _sim_stats.tcells_vasculature--;
      DBG(time_step, ": tcell ", tcell.id, " dies in vasculature\n");
    }
  } else {
    tcell.tissue_period--;
    // tcell is in the tissue
    if (tcell.tissue_period > 0) {
      // still alive
      GridCoords selected_coords;
      // FIXME: I believe that T cells should be able to detect virus in incubating cells too,
      // just with a lower probability. Note that the TCRs are binding to viral peptides, not
      // virions, so these should be detectable through MHC transport even before production of
      // complete virions within the epicell.
      if (grid_point->epicell->status == EpiCellStatus::EXPRESSING) {
        DBG(time_step, ": tcell ", tcell.id, " is inducing apoptosis at ", grid_point->coords.str(),
            "\n");
        grid_point->epicell->status = EpiCellStatus::APOPTOTIC;
        _sim_stats.expressing--;
        _sim_stats.apoptotic++;
        // we are inducing apoptosis, so don't move
        selected_coords = grid_point->coords;
      } else {
        double highest_chemokine = 0;
        for (auto &nb_coords : grid_point->neighbors) {
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
          DBG(time_step, ": tcell ", tcell.id, " at ", grid_point->coords.str(),
              " highest nb chemokine is at ", selected_coords.str(), " with ",
              highest_chemokine, "\n");
        }
      }
      DBG(time_step, ": tcell ", tcell.id, " at ", grid_point->coords.str(), " moving to ",
          selected_coords.str(), "\n");
      tcell.prev_coords = grid_point->coords;
      tissue.add_tcell(selected_coords, tcell);
    } else {
      _sim_stats.tcells_tissue--;
      DBG(time_step, ": tcell ", tcell.id, " dies in tissue at " , grid_point->coords.str(), "\n");
    }
  }
  _update_tcell_timer.stop();
}

void update_epicell(int time_step, Tissue &tissue, GridPoint *grid_point) {
  if (grid_point->epicell->status == EpiCellStatus::DEAD) return;
  _update_epicell_timer.start();
  switch (grid_point->epicell->status) {
    case EpiCellStatus::HEALTHY:
      if (grid_point->virus > 0) {
        double infection_prob = grid_point->virus * _options->infection_prob;
        if (_rnd_gen->trial_success(infection_prob)) {
          grid_point->epicell->infect();
          _sim_stats.incubating++;
        }
      }
      break;
    case EpiCellStatus::INCUBATING:
      if (grid_point->epicell->transition_to_expressing()) {
        _sim_stats.incubating--;
        _sim_stats.expressing++;
        grid_point->virus = 1.0;
      }
      break;
    case EpiCellStatus::EXPRESSING:
      if (grid_point->epicell->infection_death()) {
        grid_point->virus = 0;
        _sim_stats.dead++;
        _sim_stats.expressing--;
      } else {
        grid_point->virus = 1.0;
        grid_point->incoming_virus = 1.0;
        grid_point->chemokine = 1.0;
        grid_point->incoming_chemokine = 1.0;
        grid_point->icytokine = 1.0;
        grid_point->incoming_icytokine = 1.0;
        // apoptosis is induced directly by a tcell in update_tcell
      }
      break;
    case EpiCellStatus::APOPTOTIC:
      // FIXME: it seems that there is evidence that apoptotic cells also produce cytokines
      // This would be a necessary requirement for a feedback loop of damage caused by autoreactive
      // tcells.
      if (grid_point->epicell->apoptosis_death()) {
        grid_point->virus = 0;
        _sim_stats.dead++;
        _sim_stats.apoptotic--;
      // FIXME: no secretion of virions during apoptosis? It seems that virion secretion is
      // inhibited, but it may still happen
      /*
      } else {
        if (grid_point->epicell->is_fully_incubated()) {
          grid_point->virus = 1.0;
          grid_point->incoming_virus = 1.0;
        }
      */
      }
      break;
    default: break;
  }
  _update_epicell_timer.stop();
}

void update_concentration(int time_step, GridPoint *grid_point,
                          double &concentration, double decay_rate, double diffusion_coef,
                          Tissue &tissue, void (Tissue::*inc_incoming)(GridCoords, double)) {
  _update_concentration_timer.start();
  if (concentration > 0) {
    concentration -= decay_rate;
    if (concentration < 0) concentration = 0;
    if (concentration > 0) {
      // diffuses, reducing concentration at source
      double diffusion_amount = concentration * diffusion_coef;
      concentration -= diffusion_amount;
      // the amount diffusing spreads uniformly to all neighbors
      double amount_to_nb = diffusion_amount / grid_point->neighbors.size();
      for (auto &nb_coords : grid_point->neighbors) {
        assert(nb_coords != grid_point->coords);
        (tissue.*inc_incoming)(nb_coords, amount_to_nb);
      }
    }
  }
  _update_concentration_timer.stop();
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
    //DBG("Truncated sample file ", fname, " to ", tot_sz, " bytes\n");
  }
  upcxx::barrier();
  tot_sz = tissue.get_num_local_grid_points();
  if (!rank_me()) tot_sz += header_oss.str().length();
  // wait until rank 0 has finished setting up the file
  auto [bytes_written, grid_points_written] = tissue.dump_blocks(fname, header_oss.str(),
                                                                 view_object);
  //assert(bytes_written == tot_sz);
  if (bytes_written != tot_sz) DIE("bytes_written ", bytes_written, " != ", tot_sz, " tot_sz");
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
  _sim_stats.init();
  SLOG("# datetime     elapsed step    ", _sim_stats.header(), "\n");
  for (int time_step = 0; time_step < _options->num_timesteps; time_step++) {
    DBG("Time step ", time_step, "\n");
    if (time_step > _options->tcell_initial_delay) generate_tcells(tissue);
    // iterate through all active local grid points and update
    for (auto grid_point = tissue.get_first_active_grid_point(); grid_point;
         grid_point = tissue.get_next_active_grid_point()) {
      DBG("updating grid point ", grid_point->str(), "\n");
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
      if (grid_point->is_active()) tissue.set_active(grid_point);
    }
    barrier();
    if (time_step % _options->sample_period == 0 || time_step == _options->num_timesteps - 1) {
      chrono::duration<double> t_elapsed = NOW() - curr_t;
      if (_options->verbose) curr_t = NOW();
      SLOG_VERBOSE("[", get_current_time(true), " ", setprecision(2), fixed, setw(5), right,
                   t_elapsed.count(), "s]: ", setw(8), left, time_step, _sim_stats.to_str(), "\n");
    }
    if (!_options->verbose &&
        (time_step % five_perc == 0 || time_step == _options->num_timesteps - 1)) {
      chrono::duration<double> t_elapsed = NOW() - curr_t;
      curr_t = NOW();
      SLOG("[", get_current_time(true), " ", setprecision(2), fixed, setw(5), right,
           t_elapsed.count(), "s]: ", setw(8), left, time_step, _sim_stats.to_str(), "\n");
    }
    barrier();
    tissue.add_new_actives();
    barrier();
    _switch_round_timer.start();
    _sim_stats.virus = 0;
    _sim_stats.chemokines = 0;
    _sim_stats.icytokines = 0;
    vector<GridPoint *> to_erase = {};
    // iterate through all active local grid points and set changes
    for (auto grid_point = tissue.get_first_active_grid_point(); grid_point;
         grid_point = tissue.get_next_active_grid_point()) {
      if (grid_point->tcells) grid_point->switch_tcells_vector();
      if (grid_point->incoming_virus > 0) {
        grid_point->virus += grid_point->incoming_virus;
        if (grid_point->virus > 1) grid_point->virus = 1;
        grid_point->incoming_virus = 0;
      }
      if (grid_point->incoming_chemokine > 0) {
        grid_point->chemokine += grid_point->incoming_chemokine;
        if (grid_point->chemokine > 1) grid_point->chemokine = 1;
        grid_point->incoming_chemokine = 0;
      }
      if (grid_point->incoming_icytokine > 0) {
        grid_point->icytokine += grid_point->incoming_icytokine;
        if (grid_point->icytokine > 1) grid_point->icytokine = 1;
        grid_point->incoming_icytokine = 0;
      }
      _sim_stats.virus += grid_point->virus;
      _sim_stats.chemokines += grid_point->chemokine;
      _sim_stats.icytokines += grid_point->icytokine;
      if (!grid_point->is_active()) to_erase.push_back(grid_point);
    }
    for (auto grid_point : to_erase) tissue.erase_active(grid_point);
    _switch_round_timer.stop();
    barrier();
    sample(time_step, tissue, ViewObject::EPICELL);
    sample(time_step, tissue, ViewObject::TCELL_TISSUE);
    sample(time_step, tissue, ViewObject::VIRUS);
    sample(time_step, tissue, ViewObject::ICYTOKINE);
    sample(time_step, tissue, ViewObject::CHEMOKINE);
    _sim_stats.log();
#ifdef DEBUG
    DBG("check actives ", time_step, "\n");
    tissue.check_actives(time_step);
    barrier();
#endif
  }
  _generate_tcell_timer.done_all();
  _update_tcell_timer.done_all();
  _update_epicell_timer.done_all();
  _update_concentration_timer.done_all();
  _sample_timer.done_all();
  _switch_round_timer.done_all();

  chrono::duration<double> t_elapsed = NOW() - start_t;
  SLOG("Finished ", _options->num_timesteps, " time steps in ", setprecision(4), fixed,
       t_elapsed.count(), " s (", (double)t_elapsed.count() / _options->num_timesteps,
       " s per step)\n");
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
