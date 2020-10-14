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
  double virions = 0;

  void init() {
    if (!rank_me()) {
      log_file.open("simcov.stats");
      log_file << "# time\t" << header() << endl;
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
    auto tot_chemokines = reduce_one(chemokines, op_fast_add, 0).wait();
    auto tot_virions = reduce_one(virions, op_fast_add, 0).wait();

    tot_chemokines /= Tissue::get_num_grid_points();
    tot_virions /= Tissue::get_num_grid_points();

    ostringstream oss;
    oss << left << tot_incubating << "\t" << tot_expressing << "\t" << tot_apoptotic << "\t"
        << tot_dead << "\t" << tot_tcells_vasculature << "\t" << tot_tcells_tissue << "\t" << fixed
        << setprecision(2) << scientific << tot_chemokines << "\t" << tot_virions;
    return oss.str();
  }

  void log(int time_step) {
    string s = to_str();
    if (!rank_me()) log_file << time_step << "\t" << s << endl;
  }

};

ofstream _logstream;
bool _verbose = false;
SimStats _sim_stats;
shared_ptr<Options> _options;

IntermittentTimer generate_tcell_timer(__FILENAME__ + string(":") + "generate tcells");
IntermittentTimer update_circulating_tcells_timer(__FILENAME__ + string(":") +
                                                  "update circulating tcells");
IntermittentTimer update_tcell_timer(__FILENAME__ + string(":") + "update tcells");
IntermittentTimer update_epicell_timer(__FILENAME__ + string(":") + "update epicells");
IntermittentTimer update_concentration_timer(__FILENAME__ + string(":") + "update concentrations");
IntermittentTimer compute_updates_timer(__FILENAME__ + string(":") + "compute updates");
IntermittentTimer accumulate_concentrations_timer(__FILENAME__ + string(":") + "dispatch updates");
IntermittentTimer add_new_actives_timer(__FILENAME__ + string(":") + "add new actives");
IntermittentTimer set_active_points_timer(__FILENAME__ + string(":") + "erase inactive");
IntermittentTimer sample_timer(__FILENAME__ + string(":") + "sample");
IntermittentTimer log_timer(__FILENAME__ + string(":") + "log");


void initial_infection(Tissue &tissue) {
  BarrierTimer timer(__FILEFUNC__);
  int local_num_infections = 0;
  if (_options->infection_coords[0] != -1) {
    if (!rank_me()) {
      local_num_infections = 1;
      GridCoords coords = {_options->infection_coords[0], _options->infection_coords[1],
                           _options->infection_coords[2]};
      if (tissue.set_initial_infection(coords)) _sim_stats.incubating++;
    } else {
      local_num_infections = 0;
    }
  } else {
    // each rank generates a block of infection points
    // for large dimensions, the chance of repeat sampling for a few points is very small
    local_num_infections = _options->num_infections / rank_n();
    int remaining_infections = _options->num_infections - local_num_infections * rank_n();
    if (rank_me() < remaining_infections) local_num_infections++;
    DBG("Virions in ", local_num_infections, " locations\n");
    ProgressBar progbar(local_num_infections, "Setting initial infections");
    for (int i = 0; i < local_num_infections; i++) {
      progbar.update();
      GridCoords coords(_rnd_gen);
      DBG("infection: ", coords.str() + "\n");
      if (tissue.set_initial_infection(coords)) _sim_stats.incubating++;
      upcxx::progress();
    }
    progbar.done();
  }
  barrier();
  tissue.add_new_actives(add_new_actives_timer);
  barrier();
}

void add_tcell(Tissue &tissue, GridCoords coords, TCell tcell,
               HASH_TABLE<intrank_t, vector<pair<GridCoords, TCell>>> &tcells_to_add) {
  tcells_to_add[tissue.get_rank_for_grid_point(coords)].push_back({coords, tcell});
}

void generate_tcells(Tissue &tissue, int time_step) {
  generate_tcell_timer.start();
  double generation_rate = _options->tcell_generation_rate;
  int local_num_tcells = generation_rate / rank_n();
  int remaining_tcells = generation_rate - local_num_tcells * rank_n();
  if (rank_me() < remaining_tcells) local_num_tcells++;
  double frac = generation_rate - floor(generation_rate);
  if (rank_me() == 0 && frac > 0 && _rnd_gen->trial_success(frac)) local_num_tcells++;
  if (local_num_tcells) {
    for (int i = 0; i < local_num_tcells; i++) {
      GridCoords coords(_rnd_gen);
      string tcell_id = to_string(rank_me()) + "-" + to_string(tissue.tcells_generated);
      tissue.tcells_generated++;
      TCell tcell(tcell_id, coords);
      // choose a destination rank to add the tcell to based on the inital random coords
      tissue.add_circulating_tcell(coords, tcell);
      _sim_stats.tcells_vasculature++;
      upcxx::progress();
    }
  }
  generate_tcell_timer.stop();
}

int64_t get_rnd_coord(int64_t x, int64_t max_x) {
  int64_t new_x = x + _rnd_gen->get(0, 3) - 1;
  if (new_x < 0) new_x = 0;
  if (new_x >= max_x) new_x = max_x - 1;
  return new_x;
}

void update_circulating_tcells(int time_step, Tissue &tissue) {
  update_circulating_tcells_timer.start();
  auto circulating_tcells = tissue.get_circulating_tcells();
  for (auto it = circulating_tcells->begin(); it != circulating_tcells->end(); it++) {
    progress();
    it->vascular_period--;
    if (it->vascular_period == 0) {
      _sim_stats.tcells_vasculature--;
      circulating_tcells->erase(it--);
      continue;
    }
    GridCoords coords(_rnd_gen);
    if (tissue.try_add_tissue_tcell(coords, *it, true)) {
      _sim_stats.tcells_vasculature--;
      _sim_stats.tcells_tissue++;
      DBG(time_step, " tcell ", it->id, " extravasates at ", coords.str(), " with lifespan ",
          it->tissue_period, " will die at ", (time_step + it->tissue_period), "\n");
      circulating_tcells->erase(it--);
    }
  }
  update_circulating_tcells_timer.stop();
}

void update_tissue_tcell(int time_step, Tissue &tissue, GridPoint *grid_point,
                         HASH_TABLE<int64_t, double> &chemokines_cache) {
  update_tcell_timer.start();
  TCell *tcell = grid_point->tcell;
  if (tcell->moved) {
    // don't update tcells that were added this time step
    tcell->moved = false;
    update_tcell_timer.stop();
    return;
  }
  tcell->tissue_period--;
  if (tcell->tissue_period == 0) {
    _sim_stats.tcells_tissue--;
    DBG(time_step, " tcell ", tcell->id, " dies in tissue at ", grid_point->coords.str(), "\n");
    // not adding to a new location means this tcell is not preserved to the next time step
    delete grid_point->tcell;
    grid_point->tcell = nullptr;
    update_tcell_timer.stop();
    return;
  }
  if (tcell->binding_period != -1) {
    DBG(time_step, " tcell ", tcell->id, " is bound at ", grid_point->coords.str(), "\n");
    // this tcell is bound
    tcell->binding_period--;
    // done with binding when set to -1
    if (tcell->binding_period == 0) tcell->binding_period = -1;
  } else {
    // not bound to an epicell
    if (grid_point->epicell && (grid_point->epicell->status == EpiCellStatus::EXPRESSING ||
        grid_point->epicell->status == EpiCellStatus::INCUBATING)) {
      double binding_prob = grid_point->epicell->get_binding_prob();
      DBG(time_step, " tcell ", tcell->id, " trying to bind at ", grid_point->coords.str(), "\n");
      if (_rnd_gen->trial_success(binding_prob)) {
        DBG(time_step, " tcell ", tcell->id, " is inducing apoptosis at ", grid_point->coords.str(),
            "\n");
        if (grid_point->epicell->status == EpiCellStatus::EXPRESSING) _sim_stats.expressing--;
        if (grid_point->epicell->status == EpiCellStatus::INCUBATING) _sim_stats.incubating--;
        // as soon as this is set to apoptotic, no other tcell can bind to this epicell, so we
        // ensure that only one tcell binds to an epicell at a time
        grid_point->epicell->status = EpiCellStatus::APOPTOTIC;
        _sim_stats.apoptotic++;
        tcell->binding_period = _options->tcell_binding_period;
      }
    }
  }
  if (tcell->binding_period == -1) {
    DBG(time_step, " tcell ", tcell->id, " trying to move at ", grid_point->coords.str(), "\n");
    // didn't bind - move on chemokine gradient or at random
    GridCoords selected_coords;
    // not bound - follow chemokine gradient
    double highest_chemokine = 0;
    // get a randomly shuffled list of neighbors so the tcell doesn't always tend to move in the
    // same direction when there is a chemokine gradient
    auto nbs_shuffled = grid_point->neighbors;
    random_shuffle(nbs_shuffled.begin(), nbs_shuffled.end());
    //for (auto &nb_coords : grid_point->neighbors) {
    for (auto &nb_coords : nbs_shuffled) {
      double chemokine = 0;
      int64_t nb_idx = nb_coords.to_1d();
      auto it = chemokines_cache.find(nb_idx);
      if (it == chemokines_cache.end()) {
        chemokine = tissue.get_chemokine(nb_coords);
        chemokines_cache.insert({nb_idx, chemokine});
      } else {
        chemokine = it->second;
      }
      if (chemokine > highest_chemokine) {
        highest_chemokine = chemokine;
        selected_coords = nb_coords;
      }
      DBG(time_step, " tcell ", tcell->id, " found nb chemokine ", chemokine, " at ",
          nb_coords.str(), "\n");
    }
    if (highest_chemokine == 0) {
      // no chemokines found - move randomly
      auto rnd_nb_i = _rnd_gen->get(0, (int64_t)grid_point->neighbors.size());
      selected_coords = grid_point->neighbors[rnd_nb_i];
      DBG(time_step, " tcell ", tcell->id, " try random move to ", selected_coords.str(), "\n");
    } else {
      DBG(time_step, " tcell ", tcell->id, " - highest chemokine at ", selected_coords.str(), "\n");
    }
    // try a few times to find an open spot
    for (int i = 0; i < 5; i++) {
      if (tissue.try_add_tissue_tcell(selected_coords, *tcell, false)) {
        DBG(time_step, " tcell ", tcell->id, " at ", grid_point->coords.str(), " moves to ",
            selected_coords.str(), "\n");
        delete grid_point->tcell;
        grid_point->tcell = nullptr;
        break;
      }
      // choose another location at random
      auto rnd_nb_i = _rnd_gen->get(0, (int64_t)grid_point->neighbors.size());
      selected_coords = grid_point->neighbors[rnd_nb_i];
      DBG(time_step, " tcell ", tcell->id, " try random move to ", selected_coords.str(), "\n");
    }
  }
  update_tcell_timer.stop();
}

void update_epicell(int time_step, Tissue &tissue, GridPoint *grid_point) {
  update_epicell_timer.start();
  if (!grid_point->epicell->infectable || grid_point->epicell->status == EpiCellStatus::DEAD) {
    update_epicell_timer.stop();
    return;
  }
  if (grid_point->epicell->status != EpiCellStatus::HEALTHY)
    DBG(time_step, " epicell status ", EpiCellStatusStr[(int)grid_point->epicell->status], "\n");
  bool produce_virions = false;
  switch (grid_point->epicell->status) {
    case EpiCellStatus::HEALTHY:
      if (grid_point->virions > 0) {
        if (_rnd_gen->trial_success(_options->infectivity * grid_point->virions)) {
          grid_point->epicell->infect();
          _sim_stats.incubating++;
        }
      }
      break;
    case EpiCellStatus::INCUBATING:
      if (grid_point->epicell->transition_to_expressing()) {
        _sim_stats.incubating--;
        _sim_stats.expressing++;
      }
      break;
    case EpiCellStatus::EXPRESSING:
      if (grid_point->epicell->infection_death()) {
        _sim_stats.dead++;
        _sim_stats.expressing--;
      } else {
        produce_virions = true;
    }
      break;
    case EpiCellStatus::APOPTOTIC:
      // FIXME: it seems that there is evidence that apoptotic cells also produce cytokines
      // This would be a necessary requirement for a feedback loop of damage caused by autoreactive
      // tcells.
      if (grid_point->epicell->apoptosis_death()) {
        _sim_stats.dead++;
        _sim_stats.apoptotic--;
      } else {
        produce_virions = true;
      }
      break;
    default: break;
  }
  if (produce_virions) {
    grid_point->virions += _options->virion_production;
    grid_point->chemokine = min(grid_point->chemokine + _options->chemokine_production, 1.0);
  }
  update_epicell_timer.stop();
}

void update_chemokines(GridPoint *grid_point,
                       HASH_TABLE<int64_t, double> &chemokines_to_update) {
  update_concentration_timer.start();
  // Concentrations diffuse, i.e. the concentration at any single grid point tends to the average of
  // all the neighbors. So here we tell each neighbor what the current concentration is and later
  // those neighbors will compute their own averages. We do it in this "push" manner because then we
  // don't need to check the neighbors from every single grid point, but just push from ones with
  // concentrations > 0 (i.e. active grid points)
  if (grid_point->chemokine > 0) {
    grid_point->chemokine *= (1.0 - _options->chemokine_decay_rate);
    if (grid_point->chemokine < MIN_CONCENTRATION) grid_point->chemokine = 0;
  }
  if (grid_point->chemokine > 0) {
    for (auto &nb_coords : grid_point->neighbors) {
      assert(nb_coords != grid_point->coords);
      chemokines_to_update[nb_coords.to_1d()] += grid_point->chemokine;
    }
  }
  update_concentration_timer.stop();
}

void update_virions(GridPoint *grid_point, HASH_TABLE<int64_t, double> &virions_to_update) {
  update_concentration_timer.start();
  grid_point->virions = grid_point->virions * (1.0 - _options->virion_decay_rate);
  if (grid_point->virions < 1) grid_point->virions = 0;
  if (grid_point->virions > 0) {
    for (auto &nb_coords : grid_point->neighbors) {
      assert(nb_coords != grid_point->coords);
      virions_to_update[nb_coords.to_1d()] += grid_point->virions;
    }
  }
  update_concentration_timer.stop();
}

void diffuse(double &conc, double &nb_conc, double diffusion, int num_nbs, bool max_limit=true) {
  // set to be average of neighbors plus self
  // amount that diffuses
  double conc_diffused = diffusion * conc;
  // average out diffused amount across all neighbors
  double conc_per_point = (conc_diffused + diffusion * nb_conc) / (num_nbs + 1);
  conc = conc - conc_diffused + conc_per_point;
  if (max_limit && conc > 1.0) conc = 1.0;
  if (conc < 0) DIE("conc < 0: ", conc, " diffused ", conc_diffused, " pp ", conc_per_point);
  nb_conc = 0;
}

/*
void spread_virions(int &virions, int &nb_virions, double diffusion, int num_nbs) {
  int virions_diffused = virions * diffusion;
  int virions_left = virions - virions_diffused;
  int avg_nb_virions = (virions_diffused + nb_virions * diffusion) / (num_nbs + 1);
  virions = virions_left + avg_nb_virions;
  nb_virions = 0;
}
*/

void set_active_grid_points(Tissue &tissue) {
  set_active_points_timer.start();
  vector<GridPoint *> to_erase = {};
  // iterate through all active local grid points and set changes
  for (auto grid_point = tissue.get_first_active_grid_point(); grid_point;
       grid_point = tissue.get_next_active_grid_point()) {
    diffuse(grid_point->chemokine, grid_point->nb_chemokine, _options->chemokine_diffusion_coef,
            grid_point->neighbors.size());
    diffuse(grid_point->virions, grid_point->nb_virions, _options->virion_diffusion_coef,
            grid_point->neighbors.size(), false);
    if (grid_point->chemokine < MIN_CONCENTRATION) grid_point->chemokine = 0;
    if (grid_point->virions < 1.0) grid_point->virions = 0;
    if (grid_point->tcell) grid_point->tcell->moved = false;
    _sim_stats.chemokines += grid_point->chemokine;
    _sim_stats.virions += grid_point->virions;
    if (!grid_point->is_active()) to_erase.push_back(grid_point);
  }
  for (auto grid_point : to_erase) tissue.erase_active(grid_point);
  set_active_points_timer.stop();
}

void sample(int time_step, Tissue &tissue, ViewObject view_object) {
  // each rank writes its own blocks into the file at the appropriate locations. We compute the
  // location knowing that each position takes exactly 4 characters (up to 3 numbers and a space),
  // so we can work out how much space a block takes up, in order to position the data correctly.
  // Each point is represented by a scalar char. From 1-127 are used to indicate viral load, with a
  // color of red, going from faint to strong. From -127 to -1 are used to indicate tcell count,
  // with a color of blue, going from faint to strong. The last color (0) is used to indicate a dead
  // epicell (how do we make this a separate color?)
  // each grid point takes up a single char
  size_t tot_sz = tissue.get_num_grid_points();
  string fname_snapshot = "sample_snapshot.csv";
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
}

void run_sim(Tissue &tissue) {
  BarrierTimer timer(__FILEFUNC__);

  auto start_t = NOW();
  auto curr_t = start_t;
  auto five_perc = _options->num_timesteps / 50;
  _sim_stats.init();
  SLOG("# datetime     elapsed step    ", _sim_stats.header(),
       "\t<%active  lbln>\n");
  // store the total concentration increment updates for target grid points
  // chemokine, virions
  HASH_TABLE<int64_t, double> chemokines_to_update;
  HASH_TABLE<int64_t, double> chemokines_cache;
  HASH_TABLE<int64_t, double> virions_to_update;
  bool warned_boundary = false;
  //bool expressed = false;
  for (int time_step = 0; time_step < _options->num_timesteps; time_step++) {
    DBG("Time step ", time_step, "\n");
    if (time_step == _options->antibody_period)
      _options->virion_decay_rate *= _options->antibody_factor;
    chemokines_to_update.clear();
    virions_to_update.clear();
    chemokines_cache.clear();
    if (time_step > _options->tcell_initial_delay) {
      generate_tcells(tissue, time_step);
      barrier();
    }
    compute_updates_timer.start();
    update_circulating_tcells(time_step, tissue);
    // iterate through all active local grid points and update
    for (auto grid_point = tissue.get_first_active_grid_point(); grid_point;
         grid_point = tissue.get_next_active_grid_point()) {
      if (grid_point->chemokine > 0)
        DBG("chemokine\t", time_step, "\t", grid_point->coords.x, "\t", grid_point->coords.y, "\t",
            grid_point->coords.z, "\t", grid_point->chemokine, "\n");
      if (grid_point->virions > 0)
        DBG("virions\t", time_step, "\t", grid_point->coords.x, "\t", grid_point->coords.y, "\t",
            grid_point->coords.z, "\t", grid_point->virions, "\n");
      /*
      if (grid_point->virions > 0 && !expressed) {
        expressed = true;
        _options->virion_production = 0;
      }
      */
      if (!warned_boundary && (!grid_point->coords.x || !grid_point->coords.y ||
          (_grid_size->z > 1 && !grid_point->coords.z)) &&
          grid_point->epicell->status != EpiCellStatus::HEALTHY) {
        WARN("Hit boundary at ", grid_point->coords.str(), " ", grid_point->epicell->str(),
              " virions ", grid_point->virions, " chemokine ", grid_point->chemokine);
        warned_boundary = true;
      }
      //DBG("updating grid point ", grid_point->str(), "\n");
      upcxx::progress();
      // the tcells are moved (added to the new list, but only cleared out at the end of all
      // updates)
      if (grid_point->tcell) update_tissue_tcell(time_step, tissue, grid_point, chemokines_cache);
      update_epicell(time_step, tissue, grid_point);
      update_chemokines(grid_point, chemokines_to_update);
      update_virions(grid_point, virions_to_update);
      if (grid_point->is_active()) tissue.set_active(grid_point);
    }
    barrier();
    compute_updates_timer.stop();
    tissue.accumulate_chemokines(chemokines_to_update, accumulate_concentrations_timer);
    tissue.accumulate_virions(virions_to_update, accumulate_concentrations_timer);
    barrier();
    if (time_step % five_perc == 0 || time_step == _options->num_timesteps - 1) {
      auto num_actives = reduce_one(tissue.get_num_actives(), op_fast_add, 0).wait();
      auto perc_actives = 100.0 * num_actives / Tissue::get_num_grid_points();
      auto max_actives = reduce_one(tissue.get_num_actives(), op_fast_max, 0).wait();
      auto load_balance = max_actives ? (double)num_actives / rank_n() / max_actives : 1;
      chrono::duration<double> t_elapsed = NOW() - curr_t;
      curr_t = NOW();
      SLOG("[", get_current_time(true), " ", setprecision(2), fixed, setw(5), right,
           t_elapsed.count(), "s]: ", setw(8), left, time_step, _sim_stats.to_str());
      SLOG(setprecision(3), fixed, "\t< ", perc_actives, " ", load_balance, " >\n");
    }
    barrier();
    tissue.add_new_actives(add_new_actives_timer);
    barrier();

    _sim_stats.virions = 0;
    _sim_stats.chemokines = 0;
    set_active_grid_points(tissue);
    barrier();

    if (_options->sample_period > 0 && (time_step % _options->sample_period == 0 ||
        time_step == _options->num_timesteps - 1)) {
      sample_timer.start();
      sample(time_step, tissue, ViewObject::EPICELL);
      sample(time_step, tissue, ViewObject::TCELL_TISSUE);
      sample(time_step, tissue, ViewObject::VIRUS);
      sample(time_step, tissue, ViewObject::CHEMOKINE);
      sample_timer.stop();
    }

    log_timer.start();
    _sim_stats.log(time_step);
    barrier();
    log_timer.stop();

#ifdef DEBUG
    DBG("check actives ", time_step, "\n");
    tissue.check_actives(time_step);
    barrier();
#endif
  }

  generate_tcell_timer.done_all();
  update_circulating_tcells_timer.done_all();
  update_tcell_timer.done_all();
  update_epicell_timer.done_all();
  update_concentration_timer.done_all();
  compute_updates_timer.done_all();
  accumulate_concentrations_timer.done_all();
  add_new_actives_timer.done_all();
  set_active_points_timer.done_all();
  sample_timer.done_all();
  log_timer.done_all();

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
