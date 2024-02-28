// SimCov
//
// Steven Hofmeyr, LBNL May 2020

#include <fcntl.h>
#include <math.h>
#include <stdarg.h>
#include <unistd.h>

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <upcxx/upcxx.hpp>
#include <vector>
#include <random>
#include <unordered_set>

using namespace std;

#include "options.hpp"
#include "tissue.hpp"
#include "upcxx_utils.hpp"
#include "utils.hpp"

using namespace upcxx;
using namespace upcxx_utils;

#define NOW chrono::high_resolution_clock::now
#define STATS_COL_WIDTH 11

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
  float virions = 0;
  int64_t num_inflam_signal_all = 0;
  int64_t num_inflam_signal_cell = 0;
  //added
  int64_t num_aircell = 0;
  int64_t num_intertsitial = 0;
  int64_t num_infected = 0;
  
  int64_t num_aircell_w_inflam = 0;
  int64_t num_intertsitial_w_inflam = 0;
  int64_t num_infected_w_inflam = 0;
  
  

  void init() {
    if (!rank_me()) {
      log_file.open(_options->output_dir + "/simcov.stats");
      log_file << "# time" << header(0) << endl;
    }
  }

  string header(int width) {
    vector<string> columns = {"incb", "expr", "apop", "dead", "tvas",
                              "ttis", "inflcell", "inflall", "air" , "interstitial", "infection", "infl_air", "infl_ins", "infl_inf", "virs", "%infct"  };
    ostringstream oss;
    oss << left;
    for (auto column : columns) {
      if (width)
        oss << setw(width) << column;
      else
        oss << '\t' << column;
    }
    return oss.str();
  }

  string to_str(int width) {
    vector<int64_t> totals;
    totals.push_back(reduce_one(incubating, op_fast_add, 0).wait());
    totals.push_back(reduce_one(expressing, op_fast_add, 0).wait());
    totals.push_back(reduce_one(apoptotic, op_fast_add, 0).wait());
    totals.push_back(reduce_one(dead, op_fast_add, 0).wait());
    totals.push_back(reduce_one(tcells_vasculature, op_fast_add, 0).wait());
    totals.push_back(reduce_one(tcells_tissue, op_fast_add, 0).wait());
    totals.push_back(reduce_one(num_inflam_signal_cell, op_fast_add, 0).wait());
    totals.push_back(reduce_one(num_inflam_signal_all, op_fast_add, 0).wait());

	
    vector<float> totals_d;
    totals_d.push_back(reduce_one(virions, op_fast_add, 0).wait());
    auto perc_infected =
        100.0 * (float)(totals[0] + totals[1] + totals[2] + totals[3]) / get_num_grid_points();
    if (!_options->lung_model_dir.empty()) {
      perc_infected = 100.0
          * (float)(totals[0] + totals[1] + totals[2] + totals[3])
          / Tissue::get_num_lung_cells();
    }
	
	//added the serial matters
	totals.push_back(reduce_one(num_aircell, op_fast_add, 0).wait());
	totals.push_back(reduce_one(num_intertsitial, op_fast_add, 0).wait());
	totals.push_back(reduce_one(num_infected, op_fast_add, 0).wait());
	
	totals.push_back(reduce_one(num_aircell_w_inflam, op_fast_add, 0).wait());
	totals.push_back(reduce_one(num_intertsitial_w_inflam, op_fast_add, 0).wait());
	totals.push_back(reduce_one(num_infected_w_inflam, op_fast_add, 0).wait());


    ostringstream oss;
    oss << left;
    for (auto tot : totals) {
      if (width)
        oss << setw(width) << tot;
      else
        oss << '\t' << tot;
    }
    oss << fixed << setprecision(2) << scientific;
    for (auto tot : totals_d) {
      if (width)
        oss << setw(width) << tot;
      else
        oss << '\t' << tot;
    }
    oss << fixed << setprecision(10) << showpoint;
    if (width)
      oss << setw(width) << perc_infected;
    else
      oss << '\t' << perc_infected;
    return oss.str();
  }

  void log(int time_step) {
    string s = to_str(0);
    if (!rank_me()) log_file << time_step << s << endl;
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
IntermittentTimer sample_write_timer(__FILENAME__ + string(":") + "sample write");
IntermittentTimer log_timer(__FILENAME__ + string(":") + "log");

struct Point {
    double x, y, z;

    // Define less-than operator for set comparison
    bool operator<(const Point& other) const {
        if (x != other.x) return x < other.x;
        if (y != other.y) return y < other.y;
        return z < other.z;
    }
};

struct Point3D {
    double x;
    double y;
    double z;
};

// Function to generate all points in a sphere given the center
std::vector<Point> generateAllPointsInSphere(Point center, double radius, int numPoints) {
    std::vector<Point> points;
	set<Point> generatedPoints;  // Set to store generated points

    while (points.size() < numPoints) {
        // Generate random spherical coordinates
        double theta = 2.0 * M_PI * rand() / RAND_MAX;
        double phi = M_PI * rand() / RAND_MAX;

        // Convert spherical coordinates to Cartesian coordinates
        double x = center.x + radius * sin(phi) * cos(theta);
        double y = center.y + radius * sin(phi) * sin(theta);
        double z = center.z + radius * cos(phi);

        // Create a new point
        //Point p = {static_cast<int>(x), static_cast<int>( y ), static_cast<int> (z)};
		Point p = {x,y,z};

        // Add the point to the set
        if (generatedPoints.count(p) == 0) {
            // Add the point to the set
            generatedPoints.insert(p);
		    points.push_back(p);
        }
    }

    // Convert the set of points to a vector and return
    //std::vector<Point> result(generatedPoints.begin(), generatedPoints.end());
    return points;
}

// Function to generate unique random points on the surface of a cube with a specified center
std::vector<Point3D> generateUniqueRandomPointsOnCube(int numPoints, double sideLength, const Point3D& center) {
    std::vector<Point3D> points;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(-0.5 * sideLength, 0.5 * sideLength);
    std::unordered_set<int> usedIndices;
	// change for w narrowing
	//std::unordered_set<std::size_t> usedIndices;

    while (points.size() < numPoints) {
        Point3D point;
        int face = std::uniform_int_distribution<int>(0, 5)(gen);
        
        if (face == 0) { // front face
            point.x = center.x + dis(gen);
            point.y = center.y + dis(gen);
            point.z = center.z + 0.5 * sideLength;
        } else if (face == 1) { // back face
            point.x = center.x + dis(gen);
            point.y = center.y + dis(gen);
            point.z = center.z - 0.5 * sideLength;
        } else if (face == 2) { // top face
            point.x = center.x + dis(gen);
            point.y = center.y + 0.5 * sideLength;
            point.z = center.z + dis(gen);
        } else if (face == 3) { // bottom face
            point.x = center.x + dis(gen);
            point.y = center.y - 0.5 * sideLength;
            point.z = center.z + dis(gen);
        } else if (face == 4) { // left face
            point.x = center.x - 0.5 * sideLength;
            point.y = center.y + dis(gen);
            point.z = center.z + dis(gen);
        } else { // right face
            point.x = center.x + 0.5 * sideLength;
            point.y = center.y + dis(gen);
            point.z = center.z + dis(gen);
        }

        // Generate a unique hash for the point
        int pointHash = static_cast<int>(point.x * 10000 + point.y * 100 + point.z);
		//change for w narrowing
		//int pointHash = static_cast<std::size_t>(point.x * 10000 + point.y * 100 + point.z);

        // Check if the point is already generated
        if (usedIndices.find(pointHash) == usedIndices.end()) {
            usedIndices.insert(pointHash);
            points.push_back(point);
        }
    }
    
    return points;
}

// Function to generate unique random points inside the cube cube with a specified center
std::vector<Point3D> generateAllPointsInsideCube(double sideLength, const Point3D& center) {
    std::vector<Point3D> points;

    for (double x = center.x - 0.5 * sideLength; x <= center.x + 0.5 * sideLength; x++) {
        for (double y = center.y - 0.5 * sideLength; y <= center.y + 0.5 * sideLength; y++) {
            for (double z = center.z - 0.5 * sideLength; z <= center.z + 0.5 * sideLength; z++) {
                //points.push_back({static_cast<int>(x), static_cast<int>( y ), static_cast<int> (z)});
				points.push_back({x,y,z});
            }
        }
    }

    return points;
}


//void seed_infection(Tissue &tissue, int time_step, GridPoint *grid_point) {
void seed_infection(Tissue &tissue, int time_step) {
  // _options->infection_coords contains the coords assigned just to rank_me()
  int64_t num_success_surface = 0;
  int64_t num_failure_surface = 0;
  
  int64_t num_success_inside = 0;
  int64_t num_failure_inside = 0;
  
  for (auto it = _options->infection_coords.begin(); it != _options->infection_coords.end(); it++) {
    auto infection_coords = *it;
	int64_t num_points = 60000;
	double sideLength = 97.0;

    if (infection_coords[3] == time_step) {
	  Point3D center = {infection_coords[0], infection_coords[1], infection_coords[2]};
	  
	  //std::vector<Point> points = generateAllPointsInSphere(center, 20.0, num_points);
	  std::vector<Point3D> points = generateUniqueRandomPointsOnCube(num_points, sideLength, center);
	  std::vector<Point3D> inside_points = generateAllPointsInsideCube(sideLength - 4.0, center);
	  // Process the generated points as needed
	  
	  
	  for (const auto& point : points) {
		// Do something with each point, such as setting initial infection
		std::cout << "Processing point (" << point.x << ", " << point.y << ", " << point.z << ")" << std::endl;
		// ...
		GridCoords coords({point.x, point.y,point.z});
		auto coords_1d = coords.to_1d();
		//change for w narrowing
		//std::size_t coords_1d = coords.to_1d();
		
		
		//while (true) {
        //GridCoords new_coords(coords_1d);
        if (tissue.set_initial_infection(coords_1d)) {
          WARN("Time step ", time_step, ": SUCCESSFUL initial infection at ", coords.str());
          //ostringstream oss;
          //oss << "rank " << rank_me() << " of " << rank_n() << " time step " << time_step
          //    << ": SUCCESSFUL initial infection at " << new_coords.str() << " after " << num_tries
          //    << " tries\n";
          //cerr << oss.str();
          //break;
		  num_success_surface++;
        }
		else{
		  num_failure_surface++;
		}
		
        coords_1d++;
	
        if (coords_1d >= get_num_grid_points()) {
          WARN("Could not find epicell to match uniform initial infection coord at ", coords.str());
          break;
        }
	  //}
	 }
	 
	   for (const auto& point : inside_points) {
		// Do something with each point, such as setting initial infection
		std::cout << "Processing point (" << point.x << ", " << point.y << ", " << point.z << ")" << std::endl;
		// ...
		GridCoords coords({point.x, point.y,point.z});
		auto coords_1d = coords.to_1d();
		//change for w narrowing
		//std::size_t coords_1d = coords.to_1d();
		
		//while (true) {
        //GridCoords new_coords(coords_1d);
        if (tissue.set_initial_infection(coords_1d)) {
          WARN("Time step ", time_step, ": SUCCESSFUL initial infection at ", coords.str());
          //ostringstream oss;
          //oss << "rank " << rank_me() << " of " << rank_n() << " time step " << time_step
          //    << ": SUCCESSFUL initial infection at " << new_coords.str() << " after " << num_tries
          //    << " tries\n";
          //cerr << oss.str();
          //break;
		  num_success_inside++;
        }
		else{
		  num_failure_inside++;
		}
		
        coords_1d++;
	
        if (coords_1d >= get_num_grid_points()) {
          WARN("Could not find epicell to match uniform initial infection coord at ", coords.str());
          break;
        }
	  //}
	 }
		//GridCoords coords({point.x, point.y,point.z});
		//int64_t grid_i = coords.to_1d();
		//tissue.make_epi_dead(grid_i);
                
		//GridPoint *grid_point = tissue.get_local_grid_point(tissue.grid_points, coords.to_1d());
		//assign air to these points
		//auto grid_point = tissue.get_grid_point_sample_data(coords.to_1d());
		
		//grid_point->epicell = EpiCellStatus::DEAD
		
	//}
	
      int64_t num_points_inside = inside_points.size();
	  int64_t total_inf_points = num_success_inside + num_success_surface;
      _options->infection_coords.erase(it--);
	  
	  WARN("total points_surface ", num_points);
	  WARN("total points_inside ", num_points_inside);
	  double suc_per_sur = ((double)num_success_surface/(double)num_points)*100;
	  double fail_per_sur =  ((double)num_failure_surface/(double)num_points)*100;
	  double suc_per_ins = ((double)num_success_inside/(double)num_points_inside)*100;
	  double fail_per_ins =  ((double)num_failure_inside/(double)num_points_inside)*100;
	  WARN("Success_surface ", num_success_surface, " and Failure_surface ", num_failure_surface);
	  WARN("Success_surface ", suc_per_sur, "% and Failure_surface ", fail_per_sur,"%");
	  WARN("Success_inside ", num_success_inside, " and Failure_inside ", num_failure_inside);
	  WARN("Success_inside ", suc_per_ins, "% and Failure_inside ", fail_per_ins,"%");
	  WARN("total infection point ", total_inf_points);
    }
  }
  
  barrier();
  tissue.add_new_actives(add_new_actives_timer);
  barrier();
}

void generate_tcells(Tissue &tissue, int time_step) {
  generate_tcell_timer.start();
  int local_num = _options->tcell_generation_rate / rank_n();
  int rem = _options->tcell_generation_rate - local_num * rank_n();
  if (rank_me() < rem) local_num++;
  if (time_step == 1) WARN("rem ", rem, " local num ", local_num, "\n");
  tissue.change_num_circulating_tcells(local_num);
#ifdef DEBUG
  auto all_num = reduce_one(local_num, op_fast_add, 0).wait();
  if (!rank_me() && all_num != _options->tcell_generation_rate)
    DIE("num tcells generated ", all_num, " != generation rate ", _options->tcell_generation_rate);
#endif
  generate_tcell_timer.stop();
}

int64_t get_rnd_coord(int64_t x, int64_t max_x) {
  int64_t new_x = x + _rnd_gen->get(0, 3) - 1;
  if (new_x < 0) new_x = 0;
  if (new_x >= max_x) new_x = max_x - 1;
  return new_x;
}

void update_circulating_tcells(int time_step, Tissue &tissue, double extravasate_fraction) {
  update_circulating_tcells_timer.start();
  auto num_circulating = tissue.get_num_circulating_tcells();
  // tcells prob of dying in vasculature is 1/vascular_period
  double portion_dying = (double)num_circulating / _options->tcell_vascular_period;
  int num_dying = floor(portion_dying);
  if (_rnd_gen->trial_success(portion_dying - num_dying)) num_dying++;
  tissue.change_num_circulating_tcells(-num_dying);
  _sim_stats.tcells_vasculature -= num_dying;
  num_circulating = tissue.get_num_circulating_tcells();
  double portion_xtravasing = extravasate_fraction * num_circulating;
  int num_xtravasing = floor(portion_xtravasing);
  if (_rnd_gen->trial_success(portion_xtravasing - num_xtravasing)) num_xtravasing++;
  for (int i = 0; i < num_xtravasing; i++) {
    progress();
    GridCoords coords(_rnd_gen);
    if (tissue.try_add_new_tissue_tcell(coords.to_1d())) {
      _sim_stats.tcells_tissue++;
      DBG(time_step, " tcell extravasates at ", coords.str(), "\n");
    }
  }
  _sim_stats.tcells_vasculature = num_circulating;
  update_circulating_tcells_timer.stop();
}

void update_tissue_tcell(int time_step, Tissue &tissue, GridPoint *grid_point, vector<int64_t> &nbs,
                         HASH_TABLE<int64_t, float> &chemokines_cache) {
  update_tcell_timer.start();
  TCell *tcell = grid_point->tcell;
  if (tcell->moved) {
    // don't update tcells that were added this time step
    tcell->moved = false;
    update_tcell_timer.stop();
    return;
  }
  tcell->tissue_time_steps--;
  if (tcell->tissue_time_steps == 0) {
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
    // not bound to an epicell - try to bind first with this cell then any one of the neighbors
    auto rnd_nbs = nbs;
    // include the current location
    rnd_nbs.push_back(grid_point->coords.to_1d());
    random_shuffle(rnd_nbs.begin(), rnd_nbs.end());
    for (auto &nb_grid_i : rnd_nbs) {
      DBG(time_step, " tcell ", tcell->id, " trying to bind at ", grid_point->coords.str(), "\n");
      auto nb_epicell_status = tissue.try_bind_tcell(nb_grid_i);
      bool bound = true;
      switch (nb_epicell_status) {
        case EpiCellStatus::EXPRESSING: _sim_stats.expressing--; break;
        case EpiCellStatus::INCUBATING: _sim_stats.incubating--; break;
        case EpiCellStatus::APOPTOTIC: _sim_stats.apoptotic--; break;
        default: bound = false;
      }
      if (bound) {
        DBG(time_step, " tcell ", tcell->id, " binds at ", grid_point->coords.str(), "\n");
        tcell->binding_period = _options->tcell_binding_period;
        _sim_stats.apoptotic++;
        break; //only allow tcell to bind to one cell!
      }
    }
  }
  if (tcell->binding_period == -1) {
    DBG(time_step, " tcell ", tcell->id, " trying to move at ", grid_point->coords.str(), "\n");
    // didn't bind - move on chemokine gradient or at random
    int64_t selected_grid_i = nbs[_rnd_gen->get(0, (int64_t)nbs.size())];
    // not bound - follow chemokine gradient
    float highest_chemokine = 0;
    if (_options->tcells_follow_gradient) {
      // get a randomly shuffled list of neighbors so the tcell doesn't always tend to move in the
      // same direction when there is a chemokine gradient
      auto rnd_nbs = nbs;
      random_shuffle(rnd_nbs.begin(), rnd_nbs.end());
      for (auto nb_grid_i : rnd_nbs) {
        float chemokine = 0;
        auto it = chemokines_cache.find(nb_grid_i);
        if (it == chemokines_cache.end()) {
          chemokine = tissue.get_chemokine(nb_grid_i);
          chemokines_cache.insert({nb_grid_i, chemokine});
        } else {
          chemokine = it->second;
        }
        if (chemokine > highest_chemokine) {
          highest_chemokine = chemokine;
          selected_grid_i = nb_grid_i;
        }
        DBG(time_step, " tcell ", tcell->id, " found nb chemokine ", chemokine, " at ",
            GridCoords(selected_grid_i).str(), "\n");
      }
    }
    if (highest_chemokine == 0) {
      // no chemokines found - move randomly
      auto rnd_nb_i = _rnd_gen->get(0, (int64_t)nbs.size());
      selected_grid_i = nbs[rnd_nb_i];
      DBG(time_step, " tcell ", tcell->id, " try random move to ",
          GridCoords(selected_grid_i).str(), "\n");
    } else {
      DBG(time_step, " tcell ", tcell->id, " - highest chemokine at ",
          GridCoords(selected_grid_i).str(), "\n");
    }
    // try a few times to find an open spot
    for (int i = 0; i < 5; i++) {
      if (tissue.try_add_tissue_tcell(selected_grid_i, *tcell)) {
        DBG(time_step, " tcell ", tcell->id, " at ", grid_point->coords.str(), " moves to ",
            GridCoords(selected_grid_i).str(), "\n");
        delete grid_point->tcell;
        grid_point->tcell = nullptr;
        break;
      }
      // choose another location at random
      auto rnd_nb_i = _rnd_gen->get(0, (int64_t)nbs.size());
      selected_grid_i = nbs[rnd_nb_i];
      DBG(time_step, " tcell ", tcell->id, " try random move to ",
          GridCoords(selected_grid_i).str(), "\n");
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
    DBG(time_step, " epicell ", grid_point->epicell->str(), "\n");
  bool produce_virions = false;
  switch (grid_point->epicell->status) {
    case EpiCellStatus::HEALTHY: {
      if (grid_point->virions > 0) {
        double local_infectivity = _options->infectivity;
        if (grid_point->chemokine > 0) {
          local_infectivity *= _options->infectivity_multiplier;
        }
        if (_rnd_gen->trial_success(local_infectivity * grid_point->virions)) {
          grid_point->epicell->infect();
          _sim_stats.incubating++;
        }
      }

	   
	  
      break;
    }
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
      if (grid_point->epicell->apoptosis_death()) {
        _sim_stats.dead++;
        _sim_stats.apoptotic--;
      } else if (grid_point->epicell->was_expressing()) {
        produce_virions = true;
      }
      break;
    default: break;
  }
  if (produce_virions) {
    if (grid_point->chemokine > 0) {
      grid_point->virions += _options->virion_production * _options->virion_production_multiplier;
    } else {
      grid_point->virions += _options->virion_production;
    }
    grid_point->chemokine = min(grid_point->chemokine + _options->chemokine_production, 1.0);
  }
  update_epicell_timer.stop();
}

void update_chemokines(GridPoint *grid_point, vector<int64_t> &nbs,
                       HASH_TABLE<int64_t, float> &chemokines_to_update) {
  update_concentration_timer.start();
  // Concentrations diffuse, i.e. the concentration at any single grid point tends to the average
  // of all the neighbors. So here we tell each neighbor what the current concentration is and
  // later those neighbors will compute their own averages. We do it in this "push" manner because
  // then we don't need to check the neighbors from every single grid point, but just push from
  // ones with concentrations > 0 (i.e. active grid points)
  if (grid_point->chemokine > 0) {
    grid_point->chemokine *= (1.0 - _options->chemokine_decay_rate);
    if (grid_point->chemokine < _options->min_chemokine) grid_point->chemokine = 0;
  }
  if (grid_point->chemokine > 0) {
    for (auto &nb_grid_i : nbs) {
      chemokines_to_update[nb_grid_i] += grid_point->chemokine;
    }
  }
  update_concentration_timer.stop();
}

void update_virions(GridPoint *grid_point, vector<int64_t> &nbs,
                    HASH_TABLE<int64_t, float> &virions_to_update) {
  update_concentration_timer.start();
  grid_point->virions = grid_point->virions * (1.0 - _options->virion_clearance_rate);
  assert(grid_point->virions >= 0);
  if (grid_point->virions > 0) {
    for (auto &nb_grid_i : nbs) {
      virions_to_update[nb_grid_i] += grid_point->virions;
    }
  }
  update_concentration_timer.stop();
}

void diffuse(float &conc, float &nb_conc, double diffusion, int num_nbs) {
  // set to be average of neighbors plus self
  // amount that diffuses
  float conc_diffused = diffusion * conc;
  // average out diffused amount across all neighbors
  float conc_per_point = (conc_diffused + diffusion * nb_conc) / (num_nbs + 1);
  conc = conc - conc_diffused + conc_per_point;
  if (conc > 1.0) conc = 1.0;
  if (conc < 0) DIE("conc < 0: ", conc, " diffused ", conc_diffused, " pp ", conc_per_point);
  nb_conc = 0;
}

void spread_virions(float &virions, float &nb_virions, double diffusion, int num_nbs) {
  float virions_diffused = virions * diffusion;
  float virions_left = virions - virions_diffused;
  float avg_nb_virions = (virions_diffused + nb_virions * diffusion) / (num_nbs + 1);
  virions = virions_left + avg_nb_virions;
  nb_virions = 0;
}

void set_active_grid_points(Tissue &tissue) {
  set_active_points_timer.start();
  vector<GridPoint *> to_erase = {};
  // iterate through all active local grid points and set changes
  for (auto grid_point = tissue.get_first_active_grid_point(); grid_point;
       grid_point = tissue.get_next_active_grid_point()) {
    auto nbs = tissue.get_neighbors(grid_point->coords);
    //diffuse(grid_point->chemokine, grid_point->nb_chemokine, _options->chemokine_diffusion_coef,
            //nbs->size());
	if (grid_point->epicell->type == EpiCellType::TYPE2 || grid_point->epicell->type == EpiCellType::TYPE1 ){
		spread_virions(grid_point->virions, grid_point->nb_virions, _options->virion_diffusion_coef_infectable,
                   nbs->size());
	    diffuse(grid_point->chemokine, grid_point->nb_chemokine, _options->chemokine_diffusion_coef_infectable,
            nbs->size());
	}
	else if (grid_point->epicell->type == EpiCellType::INTERSTITIAL){
		spread_virions(grid_point->virions, grid_point->nb_virions, _options->virion_diffusion_coef_interstitial,
                   nbs->size());
	    diffuse(grid_point->chemokine, grid_point->nb_chemokine, _options->chemokine_diffusion_coef_interstitial,
            nbs->size());
	}
	else{
		spread_virions(grid_point->virions, grid_point->nb_virions, _options->virion_diffusion_coef_air,
                   nbs->size());
		diffuse(grid_point->chemokine, grid_point->nb_chemokine, _options->chemokine_diffusion_coef_air,
            nbs->size());
	}				   
    if (grid_point->chemokine < _options->min_chemokine) grid_point->chemokine = 0;
    if (grid_point->virions > MAX_VIRIONS) grid_point->virions = MAX_VIRIONS;
    if (grid_point->virions < MIN_VIRIONS) grid_point->virions = 0;
    if (grid_point->tcell) grid_point->tcell->moved = false;
    _sim_stats.virions += grid_point->virions;
    if (!grid_point->is_active()) to_erase.push_back(grid_point);
  }
  for (auto grid_point : to_erase) tissue.erase_active(grid_point);
  set_active_points_timer.stop();
}

void sample(int time_step, vector<SampleData> &samples, int64_t start_id, ViewObject view_object) {
  char cwd_buf[MAX_FILE_PATH];
  string fname = _options->output_dir + "/samples/sample_" + view_object_str(view_object) + "_" +
                 to_string(time_step) + ".vtk";
  int x_dim = _options->dimensions[0] / _options->sample_resolution;
  int y_dim = _options->dimensions[1] / _options->sample_resolution;
  int z_dim = _options->dimensions[2] / _options->sample_resolution;
  if (z_dim == 0) z_dim = 1;
  size_t tot_sz = x_dim * y_dim * z_dim;
  int spacing = 15 * _options->sample_resolution;
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
             << "SPACING " << spacing << " " << spacing << " " << spacing << "\n"
             << "ORIGIN 0 0 0\n"
             << "CELL_DATA " << (x_dim * y_dim * z_dim) << "\n"
             << "SCALARS ";
  switch (view_object) {
    case ViewObject::VIRUS: header_oss << "virus"; break;
    case ViewObject::TCELL_TISSUE: header_oss << "t-cell-tissue"; break;
    case ViewObject::EPICELL: header_oss << "epicell"; break;
    case ViewObject::CHEMOKINE: header_oss << "chemokine"; break;
    case ViewObject::INFLAM_SIGNAL_CELL: header_oss << "inflamsignalcell"; break;
    default: SDIE("unknown view object");
  }
  header_oss << " unsigned_char\n"
             << "LOOKUP_TABLE default\n";
  auto header_str = header_oss.str();
  if (!rank_me()) {
    tot_sz += header_str.size();
    // rank 0 creates the file and truncates it to the correct length
    auto fileno = open(fname.c_str(), O_WRONLY | O_CREAT, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
    if (fileno == -1) SDIE("Cannot create file ", fname, ": ", strerror(errno), "\n");
    if (ftruncate(fileno, tot_sz) == -1)
      DIE("Could not truncate ", fname, " to ", tot_sz, " bytes\n");
    close(fileno);
  }
  // wait until rank 0 has finished setting up the file
  upcxx::barrier();
  auto fileno = open(fname.c_str(), O_WRONLY, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
  if (fileno == -1) DIE("Cannot open file ", fname, ": ", strerror(errno), "\n");
  if (!upcxx::rank_me()) {
    size_t bytes_written = pwrite(fileno, header_str.c_str(), header_str.length(), 0);
    if (bytes_written != header_str.length())
      DIE("Could not write all ", header_str.length(), " bytes: only wrote ", bytes_written);
  }
  // each rank writes one portion of the dataset to the file
  unsigned char *buf = new unsigned char[samples.size()];
  double chemo_scale = 255.0 / log(1.0 / _options->min_chemokine);
  double virion_scale = 255.0 / log(MAX_VIRIONS);
  // DBG(time_step, " writing data from ", start_id, " to ", start_id + samples.size(), "\n");
  for (int64_t i = 0; i < samples.size(); i++) {
    auto &sample = samples[i];
    unsigned char val = 0;
    double scaled_chemo = 0;
    switch (view_object) {
      case ViewObject::TCELL_TISSUE:
        assert(sample.tcells >= 0);
        if (sample.tcells > 0.5)
          val = 4;
        else if (sample.tcells > 0.25)
          val = 3;
        else if (sample.tcells > 0.125)
          val = 2;
        else if (sample.tcells > 0)
          val = 1;
        break;
      case ViewObject::EPICELL:
        if (sample.has_epicell) val = static_cast<unsigned char>(sample.epicell_status) + 1;
        break;
      case ViewObject::VIRUS:
        assert(sample.virions >= 0);
        if (sample.virions > 1) val = virion_scale * log(sample.virions);
        if (sample.virions > 0 && val == 0) val = 1;
        break;
      case ViewObject::CHEMOKINE:
        assert(sample.chemokine >= 0 && sample.chemokine <= 1);
        // set chemokine to 0 to ensure we can see the tcells
        if (sample.tcells > 0) break;
        scaled_chemo = sample.chemokine / _options->min_chemokine;
        if (scaled_chemo > 1) val = chemo_scale * log(scaled_chemo);
        if (sample.chemokine > 0 && val == 0) val = 1;
        break;
      case ViewObject::INFLAM_SIGNAL_CELL:
        //if (sample.has_inflam_signal_cell) val = static_cast<unsigned char>(1);
		if (sample.has_inflam_signal_cell && sample.has_air) val = static_cast<unsigned char>(1);
		else if (sample.has_inflam_signal_cell && sample.has_interstitial) val = static_cast<unsigned char>(2);
		else if (sample.has_inflam_signal_cell && sample.has_infected) val = static_cast<unsigned char>(3);
        break;
      default: break;
    }
    buf[i] = val;
  }
  size_t fpos = header_str.length() + start_id;
  sample_write_timer.start();
  auto bytes_written = pwrite(fileno, buf, samples.size(), fpos);
  if (bytes_written != samples.size())
    DIE("Could not write all ", samples.size(), " bytes; only wrote ", bytes_written, "\n");
  delete[] buf;
  close(fileno);
  sample_write_timer.stop();
  upcxx::barrier();
}

int64_t get_samples(Tissue &tissue, vector<SampleData> &samples) {
  int64_t num_points =
      get_num_grid_points() / (_options->sample_resolution * _options->sample_resolution);
  if (_grid_size->z > 1) num_points /= _options->sample_resolution;
  int64_t num_points_per_rank = ceil((double)num_points / rank_n());
  int64_t start_id = rank_me() * num_points_per_rank;
  int64_t end_id = min((rank_me() + 1) * num_points_per_rank, num_points);
  samples.clear();
  if (end_id > start_id) {
    samples.reserve(end_id - start_id);
    int64_t i = 0;
    bool done = false;
    int block_size = _options->sample_resolution * _options->sample_resolution;
    if (_grid_size->z > 1) block_size *= _options->sample_resolution;
    vector<SampleData> block_samples;
    for (int x = 0; x < _grid_size->x && !done; x += _options->sample_resolution) {
      for (int y = 0; y < _grid_size->y && !done; y += _options->sample_resolution) {
        for (int z = 0; z < _grid_size->z; z += _options->sample_resolution) {
          if (i >= end_id) {
            done = true;
            break;
          }
          if (i >= start_id) {
            progress();
#ifdef AVERAGE_SUBSAMPLE
            float virions = 0;
            float chemokine = 0;
            int num_tcells = 0;
            bool epicell_found = false;
            bool inflam_signal_found = false;
			bool air_found = false;
			bool interstitial_found = false;
			bool infected_found = false;
            array<int, 5> epicell_counts{0};
            block_samples.clear();
            bool done_sub = false;
            for (int subx = x; subx < x + _options->sample_resolution; subx++) {
              if (subx >= _grid_size->x) break;
              for (int suby = y; suby < y + _options->sample_resolution; suby++) {
                if (suby >= _grid_size->y) break;
                for (int subz = z; subz < z + _options->sample_resolution; subz++) {
                  if (subz >= _grid_size->z) break;
                  auto sub_sd =
                      tissue.get_grid_point_sample_data(GridCoords::to_1d(subx, suby, subz));
                  num_tcells += sub_sd.tcells;
                  if (sub_sd.has_epicell) {
                    epicell_found = true;
                    switch (sub_sd.epicell_status) {
                      case EpiCellStatus::HEALTHY: epicell_counts[0]++; break;
                      case EpiCellStatus::INCUBATING: epicell_counts[1]++; break;
                      case EpiCellStatus::EXPRESSING: epicell_counts[2]++; break;
                      case EpiCellStatus::APOPTOTIC: epicell_counts[3]++; break;
                      case EpiCellStatus::DEAD: epicell_counts[4]++; break;
					  //case EpiCellStatus::TYPE2: epicell_counts[5]++; break;
					  //case EpiCellStatus::TYPE1: epicell_counts[6]++; break;
					  //case EpiCellStatus::AIR: epicell_counts[7]++; break;
					  //case EpiCellStatus::INTERSTITIAL: epicell_counts[8]++; break;
                    }
                  }
                  if (sub_sd.has_inflam_signal_cell) {
                    inflam_signal_found = true;
					if (sub_sd.has_air)
						air_found = true;
					if (sub_sd.has_interstitial)
						interstitial_found = true;
					if (sub_sd.has_infected)
						infected_found = true;
                  }
				  
                  chemokine += sub_sd.chemokine;
                  virions += sub_sd.virions;
                }
              }
            }
            EpiCellStatus epi_status = EpiCellStatus::HEALTHY;
            if (epicell_found) {
              // chose the epicell status supported by the majority of grid points
              int max_epicell_i = 0, max_count = 0;
              for (int j = 0; j < 5; j++) {
                if (max_count < epicell_counts[j]) {
                  max_count = epicell_counts[j];
                  max_epicell_i = j;
                }
              }
              switch (max_epicell_i) {
                case 0: epi_status = EpiCellStatus::HEALTHY; break;
                case 1: epi_status = EpiCellStatus::INCUBATING; break;
                case 2: epi_status = EpiCellStatus::EXPRESSING; break;
                case 3: epi_status = EpiCellStatus::APOPTOTIC; break;
                case 4: epi_status = EpiCellStatus::DEAD; break;
				//case 5: epi_status = EpiCellStatus::TYPE2; break;
				//case 6: epi_status = EpiCellStatus::TYPE1; break;
				//case 7: epi_status = EpiCellStatus::AIR; break;
				//case 8: epi_status = EpiCellStatus::INTERSTITIAL; break;
				
              }
            }
            SampleData sd = {.tcells = (double)num_tcells / block_size,
                             .has_epicell = epicell_found,
                             .has_inflam_signal_cell = inflam_signal_found,
							 .has_air = air_found,
							 .has_interstitial = interstitial_found,
							 .has_infected = infected_found,
                             .epicell_status = epi_status,
                             .virions = virions / block_size,
                             .chemokine = chemokine / block_size};
#else
            auto sd = tissue.get_grid_point_sample_data(GridCoords::to_1d(x, y, z));
#endif
            samples.push_back(sd);
          }
          i++;
        }
        if (done) break;
      }
      if (done) break;
    }
  }
  barrier();
  auto samples_written = reduce_one(samples.size(), op_fast_add, 0).wait();
  if (num_points != samples_written)
    SWARN("Number of point ", num_points, " != ", samples_written, " samples written");
  SLOG_VERBOSE("Number of samples written ", samples_written, "\n");
  return start_id;
}

void run_sim(Tissue &tissue) {
  BarrierTimer timer(__FILEFUNC__);

  auto start_t = NOW();
  auto curr_t = start_t;
  // TODO Allow for 1 timestep
  auto five_perc = (_options->num_timesteps >= 50) ? _options->num_timesteps / 50 : 1;
  _sim_stats.init();
  int64_t whole_lung_volume = (int64_t)_options->whole_lung_dims[0] *
                              (int64_t)_options->whole_lung_dims[1] *
                              (int64_t)_options->whole_lung_dims[2];
  auto sim_volume = get_num_grid_points();
  double extravasate_fraction = (double)sim_volume / whole_lung_volume;
  SLOG("Fraction of circulating T cells extravasating is ", extravasate_fraction, "\n");
  SLOG("# datetime                    step    ", _sim_stats.header(STATS_COL_WIDTH),
       "<%active  lbln>\n");
  // store the total concentration increment updates for target grid points
  // chemokine, virions
  HASH_TABLE<int64_t, float> chemokines_to_update;
  HASH_TABLE<int64_t, float> chemokines_cache;
  HASH_TABLE<int64_t, float> virions_to_update;
  bool warned_boundary = false;
  vector<SampleData> samples;
  for (int time_step = 0; time_step < _options->num_timesteps; time_step++) {
    DBG("Time step ", time_step, "\n");
    //Change
	auto grid_point = tissue.get_first_active_grid_point();
	//seed_infection(tissue, time_step, grid_point);
	seed_infection(tissue, time_step);
    barrier();
    if (time_step == _options->antibody_period)
      _options->virion_clearance_rate *= _options->antibody_factor;
    chemokines_to_update.clear();
    virions_to_update.clear();
    chemokines_cache.clear();
    if (time_step > _options->tcell_initial_delay) {
      generate_tcells(tissue, time_step);
      barrier();
    }
    compute_updates_timer.start();
    update_circulating_tcells(time_step, tissue, extravasate_fraction);
    // iterate through all active local grid points and update
    for (auto grid_point = tissue.get_first_active_grid_point(); grid_point;
         grid_point = tissue.get_next_active_grid_point()) {
      if (grid_point->chemokine > 0)
        DBG("chemokine\t", time_step, "\t", grid_point->coords.x, "\t", grid_point->coords.y, "\t",
            grid_point->coords.z, "\t", grid_point->chemokine, "\n");
      if (grid_point->virions > 0)
        DBG("virions\t", time_step, "\t", grid_point->coords.x, "\t", grid_point->coords.y, "\t",
            grid_point->coords.z, "\t", grid_point->virions, "\n");
      if (!warned_boundary && grid_point->epicell &&
          grid_point->epicell->status != EpiCellStatus::HEALTHY) {
        if (!grid_point->coords.x || grid_point->coords.x == _grid_size->x - 1 ||
            !grid_point->coords.y || grid_point->coords.y == _grid_size->y - 1 ||
            (_grid_size->z > 1 &&
             (!grid_point->coords.z || grid_point->coords.z == _grid_size->z - 1))) {
          WARN("Hit boundary at ", grid_point->coords.str(), " ", grid_point->epicell->str(),
               " virions ", grid_point->virions, " chemokine ", grid_point->chemokine);
          warned_boundary = true;
        }
      }
      // DBG("updating grid point ", grid_point->str(), "\n");
      upcxx::progress();
      auto nbs = tissue.get_neighbors(grid_point->coords);
      // the tcells are moved (added to the new list, but only cleared out at the end of all
      // updates)
      if (grid_point->tcell)
        update_tissue_tcell(time_step, tissue, grid_point, *nbs, chemokines_cache);
      if (grid_point->epicell) update_epicell(time_step, tissue, grid_point);
	 
	  
      if (!grid_point->inflam_signal_all && grid_point->chemokine > 0) {
        grid_point->inflam_signal_all = true;
        _sim_stats.num_inflam_signal_all++;
		if (grid_point->epicell && grid_point->epicell->type == EpiCellType::AIR){
		   grid_point->air = true;
		  _sim_stats.num_aircell_w_inflam++;
	    }
		if (grid_point->epicell && grid_point->epicell->type == EpiCellType::TYPE1){
		  grid_point->interstitial = true;
		  _sim_stats.num_intertsitial_w_inflam++;
		}
        if (grid_point->epicell && grid_point->epicell->type == EpiCellType::TYPE2){
		  grid_point->infected = true;
		  _sim_stats.num_infected_w_inflam++;
		}
        if (grid_point->epicell && grid_point->epicell->type != EpiCellType::AIR) {
          grid_point->inflam_signal_cell = true;
          _sim_stats.num_inflam_signal_cell++;
        }

      }
	  
	   //added
	  if (grid_point->epicell && grid_point->epicell->type == EpiCellType::AIR)
		_sim_stats.num_aircell++;
	  if (grid_point->epicell && grid_point->epicell->type == EpiCellType::TYPE1)
		_sim_stats.num_intertsitial++;
      if (grid_point->epicell && grid_point->epicell->type == EpiCellType::TYPE2)
		_sim_stats.num_infected++;
      update_chemokines(grid_point, *nbs, chemokines_to_update);
      update_virions(grid_point, *nbs, virions_to_update);
      if (grid_point->is_active()) tissue.set_active(grid_point);
    }
    barrier();
    compute_updates_timer.stop();
    tissue.accumulate_chemokines(chemokines_to_update, accumulate_concentrations_timer);
    tissue.accumulate_virions(virions_to_update, accumulate_concentrations_timer);
    barrier();
    if (time_step % five_perc == 0 || time_step == _options->num_timesteps - 1) {
      auto num_actives = reduce_one(tissue.get_num_actives(), op_fast_add, 0).wait();
      auto perc_actives = 100.0 * num_actives / get_num_grid_points();
      auto max_actives = reduce_one(tissue.get_num_actives(), op_fast_max, 0).wait();
      auto load_balance = max_actives ? (double)num_actives / rank_n() / max_actives : 1;
      chrono::duration<double> t_elapsed = NOW() - curr_t;
      curr_t = NOW();
      SLOG("[", get_current_time(), " ", setprecision(2), fixed, setw(7), right, t_elapsed.count(),
           "s]: ", setw(8), left, time_step, _sim_stats.to_str(STATS_COL_WIDTH), setprecision(3),
           fixed, "< ", perc_actives, " ", load_balance, " >\n");
    }
    barrier();
    tissue.add_new_actives(add_new_actives_timer);
    barrier();

    _sim_stats.virions = 0;

    set_active_grid_points(tissue);
    barrier();

    if (_options->sample_period > 0 &&
        (time_step % _options->sample_period == 0 || time_step == _options->num_timesteps - 1)) {
      sample_timer.start();
      samples.clear();
      int64_t start_id = get_samples(tissue, samples);
      sample(time_step, samples, start_id, ViewObject::EPICELL);
      sample(time_step, samples, start_id, ViewObject::TCELL_TISSUE);
      sample(time_step, samples, start_id, ViewObject::VIRUS);
      sample(time_step, samples, start_id, ViewObject::CHEMOKINE);
      sample(time_step, samples, start_id, ViewObject::INFLAM_SIGNAL_CELL);
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
  sample_write_timer.done_all();
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
  if (pin_thread(getpid(), local_team().rank_me()) == -1)
    WARN("Could not pin process ", getpid(), " to core ", rank_me());
  else
    SLOG_VERBOSE("Pinned processes, with process 0 (pid ", getpid(), ") pinned to core ",
                 local_team().rank_me(), "\n");
#ifdef BLOCK_PARTITION
  SLOG_VERBOSE("Using block partitioning\n");
#else
  SLOG_VERBOSE("Using linear partitioning\n");
#endif
#ifndef AVERAGE_SUBSAMPLE
  SLOG_VERBOSE("Not computing averages of subsamples\n");
#endif
  MemoryTrackerThread memory_tracker;
  memory_tracker.start();
  auto start_free_mem = get_free_mem();
  SLOG(KBLUE, "Starting with ", get_size_str(start_free_mem), " free on node 0", KNORM, "\n");
  Tissue tissue;
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
