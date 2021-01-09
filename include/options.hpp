#pragma once

#include <sys/stat.h>

#include <iostream>
#include <regex>
#include <upcxx/upcxx.hpp>

#include "CLI11.hpp"
#include "version.h"
#include "utils.hpp"

using std::array;
using std::cout;
using std::endl;
using std::vector;

#include "upcxx_utils/log.hpp"
#include "upcxx_utils/timers.hpp"

using namespace upcxx_utils;

#define YES_NO(X) ((X) ? "YES" : "NO")

class Options {
  vector<string> splitter(string in_pattern, string &content) {
    vector<string> split_content;
    std::regex pattern(in_pattern);
    copy(std::sregex_token_iterator(content.begin(), content.end(), pattern, -1),
         std::sregex_token_iterator(), back_inserter(split_content));
    return split_content;
  }

  template <typename T>
  string vec_to_str(const vector<T> &vec, const string &delimiter = ",") {
    std::ostringstream oss;
    for (auto elem : vec) {
      oss << elem;
      if (elem != vec.back()) oss << delimiter;
    }
    return oss.str();
  }

  void setup_output_dir() {
    if (!upcxx::rank_me()) {
      // create the output directory and stripe it
      if (mkdir(output_dir.c_str(), S_IRWXU) == -1) {
        // could not create the directory
        if (errno == EEXIST) {
          cerr << KLRED << "WARNING: " << KNORM << "Output directory " << output_dir
               << " already exists. May overwrite existing files\n";
        } else {
          ostringstream oss;
          oss << KLRED << "ERROR: " << KNORM << " Could not create output directory " << output_dir
              << ": " << strerror(errno) << endl;
          throw std::runtime_error(oss.str());
        }
      } else {
        // created the directory - now stripe it if possible
        auto status = std::system("which lfs 2>&1 > /dev/null");
        if (WIFEXITED(status) && WEXITSTATUS(status) == 0) {
          string cmd = "lfs setstripe -c -1 " + output_dir;
          auto status = std::system(cmd.c_str());
          if (WIFEXITED(status) && WEXITSTATUS(status) == 0)
            cout << "Set Lustre striping on the output directory\n";
          else
            cout << "Failed to set Lustre striping on output directory: " << WEXITSTATUS(status)
                 << endl;
        }
      }
    }
    upcxx::barrier();
    // now setup a samples subdirectory
    if (!upcxx::rank_me()) {
      // create the output directory and stripe it
      string samples_dir = output_dir + "/samples";
      if (mkdir(samples_dir.c_str(), S_IRWXU) == -1) {
        // could not create the directory
        if (errno == EEXIST) {
          cerr << KLRED << "WARNING: " << KNORM
               << "Samples directory already exists. May overwrite existing files\n";
        } else {
          ostringstream oss;
          oss << KLRED << "ERROR: " << KNORM
              << " Could not create samples directory : " << strerror(errno) << endl;
          throw std::runtime_error(oss.str());
        }
      }
    }
    upcxx::barrier();
  }

  void setup_log_file() {
    if (!upcxx::rank_me()) {
      string log_fname = output_dir + "/simcov.log";
      // check to see if simcov.log exists. If so, rename it
      if (file_exists(log_fname)) {
        string new_log_fname = output_dir + "/simcov-" + get_current_time(true) + ".log";
        cerr << KLRED << "WARNING: " << KNORM << log_fname << " exists. Renaming to "
             << new_log_fname << endl;
        if (rename(log_fname.c_str(), new_log_fname.c_str()) == -1)
          DIE("Could not rename ", log_fname, ": ", strerror(errno));
      }
    }
    upcxx::barrier();
  }

  void set_random_infections(int num) {
    for (int i = 0; i < num; i++) {
      if (i % rank_n() != rank_me()) continue;
      infection_coords.push_back({_rnd_gen->get(0, dimensions[0]), _rnd_gen->get(0, dimensions[1]),
                                  _rnd_gen->get(0, dimensions[2]), 0});
    }
  }

  void set_uniform_infections(int num) { SDIE("Not yet implemented"); }

  bool parse_infection_coords(vector<string> &coords_strs) {
    auto get_locations_count = [](const string &s, const string &name) -> int {
      int num = 0;
      if (s.compare(0, name.length(), name) == 0) {
        string num_str = s.substr(name.length());
        try {
          num = std::stoi(num_str);
        } catch (std::invalid_argument arg) {
          num = 0;
        }
        if (num < 1) return 0;
      }
      return num;
    };

    if (coords_strs.size() == 1) {
      int num = get_locations_count(coords_strs[0], "random:");
      if (num > 0) {
        set_random_infections(num);
        return true;
      }
      num = get_locations_count(coords_strs[0], "uniform:");
      if (num > 0) {
        set_uniform_infections(num);
        return true;
      }
    }
    for (int i = 0; i < coords_strs.size(); i++) {
      if (i % rank_n() != rank_me()) continue;
      auto coords_and_time = splitter(",", coords_strs[i]);
      if (coords_and_time.size() == 4) {
        try {
          infection_coords.push_back({std::stoi(coords_and_time[0]), std::stoi(coords_and_time[1]),
                                      std::stoi(coords_and_time[2]),
                                      std::stoi(coords_and_time[3])});
        } catch (std::invalid_argument arg) {
          coords_and_time.clear();
        }
      }
      if (coords_and_time.size() != 4) {
        ostringstream oss;
        oss << KLRED << "ERROR: " << KNORM << "incorrect specification of infection coords in "
            << "string \"" << coords_strs[i] << "\"\n"
            << " - should be four comma-separated (x,y,z,t) values or random:N or uniform:N\n";
        cerr << oss.str();
        return false;
      }
    }
    return true;
  }

 public:
  vector<int> dimensions{100, 100, 1};
  vector<int> whole_lung_dims{240, 120, 60};
  // each time step should be about 1 minute, so one day = 1440 time steps
  int num_timesteps = 2000;

  // x,y,z location and timestep
  vector<array<int, 4>> infection_coords;
  int initial_infection = 1000;
  int infectable_spacing = 1;

  // these periods are normally distributed with mean and stddev
  int incubation_period = 100;
  int apoptosis_period = 30;
  int expressing_period = 100;

  int tcell_generation_rate = 100;
  int tcell_initial_delay = 1000;
  int tcell_vascular_period = 600;
  int tcell_tissue_period = 60;
  int tcell_binding_period = 5;
  double max_binding_prob = 0.1;

  double infectivity = 0.0001;
  int virion_production = 100;
  double virion_decay_rate = 0.1;
  double virion_diffusion_coef = 0.5;

  double chemokine_production = 0.5;
  double chemokine_decay_rate = 0.01;
  double chemokine_diffusion_coef = 0.2;
  double min_chemokine = 0.001;

  double antibody_factor = 1;
  int antibody_period = 900;

  unsigned rnd_seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  string output_dir = "simcov-run-n" + to_string(upcxx::rank_n()) + "-N" +
                      to_string(upcxx::rank_n() / upcxx::local_team().rank_n()) + "-" +
                      get_current_time(true);
  int sample_period = 1;
  int sample_resolution = 1;
  int max_block_dim = 10;

  bool tcells_follow_gradient = false;

  bool show_progress = false;
  bool verbose = false;

  bool load(int argc, char **argv) {
    // SIMCOV version v0.1-a0decc6-master (Release) built on 2020-04-08T22:15:40 with g++
    string full_version_str = "SimCov version " + string(SIMCOV_VERSION) + "-" +
                              string(SIMCOV_BRANCH) + " built on " + string(SIMCOV_BUILD_DATE);
    vector<string> infection_coords_strs;

    CLI::App app(full_version_str);
    app.add_option("-d,--dim", dimensions, "Dimensions: x y z")
        ->delimiter(',')
        ->expected(3)
        ->capture_default_str();
    app.add_option("--whole-lung-dim", whole_lung_dims, "Whole lung dimensions: x y z")
        ->delimiter(',')
        ->expected(3)
        ->capture_default_str();
    app.add_option("-t,--timesteps", num_timesteps, "Number of timesteps")
        ->check(CLI::Range(1, 1000000))
        ->capture_default_str();
    app.add_option(
           "--infection-coords", infection_coords_strs,
           "Location of multiple initial infections, of form \"x1,y1,z1,t1 x2,y2,z2,t2...\"\n"
           "where x,y,z are grid coords and t is a timestep, or\n"
           "\"uniform:N\" for N uniformly distributed points or\n"
           "\"random:N\" for N randomly distributed points")
        ->delimiter(' ')
        ->capture_default_str();
    app.add_option("--initial-infection", initial_infection,
                   "Number of virions at initial infection locations")
        ->capture_default_str();
    app.add_option("--infectable-spacing", infectable_spacing,
                   "Number of grid points to space out infectable cells")
        ->capture_default_str();
    app.add_option("--incubation-period", incubation_period,
                   "Average number of time steps to expressing virions after cell is infected")
        ->capture_default_str();
    app.add_option("--apoptosis-period", apoptosis_period,
                   "Average number of time steps to death after apoptosis is induced")
        ->capture_default_str();
    app.add_option("--expressing-period", expressing_period,
                   "Average number of time steps to death after a cell starts expresssing")
        ->capture_default_str();
    app.add_option("--infectivity", infectivity,
                   "Factor multiplied by number of virions to determine probability of infection")
        ->check(CLI::Range(0.0, 1.0))
        ->capture_default_str();
    app.add_option("--virion-production", virion_production,
                   "Number of virions produced by expressing cell each time step")
        ->capture_default_str();
    app.add_option("--virion-decay", virion_decay_rate,
                   "Fraction by which virion count drops each time step")
        ->check(CLI::Range(0.0, 1.0))
        ->capture_default_str();
    app.add_option("--virion-diffusion", virion_diffusion_coef,
                   "Fraction of virions that diffuse into all neighbors each time step")
        ->check(CLI::Range(0.0, 1.0))
        ->capture_default_str();
    app.add_option("--chemokine-production", chemokine_production,
                   "Amount of chemokine produced by expressing cells each time step")
        ->check(CLI::Range(0.0, 1.0))
        ->capture_default_str();
    app.add_option("--chemokine-decay", chemokine_decay_rate,
                   "Amount by which chemokine concentration drops each time step")
        ->check(CLI::Range(0.0, 1.0))
        ->capture_default_str();
    app.add_option("--chemokine-diffusion", chemokine_diffusion_coef,
                   "Fraction of chemokine concentration that diffuses into all neighbors "
                   "each time step")
        ->check(CLI::Range(0.0, 1.0))
        ->capture_default_str();
    app.add_option("--min-chemokine", min_chemokine,
                   "Minimum chemokine concentration that triggers a T cell")
        ->check(CLI::Range(0.0, 1.0))
        ->capture_default_str();
    app.add_option("--antibody-factor", antibody_factor,
                   "Impact of antibodies; multiplier for virion decay")
        ->capture_default_str();
    app.add_option("--antibody-period", antibody_period,
                   "Number of time steps before antibodies start to be produced")
        ->capture_default_str();
    app.add_option("--tcell-generation-rate", tcell_generation_rate,
                   "Number of tcells generated at each timestep for the whole lung")
        ->capture_default_str();
    app.add_option("--tcell-initial-delay", tcell_initial_delay,
                   "Number of time steps before T cells start to be produced")
        ->capture_default_str();
    app.add_option("--tcell-vascular-period", tcell_vascular_period,
                   "Average number of time steps to death for a T cell in the vasculature")
        ->capture_default_str();
    app.add_option("--tcell-tissue-period", tcell_tissue_period,
                   "Average number of time steps to death after a T cell extravasates")
        ->capture_default_str();
    app.add_option("--tcell-binding-period", tcell_binding_period,
                   "Number of time steps a T cell is bound to an epithelial cell when inducing "
                   "apoptosis")
        ->capture_default_str();
    app.add_option("--max-binding-prob", max_binding_prob,
                   "Max probability of a T cell binding to an infected cell in one time step")
        ->capture_default_str();
    app.add_flag("--tcells-follow-gradient", tcells_follow_gradient,
                 "T cells in tissue follow the chemokine gradient")
        ->capture_default_str();
    app.add_option("-r,--seed", rnd_seed, "Random seed")->capture_default_str();
    app.add_option("--sample-period", sample_period,
                   "Number of timesteps between samples (set to 0 to disable sampling)")
        ->capture_default_str();
    app.add_option("--sample-resolution", sample_resolution, "Resolution for sampling")
        ->check(CLI::Range(1, 10000))
        ->capture_default_str();
    app.add_option("--max-block-dim", max_block_dim,
                   "Max. block dimension - larger means more locality but worse load balance. Set "
                   "to 0 for largest possible")
        ->capture_default_str();

    auto *output_dir_opt = app.add_option("-o,--output", output_dir, "Output directory");
    app.add_flag("--progress", show_progress, "Show progress");
    app.add_flag("-v, --verbose", verbose, "Verbose output");

    auto *cfg_opt = app.set_config("--config", "", "Load options from a configuration file");

    try {
      app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
      if (upcxx::rank_me() == 0) app.exit(e);
      return false;
    }

    upcxx::barrier();

    _rnd_gen = make_shared<Random>(rnd_seed + rank_me());

    if (!*output_dir_opt) {
      output_dir = "simcov-run-n" + to_string(upcxx::rank_n()) + "-N" +
                   to_string(upcxx::rank_n() / upcxx::local_team().rank_n()) + "-" +
                   get_current_time(true);
      output_dir_opt->default_val(output_dir);
    }

    if (virion_decay_rate * antibody_factor > 1.0) {
      if (!rank_me())
        cerr << "Invalid parameter settings: virion-decay * antibody_factor > 1.\n"
             << "Reduce either or both of those settings\n";
      return false;
    }

    if (!infection_coords_strs.empty() && !parse_infection_coords(infection_coords_strs))
      return false;

    if (!max_block_dim) {
      max_block_dim = min(dimensions[0], dimensions[1]);
      if (dimensions[2] > 1) max_block_dim = min(dimensions[2], max_block_dim);
    }

    if (dimensions[0] % sample_resolution || dimensions[1] % sample_resolution ||
        (dimensions[2] > 1 && dimensions[2] % sample_resolution)) {
      if (!rank_me())
        cerr << "Error: sample period " << sample_resolution
             << " must be a factor of all the dimensions\n";
      return false;
    }

    setup_output_dir();
    setup_log_file();

    init_logger(output_dir + "/simcov.log", verbose);

#ifdef DEBUG
    open_dbg("debug");
#endif

    SLOG(KLBLUE, "SimCov version ", full_version_str, KNORM, "\n");

    if (upcxx::rank_me() == 0) {
      // print out all compiler definitions
      SLOG_VERBOSE(KLBLUE, "_________________________", KNORM, "\n");
      SLOG_VERBOSE(KLBLUE, "Compiler definitions:", KNORM, "\n");
      std::istringstream all_defs_ss(ALL_DEFNS);
      vector<string> all_defs((std::istream_iterator<string>(all_defs_ss)),
                              std::istream_iterator<string>());
      for (auto &def : all_defs) SLOG_VERBOSE("  ", def, "\n");
      SLOG_VERBOSE("_________________________\n");
      SLOG(KLBLUE, "Options:", KNORM, "\n");
      auto all_opts_strings = app.config_to_str_vector(true, false);
      for (auto &opt_str : all_opts_strings) SLOG(KLBLUE, opt_str, KNORM, "\n");
      SLOG(KLBLUE, "_________________________", KNORM, "\n");
    }
    auto num_nodes = upcxx::rank_n() / upcxx::local_team().rank_n();
    SLOG("Starting run with ", upcxx::rank_n(), " processes on ", num_nodes, " node",
         (num_nodes > 1 ? "s" : ""), " at ", get_current_time(), "\n");
#ifdef DEBUG
    SWARN("Running low-performance debug mode");
#endif
    if (!upcxx::rank_me()) {
      // write out configuration file for restarts
      ofstream ofs(output_dir + "/simcov.config");
      ofs << app.config_to_str(true, true);
    }
    upcxx::barrier();
    return true;
  }
};
