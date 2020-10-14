#pragma once

#include <sys/stat.h>

#include <iostream>
#include <regex>
#include <upcxx/upcxx.hpp>

#include "CLI11.hpp"
#include "version.h"

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
    // all change to the output directory
    if (chdir(output_dir.c_str()) == -1 && !upcxx::rank_me()) {
      ostringstream oss;
      oss << KLRED << "Cannot change to output directory " << output_dir << ": " << strerror(errno)
          << KNORM << endl;
      throw std::runtime_error(oss.str());
    }
    upcxx::barrier();
    // now setup a samples subdirectory
    if (!upcxx::rank_me()) {
      // create the output directory and stripe it
      if (mkdir("samples", S_IRWXU) == -1) {
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
      // check to see if mhmxx.log exists. If so, rename it
      if (file_exists("mhmxx.log")) {
        string new_log_fname = "simcov-" + get_current_time(true) + ".log";
        cerr << KLRED << "WARNING: " << KNORM << output_dir << "/simcov.log exists. Renaming to "
             << output_dir << "/" << new_log_fname << endl;
        if (rename("simcov.log", new_log_fname.c_str()) == -1)
          DIE("Could not rename simcov.log: ", strerror(errno));
      }
    }
    upcxx::barrier();
  }

 public:
  vector<int64_t> dimensions{50, 50, 1};
  // each time step should be about 1 minute, so one day = 1440 time steps
  int num_timesteps = 200;
  int num_infections = 3;

  vector<int64_t> infection_coords{-1, -1, -1};
  int initial_infection = 10000;

  // these periods are normally distributed with mean and stddev
  vector<int> incubation_period{30, 3};
  vector<int> apoptosis_period{30, 3};
  vector<int> infection_period{80, 8};

  double tcell_generation_rate = 1;
  int tcell_initial_delay = 10;
  vector<int> tcell_vascular_period{280, 28};
  vector<int> tcell_tissue_period{6, 1};
  int tcell_binding_period = 1;

  double infectivity = 0.008;
  int virion_production = 2;
  double virion_decay_rate = 0.14;
  double virion_diffusion_coef = 1.0;

  double chemokine_production = 0.0005;
  double chemokine_decay_rate = 0.015;
  double chemokine_diffusion_coef = 0.7;

  double antibody_factor = 1;
  int antibody_period = 5760;

  unsigned rnd_seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  string output_dir = "simcov-run-n" + to_string(upcxx::rank_n()) + "-N" +
                      to_string(upcxx::rank_n() / upcxx::local_team().rank_n()) + "-" +
                      get_current_time(true);
  int sample_period = 1;
  double sample_resolution = 1.0;

  bool show_progress = false;
  bool verbose = false;

  bool load(int argc, char **argv) {
    // SIMCOV version v0.1-a0decc6-master (Release) built on 2020-04-08T22:15:40 with g++
    string full_version_str = "SimCov version " + string(SIMCOV_VERSION) + "-" +
                              string(SIMCOV_BRANCH) + " built on " + string(SIMCOV_BUILD_DATE);
    CLI::App app(full_version_str);
    app.add_option("-d,--dim", dimensions, "Dimensions: x y z")
        ->delimiter(',')
        ->expected(3)
        ->capture_default_str();
    app.add_option("-t,--timesteps", num_timesteps, "Number of timesteps")
        ->check(CLI::Range(1, 1000000))
        ->capture_default_str();
    app.add_option("--infections", num_infections, "Number of starting infections")
        ->capture_default_str();
    app.add_option("--infection-coords", infection_coords, "Location of initial infection")
        ->delimiter(',')
        ->expected(3)
        ->capture_default_str();
    app.add_option("--initial-infection", initial_infection,
                   "Number of virions at initial infection locations")
        ->capture_default_str();
    app.add_option(
           "--incubation-period", incubation_period,
           "Number of time steps to expressing virions after cell is infected (average:stddev)")
        ->delimiter(':')
        ->expected(2)
        ->check(CLI::Range(1, 50000))
        ->capture_default_str();
    app.add_option("--apoptosis-period", apoptosis_period,
                   "Number of time steps to death after apoptosis is induced (average:stddev)")
        ->delimiter(':')
        ->expected(2)
        ->check(CLI::Range(1, 50000))
        ->capture_default_str();
    app.add_option("--infected-period", infection_period,
                   "Number of time steps to death after a cell is infected (average:stddev)")
        ->delimiter(':')
        ->expected(2)
        ->check(CLI::Range(1, 50000))
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
    app.add_option("--antibody-factor", antibody_factor,
                   "Impact of antibodies; multiplier for virion decay")
        ->capture_default_str();
    app.add_option("--antibody-period", antibody_period,
                   "Number of time steps before antibodies start to be produced")
        ->capture_default_str();
    app.add_option("--tcell-generation-rate", tcell_generation_rate,
                   "Number of tcells generated at each timestep")
        ->capture_default_str();
    app.add_option("--tcell-initial-delay", tcell_initial_delay,
                   "Number of time steps before T cells start to be produced")
        ->capture_default_str();
    app.add_option("--tcell-vascular-period", tcell_vascular_period,
                   "Number of time steps to death for a T cell in the vasculature (average:stddev)")
        ->delimiter(':')
        ->expected(2)
        ->check(CLI::Range(1, 50000))
        ->capture_default_str();
    app.add_option("--tcell-tissue-period", tcell_tissue_period,
                   "Number of time steps to death after a T cell extravasates (average:stddev)")
        ->delimiter(':')
        ->expected(2)
        ->check(CLI::Range(1, 50000))
        ->capture_default_str();
    app.add_option("--tcell-binding-period", tcell_binding_period,
                   "Number of time steps a T cell is bound to an epithelial cell when inducing "
                   "apoptosis")
        ->capture_default_str();
    app.add_option("-r,--seed", rnd_seed, "Random seed")->capture_default_str();
    app.add_option("--sample-period", sample_period,
                   "Number of timesteps between samples (set to 0 to disable sampling)")
        ->capture_default_str();
    app.add_option("--sample-resolution", sample_resolution, "Resolution for sampling")
        ->check(CLI::Range(0.0, 1.0))
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
    if (infection_coords[0] != -1 && num_infections > 1) {
      num_infections = 1;
      SLOG("Initial infection coords specified, setting number of infection points to 1\n");
    }
    setup_output_dir();
    setup_log_file();

    init_logger("simcov.log", verbose);

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
      ofstream ofs("simcov.config");
      ofs << app.config_to_str(true, true);
    }
    upcxx::barrier();
    return true;
  }
};

