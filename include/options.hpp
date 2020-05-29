#ifndef __OPTIONS_H
#define __OPTIONS_H

#include <iostream>
#include <regex>
#include <sys/stat.h>
#include <upcxx/upcxx.hpp>

#include "version.h"

#include "CLI11.hpp"

using std::cout;
using std::endl;
using std::vector;

#include "upcxx_utils/log.hpp"
#include "upcxx_utils/timers.hpp"

using namespace upcxx_utils;


#define YES_NO(X) ((X) ? "YES" : "NO")

class Options {
  vector<string> splitter(string in_pattern, string& content) {
    vector<string> split_content;
    std::regex pattern(in_pattern);
    copy(std::sregex_token_iterator(content.begin(), content.end(), pattern, -1),
         std::sregex_token_iterator(), back_inserter(split_content));
    return split_content;
  }

  template <typename T>
  string vec_to_str(const vector<T> &vec, const string &delimiter=",") {
    std::ostringstream oss;
    for (auto elem : vec) {
      oss << elem;
      if (elem != vec.back()) oss << delimiter;
    }
    return oss.str();
  }

  bool extract_following_lens(const string &s, vector<unsigned> &lens, int k) {
    for (int i = 0; i < lens.size(); i++) {
      if (lens[i] == k) {
        lens.erase(lens.begin(), lens.begin() + i + 1);
        return true;
      }
    }
    return false;
  }

  void get_restart_options() {
    // read existing mhmxx.log and get most recent completed stage
    ifstream log_file("mhmxx.log");
    if (!log_file.good()) SDIE("Cannot open previous log file mhmxx.log");
    string last_stage;
    while (log_file) {
      char buf[1000];
      log_file.getline(buf, 999);
      string line(buf);
      if (line.find("Completed contig") != string::npos) last_stage = line;
      else if (line.find("Completed scaffolding") != string::npos) last_stage = line;
      if (line.find("Finished in") != string::npos) {
        SDIE("Found the end of the previous run in mhmxx.log, cannot restart");
        upcxx::barrier();
      }
    }
    if (!last_stage.empty()) {
      auto fields = split(last_stage, ' ');
      string stage_type = fields[3];
      int k = std::stoi(fields[7]);
      if (stage_type == "contig") {
        if (k == kmer_lens.back()) {
          max_kmer_len = kmer_lens.back();
          kmer_lens = {};
          stage_type = "scaffolding";
          k = scaff_kmer_lens[0];
        } else {
          if (!extract_following_lens(last_stage, kmer_lens, k))
            SDIE("Cannot find kmer length ", k, " in configuration: ", vec_to_str(kmer_lens));
          prev_kmer_len = k;
        }
        ctgs_fname = "contigs-" + to_string(k) + ".fasta";
      } else {
        max_kmer_len = kmer_lens.back();
        kmer_lens = {};
        if (!extract_following_lens(last_stage, scaff_kmer_lens, k))
            SDIE("Cannot find kmer length ", k, " in configuration: ", vec_to_str(scaff_kmer_lens));
        if (k == scaff_kmer_lens.front()) ctgs_fname = "contigs-" + to_string(k) + ".fasta";
        else ctgs_fname = "scaff-contigs-" + to_string(k) + ".fasta";
      }
      SLOG("\n*** Restarting from previous run at stage ", stage_type, " k = ",
           (stage_type == "contig" ? kmer_lens[0] : scaff_kmer_lens[0]), " ***\n\n");
      SLOG(KLBLUE, "Restart options:\n",
           "  kmer-lens =              ", vec_to_str(kmer_lens), "\n",
           "  scaff-kmer-lens =        ", vec_to_str(scaff_kmer_lens), "\n",
           "  prev-kmer-len =          ", prev_kmer_len, "\n",
           "  max-kmer-len =           ", max_kmer_len, "\n",
           "  contigs =                ", ctgs_fname, KNORM, "\n");
      if (!upcxx::rank_me() && !file_exists(ctgs_fname))
        SDIE("File ", ctgs_fname, " not found. Did the previous run have --checkpoint enabled?");
    } else {
      SLOG("No previously completed stage found, restarting from the beginning\n");
    }
  }

  void setup_output_dir() {
    if (!upcxx::rank_me()) {
      // create the output directory and stripe it if not doing a restart
      if (restart) {
        if (access(output_dir.c_str(), F_OK) == -1) {
          ostringstream oss;
          oss << KLRED << "WARNING: " << KNORM << "Output directory " << output_dir << " for restart does not exist" << endl;
          throw std::runtime_error(oss.str());
        }
      } else {
        if (mkdir(output_dir.c_str(), S_IRWXU) == -1) {
          // could not create the directory
          if (errno == EEXIST) {
            cerr << KLRED << "WARNING: " << KNORM << "Output directory " << output_dir
                 << " already exists. May overwrite existing files\n";
          } else {
            ostringstream oss;
            oss << KLRED << "ERROR: " << KNORM  << " Could not create output directory " << output_dir
                << ": " << strerror(errno) << endl;
            throw std::runtime_error(oss.str());
          }
        } else {
          // created the directory - now stripe it if possible
          auto status = std::system("which lfs 2>&1 > /dev/null");
          if (WIFEXITED(status) && WEXITSTATUS(status) == 0) {
            string cmd = "lfs setstripe -c -1 " + output_dir;
            auto status = std::system(cmd.c_str());
            if (WIFEXITED(status) && WEXITSTATUS(status) == 0) cout << "Set Lustre striping on the output directory\n";
            else cout << "Failed to set Lustre striping on output directory: " << WEXITSTATUS(status) << endl;
          }
        }
      }
    }
    upcxx::barrier();
    // after we change to the output directory, relative paths will be incorrect, so we need to fix them
    for (auto &fname : reads_fnames) {
      if (fname[0] != '/') fname = "../" + fname;
    }
    // all change to the output directory
    if (chdir(output_dir.c_str()) == -1 && !upcxx::rank_me()) {
      ostringstream oss;
      oss << KLRED << "Cannot change to output directory " << output_dir << ": " << strerror(errno) << KNORM << endl;
      throw std::runtime_error(oss.str());
    }
    upcxx::barrier();
  }

  void setup_log_file() {
    if (!upcxx::rank_me()) {
      // check to see if mhmxx.log exists. If so, and not restarting, rename it
      if (file_exists("mhmxx.log") && !restart) {
        string new_log_fname = "mhmxx-" + get_current_time(true) + ".log";
        cerr << KLRED << "WARNING: " << KNORM << output_dir << "/mhmxx.log exists. Renaming to " << output_dir << "/"
             << new_log_fname << endl;
        if (rename("mhmxx.log", new_log_fname.c_str()) == -1) DIE("Could not rename mhmxx.log: ", strerror(errno));
      } else if (!file_exists("mhmxx.log") && restart) {
        ostringstream oss;
        oss << KLRED << "ERROR: " << KNORM << " Could not restart - missing mhmxx.log in this directory" << endl;
        throw std::runtime_error(oss.str());
      }
    }
    upcxx::barrier();
  }

public:

  vector<string> reads_fnames;
  vector<unsigned> kmer_lens = {21, 33, 55, 77, 99};
  int max_kmer_len = 0;
  int prev_kmer_len = 0;
  vector<unsigned> scaff_kmer_lens = {99, 33};
  int qual_offset = 33;
  bool verbose = false;
  int max_kmer_store_mb = 50;
  int max_rpcs_in_flight = 100;
  bool force_bloom = false;
  double dynamic_min_depth = 0.9;
  int dmin_thres = 2.0;
  bool checkpoint = true;
  bool post_assm_aln = false;
  bool post_assm_only = false;
  bool dump_gfa = false;
  bool show_progress = false;
  string ctgs_fname;
#ifdef USE_KMER_DEPTHS
  string kmer_depths_fname;
#endif
  vector<int> insert_size = {0, 0};
  int min_ctg_print_len = 500;
  int break_scaff_Ns = 10;
  string output_dir = "mhmxx-run-<reads_fname[0]>-n" + to_string(upcxx::rank_n()) + "-N" +
      to_string(upcxx::rank_n() / upcxx::local_team().rank_n()) + "-" + get_current_time(true);
  bool restart = false;

  bool load(int argc, char **argv) {
    //SIMCOV version v0.1-a0decc6-master (Release) built on 2020-04-08T22:15:40 with g++
    string full_version_str = "SimCov version " + string(SIMCOV_VERSION) + "-" + string(SIMCOV_BRANCH) + " built on " +
        string(SIMCOV_BUILD_DATE);
    CLI::App app(full_version_str);
    app.add_option("-r, --reads", reads_fnames,
                   "Files containing merged and unmerged reads in FASTQ format (comma separated)")
                   ->delimiter(',') ->check(CLI::ExistingFile);
    /*
                   ->check([](const string &s) {
                     if (s == "lala") return s;
                     return string();
                   });
     */
    app.add_option("-i, --insert", insert_size,
                   "Insert size (average:stddev)")
                   ->delimiter(':') ->expected(2) ->check(CLI::Range(1, 10000));
    auto *kmer_lens_opt = app.add_option("-k, --kmer-lens", kmer_lens,
                   "kmer lengths (comma separated)")
                   ->delimiter(',') ->capture_default_str();
    app.add_option("--max-kmer-len", max_kmer_len,
                   "Maximum kmer length (need to specify if only scaffolding)")
                   ->capture_default_str() ->check(CLI::Range(0, 159));
    app.add_option("--prev-kmer-len", prev_kmer_len,
                   "Previous kmer length (need to specify if contigging and contig file is specified)")
                   ->capture_default_str() ->check(CLI::Range(0, 159));
    auto *scaff_kmer_lens_opt = app.add_option("-s, --scaff-kmer-lens", scaff_kmer_lens,
                   "kmer lengths for scaffolding (comma separated)")
                   ->delimiter(',') ->capture_default_str();
    app.add_option("-Q, --quality-offset", qual_offset,
                   "Phred encoding offset")
                   ->capture_default_str() ->check(CLI::IsMember({0, 33, 64}));
    app.add_option("-c, --contigs", ctgs_fname,
                   "File with contigs used for restart");
//                   ->check(CLI::ExistingFile);
    app.add_option("--dynamic-min-depth", dynamic_min_depth,
                   "Dynamic min. depth for DeBruijn graph traversal - set to 1.0 for a single genome")
                   ->capture_default_str() ->check(CLI::Range(0.1, 1.0));
    app.add_option("--min-depth-thres", dmin_thres,
                   "Absolute mininimum depth threshold for DeBruijn graph traversal")
                   ->capture_default_str() ->check(CLI::Range(1, 100));
    app.add_option("--max-kmer-store", max_kmer_store_mb,
                   "Maximum size for kmer store in MB")
                   ->capture_default_str() ->check(CLI::Range(1, 1000));
    app.add_option("--max-rpcs-in-flight", max_rpcs_in_flight,
                   "Maximum number of RPCs in flight, per process (0 = unlimited)")
                   ->capture_default_str() ->check(CLI::Range(0, 10000));
    app.add_option("--min-ctg-print-len", min_ctg_print_len,
                   "Minimum length required for printing a contig in the final assembly")
                   ->capture_default_str() ->check(CLI::Range(0, 100000));
    app.add_option("--break-scaff-Ns", break_scaff_Ns,
                   "Number of Ns allowed before a scaffold is broken")
                   ->capture_default_str() ->check(CLI::Range(0, 1000));
    auto *output_dir_opt = app.add_option("-o,--output", output_dir, "Output directory")
                                          ->capture_default_str();
    app.add_flag("--force-bloom", force_bloom,
                 "Always use bloom filters")
                 ->default_val(force_bloom ? "true" : "false") ->capture_default_str() ->multi_option_policy();
    app.add_flag("--checkpoint", checkpoint,
                 "Checkpoint after each contig round")
                 ->default_val(checkpoint ? "true" : "false") ->capture_default_str() ->multi_option_policy();
    app.add_flag("--restart", restart,
                 "Restart in previous directory where a run failed")
                 ->capture_default_str();
    app.add_flag("--post-assembly-align", post_assm_aln,
                 "Align reads to final assembly")
                 ->capture_default_str();
    app.add_flag("--write-gfa", dump_gfa,
                 "Dump scaffolding contig graphs in GFA2 format")
                 ->capture_default_str();
    app.add_flag("--post-assembly-only", post_assm_only,
                 "Only run post assembly")
                 ->capture_default_str();
    app.add_flag("--progress", show_progress,
                 "Show progress")
                 ->capture_default_str();
    app.add_flag("-v, --verbose", verbose,
                 "Verbose output")
                 ->capture_default_str();

    auto *cfg_opt = app.set_config("--config", "simcov.config", "Load options from a configuration file");

    try {
      app.parse(argc, argv);
    } catch(const CLI::ParseError &e) {
      if (upcxx::rank_me() == 0) app.exit(e);
      return false;
    }

    // we can get restart incorrectly set in the config file from restarted runs
    if (*cfg_opt) {
      restart = false;
      app.get_option("--restart")->default_val("false");
    }

    if (!upcxx::rank_me() && !restart && reads_fnames.empty()) {
      // FIXME: there appears to be no way to do exactly this with a validator or ->required()
      ostringstream oss;
      oss << KLRED << "Require read names if not restarting" << KNORM << endl;
      throw std::runtime_error(oss.str());
    }

    if (!upcxx::rank_me() && post_assm_only && (!max_kmer_len || ctgs_fname.empty())) {
      ostringstream oss;
      oss << KLRED << "For running only post assembly analysis, require --max-kmer_len and --contigs" << KNORM << endl;
      throw std::runtime_error(oss.str());
    }

    upcxx::barrier();

    if (!*output_dir_opt) {
      string first_read_fname = remove_file_ext(get_basename(reads_fnames[0]));
      output_dir = "simcov-run-" + first_read_fname + "-n" + to_string(upcxx::rank_n()) + "-N" +
          to_string(upcxx::rank_n() / upcxx::local_team().rank_n()) + "-" + get_current_time(true);
      output_dir_opt->default_val(output_dir);
    }

    if (restart) {
      // this mucking about is to ensure we don't get multiple failures messages if the config file does not parse
      if (!upcxx::rank_me()) app.parse_config(output_dir + "/simcov.config");
      upcxx::barrier();
      if (upcxx::rank_me()) app.parse_config(output_dir + "/simcov.config");
    }

    // make sure we only use defaults for kmer lens if none of them were set by the user
    if (*kmer_lens_opt && !*scaff_kmer_lens_opt) scaff_kmer_lens = {};
    if (*scaff_kmer_lens_opt && !*kmer_lens_opt) kmer_lens = {};

    setup_output_dir();
    setup_log_file();

    init_logger("simcov.log", verbose);

#ifdef DEBUG
    open_dbg("debug");
#endif

    SLOG(KLBLUE, "SimCov version ", full_version_str, KNORM, "\n");

    if (restart) get_restart_options();

    if (upcxx::rank_me() == 0) {
      // print out all compiler definitions
      SLOG_VERBOSE(KLBLUE, "_________________________", KNORM, "\n");
      SLOG_VERBOSE(KLBLUE, "Compiler definitions:", KNORM, "\n");
      std::istringstream all_defs_ss(ALL_DEFNS);
      vector<string> all_defs((std::istream_iterator<string>(all_defs_ss)), std::istream_iterator<string>());
      for (auto &def : all_defs) SLOG_VERBOSE("  ", def, "\n");
      SLOG_VERBOSE("_________________________\n");
      SLOG(KLBLUE, "Options:", KNORM, "\n");
      auto all_opts_strings = app.config_to_str_vector(true, false);
      for (auto &opt_str : all_opts_strings) SLOG(KLBLUE, opt_str, KNORM, "\n");
      SLOG(KLBLUE, "_________________________", KNORM, "\n");
    }
    auto num_nodes = upcxx::rank_n() / upcxx::local_team().rank_n();
    SLOG("Starting run with ", upcxx::rank_n(), " processes on ", num_nodes, " node", (num_nodes > 1 ? "s" : ""), " at ",
         get_current_time(), "\n");
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


#endif
