#pragma once
#include <fcntl.h>
#include <unistd.h>

#include <algorithm>
#include <random>

#include "upcxx_utils/log.hpp"
using namespace upcxx_utils;

using std::min;
using std::string;
using std::string_view;
using std::to_string;

#ifdef USE_BYTELL
#include "bytell_hash_map.hpp"
#define HASH_TABLE ska::bytell_hash_map
#else
#include <unordered_map>
#define HASH_TABLE std::unordered_map
#endif

int pin_thread(pid_t pid, int cid);

void dump_single_file(const string &fname, const string &out_str);

class Random {
 private:
  std::mt19937_64 generator;

  double get_prob(double max_val = 1.0) {
    return std::uniform_real_distribution<>(0, max_val)(generator);
  }

 public:
  Random(unsigned seed) : generator(seed) {}

  int get(int64_t begin, int64_t end) {
    return std::uniform_int_distribution<int64_t>(begin, end - 1)(generator);
  }

  bool trial_success(double thres) {
    assert(thres >= 0);
    if (thres > 1) return true;
    if (thres == 0) return false;
    return (get_prob() <= thres);
  }

  int get_normal(vector<int> dist_params) {
    return (int)std::normal_distribution<float>(dist_params[0], dist_params[1])(generator);
  }

  int get_poisson(int avg) { return (int)std::poisson_distribution<int>(avg)(generator); }
};

extern std::shared_ptr<Random> _rnd_gen;
