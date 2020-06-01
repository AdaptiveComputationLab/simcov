#pragma once

#include <iostream>
#include <map>
#include <fstream>
#include <stdarg.h>
#include <numeric>
#include <algorithm>
#include <upcxx/upcxx.hpp>

#include "upcxx_utils/log.hpp"
#include "upcxx_utils/progress_bar.hpp"
#include "upcxx_utils/flat_aggr_store.hpp"
#include "upcxx_utils/timers.hpp"

#include "utils.hpp"

struct Agent {
  int64_t id;
  int64_t vertex_id;
};

struct Vertex {
  int64_t id;
  vector<int64_t> neighbors;
  vector<Agent> agents;
  int64_t num_vertices;
  int64_t avg_degree;
};

class DHTGraph {
 private:
  using vertex_map_t = upcxx::dist_object<HASH_TABLE<int64_t, Vertex>>;
  vertex_map_t vertices;
  
 public:
  DHTGraph() : vertices({}) {
  }

  void clear() {
    for (auto it = vertices->begin(); it != vertices->end(); ) {
      it = vertices->erase(it);
    }
  }

  ~DHTGraph() {
    clear();
  }

  void construct(int64_t num_vertices, int avg_degree) {
  }

};
