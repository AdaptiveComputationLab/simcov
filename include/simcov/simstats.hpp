#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <upcxx/upcxx.hpp>
#include "tissue.hpp"

class SimStats {
    public:
    int64_t incubating = 0;
    int64_t expressing = 0;
    int64_t apoptotic = 0;
    int64_t dead = 0;
    int64_t tcells_vasculature = 0;
    int64_t tcells_tissue = 0;
    float chemokines = 0;
    int64_t num_chemo_pts = 0;
    float virions = 0;

    SimStats(std::shared_ptr<Options> opt);
    void init();
    std::string header(int width);
    std::string to_str(int width);
    void log(int time_step);


    private:
    std::shared_ptr<Options> options;
    std::ofstream log_file;

};
