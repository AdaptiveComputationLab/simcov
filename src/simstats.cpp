#pragma once
#include "simstats.hpp"

/**
SimStats function definitions
**/

SimStats::SimStats(std::shared_ptr<Options> opt) {
    options = opt;
}

void SimStats::init() {
    if (!upcxx::rank_me()) {
        log_file.open(options->output_dir + "/simcov.stats");
        log_file << "# time" << header(0) << endl;
    }
}

std::string SimStats::header(int width) {
    std::vector<string> columns = {"incb", "expr", "apop", "dead",    "tvas",
                              "ttis", "chem", "virs", "chempts", "%infct"};
    std::ostringstream oss;
    oss << std::left;
    for (auto column : columns) {
        if (width)
            oss << setw(width) << column;
        else
            oss << '\t' << column;
    }
    return oss.str();
}

std::string SimStats::to_str(int width) {
    std::vector<int64_t> totals;
    totals.push_back(upcxx::reduce_one(incubating, upcxx::op_fast_add, 0).wait());
    totals.push_back(upcxx::reduce_one(expressing, upcxx::op_fast_add, 0).wait());
    totals.push_back(upcxx::reduce_one(apoptotic, upcxx::op_fast_add, 0).wait());
    totals.push_back(upcxx::reduce_one(dead, upcxx::op_fast_add, 0).wait());
    totals.push_back(upcxx::reduce_one(tcells_vasculature, upcxx::op_fast_add, 0).wait());
    totals.push_back(upcxx::reduce_one(tcells_tissue, upcxx::op_fast_add, 0).wait());
    std::vector<float> totals_d;
    totals_d.push_back(upcxx::reduce_one(chemokines, upcxx::op_fast_add, 0).wait() / get_num_grid_points());
    totals_d.push_back(upcxx::reduce_one(virions, upcxx::op_fast_add, 0).wait());  // / get_num_grid_points());
    auto all_chem_pts = upcxx::reduce_one(num_chemo_pts, upcxx::op_fast_add, 0).wait();
    totals_d.push_back(all_chem_pts + totals[0] + totals[1] + totals[2] + totals[3]);
    auto perc_infected =
        100.0 * (float)(totals[0] + totals[1] + totals[2] + totals[3]) / get_num_grid_points();

    std::ostringstream oss;
    oss << std::left;
    for (auto tot : totals) {
        if (width)
            oss << std::setw(width) << tot;
        else
            oss << '\t' << tot;
    }
    oss << std::fixed << std::setprecision(2) << std::scientific;
    for (auto tot : totals_d) {
        if (width)
            oss << std::setw(width) << tot;
        else
            oss << '\t' << tot;
    }
    oss << std::fixed << std::setprecision(2) << std::showpoint;
    if (width)
        oss << std::setw(width) << perc_infected;
    else
        oss << '\t' << perc_infected;
    return oss.str();
}

void SimStats::log(int time_step) {
    std::string s = to_str(0);
    if (!upcxx::rank_me()) log_file << time_step << s << std::endl;
}