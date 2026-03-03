#pragma once

#include <vector>

struct ErrorStatistics {
    double max_abs = 0.0;
    double avg_abs = 0.0;
    double max_rel = 0.0;
    double avg_rel = 0.0;
};

struct SimulationResult {
    std::vector<double> t;
    std::vector<double> theta;
    std::vector<double> omega;
    std::vector<double> theta_errors;
    std::vector<double> omega_errors;
    ErrorStatistics theta_stats;
    ErrorStatistics omega_stats;
};
