#pragma once

#include <vector>

struct DrivenSample {
    double t = 0.0;
    double theta_numerical = 0.0;
    double theta_analytical = 0.0;
    double difference = 0.0;
};

struct DrivenSimulationResult {
    std::vector<DrivenSample> samples;
    double max_diff = 0.0;
    double avg_diff = 0.0;
    int rk4_steps = 0;
};
