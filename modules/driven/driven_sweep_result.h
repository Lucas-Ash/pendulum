#pragma once

#include <vector>

struct DrivenSweepSample {
    double drive_frequency = 0.0;
    double numerical_amplitude = 0.0;
    double analytical_amplitude = 0.0;
    double analytical_lower_stable_amplitude = 0.0;
    double analytical_upper_stable_amplitude = 0.0;
    double final_theta = 0.0;
    double final_omega = 0.0;
};

struct DrivenSweepResult {
    std::vector<DrivenSweepSample> samples;
};
