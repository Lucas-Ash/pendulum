#pragma once

#include <vector>

struct DampedSample {
    double t = 0.0;
    double analytical = 0.0;
    double numerical_linear = 0.0;
    double numerical_nonlinear = 0.0;
    double error_linear = 0.0;
    double error_nonlinear = 0.0;
    double energy_nonlinear = 0.0;
};

struct DampedSimulationResult {
    std::vector<DampedSample> samples;
    double omega0 = 0.0;
    double omega_d = 0.0;
    double zeta = 0.0;
    int rk4_steps = 0;
    double max_err_linear = 0.0;
    double max_err_nonlinear = 0.0;
};
