#pragma once

#include "modules/damped_config.h"
#include "modules/damped_result.h"

class DampedPendulumSimulator {
public:
    explicit DampedPendulumSimulator(const DampedConfig& config);
    DampedSimulationResult simulate() const;

private:
    DampedConfig config_;

    static double analytical_solution(double t, double omega_d, double gamma,
                                      double theta0, double theta_dot0);
    static double total_energy(double theta, double theta_dot, double omega0);
};
