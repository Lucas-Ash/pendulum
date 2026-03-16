#pragma once

#include "modules/damped_config.h"
#include "modules/simulation_result.h"

class DampedPendulumSimulator {
public:
    explicit DampedPendulumSimulator(const DampedConfig& config);
    SimulationResult simulate() const;

private:
    DampedConfig config_;

    static void linear_analytical_solution(double t, double omega_d, double gamma,
                                           double theta0, double theta_dot0,
                                           double& theta_exact, double& omega_exact);
    static double total_energy(double theta, double theta_dot, double omega0_sq,
                               const restoring_force::Config& restoring);
};
