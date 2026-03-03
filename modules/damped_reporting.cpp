#include "modules/damped_reporting.h"

#include <cmath>
#include <iostream>

void print_damped_simulation_summary(const DampedConfig& config,
                                     const DampedSimulationResult& result) {
    const auto& p = config.physical;
    const auto& s = config.simulation;
    const auto& settings = config.settings;

    std::cout << "============================================\n";
    std::cout << "   Damped Pendulum Simulation (Modular)\n";
    std::cout << "============================================\n";
    std::cout << "  g           = " << p.g << " m/s^2\n";
    std::cout << "  L           = " << p.L << " m\n";
    std::cout << "  omega_0     = " << result.omega0 << " rad/s\n";
    std::cout << "  gamma       = " << p.gamma << " s^-1\n";
    std::cout << "  omega_d     = " << result.omega_d << " rad/s\n";
    std::cout << "  zeta        = " << result.zeta << " (underdamped)\n";
    std::cout << "  theta_0     = " << p.theta0 << " rad ("
              << p.theta0 * 180.0 / M_PI << " deg)\n";
    std::cout << "  theta_dot_0 = " << p.theta_dot0 << " rad/s\n";
    std::cout << "  dt          = " << s.dt << " s\n";
    std::cout << "  t_start     = " << s.t_start << " s\n";
    std::cout << "  t_end       = " << s.t_end << " s\n";
    std::cout << "  output_every= " << s.output_every << "\n";
    std::cout << "  RK4 steps   = " << result.rk4_steps << "\n";
    std::cout << "  Data points = " << result.samples.size() << "\n";
    std::cout << "  Plot method = " << to_string(settings.plotting_method) << "\n";
    std::cout << "============================================\n\n";

    std::cout << "Max |error| RK4-linear vs analytical:   "
              << result.max_err_linear << "\n";
    std::cout << "Max |error| RK4-nonlinear vs analytical: "
              << result.max_err_nonlinear << "\n\n";
}
