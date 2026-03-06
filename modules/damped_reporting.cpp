#include "modules/damped_reporting.h"

#include <cmath>
#include <iostream>

void print_damped_simulation_summary(const DampedConfig& config,
                                     const SimulationResult& result) {
    const auto& p = config.physical;
    const auto& s = config.simulation;
    const auto& settings = config.settings;
    
    const double omega0_sq = (p.g / p.L) *
                             restoring_force::linearized_slope(p.restoring_force);
    const double omega0 = omega0_sq > 0.0 ? std::sqrt(omega0_sq) : 0.0;
    const double omega_d_sq = omega0 * omega0 - p.gamma * p.gamma;
    const bool underdamped = omega_d_sq > 0.0;
    const double omega_d = underdamped ? std::sqrt(omega_d_sq) : 0.0;
    const double zeta = omega0 > 0.0 ? p.gamma / omega0 : 0.0;

    std::cout << "============================================\n";
    std::cout << "   Damped Pendulum Simulation (Modular)\n";
    std::cout << "============================================\n";
    std::cout << "  g           = " << p.g << " m/s^2\n";
    std::cout << "  L           = " << p.L << " m\n";
    std::cout << "  omega_0     = " << omega0 << " rad/s\n";
    std::cout << "  gamma       = " << p.gamma << " s^-1\n";
    std::cout << "  restoring   = "
              << restoring_force::to_string(p.restoring_force.model);
    if (p.restoring_force.model == restoring_force::Model::Polynomial) {
        std::cout << " (linear=" << p.restoring_force.linear
                  << ", cubic=" << p.restoring_force.cubic << ")";
    }
    std::cout << "\n";
    if (underdamped) {
        std::cout << "  omega_d     = " << omega_d << " rad/s\n";
        std::cout << "  zeta        = " << zeta << " (underdamped)\n";
    } else {
        std::cout << "  zeta        = " << zeta << " (critical/overdamped)\n";
    }
    std::cout << "  theta_0     = " << p.theta0 << " rad ("
              << p.theta0 * 180.0 / M_PI << " deg)\n";
    std::cout << "  theta_dot_0 = " << p.theta_dot0 << " rad/s\n";
    std::cout << "  dt          = " << s.dt << " s\n";
    std::cout << "  t_start     = " << s.t_start << " s\n";
    std::cout << "  t_end       = " << s.t_end << " s\n";
    std::cout << "  output_every= " << s.output_every << "\n";
    std::cout << "  RK4 steps   = " << result.rk4_steps << "\n";
    std::cout << "  Data points = " << result.t.size() << "\n";
    std::cout << "  Plot method = " << to_string(settings.plotting_method) << "\n";
    std::cout << "  Error mode  = "
              << error_reference::to_string(settings.error_mode);
    if (settings.error_mode == error_reference::Mode::HdReference) {
        std::cout << " (factor=" << settings.error_reference_factor << ")";
    }
    std::cout << "\n";
    std::cout << "============================================\n\n";

    std::cout << "Max |error| theta vs selected reference: "
              << result.theta_stats.max_abs << "\n";
}
