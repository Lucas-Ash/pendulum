#include "modules/reporting.h"

#include <iomanip>
#include <iostream>
#include <string>

void print_accuracy_report(const ExperimentConfig& config,
                           const SimulationResult& result) {
    std::string reference;
    if (config.error_mode == error_reference::Mode::None) {
        reference = "disabled (numerical trajectory used as its own reference)";
    } else if (config.error_mode == error_reference::Mode::HdReference) {
        reference = "high-definition numerical reference (dt / " +
                    std::to_string(config.error_reference_factor) + ")";
    } else {
        reference = "small-angle analytical solution (linearized ODE)";
        if (config.analytical_model == "jacobi") {
            reference = config.restoring_force.model == restoring_force::Model::Polynomial
                            ? "Jacobi-elliptic undamped Duffing reference"
                            : "Jacobi-elliptic nonlinear pendulum reference";
        } else if (config.analytical_model == "duffing_jacobi" ||
                   config.analytical_model == "duffing") {
            reference = "Jacobi-elliptic undamped Duffing reference";
        }
    }

    std::cout << "\n=== ACCURACY ANALYSIS ===" << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Time step: " << config.dt << " s" << std::endl;
    std::cout << "Reference: " << reference << std::endl;
    std::cout << "Initial conditions: θ₀ = " << config.theta0
              << " rad, ω₀ = " << config.omega0 << " rad/s" << std::endl;
    std::cout << "Pendulum parameters: L = " << config.length
              << " m, g = " << config.gravity << " m/s²" << std::endl;
    std::cout << "Restoring force: "
              << restoring_force::to_string(config.restoring_force.model);
    if (config.restoring_force.model == restoring_force::Model::Polynomial) {
        std::cout << " (linear=" << config.restoring_force.linear
                  << ", cubic=" << config.restoring_force.cubic << ")";
    }
    std::cout << std::endl;
    std::cout << "Error analysis mode: " << error_reference::to_string(config.error_mode);
    if (config.error_mode == error_reference::Mode::HdReference) {
        std::cout << " (factor=" << config.error_reference_factor << ")";
    }
    std::cout << std::endl;

    std::cout << "\nθ error statistics:" << std::endl;
    std::cout << std::scientific << std::setprecision(12);
    std::cout << "  Maximum absolute error: " << result.theta_stats.max_abs << " rad" << std::endl;
    std::cout << "  Average absolute error: " << result.theta_stats.avg_abs << " rad" << std::endl;
    std::cout << "  Maximum relative error: " << result.theta_stats.max_rel << std::endl;
    std::cout << "  Average relative error: " << result.theta_stats.avg_rel << std::endl;

    std::cout << "\nω error statistics:" << std::endl;
    std::cout << "  Maximum absolute error: " << result.omega_stats.max_abs << " rad/s" << std::endl;
    std::cout << "  Average absolute error: " << result.omega_stats.avg_abs << " rad/s" << std::endl;
    std::cout << "  Maximum relative error: " << result.omega_stats.max_rel << std::endl;
    std::cout << "  Average relative error: " << result.omega_stats.avg_rel << std::endl;
    std::cout << std::defaultfloat;
}
