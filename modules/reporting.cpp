#include "modules/reporting.h"

#include <iomanip>
#include <iostream>

void print_accuracy_report(const ExperimentConfig& config,
                           const SimulationResult& result) {
    std::cout << "\n=== ACCURACY ANALYSIS ===" << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Time step: " << config.dt << " s" << std::endl;
    std::cout << "Reference: small-angle analytical solution (linearized ODE)" << std::endl;
    std::cout << "Initial conditions: θ₀ = " << config.theta0
              << " rad, ω₀ = " << config.omega0 << " rad/s" << std::endl;
    std::cout << "Pendulum parameters: L = " << config.length
              << " m, g = " << config.gravity << " m/s²" << std::endl;

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
