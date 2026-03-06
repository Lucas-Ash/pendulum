#include "modules/driven_reporting.h"
#include <iostream>
#include <iomanip>

void print_driven_simulation_summary(const DrivenConfig& config, const SimulationResult& result) {
    std::cout << "\nSimulation complete!\n";
    std::cout << "Total steps: " << result.rk4_steps << "\n";
    std::cout << "Output points: " << result.t.size() << "\n\n";
    
    std::cout << "Comparison Statistics:\n";
    std::cout << "  Maximum difference: " << std::scientific << std::setprecision(6) << result.theta_stats.max_abs << " rad\n";
    std::cout << "  Average difference: " << std::fixed << std::setprecision(6) << result.theta_stats.avg_abs << " rad\n\n";

    std::cout << "Restoring force: "
              << restoring_force::to_string(config.physical.restoring_force.model);
    if (config.physical.restoring_force.model == restoring_force::Model::Polynomial) {
        std::cout << " (linear=" << config.physical.restoring_force.linear
                  << ", cubic=" << config.physical.restoring_force.cubic << ")";
    }
    std::cout << "\n";
    std::cout << "Error analysis mode: "
              << error_reference::to_string(config.settings.error_mode);
    if (config.settings.error_mode == error_reference::Mode::HdReference) {
        std::cout << " (factor=" << config.settings.error_reference_factor << ")";
    }
    std::cout << "\n";
    
    if (config.settings.error_mode == error_reference::Mode::Analytical) {
        std::cout << "Note: The analytical solution is linearized about theta=0.\n";
        std::cout << "For large amplitudes or strong nonlinear terms, numerical and analytical traces diverge.\n";
    } else if (config.settings.error_mode == error_reference::Mode::HdReference) {
        std::cout << "Note: Errors are measured against a high-definition numerical reference.\n";
    } else {
        std::cout << "Note: Error analysis is disabled.\n";
    }
}
