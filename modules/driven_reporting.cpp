#include "modules/driven_reporting.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>

void print_driven_simulation_summary(const DrivenConfig& config, const SimulationResult& result) {
    std::cout << "\nSimulation complete!\n";
    std::cout << "Total steps: " << result.rk4_steps << "\n";
    std::cout << "Output points: " << result.t.size() << "\n\n";
    
    std::cout << "Comparison Statistics:\n";
    std::cout << "  Maximum difference: " << std::scientific << std::setprecision(6) << result.theta_stats.max_abs << " rad\n";
    std::cout << "  Average difference: " << std::fixed << std::setprecision(6) << result.theta_stats.avg_abs << " rad\n\n";

    std::cout << "System model: " << to_string(config.physical.system_model) << "\n";
    std::cout << "Restoring force: "
              << restoring_force::to_string(config.physical.restoring_force.model);
    if (config.physical.system_model == DrivenSystemModel::Pendulum &&
        config.physical.restoring_force.model == restoring_force::Model::Polynomial) {
        std::cout << " (linear=" << config.physical.restoring_force.linear
                  << ", cubic=" << config.physical.restoring_force.cubic << ")";
    } else if (config.physical.system_model != DrivenSystemModel::Pendulum) {
        std::cout << " (linear stiffness=" << config.physical.linear_stiffness
                  << ", cubic stiffness=" << config.physical.cubic_stiffness
                  << ", mass=" << config.physical.mass << ")";
    }
    std::cout << "\n";
    std::cout << "Damping model: "
              << damping_force::to_string(config.physical.damping_model);
    if (config.physical.damping_model == damping_force::Model::Polynomial) {
        std::cout << " (linear=" << config.physical.damping
                  << ", cubic=" << config.physical.damping_cubic << ")";
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

    if (config.mass_event.enabled) {
        std::cout << "Mass event: enabled at t=" << config.mass_event.jump_time
                  << " with delta_mass=" << config.mass_event.delta_mass << "\n";
    }
    if (config.noise.enabled) {
        std::cout << "Core noise: force_stddev=" << config.noise.force_stddev
                  << ", correlation_time=" << config.noise.correlation_time
                  << ", seed=" << config.noise.seed << "\n";
    }
}

void print_driven_sweep_summary(const DrivenConfig& config, const DrivenSweepResult& result) {
    std::cout << "\nSweep complete!\n";
    std::cout << "System model: " << to_string(config.physical.system_model) << "\n";
    std::cout << "Direction: " << to_string(config.sweep.direction) << "\n";
    std::cout << "Sweep points: " << result.samples.size() << "\n";
    if (result.samples.empty()) {
        return;
    }

    double max_amplitude = 0.0;
    for (const auto& sample : result.samples) {
        max_amplitude = std::max(max_amplitude, std::abs(sample.numerical_amplitude));
    }
    std::cout << "Max numerical amplitude: " << std::scientific << std::setprecision(6)
              << max_amplitude << "\n";
    if (config.sweep.analytical_branch_tracking) {
        std::cout << "Analytical stable-branch tracking: enabled\n";
    }
}
