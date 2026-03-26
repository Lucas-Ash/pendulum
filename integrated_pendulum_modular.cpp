#include <exception>
#include <iostream>
#include <fstream>
#include <string>

// Simple Pendulum
#include "modules/config/experiment_config.h"
#include "modules/simple/pendulum_simulator.h"
#include "modules/simple/plotting.h"
#include "modules/simple/reporting.h"
#include "modules/simple/simple_io.h"
#include "modules/runtime/runtime_env.h"

// Damped Pendulum
#include "modules/config/damped_config.h"
#include "modules/damped/damped_io.h"
#include "modules/damped/damped_plotting.h"
#include "modules/damped/damped_reporting.h"
#include "modules/damped/damped_simulator.h"

// Driven Pendulum
#include "modules/config/driven_config.h"
#include "modules/driven/driven_io.h"
#include "modules/driven/driven_plotting.h"
#include "modules/driven/driven_reporting.h"
#include "modules/driven/driven_simulator.h"

// Coupled Pendulum
#include "modules/config/coupled_config.h"
#include "modules/coupled/coupled_simulator.h"

enum class SimulationType {
    Unknown,
    Simple,
    Damped,
    Driven,
    Coupled
};

SimulationType detect_simulation_type(const std::string& config_path) {
    std::ifstream file(config_path);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open config file: " + config_path);
    }

    std::string line;
    bool has_gamma = false;
    bool has_omega_drive = false;
    bool has_explicit_driven_mode = false;
    bool has_coupled_mode = false;

    while (std::getline(file, line)) {
        if (line.find("omega_1:") != std::string::npos ||
            line.find("alpha_11:") != std::string::npos ||
            line.find("q1_0:") != std::string::npos) {
            has_coupled_mode = true;
        } else if (line.find("omega_drive:") != std::string::npos) {
            has_omega_drive = true;
        } else if ((line.find("system_model:") != std::string::npos ||
                    line.find("mode:") != std::string::npos) &&
                   line.find("duffing") != std::string::npos) {
            has_explicit_driven_mode = true;
        } else if (line.find("drive_force:") != std::string::npos ||
                   line.find("sweep:") != std::string::npos) {
            has_explicit_driven_mode = true;
        } else if (line.find("gamma:") != std::string::npos ||
                   line.find("damping:") != std::string::npos ||
                   line.find("time_damping_") != std::string::npos ||
                   line.find("lane_emden") != std::string::npos ||
                   line.find("damping_model:") != std::string::npos ||
                   line.find("damping_linear:") != std::string::npos ||
                   line.find("damping_cubic:") != std::string::npos) {
            has_gamma = true;
        }
    }

    if (has_coupled_mode) return SimulationType::Coupled;
    if (has_omega_drive || has_explicit_driven_mode) return SimulationType::Driven;
    if (has_gamma) return SimulationType::Damped;
    return SimulationType::Simple;
}

int run_simple(const std::string& config_path) {
    std::cout << "Detected Simple Pendulum simulation.\n";
    configure_pythonpath_from_venv();

    ExperimentConfig config;
    try {
        config = load_config_from_yaml(config_path);
    } catch (const std::exception& ex) {
        std::cerr << "Failed to load YAML config: " << ex.what() << std::endl;
        return 1;
    }

    std::cout << "Loaded config from: " << config_path << std::endl;
    std::cout << "length=" << config.length << ", gravity=" << config.gravity
              << ", dt=" << config.dt << ", t_max=" << config.t_max
              << ", theta0=" << config.theta0 << ", omega0=" << config.omega0
              << std::endl;

    PendulumSimulator pendulum(
        config.length, config.gravity, config.dt, config.t_max,
        config.restoring_force, config.additional_terms);
    SimulationResult result = pendulum.simulate(
        config.theta0, config.omega0, config.integrator, config.analytical_model,
        config.error_mode, config.error_reference_factor);

    write_simple_data_file(config.data_file, result);
    std::cout << "Data saved to " << config.data_file << "\n";

    print_accuracy_report(config, result);
    plot_simulation_results(config, result);
    // std::cout << "Simple simulation finished successfully. (IO/Plots temporarily disabled)\n";
    return 0;
}

int run_damped(const std::string& config_path) {
    std::cout << "Detected Damped Pendulum simulation.\n\n";

    DampedConfig config;
    try {
        config = load_damped_config_from_yaml(config_path);
    } catch (const std::exception& ex) {
        std::cerr << "Failed to load YAML config: " << ex.what() << "\n";
        return 1;
    }

    std::cout << "Loaded config from: " << config_path << "\n\n";

    try {
        DampedPendulumSimulator simulator(config);
        SimulationResult result = simulator.simulate();

        write_damped_data_file(config.settings.data_file, result);
        std::cout << "Data written to " << config.settings.data_file << "\n";

        print_damped_simulation_summary(config, result);
        render_damped_plots(config, result);
        // std::cout << "Damped simulation finished successfully. (IO/Plots temporarily disabled)\n";
    } catch (const std::exception& ex) {
        std::cerr << "Simulation failed: " << ex.what() << "\n";
        return 1;
    }

    return 0;
}

int run_driven(const std::string& config_path) {
    std::cout << "Detected Damped Driven Pendulum simulation.\n\n";

    DrivenConfig config;
    try {
        config = load_driven_config_from_yaml(config_path);
    } catch (const std::exception& ex) {
        std::cerr << "Failed to load YAML config: " << ex.what() << "\n";
        return 1;
    }

    std::cout << "Loaded config from: " << config_path << "\n\n";

    try {
        std::cout << "Running simulation...\n";
        DrivenPendulumSimulator simulator(config);
        if (config.sweep.enabled) {
            DrivenSweepResult result = simulator.simulate_sweep();
            write_driven_sweep_data_file(config.settings.sweep_data_file, result);
            std::cout << "Sweep data saved to " << config.settings.sweep_data_file << "\n";
            print_driven_sweep_summary(config, result);
            std::cout << "To view the graph, run: python3 " << config.settings.python_script << "\n";
            render_driven_sweep_plots(config, result);
        } else {
            SimulationResult result = simulator.simulate();
            write_driven_data_file(config.settings.data_file, result);
            std::cout << "Data saved to " << config.settings.data_file << "\n";

            print_driven_simulation_summary(config, result);
            std::cout << "To view the graph, run: python3 " << config.settings.python_script << "\n";
            render_driven_plots(config, result);
        }
        // std::cout << "Driven simulation finished successfully. (IO/Plots temporarily disabled)\n";
    } catch (const std::exception& ex) {
        std::cerr << "Simulation failed: " << ex.what() << "\n";
        return 1;
    }

    return 0;
}

int run_coupled(const std::string& config_path) {
    std::cout << "Detected Coupled Duffing simulation.\n\n";

    CoupledConfig config;
    try {
        config = load_coupled_config_from_yaml(config_path);
    } catch (const std::exception& ex) {
        std::cerr << "Failed to load YAML config: " << ex.what() << "\n";
        return 1;
    }

    std::cout << "Loaded config from: " << config_path << "\n\n";

    try {
        std::cout << "Running simulation...\n";
        CoupledSimulator simulator(config);
        CoupledSimulationResult result = simulator.simulate();
        write_coupled_data_file(config.settings.data_file, result);
        std::cout << "Data saved to " << config.settings.data_file << "\n";

        write_coupled_plot_script(config.settings.python_script, config.settings.data_file, config.settings.output_png);
        std::cout << "To view the graph, run: python3 " << config.settings.python_script << "\n";
        
        std::string cmd = "python3 " + config.settings.python_script;
        if (std::system(cmd.c_str()) != 0) {
            std::cerr << "Warning: Plotting script returned non-zero exit code.\n";
        }
    } catch (const std::exception& ex) {
        std::cerr << "Simulation failed: " << ex.what() << "\n";
        return 1;
    }

    return 0;
}

int main(int argc, char* argv[]) {
    std::cout << "==============================================\n";
    std::cout << "     Integrated Pendulum Simulation Runner\n";
    std::cout << "==============================================\n\n";

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <config_file.yaml>\n";
        return 1;
    }

    const std::string config_path = argv[1];

    try {
        SimulationType type = detect_simulation_type(config_path);
        
        switch (type) {
            case SimulationType::Simple:
                return run_simple(config_path);
            case SimulationType::Damped:
                return run_damped(config_path);
            case SimulationType::Driven:
                return run_driven(config_path);
            case SimulationType::Coupled:
                return run_coupled(config_path);
            default:
                std::cerr << "Error: Could not determine simulation type from config file.\n";
                return 1;
        }
    } catch (const std::exception& ex) {
        std::cerr << "Error detecting simulation type: " << ex.what() << "\n";
        return 1;
    }
}
