#include <exception>
#include <iostream>
#include <fstream>
#include <string>

// Simple Pendulum
#include "pendulum/experiment_config.h"
#include "pendulum/pendulum_simulator.h"
#include "pendulum/plotting.h"
#include "pendulum/reporting.h"
#include "pendulum/runtime_env.h"

// Damped Pendulum
#include "pendulum/damped_config.h"
#include "pendulum/damped_io.h"
#include "pendulum/damped_plotting.h"
#include "pendulum/damped_reporting.h"
#include "pendulum/damped_simulator.h"

// Driven Pendulum
#include "pendulum/driven_config.h"
#include "pendulum/driven_io.h"
#include "pendulum/driven_plotting.h"
#include "pendulum/driven_reporting.h"
#include "pendulum/driven_simulator.h"

enum class SimulationType {
    Unknown,
    Simple,
    Damped,
    Driven
};

SimulationType detect_simulation_type(const std::string& config_path) {
    std::ifstream file(config_path);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open config file: " + config_path);
    }

    std::string line;
    bool has_gamma = false;
    bool has_omega_drive = false;

    while (std::getline(file, line)) {
        if (line.find("omega_drive:") != std::string::npos) {
            has_omega_drive = true;
        } else if (line.find("gamma:") != std::string::npos || line.find("damping:") != std::string::npos) {
            has_gamma = true;
        }
    }

    if (has_omega_drive) return SimulationType::Driven;
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

    PendulumSimulator pendulum(config.length, config.gravity, config.dt, config.t_max);
    SimulationResult result = pendulum.simulate(config.theta0, config.omega0);

    print_accuracy_report(config, result);
    plot_simulation_results(result);
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
        DampedSimulationResult result = simulator.simulate();

        write_damped_data_file(config.settings.data_file, result);
        std::cout << "Data written to " << config.settings.data_file << "\n";

        print_damped_simulation_summary(config, result);
        render_damped_plots(config, result);
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
        DrivenSimulationResult result = simulator.simulate();

        write_driven_data_file(config.settings.data_file, result);
        std::cout << "Data saved to " << config.settings.data_file << "\n";

        print_driven_simulation_summary(config, result);
        
        std::cout << "To view the graph, run: python3 " << config.settings.python_script << "\n";
        render_driven_plots(config, result);
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
            default:
                std::cerr << "Error: Could not determine simulation type from config file.\n";
                return 1;
        }
    } catch (const std::exception& ex) {
        std::cerr << "Error detecting simulation type: " << ex.what() << "\n";
        return 1;
    }
}
