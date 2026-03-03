#include <exception>
#include <iostream>
#include <string>

#include "pendulum/experiment_config.h"
#include "pendulum/pendulum_simulator.h"
#include "pendulum/plotting.h"
#include "pendulum/reporting.h"
#include "pendulum/runtime_env.h"

int main(int argc, char* argv[]) {
    configure_pythonpath_from_venv();

    std::cout << "Pendulum Equation Simulator (Runge-Kutta 4)" << std::endl;
    std::cout << "=========================================" << std::endl;

    const std::string config_path = (argc > 1) ? argv[1] : "simple_pendulum.yaml";

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
