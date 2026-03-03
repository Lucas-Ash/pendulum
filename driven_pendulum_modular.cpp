#include <exception>
#include <iostream>
#include <string>

#include "pendulum/driven_config.h"
#include "pendulum/driven_io.h"
#include "pendulum/driven_plotting.h"
#include "pendulum/driven_reporting.h"
#include "pendulum/driven_simulator.h"

int main(int argc, char* argv[]) {
    std::cout << "==============================================\n";
    std::cout << "  Damped Driven Pendulum Simulation\n";
    std::cout << "  Numerical (RK4) vs Analytical Solution\n";
    std::cout << "==============================================\n\n";

    const std::string config_path =
        (argc > 1) ? argv[1] : "driven_pendulum.yaml";

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
