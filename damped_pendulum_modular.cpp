#include <exception>
#include <iostream>
#include <string>

#include "pendulum/damped_config.h"
#include "pendulum/damped_io.h"
#include "pendulum/damped_plotting.h"
#include "pendulum/damped_reporting.h"
#include "pendulum/damped_simulator.h"

int main(int argc, char* argv[]) {
    std::cout << "============================================\n";
    std::cout << "   Damped Pendulum Simulation (Modular)\n";
    std::cout << "============================================\n";

    const std::string config_path =
        (argc > 1) ? argv[1] : "damped_pendulum.yaml";

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
