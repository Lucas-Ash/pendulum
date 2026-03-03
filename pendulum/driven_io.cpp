#include "pendulum/driven_io.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

void write_driven_data_file(const std::string& path, const DrivenSimulationResult& result) {
    std::ofstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open data file for writing: " + path);
    }

    // Write header
    file << "Time,Theta_Numerical,Theta_Analytical,Difference\n";

    // Write data
    for (const auto& sample : result.samples) {
        file << std::fixed << std::setprecision(6)
             << sample.t << ","
             << sample.theta_numerical << ","
             << sample.theta_analytical << ","
             << sample.difference << "\n";
    }

    file.close();
}
