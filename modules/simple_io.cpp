#include "modules/simple_io.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

void write_simple_data_file(const std::string& path, const SimulationResult& result) {
    std::ofstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open data file for writing: " + path);
    }

    // Write header matches the unified SimulationResult struct fields as used in damped/driven
    file << "Time,Theta_Analytical,Theta,Omega,Theta_Errors,Omega_Errors,Energy\n";

    // Write data
    for (size_t i = 0; i < result.t.size(); ++i) {
        file << std::fixed << std::setprecision(6)
             << result.t[i] << ","
             << result.theta_analytical[i] << ","
             << result.theta[i] << ","
             << result.omega[i] << ","
             << result.theta_errors[i] << ","
             << result.omega_errors[i] << ","
             << result.energy[i] << "\n";
    }

    file.close();
}
