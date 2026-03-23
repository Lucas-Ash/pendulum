#include "modules/driven/driven_io.h"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

void write_driven_data_file(const std::string& path, const SimulationResult& result) {
    const std::filesystem::path output_path(path);
    if (output_path.has_parent_path()) {
        std::filesystem::create_directories(output_path.parent_path());
    }

    std::ofstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open data file for writing: " + path);
    }

    // Write header
    file << "Time,Theta_Analytical,Theta,Omega,Theta_Errors,Omega_Errors,Energy\n";

    file << std::scientific << std::setprecision(6);

    // Write data
    for (size_t i = 0; i < result.t.size(); ++i) {
        file << result.t[i] << ","
             << result.theta_analytical[i] << ","
             << result.theta[i] << ","
             << result.omega[i] << ","
             << result.theta_errors[i] << ","
             << result.omega_errors[i] << ","
             << result.energy[i] << "\n";
    }

    file.close();
}

void write_driven_sweep_data_file(const std::string& path, const DrivenSweepResult& result) {
    const std::filesystem::path output_path(path);
    if (output_path.has_parent_path()) {
        std::filesystem::create_directories(output_path.parent_path());
    }

    std::ofstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open sweep data file for writing: " + path);
    }

    file << "DriveFrequency,NumericalAmplitude,AnalyticalAmplitude,AnalyticalLowerStable,AnalyticalUpperStable,FinalTheta,FinalOmega\n";
    file << std::scientific << std::setprecision(6);

    for (const auto& sample : result.samples) {
        file << sample.drive_frequency << ","
             << sample.numerical_amplitude << ","
             << sample.analytical_amplitude << ","
             << sample.analytical_lower_stable_amplitude << ","
             << sample.analytical_upper_stable_amplitude << ","
             << sample.final_theta << ","
             << sample.final_omega << "\n";
    }

    file.close();
}
