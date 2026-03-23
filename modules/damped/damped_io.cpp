#include "modules/damped/damped_io.h"

#include <fstream>
#include <iomanip>
#include <stdexcept>

void write_damped_data_file(const std::string& path, const SimulationResult& result) {
    std::ofstream fout(path);
    if (!fout.is_open()) {
        throw std::runtime_error("Could not open output data file: " + path);
    }

    fout << std::scientific << std::setprecision(10);
    fout << "# t   theta_analytical   theta   omega   "
         << "theta_errors   omega_errors   energy\n";

    for (size_t i = 0; i < result.t.size(); ++i) {
        fout << result.t[i] << "  "
             << result.theta_analytical[i] << "  "
             << result.theta[i] << "  "
             << result.omega[i] << "  "
             << result.theta_errors[i] << "  "
             << result.omega_errors[i] << "  "
             << result.energy[i] << "\n";
    }
}
