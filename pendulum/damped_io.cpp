#include "pendulum/damped_io.h"

#include <fstream>
#include <iomanip>
#include <stdexcept>

void write_damped_data_file(const std::string& path, const DampedSimulationResult& result) {
    std::ofstream fout(path);
    if (!fout.is_open()) {
        throw std::runtime_error("Could not open output data file: " + path);
    }

    fout << std::scientific << std::setprecision(10);
    fout << "# t   analytical   numerical_linear   numerical_nonlinear   "
         << "error_linear   error_nonlinear   energy_nonlinear\n";

    for (const auto& s : result.samples) {
        fout << s.t << "  "
             << s.analytical << "  "
             << s.numerical_linear << "  "
             << s.numerical_nonlinear << "  "
             << s.error_linear << "  "
             << s.error_nonlinear << "  "
             << s.energy_nonlinear << "\n";
    }
}
