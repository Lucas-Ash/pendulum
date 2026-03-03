#include "pendulum/driven_reporting.h"
#include <iostream>
#include <iomanip>

void print_driven_simulation_summary(const DrivenConfig& config, const DrivenSimulationResult& result) {
    std::cout << "\nSimulation complete!\n";
    std::cout << "Total steps: " << result.rk4_steps << "\n";
    std::cout << "Output points: " << result.samples.size() << "\n\n";
    
    std::cout << "Comparison Statistics:\n";
    std::cout << "  Maximum difference: " << std::scientific << std::setprecision(6) << result.max_diff << " rad\n";
    std::cout << "  Average difference: " << std::fixed << std::setprecision(6) << result.avg_diff << " rad\n\n";
    
    std::cout << "Note: The analytical solution uses small-angle approximation (sin(theta) approx theta).\n";
    std::cout << "For small angles, the solutions agree closely. For larger angles, they diverge.\n";
}
