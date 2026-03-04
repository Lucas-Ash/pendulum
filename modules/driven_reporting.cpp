#include "modules/driven_reporting.h"
#include <iostream>
#include <iomanip>

void print_driven_simulation_summary(const DrivenConfig& config, const SimulationResult& result) {
    std::cout << "\nSimulation complete!\n";
    std::cout << "Total steps: " << result.rk4_steps << "\n";
    std::cout << "Output points: " << result.t.size() << "\n\n";
    
    std::cout << "Comparison Statistics:\n";
    std::cout << "  Maximum difference: " << std::scientific << std::setprecision(6) << result.theta_stats.max_abs << " rad\n";
    std::cout << "  Average difference: " << std::fixed << std::setprecision(6) << result.theta_stats.avg_abs << " rad\n\n";
    
    std::cout << "Note: The analytical solution uses small-angle approximation (sin(theta) approx theta).\n";
    std::cout << "For small angles, the solutions agree closely. For larger angles, they diverge.\n";
}
