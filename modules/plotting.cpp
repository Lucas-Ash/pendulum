#include "modules/plotting.h"

#include <vector>
#include <cstdlib>

#include "modules/add-ons/matplotlibcpp.h"
#include <iostream>
namespace plt = matplotlibcpp;

void plot_simulation_results(const ExperimentConfig& config, const SimulationResult& result) {
    if (std::getenv("QA_TEST")) return;
    
    plt::figure_size(1200, 800);

    plt::subplot2grid(2, 2, 0, 0);
    plt::named_plot("Numerical", result.t, result.theta, "b-");
    plt::named_plot("Equilibrium", result.t, std::vector<double>(result.t.size(), 0.0), "k--");
    plt::title("Pendulum Angle θ(t)");
    plt::ylabel("θ (rad)");
    plt::legend();
    plt::grid(true);

    plt::subplot2grid(2, 2, 0, 1);
    plt::named_plot("Numerical", result.t, result.omega, "r-");
    plt::title("Angular Velocity ω(t)");
    plt::ylabel("ω (rad/s)");
    plt::legend();
    plt::grid(true);

    plt::subplot2grid(2, 2, 1, 0);
    plt::named_plot("Trajectory", result.theta, result.omega, "g-");
    plt::named_plot("Equilibrium", std::vector<double>{0.0}, std::vector<double>{0.0}, "ko");
    plt::title("Phase Portrait (θ, ω)");
    plt::xlabel("θ (rad)");
    plt::ylabel("ω (rad/s)");
    plt::legend();
    plt::grid(true);
    plt::axis("equal");

    plt::subplot2grid(2, 2, 1, 1);
    plt::named_semilogy("θ error", result.t, result.theta_errors, "b-");
    plt::named_semilogy("ω error", result.t, result.omega_errors, "r-");
    plt::title("Absolute Errors (log scale)");
    plt::xlabel("Time t (s)");
    plt::ylabel("Error");
    plt::legend();
    plt::grid(true);

    plt::tight_layout();
    
    if (config.save_png) {
        plt::save(config.output_png);
        std::cout << "Plot saved to " << config.output_png << "\n";
    }
    
    if (config.show_plot) {
        plt::show();
    }
}
