#pragma once

#include <string>

enum class DrivenPlottingMethod {
    Original,  // gnuplot workflow
    New        // Python/matplotlib workflow
};

struct DrivenPhysicalConfig {
    double g = 9.81;
    double L = 1.0;
    double damping = 0.5;
    double A = 0.5;
    double omega_drive = 1.2;
    double theta0 = 0.1;
    double omega0 = 0.0;
};

struct DrivenSimulationConfig {
    double t_start = 0.0;
    double t_end = 50.0;
    double dt = 0.001;
    int output_every = 10;
};

struct DrivenSettingsConfig {
    DrivenPlottingMethod plotting_method = DrivenPlottingMethod::New;
    bool show_plot = true;
    bool save_png = true;
    bool run_plotter = true;
    std::string data_file = "driven_pendulum_data.csv";
    std::string output_png = "driven_pendulum_comparison.png";
    std::string python_script = "plot_driven_pendulum.py";
    std::string integrator = "rk4";
};

struct DrivenConfig {
    DrivenPhysicalConfig physical;
    DrivenSimulationConfig simulation;
    DrivenSettingsConfig settings;
};

DrivenConfig load_driven_config_from_yaml(const std::string& path);
std::string to_string(DrivenPlottingMethod method);
