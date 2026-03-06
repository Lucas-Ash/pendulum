#pragma once

#include <string>

enum class PlottingMethod {
    Original,  // gnuplot workflow
    New        // Python/matplotlib workflow
};

struct DampedPhysicalConfig {
    double g = 9.81;
    double L = 1.0;
    double gamma = 0.3;
    double theta0 = 0.3;
    double theta_dot0 = 0.0;
};

struct DampedSimulationConfig {
    double t_start = 0.0;
    double t_end = 15.0;
    double dt = 0.001;
    int output_every = 10;
};

struct DampedSettingsConfig {
    PlottingMethod plotting_method = PlottingMethod::Original;
    bool show_plot = true;
    bool save_png = true;
    bool plot_phase_map = false;
    bool run_plotter = true;
    std::string data_file = "damped_pendulum_data.dat";
    std::string output_png = "damped_pendulum.png";
    std::string python_script = "plot_damped_pendulum.py";
    std::string integrator = "rk4";
};

struct DampedConfig {
    DampedPhysicalConfig physical;
    DampedSimulationConfig simulation;
    DampedSettingsConfig settings;
};

DampedConfig load_damped_config_from_yaml(const std::string& path);
std::string to_string(PlottingMethod method);
