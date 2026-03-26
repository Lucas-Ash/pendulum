#pragma once

#include <string>
#include "modules/core/additional_terms.h"
#include "modules/core/error_reference.h"
#include "modules/core/restoring_force.h"

struct ExperimentConfig {
    double length = 1.0;
    double gravity = 9.81;
    double dt = 0.01;
    double t_max = 20.0;
    double theta0 = 0.5;
    double omega0 = 0.0;
    std::string data_file = "simple_pendulum_data.csv";
    bool show_plot = true;
    bool save_png = true;
    bool plot_phase_map = true;
    std::string output_png = "simple_pendulum.png";
    std::string integrator = "rk4";
    std::string analytical_model = "linear";
    error_reference::Mode error_mode = error_reference::Mode::Analytical;
    int error_reference_factor = 50;
    restoring_force::Config restoring_force;
    additional_terms::Config additional_terms;
};

ExperimentConfig load_config_from_yaml(const std::string& path);
