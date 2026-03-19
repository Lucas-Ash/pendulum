#pragma once

#include <string>

struct CoupledPhysicalConfig {
    double omega_1 = 1.0;
    double zeta_1 = 0.05;
    double alpha_11 = 0.0;
    double alpha_12 = 0.0;
    
    double omega_2 = 1.0;
    double zeta_2 = 0.05;
    double alpha_22 = 0.0;
    double alpha_21 = 0.0;

    double F = 0.0;
    double omega_drive = 1.0;
    
    double q1_0 = 0.1;
    double omega1_0 = 0.0;
    double q2_0 = 0.1;
    double omega2_0 = 0.0;
};

struct CoupledSimulationConfig {
    double t_start = 0.0;
    double t_end = 50.0;
    double dt = 0.01;
    int output_every = 1;
};

#include "modules/error_reference.h"

struct CoupledSettingsConfig {
    std::string integrator = "rk4";
    std::string data_file = "outputs/coupled_data.csv";
    std::string output_png = "outputs/coupled_plot.png";
    std::string python_script = "outputs/plot_coupled.py";
    error_reference::Mode error_mode = error_reference::Mode::None;
    int error_reference_factor = 50;
};

struct CoupledConfig {
    CoupledPhysicalConfig physical;
    CoupledSimulationConfig simulation;
    CoupledSettingsConfig settings;
};

CoupledConfig load_coupled_config_from_yaml(const std::string& path);
