#pragma once

#include "modules/damping_force.h"
#include <string>
#include "modules/error_reference.h"
#include "modules/restoring_force.h"

enum class DrivenSystemModel {
    Pendulum,
    Duffing
};

enum class DrivenPlottingMethod {
    Original,  // gnuplot workflow
    New        // Python/matplotlib workflow
};

enum class DrivenSweepDirection {
    Ascending,
    Descending
};

struct DrivenPhysicalConfig {
    double g = 9.81;
    double L = 1.0;
    double damping = 0.5;
    damping_force::Model damping_model = damping_force::Model::Linear;
    double damping_cubic = 0.0;
    double A = 0.5;
    double omega_drive = 1.2;
    double theta0 = 0.1;
    double omega0 = 0.0;
    DrivenSystemModel system_model = DrivenSystemModel::Pendulum;
    double mass = 1.0;
    double linear_stiffness = 9.81;
    double cubic_stiffness = 0.0;
    double drive_force = 0.5;
    restoring_force::Config restoring_force;
};

struct DrivenSimulationConfig {
    double t_start = 0.0;
    double t_end = 50.0;
    double dt = 0.001;
    int output_every = 10;
};

struct DrivenUnitScalesConfig {
    bool enabled = false;
    double time_scale = 1.0;
    double displacement_scale = 1.0;
    double stiffness_scale = 1.0;
};

struct DrivenMassEventConfig {
    bool enabled = false;
    double jump_time = 0.0;
    double delta_mass = 0.0;
    bool disable_drive_after_jump = false;
};

struct DrivenNoiseConfig {
    bool enabled = false;
    double force_stddev = 0.0;
    unsigned long long seed = 12345ull;
    double correlation_time = 0.0;
};

struct DrivenSweepConfig {
    bool enabled = false;
    double omega_start = 0.8;
    double omega_end = 1.6;
    int points = 25;
    double settle_time = 20.0;
    DrivenSweepDirection direction = DrivenSweepDirection::Ascending;
    bool reuse_final_state = true;
    bool analytical_branch_tracking = true;
};

struct DrivenSettingsConfig {
    DrivenPlottingMethod plotting_method = DrivenPlottingMethod::New;
    bool show_plot = true;
    bool save_png = true;
    bool plot_phase_map = false;
    bool run_plotter = true;
    std::string data_file = "driven_pendulum_data.csv";
    std::string sweep_data_file = "driven_pendulum_sweep.csv";
    std::string output_png = "driven_pendulum_comparison.png";
    std::string python_script = "plot_driven_pendulum.py";
    std::string integrator = "rk4";
    error_reference::Mode error_mode = error_reference::Mode::Analytical;
    int error_reference_factor = 50;
};

struct DrivenConfig {
    DrivenPhysicalConfig physical;
    DrivenSimulationConfig simulation;
    DrivenUnitScalesConfig unit_scales;
    DrivenMassEventConfig mass_event;
    DrivenNoiseConfig noise;
    DrivenSweepConfig sweep;
    DrivenSettingsConfig settings;
};

DrivenConfig load_driven_config_from_yaml(const std::string& path);
std::string to_string(DrivenPlottingMethod method);
std::string to_string(DrivenSystemModel model);
std::string to_string(DrivenSweepDirection direction);
