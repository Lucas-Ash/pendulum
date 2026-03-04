#pragma once

#include <string>

struct ExperimentConfig {
    double length = 1.0;
    double gravity = 9.81;
    double dt = 0.01;
    double t_max = 20.0;
    double theta0 = 0.5;
    double omega0 = 0.0;
    std::string data_file = "simple_pendulum_data.csv";
};

ExperimentConfig load_config_from_yaml(const std::string& path);
