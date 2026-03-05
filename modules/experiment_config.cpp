#include "modules/experiment_config.h"

#include <fstream>
#include <iostream>
#include <stdexcept>

#include "modules/config_utils.h"

ExperimentConfig load_config_from_yaml(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open config file: " + path);
    }

    ExperimentConfig config;
    std::string line;
    int line_number = 0;

    while (std::getline(file, line)) {
        ++line_number;
        const auto comment_pos = line.find('#');
        if (comment_pos != std::string::npos) {
            line = line.substr(0, comment_pos);
        }

        line = config_utils::trim(line);
        if (line.empty() || line == "---" || line == "...") {
            continue;
        }

        const auto colon_pos = line.find(':');
        if (colon_pos == std::string::npos) {
            throw std::runtime_error(
                "Invalid config line " + std::to_string(line_number) + " in " + path +
                ": expected key: value");
        }

        const std::string key = config_utils::trim(line.substr(0, colon_pos));
        const std::string value_text = config_utils::trim(line.substr(colon_pos + 1));

        if (value_text.empty()) {
            throw std::runtime_error("Missing value for key '" + key + "' at line " +
                                     std::to_string(line_number) + " in " + path);
        }

        if (key == "data_file") {
            config.data_file = value_text;
            // Remove surround quotes if any
            if (config.data_file.size() >= 2 && 
                ((config.data_file.front() == '"' && config.data_file.back() == '"') ||
                 (config.data_file.front() == '\'' && config.data_file.back() == '\''))) {
                config.data_file = config.data_file.substr(1, config.data_file.length() - 2);
            }
            continue;
        } else if (key == "output_png") {
            config.output_png = value_text;
            if (config.output_png.size() >= 2 && 
                ((config.output_png.front() == '"' && config.output_png.back() == '"') ||
                 (config.output_png.front() == '\'' && config.output_png.back() == '\''))) {
                config.output_png = config.output_png.substr(1, config.output_png.length() - 2);
            }
            continue;
        } else if (key == "show_plot") {
            bool bval = false;
            if (config_utils::parse_bool(value_text, bval)) {
                config.show_plot = bval;
            } else {
                throw std::runtime_error("Invalid boolean for show_plot");
            }
            continue;
        } else if (key == "save_png") {
            bool bval = false;
            if (config_utils::parse_bool(value_text, bval)) {
                config.save_png = bval;
            } else {
                throw std::runtime_error("Invalid boolean for save_png");
            }
            continue;
        } else if (key == "integrator") {
            config.integrator = value_text;
            if (config.integrator.size() >= 2 && 
                ((config.integrator.front() == '"' && config.integrator.back() == '"') ||
                 (config.integrator.front() == '\'' && config.integrator.back() == '\''))) {
                config.integrator = config.integrator.substr(1, config.integrator.length() - 2);
            }
            config.integrator = config_utils::to_lower(config.integrator);
            continue;
        } else if (key == "analytical_model") {
            config.analytical_model = value_text;
            if (config.analytical_model.size() >= 2 && 
                ((config.analytical_model.front() == '"' && config.analytical_model.back() == '"') ||
                 (config.analytical_model.front() == '\'' && config.analytical_model.back() == '\''))) {
                config.analytical_model = config.analytical_model.substr(1, config.analytical_model.length() - 2);
            }
            config.analytical_model = config_utils::to_lower(config.analytical_model);
            continue;
        }

        double value = 0.0;
        try {
            if (!config_utils::parse_double(value_text, value)) {
                throw std::runtime_error("non-numeric value");
            }
        } catch (...) {
            throw std::runtime_error("Invalid numeric value for key '" + key + "' at line " +
                                     std::to_string(line_number) + " in " + path);
        }

        if (key == "length") {
            config.length = value;
        } else if (key == "gravity") {
            config.gravity = value;
        } else if (key == "dt") {
            config.dt = value;
        } else if (key == "t_max") {
            config.t_max = value;
        } else if (key == "theta0") {
            config.theta0 = value;
        } else if (key == "omega0") {
            config.omega0 = value;
        } else {
            std::cerr << "Warning: unknown config key '" << key << "' ignored." << std::endl;
        }
    }

    if (config.length <= 0.0) {
        throw std::runtime_error("Invalid config: length must be > 0");
    }
    if (config.gravity <= 0.0) {
        throw std::runtime_error("Invalid config: gravity must be > 0");
    }
    if (config.dt <= 0.0) {
        throw std::runtime_error("Invalid config: dt must be > 0");
    }
    if (config.t_max <= 0.0) {
        throw std::runtime_error("Invalid config: t_max must be > 0");
    }

    config.data_file = config_utils::resolve_output_path(config.data_file);
    config.output_png = config_utils::resolve_output_path(config.output_png);

    return config;
}
