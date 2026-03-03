#include "modules/experiment_config.h"

#include <fstream>
#include <iostream>
#include <stdexcept>

namespace {

std::string trim(const std::string& input) {
    const auto first = input.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return "";
    }
    const auto last = input.find_last_not_of(" \t\r\n");
    return input.substr(first, last - first + 1);
}

bool parse_double(const std::string& text, double& value) {
    size_t parsed = 0;
    value = std::stod(text, &parsed);
    return parsed == text.size();
}

}  // namespace

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

        line = trim(line);
        if (line.empty() || line == "---" || line == "...") {
            continue;
        }

        const auto colon_pos = line.find(':');
        if (colon_pos == std::string::npos) {
            throw std::runtime_error(
                "Invalid config line " + std::to_string(line_number) + " in " + path +
                ": expected key: value");
        }

        const std::string key = trim(line.substr(0, colon_pos));
        const std::string value_text = trim(line.substr(colon_pos + 1));

        if (value_text.empty()) {
            throw std::runtime_error("Missing value for key '" + key + "' at line " +
                                     std::to_string(line_number) + " in " + path);
        }

        double value = 0.0;
        try {
            if (!parse_double(value_text, value)) {
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

    return config;
}
