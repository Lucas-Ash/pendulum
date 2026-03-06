#pragma once

#include <cmath>
#include <string>

namespace restoring_force {

enum class Model {
    Sine,
    Polynomial
};

struct Config {
    Model model = Model::Sine;
    double linear = 1.0;
    double cubic = 0.0;
};

inline bool parse_model_name(const std::string& normalized_name, Model& out) {
    if (normalized_name == "sine" || normalized_name == "sin") {
        out = Model::Sine;
        return true;
    }
    if (normalized_name == "polynomial" || normalized_name == "poly" || normalized_name == "duffing") {
        out = Model::Polynomial;
        return true;
    }
    return false;
}

inline const char* to_string(Model model) {
    return model == Model::Polynomial ? "polynomial" : "sine";
}

inline double linearized_slope(const Config& config) {
    return config.model == Model::Polynomial ? config.linear : 1.0;
}

inline double term(double theta, const Config& config) {
    if (config.model == Model::Polynomial) {
        const double theta_sq = theta * theta;
        return config.linear * theta + config.cubic * theta_sq * theta;
    }
    return std::sin(theta);
}

inline double potential(double theta, const Config& config) {
    if (config.model == Model::Polynomial) {
        const double theta_sq = theta * theta;
        return 0.5 * config.linear * theta_sq + 0.25 * config.cubic * theta_sq * theta_sq;
    }
    return 1.0 - std::cos(theta);
}

}  // namespace restoring_force
