#pragma once

#include <string>

namespace damping_force {

enum class Model {
    Linear,
    Polynomial
};

inline bool parse_model_name(const std::string& normalized_name, Model& out) {
    if (normalized_name == "linear" || normalized_name == "viscous") {
        out = Model::Linear;
        return true;
    }
    if (normalized_name == "polynomial" || normalized_name == "poly" ||
        normalized_name == "van_der_pol") {
        out = Model::Polynomial;
        return true;
    }
    return false;
}

inline const char* to_string(Model model) {
    return model == Model::Polynomial ? "polynomial" : "linear";
}

inline double coefficient(double theta, Model model, double linear, double cubic) {
    if (model == Model::Polynomial) {
        return linear + cubic * theta * theta;
    }
    return linear;
}

inline double term(double theta,
                   double omega,
                   Model model,
                   double linear,
                   double cubic) {
    return coefficient(theta, model, linear, cubic) * omega;
}

}  // namespace damping_force
