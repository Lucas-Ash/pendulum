#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace additional_terms {

enum class StatePowerMode {
    Signed,
    PositiveOnly
};

struct Config {
    bool inverse_cubic_enabled = false;
    double inverse_cubic_strength = 0.0;

    bool exponential_enabled = false;
    double exponential_strength = 0.0;
    double exponential_scale = 1.0;
    bool exponential_subtract_equilibrium = true;

    bool state_power_enabled = false;
    double state_power_strength = 0.0;
    double state_power_exponent = 1.0;
    StatePowerMode state_power_mode = StatePowerMode::PositiveOnly;

    bool time_damping_enabled = false;
    double time_damping_coefficient = 0.0;
    double time_damping_power = 1.0;
    double time_damping_shift = 0.0;

    double singularity_epsilon = 1e-9;
};

inline bool parse_state_power_mode(const std::string& normalized_name,
                                   StatePowerMode& out) {
    if (normalized_name == "signed" || normalized_name == "odd") {
        out = StatePowerMode::Signed;
        return true;
    }
    if (normalized_name == "positive" || normalized_name == "positive_only" ||
        normalized_name == "unsigned") {
        out = StatePowerMode::PositiveOnly;
        return true;
    }
    return false;
}

inline const char* to_string(StatePowerMode mode) {
    return mode == StatePowerMode::Signed ? "signed" : "positive_only";
}

inline bool has_velocity_dependence(const Config& config) {
    return config.time_damping_enabled &&
           std::abs(config.time_damping_coefficient) > 0.0;
}

inline double state_power_value(double theta, const Config& config) {
    if (!config.state_power_enabled ||
        std::abs(config.state_power_strength) < 1e-30) {
        return 0.0;
    }

    const double exponent = config.state_power_exponent;
    if (config.state_power_mode == StatePowerMode::Signed) {
        const double magnitude = std::pow(std::abs(theta), exponent);
        return config.state_power_strength * std::copysign(magnitude, theta);
    }

    if (theta < 0.0 && std::abs(exponent - std::round(exponent)) > 1e-12) {
        throw std::runtime_error(
            "positive_only state-power term requires theta >= 0 for non-integer exponents.");
    }
    return config.state_power_strength * std::pow(theta, exponent);
}

inline double linearized_stiffness(const Config& config) {
    double out = 0.0;
    if (config.exponential_enabled) {
        out += config.exponential_strength * config.exponential_scale;
    }
    if (config.state_power_enabled &&
        std::abs(config.state_power_exponent - 1.0) < 1e-12) {
        out += config.state_power_strength;
    }
    return out;
}

inline double acceleration(double t,
                           double theta,
                           double omega,
                           const Config& config) {
    double out = 0.0;

    if (config.inverse_cubic_enabled &&
        std::abs(config.inverse_cubic_strength) > 0.0) {
        if (std::abs(theta) <= config.singularity_epsilon) {
            throw std::runtime_error(
                "inverse_cubic term encountered a near-singular state.");
        }
        out += config.inverse_cubic_strength / (theta * theta * theta);
    }

    if (config.exponential_enabled &&
        std::abs(config.exponential_strength) > 0.0) {
        const double exp_term =
            std::exp(config.exponential_scale * theta);
        out -= config.exponential_strength *
               (exp_term -
                (config.exponential_subtract_equilibrium ? 1.0 : 0.0));
    }

    if (config.state_power_enabled &&
        std::abs(config.state_power_strength) > 0.0) {
        out -= state_power_value(theta, config);
    }

    if (config.time_damping_enabled &&
        std::abs(config.time_damping_coefficient) > 0.0) {
        const double shifted = t + config.time_damping_shift;
        if (shifted <= config.singularity_epsilon) {
            throw std::runtime_error(
                "time_damping term is singular; use t_start + time_damping_shift > 0.");
        }
        out -= config.time_damping_coefficient /
               std::pow(shifted, config.time_damping_power) * omega;
    }

    return out;
}

inline double potential(double theta, const Config& config) {
    double out = 0.0;

    if (config.inverse_cubic_enabled &&
        std::abs(config.inverse_cubic_strength) > 0.0) {
        if (std::abs(theta) <= config.singularity_epsilon) {
            throw std::runtime_error(
                "inverse_cubic term encountered a near-singular state.");
        }
        out += 0.5 * config.inverse_cubic_strength / (theta * theta);
    }

    if (config.exponential_enabled &&
        std::abs(config.exponential_strength) > 0.0) {
        if (std::abs(config.exponential_scale) <= 1e-15) {
            throw std::runtime_error(
                "exponential_scale must be non-zero when exponential terms are enabled.");
        }
        out += config.exponential_strength *
               std::exp(config.exponential_scale * theta) /
               config.exponential_scale;
        if (config.exponential_subtract_equilibrium) {
            out -= config.exponential_strength * theta;
        }
    }

    if (config.state_power_enabled &&
        std::abs(config.state_power_strength) > 0.0) {
        const double denom = config.state_power_exponent + 1.0;
        if (std::abs(denom) <= 1e-15) {
            throw std::runtime_error(
                "state_power_exponent = -1 is not supported in potential evaluation.");
        }

        if (config.state_power_mode == StatePowerMode::Signed) {
            out += config.state_power_strength *
                   std::pow(std::abs(theta), denom) / denom;
        } else if (theta >= 0.0 ||
                   std::abs(config.state_power_exponent -
                            std::round(config.state_power_exponent)) < 1e-12) {
            out += config.state_power_strength *
                   std::pow(theta, denom) / denom;
        }
    }

    return out;
}

}  // namespace additional_terms
