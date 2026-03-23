#pragma once

#include <string>

namespace error_reference {

enum class Mode {
    Analytical,
    None,
    HdReference
};

inline bool parse_mode_name(const std::string& normalized_name, Mode& out) {
    if (normalized_name == "analytical" || normalized_name == "analytic") {
        out = Mode::Analytical;
        return true;
    }
    if (normalized_name == "none" || normalized_name == "off" ||
        normalized_name == "disabled") {
        out = Mode::None;
        return true;
    }
    if (normalized_name == "hd_reference" ||
        normalized_name == "high_definition" ||
        normalized_name == "numerical_reference" ||
        normalized_name == "reference_numerical") {
        out = Mode::HdReference;
        return true;
    }
    return false;
}

inline const char* to_string(Mode mode) {
    if (mode == Mode::None) {
        return "none";
    }
    if (mode == Mode::HdReference) {
        return "hd_reference";
    }
    return "analytical";
}

}  // namespace error_reference
