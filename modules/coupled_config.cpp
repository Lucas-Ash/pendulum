#include "modules/coupled_config.h"
#include "modules/config_utils.h"

#include <fstream>
#include <iostream>
#include <stdexcept>

namespace {
void set_value(CoupledConfig& config, const std::string& key, const std::string& value_text, int line_number, const std::string& path) {
    try {
        double dval = 0.0;
        int ival = 0;

        // Physical
        if (key == "physical.omega_1" || key == "omega_1") {
            if (!config_utils::parse_double(value_text, dval)) throw std::runtime_error("numeric");
            config.physical.omega_1 = dval;
        } else if (key == "physical.zeta_1" || key == "zeta_1") {
            if (!config_utils::parse_double(value_text, dval)) throw std::runtime_error("numeric");
            config.physical.zeta_1 = dval;
        } else if (key == "physical.alpha_11" || key == "alpha_11") {
            if (!config_utils::parse_double(value_text, dval)) throw std::runtime_error("numeric");
            config.physical.alpha_11 = dval;
        } else if (key == "physical.alpha_12" || key == "alpha_12") {
            if (!config_utils::parse_double(value_text, dval)) throw std::runtime_error("numeric");
            config.physical.alpha_12 = dval;
        } else if (key == "physical.omega_2" || key == "omega_2") {
            if (!config_utils::parse_double(value_text, dval)) throw std::runtime_error("numeric");
            config.physical.omega_2 = dval;
        } else if (key == "physical.zeta_2" || key == "zeta_2") {
            if (!config_utils::parse_double(value_text, dval)) throw std::runtime_error("numeric");
            config.physical.zeta_2 = dval;
        } else if (key == "physical.alpha_22" || key == "alpha_22") {
            if (!config_utils::parse_double(value_text, dval)) throw std::runtime_error("numeric");
            config.physical.alpha_22 = dval;
        } else if (key == "physical.alpha_21" || key == "alpha_21") {
            if (!config_utils::parse_double(value_text, dval)) throw std::runtime_error("numeric");
            config.physical.alpha_21 = dval;
        } else if (key == "physical.F" || key == "F") {
            if (!config_utils::parse_double(value_text, dval)) throw std::runtime_error("numeric");
            config.physical.F = dval;
        } else if (key == "physical.omega" || key == "omega" || key == "physical.omega_drive" || key == "omega_drive") {
            if (!config_utils::parse_double(value_text, dval)) throw std::runtime_error("numeric");
            config.physical.omega_drive = dval;
        } else if (key == "physical.q1_0" || key == "q1_0") {
            if (!config_utils::parse_double(value_text, dval)) throw std::runtime_error("numeric");
            config.physical.q1_0 = dval;
        } else if (key == "physical.omega1_0" || key == "omega1_0" || key == "physical.q1_dot0" || key == "q1_dot0") {
            if (!config_utils::parse_double(value_text, dval)) throw std::runtime_error("numeric");
            config.physical.omega1_0 = dval;
        } else if (key == "physical.q2_0" || key == "q2_0") {
            if (!config_utils::parse_double(value_text, dval)) throw std::runtime_error("numeric");
            config.physical.q2_0 = dval;
        } else if (key == "physical.omega2_0" || key == "omega2_0" || key == "physical.q2_dot0" || key == "q2_dot0") {
            if (!config_utils::parse_double(value_text, dval)) throw std::runtime_error("numeric");
            config.physical.omega2_0 = dval;
        }
        // Simulation
        else if (key == "simulation.t_start" || key == "t_start") {
            if (!config_utils::parse_double(value_text, dval)) throw std::runtime_error("numeric");
            config.simulation.t_start = dval;
        } else if (key == "simulation.t_end" || key == "t_end") {
            if (!config_utils::parse_double(value_text, dval)) throw std::runtime_error("numeric");
            config.simulation.t_end = dval;
        } else if (key == "simulation.dt" || key == "dt") {
            if (!config_utils::parse_double(value_text, dval)) throw std::runtime_error("numeric");
            config.simulation.dt = dval;
        } else if (key == "simulation.output_every" || key == "output_every") {
            if (!config_utils::parse_int(value_text, ival)) throw std::runtime_error("integer");
            config.simulation.output_every = ival;
        }
        // Settings
        else if (key == "settings.integrator" || key == "integrator") {
            config.settings.integrator = config_utils::to_lower(config_utils::parse_string(value_text));
        } else if (key == "settings.data_file" || key == "data_file") {
            config.settings.data_file = config_utils::parse_string(value_text);
        } else if (key == "settings.output_png" || key == "output_png") {
            config.settings.output_png = config_utils::parse_string(value_text);
        } else if (key == "settings.python_script" || key == "python_script") {
            config.settings.python_script = config_utils::parse_string(value_text);
        } else if (key == "settings.error_analysis" || key == "error_analysis" ||
                   key == "settings.error_mode" || key == "error_mode") {
            const std::string mode_name = config_utils::to_lower(config_utils::parse_string(value_text));
            error_reference::Mode mode = error_reference::Mode::None;
            if (!error_reference::parse_mode_name(mode_name, mode)) {
                throw std::runtime_error(
                    "error_analysis must be 'analytical', 'none', or 'hd_reference'");
            }
            config.settings.error_mode = mode;
        } else if (key == "settings.error_reference_factor" || key == "error_reference_factor" ||
                   key == "settings.reference_factor" || key == "reference_factor") {
            if (!config_utils::parse_int(value_text, ival)) throw std::runtime_error("integer");
            config.settings.error_reference_factor = ival;
        } else {
            // Ignore unknown values
        }
    } catch (const std::exception&) {
        throw std::runtime_error("Invalid value for key '" + key + "' at line " + std::to_string(line_number) + " in " + path);
    }
}
} // namespace

CoupledConfig load_coupled_config_from_yaml(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open config file: " + path);
    }

    CoupledConfig config;
    std::string line;
    std::string active_section;
    int line_number = 0;

    while (std::getline(file, line)) {
        ++line_number;

        const auto comment_pos = line.find('#');
        if (comment_pos != std::string::npos) {
            line = line.substr(0, comment_pos);
        }

        if (config_utils::trim(line).empty() || config_utils::trim(line) == "---" || config_utils::trim(line) == "...") {
            continue;
        }

        size_t indent = 0;
        while (indent < line.size() && (line[indent] == ' ' || line[indent] == '\t')) {
            ++indent;
        }

        const std::string stripped = config_utils::trim(line);
        const auto colon_pos = stripped.find(':');
        if (colon_pos == std::string::npos) {
            continue; // Skip invalid lines, or maybe throw. Current codebase throws on invalid, so let's throw.
            // Wait, previous modules throw:
            // throw std::runtime_error("Invalid config line " + std::to_string(line_number) + " in " + path + ": expected key: value");
        }

        const std::string key = config_utils::trim(stripped.substr(0, colon_pos));
        const std::string value = config_utils::trim(stripped.substr(colon_pos + 1));

        if (value.empty()) {
            active_section = key;
            continue;
        }

        if (indent == 0) {
            active_section.clear();
        }

        const std::string full_key = active_section.empty() ? key : (active_section + "." + key);
        set_value(config, full_key, value, line_number, path);
    }

    if (config.settings.error_reference_factor <= 0) {
        throw std::runtime_error("Invalid config: error_reference_factor must be > 0");
    }
    if (config.settings.error_mode == error_reference::Mode::HdReference &&
        config.settings.error_reference_factor < 2) {
        throw std::runtime_error(
            "Invalid config: settings.error_reference_factor must be >= 2 for hd_reference mode");
    }

    config.settings.data_file = config_utils::resolve_output_path(config.settings.data_file);
    config.settings.output_png = config_utils::resolve_output_path(config.settings.output_png);
    config.settings.python_script = config_utils::resolve_output_path(config.settings.python_script);

    return config;
}
