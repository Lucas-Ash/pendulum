#include "modules/config/experiment_config.h"

#include <fstream>
#include <iostream>
#include <stdexcept>

#include "modules/core/config_utils.h"

ExperimentConfig load_config_from_yaml(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open config file: " + path);
    }

    ExperimentConfig config;
    std::string line;
    std::string active_section;
    int line_number = 0;

    while (std::getline(file, line)) {
        ++line_number;
        const auto comment_pos = line.find('#');
        if (comment_pos != std::string::npos) {
            line = line.substr(0, comment_pos);
        }

        const std::string trimmed = config_utils::trim(line);
        if (trimmed.empty() || trimmed == "---" || trimmed == "...") {
            continue;
        }

        size_t indent = 0;
        while (indent < line.size() && (line[indent] == ' ' || line[indent] == '\t')) {
            ++indent;
        }

        const auto colon_pos = trimmed.find(':');
        if (colon_pos == std::string::npos) {
            throw std::runtime_error(
                "Invalid config line " + std::to_string(line_number) + " in " + path +
                ": expected key: value");
        }

        const std::string key = config_utils::trim(trimmed.substr(0, colon_pos));
        const std::string value_text = config_utils::trim(trimmed.substr(colon_pos + 1));

        if (value_text.empty()) {
            const bool is_known_section =
                key == "physical" || key == "simulation" ||
                key == "settings" || key == "additional_terms";
            if (!is_known_section) {
                throw std::runtime_error("Missing value for key '" + key + "' at line " +
                                         std::to_string(line_number) + " in " + path);
            }
            active_section = key;
            continue;
        }

        if (indent == 0) {
            active_section.clear();
        }

        const std::string full_key =
            active_section.empty() ? key : (active_section + "." + key);

        if (full_key == "data_file" || full_key == "settings.data_file") {
            config.data_file = value_text;
            // Remove surround quotes if any
            if (config.data_file.size() >= 2 && 
                ((config.data_file.front() == '"' && config.data_file.back() == '"') ||
                 (config.data_file.front() == '\'' && config.data_file.back() == '\''))) {
                config.data_file = config.data_file.substr(1, config.data_file.length() - 2);
            }
            continue;
        } else if (full_key == "output_png" || full_key == "settings.output_png") {
            config.output_png = value_text;
            if (config.output_png.size() >= 2 && 
                ((config.output_png.front() == '"' && config.output_png.back() == '"') ||
                 (config.output_png.front() == '\'' && config.output_png.back() == '\''))) {
                config.output_png = config.output_png.substr(1, config.output_png.length() - 2);
            }
            continue;
        } else if (full_key == "show_plot" || full_key == "settings.show_plot") {
            bool bval = false;
            if (config_utils::parse_bool(value_text, bval)) {
                config.show_plot = bval;
            } else {
                throw std::runtime_error("Invalid boolean for show_plot");
            }
            continue;
        } else if (full_key == "save_png" || full_key == "settings.save_png") {
            bool bval = false;
            if (config_utils::parse_bool(value_text, bval)) {
                config.save_png = bval;
            } else {
                throw std::runtime_error("Invalid boolean for save_png");
            }
            continue;
        } else if (full_key == "plot_phase_map" || full_key == "settings.plot_phase_map") {
            bool bval = false;
            if (config_utils::parse_bool(value_text, bval)) {
                config.plot_phase_map = bval;
            } else {
                throw std::runtime_error("Invalid boolean for plot_phase_map");
            }
            continue;
        } else if (full_key == "integrator" || full_key == "settings.integrator") {
            config.integrator = value_text;
            if (config.integrator.size() >= 2 && 
                ((config.integrator.front() == '"' && config.integrator.back() == '"') ||
                 (config.integrator.front() == '\'' && config.integrator.back() == '\''))) {
                config.integrator = config.integrator.substr(1, config.integrator.length() - 2);
            }
            config.integrator = config_utils::to_lower(config.integrator);
            continue;
        } else if (full_key == "analytical_model" || full_key == "settings.analytical_model") {
            config.analytical_model = value_text;
            if (config.analytical_model.size() >= 2 && 
                ((config.analytical_model.front() == '"' && config.analytical_model.back() == '"') ||
                 (config.analytical_model.front() == '\'' && config.analytical_model.back() == '\''))) {
                config.analytical_model = config.analytical_model.substr(1, config.analytical_model.length() - 2);
            }
            config.analytical_model = config_utils::to_lower(config.analytical_model);
            continue;
        } else if (full_key == "error_analysis" || full_key == "error_mode" ||
                   full_key == "settings.error_analysis" || full_key == "settings.error_mode") {
            const std::string mode_name = config_utils::to_lower(config_utils::parse_string(value_text));
            error_reference::Mode mode = error_reference::Mode::Analytical;
            if (!error_reference::parse_mode_name(mode_name, mode)) {
                throw std::runtime_error(
                    "Invalid error analysis mode '" + mode_name +
                    "'. Use 'analytical', 'none', or 'hd_reference'.");
            }
            config.error_mode = mode;
            continue;
        } else if (full_key == "error_reference_factor" || full_key == "reference_factor" ||
                   full_key == "settings.error_reference_factor" ||
                   full_key == "settings.reference_factor") {
            int factor = 0;
            if (!config_utils::parse_int(value_text, factor)) {
                throw std::runtime_error(
                    "Invalid integer for error_reference_factor");
            }
            config.error_reference_factor = factor;
            continue;
        } else if (full_key == "restoring_force_model" || full_key == "restoring_model" ||
                   full_key == "physical.restoring_force_model" ||
                   full_key == "physical.restoring_model") {
            const std::string model_name = config_utils::to_lower(config_utils::parse_string(value_text));
            restoring_force::Model model = restoring_force::Model::Sine;
            if (!restoring_force::parse_model_name(model_name, model)) {
                throw std::runtime_error(
                    "Invalid restoring force model '" + model_name +
                    "'. Use 'sine' or 'polynomial'.");
            }
            config.restoring_force.model = model;
            continue;
        } else if (full_key == "additional_terms.inverse_cubic_enabled" ||
                   full_key == "inverse_cubic_enabled") {
            bool bval = false;
            if (!config_utils::parse_bool(value_text, bval)) {
                throw std::runtime_error("Invalid boolean for inverse_cubic_enabled");
            }
            config.additional_terms.inverse_cubic_enabled = bval;
            continue;
        } else if (full_key == "additional_terms.exponential_enabled" ||
                   full_key == "exponential_enabled") {
            bool bval = false;
            if (!config_utils::parse_bool(value_text, bval)) {
                throw std::runtime_error("Invalid boolean for exponential_enabled");
            }
            config.additional_terms.exponential_enabled = bval;
            continue;
        } else if (full_key == "additional_terms.exponential_subtract_equilibrium" ||
                   full_key == "exponential_subtract_equilibrium") {
            bool bval = false;
            if (!config_utils::parse_bool(value_text, bval)) {
                throw std::runtime_error(
                    "Invalid boolean for exponential_subtract_equilibrium");
            }
            config.additional_terms.exponential_subtract_equilibrium = bval;
            continue;
        } else if (full_key == "additional_terms.state_power_enabled" ||
                   full_key == "state_power_enabled") {
            bool bval = false;
            if (!config_utils::parse_bool(value_text, bval)) {
                throw std::runtime_error("Invalid boolean for state_power_enabled");
            }
            config.additional_terms.state_power_enabled = bval;
            continue;
        } else if (full_key == "additional_terms.time_damping_enabled" ||
                   full_key == "time_damping_enabled") {
            bool bval = false;
            if (!config_utils::parse_bool(value_text, bval)) {
                throw std::runtime_error("Invalid boolean for time_damping_enabled");
            }
            config.additional_terms.time_damping_enabled = bval;
            continue;
        } else if (full_key == "additional_terms.state_power_mode" ||
                   full_key == "state_power_mode") {
            const std::string mode_name =
                config_utils::to_lower(config_utils::parse_string(value_text));
            additional_terms::StatePowerMode mode =
                additional_terms::StatePowerMode::PositiveOnly;
            if (!additional_terms::parse_state_power_mode(mode_name, mode)) {
                throw std::runtime_error(
                    "Invalid state_power_mode. Use 'positive_only' or 'signed'.");
            }
            config.additional_terms.state_power_mode = mode;
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

        if (full_key == "length" || full_key == "physical.length") {
            config.length = value;
        } else if (full_key == "gravity" || full_key == "physical.gravity") {
            config.gravity = value;
        } else if (full_key == "dt" || full_key == "simulation.dt") {
            config.dt = value;
        } else if (full_key == "t_max" || full_key == "simulation.t_max" ||
                   full_key == "simulation.t_end") {
            config.t_max = value;
        } else if (full_key == "theta0" || full_key == "physical.theta0") {
            config.theta0 = value;
        } else if (full_key == "omega0" || full_key == "physical.omega0") {
            config.omega0 = value;
        } else if (full_key == "restoring_force_linear" || full_key == "restoring_linear" ||
                   full_key == "physical.restoring_force_linear" ||
                   full_key == "physical.restoring_linear") {
            config.restoring_force.linear = value;
        } else if (full_key == "restoring_force_cubic" || full_key == "restoring_cubic" ||
                   full_key == "physical.restoring_force_cubic" ||
                   full_key == "physical.restoring_cubic") {
            config.restoring_force.cubic = value;
        } else if (full_key == "additional_terms.inverse_cubic_strength" ||
                   full_key == "inverse_cubic_strength") {
            config.additional_terms.inverse_cubic_strength = value;
        } else if (full_key == "additional_terms.exponential_strength" ||
                   full_key == "exponential_strength") {
            config.additional_terms.exponential_strength = value;
        } else if (full_key == "additional_terms.exponential_scale" ||
                   full_key == "exponential_scale") {
            config.additional_terms.exponential_scale = value;
        } else if (full_key == "additional_terms.state_power_strength" ||
                   full_key == "state_power_strength") {
            config.additional_terms.state_power_strength = value;
        } else if (full_key == "additional_terms.state_power_exponent" ||
                   full_key == "state_power_exponent") {
            config.additional_terms.state_power_exponent = value;
        } else if (full_key == "additional_terms.time_damping_coefficient" ||
                   full_key == "time_damping_coefficient") {
            config.additional_terms.time_damping_coefficient = value;
        } else if (full_key == "additional_terms.time_damping_power" ||
                   full_key == "time_damping_power") {
            config.additional_terms.time_damping_power = value;
        } else if (full_key == "additional_terms.time_damping_shift" ||
                   full_key == "time_damping_shift") {
            config.additional_terms.time_damping_shift = value;
        } else if (full_key == "additional_terms.singularity_epsilon" ||
                   full_key == "singularity_epsilon") {
            config.additional_terms.singularity_epsilon = value;
        } else {
            std::cerr << "Warning: unknown config key '" << full_key << "' ignored." << std::endl;
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
    if (config.error_reference_factor <= 0) {
        throw std::runtime_error("Invalid config: error_reference_factor must be > 0");
    }
    if (config.error_mode == error_reference::Mode::HdReference &&
        config.error_reference_factor < 2) {
        throw std::runtime_error(
            "Invalid config: error_reference_factor must be >= 2 for hd_reference mode");
    }
    if (config.additional_terms.exponential_enabled &&
        std::abs(config.additional_terms.exponential_scale) <= 1e-15) {
        throw std::runtime_error(
            "Invalid config: additional_terms.exponential_scale must be non-zero");
    }
    if (config.additional_terms.singularity_epsilon <= 0.0) {
        throw std::runtime_error(
            "Invalid config: additional_terms.singularity_epsilon must be > 0");
    }

    config.data_file = config_utils::resolve_output_path(config.data_file);
    config.output_png = config_utils::resolve_output_path(config.output_png);

    return config;
}
