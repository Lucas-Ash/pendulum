#include "modules/config/damped_config.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include "modules/core/config_utils.h"

namespace {

PlottingMethod parse_plotting_method(const std::string& value_raw) {
    const std::string value = config_utils::to_lower(config_utils::parse_string(value_raw));
    if (value == "original" || value == "gnuplot") {
        return PlottingMethod::Original;
    }
    if (value == "new" || value == "python" || value == "matplotlib") {
        return PlottingMethod::New;
    }
    throw std::runtime_error("Invalid plotting_method: " + value_raw +
                             ". Use 'original' or 'new'.");
}

void set_value(DampedConfig& config, const std::string& key, const std::string& value_text,
               int line_number, const std::string& path) {
    try {
        double double_value = 0.0;
        int int_value = 0;
        bool bool_value = false;

        if (key == "physical.g" || key == "g") {
            if (!config_utils::parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.physical.g = double_value;
            return;
        }
        if (key == "physical.L" || key == "L") {
            if (!config_utils::parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.physical.L = double_value;
            return;
        }
        if (key == "physical.gamma" || key == "gamma") {
            if (!config_utils::parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.physical.gamma = double_value;
            return;
        }
        if (key == "physical.damping_model" || key == "damping_model") {
            const std::string model_name = config_utils::to_lower(config_utils::parse_string(value_text));
            damping_force::Model model = damping_force::Model::Linear;
            if (!damping_force::parse_model_name(model_name, model)) {
                throw std::runtime_error(
                    "damping_model must be 'linear' or 'polynomial'");
            }
            config.physical.damping_model = model;
            return;
        }
        if (key == "physical.damping_linear" || key == "damping_linear") {
            if (!config_utils::parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.physical.gamma = 0.5 * double_value;
            return;
        }
        if (key == "physical.damping_cubic" || key == "damping_cubic") {
            if (!config_utils::parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.physical.damping_cubic = double_value;
            return;
        }
        if (key == "physical.theta0" || key == "theta0") {
            if (!config_utils::parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.physical.theta0 = double_value;
            return;
        }
        if (key == "physical.theta_dot0" || key == "theta_dot0") {
            if (!config_utils::parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.physical.theta_dot0 = double_value;
            return;
        }
        if (key == "physical.restoring_force_model" || key == "restoring_force_model" ||
            key == "physical.restoring_model" || key == "restoring_model") {
            const std::string model_name = config_utils::to_lower(config_utils::parse_string(value_text));
            restoring_force::Model model = restoring_force::Model::Sine;
            if (!restoring_force::parse_model_name(model_name, model)) {
                throw std::runtime_error(
                    "restoring_force_model must be 'sine' or 'polynomial'");
            }
            config.physical.restoring_force.model = model;
            return;
        }
        if (key == "physical.restoring_force_linear" || key == "restoring_force_linear" ||
            key == "physical.restoring_linear" || key == "restoring_linear") {
            if (!config_utils::parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.physical.restoring_force.linear = double_value;
            return;
        }
        if (key == "physical.restoring_force_cubic" || key == "restoring_force_cubic" ||
            key == "physical.restoring_cubic" || key == "restoring_cubic") {
            if (!config_utils::parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.physical.restoring_force.cubic = double_value;
            return;
        }
        if (key == "simulation.t_start" || key == "t_start") {
            if (!config_utils::parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.simulation.t_start = double_value;
            return;
        }
        if (key == "simulation.t_end" || key == "t_end") {
            if (!config_utils::parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.simulation.t_end = double_value;
            return;
        }
        if (key == "simulation.dt" || key == "dt") {
            if (!config_utils::parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.simulation.dt = double_value;
            return;
        }
        if (key == "simulation.output_every" || key == "output_every") {
            if (!config_utils::parse_int(value_text, int_value)) throw std::runtime_error("integer");
            config.simulation.output_every = int_value;
            return;
        }
        if (key == "settings.plotting_method" || key == "plotting_method") {
            config.settings.plotting_method = parse_plotting_method(value_text);
            return;
        }
        if (key == "settings.show_plot" || key == "show_plot") {
            if (!config_utils::parse_bool(value_text, bool_value)) throw std::runtime_error("boolean");
            config.settings.show_plot = bool_value;
            return;
        }
        if (key == "settings.save_png" || key == "save_png") {
            if (!config_utils::parse_bool(value_text, bool_value)) throw std::runtime_error("boolean");
            config.settings.save_png = bool_value;
            return;
        }
        if (key == "settings.plot_phase_map" || key == "plot_phase_map") {
            if (!config_utils::parse_bool(value_text, bool_value)) throw std::runtime_error("boolean");
            config.settings.plot_phase_map = bool_value;
            return;
        }
        if (key == "settings.run_plotter" || key == "run_plotter") {
            if (!config_utils::parse_bool(value_text, bool_value)) throw std::runtime_error("boolean");
            config.settings.run_plotter = bool_value;
            return;
        }
        if (key == "settings.data_file" || key == "data_file") {
            config.settings.data_file = config_utils::parse_string(value_text);
            return;
        }
        if (key == "settings.output_png" || key == "output_png") {
            config.settings.output_png = config_utils::parse_string(value_text);
            return;
        }
        if (key == "settings.python_script" || key == "python_script") {
            config.settings.python_script = config_utils::parse_string(value_text);
            return;
        }
        if (key == "settings.integrator" || key == "integrator") {
            config.settings.integrator = config_utils::to_lower(config_utils::parse_string(value_text));
            return;
        }
        if (key == "settings.analytical_model" || key == "analytical_model") {
            config.settings.analytical_model =
                config_utils::to_lower(config_utils::parse_string(value_text));
            return;
        }
        if (key == "settings.error_analysis" || key == "error_analysis" ||
            key == "settings.error_mode" || key == "error_mode") {
            const std::string mode_name = config_utils::to_lower(config_utils::parse_string(value_text));
            error_reference::Mode mode = error_reference::Mode::Analytical;
            if (!error_reference::parse_mode_name(mode_name, mode)) {
                throw std::runtime_error(
                    "error_analysis must be 'analytical', 'none', or 'hd_reference'");
            }
            config.settings.error_mode = mode;
            return;
        }
        if (key == "settings.error_reference_factor" || key == "error_reference_factor" ||
            key == "settings.reference_factor" || key == "reference_factor") {
            if (!config_utils::parse_int(value_text, int_value)) throw std::runtime_error("integer");
            config.settings.error_reference_factor = int_value;
            return;
        }

        std::cerr << "Warning: unknown config key '" << key << "' at line "
                  << line_number << " ignored." << std::endl;
    } catch (const std::exception&) {
        throw std::runtime_error("Invalid value for key '" + key + "' at line " +
                                 std::to_string(line_number) + " in " + path);
    }
}

}  // namespace

std::string to_string(PlottingMethod method) {
    return method == PlottingMethod::Original ? "original" : "new";
}

DampedConfig load_damped_config_from_yaml(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open config file: " + path);
    }

    DampedConfig config;
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
            throw std::runtime_error("Invalid config line " + std::to_string(line_number) +
                                     " in " + path + ": expected key: value");
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

        const std::string full_key =
            active_section.empty() ? key : (active_section + "." + key);
        set_value(config, full_key, value, line_number, path);
    }

    if (config.physical.g <= 0.0) {
        throw std::runtime_error("Invalid config: physical.g must be > 0");
    }
    if (config.physical.L <= 0.0) {
        throw std::runtime_error("Invalid config: physical.L must be > 0");
    }
    if (config.physical.damping_model == damping_force::Model::Linear &&
        config.physical.gamma < 0.0) {
        throw std::runtime_error("Invalid config: physical.gamma must be >= 0");
    }
    if (config.simulation.dt <= 0.0) {
        throw std::runtime_error("Invalid config: simulation.dt must be > 0");
    }
    if (config.simulation.t_end <= config.simulation.t_start) {
        throw std::runtime_error("Invalid config: simulation.t_end must be > simulation.t_start");
    }
    if (config.simulation.output_every <= 0) {
        throw std::runtime_error("Invalid config: simulation.output_every must be > 0");
    }
    if (config.settings.data_file.empty()) {
        throw std::runtime_error("Invalid config: settings.data_file must not be empty");
    }
    if (config.settings.output_png.empty()) {
        throw std::runtime_error("Invalid config: settings.output_png must not be empty");
    }
    if (config.settings.python_script.empty()) {
        throw std::runtime_error("Invalid config: settings.python_script must not be empty");
    }
    if (config.settings.error_reference_factor <= 0) {
        throw std::runtime_error("Invalid config: settings.error_reference_factor must be > 0");
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
