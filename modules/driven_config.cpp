#include "modules/driven_config.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include "modules/config_utils.h"

namespace {

DrivenSystemModel parse_system_model(const std::string& value_raw) {
    const std::string value = config_utils::to_lower(config_utils::parse_string(value_raw));
    if (value == "pendulum") return DrivenSystemModel::Pendulum;
    if (value == "duffing") return DrivenSystemModel::Duffing;
    throw std::runtime_error("Invalid system_model: " + value_raw);
}

DrivenPlottingMethod parse_plotting_method(const std::string& value_raw) {
    const std::string value = config_utils::to_lower(config_utils::parse_string(value_raw));
    if (value == "original" || value == "gnuplot") return DrivenPlottingMethod::Original;
    if (value == "new" || value == "python" || value == "matplotlib") return DrivenPlottingMethod::New;
    throw std::runtime_error("Invalid plotting_method: " + value_raw);
}

DrivenSweepDirection parse_sweep_direction(const std::string& value_raw) {
    const std::string value = config_utils::to_lower(config_utils::parse_string(value_raw));
    if (value == "ascending" || value == "forward" || value == "up") {
        return DrivenSweepDirection::Ascending;
    }
    if (value == "descending" || value == "backward" || value == "down") {
        return DrivenSweepDirection::Descending;
    }
    throw std::runtime_error("Invalid sweep direction: " + value_raw);
}

void set_value(DrivenConfig& config, const std::string& key, const std::string& value_text,
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
        if (key == "physical.system_model" || key == "physical.mode" ||
            key == "system_model" || key == "mode") {
            config.physical.system_model = parse_system_model(value_text);
            return;
        }
        if (key == "physical.L" || key == "L") {
            if (!config_utils::parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.physical.L = double_value;
            return;
        }
        if (key == "physical.damping" || key == "damping") {
            if (!config_utils::parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.physical.damping = double_value;
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
            config.physical.damping = double_value;
            return;
        }
        if (key == "physical.damping_cubic" || key == "damping_cubic") {
            if (!config_utils::parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.physical.damping_cubic = double_value;
            return;
        }
        if (key == "physical.A" || key == "A") {
            if (!config_utils::parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.physical.A = double_value;
            return;
        }
        if (key == "physical.mass" || key == "mass" ||
            key == "duffing.mass") {
            if (!config_utils::parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.physical.mass = double_value;
            return;
        }
        if (key == "physical.linear_stiffness" || key == "linear_stiffness" ||
            key == "duffing.linear_stiffness") {
            if (!config_utils::parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.physical.linear_stiffness = double_value;
            return;
        }
        if (key == "physical.cubic_stiffness" || key == "cubic_stiffness" ||
            key == "duffing.cubic_stiffness") {
            if (!config_utils::parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.physical.cubic_stiffness = double_value;
            return;
        }
        if (key == "physical.drive_force" || key == "drive_force" ||
            key == "duffing.drive_force") {
            if (!config_utils::parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.physical.drive_force = double_value;
            return;
        }
        if (key == "physical.omega_drive" || key == "omega_drive") {
            if (!config_utils::parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.physical.omega_drive = double_value;
            return;
        }
        if (key == "physical.drive_phase" || key == "drive_phase") {
            if (!config_utils::parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.physical.drive_phase = double_value;
            return;
        }
        if (key == "physical.parametric_amplitude" || key == "parametric_amplitude") {
            if (!config_utils::parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.physical.parametric_amplitude = double_value;
            return;
        }
        if (key == "physical.parametric_frequency" || key == "parametric_frequency") {
            if (!config_utils::parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.physical.parametric_frequency = double_value;
            return;
        }
        if (key == "physical.theta0" || key == "theta0") {
            if (!config_utils::parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.physical.theta0 = double_value;
            return;
        }
        if (key == "physical.omega0" || key == "omega0") {
            if (!config_utils::parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.physical.omega0 = double_value;
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
        if (key == "unit_scales.enabled") {
            if (!config_utils::parse_bool(value_text, bool_value)) throw std::runtime_error("boolean");
            config.unit_scales.enabled = bool_value;
            return;
        }
        if (key == "unit_scales.time_scale") {
            if (!config_utils::parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.unit_scales.time_scale = double_value;
            return;
        }
        if (key == "unit_scales.displacement_scale") {
            if (!config_utils::parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.unit_scales.displacement_scale = double_value;
            return;
        }
        if (key == "unit_scales.stiffness_scale") {
            if (!config_utils::parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.unit_scales.stiffness_scale = double_value;
            return;
        }
        if (key == "mass_event.enabled") {
            if (!config_utils::parse_bool(value_text, bool_value)) throw std::runtime_error("boolean");
            config.mass_event.enabled = bool_value;
            return;
        }
        if (key == "mass_event.jump_time") {
            if (!config_utils::parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.mass_event.jump_time = double_value;
            return;
        }
        if (key == "mass_event.delta_mass") {
            if (!config_utils::parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.mass_event.delta_mass = double_value;
            return;
        }
        if (key == "mass_event.disable_drive_after_jump") {
            if (!config_utils::parse_bool(value_text, bool_value)) throw std::runtime_error("boolean");
            config.mass_event.disable_drive_after_jump = bool_value;
            return;
        }
        if (key == "noise.enabled") {
            if (!config_utils::parse_bool(value_text, bool_value)) throw std::runtime_error("boolean");
            config.noise.enabled = bool_value;
            return;
        }
        if (key == "noise.force_stddev") {
            if (!config_utils::parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.noise.force_stddev = double_value;
            return;
        }
        if (key == "noise.seed") {
            if (!config_utils::parse_int(value_text, int_value)) throw std::runtime_error("integer");
            config.noise.seed = static_cast<unsigned long long>(int_value);
            return;
        }
        if (key == "noise.correlation_time") {
            if (!config_utils::parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.noise.correlation_time = double_value;
            return;
        }
        if (key == "sweep.enabled") {
            if (!config_utils::parse_bool(value_text, bool_value)) throw std::runtime_error("boolean");
            config.sweep.enabled = bool_value;
            return;
        }
        if (key == "sweep.omega_start") {
            if (!config_utils::parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.sweep.omega_start = double_value;
            return;
        }
        if (key == "sweep.omega_end") {
            if (!config_utils::parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.sweep.omega_end = double_value;
            return;
        }
        if (key == "sweep.points") {
            if (!config_utils::parse_int(value_text, int_value)) throw std::runtime_error("integer");
            config.sweep.points = int_value;
            return;
        }
        if (key == "sweep.settle_time") {
            if (!config_utils::parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.sweep.settle_time = double_value;
            return;
        }
        if (key == "sweep.direction") {
            config.sweep.direction = parse_sweep_direction(value_text);
            return;
        }
        if (key == "sweep.reuse_final_state") {
            if (!config_utils::parse_bool(value_text, bool_value)) throw std::runtime_error("boolean");
            config.sweep.reuse_final_state = bool_value;
            return;
        }
        if (key == "sweep.analytical_branch_tracking") {
            if (!config_utils::parse_bool(value_text, bool_value)) throw std::runtime_error("boolean");
            config.sweep.analytical_branch_tracking = bool_value;
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
        if (key == "settings.sweep_data_file" || key == "sweep_data_file") {
            config.settings.sweep_data_file = config_utils::parse_string(value_text);
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

std::string to_string(DrivenPlottingMethod method) {
    return method == DrivenPlottingMethod::Original ? "original" : "new";
}

std::string to_string(DrivenSystemModel model) {
    switch (model) {
        case DrivenSystemModel::Pendulum:
            return "pendulum";
        case DrivenSystemModel::Duffing:
            return "duffing";
    }
    return "pendulum";
}

std::string to_string(DrivenSweepDirection direction) {
    return direction == DrivenSweepDirection::Descending ? "descending" : "ascending";
}

DrivenConfig load_driven_config_from_yaml(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open config file: " + path);
    }

    DrivenConfig config;
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

    if (config.physical.system_model == DrivenSystemModel::Pendulum) {
        if (config.physical.g <= 0.0) throw std::runtime_error("Invalid config: physical.g must be > 0");
        if (config.physical.L <= 0.0) throw std::runtime_error("Invalid config: physical.L must be > 0");
    } else {
        if (config.physical.mass <= 0.0) throw std::runtime_error("Invalid config: physical.mass must be > 0");
        if (config.physical.linear_stiffness <= 0.0) {
            throw std::runtime_error("Invalid config: physical.linear_stiffness must be > 0");
        }
    }
    if (config.physical.damping_model == damping_force::Model::Linear &&
        config.physical.damping < 0.0) {
        throw std::runtime_error("Invalid config: physical.damping must be >= 0");
    }
    if (config.simulation.dt <= 0.0) throw std::runtime_error("Invalid config: simulation.dt must be > 0");
    if (config.simulation.t_end <= config.simulation.t_start) throw std::runtime_error("Invalid config: simulation.t_end must be > simulation.t_start");
    if (config.simulation.output_every <= 0) throw std::runtime_error("Invalid config: simulation.output_every must be > 0");
    if (config.settings.error_reference_factor <= 0) throw std::runtime_error("Invalid config: settings.error_reference_factor must be > 0");
    if (config.settings.error_mode == error_reference::Mode::HdReference &&
        config.settings.error_reference_factor < 2) {
        throw std::runtime_error(
            "Invalid config: settings.error_reference_factor must be >= 2 for hd_reference mode");
    }
    if (config.settings.error_mode == error_reference::Mode::Analytical &&
        config.physical.system_model != DrivenSystemModel::Pendulum) {
        throw std::runtime_error(
            "Invalid config: analytical error mode is only available for pendulum mode");
    }
    if (config.unit_scales.time_scale <= 0.0) throw std::runtime_error("Invalid config: unit_scales.time_scale must be > 0");
    if (config.unit_scales.displacement_scale <= 0.0) throw std::runtime_error("Invalid config: unit_scales.displacement_scale must be > 0");
    if (config.unit_scales.stiffness_scale <= 0.0) throw std::runtime_error("Invalid config: unit_scales.stiffness_scale must be > 0");
    if (config.noise.force_stddev < 0.0) throw std::runtime_error("Invalid config: noise.force_stddev must be >= 0");
    if (config.noise.correlation_time < 0.0) throw std::runtime_error("Invalid config: noise.correlation_time must be >= 0");
    if (config.mass_event.enabled &&
        (config.mass_event.jump_time < config.simulation.t_start ||
         config.mass_event.jump_time > config.simulation.t_end)) {
        throw std::runtime_error("Invalid config: mass_event.jump_time must be within the simulated interval");
    }
    if (config.sweep.enabled) {
        if (config.sweep.points < 2) throw std::runtime_error("Invalid config: sweep.points must be >= 2");
        if (config.sweep.settle_time <= 0.0) throw std::runtime_error("Invalid config: sweep.settle_time must be > 0");
    }

    config.settings.data_file = config_utils::resolve_output_path(config.settings.data_file);
    config.settings.sweep_data_file = config_utils::resolve_output_path(config.settings.sweep_data_file);
    config.settings.output_png = config_utils::resolve_output_path(config.settings.output_png);
    config.settings.python_script = config_utils::resolve_output_path(config.settings.python_script);

    return config;
}
