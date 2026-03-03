#include "modules/driven_config.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace {

std::string trim(const std::string& input) {
    const auto first = input.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    const auto last = input.find_last_not_of(" \t\r\n");
    return input.substr(first, last - first + 1);
}

std::string to_lower(std::string text) {
    std::transform(text.begin(), text.end(), text.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return text;
}

bool parse_double(const std::string& text, double& value) {
    size_t parsed = 0;
    value = std::stod(text, &parsed);
    return parsed == text.size();
}

bool parse_int(const std::string& text, int& value) {
    size_t parsed = 0;
    value = std::stoi(text, &parsed);
    return parsed == text.size();
}

bool parse_bool(const std::string& text, bool& value) {
    const std::string lower = to_lower(trim(text));
    if (lower == "true" || lower == "yes" || lower == "on" || lower == "1") { value = true; return true; }
    if (lower == "false" || lower == "no" || lower == "off" || lower == "0") { value = false; return true; }
    return false;
}

std::string parse_string(std::string value) {
    value = trim(value);
    if (value.size() >= 2 && ((value.front() == '"' && value.back() == '"') ||
                              (value.front() == '\'' && value.back() == '\''))) {
        return value.substr(1, value.size() - 2);
    }
    return value;
}

DrivenPlottingMethod parse_plotting_method(const std::string& value_raw) {
    const std::string value = to_lower(parse_string(value_raw));
    if (value == "original" || value == "gnuplot") return DrivenPlottingMethod::Original;
    if (value == "new" || value == "python" || value == "matplotlib") return DrivenPlottingMethod::New;
    throw std::runtime_error("Invalid plotting_method: " + value_raw);
}

void set_value(DrivenConfig& config, const std::string& key, const std::string& value_text,
               int line_number, const std::string& path) {
    try {
        double double_value = 0.0;
        int int_value = 0;
        bool bool_value = false;

        if (key == "physical.g" || key == "g") {
            if (!parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.physical.g = double_value;
            return;
        }
        if (key == "physical.L" || key == "L") {
            if (!parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.physical.L = double_value;
            return;
        }
        if (key == "physical.damping" || key == "damping") {
            if (!parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.physical.damping = double_value;
            return;
        }
        if (key == "physical.A" || key == "A") {
            if (!parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.physical.A = double_value;
            return;
        }
        if (key == "physical.omega_drive" || key == "omega_drive") {
            if (!parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.physical.omega_drive = double_value;
            return;
        }
        if (key == "physical.theta0" || key == "theta0") {
            if (!parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.physical.theta0 = double_value;
            return;
        }
        if (key == "physical.omega0" || key == "omega0") {
            if (!parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.physical.omega0 = double_value;
            return;
        }
        if (key == "simulation.t_start" || key == "t_start") {
            if (!parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.simulation.t_start = double_value;
            return;
        }
        if (key == "simulation.t_end" || key == "t_end") {
            if (!parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.simulation.t_end = double_value;
            return;
        }
        if (key == "simulation.dt" || key == "dt") {
            if (!parse_double(value_text, double_value)) throw std::runtime_error("numeric");
            config.simulation.dt = double_value;
            return;
        }
        if (key == "simulation.output_every" || key == "output_every") {
            if (!parse_int(value_text, int_value)) throw std::runtime_error("integer");
            config.simulation.output_every = int_value;
            return;
        }
        if (key == "settings.plotting_method" || key == "plotting_method") {
            config.settings.plotting_method = parse_plotting_method(value_text);
            return;
        }
        if (key == "settings.show_plot" || key == "show_plot") {
            if (!parse_bool(value_text, bool_value)) throw std::runtime_error("boolean");
            config.settings.show_plot = bool_value;
            return;
        }
        if (key == "settings.save_png" || key == "save_png") {
            if (!parse_bool(value_text, bool_value)) throw std::runtime_error("boolean");
            config.settings.save_png = bool_value;
            return;
        }
        if (key == "settings.run_plotter" || key == "run_plotter") {
            if (!parse_bool(value_text, bool_value)) throw std::runtime_error("boolean");
            config.settings.run_plotter = bool_value;
            return;
        }
        if (key == "settings.data_file" || key == "data_file") {
            config.settings.data_file = parse_string(value_text);
            return;
        }
        if (key == "settings.output_png" || key == "output_png") {
            config.settings.output_png = parse_string(value_text);
            return;
        }
        if (key == "settings.python_script" || key == "python_script") {
            config.settings.python_script = parse_string(value_text);
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

        if (trim(line).empty() || trim(line) == "---" || trim(line) == "...") {
            continue;
        }

        size_t indent = 0;
        while (indent < line.size() && (line[indent] == ' ' || line[indent] == '\t')) {
            ++indent;
        }

        const std::string stripped = trim(line);
        const auto colon_pos = stripped.find(':');
        if (colon_pos == std::string::npos) {
            throw std::runtime_error("Invalid config line " + std::to_string(line_number) +
                                     " in " + path + ": expected key: value");
        }

        const std::string key = trim(stripped.substr(0, colon_pos));
        const std::string value = trim(stripped.substr(colon_pos + 1));

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

    if (config.physical.g <= 0.0) throw std::runtime_error("Invalid config: physical.g must be > 0");
    if (config.physical.L <= 0.0) throw std::runtime_error("Invalid config: physical.L must be > 0");
    if (config.physical.damping < 0.0) throw std::runtime_error("Invalid config: physical.damping must be >= 0");
    if (config.simulation.dt <= 0.0) throw std::runtime_error("Invalid config: simulation.dt must be > 0");
    if (config.simulation.t_end <= config.simulation.t_start) throw std::runtime_error("Invalid config: simulation.t_end must be > simulation.t_start");
    if (config.simulation.output_every <= 0) throw std::runtime_error("Invalid config: simulation.output_every must be > 0");

    return config;
}
