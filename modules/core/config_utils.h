#pragma once

#include <string>

namespace config_utils {

std::string trim(const std::string& input);
std::string to_lower(std::string text);
bool parse_double(const std::string& text, double& value);
bool parse_int(const std::string& text, int& value);
bool parse_bool(const std::string& text, bool& value);
std::string parse_string(std::string value);
std::string resolve_output_path(const std::string& original_path);

} // namespace config_utils
