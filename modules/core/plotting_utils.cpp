#include "modules/core/plotting_utils.h"

#include <iostream>
#include <cstdlib>

namespace plotting_utils {

bool run_python_script(const std::string& script_path) {
    const std::string command = "python3 \"" + script_path + "\"";
    const int status = std::system(command.c_str());
    if (status != 0) {
        std::cerr << "Warning: python plotting command failed: " << command << "\n";
        return false;
    }
    return true;
}

} // namespace plotting_utils
