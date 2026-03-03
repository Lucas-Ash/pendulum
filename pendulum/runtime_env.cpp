#include "pendulum/runtime_env.h"

#include <cstdlib>
#include <string>

#include <Python.h>

void configure_pythonpath_from_venv() {
    const char* venv = std::getenv("VIRTUAL_ENV");
    if (venv == nullptr || *venv == '\0') {
        return;
    }

    std::string site_packages = std::string(venv) + "/lib/python" +
                                std::to_string(PY_MAJOR_VERSION) + "." +
                                std::to_string(PY_MINOR_VERSION) + "/site-packages";

    const char* old_pythonpath = std::getenv("PYTHONPATH");
    std::string new_pythonpath = site_packages;
    if (old_pythonpath != nullptr && *old_pythonpath != '\0') {
        new_pythonpath += ":" + std::string(old_pythonpath);
    }

    setenv("PYTHONPATH", new_pythonpath.c_str(), 1);
}
