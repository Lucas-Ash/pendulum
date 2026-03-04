#pragma once

#include <string>

#include "modules/simulation_result.h"

void write_damped_data_file(const std::string& path, const SimulationResult& result);
