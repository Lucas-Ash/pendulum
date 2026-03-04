#pragma once

#include "modules/driven_config.h"
#include "modules/simulation_result.h"
#include <string>

void write_driven_data_file(const std::string& path, const SimulationResult& result);
