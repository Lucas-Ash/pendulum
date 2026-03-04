#pragma once

#include "modules/experiment_config.h"
#include "modules/simulation_result.h"
#include <string>

void write_simple_data_file(const std::string& path, const SimulationResult& result);
