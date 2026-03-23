#pragma once

#include "modules/config/driven_config.h"
#include "modules/core/simulation_result.h"
#include "modules/driven/driven_sweep_result.h"
#include <string>

void write_driven_data_file(const std::string& path, const SimulationResult& result);
void write_driven_sweep_data_file(const std::string& path, const DrivenSweepResult& result);
