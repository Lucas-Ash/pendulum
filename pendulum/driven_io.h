#pragma once

#include "pendulum/driven_config.h"
#include "pendulum/driven_result.h"
#include <string>

void write_driven_data_file(const std::string& path, const DrivenSimulationResult& result);
