#pragma once

#include <string>

#include "modules/damped_result.h"

void write_damped_data_file(const std::string& path, const DampedSimulationResult& result);
