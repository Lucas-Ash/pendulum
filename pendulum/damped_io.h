#pragma once

#include <string>

#include "pendulum/damped_result.h"

void write_damped_data_file(const std::string& path, const DampedSimulationResult& result);
