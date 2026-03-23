#pragma once

#include "modules/config/damped_config.h"
#include "modules/core/simulation_result.h"

void print_damped_simulation_summary(const DampedConfig& config,
                                     const SimulationResult& result);
