#pragma once

#include "modules/damped_config.h"
#include "modules/simulation_result.h"

void print_damped_simulation_summary(const DampedConfig& config,
                                     const SimulationResult& result);
