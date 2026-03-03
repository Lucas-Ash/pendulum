#pragma once

#include "pendulum/damped_config.h"
#include "pendulum/damped_result.h"

void print_damped_simulation_summary(const DampedConfig& config,
                                     const DampedSimulationResult& result);
