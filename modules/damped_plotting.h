#pragma once

#include "modules/damped_config.h"
#include "modules/simulation_result.h"

void render_damped_plots(const DampedConfig& config,
                         const SimulationResult& result);
