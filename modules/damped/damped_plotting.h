#pragma once

#include "modules/config/damped_config.h"
#include "modules/core/simulation_result.h"

void render_damped_plots(const DampedConfig& config,
                         const SimulationResult& result);
