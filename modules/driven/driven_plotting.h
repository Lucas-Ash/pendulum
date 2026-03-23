#pragma once

#include "modules/config/driven_config.h"
#include "modules/core/simulation_result.h"
#include "modules/driven/driven_sweep_result.h"

void render_driven_plots(const DrivenConfig& config, const SimulationResult& result);
void render_driven_sweep_plots(const DrivenConfig& config, const DrivenSweepResult& result);
