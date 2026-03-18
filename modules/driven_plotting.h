#pragma once

#include "modules/driven_config.h"
#include "modules/driven_sweep_result.h"
#include "modules/simulation_result.h"

void render_driven_plots(const DrivenConfig& config, const SimulationResult& result);
void render_driven_sweep_plots(const DrivenConfig& config, const DrivenSweepResult& result);
