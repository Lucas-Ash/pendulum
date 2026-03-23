#pragma once

#include "modules/config/driven_config.h"
#include "modules/core/simulation_result.h"
#include "modules/driven/driven_sweep_result.h"

void print_driven_simulation_summary(const DrivenConfig& config, const SimulationResult& result);
void print_driven_sweep_summary(const DrivenConfig& config, const DrivenSweepResult& result);
