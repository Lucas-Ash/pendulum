#pragma once

#include "modules/driven_config.h"
#include "modules/driven_sweep_result.h"
#include "modules/simulation_result.h"

void print_driven_simulation_summary(const DrivenConfig& config, const SimulationResult& result);
void print_driven_sweep_summary(const DrivenConfig& config, const DrivenSweepResult& result);
