#pragma once

#include "modules/config/experiment_config.h"
#include "modules/core/simulation_result.h"

void print_accuracy_report(const ExperimentConfig& config,
                           const SimulationResult& result);
