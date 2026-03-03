#pragma once

#include "modules/experiment_config.h"
#include "modules/simulation_result.h"

void print_accuracy_report(const ExperimentConfig& config,
                           const SimulationResult& result);
