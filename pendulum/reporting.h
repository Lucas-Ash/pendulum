#pragma once

#include "pendulum/experiment_config.h"
#include "pendulum/simulation_result.h"

void print_accuracy_report(const ExperimentConfig& config,
                           const SimulationResult& result);
