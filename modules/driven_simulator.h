#pragma once

#include "modules/driven_config.h"
#include "modules/driven_sweep_result.h"
#include "modules/simulation_result.h"

class DrivenPendulumSimulator {
public:
    explicit DrivenPendulumSimulator(const DrivenConfig& config);
    SimulationResult simulate() const;
    DrivenSweepResult simulate_sweep() const;

private:
    DrivenConfig config_;
};
