#pragma once

#include "modules/config/driven_config.h"
#include "modules/core/simulation_result.h"
#include "modules/driven/driven_sweep_result.h"

class DrivenPendulumSimulator {
public:
    explicit DrivenPendulumSimulator(const DrivenConfig& config);
    SimulationResult simulate() const;
    DrivenSweepResult simulate_sweep() const;

private:
    DrivenConfig config_;
};
