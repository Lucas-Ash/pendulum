#pragma once

#include "modules/driven_config.h"
#include "modules/simulation_result.h"

class DrivenPendulumSimulator {
public:
    explicit DrivenPendulumSimulator(const DrivenConfig& config);
    SimulationResult simulate() const;

private:
    DrivenConfig config_;
};
