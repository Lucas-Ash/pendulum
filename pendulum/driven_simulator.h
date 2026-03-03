#pragma once

#include "pendulum/driven_config.h"
#include "pendulum/driven_result.h"

class DrivenPendulumSimulator {
public:
    explicit DrivenPendulumSimulator(const DrivenConfig& config);
    DrivenSimulationResult simulate() const;

private:
    DrivenConfig config_;
};
