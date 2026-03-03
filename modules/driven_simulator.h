#pragma once

#include "modules/driven_config.h"
#include "modules/driven_result.h"

class DrivenPendulumSimulator {
public:
    explicit DrivenPendulumSimulator(const DrivenConfig& config);
    DrivenSimulationResult simulate() const;

private:
    DrivenConfig config_;
};
