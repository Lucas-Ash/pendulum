#pragma once

#include "modules/config/coupled_config.h"
#include "modules/coupled/coupled_state.h"
#include <string>

class CoupledSimulator {
public:
    explicit CoupledSimulator(const CoupledConfig& config);
    CoupledSimulationResult simulate() const;

private:
    CoupledConfig config_;
};

void write_coupled_data_file(const std::string& path, const CoupledSimulationResult& result);
void write_coupled_plot_script(const std::string& script_path, const std::string& data_file, const std::string& output_png);
