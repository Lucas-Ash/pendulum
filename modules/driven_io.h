#pragma once

#include "modules/driven_config.h"
#include "modules/driven_result.h"
#include <string>

void write_driven_data_file(const std::string& path, const DrivenSimulationResult& result);
