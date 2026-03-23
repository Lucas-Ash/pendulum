#!/bin/bash
g++ -O3 -std=c++17 -I. \
  -I"$(python3 -c 'import numpy; print(numpy.get_include())')" \
  integrated_pendulum_modular.cpp \
  modules/core/config_utils.cpp \
  modules/core/error_analysis.cpp \
  modules/core/plotting_utils.cpp \
  modules/config/experiment_config.cpp \
  modules/config/damped_config.cpp \
  modules/config/driven_config.cpp \
  modules/config/coupled_config.cpp \
  modules/integrators/jacobi_elliptic.cpp \
  modules/simple/pendulum_simulator.cpp \
  modules/simple/plotting.cpp \
  modules/simple/reporting.cpp \
  modules/simple/simple_io.cpp \
  modules/damped/damped_io.cpp \
  modules/damped/damped_plotting.cpp \
  modules/damped/damped_reporting.cpp \
  modules/damped/damped_simulator.cpp \
  modules/driven/driven_io.cpp \
  modules/driven/driven_plotting.cpp \
  modules/driven/driven_reporting.cpp \
  modules/driven/driven_simulator.cpp \
  modules/coupled/coupled_simulator.cpp \
  modules/runtime/runtime_env.cpp \
  -o integrated_pendulum \
  $(python3-config --cflags --embed) \
  $(python3-config --ldflags --embed)
