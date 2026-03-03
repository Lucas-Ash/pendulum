#!/bin/bash
g++ -O3 -std=c++17 -I. \
  -I"$(python3 -c 'import numpy; print(numpy.get_include())')" \
  integrated_pendulum_modular.cpp \
  modules/experiment_config.cpp \
  modules/pendulum_simulator.cpp \
  modules/plotting.cpp \
  modules/reporting.cpp \
  modules/runtime_env.cpp \
  modules/damped_config.cpp \
  modules/damped_io.cpp \
  modules/damped_plotting.cpp \
  modules/damped_reporting.cpp \
  modules/damped_simulator.cpp \
  modules/driven_config.cpp \
  modules/driven_io.cpp \
  modules/driven_plotting.cpp \
  modules/driven_reporting.cpp \
  modules/driven_simulator.cpp \
  -o integrated_pendulum \
  $(python3-config --cflags --embed) \
  $(python3-config --ldflags --embed)
