#!/bin/bash
g++ -O3 -std=c++17 -I. \
  -I"$(python3 -c 'import numpy; print(numpy.get_include())')" \
  integrated_pendulum_modular.cpp \
  pendulum/experiment_config.cpp \
  pendulum/pendulum_simulator.cpp \
  pendulum/plotting.cpp \
  pendulum/reporting.cpp \
  pendulum/runtime_env.cpp \
  pendulum/damped_config.cpp \
  pendulum/damped_io.cpp \
  pendulum/damped_plotting.cpp \
  pendulum/damped_reporting.cpp \
  pendulum/damped_simulator.cpp \
  pendulum/driven_config.cpp \
  pendulum/driven_io.cpp \
  pendulum/driven_plotting.cpp \
  pendulum/driven_reporting.cpp \
  pendulum/driven_simulator.cpp \
  -o integrated_pendulum \
  $(python3-config --cflags --embed) \
  $(python3-config --ldflags --embed)
