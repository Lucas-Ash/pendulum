#!/bin/bash
g++ -O2 -std=c++17 -I. \
  driven_pendulum_modular.cpp \
  pendulum/driven_config.cpp \
  pendulum/driven_io.cpp \
  pendulum/driven_plotting.cpp \
  pendulum/driven_reporting.cpp \
  pendulum/driven_simulator.cpp \
  -o driven_pendulum_modular
