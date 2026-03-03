g++ -O3 -std=c++17 -I. \
  -I"$(python3 -c 'import numpy; print(numpy.get_include())')" \
  simple_pendulum_modular.cpp \
  pendulum/experiment_config.cpp \
  pendulum/pendulum_simulator.cpp \
  pendulum/reporting.cpp \
  pendulum/plotting.cpp \
  pendulum/runtime_env.cpp \
  -o execute_modular \
  $(python3-config --cflags --ldflags --embed)
