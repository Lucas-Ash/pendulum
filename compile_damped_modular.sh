g++ -O2 -std=c++17 -I. \
  damped_pendulum_modular.cpp \
  pendulum/damped_config.cpp \
  pendulum/damped_io.cpp \
  pendulum/damped_plotting.cpp \
  pendulum/damped_reporting.cpp \
  pendulum/damped_simulator.cpp \
  -o damped_pendulum_modular
