#pragma once

#include "pendulum/simulation_result.h"

class PendulumSimulator {
public:
    PendulumSimulator(double length, double gravity = 9.81,
                      double timestep = 0.01, double tmax = 20.0);

    double gravity() const;
    double length() const;
    double dt() const;
    double t_max() const;

    SimulationResult simulate(double theta0 = 0.5, double omega0 = 0.0) const;

private:
    double g_;
    double L_;
    double dt_;
    double t_max_;

    void derivatives(double theta, double omega,
                     double& dtheta_dt, double& domega_dt) const;
    void rk4_step(double& theta, double& omega) const;
    void exact_linear_state(double t, double theta0, double omega0,
                            double& theta_exact, double& omega_exact) const;
};
