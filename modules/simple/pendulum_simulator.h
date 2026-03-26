#pragma once

#include <string>
#include "modules/core/additional_terms.h"
#include "modules/core/error_reference.h"
#include "modules/core/restoring_force.h"
#include "modules/core/simulation_result.h"

class PendulumSimulator {
public:
    PendulumSimulator(double length, double gravity = 9.81,
                      double timestep = 0.01, double tmax = 20.0,
                      restoring_force::Config restoring_force = {},
                      additional_terms::Config additional_terms = {});

    double gravity() const;
    double length() const;
    double dt() const;
    double t_max() const;

    SimulationResult simulate(double theta0 = 0.5, double omega0 = 0.0,
                              const std::string& integrator = "rk4",
                              const std::string& analytical_model = "linear",
                              error_reference::Mode error_mode = error_reference::Mode::Analytical,
                              int error_reference_factor = 50) const;

private:
    double g_;
    double L_;
    double dt_;
    double t_max_;
    restoring_force::Config restoring_force_;
    additional_terms::Config additional_terms_;


    void exact_linear_state(double t, double theta0, double omega0,
                            double& theta_exact, double& omega_exact) const;
    void exact_nonlinear_state(double t, double theta0, double omega0,
                               double& theta_exact, double& omega_exact) const;
    void exact_duffing_state(double t, double theta0, double omega0,
                             double& theta_exact, double& omega_exact) const;
    void exact_ermakov_pinney_state(double t, double theta0, double omega0,
                                    double& theta_exact, double& omega_exact) const;
    void exact_toda_state(double t, double theta0, double omega0,
                          double& theta_exact, double& omega_exact) const;
    void exact_toda_soliton_state(double t, double theta0, double omega0,
                                  double& theta_exact, double& omega_exact) const;
    static void exact_harmonic_state(double t, double theta0, double omega0,
                                     double omega_natural,
                                     double& theta_exact, double& omega_exact);
};
