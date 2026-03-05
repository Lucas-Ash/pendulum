#include "modules/pendulum_simulator.h"
#include "modules/rk_integrators.h"
#include "modules/error_analysis.h"

#include <algorithm>
#include <cmath>

PendulumSimulator::PendulumSimulator(double length, double gravity,
                                     double timestep, double tmax)
    : g_(gravity), L_(length), dt_(timestep), t_max_(tmax) {}

double PendulumSimulator::gravity() const { return g_; }
double PendulumSimulator::length() const { return L_; }
double PendulumSimulator::dt() const { return dt_; }
double PendulumSimulator::t_max() const { return t_max_; }



void PendulumSimulator::exact_linear_state(double t, double theta0, double omega0,
                                           double& theta_exact, double& omega_exact) const {
    const double wn = std::sqrt(g_ / L_);
    const double c = std::cos(wn * t);
    const double s = std::sin(wn * t);
    theta_exact = theta0 * c + (omega0 / wn) * s;
    omega_exact = -theta0 * wn * s + omega0 * c;
}

SimulationResult PendulumSimulator::simulate(double theta0, double omega0, const std::string& integrator) const {
    SimulationResult result;

    double t = 0.0;
    double theta = theta0;
    double omega = omega0;

    while (t <= t_max_ + 1e-12) {
        result.t.push_back(t);
        result.theta.push_back(theta);
        result.omega.push_back(omega);
        
        double theta_ref = 0.0;
        double omega_ref = 0.0;
        exact_linear_state(t, theta0, omega0, theta_ref, omega_ref);
        result.theta_analytical.push_back(theta_ref);
        result.omega_analytical.push_back(omega_ref);
        
        double energy = 0.5 * omega * omega + (g_ / L_) * (1.0 - std::cos(theta));
        result.energy.push_back(energy);

        auto derivs = [this](double t_val, const integrator::State& s) -> integrator::State {
            return {s.omega, -(g_ / L_) * std::sin(s.theta)};
        };

        integrator::State state = {theta, omega};
        if (integrator == "rk3") {
            state = integrator::rk3_step(t, dt_, state, derivs);
        } else if (integrator == "rk5") {
            state = integrator::rk5_step(t, dt_, state, derivs);
        } else {
            state = integrator::rk4_step(t, dt_, state, derivs);
        }
        theta = state.theta;
        omega = state.omega;
        t += dt_;
    }

    compute_error_statistics(result);
    result.rk4_steps = result.t.size() - 1;

    return result;
}
