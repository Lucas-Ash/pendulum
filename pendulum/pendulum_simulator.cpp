#include "pendulum/pendulum_simulator.h"

#include <algorithm>
#include <cmath>

PendulumSimulator::PendulumSimulator(double length, double gravity,
                                     double timestep, double tmax)
    : g_(gravity), L_(length), dt_(timestep), t_max_(tmax) {}

double PendulumSimulator::gravity() const { return g_; }
double PendulumSimulator::length() const { return L_; }
double PendulumSimulator::dt() const { return dt_; }
double PendulumSimulator::t_max() const { return t_max_; }

void PendulumSimulator::derivatives(double theta, double omega,
                                    double& dtheta_dt, double& domega_dt) const {
    dtheta_dt = omega;
    domega_dt = -(g_ / L_) * std::sin(theta);
}

void PendulumSimulator::rk4_step(double& theta, double& omega) const {
    double k1_theta, k1_omega, k2_theta, k2_omega;
    double k3_theta, k3_omega, k4_theta, k4_omega;
    double theta_mid, omega_mid;

    derivatives(theta, omega, k1_theta, k1_omega);

    theta_mid = theta + 0.5 * dt_ * k1_theta;
    omega_mid = omega + 0.5 * dt_ * k1_omega;
    derivatives(theta_mid, omega_mid, k2_theta, k2_omega);

    theta_mid = theta + 0.5 * dt_ * k2_theta;
    omega_mid = omega + 0.5 * dt_ * k2_omega;
    derivatives(theta_mid, omega_mid, k3_theta, k3_omega);

    theta_mid = theta + dt_ * k3_theta;
    omega_mid = omega + dt_ * k3_omega;
    derivatives(theta_mid, omega_mid, k4_theta, k4_omega);

    theta += (dt_ / 6.0) * (k1_theta + 2.0 * k2_theta + 2.0 * k3_theta + k4_theta);
    omega += (dt_ / 6.0) * (k1_omega + 2.0 * k2_omega + 2.0 * k3_omega + k4_omega);
}

void PendulumSimulator::exact_linear_state(double t, double theta0, double omega0,
                                           double& theta_exact, double& omega_exact) const {
    const double wn = std::sqrt(g_ / L_);
    const double c = std::cos(wn * t);
    const double s = std::sin(wn * t);
    theta_exact = theta0 * c + (omega0 / wn) * s;
    omega_exact = -theta0 * wn * s + omega0 * c;
}

SimulationResult PendulumSimulator::simulate(double theta0, double omega0) const {
    SimulationResult result;

    double t = 0.0;
    double theta = theta0;
    double omega = omega0;

    while (t <= t_max_ + 1e-12) {
        result.t.push_back(t);
        result.theta.push_back(theta);
        result.omega.push_back(omega);

        rk4_step(theta, omega);
        t += dt_;
    }

    result.theta_errors.resize(result.t.size());
    result.omega_errors.resize(result.t.size());

    double max_theta_error = 0.0;
    double max_omega_error = 0.0;
    double avg_theta_error = 0.0;
    double avg_omega_error = 0.0;
    double max_theta_rel_error = 0.0;
    double max_omega_rel_error = 0.0;
    double avg_theta_rel_error = 0.0;
    double avg_omega_rel_error = 0.0;

    for (size_t i = 0; i < result.t.size(); ++i) {
        double theta_ref = 0.0;
        double omega_ref = 0.0;
        exact_linear_state(result.t[i], theta0, omega0, theta_ref, omega_ref);

        const double theta_err = std::abs(result.theta[i] - theta_ref);
        const double omega_err = std::abs(result.omega[i] - omega_ref);
        const double theta_rel = theta_err / std::max(std::abs(theta_ref), 1e-12);
        const double omega_rel = omega_err / std::max(std::abs(omega_ref), 1e-12);

        result.theta_errors[i] = theta_err;
        result.omega_errors[i] = omega_err;

        max_theta_error = std::max(max_theta_error, theta_err);
        max_omega_error = std::max(max_omega_error, omega_err);
        avg_theta_error += theta_err;
        avg_omega_error += omega_err;
        max_theta_rel_error = std::max(max_theta_rel_error, theta_rel);
        max_omega_rel_error = std::max(max_omega_rel_error, omega_rel);
        avg_theta_rel_error += theta_rel;
        avg_omega_rel_error += omega_rel;
    }

    const double sample_count = static_cast<double>(result.t.size());
    result.theta_stats = {
        max_theta_error,
        avg_theta_error / sample_count,
        max_theta_rel_error,
        avg_theta_rel_error / sample_count,
    };
    result.omega_stats = {
        max_omega_error,
        avg_omega_error / sample_count,
        max_omega_rel_error,
        avg_omega_rel_error / sample_count,
    };

    return result;
}
