#include "modules/damped_simulator.h"

#include <cmath>
#include <stdexcept>

namespace {

using DerivFunc = void (*)(const double[2], double[2], double, double);

void derivs_nonlinear(const double y[2], double dydt[2], double gamma, double omega0) {
    dydt[0] = y[1];
    dydt[1] = -2.0 * gamma * y[1] - omega0 * omega0 * std::sin(y[0]);
}

void derivs_linear(const double y[2], double dydt[2], double gamma, double omega0) {
    dydt[0] = y[1];
    dydt[1] = -2.0 * gamma * y[1] - omega0 * omega0 * y[0];
}

void rk4_step(DerivFunc f, double y[2], double dt, double gamma, double omega0) {
    double k1[2], k2[2], k3[2], k4[2], ytmp[2];

    f(y, k1, gamma, omega0);
    for (int i = 0; i < 2; ++i) ytmp[i] = y[i] + 0.5 * dt * k1[i];

    f(ytmp, k2, gamma, omega0);
    for (int i = 0; i < 2; ++i) ytmp[i] = y[i] + 0.5 * dt * k2[i];

    f(ytmp, k3, gamma, omega0);
    for (int i = 0; i < 2; ++i) ytmp[i] = y[i] + dt * k3[i];

    f(ytmp, k4, gamma, omega0);
    for (int i = 0; i < 2; ++i) {
        y[i] += (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
}

}  // namespace

DampedPendulumSimulator::DampedPendulumSimulator(const DampedConfig& config)
    : config_(config) {}

void DampedPendulumSimulator::analytical_solution(double t, double omega_d, double gamma,
                                                  double theta0, double theta_dot0,
                                                  double& theta_exact, double& omega_exact) {
    const double A = theta0;
    const double B = (theta_dot0 + gamma * theta0) / omega_d;
    const double exp_term = std::exp(-gamma * t);
    const double cos_term = std::cos(omega_d * t);
    const double sin_term = std::sin(omega_d * t);

    theta_exact = exp_term * (A * cos_term + B * sin_term);
    omega_exact = -gamma * theta_exact + exp_term * (-A * omega_d * sin_term + B * omega_d * cos_term);
}

double DampedPendulumSimulator::total_energy(double theta, double theta_dot, double omega0) {
    return 0.5 * theta_dot * theta_dot + omega0 * omega0 * (1.0 - std::cos(theta));
}

SimulationResult DampedPendulumSimulator::simulate() const {
    const auto& p = config_.physical;
    const auto& s = config_.simulation;

    SimulationResult result;
    double omega0 = std::sqrt(p.g / p.L);
    if (p.gamma >= omega0) {
        throw std::runtime_error("This analytical model assumes underdamped motion: gamma < omega0");
    }
    double omega_d = std::sqrt(omega0 * omega0 - p.gamma * p.gamma);

    const int nsteps = static_cast<int>((s.t_end - s.t_start) / s.dt + 0.5);
    result.rk4_steps = nsteps;

    double y_nl[2] = {p.theta0, p.theta_dot0};

    double max_theta_error = 0.0;
    double max_omega_error = 0.0;
    double avg_theta_error = 0.0;
    double avg_omega_error = 0.0;
    double max_theta_rel_error = 0.0;
    double max_omega_rel_error = 0.0;
    double avg_theta_rel_error = 0.0;
    double avg_omega_rel_error = 0.0;

    for (int n = 0; n <= nsteps; ++n) {
        const double t = s.t_start + n * s.dt;
        double theta_exact = 0.0;
        double omega_exact = 0.0;
        analytical_solution(t, omega_d, p.gamma, p.theta0, p.theta_dot0, theta_exact, omega_exact);
        
        result.t.push_back(t);
        result.theta.push_back(y_nl[0]);
        result.omega.push_back(y_nl[1]);
        result.theta_analytical.push_back(theta_exact);
        result.omega_analytical.push_back(omega_exact);
        result.energy.push_back(total_energy(y_nl[0], y_nl[1], omega0));

        const double theta_err = std::abs(y_nl[0] - theta_exact);
        const double omega_err = std::abs(y_nl[1] - omega_exact);
        const double theta_rel = theta_err / std::max(std::abs(theta_exact), 1e-12);
        const double omega_rel = omega_err / std::max(std::abs(omega_exact), 1e-12);

        result.theta_errors.push_back(theta_err);
        result.omega_errors.push_back(omega_err);

        max_theta_error = std::max(max_theta_error, theta_err);
        max_omega_error = std::max(max_omega_error, omega_err);
        avg_theta_error += theta_err;
        avg_omega_error += omega_err;
        max_theta_rel_error = std::max(max_theta_rel_error, theta_rel);
        max_omega_rel_error = std::max(max_omega_rel_error, omega_rel);
        avg_theta_rel_error += theta_rel;
        avg_omega_rel_error += omega_rel;

        if (n < nsteps) {
            rk4_step(derivs_nonlinear, y_nl, s.dt, p.gamma, omega0);
        }
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
