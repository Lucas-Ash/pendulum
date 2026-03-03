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

double DampedPendulumSimulator::analytical_solution(double t, double omega_d, double gamma,
                                                    double theta0, double theta_dot0) {
    const double A = theta0;
    const double B = (theta_dot0 + gamma * theta0) / omega_d;
    return std::exp(-gamma * t) *
           (A * std::cos(omega_d * t) + B * std::sin(omega_d * t));
}

double DampedPendulumSimulator::total_energy(double theta, double theta_dot, double omega0) {
    return 0.5 * theta_dot * theta_dot + omega0 * omega0 * (1.0 - std::cos(theta));
}

DampedSimulationResult DampedPendulumSimulator::simulate() const {
    const auto& p = config_.physical;
    const auto& s = config_.simulation;

    DampedSimulationResult result;
    result.omega0 = std::sqrt(p.g / p.L);
    if (p.gamma >= result.omega0) {
        throw std::runtime_error("This analytical model assumes underdamped motion: gamma < omega0");
    }
    result.omega_d = std::sqrt(result.omega0 * result.omega0 - p.gamma * p.gamma);
    result.zeta = p.gamma / result.omega0;

    const int nsteps = static_cast<int>((s.t_end - s.t_start) / s.dt + 0.5);
    result.rk4_steps = nsteps;

    double y_nl[2] = {p.theta0, p.theta_dot0};
    double y_lin[2] = {p.theta0, p.theta_dot0};

    for (int n = 0; n <= nsteps; ++n) {
        const double t = s.t_start + n * s.dt;
        const double theta_exact = analytical_solution(
            t, result.omega_d, p.gamma, p.theta0, p.theta_dot0);
        const double err_lin = std::fabs(theta_exact - y_lin[0]);
        const double err_nl = std::fabs(theta_exact - y_nl[0]);
        const double E_nl = total_energy(y_nl[0], y_nl[1], result.omega0);

        if (err_lin > result.max_err_linear) result.max_err_linear = err_lin;
        if (err_nl > result.max_err_nonlinear) result.max_err_nonlinear = err_nl;

        if (n % s.output_every == 0) {
            result.samples.push_back(DampedSample{
                t,
                theta_exact,
                y_lin[0],
                y_nl[0],
                err_lin,
                err_nl,
                E_nl
            });
        }

        if (n < nsteps) {
            rk4_step(derivs_nonlinear, y_nl, s.dt, p.gamma, result.omega0);
            rk4_step(derivs_linear, y_lin, s.dt, p.gamma, result.omega0);
        }
    }

    return result;
}
