#include "modules/damped_simulator.h"
#include "modules/rk_integrators.h"
#include "modules/error_analysis.h"

#include <cmath>
#include <stdexcept>

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

        if (n < nsteps) {
            auto derivs = [gamma = p.gamma, omega0](double t_val, const integrator::State& s) -> integrator::State {
                return {s.omega, -2.0 * gamma * s.omega - omega0 * omega0 * std::sin(s.theta)};
            };
            integrator::State state = {y_nl[0], y_nl[1]};
            
            if (config_.settings.integrator == "rk3") {
                state = integrator::rk3_step(t, s.dt, state, derivs);
            } else if (config_.settings.integrator == "rk5") {
                state = integrator::rk5_step(t, s.dt, state, derivs);
            } else if (config_.settings.integrator == "semi_implicit_euler") {
                state = integrator::semi_implicit_euler_step(t, s.dt, state, derivs);
            } else if (config_.settings.integrator == "leapfrog") {
                state = integrator::leapfrog_step(t, s.dt, state, derivs);
            } else if (config_.settings.integrator == "ruth4") {
                state = integrator::ruth4_step(t, s.dt, state, derivs);
            } else if (config_.settings.integrator == "rk23") {
                state = integrator::rk23_step(t, s.dt, state, derivs);
            } else if (config_.settings.integrator == "rkf45") {
                state = integrator::rkf45_step(t, s.dt, state, derivs);
            } else {
                state = integrator::rk4_step(t, s.dt, state, derivs);
            }
            
            y_nl[0] = state.theta;
            y_nl[1] = state.omega;
        }
    }

    compute_error_statistics(result);

    return result;
}
