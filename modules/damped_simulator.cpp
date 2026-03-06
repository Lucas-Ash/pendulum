#include "modules/damped_simulator.h"
#include "modules/rk_integrators.h"
#include "modules/error_analysis.h"

#include <cmath>
#include <stdexcept>

namespace {

template <typename DerivFunc>
integrator::State advance_state(const std::string& method,
                                double t,
                                double dt,
                                const integrator::State& state,
                                DerivFunc derivs) {
    if (method == "rk3") return integrator::rk3_step(t, dt, state, derivs);
    if (method == "rk5") return integrator::rk5_step(t, dt, state, derivs);
    if (method == "semi_implicit_euler") {
        return integrator::semi_implicit_euler_step(t, dt, state, derivs);
    }
    if (method == "leapfrog") return integrator::leapfrog_step(t, dt, state, derivs);
    if (method == "ruth4") return integrator::ruth4_step(t, dt, state, derivs);
    if (method == "rk23") return integrator::rk23_step(t, dt, state, derivs);
    if (method == "rkf45") return integrator::rkf45_step(t, dt, state, derivs);
    return integrator::rk4_step(t, dt, state, derivs);
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

double DampedPendulumSimulator::total_energy(double theta, double theta_dot,
                                             double omega0_sq,
                                             const restoring_force::Config& restoring) {
    return 0.5 * theta_dot * theta_dot
         + omega0_sq * restoring_force::potential(theta, restoring);
}

SimulationResult DampedPendulumSimulator::simulate() const {
    const auto& p = config_.physical;
    const auto& s = config_.simulation;
    const auto& settings = config_.settings;

    SimulationResult result;
    const double omega0_sq = p.g / p.L;
    if (settings.error_reference_factor <= 0) {
        throw std::runtime_error("error_reference_factor must be > 0");
    }
    if (settings.error_mode == error_reference::Mode::HdReference &&
        settings.error_reference_factor < 2) {
        throw std::runtime_error(
            "hd_reference mode requires error_reference_factor >= 2");
    }

    const bool use_analytical = settings.error_mode == error_reference::Mode::Analytical;
    double omega_d = 0.0;
    if (use_analytical) {
        const double linearized_slope = restoring_force::linearized_slope(p.restoring_force);
        const double omega_lin_sq = omega0_sq * linearized_slope;
        if (omega_lin_sq <= 0.0) {
            throw std::runtime_error(
                "Damped analytical reference requires positive linearized stiffness.");
        }

        const double omega_lin = std::sqrt(omega_lin_sq);
        if (p.gamma >= omega_lin) {
            throw std::runtime_error(
                "This analytical model assumes underdamped motion: gamma < omega_linearized");
        }
        omega_d = std::sqrt(omega_lin_sq - p.gamma * p.gamma);
    }

    const int nsteps = static_cast<int>((s.t_end - s.t_start) / s.dt + 0.5);
    const double dt_ref = s.dt / static_cast<double>(settings.error_reference_factor);
    double t_ref = s.t_start;
    result.rk4_steps = nsteps;

    auto derivs = [gamma = p.gamma, omega0_sq, restoring = p.restoring_force](double t_val, const integrator::State& state) -> integrator::State {
        return {state.omega, -2.0 * gamma * state.omega
                              - omega0_sq * restoring_force::term(state.theta, restoring)};
    };

    integrator::State state = {p.theta0, p.theta_dot0};
    integrator::State state_ref = state;

    for (int n = 0; n <= nsteps; ++n) {
        const double t = s.t_start + n * s.dt;
        double theta_exact = 0.0;
        double omega_exact = 0.0;
        if (settings.error_mode == error_reference::Mode::None) {
            theta_exact = state.theta;
            omega_exact = state.omega;
        } else if (settings.error_mode == error_reference::Mode::HdReference) {
            theta_exact = state_ref.theta;
            omega_exact = state_ref.omega;
        } else {
            analytical_solution(t, omega_d, p.gamma, p.theta0, p.theta_dot0,
                                theta_exact, omega_exact);
        }
        
        result.t.push_back(t);
        result.theta.push_back(state.theta);
        result.omega.push_back(state.omega);
        result.theta_analytical.push_back(theta_exact);
        result.omega_analytical.push_back(omega_exact);
        result.energy.push_back(
            total_energy(state.theta, state.omega, omega0_sq, p.restoring_force));

        if (n < nsteps) {
            state = advance_state(settings.integrator, t, s.dt, state, derivs);
            if (settings.error_mode == error_reference::Mode::HdReference) {
                for (int k = 0; k < settings.error_reference_factor; ++k) {
                    state_ref = advance_state(settings.integrator, t_ref, dt_ref,
                                              state_ref, derivs);
                    t_ref += dt_ref;
                }
            }
        }
    }

    compute_error_statistics(result);

    return result;
}
