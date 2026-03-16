#include "modules/driven_simulator.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <functional>

#include "modules/rk_integrators.h"
#include "modules/error_analysis.h"

namespace {

    template <typename DerivFunc>
    integrator::State advance_state(const std::string& method,
                                    double t,
                                    double dt,
                                    const integrator::State& state,
                                    DerivFunc derivs,
                                    double den_gamma,
                                    double den_omega0,
                                    const std::function<double(double, const integrator::State&)>& den_residual) {
        if (method == "rk3") return integrator::rk3_step(t, dt, state, derivs);
        if (method == "rk5") return integrator::rk5_step(t, dt, state, derivs);
        if (method == "semi_implicit_euler") {
            return integrator::semi_implicit_euler_step(t, dt, state, derivs);
        }
        if (method == "leapfrog") return integrator::leapfrog_step(t, dt, state, derivs);
        if (method == "ruth4") return integrator::ruth4_step(t, dt, state, derivs);
        if (method == "rk23") return integrator::rk23_step(t, dt, state, derivs);
        if (method == "rkf45") return integrator::rkf45_step(t, dt, state, derivs);
        if (integrator::is_den_method(method)) {
            return integrator::den3_step(t, dt, state, den_gamma, den_omega0, den_residual);
        }
        return integrator::rk4_step(t, dt, state, derivs);
    }

    double damping_term(double theta, double omega, const DrivenPhysicalConfig& p) {
        return damping_force::term(
            theta, omega, p.damping_model, p.damping, p.damping_cubic);
    }

    void validate_position_only_integrator(const std::string& method,
                                           const DrivenPhysicalConfig& p) {
        if (!integrator::is_position_only_integrator(method)) {
            return;
        }
        if (std::abs(p.damping) > 1e-12 || std::abs(p.damping_cubic) > 1e-12) {
            throw std::runtime_error(
                method + " requires theta'' = a(t, theta); it is not valid with velocity-dependent damping.");
        }
    }

    void analytical_solution(double t, const DrivenPhysicalConfig& p,
                             double omega_linearized_sq,
                             double& theta_exact, double& omega_exact) {
        double denominator = std::sqrt(std::pow(omega_linearized_sq - p.omega_drive * p.omega_drive, 2)
                                     + std::pow(p.damping * p.omega_drive, 2));
        double R = p.A / denominator;
        double delta = std::atan2(p.damping * p.omega_drive,
                                  omega_linearized_sq - p.omega_drive * p.omega_drive);
        
        double theta_steady = R * std::cos(p.omega_drive * t - delta);
        double omega_steady = -R * p.omega_drive * std::sin(p.omega_drive * t - delta);
        
        double damping_half = p.damping / 2.0;
        double omega_trans_sq = omega_linearized_sq - damping_half * damping_half;
        if (omega_trans_sq <= 0) {
            theta_exact = theta_steady;
            omega_exact = omega_steady;
            return;
        }
        double omega_trans = std::sqrt(omega_trans_sq);
        
        double cos_delta = std::cos(delta);
        double sin_delta = std::sin(delta);
        
        double C1 = p.theta0 - R * cos_delta;
        double C2 = (p.omega0 + damping_half * C1 - R * p.omega_drive * sin_delta) / omega_trans;
        
        double exp_term = std::exp(-damping_half * t);
        double cos_term = std::cos(omega_trans * t);
        double sin_term = std::sin(omega_trans * t);
        
        double theta_transient = exp_term * (C1 * cos_term + C2 * sin_term);
        double omega_transient = -damping_half * theta_transient + exp_term * (-C1 * omega_trans * sin_term + C2 * omega_trans * cos_term);
                                 
        theta_exact = theta_steady + theta_transient;
        omega_exact = omega_steady + omega_transient;
    }

    double total_energy(double theta, double theta_dot, double omega0_sq,
                        const restoring_force::Config& restoring) {
        return 0.5 * theta_dot * theta_dot
             + omega0_sq * restoring_force::potential(theta, restoring);
    }
}

DrivenPendulumSimulator::DrivenPendulumSimulator(const DrivenConfig& config)
    : config_(config) {}

SimulationResult DrivenPendulumSimulator::simulate() const {
    const auto& p = config_.physical;
    const auto& s = config_.simulation;
    const auto& settings = config_.settings;
    validate_position_only_integrator(settings.integrator, p);

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
    double omega_lin_sq = 0.0;
    if (use_analytical) {
        const double linearized_slope = restoring_force::linearized_slope(p.restoring_force);
        omega_lin_sq = omega0_sq * linearized_slope;
        if (omega_lin_sq <= 0.0) {
            throw std::runtime_error(
                "Driven analytical reference requires positive linearized stiffness.");
        }
        if (std::abs(p.damping_cubic) > 1e-15) {
            throw std::runtime_error(
                "Driven analytical reference does not support damping_cubic != 0. "
                "Use error_analysis: hd_reference or none for polynomial damping.");
        }
    }

    integrator::State state = {p.theta0, p.omega0};
    integrator::State state_ref = state;
    
    double t = s.t_start;
    const double dt_ref = s.dt / static_cast<double>(settings.error_reference_factor);
    double t_ref = s.t_start;
    
    int nsteps = static_cast<int>((s.t_end - s.t_start) / s.dt + 0.5);
    result.rk4_steps = nsteps;

    auto derivs = [&p, omega0_sq](double t_val, const integrator::State& state_val) -> integrator::State {
        return {state_val.omega, -damping_term(state_val.theta, state_val.omega, p)
                               - omega0_sq * restoring_force::term(state_val.theta, p.restoring_force)
                               + p.A * std::cos(p.omega_drive * t_val)};
    };
    const double den_gamma = 0.5 * p.damping;
    const double den_omega0_sq = std::max(
        0.0, omega0_sq * restoring_force::linearized_slope(p.restoring_force));
    const double den_omega0 = std::sqrt(den_omega0_sq);
    const std::function<double(double, const integrator::State&)> den_residual =
        [&p, omega0_sq, den_omega0_sq](double t_val, const integrator::State& state_val) {
            return -damping_term(state_val.theta, state_val.omega, p)
                   - omega0_sq * restoring_force::term(state_val.theta, p.restoring_force)
                   + p.A * std::cos(p.omega_drive * t_val)
                   + p.damping * state_val.omega
                   + den_omega0_sq * state_val.theta;
        };
    auto acceleration = [&p, omega0_sq](double t_val, double theta) {
        return -omega0_sq * restoring_force::term(theta, p.restoring_force) +
               p.A * std::cos(p.omega_drive * t_val);
    };

    if (integrator::is_position_only_integrator(settings.integrator)) {
        const std::vector<integrator::State> trajectory =
            integrator::integrate_position_only_trajectory(
                settings.integrator, s.t_start, s.dt, nsteps, state, acceleration);
        std::vector<integrator::State> trajectory_ref;
        if (settings.error_mode == error_reference::Mode::HdReference) {
            trajectory_ref = integrator::integrate_position_only_trajectory(
                settings.integrator, s.t_start, dt_ref,
                nsteps * settings.error_reference_factor, state_ref, acceleration);
        }

        for (int n = 0; n <= nsteps; ++n) {
            const double t_sample = s.t_start + n * s.dt;
            const auto& sample = trajectory[static_cast<size_t>(n)];

            double theta_exact = 0.0;
            double omega_exact = 0.0;
            if (settings.error_mode == error_reference::Mode::None) {
                theta_exact = sample.theta;
                omega_exact = sample.omega;
            } else if (settings.error_mode == error_reference::Mode::HdReference) {
                const auto& ref_sample =
                    trajectory_ref[static_cast<size_t>(n * settings.error_reference_factor)];
                theta_exact = ref_sample.theta;
                omega_exact = ref_sample.omega;
            } else {
                analytical_solution(t_sample, p, omega_lin_sq, theta_exact, omega_exact);
            }

            result.t.push_back(t_sample);
            result.theta.push_back(sample.theta);
            result.omega.push_back(sample.omega);
            result.theta_analytical.push_back(theta_exact);
            result.omega_analytical.push_back(omega_exact);
            result.energy.push_back(
                total_energy(sample.theta, sample.omega, omega0_sq, p.restoring_force));
        }

        compute_error_statistics(result);
        return result;
    }

    for (int n = 0; n <= nsteps; ++n) {
        double theta_exact = 0.0;
        double omega_exact = 0.0;
        if (settings.error_mode == error_reference::Mode::None) {
            theta_exact = state.theta;
            omega_exact = state.omega;
        } else if (settings.error_mode == error_reference::Mode::HdReference) {
            theta_exact = state_ref.theta;
            omega_exact = state_ref.omega;
        } else {
            analytical_solution(t, p, omega_lin_sq, theta_exact, omega_exact);
        }
        
        result.t.push_back(t);
        result.theta.push_back(state.theta);
        result.omega.push_back(state.omega);
        result.theta_analytical.push_back(theta_exact);
        result.omega_analytical.push_back(omega_exact);
        result.energy.push_back(
            total_energy(state.theta, state.omega, omega0_sq, p.restoring_force));

        if (n < nsteps) {
            state = advance_state(settings.integrator, t, s.dt, state, derivs,
                                  den_gamma, den_omega0, den_residual);
            if (settings.error_mode == error_reference::Mode::HdReference) {
                for (int k = 0; k < settings.error_reference_factor; ++k) {
                    state_ref = advance_state(settings.integrator, t_ref, dt_ref,
                                              state_ref, derivs,
                                              den_gamma, den_omega0, den_residual);
                    t_ref += dt_ref;
                }
            }
            
            t += s.dt;
        }
    }

    compute_error_statistics(result);

    return result;
}
