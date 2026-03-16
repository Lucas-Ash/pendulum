#include "modules/damped_simulator.h"
#include "modules/rk_integrators.h"
#include "modules/error_analysis.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>

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

double linear_damping_coefficient(const DampedPhysicalConfig& physical) {
    return 2.0 * physical.gamma;
}

double damping_term(double theta, double omega, const DampedPhysicalConfig& physical) {
    return damping_force::term(
        theta,
        omega,
        physical.damping_model,
        linear_damping_coefficient(physical),
        physical.damping_cubic);
}

void validate_position_only_integrator(const std::string& method,
                                       const DampedPhysicalConfig& physical) {
    if (!integrator::is_position_only_integrator(method)) {
        return;
    }
    if (std::abs(linear_damping_coefficient(physical)) > 1e-12 ||
        std::abs(physical.damping_cubic) > 1e-12) {
        throw std::runtime_error(
            method + " requires theta'' = a(t, theta); it is not valid with velocity-dependent damping.");
    }
}

void van_der_pol_first_order_state(double t,
                                   double omega_linear,
                                   double mu,
                                   double phase,
                                   double& theta_exact,
                                   double& omega_exact) {
    const double tau = omega_linear * t + phase;
    const double epsilon = mu / omega_linear;

    theta_exact =
        2.0 * std::cos(tau) +
        epsilon * (0.75 * std::sin(tau) - 0.25 * std::sin(3.0 * tau));
    omega_exact =
        -2.0 * omega_linear * std::sin(tau) +
        mu * (0.75 * std::cos(tau) - 0.75 * std::cos(3.0 * tau));
}

double estimate_van_der_pol_phase(double omega_linear,
                                  double mu,
                                  double theta0,
                                  double theta_dot0) {
    constexpr int samples = 4096;
    constexpr double two_pi = 2.0 * M_PI;

    double best_phase = 0.0;
    double best_residual = std::numeric_limits<double>::infinity();

    for (int i = 0; i < samples; ++i) {
        const double phase = two_pi * static_cast<double>(i) / static_cast<double>(samples);
        double theta_ref = 0.0;
        double omega_ref = 0.0;
        van_der_pol_first_order_state(0.0, omega_linear, mu, phase, theta_ref, omega_ref);

        const double theta_residual = theta_ref - theta0;
        const double omega_residual = (omega_ref - theta_dot0) / std::max(omega_linear, 1e-12);
        const double residual = theta_residual * theta_residual +
                                omega_residual * omega_residual;
        if (residual < best_residual) {
            best_residual = residual;
            best_phase = phase;
        }
    }

    if (best_residual > 1e-4) {
        throw std::runtime_error(
            "Van der Pol analytical reference expects initial conditions on the weakly nonlinear limit cycle.");
    }
    return best_phase;
}

}  // namespace

DampedPendulumSimulator::DampedPendulumSimulator(const DampedConfig& config)
    : config_(config) {}

void DampedPendulumSimulator::linear_analytical_solution(double t, double omega_d, double gamma,
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
    double omega_d = 0.0;
    double omega_linear = 0.0;
    double van_der_pol_mu = 0.0;
    double van_der_pol_phase = 0.0;
    if (use_analytical) {
        const double linearized_slope = restoring_force::linearized_slope(p.restoring_force);
        const double omega_lin_sq = omega0_sq * linearized_slope;
        if (omega_lin_sq <= 0.0) {
            throw std::runtime_error(
                "Damped analytical reference requires positive linearized stiffness.");
        }
        omega_linear = std::sqrt(omega_lin_sq);

        const bool wants_van_der_pol =
            settings.analytical_model == "van_der_pol" ||
            settings.analytical_model == "van_der_pol_first_order";
        if (wants_van_der_pol) {
            if (p.damping_model != damping_force::Model::Polynomial) {
                throw std::runtime_error(
                    "Van der Pol analytical reference requires damping_model: polynomial.");
            }
            if (p.restoring_force.model != restoring_force::Model::Polynomial ||
                std::abs(p.restoring_force.cubic) > 1e-12) {
                throw std::runtime_error(
                    "Van der Pol analytical reference requires a purely linear polynomial restoring force.");
            }

            const double linear_damping = linear_damping_coefficient(p);
            van_der_pol_mu = p.damping_cubic;
            if (van_der_pol_mu <= 0.0) {
                throw std::runtime_error(
                    "Van der Pol analytical reference requires damping_cubic > 0.");
            }
            if (std::abs(linear_damping + van_der_pol_mu) > 1e-12) {
                throw std::runtime_error(
                    "Van der Pol analytical reference requires the quadratic damping coefficient "
                    "-mu*(1-theta^2), i.e. damping_linear = -mu and damping_cubic = mu "
                    "(equivalently gamma = -0.5 * damping_cubic).");
            }
            if (van_der_pol_mu / omega_linear > 0.2) {
                throw std::runtime_error(
                    "Van der Pol analytical reference is only implemented for weak nonlinearity "
                    "(damping_cubic / omega_linear <= 0.2).");
            }
            van_der_pol_phase = estimate_van_der_pol_phase(
                omega_linear, van_der_pol_mu, p.theta0, p.theta_dot0);
        } else {
            if (std::abs(p.damping_cubic) > 1e-15) {
                throw std::runtime_error(
                    "Linear damped analytical reference does not support damping_cubic != 0. "
                    "Use analytical_model: van_der_pol or error_analysis: hd_reference.");
            }
            if (p.gamma >= omega_linear) {
                throw std::runtime_error(
                    "This analytical model assumes underdamped motion: gamma < omega_linearized");
            }
            omega_d = std::sqrt(omega_lin_sq - p.gamma * p.gamma);
        }
    }

    const int nsteps = static_cast<int>((s.t_end - s.t_start) / s.dt + 0.5);
    const double dt_ref = s.dt / static_cast<double>(settings.error_reference_factor);
    double t_ref = s.t_start;
    result.rk4_steps = nsteps;

    auto derivs = [physical = p, omega0_sq](double t_val, const integrator::State& state) -> integrator::State {
        return {state.omega, -damping_term(state.theta, state.omega, physical)
                              - omega0_sq * restoring_force::term(state.theta, physical.restoring_force)};
    };
    const double den_omega0_sq = std::max(
        0.0, omega0_sq * restoring_force::linearized_slope(p.restoring_force));
    const double den_omega0 = std::sqrt(den_omega0_sq);
    const std::function<double(double, const integrator::State&)> den_residual =
        [physical = p, omega0_sq, den_omega0_sq](double t_val, const integrator::State& state) {
            return -damping_term(state.theta, state.omega, physical)
                   - omega0_sq * restoring_force::term(state.theta, physical.restoring_force)
                   + linear_damping_coefficient(physical) * state.omega
                   + den_omega0_sq * state.theta;
        };
    auto acceleration = [physical = p, omega0_sq](double t_val, double theta) {
        return -omega0_sq * restoring_force::term(theta, physical.restoring_force);
    };

    integrator::State state = {p.theta0, p.theta_dot0};
    integrator::State state_ref = state;

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
            const double t = s.t_start + n * s.dt;
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
                const bool wants_van_der_pol =
                    settings.analytical_model == "van_der_pol" ||
                    settings.analytical_model == "van_der_pol_first_order";
                if (wants_van_der_pol) {
                    van_der_pol_first_order_state(
                        t, omega_linear, van_der_pol_mu, van_der_pol_phase,
                        theta_exact, omega_exact);
                } else {
                    linear_analytical_solution(
                        t, omega_d, p.gamma, p.theta0, p.theta_dot0,
                        theta_exact, omega_exact);
                }
            }

            result.t.push_back(t);
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
            const bool wants_van_der_pol =
                settings.analytical_model == "van_der_pol" ||
                settings.analytical_model == "van_der_pol_first_order";
            if (wants_van_der_pol) {
                van_der_pol_first_order_state(
                    t, omega_linear, van_der_pol_mu, van_der_pol_phase,
                    theta_exact, omega_exact);
            } else {
                linear_analytical_solution(
                    t, omega_d, p.gamma, p.theta0, p.theta_dot0,
                    theta_exact, omega_exact);
            }
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
                                  p.gamma, den_omega0, den_residual);
            if (settings.error_mode == error_reference::Mode::HdReference) {
                for (int k = 0; k < settings.error_reference_factor; ++k) {
                    state_ref = advance_state(settings.integrator, t_ref, dt_ref,
                                              state_ref, derivs,
                                              p.gamma, den_omega0, den_residual);
                    t_ref += dt_ref;
                }
            }
        }
    }

    compute_error_statistics(result);

    return result;
}
