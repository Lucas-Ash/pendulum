#include "modules/damped/damped_simulator.h"
#include "modules/integrators/simulation_runner.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace {

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
        std::abs(physical.damping_cubic) > 1e-12 ||
        additional_terms::has_velocity_dependence(physical.additional_terms)) {
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

    const double omega0_sq = p.g / p.L;
    const bool use_analytical = settings.error_mode == error_reference::Mode::Analytical;
    double omega_d = 0.0;
    double omega_linear = 0.0;
    double omega_lin_sq = 0.0;
    double van_der_pol_mu = 0.0;
    double van_der_pol_phase = 0.0;
    const bool wants_lane_emden =
        settings.analytical_model == "lane_emden";
    if (use_analytical) {
        const bool wants_van_der_pol =
            settings.analytical_model == "van_der_pol" ||
            settings.analytical_model == "van_der_pol_first_order";
        if (wants_lane_emden) {
            if (std::abs(linear_damping_coefficient(p)) > 1e-12 ||
                std::abs(p.damping_cubic) > 1e-12) {
                throw std::runtime_error(
                    "Lane-Emden analytical reference requires classical damping to be disabled.");
            }
            if (p.restoring_force.model != restoring_force::Model::Polynomial ||
                std::abs(p.restoring_force.linear) > 1e-12 ||
                std::abs(p.restoring_force.cubic) > 1e-12) {
                throw std::runtime_error(
                    "Lane-Emden analytical reference requires the base restoring force to be disabled.");
            }
            if (!p.additional_terms.time_damping_enabled ||
                std::abs(p.additional_terms.time_damping_coefficient - 2.0) > 1e-12 ||
                std::abs(p.additional_terms.time_damping_power - 1.0) > 1e-12) {
                throw std::runtime_error(
                    "Lane-Emden analytical reference requires time_damping_coefficient = 2 and time_damping_power = 1.");
            }
            if (!p.additional_terms.state_power_enabled ||
                std::abs(p.additional_terms.state_power_strength - 1.0) > 1e-12) {
                throw std::runtime_error(
                    "Lane-Emden analytical reference requires state_power_strength = 1.");
            }
            const double lane_index = p.additional_terms.state_power_exponent;
            const bool supported_index =
                std::abs(lane_index) < 1e-12 ||
                std::abs(lane_index - 1.0) < 1e-12 ||
                std::abs(lane_index - 5.0) < 1e-12;
            if (!supported_index) {
                throw std::runtime_error(
                    "Lane-Emden analytical reference is implemented for polytropic indices n = 0, 1, or 5.");
            }
            if (s.t_start + p.additional_terms.time_damping_shift <=
                p.additional_terms.singularity_epsilon) {
                throw std::runtime_error(
                    "Lane-Emden analytical reference requires t_start + time_damping_shift > 0.");
            }
        } else {
            const double linearized_slope =
                restoring_force::linearized_slope(p.restoring_force) +
                additional_terms::linearized_stiffness(p.additional_terms) / omega0_sq;
            omega_lin_sq = omega0_sq * linearized_slope;
            if (omega_lin_sq <= 0.0) {
                throw std::runtime_error(
                    "Damped analytical reference requires positive linearized stiffness.");
            }
            omega_linear = std::sqrt(omega_lin_sq);
        }

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
        } else if (!wants_lane_emden) {
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
    const integrator::State initial_state = {p.theta0, p.theta_dot0};

    auto derivs = [physical = p, omega0_sq](double t_val, const integrator::State& state) -> integrator::State {
        return {state.omega, -damping_term(state.theta, state.omega, physical)
                              - omega0_sq * restoring_force::term(state.theta, physical.restoring_force)
                              + additional_terms::acceleration(
                                    t_val, state.theta, state.omega, physical.additional_terms)};
    };
    const double den_omega0_sq = std::max(
        0.0, omega0_sq * restoring_force::linearized_slope(p.restoring_force) +
                 additional_terms::linearized_stiffness(p.additional_terms));
    const double den_omega0 = std::sqrt(den_omega0_sq);
    auto den_residual =
        [physical = p, omega0_sq, den_omega0_sq](double t_val, const integrator::State& state) {
            return -damping_term(state.theta, state.omega, physical)
                   - omega0_sq * restoring_force::term(state.theta, physical.restoring_force)
                   + additional_terms::acceleration(
                         t_val, state.theta, state.omega, physical.additional_terms)
                   + linear_damping_coefficient(physical) * state.omega
                   + den_omega0_sq * state.theta;
        };
    auto acceleration = [physical = p, omega0_sq](double t_val, double theta) {
        return -omega0_sq * restoring_force::term(theta, physical.restoring_force) +
               additional_terms::acceleration(
                   t_val, theta, 0.0, physical.additional_terms);
    };

    auto reference_state =
        [settings, omega_linear, van_der_pol_mu, van_der_pol_phase, omega_d, p](
            double t,
            const integrator::State& sample,
            const integrator::State* hd_reference_state,
            double& theta_reference,
            double& omega_reference) {
            if (settings.error_mode == error_reference::Mode::None) {
                theta_reference = sample.theta;
                omega_reference = sample.omega;
                return;
            }
            if (settings.error_mode == error_reference::Mode::HdReference) {
                theta_reference = hd_reference_state->theta;
                omega_reference = hd_reference_state->omega;
                return;
            }

            const bool wants_van_der_pol =
                settings.analytical_model == "van_der_pol" ||
                settings.analytical_model == "van_der_pol_first_order";
            if (settings.analytical_model == "lane_emden") {
                const double xi = t + p.additional_terms.time_damping_shift;
                const double index = p.additional_terms.state_power_exponent;
                if (std::abs(index) < 1e-12) {
                    theta_reference = 1.0 - xi * xi / 6.0;
                    omega_reference = -xi / 3.0;
                } else if (std::abs(index - 1.0) < 1e-12) {
                    theta_reference = std::sin(xi) / xi;
                    omega_reference = (xi * std::cos(xi) - std::sin(xi)) / (xi * xi);
                } else {
                    const double scale = 1.0 + xi * xi / 3.0;
                    theta_reference = 1.0 / std::sqrt(scale);
                    omega_reference = -xi / (3.0 * std::pow(scale, 1.5));
                }
                return;
            }
            if (wants_van_der_pol) {
                van_der_pol_first_order_state(
                    t, omega_linear, van_der_pol_mu, van_der_pol_phase,
                    theta_reference, omega_reference);
                return;
            }

            linear_analytical_solution(
                t, omega_d, p.gamma, p.theta0, p.theta_dot0,
                theta_reference, omega_reference);
        };

    auto energy_of = [omega0_sq, p](const integrator::State& sample) {
        return total_energy(sample.theta, sample.omega, omega0_sq, p.restoring_force) +
               additional_terms::potential(sample.theta, p.additional_terms);
    };

    return simulation_runner::run(
        {
            settings.integrator,
            s.t_start,
            s.dt,
            nsteps,
            settings.error_mode,
            settings.error_reference_factor,
            p.gamma,
            den_omega0,
        },
        initial_state,
        derivs,
        acceleration,
        den_residual,
        reference_state,
        energy_of);
}
