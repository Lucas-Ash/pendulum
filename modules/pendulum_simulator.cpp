#include "modules/pendulum_simulator.h"
#include "modules/rk_integrators.h"
#include "modules/error_analysis.h"
#include "modules/jacobi_elliptic.h"

#include <algorithm>
#include <cmath>
#include <functional>
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
    if (method == "rk3") {
        return integrator::rk3_step(t, dt, state, derivs);
    }
    if (method == "rk5") {
        return integrator::rk5_step(t, dt, state, derivs);
    }
    if (method == "semi_implicit_euler") {
        return integrator::semi_implicit_euler_step(t, dt, state, derivs);
    }
    if (method == "leapfrog") {
        return integrator::leapfrog_step(t, dt, state, derivs);
    }
    if (method == "ruth4") {
        return integrator::ruth4_step(t, dt, state, derivs);
    }
    if (method == "rk23") {
        return integrator::rk23_step(t, dt, state, derivs);
    }
    if (method == "rkf45") {
        return integrator::rkf45_step(t, dt, state, derivs);
    }
    if (integrator::is_den_method(method)) {
        return integrator::den3_step(t, dt, state, den_gamma, den_omega0, den_residual);
    }
    return integrator::rk4_step(t, dt, state, derivs);
}

}  // namespace

PendulumSimulator::PendulumSimulator(double length, double gravity,
                                     double timestep, double tmax,
                                     restoring_force::Config restoring_force)
    : g_(gravity), L_(length), dt_(timestep), t_max_(tmax),
      restoring_force_(restoring_force) {}

double PendulumSimulator::gravity() const { return g_; }
double PendulumSimulator::length() const { return L_; }
double PendulumSimulator::dt() const { return dt_; }
double PendulumSimulator::t_max() const { return t_max_; }


void PendulumSimulator::exact_harmonic_state(double t, double theta0, double omega0,
                                             double omega_natural,
                                             double& theta_exact,
                                             double& omega_exact) {
    if (omega_natural < 1e-15) {
        theta_exact = theta0 + omega0 * t;
        omega_exact = omega0;
        return;
    }

    const double c = std::cos(omega_natural * t);
    const double s = std::sin(omega_natural * t);
    theta_exact = theta0 * c + (omega0 / omega_natural) * s;
    omega_exact = -theta0 * omega_natural * s + omega0 * c;
}

void PendulumSimulator::exact_linear_state(double t, double theta0, double omega0,
                                           double& theta_exact, double& omega_exact) const {
    const double wn_sq = (g_ / L_) * restoring_force::linearized_slope(restoring_force_);
    if (wn_sq < 0.0) {
        throw std::runtime_error(
            "Linear analytical model requires positive linearized restoring stiffness.");
    }
    exact_harmonic_state(t, theta0, omega0, std::sqrt(std::max(wn_sq, 0.0)),
                         theta_exact, omega_exact);
}

void PendulumSimulator::exact_nonlinear_state(double t, double theta0, double omega0,
                                              double& theta_exact, double& omega_exact) const {
    if (restoring_force_.model == restoring_force::Model::Polynomial) {
        exact_duffing_state(t, theta0, omega0, theta_exact, omega_exact);
        return;
    }

    const double wn = std::sqrt(g_ / L_);
    
    // Total dimensionless energy E / (m*g*L)
    const double E = 0.5 * (omega0 / wn) * (omega0 / wn) + (1.0 - std::cos(theta0));
    
    // For bounded motion (oscillation), E < 2
    if (E >= 2.0) {
        // Fallback to linear if the energy causes continuous rotation (unbounded).
        // A full treatment handles m > 1 via reciprocal modulus conversions, 
        // but for a "simple pendulum" demonstration, let's bound to oscillation.
        exact_linear_state(t, theta0, omega0, theta_exact, omega_exact);
        return;
    }

    const double k_sq = E / 2.0;  // m parameter for Jacobi functions
    const double k = std::sqrt(k_sq);
    
    // Modulus bounds
    if (k_sq < 1e-15) {
        exact_linear_state(t, theta0, omega0, theta_exact, omega_exact);
        return;
    }

    // Determine phase offset phi_0
    double sin_theta_half = std::sin(theta0 / 2.0) / k;
    sin_theta_half = std::max(-1.0, std::min(1.0, sin_theta_half));
    
    // Fix the initial v0 determination.
    // Given theta_0 and omega_0, we can define the initial state directly via sn, cn, dn.
    // Instead of wrestling with inverse elliptic integrals over branches,
    // we solve for the time offset t_0 equivalent to the state.
    // Actually, std::ellint_1(k, phi) is F(phi, k). The k parameter is the first argument in C++17.
    // However, it's easier to use the fact that if omega_0 = 0 (like our default cases), 
    // the pendulum is released from rest at theta0.
    // Hence, E = 1 - cos(theta0) and k = sin(theta0/2).
    // The maximum angle is exactly theta0, so the phase starts at F(pi/2, k) which is the quarter period K(k).
    
    double phi_0 = 0.0;
    if (std::abs(sin_theta_half) >= 1.0) {
        phi_0 = M_PI / 2.0 * (sin_theta_half > 0 ? 1.0 : -1.0);
    } else {
        phi_0 = std::asin(sin_theta_half);
    }
    
    // Calculate complete elliptic integral of the first kind K(k)
    double K = std::comp_ellint_1(k);
    
    // Our starting phase v0.
    // If omega0 == 0, we start at amplitude, so v0 = K.
    // If we have an initial velocity, we need F(phi_0, k).
    // Note C++17 signature is ellint_1(k, phi).
    double v0 = std::ellint_1(k, phi_0);
    if (omega0 < 0) {
        v0 = K + (K - v0);
    }

    // Current phase variable v_t
    double v_t = wn * t + v0;

    double sn, cn, dn;
    math_utils::jacobi_sn_cn_dn(v_t, k_sq, sn, cn, dn);

    theta_exact = 2.0 * std::asin(k * sn);
    omega_exact = 2.0 * wn * k * cn;
}

void PendulumSimulator::exact_duffing_state(double t, double theta0, double omega0,
                                            double& theta_exact, double& omega_exact) const {
    if (restoring_force_.model != restoring_force::Model::Polynomial) {
        throw std::runtime_error(
            "Duffing analytical solution requires polynomial restoring force.");
    }

    const double alpha = (g_ / L_) * restoring_force_.linear;
    const double beta = (g_ / L_) * restoring_force_.cubic;

    if (std::abs(beta) < 1e-15) {
        if (alpha < 0.0) {
            throw std::runtime_error(
                "Duffing analytical fallback requires non-negative linear coefficient.");
        }
        exact_harmonic_state(t, theta0, omega0, std::sqrt(std::max(alpha, 0.0)),
                             theta_exact, omega_exact);
        return;
    }
    if (beta < 0.0) {
        throw std::runtime_error(
            "Duffing Jacobi analytical solution currently requires restoring_force_cubic > 0.");
    }
    if (alpha < 0.0) {
        throw std::runtime_error(
            "Duffing Jacobi analytical solution currently requires restoring_force_linear >= 0.");
    }

    const double theta0_sq = theta0 * theta0;
    const double energy = 0.5 * omega0 * omega0
                        + 0.5 * alpha * theta0_sq
                        + 0.25 * beta * theta0_sq * theta0_sq;

    const double discr = alpha * alpha + 4.0 * beta * energy;
    if (discr < 0.0) {
        throw std::runtime_error(
            "Duffing analytical solution failed: negative turning-point discriminant.");
    }

    const double amplitude_sq = (-alpha + std::sqrt(discr)) / beta;
    if (amplitude_sq < -1e-12) {
        throw std::runtime_error(
            "Duffing analytical solution failed: computed negative amplitude^2.");
    }

    const double amplitude = std::sqrt(std::max(0.0, amplitude_sq));
    if (amplitude < 1e-15) {
        theta_exact = 0.0;
        omega_exact = 0.0;
        return;
    }

    const double omega_sq = alpha + beta * amplitude_sq;
    if (omega_sq <= 0.0) {
        throw std::runtime_error(
            "Duffing analytical solution requires alpha + beta*A^2 > 0.");
    }
    const double omega = std::sqrt(omega_sq);

    const double m = (beta * amplitude_sq) / (2.0 * omega_sq);
    if (m < -1e-12 || m > 1.0 + 1e-12) {
        throw std::runtime_error(
            "Duffing analytical solution produced Jacobi parameter outside [0, 1].");
    }
    const double m_clamped = std::max(0.0, std::min(1.0, m));

    double cn0 = theta0 / amplitude;
    cn0 = std::max(-1.0, std::min(1.0, cn0));
    const double phi0 = std::acos(cn0);
    const double u_abs = std::ellint_1(std::sqrt(m_clamped), phi0);
    const double u0 = (omega0 >= 0.0) ? -u_abs : u_abs;

    double sn = 0.0;
    double cn = 1.0;
    double dn = 1.0;
    math_utils::jacobi_sn_cn_dn(omega * t + u0, m_clamped, sn, cn, dn);

    theta_exact = amplitude * cn;
    omega_exact = -amplitude * omega * sn * dn;
}

SimulationResult PendulumSimulator::simulate(
    double theta0,
    double omega0,
    const std::string& integrator,
    const std::string& analytical_model,
    error_reference::Mode error_mode,
    int error_reference_factor) const {
    SimulationResult result;

    if (error_reference_factor <= 0) {
        throw std::runtime_error("error_reference_factor must be > 0");
    }
    if (error_mode == error_reference::Mode::HdReference &&
        error_reference_factor < 2) {
        throw std::runtime_error(
            "hd_reference mode requires error_reference_factor >= 2");
    }

    const int nsteps = static_cast<int>(t_max_ / dt_ + 0.5);
    const double dt_ref = dt_ / static_cast<double>(error_reference_factor);
    double t_ref = 0.0;

    integrator::State state = {theta0, omega0};
    integrator::State state_ref = state;

    auto derivs = [this](double t_val, const integrator::State& s) -> integrator::State {
        return {s.omega,
                -(g_ / L_) * restoring_force::term(s.theta, restoring_force_)};
    };
    const double den_omega0_sq = std::max(
        0.0, (g_ / L_) * restoring_force::linearized_slope(restoring_force_));
    const double den_omega0 = std::sqrt(den_omega0_sq);
    const std::function<double(double, const integrator::State&)> den_residual =
        [this, den_omega0_sq](double, const integrator::State& s) {
            return -(g_ / L_) * restoring_force::term(s.theta, restoring_force_) +
                   den_omega0_sq * s.theta;
        };
    auto acceleration = [this](double, double theta) {
        return -(g_ / L_) * restoring_force::term(theta, restoring_force_);
    };

    if (integrator::is_position_only_integrator(integrator)) {
        const std::vector<integrator::State> trajectory =
            integrator::integrate_position_only_trajectory(
                integrator, 0.0, dt_, nsteps, state, acceleration);
        std::vector<integrator::State> trajectory_ref;
        if (error_mode == error_reference::Mode::HdReference) {
            trajectory_ref = integrator::integrate_position_only_trajectory(
                integrator, 0.0, dt_ref, nsteps * error_reference_factor, state_ref,
                acceleration);
        }

        for (int n = 0; n <= nsteps; ++n) {
            const double t = n * dt_;
            const auto& sample = trajectory[static_cast<size_t>(n)];

            result.t.push_back(t);
            result.theta.push_back(sample.theta);
            result.omega.push_back(sample.omega);

            double theta_ref_val = 0.0;
            double omega_ref_val = 0.0;
            if (error_mode == error_reference::Mode::None) {
                theta_ref_val = sample.theta;
                omega_ref_val = sample.omega;
            } else if (error_mode == error_reference::Mode::HdReference) {
                const auto& ref_sample =
                    trajectory_ref[static_cast<size_t>(n * error_reference_factor)];
                theta_ref_val = ref_sample.theta;
                omega_ref_val = ref_sample.omega;
            } else {
                const bool wants_duffing_analytical =
                    analytical_model == "duffing_jacobi" ||
                    analytical_model == "duffing";
                if (wants_duffing_analytical ||
                    (analytical_model == "jacobi" &&
                     restoring_force_.model == restoring_force::Model::Polynomial)) {
                    exact_duffing_state(t, theta0, omega0, theta_ref_val, omega_ref_val);
                } else if (analytical_model == "jacobi") {
                    exact_nonlinear_state(t, theta0, omega0, theta_ref_val, omega_ref_val);
                } else {
                    exact_linear_state(t, theta0, omega0, theta_ref_val, omega_ref_val);
                }
            }

            result.theta_analytical.push_back(theta_ref_val);
            result.omega_analytical.push_back(omega_ref_val);
            result.energy.push_back(
                0.5 * sample.omega * sample.omega +
                (g_ / L_) * restoring_force::potential(sample.theta, restoring_force_));
        }

        compute_error_statistics(result);
        result.rk4_steps = nsteps;
        return result;
    }

    for (int n = 0; n <= nsteps; ++n) {
        const double t = n * dt_;
        const double theta = state.theta;
        const double omega = state.omega;

        result.t.push_back(t);
        result.theta.push_back(theta);
        result.omega.push_back(omega);
        
        double theta_ref = 0.0;
        double omega_ref = 0.0;
        if (error_mode == error_reference::Mode::None) {
            theta_ref = theta;
            omega_ref = omega;
        } else if (error_mode == error_reference::Mode::HdReference) {
            theta_ref = state_ref.theta;
            omega_ref = state_ref.omega;
        } else {
            const bool wants_duffing_analytical =
                analytical_model == "duffing_jacobi" ||
                analytical_model == "duffing";
            if (wants_duffing_analytical ||
                (analytical_model == "jacobi" &&
                 restoring_force_.model == restoring_force::Model::Polynomial)) {
                exact_duffing_state(t, theta0, omega0, theta_ref, omega_ref);
            } else if (analytical_model == "jacobi") {
                exact_nonlinear_state(t, theta0, omega0, theta_ref, omega_ref);
            } else {
                exact_linear_state(t, theta0, omega0, theta_ref, omega_ref);
            }
        }
        result.theta_analytical.push_back(theta_ref);
        result.omega_analytical.push_back(omega_ref);
        
        double energy = 0.5 * omega * omega +
                        (g_ / L_) *
                        restoring_force::potential(theta, restoring_force_);
        result.energy.push_back(energy);

        if (n < nsteps) {
            state = advance_state(integrator, t, dt_, state, derivs,
                                  0.0, den_omega0, den_residual);
            if (error_mode == error_reference::Mode::HdReference) {
                for (int k = 0; k < error_reference_factor; ++k) {
                    state_ref = advance_state(
                        integrator, t_ref, dt_ref, state_ref, derivs,
                        0.0, den_omega0, den_residual);
                    t_ref += dt_ref;
                }
            }
        }
    }

    compute_error_statistics(result);
    result.rk4_steps = nsteps;

    return result;
}
