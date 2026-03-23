#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace integrator {

struct State {
    double theta;
    double omega;

    State operator+(const State& other) const {
        return {theta + other.theta, omega + other.omega};
    }

    State operator*(double scalar) const {
        return {theta * scalar, omega * scalar};
    }
};

inline bool is_velocity_verlet_method(const std::string& method) {
    return method == "velocity_verlet" || method == "velocity-verlet" ||
           method == "verlet";
}

inline bool is_runge_kutta_nystrom_method(const std::string& method) {
    return method == "runge_kutta_nystrom" || method == "runge-kutta-nystrom" ||
           method == "rkn4";
}

inline bool is_numerov_method(const std::string& method) {
    return method == "numerov";
}

inline bool is_den_method(const std::string& method) {
    return method == "den3" || method == "den";
}

inline bool is_position_only_integrator(const std::string& method) {
    return is_velocity_verlet_method(method) ||
           is_runge_kutta_nystrom_method(method) ||
           is_numerov_method(method);
}

struct Propagator {
    double phi11;
    double phi12;
    double phi21;
    double phi22;
};

inline Propagator den_propagator(double tau, double gamma, double omega0) {
    const double omega0_sq = omega0 * omega0;
    const double discriminant = omega0_sq - gamma * gamma;
    const double threshold = 1e-10 * std::max(omega0_sq, 1.0);

    if (discriminant > threshold) {
        const double omega_d = std::sqrt(discriminant);
        const double exp_factor = std::exp(-gamma * tau);
        const double cos_term = std::cos(omega_d * tau);
        const double sin_term = std::sin(omega_d * tau);
        return {
            exp_factor * (cos_term + (gamma / omega_d) * sin_term),
            exp_factor * (sin_term / omega_d),
            exp_factor * (-(omega0_sq / omega_d) * sin_term),
            exp_factor * (cos_term - (gamma / omega_d) * sin_term),
        };
    }

    if (discriminant < -threshold) {
        const double omega_h = std::sqrt(-discriminant);
        const double exp_factor = std::exp(-gamma * tau);
        const double cosh_term = std::cosh(omega_h * tau);
        const double sinh_term = std::sinh(omega_h * tau);
        return {
            exp_factor * (cosh_term + (gamma / omega_h) * sinh_term),
            exp_factor * (sinh_term / omega_h),
            exp_factor * (-(omega0_sq / omega_h) * sinh_term),
            exp_factor * (cosh_term - (gamma / omega_h) * sinh_term),
        };
    }

    const double exp_factor = std::exp(-gamma * tau);
    return {
        exp_factor * (1.0 + gamma * tau),
        exp_factor * tau,
        exp_factor * (-omega0_sq * tau),
        exp_factor * (1.0 - gamma * tau),
    };
}

inline State den_constant_residual_response(double tau, double gamma, double omega0) {
    if (std::abs(omega0) < 1e-14) {
        if (std::abs(gamma) < 1e-14) {
            return {0.5 * tau * tau, tau};
        }
        const double exp_term = std::exp(-2.0 * gamma * tau);
        return {
            tau / (2.0 * gamma) - (1.0 - exp_term) / (4.0 * gamma * gamma),
            (1.0 - exp_term) / (2.0 * gamma),
        };
    }

    const Propagator phi = den_propagator(tau, gamma, omega0);
    return {
        (1.0 - phi.phi22 - 2.0 * gamma * phi.phi12) / (omega0 * omega0),
        phi.phi12,
    };
}

// A generic Runge-Kutta 4th Order integrator step function.
// DerivFunc is a callable that takes (time t, const State& current_state) 
// and returns the State derivative.
template <typename StateType, typename DerivFunc>
StateType rk4_step(double t, double dt, const StateType& state, DerivFunc derivatives) {
    StateType k1 = derivatives(t, state);
    StateType k2 = derivatives(t + dt * 0.5, state + k1 * (dt * 0.5));
    StateType k3 = derivatives(t + dt * 0.5, state + k2 * (dt * 0.5));
    StateType k4 = derivatives(t + dt, state + k3 * dt);
    
    return state + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt / 6.0);
}

// A generic Runge-Kutta 3rd Order (Standard) integrator step function.
template <typename StateType, typename DerivFunc>
StateType rk3_step(double t, double dt, const StateType& state, DerivFunc derivatives) {
    StateType k1 = derivatives(t, state);
    StateType k2 = derivatives(t + dt * 0.5, state + k1 * (dt * 0.5));
    StateType k3 = derivatives(t + dt, state + k1 * (-dt) + k2 * (2.0 * dt));
    
    return state + (k1 + k2 * 4.0 + k3) * (dt / 6.0);
}

// A generic Runge-Kutta 5th Order (Butcher's) integrator step function.
template <typename StateType, typename DerivFunc>
StateType rk5_step(double t, double dt, const StateType& state, DerivFunc derivatives) {
    StateType k1 = derivatives(t, state);
    StateType k2 = derivatives(t + dt * 0.25, state + k1 * (dt * 0.25));
    StateType k3 = derivatives(t + dt * 0.25, state + k1 * (dt * 0.125) + k2 * (dt * 0.125));
    StateType k4 = derivatives(t + dt * 0.5, state + k2 * (-dt * 0.5) + k3 * dt);
    StateType k5 = derivatives(t + dt * 0.75, state + k1 * (dt * 3.0/16.0) + k4 * (dt * 9.0/16.0));
    StateType k6 = derivatives(t + dt, state + k1 * (-dt * 3.0/7.0) + k2 * (dt * 2.0/7.0) + 
                                   k3 * (dt * 12.0/7.0) + k4 * (-dt * 12.0/7.0) + k5 * (dt * 8.0/7.0));
    
    return state + (k1 * 7.0 + k3 * 32.0 + k4 * 12.0 + k5 * 32.0 + k6 * 7.0) * (dt / 90.0);
}

// Semi-Implicit Euler (Symplectic 1st Order)
template <typename DerivFunc>
State semi_implicit_euler_step(double t, double dt, const State& state, DerivFunc derivatives) {
    State k1 = derivatives(t, state);
    State next_state = state;
    next_state.omega += k1.omega * dt;
    next_state.theta += next_state.omega * dt;
    return next_state;
}

// Leapfrog (Symplectic 2nd Order)
template <typename DerivFunc>
State leapfrog_step(double t, double dt, const State& state, DerivFunc derivatives) {
    State state_half = state;
    state_half.theta += state.omega * (dt * 0.5);
    
    State k_mid = derivatives(t + dt * 0.5, state_half);
    State next_state = state_half;
    next_state.omega += k_mid.omega * dt;
    
    next_state.theta += next_state.omega * (dt * 0.5);
    return next_state;
}

// Velocity-Verlet for theta'' = a(t, theta).
template <typename AccelFunc>
State velocity_verlet_step(double t, double dt, const State& state, AccelFunc acceleration) {
    const double a0 = acceleration(t, state.theta);
    const double theta_next = state.theta + state.omega * dt + 0.5 * a0 * dt * dt;
    const double a1 = acceleration(t + dt, theta_next);

    return {theta_next, state.omega + 0.5 * (a0 + a1) * dt};
}

// Classical 4th-order Runge-Kutta-Nystrom for theta'' = a(t, theta).
template <typename AccelFunc>
State runge_kutta_nystrom_step(double t, double dt, const State& state, AccelFunc acceleration) {
    const double k1_theta = state.omega;
    const double k1_omega = acceleration(t, state.theta);

    const double k2_theta = state.omega + k1_omega * (dt * 0.5);
    const double k2_omega =
        acceleration(t + dt * 0.5, state.theta + k1_theta * (dt * 0.5));

    const double k3_theta = state.omega + k2_omega * (dt * 0.5);
    const double k3_omega =
        acceleration(t + dt * 0.5, state.theta + k2_theta * (dt * 0.5));

    const double k4_theta = state.omega + k3_omega * dt;
    const double k4_omega =
        acceleration(t + dt, state.theta + k3_theta * dt);

    return {
        state.theta + (k1_theta + 2.0 * k2_theta + 2.0 * k3_theta + k4_theta) * (dt / 6.0),
        state.omega + (k1_omega + 2.0 * k2_omega + 2.0 * k3_omega + k4_omega) * (dt / 6.0),
    };
}

template <typename ResidualFunc>
State den3_step(double t,
                double dt,
                const State& state,
                double gamma,
                double omega0,
                ResidualFunc residual) {
    const Propagator h = den_propagator(dt, gamma, omega0);
    const Propagator h2 = den_propagator(dt * 0.5, gamma, omega0);
    const State response_h = den_constant_residual_response(dt, gamma, omega0);
    const State response_h2 = den_constant_residual_response(dt * 0.5, gamma, omega0);

    const double g1 = residual(t, state);

    const double y2 = h2.phi11 * state.theta + h2.phi12 * state.omega +
                      response_h2.theta * g1;
    const double v2 = h2.phi21 * state.theta + h2.phi22 * state.omega +
                      response_h2.omega * g1;
    const State stage2{y2, v2};
    const double g2 = residual(t + dt * 0.5, stage2);

    const double y3 = h.phi11 * state.theta + h.phi12 * state.omega +
                      response_h.theta * g2;
    const double v3 = h.phi21 * state.theta + h.phi22 * state.omega +
                      response_h.omega * g2;
    const State stage3{y3, v3};
    const double g3 = residual(t + dt, stage3);

    return {
        h.phi11 * state.theta + h.phi12 * state.omega +
            dt * (h.phi12 * (g1 / 6.0) + h2.phi12 * (2.0 * g2 / 3.0)),
        h.phi21 * state.theta + h.phi22 * state.omega +
            dt * (h.phi22 * (g1 / 6.0) + h2.phi22 * (2.0 * g2 / 3.0) + g3 / 6.0),
    };
}

template <typename AccelFunc>
long double numerov_next_theta(long double t,
                               long double dt,
                               long double prev_theta,
                               long double current_theta,
                               AccelFunc acceleration) {
    const auto eval_accel = [&](long double time, long double theta) {
        return static_cast<long double>(
            acceleration(static_cast<double>(time), static_cast<double>(theta)));
    };

    const long double dt_sq = dt * dt;
    const long double a_prev = eval_accel(t - dt, prev_theta);
    const long double a_curr = eval_accel(t, current_theta);
    const long double rhs = 2.0L * current_theta - prev_theta +
                            (dt_sq / 12.0L) * (a_prev + 10.0L * a_curr);

    long double theta_next = 2.0L * current_theta - prev_theta + a_curr * dt_sq;
    const long double scale =
        std::max({1.0L, std::fabs(prev_theta), std::fabs(current_theta)});
    const long double tol =
        std::max(128.0L * std::numeric_limits<long double>::epsilon() * scale,
                 1e-3L * dt_sq * dt_sq * dt_sq * scale);

    for (int iter = 0; iter < 20; ++iter) {
        const long double a_next = eval_accel(t + dt, theta_next);
        const long double residual = theta_next - rhs - (dt_sq / 12.0L) * a_next;
        if (std::fabs(residual) <= tol) {
            return theta_next;
        }

        const long double eps =
            std::sqrt(std::numeric_limits<long double>::epsilon()) *
            std::max(1.0L, std::fabs(theta_next));
        const long double a_plus = eval_accel(t + dt, theta_next + eps);
        const long double a_minus = eval_accel(t + dt, theta_next - eps);
        const long double daccel_dtheta = (a_plus - a_minus) / (2.0L * eps);
        const long double jacobian = 1.0L - (dt_sq / 12.0L) * daccel_dtheta;

        if (std::fabs(jacobian) < 1e-18L) {
            break;
        }

        const long double delta = residual / jacobian;
        theta_next -= delta;
        if (std::fabs(delta) <= tol) {
            return theta_next;
        }
    }

    throw std::runtime_error("Numerov solver failed to converge for theta_{n+1}.");
}

template <typename AccelFunc>
std::vector<double> numerov_consistent_velocities(double t_start,
                                                  double dt,
                                                  const std::vector<State>& states,
                                                  AccelFunc acceleration) {
    std::vector<double> omega(states.size(), 0.0);
    if (states.empty()) {
        return omega;
    }

    omega[0] = states[0].omega;
    if (states.size() == 1u) {
        return omega;
    }

    std::vector<double> accel(states.size(), 0.0);
    for (size_t i = 0; i < states.size(); ++i) {
        const double t = t_start + static_cast<double>(i) * dt;
        accel[i] = acceleration(t, states[i].theta);
    }

    omega[1] = states[1].omega;
    for (size_t i = 1; i + 1 < states.size(); ++i) {
        omega[i] = (states[i + 1].theta - states[i - 1].theta) / (2.0 * dt) -
                   dt * (accel[i + 1] - accel[i - 1]) / 12.0;
    }

    if (states.size() >= 5u) {
        const size_t n = states.size() - 1u;
        omega[n] = (25.0 * states[n].theta - 48.0 * states[n - 1u].theta +
                    36.0 * states[n - 2u].theta - 16.0 * states[n - 3u].theta +
                    3.0 * states[n - 4u].theta) / (12.0 * dt);
    } else {
        const size_t n = states.size() - 1u;
        omega[n] = states[n].omega;
    }
    return omega;
}

template <typename AccelFunc>
std::vector<State> integrate_position_only_trajectory(const std::string& method,
                                                      double t_start,
                                                      double dt,
                                                      int nsteps,
                                                      const State& initial,
                                                      AccelFunc acceleration) {
    if (nsteps < 0) {
        throw std::runtime_error("nsteps must be non-negative.");
    }

    std::vector<State> states(static_cast<size_t>(nsteps) + 1u);
    states[0] = initial;
    if (nsteps == 0) {
        return states;
    }

    if (is_velocity_verlet_method(method)) {
        double t = t_start;
        for (int i = 0; i < nsteps; ++i) {
            states[static_cast<size_t>(i) + 1u] =
                velocity_verlet_step(t, dt, states[static_cast<size_t>(i)], acceleration);
            t += dt;
        }
        return states;
    }

    if (is_runge_kutta_nystrom_method(method)) {
        double t = t_start;
        for (int i = 0; i < nsteps; ++i) {
            states[static_cast<size_t>(i) + 1u] =
                runge_kutta_nystrom_step(t, dt, states[static_cast<size_t>(i)], acceleration);
            t += dt;
        }
        return states;
    }

    if (!is_numerov_method(method)) {
        throw std::runtime_error("Unsupported position-only integrator: " + method);
    }

    if (nsteps < 4) {
        return integrate_position_only_trajectory(
            "runge_kutta_nystrom", t_start, dt, nsteps, initial, acceleration);
    }

    states[1] = runge_kutta_nystrom_step(t_start, dt, initial, acceleration);
    if (nsteps >= 2) {
        states[2] = runge_kutta_nystrom_step(t_start + dt, dt, states[1], acceleration);
    }
    std::vector<long double> theta_values(states.size(), 0.0L);
    for (size_t i = 0; i < states.size() && i < 3u; ++i) {
        theta_values[i] = static_cast<long double>(states[i].theta);
    }
    for (int i = 2; i < nsteps; ++i) {
        const long double t = static_cast<long double>(t_start) +
                              static_cast<long double>(i) * static_cast<long double>(dt);
        theta_values[static_cast<size_t>(i) + 1u] = numerov_next_theta(
            t,
            static_cast<long double>(dt),
            theta_values[static_cast<size_t>(i) - 1u],
            theta_values[static_cast<size_t>(i)],
            acceleration);
        states[static_cast<size_t>(i) + 1u] = {
            static_cast<double>(theta_values[static_cast<size_t>(i) + 1u]),
            0.0,
        };
    }

    const std::vector<double> omega =
        numerov_consistent_velocities(t_start, dt, states, acceleration);
    for (size_t i = 0; i < states.size(); ++i) {
        states[i].omega = omega[i];
    }

    return states;
}

// Ruth's 4th Order Symplectic
template <typename DerivFunc>
State ruth4_step(double t, double dt, const State& state, DerivFunc derivatives) {
    const double w = std::cbrt(2.0);
    const double w_frac = 1.0 / (2.0 - w);
    
    const double c1 = w_frac / 2.0;
    const double c2 = (1.0 - w) * w_frac / 2.0;
    const double c3 = c2;
    const double c4 = c1;
    
    const double d1 = w_frac;
    const double d2 = -w * w_frac;
    const double d3 = d1;
    
    State s = state;
    
    s.theta += s.omega * (c1 * dt);
    s.omega += derivatives(t + c1 * dt, s).omega * (d1 * dt);
    
    s.theta += s.omega * (c2 * dt);
    s.omega += derivatives(t + (c1 + c2) * dt, s).omega * (d2 * dt);
    
    s.theta += s.omega * (c3 * dt);
    s.omega += derivatives(t + (c1 + c2 + c3) * dt, s).omega * (d3 * dt);
    
    s.theta += s.omega * (c4 * dt);
    
    return s;
}

// Bogacki-Shampine (RK23) explicit step
template <typename StateType, typename DerivFunc>
StateType rk23_step(double t, double dt, const StateType& state, DerivFunc derivatives) {
    StateType k1 = derivatives(t, state);
    StateType k2 = derivatives(t + dt * 0.5, state + k1 * (dt * 0.5));
    StateType k3 = derivatives(t + dt * 0.75, state + k2 * (dt * 0.75));
    
    // We use the 3rd-order step to advance the solution
    return state + (k1 * (2.0/9.0) + k2 * (1.0/3.0) + k3 * (4.0/9.0)) * dt;
}

// Runge-Kutta-Fehlberg (RKF45) explicit step
template <typename StateType, typename DerivFunc>
StateType rkf45_step(double t, double dt, const StateType& state, DerivFunc derivatives) {
    StateType k1 = derivatives(t, state);
    StateType k2 = derivatives(t + dt * 0.25, state + k1 * (dt * 0.25));
    StateType k3 = derivatives(t + dt * (3.0/8.0), state + k1 * (dt * 3.0/32.0) + k2 * (dt * 9.0/32.0));
    StateType k4 = derivatives(t + dt * (12.0/13.0), state + k1 * (dt * 1932.0/2197.0) + k2 * (-dt * 7200.0/2197.0) + k3 * (dt * 7296.0/2197.0));
    StateType k5 = derivatives(t + dt, state + k1 * (dt * 439.0/216.0) + k2 * (-dt * 8.0) + k3 * (dt * 3680.0/513.0) + k4 * (-dt * 845.0/4104.0));
    StateType k6 = derivatives(t + dt * 0.5, state + k1 * (-dt * 8.0/27.0) + k2 * (dt * 2.0) + k3 * (-dt * 3544.0/2565.0) + k4 * (dt * 1859.0/4104.0) + k5 * (-dt * 11.0/40.0));
    
    // We use the 5th-order estimate to advance the solution
    return state + (k1 * (16.0/135.0) + k3 * (6656.0/12825.0) + k4 * (28561.0/56430.0) + k5 * (-9.0/50.0) + k6 * (2.0/55.0)) * dt;
}

} // namespace integrator
