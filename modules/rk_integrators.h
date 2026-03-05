#pragma once

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

// A generic Runge-Kutta 4th Order integrator step function.
// DerivFunc is a callable that takes (time t, const State& current_state) 
// and returns the State derivative.
template <typename DerivFunc>
State rk4_step(double t, double dt, const State& state, DerivFunc derivatives) {
    State k1 = derivatives(t, state);
    State k2 = derivatives(t + dt * 0.5, state + k1 * (dt * 0.5));
    State k3 = derivatives(t + dt * 0.5, state + k2 * (dt * 0.5));
    State k4 = derivatives(t + dt, state + k3 * dt);
    
    return state + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt / 6.0);
}

// A generic Runge-Kutta 3rd Order (Standard) integrator step function.
template <typename DerivFunc>
State rk3_step(double t, double dt, const State& state, DerivFunc derivatives) {
    State k1 = derivatives(t, state);
    State k2 = derivatives(t + dt * 0.5, state + k1 * (dt * 0.5));
    State k3 = derivatives(t + dt, state + k1 * (-dt) + k2 * (2.0 * dt));
    
    return state + (k1 + k2 * 4.0 + k3) * (dt / 6.0);
}

// A generic Runge-Kutta 5th Order (Butcher's) integrator step function.
template <typename DerivFunc>
State rk5_step(double t, double dt, const State& state, DerivFunc derivatives) {
    State k1 = derivatives(t, state);
    State k2 = derivatives(t + dt * 0.25, state + k1 * (dt * 0.25));
    State k3 = derivatives(t + dt * 0.25, state + k1 * (dt * 0.125) + k2 * (dt * 0.125));
    State k4 = derivatives(t + dt * 0.5, state + k2 * (-dt * 0.5) + k3 * dt);
    State k5 = derivatives(t + dt * 0.75, state + k1 * (dt * 3.0/16.0) + k4 * (dt * 9.0/16.0));
    State k6 = derivatives(t + dt, state + k1 * (-dt * 3.0/7.0) + k2 * (dt * 2.0/7.0) + 
                                   k3 * (dt * 12.0/7.0) + k4 * (-dt * 12.0/7.0) + k5 * (dt * 8.0/7.0));
    
    return state + (k1 * 7.0 + k3 * 32.0 + k4 * 12.0 + k5 * 32.0 + k6 * 7.0) * (dt / 90.0);
}

} // namespace integrator
