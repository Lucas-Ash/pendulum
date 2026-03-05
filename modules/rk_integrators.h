#pragma once

#include <cmath>

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
template <typename DerivFunc>
State rk23_step(double t, double dt, const State& state, DerivFunc derivatives) {
    State k1 = derivatives(t, state);
    State k2 = derivatives(t + dt * 0.5, state + k1 * (dt * 0.5));
    State k3 = derivatives(t + dt * 0.75, state + k2 * (dt * 0.75));
    
    // We use the 3rd-order step to advance the solution
    return state + (k1 * (2.0/9.0) + k2 * (1.0/3.0) + k3 * (4.0/9.0)) * dt;
}

// Runge-Kutta-Fehlberg (RKF45) explicit step
template <typename DerivFunc>
State rkf45_step(double t, double dt, const State& state, DerivFunc derivatives) {
    State k1 = derivatives(t, state);
    State k2 = derivatives(t + dt * 0.25, state + k1 * (dt * 0.25));
    State k3 = derivatives(t + dt * (3.0/8.0), state + k1 * (dt * 3.0/32.0) + k2 * (dt * 9.0/32.0));
    State k4 = derivatives(t + dt * (12.0/13.0), state + k1 * (dt * 1932.0/2197.0) + k2 * (-dt * 7200.0/2197.0) + k3 * (dt * 7296.0/2197.0));
    State k5 = derivatives(t + dt, state + k1 * (dt * 439.0/216.0) + k2 * (-dt * 8.0) + k3 * (dt * 3680.0/513.0) + k4 * (-dt * 845.0/4104.0));
    State k6 = derivatives(t + dt * 0.5, state + k1 * (-dt * 8.0/27.0) + k2 * (dt * 2.0) + k3 * (-dt * 3544.0/2565.0) + k4 * (dt * 1859.0/4104.0) + k5 * (-dt * 11.0/40.0));
    
    // We use the 5th-order estimate to advance the solution
    return state + (k1 * (16.0/135.0) + k3 * (6656.0/12825.0) + k4 * (28561.0/56430.0) + k5 * (-9.0/50.0) + k6 * (2.0/55.0)) * dt;
}

} // namespace integrator
