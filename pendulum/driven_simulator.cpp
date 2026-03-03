#include "pendulum/driven_simulator.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace {
    struct State {
        double theta;
        double omega;
        State operator+(const State& other) const { return {theta + other.theta, omega + other.omega}; }
        State operator*(double scalar) const { return {theta * scalar, omega * scalar}; }
    };

    State derivatives(double t, const State& state, const DrivenPhysicalConfig& p) {
        State d_state;
        d_state.theta = state.omega;
        d_state.omega = -p.damping * state.omega
                        - (p.g / p.L) * std::sin(state.theta)
                        + p.A * std::cos(p.omega_drive * t);
        return d_state;
    }

    State rk4_step(double t, double dt, const State& state, const DrivenPhysicalConfig& p) {
        State k1 = derivatives(t, state, p);
        State k2 = derivatives(t + dt * 0.5, state + k1 * (dt * 0.5), p);
        State k3 = derivatives(t + dt * 0.5, state + k2 * (dt * 0.5), p);
        State k4 = derivatives(t + dt, state + k3 * dt, p);
        return state + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt / 6.0);
    }

    double analytical_solution(double t, const DrivenPhysicalConfig& p) {
        double omega_0 = std::sqrt(p.g / p.L);
        double denominator = std::sqrt(std::pow(omega_0 * omega_0 - p.omega_drive * p.omega_drive, 2) + std::pow(p.damping * p.omega_drive, 2));
        double R = p.A / denominator;
        double delta = std::atan2(p.damping * p.omega_drive, omega_0 * omega_0 - p.omega_drive * p.omega_drive);
        
        double theta_steady = R * std::cos(p.omega_drive * t - delta);
        
        double damping_half = p.damping / 2.0;
        double omega_trans_sq = omega_0 * omega_0 - damping_half * damping_half;
        if (omega_trans_sq <= 0) {
            // Overdamped or critically damped edge case (the original code assumes underdamped)
            return theta_steady; // simplified fallback
        }
        double omega_trans = std::sqrt(omega_trans_sq);
        
        double cos_delta = std::cos(delta);
        double sin_delta = std::sin(delta);
        
        double C1 = p.theta0 - R * cos_delta;
        double C2 = (p.omega0 + damping_half * C1 - R * p.omega_drive * sin_delta) / omega_trans;
        
        double theta_transient = C1 * std::exp(-damping_half * t) * std::cos(omega_trans * t) +
                                 C2 * std::exp(-damping_half * t) * std::sin(omega_trans * t);
                                 
        return theta_steady + theta_transient;
    }
}

DrivenPendulumSimulator::DrivenPendulumSimulator(const DrivenConfig& config)
    : config_(config) {}

DrivenSimulationResult DrivenPendulumSimulator::simulate() const {
    const auto& p = config_.physical;
    const auto& s = config_.simulation;

    DrivenSimulationResult result;
    State state = {p.theta0, p.omega0};
    
    double t = s.t_start;
    
    double max_diff = 0.0;
    double sum_diff = 0.0;
    int sample_count = 0;

    int nsteps = static_cast<int>((s.t_end - s.t_start) / s.dt + 0.5);
    result.rk4_steps = nsteps;

    for (int n = 0; n <= nsteps; ++n) {
        if (n % s.output_every == 0) {
            double theta_exact = analytical_solution(t, p);
            double diff = std::abs(state.theta - theta_exact);
            
            result.samples.push_back({
                t,
                state.theta,
                theta_exact,
                diff
            });
            
            max_diff = std::max(max_diff, diff);
            sum_diff += diff;
            sample_count++;
        }

        if (n < nsteps) {
            state = rk4_step(t, s.dt, state, p);
            t += s.dt;
        }
    }

    result.max_diff = max_diff;
    result.avg_diff = sample_count > 0 ? (sum_diff / sample_count) : 0.0;

    return result;
}
