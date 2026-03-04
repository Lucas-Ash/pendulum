#include "modules/driven_simulator.h"
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

    void analytical_solution(double t, const DrivenPhysicalConfig& p, double& theta_exact, double& omega_exact) {
        double omega_0 = std::sqrt(p.g / p.L);
        double denominator = std::sqrt(std::pow(omega_0 * omega_0 - p.omega_drive * p.omega_drive, 2) + std::pow(p.damping * p.omega_drive, 2));
        double R = p.A / denominator;
        double delta = std::atan2(p.damping * p.omega_drive, omega_0 * omega_0 - p.omega_drive * p.omega_drive);
        
        double theta_steady = R * std::cos(p.omega_drive * t - delta);
        double omega_steady = -R * p.omega_drive * std::sin(p.omega_drive * t - delta);
        
        double damping_half = p.damping / 2.0;
        double omega_trans_sq = omega_0 * omega_0 - damping_half * damping_half;
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

    double total_energy(double theta, double theta_dot, double omega0) {
        return 0.5 * theta_dot * theta_dot + omega0 * omega0 * (1.0 - std::cos(theta));
    }
}

DrivenPendulumSimulator::DrivenPendulumSimulator(const DrivenConfig& config)
    : config_(config) {}

SimulationResult DrivenPendulumSimulator::simulate() const {
    const auto& p = config_.physical;
    const auto& s = config_.simulation;

    SimulationResult result;
    State state = {p.theta0, p.omega0};
    double omega0 = std::sqrt(p.g / p.L);
    
    double t = s.t_start;
    
    double max_theta_error = 0.0;
    double max_omega_error = 0.0;
    double avg_theta_error = 0.0;
    double avg_omega_error = 0.0;
    double max_theta_rel_error = 0.0;
    double max_omega_rel_error = 0.0;
    double avg_theta_rel_error = 0.0;
    double avg_omega_rel_error = 0.0;

    int nsteps = static_cast<int>((s.t_end - s.t_start) / s.dt + 0.5);
    result.rk4_steps = nsteps;

    for (int n = 0; n <= nsteps; ++n) {
        double theta_exact = 0.0;
        double omega_exact = 0.0;
        analytical_solution(t, p, theta_exact, omega_exact);
        
        result.t.push_back(t);
        result.theta.push_back(state.theta);
        result.omega.push_back(state.omega);
        result.theta_analytical.push_back(theta_exact);
        result.omega_analytical.push_back(omega_exact);
        result.energy.push_back(total_energy(state.theta, state.omega, omega0));

        const double theta_err = std::abs(state.theta - theta_exact);
        const double omega_err = std::abs(state.omega - omega_exact);
        const double theta_rel = theta_err / std::max(std::abs(theta_exact), 1e-12);
        const double omega_rel = omega_err / std::max(std::abs(omega_exact), 1e-12);

        result.theta_errors.push_back(theta_err);
        result.omega_errors.push_back(omega_err);

        max_theta_error = std::max(max_theta_error, theta_err);
        max_omega_error = std::max(max_omega_error, omega_err);
        avg_theta_error += theta_err;
        avg_omega_error += omega_err;
        max_theta_rel_error = std::max(max_theta_rel_error, theta_rel);
        max_omega_rel_error = std::max(max_omega_rel_error, omega_rel);
        avg_theta_rel_error += theta_rel;
        avg_omega_rel_error += omega_rel;

        if (n < nsteps) {
            state = rk4_step(t, s.dt, state, p);
            t += s.dt;
        }
    }

    const double sample_count = static_cast<double>(result.t.size());
    result.theta_stats = {
        max_theta_error,
        avg_theta_error / sample_count,
        max_theta_rel_error,
        avg_theta_rel_error / sample_count,
    };
    result.omega_stats = {
        max_omega_error,
        avg_omega_error / sample_count,
        max_omega_rel_error,
        avg_omega_rel_error / sample_count,
    };

    return result;
}
