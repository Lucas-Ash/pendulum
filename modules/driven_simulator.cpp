#include "modules/driven_simulator.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>

#include "modules/rk_integrators.h"
#include "modules/error_analysis.h"

namespace {

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
    integrator::State state = {p.theta0, p.omega0};
    double omega0 = std::sqrt(p.g / p.L);
    
    double t = s.t_start;
    
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

        if (n < nsteps) {
            auto derivs = [&p](double t_val, const integrator::State& s) -> integrator::State {
                return {s.omega, -p.damping * s.omega - (p.g / p.L) * std::sin(s.theta) + p.A * std::cos(p.omega_drive * t_val)};
            };
            
            if (config_.settings.integrator == "rk3") {
                state = integrator::rk3_step(t, s.dt, state, derivs);
            } else if (config_.settings.integrator == "rk5") {
                state = integrator::rk5_step(t, s.dt, state, derivs);
            } else if (config_.settings.integrator == "semi_implicit_euler") {
                state = integrator::semi_implicit_euler_step(t, s.dt, state, derivs);
            } else if (config_.settings.integrator == "leapfrog") {
                state = integrator::leapfrog_step(t, s.dt, state, derivs);
            } else if (config_.settings.integrator == "ruth4") {
                state = integrator::ruth4_step(t, s.dt, state, derivs);
            } else if (config_.settings.integrator == "rk23") {
                state = integrator::rk23_step(t, s.dt, state, derivs);
            } else if (config_.settings.integrator == "rkf45") {
                state = integrator::rkf45_step(t, s.dt, state, derivs);
            } else {
                state = integrator::rk4_step(t, s.dt, state, derivs);
            }
            
            t += s.dt;
        }
    }

    compute_error_statistics(result);

    return result;
}
