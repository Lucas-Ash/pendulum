#pragma once

#include <stdexcept>
#include <string>
#include <vector>

#include "modules/error_analysis.h"
#include "modules/error_reference.h"
#include "modules/rk_integrators.h"
#include "modules/simulation_result.h"

namespace simulation_runner {

struct RunOptions {
    std::string integrator_method;
    double t_start = 0.0;
    double dt = 0.0;
    int nsteps = 0;
    error_reference::Mode error_mode = error_reference::Mode::Analytical;
    int error_reference_factor = 50;
    double den_gamma = 0.0;
    double den_omega0 = 0.0;
};

inline void validate_options(const RunOptions& options) {
    if (options.dt <= 0.0) {
        throw std::runtime_error("dt must be > 0");
    }
    if (options.nsteps < 0) {
        throw std::runtime_error("nsteps must be non-negative");
    }
    if (options.error_reference_factor <= 0) {
        throw std::runtime_error("error_reference_factor must be > 0");
    }
    if (options.error_mode == error_reference::Mode::HdReference &&
        options.error_reference_factor < 2) {
        throw std::runtime_error(
            "hd_reference mode requires error_reference_factor >= 2");
    }
}

template <typename DerivFunc, typename ResidualFunc>
integrator::State advance_state(const std::string& method,
                                double t,
                                double dt,
                                const integrator::State& state,
                                DerivFunc derivs,
                                double den_gamma,
                                double den_omega0,
                                ResidualFunc den_residual) {
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

template <typename ReferenceFunc, typename EnergyFunc>
inline void append_sample(SimulationResult& result,
                          double t,
                          const integrator::State& sample,
                          const integrator::State* hd_reference_state,
                          ReferenceFunc reference_state,
                          EnergyFunc energy_of) {
    result.t.push_back(t);
    result.theta.push_back(sample.theta);
    result.omega.push_back(sample.omega);

    double theta_reference = 0.0;
    double omega_reference = 0.0;
    reference_state(t, sample, hd_reference_state, theta_reference, omega_reference);
    result.theta_analytical.push_back(theta_reference);
    result.omega_analytical.push_back(omega_reference);
    result.energy.push_back(energy_of(sample));
}

template <typename DerivFunc,
          typename AccelFunc,
          typename ResidualFunc,
          typename ReferenceFunc,
          typename EnergyFunc>
SimulationResult run(const RunOptions& options,
                     const integrator::State& initial_state,
                     DerivFunc derivs,
                     AccelFunc acceleration,
                     ResidualFunc den_residual,
                     ReferenceFunc reference_state,
                     EnergyFunc energy_of) {
    validate_options(options);

    SimulationResult result;
    result.rk4_steps = options.nsteps;

    if (integrator::is_position_only_integrator(options.integrator_method)) {
        const std::vector<integrator::State> trajectory =
            integrator::integrate_position_only_trajectory(
                options.integrator_method,
                options.t_start,
                options.dt,
                options.nsteps,
                initial_state,
                acceleration);

        std::vector<integrator::State> trajectory_ref;
        if (options.error_mode == error_reference::Mode::HdReference) {
            trajectory_ref = integrator::integrate_position_only_trajectory(
                options.integrator_method,
                options.t_start,
                options.dt / static_cast<double>(options.error_reference_factor),
                options.nsteps * options.error_reference_factor,
                initial_state,
                acceleration);
        }

        for (int n = 0; n <= options.nsteps; ++n) {
            const double t = options.t_start + static_cast<double>(n) * options.dt;
            const auto& sample = trajectory[static_cast<size_t>(n)];
            const integrator::State* hd_reference_state = nullptr;
            if (options.error_mode == error_reference::Mode::HdReference) {
                hd_reference_state =
                    &trajectory_ref[static_cast<size_t>(n * options.error_reference_factor)];
            }
            append_sample(
                result, t, sample, hd_reference_state, reference_state, energy_of);
        }

        compute_error_statistics(result);
        return result;
    }

    integrator::State state = initial_state;
    integrator::State state_ref = initial_state;
    const double dt_ref = options.dt / static_cast<double>(options.error_reference_factor);
    double t_ref = options.t_start;

    for (int n = 0; n <= options.nsteps; ++n) {
        const double t = options.t_start + static_cast<double>(n) * options.dt;
        const integrator::State* hd_reference_state = nullptr;
        if (options.error_mode == error_reference::Mode::HdReference) {
            hd_reference_state = &state_ref;
        }

        append_sample(result, t, state, hd_reference_state, reference_state, energy_of);

        if (n >= options.nsteps) {
            continue;
        }

        state = advance_state(
            options.integrator_method,
            t,
            options.dt,
            state,
            derivs,
            options.den_gamma,
            options.den_omega0,
            den_residual);

        if (options.error_mode != error_reference::Mode::HdReference) {
            continue;
        }

        for (int k = 0; k < options.error_reference_factor; ++k) {
            state_ref = advance_state(
                options.integrator_method,
                t_ref,
                dt_ref,
                state_ref,
                derivs,
                options.den_gamma,
                options.den_omega0,
                den_residual);
            t_ref += dt_ref;
        }
    }

    compute_error_statistics(result);
    return result;
}

}  // namespace simulation_runner
