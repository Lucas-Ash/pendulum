#include "modules/driven/driven_simulator.h"
#include "modules/integrators/simulation_runner.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

constexpr double kPi = 3.14159265358979323846;

struct ResolvedDrivenSystem {
    DrivenSystemModel model = DrivenSystemModel::Pendulum;
    double base_mass = 1.0;
    double linear_stiffness = 0.0;
    double cubic_stiffness = 0.0;
    double drive_force = 0.0;
    double omega_drive = 0.0;
    double pendulum_omega0_sq = 0.0;
    restoring_force::Config restoring_force;
};

double mass_scale_denominator(const DrivenUnitScalesConfig& scales) {
    return scales.stiffness_scale * scales.time_scale * scales.time_scale;
}

double force_scale_denominator(const DrivenUnitScalesConfig& scales) {
    return scales.stiffness_scale * scales.displacement_scale;
}

double cubic_stiffness_scale_factor(const DrivenUnitScalesConfig& scales) {
    return scales.displacement_scale * scales.displacement_scale /
           scales.stiffness_scale;
}

double cubic_damping_scale_factor(const DrivenUnitScalesConfig& scales) {
    return scales.displacement_scale * scales.displacement_scale /
           (scales.stiffness_scale * scales.time_scale * scales.time_scale *
            scales.time_scale);
}

DrivenConfig apply_unit_scales(const DrivenConfig& config) {
    if (!config.unit_scales.enabled) {
        return config;
    }

    if (config.physical.system_model == DrivenSystemModel::Pendulum) {
        throw std::runtime_error(
            "unit_scales are only supported for Duffing mode.");
    }

    const auto& scales = config.unit_scales;
    DrivenConfig scaled = config;

    scaled.physical.omega_drive *= scales.time_scale;
    scaled.physical.theta0 /= scales.displacement_scale;
    scaled.physical.omega0 *= scales.time_scale / scales.displacement_scale;
    scaled.physical.damping /= (scales.stiffness_scale * scales.time_scale);
    scaled.physical.damping_cubic *= cubic_damping_scale_factor(scales);

    scaled.physical.mass /= mass_scale_denominator(scales);
    scaled.physical.linear_stiffness /= scales.stiffness_scale;
    scaled.physical.cubic_stiffness *= cubic_stiffness_scale_factor(scales);
    scaled.physical.drive_force /= force_scale_denominator(scales);

    scaled.simulation.t_start /= scales.time_scale;
    scaled.simulation.t_end /= scales.time_scale;
    scaled.simulation.dt /= scales.time_scale;

    scaled.mass_event.jump_time /= scales.time_scale;
    scaled.mass_event.delta_mass /= mass_scale_denominator(scales);

    scaled.noise.force_stddev /= force_scale_denominator(scales);
    scaled.noise.correlation_time /= scales.time_scale;

    scaled.sweep.omega_start *= scales.time_scale;
    scaled.sweep.omega_end *= scales.time_scale;
    scaled.sweep.settle_time /= scales.time_scale;

    return scaled;
}

double splitmix64_unit(unsigned long long value) {
    unsigned long long z = value + 0x9e3779b97f4a7c15ull;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
    z ^= (z >> 31);
    return static_cast<double>((z >> 11) & ((1ull << 53) - 1ull)) /
           static_cast<double>(1ull << 53);
}

double gaussian_from_index(unsigned long long seed, unsigned long long index) {
    const double u1 = std::max(1e-12, splitmix64_unit(seed + 2ull * index));
    const double u2 = splitmix64_unit(seed + 2ull * index + 1ull);
    return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * kPi * u2);
}

ResolvedDrivenSystem resolve_system(const DrivenConfig& config) {
    ResolvedDrivenSystem system;
    const auto& p = config.physical;
    system.model = p.system_model;
    system.omega_drive = p.omega_drive;

    if (p.system_model == DrivenSystemModel::Pendulum) {
        system.base_mass = 1.0;
        system.pendulum_omega0_sq = p.g / p.L;
        system.linear_stiffness =
            system.pendulum_omega0_sq *
            restoring_force::linearized_slope(p.restoring_force);
        system.cubic_stiffness =
            p.restoring_force.model == restoring_force::Model::Polynomial
                ? system.pendulum_omega0_sq * p.restoring_force.cubic
                : 0.0;
        system.drive_force = p.A;
        system.restoring_force = p.restoring_force;
        return system;
    }

    if (p.system_model == DrivenSystemModel::Duffing) {
        system.base_mass = p.mass;
        system.linear_stiffness = p.linear_stiffness;
        system.cubic_stiffness = p.cubic_stiffness;
        system.drive_force = p.drive_force;
        system.restoring_force.model = restoring_force::Model::Polynomial;
        system.restoring_force.linear = p.linear_stiffness;
        system.restoring_force.cubic = p.cubic_stiffness;
        return system;
    }
    return system;
}

bool has_dynamic_mass(const DrivenConfig& config) {
    return config.mass_event.enabled &&
           config.physical.system_model != DrivenSystemModel::Pendulum &&
           std::abs(config.mass_event.delta_mass) > 0.0;
}

bool has_core_noise(const DrivenConfig& config) {
    return config.noise.enabled && config.noise.force_stddev > 0.0;
}

double damping_term(double theta, double omega, const DrivenPhysicalConfig& p) {
    return damping_force::term(
        theta, omega, p.damping_model, p.damping, p.damping_cubic);
}

void validate_position_only_integrator(const std::string& method,
                                       const DrivenPhysicalConfig& p,
                                       const DrivenConfig& config) {
    if (!integrator::is_position_only_integrator(method)) {
        return;
    }
    if (std::abs(p.damping) > 1e-12 || std::abs(p.damping_cubic) > 1e-12) {
        throw std::runtime_error(
            method + " requires theta'' = a(t, theta); it is not valid with velocity-dependent damping.");
    }
    if (has_core_noise(config)) {
        throw std::runtime_error(
            method + " is not supported with stochastic force noise enabled.");
    }
}

void validate_special_features(const DrivenConfig& config) {
    if (config.mass_event.enabled &&
        config.physical.system_model == DrivenSystemModel::Pendulum) {
        throw std::runtime_error(
            "Mass events are only supported for explicit Duffing mode.");
    }
    if (integrator::is_den_method(config.settings.integrator) &&
        (has_dynamic_mass(config) || has_core_noise(config))) {
        throw std::runtime_error(
            "den3 is not supported with dynamic mass events or stochastic force noise.");
    }
    if (config.unit_scales.enabled &&
        config.physical.system_model == DrivenSystemModel::Pendulum) {
        throw std::runtime_error(
            "unit_scales are only supported for Duffing mode.");
    }
}

double mass_at_time(double t,
                    const DrivenConfig& config,
                    const ResolvedDrivenSystem& system) {
    if (!has_dynamic_mass(config) || t < config.mass_event.jump_time) {
        return system.base_mass;
    }
    return system.base_mass + config.mass_event.delta_mass;
}

double drive_force_at_time(double t,
                           const DrivenConfig& config,
                           const ResolvedDrivenSystem& system) {
    if (config.mass_event.enabled &&
        config.mass_event.disable_drive_after_jump &&
        t >= config.mass_event.jump_time) {
        return 0.0;
    }
    return system.drive_force;
}

double noise_force_at_time(double t,
                           const DrivenConfig& config) {
    if (!has_core_noise(config)) {
        return 0.0;
    }
    const double correlation_time =
        config.noise.correlation_time > 0.0
            ? config.noise.correlation_time
            : config.simulation.dt;
    const double shifted_time = std::max(0.0, t - config.simulation.t_start);
    const auto bucket = static_cast<unsigned long long>(
        std::floor(shifted_time / correlation_time + 1e-12));
    return config.noise.force_stddev *
           gaussian_from_index(config.noise.seed, bucket);
}

double acceleration_of(double t,
                       double theta,
                       double omega,
                       const DrivenConfig& config,
                       const ResolvedDrivenSystem& system,
                       double drive_frequency) {
    const double mass = mass_at_time(t, config, system);
    const double drive_force = drive_force_at_time(t, config, system);
    const double drive_term = drive_force * std::cos(drive_frequency * t + config.physical.drive_phase);
    
    const double parametric_freq = config.physical.parametric_frequency == 0.0 ? 
                                   2.0 * drive_frequency : config.physical.parametric_frequency;
    const double parametric_stiffness = config.physical.parametric_amplitude * std::cos(parametric_freq * t);
    const double parametric_excitation_term = -parametric_stiffness * theta;

    const double noise_term = noise_force_at_time(t, config);

    if (system.model == DrivenSystemModel::Pendulum) {
        return -damping_term(theta, omega, config.physical) -
               system.pendulum_omega0_sq *
                   restoring_force::term(theta, system.restoring_force) +
               parametric_excitation_term +
               drive_term + noise_term;
    }

    return (-damping_term(theta, omega, config.physical) -
            system.linear_stiffness * theta -
            system.cubic_stiffness * theta * theta * theta +
            parametric_excitation_term +
            drive_term + noise_term) /
           mass;
}

double total_energy(double theta,
                    double omega,
                    const DrivenConfig& config,
                    const ResolvedDrivenSystem& system) {
    if (system.model == DrivenSystemModel::Pendulum) {
        return 0.5 * omega * omega +
               system.pendulum_omega0_sq *
                   restoring_force::potential(theta, system.restoring_force);
    }
    return 0.5 * system.base_mass * omega * omega +
           0.5 * system.linear_stiffness * theta * theta +
           0.25 * system.cubic_stiffness * theta * theta * theta * theta;
}

void analytical_solution(double t,
                         const DrivenPhysicalConfig& p,
                         double drive_frequency,
                         double omega_linearized_sq,
                         double& theta_exact,
                         double& omega_exact) {
    const double denominator = std::sqrt(
        std::pow(omega_linearized_sq - drive_frequency * drive_frequency, 2) +
        std::pow(p.damping * drive_frequency, 2));
    const double response = p.A / denominator;
    const double delta = std::atan2(
        p.damping * drive_frequency,
        omega_linearized_sq - drive_frequency * drive_frequency);

    const double theta_steady = response * std::cos(drive_frequency * t - delta);
    const double omega_steady = -response * drive_frequency *
                                std::sin(drive_frequency * t - delta);

    const double damping_half = p.damping / 2.0;
    const double omega_trans_sq = omega_linearized_sq - damping_half * damping_half;
    if (omega_trans_sq <= 0.0) {
        theta_exact = theta_steady;
        omega_exact = omega_steady;
        return;
    }

    const double omega_trans = std::sqrt(omega_trans_sq);
    const double cos_delta = std::cos(delta);
    const double sin_delta = std::sin(delta);
    const double c1 = p.theta0 - response * cos_delta;
    const double c2 =
        (p.omega0 + damping_half * c1 - response * drive_frequency * sin_delta) /
        omega_trans;

    const double exp_term = std::exp(-damping_half * t);
    const double cos_term = std::cos(omega_trans * t);
    const double sin_term = std::sin(omega_trans * t);
    const double theta_transient = exp_term * (c1 * cos_term + c2 * sin_term);
    const double omega_transient =
        -damping_half * theta_transient +
        exp_term * (-c1 * omega_trans * sin_term + c2 * omega_trans * cos_term);

    theta_exact = theta_steady + theta_transient;
    omega_exact = omega_steady + omega_transient;
}

SimulationResult simulate_window(const DrivenConfig& config,
                                 const ResolvedDrivenSystem& system,
                                 const integrator::State& initial_state,
                                 double t_start,
                                 double duration,
                                 double drive_frequency) {
    const auto& p = config.physical;
    const auto& settings = config.settings;
    const int nsteps = static_cast<int>(duration / config.simulation.dt + 0.5);
    const double omega_lin_sq = system.linear_stiffness;

    auto derivs = [&config, &system, drive_frequency](double t_val,
                                                      const integrator::State& state_val)
        -> integrator::State {
        return {
            state_val.omega,
            acceleration_of(t_val, state_val.theta, state_val.omega, config, system, drive_frequency)
        };
    };

    const double damping_over_mass =
        system.model == DrivenSystemModel::Pendulum
            ? p.damping
            : p.damping / system.base_mass;
    const double den_gamma = 0.5 * damping_over_mass;
    const double den_omega0 = std::sqrt(std::max(0.0, system.linear_stiffness / system.base_mass));
    auto den_residual =
        [&config, &system, drive_frequency, damping_over_mass](double t_val,
                                                               const integrator::State& state_val) {
            return acceleration_of(
                       t_val, state_val.theta, state_val.omega, config, system, drive_frequency) +
                   damping_over_mass * state_val.omega +
                   (system.linear_stiffness / system.base_mass) * state_val.theta;
        };

    auto acceleration = [&config, &system, drive_frequency](double t_val, double theta) {
        return acceleration_of(t_val, theta, 0.0, config, system, drive_frequency);
    };

    auto reference_state = [&settings, &p, omega_lin_sq, drive_frequency, &system](
                               double t_sample,
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
        if (system.model != DrivenSystemModel::Pendulum) {
            theta_reference = sample.theta;
            omega_reference = sample.omega;
            return;
        }
        analytical_solution(t_sample, p, drive_frequency, omega_lin_sq,
                            theta_reference, omega_reference);
    };

    auto energy_of = [&config, &system](const integrator::State& sample) {
        return total_energy(sample.theta, sample.omega, config, system);
    };

    return simulation_runner::run(
        {
            settings.integrator,
            t_start,
            config.simulation.dt,
            nsteps,
            settings.error_mode,
            settings.error_reference_factor,
            den_gamma,
            den_omega0,
        },
        initial_state,
        derivs,
        acceleration,
        den_residual,
        reference_state,
        energy_of);
}

double steady_state_amplitude(const SimulationResult& result,
                              double drive_frequency) {
    if (result.theta.empty() || result.t.empty()) {
        return 0.0;
    }
    const size_t start = static_cast<size_t>(0.65 * static_cast<double>(result.theta.size()));
    if (start >= result.theta.size()) {
        return 0.0;
    }

    if (std::abs(drive_frequency) < 1e-15) {
        double amplitude = 0.0;
        for (size_t i = start; i < result.theta.size(); ++i) {
            amplitude = std::max(amplitude, std::abs(result.theta[i]));
        }
        return amplitude;
    }

    size_t window_start = start;
    if (result.t.size() >= 2u) {
        const double dt = result.t[1] - result.t[0];
        const double period = 2.0 * kPi / std::abs(drive_frequency);
        if (dt > 0.0 && period > dt) {
            const size_t samples_per_period = static_cast<size_t>(
                std::max(1ll, std::llround(period / dt)));
            const size_t available_samples = result.theta.size() - start;
            const size_t full_periods = available_samples / samples_per_period;
            if (full_periods >= 1u) {
                const size_t measurement_count = full_periods * samples_per_period;
                window_start = result.theta.size() - measurement_count;
            }
        }
    }

    const size_t count = result.theta.size() - window_start;
    if (count == 0u) {
        return 0.0;
    }

    double mean = 0.0;
    for (size_t i = window_start; i < result.theta.size(); ++i) {
        mean += result.theta[i];
    }
    mean /= static_cast<double>(count);

    double cosine_projection = 0.0;
    double sine_projection = 0.0;
    for (size_t i = window_start; i < result.theta.size(); ++i) {
        const double centered = result.theta[i] - mean;
        const double phase = drive_frequency * result.t[i];
        cosine_projection += centered * std::cos(phase);
        sine_projection += centered * std::sin(phase);
    }

    const double scale = 2.0 / static_cast<double>(count);
    const double a = scale * cosine_projection;
    const double b = scale * sine_projection;
    return std::sqrt(a * a + b * b);
}

std::vector<double> harmonic_balance_roots(const DrivenConfig& config,
                                           const ResolvedDrivenSystem& system,
                                           double drive_frequency) {
    if (system.model == DrivenSystemModel::Pendulum ||
        config.physical.damping_model != damping_force::Model::Linear ||
        std::abs(config.physical.damping_cubic) > 1e-15) {
        return {};
    }

    const double mass = system.base_mass;
    const double c = config.physical.damping;
    const double k1 = system.linear_stiffness;
    const double k3 = system.cubic_stiffness;
    const double force = system.drive_force;
    const double delta = k1 - mass * drive_frequency * drive_frequency;

    if (std::abs(k3) < 1e-30) {
        const double denominator = std::sqrt(delta * delta + c * c * drive_frequency * drive_frequency);
        return {denominator == 0.0 ? 0.0 : force / denominator};
    }

    const double a3 = (9.0 / 16.0) * k3 * k3;
    const double a2 = 1.5 * k3 * delta;
    const double a1 = delta * delta + c * c * drive_frequency * drive_frequency;
    const double a0 = -force * force;

    const double b = a2 / a3;
    const double c2 = a1 / a3;
    const double d2 = a0 / a3;
    const double q = (3.0 * c2 - b * b) / 9.0;
    const double r = (9.0 * b * c2 - 27.0 * d2 - 2.0 * b * b * b) / 54.0;
    const double d = q * q * q + r * r;

    std::vector<double> roots_y;
    if (d >= 0.0) {
        const double s = std::cbrt(r + std::sqrt(d));
        const double t = std::cbrt(r - std::sqrt(d));
        roots_y.push_back(-b / 3.0 + s + t);
    } else {
        const double theta = std::acos(r / std::sqrt(-(q * q * q)));
        const double radius = 2.0 * std::sqrt(-q);
        for (int k = 0; k < 3; ++k) {
            roots_y.push_back(
                radius * std::cos((theta + 2.0 * kPi * static_cast<double>(k)) / 3.0) -
                b / 3.0);
        }
    }

    std::vector<double> amplitudes;
    for (double y : roots_y) {
        if (y >= 0.0) {
            amplitudes.push_back(std::sqrt(y));
        }
    }
    std::sort(amplitudes.begin(), amplitudes.end());
    amplitudes.erase(
        std::unique(amplitudes.begin(), amplitudes.end(),
                    [](double lhs, double rhs) { return std::abs(lhs - rhs) < 1e-18; }),
        amplitudes.end());
    return amplitudes;
}

double choose_stable_branch(const std::vector<double>& roots,
                            DrivenSweepDirection direction,
                            double cubic_stiffness) {
    if (roots.empty()) {
        return 0.0;
    }

    if (roots.size() < 3) {
        return roots.front();
    }

    const bool hardening = cubic_stiffness >= 0.0;
    const bool choose_lower =
        (hardening && direction == DrivenSweepDirection::Ascending) ||
        (!hardening && direction == DrivenSweepDirection::Descending);
    if (choose_lower) {
        return roots.front();
    }
    return roots.back();
}

}  // namespace

DrivenPendulumSimulator::DrivenPendulumSimulator(const DrivenConfig& config)
    : config_(config) {}

SimulationResult DrivenPendulumSimulator::simulate() const {
    const DrivenConfig scaled_config = apply_unit_scales(config_);
    validate_special_features(scaled_config);
    validate_position_only_integrator(
        scaled_config.settings.integrator, scaled_config.physical, scaled_config);
    const ResolvedDrivenSystem system = resolve_system(scaled_config);
    const integrator::State initial_state = {
        scaled_config.physical.theta0,
        scaled_config.physical.omega0,
    };
    const double duration =
        scaled_config.simulation.t_end - scaled_config.simulation.t_start;
    return simulate_window(
        scaled_config,
        system,
        initial_state,
        scaled_config.simulation.t_start,
        duration,
        scaled_config.physical.omega_drive);
}

DrivenSweepResult DrivenPendulumSimulator::simulate_sweep() const {
    const DrivenConfig scaled_config = apply_unit_scales(config_);
    validate_special_features(scaled_config);
    validate_position_only_integrator(
        scaled_config.settings.integrator, scaled_config.physical, scaled_config);
    const ResolvedDrivenSystem system = resolve_system(scaled_config);

    std::vector<double> frequencies(static_cast<size_t>(scaled_config.sweep.points), 0.0);
    const double span = scaled_config.sweep.omega_end - scaled_config.sweep.omega_start;
    for (int i = 0; i < scaled_config.sweep.points; ++i) {
        const double alpha =
            static_cast<double>(i) / static_cast<double>(scaled_config.sweep.points - 1);
        frequencies[static_cast<size_t>(i)] =
            scaled_config.sweep.omega_start + alpha * span;
    }
    if (scaled_config.sweep.direction == DrivenSweepDirection::Descending) {
        std::reverse(frequencies.begin(), frequencies.end());
    }

    integrator::State state = {scaled_config.physical.theta0, scaled_config.physical.omega0};
    double segment_start = scaled_config.simulation.t_start;
    DrivenSweepResult result;
    result.samples.reserve(frequencies.size());

    for (double drive_frequency : frequencies) {
        const SimulationResult segment = simulate_window(
            scaled_config,
            system,
            state,
            segment_start,
            scaled_config.sweep.settle_time,
            drive_frequency);

        DrivenSweepSample sample;
        sample.drive_frequency = drive_frequency;
        sample.numerical_amplitude = steady_state_amplitude(segment, drive_frequency);
        sample.final_theta = segment.theta.empty() ? 0.0 : segment.theta.back();
        sample.final_omega = segment.omega.empty() ? 0.0 : segment.omega.back();

        const std::vector<double> hb_roots =
            harmonic_balance_roots(scaled_config, system, drive_frequency);
        if (!hb_roots.empty()) {
            sample.analytical_lower_stable_amplitude = hb_roots.front();
            sample.analytical_upper_stable_amplitude =
                hb_roots.size() >= 3 ? hb_roots.back() : hb_roots.front();
            sample.analytical_amplitude =
                scaled_config.sweep.analytical_branch_tracking
                    ? choose_stable_branch(
                          hb_roots,
                          scaled_config.sweep.direction,
                          system.cubic_stiffness)
                    : hb_roots.front();
        }

        result.samples.push_back(sample);

        if (scaled_config.sweep.reuse_final_state) {
            state.theta = sample.final_theta;
            state.omega = sample.final_omega;
        } else {
            state.theta = scaled_config.physical.theta0;
            state.omega = scaled_config.physical.omega0;
        }
        segment_start += scaled_config.sweep.settle_time;
    }

    return result;
}
