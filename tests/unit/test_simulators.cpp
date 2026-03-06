#include "tests/test_framework.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "modules/damped_simulator.h"
#include "modules/driven_simulator.h"
#include "modules/pendulum_simulator.h"

namespace {

double max_in_vector(const std::vector<double>& values) {
    double v = values.empty() ? 0.0 : values.front();
    for (double x : values) {
        v = std::max(v, x);
    }
    return v;
}

double segment_abs_max(const std::vector<double>& values, size_t start, size_t end) {
    double out = 0.0;
    for (size_t i = start; i < end; ++i) {
        out = std::max(out, std::fabs(values[i]));
    }
    return out;
}

double correlation_with_driver(
    const std::vector<double>& t,
    const std::vector<double>& theta,
    double omega_drive,
    size_t start) {
    double num = 0.0;
    double den_theta = 0.0;
    double den_drive = 0.0;
    for (size_t i = start; i < t.size(); ++i) {
        const double drive = std::cos(omega_drive * t[i]);
        num += theta[i] * drive;
        den_theta += theta[i] * theta[i];
        den_drive += drive * drive;
    }
    const double denom = std::sqrt(den_theta * den_drive);
    return denom > 1e-15 ? num / denom : 0.0;
}

}  // namespace

TEST(PendulumSimulatorLinearAnalyticalCloseForSmallAngle) {
    PendulumSimulator sim(1.0, 9.81, 0.001, 2.0);
    const SimulationResult out = sim.simulate(0.01, 0.0, "rk5", "linear");

    EXPECT_FALSE(out.t.empty());
    EXPECT_TRUE(out.theta.size() == out.theta_analytical.size());
    EXPECT_TRUE(out.theta_stats.max_abs < 1e-4);
    EXPECT_TRUE(out.omega_stats.max_abs < 1e-3);
}

TEST(PendulumSimulatorJacobiAnalyticalProducesFiniteTrajectories) {
    PendulumSimulator sim(1.0, 9.81, 0.001, 3.0);
    const SimulationResult out = sim.simulate(1.0, 0.0, "rk4", "jacobi");

    EXPECT_FALSE(out.t.empty());
    EXPECT_NEAR(out.theta.front(), out.theta_analytical.front(), 1e-12);
    EXPECT_NEAR(out.omega.front(), out.omega_analytical.front(), 1e-12);
    EXPECT_TRUE(out.theta_stats.max_abs < 0.2);
    for (double v : out.theta) {
        EXPECT_FINITE(v);
    }
}

TEST(PendulumSimulatorEnergyConservationWithoutDamping) {
    PendulumSimulator sim(1.0, 9.81, 0.001, 10.0);
    const SimulationResult out = sim.simulate(0.4, 0.0, "rk5", "jacobi");

    const double min_e = *std::min_element(out.energy.begin(), out.energy.end());
    const double max_e = *std::max_element(out.energy.begin(), out.energy.end());
    EXPECT_TRUE((max_e - min_e) < 1e-3);
}

TEST(PendulumSimulatorIntegratorChoicesAllSane) {
    PendulumSimulator sim(1.0, 9.81, 0.002, 4.0);

    const SimulationResult rk3 = sim.simulate(0.4, -0.1, "rk3", "jacobi");
    const SimulationResult rk4 = sim.simulate(0.4, -0.1, "rk4", "jacobi");
    const SimulationResult rk5 = sim.simulate(0.4, -0.1, "rk5", "jacobi");

    EXPECT_EQ(rk3.t.size(), rk4.t.size());
    EXPECT_EQ(rk4.t.size(), rk5.t.size());
    EXPECT_TRUE(std::fabs(rk4.theta.back() - rk5.theta.back()) < 5e-3);
    EXPECT_TRUE(std::fabs(rk3.theta.back() - rk5.theta.back()) < 2e-2);
}

TEST(PendulumSimulatorDuffingJacobiReferenceTracksNumericalSolution) {
    restoring_force::Config restoring;
    restoring.model = restoring_force::Model::Polynomial;
    restoring.linear = 1.0;
    restoring.cubic = 0.35;

    PendulumSimulator sim(1.0, 9.81, 0.0005, 6.0, restoring);
    const SimulationResult out = sim.simulate(0.45, 0.0, "rk5", "duffing_jacobi");

    EXPECT_FALSE(out.t.empty());
    EXPECT_NEAR(out.theta.front(), out.theta_analytical.front(), 1e-12);
    EXPECT_NEAR(out.omega.front(), out.omega_analytical.front(), 1e-12);
    EXPECT_TRUE(out.theta_stats.max_abs < 5e-3);
    EXPECT_TRUE(out.omega_stats.max_abs < 5e-2);

    const double min_e = *std::min_element(out.energy.begin(), out.energy.end());
    const double max_e = *std::max_element(out.energy.begin(), out.energy.end());
    EXPECT_TRUE((max_e - min_e) < 2e-3);
}

TEST(PendulumSimulatorCanDisableErrorAnalysis) {
    PendulumSimulator sim(1.0, 9.81, 0.002, 3.0);
    const SimulationResult out = sim.simulate(
        0.35, 0.0, "rk4", "jacobi", error_reference::Mode::None, 50);

    EXPECT_FALSE(out.t.empty());
    EXPECT_TRUE(out.theta_stats.max_abs < 1e-15);
    EXPECT_TRUE(out.omega_stats.max_abs < 1e-15);
}

TEST(PendulumSimulatorSupportsHdReferenceErrorAnalysis) {
    PendulumSimulator sim(1.0, 9.81, 0.01, 2.0);
    const SimulationResult out = sim.simulate(
        0.7, -0.1, "rk4", "jacobi", error_reference::Mode::HdReference, 50);

    EXPECT_FALSE(out.t.empty());
    EXPECT_TRUE(out.theta_stats.max_abs < 5e-2);
    EXPECT_TRUE(out.theta_stats.max_abs > 1e-10);
    EXPECT_TRUE(out.omega_stats.max_abs < 2e-1);
}

TEST(DampedSimulatorShowsDecayAndEnergyDissipation) {
    DampedConfig cfg;
    cfg.physical.gamma = 0.2;
    cfg.physical.theta0 = 0.6;
    cfg.physical.theta_dot0 = 0.0;
    cfg.simulation.t_start = 0.0;
    cfg.simulation.t_end = 10.0;
    cfg.simulation.dt = 0.002;
    cfg.settings.integrator = "rk5";

    DampedPendulumSimulator sim(cfg);
    const SimulationResult out = sim.simulate();

    EXPECT_FALSE(out.t.empty());
    const size_t quarter = out.theta.size() / 4;
    const double first_peak = segment_abs_max(out.theta, 0, quarter);
    const double last_peak = segment_abs_max(out.theta, out.theta.size() - quarter, out.theta.size());
    EXPECT_TRUE(last_peak < first_peak);

    EXPECT_TRUE(out.energy.back() < out.energy.front());
    EXPECT_TRUE(out.theta_stats.avg_abs < 0.1);
}

TEST(DampedSimulatorRejectsNonUnderdampedParameters) {
    DampedConfig cfg;
    cfg.physical.g = 9.81;
    cfg.physical.L = 1.0;
    cfg.physical.gamma = 4.0;
    cfg.simulation.t_end = 1.0;
    cfg.simulation.dt = 0.01;

    DampedPendulumSimulator sim(cfg);
    EXPECT_THROW(sim.simulate());
}

TEST(DampedSimulatorPolynomialRestoringRuns) {
    DampedConfig cfg;
    cfg.physical.gamma = 0.2;
    cfg.physical.theta0 = 0.4;
    cfg.physical.theta_dot0 = 0.0;
    cfg.physical.restoring_force.model = restoring_force::Model::Polynomial;
    cfg.physical.restoring_force.linear = 1.0;
    cfg.physical.restoring_force.cubic = 0.3;
    cfg.simulation.t_start = 0.0;
    cfg.simulation.t_end = 6.0;
    cfg.simulation.dt = 0.002;

    DampedPendulumSimulator sim(cfg);
    const SimulationResult out = sim.simulate();
    EXPECT_FALSE(out.t.empty());
    EXPECT_FINITE(max_in_vector(out.theta_errors));
}

TEST(DampedSimulatorHdReferenceErrorAnalysisRuns) {
    DampedConfig cfg;
    cfg.physical.gamma = 0.1;
    cfg.physical.theta0 = 0.5;
    cfg.simulation.t_end = 4.0;
    cfg.simulation.dt = 0.01;
    cfg.settings.error_mode = error_reference::Mode::HdReference;
    cfg.settings.error_reference_factor = 50;
    cfg.settings.integrator = "rk4";

    DampedPendulumSimulator sim(cfg);
    const SimulationResult out = sim.simulate();
    EXPECT_FALSE(out.t.empty());
    EXPECT_TRUE(out.theta_stats.max_abs < 5e-2);
}

TEST(DrivenSimulatorSteadyStateAndPhaseLocking) {
    DrivenConfig cfg;
    cfg.physical.g = 9.81;
    cfg.physical.L = 1.0;
    cfg.physical.damping = 0.4;
    cfg.physical.A = 0.5;
    cfg.physical.omega_drive = 1.2;
    cfg.physical.theta0 = 0.0;
    cfg.physical.omega0 = 0.0;
    cfg.simulation.t_start = 0.0;
    cfg.simulation.t_end = 25.0;
    cfg.simulation.dt = 0.002;
    cfg.settings.integrator = "rk5";

    DrivenPendulumSimulator sim(cfg);
    const SimulationResult out = sim.simulate();

    EXPECT_FALSE(out.t.empty());
    const size_t half = out.theta.size() / 2;
    const size_t quarter = out.theta.size() / 4;
    const double mid_amp = segment_abs_max(out.theta, half - quarter, half);
    const double late_amp = segment_abs_max(out.theta, out.theta.size() - quarter, out.theta.size());
    EXPECT_TRUE(late_amp > 0.0);
    EXPECT_TRUE(std::fabs(late_amp - mid_amp) < 0.2);

    const double corr = correlation_with_driver(out.t, out.theta, cfg.physical.omega_drive, half);
    EXPECT_TRUE(corr > 0.2);
}

TEST(DrivenSimulatorParameterSweepRunsAndStaysFinite) {
    std::vector<std::pair<double, double>> cases = {
        {0.2, 0.8},
        {0.5, 1.2},
        {1.0, 1.8},
    };

    for (const auto& [A, omega_drive] : cases) {
        DrivenConfig cfg;
        cfg.physical.A = A;
        cfg.physical.omega_drive = omega_drive;
        cfg.physical.damping = 0.3;
        cfg.simulation.t_end = 6.0;
        cfg.simulation.dt = 0.002;
        cfg.settings.integrator = "rk4";

        DrivenPendulumSimulator sim(cfg);
        const SimulationResult out = sim.simulate();
        EXPECT_FALSE(out.t.empty());
        EXPECT_FINITE(max_in_vector(out.theta_errors));
    }
}

TEST(DrivenSimulatorPolynomialRestoringRuns) {
    DrivenConfig cfg;
    cfg.physical.damping = 0.25;
    cfg.physical.A = 0.3;
    cfg.physical.omega_drive = 1.1;
    cfg.physical.restoring_force.model = restoring_force::Model::Polynomial;
    cfg.physical.restoring_force.linear = 1.0;
    cfg.physical.restoring_force.cubic = 0.2;
    cfg.simulation.t_end = 5.0;
    cfg.simulation.dt = 0.002;

    DrivenPendulumSimulator sim(cfg);
    const SimulationResult out = sim.simulate();
    EXPECT_FALSE(out.t.empty());
    EXPECT_FINITE(max_in_vector(out.theta_errors));
}

TEST(DrivenSimulatorCanDisableErrorAnalysis) {
    DrivenConfig cfg;
    cfg.physical.damping = 0.3;
    cfg.physical.A = 0.4;
    cfg.physical.omega_drive = 1.1;
    cfg.simulation.t_end = 3.0;
    cfg.simulation.dt = 0.01;
    cfg.settings.error_mode = error_reference::Mode::None;
    cfg.settings.integrator = "rk4";

    DrivenPendulumSimulator sim(cfg);
    const SimulationResult out = sim.simulate();
    EXPECT_FALSE(out.t.empty());
    EXPECT_TRUE(out.theta_stats.max_abs < 1e-15);
}
