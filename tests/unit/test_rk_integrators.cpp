#include "tests/test_framework.h"

#include <cmath>
#include <functional>
#include <string>
#include <vector>

#include "modules/rk_integrators.h"

namespace {

integrator::State integrate(
    const std::string& method,
    double dt,
    double t_end,
    integrator::State initial,
    const std::function<integrator::State(double, const integrator::State&)>& derivs) {
    int steps = static_cast<int>(std::round(t_end / dt));
    double t = 0.0;
    integrator::State state = initial;
    for (int i = 0; i < steps; ++i) {
        if (method == "rk3") {
            state = integrator::rk3_step(t, dt, state, derivs);
        } else if (method == "rk5") {
            state = integrator::rk5_step(t, dt, state, derivs);
        } else if (method == "semi_implicit_euler") {
            state = integrator::semi_implicit_euler_step(t, dt, state, derivs);
        } else if (method == "leapfrog") {
            state = integrator::leapfrog_step(t, dt, state, derivs);
        } else if (method == "ruth4") {
            state = integrator::ruth4_step(t, dt, state, derivs);
        } else if (method == "rk23") {
            state = integrator::rk23_step(t, dt, state, derivs);
        } else if (method == "rkf45") {
            state = integrator::rkf45_step(t, dt, state, derivs);
        } else {
            state = integrator::rk4_step(t, dt, state, derivs);
        }
        t += dt;
    }
    return state;
}

double observed_order_from_halving(const std::vector<double>& errors) {
    std::vector<double> local_orders;
    for (size_t i = 1; i < errors.size(); ++i) {
        if (errors[i - 1] <= 0.0 || errors[i] <= 0.0) {
            continue;
        }
        const double p = std::log(errors[i - 1] / errors[i]) / std::log(2.0);
        local_orders.push_back(p);
    }

    double sum = 0.0;
    for (double p : local_orders) {
        sum += p;
    }
    return local_orders.empty() ? 0.0 : sum / static_cast<double>(local_orders.size());
}

integrator::State integrate_position_only(
    const std::string& method,
    double dt,
    double t_end,
    integrator::State initial,
    const std::function<double(double, double)>& acceleration) {
    const int steps = static_cast<int>(std::round(t_end / dt));
    const auto states = integrator::integrate_position_only_trajectory(
        method, 0.0, dt, steps, initial, acceleration);
    return states.back();
}

integrator::State integrate_den(
    double dt,
    double t_end,
    double gamma,
    double omega0,
    integrator::State initial,
    const std::function<double(double, const integrator::State&)>& residual) {
    const int steps = static_cast<int>(std::round(t_end / dt));
    double t = 0.0;
    integrator::State state = initial;
    for (int i = 0; i < steps; ++i) {
        state = integrator::den3_step(t, dt, state, gamma, omega0, residual);
        t += dt;
    }
    return state;
}

}  // namespace

TEST(RkIntegratorsStateArithmetic) {
    const integrator::State a{1.0, -2.0};
    const integrator::State b{3.5, 4.0};

    const integrator::State c = a + b;
    EXPECT_NEAR(c.theta, 4.5, 1e-12);
    EXPECT_NEAR(c.omega, 2.0, 1e-12);

    const integrator::State d = a * 2.0;
    EXPECT_NEAR(d.theta, 2.0, 1e-12);
    EXPECT_NEAR(d.omega, -4.0, 1e-12);
}

TEST(RkIntegratorsConstantDerivativeIsExact) {
    auto derivs = [](double, const integrator::State&) -> integrator::State {
        return {2.0, -3.0};
    };

    const integrator::State initial{1.0, 5.0};
    const double dt = 0.1;

    const auto s3 = integrator::rk3_step(0.0, dt, initial, derivs);
    const auto s4 = integrator::rk4_step(0.0, dt, initial, derivs);
    const auto s5 = integrator::rk5_step(0.0, dt, initial, derivs);
    const auto s23 = integrator::rk23_step(0.0, dt, initial, derivs);
    const auto s45 = integrator::rkf45_step(0.0, dt, initial, derivs);

    const double theta_expected = 1.0 + 2.0 * dt;
    const double omega_expected = 5.0 - 3.0 * dt;

    EXPECT_NEAR(s3.theta, theta_expected, 1e-12);
    EXPECT_NEAR(s3.omega, omega_expected, 1e-12);
    EXPECT_NEAR(s4.theta, theta_expected, 1e-12);
    EXPECT_NEAR(s4.omega, omega_expected, 1e-12);
    EXPECT_NEAR(s5.theta, theta_expected, 1e-12);
    EXPECT_NEAR(s5.omega, omega_expected, 1e-12);
    EXPECT_NEAR(s23.theta, theta_expected, 1e-12);
    EXPECT_NEAR(s23.omega, omega_expected, 1e-12);
    EXPECT_NEAR(s45.theta, theta_expected, 1e-12);
    EXPECT_NEAR(s45.omega, omega_expected, 1e-12);
}

TEST(PositionOnlyIntegratorsConstantAccelerationIsExact) {
    auto acceleration = [](double, double) { return -3.0; };

    const integrator::State initial{1.0, 5.0};
    const double dt = 0.1;
    const double t_end = 0.5;

    const auto verlet = integrate_position_only("velocity_verlet", dt, t_end, initial, acceleration);
    const auto rkn = integrate_position_only("runge_kutta_nystrom", dt, t_end, initial, acceleration);
    const auto numerov = integrate_position_only("numerov", dt, t_end, initial, acceleration);

    const double theta_expected = initial.theta + initial.omega * t_end - 1.5 * t_end * t_end;
    const double omega_expected = initial.omega - 3.0 * t_end;

    EXPECT_NEAR(verlet.theta, theta_expected, 1e-12);
    EXPECT_NEAR(verlet.omega, omega_expected, 1e-12);
    EXPECT_NEAR(rkn.theta, theta_expected, 1e-12);
    EXPECT_NEAR(rkn.omega, omega_expected, 1e-12);
    EXPECT_NEAR(numerov.theta, theta_expected, 1e-12);
    EXPECT_NEAR(numerov.omega, omega_expected, 1e-12);
}

TEST(DenIntegratorExactForHomogeneousDampedOscillator) {
    const double gamma = 0.2;
    const double omega0 = 1.5;
    const double dt = 0.1;
    const double t_end = 2.0;
    const integrator::State initial{0.7, -0.15};
    auto residual = [](double, const integrator::State&) { return 0.0; };

    const integrator::State out = integrate_den(dt, t_end, gamma, omega0, initial, residual);
    const integrator::Propagator phi = integrator::den_propagator(t_end, gamma, omega0);
    const double theta_exact = phi.phi11 * initial.theta + phi.phi12 * initial.omega;
    const double omega_exact = phi.phi21 * initial.theta + phi.phi22 * initial.omega;

    EXPECT_NEAR(out.theta, theta_exact, 1e-12);
    EXPECT_NEAR(out.omega, omega_exact, 1e-12);
}

TEST(DenIntegratorObservedOrderMatchesFourthOrder) {
    const double gamma = 0.15;
    const double omega0 = 1.2;
    const double t_end = 1.0;
    const integrator::State initial{0.4, -0.3};
    auto residual = [](double t, const integrator::State&) {
        return std::sin(2.0 * t);
    };

    auto error_at = [&](double dt) {
        const integrator::State out = integrate_den(dt, t_end, gamma, omega0, initial, residual);
        const integrator::State ref = integrate_den(dt / 32.0, t_end, gamma, omega0, initial, residual);
        return std::max(std::fabs(out.theta - ref.theta), std::fabs(out.omega - ref.omega));
    };

    const std::vector<double> dts = {0.2, 0.1, 0.05, 0.025};
    std::vector<double> errors;
    for (double dt : dts) {
        errors.push_back(error_at(dt));
    }

    const double order = observed_order_from_halving(errors);
    EXPECT_TRUE(order > 3.5 && order < 4.5);
}

TEST(RkIntegratorsConvergenceOnExponentialDecay) {
    auto derivs = [](double, const integrator::State& s) -> integrator::State {
        return {-s.theta, -s.omega};
    };

    const integrator::State initial{1.0, -2.0};
    const double exact_theta = std::exp(-1.0) * initial.theta;
    const double exact_omega = std::exp(-1.0) * initial.omega;

    auto error_at = [&](const std::string& method, double dt) {
        const auto out = integrate(method, dt, 1.0, initial, derivs);
        const double e_theta = std::fabs(out.theta - exact_theta);
        const double e_omega = std::fabs(out.omega - exact_omega);
        return std::max(e_theta, e_omega);
    };

    const double e3_dt = error_at("rk3", 0.1);
    const double e3_half = error_at("rk3", 0.05);
    const double e4_dt = error_at("rk4", 0.1);
    const double e4_half = error_at("rk4", 0.05);
    const double e5_dt = error_at("rk5", 0.1);
    const double e5_half = error_at("rk5", 0.05);
    const double e23_dt = error_at("rk23", 0.1);
    const double e23_half = error_at("rk23", 0.05);
    const double e45_dt = error_at("rkf45", 0.1);
    const double e45_half = error_at("rkf45", 0.05);

    EXPECT_TRUE(e3_half < e3_dt);
    EXPECT_TRUE(e4_half < e4_dt);
    EXPECT_TRUE(e5_half < e5_dt);
    EXPECT_TRUE(e23_half < e23_dt);
    EXPECT_TRUE(e45_half < e45_dt);

    EXPECT_TRUE(e3_dt / e3_half > 4.0);
    EXPECT_TRUE(e4_dt / e4_half > 8.0);
    EXPECT_TRUE(e5_dt / e5_half > 12.0);
    EXPECT_TRUE(e23_dt / e23_half > 4.0);
    EXPECT_TRUE(e45_dt / e45_half > 12.0);
}

TEST(RkIntegratorsObservedOrderMatchesMethodOrder) {
    auto derivs = [](double, const integrator::State& s) -> integrator::State {
        return {-2.0 * s.theta, -2.0 * s.omega};
    };

    const integrator::State initial{1.0, -0.5};
    const double t_end = 1.0;
    const double exact_scale = std::exp(-2.0 * t_end);

    auto max_error = [&](const std::string& method, double dt) {
        const integrator::State out = integrate(method, dt, t_end, initial, derivs);
        const double theta_exact = initial.theta * exact_scale;
        const double omega_exact = initial.omega * exact_scale;
        const double err_theta = std::fabs(out.theta - theta_exact);
        const double err_omega = std::fabs(out.omega - omega_exact);
        return std::max(err_theta, err_omega);
    };

    const std::vector<double> dts = {0.2, 0.1, 0.05, 0.025};
    std::vector<double> e3;
    std::vector<double> e4;
    std::vector<double> e5;
    std::vector<double> e23;
    std::vector<double> e45;
    for (double dt : dts) {
        e3.push_back(max_error("rk3", dt));
        e4.push_back(max_error("rk4", dt));
        e5.push_back(max_error("rk5", dt));
        e23.push_back(max_error("rk23", dt));
        e45.push_back(max_error("rkf45", dt));
    }

    const double p3 = observed_order_from_halving(e3);
    const double p4 = observed_order_from_halving(e4);
    const double p5 = observed_order_from_halving(e5);
    const double p23 = observed_order_from_halving(e23);
    const double p45 = observed_order_from_halving(e45);

    EXPECT_TRUE(p3 > 2.6 && p3 < 3.4);
    EXPECT_TRUE(p4 > 3.6 && p4 < 4.4);
    EXPECT_TRUE(p5 > 4.4 && p5 < 5.6);
    EXPECT_TRUE(p23 > 2.6 && p23 < 3.4); // RK23 is 3rd order step
    EXPECT_TRUE(p45 > 4.4 && p45 < 5.6); // RKF45 is 5th order step

    EXPECT_TRUE(e3.back() < e3.front());
    EXPECT_TRUE(e4.back() < e4.front());
    EXPECT_TRUE(e5.back() < e5.front());
    EXPECT_TRUE(e23.back() < e23.front());
    EXPECT_TRUE(e45.back() < e45.front());
}

TEST(RkIntegratorsSymplecticHarmonicOscillator) {
    auto derivs = [](double, const integrator::State& s) -> integrator::State {
        return {s.omega, -s.theta};
    };

    const integrator::State initial{1.0, 0.0};
    const double t_end = 2.0 * M_PI; // One full period
    auto error_at = [&](const std::string& method, double dt) {
        int steps = std::round(t_end / dt);
        double t_actual_end = steps * dt;
        const double exact_theta = initial.theta * std::cos(t_actual_end) + initial.omega * std::sin(t_actual_end);
        const double exact_omega = -initial.theta * std::sin(t_actual_end) + initial.omega * std::cos(t_actual_end);
        const auto out = integrate(method, dt, t_actual_end, initial, derivs);
        const double e_theta = std::fabs(out.theta - exact_theta);
        const double e_omega = std::fabs(out.omega - exact_omega);
        return std::max(e_theta, e_omega);
    };

    const double e_euler_dt = error_at("semi_implicit_euler", 0.01);
    const double e_euler_half = error_at("semi_implicit_euler", 0.005);
    
    const double e_leapfrog_dt = error_at("leapfrog", 0.01);
    const double e_leapfrog_half = error_at("leapfrog", 0.005);
    
    const double e_ruth_dt = error_at("ruth4", 0.1);
    const double e_ruth_half = error_at("ruth4", 0.05);

    EXPECT_TRUE(e_euler_half < e_euler_dt);
    EXPECT_TRUE(e_leapfrog_half < e_leapfrog_dt);
    EXPECT_TRUE(e_ruth_half < e_ruth_dt);

    // Euler is 1st order (~2x better)
    EXPECT_TRUE(e_euler_dt / e_euler_half > 1.8);
    // Leapfrog is 2nd order (~4x better)
    EXPECT_TRUE(e_leapfrog_dt / e_leapfrog_half > 3.8);
    // Ruth is 4th order (~16x better)
    EXPECT_TRUE(e_ruth_dt / e_ruth_half > 14.0);
}

TEST(PositionOnlyIntegratorsObservedOrderMatchesMethodOrder) {
    auto acceleration = [](double, double theta) {
        return -theta;
    };

    const integrator::State initial{1.0, -0.25};
    const double t_end = 1.0;
    const std::vector<double> dts = {0.2, 0.1, 0.05, 0.025};

    auto error_at = [&](const std::string& method, double dt) {
        const integrator::State out =
            integrate_position_only(method, dt, t_end, initial, acceleration);
        const double exact_theta =
            initial.theta * std::cos(t_end) + initial.omega * std::sin(t_end);
        const double exact_omega =
            -initial.theta * std::sin(t_end) + initial.omega * std::cos(t_end);
        return std::max(std::fabs(out.theta - exact_theta),
                        std::fabs(out.omega - exact_omega));
    };

    std::vector<double> verlet_errors;
    std::vector<double> rkn_errors;
    std::vector<double> numerov_errors;
    for (double dt : dts) {
        verlet_errors.push_back(error_at("velocity_verlet", dt));
        rkn_errors.push_back(error_at("runge_kutta_nystrom", dt));
        numerov_errors.push_back(error_at("numerov", dt));
    }

    const double verlet_order = observed_order_from_halving(verlet_errors);
    const double rkn_order = observed_order_from_halving(rkn_errors);
    const double numerov_order = observed_order_from_halving(numerov_errors);

    EXPECT_TRUE(verlet_order > 1.8 && verlet_order < 2.2);
    EXPECT_TRUE(rkn_order > 3.6 && rkn_order < 4.4);
    EXPECT_TRUE(numerov_order > 3.4 && numerov_order < 4.4);

    EXPECT_TRUE(verlet_errors.back() < verlet_errors.front());
    EXPECT_TRUE(rkn_errors.back() < rkn_errors.front());
    EXPECT_TRUE(numerov_errors.back() < numerov_errors.front());
}
