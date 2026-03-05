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

    const double theta_expected = 1.0 + 2.0 * dt;
    const double omega_expected = 5.0 - 3.0 * dt;

    EXPECT_NEAR(s3.theta, theta_expected, 1e-12);
    EXPECT_NEAR(s3.omega, omega_expected, 1e-12);
    EXPECT_NEAR(s4.theta, theta_expected, 1e-12);
    EXPECT_NEAR(s4.omega, omega_expected, 1e-12);
    EXPECT_NEAR(s5.theta, theta_expected, 1e-12);
    EXPECT_NEAR(s5.omega, omega_expected, 1e-12);
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

    EXPECT_TRUE(e3_half < e3_dt);
    EXPECT_TRUE(e4_half < e4_dt);
    EXPECT_TRUE(e5_half < e5_dt);

    EXPECT_TRUE(e3_dt / e3_half > 4.0);
    EXPECT_TRUE(e4_dt / e4_half > 8.0);
    EXPECT_TRUE(e5_dt / e5_half > 12.0);
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
    for (double dt : dts) {
        e3.push_back(max_error("rk3", dt));
        e4.push_back(max_error("rk4", dt));
        e5.push_back(max_error("rk5", dt));
    }

    const double p3 = observed_order_from_halving(e3);
    const double p4 = observed_order_from_halving(e4);
    const double p5 = observed_order_from_halving(e5);

    EXPECT_TRUE(p3 > 2.6 && p3 < 3.4);
    EXPECT_TRUE(p4 > 3.6 && p4 < 4.4);
    EXPECT_TRUE(p5 > 4.4 && p5 < 5.6);

    EXPECT_TRUE(e3.back() < e3.front());
    EXPECT_TRUE(e4.back() < e4.front());
    EXPECT_TRUE(e5.back() < e5.front());
}
