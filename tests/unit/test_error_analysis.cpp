#include "tests/test_framework.h"

#include "modules/error_analysis.h"
#include "modules/simulation_result.h"

TEST(ErrorAnalysisEmptyResultDoesNotCrash) {
    SimulationResult result;
    compute_error_statistics(result);

    EXPECT_TRUE(result.theta_errors.empty());
    EXPECT_TRUE(result.omega_errors.empty());
    EXPECT_NEAR(result.theta_stats.max_abs, 0.0, 1e-15);
    EXPECT_NEAR(result.omega_stats.max_abs, 0.0, 1e-15);
}

TEST(ErrorAnalysisSinglePointZeroError) {
    SimulationResult result;
    result.t = {0.0};
    result.theta = {0.1};
    result.omega = {-0.2};
    result.theta_analytical = {0.1};
    result.omega_analytical = {-0.2};

    compute_error_statistics(result);

    EXPECT_EQ(result.theta_errors.size(), 1u);
    EXPECT_EQ(result.omega_errors.size(), 1u);
    EXPECT_NEAR(result.theta_errors[0], 0.0, 1e-15);
    EXPECT_NEAR(result.omega_errors[0], 0.0, 1e-15);
    EXPECT_NEAR(result.theta_stats.max_abs, 0.0, 1e-15);
    EXPECT_NEAR(result.theta_stats.avg_abs, 0.0, 1e-15);
    EXPECT_NEAR(result.omega_stats.max_abs, 0.0, 1e-15);
    EXPECT_NEAR(result.omega_stats.avg_abs, 0.0, 1e-15);
}

TEST(ErrorAnalysisManualErrorsMatchExpectedStats) {
    SimulationResult result;
    result.t = {0.0, 0.5, 1.0};
    result.theta = {1.0, 0.8, 0.6};
    result.theta_analytical = {1.0, 1.0, 1.0};
    result.omega = {0.0, 0.2, 0.4};
    result.omega_analytical = {0.0, 0.1, 0.3};

    compute_error_statistics(result);

    EXPECT_NEAR(result.theta_stats.max_abs, 0.4, 1e-12);
    EXPECT_NEAR(result.theta_stats.avg_abs, (0.0 + 0.2 + 0.4) / 3.0, 1e-12);
    EXPECT_NEAR(result.theta_stats.max_rel, 0.4, 1e-12);
    EXPECT_NEAR(result.theta_stats.avg_rel, (0.0 + 0.2 + 0.4) / 3.0, 1e-12);

    EXPECT_NEAR(result.omega_stats.max_abs, 0.1, 1e-12);
    EXPECT_NEAR(result.omega_stats.avg_abs, (0.0 + 0.1 + 0.1) / 3.0, 1e-12);
    EXPECT_NEAR(result.omega_stats.max_rel, 1.0, 1e-12);
    EXPECT_NEAR(result.omega_stats.avg_rel, (0.0 + 1.0 + (0.1 / 0.3)) / 3.0, 1e-12);
}

TEST(ErrorAnalysisZeroAnalyticalUsesGuardDenominator) {
    SimulationResult result;
    result.t = {0.0};
    result.theta = {1e-9};
    result.theta_analytical = {0.0};
    result.omega = {-2e-9};
    result.omega_analytical = {0.0};

    compute_error_statistics(result);

    EXPECT_NEAR(result.theta_stats.max_rel, 1e3, 1e-9);
    EXPECT_NEAR(result.omega_stats.max_rel, 2e3, 1e-9);
}
