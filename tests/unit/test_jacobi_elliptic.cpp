#include "tests/test_framework.h"

#include <cmath>

#include "modules/jacobi_elliptic.h"

TEST(JacobiSpecialCaseMZeroMatchesTrig) {
    const double u = 1.2345;
    double sn = 0.0;
    double cn = 0.0;
    double dn = 0.0;
    math_utils::jacobi_sn_cn_dn(u, 0.0, sn, cn, dn);

    EXPECT_NEAR(sn, std::sin(u), 1e-12);
    EXPECT_NEAR(cn, std::cos(u), 1e-12);
    EXPECT_NEAR(dn, 1.0, 1e-12);
}

TEST(JacobiSpecialCaseMOneMatchesHyperbolic) {
    const double u = 0.87;
    double sn = 0.0;
    double cn = 0.0;
    double dn = 0.0;
    math_utils::jacobi_sn_cn_dn(u, 1.0, sn, cn, dn);

    const double sech = 1.0 / std::cosh(u);
    EXPECT_NEAR(sn, std::tanh(u), 1e-12);
    EXPECT_NEAR(cn, sech, 1e-12);
    EXPECT_NEAR(dn, sech, 1e-12);
}

TEST(JacobiBoundaryInsideIntervalStaysFiniteAndConsistent) {
    const double u = 0.5;
    for (double m : {1e-12, 1.0 - 1e-12}) {
        double sn = 0.0;
        double cn = 0.0;
        double dn = 0.0;
        math_utils::jacobi_sn_cn_dn(u, m, sn, cn, dn);

        EXPECT_FINITE(sn);
        EXPECT_FINITE(cn);
        EXPECT_FINITE(dn);
        EXPECT_NEAR(sn * sn + cn * cn, 1.0, 1e-8);
        EXPECT_NEAR(dn * dn + m * sn * sn, 1.0, 1e-8);
    }
}

TEST(JacobiKnownValuesAtZero) {
    double sn = 99.0;
    double cn = 99.0;
    double dn = 99.0;
    math_utils::jacobi_sn_cn_dn(0.0, 0.5, sn, cn, dn);

    EXPECT_NEAR(sn, 0.0, 1e-15);
    EXPECT_NEAR(cn, 1.0, 1e-15);
    EXPECT_NEAR(dn, 1.0, 1e-15);
}

TEST(JacobiOutOfBoundsMStillReturnsFiniteValues) {
    for (double m : {-0.5, 1.5}) {
        double sn = 0.0;
        double cn = 0.0;
        double dn = 0.0;
        math_utils::jacobi_sn_cn_dn(0.7, m, sn, cn, dn);
        EXPECT_FINITE(sn);
        EXPECT_FINITE(cn);
        EXPECT_FINITE(dn);
    }
}
