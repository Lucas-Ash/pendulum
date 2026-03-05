#include "modules/jacobi_elliptic.h"
#include <cmath>
#include <vector>

namespace math_utils {

void jacobi_sn_cn_dn(double u, double m, double& sn, double& cn, double& dn) {
    if (m < 0.0 || m > 1.0) {
        // Fallbacks for out-of-bounds (though physical pendulums should have m in [0, 1])
        if (m < 0.0) {
            double mu = -m / (1.0 - m);
            double sn_mu, cn_mu, dn_mu;
            jacobi_sn_cn_dn(u * std::sqrt(1.0 - m), mu, sn_mu, cn_mu, dn_mu);
            sn = sn_mu / dn_mu;
            cn = cn_mu / dn_mu;
            dn = 1.0 / dn_mu;
            return;
        }
        if (m > 1.0) {
            double mu = 1.0 / m;
            double sn_mu, cn_mu, dn_mu;
            jacobi_sn_cn_dn(u * std::sqrt(m), mu, sn_mu, cn_mu, dn_mu);
            sn = sn_mu * std::sqrt(mu);
            cn = dn_mu;
            dn = cn_mu;
            return;
        }
    }

    if (m == 0.0) {
        sn = std::sin(u);
        cn = std::cos(u);
        dn = 1.0;
        return;
    }

    if (m == 1.0) {
        sn = std::tanh(u);
        double sech = 1.0 / std::cosh(u);
        cn = sech;
        dn = sech;
        return;
    }

    const int max_iter = 10;
    std::vector<double> a(max_iter);
    std::vector<double> c(max_iter);

    a[0] = 1.0;
    c[0] = std::sqrt(m);

    int n = 0;
    while (n < max_iter - 1) {
        if (c[n] < 1e-15) {
            break;
        }
        double b = std::sqrt(a[n] * a[n] - c[n] * c[n]);
        a[n + 1] = 0.5 * (a[n] + b);
        c[n + 1] = 0.5 * (a[n] - b);
        n++;
    }

    double phi = std::pow(2.0, n) * a[n] * u;
    for (int i = n; i > 0; --i) {
        phi = 0.5 * (std::asin(c[i] * std::sin(phi) / a[i]) + phi);
    }

    sn = std::sin(phi);
    cn = std::cos(phi);
    dn = std::sqrt(1.0 - m * sn * sn);
}

} // namespace math_utils
