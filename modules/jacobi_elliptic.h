#pragma once

namespace math_utils {

// Computes the Jacobi elliptic functions sn(u, m), cn(u, m), and dn(u, m)
// where m is the parameter (m = k^2), NOT the elliptic modulus k.
// Uses the Arithmetic-Geometric Mean (AGM) / descending Landen transformation.
void jacobi_sn_cn_dn(double u, double m, double& sn, double& cn, double& dn);

} // namespace math_utils
