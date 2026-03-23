#include "modules/core/error_analysis.h"
#include <cmath>
#include <algorithm>

void compute_error_statistics(SimulationResult& result) {
    double max_theta_error = 0.0;
    double max_omega_error = 0.0;
    double avg_theta_error = 0.0;
    double avg_omega_error = 0.0;
    double max_theta_rel_error = 0.0;
    double max_omega_rel_error = 0.0;
    double avg_theta_rel_error = 0.0;
    double avg_omega_rel_error = 0.0;

    result.theta_errors.clear();
    result.omega_errors.clear();
    result.theta_errors.reserve(result.t.size());
    result.omega_errors.reserve(result.t.size());

    for (size_t i = 0; i < result.t.size(); ++i) {
        const double theta_err = std::abs(result.theta[i] - result.theta_analytical[i]);
        const double omega_err = std::abs(result.omega[i] - result.omega_analytical[i]);
        
        const double theta_rel = theta_err / std::max(std::abs(result.theta_analytical[i]), 1e-12);
        const double omega_rel = omega_err / std::max(std::abs(result.omega_analytical[i]), 1e-12);

        result.theta_errors.push_back(theta_err);
        result.omega_errors.push_back(omega_err);

        max_theta_error = std::max(max_theta_error, theta_err);
        max_omega_error = std::max(max_omega_error, omega_err);
        avg_theta_error += theta_err;
        avg_omega_error += omega_err;
        max_theta_rel_error = std::max(max_theta_rel_error, theta_rel);
        max_omega_rel_error = std::max(max_omega_rel_error, omega_rel);
        avg_theta_rel_error += theta_rel;
        avg_omega_rel_error += omega_rel;
    }

    const double sample_count = static_cast<double>(result.t.size());
    if (sample_count > 0.0) {
        result.theta_stats = {
            max_theta_error,
            avg_theta_error / sample_count,
            max_theta_rel_error,
            avg_theta_rel_error / sample_count,
        };
        result.omega_stats = {
            max_omega_error,
            avg_omega_error / sample_count,
            max_omega_rel_error,
            avg_omega_rel_error / sample_count,
        };
    }
}
