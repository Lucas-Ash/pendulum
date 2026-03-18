#pragma once

#include <algorithm>
#include <numeric>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace frequency_estimation {

inline double estimate_dominant_frequency_hz(const std::vector<double>& signal,
                                             double dt,
                                             double expected_hz,
                                             double relative_window = 0.2) {
    if (signal.size() < 8) {
        throw std::runtime_error("frequency estimation requires at least 8 samples");
    }
    if (dt <= 0.0) {
        throw std::runtime_error("frequency estimation requires dt > 0");
    }

    const double mean =
        std::accumulate(signal.begin(), signal.end(), 0.0) /
        static_cast<double>(signal.size());
    std::vector<double> centered(signal.size(), 0.0);
    for (size_t i = 0; i < signal.size(); ++i) {
        const double window = 0.5 - 0.5 * std::cos(
            2.0 * 3.14159265358979323846 * static_cast<double>(i) /
            static_cast<double>(signal.size() - 1));
        centered[i] = (signal[i] - mean) * window;
    }

    const size_t nfreq = signal.size() / 2 + 1;
    std::vector<double> magnitudes(nfreq, 0.0);
    std::vector<double> freqs(nfreq, 0.0);
    for (size_t k = 0; k < nfreq; ++k) {
        double real = 0.0;
        double imag = 0.0;
        for (size_t n = 0; n < centered.size(); ++n) {
            const double angle =
                -2.0 * 3.14159265358979323846 *
                static_cast<double>(k) * static_cast<double>(n) /
                static_cast<double>(centered.size());
            real += centered[n] * std::cos(angle);
            imag += centered[n] * std::sin(angle);
        }
        magnitudes[k] = std::sqrt(real * real + imag * imag);
        freqs[k] = static_cast<double>(k) / (dt * static_cast<double>(centered.size()));
    }

    const double low = std::max(0.0, expected_hz * (1.0 - relative_window));
    const double high = expected_hz * (1.0 + relative_window);
    size_t best_index = 0;
    double best_magnitude = -1.0;
    for (size_t k = 0; k < nfreq; ++k) {
        if (freqs[k] < low || freqs[k] > high) {
            continue;
        }
        if (magnitudes[k] > best_magnitude) {
            best_magnitude = magnitudes[k];
            best_index = k;
        }
    }
    if (best_magnitude < 0.0) {
        best_index = static_cast<size_t>(
            std::distance(magnitudes.begin(),
                          std::max_element(magnitudes.begin(), magnitudes.end())));
    }
    if (best_index == 0 || best_index + 1 >= nfreq) {
        return freqs[best_index];
    }

    const double y1 = std::log(std::max(magnitudes[best_index - 1], 1e-30));
    const double y2 = std::log(std::max(magnitudes[best_index], 1e-30));
    const double y3 = std::log(std::max(magnitudes[best_index + 1], 1e-30));
    const double denominator = y1 - 2.0 * y2 + y3;
    if (std::abs(denominator) < 1e-30) {
        return freqs[best_index];
    }
    const double offset = 0.5 * (y1 - y3) / denominator;
    const double bin_width = freqs[1] - freqs[0];
    return freqs[best_index] + offset * bin_width;
}

inline std::vector<double> slice(const std::vector<double>& values,
                                 size_t begin,
                                 size_t end) {
    if (begin > end || end > values.size()) {
        throw std::runtime_error("invalid slice bounds");
    }
    return std::vector<double>(values.begin() + static_cast<std::ptrdiff_t>(begin),
                               values.begin() + static_cast<std::ptrdiff_t>(end));
}

}  // namespace frequency_estimation
