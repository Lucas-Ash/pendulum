#pragma once

#include <vector>

struct CoupledState {
    double q1;
    double omega1;
    double q2;
    double omega2;

    CoupledState operator+(const CoupledState& other) const {
        return {q1 + other.q1, omega1 + other.omega1, q2 + other.q2, omega2 + other.omega2};
    }

    CoupledState operator*(double scalar) const {
        return {q1 * scalar, omega1 * scalar, q2 * scalar, omega2 * scalar};
    }
};

struct CoupledSimulationResult {
    std::vector<double> t;
    std::vector<double> q1;
    std::vector<double> omega1;
    std::vector<double> q2;
    std::vector<double> omega2;

    std::vector<double> q1_reference;
    std::vector<double> omega1_reference;
    std::vector<double> q2_reference;
    std::vector<double> omega2_reference;
};
