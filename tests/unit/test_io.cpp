#include "tests/test_framework.h"

#include <fstream>

#include "modules/damped_io.h"
#include "modules/driven_io.h"
#include "modules/simple_io.h"
#include "modules/simulation_result.h"
#include "tests/test_helpers.h"

namespace {

SimulationResult sample_result() {
    SimulationResult result;
    result.t = {0.0, 0.1, 0.2};
    result.theta_analytical = {1.0, 0.9, 0.8};
    result.theta = {1.0, 0.89, 0.79};
    result.omega = {0.0, -0.1, -0.2};
    result.theta_errors = {0.0, 0.01, 0.01};
    result.omega_errors = {0.0, 0.01, 0.02};
    result.energy = {1.2, 1.1, 1.0};
    return result;
}

}  // namespace

TEST(SimpleIoWriteThenReadRoundTrip) {
    TempDir temp;
    const auto out_path = temp.child("simple.csv");

    const SimulationResult input = sample_result();
    write_simple_data_file(out_path.string(), input);

    const CsvData csv = read_csv(out_path);
    EXPECT_EQ(csv.header, "Time,Theta_Analytical,Theta,Omega,Theta_Errors,Omega_Errors,Energy");
    EXPECT_EQ(csv.rows.size(), input.t.size());
    EXPECT_EQ(csv.rows[0].size(), 7u);
    EXPECT_NEAR(csv.rows[2][2], input.theta[2], 1e-6);
}

TEST(DampedIoWriteThenReadRoundTrip) {
    TempDir temp;
    const auto out_path = temp.child("damped.dat");

    const SimulationResult input = sample_result();
    write_damped_data_file(out_path.string(), input);

    const DatData dat = read_dat(out_path);
    EXPECT_TRUE(dat.header.find("theta_analytical") != std::string::npos);
    EXPECT_EQ(dat.rows.size(), input.t.size());
    EXPECT_EQ(dat.rows[0].size(), 7u);
    EXPECT_NEAR(dat.rows[1][3], input.omega[1], 1e-10);
}

TEST(DrivenIoWriteThenReadRoundTrip) {
    TempDir temp;
    const auto out_path = temp.child("driven.csv");

    const SimulationResult input = sample_result();
    write_driven_data_file(out_path.string(), input);

    const CsvData csv = read_csv(out_path);
    EXPECT_EQ(csv.header, "Time,Theta_Analytical,Theta,Omega,Theta_Errors,Omega_Errors,Energy");
    EXPECT_EQ(csv.rows.size(), input.t.size());
    EXPECT_EQ(csv.rows[0].size(), 7u);
    EXPECT_NEAR(csv.rows[1][6], input.energy[1], 1e-6);
}

TEST(DampedIoUsesScientificNotation) {
    TempDir temp;
    const auto out_path = temp.child("damped_scientific.dat");

    const SimulationResult input = sample_result();
    write_damped_data_file(out_path.string(), input);

    std::ifstream file(out_path);
    EXPECT_TRUE(file.is_open());
    std::string line;
    std::getline(file, line);  // header
    std::getline(file, line);  // first data row
    EXPECT_TRUE(line.find('e') != std::string::npos || line.find('E') != std::string::npos);
}

TEST(DrivenIoUsesScientificNotation) {
    TempDir temp;
    const auto out_path = temp.child("driven_scientific.csv");

    const SimulationResult input = sample_result();
    write_driven_data_file(out_path.string(), input);

    std::ifstream file(out_path);
    EXPECT_TRUE(file.is_open());
    std::string line;
    std::getline(file, line);  // header
    std::getline(file, line);  // first data row
    EXPECT_TRUE(line.find('e') != std::string::npos || line.find('E') != std::string::npos);
}

TEST(DrivenSweepIoWritesExpectedCsv) {
    TempDir temp;
    const auto out_path = temp.child("driven_sweep.csv");

    DrivenSweepResult result;
    result.samples.push_back({1.0, 0.1, 0.09, 0.08, 0.11, 0.02, -0.03});
    result.samples.push_back({1.2, 0.2, 0.18, 0.16, 0.21, 0.01, -0.02});

    write_driven_sweep_data_file(out_path.string(), result);

    const CsvData csv = read_csv(out_path);
    EXPECT_EQ(csv.header,
              "DriveFrequency,NumericalAmplitude,AnalyticalAmplitude,AnalyticalLowerStable,AnalyticalUpperStable,FinalTheta,FinalOmega");
    EXPECT_EQ(csv.rows.size(), 2u);
    EXPECT_NEAR(csv.rows[1][1], 0.2, 1e-12);
    EXPECT_NEAR(csv.rows[0][4], 0.11, 1e-12);
}
