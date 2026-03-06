#include "tests/test_framework.h"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "modules/damped_config.h"
#include "modules/damped_io.h"
#include "modules/damped_simulator.h"
#include "modules/driven_config.h"
#include "modules/driven_io.h"
#include "modules/driven_simulator.h"
#include "modules/experiment_config.h"
#include "modules/pendulum_simulator.h"
#include "modules/simple_io.h"
#include "tests/test_helpers.h"

namespace {

struct Trajectory {
    double dt = 0.0;
    std::vector<double> theta;
    std::vector<double> omega;
};

bool env_truthy(const char* name) {
    const char* raw = std::getenv(name);
    if (raw == nullptr) {
        return false;
    }
    const std::string value(raw);
    return value == "1" || value == "true" || value == "TRUE" || value == "yes" || value == "on";
}

std::string to_yaml_num(double value) {
    std::ostringstream oss;
    oss << std::setprecision(17) << value;
    return oss.str();
}

Trajectory run_simple_serial_case(
    const TempDir& temp,
    const std::string& tag,
    const std::string& integrator,
    double dt,
    double t_max,
    double theta0,
    double omega0) {
    const std::string yaml_name = "serial_" + tag + "_" + integrator + "_" + std::to_string(dt) + ".yaml";
    const std::string data_rel = "serial/" + tag + "_" + integrator + "_" + std::to_string(dt) + ".csv";

    const auto config_path = temp.write_file(
        yaml_name,
        "length: 1.0\n"
        "gravity: 9.81\n"
        "dt: " + to_yaml_num(dt) + "\n"
        "t_max: " + to_yaml_num(t_max) + "\n"
        "theta0: " + to_yaml_num(theta0) + "\n"
        "omega0: " + to_yaml_num(omega0) + "\n"
        "integrator: \"" + integrator + "\"\n"
        "analytical_model: jacobi\n"
        "data_file: \"" + data_rel + "\"\n"
        "show_plot: false\n"
        "save_png: false\n");

    const ExperimentConfig cfg = load_config_from_yaml(config_path.string());
    const auto output_path = temp.child(cfg.data_file);
    std::filesystem::create_directories(output_path.parent_path());

    PendulumSimulator simulator(
        cfg.length, cfg.gravity, cfg.dt, cfg.t_max, cfg.restoring_force);
    const SimulationResult result = simulator.simulate(
        cfg.theta0, cfg.omega0, cfg.integrator, cfg.analytical_model,
        cfg.error_mode, cfg.error_reference_factor);
    write_simple_data_file(output_path.string(), result);

    const CsvData csv = read_csv(output_path);
    EXPECT_EQ(csv.rows.size(), result.t.size());

    Trajectory out;
    out.dt = cfg.dt;
    out.theta = result.theta;
    out.omega = result.omega;
    return out;
}

double trajectory_rms_error(const Trajectory& candidate, const Trajectory& reference) {
    const double ratio_real = candidate.dt / reference.dt;
    const size_t ratio = static_cast<size_t>(std::llround(ratio_real));
    EXPECT_TRUE(ratio > 0u);
    EXPECT_NEAR(static_cast<double>(ratio), ratio_real, 1e-9);

    EXPECT_TRUE(reference.theta.size() >= (candidate.theta.size() - 1u) * ratio + 1u);
    EXPECT_EQ(candidate.theta.size(), candidate.omega.size());
    EXPECT_EQ(reference.theta.size(), reference.omega.size());

    double sq_sum = 0.0;
    for (size_t i = 0; i < candidate.theta.size(); ++i) {
        const size_t ref_i = i * ratio;
        const double d_theta = candidate.theta[i] - reference.theta[ref_i];
        const double d_omega = candidate.omega[i] - reference.omega[ref_i];
        sq_sum += d_theta * d_theta + d_omega * d_omega;
    }
    return std::sqrt(sq_sum / static_cast<double>(candidate.theta.size()));
}

double loglog_slope(const std::vector<double>& dts, const std::vector<double>& errors) {
    EXPECT_EQ(dts.size(), errors.size());
    const size_t n = dts.size();
    EXPECT_TRUE(n >= 2u);

    double sum_x = 0.0;
    double sum_y = 0.0;
    double sum_xx = 0.0;
    double sum_xy = 0.0;

    for (size_t i = 0; i < n; ++i) {
        const double x = std::log(dts[i]);
        const double y = std::log(std::max(errors[i], 1e-30));
        sum_x += x;
        sum_y += y;
        sum_xx += x * x;
        sum_xy += x * y;
    }

    const double denom = static_cast<double>(n) * sum_xx - sum_x * sum_x;
    EXPECT_TRUE(std::fabs(denom) > 1e-18);
    return (static_cast<double>(n) * sum_xy - sum_x * sum_y) / denom;
}

std::string format_dt_for_path(double dt) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(9) << dt;
    std::string s = oss.str();
    while (!s.empty() && s.back() == '0') {
        s.pop_back();
    }
    if (!s.empty() && s.back() == '.') {
        s.pop_back();
    }
    for (char& c : s) {
        if (c == '.') {
            c = 'p';
        }
    }
    return s.empty() ? "0" : s;
}

void write_convergence_csv(
    const std::vector<double>& dts,
    const std::map<std::string, std::vector<double>>& errors_by_method) {
    const std::filesystem::path out_dir = "tests/artifacts/convergence";
    std::filesystem::create_directories(out_dir);
    const std::filesystem::path out_csv = out_dir / "convergence_results.csv";

    std::ofstream file(out_csv);
    if (!file.is_open()) {
        throw std::runtime_error("Could not write convergence CSV: " + out_csv.string());
    }

    file << "integrator,dt,rms_error\n";
    for (const auto& [method, errors] : errors_by_method) {
        EXPECT_EQ(errors.size(), dts.size());
        for (size_t i = 0; i < dts.size(); ++i) {
            file << method << ","
                 << std::setprecision(17) << dts[i] << ","
                 << std::setprecision(17) << errors[i] << "\n";
        }
    }
}

void maybe_render_convergence_plot() {
    if (!env_truthy("PENDULUM_CONVERGENCE_PLOT")) {
        return;
    }

    const bool show = env_truthy("PENDULUM_CONVERGENCE_SHOW");
    std::string command =
        "python3 tests/plot_convergence.py "
        "--input tests/artifacts/convergence/convergence_results.csv "
        "--output tests/artifacts/convergence/convergence_results.png";
    if (show) {
        command += " --show";
    }

    const int rc = std::system(command.c_str());
    EXPECT_EQ(rc, 0);
}

}  // namespace

TEST(SerialIntegrationSimplePendulumFlowWritesExpectedCsv) {
    EnvVarGuard guard("QA_TEST");
    guard.set("1");

    TempDir temp;
    const auto config_path = temp.write_file(
        "simple.yaml",
        "length: 1.0\n"
        "gravity: 9.81\n"
        "dt: 0.01\n"
        "t_max: 1.0\n"
        "theta0: 0.1\n"
        "omega0: 0.0\n"
        "integrator: rk4\n"
        "analytical_model: linear\n"
        "data_file: \"serial/simple.csv\"\n"
        "show_plot: false\n"
        "save_png: false\n");

    const ExperimentConfig cfg = load_config_from_yaml(config_path.string());
    const auto output_path = temp.child(cfg.data_file);
    std::filesystem::create_directories(output_path.parent_path());

    PendulumSimulator simulator(
        cfg.length, cfg.gravity, cfg.dt, cfg.t_max, cfg.restoring_force);
    const SimulationResult result = simulator.simulate(
        cfg.theta0, cfg.omega0, cfg.integrator, cfg.analytical_model,
        cfg.error_mode, cfg.error_reference_factor);
    write_simple_data_file(output_path.string(), result);

    const CsvData csv = read_csv(output_path);
    const size_t expected_rows = static_cast<size_t>(std::round(cfg.t_max / cfg.dt)) + 1u;
    EXPECT_EQ(csv.rows.size(), expected_rows);
    EXPECT_NEAR(csv.rows.front()[0], 0.0, 1e-9);
    EXPECT_NEAR(csv.rows.back()[0], cfg.t_max, 1e-6);
    EXPECT_TRUE(std::fabs(csv.rows.back()[2]) < 1.0);
}

TEST(SerialIntegrationDuffingUndampedFlowWritesExpectedCsv) {
    EnvVarGuard guard("QA_TEST");
    guard.set("1");

    TempDir temp;
    const auto config_path = temp.write_file(
        "duffing.yaml",
        "length: 1.0\n"
        "gravity: 9.81\n"
        "dt: 0.001\n"
        "t_max: 3.0\n"
        "theta0: 0.45\n"
        "omega0: 0.0\n"
        "integrator: rk5\n"
        "analytical_model: duffing_jacobi\n"
        "restoring_force_model: polynomial\n"
        "restoring_force_linear: 1.0\n"
        "restoring_force_cubic: 0.35\n"
        "data_file: \"serial/duffing.csv\"\n"
        "show_plot: false\n"
        "save_png: false\n");

    const ExperimentConfig cfg = load_config_from_yaml(config_path.string());
    const auto output_path = temp.child(cfg.data_file);
    std::filesystem::create_directories(output_path.parent_path());

    PendulumSimulator simulator(
        cfg.length, cfg.gravity, cfg.dt, cfg.t_max, cfg.restoring_force);
    const SimulationResult result =
        simulator.simulate(cfg.theta0, cfg.omega0, cfg.integrator, cfg.analytical_model,
                           cfg.error_mode, cfg.error_reference_factor);
    write_simple_data_file(output_path.string(), result);

    const CsvData csv = read_csv(output_path);
    const size_t expected_rows = static_cast<size_t>(std::round(cfg.t_max / cfg.dt)) + 1u;
    EXPECT_EQ(csv.rows.size(), expected_rows);
    EXPECT_TRUE(result.theta_stats.max_abs < 1e-2);
}

TEST(SerialIntegrationSimplePendulumHdReferenceFlowWritesExpectedCsv) {
    EnvVarGuard guard("QA_TEST");
    guard.set("1");

    TempDir temp;
    const auto config_path = temp.write_file(
        "simple_hd.yaml",
        "length: 1.0\n"
        "gravity: 9.81\n"
        "dt: 0.01\n"
        "t_max: 2.0\n"
        "theta0: 0.6\n"
        "omega0: 0.0\n"
        "integrator: rk4\n"
        "analytical_model: jacobi\n"
        "error_analysis: hd_reference\n"
        "error_reference_factor: 50\n"
        "data_file: \"serial/simple_hd.csv\"\n"
        "show_plot: false\n"
        "save_png: false\n");

    const ExperimentConfig cfg = load_config_from_yaml(config_path.string());
    const auto output_path = temp.child(cfg.data_file);
    std::filesystem::create_directories(output_path.parent_path());

    PendulumSimulator simulator(
        cfg.length, cfg.gravity, cfg.dt, cfg.t_max, cfg.restoring_force);
    const SimulationResult result = simulator.simulate(
        cfg.theta0, cfg.omega0, cfg.integrator, cfg.analytical_model,
        cfg.error_mode, cfg.error_reference_factor);
    write_simple_data_file(output_path.string(), result);

    const CsvData csv = read_csv(output_path);
    EXPECT_FALSE(csv.rows.empty());
    EXPECT_TRUE(result.theta_stats.max_abs < 5e-2);
}

TEST(SerialIntegrationDampedPendulumFlowWritesExpectedDat) {
    EnvVarGuard guard("QA_TEST");
    guard.set("1");

    TempDir temp;
    const auto config_path = temp.write_file(
        "damped.yaml",
        "physical:\n"
        "  g: 9.81\n"
        "  L: 1.0\n"
        "  gamma: 0.2\n"
        "  theta0: 0.4\n"
        "  theta_dot0: 0.0\n"
        "simulation:\n"
        "  t_start: 0.0\n"
        "  t_end: 2.0\n"
        "  dt: 0.01\n"
        "settings:\n"
        "  data_file: \"serial/damped.dat\"\n"
        "  output_png: \"serial/damped.png\"\n"
        "  python_script: \"serial/plot.py\"\n"
        "  run_plotter: false\n"
        "  show_plot: false\n"
        "  save_png: false\n");

    const DampedConfig cfg = load_damped_config_from_yaml(config_path.string());
    const auto output_path = temp.child(cfg.settings.data_file);
    std::filesystem::create_directories(output_path.parent_path());

    DampedPendulumSimulator simulator(cfg);
    const SimulationResult result = simulator.simulate();
    write_damped_data_file(output_path.string(), result);

    const DatData dat = read_dat(output_path);
    const size_t expected_rows =
        static_cast<size_t>(std::round((cfg.simulation.t_end - cfg.simulation.t_start) / cfg.simulation.dt)) + 1u;
    EXPECT_EQ(dat.rows.size(), expected_rows);
    EXPECT_NEAR(dat.rows.front()[0], cfg.simulation.t_start, 1e-9);
    EXPECT_NEAR(dat.rows.back()[0], cfg.simulation.t_end, 1e-6);
    EXPECT_TRUE(dat.rows.back()[6] < dat.rows.front()[6]);
}

TEST(SerialIntegrationDrivenPendulumFlowWritesExpectedCsv) {
    EnvVarGuard guard("QA_TEST");
    guard.set("1");

    TempDir temp;
    const auto config_path = temp.write_file(
        "driven.yaml",
        "physical:\n"
        "  g: 9.81\n"
        "  L: 1.0\n"
        "  damping: 0.4\n"
        "  A: 0.3\n"
        "  omega_drive: 1.0\n"
        "  theta0: 0.0\n"
        "  omega0: 0.0\n"
        "simulation:\n"
        "  t_start: 0.0\n"
        "  t_end: 2.5\n"
        "  dt: 0.01\n"
        "settings:\n"
        "  data_file: \"serial/driven.csv\"\n"
        "  output_png: \"serial/driven.png\"\n"
        "  python_script: \"serial/plot.py\"\n"
        "  run_plotter: false\n"
        "  show_plot: false\n"
        "  save_png: false\n");

    const DrivenConfig cfg = load_driven_config_from_yaml(config_path.string());
    const auto output_path = temp.child(cfg.settings.data_file);
    std::filesystem::create_directories(output_path.parent_path());

    DrivenPendulumSimulator simulator(cfg);
    const SimulationResult result = simulator.simulate();
    write_driven_data_file(output_path.string(), result);

    const CsvData csv = read_csv(output_path);
    const size_t expected_rows =
        static_cast<size_t>(std::round((cfg.simulation.t_end - cfg.simulation.t_start) / cfg.simulation.dt)) + 1u;
    EXPECT_EQ(csv.rows.size(), expected_rows);
    EXPECT_NEAR(csv.rows.front()[0], cfg.simulation.t_start, 1e-9);
    EXPECT_NEAR(csv.rows.back()[0], cfg.simulation.t_end, 1e-6);
    EXPECT_TRUE(std::fabs(csv.rows.back()[2]) < 5.0);
}

TEST(SerialIntegrationSimplePendulumConvergesWithIncreasingResolutionAcrossIntegrators) {
    EnvVarGuard guard("QA_TEST");
    guard.set("1");

    TempDir temp;
    const double t_max = 2.0;
    const double theta0 = 0.7;
    const double omega0 = -0.2;

    const std::vector<double> dts = {
        0.08,
        0.04,
        0.02,
        0.01,
        0.005,
        0.0025,
        0.00125,
        0.000625,
        0.0003125,
        0.00015625,
        0.000078125,
    };

    const double reference_dt = 0.0000390625;
    const Trajectory reference =
        run_simple_serial_case(temp, "reference", "rk5", reference_dt, t_max, theta0, omega0);

    std::map<std::string, std::vector<double>> errors_by_method;

    for (const std::string& method : std::array<std::string, 8>{"rk3", "rk4", "rk5", "rk23", "rkf45", "semi_implicit_euler", "leapfrog", "ruth4"}) {
        std::vector<double> errors;
        errors.reserve(dts.size());
        for (double dt : dts) {
            const Trajectory traj = run_simple_serial_case(
                temp,
                "scan_" + format_dt_for_path(dt),
                method,
                dt,
                t_max,
                theta0,
                omega0);
            errors.push_back(trajectory_rms_error(traj, reference));
        }

        const double slope = loglog_slope(dts, errors);
        EXPECT_TRUE(slope > 1.0);
        EXPECT_TRUE(errors.back() < errors.front() * 0.2);
        errors_by_method[method] = std::move(errors);
    }

    write_convergence_csv(dts, errors_by_method);
    maybe_render_convergence_plot();
}
