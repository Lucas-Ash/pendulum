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

#include "modules/config/damped_config.h"
#include "modules/config/driven_config.h"
#include "modules/config/experiment_config.h"
#include "modules/damped/damped_io.h"
#include "modules/damped/damped_simulator.h"
#include "modules/driven/driven_io.h"
#include "modules/driven/driven_simulator.h"
#include "modules/simple/pendulum_simulator.h"
#include "modules/simple/simple_io.h"
#include "tests/test_helpers.h"

namespace {

struct Trajectory {
    double dt = 0.0;
    std::vector<double> theta;
    std::vector<double> omega;
    std::vector<double> energy;
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

Trajectory result_to_trajectory(const SimulationResult& result, double dt) {
    Trajectory out;
    out.dt = dt;
    out.theta = result.theta;
    out.omega = result.omega;
    out.energy = result.energy;
    return out;
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

    return result_to_trajectory(result, cfg.dt);
}

Trajectory run_damped_serial_case(const TempDir& temp,
                                  const std::string& integrator,
                                  double dt,
                                  double t_end,
                                  double gamma,
                                  double theta0,
                                  double theta_dot0) {
    const std::string yaml = "physical:\n"
        "  g: 9.81\n"
        "  L: 1.0\n"
        "  gamma: " + to_yaml_num(gamma) + "\n"
        "  theta0: " + to_yaml_num(theta0) + "\n"
        "  theta_dot0: " + to_yaml_num(theta_dot0) + "\n"
        "simulation:\n"
        "  t_start: 0.0\n"
        "  t_end: " + to_yaml_num(t_end) + "\n"
        "  dt: " + to_yaml_num(dt) + "\n"
        "  output_every: 1\n"
        "settings:\n"
        "  data_file: \"serial/damped_" + integrator + ".dat\"\n"
        "  output_png: \"serial/damped.png\"\n"
        "  python_script: \"serial/plot.py\"\n"
        "  integrator: \"" + integrator + "\"\n"
        "  error_analysis: hd_reference\n"
        "  error_reference_factor: 50\n"
        "  run_plotter: false\n"
        "  show_plot: false\n"
        "  save_png: false\n";
    const auto config_path = temp.write_file("damped_" + integrator + ".yaml", yaml);
    const DampedConfig cfg = load_damped_config_from_yaml(config_path.string());
    DampedPendulumSimulator simulator(cfg);
    const SimulationResult result = simulator.simulate();
    return result_to_trajectory(result, cfg.simulation.dt);
}

Trajectory run_driven_serial_case(const TempDir& temp,
                                  const std::string& integrator,
                                  double dt,
                                  double t_end,
                                  double damping,
                                  double A,
                                  double omega_drive,
                                  double theta0,
                                  double omega0) {
    const std::string yaml = "physical:\n"
        "  g: 9.81\n"
        "  L: 1.0\n"
        "  damping: " + to_yaml_num(damping) + "\n"
        "  A: " + to_yaml_num(A) + "\n"
        "  omega_drive: " + to_yaml_num(omega_drive) + "\n"
        "  theta0: " + to_yaml_num(theta0) + "\n"
        "  omega0: " + to_yaml_num(omega0) + "\n"
        "simulation:\n"
        "  t_start: 0.0\n"
        "  t_end: " + to_yaml_num(t_end) + "\n"
        "  dt: " + to_yaml_num(dt) + "\n"
        "  output_every: 1\n"
        "settings:\n"
        "  data_file: \"serial/driven_" + integrator + ".csv\"\n"
        "  output_png: \"serial/driven.png\"\n"
        "  python_script: \"serial/plot.py\"\n"
        "  integrator: \"" + integrator + "\"\n"
        "  error_analysis: hd_reference\n"
        "  error_reference_factor: 50\n"
        "  run_plotter: false\n"
        "  show_plot: false\n"
        "  save_png: false\n";
    const auto config_path = temp.write_file("driven_" + integrator + ".yaml", yaml);
    const DrivenConfig cfg = load_driven_config_from_yaml(config_path.string());
    DrivenPendulumSimulator simulator(cfg);
    const SimulationResult result = simulator.simulate();
    return result_to_trajectory(result, cfg.simulation.dt);
}

size_t checked_ratio(const Trajectory& candidate, const Trajectory& reference) {
    const double ratio_real = candidate.dt / reference.dt;
    const size_t ratio = static_cast<size_t>(std::llround(ratio_real));
    EXPECT_TRUE(ratio > 0u);
    EXPECT_NEAR(static_cast<double>(ratio), ratio_real, 1e-9);

    EXPECT_TRUE(reference.theta.size() >= (candidate.theta.size() - 1u) * ratio + 1u);
    EXPECT_EQ(candidate.theta.size(), candidate.omega.size());
    EXPECT_EQ(reference.theta.size(), reference.omega.size());
    EXPECT_EQ(candidate.theta.size(), candidate.energy.size());
    EXPECT_EQ(reference.theta.size(), reference.energy.size());
    return ratio;
}

double component_rms_error(const std::vector<double>& candidate,
                           const std::vector<double>& reference,
                           size_t ratio) {
    EXPECT_TRUE(reference.size() >= (candidate.size() - 1u) * ratio + 1u);

    double sq_sum = 0.0;
    for (size_t i = 0; i < candidate.size(); ++i) {
        const size_t ref_i = i * ratio;
        const double delta = candidate[i] - reference[ref_i];
        sq_sum += delta * delta;
    }
    return std::sqrt(sq_sum / static_cast<double>(candidate.size()));
}

double trajectory_rms_error(const Trajectory& candidate, const Trajectory& reference) {
    const size_t ratio = checked_ratio(candidate, reference);

    double sq_sum = 0.0;
    for (size_t i = 0; i < candidate.theta.size(); ++i) {
        const size_t ref_i = i * ratio;
        const double d_theta = candidate.theta[i] - reference.theta[ref_i];
        const double d_omega = candidate.omega[i] - reference.omega[ref_i];
        sq_sum += d_theta * d_theta + d_omega * d_omega;
    }
    return std::sqrt(sq_sum / static_cast<double>(candidate.theta.size()));
}

double trajectory_theta_rms_error(const Trajectory& candidate, const Trajectory& reference) {
    return component_rms_error(candidate.theta, reference.theta, checked_ratio(candidate, reference));
}

double trajectory_omega_rms_error(const Trajectory& candidate, const Trajectory& reference) {
    return component_rms_error(candidate.omega, reference.omega, checked_ratio(candidate, reference));
}

double trajectory_energy_rms_error(const Trajectory& candidate, const Trajectory& reference) {
    return component_rms_error(candidate.energy, reference.energy, checked_ratio(candidate, reference));
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

double filtered_loglog_slope(const std::vector<double>& dts, const std::vector<double>& errors) {
    EXPECT_EQ(dts.size(), errors.size());
    EXPECT_TRUE(dts.size() >= 3u);

    std::vector<double> sorted_errors = errors;
    std::sort(sorted_errors.begin(), sorted_errors.end());
    const size_t tail_count = std::min<size_t>(3u, sorted_errors.size());
    const double floor_est = sorted_errors[tail_count / 2u];
    const double threshold = floor_est * 8.0;

    size_t fit_count = errors.size();
    for (size_t i = 0; i < errors.size(); ++i) {
        if (errors[i] <= threshold && i >= 3u) {
            fit_count = i;
            break;
        }
    }
    if (fit_count < 3u) {
        fit_count = errors.size();
    }

    std::vector<double> fit_dts(dts.begin(), dts.begin() + static_cast<std::ptrdiff_t>(fit_count));
    std::vector<double> fit_errors(errors.begin(), errors.begin() + static_cast<std::ptrdiff_t>(fit_count));
    return loglog_slope(fit_dts, fit_errors);
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

void write_convergence_csv_to_path(
    const std::vector<double>& dts,
    const std::map<std::string, std::vector<double>>& errors_by_method,
    const std::filesystem::path& out_csv) {
    std::filesystem::create_directories(out_csv.parent_path());
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

void write_convergence_csv(
    const std::vector<double>& dts,
    const std::map<std::string, std::vector<double>>& errors_by_method) {
    write_convergence_csv_to_path(dts, errors_by_method,
                                  "tests/artifacts/convergence/convergence_results.csv");
}

void write_numerov_convergence_report(
    const std::vector<double>& dts,
    const std::vector<double>& theta_errors,
    const std::vector<double>& omega_errors,
    const std::vector<double>& state_errors,
    const std::vector<double>& energy_errors) {
    const std::filesystem::path out_dir = "tests/artifacts/convergence";
    std::filesystem::create_directories(out_dir);

    const std::filesystem::path report_csv = out_dir / "numerov_convergence_report.csv";
    std::ofstream report(report_csv);
    if (!report.is_open()) {
        throw std::runtime_error("Could not write Numerov convergence report: " + report_csv.string());
    }

    EXPECT_EQ(theta_errors.size(), dts.size());
    EXPECT_EQ(omega_errors.size(), dts.size());
    EXPECT_EQ(state_errors.size(), dts.size());
    EXPECT_EQ(energy_errors.size(), dts.size());

    report << "dt,rms_theta_error,rms_omega_error,rms_state_error,rms_energy_error\n";
    for (size_t i = 0; i < dts.size(); ++i) {
        report << std::setprecision(17) << dts[i] << ","
               << theta_errors[i] << ","
               << omega_errors[i] << ","
               << state_errors[i] << ","
               << energy_errors[i] << "\n";
    }

    const std::filesystem::path summary_csv = out_dir / "numerov_convergence_summary.csv";
    std::ofstream summary(summary_csv);
    if (!summary.is_open()) {
        throw std::runtime_error("Could not write Numerov convergence summary: " + summary_csv.string());
    }

    summary << "metric,filtered_gradient\n";
    summary << "theta," << std::setprecision(12) << filtered_loglog_slope(dts, theta_errors) << "\n";
    summary << "omega," << std::setprecision(12) << filtered_loglog_slope(dts, omega_errors) << "\n";
    summary << "state," << std::setprecision(12) << filtered_loglog_slope(dts, state_errors) << "\n";
    summary << "energy," << std::setprecision(12) << filtered_loglog_slope(dts, energy_errors) << "\n";
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

    int rc = std::system(command.c_str());
    EXPECT_EQ(rc, 0);

    command =
        "python3 tests/plot_performance_by_order.py "
        "--input tests/artifacts/convergence/convergence_results.csv "
        "--output tests/artifacts/convergence/performance_by_order.png "
        "--t-max 2.0";
    rc = std::system(command.c_str());
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

TEST(SerialIntegrationVanDerPolFlowWritesExpectedDat) {
    EnvVarGuard guard("QA_TEST");
    guard.set("1");

    TempDir temp;
    const auto config_path = temp.write_file(
        "van_der_pol.yaml",
        "physical:\n"
        "  g: 1.0\n"
        "  L: 1.0\n"
        "  damping_model: polynomial\n"
        "  damping_linear: -0.05\n"
        "  damping_cubic: 0.05\n"
        "  theta0: 2.0\n"
        "  theta_dot0: 0.0\n"
        "  restoring_force_model: polynomial\n"
        "  restoring_force_linear: 1.0\n"
        "  restoring_force_cubic: 0.0\n"
        "simulation:\n"
        "  t_start: 0.0\n"
        "  t_end: 6.0\n"
        "  dt: 0.01\n"
        "settings:\n"
        "  data_file: \"serial/van_der_pol.dat\"\n"
        "  output_png: \"serial/van_der_pol.png\"\n"
        "  python_script: \"serial/plot.py\"\n"
        "  analytical_model: van_der_pol\n"
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
    EXPECT_FALSE(dat.rows.empty());
    EXPECT_TRUE(result.theta_stats.max_abs < 1e-2);
    EXPECT_TRUE(result.omega_stats.max_abs < 2e-2);
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

TEST(SerialIntegrationDrivenDuffingSweepWritesExpectedCsv) {
    EnvVarGuard guard("QA_TEST");
    guard.set("1");

    TempDir temp;
    const auto config_path = temp.write_file(
        "driven_duffing_sweep.yaml",
        "physical:\n"
        "  system_model: duffing\n"
        "  mass: 1.0\n"
        "  linear_stiffness: 1.0\n"
        "  cubic_stiffness: 50.0\n"
        "  drive_force: 0.4\n"
        "  damping: 0.05\n"
        "simulation:\n"
        "  t_start: 0.0\n"
        "  t_end: 1.0\n"
        "  dt: 0.01\n"
        "sweep:\n"
        "  enabled: true\n"
        "  omega_start: 0.8\n"
        "  omega_end: 2.5\n"
        "  points: 10\n"
        "  settle_time: 8.0\n"
        "settings:\n"
        "  error_mode: none\n"
        "  run_plotter: false\n"
        "  save_png: false\n"
        "  show_plot: false\n"
        "  sweep_data_file: \"serial/driven_sweep.csv\"\n");

    const DrivenConfig cfg = load_driven_config_from_yaml(config_path.string());
    const auto output_path = temp.child(cfg.settings.sweep_data_file);
    std::filesystem::create_directories(output_path.parent_path());

    const DrivenSweepResult result = DrivenPendulumSimulator(cfg).simulate_sweep();
    write_driven_sweep_data_file(output_path.string(), result);

    const CsvData csv = read_csv(output_path);
    EXPECT_EQ(csv.rows.size(), static_cast<size_t>(cfg.sweep.points));
    EXPECT_EQ(csv.rows.front().size(), 7u);
    EXPECT_TRUE(csv.rows[5][1] > 0.0);
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
    std::vector<double> numerov_theta_errors;
    std::vector<double> numerov_omega_errors;
    std::vector<double> numerov_state_errors;
    std::vector<double> numerov_energy_errors;

    const std::map<std::string, std::pair<double, double>> expected_orders = {
        {"rk3", {2.5, 3.5}},
        {"rk4", {3.5, 4.5}},
        {"rk5", {4.3, 5.7}},
        {"rk23", {2.5, 3.5}},
        {"rkf45", {4.3, 5.7}},
        {"den3", {3.5, 4.5}},
        {"semi_implicit_euler", {0.8, 1.3}},
        {"leapfrog", {1.7, 2.3}},
        {"ruth4", {3.3, 4.7}},
        {"velocity_verlet", {1.7, 2.3}},
        {"runge_kutta_nystrom", {3.5, 4.5}},
        {"numerov", {3.3, 4.5}},
    };

    for (const auto& [method, expected_order] : expected_orders) {
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
            const double state_error = trajectory_rms_error(traj, reference);
            errors.push_back(state_error);

            if (method == "numerov") {
                numerov_theta_errors.push_back(trajectory_theta_rms_error(traj, reference));
                numerov_omega_errors.push_back(trajectory_omega_rms_error(traj, reference));
                numerov_state_errors.push_back(state_error);
                numerov_energy_errors.push_back(trajectory_energy_rms_error(traj, reference));
            }
        }

        const double slope = filtered_loglog_slope(dts, errors);
        EXPECT_TRUE(slope > expected_order.first);
        EXPECT_TRUE(slope < expected_order.second);
        EXPECT_TRUE(errors.back() < errors.front() * 0.2);
        errors_by_method[method] = std::move(errors);
    }

    write_convergence_csv(dts, errors_by_method);
    write_numerov_convergence_report(
        dts,
        numerov_theta_errors,
        numerov_omega_errors,
        numerov_state_errors,
        numerov_energy_errors);
    maybe_render_convergence_plot();
}

TEST(SerialIntegrationDampedPendulumConvergesDen3AdaptedBenchmark) {
    EnvVarGuard guard("QA_TEST");
    guard.set("1");

    TempDir temp;
    const double t_end = 2.0;
    const double gamma = 0.2;
    const double theta0 = 0.4;
    const double theta_dot0 = 0.0;

    const std::vector<double> dts = {
        0.08, 0.04, 0.02, 0.01, 0.005, 0.0025, 0.00125, 0.000625, 0.0003125,
    };

    const double reference_dt = 0.00015625;
    const Trajectory reference =
        run_damped_serial_case(temp, "rk5", reference_dt, t_end, gamma, theta0, theta_dot0);

    const std::vector<std::string> methods = {"den3", "rk3", "rk4", "rk5", "rk23", "rkf45"};
    std::map<std::string, std::vector<double>> errors_by_method;

    for (const auto& method : methods) {
        std::vector<double> errors;
        errors.reserve(dts.size());
        for (double dt : dts) {
            const Trajectory traj =
                run_damped_serial_case(temp, method, dt, t_end, gamma, theta0, theta_dot0);
            errors.push_back(trajectory_rms_error(traj, reference));
        }
        errors_by_method[method] = std::move(errors);
    }

    write_convergence_csv_to_path(
        dts, errors_by_method,
        "tests/artifacts/convergence/damped_convergence_results.csv");
}

TEST(SerialIntegrationDrivenPendulumConvergesDen3AdaptedBenchmark) {
    EnvVarGuard guard("QA_TEST");
    guard.set("1");

    TempDir temp;
    const double t_end = 2.0;
    const double damping = 0.4;
    const double A = 0.3;
    const double omega_drive = 1.0;
    const double theta0 = 0.0;
    const double omega0 = 0.0;

    const std::vector<double> dts = {
        0.08, 0.04, 0.02, 0.01, 0.005, 0.0025, 0.00125, 0.000625, 0.0003125,
    };

    const double reference_dt = 0.00015625;
    const Trajectory reference = run_driven_serial_case(
        temp, "rk5", reference_dt, t_end, damping, A, omega_drive, theta0, omega0);

    const std::vector<std::string> methods = {"den3", "rk3", "rk4", "rk5", "rk23", "rkf45"};
    std::map<std::string, std::vector<double>> errors_by_method;

    for (const auto& method : methods) {
        std::vector<double> errors;
        errors.reserve(dts.size());
        for (double dt : dts) {
            const Trajectory traj = run_driven_serial_case(
                temp, method, dt, t_end, damping, A, omega_drive, theta0, omega0);
            errors.push_back(trajectory_rms_error(traj, reference));
        }
        errors_by_method[method] = std::move(errors);
    }

    write_convergence_csv_to_path(
        dts, errors_by_method,
        "tests/artifacts/convergence/driven_convergence_results.csv");
}
