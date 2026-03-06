#include "tests/test_framework.h"

#include "modules/damped_config.h"
#include "modules/driven_config.h"
#include "modules/experiment_config.h"
#include "tests/test_helpers.h"

TEST(ExperimentConfigLoadsValidYaml) {
    EnvVarGuard guard("QA_TEST");
    guard.set("1");

    TempDir temp;
    const auto path = temp.write_file(
        "simple.yaml",
        "length: 2.0\n"
        "gravity: 9.5\n"
        "dt: 0.01\n"
        "t_max: 5.0\n"
        "theta0: 0.3\n"
        "omega0: -0.1\n"
        "data_file: \"qa/out.csv\"\n"
        "show_plot: false\n"
        "save_png: true\n"
        "plot_phase_map: false\n"
        "output_png: \"qa/out.png\"\n"
        "integrator: \"RK5\"\n"
        "analytical_model: \"Jacobi\"\n"
        "error_analysis: hd_reference\n"
        "error_reference_factor: 80\n"
        "restoring_force_model: polynomial\n"
        "restoring_force_linear: 1.2\n"
        "restoring_force_cubic: 0.35\n");

    const ExperimentConfig cfg = load_config_from_yaml(path.string());
    EXPECT_NEAR(cfg.length, 2.0, 1e-12);
    EXPECT_NEAR(cfg.gravity, 9.5, 1e-12);
    EXPECT_NEAR(cfg.dt, 0.01, 1e-12);
    EXPECT_NEAR(cfg.t_max, 5.0, 1e-12);
    EXPECT_NEAR(cfg.theta0, 0.3, 1e-12);
    EXPECT_NEAR(cfg.omega0, -0.1, 1e-12);
    EXPECT_EQ(cfg.data_file, "qa/out.csv");
    EXPECT_FALSE(cfg.show_plot);
    EXPECT_TRUE(cfg.save_png);
    EXPECT_FALSE(cfg.plot_phase_map);
    EXPECT_EQ(cfg.output_png, "qa/out.png");
    EXPECT_EQ(cfg.integrator, "rk5");
    EXPECT_EQ(cfg.analytical_model, "jacobi");
    EXPECT_TRUE(cfg.error_mode == error_reference::Mode::HdReference);
    EXPECT_EQ(cfg.error_reference_factor, 80);
    EXPECT_TRUE(cfg.restoring_force.model == restoring_force::Model::Polynomial);
    EXPECT_NEAR(cfg.restoring_force.linear, 1.2, 1e-12);
    EXPECT_NEAR(cfg.restoring_force.cubic, 0.35, 1e-12);
}

TEST(ExperimentConfigMissingOptionalUsesDefaults) {
    EnvVarGuard guard("QA_TEST");
    guard.set("1");

    TempDir temp;
    const auto path = temp.write_file(
        "minimal.yaml",
        "length: 1.2\n"
        "gravity: 9.81\n"
        "dt: 0.02\n"
        "t_max: 2.0\n"
        "theta0: 0.1\n"
        "omega0: 0.0\n");

    const ExperimentConfig cfg = load_config_from_yaml(path.string());
    EXPECT_EQ(cfg.data_file, "simple_pendulum_data.csv");
    EXPECT_TRUE(cfg.show_plot);
    EXPECT_TRUE(cfg.save_png);
    EXPECT_TRUE(cfg.plot_phase_map);
    EXPECT_EQ(cfg.output_png, "simple_pendulum.png");
    EXPECT_EQ(cfg.integrator, "rk4");
    EXPECT_EQ(cfg.analytical_model, "linear");
    EXPECT_TRUE(cfg.error_mode == error_reference::Mode::Analytical);
    EXPECT_EQ(cfg.error_reference_factor, 50);
    EXPECT_TRUE(cfg.restoring_force.model == restoring_force::Model::Sine);
    EXPECT_NEAR(cfg.restoring_force.linear, 1.0, 1e-12);
    EXPECT_NEAR(cfg.restoring_force.cubic, 0.0, 1e-12);
}

TEST(ExperimentConfigInvalidYamlThrows) {
    TempDir temp;
    const auto bad_syntax = temp.write_file("bad_syntax.yaml", "length 1.0\n");
    const auto bad_type = temp.write_file(
        "bad_type.yaml",
        "length: 1.0\n"
        "gravity: abc\n"
        "dt: 0.01\n"
        "t_max: 1.0\n");
    const auto missing_value = temp.write_file("missing_value.yaml", "length:\n");

    EXPECT_THROW(load_config_from_yaml(bad_syntax.string()));
    EXPECT_THROW(load_config_from_yaml(bad_type.string()));
    EXPECT_THROW(load_config_from_yaml(missing_value.string()));
}

TEST(ExperimentConfigOutputPathResolvesWithoutQaEnv) {
    EnvVarGuard guard("QA_TEST");
    guard.unset();

    TempDir temp;
    WorkingDirGuard wd(temp.path());
    const auto path = temp.write_file(
        "config.yaml",
        "length: 1.0\n"
        "gravity: 9.81\n"
        "dt: 0.01\n"
        "t_max: 1.0\n"
        "theta0: 0.2\n"
        "omega0: 0.0\n"
        "data_file: \"nested/a.csv\"\n"
        "output_png: \"nested/b.png\"\n");

    const ExperimentConfig cfg = load_config_from_yaml(path.string());
    EXPECT_EQ(cfg.data_file, "outputs/a.csv");
    EXPECT_EQ(cfg.output_png, "outputs/b.png");
}

TEST(DampedConfigLoadsValidYamlAndNestedSections) {
    EnvVarGuard guard("QA_TEST");
    guard.set("1");

    TempDir temp;
    const auto path = temp.write_file(
        "damped.yaml",
        "physical:\n"
        "  g: 9.81\n"
        "  L: 2.0\n"
        "  gamma: 0.3\n"
        "  theta0: 0.4\n"
        "  theta_dot0: -0.2\n"
        "  restoring_force_model: polynomial\n"
        "  restoring_force_linear: 0.9\n"
        "  restoring_force_cubic: 0.2\n"
        "simulation:\n"
        "  t_start: 0.0\n"
        "  t_end: 3.0\n"
        "  dt: 0.01\n"
        "  output_every: 5\n"
        "settings:\n"
        "  plotting_method: new\n"
        "  show_plot: false\n"
        "  save_png: false\n"
        "  plot_phase_map: true\n"
        "  error_analysis: none\n"
        "  error_reference_factor: 64\n"
        "  run_plotter: false\n"
        "  data_file: \"damped.dat\"\n"
        "  output_png: \"damped.png\"\n"
        "  python_script: \"plot.py\"\n"
        "integrator: \"RK3\"\n");

    const DampedConfig cfg = load_damped_config_from_yaml(path.string());
    EXPECT_NEAR(cfg.physical.L, 2.0, 1e-12);
    EXPECT_NEAR(cfg.physical.gamma, 0.3, 1e-12);
    EXPECT_NEAR(cfg.simulation.t_end, 3.0, 1e-12);
    EXPECT_EQ(cfg.simulation.output_every, 5);
    EXPECT_EQ(cfg.settings.integrator, "rk3");
    EXPECT_EQ(cfg.settings.data_file, "damped.dat");
    EXPECT_EQ(to_string(cfg.settings.plotting_method), "new");
    EXPECT_TRUE(cfg.settings.plot_phase_map);
    EXPECT_TRUE(cfg.settings.error_mode == error_reference::Mode::None);
    EXPECT_EQ(cfg.settings.error_reference_factor, 64);
    EXPECT_TRUE(cfg.physical.restoring_force.model == restoring_force::Model::Polynomial);
    EXPECT_NEAR(cfg.physical.restoring_force.linear, 0.9, 1e-12);
    EXPECT_NEAR(cfg.physical.restoring_force.cubic, 0.2, 1e-12);
}

TEST(DampedConfigDefaultsAndValidation) {
    EnvVarGuard guard("QA_TEST");
    guard.set("1");

    TempDir temp;
    const auto minimal = temp.write_file("damped_minimal.yaml", "");
    const auto invalid = temp.write_file(
        "damped_invalid.yaml",
        "physical:\n"
        "  g: -1\n");

    const DampedConfig cfg = load_damped_config_from_yaml(minimal.string());
    EXPECT_NEAR(cfg.physical.g, 9.81, 1e-12);
    EXPECT_EQ(cfg.settings.integrator, "rk4");
    EXPECT_FALSE(cfg.settings.plot_phase_map);
    EXPECT_TRUE(cfg.settings.error_mode == error_reference::Mode::Analytical);
    EXPECT_EQ(cfg.settings.error_reference_factor, 50);
    EXPECT_TRUE(cfg.physical.restoring_force.model == restoring_force::Model::Sine);
    EXPECT_EQ(to_string(PlottingMethod::Original), "original");
    EXPECT_EQ(to_string(PlottingMethod::New), "new");

    EXPECT_THROW(load_damped_config_from_yaml(invalid.string()));
}

TEST(DrivenConfigLoadsValidYamlAndValidation) {
    EnvVarGuard guard("QA_TEST");
    guard.set("1");

    TempDir temp;
    const auto path = temp.write_file(
        "driven.yaml",
        "physical:\n"
        "  g: 9.81\n"
        "  L: 1.3\n"
        "  damping: 0.2\n"
        "  A: 0.4\n"
        "  omega_drive: 1.1\n"
        "  theta0: 0.05\n"
        "  omega0: 0.01\n"
        "  restoring_force_model: polynomial\n"
        "  restoring_force_linear: 1.05\n"
        "  restoring_force_cubic: 0.15\n"
        "simulation:\n"
        "  t_start: 0.0\n"
        "  t_end: 4.0\n"
        "  dt: 0.01\n"
        "  output_every: 2\n"
        "settings:\n"
        "  plotting_method: original\n"
        "  show_plot: false\n"
        "  save_png: false\n"
        "  plot_phase_map: true\n"
        "  error_analysis: hd_reference\n"
        "  error_reference_factor: 55\n"
        "  run_plotter: false\n"
        "  data_file: \"driven.csv\"\n"
        "  output_png: \"driven.png\"\n"
        "  python_script: \"driven_plot.py\"\n"
        "integrator: \"RK5\"\n");

    const DrivenConfig cfg = load_driven_config_from_yaml(path.string());
    EXPECT_NEAR(cfg.physical.L, 1.3, 1e-12);
    EXPECT_NEAR(cfg.physical.A, 0.4, 1e-12);
    EXPECT_EQ(cfg.simulation.output_every, 2);
    EXPECT_EQ(cfg.settings.integrator, "rk5");
    EXPECT_EQ(cfg.settings.data_file, "driven.csv");
    EXPECT_EQ(to_string(cfg.settings.plotting_method), "original");
    EXPECT_TRUE(cfg.settings.plot_phase_map);
    EXPECT_TRUE(cfg.settings.error_mode == error_reference::Mode::HdReference);
    EXPECT_EQ(cfg.settings.error_reference_factor, 55);
    EXPECT_TRUE(cfg.physical.restoring_force.model == restoring_force::Model::Polynomial);
    EXPECT_NEAR(cfg.physical.restoring_force.linear, 1.05, 1e-12);
    EXPECT_NEAR(cfg.physical.restoring_force.cubic, 0.15, 1e-12);

    const auto invalid = temp.write_file(
        "driven_invalid.yaml",
        "simulation:\n"
        "  dt: 0\n");
    EXPECT_THROW(load_driven_config_from_yaml(invalid.string()));
}

TEST(DrivenConfigOutputPathResolvesWithoutQaEnv) {
    EnvVarGuard guard("QA_TEST");
    guard.unset();

    TempDir temp;
    WorkingDirGuard wd(temp.path());
    const auto path = temp.write_file(
        "driven_paths.yaml",
        "settings:\n"
        "  data_file: \"a/b/c.csv\"\n"
        "  output_png: \"a/b/c.png\"\n"
        "  python_script: \"a/b/c.py\"\n");

    const DrivenConfig cfg = load_driven_config_from_yaml(path.string());
    EXPECT_EQ(cfg.settings.data_file, "outputs/c.csv");
    EXPECT_EQ(cfg.settings.output_png, "outputs/c.png");
    EXPECT_EQ(cfg.settings.python_script, "outputs/c.py");
    EXPECT_FALSE(cfg.settings.plot_phase_map);
}
