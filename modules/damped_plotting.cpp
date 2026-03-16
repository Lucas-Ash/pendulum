#include "modules/damped_plotting.h"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include "modules/plotting_utils.h"

namespace {

std::string reference_label(const DampedConfig& config) {
    if (config.settings.error_mode == error_reference::Mode::HdReference) {
        return "High-resolution reference";
    }
    if (config.settings.analytical_model == "van_der_pol" ||
        config.settings.analytical_model == "van_der_pol_first_order") {
        return "Van der Pol asymptotic reference";
    }
    return "Analytical reference";
}

std::string numerical_label(const DampedConfig& config) {
    return config.settings.integrator + " numerical trajectory";
}

bool write_python_plot_script(const DampedConfig& config,
                              const SimulationResult& result) {
    const auto& settings = config.settings;
    const std::string ref_label = reference_label(config);
    const std::string num_label = numerical_label(config);

    std::ofstream script(settings.python_script);
    if (!script.is_open()) {
        std::cerr << "Warning: could not create plotting script: "
                  << settings.python_script << "\n";
        return false;
    }

    script << "#!/usr/bin/env python3\n";
    script << "import numpy as np\n";
    script << "import matplotlib.pyplot as plt\n\n";
    script << "data = np.loadtxt('" << settings.data_file << "', comments='#')\n";
    script << "t = data[:,0]\n";
    script << "theta_exact = data[:,1]\n";
    script << "theta_nl = data[:,2]\n";
    script << "omega_nl = data[:,3]\n";
    script << "err_nl = data[:,4]\n";
    script << "err_omega = data[:,5]\n";
    script << "energy_nl = data[:,6]\n\n";
    if (settings.plot_phase_map) {
        script << "fig, axes = plt.subplots(2, 2, figsize=(13, 10))\n";
        script << "fig.suptitle('Damped Pendulum "
               << "(theta0=" << config.physical.theta0
               << " rad, gamma=" << config.physical.gamma << ")', fontsize=14)\n\n";
        script << "ax_theta = axes[0, 0]\n";
        script << "ax_err = axes[0, 1]\n";
        script << "ax_energy = axes[1, 0]\n";
        script << "ax_phase = axes[1, 1]\n\n";

        script << "ax_theta.plot(t, theta_exact, label='" << ref_label << "', linewidth=2.2)\n";
        script << "ax_theta.plot(t, theta_nl, '-.', label='" << num_label << "', linewidth=1.8)\n";
        script << "ax_theta.set_xlabel('Time t (s)')\n";
        script << "ax_theta.set_ylabel('theta (rad)')\n";
        script << "ax_theta.set_title('Angular Displacement')\n";
        script << "ax_theta.grid(True, alpha=0.3)\n";
        script << "ax_theta.legend(loc='best')\n\n";

        script << "ax_err.semilogy(t, err_nl, label='|theta - reference|', linewidth=1.6)\n";
        script << "ax_err.set_xlabel('Time t (s)')\n";
        script << "ax_err.set_ylabel('|Error| (rad)')\n";
        script << "ax_err.set_title('Absolute Error vs Reference')\n";
        script << "ax_err.grid(True, alpha=0.3)\n";
        script << "ax_err.legend(loc='best')\n\n";

        script << "ax_energy.plot(t, energy_nl, label='E(t) nonlinear', linewidth=1.6)\n";
        script << "ax_energy.set_xlabel('Time t (s)')\n";
        script << "ax_energy.set_ylabel('E / (m L^2)')\n";
        script << "ax_energy.set_title('Mechanical Energy (nonlinear pendulum)')\n";
        script << "ax_energy.grid(True, alpha=0.3)\n";
        script << "ax_energy.legend(loc='best')\n\n";

        script << "ax_phase.plot(theta_nl, omega_nl, label='Phase trajectory', linewidth=1.4)\n";
        script << "ax_phase.plot([0.0], [0.0], 'ko', label='Equilibrium')\n";
        script << "ax_phase.set_xlabel('theta (rad)')\n";
        script << "ax_phase.set_ylabel('omega (rad/s)')\n";
        script << "ax_phase.set_title('Phase Map (theta, omega)')\n";
        script << "ax_phase.grid(True, alpha=0.3)\n";
        script << "ax_phase.legend(loc='best')\n";
        if (config.physical.damping_model == damping_force::Model::Linear) {
            script << "ax_phase.set_aspect('equal', adjustable='box')\n";
        }
        script << "\n";
    } else {
        script << "fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)\n";
        script << "fig.suptitle('Damped Pendulum "
               << "(theta0=" << config.physical.theta0
               << " rad, gamma=" << config.physical.gamma << ")', fontsize=14)\n\n";

        script << "axes[0].plot(t, theta_exact, label='" << ref_label << "', linewidth=2.2)\n";
        script << "axes[0].plot(t, theta_nl, '-.', label='" << num_label << "', linewidth=1.8)\n";
        script << "axes[0].set_ylabel('theta (rad)')\n";
        script << "axes[0].set_title('Angular Displacement')\n";
        script << "axes[0].grid(True, alpha=0.3)\n";
        script << "axes[0].legend(loc='best')\n\n";

        script << "axes[1].semilogy(t, err_nl, label='|theta - reference|', linewidth=1.6)\n";
        script << "axes[1].set_ylabel('|Error| (rad)')\n";
        script << "axes[1].set_title('Absolute Error vs Reference')\n";
        script << "axes[1].grid(True, alpha=0.3)\n";
        script << "axes[1].legend(loc='best')\n\n";

        script << "axes[2].plot(t, energy_nl, label='E(t) nonlinear', linewidth=1.6)\n";
        script << "axes[2].set_xlabel('Time t (s)')\n";
        script << "axes[2].set_ylabel('E / (m L^2)')\n";
        script << "axes[2].set_title('Mechanical Energy (nonlinear pendulum)')\n";
        script << "axes[2].grid(True, alpha=0.3)\n";
        script << "axes[2].legend(loc='best')\n\n";
    }

    script << "plt.tight_layout()\n";
    if (settings.save_png) {
        script << "plt.savefig('" << settings.output_png
               << "', dpi=150, bbox_inches='tight')\n";
        script << "print('Plot saved to " << settings.output_png << "')\n";
    }
    if (settings.show_plot) {
        script << "plt.show()\n";
    } else {
        script << "plt.close(fig)\n";
    }

    return true;
}



bool render_with_gnuplot(const DampedConfig& config,
                         const SimulationResult& result) {
    const auto& settings = config.settings;
    const auto& p = config.physical;
    const std::string ref_label = reference_label(config);
    const std::string num_label = numerical_label(config);
    const char* gnuplot_cmd = settings.show_plot ? "gnuplot -persistent" : "gnuplot";

    FILE* gp = popen(gnuplot_cmd, "w");
    if (!gp) {
        std::cerr << "Warning: could not start gnuplot.\n";
        return false;
    }

    fprintf(gp,
        "set encoding utf8\n"
        "set grid lw 0.5\n"
        "set style line 1 lc rgb '#0060ad' lw 2.5 dt 1\n"
        "set style line 2 lc rgb '#dd181f' lw 2.0 dt 2\n"
        "set style line 3 lc rgb '#228B22' lw 2.0 dt 4\n"
        "set style line 4 lc rgb '#dd181f' lw 1.5\n"
        "set style line 5 lc rgb '#228B22' lw 1.5\n"
        "set style line 6 lc rgb '#8B008B' lw 1.5\n"
        "set style line 7 lc rgb '#444444' lw 1.8\n");

    if (settings.show_plot) {
        fprintf(gp,
            "set terminal qt size 1100,850 enhanced font 'Arial,11'\n");

        if (settings.plot_phase_map) {
            fprintf(gp,
                "set multiplot layout 2,2 title "
                "'Damped Pendulum  (theta0=%.2f rad, gamma=%.2f)' font ',14'\n",
                p.theta0, p.gamma);
        } else {
            fprintf(gp,
                "set multiplot layout 3,1 title "
                "'Damped Pendulum  (theta0=%.2f rad, gamma=%.2f)' font ',14'\n",
                p.theta0, p.gamma);
        }

        fprintf(gp,
            "set xlabel 'Time t (s)'\n"
            "set ylabel 'theta (rad)'\n"
            "set title 'Angular Displacement' font ',12'\n"
            "set key top right\n"
            "plot '%s' u 1:2 w l ls 1 title '%s',"
            "     ''   u 1:3 w l ls 3 title '%s'\n",
            settings.data_file.c_str(), ref_label.c_str(), num_label.c_str());

        fprintf(gp,
            "set ylabel '|Error| (rad)'\n"
            "set title 'Absolute Error vs Reference' font ',12'\n"
            "set format y '10^{%%L}'\n"
            "set logscale y\n"
            "plot '%s' u 1:5 w l ls 4 title '|theta - reference|',"
            "     ''   u 1:6 w l ls 5 title '|omega - reference|'\n"
            "unset logscale y\n"
            "set format y '%%g'\n",
            settings.data_file.c_str());

        fprintf(gp,
            "set ylabel 'E / (m L^2)'\n"
            "set title 'Mechanical Energy (nonlinear pendulum)' font ',12'\n"
            "plot '%s' u 1:7 w l ls 6 title 'E(t) nonlinear'\n",
            settings.data_file.c_str());

        if (settings.plot_phase_map) {
            fprintf(gp,
                "set xlabel 'theta (rad)'\n"
                "set ylabel 'omega (rad/s)'\n"
                "set title 'Phase Map (theta, omega)' font ',12'\n"
                "plot '%s' u 3:4 w l ls 7 title 'Phase trajectory',"
                "     ''   u (0):(0) w p pt 7 ps 1.1 lc rgb '#000000' title 'Equilibrium'\n",
                settings.data_file.c_str());
        }

        fprintf(gp, "unset multiplot\n");
        fflush(gp);
    }

    if (settings.save_png) {
        fprintf(gp,
            "set terminal pngcairo size 1200,900 enhanced font 'Arial,12'\n"
            "set output '%s'\n",
            settings.output_png.c_str());

        if (settings.plot_phase_map) {
            fprintf(gp,
                "set multiplot layout 2,2 title "
                "'Damped Pendulum  (theta0=%.2f rad, gamma=%.2f)' font ',14'\n",
                p.theta0, p.gamma);
        } else {
            fprintf(gp,
                "set multiplot layout 3,1 title "
                "'Damped Pendulum  (theta0=%.2f rad, gamma=%.2f)' font ',14'\n",
                p.theta0, p.gamma);
        }

        fprintf(gp,
            "set xlabel 'Time t (s)'\n"
            "set ylabel 'theta (rad)'\n"
            "set title 'Angular Displacement' font ',12'\n"
            "set key top right\n"
            "plot '%s' u 1:2 w l ls 1 title '%s',"
            "     ''   u 1:3 w l ls 3 title '%s'\n",
            settings.data_file.c_str(), ref_label.c_str(), num_label.c_str());

        fprintf(gp,
            "set ylabel '|Error| (rad)'\n"
            "set title 'Absolute Error vs Reference' font ',12'\n"
            "set format y '10^{%%L}'\n"
            "set logscale y\n"
            "plot '%s' u 1:5 w l ls 4 title '|theta - reference|',"
            "     ''   u 1:6 w l ls 5 title '|omega - reference|'\n"
            "unset logscale y\n"
            "set format y '%%g'\n",
            settings.data_file.c_str());

        fprintf(gp,
            "set ylabel 'E / (m L^2)'\n"
            "set title 'Mechanical Energy (nonlinear pendulum)' font ',12'\n"
            "plot '%s' u 1:7 w l ls 6 title 'E(t) nonlinear'\n",
            settings.data_file.c_str());

        if (settings.plot_phase_map) {
            fprintf(gp,
                "set xlabel 'theta (rad)'\n"
                "set ylabel 'omega (rad/s)'\n"
                "set title 'Phase Map (theta, omega)' font ',12'\n"
                "plot '%s' u 3:4 w l ls 7 title 'Phase trajectory',"
                "     ''   u (0):(0) w p pt 7 ps 1.1 lc rgb '#000000' title 'Equilibrium'\n",
                settings.data_file.c_str());
        }

        fprintf(gp, "unset multiplot\n");
        fflush(gp);
    }

    pclose(gp);
    return true;
}

}  // namespace

void render_damped_plots(const DampedConfig& config,
                         const SimulationResult& result) {
    if (config.settings.plotting_method == PlottingMethod::Original) {
        if (render_with_gnuplot(config, result)) {
            if (config.settings.save_png) {
                std::cout << "Plot saved to " << config.settings.output_png << "\n";
            }
            return;
        }
        std::cerr << "Falling back to new plotting method.\n";
    }

    if (!write_python_plot_script(config, result)) {
        return;
    }
    std::cout << "Plotting script created: " << config.settings.python_script << "\n";

    if (config.settings.run_plotter) {
        plotting_utils::run_python_script(config.settings.python_script);
    } else {
        std::cout << "Run manually: python3 " << config.settings.python_script << "\n";
    }
}
