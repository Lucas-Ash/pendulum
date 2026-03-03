#include "pendulum/damped_plotting.h"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

namespace {

bool write_python_plot_script(const DampedConfig& config,
                              const DampedSimulationResult& result) {
    const auto& settings = config.settings;

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
    script << "theta_lin = data[:,2]\n";
    script << "theta_nl = data[:,3]\n";
    script << "err_lin = data[:,4]\n";
    script << "err_nl = data[:,5]\n";
    script << "energy_nl = data[:,6]\n\n";
    script << "fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)\n";
    script << "fig.suptitle('Damped Pendulum "
           << "(theta0=" << config.physical.theta0
           << " rad, gamma=" << config.physical.gamma
           << ", omega0=" << result.omega0 << ")', fontsize=14)\n\n";

    script << "axes[0].plot(t, theta_exact, label='Analytical (linearized)', linewidth=2.2)\n";
    script << "axes[0].plot(t, theta_lin, '--', label='RK4 linearized', linewidth=1.8)\n";
    script << "axes[0].plot(t, theta_nl, '-.', label='RK4 nonlinear (sin theta)', linewidth=1.8)\n";
    script << "axes[0].set_ylabel('theta (rad)')\n";
    script << "axes[0].set_title('Angular Displacement')\n";
    script << "axes[0].grid(True, alpha=0.3)\n";
    script << "axes[0].legend(loc='best')\n\n";

    script << "axes[1].semilogy(t, err_lin, label='|RK4 linear - analytical|', linewidth=1.6)\n";
    script << "axes[1].semilogy(t, err_nl, label='|RK4 nonlinear - analytical|', linewidth=1.6)\n";
    script << "axes[1].set_ylabel('|Error| (rad)')\n";
    script << "axes[1].set_title('Absolute Error vs Analytical Solution')\n";
    script << "axes[1].grid(True, alpha=0.3)\n";
    script << "axes[1].legend(loc='best')\n\n";

    script << "axes[2].plot(t, energy_nl, label='E(t) nonlinear', linewidth=1.6)\n";
    script << "axes[2].set_xlabel('Time t (s)')\n";
    script << "axes[2].set_ylabel('E / (m L^2)')\n";
    script << "axes[2].set_title('Mechanical Energy (nonlinear pendulum)')\n";
    script << "axes[2].grid(True, alpha=0.3)\n";
    script << "axes[2].legend(loc='best')\n\n";

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

bool run_python_plotter(const DampedConfig& config) {
    const std::string command = "python3 \"" + config.settings.python_script + "\"";
    const int status = std::system(command.c_str());
    if (status != 0) {
        std::cerr << "Warning: python plotting command failed: " << command << "\n";
        return false;
    }
    return true;
}

bool render_with_gnuplot(const DampedConfig& config,
                         const DampedSimulationResult& result) {
    const auto& settings = config.settings;
    const auto& p = config.physical;
    const char* gnuplot_cmd = settings.show_plot ? "gnuplot -persistent" : "gnuplot";

    FILE* gp = popen(gnuplot_cmd, "w");
    if (!gp) {
        std::cerr << "Warning: could not start gnuplot.\n";
        return false;
    }

    if (settings.show_plot) {
        fprintf(gp,
            "set terminal qt size 1100,850 enhanced font 'Arial,11'\n"
            "set encoding utf8\n"
            "set grid lw 0.5\n"
            "set style line 1 lc rgb '#0060ad' lw 2.5 dt 1\n"
            "set style line 2 lc rgb '#dd181f' lw 2.0 dt 2\n"
            "set style line 3 lc rgb '#228B22' lw 2.0 dt 4\n"
            "set style line 4 lc rgb '#dd181f' lw 1.5\n"
            "set style line 5 lc rgb '#228B22' lw 1.5\n"
            "set style line 6 lc rgb '#8B008B' lw 1.5\n");

        fprintf(gp,
            "set multiplot layout 3,1 title "
            "'Damped Pendulum  (theta0=%.2f rad, gamma=%.2f, omega0=%.4f)' font ',14'\n",
            p.theta0, p.gamma, result.omega0);

        fprintf(gp,
            "set xlabel 'Time t (s)'\n"
            "set ylabel 'theta (rad)'\n"
            "set title 'Angular Displacement' font ',12'\n"
            "set key top right\n"
            "plot '%s' u 1:2 w l ls 1 title 'Analytical (linearized)',"
            "     ''   u 1:3 w l ls 2 title 'RK4 linearized',"
            "     ''   u 1:4 w l ls 3 title 'RK4 nonlinear (sin theta)'\n",
            settings.data_file.c_str());

        fprintf(gp,
            "set ylabel '|Error| (rad)'\n"
            "set title 'Absolute Error vs Analytical Solution' font ',12'\n"
            "set format y '10^{%%L}'\n"
            "set logscale y\n"
            "plot '%s' u 1:5 w l ls 4 title '|RK4 linear - analytical|',"
            "     ''   u 1:6 w l ls 5 title '|RK4 nonlinear - analytical|'\n"
            "unset logscale y\n"
            "set format y '%%g'\n",
            settings.data_file.c_str());

        fprintf(gp,
            "set ylabel 'E / (m L^2)'\n"
            "set title 'Mechanical Energy (nonlinear pendulum)' font ',12'\n"
            "plot '%s' u 1:7 w l ls 6 title 'E(t) nonlinear'\n",
            settings.data_file.c_str());

        fprintf(gp, "unset multiplot\n");
        fflush(gp);
    }

    if (settings.save_png) {
        fprintf(gp,
            "set terminal pngcairo size 1200,900 enhanced font 'Arial,12'\n"
            "set output '%s'\n",
            settings.output_png.c_str());

        fprintf(gp,
            "set multiplot layout 3,1 title "
            "'Damped Pendulum  (theta0=%.2f rad, gamma=%.2f, omega0=%.4f)' font ',14'\n",
            p.theta0, p.gamma, result.omega0);

        fprintf(gp,
            "set xlabel 'Time t (s)'\n"
            "set ylabel 'theta (rad)'\n"
            "set title 'Angular Displacement' font ',12'\n"
            "set key top right\n"
            "plot '%s' u 1:2 w l ls 1 title 'Analytical (linearized)',"
            "     ''   u 1:3 w l ls 2 title 'RK4 linearized',"
            "     ''   u 1:4 w l ls 3 title 'RK4 nonlinear (sin theta)'\n",
            settings.data_file.c_str());

        fprintf(gp,
            "set ylabel '|Error| (rad)'\n"
            "set title 'Absolute Error vs Analytical Solution' font ',12'\n"
            "set format y '10^{%%L}'\n"
            "set logscale y\n"
            "plot '%s' u 1:5 w l ls 4 title '|RK4 linear - analytical|',"
            "     ''   u 1:6 w l ls 5 title '|RK4 nonlinear - analytical|'\n"
            "unset logscale y\n"
            "set format y '%%g'\n",
            settings.data_file.c_str());

        fprintf(gp,
            "set ylabel 'E / (m L^2)'\n"
            "set title 'Mechanical Energy (nonlinear pendulum)' font ',12'\n"
            "plot '%s' u 1:7 w l ls 6 title 'E(t) nonlinear'\n",
            settings.data_file.c_str());

        fprintf(gp, "unset multiplot\n");
        fflush(gp);
    }

    pclose(gp);
    return true;
}

}  // namespace

void render_damped_plots(const DampedConfig& config,
                         const DampedSimulationResult& result) {
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
        run_python_plotter(config);
    } else {
        std::cout << "Run manually: python3 " << config.settings.python_script << "\n";
    }
}
