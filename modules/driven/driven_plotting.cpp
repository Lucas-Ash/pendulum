#include "modules/driven/driven_plotting.h"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include "modules/core/plotting_utils.h"

namespace {

void generate_plot_script(const DrivenConfig& config, const SimulationResult& result) {
    const auto& s = config.settings;
    const auto& p = config.physical;
    const auto& sim = config.simulation;

    const std::filesystem::path script_path(s.python_script);
    if (script_path.has_parent_path()) {
        std::filesystem::create_directories(script_path.parent_path());
    }
    
    std::ofstream script(s.python_script);
    if (!script.is_open()) {
        throw std::runtime_error("Could not create plotting script: " + s.python_script);
    }

    script << "#!/usr/bin/env python3\n";
    script << "\"\"\"\n";
    script << "Pendulum Simulation Visualization\n";
    script << "Compares Numerical (RK4) vs Analytical (Small-Angle) Solutions\n";
    script << "\"\"\"\n\n";
    script << "import numpy as np\n";
    script << "import matplotlib.pyplot as plt\n\n";
    script << "# Load data from CSV\n";
    script << "data = np.loadtxt('" << s.data_file << "', delimiter=',', skiprows=1)\n";
    script << "time = data[:, 0]\n";
    script << "theta_numerical = data[:, 2]\n";
    script << "omega_numerical = data[:, 3]\n";
    script << "theta_analytical = data[:, 1]\n";
    script << "difference = data[:, 4]\n\n";
    if (s.plot_phase_map) {
        script << "# Create figure with subplots\n";
        script << "fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 12))\n\n";
        script << "# Plot 1: Theta vs Time\n";
        script << "ax1.plot(time, theta_numerical, 'b-', linewidth=1.5, label='Numerical (RK4)', alpha=0.8)\n";
        script << "ax1.plot(time, theta_analytical, 'r--', linewidth=1.5, label='Analytical (Small-Angle)', alpha=0.8)\n";
        script << "ax1.set_xlabel('Time (s)', fontsize=12)\n";
        script << "ax1.set_ylabel('Angle theta (rad)', fontsize=12)\n";
        script << "ax1.set_title('Damped Driven Pendulum: Numerical vs Analytical Solution', fontsize=14, fontweight='bold')\n";
        script << "ax1.legend(fontsize=11)\n";
        script << "ax1.grid(True, alpha=0.3)\n";
        script << "ax1.set_xlim([0, " << sim.t_end << "])\n\n";
        script << "# Plot 2: Difference between solutions\n";
        script << "ax2.plot(time, difference, 'g-', linewidth=1, label='|Numerical - Analytical|')\n";
        script << "ax2.set_xlabel('Time (s)', fontsize=12)\n";
        script << "ax2.set_ylabel('Absolute Difference (rad)', fontsize=12)\n";
        script << "ax2.set_title('Difference Between Numerical and Analytical Solutions', fontsize=14, fontweight='bold')\n";
        script << "ax2.legend(fontsize=11)\n";
        script << "ax2.grid(True, alpha=0.3)\n";
        script << "ax2.set_xlim([0, " << sim.t_end << "])\n\n";
        script << "# Plot 3: Phase map\n";
        script << "ax3.plot(theta_numerical, omega_numerical, 'm-', linewidth=1.2, label='Phase trajectory')\n";
        script << "ax3.plot([0.0], [0.0], 'ko', label='Equilibrium')\n";
        script << "ax3.set_xlabel('theta (rad)', fontsize=12)\n";
        script << "ax3.set_ylabel('omega (rad/s)', fontsize=12)\n";
        script << "ax3.set_title('Phase Map (theta, omega)', fontsize=14, fontweight='bold')\n";
        script << "ax3.legend(fontsize=11)\n";
        script << "ax3.grid(True, alpha=0.3)\n";
        script << "ax3.set_aspect('equal', adjustable='box')\n\n";
    } else {
        script << "# Create figure with subplots\n";
        script << "fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))\n\n";
        script << "# Plot 1: Theta vs Time\n";
        script << "ax1.plot(time, theta_numerical, 'b-', linewidth=1.5, label='Numerical (RK4)', alpha=0.8)\n";
        script << "ax1.plot(time, theta_analytical, 'r--', linewidth=1.5, label='Analytical (Small-Angle)', alpha=0.8)\n";
        script << "ax1.set_xlabel('Time (s)', fontsize=12)\n";
        script << "ax1.set_ylabel('Angle theta (rad)', fontsize=12)\n";
        script << "ax1.set_title('Damped Driven Pendulum: Numerical vs Analytical Solution', fontsize=14, fontweight='bold')\n";
        script << "ax1.legend(fontsize=11)\n";
        script << "ax1.grid(True, alpha=0.3)\n";
        script << "ax1.set_xlim([0, " << sim.t_end << "])\n\n";
        script << "# Plot 2: Difference between solutions\n";
        script << "ax2.plot(time, difference, 'g-', linewidth=1, label='|Numerical - Analytical|')\n";
        script << "ax2.set_xlabel('Time (s)', fontsize=12)\n";
        script << "ax2.set_ylabel('Absolute Difference (rad)', fontsize=12)\n";
        script << "ax2.set_title('Difference Between Numerical and Analytical Solutions', fontsize=14, fontweight='bold')\n";
        script << "ax2.legend(fontsize=11)\n";
        script << "ax2.grid(True, alpha=0.3)\n";
        script << "ax2.set_xlim([0, " << sim.t_end << "])\n\n";
    }

    script << "# Add parameter annotation\n";
    script << "params_text = 'Parameters:\\n";
    script << "g = " << p.g << " m/s^2\\n";
    script << "L = " << p.L << " m\\n";
    script << "damping = " << p.damping << "\\n";
    script << "A = " << p.A << "\\n";
    script << "omega_drive = " << p.omega_drive << " rad/s\\n";
    script << "theta0 = " << p.theta0 << " rad\\n";
    script << "omega0 = " << p.omega0 << " rad/s'\n";
    script << "fig.text(0.02, 0.02, params_text, fontsize=9, verticalalignment='bottom',\n";
    script << "         fontfamily='monospace', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))\n\n";
    script << "plt.tight_layout()\n";
    script << "plt.subplots_adjust(bottom=0.15)\n\n";
    
    if (s.save_png) {
        script << "# Save figure\n";
        script << "plt.savefig('" << s.output_png << "', dpi=150, bbox_inches='tight')\n";
    }
    
    if (s.show_plot) {
        script << "# Display plot\n";
        script << "plt.show()\n";
    }
    script.close();
}

void generate_sweep_plot_script(const DrivenConfig& config,
                                const DrivenSweepResult& result) {
    const auto& s = config.settings;
    const char* title =
        config.physical.system_model == DrivenSystemModel::Duffing
            ? "Driven Duffing sweep"
            : "Driven pendulum sweep";

    const std::filesystem::path script_path(s.python_script);
    if (script_path.has_parent_path()) {
        std::filesystem::create_directories(script_path.parent_path());
    }

    std::ofstream script(s.python_script);
    if (!script.is_open()) {
        throw std::runtime_error("Could not create plotting script: " + s.python_script);
    }

    script << "#!/usr/bin/env python3\n";
    script << "from pathlib import Path\n";
    script << "import numpy as np\n";
    script << "import matplotlib.pyplot as plt\n\n";
    script << "root = Path(__file__).resolve().parents[3]\n";
    script << "data = np.loadtxt(root / '" << s.sweep_data_file
           << "', delimiter=',', skiprows=1)\n";
    script << "data = np.atleast_2d(data)\n";
    script << "freq = data[:, 0]\n";
    script << "num = data[:, 1]\n";
    script << "ana = data[:, 2]\n";
    script << "low = data[:, 3]\n";
    script << "high = data[:, 4]\n\n";
    script << "fig, ax = plt.subplots(figsize=(10, 6))\n";
    script << "ax.plot(freq, num, marker='o', markersize=3.0, linewidth=1.8, "
              "label='Numerical amplitude')\n";
    script << "if np.any(ana > 0):\n";
    script << "    ax.plot(freq, ana, '--', linewidth=1.4, label='Analytical amplitude')\n";
    script << "if np.any(high > low):\n";
    script << "    ax.plot(freq, low, ':', linewidth=1.0, label='Analytical lower stable')\n";
    script << "    ax.plot(freq, high, ':', linewidth=1.0, label='Analytical upper stable')\n";
    script << "ax.set_xlabel('Drive frequency (rad/s)')\n";
    script << "ax.set_ylabel('Steady-state amplitude')\n";
    script << "ax.set_title('" << title << "')\n";
    script << "ax.grid(True, alpha=0.3)\n";
    script << "ax.legend()\n";
    script << "plt.tight_layout()\n";
    if (s.save_png) {
        script << "plt.savefig(root / '" << s.output_png << "', dpi=150, bbox_inches='tight')\n";
    }
    if (s.show_plot) {
        script << "plt.show()\n";
    } else {
        script << "plt.close(fig)\n";
    }
    script.close();
    (void)result;
}

} // namespace

void render_driven_plots(const DrivenConfig& config, const SimulationResult& result) {
    if (config.settings.plotting_method == DrivenPlottingMethod::New) {
        generate_plot_script(config, result);
        
        if (config.settings.run_plotter) {
            plotting_utils::run_python_script(config.settings.python_script);
        }
    } else {
        std::cerr << "Warning: 'original' plotting method (gnuplot) is not implemented for driven pendulum.\n";
    }
}

void render_driven_sweep_plots(const DrivenConfig& config, const DrivenSweepResult& result) {
    if (config.settings.plotting_method == DrivenPlottingMethod::New) {
        generate_sweep_plot_script(config, result);
        if (config.settings.run_plotter) {
            plotting_utils::run_python_script(config.settings.python_script);
        }
    } else {
        std::cerr << "Warning: 'original' plotting method (gnuplot) is not implemented for driven sweeps.\n";
    }
}
