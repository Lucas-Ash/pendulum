#include "modules/driven_plotting.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace {

void generate_plot_script(const DrivenConfig& config, const SimulationResult& result) {
    const auto& s = config.settings;
    const auto& p = config.physical;
    const auto& sim = config.simulation;
    
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
    script << "theta_analytical = data[:, 1]\n";
    script << "difference = data[:, 4]\n\n";
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

} // namespace

void render_driven_plots(const DrivenConfig& config, const SimulationResult& result) {
    if (config.settings.plotting_method == DrivenPlottingMethod::New) {
        generate_plot_script(config, result);
        
        if (config.settings.run_plotter) {
            std::string cmd = "python3 " + config.settings.python_script;
            int ret = std::system(cmd.c_str());
            if (ret != 0) {
                std::cerr << "Warning: Python plotting script returned " << ret << "\n";
            }
        }
    } else {
        std::cerr << "Warning: 'original' plotting method (gnuplot) is not implemented for driven pendulum.\n";
    }
}
