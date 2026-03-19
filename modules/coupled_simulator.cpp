#include "modules/coupled_simulator.h"
#include "modules/rk_integrators.h"
#include <cmath>
#include <fstream>
#include <iostream>

CoupledSimulator::CoupledSimulator(const CoupledConfig& config)
    : config_(config) {}

CoupledSimulationResult CoupledSimulator::simulate() const {
    CoupledSimulationResult result;
    const auto& c = config_.physical;
    const auto& sim = config_.simulation;

    int nsteps = static_cast<int>((sim.t_end - sim.t_start) / sim.dt + 0.5);
    int output_every = sim.output_every > 0 ? sim.output_every : 1;

    result.t.reserve(nsteps / output_every + 1);
    result.q1.reserve(result.t.capacity());
    result.omega1.reserve(result.t.capacity());
    result.q2.reserve(result.t.capacity());
    result.omega2.reserve(result.t.capacity());
    
    result.q1_reference.reserve(result.t.capacity());
    result.omega1_reference.reserve(result.t.capacity());
    result.q2_reference.reserve(result.t.capacity());
    result.omega2_reference.reserve(result.t.capacity());

    auto derivatives = [&c](double t, const CoupledState& state) -> CoupledState {
        double dq1_dt = state.omega1;
        double dq2_dt = state.omega2;

        double domega1_dt = c.F * std::cos(c.omega_drive * t)
                            - 2.0 * c.zeta_1 * c.omega_1 * state.omega1
                            - c.omega_1 * c.omega_1 * state.q1
                            - c.alpha_11 * state.q1 * state.q1 * state.q1
                            - c.alpha_12 * state.q1 * state.q2 * state.q2;

        double domega2_dt = - 2.0 * c.zeta_2 * c.omega_2 * state.omega2
                            - c.omega_2 * c.omega_2 * state.q2
                            - c.alpha_22 * state.q2 * state.q2 * state.q2
                            - c.alpha_21 * state.q1 * state.q1 * state.q2;

        return {dq1_dt, domega1_dt, dq2_dt, domega2_dt};
    };

    CoupledState state = {c.q1_0, c.omega1_0, c.q2_0, c.omega2_0};
    CoupledState state_ref = state;
    double t = sim.t_start;
    double t_ref = sim.t_start;
    double dt_ref = sim.dt / static_cast<double>(config_.settings.error_reference_factor);

    auto step_forward = [this, &derivatives](double curr_t, double curr_dt, const CoupledState& curr_state) {
        if (config_.settings.integrator == "rk4") return integrator::rk4_step(curr_t, curr_dt, curr_state, derivatives);
        if (config_.settings.integrator == "rk5") return integrator::rk5_step(curr_t, curr_dt, curr_state, derivatives);
        if (config_.settings.integrator == "rk3") return integrator::rk3_step(curr_t, curr_dt, curr_state, derivatives);
        if (config_.settings.integrator == "rk23") return integrator::rk23_step(curr_t, curr_dt, curr_state, derivatives);
        if (config_.settings.integrator == "rkf45") return integrator::rkf45_step(curr_t, curr_dt, curr_state, derivatives);
        return integrator::rk4_step(curr_t, curr_dt, curr_state, derivatives); // fallback
    };

    for (int i = 0; i <= nsteps; ++i) {
        if (i % output_every == 0) {
            result.t.push_back(t);
            result.q1.push_back(state.q1);
            result.omega1.push_back(state.omega1);
            result.q2.push_back(state.q2);
            result.omega2.push_back(state.omega2);
            
            double q1_ref_val = state.q1;
            double omega1_ref_val = state.omega1;
            double q2_ref_val = state.q2;
            double omega2_ref_val = state.omega2;

            if (config_.settings.error_mode == error_reference::Mode::HdReference) {
                q1_ref_val = state_ref.q1;
                omega1_ref_val = state_ref.omega1;
                q2_ref_val = state_ref.q2;
                omega2_ref_val = state_ref.omega2;
            }
            
            result.q1_reference.push_back(q1_ref_val);
            result.omega1_reference.push_back(omega1_ref_val);
            result.q2_reference.push_back(q2_ref_val);
            result.omega2_reference.push_back(omega2_ref_val);
        }

        if (i < nsteps) {
            state = step_forward(t, sim.dt, state);
            if (config_.settings.error_mode == error_reference::Mode::HdReference) {
                for (int k = 0; k < config_.settings.error_reference_factor; ++k) {
                    state_ref = step_forward(t_ref, dt_ref, state_ref);
                    t_ref += dt_ref;
                }
            }
            t += sim.dt;
        }
    }

    return result;
}

void write_coupled_data_file(const std::string& path, const CoupledSimulationResult& result) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Could not open data file for writing: " + path);
    }
    out << "t,q1,omega1,q2,omega2,q1_ref,omega1_ref,q2_ref,omega2_ref\n";
    for (size_t i = 0; i < result.t.size(); ++i) {
        out << result.t[i] << "," << result.q1[i] << "," << result.omega1[i] 
            << "," << result.q2[i] << "," << result.omega2[i] << ","
            << result.q1_reference[i] << "," << result.omega1_reference[i] << ","
            << result.q2_reference[i] << "," << result.omega2_reference[i] << "\n";
    }
}

void write_coupled_plot_script(const std::string& script_path, const std::string& data_file, const std::string& output_png) {
    std::ofstream script(script_path);
    if (!script) {
        throw std::runtime_error("Could not write coupled python script: " + script_path);
    }
    script << "import pandas as pd\n";
    script << "import matplotlib.pyplot as plt\n";
    script << "import sys\n\n";
    script << "def main():\n";
    script << "    try:\n";
    script << "        df = pd.read_csv('" << data_file << "')\n";
    script << "    except Exception as e:\n";
    script << "        print(f'Error reading data: {e}')\n";
    script << "        sys.exit(1)\n\n";
    script << "    has_error = 'q1_ref' in df.columns and (df['q1'] != df['q1_ref']).any()\n";
    script << "    fig, axes = plt.subplots(3 if has_error else 2, 1, figsize=(10, 12 if has_error else 8))\n";
    script << "    if not has_error:\n";
    script << "        ax1, ax2 = axes\n";
    script << "    else:\n";
    script << "        ax1, ax2, ax3 = axes\n\n";
    script << "    ax1.plot(df['t'], df['q1'], label='q1')\n";
    script << "    ax1.plot(df['t'], df['q2'], label='q2')\n";
    script << "    ax1.set_title('Time Series: q1 and q2')\n";
    script << "    ax1.set_xlabel('Time (t)')\n";
    script << "    ax1.set_ylabel('Displacement')\n";
    script << "    ax1.legend()\n";
    script << "    ax1.grid(True)\n\n";
    script << "    ax2.plot(df['q1'], df['q2'], color='purple')\n";
    script << "    ax2.set_title('Phase Space: q1 vs q2')\n";
    script << "    ax2.set_xlabel('q1')\n";
    script << "    ax2.set_ylabel('q2')\n";
    script << "    ax2.grid(True)\n\n";
    script << "    if has_error:\n";
    script << "        ax3.plot(df['t'], abs(df['q1_ref'] - df['q1']), label='Error q1')\n";
    script << "        ax3.plot(df['t'], abs(df['q2_ref'] - df['q2']), label='Error q2')\n";
    script << "        ax3.set_title('Absolute Error vs HD Reference')\n";
    script << "        ax3.set_xlabel('Time (t)')\n";
    script << "        ax3.set_ylabel('Absolute Error')\n";
    script << "        ax3.legend()\n";
    script << "        ax3.grid(True)\n\n";
    script << "    plt.tight_layout()\n";
    script << "    plt.savefig('" << output_png << "')\n";
    script << "    print('Saved plot to " << output_png << "')\n\n";
    script << "if __name__ == '__main__':\n";
    script << "    main()\n";
}
