#!/usr/bin/env python3
"""
plot_poincare.py
Generates a Poincare section for the driven pendulum simulation by sampling
the state (theta, omega) exactly once every driving period T = 2*pi / omega_drive.

Usage:
  python scripts/plot_poincare.py [--omega-drive OMEGA] [--t-end T_END] [--base-config YAML]
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import tempfile
import yaml
import os
from pathlib import Path

def run_simulation(exec_path: Path, config_path: Path, cwd: Path) -> subprocess.CompletedProcess:
    env = os.environ.copy()
    env["QA_TEST"] = "1"
    return subprocess.run(
        [str(exec_path.resolve()), str(config_path.resolve())],
        cwd=str(cwd),
        env=env,
        capture_output=True,
        text=True,
        timeout=120,
    )

def normalize_angle(theta):
    """Normalize angle to [-pi, pi] interval."""
    return (theta + np.pi) % (2 * np.pi) - np.pi

def plot_poincare(omega_drive, t_end, base_config_path, project_root):
    exec_path = project_root / "integrated_pendulum"
    if not exec_path.exists():
        print(f"Error: Executable not found at {exec_path}")
        return

    # Load base config
    with open(base_config_path, 'r') as f:
        base_config = yaml.safe_load(f)
    
    # We don't want the simulation to plot
    if 'settings' in base_config:
        base_config['settings']['run_plotter'] = False
        base_config['settings']['show_plot'] = False

    if omega_drive <= 0:
        print("Error: omega_drive must be strictly positive.")
        return

    print(f"Generating Poincare section for omega_drive = {omega_drive}")
    print(f"Running simulation until t = {t_end}")
    
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir_path = Path(temp_dir)
        
        # Modify config
        config = yaml.safe_load(yaml.dump(base_config))
        config['physical']['omega_drive'] = float(omega_drive)
        config['simulation']['t_end'] = float(t_end)
        
        # We need a small dt to accurately interpolate periods
        if float(config['simulation'].get('dt', 0.01)) > 0.01:
            config['simulation']['dt'] = 0.005
            
        # Optional: increase output_every to not create massive files, but keep it small enough for interpolation
        config['simulation']['output_every'] = 1 
        
        data_file = temp_dir_path / "poincare_data.csv"
        config['settings']['data_file'] = str(data_file.resolve())
        
        # Write temp config
        temp_config_path = temp_dir_path / "poincare_config.yaml"
        with open(temp_config_path, 'w') as f:
            yaml.dump(config, f)
        
        # Run simulation
        print("Running simulation... this may take a moment.")
        proc = run_simulation(exec_path, temp_config_path, project_root)
        if proc.returncode != 0:
            print("Simulation FAILED")
            print(proc.stderr)
            return
            
        if not data_file.exists():
            print("Error: Output data file not generated.")
            return
            
        print("Processing data...")
        try:
            data = np.loadtxt(data_file, delimiter=',', skiprows=1)
            time = data[:, 0]
            theta = data[:, 2] # numerical theta
            omega = data[:, 3] # numerical omega
            
            # Driving period
            T = 2 * np.pi / omega_drive
            
            # Find times that are multiples of T
            # (skip the first few periods to drop transients)
            max_period = int(np.floor(time[-1] / T))
            
            # Discard first 20% of periods
            start_period = int(max_period * 0.2)
            if max_period - start_period < 10:
                print("Warning: Might not be enough periods for a good Poincare section.")
                start_period = 0
            
            poincare_theta = []
            poincare_omega = []
            
            for p in range(start_period, max_period + 1):
                target_t = p * T
                # Find closest index
                idx = np.searchsorted(time, target_t)
                
                # Check bounds
                if idx == 0:
                    best_idx = 0
                elif idx == len(time):
                    best_idx = len(time) - 1
                else:
                    # Linear interpolation would be better, but closest point is simpler and fine for small dt
                    if (time[idx] - target_t) < (target_t - time[idx-1]):
                        best_idx = idx
                    else:
                        best_idx = idx - 1
                
                # Interpolation approach for slightly more accuracy (optional):
                # We can linearly interpolate theta and omega at exactly target_t
                if 0 < idx < len(time):
                    t0, t1 = time[idx-1], time[idx]
                    th0, th1 = theta[idx-1], theta[idx]
                    om0, om1 = omega[idx-1], omega[idx]
                    
                    factor = (target_t - t0) / (t1 - t0)
                    interp_th = th0 + factor * (th1 - th0)
                    interp_om = om0 + factor * (om1 - om0)
                    
                    poincare_theta.append(normalize_angle(interp_th))
                    poincare_omega.append(interp_om)

            print(f"Extracted {len(poincare_theta)} points for Poincare section.")
            
            # Plot
            plt.figure(figsize=(8, 8))
            plt.scatter(poincare_theta, poincare_omega, s=2, c='k', alpha=0.5)
            plt.xlabel('$\\theta$ (mod $2\\pi$) (rad)', fontsize=14)
            plt.ylabel('$\\omega$ (rad/s)', fontsize=14)
            plt.title(f'Poincare Section ($\\omega_{{drive}} = {omega_drive:.3f}$, $A = {config["physical"].get("A", "unknown")}$)', fontsize=16)
            plt.grid(True, alpha=0.3)
            plt.xlim([-np.pi, np.pi])
            
            # Save
            out_dir = project_root / "outputs"
            out_dir.mkdir(exist_ok=True)
            out_file = out_dir / f"poincare_section_w{omega_drive:.3f}.png"
            
            plt.tight_layout()
            plt.savefig(out_file, dpi=150)
            print(f"Poincare section plot saved to {out_file}")
            
        except Exception as e:
            print(f"Error analyzing data: {e}")

def main():
    parser = argparse.ArgumentParser(description="Generate Poincare section for driven pendulum.")
    parser.add_argument("--omega-drive", type=float, default=0.66667, help="Driving frequency (rad/s)")
    parser.add_argument("--t-end", type=float, default=1500.0, help="Simulation end time (longer = more points)")
    parser.add_argument("--base-config", type=str, default="driven_pendulum.yaml", help="Path to base YAML config file")
    
    args = parser.parse_args()
    
    project_root = Path(__file__).resolve().parent.parent
    base_config_path = (project_root / args.base_config).resolve()
    
    if not base_config_path.exists():
        print(f"Error: Base config file not found: {base_config_path}")
        return
        
    plot_poincare(args.omega_drive, args.t_end, base_config_path, project_root)

if __name__ == "__main__":
    main()
