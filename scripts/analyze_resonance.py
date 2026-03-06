#!/usr/bin/env python3
"""
analyze_resonance.py
Sweeps across driving frequencies for the driven pendulum, running the simulation
for each, and plots the maximum steady-state amplitude vs. driving frequency.

Usage:
  python scripts/analyze_resonance.py [--omega-min MIN] [--omega-max MAX] [--steps N] [--base-config YAML]
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import tempfile
import yaml
import os
import shutil
from pathlib import Path

def run_simulation(exec_path: Path, config_path: Path, cwd: Path) -> subprocess.CompletedProcess:
    env = os.environ.copy()
    env["QA_TEST"] = "1"
    # Disable run_plotter in case it wasn't disabled in the modified YAML
    return subprocess.run(
        [str(exec_path.resolve()), str(config_path.resolve())],
        cwd=str(cwd),
        env=env,
        capture_output=True,
        text=True,
        timeout=120,
    )

def analyze_resonance(omega_min, omega_max, steps, base_config_path, project_root):
    exec_path = project_root / "integrated_pendulum"
    if not exec_path.exists():
        print(f"Error: Executable not found at {exec_path}")
        return

    # Load base config
    with open(base_config_path, 'r') as f:
        base_config = yaml.safe_load(f)
    
    # We don't want the simulation to plot during the sweep
    if 'settings' in base_config:
        base_config['settings']['run_plotter'] = False
        base_config['settings']['show_plot'] = False

    try:
        g = float(base_config.get('physical', {}).get('g', 9.81))
        L = float(base_config.get('physical', {}).get('L', 1.0))
        omega0 = np.sqrt(g / L)
    except Exception:
        omega0 = None

    if omega0 is not None and omega_min < omega0 < omega_max:
        # Concentrate points around theoretical natural frequency (omega0) using cubic spacing
        n_lower = max(2, int(steps * (omega0 - omega_min) / (omega_max - omega_min)))
        n_upper = max(2, steps - n_lower)
        
        t_lower = np.linspace(1, 0, n_lower, endpoint=False)
        omegas_lower = omega0 - (omega0 - omega_min) * (t_lower ** 3)
        
        t_upper = np.linspace(0, 1, n_upper)
        omegas_upper = omega0 + (omega_max - omega0) * (t_upper ** 3)
        
        omegas = np.concatenate([omegas_lower, omegas_upper])
    else:
        omegas = np.linspace(omega_min, omega_max, steps)

    amplitudes = []
    analytical_amplitudes = []

    print(f"Sweeping omega_drive from {omega_min} to {omega_max} in {steps} steps...")
    
    # Create a temporary directory for config files and outputs
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir_path = Path(temp_dir)
        
        for i, omega in enumerate(omegas):
            print(f"  Step {i+1}/{steps}: omega_drive = {omega:.3f} ...", end=" ", flush=True)
            
            # Modify config
            config = yaml.safe_load(yaml.dump(base_config)) # Deep copy basically
            config['physical']['omega_drive'] = float(omega)
            
            # Direct output to temp dir
            data_file = temp_dir_path / f"data_{i}.csv"
            config['settings']['data_file'] = str(data_file.resolve())
            
            # Write temp config
            temp_config_path = temp_dir_path / f"config_{i}.yaml"
            with open(temp_config_path, 'w') as f:
                yaml.dump(config, f)
            
            # Run simulation
            proc = run_simulation(exec_path, temp_config_path, project_root)
            if proc.returncode != 0:
                print("FAILED")
                print(proc.stderr)
                amplitudes.append(0.0)
                analytical_amplitudes.append(0.0)
                continue
            
            # Analyze output
            if not data_file.exists():
                print("NO DATA FILE")
                amplitudes.append(0.0)
                analytical_amplitudes.append(0.0)
                continue
            
            try:
                data = np.loadtxt(data_file, delimiter=',', skiprows=1)
                
                # data shape is (N, 5): time, theta_analytical, theta_numerical, omega_numerical, diff
                times = data[:, 0]
                thetas = data[:, 2] # numerical theta
                theta_analytical = data[:, 1] # analytical theta
                
                # Check for NaNs or infs (instability)
                if np.any(np.isnan(thetas)) or np.any(np.isinf(thetas)):
                    print("UNSTABLE (NaN/Inf)")
                    amplitudes.append(np.nan)
                    analytical_amplitudes.append(np.nan)
                    continue
                
                # Discard the first 20% of data to remove transients (steady state)
                n_points = len(thetas)
                steady_state_idx = int(0.2 * n_points)
                
                if steady_state_idx >= n_points:
                    print("NOT ENOUGH DATA")
                    amplitudes.append(0.0)
                    analytical_amplitudes.append(0.0)
                    continue
                
                steady_thetas = thetas[steady_state_idx:]
                steady_thetas_analytical = theta_analytical[steady_state_idx:]
                
                # Maximum amplitude in steady state
                max_amp = np.max(np.abs(steady_thetas))
                amplitudes.append(max_amp)
                
                max_amp_analytical = np.max(np.abs(steady_thetas_analytical))
                analytical_amplitudes.append(max_amp_analytical)
                
                print(f"Amp = {max_amp:.4f} (Analytic = {max_amp_analytical:.4f})")
            except Exception as e:
                print(f"ERROR processing data: {e}")
                amplitudes.append(0.0)
                analytical_amplitudes.append(0.0)

    # Plot results
    plt.figure(figsize=(10, 6))
    plt.plot(omegas, amplitudes, 'bo-', linewidth=2, markersize=5, label='Numerical Maximum')
    plt.plot(omegas, analytical_amplitudes, 'r--', linewidth=2, label='Analytical Maximum')
    plt.xlabel('Driving Frequency $\\omega_{drive}$ (rad/s)', fontsize=14)
    plt.ylabel('Maximum Steady-State Amplitude $\\theta_{max}$ (rad)', fontsize=14)
    plt.title('Resonance Curve for Driven Pendulum', fontsize=16)
    plt.grid(True)
    
    # Add vertical line at natural frequency (sqrt(g/L)) if possible
    try:
        g = float(base_config.get('physical', {}).get('g', 9.81))
        L = float(base_config.get('physical', {}).get('L', 1.0))
        omega0 = np.sqrt(g / L)
        plt.axvline(x=omega0, color='g', linestyle=':', label=f'Natural Freq $\\omega_0$ \\approx {omega0:.2f}')
    except:
        pass
        
    plt.legend()

    out_dir = project_root / "outputs"
    out_dir.mkdir(exist_ok=True)
    out_file = out_dir / "resonance_plot.png"
    
    plt.tight_layout()
    plt.savefig(out_file, dpi=150)
    print(f"\nResonance plot saved to {out_file}")

def main():
    parser = argparse.ArgumentParser(description="Analyze resonance of driven pendulum.")
    parser.add_argument("--omega-min", type=float, default=0.5, help="Minimum driving frequency (rad/s)")
    parser.add_argument("--omega-max", type=float, default=5.0, help="Maximum driving frequency (rad/s)")
    parser.add_argument("--steps", type=int, default=80, help="Number of frequency steps")
    parser.add_argument("--base-config", type=str, default="driven_pendulum.yaml", help="Path to base YAML config file")
    
    args = parser.parse_args()
    
    project_root = Path(__file__).resolve().parent.parent
    base_config_path = (project_root / args.base_config).resolve()
    
    if not base_config_path.exists():
        print(f"Error: Base config file not found: {base_config_path}")
        return
        
    analyze_resonance(args.omega_min, args.omega_max, args.steps, base_config_path, project_root)

if __name__ == "__main__":
    main()
