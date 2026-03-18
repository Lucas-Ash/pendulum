#!/usr/bin/env python3
"""
Build a resonance curve for the driven pendulum and compare it against the
linear damped-driven analytical response.

This script supports two workflows:
1. Native sweep mode: run one YAML with ``sweep.enabled: true`` and analyze the
   generated sweep CSV.
2. Legacy point-by-point mode: vary ``omega_drive`` across repeated single runs.

In both cases the script writes:
- ``outputs/<prefix>.csv``: resonance data by drive frequency
- ``outputs/<prefix>_drift.csv``: drift/summary metrics versus analytical theory
- ``outputs/<prefix>.png``: comparison plot
"""

import argparse
import csv
import math
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

def nested_or_top_level(config: dict, section: str, key: str, default):
    section_value = config.get(section, {})
    if isinstance(section_value, dict) and key in section_value:
        return section_value[key]
    return config.get(key, default)


def restoring_linearized_slope(config: dict) -> float:
    model = str(
        nested_or_top_level(config, "physical", "restoring_force_model", "sine")
    ).strip().lower()
    if model in {"polynomial", "duffing"}:
        return float(nested_or_top_level(config, "physical", "restoring_force_linear", 1.0))
    return 1.0


def driven_parameters(config: dict) -> tuple[float, float, float, float]:
    g = float(nested_or_top_level(config, "physical", "g", 9.81))
    length = float(nested_or_top_level(config, "physical", "L", 1.0))
    damping = float(nested_or_top_level(config, "physical", "damping", 0.0))
    drive_amplitude = float(nested_or_top_level(config, "physical", "A", 0.0))
    return g, length, damping, drive_amplitude


def linearized_omega0_sq(config: dict) -> float:
    g, length, _, _ = driven_parameters(config)
    return (g / length) * restoring_linearized_slope(config)


def analytical_amplitude(config: dict, drive_frequency: float) -> float:
    omega0_sq = linearized_omega0_sq(config)
    _, _, damping, drive_amplitude = driven_parameters(config)
    denominator = math.sqrt(
        (omega0_sq - drive_frequency * drive_frequency) ** 2 +
        (damping * drive_frequency) ** 2
    )
    if denominator <= 0.0:
        return math.inf
    return drive_amplitude / denominator


def theoretical_peak_frequency(config: dict) -> float:
    omega0_sq = linearized_omega0_sq(config)
    _, _, damping, _ = driven_parameters(config)
    peak_sq = omega0_sq - 0.5 * damping * damping
    return math.sqrt(max(0.0, peak_sq))


def build_omega_grid(omega_min: float, omega_max: float, steps: int, omega0: float | None) -> np.ndarray:
    if omega0 is not None and omega_min < omega0 < omega_max:
        n_lower = max(2, int(steps * (omega0 - omega_min) / (omega_max - omega_min)))
        n_upper = max(2, steps - n_lower)

        t_lower = np.linspace(1.0, 0.0, n_lower, endpoint=False)
        omegas_lower = omega0 - (omega0 - omega_min) * (t_lower ** 3)

        t_upper = np.linspace(0.0, 1.0, n_upper)
        omegas_upper = omega0 + (omega_max - omega0) * (t_upper ** 3)

        return np.concatenate([omegas_lower, omegas_upper])
    return np.linspace(omega_min, omega_max, steps)


def steady_state_amplitude(times: np.ndarray, theta_values: np.ndarray, drive_frequency: float) -> float:
    if len(theta_values) == 0:
        return 0.0
    start = int(0.65 * len(theta_values))
    if start >= len(theta_values):
        return 0.0

    window_start = start
    if len(times) >= 2 and abs(drive_frequency) > 1e-15:
        dt = float(times[1] - times[0])
        period = 2.0 * math.pi / abs(drive_frequency)
        if dt > 0.0 and period > dt:
            samples_per_period = max(1, int(round(period / dt)))
            available_samples = len(theta_values) - start
            full_periods = available_samples // samples_per_period
            if full_periods >= 1:
                measurement_count = full_periods * samples_per_period
                window_start = len(theta_values) - measurement_count

    steady_times = times[window_start:]
    steady_theta = theta_values[window_start:]
    mean = float(np.mean(steady_theta))
    centered = steady_theta - mean
    phases = drive_frequency * steady_times
    cosine_projection = float(np.sum(centered * np.cos(phases)))
    sine_projection = float(np.sum(centered * np.sin(phases)))
    scale = 2.0 / len(steady_theta)
    a = scale * cosine_projection
    b = scale * sine_projection
    return math.sqrt(a * a + b * b)


def write_csv(path: Path, header: list[str], rows: list[list[float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def run_native_sweep(base_config: dict,
                     base_config_path: Path,
                     exec_path: Path,
                     project_root: Path) -> tuple[np.ndarray, np.ndarray]:
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir_path = Path(temp_dir)
        config = yaml.safe_load(yaml.dump(base_config))
        settings = config.setdefault("settings", {})
        settings["run_plotter"] = False
        settings["show_plot"] = False
        settings["save_png"] = False
        sweep_csv = temp_dir_path / "native_sweep.csv"
        settings["sweep_data_file"] = str(sweep_csv.resolve())
        settings["python_script"] = str((temp_dir_path / "plot_sweep.py").resolve())
        settings["output_png"] = str((temp_dir_path / "sweep.png").resolve())
        temp_config_path = temp_dir_path / base_config_path.name
        with temp_config_path.open("w", encoding="utf-8") as handle:
            yaml.safe_dump(config, handle, sort_keys=False)

        proc = run_simulation(exec_path, temp_config_path, project_root)
        if proc.returncode != 0:
            raise RuntimeError(proc.stderr or proc.stdout or "native sweep execution failed")
        if not sweep_csv.exists():
            raise FileNotFoundError(f"Expected sweep output was not generated: {sweep_csv}")

        sweep_data = np.loadtxt(sweep_csv, delimiter=",", skiprows=1)
        sweep_data = np.atleast_2d(sweep_data)
        return sweep_data[:, 0], sweep_data[:, 1]


def run_pointwise_sweep(omega_min: float,
                        omega_max: float,
                        steps: int,
                        base_config: dict,
                        base_config_path: Path,
                        exec_path: Path,
                        project_root: Path) -> tuple[np.ndarray, np.ndarray]:
    omega0 = math.sqrt(max(0.0, linearized_omega0_sq(base_config)))
    omegas = build_omega_grid(omega_min, omega_max, steps, omega0)
    amplitudes: list[float] = []

    print(f"Sweeping omega_drive from {omega_min} to {omega_max} in {steps} steps...")

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir_path = Path(temp_dir)
        for index, omega in enumerate(omegas):
            print(f"  Step {index + 1}/{len(omegas)}: omega_drive = {omega:.6f} ...", end=" ", flush=True)
            config = yaml.safe_load(yaml.dump(base_config))
            settings = config.setdefault("settings", {})
            settings["run_plotter"] = False
            settings["show_plot"] = False
            settings["save_png"] = False
            settings["python_script"] = str((temp_dir_path / f"plot_{index}.py").resolve())
            settings["output_png"] = str((temp_dir_path / f"plot_{index}.png").resolve())
            config.setdefault("physical", {})["omega_drive"] = float(omega)
            data_file = temp_dir_path / f"point_{index}.csv"
            settings["data_file"] = str(data_file.resolve())
            temp_config_path = temp_dir_path / f"{base_config_path.stem}_{index}.yaml"
            with temp_config_path.open("w", encoding="utf-8") as handle:
                yaml.safe_dump(config, handle, sort_keys=False)

            proc = run_simulation(exec_path, temp_config_path, project_root)
            if proc.returncode != 0:
                raise RuntimeError(proc.stderr or proc.stdout or f"simulation failed for omega={omega}")
            if not data_file.exists():
                raise FileNotFoundError(f"Expected data file was not generated: {data_file}")

            data = np.loadtxt(data_file, delimiter=",", skiprows=1)
            data = np.atleast_2d(data)
            times = data[:, 0]
            thetas = data[:, 2]
            if np.any(np.isnan(thetas)) or np.any(np.isinf(thetas)):
                raise RuntimeError(f"Encountered unstable solution at omega={omega}")

            amplitude = steady_state_amplitude(times, thetas, float(omega))
            amplitudes.append(amplitude)
            print(f"Amp = {amplitude:.6f}")

    return omegas, np.asarray(amplitudes)


def summarise_drift(omegas: np.ndarray,
                    numerical_amplitudes: np.ndarray,
                    analytical_amplitudes: np.ndarray,
                    config: dict) -> dict[str, float]:
    abs_drift = np.abs(numerical_amplitudes - analytical_amplitudes)
    rel_drift = abs_drift / np.maximum(np.abs(analytical_amplitudes), 1e-30)

    peak_index_numerical = int(np.argmax(numerical_amplitudes))
    peak_index_analytical = int(np.argmax(analytical_amplitudes))

    return {
        "omega0": math.sqrt(max(0.0, linearized_omega0_sq(config))),
        "theoretical_peak_frequency": theoretical_peak_frequency(config),
        "numerical_peak_frequency": float(omegas[peak_index_numerical]),
        "analytical_peak_frequency": float(omegas[peak_index_analytical]),
        "peak_frequency_abs_drift": abs(
            float(omegas[peak_index_numerical]) - float(omegas[peak_index_analytical])
        ),
        "numerical_peak_amplitude": float(numerical_amplitudes[peak_index_numerical]),
        "analytical_peak_amplitude": float(analytical_amplitudes[peak_index_analytical]),
        "peak_amplitude_abs_drift": abs(
            float(numerical_amplitudes[peak_index_numerical]) -
            float(analytical_amplitudes[peak_index_analytical])
        ),
        "mean_abs_drift": float(np.mean(abs_drift)),
        "max_abs_drift": float(np.max(abs_drift)),
        "rmse_abs_drift": float(np.sqrt(np.mean(abs_drift ** 2))),
        "mean_rel_drift": float(np.mean(rel_drift)),
        "max_rel_drift": float(np.max(rel_drift)),
    }


def enforce_thresholds(summary: dict[str, float],
                       max_rmse: float | None,
                       max_abs_drift: float | None,
                       max_peak_frequency_drift: float | None,
                       max_peak_amplitude_drift: float | None) -> None:
    failures: list[str] = []
    if max_rmse is not None and summary["rmse_abs_drift"] > max_rmse:
        failures.append(
            f"rmse_abs_drift={summary['rmse_abs_drift']:.6e} exceeds {max_rmse:.6e}"
        )
    if max_abs_drift is not None and summary["max_abs_drift"] > max_abs_drift:
        failures.append(
            f"max_abs_drift={summary['max_abs_drift']:.6e} exceeds {max_abs_drift:.6e}"
        )
    if (max_peak_frequency_drift is not None and
            summary["peak_frequency_abs_drift"] > max_peak_frequency_drift):
        failures.append(
            "peak_frequency_abs_drift="
            f"{summary['peak_frequency_abs_drift']:.6e} exceeds {max_peak_frequency_drift:.6e}"
        )
    if (max_peak_amplitude_drift is not None and
            summary["peak_amplitude_abs_drift"] > max_peak_amplitude_drift):
        failures.append(
            "peak_amplitude_abs_drift="
            f"{summary['peak_amplitude_abs_drift']:.6e} exceeds {max_peak_amplitude_drift:.6e}"
        )

    if failures:
        raise RuntimeError("Resonance drift checks failed:\n  " + "\n  ".join(failures))


def analyze_resonance(args: argparse.Namespace, project_root: Path) -> None:
    exec_path = project_root / "integrated_pendulum"
    if not exec_path.exists():
        raise FileNotFoundError(f"Executable not found at {exec_path}")

    with args.base_config.open("r", encoding="utf-8") as handle:
        base_config = yaml.safe_load(handle)

    sweep_config = base_config.get("sweep", {})
    use_native_sweep = bool(
        args.use_native_sweep or
        (isinstance(sweep_config, dict) and bool(sweep_config.get("enabled", False)))
    )

    if use_native_sweep:
        omegas, numerical_amplitudes = run_native_sweep(
            base_config,
            args.base_config,
            exec_path,
            project_root,
        )
    else:
        omega_min = args.omega_min
        omega_max = args.omega_max
        steps = args.steps
        if omega_min is None or omega_max is None or steps is None:
            raise ValueError("Point-by-point mode requires --omega-min, --omega-max, and --steps.")
        omegas, numerical_amplitudes = run_pointwise_sweep(
            omega_min,
            omega_max,
            steps,
            base_config,
            args.base_config,
            exec_path,
            project_root,
        )

    analytical_amplitudes = np.asarray(
        [analytical_amplitude(base_config, float(omega)) for omega in omegas]
    )
    summary = summarise_drift(omegas, numerical_amplitudes, analytical_amplitudes, base_config)
    enforce_thresholds(
        summary,
        args.max_rmse,
        args.max_abs_drift,
        args.max_peak_frequency_drift,
        args.max_peak_amplitude_drift,
    )

    abs_drift = np.abs(numerical_amplitudes - analytical_amplitudes)
    rel_drift = abs_drift / np.maximum(np.abs(analytical_amplitudes), 1e-30)

    output_dir = project_root / "outputs"
    output_dir.mkdir(exist_ok=True)
    output_prefix = args.output_prefix or "resonance_plot"
    curve_csv = output_dir / f"{output_prefix}.csv"
    drift_csv = output_dir / f"{output_prefix}_drift.csv"
    plot_path = output_dir / f"{output_prefix}.png"

    curve_rows = [
        [
            float(omega),
            float(num_amp),
            float(ana_amp),
            float(abs_err),
            float(rel_err),
        ]
        for omega, num_amp, ana_amp, abs_err, rel_err in zip(
            omegas,
            numerical_amplitudes,
            analytical_amplitudes,
            abs_drift,
            rel_drift,
        )
    ]
    write_csv(
        curve_csv,
        [
            "DriveFrequency",
            "NumericalAmplitude",
            "AnalyticalAmplitude",
            "AbsoluteDrift",
            "RelativeDrift",
        ],
        curve_rows,
    )
    write_csv(
        drift_csv,
        list(summary.keys()),
        [[summary[key] for key in summary.keys()]],
    )

    plt.figure(figsize=(10, 6))
    plt.plot(omegas, numerical_amplitudes, "o-", linewidth=2, markersize=5, label="Numerical")
    plt.plot(omegas, analytical_amplitudes, "--", linewidth=2, label="Analytical")
    plt.axvline(
        x=summary["omega0"],
        color="g",
        linestyle=":",
        linewidth=1.2,
        label=f"Linearized $\\omega_0$ = {summary['omega0']:.3f}",
    )
    plt.xlabel("Driving Frequency $\\omega_{drive}$ (rad/s)", fontsize=14)
    plt.ylabel("Steady-State Amplitude", fontsize=14)
    plt.title("Driven Pendulum Resonance Curve", fontsize=16)
    plt.grid(True, alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(plot_path, dpi=150)

    print(f"Saved resonance data to {curve_csv}")
    print(f"Saved resonance drift summary to {drift_csv}")
    print(f"Saved resonance plot to {plot_path}")
    print(
        "Drift summary: "
        f"rmse={summary['rmse_abs_drift']:.6e}, "
        f"max_abs={summary['max_abs_drift']:.6e}, "
        f"peak_freq_drift={summary['peak_frequency_abs_drift']:.6e}, "
        f"peak_amp_drift={summary['peak_amplitude_abs_drift']:.6e}"
    )

    if args.show:
        plt.show()
    else:
        plt.close()


def main() -> int:
    parser = argparse.ArgumentParser(description="Analyze resonance of a driven pendulum.")
    parser.add_argument("--omega-min", type=float, default=None, help="Minimum driving frequency (rad/s)")
    parser.add_argument("--omega-max", type=float, default=None, help="Maximum driving frequency (rad/s)")
    parser.add_argument("--steps", type=int, default=None, help="Number of frequency steps")
    parser.add_argument(
        "--base-config",
        type=Path,
        default=Path("QA/driven_pendulum/driven_pendulum.yaml"),
        help="Path to base YAML config file",
    )
    parser.add_argument(
        "--output-prefix",
        default=None,
        help="Output file stem under outputs/",
    )
    parser.add_argument(
        "--use-native-sweep",
        action="store_true",
        help="Run the YAML's native sweep mode instead of repeated single-frequency simulations",
    )
    parser.add_argument("--max-rmse", type=float, default=None, help="Fail if amplitude RMSE exceeds this value")
    parser.add_argument(
        "--max-abs-drift",
        type=float,
        default=None,
        help="Fail if the maximum amplitude drift exceeds this value",
    )
    parser.add_argument(
        "--max-peak-frequency-drift",
        type=float,
        default=None,
        help="Fail if the numerical and analytical peak frequencies diverge by more than this value",
    )
    parser.add_argument(
        "--max-peak-amplitude-drift",
        type=float,
        default=None,
        help="Fail if the numerical and analytical peak amplitudes diverge by more than this value",
    )
    parser.add_argument("--show", action="store_true", help="Display the generated plot")

    args = parser.parse_args()

    project_root = Path(__file__).resolve().parent.parent
    args.base_config = (project_root / args.base_config).resolve()
    if not args.base_config.exists():
        raise FileNotFoundError(f"Base config file not found: {args.base_config}")

    analyze_resonance(args, project_root)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
