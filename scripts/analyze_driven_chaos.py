#!/usr/bin/env python3
"""
Build bifurcation diagrams and estimate maximal Lyapunov exponents
for the driven pendulum using the repository's driven-YAML format.

Examples:
  python scripts/analyze_driven_chaos.py bifurcation   
  --base-config QA/driven_bifurcation.yaml   --parameter damping   
  --min 0.725 --max 0.746   --steps 1800   --warmup-periods 300   
  --sample-periods 100   --steps-per-period 300   --workers 0   
  --output-prefix driven_bifurcation_d

  python scripts/analyze_driven_chaos.py lyapunov \
      --base-config QA/driven_pendulum/driven_pendulum.yaml
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import math
import os
from dataclasses import dataclass, replace
from pathlib import Path

import matplotlib.pyplot as plt
import yaml


@dataclass(frozen=True)
class DrivenAnalysisConfig:
    g: float
    length: float
    damping: float
    damping_model: str
    damping_cubic: float
    drive_amplitude: float
    omega_drive: float
    theta0: float
    omega0: float
    t_start: float
    base_dt: float
    restoring_model: str
    restoring_linear: float
    restoring_cubic: float


@dataclass(frozen=True)
class BifurcationSample:
    task_index: int
    parameter_value: float
    rows: list[tuple[float, int, float, float, float, float, float]]


PARAMETER_ALIASES = {
    "g": "g",
    "gravity": "g",
    "L": "length",
    "length": "length",
    "damping": "damping",
    "damping_cubic": "damping_cubic",
    "A": "drive_amplitude",
    "amplitude": "drive_amplitude",
    "drive_amplitude": "drive_amplitude",
    "omega_drive": "omega_drive",
    "theta0": "theta0",
    "omega0": "omega0",
    "restoring_force_linear": "restoring_linear",
    "restoring_linear": "restoring_linear",
    "restoring_force_cubic": "restoring_cubic",
    "restoring_cubic": "restoring_cubic",
}


def load_yaml(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        loaded = yaml.safe_load(handle)
    if not isinstance(loaded, dict):
        raise ValueError(f"Expected a mapping in {path}")
    return loaded


def nested_or_top_level(config: dict, section: str, key: str, default):
    section_value = config.get(section, {})
    if isinstance(section_value, dict) and key in section_value:
        return section_value[key]
    return config.get(key, default)


def load_analysis_config(path: Path) -> DrivenAnalysisConfig:
    raw = load_yaml(path)

    restoring_model = str(
        nested_or_top_level(raw, "physical", "restoring_force_model", "sine")
    ).strip().lower()
    if restoring_model == "duffing":
        restoring_model = "polynomial"
    if restoring_model not in {"sine", "polynomial"}:
        raise ValueError(
            "restoring_force_model must be 'sine' or 'polynomial' for driven chaos analysis"
        )

    config = DrivenAnalysisConfig(
        g=float(nested_or_top_level(raw, "physical", "g", 9.81)),
        length=float(nested_or_top_level(raw, "physical", "L", 1.0)),
        damping=float(nested_or_top_level(raw, "physical", "damping", 0.5)),
        damping_model=str(
            nested_or_top_level(raw, "physical", "damping_model", "linear")
        ).strip().lower(),
        damping_cubic=float(
            nested_or_top_level(raw, "physical", "damping_cubic", 0.0)
        ),
        drive_amplitude=float(nested_or_top_level(raw, "physical", "A", 0.5)),
        omega_drive=float(nested_or_top_level(raw, "physical", "omega_drive", 1.2)),
        theta0=float(nested_or_top_level(raw, "physical", "theta0", 0.1)),
        omega0=float(nested_or_top_level(raw, "physical", "omega0", 0.0)),
        t_start=float(nested_or_top_level(raw, "simulation", "t_start", 0.0)),
        base_dt=float(nested_or_top_level(raw, "simulation", "dt", 0.01)),
        restoring_model=restoring_model,
        restoring_linear=float(
            nested_or_top_level(raw, "physical", "restoring_force_linear", 1.0)
        ),
        restoring_cubic=float(
            nested_or_top_level(raw, "physical", "restoring_force_cubic", 0.0)
        ),
    )
    validate_config(config)
    return config


def validate_config(config: DrivenAnalysisConfig) -> None:
    if config.g <= 0.0:
        raise ValueError("g must be > 0")
    if config.length <= 0.0:
        raise ValueError("L must be > 0")
    if config.damping_model not in {"linear", "polynomial", "van_der_pol"}:
        raise ValueError("damping_model must be 'linear' or 'polynomial'")
    if config.damping_model == "linear" and config.damping < 0.0:
        raise ValueError("damping must be >= 0")
    if config.omega_drive <= 0.0:
        raise ValueError("omega_drive must be > 0 for stroboscopic analysis")


def damping_term(theta: float, omega: float, config: DrivenAnalysisConfig) -> float:
    coefficient = config.damping
    if config.damping_model in {"polynomial", "van_der_pol"}:
        coefficient += config.damping_cubic * theta * theta
    return coefficient * omega


def resolve_parameter_name(name: str) -> str:
    if name not in PARAMETER_ALIASES:
        valid = ", ".join(sorted(PARAMETER_ALIASES))
        raise ValueError(f"Unsupported parameter '{name}'. Choose from: {valid}")
    return PARAMETER_ALIASES[name]


def update_parameter(config: DrivenAnalysisConfig, parameter: str, value: float) -> DrivenAnalysisConfig:
    field = resolve_parameter_name(parameter)
    updated = replace(config, **{field: float(value)})
    validate_config(updated)
    return updated


def drive_period(config: DrivenAnalysisConfig) -> float:
    return 2.0 * math.pi / config.omega_drive


def derive_steps_per_period(config: DrivenAnalysisConfig) -> int:
    period = drive_period(config)
    if config.base_dt > 0.0:
        inferred = math.ceil(period / config.base_dt)
        return max(250, min(1200, inferred))
    return 600


def restoring_term(theta: float, config: DrivenAnalysisConfig) -> float:
    if config.restoring_model == "polynomial":
        return config.restoring_linear * theta + config.restoring_cubic * theta * theta * theta
    return math.sin(theta)


def restoring_derivative(theta: float, config: DrivenAnalysisConfig) -> float:
    if config.restoring_model == "polynomial":
        return config.restoring_linear + 3.0 * config.restoring_cubic * theta * theta
    return math.cos(theta)


def rhs(theta: float, omega: float, time_value: float, config: DrivenAnalysisConfig) -> tuple[float, float]:
    omega0_sq = config.g / config.length
    return (
        omega,
        -damping_term(theta, omega, config)
        - omega0_sq * restoring_term(theta, config)
        + config.drive_amplitude * math.cos(config.omega_drive * time_value),
    )


def tangent_rhs(
    theta: float,
    omega: float,
    delta_theta: float,
    delta_omega: float,
    time_value: float,
    config: DrivenAnalysisConfig,
) -> tuple[float, float, float, float]:
    theta_dot, omega_dot = rhs(theta, omega, time_value, config)
    omega0_sq = config.g / config.length
    damping_coefficient = config.damping
    damping_theta_derivative = 0.0
    if config.damping_model in {"polynomial", "van_der_pol"}:
        damping_coefficient += config.damping_cubic * theta * theta
        damping_theta_derivative = 2.0 * config.damping_cubic * theta
    return (
        theta_dot,
        omega_dot,
        delta_omega,
        -damping_theta_derivative * omega * delta_theta
        -damping_coefficient * delta_omega
        - omega0_sq * restoring_derivative(theta, config) * delta_theta,
    )


def rk4_step(theta: float, omega: float, time_value: float, dt: float, config: DrivenAnalysisConfig) -> tuple[float, float]:
    k1_theta, k1_omega = rhs(theta, omega, time_value, config)
    k2_theta, k2_omega = rhs(
        theta + 0.5 * dt * k1_theta,
        omega + 0.5 * dt * k1_omega,
        time_value + 0.5 * dt,
        config,
    )
    k3_theta, k3_omega = rhs(
        theta + 0.5 * dt * k2_theta,
        omega + 0.5 * dt * k2_omega,
        time_value + 0.5 * dt,
        config,
    )
    k4_theta, k4_omega = rhs(
        theta + dt * k3_theta,
        omega + dt * k3_omega,
        time_value + dt,
        config,
    )

    theta_next = theta + (dt / 6.0) * (k1_theta + 2.0 * k2_theta + 2.0 * k3_theta + k4_theta)
    omega_next = omega + (dt / 6.0) * (k1_omega + 2.0 * k2_omega + 2.0 * k3_omega + k4_omega)
    return theta_next, omega_next


def rk4_tangent_step(
    theta: float,
    omega: float,
    delta_theta: float,
    delta_omega: float,
    time_value: float,
    dt: float,
    config: DrivenAnalysisConfig,
) -> tuple[float, float, float, float]:
    k1 = tangent_rhs(theta, omega, delta_theta, delta_omega, time_value, config)
    k2 = tangent_rhs(
        theta + 0.5 * dt * k1[0],
        omega + 0.5 * dt * k1[1],
        delta_theta + 0.5 * dt * k1[2],
        delta_omega + 0.5 * dt * k1[3],
        time_value + 0.5 * dt,
        config,
    )
    k3 = tangent_rhs(
        theta + 0.5 * dt * k2[0],
        omega + 0.5 * dt * k2[1],
        delta_theta + 0.5 * dt * k2[2],
        delta_omega + 0.5 * dt * k2[3],
        time_value + 0.5 * dt,
        config,
    )
    k4 = tangent_rhs(
        theta + dt * k3[0],
        omega + dt * k3[1],
        delta_theta + dt * k3[2],
        delta_omega + dt * k3[3],
        time_value + dt,
        config,
    )

    return (
        theta + (dt / 6.0) * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]),
        omega + (dt / 6.0) * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]),
        delta_theta + (dt / 6.0) * (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2]),
        delta_omega + (dt / 6.0) * (k1[3] + 2.0 * k2[3] + 2.0 * k3[3] + k4[3]),
    )


def integrate_period(
    theta: float,
    omega: float,
    time_value: float,
    config: DrivenAnalysisConfig,
    steps_per_period: int,
) -> tuple[float, float, float]:
    if steps_per_period <= 0:
        raise ValueError("steps_per_period must be > 0")
    dt = drive_period(config) / steps_per_period
    for _ in range(steps_per_period):
        theta, omega = rk4_step(theta, omega, time_value, dt, config)
        time_value += dt
    return theta, omega, time_value


def integrate_period_with_tangent(
    theta: float,
    omega: float,
    delta_theta: float,
    delta_omega: float,
    time_value: float,
    config: DrivenAnalysisConfig,
    steps_per_period: int,
) -> tuple[float, float, float, float, float]:
    if steps_per_period <= 0:
        raise ValueError("steps_per_period must be > 0")
    dt = drive_period(config) / steps_per_period
    for _ in range(steps_per_period):
        theta, omega, delta_theta, delta_omega = rk4_tangent_step(
            theta, omega, delta_theta, delta_omega, time_value, dt, config
        )
        time_value += dt
    return theta, omega, delta_theta, delta_omega, time_value


def wrap_angle(theta: float) -> float:
    return (theta + math.pi) % (2.0 * math.pi) - math.pi


def write_csv(path: Path, header: list[str], rows: list[tuple]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def build_output_paths(project_root: Path, prefix: str) -> tuple[Path, Path]:
    out_dir = project_root / "outputs"
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir / f"{prefix}.csv", out_dir / f"{prefix}.png"


def choose_scatter_style(num_points: int, point_size: float | None, alpha: float | None) -> tuple[float, float]:
    if point_size is None:
        if num_points <= 2_000:
            point_size = 10.0
        elif num_points <= 20_000:
            point_size = 3.0
        else:
            point_size = 1.0
    if alpha is None:
        if num_points <= 2_000:
            alpha = 0.85
        elif num_points <= 20_000:
            alpha = 0.55
        else:
            alpha = 0.35
    return point_size, alpha


def padded_limits(values: list[float], explicit_min: float | None, explicit_max: float | None) -> tuple[float, float]:
    if explicit_min is not None and explicit_max is not None:
        if explicit_max <= explicit_min:
            raise ValueError("explicit axis max must be greater than explicit axis min")
        return explicit_min, explicit_max

    vmin = min(values)
    vmax = max(values)
    span = vmax - vmin
    if span <= 1e-12:
        pad = 1e-3
    else:
        pad = max(0.05 * span, 1e-4)

    lower = explicit_min if explicit_min is not None else vmin - pad
    upper = explicit_max if explicit_max is not None else vmax + pad
    if upper <= lower:
        raise ValueError("computed axis limits are invalid")
    return lower, upper


def resolve_worker_count(requested: int) -> int:
    if requested < 0:
        raise ValueError("workers must be >= 0")
    if requested == 0:
        detected = os.cpu_count() or 1
        return max(1, detected)
    return max(1, requested)


def resolve_section_coordinate(args: argparse.Namespace) -> str:
    if args.section_coordinate is not None:
        return args.section_coordinate
    if args.theta_mode == "raw":
        return "theta_raw"
    if args.theta_mode == "wrapped":
        return "theta_wrapped"
    return "v"


def bifurcation_coordinate_metadata(coordinate: str) -> tuple[int, str, str]:
    if coordinate == "theta_raw":
        return 3, "theta", "theta span"
    if coordinate == "theta_wrapped":
        return 4, "theta (wrapped to [-pi, pi])", "wrapped-theta span"
    if coordinate in {"omega", "v"}:
        return 6, "v (stroboscopic angular velocity)", "v span"
    raise ValueError(f"Unsupported section coordinate: {coordinate}")


def sample_bifurcation_parameter(
    task_index: int,
    base_config: DrivenAnalysisConfig,
    parameter: str,
    parameter_value: float,
    warmup_periods: int,
    sample_periods: int,
    steps_per_period: int,
) -> BifurcationSample:
    current = update_parameter(base_config, parameter, parameter_value)
    theta = current.theta0
    omega = current.omega0
    time_value = current.t_start

    for _ in range(warmup_periods):
        theta, omega, time_value = integrate_period(
            theta, omega, time_value, current, steps_per_period
        )

    rows: list[tuple[float, int, float, float, float, float, float]] = []
    for sample_index in range(sample_periods):
        theta, omega, time_value = integrate_period(
            theta, omega, time_value, current, steps_per_period
        )
        rows.append(
            (
                parameter_value,
                sample_index,
                time_value,
                theta,
                wrap_angle(theta),
                omega,
                omega,
            )
        )

    return BifurcationSample(
        task_index=task_index,
        parameter_value=parameter_value,
        rows=rows,
    )


def run_bifurcation(args: argparse.Namespace, project_root: Path) -> None:
    config = load_analysis_config(args.base_config)
    steps_per_period = args.steps_per_period or derive_steps_per_period(config)
    parameter = resolve_parameter_name(args.parameter)
    workers = resolve_worker_count(args.workers)
    if args.steps < 2:
        raise ValueError("steps must be >= 2 for bifurcation mode")
    if args.warmup_periods < 0:
        raise ValueError("warmup_periods must be >= 0")
    if args.sample_periods <= 0:
        raise ValueError("sample_periods must be > 0")
    parameter_values = [
        args.min + (args.max - args.min) * index / (args.steps - 1)
        for index in range(args.steps)
    ]

    rows: list[tuple] = []
    total = len(parameter_values)
    print(f"Using {workers} worker process{'es' if workers != 1 else ''} for the sweep")

    if workers == 1:
        for index, parameter_value in enumerate(parameter_values, start=1):
            sample = sample_bifurcation_parameter(
                task_index=index - 1,
                base_config=config,
                parameter=parameter,
                parameter_value=parameter_value,
                warmup_periods=args.warmup_periods,
                sample_periods=args.sample_periods,
                steps_per_period=steps_per_period,
            )
            rows.extend(sample.rows)
            print(f"[{index}/{total}] {parameter}={parameter_value:.6f}")
    else:
        completed: list[BifurcationSample] = []
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = [
                executor.submit(
                    sample_bifurcation_parameter,
                    index,
                    config,
                    parameter,
                    parameter_value,
                    args.warmup_periods,
                    args.sample_periods,
                    steps_per_period,
                )
                for index, parameter_value in enumerate(parameter_values)
            ]
            for done_count, future in enumerate(as_completed(futures), start=1):
                sample = future.result()
                completed.append(sample)
                print(f"[{done_count}/{total}] {parameter}={sample.parameter_value:.6f}")

        completed.sort(key=lambda item: item.task_index)
        for sample in completed:
            rows.extend(sample.rows)

    prefix = args.output_prefix or f"driven_bifurcation_{parameter}"
    csv_path, png_path = build_output_paths(project_root, prefix)
    write_csv(
        csv_path,
        ["parameter", "sample_index", "time", "theta", "theta_wrapped", "omega", "v"],
        rows,
    )

    coordinate = resolve_section_coordinate(args)
    y_column, y_label, span_label = bifurcation_coordinate_metadata(coordinate)
    x_values = [row[0] for row in rows]
    y_values = [row[y_column] for row in rows]
    x_min, x_max = padded_limits(x_values, args.x_min, args.x_max)
    y_min, y_max = padded_limits(y_values, args.y_min, args.y_max)
    point_size, alpha = choose_scatter_style(len(rows), args.point_size, args.alpha)

    y_span = max(y_values) - min(y_values)
    print(
        f"Plotting {len(rows)} points with {span_label} {y_span:.6e} "
        f"and y-limits [{y_min:.6f}, {y_max:.6f}]"
    )
    if y_span < 1e-2:
        print(
            "Warning: the sampled orbit occupies a very narrow band. "
            "This usually means the system is on a low-period attractor, not a broad chaotic branch."
        )

    plt.figure(figsize=(10, 6))
    plt.scatter(x_values, y_values, s=point_size, c="#111111", alpha=alpha, linewidths=0, rasterized=True)
    plt.xlabel(parameter)
    plt.ylabel(y_label)
    plt.title("Driven Pendulum Bifurcation Diagram")
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    plt.grid(True, alpha=0.2)
    plt.tight_layout()
    plt.savefig(png_path, dpi=300)

    print(f"Saved bifurcation samples to {csv_path}")
    print(f"Saved bifurcation plot to {png_path}")

    if args.show:
        plt.show()
    else:
        plt.close()


def run_lyapunov(args: argparse.Namespace, project_root: Path) -> None:
    config = load_analysis_config(args.base_config)
    steps_per_period = args.steps_per_period or derive_steps_per_period(config)
    if args.transient_periods < 0:
        raise ValueError("transient_periods must be >= 0")
    if args.sample_periods <= 0:
        raise ValueError("sample_periods must be > 0")

    theta = config.theta0
    omega = config.omega0
    time_value = config.t_start

    for _ in range(args.transient_periods):
        theta, omega, time_value = integrate_period(
            theta, omega, time_value, config, steps_per_period
        )

    delta_theta = 1.0
    delta_omega = 0.0
    period = drive_period(config)
    accumulated_log = 0.0
    history_rows: list[tuple] = []

    for sample_index in range(1, args.sample_periods + 1):
        theta, omega, delta_theta, delta_omega, time_value = integrate_period_with_tangent(
            theta,
            omega,
            delta_theta,
            delta_omega,
            time_value,
            config,
            steps_per_period,
        )

        norm = math.hypot(delta_theta, delta_omega)
        if not math.isfinite(norm) or norm <= 0.0:
            raise RuntimeError("Lyapunov tangent vector collapsed or became non-finite")

        accumulated_log += math.log(norm)
        elapsed_time = sample_index * period
        lambda_estimate = accumulated_log / elapsed_time
        history_rows.append(
            (
                sample_index,
                time_value,
                elapsed_time,
                lambda_estimate,
                norm,
                theta,
                omega,
            )
        )

        delta_theta /= norm
        delta_omega /= norm

        if sample_index % max(1, args.sample_periods // 10) == 0 or sample_index == args.sample_periods:
            print(f"[{sample_index}/{args.sample_periods}] lambda_max ~= {lambda_estimate:.8f}")

    lambda_max = history_rows[-1][3]
    prefix = args.output_prefix or "driven_lyapunov"
    csv_path, png_path = build_output_paths(project_root, prefix)
    write_csv(
        csv_path,
        [
            "period_index",
            "time",
            "elapsed_time",
            "lambda_estimate",
            "stretch_norm",
            "theta",
            "omega",
        ],
        history_rows,
    )

    times = [row[2] for row in history_rows]
    estimates = [row[3] for row in history_rows]

    plt.figure(figsize=(10, 5))
    plt.plot(times, estimates, color="#0b5394", linewidth=1.8)
    plt.axhline(0.0, color="#666666", linestyle="--", linewidth=1.0)
    plt.xlabel("Elapsed analysis time")
    plt.ylabel("Largest Lyapunov exponent")
    plt.title("Driven Pendulum Lyapunov Estimate")
    plt.grid(True, alpha=0.2)
    plt.tight_layout()
    plt.savefig(png_path, dpi=180)

    print(f"Maximal Lyapunov exponent: {lambda_max:.10f}")
    print(f"Saved Lyapunov history to {csv_path}")
    print(f"Saved Lyapunov plot to {png_path}")

    if args.show:
        plt.show()
    else:
        plt.close()


def build_parser(project_root: Path) -> argparse.ArgumentParser:
    default_config = project_root / "QA" / "driven_pendulum" / "driven_pendulum.yaml"

    parser = argparse.ArgumentParser(
        description="Driven pendulum chaos analysis tools",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    bifurcation = subparsers.add_parser(
        "bifurcation",
        help="Sweep one parameter and build a stroboscopic bifurcation diagram",
    )
    bifurcation.add_argument(
        "--base-config",
        type=Path,
        default=default_config,
        help="Driven YAML config to use as the baseline",
    )
    bifurcation.add_argument("--parameter", default="A", help="Parameter to sweep")
    bifurcation.add_argument("--min", type=float, required=True, help="Sweep lower bound")
    bifurcation.add_argument("--max", type=float, required=True, help="Sweep upper bound")
    bifurcation.add_argument("--steps", type=int, default=250, help="Number of sweep points")
    bifurcation.add_argument(
        "--warmup-periods",
        type=int,
        default=300,
        help="Driving periods to discard before sampling",
    )
    bifurcation.add_argument(
        "--sample-periods",
        type=int,
        default=80,
        help="Driving periods to record per parameter value",
    )
    bifurcation.add_argument(
        "--steps-per-period",
        type=int,
        default=None,
        help="RK4 steps per driving period; derived from config when omitted",
    )
    bifurcation.add_argument(
        "--section-coordinate",
        choices=("v", "omega", "theta_wrapped", "theta_raw"),
        default=None,
        help="Poincare-section coordinate to plot on the y-axis (default: v)",
    )
    bifurcation.add_argument(
        "--theta-mode",
        choices=("wrapped", "raw"),
        default=None,
        help="Deprecated compatibility option: plot wrapped/raw theta when --section-coordinate is omitted",
    )
    bifurcation.add_argument(
        "--output-prefix",
        default=None,
        help="Output file stem under outputs/",
    )
    bifurcation.add_argument("--point-size", type=float, default=None, help="Scatter marker size")
    bifurcation.add_argument("--alpha", type=float, default=None, help="Scatter alpha")
    bifurcation.add_argument("--x-min", type=float, default=None, help="Manual x-axis minimum")
    bifurcation.add_argument("--x-max", type=float, default=None, help="Manual x-axis maximum")
    bifurcation.add_argument("--y-min", type=float, default=None, help="Manual y-axis minimum")
    bifurcation.add_argument("--y-max", type=float, default=None, help="Manual y-axis maximum")
    bifurcation.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Worker processes for the sweep; use 0 to use all available CPU cores",
    )
    bifurcation.add_argument("--show", action="store_true", help="Display the generated plot")

    lyapunov = subparsers.add_parser(
        "lyapunov",
        help="Estimate the maximal Lyapunov exponent for one driven configuration",
    )
    lyapunov.add_argument(
        "--base-config",
        type=Path,
        default=default_config,
        help="Driven YAML config to analyze",
    )
    lyapunov.add_argument(
        "--transient-periods",
        type=int,
        default=400,
        help="Driving periods to discard before Lyapunov accumulation",
    )
    lyapunov.add_argument(
        "--sample-periods",
        type=int,
        default=1000,
        help="Driving periods used in the Lyapunov average",
    )
    lyapunov.add_argument(
        "--steps-per-period",
        type=int,
        default=None,
        help="RK4 steps per driving period; derived from config when omitted",
    )
    lyapunov.add_argument(
        "--output-prefix",
        default=None,
        help="Output file stem under outputs/",
    )
    lyapunov.add_argument("--show", action="store_true", help="Display the generated plot")

    return parser


def main() -> int:
    project_root = Path(__file__).resolve().parent.parent
    parser = build_parser(project_root)
    args = parser.parse_args()

    if not args.base_config.exists():
        raise FileNotFoundError(f"Base config file not found: {args.base_config}")

    if args.command == "bifurcation":
        run_bifurcation(args, project_root)
        return 0
    if args.command == "lyapunov":
        run_lyapunov(args, project_root)
        return 0
    raise ValueError(f"Unknown command: {args.command}")


if __name__ == "__main__":
    raise SystemExit(main())
