#!/usr/bin/env python3
"""
Plot physics-centered benchmarks for specialized integrator schemes.

Each benchmark is centered on a scheme adapted to a specific problem:
- Damped oscillator: DEN3 (exponential integrator for θ̈ + 2γθ̇ + ω₀²θ = g)
- Driven damped: DEN3 (same adaptation)
- Conservative: symplectic (leapfrog, ruth4) and position-only (VV, RKN, Numerov)

Reads convergence CSVs from tests/artifacts/convergence/, generates comparison
plots showing how well adapted each niche scheme is on its native problem.
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt


# integrator -> function evaluations per step
EVALS_PER_STEP = {
    "rk3": 3,
    "rk23": 3,
    "den3": 3,
    "rk4": 4,
    "rk5": 6,
    "rkf45": 6,
}


def load_convergence_csv(csv_path: Path):
    """Load integrator, dt, rms_error from convergence CSV."""
    if not csv_path.exists():
        return None
    series = defaultdict(list)
    with open(csv_path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            method = row["integrator"].strip()
            dt = float(row["dt"])
            err = float(row["rms_error"])
            series[method].append((dt, err))
    return dict(series) if series else None


def plot_den3_adaptation(
    series: dict,
    output_path: Path,
    title: str,
    t_max: float,
):
    """Plot DEN3 vs RK schemes: convergence and performance."""
    if not series or "den3" not in series:
        print(f"Skipping {output_path.name}: no DEN3 data in {title}")
        return False

    methods = ["den3", "rk3", "rk4", "rk5", "rk23", "rkf45"]
    available = [m for m in methods if m in series]
    if len(available) < 2:
        return False

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Left: convergence (error vs dt)
    ax = axes[0]
    for method in available:
        points = sorted(series[method], key=lambda p: p[0], reverse=True)
        dts = [p[0] for p in points]
        errs = [p[1] for p in points]
        label = f"{method.upper()}" + (" (adapted)" if method == "den3" else "")
        ax.loglog(dts, errs, marker="o", linewidth=2, markersize=6, label=label)
    ax.invert_xaxis()
    ax.set_xlabel("Time Step (dt)")
    ax.set_ylabel("RMS Trajectory Error")
    ax.set_title(f"{title}: Convergence")
    ax.grid(True, which="both", linestyle="--", alpha=0.4)
    ax.legend(loc="best")

    # Right: performance (error vs total function evaluations)
    ax = axes[1]
    for method in available:
        evals = EVALS_PER_STEP.get(method, 4)
        points = sorted(series[method], key=lambda p: p[0], reverse=True)
        dts = [p[0] for p in points]
        errs = [p[1] for p in points]
        nsteps = [max(1, int(t_max / dt + 0.5)) for dt in dts]
        total_evals = [n * evals for n in nsteps]
        label = f"{method.upper()}" + (" (adapted)" if method == "den3" else "")
        ax.loglog(total_evals, errs, marker="o", linewidth=2, markersize=6, label=label)
    ax.set_xlabel("Total Function Evaluations")
    ax.set_ylabel("RMS Trajectory Error")
    ax.set_title(f"{title}: Performance (DEN3-adapted)")
    ax.grid(True, which="both", linestyle="--", alpha=0.4)
    ax.legend(loc="best")

    fig.suptitle(
        f"{title}\nDEN3 (exponential integrator) adapted for damped systems",
        fontsize=11,
    )
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=160, bbox_inches="tight")
    print(f"Saved {output_path}")
    return True


def plot_conservative_adaptation(
    series: dict,
    output_path: Path,
    t_max: float = 2.0,
):
    """Plot symplectic and position-only schemes on conservative problem."""
    symplectic = ["semi_implicit_euler", "leapfrog", "ruth4"]
    position_only = ["velocity_verlet", "runge_kutta_nystrom", "numerov"]
    evals = {
        "semi_implicit_euler": 1,
        "leapfrog": 1,
        "ruth4": 3,
        "velocity_verlet": 2,
        "runge_kutta_nystrom": 4,
        "numerov": 4,
    }

    available = [m for m in symplectic + position_only if m in series]
    if not available:
        print(f"Skipping {output_path.name}: no conservative scheme data")
        return False

    fig, ax = plt.subplots(figsize=(8, 6))
    for method in available:
        e = evals.get(method, 2)
        points = sorted(series[method], key=lambda p: p[0], reverse=True)
        dts = [p[0] for p in points]
        errs = [p[1] for p in points]
        nsteps = [max(1, int(t_max / dt + 0.5)) for dt in dts]
        total_evals = [n * e for n in nsteps]
        family = "symplectic" if method in symplectic else "position-only"
        label = f"{method.replace('_', '-')} ({family})"
        ax.loglog(total_evals, errs, marker="o", linewidth=2, markersize=6, label=label)
    ax.set_xlabel("Total Function Evaluations")
    ax.set_ylabel("RMS Trajectory Error")
    ax.set_title("Conservative pendulum: symplectic & position-only schemes")
    ax.grid(True, which="both", linestyle="--", alpha=0.4)
    ax.legend(loc="best")
    fig.suptitle(
        "Simple undamped pendulum (θ̈ = a(t,θ))\nAdapted for symplectic and Nyström-style integrators",
        fontsize=11,
    )
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=160, bbox_inches="tight")
    print(f"Saved {output_path}")
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Plot physics-centered benchmarks for specialized schemes."
    )
    parser.add_argument(
        "--artifacts-dir",
        default="tests/artifacts/convergence",
        help="Directory containing convergence CSVs.",
    )
    args = parser.parse_args()

    artifacts = Path(args.artifacts_dir)
    saved = 0

    # Damped (DEN3-adapted)
    damped_csv = artifacts / "damped_convergence_results.csv"
    damped_series = load_convergence_csv(damped_csv)
    if damped_series:
        out = artifacts / "physics_benchmark_damped_den3.png"
        if plot_den3_adaptation(
            damped_series, out, "Damped oscillator", t_max=2.0
        ):
            saved += 1

    # Driven (DEN3-adapted)
    driven_csv = artifacts / "driven_convergence_results.csv"
    driven_series = load_convergence_csv(driven_csv)
    if driven_series:
        out = artifacts / "physics_benchmark_driven_den3.png"
        if plot_den3_adaptation(
            driven_series, out, "Driven damped oscillator", t_max=2.0
        ):
            saved += 1

    # Conservative (symplectic / position-only)
    simple_csv = artifacts / "convergence_results.csv"
    simple_series = load_convergence_csv(simple_csv)
    if simple_series:
        out = artifacts / "physics_benchmark_conservative.png"
        if plot_conservative_adaptation(simple_series, out):
            saved += 1

    if saved == 0:
        print(
            "No physics benchmark plots generated. "
            "Run ./tests/run_cpp_tests.sh --plot-convergence first to produce "
            "damped_convergence_results.csv, driven_convergence_results.csv, "
            "and convergence_results.csv."
        )
    else:
        print(f"Generated {saved} physics benchmark plot(s).")


if __name__ == "__main__":
    main()
