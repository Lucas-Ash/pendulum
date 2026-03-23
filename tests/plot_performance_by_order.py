#!/usr/bin/env python3
"""
Plot performance (error vs computational cost) for integrators grouped by order.

Compares schemes of the same convergence order: at a given accuracy target,
which scheme achieves it with least cost (function evaluations)? This supplements
the convergence test by testing efficiency within each order class.

Reads: tests/artifacts/convergence/convergence_results.csv (from run_cpp_tests.sh --plot-convergence)
Output: tests/artifacts/convergence/performance_by_order.png
        tests/artifacts/convergence/performance_by_order.csv
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt


# Integrator -> (convergence order, function evaluations per step).
# evals_per_step: derivative/acceleration evaluations per integration step.
ORDER_AND_COST = {
    "semi_implicit_euler": (1, 1),
    "leapfrog": (2, 1),
    "velocity_verlet": (2, 2),
    "rk3": (3, 3),
    "rk23": (3, 3),
    "den3": (3, 3),
    "rk4": (4, 4),
    "runge_kutta_nystrom": (4, 4),
    "numerov": (4, 4),   # nominal (1 base + ~3 Newton iters typical)
    "ruth4": (4, 3),
    "rk5": (5, 6),
    "rkf45": (5, 6),
}


def load_series(csv_path: str):
    """Load integrator, dt, rms_error from convergence CSV."""
    series = defaultdict(list)
    with open(csv_path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            method = row["integrator"].strip()
            dt = float(row["dt"])
            err = float(row["rms_error"])
            series[method].append((dt, err))
    return series


def main():
    parser = argparse.ArgumentParser(
        description="Plot performance (error vs cost) for integrators by order."
    )
    parser.add_argument(
        "--input",
        default="tests/artifacts/convergence/convergence_results.csv",
        help="Input convergence CSV path.",
    )
    parser.add_argument(
        "--output",
        default="tests/artifacts/convergence/performance_by_order.png",
        help="Output PNG path.",
    )
    parser.add_argument(
        "--csv",
        default="tests/artifacts/convergence/performance_by_order.csv",
        help="Output CSV path for performance data.",
    )
    parser.add_argument(
        "--t-max",
        type=float,
        default=2.0,
        help="Simulation t_max used to generate convergence data (for nsteps = t_max/dt).",
    )
    args = parser.parse_args()

    series = load_series(args.input)
    if not series:
        raise RuntimeError(f"No convergence data found in {args.input}")

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    csv_path = Path(args.csv)
    csv_path.parent.mkdir(parents=True, exist_ok=True)

    # Group by order; only include groups with 2+ schemes
    by_order = defaultdict(list)
    for method, points in series.items():
        if method not in ORDER_AND_COST:
            continue
        order, evals = ORDER_AND_COST[method]
        by_order[order].append((method, evals, sorted(points, key=lambda p: p[0], reverse=True)))

    # Sort orders and filter to groups with multiple schemes
    orders_with_compare = [
        (order, methods) for order, methods in sorted(by_order.items())
        if len(methods) >= 2
    ]

    if not orders_with_compare:
        raise RuntimeError(
            "No order groups with 2+ schemes found. "
            "Ensure convergence_results.csv includes multiple schemes per order."
        )

    n_panels = len(orders_with_compare)
    ncols = 2 if n_panels >= 2 else 1
    nrows = (n_panels + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(7 * ncols, 5 * nrows))
    if n_panels == 1:
        axes = [axes]
    else:
        axes = axes.flatten()

    # Hide unused subplots
    for j in range(n_panels, len(axes)):
        axes[j].set_visible(False)

    csv_rows = []
    csv_rows.append(["order", "integrator", "evals_per_step", "dt", "nsteps", "total_evals", "rms_error"])

    for idx, (order, methods) in enumerate(orders_with_compare):
        ax = axes[idx]
        for method, evals_per_step, points in methods:
            dts = [p[0] for p in points]
            errs = [p[1] for p in points]
            nsteps = [max(1, int(args.t_max / dt + 0.5)) for dt in dts]
            total_evals = [n * evals_per_step for n in nsteps]

            ax.loglog(
                total_evals,
                errs,
                marker="o",
                linewidth=2,
                markersize=6,
                label=method.replace("_", "-"),
            )

            for i, (dt, err, ne, te) in enumerate(zip(dts, errs, nsteps, total_evals)):
                csv_rows.append([order, method, evals_per_step, dt, ne, te, err])

        ax.set_xlabel("Total Function Evaluations")
        ax.set_ylabel("RMS Trajectory Error")
        ax.set_title(f"Order {order} schemes: Error vs Cost")
        ax.grid(True, which="both", linestyle="--", alpha=0.4)
        ax.legend(loc="upper right")

    fig.suptitle(
        "Performance comparison: schemes of the same order\n"
        "(Lower-left curve = more efficient at given accuracy)",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(output_path, dpi=160, bbox_inches="tight")
    print(f"Saved performance plot to {output_path}")

    with open(csv_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerows(csv_rows)
    print(f"Saved performance data to {csv_path}")


if __name__ == "__main__":
    main()
