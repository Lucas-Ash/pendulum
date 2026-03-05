#!/usr/bin/env python3
"""Plot integrator convergence curves from convergence_results.csv."""

import argparse
import csv
from collections import defaultdict
import math
from pathlib import Path
import statistics

import matplotlib.pyplot as plt


def load_series(csv_path: str):
    series = defaultdict(list)
    with open(csv_path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            method = row["integrator"].strip()
            dt = float(row["dt"])
            err = float(row["rms_error"])
            series[method].append((dt, err))
    return series


def fit_loglog_gradient(points):
    x_vals = [math.log(p[0]) for p in points]
    y_vals = [math.log(max(p[1], 1e-300)) for p in points]
    n = len(points)
    sum_x = sum(x_vals)
    sum_y = sum(y_vals)
    sum_xx = sum(x * x for x in x_vals)
    sum_xy = sum(x * y for x, y in zip(x_vals, y_vals))
    denom = n * sum_xx - sum_x * sum_x
    if abs(denom) < 1e-18:
        raise RuntimeError("Degenerate data for gradient fit.")
    slope = (n * sum_xy - sum_x * sum_y) / denom
    return slope


def detect_saturation(points, floor_window=3, saturation_factor=8.0, min_fit_points=3):
    """
    Split points into pre-saturation (fit region) and saturation region.

    points: sorted by dt descending (coarse->fine).
    """
    errors = [p[1] for p in points]
    tail_count = min(floor_window, len(errors))
    floor_est = statistics.median(sorted(errors)[:tail_count])
    threshold = floor_est * saturation_factor

    sat_start = None
    for i, err in enumerate(errors):
        if err <= threshold and i >= min_fit_points:
            sat_start = i
            break

    if sat_start is None:
        fit_points = points
        sat_points = []
    else:
        fit_points = points[:sat_start]
        sat_points = points[sat_start:]

    if len(fit_points) < min_fit_points:
        fit_points = points
        sat_points = []
        sat_start = None

    return {
        "fit_points": fit_points,
        "sat_points": sat_points,
        "floor_est": floor_est,
        "threshold": threshold,
        "sat_start": sat_start,
    }


def write_summary(summary_path, rows):
    with open(summary_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "integrator",
                "gradient_filtered",
                "noise_floor_estimate",
                "saturation_threshold",
                "fit_points",
                "saturated_points",
                "saturation_dt_start",
            ]
        )
        for row in rows:
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(description="Plot RK convergence across resolutions.")
    parser.add_argument(
        "--input",
        default="tests/artifacts/convergence/convergence_results.csv",
        help="Input CSV path.",
    )
    parser.add_argument(
        "--output",
        default="tests/artifacts/convergence/convergence_results.png",
        help="Output PNG path.",
    )
    parser.add_argument(
        "--summary",
        default="tests/artifacts/convergence/convergence_summary.csv",
        help="Output CSV path for filtered convergence gradients.",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display the figure interactively in addition to saving PNG.",
    )
    args = parser.parse_args()

    series = load_series(args.input)
    if not series:
        raise RuntimeError(f"No convergence data found in {args.input}")

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path = Path(args.summary)
    summary_path.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(10, 6))
    summary_rows = []
    for method in sorted(series):
        points = sorted(series[method], key=lambda item: item[0], reverse=True)

        sat = detect_saturation(points)
        fit_points = sat["fit_points"]
        sat_points = sat["sat_points"]
        slope = fit_loglog_gradient(fit_points)

        dts_all = [p[0] for p in points]
        errs_all = [p[1] for p in points]
        dts_fit = [p[0] for p in fit_points]
        errs_fit = [p[1] for p in fit_points]
        dts_sat = [p[0] for p in sat_points]
        errs_sat = [p[1] for p in sat_points]

        line = plt.loglog(
            dts_all,
            errs_all,
            marker="o",
            linewidth=1.2,
            alpha=0.35,
            linestyle="--",
            label=f"{method.upper()} raw",
        )[0]
        color = line.get_color()

        plt.loglog(
            dts_fit,
            errs_fit,
            marker="o",
            linewidth=2.4,
            color=color,
            label=f"{method.upper()} filtered (grad={slope:.3f})",
        )
        if dts_sat:
            plt.loglog(
                dts_sat,
                errs_sat,
                marker="x",
                linewidth=0,
                markersize=8,
                color=color,
                label=f"{method.upper()} saturation",
            )

        sat_dt = dts_sat[0] if dts_sat else ""
        summary_rows.append(
            [
                method,
                f"{slope:.12g}",
                f"{sat['floor_est']:.12g}",
                f"{sat['threshold']:.12g}",
                len(fit_points),
                len(sat_points),
                sat_dt,
            ]
        )
        print(
            f"{method.upper()}: gradient(filtered)={slope:.6f}, "
            f"floor~{sat['floor_est']:.3e}, fit_points={len(fit_points)}, "
            f"saturated_points={len(sat_points)}"
        )

    plt.gca().invert_xaxis()
    plt.xlabel("Time Step (dt)")
    plt.ylabel("RMS Trajectory Error")
    plt.title("Simple Pendulum Convergence (Saturation-Filtered)")
    plt.grid(True, which="both", linestyle="--", alpha=0.4)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=160)
    print(f"Saved convergence plot to {output_path}")

    write_summary(summary_path, summary_rows)
    print(f"Saved convergence summary to {summary_path}")

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
