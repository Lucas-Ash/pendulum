#!/usr/bin/env python3
"""
Run QA YAML scripts and optionally compare results to baseline outputs.

Usage:
  python run_qa_tests.py [--compare] [--exec PATH] [--qa-dir PATH]

Options:
  --compare, -c    Compare new results to existing outputs; report differences
  --exec PATH      Path to integrated_pendulum executable (default: ./integrated_pendulum)
  --qa-dir PATH    QA directory (default: QA)
  --tolerance TOL  Max absolute difference to consider OK when comparing (default: 1e-10)
"""

import argparse
import csv
from dataclasses import dataclass
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def find_qa_yamls(qa_dir: Path) -> list[Path]:
    """Find all .yaml config files in the QA directory."""
    return sorted(qa_dir.rglob("*.yaml"))


def get_output_paths(config_path: Path) -> list[Path]:
    """Extract output file paths from a YAML config. Returns paths relative to project root."""
    paths = []
    data_file_re = re.compile(
        r'^\s*(?:data_file|sweep_data_file):\s*["\']?([^"\'\s]+)["\']?\s*(?:#.*)?$'
    )
    with open(config_path) as f:
        for line in f:
            m = data_file_re.match(line.strip())
            if m:
                paths.append(Path(m.group(1)))
    return paths


def load_data_file(path: Path) -> tuple[list[str], list[list[float]]]:
    """
    Load a data file (CSV or space-separated .dat).
    Returns (column_names, rows) where rows are lists of floats.
    """
    path = Path(path)
    if not path.exists():
        return [], []

    rows = []
    columns = []

    with open(path) as f:
        if path.suffix.lower() == ".csv":
            reader = csv.reader(f)
            lines = list(reader)
            if lines:
                columns = lines[0]
                for line in lines[1:]:
                    rows.append([float(x) for x in line])
        else:
            first_data_line = True
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith("#"):
                    columns = line.lstrip("#").split()
                    continue
                if "," in line:
                    parts = [part.strip() for part in line.split(",")]
                    if first_data_line and any(not part.replace(".", "", 1).replace("-", "", 1).replace("+", "", 1).replace("e", "", 1).replace("E", "", 1).isdigit() for part in parts):
                        columns = parts
                        first_data_line = False
                        continue
                    rows.append([float(x) for x in parts])
                    first_data_line = False
                    continue
                parts = line.split()
                if parts:
                    rows.append([float(x) for x in parts])
                    first_data_line = False

    return columns, rows


def compare_data(
    baseline_path: Path,
    new_path: Path,
    name: str,
) -> dict:
    """
    Compare two data files. Returns summary statistics.
    """
    bl_cols, bl_rows = load_data_file(baseline_path)
    new_cols, new_rows = load_data_file(new_path)

    result = {
        "name": name,
        "file": str(new_path),
        "baseline_rows": len(bl_rows),
        "new_rows": len(new_rows),
        "row_match": len(bl_rows) == len(new_rows),
        "columns": bl_cols or new_cols,
        "max_abs_diff": {},
        "mean_abs_diff": {},
        "rmse": {},
        "max_rel_diff": {},
    }

    if not bl_rows or not new_rows:
        result["error"] = "Empty baseline or new data"
        return result

    ncols = min(len(bl_cols) if bl_cols else len(bl_rows[0]), len(bl_rows[0]), len(new_rows[0]))
    nrows = min(len(bl_rows), len(new_rows))
    cols = bl_cols[:ncols] if bl_cols else [f"col_{i}" for i in range(ncols)]

    for i, col in enumerate(cols):
        if i >= len(bl_rows[0]) or i >= len(new_rows[0]):
            break
        bl_vals = [r[i] for r in bl_rows[:nrows]]
        new_vals = [r[i] for r in new_rows[:nrows]]
        diffs = [abs(a - b) for a, b in zip(bl_vals, new_vals)]
        result["max_abs_diff"][col] = max(diffs) if diffs else 0
        result["mean_abs_diff"][col] = sum(diffs) / len(diffs) if diffs else 0
        result["rmse"][col] = (sum(d * d for d in diffs) / len(diffs)) ** 0.5 if diffs else 0
        # Relative diff (avoid div by zero)
        rel_diffs = []
        for a, b in zip(bl_vals, new_vals):
            denom = abs(a) if abs(a) > 1e-30 else 1e-30
            rel_diffs.append(abs(a - b) / denom)
        result["max_rel_diff"][col] = max(rel_diffs) if rel_diffs else 0

    return result


def run_simulation(exec_path: Path, config_path: Path, cwd: Path) -> subprocess.CompletedProcess:
    """Run the pendulum simulation for a config file."""
    env = os.environ.copy()
    env["QA_TEST"] = "1"
    return subprocess.run(
        [str(exec_path.resolve()), str(config_path.resolve())],
        cwd=str(cwd),
        env=env,
        capture_output=True,
        text=True,
        timeout=300,
    )


@dataclass(frozen=True)
class AnalysisArtifact:
    name: str
    baseline_path: Path
    generated_path: Path


@dataclass(frozen=True)
class AnalysisRegression:
    name: str
    artifacts: tuple[AnalysisArtifact, ...]
    command: list[str]
    timeout_seconds: int = 300


def build_analysis_regressions(project_root: Path) -> list[AnalysisRegression]:
    return [
        AnalysisRegression(
            name="driven_bifurcation_v",
            artifacts=(
                AnalysisArtifact(
                    name="driven_bifurcation_v",
                    baseline_path=project_root / "QA" / "driven_pendulum" / "outputs" / "driven_bifurcation_v_reference.csv",
                    generated_path=project_root / "outputs" / "qa_driven_bifurcation_v.csv",
                ),
            ),
            command=[
                sys.executable,
                str(project_root / "scripts" / "analyze_driven_chaos.py"),
                "bifurcation",
                "--base-config",
                str(project_root / "QA" / "driven_bifurcation.yaml"),
                "--parameter",
                "damping",
                "--min",
                "0.725",
                "--max",
                "0.746",
                "--steps",
                "192",
                "--warmup-periods",
                "160",
                "--sample-periods",
                "24",
                "--steps-per-period",
                "240",
                "--workers",
                "1",
                "--section-coordinate",
                "v",
                "--output-prefix",
                "qa_driven_bifurcation_v",
            ],
        ),
        AnalysisRegression(
            name="driven_resonance_sweep",
            artifacts=(
                AnalysisArtifact(
                    name="driven_resonance_sweep_curve",
                    baseline_path=project_root / "QA" / "driven_pendulum" / "outputs" / "driven_resonance_reference.csv",
                    generated_path=project_root / "outputs" / "qa_driven_resonance.csv",
                ),
                AnalysisArtifact(
                    name="driven_resonance_sweep_drift",
                    baseline_path=project_root / "QA" / "driven_pendulum" / "outputs" / "driven_resonance_drift_reference.csv",
                    generated_path=project_root / "outputs" / "qa_driven_resonance_drift.csv",
                ),
            ),
            command=[
                sys.executable,
                str(project_root / "scripts" / "analyze_resonance.py"),
                "--base-config",
                str(project_root / "QA" / "driven_pendulum" / "driven_pendulum_sweep.yaml"),
                "--output-prefix",
                "qa_driven_resonance",
                "--use-native-sweep",
                "--max-rmse",
                "2e-3",
                "--max-abs-drift",
                "5e-3",
                "--max-peak-frequency-drift",
                "8e-2",
                "--max-peak-amplitude-drift",
                "5e-3",
            ],
        ),
        AnalysisRegression(
            name="driven_lyapunov",
            artifacts=(
                AnalysisArtifact(
                    name="driven_lyapunov",
                    baseline_path=project_root / "QA" / "driven_pendulum" / "outputs" / "driven_lyapunov_reference.csv",
                    generated_path=project_root / "outputs" / "qa_driven_lyapunov.csv",
                ),
            ),
            command=[
                sys.executable,
                str(project_root / "scripts" / "analyze_driven_chaos.py"),
                "lyapunov",
                "--base-config",
                str(project_root / "QA" / "driven_pendulum" / "driven_pendulum.yaml"),
                "--transient-periods",
                "120",
                "--sample-periods",
                "240",
                "--steps-per-period",
                "300",
                "--output-prefix",
                "qa_driven_lyapunov",
            ],
        ),
    ]


def run_analysis_regression(case: AnalysisRegression, cwd: Path) -> subprocess.CompletedProcess:
    env = os.environ.copy()
    env.setdefault("MPLBACKEND", "Agg")
    env.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "pendulum_mplconfig"))
    return subprocess.run(
        case.command,
        cwd=str(cwd),
        env=env,
        capture_output=True,
        text=True,
        timeout=case.timeout_seconds,
    )


def print_comparison_result(result: dict, relative_path: str, tolerance: float) -> bool:
    if "error" in result:
        print(f"\n  {result['name']} / {result['file']}: {result['error']}")
        return False

    row_issue = ""
    if not result["row_match"]:
        row_issue = (
            f" (ROW COUNT: baseline={result['baseline_rows']} vs new={result['new_rows']})"
        )
    print(f"\n  {result['name']}: {relative_path}{row_issue}")

    ok = result["row_match"]
    for col in (result["columns"] or list(result["max_abs_diff"].keys())):
        if col not in result["max_abs_diff"]:
            continue
        max_d = result["max_abs_diff"][col]
        mean_d = result["mean_abs_diff"][col]
        rmse = result["rmse"][col]
        max_r = result["max_rel_diff"][col]
        status = "OK" if max_d < tolerance else "DIFF"
        if max_d >= tolerance:
            ok = False
        print(
            f"    {col}: max_abs={max_d:.2e} mean_abs={mean_d:.2e} "
            f"rmse={rmse:.2e} max_rel={max_r:.2e} [{status}]"
        )

    return ok


def main():
    parser = argparse.ArgumentParser(
        description="Run QA YAML scripts and optionally compare to baseline outputs"
    )
    parser.add_argument(
        "--compare", "-c",
        action="store_true",
        help="Compare new results to existing outputs; report differences",
    )
    parser.add_argument(
        "--exec",
        default="./integrated_pendulum",
        help="Path to integrated_pendulum executable",
    )
    parser.add_argument(
        "--qa-dir",
        default="QA",
        help="QA directory path",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=1e-10,
        help="Max absolute difference to consider OK when comparing (default: 1e-10)",
    )
    args = parser.parse_args()

    project_root = Path(__file__).resolve().parent.parent
    qa_dir = (project_root / args.qa_dir).resolve()
    exec_path = (project_root / args.exec).resolve()

    if not qa_dir.exists():
        print(f"Error: QA directory not found: {qa_dir}", file=sys.stderr)
        sys.exit(1)
    if not exec_path.exists():
        print(f"Error: Executable not found: {exec_path}", file=sys.stderr)
        print("Run ./compile_integrated.sh first.", file=sys.stderr)
        sys.exit(1)

    yamls = find_qa_yamls(qa_dir)
    if not yamls:
        print(f"No YAML configs found in {qa_dir}", file=sys.stderr)
        sys.exit(1)

    # Collect all output paths from configs
    output_paths = []
    for yp in yamls:
        for op in get_output_paths(yp):
            full_path = (project_root / op).resolve()
            if full_path not in [p[1] for p in output_paths]:
                output_paths.append((yp.stem, full_path))

    # Baseline backup when comparing
    baseline_dir = None
    if args.compare and output_paths:
        baseline_dir = Path(tempfile.mkdtemp(prefix="qa_baseline_"))
        for _, op in output_paths:
            if op.exists():
                rel = op.relative_to(project_root)
                dest = baseline_dir / rel
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(op, dest)

    # Run simulations
    print("Running QA simulations...\n")
    failed = []
    for yp in yamls:
        rel = yp.relative_to(project_root)
        print(f"  [{rel}]", end=" ")
        proc = run_simulation(exec_path, yp, project_root)
        if proc.returncode != 0:
            print("FAILED")
            failed.append((rel, proc.stderr or proc.stdout or ""))
        else:
            print("OK")

    if failed:
        print("\nFailures:")
        for rel, err in failed:
            print(f"  {rel}:")
            for line in (err or "").strip().split("\n")[-5:]:
                print(f"    {line}")
        # Exclude failed tests' outputs from comparison
        failed_stems = {Path(rel).stem for rel, _ in failed}
        output_paths = [(n, p) for n, p in output_paths if n not in failed_stems]

    # Compare if requested
    all_ok = True
    if args.compare and baseline_dir:
        print("\n" + "=" * 60)
        print("Comparison: new results vs baseline (existing outputs)")
        print("=" * 60)

        for name, new_path in output_paths:
            rel = new_path.relative_to(project_root)
            baseline_path = baseline_dir / rel
            r = compare_data(baseline_path, new_path, name)
            all_ok = print_comparison_result(r, str(rel), args.tolerance) and all_ok

        shutil.rmtree(baseline_dir, ignore_errors=True)

    analysis_failures = []
    analysis_regressions = build_analysis_regressions(project_root)
    analysis_results: list[tuple[AnalysisArtifact, dict]] = []

    if analysis_regressions:
        print("\nRunning analysis regressions...\n")
        for case in analysis_regressions:
            print(f"  [{case.name}]", end=" ")
            proc = run_analysis_regression(case, project_root)
            if proc.returncode != 0:
                print("FAILED")
                analysis_failures.append((case.name, proc.stderr or proc.stdout or ""))
                continue

            missing_outputs = [
                artifact.generated_path
                for artifact in case.artifacts
                if not artifact.generated_path.exists()
            ]
            if missing_outputs:
                print("FAILED")
                analysis_failures.append(
                    (
                        case.name,
                        "Expected analysis output was not generated: "
                        + ", ".join(str(path) for path in missing_outputs),
                    )
                )
                continue

            print("OK")
            if args.compare:
                for artifact in case.artifacts:
                    analysis_results.append(
                        (
                            artifact,
                            compare_data(
                                artifact.baseline_path,
                                artifact.generated_path,
                                artifact.name,
                            ),
                        )
                    )

    if analysis_failures:
        print("\nAnalysis regression failures:")
        for name, err in analysis_failures:
            print(f"  {name}:")
            for line in (err or "").strip().split("\n")[-8:]:
                print(f"    {line}")
        all_ok = False

    if args.compare and analysis_results:
        print("\n" + "=" * 60)
        print("Comparison: analysis regressions vs reference values")
        print("=" * 60)
        for artifact, result in analysis_results:
            rel = artifact.generated_path.relative_to(project_root)
            all_ok = print_comparison_result(result, str(rel), args.tolerance) and all_ok

    if args.compare:
        if all_ok:
            print(f"\n  All comparisons within numerical tolerance (max_abs < {args.tolerance}).")
        else:
            print("\n  Some differences detected. Review the metrics above.")

    print("\nDone.")
    if failed or analysis_failures or (args.compare and not all_ok):
        sys.exit(1)


if __name__ == "__main__":
    main()
