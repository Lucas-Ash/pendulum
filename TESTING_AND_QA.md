# Testing and QA Guide

This repository has two main verification layers:

1. C++ unit + integration tests
2. QA regression runs against YAML simulation configs

Use this guide as the single reference for how to run both.

## What is covered

The C++ suite (in `tests/`) covers:

- Unit tests:
  - config parsing/utilities
  - Jacobi elliptic functions
  - RK3/RK4/RK5 integrator behavior and observed order
  - error statistics
  - simulator behavior checks
  - IO format checks
- Serial integration tests:
  - simple/damped/driven config load -> simulate -> write -> read
  - convergence vs resolution across integrators

The QA regression script (`scripts/run_qa_tests.py`) runs end-to-end simulations from `QA/**/*.yaml`, and it also runs a fixed driven-chaos bifurcation sweep whose generated `v`-coordinate samples are compared against a checked-in reference CSV.

## Main entry points

Run from repo root (`/home/user/pendulum`).

### 1) C++ tests only

```bash
./tests/run_cpp_tests.sh
```

What it does:

- Compiles the C++ test binary (`tests/build/pendulum_tests`)
- Runs all tests
- Shows progress for each test with:
  - current index (`[RUN i/N]`)
  - per-test duration
  - overall percent
  - elapsed time
  - ETA

### 2) C++ tests + optional convergence plot

```bash
./tests/run_cpp_tests.sh --plot-convergence
```

Also available:

```bash
./tests/run_cpp_tests.sh --show-convergence
```

Notes:

- `--plot-convergence` generates convergence artifacts.
- `--show-convergence` also opens the plot window.
- The script activates `open_env` automatically for plotting.

### 3) All tests runner

```bash
./tests/run_all_tests.sh
```

This forwards all CLI flags to `run_cpp_tests.sh`, so these also work:

```bash
./tests/run_all_tests.sh --plot-convergence
./tests/run_all_tests.sh --show-convergence
```

## QA regression layer

To include QA baseline comparisons:

```bash
RUN_QA_SCRIPT=1 ./tests/run_all_tests.sh
```

With convergence plotting and QA in one run:

```bash
RUN_QA_SCRIPT=1 ./tests/run_all_tests.sh --plot-convergence
```

When `RUN_QA_SCRIPT=1`:

- `open_env` is activated
- `integrated_pendulum` is rebuilt (`./compile_integrated.sh`)
- `scripts/run_qa_tests.py --compare` runs over `QA/`
- the driven-chaos QA regression rebuilds `outputs/qa_driven_bifurcation_v.csv` and compares it to `QA/driven_pendulum/outputs/driven_bifurcation_v_reference.csv`
- the QA YAML set now includes `QA/van_der_pol/van_der_pol.yaml`, which checks a weakly nonlinear Van der Pol case against a first-order analytical reference
- `QA/van_der_pol/van_der_pol_strong_mu.yaml` adds a large-`mu` relaxation-oscillation case checked against a higher-resolution numerical reference

## Convergence artifacts

When convergence plotting is enabled, outputs are written to:

- `tests/artifacts/convergence/convergence_results.csv`
- `tests/artifacts/convergence/convergence_results.png`
- `tests/artifacts/convergence/convergence_summary.csv`

`convergence_summary.csv` includes, per integrator:

- filtered convergence gradient (`gradient_filtered`)
- estimated numerical floor (`noise_floor_estimate`)
- saturation threshold
- number of fit points vs saturated points
- `dt` where saturation starts

## Saturation and filtered gradients

At very fine `dt`, errors stop decreasing due to floating-point floor (often around `1e-14` for doubles in this setup).

The convergence plot script:

- detects saturation region near the floor
- excludes saturated points from the fitted gradient
- reports a filtered slope closer to the true method order

## Driven chaos analysis

For driven-pendulum post-processing beyond the C++ test/QA flow, use:

```bash
python scripts/analyze_driven_chaos.py bifurcation   --base-config QA/driven_bifurcation.yaml   --parameter A   --min 1.42 --max 1.5   --steps 700   --warmup-periods 1500   --sample-periods 180   --steps-per-period 800   --workers 0   --output-prefix driven_bifurcation_A
```

This writes a stroboscopic bifurcation CSV and PNG to `outputs/`. By default the bifurcation plot now uses the Poincare-section velocity `v` on the y-axis; use `--section-coordinate theta_wrapped` or `--section-coordinate theta_raw` to recover the older theta-based views.
Set `--workers` above `1` to split the parameter sweep across processes;
`--workers 0` means "use all detected CPU cores".

To estimate the maximal Lyapunov exponent for one driven configuration:

```bash
python scripts/analyze_driven_chaos.py lyapunov \
  --base-config QA/driven_pendulum/driven_pendulum.yaml
```

This writes the Lyapunov-history CSV/PNG to `outputs/` and prints the final
largest-exponent estimate in inverse-time units.

## Typical workflow

Fast verification:

```bash
./tests/run_cpp_tests.sh
```

Numerics-focused verification:

```bash
./tests/run_cpp_tests.sh --plot-convergence
cat tests/artifacts/convergence/convergence_summary.csv
```

Full pre-CI check:

```bash
RUN_QA_SCRIPT=1 ./tests/run_all_tests.sh --plot-convergence
```

## Exit codes

- `0`: all executed checks passed
- non-zero: at least one test/check failed

## Common notes

- You may see compile warnings from `config_utils.cpp` (`std::system` return value not used). This is a warning, not a test failure.
- QA compare mode may show differences if baseline outputs intentionally changed.
- Plotting YAML now supports `plot_phase_map`:
  - simple config: top-level `plot_phase_map` (default `true`)
  - damped/driven config: `settings.plot_phase_map` (default `false`, with top-level alias `plot_phase_map`)
- Restoring-force model options are now supported in all simulation configs:
  - `restoring_force_model`: `sine` (default) or `polynomial`
  - `restoring_force_linear`: linear term coefficient (default `1.0`)
  - `restoring_force_cubic`: cubic term coefficient (default `0.0`)
  - damped/driven support `physical.*` and top-level aliases for these keys
- Damping-model options are supported in damped/driven configs:
  - `damping_model`: `linear` (default) or `polynomial`
  - damped uses `gamma` as the legacy linear-damping parameter, with `damping_linear` as an alias for the full coefficient multiplying `theta_dot`
  - driven uses `damping` as the legacy linear-damping parameter, with `damping_linear` as an alias
  - `damping_cubic` adds the quadratic state-dependent contribution to the damping coefficient, so the damping term becomes `(damping_linear + damping_cubic * theta^2) * theta_dot`
- Damped configs also support `analytical_model`:
  - `linear` (default): the existing underdamped linear reference
  - `van_der_pol`: a first-order weakly nonlinear Van der Pol limit-cycle approximation
  - the Van der Pol analytical reference requires the standard quadratic coefficient form `-mu * (1 - theta^2)`, implemented as `damping_model: polynomial`, `damping_linear = -mu`, `damping_cubic = mu`, plus a purely linear polynomial restoring force
- Error-analysis options are supported in all simulation configs:
  - `error_analysis` / `error_mode`: `analytical` (default), `none`, or `hd_reference`
  - `error_reference_factor`: refinement factor for `hd_reference` mode (default `50`)
  - simple supports top-level keys; damped/driven support `settings.*` and top-level aliases
- For simple (undamped) runs, `analytical_model: duffing_jacobi` enables a Jacobi-elliptic analytical reference for the polynomial Duffing case.
- QA now includes `QA/duffing_oscillator/duffing_undamped.yaml`.
- QA now also includes `QA/van_der_pol/van_der_pol.yaml`.
