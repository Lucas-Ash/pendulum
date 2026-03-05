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

The QA regression script (`scripts/run_qa_tests.py`) runs end-to-end simulations from `QA/**/*.yaml` and can compare outputs against existing baselines.

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
