#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="$ROOT_DIR/tests/build"
BIN="$BUILD_DIR/pendulum_tests"
PLOT_CONVERGENCE=0
SHOW_CONVERGENCE=0

open_env() {
  if [[ -f "$ROOT_DIR/open_env/bin/activate" ]]; then
    # shellcheck disable=SC1091
    source "$ROOT_DIR/open_env/bin/activate"
  else
    echo "open_env virtual environment not found at $ROOT_DIR/open_env" >&2
    return 1
  fi
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --plot-convergence)
      PLOT_CONVERGENCE=1
      shift
      ;;
    --show-convergence)
      PLOT_CONVERGENCE=1
      SHOW_CONVERGENCE=1
      shift
      ;;
    *)
      echo "Unknown option: $1" >&2
      echo "Usage: $0 [--plot-convergence] [--show-convergence]" >&2
      exit 2
      ;;
  esac
done

mkdir -p "$BUILD_DIR"

echo "[1/2] Compiling C++ test suite..."
g++ -O2 -std=c++17 -I"$ROOT_DIR" \
  "$ROOT_DIR/tests/test_main.cpp" \
  "$ROOT_DIR/tests/test_helpers.cpp" \
  "$ROOT_DIR/tests/unit/test_config_utils.cpp" \
  "$ROOT_DIR/tests/unit/test_jacobi_elliptic.cpp" \
  "$ROOT_DIR/tests/unit/test_rk_integrators.cpp" \
  "$ROOT_DIR/tests/unit/test_error_analysis.cpp" \
  "$ROOT_DIR/tests/unit/test_configs.cpp" \
  "$ROOT_DIR/tests/unit/test_simulators.cpp" \
  "$ROOT_DIR/tests/unit/test_io.cpp" \
  "$ROOT_DIR/tests/integration/test_serial_integration.cpp" \
  "$ROOT_DIR/modules/config_utils.cpp" \
  "$ROOT_DIR/modules/jacobi_elliptic.cpp" \
  "$ROOT_DIR/modules/error_analysis.cpp" \
  "$ROOT_DIR/modules/experiment_config.cpp" \
  "$ROOT_DIR/modules/damped_config.cpp" \
  "$ROOT_DIR/modules/driven_config.cpp" \
  "$ROOT_DIR/modules/pendulum_simulator.cpp" \
  "$ROOT_DIR/modules/damped_simulator.cpp" \
  "$ROOT_DIR/modules/driven_simulator.cpp" \
  "$ROOT_DIR/modules/simple_io.cpp" \
  "$ROOT_DIR/modules/damped_io.cpp" \
  "$ROOT_DIR/modules/driven_io.cpp" \
  -o "$BIN"

echo "[2/2] Running C++ test suite..."
if [[ "$PLOT_CONVERGENCE" == "1" ]]; then
  open_env
  export PENDULUM_CONVERGENCE_PLOT=1
  if [[ "$SHOW_CONVERGENCE" == "1" ]]; then
    export PENDULUM_CONVERGENCE_SHOW=1
  else
    unset PENDULUM_CONVERGENCE_SHOW
  fi
else
  unset PENDULUM_CONVERGENCE_PLOT
  unset PENDULUM_CONVERGENCE_SHOW
fi

"$BIN"
