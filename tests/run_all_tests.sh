#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

open_env() {
  if [[ -f "$ROOT_DIR/open_env/bin/activate" ]]; then
    # shellcheck disable=SC1091
    source "$ROOT_DIR/open_env/bin/activate"
  else
    echo "open_env virtual environment not found at $ROOT_DIR/open_env" >&2
    return 1
  fi
}

"$ROOT_DIR/tests/run_cpp_tests.sh" "$@"

if [[ "${RUN_QA_SCRIPT:-0}" == "1" ]]; then
  open_env
  "$ROOT_DIR/compile_integrated.sh"
  python3 "$ROOT_DIR/scripts/run_qa_tests.py" --compare --exec "./integrated_pendulum" --qa-dir "QA"
fi
