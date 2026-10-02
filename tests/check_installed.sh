#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
PY="${PYTHON:-python3}"
CC="${CC:-cc}"
prefix="$PWD/.coinor-deps/solnp-install"
mkdir -p tests/logs/installed-api
export LD_LIBRARY_PATH="$prefix/lib:${OSQP_HOME:-$PWD/.coinor-deps/osqp}/lib:${LD_LIBRARY_PATH:-}"
export DYLD_LIBRARY_PATH="$prefix/lib:${OSQP_HOME:-$PWD/.coinor-deps/osqp}/lib:${DYLD_LIBRARY_PATH:-}"
"$CC" -O2 tests/test_nlp.c -I"$prefix/include/solnp" -L"$prefix/lib" -lsolnp -lm \
  -o tests/logs/installed-api/test_installed
"$PY" tests/run.py --binary "$PWD/tests/logs/installed-api/test_installed" \
  --output tests/logs/installed-api/tiny
