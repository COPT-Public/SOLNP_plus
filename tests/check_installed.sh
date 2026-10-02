#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
PY="${PYTHON:-python3}"
CC="${CC:-cc}"
mkdir -p tests/logs/installed-api
if [ -d .coinor-deps/hdsdp-install ]; then
  prefix="$PWD/.coinor-deps/hdsdp-install"
  mkl_lib="${MKLROOT:?Set MKLROOT}/lib/intel64"
  [ -d "$mkl_lib" ] || mkl_lib="$MKLROOT/lib"
  "$CC" -O2 tests/test_sdp.c -I"$prefix/include/hdsdp" -L"$prefix/lib" -lhdsdp \
    -L"$mkl_lib" -lmkl_rt -lm -lpthread -ldl -o tests/logs/installed-api/test_installed
  "$PY" tests/run.py --project HDSDP --binary "$PWD/tests/logs/installed-api/test_installed" --output tests/logs/installed-api/tiny
  "$PY" tests/run.py --project HDSDP --case trace --binary "$PWD/tests/logs/installed-api/test_installed" --output tests/logs/installed-api/trace
else
  prefix="$PWD/.coinor-deps/solnp-install"
  export LD_LIBRARY_PATH="$prefix/lib:${OSQP_HOME:-$PWD/.coinor-deps/osqp}/lib:${LD_LIBRARY_PATH:-}"
  "$CC" -O2 tests/test_nlp.c -I"$prefix/include/solnp" -L"$prefix/lib" -lsolnp -lm \
    -o tests/logs/installed-api/test_installed
  "$PY" tests/run.py --project SOLNP_plus --binary "$PWD/tests/logs/installed-api/test_installed" --output tests/logs/installed-api/tiny
fi
