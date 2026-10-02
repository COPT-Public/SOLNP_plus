#!/usr/bin/env bash
set -euxo pipefail
cd "$(dirname "$0")"
CM="${CMAKE:-cmake}"
PY="${PYTHON:-python3}"
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MKL_THREADING_LAYER=SEQUENTIAL
"$PY" tests/generate_problem.py
if [ -z "${OSQP_HOME:-}" ]; then
  "$CM" -S source/ThirdParty/osqp -B .coinor-deps/osqp-build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PWD/.coinor-deps/osqp" -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DUNITTESTS=OFF
  "$CM" --build .coinor-deps/osqp-build -j4
  "$CM" --install .coinor-deps/osqp-build
  export OSQP_HOME="$PWD/.coinor-deps/osqp"
fi
"$CM" -S tests -B build-coinor -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PWD/.coinor-deps/solnp-install"
"$CM" --build build-coinor -j4
"$CM" --install build-coinor
"$CM" --build build-coinor --target test
"$PY" tests/run.py --project SOLNP_plus

bash tests/check_installed.sh
