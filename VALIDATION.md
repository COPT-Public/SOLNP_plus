# SOLNP_plus validation - 2026-10-01

Build/install/test script exit: **0**. SOLNP+ was built on the user-authorized
Linux validation server. The source files are recorded in
`SOURCE_MANIFEST.sha256`.

Environment: Ubuntu 22.04.5 x86_64; GCC/G++ 11.4.0; CMake 3.30.5;
Python 3.10.12; NumPy 1.26.4. SOLNP+ uses BLAS/LAPACK 3.10.0 and OSQP 0.6.3;
dependencies were installed in user prefixes. This validation covers the CPU
build and installation path. Optional language interfaces and native Windows
were not tested. UTC timestamps on 2026-09-30 in raw logs are 2026-10-01 in
Shanghai.

| Case | Objective | Primal residual | Result |
|---|---:|---:|---|
| tiny | 3.94430452611e-31 | 8.882e-16 | PASS |
| installed-api/tiny | 3.94430452611e-31 | 8.882e-16 | PASS |

Additional checks: {"ctest": "1/1", "installed_api_compile": "PASS", "installed_api_numerical_cases": "1/1"}.
See `validation/SUMMARY.json` and `validation/build-and-tests.log` for full evidence.

These are installation correctness tests, not performance benchmarks or proof
for every solver feature.
