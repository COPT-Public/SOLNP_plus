# Installation and checker tests

Run `bash build.sh` after following INSTALL. All tests use bundled small models
with analytic solutions and return process exit 0 only if accepted.

problem.txt: min (x-1)^2+(y-2)^2, x+y=3, -5<=x,y<=5.
Expected x=(1,2), objective 0, status=1. Independent equality residual <=1e-6,
objective <=1e-6, coordinate error <=1e-3; bounds and solver objective are checked.

`run.py` uses NumPy to independently recompute the specified feasibility and
objective conditions from returned solution data. Every invocation creates a
new run-* output folder. The parent verification.json is reset before execution
and replaced with the result; no earlier solution is read. Missing output,
timeout, nonzero process exit, nonfinite values or failed checks cause failure.
The run folder contains command.json, solver.log and raw solver outputs.

`tests/check_installed.sh` compiles the same C consumer using only headers
and libraries in the installation prefix, then reruns numerical verification.
This detects missing installed headers or unusable installed libraries.

To add a case:
1. Add a small model with an analytic expected answer to examples/ and describe
   its mathematics/expected answer in examples/README.md.
2. Extend generate_problem.py if the model is generated. In run.py add an explicit
   --case choice and the corresponding model, expected optimum and residual checks.
3. Add its command to build.sh (and add_test in the relevant CMakeLists if using
   CTest). Use the installed executable/consumer where applicable.
4. Confirm a real solve passes and deliberately bad output fails; record tolerances.
For cuPDLP-C's four-case suite, add a model in tests/data/ and an entry in
tests/cases.json, then register it in tests/CMakeLists.txt. Compressed input is
generated temporarily; no external benchmark downloads are required.

These tests do not benchmark performance or validate optional language interfaces.
