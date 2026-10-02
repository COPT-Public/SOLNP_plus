"""Run the SOLNP+ installation test with independent numerical checks."""
import argparse
import datetime
import json
import os
import platform
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--project', default='SOLNP_plus', choices=['SOLNP_plus'])
parser.add_argument('--binary', type=Path)
parser.add_argument('--output', type=Path, default=Path('tests/logs/tiny'))
args = parser.parse_args()

root = Path.cwd()
output_root = args.output.resolve()
output_root.mkdir(parents=True, exist_ok=True)
run_dir = Path(tempfile.mkdtemp(prefix='run-', dir=output_root))
(output_root / 'verification.json').write_text(
    json.dumps({'pass': False, 'state': 'running', 'run_directory': str(run_dir)})
)

env = os.environ.copy()
env.update(
    OMP_NUM_THREADS='1',
    OPENBLAS_NUM_THREADS='1',
    MKL_NUM_THREADS='1',
    MKL_THREADING_LAYER='SEQUENTIAL',
)


def execute(command):
    (run_dir / 'command.json').write_text(json.dumps(command, indent=2))
    with (run_dir / 'solver.log').open('w') as log:
        process = subprocess.run(
            command,
            cwd=root,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            timeout=240,
        )
    return process.returncode, (run_dir / 'solver.log').read_text(errors='replace')


def emit(report):
    report['run_directory'] = str(run_dir)
    encoded = json.dumps(report, indent=2, allow_nan=False)
    (run_dir / 'verification.json').write_text(encoded)
    (output_root / 'verification.json').write_text(encoded)
    print(encoded, flush=True)


report = {
    'project': 'SOLNP_plus',
    'timestamp_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'host': platform.node(),
    'primal_feasibility_tolerance': 1e-6,
    'pass': False,
}

try:
    binary = args.binary or root / 'build-coinor/test_nlp'
    return_code, log = execute([str(binary), str(root / 'examples/problem.txt')])
    records = [
        line.split('COINOR_RESULT ', 1)[1]
        for line in log.splitlines()
        if line.startswith('COINOR_RESULT ')
    ]
    result = json.loads(records[-1])
    x = np.array(result['x'])
    objective = float((x[0] - 1) ** 2 + (x[1] - 2) ** 2)
    residual = float(abs(x.sum() - 3))
    passed = (
        return_code == 0
        and result['status'] == 1
        and np.isfinite(x).all()
        and residual <= 1e-6
        and objective <= 1e-6
        and np.max(np.abs(x - [1, 2])) <= 1e-3
    )
    result.update(objective=objective, primal_residual=residual, solver_exit=return_code)
    result['pass'] = bool(passed)
    report.update(result)
except Exception as error:
    report['error'] = f'{type(error).__name__}: {error}'

emit(report)
sys.exit(0 if report['pass'] else 1)
