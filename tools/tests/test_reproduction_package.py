#!/usr/bin/env python3
from pathlib import Path
import importlib.util
import subprocess
import sys

root = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('verify_inputs', root / 'tools/verify_reproduction_inputs.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
report = module.verify()
assert report['status'] == 'INPUTS_AND_RUN_PLAN_PASS'
result = subprocess.run([sys.executable, str(root / 'tools/reproduce_production.py'), '--help'], capture_output=True, text=True)
assert result.returncode == 0, result.stderr
assert 'operator' in result.stdout and 'aggregate' in result.stdout
print('PASS: released production inputs and 38-task / 1200-case run plans')
