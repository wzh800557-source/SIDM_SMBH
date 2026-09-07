#!/usr/bin/env python3
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    tests = sorted((ROOT / "fp_solver" / "tests").glob("test_*.py"))
    tests += sorted((ROOT / "parameter_scan" / "tests").glob("test_*.py"))
    env = dict(os.environ)
    env["MPLBACKEND"] = "Agg"
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env["PYTHONPATH"] = os.pathsep.join(
        [str(ROOT / "fp_solver"), str(ROOT / "parameter_scan"), env.get("PYTHONPATH", "")]
    )
    for test in tests:
        print(f"==> {test.relative_to(ROOT)}", flush=True)
        subprocess.run([sys.executable, str(test)], cwd=ROOT, env=env, check=True)
    print(f"PASS: {len(tests)} test programs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
