#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Run the black-hole-aware fluid response with the bundled fluid engine."""
from pathlib import Path
import runpy
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "fp_solver"))
if __name__ == "__main__":
    runpy.run_path(str(ROOT / "fp_solver" / "measured_fluid_feedback.py"), run_name="__main__")
