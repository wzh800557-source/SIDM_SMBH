"""Locate the fluid dependency without adding installation paths to run metadata."""
from __future__ import annotations

import importlib
import os
from pathlib import Path
import sys


def load_fluid_engine():
    bundled = Path(__file__).resolve().parents[1] / "fluid_solver" / "GravothermalSIDM"
    root = Path(os.environ.get("GRAVOTHERMAL_ROOT") or bundled).expanduser().resolve()
    for name in ("evolve.py", "record.py"):
        if not (root / "SourcePy" / name).is_file():
            raise ImportError("Fluid engine missing. Restore fluid_solver/GravothermalSIDM or set GRAVOTHERMAL_ROOT.")
    sys.path.insert(0, str(root))
    modules = tuple(importlib.import_module(f"SourcePy.{name}") for name in ("evolve", "record"))
    for module in modules:
        if root not in Path(module.__file__).resolve().parents:
            raise ImportError("A different SourcePy engine was already imported. Use a fresh Python process.")
    return modules
