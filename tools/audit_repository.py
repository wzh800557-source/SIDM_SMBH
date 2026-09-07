#!/usr/bin/env python3
from __future__ import annotations

import ast
import json
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FORBIDDEN = ("/" + "Users/", "/home/" + "wzh557", "orcd" + "/pool")


def main() -> int:
    readme = ROOT / "README.md"
    assert readme.is_file() and readme.stat().st_size == 0
    excluded = {".git", ".venv", "venv", "build", "dist", "results", "logs"}
    files = [path for path in ROOT.rglob("*") if path.is_file()
             and not excluded.intersection(path.relative_to(ROOT).parts)]
    assert files
    for path in files:
        relative = path.relative_to(ROOT)
        assert "__pycache__" not in relative.parts
        assert path.name != ".DS_Store"
        assert path.stat().st_size < 5_000_000, relative
        if path.suffix == ".py":
            text = path.read_text()
            ast.parse(text, filename=str(relative))
            assert not any(item in text for item in FORBIDDEN), relative
        elif path.suffix in {".sh", ".sbatch"}:
            text = path.read_text()
            assert not any(item in text for item in FORBIDDEN), relative
            subprocess.run(["bash", "-n", str(path)], check=True)
        elif path.suffix == ".json":
            json.loads(path.read_text())
    print(f"PASS: repository audit ({len(files)} files)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
