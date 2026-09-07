#!/usr/bin/env python3
"""Report dense-scan progress without waiting for the aggregate job."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path


def status(path: Path) -> str:
    if not path.is_file():
        return "MISSING"
    try:
        return str(json.loads(path.read_text()).get("status", "UNKNOWN"))
    except (OSError, json.JSONDecodeError):
        return "UNREADABLE"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid", type=Path, required=True)
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    grid = json.loads(args.grid.read_text())
    prep: Counter[str] = Counter()
    fluid: Counter[str] = Counter()
    for case in grid["cases"]:
        root = args.root / "cases" / case["tag"]
        prep[status(root / "PREP_STATUS.json")] += 1
        fluid[status(root / "all_halo_flux" / "STATUS.json")] += 1
    print(
        json.dumps(
            {
                "case_count": grid["case_count"],
                "preparation": dict(prep),
                "fluid": dict(fluid),
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
