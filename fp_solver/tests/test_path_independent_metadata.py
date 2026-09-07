#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path


def main() -> int:
    package = Path(__file__).resolve().parents[1]
    scripts = (
        "prepare_bh_remapped_profile.py",
        "hydrostatic_bridge.py",
        "boundary_criterion_scan.py",
        "finite_angle_capture.py",
        "combine_ej_operators.py",
        "solve_ej_steady_state.py",
        "assess_ej_convergence.py",
        "assemble_absolute_closure.py",
        "measured_fluid_feedback.py",
        "analyze_fluid_response.py",
    )
    for name in scripts:
        text = (package / name).read_text()
        assert ".resolve()" not in text, name
    bridge = (package / "hydrostatic_bridge.py").read_text()
    assert '"input_profile": provenance_name(args.profile)' in bridge
    assert '"output_profile": provenance_name(args.out_profile)' in bridge
    prep = (package.parent / "slurm" / "fp_solver" / "run_production_ej_prep.sbatch").read_text()
    assert "'root':root.name" in prep
    assert "'package':pkg.name" in prep
    combine = (package / "combine_ej_operators.py").read_text()
    assert '"input_paths": [path.name for path in paths]' in combine
    assert '"input_paths": [str(path) for path in paths]' not in combine
    print("PASS: production metadata are independent of execution-root paths")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
