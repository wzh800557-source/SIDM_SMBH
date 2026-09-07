#!/usr/bin/env python3
"""Construct and preflight an exact-r_in, physically normalized GNC run deck."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import shutil
from pathlib import Path
from typing import Iterable

import numpy as np


def replace_assignment(text: str, key: str, value: str) -> str:
    pattern = re.compile(
        rf"^(?P<prefix>\s*{re.escape(key)}\s*=\s*).*$", re.MULTILINE | re.IGNORECASE
    )
    matches = list(pattern.finditer(text))
    if len(matches) != 1:
        raise ValueError(f"expected one assignment for {key!r}, found {len(matches)}")
    return pattern.sub(lambda m: m.group("prefix") + value, text, count=1)


def replace_first_data_line(text: str, value: str) -> str:
    lines = text.splitlines()
    for i, line in enumerate(lines):
        s = line.strip()
        if s and not s.startswith("#"):
            if not re.fullmatch(r"[+\-]?\d+", s):
                raise ValueError(f"first model data line is not the MPI task count: {line!r}")
            lines[i] = value
            return "\n".join(lines) + "\n"
    raise ValueError("model has no data line")


def replace_grid_pair(text: str, gx: int, dc: int) -> str:
    lines = text.splitlines()
    marker = None
    for i, line in enumerate(lines):
        if "gx_bins dc_bins" in line:
            marker = i
            break
    if marker is None:
        raise ValueError("could not locate gx_bins/dc_bins marker")
    for i in range(marker + 1, len(lines)):
        s = lines[i].strip()
        if s and not s.startswith("#"):
            vals = s.split()
            if len(vals) < 2:
                raise ValueError("malformed gx/dc grid line")
            lines[i] = f"{gx:d}  {dc:d}"
            return "\n".join(lines) + "\n"
    raise ValueError("no grid line follows gx_bins/dc_bins marker")


def assignment_float(text: str, key: str) -> float:
    pattern = re.compile(
        rf"^\s*{re.escape(key)}\s*=\s*([^#\s]+)", re.MULTILINE | re.IGNORECASE
    )
    m = pattern.search(text)
    if not m:
        raise ValueError(f"missing assignment {key}")
    return float(re.sub(r"[dD]", "e", m.group(1)))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_profile_norm(path: Path) -> dict:
    values = np.loadtxt(path, comments="#", max_rows=1)
    if values.size != 8:
        raise ValueError("profile_norm.in must contain exactly eight numeric values")
    keys = ("rh", "n0", "rb", "mpart", "beta", "rhob", "sigb", "asymp")
    return dict(zip(keys, map(float, np.ravel(values))))


def main(argv: Iterable[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--base-model", type=Path, required=True)
    p.add_argument("--normalized-dir", type=Path, required=True)
    p.add_argument("--run-dir", type=Path, required=True)
    p.add_argument("--ranks", type=int, default=12)
    p.add_argument("--gx-bins", type=int, default=48)
    p.add_argument("--dc-bins", type=int, default=72)
    p.add_argument("--dt-tnr", type=float, default=0.1)
    p.add_argument("--total-tnr", type=float, default=6.0)
    p.add_argument("--updates-per-snapshot", type=int, default=10)
    p.add_argument(
        "--sidm-kernel", choices=("yukawa-tchannel",),
        default="yukawa-tchannel",
    )
    p.add_argument(
        "--sigma-over-m", type=float, default=100.0,
        help="zero-velocity total cross section sigma0/m [cm^2/g]",
    )
    p.add_argument("--w-kms", type=float, default=80.0)
    args = p.parse_args(argv)

    if args.sigma_over_m <= 0.0 or args.w_kms <= 0.0:
        raise ValueError("SIDM kernel parameters must be positive")

    if args.dc_bins % args.ranks:
        raise ValueError("dc-bins must be divisible by MPI ranks")
    nsnap_float = args.total_tnr / args.dt_tnr
    nsnap = int(round(nsnap_float))
    if not math.isclose(nsnap_float, nsnap, rel_tol=0.0, abs_tol=1e-10):
        raise ValueError("total-tnr must be an integer multiple of dt-tnr")

    nd = args.normalized_dir
    required = [
        nd / "df_normalization.json",
        nd / "profile_norm.in",
        nd / "mfrac.normalized.in",
        nd / "f_evolved_normalized.txt",
    ]
    required += [nd / f"initial_xj_rank{k:04d}.dat" for k in range(args.ranks)]
    missing = [str(x) for x in required if not x.is_file()]
    if missing:
        raise FileNotFoundError("missing normalized inputs: " + ", ".join(missing))
    dfdiag = json.loads((nd / "df_normalization.json").read_text())
    if dfdiag.get("status") != "PASS":
        raise ValueError(f"DF normalization did not pass: {dfdiag.get('status')}")
    norm = read_profile_norm(nd / "profile_norm.in")
    xb = norm["rh"] / (2.0 * norm["rb"])
    xmin = float(dfdiag["xmin"])
    xmax = float(dfdiag["xmax"])
    if not xmin < xb < xmax:
        raise ValueError("exact boundary energy is not inside normalized DF grid")
    if dfdiag.get("capture_species") != "SBH point-particle proxy":
        raise ValueError("normalized DF is not configured for relativistic point-particle capture")
    if not math.isclose(
        norm["asymp"],
        float(dfdiag["asymp_g_at_x_boundary"]),
        rel_tol=5.0e-13,
        abs_tol=0.0,
    ):
        raise ValueError("profile_norm.in and DF diagnostics disagree on g(x_boundary)")

    # The existing cluster loader skips exactly one header record and then reads
    # ``nx nj``.  Validate that byte-level contract before submitting an MPI job.
    table_lines = (nd / "f_evolved_normalized.txt").read_text().splitlines()
    if len(table_lines) < 4 or not table_lines[0].startswith("#"):
        raise ValueError("normalized DF table has no single header record")
    try:
        nx_table, nj_table = map(int, table_lines[1].split())
    except Exception as exc:
        raise ValueError("DF loader cannot read nx,nj from the second record") from exc
    if (nx_table, nj_table) != (args.gx_bins, args.gx_bins):
        raise ValueError(
            f"DF table grid {(nx_table, nj_table)} does not match gx-bins={args.gx_bins}"
        )

    model = args.base_model.read_text()
    model = replace_first_data_line(model, str(args.ranks))
    model = replace_assignment(model, "mbh", f"{float(dfdiag['mbh_msun']):.16e}")
    model = replace_assignment(model, "emin_factor", f"{xmin:.16e}")
    model = replace_assignment(model, "emax_factor", f"{xmax:.16e}")
    model = replace_assignment(model, "eboundary", f"{xb:.16e}")
    model = replace_assignment(
        model, "num of ge update per snap", str(args.updates_per_snapshot)
    )
    model = replace_assignment(model, "timestep_snapshot (output)", f"{args.dt_tnr:.16e}")
    model = replace_assignment(model, "num_of_snapshot", str(nsnap))
    model = replace_assignment(model, "clone x0", f"{10.0 * xb:.16e}")
    model = replace_grid_pair(model, args.gx_bins, args.dc_bins)

    # Deck round-trip checks catch the sequential-parser failure that previously
    # turned a blank clone value into -Infinity/NaN later in initialization.
    parsed = {
        key: assignment_float(model, key)
        for key in ("mbh", "emin_factor", "emax_factor", "eboundary", "clone x0")
    }
    if not all(np.isfinite(list(parsed.values()))):
        raise ValueError("model deck contains non-finite numeric values")
    r_from_deck = norm["rh"] / (2.0 * parsed["eboundary"])
    if abs(r_from_deck / norm["rb"] - 1.0) > 1.0e-12:
        raise ValueError("deck does not map back to exact r_in")
    if not math.isclose(parsed["clone x0"] / 10.0, xb, rel_tol=1e-12):
        raise ValueError("clone threshold is not the boundary energy")
    guard_cells = (
        (math.log10(xb) - math.log10(xmin))
        / (math.log10(xmax) - math.log10(xmin))
        * (args.gx_bins - 1)
    )
    if guard_cells < 2.0:
        raise ValueError("fewer than two low-energy guard cells")

    # Check every rank's input count and strict grid containment.  The proposal
    # weights describe one distributed Monte-Carlo ensemble, so their unit-mean
    # constraint is global rather than rank-local.
    sample_summary = []
    total_weight = 0.0
    total_weight_squared = 0.0
    total_samples = 0
    for rank in range(args.ranks):
        fp = nd / f"initial_xj_rank{rank:04d}.dat"
        with fp.open() as f:
            n = int(f.readline().strip())
        arr = np.loadtxt(fp, skiprows=1)
        if arr.shape != (n, 3):
            raise ValueError(f"rank {rank}: declared {n} samples, read {arr.shape}")
        if not (
            np.all(np.isfinite(arr))
            and np.all((arr[:, 0] > xmin) & (arr[:, 0] < xmax))
            and np.all((arr[:, 1] > float(dfdiag["jmin"])) & (arr[:, 1] < float(dfdiag["jmax"])))
            and np.all(arr[:, 2] > 0.0)
        ):
            raise ValueError(f"rank {rank}: invalid x, j, or importance weight")
        rank_weight_sum = float(math.fsum(float(w) for w in arr[:, 2]))
        rank_weight_mean = rank_weight_sum / n
        total_weight += rank_weight_sum
        total_weight_squared += float(np.sum(arr[:, 2] ** 2))
        total_samples += n
        sample_summary.append(
            {
                "rank": rank,
                "n": n,
                "active": int(np.count_nonzero(arr[:, 0] >= xb)),
                "importance_weight_sum": rank_weight_sum,
                "importance_weight_mean": rank_weight_mean,
            }
        )
    global_weight_mean = total_weight / total_samples
    if abs(global_weight_mean - 1.0) > 5.0e-10:
        raise ValueError(
            "importance weights do not have unit global mean: "
            f"{global_weight_mean:.16e}"
        )
    global_effective_sample_size = total_weight**2 / total_weight_squared

    # This fingerprint excludes boundary-specific metadata and g(x_boundary).
    # Runs with the same numerical DF table and the same weighted sample ensemble
    # must therefore have the same value even when their reporting boundary moves.
    common_reservoir_hash = hashlib.sha256()
    for line in table_lines[1:]:
        common_reservoir_hash.update(line.encode("ascii"))
        common_reservoir_hash.update(b"\n")
    for rank in range(args.ranks):
        common_reservoir_hash.update(
            (nd / f"initial_xj_rank{rank:04d}.dat").read_bytes()
        )
    common_reservoir_fingerprint = common_reservoir_hash.hexdigest()

    args.run_dir.mkdir(parents=True, exist_ok=True)
    (args.run_dir / "model.in").write_text(model)
    (args.run_dir / "sidm_kernel.in").write_text(
        f"{args.sigma_over_m:.16e} {args.w_kms:.16e}\n"
    )
    shutil.copy2(nd / "mfrac.normalized.in", args.run_dir / "mfrac.in")
    for src in required:
        if src.name == "mfrac.normalized.in":
            continue
        shutil.copy2(src, args.run_dir / src.name)
    # ``main/ini.f90`` on the production tree calls load_df_from_table with this
    # fixed filename before the weighted samples reconstruct the same DF.
    shutil.copy2(
        nd / "f_evolved_normalized.txt", args.run_dir / "f_evolved.txt"
    )

    manifest = {
        "schema": "gnc-absolute-run-v2",
        "status": "PREFLIGHT_PASS",
        "exact_boundary_relation": "x_boundary = R_h / (2 r_in)",
        "rh_pc": norm["rh"],
        "r_in_pc": norm["rb"],
        "x_boundary": xb,
        "normalization_radius_pc": float(
            dfdiag.get("normalization_radius_pc", norm["rb"])
        ),
        "normalization_x": float(dfdiag.get("normalization_x", xb)),
        "proposal_boundary_radius_pc": float(
            dfdiag.get("proposal_boundary_radius_pc", norm["rb"])
        ),
        "proposal_x_boundary": float(dfdiag.get("proposal_x_boundary", xb)),
        "common_reservoir_fingerprint": common_reservoir_fingerprint,
        "emin_factor": xmin,
        "emax_factor": xmax,
        "guard_cells_below_boundary": guard_cells,
        "clone_x0_deck": 10.0 * xb,
        "df_asymp": norm["asymp"],
        "cusp_slope_beta": dfdiag["beta"],
        "n0_pc3": norm["n0"],
        "rho_boundary_msun_pc3": norm["rhob"],
        "sigma_boundary_kms": norm["sigb"],
        "particle_mass_msun": norm["mpart"],
        "capture_species": dfdiag["capture_species"],
        "capture_radius_model": dfdiag["capture_radius_model"],
        "sidm_kernel": args.sidm_kernel,
        "sigma0_over_m_cm2_g": args.sigma_over_m,
        "sigma0_convention": "zero-velocity total cross section in the t-channel Born kernel",
        "yukawa_w_kms": args.w_kms,
        "ranks": args.ranks,
        "gx_bins": args.gx_bins,
        "dc_bins": args.dc_bins,
        "dt_tnr": args.dt_tnr,
        "snapshots": nsnap,
        "total_tnr": nsnap * args.dt_tnr,
        "samples": sample_summary,
        "importance_weight_global_mean": global_weight_mean,
        "importance_weight_global_effective_size": global_effective_sample_size,
        "importance_weight_normalization": "one global mean over all MPI ranks",
        "files": {
            p.name: sha256(p)
            for p in sorted(args.run_dir.iterdir())
            if p.is_file() and p.name != "run_manifest.json"
        },
    }
    (args.run_dir / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
