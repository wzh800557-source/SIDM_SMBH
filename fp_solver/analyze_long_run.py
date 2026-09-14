#!/usr/bin/env python3
"""Measure an absolute GNC capture closure from a multi-relaxation-time run.

The weighted event table supplies the physical mass captured per output interval.
The patched ``pro`` executable also writes every compact-object plunge, retaining
its final and initial orbital energy and its native Monte-Carlo weight.  Combining
these two records gives a directly measured mass current, capture-weighted energy
current, and the dimensionless coefficients required by the fluid boundary.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Iterable, Tuple

import numpy as np

G_PC_KMS2_MSUN = 4.30091e-3
KMS_TO_PCMYR = 1.0227121650537077
CM2_G_TO_PC2_MSUN = 2.0884205246637706e-4


def read_named_table(path: Path) -> np.ndarray:
    arr = np.genfromtxt(path, names=True, dtype=float, autostrip=True)
    if arr.dtype.names is None:
        raise ValueError(f"no named columns in {path}")
    return np.atleast_1d(arr)


def find_event_table(root: Path, weighted: bool) -> Path:
    suffix = "_event_Nweight.txt" if weighted else "_event_N.txt"
    candidates = sorted(root.rglob(f"*{suffix}"))
    if not candidates:
        raise FileNotFoundError(f"no {suffix} beneath {root}")
    for species in ("BH", "NS", "WD", "BD", "MS"):
        preferred = [p for p in candidates if f"/{species}/" in p.as_posix()]
        if preferred:
            return preferred[0]
    return candidates[0]


def get_column(arr: np.ndarray, alternatives: Tuple[str, ...]) -> np.ndarray:
    names = {n.lower(): n for n in arr.dtype.names or ()}
    for alt in alternatives:
        if alt.lower() in names:
            return np.asarray(arr[names[alt.lower()]], float)
    raise KeyError(f"none of {alternatives} found in columns {arr.dtype.names}")


def relative_change(a: float, b: float) -> float:
    scale = max(abs(a), abs(b), np.finfo(float).tiny)
    return abs(a - b) / scale


def load_manifest(run: Path) -> dict:
    p = run / "run_manifest.json"
    if not p.is_file():
        raise FileNotFoundError(p)
    out = json.loads(p.read_text())
    if out.get("status") != "PREFLIGHT_PASS":
        raise ValueError("run manifest did not pass preflight")
    return out


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_plunge_records(run: Path) -> np.ndarray:
    files = sorted(
        run.rglob("plunge_records_*.txt"),
        key=lambda p: int(re.findall(r"(\d+)", p.stem)[-1]),
    )
    bh_files = [path for path in files if "/BH/" in path.as_posix()]
    if bh_files:
        files = bh_files
    rows = []
    for p in files:
        data_lines = [
            line for line in p.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
        if not data_lines:
            continue
        a = np.loadtxt(data_lines, ndmin=2)
        if a.shape[1] not in (9, 10):
            raise ValueError(f"{p} has {a.shape[1]} plunge columns, expected 9 or 10")
        if a.shape[1] == 9:
            # Older patched runs did not record particle creation time.  Keep
            # them readable, but the production gate records this provenance.
            a = np.column_stack([a, np.full(a.shape[0], np.nan)])
        rows.append(a)
    return np.vstack(rows) if rows else np.empty((0, 10), float)


def load_inner_inventory(run: Path) -> np.ndarray:
    """Load per-snapshot mass and binding-energy inventories from patched GNC."""

    files = sorted(
        run.rglob("inner_inventory_*.txt"),
        key=lambda p: int(re.findall(r"(\d+)", p.stem)[-1]),
    )
    bh_files = [path for path in files if "/BH/" in path.as_posix()]
    if bh_files:
        files = bh_files
    rows = []
    for path in files:
        data = np.loadtxt(path, comments="#", ndmin=2)
        if data.shape != (1, 6):
            raise ValueError(f"{path} has shape {data.shape}, expected one row and six columns")
        rows.append(data[0])
    return np.asarray(rows, dtype=float) if rows else np.empty((0, 6), float)


def effective_count(weights: np.ndarray) -> float:
    weights = np.asarray(weights, dtype=float)
    total = float(weights.sum())
    squared = float(np.square(weights).sum())
    return total * total / squared if squared > 0.0 else 0.0


def weighted_moments(values: np.ndarray, weights: np.ndarray) -> dict:
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    if values.size == 0 or weights.size != values.size or float(weights.sum()) <= 0.0:
        return {
            "mean": math.nan,
            "variance": math.nan,
            "standard_error": math.nan,
            "weight_sum": 0.0,
            "weight_squared_sum": 0.0,
            "weighted_value_sum": 0.0,
            "weighted_value_squared_sum": 0.0,
            "effective_count": 0.0,
        }
    mean = float(np.average(values, weights=weights))
    variance = float(np.average(np.square(values - mean), weights=weights))
    neff = effective_count(weights)
    return {
        "mean": mean,
        "variance": variance,
        "standard_error": math.sqrt(variance / neff) if neff > 0.0 else math.nan,
        "weight_sum": float(weights.sum()),
        "weight_squared_sum": float(np.square(weights).sum()),
        "weighted_value_sum": float(np.dot(weights, values)),
        "weighted_value_squared_sum": float(np.dot(weights, np.square(values))),
        "effective_count": neff,
    }


def half_plateau(numerators: np.ndarray, durations: np.ndarray) -> dict:
    """Compare time-averaged currents in the first and second retained halves."""

    numerators = np.asarray(numerators, dtype=float)
    durations = np.asarray(durations, dtype=float)
    midpoint = numerators.size // 2
    if midpoint < 1 or numerators.size - midpoint < 1:
        return {
            "available": False,
            "first_half": math.nan,
            "second_half": math.nan,
            "relative_change": math.inf,
        }
    first = float(numerators[:midpoint].sum() / durations[:midpoint].sum())
    second = float(numerators[midpoint:].sum() / durations[midpoint:].sum())
    return {
        "available": True,
        "first_half": first,
        "second_half": second,
        "relative_change": relative_change(first, second),
    }


def main(argv: Iterable[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("run", type=Path)
    p.add_argument("--dt-tnr", type=float, default=0.1)
    p.add_argument("--window-tnr", type=float, default=1.0)
    p.add_argument("--burn-in-tnr", type=float, default=2.0)
    p.add_argument("--plateau-tolerance", type=float, default=0.25)
    p.add_argument("--minimum-raw-captures", type=int, default=20)
    p.add_argument("--minimum-effective-captures", type=float, default=50.0)
    p.add_argument(
        "--minimum-boundary-effective-captures",
        type=float,
        default=10.0,
        help="per-run diagnostic floor; the ensemble gate applies its own pooled floor",
    )
    p.add_argument(
        "--inventory-tolerance",
        type=float,
        default=0.25,
        help="maximum first-half/second-half change in inner mass and binding inventory",
    )
    time_group = p.add_mutually_exclusive_group()
    time_group.add_argument(
        "--physical-tnr-myr", type=float, default=None,
        help="physical duration assigned to one dimensionless GNC TNR"
    )
    time_group.add_argument(
        "--sigma-over-m-cm2-g", type=float, default=None,
        help=("derive the physical TNR from the orbital scattering time at "
              "R_h: [rho(R_h) (sigma/m) sqrt(G M_bh/R_h)]^-1")
    )
    p.add_argument(
        "--kernel-compatible", action="store_true",
        help=("assert that the SIDM differential-scattering kernel has the same "
              "dimensionless velocity and angular dependence as the GNC tables")
    )
    # The native event summary is written with limited decimal precision.  A
    # one-per-mille gate is still much tighter than any statistical error while
    # not rejecting the direct plunge ledger solely because of text rounding.
    p.add_argument("--mass-consistency-tolerance", type=float, default=1.0e-3)
    p.add_argument("--out-json", type=Path, default=None)
    p.add_argument("--out-csv", type=Path, default=None)
    p.add_argument("--out-plunge-csv", type=Path, default=None)
    args = p.parse_args(argv)

    manifest = load_manifest(args.run)
    build_status_path = args.run / "build_status.json"
    build_status = (
        json.loads(build_status_path.read_text())
        if build_status_path.is_file() else {}
    )
    required_build_flags = (
        "compiled",
        "absolute_normalization_patch",
        "plunge_records_patch",
        "inner_inventory_patch",
        "snapshot_terminal_coefficients_only",
        "weighted_xj_loader",
        "born_kernel",
    )
    validated_build = (
        build_status.get("status") == "COMPILED"
        and all(build_status.get(key) is True for key in required_build_flags)
    )
    executable_hashes = {
        name: sha256(args.run / name)
        for name in ("ini", "main", "pro")
        if (args.run / name).is_file()
    }
    solver_fingerprint = (
        hashlib.sha256(
            json.dumps(executable_hashes, sort_keys=True).encode("utf-8")
        ).hexdigest()
        if len(executable_hashes) == 3 else None
    )
    ranks = int(manifest["ranks"])
    particle_mass = float(manifest["particle_mass_msun"])
    mbh = float(json.loads((args.run / "df_normalization.json").read_text())["mbh_msun"])
    rh = float(manifest["rh_pc"])
    rb = float(manifest["r_in_pc"])
    rho_b = float(manifest["rho_boundary_msun_pc3"])
    sigma_b = float(manifest["sigma_boundary_kms"])
    rho_h = float(manifest["n0_pc3"]) * particle_mass

    weighted_path = find_event_table(args.run, True)
    raw_path = find_event_table(args.run, False)
    weighted = read_named_table(weighted_path)
    raw = read_named_table(raw_path)
    if weighted.size != raw.size:
        raise ValueError("weighted and raw event tables have different lengths")

    capture_names = ("N_plunge", "N_td", "N_tdfull", "N_tdempty", "N_tmpfull")
    cap_w = get_column(weighted, capture_names)
    cap_n_average = get_column(raw, capture_names)
    boundary_w = get_column(weighted, ("N_emax",))
    boundary_n_average = get_column(raw, ("N_emax",))
    inner_inventory_event_w = get_column(weighted, ("N_norm_bd",))
    inner_inventory_event_n_average = get_column(raw, ("N_norm_bd",))
    tsnap_gnc_native_myr = get_column(weighted, ("Tsnap", "TsnapMyr"))
    if not np.all(np.isfinite(np.r_[cap_w, cap_n_average, boundary_w,
                                    boundary_n_average, inner_inventory_event_w,
                                    inner_inventory_event_n_average,
                                    tsnap_gnc_native_myr])):
        raise ValueError("event table contains non-finite values")
    if np.any(np.diff(tsnap_gnc_native_myr) <= 0.0):
        raise ValueError("GNC snapshot times are not strictly increasing")

    spw_float = args.window_tnr / args.dt_tnr
    spw = int(round(spw_float))
    if not math.isclose(spw_float, spw, abs_tol=1.0e-10):
        raise ValueError("window-tnr must be an integer multiple of dt-tnr")
    nwin = cap_w.size // spw
    if nwin < 2:
        raise ValueError("run is shorter than two analysis windows")
    nuse = nwin * spw
    cap_w = cap_w[:nuse]
    cap_n_average = cap_n_average[:nuse]
    boundary_w = boundary_w[:nuse]
    boundary_n_average = boundary_n_average[:nuse]
    inner_inventory_event_w = inner_inventory_event_w[:nuse]
    inner_inventory_event_n_average = inner_inventory_event_n_average[:nuse]
    tsnap_gnc_native_myr = tsnap_gnc_native_myr[:nuse]

    gnc_native_tnr_myr = float(
        tsnap_gnc_native_myr[-1] / (nuse * args.dt_tnr)
    )
    vh_kms = math.sqrt(G_PC_KMS2_MSUN * mbh / rh)
    if args.physical_tnr_myr is not None:
        if not np.isfinite(args.physical_tnr_myr) or args.physical_tnr_myr <= 0.0:
            raise ValueError("physical-tnr-myr must be finite and positive")
        physical_tnr_myr = float(args.physical_tnr_myr)
        time_mapping = "explicit physical TNR supplied by caller"
    elif args.sigma_over_m_cm2_g is not None:
        if not np.isfinite(args.sigma_over_m_cm2_g) or args.sigma_over_m_cm2_g <= 0.0:
            raise ValueError("sigma-over-m-cm2-g must be finite and positive")
        sigma_pc2_msun = args.sigma_over_m_cm2_g * CM2_G_TO_PC2_MSUN
        physical_tnr_myr = 1.0 / (
            rho_h * sigma_pc2_msun * vh_kms * KMS_TO_PCMYR
        )
        time_mapping = (
            "SIDM orbital scattering time at R_h, using v_h=sqrt(G M_bh/R_h)"
        )
    else:
        physical_tnr_myr = gnc_native_tnr_myr
        time_mapping = "native relaxation clock recorded by the GNC executable"

    tsnap_myr = (
        np.arange(1, nuse + 1, dtype=float) * args.dt_tnr * physical_tnr_myr
    )
    dt_snapshot_myr = np.full(nuse, args.dt_tnr * physical_tnr_myr)
    cap_w_win = cap_w.reshape(nwin, spw).sum(axis=1)
    cap_n_win_actual = cap_n_average.reshape(nwin, spw).sum(axis=1) * ranks
    bd_w_win = boundary_w.reshape(nwin, spw).sum(axis=1)
    bd_n_win_actual = boundary_n_average.reshape(nwin, spw).sum(axis=1) * ranks
    dt_win_myr = dt_snapshot_myr.reshape(nwin, spw).sum(axis=1)
    t1_myr = tsnap_myr[spw - 1:nuse:spw]
    t0_myr = t1_myr - dt_win_myr
    t0_tnr = np.arange(nwin, dtype=float) * args.window_tnr
    t1_tnr = t0_tnr + args.window_tnr
    keep = t0_tnr >= args.burn_in_tnr - 1.0e-12
    if np.count_nonzero(keep) < 2:
        raise ValueError("need at least two complete post-burn-in windows")

    mass_win = cap_w_win * particle_mass
    mdot_win = mass_win / dt_win_myr
    post_time_myr = float(dt_win_myr[keep].sum())
    post_mass_table = float(mass_win[keep].sum())
    mdot = post_mass_table / post_time_myr
    last_change = relative_change(float(mdot_win[-2]), float(mdot_win[-1]))
    post_window_rates = np.asarray(mdot_win[keep], float)
    window_rate_mean = float(np.mean(post_window_rates))
    window_rate_std = float(np.std(post_window_rates, ddof=1)) if post_window_rates.size > 1 else 0.0
    window_rate_sem = window_rate_std / math.sqrt(post_window_rates.size)

    plunge = load_plunge_records(args.run)
    # columns: isnap, x_final, x_initial, j_final, j_initial, weight_real,
    #          particle mass, exit time, periapse, creation time
    burn_snapshot = int(round(args.burn_in_tnr / args.dt_tnr))
    if plunge.size:
        in_range = (plunge[:, 0] > burn_snapshot) & (plunge[:, 0] <= nuse)
        post_plunge = plunge[in_range]
    else:
        post_plunge = plunge

    energy_unit_kms2 = G_PC_KMS2_MSUN * mbh / rh
    mass_flux_scale = 4.0 * math.pi * rb**2 * rho_b * sigma_b * KMS_TO_PCMYR
    energy_flux_scale = mass_flux_scale * sigma_b**2
    x_boundary = float(manifest["x_boundary"])
    boundary_binding_kms2 = x_boundary * energy_unit_kms2

    if post_plunge.size:
        physical_number_weight = post_plunge[:, 5] / ranks
        physical_mass_weight = physical_number_weight * post_plunge[:, 6]
        x_final_binding = np.abs(post_plunge[:, 1])
        x_initial_binding = np.abs(post_plunge[:, 2])
        all_final = weighted_moments(x_final_binding, physical_mass_weight)
        all_initial = weighted_moments(x_initial_binding, physical_mass_weight)
        post_mass_records = all_final["weight_sum"]
        mean_x_final = all_final["mean"]
        mean_x_initial = all_initial["mean"]
        capture_effective_count = all_final["effective_count"]
        mean_x_final_standard_error = all_final["standard_error"]
        raw_capture_records = int(post_plunge.shape[0])
    else:
        physical_mass_weight = np.empty(0, float)
        x_final_binding = np.empty(0, float)
        x_initial_binding = np.empty(0, float)
        all_final = weighted_moments(x_final_binding, physical_mass_weight)
        all_initial = weighted_moments(x_initial_binding, physical_mass_weight)
        post_mass_records = 0.0
        mean_x_final = math.nan
        mean_x_initial = math.nan
        capture_effective_count = 0.0
        mean_x_final_standard_error = math.nan
        raw_capture_records = 0

    mass_consistency = (
        abs(post_mass_records - post_mass_table)
        / max(abs(post_mass_records), abs(post_mass_table), np.finfo(float).tiny)
    )
    mass_consistency_pass = mass_consistency <= args.mass_consistency_tolerance
    legacy_plateau_pass = last_change <= args.plateau_tolerance
    raw_count_pass = raw_capture_records >= args.minimum_raw_captures
    effective_count_pass = capture_effective_count >= args.minimum_effective_captures
    count_pass = raw_count_pass and effective_count_pass

    mean_capture_binding_kms2 = mean_x_final * energy_unit_kms2
    mean_initial_binding_kms2 = mean_x_initial * energy_unit_kms2
    c_m = mdot / mass_flux_scale
    l_capture = mdot * mean_capture_binding_kms2
    l_sink = -1.5 * mdot * sigma_b**2
    l_source_ceiling = l_capture
    c_e_capture = l_capture / energy_flux_scale
    c_e_sink = l_sink / energy_flux_scale
    c_e_source_ceiling = l_source_ceiling / energy_flux_scale

    # A late plunge need not represent a boundary-fed steady current.  The
    # initial-energy label is inherited by clones and boundary replacements, so
    # x_initial < x_boundary selects trajectories supplied from outside r_in and
    # removes depletion of the initially populated inner cusp.
    supplied_mask = x_initial_binding < x_boundary
    supplied = post_plunge[supplied_mask] if post_plunge.size else post_plunge
    supplied_mass_weight = physical_mass_weight[supplied_mask]
    supplied_x_final = x_final_binding[supplied_mask]
    supplied_x_initial = x_initial_binding[supplied_mask]
    supplied_release_x = supplied_x_final - x_boundary
    supplied_final_stats = weighted_moments(supplied_x_final, supplied_mass_weight)
    supplied_initial_stats = weighted_moments(supplied_x_initial, supplied_mass_weight)
    supplied_release_stats = weighted_moments(supplied_release_x, supplied_mass_weight)

    supplied_mass_win = np.zeros(nwin, float)
    supplied_raw_win = np.zeros(nwin, float)
    supplied_capture_mass_x_win = np.zeros(nwin, float)
    supplied_release_mass_x_win = np.zeros(nwin, float)
    if supplied.size:
        supplied_window_index = ((supplied[:, 0].astype(int) - 1) // spw).astype(int)
        valid_window = (supplied_window_index >= 0) & (supplied_window_index < nwin)
        supplied_window_index = supplied_window_index[valid_window]
        supplied_weights_in_window = supplied_mass_weight[valid_window]
        supplied_final_in_window = supplied_x_final[valid_window]
        supplied_release_in_window = supplied_release_x[valid_window]
        np.add.at(supplied_mass_win, supplied_window_index, supplied_weights_in_window)
        np.add.at(supplied_raw_win, supplied_window_index, 1.0)
        np.add.at(
            supplied_capture_mass_x_win,
            supplied_window_index,
            supplied_weights_in_window * supplied_final_in_window,
        )
        np.add.at(
            supplied_release_mass_x_win,
            supplied_window_index,
            supplied_weights_in_window * supplied_release_in_window,
        )

    supplied_mdot_win = supplied_mass_win / dt_win_myr
    supplied_capture_current_win = supplied_capture_mass_x_win * energy_unit_kms2 / dt_win_myr
    supplied_boundary_advective_current_win = (
        supplied_mass_win * x_boundary * energy_unit_kms2 / dt_win_myr
    )
    supplied_outward_current_win = supplied_release_mass_x_win * energy_unit_kms2 / dt_win_myr
    supplied_time = float(dt_win_myr[keep].sum())
    supplied_mass = float(supplied_mass_win[keep].sum())
    supplied_mdot = supplied_mass / supplied_time
    supplied_capture_current = float(
        supplied_capture_mass_x_win[keep].sum() * energy_unit_kms2 / supplied_time
    )
    supplied_boundary_advective_current = float(
        supplied_mass * x_boundary * energy_unit_kms2 / supplied_time
    )
    supplied_outward_current = float(
        supplied_release_mass_x_win[keep].sum() * energy_unit_kms2 / supplied_time
    )
    supplied_identity_residual = (
        supplied_capture_current
        - supplied_boundary_advective_current
        - supplied_outward_current
    )
    supplied_identity_relative = abs(supplied_identity_residual) / max(
        abs(supplied_capture_current),
        abs(supplied_boundary_advective_current),
        abs(supplied_outward_current),
        np.finfo(float).tiny,
    )
    supplied_mass_plateau = half_plateau(supplied_mass_win[keep], dt_win_myr[keep])
    supplied_energy_plateau = half_plateau(
        supplied_release_mass_x_win[keep] * energy_unit_kms2,
        dt_win_myr[keep],
    )
    supplied_plateau_pass = (
        supplied_mass_plateau["available"]
        and supplied_energy_plateau["available"]
        and supplied_mass_plateau["relative_change"] <= args.plateau_tolerance
        and supplied_energy_plateau["relative_change"] <= args.plateau_tolerance
    )
    supplied_effective_pass = (
        supplied_final_stats["effective_count"]
        >= args.minimum_boundary_effective_captures
    )
    supplied_energy_physical = bool(
        supplied.size
        and np.all(np.isfinite(supplied_release_x))
        and np.all(supplied_release_x >= 0.0)
        and supplied_outward_current >= 0.0
    )

    # Directly track the live kinetic reservoir.  A constant current measured
    # while this inventory drains is a transient, not a stationary closure.
    inventory = load_inner_inventory(args.run)
    if inventory.size:
        inventory = inventory[
            (inventory[:, 0] > burn_snapshot) & (inventory[:, 0] <= nuse)
        ]
    expected_inventory_snapshots = max(nuse - burn_snapshot, 0)
    inventory_completeness = (
        float(inventory.shape[0]) / expected_inventory_snapshots
        if expected_inventory_snapshots else 0.0
    )
    inventory_mass_plateau = half_plateau(
        inventory[:, 2] if inventory.size else np.empty(0),
        np.ones(inventory.shape[0], float),
    )
    inventory_energy_plateau = half_plateau(
        inventory[:, 3] * energy_unit_kms2 if inventory.size else np.empty(0),
        np.ones(inventory.shape[0], float),
    )
    inventory_stationary = bool(
        inventory_completeness >= 0.95
        and inventory_mass_plateau["available"]
        and inventory_energy_plateau["available"]
        and inventory_mass_plateau["relative_change"] <= args.inventory_tolerance
        and inventory_energy_plateau["relative_change"] <= args.inventory_tolerance
    )

    # The capture records carry unequal physical weights.  Their effective count,
    # rather than the raw row count, sets the counting error.  The legacy all-
    # capture quantities remain in the output for comparison with earlier runs.
    counting_fractional_error = (
        1.0 / math.sqrt(capture_effective_count)
        if capture_effective_count > 0.0 else math.inf
    )
    window_fractional_error = (
        window_rate_sem / abs(mdot) if mdot != 0.0 else math.inf
    )
    mdot_fractional_error = max(counting_fractional_error, window_fractional_error)
    mdot_standard_error = abs(mdot) * mdot_fractional_error
    mean_capture_binding_standard_error_kms2 = (
        mean_x_final_standard_error * energy_unit_kms2
        if np.isfinite(mean_x_final_standard_error) else math.nan
    )
    capture_energy_fractional_error = math.sqrt(
        mdot_fractional_error**2
        + (mean_x_final_standard_error / mean_x_final)**2
    ) if np.isfinite(mean_x_final_standard_error) and mean_x_final != 0.0 else math.inf
    capture_energy_current_standard_error = abs(l_capture) * capture_energy_fractional_error
    sink_current_standard_error = abs(l_sink) * mdot_fractional_error

    supplied_window_rates = supplied_mdot_win[keep]
    supplied_window_mean = float(np.mean(supplied_window_rates))
    supplied_window_std = (
        float(np.std(supplied_window_rates, ddof=1))
        if supplied_window_rates.size > 1 else 0.0
    )
    supplied_window_sem = supplied_window_std / math.sqrt(supplied_window_rates.size)
    supplied_counting_fractional_error = (
        1.0 / math.sqrt(supplied_final_stats["effective_count"])
        if supplied_final_stats["effective_count"] > 0.0 else math.inf
    )
    supplied_window_fractional_error = (
        supplied_window_sem / abs(supplied_mdot)
        if supplied_mdot != 0.0 else math.inf
    )
    supplied_mdot_fractional_error = max(
        supplied_counting_fractional_error, supplied_window_fractional_error
    )
    supplied_mdot_standard_error = abs(supplied_mdot) * supplied_mdot_fractional_error
    supplied_release_mean_fractional_error = (
        supplied_release_stats["standard_error"] / abs(supplied_release_stats["mean"])
        if np.isfinite(supplied_release_stats["standard_error"])
        and supplied_release_stats["mean"] != 0.0 else math.inf
    )
    supplied_outward_fractional_error = math.sqrt(
        supplied_mdot_fractional_error**2
        + supplied_release_mean_fractional_error**2
    )
    supplied_outward_standard_error = (
        abs(supplied_outward_current) * supplied_outward_fractional_error
    )

    physical_kernel_compatible = args.kernel_compatible
    diagnostic_failures = []
    if not mass_consistency_pass:
        diagnostic_failures.append("event and plunge mass budgets disagree")
    if raw_capture_records == 0:
        diagnostic_failures.append("no post-burn-in plunge records")
    if inventory_completeness < 0.95:
        diagnostic_failures.append("inner mass-energy inventory records are incomplete")
    if not np.isfinite(supplied_outward_current):
        diagnostic_failures.append("boundary-supplied energy current is unavailable")
    if not validated_build or solver_fingerprint is None:
        diagnostic_failures.append("patched production FP build provenance is incomplete")
    status = (
        "DIRECT_FP_ENERGY_DIAGNOSTIC_COMPLETE"
        if not diagnostic_failures else "DIRECT_FP_ENERGY_DIAGNOSTIC_INCOMPLETE"
    )
    acceptance_blockers = []
    if not supplied_effective_pass:
        acceptance_blockers.append("too few effective boundary-supplied captures in this run")
    if not supplied_plateau_pass:
        acceptance_blockers.append("boundary-supplied mass or energy current has not plateaued")
    if not inventory_stationary:
        acceptance_blockers.append("inner mass or binding-energy inventory is not stationary")
    if not supplied_energy_physical:
        acceptance_blockers.append("released binding-energy current is not non-negative")
    if not physical_kernel_compatible:
        acceptance_blockers.append("SIDM kernel compatibility was not asserted")
    statistical_acceptance_pass = not acceptance_blockers and not diagnostic_failures
    physical_blockers = ([] if physical_kernel_compatible else [
        "The requested physical time mapping must use the same Born/Yukawa "
        "kernel as the GNC diffusion tables."
    ])

    rows = np.column_stack([
        t0_tnr, t1_tnr, t0_myr, t1_myr, cap_w_win, cap_n_win_actual,
        mass_win, mdot_win, bd_w_win, bd_n_win_actual,
        supplied_raw_win, supplied_mass_win, supplied_mdot_win,
        supplied_capture_current_win, supplied_boundary_advective_current_win,
        supplied_outward_current_win,
    ])
    out_csv = args.out_csv or args.run / "absolute_closure_windows.csv"
    out_json = args.out_json or args.run / "absolute_closure_diagnostics.json"
    out_plunge = args.out_plunge_csv or args.run / "absolute_plunge_records.csv"
    np.savetxt(
        out_csv,
        rows,
        delimiter=",",
        header=("t0_TNR,t1_TNR,t0_Myr,t1_Myr,captures_weighted_number,"
                "captures_raw_actual,captured_mass_Msun,Mdot_Msun_per_Myr,"
                "boundary_exits_weighted,boundary_exits_raw_actual,"
                "boundary_supplied_raw_records,boundary_supplied_mass_Msun,"
                "boundary_supplied_Mdot_Msun_per_Myr,"
                "boundary_supplied_capture_binding_current_Msun_kms2_per_Myr,"
                "boundary_advected_binding_current_Msun_kms2_per_Myr,"
                "outward_released_binding_current_Msun_kms2_per_Myr"),
        comments="",
        fmt="%.12e",
    )
    np.savetxt(
        out_plunge,
        plunge,
        delimiter=",",
        header=("isnap,x_final,x_initial,j_final,j_initial,weight_real_per_task,"
                "particle_mass_Msun,exit_time_GNC_native_Myr,rp_AU,"
                "create_time_GNC_native"),
        comments="",
        fmt="%.16e",
    )

    diag = {
        "schema": "gnc-direct-energy-current-diagnostic-v3",
        "status": status,
        "failures": diagnostic_failures,
        "acceptance_blockers": acceptance_blockers,
        "statistical_acceptance_pass": statistical_acceptance_pass,
        "physical_kernel_compatible": physical_kernel_compatible,
        "physical_blockers": physical_blockers,
        "weighted_event_table": weighted_path.name,
        "raw_event_table": raw_path.name,
        "capture_column": next(n for n in capture_names if n.lower() in
                               {x.lower() for x in weighted.dtype.names or ()}),
        "snapshots_analyzed": int(nuse),
        "total_tnr_analyzed": float(nwin * args.window_tnr),
        "total_time_myr": float(tsnap_myr[-1]),
        "tnr_myr_used_for_physical_flux": physical_tnr_myr,
        "physical_time_mapping": time_mapping,
        "gnc_native_tnr_myr": gnc_native_tnr_myr,
        "gnc_native_total_time_myr": float(tsnap_gnc_native_myr[-1]),
        "sidm_sigma_over_m_cm2_g_for_time_mapping": args.sigma_over_m_cm2_g,
        "rho_at_rh_msun_pc3": rho_h,
        "orbital_speed_at_rh_kms": vh_kms,
        "window_tnr": args.window_tnr,
        "burn_in_tnr": args.burn_in_tnr,
        "post_burn_in_time_myr": post_time_myr,
        "post_burn_in_captured_mass_event_table_msun": post_mass_table,
        "post_burn_in_captured_mass_plunge_records_msun": post_mass_records,
        "event_vs_plunge_mass_relative_difference": mass_consistency,
        "post_burn_in_direct_plunge_count": raw_capture_records,
        "post_burn_in_capture_effective_count": capture_effective_count,
        "post_burn_in_weighted_boundary_exits": float(bd_w_win[keep].sum()),
        "post_burn_in_raw_boundary_exits_actual": float(bd_n_win_actual[keep].sum()),
        "measured_mdot_msun_per_myr": mdot,
        "measured_mdot_standard_error_msun_per_myr": mdot_standard_error,
        "measured_mdot_fractional_error": mdot_fractional_error,
        "capture_counting_fractional_error": counting_fractional_error,
        "capture_window_fractional_standard_error": window_fractional_error,
        "post_burn_in_window_rate_mean_msun_per_myr": window_rate_mean,
        "post_burn_in_window_rate_std_msun_per_myr": window_rate_std,
        "last_two_capture_rate_relative_change": last_change,
        "plateau_tolerance": args.plateau_tolerance,
        "minimum_raw_captures": args.minimum_raw_captures,
        "minimum_effective_captures": args.minimum_effective_captures,
        "plateau_pass": legacy_plateau_pass,
        "raw_count_pass": raw_count_pass,
        "effective_count_pass": effective_count_pass,
        "mass_consistency_pass": mass_consistency_pass,
        "mean_capture_x_final": mean_x_final,
        "mean_capture_x_initial": mean_x_initial,
        "energy_unit_kms2": energy_unit_kms2,
        "mean_capture_binding_energy_kms2": mean_capture_binding_kms2,
        "mean_capture_binding_energy_standard_error_kms2": mean_capture_binding_standard_error_kms2,
        "mean_initial_binding_energy_kms2": mean_initial_binding_kms2,
        "capture_energy_current_msun_kms2_per_myr": l_capture,
        "capture_energy_current_standard_error_msun_kms2_per_myr": capture_energy_current_standard_error,
        "thermal_sink_current_msun_kms2_per_myr": l_sink,
        "thermal_sink_current_standard_error_msun_kms2_per_myr": sink_current_standard_error,
        "returned_energy_source_ceiling_msun_kms2_per_myr": l_source_ceiling,
        "returned_energy_source_ceiling_standard_error_msun_kms2_per_myr": capture_energy_current_standard_error,
        "mass_flux_scale_msun_per_myr": mass_flux_scale,
        "energy_flux_scale_msun_kms2_per_myr": energy_flux_scale,
        "C_M": c_m,
        "C_E_capture_inward": c_e_capture,
        "C_E_thermal_sink": c_e_sink,
        "C_E_returned_source_ceiling": c_e_source_ceiling,
        "boundary_supplied_current": {
            "selection": "post-burn-in plunge with abs(x_initial) < x_boundary",
            "interpretation": (
                "This excludes particles drawn from the initially populated "
                "inner cusp. The selected current is supplied through the outer "
                "FP reservoir rather than by transient cusp depletion."
            ),
            "raw_capture_count": int(supplied.shape[0]),
            "effective_capture_count": supplied_final_stats["effective_count"],
            "minimum_per_run_effective_capture_count": (
                args.minimum_boundary_effective_captures
            ),
            "effective_count_pass": supplied_effective_pass,
            "captured_mass_msun": supplied_mass,
            "time_myr": supplied_time,
            "mdot_msun_per_myr": supplied_mdot,
            "mdot_standard_error_msun_per_myr": supplied_mdot_standard_error,
            "mdot_fractional_error": supplied_mdot_fractional_error,
            "mean_x_initial": supplied_initial_stats["mean"],
            "mean_x_capture": supplied_final_stats["mean"],
            "mean_delta_x_released": supplied_release_stats["mean"],
            "boundary_x": x_boundary,
            "boundary_binding_energy_kms2": boundary_binding_kms2,
            "mean_capture_binding_energy_kms2": (
                supplied_final_stats["mean"] * energy_unit_kms2
            ),
            "mean_released_binding_energy_kms2": (
                supplied_release_stats["mean"] * energy_unit_kms2
            ),
            "capture_binding_current_msun_kms2_per_myr": supplied_capture_current,
            "boundary_advected_binding_current_msun_kms2_per_myr": (
                supplied_boundary_advective_current
            ),
            "outward_released_binding_current_msun_kms2_per_myr": (
                supplied_outward_current
            ),
            "outward_released_binding_current_standard_error_msun_kms2_per_myr": (
                supplied_outward_standard_error
            ),
            "current_identity": (
                "capture binding current = boundary advected binding current "
                "+ outward released binding current"
            ),
            "current_identity_residual_msun_kms2_per_myr": (
                supplied_identity_residual
            ),
            "current_identity_relative_residual": supplied_identity_relative,
            "nonnegative_released_energy_pass": supplied_energy_physical,
            "mass_current_plateau": supplied_mass_plateau,
            "released_energy_current_plateau": supplied_energy_plateau,
            "plateau_tolerance": args.plateau_tolerance,
            "plateau_pass": supplied_plateau_pass,
            "C_M": supplied_mdot / mass_flux_scale,
            "C_E_capture_binding": supplied_capture_current / energy_flux_scale,
            "C_E_boundary_advected_binding": (
                supplied_boundary_advective_current / energy_flux_scale
            ),
            "C_E_outward_released_binding": (
                supplied_outward_current / energy_flux_scale
            ),
            "sufficient_statistics": {
                "physical_mass_weight_sum_msun": supplied_final_stats["weight_sum"],
                "physical_mass_weight_squared_sum_msun2": supplied_final_stats[
                    "weight_squared_sum"
                ],
                "weighted_x_capture_sum_msun": supplied_final_stats[
                    "weighted_value_sum"
                ],
                "weighted_x_capture_squared_sum_msun": supplied_final_stats[
                    "weighted_value_squared_sum"
                ],
                "weighted_x_initial_sum_msun": supplied_initial_stats[
                    "weighted_value_sum"
                ],
                "weighted_delta_x_sum_msun": supplied_release_stats[
                    "weighted_value_sum"
                ],
                "weighted_delta_x_squared_sum_msun": supplied_release_stats[
                    "weighted_value_squared_sum"
                ],
            },
        },
        "inner_reservoir_stationarity": {
            "source": "direct per-snapshot FP mass and binding-energy inventory",
            "records_found": int(inventory.shape[0]),
            "records_expected": expected_inventory_snapshots,
            "completeness_fraction": inventory_completeness,
            "mass_inventory_plateau": inventory_mass_plateau,
            "binding_energy_inventory_plateau": inventory_energy_plateau,
            "tolerance": args.inventory_tolerance,
            "stationary": inventory_stationary,
            "event_table_post_burn_weighted_inner_number_mean": float(
                np.mean(inner_inventory_event_w[burn_snapshot:nuse])
            ),
            "event_table_post_burn_raw_inner_number_mean_actual": float(
                np.mean(inner_inventory_event_n_average[burn_snapshot:nuse]) * ranks
            ),
        },
        "run_provenance": {
            "gnc_seed": manifest.get("gnc_seed"),
            "same_initialization_seed": manifest.get("same_initialization_seed"),
            "same_evolution_seed": manifest.get("same_evolution_seed"),
            "gx_bins": manifest.get("gx_bins"),
            "dc_bins": manifest.get("dc_bins"),
            "ranks": ranks,
            "common_reservoir_fingerprint": manifest.get(
                "common_reservoir_fingerprint"
            ),
            "physical_configuration_fingerprint": manifest.get(
                "physical_configuration_fingerprint"
            ),
            "solver_executable_fingerprint": solver_fingerprint,
            "solver_executable_sha256": executable_hashes,
            "validated_patched_build": validated_build,
            "build_status": build_status,
            "cfs_sha256": manifest.get("cfs_sha256"),
            "sidm_kernel": manifest.get("sidm_kernel"),
            "sigma0_over_m_cm2_g": manifest.get("sigma0_over_m_cm2_g"),
            "yukawa_w_kms": manifest.get("yukawa_w_kms"),
        },
        "r_in_pc": rb,
        "rho_boundary_msun_pc3": rho_b,
        "sigma_boundary_kms": sigma_b,
        "mbh_msun": mbh,
        "rh_pc": rh,
        "kernel_scope": (
            "The production current uses the fixed t-channel Born/Yukawa kernel. "
            "The --kernel-compatible assertion records that the physical time "
            "mapping and the coupled halo use that same kernel."
        ),
        "window_csv": out_csv.name,
        "plunge_csv": str(out_plunge.resolve()),
        "sign_convention": (
            "Mass and binding-energy capture currents are positive inward "
            "magnitudes. Released binding energy is positive outward. Signed "
            "orbital energies at the interface and capture surface are negative."
        ),
    }
    out_json.write_text(json.dumps(diag, indent=2, sort_keys=True) + "\n")
    print(json.dumps(diag, indent=2, sort_keys=True))
    return 0 if status == "DIRECT_FP_ENERGY_DIAGNOSTIC_COMPLETE" else 4


if __name__ == "__main__":
    raise SystemExit(main())
