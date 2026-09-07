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


def load_plunge_records(run: Path) -> np.ndarray:
    files = sorted(
        run.rglob("plunge_records_*.txt"),
        key=lambda p: int(re.findall(r"(\d+)", p.stem)[-1]),
    )
    rows = []
    for p in files:
        data_lines = [
            line for line in p.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
        if not data_lines:
            continue
        a = np.loadtxt(data_lines, ndmin=2)
        if a.shape[1] != 9:
            raise ValueError(f"{p} has {a.shape[1]} plunge columns, expected 9")
        rows.append(a)
    return np.vstack(rows) if rows else np.empty((0, 9), float)


def main(argv: Iterable[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("run", type=Path)
    p.add_argument("--dt-tnr", type=float, default=0.1)
    p.add_argument("--window-tnr", type=float, default=1.0)
    p.add_argument("--burn-in-tnr", type=float, default=2.0)
    p.add_argument("--plateau-tolerance", type=float, default=0.25)
    p.add_argument("--minimum-raw-captures", type=int, default=20)
    p.add_argument("--minimum-effective-captures", type=float, default=50.0)
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
    tsnap_gnc_native_myr = get_column(weighted, ("Tsnap", "TsnapMyr"))
    if not np.all(np.isfinite(np.r_[cap_w, cap_n_average, boundary_w,
                                        boundary_n_average, tsnap_gnc_native_myr])):
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
    #          particle mass, exit time, periapse
    burn_snapshot = int(round(args.burn_in_tnr / args.dt_tnr))
    post_plunge = plunge[plunge[:, 0] > burn_snapshot] if plunge.size else plunge
    if post_plunge.size:
        physical_number_weight = post_plunge[:, 5] / ranks
        physical_mass_weight = physical_number_weight * post_plunge[:, 6]
        post_mass_records = float(physical_mass_weight.sum())
        # GNC stores bound orbital energies with the conventional negative sign.
        # The closure uses positive binding energy x=|E|/(GM/R_h).
        x_final_binding = np.abs(post_plunge[:, 1])
        x_initial_binding = np.abs(post_plunge[:, 2])
        mean_x_final = float(np.average(x_final_binding, weights=physical_mass_weight))
        mean_x_initial = float(np.average(x_initial_binding, weights=physical_mass_weight))
        weight_sum = float(physical_mass_weight.sum())
        weight_sq_sum = float(np.square(physical_mass_weight).sum())
        capture_effective_count = weight_sum**2 / weight_sq_sum if weight_sq_sum > 0.0 else 0.0
        x_final_variance = float(np.average(
            np.square(x_final_binding - mean_x_final), weights=physical_mass_weight
        ))
        mean_x_final_standard_error = (
            math.sqrt(x_final_variance / capture_effective_count)
            if capture_effective_count > 0.0 else math.nan
        )
        raw_capture_records = int(post_plunge.shape[0])
    else:
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
    plateau_pass = last_change <= args.plateau_tolerance
    raw_count_pass = raw_capture_records >= args.minimum_raw_captures
    effective_count_pass = capture_effective_count >= args.minimum_effective_captures
    count_pass = raw_count_pass and effective_count_pass

    energy_unit_kms2 = G_PC_KMS2_MSUN * mbh / rh
    mean_capture_binding_kms2 = mean_x_final * energy_unit_kms2
    mean_initial_binding_kms2 = mean_x_initial * energy_unit_kms2
    mass_flux_scale = 4.0 * math.pi * rb**2 * rho_b * sigma_b * KMS_TO_PCMYR
    energy_flux_scale = mass_flux_scale * sigma_b**2
    c_m = mdot / mass_flux_scale
    l_capture = mdot * mean_capture_binding_kms2
    l_sink = -1.5 * mdot * sigma_b**2
    l_source_ceiling = l_capture
    c_e_capture = l_capture / energy_flux_scale
    c_e_sink = l_sink / energy_flux_scale
    c_e_source_ceiling = l_source_ceiling / energy_flux_scale

    # The capture records carry unequal physical weights.  Their effective count,
    # rather than the raw Monte-Carlo row count, sets the counting error.  The
    # window-to-window standard error additionally captures residual temporal
    # variation.  We take the larger of the two for a conservative mass-current
    # uncertainty and propagate the capture-energy mean independently.
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

    status = "ABSOLUTE_CLOSURE_MEASURED"
    failures = []
    if not plateau_pass:
        failures.append("capture rate has not plateaued")
    if not raw_count_pass:
        failures.append("too few direct plunge records")
    if not effective_count_pass:
        failures.append("too few effective weighted captures")
    if not mass_consistency_pass:
        failures.append("event and plunge mass budgets disagree")
    if not np.isfinite(c_e_capture):
        failures.append("capture-weighted energy is unavailable")
    if failures:
        status = "MORE_RUNTIME_OR_DEBUG_REQUIRED"
    statistical_acceptance_pass = not failures
    physical_kernel_compatible = (
        args.sigma_over_m_cm2_g is None or args.kernel_compatible
    )
    physical_blockers = []
    if statistical_acceptance_pass and not physical_kernel_compatible:
        status = "STATISTICALLY_MEASURED_KERNEL_NOT_VALIDATED"
        physical_blockers.append(
            "GNC uses its forward-peaked small-angle SIDM kernel, which has not "
            "been shown to match the isotropic constant-cross-section fluid run"
        )

    rows = np.column_stack([
        t0_tnr, t1_tnr, t0_myr, t1_myr, cap_w_win, cap_n_win_actual,
        mass_win, mdot_win, bd_w_win, bd_n_win_actual,
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
                "boundary_exits_weighted,boundary_exits_raw_actual"),
        comments="",
        fmt="%.12e",
    )
    np.savetxt(
        out_plunge,
        plunge,
        delimiter=",",
        header=("isnap,x_final,x_initial,j_final,j_initial,weight_real_per_task,"
                "particle_mass_Msun,exit_time_GNC_native_Myr,rp_AU"),
        comments="",
        fmt="%.16e",
    )

    diag = {
        "schema": "gnc-absolute-closure-diagnostic-v2",
        "status": status,
        "failures": failures,
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
        "plateau_pass": plateau_pass,
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
        "r_in_pc": rb,
        "rho_boundary_msun_pc3": rho_b,
        "sigma_boundary_kms": sigma_b,
        "mbh_msun": mbh,
        "rh_pc": rh,
        "kernel_scope": (
            "This GNC build uses its forward-peaked, small-angle SIDM diffusion "
            "kernel, including the energy-dependent effective Coulomb logarithm "
            "and SIDM normalization. It is not an isotropic hard-scattering "
            "operator and therefore is not automatically compatible with the "
            "constant-cross-section fluid calculation."
        ),
        "window_csv": out_csv.name,
        "plunge_csv": str(out_plunge.resolve()),
        "sign_convention": (
            "capture current is a positive inward magnitude; the fluid sink is "
            "negative and the hypothetical fully returned source is positive"
        ),
    }
    out_json.write_text(json.dumps(diag, indent=2, sort_keys=True) + "\n")
    print(json.dumps(diag, indent=2, sort_keys=True))
    if status == "ABSOLUTE_CLOSURE_MEASURED":
        return 0
    return 5 if status == "STATISTICALLY_MEASURED_KERNEL_NOT_VALIDATED" else 4


if __name__ == "__main__":
    raise SystemExit(main())
