#!/usr/bin/env python3
"""Combine independent sparse ``(E,J)`` jump operators without selection bias."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from finite_angle_capture import aggregate_ej_jump_records
import ej_energy_ledger
import ej_energy_dg


GRID_ROUNDOFF_RTOL = 1.0e-13


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load(path: Path) -> tuple[dict[str, np.ndarray], dict]:
    with np.load(path, allow_pickle=False) as archive:
        arrays = {key: archive[key] for key in archive.files if key != "metadata_json"}
        metadata = json.loads(str(archive["metadata_json"]))
    return arrays, metadata


def edge_max_relative_difference(
    candidate: np.ndarray, reference: np.ndarray, key: str,
) -> float:
    """Compare finite grid edges while requiring identical infinite sentinels."""
    candidate = np.asarray(candidate, dtype=float)
    reference = np.asarray(reference, dtype=float)
    if candidate.shape != reference.shape:
        raise RuntimeError(f"EJ operators use different {key} shapes")
    if np.any(np.isnan(candidate)) or np.any(np.isnan(reference)):
        raise RuntimeError(f"EJ operators contain NaN in {key}")
    candidate_finite = np.isfinite(candidate)
    reference_finite = np.isfinite(reference)
    if not np.array_equal(candidate_finite, reference_finite):
        raise RuntimeError(f"EJ operators use different finite/infinite {key} edges")
    infinite = ~reference_finite
    if np.any(infinite) and not np.array_equal(
        np.signbit(candidate[infinite]), np.signbit(reference[infinite])
    ):
        raise RuntimeError(f"EJ operators use different infinite {key} sentinels")
    if not np.any(reference_finite):
        return 0.0
    candidate_values = candidate[reference_finite]
    reference_values = reference[reference_finite]
    scale = np.maximum(np.abs(reference_values), np.finfo(float).tiny)
    return float(np.max(np.abs(candidate_values - reference_values) / scale))


def combine(paths: list[Path]) -> tuple[dict[str, np.ndarray], dict]:
    if len(paths) < 2:
        raise ValueError("combine at least two independent EJ operators")
    loaded = [load(path) for path in paths]
    reference_arrays, reference = loaded[0]
    seeds = [metadata.get("seed") for _, metadata in loaded]
    if any(seed is None for seed in seeds):
        raise RuntimeError("every EJ operator needs an independent seed label")
    if len(set(seeds)) != len(seeds):
        raise RuntimeError("cannot combine repeated EJ random seeds")
    invariant_metadata = (
        "schema", "profile_sha256", "df_table_sha256", "mbh_msun", "rh_pc",
        "r_inner_pc", "reservoir_radius_pc", "sigma0_over_m_cm2_g", "w_kms",
        "energy_bins", "angular_bins", "domain_states_per_cell", "state_count",
        "post_capture_sentinel", "post_reservoir_sentinel",
        "reservoir_state_rule", "nonkepler_exchange",
    )
    optional_invariant_metadata = ("angular_grid", "energy_ledger_schema", "energy_ledger_includes_self_cell_collisions", "energy_occupation_order", "energy_dg_schema")
    edge_max_relative_spread = {
        "x_edges": 0.0,
        "j_over_jlc_edges": 0.0,
    }
    for arrays, metadata in loaded[1:]:
        for key in invariant_metadata:
            if metadata.get(key) != reference.get(key):
                raise RuntimeError(f"EJ operators differ in {key}")
        for key in optional_invariant_metadata:
            if metadata.get(key) != reference.get(key):
                raise RuntimeError(f"EJ operators differ in {key}")
        for key in ("x_edges", "j_over_jlc_edges"):
            candidate = np.asarray(arrays[key], float)
            reference_edge = np.asarray(reference_arrays[key], float)
            maximum_relative = edge_max_relative_difference(
                candidate, reference_edge, key
            )
            edge_max_relative_spread[key] = max(
                edge_max_relative_spread[key], maximum_relative
            )
            if maximum_relative > GRID_ROUNDOFF_RTOL:
                raise RuntimeError(
                    f"EJ operators use physically different {key}: maximum "
                    f"relative difference {maximum_relative:.6e} exceeds "
                    f"{GRID_ROUNDOFF_RTOL:.1e}"
                )

    statistical_weights = np.array([
        float(metadata["radial_bins"] * metadata["pairs_per_bin"])
        for _, metadata in loaded
    ])
    total_weight = float(np.sum(statistical_weights))
    energy_keys = (
        "rate_direct_binding_msun_kms2_per_myr",
        "rate_immediate_binding_msun_kms2_per_myr",
    )
    energy_presence = [
        tuple(key in arrays for key in energy_keys) for arrays, _ in loaded
    ]
    if any(any(item) and not all(item) for item in energy_presence):
        raise RuntimeError("an EJ operator contains only one capture-energy column")
    has_capture_energy = all(all(item) for item in energy_presence)
    if any(all(item) for item in energy_presence) and not has_capture_energy:
        raise RuntimeError("cannot mix legacy and capture-energy EJ operators")
    chunks = {key: [] for key in (
        "source1", "source2", "pre", "post",
        "rate_direct_msun_per_myr", "rate_immediate_msun_per_myr",
        "sample_count",
    )}
    if has_capture_energy:
        chunks.update({key: [] for key in energy_keys})
    for (arrays, _), weight in zip(loaded, statistical_weights):
        for key in ("source1", "source2", "pre", "post", "sample_count"):
            chunks[key].append(arrays[key])
        chunks["rate_direct_msun_per_myr"].append(
            arrays["rate_direct_msun_per_myr"] * weight
        )
        chunks["rate_immediate_msun_per_myr"].append(
            arrays["rate_immediate_msun_per_myr"] * weight
        )
        if has_capture_energy:
            for key in energy_keys:
                chunks[key].append(arrays[key] * weight)
    combined = aggregate_ej_jump_records(
        *(np.concatenate(chunks[key]) for key in (
            "source1", "source2", "pre", "post",
            "rate_direct_msun_per_myr", "rate_immediate_msun_per_myr",
            "sample_count",
        )),
        direct_binding_energy_rate=(
            np.concatenate(chunks[energy_keys[0]])
            if has_capture_energy else None
        ),
        immediate_binding_energy_rate=(
            np.concatenate(chunks[energy_keys[1]])
            if has_capture_energy else None
        ),
        nstates=int(reference["state_count"]),
    )
    combined["rate_direct_msun_per_myr"] /= total_weight
    combined["rate_immediate_msun_per_myr"] /= total_weight
    if has_capture_energy:
        for key in energy_keys:
            combined[key] /= total_weight
    else:
        # Preserve the legacy schema rather than presenting fabricated zero
        # capture-energy measurements.
        for key in energy_keys:
            combined.pop(key, None)
    ledger_presence = [any(k in row for k in ej_energy_ledger.KEYS) for row, _ in loaded]
    if any(ledger_presence):
        combined.update(ej_energy_ledger.combine(
            [row for row, _ in loaded], int(reference["state_count"]),
            statistical_weights / total_weight,
        ))
    dg_presence = [any(k in row for k in ej_energy_dg.KEYS) for row, _ in loaded]
    if any(dg_presence):
        if not all(all(k in row for k in ej_energy_dg.KEYS) for row, _ in loaded):
            raise RuntimeError("cannot mix complete DG1, partial, and legacy operators")
        dg_rows = [{k: row[k] * (w if k in ej_energy_dg.RATE_KEYS else 1.)
                    for k in ej_energy_dg.KEYS}
                   for (row, _), w in zip(loaded, statistical_weights / total_weight)]
        combined.update(ej_energy_dg.aggregate(dg_rows, int(reference["state_count"])))
    combined["x_edges"] = np.asarray(reference_arrays["x_edges"])
    combined["j_over_jlc_edges"] = np.asarray(
        reference_arrays["j_over_jlc_edges"]
    )
    capture = combined["post"] == int(reference["post_capture_sentinel"])
    reservoir = combined["post"] == int(reference["post_reservoir_sentinel"])
    metadata = {
        **{key: reference[key] for key in invariant_metadata},
        "schema": "finite-angle-ej-jump-operator-v1",
        "combination_schema": "independent-zero-filled-weighted-mean-v1",
        "input_count": len(paths),
        "input_paths": [path.name for path in paths],
        "input_sha256": [sha256_file(path) for path in paths],
        "seeds": seeds,
        "statistical_weights": statistical_weights.tolist(),
        "combined_radial_pair_draws": int(total_weight),
        "edge_grid_roundoff_equivalence_rtol": GRID_ROUNDOFF_RTOL,
        "input_edge_grid_max_relative_spread": edge_max_relative_spread,
        "radial_bins": int(reference["radial_bins"]),
        "pairs_per_bin": int(total_weight / reference["radial_bins"]),
        "record_count": int(combined["source1"].size),
        "represented_event_records": int(np.sum(combined["sample_count"])),
        "direct_capture_roundtrip_msun_per_myr": float(np.sum(
            combined["rate_direct_msun_per_myr"][capture]
        )),
        "immediate_capture_roundtrip_msun_per_myr": float(np.sum(
            combined["rate_immediate_msun_per_myr"][capture]
        )),
        "capture_energy_columns_present": has_capture_energy,
        "reservoir_exchange_rate_isotropic_msun_per_myr": float(np.sum(
            combined["rate_direct_msun_per_myr"][reservoir]
        )),
        "note": (
            "Rates are weighted by radial_bins times pairs_per_bin. A sparse "
            "transition absent from one realization contributes zero in that "
            "realization rather than being omitted from the mean."
        ),
    }
    if has_capture_energy:
        metadata.update({
            "direct_capture_binding_roundtrip_msun_kms2_per_myr": float(
                np.sum(combined[energy_keys[0]][capture])
            ),
            "immediate_capture_binding_roundtrip_msun_kms2_per_myr": float(
                np.sum(combined[energy_keys[1]][capture])
            ),
        })
    for key in optional_invariant_metadata:
        if key in reference:
            metadata[key] = reference[key]
    return combined, metadata


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("operators", nargs="+", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--out-json", required=True, type=Path)
    args = parser.parse_args()
    if args.out.suffix != ".npz":
        raise ValueError("combined operator output must have suffix .npz")
    arrays, metadata = combine(args.operators)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.out,
        **arrays,
        metadata_json=np.asarray(json.dumps(metadata, sort_keys=True)),
    )
    report = {
        **metadata,
        "path": args.out.name,
        "sha256": sha256_file(args.out),
    }
    args.out_json.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
