"""Exact sampled orbital-energy moments, separate from the population operator.

Specific energy is epsilon = v^2/2 - GM/r = -binding, in (km/s)^2.
This ledger includes self-cell collisions. It audits an orbit-space reservoir,
not a spatial conductive luminosity, and never authorizes a fluid heat source.
"""
from __future__ import annotations

import numpy as np

SCHEMA = "exact-orbital-energy-moments-v1"
INDEX_KEYS = ("ledger_source1", "ledger_source2", "ledger_pre", "ledger_post")
RATE_KEYS = tuple(
    f"ledger_rate_{model}_{moment}_msun_kms2_per_myr"
    for model in ("direct", "immediate") for moment in ("pre", "post")
) + ("ledger_rate_physical_delta_msun_kms2_per_myr",)
KEYS = INDEX_KEYS + RATE_KEYS


def aggregate(records, nstates):
    """Sum exact moments using the same occupation-pair keys as mass rates."""
    arrays = {k: np.asarray(records[k]).reshape(-1) for k in KEYS}
    if len({v.size for v in arrays.values()}) != 1:
        raise ValueError("inconsistent energy-ledger lengths")
    a, b, pre, post = (arrays[k].astype(np.int64) for k in INDEX_KEYS)
    if not (np.all((a >= 0) & (a < nstates)) and
            np.all((b >= 0) & (b < nstates)) and
            np.all((pre >= 0) & (pre < nstates)) and
            np.all((post >= -2) & (post < nstates))):
        raise ValueError("energy-ledger state outside grid")
    for k in RATE_KEYS:
        if not np.all(np.isfinite(arrays[k])):
            raise ValueError("nonfinite energy moment")
    # Production grids fit int64, checked explicitly rather than wrapping.
    if int(nstates)**3 * (int(nstates) + 2) > np.iinfo(np.int64).max:
        raise ValueError("energy-ledger sparse key exceeds int64")
    key = ((a * nstates + b) * nstates + pre) * (nstates + 2) + post + 2
    unique, inverse = np.unique(key, return_inverse=True)
    reduced = unique // (nstates + 2)
    result = dict(zip(INDEX_KEYS, (
        reduced // nstates // nstates,
        (reduced // nstates) % nstates,
        reduced % nstates,
        unique % (nstates + 2) - 2,
    )))
    result.update({k: np.bincount(inverse, weights=arrays[k]) for k in RATE_KEYS})
    return result


def build_records(source1, source2, pre_bindings, post_bindings, post_states,
                  captures, survivals, weight, nstates, nonkepler_exchange=True):
    """Retain both partners, including discarded/self-cell outcomes.

    For direct capture the modeled moments are survival weighted. The full
    physical collision delta is also kept, so omission of the rejected
    loss-cone branch appears as an explicit energy defect, not invented heat.
    """
    weight = np.asarray(weight, float).reshape(-1)
    if np.any(weight < 0) or not np.all(np.isfinite(weight)):
        raise ValueError("invalid event weight")
    chunks = []
    for pre, before, after, post, cap, survival in zip(
            (source1, source2), pre_bindings, post_bindings, post_states,
            captures, survivals):
        before, after = -np.asarray(before, float), -np.asarray(after, float)
        cap = np.asarray(cap, bool)
        survival = np.asarray(survival, float)
        if not np.all(np.isfinite(survival)) or np.any((survival < 0) | (survival > 1)):
            raise ValueError("survival probability outside [0,1]")
        post = np.where(cap, -1, np.asarray(post, np.int64))
        direct_weight = weight * np.where(cap, survival, 1.0)
        immediate_weight = weight.copy()
        # Regression-only cancellation of non-Kepler outcomes must remain
        # visible in the full physical ledger, even when omitted by the model.
        if not nonkepler_exchange:
            omitted = after >= 0
            direct_weight = np.where(omitted, 0., direct_weight)
            immediate_weight = np.where(omitted, 0., immediate_weight)
        row = dict(zip(INDEX_KEYS, (source1, source2, pre, post)))
        for model, w in (("direct", direct_weight), ("immediate", immediate_weight)):
            row[f"ledger_rate_{model}_pre_msun_kms2_per_myr"] = w * before
            row[f"ledger_rate_{model}_post_msun_kms2_per_myr"] = w * after
        row[RATE_KEYS[-1]] = weight * (after - before)
        chunks.append(row)
    return combine(chunks, nstates)


def combine(chunks, nstates, weights=None):
    if not chunks:
        raise ValueError("no energy records")
    presence = [all(k in row for k in KEYS) for row in chunks]
    if not all(presence):
        raise ValueError("missing energy ledger; do not mix old and new operators")
    if weights is None:
        weights = np.ones(len(chunks))
    weights = np.asarray(weights, float)
    if weights.size != len(chunks) or not np.all(np.isfinite(weights)) or np.any(weights < 0):
        raise ValueError("invalid ledger combination weights")
    records = {}
    for k in KEYS:
        records[k] = np.concatenate([
            np.asarray(row[k]) * (w if k in RATE_KEYS else 1)
            for row, w in zip(chunks, weights)
        ])
    return aggregate(records, nstates)


def audit(arrays, occupation, model, moment_tolerance=0.01):
    """Evaluate exact first moments without replacing them by bin centroids.

    The sum of state energy changes plus escaped/captured energy equals the
    modeled collision delta. In a conservative full collision operator this
    delta vanishes. Steady cell populations need not imply steady first moments.
    """
    if model not in ("direct", "immediate"):
        raise ValueError("unknown capture model")
    present = [k in arrays for k in KEYS]
    if not any(present):
        return {"status": "NOT_MEASURED", "fluid_coupling_authorized": False}
    if not all(present):
        raise ValueError("partial energy ledger")
    occupation = np.asarray(occupation, float)
    if np.any(occupation < 0) or not np.all(np.isfinite(occupation)):
        raise ValueError("invalid occupations")
    a, b, pre, post = (np.asarray(arrays[k], np.int64) for k in INDEX_KEYS)
    pair = occupation[a] * occupation[b]
    before = pair * arrays[f"ledger_rate_{model}_pre_msun_kms2_per_myr"]
    after = pair * arrays[f"ledger_rate_{model}_post_msun_kms2_per_myr"]
    physical = pair * arrays[RATE_KEYS[-1]]
    represented = post >= 0
    dot = -np.bincount(pre, weights=before, minlength=occupation.size)
    dot += np.bincount(post[represented], weights=after[represented], minlength=occupation.size)
    connected = np.arange(occupation.size) % 2 == 1
    capture = float(np.sum(after[post == -1]))
    returned = float(np.sum(after[post == -2]))
    internal = float(np.sum(dot[~connected]))
    reservoir = float(np.sum(dot[connected]))
    defect = float(np.sum(after - before))
    physical_defect = float(np.sum(physical))
    scale = max(float(np.sum(np.abs(before)) + np.sum(np.abs(after))), 1e-300)
    flux_scale = max(abs(capture), abs(returned), abs(reservoir), 1e-300)
    internal_l1 = float(np.sum(np.abs(dot[~connected])))
    identity = internal + reservoir + capture + returned - defect
    # Detailed per-cell nonstationarity, not just the easily cancelled global sum.
    gates = {
        "physical_pair_energy_conserved": abs(physical_defect) / scale < 1e-10,
        "ledger_identity": abs(identity) / scale < 1e-10,
        "modeled_collision_defect_small_vs_boundary_currents": abs(defect) / flux_scale < moment_tolerance,
        "internal_energy_moments_stationary": internal_l1 / flux_scale < moment_tolerance,
    }
    return {
        "schema": SCHEMA,
        "status": "MOMENT_AUDIT_PASS" if all(gates.values()) else "MOMENT_AUDIT_FAIL",
        "gates": gates,
        "energy_convention": "epsilon=v^2/2-GM/r; negative for bound particles",
        "rate_unit": "Msun (km/s)^2 / Myr",
        "capture_orbital_energy_rate": capture,
        "unrepresented_reservoir_return_energy_rate": returned,
        "connected_states_energy_change_rate": reservoir,
        "internal_states_energy_change_rate": internal,
        "internal_states_energy_change_l1": internal_l1,
        "modeled_collision_energy_defect": defect,
        "physical_pair_energy_defect": physical_defect,
        "omitted_history_energy_delta": physical_defect - defect,
        "identity_residual": identity,
        "identity_relative": abs(identity) / scale,
        "physical_pair_defect_relative": abs(physical_defect) / scale,
        "modeled_defect_over_boundary_scale": abs(defect) / flux_scale,
        "internal_l1_over_boundary_scale": internal_l1 / flux_scale,
        "moment_tolerance": moment_tolerance,
        "fluid_coupling_authorized": False,
        "note": "Orbit-space moments, not a measured spatial heat flux. Passing this audit still requires flux matching, work terms, and convergence before fluid coupling.",
    }
