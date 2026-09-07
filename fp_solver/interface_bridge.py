#!/usr/bin/env python3
"""Validation gates for the hydrostatic GNC-to-fluid subgrid bridge."""

from __future__ import annotations

import math


def validate_interface_bridge(
    bridge: dict,
    *,
    mbh_msun: float,
    sigma_over_m_cm2_g: float,
    w_kms: float,
    fluid_profile_sha256: str,
    fluid_inner_shell_mass_msun: float,
    closure_radius_pc: float | None,
    radius_tolerance: float = 1.0e-3,
) -> dict:
    """Verify that the unresolved GNC domain is a controlled shell fraction."""

    if bridge.get("status") != "OK":
        raise RuntimeError(f"interface bridge was not accepted: {bridge.get('status')}")
    if bridge.get("schema") != "absolute-closure-hydrostatic-bridge-v3":
        raise RuntimeError("interface bridge metadata schema is stale")
    if bridge.get("input_profile_sha256") != fluid_profile_sha256:
        raise RuntimeError("fluid profile does not match the bridge input")
    for key, expected in (
        ("mbh_msun", mbh_msun),
        ("sigma_over_m_cm2_g", sigma_over_m_cm2_g),
        ("yukawa_w_kms", w_kms),
    ):
        measured = float(bridge[key])
        if not math.isclose(measured, expected, rel_tol=1.0e-10, abs_tol=0.0):
            raise RuntimeError(f"interface bridge metadata mismatch for {key}")
    if bridge.get("fluid_subgrid_bridge_gate") != "PASS":
        raise RuntimeError("GNC-owned mass exceeds the hydrostatic bridge gate")
    mass_fraction = float(bridge["mass_inside_r_in_fraction_of_first_shell"])
    maximum_fraction = float(bridge["maximum_allowed_subgrid_mass_fraction"])
    if not 0.0 <= mass_fraction <= maximum_fraction:
        raise RuntimeError("invalid GNC-owned first-shell mass fraction")
    volume_fraction = float(bridge["volume_inside_r_in_fraction_of_first_shell"])
    if not 0.0 <= volume_fraction < 1.0:
        raise RuntimeError("invalid GNC-owned first-shell volume fraction")
    target_mass = float(bridge["target_first_shell_mass_msun"])
    reconstructed_mass = float(bridge["reconstructed_first_shell_mass_msun"])
    profile_mass_error = abs(fluid_inner_shell_mass_msun / target_mass - 1.0)
    bridge_mass_error = abs(reconstructed_mass / target_mass - 1.0)
    pressure_error = abs(
        float(bridge["pressure_edge_inside"])
        / float(bridge["pressure_edge_outside"])
        - 1.0
    )
    if profile_mass_error > 2.0e-8:
        raise RuntimeError("fluid first-shell mass differs from the bridged shell")
    if bridge_mass_error > 2.0e-8 or pressure_error > 2.0e-8:
        raise RuntimeError("hydrostatic bridge fails mass or pressure continuity")
    radius_error = None
    if closure_radius_pc is not None:
        bridge_radius = float(bridge["r_in_pc"])
        radius_error = abs(closure_radius_pc / bridge_radius - 1.0)
        if radius_error > radius_tolerance:
            raise RuntimeError("capture current and fluid bridge use different boundaries")
    return {
        "status": "PASS",
        "method": (
            "mass-matched hydrostatic continuation of the first Lagrangian "
            "shell to the orbit-defined GNC boundary"
        ),
        "r_in_pc": float(bridge["r_in_pc"]),
        "closure_to_bridge_radius_relative_error": radius_error,
        "GNC_owned_mass_msun": float(bridge["mass_inside_r_in_msun"]),
        "GNC_owned_mass_fraction_of_first_shell": mass_fraction,
        "GNC_owned_volume_fraction_of_first_shell": volume_fraction,
        "fluid_first_shell_mass_relative_error": profile_mass_error,
        "bridge_mass_relative_error": bridge_mass_error,
        "pressure_continuity_relative_error": pressure_error,
    }
