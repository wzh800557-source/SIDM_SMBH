#!/usr/bin/env python3
"""Project the solved even-in-v_r occupation into density and pressure.

Uses the same fixed Kepler potential and tabulated DF as the jump operator.
Angular cell edges, the loss cone, and the reservoir boundary are integrated
piecewise. A DF even in v_r cannot supply a spatial odd velocity moment (heat
or advective flux), so this reconstruction does not claim to measure one.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import numpy as np

from finite_angle_capture import (G_PC_KMS2_MSUN, C_KMS, read_df_table,
                                 ej_state_indices, sha256_file)
from normalize_gnc_df import df_forward_matrix, load_profile, log_interp
from solve_ej_steady_state import load_operator


def angular_moments(r, binding, speed, gm, v0, jlc, rb, xedges, yedges, occupation):
    jmax = r * speed
    mumax = np.sqrt(np.maximum(0., 1. - (jlc / np.maximum(jmax, 1e-300))**2))
    # Each discontinuity in D(E,J) is an exact integration boundary in mu.
    jcuts = jlc * yedges[np.isfinite(yedges)]
    ratios = jcuts[None, :] / np.maximum(jmax[:, None], 1e-300)
    cuts = np.sqrt(np.maximum(0., 1. - ratios**2))
    # The connected-state flag changes when the apocentre reaches rb.
    jboundary = np.sqrt(np.maximum(0., 2. * rb**2 * (gm / rb - binding)))
    mucross = np.sqrt(np.maximum(0., 1. - (jboundary / np.maximum(jmax, 1e-300))**2))
    cuts = np.sort(np.minimum(np.column_stack((np.zeros_like(mumax), cuts,
                                               mucross, mumax)), mumax[:, None]), axis=1)
    lo, hi = cuts[:, :-1], cuts[:, 1:]
    mid = .5 * (lo + hi)
    J = jmax[:, None] * np.sqrt(1. - mid**2)
    states = ej_state_indices(np.broadcast_to(binding[:, None], J.shape), J,
                              gm, v0, jlc, rb, xedges, yedges)
    if np.any((states < 0) & ((hi-lo) > 1e-14)):
        raise ValueError("projection sampled outside represented energy grid")
    D = occupation[np.maximum(states, 0)]
    integral0 = np.sum(D * (hi-lo), axis=1)
    integral2 = np.sum(D * (hi**3-lo**3) / 3., axis=1)
    return integral0, integral2, mumax


def project_radius(r, xgrid, gx, xedges, yedges, occupation, gm, rh, rb,
                   n0, particle_mass, order=8):
    v0 = np.sqrt(gm/rh)
    psi = rh/r
    jlc = 4.*gm/C_KMS
    # Energy states with Jc < Jlc contain no outside-loss-cone particles.
    upper = min(psi, float(xedges[-1]), float(xgrid[-1]))
    lower = max(float(xgrid[0]), float(xedges[0]))
    if upper <= lower:
        raise ValueError("no represented velocities at requested radius")
    edges = np.unique(np.r_[lower, upper, xgrid[(xgrid>lower)&(xgrid<upper)],
                            xedges[(xedges>lower)&(xedges<upper)]])
    z, w = np.polynomial.legendre.leggauss(order)
    x = (.5*(edges[1:]+edges[:-1])[:, None] + .5*np.diff(edges)[:, None]*z).ravel()
    weights = (.5*np.diff(edges)[:, None]*w).ravel()
    g = np.interp(x, xgrid, gx)
    q = psi-x
    speed = v0*np.sqrt(2.*q)
    a0, a2, outside = angular_moments(r, x*v0**2, speed, gm, v0, jlc, rb,
                                    xedges, yedges, occupation)
    weight = weights*g*np.sqrt(q)
    den = float(np.sum(weight*a0))
    factor = 2./np.sqrt(np.pi)*n0*particle_mass
    rho = factor*den
    pr = factor*float(np.sum(weight*2.*v0**2*q*a2))
    pt = factor*float(np.sum(weight*2.*v0**2*q*(a0-a2)))
    iso_rho = particle_mass*n0*float((df_forward_matrix(np.array([psi]), xgrid)@gx)[0])
    outside_rho = factor*float(np.sum(weight*outside))
    return {
        "r_pc": float(r), "rho_msun_pc3": rho,
        "pressure_radial_msun_kms2_pc3": pr,
        "pressure_tangential_total_msun_kms2_pc3": pt,
        "sigma_1d_kms": np.sqrt((pr+pt)/(3.*rho)) if rho>0 else None,
        "beta": 1.-pt/(2.*pr) if pr>0 else None,
        "rho_undepleted_full_msun_pc3": iso_rho,
        "rho_undepleted_outside_lc_msun_pc3": outside_rho,
        "rho_over_undepleted_outside_lc": rho/outside_rho if outside_rho>0 else None,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("operator", "steady-json", "steady-csv", "df-table", "normalization-json", "profile", "out-json", "out-csv"):
        p.add_argument("--"+name, type=Path, required=True)
    p.add_argument("--model", choices=("direct", "immediate"), default="immediate")
    p.add_argument("--radii", type=int, default=48)
    p.add_argument("--quadrature-order", type=int, default=8)
    args=p.parse_args()
    arrays, meta=load_operator(args.operator)
    steady=json.loads(args.steady_json.read_text())
    norm=json.loads(args.normalization_json.read_text())
    xgrid, _, table, tablemeta=read_df_table(args.df_table)
    if steady["operator_sha256"] != sha256_file(args.operator):
        raise ValueError("steady-state/operator mismatch")
    for key, path in (("df_table_sha256", args.df_table), ("profile_sha256", args.profile)):
        if meta[key] != sha256_file(path) or norm[key] != meta[key]:
            raise ValueError("projection input mismatch: "+key)
    if norm.get("status") != "PASS":
        raise ValueError("unaccepted DF normalization")
    rows=list(csv.DictReader(args.steady_csv.open()))
    if [int(row["state"]) for row in rows] != list(range(meta["state_count"])):
        raise ValueError("missing or reordered occupation states")
    occ=np.array([float(row["occupation_"+args.model]) for row in rows])
    if np.any(occ<0) or not np.all(np.isfinite(occ)):
        raise ValueError("invalid occupation")
    # The solver now binds its CSV in new output. Legacy occupations may be
    # inspected, but cannot pass production input validation here.
    csv_verified=steady.get("occupation_csv_sha256")==sha256_file(args.steady_csv)
    rprof, rhoprof, sigprof=load_profile(args.profile)
    rr=np.geomspace(meta["r_inner_pc"]*1.001, meta["reservoir_radius_pc"], args.radii)
    kwargs=dict(xgrid=xgrid,gx=np.mean(table,axis=1),xedges=arrays["x_edges"],
                yedges=arrays["j_over_jlc_edges"],occupation=occ,
                gm=G_PC_KMS2_MSUN*meta["mbh_msun"],rh=meta["rh_pc"],
                rb=meta["reservoir_radius_pc"],n0=norm["n0_pc3"],
                particle_mass=norm["particle_mass_msun"])
    projected=[project_radius(r, **kwargs, order=args.quadrature_order) for r in rr]
    refined=[project_radius(r, **kwargs, order=2*args.quadrature_order) for r in rr]
    for row in refined:
        r=row["r_pc"]
        row["rho_bridge_msun_pc3"]=float(log_interp(rprof,rhoprof,r))
        row["sigma_bridge_kms"]=float(log_interp(rprof,sigprof,r))
        row["rho_over_bridge"]=row["rho_msun_pc3"]/row["rho_bridge_msun_pc3"]
        row["pressure_trace_over_bridge"]=((row["pressure_radial_msun_kms2_pc3"]+row["pressure_tangential_total_msun_kms2_pc3"])/(3.*row["rho_bridge_msun_pc3"]*row["sigma_bridge_kms"]**2))
    errors={key:float(max(abs(a[key]/b[key]-1.) for a,b in zip(projected,refined) if b[key]>0))
            for key in ("rho_msun_pc3","pressure_radial_msun_kms2_pc3","pressure_tangential_total_msun_kms2_pc3")}
    gates={"occupation_csv_verified":csv_verified,
           "population_solver_pass":steady["models"][args.model]["status"]=="EJ_STEADY_SMOKE_PASS",
           "quadrature_refinement_le_1_percent":max(errors.values())<=.01}
    out={"status":"PROJECTION_DIAGNOSTIC_PASS" if all(gates.values()) else "PROJECTION_DIAGNOSTIC_FAIL",
         "gates":gates,"model":args.model,"maximum_relative_refinement":errors,
         "rho_over_bridge_range":[min(r["rho_over_bridge"] for r in refined),max(r["rho_over_bridge"] for r in refined)],
         "pressure_trace_over_bridge_range":[min(r["pressure_trace_over_bridge"] for r in refined),max(r["pressure_trace_over_bridge"] for r in refined)],
         "fluid_coupling_authorized":False,
         "note":"Fixed Kepler potential, even-in-radial-velocity DF. These are density and pressure moments, not a measured radial conductive or advective flux."}
    args.out_json.parent.mkdir(parents=True,exist_ok=True)
    args.out_csv.parent.mkdir(parents=True,exist_ok=True)
    with args.out_csv.open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(refined[0]));writer.writeheader();writer.writerows(refined)
    args.out_json.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
    return 0 if all(gates.values()) else 2

if __name__=="__main__":
    raise SystemExit(main())
