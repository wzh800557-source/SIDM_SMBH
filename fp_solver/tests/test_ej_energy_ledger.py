#!/usr/bin/env python3
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import ej_energy_ledger as el
from finite_angle_capture import scatter_equal_mass, ej_state_indices
from reproject_depleted_df import project_radius

def main():
    # Exact two-particle example: total signed energy -10 before and after.
    # Partner 2 remains in one cell but receives 4 energy units.
    row=el.build_records(np.array([0]),np.array([2]),
        (np.array([4.]),np.array([6.])),(np.array([8.]),np.array([2.])),
        (np.array([-1]),np.array([2])),(np.array([True]),np.array([False])),
        (np.array([.25]),np.array([0.])),np.array([1.]),4)
    imm=el.audit(row,np.ones(4),'immediate');direct=el.audit(row,np.ones(4),'direct')
    assert abs(imm['modeled_collision_energy_defect'])<1e-14
    assert abs(imm['capture_orbital_energy_rate']+8)<1e-14
    assert abs(direct['modeled_collision_energy_defect']-3)<1e-14
    assert abs(direct['omitted_history_energy_delta']+3)<1e-14
    assert not direct['gates']['modeled_collision_defect_small_vs_boundary_currents']
    assert not imm['fluid_coupling_authorized']
    assert el.audit({},np.ones(4),'direct')['status']=='NOT_MEASURED'
    assert imm['identity_relative']<1e-14
    occ=np.array([.7,1.,.2,1.]); weighted=el.audit(row,occ,'immediate')
    assert abs(weighted['capture_orbital_energy_rate']+8*.7*.2)<1e-14
    assert weighted['physical_pair_defect_relative']<1e-14
    merged=el.combine([row,row],4,weights=[.25,.75])
    for key in el.KEYS: np.testing.assert_allclose(merged[key],row[key])
    broken=dict(row);del broken[el.RATE_KEYS[0]]
    try: el.audit(broken,occ,'direct');raise AssertionError('accepted partial ledger')
    except ValueError: pass
    # Independent physical collision benchmark: exact momentum and kinetic
    # energy, Maxwellian pair detailed balance before binning, arbitrary D_aD_b.
    rng=np.random.default_rng(77321);n=20000;sigma=3.
    v1=rng.normal(0,sigma,(n,3));v2=rng.normal(0,sigma,(n,3))
    a,b,_=scatter_equal_mass(v1,v2,2.,rng)
    np.testing.assert_allclose(a+b,v1+v2,rtol=1e-12,atol=1e-12)
    before=np.sum(v1*v1+v2*v2,axis=1);after=np.sum(a*a+b*b,axis=1)
    assert np.max(abs(after-before))/np.max(before)<1e-14
    assert np.max(abs(np.expm1((before-after)/(2*sigma**2))))<1e-12
    psi=400.;binding=(psi-.5*np.sum(v1*v1,axis=1),psi-.5*np.sum(v2*v2,axis=1))
    binding_after=(psi-.5*np.sum(a*a,axis=1),psi-.5*np.sum(b*b,axis=1))
    # Arbitrary source labels test cancellation under depletion, including
    # pairs in the same cell, across cells, and fixed reservoir cells.
    s1=rng.integers(0,8,n);s2=rng.integers(0,8,n)
    er=el.build_records(s1,s2,binding,binding_after,
        (rng.integers(0,8,n),rng.integers(0,8,n)),
        (np.zeros(n,bool),np.zeros(n,bool)),(np.zeros(n),np.zeros(n)),
        rng.uniform(.1,2,n),8)
    audit=el.audit(er,rng.uniform(.1,2,8),'immediate')
    assert audit['physical_pair_defect_relative']<1e-13
    assert abs(audit['modeled_collision_energy_defect'])<1e-8
    # Project constant isotropic g in a Kepler potential. Without the tiny
    # loss cone rho proportional psi^(3/2), sigma_1d^2=2 psi v0^2/5.
    x=np.geomspace(1e-9,1e4,400);xe=np.geomspace(1e-9,1e4,24)
    ye=np.array([1.,3.,100.,np.inf]);D=np.ones(2*(len(xe)-1)*(len(ye)-1))
    for r in (.1,1.,10.):
        p=project_radius(r,x,np.ones_like(x),xe,ye,D,1.,1.,2.,1.,1.,order=24)
        analytic=4/(3*np.sqrt(np.pi))*r**-1.5
        assert abs(p['rho_msun_pc3']/analytic-1)<2e-5
        assert abs(p['sigma_1d_kms']**2/(.4/r)-1)<2e-5
        assert abs(p['beta'])<1e-6
        half=project_radius(r,x,np.ones_like(x),xe,ye,D*.5,1.,1.,2.,1.,1.,order=24)
        assert abs(half['rho_msun_pc3']/p['rho_msun_pc3']-.5)<1e-12
        assert abs(half['sigma_1d_kms']/p['sigma_1d_kms']-1)<1e-12
    print('PASS: paired energy conservation, self-cell exchange, rejected-history defect, depletion weighting, Maxwellian pair balance, analytic density/pressure projection')
    return 0
if __name__=='__main__': raise SystemExit(main())
