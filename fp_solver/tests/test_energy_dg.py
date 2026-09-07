"""Regression tests for sampled first moments and positive energy-basis solves."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
import ej_energy_dg as dg
import ej_energy_ledger as old
from solve_ej_energy_dg import solve,weak_coefficients


def manufactured():
    # Open, algebraic collision network, not a physical halo benchmark. Supply
    # and capture act at two distinct energies and determine both endpoint D's.
    lo=np.array([1.,1.]);hi=np.array([3.,3.]);target=np.array([.4,1.2,1.,1.])
    b=np.array([1.2,1.2,2.7,2.7]);other=np.array([2.,2.,1.5,1.5])
    s1=np.array([0,1,0,1]);s2=np.ones(4,dtype=int)
    p1=np.array([-1,0,-1,0]);p2=np.ones(4,dtype=int)
    cap=np.array([True,False,True,False]);z=np.zeros(4,bool);surv=np.ones(4)
    injected=(1.-(b-1)/2)*.4+(b-1)/2*1.2
    w=np.where(cap,1.,injected)
    args=(s1,s2,(b,other),(b,other),(p1,p2),(cap,z),(surv,surv),w)
    a=dg.build_records(*args,lo,hi)
    a['x_edges']=np.array([1.,3.]);a['j_over_jlc_edges']=np.array([1.,np.inf])
    meta={'state_count':2,'angular_bins':1,'mbh_msun':1/4.30091e-3,'rh_pc':1.}
    return a,meta,args,lo,hi,target


def main():
    a,m,args,lo,hi,target=manufactured()
    # Constant endpoint occupations reduce exactly to the old cellwise model.
    ledger=old.build_records(*args,nstates=2)
    rng=np.random.default_rng(3109)
    for model in ('direct','immediate'):
        for _ in range(10):
            occ=rng.uniform(.05,2.,2)
            dm,du,result=dg.moments(a,np.repeat(occ,2),model,2)
            ref=old.audit(ledger,occ,model)
            for key in ('capture_orbital_energy_rate','connected_states_energy_change_rate','internal_states_energy_change_rate'):
                assert np.isclose(result[key],ref[key],rtol=1e-12,atol=1e-12),(key,result,ref)
        dm,du,r=dg.moments(a,target,model,2)
        assert max(abs(dm[0]),abs(du[0]))<1e-14
        out,d=solve(a,m,model,max_nfev=100)
        assert out['status']=='DG_MOMENT_SOLVE_PASS',out
        assert np.allclose(d,target,rtol=1e-9,atol=1e-9),(d,target)
    # One constant D balances the particle population but cannot balance both
    # distinct supplied/captured energies in this manufactured example.
    cap_old=(args[-1][1]+args[-1][3])/2
    _,_,r=dg.moments(a,np.array([cap_old,cap_old,1.,1.]),'immediate',2)
    assert r['internal_mass_l1_over_capture']<1e-12
    assert r['internal_l1_over_boundary_scale']>.01
    # Both particles can stay in their parent cells yet exchange energy. These
    # outcomes must not be discarded merely because the population is unchanged.
    src=np.array([0]);res=np.array([1]);w=np.array([2.]);z=np.array([False]);one=np.ones(1)
    ar=dg.build_records(src,res,(np.array([1.2]),np.array([2.7])),
       (np.array([1.8]),np.array([2.1])),(src,res),(z,z),(one,one),w,lo,hi)
    occ=np.array([.2,.7,1.,1.]);dm,du,r=dg.moments(ar,occ,'immediate',2)
    assert np.max(abs(dm))<1e-14 and abs(du[0])>.1
    assert r['physical_pair_defect_relative']<1e-14
    assert abs(sum(du))<1e-14
    s1,s2,row,v=weak_coefficients(ar,lo,hi,'immediate')
    weak=np.bincount(row,weights=v*occ[s1]*occ[s2],minlength=4).reshape(2,2)
    assert np.allclose(weak.sum(axis=1),dm,atol=1e-14)
    assert np.allclose(-lo*weak[:,0]-hi*weak[:,1],du,atol=1e-14)
    # Structured aggregation remains correct beyond the old integer radix limit.
    n=100000;big={k:a[k].copy() for k in dg.KEYS}
    big['dg_source1']+=2*(n-2);big['dg_source2']+=2*(n-2)
    big['dg_pre']+=n-2;mask=big['dg_post']>=0;big['dg_post'][mask]+=n-2
    again=dg.aggregate([big],n)
    assert all(np.allclose(big[k],again[k]) for k in dg.KEYS)
    try:dg.basis(np.array([5.]),np.array([0]),lo,hi)
    except ValueError:pass
    else:raise AssertionError('out-of-cell energy accepted')
    print('PASS: DG1 manufactured solve, old-model recovery, nonstationary mass-only counterexample, same-cell exchange, exact weak moment identities, large-state indexing, support rejection')
if __name__=='__main__':main()
