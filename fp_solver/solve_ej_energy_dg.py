#!/usr/bin/env python3
"""Solve simultaneous mass and exact energy moments with positive DG1 occupations.

Research implementation. A passed solve does not establish resolution convergence,
correct loss-cone rescattering, or a spatial heat flux for fluid coupling.
"""
from __future__ import annotations
import argparse,csv,json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.sparse import coo_matrix, diags
import ej_energy_dg as dg
from solve_ej_steady_state import load_operator, sha256_file


def weak_coefficients(a, lo, hi, model):
    """Mass and binding-energy hat-function tests span {1, epsilon} per cell."""
    n=len(lo); pre=a['dg_pre'];post=a['dg_post'];represented=post>=0
    mass=a['dg_'+model+'_mass'];before=a['dg_'+model+'_pre_energy'];after=a['dg_'+model+'_post_energy']
    source1=[];source2=[];row=[];val=[]
    for state,mask,m,u,sign in ((pre,np.ones(pre.size,bool),mass,before,-1.),
                              (post,represented,mass,after,1.)):
        ids=state[mask];width=hi[ids]-lo[ids]
        # epsilon is negative. Phi_lo=(Bhi+epsilon)/width.
        vl=sign*(hi[ids]*m[mask]+u[mask])/width
        vh=sign*(-lo[ids]*m[mask]-u[mask])/width
        for endpoint,v in enumerate((vl,vh)):
            source1.append(a['dg_source1'][mask]);source2.append(a['dg_source2'][mask])
            row.append(2*ids+endpoint);val.append(v)
    return tuple(np.concatenate(v) for v in (source1,source2,row,val))


def solve(a,meta,model,max_nfev=1000,initial=None,continuation=False):
    n=meta['state_count'];nd=2*n
    lo,hi=dg.state_edges(a['x_edges'],meta['angular_bins'],4.30091e-3*meta['mbh_msun']/meta['rh_pc'])
    s1,s2,row,val=weak_coefficients(a,lo,hi,model)
    scale=np.bincount(row,weights=abs(val),minlength=nd)
    connected=(np.arange(nd)//2)%2==1
    floor=max(float(np.max(scale))*1e-12,1e-300)
    # Keep two equations/unknowns together. A weakly sampled slope must fail
    # the independent rank/energy checks rather than be silently suppressed.
    parents=(scale.reshape(n,2).max(axis=1)>floor)&(np.arange(n)%2==0)
    active=np.repeat(parents,2);idx=np.flatnonzero(active)
    if not len(idx): raise ValueError('no active internal moments')
    eqscale=np.maximum(scale[idx],floor)
    col=np.full(nd,-1,dtype=int);col[idx]=np.arange(len(idx))
    eq=col[row];keep=eq>=0
    sr1=s1[keep];sr2=s2[keep];eq=eq[keep];coef=val[keep]/eqscale[eq]
    dense=len(idx)<=2500
    def occupation(x):
        d=np.ones(nd);d[idx]=x;return d
    def residual(x):
        d=occupation(x)
        return np.bincount(eq,weights=coef*d[sr1]*d[sr2],minlength=len(idx))
    def jac(x):
        d=occupation(x); rr=[];cc=[];vv=[]
        for src,other in ((sr1,sr2),(sr2,sr1)):
            c=col[src];ok=c>=0
            rr.append(eq[ok]);cc.append(c[ok]);vv.append(coef[ok]*d[other[ok]])
        j=coo_matrix((np.concatenate(vv),(np.concatenate(rr),np.concatenate(cc))),shape=(len(idx),len(idx))).tocsr()
        return j.toarray() if dense else j
    x0=np.ones(len(idx)) if initial is None else np.maximum(np.asarray(initial)[idx],1e-25)
    upper=np.exp(5.)
    x=np.minimum(x0,upper*.99); stages=[]
    # Sum and difference of endpoint weak tests span the same mass and energy
    # constraints. Continuation starts with flat-in-cell D and gradually replaces
    # the slope constraint by its measured energy equation. Only lambda=1 is
    # eligible for acceptance, and the original exact audit remains unchanged.
    parent_scale=np.maximum(eqscale.reshape(-1,2).max(axis=1),floor)
    rr=np.repeat(np.arange(len(idx)),2)
    cc=np.tile(np.arange(len(idx)).reshape(-1,2),(1,2)).ravel()
    weights=np.tile([1.,1.,-1.,1.],len(idx)//2)
    transform=coo_matrix((weights,(rr,cc)),shape=(len(idx),len(idx))).tocsr()
    transform=diags(np.repeat(1./parent_scale,2))@transform@diags(eqscale)
    slope=coo_matrix((np.tile([-1.,1.],len(idx)//2),
                     (np.repeat(np.arange(1,len(idx),2),2),np.arange(len(idx)))),
                     shape=(len(idx),len(idx))).tocsr()
    if dense: transform=transform.toarray();slope=slope.toarray()
    for lam in ((0.,.1,.25,.5,.75,1.) if continuation else (1.,)):
        def stage_r(y):
            r=transform@residual(y)
            r[1::2]*=lam
            r[1::2]+=(1.-lam)*(y[1::2]-y[::2])
            return r
        def stage_j(y):
            base=transform@jac(y)
            fac=np.ones(len(idx));fac[1::2]=lam
            return fac[:,None]*base+(1.-lam)*slope if dense else diags(fac)@base+(1.-lam)*slope
        fit=least_squares(stage_r,x,jac=stage_j,bounds=(0.,upper),
                          method='trf',x_scale='jac',ftol=1e-12,xtol=1e-12,
                          gtol=1e-12,max_nfev=max_nfev)
        x=fit.x
        stages.append({'lambda':lam,'nfev':fit.nfev,'success':bool(fit.success),
                       'maximum_stage_residual':float(np.max(abs(stage_r(x))))})
    d=occupation(x); dm,du,audit=dg.moments(a,d,model,n)
    # Independent validation of both parent-cell moments, not only weighted fits.
    maxres=float(np.max(abs(residual(x))))
    j=jac(x)
    if dense:
        sing=np.linalg.svd(j,compute_uv=False)
        rank=int(np.sum(sing>sing[0]*max(j.shape)*np.finfo(float).eps)) if sing.size and sing[0]>0 else 0
        rank_pass=rank==len(idx)
    else:
        rank=None;rank_pass=False # explicitly unverified, never guessed from solver success
    gates={
        'optimizer_success':bool(fit.success),'scaled_moment_residual_le_1e_7':maxres<=1e-7,
        'full_jacobian_rank':rank_pass,'nonnegative_occupation':bool(np.all(d>=0)),
        'no_upper_occupation_bound_hit':bool(np.all(x<upper*(1.-1e-7))),
        'positive_capture':audit['capture_mdot_msun_per_myr']>0,
        'internal_mass_l1_over_capture_le_1_percent':audit['internal_mass_l1_over_capture']<=.01,
        'internal_energy_l1_over_boundary_le_1_percent':audit['internal_l1_over_boundary_scale']<=.01,
        'energy_defect_over_boundary_le_1_percent':audit['modeled_defect_over_boundary_scale']<=.01,
        'physical_pair_energy_conserved':audit['physical_pair_defect_relative']<=1e-10,
        'mass_identity':audit['mass_identity_relative']<=1e-10,
        'energy_identity':audit['energy_identity_relative']<=1e-10,
    }
    result={'status':'DG_MOMENT_SOLVE_PASS' if all(gates.values()) else 'DG_MOMENT_SOLVE_FAIL',
            'gates':gates,'model':model,'active_dofs':len(idx),'jacobian_rank':rank,
            'rank_method':'dense SVD' if dense else 'not verified for large sparse solve',
            'maximum_energy_bin_edge_ratio':float(np.max(hi/lo)),
            'maximum_scaled_residual':maxres,'nfev':sum(v['nfev'] for v in stages),'optimizer_message':fit.message,
            'continuation':continuation,'solver_stages':stages,
            'minimum_active_occupation':float(min(x)),'maximum_active_occupation':float(max(x)),
            'audit':audit,'fluid_coupling_authorized':False}
    return result,d


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('operator','out-json','out-csv'):p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--maximum-evaluations',type=int,default=1000)
    p.add_argument('--model',choices=('direct','immediate','both'),default='both')
    p.add_argument('--initial-population-csv',type=Path)
    p.add_argument('--continuation',action='store_true')
    args=p.parse_args();a,m=load_operator(args.operator)
    if m.get('energy_dg_schema')!=dg.SCHEMA or any(k not in a for k in dg.KEYS):
        raise ValueError('requires freshly sampled DG1 source-energy basis, not legacy aggregate moments')
    models={};occup={}
    for model in ('direct','immediate') if args.model=='both' else (args.model,):
        initial=None
        if args.initial_population_csv:
            rows=list(csv.DictReader(args.initial_population_csv.open()))
            key='occupation_'+model
            if len(rows)!=m['state_count'] or key not in rows[0]:
                raise ValueError('initial CSV does not match parent states')
            initial=np.repeat([float(r[key]) for r in rows],2)
        models[model],occup[model]=solve(a,m,model,args.maximum_evaluations,initial,args.continuation)
        print(json.dumps(models[model]),flush=True)
    args.out_csv.parent.mkdir(parents=True,exist_ok=True)
    with args.out_csv.open('w') as f:
        writer=csv.writer(f);writer.writerow(['state','endpoint']+['occupation_'+k for k in models])
        for i in range(2*m['state_count']):writer.writerow([i//2,i%2]+[occup[k][i] for k in models])
    result={'schema':dg.SCHEMA,'status':'DG_MOMENT_SOLVE_PASS' if all(v['status']=='DG_MOMENT_SOLVE_PASS' for v in models.values()) else 'DG_MOMENT_SOLVE_FAIL',
            'operator_sha256':sha256_file(args.operator),'occupation_csv_sha256':sha256_file(args.out_csv),
            'models':models,'fluid_coupling_authorized':False,
            'note':'Two-moment discretization test. No spatial heat flux or fluid response is inferred.'}
    args.out_json.write_text(json.dumps(result,indent=2)+'\n')
    return 0 if result['status']=='DG_MOMENT_SOLVE_PASS' else 2
if __name__=='__main__':raise SystemExit(main())
