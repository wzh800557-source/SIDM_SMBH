"""Positive piecewise-linear energy occupation and exact weak collision moments.

D_s(B) = (1-t) D_s,lo + t D_s,hi, t=(B-Blo)/(Bhi-Blo).
Here B=-epsilon is positive binding energy. Two equations per internal state
set its mass and signed orbital-energy rates to zero. No energy correction or
replacement by a cell-centre energy is applied to sampled collisions.
"""
from __future__ import annotations
import numpy as np

SCHEMA = 'finite-angle-energy-dg1-v1'
INDEX_KEYS = ('dg_source1', 'dg_source2', 'dg_pre', 'dg_post')
RATE_KEYS = tuple('dg_'+m+'_'+q for m in ('direct','immediate')
                  for q in ('mass','pre_energy','post_energy')) + ('dg_physical_delta',)
KEYS = INDEX_KEYS + RATE_KEYS


def state_edges(xedges, angular_bins, v0_squared):
    b = np.asarray(xedges, float)*v0_squared
    if np.any(~np.isfinite(b)) or np.any(b<=0) or np.any(np.diff(b)<=0):
        raise ValueError('invalid energy edges')
    return np.repeat(b[:-1], 2*angular_bins), np.repeat(b[1:], 2*angular_bins)


def basis(binding, state, lo, hi):
    t = (np.asarray(binding)-lo[state])/(hi[state]-lo[state])
    if np.any(~np.isfinite(t)) or np.any(t < -1e-9) or np.any(t > 1+1e-9):
        raise ValueError('source energy outside its parent cell')
    t = np.clip(t, 0., 1.)  # roundoff at the already validated endpoints only
    return np.stack((1.-t, t), axis=-1)


def aggregate(rows, nstates):
    if not rows:
        raise ValueError('no DG records')
    a = {k:np.concatenate([r[k] for r in rows]) for k in KEYS}
    if len({v.size for v in a.values()}) != 1:
        raise ValueError('inconsistent DG record lengths')
    for k in INDEX_KEYS:
        a[k]=np.asarray(a[k], np.int64)
    for k in RATE_KEYS:
        if np.any(~np.isfinite(a[k])):
            raise ValueError('nonfinite DG moment')
    if any(np.any((a[k]<0)|(a[k]>=2*nstates)) for k in INDEX_KEYS[:2]):
        raise ValueError('DG source outside basis')
    if np.any((a['dg_pre']<0)|(a['dg_pre']>=nstates)) or np.any((a['dg_post'] < -2)|(a['dg_post']>=nstates)):
        raise ValueError('DG parent state outside grid')
    if any(np.any(a['dg_'+m+'_mass']<0) for m in ('direct','immediate')):
        raise ValueError('negative collision weight')
    # Structured keys avoid mixed-radix int64 overflow at production resolution.
    keys=np.empty(a[INDEX_KEYS[0]].size,dtype=[(k,'<i8') for k in INDEX_KEYS])
    for k in INDEX_KEYS: keys[k]=a.pop(k)
    unique,inverse=np.unique(keys,return_inverse=True)
    result={k:unique[k].copy() for k in INDEX_KEYS}
    result.update({k:np.bincount(inverse,weights=a[k]) for k in RATE_KEYS})
    return result


def build_records(source1, source2, pre_bindings, post_bindings, post_states,
                  captures, survivals, weight, lo, hi, nonkepler_exchange=True):
    nstates=len(lo)
    weight=np.asarray(weight,float)
    if np.any(~np.isfinite(weight)) or np.any(weight<0):
        raise ValueError('invalid event weight')
    f1=basis(pre_bindings[0],source1,lo,hi)
    f2=basis(pre_bindings[1],source2,lo,hi)
    reduced=[]
    for i in (0,1):
        for j in (0,1):
            chunks=[]
            w=weight*f1[:,i]*f2[:,j]
            for pre,before,after,post,cap,surv in zip(
                    (source1,source2),pre_bindings,post_bindings,post_states,captures,survivals):
                before=-np.asarray(before,float); after=-np.asarray(after,float)
                cap=np.asarray(cap,bool);surv=np.asarray(surv,float)
                if np.any(~np.isfinite(surv)) or np.any((surv<0)|(surv>1)):
                    raise ValueError('invalid capture survival')
                row=dict(zip(INDEX_KEYS,(2*source1+i,2*source2+j,pre,np.where(cap,-1,post))))
                for model in ('direct','immediate'):
                    wm=w*np.where(cap,surv,1.) if model=='direct' else w.copy()
                    if not nonkepler_exchange: wm=np.where(after>=0,0.,wm)
                    row['dg_'+model+'_mass']=wm
                    row['dg_'+model+'_pre_energy']=wm*before
                    row['dg_'+model+'_post_energy']=wm*after
                row['dg_physical_delta']=w*(after-before)
                chunks.append(row)
            reduced.append(aggregate(chunks,nstates))
    return aggregate(reduced,nstates)


def moments(arrays, occupation, model, nstates):
    d=np.asarray(occupation,float)
    if d.shape != (2*nstates,) or np.any(~np.isfinite(d)) or np.any(d<0):
        raise ValueError('invalid positive DG occupation')
    src1,src2,pre,post=(arrays[k] for k in INDEX_KEYS)
    pair=d[src1]*d[src2]
    mass=pair*arrays['dg_'+model+'_mass']
    before=pair*arrays['dg_'+model+'_pre_energy']
    after=pair*arrays['dg_'+model+'_post_energy']
    mask=post>=0
    dm=-np.bincount(pre,weights=mass,minlength=nstates)
    du=-np.bincount(pre,weights=before,minlength=nstates)
    dm+=np.bincount(post[mask],weights=mass[mask],minlength=nstates)
    du+=np.bincount(post[mask],weights=after[mask],minlength=nstates)
    internal=np.arange(nstates)%2==0
    cap=float(np.sum(mass[post==-1])); ret=float(np.sum(mass[post==-2]))
    capu=float(np.sum(after[post==-1])); retu=float(np.sum(after[post==-2]))
    connectu=float(np.sum(du[~internal])); defect=float(np.sum(after-before))
    physical=float(np.sum(pair*arrays['dg_physical_delta']))
    flux_scale=max(abs(capu),abs(retu),abs(connectu),1e-300)
    interaction=max(float(np.sum(abs(before))+np.sum(abs(after))),1e-300)
    mass_scale=max(float(np.sum(mass)),1e-300)
    result={
        'capture_mdot_msun_per_myr':cap,'reservoir_return_msun_per_myr':ret,
        'capture_orbital_energy_rate':capu,'unrepresented_reservoir_return_energy_rate':retu,
        'connected_states_energy_change_rate':connectu,
        'internal_states_energy_change_rate':float(np.sum(du[internal])),
        'internal_l1_over_boundary_scale':float(np.sum(abs(du[internal])))/flux_scale,
        'internal_mass_l1_over_capture':float(np.sum(abs(dm[internal])))/max(cap,1e-300),
        'modeled_defect_over_boundary_scale':abs(defect)/flux_scale,
        'physical_pair_defect_relative':abs(physical)/interaction,
        'energy_identity_relative':abs(float(np.sum(du))+capu+retu-defect)/interaction,
        'mass_identity_relative':abs(float(np.sum(dm))+cap+ret)/mass_scale,
        'modeled_collision_energy_defect':defect,'physical_pair_energy_defect':physical,
        'rate_unit':'Msun (km/s)^2 / Myr','fluid_coupling_authorized':False,
    }
    return dm,du,result
