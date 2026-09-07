#!/usr/bin/env python3
"""Compare local and orbit-integrated GNC/fluid boundary definitions."""
from __future__ import annotations
import argparse, hashlib, json, math
from pathlib import Path
import numpy as np

from sidm_born_kernel import tchannel_viscosity_over_sigma0

G=4.30091e-3
C_KMS=2.99792458e5
MSUN_G=1.98847e33
PC_CM=3.0856775814913673e18
CM2_G_TO_PC2_MSUN=MSUN_G/PC_CM**2


def sha256_file(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):
            digest.update(block)
    return digest.hexdigest()


def load(path):
    a=np.loadtxt(path,comments='#',ndmin=2)
    r,rho,sig=(np.asarray(a[:,i],float) for i in range(3))
    m=np.isfinite(r)&np.isfinite(rho)&np.isfinite(sig)&(r>0)&(rho>0)&(sig>0)
    o=np.argsort(r[m]); return r[m][o],rho[m][o],sig[m][o]

def logint(r,y,x):
    x=np.asarray(x,float)
    return np.exp(np.interp(np.log(x),np.log(r),np.log(y)))

def relative_speed(v_test, sigma_bg):
    return np.sqrt(np.asarray(v_test)**2 + 3.0*np.asarray(sigma_bg)**2)

def som_effective(som0, vrel, kernel, w_kms):
    if kernel == 'constant': return som0*np.ones_like(np.asarray(vrel),dtype=float)
    if kernel == 'yukawa-tchannel': return som0*tchannel_viscosity_over_sigma0(np.asarray(vrel)/w_kms)
    raise ValueError(f'unknown kernel {kernel}')

def large_angle_cross_section(som0,vrel,kernel,w_kms):
    """Cross section for individual deflections theta>pi/2."""
    vrel=np.asarray(vrel,float)
    if kernel=='constant': return 0.5*som0*np.ones_like(vrel)
    y=vrel/w_kms; aa=y*y
    # sigma_tot/sigma0=1/(1+a), and the fraction above pi/2 is 1/(a+2).
    return som0/((1.0+aa)*(2.0+aa))

def total_cross_section(som0,vrel,kernel,w_kms):
    """Total event cross section for the same microscopic kernel.

    The FP current is controlled by an angular-transport moment, whereas the
    number of discrete collisions is controlled by the total cross section.
    Their ratio determines whether loss-cone refilling is genuinely diffusive
    or is instead supplied by rare, non-local angular jumps.
    """
    vrel=np.asarray(vrel,float)
    if kernel=='constant': return som0*np.ones_like(vrel)
    if kernel=='yukawa-tchannel':
        aa=(vrel/w_kms)**2
        return som0/(1.0+aa)
    raise ValueError(f'unknown kernel {kernel}')

def local_N(r,rho,sig,som,mbh,a,kernel,w_kms):
    den=logint(r,rho,a); sbg=logint(r,sig,a)
    vc=math.sqrt(G*mbh/a)
    vrel=float(relative_speed(vc,sbg))
    return 2*math.pi*a*den*float(som_effective(som,vrel,kernel,w_kms))*vrel/vc

def knudsen(r,rho,sig,som,a,kernel,w_kms):
    den=logint(r,rho,a); vel=logint(r,sig,a)
    # The control uses the same angular-transport cross section as the orbital
    # criterion, evaluated at the background relative RMS speed.
    vrel=math.sqrt(6.0)*vel
    mfp=1.0/(den*float(som_effective(som,vrel,kernel,w_kms)))
    h=vel/math.sqrt(4.0*math.pi*G*den)
    return mfp/h

def Norb_all_components(r,rho,sig,som,mbh,a,e,kernel,w_kms,nquad=4096):
    """Return orbit-integrated (transport, large-angle, total-event) depths."""
    # Eccentric-anomaly quadrature.  The relative speed includes the local
    # isotropic background dispersion and the cross section is evaluated from
    # the same angular-transport kernel used to define the local criterion.
    u=(np.arange(nquad)+0.5)*(2*math.pi/nquad)
    q=1-e*np.cos(u)
    rr=a*q
    if rr.min()<r[0] or rr.max()>r[-1]: return np.nan, np.nan, np.nan
    dens=logint(r,rho,rr)
    sbg=logint(r,sig,rr)
    n=math.sqrt(G*mbh/a**3)          # (km/s)/pc
    v=math.sqrt(G*mbh/a)*np.sqrt(2/q-1)
    vrel=relative_speed(v,sbg)
    someff=som_effective(som,vrel,kernel,w_kms)
    somlarge=large_angle_cross_section(som,vrel,kernel,w_kms)
    somtotal=total_cross_section(som,vrel,kernel,w_kms)
    dt=q/n*(2*math.pi/nquad)         # pc/(km/s)
    return (float(np.sum(dens*someff*vrel*dt)),
            float(np.sum(dens*somlarge*vrel*dt)),
            float(np.sum(dens*somtotal*vrel*dt)))

def Norb_components(r,rho,sig,som,mbh,a,e,kernel,w_kms,nquad=4096):
    values=Norb_all_components(r,rho,sig,som,mbh,a,e,kernel,w_kms,nquad)
    return values[0],values[1]

def Norb(r,rho,sig,som,mbh,a,e,kernel,w_kms,nquad=4096):
    return Norb_components(r,rho,sig,som,mbh,a,e,kernel,w_kms,nquad)[0]

def loss_cone_R(mbh,a):
    """Keplerian loss-cone size R_lc=J_lc^2/J_c^2.

    We use J_lc=4GM/c, consistent with the plunge boundary used by GNC.
    The value is capped at unity when no non-plunging angular momentum range
    remains at the specified semimajor axis.
    """
    return min(1.0,16.0*G*mbh/(C_KMS*C_KMS*a))

def loss_cone_orbit_components(r,rho,sig,som,mbh,a,kernel,w_kms,nquad=4096):
    """Return (q, N_orb, N_large, R_lc, e_lc) at the loss-cone edge.

    The Cohn--Kulsrud filling parameter is q approximately N_orb/R_lc when
    N_orb is formed with the angular-transport cross section.  Evaluating the
    scattering moments on the orbit at J=J_lc identifies the radius that
    controls loss-cone refilling, rather than testing the collision kernel at
    the much larger FP/fluid handoff radius where q is already very large.
    """
    rlc=loss_cone_R(mbh,a)
    elc=math.sqrt(max(0.0,1.0-rlc))
    nv,nl=Norb_components(r,rho,sig,som,mbh,a,elc,kernel,w_kms,nquad)
    return nv/rlc,nv,nl,rlc,elc

def loss_cone_jump_components(r,rho,sig,som,mbh,a,kernel,w_kms,nquad=4096):
    """Return the loss-cone diffusion and discrete-event diagnostics.

    ``events_per_diffusion_time`` estimates the number of physical collisions
    during the time required for the transport moment to move a particle
    across the loss-cone width.  A Fokker--Planck absorbing boundary requires
    this number to be much larger than unity.  If it is small, capture is a
    finite-jump Boltzmann problem even when deflections larger than pi/2 carry
    only a small fraction of the viscosity current.
    """
    rlc=loss_cone_R(mbh,a)
    elc=math.sqrt(max(0.0,1.0-rlc))
    nv,nl,nt=Norb_all_components(
        r,rho,sig,som,mbh,a,elc,kernel,w_kms,nquad
    )
    q=nv/rlc
    events_per_diffusion_time=nt/q
    rms_event_over_loss_cone=(
        math.sqrt((nv/nt)/rlc) if nt>0.0 and rlc>0.0 else math.inf
    )
    return q,nv,nl,nt,rlc,elc,events_per_diffusion_time,rms_event_over_loss_cone

def root_log(fun,lo,hi,ngrid=500):
    grid=np.geomspace(lo,hi,ngrid)
    vals=np.array([fun(x) for x in grid])-1
    ok=np.isfinite(vals)
    roots=[]
    for i in range(len(grid)-1):
        if ok[i] and ok[i+1] and vals[i]*vals[i+1]<0:
            zl, zh = math.log(grid[i]), math.log(grid[i+1])
            fl, fh = vals[i], vals[i+1]
            for _ in range(80):
                zm = 0.5*(zl+zh)
                fm = fun(math.exp(zm))-1
                if fl*fm <= 0:
                    zh, fh = zm, fm
                else:
                    zl, fl = zm, fm
            roots.append(math.exp(0.5*(zl+zh)))
    return roots

def thermal_N(r,rho,sig,som,mbh,a,kernel,w_kms,emax=0.99,ne=48,nquad=2048):
    # Isotropic angular momenta imply p(e)=2e. Renormalize after excluding e>emax.
    x,w=np.polynomial.legendre.leggauss(ne)
    es=0.5*(x+1)*emax; ws=0.5*emax*w
    vals=np.array([Norb(r,rho,sig,som,mbh,a,float(e),kernel,w_kms,nquad) for e in es])
    return float(np.sum(ws*2*es*vals)/(emax**2))

def thermal_components(r,rho,sig,som,mbh,a,kernel,w_kms,emax=0.99,ne=48,nquad=2048):
    x,w=np.polynomial.legendre.leggauss(ne)
    es=0.5*(x+1)*emax; ws=0.5*emax*w
    vals=np.array([Norb_components(r,rho,sig,som,mbh,a,float(e),kernel,w_kms,nquad) for e in es])
    avg=np.sum((ws*2*es)[:,None]*vals,axis=0)/(emax**2)
    return float(avg[0]),float(avg[1])

def local_kernel_diagnostic(r,rho,sig,som,mbh,rb,kernel,w_kms):
    sbg=logint(r,sig,rb); vc=math.sqrt(G*mbh/rb)
    vrel=float(relative_speed(vc,sbg)); sv=float(som_effective(som,vrel,kernel,w_kms))
    sl=float(large_angle_cross_section(som,vrel,kernel,w_kms))
    return {'vrel_kms':vrel,'vrel_over_w':vrel/w_kms if kernel=='yukawa-tchannel' else None,
            'large_angle_cross_section_fraction_of_viscosity':sl/sv}

def classify_loss_cone_kernel(fraction, events, maximum_fraction=0.1, minimum_events=10.0):
    """Missing or non-finite critical-orbit diagnostics do not establish failure."""
    if fraction is None or events is None:
        return 'DIAGNOSTIC_UNAVAILABLE'
    if not np.isfinite(fraction) or not np.isfinite(events):
        return 'DIAGNOSTIC_UNAVAILABLE'
    if fraction < 0 or events < 0:
        return 'DIAGNOSTIC_UNAVAILABLE'
    return ('PASS' if fraction <= maximum_fraction and events >= minimum_events
            else 'NONLOCAL_JUMP_OPERATOR_REQUIRED')


def main():
    p=argparse.ArgumentParser()
    p.add_argument('profile',type=Path); p.add_argument('--mbh',type=float,default=4e6)
    p.add_argument('--sigma-over-m',type=float,default=100.); p.add_argument('--rh',type=float)
    p.add_argument('--kernel',choices=('constant','yukawa-tchannel'),default='constant')
    p.add_argument('--w-kms',type=float,default=80.0)
    p.add_argument('--max-large-angle-viscosity-fraction',type=float,default=0.10)
    p.add_argument('--minimum-events-per-loss-cone-diffusion-time',type=float,default=10.0)
    p.add_argument('--bridge-json',type=Path)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args(); r,rho,sig=load(a.profile); som=a.sigma_over_m*CM2_G_TO_PC2_MSUN
    bridge=json.loads(a.bridge_json.read_text()) if a.bridge_json else {}
    profile_sha256=sha256_file(a.profile)
    if bridge:
        if bridge.get('status')!='OK' or bridge.get('schema')!='absolute-closure-hydrostatic-bridge-v3':
            raise RuntimeError('hydrostatic bridge did not pass the v3 acceptance gates')
        if bridge.get('output_profile_sha256')!=profile_sha256:
            raise RuntimeError('boundary profile does not match the hydrostatic bridge output')
        for key,measured,expected in (
            ('mbh_msun',float(bridge['mbh_msun']),a.mbh),
            ('sigma_over_m_cm2_g',float(bridge['sigma_over_m_cm2_g']),a.sigma_over_m),
        ):
            if not math.isclose(measured,expected,rel_tol=1e-10,abs_tol=0.0):
                raise RuntimeError(f'boundary metadata mismatch for {key}')
        if bridge.get('collision_kernel')!=a.kernel:
            raise RuntimeError('boundary and bridge collision kernels differ')
        if a.kernel=='yukawa-tchannel' and not math.isclose(float(bridge['yukawa_w_kms']),a.w_kms,rel_tol=1e-10,abs_tol=0.0):
            raise RuntimeError('boundary and bridge Yukawa velocity scales differ')
    rh=a.rh or bridge.get('rh_pc') or G*a.mbh/sig[0]**2
    if bridge and not math.isclose(float(bridge['rh_pc']),float(rh),rel_tol=1e-10,abs_tol=0.0):
        raise RuntimeError('boundary and bridge influence radii differ')
    lo=max(r[0]*1.01,6*G*a.mbh/(2.99792458e5**2)*1.01); hi=min(rh/1.01,r[-1]/2)
    criteria={
      'local_circular': lambda x: local_N(r,rho,sig,som,a.mbh,x,a.kernel,a.w_kms),
      'orbit_e0p7': lambda x: Norb(r,rho,sig,som,a.mbh,x,0.7,a.kernel,a.w_kms),
      'orbit_e0p9': lambda x: Norb(r,rho,sig,som,a.mbh,x,0.9,a.kernel,a.w_kms),
      'orbit_thermal_e_lt_0p99': lambda x: thermal_N(r,rho,sig,som,a.mbh,x,a.kernel,a.w_kms),
    }
    if a.minimum_events_per_loss_cone_diffusion_time <= 1.0:
        raise ValueError('minimum diffusive-event count must exceed unity')
    out={'schema':'gnc-fluid-boundary-criterion-scan-v4','profile':a.profile.name,'profile_sha256':profile_sha256,'bridge_json_sha256':sha256_file(a.bridge_json) if a.bridge_json else None,'mbh_msun':a.mbh,'sigma_over_m_cm2_g':a.sigma_over_m,'collision_kernel':a.kernel,'yukawa_w_kms':a.w_kms if a.kernel=='yukawa-tchannel' else None,'rh_pc':rh,'criteria':{}}
    for name,f in criteria.items():
        roots=root_log(f,lo,hi)
        rin=roots[0] if roots else None
        row={'r_boundary_pc':rin,'x_boundary':rh/(2*rin) if rin else None,'N_local_at_boundary':local_N(r,rho,sig,som,a.mbh,rin,a.kernel,a.w_kms) if rin else None,'all_roots_pc':roots}
        if rin:
            row.update(local_kernel_diagnostic(r,rho,sig,som,a.mbh,rin,a.kernel,a.w_kms))
            if name=='orbit_thermal_e_lt_0p99':
                nv,nl=thermal_components(r,rho,sig,som,a.mbh,rin,a.kernel,a.w_kms)
            elif name.startswith('orbit_e'):
                ee=float(name.split('e',1)[1].replace('p','.'))
                nv,nl=Norb_components(r,rho,sig,som,a.mbh,rin,ee,a.kernel,a.w_kms)
            else:
                nv=row['N_local_at_boundary']; nl=nv*row['large_angle_cross_section_fraction_of_viscosity']
            row['orbit_large_angle_probability']=nl
            row['orbit_large_angle_fraction_of_viscosity_current']=nl/nv
            row['fp_angular_kernel_status']=('PASS' if nl/nv<=a.max_large_angle_viscosity_fraction else 'ANGULAR_JUMP_CORRECTION_REQUIRED')
        out['criteria'][name]=row
    knroots=root_log(lambda x:knudsen(r,rho,sig,som,x,a.kernel,a.w_kms),lo,min(r[-1]/1.01,1.0e4))
    out['controls']={
      'Knudsen_equals_1': {'all_roots_pc':knroots,'note':'conductivity transition, not a loss-cone validity boundary'},
      'influence_radius': {'r_boundary_pc':rh,'N_local_at_boundary':local_N(r,rho,sig,som,a.mbh,rh,a.kernel,a.w_kms),'note':'common geometric hand-off, but invalid wherever N_local>1'},
    }
    fractions=[v.get('orbit_large_angle_fraction_of_viscosity_current') for v in out['criteria'].values()]
    fractions=[x for x in fractions if x is not None]
    boundary_gate={
      'maximum_allowed_fraction':a.max_large_angle_viscosity_fraction,
      'maximum_measured_fraction':max(fractions) if fractions else None,
      'status':('PASS' if fractions and max(fractions)<=a.max_large_angle_viscosity_fraction else 'ANGULAR_JUMP_CORRECTION_REQUIRED'),
      'definition':'fraction of the handoff-orbit viscosity current supplied by individual deflections theta>pi/2',
      'role':'conservative diagnostic for general phase-space diffusion; it is not the loss-cone capture gate when q is already much greater than unity at the handoff',
    }
    out['fp_boundary_angular_kernel_diagnostic']=boundary_gate

    # Loss-cone capture is controlled near q=1, not at the handoff where
    # N_orb=1 but R_lc is tiny and hence q>>1.  Evaluate the actual
    # loss-cone-edge orbit and test the small-angle approximation there.
    rin_ref=out['criteria']['local_circular']['r_boundary_pc']
    rg=G*a.mbh/C_KMS**2
    qlo=max(16.0*rg*1.01,r[0]*1.01)
    qhi=min((rin_ref/1.01 if rin_ref else rh/1.01),rh/1.01,r[-1]/2.0)
    qfun=lambda x: loss_cone_orbit_components(
        r,rho,sig,som,a.mbh,x,a.kernel,a.w_kms,4096)[0]
    qroots=root_log(qfun,qlo,qhi,ngrid=240) if qhi>qlo else []
    qrows=[]
    for rc in qroots:
        (qval,nv,nl,nt,rlc,elc,events_per_diffusion_time,
         rms_event_over_loss_cone)=loss_cone_jump_components(
            r,rho,sig,som,a.mbh,rc,a.kernel,a.w_kms,16384)
        qrows.append({
          'semimajor_axis_pc':rc,
          'over_local_handoff_radius':rc/rin_ref if rin_ref else None,
          'x_energy':rh/(2.0*rc),
          'R_lc':rlc,
          'e_lc':elc,
          'N_orb_at_J_lc':nv,
          'q_estimate':qval,
          'large_angle_probability_per_orbit':nl,
          'large_angle_fraction_of_viscosity_current':nl/nv,
          'total_collision_probability_per_orbit':nt,
          'total_to_transport_cross_section_orbit_average':nt/nv,
          'events_per_loss_cone_diffusion_time':events_per_diffusion_time,
          'rms_single_event_transport_kick_over_loss_cone_width':rms_event_over_loss_cone,
        })
    critical_fraction=(qrows[0]['large_angle_fraction_of_viscosity_current']
                       if qrows else None)
    critical_events=(qrows[0]['events_per_loss_cone_diffusion_time']
                     if qrows else None)
    angular_pass=(critical_fraction is not None and
                  critical_fraction<=a.max_large_angle_viscosity_fraction)
    event_pass=(critical_events is not None and
                critical_events>=a.minimum_events_per_loss_cone_diffusion_time)
    critical_pass=angular_pass and event_pass
    capture_gate={
      'maximum_allowed_fraction':a.max_large_angle_viscosity_fraction,
      'minimum_events_per_loss_cone_diffusion_time':a.minimum_events_per_loss_cone_diffusion_time,
      'critical_points':qrows,
      'large_angle_tail_status':('UNAVAILABLE' if critical_fraction is None else ('PASS' if angular_pass else 'FAIL')),
      'many_small_events_status':('UNAVAILABLE' if critical_events is None else ('PASS' if event_pass else 'FAIL')),
      'status':classify_loss_cone_kernel(critical_fraction, critical_events, a.max_large_angle_viscosity_fraction, a.minimum_events_per_loss_cone_diffusion_time),
      'definition':'diffusive loss-cone validity at q=N_transport/R_lc=1: a small large-angle tail and many physical collisions per loss-cone diffusion time are both required',
      'handoff_q_estimate':(
          loss_cone_orbit_components(r,rho,sig,som,a.mbh,rin_ref,a.kernel,a.w_kms,4096)[0]
          if rin_ref else None),
      'note':'The theta>pi/2 tail alone is insufficient. If fewer than many collisions occur while diffusion crosses J_lc, the FP sink misses direct finite-angle transitions and a non-local Boltzmann or Monte-Carlo jump operator is required.',
    }
    out['fp_loss_cone_angular_kernel_gate']=capture_gate
    # Backward-compatible key consumed by production scripts.  It now carries
    # the physically relevant capture gate rather than the over-conservative
    # handoff-radius diagnostic.
    out['fp_angular_kernel_gate']=capture_gate
    a.out.parent.mkdir(parents=True,exist_ok=True); a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps(out,indent=2,sort_keys=True))
if __name__=='__main__': main()
