"""Synthetic regression tests. These are not production halo results."""
import sys, math, json, csv, tempfile, hashlib
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from publication_inputs import load_validated_response, measured_response_endpoint
from boundary_energy_budget import transfer_capture
from velocity_average import maxwellian_rate_average
from sidm_born_kernel import tchannel_viscosity_over_sigma0
from boundary_criterion_scan import classify_loss_cone_kernel
from fluid_bh_scaleheight import bh_lmfp_scaleheight_factor, bh_lmfp_scaleheight_factor_physical
from solve_ej_steady_state import population_rate

checks=[]
def check(name,condition):
    assert condition,name
    checks.append(name)
def rejects(fn):
    try: fn()
    except (RuntimeError,ValueError): return True
    return False

with tempfile.TemporaryDirectory() as td:
 d=Path(td)
 (d/'fluid_response_analysis_old.json').write_text('{"status":"PASS"}')
 check('no historical fallback',rejects(lambda:load_validated_response(d)))
 ap=d/'fluid_response_analysis.json';cp=d/'fluid_response_comparison.csv'
 fields=['tau_relax','sink_over_control_rho_inner_mean_msun_pc3','sink_over_control_r_inner_pc','sink_over_control_sigma_inner_1d_kms']
 with cp.open('w') as f:
  w=csv.writer(f);w.writerow(fields);w.writerow([0,1,1,1]);w.writerow([1,1.01,.99,1.001])
 a={'status':'PASS','gates':{'control_completed_physically':True,'summary_trajectory_endpoints_match':True},'common_tau_relax':1,'inner_mean_density_sink_over_control_at_common_end':1.01,'inner_radius_sink_over_control_at_common_end':.99,'inner_dispersion_sink_over_control_at_common_end':1.001,'identity':{'fluid_comparison_csv_sha256':hashlib.sha256(cp.read_bytes()).hexdigest()},'common_time_myr':2,'captured_mass_at_common_end_msun':3}
 ap.write_text(json.dumps(a));aa,rr=load_validated_response(d)
 check('synthetic consistent ledger loads',len(rr)==2)
 check('measured endpoint used',measured_response_endpoint(aa)==(2,3))
 for label,change in [('failed status',{'status':'FAIL'}),('failed gate',{'gates':{'a':False}}),('endpoint mismatch',{'common_tau_relax':2}),('fixture rejection',{'is_fixture':True})]:
  bad=dict(a);bad.update(change);ap.write_text(json.dumps(bad));check(label,rejects(lambda:load_validated_response(d)))
 check('no endpoint from nominal stopping limit',rejects(lambda:measured_response_endpoint({})))
 ap.write_text(json.dumps(a));cp.write_text(cp.read_text()+'2,1,1,1\n');check('altered table rejected',rejects(lambda:load_validated_response(d)))

base=transfer_capture(100,300,1000,2,3)
check('mean energy capture leaves specific energy unchanged',abs(base['specific_energy']-3)<1e-14)
check('mass conservation',abs(base['mass_residual'])<1e-12)
check('thermal budget conservation',abs(base['thermal_budget_residual'])<1e-12)
check('preferentially hot capture cools',transfer_capture(100,300,1000,2,6)['specific_energy']<3)
check('preferentially cold capture raises specific energy',transfer_capture(100,300,1000,2,1)['specific_energy']>3)
check('outward heat and work balanced',abs(transfer_capture(100,300,1000,2,3,10,4)['thermal_budget_residual'])<1e-12)
check('reject exhausted shell',rejects(lambda:transfer_capture(100,300,1000,100,3)))
# Physical conversion must agree exactly with the existing dimensionless implementation.
r=np.array([.1,1.,10.]);rho=np.array([100.,3.,.2]);rs=4.;rhos=2.;mbh=17.
physical=bh_lmfp_scaleheight_factor_physical(r,rho,mbh,.03)
code=bh_lmfp_scaleheight_factor(r/rs,rho/rhos,mbh/(4*np.pi*rhos*rs**3),.03/rs)
check('4pi physical/code unit equivalence',np.allclose(physical,code,rtol=1e-14,atol=0))
check('absent critical root is unavailable',classify_loss_cone_kernel(None,None)=='DIAGNOSTIC_UNAVAILABLE')
check('nonfinite critical diagnostic is unavailable',classify_loss_cone_kernel(float('nan'),2)=='DIAGNOSTIC_UNAVAILABLE')
check('evaluated local pass',classify_loss_cone_kernel(.05,20)=='PASS')
check('evaluated local fail',classify_loss_cone_kernel(.05,.2)=='NONLOCAL_JUMP_OPERATOR_REQUIRED')
# Independent analytic noncentral Maxwell first moment, constant cross section.
for v,sig in ((0.,2.),(.1,2.),(2.,2.),(20.,2.)):
 exact=2*sig*math.sqrt(2/math.pi) if v==0 else sig*math.sqrt(2/math.pi)*math.exp(-v*v/(2*sig*sig))+(v+sig*sig/v)*math.erf(v/(math.sqrt(2)*sig))
 val=maxwellian_rate_average(v,sig,lambda x:np.ones_like(x))
 check(f'constant-cross-section analytic mean v={v}',abs(val/exact-1)<1e-10)
check('cold limit',maxwellian_rate_average(2,0,lambda x:np.ones_like(x))==2)
fun=lambda v:tchannel_viscosity_over_sigma0(v/80.)
for v,sig in ((0.,80.),(115.25,79.611),(400.,30.)):
 a=maxwellian_rate_average(v,sig,fun,192);b=maxwellian_rate_average(v,sig,fun,384)
 check(f'nonlinear quadrature convergence v={v}',abs(a/b-1)<1e-7)
# A synthetic reversible operator has equilibrium and conserves represented mass.
occ=np.ones(4);src1=np.array([0,2,0,2]);src2=np.array([1,1,3,3]);pre=np.array([0,2,0,2]);post=np.array([2,0,2,0]);rate=np.array([2.,2.,1.,1.])
dot,_=population_rate(occ,src1,src2,pre,post,rate)
check('synthetic no-sink detailed-balance equilibrium',np.max(np.abs(dot))<1e-14)
for occupation in (np.array([.2,1,.7,1]),np.array([2.,1,.1,1])):
 dot,_=population_rate(occupation,src1,src2,pre,post,rate)
 check('closed operator mass conservation',abs(dot.sum())<1e-14)
print(json.dumps({'status':'PASS','test_count':len(checks),'checks':checks,'scope':'synthetic regression and analytic limits, not a production halo or energy-closure validation'},indent=2))
