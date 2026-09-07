#!/usr/bin/env python3
"""Manufactured steady-cusp transport benchmark, not a halo evolution.

Reference: Shapiro & Paschalidis 2014, https://arxiv.org/abs/1402.0005.
For sigma proportional v^-a the weakly collisional cusp has
rho proportional r^(-(3+a)/4), v^2=GM/[(1+beta)r]. This tests whether the
implemented BH-limited conductivity gives a radially constant positive L,
including refinement of the stock finite-difference flux. Constants are scaled
out. It does not calibrate the conduction amplitude or the Yukawa transition.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from fluid_bh_scaleheight import bh_corrected_lmfp_inverse, shell_luminosity


def run_benchmark():
    rows=[]
    for a in (0.,4.):
        beta=(3.+a)/4.
        for n in (64,128,256):
            r=np.geomspace(1e-4,.1,n);rho=1e-10*r**(-beta)
            v2=1./((beta+1.)*r);v=np.sqrt(v2)
            # Stock dimensionless LMFP k proportional rho sigma v^3.
            k_stock=rho*v**(-a)*v**3
            inverse,_=bh_corrected_lmfp_inverse(1./k_stock,r,rho,1.)
            k=1./inverse
            flux=shell_luminosity(r,1.5*v2,k)
            analytic=1.5*k/(beta+1.)
            sl=slice(3,-3)
            row={'cross_section_velocity_exponent':a,'density_slope':beta,'shells':n,
                 'maximum_relative_flux_error':float(np.max(abs(flux[sl]/analytic[sl]-1.))),
                 'relative_luminosity_spread':float(np.ptp(flux[sl])/np.mean(flux[sl])),
                 'luminosity_outward':bool(np.all(flux[sl]>0))}
            rows.append(row)
    gates={'outward_luminosity':all(x['luminosity_outward'] for x in rows),
           'constant_cusp_luminosity':max(x['relative_luminosity_spread'] for x in rows)<1e-6,
           'flux_error_decreases_with_refinement':all(rows[j+1]['maximum_relative_flux_error']<rows[j]['maximum_relative_flux_error'] for j in (0,1,3,4)),
           'finest_flux_error_below_2_percent':max(rows[j]['maximum_relative_flux_error'] for j in (2,5))<.02}
    return {'status':'MANUFACTURED_CUSP_TRANSPORT_PASS' if all(gates.values()) else 'FAIL',
            'gates':gates,'rows':rows,'source':'https://arxiv.org/abs/1402.0005',
            'scope':'Manufactured BH-dominated LMFP power laws. Not a time-evolved cusp, production Yukawa solution, or calibrated luminosity.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out-json',type=Path);args=p.parse_args()
    result=run_benchmark();s=json.dumps(result,indent=2)+'\n'
    if args.out_json:
        args.out_json.parent.mkdir(parents=True,exist_ok=True);args.out_json.write_text(s)
    print(s);raise SystemExit(0 if result['status']=='MANUFACTURED_CUSP_TRANSPORT_PASS' else 2)
