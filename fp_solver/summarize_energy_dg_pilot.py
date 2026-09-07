#!/usr/bin/env python3
"""Gather pilot diagnostics without promoting process completion to acceptance."""
import argparse,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args()
rows=[];missing=[]
for e in (32,64):
 for seed in (314159,271828):
  path=a.root/f'e{e}_fine_s{seed}'/'diagnostic_summary.json'
  if not path.exists():missing.append(str(path));continue
  v=json.loads(path.read_text())
  for m in ('direct','immediate'):
   old=v['population_control'][m];new=v['dg_models'][m]
   rows.append({'energy_bins':e,'seed':seed,'model':m,
    'population_status':old['status'],'population_mdot':old['steady_capture_mdot_msun_per_myr'],
    'population_energy_residual':old['energy_moment_audit']['internal_l1_over_boundary_scale'],
    'dg':new})
out={'status':'DG_PILOT_DIAGNOSTICS_COMPLETE' if not missing else 'INCOMPLETE_OR_FAILED',
 'missing_runs':missing,'rows':rows,'fluid_coupling_authorized':False,
 'note':'Numerical moment solves still require seed, energy and angular refinement, self-consistent loss-cone transport, spatial flux and work-term matching before fluid evolution.'}
(a.root/'dg_pilot_summary.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
