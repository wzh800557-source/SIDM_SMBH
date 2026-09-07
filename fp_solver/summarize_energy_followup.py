#!/usr/bin/env python3
"""Summarize the four energy-ledger runs without promoting diagnostics to closure."""
import argparse
import json
from pathlib import Path

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);args=p.parse_args()
    rows=[];missing=[]
    for E in (64,96):
        for seed in (314159,271828):
            tag=f'e{E}_j257_s{seed}';f=args.root/tag/'diagnostic_summary.json'
            if not f.exists():missing.append(tag);continue
            x=json.loads(f.read_text());s=json.loads((f.parent/'steady.json').read_text())
            for m in ('direct','immediate'):
                row={'energy_bins':E,'seed':seed,'model':m,
                     'population_status':s['models'][m]['status'],
                     'mdot_msun_per_myr':s['models'][m]['steady_capture_mdot_msun_per_myr'],
                     'energy_audit':x['energy_moments'][m],
                     'projection':x['projections'][m]}
                rows.append(row)
    out={'status':'FOLLOWUP_DIAGNOSTICS_COMPLETE' if not missing else 'INCOMPLETE_OR_FAILED',
         'missing_runs':missing,'rows':rows,'fluid_coupling_authorized':False,
         'note':'Compare seeds and energy grids before interpreting any current. A moment audit is not a spatial conductive-flux measurement.'}
    (args.root/'followup_summary.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'missing_runs':missing,'rows':len(rows)},indent=2))
    return 0 if not missing else 2
if __name__=='__main__':raise SystemExit(main())
