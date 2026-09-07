#!/usr/bin/env python3
"""Run the fixed-snapshot operator ladder without historical Slurm dependencies."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import subprocess
import sys

REPO=Path(__file__).resolve().parents[1]
FP=REPO/'fp_solver'

def run(script, *args, allowed=(0,)):
    command=[sys.executable,str(FP/script),*map(str,args)]
    result=subprocess.run(command,check=False)
    if result.returncode not in allowed:
        raise RuntimeError(f'{script} failed with exit code {result.returncode}')
    return result.returncode

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=['operator','steady','aggregate'])
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--plan',type=Path,default=REPO/'configs/production_operators.json')
    p.add_argument('--index',type=int,help='operator task index, 0 through 37')
    args=p.parse_args();root=args.root.resolve();plan=json.loads(args.plan.read_text())
    criteria=json.loads((root/'bridge/boundary_criteria.json').read_text())['criteria']
    if args.stage=='operator':
        if args.index is None or not 0<=args.index<len(plan['tasks']):p.error('valid --index required')
        task=plan['tasks'][args.index];out=root/'operators'/f"{task['group']}_s{task['seed']}"
        out.mkdir(parents=True,exist_ok=False)
        run('finite_angle_capture.py','--profile',root/'bridge/bridged_profile.txt',
            '--df-table',root/'normalized/f_evolved_normalized.txt',
            '--normalization-json',root/'normalized/df_normalization.json',
            '--out-json',out/'run.json','--out-csv',out/'run.csv','--out-ej-operator',out/'operator.npz',
            '--mbh',plan['mbh_msun'],'--sigma0-over-m',plan['sigma0_over_m_cm2_g'],
            '--w-kms',plan['w_kms'],'--r-outer',criteria[task['criterion']]['r_boundary_pc'],
            '--radial-bins',plan['radial_bins'],'--pairs-per-bin',plan['pairs_per_radial_bin'],
            '--cap-importance-fraction',0.5,'--direction-edge-fraction',0.3,
            '--direction-edge-power',4,'--survival-anomaly-order',16,'--survival-collision-order',48,
            '--ej-energy-bins',task['energy_bins'],'--ej-angular-grid',task['angular_grid'],'--seed',task['seed'])
        return 0
    groups=list(dict.fromkeys(t['group'] for t in plan['tasks']))
    if args.stage=='steady':
        for group in groups:
            out=root/'steady_reproduction'/group;out.mkdir(parents=True,exist_ok=False)
            seeds=[t['seed'] for t in plan['tasks'] if t['group']==group]
            operators=[root/'operators'/f'{group}_s{s}'/'operator.npz' for s in seeds]
            run('combine_ej_operators.py',*operators,'--out',out/'operator_combined.npz','--out-json',out/'operator_combined.json')
            for label,op in [*[(f's{s}',o) for s,o in zip(seeds,operators)],('combined',out/'operator_combined.npz')]:
                run('solve_ej_steady_state.py','--operator',op,'--out-json',out/f'steady_{label}.json',
                    '--out-csv',out/f'steady_{label}.csv','--residual-tolerance',1e-7,'--maximum-evaluations',10000)
        return 0
    out=root/'aggregate_reproduction';out.mkdir(parents=True,exist_ok=False)
    results=[root/'steady_reproduction'/g/f'steady_{s}.json' for g in groups for s in ('s314159','s271828','combined')]
    run('assess_ej_convergence.py',*results,'--out-json',out/'ej_convergence.json',
        '--out-csv',out/'ej_convergence.csv','--fractional-tolerance',0.1,'--minimum-seeds',2,
        '--minimum-boundaries',3,'--primary-boundary-radius',criteria['local_circular']['r_boundary_pc'],
        '--require-capture-energy',allowed=(0,2))
    # The assembler checks every mass gate independently. A failed energy test is
    # recorded and cannot authorize an absolute energy closure.
    run('assemble_fluid_mass_closure.py','--remap-json',root/'remap/remap.json',
        '--bridge-json',root/'bridge/bridge.json','--normalization-json',root/'normalized/df_normalization.json',
        '--convergence-json',out/'ej_convergence.json',
        '--steady-json',root/'steady_reproduction/local_circular_e64_j257/steady_combined.json',
        '--out',out/'fluid_mass_closure.json')
    return 0

if __name__=='__main__':
    raise SystemExit(main())
