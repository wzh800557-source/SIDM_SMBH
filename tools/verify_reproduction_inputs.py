#!/usr/bin/env python3
"""Check released inputs and run-plan coverage without launching simulations."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import re

REPO = Path(__file__).resolve().parents[1]


def verify() -> dict:
    manifest = json.loads((REPO / 'data/manifest.json').read_text())
    paths = set()
    for record in manifest['files']:
        path = REPO / record['path']
        assert path.is_relative_to(REPO)
        assert path.is_file(), record['path']
        assert hashlib.sha256(path.read_bytes()).hexdigest() == record['public_sha256'], record['path']
        paths.add(record['path'])
    actual = {str(p.relative_to(REPO)) for folder in ('production', 'reference')
              for p in (REPO / 'data' / folder).rglob('*') if p.is_file()}
    assert actual == paths, actual.symmetric_difference(paths)
    plan = json.loads((REPO / 'configs/production_operators.json').read_text())
    tasks = plan['tasks']
    assert len(tasks) == 38
    assert len({(t['group'], t['seed']) for t in tasks}) == 38
    groups = {t['group'] for t in tasks}
    assert len(groups) == 19
    for group in groups:
        assert {t['seed'] for t in tasks if t['group'] == group} == {314159, 271828}
    criteria = json.loads((REPO / 'data/production/bridge/boundary_criteria.json').read_text())['criteria']
    assert all(t['criterion'] in criteria for t in tasks)
    with (REPO / 'data/reference/scan20/scan20_results.csv').open(newline='') as stream:
        rows = list(csv.DictReader(stream))
    summary = json.loads((REPO / 'data/reference/scan20/scan20_summary.json').read_text())
    assert len(rows) == summary['case_count'] == 1200
    assert summary['grid_shape'] == [3, 20, 20]
    grid = json.loads((REPO / 'parameter_scan20/scan_grid.json').read_text())
    assert len(grid['cases']) == 1200
    # Every explicitly named Python entry point in the submission scripts must
    # be present. This is a file-coverage check, not an execution test.
    sources = {p.name for p in REPO.rglob('*.py') if '.git' not in p.parts}
    dependencies = set()
    for path in (REPO / 'slurm').rglob('*'):
        if path.suffix in {'.sh', '.sbatch'}:
            dependencies.update(re.findall(r'\b([A-Za-z_][A-Za-z_0-9]*\.py)\b', path.read_text()))
    assert not dependencies - sources, sorted(dependencies - sources)
    assert (REPO / 'README.md').stat().st_size == 0
    return {
        'status': 'INPUTS_AND_RUN_PLAN_PASS',
        'released_data_files': len(paths),
        'operator_tasks': len(tasks),
        'operator_groups': len(groups),
        'dense_scan_cases': len(rows),
        'submission_python_entry_points': len(dependencies),
        'scope': 'File integrity and run-plan coverage only. No simulations are run by this check.'
    }


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
