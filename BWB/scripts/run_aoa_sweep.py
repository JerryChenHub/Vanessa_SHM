"""Run both BWB parts at 3, 4, 5 and 6 degrees, preserving actual solver evidence.

Run from the repository root. This is a finite foreground batch, not a scheduled
task. The existing 2-degree input and results remain the comparison baseline.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import zipfile

from BWB.paths import REPOSITORY_ROOT as ROOT, WORKSPACE_DIR
from .build_sol144_demo import build
from .build_bwb_aero_panels import provenance
from .run_sol144 import run
from .export_sol144_strain import export_case
from .check_sol144_reasonableness import check

ANGLES = (3, 4, 5, 6)
EXPECTED_DIAGNOSTICS = {
    'Intel MKL ERROR: Parameter 13 was incorrect on entry to DGEMM .',
    'Launcher returned exit code 1',
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def archive_evidence(case, output):
    paths = sorted(p for p in case.iterdir() if p.is_file() and
                   (p.suffix.lower() in ('.bdf', '.inc', '.f06', '.f04', '.log', '.op2', '.json')
                    or p.name.startswith('launcher.')))
    archive = output / 'solver_evidence.zip'
    with zipfile.ZipFile(archive, 'w', zipfile.ZIP_DEFLATED, compresslevel=9) as stream:
        for path in paths:
            stream.write(path, path.name)
    with zipfile.ZipFile(archive) as stream:
        if stream.testzip() is not None:
            raise ValueError('Archive integrity check failed')
        for path in paths:
            if hashlib.sha256(stream.read(path.name)).hexdigest() != sha(path):
                raise ValueError(f'Archive byte mismatch: {path.name}')
    (output / 'solver_evidence_manifest.json').write_text(json.dumps({
        'archive': archive.name, 'sha256': sha(archive),
        'source_run': case.relative_to(ROOT).as_posix(),
        'members_sha256': {p.name: sha(p) for p in paths},
    }, indent=2) + '\n', encoding='utf-8')


def execute(run_root, solver, license_server):
    run_root = Path(run_root).resolve()
    if not run_root.is_relative_to(WORKSPACE_DIR):
        raise ValueError('Run directory must be in workspace/BWB')
    input_root = ROOT / 'BWB/cases/aoa_sweep'
    output_root = ROOT / 'BWB/results/aoa_sweep'
    baseline_review = ROOT / 'BWB/results/sol144_demo/runtime_diagnostic_review.json'
    jobs = [(alpha, part, f'aoa_{alpha:02d}') for alpha in ANGLES for part in (1, 2)]
    for _, part, name in jobs:
        if (run_root / name / f'part_{part}' / 'bwb_sol144.f06').exists():
            raise FileExistsError('Use a fresh run directory; previous evidence is preserved')
        if (output_root / name / f'part_{part}' / 'case_review.json').exists():
            raise FileExistsError('Sweep results already exist; do not silently replace them')
    for alpha in ANGLES:
        build(input_root / f'aoa_{alpha:02d}', alpha_deg=alpha)
    for alpha, part, name in jobs:
        source = input_root / name / f'part_{part}'
        case = run_root / name / f'part_{part}'
        output = output_root / name / f'part_{part}'
        print(f'RUNNING AOA {alpha} deg / part {part}', flush=True)
        run(source, case, solver, license_server)
        summary = json.loads((case / 'run-summary.json').read_text())
        errors = set(summary['runtime_errors'])
        if (not summary['solver_completed'] or not summary['fresh_f06'] or
                summary['fatal_messages'] or summary['warning_messages'] or
                (errors and errors != EXPECTED_DIAGNOSTICS)):
            raise RuntimeError(f'Unreviewed solver failure or warning: {case}')
        if errors and not baseline_review.is_file():
            raise RuntimeError('Missing prior empty-matrix diagnostic evidence')
        strain = export_case(case, output)
        physics = check(case, output, part)
        (output / 'reasonableness.json').write_text(json.dumps(physics, indent=2) + '\n', encoding='utf-8')
        f04 = (case / 'bwb_sol144.f04').read_text(errors='replace')
        feedback = [line.strip() for line in f04.splitlines()
                    if re.search(r'\bKALL\s+17346\s+17346', line)]
        f06 = (case / 'bwb_sol144.f06').read_text(errors='replace')
        monitors = re.findall(r'^\s+CZ\s+([\d.E+-]+)\s+([\d.E+-]+)\s*$', f06, re.MULTILINE)
        rigid, elastic = map(float, monitors[-1]) if monitors else (0., 0.)
        accepted = (physics['all_coarse_screens_pass'] and bool(feedback) and
                    abs(elastic - rigid) > 1e-4 and strain['max_abs_root_dof'] == 0.)
        review = {
            'status': 'accepted_for_workflow_demo_with_documented_runtime_diagnostic' if accepted and errors
                      else 'accepted_for_workflow_demo' if accepted else 'review_needed',
            'accepted_for_demo': accepted, 'aoa_deg': alpha, 'part': part,
            'raw_solver_status': summary['status'], 'runtime_errors_retained': sorted(errors),
            'prior_runtime_diagnostic_evidence': provenance(baseline_review) if errors else None,
            'runtime_review_scope': 'Same solver, geometry, Mach, q, fixed trim-variable count and error signature as the traced 2-degree baseline; these angle cases are not individually verbose-traced.',
            'source_f06_sha256': sha(case / 'bwb_sol144.f06'),
            'nonzero_aero_stiffness_evidence': feedback,
            'rigid_body_up_force_N': rigid, 'elastic_body_up_force_N': elastic,
            'elastic_increment_percent': (elastic / rigid - 1.) * 100 if rigid else None,
            'reasonableness_sha256': sha(output / 'reasonableness.json'),
            'strain_summary_sha256': sha(output / 'strain_summary.json'),
            'driver': provenance(Path(__file__)),
        }
        (output / 'case_review.json').write_text(json.dumps(review, indent=2) + '\n', encoding='utf-8')
        archive_evidence(case, output)
        if not accepted:
            raise RuntimeError(f'Reasonableness or coupling check failed: {case}')
        print(f'REVIEWED AOA {alpha} deg / part {part}: force={elastic:.6f} N', flush=True)
    print('All eight actual SOL144 runs exported and checked.', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-root', type=Path, default=WORKSPACE_DIR / 'aoa_sweep/run_01')
    parser.add_argument('--solver', default=r'D:\Siemens\Femap 2406\nastran\bin\nastran.exe')
    parser.add_argument('--license-server', default='29000@license-siemens.eng.uci.edu')
    args = parser.parse_args()
    execute(args.run_root, args.solver, args.license_server)
