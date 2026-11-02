"""Solve and inspect the independent BWB components at 200 m/s and AOA 3--6.

Physical screen failures are reported without discarding complete linear results.
This finite foreground batch creates no scheduled job and preserves old cases.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re

from BWB.paths import REPOSITORY_ROOT as ROOT, WORKSPACE_DIR, TRIAL_CASES_DIR, TRIAL_RESULTS_DIR
from .build_sol144_demo import build
from .build_bwb_aero_panels import provenance
from .run_sol144 import run
from .run_aoa_sweep import archive_evidence, sha, EXPECTED_DIAGNOSTICS

ANGLES = (3, 4, 5, 6)
INPUT_ROOT = TRIAL_CASES_DIR / 'speed_200'
OUTPUT_ROOT = TRIAL_RESULTS_DIR / 'speed_200'


def solve(run_root, solver, license_server):
    jobs = [(a, p) for a in ANGLES for p in (1, 2)]
    for alpha, part in jobs:
        case = run_root / f'aoa_{alpha:02}/part_{part}'
        if (case / 'bwb_sol144.f06').exists():
            raise FileExistsError(f'Previous solver evidence must be preserved: {case}')
    for alpha in ANGLES:
        destination = build(INPUT_ROOT / f'aoa_{alpha:02}', velocity=200., alpha_deg=alpha)
        for part in (1, 2):
            folder = destination / f'part_{part}'
            deck = folder / 'bwb_sol144.bdf'
            text = deck.read_text().replace('ASSUMED LOW SPEED', 'ASSUMED 200 M/S').replace(
                'ASSUMED_LOW_SPEED_ELASTIC', 'ASSUMED_200MPS_ELASTIC')
            deck.write_text(text, encoding='ascii')
            manifest_path = folder / 'case_manifest.json'
            manifest = json.loads(manifest_path.read_text())
            manifest['inputs'] = [provenance(folder / Path(item['path']).name) for item in manifest['inputs']]
            manifest['speed_sweep_driver'] = provenance(Path(__file__))
            manifest_path.write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
    for alpha, part in jobs:
        relative = f'aoa_{alpha:02}/part_{part}'
        print(f'RUNNING 200 m/s, AOA {alpha}, part {part}', flush=True)
        run(INPUT_ROOT / relative, run_root / relative, solver, license_server)
    print('All eight solver invocations finished; inspect raw status before using results.', flush=True)


def postprocess(run_root):
    from .export_sol144_strain import export_case
    from .check_sol144_reasonableness import check

    failures = []
    diagnostic_review = OUTPUT_ROOT / 'runtime_diagnostic_review.json'
    for alpha in ANGLES:
        for part in (1, 2):
            relative = f'aoa_{alpha:02}/part_{part}'
            case, output = run_root / relative, OUTPUT_ROOT / relative
            output.mkdir(parents=True, exist_ok=True)
            summary = json.loads((case / 'run-summary.json').read_text())
            errors = set(summary['runtime_errors'])
            complete = (summary['solver_completed'] and summary['fresh_f06'] and
                        not summary['fatal_messages'] and (not errors or errors == EXPECTED_DIAGNOSTICS))
            if not complete:
                archive_evidence(case, output)
                (output / 'case_review.json').write_text(json.dumps({
                    'status': 'solver_failed_or_unreviewed_error', 'accepted_for_demo': False,
                    'aoa_deg': alpha, 'part': part, 'raw_run_summary': summary,
                }, indent=2) + '\n', encoding='utf-8')
                failures.append(relative)
                continue
            strain = export_case(case, output)
            physics = check(case, output, part)
            (output / 'reasonableness.json').write_text(json.dumps(physics, indent=2) + '\n', encoding='utf-8')
            f04 = (case / 'bwb_sol144.f04').read_text(errors='replace')
            feedback = [line.strip() for line in f04.splitlines()
                        if re.search(r'\bKALL\s+17346\s+17346', line)]
            f06 = (case / 'bwb_sol144.f06').read_text(errors='replace')
            monitors = re.findall(r'^\s+CZ\s+([\d.E+-]+)\s+([\d.E+-]+)\s*$', f06, re.MULTILINE)
            rigid, elastic = map(float, monitors[-1]) if monitors else (0., 0.)
            screens = physics['subcases'][0]['screens']
            numerical = (all(screens[key] for key in ('force_balance', 'moment_balance', 'root_fixed'))
                         and bool(feedback) and abs(elastic - rigid) > 1e-4)
            limits = [key for key, passed in screens.items() if not passed]
            reviewed = not errors or diagnostic_review.is_file()
            accepted = numerical and reviewed and physics['all_coarse_screens_pass'] and not summary['warning_messages']
            review = {
                'status': ('linear_result_exceeds_model_screens' if limits else
                           'linear_result_with_runtime_diagnostic' if errors else 'linear_result_complete'),
                'accepted_for_demo': accepted, 'result_integrity_checks_pass': numerical,
                'linear_physics_screens_pass': physics['all_coarse_screens_pass'],
                'failed_physical_screens': limits, 'aoa_deg': alpha, 'part': part, 'velocity_m_s': 200.,
                'raw_solver_status': summary['status'], 'runtime_errors_retained': sorted(errors),
                'warning_messages_retained': summary['warning_messages'],
                'runtime_diagnostic_review': provenance(diagnostic_review) if diagnostic_review.is_file() else None,
                'runtime_review_scope': 'At-speed trace documented separately; raw status and all diagnostics are retained.',
                'source_f06_sha256': sha(case / 'bwb_sol144.f06'),
                'nonzero_aero_stiffness_evidence': feedback,
                'rigid_body_up_force_N': rigid, 'elastic_body_up_force_N': elastic,
                'elastic_increment_percent': (elastic / rigid - 1.) * 100 if rigid else None,
                'reasonableness_sha256': sha(output / 'reasonableness.json'),
                'strain_summary_sha256': sha(output / 'strain_summary.json'),
                'driver': provenance(Path(__file__)),
                'limitations': ['A finite linear solution does not prove stability below divergence.',
                                'Yield or large-deformation screen failures invalidate physical predictions of this linear elastic model.',
                                'No buckling, plasticity, geometrically nonlinear equilibrium or flow separation is simulated.'],
            }
            (output / 'case_review.json').write_text(json.dumps(review, indent=2) + '\n', encoding='utf-8')
            archive_evidence(case, output)
            print(f'200 m/s AOA {alpha} part {part}: numerical={numerical}, failed physical screens={limits}', flush=True)
    if failures:
        raise RuntimeError(f'Incomplete or unreviewed solver results: {failures}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-root', type=Path, default=WORKSPACE_DIR / 'speed_200/run_01')
    parser.add_argument('--solver', default=r'D:\Siemens\Femap 2406\nastran\bin\nastran.exe')
    parser.add_argument('--license-server', default='29000@license-siemens.eng.uci.edu')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--solve-only', action='store_true')
    group.add_argument('--postprocess-only', action='store_true')
    args = parser.parse_args()
    run_root = args.run_root.resolve()
    if not run_root.is_relative_to(WORKSPACE_DIR):
        parser.error('Run directory must be in workspace/BWB')
    if not args.postprocess_only:
        solve(run_root, args.solver, args.license_server)
    if not args.solve_only:
        postprocess(run_root)


if __name__ == '__main__':
    main()
