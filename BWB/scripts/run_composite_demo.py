"""Build, solve, export and check one fresh two-part woven-composite trial.

Run from the repository root, for example::

    python -m BWB.scripts.run_composite_demo --output workspace/BWB/trials/composite_skin/new_baseline
    python -m BWB.scripts.run_composite_demo --output workspace/BWB/trials/composite_skin/new_angle_30 --angle-deg 30

Each part is launched once through run_sol144.run with its normal solver and
license settings. Complete results may be exported for diagnosis even when
runtime errors are recorded; those errors are never accepted automatically.
Exit 0 means clean execution plus passing coarse screens, 2 means review is
required, and an exception reports an interrupted/failed workflow. No retry,
scheduled job, Femap session or permission workaround is created.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path

from BWB.paths import REPOSITORY_ROOT as ROOT, TRIALS_DIR

DEFAULT_CONFIG = ROOT / 'BWB/config/composite_skin.json'
DEFAULT_SOLVER = r'D:\Siemens\Femap 2406\nastran\bin\nastran.exe'
DEFAULT_LICENSE_SERVER = '29000@license-siemens.eng.uci.edu'


def timestamp():
    return datetime.now(timezone.utc).isoformat()


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run_trial(config, output, angle_override=None, solver=DEFAULT_SOLVER,
              license_server=DEFAULT_LICENSE_SERVER, threads=2):
    """Create a fresh trial; return (aggregate summary, process exit code)."""
    output = Path(output).resolve()
    if not output.is_relative_to(TRIALS_DIR.resolve()) or output == TRIALS_DIR.resolve():
        raise ValueError('Choose a new trial directory under workspace/BWB/trials')
    if output.exists():
        raise FileExistsError(f'Existing trial must be preserved; choose another --output: {output}')
    config, solver = Path(config).resolve(strict=True), Path(solver).resolve(strict=True)
    if not config.is_file() or not solver.is_file():
        raise ValueError('Config and solver must both be existing files')
    if angle_override is not None and not math.isfinite(angle_override):
        raise ValueError('--angle-deg must be finite')
    if threads not in range(1, 17) or not str(license_server).strip():
        raise ValueError('Use 1..16 threads and a nonempty license server')

    # Delay imports so --help remains usable while optional workflow modules
    # are unavailable. Resolve all stages before creating files or launching.
    from .build_composite_skin import build
    from .run_sol144 import run
    from .export_composite_sol144 import export_case
    from .check_composite_sol144 import check

    summary = {
        'schema_version': 1, 'workflow_status': 'building', 'started_utc': timestamp(),
        'trial_root': output.relative_to(ROOT).as_posix(),
        'config': str(config), 'config_sha256': file_hash(config),
        'angle_override_deg': angle_override,
        'solver': str(solver), 'license_server': license_server, 'threads': threads,
        'driver': 'BWB/scripts/run_composite_demo.py', 'driver_sha256': file_hash(__file__),
        'runtime_diagnostics_autoaccepted': False,
        'diagnostic_policy': ('Raw solver status, launcher return code, errors and warnings are preserved. '
                              'Exporting a complete field for inspection is not acceptance of its diagnostics. '
                              'Coarse-screen results are reported independently of runtime cleanliness.'),
        'limitations': ['Passing coarse screens is not material-strength, mesh, experimental or aeroelastic-stability validation.',
                        'This workflow does not assign or certify composite strength allowables or failure indices.'],
        'cases': [],
    }
    output.mkdir(parents=True, exist_ok=False)
    summary_path = output / 'trial_summary.json'

    def save():
        summary['updated_utc'] = timestamp()
        summary_path.write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n', encoding='utf-8')

    save()
    try:
        cases_root = Path(build(config, output, angle_override)).resolve()
        if cases_root != output / 'cases':
            raise RuntimeError('Builder returned an unexpected cases directory')
        summary['workflow_status'] = 'running'
        save()
        for part in (1, 2):
            source = cases_root / f'part_{part}'
            run_dir, export_dir = output / f'runs/part_{part}', output / f'results/part_{part}'
            if run_dir.exists() or export_dir.exists():
                raise FileExistsError('Refusing to overwrite an existing run or result directory')
            record = {
                'part': part, 'stage': 'solver',
                'case_dir': source.relative_to(ROOT).as_posix(),
                'run_dir': run_dir.relative_to(ROOT).as_posix(),
                'export_dir': export_dir.relative_to(ROOT).as_posix(),
                'runtime_diagnostics_autoaccepted': False,
            }
            summary['cases'].append(record)
            save()
            print(f'Part {part}: invoking SOL144 once; run directory {run_dir}', flush=True)
            record['runner_returncode'] = run(source, run_dir, str(solver), license_server, threads)
            run_summary_path = run_dir / 'run-summary.json'
            raw = json.loads(run_summary_path.read_text(encoding='utf-8'))
            record['raw_solver_status'] = raw.get('status')
            record['raw_run_summary'] = raw
            record['run_summary_sha256'] = file_hash(run_summary_path)
            record['solver_outputs_complete'] = bool(raw.get('solver_completed') and raw.get('fresh_f06')
                                                     and not raw.get('fatal_messages'))
            record['solver_clean'] = bool(
                record['solver_outputs_complete'] and raw.get('status') == 'completed'
                and record['runner_returncode'] == 0 and raw.get('launcher_returncode') == 0
                and not raw.get('runtime_errors') and not raw.get('warning_messages'))
            record['requires_diagnostic_review'] = not record['solver_clean']
            save()
            if not record['solver_outputs_complete']:
                record.update(stage='postprocessing_skipped',
                              reason='Solver did not provide fresh complete output without fatal messages',
                              all_coarse_screens_pass=None)
                save()
                continue

            record['stage'] = 'export'
            save()
            exported = export_case(run_dir, export_dir)
            record['exporter_result'] = exported
            diagnostic = exported.get('diagnostics', {})
            if diagnostic.get('runtime_errors') or diagnostic.get('warning_messages'):
                record['requires_diagnostic_review'] = True
                record['solver_clean'] = False
            record['stage'] = 'check'
            save()
            checked = check(run_dir, export_dir)
            if not isinstance(checked.get('all_coarse_screens_pass'), bool):
                raise ValueError('Checker must return a boolean all_coarse_screens_pass')
            check_path = export_dir / 'composite_check.json'
            check_path.write_text(json.dumps(checked, indent=2, allow_nan=False) + '\n', encoding='utf-8')
            record.update(stage='complete', checker_result=checked,
                          checker_output=check_path.relative_to(ROOT).as_posix(),
                          checker_output_sha256=file_hash(check_path),
                          all_coarse_screens_pass=checked['all_coarse_screens_pass'])
            print(f"Part {part}: raw status={record['raw_solver_status']}; "
                  f"coarse screens={record['all_coarse_screens_pass']}; "
                  f"diagnostic review required={record['requires_diagnostic_review']}", flush=True)
            save()

        records = summary['cases']
        all_exported = len(records) == 2 and all(record['stage'] == 'complete' for record in records)
        all_clean = all_exported and all(record['solver_clean'] for record in records)
        all_screens = all(record['all_coarse_screens_pass'] for record in records) if all_exported else None
        review_required = not (all_clean and all_screens)
        summary.update(
            workflow_status=('incomplete_solver_outputs' if not all_exported else
                             'complete_with_review_required' if review_required else
                             'complete_clean_and_coarse_screens_pass'),
            all_results_exported=all_exported, all_solver_runs_clean=all_clean,
            all_coarse_screens_pass=all_screens, review_required=review_required,
            finished_utc=timestamp(), exit_code=2 if review_required else 0)
        save()
        return summary, summary['exit_code']
    except BaseException as error:
        # Keep failed/interrupt evidence, then re-raise visibly. No retry or
        # suppression of permission, license, solver or postprocessor errors.
        summary.update(workflow_status='interrupted_or_failed', review_required=True,
                       exception_type=type(error).__name__, exception=str(error),
                       finished_utc=timestamp())
        save()
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', type=Path, default=DEFAULT_CONFIG)
    parser.add_argument('--output', type=Path, required=True,
                        help='New, nonexistent directory under workspace/BWB/trials; contains cases/runs/results')
    parser.add_argument('--angle-deg', type=float,
                        help='Override regional warp reference angles, preserving configured ply offsets')
    parser.add_argument('--solver', default=DEFAULT_SOLVER)
    parser.add_argument('--license-server', default=DEFAULT_LICENSE_SERVER)
    parser.add_argument('--threads', type=int, choices=range(1, 17), default=2)
    args = parser.parse_args()
    summary, exit_code = run_trial(args.config, args.output, args.angle_deg,
                                   args.solver, args.license_server, args.threads)
    print(json.dumps({key: value for key, value in summary.items() if key != 'cases'}, indent=2))
    return exit_code


if __name__ == '__main__':
    raise SystemExit(main())
