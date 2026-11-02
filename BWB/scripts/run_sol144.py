"""Run a BWB case in BWB scratch space; reuse only NASA's read-only log inspector."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

from BWB.paths import REPOSITORY_ROOT as ROOT, WORKSPACE_DIR
from NASA_SHM.core.run_local_nastran import inspect_results
from .build_bwb_aero_panels import provenance


def run(source, destination, solver, license_server, threads=2):
    source, destination = Path(source).resolve(strict=True), Path(destination).resolve()
    if not destination.is_relative_to(WORKSPACE_DIR):
        raise ValueError('Solver outputs must be in workspace/BWB')
    if destination == source:
        raise ValueError('Keep permanent input separate from runtime files')
    destination.mkdir(parents=True, exist_ok=True)
    if (destination / 'bwb_sol144.f06').exists():
        raise FileExistsError('Use a fresh run directory to preserve prior solver evidence')
    for path in source.iterdir():
        if path.suffix.lower() in ('.bdf', '.inc', '.json'):
            shutil.copyfile(path, destination / path.name)
    scratch = destination / 'scratch'
    scratch.mkdir(exist_ok=True)
    solver = Path(solver).resolve(strict=True)
    command = [str(solver), 'bwb_sol144.bdf', f'auth={license_server}', 'program=NASTRAN',
               'memory=1GB', f'parallel={threads}', f'sdirectory={scratch}',
               'scr=yes', 'notify=no', 'old=yes']
    env = os.environ.copy()
    env.update(SPLM_LICENSE_SERVER=license_server, OMP_NUM_THREADS=str(threads), MKL_NUM_THREADS=str(threads))
    env.pop('NOLICENSECHECK', None)
    started = time.time()
    metadata = {'command': command, 'started_unix': started,
                'source_case': source.relative_to(ROOT).as_posix(),
                'input_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                 for p in destination.iterdir() if p.suffix.lower() in ('.bdf', '.inc')},
                'reused_inspector': provenance(ROOT / 'NASA_SHM/core/run_local_nastran.py'),
                'runner': provenance(Path(__file__))}
    (destination / 'run-invocation.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')
    with (destination / 'launcher.stdout.txt').open('wb') as stdout, (destination / 'launcher.stderr.txt').open('wb') as stderr:
        result = subprocess.run(command, cwd=destination, env=env, stdout=stdout, stderr=stderr,
                                creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0))
    metadata['launcher_returncode'] = result.returncode
    f06 = destination / 'bwb_sol144.f06'
    metadata['fresh_f06'] = f06.exists() and f06.stat().st_mtime >= started - 2
    summary = inspect_results(destination, 'bwb_sol144')
    # NASA's inspector recognizes unlabelled headers; BWB uses descriptive LABELs.
    if f06.exists():
        summary['subcases_in_output'] = sorted({int(value) for value in re.findall(
            r'^0[^\r\n]*\bSUBCASE\s+(\d+)\s*$', f06.read_text(errors='replace'), re.MULTILINE)})
    if not metadata['fresh_f06']:
        summary.update(status='failed_or_incomplete', solver_completed=False)
    if result.returncode:
        summary['runtime_errors'].append(f'Launcher returned exit code {result.returncode}')
        summary['status'] = 'completed_with_runtime_errors' if summary['solver_completed'] else 'failed_or_incomplete'
    summary.update(metadata)
    (destination / 'run-summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(summary, indent=2))
    return 0 if summary['status'] == 'completed' else 2


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('destination', type=Path)
    parser.add_argument('--solver', default=r'D:\Siemens\Femap 2406\nastran\bin\nastran.exe')
    parser.add_argument('--license-server', default='29000@license-siemens.eng.uci.edu')
    parser.add_argument('--threads', type=int, choices=range(1, 17), default=2)
    args = parser.parse_args()
    raise SystemExit(run(args.source, args.destination, args.solver, args.license_server, args.threads))
