"""Extract this project's F06 displacement and shell membrane-strain tables.

    python -m core.export_results workspace/cases/sample_1

No extrapolation, averaging between elements, or filling missing results with
zero is performed. This parser intentionally rejects fiber-distance output;
the project's default STRAIN,CENTER request outputs membrane strain/curvature.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re

from .run_local_nastran import inspect_results
from .paths import require_working_output


def read_tables(path):
    displacements, strains = {}, {}
    section, subcase, strain_mode = None, None, None
    element_type = None
    ended = False
    with Path(path).open(errors="replace") as source:
        for line in source:
            ended |= '* * * END OF JOB * * *' in line
            match = re.search(r'^0\s+SUBCASE\s+(\d+)\s*$', line)
            if match:
                subcase = int(match[1])
            compact = ''.join(line.split())
            if 'DISPLACEMENTVECTOR' in compact:
                section = 'displacement'
            elif 'STRAINSINQUADRILATERALELEMENTS(QUAD4)' in compact:
                section, element_type, strain_mode = 'strain', 'CQUAD4', None
            elif 'STRAINSINTRIANGULARELEMENTS(TRIA3)' in compact:
                section, element_type, strain_mode = 'strain', 'CTRIA3', None
            elif any(title in compact for title in ('FORCESOFSINGLE', 'STRESSESIN', 'STRAINSINBAR', 'ANALYSISSUMMARY')):
                section = None
            if section == 'strain':
                if 'CURVATURE' in line:
                    strain_mode = 'membrane_curvature'
                if 'DISTANCE' in line:
                    raise ValueError('Fiber-distance strain output requires a different extraction rule')
            fields = line.split()
            if section == 'displacement' and len(fields) == 8 and fields[1] == 'G' and fields[0].isdigit():
                values = tuple(map(float, fields[2:]))
                key = (subcase, int(fields[0]))
                if key in displacements:
                    raise ValueError(f'Duplicate displacement {key}')
                displacements[key] = values
            elif section == 'strain' and len(fields) == 10 and fields[0] == '0' and fields[1].isdigit():
                if strain_mode != 'membrane_curvature' or float(fields[2]) != 0:
                    raise ValueError('Unrecognized shell strain row; refusing to assume membrane strain')
                key = (subcase, int(fields[1]))
                if key in strains:
                    raise ValueError(f'Duplicate shell strain {key}')
                strains[key] = (element_type, *map(float, fields[3:6]))
    if not ended or not displacements or not strains:
        raise ValueError('Incomplete F06: end marker, displacements and shell strains are required')
    if not all(math.isfinite(v) for row in displacements.values() for v in row):
        raise ValueError('Nonfinite displacement')
    if not all(math.isfinite(v) for row in strains.values() for v in row[1:]):
        raise ValueError('Nonfinite strain')
    return displacements, strains


def export_case(case_dir):
    case_dir = require_working_output(case_dir)
    summary = inspect_results(case_dir)
    if not summary['solver_completed']:
        raise ValueError('Solver did not complete; no results exported')
    input_hashes = {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in case_dir.iterdir() if p.suffix.lower() in {'.bdf','.inc'}}
    manifest = case_dir/'case_manifest.json'
    labels = json.loads(manifest.read_text(encoding='utf-8')) if manifest.exists() else None
    if labels is not None and labels['case_input_sha256'] != input_hashes:
        raise ValueError('Damage labels do not match the solved inputs')
    saved_summary = case_dir/'run-summary.json'
    if saved_summary.exists():
        recorded = json.loads(saved_summary.read_text(encoding='utf-8'))
        if recorded.get('fresh_f06') is False or (recorded.get('input_sha256') and recorded['input_sha256'] != input_hashes):
            raise ValueError('Inputs or output freshness do not match the recorded solver run')
        summary['runtime_errors'] = sorted(set(summary['runtime_errors']+recorded.get('runtime_errors',[])))
        summary['launcher_returncode'] = recorded.get('launcher_returncode')
    displacement, strain = read_tables(case_dir/'sol144.f06')
    subcases = sorted({key[0] for key in displacement})
    case_control = case_dir/'subcases.inc'
    if case_control.exists():
        expected = {int(s) for s in re.findall(r'^\s*SUBCASE\s+(\d+)',case_control.read_text(),re.M|re.I)}
        if expected and expected != set(subcases):
            raise ValueError('Result subcases do not match the requested subcases.inc')
    element_ids = sorted({key[1] for key in strain})
    node_ids = {key[1] for key in displacement}
    expected_nodes, expected_shells = set(), set()
    with (case_dir/'FEM_only.bdf').open() as mesh:
        for line in mesh:
            card = line[:8].strip().rstrip('*')
            width = 16 if line[:8].strip().endswith('*') else 8
            if card == 'GRID':
                expected_nodes.add(int(line[8:8+width]))
            elif card in {'CQUAD4', 'CTRIA3'}:
                expected_shells.add(int(line[8:8+width]))
    if node_ids != expected_nodes or set(element_ids) != expected_shells:
        raise ValueError('Result IDs do not match all GRID and shell element IDs in FEM_only.bdf')
    if set(strain) != {(s, e) for s in subcases for e in element_ids}:
        raise ValueError('Missing shell strain rows between subcases')
    if set(displacement) != {(s, n) for s in subcases for n in node_ids}:
        raise ValueError('Missing displacement rows between subcases')
    output = case_dir/'results'
    output.mkdir(exist_ok=True)
    with (output/'displacements.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['subcase', 'node_id', 'ux', 'uy', 'uz', 'rx', 'ry', 'rz'])
        writer.writerows([*key, *values] for key, values in sorted(displacement.items()))
    with (output/'shell_membrane_strains.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['subcase', 'element_id', 'element_type', 'exx_element', 'eyy_element', 'gamma_xy_element'])
        writer.writerows([*key, *values] for key, values in sorted(strain.items()))
    # ID + 3 columns per subcase can be passed directly to report_figure.py.
    with (output/'strain_wide.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['ID']+[f'subcase_{s}_{v}' for s in subcases for v in ['exx_element', 'eyy_element', 'gamma_xy_element']])
        writer.writerows([eid]+[value for s in subcases for value in strain[(s,eid)][1:]] for eid in element_ids)
    with (case_dir/'sol144.f06').open('rb') as source:
        digest = hashlib.file_digest(source, 'sha256').hexdigest()
    metadata = {
        'source_f06':str(case_dir/'sol144.f06'), 'source_sha256':digest,
        'input_sha256':input_hashes,
        'solver_status':summary['status'], 'runtime_errors':summary['runtime_errors'],
        'launcher_returncode':summary.get('launcher_returncode'),
        'warning_messages':summary['warning_messages'],
        'subcases':subcases, 'nodes_per_subcase':len(node_ids), 'shell_elements_per_subcase':len(element_ids),
        'coordinates':'displacement: GRID output coordinate system; strain: element coordinate system',
        'strain_definition':'First row of STRAIN/CURVATURE output: membrane strain; gamma_xy is engineering shear strain as printed by Nastran.',
        'units':'Inherited from BDF. Strains dimensionless, rotations radians. No unit conversion.',
        'precision':'F06 printed precision, typically 7 significant digits; use OP2 for full precision in future batches.',
        'missing_data':'Only CQUAD4 and CTRIA3 membrane strains exported. CBAR absence is not zero strain.',
    }
    if labels is not None:
        (output/'damage_labels.json').write_text(json.dumps(labels,indent=2),encoding='utf-8')
        metadata['damage_labels'] = str(output/'damage_labels.json')
    (output/'metadata.json').write_text(json.dumps(metadata,indent=2),encoding='utf-8')
    return metadata


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case_dir', type=Path)
    args = parser.parse_args()
    print(json.dumps(export_case(args.case_dir), indent=2))
