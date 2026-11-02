"""Export actual PCOMP/MAT8 ply results from a completed BWB SOL144 F06.

Request STRESS/STRAIN(PRINT,PLOT,SORT1,REAL,CENTER)=ALL, NOCOMPS=1,
and SOUT=YES for every explicitly listed PCOMP ply. No solver is started.
In-plane results are at ply midpoints; interlaminar columns describe the
interface above that ply. No isotropic stress reconstruction is performed.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import re
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from BWB.scripts.inspect_shell_bdf import cards, load_model, nastran_float
from BWB.scripts.export_sol144_strain import (
    END_MARKER, SUBCASE_HEADER, _expected_subcases, file_hash,
    principal_strains, read_displacements, require,
)


def read_composite_model(deck):
    """Read explicit full-stack PCOMP/MAT8 properties and element references."""
    deck = Path(deck)
    model = load_model(deck)
    require(model['nodes'] and model['elements'], 'No structural GRID/shell model')
    require(not any(model['duplicate_ids'].values()), 'Duplicate structural IDs')
    require(all(n['cp'] == n['cd'] == 0 for n in model['nodes'].values()),
            'Structural GRID input/output coordinates must be basic')
    properties, materials, orientations = {}, {}, {}
    for name, fields, location in cards(deck):
        f = lambda index: fields[index] if index < len(fields) else ''
        if name == 'MAT8':
            mid = int(f(1))
            require(mid not in materials, f'Duplicate MAT8 {mid}')
            materials[mid] = {key: nastran_float(f(index)) for index, key in (
                (2, 'E1_Pa'), (3, 'E2_Pa'), (4, 'nu12'), (5, 'G12_Pa'),
                (6, 'G1z_Pa'), (7, 'G2z_Pa'), (8, 'rho_kg_m3'))}
            # Simcenter 2406 QRG MAT8, printed pp. 15-20/21 (PDF 1982/1983).
            # Keep specified blanks distinct from Nastran's Xc=Xt, Yc=Yt
            # defaults. STRN=1 changes allowable units to strain, not Pa.
            specified = {key: nastran_float(f(index)) if f(index) else None
                         for index, key in ((12, 'Xt'), (13, 'Xc'), (14, 'Yt'),
                                            (15, 'Yc'), (16, 'S'))}
            require(all(value is None or math.isfinite(value) and value > 0
                        for value in specified.values()), f'Invalid MAT8 allowable: {mid}')
            if any(value is not None for value in specified.values()):
                require(all(specified[key] is not None for key in ('Xt', 'Yt', 'S')),
                        f'MAT8 {mid}: specifying any allowable requires Xt, Yt and S')
            require(not f(19) or nastran_float(f(19)) == 1.,
                    f'MAT8 {mid}: STRN must be blank (stress) or 1 (strain)')
            mode = 'strain' if f(19) else 'stress'
            effective = dict(specified)
            effective['Xc'] = specified['Xc'] if specified['Xc'] is not None else specified['Xt']
            effective['Yc'] = specified['Yc'] if specified['Yc'] is not None else specified['Yt']
            materials[mid].update(allowable_mode=mode, specified_allowables=specified,
                                  effective_allowables=effective,
                                  strengths_Pa=specified if mode == 'stress' else None,
                                  effective_strengths_Pa=effective if mode == 'stress' else None,
                                  strength_failure_assessment_performed=False)
            require(all(materials[mid][key] > 0 for key in ('E1_Pa', 'E2_Pa', 'G12_Pa')),
                    f'Nonpositive MAT8 in-plane stiffness: {mid}')
        elif name == 'PCOMP':
            pid = int(f(1))
            require(pid not in properties and pid not in model['properties'], f'Duplicate property {pid}')
            require(not f(8), f'PCOMP {pid}: exporter requires an explicit full stack, LAM blank')
            plies, previous_mid, previous_thickness = [], None, None
            for index in range(9, len(fields), 4):
                ply_fields = [f(index + offset) for offset in range(4)]
                if not any(ply_fields):
                    continue
                mid = int(ply_fields[0]) if ply_fields[0] else previous_mid
                thickness = nastran_float(ply_fields[1]) if ply_fields[1] else previous_thickness
                require(mid is not None and thickness is not None and thickness > 0,
                        f'PCOMP {pid}: missing material or positive ply thickness')
                require(ply_fields[3].upper() == 'YES', f'PCOMP {pid}: every ply needs SOUT=YES')
                plies.append(dict(ply_id=len(plies) + 1, material_id=mid, thickness_m=thickness,
                                  angle_deg=nastran_float(ply_fields[2])))
                previous_mid, previous_thickness = mid, thickness
            require(plies, f'No PCOMP plies: {location}')
            thickness = sum(p['thickness_m'] for p in plies)
            z = nastran_float(f(2)) if f(2) else -thickness / 2
            for ply in plies:
                ply.update(bottom_z_m=z, mid_z_m=z + ply['thickness_m'] / 2,
                           interface_above_z_m=z + ply['thickness_m'])
                z += ply['thickness_m']
            properties[pid] = dict(thickness_m=thickness, plies=plies, failure_theory=f(5))
        elif name in ('CQUAD4', 'CTRIA3'):
            orientation = f(7 if name == 'CQUAD4' else 6)
            orientations[int(f(1))] = orientation or '0.0'
    require(properties and materials, 'PCOMP and MAT8 input required')
    for eid, element in model['elements'].items():
        require(element['pid'] in properties, f'Shell {eid} does not reference a PCOMP')
        require(all(n in model['nodes'] for n in element['nodes']), f'Shell {eid} references missing GRID')
        for ply in properties[element['pid']]['plies']:
            require(ply['material_id'] in materials, f'Shell {eid} references a non-MAT8 ply material')
    return model, properties, materials, orientations


def read_ply_tables(path):
    """Parse real SORT1 layered-composite tables; reject duplicate/nonfinite rows."""
    output = {'stress': {}, 'strain': {}}
    subcase = section = element_type = None
    ended, last_eid = False, None
    with Path(path).open(errors='replace') as stream:
        for line_number, line in enumerate(stream, 1):
            require('FATAL' not in line.upper(), f'F06 fatal at line {line_number}')
            ended |= END_MARKER in line
            match = SUBCASE_HEADER.search(line)
            if match:
                subcase = int(match[1])
            compact = ''.join(line.split()).upper()
            if line.lstrip().startswith('***'):
                section = None
            if 'LAYEREDCOMPOSITEELEMENTS' in compact:
                section = 'strain' if 'STRAINSIN' in compact else 'stress' if 'STRESSESIN' in compact else None
                element_type = 'CQUAD4' if 'QUAD4' in compact else 'CTRIA3' if 'TRIA3' in compact else None
                require(section is not None and element_type is not None, 'Unsupported composite output table')
                last_eid = None
                continue
            if any(title in compact for title in ('DISPLACEMENTVECTOR', 'FORCESOFSINGLE',
                    'AERODYNAMIC', 'FAILUREINDICES', 'FAILUREINDEX')):
                section = None
            if not section:
                continue
            fields = line.split()
            if len(fields) == 12 and fields[0] == '0':
                fields = fields[1:]
            if len(fields) == 10 and fields[0].isdigit() and last_eid is not None:
                fields = [str(last_eid)] + fields
            if len(fields) < 2 or not fields[0].isdigit() or not fields[1].isdigit():
                continue
            require(len(fields) == 11, f'Unexpected composite numeric row at F06 line {line_number}')
            require(subcase is not None, 'Composite table lacks a subcase header')
            eid, ply = map(int, fields[:2])
            values = tuple(nastran_float(token) for token in fields[2:])
            require(all(math.isfinite(value) for value in values), 'Nonfinite ply result')
            key = (subcase, eid, ply)
            require(key not in output[section], f'Duplicate {section} result {key}')
            output[section][key] = dict(element_type=element_type, values=values)
            last_eid = eid
    require(ended and output['stress'] and output['strain'], 'Incomplete composite stress/strain/end-of-job output')
    return output


def validate_case(case_dir, job_name='bwb_sol144', root_tolerance=1e-10):
    case_dir = Path(case_dir).resolve()
    require(Path(job_name).name == job_name, 'Job name must not contain a directory')
    require(math.isfinite(root_tolerance) and root_tolerance >= 0, 'Invalid root tolerance')
    deck, f06 = case_dir / f'{job_name}.bdf', case_dir / f'{job_name}.f06'
    require(deck.is_file() and f06.is_file(), 'Need completed F06 and full input deck')
    model, properties, materials, orientations = read_composite_model(deck)
    results = read_ply_tables(f06)
    displacements = read_displacements(f06)
    subcases = sorted({s for s, _, _ in results['strain']})
    require(_expected_subcases(deck) == set(subcases), 'Composite results do not match requested subcases')
    expected = {(s, eid, ply['ply_id']) for s in subcases for eid, e in model['elements'].items()
                for ply in properties[e['pid']]['plies']}
    for quantity in ('stress', 'strain'):
        require(set(results[quantity]) == expected, f'Missing/extra composite {quantity} element/ply rows')
        for (_, eid, _), row in results[quantity].items():
            require(row['element_type'] == model['elements'][eid]['type'], 'Composite element type mismatch')
    expected_grid = {(s, nid) for s in subcases for nid in model['nodes']}
    structural = {key: value for key, value in displacements.items() if key in expected_grid}
    extra = {key: value for key, value in displacements.items() if key not in expected_grid}
    require(set(structural) == expected_grid, 'Missing structural displacement rows')
    if extra:
        aero_count = 0
        for name, fields, _ in cards(deck):
            if name == 'CAERO1':
                ns, nc = int(fields[4] or 0), int(fields[5] or 0)
                require(ns > 0 and nc > 0, 'Nonuniform CAERO1 extra-point output is unsupported')
                aero_count += ns * nc + (ns + 1) * (nc + 1)
        first = max(model['nodes']) + 1
        expected_extra = {(s, n) for s in subcases for n in range(first, first + aero_count)}
        require(aero_count > 0 and set(extra) == expected_extra, 'Unexpected nonstructural GRID output')
        require(all(abs(v[i]) <= root_tolerance for v in extra.values() for i in (0, 1, 3, 5)),
                'Unexpected DOFs on generated aerodynamic points')
    roots = model['constraints'].get(1, {})
    require(roots and all(set(dofs) == set(range(1, 7)) and all(v == 0 for v in dofs.values())
                          for dofs in roots.values()), 'SPC set 1 must fully clamp root nodes')
    require(all(n in model['nodes'] for n in roots), 'Root SPC references a missing node')
    root_max = max(abs(v) for s in subcases for n in roots for v in structural[(s, n)])
    require(root_max <= root_tolerance, 'Fixed root has nonzero displacement/rotation')
    for s in subcases:
        for quantity in ('stress', 'strain'):
            require(any(abs(v) > 0 for key, row in results[quantity].items() if key[0] == s
                        for v in row['values'][:3]), f'All-zero in-plane ply {quantity}: subcase {s}')
        require(any(abs(v) > 0 for key, values in structural.items() if key[0] == s for v in values),
                f'All-zero structural response: subcase {s}')
    hashes = {p.name: file_hash(p) for p in sorted(case_dir.iterdir()) if p.suffix.lower() in ('.bdf', '.inc')}
    run_path = case_dir / 'run-summary.json'
    run = json.loads(run_path.read_text(encoding='utf-8')) if run_path.is_file() else {}
    require(run.get('fresh_f06') is not False, 'Solver output was recorded as stale')
    require(not run.get('input_sha256') or run['input_sha256'] == hashes, 'Inputs changed after solve')
    diagnostics = dict(runtime_errors=list(run.get('runtime_errors', [])), warning_messages=[])
    for path in (f06, case_dir / f'{job_name}.log'):
        if path.is_file():
            for line in path.read_text(errors='replace').splitlines():
                require('FATAL' not in line.upper(), f'Fatal diagnostic in {path.name}')
                if 'WARNING MESSAGE' in line:
                    diagnostics['warning_messages'].append(line.strip())
                if path.suffix == '.log' and re.search(r'\bERROR\b', line, re.I):
                    diagnostics['runtime_errors'].append(line.strip())
    diagnostics = {key: sorted(set(values)) for key, values in diagnostics.items()}
    return dict(case_dir=case_dir, f06=f06, model=model, properties=properties, materials=materials,
                orientations=orientations, results=results, subcases=subcases, displacement=structural,
                extra=extra, roots=roots, root_max=root_max, root_tolerance=root_tolerance,
                input_sha256=hashes, run=run, diagnostics=diagnostics)


def write_csv(path, rows):
    with Path(path).open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def ply_extrema(rows):
    """Summarize actual ply results without inferring a failure criterion."""
    critical = max(rows, key=lambda row: abs(row['dominant_principal_microstrain']))
    return dict(ply_rows=len(rows), shell_elements=len({row['element_id'] for row in rows}),
                material_ids=sorted({row['material_id'] for row in rows}),
                property_ids=sorted({row['property_id'] for row in rows}),
                max_abs_principal_microstrain=abs(critical['dominant_principal_microstrain']),
                max_abs_principal_element_id=critical['element_id'],
                max_abs_principal_ply_id=critical['ply_id'],
                max_abs_epsilon_11=max(abs(row['epsilon_11']) for row in rows),
                max_abs_epsilon_22=max(abs(row['epsilon_22']) for row in rows),
                max_abs_gamma_12=max(abs(row['gamma_12']) for row in rows),
                max_abs_sigma_11_Pa=max(abs(row['sigma_11_Pa']) for row in rows),
                max_abs_sigma_22_Pa=max(abs(row['sigma_22_Pa']) for row in rows),
                max_abs_tau_12_Pa=max(abs(row['tau_12_Pa']) for row in rows))


def export_case(case_dir, output_dir, job_name='bwb_sol144', root_tolerance=1e-10):
    data = validate_case(case_dir, job_name, root_tolerance)
    output_dir = Path(output_dir).resolve()
    require(output_dir != data['case_dir'], 'Keep exported results separate from runtime inputs')
    model, rows = data['model'], []
    for key in sorted(data['results']['strain']):
        subcase, eid, ply_id = key
        element = model['elements'][eid]
        ply = data['properties'][element['pid']]['plies'][ply_id - 1]
        stress, strain = (data['results'][q][key]['values'] for q in ('stress', 'strain'))
        first, second, _ = principal_strains(*strain[:3])
        centroid = np.mean([model['nodes'][n]['xyz'] for n in element['nodes']], axis=0)
        row = dict(subcase=subcase, element_id=eid, element_type=element['type'], property_id=element['pid'],
                   ply_id=ply_id, material_id=ply['material_id'], ply_angle_deg=ply['angle_deg'],
                   element_theta_or_mcid=data['orientations'][eid], ply_thickness_m=ply['thickness_m'],
                   ply_mid_z_m=ply['mid_z_m'], interface_above_z_m=ply['interface_above_z_m'],
                   centroid_x_m=float(centroid[0]), centroid_y_m=float(centroid[1]), centroid_z_m=float(centroid[2]),
                   sigma_11_Pa=stress[0], sigma_22_Pa=stress[1], tau_12_Pa=stress[2],
                   interface_tau_xz_material_Pa=stress[3], interface_tau_yz_material_Pa=stress[4],
                   stress_principal_angle_deg=stress[5], stress_principal_major_Pa=stress[6],
                   stress_principal_minor_Pa=stress[7], stress_max_shear_Pa=stress[8],
                   epsilon_11=strain[0], epsilon_22=strain[1], gamma_12=strain[2],
                   interface_gamma_xz_material=strain[3], interface_gamma_yz_material=strain[4],
                   strain_principal_angle_deg=strain[5], strain_principal_major_reported=strain[6],
                   strain_principal_minor_reported=strain[7], strain_max_shear_reported=strain[8],
                   epsilon_11_microstrain=strain[0] * 1e6, epsilon_22_microstrain=strain[1] * 1e6,
                   gamma_12_microstrain=strain[2] * 1e6, principal_max_microstrain=first * 1e6,
                   principal_min_microstrain=second * 1e6,
                   dominant_principal_microstrain=(first if abs(first) >= abs(second) else second) * 1e6)
        rows.append(row)
    displacements = []
    for (s, nid), values in sorted(data['displacement'].items()):
        xyz = model['nodes'][nid]['xyz']
        displacements.append(dict(subcase=s, grid_id=nid, x_m=xyz[0], y_m=xyz[1], z_m=xyz[2],
                                  ux_m=values[0], uy_m=values[1], uz_m=values[2],
                                  rx_rad=values[3], ry_rad=values[4], rz_rad=values[5],
                                  translation_m=math.sqrt(sum(v * v for v in values[:3]))))
    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(output_dir / 'ply_stress_strain.csv', rows)
    write_csv(output_dir / 'grid_displacements.csv', displacements)
    files = ['ply_stress_strain.csv', 'grid_displacements.csv']
    if data['extra']:
        extra_rows = [dict(subcase=s, solver_output_point_id=n, t1_local_m=v[0], t2_local_m=v[1],
                           t3_local_m=v[2], r1_local_rad=v[3], r2_local_rad=v[4], r3_local_rad=v[5])
                      for (s, n), v in sorted(data['extra'].items())]
        write_csv(output_dir / 'solver_generated_aero_displacements.csv', extra_rows)
        files.append('solver_generated_aero_displacements.csv')
    cases = []
    for subcase in data['subcases']:
        selected = [r for r in rows if r['subcase'] == subcase]
        case = dict(subcase=subcase, **ply_extrema(selected),
                    max_translation_m=max(r['translation_m'] for r in displacements if r['subcase'] == subcase))
        for key, row_key in (('per_property_results', 'property_id'), ('per_material_results', 'material_id')):
            case[key] = {value: ply_extrema([row for row in selected if row[row_key] == value])
                         for value in sorted({row[row_key] for row in selected})}
        cases.append(case)
    diagnostic = bool(data['diagnostics']['runtime_errors'] or data['diagnostics']['warning_messages'])
    summary = dict(schema_version=2, status='exported_with_diagnostics' if diagnostic else 'exported',
                   diagnostic_only=diagnostic, solver_end_of_job=True, source_f06=str(data['f06']),
                   source_f06_sha256=file_hash(data['f06']), input_sha256=data['input_sha256'],
                   exporter_sha256=file_hash(__file__), structural_nodes=len(model['nodes']),
                   shell_elements=len(model['elements']), ply_rows=len(rows), subcases=cases,
                   pcomp_properties=data['properties'], mat8_materials=data['materials'],
                   root_constraint_set=1, root_fixed_nodes=len(data['roots']), max_abs_root_dof=data['root_max'],
                   root_tolerance=root_tolerance, all_expected_element_ply_results_present=True,
                   all_results_finite=True, nonzero_stress_strain_and_displacement_per_subcase=True,
                   units='m-N-kg-s; stresses Pa, strains dimensionless, rotations radians; microstrain = strain * 1e6.',
                   output_location='Element centroid; in-plane values at each ply midpoint. Interlaminar columns at the interface above that ply; final-ply interface values are nominally zero, subject to numerical roundoff.',
                   coordinates='11/22/12 in ply coordinates, with engineering shear strain. Interlaminar stress and strain xz/yz in element material coordinates, as labeled XZ-MAT/YZ-MAT in the actual F06.',
                   interpretation='Direct PCOMP ply output, no MAT1 or isotropic stress reconstruction; no composite strength/failure certification.',
                   strength_interpretation=('MAT8 allowables are preserved as specified, with optional compression defaults reported separately. '
                                            'All-null allowables mean unassigned, not zero strength. No failure index, utilization, '
                                            'isotropic yield criterion or strength pass/fail assessment is inferred.'),
                   recorded_launcher_returncode=data['run'].get('launcher_returncode'), diagnostics=data['diagnostics'],
                   files=files, csv_sha256={name: file_hash(output_dir / name) for name in files})
    (output_dir / 'composite_summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--job-name', default='bwb_sol144')
    parser.add_argument('--root-tolerance', type=float, default=1e-10)
    args = parser.parse_args()
    try:
        result = export_case(args.case_dir, args.output_dir, args.job_name, args.root_tolerance)
    except (ValueError, OSError, KeyError) as exc:
        parser.exit(1, f'Composite export rejected: {exc}\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
