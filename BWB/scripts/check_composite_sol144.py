"""Check actual composite SOL144 load balance and scales without a metal yield rule."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import re

import numpy as np

from .check_aero_panels import grouped_cards, read_aero_geometry
from .check_sol144_reasonableness import (digest, pressure_consistency, reaction_consistency,
                                        read_load_tables, require, shift_moment)
from .inspect_shell_bdf import cards, load_model, nastran_float


def constitutive_residual(rows, materials):
    """Compare independent ply stress/strain tables against orthotropic compliance."""
    error = scale = 0.
    for row in rows:
        m = materials[str(row['material_id'])]
        sx, sy, tau = (float(row[k]) for k in ('sigma_11_Pa', 'sigma_22_Pa', 'tau_12_Pa'))
        predicted = np.array([(sx-m['nu12']*sy)/m['E1_Pa'],
                              sy/m['E2_Pa']-m['nu12']*sx/m['E1_Pa'], tau/m['G12_Pa']])
        actual = np.array([float(row[k]) for k in ('epsilon_11','epsilon_22','gamma_12')])
        require(np.all(np.isfinite(predicted)) and np.all(np.isfinite(actual)), 'Nonfinite ply data')
        error = max(error, float(np.max(abs(predicted-actual))))
        scale = max(scale, float(np.max(abs(actual))))
    require(scale > 0, 'Zero strain field')
    return dict(max_abs_strain_difference=error, max_strain_scale=scale, relative_error=error/scale,
                tolerance=1e-5, passed=error/scale < 1e-5)


def coupling_evidence(case_dir, job_name):
    f04 = (case_dir/f'{job_name}.f04').read_text(errors='replace')
    matrix = [line.strip() for line in f04.splitlines() if re.search(r'\bKALL\s+\d+\s+\d+', line)]
    nonzero = False
    for line in matrix:
        fields = line.split(); offset = fields.index('KALL')
        # F04 matrix synopsis: name, rows, columns, form, type, words, density.
        nonzero |= (int(fields[offset+1]) > 0 and int(fields[offset+2]) > 0
                    and nastran_float(fields[offset+6]) > 0)
    f06 = (case_dir/f'{job_name}.f06').read_text(errors='replace')
    marker = 'AERODYNAMIC MONITOR POINT TOTAL VEHICLE COEFFICIENTS:'
    require(marker in f06, 'Missing aeroelastic monitor table')
    match = re.search(r'^\s*WIND\s+CL\s+(\S+)\s+(\S+).*?(\S+)\s*$',
                      f06.split(marker, 1)[1], re.MULTILINE)
    require(match is not None, 'Missing wind-axis CL in aeroelastic monitor')
    rigid, increment, total = map(nastran_float, match.groups())
    return dict(nonzero_KALL=bool(nonzero), KALL_matrix_records=matrix,
                rigid_CL=rigid, elastic_CL_increment=increment, total_CL=total,
                nonzero_elastic_lift_increment=abs(increment) > 1e-10,
                scope='The current builder creates one TRIM subcase per case.')


def check(case_dir, export_dir, job_name='bwb_sol144'):
    case_dir, export_dir = Path(case_dir).resolve(), Path(export_dir).resolve()
    summary_path = export_dir/'composite_summary.json'
    exported = json.loads(summary_path.read_text(encoding='utf-8'))
    f06 = case_dir/f'{job_name}.f06'
    require(exported['source_f06_sha256'] == digest(f06), 'Export belongs to another F06')
    inputs = {p.name:digest(p) for p in case_dir.iterdir() if p.suffix.lower() in ('.bdf','.inc')}
    require(inputs == exported['input_sha256'], 'Case inputs changed after export')
    for name, expected in exported['csv_sha256'].items():
        require(digest(export_dir/name) == expected, f'Exported CSV changed: {name}')
    model = load_model(case_dir/'structural_cards.bdf')
    roots = sorted(model['constraints'].get(1, {}))
    require(roots and all(n['cp'] == n['cd'] == 0 for n in model['nodes'].values()),
            'Expected fixed roots and basic structural coordinates')
    groups = grouped_cards(case_dir/'aero_cards.bdf')
    aero = read_aero_geometry(groups); boxes = aero['boxes']
    af = groups['AEROS'][0][0]; origin = aero['frames'][int(af[2] or 0)][0]
    trim = [f for name,f,_ in cards(case_dir/'flight_condition.bdf') if name == 'TRIM']
    require(len(trim) == 1, 'Expected one TRIM condition')
    q = nastran_float(trim[0][3]); alpha = nastran_float(trim[0][trim[0].index('ANGLEA')+1])
    require(q > 0, 'Nonpositive dynamic pressure')
    pressures, resultants, reactions, precision = read_load_tables(f06, include_precision=True)
    subcases = sorted(resultants)
    require(len(subcases) == 1, 'This checker expects one subcase per trial')
    require(set(pressures) == {(s,b) for s in subcases for b in boxes}, 'Incomplete pressure coverage')
    require(set(reactions) == {(s,n) for s in subcases for n in roots}, 'Incomplete root reactions')
    with (export_dir/'grid_displacements.csv').open(newline='') as stream:
        displacements = list(csv.DictReader(stream))
    with (export_dir/'ply_stress_strain.csv').open(newline='') as stream:
        plies = list(csv.DictReader(stream))
    constitutive = constitutive_residual(plies, exported['mat8_materials'])
    coupling = coupling_evidence(case_dir, job_name)
    area, span = (aero['reference'][k] for k in ('area_m2','span_m'))
    aspect = span*span/area
    slope = 2*math.pi*aspect/(2+math.sqrt(4+aspect*aspect))
    rough_lift = q*area*slope*alpha
    reports = []
    for subcase in subcases:
        force, moment = np.zeros(3), np.zeros(3); pressure_ratios = []
        for bid, box in boxes.items():
            cp, pressure = pressures[(subcase,bid)]
            residual, tolerance = pressure_consistency(cp,pressure,q,*precision['pressures'][(subcase,bid)])
            pressure_ratios.append(residual/tolerance)
            corners = box['corners']
            point = .375*(corners[0]+corners[3])+.125*(corners[1]+corners[2])
            load = pressure*box['area']*box['normal']
            force += load; moment += np.cross(point,load)
        listed = resultants[subcase]
        difference, tolerance = reaction_consistency(
            [model['nodes'][nid]['xyz'] for nid in roots],
            [reactions[(subcase,nid)] for nid in roots],
            [precision['reactions'][(subcase,nid)] for nid in roots],
            listed, precision['resultants'][subcase])
        root_moment = shift_moment(moment,force,origin)
        reaction_moment = shift_moment(listed[3:],listed[:3],origin)
        force_error = float(np.linalg.norm(force+listed[:3]))
        moment_error = float(np.linalg.norm(root_moment+reaction_moment))
        force_relative = force_error/max(float(np.linalg.norm(force)),q*area*1e-6)
        moment_relative = moment_error/max(float(np.linalg.norm(root_moment)),q*area*span*1e-6)
        selected = [row for row in displacements if int(row['subcase']) == subcase]
        require({int(row['grid_id']) for row in selected} == set(model['nodes']), 'Incomplete displacement IDs')
        dofs = np.array([[float(row[k]) for k in ('ux_m','uy_m','uz_m','rx_rad','ry_rad','rz_rad')] for row in selected])
        require(np.all(np.isfinite(dofs)), 'Nonfinite displacement field')
        translation = float(np.max(np.linalg.norm(dofs[:,:3],axis=1)))
        rotation = float(np.max(np.linalg.norm(dofs[:,3:],axis=1)))
        root_dof = max(float(np.max(abs(value))) for row,value in zip(selected,dofs) if int(row['grid_id']) in roots)
        lift = float(force@aero['up']); lift_ratio = lift/rough_lift if rough_lift else None
        screens = dict(force_balance=force_relative < 1e-4, moment_balance=moment_relative < 1e-4,
                       pressure_equals_qCp=max(pressure_ratios) <= 1,
                       printed_reactions_consistent=bool(np.all(abs(difference) <= tolerance)),
                       root_fixed=root_dof <= 1e-10, small_translation=translation/span < .01,
                       small_rotation=rotation < .05, orthotropic_stress_strain_consistent=constitutive['passed'],
                       aeroelastic_stiffness_active=coupling['nonzero_KALL'],
                       aeroelastic_lift_changes=coupling['nonzero_elastic_lift_increment'],
                       same_lift_order_as_finite_wing=lift_ratio is not None and .25 < lift_ratio < 4)
        reports.append(dict(subcase=subcase,screens=screens, aerodynamic_force_CAD_N=force.tolist(),
            aerodynamic_moment_root_Nm=root_moment.tolist(),root_origin_CAD_m=origin.tolist(),
            root_reaction_force_CAD_N=listed[:3].tolist(),root_reaction_moment_Nm=reaction_moment.tolist(),
            force_balance_relative=force_relative,moment_balance_relative=moment_relative,
            force_balance_error_N=force_error,moment_balance_error_Nm=moment_error,
            pressure_max_error_over_print_tolerance=max(pressure_ratios),
            reactions_max_error_over_print_tolerance=float(np.max(abs(difference)/tolerance)),
            max_translation_m=translation,max_translation_over_span=translation/span,
            max_rotation_rad=rotation,max_abs_root_dof=root_dof,body_up_lift_N=lift,
            body_up_CL=lift/(q*area),lift_to_rough_finite_wing_ratio=lift_ratio))
    passed = all(all(r['screens'].values()) for r in reports)
    return dict(status='coarse_checks_pass' if passed else 'review_needed',all_coarse_screens_pass=passed,
                diagnostic_only=exported['diagnostic_only'],runtime_diagnostics=exported['diagnostics'],
                source_f06_sha256=digest(f06),export_summary_sha256=digest(summary_path),checker_sha256=digest(__file__),
                pressure_boxes=len(boxes),root_reaction_nodes=len(roots),constitutive_consistency=constitutive,
                aeroelastic_coupling=coupling,subcases=reports,
                rough_finite_wing=dict(lift_N=rough_lift,formula='q*S*2*pi*AR/(2+sqrt(4+AR^2))*alpha',
                                      scope='Isolated unswept uncambered finite-wing scale; ignores winglet, camber and sweep.'),
                limitations=['Physical checks do not clear raw solver runtime diagnostics.',
                             'No mesh convergence, buckling, strength or experimental validation is claimed.',
                             'Strain output is at ply midpoints, not an outer-fiber failure envelope.'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case-dir',type=Path,required=True)
    parser.add_argument('--export-dir',type=Path,required=True)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    result = check(args.case_dir,args.export_dir)
    text = json.dumps(result,indent=2)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(text,encoding='utf-8')
    print(text,end='')
    return 0 if result['all_coarse_screens_pass'] else 2


if __name__ == '__main__':
    raise SystemExit(main())
