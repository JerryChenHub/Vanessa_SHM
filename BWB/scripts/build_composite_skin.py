"""Assign four parameterized woven carbon/epoxy regions and build trial SOL144 decks."""
from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from pathlib import Path

import numpy as np

from BWB.paths import REPOSITORY_ROOT as ROOT, TRIALS_DIR, WORKSPACE_DIR
from .build_bwb_aero_panels import provenance, real, write_card
from .build_sol144_demo import build as build_flight
from .inspect_shell_bdf import cards, load_model

REGIONS = {'upper': 101, 'lower': 102, 'leading_edge': 103, 'trailing_edge': 104}
DEFAULT_CONFIG = ROOT / 'BWB/config/composite_skin.json'


def material_values(config):
    m = config['material']
    scale = float(m.get('stiffness_scale', 1.))
    names = ('E1_Pa', 'E2_Pa', 'G12_Pa', 'G13_Pa', 'G23_Pa')
    values = {key: float(m[key]) * scale for key in names}
    values.update(nu12=float(m['nu12']), density_kg_m3=float(m['density_kg_m3']))
    if not all(math.isfinite(v) for v in values.values()) or min(values[k] for k in names) <= 0:
        raise ValueError('Finite positive elastic moduli are required')
    if values['density_kg_m3'] <= 0 or not 0 <= values['nu12'] < 1:
        raise ValueError('Invalid material density or Poisson ratio')
    if 1 - values['nu12']**2 * values['E2_Pa'] / values['E1_Pa'] <= 0:
        raise ValueError('Orthotropic plane stress compliance is not positive definite')
    return values


def transformed_stiffness(material, angle_deg):
    e1, e2, nu, g = (material[k] for k in ('E1_Pa', 'E2_Pa', 'nu12', 'G12_Pa'))
    denominator = 1 - nu * nu * e2 / e1
    q11, q22, q12, q66 = e1 / denominator, e2 / denominator, nu * e2 / denominator, g
    m, n = math.cos(math.radians(angle_deg)), math.sin(math.radians(angle_deg))
    a = q11*m**4 + 2*(q12+2*q66)*m*m*n*n + q22*n**4
    b = q11*n**4 + 2*(q12+2*q66)*m*m*n*n + q22*m**4
    c = (q11+q22-4*q66)*m*m*n*n + q12*(m**4+n**4)
    d = (q11-q12-2*q66)*m**3*n - (q22-q12-2*q66)*m*n**3
    e = (q11-q12-2*q66)*m*n**3 - (q22-q12-2*q66)*m**3*n
    f = (q11+q22-2*q12-2*q66)*m*m*n*n + q66*(m**4+n**4)
    return np.array([[a,c,d],[c,b,e],[d,e,f]])


def laminate(region, material, sign):
    """Resolve a complete stack without silently overriding thickness inputs.

    A ply-angle list defines the count; otherwise ply_count or a per-ply
    thickness list defines it and all angle offsets default to zero. Redundant
    thickness specifications must agree within numerical roundoff.
    """
    if sign not in (-1, 1):
        raise ValueError('Part angle sign must be +1 or -1')
    count = None
    if 'ply_count' in region:
        value = region['ply_count']
        try:
            number = float(value)
        except (TypeError, ValueError) as exc:
            raise ValueError('ply_count must be a positive integer') from exc
        if isinstance(value, bool) or not math.isfinite(number) or number < 1 or not number.is_integer():
            raise ValueError('ply_count must be a positive integer')
        count = int(number)
    if 'ply_angles_deg' in region:
        offsets = region['ply_angles_deg']
        if not isinstance(offsets, (list, tuple)) or not offsets:
            raise ValueError('At least one woven ply angle is required')
        if count is not None and count != len(offsets):
            raise ValueError('ply_count conflicts with ply_angles_deg')
        count = len(offsets)
    else:
        if count is None and isinstance(region.get('ply_thicknesses_m'), (list, tuple)):
            count = len(region['ply_thicknesses_m'])
        if not count:
            raise ValueError('Provide ply_angles_deg, ply_count or nonempty ply_thicknesses_m')
        offsets = [0.] * count
    try:
        angles = [sign*(float(region.get('angle_deg', 0.)) + float(a)) for a in offsets]
    except (TypeError, ValueError) as exc:
        raise ValueError('Finite ply angles are required') from exc
    if not all(math.isfinite(value) for value in angles):
        raise ValueError('Finite ply angles are required')

    def positive(value):
        try:
            thickness = float(value)
        except (TypeError, ValueError) as exc:
            raise ValueError('Finite positive ply thicknesses are required') from exc
        if isinstance(value, bool) or not math.isfinite(thickness) or thickness <= 0:
            raise ValueError('Finite positive ply thicknesses are required')
        return thickness

    scalar = positive(region['ply_thickness_m']) if 'ply_thickness_m' in region else None
    requested_total = positive(region['total_thickness_m']) if 'total_thickness_m' in region else None
    if 'ply_thicknesses_m' in region:
        specified = region['ply_thicknesses_m']
        if not isinstance(specified, (list, tuple)) or len(specified) != count:
            raise ValueError('ply_thicknesses_m must contain one thickness per ply')
        thicknesses = [positive(value) for value in specified]
        if scalar is not None and not all(math.isclose(value, scalar, rel_tol=1e-12, abs_tol=1e-15)
                                          for value in thicknesses):
            raise ValueError('ply_thickness_m conflicts with ply_thicknesses_m')
    elif scalar is not None:
        thicknesses = [scalar] * count
    elif requested_total is not None:
        thicknesses = [requested_total / count] * count
    else:
        raise ValueError('Provide ply_thickness_m, ply_thicknesses_m or total_thickness_m')
    uniform = all(value == thicknesses[0] for value in thicknesses)
    # Preserve the original scalar multiply and integration order for defaults.
    total = thicknesses[0] * count if uniform else math.fsum(thicknesses)
    if not math.isfinite(total):
        raise ValueError('Finite total laminate thickness is required')
    if requested_total is not None and not math.isclose(total, requested_total, rel_tol=1e-12, abs_tol=1e-15):
        raise ValueError('total_thickness_m conflicts with the sum of ply thicknesses')
    a = np.zeros((3,3)); b = a.copy(); d = a.copy()
    z = -total/2
    for angle, thickness in zip(angles, thicknesses):
        top = z + thickness
        q = transformed_stiffness(material, angle)
        a += q*(top-z); b += q*(top**2-z**2)/2; d += q*(top**3-z**3)/3
        z = top
    result = dict(angles_deg=angles, ply_thicknesses_m=thicknesses, total_thickness_m=total,
                  A_N_per_m=a.tolist(), B_N=b.tolist(), D_Nm=d.tolist(),
                  areal_mass_kg_m2=material['density_kg_m3']*total,
                  symmetric=bool(np.allclose(angles, angles[::-1], atol=1e-12, rtol=0)
                                 and np.allclose(thicknesses, thicknesses[::-1], atol=1e-15, rtol=1e-12)))
    if uniform:
        result['ply_thickness_m'] = thicknesses[0]
    return result


def classify(xi, side, config):
    le, te = config['leading_edge_chord_fraction'], config['trailing_edge_chord_fraction']
    if xi <= le:
        return 'leading_edge'
    if xi >= 1-te:
        return 'trailing_edge'
    return side


def element_assignment(model, layout, config):
    nodes = {nid: np.array(n['xyz']) for nid, n in model['nodes'].items()}
    grids = {key: np.array(value) for key, value in layout['patch_node_ids'].items()}
    positions = {name: {int(n): (j,i) for j,row in enumerate(grid) for i,n in enumerate(row)}
                 for name,grid in grids.items()}
    result = {}
    for eid, element in model['elements'].items():
        patch = layout['element_patches'][str(eid)]
        xyz = np.array([nodes[n] for n in element['nodes']]); centre = xyz.mean(axis=0)
        chord_fractions = None
        if patch.endswith(('_upper', '_lower')):
            grid = grids[patch]
            loc = [positions[patch][nid] for nid in element['nodes']]
            row = float(np.mean([p[0] for p in loc]))
            j = min(int(row), len(grid)-2); t = row-j
            le = (1-t)*nodes[int(grid[j,0])] + t*nodes[int(grid[j+1,0])]
            te = (1-t)*nodes[int(grid[j,-1])] + t*nodes[int(grid[j+1,-1])]
            tangent = (nodes[int(grid[j+1,0])]+nodes[int(grid[j+1,-1])]
                       -nodes[int(grid[j,0])]-nodes[int(grid[j,-1])])/2
            xi = float((centre[0]-le[0])/(te[0]-le[0]))
            chord_fractions = [float((nodes[nid][0] - nodes[int(grid[r, 0])][0])
                                     / (nodes[int(grid[r, -1])][0] - nodes[int(grid[r, 0])][0]))
                               for nid, (r, _) in zip(element['nodes'], loc)]
            side = patch.rsplit('_',1)[1]
            region = classify(xi, side, config)
        elif patch.endswith('_trailing_edge'):
            zone = patch.removesuffix('_trailing_edge')
            grid = grids[zone+'_lower']
            indices = [positions[zone+'_lower'][n][0] for n in element['nodes'] if n in positions[zone+'_lower']]
            j = min(indices)
            tangent = nodes[int(grid[j+1,-1])] - nodes[int(grid[j,-1])]
            xi, region = 1., 'trailing_edge'
        elif patch == 'tip_cap':
            grid = grids['winglet_lower']
            le, te = nodes[int(grid[-1,0])], nodes[int(grid[-1,-1])]
            xi = float((centre[0]-le[0])/(te[0]-le[0]))
            tangent = te-le  # span is normal to the tip cap; use the chord reference.
            region = classify(xi, config['tip_cap_middle_region'], config)
        else:
            raise ValueError(f'Unclassified patch: {patch}')
        if not -.02 <= xi <= 1.02:
            raise ValueError(f'Invalid chord fraction: {eid}, {xi}')
        normal = np.cross(xyz[1]-xyz[0], xyz[2]-xyz[0]) if len(xyz)==3 else np.cross(xyz[2]-xyz[0], xyz[3]-xyz[1])
        normal /= np.linalg.norm(normal)
        tangent = tangent - normal*np.dot(tangent, normal)
        if np.linalg.norm(tangent) < 1e-10:
            raise ValueError(f'Degenerate material reference direction: {eid}')
        tangent /= np.linalg.norm(tangent)
        area = sum(np.linalg.norm(np.cross(xyz[k]-xyz[0],xyz[k+1]-xyz[0]))/2 for k in range(1,len(xyz)-1))
        result[eid] = dict(region=region, pid=REGIONS[region], patch=patch, chord_fraction=xi,
                           centroid_m=centre.tolist(), normal=normal.tolist(), reference_axis=tangent.tolist(),
                           reference_kind='tip_chord' if patch=='tip_cap' else 'outward_span',
                           area_m2=float(area), mcid=50000+eid)
        if chord_fractions is not None:
            result[eid].update(chord_node_fractions=chord_fractions, source_span_row=j)
    if set(result) != set(model['elements']) or set(v['region'] for v in result.values()) != set(REGIONS):
        raise ValueError('Four-region coverage is incomplete')
    return result


def report_band_boundaries(assignment, config):
    """Report actual band endpoints per surface patch, retaining centroid labels."""
    patches = {}
    for entry in assignment.values():
        if 'chord_node_fractions' in entry:
            rows = patches.setdefault(entry['patch'], {})
            rows.setdefault(entry['source_span_row'], []).append(entry)
    report = {}
    for patch, rows in sorted(patches.items()):
        boundaries = {}
        for region, target in (('leading_edge', config['leading_edge_chord_fraction']),
                               ('trailing_edge', config['trailing_edge_chord_fraction'])):
            samples = []
            missing = []
            for row, entries in sorted(rows.items()):
                values = [xi for entry in entries if entry['region'] == region
                          for xi in entry['chord_node_fractions']]
                if not values:
                    missing.append(row)
                    continue
                xi = max(values) if region == 'leading_edge' else min(values)
                fraction = xi if region == 'leading_edge' else 1 - xi
                samples.append(dict(source_span_row=row, boundary_xi=xi, achieved_fraction=fraction))
            fractions = [sample['achieved_fraction'] for sample in samples]
            positions = [sample['boundary_xi'] for sample in samples]
            boundaries[region] = dict(
                boundary='aft_edge' if region == 'leading_edge' else 'forward_edge',
                target_fraction=target, represented_span_rows=len(samples), missing_span_rows=missing,
                boundary_xi_min=min(positions) if positions else None,
                boundary_xi_max=max(positions) if positions else None,
                achieved_fraction_min=min(fractions) if fractions else None,
                achieved_fraction_max=max(fractions) if fractions else None,
                fraction_error_min=min(fractions)-target if fractions else None,
                fraction_error_max=max(fractions)-target if fractions else None,
                max_abs_fraction_error=max(abs(value-target) for value in fractions) if fractions else None)
        report[patch] = dict(source_span_rows=len(rows), **boundaries)
    return report


def build(config_path=DEFAULT_CONFIG, output=TRIALS_DIR/'composite_skin', angle_override=None):
    output = Path(output).resolve()
    if not output.is_relative_to(WORKSPACE_DIR):
        raise ValueError('Generated composite trial files must be in workspace/BWB')
    config = json.loads(Path(config_path).read_text(encoding='utf-8-sig'))
    le, te = (float(config[key]) for key in ('leading_edge_chord_fraction','trailing_edge_chord_fraction'))
    if not (0 < le < 1 and 0 < te < 1 and le+te < 1):
        raise ValueError('LE and TE fractions must be positive and sum to less than one')
    if set(config['regions']) != set(REGIONS) or config['tip_cap_middle_region'] not in ('upper','lower'):
        raise ValueError('Invalid region definitions')
    if angle_override is not None:
        for region in config['regions'].values():
            region['angle_deg'] = float(angle_override)
    material = material_values(config)
    stacks_by_part = {}
    for part in (1, 2):
        sign = float(config['part_angle_sign'][str(part)])
        stacks_by_part[part] = {name: laminate(value, material, sign) for name, value in config['regions'].items()}
    flight = config['flight']
    cases_root = output/'cases'
    if any(cases_root.glob('part_*/bwb_sol144.bdf')):
        raise FileExistsError('Choose a new trial output directory to preserve existing assignments')
    build_flight(cases_root, velocity=flight['velocity_m_s'], density=flight['density_kg_m3'],
                 sound_speed=flight['sound_speed_m_s'], alpha_deg=flight['alpha_deg'])
    for part in (1,2):
        source = ROOT/f'BWB/models/aero_panels/part_{part}/structural_cards.bdf'
        layout_path = ROOT/f'BWB/models/skin_structured/part_{part}/mesh_layout.json'
        model = load_model(source); layout = json.loads(layout_path.read_text())
        assignment = element_assignment(model, layout, config)
        stacks = stacks_by_part[part]
        chunks = ['$ Four-zone cured woven carbon/epoxy trial; SI units.\n']
        for name, fields, _ in cards(source):
            if name in ('MAT1','PSHELL'):
                continue
            values = fields[1:]
            if name in ('CQUAD4','CTRIA3'):
                entry = assignment[int(fields[1])]
                count = 4 if name=='CQUAD4' else 3
                values = [fields[1], entry['pid'], *fields[3:3+count], entry['mcid']]
            chunks.append(write_card(name, values))
        chunks.append(write_card('PARAM',['NOCOMPS',1]))
        chunks.append(write_card('MAT8',[1,real(material['E1_Pa']),real(material['E2_Pa']),real(material['nu12']),
                                        real(material['G12_Pa']),real(material['G13_Pa']),real(material['G23_Pa']),
                                        real(material['density_kg_m3'])]))
        for region,pid in REGIONS.items():
            stack=stacks[region]
            values=[pid,real(-stack['total_thickness_m']/2),0.,'', '',0.,0.,'']
            for angle, thickness in zip(stack['angles_deg'], stack['ply_thicknesses_m']):
                values.extend([1,real(thickness),real(angle),'YES'])
            chunks.append(write_card('PCOMP',values))
        for entry in assignment.values():
            centre=np.array(entry['centroid_m']); normal=np.array(entry['normal']); tangent=np.array(entry['reference_axis'])
            chunks.append(write_card('CORD2R',[entry['mcid'],0,*map(real,centre),*map(real,centre+normal),*map(real,centre+tangent)]))
        folder=cases_root/f'part_{part}'
        (folder/'structural_cards.bdf').write_text(''.join(chunks),encoding='ascii')
        deck=folder/'bwb_sol144.bdf'
        deck.write_text(deck.read_text().replace('ASSUMED LOW SPEED AEROELASTIC STRAIN DEMO','FOUR ZONE WOVEN CARBON EPOXY TRIAL')
                        .replace('REAL,VONMISES,CENTER','REAL,CENTER').replace('REAL,FIBER,CENTER','REAL,CENTER'),encoding='ascii')
        counts=Counter(entry['region'] for entry in assignment.values())
        mass=sum(entry['area_m2']*stacks[entry['region']]['areal_mass_kg_m2'] for entry in assignment.values())
        report=dict(part=part,config=config,material_effective=material,laminates=stacks,element_counts=dict(counts),
                    shell_mass_kg=mass,elements=assignment,source_structure=provenance(source),source_layout=provenance(layout_path),
                    band_boundaries=report_band_boundaries(assignment, config),
                    generator=provenance(Path(__file__)),
                    discretization='Element-centroid selection; bands follow existing elements, not exact CAD cuts. Refinement is needed to resolve arbitrary narrow widths.',
                    coordinate_convention='Each element MCID X follows outward span projected onto shell. Tip cap uses chord. PCOMP angles rotate woven warp; fill is +90 degrees. Part 2 angle sign mirrors the fabric path.',
                    inputs=[provenance(folder/name) for name in ('structural_cards.bdf','aero_cards.bdf','spline_cards.bdf','dmicards.bdf','flight_condition.bdf','bwb_sol144.bdf')])
        (folder/'composite_assignment.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        manifest=json.loads((folder/'case_manifest.json').read_text())
        manifest['assumptions']=[s for s in manifest['assumptions'] if 'Aluminium' not in s and 'Strains are linear' not in s]
        manifest['assumptions'] += ['MAT8 effective cured woven carbon/epoxy, four explicit PCOMP regions, no added epoxy layer.',
                                    'In-plane ply stress and strain are recovered at ply midpoints; no aluminium yield criterion applies.']
        manifest['inputs']=report['inputs']; manifest['composite_assignment']=provenance(folder/'composite_assignment.json')
        (folder/'case_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
        print(f'Part {part}: {dict(counts)}, shell mass {mass:.6f} kg',flush=True)
    return cases_root


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,default=DEFAULT_CONFIG)
    parser.add_argument('--output',type=Path,default=TRIALS_DIR/'composite_skin')
    parser.add_argument('--angle-deg',type=float,help='Override reference warp angle for all four regions; preserve each ply offset')
    args=parser.parse_args()
    print(build(args.config,args.output,args.angle_deg))
