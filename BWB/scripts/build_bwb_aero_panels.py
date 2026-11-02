"""Build independent BWB aerodynamic geometry from paired structured skin rows.

NASA_SHM bulk-card organization is a read-only reference. No NASA aerodynamic
numbers, loads, units, or source code are reused. The combined deck is intended
for preprocessor import; flight conditions and a SOL144 solution are not added.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

try:
    from .inspect_shell_bdf import load_model
except ImportError:
    from inspect_shell_bdf import load_model

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE = ROOT / 'BWB/models/skin_structured'
DEFAULT_OUTPUT = ROOT / 'BWB/models/aero_panels'
ZONES = ('main', 'transition', 'bend', 'winglet')
FLOW = np.array([1., 0., 0.])


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def provenance(path):
    path = Path(path).resolve()
    raw = path.read_bytes()
    # Git can legitimately change CRLF to LF on checkout. Compare text content
    # across hosts, while retaining the exact generation-time bytes as evidence.
    return {'path': path.relative_to(ROOT).as_posix(),
            'sha256': hashlib.sha256(raw.replace(b'\r\n', b'\n')).hexdigest(),
            'hash_mode': 'crlf_to_lf', 'raw_sha256_at_generation': digest(path)}


def real(value):
    """Ten significant digits, fitting even signed reals into a 16-column field."""
    return f'{float(value):.9E}'


def write_card(name, values):
    """Use legal small free fields or paired 16-column large fixed fields.

    Simcenter User Guide 2-5..2-8 limits physical lines to 80 columns and small
    fields to eight characters. Femap accepts longer free fields, but that is
    not proof of solver compatibility. Large fields are always paired.
    """
    fields = [str(v) if v is not None else '' for v in values]
    if any(len(f) > 8 for f in fields):
        if any(len(f) > 16 for f in fields):
            raise ValueError(f'{name}: a field exceeds 16 characters')
        fields += [''] * ((-len(fields)) % 8)
        return '\n'.join(f'{name + "*" if i == 0 else "*":<8}'
                         + ''.join(f'{f:>16}' for f in fields[i:i+4])
                         for i in range(0,len(fields),4)) + '\n'
    return '\n'.join(','.join([name if i == 0 else '+'] + fields[i:i+8])
                     for i in range(0, len(fields), 8)) + '\n'


def panel_normal(p1, p4):
    vector = np.cross(FLOW, np.asarray(p4) - p1)
    magnitude = np.linalg.norm(vector)
    if magnitude < 1e-10:
        raise ValueError('CAERO1 span is degenerate or parallel to flow')
    return vector / magnitude


def normalwash_from_mean_section(points, normal, xi):
    """Negative local-normal displacement slope at a normalized chord station.

PCHIP is differentiated with respect to normalized mean-section X and then
divided by the physical X chord. This includes camber and geometric incidence.
"""
    points = np.asarray(points, dtype=float)
    chord = float(points[-1, 0] - points[0, 0])
    if chord <= 0:
        raise ValueError('Mean section must run from leading to trailing edge')
    u = (points[:, 0] - points[0, 0]) / chord
    if np.any(np.diff(u) <= 0):
        raise ValueError('Mean-section X coordinates must increase strictly')
    if not 0 <= xi <= 1:
        raise ValueError('Chord sample outside [0, 1]')
    h = (points - points[0]) @ np.asarray(normal)
    return -float(PchipInterpolator(u, h).derivative()(xi)) / chord


def extract_structure(source, destination):
    """Preserve the original structural bulk entries without its static executive."""
    raw = Path(source).read_bytes()
    lines = raw.splitlines(keepends=True)
    start = next((i+1 for i, s in enumerate(lines) if s.strip().upper() == b'BEGIN BULK'), None)
    if start is None:
        raise ValueError(f'No BEGIN BULK in {source}')
    end = next((i for i in range(start, len(lines)) if lines[i].strip().upper().startswith(b'ENDDATA')), len(lines))
    Path(destination).write_bytes(b''.join(lines[start:end]))


def region_support(lower, upper, lower_ids, upper_ids, rows, normal):
    """Use a distributed set from one local-normal-positive physical skin.

The selection avoids almost coincident upper/lower projections. It adds no
rigid links or artificial stiffness. Local skin stress sensitivity still needs
assessment before using this initial spline selection for SHM predictions.
"""
    a, b = rows
    use_upper = float(np.mean((upper[a:b+1]-lower[a:b+1]) @ normal)) >= 0
    points, node_ids = (upper, upper_ids) if use_upper else (lower, lower_ids)
    row_stride = max(1, math.ceil((b-a)/10))
    chosen_rows = sorted(set(range(a, b+1, row_stride)) | {a, b})
    columns = sorted(set(range(0, points.shape[1], 4)) | {points.shape[1]-1})
    selected = list(dict.fromkeys(int(node_ids[j, i]) for j in chosen_rows for i in columns))
    coord = np.array([points[j, i] for j in chosen_rows for i in columns])
    yaxis = np.cross(normal, FLOW)
    projection = np.column_stack((coord @ FLOW, coord @ yaxis))
    if np.linalg.matrix_rank(projection-projection.mean(axis=0), tol=1e-9) != 2:
        raise ValueError('Spline support does not span a local surface')
    rounded = {tuple(np.round(p/1e-8).astype(np.int64)) for p in projection}
    if len(rounded) != len(selected):
        raise ValueError('Coincident spline support projections')
    return selected, ('upper' if use_upper else 'lower'), chosen_rows, columns


def build_part(part, source_root=DEFAULT_SOURCE, output_root=DEFAULT_OUTPUT, nchord=10):
    if nchord < 4:
        raise ValueError('Use at least four chordwise boxes for this demonstration')
    source = Path(source_root) / f'part_{part}'
    output = Path(output_root) / f'part_{part}'
    output.mkdir(parents=True, exist_ok=True)
    bdf = source / f'part_{part}_skin.bdf'
    layout_path = source / 'mesh_layout.json'
    model = load_model(bdf)
    layout = json.loads(layout_path.read_text(encoding='utf-8'))
    if any(n['cp'] != 0 for n in model['nodes'].values()):
        raise ValueError('Transform structural node coordinates before aerodynamic extraction')
    if any(model['duplicate_ids'].values()):
        raise ValueError('Duplicate structural IDs in source model')
    nodes = {nid: np.array(n['xyz'], dtype=float) for nid, n in model['nodes'].items()}
    base = int(part) * 1_000_000
    cid, pid = base+900001, base+900002
    root_z = float(nodes[layout['patch_node_ids']['main_lower'][0][0]][2])
    panels, boxes = [], []
    developed_area = projected_area = chord_second_moment = 0.

    for zone in ZONES:
        lower_ids = np.asarray(layout['patch_node_ids'][zone+'_lower'])
        upper_ids = np.asarray(layout['patch_node_ids'][zone+'_upper'])
        lower = np.array([[nodes[int(n)] for n in row] for row in lower_ids])
        upper = np.array([[nodes[int(n)] for n in row] for row in upper_ids])
        if lower.shape != upper.shape:
            raise ValueError(f'Unpaired upper/lower section rows in {zone}')
        mean = (lower+upper)/2
        le = mean[:, 0]
        chord = mean[:, -1, 0] - le[:, 0]
        cuts = list(range(len(le))) if zone == 'bend' else [0, len(le)-1]
        for subindex, (a, b) in enumerate(zip(cuts[:-1], cuts[1:])):
            # A positive local CAD-Y normal gives consistent conventions across
            # the two separate mirrored models; reverse span indexing as well.
            reverse = panel_normal(le[a], le[b])[1] < 0
            i1, i4 = (b, a) if reverse else (a, b)
            p1, p4 = le[i1], le[i4]
            c1, c4 = float(chord[i1]), float(chord[i4])
            if min(c1,c4) <= 0:
                raise ValueError('Nonpositive projected chord')
            normal = panel_normal(p1, p4)
            eid = base+10_000*len(panels)+1
            sid = eid+5000
            ns = b-a
            ids, side, support_rows, support_columns = region_support(
                lower, upper, lower_ids, upper_ids, (a,b), normal)
            name = f'{zone}_{subindex+1:02}' if zone == 'bend' else zone
            span = p4-p1
            distance = float(np.linalg.norm(span[1:]))
            developed_area += distance*(c1+c4)/2
            projected_area += abs(float(span[2]))*(c1+c4)/2
            chord_second_moment += abs(float(span[2]))*(c1*c1+c1*c4+c4*c4)/3
            panel = dict(eid=eid, name=name, zone=zone, source_row_range=[a,b],
                         span_reversed=bool(reverse), pid=pid, cp=0, igid=1,
                         nspan=ns, nchord=nchord, p1=p1.tolist(), p4=p4.tolist(),
                         x12=c1, x43=c4, normal=normal.tolist(),
                         box_ids=list(range(eid,eid+ns*nchord)), spline_id=sid, set_id=sid,
                         structural_node_ids=ids, support_skin=side,
                         support_rows=support_rows, support_columns=support_columns)
            panels.append(panel)

            span_positions = ((le[a:b+1]-p1) @ span) / float(span @ span)
            order = np.argsort(span_positions)
            span_positions = span_positions[order]
            if np.any(np.diff(span_positions) <= 0):
                raise ValueError('Nonmonotonic mean-section span positions')

            def point(eta, xi):
                return p1+eta*span+xi*((1-eta)*c1+eta*c4)*FLOW

            for j in range(ns):
                eta0, eta1, eta = j/ns, (j+1)/ns, (j+.5)/ns
                for i in range(nchord):
                    xi0, xi1, xi = i/nchord, (i+1)/nchord, (i+.75)/nchord
                    slopes = np.array([normalwash_from_mean_section(row, normal, xi)
                                       for row in mean[a:b+1]])[order]
                    w = float(np.interp(eta, span_positions, slopes))
                    corners = np.array([point(eta0,xi0),point(eta0,xi1),
                                        point(eta1,xi1),point(eta1,xi0)])
                    area = .5*np.linalg.norm(np.cross(corners[1]-corners[0],corners[2]-corners[0]))
                    area += .5*np.linalg.norm(np.cross(corners[2]-corners[0],corners[3]-corners[0]))
                    box_id = eid+i+nchord*j
                    boxes.append(dict(id=box_id, j_index=len(boxes)+1, panel_id=eid,
                                      span_index=j, chord_index=i, corners=corners.tolist(),
                                      normal=normal.tolist(), control_point=point(eta,xi).tolist(),
                                      pressure_point=point(eta,(i+.25)/nchord).tolist(),
                                      w2gj_rad=w, area_m2=float(area)))

    if len(set(b['id'] for b in boxes)) != len(boxes):
        raise ValueError('Aerodynamic box IDs overlap')
    all_le = np.array([p[k] for p in panels for k in ['p1','p4']])
    refs = dict(area_convention='horizontal_projection_including_winglet',
                sref_m2=projected_area, bref_m=float(np.ptp(all_le[:,2])),
                cref_m=chord_second_moment/projected_area, developed_area_m2=developed_area,
                scope='single supplied component, not full assembled aircraft')
    aero = ['$ BWB aerodynamic geometry; m-N-kg-s; no flight condition.\n',
            '$ ACS: X=CAD X, Y=-CAD Z, Z=CAD Y. No symmetry images.\n',
            write_card('CORD2R',[cid,0,*map(real,[0,0,root_z,0,1,root_z,1,0,root_z])]),
            write_card('AEROS',[cid,cid,real(refs['cref_m']),real(refs['bref_m']),real(refs['sref_m2']),0,0]),
            write_card('PAERO1',[pid])]
    spline = ['$ BWB local surface splines: initial support selection, no added rigid links.\n']
    for p in panels:
        aero += [f"$ REGION {p['name']} | boxes {p['box_ids'][0]}..{p['box_ids'][-1]}\n",
                 write_card('CAERO1',[p['eid'],pid,0,p['nspan'],nchord,None,None,1,
                                     *map(real,p['p1']),real(p['x12']),*map(real,p['p4']),real(p['x43'])])]
        spline += [f"$ {p['name']}: {p['support_skin']} skin, {len(p['structural_node_ids'])} support GRID points\n",
                   write_card('SPLINE1',[p['spline_id'],p['eid'],p['box_ids'][0],p['box_ids'][-1],
                                         p['set_id'],real(0),'IPS','BOTH']),
                   write_card('SET1',[p['set_id'],*p['structural_node_ids']])]
    dmi = ['$ W2GJ: negative mean-surface normal slope at each box 3/4-chord control point.\n',
           '$ Sorted by CAERO1 EID, then chord-fast box ID. Values in linearized radians.\n',
           write_card('DMI',['W2GJ',0,2,2,0,None,len(boxes),1])]
    for box in boxes:
        dmi += [f"$ j={box['j_index']} box={box['id']}\n",
                write_card('DMI',['W2GJ',1,box['j_index'],real(box['w2gj_rad'])])]
    files = dict(aero='aero_cards.bdf',spline='spline_cards.bdf',camber='dmicards.bdf',
                 structure='structural_cards.bdf',combined=f'part_{part}_aero_structure.bdf')
    for key, content in [('aero',aero),('spline',spline),('camber',dmi)]:
        (output/files[key]).write_text(''.join(content),encoding='ascii')
    extract_structure(bdf, output/files['structure'])
    combined = ('$ PREPROCESSOR GEOMETRY PACKAGE ONLY. No SOL, TRIM, load case or aerodynamic solution.\n'
                '$ Keep these authoritative bulk includes when preparing a future SOL144 deck.\n'
                'CEND\nBEGIN BULK\n'
                + ''.join(f"INCLUDE '{files[k]}'\n" for k in ['structure','aero','spline','camber'])
                + 'ENDDATA\n')
    (output/files['combined']).write_text(combined,encoding='ascii')
    nasa_files = [ROOT/'NASA_SHM/core/base'/name for name in
                  ['aero_cards.bdf','spline_cards.bdf','dmicards.bdf','sol144.bdf','trim_cards.bdf']]
    manifest = dict(schema_version=1,part=str(part),units='m-N-kg-s',
        status='generated aerodynamic input geometry; no SOL144 solution',
        source_files=[provenance(bdf),provenance(layout_path)],
        nasa_reference_files=[provenance(p) for p in nasa_files],
        nasa_reuse='Read-only bulk-card organization; no NASA numerical geometry, units, loads, or correction values copied.',
        generator=provenance(Path(__file__)),
        structural_counts=dict(nodes=len(model['nodes']),elements=len(model['elements'])),
        aerodynamic_coordinate_system=dict(id=cid,origin=[0.,0.,root_z],
            x_axis=[1.,0.,0.],y_axis=[0.,0.,-1.],z_axis=[0.,1.,0.]),
        reference_geometry=refs,panels=panels,boxes=boxes,files=files,
        w2gj_metadata=dict(sampling='box_three_quarter_chord_midspan',
            formula='-d(mean_normal_height)/dx',
            chord_interpolation='PCHIP of paired mean section normal height versus normalized mean X; analytic derivative / physical X chord',
            span_interpolation='linear interpolation by projection of source-row leading edges onto CAERO1 leading-edge segment',
            sorting='single IGID; increasing CAERO1 EID; chord index varies fastest',
            reference='Installed Simcenter Nastran 2406 aero.pdf PDF pages 21,24,43,272; QRG.pdf PDF page 1378',
            documentation_note='The generic aero.pdf PDF147 midchord wording conflicts with the j-set definition and SOL144 example; use the 3/4-chord control point.'),
        assumptions=['Separate components in original CAD coordinates, not an assembled full wing.',
            'One layer of piecewise planar CAERO1 reference surfaces; camber and incidence supplied through W2GJ.',
            'Uniform chordwise boxes; source span row counts retained. Density is not convergence-validated.',
            'IPS SPLINE1 uses distributed nodes on one local-normal-positive skin; no artificial stiffness is added.',
            'Spline load distribution and local skin strains require solver verification before SHM interpretation.',
            'No AERO velocity/density card, AESTAT, TRIM, gravity, 1g target or flight-case values are assumed.',
            'SREF includes the horizontal projection of this component including its winglet; coefficients are not whole-aircraft coefficients.'])
    (output/'aero_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(json.dumps({'part':str(part),'regions':len(panels),'boxes':len(boxes),
                      'spline_nodes':[len(p['structural_node_ids']) for p in panels],
                      'reference_geometry':refs,'w2gj_range_rad':[min(b['w2gj_rad'] for b in boxes),max(b['w2gj_rad'] for b in boxes)]},indent=2))
    return manifest


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parts',nargs='+',choices=['1','2'],default=['1','2'])
    parser.add_argument('--source',type=Path,default=DEFAULT_SOURCE)
    parser.add_argument('--output',type=Path,default=DEFAULT_OUTPUT)
    parser.add_argument('--chord-divisions',type=int,default=10)
    args=parser.parse_args()
    for part in args.parts:
        build_part(part,args.source,args.output,args.chord_divisions)


if __name__=='__main__':
    main()
