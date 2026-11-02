"""Plot actual four-region composite assignment and first-ply material directions.

Reads part 1 of the baseline and angle_30 temporary cases. Isotropic regions
have no directional arrows. The only output is assignment_overview.png in the
selected temporary trial directory. No Femap is used.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

try:
    from .inspect_shell_bdf import load_model
except ImportError:
    from inspect_shell_bdf import load_model

ROOT = Path(__file__).resolve().parents[2]
TRIAL_ROOT = ROOT / 'workspace/BWB/trials/composite_skin'
REGIONS = ('upper', 'lower', 'leading_edge', 'trailing_edge')
COLORS = {'upper': '#4b91b6', 'lower': '#e5a664',
          'leading_edge': '#65b699', 'trailing_edge': '#ab8ac5'}
WARP, FILL = '#234e7a', '#d26936'


def material_model(laminate):
    """Legacy reports described four orthotropic woven-fabric regions."""
    model = laminate.get('material_model', 'orthotropic')
    if model not in ('orthotropic', 'plane_isotropic'):
        raise ValueError(f'Unsupported material model in assignment report: {model!r}')
    return model


def directional_regions(assignment):
    return {region for region, laminate in assignment['laminates'].items()
            if material_model(laminate) == 'orthotropic'}


def camera_matrix(elevation=45.0, azimuth=-64.0):
    """Orthographic side/upper view; convert CAD XYZ to span/chord/height first."""
    elevation, azimuth = np.deg2rad([elevation, azimuth])
    toward_viewer = np.array([np.cos(elevation) * np.cos(azimuth),
                              np.cos(elevation) * np.sin(azimuth), np.sin(elevation)])
    right = np.array([-np.sin(azimuth), np.cos(azimuth), 0.0])
    up = np.cross(toward_viewer, right)
    # Each row maps CAD [X, Y, Z] into screen-right, screen-up, camera depth.
    return np.array([right, up, toward_viewer])[:, [1, 2, 0]]


def read_case(folder: Path):
    assignment = json.loads((folder / 'composite_assignment.json').read_text(encoding='utf-8'))
    model = load_model(folder / 'structural_cards.bdf')
    if any(node['cp'] != 0 for node in model['nodes'].values()) or any(model['duplicate_ids'].values()):
        raise ValueError('Plot requires unique entities and basic GRID coordinates')
    elements = {int(eid): record for eid, record in assignment['elements'].items()}
    if set(elements) != set(model['elements']):
        raise ValueError('Assignment and structural BDF do not contain the same shell elements')
    polygons = {eid: np.array([model['nodes'][nid]['xyz'] for nid in element['nodes']])
                for eid, element in model['elements'].items()}
    for eid, record in elements.items():
        normal, reference = np.array(record['normal']), np.array(record['reference_axis'])
        if (record['region'] not in REGIONS or not np.isclose(np.linalg.norm(normal), 1.0)
                or not np.isclose(np.linalg.norm(reference), 1.0)
                or abs(normal @ reference) > 1e-7):
            raise ValueError(f'Invalid material tangent frame for element {eid}')
        if not np.allclose(polygons[eid].mean(axis=0), record['centroid_m'], atol=1e-7, rtol=0):
            raise ValueError(f'Assignment centroid differs from BDF element {eid}')
    return assignment, elements, polygons


def fabric_directions(record, angle_deg):
    x = np.asarray(record['reference_axis'])
    normal = np.asarray(record['normal'])
    y = np.cross(normal, x)
    angle = np.deg2rad(angle_deg)
    warp = np.cos(angle) * x + np.sin(angle) * y
    fill = -np.sin(angle) * x + np.cos(angle) * y
    return warp, fill


def visible_samples(elements, polygons, camera, maximum=22, minimum_spacing=.075,
                    eligible_regions=None):
    """Select sparse front-surface anchors, checking occlusion in the actual mesh."""
    projected = {eid: points @ camera.T for eid, points in polygons.items()}
    triangles = []
    for points in projected.values():
        triangles.append(points[[0, 1, 2]])
        if len(points) == 4:
            triangles.append(points[[0, 2, 3]])
    triangles = np.asarray(triangles)
    a, b, c = triangles[:, 0], triangles[:, 1], triangles[:, 2]
    v0, v1 = b[:, :2] - a[:, :2], c[:, :2] - a[:, :2]
    denominator = v0[:, 0] * v1[:, 1] - v1[:, 0] * v0[:, 1]
    nondegenerate = np.abs(denominator) > 1e-12
    denominator = np.where(nondegenerate, denominator, 1.0)
    candidates = []
    for eid, record in sorted(elements.items()):
        if eligible_regions is not None and record['region'] not in eligible_regions:
            continue
        if np.asarray(record['normal']) @ camera[2] <= .12:
            continue
        point = np.asarray(record['centroid_m']) @ camera.T
        v2 = point[:2] - a[:, :2]
        u = (v2[:, 0] * v1[:, 1] - v1[:, 0] * v2[:, 1]) / denominator
        v = (v0[:, 0] * v2[:, 1] - v2[:, 0] * v0[:, 1]) / denominator
        covers = nondegenerate & (u >= -1e-7) & (v >= -1e-7) & (u + v <= 1 + 1e-7)
        depths = a[:, 2] + u * (b[:, 2] - a[:, 2]) + v * (c[:, 2] - a[:, 2])
        if np.any(covers) and np.max(depths[covers]) - point[2] < .0007:
            candidates.append((eid, point[:2]))
    if not candidates:
        return [], projected
    points = np.asarray([point for _, point in candidates])
    selected = [int(np.argmin(np.linalg.norm(points - points.mean(axis=0), axis=1)))]
    while len(selected) < maximum:
        nearest = np.linalg.norm(points[:, None] - points[selected][None], axis=2).min(axis=1)
        next_index = int(np.argmax(nearest))
        if nearest[next_index] < minimum_spacing:
            break
        selected.append(next_index)
    return [candidates[index][0] for index in selected], projected


def render(trial_root=TRIAL_ROOT):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import to_rgb
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    trial_root = Path(trial_root).resolve(strict=True)
    if not trial_root.is_relative_to((ROOT / 'workspace/BWB/trials').resolve()):
        raise ValueError('Composite assignment plots belong under workspace/BWB/trials')
    baseline, elements, polygons = read_case(trial_root / 'baseline/cases/part_1')
    rotated, rotated_elements, rotated_polygons = read_case(trial_root / 'angle_30/cases/part_1')
    if set(elements) != set(rotated_elements):
        raise ValueError('Comparison cases have different element IDs')
    for eid in elements:
        if (elements[eid]['region'] != rotated_elements[eid]['region']
                or not np.array_equal(polygons[eid], rotated_polygons[eid])
                or not np.allclose(elements[eid]['reference_axis'], rotated_elements[eid]['reference_axis'])
                or not np.allclose(elements[eid]['normal'], rotated_elements[eid]['normal'])):
            raise ValueError('Geometry, regions or reference frames differ between comparison cases')

    camera = camera_matrix()
    directional = directional_regions(baseline) | directional_regions(rotated)
    sampled, projected = visible_samples(elements, polygons, camera,
                                         eligible_regions=directional)
    drawing_order = sorted(projected, key=lambda eid: projected[eid][:, 2].mean())
    all_points = np.concatenate(list(projected.values()))[:, :2] * 1000
    low, high = all_points.min(axis=0), all_points.max(axis=0)
    extent = high - low
    margin = np.array([.07, .12]) * max(extent)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'text.color': '#263c50', 'axes.titlecolor': '#263c50'})
    fig, axes = plt.subplots(1, 3, figsize=(18.0, 7.0), facecolor='#f6f8fb')
    fig.subplots_adjust(left=.025, right=.985, bottom=.29, top=.81, wspace=.035)
    fig.text(.035, .94, 'BWB composite skin | Regions and material directions',
             fontsize=24, weight='bold')
    fig.text(.035, .89, f'Part 1  /  Actual BDF shell geometry  /  {len(elements):,} elements  /  '
             'Same upper-side orthographic view in all panels', fontsize=12, color='#64758a')

    def angle_label(assignment):
        angles = {data['angles_deg'][0] for data in assignment['laminates'].values()
                  if material_model(data) == 'orthotropic'}
        return '/'.join(f'{value:g}°' for value in sorted(angles)) or 'isotropic'

    titles = ('01  Four material regions', f'02  Baseline: directional plies {angle_label(baseline)}',
              f'03  Rotated: directional plies {angle_label(rotated)}')
    for index, axis in enumerate(axes):
        axis.set_facecolor('#f6f8fb')
        faces = []
        for eid in drawing_order:
            base_color = COLORS[elements[eid]['region']] if index == 0 else '#dce6ee'
            facing = abs(np.asarray(elements[eid]['normal']) @ camera[2])
            faces.append(np.clip(np.asarray(to_rgb(base_color)) * (.82 + .18 * facing), 0, 1))
        axis.add_collection(PolyCollection([projected[eid][:, :2] * 1000 for eid in drawing_order],
                                           facecolors=faces, edgecolors='#62788a',
                                           linewidths=.16, antialiased=True, zorder=1))
        axis.set_xlim(low[0] - margin[0], high[0] + margin[0])
        axis.set_ylim(low[1] - margin[1], high[1] + margin[1])
        axis.set_aspect('equal')
        axis.axis('off')
        axis.set_title(titles[index], loc='left', fontsize=15, weight='bold', pad=16)
        if index:
            data, assignment = (elements, baseline) if index == 1 else (rotated_elements, rotated)
            for eid in sampled:
                record = data[eid]
                laminate = assignment['laminates'][record['region']]
                if material_model(laminate) == 'plane_isotropic':
                    continue
                angle = laminate['angles_deg'][0]
                center = np.asarray(record['centroid_m']) @ camera.T * 1000
                for direction, color in zip(fabric_directions(record, angle), (WARP, FILL)):
                    vector = direction @ camera.T * 55.0
                    start, end = center[:2] - vector[:2] / 2, center[:2] + vector[:2] / 2
                    axis.annotate('', xy=end, xytext=start,
                                  arrowprops={'arrowstyle': '-|>', 'color': color, 'lw': 1.35,
                                              'mutation_scale': 6, 'shrinkA': 0, 'shrinkB': 0}, zorder=3)
        # A physical scale applies to this orthographic projection, without stretching.
        bar_x, bar_y = low[0], low[1] - .075 * max(extent)
        axis.plot([bar_x, bar_x + 100], [bar_y, bar_y], color='#52687b', lw=2)
        axis.text(bar_x + 50, bar_y - 16, '100 mm', fontsize=9, color='#64758a', ha='center', va='top')

    region_labels = {'upper': 'Upper', 'lower': 'Lower (mostly hidden)',
                     'leading_edge': 'Leading edge', 'trailing_edge': 'Trailing edge'}
    handles = []
    for region in REGIONS:
        laminate = baseline['laminates'][region]
        description = ('equivalent plane-isotropic'
                       if material_model(laminate) == 'plane_isotropic'
                       else f"{len(laminate['angles_deg'])} orthotropic plies")
        text = (f"{region_labels[region]}: {description}, "
                f"{laminate['total_thickness_m'] * 1000:.3f} mm")
        handles.append(Patch(facecolor=COLORS[region], edgecolor='none', label=text))
    fig.legend(handles=handles, loc='lower left', bbox_to_anchor=(.030, .18), ncol=2,
               frameon=False, fontsize=10.5, handlelength=1.4, columnspacing=1.7)
    fig.legend(handles=[Line2D([0], [0], color=WARP, lw=2, label='Material axis 1: θ'),
                        Line2D([0], [0], color=FILL, lw=2, label='Material axis 2: θ + 90°')],
               loc='lower right', bbox_to_anchor=(.974, .21), ncol=2, frameon=False, fontsize=11)
    fig.text(.035, .15, 'Arrows show orthotropic material axes only. Plane-isotropic regions have no preferred in-plane direction.',
             fontsize=11, weight='bold')
    fig.text(.035, .083,
             'θ rotates from each shell reference axis (projected outward span; chord on the tip cap). '
             'Arrows show the first ply only, at sparse visible locations.\n'
             'Lower skin is mostly occluded; assignment includes both surfaces. '
             'Arrow lengths are illustrative. Geometry uses true scale; ply thickness is listed, not drawn.',
             fontsize=10.5, color='#64758a', linespacing=1.65)
    fig.text(.035, .026,
             'Source: baseline / angle_30 cases/part_1/composite_assignment.json + structural_cards.bdf  '
             '|  Temporary assignment illustration, not a stress or strength result.',
             fontsize=9.2, color='#7b8b9e')
    output = trial_root / 'assignment_overview.png'
    fig.savefig(output, dpi=190, facecolor=fig.get_facecolor())
    plt.close(fig)
    print(json.dumps({'output': str(output), 'elements': len(elements),
                      'visible_arrow_pairs_per_panel': {
                          name: sum(elements[eid]['region'] in directional_regions(assignment)
                                    for eid in sampled)
                          for name, assignment in [('baseline', baseline), ('angle_30', rotated)]},
                      'material_models': {region: material_model(data)
                                          for region, data in baseline['laminates'].items()},
                      'first_ply_angles_deg': {name: {region: data['angles_deg'][0]
                                                    for region, data in assignment['laminates'].items()
                                                    if material_model(data) == 'orthotropic'}
                                              for name, assignment in [('baseline', baseline), ('angle_30', rotated)]}}, indent=2))
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--trial-root', type=Path, default=TRIAL_ROOT)
    args = parser.parse_args()
    render(args.trial_root)


if __name__ == '__main__':
    main()
