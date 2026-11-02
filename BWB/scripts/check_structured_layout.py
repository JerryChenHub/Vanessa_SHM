"""Independently check exported BDF coordinates against structured grid metadata.

No Femap session or CAD projection is used. These checks establish layout and
connectivity, not solver accuracy or continuous CAD-surface approximation error.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path

import numpy as np

from inspect_shell_bdf import load_model


SKIN_PATCHES = tuple(f'{zone}_{side}' for zone in
                     ('main', 'transition', 'bend', 'winglet')
                     for side in ('lower', 'upper'))


def summary(values):
    a = np.asarray(values, dtype=float)
    return dict(count=len(a), min=float(a.min()), median=float(np.median(a)),
                p95=float(np.percentile(a, 95)), max=float(a.max()))


def quality(element_ids, elements, xyz):
    ratios, angles, warp, areas = [], [], [], []
    reversed_corners, nonpositive_area = [], []
    for eid in element_ids:
        p = np.asarray([xyz[n] for n in elements[eid]['nodes']])
        edge = np.roll(p, -1, axis=0) - p
        lengths = np.linalg.norm(edge, axis=1)
        ratios.append(float(lengths.max() / lengths.min()))
        corner_normals = np.cross(edge, -np.roll(edge, 1, axis=0))
        reference = corner_normals.sum(axis=0)
        if np.any(corner_normals @ reference <= 0):
            reversed_corners.append(eid)
        cosines = np.sum(edge * -np.roll(edge, 1, axis=0), axis=1)
        cosines /= lengths * np.roll(lengths, 1)
        angles.extend(np.degrees(np.arccos(np.clip(cosines, -1, 1))).tolist())
        n1 = np.cross(p[1] - p[0], p[2] - p[0])
        area = np.linalg.norm(n1) / 2
        if len(p) == 4:
            n2 = np.cross(p[2] - p[0], p[3] - p[0])
            area += np.linalg.norm(n2) / 2
            cosine = np.dot(n1, n2) / (np.linalg.norm(n1) * np.linalg.norm(n2))
            warp.append(float(np.degrees(np.arccos(np.clip(cosine, -1, 1)))))
        areas.append(float(area))
        if area <= 0:
            nonpositive_area.append(eid)
    return dict(element_count=len(element_ids),
                element_types=dict(Counter(elements[e]['type'] for e in element_ids)),
                edge_length_ratio=summary(ratios), corner_angle_deg=summary(angles),
                diagonal_warp_deg=summary(warp) if warp else None,
                area_m2=summary(areas), total_area_m2=sum(areas),
                reversed_or_flat_corner_element_ids=reversed_corners,
                nonpositive_area_element_ids=nonpositive_area)


def check_layout(folder: Path):
    candidates = list(folder.glob('part_*_skin.bdf'))
    if len(candidates) != 1:
        raise ValueError(f'Expected one skin BDF in {folder}')
    path = candidates[0]
    model = load_model(path)
    layout = json.loads((folder / 'mesh_layout.json').read_text(encoding='utf-8'))
    elements = model['elements']
    xyz = {n: np.asarray(v['xyz']) for n, v in model['nodes'].items()}
    patch_of = {int(e): patch for e, patch in layout['element_patches'].items()}
    assert set(patch_of) == set(elements), 'Patch labels do not cover exported elements'
    assert all(v['cp'] == 0 for v in model['nodes'].values()), 'Non-basic GRID coordinates'
    assert not any(model['duplicate_ids'].values()), 'Duplicate BDF entity IDs'
    neighbours = defaultdict(set)
    attached_elements = defaultdict(set)
    edge_uses = defaultdict(list)
    cells = defaultdict(list)
    for eid, e in elements.items():
        conn = e['nodes']
        assert len(set(conn)) == len(conn), f'Repeated corner in element {eid}'
        cells[frozenset(conn)].append(eid)
        for a, b in zip(conn, conn[1:] + conn[:1]):
            neighbours[a].add(b)
            neighbours[b].add(a)
            edge_uses[tuple(sorted((a, b)))].append((eid, a, b))
        for n in conn:
            attached_elements[n].add(eid)
    assert all(len(v) == 1 for v in cells.values()), 'Overlapping duplicate elements'
    assert all(len(v) <= 2 for v in edge_uses.values()), 'Nonmanifold edge'
    assert all(v[0][1:] != v[1][1:] for v in edge_uses.values() if len(v) == 2), 'Inconsistent normals'

    patch_reports = {}
    matched = set()
    grids = {}
    # Small-field BDF export rounds metre coordinates at roughly micrometre scale.
    fraction_tolerance = 1e-5
    coordinate_tolerance = 2e-6
    for name in SKIN_PATCHES:
        grid = np.asarray(layout['patch_node_ids'][name], dtype=int)
        grids[name] = grid
        assert len(set(grid.ravel())) == grid.size, f'{name}: collapsed grid nodes'
        actual = []
        for j in range(grid.shape[0] - 1):
            for i in range(grid.shape[1] - 1):
                corners = frozenset(grid[j:j + 2, i:i + 2].ravel())
                assert len(cells[corners]) == 1, f'{name}: missing or duplicated cell {j},{i}'
                eid = cells[corners][0]
                assert patch_of[eid] == name and elements[eid]['type'] == 'CQUAD4'
                actual.append(eid)
        expected = {e for e, p in patch_of.items() if p == name}
        assert set(actual) == expected, f'{name}: extra elements outside logical cells'
        matched.update(actual)
        inside = grid[1:-1, 1:-1].ravel()
        assert all(len(neighbours[n]) == 4 for n in inside), f'{name}: non-four-valent interior'
        assert all(len(attached_elements[n]) == 4 for n in inside), f'{name}: interior incidence not four'
        p = np.asarray([[xyz[n] for n in row] for row in grid])
        report = dict(rows=grid.shape[0], columns=grid.shape[1], quad_count=len(actual),
                      interior_nodes_checked=len(inside), interior_edge_valence=4,
                      interior_element_incidence=4)
        if name.startswith('main_'):
            row_spread = float(np.ptp(p[:, :, 2], axis=1).max())
            width = p[:, -1, 0] - p[:, 0, 0]
            fractions = (p[:, :, 0] - p[:, 0, None, 0]) / width[:, None]
            fraction_error = float(np.max(abs(fractions - np.linspace(0, 1, grid.shape[1]))))
            spacing = np.diff(p[:, 0, 2])
            span_error = float(np.max(abs(spacing - spacing.mean())))
            assert np.all(width > 0) and np.all(np.diff(p[:, :, 0], axis=1) > 0)
            assert row_spread <= coordinate_tolerance, f'{name}: chord rows leave common Z plane'
            assert fraction_error <= fraction_tolerance, f'{name}: chord fractions are not uniform'
            assert span_error <= coordinate_tolerance, f'{name}: span stations not regular'
            report.update(max_row_z_spread_m=row_spread,
                          max_normalized_chord_fraction_error=fraction_error,
                          mean_signed_span_step_m=float(spacing.mean()),
                          max_span_step_error_m=span_error)
        patch_reports[name] = report

    for side in ('lower', 'upper'):
        for previous, following in zip(('main', 'transition', 'bend'), ('transition', 'bend', 'winglet')):
            assert np.array_equal(grids[f'{previous}_{side}'][-1], grids[f'{following}_{side}'][0]), 'Interface not conforming'
    for zone in ('main', 'transition', 'bend', 'winglet'):
        assert np.array_equal(grids[f'{zone}_lower'][:, 0], grids[f'{zone}_upper'][:, 0]), 'Leading edge not conforming'

    root_nodes = set(grids['main_lower'][0]) | set(grids['main_upper'][0])
    free_edges = [edge for edge, uses in edge_uses.items() if len(uses) == 1]
    assert all(a in root_nodes and b in root_nodes for a, b in free_edges), 'Unintended opening outside root'
    boundary_neighbours = defaultdict(set)
    for a, b in free_edges:
        boundary_neighbours[a].add(b)
        boundary_neighbours[b].add(a)
    assert set(boundary_neighbours) == root_nodes and all(len(v) == 2 for v in boundary_neighbours.values())
    visited = {next(iter(xyz))}
    pending = list(visited)
    while pending:
        for n in neighbours[pending.pop()] - visited:
            visited.add(n)
            pending.append(n)
    assert visited == set(xyz), 'Disconnected nodes or shell pieces'
    euler = len(xyz) - len(edge_uses) + len(elements)
    assert euler == 1, f'Expected disk topology (one open root), got Euler characteristic {euler}'
    skin_quality = quality(sorted(matched), elements, xyz)
    assert not skin_quality['reversed_or_flat_corner_element_ids'], 'Invalid skin quad corner geometry'
    assert not skin_quality['nonpositive_area_element_ids'], 'Degenerate skin element'
    groups = dict(eight_surface_skin_patches=skin_quality)
    for label, selector in (
        ('thin_trailing_edge_closure', lambda p: p.endswith('_trailing_edge')),
        ('thin_tip_cap_closure', lambda p: p == 'tip_cap')):
        groups[label] = quality([e for e, p in patch_of.items() if selector(p)], elements, xyz)
    report = dict(status='PASS', bdf=str(path.resolve()),
                  checks='Exported BDF coordinates and connectivity independently compared with layout metadata',
                  tolerance=dict(coordinate_m=coordinate_tolerance, normalized_chord_fraction=fraction_tolerance),
                  node_count=len(xyz), element_count=len(elements), patches=patch_reports,
                  topology=dict(connected_components=1, euler_characteristic=euler,
                                free_edges=len(free_edges), free_edges_only_on_root=True,
                                conforming_patch_interfaces=True), quality_groups=groups,
                  limitations=[
                      'Quality metrics are geometric screens, not Nastran acceptance criteria or convergence verification.',
                      'CAD distance and projection validity are not independently rechecked without a CAD kernel.',
                      'Eight skin patches use CAD-projected interior nodes; trailing edge is an explicitly simplified ruled closure.',
                      'Straight-edged shell facets approximate curved CAD between nodes.',
                      'Parallel chord rows are constant-Z sections in the main wing; their three-dimensional lines follow airfoil curvature.'
                  ])
    output = folder / 'structured_layout_qa.json'
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folders', nargs='+', type=Path)
    args = parser.parse_args()
    for folder in args.folders:
        result = check_layout(folder)
        print(f"PASS {folder}: {result['node_count']} nodes, {result['element_count']} elements; "
              f"eight regular quad patches and {result['topology']['free_edges']} root-only free edges")
