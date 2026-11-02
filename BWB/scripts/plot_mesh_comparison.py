"""Compare actual free and structured upper-skin BDF meshes in CAD Z/X.

Default view: part 1, central main-wing span Z = 0.20..0.65 m. The old
upper skin is selected by positive CAD-Y element normals on the main wing;
the new upper skin uses element_patches in its mesh_layout.json. Geometry
and element connectivity are read from the BDFs without reshaping them.
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

PROJECT = Path(__file__).resolve().parents[1]


def polygons_from_bdf(path: Path):
    model = load_model(path)
    if any(node['cp'] != 0 for node in model['nodes'].values()):
        raise ValueError(f'Nonzero GRID coordinate system in {path}')
    if any(model['duplicate_ids'].values()):
        raise ValueError(f'Duplicate entity IDs in {path}')
    return {eid: np.asarray([model['nodes'][nid]['xyz'] for nid in element['nodes']])
            for eid, element in model['elements'].items()}


def render(old_path: Path, new_path: Path, layout_path: Path, output: Path,
           span_min=.20, span_max=.65):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection

    if span_min >= span_max:
        raise ValueError('span_min must be smaller than span_max')
    old_all = polygons_from_bdf(old_path)
    new_all = polygons_from_bdf(new_path)
    layout = json.loads(layout_path.read_text(encoding='utf-8'))

    def is_old_upper(points):
        normal = np.cross(points[1] - points[0], points[2] - points[0])
        length = np.linalg.norm(normal)
        return (length > 0 and normal[1] / length > .5
                and -.00001 <= points[:, 2].mean() <= .78001)

    old = [points for points in old_all.values() if is_old_upper(points)]
    new = [new_all[int(eid)] for eid, patch in layout['element_patches'].items()
           if patch == 'main_upper']

    def overlaps_view(points):
        return points[:, 2].max() > span_min and points[:, 2].min() < span_max

    selected = [[points for points in polygons if overlaps_view(points)]
                for polygons in (old, new)]
    if not all(selected):
        raise ValueError('No upper-main-wing elements selected in the requested span window')
    # Shared limits and equal physical aspect make the two topology patterns comparable.
    chord_values = np.concatenate([points[:, 0] for polygons in selected for points in polygons])
    chord_low, chord_high = chord_values.min(), chord_values.max()
    chord_margin = .025 * (chord_high - chord_low)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 8.5), sharex=True, sharey=True,
                             facecolor='#f5f8fa')
    fig.suptitle('BWB | Upper-skin mesh pattern', x=.055, y=.96, ha='left',
                 fontsize=23, fontweight='bold', color='#203c49')
    fig.text(.055, .912, 'Part 1  ·  Same CAD coordinates and scale  ·  Central main-wing region',
             fontsize=12, color='#526975')

    for index, (axis, polygons, label, path) in enumerate(zip(
            axes, selected, ('Previous: free surface mesh', 'Updated: structured span / chord grid'),
            (old_path, new_path))):
        projected = [points[:, [2, 0]] * 1000 for points in polygons]
        axis.add_collection(PolyCollection(projected, facecolors='#a9cfce', edgecolors='#315c68',
                                           linewidths=.65, antialiased=True))
        axis.set_facecolor('#ffffff')
        axis.set_xlim(span_min * 1000, span_max * 1000)
        axis.set_ylim((chord_low - chord_margin) * 1000, (chord_high + chord_margin) * 1000)
        axis.set_aspect('equal', adjustable='box')
        axis.set_title(label, loc='left', fontsize=13, fontweight='bold', color='#203c49', pad=14)
        axis.set_xlabel('CAD Z / span (mm)', labelpad=10, color='#526975')
        axis.tick_params(colors='#607783', labelsize=9)
        for spine in axis.spines.values():
            spine.set_color('#cedae0')
        try:
            source_label = path.resolve().relative_to(PROJECT.parent).as_posix()
        except ValueError:
            source_label = str(path.resolve())
        fig.text(.055 if index == 0 else .535, .112, f'Source: {source_label}',
                 fontsize=8.1, color='#526975')
    axes[0].set_ylabel('CAD X / chord (mm)', labelpad=10, color='#526975')
    fig.text(.055, .065,
             'Actual exported shell elements; orthographic Z/X projection. No geometry or connectivity was altered for this plot.',
             fontsize=10, color='#526975')
    fig.text(.055, .035,
             'Selection: previous mesh = main-wing normals toward +Y; updated mesh = mesh_layout.json main_upper patch.',
             fontsize=9, color='#728691')
    fig.subplots_adjust(left=.065, right=.975, top=.855, bottom=.195, wspace=.13)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=190, facecolor=fig.get_facecolor())
    plt.close(fig)
    return {'output': str(output.resolve()), 'selected_elements_old': len(selected[0]),
            'selected_elements_new': len(selected[1]), 'span_window_m': [span_min, span_max]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--old', type=Path, default=PROJECT / 'models/skin_demo/part_1/part_1_skin.bdf')
    parser.add_argument('--new', type=Path, default=PROJECT / 'models/skin_structured/part_1/part_1_skin.bdf')
    parser.add_argument('--layout', type=Path, default=PROJECT / 'models/skin_structured/part_1/mesh_layout.json')
    parser.add_argument('--output', type=Path, default=PROJECT / 'models/skin_structured/mesh_pattern_comparison.png')
    parser.add_argument('--span-min', type=float, default=.20)
    parser.add_argument('--span-max', type=float, default=.65)
    args = parser.parse_args()
    print(json.dumps(render(args.old, args.new, args.layout, args.output, args.span_min, args.span_max), indent=2))


if __name__ == '__main__':
    main()
