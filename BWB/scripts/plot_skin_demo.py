"""Render two existing BDF shell meshes as separate models in original CAD coordinates.

Run from the repository root::

    python -m BWB.scripts.plot_skin_demo part_1.bdf part_2.bdf

Input coordinates are meters; only plotted coordinates are converted to mm.
No geometry is generated, moved, mirrored or assembled by this script.
"""
from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

try:
    from .inspect_shell_bdf import load_model
except ImportError:  # Also allow direct execution by absolute script path.
    from inspect_shell_bdf import load_model


DEFAULT_OUTPUT = Path(__file__).resolve().parents[1] / "models" / "skin_structured" / "skin_mesh_overview.png"


def full_clamp_ids(model, spc_set):
    constraints = {nid: dofs.copy() for nid, dofs in model["permanent"].items()}
    selected_sets = [spc_set] if spc_set is not None else sorted(model["constraints"])
    for sid in selected_sets:
        for nid, dofs in model["constraints"].get(sid, {}).items():
            constraints.setdefault(nid, {}).update(dofs)
    return [nid for nid in model["nodes"]
            if all(constraints.get(nid, {}).get(dof) == 0.0 for dof in range(1, 7))]


def render(first_path: Path, second_path: Path, output: Path, spc_set=1):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    paths = [Path(first_path), Path(second_path)]
    models = [load_model(path) for path in paths]
    for path, model in zip(paths, models):
        if not model["nodes"] or not model["elements"]:
            raise ValueError(f"No supported shell mesh in {path}")
        if any(node["cp"] != 0 for node in model["nodes"].values()):
            raise ValueError(f"Nonzero GRID CP in {path}; transform coordinates before plotting")
        if any(model["duplicate_ids"].values()):
            raise ValueError(f"Duplicate entity IDs in {path}; inspect the BDF before plotting")
        if not full_clamp_ids(model, spc_set):
            raise ValueError(f"No zero-displacement DOF 123456 root nodes found in {path} for SPC set {spc_set}")

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.labelcolor": "#374151", "text.color": "#172b3a",
                         "xtick.color": "#64748b", "ytick.color": "#64748b"})
    figure = plt.figure(figsize=(16, 8.1), facecolor="#f7f9fc")
    colors = ["#93bad8", "#9bc7c0"]
    figure.suptitle("BWB | Shell skin mesh", x=.055, y=.965, ha="left", fontsize=23, fontweight="bold")
    figure.text(.055, .917, "Two separate models  ·  Skin only  ·  No ribs", fontsize=12, color="#526374")

    def plotted_xyz(xyz):
        # Horizontal plot axis = span (CAD Z); vertical = height (CAD Y).
        # Values retain their original coordinate origin and signs.
        return (1000 * xyz[2], 1000 * xyz[0], 1000 * xyz[1])

    for index, (path, model, color) in enumerate(zip(paths, models, colors)):
        axis = figure.add_subplot(1, 2, index + 1, projection="3d")
        axis.set_facecolor("#f7f9fc")
        coordinates = {nid: plotted_xyz(node["xyz"]) for nid, node in model["nodes"].items()}
        polygons = []
        for eid, element in model["elements"].items():
            try:
                polygons.append([coordinates[nid] for nid in element["nodes"]])
            except KeyError as exc:
                raise ValueError(f"Element {eid} in {path} references absent GRID {exc.args[0]}") from exc
        skin = Poly3DCollection(polygons, facecolors=color, edgecolors="#344c5e",
                               linewidths=.20, alpha=.94, antialiased=True)
        axis.add_collection3d(skin)
        clamp_ids = full_clamp_ids(model, spc_set)
        clamp = [coordinates[nid] for nid in clamp_ids]
        axis.scatter(*zip(*clamp), c="#d53f45", s=14, depthshade=False, edgecolors="white", linewidths=.25)

        points = list(coordinates.values())
        low = [min(point[dimension] for point in points) for dimension in range(3)]
        high = [max(point[dimension] for point in points) for dimension in range(3)]
        sizes = [high[dimension] - low[dimension] for dimension in range(3)]
        largest = max(sizes)
        for dimension, setter in enumerate((axis.set_xlim, axis.set_ylim, axis.set_zlim)):
            margin = max(sizes[dimension] * .065, largest * .025)
            setter(low[dimension] - margin, high[dimension] + margin)
        # Independent centers; preserve physical aspect rather than stretching height.
        axis.set_box_aspect([max(size, largest * .06) for size in sizes])
        axis.view_init(elev=24, azim=-64 if index == 0 else -116)
        axis.set_proj_type("ortho")
        axis.set_xlabel("CAD Z / span (mm)", labelpad=10)
        axis.set_ylabel("CAD X / chord (mm)", labelpad=10)
        axis.set_zlabel("CAD Y / height (mm)", labelpad=8)
        for pane in (axis.xaxis.pane, axis.yaxis.pane, axis.zaxis.pane):
            pane.set_facecolor((1, 1, 1, .30))
            pane.set_edgecolor((.8, .84, .89, .5))
        for dimension_axis in (axis.xaxis, axis.yaxis, axis.zaxis):
            dimension_axis._axinfo["grid"].update(color=(.74, .79, .85, .35), linewidth=.6)
        axis.tick_params(labelsize=8, pad=1)
        types = Counter(element["type"] for element in model["elements"].values())
        counts = "  ·  ".join(f"{count:,} {kind}" for kind, count in sorted(types.items()))
        axis.set_title(f"PART {index + 1}\n{len(model['nodes']):,} nodes  ·  {counts}", loc="left",
                       y=.91, fontsize=12, linespacing=1.7, fontweight="normal")
        x = .06 if index == 0 else .54
        thicknesses = sorted({prop["thickness"] for prop in model["properties"].values()
                              if prop["thickness"] is not None})
        thickness_text = ", ".join(f"{1000 * value:g}" for value in thicknesses) or "unspecified"
        figure.text(x, .167, f"Shell thickness: {thickness_text} mm  ·  Clamped nodes: {len(clamp_ids):,}",
                    fontsize=10, color="#526374")
        figure.text(x, .137, f"Source: {path.parent.name}/{path.name}", fontsize=8.6, color="#718096")

    figure.legend(handles=[Patch(facecolor="#a7c5d9", edgecolor="#344c5e", label="Actual BDF shell elements"),
                           Line2D([0], [0], marker="o", linestyle="", color="#d53f45", markersize=6,
                                  label="Root clamp: DOF 1–6 = 0")],
                  loc="lower left", bbox_to_anchor=(.048, .058), frameon=False, ncol=2, fontsize=10)
    figure.text(.055, .027,
                "Original CAD coordinates; each view centered independently. Display units: mm (source data: m). These are not an assembly.",
                fontsize=9, color="#718096")
    figure.subplots_adjust(left=.025, right=.975, top=.88, bottom=.20, wspace=.06)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=190, facecolor=figure.get_facecolor())
    plt.close(figure)
    return output.resolve()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("part_1_bdf", type=Path)
    parser.add_argument("part_2_bdf", type=Path)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--spc-set", type=int, default=1, help="Direct SPC/SPC1 set used to mark the root nodes")
    args = parser.parse_args()
    try:
        output = render(args.part_1_bdf, args.part_2_bdf, args.output, args.spc_set)
    except (OSError, ValueError) as exc:
        parser.exit(2, f"Mesh preview failed: {exc}\n")
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
