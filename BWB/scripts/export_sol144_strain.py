"""Export verified SOL144 GRID displacements and shell fiber strains for BWB.

Request DISPLACEMENT(PRINT,SORT1,REAL)=ALL and
STRAIN(PRINT,SORT1,REAL,FIBER,CENTER)=ALL. No solver is started by this script.
NASA_SHM's fiber-table reader is imported explicitly, read-only; its export
workflow and output-directory policy are not used.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import sys
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from NASA_SHM.core.plate_strain import read_fiber_strains

try:
    from .inspect_shell_bdf import cards, load_model, nastran_float
except ImportError:
    from inspect_shell_bdf import cards, load_model, nastran_float


NASA_READER = ROOT / "NASA_SHM/core/plate_strain.py"
END_MARKER = "* * * END OF JOB * * *"
SUBCASE_HEADER = re.compile(r"^0\s+.*?\bSUBCASE\s+(\d+)\s*$")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def principal_strains(exx, eyy, gamma_xy):
    """In-plane strain tensor eigenvalues; gamma_xy is engineering shear."""
    require(all(math.isfinite(v) for v in (exx, eyy, gamma_xy)), "Nonfinite strain component")
    average = (exx + eyy) / 2
    radius = math.hypot((exx - eyy) / 2, gamma_xy / 2)
    angle = math.degrees(0.5 * math.atan2(gamma_xy, exx - eyy))
    return average + radius, average - radius, angle


def read_displacements(path):
    """Read real SORT1 GRID translations/rotations with explicit subcase IDs."""
    rows, subcase, active, ended = {}, None, False, False
    with Path(path).open(errors="replace") as stream:
        for line in stream:
            require("FATAL" not in line.upper(), f"Fatal diagnostic in F06: {line.strip()}")
            ended |= END_MARKER in line
            match = SUBCASE_HEADER.search(line)
            if match:
                subcase = int(match[1])
            compact = "".join(line.split()).upper()
            if "DISPLACEMENTVECTOR" in compact:
                active = True
            elif any(title in compact for title in (
                    "STRAINSIN", "STRESSESIN", "FORCESOFSINGLE", "FORCESOFMULTI",
                    "ANALYSISSUMMARY", "AERODYNAMICMONITOR", "GRIDPOINTFORCEBALANCE")):
                active = False
            if not active:
                continue
            fields = line.split()
            if len(fields) == 9 and fields[0] == "0":
                fields = fields[1:]
            if len(fields) == 8 and fields[0].isdigit() and fields[1] == "G":
                require(subcase is not None, "Displacement row has no subcase header")
                values = tuple(nastran_float(v) for v in fields[2:])
                require(all(math.isfinite(v) for v in values), "Nonfinite displacement")
                key = (subcase, int(fields[0]))
                require(key not in rows, f"Duplicate displacement row {key}")
                rows[key] = values
    require(ended, "F06 has no END OF JOB marker")
    require(rows, "F06 has no real GRID displacement table")
    return rows


def read_bwb_fiber_strains(path):
    """Adapt labeled page headers in a temporary copy for NASA's pure reader."""
    normalized, changed = [], 0
    with Path(path).open(errors="replace") as stream:
        for line in stream:
            match = SUBCASE_HEADER.search(line)
            if match:
                replacement = f"0 SUBCASE {match[1]}\n"
                changed += replacement != line
                line = replacement
            normalized.append(line)
    text = "".join(normalized)
    with tempfile.TemporaryDirectory(prefix="bwb_f06_reader_") as temporary:
        adapted = Path(temporary) / "reader_input.f06"
        adapted.write_text(text, encoding="utf-8")
        rows = read_fiber_strains(adapted)
    return rows, dict(normalized_subcase_headers=changed,
                      normalized_text_sha256=hashlib.sha256(text.encode("utf-8")).hexdigest(),
                      original_f06_modified=False,
                      procedure="Only labeled/page SUBCASE header lines normalized in a temporary parser input.")


def _expected_subcases(deck):
    """Read explicit case-control subcases, including relative INCLUDE files."""
    found, visited = set(), set()

    def visit(path):
        path = path.resolve()
        if path in visited:
            return
        visited.add(path)
        text = path.read_text(errors="replace")
        found.update(int(value) for value in re.findall(r"^\s*SUBCASE\s+(\d+)\s*$", text, re.M | re.I))
        for value in re.findall(r"^\s*INCLUDE\s+['\"]([^'\"]+)['\"]", text, re.M | re.I):
            child = path.parent / value
            require(child.is_file(), f"Missing INCLUDE input: {child}")
            visit(child)

    if deck.is_file():
        visit(deck)
    return found


def validate_case(case_dir, job_name="bwb_sol144", root_tolerance=1e-10):
    case_dir = Path(case_dir).resolve()
    require(Path(job_name).name == job_name, "Job name must not contain a directory")
    require(math.isfinite(root_tolerance) and root_tolerance >= 0, "Invalid root tolerance")
    f06, structural_path = case_dir / f"{job_name}.f06", case_dir / "structural_cards.bdf"
    require(f06.is_file() and structural_path.is_file(), "Need F06 and structural_cards.bdf")
    model = load_model(structural_path)
    require(model["nodes"] and model["elements"], "Structural model contains no GRID/shell data")
    require(not any(model["duplicate_ids"].values()), "Duplicate structural IDs")
    require(all(n["cp"] == n["cd"] == 0 for n in model["nodes"].values()),
            "This exporter requires basic GRID input/output coordinate systems")
    displacement = read_displacements(f06)
    strain, parser_adaptation = read_bwb_fiber_strains(f06)
    subcases = sorted({s for s, _ in displacement})
    requested = _expected_subcases(case_dir / f"{job_name}.bdf")
    require(not requested or requested == set(subcases), "F06 subcases do not match requested case control")
    expected_grid = {(s, n) for s in subcases for n in model["nodes"]}
    expected_shell = {(s, e) for s in subcases for e in model["elements"]}
    extra_displacements = {key: value for key, value in displacement.items() if key not in expected_grid}
    generated_aero_count = 0
    if extra_displacements:
        # Simcenter's DISPLACEMENT=ALL may include internally generated CAERO1
        # box centers and drawing corner points. These are NOT basic-coordinate
        # structural GRID results. Validate this specific output pattern and
        # preserve the extra table separately instead of silently discarding it.
        aero_path = case_dir / "aero_cards.bdf"
        require(aero_path.is_file(), "Extra GRID results without aerodynamic input provenance")
        for name, fields, _ in cards(aero_path):
            if name == "CAERO1":
                ns, nc = int(fields[4] or 0), int(fields[5] or 0)
                require(ns > 0 and nc > 0, "Extra aerodynamic GRID validation requires uniform CAERO1 divisions")
                generated_aero_count += ns * nc + (ns + 1) * (nc + 1)
        require(generated_aero_count > 0, "Extra GRID results without CAERO1 geometry")
        first_extra = max(model["nodes"]) + 1
        expected_extra = {(s, n) for s in subcases
                          for n in range(first_extra, first_extra + generated_aero_count)}
        require(set(extra_displacements) == expected_extra,
                "Unexpected nonstructural GRID IDs/count; cannot classify as generated CAERO1 output")
        require(all(abs(values[index]) <= root_tolerance for values in extra_displacements.values()
                    for index in (0, 1, 3, 5)),
                "Unexpected DOFs on generated aerodynamic GRID output")
        displacement = {key: value for key, value in displacement.items() if key in expected_grid}
    require(set(displacement) == expected_grid, "Missing/extra GRID displacement rows or subcases")
    require(set(strain) == expected_shell, "Missing/extra shell fiber-strain rows or subcases")
    for (subcase, eid), row in strain.items():
        shell = model["elements"][eid]
        require(row["type"] == shell["type"], f"Shell result type mismatch: {subcase}/{eid}")
        require(shell["pid"] in model["properties"], f"Missing PSHELL for element {eid}")
        prop = model["properties"][shell["pid"]]
        thickness = prop["thickness"]
        require(thickness is not None and thickness > 0, f"Nonpositive shell thickness: {eid}")
        z1 = -thickness / 2 if prop["z1"] is None else prop["z1"]
        z2 = thickness / 2 if prop["z2"] is None else prop["z2"]
        for actual, expected in zip((row["z_bottom"], row["z_top"]), sorted((z1, z2))):
            require(math.isclose(actual, expected, rel_tol=2e-6, abs_tol=1e-10),
                    f"Fiber distance differs from PSHELL: {subcase}/{eid}")
        for surface in ("top", "bottom", "mid"):
            require(all(math.isfinite(v) for v in row[surface]), f"Nonfinite strain: {subcase}/{eid}/{surface}")
    root_constraints = model["constraints"].get(1, {})
    require(root_constraints, "No root constraints in SPC set 1")
    require(all(set(dofs) == set(range(1, 7)) and all(value == 0 for value in dofs.values())
                for dofs in root_constraints.values()), "SPC set 1 must fully clamp each root GRID")
    root_max = 0.0
    for subcase in subcases:
        for nid, dofs in root_constraints.items():
            require(nid in model["nodes"], f"SPC references missing GRID {nid}")
            for dof, expected in dofs.items():
                difference = abs(displacement[(subcase, nid)][dof - 1] - expected)
                root_max = max(root_max, difference)
                require(difference <= root_tolerance,
                        f"Fixed-root displacement exceeds tolerance: subcase {subcase}, GRID {nid}, DOF {dof}")
    input_hashes = {path.name: file_hash(path) for path in sorted(case_dir.iterdir())
                    if path.suffix.lower() in (".bdf", ".inc")}
    recorded_path = case_dir / "run-summary.json"
    recorded = json.loads(recorded_path.read_text(encoding="utf-8")) if recorded_path.is_file() else {}
    require(recorded.get("fresh_f06") is not False, "Recorded solver output is stale")
    if recorded.get("input_sha256"):
        require(recorded["input_sha256"] == input_hashes, "Solved input hashes differ from current case inputs")
    diagnostics = {"warning_messages": [], "runtime_error_lines": []}
    for path in (f06, case_dir / f"{job_name}.log"):
        if path.is_file():
            for line in path.read_text(errors="replace").splitlines():
                require("FATAL" not in line.upper(), f"Fatal diagnostic in {path.name}: {line.strip()}")
                if "WARNING MESSAGE" in line:
                    diagnostics["warning_messages"].append(line.strip())
                if path.suffix == ".log" and re.search(r"\bERROR\b", line, re.I):
                    diagnostics["runtime_error_lines"].append(line.strip())
    for key in diagnostics:
        diagnostics[key] = sorted(set(diagnostics[key]))
    diagnostics["recorded_runtime_errors"] = recorded.get("runtime_errors", [])
    return dict(model=model, displacement=displacement, strain=strain, subcases=subcases,
                f06=f06, structural_path=structural_path, input_sha256=input_hashes,
                root_nodes=sorted(root_constraints), root_max_abs_dof=root_max,
                root_tolerance=root_tolerance, recorded_run=recorded, diagnostics=diagnostics,
                parser_adaptation=parser_adaptation, extra_displacements=extra_displacements,
                generated_aero_count=generated_aero_count)


def _plot_strain(path, model, rows, subcase, diagnostic_only=False):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import TwoSlopeNorm
    from matplotlib.ticker import FormatStrFormatter, MaxNLocator
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    ids = sorted(model["elements"])
    polygons = [np.array([model["nodes"][n]["xyz"] for n in model["elements"][eid]["nodes"]]) for eid in ids]
    xyz = np.array([n["xyz"] for n in model["nodes"].values()])
    low, high = xyz.min(axis=0), xyz.max(axis=0)
    dimensions = np.maximum(high - low, 1e-6)
    by_key = {(r["element_id"], r["surface"]): r for r in rows if r["subcase"] == subcase}
    magnitude = max(abs(r["dominant_principal_microstrain"]) for r in by_key.values())
    magnitude = max(magnitude, 1e-9)
    norm, cmap = TwoSlopeNorm(vmin=-magnitude, vcenter=0, vmax=magnitude), plt.get_cmap("coolwarm")
    fig = plt.figure(figsize=(15, 6.5), facecolor="white")
    axes = []
    for index, (surface, label) in enumerate((("top", "+element-z fiber"),
                                             ("bottom", "-element-z fiber"),
                                             ("mid", "Mid: component mean")), 1):
        ax = fig.add_subplot(1, 3, index, projection="3d")
        axes.append(ax)
        values = [by_key[(eid, surface)]["dominant_principal_microstrain"] for eid in ids]
        collection = Poly3DCollection(polygons, facecolors=cmap(norm(values)), linewidths=0.10,
                                      edgecolors=(0.1, 0.1, 0.1, 0.25), rasterized=True)
        ax.add_collection3d(collection)
        ax.set(xlim=(low[0], high[0]), ylim=(low[1], high[1]), zlim=(low[2], high[2]),
               xlabel="CAD X [m]", ylabel="CAD Y [m]", zlabel="CAD Z [m]", title=label)
        ax.set_box_aspect(dimensions)
        ax.view_init(elev=24, azim=-56)
        ax.tick_params(labelsize=7, pad=0)
        ax.xaxis.set_major_locator(MaxNLocator(5))
        ax.yaxis.set_major_locator(MaxNLocator(3))
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
        ax.zaxis.set_major_locator(MaxNLocator(6))
        ax.grid(False)
    title = f"BWB SOL144 | Subcase {subcase} | shell centroid strain"
    if diagnostic_only:
        title += "\nDIAGNOSTIC RUN: runtime errors reported"
    fig.suptitle(title, fontsize=13 if diagnostic_only else 15, y=.98)
    fig.text(.5, .08, "Signed principal strain with largest magnitude; original geometry; no nodal averaging.",
             ha="center", fontsize=10)
    fig.subplots_adjust(left=.02, right=.88, top=.9, bottom=.15, wspace=.07)
    colorbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=fig.add_axes([.91, .23, .015, .53]))
    colorbar.set_label("Principal strain [microstrain]")
    fig.savefig(path, dpi=200, facecolor="white")
    plt.close(fig)


def export_case(case_dir, output_dir, job_name="bwb_sol144", root_tolerance=1e-10, make_plot=True):
    validated = validate_case(case_dir, job_name, root_tolerance)
    model, strain, displacement = (validated[key] for key in ("model", "strain", "displacement"))
    diagnostic_only = bool(validated["diagnostics"]["runtime_error_lines"] or
                           validated["diagnostics"]["recorded_runtime_errors"])
    output_dir = Path(output_dir).resolve()
    rows = []
    for (subcase, eid), fibers in sorted(strain.items()):
        shell = model["elements"][eid]
        centroid = np.mean([model["nodes"][n]["xyz"] for n in shell["nodes"]], axis=0)
        for surface in ("top", "bottom", "mid"):
            exx, eyy, gamma = fibers[surface]
            first, second, angle = principal_strains(exx, eyy, gamma)
            z = (fibers["z_top"] if surface == "top" else fibers["z_bottom"] if surface == "bottom"
                 else (fibers["z_top"] + fibers["z_bottom"]) / 2)
            rows.append(dict(subcase=subcase, element_id=eid, element_type=shell["type"],
                property_id=shell["pid"], surface=surface, fiber_z_m=z,
                centroid_x_m=float(centroid[0]), centroid_y_m=float(centroid[1]), centroid_z_m=float(centroid[2]),
                exx_element=exx, eyy_element=eyy, gamma_xy_element=gamma,
                exx_microstrain=exx * 1e6, eyy_microstrain=eyy * 1e6, gamma_xy_microstrain=gamma * 1e6,
                principal_max_microstrain=first * 1e6, principal_min_microstrain=second * 1e6,
                principal_angle_deg=angle,
                dominant_principal_microstrain=(first if abs(first) >= abs(second) else second) * 1e6))
    output_dir.mkdir(parents=True, exist_ok=True)
    strain_path = output_dir / "shell_strain_long.csv"
    with strain_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    displacement_path = output_dir / "grid_displacements.csv"
    displacement_rows = []
    for (subcase, nid), values in sorted(displacement.items()):
        xyz = model["nodes"][nid]["xyz"]
        displacement_rows.append(dict(subcase=subcase, grid_id=nid, x_m=xyz[0], y_m=xyz[1], z_m=xyz[2],
            ux_m=values[0], uy_m=values[1], uz_m=values[2], rx_rad=values[3], ry_rad=values[4], rz_rad=values[5],
            translation_m=math.sqrt(sum(v * v for v in values[:3]))))
    with displacement_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(displacement_rows[0]))
        writer.writeheader()
        writer.writerows(displacement_rows)
    extra_path = None
    if validated["extra_displacements"]:
        extra_path = output_dir / "solver_generated_aero_displacements.csv"
        with extra_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(["subcase", "solver_output_point_id", "t1_local_m", "t2_local_m", "t3_local_m",
                             "r1_local_rad", "r2_local_rad", "r3_local_rad"])
            writer.writerows([*key, *values] for key, values in sorted(validated["extra_displacements"].items()))
    cases = []
    for subcase in validated["subcases"]:
        case_rows = [row for row in rows if row["subcase"] == subcase]
        surface_summary = {}
        for surface in ("top", "bottom", "mid"):
            selected = [row for row in case_rows if row["surface"] == surface]
            critical = max(selected, key=lambda row: abs(row["dominant_principal_microstrain"]))
            surface_summary[surface] = dict(
                max_principal_microstrain=max(r["principal_max_microstrain"] for r in selected),
                min_principal_microstrain=min(r["principal_min_microstrain"] for r in selected),
                max_abs_principal_microstrain=abs(critical["dominant_principal_microstrain"]),
                max_abs_element_id=critical["element_id"])
        largest_displacement = max((r for r in displacement_rows if r["subcase"] == subcase),
                                   key=lambda row: row["translation_m"])
        item = dict(subcase=subcase, surfaces=surface_summary,
                    max_translation_m=largest_displacement["translation_m"],
                    max_translation_grid_id=largest_displacement["grid_id"])
        if make_plot:
            image_name = "strain_overview.png" if len(validated["subcases"]) == 1 else f"strain_subcase_{subcase}.png"
            _plot_strain(output_dir / image_name, model, rows, subcase, diagnostic_only)
            item["plot"] = image_name
        cases.append(item)
    source_bytes = NASA_READER.read_bytes()
    summary = dict(status="exported_with_runtime_diagnostics" if diagnostic_only else "exported",
        diagnostic_only=diagnostic_only, solver_end_of_job=True, fatal_messages=[],
        job_name=job_name, case_dir=str(Path(case_dir).resolve()),
        source_f06=str(validated["f06"]), source_f06_sha256=file_hash(validated["f06"]),
        f06_parser_adaptation=validated["parser_adaptation"],
        input_sha256=validated["input_sha256"],
        nodes_per_subcase=len(model["nodes"]), shells_per_subcase=len(model["elements"]),
        solver_generated_aero_points_per_subcase=validated["generated_aero_count"],
        nonstructural_grid_policy="CAERO1 box-center plus corner-point count and local DOF pattern checked; preserved separately in local aerodynamic coordinates, not treated as structural GRID data.",
        strain_csv_rows=len(rows), displacement_csv_rows=len(displacement_rows), subcases=cases,
        root_constraint_set=1, root_fixed_nodes=len(validated["root_nodes"]),
        max_abs_root_dof=validated["root_max_abs_dof"], root_tolerance=root_tolerance,
        coordinates="Displacements: basic CAD coordinates. Strains: local element axes, centroid values.",
        strain_definition="gamma_xy is engineering shear; principal strains use gamma_xy/2 in the tensor.",
        surfaces={"top": "+element-z fiber", "bottom": "-element-z fiber",
                  "mid": "Arithmetic mean of top/bottom components; principal strains recomputed from that mean."},
        units="m-N-kg-s; strains dimensionless and microstrain; rotations radians.",
        precision="Nastran F06 printed precision; no interpolation or nodal averaging of results.",
        nasa_reuse=dict(path=NASA_READER.relative_to(ROOT).as_posix(), symbol="read_fiber_strains",
                       mode="explicit read-only Python import; no NASA output workflow called",
                       sha256=hashlib.sha256(source_bytes.replace(b"\r\n", b"\n")).hexdigest(),
                       hash_mode="crlf_to_lf"),
        recorded_output_freshness=validated["recorded_run"].get("fresh_f06"),
        recorded_launcher_returncode=validated["recorded_run"].get("launcher_returncode"),
        diagnostics=validated["diagnostics"],
        files={"strain": strain_path.name, "displacements": displacement_path.name},
        csv_sha256={strain_path.name: file_hash(strain_path), displacement_path.name: file_hash(displacement_path)})
    if extra_path:
        summary["files"]["solver_generated_aero_displacements"] = extra_path.name
        summary["csv_sha256"][extra_path.name] = file_hash(extra_path)
    (output_dir / "strain_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--job-name", default="bwb_sol144")
    parser.add_argument("--root-tolerance", type=float, default=1e-10,
                        help="Absolute fixed DOF tolerance in m for translations and rad for rotations")
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args()
    try:
        result = export_case(args.case_dir, args.output_dir, args.job_name, args.root_tolerance, not args.no_plot)
    except (ValueError, OSError, KeyError) as exc:
        parser.exit(1, f"Result export rejected: {exc}\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
