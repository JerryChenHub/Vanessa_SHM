"""Read-only shell-mesh QA for free, small-field and large-field Nastran BDF.

Example (SI model, two separate roots)::

    python -m BWB.scripts.inspect_shell_bdf model.bdf --units m \
        --root-plane z=0 --root-plane z=1.56 --output qa.json

The geometric measures are screening metrics, not Nastran element-quality
criteria. Coordinates in a nonzero GRID CP are not transformed. INCLUDE files
are followed relative to the file containing the INCLUDE statement.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import itertools
import json
import math
from pathlib import Path
import re
import sys


def nastran_float(value: str, default: float = 0.0) -> float:
    text = value.strip().replace("D", "E").replace("d", "e")
    if not text:
        return default
    if "e" not in text.lower():
        text = re.sub(r"(?<=[\d.])([+-]\d+)$", r"e\1", text)
    number = float(text)
    if not math.isfinite(number):
        raise ValueError(f"Nonfinite real: {value!r}")
    return number


def cards(path: Path, ancestors: tuple[Path, ...] = ()):
    """Yield (name, logical fields, location), preserving blank fields."""
    path = path.resolve()
    if path in ancestors:
        raise ValueError(f"Cyclic INCLUDE: {path}")
    lines = path.read_text(encoding="utf-8-sig", errors="replace").splitlines()
    current = None
    continuation_id = ""
    index = 0
    while index < len(lines):
        line_number = index + 1
        line = lines[index].split("$", 1)[0].rstrip()
        index += 1
        if not line.strip():
            continue
        if line.lstrip().upper().startswith("INCLUDE"):
            if current is not None:
                yield current
                current = None
            include = line.strip()[7:].strip()
            quote = include[:1]
            if quote in ("'", '"'):
                while include.count(quote) < 2 and index < len(lines):
                    include += lines[index].split("$", 1)[0].strip()
                    index += 1
                match = re.fullmatch(r"(['\"])(.*?)\1\s*", include)
                if not match:
                    raise ValueError(f"Invalid INCLUDE at {path}:{line_number}")
                include = match[2]
            if not include:
                raise ValueError(f"Empty INCLUDE at {path}:{line_number}")
            yield from cards(path.parent / include, ancestors + (path,))
            continuation_id = ""
            continue
        free = "," in line[:8]
        first = line.split(",", 1)[0].strip() if free else line[:8].strip()
        is_continuation = (not first or first.startswith(("+", "*"))
                           or bool(continuation_id and first == continuation_id))
        if free:
            parts = [part.strip() for part in line.split(",")]
            has_continuation_field = len(parts) > 9 and (not parts[9] or parts[9].startswith(("+", "*")))
            next_id = parts[9] if has_continuation_field else ""
            fields = parts[1:9] if has_continuation_field else parts[1:]
            # Fixed-width continuation positions apply to free-field records too.
            if next_id or (index < len(lines) and lines[index].lstrip().startswith(("+", "*", ","))):
                fields += [""] * (8 - len(fields))
        else:
            width = 16 if first.endswith("*") or first.startswith("*") else 8
            fields = [line[start:start + width].strip() for start in range(8, 72, width)]
            next_id = line[72:80].strip()
        if is_continuation:
            if current is None:
                raise ValueError(f"Orphan continuation at {path}:{line_number}")
            current[1].extend(fields)
        else:
            if current is not None:
                yield current
            name = first.rstrip("*").upper()
            current = [name, [name] + fields, f"{path}:{line_number}"]
        continuation_id = next_id
    if current is not None:
        yield current


def expanded_ids(tokens: list[str]) -> list[int]:
    tokens = [value.upper() for value in tokens if value]
    result = []
    i = 0
    while i < len(tokens):
        start = int(tokens[i])
        if i + 2 < len(tokens) and tokens[i + 1] == "THRU":
            end = int(tokens[i + 2])
            i += 3
            step = 1
            if i + 1 < len(tokens) and tokens[i] == "BY":
                step = int(tokens[i + 1])
                i += 2
            if step <= 0 or end < start:
                raise ValueError("Invalid THRU/BY range")
            result.extend(range(start, end + 1, step))
        else:
            result.append(start)
            i += 1
    return result


def _components(text: str) -> list[int]:
    result = [int(char) for char in text if char.strip()]
    if any(value not in range(1, 7) for value in result):
        raise ValueError(f"Expected structural DOF digits 1..6, got {text!r}")
    return result


def load_model(path: Path) -> dict:
    model = {name: {} for name in ("nodes", "elements", "materials", "properties")}
    model.update(constraints=defaultdict(lambda: defaultdict(dict)), permanent={},
                 card_counts=Counter(), warnings=[], duplicate_ids=defaultdict(list))
    def put(kind, key, value):
        if key in model[kind]:
            model["duplicate_ids"][kind].append(key)
        model[kind][key] = value
    for name, fields, location in cards(path):
        model["card_counts"][name] += 1
        field = lambda i: fields[i] if i < len(fields) else ""
        try:
            if name == "GRID":
                nid = int(field(1))
                put("nodes", nid, {"xyz": tuple(nastran_float(field(i)) for i in (3, 4, 5)),
                                   "cp": int(field(2) or 0), "cd": int(field(6) or 0)})
                if field(7):
                    model["permanent"][nid] = {dof: 0.0 for dof in _components(field(7))}
            elif name in ("CQUAD4", "CTRIA3"):
                count = 4 if name == "CQUAD4" else 3
                put("elements", int(field(1)), {"type": name, "pid": int(field(2)),
                    "nodes": tuple(int(field(i)) for i in range(3, 3 + count))})
            elif name == "MAT1":
                put("materials", int(field(1)), {key: nastran_float(field(i)) if field(i) else None
                    for i, key in ((2, "E"), (3, "G"), (4, "nu"), (5, "rho"))})
            elif name == "PSHELL":
                values = {key: int(field(i)) if field(i) else None
                          for i, key in ((2, "mid1"), (4, "mid2"), (6, "mid3"), (11, "mid4"))}
                values.update({key: nastran_float(field(i)) if field(i) else None
                    for i, key in ((3, "thickness"), (5, "bending_inertia_ratio"),
                                   (7, "transverse_shear_ratio"), (8, "nsm"),
                                   (9, "z1"), (10, "z2"))})
                put("properties", int(field(1)), values)
            elif name == "SPC1":
                sid, dofs = int(field(1)), _components(field(2))
                for nid in expanded_ids(fields[3:]):
                    model["constraints"][sid][nid].update({dof: 0.0 for dof in dofs})
            elif name == "SPC":
                sid = int(field(1))
                for i in range(2, len(fields), 3):
                    if field(i):
                        nid = int(field(i))
                        model["constraints"][sid][nid].update(
                            {dof: nastran_float(field(i + 2)) for dof in _components(field(i + 1))})
        except (ValueError, IndexError) as exc:
            raise ValueError(f"{location}: {name}: {exc}") from exc
    if model["card_counts"]["SPCADD"]:
        model["warnings"].append("SPCADD is not expanded; select a direct SPC/SPC1 set or inspect the union explicitly.")
    if any(node["cp"] != 0 for node in model["nodes"].values()):
        model["warnings"].append("Nonzero GRID CP: raw coordinate triples are not transformed; all geometric metrics and root-plane matches may be invalid.")
    if any(model["duplicate_ids"].values()):
        model["warnings"].append("Duplicate IDs found; metrics use the final definition of each ID.")
    return model


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def norm(a):
    return math.sqrt(sum(value * value for value in a))


def angle(a, b):
    denominator = norm(a) * norm(b)
    if denominator == 0:
        return None
    cosine = sum(x * y for x, y in zip(a, b)) / denominator
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def stats(values):
    values = sorted(value for value in values if value is not None and math.isfinite(value))
    if not values:
        return {"count": 0, "min": None, "median": None, "p95": None, "max": None, "mean": None}
    def percentile(fraction):
        position = (len(values) - 1) * fraction
        lower = int(position)
        upper = min(lower + 1, len(values) - 1)
        return values[lower] + (values[upper] - values[lower]) * (position - lower)
    return {"count": len(values), "min": values[0], "median": percentile(.5),
            "p95": percentile(.95), "max": values[-1], "mean": sum(values) / len(values)}


def near_coincident(nodes, tolerance):
    if tolerance <= 0:
        return {"enabled": False, "tolerance": tolerance}
    cells = defaultdict(list)
    count, examples = 0, []
    for nid, node in nodes.items():
        xyz = node["xyz"]
        cell = tuple(math.floor(value / tolerance) for value in xyz)
        for delta in itertools.product((-1, 0, 1), repeat=3):
            other_cell = tuple(a + b for a, b in zip(cell, delta))
            for other_id in cells.get(other_cell, ()):
                distance = norm(sub(xyz, nodes[other_id]["xyz"]))
                if distance <= tolerance:
                    count += 1
                    if len(examples) < 100:
                        examples.append({"nodes": [other_id, nid], "distance": distance})
        cells[cell].append(nid)
    return {"enabled": True, "tolerance": tolerance, "pair_count": count,
            "pair_examples": examples, "examples_truncated": count > len(examples)}


def inspect(path, units="unspecified", coincident_tol=1e-7, root_planes=(), root_tol=1e-6, spc_set=None):
    model = load_model(Path(path))
    nodes, elements = model["nodes"], model["elements"]
    xyz = [node["xyz"] for node in nodes.values()]
    bounds = {"min": [min(p[i] for p in xyz) for i in range(3)],
              "max": [max(p[i] for p in xyz) for i in range(3)]} if xyz else None
    if bounds:
        bounds["size"] = sub(bounds["max"], bounds["min"])
    parent = {nid: nid for nid in nodes}
    def find(nid):
        while parent[nid] != nid:
            parent[nid] = parent[parent[nid]]
            nid = parent[nid]
        return nid
    def union(a, b):
        parent[find(b)] = find(a)
    edges = defaultdict(list)
    used_nodes, valid_elements = set(), {}
    missing_nodes, missing_properties, repeated_nodes, zero_area = [], [], [], []
    areas, ratios, lengths, warp_angles, corner_angles = [], [], [], [], []
    per_property = defaultdict(lambda: {"elements": 0, "area": 0.0})
    element_metrics = []
    for eid, element in elements.items():
        ids = element["nodes"]
        absent = [nid for nid in ids if nid not in nodes]
        if element["pid"] not in model["properties"]:
            missing_properties.append(eid)
        if absent:
            missing_nodes.append({"element": eid, "node_ids": absent})
            continue
        if len(set(ids)) != len(ids):
            repeated_nodes.append(eid)
        valid_elements[eid] = element
        used_nodes.update(ids)
        for nid in ids[1:]:
            union(ids[0], nid)
        pts = [nodes[nid]["xyz"] for nid in ids]
        edge_lengths = []
        for i, nid in enumerate(ids):
            other = ids[(i + 1) % len(ids)]
            edges[tuple(sorted((nid, other)))].append((eid, nid, other))
            edge_lengths.append(norm(sub(pts[i], pts[(i + 1) % len(ids)])))
            corner_angles.append(angle(sub(pts[i - 1], pts[i]), sub(pts[(i + 1) % len(ids)], pts[i])))
        normals = [cross(sub(pts[1], pts[0]), sub(pts[2], pts[0]))]
        if len(ids) == 4:
            normals.append(cross(sub(pts[2], pts[0]), sub(pts[3], pts[0])))
            warp_angles.append(angle(*normals))
        area = sum(norm(value) / 2 for value in normals)
        if area == 0:
            zero_area.append(eid)
        ratio = max(edge_lengths) / min(edge_lengths) if min(edge_lengths) > 0 else None
        areas.append(area)
        ratios.append(ratio)
        lengths.extend(edge_lengths)
        per_property[element["pid"]]["elements"] += 1
        per_property[element["pid"]]["area"] += area
        element_metrics.append({"id": eid, "edge_length_ratio": ratio, "area": area})
    component_nodes = defaultdict(list)
    component_elements = Counter()
    for nid in used_nodes:
        component_nodes[find(nid)].append(nid)
    for element in valid_elements.values():
        component_elements[find(element["nodes"][0])] += 1
    components = [{"node_count": len(ids), "element_count": component_elements[root],
                   "first_node_id": min(ids)} for root, ids in component_nodes.items()]
    components.sort(key=lambda value: (-value["element_count"], value["first_node_id"]))
    constraints = {nid: dofs.copy() for nid, dofs in model["permanent"].items()}
    selected_sets = [spc_set] if spc_set is not None else sorted(model["constraints"])
    for sid in selected_sets:
        for nid, dofs in model["constraints"].get(sid, {}).items():
            constraints.setdefault(nid, {}).update(dofs)
    def clamped(nid):
        return all(constraints.get(nid, {}).get(dof) == 0.0 for dof in range(1, 7))
    root_reports, root_ids, root_node_sets = [], set(), []
    for axis, value in root_planes:
        ids = {nid for nid, node in nodes.items() if abs(node["xyz"][axis] - value) <= root_tol}
        root_ids.update(ids)
        root_node_sets.append(ids)
        full = sum(clamped(nid) for nid in ids)
        partial = sum(bool(constraints.get(nid)) and not clamped(nid) for nid in ids)
        root_reports.append({"axis": "xyz"[axis], "coordinate": value, "node_count": len(ids),
            "shell_node_count": len(ids & used_nodes), "zero_displacement_123456_count": full,
            "partial_or_nonzero_constraint_count": partial, "unconstrained_count": len(ids) - full - partial,
            "fully_clamped_fraction": full / len(ids) if ids else None,
            "not_fully_clamped_node_examples": sorted(nid for nid in ids if not clamped(nid))[:100]})
    free_edges = [edge for edge, uses in edges.items() if len(uses) == 1]
    nonmanifold = [edge for edge, uses in edges.items() if len(uses) > 2]
    inconsistent = [edge for edge, uses in edges.items()
                    if len(uses) == 2 and uses[0][1:] == uses[1][1:]]
    absent_materials = []
    for pid, prop in model["properties"].items():
        for key in ("mid1", "mid2", "mid3", "mid4"):
            if prop[key] and prop[key] > 0 and prop[key] not in model["materials"]:
                absent_materials.append({"pid": pid, "field": key, "mid": prop[key]})
    if not nodes or not elements:
        model["warnings"].append("No GRID or supported CQUAD4/CTRIA3 shell elements found.")
    if spc_set is not None and spc_set not in model["constraints"]:
        model["warnings"].append(f"Requested direct SPC set {spc_set} was not found.")
    return {
        "source": str(Path(path).resolve()), "coordinate_units": units,
        "limitations": ["Only CQUAD4/CTRIA3 shells are included in connectivity and surface metrics.",
            "Edge-length ratio and diagonal-triangle warp are simple screening measures, not Nastran criteria.",
            "No coordinate-system transformations or SPCADD expansion are performed.",
            "Without --spc-set, reported constraints are the union of all direct sets, not an inferred active case.",
            "Root planes are user-specified geometric selections; they are not inferred from a bounding box."],
        "warnings": model["warnings"], "card_counts": dict(model["card_counts"]),
        "node_count": len(nodes), "shell_element_count": len(elements),
        "shell_types": dict(Counter(e["type"] for e in elements.values())),
        "bbox": bounds, "coordinates_are_basic": all(n["cp"] == 0 for n in nodes.values()),
        "duplicate_ids": dict(model["duplicate_ids"]),
        "missing_node_references": missing_nodes, "element_ids_missing_pshell": missing_properties,
        "missing_material_references": absent_materials,
        "elements_with_repeated_node_ids": repeated_nodes, "zero_area_element_ids": zero_area,
        "connected_component_count": len(components), "connected_components": components,
        "unused_grid_count": len(nodes.keys() - used_nodes),
        "unused_grid_examples": sorted(nodes.keys() - used_nodes)[:100],
        "edges": {"unique_count": len(edges), "free_count": len(free_edges),
            "nonmanifold_count": len(nonmanifold), "inconsistent_orientation_count": len(inconsistent),
            "free_edge_examples": free_edges[:100], "nonmanifold_edge_examples": nonmanifold[:100],
            "inconsistent_orientation_examples": inconsistent[:100],
            "free_edges_both_nodes_clamped": sum(clamped(a) and clamped(b) for a, b in free_edges),
            "free_edges_on_selected_roots": sum(any(a in ids and b in ids for ids in root_node_sets)
                                                 for a, b in free_edges)},
        "near_coincident_nodes": near_coincident(nodes, coincident_tol),
        "geometry_metrics": {"total_area": sum(areas), "element_area": stats(areas),
            "edge_length": stats(lengths), "longest_to_shortest_edge_ratio": stats(ratios),
            "corner_angle_deg": stats(corner_angles), "quad_diagonal_warp_angle_deg": stats(warp_angles),
            "worst_edge_ratio_elements": sorted(element_metrics,
                key=lambda item: (item["edge_length_ratio"] is None, item["edge_length_ratio"] or 0), reverse=True)[:20]},
        "materials": model["materials"], "plate_properties": model["properties"],
        "property_usage": dict(per_property),
        "constraints": {"direct_set_node_counts": {sid: len(values) for sid, values in model["constraints"].items()},
            "selected_direct_sets": selected_sets, "permanently_constrained_grid_count": len(model["permanent"]),
            "constrained_existing_node_count": sum(nid in nodes for nid in constraints),
            "fully_clamped_existing_node_count": sum(clamped(nid) for nid in nodes),
            "missing_node_ids": sorted(set(constraints) - nodes.keys()), "root_tolerance": root_tol,
            "root_planes": root_reports,
            "fully_clamped_nodes_outside_selected_roots": sum(clamped(nid) and nid not in root_ids for nid in nodes)
                if root_planes else None},
    }


def parse_root(text):
    match = re.fullmatch(r"([xyzXYZ])\s*=\s*(.+)", text)
    if not match:
        raise argparse.ArgumentTypeError("Expected x=number, y=number or z=number")
    try:
        return "xyz".index(match[1].lower()), nastran_float(match[2])
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("bdf", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--units", default="unspecified", help="Coordinate unit label only; no conversion")
    parser.add_argument("--coincident-tol", type=float, default=1e-7, help="Coordinate units; zero disables check")
    parser.add_argument("--root-plane", action="append", type=parse_root, default=[])
    parser.add_argument("--root-tol", type=float, default=1e-6, help="Coordinate units")
    parser.add_argument("--spc-set", type=int, help="Direct SPC/SPC1 SID; default is union of all direct sets")
    args = parser.parse_args()
    if any(not math.isfinite(value) or value < 0 for value in (args.coincident_tol, args.root_tol)):
        parser.error("Tolerances must be finite and nonnegative")
    try:
        result = inspect(args.bdf, args.units, args.coincident_tol, args.root_plane, args.root_tol, args.spc_set)
    except (OSError, ValueError) as exc:
        parser.exit(2, f"BDF inspection failed: {exc}\n")
    output = json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(output + "\n", encoding="utf-8")
        print(f"Wrote {args.output.resolve()}")
    else:
        print(output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
