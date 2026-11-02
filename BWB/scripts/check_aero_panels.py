"""Independent, read-only checks of generated BWB aerodynamic bulk data.

The checker reads the exported cards, never the generator's implementation.
Passing these checks establishes consistent input data, not aerodynamic fidelity
or a successfully completed Nastran solution.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

try:
    from .inspect_shell_bdf import cards, expanded_ids, load_model, nastran_float
except ImportError:
    from inspect_shell_bdf import cards, expanded_ids, load_model, nastran_float


class ValidationError(ValueError):
    """An exported card or its provenance is inconsistent."""


def require(condition, message):
    if not condition:
        raise ValidationError(message)


def sha256(path, mode="raw"):
    require(mode in ("raw", "crlf_to_lf"), f"Unknown provenance hash mode: {mode}")
    content = Path(path).read_bytes()
    if mode == "crlf_to_lf":
        content = content.replace(b"\r\n", b"\n")
    return hashlib.sha256(content).hexdigest()


def grouped_cards(path):
    grouped = defaultdict(list)
    for name, fields, location in cards(Path(path)):
        grouped[name].append((fields, location))
    return grouped


def check_physical_records(path):
    """Enforce the Simcenter Nastran 80-column input record convention."""
    for number, raw in enumerate(Path(path).read_text(encoding="ascii").splitlines(), 1):
        line = raw.split("$", 1)[0]
        if not line.strip():
            continue
        require(len(line) <= 80, f"{Path(path).name}:{number}: bulk input exceeds 80 columns")
        if "," in line[:8]:
            fields = line.split(",")
            require(all(len(value.strip()) <= 8 for value in fields[1:]),
                    f"{Path(path).name}:{number}: free-field data exceeds eight characters")


def field(fields, index, default=""):
    return fields[index] if index < len(fields) and fields[index] else default


def coordinate_systems(grouped):
    """Resolve right-handed CORD2R frames from the actual bulk data."""
    definitions = {}
    for f, location in grouped.get("CORD2R", []):
        cid = int(f[1])
        require(cid > 0 and cid not in definitions, f"Duplicate/invalid CORD2R {cid}: {location}")
        definitions[cid] = (int(field(f, 2, "0")),
                            np.array([nastran_float(field(f, i)) for i in range(3, 12)]).reshape(3, 3))
    frames = {0: (np.zeros(3), np.eye(3))}

    def resolve(cid, pending=()):
        if cid in frames:
            return frames[cid]
        require(cid in definitions, f"Missing CORD2R {cid}")
        require(cid not in pending, f"Cyclic coordinate reference involving {cid}")
        rid, abc = definitions[cid]
        origin, basis = resolve(rid, pending + (cid,))
        a, b, c = origin + abc @ basis.T
        z = b - a
        require(np.linalg.norm(z) > 1e-12, f"Degenerate CORD2R {cid} z axis")
        z /= np.linalg.norm(z)
        y = np.cross(z, c - a)
        require(np.linalg.norm(y) > 1e-12, f"Degenerate CORD2R {cid} x-z plane")
        y /= np.linalg.norm(y)
        x = np.cross(y, z)
        basis = np.column_stack((x, y, z))
        require(np.linalg.det(basis) > 0.999999999, f"Non-right-handed CORD2R {cid}")
        frames[cid] = (a, basis)
        return frames[cid]

    for cid in definitions:
        resolve(cid)
    return frames


def read_aero_geometry(grouped):
    """Expand uniform CAERO1 entries independently into their numbered boxes."""
    frames = coordinate_systems(grouped)
    require(len(grouped.get("AEROS", [])) == 1, "Expected exactly one AEROS")
    af = grouped["AEROS"][0][0]
    acsid, rcsid = int(field(af, 1, "0")), int(field(af, 2, "0"))
    require(acsid in frames and rcsid in frames, "Undefined AEROS coordinate system")
    flow = frames[acsid][1][:, 0]
    up = frames[acsid][1][:, 2]
    reference = dict(chord_m=nastran_float(field(af, 3)), span_m=nastran_float(field(af, 4)),
                     area_m2=nastran_float(field(af, 5)),
                     symxz=int(field(af, 6, "0")), symxy=int(field(af, 7, "0")))
    require(all(reference[k] > 0 for k in ("chord_m", "span_m", "area_m2")),
            "Nonpositive AEROS reference dimension")
    properties = [int(f[1]) for f, _ in grouped.get("PAERO1", [])]
    require(len(properties) == len(set(properties)), "Duplicate PAERO1 ID")
    macro, boxes = {}, {}
    for f, location in grouped.get("CAERO1", []):
        eid, pid, cp = int(f[1]), int(f[2]), int(field(f, 3, "0"))
        ns, nc = int(field(f, 4, "0")), int(field(f, 5, "0"))
        require(eid > 0 and eid not in macro, f"Duplicate/invalid CAERO1 {eid}")
        require(pid in properties, f"CAERO1 {eid} references missing PAERO1 {pid}")
        require(cp in frames, f"CAERO1 {eid} references missing CP {cp}")
        require(ns > 0 and nc > 0, f"CAERO1 {eid} must use positive uniform divisions")
        require(int(field(f, 6, "0")) == 0 and int(field(f, 7, "0")) == 0,
                f"CAERO1 {eid}: AEFACT spacing is not supported by this checker")
        p1 = np.array([nastran_float(field(f, i)) for i in (9, 10, 11)])
        p4 = np.array([nastran_float(field(f, i)) for i in (13, 14, 15)])
        origin, basis = frames[cp]
        p1, p4 = origin + basis @ p1, origin + basis @ p4
        c1, c4 = nastran_float(field(f, 12)), nastran_float(field(f, 16))
        require(c1 > 0 and c4 > 0, f"CAERO1 {eid}: nonpositive chord")
        normal = np.cross(flow, p4 - p1)
        height = float(np.linalg.norm(normal))
        require(height > 1e-12, f"CAERO1 {eid}: zero aerodynamic span")
        normal /= height
        points = []
        for j in range(ns + 1):
            t = j / ns
            le, chord = p1 * (1 - t) + p4 * t, c1 * (1 - t) + c4 * t
            points.append([le + flow * chord * (i / nc) for i in range(nc + 1)])
        points = np.asarray(points)
        ids = []
        for j in range(ns):
            for i in range(nc):
                bid = eid + j * nc + i
                require(bid not in boxes, f"Overlapping aerodynamic box ID {bid}")
                quad = points[[j, j, j + 1, j + 1], [i, i + 1, i + 1, i]]
                area_vector = (np.cross(quad[1] - quad[0], quad[2] - quad[0])
                               + np.cross(quad[2] - quad[0], quad[3] - quad[0])) / 2
                area = float(np.linalg.norm(area_vector))
                require(area > 1e-14, f"Box {bid}: zero area")
                require(float(np.dot(area_vector, normal)) > 0, f"Box {bid}: reversed corner order")
                boxes[bid] = dict(caero=eid, span_index=j, chord_index=i, corners=quad,
                                  normal=normal, area=area,
                                  projected_area=abs(float(np.dot(area_vector, up))))
                ids.append(bid)
        macro[eid] = dict(p1=p1, p4=p4, c1=c1, c4=c4, points=points,
                          pid=pid, cp=cp, igid=int(field(f, 8, "0")),
                          nspan=ns, nchord=nc, normal=normal, ids=ids,
                          area=height * (c1 + c4) / 2, location=location)
    require(macro, "No CAERO1 geometry")
    return dict(macros=macro, boxes=boxes, reference=reference, flow=flow, up=up,
                frames=frames)


def close(actual, expected, description, atol=2e-8, rtol=2e-8):
    require(np.allclose(actual, expected, atol=atol, rtol=rtol),
            f"{description}: actual {actual!r}, expected {expected!r}")


def read_w2gj(grouped, count):
    header, entries = [], {}
    for f, location in grouped.get("DMI", []):
        if field(f, 1).upper() != "W2GJ":
            continue
        column = int(field(f, 2))
        if column == 0:
            header.append(f)
            continue
        require(column == 1, f"W2GJ has unexpected column {column}: {location}")
        row = int(field(f, 3))
        require(row not in entries, f"Duplicate W2GJ row {row}")
        require(1 <= row <= count, f"W2GJ row {row} is out of bounds")
        # The generator deliberately emits one scalar per DMI data entry.
        require(not any(f[5:]), f"Unexpected compressed W2GJ entry: {location}")
        entries[row] = nastran_float(field(f, 4))
    require(len(header) == 1, "Expected exactly one W2GJ DMI header")
    f = header[0]
    require(int(field(f, 3)) == 2, "W2GJ must be a rectangular matrix")
    require(int(field(f, 4)) in (1, 2), "W2GJ must have real entries")
    require(not field(f, 6), "W2GJ DMI header reserved field 6 must be blank")
    require(int(field(f, 7)) == count and int(field(f, 8)) == 1,
            "W2GJ dimensions do not match aerodynamic box count")
    require(set(entries) == set(range(1, count + 1)), "Missing W2GJ rows")
    return entries


def source_path(repo_root, record):
    path = (repo_root / record["path"]).resolve()
    require(path.is_relative_to(repo_root.resolve()), f"Source escapes repository: {record['path']}")
    require(path.is_file(), f"Missing provenance source: {path}")
    require(sha256(path, record.get("hash_mode", "raw")) == record["sha256"],
            f"Source SHA256 mismatch: {record['path']}")
    return path


def check_aero(folder, repo_root=None):
    folder = Path(folder).resolve()
    repo_root = Path(repo_root).resolve() if repo_root else Path(__file__).resolve().parents[2]
    manifest = json.loads((folder / "aero_manifest.json").read_text(encoding="utf-8"))
    require(manifest["schema_version"] == 1, "Unsupported manifest schema")
    require(manifest["units"] == "m-N-kg-s", "Expected SI units")
    sources = [source_path(repo_root, record) for record in manifest["source_files"]]
    references = [source_path(repo_root, record) for record in manifest.get("nasa_reference_files", [])]
    if "generator" in manifest:
        source_path(repo_root, manifest["generator"])
    bdfs = [path for path in sources if path.suffix.lower() == ".bdf"]
    layouts = [path for path in sources if path.name == "mesh_layout.json"]
    require(len(bdfs) == len(layouts) == 1, "Expected one source structural BDF and layout")
    source = load_model(bdfs[0])
    layout = json.loads(layouts[0].read_text(encoding="utf-8"))
    files = {name: folder / filename for name, filename in manifest["files"].items()}
    require(all(path.is_file() for path in files.values()), "Missing generated bulk data file")
    for key in ("aero", "spline", "camber"):
        check_physical_records(files[key])
    groups = grouped_cards(files["aero"])
    aero = read_aero_geometry(groups)
    macros, boxes = aero["macros"], aero["boxes"]
    require(np.allclose(aero["flow"], (1, 0, 0)), "Expected CAD +X flow")
    require(np.allclose(aero["up"], (0, 1, 0)), "Expected CAD +Y aerodynamic up")
    require(aero["reference"]["symxz"] == aero["reference"]["symxy"] == 0,
            "Independent clamped parts must not silently enable aerodynamic symmetry")
    panel_records = manifest["panels"]
    require(len(panel_records) == len(macros), "Manifest panel count mismatch")
    require({p["eid"] for p in panel_records} == set(macros), "Manifest CAERO IDs mismatch")
    for p in panel_records:
        eid, zone = p["eid"], p["zone"]
        actual = macros[eid]
        for key in ("nspan", "nchord", "pid", "cp", "igid"):
            require(actual[key] == p[key], f"CAERO1 {eid}: manifest {key} mismatch")
        close(actual["p1"], p["p1"], f"CAERO1 {eid} P1")
        close(actual["p4"], p["p4"], f"CAERO1 {eid} P4")
        close(actual["c1"], p["x12"], f"CAERO1 {eid} X12")
        close(actual["c4"], p["x43"], f"CAERO1 {eid} X43")
        close(actual["normal"], p["normal"], f"CAERO1 {eid} normal")
        require(actual["normal"][1] > 0, f"CAERO1 {eid}: downward/inconsistent normal")
        require(actual["ids"] == p["box_ids"], f"CAERO1 {eid}: wrong box numbering")
        lower, upper = (layout["patch_node_ids"][f"{zone}_{side}"] for side in ("lower", "upper"))
        a, b = p["source_row_range"]
        require(0 <= a < b < len(lower), f"CAERO1 {eid}: invalid source row range")
        order = (b, a) if p["span_reversed"] else (a, b)
        for row, point, chord in zip(order, (actual["p1"], actual["p4"]), (actual["c1"], actual["c4"])):
            le = np.mean([source["nodes"][grid[row][0]]["xyz"] for grid in (lower, upper)], axis=0)
            te = np.mean([source["nodes"][grid[row][-1]]["xyz"] for grid in (lower, upper)], axis=0)
            close(point, le, f"CAERO1 {eid}: source leading-edge location")
            close(chord, te[0] - le[0], f"CAERO1 {eid}: source projected chord")

    # Independently find complete shared chord edges; the eleven macro faces
    # must form one open chain with ten joints and no branches.
    ends = [(eid, end, values) for eid, p in macros.items()
            for end, values in ((0, p["points"][0]), (1, p["points"][-1]))]
    adjacency, end_uses, joints = defaultdict(set), Counter(), []
    for i, (eid, end, points) in enumerate(ends):
        for other, other_end, other_points in ends[i + 1:]:
            if eid != other and points.shape == other_points.shape and np.allclose(points, other_points, atol=2e-8, rtol=0):
                adjacency[eid].add(other)
                adjacency[other].add(eid)
                end_uses[(eid, end)] += 1
                end_uses[(other, other_end)] += 1
                joints.append((eid, other))
                require(np.dot(macros[eid]["normal"], macros[other]["normal"]) > 0,
                        f"Opposite normals at joint {eid}/{other}")
    require(all(n == 1 for n in end_uses.values()), "Branching/overlapping aerodynamic edges")
    if len(macros) > 1:
        require(len(joints) == len(macros) - 1, "Missing or duplicated aerodynamic interfaces")
        require(sorted(len(adjacency[e]) for e in macros) == [1, 1] + [2] * (len(macros) - 2),
                "Aerodynamic macros do not form a single open chain")
        seen, pending = set(), [next(iter(macros))]
        while pending:
            eid = pending.pop()
            if eid not in seen:
                seen.add(eid)
                pending.extend(adjacency[eid] - seen)
        require(seen == set(macros), "Disconnected aerodynamic macros")

    developed = sum(box["area"] for box in boxes.values())
    projected = sum(box["projected_area"] for box in boxes.values())
    z = [p[key][2] for p in macros.values() for key in ("p1", "p4")]
    span = max(z) - min(z)
    chord = sum(abs(p["p4"][2] - p["p1"][2]) *
                (p["c1"] ** 2 + p["c1"] * p["c4"] + p["c4"] ** 2) / 3
                for p in macros.values()) / projected
    reference = manifest["reference_geometry"]
    require(reference["area_convention"] == "horizontal_projection_including_winglet",
            "Unsupported reference area convention")
    for actual, manifest_key, card_key in ((projected, "sref_m2", "area_m2"),
                                          (span, "bref_m", "span_m"),
                                          (chord, "cref_m", "chord_m")):
        close(actual, reference[manifest_key], f"Manifest {manifest_key}")
        close(actual, aero["reference"][card_key], f"AEROS {card_key}")
    close(developed, reference["developed_area_m2"], "Developed aerodynamic area")

    structural = load_model(files["structure"])
    require(not any(structural["duplicate_ids"].values()), "Duplicate structural entity IDs")
    require(len(structural["nodes"]) == manifest["structural_counts"]["nodes"], "Structural node count mismatch")
    require(len(structural["elements"]) == manifest["structural_counts"]["elements"], "Structural element count mismatch")
    for key in ("nodes", "elements", "materials", "properties", "constraints", "permanent"):
        require(structural[key] == source[key], f"Structural {key} changed while preparing aerodynamic inputs")
    splines = grouped_cards(files["spline"])
    sets = {}
    for f, _ in splines.get("SET1", []):
        sid = int(f[1])
        require(sid not in sets, f"Duplicate SET1 {sid}")
        ids = expanded_ids(f[2:])
        require(len(ids) == len(set(ids)), f"Duplicate node in SET1 {sid}")
        require(set(ids) <= structural["nodes"].keys(), f"SET1 {sid}: unknown structural nodes")
        require(len(ids) >= 3, f"SET1 {sid}: insufficient nodes for surface spline")
        xyz = np.array([structural["nodes"][n]["xyz"] for n in ids])
        require(np.linalg.matrix_rank(xyz - xyz.mean(axis=0), tol=1e-9) >= 2,
                f"SET1 {sid}: collinear support nodes")
        sets[sid] = ids
    spline_ids, covered = set(), Counter()
    by_eid = {p["eid"]: p for p in panel_records}
    for f, _ in splines.get("SPLINE1", []):
        sid, caero, box1, box2, setg = (int(f[i]) for i in range(1, 6))
        require(sid not in spline_ids, f"Duplicate SPLINE1 {sid}")
        spline_ids.add(sid)
        require(caero in macros and setg in sets, f"SPLINE1 {sid}: missing CAERO1 or SET1")
        require(box1 <= box2, f"SPLINE1 {sid}: reversed box range")
        ids = set(range(box1, box2 + 1))
        require(ids == set(macros[caero]["ids"]), f"SPLINE1 {sid}: range does not match its CAERO1")
        covered.update(ids)
        record = by_eid[caero]
        require(sid == record["spline_id"] and setg == record["set_id"], f"SPLINE1 {sid}: manifest ID mismatch")
        require(set(sets[setg]) == set(record["structural_node_ids"]), f"SPLINE1 {sid}: manifest support mismatch")
        xyz = np.array([structural["nodes"][n]["xyz"] for n in sets[setg]])
        span_axis = np.cross(macros[caero]["normal"], aero["flow"])
        projected_nodes = np.column_stack((xyz @ aero["flow"], xyz @ span_axis))
        require(np.linalg.matrix_rank(projected_nodes - projected_nodes.mean(axis=0), tol=1e-9) == 2,
                f"SPLINE1 {sid}: support is collinear in aerodynamic plane")
        rounded = {tuple(np.rint(point / 1e-8).astype(np.int64)) for point in projected_nodes}
        require(len(rounded) == len(sets[setg]), f"SPLINE1 {sid}: coincident projected support points")
        require(nastran_float(field(f, 6)) >= 0, f"SPLINE1 {sid}: negative attachment flexibility")
        require(field(f, 7, "IPS") == "IPS" and field(f, 8, "BOTH") == "BOTH",
                f"SPLINE1 {sid}: unexpected interpolation or load/displacement usage")
    require(set(covered) == set(boxes) and all(n == 1 for n in covered.values()),
            "Spline coupling must cover every aerodynamic box exactly once")
    require(len(spline_ids) == len(macros) == len(sets), "Unexpected spline/set count")

    w2gj = read_w2gj(grouped_cards(files["camber"]), len(boxes))
    # Independently recover the documented initial normalwash from source
    # section means. This does not call the generator or trust its sampled data.
    sampled_wash = {}
    for p in panel_records:
        macro = macros[p["eid"]]
        lower, upper = (layout["patch_node_ids"][f"{p['zone']}_{side}"] for side in ("lower", "upper"))
        first, last = p["source_row_range"]
        vector = macro["p4"] - macro["p1"]
        ts, row_slopes = [], []
        xi = (np.arange(macro["nchord"]) + 0.75) / macro["nchord"]
        for row in range(first, last + 1):
            xyz = np.array([(np.array(source["nodes"][a]["xyz"]) +
                             np.array(source["nodes"][b]["xyz"])) / 2
                            for a, b in zip(lower[row], upper[row])])
            chord_length = xyz[-1, 0] - xyz[0, 0]
            require(chord_length > 0, "Nonpositive source section chord")
            u = (xyz[:, 0] - xyz[0, 0]) / chord_length
            require(np.all(np.diff(u) > 0), "Nonmonotone source chord coordinates")
            h = (xyz - xyz[0]) @ macro["normal"]
            derivative = PchipInterpolator(u, h).derivative()(xi)
            row_slopes.append(-derivative / chord_length)
            ts.append(float(np.dot(xyz[0] - macro["p1"], vector) / np.dot(vector, vector)))
        order = np.argsort(ts)
        ts, row_slopes = np.array(ts)[order], np.array(row_slopes)[order]
        require(np.all(np.diff(ts) > 0), "Source span coordinates are not strictly increasing")
        for j in range(macro["nspan"]):
            eta = (j + 0.5) / macro["nspan"]
            for i in range(macro["nchord"]):
                bid = p["eid"] + j * macro["nchord"] + i
                sampled_wash[bid] = float(np.interp(eta, ts, row_slopes[:, i]))
    records = manifest["boxes"]
    require([r["id"] for r in records] == sorted(boxes), "Manifest box order is not ascending aerodynamic ID")
    for row, record in enumerate(records, 1):
        box, bid = boxes[record["id"]], record["id"]
        require(box["caero"] == record["panel_id"], f"Box {bid}: wrong parent CAERO1")
        require(box["span_index"] == record["span_index"] and box["chord_index"] == record["chord_index"],
                f"Box {bid}: wrong span/chord indexing")
        close(box["corners"], record["corners"], f"Box {bid} corners")
        close(box["normal"], record["normal"], f"Box {bid} normal")
        close(box["area"], record["area_m2"], f"Box {bid} area", atol=1e-11)
        q = box["corners"]
        le, te = (q[0] + q[3]) / 2, (q[1] + q[2]) / 2
        close(le + 0.75 * (te - le), record["control_point"], f"Box {bid} control point")
        close(le + 0.25 * (te - le), record["pressure_point"], f"Box {bid} pressure point")
        close(w2gj[row], record["w2gj_rad"], f"Box {bid} W2GJ row {row}", atol=1e-11)
        close(w2gj[row], sampled_wash[bid], f"Box {bid} source-derived W2GJ", atol=1e-10)

    combined = grouped_cards(files["combined"])
    combined_geometry = read_aero_geometry(combined)
    require(set(combined_geometry["boxes"]) == set(boxes), "Combined BDF aerodynamic boxes differ")
    for name, expected in (("CAERO1", len(macros)), ("SPLINE1", len(macros)),
                           ("SET1", len(sets)), ("GRID", len(structural["nodes"])),
                           ("AEROS", 1)):
        require(len(combined.get(name, [])) == expected, f"Combined BDF {name} count mismatch")
    require(read_w2gj(combined, len(boxes)) == w2gj, "Combined BDF W2GJ differs")
    return dict(status="PASS", part=manifest["part"], units=manifest["units"],
                caero1_macros=len(macros), aerodynamic_boxes=len(boxes),
                shared_chord_interfaces=len(joints), spline1_count=len(spline_ids),
                structural_nodes=len(structural["nodes"]), structural_elements=len(structural["elements"]),
                w2gj_rows=len(w2gj), w2gj_min_rad=min(w2gj.values()), w2gj_max_rad=max(w2gj.values()),
                developed_area_m2=developed, projected_reference_area_m2=projected,
                reference_span_m=span, reference_chord_m=chord,
                source_hashes_checked=len(sources), nasa_reference_hashes_checked=len(references),
                generator_hash_checked="generator" in manifest,
                checks=["80-column records and fixed/free field widths", "exported CAERO1 geometry and box ordering", "CAD source leading edges and projected chords",
                        "shared edges and upward continuous normals", "projected reference area/span/chord",
                        "unchanged structural cards and existing spline support nodes",
                        "complete single spline coverage", "W2GJ matrix dimensions and complete row ordering",
                        "independent source-section PCHIP normalwash reconstruction",
                        "combined INCLUDE expansion", "source SHA256 provenance"],
                limitations=["No Nastran solve, aerodynamic convergence or physical validation is implied.",
                             "W2GJ values are checked against the documented sampled correction, not against measured aerodynamic data."])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path, help="Generated part folder containing aero_manifest.json")
    parser.add_argument("--repo-root", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        result = check_aero(args.folder, args.repo_root)
    except (ValueError, KeyError, OSError, TypeError) as exc:
        result = dict(status="FAIL", error=str(exc))
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
