"""Import authoritative BWB Nastran cards as native Femap aero entities.

This creates a separate model; it never opens or edits the skin .modfem files.
CAERO1 entities remain Aero Panels (type 95), not CQUAD4 structural elements.
The input BDF is the authority for analysis. Saving this viewing database does
not run a solver or certify a subsequent Femap re-export of every input card.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from femap_connection import ROOT, check, connect, ids


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def entity_ids(app, kind: int) -> list[int]:
    selected = app.feSet
    check(selected.AddAll(kind), f"select Femap entity type {kind}")
    return ids(selected)


def imported_panels(app) -> list[dict]:
    panel = app.feAeroPanel
    result = []
    for eid in entity_ids(app, 95):
        check(panel.Get(eid), f"read native Aero Panel {eid}")
        result.append({
            "eid": eid, "title": panel.title, "native_type": panel.type,
            "pid": panel.pid, "cp": panel.icp, "igid": panel.nIgid,
            "nspan": panel.PanelSpanCount(), "nchord": panel.PanelChordCount(),
            "p1": list(panel.vPt1), "p4": list(panel.vPt4),
            "x12": panel.dChord12, "x43": panel.dChord34,
        })
    return result


def audit_camber_storage(app, source: Path) -> dict:
    """Report observed storage, without claiming an untested DMI round trip."""
    matrices = []
    matrix = app.feMatrixInput
    for entity_id in entity_ids(app, 150):
        check(matrix.Get(entity_id), "read imported direct matrix")
        matrices.append({"id": entity_id, "name": matrix.Name})
    text_blocks = []
    manager = app.feAnalysisMgr
    for entity_id in entity_ids(app, 60):
        check(manager.Get(entity_id), "read imported analysis settings")
        for field in ("NasExecStartText", "NasExecEndText", "NasBulkStartText", "NasBulkEndText"):
            text_id = getattr(manager, field)
            if text_id > 0:
                rc, value = manager.GetText(text_id)
                check(rc, f"read imported {field}")
                if "W2GJ" in value.upper():
                    text_blocks.append({"analysis_id": entity_id, "field": field, "text": value})
    source_lines = [line.strip() for line in source.read_text(encoding="utf-8").splitlines()
                    if line.strip() and not line.lstrip().startswith("$")]
    stored_lines = {"".join(line.upper().split()) for block in text_blocks
                    for line in block["text"].splitlines()}
    all_text_retained = bool(source_lines) and all("".join(line.upper().split()) in stored_lines
                                                for line in source_lines)
    return {
        "source": str(source.relative_to(ROOT)), "source_sha256": sha256(source),
        "native_matrix_entities": matrices,
        "W2GJ_found_as_native_matrix": any(item["name"].strip().upper() == "W2GJ" for item in matrices),
        "W2GJ_found_in_analysis_text": bool(text_blocks),
        "all_source_noncomment_lines_retained_in_analysis_text": all_text_retained,
        "analysis_text_locations": [{key: value for key, value in block.items() if key != "text"}
                                    for block in text_blocks],
        "solver_round_trip_verified": False,
        "instruction": "Use authoritative dmicards.bdf; native viewing model is not a replacement analysis deck",
    }


def audit_splines(app, expected: dict[int, dict]) -> list[dict]:
    spline_ids = entity_ids(app, 98)
    if set(spline_ids) != {panel["spline_id"] for panel in expected.values()}:
        raise RuntimeError("Imported native Aero Spline IDs differ from the manifest")
    native = app.feAeroSpline
    result = []
    for panel in expected.values():
        check(native.Get(panel["spline_id"]), "read imported aerodynamic spline")
        box_ids = ids(native.GetBoxSet())
        node_ids = ids(native.GetNodeSet())
        if (native.icaero != panel["eid"] or native.isetg != panel["set_id"]
                or set(box_ids) != set(panel["box_ids"])
                or set(node_ids) != set(panel["structural_node_ids"])):
            raise RuntimeError(f"Spline {panel['spline_id']} imported references differ from manifest")
        result.append({"id": panel["spline_id"], "caero": native.icaero,
                       "set_id": native.isetg, "box_count": len(box_ids),
                       "structural_node_count": len(node_ids),
                       "box_ids_match": True, "structural_node_ids_match": True})
    return result


def prepare_view(app, show_structure: bool) -> dict:
    rc, view_id = app.feAppGetActiveView()
    check(rc, "get active view")
    view = app.feView
    check(view.Get(view_id), "read active view")
    view.title = "BWB native CAERO1 - geometry review, not solved"
    view.Mode, view.Contour, view.Deformed, view.BC = 0, 0, 0, 0
    # FVI_AERO_PANEL / FVI_AERO_MESH / FVI_ELEMENT from Femap 2406 typelib.
    for option in (102, 103):
        view.SetDraw(option, True)
        view.SetLabel(option, 0)
    view.SetDraw(9, show_structure)
    view.SetLabel(9, 0)
    for option in (2, 3, 4, 5, 7, 8, 104, 105, 108):
        view.SetDraw(option, False)
    check(view.Put(view_id), "save aerodynamic viewing options")
    orientation = app.feViewOrient
    check(orientation.Get(view_id), "read view orientation")
    orientation.Perspective = False
    orientation.SetRotationAngles(-65, 20, -25)
    check(orientation.Put(view_id), "set view orientation")
    app.feAppVisible(True)
    app.feViewAutoscaleAll(view_id, True)
    app.feViewRegenerate(view_id)
    return {"view_id": view_id, "aero_mesh_visible": True,
            "structural_mesh_visible": show_structure,
            "structural_mesh_retained_in_database": True}


def build(app, deck: Path, output: Path, manifest: dict, show_structure: bool) -> dict:
    sources = {"combined": deck}
    sources.update({key: deck.parent / value for key, value in manifest["files"].items()})
    source_hashes = {key: sha256(path) for key, path in sources.items()}
    check(app.feFileNew(), "create separate aerodynamic database")
    check(app.feFileReadNastran(8, str(deck)), "import authoritative NX Nastran BDF")
    actual = imported_panels(app)
    if not actual:
        raise RuntimeError("Femap imported no native Aero Panels; refusing a misleading view")
    if any(panel["native_type"] != 0 for panel in actual):
        raise RuntimeError("Unexpected aero body type; this model requires CAERO1 panels")

    expected = {int(panel["eid"]): panel for panel in manifest["panels"]}
    if set(expected) != {panel["eid"] for panel in actual}:
        raise RuntimeError("Imported native Aero Panel IDs differ from the manifest")
    zone_colors = {"main": 110, "transition": 17, "bend": 39, "winglet": 103}
    zone_layers = {zone: 201 + index for index, zone in enumerate(zone_colors)}
    layer = app.feLayer
    for zone, layer_id in zone_layers.items():
        layer.title = f"AERO {zone} - native CAERO1"
        check(layer.Put(layer_id), "create aerodynamic zone layer")
    for panel in actual:
        spec = expected[panel["eid"]]
        for field in ("nspan", "nchord", "cp", "pid", "igid"):
            if field in spec and panel[field] != spec[field]:
                raise RuntimeError(f"Panel {panel['eid']} {field}: {panel[field]} != {spec[field]}")
        for field in ("p1", "p4"):
            if field in spec and max(abs(a-b) for a, b in zip(panel[field], spec[field])) > 1e-8:
                raise RuntimeError(f"Panel {panel['eid']} imported {field} coordinates differ")
        for field in ("x12", "x43"):
            if field in spec and abs(panel[field] - spec[field]) > 1e-8:
                raise RuntimeError(f"Panel {panel['eid']} imported {field} differs")
        native = app.feAeroPanel
        check(native.Get(panel["eid"]), "read panel for descriptive title")
        native.title = spec.get("name", f"CAERO1 {panel['eid']}")
        if spec.get("zone") in zone_colors:
            native.color = zone_colors[spec["zone"]]
            native.layer = zone_layers[spec["zone"]]
        check(native.Put(panel["eid"]), "set native aerodynamic title")
        panel["title"] = native.title

    node_count, element_count = len(entity_ids(app, 7)), len(entity_ids(app, 8))
    for field, actual_count in (("nodes", node_count), ("elements", element_count)):
        if actual_count != manifest["structural_counts"][field]:
            raise RuntimeError(f"Structural {field}: {actual_count} != {manifest['structural_counts'][field]}")
    panel_boxes = sum(panel["nspan"] * panel["nchord"] for panel in actual)
    if panel_boxes != sum(panel["nspan"] * panel["nchord"] for panel in expected.values()):
        raise RuntimeError("Imported aerodynamic box count differs from manifest")
    spline_audit = audit_splines(app, expected)
    camber_storage = audit_camber_storage(app, sources["camber"])
    view = prepare_view(app, show_structure)
    output.parent.mkdir(parents=True, exist_ok=True)
    check(app.feFileSaveAs(False, str(output)), "save independent native aero model")
    if {key: sha256(path) for key, path in sources.items()} != source_hashes:
        raise RuntimeError("Source input BDF changed during viewing-model generation")
    report = {
        "input_bdf": str(deck.relative_to(ROOT)), "input_sha256": source_hashes["combined"],
        "source_hashes": source_hashes,
        "output_model": str(output.relative_to(ROOT)),
        "purpose": "Native aerodynamic entity review; BDF remains analysis authority",
        "native_aero_entity_type": 95, "aero_panels": len(actual),
        "aero_boxes": panel_boxes, "structural_nodes": node_count,
        "structural_elements": element_count, "pseudo_structural_aero_elements_created": 0,
        "native_aero_splines": len(spline_audit), "spline_audit": spline_audit,
        "camber_storage_audit": camber_storage,
        "panels": actual, "view": view, "solved": False,
        "source_bdf_unchanged": True,
        "validation": "Imported native panel IDs, divisions and supplied geometry match manifest",
    }
    (output.parent / "femap_import_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deck", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--show-structure", action="store_true")
    args = parser.parse_args()
    deck, manifest_path, output = args.deck.resolve(), args.manifest.resolve(), args.output.resolve()
    if output.suffix.lower() != ".modfem":
        parser.error("--output must be a .modfem file")
    if not output.is_relative_to(ROOT / "BWB/models") or output.is_relative_to(ROOT / "BWB/models/skin_structured"):
        parser.error("Save to a dedicated BWB/models aero folder, separate from the source skin models")
    if not deck.is_file() or not manifest_path.is_file():
        parser.error("The input deck and manifest must already exist")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    app = connect(True)
    previous = (app.DialogAutoAnswer, app.DialogAutoSkip, app.DialogAutoSkipMsg)
    try:
        app.DialogAutoAnswer = app.DialogAutoSkip = app.DialogAutoSkipMsg = 1
        report = build(app, deck, output, manifest, args.show_structure)
        print(json.dumps({key: value for key, value in report.items() if key != "panels"}, indent=2))
    finally:
        app.DialogAutoAnswer, app.DialogAutoSkip, app.DialogAutoSkipMsg = previous


if __name__ == "__main__":
    main()
