"""Coarse equilibrium, strain and small-displacement checks of a BWB SOL144 run.

This reads actual APRES/SPCFORCE output and exported CSVs. The finite-wing and
thin-wall estimates are order-of-magnitude comparisons, not convergence or
strength certification. No solver is started and no result files are modified.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import re

import numpy as np

try:
    from .inspect_shell_bdf import cards, load_model, nastran_float
    from .check_aero_panels import grouped_cards, read_aero_geometry
except ImportError:
    from inspect_shell_bdf import cards, load_model, nastran_float
    from check_aero_panels import grouped_cards, read_aero_geometry

ROOT = Path(__file__).resolve().parents[2]
PRINT_ABSOLUTE_FLOOR = 1e-10
PRINT_ARITHMETIC_RTOL = 128 * np.finfo(float).eps


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def plane_stress_von_mises(exx, eyy, gamma, young=69e9, poisson=.33):
    factor = young / (1 - poisson * poisson)
    sx, sy = factor * (exx + poisson * eyy), factor * (eyy + poisson * exx)
    shear = young / (2 * (1 + poisson)) * gamma
    return math.sqrt(max(0., sx * sx - sx * sy + sy * sy + 3 * shear * shear))


def shift_moment(moment_about_zero, force, new_origin):
    return np.asarray(moment_about_zero) - np.cross(new_origin, force)


def printed_roundoff(token):
    """Half the last printed decimal place, retaining each F06 token's exponent."""
    token = token.strip().upper().replace("D", "E")
    if "E" not in token:
        token = re.sub(r"(?<=[0-9.])([+-]\d+)$", r"E\1", token)
    match = re.fullmatch(r"[+-]?(\d*\.?\d+|\d+\.)(?:E([+-]?\d+))?", token)
    require(match is not None and math.isfinite(nastran_float(token)),
            f"Unrecognized finite printed number: {token}")
    mantissa, exponent = match.group(1), int(match.group(2) or 0)
    decimals = len(mantissa.split(".")[1]) if "." in mantissa else 0
    return .5 * 10. ** (exponent - decimals)


def print_tolerance(propagated_roundoff, magnitude):
    """Printed-data uncertainty plus a small absolute and arithmetic-relative floor."""
    return (np.asarray(propagated_roundoff) + PRINT_ABSOLUTE_FLOOR
            + PRINT_ARITHMETIC_RTOL * np.asarray(magnitude))


def cross_magnitude_bound(a, b):
    """Componentwise |a cross b| bound without cancellation; also propagates errors."""
    a, b = np.abs(a), np.abs(b)
    return np.array([a[1] * b[2] + a[2] * b[1],
                     a[2] * b[0] + a[0] * b[2],
                     a[0] * b[1] + a[1] * b[0]])


def pressure_consistency(coefficient, pressure, q, coefficient_roundoff, pressure_roundoff):
    residual = abs(pressure - q * coefficient)
    tolerance = float(print_tolerance(pressure_roundoff + abs(q) * coefficient_roundoff,
                                      abs(pressure) + abs(q * coefficient)))
    return residual, tolerance


def reaction_consistency(points, values, roundoff, listed, listed_roundoff):
    """Propagate printed nodal force/moment uncertainty through their resultant sum."""
    total, budget, magnitude = np.zeros(6), np.zeros(6), np.zeros(6)
    for point, row, error in zip(points, values, roundoff):
        total[:3] += row[:3]
        total[3:] += row[3:] + np.cross(point, row[:3])
        budget[:3] += error[:3]
        budget[3:] += error[3:] + cross_magnitude_bound(point, error[:3])
        magnitude[:3] += np.abs(row[:3])
        magnitude[3:] += np.abs(row[3:]) + cross_magnitude_bound(point, row[:3])
    tolerance = print_tolerance(budget + listed_roundoff, magnitude + np.abs(listed))
    return total - listed, tolerance


def read_load_tables(path, include_precision=False):
    pressures, resultants, reactions = {}, {}, {}
    precision = {"pressures": {}, "resultants": {}, "reactions": {}}
    section, subcase, ended = None, None, False
    with Path(path).open(errors="replace") as stream:
        for line in stream:
            require("FATAL" not in line.upper(), "F06 contains a FATAL diagnostic")
            ended |= "* * * END OF JOB * * *" in line
            match = re.search(r"^0\s+.*?\bSUBCASE\s+(\d+)\s*$", line)
            if match:
                subcase = int(match[1])
            compact = "".join(line.split()).upper().replace("-", "")
            if "AERODYNAMICPRESSURESONTHEAERODYNAMICELEMENTS" in compact:
                section = "pressure"
            elif "SPCFORCERESULTANT" in compact:
                section = "resultant"
            elif "FORCESOFSINGLEPOINTCONSTRAINT" in compact:
                section = "spc"
            elif any(title in compact for title in ("AERODYNAMICFORCES", "STRAINSIN", "STRESSESIN", "DISPLACEMENTVECTOR")):
                section = None
            values = line.split()
            if section == "pressure" and len(values) == 5 and values[0].isdigit() and values[1] == "LS":
                require(subcase is not None, "APRES table has no subcase header")
                key = (subcase, int(values[4]))
                require(key not in pressures, f"Duplicate aerodynamic pressure {key}")
                pressures[key] = (nastran_float(values[2]), nastran_float(values[3]))
                precision["pressures"][key] = tuple(printed_roundoff(v) for v in values[2:4])
            elif section == "resultant" and values[:1] == ["TOTALS"]:
                require(subcase is not None and subcase not in resultants and len(values) == 7,
                        "Unrecognized or duplicated SPCFORCE resultant")
                resultants[subcase] = np.array([nastran_float(v) for v in values[1:]])
                precision["resultants"][subcase] = np.array([printed_roundoff(v) for v in values[1:]])
                section = None
            elif section == "spc" and len(values) == 8 and values[0].isdigit() and values[1] == "G":
                key = (subcase, int(values[0]))
                require(subcase is not None and key not in reactions, f"Duplicate/unidentified SPCFORCE row {key}")
                reactions[key] = np.array([nastran_float(v) for v in values[2:]])
                precision["reactions"][key] = np.array([printed_roundoff(v) for v in values[2:]])
    require(ended and pressures and resultants and reactions, "Incomplete APRES/SPCFORCE/end-of-job output")
    if include_precision:
        return pressures, resultants, reactions, precision
    return pressures, resultants, reactions


def root_section_estimate(model, root_nodes, aero_force, aero_moment_zero, young):
    edge_uses = defaultdict(list)
    for eid, element in model["elements"].items():
        conn = element["nodes"]
        for a, b in zip(conn, conn[1:] + conn[:1]):
            if a in root_nodes and b in root_nodes:
                edge_uses[tuple(sorted((a, b)))].append(eid)
    pieces = []
    for (a, b), attached in edge_uses.items():
        if len(attached) != 1:
            continue
        p, q = (np.array(model["nodes"][n]["xyz"]) for n in (a, b))
        thickness = model["properties"][model["elements"][attached[0]]["pid"]]["thickness"]
        pieces.append((p, q, np.linalg.norm(q - p) * thickness))
    require(len(pieces) == len(root_nodes), "Root boundary is not a complete single closed contour")
    area = sum(a for _, _, a in pieces)
    centroid = sum((p + q) / 2 * a for p, q, a in pieces) / area
    ixx = iyy = ixy = 0.
    for p, q, a in pieces:
        mid = (p + q) / 2 - centroid
        delta = q - p
        ixx += a * (mid[1] ** 2 + delta[1] ** 2 / 12)
        iyy += a * (mid[0] ** 2 + delta[0] ** 2 / 12)
        ixy += a * (mid[0] * mid[1] + delta[0] * delta[1] / 12)
    moment = shift_moment(aero_moment_zero, aero_force, centroid)
    coefficients = np.linalg.solve([[iyy, ixy], [ixy, ixx]], [-moment[1], moment[0]])
    section_stress = [float(coefficients @ (np.array(model["nodes"][n]["xyz"]) - centroid)[:2]) for n in root_nodes]
    return dict(thin_wall_area_m2=float(area), centroid_cad_m=centroid.tolist(),
                I_about_CAD_X_m4=float(ixx), I_about_CAD_Y_m4=float(iyy), Ixy_m4=float(ixy),
                bending_moment_about_section_centroid_Nm=moment.tolist(),
                estimated_max_axial_bending_microstrain=max(abs(v) for v in section_stress) / young * 1e6,
                assumptions="Thin-wall contour integration and unsymmetric beam bending along CAD Z; ignores shear, torsion, local shell bending and concentrations.")


def check(case_dir, export_dir, part_index, job_name="bwb_sol144", yield_mpa=200.):
    case_dir, export_dir = Path(case_dir).resolve(), Path(export_dir).resolve()
    require(part_index in (1, 2), "Part index must be 1 or 2")
    require(math.isfinite(yield_mpa) and yield_mpa > 0, "Positive assumed yield strength required")
    f06 = case_dir / f"{job_name}.f06"
    summary_path = export_dir / "strain_summary.json"
    exported = json.loads(summary_path.read_text(encoding="utf-8"))
    require(exported["source_f06_sha256"] == digest(f06), "Export belongs to a different F06")
    current_inputs = {p.name: digest(p) for p in case_dir.iterdir() if p.suffix.lower() in (".bdf", ".inc")}
    require(current_inputs == exported["input_sha256"], "Case inputs changed after strain export")
    for name, expected in exported["csv_sha256"].items():
        require(digest(export_dir / name) == expected, f"Exported CSV hash mismatch: {name}")
    model = load_model(case_dir / "structural_cards.bdf")
    root_nodes = set(model["constraints"].get(1, {}))
    require(root_nodes, "Missing root SPC set 1")
    require(all(n["cp"] == n["cd"] == 0 for n in model["nodes"].values()), "Expected basic structural coordinates")
    manifest_path = ROOT / f"BWB/models/aero_panels/part_{part_index}/aero_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    boxes = {b["id"]: b for b in manifest["boxes"]}
    groups = grouped_cards(case_dir / "aero_cards.bdf")
    aero = read_aero_geometry(groups)
    require(set(boxes) == set(aero["boxes"]), "Actual CAERO1 box IDs differ from manifest")
    for bid, box in boxes.items():
        actual = aero["boxes"][bid]
        require(np.allclose(actual["corners"], box["corners"], atol=2e-8, rtol=2e-8),
                f"Actual CAERO1 box geometry differs from manifest: {bid}")
        require(np.allclose(actual["normal"], box["normal"], atol=2e-8, rtol=2e-8)
                and math.isclose(actual["area"], box["area_m2"], rel_tol=2e-8, abs_tol=1e-12),
                f"Manifest box normal/area mismatch: {bid}")
        qcorners = actual["corners"]
        point = .375 * (qcorners[0] + qcorners[3]) + .125 * (qcorners[1] + qcorners[2])
        require(np.allclose(point, box["pressure_point"], atol=2e-8, rtol=2e-8),
                f"Manifest quarter-chord pressure point mismatch: {bid}")
    af = groups["AEROS"][0][0]
    root_origin = aero["frames"][int(af[2] or 0)][0]
    root_z = float(np.mean([model["nodes"][nid]["xyz"][2] for nid in root_nodes]))
    require(np.allclose(root_origin, (0, 0, root_z), atol=1e-8), "Unexpected aerodynamic moment reference origin")
    trim = [fields for name, fields, _ in cards(case_dir / "flight_condition.bdf") if name == "TRIM"]
    require(len(trim) == 1, "Expected one assumed TRIM flight condition")
    dynamic_pressure = nastran_float(trim[0][3])
    alpha_index = trim[0].index("ANGLEA")
    alpha = nastran_float(trim[0][alpha_index + 1])
    require(dynamic_pressure > 0, "Nonpositive dynamic pressure")
    pressures, resultants, reactions, precision = read_load_tables(f06, include_precision=True)
    subcases = sorted(resultants)
    require(set(pressures) == {(s, b) for s in subcases for b in boxes}, "APRES box/subcase coverage is incomplete")
    require(set(reactions) == {(s, n) for s in subcases for n in root_nodes}, "SPCFORCE root/subcase coverage differs")
    with (export_dir / "shell_strain_long.csv").open(newline="") as stream:
        strain_rows = list(csv.DictReader(stream))
    with (export_dir / "grid_displacements.csv").open(newline="") as stream:
        displacement_rows = list(csv.DictReader(stream))
    young, poisson = 69e9, .33
    require(all(math.isclose(m["E"], young, rel_tol=1e-8) and math.isclose(m["nu"], poisson, rel_tol=1e-8)
                for m in model["materials"].values()), "Stress reconstruction assumes MAT1 E=69 GPa, nu=0.33")
    area, span = aero["reference"]["area_m2"], aero["reference"]["span_m"]
    aspect_ratio = span * span / area
    slope = 2 * math.pi * aspect_ratio / (2 + math.sqrt(4 + aspect_ratio * aspect_ratio))
    rough_lift = dynamic_pressure * area * slope * alpha
    reports = []
    for subcase in subcases:
        force, moment = np.zeros(3), np.zeros(3)
        cp_checks = []
        for bid, box in boxes.items():
            coefficient, pressure = pressures[(subcase, bid)]
            error, tolerance = pressure_consistency(coefficient, pressure, dynamic_pressure,
                                                     *precision["pressures"][(subcase, bid)])
            cp_checks.append((error, tolerance, bid))
            load = pressure * box["area_m2"] * np.array(box["normal"])
            force += load
            moment += np.cross(box["pressure_point"], load)
        max_cp_error, max_cp_tolerance, max_cp_box = max(cp_checks, key=lambda row: row[0])
        worst_cp = max(cp_checks, key=lambda row: row[0] / row[1])
        require(worst_cp[0] <= worst_cp[1],
                f"APRES pressure is inconsistent with q times Cp: box {worst_cp[2]}, "
                f"error {worst_cp[0]:.8g} Pa exceeds printed-roundoff tolerance {worst_cp[1]:.8g} Pa")
        listed = resultants[subcase]
        node_order = sorted(root_nodes)
        nodal_difference, nodal_tolerance = reaction_consistency(
            [model["nodes"][nid]["xyz"] for nid in node_order],
            [reactions[(subcase, nid)] for nid in node_order],
            [precision["reactions"][(subcase, nid)] for nid in node_order],
            listed, precision["resultants"][subcase])
        nodal_vs_listed = float(np.max(abs(nodal_difference)))
        require(np.all(abs(nodal_difference) <= nodal_tolerance),
                "Nodal root reactions do not reproduce printed SPCFORCE resultant: "
                f"component residuals {nodal_difference.tolist()}, tolerances {nodal_tolerance.tolist()} (N, N m)")
        root_moment = shift_moment(moment, force, root_origin)
        reaction_root_moment = shift_moment(listed[3:], listed[:3], root_origin)
        force_error = float(np.linalg.norm(force + listed[:3]))
        moment_error = float(np.linalg.norm(root_moment + reaction_root_moment))
        force_relative = force_error / max(float(np.linalg.norm(force)), dynamic_pressure * area * 1e-6)
        moment_relative = moment_error / max(float(np.linalg.norm(root_moment)), dynamic_pressure * area * span * 1e-6)
        vm = {surface: (0., None) for surface in ("top", "bottom", "mid")}
        strain_count = 0
        actual_root_peak = 0.
        root_elements = {eid for eid, e in model["elements"].items() if set(e["nodes"]) & root_nodes}
        for row in strain_rows:
            if int(row["subcase"]) != subcase:
                continue
            strain_count += 1
            components = [float(row[key]) for key in ("exx_element", "eyy_element", "gamma_xy_element")]
            require(all(math.isfinite(v) for v in components), "Nonfinite exported strain")
            value = plane_stress_von_mises(*components, young, poisson)
            surface, eid = row["surface"], int(row["element_id"])
            if value > vm[surface][0]:
                vm[surface] = (value, eid)
            if eid in root_elements:
                actual_root_peak = max(actual_root_peak, abs(float(row["dominant_principal_microstrain"])))
        require(strain_count == 3 * len(model["elements"]), "Unexpected exported strain count")
        max_displacement = max_rotation = root_dof = 0.
        displacement_count = 0
        for row in displacement_rows:
            if int(row["subcase"]) != subcase:
                continue
            displacement_count += 1
            values = np.array([float(row[key]) for key in ("ux_m", "uy_m", "uz_m", "rx_rad", "ry_rad", "rz_rad")])
            require(np.all(np.isfinite(values)), "Nonfinite exported displacement")
            max_displacement = max(max_displacement, float(np.linalg.norm(values[:3])))
            max_rotation = max(max_rotation, float(np.linalg.norm(values[3:])))
            if int(row["grid_id"]) in root_nodes:
                root_dof = max(root_dof, float(np.max(abs(values))))
        require(displacement_count == len(model["nodes"]), "Unexpected exported structural displacement count")
        peak_vm = max(v[0] for v in vm.values())
        body_lift = float(force @ aero["up"])
        ratio = body_lift / rough_lift if rough_lift else None
        section = root_section_estimate(model, root_nodes, force, moment, young)
        section["actual_max_abs_root_adjacent_principal_microstrain"] = actual_root_peak
        screens = dict(force_balance=force_relative < 1e-4, moment_balance=moment_relative < 1e-4,
                       root_fixed=root_dof <= 1e-10, small_translation=max_displacement / span < .01,
                       small_rotation=max_rotation < .05, below_assumed_yield=peak_vm < yield_mpa * 1e6,
                       same_lift_order_as_isolated_finite_wing=ratio is not None and .25 < ratio < 4.)
        reports.append(dict(subcase=subcase, screens=screens,
            aerodynamic_force_CAD_N=force.tolist(), aerodynamic_moment_basic_origin_Nm=moment.tolist(),
            aerodynamic_moment_root_origin_Nm=root_moment.tolist(), root_origin_CAD_m=root_origin.tolist(),
            spc_force_CAD_N=listed[:3].tolist(), spc_moment_basic_origin_Nm=listed[3:].tolist(),
            spc_moment_root_origin_Nm=reaction_root_moment.tolist(),
            force_balance_error_N=force_error, force_balance_relative=force_relative,
            moment_balance_error_Nm=moment_error, moment_balance_relative=moment_relative,
            nodal_vs_listed_spc_max_abs_component_difference=nodal_vs_listed,
            max_pressure_vs_qCp_error_Pa=max_cp_error,
            printed_output_consistency=dict(
                pressure_at_max_error_box_id=max_cp_box,
                pressure_tolerance_at_max_error_Pa=max_cp_tolerance,
                pressure_max_error_over_tolerance=worst_cp[0] / worst_cp[1],
                pressure_worst_ratio_box_id=worst_cp[2],
                pressure_worst_ratio_error_Pa=worst_cp[0],
                pressure_worst_ratio_tolerance_Pa=worst_cp[1],
                nodal_minus_listed_force_N=nodal_difference[:3].tolist(),
                nodal_force_tolerance_N=nodal_tolerance[:3].tolist(),
                nodal_minus_listed_moment_Nm=nodal_difference[3:].tolist(),
                nodal_moment_tolerance_Nm=nodal_tolerance[3:].tolist(),
                nodal_max_error_over_tolerance=float(np.max(abs(nodal_difference) / nodal_tolerance))),
            max_translation_m=max_displacement, max_translation_over_reference_span=max_displacement / span,
            max_rotation_rad=max_rotation, max_rotation_deg=math.degrees(max_rotation), max_abs_root_dof=root_dof,
            von_mises_by_surface={surface: dict(peak_MPa=value / 1e6, element_id=eid) for surface, (value, eid) in vm.items()},
            peak_von_mises_MPa=peak_vm / 1e6, assumed_yield_MPa=yield_mpa, yield_utilization=peak_vm / (yield_mpa * 1e6),
            body_up_lift_N=body_lift, body_up_CL=body_lift / (dynamic_pressure * area),
            lift_to_rough_finite_wing_ratio=ratio, root_thin_wall_estimate=section))
    all_screens = all(all(report["screens"].values()) for report in reports)
    diagnostic = bool(exported.get("diagnostic_only"))
    return dict(status="diagnostic_only" if diagnostic else "coarse_checks_pass" if all_screens else "review_needed",
        diagnostic_only=diagnostic, all_coarse_screens_pass=all_screens, part=part_index,
        source_f06_sha256=digest(f06), export_summary_sha256=digest(summary_path),
        aero_manifest_sha256=digest(manifest_path), pressure_boxes_per_subcase=len(boxes), root_reaction_nodes=len(root_nodes),
        material=dict(young_modulus_Pa=young, poisson_ratio=poisson, assumed_yield_MPa=yield_mpa),
        printed_output_tolerance_policy=dict(
            method="Half each F06 token's last printed digit, propagated through q*Cp and sums of forces/moments including moment arms; absolute and floating-point arithmetic floors added.",
            absolute_floor=PRINT_ABSOLUTE_FLOOR, arithmetic_relative_floor=PRINT_ARITHMETIC_RTOL,
            units="Comparison-specific: Pa for pressure, N for force, N m for moment.",
            physical_screen_thresholds_unchanged=True),
        rough_finite_wing=dict(reference_area_m2=area, reference_span_m=span, aspect_ratio=aspect_ratio,
            slope_per_rad=slope, formula="a=2*pi*AR/(2+sqrt(4+AR^2)); L=q*S*a*alpha",
            dynamic_pressure_Pa=dynamic_pressure, alpha_deg=math.degrees(alpha), estimated_lift_N=rough_lift,
            scope="Isolated uncambered unswept finite-wing scale only; ignores CAD camber, sweep, winglet interaction and elastic redistribution."),
        subcases=reports,
        limitations=["A force-balanced result may still contain an incorrect modeling assumption.",
                     "Coarse magnitude checks are not mesh convergence, experimental validation or strength certification.",
                     "Plane-stress von Mises is reconstructed from printed linear isotropic shell strains."])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-dir", type=Path, required=True)
    parser.add_argument("--export-dir", type=Path, required=True)
    parser.add_argument("--part-index", type=int, choices=(1, 2), required=True)
    parser.add_argument("--job-name", default="bwb_sol144")
    parser.add_argument("--assumed-yield-mpa", type=float, default=200.)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        result = check(args.case_dir, args.export_dir, args.part_index, args.job_name, args.assumed_yield_mpa)
    except (ValueError, OSError, KeyError) as exc:
        parser.exit(1, f"Reasonableness check rejected: {exc}\n")
    text = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")
    return 0 if result["all_coarse_screens_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
