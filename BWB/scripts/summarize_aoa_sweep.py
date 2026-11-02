"""Summarize temporary BWB 3/4/5/6-degree trials against the existing 2-degree run.

No solver is started. All ten result folders must exist before any aggregate
artifact is written. Component-wise responses are checked for affine dependence
on angle of attack; extrema and vector magnitudes are not required to be linear.
Trial inputs/results belong under workspace/BWB/trials, not the project models.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RESULTS = ROOT / "workspace/BWB/trials/results/aoa_sweep"
ANGLES = (2, 3, 4, 5, 6)
DISPLACEMENT_COMPONENTS = ("ux_m", "uy_m", "uz_m", "rx_rad", "ry_rad", "rz_rad")
STRAIN_COMPONENTS = ("exx_element", "eyy_element", "gamma_xy_element")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify_recorded_hash(path, expected):
    """Permit only LF/CRLF checkout conversion, never a data-content change."""
    raw = Path(path).read_bytes()
    lf = raw.replace(b"\r\n", b"\n")
    variants = (raw, lf, lf.replace(b"\n", b"\r\n"))
    require(any(hashlib.sha256(value).hexdigest() == expected for value in variants),
            f"Recorded result hash mismatch: {path}")
    return "exact" if hashlib.sha256(raw).hexdigest() == expected else "LF/CRLF-equivalent"


def read_csv(path):
    with Path(path).open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def load_case(folder, alpha, part):
    paths = {name: folder / filename for name, filename in (
        ("summary", "strain_summary.json"), ("reason", "reasonableness.json"),
        ("strain", "shell_strain_long.csv"), ("displacement", "grid_displacements.csv"))}
    require(all(path.is_file() for path in paths.values()), f"Results incomplete: {folder}")
    summary = json.loads(paths["summary"].read_text(encoding="utf-8-sig"))
    reason = json.loads(paths["reason"].read_text(encoding="utf-8-sig"))
    require(summary["solver_end_of_job"] and not summary["fatal_messages"], f"Incomplete/fatal solver case: {folder}")
    require(len(summary["subcases"]) == len(reason["subcases"]) == 1, f"Expected one subcase: {folder}")
    require(reason["part"] == part, f"Wrong part in {folder}")
    actual_alpha = reason["rough_finite_wing"]["alpha_deg"]
    require(abs(actual_alpha - alpha) < 1e-7, f"Wrong angle of attack in {folder}: {actual_alpha}")
    require(reason["source_f06_sha256"] == summary["source_f06_sha256"], f"F06 provenance differs: {folder}")
    hash_checks = {"summary": verify_recorded_hash(paths["summary"], reason["export_summary_sha256"])}
    for name in ("strain", "displacement"):
        hash_checks[name] = verify_recorded_hash(paths[name], summary["csv_sha256"][paths[name].name])
    review_path = folder / "case_review.json"
    review = None
    if review_path.is_file():
        review = json.loads(review_path.read_text(encoding="utf-8-sig"))
        require(review["source_f06_sha256"] == summary["source_f06_sha256"], f"Case review belongs to a different F06: {folder}")
        verify_recorded_hash(paths["reason"], review["reasonableness_sha256"])
        verify_recorded_hash(paths["summary"], review["strain_summary_sha256"])
        paths["review"] = review_path
    displacements = read_csv(paths["displacement"])
    strains = read_csv(paths["strain"])
    require(len(displacements) == summary["nodes_per_subcase"], f"Displacement count differs: {folder}")
    require(len(strains) == 3 * summary["shells_per_subcase"], f"Strain count differs: {folder}")
    displacements.sort(key=lambda row: int(row["grid_id"]))
    strains.sort(key=lambda row: (int(row["element_id"]), row["surface"]))
    node_ids = [int(r["grid_id"]) for r in displacements]
    strain_keys = [(int(r["element_id"]), r["surface"]) for r in strains]
    require(len(node_ids) == len(set(node_ids)) and len(strain_keys) == len(set(strain_keys)),
            f"Duplicate exported result keys: {folder}")
    disp = np.array([[float(r[key]) for key in DISPLACEMENT_COMPONENTS] for r in displacements])
    strain = np.array([[float(r[key]) for key in STRAIN_COMPONENTS] for r in strains])
    principal = np.array([[float(r[key]) for key in ("principal_max_microstrain", "principal_min_microstrain")]
                          for r in strains])
    xyz = np.array([[float(r[key]) for key in ("x_m", "y_m", "z_m")] for r in displacements])
    require(all(np.all(np.isfinite(values)) for values in (disp, strain, principal, xyz)), f"Nonfinite result in {folder}")
    strain_case, mechanics = summary["subcases"][0], reason["subcases"][0]
    require(strain_case["subcase"] == mechanics["subcase"], f"Subcase mismatch: {folder}")
    require({int(r["subcase"]) for r in displacements + strains} == {strain_case["subcase"]},
            f"CSV subcase mismatch: {folder}")
    fiber_surface = max(("top", "bottom"), key=lambda surface: strain_case["surfaces"][surface]["max_abs_principal_microstrain"])
    metrics = dict(alpha_deg=alpha, part=part, body_up_lift_N=mechanics["body_up_lift_N"],
        body_up_CL=mechanics["body_up_CL"], lateral_CAD_Fz_N=mechanics["aerodynamic_force_CAD_N"][2],
        max_translation_mm=mechanics["max_translation_m"] * 1e3,
        max_translation_grid_id=strain_case["max_translation_grid_id"],
        max_rotation_deg=mechanics["max_rotation_deg"],
        top_peak_abs_principal_microstrain=strain_case["surfaces"]["top"]["max_abs_principal_microstrain"],
        bottom_peak_abs_principal_microstrain=strain_case["surfaces"]["bottom"]["max_abs_principal_microstrain"],
        mid_peak_abs_principal_microstrain=strain_case["surfaces"]["mid"]["max_abs_principal_microstrain"],
        max_fiber_abs_principal_microstrain=strain_case["surfaces"][fiber_surface]["max_abs_principal_microstrain"],
        max_fiber_surface=fiber_surface, max_fiber_element_id=strain_case["surfaces"][fiber_surface]["max_abs_element_id"],
        peak_von_mises_MPa=mechanics["peak_von_mises_MPa"], yield_utilization=mechanics["yield_utilization"],
        force_balance_relative=mechanics["force_balance_relative"], moment_balance_relative=mechanics["moment_balance_relative"],
        max_abs_root_dof=mechanics["max_abs_root_dof"], coarse_screens_pass=reason["all_coarse_screens_pass"],
        accepted_for_workflow_demo=review.get("accepted_for_demo", False) if review else None,
        diagnostic_only=bool(summary.get("diagnostic_only") or reason.get("diagnostic_only")))
    return dict(folder=folder, alpha=actual_alpha, part=part, summary=summary, reason=reason, metrics=metrics,
                node_ids=node_ids, strain_keys=strain_keys, xyz=xyz, displacements=disp,
                strains=strain, principal=principal, hash_checks=hash_checks,
                source_files=[dict(path=str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path),
                                   sha256=digest(path)) for path in paths.values()])


def affine_check(angles, values, labels, unit, absolute_tolerance, relative_tolerance=5e-6, roundoff_scale_floor=None):
    """Fit y=a*alpha+b independently for every physical component."""
    values = np.asarray(values, dtype=float)
    require(values.ndim == 2 and values.shape[0] == len(angles), "Invalid affine-check data")
    design = np.column_stack((angles, np.ones(len(angles))))
    coefficients = np.linalg.lstsq(design, values, rcond=None)[0]
    residual = values - design @ coefficients
    largest = np.max(abs(residual), axis=0)
    scale = np.max(abs(values), axis=0)
    if roundoff_scale_floor is not None:
        scale = np.maximum(scale, np.asarray(roundoff_scale_floor))
    tolerance = absolute_tolerance + relative_tolerance * scale
    normalized = largest / tolerance
    worst = int(np.argmax(normalized))
    return dict(form="component(alpha_deg) = slope * alpha_deg + intercept", unit=unit,
        scalar_components_checked=values.shape[1], components_outside_tolerance=int(np.sum(normalized > 1)),
        passed=bool(np.all(normalized <= 1)), absolute_tolerance=absolute_tolerance,
        relative_tolerance=relative_tolerance, max_abs_fit_residual=float(largest.max()),
        max_normalized_residual=float(normalized[worst]), worst_component=labels[worst],
        worst_component_slope_per_deg=float(coefficients[0, worst]),
        worst_component_intercept_at_zero_deg=float(coefficients[1, worst]),
        max_abs_intercept_at_zero_deg=float(np.max(abs(coefficients[1]))),
        roundoff_scale_floor_applied=roundoff_scale_floor is not None,
        note="Nonzero intercept is allowed because initial camber/twist W2GJ remains active. No fit through zero is imposed.")


def compare_parts(first, second):
    require(first["node_ids"] == second["node_ids"] and first["strain_keys"] == second["strain_keys"],
            "Mirrored parts do not have matching structural result IDs")
    reflected_xyz = first["xyz"].copy()
    reflected_xyz[:, 2] = 1.56 - reflected_xyz[:, 2]
    geometry_difference = float(np.max(np.linalg.norm(second["xyz"] - reflected_xyz, axis=1)))
    require(geometry_difference < 2e-6, "Parts do not match the documented CAD reflection Z2=1.56-Z1")
    # Displacement/force are polar vectors; rotation/moment are axial vectors.
    translation_difference = second["displacements"][:, :3] - first["displacements"][:, :3] * [1, 1, -1]
    rotation_difference = second["displacements"][:, 3:] - first["displacements"][:, 3:] * [-1, -1, 1]
    principal_difference = abs(second["principal"] - first["principal"])
    mechanics1, mechanics2 = (case["reason"]["subcases"][0] for case in (first, second))
    force_difference = np.array(mechanics2["aerodynamic_force_CAD_N"]) - np.array(mechanics1["aerodynamic_force_CAD_N"]) * [1, 1, -1]
    moment_difference = np.array(mechanics2["aerodynamic_moment_root_origin_Nm"]) - np.array(mechanics1["aerodynamic_moment_root_origin_Nm"]) * [-1, -1, 1]
    return dict(alpha_deg=first["metrics"]["alpha_deg"], max_geometry_reflection_error_m=geometry_difference,
        max_reflected_translation_vector_difference_mm=float(np.max(np.linalg.norm(translation_difference, axis=1))) * 1e3,
        max_reflected_rotation_vector_difference_rad=float(np.max(np.linalg.norm(rotation_difference, axis=1))),
        max_matched_principal_strain_difference_microstrain=float(principal_difference.max()),
        reflected_aerodynamic_force_difference_norm_N=float(np.linalg.norm(force_difference)),
        reflected_root_moment_difference_norm_Nm=float(np.linalg.norm(moment_difference)),
        lift_difference_N=abs(second["metrics"]["body_up_lift_N"] - first["metrics"]["body_up_lift_N"]),
        max_translation_peak_difference_mm=abs(second["metrics"]["max_translation_mm"] - first["metrics"]["max_translation_mm"]),
        max_fiber_principal_peak_difference_microstrain=abs(second["metrics"]["max_fiber_abs_principal_microstrain"] - first["metrics"]["max_fiber_abs_principal_microstrain"]),
        peak_von_mises_difference_MPa=abs(second["metrics"]["peak_von_mises_MPa"] - first["metrics"]["peak_von_mises_MPa"]))


def trend(values):
    differences = np.diff(values)
    tolerance = 1e-9 * max(1., max(abs(value) for value in values))
    if np.all(differences >= -tolerance):
        return "nondecreasing"
    if np.all(differences <= tolerance):
        return "nonincreasing"
    return "mixed"


def write_comparison_plot(path, metrics, has_diagnostics):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 7.7), constrained_layout=True)
    specifications = (("body_up_lift_N", "Body-up aerodynamic force [N]"),
                      ("max_translation_mm", "Maximum displacement magnitude [mm]"),
                      ("max_fiber_abs_principal_microstrain", "Peak |fiber principal strain| [microstrain]"),
                      ("peak_von_mises_MPa", "Peak reconstructed von Mises stress [MPa]"))
    for ax, (key, ylabel) in zip(axes.ravel(), specifications):
        for part, color, style, marker in ((1, "#147d92", "-", "o"), (2, "#c46d24", "--", "x")):
            selected = [row for row in metrics if row["part"] == part]
            ax.plot([r["alpha_deg"] for r in selected], [r[key] for r in selected],
                    color=color, linestyle=style, marker=marker, markersize=6,
                    linewidth=1.6, label=f"Part {part}")
        ax.set(xlabel="Angle of attack [deg]", ylabel=ylabel, xticks=ANGLES)
        ax.grid(alpha=.22)
        ax.legend(frameon=False)
    subtitle = "Fixed geometry, dynamic pressure and W2GJ; points are independent solved jobs"
    if has_diagnostics:
        subtitle += "\nRuntime diagnostics retained; see README"
    fig.suptitle("BWB | Angle-of-attack sweep\n" + subtitle, fontsize=13)
    fig.savefig(path, dpi=200, facecolor="white")
    plt.close(fig)


def summarize(results=DEFAULT_RESULTS, baseline=None):
    results = Path(results).resolve()
    baseline = Path(baseline).resolve() if baseline else results.parent / "sol144_demo"
    cases = {}
    # Complete collection and verification precede all aggregate output writes.
    for angle in ANGLES:
        for part in (1, 2):
            folder = (baseline if angle == 2 else results / f"aoa_{angle:02}") / f"part_{part}"
            cases[(angle, part)] = load_case(folder, angle, part)
    dynamic_pressures = [case["reason"]["rough_finite_wing"]["dynamic_pressure_Pa"] for case in cases.values()]
    require(np.allclose(dynamic_pressures, dynamic_pressures[0], atol=1e-8, rtol=0), "Dynamic pressure changed during sweep")
    affine, trends = {}, {}
    for part in (1, 2):
        selected = [cases[(angle, part)] for angle in ANGLES]
        reference = selected[0]
        for case in selected[1:]:
            require(case["node_ids"] == reference["node_ids"] and case["strain_keys"] == reference["strain_keys"],
                    f"Structural result indexing changed for part {part}")
            for name in ("structural_cards.bdf", "aero_cards.bdf", "spline_cards.bdf", "dmicards.bdf"):
                require(case["summary"]["input_sha256"][name] == reference["summary"]["input_sha256"][name],
                        f"Fixed model input changed for part {part}: {name}")
        angles = [case["alpha"] for case in selected]
        transl_labels = [f"GRID {nid} {component}" for nid in reference["node_ids"] for component in DISPLACEMENT_COMPONENTS[:3]]
        rot_labels = [f"GRID {nid} {component}" for nid in reference["node_ids"] for component in DISPLACEMENT_COMPONENTS[3:]]
        strain_labels = [f"element {eid} {surface} {component}" for eid, surface in reference["strain_keys"] for component in STRAIN_COMPONENTS]
        # Mid-components average printed top/bottom values. When those cancel,
        # their roundoff is set by the parent magnitudes, not the near-zero mean.
        all_strains = np.array([case["strains"] for case in selected])
        strain_roundoff_scale = np.max(abs(all_strains), axis=0)
        strain_key_index = {key: index for index, key in enumerate(reference["strain_keys"])}
        for (eid, surface), index in strain_key_index.items():
            if surface == "mid":
                strain_roundoff_scale[index] = np.maximum(strain_roundoff_scale[strain_key_index[(eid, "top")]],
                                                          strain_roundoff_scale[strain_key_index[(eid, "bottom")]])
        affine[str(part)] = dict(
            translations=affine_check(angles, np.array([c["displacements"][:, :3].ravel() for c in selected]), transl_labels, "m", 1e-12),
            rotations=affine_check(angles, np.array([c["displacements"][:, 3:].ravel() for c in selected]), rot_labels, "rad", 1e-12),
            local_strain_components=affine_check(angles, all_strains.reshape(len(selected), -1), strain_labels, "dimensionless", 1e-13,
                                                 roundoff_scale_floor=strain_roundoff_scale.ravel()),
            aerodynamic_force=affine_check(angles, np.array([c["reason"]["subcases"][0]["aerodynamic_force_CAD_N"] for c in selected]),
                                           ["CAD Fx", "CAD Fy", "CAD Fz"], "N", 1e-7),
            aerodynamic_root_moment=affine_check(angles, np.array([c["reason"]["subcases"][0]["aerodynamic_moment_root_origin_Nm"] for c in selected]),
                                                 ["CAD Mx", "CAD My", "CAD Mz"], "N*m", 1e-7))
        trends[str(part)] = {key: dict(values=[c["metrics"][key] for c in selected],
                                      direction=trend([c["metrics"][key] for c in selected]))
                            for key in ("body_up_lift_N", "max_translation_mm", "max_fiber_abs_principal_microstrain", "peak_von_mises_MPa")}
    metrics = [cases[(angle, part)]["metrics"] for angle in ANGLES for part in (1, 2)]
    mirror = [compare_parts(cases[(angle, 1)], cases[(angle, 2)]) for angle in ANGLES]
    diagnostic_cases = []
    for (angle, part), case in cases.items():
        if case["metrics"]["diagnostic_only"]:
            diagnostic_cases.append(dict(alpha_deg=angle, part=part, diagnostics=case["summary"].get("diagnostics", {})))
    affine_pass = all(check["passed"] for part in affine.values() for check in part.values())
    report = dict(status="summary_with_runtime_diagnostics" if diagnostic_cases else "summarized",
        angles_deg=list(ANGLES), new_angles_deg=[3, 4, 5, 6], parts=[1, 2], total_jobs=10, new_jobs=8,
        baseline=str(baseline), dynamic_pressure_Pa=dynamic_pressures[0],
        all_coarse_screens_pass=all(row["coarse_screens_pass"] for row in metrics),
        all_new_jobs_accepted_for_workflow_demo=all(row["accepted_for_workflow_demo"] for row in metrics if row["alpha_deg"] > 2),
        all_component_affine_checks_pass=affine_pass, component_affine_checks=affine,
        trends=trends, mirror_comparison=mirror,
        maximum_mirror_differences={key: max(row[key] for row in mirror) for key in mirror[0] if key != "alpha_deg"},
        cases=metrics, diagnostic_cases=diagnostic_cases,
        provenance=[dict(alpha_deg=angle, part=part, hash_checks=case["hash_checks"], source_files=case["source_files"],
                         source_f06_sha256=case["summary"]["source_f06_sha256"]) for (angle, part), case in cases.items()],
        mirror_convention="CAD Z2=1.56-Z1. Polar vector sign transform (1,1,-1); axial vector sign transform (-1,-1,1). Same shell IDs/surfaces compared using principal-strain invariants, not unlike local axes.",
        limitations=["Maximal principal strains, stress extrema and displacement magnitudes are not required to be affine.",
                     "Within-part local strain components and GRID DOFs are fit with a free intercept; W2GJ is retained.",
                     "All results retain the underlying solver runtime diagnostics and assumed material/flight conditions.",
                     "A sweep consistency check is not grid convergence, experimental validation or certification."])
    results.mkdir(parents=True, exist_ok=True)
    for filename, rows in (("aoa_summary.csv", metrics), ("mirror_comparison.csv", mirror)):
        with (results / filename).open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    (results / "aoa_summary.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    write_comparison_plot(results / "aoa_comparison.png", metrics, bool(diagnostic_cases))
    table = ["| 迎角 | 向上气动力 1 / 2（N） | 最大位移 1 / 2（mm） | 最大纤维主应变绝对值 1 / 2（µε） | 峰值 von Mises 1 / 2（MPa） |",
             "|---:|---:|---:|---:|---:|"]
    for angle in ANGLES:
        pair = [cases[(angle, part)]["metrics"] for part in (1, 2)]
        entries = [f"{pair[0][key]:.6g} / {pair[1][key]:.6g}" for key in
                   ("body_up_lift_N", "max_translation_mm", "max_fiber_abs_principal_microstrain", "peak_von_mises_MPa")]
        table.append(f"| {angle}° | " + " | ".join(entries) + " |")
    maxima = report["maximum_mirror_differences"]
    lines = ["# BWB 临时迎角试算", "",
        "试算结果存放于 `workspace/BWB/trials/results/aoa_sweep/`，不属于正式模型成果。", "",
        "两个零件分别计算 3°、4°、5°、6°，并保留此前 2° 结果作为基线，共 10 组结果。各组是独立 SOL144 作业；根部固定、材料、结构、气动面板、SPLINE、W2GJ 和动压保持一致。",
        f"动压为 {dynamic_pressures[0]:g} Pa。此处向上气动力为各独立零件的 CAD +Y 气动力；该方向与严格的风轴升力方向存在迎角差异。", "",
        *table, "", "![迎角比较](aoa_comparison.png)", "",
        f"逐节点平移、转角、逐单元三个表面的局部 exx/eyy/γxy 以及总气动力/根部力矩，均按 `a × α + b` 检查：{'全部通过' if affine_pass else '存在超出打印精度容差的分量，见 JSON'}。拟合保留截距，因为 W2GJ 包含初始弯度/扭转。最大位移模长、主应变和应力极值只汇总趋势，不强制线性或穿过原点。", "",
        "mid 分量由已打印的上下纤维分量平均得到；当上下分量相消时，其打印舍入容差按母纤维幅值传播，不能按接近零的平均值缩小。", "",
        f"镜像变换后，10 组中的最大同节点位移差为 {maxima['max_reflected_translation_vector_difference_mm']:.6g} mm；同单元主应变最大差为 {maxima['max_matched_principal_strain_difference_microstrain']:.6g} µε；两件向上气动力最大差为 {maxima['lift_difference_N']:.6g} N。位移和力使用极向量变换，转角和力矩使用轴向量变换；主应变比较不受局部单元轴差异影响。", "",
        f"各工况的力/力矩平衡、根部约束、小变形及假设屈服强度粗筛：{'全部通过' if report['all_coarse_screens_pass'] else '至少一项需复查'}。详细值见各零件目录的 `reasonableness.json`。"]
    six_top = cases[(6, 1)]["summary"]["subcases"][0]["surfaces"]["top"]["max_principal_microstrain"]
    six_abs = cases[(6, 1)]["metrics"]["max_fiber_abs_principal_microstrain"]
    lines += ["", "Femap 模型：", "",
              "- [零件 1：3°–6°结果](part_1_aoa_3to6.modfem)",
              "- [零件 2：3°–6°结果](part_2_aoa_3to6.modfem)", "",
              f"每件 Output Set 1–4 对应 3°、4°、5°、6°，默认显示 6°。Femap 默认 7071 是 top-major（正纤维侧最大主应变），6° 峰值约 {six_top:.6g} µε；表中的 {six_abs:.6g} µε 是跨上下两纤维表面的最大主应变绝对值，来自 bottom 压应变。两者统计量不同。"]
    if diagnostic_cases:
        errors = sorted({message for case in diagnostic_cases for key in ("runtime_error_lines", "recorded_runtime_errors")
                         for message in case["diagnostics"].get(key, [])})
        accepted_text = "8 个新工况的 case_review 均接受用于流程演示。" if report["all_new_jobs_accepted_for_workflow_demo"] else "个别新工况尚无接受用于流程演示的记录。"
        lines += ["", accepted_text,
                  f"**运行诊断仍保留：{len(diagnostic_cases)} 组带有原始运行错误记录。** 数据完成标记和结果表完整，本汇总保留已审阅的诊断信息，不将这些作业标记为无错误求解。", "",
                  *[f"- `{error}`" for error in errors]]
    lines += ["", "文件：`aoa_summary.csv` 为工况指标，`mirror_comparison.csv` 为两件差异，`aoa_summary.json` 包含分量仿射性、输入校验和来源记录。这里只检查当前假设下的结果一致性与量级，不替代网格收敛、试验验证或强度认证。", "",
              "复现汇总：", "", "```powershell", "python BWB/scripts/summarize_aoa_sweep.py", "```", ""]
    (results / "README.md").write_text("\n".join(lines), encoding="utf-8")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--baseline", type=Path)
    args = parser.parse_args()
    try:
        report = summarize(args.results, args.baseline)
    except (ValueError, OSError, KeyError) as exc:
        parser.exit(1, f"Sweep summary rejected: {exc}\n")
    print(json.dumps({key: report[key] for key in ("status", "total_jobs", "new_jobs", "all_coarse_screens_pass",
                                                   "all_component_affine_checks_pass", "maximum_mirror_differences")}, indent=2))


if __name__ == "__main__":
    main()
