"""Compare completed 200 m/s SOL144 jobs with the actual 20 m/s AOA sweep.

Reads existing result tables only. Physical-screen failures are reported as
model-applicability limitations, independently of result/solver completion.
Nothing is inferred by multiplying the low-speed response by 100.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

try:
    from .inspect_shell_bdf import cards, nastran_float
    from .summarize_aoa_sweep import load_case, require, verify_recorded_hash
except ImportError:
    from inspect_shell_bdf import cards, nastran_float
    from summarize_aoa_sweep import load_case, require, verify_recorded_hash

ROOT = Path(__file__).resolve().parents[2]
ANGLES = (3, 4, 5, 6)
MODEL_FILES = ("structural_cards.bdf", "aero_cards.bdf", "spline_cards.bdf", "dmicards.bdf")
YIELD_MPA = 200.0


def source_record(path):
    path = Path(path)
    return dict(path=str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path),
                sha256=hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest(),
                hash_mode="crlf_to_lf")


def verify_source_record(record):
    path = ROOT / record["path"]
    mode = record.get("hash_mode", "raw")
    require(mode in ("raw", "crlf_to_lf"), f"Unknown source hash mode: {mode}")
    raw = path.read_bytes()
    if mode == "crlf_to_lf":
        raw = raw.replace(b"\r\n", b"\n")
    require(hashlib.sha256(raw).hexdigest() == record["sha256"], f"Review source hash differs: {path}")
    return path


def read_flight(case, input_dir, expected_speed):
    path = input_dir / "flight_condition.bdf"
    verify_recorded_hash(path, case["summary"]["input_sha256"][path.name])
    trim = [fields for name, fields, _ in cards(path) if name == "TRIM"]
    require(len(trim) == 1, f"Expected one TRIM card: {path}")
    fields = trim[0]
    mach, q = nastran_float(fields[2]), nastran_float(fields[3])
    alpha = math.degrees(nastran_float(fields[fields.index("ANGLEA") + 1]))
    manifest_path = input_dir / "case_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
    flight = manifest["flight"]
    speed, density, sound = (flight[key] for key in ("velocity_m_s", "density_kg_m3", "sound_speed_m_s"))
    require(math.isclose(speed, expected_speed, abs_tol=1e-8), f"Wrong speed: {input_dir}")
    require(math.isclose(mach, speed / sound, rel_tol=1e-8), f"TRIM Mach differs from manifest: {input_dir}")
    require(math.isclose(q, .5 * density * speed * speed, rel_tol=1e-8), f"TRIM q differs from manifest: {input_dir}")
    require(math.isclose(q, case["reason"]["rough_finite_wing"]["dynamic_pressure_Pa"], rel_tol=1e-8),
            f"Result q differs from input: {input_dir}")
    require(abs(alpha - case["alpha"]) < 1e-7, f"Result AOA differs from input: {input_dir}")
    model_hashes = {}
    for name in MODEL_FILES:
        verify_recorded_hash(input_dir / name, case["summary"]["input_sha256"][name])
        model_hashes[name] = source_record(input_dir / name)["sha256"]
    return dict(speed_m_s=speed, mach=mach, dynamic_pressure_Pa=q, density_kg_m3=density,
                sound_speed_m_s=sound, model_hashes=model_hashes,
                provenance=[source_record(path), source_record(manifest_path)])


def stress_statistics(case):
    """Count unique shell centroids whose top or bottom VM exceeds 200 MPa."""
    material = case["reason"]["material"]
    young, poisson = material["young_modulus_Pa"], material["poisson_ratio"]
    require(young > 0 and -1 < poisson < .5, "Invalid isotropic material")
    require(math.isclose(material["assumed_yield_MPa"], YIELD_MPA, rel_tol=1e-8),
            "The requested 200 MPa threshold differs from the case material screen")
    fiber = np.array([surface in ("top", "bottom") for _, surface in case["strain_keys"]])
    keys = [key for key, selected in zip(case["strain_keys"], fiber) if selected]
    strain = case["strains"][fiber]
    exx, eyy, gamma = strain.T
    factor = young / (1 - poisson * poisson)
    sx, sy = factor * (exx + poisson * eyy), factor * (eyy + poisson * exx)
    tau = young / (2 * (1 + poisson)) * gamma
    vm = np.sqrt(np.maximum(0., sx * sx - sx * sy + sy * sy + 3 * tau * tau)) / 1e6
    require(np.all(np.isfinite(vm)), "Nonfinite reconstructed stress")
    peak = int(np.argmax(vm))
    expected = case["reason"]["subcases"][0]["peak_von_mises_MPa"]
    require(math.isclose(float(vm[peak]), expected, rel_tol=1e-8, abs_tol=1e-9),
            f"Reconstructed stress differs from reasonableness report: {case['folder']}")
    exceed_ids = sorted({eid for (eid, _), value in zip(keys, vm) if value > YIELD_MPA})
    return dict(peak_von_mises_MPa=float(vm[peak]), peak_stress_element_id=keys[peak][0],
                peak_stress_surface=keys[peak][1], elements_above_200MPa=len(exceed_ids),
                fiber_samples_above_200MPa=int(np.sum(vm > YIELD_MPA)),
                element_ids_above_200MPa=exceed_ids,
                elements_checked=len({eid for eid, _ in keys}),
                threshold_MPa=YIELD_MPA, young_modulus_Pa=young, poisson_ratio=poisson)


def collect_case(folder, input_dir, angle, part, speed):
    case = load_case(folder, angle, part)
    flight = read_flight(case, input_dir, speed)
    stresses = stress_statistics(case)
    mechanics = case["reason"]["subcases"][0]
    screens = mechanics["screens"]
    review_path = folder / "case_review.json"
    review = json.loads(review_path.read_text(encoding="utf-8-sig")) if review_path.is_file() else {}
    if speed == 200:
        require("runtime_diagnostic_review" in review, f"New-speed runtime review is missing: {folder}")
        verify_source_record(review["runtime_diagnostic_review"])
    m = case["metrics"]
    applicable = all(screens[key] for key in ("small_translation", "small_rotation", "below_assumed_yield"))
    row = dict(speed_m_s=speed, alpha_deg=angle, part=part, mach=flight["mach"],
        dynamic_pressure_Pa=flight["dynamic_pressure_Pa"], body_up_force_N=m["body_up_lift_N"],
        body_up_force_coefficient=m["body_up_CL"], max_translation_mm=m["max_translation_mm"],
        max_translation_grid_id=m["max_translation_grid_id"], max_rotation_deg=m["max_rotation_deg"],
        max_translation_over_reference_span=mechanics["max_translation_over_reference_span"],
        max_fiber_abs_principal_microstrain=m["max_fiber_abs_principal_microstrain"],
        peak_strain_element_id=m["max_fiber_element_id"], peak_strain_surface=m["max_fiber_surface"],
        **{key: value for key, value in stresses.items() if key != "element_ids_above_200MPa"},
        yield_utilization=stresses["peak_von_mises_MPa"] / YIELD_MPA,
        small_translation_screen=bool(screens["small_translation"]),
        small_rotation_screen=bool(screens["small_rotation"]),
        below_assumed_yield_screen=bool(screens["below_assumed_yield"]),
        within_small_deflection_and_elastic_screens=applicable,
        interpretation="within_coarse_model_screens" if applicable else "linear_solution_beyond_model_screens",
        force_balance_screen=bool(screens["force_balance"]), moment_balance_screen=bool(screens["moment_balance"]),
        root_fixed_screen=bool(screens["root_fixed"]), force_balance_relative=m["force_balance_relative"],
        moment_balance_relative=m["moment_balance_relative"], max_abs_root_dof=m["max_abs_root_dof"],
        solver_end_of_job=case["summary"]["solver_end_of_job"],
        raw_solver_status=review.get("raw_solver_status", case["summary"]["status"]),
        runtime_diagnostics_retained=m["diagnostic_only"],
        accepted_for_workflow_demo=review.get("accepted_for_demo"), review_status=review.get("status"))
    require(row["below_assumed_yield_screen"] == (row["peak_von_mises_MPa"] < YIELD_MPA),
            f"Yield-screen disagreement: {folder}")
    case.update(flight=flight, stress_statistics=stresses, row=row, review=review)
    return case


def compare_speed(low, high):
    require(low["flight"]["model_hashes"] == high["flight"]["model_hashes"],
            f"Structural/aero/SPLINE/W2GJ input changed between speeds for part {low['part']}")
    require(low["node_ids"] == high["node_ids"] and low["strain_keys"] == high["strain_keys"],
            "Structural result indexing differs between speeds")
    require(np.allclose(low["xyz"], high["xyz"], atol=1e-12, rtol=0), "Structural coordinates differ between speeds")
    row = dict(alpha_deg=high["row"]["alpha_deg"], part=high["part"],
               dynamic_pressure_ratio=high["row"]["dynamic_pressure_Pa"] / low["row"]["dynamic_pressure_Pa"])
    for key in ("body_up_force_N", "max_translation_mm", "max_fiber_abs_principal_microstrain", "peak_von_mises_MPa"):
        a, b = low["row"][key], high["row"][key]
        row[f"{key}_20m_s"], row[f"{key}_200m_s"] = a, b
        row[f"{key}_actual_ratio"] = b / a if a != 0 else None
    return row


def write_plot(path, rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8), constrained_layout=True)
    specs = (("body_up_force_N", "Body-up aerodynamic force [N]"),
             ("max_translation_mm", "Maximum displacement [mm]"),
             ("max_fiber_abs_principal_microstrain", "Peak |fiber principal strain| [microstrain]"),
             ("peak_von_mises_MPa", "Peak reconstructed von Mises stress [MPa]"))
    for ax, (key, label) in zip(axes.ravel(), specs):
        for speed, color in ((20, "#177d91"), (200, "#b3472d")):
            for part, style, marker in ((1, "-", "o"), (2, "--", "x")):
                data = [r for r in rows if r["speed_m_s"] == speed and r["part"] == part]
                ax.plot([r["alpha_deg"] for r in data], [r[key] for r in data], color=color,
                        linestyle=style, marker=marker, markersize=5, label=f"{speed} m/s, part {part}")
        # A logarithmic response axis keeps both actual speeds readable.
        if all(r[key] > 0 for r in rows):
            ax.set_yscale("log")
        else:
            ax.set_yscale("symlog", linthresh=max(1e-9, max(abs(r[key]) for r in rows) * 1e-3))
        if key == "peak_von_mises_MPa":
            ax.axhline(YIELD_MPA, color="0.35", linewidth=1, linestyle=":", label="Assumed yield: 200 MPa")
        ax.set(xlabel="Angle of attack [deg]", ylabel=label, xticks=ANGLES)
        ax.grid(alpha=.2, which="both")
        ax.legend(frameon=False, fontsize=8)
    fig.suptitle("BWB | Actual solved responses at 20 and 200 m/s\n"
                 "Linear model; applicability flags and runtime diagnostics are retained", fontsize=13)
    fig.savefig(path, dpi=200, facecolor="white")
    plt.close(fig)


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_readme(path, report):
    high = [r for r in report["cases"] if r["speed_m_s"] == 200]
    table = ["| 迎角 | 零件 | 向上气动力 N | 最大位移 mm | 最大绝对主应变 µε | 峰值 VM MPa | >200 MPa 单元数 |",
             "|---:|---:|---:|---:|---:|---:|---:|"]
    for r in high:
        table.append(f"| {r['alpha_deg']}° | {r['part']} | {r['body_up_force_N']:.6g} | {r['max_translation_mm']:.6g} | "
                     f"{r['max_fiber_abs_principal_microstrain']:.6g} | {r['peak_von_mises_MPa']:.6g} | {r['elements_above_200MPa']} |")
    limits = ["| 迎角 / 零件 | 位移/展长 | 最大转角 ° | 小位移 / 小转角 / 弹性粗筛 | 应变峰值单元/侧 | 应力峰值单元/侧 |",
              "|---:|---:|---:|---|---|---|"]
    for r in high:
        flags = " / ".join("通过" if r[key] else "超出" for key in
                           ("small_translation_screen", "small_rotation_screen", "below_assumed_yield_screen"))
        limits.append(f"| {r['alpha_deg']}° / {r['part']} | {r['max_translation_over_reference_span']:.5g} | "
                      f"{r['max_rotation_deg']:.5g} | {flags} | {r['peak_strain_element_id']} / {r['peak_strain_surface']} | "
                      f"{r['peak_stress_element_id']} / {r['peak_stress_surface']} |")
    ratio_table = ["| 迎角 / 零件 | 向上气动力比 | 最大位移比 | 最大绝对主应变比 | 峰值应力比 |",
                   "|---:|---:|---:|---:|---:|"]
    for r in report["speed_comparison"]:
        values = [r[f"{key}_actual_ratio"] for key in
                  ("body_up_force_N", "max_translation_mm", "max_fiber_abs_principal_microstrain", "peak_von_mises_MPa")]
        ratio_table.append(f"| {r['alpha_deg']}° / {r['part']} | " +
                           " | ".join("未定义" if value is None else f"{value:.5g}" for value in values) + " |")
    restricted = sum(not r["within_small_deflection_and_elastic_screens"] for r in high)
    lines = ["# BWB：200 m/s 与 20 m/s 求解结果比较", "",
        f"200 m/s 工况使用 Mach {high[0]['mach']:.7f}、动压 {high[0]['dynamic_pressure_Pa']:.0f} Pa；"
        "密度 1.225 kg/m³、声速 340.3 m/s。两个零件各自固定根部，分别计算 3°、4°、5°、6°。"
        "结构仍为假设铝材、1 mm 蒙皮，无 rib/spar。", "",
        "**8 组均超过假设 200 MPa 屈服及小变形筛查界限，以下为超出适用范围的线性输出。**", "",
        "本汇总读取两个速度各自的实际 SOL144 输出。动压增加 100 倍，但 Mach 和气动弹性反馈也改变；"
        "表中响应没有按 100 倍推算。向上气动力是 CAD +Y 分量，方向与风轴升力存在迎角差异。", "",
        *table, "", "![实际速度比较](speed_comparison.png)", "",
        "图的响应纵轴使用对数刻度，以同时显示两个速度；应力图虚线为假设的 200 MPa 屈服强度。", "",
        f"200 m/s 的 {restricted}/8 组结果超出至少一项小变形/弹性粗筛。"
        "粗筛采用最大位移/参考展长 <1%、最大转角 <0.05 rad（约 2.865°）以及峰值 VM <200 MPa。"
        "这些是当前模型的适用性检查，超出并不等同于求解作业失败，通过也不证明结构安全。", "",
        *limits, "",
        "最大主应变统计跨 top/bottom 两个纤维面的绝对值；top/bottom 指单元局部 ±z 侧，"
        "不等同于几何上蒙皮/下蒙皮。von Mises 由实际局部 exx、eyy、工程剪应变重构平面应力。"
        "超限数按任一 top/bottom 纤维超过 200 MPa 的唯一单元计数，仅代表单元中心采样，"
        "不代表非线性屈服区面积或实际塑性预测。", "",
        "下表为实际结果的 200 m/s / 20 m/s 比值，分别比较各自峰值；峰值位置可能变化。", "",
        *ratio_table, "",
        "线性 SOL144 不包含大变形几何非线性、塑性、屈曲或失速。超过屈服的线性应力只是"
        "弹性外推，超过小变形范围的位移也不能直接当作真实变形。Mach 约 0.588 的当前气动模型"
        "仍需确认局部压缩性/流动适用性；这里没有补做黏性或激波 CFD。", "",
        "[稳定性审查](stability_review.json) 确认存在非零气动刚度并完成有效刚度矩阵分解，"
        "但未进行 DIVERG 特征值分析，没有计算临界发散动压或稳定裕度。"
        "静态线性方程有数值解不能证明 200 m/s 低于发散速度；大幅响应也不能单独证明已经发散。", "",
        f"运行记录：{sum(r['solver_end_of_job'] for r in high)}/8 组有 END OF JOB，"
        f"{sum(r['runtime_diagnostics_retained'] for r in high)}/8 组保留原始运行诊断。"
        "每组 case_review 的原始状态与适用性检查分别保存，未将带诊断的作业改写为无错误求解。", "",
        "[运行诊断审查](runtime_diagnostic_review.json) 在 200 m/s、6°、零件 1 直接跟踪到 "
        "DGEMM 的零输出行调用（m=0、ldc=0），普通运行和跟踪运行的三份导出 CSV 逐字节一致。"
        "其余七组共享同一诊断特征，未逐组跟踪。原始 MKL 错误和 launcher 返回码 1 继续保留；"
        "此数值复现性检查不消除上述物理适用性限制。", "",
        "`speed_summary.csv` 包含两速度共 16 组的峰值、单元、约束、平衡和适用性标记；"
        "`speed_comparison.csv` 保存实际响应比；`speed_summary.json` 保存超限单元 ID、输入/结果来源和审查记录。", "",
        "Femap 结果文件：", "",
        "- [零件 1：200 m/s，3°–6°](part_1_200mps_aoa_3to6.modfem)",
        "- [零件 2：200 m/s，3°–6°](part_2_200mps_aoa_3to6.modfem)", "",
        "每件有 4 个 Output Set，依次对应 3°、4°、5°、6°，均可切换；默认显示 6° 的 bottom von Mises 应力。"
        "Femap 应力单位为 Pa，显示数值除以 1e6 得到 MPa。", "",
        "复现汇总：", "", "```powershell", "python BWB/scripts/summarize_speed_200.py", "```", ""]
    path.write_text("\n".join(lines), encoding="utf-8")


def summarize(results, baseline, inputs, baseline_inputs, stability_review=None):
    cases = {}
    for speed, output_root, input_root in ((20, baseline, baseline_inputs), (200, results, inputs)):
        for angle in ANGLES:
            for part in (1, 2):
                relative = Path(f"aoa_{angle:02}") / f"part_{part}"
                cases[(speed, angle, part)] = collect_case(output_root / relative, input_root / relative, angle, part, speed)
    comparisons = [compare_speed(cases[(20, a, p)], cases[(200, a, p)]) for a in ANGLES for p in (1, 2)]
    rows = [case["row"] for case in cases.values()]
    high = [row for row in rows if row["speed_m_s"] == 200]
    stability_review = stability_review or results / "stability_review.json"
    stability_data = json.loads(stability_review.read_text(encoding="utf-8-sig"))
    reviewed_cases = {(entry["aoa_deg"], entry["part"]): entry for entry in stability_data["cases"]}
    require(set(reviewed_cases) == {(angle, part) for angle in ANGLES for part in (1, 2)},
            "Stability review does not cover all eight 200 m/s cases")
    for (angle, part), entry in reviewed_cases.items():
        require(entry["f06_sha256"] == cases[(200, angle, part)]["summary"]["source_f06_sha256"],
                f"Stability review F06 differs for {angle} degrees, part {part}")
    stability = dict(assessment="not_established_by_this_summary",
                     note="A completed static solve does not establish freedom from aeroelastic divergence.",
                     independent_review=stability_data, provenance=source_record(stability_review))
    runtime_path = results / "runtime_diagnostic_review.json"
    runtime_data = json.loads(runtime_path.read_text(encoding="utf-8-sig"))
    traced = cases[(200, runtime_data["angle_of_attack_deg"], runtime_data["part"])]
    require(runtime_data["f06_sha256"]["normal"] == traced["summary"]["source_f06_sha256"],
            "Runtime diagnostic review belongs to a different F06")
    for name, comparison in runtime_data["result_comparison"].items():
        require(comparison["normal_sha256"] == traced["summary"]["csv_sha256"][name],
                f"Runtime review result hash differs: {name}")
    report = dict(schema_version=1, status="summarized_completed_linear_solutions", angles_deg=list(ANGLES),
        speeds_m_s=[20, 200], new_jobs=8, baseline_jobs=8,
        all_200m_s_jobs_have_end_of_job=all(row["solver_end_of_job"] for row in high),
        all_200m_s_equilibrium_and_root_screens_pass=all(row[key] for row in high for key in
                                                       ("force_balance_screen", "moment_balance_screen", "root_fixed_screen")),
        all_200m_s_within_small_deflection_and_elastic_screens=all(row["within_small_deflection_and_elastic_screens"] for row in high),
        screen_thresholds=dict(max_translation_over_span=.01, max_rotation_rad=.05, assumed_yield_MPa=YIELD_MPA),
        cases=rows, speed_comparison=comparisons, aeroelastic_stability=stability,
        runtime_diagnostic_review=dict(review=runtime_data, provenance=source_record(runtime_path)),
        element_ids_above_200MPa=[dict(speed_m_s=s, alpha_deg=a, part=p,
            element_ids=case["stress_statistics"]["element_ids_above_200MPa"]) for (s, a, p), case in cases.items()],
        reviews=[dict(speed_m_s=s, alpha_deg=a, part=p, review=case["review"],
                      diagnostics=case["summary"].get("diagnostics", {})) for (s, a, p), case in cases.items()],
        provenance=[dict(speed_m_s=s, alpha_deg=a, part=p, source_f06_sha256=case["summary"]["source_f06_sha256"],
            result_files=case["source_files"], input_files=case["flight"]["provenance"],
            fixed_model_hashes=case["flight"]["model_hashes"], hash_checks=case["hash_checks"])
            for (s, a, p), case in cases.items()],
        limitations=["All ratios use actual independent solver results; no 100x response scaling is performed.",
                     "Small-deflection or elastic-screen exceedance is a physical interpretation limit, not a solver-completion flag.",
                     "Stress above the assumed yield is linear elastic extrapolation, not plastic stress prediction.",
                     "Centroid fiber checks do not prove mesh convergence, local buckling resistance or structural safety.",
                     "Static solution completion does not establish stability against divergence."])
    results.mkdir(parents=True, exist_ok=True)
    write_csv(results / "speed_summary.csv", rows)
    write_csv(results / "speed_comparison.csv", comparisons)
    (results / "speed_summary.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    write_plot(results / "speed_comparison.png", rows)
    write_readme(results / "README.md", report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=ROOT / "BWB/results/speed_200")
    parser.add_argument("--baseline", type=Path, default=ROOT / "BWB/results/aoa_sweep")
    parser.add_argument("--cases", type=Path, default=ROOT / "BWB/cases/speed_200")
    parser.add_argument("--baseline-cases", type=Path, default=ROOT / "BWB/cases/aoa_sweep")
    parser.add_argument("--stability-review", type=Path)
    args = parser.parse_args()
    try:
        report = summarize(args.results.resolve(), args.baseline.resolve(), args.cases.resolve(),
                           args.baseline_cases.resolve(), args.stability_review.resolve() if args.stability_review else None)
    except (OSError, ValueError, KeyError) as exc:
        parser.exit(1, f"Speed comparison rejected: {exc}\n")
    print(json.dumps({key: report[key] for key in ("status", "new_jobs", "baseline_jobs",
        "all_200m_s_jobs_have_end_of_job", "all_200m_s_equilibrium_and_root_screens_pass",
        "all_200m_s_within_small_deflection_and_elastic_screens")}, indent=2))


if __name__ == "__main__":
    main()
