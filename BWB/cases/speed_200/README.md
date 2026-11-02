# 200 m/s 高动压工况

保持既有结构、材料、厚度、根部约束、CAERO1、SPLINE1 和 W2GJ 不变，将速度由 20 m/s 改为 200 m/s；两个零件各自计算迎角 3°、4°、5°、6°。密度仍假设为 1.225 kg/m³、声速为 340.3 m/s，因此 TRIM 使用 Mach 0.5877167205、动压 24,500 Pa。侧滑为零，AEQR=1，均实际调用 SOL144 求解。

```text
speed_200/
├── aoa_03/part_1/、part_2/
├── aoa_04/part_1/、part_2/
├── aoa_05/part_1/、part_2/
└── aoa_06/part_1/、part_2/
```

每件目录含完整 BDF 入口、飞况、四组结构/气动输入和来源清单。各自独立固定根部，没有整机配平、惯性载荷、rib 或 spar。材料仍为 E=69 GPa、ν=0.33、1 mm 铝蒙皮；200 MPa 只是后处理假设屈服阈值，没有塑性材料本构。

本次八个算例均完成结果恢复与载荷平衡检查，但全部超过假设屈服强度及原有小变形筛查界限。因此保留线性输出用于识别当前建模假设的限制，不能将应力、应变、位移解释为屈服后的真实响应。未计算静气动弹性发散临界动压，也未模拟屈曲、大变形、塑性或流动分离。详见 [实际结果](../../results/speed_200/README.md)。

在仓库根目录运行一次前台批处理（需要本机 Nastran 及许可证）：

```powershell
python -m BWB.scripts.run_speed_200 --solve-only
python -m BWB.scripts.run_speed_200 --postprocess-only
```

脚本拒绝覆盖已有求解 F06；本仓库保留本次结果，重跑应指定新的 `--run-root workspace/BWB/speed_200/新目录`。默认后处理会更新该扫描的正式导出，需先保留需要追溯的旧成果。旧 20 m/s 的输入和结果不被修改。

本次额外对 200 m/s、6°、零件 1 运行 `MKL_VERBOSE=1` 的独立重复作业，保存在 `workspace/BWB/speed_200/aoa_06_part_1_trace/`。详细跟踪将同一运行诊断定位为空输出矩阵调用，所有导出的应变和结构/气动位移 CSV 与普通运行逐字节相同；原始错误与非零退出码仍保留。运行 `python -m BWB.scripts.review_speed_200_runtime` 可核对该现有跟踪，不启动新作业。大型跟踪日志保留在 workspace，Git 中记录其散列和关键原文。

汇总及 Femap 导入：

```powershell
python BWB/scripts/summarize_speed_200.py --stability-review BWB/results/speed_200/stability_review.json
python BWB/scripts/view_femap_speed_sweep.py --run-root workspace/BWB/speed_200/run_01 --speed-mps 200 --default-result bottom-vonmises-stress --allow-reviewed-empty-matrix-diagnostic
```

Femap 导入器拒绝覆盖同名结果模型，重建时用 `--output-dir` 指定 `BWB/results/` 下新的目录。审查依据通过独立的 `--review-dir` 指定，默认读取本次 `BWB/results/speed_200`；需与实际运行的速度和 F06 对应。每件四组结果按 AOA 3°–6°保存，默认 6° bottom von Mises，应力原生单位为 Pa（除以 10⁶ 为 MPa）。模型仅用于查看实际 OP2；求解仍以外部完整 BDF/W2GJ 为准。

高动压放大了 F06 打印舍入的绝对差值。后处理按每个实际打印数值的末位精度传播 `pressure=q×Cp` 和节点反力合成误差，记录残差与容差，并用测试确认真实不一致仍会被拒绝；载荷平衡、根部、小变形与屈服筛查阈值保持原值。
