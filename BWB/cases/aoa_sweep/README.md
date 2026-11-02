# BWB 迎角 3°–6° 扫描输入

延续 [2° 基准工况](../sol144_demo/README.md)，速度 20 m/s、密度 1.225 kg/m³、Mach 0.05877167205、动压 245 Pa、侧滑 0°、1 mm 铝蒙皮及各自根部固定均保持不变。仅修改迎角，并将度转换为 TRIM 使用的弧度。

```text
aoa_sweep/
├── aoa_03/part_1/、part_2/   3°，两件独立求解
├── aoa_04/part_1/、part_2/   4°
├── aoa_05/part_1/、part_2/   5°
└── aoa_06/part_1/、part_2/   6°
```

每件目录都是完整独立输入，含 `bwb_sol144.bdf`、结构/气动/样条/W2GJ、`flight_condition.bdf` 和来源清单。8 个工况均实际求解，不由 2° 结果按比例生成。保留初始弯度 W2GJ，因此响应不要求经过迎角为零时的原点。

本次命令（仓库根目录执行，需要已安装的 Nastran 和可用许可证）：

```powershell
python -m BWB.scripts.run_aoa_sweep --run-root workspace/BWB/aoa_sweep/run_01
python -m BWB.scripts.summarize_aoa_sweep
```

第一条是一次性前台批处理，没有创建定时任务。它生成输入、顺序求解、提取应变、检查力/力矩和量级，并逐件归档求解证据。为保留旧结果，脚本拒绝覆盖既有运行 F06 或已完成的扫描成果；重跑单个算例可用 `BWB.scripts.run_sol144` 指定新的临时目录，再调用现有导出与检查入口。

实际结果、对比图、所有应变 CSV 与 Femap 多工况模型位于 [扫描结果目录](../../results/aoa_sweep/README.md)。沿用基准工况的数学库诊断记录，所有原始错误保留；接受为演示仍需每个工况的完整性、平衡、弹性反馈和合理性检查通过。线性势流模型不能验证这些迎角下真实流动是否分离或失速。

在保留上述运行目录、安装 Femap 的本机上，可另存一个包含四组 OP2 结果的模型：

```powershell
python BWB/scripts/view_femap_aoa_sweep.py --run-root workspace/BWB/aoa_sweep/run_01 --output-dir BWB/results/aoa_sweep_rebuild --allow-reviewed-empty-matrix-diagnostic
```

输出目录应尚未包含同名模型。该选项仅允许此次已审阅的空矩阵诊断组合，原始错误仍写入各结果集备注；不能用来跳过其他求解失败。
