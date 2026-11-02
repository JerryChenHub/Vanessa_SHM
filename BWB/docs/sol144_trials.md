# SOL144 临时试算

飞况探索与流程试跑产生的文件统一归入 `workspace/BWB/`，不作为仓库成果。可复用脚本、代码回归测试及当前 CAD/结构/气动模型仍放在 `BWB/`。

```text
workspace/BWB/
├── trials/
│   ├── cases/
│   │   ├── sol144_demo/       假设基准输入
│   │   ├── aoa_sweep/         20 m/s迎角扫描输入
│   │   └── speed_200/         200 m/s试算输入
│   ├── results/              同名分组的CSV、图、Femap及诊断/证据
│   └── relocation_manifest.json  本机迁移前后路径及逐文件SHA256
├── sol144_demo/               已有求解器工作目录
├── aoa_sweep/                 已有扫描工作目录
└── speed_200/                 已有高速试算工作目录
```

这些目录由运行过程产生，干净克隆不包含历史试算数据。旧诊断和来源清单保留生成时的路径及散列；本机迁移清单记录新位置。历史文件中的旧路径是来源记录，不能据此重新向正式目录写入试算。求解工作目录未迁移。

## 生成、求解和提取

在仓库根目录执行，需要本机安装 Nastran、有效许可证及 Python 依赖：

```powershell
# 生成临时输入；输出目录必须位于workspace/BWB
python -m BWB.scripts.build_sol144_demo --output workspace/BWB/trials/cases/example --velocity 20 --alpha-deg 3

# 每次使用新的求解目录，保留已有F06
python -m BWB.scripts.run_sol144 workspace/BWB/trials/cases/example/part_1 workspace/BWB/example/run_01

# 导出和检查的参数见帮助；输出指定到workspace/BWB/trials/results
python -m BWB.scripts.export_sol144_strain --help
python -m BWB.scripts.check_sol144_reasonableness --help
```

`run_aoa_sweep.py`、`run_speed_200.py` 保留作为批处理代码，其输入与结果默认写入上述临时目录，原始求解工作目录也在 `workspace/BWB/`。已存在的求解 F06 不被覆盖。脚本对诊断和历史基准的依赖需在本机准备；没有相应审查依据时，应保留原始错误并完成审查，不能把退出码改写为成功。

汇总脚本的默认输入和输出同样位于 `trials`。Femap 导入器默认从 `trials/results` 读取审查并写入独立结果模型；重建用新的 `--output-dir`，必要时另设 `--review-dir` 指向对应的已审查结果。结构和气动源模型仍从 `BWB/models/` 只读获取。

## 保留的实现要点

- 工况以 SI 单位组织，Mach、动压和迎角弧度由指定条件计算；不继承 NASA 模型的单位换算或飞况数值。
- TRIM 显式提供 `ANGLEA`、`SIDES=0` 和 `AEQR=1`，并检查非零气动刚度与变形反馈。W2GJ 以完整 DMI 列保留，不能由几何查看模型的重导出代替。
- 完整性检查覆盖结构节点、壳单元两纤维面、根部约束、压力/系数一致性及反力平衡。打印舍入容差来自 F06 的有效位数，物理阈值独立保留。
- 应力限值是工况中明确记录的材料假设；数值完成、平衡通过或误差可复现不能替代小变形、屈服、屈曲、发散及气动适用性判断。
- 原始错误和来源哈希不删改。运行诊断的审查范围只适用于其证据覆盖的模型、求解器和工况。

代码回归测试会在 `workspace/BWB/` 中生成隔离输入，结束后清理，不要求历史试算存在：

```powershell
python -m unittest discover -s BWB/tests -q
```
