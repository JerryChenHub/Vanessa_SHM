# NASA_SHM 与 BWB

本仓库承载两个独立项目，沿用一个 Git 仓库和一个 `workspace/`。外层目录名 `Vanessa_SHM` 暂时保留，项目归属以子目录为准。

```text
Vanessa_SHM/
├── NASA_SHM/                 原项目的代码、模型和正式数据
│   ├── core/                 代码、基础模型、网格参考数据和图片
│   ├── tests/                原项目回归测试
│   ├── data_2026-10-17/       已有正式数据集及计算记录
│   └── README.md             NASA_SHM 使用说明
├── BWB/                      BWB 蒙皮、气动面板与后续分析
│   ├── paths.py              BWB 自身路径
│   ├── data/                 BWB 自有输入数据
│   ├── scripts/              CAD 导出、蒙皮/气动生成、检查和绘图
│   ├── tests/                气动卡片与几何约定的回归测试
│   ├── models/skin_demo/     CAD 输入与首版自由网格参考
│   ├── models/skin_structured/ 当前规则蒙皮网格、检查及对比图
│   ├── models/aero_panels/   两件独立的气动/结构输入与 Femap 查看模型
│   ├── cases/sol144_demo/    假设飞况下可重建的 SOL144 输入
│   ├── results/sol144_demo/  实际应变、云图、Femap 结果与求解证据
│   └── README.md             项目边界与复用约定
└── workspace/                共享临时工作区，不提交计算文件
    ├── NASA_SHM/             原有 case、结果、诊断和历史归档
    └── BWB/                  BWB 临时文件、计算和中间结果
```

## 项目边界

- 各项目的代码、固定输入、测试和正式成果归各自目录；临时输出放入各自的 `workspace/<项目名>/`。
- BWB 可以显式导入 `NASA_SHM.core` 或读取 NASA_SHM 数据。调用时传入 BWB 自己的输入和输出路径，避免使用 NASA_SHM 专属的模型与输出默认值。
- NASA_SHM 不依赖 BWB。BWB 专用改动放入 BWB；复用现有数据时保留来源、版本或文件哈希、单位和变换说明，派生成果写入 BWB。
- NASA_SHM 原始计算数据、求解日志和源码快照保留；共享 workspace 不改变数据所有权。

## 从仓库根目录使用

```powershell
# NASA_SHM 生成与检查入口
python -m NASA_SHM.core.data_generation --help
python -m NASA_SHM.core.batch_run --help

# 运行已有回归测试
Push-Location NASA_SHM
python -m unittest discover -s tests -v
Pop-Location
```

NASA_SHM 默认生成位置为 `workspace/NASA_SHM/cases/`。显式传入的相对路径仍以当前工作目录为基准。

## 两个项目的数据流

NASA_SHM 的流程为：基础 BDF 与网格索引 → 刚度损伤样本 → SOL144 → 位移/应变结果与标签。已有正式数据集和求解记录留在 NASA_SHM 或所属 workspace 中，不复制成 BWB 的计算结果。

BWB 的流程为：原始 CAD 的只读副本 → SolidWorks API 导出 Parasolid → Femap API 生成规则蒙皮 → 从蒙皮行列提取参考升力面 → 生成 CAERO1、SPLINE1 和 W2GJ → 独立检查与 Femap 查看。气动文件沿用 NASA_SHM 的卡片分文件组织，但几何、插值节点和初始弯度/扭转输入由 BWB 自己计算；没有移用 NASA 的气动数值、工况或应变结果。

当前 BWB 气动包定义每件 11 个 CAERO1 大面、590 个小面板及对应结构连接。独立的 `cases/sol144_demo/` 已补齐假设飞况：20 m/s、迎角 2°、侧滑 0°、1 mm 铝蒙皮、根部固定。两件均实际运行 SOL144，得到完整的位移与蒙皮两表面应变，详见 [结果与合理性判断](BWB/results/sol144_demo/README.md)。两件保持各自 CAD 坐标并独立求解，不代表整机装配或自由飞行配平。

两件的输入检查、结果覆盖、力/力矩平衡和量级筛查均通过。求解器正常到达 END OF JOB，仍保留一条已追踪到零行矩阵的 MKL 诊断和非零启动器退出码；没有将它们改写为无错误运行。当前成果适合流程演示，尚未完成气动或网格精度验证，也未生成 BWB 损伤数据集。

详细说明：[NASA_SHM](NASA_SHM/README.md)、[BWB](BWB/README.md)、[气动生成与文件说明](BWB/models/aero_panels/README.md)。Windows、Femap 和 Python 依赖配置见 BWB 说明。

迁移前文件校验清单、旧索引备份及迁移校验记录位于 `workspace/NASA_SHM/archive/project-split/`。
