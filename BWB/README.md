# BWB

本项目从两个 SolidWorks 零件生成规则蒙皮，再由同一套蒙皮行列建立气动参考面、结构插值连接和弯度/扭转修正。已在假设的 20 m/s、2° 迎角工况下分别运行两件的 SOL144，获得完整应变；两件尚未装配为整机。详见 [求解输入与复现](cases/sol144_demo/README.md) 和 [实际结果与合理性判断](results/sol144_demo/README.md)。

- [规则蒙皮模型](models/skin_structured/README.md)：主翼采用平行弦向截面和连续展向列，保留材料、1 mm 厚度及根部固定的演示假设。
- [气动面板模型](models/aero_panels/README.md)：每件 11 个 CAERO1、590 个气动小面板，附 SPLINE1/SET1、W2GJ、独立检查报告和 Femap 查看模型。
- `models/skin_demo/`：原始 CAD 副本、Parasolid 导出及首版自由网格，保留作来源记录和对比。
- `cases/sol144_demo/`：假设工况、完整 SOL144 入口、输入来源与求解诊断说明。
- `results/sol144_demo/`：应变 CSV、位移、合理性检查、云图、Femap 结果模型及压缩求解证据。

## 目录与生成逻辑

```text
BWB/
├── paths.py                          BWB 路径定义
├── data/                             BWB 自有输入数据
├── scripts/
│   ├── ExportSkinGeometry.cs         SolidWorks API 导出实现
│   ├── export_skin_geometry.ps1      CAD 导出入口
│   ├── build_femap_skin_demo.py      规则蒙皮生成入口
│   ├── structured_skin_grid.py       CAD 投影与规则行列拓扑
│   ├── build_bwb_aero_panels.py       气动面、SPLINE 与 W2GJ 生成
│   ├── check_aero_panels.py          从导出卡片独立重建并检查
│   ├── view_femap_aero_panels.py     导入原生 Aero Panel 与查看模型
│   └── plot_aero_panels.py           结构/气动叠加图
├── tests/test_aero_panels.py          几何、符号、编号与输入检查
└── models/
    ├── skin_demo/                    源 CAD 与历史自由网格
    ├── skin_structured/              当前规则结构网格
    └── aero_panels/                  当前气动输入包及可追溯报告
        ├── part_1/
        ├── part_2/
        ├── aero_panel_overview.png
        └── aero_spline_support.png
```

气动生成以 `skin_structured/part_N/part_N_skin.bdf` 和 `mesh_layout.json` 为输入。上下蒙皮配对提取平均线，主翼、短过渡、弯折、翼梢分别建立顺流的平面 CAERO1；弯折分成 8 个小平面。SPLINE1 将每个大面连接到附近单侧蒙皮上分布的结构节点；当前支撑选取是建模方案，尚未验证其载荷传递和局部蒙皮应变精度。真实平均线相对参考面的弯度与扭转以 W2GJ 初始下洗补入，气动网格本身不包裹上下蒙皮。

`aero_panels/part_N/part_N_aero_structure.bdf` 汇总结构与气动 bulk 卡片用于导入；它没有 SOL、TRIM 或求解输出请求。实际求解入口是 `cases/sol144_demo/part_N/bwb_sol144.bdf`，保留全部四组输入并增加工况。原生气动查看 `.modfem` 没有保存 W2GJ；求解必须包含外部 `dmicards.bdf`，不能盲目重导出查看模型。`results/` 中的 `.modfem` 已导入实际 OP2 应变供查看。气动精度、结构收敛和 BWB 损伤数据生成仍未完成。

- `BWB/`：本项目代码、说明和需要长期保留的成果。
- `BWB/data/`：本项目输入数据，记录来源、单位和版本。
- `workspace/BWB/`：本项目临时计算文件、日志和中间结果，与 NASA_SHM 的临时文件分开存放。
- `BWB/paths.py`：本项目路径定义；导入时不创建目录，也不加载 NASA_SHM。

蒙皮生成入口为 `BWB/scripts/build_femap_skin_demo.py`；在仓库根目录执行，需 Windows 桌面 COM 和 Femap。气动卡片生成和独立检查直接处理 BDF/JSON，不需要启动 Femap；原生查看模型由 `view_femap_aero_panels.py` 通过 Femap COM 生成。这些建模入口均不调用求解器。具体重建命令见各模型目录的 README。

当前两件气动输入的独立检查均为 PASS，原生面板和样条与输入一致。SOL144 已完成完整结果恢复，力/力矩平衡和粗略量级检查通过；保留空矩阵 MKL 诊断、9 个封口单元的几何检查超限，并未把流程演示称为应变精度验证。详细数字和处理边界见结果报告。

## 本机环境

安装与当前 Python 版本兼容的 `pywin32`，以及网格处理和绘图所需的 `numpy`、`scipy`、`matplotlib`：

```powershell
python -m pip install pywin32 numpy scipy matplotlib
```

Femap 引导脚本优先使用当前 Python 环境；需要单独的 pywin32 运行库时，可设置 `BWB_PYWIN32_PATH` 指向其目录。Femap 的 `femap.tlb` 默认从 COM 注册信息或标准安装目录查找；非标准安装可用环境变量 `FEMAP_TLB` 指定完整文件路径。首次 CAD 导出还需要 SolidWorks；导出入口见 `scripts/export_skin_geometry.ps1`。

只导入 `scripts/femap_connection.py` 会读取类型库并生成本项目的 COM 包装缓存，不会连接或启动 Femap；运行建模入口时需要能够访问 Windows 桌面 COM 的会话。缓存写入 `workspace/BWB/femap_demo`，不提交 Git。

## 跨项目复用

如需复用 NASA_SHM，显式导入 `NASA_SHM.core` 中需要的模块，或读取 `NASA_SHM/data_2026-10-17/` 中的数据。将这些来源作为只读依赖；不修改 NASA_SHM 的代码、原始输入和已有结果。调用复用函数时显式传入输入和输出路径，不沿用 NASA_SHM 的默认输出位置。BWB 生成的文件应写入本项目目录或 `workspace/BWB/`。

记录所用文件的源路径、SHA256、单位和工况；复用代码时记录代码版本或文件 SHA256，并核对其模型假设是否适用于 BWB。需要改变已有方法时，在 BWB 中编写适配代码。

本次气动建模参考 `NASA_SHM/core/base/aero_cards.bdf`、`spline_cards.bdf`、`dmicards.bdf` 的组织方式，并在清单中记录来源。NASA 气动输入含 500 个小面板；其文件注释提到的 MATLAB 生成脚本不在当前仓库，不能声称复刻了未知的几何插值算法。BWB 不复制 NASA 的坐标、W2GJ 数值、飞行工况、对称条件或历史结果，也不修改 NASA_SHM。

仅在当前 Python 未安装 pywin32、且没有设置运行库路径时，COM 引导脚本可只读使用本机已有的 `workspace/NASA_SHM/femap_preview/python_deps` 作为兼容回退。该临时目录不随 Git 仓库分发，新环境应自行安装依赖。
