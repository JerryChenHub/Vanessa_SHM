# NASA_SHM

以基础 Nastran 模型生成刚度损伤 case，在本机运行包含气动力的 SOL144，并提取位移及壳膜内应变。

## 目录

```text
NASA_SHM/
├── core/                       核心 Python 包及固定输入
│   ├── paths.py                统一默认路径，以项目位置为基准
│   ├── base/                   原 base_file；基础 BDF / INCLUDE 文件
│   ├── data/                   原 _data；网格索引及历史参考 CSV
│   ├── assets/                 原 img；参考图片
│   ├── bdf_damage.py           BDF 读取与损伤材料/属性修改
│   ├── nodal_info.py           网格读取、选区及查看工具
│   ├── data_generation.py      从基础模型生成损伤 case
│   ├── run_local_nastran.py    本机求解器调用与完成状态检查
│   ├── export_results.py       结果校验、CSV 和损伤标签导出
│   └── report_figure.py        绘图工具
├── tests/                      回归测试
└── data_2026-10-17/             正式数据集及计算记录；不提交 Git

../workspace/NASA_SHM/           仓库共享工作区内的本项目区域；不提交 Git
    ├── cases/                  新生成 case 的默认位置
    ├── figures/                新生成图片的默认位置
    ├── local-validated-20261017/ 已完成的新损伤样本及结果索引
    ├── local-sol144-20261017-155146/ 基准结果和早期试算
    └── archive/                旧一次性脚本、缓存及迁移前元数据备份
```

`core/data` 中的 `1.csv`、`9.csv`、`undamged.csv` 是历史参考文件，保留原名；新结果写入各 case 的 `results/`。
原 `Temp/` 后来整理为 `workspace/`，现整体位于仓库根的 `workspace/NASA_SHM/`。原始求解日志保留当时的路径，活动结果索引中的路径已更新。

## 使用

使用 Python 3.11 或更新版本。生成及绘图使用 `numpy`、`pandas`、`matplotlib`；图片画圈工具额外使用 `opencv-python`。求解和结果提取模块使用 Python 标准库。本机需要已安装的 Simcenter Nastran；使用 UCI 许可时先连接 VPN。

以下命令在 `NASA_SHM/` 内运行（先从仓库根执行 `cd NASA_SHM`）；Python 内调用可使用 `from core.<模块> import ...`。从仓库根或 BWB 复用时使用 `from NASA_SHM.core.<模块> import ...`，模块命令相应写为 `python -m NASA_SHM.core.<模块>`。

```powershell
# 1. 生成一个损伤状态，包含模板中的全部 10 个工况
python -m core.data_generation --samples 1 --seed 20261017

# 2. 调用本机求解器
python -m core.run_local_nastran ../workspace/NASA_SHM/cases/sample_1

# 3. 校验结果并提取 CSV
python -m core.export_results ../workspace/NASA_SHM/cases/sample_1

# 4. 绘图（默认输出 ../workspace/NASA_SHM/figures/strain_presence.png）
python -m core.report_figure --strain-path ../workspace/NASA_SHM/cases/sample_1/results/strain_wide.csv

# 仅检查已有求解结果
python -m core.run_local_nastran ../workspace/NASA_SHM/cases/sample_1 --inspect-only

# 回归测试
python -m unittest discover -s tests -v
```

使用 `--output-dir ../workspace/NASA_SHM/cases/another_batch` 创建另一批 case。生成器拒绝覆盖已有 `sol144.f06` 的目录。
默认输入和输出位置由 `core/paths.py` 决定，不依赖当前工作目录；显式传入的相对路径仍相对于调用者的工作目录。

## 数据流和约定

`core/base` + `core/data` → 复制及损伤修改 → `../workspace/NASA_SHM/cases/sample_i` → Nastran SOL144 → F06/XDB → `results/*.csv`。

- 一个 sample 是一种损伤状态；当前模板有 10 个飞行工况。`--samples 10` 表示 10 种损伤状态、100 个工况结果。
- `damage_ratio` / `stiffness_scale` 表示**剩余刚度比例**。壳膜内、弯曲、剪切及耦合材料引用和梁的 E/G 同比降低；密度、几何保持不变。
- `case_manifest.json` 保存随机种子、选区、实际修改单元及输入哈希。
- `run-summary.json` 记录求解完成状态、警告及运行诊断；不能仅根据启动器退出码判断计算是否完成。
- `results/strain_wide.csv` 为单元 ID 加每工况 3 个膜内应变分量；坐标系是单元坐标系，剪切量为工程剪应变。
- CSV 来自 F06 打印精度；壳膜内应变与曲率分别处理。梁没有出现在壳应变 CSV 中，不代表梁应变为零。
- 目录整理不改变损伤、载荷或结果提取的物理定义。旧诊断脚本在 `../workspace/NASA_SHM/archive/legacy_scripts` 仅作记录，日常使用 `core` 中的入口。

已完成结果入口：[RESULTS.json](../workspace/NASA_SHM/local-validated-20261017/RESULTS.json)。

## 后台批量数据

```powershell
python -m core.batch_run prepare --output-dir data_2026-10-17 --count 100 --seed 20261017
python -m core.batch_run run --output-dir data_2026-10-17
```

第二条命令可以作为独立后台进程运行；恢复中断的批次仍使用同一条命令。`NASA_SHM/data_<日期>/` 保留三层应变 CSV、每个 case 的位置和损伤参数、整批抽样规则、输入压缩包、校验信息及运行进度。`status.json` 显示完成数和预计结束时间。批处理临时目录位于仓库根的 `workspace/NASA_SHM/batch_data_<日期>/`。

`core/plate_strain.py` 读取 `STRAIN(FIBER)`：top 为壳局部 +z，bottom 为局部 -z，mid 为两个表面应变分量的算术平均。每个表为 12,475 个壳单元、10 工况 × 3 分量。批处理在成功导出并校验后清理对应 workspace 临时目录；失败 case 留存诊断文件。工作期间阻止 Windows 空闲休眠，退出时恢复；计算仍需要本机运行及 VPN/许可连接。

进度文件遇到 Windows 文件占用时自动重试；持续占用只延迟进度视图更新，运行状态仍保存到 `batch.json`。每次求解使用独立的 `solve_*` 临时目录，避免重试复用旧 XDB；异常退出会先终止并回收本次求解进程。恢复运行会校验并跳过已完成的数据，成功补算后将旧失败记录移入 `failure-history.jsonl`。

启动求解前检查 Windows 可提交内存（求解器额度加 1 GiB 余量）及可用物理内存；不足时显示 `waiting_resources`，每 30 秒检查并自动恢复。内存分配失败会保留诊断并等待重试，不消耗普通 case 失败次数。ETA 优先采用本次运行的速度；等待资源期间不显示确定完成时间。

## 项目拆分记录

本项目已从仓库根迁入 `NASA_SHM/`，临时文件位于 `../workspace/NASA_SHM/`。BWB 可显式复用本项目代码与数据，专用改动和派生结果归 BWB。

现有批次的路径和当前可执行源码校验值已完成迁移。`batch.json` 的 `layout_migration` 记录迁移前源码哈希及新旧位置；原始 `shared/pipeline_source.zip` 和各 case 的计算元数据保持不变，新布局代码另存为 `shared/pipeline_source_after_project_split.zip`。迁移前批次索引另存为 `shared/batch_before_project_split.json`。其他原索引备份和文件完整性报告见 `../workspace/NASA_SHM/archive/project-split/`。历史日志、运行命令及旧诊断脚本仍记录当时的路径。
