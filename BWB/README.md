# BWB

已建立基于两个 SolidWorks 零件的 Femap 结构化蒙皮面网格，当前说明与原生模型见 `models/skin_structured/README.md`。主翼采用平行弦向截面和连续展向列；首版自由网格与原始 CAD 输入保留在 `models/skin_demo/` 供对比。尚未建立气动面板或运行 SOL144。

- `BWB/`：本项目代码、说明和需要长期保留的成果。
- `BWB/data/`：本项目输入数据，记录来源、单位和版本。
- `workspace/BWB/`：本项目临时计算文件、日志和中间结果，与 NASA_SHM 的临时文件分开存放。
- `BWB/paths.py`：本项目路径定义；导入时不创建目录，也不加载 NASA_SHM。

网格生成入口为 `BWB/scripts/build_femap_skin_demo.py`；在仓库根目录执行，需 Windows 桌面 COM 和 Femap。该脚本不调用求解器。

## 本机环境

安装与当前 Python 版本兼容的 `pywin32`，以及网格处理和绘图所需的 `numpy`、`scipy`、`matplotlib`：

```powershell
python -m pip install pywin32 numpy scipy matplotlib
```

Femap 引导脚本优先使用当前 Python 环境；需要单独的 pywin32 运行库时，可设置 `BWB_PYWIN32_PATH` 指向其目录。Femap 的 `femap.tlb` 默认从 COM 注册信息或标准安装目录查找；非标准安装可用环境变量 `FEMAP_TLB` 指定完整文件路径。首次 CAD 导出还需要 SolidWorks；导出入口见 `scripts/export_skin_geometry.ps1`。

只导入 `scripts/femap_connection.py` 会读取类型库并生成本项目的 COM 包装缓存，不会连接或启动 Femap；运行建模入口时需要能够访问 Windows 桌面 COM 的会话。缓存写入 `workspace/BWB/femap_demo`，不提交 Git。

## 跨项目复用

如需复用 NASA_SHM，显式导入 `NASA_SHM.core` 中需要的模块，或读取 `NASA_SHM/data_2026-10-17/` 中的数据。将这些来源作为只读依赖；不修改 NASA_SHM 的代码、原始输入和已有结果。调用复用函数时显式传入输入和输出路径，不沿用 NASA_SHM 的默认输出位置。BWB 生成的文件应写入本项目目录或 `workspace/BWB/`。

记录所用文件的源路径、SHA256、单位和工况；复用代码时记录代码版本或文件 SHA256，并核对其模型假设是否适用于 BWB。需要改变已有方法时，在 BWB 中编写适配代码。当前没有加载 NASA_SHM 分析代码或数据；仅在当前 Python 未安装 pywin32、且没有设置运行库路径时，COM 引导脚本可只读使用本机已有的 `workspace/NASA_SHM/femap_preview/python_deps` 作为兼容回退。该临时目录不随 Git 仓库分发，新环境应自行安装依赖。
