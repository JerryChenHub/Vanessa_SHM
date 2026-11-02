# BWB SOL144 演示工况

两个零件独立固定根部，使用同一组假设；不作左右翼装配或自由飞行配平。完整输入均为 m–N–kg–s。

| 项目 | 假设 |
|---|---:|
| 速度 / 密度 / 声速 | 20 m/s / 1.225 kg/m³ / 340.3 m/s |
| Mach / 动压 | 0.05877167205 / 245 Pa |
| 迎角 / 侧滑角 | 2° / 0° |
| 蒙皮 | 各向同性铝，E=69 GPa，ν=0.33，ρ=2700 kg/m³，厚 1 mm |
| 约束 | 每件 49 个根部节点，DOF 123456 固定 |
| 载荷 | 仅气动力，无重力、惯性载荷、附加质量、残余应力 |
| 结构 | 2940 GRID，2914 CQUAD4 + 1 CTRIA3，无 rib/spar |
| 气动 | 11 CAERO1，590 小面板；线性亚声速势流，保留 W2GJ |

不输入 CL。Nastran 根据气动几何、迎角、初始弯度/扭转计算压力，通过 SPLINE1 向结构传递载荷，并把结构位移反馈给气动模型。TRIM 的 AEQR=1；显式 SIDES=0 规避本机版本在第二组变量空白时未启用弹性反馈的实测问题。

每件目录内：

```text
bwb_sol144.bdf         SOL144、SPC、TRIM及输出请求
flight_condition.bdf  AESTAT ANGLEA/SIDES、Mach、q、角度、AEQR
structural_cards.bdf  现有蒙皮、材料、厚度、根部约束
aero_cards.bdf        气动坐标、参考尺度、PAERO1、CAERO1
spline_cards.bdf      SPLINE1与SET1
dmicards.bdf          单一连续列的590个W2GJ初始下洗
case_manifest.json   假设、源文件与SHA256
```

从仓库根目录重建和运行（已安装 Nastran，许可证服务器可达）：

```powershell
python -m BWB.scripts.build_sol144_demo
python -m BWB.scripts.run_sol144 BWB/cases/sol144_demo/part_1 workspace/BWB/sol144_demo/my_part_1
python -m BWB.scripts.export_sol144_strain --case-dir workspace/BWB/sol144_demo/my_part_1 --output-dir BWB/results/sol144_demo/my_part_1
python -m BWB.scripts.check_sol144_reasonableness --case-dir workspace/BWB/sol144_demo/my_part_1 --export-dir BWB/results/sol144_demo/my_part_1 --part-index 1 --output BWB/results/sol144_demo/my_part_1/reasonableness.json
```

part_2 同理。使用新的运行目录，脚本会拒绝覆盖已有 F06。`build_sol144_demo --help` 可修改速度、密度、声速和迎角；这会重建输入，不会自动执行求解。检查必须针对新的输入/结果重做，不能沿用本次合理性结论。

运行器对任何运行时错误保持非零返回码；本次已知的空矩阵诊断也不例外。随后应检查 `run-summary.json`、完整结果与 [求解诊断记录](solver_notes.md)，不得只用进程退出码宣告成功或失败。最终成果见 [结果说明](../../results/sol144_demo/README.md)。

NASA_SHM 仅作为 SOL144 卡片组织参考，并只读复用日志检查和 F06 应变读取函数；没有修改 NASA_SHM 或套用它的英制参数、几何、载荷和应变。
