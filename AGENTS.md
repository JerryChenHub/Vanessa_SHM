# 工作约定

- `NASA_SHM/` 与 `BWB/` 是两个独立项目。代码、固定输入、测试和正式成果放在所属项目内。
- 共用根目录 `workspace/`，临时内容分别放入 `workspace/NASA_SHM/` 与 `workspace/BWB/`。
- BWB 复用 NASA_SHM 时使用显式导入与路径，记录数据来源和处理过程；BWB 专用改动不得混入 NASA_SHM，NASA_SHM 不反向依赖 BWB。
- SolidWorks 建模调用 API 完成；Fillet 倒角等在最后添加；螺丝等优先使用 TOOLBOX。建模中关闭暂时不需要的窗口。
- 禁止开启定时任务。
