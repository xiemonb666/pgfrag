# 独立源码与结果审阅

结论：指定范围内未发现可复现的误报物理 PASS。薄壁工具把局部几何前置与物理/全局资格分开；bias 报告把样例与物理结论分开。现有结果仍是解析夹具，不能据此声称真实 mesh bias 已消除。

薄壁检查使用实际输入的 annular-wedge 单元区间、射线连续覆盖和加权质心投影层数；仅 `h/t` 不会给 PASS。解析结果中，t=2 mm、h=0.5 mm 为 4 区间/4 层、`LOCAL_RESOLUTION_PRECONDITION_MET`；h=2 mm 为 1/1、`UNDERRESOLVED`。单位按每个量显式换算；pointweight 与解析环扇体积逐单元对账，质量按密度×体积核对。段检查在圆柱内外半径与 z 平面的断点间检查整段，能检出端点在材料内、连线穿空腔的情况。代码固定输出 `physical_pass=false`、`global_pass=false`。`actual_section_samples.png` 展示解析 fixture，不是生产 mesh 截面。

bias 配对固定几何、load、材料参数/场 hash、随机种子、求解器与公共 rest-frame ID；路径/表面先按 `local_to_common` 变换。路径量是离散样点的对称最近邻 p95/Hausdorff，面积为三角面重算，能量逐项按单位换算，起裂时间比较绝对差；这些定义在 README 中明示。无断裂返回 UNKNOWN，图连通分量不冒充几何碎片，缺指标/资格也不能得到 within-tolerance 状态；样例报告为 `FIXTURE_ONLY`、`NOT_ADMITTED`、`bias_eliminated=false`。

限制：最近样点距离不等同连续曲线距离，结果依赖离散采样；native qualification 是描述符自报 PASS，需后续来源审查。没有调用测试、实现、solver 或 CAD；对既有 JSON 只做一次 PowerShell 读取汇总（<1 s）。
