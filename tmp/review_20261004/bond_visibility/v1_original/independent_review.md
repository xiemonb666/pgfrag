# Bond visibility v1 独立只读审阅

审阅对象固定为下列 SHA-256；本报告只对这些字节成立。未运行 Python、CAD、pytest、STEP、Julia、solver 或 Git；只读源码/设计和已安装 CadQuery/OCP 文件。

| 文件 | SHA-256 |
|---|---|
| `pgfrag/src/pgfrag/bond_visibility.py` | `7b87ffd1bb94f447de4ed379b9e047aa0da3cf02626629bb86b6ba8fe97eb86` |
| `pgfrag/tests/test_bond_visibility.py` | `533702781d8d27295260b108e2783949cdac1691709704fa96edfed6fccfa84a` |
| `pgfrag/tmp/review_20261004/bond_visibility/design_and_freeze.md` | `1e943da11387344169d66350dbdb9740db18ed6b245f8e5643c318141884fe2` |
| `pgfrag/tmp/review_20261004/cavity_bond_safety/design.md` | `6035dfa7c3fdeade43a61ff653b03f55f002326ddebe95b5228db0d0b245e0b` |

## 结论

核心方法选对了：在有限线段上取 B-rep 与面交点，再逐开区间用 OCC solid classifier 判断占据；端点、区间、相反域和不确定几何按三态处理。实现还不能作为入 native 前的准入证明。v1 有一个可能把短空腔吞掉并返回 `IN_MATERIAL` 的聚类缺陷；precision/provenance 的可信度依赖调用者自述；UNKNOWN 的 asset 级拒绝门禁及 native 接入均不存在。当前 `whole_asset_physics_status` 正确写为 `NOT_ASSESSED`，但这也意味着没有 asset admission 结果。

## P1 — 聚类链可能吞掉空腔并误放行

`bond_visibility.py:446-497` 使用相邻 hit 距离做 single-linkage 聚类（`parameter - previous <= tolerance`），但仅当整个组的总跨度 `<= tolerance` 才写 `UNRESOLVED_NEARBY_BOUNDARY_HITS`。如果同一 solid、同一 face 上一串多次交点的每个相邻间距都不超过 tolerance、而组总跨度超过 tolerance，`solid_faces` 只有一个元素，`hit_reasons` 为空；第 481 行把整串压成平均事件。`_pair_intervals` 随后只在这个均值两侧分类（568-626），原本应因宽度落入误差带而为 UNKNOWN 的小 VOID 开区间已经消失；两侧中点恰好都在实体时，`classify_candidate_pairs` 可走到 `IN_MATERIAL`（750-773）。

纯数学反例足以锁定此分支：`tol=1`，同一 face 命中参数 `[10.0, 10.9, 11.8, 12.7]`；单链合并总跨度 `2.7 > tol`，却没有 ambiguity reason。若交点交替进出实体，两个宽 `0.9` 的 VOID 区间本应是薄区间 UNKNOWN；平均点 11.35 变成唯一切点。把候选段其余部分取长、且该平均切点两侧中点落在实体时，当前结果可能变成 `IN_MATERIAL`。原始 hits 虽留在 JSON 中，决策却未消费这些独立交点。应对任意同组 distinct parameters fail closed（或保留各事件并证明每一间隔），并增加纯数学回归测试。此缺陷覆盖本审阅中的“全交点/切点/薄层/不得以采样或合并假 PASS”要求。

## P1 — asset UNKNOWN 准入门禁及 native 路径尚未实现

`classify_candidate_pairs` 只输出 pair-level `admissible_static_candidate`，总返回状态固定为 `STATIC_CLASSIFICATION_COMPLETE`，ledger 明示 `whole_asset_physics_status: NOT_ASSESSED`（814-871）。源码搜索未发现 `pgfrag` backend/assets/native 构造路径调用 `load_canonical_step` 或 `classify_candidate_pairs`；bond visibility 仅存在于该新模块。因而现在没有机制能在 UNKNOWN 落到关键壁/路径或超过预声明限值时拒绝整个 asset，也没有证明筛选先于 neighbor count、bond IDs、correction、kernel 和 storage。v1 作为隔离原型是清楚标注的；在这个门禁与拓扑重算完成前，不能据 pair ledger 把 vessel 标成 admitted 或启动真实求解。

## P2 — 精度与几何 provenance 是未验证的 caller assertion

`_source_from_shape` 只要求 `source_accuracy_mm` 是有限非负数、`source_accuracy_basis` 是非空文本（189-257）；没有把声明同 OCC 顶点/边/面容差、STEP precision、形状坐标量级或机器精度核对。于是错误地声明较小误差带时，宽于代码算出的 `2*(accuracy+tolerance)`、但窄于 B-rep 实际 tolerance 的结构可能被分类为 IN。应把读取到的 OCC topology tolerance 纳入误差上界；没有可信的来源精度证据则返回 UNKNOWN。现有 thin test 是调用者明确传入 0.01 mm 后验证该区间变 UNKNOWN，未覆盖“声明比实际 B-rep tolerance 小”的路径。

`load_canonical_step` 确实对文件字节求 SHA-256 并和调用者传入的 digest 比较（261-292），这比只验 64 位格式好；但没有仓内可信 manifest/source pin 绑定 canonical STEP 与真实候选点集。调用方若对错误 STEP 自己计算 digest 再传入，也会匹配。另有窄的 TOCTOU：第 287 行先哈希路径，第 303 行 CadQuery 再次按路径打开；并发替换文件会让输出记录一个 digest、实际导入另一个文件。将已哈希的不可变字节快照交给导入器，且显式绑定点集来源/资产 ID 与几何 hash，才能形成可审计来源链。当前测试只验证“非 STEP 内容 + 错误 digest 在 CAD 导入前拒绝”，没有有效 STEP 匹配导入或 asset-to-candidate 绑定正控。

## API 与几何判定核对

检查了仓内实际安装的 CadQuery 2.8.0 源码 `third_party/vesselforge-runtime/venv/Lib/site-packages/cadquery/occ_impl/importers/__init__.py`。其 `importStep(fileName, unit="MM")` 在第 97-99 行设置 `xstep.cascade.unit`；源码 docstring（92-94 行）明确这是目标 OCC 单位，OCCT 按 STEP 文件声明单位缩放至目标。因此传 `unit="MM"`、将参考点声明为 mm、并显式输出 `m = mm * 1e-3`（代码 303、835、853-855）这条单位链一致；需要注意 importer 改写 OCC process-global `Interface_Static`，并发 import 需要串行化或保护。

同一环境 OCP 包 metadata 为 `cadquery-ocp 7.9.3.1.1`，Python 类绑定由二进制扩展提供，安装目录没有这些类的可读 `.pyi`/C++ 源。检查了项目现有 OCC 调用模式，并对照 [OCCT IntCurvesFace_ShapeIntersector API](https://dev.opencascade.org/doc/occt-7.6.0/refman/html/class_int_curves_face___shape_intersector.html)：`Load(shape,tol)`、`Perform(gp_Lin,pinf,psup)`、`IsDone/NbPnt/WParameter/Face` 的调用形状正确；它与单独的 `PerformNearest` 是不同入口，当前用的是全结果入口。`gp_Dir` 将方向单位化，所以 `WParameter` 参数与长度 mm 对应。`BRepClass3d_SolidClassifier.Perform(point,tolerance)` 与 `State` 用法同仓既有调用一致。这里确认了调用面与 API 契约，没执行绑定，也不能据此证明 OCP 7.9.3.1.1 对本项目几何的运行行为。

`_occ_topology` 用 `BRepCheck_Analyzer.IsValid()`、solid leaf 检查、shell closed 标志筛形状（127-174），开放/非 solid/检查异常都会降为 UNKNOWN；单一 solid 默认归一个 domain，多 solid 要显式给 ID。这个 fail-closed 骨架合理。端点和 interval 的重叠 solid / ON / 未知状态也不会 eligible。`raw_hits`、参数、面/solid 索引保留（376-443、779-808）；两个有向记录从 canonical 无向 pair 同步生成，重复反向输入合并（667-688、775-785），id 不存在/相同/重复点坐标索引会报错。它们是有用的审计及 reciprocity 设计，但不能抵消聚类之后的决策丢交点问题。

## 测试覆盖与可执行下一步

现有测试代码有 hollow-shell air path（63-104）、实体路径（106-128）、分域 handle gap（130-143）、球切线与 ON 面（145-166）、薄壁精度阈值（168-181）、非 solid / precision missing（183-204）、pair reciprocity/id canonicalization（206-221）。这些都是**未执行**的 in-memory fixtures。设计中要求的 coupon negative control（mask 后整个 eligible pair set 与原 native graph 相同、blocked=0、step-zero damage不变）仍只是计划；测试只检查单个 coupon pair 为 IN，不能叫实测负控。测试也没有实 STEP hash-match/import、coincident/edge/vertex hits、同 face 近邻交点链、交点不足/失败、actual OCC topology tolerance 大于声明、asset UNKNOWN admission 或 native degree/correction/storage 顺序检查。

推荐可执行下一步：修复并锁定上面的纯数学聚类回归；随后只运行 hollow shell 和 solid coupon 的最小 CadQuery/pytest 门槛，验证 hollow shell 的四个交点与 `IN/OUT/IN` 区间，coupon 的 pair eligible 且 UNKNOWN/OUT 计数为零。这个小测试足以判断 OCC API 调用及 v2 核心分段路径是否可用，但不足以批准真实 vessel 入 solver。继续前至少还需真实 canonical STEP 的可信 digest/unit/frame/precision manifest、针对实际 vessel 的未知关键区域拒绝门禁、有效 STEP hash-match 的读入测试、薄层/切点/重合/同面交点链负控、以及 native integration 验证 `n_neighbors`、bond IDs、correction/kernel/storage 都基于筛后的同一张图。不要把这份源码审阅当作求解器验收。

## 数值边界与“无点/切线”检查

在 v1 固定源码中，`source_accuracy_mm + distance_tolerance_mm`（691-695）可能由两个有限数相加溢出为 `inf`。后续每个正宽 interval 都会满足 `width <= 2*inf` 并强制 UNKNOWN，因此按可见 Python 控制流它不会直接产生 `IN_MATERIAL`；但 endpoint classifier 仍先收到无穷 tolerance（717-718），没有在调用 OCC 前变成结构化 UNKNOWN。应显式拒绝非有限合并带，避免把无效容差交给几何内核。

`math.sqrt(math.fsum(component * component ...))`（701-706）可能在极大有限坐标差下得到非有限长度，也可能由 `fsum` 溢出异常；只有返回值已是非有限时才转成 UNKNOWN，异常没有捕获。cluster 参数均为非负，`math.fsum(parameters) / len(parameters)`（481）也可能在平均值本可有限时因求和中间值溢出而抛出异常。两类目前更像 fail-stop 而非直接假 PASS；使用 `math.hypot` 和稳定均值、并把数值异常映射成 UNKNOWN 会更符合 ledger 契约。

对 tangent/no-hit 分支，没有从静态源码证明当前 OCC 会漏掉一个有效相交。OCCT `Perform` 是全交点入口，`PerformNearest` 是另一个 API；代码也检查 `IsDone`。但若 `IsDone=True, NbPnt=0`，算法会把整段当一个区间、只测一个中点（568-626）；两端 IN 且中点 IN 就可放行。故安全性依赖 OCC 对共面、切点、共享边/顶点和重合面的完整交点结果。代码没用 `Transition()`/`State()` 作独立校验；球的外切测试两端在球外，不能单独证明“端点在材料内的擦边路径一定 UNKNOWN”。这些需在后续可执行 CAD 负控里实测。当前未把它们认定为已复现的源码假 PASS。
