# 下一阶段断裂路线：最小决策备忘录（2026-10-04）

**结论。** 下一步不迁移 solver：保留 Akantu 固定界面结果作失败基线，围绕现有 Peridynamics.jl/OSB 做一个隔离的“任意材料点集输入 + 每步 native 断键事件”实验适配器。它是当前最接近可验证自由路径的本地力学后端：native Body 可接收任意 3×N 位置和 N 个体积数组（构造器只检查维度与 NaN；正值、有限性和物理求积意义须由调用方验证），OSB native bond 依临界伸长断裂；PG wrapper 目前却只构造 uniform box，并只写终态断键。裂面重建、接触、能量闭合、真实陶瓷标定和 R2 血缘均未完成。**当前合格真实器物理碎裂样本为 0**；现有玩具能力不缩短这条资格链。

## 真实完成度与方法边界

| 现有路线 | 已实现且有证据 | 不能声称 / 当前距离 |
|---|---|---|
| Akantu cohesive FEM | 3D extrinsic 动态插入和指定中央界面双侧 native 连通性通过；仅在冻结 toy 网格上证明。 | 插入沿 bulk 三角面候选，方向仍受面网格限制（由插入接口推断）。两例 signed cohesive dissipation 负、Gc·A 未到，能量 gate **FAIL / NOT_REACHED**；contact energy selector 未验证、contact **NOT_RUN**。不应用来预测自由陶瓷裂纹。 |
| PG Peridynamics.jl OSB | 3D mechanics 已有；Body(material, position, volume) 接受通用点/体积；native 每步按 critical stretch 更新 bond activity。 | PG adapter 只用 uniform box、均质 toy 参数与轴向速度边界；broken_bonds.csv 是末态子集，无断键时刻。无裂面网格/paired positive-area lineage、无 contact case、外功边界下总能量守恒 UNKNOWN。PD 不受 tetra 面路径约束，但仍依赖点布局、horizon、求积/边界修正；自由裂面不能由断键端点唯一推出。 |
| SynLab shared cutter / implicit | 已有共同切割、碎片事件和正面积配对接口；可接既有几何/R2审计。 | cutter 是输入的观察/程序化表面，不由力学求解产生。观察扫描路线可承载**实测裂面**证据；程序化表面只能作为单列 SYNTHETIC 诊断层，不能冒充真实物理。 |

若未来重开 cohesive 自由路径，至少须另报薄壁 t/h 与 h/l_process-zone，并通过独立重网格方向测试；当前固定界面 toy 既无自由裂纹、也没有这些尺度证据。

Akantu intrinsic bilinear 201 点 run02 仅是已知中央界面、全 DOF 规定加载的 toy law 检查；300.046 s 超时，仅 168/201 点，端点能量 gate NOT_EVALUATED_INCOMPLETE。不得把它计为自由路径验证。已记录证据见 source hash 清单。

## 建议只做的下一实现单元

在 pgfrag/experiments/ 增加一个**实验专用、非 production** Julia 驱动器：输入冻结的 rest-space xyz + volume，并把 point/material IDs 以旁路清单绑定每列；通过现有 Body(mat, position, volume) 建 native 3D OSB body；在每个 native velocity-Verlet 力更新后，对 bond_active 与上步快照作差，写规范化 (point_i, point_j, first_break_step/time)。calc_failure! 已在 native loop 中更新该数组；现有 run_fixture.jl 已能导出 native 全局点号和终态 bond 状态。因此这是复用 mechanics、补齐 event adapter，不是重写 PD。native loop 没有已证实的用户回调；驱动器若调用内部 timestep API，须先用同一无断键短例证明其状态与原 native submit 等价。拒绝最近点、坐标猜 ancestry；点号必须绑定不可变输入清单。

下一步单元的 gate 仅是：同一冻结点集在短无断键例中与 native 状态一致；加载例记录每步单调、互反一致的断键 pair，首断步可复算，起始/终态 bond 数及质量闭合。任何 ID 丢失、非互反或 timestep 路径改写即停止；只报告 adapter FAIL/UNKNOWN，不导出 R2。暂不实现裂面重建或杯体运行。

## 自由路径方向偏差的最小验证（后续由 root 冻结门槛）

1. 用同一固定坐标系的 3D 各向同性预缺口试件，预缺口和载荷轴相对离散方向为斜交。几何、缺口、材料场、边界、加载历程、终止时间和随机场逐项冻结。
2. 固定物理 horizon δ 与材料/载荷，构造**三个独立生成**的点布局方向/seed；在每个布局上做两个 h 水平（共至少 6 个 case）。空间材料场必须在固定 rest 坐标上采样。分别记录 h/t、δ/h、δ/t。随后才独立做 δ 敏感性组，不能将 h、δ、材料和载荷一起调。
3. 把整个试件连同载荷一起刚体旋转只检验协变，不算独立布局。有效方向对照是在不动试件、缺口、载荷、材料坐标系时，重新生成/旋转采样格点并重新裁剪到相同实体。
4. 运行前由 root 冻结路径/起裂、反力—位移、断裂功、质量/能量残差、裂面面积/法向/拓扑的容差和终止条件。当前账本没有证据支持数值阈值；未冻结即 NOT_EVALUATED。每个配置需足够的事件时间与自由断面真值，不能只比末态组件或图像。
5. 任一布局方向差在细化后仍超预注册界、能量/材料 gate 失败，或真实裂面 ancestry 无法通过 R2 源唯一性合同，就停止该路线对该尺度的物理资格判断，标 NumericPhysical FAIL/UNKNOWN；不能通过改善外观补救。

此实验即便通过也只验证一个 toy 试件上的离散公平性，不是粗糙陶瓷或器物断裂认证。薄壁实物尚须实测 t/h、标定 δ/t/材料参数和冲击/接触；真实陶瓷不同材质需分别标定强度/断裂能及不确定度。R2 仍需 parent、fragment、event、paired positive-area 新断面、boundary 与 pose 血缘，并通过现有 fracture_truth.md 的来源唯一、材料内部采样、面积、gap 与法向门。 若未来满足物理门，最小 annotation 接入预期是隔离的 PG→SynLab FractureResult adapter 为既有对象字段供证据；不改 R2 schema/阈值。没有可追溯 ancestry 和面积裂面时不接 label。

## 止损备用路线

若 native event adapter 或后续 PD 数值 gate 不过，不切到本轮已能量失败的 Akantu，也不另装 solver。用已有 SynLab shared-cutter 只接**实测断口扫描/观察面**，保留母体到 fragment/event/配对面的可审 lineage，作为观察证据→历史相邻实验路径；若只有程序化 cutter，则 SYNTHETIC 隔离。没有实测表面且无物理通过样本时，历史 R2 留 UNKNOWN，不产出物理标签。此止损只保护数据路线，不替代用户要求的真实物理目标。

## 本地证据索引

源路径、字节数、SHA-256 与支撑行号见同目录 source_hashes.json。本审计未运行 solver、未改生产/已有报告、未读动态 run03，也未用外部技术资料；mesh/horizon 偏差比较依本地研究 ledger 的既有预案，具体验收数值须 root 另行冻结。



