# PD 加载 coupon：最小判别测试设计（2026-10-04）

**结论。** 现有 tinyload 配置几乎不可能触发 native OSB 断键：在冻结材料数值下，临界伸长约 0.161%，而它 5 步仅给 4 mm 加载点中心距约 0.000023% 名义应变，约为阈值的 1/7,000。它只能作极小弹性/事件观察器冒烟测试，零事件不证明断裂实现。建议 root 从下列两个独立加载 gate 里选择、重新冻结新 case/common/manifest 后再决定是否放行；本报告不运行或改变任何 solver/config。

**已知状态与范围。** 主审已逐点验收同一 125 点、thread=1 的零载 case 与实际 PNG；该结果只证明 native observer 流程及无载活动不变。当前读取的 runner/common/tinyload 三文件与 pgfrag/tmp/review_20261004/pd_point_event_probe/finalfreeze_v3 副本字节一致，SHA 见 JSON。tinyload 配置自身仍是 5 步来源准备态；不改写它。本文只读配置、初始点集及 pinned Peridynamics 源，不运行 Julia/求解器/CAD/Git，也不改代码。

## native 阈值、dt 与尺度

native OSB 把 E、ν 转成 K=E/[3(1−2ν)] 与 G=E/[2(1+ν)]，再按 Gc 和 horizon δ 算 CriticalStretch 的 εc=√(5Gc/(9Kδ))；每条允许失效的键仅在当前伸长 (l−L)/L > εc 时被置为 inactive。代入 UNVALIDATED_PROTOTYPE 参数 E=5 GPa、ν=0.24、Gc=45 J/m²、δ=3 mm、ρ=2100 kg/m³：

- K=3.2051 GPa，G=2.0161 GPa，εc≈1.61245×10⁻³（无量纲，0.161245%）。
- 以加载点中心面间距 Lg=4 mm 为均匀轴向估算，达到 εc 需总相对位移 εc·Lg≈6.45 μm。若 top/bottom 速度分别为 +v/−v，名义 ε≈2vT/Lg；这只估计总体加载量，不预测每个非局部键的真实局部伸长。
- 由这些各向同性常数算出的连续体参考纵波/横波速约 1675/980 m/s；5 mm 域的纵波通过时间约 2.99 μs，4 mm 加载面间约 2.39 μs。它们是工程波时标，不是已验证的 OSB 色散或准静态判据。
- native 稳定步公式是 dt_i=√(2ρ / Σ_j(V_j·bc/L_ij))，bc=18K/(πδ⁴)，每 chunk 取最小再乘 safety=0.8。按 1 mm 立方 cell-center 点集、半径 δ=3 mm、最少邻居角点作独立求和，若包含 L=δ 的邻居，估算 accepted dt≈0.915 μs；若邻域搜索排除正好等于 δ 的点，约 0.987 μs。实际运行必须以 native preflight/accepted dt 为准，不能把这一区间替代运行记录。
- 点中心坐标是 0.5–4.5 mm，top/bottom 两个 25 点加载层相距 4 mm；5 mm 是带 cell 体积的实体包络，不是 gauge length。δ/h=3，δ/5 mm=0.6。5 mm 立方体内切半径只有 2.5 mm，小于 δ，因此每个点的 horizon 都被外表截断。这是强 surface-effect API/coupon 诊断，不能代表 bulk 材料或薄壁器物。

当前 tinyload 的 v=1e−4 m/s、5 步对应 T≈4.57–4.94 μs、总面间增量仅 0.915–0.987 nm、名义应变约 2.29–2.47×10⁻⁷（εc 的 0.014–0.015%）。按同一名义应变要走到阈值约需 32,700–35,300 步，而 common_case 仅允许 5 步。因此 tinyload 无 transitions 是预期的无载/微小响应链路观察，分类应为“弹性诊断、断裂未到达”，不是 crack PASS，也不能因无事件去调大现有 case。

## 两项相互独立的加载 gate

两项共用同一 125 点均匀 voxel 输入、OSB、SI、δ=3 mm、ρ/E/ν/Gc、native AUTO dt、thread=1、single chunk。source runner 已逐步记录 step 0 与每步 point displacement/velocity/acceleration、b_int/b_ext force density、damage/active degree，以及每条 directed 和 canonical bond activity/transition。建议新 common 的 max_steps 仅增到 320，保留 120 s timeout；分别新建 case TOML 与 manifest 锁定参数/输入哈希。现有 common/tinyload 仍保持原样。启动编译时间和加载用时未测；零载的单步成功不能外推耗时。

| 独立 gate | 可执行参数建议（写入新的冻结 case，root 决定后执行） | 可证伪预期与边界 |
|---|---|---|
| A. 加载弹性/no-break 阳性响应 + 零载阴性控 | 阴性控复用主审已验收的 1-step zero-load 证据，不重跑。新增阳性加载用 axis=3、top/bottom ±0.005 m/s、64 步、native AUTO dt；当前速度 BC 是恒速函数，无 ramp。 | 均匀估算面间增量 0.586–0.632 μm、名义应变 0.0146–0.0158%（εc 的约 9.1–9.8%）。期待非零位移、速度与 b_int 响应，所有键仍 active；从 rest 坐标和逐步点位置另算最大键伸长，运行前冻结其上限（可先用 0.5εc 作保守诊断门，不作物理常数）。若响应全零/映射或 dt 错，observer/load 路径 FAIL；若有任何断键或最大键伸长越界，no-break 控制 FAIL，停止用该载荷诊断阈值。它仍是骤启后的显式动力测试，不是已证明的准静态弹性模量测试。 |
| B. native critical-stretch crossing event diagnostic | 同一 fixture；axis=3、top/bottom ±0.02 m/s、320 步、native AUTO dt。不预删键、不加固定断面/碎片、不做参数搜寻。 | 均匀估算末态名义应变 0.293–0.316%（约 1.82–1.96εc）；名义阈值到达约在 163–176 步，1.5εc 约 245–264 步。真实首次事件依 native 局部键伸长，不保证落在名义步附近。必须留下 step0 至 320 的完整 force/displacement/velocity/activity 与所有 transitions，复核首个真 active→inactive pair、时间戳、互反键状态和阈值伸长。无事件则 NOT_REACHED，停在本冻结预算；早期突发事件、质量/映射错误、非有限状态或单步大量键塌陷则 STOP/FAIL。通过只证明指定 toy 参数下 native threshold crossing 被观察到，不证明物理断裂。 |

constant velocity BC 在 step 0 之前为零，native 每步却把 velocity_half 设置为 ±v，故 t=0 开始存在速度阶跃和应力波。B 即便有约 50 个 gauge 波通过时间才达到名义 εc，也仍是显式动力事件诊断，不能按步数或小速度直接叫 quasi-static；无阻尼/平衡证据时不能把 hold 的短暂稳定当静力平衡。若 root 要测准静态路径，应另行授权并冻结一版平滑 ramp 的实现与 hash；当前 runner 只从 TOML 读取常速，不支持 ramp，不在本审计扩改。

## 成本、反力与能量边界

按 125 点规则立方格对 L≤3 mm 组合求和，初始图约 6,818 directed / 3,409 canonical 键（正好位于 δ 的边由 native 搜索语义决定，step-0 native 计数才是权威）。observer 每个样本保存 125 点行和约 10,227 条键状态行；A 的 65 个样本约 0.665 百万键状态行，B 的 321 个样本约 3.283 百万，且 runner 在内存中暂存向量行。它比 5 步例扩大几十倍至数百倍，存储和峰值内存不可忽略；没有实测时间或 ETA。120 s 到时就记录 timeout/NOT_REACHED，不自动重试、不扩预算。

b_int/b_ext 字段单位是 N/m³，必须乘各自 point volume 后才是 N。恒速受控边界的 actuator reaction 没有直接 native 列；动态估计需要以受控点离散动量变化/Δt 减去积分后的 b_int+b_ext，并核对 Verlet 时间错位和符号。不可把某一面的 Σb_intV 直接称为 actuator force。由该反力与边界速度积分的外功、Σ½ρV|v|² 动能、以及每步 native strain-energy density 积分，才可能组成加载能量账；当前 probe 不写总反力、外功或 strain-energy density，能量状态仍 UNKNOWN。零载能量平衡是退化检查，不代表加载能量 gate。

另一个独立本构风险：native OSB 每次 force_density 更新都会以当前 bond_active 重新求 weighted volume m_i=Σ_active(ωβL²V_j)；dilatation、力和 strain energy 随活动族与 m_i 一起改变。断键后不能以固定初始 m 解释响应，也不能把 native strain-energy 的跳变直接叫作断裂耗散。任何能量 gate 前，须证明逐步 m_i/能量定义、驱动功与动量账本闭合；现有输出不足以判 PASS。即使阈值公式把 Gc=45 映射到 εc，该 Gc 也是 UNVALIDATED_PROTOTYPE 的输入，不是陶瓷试验标定，也不是模型自证真实韧度。

## 止损与后续门

本轮合理终点只有 observer/load/event 诊断。A 失败则先不跑 B；A 成功后 root 可独立决定是否冻结 B。B 未在 320 步内触发则记 NOT_REACHED，不在同 manifest 改 v/Gc/δ 找成功；触发但能量、表面或材料未知仍保持物理资格 UNKNOWN。事件通过不代表 mesh/point orientation 无偏：方向/空间 layout、cross-resolution、δ/h 和不同材料标定另设 gate。也不代表 contact、3D 裂面、fragment pose 或 R2 真邻关系已通过。
