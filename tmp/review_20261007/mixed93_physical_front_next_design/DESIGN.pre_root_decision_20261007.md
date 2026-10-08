# 93-DOF 固定裂面到首个前沿候选：只读设计

状态：DESIGN_ONLY；本轮只读源与 JSON/Markdown 元数据，只写本目录的 DESIGN.md、PLAN.json。没有运行 solver、provider、NPZ 或测试。本设计不依赖一项尚未产生的 93-DOF 数值成功结果，也不扩大 root 当前固定裂面任务。

## 结论与当前边界

建议把 root 当前工作与前沿推进分成两张实验卡。当前卡仍是单一固定直裂面、virgin-elastic Mode-I CZM 耦合：最多 16 个新 MLS 查询、3 次 cohesive primitive assembly、1 次 99 阶 KKT solve。它可以建立固定几何下的真实耦合输入/输出契约，但不能同时证明旧面与新面都已平衡，也不能给出裂纹增长。

后续的单步卡必须多一个候选几何的机械重分析；若没有可复用且完全匹配的旧面平衡 receipt，就需要旧面、候选面各一次平衡求解。新卡应明确重新预算 MLS/装配/求解，不得把当前 16/3/1 当成两态增长预算。提交对象始终是一个有界几何 candidate，而非自动裂纹速度或下一步循环。

建议只提供两个物理计账分支，且一次运行只能选择一个：

1. Griffith/LEFM 分支：旧、新态无惯性、无其他不可逆耗散；新增裂面能量恰为 GcΔA。能量门使用 v2 已声明的公式。
2. Cohesive-history 分支：新增损耗由实际 cohesive history 的 ΔDcoh 给出，能量门应是 Wext − ΔUstored − ΔDcoh（再显式列出任何其他耗散）。当前 v2 的 GcΔA 门不支持该分支；不得同时把同一 cohesive 断裂耗散再加一次 GcΔA。

目前 93-DOF cohesive case 的参数 K0=1e8 Pa/m、tpeak=1000 Pa、Gc=0.01 J/m² 是 toy bridge 量，不是材料数据；不能用它们接受真实裂面推进。没有真实材料 Gc、真实物理外载、current J 字段及机械重分析 receipt 时，终态只能 UNKNOWN/NEEDS_RECOMPUTE，active old surface 保持不变。

## 只读证据与含义

- 93-DOF 缓存试验说明书为 pgfrag/experiments/mixed93_cached_lefm_solve_20261006/README.md，case 为同目录 case.json。布局是 27 个 BASE source × 3 位移分量加 12 个 GLOBAL_FRONT 分支自由度，总计 93；当前几何是 y=0 平面、x∈[-R,0]、z 方向有限直前沿。E=1 MPa、ν=0.3 与 KI/KII/KIII manufactured fields 被明确标为 toy numerical benchmark。该源唯一执行曾以 JSON 序列化错误 rc=1 结束；JSON 独立恢复审阅报告缓存数组 internally consistent、physics unqualified。这里仅读取文字和 JSON 元数据，没有碰 NPZ 内容。
- 当前 cohesive bridge case 为 pgfrag/experiments/mixed93_actual_cohesive_bridge_20261006/case.json，case SHA-256 为 fcc4b46ed101b0f6067a9f3ec10816909c9fa931edd5a8cdde4c7f46cb3f3a71，driver SHA-256 为 023d730454802ec755df41c9878744aa086cd5b97a8a88e9d5bbd4181b6f1412。它是 SOURCE_READY_EXECUTION_DISABLED，真实数值桥未实现/未审阅；16 个 y=0 界面点、最多 16 个新 MLS 查询、最多 3 次 primitive assembly、一次 KKT solve。它的四支 GLOBAL_FRONT 分支在 F1 上给出当前直前沿跳跃，且法向投影仅命中指定 F1-uy 系数。这个 GLOBAL_FRONT12 布局是一个已声明单平直 front 的固定前沿 benchmark；它不包含任意工件多裂纹、多 front、分叉/交汇的 basis、所有权或更新算法。拓扑变成多 front 时不能保持 93 系数声明不变却宣称支持增长。
- scalar cohesive interface assembly 的 README 和源模块 pgfrag/experiments/cohesive_interface_assembly_20261006/ 说明：给定同一 QP 的 H+、H−、跳量和 immutable committed history，primitive 只做 read-only trial，返回力、切线、stored/dissipated/work 累计量；每侧面积只计一次，历史键须和 QP IDs 完全匹配。模块仅支持单法向 Mode I，压缩、切向/混合模式是 UNKNOWN；softening tangent 可为 indefinite，kink 为 nondifferentiable；没有全局 solver、front coupling 或材料校准。接口给出的 trial energy 不能单独代替旧/新平衡态间的外力增量功。
- pgfrag/experiments/explicit_surface_state_20261005/README.md 与源契约要求 CrackPatch 两侧有数值位移数组，surface state 固定 NOT_COUPLED / mechanically_active=false；surface/patch hash 和 parent ID 是几何 lineage 元数据，不证明几何包含关系、全局不交、平衡或物理资格。提议状态可先存 geometry-only strip，但不能补零伪造 candidate 位移；只有候选面重解产生的两侧 displacement arrays 才可用来构造重新平衡的 ExplicitSurfaceState。
- pgfrag/experiments/griffith_front_step_20261006/README.md 的单步核只负责 straight finite Mode III、预声明 +x 方向、bounds、一次 propose/confirm/reject；G=Gc 持平，G>Gc 只是提出候选；source/status/hash 是调用者 receipt 声明，不验证 solver 真实性。其 v2 能量门 pgfrag/experiments/griffith_front_energy_gate_v2_20261006/README.md 明确不桥接含义不清的旧 delta_potential_J。唯一 supported mode 是无惯性、无其他耗散的 quasi-static elastic；公式是 residual = Wext − ΔUstored − GcΔA，ΔΠ = ΔUstored − Wext，ΔΠ 只作审计量，不可在 residual 中再减一次。abs_tol、relative_tol、正 reference_energy、单位和 parent/candidate/reanalysis/hash binding 均需实际数值输入；缺项或超限保留 active old surface。外部 current reanalysis 可得到 COMMITTED_EXPERIMENTAL_REANALYSIS，但这不是 physical validation；fixture 仅能 COMMITTED_FIXTURE_ONLY。
- 已存的固定载荷有限变化记录都只是合成诊断：tip-aligned Mode-III 记录与独立 JSON 审阅给出 δa=0.5 mm 时的 G≈0.01411668 / 0.01411683 J/m²（theta32/64），但物理断裂仍 UNKNOWN、G 未与材料阈值比较，角向/径向积分没有认证误差界；移动 tip 后旧 tip manufactured field 已非 exact。fixed_load_front_reanalysis_20261006 的 result SHA-256 为 a60e51c85c26510c358ee57cafa4bfdda6ed89ea463d99bcc86551f64220732c，status FINITE_CHANGE_RELEASE_DIAGNOSTIC_ONLY_NOT_COMMITTED，δa=0.5 mm、单面投影 ΔA=5e-6 m²、ΔΠ=-2.7331683517067647e-8 J、G diagnostic=0.005466336703413529 J/m²、未比 Gc/未提交；记录注明候选 tip/branch 未采样，独立评论提到约 0.013 的 finite-G 差异与旧 tip quadrature 风险。这些不同构造的数字不应互相拟合或用来给 93 toy Gc 定标，也不是当前 93 case 的真实 J。

## 建议的最小接口与单步主线

令 parent surface revision 带 immutable geometry hash、实际旧面两侧位移、物理外载和材料来源绑定。将以下三类 receipt 明确分开：驱动力证据、候选 reanalysis、能量门。每个数值都带 SI 单位；hash/status 是 lineage 与 provenance，不是物理通过证明。

1. **当前驱动力 receipt**：实际 J 组件对 parent exact field 产出 G[J/m²]、J result ID、J 源 hash、solver-result hash、field/source hash、parent surface hash、当前 coupled 状态，以及 contour 与 path sensitivity 的独立状态/证据 hash；同时传真实材料 Gc[J/m²]、来源/校准 hash。source component 必须确认 mode、前沿和 +x 方向与本模型一致。缺失、未耦合、过期字段或未检查敏感性 → UNKNOWN/NEEDS_RECOMPUTE。不能由字符串 REVIEWED 或 qexact 构成数值驱动力。
2. **受限 propose**：只对已验证的一条有限直 front 沿 +x 比较 G 和 Gc。G<Gc → ARRESTED；G=Gc → HELD_AT_THRESHOLD；G>Gc → 只产生一个 geometry-only strip candidate，并进入 PROPOSED_WAITING_REANALYSIS。δa 是调用者要求的几何试步，经 finite-positive、step cap、seed frame、domain bounds 限制；它不是由本控制器猜出的速度/dynamic law。候选 ID 绑定 parent hash、geometry hash 与 step sequence。此时旧 active surface 不变。
3. **前沿改变后重建模型**：以候选 tip 重新定义 branch origin/方向，重算所有受影响体/外边界/界面采样处的 basis values 和 gradients、MLS supports、积分与 KKT 装配，记录 basis generation/hash。旧解系数不能复制到新 basis 当作解，q_new 与 q_old 也不能逐项相减，因为两态语义基可能不同。旧位移只可经声明的物理空间映射作为迭代初值，最终必须重新满足候选几何平衡；旧面/新 strip 上的 side displacement 都由候选解重建。历史只按精确稳定 QP ID 和几何 lineage 映射；新 QP 只有有明确本构规则时才初始化 virgin，否则 UNKNOWN。
4. **相同真实外载的两态机械再平衡**：必须绑定同一个真实物理 load definition/hash、同一个物理外部边界/控制目标与材料版本。候选前后实际边界牵引不重标定、不按新 KI 缩放。先使用可复用的 exact parent reanalysis，或在缺失时对旧几何求解；再对 candidate geometry 单独求解。返回两态真实 displacement/strain/traction/history、线性/非线性收敛和稳定支路信息、候选 state hash、solver/source/result/field hashes、geometry/basis/load/history hashes。若只重算一个面、解不收敛、branch 不稳定或无法验证 hash binding，candidate 留待重算而旧面保持 active。
5. **完整能量/面积账本**：在物理空间而非不同系数向量上求外功增量；至少记录 ΔUbulk、ΔUcohesive,recoverable、ΔDcoh、Wboundary、可能的 Wgauge/support、ΔA、GcΔA、units、积分规则/误差依据及每项来源 hash。ΔA 对单个内部裂面按裂面参考面积计一次；直 front 试步的几何值是 front_length×δa，不因 H+、H− 两张位移侧而乘二。区分旧/新端点总 potential 与路径 work，禁止用歧义 delta_potential_J 顶替 ΔUstored 或 W。
6. **按所选物理模型确认**：
   - Griffith v2 分支只在无惯性、无 cohesive/其他不可逆耗散（或这些耗散已独立证明为零）时允许确认。令 Ustored 含全部可恢复的 bulk/interface elastic 能；门为 r=Wext−ΔUstored−GcΔA，ΔΠ=ΔUstored−Wext 仅派生记录。实际 abs_tol_J≥0、relative_tol≥0、reference_energy_J>0 必须有限，且依据/哈希对应当前离散和误差估计；超出 abs_tol+relative_tol×reference_energy 或缺项即 NEEDS_RECOMPUTE。
   - 若 cohesive history 发展出 ΔDcoh，不能在该耗散之外再减同一新面的 GcΔA；只有经独立推导的 cohesive ledger 才可能启用另一个能量核，当前 v2 直接 UNKNOWN。
   - 动态/惯性、mixed-mode cohesive、软化 snapback、非凸解未证稳定、非线性未收敛均为 UNKNOWN/NEEDS_RECOMPUTE；不能用“G>Gc”跳过能量门。
   - 新旧态若使用非零约束目标，93 KKT 恒等式中有 g·lambda：2U−W+g·lambda=0。它会影响能量闭合。必须预先选择同一零目标纯 gauge 并证明无约束 work，或把约束/支承反力沿实际加载路径作成有符号 Wgauge，纳入 ledger；不能默默丢弃 g·lambda，也不能把它作为数值残差阈值来“通过”。当前 bridge 的 gauge target 来自缩放 manufactured q_solution，不能直接用作真实断裂能量控制。
7. **终态**：完整、同候选绑定、外部 REANALYZED_CURRENT reanalysis 与纯 Griffith gate 通过时，只允许 immutable state update 标记 COMMITTED_EXPERIMENTAL_REANALYSIS / EXPERIMENTAL_UNCALIBRATED；这只说明外部 receipt 数值符合声明的实验门。外部 receipt 由独立 adapter 实际验证前，状态机中的 hash 和状态词仍是声明。任何 rejection、UNKNOWN、NEEDS_RECOMPUTE 都保留 parent active。一个 turn 后停止，无自动 propose 第二步。

## 状态图

ACTIVE_OLD → ARRESTED（当前 G<Gc）或 HELD（G=Gc）或 PROPOSED（G>Gc，只生成一个有界 geometry）或 UNKNOWN（驱动缺失、过期、未耦合）。
PROPOSED → ACTIVE_OLD（明确 reject）或 NEEDS_RECOMPUTE（候选重解/账本不完整或超限）或 UNKNOWN（dynamic、mixed mode、软化不稳定）或 COMMITTED_EXPERIMENTAL_REANALYSIS（同候选 current reanalysis 与所选能量门通过）。所有非 commit 分支都保留 parent active；一次后停止。

## 可独立并行准备的工作与 Astra 决策

可单独准备而不写现有生产代码的工作包：
- J 输入适配契约：把真实 J solver field、路径/轮廓 sensitivity 及状态 hash 接到 parent geometry，并拒绝旧字段/错 mode。
- 移动 tip basis 构造预算：枚举候选后必须重新求值的 BASE MLS 与 GLOBAL_FRONT branch 体/边界/界面坐标，量化新查询/装配/存储上限；说明什么缓存只可复用。
- energy adapter：外部真实 load 物理功、stored energy、cohesive dissipation、单次面积 convention、quadrature uncertainty 与容差来源的独立账本及 sign review。
- reanalysis binding：几何/历史变更后解出全部 old/new patch displacement fields 的 external receipt 验证，以及稳定支路/solver convergence 证据。

实施前需要 Astra/root 决定的真实困难是一个不应靠代码猜测的模型选择：第一步究竟采用 (A) 新增 traction-free 面的 LEFM/Griffith GcΔA 门，还是 (B) cohesive damage history 耗散门；现有 93 case 是 toy、elastic-only cohesive bridge，而旧 v2 明确只实现了 (A)。同时必须裁定 93 的非零 gauge target：改成同一零目标纯 gauge并重做真实 fixed-load baseline，还是证明并记录真实支承/约束路径功 Wgauge。二者未定前，机械再平衡即使数值残差小，也不足以作裂纹推进能量判决。

明确未实现：分叉、转向、起裂、任意多裂纹/多前沿、fragment-count 推断、速度或 dynamic law、mixed cohesive、snapback continuation、真实材料校准、production surface promotion。任何受限 mainline 的输出都不能预设 fragment count。
