# Energy ledger 与 bias protocol 独立只读审阅

## 结论

保存的 energy-ledger 用例为 7 passed，bias-protocol 为 6 passed；本次没有重跑。energy-ledger 的 signed D、正负外功公式、首样本增量基线、缺 contact/clock 时 UNKNOWN、all-DOF N/A、GcA 只作 constitutive reference 且 calibration 恒 UNKNOWN 均有直接实现。bias protocol 会把旋转、相同 artifact、相同旋转不变量排除出独立候选；哈希差异和 seed/direction 差异只给 UNVERIFIED_INDEPENDENT_CANDIDATE，未升级为实际独立布局或运行结果。

发现三个可静态复现的 energy false-pass 输入（转换/汇总溢出、bool 容差、重复 required key），一个可复现的 bias 假“不具可比性”（SHA 大小写），以及单样本自由能平衡的空区间 PASS。下面给出行号和输入。

## 审阅证据

- energy-ledger freeze 的 integrity.json 为 FROZEN_MATCH；live source/test/raw audit input 的 SHA-256 与冻结副本和 pytest 前后 receipt 一致。保存输出为 7 passed in 0.02s，退出码 0。
- bias-protocol 的 test_receipt.json 记录 source、test、protocol 文件在 pytest 前后未变；当前 source/test 哈希与 receipt 相同。保存输出为 6 passed in 0.04s，退出码 0。
- 只读检查了两个模块、tests、energy frozen input/receipt、bias protocol 与 test receipt；没有运行测试/solver/Git。没有写生产或测试文件。

## 实质问题

### P1：有限能量输入可在 J 转换或 balance 算术后变成 inf/NaN，并误报 PASS

位置：src/pgfrag/energy_ledger.py:167–178 对原单位浮点值做有限检查；202 行换算成 J 后没有再检查；539–556 行做角色求和、首样本差分、残差、scale、limit 时也没有 finite guard。abs(residual) > limit 对 NaN 为 false，对 inf > inf 也为 false。

全有限输入复现（其余 metadata、EvidenceRef、physical clock、boundary-work definition 按 tests/test_energy_ledger.py 的 analytic positive control 填齐，并使用正常 finite tolerance）：clock [0,1] s；K [0,0] J，stored [0,1e308] J，D [0,1e308] J 且 cumulative=True，外功 [0,0] J 且 cumulative=True。末点实际 residual 为 2e308 J，float 汇总得到 inf；scale 和 relative limit 也为 inf。因 inf > inf 为 false，balance 标 PASS，D audit 也 PASS，overall 可标 PASS_ENERGY_AUDIT_ONLY，尽管外功为零、储能与耗散增量均为正。

转换后产生 NaN 的变体：stored 使用合法单位 kWh、原始有限值 [1e308,1e308]，D 使用 [0,1] J，其他项为零。乘 3.6e6 后 stored 两个样本均为 inf，差分 inf-inf 为 NaN；末点 residual/scale 为 NaN，失败比较仍为 false。应在单位换算后、聚合、差分、残差、scale 和 tolerance 构造各阶段拒绝非有限结果；任何非有限 residual/limit 均不能得到 PASS。

### P2：bool 被当成能量和 tolerance 数值；bool 相对容差会让任意 balance 通过

位置：_as_finite_tuple 的 float(value)（171–177 行）接受 True、False 和 "1.0" 这样的数值字符串；_validate_tolerance（215–217 行）同样只转 float，因此 balance_rel=True 被接受为 1.0。balance 的绝对 residual 按三角不等式必不大于 sum(abs(delta terms))；相对容差为 1 时允许值至少等于该 scale，故任意不闭合轨迹都能通过。

复现：两样本、其他 qualification/evidence 完整的自由动力学 ledger，K [0,0] J，stored [0,2] J，D [0,0] J，外功 [0,1] J；正确 residual 为 1 J，应失败。只将 EnergyTolerance.balance_rel=True，实现允许值变成 3 J 并输出 PASS。原始 bool 通道值也会转换后被读取为 0/1 J。应拒绝所有 bool 数值，并明确决定是否拒绝字符串数值，而不是静默转换。

### P2：非 bool 的 cumulative 值按 Python truthiness 冒充累计历史

位置：EnergyChannel.cumulative 声明为 bool（84 行），但没有运行时类型校验；耗散下降审计用 if channel.cumulative（304 行），自由平衡累计性门控用 not history.channels[name].cumulative（519–525 行）。因此 cumulative="false" 是 truthy，会被当作累计 D/外功接受；输出也原样写出字符串。JSON 配置把 "false" 用作布尔字面值时会把非累计历史错误升级为 balance 输入。应要求严格 bool 或显式解析并拒绝非布尔值。

### P2：required balance 通道重复列出时会重复累计该角色，能把失衡算成 PASS

位置：required-role 列表从 required_balance_channels 原样迭代（495–503 行），并在 539–542 行按该列表求和；Qualification 的 Sequence 没有去重校验。

复现：K [0,0]，stored [0,2]，D [0,0]，外功 [0,1]（D/外功均正确标累计），required 为 [kinetic, stored, dissipated, external_work, external_work]。无重复时 residual 为 +1 J，应失败；重复外功 key 后 work_names 含两次外功，增量变成 2 J，残差被抵消成 0 并 PASS。应拒绝重复 required key，而不是静默重复计量。

### P2：只有一个样本时，自由能 balance 产生空区间 PASS

位置：_as_finite_tuple 只要求至少一个样本（175 行）；audit_energy_history 允许 balance 用单样本计算（539–556 行）。给自由动力学所有必需角色各提供一个有限样本，完整 caller EvidenceRef、physical clock 和 boundary work metadata；所有首样本差分按 baseline 定义为 0，residual 为 0，结果就是 PASS / PASS_ENERGY_AUDIT_ONLY。这证明不了任何时间区间中的能量闭合。自由能 balance 建议至少需要两个时间样本，否则返回 UNKNOWN；单样本仍可供端点 D / constitutive-only 比较使用。

### P2：SHA-256 大小写有效但物理案例签名比较区分大小写，导致假 NOT_COMPARABLE

位置：bias_protocol.py:46–47 接受大写十六进制；_physical_signature / _same_physical_case（513–527 行）直接比较原字符串，未将 hash 小写规范化。复现：取 tests/test_bias_protocol.py:63–100 的 manifest 并复制为第二 study，只把第二份 physical_case.boundary_load_sha256 从 lowercase 改为同一 digest 的 uppercase。两份 hash 均合法且表示同一 SHA-256，但 _same_physical_case 返回 MISMATCH，整体被标 NOT_COMPARABLE。应在物理字段签名前规范化 SHA 字段大小写。

## 重要边界与接入风险

- 外功与基线：balance 用 ΔK + Δstored + ΔD - ΔW，外功正方向为 positive_into_system；正负外功保持有符号。所有 *_abs_j 阈值按 J 使用，channel 值先按单位转成 J；正常 finite 区间内公式和量纲方向一致。基线明确是窗口首行，没有偷偷假设初值为零。
- 来源、接触和校准：独立平衡 basis、channel EvidenceRef 与自由动力学 source/method 只标 VERIFIED_BY_CALLER_REFERENCE，不能核验实际文件或测量方法；BoundaryWorkDefinition.source 仅需非空字符串。required channel 列表决定哪些项被视为完整，包括 cohesive contact，模块不会从模型自动推导。材料 calibration 固定输出 UNKNOWN；MATCHED_REFERENCE 只表示 caller-supplied final D 与 GcA 的算术差在指定容差内。source-grounded 两个负 D 用例仍保留符号并 FAIL；小负 roundoff 保留并标 PASS_WITH_ROUNDOFF。
- Bias readiness 的含义：rotations、相同 artifact SHA、相同 rotation-invariant digest 不会算独立生成；不同 SHA/seed/direction 的候选标 UNVERIFIED_INDEPENDENT_CANDIDATE。READY_NOT_RUN 是预注册设计矩阵状态；出现 run evidence 后转 RUN_EVIDENCE_PRESENT_ANALYSIS_PENDING。这些哈希是声明，模块不从 mesh/point payload 重算。
- PD 尺度：horizon_m 与 process_zone.length_m 使用独立字段；任何缺失都令 peridynamic_resolution_qualification 为 UNKNOWN，不会从 dx、t 或最大边推算。两者值和 process-zone 来源仍是调用方声明，known 状态最多为 DECLARED_KNOWN_UNVERIFIED。另需主代理按实验验收条件判断：direction/spatial pair readiness（769–817 行）只对比 dx/h、dt、direction、generator，不要求 PD runs 的 horizon_m/process-zone 相同；dt pair（840–843 行）才要求二者相同。若方向或空间比较必须固定这两个物理尺度，则不同尺度的 pair 仍能标 READY，属于混杂设计。
- 裂面指标：没有 positive-area crack-surface metric；metric 名称包含 crack_direction，surface 只是泛型 qualified_surface。surface 缺失或 receipt/hash 链断裂时值为 null / NOT_RUN。若调用方自填 QUALIFIED 和匹配哈希，结果最多为 PROVENANCE_LINKED_UNVERIFIED；模块不读/重哈希 surface bytes，也不检查它代表 positive-area fracture surface。不能据此声称获得正面积裂面或可信裂纹方向。
- NaN：energy 通道的原始 NaN/Inf 会在初始 finite 检查中拒绝；缺口是在单位变换或算术产生非有限数后没有复检。bias 的 canonical_sha256(... allow_nan=False) 会拒绝 manifest 中的 NaN/Inf 并令 manifest invalid；_positive_number 也拒绝 bool、NaN/Inf。

## 保存测试覆盖

energy tests 覆盖 signed 大负 D、小负 roundoff、contact/clock 缺失、单位拒绝、全约束 GcA 参照、合成三点分析闭合、status-only evidence 拒绝；没有覆盖有限输入算术溢出/NaN、bool 数值或 tolerance、累计标记类型、重复 required key、单样本 free-balance。

bias tests 覆盖旋转/相同 artifact 不冒充独立生成、物理身份缺失和材料不匹配、geometry-only metric 和零值不等于 PASS、hash-linked metric 仍为 unverified、预注册矩阵与 PD resolution UNKNOWN。没有覆盖 physical hash 大小写语义、PD paired contrasts 的 horizon/process-zone 控制、positive-area surface 类型。
