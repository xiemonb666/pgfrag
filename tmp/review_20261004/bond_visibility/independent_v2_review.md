# Bond visibility v2 独立复审

只读审阅冻结字节；没有运行 Python/CAD/STEP/pytest/Julia/solver/Git。v2 实现与测试 pin：

| 文件 | SHA-256 |
|---|---|
| `src/pgfrag/bond_visibility.py` | `3dfa2827a6e5b2e8589c0381e8729b2109092aaa7e30486777e4dfb927b21216` |
| `tests/test_bond_visibility.py` | `ef2820f5813e7b83e433eac443f4063317c828152c85790dd3b98b85eb6e2307` |
| `v2_diff.patch` | `3fd0b7231ee6f28d906f768062c1bec271d4f7cef22ea5f33dd076a8a6512850` |
| `v2_freeze_manifest.json` | `1b1cf3e7e1df4a8189266a2b4b78f61b0d0832abc332a7e715bbfff611a0e8b9` |
| 当前 `design_and_freeze.md` | `b4d22bbc7079f0ef75a9e501f4e35f7445bdf9439af356a7da7b8379378c39c9` |

以上 source/test/patch/design 值与冻结 manifest 一致。

## 复审结论

**v2 修复了 v1 的交点链假放行。** `_cluster_hits` 现在将任意不同参数数目的交点组标成 unresolved；主分类在读区间 occupancy 前因 `hit_reasons` 返回 `UNKNOWN`，所有 interval 标成 diagnostic-only，`admissible_static_candidate` 只由最终 `IN_MATERIAL` 决定（代码 498-577、897-965）。纯数学测试用模拟的 IN 中点覆盖单链和短空腔，并断言 `UNKNOWN` / 不 eligible（测试 266-339）；这直接覆盖旧缺陷，不依赖 CAD 样本。静态源码未发现 diagnostic midpoint 能绕过 unresolved 状态成为 bond。

数字防护也有实质改善：segment vector 先检查有限，再用 `math.hypot`；不有限、零长度、下溢长度都有独立 UNKNOWN reason（800-852）。合并精度带及 `2*band` 都必须有限，失败时 endpoint classifier 收到 `None` 并 UNKNOWN；interval/hit 中点和线段插值有有限检查，raw hit 点失败仍保留 raw 参数并使 pair UNKNOWN（73-120、426-495、606-728）。新增 numeric 测试用 `allow_nan=False` 编码结果并覆盖带宽溢出、端点差溢出及大 hit 参数（测试 341-437）。这足以支持一次有界的下一步运行检查，但不是物理或资产验收。

## 尚存问题

- **P1 接入门禁仍缺失。** v2 仍明确声明 `native graph is not connected` 和 `whole_asset_physics_status=NOT_ASSESSED`；对 `pgfrag/src/pgfrag` 的静态调用搜索没有 backend/assets 调用。没有根据关键区域/阈值 UNKNOWN 拒绝整个 asset，也没有证明筛选先于 neighbor counts、bond IDs、correction、kernel、storage。实际 vessel 仍不得据此准入求解。
- **P2 来源精度仍是 caller trust。** `_source_from_shape` 仅验数值有限非负、basis 非空；未读取 OCC 顶点/边/面容差或将其同申报值比较。STEP expected SHA 由 caller 提供，且候选坐标没有独立 asset/source hash 绑定。这个实现能核对“读到的字节是否等于 caller 给的摘要”，不能自行证明摘要、精度声明或点集来源权威。真实资产之前仍需可信 manifest 与几何 tolerance 上界；无证据就 UNKNOWN。
- **P2 极小浮点中点边缘。** `_safe_midpoint(low, high)` 用两个半值求和，虽可避免大同号坐标相加溢出，但最小正 IEEE-754 subnormal `x=5e-324` 满足 `x/2 == 0`，所以 `_safe_midpoint(x,x)` 会落到 0。若单个精确切点在此参数、classifier 将两端和整段诊断中点都判 IN，事件被移到 segment 起点；区间构造会丢弃 `t=0` cut，切线检测找不到 before interval，理论上可能让切线候选成为 IN。真实 OCC 内建 tolerance 可能会把这种几何判为 ON，但代码没有把该下界纳入精度门禁。现有极值测试只覆盖大数中点，不覆盖 subnormal/相邻浮点值。建议对 `low == high` 返回原值，并对没有可表示严格内部点的区间 fail UNKNOWN，或明确要求 tolerance 高于数值分辨率。
- **P2 无效精度 basis 可破坏严格 JSON。** 精度 basis 类型不正确时会加 UNKNOWN reason，但 `_source_from_shape` 仍把原值存进 `source_accuracy_basis`（239-307）；若调用者传 `float('nan')`，结果对象仍含 NaN，外部 `json.dumps(..., allow_nan=False)` 会失败。正常合法 contract 要求字符串，所以这不是 eligible 假放行；但“所有拒绝结果可有限 JSON 序列化”还未覆盖。无效 basis 应规范成 `None`，并增加一个纯数据回归。
- **切线/共面/共享边 OCC 行为仍未实测。** 当前结构在 `IsDone=True,NbPnt=0` 时把整段变成一个 midpoint interval；端点和中点均 IN 时允许 eligible，正确性依赖 OCC 全命中返回完整。测试中的球外切、ON 面用例仍未运行；没有端点均在材料内的切线、共面重合、共享 edge/vertex 的 CAD 控制。

## 最小运行建议

v2 的 pure-math regression 已经是有效的软件负控，能证明 diagnostic interval 即便全 IN 也不能把 unresolved hits 变成 eligible。建议 root 可放行设计中列出的**一次有界 pytest 命令**：pure-math regression 类，加 hollow-shell air-path、hollow-shell all-hits、solid-coupon 三项。它只能验证这条 OCP 调用路径和首批分段行为；不验证真实 STEP hash-match、precision 可信度、切点/重合/边界退化，也不接 native，因此不能标记真实 asset admitted 或物理 PASS。若要宣称有限范围完整，另补上 subnormal 中点保护测试。

## v1 design SHA 差异

v1 报告记下当时读取的 live design SHA 为 `1e943da11387344169d66350dbdb9740db18ed6b245f8e5643c318141884fe2`；v1 原档 manifest 中保存的 `design_and_freeze.md` 是 `f95c86eada924153ae0c0467e928d765d4d23750dc6abc5759215e5a3a3d227e`。同一 manifest 说明 v1 source/test 的 SHA 与原档一致，而 design 不一致；当前活动设计已是 v2 的 `b4d22...`。最稳妥的归因是：v1 审阅读取的是归档前的 live draft，design 在形成原档时已有不同字节。原 `1e943...` 字节没有留在 archive，故无法重建其逐字差异；本报告不改写 v1 归档，也不把 v1 报告说成审阅过 f95 版本。
