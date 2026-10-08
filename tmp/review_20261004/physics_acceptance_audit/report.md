# Toy cohesive law 到真实器物自由断裂：五项物理证据缺口

本审计只读 PG 源码、冻结 manifest、留存的 case 结果和现有 research ledger；不运行 solver，不裁决后续架构路线。未读取或引用动态 run03 结果。**当前物理 QA = FAILED/PENDING；没有证据支持 physical R2 或物理预训练准入。**

## 已知证据边界

- Akantu extrinsic 96-tet toy 在预指定中央面动态插入 8 个 cohesive facets，能证明固定合同面的 native 插入/成对连通性；它不测试自由裂纹。mode-I 两例的原生 `dissipated` 是实质负值（约 `-5.57e-5 J` 与 `-9.73e-6 J`），超过运行前冻结的数值误差容忍范围，故为**已有 FAIL**；`Gc·A` 能量目标未达到。保留有符号原值，不用绝对值修正。run02 部分结果中的约 `-1.4e-23 J` 属舍入量级，须保留 raw 值但本身不构成 FAIL；该 case 仍因未完成而是 `NOT_EVALUATED_INCOMPLETE`。详见 `relicweave-synlab/docs/a12/COHESIVE_CAPABILITY_20261004.md` 及对应 frozen manifest、status 和 trace。
- intrinsic bilinear run02 是全 DOF 规定纯 Mode-I 的本构小片，外部 300.046 s 超时，只留下 168/201 个闭合样本，索引 0–167；没有完整终态，final `GcA` 与候选本构 gate 为 **NOT_EVALUATED_INCOMPLETE**，不是通过也不是完整失效失败。它不验证自由裂纹路径。细节在 `pgfrag/experiments/cohesive_bilinear_patch_20261004/manifest_bilinear_modeI_201_v2.json` 和 `relicweave-synlab/tmp/review_20261004/implementation/cohesive_bilinear_patch_20261004/controlled_intrinsic_bilinear_modeI_201_run02/`。
- 现有 h0.5 结果仅为有限静态几何 QA；没有被接受的 canonical cup-wall 动态自由裂纹独立重网格比较，也没有真实 ROI 的 `t/h` 或 cohesive `h/l_process-zone` 验证。CONTINUATION 中 Mesh32/Mesh64/Boolean 是 rough geometry / wall-mesh 诊断，不能证明裂纹方向偏差或收敛。
- Cohesive case 的 contact `NOT_RUN`。执行器请求 `getEnergy("contact")`，而所审 native cohesive energy key 是 `"cohesive contact"`；因此旧 trace 中 `contact=0` 不能解释为接触能为零，相关 ledger 可能漏项。真实陶瓷材料参数、粗糙断面的盲测也均未完成。

## 五项最小可证伪验收

| 证据缺口与当前状态 | 最小验收设计及可证伪判据 | 不通过时的路线条件 |
|---|---|---|
| **1. 网格锁定 / 自由路径。** 已有中央指定面插入；自由起裂路径、独立 remesh/orientation 偏差均 **NOT_RUN**。 | 先冻结几何 QA 通过的 canonical cup-wall 小 ROI、预缺口、载荷/边界及同一 rest-space 材料场。用至少 3 组独立网格方向或 seed，并在每组做至少一档 h 加密；固定 law 和材料，只另做独立 dt 检查。对齐材料坐标，比较起裂载荷/时间、力–位移、断裂功、路径/分叉拓扑、裂面面积/法向与质量。root 在跑前冻结各指标容差；独立方向差或粗细差超界即否证无网格偏差，通过同网格旋转不能替代此试验。没有预冻结容差则只能 `NOT_EVALUATED`。 | 本轮止损建议最多安排一次有依据的分辨率修正；仍超界则登记 `NumericPhysical FAIL`，停止 cohesive 物理晋级。后续是否评估有限 PD 共同裂面重建由 root 决定，且 PD 必须单独通过裂面歧义/R2 门；不得退回任意共享面标签。 |
| **2. 薄壁厚度与过程区分辨率。** 已有几何网格不等于动态物理分辨率；实测壁厚尺度和过程区收敛 **NOT_RUN**。 | 对 canonical ROI 实测局部壁厚 `t` 和局部单元尺度 `h`，报告分布而非单点 `t/h`；从已校准 cohesive law 求/测特征过程区 `l_p`，报告 `h/l_p` 与 `t/l_p`。在固定物理参数下同时加密 h，并比较起裂、耗散、路径和厚度方向应力/损伤场；若评估 PD，另报告 `δ/h`、`δ/t`，不与 CZM 收敛合并。先冻结所需分辨率及输出容差，再验证细化结果落入界限。 | 过程区/厚度未解析或指标随 h 显著变化时为 `NumericPhysical FAIL/NOT_EVALUATED`；本轮止损建议只安排一次有依据的细化诊断，仍不满足就停止物理标签晋级并交 root 决定是否转用另一种已冻结表示。不得用静态几何 QA、全局 CFL 或更小 dt 代替空间解析。 |
| **3. 接触与能量账本。** Extrinsic 实质负耗散为**已有 FAIL**；完整 GcA 未达到。 Contact 是 **NOT_RUN**，且 `contact` key 未验证，历史 contact 能量值应为 **UNKNOWN**， ledger 可能漏项。Intrinsic run02 因截断而能量门 **NOT_EVALUATED_INCOMPLETE**。 | 在文档确认的正确 cohesive-contact key 上做独立 key/preflight 断言；执行已知分离→闭合/接触小例和 no-contact 负控，导出界面 traction/opening/damage、接触力/穿透、冲量/动量，以及反力外功、动能、bulk/cohesive 可逆能、损耗、接触/摩擦/阻尼能的逐步账本。闭合后的穿透、能量残差和 GcA 误差需在预注册容差内；超过运行前冻结数值误差容忍范围的实质负耗散或缺账项可否证物理账本，舍入量级负值应保留 raw 但不单独判失败；不把 unknown selector 的 0 当物理结果。 | 先修正/验证测量接口后，若原生损耗仍超过冻结误差容忍范围为负或账本无法闭合，当前 cohesive 物理路线保持拒收，交 root 作路线决定；contact `NOT_RUN` 不可由 cohesive 接触惩罚接口替代为 PASS。旧 trace/status 保持不改，后续任何修正必须有新 manifest 与新 case。 |
| **4. 粗糙断面形貌。** Toy law 与固定平面界面没有真实粗糙断口对照；Morphology/roughness 是 **NOT_VALIDATED**。Mesh32/Mesh64 Boolean 是几何/磨损 QA，不是裂纹方向或断口形貌验证。 | 对目标材料制作独立真实断裂试样/器物 ROI 的高分辨率 3D 扫描；预先划分校准与盲测样本。以相同实测尺度和材料坐标，比较预测/测量表面的多尺度高度谱或相关长度、RMS 坡度/粗糙度、法向各向异性、面积与分叉/边界统计，并用扫描不确定度冻结容差。盲测统计超界，或只靠求解后加噪/平滑/切分才接近，都否证“物理粗糙度预测”。 | 若粗糙度缺乏可验证微结构机制或盲测失败，模拟断面只能标为 synthetic/unvalidated；转为采集/重建观察到的表面作为独立证据类型，或等待经 root 批准的材料形貌模型验证，不把美化结果当物理 R2。 |
| **5. 分材质本构标定。** 当前 `E, ν, ρ, σc, Gc, β, δc` 均是 toy 值；无材料来源/不确定度/独立实验验证，属于 **UNVALIDATED_TOY / NOT_RUN**。 | 对每个目标陶瓷配方/烧成批次建立来源可追溯、含不确定度的 `E, ν, ρ`、强度、I/II/混合模态断裂能、率效应、各向异性及空间相关尺度；用一组试样标定，再以独立缺口纯 I、剪切/混合模态及目标冲击率试样盲测力–位移、失效载荷、断裂功与形貌。预测误差需在事先冻结的实验和参数不确定度范围内。 | 未标定或盲测不符时，该材料及相邻未验证材质均保持 `CalibratedMaterial=UNKNOWN/FAIL`，不得材料间迁移、不得放入 physical R2 roster；只可限定到已通过的具体批次/条件，或用实测断口证据。通过 toy law 或单一 mode-I 曲线不代表陶瓷资格。 |

## 结论

已失败的是 extrinsic toy 的 signed cohesive dissipation gate；尚未评价的是其完整脱粘、动态 contact、自由路径方向偏差、薄壁/过程区收敛、粗糙断口盲测及真实材料校准。Intrinsic run02 因运行不完整只可作部分本构观察。五项任一为 FAIL、UNKNOWN 或 NOT_RUN 时都不产出 physical R2；mesh/build/API 成功或 Boolean/表面视觉 QA 不能替代物理验收。路线裁决留给 root。

文中的“本轮最多/只安排一次修正”仅是本轮资源与止损 advisory，不构成永久方法约束；root 可根据后续新证据重新冻结路线和预算。

本审计来源文件的 SHA-256 和字节数见同目录 `source_hashes.json`。
