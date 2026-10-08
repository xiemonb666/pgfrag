# 93-DOF 前沿有限扰动诊断设计（root决定版）

状态：SOURCE_ONLY_PREPARATION。保留上一决策版于 DESIGN.pre_root_decision_20261007.md、PLAN.pre_root_decision_20261007.json。本轮只修改本目录文档并新增短 source contract；没有运行 solver、provider、NPZ、数值源或测试，没有改已有 case。

## root 已定架构

阶段 1 是单一 93-DOF toy Mode-I 纯 LEFM 诊断：旧前沿 a=0 与有限扰动 candidate a=0.5 mm 沿 +x，在完全相同的 old-tip 外部 traction、相同物理刚体 moment 定义及零目标 gauge 下，各自重新求平衡。裂面两侧 traction-free；不装配 current virgin CZM，不加 cohesive stored energy，不引入材料 Gc，不作 fracture threshold，不提交 ExplicitSurfaceState，不宣称增长。结果仅为 FINITE_PERTURBATION_DIAGNOSTIC_ONLY。

这不是之前建议的 162-DOF Mode-III/384点方案；该方案与当前 mixed93 operator/boundary load 不匹配，禁止复用。也不能拿已有 mixed93 manufactured q_solution 当旧态 baseline：它的 physical gauge target 是从 manufactured 解生成的非零 g。旧、新态都必须分别在 g=0 下重解。

阶段 2 才进入真实裂纹增长：显式 cohesive history 随平衡路径渐进损伤，ledger 使用 stored-energy change 与累计/增量 cohesive dissipation ΔDcoh；cohesive law 的 Gc 参数用于本构，不再另减同一新裂面的 GcΔA。阶段 1 的有限 ΔΠ/面积数字不能替代阶段 2 的损伤、耗散、校准或 commit 门。

## 固定数据与物理合同

载荷来源是冻结 mixed93 pure-I 边界 traction 数组：缓存参考 KI=100 Pa√m，按 root/fcc4 的比例 0.01 缩放到 KI=1 Pa√m；它仍是 toy old-tip manufactured traction，不是已校准实物载荷。保持每个 boundary quadrature row 的 xyz、单位法向、单侧面积权重和缩放后的 traction_Pa 不变。对旧/新 basis 分别以这些同一物理样本积分广义外力，故 f_old 与 f_new 可以不同；不得把系数向量相等误当成物理载荷不变，也不得从 candidate 解重新计算 traction 或调整 KI。

冻结物理边界顺序和切片为 side [0,256)、cap_bottom [256,512)、cap_top [512,768)。数量为 side 256、两 cap 共 512、boundary 共 768；来源是 mixed93 a0/theta32 cache。体积积分 512 点；另 128 off-grid 点只用于场诊断。source grid 为 27 个固定 xyz、93 semantic coefficients 为 81 BASE + 12 GLOBAL_FRONT；域 R=L=0.01 m，z∈[-0.005,0.005] m，裂面 y=0、旧 tip x=0、front 平行 z。

使用固定的 a0/theta32 几何采样做同网格有限扰动对照：基准和 candidate 都在相同 512 volume + 768 boundary 坐标与权重上重新计算各自 branch values/gradients、K、physical-moment C 和广义 f。若 source IDs、MLS 配置、query coordinates、row IDs、weights、normals 或 cache hash 任一不精确匹配，则不得声称零重查询；本 scope 直接阻断，不暗中补 query。通过精确缓存映射时，BASE MLS values/gradients 可复用，provider 新 query=0；candidate branch 项在所有积分行解析重算。固定点规则只形成有限积分诊断，不证明候选 tip 周围的积分收敛。

每态求解约束为 K(a)q + C(a)^T lambda=f(a)、C(a)q=0。C(a) 在改变 basis 后要重新表达同一组 volume-average 三平移和归一化体积 rotation moments，中心与尺度保持固定；zero target 避免 manufactured boundary-motion gauge work。至少两次 99×99 KKT 解（old、candidate），每态一次 93×93 bulk assembly；cohesive primitive calls=0；不重试、不做多步循环。旧 manufactured nonzero-g cache 不算其中一次。

各态报告 U=0.5qᵀKq、W=fᵀq（同时核对同一外部 traction 下的边界场积分）、Π=U−W、g·lambda=0、平衡/约束/谱残差。可报告 ΔΠ 和 −ΔΠ/(L·δa) 作为有限差分数字，但不能叫经验证 J、不能比较 Gc、不能驱动状态提交。δa=0.5 mm 时单面几何面积差为 5e-6 m²；双侧 displacement traces 不把裂面面积乘二。

## 证据与边界

具体采样、单位、traction 缩放、重查询预算在同目录 SOURCE_CONTRACT.md / SOURCE_CONTRACT.json。上游 mixed93 case SHA-256 为 67ff4df51a0269e88b8d107a2a95929c3990a7cd9605585a31a2df583434cc1c；fcc4 case SHA-256 为 fcc4b46ed101b0f6067a9f3ec10816909c9fa931edd5a8cdde4c7f46cb3f3a71；root 授权缓存数组 SHA-256 元数据为 3e1e6fd17751e4e06113010f20ea33c7b69914ca90431e26bbf76124951b5d41。本轮仅读这些 case/source/receipt 的 JSON/Markdown；没有打开 NPZ。mixed93 实际 Cohesive bridge 仍为 SOURCE_READY_EXECUTION_DISABLED，且当前 fixed bridge 为 virgin elastic toy-only，不作为本阶段 operator 或能量输入。

generic cohesive assembly 的历史接口可保留给阶段 2：trial 是 immutable/read-only，汇总 stored/dissipated/work；压缩、切向/混合模式 UNKNOWN；softening tangent 可能 indefinite、拐点 nondifferentiable；没有全局非线性求解器或稳定支路跟踪。若实际历史耗散出现，阶段 1 立即不再适用，必须走阶段 2 cohesive-history ledger。snapback、非凸状态未证稳定或求解不收敛一律 UNKNOWN，不强推。

GLOBAL_FRONT12 仍只描述单条直前沿的 branch jump。多裂纹、多 front、转向、分叉、起裂、碎片数量推断、dynamic crack speed 均未实现。有限扰动 candidate 永不改写 active surface；此阶段没有 commit 状态。

## 合作交付

full3d operator agent 应据 SOURCE_CONTRACT 实现其另有目录所有权的数值源和 receipts，本目录仅提供冻结载荷、query、gauge、求解和 nonclaim 合同。实施入口必须区分 physical traction data 与 basis-dependent generalized force vectors，并分别把 old/candidate zero-gauge reanalysis 绑定到 geometry/basis/load/source hashes。
