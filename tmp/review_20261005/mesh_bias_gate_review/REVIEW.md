# Native fracture 网格偏置门审阅

**日期：** 2026-10-05  
**范围：** 只读检查 `pgfrag/docs` 与 `pgfrag/experiments`；先用 `rg` 定位。未运行 solver，未改生产代码或现有实验。

## 结论与现有证据

目前**没有已执行且可接受的 native fracture 网格偏置验证**。现有结果证明了有限几何 QA 和窄 solver/API 能力；没有证明裂纹路径、断裂功或碎片证据对独立四面体重网格稳定。

- `experiments/wall_resolution/run_roi_trial.py` 对 VF 母体的一块杯壁 CAD ROI 做静态网格、实体内点、体积和壁厚取样。ROI 框为 x=[−23,−14]、y=[−1.4,1.4]、z=[23.1,25.9] mm。脚本明言未建 solver/bond graph/crack；h=0.5 mm 时的 t/h≈8 只是几何采样比，δ=1 mm 只是邻域诊断值。旋转原网格是 covariance 检查；独立 remesh 仍是待做事项。
- `docs/research/INNOVATION_LEDGER.md` 明记：无标定材料、无独立 remesh 方向结果、无 physical R2；h0.5 只有有限静态几何 QA，旧 h8/h6 whole-solid 网格被拒。方向、能量验收容差没有本地证据基础。
- PD loaded diagnostic 的窄 pass 没有 bond event。B 是对所有点施加仿射变形后调用一次力/失效评估，不是积分裂纹。原 5 mm 立方体取 h=1、δ=3 mm，每个点的 horizon 都碰到外表面。
- `free_evolution_coupon_20261005` 标记为 `SOURCE_ONLY_NOT_RUN`：解析盒 4,500 点、h=0.6、δ=1.8 mm、δ/h=3、prototype material 未校准、NoCorrection 且禁用失效。它是可逆响应设计，不是自由裂纹或网格收敛证据。
- cohesive 201 点样例是受控本构面、8 个界面单元、全自由度指定加载；其 manifest 声明不提供 mesh-bias、R2 或训练资格，run 也不完整。现有陶瓷基准没有过程区长度或 δ 映射。

## 最小可辨识 native 对照

样本选同一 VF 母体的 **3D 杯壁缺口试样**，以现有 cup-wall ROI 定位和测壁厚；ROI 裁切面不能作为裂面。应恢复足够母体几何，或让裁切面和夹持边远离裂纹评估区，距离由已标定的 δ/过程区长度决定。用 native grip velocity BC（或固定的实测冲击波形）加载，只约束物理夹持点，让内部自由演化；不做全点仿射加载。

开始 fracture 比较前须冻结材料卡、物理 δ（PD）或 law-derived 过程区长度（cohesive）、载荷历程和表面模型。若有强度/缺陷随机场，按共同 rest-material 坐标和固定 seed 生成一次，再插值到所有网格；不得按点 ID 重新抽样。缺少标定时只能报告数值敏感性，physical label 保持 UNKNOWN。

| 阶段 | Native case | 变量与测量 |
|---|---:|---|
| 网格/邻域预检，不跑 solver | 3 个 placement × 每个 h | baseline、mesher-frame rotation、半格平移；映回同一 CAD/material frame。验证点集与 tet 方向分布确实改变。测局部 h、壁厚 t、t/h、δ/h、δ/t（cohesive 则 h/ℓₚ、t/ℓₚ）、邻域体积/方向张量、壁厚层数，并逐 bond 检查连线是否穿过 CAD 外部；区分真实自由表面和人工裁切面。 |
| 完整表面响应预检，无失效 | 2 | 同一网格/载荷比较 NoCorrection 和可用的 energy surface correction；比较表面应力/应变、刚度、反力、能量，选定后在裂纹比较中固定。surface correction 与跨空腔 bond clipping 分开审。 |
| 固定物理 δ 的空间/方向矩阵 | 9 个 fracture case | δ 固定，h 按 δ/h=3、4、5 细化；每级跑 baseline、旋转 remesh、平移 remesh。统一 CAD、缺口、rest-field、载荷时长、材料、失效律、修正和 extractor；每次重建 bond graph。 |
| 时间步检查 | +1 | 最细 h 的 baseline 单独将 Δt 减半，保持物理载荷时程不变；不与 h 收敛混算。 |
| δ/过程区敏感性 | +1 | 最细 h、baseline mesh 固定，只换一个由独立标定不确定度支持的第二尺度；若没有合理区间，记 UNKNOWN。 |

合计最小 **13 个 native cases**，另加网格/邻域预检。δ/h 的 3、4、5 是收敛序列候选，不是通用合格标准。PD-centroid 离散应由实际点体积/邻距给出 h，不能只报 Gmsh 目标边长。旋转原网格并同步旋转载荷只检验 covariance；变换 mesher frame 后必须 inverse-map 回同一物理坐标，且检查 mesh/point set 实际不同。始终报告共同 rest frame 下的力—位移、起裂负荷/时间/位置、裂纹路径/表面面积及法向不确定度、有效断裂功、能量/动量账、组件质量；不设目标碎片数。PD bond failures 对裂面只是约束，不能自动当唯一、有面积的 CrackSurface。

## 预定门槛与预算

- **FAIL（几何）：** 裂纹评估区任一处 t/h≤1（壁厚不及一个实际采样间距）；有效 bond 穿过实体外部/空腔；或裁切/夹持边驱动裂纹。更大的 t/h 也不直接通过，仍须看 h 序列。
- **FAIL（数值）：** root 在确认试验前冻结的方向/能量/路径容差被超过且 h 细化不收敛；或半 Δt 结果超过独立时间容差；或能量、动量、质量账超出独立预算。方向重网格后裂纹在母体坐标改变并持续贴合 tet 方向，是直接偏置证据。裂面不唯一则 CrackSurface 门失败。
- **UNKNOWN：** 未标定材料/载荷/δ 或过程区；未冻结容差；修正、horizon 支持/穿空腔 bond 未核实；三组 mesh 实际没有变化；或 run/output 不完整。t/h≈8、solver exit 0、静态 CAD QA 都不能单独使门通过。
- **数值通过不等于物理或 A12 通过。** 要进入 labels 仍需事件支持的成对裂面证据、材料/载荷校准及后续 A12 degradation/full-label 检查；拓扑歧义保留 UNKNOWN。

现有文件没有可依据的方向、能量或拓扑百分比阈值；PG-RC-002 要求由 root 在 pilot 后、确认集前冻结。此处不造一个“5%”之类的门槛。

资源按一个试样、一个冻结 field realization、单线程规划。先预检全部网格和点数，再决定是否执行；已有 reversible coupon 源设计上限为每例 5,000 点/192 步，约 150 MiB/193 个全历史状态。若 13 例都能落在该预算内，全历史粗估约 1.9 GiB；192 步可能不足以观察完整 fracture，故需围绕起裂保留全速 bond/point history 并记录到事件完成。超过预算或没观察到完整事件就停为 UNKNOWN，不得暗中加粗 h/缩短物理载荷。wall-resolution 的 10,000-tet 上限是几何试验上限，不是 native solver 预算；现有材料也不足以估运行时间。

**出处：** `experiments/wall_resolution/run_roi_trial.py`；`docs/research/INNOVATION_LEDGER.md` PG-RC-002；`experiments/native_resume_20261005/next_free_evolution_coupon_source_plan.md`；`experiments/free_evolution_coupon_20261005/SOURCE_REVIEW.md`；`experiments/pd_loaded_diagnostics_20261004/manifest.json`；`experiments/cohesive_bilinear_patch_20261004/manifest_bilinear_modeI_201_run03.json`；`docs/references/CERAMIC_FRACTURE_BENCHMARK_20261004.md`。

