# PD loaded diagnostics：决策与A结果 v4（2026-10-04）

**修订状态。** 本稿纠正 v1 的稳定步估算错误，并记录主审授权后唯一一次 A 加载 observer run 的结果。v1 原文、v2 纠错前稿和 v3 运行前设计稿分别保存在 `decision_v1.md`、`decision_v2_pre_wavecorrection.md`、`decision_v3_pre_a_run.md`；修订轨迹见 `decision_history.diff.md`。v1 把角点局部步长误当全局步长，导致 320 步预测高估约 2 倍。A 使用冻结且经审阅的 runner/case/common/launcher 哈希单次执行并通过窄响应 gate；B 仍为未运行设计。

**范围与已知证据。** 审阅了 PG-Frag frozen point-probe、run01 原始 summary/receipt 和 pinned Peridynamics 源；按授权只执行了一次新的 8 步 A loaded observer run，未运行 B、CAD 或 Git，未更改冻结 runner、输入或旧 case。零载 `zero_load_control_run01` 已由主审验收：125 点、1 thread/1 chunk，native accepted dt `4.925124976599115e-7 s`，6,124 directed / 3,062 canonical bonds，2 个采样态、零 transitions，数据/质量/互反状态 gate PASS；总能量账本 UNKNOWN。`common_case.toml` 仍是 `max_steps=5`、timeout 120 s；tinyload 是 5 步、±`1e-4 m/s` 的准备配置，未运行。其实际 dt 下名义应变仅 `1.2313e-7`，约 `0.00764% εc`，零事件不能验证 fracture。

点集为 5 mm 立方体的 cell-center 离散，坐标中心面 `z=0.5…4.5 mm`，加载面 gauge 为 4 mm；horizon `δ=3 mm`，全域强 surface-effect。材料 `E=5 GPa, ν=.24, Gc=45 J/m², ρ=2100 kg/m³` 都是 `UNVALIDATED_PROTOTYPE`，不是陶瓷实测/标定。native 源给出 `K=E/[3(1−2ν)]`、`εc=√(5Gc/(9Kδ))`，故 toy `εc≈1.61245155e-3`。不可由此推称真实断裂资格。

## 推荐的两个独立小 gate

A 与 B 独立：A 复用 frozen runner 的原始 SHA，B 是全新隔离的 native constitutive branch diagnostic；不改旧 runner、common、tinyload 或 run01 证据。A 配置用独立 `max_steps=8` common/case、同一 frozen point CSV/metadata、SI、OSB/CriticalStretch、`AUTO_NATIVE_STABLE_BOUND`、thread=1、单 chunk、±`0.005 m/s` 的 z 向 top/bottom 恒速 BC。沿用 run01 实测 dt，8 步 `T=3.940099981279292 μs`、名义 gauge 应变 `9.85024995319823e-6≈0.6109% εc`。按连续体参考公式 cL=√((K+4G/3)/ρ)≈1675 m/s，4 mm gauge 波时约 2.39 μs；8 步观察窗约 3.94 μs，即约 1.65 个 gauge transit。native VelocityVerlet 对全部 point velocity 与 velocity_half 字段默认置零，loaded runner 分支未设初始速度；在首个积分步，边界条件再将 top/bottom 半步速度直接设为 ±v，所以载荷是骤启显式动力响应，未做 hold/平衡检查，不能称准静态，也不测能量闭合或自由断裂。A 单次实际运行与原始数据核验见下节。A source pins 与 argv/exit 记录在实验 `manifest.json` 和 launch receipt；启动及运行约 8.51 s，timeout 冻结为 120 s，没有重试。

### A loaded-response observer：实际窄 gate 通过

单次 run `loaded_elastic_a8_run01` 于 2026-10-04 完成：wrapper exit 0、Julia child exit 0、未 timeout。实际 native dt 与 run01 一致，为 `4.925124976599115e-7 s`；总时长 `3.940099981279292e-6 s`。原始输出有 1,125 point-state rows、55,116 directed bond-state rows 和 27,558 canonical bond-state rows（共 9 个采样态）。每步 125 点、6,124 directed / 3,062 canonical 键；每步 mass sum 为 `0.00026249999999999966 kg`。坐标满足 `current = reference + U`，最大分量残差 `1.608322006572579e-18 m`；互反及 canonical activity 每步一致。9 态中 U、V 和内部力密度均出现非零值，两个 transition CSV 均只有表头，未观测到断键。只读逐步极值记录在 `pd_loaded_diagnostics/reports/loaded_elastic_a8_readonly_audit.json`；以相同原始 CSV 绘制的科学图为 `pd_loaded_diagnostics/reports/loaded_elastic_a8_science.png`（SHA-256 `b3df2be1cff05a4060d4b81112686287f9d925aa8d5511bb5af9ad9ac317764e`），几何位移倍率为 1×，U_z/位移单位和内部力密度极值标在图中。运行原始数据与 receipt 在 `pd_loaded_diagnostics/snapshots/loaded_elastic_a8_run01/` 和 `pd_loaded_diagnostics/receipts/`。

| step | max displacement norm (m) | max velocity norm (m/s) | max internal body-force density norm (N/m³) | active canonical bonds |
|---:|---:|---:|---:|---:|
| 0 | 0 | 0 | 0 | 3062 |
| 1 | 2.46256e-9 | 4.02814e-3 | 1.15057e7 | 3062 |
| 2 | 4.95156e-9 | 3.53930e-3 | 2.01164e7 | 3062 |
| 3 | 7.57403e-9 | 3.87061e-3 | 2.63545e7 | 3062 |
| 4 | 1.03513e-8 | 4.33615e-3 | 3.14637e7 | 3062 |
| 5 | 1.31454e-8 | 4.39778e-3 | 3.72705e7 | 3062 |
| 6 | 1.58362e-8 | 3.94557e-3 | 4.50659e7 | 3062 |
| 7 | 1.84288e-8 | 2.51464e-3 | 5.43165e7 | 3062 |
| 8 | 2.09902e-8 | 2.81092e-3 | 6.35960e7 | 3062 |

这只通过多步 observer/加载响应链路 gate。时间窗短且为骤启动态响应，所有材料参数仍是原型值；结果不验证弹性模量、准静态平衡、能量闭合、断裂阈值、自由裂纹、接触、粗糙断面、材料真实性或 RelicWeave R2 血缘。

B 不是数百步恒速加载。用全点仿射单轴映射 `F=diag(1,1,1+η)`，固定 `η∈{0, 0.5εc, 1.5εc}`，分别从同一 pristine native 状态构造每个状态，逐状态只调用一次 native force/failure evaluation，不进入 VelocityVerlet 时间积分。`η=0` 和 `.5εc` 是无断键控，`1.5εc` 是正控：按 native 几何必有伸长超过阈值的 z 向邻键。输出每个 directed bond 的 source/neighbor point id、参考长 `L`、native 当前长 `l`、`(l−L)/L`、point εc、fail-permit、调用前后 active，以及 canonical pair 汇总；不得手改 `bond_active`。检查低于阈值无 transition、超阈值确有 native active→inactive、互反方向与 canonical pair 一致。通过只说明 pinned CriticalStretch 失效分支响应真实仿射几何；全点同一仿射运动不是裂纹、自由动力学、材料标定、断面或 fragment 证据。

静态源审阅显示该 B 路径有原生调用依据：`threads_data_handler` 构造一个单 chunk，OSB `calc_force_density!` 对每点依次执行 `calc_failure!`、`calc_damage!`、`force_density_point!`；failure 使用 native `storage.position` 对原生 bond reference length 算伸长。Body 默认允许失效。新 driver 可复制 pristine state 的当前位置后应用 F，再调用 native handler 的 force evaluation 并比较真实前后键数组。现有冻结 runner 只支持 `zero_load`/top-bottom symmetric velocity，不能直接表达 B；`displacement_bc!` 源文档说明只用于 NewtonKrylov，不能冒充 VelocityVerlet 支持。当前证据仅为源码静态可行性，B driver 尚未执行；root曾读取早期草稿哈希前缀 c1855a8c…，但精确字节未保留、不重造。现保留的实际修前字节是 d7c8c11fe9d24ff3bd395240942d625bae0e17e6cbf70619b2d87aade09f03aa，并仅针对作用域读值及 volume finite 检查做了修订。

## 止损与物理边界

A 的设计 gate 已通过，且实际未发现断键；不因此把 A 解释成断裂 gate。B 尚未运行。B 在 native call 不能稳定构造、没有阈值下/上正确区分、pair 互反不一致或 provenance 不全时记 BLOCKED/FAIL；不手动改键、不改 E/Gc/δ 找成功。两 gate 都不是材料真实性验收。

native OSB 每次 force evaluation 都用当前 `bond_active` 重算 weighted volume `m_i`，再用于 dilatation/force/strain energy；破键会改变后续刚度/能量定义。故 B 的力值和断键后能量不作通过标准。现有 observer 缺少 actuator reaction、加载功与完整能量账本；energy、contact、真实 rough fracture surface、不同网格方向偏差、跨分辨率及陶瓷材料参数标定仍为 UNKNOWN/未做。不能把 toy PASS 升级为真实器物物理或 RelicWeave R2 真邻训练证据。
