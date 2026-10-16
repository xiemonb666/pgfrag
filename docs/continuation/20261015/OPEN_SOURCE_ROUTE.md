# 开源物理路线：一个轴对称纹样器物的动态断裂数据

**研究日期：** 2026-10-15。**范围：** 只读核对项目状态与官方开源来源；未安装依赖、未运行求解器、未重跑旧 one-shot。本文工期是工程估计，不是测得的运行性能。

## 结论

若“自由裂纹”和大旋转客观性是硬要求，短期主路线应是项目已锁定的 Peridynamics.jl 0.5.4，使用 **OSBMaterial + CriticalStretch + VelocityVerlet**；不选 3D bond-based BB，因为它把泊松比固定在 0.25。OSB 的失效取决于参考/当前键长比，纯刚体旋转不改变伸长，力向量随旋转协变，因而这条路径具备原生大旋转客观性的结构基础。仍须做刚体旋转回归，不能把理论不变性当作本工程的数值验收。

如果“尽快取得可测面积的裂面/碎片数据”优先，且接受裂纹沿有限元面网格扩展，Akantu 5.0.7 cohesive 路径更直接：它有显式动力、内禀/外插 cohesive 面、面损伤/开口/牵引输出以及 FragmentManager 的质量、惯量、速度统计。用户环境说明有 WSL Akantu 5.0.7；仓库盘点同时看到 v5.0.7.post1 源码与 venv 材料，但本轮未核实导入/运行。它不是任意几何中的无网格自由裂纹。Peridynamics.jl 则保留自由断键，但键事件不会自动变成唯一的正面积裂面。

当前工程没有一条现成路线能直接产出“经实测校准的真实器物训练样本”。最快可信实现是有限扩展已有 PG-Frag/Peridynamics 路线，不加新的大依赖：器物体积材料点及正体积权重、内腔可见性键筛选、逐步断键记录、裂面/碎片重建和材料/网格验证都需要补齐。Akantu 可作为网格约束裂面基线，不是默认替代。

## 按此目标比较

| 路线 | 已有能力与适合点 | 关键缺口/限制 | 建议 |
|---|---|---|---|
| Peridynamics.jl `BBMaterial` | 最简单的中心键力；`CriticalStretch` 删除失效键；原生 `VelocityVerlet` 和预裂纹 API。键长度只依赖当前/参考距离，刚体旋转客观。官方 Kalthoff–Winkler 示例用 3D 动态断裂。 | 3D bond-based 弹性限制 `ν=0.25`；不匹配未必恰为 0.25 的陶瓷实测弹性。裂纹事件只是一组断键，非裂面几何。 | 仅用于 ν≈0.25 的材料对照或上游案例复现。 |
| Peridynamics.jl `OSBMaterial` | 保留中心键力，但力幅由局部加权体积/膨胀率状态决定；允许以 `E,ν` 表征各向同性线弹性并使用 `Gc` 或临界伸长；刚体旋转不改变键长、膨胀率和断键判断。项目已锁定 Julia 1.10 与该上游 physics source。 | `CriticalStretch` 是理想不可逆脆性拉伸删键；当前 PG-Frag 表面资产仍 fail-closed，孔腔穿越键需在原生 BondSystem 构造前剔除；生产 runner 没有逐步裂面/碎片输出或已证明的碎后自接触。一次历史 coupon run 的 6,810 个断键全在夹持边界，不能算内部裂纹通过。 | **满足“自由裂纹 + 大旋转”的首选**，先保持 OSB，不先换 correspondence。 |
| Peridynamics.jl `CRMaterial` / correspondence | 对应力模型提供 deformation-gradient 本构；`CRMaterial` 明确做大变形应力旋转，0.5.4 仅支持 `SaintVenantKirchhoff` 和 `LinearElastic`。 | 位移梯度/逆矩阵和零能模态稳定化增加了首例风险；大损伤下局部形状张量退化是额外门。并非最快自由碎裂入口。 | 除非之后需要与实验相符的有限变形本构且 OSB 不足，否则不放进首个器物例。 |
| Akantu 5.0.7 cohesive | 仓库里有 v5.0.7.post1 源码/venv 候选；3D 连续体与 cohesive facets；显式 lumped-mass 时间积分；外插插入按有效牵引 `σ_eff > σ_c` 触发，可限定允许区域；支持 cohesive 损伤、开口、牵引及耗散能读取，亦有破裂碎片的质量/惯量/速度统计。 | 裂纹仍沿被插入的网格面发展，路径受面方向与网格影响；3D vessel CAD→高质量体网格、cohesive 参数、接触及结果适配尚未接入 PG-Frag。`finite_deformation=true` 是大变形选项，不等于对当前材料与裂面算法已通过刚体转动客观性测试。 | 若接受预期裂纹带/mesh-face 路径、需较快的面积型裂纹基线，值得做一次固定面参考案；不替代 PD 自由裂纹。 |
| Peridigm | Sandia 的成熟开源 PD 代码，官方示例覆盖显式动态、LPS/OSB、impact/contact 和脆性破裂。 | 依赖 Trilinos/MPI/CMake 的独立构建和数据适配；当前已有 Peridynamics.jl 锁定运行时，换引擎不会解决几何/裂面/校准问题。 | 目前不迁移，也不安装大型依赖。 |

## 方程、薄壁与裂面解释

对任意材料点对，Peridynamics.jl `CriticalStretch` 使用

\[
\varepsilon_{ij}=\frac{l_{ij}-L_{ij}}{L_{ij}},\qquad
l_{ij}=|\mathbf x_j-\mathbf x_i|,\quad L_{ij}=|\mathbf X_j-\mathbf X_i|,
\]

并在允许失效且 `εij > εc` 时不可逆地删除键。因 `|Q r|=|r|`，刚体变换 `x'=Qx+c` 不改变该阈值；OSB 状态力也只由这些键长度、核权和当前键方向构造，故纯刚体旋转时 `l=L`、伸长/膨胀率为零，内力为零，叠加旋转下力按 `Q` 旋转。这是模型结构推论，不是器物算例结果。BB 的键力是 `f=c ε n`，3D 固定 `ν=1/4`；OSB 取 `K=E/[3(1−2ν)]`，保留独立的 `E,ν`。

均匀 OSB 在 0.5.4 源码中的核心量为（`ω` 是 kernel 权、`β` 是 surface factor）：

\[
m_i=\sum_{j\in H_i}\omega_{ij}\beta_{ij}L_{ij}^2V_j,\qquad
\theta_i=\frac{3}{m_i}\sum_{j\in H_i}\omega_{ij}\beta_{ij}L_{ij}(l_{ij}-L_{ij})V_j,
\]
\[
\mathbf p_{ij}=\omega_{ij}\beta_{ij}\frac{\theta_i(3K/m_i-5G/m_i)L_{ij}+(15G/m_i)(l_{ij}-L_{ij})}{l_{ij}}\,\Delta\mathbf x_{ij},\qquad
\mathbf b_i^{int}=\sum_j(\mathbf p_{ij}-\mathbf p_{ji})V_j.
\]
源码按有向键分别计算两点的状态力，并将每一项的相反体积加权贡献写回两端；故合力需合并 `ij` 与 `ji` 两个方向。超叠刚体旋转时所有 `l,L` 与内禀状态量保持不变，当前 bond 向量随 `Q` 变换；无变形参考态的内力仍为零。

3D `CriticalStretch` 在该固定实现中由 `Gc` 换算为

\[
\varepsilon_c=\sqrt{\frac{5G_c}{9K\delta}},
\]

因此 `Gc`、`K`、horizon `δ` 不是互不相关的独立输入；API 要求 `Gc` 与 `εc` 只给一个。这个关系用于该模型的阈值，不表示仿真中观测到的裂纹耗散恰好等于 `Gc A`，也不包含陶瓷的速率/缺陷统计。

### 薄壁的实际门

Peridynamics.jl 官方例子里 `δ/Δx=3.015` 和 `4.015` 是已展示的离散设置，不是器壁合格阈值。其 Kalthoff 样例厚度为 9 个点距、horizon 约 4 个点距；这只能作为首轮 3D 壁厚采样的参考。针对 VesselForge 实体，入场前须记录局部壁厚 `t`、点距 `Δx`、horizon `δ`、点层数 `t/Δx`、比值 `t/δ`、局部纹样波长/幅度与曲率，并在同一闭合实体上做点距/horizon 对照。若仅一两层点跨壁，或纹样与壁厚低于分辨尺度，应先停止，不把看起来连续的点云当作有效 3D 壳模型。

`EnergySurfaceCorrection` 是针对缺失 horizon 邻点的能量型自由表面修正：代码分别施加三个坐标方向小均匀伸长，取 `m_{d,i}=Ψ^{theory}_d/Ψ^{PD}_{d,i}`，再按 bond 方向组合邻点修正因子。它不是对薄壳、弯曲、孔腔、纹样或裂纹面的通用边界条件；代码初始化时按未变形键网计算。对 3D bond 方向 `n`，实现把各向因子端点平均为 `s̄_d=(s_{d,i}+s_{d,j})/2`，再取 `β_ij=[Σ_d n_d²/s̄_d²]^{-1/2}`。Le–Bobaru 原论文专门比较 PD 表面修正的弹性与断裂效果；PG 仍需用实际器壁厚度和点云做独立对照，不能凭启用此选项晋级。首例把 `NoCorrection` 留作诊断对照，把 `EnergySurfaceCorrection` 作为待验证候选；只有闭合实体键筛选完成后再建键、计算修正系数，并与解析/实测弹性壁响应对比。

尤其是 hollow vessel：上游 `Body(position, volume)` 按 Euclidean horizon 搜索生成候选键，不接收几何可见性策略。若只凭端点在材料内，键可能穿过器物内腔，把内外壁错误耦合。项目已有 CAD 可见性分类器和 native 初始键筛选的来源/设计稿，但它们仍是未接入/未验证路径。筛选必须发生在 native halo、键索引、表面修正、OSB 加权体积构造前；跨空腔几何排除不得记为物理断裂。

### 键图不是裂面

PG-Frag 的 [`bond_events.py`](../../../src/pgfrag/bond_events.py) 事件导出合同明确把连通分量标为 `intact_horizon_bond_graph_component_only`，`positive_area_crack_surface=UNKNOWN_NOT_DERIVED`。断键可以约束裂面的位置，却一般不唯一决定裂面三角形、面积、分支与双侧祖先；键图分成两个连通分量也不自动说明有可用面积裂纹或刚体碎片。Akantu 的 cohesive facets 自带参考面与面积，因而容易给出“面级标签”，但仍须检验它们确实完全脱粘、属于不同后继实体，并能与实物断口配准。任何训练用裂纹真值应与 solver 预测分开存。

## 最小验证与最快实际步骤

1. **材性先验。** 从同批/同条件实验给 `ρ,E,ν`；用动态拉伸/弯曲/预裂纹冲击数据约束强度与断裂能。OSB 首案固定一组 `E,ν,ρ,Gc`，由公式算 `εc`，不要用“让器物裂开”反向调参。无实测数据时只叫诊断数据，不能叫校准陶瓷结果。
2. **实现 sanity gate。** 用无损均匀点集通过力/能量为零、线弹性响应、总内力反对称/动量守恒、零损伤的 0°/90°/任意轴刚体旋转；再用空心圆环/短筒验证 cavity-bond 全部正确剔除、内外表面支持量、正体积与质量闭合。
3. **断裂与数值敏感性。** 先复现 Peridynamics.jl Kalthoff–Winkler 动态预裂纹基准（软件/时间步/裂纹事件链自检，不是陶瓷校准）；再测一组未用于拟合的陶瓷裂纹路径/载荷历史。快速筛查至少要 `Δt` 减半和两个点距、两个 `δ/Δx` 候选；若要宣称离散收敛，固定物理 `δ` 做三层 `Δx` 细化，并单独报告改变 `δ` 的敏感性。比较峰值力/裂纹起始时刻/路径距离/面积、动能-弹性能-边界功-断裂耗散和碎片质量/动量，不以碎片数作目标。
4. **训练资格。** 第一条 vessel run 固定 asset/STEP SHA、材料批号/参数来源、纹样进入几何或材料场的定义、支承/冲击、时间步、点云/正体积权重、horizon 与键策略、solver/driver hashes。预测裂纹/碎片标签与实测高速影像/扫描/传感器标签分字段记录；只有实测和复核数据才叫实验标签。缺实测样本时可留作 synthetic diagnostic，但不能冒充 experimental training data。

### 哪条最快

- **最快得到一个有意义的无网格断键诊断：** 继续复用当前 Peridynamics.jl OSB。它已有本地 Julia 1.10、固定 0.5.4 physics source、VelocityVerlet 和实验事件记录器；旧 4,500 点 coupon 的历史 launcher 为 96 步记录约 40 秒，但该案未校准且断裂集中在夹持边界，故只说明算力路径能动，不说明器物物理已通。若补齐 STEP/闭合实体到体积材料点的入口、并有实测载荷，工程估计数个工作日可做出第一个**诊断** case；真正裂面、校准及 3-level 数值门至少按 1–3 周计划。若材料实验和实测标签未准备好，资格时间未知。
- **最快得到带明确面面积的对照裂纹：** Akantu 外插 cohesive，适用于已选定裂纹带且允许 mesh-face 限制；且本地 venv 确认可导入后，若器物 tetra mesh 已准备，预计可在数个工作日内做出受限的 explicit diagnostic。若还需核实环境、CAD 体网格、接触/碎片导出和 PG adapter，优势会缩小。该估计不等同无网格裂纹。
- **当前建议：** 不换大依赖；先让 OSB 首例以“来源绑定诊断”通过物体/孔腔/薄壁键图及旋转门。若真实器物要求的界面面积训练标签是首要目标，再补一组 Akantu mesh-face baseline；不能用 cohesive 通过替代 PD 自由路径结论。完整 physical/training-ready 结果须等材料和实测标签门通过。

## 证据边界与来源

[根 handover（2026-10-07）](../../../../docs/handover/20261007/HANDOVER.md) 仍把物理合格包记为 0；[pgfrag README](../../../README.md) 明说生产 runner 当前只有有界均质 bar 片段、表面资产返回 `UNSUPPORTED`。[headless backend 说明](../../headless_backend.md) 说明 Julia fork 的 I/O patch 不改 20 个 physics/time-solver 源文件，但只验了来源身份和 import/run viability，不是上游回归、校准或物理验证。[历史 pilot summary](../../../experiments/free_fracture_pilot_20261005/run_attempt_01/results/summary.toml) 与 [launcher receipt](../../../experiments/free_fracture_pilot_20261005/run_attempt_01/launcher_receipt.json) 记录一例实跑：4,500 点、96 步、6,810 个断键全在夹持边界；不重跑、不晋级。

- [Peridynamics.jl v0.5.4 OSB 源码（官方标签）](https://github.com/kaipartmann/Peridynamics.jl/blob/v0.5.4/src/physics/ordinary_state_based.jl)；[BB 源码](https://github.com/kaipartmann/Peridynamics.jl/blob/v0.5.4/src/physics/bond_based.jl)；[临界伸长公式源码](https://github.com/kaipartmann/Peridynamics.jl/blob/v0.5.4/src/physics/fracture.jl)；[CR 应力旋转源码](https://github.com/kaipartmann/Peridynamics.jl/blob/v0.5.4/src/physics/correspondence_rotated.jl)。
- [Peridynamics.jl 官方动态 Kalthoff–Winkler 教程](https://kaipartmann.github.io/Peridynamics.jl/stable/generated/tutorial_kalthoff-winkler_dynfrac/)；[Peridynamics.jl JuliaCon 2024 论文](https://proceedings.juliacon.org/papers/10.21105/jcon.00165.pdf)。教程网页是当前 stable 文档；具体公式/API 以项目固定的 v0.5.4 源码为准。
- [Akantu v5.0.7-post1 官方 solid mechanics/cohesive 手册](https://akantu.readthedocs.io/en/v5.0.7-post1/manual/solidmechanicsmodel.html)；[官方 cohesive 示例](https://akantu.readthedocs.io/en/v5.0.7-post1/examples.html)；[Akantu 官方 I/O 文档](https://akantu.readthedocs.io/en/v5.0.7-post1/manual/io.html)；[Akantu JOSS 论文](https://joss.theoj.org/papers/10.21105/joss.05253)。
- [Le & Bobaru, *Surface corrections for peridynamic models in elasticity and fracture*, Computational Mechanics (2018)](https://doi.org/10.1007/s00466-017-1469-1)。
- [Sandia/Peridigm 官方仓库](https://github.com/peridigm/peridigm)；其 README 说明依赖 Trilinos，并以显式 impact/fracture、OSB、contact 为能力例。








