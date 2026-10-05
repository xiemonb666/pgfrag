# Finite-front shadow 决策（2026-10-05）

**结论：当前 sharp-bool visibility + 分区内 MLS 梯度 + 纯体积分 K，不是已成立的 conforming 裂纹 Galerkin 方法。首选连续 MLS 基底上的 partition-of-unity 裂面/前沿增广；保留显式裂面路线 B，但停止把 bool-mask 小修补当作合法耦合路径。** 本次仅只读审查和解析推导，未运行 solver、CAD、pytest，未改变实现。当前仍为 `GEOMETRIC_PRIMITIVE / NOT_COUPLED`，不能产生 R2 或 GT。

## 1. 一个不依赖舍入误差的反例

物理裂面为有限三角片，顶点 `(0,-1,0),(0,1,0),(-2,0,0)`；取源点 `s=(-1,0,-1)`、查询 `q(t)=(1+t,0,1)`。连线穿过 `z=0` 于 `(t/2,0,0)`：小负 t 遮挡，小正 t 可见。切换点 `q0=(1,0,1)` 明显不在裂面上。

另设四个一直可见的源点 `q0+a(1,1,1), q0+a(1,-1,-1), q0+a(-1,1,-1), q0+a(-1,-1,1)`，`a=1/4`，全部位于 `z>0`；共同支持半径 `R=4`。以 q0 为线性基底原点，这四点在极限处同权 `w0>0`，其 moment matrix 为 `A0=4w0 diag(1,a²,a²,a²)`，两侧均满秩。给 s 的系数 1，其余 0。令 `ws=w(||q0-s||/R)>0`，Sherman–Morrison 直接给出：

```text
u(q0−) = 0,
u(q0+) = φs(q0+) = ws / [4w0 + ws(1+8/a²)] > 0.
```

这是现有 Wendland C2 也适用的解析反例，不是已执行测试。非零迹差在邻域中持续，形成 Γ 外的二维 shadow 面。其分布梯度含 `[[u]] ⊗ n δS`，故一般不属于 `H¹(Ω\Γ)`；只积分分片经典梯度不能代表原弱式。按分片积分部件，真实光滑应力产生的内部项 `∫S (σn)·[[v]]` 一般不消失。结果可能虚增柔度或失去一致性；因离散跳跃并不完全独立，不能夸称等价于逐面新增自由裂纹。

逐 quad 的 exact Fraction 只证明各开区间导数。把切换点标 UNKNOWN、避开它积分、提高 quadrature 或调小 offset 均不补这个空间缺陷。仿射系数在所有满秩 mask 下重现同一场，因此 affine patch 恰好会漏检。无裂纹 K 小测仍有独立价值。

## 2. 两个选项，首选增广

**首选：连续、无可见性删源的 MLS 基底 φ，加受支持域拓扑限制的 Heaviside 和前沿分支函数：**

```text
uh = Σi φi di + Σj∈J φj Hj aj + Σk∈T Σα φk Fα bkα.
F = √r {sin(θ/2), cos(θ/2), sin(θ/2)sinθ, cos(θ/2)sinθ}.
```

J 只包含支持域被真实 Γ 完全切成两部分、且不含前沿的基函数；T 覆盖前沿。严禁给前沿支持套无限平面 H。先用同一个解析直前沿坐标系，`θ∈(-π,π)` 的支割仅在真实裂面上。`∇(φF)=F∇φ+φ∇F`；`|∇F|=O(r^-1/2)` 在前沿法向截面可积平方，因此这个构造在 Γ 外有 H¹ 合法性。不得按“最近三角片/线段”随意切换局部坐标制造新跳面。曲前沿需要连续坐标/图册及一致支割，另行审查。

这借鉴已发表 PU/XFEM 构造，应用到当前 MLS 是明确的工程方法变更，不能继续称 Ai–Augarde 2018 的精确复现。代价是额外 DOF、增广相关性和裂面/奇异积分；收益是合法性由基函数结构给出。H 与 tip 集合按完整支持域确定；还需检查过渡区逼近误差与线性相关，不能靠正则化掩盖额外零模。

**备选：真正连续的 diffraction support。** 连续绕射距离及其导数可消除突删源，较少新增 DOF；但必须明确 3D 有限/曲前沿路径、路径切换处连续性、裂面两侧迹和满秩性。不能只对 bool mask 加 sigmoid，也不能仅在原错误基底上加 tip enrichment。本次未获得可直接复现的完整 3D 绕射规格，故不优先启动；不把论文摘要当实施证明。

## 3. 首次最小合法耦合

先做**冻结直前沿的三维挤出 Mode III 解析边值问题**，再回到既定倾斜 penny crack。域为半径 R、长度 L 的圆柱，裂面是负 x 径向半平面在柱内的部分，前沿为轴线；外边界施加下式精确应力的牵引，包括两端盖，并以六个积分约束去刚体模态：

```text
u=(0,0,A√r sin(θ/2)), A=(KIII/μ)√(2/π), b=0.
σxz=μ∂x uz, σyz=μ∂y uz; crack-face traction=0.
Uexact = μπA²RL/4 = KIII²RL/(2μ).
```

这是由调和场直接推导的工程验证题，不是论文复现，也不验证闭合曲前沿。验证能量范数/位移、裂面迹、力矩/力平衡、能量和纯 III 系数；负 x 两侧 jump 非零且按 √r 闭合，正 x 两侧连续。解析应力边载避免把 MLS 系数误当物理节点位移。先独立验证 quadrature，再固定几何做 h 收敛；不要只用“所有节点全局增广后精确包含答案”的一次代数通过代替局部增广验证。

复用：无 mask 的 MLS 值/解析梯度、Wendland、D/engineering-strain B 和 `BᵀDB` 内核、物理位移映射原理、显式 Γ/前沿/lineage、可见性作为几何诊断。必须新建增广 DOF 注册与 provider：现有 `3*node_count` 和“每行等于节点数”的接口不能伪装承载新增系数。补全裂面分侧与前沿奇异 quadrature、边载、积分约束。现有 visibility offset 不参与新基底，也不复制实物质量。

## 4. 极限与拒绝标准

* **固定 h 的函数空间门槛：** 对每个基函数/非仿射系数，沿 shadow 面、支割延长线、图册边界、support 边界取 ε→0 双侧迹。Γ 外差应趋零；有非零极限即拒绝，不能以“数值较小”豁免。Γ 上只允许声明的跳跃；tip 场检查有限能量。
* **积分与网格分开：** 固定点集/Γ 先细化 quadrature，之后细化 h，并旋转非对称点布局；比较能量、位移、平衡、裂面 traction 残差。拒绝额外零模、负能量、条件恶化被静默截断或结果随积分改变不收敛。牵引自由是弱式/收敛结果，不能要求粗网格每一点精确为零。
* **有限曲前沿：** 直前沿通过后做冻结倾斜 penny；独立细化圆周/裂面 tessellation 与 h，检查一圈局部 frame/trace 连续、几何扰动 δ→0 的响应极限、SIF 提取域变化后的稳定性（采用与曲前沿项一致的公式）。不得把无限体 SIF 当有限立方体逐点精确解。真实传播路径的 Δa→0 检查留到静态门槛后；不凭一个更新宣称路径收敛。
* 本审查不发明百分比/容差。主审须先冻结无量纲误差指标、独立参考、数值精度预算和运行上限；预算不足输出 UNKNOWN。任何 ghost jump、未经定义的曲前沿图册或未收敛积分均禁止 promotion；静态内裂面亦不产生两个碎片的 R2。

## 证据及主审待决

[Ai 官方 thesis](https://etheses.dur.ac.uk/12679/1/WAi_PhDthesis.pdf)，§6.2、印刷 pp.128–130（本地选页文本 pp.153–155）：确实规定有限 visibility、前沿点不 split、近前沿 cracking-particle support 全截断和 visibility-only offset；这些段落没有消除 shadow jump 的条件。此结论针对当前数学构造，不是声称作者全部结果无效。

[Moës–Dolbow–Belytschko 原论文全文镜像](https://paperzz.com/doc/9161143/a-finite-element-method-for-crack-growth-without-remeshing)，§3.1 pp.135–138/Eqs.(15)–(20)：支持域切分决定 H，含 tip 支持用分支函数；§3.3 p.139 专门处理裂面切割积分。该文是 2D FE，不能据此声称当前 3D EFGM 已复现；上述 H¹ 判断与 Mode III 题为本报告推导。

[Belytschko–Fleming 1999](https://www.sciencedirect.com/science/article/pii/S0045794998002053)，Introduction 对 §4 的说明，明确原 visibility 在非凸边界引起内部不连续，并描述 diffraction/transparency/see-through。当前仅检索到 abstract/introduction，不足以提供其完整算法。

主审只需决定：①接受从 bool-masked CPM 实施草案转为 PU-enriched EFGM（路线 B 不变）；②批准直前沿 Mode III 为 penny 前置耦合门槛；③冻结新 provider/DOF 合同、误差指标与运行预算。不得以本报告本身批准 solver run 或 GT。输入哈希及未执行事项见 `DECISION.json`。
