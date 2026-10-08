# 方法修订：连续 MLS × front PU，Mode III

**仅设计，未运行。** 接受 root 决策：停止 sharp-bool visibility 力学基底，保留显式裂面 B，改用连续 MLS 加合法前沿 PU 增广。旧 `COUPLING_PILOT_DESIGN.*` 原样保留作历史，不再作为实施路线。

## 固定制造解

圆柱 `r≤R=0.01m, 0≤z≤L=0.01m`，单一直前沿为 z 轴；`Γ: y=0,x<0` 径向裂片延至圆柱外壁和两端盖。全域用同一 z 轴前沿坐标 `(r,θ)`、`θ∈(-π,π)`，禁止 nearest-triangle frame 切换。固定 `E=1e6Pa, ν=.3, KIII=100Pa√m`，`μ=E/[2(1+ν)]`：

```text
uz=A√r sin(θ/2), ux=uy=0, A=(KIII/μ)√(2/π), b=0
σxz=μ∂xuz, σyz=μ∂yuz, Uexact=KIII²RL/(2μ)
```

外圆柱侧和两端盖均施加解析 `t=σn`（端盖含 xz、yz 剪力）；裂面零牵引；六个积分约束去刚体模态。此为制造解，`r→0` 应变发散，不能称真实材料或整域小应变物理。记录实际 `r_min` 及任意数值截断；禁止静默clip或称其材料校准。正半径排除区能量占比 `r_min/R`，须补解析值或报告漏项。

## 共享接口与验收

Provider 返回独立 `BasisRegistry`：连续MLS、`J_H`、`T_F` 的全局唯一 DOF、基函数类型、分量与来源；`ShapeRow(x)` 返回物理基值、材料梯度、DOF IDs、成员证据及状态。新增系数不得伪装粒子或混入 `3*node_count`。`J_H` 仅收全support被真实Γ切成两非空部分且support不含front的基函数；否则禁用。首阶front仅 `F1=√r sin(θ/2)e_z`，`T_F` 覆盖front；梯度用乘积法则。若小布局自然导致所有site均在 `T_F`，精确复现只算代数桥接，不能证明局部增广收敛。

使用front-aware体积/边界积分并分辨 `r^-1/2` 梯度；root须先冻结点布局、support、`J_H/T_F`、quadrature、`r_min`、误差指标及DOF/内存/时间上限，未冻结则 `UNKNOWN`。验收位移/应力、裂面迹与零牵引、两端盖及圆柱侧加载、力矩平衡、仅六个刚体零模、无额外零模、能量趋近 `Uexact`；先固定几何收敛积分，再细化点云并旋转非对称布局。曲前沿、混合模态、SIF及传播留后续门槛。不得由此宣称真实破裂、R2或GT。

决策文件及SHA-256见JSON；未批准运行。
