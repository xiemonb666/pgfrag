# Tip-aligned quadrature math review

只读审阅全3D算子、静力桥、fixed-front case和已冻结结果；未运行数值代码。首版映射与最终固定c端盖映射在几何和Jacobian上成立，source review待main给出READY地址。

体积：`x=a+ρcosθ, y=ρsinθ, ρmax=−a cosθ+sqrt(R²−a²sin²θ)`，`ρ=ρmax s²`，故`J_V=ρ dρ/ds=2ρmax²s³`（m²）；乘`dz`后为m³。对`a<R`沿各θ从tip到圆周，覆盖整个圆盘。新tip`gradF1~ρ^-1/2`，映射后的刚度积可积；θ内点不踩`±π`支路。

端盖固定分界`c=0.0005m`对两种a及所有θ级复用同一物理点集：内环`r=c sin²(πs/2)`，`J_A=r dr/ds=(πc²/2)sin²(πs/2)sinπs`；外环`r=c+(R−c)s²`，`J_A=2(R−c)s[c+(R−c)s²]`，均为m²。二者依次覆盖`[0,c]`、`[c,R]`，仅在零测度环相接。内环在s=0使旧tip traction的面积加权量为O(s²)；内外环两端均聚焦r=c的新tip端盖尖点。上下半圆覆盖全角域并分开y=0两支，不查询端点。旧tip边载始终用冻结的`σ_old n`，a变化只令新基下的`f_a`变化。Side/ cap geometry不随tip平移。

总4992行分20块、每块≤256可容纳（20×256=5120）；case须记录逐块行数/唯一点、边界载荷复用规则及数组峰值。通过守恒或加密稳定只证明离散一致性，不证明材料或裂纹物理资格。
