# Mixed branch enrichment rank review

**范围：**只解析审阅 `mixed_enrichment_algebra.py`、`mixed_lefm_field.py`；未调用 MLS、组装、积分或求解器。令局部坐标为 x,y，a=theta/2，F_j=sqrt(r/R)h_j，h=(sin a,cos a,sin a sin(theta),cos a sin(theta)。

## 全局每源分支

若同一完整节点集再现局部坐标（sum_i x_iN_i=x、sum_i y_iN_i=y），则每个位移分量有至少两条精确零系数方向：

\[
\sum_i(x_iN_iF_3+y_iN_iF_1-y_iN_iF_4)=xF_3+y(F_1-F_4)=0,
\quad
\sum_i(x_iN_iF_4+y_iN_iF_3-y_iN_iF_2)=xF_4+y(F_3-F_2)=0.
\]

由 sin(theta)=2sin(a)cos(a) 得。非退化节点下三个位移分量共至少六个零系数模，产生零位移、零应变，精确积分也无法消除。

仅取前平面仿射系数 p_j=a_j+b_jx+c_jy 时，12个 branch 函数维数为10：四个F独立，八个xF_j,yF_j恰有上面两条关系；展开 h_3=(cos(a)-cos(3a))/2、h_4=(sin(a)+sin(3a))/2，比对半角1/2、3/2、5/2谐波可验证无其他关系。完整3D仿射系数再加四个独立 zF_j，维数为14。此计数不等于一般MLS系数场的秩。

只富集局部节点子集时，完整节点集再现坐标不保证子集仍再现；若局部活动族本身再现x,y，则相同零模仍存在。shifted enrichment 是不同基，需另析。

## 能量核与LEFM

连通、正定弹性体的物理零应变场只有六刚体运动，三维线性再现使BASE可表示它们。系数核还含上述至少六个零场方向，故六刚体不是系数空间唯一零模。所读源未定义MLS、组装和边界，故不能确定矩阵总秩。

所给三个直裂纹前沿各向同性 leading LEFM 位移均可由四分支的常系数精确表示。令 kappa=3-4nu，A_I=K_I sqrt(R)/(2mu sqrt(2pi))，A_II=K_II sqrt(R)/(2mu sqrt(2pi))，A_III=K_III sqrt(2R/pi)/mu：

\[
I:(u_x,u_y)=A_I((\kappa-1)F_2+F_3,(\kappa+1)F_1-F_4),\quad
II:(u_x,u_y)=A_{II}((\kappa+1)F_1+F_4,(1-\kappa)F_2+F_3),\quad
III:u_z=A_{III}F_1.
\]

因此制造场首链用27个BASE源加每分量4个全局常系数branch（93主系数）是紧凑且足够的；它保留四函数并精确表达混合I/II/III，不声称通用XFEM完备。全牵引边界保留六个物理刚体模，故93主变量加六个KKT规范乘子为99阶合理。不要用ridge、伪逆或暗截断隐藏原405向量参数化的已知零模。

**未决：**未给出MLS支持、节点几何、其他依赖与完整405系数总秩；需查实际形函数/组装。

