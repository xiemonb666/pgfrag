# Tip-aligned study metadata review

读过`study_case.json`和`README.md`；未读在写中的Python源码，未运行数值代码。

metadata已正确区分：固定物理圆柱/旧tip全应力牵引与变化的广义载荷`f_a`；旧tipexact只用于a=0场误差，a=.5只报告相同点的解场差；`G_delta_a`与旧tip局部解析`0.013 J/m²`不作直接pass对比。fixed-c端盖映射、θ支路分段、4组体积数和4992/20块预算自洽。

需改三处措辞：`geometry_and_physics.crack_and_branch_warning`、`limitations[0]`和README第2段称`0<x<a`为trial basis与old physical crack cut不匹配。按候选几何定义，`S_a={(x,0,z):x<a, -L/2≤z≤L/2}`，试函数F_a的支路正对应`Ω\\S_a`；没有候选空间内部的裂面不匹配。应明确a>0仅是未提交的候选前沿再分析：旧tip场只在`S_0`是exact解，在a>0仅生成冻结的外部牵引，不存在被声称的shifted exact误差或自动裂纹状态。

READY源码审阅重点：确认F_a在`S_a`两侧提供允许跳量且弱式不跨裂面强制连续；volume basis使用当前a；外载traction始终由旧tip场及原物理法向产生，`f_a`用新基重积；a>0不把旧场误差/旧基梯度当新tip exact；旧场仅按case规定用于a=0对比及所有物理边界traction。
