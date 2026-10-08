# 2026-10-06 物理核心补记

**主线与资格。**RelicWeave核心仍是“观察碎片→候选排序→证据→接受/拒绝/补采/STOP”，GT仅离线审计；完整物理包0。A12 30–40%、VF 30–40%、PG 25–35%是主观工程进度区间，不是物理完成度或测试通过率。物理突破前总DOCX不更新；已记录工程方法改进，论文新颖性尚待文献对照。Git已推送119文件（PG c13770…；VF f34800…），独立blob复核进行中；本补记为source-only draft，待root最终审阅。

**VF决策。**Jar外轮廓10个非平直span的R9目标及独立C4/feature/regular证书已接受；优化器80次迭代未收敛，producer未晋级。Bowl拒绝图保留−300 mm反例及6个精确R9反例；停止继续SLSQP补丁。新保守LP的bowl、jar均为solver status 2、无候选；这不证明所有器形不可行，停止此LP。

**内壁offset。**旧run把z传给需要z′的算子，rc2/0.078s；2,890 B失败产物与source freeze保留。新z′修复单次rc0/0.187s，12-span内壁local candidate的局部证书PASS：37/64 nodes、300次display、Q degree18且每span 19系数。独立saved-Q与视觉复核待审，不能先标accept；global strip/wall/lip/native/physics均UNKNOWN。

**PG与接口。**CZM通用接口toy含21次assembly、3分支/24 DOF、非原点全交叉FD相对误差1.86e−11；仅接口数学验证，不代表实际93 DOF求解或器物耦合。生产侧当前仅static checker状态SOURCE_READY_METADATA_ONLY；实际数值bridge由另一agent编写，未数值运行。连续折线v2单次5测试产生7条返回：同顶点集有正连续距离；平行偏移4样本精确0.5 m，刚体pose 6样本精确0.25 m；极小atol到128样本仍UNKNOWN。范围只覆盖给定折线，真实mesh bias仍未知。

93缓存producer在0.6208 s、rc1于JSON bool序列化失败，arrays.npz已保存，原result.json缺失；serializer-only 0.129 s检查PASS。独立saved-93 reader rc0/0.296708 s，接受为缓存制造解的数学恢复：6零/87正/无显著负谱，公式复现≤7.08e−16，f/K独立重算误差0、KKT恢复0；工程场误差0.005431783（0.54318%，纯I体积strain），off-grid stress 0.1658–0.2294%。不把原失败运行改标成功，不认定物理碎裂。mass saved-audit数学检查已接受（rc0/0.1681s）：M93全正、零/负特征值各0且Cholesky成功，六个刚体速度点误差8.88e−16，6×6动能/惯量归一化最大1.728e−15。ρ=1 toy、NOT_MATERIAL_CALIBRATED，不是动力学冲击/接触证据，物理包仍0。

**下一前沿。**需要真实裂纹前沿与能量耦合、动力学及碎片接触、材料标定、薄壁和网格方向对照、粗糙断口、纹样/材质退化的详细标签，以及A12资格验证。接受的离散接口、折线度量和几何证书都不替代这些物理证据。