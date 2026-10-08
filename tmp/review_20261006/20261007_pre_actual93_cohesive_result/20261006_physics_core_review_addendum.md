# 2026-10-06 物理核心补记（2026-10-07 DRAFT）

**主线与资格。**RelicWeave核心仍是“观察碎片→候选排序→证据→接受/拒绝/补采/STOP”，GT仅离线审计；完整物理包0。A12 30–40%、VF 30–40%、PG 25–35%为主观工程进度区间，不是物理完成度或测试通过率。物理突破前总DOCX不更新；工程方法改进已记录，论文新颖性尚待文献对照。Git独立postverify确认58文件/959,548规划字节逐字节复核通过：PG提交6aabf6191eb232b4ad446334189d15e894bc3bff，VF提交8acff8ac7db9cf212ed87ede3054a8311b69cc3a；远端匹配，独立blob复核通过。

**VF决策与内壁offset。**Jar外轮廓10个非平直span的R9目标及独立C4/feature/regular证书已接受；优化器80次迭代未收敛，producer未晋级。Bowl拒绝图保留−300 mm反例及6个精确R9反例；停止继续SLSQP补丁。新保守LP的bowl、jar均为solver status 2、无候选；这不证明所有器形不可行，停止此LP。旧offset run把z传给需要z′的算子，rc2/0.078s；2,890 B失败产物与source freeze保留。新z′修复run rc0/0.187s；独立saved-offset检查rc0/0.175774s，12 spans、11 joins及300个display点均通过且误差0。主审接受Jar局部二维法线offset候选；新布局图等比例坐标、端点与图例可辨。范围仅jar local 2D；全局strip/wall/lip/native/physics仍UNKNOWN。

**PG、接口与93桥。**CZM通用toy接口含21次assembly、3分支/24 DOF、非原点全交叉FD相对误差1.86e−11；仅接口数学验证，不代表实际器物耦合。连续折线v2单次5测试产生7条返回：同顶点集有正连续距离；平行偏移4样本精确0.5 m，刚体pose 6样本精确0.25 m；极小atol到128样本仍UNKNOWN，真实mesh bias未验证。当前93 cohesive bridge worker/supervisor源审已接受为SOURCE_READY_WITH_EXTERNAL_ATTEMPT_WATCHDOG_REQUIRED；外部watchdog单次执行正在准备，数值结果待验收。旧metadata-only review仅记录早期driver，不代表当前runner数值PASS。

MLS rotation focused run rc0/0.188696s：64 queries、5项检查；4组各16项且support=27，最大condition=18.9585。proper rotation权重误差1.665e−16、归一化梯度误差6.395e−16；mm权重/梯度误差1.11e−16/2.665e−16；anisotropic负对照delta=0.0892173，大于1e−8。仅接受为provider covariance prerequisite；不证明四面体/裂纹/网格方向偏差消除，力学多方向对照仍UNKNOWN。

**缓存数学检查与待补证。**93缓存producer曾在0.6208 s、rc1于JSON bool序列化失败，arrays.npz已保存而原result.json缺失；serializer-only 0.129 s检查PASS。独立saved-93 reader rc0/0.296708 s，接受为缓存制造解数学恢复：6零/87正/无显著负谱，公式复现≤7.08e−16，f/K独立重算误差0、KKT恢复0；工程场误差0.005431783（0.54318%，纯I体积strain），off-grid stress 0.1658–0.2294%。不把原失败运行改标成功，不认定物理碎裂。mass saved-audit数学检查rc0/0.1681s：M93全正、零/负特征值各0且Cholesky成功；六个刚体速度点误差8.88e−16，6×6动能/惯量归一化最大1.728e−15。rho=1 toy、NOT_MATERIAL_CALIBRATED，不是动力学冲击/接触证据，物理包仍0。移动前沿与材料标定方案仍是source-only设计，无数值PASS。

**下一前沿。**需要真实裂纹前沿与能量耦合、动力学及碎片接触、材料标定、薄壁和网格方向对照、粗糙断口、纹样/材质退化详细标签及A12资格验证。局部offset、协变前置检查、缓存制造解与离散接口都不替代这些证据。本补记保持DRAFT，待本轮PG cohesive外部watchdog数值结果与主审再验收后冻结。