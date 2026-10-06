# 2026-10-06 物理主线状态

**决策。**优先物理链：full-3D弹性→校准裂纹扩展→断面/碎片包。核心仍是“观察碎片→排序候选→证据→接受/拒绝/补采/STOP”，GT离线审计。合格完整包0；A12 focused25、退化primitive12 PASS；P02/P03端到端未通过。

**已运行。**full3D静力桥run02退出0、0.754秒，pins前后不变；27点、162 DOF（6刚体零模态、156正模态）。位移/应变误差1.90%/5.31%；U=1.31335e−6 J，对解析1.3e−6 J为1.03%。独立审计质量/惯量一致，质量矩阵相对差1.84e−16、条件估计2.29e7、应力L2误差8.09%。ωmax=2.187e6 rad/s、dt上界9.144e−7 s仅适用玩具ρ=1，非物理时间步。主算法向σyy=244.306 Pa，对照剪应力1513.807 Pa（比值16.14%）。制造载荷包含两端盖面内牵引；该应力差异本身不证明载荷错误。

fixed-front run01 rc=0、0.551s，固定牵引，δa=0.5mm；新增1次K+KKT，4残差≤1e−9。Gfinite=0.005466 J/m² vs同一材料/制造基准局部解析0.013，可比诊断但非同一极限；δa/R=5%，旧tip积分与场误差未收敛，约58%差距尚不能归因物理误差。Gc未比、增长未提交。mixed-LEFM源ed9006…的54点解析/旋转独立审计通过，无solver耦合、标定或裂纹扩展。

**缺口与下个物理门。**仍缺裂尖积分收敛、h/dt/方向偏差研究、材料实验标定、真实裂口、自动扩展、接触/内聚耦合、求解所得碎片数/断面/事件谱系及完整标签包。本轮优先重做裂尖对齐奇异积分并验证场/能量精度，再补材料基准与裂纹起始/扩展/断面验证；随后核验接触与碎片质量、惯量、谱系，才供A12形成R2及候选证据。不得把bond components当物理碎片。两片箱体OBJ为FIXTURE_GEOMETRY_ONLY/NOT_COUPLED。

**器物与创新边界。**VF仅一族生产处理，smooth CAD实验。bowl/jar闭合STEP非物理断裂母体；壁厚/自交未知，两帽网格仅表面级。精确曲率半径目标12段4过8败；绝对曲率充分屏7过1败，失败不证明尖点，signed regularity未知。全局五次B-spline仅初始约束源就绪，确定性fairness源未完；未求控制网/出曲线。MLS替代面追踪、mesh-bias门、signed曲率区分、全局B-spline均为工程选项/待论文验证，未确认新颖性。总DOCX待物理突破后更新。

见：pgfrag/tmp/review_20261006/{full3d_static_solve_bridge_20261006,bridge_mass_independent,fixed_load_front_reanalysis_20261006,mixed_lefm_independent,fracture_artifact_contract_20261006}；vesselforge/tmp/review_20261006/{exact_curvature_certificate_20261006,catalog_native_adapter_20261006}；tmp/review_20261006/global_bspline_meridian_20261006。