# 2026-10-06 物理主线状态

**决策。**优先物理链：full-3D弹性→校准裂纹扩展→断面/碎片包。核心仍是“观察碎片→排序候选→证据→接受/拒绝/补采/STOP”，GT离线审计。合格完整包0；A12 focused25、退化primitive12 PASS；P02/P03端到端未通过。

**已运行。**full3D静力桥run02退出0、0.754秒，pins不变；27点、162 DOF（6刚体零模态、156正模态）。位移/应变误差1.90%/5.31%；U=1.31335e−6 J，对解析1.3e−6 J为1.03%。独立审计质量/惯量一致，质量矩阵相对差1.84e−16、条件估计2.29e7、应力L2误差8.09%。ωmax=2.187e6、dt上界9.144e−7 s仅适用玩具ρ=1。主算σyy=244.306 Pa，对照剪应力1513.807 Pa；制造载荷含两端盖面内牵引，该差异不证明载荷错误。

fixed-front有限扩展G=0.005466 J/m²与同一材料/制造基准局部解析0.013可作诊断比较但非同一极限；δa/R=5%，旧tip积分/场误差未收敛，约58%差距不能归因物理误差。独立复核缓存数组20/20代数项PASS（0.207s），仅审旧数组、无物理解算；渲染0.154s。Gc未比、增长未提交。mixed-LEFM字段54点独立审计通过；新四分支helper 7项算术检查PASS（0.107323s），11坐标、18次helper调用、405语义ID清点；现有算子上限375 DOF，故混合solver未运行。

**当前物理门。**优先重做裂尖对齐奇异积分与场/能量精度验证。约千行tip-quadrature source已root审阅；trust gate、early overwrite、adapter basis-row修正进行中，数值NOT RUN。仍缺h/dt/方向偏差研究、材料标定、真实裂口、自动扩展、接触/内聚耦合、物理碎片数/断面/事件谱系及完整标签包。两片箱体OBJ仅FIXTURE_GEOMETRY_ONLY/NOT_COUPLED；不得用bond components代替碎片。

**VesselForge候选。**11控制点global min-jerk两家族实际4分段KKT，0.208764s；9等式/4严格约束、C4、hull/regularity及独立24点de Boor抽查通过。R9目标12段9过3个端点反例：bowl轴R4.699 mm，κ=+0.21280775646/mm；bowl口R5.435、jar口R2.131 mm均负κ，三处1−3κ均正；不降低R9。现有轮廓缺明确平底/neck。17控制semantic flat-base+cyl-neck版仅SOURCE_ONLY、未运行。旧catalog曲率证书另为12段4过8败、偏移充分屏7过1败；失败不证明尖点，signed regularity及整体壁厚/自交仍UNKNOWN。VF仅一族生产处理，smooth CAD实验；两帽网格仅表面级，bowl/jar闭合STEP不构成物理断裂母体。

MLS替代面追踪、mesh-bias资格门、signed curvature区分、全局B-spline均是工程选项/待论文验证，未确认新颖性。合格完整物理包经标定和独立复现后再更新总DOCX。

证据：pgfrag/tmp/review_20261006/{full3d_static_solve_bridge_20261006,bridge_mass_independent,fixed_load_front_saved_independent,mixed_lefm_enrichment_algebra_20261006,tip_quadrature_independent}；vesselforge/tmp/review_20261006/{global_bspline_independent,catalog_corner_fairing_20261006,catalog_corner_fairing_independent}；tmp/review_20261006/tip_aligned_quadrature_study_20261006。