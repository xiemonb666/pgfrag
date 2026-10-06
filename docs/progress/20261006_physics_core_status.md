# 2026-10-06 物理主线状态

**方向。**核心链为“观察碎片→排序候选→证据→接受/拒绝/补采/STOP”，GT离线审计；合格完整物理包0。A12 focused25、退化primitive12 PASS，P02/P03端到端未通过。本轮非物理突破，总DOCX不更新。

**弹性/裂尖。**Full-3D桥run02 rc0/0.754s，27点、162 DOF（6刚体+156正模态），位移/应变误差1.90%/5.31%，U=1.31335e−6 J对1.3e−6差1.03%；独立质量/惯量一致，质量矩阵差1.84e−16、条件2.29e7、应力L2误差8.09%。载荷含两端盖面内牵引；σyy=244.306 Pa对剪应力1513.807 Pa不证明载荷错。

Tip-aligned四case outer 3.531s（worker 2.953s），raw 4992/unique4224/17批；每谱6零、156正、无显著负（阈值64ε·162，非对称门1e−12）。a=0体积应力误差1.345%/1.347%，128缓存独立位置0.51265%/0.51279%；更高θ无改善非收敛证明。U=1.300233e−6对1.3e−6 J约0.018%。外牵引沿用a=0制造解解析完整牵引；a=0.5 mm为手动候选前沿，精确场未知。新G=0.01411668/0.01411683 J/m²（一侧L·δa），不以局部0.013作PASS。独立raw K相对差4.6–6.09e−16，f/C0/RBC/能量及各门一致，rc0/约1.17s。

内聚toy assembly单次运行含17次assembly调用，rc0，outer 0.336s/internal 0.00435s；弹性/软化/卸载FD相对最大2.33e−12，耗散200 J=2 m²×Gc 100 J/m²，ledger误差2.84e−14 J，力/矩残差0。源码/手算/图审通过；仅节点toy，非实际器物耦合。裂尖精度待验证；mesh bias、h/dt、方向、真实材料/裂纹/扩展、碎片及真实接触耦合仍UNKNOWN；两帽OBJ仅表面fixture。

**混合场。**4分支helper 7算术项PASS（0.107323s），11坐标/18调用/405语义ID；每源增广至少6额外系数零场（每分量2个跨分支依赖）。93/405语义raw operator source通过；独立K tiny 24/60 rc0、0.1076s，24×24 raw K逐项误差0、非对称相对4.84e−17，6跨分支零场候选最大Knull 1.97e−16。27源仅IDs，无93谱/求解；93 cached仅source preparation，现存a0 θ32：volume512+boundary768+offgrid128，0新MLS。27 BASE+4 GLOBAL/组件、93 DOF/99 KKT仍未解算；保留405诊断，不用ridge/pinv。mixed-LEFM 54点独立审计通过；内聚未耦合器物。

**VesselForge。**17控制flat-base/cyl-neck两族共4次KKT，0.413909s，特征/严格约束/C4/hull/regularity过；新20非平直目标13 PROVED/4 FAIL（2独立bowl点）/3 UNKNOWN，另计4直段。旧R9反例：bowl轴R4.699（κ=+0.21280775646/mm）、bowl口R5.435及jar口R2.131（后两者κ负），三处inward factor正。旧bowl apex→neck Δr=3、H=8.6 mm，小于R9单调/竖切线必要界√99=9.95 mm，总变差证明成立。新apex z=72/H=14仅必要界/source，非优化PASS；源码修正待最终审阅，优化器未运行。全局fairness未完成；壁厚/自交/signed regularity UNKNOWN。

MLS、mesh-bias门、signed curvature区分、全局B-spline均为待论文验证的工程选择。当前qualified packs仍0。

证据：pgfrag/tmp/review_20261006/{fixed_load_front_saved_independent,tip_aligned_saved_independent,cohesive_interface_assembly_20261006,mixed_enriched_operator_20261006,mixed_enrichment_rank_critic}；vesselforge/tmp/review_20261006/{semantic_base_neck_bspline_20261006,semantic_base_neck_independent}。