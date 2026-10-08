# Mixed93 最终独立源码审阅

**SOURCE_READY；未发现剩余源码阻断。** Driver SHA256 `61e88d32a6d5543c891e92c8996b66f1562db6a20781db061c6f266d0dffbfa8`；case `ab57f5969a22218f4c347633ca312029ffd8ac6446200eb6c9481529e8a1c5d0`；README `a7083928fdb77fdd42a19539a7f392af4287f7faeb5c906540d10186dfdd4ba9`；SOURCE_PINS `b93cce82be3b4023a4de333b7f8755b55a44275e052751a9b38cc37fdcccb40f`。

语义布局为27个排序BASE×3加4个global branch×3；F1–F4系数符号、归一化和branch-major槽位与冻结公式一致。B列、全6×6各向同性D、engineering shear、应力及侧面/两端盖全量`σn`一致；raw K不对称化、不截断，无ridge/pinv/载荷投影。缓存result元数据为4992 raw、4224 unique、17批，MLS/旧场各4224行均OK；选点512体积、256侧面、512端盖、128 offgrid，case/hash/result内array SHA与row-ID映射相符；未读NPZ数据。

物理矩阵实现`V⁻¹∫[u,(x×u)/R]dV`，分量符号和原点中心正确；`g=Cphys q_exact`与解析场直接积分核对。对角缩放谱/非对称门、刚体秩空域、99×99缩放KKT的q/λ还原、物理弱式/约束残差及`2U−q·f+g·λ`一致。B和四个q_exact在体积、边界、offgrid三云均分块≤256并先于首个solve复核。2%仅是工程诊断，不代表物理资格。

当前case SHA与SOURCE_PINS current/freeze一致，执行授权均为false；默认源检查路径先于NumPy/cache读取。执行前会核当前case、授权位、门值、entry/依赖/缓存哈希，输出限于本study新目录。未来GO需按README人工审阅仅授权字段及current pin更新。未运行数值代码/求解/CAD/Git，也未改source/case/pins；只读缓存result元数据JSON。