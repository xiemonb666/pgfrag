# 最终独立窄审：能量门与 mesh binding

结论：窄接受为 fixture 状态机与几何绑定 plumbing；没有实际断裂力学/J-integral 耦合，也没有物理断裂面或粗糙度识别证据。下一步应让真实 solver 的外功、储能变化、裂面能和 J/前沿证据带可验证来源进入门控，并校准数值容差；mesh 通用路径还需逐连通壳方向检查及相交/覆盖验证。

v2 实现按 R = W_increment − ΔU_stored − Gc·ΔA 检查，容差为 abs_tol + relative_tol·reference_energy；派生 Δtotal potential 只报告，不重复扣减。11 项收据全通过，9 个保存状态的算术可复核：两类大残差拒绝为 NEEDS_RECOMPUTE，旧 active surface 保留；wrong-parent、缺容差及 legacy 输入也拒绝并保留父面，动态模式为 UNKNOWN。balanced 状态仅 COMMITTED_FIXTURE_ONLY；synthetic external 分支仅 COMMITTED_EXPERIMENTAL_REANALYSIS。caller 提交的哈希/来源只是绑定声明，未独立验证；容差 basis/hash 未经校准。原 v1 当前源 SHA256 与 v1 旧收据一致：3576a2aa6a68771578ed93e70d618f893b3465db476006f98e811ab12b44fba2。

提示 JSON 中缺容差、动态、wrong-parent 三行的 agrees=false：这表示门控提前停止且 saved residual=null，脚本拿现有能量项算出的 0 不可与未计算 residual 比较；应读作 NOT_APPLICABLE_BLOCKED_EARLY，不是数值不一致。未重跑原件、未覆盖 JSON。

冻结 planar BOX 的两 OBJ 哈希与内嵌顶点/面、接口 pair ledger 绑定；12 三角面/件、封闭边、正有符号体积 4.44 与 3.56 m³，和为 8 m³。两对 cut 面各 2.1774296774 m²、法向点积约 −1；该盒体 fixture 的连通/朝向由冻结拓扑读回确认。体积和不证明一般分片无重叠或完整覆盖；通用 validator 的整体正体积也不能保证多个连通壳各自外向或排除自交。输出限定 FIXTURE_GEOMETRY_ONLY、NOT_COUPLED、physical_pack_count=0，平面切割不能作为真实粗糙裂面 GT。

独立数值与 SHA256 清单：[independent_readback.json](independent_readback.json)，SHA256 81215da6a5999059f33dd56d2ea4cf46e4e85f9ca86b32d0d7d22d47e6badbb8。既有 mesh 审阅/实际收据 SHA256：d94278c391a9afc3d9ec3bb7c668d610a399cd8df18b7341f68b2f38cebb73c8 / e14f8eae0ccdeab5b8e909d7748dd225187b122d5cdead2b967129e4a3db81f1。
