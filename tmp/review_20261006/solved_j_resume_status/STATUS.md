# Solved J 状态

Primary `39cb…24d3`；stage `1377…fac9`。Stage 未限定：8 轮廓、224 点/解档、448 梯度样本；两档 54 项 ID、d 与主结果相同。

算术 PASS：独立 `math.fsum` 与保存 J 最大差 `1.734723475976807e−18 J/m²`；解析 `0.013` 仅比较。更正值 `0.0006968643143483756` 是两档静力解差，不是轮廓 GL12/16 差。

Receipt 未通过：父结果仍 `UNKNOWN_WORKER_RESULT_CONTRACT`，worker rc=0。mapper 查缺失的 `quadrature_rule`；源字段为 `theta_order`/`volume_point_count`。候选映射未应用；stage 记 `contract_evaluation_count=0`。两次 JSON 复核失败；fsum 输入未变。

物理边界：中面 `|grad_z|≈1e−17`，镜像控制粗/细档 `2.3473e−5`/`2.2475e−7`，不能证明整域 z 不变；标量 `u_z` 加中面二维 J 不是全 3D 断裂能/G。

Full-3D 缺三接口：矢量位移/应力场、前缘局部 3D domain-J、裂纹/材料本构与能量平衡。
