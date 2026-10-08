# 独立静态与保存数组审阅

## Local-front saved K

输入保存件 SHA-256 为 `EBDEA44857AB2B1A358FC679D0B1E4996943905A5996BA83F8F8F1BA68340762`；冻结 operator 源 SHA 为 `77B0DB5F26855FFB6D444154D42B04C5644FC99737BE177420A3EFAD76AC9E59`，freeze JSON SHA 为 `8B7E29A701A53EC9F00F24CB5ED2B45D80FFC470DF2CE225DBBF4380026424AF`。只读重算保存的 140×140 K，没有导入 provider/operator、重组 K 或求解。

窄数值回读通过：K 有限且对称，重算矩阵 SHA 与收据相同；`eigvalsh` 与保存谱最大差为 0。谱有 1 个舍入级负特征值（最小 `−7.50e−13 N/m`），0 个超出 `1.4091e−7 N/m` 阈值的负值，恰有 1 个预期平移零模。由 canonical registry 重建的 100/30/10 个 BASE/F1/Heaviside DOF 顺序唯一；平移缩放残差 `3.52e−15`，归一化均值行对 1 m 平移给出 1 m。256 个权重均为正，总量 `6.4e−6 m³`。从固定晶格及语义键重建 `Uz=x` 后，`0.5 dᵀKd=1.2307692307692366 J`，与解析 `1.2307692307692308 J` 相对差 `4.69e−15`。

这只接受保存矩阵的窄代数回读。**原 OUTER125 仍是失败且未晋升**：worker 发出的小写 basis 计数字段与父校验器期待的枚举名不一致；该失败必须保留。本次未运行 MLS、operator、solver、测试或渲染。审计 JSON和谱图分别为 [local_saved_array_audit.json](local_saved_array_audit.json) 与 [local_saved_K_spectrum.png](local_saved_K_spectrum.png)。外层命令的子进程退出码和耗时没有保存在回执中，本报告不推断这两个值。

## Full3D

最终 source freeze 尚未提供，因此**暂不作 READY 判定**。未冻结快照中发现一个待核 promotion blocker：完整 `_promotable_success_result` 校验器有定义，但未看到 parent `main()` 调用；父路径读取 staging JSON 后直接写 primary，且 finally 删除 staging。作者已在补此路径；收到最终 freeze 后只做一次窄核对，确认 rc0 才提升、失败谱保留 staging、且实际发射的 status/eligibility 值与校验器匹配。此前发现的端盖 `Mz` 公式和数值矩仍需以冻结版本确认。
