# Mixed93 saved cohesive result：只读独立审阅

结论：`READONLY_DIAGNOSTIC_CONSISTENT_PHYSICAL_SCOPE_UNQUALIFIED`。保存的实际包与冻结 93-DOF Mode-I 固定裂面合同在我检查的几何、符号、history、能量配对和调用计数上自洽；结果仍只属 toy diagnostic，不是材料、裂纹扩展或物理通过。

审阅输入身份：worker SHA-256 `36ac1710286b0abf607bda5b1c6ffab893f8191830eae091517e2acd41d62236`；case `fcc4b46ed101b0f6067a9f3ec10816909c9fa931edd5a8cdde4c7f46cb3f3a71`；run01 `result.json` 1,440,429 B，SHA `e0819a792789fe8fdada86a75983e7cb8465c5d0aa78140ca93a84c4440159f2`。它与 inner receipt 所记结果 SHA/长度相符；inner receipt SHA `eccd3510fa0c1281484b7ea97a57e226551e5c7d400dbd5da60b6d9772aa6474`。外层 attempt receipt SHA `955250d7cf20c8bc9e467b1c272e630c7e6f9b31270a07587dd4d0ce7cb8a187`；全程 capture SHA `5851460bcdd45dd2f03e605292e9da14aa325e92353e8b70a3e2cabfcf1db403`。inner 与外层记录均为 rc 0、无 timeout、结果 `ACCEPTED_DIAGNOSTIC_ONLY`、前后 source pins 相同。capture 记完整外层 elapsed 0.972552 s、host command wall 1.164988 s（均低于 60 s）。NPZ 未由本审阅打开；仅沿用 receipt 中 raw-byte SHA `3e1e6fd17751e4e06113010f20ea33c7b69914ca90431e26bbf76124951b5d41`。

结果包含 16 条 MLS 查询、16 条 interface H+/H−/J、16 条 baseline 与 16 条 final point histories。MLS 保存为 16/16 `OK`、每点 support 27；归一化值矩、partition 和梯度矩误差分别约 `1.73e-16`、`4.44e-16`、`5.55e-16`。界面从 `x=-R..0, z=-L/2..L/2` 映射；首点 `x=-0.00930568 m, z=-0.00430568 m`，`n=+y`，单侧权重 `3.02507e-6 m²`。保存的首点 F1 `H+_uy82=+0.9646596062, H-_uy82=-0.9646596062`，故 `J_y82=1.9293192124`；F2 的 uy 跳跃为 0。相应 `q0_82=1.4522648224e-7 m` 与 `qc_82=9.1799132258e-9 m` 对应的保存 opening 分别为 `2.8018824e-7 m` 与 `1.7710983e-8 m`。这与 `KG[82,82]=20000 N/m`、其余为零及 +y 开口一致。全界面单侧面积记录和为 `9.999999999999999e-5 m²`，没有把两侧面积重复累计。

逐点 ledger 显示 baseline 与 final 均为正开口、`damage=0`、耗散为 0；开口范围分别为 `[7.6534e-8, 2.8019e-7] m`、`[4.8378e-9, 1.7711e-8] m`，低于 case 的 `delta0=1e-5 m`。首点 baseline primitive 内力 y 分量 `+8.4759040e-5 N`，物理上表面牵引 y 分量 `-28.018824 Pa`；它们的 units 与一次面积权重相配。下表面记录为相反物理牵引。final 仍为同号且弹性可微；剪切、压缩没有被投影处理。三次 primitive assembly（含已丢弃的 fully-damaged 零切线 control）与一次 99 阶 KKT solve 的尝试/成功计数一致；control 汇总力、切线范数均为 0，且明确不进物理解。

缩放为 `q_solution`、外力、物理 gauge target 与乘子乘 `0.01`，`q_exact` 只列诊断。保存的非零解变化范数 `3.5866e-7 m`，cohesive force 范数 `1.8360e-4 N`，远大于所记舍入 floor；这支持该 toy 接口确实改变本次固定裂面平衡解。保存的 KKT/弱式/gauge 相对误差分别为 `4.43e-17 / 1.51e-17 / 2.40e-17`；输出还报告 K、K+KG 均有 6 零/87 正/0 负。**这些 spectrum 与残差数字是读取已保存输出，本审阅未重新组装、求解或求特征值。**

能量语义检查通过逐字段辨认：`Ubulk=2.80764814e-11 J`、`UΓ=8.42708068e-13 J`、`qc·fext=5.78383789e-11 J`、`g·λ=1.08415e-26 J`，保存的能量恒等式相对误差 `3.71e-16`；这些与 `2(Ubulk+UΓ)-qc·fext+g·λ=0` 的 KKT 驻值配对相符，**qc·fext 不是载荷路径增量外功**。primitive `external_work_j` 另记 virgin 弹性 cohesive 路径积分 `8.42708068e-13 J`，与该 primitive 的 stored energy 相符，没有冒称全局外功。

outer receipt 报告 combined output `1,472,018 / 12,582,912 B`，pins 前后相等；inner raw stdout 为 316 B、stderr 0 B。结果、receipt、日志与锁均被绑定。无 solver/eigensolve 的独立复算，故不提升已保存 gate 为二次独立数值验证；本审阅也没有运行 saved checker。

物理边界：这是 27 个 source、81 BASE + 12 global front 系数下的固定平面 Mode-I cohesive 弹性加载案例。没有裂面推进、混合模态、过程区、真实 (K_0/t_{peak}/G_c) 标定或陶瓷资格；`ACCEPTED_DIAGNOSTIC_ONLY` 只指声明的代数和接口门。
