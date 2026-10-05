# PG-Frag 原生自由演化断裂小试：源码审阅

**冻结状态：** `FROZEN_SOURCE_ONLY`。fixed case 的 prelaunch `run_status` 保留 `SOURCE_ONLY_NOT_RUN`；本次没有运行 Julia、测试或 CAD。以后若执行，结果中的 `source_freeze_status=SOURCE_PINNED`、`source_execution_status=EXECUTED` 只说明源已钉住且 driver 实际执行过；结果摘要是 provisional，最终过程和预算接受以 launcher receipt 的 exit、`cap_exceeded`、`native_process_killed`、`final_budget_ok` 与 `result_authoritative` 为准。

## 边界与 intact 数值路径门

小试采用完整解析 uniform box：12 × 9 × 9 mm，`h=0.6 mm`，`δ=1.8 mm`，20 × 15 × 15 = 4,500 点，点体积 `h³`，总质量 0.0020412 kg。没有点删除 mask、缺口、预裂面、手动断 bond、全点仿射运动或固定碎片数。几何缺口暂缓：native horizon 仍可能跨越被删点形成跨空域 bond，只有 immutable point/bond policy 实际 readback 后才能评估缺口。

完整 intact coupon 的独立报告与固定 case 被作为两项只读 baseline pin：报告 SHA-256 `96eb441f9afb7b07d5cc672c4be4681a7b9bf587f7a4ad291785e196e15178ec`；case SHA-256 `e6aa57f3edf136970e9f7108a1df2535a1e7c0d4cb56ee89fee8d858cd3c603d`。报告的 20 项独立算术检查覆盖 native state、grip path、boundary impulse、momentum 和输出完整性。Root 将此结果接受为 `numerical_loading_gate=INTACT_NUMERICAL_PATH_REVIEWED`，仅表示工程状态/边界路径算术已被审阅；`physical_loading_validity=UNKNOWN_UNCALIBRATED`。它没有证明科学加载有效、能量闭合、材料参数真实或断裂物理合格。报告写明能量 residual 仍未资格化、dwell 未达静止，且审阅没有使用 intact `summary.toml` 作独立验算。

pilot 固定使用 root 审阅的 intact native raw stability bound `3.4242302121147243e-7 s`，以及 `dt_fraction=0.10951366490292391`，目标 `dt=3.75e-8 s`。driver 在 native 稳定界与 baseline 不匹配、固定乘积不为 `3.75e-8 s`、VelocityVerlet 接受的 `Δt` 被改动，或相位分辨率不满足时直接阻塞；没有自动裁剪、自动重试或按断裂结果调 dt。此前固定物理时间程保持不变：ramp 1.0 μs、hold 0.6 μs、unload 1.0 μs、fixed-grip dwell 1.0 μs；预注册驱动幅值固定为 native `εc` 的 1.5 倍。最多 96 步，实际是否产生断键由 native 状态决定，不通过后处理补裂。

## 原生失效规则与材料输入

冻结的 Peridynamics.jl OSB `CriticalStretch` 对每条 `fail_permit=true` 的有向 bond 计算 `ε=(l−L)/L`，且只在严格 `ε > εc` 时将该方向 inactive。native 机制不会把 inactive bond 重新激活；pilot 也逐状态检测重激活并失败。参考图须为 reciprocal、全 active；每个时刻检查 directed active 数是 active canonical pair 数的两倍、每点 active 邻居与 canonical 图吻合、累计事件数等于 broken-pair 状态。

只提供 `Gc` 时，冻结的 native OSB 参数转换为 `εc = sqrt(5 Gc / (9 K δ))`，`K = E/[3(1−2ν)]`。当前诊断值 `E=5 GPa`、`ν=0.24`、`ρ=2100 kg/m³`、`Gc=45 J/m²` 得 `K≈3.205128205 GPa`、`εc≈0.002081666`，目标应变约 0.0031225，对应 11.4 mm grip span 上约 35.60 μm。它们是未校准原型数值；幅值只为固定阈值机制诊断预注册，不代表测得强度。`Gc` 只进入 native 阈值，不能由此结果声称已测得断裂功或耗散。

native bond force path 每步依次运行 failure、damage 与 force 更新；OSB damage 是邻域 inactive-bond 比例，不是裂纹面积或耗散能。关键冻结源码：`physics/fracture.jl`（`Gc`→`εc`）、`physics/ordinary_state_based.jl`（strict threshold 与 OSB 能量）、`discretization/bond_system.jl`（逐步失效顺序）、`time_solvers/velocity_verlet.jl` 和 `physics/boundary_conditions.jl`（时间层及 BC）。这些路径随 common-source pins 一并冻结。

## 驱动与动力学含义

driver 只用 native `velocity_bc!` 驱动 x 两端 grips：左端 xyz 固定，右端 yz 固定、x 方向沿连续物理时间 cubic-smoothstep ramp/hold/unload；随后 grip 位置保持，内部点继续由 native VelocityVerlet 自由更新。固定输入估算的纵波速约 1,675 m/s、12 mm 穿越时间约 7.16 μs；1.0 μs ramp 约为其 0.14 倍，因此响应具有明显惯性，明确不是准静态。

VelocityVerlet 的 native 边界顺序需按时间层解释：约束把 `velocity_half` 覆写为规定值后发生 drift，力更新再生成存储的 full-step velocity。因而 dwell 中 grip 位置保持、受约束 `velocity_half=0`，但记录的 grip full-step velocity 仍是 raw post-force 值，可能非零。energy CSV 将自由 DOF full-step 动能与受约束 DOF raw full-step 动能分开列出；两者仍是 solver 诊断，不称为物理总动能。native BC impulse 与 half-step kinetic jump 也仅是数值账本。

## 自由表面、事件、质量及图标签

固定 OSB / `NoCorrection`，自由表面截断不会被修正。按球形 horizon containment，只有 1,134/4,500 = 25.2% 点有完整球形 support；其余点的 weighted volume、dilatation 与力响应可能受边界 support 截断影响，进而改变伸长、起裂和裂纹扩展。`δ/h=3` 不构成收敛证明。本实验没有内缺口，避免将几何排除误当作物理断裂，也不宣称表面误差、裂面或路径偏差已受控。

结果保存完整点状态、native strain-energy density、damage、精确 per-point active-bond counts；另保存 reference canonical bond roster、仅由 native reciprocal active→inactive 转换产生的事件、active-bond graph 连通分量标签和质量统计。事件 ceiling 由实际 canonical pair 数和剩余固定输出预算共同决定，始终不高于 canonical pair 数；检测到超过 ceiling 或任何状态计数矛盾即失败。零事件会如实报告 `NO_BOND_FAILURE`。连通分量数只是图诊断，不代表可见物理碎片、共同裂面或目标碎片数量。

`energy_reaction_history.csv` 包含 native strain energy、raw full-step kinetic energy（按自由/约束 DOF 分项）、BC impulse、半步约束动能跳变、动量更新 residual 和 surrogate actuator-work。断键改变 OSB active support，且 VelocityVerlet 状态时间交错；能量账本不构成 closure，下降量不能解释为 `Gc` 耗散。未计算 fracture-dissipation history。

## 资源估算与冻结范围

case 限制 5,000 点、180 步、190 MiB 结果目录、200 MiB 总 attempt、180 秒。固定 `dt=3.75e-8 s` 对 3.6 μs 物理时程为 96 步、97 个 state。按 21-field Float64 state 与 Int32 graph labels 估算二进制分别为 73,332,000 与 1,746,000 bytes；按 204,778 canonical pairs，96-byte canonical 行上限估计 19,658,688 bytes；16 MiB 预留其他固定输出；256-byte/event 行上限且 event ceiling 为 204,778 时最多估计 52,423,168 bytes。保守结果输出估计合计 163,937,072 bytes，低于 199,229,440-byte 结果目录预算；每个实际文件写入后仍会再做总字节门，launcher 对 attempt 目录及最终 receipt 也做最终门控。摘要本身只 provisional，不可覆盖最终 receipt 的超限状态。

source freeze 由已接受 intact freeze 的 78 个 common runtime/backend/vendor/dependency-helper pins、4 个 pilot 专属 driver/case/launcher/review pins，以及 2 个 intact baseline report/case pins 组成；78 个 common pin 路径与其已接受 SHA-256 保持一致，bytes 与 SHA-256 均按当前文件重新核对。冻结文件不 pin 自身。原 intact coupon、baseline 审阅文件、vendor、production、Project、Manifest、Git 与 DOCX 均保持只读。此 source preparation 不包含运行结果、真实材料资格、准静态资格、能量耗散资格、网格/方向/dt 收敛或物理碎片结论。
