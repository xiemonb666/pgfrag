# Mixed93 cohesive bridge：独立源码审阅

结论：`SOURCE_READY_WITH_EXTERNAL_ATTEMPT_WATCHDOG_REQUIRED`。在冻结输入下，数值路径的源码接口和方程装配相容；未发现阻止一次受控成功路径的源码问题。此结论只审源码与元数据，不是数值结果、物理资格或执行授权。

冻结身份核验：`run_actual_cohesive_bridge.py` 为 479 行，SHA-256 `36ac1710286b0abf607bda5b1c6ffab893f8191830eae091517e2acd41d62236`；supervisor 为 160 行，SHA-256 `99341a82fa1cbe4192fcac3bf76030e7d094c254606414753c80de6e31a3fba9`。`NUMERIC_SOURCE_PINS.json` SHA `090d6fc5f8c695765a7bb6fe91490d0e7c54aeca337a6f9f318734de4cae3a14`，配套 `NUMERIC_SOURCE_README.md` SHA `124a012d43cd643a7e490ae3d6fb1e2399f17368182988bfbf5fb07297783f76`。case、static checker、原 SOURCE_PINS 分别为 `fcc4b46e…b3f3a71`、`023d7304…b6f1412`、`9cb1fee4…c4554f3`。numeric pins 中 20 个源码/元数据引用全部匹配；接受缓存仅做原始字节 SHA，匹配 `3e1e6fd17751e4e06113010f20ea33c7b69914ca90431e26bbf76124951b5d41`，未解析 NPZ 内容。case 与两套 pins 均保持执行关闭；`numeric_runs` 不存在。

API 静态核对通过。MLS 配置参数与 `MLS3DConfig` 一致；`shape_functions_3d(source, queries, config, visibility_mask=None)` 返回的 `weights/gradients/statuses/support_counts/condition_numbers` 与 worker 使用吻合，16 行在一次批量调用中请求。93 语义序列由 81 个 source BASE 后接 12 个 branch-major global 分量构成；`H+ / H-` 共用每点 BASE 值，F1 两侧为 `±sqrt(-x/R)`，法向选中 `[82,82]` 的唯一 `KG=2K0RL=20000 N/m` 项。trace 构造器、`EnrichedFacePoint.point_id/jump_operator`、rigid 候选和 `EnrichedTrialResult/PointTrialResult` 访问字段与冻结源一致。trace 复用其 pin 过的 derivative 模块及 `_SCALAR` 类对象；worker 取同一 `_DERIVATIVE._SCALAR`，没有重复定义材料/状态类。mode、剪切/压缩拒绝、positive internal force 与物理上表面 `-tn` 符号、单侧面积计数和 history 字段的用法均相符。

成功路径会保存 16 个 MLS 查询行、逐点 H+/H−/J 和刚体证据，以及 baseline/final 16 点 history、primitive 力与切线；3 次 assembly 和 1 次 solve 在调用前递增并另记成功返回数。`q_solution`、载荷、gauge、物理乘子按 0.01 缩放，`q_exact` 仅诊断；K 保持原值。KKT 使用原始 K+正 KG、物理六个 gauge，并恢复物理乘子，没有 ridge/pinv。能量式 `2(Ubulk+UΓ)-q·f+g·λ` 明确标为代数配对；primitive `external_work_j` 另列为 cohesive path integral，未冒称全局增量外功。worker 结果与 receipt 的数值序列化用 `allow_nan=False`。

外层 60 秒门仍为必需：supervisor 从 preflight 计时并把剩余时限给 child，但 `total_seconds_from_preflight` 在 child 返回时取样，未覆盖后续 post-pin 核验、日志及 receipt 序列化。由独立 external watchdog 记录完整外层进程的真实 elapsed 与前后 pins；超过 60 秒即拒收，即使内部诊断为 accepted。64 MiB 是源码中的 working-memory estimate 门，不是实测 process RSS。另 `interface_qps` 仅在完整构造循环后赋值；任一循环异常导致 `NOT_ACCEPTED` 时，0 表示没有完整 ledger，不能解释为零次尝试/零个已构造点。只有完整成功包的 16 行证据可接受。supervisor 的禁用态字符串及 lock 元数据 JSON 未显式 `allow_nan=False`，但仅包含固定字符串、布尔、PID 与 `time.time()`，不是数值结果路径。

未导入 NumPy，未打开/解析 NPZ，未调用 MLS/provider、primitive、solver 或测试。当前仍为 `SOURCE_READY_EXECUTION_DISABLED`；任何实际数值成功及其物理解释均未知。
