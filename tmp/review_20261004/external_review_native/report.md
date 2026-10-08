# 外部报告原生PG证据独立核验

审阅对象为附件《已粘贴的文本.txt》，SHA-256 为 0b3544d4b2a02d8132660834ea5a3d44c71a2054836d3ffeb8ac29756a771bad。附件中的 #1–41 位于物理行 3–83，#153–192 位于物理行 307–385，#202–210 位于物理行 405–421；下文按报告自己的条目编号处理。

## 给 root 的关键结论

报告把被捕获的包预编译 workload/headless 文件输出异常称为“recompile error”，并据此怀疑零载荷未运行；冻结源码和收据显示该异常发生在包预编译 workload 且被 catch，随后原生 submit→VelocityVerlet 一步与 step0/1 hooks 实际完成，点/键数据门通过。它只证明零载荷接口路径，不证明非零构成行为、断键/裂纹、能量闭合、材料校准、真实 A12 链或训练资格。报告声称缺失的力密度、体积、材料值、dt、SI单位、版本与哈希，在这次原始产物里已有记录；图是否展示这些信息未由附件文字核验。5 mm/3 mm horizon、125 点是明显的边界/分辨率限制，不是 P0 bug。rest/current 直连缓冲区的担忧可由静态源码排除；非零多步快照仍未测试。本报告只供主代理裁决，不作验收。

## 口径与限制

本审阅只读取附件文本、冻结的 native run01 源码、进程收据、CSV、summary，以及已冻结的 bias/bond/energy 源码或其哈希收据。未运行测试、solver、Julia、CAD、Git；没有委派。附件提及的三张图未在本轮作为图像输入，因此凡陈述“图中没显示”的内容，只有在 raw source/artifact 有反证时才评价“数据层是否存在”，不推断图面事实。

类别含义：CONFIRMED_CURRENT_GAP 表示这份实验/产物确实缺少该能力或观察；RISK_NOT_TESTED 表示需要另外的输入、工况或比较试验；EXISTING_RAW_EVIDENCE 表示冻结原始文件已有对应记录；DIAGNOSTIC_SCOPE_NOT_BUG 表示能力边界与本零载荷 smoke probe 一致；INCORRECT_INFERENCE 表示报告把已知异常或证据作了错误推断。外部报告的 P0/P1 标签不作为本审阅严重度。

覆盖与分类计数：本轮覆盖所请求的 90 个条目 ID（#1–41、#153–192、#202–210），归并成 37 组；其余 120 个编号未在本轮核验。37 组按类别计数为 INCORRECT_INFERENCE 11、CONFIRMED_CURRENT_GAP 14、EXISTING_RAW_EVIDENCE 9、RISK_NOT_TESTED 13、DIAGNOSTIC_SCOPE_NOT_BUG 8。类别可重叠，故这些数不相加为 37。外部报告的 210 条断言不等于 210 个已确认 bug。

source pin 核对：native freeze receipt 声明 31 个冻结源文件，另记录 Julia runtime SHA；artifact_hashes.json 列 11 项产物。报告 JSON 中所有 SHA-256 字段均为 64 位小写 hex，未发现格式错误。

## #1–41：PG、单位、离散与求解路径

- **#1–2（附件行3、5）— INCORRECT_INFERENCE。** stderr 确实有“precompilation errored”，但不是 recompile 失败：日志将异常定位为 headless backend 拒绝预编译 workload 的 file-output Job，源码在模块 include 周围 catch 并记录异常。之后 child 与 wrapper 都观察到 exit 0、未超时，stdout 明确报告 125 点、1 步、0 个 transition；execution receipt 的状态是 DATA_GATE_PASS_WITH_NONFATAL_PACKAGE_PRECOMPILE_ERROR_RECORDED。异常应如实保留，但它不与实际 run 完成矛盾，也不构成物理正确的证明。（E1、E5、E6）

- **#3–4（附件行7、9）— INCORRECT_INFERENCE；物理正确性仍为 RISK_NOT_TESTED。** 冻结 runner 创建 Job 并调用 Peridynamics.submit；submit_threads 初始化 data handler 后调用 solve!；VelocityVerlet 的 solve! 调用 reference hook、每步调用 force-density calculation 和结果 hook。运行收据记录实际 Julia child、argv、step0/1 两组 point rows 和完整 bond rows。因此“求解器未实际参与/返回初始态所以无法区分”与当前代码路径和产物不符。零载荷只能证明此具体原生接口路径走完，不能证明非零响应、材料定律或求解数值正确。（E1、E4、E5）

- **#5–8（附件行11、13、15、17）— CONFIRMED_CURRENT_GAP；DIAGNOSTIC_SCOPE_NOT_BUG。** 这次只有 zero_load、1 native step；最大速度与记录的力密度均为 0，damage 没变化，directed/canonical transition 为 0。完整图状态仍全部保留；summary 把物理 fragment count 标为 UNKNOWN。故没有观测到断键/裂纹成立，但这不是对断裂实现的反例，也不能把 UNKNOWN 解释成“永远无破裂”。（E2、E3、E4）

- **#9（附件行19）— EXISTING_RAW_EVIDENCE；RISK_NOT_TESTED。** callback 的字段不是由图像推断：runner 从 native ThreadsBodyDataHandler 的 chunk/storage 读取位置、速度、加速度、b_int/b_ext、damage、active bond state，并在原生 reference/step hook 中采样。该证据能建立观测点来自 native state hooks，不能建立 constitutive response 正确。（E1、E5）

- **#10（附件行21）— INCORRECT_INFERENCE；非零快照行为仍为 RISK_NOT_TESTED。** 静态源码显示 observer 的 rest_position 引用 input.position；Body 构造器会保存传入数组，但 native chunk 用 body.position[:, point_ids] 切片，VelocityVerlet 初始化 solver position 时又 copy(system.position)。采样读取 options.rest_position 与 solver storage.position，并把标量值放进新 CSV 行。因此当前 current 与 rest 不是同一变形缓冲区。zero-load 下二者相等也确实不能以图面单独证明这一点；非零多步变形快照仍未测试。（E1、E5）

- **#11–12（附件行23、25）— #11 INCORRECT_INFERENCE；#12 CONFIRMED_CURRENT_GAP 且 DIAGNOSTIC_SCOPE_NOT_BUG。** point_state_history 的物理状态量在两帧不变符合零载荷；但所有 CSV 字段并非完全相同，step 从 0 到 1、time_s 从 0 到 4.925124976599115e-7 s。只运行 1 步无法观察多步漂移或累计积分误差，这正是本 gate 的范围限制。（E2、E3）

- **#13（附件行27）— CONFIRMED_CURRENT_GAP。** 本产物没有完整 work/energy ledger；summary 将 energy_ledger 明确置为 UNKNOWN，缺项没有补零。这是对该 run 的准确限制，不是已检测到能量错误。（E3、E8）

- **#14（附件行29）— EXISTING_RAW_EVIDENCE（限定为力密度）。** point_state_history CSV 含 b_int 与 b_ext 三个分量，列名单位为 N/m³；收据确认已审计所有输出 force-density 值且本零载荷样本均为 0。没有逐 bond force、reaction force 或能量记录；原始 CSV 有力密度不代表图中显示了力，也不把力密度等同于合力。（E1、E3、E4）

- **#15（附件行31）— CONFIRMED_CURRENT_GAP；RISK_NOT_TESTED。** 没有 strain/stress 字段或非零载荷下的 constitutive 结果，因此无法从本 run 验证材料响应。native force-density 字段和模型调用存在，但全为零。（E1、E3、E4）

- **#16–17（附件行33、35）— EXISTING_RAW_EVIDENCE（本次参数传递）；材料校准仍为 CONFIRMED_CURRENT_GAP。** case/config 与 summary 记录 OSB、rho=2100 kg/m³、E=5.0 GPa、nu=0.24、Gc=45 J/m²、horizon=0.003 m；runner 将这些值传给 material!(body; horizon, rho, E, nu, Gc)，并提交 native Job。它证明本次传参路径与参数值，不证明这些值是实测、已校准，或后端材料定律正确。（E1、E2、E3）

- **#18（附件行37）— RISK_NOT_TESTED；DIAGNOSTIC_SCOPE_NOT_BUG。** 输入是 5 mm 实心立方体、1 mm 单元中心点，horizon=3 mm，即 δ/L=0.6、约 3 个点间距；立方体中心到最近表面仅 2.5 mm，所以没有点能拥有完全位于体内的 3 mm 邻域。它确实不适合据此评价 bulk 物理；“边界效应强”是规模/校准限制，不是已观测到的 solver bug。（E2、E3）

- **#19（附件行39）— DIAGNOSTIC_SCOPE_NOT_BUG。** 125 点属于明确声明的接口/事件 smoke fixture，不是空间收敛或数值准确性基准。（E2）

- **#20（附件行41）— EXISTING_RAW_EVIDENCE。** 元数据声明每个 cell-center quadrature point 的体积是 1 mm³，即 1e-9 m³；CSV 每点有 volume_m3，summary 记录总体积 1.25e-7 m³。元数据也明确禁止把它说成 Voronoi 体积或真实壳体体积；图是否显示该字段未知。（E2、E3）

- **#21–27（附件行43–55）— EXISTING_RAW_EVIDENCE（已捕获图状态）。** point CSV 每点有 n_active_bonds，canonical/directed 状态文件完整保留两帧，summary 明确初始 directed=6124、canonical undirected=3062；receipt 记录 12248 条 directed 与 6124 条 canonical sample rows、 reciprocal directions present and equal。observer 对观察到的每对 bond 强制两向互反、状态一致，拒绝 self-bond、duplicate direction，并要求 directed 数量等于 canonical 数量的两倍。故“仅抽查 99/3062”不适用于 raw CSV；“3062 语义不明”也与字段名/summary 不符。范围限制：这些检查验证已观测 native graph 的互反与计数，不独立证明 bond 邻接集合与预期几何/物理域完全无遗漏。（E1、E3、E4）

- **#28–29（附件行57、59）— CONFIRMED_CURRENT_GAP。** raw observer 不输出 bond-length histogram 或 stretch field。坐标和完整键列表使离线重建成为可能，但这次没有相应分布或预应变判据产物；不能由此认定键长错误。（E1、E3）

- **#30（附件行61）— EXISTING_RAW_EVIDENCE（常见坏数 gate）。** 输入读入时显式拒绝 NaN/Inf 坐标或体积及非正体积；execution receipt 对 250 个 point samples 记录所有位置、速度、加速度、力密度、damage、volume、mass、time、dt finite，并记 point audit PASS。代码未声明专门 subnormal 检查；图是否展示坏数面板未知。（E1、E4）

- **#31–33（附件行63、65、67）— INCORRECT_INFERENCE（对原始 run）。** dt 并非不透明硬编码：配置指定 AUTO_NATIVE_STABLE_BOUND 和 safety factor 0.8；runner 先由 Peridynamics.calc_stable_timestep(native handler, safety) 计算界，再要求 submit 后实际 solver.Δt 与该值相符；summary/CSV 记录 dt 与时间。源码明确 VelocityVerlet(steps=1, safety_factor=0.8)，freeze receipt pin 该源码和项目。此证据证明来源/积分器与原生稳定步长调用，不替代物理单位、离散收敛或经验精度验证。（E1、E2、E3、E5、E7）

- **#34–35（附件行69、71）— EXISTING_RAW_EVIDENCE。** zero_load 分支给所有点在三轴设置零速度初始条件，不创建速度边界条件或固定点集合；case 明确 speed=0。它不是“空白配置”，也没有证据显示默认隐藏约束。该设置没有检验有载荷边界条件设计。（E1、E2）

- **#36–37（附件行73、75）— CONFIRMED_CURRENT_GAP；DIAGNOSTIC_SCOPE_NOT_BUG。** 本次点集元数据明确说是内部生成的 cell-center fixture，不是实测试件、导入 shell 或 A12 器物；没有 CAD/tetra/surface 到 particles/bonds 的上游链，因此报告对真实器物链的资格保留成立。（E2、E3）

- **#38–40（附件行77、79、81）— 对原生路径的混单位推断为 INCORRECT_INFERENCE；上游转换仍 RISK_NOT_TESTED。** native point CSV/metadata 坐标单位 m、体积 m³；runner 要求这两种 SI 单位，material 值带 SI 后缀，horizon 以 m 传给 material!，输出含时间 s、力密度 N/m³。图用 mm 轴而另列 m 位移本身不证明转换错；附件没提供图像，也没有 mm CAD 输入进入此次 probe。此次数据路径单位契约有原始证据，但 relief/CAD 到 PG 的跨模块转换没在此 run 中经过。
- **#41（附件行83）— CONFIRMED_CURRENT_GAP。** relief/CAD/PG 的端到端统一单位审计不在该隔离 point-set probe 内，故全链结论未知。（E1、E2、E3）

## #153–192：真实几何链、裂纹与数据资格

- **#153–156（附件行307–313）— CONFIRMED_CURRENT_GAP。** 运行输入不是 Bowl/Beaker/Jar、不是 surface/tetra/shell，也没有 cavity；所以真实 surface→volume→particles/bonds 以及 cavity 粒子化后保空的链条未演示。（E2、E3）
- **#157–162（附件行315–325）— RISK_NOT_TESTED。** 薄壁两侧跨空腔误连、近邻但拓扑不相接的壁/把手/肩部、邻域的可见性或材料域约束，都是 hollow vessel input 才能暴露的问题；此处只有单一连通实心 cube，没有证据表明已出错，也没有测试排除这些风险。（E2、E3、E5）
- **#163–165（附件行327–331）— CONFIRMED_CURRENT_GAP。** native run 对整块 Body 使用一个 OSBMaterial，material! 仅传全局 rho/E/nu/Gc/horizon；未输入分区材料、陶瓷异质性或纹样区域结构弱化字段。（E1、E2）
- **#166–176（附件行333–353）— CONFIRMED_CURRENT_GAP；DIAGNOSTIC_SCOPE_NOT_BUG。** 没有裂纹起始/扩展/分叉、碎片分离、自由表面、碰撞、fragment 分布或损伤能量证据；native 图状态 transition=0，physical_fragment_count 明确 UNKNOWN。Gc 是配置材料参数，不是已测得的 fracture-energy history；CriticalStretch 被调用也不代表 threshold 在零载荷下触发。（E1、E2、E3、E4）
- **#177–182（附件行355–365）— RISK_NOT_TESTED。** 本次只有一个 point layout、一个 dt、一个 horizon，且没有 benchmark。#178 旋转整件几何、载荷和材料场检验 objectivity；它不能替代独立生成/重排点布局来测试方向偏差。#179 是空间点距/离散比较；#180 改 horizon 会改变非局部作用尺度/模型参数，应作为敏感性/校准问题解释，不等同于同一模型的 dx 收敛，也没有普适“裂纹必须完全不变”判据。bias protocol v2 明确把 horizon 与 process-zone 分开比较，未声明其设计检查能证明独立性或移除偏差。（E2、E5、E9）
- **#183–185（附件行367–373）— CONFIRMED_CURRENT_GAP（批量资格）；但每-run QC 已有 EXISTING_RAW_EVIDENCE。** 单个 native run receipt 有 point audit、bond audit 和哈希；它不构成批量训练样本的失败隔离、逐样本 QC 或数据集认证。summary 将 R2/training eligibility 标为 UNKNOWN。（E3、E4）
- **#186（附件行375）— 对这一个 fixture 为 INCORRECT_INFERENCE；真实几何数据 provenance 仍缺。** 冻结配置、生成器 hash、point CSV/metadata hash 和 native 源 manifest 可重建本次输入设置；其来源明确是生成的内部 fixture。它不包含真实 CAD/ornament/试件 provenance。（E1、E2、E4）
- **#187–188（附件行377、379）— 对本次 native input 为 INCORRECT_INFERENCE；跨数据集重复性 RISK_NOT_TESTED。** 使用的 primary uniform layout 是无随机数的确定性 cell-center 构造；同一生成器还定义了另一个 seed=17 的 jitter candidate，但零载荷 run 实际输入是 uniform layout。Julia Manifest 记录 Julia 1.10.12 和 Peridynamics 0.5.4+pgfrag.1；run receipt 固定 runtime executable SHA、31 个源文件 SHA，artifact hash manifest 列 11 项。receipt 未独立探测 Julia binary 版本字符串，但路径、manifest 与二进制 hash 均被记录。重复 run 的字节一致性仍未在本审阅中验证。（E2、E4、E7）
- **#189（附件行379的错误前提，及报告行380）— INCORRECT_INFERENCE。** 非致命 package precompile workload 异常不能推出原生 run 不可复现；真实 child 观察到 exit 0、非 timeout，输入 pre/post SHA 相同且输出 artifacts 被 hash pin。重复执行得到相同结果没有证据，因此只能说 repeat determinism 未测试。（E4、E5）
- **#190–192（附件行381–385）— RISK_NOT_TESTED。** 真实几何/纹样同 seed 重建未在 native fixture 中执行；同一物理 seed fracture 一致性因本次无断裂行为也未测试。uniform points 的确定性生成器并不补足 vessel/ornament/fracture 的重复性证据。（E2、E3）

## #202–210：诊断状态与端到端凭证

- **#202–204（附件行405–409）— RISK_NOT_TESTED / DIAGNOSTIC_SCOPE_NOT_BUG。** 对 Geometry/Ornament 图的面板指标没有相应图像或上游 artifacts 可供本次核验。就 PG native artifacts 而言，数值 gate、单位、数据 counts、unknown scopes 和 hashes 都有 raw receipt/summary；是否把它们画进原图不能由文本附件确定。（E3、E4）
- **#205–206（附件行411–413）— #205 INCORRECT_INFERENCE；#206 DIAGNOSTIC_SCOPE_NOT_BUG。** package precompilation error 与 solver run 并存不代表“假完成”；前者被包模块捕获，实际 run 的 DATA_GATE 状态显式注明该非致命异常。summary 的 COMPLETED 只表示 isolated observer experiment 完成；同一 summary 同时列 calibration UNVALIDATED、energy/contact/physical fragments/R2 UNKNOWN。因此报告关于“完成不等于物理正确”的提醒是对的，若推断这次没有执行 solver 或 summary 宣告了物理资格则不成立。（E3、E4、E5）
- **#207（附件行415）— CONFIRMED_CURRENT_GAP（端到端集成资格）。** 本实验是 isolated native observer；bias 模块只检查 study design/provenance、不调用 solver，energy ledger 只审 supplied histories、不校准材料，bond exporter 明确把 caller-supplied lineage 标为 unverified 且不把 graph components 当 verified fragments。它们目前没有一个共同的端到端资格结果连接 Geometry→Ornament→PD fracture。（E9）
- **#208（附件行417）— CONFIRMED_CURRENT_GAP（跨模块）；#209（附件行419）— 对本次 PG run 为 INCORRECT_INFERENCE；#210（附件行421）— RISK_NOT_TESTED。** 本次已有 run_id、case_id、input_id 与多级 SHA，但它们只标记 isolated fixture/run，无法关联 geometry、ornament、mesh、fracture 的同一个物理样本。hash 存在不等于旧样本失效/重认证策略存在；dataset lifecycle 没在当前证据中验证。（E2、E4、E7、E9）

## 证据索引（实际文件、行与 SHA-256）

- **E0 外部报告**：C:/Users/xiemo/.codex/attachments/91ec634d-1d25-4a0f-b3d6-8541952e5011/已粘贴的文本.txt，SHA 0b3544d4b2a02d8132660834ea5a3d44c71a2054836d3ffeb8ac29756a771bad。目标条目物理行映射见报告开头。
- **E1 frozen observer source**：pgfrag/tmp/review_20261004/pd_point_event_probe/finalfreeze_v3/run_pd_point_event_probe.jl，SHA cb3ba390b3223d588fd203261a363e9c80b6cd6830583780ff195dcacd5be290。输入契约行55–135；native state/bond checks 行138–263；hooks 行266–280；配置、material、dt、submit 行292–385；CSV/summary 行388–430。
- **E2 frozen configs/input contract**：case_zero_load_control.toml 行1–9，SHA a2e6694d047e70ed445ee9037bbcbabbd99095dd060ebe1c78bed674199b2178；common_case.toml 行1–20，SHA ca2e4a64946d32f2ba6abdc1bd5b4d1768f5f7587812f3bcce490c0713a0f3bf；coupon_uniform_voxel_5cube_125.metadata.toml 行1–24，SHA 8f6235661f7ba1ff4a8ac03e08be84284b214ce089118f83381abeb5abda10d0；points CSV SHA 9b9207727364c667ecc0b592ccff34405eaab75a0992e057a2d8a7a121b706da；make_point_inputs.ps1 行5–20、31–75，SHA f218b55964634631fccdaf33ca24649c5c64bc0c2c89c38e7800942fa531b0a0。
- **E3 native raw outputs**：summary.toml 行1–39，SHA e69c85dd07a6e42da4efc7f48520d81e57ea79524e87919b7318e683488ae9c7；point_state_history.csv 行1–251，SHA 9726372d415bbf8201b3dcde659dcacdff06b43777cee20705ec03019fbe6175；directed_bond_state_history.csv SHA 416c6dfb33e500c198612f2c793250fa5eebb9d607bc09d1a3b6d229f99c2845；canonical_bond_state_history.csv SHA 8596a6739608c2ba06d5d8d3b92997c317a45a7056c9b9cb5a31018e9ab5dba3；directed_bond_transitions.csv 仅 header 行1、SHA 2bf1e245bd8d29e1fb69c631a7a5101943690176eccaa6ec20c84999bf983137；canonical_bond_transitions.csv 仅 header 行1、SHA 9add03bb7e2c69d5e20543e6fedc0b385f240b20f43445abacc9769032589692。
- **E4 receipts**：execution_receipt.json SHA cf16573ef82dc94b70de47c1c0919f094808ded2213bae8f28d5c1223fef471b，关键字段 JSON pointers $.status、$.julia_child、$.point_state_audit、$.bond_state_audit、$.freeze，含 31 source_file_count 与输入 pre/post hash；launch_receipt.json SHA dcdfed856d6d78eb8d0a97be44bac504a275832505836af8e2147e35607c5bfa；stdout.log 行1 SHA 8d48230adaa75446bc8069848e0e00eb49deff8750efc62ee0a4ad0657b9fff8；stderr.log 行1–5 SHA 212e034693d4035986be206b89e1c89de402bc2b8ee0f0a0ce326451d0853bee；artifact_hashes.json SHA 19f8533c414152c1318b516bfd7535a95ca1627137d6492b08f7f24eef3d0299，列出 11 个 artifact hashes。
- **E5 frozen native Julia source**：Peridynamics.jl 行157–161 catch，SHA 195127ad9eaaafea064c9b7087c83d7b80b3330b3a46f62b38e3838b00b1f4cb；precompile_workload.jl 行1–71（带 path 的 workload submit 行64–70），SHA c1caeee907c1e9279bc32721c4e5a11e4f20e655413acdaeb5175938acd29dcb；submit.jl 行14–26、99–114，SHA c23e6ae68fa8331023325723457a88945242d49946d39e9ad04508e0dff68ffa；velocity_verlet.jl 行165–195、297–303，SHA 584f0afbb1d0cd9808bfeb56d128327d28f7eca486666ffe2a7bbc6feeea48d6。
- **E6 buffer relationship**：body.jl 行101–124，SHA 60762a4f52c9b6c0b8f23d295b266c15bc697b1e2a356ddd759a561236bd3421；bond_system.jl 行178–182，SHA c98bed960c6b6789f015f02e400d6b3a99e09df22c529f386c384a1564bb69aa；配合 E1 行164–180、266–280 和 E5 VV position copy 行297–303。
- **E7 native runtime pins**：Project.toml SHA `eebff26de1f99a56b29a43ad767d2877bba5b6a0e4a03de1f909eab00f59f4ae`；Manifest.toml 行3–5、417–421，SHA `eb7d8902b84cfc64faaa0d45b0e6a3ae8534a252d06399f3db074da5a738e594`，记 Julia 1.10.12 与 Peridynamics 0.5.4+pgfrag.1；当前 runtime executable `third_party/pgfrag-runtime/extract/julia-1.10.12/bin/julia.exe` 的 GetFileHash 为 `29a5cb5bde6a3f4a151f6e714591ca9344638303d9ed2a893ab31ada3de5d103`，与 frozen execution receipt 的 runtime_sha256 一致。
- **E8 frozen energy/bond sources 与后续收据分开记载**：energy_ledger.py frozen v2 行1–4、28–31、846–850，SHA 6d7cea82c0e1bd105e143e3a3f5324595d2594c99201d9de7a0003d28496ad0e；freeze_v2/integrity.json 的 FROZEN_PRE_REVIEW_NOT_TESTED 是冻结当时状态，SHA 24a65fdfbfdf73185907f74aeecdc7b791aff6dd9e51f3ae821f866913f230f5；后续 v2_run/receipt.json 报告 14 tests PASS、exit 0，SHA a9d22f28044cb834fc37630957295048f949e25c3dcfed4462b021fba8ce0e19；captured stdout SHA 874638d648852fb4725495a6790b404d715c1088f6ddb8a91c2af9f8782f3302 报告 14 passed。bond_events.py.v2-frozen 行47–49、144–151、211–219，SHA 48f2c8785527a84e5237270306fe68427be37f1e643780dd3331da51099a05db；freeze_manifest_v2.json 的 NOT_RUN_WAITING_FOR_ROOT_REVIEW 是冻结当时状态，SHA 28415ad2d7e487324df3eeb665a41077a6c63298969f1f2b0c5cae282539ef34；后续 v2_pytest_receipt.json 与 stdout 报告 19 tests PASS、exit 0，且明确 adapter integration NOT_CONNECTED_TO_PGFRAG_RUNNER，receipt SHA 88a8783379122d6798161212aaf2eb61998d698d347e021fc58ed268b2cac41b；stdout SHA b88b69092c878a260d0a624793773ab4a11b90f4b9873e8a443b6043c414dd39 报告 19 passed；JUnit SHA 1ccc3bb73609c7666c76b2549d8067563521058e50027efe6a25d2842dc022c6 记录 19 tests、无 failures/errors/skips。测试状态不代表物理资格或原生链集成。
- **E9 frozen bias source/receipt**：pgfrag/src/pgfrag/bias_protocol.py 行1–5、921–923、1060–1075，SHA 95dd958818cfdf44cc40f3c61b2a6eb60a68c10f99121d749876e6c1a4f26779，与 bias_protocol/v2_freeze.json 的登记 hash 相同；代码声明只审研究设计/provenance、不运行力学 solver，objectivity 可是 DESIGNED_NOT_RUN，directional bias / point-cloud independence 未建立。v2 test receipt 为 11 tests pass、solver_invoked=false，SHA 266eaeb93c7c2778ec9a30bccf68e02cda92dca20719f8b5ed445fd14009e85e。该模块不认证当前 native run 的材料或数据资格。

## root-only note：非外部报告条目

loaded-gate dt 是 root-only 事项，不计入外部报告 #1–210。v1 角点估算错误是已确认、已修订的历史发现。保存的 `decision_v1.md:14`（SHA `d29e2037905709164244dcb9ced0f29acd5d5a0702d5f331a6ce7feb56ad5e60`）原句称“最少邻居角点”估算 accepted dt 约 0.915–0.987 μs；这不是全点 native 稳定步长的全局最小值。原生 VelocityVerlet 对每点算 `dt_i` 后取全点最小，再乘 safety（`velocity_verlet.jl:136-162`，SHA `584f0afbb1d0cd9808bfeb56d128327d28f7eca486666ffe2a7bbc6feeea48d6`）。保存的 `decision_v2_pre_wavecorrection.md:3,11`（SHA `b074ebafb34783a2056ec2bdc6e24483b35af6145b6e7fb9d9948f92c898da96`）明确称 v1 估算错误、320 步预测高估约 2 倍，改用 run01 实测 dt 作预估并要求复核 loaded run 的实际 accepted dt；该文同时说明新增实验尚未运行。`decision_history.diff.md:7`（SHA `4216eb0de0b072eec7dbe46a7cc77e725d6be00ceb0511f3978d20c64a7b1abd`）也记录以 run01 native accepted dt 替代角点估算。故历史错误已在设计中修订；新 loaded 方案的实际 accepted dt 仍为 NOT_RUN。zero-load 实测 dt 为 4.925124976599115e-7 s（`summary.toml:12,27`，SHA `e69c85dd07a6e42da4efc7f48520d81e57ea79524e87919b7318e683488ae9c7`）。前次审阅记录的 live `decision.md` 版本 SHA `2de9ac4fdbb5a2cd73b7e48cc6c7d3d4aa302c46be85c1c6f65cf5701bacf4b0` 仅保留为当时观察值；该文件此后继续变化，该值不是当前 hash，也不是历史结论的唯一来源。

另核对 cavity design 引用：本 external-review 报告没有引用 cavity design SHA。实际 `cavity_bond_safety/design.md` 的 GetFileHash 为 `6035fa7c3fdead24a43a61ff653b03f55f002326ddebe95b5228db0d0b245e0b`；但 `bond_visibility/design_and_freeze.md:7` 与 `independent_review.md/.json:10` 引用文本 `6035dfa7c3fdeade43a61ff653b03f55f002326ddebe95b5228db0d0b245e0b`，它只有 63 位十六进制字符，并且不同于实际 64 位 SHA-256。这是实际存在的引用差异。被核对文档自身 SHA：design_and_freeze `f95c86eada924153ae0c0467e928d765d4d23750dc6abc5759215e5a3a3d227e`；independent_review.md `270f623bef711e94f21251ce047cd651521ce309c3a87d38ed16257095ce6f58`；independent_review.json `ab172121cba83a0d637031d67c4d7207c82f4d0ccc437649b90e5a7889d3de6f`。未修改这些被审阅文件。

## 尚未知

- 三张原始诊断图是否把 CSV/config 中现有指标实际展示出来。
- 在非零多步载荷下的独立状态快照、非零 constitutive response、真实损伤/断键/裂纹与自由表面。
- 真实 A12 CAD/surface/tetra 到颗粒/键的完整 lineage、空腔与薄壁邻接质量、接触、裂纹表面及 material-region 输入。
- 空间点距/独立生成布局方向测试、horizon/过程区敏感性、dt 敏感性、解析/实验 benchmark、repeat-run deterministic output。
- 批量数据集失败隔离、逐样本 QC、全链 artifact ID 与样本版本失效/重认证策略。

这些是独立审阅的证据边界；最终接受与否由 root 决定。

