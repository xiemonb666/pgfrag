# PG-Frag 薄壁空腔接入前：键可见性设计审查

**结论：whole-vessel PD 接入应被此项阻断，直到静态材料拓扑筛键和相应力学记账通过验收。** 当前 native Peridynamics.jl 只按参考坐标的欧氏距离建立 horizon 邻域；它没有检查两点之间的线段是否留在材料内。125 点实心 coupon 不会覆盖空腔或把手间隙的风险。此审查只读完成；没有改代码或实验输出，也没有运行 Julia、solver、CAD 或测试。正在运行的实心 zero-load 实验可以继续。

## 源码确认

PG-Frag 后端优先从 `pgfrag/src/pgfrag/resources/julia` 载入 Julia 项目（`backend.py:122-134`，启动参数在 `backend.py:563-567`）。该项目 Manifest 将本地 Peridynamics.jl 固定为 `0.5.4+pgfrag.1`（Manifest `:417-421`，包 Project `:1-3`）。下面引用的是这个实际资源副本；关键源文件与仓库根目录的 vendor 副本逐文件 SHA-256 相同。

| 判断 | 代码位置 | SHA-256 |
|---|---|---|
| horizon 键枚举和预裂纹筛键 | `pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/discretization/bond_system.jl:55-164` | `c98bed960c6b6789f015f02e400d6b3a99e09df22c529f386c384a1564bb69aa` |
| 激活的 Julia Manifest（Peridynamics.jl 与 PointNeighbors pin） | `pgfrag/src/pgfrag/resources/julia/Manifest.toml:417-421,434-438` | `eb7d8902b84cfc64faaa0d45b0e6a3ae8534a252d06399f3db074da5a738e594` |
| Peridynamics.jl 导出接口 | `.../vendor/Peridynamics.jl/src/Peridynamics.jl:27-31` | `195127ad9eaaafea064c9b7087c83d7b80b3330b3a46f62b38e3838b00b1f4cb` |
| Body 输入和保存字段 | `.../discretization/body.jl:2-27,59-124,201-207` | `60762a4f52c9b6c0b8f23d295b266c15bc697b1e2a356ddd759a561236bd3421` |
| 预裂纹语义 | `.../discretization/predefined_cracks.jl:12-22,31-60,63-111` | `6b2d67e1af7934016e2f40d0f695f42c3cee899229a3a65d422ec31a9341d9d2` |
| OSB 力、weighted volume 和 dilatation | `.../physics/ordinary_state_based.jl:137-235` | `65ffd614270c8f409a28f83e5def1ea1abc5f94ddcb29eac7b2cd7320e720846` |
| correction 初始化和 bond graph 使用 | `.../discretization/bond_system_corrections.jl:43-73,77-100,140-179` | `63c3a3c458f84dcaa0ad3345c141f67fad2ae42f165b1df3c32ecc6e9c508186` |
| Abaqus Tet4/Hex8 到 material points 的转换 | `.../src/AbaqusMeshConverter/AbaqusMeshConverter.jl:7-90` | `b3d09a00fba94ea3a66e132b2bc61573cafea3fcbdf8d35d35dde1e674c836b1` |
| 当前 solver fixture | `pgfrag/src/pgfrag/resources/julia/run_fixture.jl:64-70,126-143` | `2468586567a2a777de155febd8b1ad9829016f682ede7b69c5f2f7a5c9952534` |
| PG-Frag 资产和后端门禁 | `pgfrag/src/pgfrag/assets.py:145-180,230-232`; `backend.py:461-487` | `assets.py`: `9e56b9f1f2ef289b061cf9a1f4241d473990f895d330e40cfc7b100caeaad889`; `backend.py`: `8b3577201e0e8e94e747b3c5c3529c2d117d8ceccb729f0d2d2510dd3fd834dc` |
| endpoint fragment/event 语义 | `pgfrag/src/pgfrag/backend.py:40-74`; `bond_events.py:211-220` | `backend.py`: `8b3577201e0e8e94e747b3c5c3529c2d117d8ceccb729f0d2d2510dd3fd834dc`; `bond_events.py`: `48f2c8785527a84e5237270306fe68427be37f1e643780dd3331da51099a05db` |
| PointNeighbors 距离筛选 | Manifest `:434-438`; pinned depot `PointNeighbors/STFG7/src/neighborhood_search.jl:390-417` | tree SHA-1 `e1273088c75fe6f79e34dfdd6cb67a1cd565619b`; source file SHA-256 `040525485b7b4d4eceaf57f7e55f8e14c67ced570668989843c3baa4616b813f` |
| OCC clipped-cell 候选接口 | `pgfrag/experiments/cad_clipped_cells_spike.py:850-955` | `e9bd8ef98824efcdd7bf47589196d439a969603ab7da8e46b59d47d648084914` |

`BondSystem` 在 `bond_system.jl:96-123` 用 `GridNeighborhoodSearch` 和每点 horizon `δ` 枚举邻居；PointNeighbors 对位移向量计算平方欧氏距离，并以 `distance² <= search_radius²` 纳入候选（固定版本 `neighborhood_search.jl:398-417`）。Peridynamics 随后只调用 point-set precrack 过滤（`bond_system.jl:134-164`）。这两处没有固体占据、材料连通性或 segment-in-material 判定；导出的 API 列表包含 `precrack!`，未提供可见性/遮挡筛键 API（`Peridynamics.jl:27-31`）。

Native `Body` 接受位置和点体积，并在结构中保存 `position`、`volume`、点集等字段，没有原始 CAD 或体单元面邻接（`body.jl:59-100`）。此外，`Body(mat, inp_file)` 支持 Abaqus Tet4/Hex8：`read_inp` 将每个单元转成中心点及单元体积并保留 element sets（`AbaqusMeshConverter.jl:63-90`），但转换结果不包含原四面体边界/邻接，不能在该点云上可靠遍历材料路径。若采用 tet 路线，筛线段应在转换丢弃 connectivity 之前执行，或让源体网格邻接随点保留下来。

当前 solver fixture 是 `uniform_box` 点集（`run_fixture.jl:64-70`）。PG-Frag 的 forge asset 包含 surface vertices/faces 与外观/材料字段，但声明表面到体点离散及 void topology 尚未实现，后端因此在不启动 solver 的情况下返回 UNSUPPORTED（`assets.py:145-180,230-232`; `backend.py:461-487`）。另有 OCC clipped-cell 候选可输出中心和精确正体积，但源文件明确标为静态几何候选、不是 PG job/solver 输入/键图（`cad_clipped_cells_spike.py:2-5,899-918`）。这些接口表明有可能获得材料体积；它们尚未提供 PD 键的材料路径筛选。

## 最小设计

对参考态中每一条欧氏 horizon 候选边，运行一次静态三态谓词：`IN_MATERIAL`、`OUT_OR_VOID`、`UNKNOWN`。只允许 `IN_MATERIAL` 成为 native bond。`OUT_OR_VOID` 以及 `UNKNOWN` 都不入力学键图；`UNKNOWN` 必须保留单独计数和原因，不能按“可能没问题”强连，也不能悄悄当作 `OUT_OR_VOID` 删除后宣称输入合格。若 UNKNOWN 数量超过预先声明的限值，或影响关键薄壁、连通路径、局部 horizon 支持/图连通，则整个 asset 标为 `UNKNOWN/NOT_ADMITTED` 并停止进入物理求解。该限值和“关键”区域须在实现前定义，不能事后放宽。端点需明确位于同一材料体/连通材料域内部，且整条开线段必须在该材料内。此约束应基于未变形参考坐标计算一次并固定，而非每步随位移重算。

实现时把可见性筛选放在 native candidate graph 建好后、`BondSystem` 完成初始化前：也就是类似 `filter_bonds!` 的图构造阶段，按 canonical pair 一次分类并对两个有向副本对称执行；同时更新 `n_neighbors`，再生成 `bond_ids`、surface correction、kernel 与 storage。当前构造顺序先拿到 bonds/count，再创建 correction/kernel（`bond_system.jl:55-74`）。几何 mask 应从一开始定义模型的初始 bond domain；不要在末态结果上删键。

可靠的线段判定至少应执行以下步骤：

1. 确认两个端点都明确在同一材料体内；边界距离落在源几何误差带内时返回 `UNKNOWN`。
2. 求线段与该材料体全部边界的参数交点，包含内腔边界。合并数值上无法分辨的 hit 时不得把小间隔吞成“连续材料”。
3. 按交点把线段分成开区间，对每段中点做鲁棒体内/体外分类；每一段必须属于同一材料连通域。仅两个端点为 IN 不够。
4. 切线擦边、与面近共面、交点次序不确定、薄层厚度低于误差带、非流形/自交/开口面、float32 坐标导致的占据不确定，一律返回 `UNKNOWN` 并不建键。
5. 按无向 pair 缓存判定和来源 hash/精度/算法版本；输出 `AIR_PATH`、`DISCONNECTED_MATERIAL`、`GEOMETRY_UNKNOWN` 等静态拓扑原因。汇总 UNKNOWN 影响的点、局部邻域支持和关键几何区；触发 asset 级门禁时不启动 solver。不得将其写成已断键事件或静默等同空气。

若模型有多个材料体，应先按 body/material ID 分区；几何相邻不代表材料连续。真实接触面应由 contact 接口处理，不能用跨空气的 bond 替代接触。

## 可选实现路线

| 路线 | 判定方法 | 权衡 |
|---|---|---|
| A. 源 CAD/B-rep 全交点加区间占据分类（推荐） | 对候选线段求与封闭 solid 的全部交点，逐区间用 CAD solid classifier 分类。若源是三角 surface，则先要求封闭、定向一致、无自交，再用 all-hits/winding occupancy；不满足返回 UNKNOWN。 | 对内腔和多次出入最直接，能避免端点测试漏掉中间 void。必须保留 authoritative source 几何与高精度；surface 拓扑坏或擦边情况 fail-closed。当前 asset 合同只有表面网格，所以需要上游提供/验证可判定体域。 |
| B. 原 tet volume 网格的线段穿越 | 用 BVH 找出与线段相交的 tetra，求各自的参数区间，要求同一材料的 tet 内部区间覆盖完整 `[0,1]`；穿过边/面时只在经验证的几何容差内合并。 | 与现有 `read_inp` Tet4 点体积路线相配，体域来源清楚；但 native 转换会丢掉连接关系，必须先筛或保留单元拓扑。结果受体网格边界质量、退化 tet 和网格误差影响。 |
| C. SDF/体素的保守 segment supercover | 遍历线段经过的全部体素，并用有误差界的 signed-distance/Lipschitz 判定每段确实在体内；边界体素或界不充分的区间返回 UNKNOWN。 | 批量查询和加速方便；内存及分辨率成本高，边界误差直接影响薄壁。只在有保守误差界时可采用，简单等距采样会漏过比采样间隔更薄的穿越区域。 |

推荐先把数据契约改成携带一个权威、封闭且可分类的源材料体域，再实现 A；若用户指定的输入本来是体 tet 网格，可用 B，但必须在 `Body` 点云转换前保留 traversal 所需拓扑。当前 surface-only 资产若无法证明闭合占据，不能靠 vertex/face endpoints 推断空腔内外。C 可作为有界近似，不作为默认真值来源。

## 初始键、OSB 与损伤语义

`bond_active` 是 fracture data field；`init_field_system(..., Val(:bond_active))` 默认把所有系统 bond 初始化为 `true`（`bond_system.jl:349-367`）。`precrack!` 适合解释两种初始裂纹状态，但不等于通用材料可见性筛选。`update_dmg=false` 会让 set-cross bonds 在建系统时被过滤并从 `n_neighbors` 移除；默认 `update_dmg=true` 则把它们保留在图中、初始化时置 inactive 并计算初始 damage（`predefined_cracks.jl:12-14,31-60,63-81,106-111`; `bond_system.jl:142-164,228-245,317-367`）。几何拓扑 no-bond 应仿照前一种“在分母和系统图之外”的语义，但要用 segment predicate；point-set precrack 只能表达集合两侧的交叉关系，不能代表任意 pair 几何遮挡。

OSB 每点 weighted volume `m` 的累加含 `bond_active`（`ordinary_state_based.jl:207-219`）；dilatation 的累加也含同一 activity mask，且用同一 `m` 作归一化（`:222-234`）。因此，初始化前屏蔽键后这两项会共同排除它。可是 native damage 是 `1 - n_active_bonds / system.n_neighbors`（`bond_system.jl:317-321`）。如果只在 storage 上将遮挡边的 `bond_active=false`，却保留原 `system.n_neighbors`，会把空气路径误报成 step-zero material damage。`bond_active` 是 constitutive break 状态，不应用作永久拓扑类别。

此外，`EnergySurfaceCorrection` 的 `mfactor/scfactor` 在系统初始化时从 bond graph 和活跃键计算（`bond_system_corrections.jl:43-73,77-100,140-179`）；native 构造 correction 时也已经传入 bond 数量。故筛选必须早于 correction/storage 初始化并让其基于最终 admissible graph 重算。事后删末态键无法撤销此前错误的 OSB force、`m`、dilatation、damage、surface correction 和碎片连通分量。

在报表中把候选 pair 分成“有欧氏 horizon 候选但几何未成键”与“已成键后在 step 0/之后断裂”两种。前者不进入 `bond_states.csv` 的失效集合、`broken_bonds.csv`、断键事件时间线、step-zero damage 或 R2/裂纹统计；为可审计性另存静态 geometry-mask ledger。现有输出将 `bond_states.csv` 定义为 initial native bonds 的终态活动状态（`backend.py:48-64`），fragment 标签也仅是 endpoint intact-bond graph components、不是经验证的刚体碎片（`backend.py:40-46`; `bond_events.py:211-220`），新 ledger 必须避免改变这些字段的含义。

## horizon 与校准边界

该可见性规则是要显式声明和验证的材料相互作用域规则；Peridynamics.jl 默认模型只按距离，不隐含“线段须在材料内”。筛掉跨空腔边会改变 horizon family、OSB 积分支持、表面修正和有效刚度。因此，不用缩小 `δ` 或增大点间距来隐藏假连键；也不据此重选过程区长度或材料参数。保持已选/待校准的物理 horizon `δ`、点距 `h`、断裂/过程区参数独立，报告 `δ/h` 与薄壁局部离散分辨率，并在相同筛选规则下做空间细化、薄壁响应和材料参数校准。没有通用的 wall-thickness/horizon 比率可在这项审查中替代验证。若精度不足以判定薄壁或空腔路径，结果是 UNKNOWN/拒绝该 bond 或阻断该几何输入。

## 最小正负控设计（只设计，未执行）

- **实心 coupon 负控：** 对现有 125 点实心 zero-load 输入，使用同一 horizon 候选和可见性谓词。要求 geometry-blocked 数量为 0、eligible bond pair 集合与原图相同、mask 不引入 damage/破键/R2 事件；已有 zero-load 实验照常继续。
- **薄壁空腔正控：** 建封闭、已分类的空腔壳体，在两侧薄壁内布点并保证一对相对内壁点的端点都明确 IN、欧氏距离小于 `δ`，但连线穿过空腔。要求该 pair 为 `AIR_PATH` / 非 bond；同一壁厚方向和壁面内经实体的有效 segment 仍保留。设置低于几何精度/薄壁分辨率的近边界 pair，要求 UNKNOWN 而非强连；若这些 UNKNOWN 影响关键壁支撑或超过预声明限值，整个 asset 应保持 UNKNOWN/不准入。
- **把手间隙正控：** 准备端点各自位于实体、彼此距离小于 `δ`、中间连线穿越外部空气的把手与壶体邻近段（或开环手柄两臂）；要求无跨隙 bond。另放一对沿真实实体连接/融合部位的控制 segment，要求通过；若 CAD 将二者定义为独立材料体，则接触关系另测，不能强造材料 bond。

验收需同时核对 pair 级筛选 reciprocity、每点 `n_neighbors` 等于 admissible 图 degree、`m`/dilatation 使用同一静态图、correction 在该图上初始化、几何排除 ledger 与动态断键 ledger 无重叠，以及负控初始图完全不变。先通过静态拓扑门，再进入力学比较；本次审查未运行这些检查。
