# PG-Frag continuation inventory — 2026-10-15

只读盘点；未启动服务、未重跑历史任务、未打开大数组、未改源文件。当前用户授权已覆盖本轮继续工作；10-07 handover 中“需用户重新启动”的旧停工边界不再适用。

## 仓库与运行时

- `pgfrag` branch `backup/forge-work-20261003`，HEAD `21aa69dd131031944a3139b482abf6a0069da023`（2026-10-08）。tracked tree clean；盘点时有 608 条 untracked 路径，多为实验/审阅产物。pgfrag 与父目录未发现额外 `AGENTS.md`。
- 已固定并随项目带有 Peridynamics.jl 0.5.4 源码，见 [`docs/headless_backend.md`](../headless_backend.md)：Julia 1.10.12 可执行文件已解压；锁定环境从 `vendor/Peridynamics.jl` 加载，headless patch 只去掉 VTK I/O。未在本盘点中载入或运行求解器。
- Python 3.12.14 项目环境存在于 `third_party/vesselforge-runtime/venv/Scripts/python.exe`；当前 `python` 命令是 Microsoft Store alias。
- 备选 Akantu 5.0.7.post1 有 CPython 3.12 Linux x86_64 wheel、Linux venv 和部分 cohesive 源码（`third_party/akantu-runtime`）；未接入 PG-Frag，非当前 Windows 原生后端，实际执行环境未核实。相比已有 Peridynamics 路径，需要另建网格、cohesive 插入/本构及结果适配。

## 已有路径和边界

- 生产 job API 在 `src/pgfrag/jobs.py`，求解编排在 `src/pgfrag/backend.py`，Julia runner 在 `src/pgfrag/resources/julia/run_fixture.jl`（仓库根也有 runner 副本）。job 目前只接受 `uniform_box_bar_fixture` 或 `forge_asset_surface`；后者只做资产合同检查，然后 fail-closed 为 `UNSUPPORTED`，因为没有表面到体积材料点离散及空腔/拓扑检查。
- 可求解材料是原生 BB/OSB（所有材料字段均强制标为 `UNVALIDATED_PROTOTYPE`），加载只支持零载荷或对称轴向速度。唯一生产资源档 `local_smoke` 上限为 1,000 点、500 步、单 Julia 线程、120 秒。
- Peridynamics.jl 有 `CriticalStretch` + `material!`、`VelocityVerlet`、`precrack!`，也有多体 short-range-force `contact!`。现有 runner 没有接入预裂纹、器物几何或接触。生产结果仅导出端点点状态、断键对、损伤和完整键图连通分量；不输出逐步事件、裂面、带姿态/接触的刚体碎片。连通分量只表示端点键图，不等同已验证刚性碎片。
- `experiments/free_fracture_pilot_20261005/run_free_fracture_pilot.jl` 是较完整的实验级动态断键/事件记录器，保存 `native_failure_events.csv` 与图状态历史；handover 审阅确认其 6,810 个 first-break 全在右夹持边界并跨过 grip boundary，故没有证明内部裂纹或碎片。不要重跑该 one-shot。
- 其余已有证据也不构成器物破裂：4,500 点 free-evolution coupon 明确禁用断裂、damage 为零；93-DOF cohesive bridge 与 opening-history 是固定界面 toy，opening-history 后期 D=1 但 bulk 仍承载且未断开成碎片。`fracture_artifact_contract_20261006` 有母体/碎片/观察件 fixture，但只是数据合同样例。
- `experiments/forge_volume_ingress_20261005/worker.py` 是未运行、未接生产 backend 的 SOURCE_ONLY STEP→材料点候选；其 256 cells/500 points、180 秒、20 MiB 只适合入口诊断，不是器物网格求解预算。

## 最快可信的最小实现

1. 选一个 VesselForge 生成的空心器型，带封闭实体、内腔、壁厚证据和一套明确的纹样参数；纹样必须进入真实表面/厚度几何或被材料场明确表达，不能只作为渲染贴图。用 `forge.asset.v1` 绑定几何、单位、坐标、材料来源和所有哈希。
2. 在 PG-Frag 包内实现正式的 STEP/闭合三角面到 3D 体积材料点适配器；检查 watertight/orientation、连通实体、内腔、局部壁厚、点体积守恒与最小壁厚采样。不得把边界表面顶点当成实体材料点。然后扩展 `jobs.py`、`backend.py` 和两份 runner，让 runner 接收真实点/体积及可审计加载，而不是固定长方体。
3. 首个力学候选沿用 OSB + CriticalStretch + VelocityVerlet，保留材料/几何输入 pins；定义单一、可复现实验夹具和加载路径。输出每步断键事件和检查点状态，独立重建键图组件；再生成可追溯的碎片点云/表面候选。首版可报告键图分裂，不能声称已解决破裂后自接触、摩擦或碎片碰撞。
4. 将合成样本与实测样本贯通：同一 case 记录 asset/pattern/material-batch/fixture/solver/source hashes、加载和离散参数、结果、观察到的裂纹/碎片标签及其证据来源；模型产出的组件数不能充当训练真值。Akantu 暂列候选，不是此最短路线。

## 建议预算与数值边界

- 已有 500 点体积入口和 1,000 点/500 步 smoke cap 无法代表真实器物。首个**计划性**一器型诊断可先以 1 个器型、1 个纹样、约 10k–20k 点、500–1,000 步、单 CPU 进程、10 分钟墙钟、256 MiB 输出作为硬预算；这只是待测起点，不是性能承诺。实际点数必须先按器壁厚度/视界分辨率核定；若壁厚无法由数层点表示，应停止并加粗模型或缩小器型。全逐点逐步历史可能迅速超过预算，优先保留断键事件、固定检查点和可复验状态哈希。
- 需要后续独立验证：同批陶瓷密度、E、ν、强度、断裂能及其不确定度；材料模型参数反演；时间步减半、点距/视界与方向收敛；自由表面/薄壁偏差；能量与边界功核对；真实裂纹路径、断面与实测碎片对照。现有可逆弹性、静态弹性和 cohesive toy 的通过项不能代替这些门。
- CriticalStretch 是理想脆性、不可逆删键模型；几何离散、视界偏差、加载/时间步、材料校准会影响裂纹路线和损伤。当前 PG 还缺器物级体积离散、断后自接触、裂面历史映射和经实测确认的碎片姿态，因此第一版应称“来源绑定的动态破裂诊断”，不能称合格材料预测或真实撞击碎片仿真。

主要事实来源：[`README.md`](../../README.md)、[`headless_backend.md`](../headless_backend.md)、[`20261006_physics_core_status.md`](../progress/20261006_physics_core_status.md)、根 `docs/handover/20261007/HANDOVER.md` 的 PG-Frag 段，以及上述源文件/实验摘要。
