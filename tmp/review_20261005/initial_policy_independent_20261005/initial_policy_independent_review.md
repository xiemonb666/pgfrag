# InitialBondPolicy 独立静态审阅

结论：冻结补丁的核心 admission seam 与数据重建顺序符合隔离 contract；发现一个可能绕过“CPU-only”的具体入口：单 rank MPI 单体会被误当成 one-chunk 而放行，需由主审判定该边界是否严格禁止。

审阅只读完成。git apply --check 对冻结 patch 返回 0；manifest 所列路径均存在且 SHA-256 匹配。没有应用 patch、运行 Julia 或执行测试。冻结身份：patch dc2803eda64e796901a5ad943b3ecf4e9f25da204f99658fdd1e32d3353bdbbe；helper 80a456f86adf7a02e706a90146928f83af36bad0f66cce60e7f9c9bb2ed7749a；test source 8636f20aa8e125fe3923723e94f31f2693d018f30e1729925a6c9d9aa1fd0889；manifest d96182e68b9da91062069e96282d8abbbf824b55e5e1e79516275154ef6cccce；architecture report 52ba1a835fb81b5ce5f6d6c284bc1a458284ddd08e865d3002ef737d05d669c3。

具体 launch 边界：补丁只在 MPI multibody entry 拒绝带 policy 的 body。未修改的单体入口 src/core/mpi_body_data_handler.jl:81-85 用 mpi_nranks() 建 decomposition；force_mpi_run!() 允许 1 rank。此时 decomposition 是完整单 chunk，helper 的 _validate_initial_policy_scope 在 initial_bond_policy.jl:789-792 会接受它，因而 submit 的 MPI 路径能绕过严格的 Threads/CPU-only 限制。若 contract 要求仅 CPU Threads，应在单体 MPI 入口也拒绝 policy。多 rank 单体会因 pd.n_chunks != 1 fail closed；多体 MPI 的补丁入口也先检查 policy 再抛出未实现错误。

其余关键 mechanics 静态吻合：真实 JSON/TOML 字节分别与外部传入 SHA 比对；CAD、CSV、classifier JSON/source、raw-hit sidecar 均从实际文件双读并与 pinned SHA 比对（initial_bond_policy.jl:694-773）；CSV 的 ID、行序和 Float64 坐标须与 TOML 精确一致，Body 坐标用 1e-12 m、rtol=0 比较（:681-689, :794-819）。完整 native directed raw 搜索在 precrack filter 前逐方向等于所有 TOML canonical pairs，并检查无漏项、额外项、重复项和 reciprocal 缺失；mask 后重建 n_neighbors，拒绝 UNKNOWN、空候选、空 mask 与零度点（:907-955, :781-786）。后续 halo/localization、bond IDs、correction、kernels 和 BodyChunk storage 均消费过滤后的 graph；OSB weighted-volume/dilatation 遍历同一 system.bonds。Geometry exclusion 不更改 failure/damage 或物理事件字段，receipt 保持 NOT_RUN/unauthorized。

依赖与 include 顺序静态通过：patch 的唯一 Project.toml 变化是在 fork overlay deps 增加 stdlib SHA/TOML；helper include 时 AbstractBody 已定义，method signatures 不提前引用后定义的 PointDecomposition 或 Bond。Classifier 的科学 truth/calibration 与 native readback 仍未知，不影响此处隔离 contract 静态接受。
