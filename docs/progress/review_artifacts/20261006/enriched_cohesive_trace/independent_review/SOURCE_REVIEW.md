# Enriched cohesive trace：source-only 独立审阅

判断：源码 `READY_FOR_NARROW_TOY_MICROCHECK`；没有执行结果，本模块仍为 `SOURCE_ONLY_NOT_COUPLED`，不构成物理或真实谱资格。冻结读取：主源 315 行，SHA-256 `13050baf07a896396d1e6fc809a4eb4a309865900f4a30c42fa23e2bf7184b29`；future checker 154 行，SHA-256 `1bbae06cdd1c74e49e24680024c1113e9934c03e1d4645f0af09baccc75a2142`。只读源码、README、契约及冻结节点导数；未执行 checker/场/assembler/solver。

主公式静态一致：`J=H+−H−`、`δ=n·Jq`、`g=Jᵀn`，故 `f=Agt`、`K=A k_alg ggᵀ`；冻结导数复用 history 状态及单侧切线。单侧面积只计一次，重复点/历史键/系数序拒绝；共点、单位法向、正面积、两侧刚体候选迹重现及 `J` 消去均有检查。16 QP/93系数上限早于稠密切线分配；软化负切线不被截正，压缩与切向/混合响应明确拒绝。重复试算不提交状态；拐点返回 nondifferentiable 与所选单侧branch。

最终 checker 源为24维toy夹具：非原点检验旋转迹，bad-rigid 注入检验validator；三条front系数在弹性、软化、卸载各作列有限差分，并对front子块完整 `A k aaᵀ`（含非对角交叉项）作独立比较。6个候选仅占BASE列、BASE两侧迹相同，front力与rigid候选点积为零符合该fixture语义；不把它外推为一般93维谱证明。checker声明21次assembly calls，在README所写24次上限内；声明15秒预算。源码检查不等于执行PASS，当前无据称数值通过的结果。拐点本身不在FD点上测试，断点状态只经静态调用链确认。

六个候选向量由调用方提供；源码不检其线性独立性，也不从真实GLOBALFRONT基底生成物理刚体模态。README已限定刚体证据只覆盖所给候选，`BASE/GLOBALFRONT` 也是标签。未发现可复现公式/接口错误。仍未知：实际93维刚体谱、真实provider/front耦合、mixed-mode物理响应及断裂资格。
