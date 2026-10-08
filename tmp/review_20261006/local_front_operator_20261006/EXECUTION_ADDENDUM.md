# One-shot execution addendum

单次授权 VF `-B` 运行针对固定 256 查询、140 DOF 标量算子，外部 rc=125，耗时 0.765 s。未晋级因父校验器期待大写 basis 标签，而 registry 为小写 wire 值，且含已知未用项 `front_f1_unshifted: 0`；这是映射校验失败，不是 MLS/谱失败。staging 与 UNKNOWN receipt 原字节保留，12 pins 前后匹配；未求解、加载、改源或重跑。

[basis_count_mapping_fix.py](basis_count_mapping_fix.py) 只核对已存 staging SHA、140 个 ordered keys 和精确映射 `base:100 / front_f1_shifted:30 / heaviside:10 / front_f1_unshifted:0`；未知或额外 kind 拒绝，不调用 provider/operator。`max_gradient_error_per_m` 是重建 Uz=x 的梯度误差，单位 m/m（无量纲）；基函数梯度 dN/dx 单位 1/m。
