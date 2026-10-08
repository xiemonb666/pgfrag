# Fixed-load 前沿重分析 source-only 审阅

**判定：核心诊断路径 `READY_FOR_DIAGNOSTIC_AFTER_PIN`；目前不得执行。** 仅静态读了 module、case、README，并核对 SHA/存在性；没有导入实现、调用 provider、装配、求解、CAD、测试或 Git。

BASE φ/grad 用显式 `(source_id,basis_kind)` 映射提取，目标 registry 的 162 个语义 slot 再映射到新 basis 行；`BASE` 值和梯度原样保留。`front_f1_unshifted` 以 `φ·F̄` 替换，梯度含 `F̄∇φ+φ∇F̄` 的完整乘积项；tip/负 x 分支冻结采样点为 `UNKNOWN` 并阻断。哈希纳入算子 SHA、位移、canonical ID 和新值/梯度字节。

固定 27 个 rest sites、162 DOF、384 原 quadrature（256 volume/128 boundary）；新查询为 0。经新基底重新积分原 128 行 traction；旧载荷只用于回构误差校验，未更新 KIII traction、投影载荷或迁移旧系数。模块把新值/梯度交给 full-3D 各向同性 assembler，源码声称 `λ+μ BᵀDB`；本窄审未读取依赖实现，故只确认接线，未独立审证依赖内部 B/D。

源中构造六个体积物理矩并以同一体积质心归一转动矩；g 固定取 baseline 六分量，且要求 `g[3]>0`。先门控缩放 K 为 6 个近零、156 个正特征值、无负值，再构造唯一一次对角缩放 168×168 直接 KKT solve；无伪逆、ridge、载荷修复或系数 pin。`Π=U−f·d`，`−ΔΠ/(LΔa)` 用单侧面积；`L=.01 m, Δa=.0005 m`，故 `ΔA=5e-6 m²`，双面总面积另报 `1e-5 m²`。不比较 Gc、不提交增长，也不宣称自主碎裂。

算子结果实算 SHA 与指定 `3b31e7…f89cd2` 相同，source-pins SHA 也匹配。baseline JSON 当前不存在且 case SHA 为 null；执行会在 pin 校验阶段阻断，早于 NumPy 导入。**具体预算守卫缺口：** 60 秒计时从 `_run` 开始，未覆盖 source 校验、JSON 读取及数值模块导入；另外改用另一个允许目录内 `--output` 可再次执行，因此 `retry:false` 没有跨调用强制。若预算是整次唯一尝试，执行前应补守卫或由主审明确限定单次调用。无此点外的前沿公式接线 BUG。
