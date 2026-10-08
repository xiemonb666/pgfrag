# Mixed enriched operator：独立 source-only 审阅

## READY（静态源审）

审阅对象：`mixed_enriched_operator.py` SHA-256 `f28f3b55e1ee15212c72ffba75bfe5a4b4f324a7fd4921d63dd47cbd840d7e6f`；checker 实际位于 `pgfrag/tmp/review_20261006/mixed_enriched_operator_20261006/operator_arithmetic_check.py`，SHA-256 `9a12def77c85e28b7b642dbf431630422b092414b8cabe77c208023307648286`。

静态映射一致：源ID规范化排序；BASE按 source/component 排列，93布局末尾只放一次全局 F1…F4×三分量，27源为93；诊断布局保留每源分支，为405。adapter逐项核对旧ID/顺序、scope/owner/component/index、B列与梯度及分支乘积法则；operator入口再次核对布局ID、各字段和同一标量场跨分量一致。输入source顺序改变时，位置按ID配对、输出按排序ID生成。

6分量B符合工程剪切顺序；各向同性C完整包含法向Lamé块及三个mu剪切项。原始K按正体积权重累加B_i^T C B_j，f按正面积权重累加形函数乘给定牵引分量；声明单位N/m、N一致。刚体候选只占BASE槽，平移常数、转动为 omega×(x-center)；405诊断额外模式系数对应两条 x/y 分支恒等式，每分量两条。它们是否精确零模仍以形函数的线性再现为前提。

27源上限及每个chunk最多256行、最多256个chunk，在K的n²数组分配前检查；不见ridge、pinv、秩截断、provider或solver调用。模块只依赖Python标准库。

## 阻塞 / UNKNOWN

若目标是**独立核验原始K所有矩阵项**，checker不足：它独立重建B、分支值/梯度、仿射场和牵引力，但K本身由被审operator组装，再只检查K乘已知零模；未独立重算参考K。其`93`检查是DOF/ID数量与唯一性，明确不是93谱或物理资格证明。checker文件头要求另行GO，本轮未执行任何K、field、checker、provider或solver。

因此源级判断为 **READY：bounded raw operator与93/405语义映射可进入获准的算术检查**；**UNKNOWN：93矩阵完整谱/秩、实际MLS额外依赖、具体积分与边界条件下的物理可解性**。本轮无可静态复现的字段错配或公式错误；若要把独立性提升到K项级，需另备独立装配参考。
