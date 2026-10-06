# Saved-array J 独立复核

按stage现成8条contour的坐标、法向、权重和梯度，每档224点依给定式以math.fsum重算；两档共448样点。样点字段与保存数组逐项一致；54项IDs和系数binding均与primary逐项相同。

对照解析值0.013 J/m²：每弧12点绝对误差0.085–0.464%，每弧16点0.094–0.460%；逐路径重算与保存J最大差1.73e-18 J/m²。数值和8路径图见附件。

原coupling receipt保持UNKNOWN_WORKER_RESULT_CONTRACT，原件未改并逐字节副本归档。此前审阅器两次JSON解析失败，2.4737s未产数值；这是解析失败，不是物理失败。细档控制点|grad_z|max=2.24750839619e-07非零，因此这里只报告midsection functional，不据此称全域二维或实际G。
