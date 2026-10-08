# Bond-event export 独立数据门审阅

## 关键结论

冻结清单对应的生产文件、测试文件和手工 raw fixture 均与工作树内容一致；保存的 pytest 记录为 `15 passed`，本次未重跑。状态验证逻辑对 bool/数值类型、互反方向、无 healing、step 0 初始断键与后续事件、孤点质量和输入只读处理均有明确控制。输出把 lineage 标成调用方提供且未验证；fixture 也明确是手工数据，因此它只能证明接口控制，不证明真实 solver 字段来源。发现一个可复现的 P2：常量 dt 检查的容差随绝对时间戳 ULP 增长，在时间分辨率不足时仍可能输出 `KNOWN`。没有发现 P1。

## 证据与范围

- 比对 `freeze_manifest.json` 中的 SHA-256 与当前 `src/pgfrag/bond_events.py`、`tests/test_bond_events.py`、`native_shaped_fixture.json`；三者均与对应 `.frozen` 文件一致。
- `pytest_run.txt` 保存的是 `tests/test_bond_events.py` 的退出码 0、`15 passed in 0.06s`。本次只读该记录，没有运行 pytest、solver 或生成流程。
- raw fixture 的 `fixture_kind` 明确写为 “hand-authored native-shaped data … not solver output”；测试 provenance 同样填了 “hand-authored … not solver output”。
- 没有读取仍动态的 `energy_ledger.py` / `bias_protocol.py`，也没有审阅 A12 完整链 audit 模块。

## 发现

### P2：常量 dt 时间比较会在时间戳精度不足时错误声明 `KNOWN`

位置：`src/pgfrag/bond_events.py:525`、`:555–558`。`_times_close` 的容差取 `32 * ulp(max(|actual|, |expected|, |scale|))`；在常量 dt 路径，`actual` 和 `expected` 是绝对时间戳，因此较大的原点可以让容差远大于 `expected_dt_s`。可复现方式：以 `tmp/review_20261004/bond_events/native_shaped_fixture.json` 为基准，保留其余输入和 directed-bond 状态，将 `expected_dt_s` 覆盖为 `1e-16`，并将各帧 `time_s` 覆盖为下面的序列；step 2 原有的 active→inactive 转换仍会输出断键事件：

```text
expected_final_step = 2
expected_dt_s = 1e-16
frames: (step, time_s) =
  (0, 1.0)
  (1, 1.0000000000000002)
  (2, 1.0000000000000004)
```

在 IEEE double 下，相邻记录间隔各为约 `2.22e-16 s`，是声明 dt 的约 2.22 倍；但 `ulp(1.0) ≈ 2.22e-16`，容差约 `7.1e-15 s`，足以让 step 1 和 step 2 的时间误差都通过。因此这组记录会被标为 `CONSTANT_DT_VALIDATED`，从而允许输出断键事件及其时间，尽管给定时间戳不能支持该 dt 声明。建议时间分辨率不足以区分 dt 时将时序标为 `UNKNOWN`（或拒绝该声明）；也可用相邻帧时间差与 dt 比较，容差按间隔尺度限幅。

这是针对“时间戳原点显著大于 dt”的数值边界问题。通常从 0 开始、dt 明显大于时间戳 ULP 的轨迹不触发它。可见代码行没有对应测试；本报告仅依据比较逻辑和 double 精度计算，没有执行复现。

## 按关注项检查

- **Bool / numeric 类型：** point IDs、bond endpoints、step、`expected_final_step` 拒绝 bool；几何量、密度、时间和 dt 通过 `_finite_real` 拒绝 bool、非 `Real` 及非有限值；`active` 要求内置 `bool`。没有发现 bool 以 `1` / `0` 身份混入点 ID 或数值字段的路径。numpy 标量/数组是否可直接接入取决于具体类型与容器协议，当前手工 JSON adapter 已显式转成 Python 整数；尚无真实 adapter 可审。
- **缺帧 / 末端截断：** 多帧轨迹必须从 step 0 连续到 `expected_final_step`；内部缺帧和多帧尾部未到预期末步都会在 `_validate_frames` 抛 `BondEventValidationError`（行 492–504），不会生成伪事件或部分组件。单帧 endpoint 则返回 `UNKNOWN`、无事件并保留该帧组件（测试 131–148 行）。这是明确的严格输入边界；若实际运行可能中途终止且需要保留有效前缀组件，调用方必须知道带着预期末步传入时会整次失败。
- **dt 单位 / 时序：** 名称、文档都声明秒，代码无法识别调用方将 ms/μs 数值错误贴到 `_s` 字段；接入时必须先统一转 SI。常量 dt 允许任意非负首帧时间并校验后续时间；变 dt 路径校验每个间隔与该帧 `step_dt_s`，并要求每个非零帧提供间隔 dt。非递增时间拒绝，缺任一帧时间或 dt 证据则 suppress events 并令 timeline `UNKNOWN`。
- **互反、healing、初始断键与事件：** 每个 initial bond 每帧都要求两方向且 active 值相同；`false→true` 直接拒绝。只有 step 0 的 false 会进 `initial_state_broken_bonds`；事件只从后续相邻完整帧中的 `true→false` 生成。因此未见重复首断事件或 healing 后重报。
- **孤点质量、输入不变：** 组件邻接从全部点初始化，未连接点会作为单点组件；质量按所有点的 `density * volume` 汇总。函数构造本地 dict/tuple 并读取输入，没有对调用方容器的写操作；本项是代码审阅判断，测试未显式快照验证输入不变。
- **UNKNOWN 含义与 lineage 资格：** 时序未知时 `break_events` 保持空并编码为 `UNKNOWN_NO_INFERENCE`，而不是“无破裂”。provenance 明确 `SUPPLIED_NATIVE_STATE_CALLER_PROVIDED_UNVERIFIED`，field 字符串也只作 lineage label。当前 fixture/test 不具备 native 来源证明；生产接入仍须提供真实 solver 帧、映射后的稳定全局点 ID、方向及字段位置，并按声明单位记录时间。该资格限制已在输出显式表达，这里不把它另报为代码缺陷。

## 保存测试覆盖缺口

保存的 15 项测试覆盖事件、组件、初始断键、单 endpoint/missing time/缺 dt、变 dt、互反不一致、重复或缺失方向、healing、内部缺帧/乱序、时间不符、self-bond/重复初始 bond/NaN 及孤点。没有看到专门测试：bool 数值输入、尾部截断多帧、dt 与绝对时间 ULP 的比例边界、输入对象不变性。源码已有 bool 和 healing 等检查；上述时间 ULP 项是本次唯一列出的可复现缺陷。
