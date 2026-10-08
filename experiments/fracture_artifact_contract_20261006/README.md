# PG fracture artifact contract (2026-10-06)

本目录定义一个隔离的实验导出格式 `pgfrag.fracture_artifact.v1`，供后续 A12 surface 组件接入。它不调用 solver，不重建裂面，也不改动 `pgfrag` production runner、SynLab PG source reader 或 chain dispatch。

## 安全边界

母体、材料、load/solver/discretization provenance、闭合 fragment mesh 的相对文件路径/hash/mass/volume、双侧 reciprocal fracture surfaces、母坐标/反向单位法线/三角形、事件时间/step/lineage 和 `clean_ground_truth` 都被严格要求显式表达。surface pair 必须引用两个不同的闭合 fragment；patch 两侧 surface ID 必须互为 reciprocal，片段归属要和 pair 一致，母坐标逐点匹配、法线单位长度且反向，两侧三角形绕序相反，声明面积须与母坐标三角面积吻合。Bond/point graph component 只可作诊断，不能替代闭合 mesh partition。

`observation` 和 `clean_ground_truth` 是两个独立通道。`read_observation()` 只返回 `a12.fracture_observation.v1` 的观测片段、sha 已核验的观测网格/mask bytes（base64）和观测处理 provenance；它不返回母体 ID/hash、文件系统路径、事件、surface pair、lineage 或任何 label。观测 payload 当前只接受最小 ASCII PLY（三维顶点/三角面，无 comments 或 metadata）和只含 `visible_vertex_indices` 的 JSON mask，避免在几何 payload 里夹带母体答案字段。这样 A12 可以只用投影中的观测 payload 为观测碎片产生 candidate 与 evidence。`read_clean_ground_truth_for_audit()` 要求明确的 `purpose="offline_audit"`，输出也带 `AUDIT_ONLY_CLEAN_GROUND_TRUTH` 角色。两通道没有自动合并方法。

材料 `UNCALIBRATED`、fixture-only 几何、checkbox 或 caller 自报 PASS 都不会变成预训练资格。`FIXTURE_GEOMETRY_ONLY` 和 `SIMULATION_UNCALIBRATED` 的 `pretraining_qualified` 恒为 false。`PHYSICALLY_VALIDATED` 必须含来源、标准、验证者及签名 receipt，并由调用方传入的可信 `physical_receipt_verifier` 独立验证；本实验没有内置信任根，所以默认验证会拒绝 physical promotion。schema 只允许指定字段，`physical_validation_passed` 这类 checkbox 会因未知字段被拒绝。

## 不确定值

- 有数值（包括真实零）写作 `{"state":"KNOWN","value":0.0,"unit":"N"}`。
- 已知未知值写作 `{"state":"UNKNOWN","reason":"..."}`；不适用则用 `NOT_APPLICABLE` 并说明原因。未知不能带 `value` 或 `unit`。
- 缺必填字段、只有 `KNOWN` 却缺 `value`、单位不符、非有限数和重复 JSON key 都失败关闭。

闭合 fragment 必须有单独的 mesh reference、sha256、参考 frame、SI 单位和 closure method。`fragmentation_basis` 只接受闭合 mesh partition；只给 graph component 会被拒绝。mesh `watertight=true` 仅描述来源提供的 closure 证据；fixture 的 `FIXTURE_ASSERTION_ONLY` 不构成几何或物理 QA。

## 使用

```python
from fracture_artifact_contract import (
    export_artifact,
    read_observation,
    read_clean_ground_truth_for_audit,
    validate_artifact,
)

report = validate_artifact(artifact)
export_artifact(artifact, "fracture_artifact.json")
inference_input = read_observation("fracture_artifact.json")
audit_labels = read_clean_ground_truth_for_audit(
    "fracture_artifact.json", purpose="offline_audit"
)
```

`fracture_artifact.schema.json` 给出 Draft 2020-12 的严格形状；Python validator 额外执行跨字段 reciprocal、lineage、units、normal、triangle winding/面积、label pairing、leakage 和 qualification 检查。payload 路径必须是包目录内规范 POSIX 相对路径；export 会核验母体、fragment、observed mesh 与 mask 的真实文件 SHA，observation reader 只读取并核验观测 mesh/mask，offline GT reader 另外核验母体/fragment meshes。JSON 与 payload 读取上限分别为 32 MB、单文件 8 MB/总计 30 MB；实体数组有 100,000 行上限、每个 surface 有 200,000 个 triangle 上限。这里的 fixture 与 focused test records 都小于 64，不会启动 solver 或批量造数据。

## 固定 fixtures 与验证

`fixtures/sample_fracture_artifact.json` 是可实际读取的几何完整小样例，并明确为 `FIXTURE_GEOMETRY_ONLY`；其中母体/fragment/observed mesh 与 visibility mask 的 digest 对应 `fixtures/payload/` 下的迷你文件。load/material/discretization 中的 synthetic provenance hashes 仅是测试标记，不是实际 solver、材料实验或产品几何证据。其余固定 case 分别演示 missing labels、跨 fragment 错配、单位 mismatch、非标准 JSON NaN、GT 字段泄漏，以及显式 UNKNOWN 与 genuine zero 的区分。

在 repo 根目录运行：

```powershell
& 'relicweave-synlab/.venv/Scripts/python.exe' -m unittest -v pgfrag.experiments.fracture_artifact_contract_20261006.test_fracture_artifact_contract
```

已有 `synlab.qa.source_readers.pgfrag_result.read_pgfrag_result()` / chain dispatch 仍只验证 `pgfrag.result.v1` endpoint rows，并保持 `PG_DIAGNOSTIC_ONLY`。它不会因为这个新增契约就推断出物理 fragment、裂面、event 或 R2。本目录是后续 A12 adapter 的输入约定，不宣称已集成或已获物理验证。
