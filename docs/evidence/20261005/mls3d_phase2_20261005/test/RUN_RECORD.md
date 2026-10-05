# Phase2 focused MLS shape pytest

One invocation only, using the Syn venv. No solver, operator, crack benchmark, or geometry workflow was run.

Artifact status: `MATHEMATICAL_PRIMITIVE_ONLY` / `NOT_COUPLED`; no fracture or R2 claim.

## Command

Actual argv is recorded as an array in `RUN_METADATA.json`.

```text
C:\Users\xiemo\Desktop\ShardForge\relicweave-synlab\.venv\Scripts\python.exe -m pytest -q pgfrag/experiments/mls3d_kernel_20261005/test_mls3d_shape.py --junitxml=pgfrag/tmp/review_20261005/mls3d_kernel_20261005/shape_api_phase2/phase2_focused_test/junit.xml
```

Working directory: `C:\Users\xiemo\Desktop\ShardForge`
Timeout enforced on pytest process: 30,000 ms
Return code: 0
Elapsed wall time: 1,465 ms
Actual testcase count from JUnit //testcase: 7

## Captured output

stdout (`stdout.txt`):
```text
.......                                                                  [100%]
7 passed in 0.21s
```

stderr (`stderr.txt`):
```text
(empty)
```

JUnit XML: `junit.xml`

## Source identity

`PRE_RUN_PHASE2_SOURCE_SHA256.txt` and `POST_RUN_PHASE2_SOURCE_SHA256.txt` cover the frozen API/test and phase1 core dependency. The phase1 module/test/README pre/post hashes are in the corresponding `*_PHASE1_ORIGINAL_SHA256.txt` files. `POST_RUN_SNAPSHOT_CHECK.txt` verifies both phase2 snapshot identity and all three phase1 original snapshots.
