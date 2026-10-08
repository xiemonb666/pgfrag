# Focused MLS3D acceptance pytest

One invocation only; run with the existing Syn venv. No solver, benchmark, or geometry workflow was invoked.

Artifact status remains `MATHEMATICAL_PRIMITIVE_ONLY` / `NOT_COUPLED`; no fracture or R2 claim is made.

## Command

Actual argv is recorded as an array in `RUN_METADATA.json`.

```text
C:\Users\xiemo\Desktop\ShardForge\relicweave-synlab\.venv\Scripts\python.exe -m pytest -q pgfrag/experiments/mls3d_kernel_20261005/test_mls3d_mls.py --junitxml=pgfrag/tmp/review_20261005/mls3d_kernel_20261005/focused_test/junit.xml
```

Working directory: `C:\Users\xiemo\Desktop\ShardForge`
Timeout enforced on pytest process: 30,000 ms
Return code: 0
Elapsed wall time: 1,452 ms
Actual testcase count from JUnit `//testcase`: 6

## Captured output

stdout (`stdout.txt`):
```text
......                                                                   [100%]
6 passed in 0.17s
```

stderr (`stderr.txt`):
```text
(empty)
```

JUnit XML: `junit.xml`

## Source freeze

Before-run snapshots are in `source_snapshot/`; SHA-256 manifests are `PRE_RUN_SHA256.txt` and `POST_RUN_SHA256.txt`. `SOURCE_IDENTITY_CHECK.txt` compares pre-run source, post-run source, and snapshot hashes.
