#!/usr/bin/env bash
set +e
PY=/mnt/c/Users/xiemo/Desktop/ShardForge/third_party/akantu-runtime/venv/bin/python
DRIVER=/mnt/c/Users/xiemo/Desktop/ShardForge/pgfrag/experiments/cohesive_bilinear_patch_20261004/init_preflight_dump_v3.py
OUT=/mnt/c/Users/xiemo/Desktop/ShardForge/relicweave-synlab/tmp/review_20261004/implementation/cohesive_bilinear_patch_20261004/preflight_dump_v3
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$PY" "$DRIVER"
rc=$?
printf '%s\n' "$rc" > "$OUT/python_process_exit_code_v3.txt"
exit "$rc"