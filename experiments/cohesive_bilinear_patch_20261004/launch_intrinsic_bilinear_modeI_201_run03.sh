#!/usr/bin/env bash
set +e
PY=/mnt/c/Users/xiemo/Desktop/ShardForge/third_party/akantu-runtime/venv/bin/python
DRIVER=/mnt/c/Users/xiemo/Desktop/ShardForge/pgfrag/experiments/cohesive_bilinear_patch_20261004/run_intrinsic_bilinear_modeI_201_run03.py
OUT=/home/xiemo/cohesive_bilinear_patch_20261004/controlled_intrinsic_bilinear_modeI_201_run03
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$PY" "$DRIVER"
rc=$?
printf '%s\n' "$rc" > "$OUT/python_process_exit_code.txt"
exit "$rc"