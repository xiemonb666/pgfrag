#!/usr/bin/env bash
set -euo pipefail

ROOT=/mnt/c/Users/xiemo/Desktop/ShardForge
PY="$ROOT/third_party/akantu-runtime/venv/bin/python"
EXP="$ROOT/pgfrag/experiments/free_cohesive_patch_20261015"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "$EXP"
exec "$PY" "$EXP/run_native_free_patch.py"
