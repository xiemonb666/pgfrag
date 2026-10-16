#!/usr/bin/env bash
set -euo pipefail

ROOT=/mnt/c/Users/xiemo/Desktop/ShardForge
EXPERIMENT="$ROOT/pgfrag/experiments/continuation_20261016/reference_v6_preflight"
PYTHON="$ROOT/third_party/akantu-runtime/venv/bin/python"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_DYNAMIC=FALSE

cd "$EXPERIMENT"
exec timeout --foreground --signal=TERM --kill-after=5s 300s \
  taskset -c 0 "$PYTHON" \
  "$EXPERIMENT/inputs/cohesive_vessel.py" \
  --case "$EXPERIMENT/inputs/case.json" \
  --out "$EXPERIMENT/native_preflight/run01" \
  --preflight-only
