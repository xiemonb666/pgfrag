# Extrinsic capture preflight

`capture_preflight.py` performs a bounded native API capture against the frozen `coarse654_20261004_01` 96-T4 coupon. It has two explicit modes:

- `extrinsic` repeats the old probe's model initialization, mass assembly, interface registration, and automatic-insertion eligibility update. It must report eight eligible facets, zero already inserted facets, and zero cohesive elements.
- `intrinsic-capability` starts a separate fresh model with only the saved `interface` surface registered before initialization. Its eight cohesive elements exist only to establish native field getter shapes; this run is not an extrinsic insertion or event replay.

Neither mode applies boundary conditions, calls `checkCohesiveStress`, calls `solveStep`, advances time, or assembles internal forces. The setup calls `assembleMassLumped()` once, matching the frozen probe setup. Native arrays are copied into the paired data directory as compressed NPZ files. Getter errors and empty arrays are retained as such.

Run each mode from PowerShell at the repository root. The wrapper fixes the thread count, limits wall time to 60 seconds, limits each created file to 16 MiB using the Linux 512-byte file-size blocks, and captures the bounded console log:

```powershell
wsl.exe -d Ubuntu -- bash -lc 'cd /mnt/c/Users/xiemo/Desktop/ShardForge && ulimit -f 32768 && export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 && timeout --signal=TERM --kill-after=2s 60s third_party/akantu-runtime/venv/bin/python pgfrag/experiments/continuation_20261016/extrinsic_capture_preflight/capture_preflight.py --mode extrinsic > pgfrag/data/continuation_20261016/extrinsic_capture_preflight/extrinsic_stdout.txt 2>&1; rc=$?; wc -c pgfrag/data/continuation_20261016/extrinsic_capture_preflight/extrinsic_stdout.txt; cat pgfrag/data/continuation_20261016/extrinsic_capture_preflight/extrinsic_stdout.txt; exit $rc'
```

```powershell
wsl.exe -d Ubuntu -- bash -lc 'cd /mnt/c/Users/xiemo/Desktop/ShardForge && ulimit -f 32768 && export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 && timeout --signal=TERM --kill-after=2s 60s third_party/akantu-runtime/venv/bin/python pgfrag/experiments/continuation_20261016/extrinsic_capture_preflight/capture_preflight.py --mode intrinsic-capability > pgfrag/data/continuation_20261016/extrinsic_capture_preflight/intrinsic_stdout.txt 2>&1; rc=$?; wc -c pgfrag/data/continuation_20261016/extrinsic_capture_preflight/intrinsic_stdout.txt; cat pgfrag/data/continuation_20261016/extrinsic_capture_preflight/intrinsic_stdout.txt; exit $rc'
```

The frozen paths, input/source/wheel hashes, public binding documentation, getter results, array shapes and hashes, no-step contract, and next-step capture ordering are in [the paired report](../../../data/continuation_20261016/extrinsic_capture_preflight/REPORT.md) and the two machine records.
