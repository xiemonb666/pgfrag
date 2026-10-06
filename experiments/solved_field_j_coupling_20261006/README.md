# Solved-field to Mode III J readback

`read_solved_field_j.py` reads the frozen two-level static-solve JSON, binds
each coefficient to its semantic canonical ID, and samples the existing
54-DOF continuous MLS/front-PU provider. It does not assemble a stiffness
matrix or invoke the static solver.

The planned provider batch contains 224 midsection points from two ellipse
scales, two rotations, and GL(12/16) per smooth arc. It adds 24 reflected
controls at `z = +/- 0.0025 m`, for 248 provider queries total. The worker has
a 60-second process timeout and the result JSON has an 8 MiB cap.

The full solved gradient is saved and passed unchanged to the pinned
`j_integral_precopy_guard_v3.py`. An integrator `UNKNOWN` status and reason are
preserved. Any numeric result is a midsection diagnostic only; it does not
establish full-domain 2D behavior, physical `G`, or crack-growth readiness.
The closed-form `KIII^2/(2 mu)` value appears only as a comparison reference.
The report keeps GL(12) versus GL(16) sensitivity within each fixed ellipse
separate from size/rotation comparisons made at a fixed GL order. It does not
combine those effects into one contour spread, and neither comparison is
treated as a path-independence qualification.

The provider batch reports 248 sampled rows. The adapter contracts its returned
value and gradient matrices directly; `contract_evaluation_count` records the
provider's actual `evaluate` count, which is zero for this matrix path.

The source and dependency hashes are listed in
`pgfrag/tmp/review_20261006/solved_field_j_coupling_20261006/source_freeze.json`.
The source-freeze record is marked as awaiting main-agent review. Do not run the
adapter before approval. After approval, the single intended command is:

```powershell
third_party/vesselforge-runtime/venv/Scripts/python.exe -B pgfrag/experiments/solved_field_j_coupling_20261006/read_solved_field_j.py --execute-readback
```

The command refuses to overwrite an existing run output. Its default output
directory is `pgfrag/tmp/review_20261006/solved_field_j_coupling_20261006/runs/run01/`.
