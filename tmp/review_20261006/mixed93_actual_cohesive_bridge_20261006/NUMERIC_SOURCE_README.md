# Mixed93 actual cohesive bridge source package

This package adds a root-go-gated worker for one fixed Mode-I crack-sheet coupling case. The frozen `case.json`, 281-line static checker, and original `SOURCE_PINS.json` remain unchanged. The case is an illustrative elastic cohesive bridge; it is not a material calibration, process-zone qualification, crack-growth result, or physical fracture claim.

`run_actual_cohesive_bridge_once.py` defaults to SOURCE_READY. That path runs the existing source-only checker, parses the new Python sources, verifies the source hashes, and hashes the accepted NPZ bytes. It does not import NumPy or parse/open the NPZ archive, call an MLS provider, or solve. The numerical worker has no top-level numerical imports and refuses execution unless `NUMERIC_SOURCE_PINS.json` contains an explicit root-go record and the matching command token.

The proposed authorized command is:

```powershell
& 'C:\Users\xiemo\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B pgfrag/experiments/mixed93_actual_cohesive_bridge_20261006/run_actual_cohesive_bridge_once.py
```

No numerical execution is enabled by this source package. After root review, root may set `execution_authorized=true` and fill `root_go.authorized`, `root_go.token`, author, and UTC in the numeric pins. Only then may the supervisor be invoked with `--execute --root-go <exact-token>`. It creates one exclusive `numeric_runs/run01`, keeps a one-attempt lock, applies one 60-second deadline from preflight through the child process, captures stdout/stderr as raw bytes, and writes an outer receipt with exact file sizes and hashes. The 12 MiB combined-output cap includes result JSON, raw stdout/stderr, lock, and receipt, with receipt byte accounting solved to a serialized-size fixed point. It never retries or rewrites the worker result.

The worker is designed to use the cached `pure_I_q_solution`, `pure_I_boundary_force_N`, `pure_I_physical_gauge_target_m`, and `pure_I_lambda_physical`, each scaled by 0.01; it keeps `pure_I_q_exact` diagnostic-only. It forms 16 new direct continuous BASE MLS queries and 16 finite-sheet Gauss points, builds shared-point upper/lower traces in the frozen 93 semantic DOF registry, calls the pinned cohesive primitive at most three times, and forms one 99-order coupled KKT solve with the analytic positive `KG[82,82]=2 K0 R L=20000 N/m`. Result JSON retains query shapes, per-point `H+`, `H-`, `J`, six rigid modes, cached and coupled matrices/vectors, and baseline/final cohesive point ledgers so reviewers can reconstruct the coupling. Failed or unsupported gates produce `NOT_ACCEPTED` with details. No raw-stiffness repair, ridge, pseudoinverse, coefficient pinning, load projection, or discarded shear response is permitted.

The internal generalized cohesive force is the derivative of cohesive potential and enters the left-hand side with positive `KG`. The primitive's `+A t n` positive-face force is internal bookkeeping; the physical traction applied to the upper face is `-t n`. The sign audit checks this directly and also verifies the energy derivative, so action/reaction by itself cannot pass the direction gate.
