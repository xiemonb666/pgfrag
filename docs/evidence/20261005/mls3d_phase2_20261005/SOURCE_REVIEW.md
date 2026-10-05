# Phase 2 source review — direct MLS shape functions

**Artifact status:** `MATHEMATICAL_PRIMITIVE_ONLY`  
**Coupling status:** `NOT_COUPLED`  
**Fracture status:** `NO_ACTUAL_FRACTURE`  
**R2 status:** `NO_R2_CLAIM`  
**Execution status:** frozen phase2 API acceptance source passed one focused
pytest run (7 JUnit testcases, exit 0). No solver, operator, benchmark, or
geometry workflow was run.

## Phase 1 preservation

Before adding this API, phase1's already accepted `mls3d_mls.py`,
`test_mls3d_mls.py`, and `README.md` were copied into `phase1_snapshot/`. Their
hashes match phase1's focused-test snapshots and run record: exactly one
Syn-venv pytest run, exit code 0, six testcase elements in JUnit. No phase1
source file was edited for phase2. The snapshot digests are in
`PHASE1_SNAPSHOT_SHA256.txt`; source identity for phase2 can be checked against
`PHASE2_SOURCE_SHA256.txt`.

## Scope and source reuse

The new API is in
`pgfrag/experiments/mls3d_kernel_20261005/mls3d_shape.py`. It imports the
reviewed phase1 `MLS3DConfig`, `MLSStatus`, `wendland_c2`, hard caps, and strict
position/mask validators. It does not change the phase1 scalar/vector evaluator
or its passed-test sources.

The target use is an EFGM-style assembly caller that needs a nodal shape row
without interpolating an `N`-column identity matrix. The assembly agent
confirmed that dense source rows in canonical global-ID order are compatible
with its adapter; successful masked/out-of-support columns are zero and the
gradient orientation is `(query, source, axis)`.

For each query, source basis values `p_i` are fixed, while the query basis and
radial weights vary. The moment matrix is `A=sum_i w_i p_i p_i.T`. Solving the
dual system `A.T c=p(x)` yields `Phi_i=w_i p_i.T c`. Differentiating that solve
gives `A.T dc/dx_j=dp/dx_j-(dA.T/dx_j)c`; the code then differentiates both
the weight and dual coefficients. No `N-by-N` identity field is materialized.

Per-query outputs are dense and bounded at `Q<=500`, `N<=500`. The maximum
weight array has 250,000 entries and the maximum gradient array has 750,000
entries. Input shape/cap/dtype/finite checks run before conversion and result
allocation. Position inputs are copied after validation; visibility is copied
and made read-only. Result arrays use immutable bytes-backed storage, so
writeability cannot be re-enabled. The same no-pseudoinverse SVD rank and condition checks,
duplicate rejection, poor-support result, fixed-mask stability rule, and
support-boundary UNKNOWN guard apply.

## Review gates

The acceptance source is
`pgfrag/experiments/mls3d_kernel_20261005/test_mls3d_shape.py`. It covers:

- partition of unity and zero gradient sum;
- exact affine coordinate reconstruction and identity coordinate gradient;
- scalar/vector nonlinear value and gradient equivalence with phase1's
  established interpolation interface;
- proper-rotation covariance and independent finite differences of each
  shape-gradient axis;
- UNKNOWN outcomes for poor support, rank loss, duplicate support,
  ill-conditioning, radial and visibility switches;
- fail-closed source/query budgets, including a huge zero-stride view rejected
  from shape metadata;
- immutable result buffers, source-input snapshot isolation, and a successful
  locally stable mask row with hidden-source weights and gradients set to zero.

The focused run receipt is
`pgfrag/tmp/review_20261005/mls3d_kernel_20261005/shape_api_phase2/phase2_focused_test/RUN_RECORD.md`.
The run passed all seven testcase elements. Passing these checks verifies only
this shape-function primitive. It does not validate EFGM stiffness assembly,
mechanics coupling, fracture evolution, surfaces, fragment labels, R2, or
production readiness.
