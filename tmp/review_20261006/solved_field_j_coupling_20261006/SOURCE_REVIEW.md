# Source review: solved-field J coupling

Status: `SOURCE_IMPLEMENTATION_READY_AWAITING_ROOT_REVIEW`

This source-only change adds the readback adapter and experiment notes. It does
not edit the frozen static solve, J-integrator, provider, fixture, production
code, or documentation artifacts.

The adapter checks the exact existing solve-result hash, solver-source hash,
source-freeze hash, provider dependencies, and J-integrator v3 hash before
loading numerical dependencies. It rebuilds the frozen semantic registry,
creates one provider batch, and maps the stored coefficient vector through
canonical IDs before contracting `phi` and all three components of `dphi`.
The adapter never assembles `K` and never launches the static solver.

The 224 J-contour samples are the established two-scale, two-rotation,
GL(12/16)-per-arc set. Twenty-four additional points use twelve mirrored
`z=+/-0.0025 m` pairs. All 248 points are submitted in one bounded provider
batch. A subprocess timeout enforces 60 seconds. Result serialization enforces
8 MiB. The result stores contour coordinates, normals, `ds` weights, solved
`u_z`, the unmodified `(grad_x, grad_y, grad_z)`, the integrator result and
reason, comparisons with the analytic scalar, GL(12/16) sensitivity within a
fixed ellipse, and separate cross-size/rotation contour sensitivity at each
fixed GL order. No combined spread mixes geometry and quadrature-order effects.

The adapter passes exact contiguous NumPy arrays to J-integrator v3. If v3
returns `UNKNOWN`, including `field_must_be_z_invariant`, the adapter preserves
that status and reason. There is no tolerance override or gradient projection.
Every numerical J is labeled a midsection diagnostic and cannot qualify full
3D behavior, physical `G`, or crack growth. Parent promotion requires a zero
child exit plus the complete schema, frozen input identity, 248-query batch,
eight contour slices, two solved levels, semantic coefficient bindings,
finite numeric payload, and separate sensitivity records. Failed or
incomplete results remain `UNKNOWN` with the worker staging artifact kept.
The provider reports 248 rows sampled in one factory batch. The adapter uses
the returned `values` and `gradients` matrices directly, so
`contract_evaluation_count` records the actual per-point `evaluate` count of
zero and never calls `evaluate` to manufacture a query count.

Budget declaration: 248/256 provider queries; 60.0-second hard wall cap with a
59.0-second worker timeout; 8 MiB maximum result JSON. The frozen provider
`values` and `gradients` arrays occupy about 0.41 MiB together; query
coordinates add about 0.006 MiB. The adapter also records a tracked-array
estimate. Provider MLS/exact-field internals and Python object overhead are
excluded; serialized output size is capped at runtime.

Execution and numerical validation have not been run. Wait for main-agent
approval before one readback execution.
