# Phase 2 direct MLS shape-function contract

**Artifact status:** `MATHEMATICAL_PRIMITIVE_ONLY`  
**Coupling status:** `NOT_COUPLED`  
**Fracture status:** `NO_ACTUAL_FRACTURE`  
**R2 status:** `NO_R2_CLAIM`

## API

Implementation: `pgfrag/experiments/mls3d_kernel_20261005/mls3d_shape.py`.

```python
shape_functions_3d(
    source_positions: np.ndarray,  # (N, 3)
    query_positions: np.ndarray,   # (Q, 3)
    config: MLS3DConfig,
    *,
    visibility_mask: np.ndarray | None = None,  # (N,), bool
    mask_locally_constant: bool | np.ndarray | None = None,  # scalar or (Q,), bool
) -> ShapeFunctions3DResult
```

The result has `weights[Q,N]` for `Phi[query,source]`, `gradients[Q,N,3]` for
`dPhi[query,source]/dX[axis]`, `statuses[Q]`, a computed `reasons[Q]`,
`support_counts[Q]`, and `condition_numbers[Q]`. `reasons[q]` is `None` when
`statuses[q] == MLSStatus.OK`; otherwise it is the exact status string. `ok[Q]`
provides boolean success flags. All returned arrays are backed by immutable
bytes; callers cannot re-enable writes with `setflags(write=True)`.

Rows use source input order. For the EFGM adapter, pass source positions in
canonical global-ID order `0..N-1`; then every result column already maps to
that global source ID. Successful rows have zero entries for hidden and
out-of-support sources. Any UNKNOWN row has NaN weights and gradients for all
N sources so downstream assembly cannot treat missing support as zero.

## Bounded input and snapshot behavior

Inputs must be plain NumPy ndarrays, not array-like objects or ndarray
subclasses. `N<=500`, `Q<=500`, and both position arrays have shape `(count,3)`.
Shape, dtype, and finite-value validation happens before dtype conversion and
before result arrays are allocated. A cap violation raises `ValueError` with a
`budgetRejected` marker; it does not create an oversized UNKNOWN result.

After validation, the API takes private float64 copies of source and query
positions and a private read-only copy of the visibility mask. The caller may
then reuse its input buffers without changing the current calculation. The
fixed mask is not a crack state or visibility law. With a mask, the caller must
explicitly assert a locally constant region for each query; unspecified or
false entries return `UNKNOWN_VISIBILITY_SWITCH`. Queries near an active radial
support boundary, including exact `q=1`, return `UNKNOWN_SUPPORT_BOUNDARY`.

## Shape solve and derivative

For each query, with the same fixed reference-frame basis, compact Wendland C2
weight, radius, and guards as phase 1, assemble

`A=sum_i w_i p_i p_i.T`, `q(x)=p(x)`, `A.T c=q`, `Phi_i=w_i p_i.T c`.

For spatial axis `j`,

`A.T dc/dx_j = dp/dx_j - (dA.T/dx_j)c`,
`dPhi_i/dx_j = (dw_i/dx_j)(p_i.T c) + w_i(p_i.T dc/dx_j)`.

The implementation uses four-by-four solves per query and expands one dense
`Q-by-N` nodal shape matrix and `Q-by-N-by-3` gradient array. It does not form
an `N-by-N` identity field or call the scalar interpolator N times. Gradients
are ordered `(query, source, reference axis)` and are directly consumable by
an EFGM-style assembly adapter.

## UNKNOWN behavior and claims

Less than four active sources, duplicate local source positions, rank loss,
excessive condition number, a failed linear solve, non-finite shape output,
visibility switching, and support-boundary switching return an UNKNOWN status
with NaN nodal rows. No pseudo-inverse or regularization masks missing support.

This function only returns interpolation weights and spatial derivatives. It
does not assemble an elastic operator, evolve a displacement, create a crack,
or produce a fracture surface, event, fragment label, or R2 pair. Its status is
`NOT_COUPLED`; no mechanics or fracture truth is claimed.

## Acceptance source status

`pgfrag/experiments/mls3d_kernel_20261005/test_mls3d_shape.py` contains source
checks for partition unity, zero sum-gradient, affine coordinate and identity
gradient reproduction, equivalence to the existing scalar/vector MLS path,
rotation covariance, independent central finite differences of shape
gradients, UNKNOWN support/conditioning/switch paths, and pre-allocation
budget rejection. It also checks output write protection, source-snapshot
isolation, and fixed-mask success with hidden-node weights and gradients equal
to zero. All seven testcase elements passed in the one authorized focused
Syn-venv run; raw streams, JUnit, argv, timing, and source identity receipts are
in `phase2_focused_test/`.
