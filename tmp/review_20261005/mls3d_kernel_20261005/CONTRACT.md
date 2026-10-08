# MLS3D primitive contract

**Artifact status:** `MATHEMATICAL_PRIMITIVE_ONLY`  
**Fracture status:** `NO_ACTUAL_FRACTURE`  
**Production status:** `NOT_PRODUCTION`

## API boundary

The source is `pgfrag/experiments/mls3d_kernel_20261005/mls3d_mls.py`.
`evaluate_mls3d(source_positions, source_values, query_positions, config, ...)`
returns one interpolation result and one spatial gradient per query. Positions
are three-dimensional coordinates in one material/reference frame. Source row
`i` is paired with value row `i`; a one-dimensional value array is a scalar
field and a `(n_source, n_component)` array is a vector or multi-field.

The polynomial basis is `p(X)=[1,(X-o)/s]`, where the fixed configuration
supplies one reference origin `o`, one positive scalar length scale `s`, and
one positive compact support radius `R`. Reuse the same origin and scale for a
single reference domain. `s` is isotropic so a proper rigid rotation does not
change the interpolation metric. The compact Wendland C2 weight uses
`q=|X_i-x|/R` and contributes only for `q<1`.

## Output and failure behavior

Successful rows have `status=OK`, finite `values`, finite `gradients`, a support
count, and a moment-matrix condition number. Gradient entries have layout
`gradients[query, field_component, reference_axis]` and mean
`d(value_component)/dX_axis`. Returned arrays are read-only.

Every failed query has a non-OK status and NaN values and gradients. The
primitive does not return a least-squares pseudo-inverse or regularized answer
for insufficient support. Fewer than four visible in-radius sources, duplicate
local source positions, rank loss, excessive condition number, or a failed
linear solve remain explicit UNKNOWN outcomes. Exceeding 500 sources, 500
queries, or three field components raises `ValueError` with a `budgetRejected`
marker before numeric conversion or result allocation; no oversized result is
constructed.

Inputs must be plain `numpy.ndarray` objects. Array-like objects and ndarray
subclasses are rejected before invoking conversion hooks. Source and query
positions are capped at 500 rows each; field values are capped at 500 rows and
1–3 components. Exceeding a cap raises `ValueError` with a `budgetRejected`
marker before finite scanning, dtype conversion, or output allocation. Shape
metadata is inspected before conversion, so even a zero-stride view with a
huge logical component axis is rejected cheaply. Invalid shapes, non-finite
input fields, malformed masks, and invalid configuration values also raise
`ValueError`. A valid zero-query batch returns an empty result.

## Visibility-mask semantics

An optional plain NumPy boolean `visibility_mask` selects sources before support assembly.
The function snapshots the mask and treats it as fixed; it does not calculate
visibility and does not interpret a mask as a crack. If a mask is supplied,
the caller must also set `mask_locally_constant=True` or pass a plain NumPy
boolean vector with one entry per query. An unspecified or false entry returns
`UNKNOWN_VISIBILITY_SWITCH`.
Callers must mark queries on any externally known mask boundary or switching
surface false. This primitive cannot infer an external mask's boundary from
the mask vector alone.

Queries within the configured tolerance of the radial support boundary return
`UNKNOWN_SUPPORT_BOUNDARY`, including exact `q=1`. This guard is conservative:
the selected C2 kernel itself has a continuous first derivative there, while
neighborhood/mask transition validation is outside this primitive.

## Analytic differentiation

For a fixed visible support, the code differentiates the weighted moment
matrix and right-hand side, solves the differentiated system, and includes the
query-basis derivative. Gradients are valid only while the source set,
reference frame, configuration, and visibility mask stay fixed locally.

`green_lagrange_strain` is a separate kinematic helper with
`E=0.5(F^T F-I)`, `F=I+grad(u)`. It provides a finite-strain reproduction
check for rigid translations and rotations. It does not create mechanics or
advance a crack.

## Required acceptance checks

`test_mls3d_mls.py` lists checks for constant and affine reproduction,
translation and finite-rotation zero strain, generic proper-rotation
covariance, an independent central finite-difference gradient comparison away
from support switches, and explicit UNKNOWN behavior for sparse, rank-deficient,
duplicate, radial-boundary, and visibility-switch cases. It also records hard
budget rejection for 501 sources/queries and more than three field components,
a huge zero-stride broadcast field view rejected before copying, and a hostile
array-like whose `__array__`/`tolist` hooks must remain untouched. Tests are
present but have **not** been executed pending root review.
