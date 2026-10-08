# 3-D linear MLS primitive — source review

**Artifact status:** `MATHEMATICAL_PRIMITIVE_ONLY`  
**Fracture status:** `NO_ACTUAL_FRACTURE`  
**Production status:** `NOT_PRODUCTION`  
**Execution status:** source and tests frozen for review; no test, CPM solver,
or fracture run was executed.

## Scope

The implementation at
`pgfrag/experiments/mls3d_kernel_20261005/mls3d_mls.py` evaluates bounded,
three-dimensional, linear moving-least-squares (MLS) interpolation and analytic
spatial gradients. It supplies a small numerical primitive that an explicit
surface mechanics experiment could call. It does not implement the rest of the
Ai–Augarde adaptive cracking-particle method, particle visibility, meshfree
support modification, crack-front state, fracture criterion, splitting,
quadrature, or a solver.

The root architecture decision selects a mechanics-coupled explicit crack
surface as the future geometric-truth producer. The chosen branch is a new
formulation, rather than an interpretation of the existing ordinary
state-based peridynamic state. Nothing in this primitive records crack
geometry or mechanics-grounded surface labels.

## Primary-source alignment and limit

The source extraction for Ai and Augarde's official doctoral thesis reports
that §2.2 (PDF pages 57–58, print pages 32–33) gives the 3-D linear basis
`p=[1,x,y,z]^T`, the MLS form `Phi=p^T A^-1 B`, and differentiation of the
matrix inverse. This implementation uses that 3-D polynomial space and
analytic differentiation of the equivalent linear solve `A a=b`:

`A (da/dx_j) = (db/dx_j) - (dA/dx_j) a`.

The thesis gives a quartic compact spline as an example/general recommendation.
The exact kernel and support radius used in Ai and Augarde's 2018 §4.2 3-D
benchmark were not established by the inspected material. As recorded in the
root architecture decision, do not call this an exact benchmark reproduction.
The root-selected Wendland C2 kernel and support radius are engineering inputs
for this standalone primitive; they are not attributed to the paper. No
support radius from a related 2-D paper is transferred to the 3-D example.

Primary-source record: W. Ai, W. G. (Bill) Price, and I. S. Howard, *An
adaptive cracking-particle method providing explicit and accurate description
of 3D crack surfaces*, Durham repository primary record:
<https://durham-repository.worktribe.com/output/1335971/an-adaptive-cracking-particle-method-providing-explicit-and-accurate-description-of-3d-crack-surfaces>.
The official thesis identified in the parallel source extraction is
<https://etheses.dur.ac.uk/12679/1/WAi_PhDthesis.pdf>.

## Formulation implemented

For source locations `X_i`, query `x`, fixed reference origin `o`, and fixed
positive scalar scale `s`, the basis is

`p_i=[1,(X_i-o)/s]`, `p(x)=[1,(x-o)/s]`.

For the visible source set within support radius `R`, define

`q_i=||X_i-x||/R`,  
`w_i=(1-q_i)^4(4q_i+1)` for `0<=q_i<1`, otherwise zero,  
`A(x)=sum_i w_i p_i p_i^T`,  
`b(x)=sum_i w_i p_i u_i`,  
`a(x)=solve(A(x),b(x))`,  
`u_hat(x)=p(x)^T a(x)`.

The selected radial derivative is `dw/dq=-20q(1-q)^3` on support. For each
spatial axis, the code computes

`dA/dx_j=sum_i (dw_i/dx_j) p_i p_i^T`,  
`db/dx_j=sum_i (dw_i/dx_j) p_i u_i`,  
`da/dx_j=solve(A, db/dx_j-(dA/dx_j)a)`,  
`du_hat/dx_j=(dp(x)/dx_j)^T a+p(x)^T(da/dx_j)`.

Thus the gradient includes weight changes, the dependency of the moment matrix
and right-hand side, the derivative of the solved coefficients (equivalent to
the inverse derivative), and the query-basis derivative. The implementation
uses `numpy.linalg.solve` after SVD rank and condition checks and never
substitutes a pseudo-inverse.

The source uses fixed reference coordinates rather than a query-centered
changing polynomial frame. This preserves a consistent differentiable basis
across queries. A single isotropic scale is required for rigid-rotation
covariance. Ill-scaled user coordinate choices may be rejected as
ill-conditioned; the caller chooses and holds the origin, scale, and support
radius fixed for its reference domain.

## Boundaries and unknowns

Only fixed-mask local gradients are supported. A supplied boolean visibility
mask is snapshotted. Queries without an explicit local-stability assertion
return UNKNOWN; callers must likewise label known mask-boundary or
mask-switching queries UNKNOWN. The primitive has no visibility boundary
detector because it receives only a fixed mask vector.

Radial support-boundary queries within a configured tolerance (including exact
`q=1`) return UNKNOWN. Below-four support, local duplicate source positions,
rank loss, excessive condition number, and failed solves also return UNKNOWN
with NaNs. Inputs must be plain NumPy ndarrays. Shape and finite checks precede
dtype conversion; source/query rows are capped at 500 and value components at
three. Cap violations raise `ValueError` with `budgetRejected` before casts,
finite scans, or result allocation. Array-like conversion hooks are not called.
The acceptance source includes 501-row rejects, a huge zero-stride value view,
and a hostile `__array__`/`tolist` counterexample. These checks avoid
manufacturing gradients where support is insufficient or resource limits are
exceeded.

## Review status

The source-only acceptance file is
`pgfrag/experiments/mls3d_kernel_20261005/test_mls3d_mls.py`. It covers:

- exact constant and affine field reproduction and gradient;
- rigid translation and finite proper rotation with zero Green–Lagrange strain;
- covariance under a generic proper rotation of layout, query, values, and
  reference origin;
- an independent central finite-difference check away from switching surfaces;
- explicit UNKNOWN outcomes for poor support, rank deficiency, duplicates,
  support boundary, and visibility-mask switch;
- hard `budgetRejected` errors for 501-row sources/queries, more than three
  components, a huge zero-stride broadcast field view, and a non-ndarray object
  with guarded `__array__` and `tolist` hooks.

These tests are authored but not executed per the review freeze. Passing them
would verify only this numerical primitive; it would not establish an explicit
crack surface, the 2018 benchmark, fracture propagation, surface export,
material calibration, or production readiness.
