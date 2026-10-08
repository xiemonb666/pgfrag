# Source freeze and deferred verification plan

The operator and its focused acceptance tests are
`pgfrag/experiments/efgm_elastic_operator_20261005/efgm_elastic_operator.py`
and `pgfrag/experiments/efgm_elastic_operator_20261005/test_efgm_elastic_operator.py`.
Both SHA-256 values are recorded in `source_freeze.json`. One pinned pytest
invocation ran all 11 tests in 0.52 seconds under a 30-second timeout; the
JUnit report, stdout/stderr, return code, and pre/post hashes are under
`focused_test/`. No solver or benchmark was run.

## Bounded fixture

Use a cube `[-1,1]^3`, with 27 coefficient sites at the tensor product of
`{-1,0,1}`. Divide each axis into two unit cells and use two-point Gauss
quadrature on each axis in each cell. This gives eight points per cell, 64
points total, with volume weight `1/8` at each point and total volume 8. Use
`E=1`, `nu=0.3`, and dimensionless tension scale `sigma=1e-3` if a later load
path is added. Maximum fixture size stays far below 375 DOFs and 256 volume
points. The focused source tests completed in 0.52 seconds; this budget does
not include a solver or mechanics benchmark.

## Checks to run after review

1. **Input and allocation guards:** reject ndarray subclasses and conversion
   hooks; reject 126 coefficient sites and 257 quadrature points before any
   conversion copy or operator allocation. Accept the 27-site, 64-point
   fixture. Check malformed dense all-node shape ID ordering/length and ensure
   a non-OK support status yields no matrix while preserving the provider's
   original UNKNOWN reason code in the diagnostic.
2. **Material positivity:** for representative stable values including
   `(E,nu)=(1,0.3)`, verify the 6x6 engineering-shear `D` is symmetric positive
   definite. Reject `E<=0` and `nu` outside `(-1,0.5)`.
3. **Galerkin symmetry and six modes:** assemble the independent global
   affine provider on the fixture; require
   `||K-K.T||_F/max(||K||_F,1) <= 1e-12` and
   `||K R||_F/max(||K||_F*||R||_F,1) <= 1e-10` for the matrix `R` of three
   translation and three infinitesimal rotation vectors. Check all six are
   independent on the nondegenerate grid. This provider spans only four global
   polynomial modes and can induce other null modes; this check verifies that
   the required rigid modes are present, not that these are the only modes or
   that a local MLS/EFGM support is stable.
4. **Affine constant-strain energy:** set coefficient values to `A X + b` for
   a fixed nonsymmetric `A`. Compare `0.5*d.T@K@d` with the independent tensor
   result `0.5*V*epsilon:C:epsilon`, where `epsilon=sym(A)` and `V=8`. Include
   normal and shear components so the engineering-shear convention is
   exercised; require relative energy difference at most `1e-10`.
5. **Load linearity belongs to a later solver slice:** this operator accepts
   no load and computes no reactions. When a static system is added, hold K
   and physical constraint rows fixed, scale prescribed traction by two, and
   check displacement and reactions scale by two while strain energy scales
   by four. Do not claim this check from operator-only assembly.
6. **Finite-input overflow counterexamples:** use a bounded one-point case
   with finite gradients of magnitude `1e200` and volume weight `1e308`; require
   `ElasticOperatorNumericsError` before any result can be returned. Use two
   finite per-point local matrices near `1.3e308` each to force overflow only
   during global accumulation and require the same exception. Separately use
   finite `lambda` and `mu` whose diagonal sum `lambda+2*mu` overflows and
   require rejection of D. Never clip values or replace them with zero.

These are acceptance checks written in
`test_efgm_elastic_operator.py` for the operator primitive, not a crack
benchmark. All 11 passed in the one pinned focused invocation. The test source
does not claim MLS support stability or uniqueness of the six null modes.
The MLS adapter maps `MLSStatus.UNKNOWN_POOR_SUPPORT` to the specific
`UNKNOWN_POOR_SUPPORT` status. Other MLS UNKNOWN statuses map to generic
`UNKNOWN` and preserve the exact MLS reason code in `diagnostic`. The
shape-provider fixture must be independently affine-complete; it must not
reuse this assembly's B/D logic to compute its reference energy. Root must
review this provider name/status vocabulary before connecting the MLS weights
API. A later crack slice must show stiffness changes only through explicit
surface visibility and the resulting shape/support rows.

## Rigid gauge follow-up

The K null modes are diagnostic vectors, not sufficient gauge constraints. A
future static solve must impose three integral mean-translation constraints
and three integral mean-rotation constraints with Lagrange multipliers, using
physical `N_i(X)` and `grad N_i(X)`. A candidate independent definition is
`integral(u)dV=0` for translation and `integral(0.5 curl(u))dV=0` for rotation;
compare the exact source equations (6.7)-(6.8) before adopting it. No
coefficient pinning is allowed as a substitute for physical Dirichlet rows or
those gauges.
