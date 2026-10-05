# Bounded MLS-to-elastic operator diagnostic

This package connects the frozen 3-D MLS shape-function API to the frozen
small-strain elastic Galerkin operator for one fixed SI cube. It is source-only
until a separate review authorizes its one diagnostic invocation. It is not a
solver, crack benchmark, fracture simulation, or EFGM stability claim. The
broader explicit-surface/visibility formulation is governed by
`../explicit_crack_formulation_20261005/FORMULATION_CONTRACT.md`; that
contract records the published CPM/EFGM adoption and is provenance, not a
novelty claim for this adapter.

## Fixed case

`fixed_case.json` is the case record. It fixes the domain to `[-1,1]^3 m`,
places 27 MLS coefficient sites on the tensor grid `{-1,0,1}^3 m`, and uses a
single 4-by-4-by-4 tensor-product Gauss-Legendre volume rule (64 points,
volume 8 m^3). The isotropic small-strain material is E = 1e6 Pa and nu = 0.3.
The MLS reference origin is 0 m, reference scale is 1 m, and the fixed support
radius is 2.5 m. That radius is a numerical support length, not a material
parameter; it is not tuned by this diagnostic.

The identity-visibility path considers all candidate sources and passes no
visibility mask and no `mask_locally_constant` assertion. There is no crack
surface in this case. A future crack adapter must compute visibility from an
explicit surface and must establish query-local differentiability/stability;
it cannot reuse this identity path to claim crack coupling.

## Adapter contract

`MLSShapeProvider` makes one batched call to `shape_functions_3d` with source
positions `(27,3)` and the exact ordered quadrature positions `(64,3)`. MLS
weights `(64,27)` and analytic gradients `(64,27,3)` are returned in source
input order, which this private adapter maps directly to dense global IDs
`0..26`. Weights are dimensionless and gradient entries have units `1/m`.
Every point query from the operator must match the next frozen quadrature row.

An MLS `OK` row is exposed as operator `OK`. `UNKNOWN_POOR_SUPPORT` retains
that status; other MLS UNKNOWN statuses become generic operator UNKNOWN with
the exact originating reason, support count, and condition diagnostic. An
UNKNOWN row is never replaced with zeros or retried, and the operator must
return a blocked result without a partial stiffness matrix. The adapter stores
its coordinates and ID vector in bytes-backed immutable arrays so the cached
row-to-ID association cannot be made writable again.

For this fixed input the intended result orientation is
`Phi[q,i]` and `dPhi[q,i,j] = d(Phi_i)/dX_j`. The operator uses engineering
shear strain order `[eps_xx, eps_yy, eps_zz, gamma_xy, gamma_yz, gamma_xz]`,
where `gamma_ij = du_i/dX_j + du_j/dX_i`, with conjugate shear stresses
`[tau_xy, tau_yz, tau_xz]`. Thus D has shear diagonals mu. The assembled
stiffness has units N/m and `0.5 d.T K d` has units J.

MLS coefficients are not generally physical nodal displacements. This package
does not clamp coefficient entries to prescribe displacement. There is no
static solve or gauge fixing; the six translations/infinitesimal rotations are
only diagnostic vectors. Any future physical Dirichlet condition must use
shape rows `H(X) d = u(X)`, and a static solve must handle the six rigid-body
gauge modes independently.

## Single diagnostic plan

`operator_diagnostic.py` is the sole planned invocation. It estimates explicit
numeric array storage before MLS/operator work, rejects estimates over 32 MiB,
and runs the numerical work in a child process with a hard 30-second timeout.
The fixed input is much smaller than the operator's 375-DOF and 256-quadrature
caps. The estimate is not a bound on Python object overhead or vendor BLAS
workspace; those are stated separately in the JSON.

The JSON report is intended to include:

- MLS status, support count, and condition information for all 64 points;
- actual MLS affine interpolation/displacement-gradient errors if all MLS
  rows are `OK` (otherwise an explicit not-evaluated reason);
- the single assembled 81-by-81 K, symmetry, six rigid-mode residual, all 81
  eigenvalues, raw and beyond-roundoff negative counts, and extra near-zero
  modes after six;
- operator quadratic energy for an affine displacement coefficient field,
  compared with an independently evaluated isotropic tensor energy over the
  8 m^3 cube;
- input, timing, and source hashes.

Roundoff references are derived from binary64 unit roundoff, an operation
count covering two length-six constitutive contractions, quadrature
accumulation, and DOF-scale diagnostic products, and the observed finite MLS
condition number. They are
diagnostic labels, not pass/fail thresholds. The six rigid checks establish
only that those modes are present; any additional null modes remain possible.
The symmetric part of K is used only for eigenvalue reporting; K is not
symmetrized, regularized, clipped, or forced positive definite.

The runner applies no loads, invokes no solver, loads no crack geometry, and
does no SIF, crack-growth, calibration, timing benchmark, or rendering work.
The chosen material is an engineering fixture, not ceramic calibration. A
successful run would verify this bounded operator wiring and its stated
diagnostics only; it would not establish EFGM stability or actual fracture
capability.
