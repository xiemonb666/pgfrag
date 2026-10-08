# Bounded 3-D EFGM elastic operator

This folder holds a source-only small-strain isotropic elastic Galerkin
operator for the explicit-crack mechanics branch. It assembles a dense,
symmetric stiffness matrix from caller-supplied MLS shape values, analytic
gradients, global support IDs, and positive volume quadrature weights. The
hard limits are 375 displacement DOFs and 256 quadrature points.

The operator is intentionally independent of a particular MLS implementation.
An adapter implements `ShapeProvider.evaluate(material_position)` and returns
`ShapeSample(weights, gradients, support_ids, status)` as dense all-node rows
in canonical global order `0..N-1`; hidden or out-of-support contributions
are zero. Non-OK support status fails closed with no stiffness matrix. The
adapter can later obtain these values from the root MLS kernel after its
shape-weight API is added.

Finite inputs that overflow the constitutive matrix or a local/global
stiffness entry raise `ElasticOperatorNumericsError`; the operator returns no
partial or clipped matrix.

The strain vector uses engineering shear
`[eps_xx, eps_yy, eps_zz, gamma_xy, gamma_yz, gamma_xz]`, with
`gamma_ij = du_i/dX_j + du_j/dX_i`; the constitutive shear diagonals are `mu`.
Coordinates and quadrature volumes must use consistent length units, while
`E` uses the matching pressure unit. For the normalized formulation fixture,
`E=1`, `nu=0.3`, and `sigma=1e-3`; stress is not an input to stiffness assembly.

`physical_displacement_row` returns the physical field map `H(X)` for a
Dirichlet equation `H(X) d = u_bar`. It must not be replaced by setting an MLS
coefficient to a prescribed displacement. `rigid_mode_coefficient_vectors`
returns the six coefficient vectors used to check the three translations and
three infinitesimal rotations. It does not create the integral Lagrange
multiplier gauges needed by a static solve.

There is no crack state, finite-surface visibility, external load, boundary
integration, solver, SIF, or propagation code here. A future visibility-aware
shape provider can alter the support and therefore stiffness across a crack;
this module does not fabricate that coupling. This is not a full ceramic
shatter model or a reproduction of the published benchmark.

See `pgfrag/tmp/review_20261005/efgm_elastic_operator_20261005/` for the
provider contract, source review, freeze hashes, and focused test receipt. The
frozen test source passed once: 11 focused tests in 0.52 seconds from a pinned
snapshot. No solver or benchmark was run. Its global affine provider has only
four polynomial modes and can induce additional K null modes. The six-mode
check verifies the required rigid modes are present; it does not show they are
the only null modes or establish local MLS/EFGM support stability.
