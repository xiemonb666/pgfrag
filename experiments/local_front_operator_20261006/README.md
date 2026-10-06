# Local-front scalar antiplane operator diagnostic

This source package assembles one fixed 140-DOF scalar Poisson stiffness matrix from the pinned local-front MLS provider. It is an antiplane scalar bridge for the finite crack fixture, not full 3-D isotropic elasticity, a boundary-value solve, or fracture qualification.

The source fixture has 100 sites in the box `[-0.01, 0.01] × [-0.01, 0.01] × [-0.008, 0.008] m`, a finite crack on `y=0, x<0`, 100 BASE coefficients, 10 Heaviside coefficients, and 30 shifted-front coefficients. The provider and its inputs are read-only dependencies pinned by SHA-256 in `operator_diagnostic.py`. The local case and operator implementation are owned by this directory.

The assembled scalar form is

`K_ab = integral_Omega mu * grad(N_a) dot grad(N_b) dV`, with `mu = E/(2(1+nu))`.

All three scalar gradient components, including `dUz/dz`, carry coefficient `mu`. In full isotropic elasticity, the `epsilon_zz` contribution has coefficient `lambda + 2*mu`; this bridge does not implement that vector constitutive operator. Matrix entries have units N/m, coefficient values m, and quadratic energy J.

The volume rule splits x and y at zero into four open quadrants. Each quadrant uses tensor-product GL4 in x, y, and z, giving 64 positive-weight points per quadrant and 256 total. No query is clipped, nudged, or placed on the crack or front. The mapped weights sum to the box volume `6.4e-6 m^3`. GL4 exactly integrates one-dimensional polynomials through degree seven, but the MLS/enrichment integrand is not polynomial and the front gradients are singular; this single rule is a bounded diagnostic with no convergence or singular-integration claim.

The source lattice axes are x,z = `[-.008,-.004,0,.004,.008] m` and y = `[-.008,-.004,.004,.008] m`. The maximum coordinatewise nearest-site offsets over the box are `[.002,.004,.002] m`; the complete Cartesian lattice therefore gives a maximum nearest-site distance `0.004898979485566356 m`. With fixed support radius `h=.007 m`, the geometric margin is `0.002101020514433644 m`. This proves only that each point lies strictly inside at least one geometric support ball. The real per-query support counts, MLS rank/condition numbers, and UNKNOWN rows remain runtime facts.

The diagnostic reports the complete K and its symmetric-part spectrum, a BASE-only constant translation, a BASE affine `Uz=x` patch, independent scalar tensor energy, a normalized volume-mean row, canonical registry order, and the explicit numeric-array estimate. A significant negative eigenvalue or any nullity other than the single scalar translation leaves `assembly_status=OK` but makes mechanical eligibility `BLOCKED_INVALID_SPECTRUM` and `solve_allowed=false`. Such a complete diagnostic returns rc=3 and remains in staging with an UNKNOWN receipt, never the accepted primary. Only a spectrum with exactly the one expected translation null mode and no significant negative eigenvalues has `solve_allowed=true` for a later reviewed step; this diagnostic still does not invoke a solver. There is no clipping, pseudoinverse, ridge, or truncation. An UNKNOWN provider/MLS row blocks assembly and retains its original coordinate and reason.

The declared array estimate is 4,517,536 bytes (about 4.31 MiB) against an 8 MiB cap for explicit numeric arrays; Python objects and vendor eigensolver workspace are excluded. The hard run deadline is 60 seconds. The wrapper requires a new output path, captures the actual argv/return code/stdout/stderr, validates the complete JSON and hashes, and promotes only a complete rc=0 result. Nonzero, partial, or invalid results are kept as a staging file with a separate UNKNOWN receipt. No diagnostic was run while preparing this source package; execution remains subject to root review and authorization.

The fixed inputs are in `fixed_case.json`; the bounded implementation and no-overwrite worker wrapper are in `operator_diagnostic.py`.
