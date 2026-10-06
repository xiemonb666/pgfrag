# Source review: complete 3-D cracked elastic operator

Status: `SOURCE_READY_FOR_MAIN_AGENT_REVIEW`; no provider, operator, solver,
CAD, or Git command was executed for this package.

## Scope and fixed size

- 27 physical source sites on the frozen `3x3x3` grid.
- 54 scalar basis rows: source-major BASE and FRONT_F1 for each source.
- 162 vector coefficients: UX/UY/UZ for both basis families at all sources.
- One 256-point volume rule (GL4 x GL16 x GL4) and 128 outer-boundary points,
  for 384 predeclared unique locations split into provider batches of 256 and
  128 queries.
- No Heaviside enrichment, visibility mask, solver, crack growth, or material
  calibration.

## Operator/API path

1. `full3d_cracked_elastic_operator.py::build_full3d_registry(source_ids)`
   freezes semantic `DofKey(source, component, basis, index)` identities.
2. `mls_full3d_enriched_provider.py::ContinuousMLSFull3DFrontProvider`
   calls the existing `shape_functions_3d` once per fixed batch: 256 volume
   rows and 128 boundary rows. The two batches share the same 27 sources and
   162-key registry; they form BASE/F1 values and product-rule gradients.
3. `assemble_full3d_enriched_operator(...)` lifts those 54 scalar basis rows
   to 162 vector DOFs, calls the frozen `efgm_elastic_operator` assembly, and
   permutes its basis-major/component-minor matrix into the registry's
   canonical semantic order. The duplicated basis coordinates are only
   legacy assembler bookkeeping; the physical source set stays at 27.
4. `rigid_mode_coefficient_vectors(...)` builds six registry-aligned modes:
   three unit translations and three rotations `omega x (x-centroid)`. Only
   BASE coefficients carry values; all FRONT_F1 entries are zero.

The reused primitive builds the complete engineering-strain B and full
isotropic D from lambda and mu. Its convention is
`[eps_xx,eps_yy,eps_zz,gamma_xy,gamma_yz,gamma_xz]`; no scalar `mu*I`
constitutive substitute is present. For the manufactured field, the
reference has only `sigma_xz=mu*duz/dx` and `sigma_yz=mu*duz/dy` (plus symmetric
counterparts).

## Boundary traction check in the source

The lateral load uses the full reference traction `sigma*n`. The top and
bottom caps use respectively `(mu*duz/dx, mu*duz/dy, 0)` and
`(-mu*duz/dx, -mu*duz/dy, 0)`. The diagnostic integrates the three-component
boundary load vector, records each cap's x/y traction and resultant, and
computes exact and reconstructed work. The cap in-plane traction is retained
even though the exact manufactured displacement is purely UZ and `tz=0` on
both caps. The cap rule uses GL4 in `s` on `(0,1)`, with `r=R*s^2` and
`dA=2*R^2*s^3 ds dtheta`; volume quadrature stays at the frozen GL4 rule in
`r` over `[0,R]`. The transformed cap integrals make the singular traction's
rigid translation and rotation radial factors low-degree polynomials while
retaining 32 samples per cap. The branch faces are assigned analytic zero
traction without sampling an exact field value on the branch cut.

The source runtime plan checks both the six integrated force/moment values and
`rigid_vectors.T @ load_vector`. With `Q=mu*A*R^(3/2)`, the analytic side
moment is `Mx=4/3*Q*L`; top cap force is `Fy=4/3*Q` with moments
`Mx=-2/3*Q*L, Mz=4/5*Q*R`; bottom cap force is `Fy=-4/3*Q` with moments
`Mx=-2/3*Q*L, Mz=-4/5*Q*R`. Side force and the other side moments are zero.
The two caps and side sum to zero force and zero moment. Runtime records the
quadrature imbalance directly, with no load projection or equilibrium repair.

## Frozen benchmark and qualification limits

`fixed_case.json` keeps `R=L=0.01 m`, `E=1 MPa`, `nu=0.3`, and
`KIII=100 Pa*sqrt(m)`. The analytic energy is `1.3e-6 J`; all these values
are explicitly marked as a synthetic benchmark, not ceramic material data.
The front strain is singular like `r^-1/2`, so small-strain linear elasticity
is interpreted as a near-tip LEFM limit away from the front. Positive radial
nodes avoid the singular point but do not remove a hole from the domain. One
finite quadrature rule is not a convergence proof.

The future spectrum gate expects exactly six near-zero rigid eigenvalues and
156 positive eigenvalues. Extra near-zero modes and significant negative
modes fail the diagnostic without ridge, pseudoinverse, truncation, or matrix
repair. A spectrum or benchmark-energy match does not qualify mechanics.

The successful worker payload includes a `static_solve_replay.v1` block. Its
stable fields include the semantic `operator.K_N_per_m`,
`operator.boundary_load_N`, `operator.rigid_mode_coefficients` (162 by 6),
registry canonical IDs and ordered key records, source positions and centroid,
all 384 cached query/basis/exact-field rows, quadrature weights and boundary
normals, and the 128 full-vector boundary traction rows. Basis gradients are
in `1/m`; exact displacement gradients are `m/m` (dimensionless). A bridge can
reuse those rows without another MLS or operator evaluation. The pure schema
fixture checked fixed-case fields and actual lowercase enum wire values; it
did not import a numerical module.

The worker pins every numerical dependency before import. An assembly or bad
spectrum result remains a stage file; its parent path contains only a failure
receipt with the stage path and SHA-256. The parent promotes a result only for
worker return code zero, exact six/156 spectrum counts, finite complete arrays,
and exact source pins. Promotion leaves `solve_allowed=false`.

## Run envelope, not yet exercised

The diagnostic source pins frozen dependencies before importing numerical
modules. The planned single worker attempt has a 60-second parent timeout and
an 8,388,608-byte explicit-array cap. Its static estimate is 6,845,272 bytes
(about 6.53 MiB), leaving 1,543,336 bytes below the cap; Python object
overhead and vendor BLAS workspace are excluded. Timeout is a final recorded
failure with no retry. The complete replay JSON is estimated below 5 MiB and
the parent refuses promotion if its actual stage file exceeds 8 MiB. The
prepared invocation and exact case are documented in the package README and
`fixed_case.json`.

The entrypoint SHA-256 is
`e0f0cbbf0e60f65ce573ae8eed163e56c9b40a58be85b251445fd988e251936a`; all
upstream hashes are listed in `SOURCE_PINS.json` and are checked before any
numerical import. The fixed-case and semantic wire-schema fixture passed,
including the registry's actual lowercase enum values. The pure-stdlib fixture
also passed the complete rc0 promotion gate using finite zero-filled arrays of
the frozen shapes and a six-zero/156-positive synthetic spectrum; it performs
no operator or solver operation.

No numerical row, stiffness, traction vector, or spectrum has been generated
in this source phase. Main-agent authorization is required before the
prepared run is invoked.
