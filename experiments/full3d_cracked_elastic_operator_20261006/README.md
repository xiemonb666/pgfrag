# Full 3-D cracked elastic operator source package

This package prepares the complete 3-D small-strain elastic stiffness for the
fixed straight-front Mode III manufactured-field case. The scalar 54-DOF
UZ bridge remains untouched and is not used as the mechanics operator here.
The new matrix contains UX, UY, and UZ coefficients for both BASE and
unshifted normalized FRONT_F1 at each of the 27 source sites, for 162 semantic
DOFs total.

## Operator and coefficient identity

The existing reviewed 3-D elastic primitive supplies engineering strain in
the order

```text
[eps_xx, eps_yy, eps_zz, gamma_xy, gamma_yz, gamma_xz]
```

and the full isotropic constitutive matrix built from
`lambda=E*nu/((1+nu)*(1-2*nu))` and `mu=E/(2*(1+nu))`. The resulting operator
is `K = integral(B.T @ D @ B) dV`. In particular, the `eps_zz` and normal
coupling terms remain in D; the operator does not substitute `mu*I`.

The 54 scalar rows are `phi_i` and `phi_i*Fbar`, where `Fbar=F1/sqrt(R)`.
The F1 gradients use the full product rule. Each scalar row is lifted to all
three displacement components, then the 162-by-162 matrix is permuted into
the frozen `DofRegistry.canonical_ids` order. The repeated source coordinates
used in the legacy assembler call are only coefficient-basis bookkeeping
centers; they do not describe 54 physical particles or change the 27-site
material geometry.

`build_full3d_registry`, `assemble_full3d_enriched_operator`, and
`rigid_mode_coefficient_vectors` are the public mechanics helpers in
`full3d_cracked_elastic_operator.py`. The rigid-mode columns are three unit
translations and three infinitesimal rotations. BASE entries use the source
coordinate displacement pattern required by MLS linear completeness; every
FRONT_F1 entry is zero.

## Fixed manufactured case and boundary load

The cylinder has `R=L=0.01 m`, the negative-x plane `y=0` is the crack cut,
and the straight front follows the z-axis. The 27 sources are the product
grid `{-0.005,0,0.005}^3`, with continuous unmasked linear MLS and Wendland C2
support radius `h=0.025 m`. All sources receive the global front enrichment;
there is no Heaviside term.

The benchmark field is `u=(0,0,A*sqrt(r)*sin(theta/2))`, with
`A=(KIII/mu)*sqrt(2/pi)`. Its nonzero stresses are
`sigma_xz=sigma_zx=mu*duz/dx` and
`sigma_yz=sigma_zy=mu*duz/dy`. The free outer-cylinder load is the complete
traction `sigma*n`. On the top cap it is
`(mu*duz/dx, mu*duz/dy, 0)`; on the bottom cap it is
`(-mu*duz/dx, -mu*duz/dy, 0)`. The full-vector boundary load diagnostic
integrates all three components. Both crack faces have zero traction by the
analytic Mode III solution and are not queried across the branch cut.

For the cap quadrature, `s` uses GL4 on `(0,1)`, `r=R*s^2`, and
`dA=2*R^2*s^3 ds dtheta`; this keeps 32 points per cap while resolving the
`1/sqrt(r)` traction in the rigid force and moment integrals. The diagnostic
checks both the integrated six-component boundary wrench and
`rigid_modes.T @ boundary_load`. Analytically, with `Q=mu*A*R^(3/2)`, the
side moment is `(4/3)*Q*L` about x; the top cap has force `(0,(4/3)*Q,0)` and
moment `(-(2/3)*Q*L,0,(4/5)*Q*R)`; the bottom cap has the opposite force and
moment `(-(2/3)*Q*L,0,-(4/5)*Q*R)`. Their total force and moment are zero. The
diagnostic records quadrature differences as-is and never projects or repairs
the load.

For exact reproduction, all BASE coefficients are zero; FRONT_F1 coefficients
are zero in UX/UY and `A*sqrt(R)` in UZ. `E=1 MPa`, `nu=0.3`, and
`KIII=100 Pa*sqrt(m)` produce the reference `1.3e-6 J` energy. These are
synthetic engineering benchmark values, not ceramic material data.

The radial interval still covers `[0,R]`. The positive GL4 nodes simply avoid
sampling the singular front. The strain grows as `r^-1/2`, so this is a
small-strain LEFM limit away from the mathematical front; no finite-strain or
tip-core model is included. The one fixed finite quadrature rule cannot
establish convergence.

## Prepared diagnostic and execution limits

`fixed_case.json` freezes the geometry, semantic DOFs, full cap traction,
Gauss rules, and limits. It predeclares all 384 unique locations: one 256-point
volume batch using GL4(r) x GL16(theta) x GL4(z), plus one 128-point boundary
batch containing 64 side points and 32 points per cap. Each provider batch is
at most 256 queries. The diagnostic in `full3d_operator_diagnostic.py` is
prepared for one worker attempt with a 60-second wall timeout and an 8 MiB
explicit numeric-array budget. The source estimate is 6,845,272 bytes (about
6.53 MiB), leaving 1,543,336 bytes below the cap; Python object overhead and
BLAS-internal workspace are excluded. A late timeout is recorded and never
retried.

The planned spectrum contains 162 eigenvalues: six BASE-only rigid modes and
156 positive modes are expected. Extra near-zero modes or significant
negative modes fail the diagnostic. No ridge, pseudoinverse, truncation, or
matrix repair is allowed. The code does not invoke a solver or apply the load
to a solver. An observed spectrum, energy match, or exact partition-of-unity
field is not mechanics qualification.

The worker saves its complete numerical arrays under
`static_solve_replay.schema = "static_solve_replay.v1"`. This includes the
162-by-162 semantic stiffness, 162-entry boundary load, 162-by-6 rigid-mode
coefficient matrix, registry canonical IDs and ordered semantic keys, source
coordinates and the source-centroid rotation center, the 384 query positions,
volume and boundary weights/normals/region slices, cached 384-by-54 basis
values and gradients, exact displacement/gradient/stress/status rows, and the
128-by-3 boundary traction. The solve bridge can consume these arrays without
re-running MLS or the operator. Rigid rotations are `omega x (x-centroid)` in
radians without division by `R`. JSON keys, units, shapes, query order, and
the lowercase semantic enum values are checked before a successful result can
be promoted.

The nested replay fields used by the bridge are `operator.K_N_per_m`,
`operator.boundary_load_N`, `operator.rigid_mode_coefficients`,
`registry.ordered_canonical_ids`, `registry.ordered_keys`,
`sources.positions_xyz_m`, `sources.rigid_mode_center_xyz_m`,
`basis_row_ids`, `quadrature.positions_xyz_m`,
`quadrature.volume_weights_m3`, `quadrature.boundary_normals`,
`quadrature.boundary_area_weights_m2`, `quadrature.region_slices_half_open`,
`basis.values`, `basis.gradients_per_m`, `exact_field.displacement_m`,
`exact_field.displacement_gradient_dimensionless`, `exact_field.statuses`,
and `boundary.traction_Pa`. The nested JSON holds roughly 123,000 numeric leaves;
using up to 35 bytes per leaf gives a conservative serialized estimate below
5 MiB before small metadata. The parent checks the actual stage file size
against 8 MiB before promotion.

The saved MLS/enriched basis gradients have units `1/m`; saved displacement
gradients `du_i/dx_j` are `m/m`, hence dimensionless. The source centroid is
the rigid-rotation origin, and the recorded rotations are not normalized by
the cylinder radius.

An assembly failure or a bad spectrum remains a stage diagnostic. The primary
path receives only a failure receipt with the stage path and SHA-256; only a
zero-return worker with the exact six/156 spectrum, complete finite replay
arrays, and matching source pins can be promoted. Even that pass keeps
`solve_allowed=false`.

The prepared command, reserved for the main agent after source review and
authorization, is:

```powershell
python -B pgfrag/experiments/full3d_cracked_elastic_operator_20261006/full3d_operator_diagnostic.py --output tmp/review_20261006/full3d_cracked_elastic_operator_20261006/run01/result.json
```

The pure schema fixture checks the fixed JSON structure and semantic wire
values without importing MLS or operator modules. Source preparation did not
run the provider, operator, solver, CAD, or Git. No numerical run result
exists yet.
