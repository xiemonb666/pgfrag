# Elastic operator source contract

Status: source-only, not executed. This module is a bounded EFGM linear
elasticity primitive. It adopts explicit crack surfaces as the future
mechanical state, while keeping this first operator independent of crack
visibility and fracture evolution. It is not the existing OSB/CriticalStretch
mechanics and does not claim a full ceramic shatter formulation.

## Input and output

The public assembly call is:

```text
assemble_elastic_operator(
    coefficient_positions,
    quadrature_positions,
    quadrature_volume_weights,
    IsotropicMaterial(E, nu),
    ShapeProvider,
) -> ElasticAssemblyResult
```

Each point from `ShapeProvider.evaluate(X)` returns a `ShapeSample` containing
plain NumPy arrays `weights[N]`, `gradients[N,3]`, and integer
`support_ids[N]`, plus a `ShapeStatus`. This contract uses dense all-node rows
in canonical order: IDs must equal `0..N-1`; in an `OK` sample,
hidden/out-of-support values and gradients are zero. Gradient column `j` is `d(N_i)/dX_j` in the
reference/material frame. The shape-value row is retained for physical
boundary constraints; only analytic gradients enter strain assembly. A
generic `UNKNOWN` must carry the provider's original reason/status in its
diagnostic field; `UNKNOWN_POOR_SUPPORT` retains its specific status.

The accepted hard limits are 125 coefficient sites / 375 vector DOFs and 256
quadrature points. Input type, dimensions, row caps, real dtype, and finite
values are checked before coordinate conversion copies or dense stiffness
allocation. Provider support count is capped by the already-validated site
count before provider arrays are copied. Dense output therefore has at most
375 by 375 entries. Positive quadrature volume weights are required.

An `UNKNOWN_POOR_SUPPORT` or `UNKNOWN` provider result returns
`BLOCKED_UNKNOWN_SUPPORT` with `stiffness=None`. A provider `BLOCKED` result
returns `BLOCKED_SHAPE_STATUS`. Any non-OK quadrature stops assembly before
the global matrix is allocated; no failed point is replaced with zeros. A
malformed `OK` sample raises
`ShapeProviderContractError`. A non-finite constitutive matrix, local matrix,
or accumulated global matrix raises `ElasticOperatorNumericsError`; there is
no clipping or `nan_to_num` fallback. No failure path returns a partial
stiffness matrix.

## Weak form and convention

For each quadrature point, the module assembles

```text
K += B.T @ D @ B * volume_weight
```

with `B` using node-major vector coefficients and the engineering-shear
strain order `[eps_xx, eps_yy, eps_zz, gamma_xy, gamma_yz, gamma_xz]`, where
`gamma_ij = du_i/dX_j + du_j/dX_i`. The stress vector is
`[sig_xx, sig_yy, sig_zz, tau_xy, tau_yz, tau_xz]`, so each shear diagonal in
`D` is `mu = E/[2(1+nu)]`. The normal block uses the Lamé coefficient
`lambda = E nu/[(1+nu)(1-2nu)]`. Material validation requires finite `E>0` and
`-1<nu<0.5`.

With consistent length and pressure units, stiffness units are pressure times
length, and `0.5*d.T@K@d` has energy units. `sigma=1e-3` is a normalized load
scale for the later mechanics fixture and is not consumed here.

## Physical boundary rows and rigid modes

MLS coefficient entries are not generally physical displacements at the
coefficient sites. A prescribed displacement at physical point `X` must be
represented as `H(X) d = u_bar`, with `H` assembled from the provider's shape
values by `physical_displacement_row`. This module does not clamp coefficient
entries or solve constrained equations.

On a connected, linearly complete free body, the stiffness should have six
infinitesimal rigid modes: three uniform translations and three rotations
`u_i = omega cross (X_i - X_centroid)` in the coefficient vector. The helper
`rigid_mode_coefficient_vectors` supplies those six columns for a later null
mode check. A static solve still needs six gauge constraints. The formulation
contract calls for integral Lagrange-multiplier constraints on mean
translations and mean rotations; those constraints must be formed from the
physical field shape rows and gradients, not from pinned MLS coefficients.
The exact implementation should be checked against the source's Eqs. (6.7)
and (6.8) before any solver is added. The current source intentionally has no
constraint matrix or solver.

## Crack and benchmark boundary

The next adapter can ask the explicit surface/visibility layer for one
side-specific support and pass its shape rows into this assembly. Zeroing
opposite-side contributions then changes `B` and hence `K` through the weak form.
This source has no crack, visibility query, support split, surface state,
cohesive law, contact, growth, or SIF path. The formulation contract's
inclined penny-crack example remains a future benchmark, not a current result.

The first frozen mathematical fixture is only 27 coefficient sites on a
`3x3x3` grid over a side-2 cube and 64 volume points from `2x2x2` background
cells with `2x2x2` Gauss points per cell. Its material values are `E=1`,
`nu=0.3`; a later load path may use `sigma=1e-3`. It is a primitive check and
not a published full-resolution benchmark.
