# Mode III continuous-MLS PU scalar operator

This source-only package defines the first enriched scalar `u_z` operator
bridge for the frozen straight-front Mode III manufactured field. It uses the
frozen `DofRegistry`/antiplane provider contract, phase-2 continuous MLS
shape values and gradients, the frozen isotropic constitutive matrix, and the
root-reviewed analytic field reference. Shared primitive files remain
read-only. No test, solver, or numerical evaluation is part of this source
freeze.

## Fixed construction

The case is a cylinder `r <= R = 0.01 m`, `z in [-L/2,L/2]`, with `L=0.01 m`
and a straight z-axis front. The crack cut is the negative-x radial half-plane
`y=0`, `x<0`. The source grid is the 27-point product of
`{-0.005,0,0.005} m`; its `0.005 m` grid spacing is separate from the fixed
MLS support radius `h=0.025 m`. MLS is continuous with Wendland C2, origin
`(0,0,0) m`, reference scale `R=0.01 m`, and no visibility mask.

Every source has exactly two coefficients in the root-frozen registry: UZ
`BASE` and UZ `FRONT_F1_UNSHIFTED`, for 54 coefficients total. The unshifted
front function is normalized as

```text
Fbar = F1 / sqrt(R)                         dimensionless
psi_i = phi_i Fbar
grad(psi_i) = Fbar grad(phi_i) + phi_i grad(Fbar)
```

Both coefficient families therefore have units metres. The manufactured
coefficient vector sets all BASE entries to zero and every front entry to
`A sqrt(R) m`, where `A=(KIII/mu)sqrt(2/pi)` and `KIII=100 Pa sqrt(m)`.
No source value of `F1` is needed, so sources on the crack cut and front need
no invented side value. Volume and boundary queries avoid the front and exact
branch cut through positive radial and interior Gauss nodes; the adapter has
no epsilon-offset or masking path.

All 27 sites belong to the front set `T_F`, because each fixed support ball
contains the complete straight front segment. `J_H` is empty because these
same supports contain the front and do not meet the root's front-exclusion
condition for Heaviside enrichment. This all-site choice is only a global
algebraic bridge for this small manufactured field. It does not demonstrate
local enrichment qualification, a front transition, or convergence.

## Restricted scalar stiffness and quadrature

The scalar column for a UZ shape gradient is mapped through the existing
engineering-strain convention as
`[0,0,dphi/dz,0,dphi/dy,dphi/dx]`. The assembled scalar bilinear form is

```text
K_ij = integral [mu (d_i,x d_j,x + d_i,y d_j,y)
                 + (lambda + 2 mu) d_i,z d_j,z] dV
```

The longitudinal term is kept because MLS products may vary in z even though
the exact Mode III field does not. The implementation uses the frozen
isotropic `D` and antiplane strain mapping; it does not pass 54 coefficients
through the old 3-DOF-per-node assembler.

The frozen one-shot quadrature has 128 volume points (`r`: GL4 over the full
`[0,R]` interval with positive nodes and cylindrical Jacobian `r`; `theta`:
GL8 over the open branch interval `(-pi,pi)`; `z`: GL4 over `[-L/2,L/2]`), 64
lateral points (theta GL16, z GL4), and 32 points on each end cap (r GL4,
theta GL8). The 256 unique positions are submitted to one MLS batch in the
order volume, lateral side, top cap, bottom cap. `r_min` records the smallest
positive sample only; the radial interval still covers all of `[0,R]`, with
no hole, cutoff, or clipping.

The planned JSON diagnostic computes the 54-by-54 scalar K once, coefficient
field and gradient errors against the analytic reference, independent exact
tensor energy, traction work, symmetry, the one representable scalar UZ
translation residual, all eigenvalues of the symmetric part, negative values,
and additional near-zero values after that scalar rigid translation. The
roundoff scales are descriptive only. Six full-vector rigid modes cannot be
represented in this UZ-only coefficient space.

The external side traction contributes to scalar UZ virtual work. On end caps,
the scalar traction is `t_z=0`; full vector x/y end-cap tractions are recorded
as reference boundary data but are excluded from the restricted scalar
virtual work. No traction is applied to a solver because there is no solver.
The restricted scalar model cannot establish full 3-D equilibrium.

## Qualification boundary

An exact partition-of-unity field, small energy discrepancy, symmetric K,
small scalar translation residual, or one observed spectrum is insufficient
to qualify the mechanics field. Extra null modes and negative eigenvalues are
reported without pseudoinverse, ridge, truncation, or positive-definiteness
enforcement. This remains a fixed uncracked operator diagnostic, not an
EFGM stability result, finite-front fracture result, SIF validation, solver
coupling, or material calibration.
