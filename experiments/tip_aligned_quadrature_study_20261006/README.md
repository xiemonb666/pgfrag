# Tip-aligned full-3D quadrature study

Status: `SOURCE_ONLY_AWAITING_ROOT_REVIEW`. This source does not modify the
frozen full-3D operator/provider, prior result, or static-solve bridge. The
source-check path uses only Python's standard library. Numerical imports,
quadrature generation, field/provider calls, assembly, and solves are confined
to the future `--execute` path and have not been run for this study.

The study has four cases: trial F1 center `a=0` or `a=0.5 mm`, each with
volume theta GL32 or GL64. All use the same physical cylinder `R=L=0.01 m`,
the 27 frozen source sites, continuous MLS, 162 semantic UX/UY/UZ DOFs, and
full isotropic small-strain `B^T D B` with lambda-mu coupling. The trial offset
changes the global F1 center and specifies a manually selected candidate
crack-front offset. The exact Mode III field, KIII, exterior traction, source
grid, and physical boundary remain fixed at the old tip; the old-tip
manufactured field is an exact comparison only for `a=0`. At `a=0.5 mm`, the
candidate crack surface extends to `x<a` and its true solution is unknown.
This is a fixed-load candidate reanalysis, with no committed crack state,
autonomous growth claim, or material qualification. No bulk quadrature rows
lie on the zero-measure crack surface.

For each volume case, integration is centered at `(a,0)` in the disk plane:

```text
rho_max(theta,a) = -a*cos(theta) + sqrt(R^2 - a^2*sin(theta)^2)
rho = rho_max*s^2
dV = 2*rho_max^2*s^3 ds dtheta dz
```

The rule is GL4 in `s`, GL4 in `z`, and GL32 or GL64 in theta over `(-pi,pi)`.
It covers the full disk with open Gauss nodes and creates no radial hole.
The four volume clouds contain 3,072 requested rows total. Assembly uses the
unchanged frozen assembler in deterministic chunks of at most 256 rows and
sums raw 162x162 stiffness matrices; no frozen point limits are raised.

The side rule is fixed at GL64(theta) x GL4(z), 256 rows. The cap rule uses
one shared physical point cloud for both trial offsets, so side and cap
traction discretization is fixed as `a` changes. On each cap, theta is split
into `(-pi,0)` and `(0,pi)`, with GL8 and GL16 on each half for coarse/fine
load diagnostics. The radial coordinate is split at `c=0.5 mm` for both
offsets:

```text
inner: r = c*sin(pi*s/2)^2, dr/ds = (pi*c/2)*sin(pi*s)
outer: r = c + (R-c)*s^2, dr/ds = 2*(R-c)*s
dA = r*(dr/ds)*ds*dtheta
```

GL4 is used on each radial segment. Thus each cap has 128 coarse or 256 fine
points; both caps have 256 or 512. The physical cap points, normals, area
weights, and old-tip exact traction are shared by both `a` values. Each case
records the coarse/fine cap generalized load, resultant, moment, and work;
the solve always uses the fine GL16-per-semicircle cap rule. Cap `tx` and `ty`
are retained from the full `sigma*n`. The split focuses quadrature around the
old traction singularity and the shifted F1 singularity radius, but finite
GL4/GL8/GL16 rules do not certify integration error.

The declared maximum is 4,992 query rows: 3,072 volume, 256 side, 1,536 cap
rows counting each trial-offset load evaluation, and 128 shared off-grid
rows. The raw group plan is at most 20 chunks of 256. All physical xyz rows
are exact-deduplicated before the first MLS/reference call; actual requested,
unique, and provider-batch counts are recorded separately. Each provider call
stays at or below 256 rows. The prior 384-location replay is not used as a
cache. The off-grid rule has four old-tip radii, eight angles away from the
branch cut, and four z levels.

The fixed gauge is the analytic old-tip volume moment
`g=[0,0,0,A*sqrt(R)*16/(21*pi),0,0] m` with
`g_rx=5.031114771170333e-6 m` for all four cases. This differs from the old
finite-theta saved value `5.031213792347516e-6 m` by relative
`1.9681756765e-5`; it changes only the numerical gauge target, not KIII or
traction. The six constraints are physical volume moments computed from each
case's own integration points. Every stiffness must pass the six-rigid/156
positive spectrum gate before the direct six-constraint saddle solve. Extra
zero or significant negative modes fail closed; there is no ridge,
pseudoinverse, coefficient pinning, load projection, or spectral truncation.
The scaled raw stiffness must be finite and have relative asymmetry at most
`1e-12`, measured as `||Kscaled-Kscaled^T||_F/max(||Kscaled||_F,tiny)`; its
eigenvalues must also be finite. The engineering zero/negative
threshold is `64*eps_float64*162*spectral_scale`, and the rigid-mode relative
residual limit is `128*eps_float64*162`. These numerical thresholds are not a
strict eigenvalue proof. MLS condition numbers are recorded and governed by the
provider's own qualification limit; they never relax mechanical thresholds.

The solve also records and gates the weak-equilibrium, gauge, saddle backward,
and constrained-energy identity residuals at `1e-9`; the reconstructed-strain
versus raw-`K` energy defect has the same `1e-9` limit. These are fixed-front
algebraic consistency checks, not physical-accuracy or fracture acceptance
criteria. The single-attempt output defaults to
`pgfrag/tmp/review_20261006/tip_aligned_quadrature_study_20261006/run01` and
must remain under this study's own `pgfrag/tmp` review subtree.

At `a=0`, weighted displacement, gradient, full-stress, and spurious-normal
stress errors are compared with the old-tip manufactured field over the
volume and shared off-grid rows. At `a=0.5 mm`, report the solved-field change
from the `a=0` fixed-load solution on common physical points and the numerical normal to
shear stress ratio; do not label it as an exact-new-tip error. Report
`U=0.5*d^T*K*d`, `W=d^T*f`, `Pi=U-W`, the constrained identity
`2U-W=-g^T*lambda`, six rigid virtual works, boundary forces/moments, and
KTT-versus-reconstructed-strain energy consistency.

For each volume theta order, the finite offset value is

```text
G_delta_a = -[Pi(a=0.5 mm)-Pi(a=0)]/(L*0.5 mm)
```

The prior `0.005466 J/m^2` result is context. The exact `0.013 J/m^2` value
is an old-tip local/asymptotic reference; it is not a percent-error criterion
for finite `delta-a` energy differences. No J-contour points or J-integral
path are added in this first study.

Future execution limits are one local attempt capped at 60 seconds, 32 MiB of
explicit numeric arrays, and 20 MiB combined NPZ/JSON output. The NPZ stores
numerical arrays; JSON stores metric summaries and the NPZ hash. The JSON does
not contain bulk arrays. This is a finite quadrature comparison, not a proof
of convergence. Small-strain LEFM remains a near-tip idealization without a
finite-strain/core model, and the synthetic E, nu, KIII, and energies are not
ceramic material calibration.
