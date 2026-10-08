# Mixed93 objectivity design, adapter revision

## Review response

The independent source review identified two frame-sensitive defaults in the
frozen worker: branch functions were fed world x/y, and physical moments used
world x with an implicit zero center. The adapter now uses the actual
`CrackFrame` and local analytic field API; the moment helper receives x-c for
both its coefficient matrix and direct-field audit. The rotated center is
`c'=Qc+t`, so centered wrench moments rotate by Q and the world-origin moment
uses `Q M0+t cross F'`.

The public `mixed_enriched_operator` raw stiffness and force assemblers are
reused with its canonical 93 semantic IDs. The frozen worker's private
`_physical_moment_matrix`, `_direct_physical_moments`, `_scaled_spectrum_gate`,
`_solve_case`, and `_reconstruct_field` are reused behind a frame-aware row
builder. Every state forms fresh 31-scalar values and world gradients, raw K/C/f,
scaled spectrum and KKT solve. Both targets are g=0. No projection, ridge,
pseudoinverse, eigenvalue repair, or spectrum reuse is allowed.

## State map and checks

For Q=Rz(0.37)Rx(-0.23) and t=(0.003,-0.007,0.004) m, transform all source,
volume, and boundary points, crack origin/tip/frame, MLS origin, normals,
traction/stress/strain tensors, and vector coefficients. Keep measures, R/L,
E/nu, local KI/KII/KIII, reference scale, and support radius unchanged. Crack
local coordinates are xi=H^T(x-tip); branch values are evaluated in xi and
local gradients map to world by H. The same baseline xi rows feed analytic field
evaluation in both frames, avoiding ambiguous branch-side flips. A crack-cut
query without an explicit trace side blocks.

With A=I31 kron Q and B=diag(Q,Q), compare independently assembled
K'=A K A^T, C'=B C A^T, f'=A f, q'=Aq, and lambda'=B lambda. Convert engineering
strain shear by gamma/2 and stress shear directly when moving to tensors; rotate
T'=Q T Q^T, then compare volume reconstructed fields. Use the public raw
operator in 256-row chunks and separately recompute each fresh spectrum; require
6 zero, 87 positive, no significant negative modes before solving each 99x99
system. Units are K N/m, C dimensionless, q/g m, f/lambda N, strain unitless,
stress Pa, and work J.

The center wrench checks F'=QF, M_c'=Q M_c; the origin wrench checks
M0'=Q M0+t cross F'. Scale force residual by integral |traction| dA and moment
by R times that quantity, avoiding near-zero net resultants. The proposed capture
uses 1,280 MLS queries/state in five 256-row chunks, 2,560 total, 1,536 analytic
boundary field calls, two raw assemblies, two spectrum checks, and two solves.
The 30 s, 8 MiB, 64 MiB and numerical tolerances are provisional.

## Scope boundary

This first adapter checks coordinate covariance only. Co-rotating the same
source cloud and quadrature does not test relative cloud, quadrature-rule,
tetrahedral, or crack-front orientation bias. Those studies are excluded here;
quadrature sensitivity remains unresolved, and the moving-front antiplane
provider is not a mixed93 fracture solver. No numerical work has been run.
