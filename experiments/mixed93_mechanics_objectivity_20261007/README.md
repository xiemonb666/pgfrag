# Mixed93 frame-aware mechanics objectivity adapter

Status: source only. `capture_arrays` is guarded by `numeric_execution_authorized`,
`execution_authorization.authorized`, and an exact root-go token. This revision
has all authorization fields false/null. A future runner must verify the case,
all source pins, and the cached archive digest before opening the NPZ; it then
passes the six arrays and canonical source IDs to `capture_arrays`. This module
does not open or load the archive.

The adapter reuses `mixed_enriched_operator` public raw K/f assembly and the
pinned cached worker's private physical-moment, fresh spectrum, KKT solve, and
field-reconstruction helpers. It builds every 31-scalar semantic row itself.
For each query, it evaluates the four branch functions at
`xi=H^T(x-tip)` and maps the local gradient to world coordinates with H. BASE
MLS values/gradients come from the actual `shape_functions_3d` API. The rotated
state transforms source/query points and MLS reference origin, while scalar MLS
lengths and all physical measures/material/load scalars stay fixed.

The analytic traction uses `evaluate_local_point` with the co-rotated
`CrackFrame`, which returns world displacement/tensors. Thus exact branch-side
selection uses the same local coordinates for both states. An exact crack-cut
row without a pinned trace side blocks. Engineering strain shear components are
halved when making a tensor; stress shear components are not. Both tensors are
rotated before any Voigt comparison.

The gauge helper receives `x-c` in both states, with `c'=Qc+t`; both targets
are zero. Centered moment vectors rotate by Q. The returned world-origin wrench
also checks `M0'=Q M0+t cross F'`, with force normalized by integral absolute
traction and moment by R times that integral. The record includes both
geometric frames/wrenches and raw K/C/f/q/lambda and volume displacement,
engineering strain, and stress arrays. Each state independently assembles
fresh K/C/f, recomputes the spectrum, and solves its own 99x99 KKT system.
There is no spectral reuse or matrix repair. Volume field comparisons reuse
the 512 queried volume rows.

The fixed-BASE provider check remains limited to 16 interface queries. This
adapter proposes 2,560 real MLS queries, 1,536 analytic boundary evaluations,
two fresh raw assemblies, two spectrum checks, and two solves. Budgets and
error gates remain provisional until ROOT reviews this source and authorizes a
single capture. Its classification is coordinate covariance of the same toy
mechanics representation. Relative source-cloud bias, quadrature sensitivity,
tetrahedral bias, crack growth, and physical material qualification remain
unresolved; the actual tetra/front experiment requires a separate solver path.
