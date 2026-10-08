# Mode III scalar static-solve v2 source review

**Status: distinct output-evidence source freeze; no solver, tests, or numerical assembly run.** V1 at `pgfrag/experiments/modeiii_pu_static_solve_20261005/` is preserved byte-for-byte. This bundle changes the result evidence and serialization safeguards only. It does not change the coefficient registry, MLS/front-PU provider, quadrature, stiffness/load/mean assembly, diagonal scaling, spectrum tolerances, saddle system, or direct solve.

## Output additions

The result publishes the ordered canonical coefficient IDs and states that they index all vectors and both axes of each raw stiffness matrix. Each volume level retains the complete 54×54 `K` matrix in N/m and its normalized physical-mean row in dimensionless form. The applied 54-entry scalar load vector is written once in N. Existing coefficient vectors and weak defects remain present, allowing independent reconstruction of `K*d + c*lambda - f`, gauge residuals, energy, and work from saved data.

The off-grid record now includes all 256 query coordinates in m, exact `u_z` in m, and exact `grad(u_z)` as a 256×3 dimensionless strain-like array. Each quadrature level also stores the solved `u_z` and solved gradient on those same coordinates. Basis gradients remain in 1/m; multiplying by coefficient displacements in m gives dimensionless `grad(u_z)`. V2 renames the prior mislabelled absolute gradient-error unit from per-meter to dimensionless. The prior v1 key remains unchanged in its immutable snapshot.

## Finite and output-size guards

Before successful result serialization, the worker checks the retained matrices, mean/load rows, analytic coefficient defect, off-grid query/reference/solution/error arrays, coefficients, energy, external work, weak residual, gauge values, and energy/work errors for finiteness. `_json_safe` now raises on any nonfinite float instead of mapping it to `null`. A serialization failure produces a compact UNKNOWN receipt before the worker returns failure. Result JSON has a hard 4 MiB cap; an oversized candidate writes `UNKNOWN_RESULT_JSON_BUDGET` without promoting the large evidence payload.

The explicit numeric-array estimate adds the retained evidence arrays using fixed counts: two 54×54 matrices, three 54-vectors (two mean rows and one load), and 256×(3 coordinate + 4 reference + 2×4 solved values/gradient components) float64 entries. This is approximately 78.7 kB of added numeric payload, in addition to the existing conservative 32 MiB bounded estimate; Python list overhead remains excluded. The worker enforces the existing 30-second parent deadline and keeps all prior fail-closed UNKNOWN behavior.

## Preserved numerical scope

The fixed case remains `R=0.01 m`, `h=0.025 m`, `L=0.01 m`, `E=10^6 Pa`, `nu=0.3`, `KIII=100 Pa sqrt(m)`, and the 27-source grid with 54 scalar UZ coefficients. The 128-to-256 change is GL8-to-GL16 theta only; radial and axial orders remain GL4. It diagnoses angular sensitivity for the half-angle front term and rational MLS rows without claiming convergence or general qualification. Full cap tractions are still recorded, scalar cap `t_z` must remain zero, and x/y cap work stays outside the restricted UZ RHS.

V1 source SHA-256 is preserved as context: solver `D9294D9A90D20505D6ECD01354A2D9E0897F56B5FF0BB98AAD01A4FD31513787`, case `93C51463D717ADAA8E78BC163FB8915AAF26164F411C2B46A20D80492BDA42D4`. No original v1 source, case, review, or freeze was overwritten. Root must review this distinct source freeze before authorizing a solve.
