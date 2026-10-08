# Mode III scalar static-solve source review

**Status: source freeze ready for root review; no solver, tests, or numerical assembly run.** This is a restricted static traction solve layered over the frozen continuous MLS/front-PU provider. The source does not alter the provider, field reference, registry contract, operator case, or any production code.

## Mechanics and assembly

The experiment consumes the root-frozen operator case with `R=0.01 m`, `h=0.025 m`, `L=0.01 m`, `E=10^6 Pa`, `nu=0.3`, and `KIII=100 Pa sqrt(m)`. It builds the fixed 3×3×3 source grid, freezes 27 BASE and 27 unshifted FRONT_F1 UZ keys, and reuses the provider’s 54-column canonical rows. All-site front support (`T_F` all 27, `J_H` empty) remains an algebraic bridge; this does not establish local enrichment quality.

For each quadrature point it maps `grad(u_z)` through the pinned 6-component engineering-strain adapter and the pinned isotropic `D`. This yields the restricted scalar stiffness `K = ∫BᵀDB dV`, retaining the `lambda+2mu` axial-gradient term, and the side load `f = ∫side psi*t_z dS`. The analytic field reference supplies `sigma*n` on the 64 side points and both 32-point caps. It checks cap `t_z` is exactly zero, records all cap traction components, and omits cap x/y work from this UZ-only right-hand side. The branch cut remains the prescribed negative-x sheet; no front/branch-face query, clipping, radial hole, or epsilon shift is introduced.

The physical mean row is `c=(1/V)∫psi dV`; the only intended removed scalar mode is uniform UZ translation. The script equilibrates with `S=diag(K_aa^-1/2)`, forms the symmetric scaled 55×55 Lagrange system, and calls `numpy.linalg.solve` directly. It checks for exactly one near-zero eigenvalue, a small expected-translation residual, and unit physical mean for that translation before attempting a solve. An additional near-zero eigenvalue, missing translation, nonpositive diagonal, or unavailable gauge writes an UNKNOWN record and stops. The full spectrum and roundoff reference remain diagnostic only. No coefficient pin, pseudoinverse, ridge, penalty, or silent null deletion exists.

## Fixed quadrature and interpretation

The 128-to-256 volume comparison changes theta order only: GL4 radial × GL8 theta × GL4 z versus GL4 radial × GL16 theta × GL4 z. Radial and z orders stay four. This targets angular sensitivity from the half-angle front enrichment and rational MLS rows. Angular error can remain at GL16; two levels expose fixed-rule sensitivity but do not establish convergence, and no post-result threshold tuning or general qualification is allowed. Both levels share GL16 theta × GL4 z side loading (64 points), with GL4 radial × GL8 theta on each cap (32 points). An independent frozen 256-point set avoids the crack sheet and front.

The result records the analytic-coefficient defect `K*d_star-f` separately from the actual solve residual. For solved levels it also records the mean constraint, side net-force balance, analytic traction work, energy and external-work errors, and off-grid displacement/gradient max and RMS errors using fixed `A*sqrt(R)` and `A/sqrt(R)` scales. The classification remains diagnostic-only: no full-vector equilibrium, finite crack-end solution, material calibration, fracture, SIF, propagation, or local-enrichment convergence claim follows.

## Fail-closed behavior and scope

Before NumPy or numerical imports, the worker checks exact hashes for its case and immutable dependencies; it then preflights 125 source, 375 total DOF, 256-row query, 32 MiB explicit numeric-array, and 30-second parent-timeout limits. It uses two fixed volume assemblies, a direct solve, and one off-grid batch only. Any source mismatch, non-OK provider/reference row, nonfinite assembly, extra null, failed direct solve, or budget/deadline failure produces an UNKNOWN JSON before the worker raises its failure status. Existing outputs are not overwritten. Python object overhead and vendor BLAS workspace are excluded from the explicitly reported numeric-array estimate.

The script intentionally has not been launched. Root reported the preceding operator primary diagnostic hash as `8a62d61384156921fecc585b2336f161dd572988eadae8bed46afadf105c22c2` and noted that its scalar matrix had one translation mode and no additional null candidate; that evidence is context only and is not consumed by this source. Review this source and its hashes before deciding whether to authorize any numerical solve.

## Review outcome requested

No unresolved source-design choice is being escalated. The root should inspect the complete source freeze and decide whether the restricted mechanics assembly and fail-closed numerical gates are fit for one diagnostic run. No conclusion about physical or general method validity is made here.
