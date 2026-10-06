# Tip-aligned quadrature study source review

Status: `SOURCE_ONLY_AWAITING_ROOT_REVIEW`; numerical execution is explicitly
unauthorized in the case file. No quadrature points, MLS/reference calls,
operator assembly, eigensolve, or saddle solve were run for this study.

The prepared entrypoint builds the four fixed `(a, theta_order)` cases,
deduplicates declared xyz rows before one MLS/reference pass split into
<=256-row calls,
reconstructs shifted normalized F1 rows with the product rule, and sums raw
full-vector stiffness matrices through the frozen <=256-row assembler. It
preflights the six-rigid/156-positive spectrum and solves only when that gate
passes. Its cached-row adapter exposes the frozen source-major scalar
`basis_row_ids`, matching the assembler contract without querying the provider
again. The solve gates weak equilibrium, gauge residual, saddle backward
error, and constrained-energy identity at `1e-9`; reconstructed-strain versus
raw-K energy has the same limit. These are algebraic consistency checks, not
physical-accuracy or fracture-acceptance criteria. The reviewable outputs include raw K, load/gauge/rigid arrays,
registry IDs, quadrature maps/weights/normals, base MLS/reference arrays,
per-case volume/boundary/off-grid basis arrays, and solved metrics in a
bounded NPZ plus a small JSON summary.

Mechanical spectrum checks are independent of MLS conditioning: non-finite
scaled matrices/eigenvalues fail closed, scaled relative asymmetry must be at
most `1e-12` using `||Kscaled-Kscaled^T||_F/max(||Kscaled||_F,tiny)`,
zero/negative classification uses
`64*eps_float64*162*spectral_scale`, and the relative rigid residual gate is
`128*eps_float64*162`. These are engineering numerical tolerances, not a
strict eigenvalue proof; MLS condition numbers remain subject to the provider's
own qualification limit.

Both candidate offsets share exactly the same cap coordinates, normals,
weights, and old-tip full-stress traction. The cap radius is split at
`c=0.0005 m` using the frozen case's two-level GL8/GL16 semicircle rules and
GL4 on each radial segment. Each case solves with the fine GL16 cap load.
The full requested-row upper bound is 4,992 and the pre-dedup batch bound is
20 batches of at most 256 rows. The volume rules are GL4(radial) x
GL32/GL64(theta) x GL4(z), with `rho=rho_max*s^2` centered at the trial tip.

This is a manually specified fixed-load candidate reanalysis. The old-tip
manufactured field is exact only for `a=0`; the true shifted candidate
solution is unknown. The source makes no committed crack-state, autonomous
growth, convergence-proof, or material-qualification claim. It adds no J
contour queries. The existing full-3D operator/provider/result and previous
receipts remain pinned read-only inputs.

The pure source check passed with the bundled Python runtime:

```text
SOURCE_CHECK_OK_NO_NUMERICAL_IMPORTS
requested_rows_max=4992
raw_batch_upper_bound=20
fixed_shared_cap_split_m=0.0005
cached_row_basis_contract=AST_VERIFIED
future_execution=NOT_RUN
```

It used Python's standard library only. No Git command was run. Before any
future numerical execution, root review and separate explicit authorization
are required; the `study_case.json` gate currently blocks execution.
The future output path is constrained to this study's own
`pgfrag/tmp/review_20261006/tip_aligned_quadrature_study_20261006` subtree,
with default `run01`; existing workspace-root `tmp/review` artifacts remain
source-review records and are not output targets.
