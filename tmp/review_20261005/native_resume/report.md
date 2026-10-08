# PG-Frag native dynamics and bond-break continuation review

**Review date:** 2026-10-05  
**Status:** Existing B native threshold run completed; read-only audit passed. No solver was rerun.  
**Scope:** Continue the 2026-10-04 A/B work by checking receipts, all pinned sources, every directed/canonical B bond row, the native failure path, and both scientific images. Keep the old receipts, raw CSVs, design notes, and trial records unchanged.

## What the run actually did

The B run was not left in an interrupted state. Its launch receipt says `COMPLETED`, wrapper exit `0`, child exit `0`, child exit observed, wait completed, no timeout, no launch error, and 4.4945669 seconds elapsed. The original stdout contains `PASS_NATIVE_BRANCH_DIAGNOSTIC_ONLY`; stderr is empty. There was no Julia process at the start of this continuation audit. The existing files are sufficient, so this review did not launch Julia or rerun the native case.

The frozen manifest and design note still carry their pre-run labels (`A_NARROW_GATE_PASS__B_NOT_RUN` / `FROZEN_PENDING_SINGLE_RUN`). They remain unmodified as history. The completed receipt and raw snapshot, independently checked below, establish B's observed status; this report supplements the older design record.

The pinned driver built three fresh native OSB handlers from the same 125-point reference and applied `F = diag(1,1,1+η)` to every point for `η/εc = 0, 0.5, 1.5`. For each branch it called `Peridynamics.calc_force_density!` once, then exported the native before/after activity, native bond lengths, stretch, and failure-permission fields. It did not call VelocityVerlet integration and did not manually edit `bond_active`. The correct conclusion is a native CriticalStretch branch diagnostic under an imposed affine geometry, not a time-evolving fracture simulation.

## Read-only checks and results

The local assertion audit is [audit_native_run.py](../../../experiments/native_resume_20261005/audit_native_run.py). It passed against the immutable old receipt and raw files. It rehashed all 26 sources in the receipt both before/after lists against current files, including byte counts; all 26 matched. The audit also matched the current Julia executable, B driver, manifest, preflight receipt, stdout/stderr, and all three output artifacts to their recorded hashes. Full pin details are in [native_resume_audit.json](native_resume_audit.json).

The original [directed bond CSV](../../review_20261004/pd_loaded_diagnostics/snapshots/affine_threshold_b_run01/affine_native_bond_states.csv) is 2,429,849 bytes, SHA-256 `e55b16f61a296ae092ce4e1c01166f7a5eff49735338b19e6364d25178a4af8d`; the [canonical pair CSV](../../review_20261004/pd_loaded_diagnostics/snapshots/affine_threshold_b_run01/affine_canonical_pair_states.csv) is 1,158,067 bytes, SHA-256 `b4c5ae70126bb9997f28b4a766c4ae7986360246e335d6f904af0e8cafafd140`. The [native summary](../../review_20261004/pd_loaded_diagnostics/snapshots/affine_threshold_b_run01/summary.toml) is 1,497 bytes, SHA-256 `5a9f3b53eac2fe122ac820e7e43c5239d4fae490fcc4c73aa42fbd9af1963efc`. The immutable [launch receipt](../../review_20261004/pd_loaded_diagnostics/receipts/affine_threshold_b_run01.launch_receipt.json) SHA-256 is `8ab04d7c70fc5a35acac70ed0907a4939a238af83e4f949435946a06bab53e35`; preflight receipt SHA-256 is `4a7738f0a2440ad392ad01b50354fc8d7f0040f7fffb3da2992dde300e81cc1d`.

The point fixture SHA is `9b9207727364c667ecc0b592ccff34405eaab75a0992e057a2d8a7a121b706da`; it contains 125 cell centers in a 5×5×5, 1 mm-spaced grid, each with a 1 mm³ quadrature volume. The material values are `E=5 GPa`, `ν=0.24`, `ρ=2100 kg/m³`, `Gc=45 J/m²`; the declared horizon is `δ=3 mm`. Using the pinned OSB bulk modulus `K = E/[3(1−2ν)]`, the source CriticalStretch formula `εc = sqrt(5 Gc/(9 K δ))` gives `0.0016124515496597099`, exactly matching the native summary within floating precision. These remain `UNVALIDATED_PROTOTYPE` values.

The raw directed CSV has 18,372 data rows, or 6,124 directed bonds in each of the three states. The canonical CSV has 9,186 data rows, or 3,062 pairs per state. The assertion audit checked each row's unique directed endpoint key, canonical endpoint IDs, finite positive lengths, finite stretch and threshold, fixture SHA, reconstructed reference length, and reconstructed affine current length. For every row, its recorded stretch agreed with `(current_length − reference_length)/reference_length` and the independent point-coordinate calculation. Every canonical pair had exactly two reversed directed rows; lengths, stretch, epsilon, permission, and activity were identical across directions. Maximum forward/reverse stretch difference was zero.

Every row began active and all 6,124 directed rows per state were failure-permitted. Per directed row, `active_after` and `transition` agreed with the native rule `stretch > εc && fail_permit`. The independent check found no threshold or reciprocity mismatches:

| Prescribed state | Directed rows | Canonical pairs | Stretch/εc range | Broken canonical pairs | Broken directed rows |
|---|---:|---:|---:|---:|---:|
| rest (`η/εc=0`) | 6,124 | 3,062 | 0 to 0 | 0 | 0 |
| below (`η/εc=0.5`) | 6,124 | 3,062 | 0 to 0.5000000000000498 | 0 | 0 |
| above (`η/εc=1.5`) | 6,124 | 3,062 | 0 to 1.500000000000283 | 632 | 1,264 |

For the 632 broken pairs, the smallest recorded stretch was `1.0004024638033726 εc`; the greatest retained stretch in that above-threshold branch was `0.7504529540581504 εc`. Thus the strict native comparison separates all observed transitions from all retained rows. B's narrow gate passed as recorded: `PASS_NATIVE_BRANCH_DIAGNOSTIC_ONLY`.

The native source confirms the audited rule. `calc_failure!` computes `ε = (l-L)/L` and sets a permitted bond inactive only when `ε > εc` in [ordinary_state_based.jl](../../../src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/physics/ordinary_state_based.jl:138). The bond-system force path calls native failure, damage, then force evaluation in [bond_system.jl](../../../src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/discretization/bond_system.jl:264). The B driver's SHA is `2ea0ce90831ec1a9d299a02fc7a223b8bbff1cd411991ce80a0ddbf40e408c2c`; all source hashes and raw artifact hashes are indexed in the JSON audit.

## Image review

The existing A image at [loaded_elastic_a8_science.png](../../review_20261004/pd_loaded_diagnostics/reports/loaded_elastic_a8_science.png) is still hash-matched to its frozen manifest entry (`b3df2be1cff05a4060d4b81112686287f9d925aa8d5511bb5af9ad9ac317764e`). Its rest/current overlay uses true 1× scale; the expected overlap follows from about 21 nm maximum displacement on a 5 mm body. The plotted axial field, displacement trace, and internal-force-density trace are labeled. The image states that it shows no crack, contact, or energy closure. Root reports that they have also viewed and accepted this representation as a narrow A observer figure.

The new B figure [native_b_affine_threshold_audit.png](native_b_affine_threshold_audit.png) was generated only from the existing B CSVs. Its left panel shows normalized stretch with the strict `εc` threshold and separates the 632 broken pairs from 2,430 retained pairs in the `above` branch. Its right panel pairs the low-ID→high-ID and reverse directed stretches; all points lie on `y=x`, with reported maximum difference zero. The title and footer explicitly say the map was prescribed to every point, with no time integration and no free fracture. Root reports independently viewing and accepting the same representation.

## Facts, interpretation, and unknowns

**Facts from receipts and raw data:** The 125-point native graph has 6,124 directed/3,062 canonical bonds. B has native failure transitions on 632 reciprocal pairs only in its `1.5 εc` branch. A's separate eight-step loaded observer run has nonzero motion/force response and no bond event. Both figures are sourced from actual output data. All 26 receipt-pinned sources still match.

**Interpretation:** B verifies that this pinned OSB/CriticalStretch failure path reacts to a prescribed geometric threshold crossing, with state reciprocity and point-to-bond length mapping intact. It does not verify dynamic crack initiation or propagation. A verifies a narrow, very short loaded-response observer path, not fracture or quasistatic material response. Neither case evaluates rotation objectivity or independent remeshing/orientation sensitivity.

**Geometry and model limits:** Here `δ/h=3`, and the cube is only `5 mm = 1.67δ` across. The most interior point is 2.5 mm from its nearest exterior face, less than `δ=3 mm`; no point has a complete spherical horizon inside the specimen. The driver uses `OSBMaterial()` and the pinned constructor defaults to `NoCorrection` ([source line 94](../../../src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/physics/ordinary_state_based.jl:94)). The geometry is therefore strongly surface-dominated. This does not invalidate B's branch diagnostic, but it sharply limits transfer to a coupon or artifact.

Energy ledger is `UNKNOWN_NOT_A_GATE`; contact is `NOT_RUN`; physical material qualification is none; free fracture, fragment generation, and artifact-scale behavior are not established. The native OSB strain-energy routine recalculates weighted volume and dilatation using the active graph ([ordinary_state_based.jl](../../../src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/physics/ordinary_state_based.jl:207)), so an instantaneous strain-energy drop alone is not a verified fracture-dissipation ledger.

## Next source-only gate proposal

The next source review and benchmark route is recorded in [next_free_evolution_coupon_source_plan.md](../../../experiments/native_resume_20261005/next_free_evolution_coupon_source_plan.md). It separates (1) objectivity under rigid translation/rotation of one fixed point set from (2) mesh orientation and `h`/`δ` convergence; chooses a coupon with meaningful interior support and a tested surface correction; distinguishes all-point affine prescription from grip-only time-dependent velocity conditions integrated by VelocityVerlet; adds reaction-force/work and a closed energy ledger before a free-path fracture run; and proposes a geometric notch with no manually deleted bonds so the continuation path remains unconstrained. This is source-only design, not authorization to edit or run the next case.

One primary quasi-static glass-fiber coupon study reports horizon/spacing sensitivity and favors three particle spacings in that studied setup; treat it as context, not a universal setting. A separate 2D/3D dynamic-fracture study explicitly examines horizon, specimen size/thickness, surface effects, crack observables, and elastic/kinetic/dissipated energy. [Mesh sensitivity study](https://doi.org/10.1016/j.proeng.2017.02.116), [dynamic-fracture study](https://doi.org/10.1007/s00466-021-02017-1). The official API documents `OSBMaterial()` as using `NoCorrection` by default, exposes `EnergySurfaceCorrection`, describes time-dependent `velocity_bc!`, and limits `displacement_bc!` to NewtonKrylov; VelocityVerlet is the time-integrating solver path. [Peridynamics.jl public API](https://kaipartmann.github.io/Peridynamics.jl/stable/public_api_reference/).

All older design notes, receipts, and raw artifacts remain unmodified. The old B-not-run manifest status is retained as historical design metadata; consult this post-run report and its hash-indexed audit for B's actual observed result. The new plan does not repeat an old negative trial and does not infer a physical pass from the B threshold event.
