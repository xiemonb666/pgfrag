# Final free-fracture pilot source audit

Scope: final fixed case, Julia driver, one-shot PowerShell launcher, 84-entry `source_freeze.toml`, and the referenced pinned Peridynamics sources. Read-only source/hash review only; no Julia, tests, Python, CAD, or Git execution.

## Freeze verification

The current freeze SHA-256 is `1f4c3e9b10a607d303bbc3435bf22ae8b591ce8f6b5651aacb6cba3c3bfd7217`, exactly matching the requested digest. All **84/84** frozen paths exist and match both their declared byte count and SHA-256.

| File | Bytes | SHA-256 |
|---|---:|---|
| `run_free_fracture_pilot.jl` | 82117 | `9922d52c13468df5152d88c1fdb1b7c9fc6ce473a1a3a9e56507aaf37479f50b` |
| `fixed_case.toml` | 3276 | `45b9f8bb5d10522fa2e6fa91ca8d944e2d6af27396034b489e27a026783bbedd` |
| `run_once.ps1` | 17288 | `97434a0f5b17dd93ae7a0064f98515b4651aad583b0a039b7907c646c89ba69d` |
| Pinned `SOURCE_REVIEW.md` | 7908 | `30f12cfd0c2235bd5c0bc34d8b79c8bea391022c71f11ff4f9fac764a088ce42` |
| Pinned intact-coupon independent review | — | `96eb441f9afb7b07d5cc672c4be4681a7b9bf587f7a4ad291785e196e15178ec` |

The two baseline references in the fixed case are also pins, and the driver checks that their case-declared hashes agree with the source-freeze hashes (`run_free_fracture_pilot.jl:842-850`). The driver validates every freeze entry (`:17-39`); the launcher validates its driver, case, Julia executable, and itself (`run_once.ps1:107-119`). The freeze file is not self-pinned and the launcher records, rather than hard-codes, its expected digest. The provided root digest supplies that external approval check; this is not a mismatch in the reviewed snapshot.

## Requested corrections verified

- **Numerical and physical status are separate.** The fixed case requires `INTACT_NUMERICAL_PATH_REVIEWED` and `UNKNOWN_UNCALIBRATED` as separate fields (`fixed_case.toml:41-52`); the driver requires both and keeps `load_validity_claim=false` (`run_free_fracture_pilot.jl:887-897, 1241-1247`). It pins the intact case and independent arithmetic review, without presenting them as calibration or physical validity.
- **The time step is fixed against the baseline.** The driver computes the native raw bound with one chunk, compares it to the pinned baseline, requires the exact predeclared fraction and 37.5 ns result, then rejects any phase-resolution need for clipping (`:968-1009, 1020-1027`). It verifies the constructed VelocityVerlet retained that step before submission (`:1075-1078`). No automatic retuning path remains.
- **Native event accounting is internally cross-checked.** The recorder rejects reactivation and checks each active-to-inactive transition against the strict permitted `stretch > εc` rule (`:307-319`). Every state checks reciprocal pair state, directed/canonical totals, per-point active-neighbor counts, total active counts, and cumulative events equal broken canonical pairs (`:327-349`). This agrees with the pinned native rule, which only changes active bonds to inactive when permitted stretch exceeds the critical value (`bond_system.jl:288-305`; OSB override `ordinary_state_based.jl:137-156`).
- **Event output is budgeted before solving.** The driver estimates point and graph binaries, bounds canonical and event CSV rows, reserves 16 MiB for other files, derives an event-count ceiling from the remaining result budget, and refuses the run if the envelope does not fit (`run_free_fracture_pilot.jl:1043-1063`). It checks the ceiling before each event append (`:314-315`). Writers still check actual output size after each file (`:1132-1159`), so a conservative-envelope miss fails the native process and is reflected in the receipt rather than being accepted.
- **Constrained kinetic energy is qualified.** The ledger now separates free and constrained-DOF contributions to raw full-step kinetic energy (`:654-685`) and explains that native velocity BCs set `velocity_half`, while the post-force full-step grip velocity may be nonzero (`:1234-1235, 1283-1289`). This is consistent with pinned native BC dispatch to `:velocity_half` (`boundary_conditions.jl:422-432`) and the Verlet full-velocity update (`velocity_verlet.jl:257-294`). The ledger remains explicitly diagnostic; `Gc` is only the threshold input and fracture dissipation remains uncomputed (`:1288-1290`).
- **The final launcher receipt governs acceptance.** The summary identifies itself as provisional (`run_free_fracture_pilot.jl:1160-1164`). The launcher records the actual exit code, whether it killed the child, output and attempt cap results, and makes `result_authoritative` true only after completed status, exit 0, no cap/timeout/kill, and the final byte checks pass (`run_once.ps1:206-240, 302-344, 349-388`). Cap status is no longer labeled as killed when the process already exited.

## Mechanics and launch assessment

The native lifecycle remains coherent: initialization computes reference forces before the reference-state hook; each completed Verlet step applies the half-step BC, updates positions, evaluates native failure/forces, updates acceleration and full velocity, then reaches the recorder (`data_handler.jl:1-4`; `submit.jl:99-112`; `velocity_verlet.jl:165-195`). The driver adds no artificial external force or manual bond edits, retains the full unnotched point set, and uses the native OSB energy exporter. Point IDs are checked as a one-to-one 1..n mapping with one chunk and no halo (`run_free_fracture_pilot.jl:232-241`). The native critical-stretch formula still matches `sqrt(5Gc/(9Kδ))` (`fracture.jl:110-129`; driver `:992-1001`).

**No mechanics or launch blocker found in the final source snapshot.** The known scientific limits remain accurately labeled: material uncalibrated, `NoCorrection` surface truncation open, dynamic short dwell, no physical-fragment claim, and no fracture-energy closure. These are interpretation limits, not source launch blockers for the bounded diagnostic. The final source is still marked not run and the launcher retains the parent-authorization comment; this audit did not execute it.
