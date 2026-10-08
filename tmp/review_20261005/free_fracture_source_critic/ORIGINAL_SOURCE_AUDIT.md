# Original free-fracture pilot source audit

Audit scope: the original source snapshot for `analytic_osb_free_fracture_pilot_h0600_delta1800`, using `fixed_case.toml` (the experiment directory does not contain `free_fracture_pilot.toml`), `run_free_fracture_pilot.jl`, `run_once.ps1`, `source_freeze.toml`, and only the Peridynamics files named by that freeze. This was a read-only source audit. No Julia, Python, CAD, test, or Git command was run.

## Snapshot integrity

At the audit hash pass, all 19 entries in the original `source_freeze.toml` matched both the frozen byte count and SHA-256. The freeze file's own current SHA-256 was `a63b5d1a2fb0494975c0f510f444ff9b5e8c6ed23a6ababb1808138937bc743e`; the freeze file itself is not one of its pins. The contents below describe that original snapshot; a later authored revision must be reviewed against its new freeze hash.

| Pinned path | Bytes | Current SHA-256 | Match |
|---|---:|---|---|
| `pgfrag/experiments/free_fracture_pilot_20261005/run_free_fracture_pilot.jl` | 71872 | `a82fbf0afb459ff992251670687839f9d93f1ff1fd0e34d3949cb36121da8f3d` | yes |
| `pgfrag/experiments/free_fracture_pilot_20261005/fixed_case.toml` | 2051 | `16645c14fdd36a0cdf5c3b03decc1e4e7286e09e7e5fedeba1efdfef2e9d7ff8` | yes |
| `pgfrag/experiments/free_fracture_pilot_20261005/run_once.ps1` | 15914 | `0543a95833af1d55218fe8fddb4ddac032b916bdc452b4778535069c14f4b7d6` | yes |
| `pgfrag/tmp/review_20261005/free_fracture_pilot/SOURCE_REVIEW.md` | 7438 | `160a477ca1b1c7cdc413ec9d96eddd3e8fcfdc11248b02b84e2ad1e87910073a` | yes |
| `third_party/pgfrag-runtime/extract/julia-1.10.12/bin/julia.exe` | 167776 | `29a5cb5bde6a3f4a151f6e714591ca9344638303d9ed2a893ab31ada3de5d103` | yes |
| `pgfrag/src/pgfrag/resources/julia/Project.toml` | 107 | `eebff26de1f99a56b29a43ad767d2877bba5b6a0e4a03de1f909eab00f59f4ae` | yes |
| `pgfrag/src/pgfrag/resources/julia/Manifest.toml` | 22753 | `eb7d8902b84cfc64faaa0d45b0e6a3ae8534a252d06399f3db074da5a738e594` | yes |
| `pgfrag/src/pgfrag/resources/julia/backend_provenance.json` | 8223 | `7e8a4638224ca61d9fca5153a60a47c217989b6070e414229ec8e9bc0afe942d` | yes |
| `pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/physics/fracture.jl` | 6507 | `d8b0381584ca33a5ef33071d2e31abb66b96ed3c639facda2cd3bec80ca67075` | yes |
| `pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/physics/ordinary_state_based.jl` | 11286 | `65ffd614270c8f409a28f83e5def1ea1abc5f94ddcb29eac7b2cd7320e720846` | yes |
| `pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/discretization/bond_system.jl` | 14917 | `c98bed960c6b6789f015f02e400d6b3a99e09df22c529f386c384a1564bb69aa` | yes |
| `pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/physics/boundary_conditions.jl` | 26235 | `eca99fbf67cf5169a0ea8e34f72b949f3b607155b4a5837cfea1b2a000d2649b` | yes |
| `pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/time_solvers/velocity_verlet.jl` | 12397 | `584f0afbb1d0cd9808bfeb56d128327d28f7eca486666ffe2a7bbc6feeea48d6` | yes |
| `pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/core/threads_body_data_handler.jl` | 8327 | `38d01769f5d75eb594a984a74d1feea1ecd553f762876e25e1e4c2ce3d19f116` | yes |
| `pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/core/data_handler.jl` | 426 | `de1dba3e5d4425f181fb478692ac654077f1c55674236aedb45d963a59a93388` | yes |
| `pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/core/submit.jl` | 3674 | `c23e6ae68fa8331023325723457a88945242d49946d39e9ad04508e0dff68ffa` | yes |
| `pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/discretization/body.jl` | 18152 | `60762a4f52c9b6c0b8f23d295b266c15bc697b1e2a356ddd759a561236bd3421` | yes |
| `pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/physics/material_parameters.jl` | 8596 | `4a032a539214e3d1f1291d6e63b19ff6254533f849950c10d5a8e4b98fde96ab` | yes |
| `pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/core/parameter_handler.jl` | 2788 | `2630c4ba6651152071d6a4b7528a88e651e02a9c37337f39bde26cc0dd30264c` | yes |

The launcher checks the case, driver, launcher, and Julia executable against pins (`run_once.ps1:107-119`); the Julia driver checks every entry in the freeze (`run_free_fracture_pilot.jl:17-39`). Since `source_freeze.toml` is not externally pinned, its digest is recorded but not compared with an independent expected digest. A coordinated edit to the freeze and its pinned files could therefore create a new self-consistent set. This is a trust-anchor limitation, not a mismatch in this snapshot.

## Findings

### Run gate: original driver accepts the still-pending coupon state

The fixed case sets `load_validity_status = "PENDING_INTACT_COUPON_ACCEPTANCE"` (`fixed_case.toml:32-37`). The driver requires that exact pending string (`run_free_fracture_pilot.jl:833-839`) and then continues toward `submit`; it does not require a separate affirmative numerical acceptance result. If accepted-coupon evidence is a prerequisite to launching the fracture pilot, this is a **run blocker in the original snapshot**. The source also correctly keeps the result's load-validity claim pending (`run_free_fracture_pilot.jl:1128-1129`).

### Energy interpretation: held grip positions do not imply zero stored full-step velocity

Native `velocity_bc!` attaches the condition to `:velocity_half` (`boundary_conditions.jl:422-432`). Velocity Verlet then computes full velocity as `velocity_half + acceleration*dt/2` after the force update (`velocity_verlet.jl:257-294`). During a fixed-grip dwell the prescribed half-step velocity keeps positions fixed, but reaction acceleration can leave a nonzero stored full-step grip velocity. The driver validates fixed grip displacement (`run_free_fracture_pilot.jl:999-1024`) and sums kinetic energy over every point's full-step velocity (`:638-652`), while the summary says both grips remain at zero velocity in dwell (`:1120-1122`). This does **not** block the native position-controlled fracture evolution, but it blocks a literal zero-full-velocity or physical total-kinetic-energy interpretation. Keep that energy ledger diagnostic, qualify the grip statement, and/or report constrained-DOF kinetic energy separately.

### Event and output budget: bounded by this native graph, but no explicit event ceiling

The recorder stores every transition in a growable `FractureEvent[]` (`run_free_fracture_pilot.jl:73-97, 303-312`). The pinned native failure code only changes active bonds to inactive when `stretch > εc && fail_permit` (`bond_system.jl:288-305`; OSB override `ordinary_state_based.jl:137-156`), so the original implementation should have at most one event per canonical pair. However, the driver does not assert no reactivation or reconcile event counts against active-pair deltas. It reserves a fixed 48 MiB for bonds/events/CSV (`:963-966`) and writes the complete canonical-bond and event CSVs before checking the output budget (`:1046-1057`). At the launcher's maximum 180 steps, point-state bytes are 136,836,000 and graph labels 3,258,000; adding 48 MiB gives 190,425,648 bytes against the 199,229,440-byte output cap, leaving 8,803,792 bytes beyond that reserve (`run_once.ps1:94-101`). The reserve is plausible for the fixed geometry, and the process-level 200 MiB cap exists (`run_once.ps1:182-199`), but it is not an explicit event-size proof; a single full event CSV is written before either budget check can stop it. Add an event count/byte ceiling and preserve the native monotonicity/count checks. The planned authored revision reportedly adds these checks.

The recorder itself also allocates full point histories before solving (`run_free_fracture_pilot.jl:99-109`): six 4500×3×(steps+1) Float64 tensors plus scalar histories and graph labels. The case bounds point count and steps, but neither launcher nor driver sets a RAM limit. This is a fixed-size resource risk rather than an unbounded state-history allocation; the event vector is the data-dependent part.

### Completion status can disagree with final budget result

The driver writes `summary.toml` with status `ENGINEERING_DIAGNOSTIC_COMPLETED` and only then performs the final output-budget check (`run_free_fracture_pilot.jl:1191-1193`). If the summary itself pushes output over the cap, Julia exits nonzero but the already-written summary still says completed. The launcher correctly recognizes a nonzero native exit (`run_once.ps1:221-235`), but downstream readers must prefer the receipt/exit code over the summary in that case. The launcher also sets `diskCapExceeded` after observing an over-cap tree, then only calls `Kill()` if the child is still alive (`run_once.ps1:192-199`); if the child has just exited, the receipt can say `DISK_CAP_KILLED` although no kill occurred. These are receipt semantics issues, not solver blockers.

## Native mechanics checked

- Initial state and output chronology are coherent: thread submission builds the handler and initializes native forces before `solve!` (`submit.jl:99-112`; `data_handler.jl:1-4`); `solve!` emits the reference state, then Velocity Verlet applies half-step BCs, drifts, computes forces/failure, updates acceleration/full velocity, and exports each completed state (`velocity_verlet.jl:165-195`). The pilot's exact options overload observes those hooks (`run_free_fracture_pilot.jl:333-345`).
- The driver adds no external/artificial force and manually edits no bonds. `Body` initializes failure permission true (`body.jl:101-123`); `material!` sets permissions based on nonzero fracture parameters (`material_parameters.jl:64-85`; `fracture.jl:147-175`); the pilot checks permission and reference activity (`run_free_fracture_pilot.jl:898-908, 267-276`). Its only loading is native velocity boundary conditions on the two x-end grip sets. It does not use precracks, a point-removal mask, affine displacement replacement, or contact (`fixed_case.toml:6-18, 32-42`; driver `:898-908, 973-976`).
- The native critical-stretch rule is strict and the driver checks a transition's measured stretch against it (`bond_system.jl:288-305`; driver `:303-312`). Event time is the completed state time `step*dt` from the hook (`:222-228, 303-312, 340-345`).
- The point-ID mapping is explicitly one-to-one: one chunk, IDs exactly 1..n, and no halo points (`run_free_fracture_pilot.jl:221-239`). Per-point output, reference coordinates, canonical bonds, and graph labels are consistently indexed by those IDs (`:242-252, 123-164, 176-212`).
- The native OSB exporter recomputes bond lengths and energy density using the currently active graph (`ordinary_state_based.jl:238-286`); the driver zeros the storage first so unsupported points cannot retain a prior state's value (`run_free_fracture_pilot.jl:254-265`). The threshold formula matches pinned native code `εc = sqrt(5Gc/(9Kδ))` (`fracture.jl:110-129`; driver `:918-927`). The case correctly labels material values uncalibrated, `NoCorrection` surface truncation as open, graph components as connectivity diagnostics, and fracture dissipation as not computed (`fixed_case.toml:20-30, 59-64`; driver `:1094-1103, 1157-1169`).
- The secant function is consistent with native half-step velocity BCs: it sets `(u(t_n)-u(t_{n-1}))/dt` (`run_free_fracture_pilot.jl:348-379`), and the VV drift multiplies that half-step velocity by `dt` (`velocity_verlet.jl:181-195, 247-254`). This reproduces the prescribed grip displacement at sampled solver times; the code explicitly allows the sampled maximum to miss the continuous hold peak (`driver:1016-1024, 1135`). The trajectory remains a short explicit dynamic diagnostic, not quasistatic.
- `NoCorrection` with only 25.2% complete-horizon interior is a declared truncation limitation, not hidden geometry mutation (`fixed_case.toml:6-18, 20-30`; driver `:856-880, 1089-1096`). No physical-fragment, real-object, or calibrated-material claim is made.

## Acceptance view for the original snapshot

The pinned native event lifecycle, strict threshold, state chronology, no-manual-bond-edit policy, point mapping, and prescribed position path are internally coherent. The original case's pending coupon gate is the only clear prelaunch blocker if the root's stated acceptance prerequisite applies. The energy result must remain a discrete diagnostic because of the constrained full-step velocity issue and the explicitly unqualified residual. Output/event ceilings and completion receipt semantics deserve the planned small hardening before launch; none indicates a deterministic failure of the native solver in this fixed 4500-point case.
