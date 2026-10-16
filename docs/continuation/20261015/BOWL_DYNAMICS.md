# Native spherical-bowl dynamics adapter

This continuation adds a bounded, source-pinned native dynamics path for the two
VesselForge `forge.asset.v1` spherical-bowl assets. The producer is
`src/pgfrag/bowl_dynamics.py`; its Julia runner is
`src/pgfrag/resources/julia/run_bowl_dynamics.jl`. Preparation creates a fresh
attempt and immutable source snapshot. It does not launch Julia unless the caller
explicitly selects `--execute`.

## Frozen case and limits

The first comparison uses concentric lower hemispheres with `Ro=12 mm`,
`Ri=8 mm`, and center `(0, 0, 12 mm)`, sampled at cell centroids with `h=1 mm`
and horizon `delta=2.5 mm`. A separate `h=0.8 mm` case is a spatial-sensitivity
probe at the same horizon; it is not a convergence claim. Selected cells retain
the full `h^3` volume, so reported lattice mass is a discrete model mass and
the curved wall is not cut-cell corrected.

Mechanics are copied from and bound to `manifest.material.mechanics`:
`rho=2300 kg/m^3`, `E=70 GPa`, `nu=0.25`, and `Gc=7.6 J/m^2`. Every case is
marked `UNVALIDATED_PROTOTYPE`; these are illustrative values without specimen
measurement or calibration. Downstream case or runner overrides are rejected.

The native model is `BBMaterial{NoCorrection}` with `CriticalStretch`, advanced
by Peridynamics.jl 0.5.4 `VelocityVerlet`, with one Julia thread and free
release. The frozen loads are `alpha=0, 100, 300 s^-1`; for nonzero alpha the
initial affine rate tensor is `diag(alpha, -0.25 alpha, 0)`. Velocity is
projected to zero total linear and angular momentum about the mass center. No
impact, drop, contact, external grip, precrack, or boundary condition is used.
The baseline is 400 steps at safety factor 0.5. The half-dt sensitivity uses
800 steps at safety factor 0.25, preserving the same modeled duration.

The pre-solver budgets are 12,000 material points, 2,000,000 canonical pairs,
1,000 steps, 120 seconds, and a 256 MiB per-attempt output cap. The cap was
raised from 128 MiB before the first native solve so the same-duration 800-step
half-dt comparison could preserve its complete audit history. This changed an
I/O allocation limit only; no physical, numerical, or data-qualification
threshold changed. Receipts report peak materialized and persistent output
separately. Attempts are never reused or overwritten.

## Cavity filter and native integration

`GeometryFilteredBody <: Peridynamics.AbstractBody` forwards the wrapped native
`Body` properties. Its type-specific `find_bonds` first calls the native bond
builder, then removes any pair whose closed reference segment passes through the
open inner sphere. It updates `n_neighbors` before Peridynamics derives halo
points, bond ranges, kernels, or storage. Each excluded canonical pair is
recorded as an `INITIAL_GEOMETRIC_EXCLUSION`; it is not a fracture event. The
adapter's `pre_submission_check` delegates to the native body check. There is no
global override of Peridynamics methods.

The runner uses the native `VelocityVerlet` time-stepper operations and native
bond force/damage/energy routines. It checks the no-safety CFL limit and the
accepted safety-scaled time step before stepping. At the initial state it sets
acceleration from the initialized native internal/external force density while
leaving the imposed IC velocities unchanged, so the first native half-step has
the corresponding initial acceleration.

For the native BB material,

```text
K       = E / (3 (1 - 2 nu))
c       = 18 K / (pi delta^4)
epsilon_c = sqrt(5 Gc / (9 K delta))
U_ij    = 0.5 c V_i V_j L_ij epsilon_ij^2
```

CriticalStretch failure is the native strict test `stretch > epsilon_c`.
Events are recorded after the current position update, with stretch and pair
potential evaluated at the same state immediately before the pair's removal.
The cumulative break ledger sums these discrete `U_ij` values; it is not
`Gc * crack_area`.

## Result and provenance

`result.json` uses schema `pgfrag.native_bowl_result.v1`, reports `status=PASS`
only after native completion and packaging, binds the raw Forge manifest bytes
through top-level `source_forge_manifest_sha256`, and provides the compressed
NPZ descriptor at top-level `arrays`. The array contract stores the sorted
canonical zero-based pair roster (`i < j`), little-bit-order packed active
states after each state's native failure evaluation, every point state and
damage state, component representatives (minimum point ID, including
isolates), first-break rows, native strain-energy density, and a native
per-state `K/U/D/P/L/force` ledger. Checkpoint acceleration and force are
included as separate arrays.

`argv.json` freezes the prepared command and source receipts. The exact source
manifest is copied as `source_snapshot/forge_manifest.json`; the runner, Python
producer, Julia project/manifest, backend provenance and vendored
Peridynamics.jl bytes are also snapshotted and hashed. An execution uses the
snapshot runner, the isolated `third_party/pgfrag-runtime/depot`, one Julia
thread, and a separate `run_receipt.json`; canonical and snapshot pins are
checked before and after native execution. Raw full-trajectory state files are
packaged into `arrays.npz` in the new attempt and removed after packaging;
native CSV ledgers remain available for independent readback.

The result is a bounded dynamics diagnostic on an approximate centroid
discretization. It does not establish material calibration, spatial or temporal
convergence, experimental fracture prediction, contact response, or validated
fragment mechanics.

## Executed native results and independent readback

The following source-pinned attempts completed from their frozen snapshots under
`pgfrag/data/continuation_20261015/`. The independent audits are saved under
`pgfrag/experiments/bowl_independent_audit_20261015/`.

| Case | Steps / `dt` | First-break events | Independent final active-graph components | Final discrete break energy `D` | Independent energy result | Wall / peak materialized / persistent output |
|---|---:|---:|---:|---:|---|---:|
| `zero/run_attempt_03` | 400 / `6.4072827964e-8 s` | 0 | 1 | `0 J` | `PASS`; zero-load state stayed at rest and total-energy drift was 0 | `16.109 s` / `69,587,124 B` / `4,786,325 B` |
| `alpha100/run_attempt_01` | 400 / `6.4072827964e-8 s` | 314 | 1 | `7.1189775856e-5 J` | `PASS`; maximum absolute total-energy drift was `0.2429%` under the auditor's peak-total-energy normalization | `18.328 s` / `114,128,038 B` / `49,327,239 B` |
| `alpha300/run_attempt_01` | 400 / `6.4072827964e-8 s` | 21,220 | not independently established | `5.6602573956e-3 J` | Audit attempt is `INPUT_ERROR`: the plotting dependency could not be imported to render the requested audit figures. This is not a physics failure or a passed independent audit. | `18.266 s` / `116,184,458 B` / `51,383,659 B` |
| `alpha300_halfdt/run_attempt_01` | 800 / `3.2036413982e-8 s` | 20,768 | 58 | `5.3156816325e-3 J` | `PASS`; final total mechanical energy is `2.1761%` above the initial energy. The audit's maximum absolute drift is `2.1298%` because its metric divides by the maximum absolute total energy over the trajectory. | `25.828 s` / `224,748,773 B` / `96,031,574 B` |

Each completed case used 2,516 points, 70,648 retained canonical pairs, and 24
initial geometric exclusions. The half-dt attempt preserves the same modeled
duration as the 400-step cases. Its 58 graph components include isolated points
and are connectivity labels only; they do not establish 58 rigid fragments,
positive-area crack faces, or R2. Its largest graph component still carries
`0.0053544 kg` of the `0.0057868 kg` discrete model mass. The recorded break
energy is the sum of pair potential removed at discrete break states, not
`Gc * crack_area`.

Independent audits `zero_03_root01`, `alpha100_root01`, and
`alpha300_halfdt_root01` completed with `PASS`. The alpha300 independent audit
directory contains the input error described above, so no independent audit
result is claimed for that baseline. The `h=0.8 mm` fixed-delta spatial
sensitivity case remains unrun in this snapshot; the available half-dt run is a
time-step sensitivity study and does not establish convergence.

The earlier failed zero-load attempts are preserved. `zero/run_attempt_01`
stopped because the isolated Julia depot lacked `StaticArrays`; `zero/run_attempt_02`
stopped because the adapter lacked a native `parameter_setup_type` dispatch for
`GeometryFilteredBody`. The successful `zero/run_attempt_03` uses a fresh
attempt and its own source receipts. None of these runs turns the uncalibrated
prototype values into measured ceramic properties.
