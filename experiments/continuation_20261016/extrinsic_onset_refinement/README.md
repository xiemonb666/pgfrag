# Extrinsic onset refinement protocol

This folder contains the root-gated native explicit driver and the offline
comparison plan for the frozen 96-T4 coupon. It does not contain a campaign
result. `--prepare` is a hash, mesh-area, and output-budget check; it does not
construct an Akantu model. A native run requires `--launch` and the exact
SHA-256 of the reviewed data-folder protocol, after root review.

The run matrix uses the same source MSH, `material.dat`, beta=3, eligibility
surface, and A/B/C/D order as the old probe. The common endpoint is
`666 * 5.8582308142940834e-8 s`, calculated with decimal arithmetic. Current
requested runs are `h`, `h/2`, and `h/4` (666, 1332, and 2664 steps). The
`h/8` run has 5328 steps and stays reserved for a separate root confirmation.
Each resolution is a fresh initialization; no historical NPZ is loaded as an
initial state.

The executable preserves one attempt per resolution. It refuses an existing
run directory, records a failure without retrying, caps each run at 300 seconds,
the campaign at 20 minutes, and all generated data at 256 MiB. The frozen
four-resolution forecast is in
[`frozen_protocol.json`](../../../../data/continuation_20261016/extrinsic_onset_refinement/frozen_protocol.json).

## Preparation commands

Run only these two commands to materialize the protocol and API inventory. They
do not initialize a model or call a solver step:

```powershell
wsl.exe --exec env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 /mnt/c/Users/xiemo/Desktop/ShardForge/third_party/akantu-runtime/venv/bin/python /mnt/c/Users/xiemo/Desktop/ShardForge/pgfrag/experiments/continuation_20261016/extrinsic_onset_refinement/refinement_driver.py --inspect-wheel-api
wsl.exe --exec env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 /mnt/c/Users/xiemo/Desktop/ShardForge/third_party/akantu-runtime/venv/bin/python /mnt/c/Users/xiemo/Desktop/ShardForge/pgfrag/experiments/continuation_20261016/extrinsic_onset_refinement/refinement_driver.py --prepare
```

The protocol builder prints the exact protocol SHA-256 and three launch
commands. Those launch commands are for root review only; they must not be run
until root explicitly authorizes the native campaign. Each command runs only
one resolution and refuses to overwrite or retry its result directory. The
reserved h/8 command is not emitted by the current launcher.

## Capture behavior

The driver repeats the old sequence exactly: extrinsic explicit initialization,
one initial lumped-mass assembly, registration of the physical `interface`,
the bottom zero-z condition, then for each step A before `applyBC`, B after the
top ramp BC and before `checkCohesiveStress`, C immediately after insertion,
`solveStep`, and D immediately afterward. It issues no extra explicit force or
mass assembly for observation. The cached v5.0.7 `getEnergy("kinetic")`
implementation calls lumped-mass assembly, and `onNodesAdded` marks mass for
reassembly. Thus the original C order—copy arrays, then call energy getters—may
contain a native mass-assembly barrier after the copied C arrays. The installed
5.0.7.post1 wheel has distinct provenance, so equivalence to the cached source
is unknown. On insertion steps, the driver saves separate C pre-energy and
post-energy getter copies and compares mass, forces, kinematics, blocked DOFs,
coordinates, and connectivity. D arrays are also copied before the legacy
energy getters; a separate D post-getter copy is recorded on insertion steps.
Those extra public getter reads request no second energy getter or explicit
assembly. Energy values stay labeled with their post-array getter phase.

Each A and B getter result is copied while the run is at that phase. The driver
keeps the current A and the previous four no-insertion B states in memory; on
an insertion it writes dense A/B/C/D snapshots and the preceding B history. It
also preserves the final four B states if no later insertion flushes them. D
captures are written in 128-step NPZ chunks through the common endpoint.

Every state read copies displacement, velocity, acceleration, mass,
internal/external force, blocked DOFs, bulk material stress, native cohesive
internal fields, native energy getters, rest coordinates and connectivity.
All D nodal arrays, bulk stress, and cohesive fields are saved. Topology and
rest-coordinate hashes are saved in the every-step scalar trace; full
rest-coordinate/connectivity arrays are persisted in the initial and dense
event snapshots. The generic `facet_ancestry` module checks every sampled
snapshot against the frozen source cell-row/local-corner lineage and exact rest
coordinates. This depends on unchanged native bulk-row ordering; errors remain
`UNKNOWN` and are never repaired by nearest-coordinate matching.

At B, pre-check candidate facet stress and geometric normals are copied from
the native `FacetsFEEngine` using the native check-facet mask. The archive
retains native slot indices and the within-snapshot QP grouping. No candidate
slot is joined to a cohesive FID, and no QP order is carried across insertion.
After insertion, the native `normal`, opening, traction, effective strengths,
insertion stress, damage/history, contact fields, and energy fields are copied
again from fresh native views. Cohesive row/FID mapping is snapshot-local.

The actual wheel exposes public FEEngine shape functions and integration-point
counts. The one permitted zero-step shape-only preflight verified the raw
`getShapes` layout for this coupon: 696 rows group into 232 triangle slots with
three QPs each, using the same symmetric three-point Tri3 shape rule. Exact
native facet connectivity and coordinates recover each selected triangle's
area as 0.5e-6 m², and native normals align with the geometric normals. The
driver therefore records `area/3` for each QP only for these eight planar
affine coupon facets, retaining raw shapes and native QP order. No public
per-QP weight/Jacobian getter is claimed, and this geometric measure is not
generalized to other meshes or curved facets. Stable facet identity and
within-snapshot cohesive ordering still must be verified before signed work is
reconstructed. The top internal reaction sign remains an estimate. The wheel
also has no public `getTime()`; event time is recorded as exact driver `step*dt`
and the scheduled BC value, not as a native clock reading.

## Offline comparison

Read only the generated run folders. For each resolution, validate source,
protocol, runtime, run status, step count, archive row counts, and the exact
mesh/field hashes before interpreting histories. Preserve failures and partial
captures; do not restart from any saved array.

Compare the first insertion time and stable source facet identity set. For
each inserted cohesive facet, calculate `q = sqrt(delta.T @ B @ delta)` only
when the native opening and actual native normal are populated, with the frozen
`B = n*n.T + 9*(I-n*n.T)` for beta=3. Report the first positive q and q/h².
Compare `d ~ B^-1*t_insertion` only if the captured `insertion_stress` semantics
are confirmed from the native law. Never infer sigma from traction and then use
that sigma to reproduce the same traction. Contact opening/traction remains a
separate branch.

Compare post-insertion QP histories only after mapping each snapshot’s
cohesive connectivity to stable FIDs and confirming the within-element
snapshot ordering. The pre-insertion candidate-facet slots and post-insertion
cohesive QPs are independent orderings. Reconstruct signed QP work only when
the coupon-specific `area/3` measure, state timing, field semantics, and those
maps are supported; never clamp or take absolute values of energy increments.
Keep canonical cohesive contact energy separate from its `contact` alias.
Report native energy getters alongside raw nodal forces, mass, and blocked
DOFs, and keep C pre/post-energy mass states distinct.

For the provisional 2% matched-time screen, compare fixed observables at the
common endpoint and matched physical times using an L-infinity absolute
difference divided by a frozen scale: force estimate 280 N, cohesive opening
`delta_c=2.1714285714e-7 m`, native energies `Gc*A=3.04e-5 J`, and total mass
`rho*V=4.0e-5 kg`. Label reaction-force histories as estimates. The provisional
energy-defect threshold is `0.01 * Gc*A = 3.04e-7 J`; it cannot be called a
ledger pass until the actuator force and prescribed-kinematic phase are
verified. The frozen-state tangent/eigenvalue calculation is a diagnostic and
does not certify an onset CFL step.

The local onset ledger is anchored at the last verified pre-event D state in
the linear ramp. It computes `K_phys = 0.5*sum_free(m*v^2) +
0.5*sum_blocked(m*(g')^2)` from captured nodal mass and velocity and reports
native kinetic energy separately; the unresolved startup impulse remains
explicit. Boundary reaction uses `R_b = m_b*g'' - Fint_b - Fext_b` only for a
verified force phase. During a verified constant-velocity ramp with zero
external force this reduces to `-Fint_b`; trapezoidal boundary work remains an
estimate until actuator force and the prescribed-kinematic phase are proven.
The global energy ledger stays `UNKNOWN` when those quantities are unavailable.

The `q_first` trend is also checked for the possible initial layer: if
`q_first ~ h^p`, then `h²*lambda_max` scales like `h^(2-p)`. A shrinking
near-zero-opening layer requires physical history evidence; under-resolution
is not reported as the sole cause of the historical negative signed work.
