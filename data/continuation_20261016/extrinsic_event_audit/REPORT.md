# Extrinsic event offline audit

This is a read-only reconstruction of the frozen `coarse654_20261004_01` evidence. The script does not import Akantu, start a solver, edit historical evidence, or evolve a native state. The script and machine-readable output are in this directory's paired experiment/data folders.

## Signed work and first negative getter

Areas were recomputed from the paired Tri3 faces in the saved post-insertion mesh: 8 facets, 24 quadrature points, and 3.9999999999999998e-06 m² total. Each QP receives one third of its actual triangle area. The recorded 650D→651D vector path gives signed ΔW = 3.6823662360605853e-07 J, R = 1.956104102232648e-07 J, and D = 1.7262621338279373e-07 J; the native D getter is 1.7262621338279373e-07 J. The 651 reconstruction matches the native getter within the stated roundoff tolerance.

For 651D→652D, signed ΔW = -5.6097918000056439e-05 J, cumulative W = -5.5729681376450381e-05 J, R = 0 J, and D = -5.5729681376450381e-05 J; the native D getter is -5.5729681376450387e-05 J. 22 of 24 QPs have negative work increments at step 652. The first negative native getter is 652D. No energy values were clamped or replaced by magnitudes.

The saved maximum Euclidean opening increases from 2.17288819e-09 m at 651D to 5.34502098e-07 m at 652D, while 24 QP traction vectors at 652D are exactly zero. The signed saved vector path explains the first negative dissipated-energy value. It supports a leading unresolved near-zero-opening/directional-stiffness mechanism; the offline record does not establish it as the sole cause.

## Insertion state and operative mass

At step 650, the recorded B state has 45 nodes, while C and D have 54 nodes and 9 appended nodes. In C, the duplicated nodes have equal displacement and velocity to their source nodes, and their appended mass rows are zero. D is the post-solve state: its appended-node masses sum to 5.0000000000000004e-06 kg per component, total mass is 4.0000000000000003e-05 kg per component, and total mass agrees with B to the recorded numerical precision. The zero appended masses in C therefore are not used as the operative post-insertion mass.

P/L/K are recorded for B, C, and D in `audit.json`; B and C values use the mass arrays saved at those phases. C's appended zero masses make its pre-solve P/L/K a stale-mass diagnostic, while D's reassembled mass is the operative post-solve result.

Using D's saved post-solve mass, velocity, displacement and rest coordinates gives P = [4.963084e-24, 5.376674e-23, 2.869928e-06] kg·m/s, L about the mass center = [3.258181e-11, -3.258181e-11, 2.100263e-26] kg·m²/s, and K = 1.648483312076419e-07 J. Acceleration values are not in the NPZ snapshots; the tagged source says `onNodesAdded` copies acceleration to duplicates, but this evidence set cannot verify its numeric values. Duplicated displacement remains equal at the saved D sample; velocity can differ after solve.

## Force and constitutive-data limits

The snapshots save only the assembled `internal_force` vector from `model.getInternalForce()`. Tagged source ancestry shows the cohesive override computes cohesive tractions and then calls the base internal-force assembly, so the saved force is a combined vector, not separately attributable bulk and cohesive vectors. `external_force` and acceleration are not among the saved arrays; the trace has only a scalar top reaction estimate. Tagged `checkInsertion` computes per-QP `new_sigma` from the effective facet stress, then uses `delta_c_eff = sigma_limit / new_sigma * delta_c` when the nominal `delta_c` is explicit; it stores those values in internal fields. This input has an explicit nominal `delta_c`, but the phase capture omitted `sigma_c_eff`, `delta_c_eff`, `delta_max`, damage, normals, and stress tensors. Mesh geometry can recover geometric normals, while the recorded wheel is 5.0.7.post1 and the local inspected source is tagged v5.0.7. Accordingly this audit does not claim a wheel-specific 651 traction-law reproduction or a finite 651 Jacobian/eigenvalue result. It does not infer a parameter from 651 traction and then use that same traction as its validation target.

The cached source manifest's files and the installed wheel are separately hashed in `audit.json`. The wheel/source version distinction is retained; the local source cache does not contain `material_cohesive_linear_inline_impl.hh`.

## Reproduction

From the repository root in PowerShell, run:

```powershell
& 'C:\Users\xiemo\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' 'pgfrag/experiments/continuation_20261016/extrinsic_event_audit/audit_event.py'
```

The script only reads the frozen input/evidence files and writes `audit.json` and this report under `pgfrag/data/continuation_20261016/extrinsic_event_audit/`. Its own SHA-256 and all required input/source hashes are in the JSON.
