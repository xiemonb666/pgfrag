# Free-fracture pilot — execution review

Date: 2026-10-05 (UTC execution window: 09:49:26–09:50:06)

## Authority and source identity

This report records the single fixed execution explicitly authorized by root after review of source freeze `1f4c3e9b10a607d303bbc3435bf22ae8b591ce8f6b5651aacb6cba3c3bfd7217`. The approved launcher was invoked once. There was no retry, source change, threshold change, fracture post-edit, or second run. The freeze and all 84 pinned files were rechecked after execution: 84 pins, zero byte/hash mismatches. The run receipt is authoritative for process and budget status.

Correction to an earlier progress message: I incorrectly summarized the `SOURCE_REVIEW.md` hash as `d496…`. That was a message-level hash error. The review file on disk and its actual freeze pin both have SHA-256 `30f12cfd0c2235bd5c0bc34d8b79c8bea391022c71f11ff4f9fac764a088ce42`; no frozen file was changed to make this correction. The pilot source identities are driver SHA-256 `9922d52c13468df5152d88c1fdb1b7c9fc6ce473a1a3a9e56507aaf37479f50b`, fixed case SHA-256 `45b9f8bb5d10522fa2e6fa91ca8d944e2d6af27396034b489e27a026783bbedd`, and launcher SHA-256 `97434a0f5b17dd93ae7a0064f98515b4651aad583b0a039b7907c646c89ba69d`. The accepted intact baseline references are fixed case SHA-256 `e6aa57f3edf136970e9f7108a1df2535a1e7c0d4cb56ee89fee8d858cd3c603d` and independent review SHA-256 `96eb441f9afb7b07d5cc672c4be4681a7b9bf587f7a4ad291785e196e15178ec`.

## Execution result

`run_once.ps1` returned `COMPLETED_EXIT_0`; the pinned native Julia process exited 0 in 39,962 ms. It was not killed and did not time out. Disk, output, and attempt caps were clear; `final_budget_ok=true` and `result_authoritative=true`. The result directory is 88,921,168 bytes against a 199,229,440-byte cap; the final attempt directory is 88,924,215 bytes against a 209,715,200-byte cap. The frozen preflight worst-case output estimate was 163,937,072 bytes. The run covered 4,500 points, 96 steps at fixed `dt=3.75e-8 s`, and ended at `3.6e-6 s`. The actual native raw stability bound exactly matched the accepted intact baseline at `3.4242302121147243e-7 s`; the fixed time step was not autotuned.

Stdout reports `points=4500`, `interior=1134/4500`, `steps=96`, `canonical_bonds=204778`, and `events=6810`. Stderr contains only the pinned solver's warning that an explicitly specified step size should be checked against its CFD condition. This pilot's numerical path gate remains `INTACT_NUMERICAL_PATH_REVIEWED`; its physical loading validity remains `UNKNOWN_UNCALIBRATED`.

## Native damage and graph observations

The coupon is the complete analytic uniform box, with no notch, point-removal mask, precrack, or manually edited bond state. Every native canonical bond was initially active and permitted to fail. The pinned OSB `CriticalStretch` criterion uses native stretch exceeding

`epsilon_c = sqrt(5*Gc / (9*K*delta)) = 0.002081665999466133`.

Here the uncalibrated prototype values set `Gc=45 J/m²`, `K=3.205128205128205e9 Pa`, and `delta=1.8 mm`; the registered loading target was a fixed 1.5× threshold (`0.0031224989991991995` engineering strain). These are diagnostic inputs, not measured ceramic properties. In this implementation, `Gc` parameterizes the failure threshold; a dissipated-fracture-energy history was not computed, so this run does not establish that `Gc` was dissipated or validate a fracture-energy balance.

Native failure events were observed. The first was at step 4 (`1.5e-7 s`), pair 19–20, with stretch `0.003587303573499008` and stretch/threshold ratio `1.7232848950883632`. The run recorded 6,810 cumulative native pair transitions; the final broken canonical-pair count is also 6,810, leaving 197,968 of 204,778 canonical pairs active (395,936 active directed bonds). The driver's run-time checks for no pair reactivation, event/broken-pair agreement, `n_active = 2 × active canonical pairs`, and agreement between stored and summed per-point active counts completed without aborting.

A row-level crosswalk of all 6,810 events against point_reference.csv shows that every failed pair has a right-grip endpoint; there are zero events without one. The only event-pair reference-x combinations are 4.5↔5.7 mm and 5.1↔5.7 mm. This localizes all observed failure to bonds joining the right-grip layer to its neighboring x-layers; it is a grip-end separation, not internal crack growth.

At step 96, the active-bond graph has two connected components. The minimum-point-ID labels identify 4,275 points under label 1 (mass `0.00193914 kg`, about 95%) and 225 points under label 20 (mass `0.00010206 kg`, about 5%), out of total mass `0.0020412000000000026 kg`. The 225-point component is the right-grip layer. Thus the two-component result records right-grip-layer detachment; it is not evidence of an internal crack, crack-front propagation through the coupon, or physical shard fragmentation. Component count was not a target, and the native event/graph observation must not be presented as fracture-physics success. Because `NoCorrection` leaves free-surface support truncation uncorrected and only 1,134 points (25.2%) have complete spherical horizons, the result carries known boundary-support limitations and does not qualify free-surface response.

## Loading and energy limits

The case uses a fixed 1 μs ramp, 0.6 μs hold, 1 μs unload, and 1 μs dwell. The estimated longitudinal transit time is `7.163e-6 s`; the ramp is about 0.140 transit times. This is a short explicit dynamic diagnostic, not a quasi-static test. The intact baseline's independent review covers numerical path/state and arithmetic checks, not scientific load validity, calibration, or fracture behavior. No claim of valid physical loading is made here.

Energy output is raw solver diagnostic data. At the final step, native OSB strain energy is `1.9385e-5 J`, full-step kinetic energy is `2.01994e-5 J` (free-DOF portion `2.01994e-5 J`; constrained-DOF raw portion `1.67e-14 J`), accumulated discrete BC kinetic work is `0.00309608 J`, and the recorded diagnostic residual is `-0.00305650 J`. The constrained grip positions are held and constrained `velocity_half` is zeroed during dwell, while stored full-step grip velocity can remain nonzero because of VelocityVerlet time staggering. These conventions and the large residual preclude energy closure or material qualification; the actuator-work trapezoid is separately marked surrogate-only.

## Run artifacts

All hashes below are SHA-256 of actual post-run files. Paths are relative to the workspace.

| Artifact | Bytes | SHA-256 |
|---|---:|---|
| `pgfrag/experiments/free_fracture_pilot_20261005/run_attempt_01/launcher_receipt.json` | 2,755 | `0bae8302c765effdd2614f72a2e1802d4779aa9ea4910829a1281f063c4eb557` |
| `pgfrag/experiments/free_fracture_pilot_20261005/run_attempt_01/native_stdout.log` | 118 | `84478924ac9ba53b17a977e36b87a1bfa12ede462605963fd5641070e7bebab8` |
| `pgfrag/experiments/free_fracture_pilot_20261005/run_attempt_01/native_stderr.log` | 227 | `b87722a6adc7800d3eea59db8e3e633b3e2cb9e160eaf980b5772bce56980937` |
| `pgfrag/experiments/free_fracture_pilot_20261005/run_attempt_01/results/canonical_bonds.csv` | 12,255,264 | `cc0f33ecdc350fe28c6f3bcffe506bbfa327a4df2ad84494e5aa6c38f4d8aa47` |
| `pgfrag/experiments/free_fracture_pilot_20261005/run_attempt_01/results/energy_reaction_history.csv` | 65,563 | `8287c68ca0cc612b406d87759c4792e15eb8836722273a5e14df0aec0f07f667` |
| `pgfrag/experiments/free_fracture_pilot_20261005/run_attempt_01/results/graph_component_labels_layout.toml` | 684 | `7a6a0e137dcccd807f6942b807699f65e9a3f062d26c8e817296a90ebc195899` |
| `pgfrag/experiments/free_fracture_pilot_20261005/run_attempt_01/results/graph_component_labels.i32le` | 1,746,000 | `7575f682d5e5efbb96f3128102550b14f0a906b5f92cc23df1e4a3976cff4d16` |
| `pgfrag/experiments/free_fracture_pilot_20261005/run_attempt_01/results/graph_state_history.csv` | 17,567 | `37c12ca66913c1d192b68364adfe40ce70bae01779f2cc90096788673c645f62` |
| `pgfrag/experiments/free_fracture_pilot_20261005/run_attempt_01/results/native_failure_events.csv` | 957,282 | `e37e0785da631d62be4413d54577ba14439d24b347fa728e832124ea0a5e30b6` |
| `pgfrag/experiments/free_fracture_pilot_20261005/run_attempt_01/results/point_reference.csv` | 520,394 | `7194fadbcd2af23b3205aa4349b471d6cb9457d24ee84b88e91404af35579257` |
| `pgfrag/experiments/free_fracture_pilot_20261005/run_attempt_01/results/point_state_history.f64le` | 73,332,000 | `2512eb66bebf85fed1d8d830b504f994fd5a01d6beab4af5ff3872543ecd7edf` |
| `pgfrag/experiments/free_fracture_pilot_20261005/run_attempt_01/results/point_state_layout.toml` | 1,816 | `6e5f1ca4ccf7c8680bf693a6f96236dca1e52382bf58126c2cfeaa276074959d` |
| `pgfrag/experiments/free_fracture_pilot_20261005/run_attempt_01/results/summary.toml` | 24,598 | `1e3396d30aa124a3371b47bcd6cc40eea4118312a18d685bbfbe3fdaa60a9c8c` |

The result artifacts remain under `pgfrag/experiments/free_fracture_pilot_20261005/run_attempt_01/`. This report is a post-run review record and was not part of, nor added to, the source freeze.


