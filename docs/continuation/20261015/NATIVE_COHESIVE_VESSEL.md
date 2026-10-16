# Native cohesive vessel experiment — 2026-10-15

**Status:** bounded native engineering diagnostic complete; the overall goal is
open. Both principal runs independently pass arithmetic/source auditing, while
physical qualification remains **NOT_PHYSICAL_QUALIFIED**. The count of
qualified physical packages remains **0**.

## Frozen case and native runs

The experiment uses the conforming tetrahedral bowl partition recorded in
[`mesh_result.json`](../../../experiments/cohesive_vessel_20261015/weak_seam_v4/mesh_result.json):
1,212 source nodes, 4,342 linear T4 elements, and 158 conformal x=0 interface
triangles. The baseline and half-step source cases are frozen under
[`weak_seam_v4`](../../../experiments/cohesive_vessel_20261015/weak_seam_v4/case.json)
and
[`weak_seam_halfdt_v2`](../../../experiments/cohesive_vessel_20261015/weak_seam_halfdt_v2/case.json).
Native cohesive insertion produces 1,317 nodes and 474 quadrature points. The
independent audits confirm that native coordinates, T4 ancestry, cohesive
pairing, nodal mass, and frozen source inputs match the declared sources.

The bulk material is `rho=2300 kg/m^3`, `E=70 GPa`, `nu=0.25`; it is an
uncalibrated prototype. The prescribed seam uses `sigma_c=10 MPa`,
`Gc=7.6 J/m^2`, `delta0=10 nm`, `beta=kappa=1`, and zero compression penalty.
All native DOFs are free. Opposed exterior patches receive the frozen
`1600*sin^2(pi*t/25 us) N` load through 25 us, followed by free evolution to
50 us. This is a prescribed weak-plane experiment, without contact or a
homogeneous free crack path.

| Frozen run | Steps / dt | Native nodal states | Independent audit |
| --- | --- | --- | --- |
| [`weak_seam_v4`](../../../data/continuation_20261015/cohesive_weak_seam_v4/result.json) | 3,244 / `1.5413070283600494e-8 s` | 3,245 | [`AUDIT_ARITHMETIC_PASS`](../../../experiments/cohesive_vessel_20261015/outaudit/weak_seam_v4_run04/audit.json) |
| [`weak_seam_halfdt_v2`](../../../data/continuation_20261015/cohesive_weak_seam_halfdt_v2/result.json) | 6,488 / `7.706535141800247e-9 s` | 6,489 | [`AUDIT_ARITHMETIC_PASS`](../../../experiments/cohesive_vessel_20261015/outaudit/weak_seam_halfdt_v2_run01/audit.json) |

Each audit covers every solver step and every native node, verifies exact
source hashes, and reports `qualification_status=NOT_PHYSICAL_QUALIFIED`.
Both runs reach all-QP full damage. The first all-QP failure occurs at
`11.12824 us` in the baseline and `11.12053 us` at half timestep. These are
native constitutive events for this frozen prescribed interface; they do not
establish a calibrated physical crack path.

For the measured discrete interface area `A_h=1.2562331510181291e-4 m^2`, the
reference `Gc*A_h` is `9.547371947737782e-4 J`. The baseline final native
dissipation is `9.545672102122966e-4 J`, `-0.017804%` relative to that scale;
the half-step result is `9.547046732062301e-4 J`, `-0.0034063%`. Maximum
algorithmic external-energy balance residuals are `2.19016e-7 J`
(`2.2940e-4 * Gc*A_h`) and `5.87063e-8 J` (`6.1489e-5 * Gc*A_h`). The
[`two-level fixed-mesh sensitivity record`](../../../experiments/cohesive_vessel_20261015/time_sensitivity_v4_halfdt_v2.json)
reports a final mass-weighted
displacement difference of `1.46764e-9 m`, maximum nodal displacement
difference of `1.08870e-8 m`, relative mass-weighted velocity difference of
`0.14754%`, and dissipation difference of `1.43980e-4 * Gc*A_h`. Two timestep
levels do not establish asymptotic convergence; the mesh, orientation, and
material remain uncalibrated.

## Closed native components and v5 packet

[`native_cohesive_fragments_v5`](../../../data/continuation_20261015/native_cohesive_fragments_v5/manifest.json)
is derived from the baseline run. Its separation record places the first
all-QP full-damage and zero-traction state at step 722. All saved subsequent
states satisfy the joint failure gate and the global x-plane separation test;
the minimum saved global x gap after debonding is `1.22001e-6 m`, and the final
gap is `3.84005e-4 m`. There is no continuous-time self-intersection
certificate and no contact dynamics.

The packet covers all 4,342 T4s and all 158 failed source facets. Its component
audit reports full rest/current volume and mass coverage, closed native
boundaries, verified source-exterior ancestry, and no shared native nodes:

| Component | T4s | Nodes | Boundary faces | Rest volume (m^3) | Current volume (m^3) | Rest mass (kg) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 2,239 | 670 | 1,116 | `1.26956100956746e-6` | `1.26957764400451e-6` | `0.00291999032200516` |
| 1 | 2,103 | 647 | 1,092 | `1.27032118833526e-6` | `1.27033638131261e-6` | `0.00292173873317110` |

The v5 `state.npz` adds native cohesive connectivity and the QP damage,
traction, and opening arrays; the earlier
[`native_cohesive_fragments_v4`](../../../data/continuation_20261015/native_cohesive_fragments_v4/manifest.json)
packet is retained. The v5 packet marks `physical_pass=false`,
`material_calibration_status=UNCALIBRATED`, `mesh_sensitivity=NOT_ESTABLISHED`,
`typed_R2_qualification=NOT_ESTABLISHED`, and `pretrain_eligible=false`. Its
closed-component result is topological eligibility for this prescribed seam,
not physical fragment qualification.

The existing [two-panel PNG preview](../../../experiments/cohesive_vessel_20261015/visualization_v4/native_fragments_v4_preview_01/native_t50us_weak_seam_two_panel.png)
is an older v4 visualization from the penalty-zero model. It is not a view of
the pending compression/contact implementation and should not be presented as
one.

## Bonded controls and failed acceptance gate

The independently executed [perfect-bond bulk control](../../../data/control_bonded_bulk/run_attempt_01/diagnostic.json)
uses the same source mesh, bulk material, 100 N patch load, pulse duration, and
50 us physical window as the nonsoft cohesive comparison. It has 1,212 nodes,
4,342 T4s, no cohesive elements, exact source node/T4 ordering, a per-node mass
audit over all three components, and the same patch IDs/weights as the cohesive
control (maximum weight difference `1.46e-16`). Its peak-load response is
`6.442884791215271e-7 m`, or `6.442884791215271e-9 m/N`, at 12.5 us.

[`control_comparison_v1.json`](../../../experiments/cohesive_vessel_20261015/control_comparison_v1.json)
compares that bulk control with `bonded_cohesive_control_v1`. The cohesive
control has no QP damage at 100 N and agrees with the bulk response within
`0.01314%` at the peak-load instant, but the declared 2% transient-compliance
gate **fails**: pulse-window RMS relative difference is `2.9908%`, and
full-window RMS relative difference is `308.34%`. The peak-only agreement does
not override the failed time-history gate.

The separate [`current_gap_diagnosis_v1.json`](../../../experiments/cohesive_vessel_20261015/current_gap_diagnosis_v1.json)
finds a minimum sampled current-normal gap of `-5.20205 um` in the nonsoft
cohesive control and `-0.420158 um` in `weak_seam_v4`. These are compressive
overlaps in a model whose old material card has zero penalty and no compression
contact branch. Neither result is a contact-qualified bonded control.

The other controls remain distinct and are not interchangeable: the zero-load
control has zero displacement and damage; the 100 N weak-seam low control has
maximum damage `0.00389715` and is therefore not an undamaged elastic control;
the 1 N low-elastic control has zero damage at a different load; and the
100 N high-strength cohesive control has zero damage but still exhibits the
negative sampled gap and fails the transient RMS gate. The new
`Kn=1e15` compression/contact branch is still under implementation and has no
run result or receipt here. It is not counted as a pass, and source changes
alone cannot upgrade qualification.

## A12 research package and remaining work

[`native_cohesive_v5_dataset`](../../../../relicweave-synlab/data/continuation_20261015/native_cohesive_v5_dataset/manifest.json)
is present as a research-only package generated from the v5 native fragment
packet. Its manifest records `physical_control_status=FAILED_TRANSIENT_BONDED_COMPLIANCE`,
`physical_pass=false`, `pretrain_eligible=false`,
`deployment_eligible=false`, `production_a12_complete=false`, no generated
relation labels, an unassigned research split, and an ideal software RGB-D
proxy that is not rig-validated. The new A12 bridge
[`native_cohesive.py`](../../../../relicweave-synlab/synlab/integrations/native_cohesive.py)
and v5 package still await root review; this record does not claim production
A12 completion or training eligibility.

The bounded engineering path now reaches native separation, closed
tet-derived components, and software-observation research data for one
uncalibrated preinserted weak plane. The overall objective remains open because
the matched transient bonded-compliance gate fails and compression/contact
behavior is unresolved. Required work still includes a reviewed and run
compression/contact control with energy checks, mesh/orientation sensitivity,
homogeneous crack-path and material calibration studies, continuous-time
separation/contact validation, and A12 bridge review plus production admission.
No measured ceramic response, physical-fragment qualification, or qualified
physical package is claimed.

Older frozen cases and failures—including `weak_seam_v1` through
`weak_seam_v3`, `weak_seam_halfdt_v1`, and the
`outaudit/weak_seam_v3_input_error` record—remain preserved by case and were
not overwritten or promoted by this update.
