# Free intrinsic cohesive patch source review

This source bundle implements the frozen 2 × 2 × 4 mm, 96-tetrahedron mesh as a bounded intrinsic-cohesive coupling observation. It prescribes only the top and bottom z components plus the three in-plane rigid-mode pins. All 54 cohesive-side components and the geometric interior bulk components remain free. The source freeze used for the native run is preserved in `experiments/free_cohesive_patch_20261015/frozen_run03/`; its completed run and independent readback are preserved in `experiments/free_cohesive_patch_20261015/coupled_readback_run02/`.

## Stability gate

The inspected `SolidMechanicsModel::getStableTimeStep()` loops over `_ek_regular` elements and does not include intrinsic cohesive stiffness. Its native value is therefore only a bulk sanity check. Before the saved run started, the frozen source received a parent-reviewed initial stiffness-to-lumped-mass bound over the free bulk and interface DOFs. `run_started.json` and `run_status.json` record `lambda_max = 6.148965883903991e16 s^-2`, the derived total stable step `8.0654565130e-9 s`, and solver `dt = 2e-9 s`; the runner passed this total-bound gate before step 1. The bound covers the initial high-frequency stiffness. It does not establish evolving-geometry, softening-accuracy, or contact stability.

The bilinear material input has no separate normal or shear stiffness parameter. The version-pinned law defines the elastic limit `delta_0`, sets `delta_c = 2 G_c / sigma_c`, rescales its internal peak value, and initializes `delta_max = delta_0`. For the frozen values, the pre-peak slopes are `K_n = sigma_c / delta_0 = 3.2236842105e15 Pa/m` and `K_t = (beta²/kappa) K_n = 2.9013157895e16 Pa/m` (`beta=3`, `kappa=1`). The configured compression penalty is zero, so this case has no compression penalty stiffness. The reviewed parent bound assembled these source-derived slopes directly; it did not rely on the native tangent API while its `sigma_c` versus `sigma_c_eff` behavior remained unresolved. These semantics are pinned to the [Akantu v5.0.7 bilinear law](https://gitlab.com/akantu/akantu/-/raw/v5.0.7/src/model/solid_mechanics/solid_mechanics_model_cohesive/materials/constitutive_laws/material_cohesive_bilinear.cc) and [linear-law traction/tangent implementation](https://gitlab.com/akantu/akantu/-/raw/v5.0.7/src/model/solid_mechanics/solid_mechanics_model_cohesive/materials/constitutive_laws/material_cohesive_linear_inline_impl.hh).

## Traction sign and native dissipation

The frozen case previously described the native traction as a resisting traction and proposed negating its opening work. The inspected code and prior pure-opening control do not support that convention. In the checked-in Akantu v5.0.7 source, `MaterialCohesive::updateEnergies` accumulates

```text
total_energy += 0.5 * (opening - previous_opening) · (traction + previous_traction)
reversible_energy = 0.5 * traction · opening
dissipated_energy = total_energy - reversible_energy
```

The native traction is therefore energy-conjugate to native opening for this API. Reversing the cohesive side orientation changes the paired opening and traction signs together and leaves their dot product unchanged. The older intrinsic, fully prescribed pure-opening trace shows positive z opening with positive z traction and a positive signed interface-work integral; its final dissipated energy reached `3.04e-5 J`, while the smallest raw value was only `-1.3877787807814456e-23 J`.

That sign result does **not** explain the negative dissipation reported by the separate extrinsic dynamic runs. Those traces first went negative two steps after insertion. The checked-in source confirms `SolidMechanicsModel::afterSolveStep` forwards to constitutive materials and the energy getters only read or copy stored energy fields, so querying `getEnergy` once per state cannot itself advance the history. However, the source bundle omits the full `ModelSolver` / `TimeStepSolver` and `savePreviousState` implementation needed to establish the precise snapshot order around `solveStep()`. It also does not establish byte-for-byte correspondence between the installed `5.0.7.post1` wheel and the reviewed `v5.0.7` tag. The earlier negative values therefore remain unresolved; no traction flip or absolute-value cleanup is allowed.

For the completed intrinsic run, the independent analyzer integrated raw, positive-area quadrature work from saved frames and compared it with native `dissipated + reversible`. The diagnostic gates passed: no sample was below the predeclared `-1e-12 J` roundoff floor and the history identity closed within `1e-12 J`. When a negative value was observed, the report distinguished a closed source-formula identity from a work/history mismatch. Each energy key was queried once per state; the inspected `getEnergy` and `getDissipatedEnergy` implementations read/copy stored arrays and do not advance the history themselves. The state-order question remains unresolved for the separate extrinsic episode because the reviewed source slice omits the solver lifecycle.

The top actuator-work estimate records `mass × acceleration`, internal force, and external force separately, then applies the residual convention declared in `case.json`. Its energy balance remains diagnostic because the constrained boundary response has inertia. It is independent of the cohesive signed-work identity.

## Files and checks

- `pgfrag/experiments/free_cohesive_patch_20261015/run_native_free_patch.py` — single-attempt native run source; saves the initialized state and every post-solve 24-point cohesive field to raw NPY arrays.
- `pgfrag/experiments/free_cohesive_patch_20261015/free_patch_contract.py` — MSH 2.2 parser, exact bulk-cell/local-ordinal ancestry, positive-area parent-facet pairing, and native signed-work helper.
- `pgfrag/experiments/free_cohesive_patch_20261015/analyze_native_free_patch.py` — independent raw-history, BC, topology, lineage, monotonicity, energy, and sign diagnostics.
- `pgfrag/experiments/free_cohesive_patch_20261015/test_source_contract.py` — native-free syntax/contract tests; does not construct an Akantu model or call the solver.
- `pgfrag/experiments/free_cohesive_patch_20261015/launch_free_patch.sh` — fixed WSL launcher using the existing isolated Akantu runtime.

## Completed native run and independent readback

The frozen source in `frozen_run03/` produced one native run, identified in
`coupled_readback_run02/run_started.json` and `run_status.json`. The run
completed 10,000 solver steps and saved all 10,001 expected frames. The saved
analysis reports `ANALYSIS_COMPLETE_DIAGNOSTIC; NO_PHYSICAL_OR_FRAGMENT_QUALIFICATION`
and all 15 declared gates as true. It verified unchanged 54-node / 96-bulk-tet
/ 8-cohesive-element topology, the free-DOF contract, exact source-cell and
local-ordinal ancestry, positive parent-facet areas, finite native histories,
the total stability bound, and signed interface-work closure.

All eight input interface facets have two source tetrahedron parents and
positive rest area; their summed area is `4e-6 m^2`. All 24 cohesive
quadrature histories were saved without cell averaging. The independent signed
interface work ends at `3.0410667858528943e-5 J`, compared with native
dissipated plus reversible energy of `3.0410667858528932e-5 J`; the largest
history residual is `1.0164e-20 J`. The reference `Gc * A` is `3.04e-5 J`, a
comparison within this toy law and geometry, not a material calibration. The
declared full-debond-plus-hold rule first registers at history index 4,130 and
remains satisfied through its hold window, whose end index is exclusive at
6,131.

The negative native dissipation remains at roundoff scale only in this
intrinsic run (`-7.2214e-15 J` minimum, with no sample below the declared
`-1e-12 J` tolerance), and the independent signed-work identity closes. This
does not resolve the separate earlier extrinsic dynamic run's negative-
dissipation episode; its solver history-order question remains open.

The lineages are based on source Gmsh face IDs, source tetrahedron element IDs and local vertex ordinals, native cohesive connectivity, and post-insertion native bulk adjacency. Coordinates are used as a rest-state consistency assertion, never to pair nodes. Each of the eight declared faces must have two input parents, exactly one native bulk parent per cohesive side, matching source node sets, positive equal side areas, and total area `4e-6 m²`.

## Scope

This is a free-DOF native coupling spike on one predeclared planar interface.
The run verifies lineage and positive area for that inserted interface; it does
not establish a propagating crack surface or fragment lineage. The configured
characteristic process-zone length is about `0.109 mm`, while the coarsest
relevant edge is `1 mm`; it cannot support a crack-path or process-zone
resolution claim. The toy material is not calibrated and must not be
transferred to the vessel. No physical fracture qualification, R2
qualification, mesh convergence, or pretraining eligibility follows from this
case.
