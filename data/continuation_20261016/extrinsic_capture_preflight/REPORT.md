# Native zero-step capture preflight

The actual Akantu 5.0.7.post1 runtime loaded the frozen coupon and exposed the native capture interfaces needed for a future, separately authorized replay. The extrinsic setup produced 96 tetrahedra, 45 nodes, eight eligible interface facets, and zero cohesive elements. It took about 0.54 seconds inside Python. No boundary condition, stress check, force assembly, explicit solve step, or model evolution occurred. The separate intrinsic capability run created eight cohesive elements during `initFull()` and then called `assembleMassLumped()` once as the initial mass assembly.

## Field and getter availability

`model.flattenInternalReal(name, kind, ghost_type)` is public in the wheel. Its binding documentation returns an element-type map, which the script reads with `model.flattenInternalReal(name, aka._ek_cohesive, aka._not_ghost)(aka._cohesive_3d_6, aka._not_ghost)`. The requested native field names all read successfully. In the extrinsic run, the cohesive arrays have zero rows because no cohesive elements exist. In the separate intrinsic capability initialization, they have 24 quadrature rows for eight Tri3 elements and three quadrature points per facet.

| Native field name | Extrinsic shape | Intrinsic capability shape | Initial intrinsic reading |
|---|---:|---:|---|
| `sigma_c_eff` | `(0, 1)` | `(24, 1)` | zero; no extrinsic `checkInsertion` occurred |
| `delta_c_eff` | `(0, 1)` | `(24, 1)` | the explicit nominal input `delta_c` default, `2.1714285714285715e-7 m`; not insertion-derived |
| `insertion_stress` | `(0, 3)` | `(24, 3)` | zero; no extrinsic insertion stress was captured |
| `normal` | `(0, 3)` | `(24, 3)` | zero cohesive field before traction evaluation |
| `delta max` | `(0, 1)` | `(24, 1)` | zero |
| `damage` | `(0, 1)` | `(24, 1)` | zero |
| `opening` | `(0, 3)` | `(24, 3)` | zero |
| `tractions` | `(0, 3)` | `(24, 3)` | zero |
| `contact_tractions` | `(0, 3)` | `(24, 3)` | zero |
| `contact_opening` | `(0, 3)` | `(24, 3)` | zero |
| `reversible_energy`, `total_energy` | `(0, 1)` | `(24, 1)` | zero |
| bulk `stress` | `(96, 9)` | `(96, 9)` | zero in this unloaded initialization |

The material object exposes public `getOpening` and `getTraction` getters, both with shape `(24, 3)` in the intrinsic capability run. Its direct Python getters `getSigmaCEff`, `getDeltaCEff`, `getInsertionStress`, `getNormal`, `getDeltaMax`, `getDamage`, `getContactTraction`, and `getContactOpening` are absent. The `flattenInternalReal` path is available for all of these fields. The bulk `MaterialLinearElastic3D.getStress(_tetrahedron_4, _not_ghost)` getter returns `(96, 9)`. The cohesive model's `getStressOnFacets(_triangle_3)` getter returns `(696, 18)` in the extrinsic initialization; both bulk and facet stress values are zero because no force assembly or solve was done.

The cohesive model's `getTime()` is absent from the actual wheel. The available time methods are `getStableTimeStep`, `getTimeStepSolver`, and `setTimeStep`; none returns the current model clock. Future event capture must retain the driver's step index, `dt`, prescribed time `step * dt`, and actual `applyBC` value. The frozen driver checks `checkCohesiveStress` after applying the next boundary value. The check interpolates existing bulk stress to candidate facets; the boundary setter changes displacement but does not recompute that stress.

The public DOF and state getters all work. In the extrinsic zero-step capture, acceleration, external force, internal force, mass, displacement, velocity, and blocked DOFs are copied with shape `(45, 3)`; the intrinsic capability initialization returns `(54, 3)`. Floating state and force arrays are zero; blocked DOFs are Boolean and all false. Mass is nonzero. Exact copied arrays are in the NPZ archives.

## Native normal and actual mesh geometry

The saved Gmsh file contains eight physical Tri3 interface faces. Areas recomputed from their actual node coordinates are `0.5e-6 m²` each and `4.0e-6 m²` total. The input triangle winding alternates between `+z` and `-z`, so the raw surface element winding is not a reliable proxy for the cohesive material orientation.

In the extrinsic setup, `FacetsFEEngine.getNormalsOnIntegrationPoints(_triangle_3, _not_ghost)` returns `(696, 3)`: 232 candidate facets with three quadrature points per facet. Selecting the eight `getCheckFacets()` entries in the same native triangle ElementTypeMap ordering yields `(8, 3, 3)`, with native normal `+z` at every selected quadrature point. This is a shape-aligned native candidate-facet check; no exact native facet-to-Gmsh-node-ID join is claimed. In the intrinsic capability run, the facet normal getter reports no remaining non-ghost triangle facet entry after intrinsic insertion. The native cohesive `normal` internal field is readable in both modes but stays empty/zero until a traction evaluation.

## Capture ordering for any later replay

The frozen probe orders its phases as A before `applyBC`, B immediately after `applyBC` and before `checkCohesiveStress`, C after the check/insertion and before `solveStep`, then D after `solveStep`. Capture the driver's prescribed time and applied BC at A/B; at B also copy displacement and the already-existing bulk and facet stress. At C capture new node IDs, public parent mapping if exposed, duplicated displacement/velocity/acceleration/blocked DOFs, and the pre-reassembly appended mass rows. At D repeat the arrays after `solveStep` and its operative mass update. Keep the phase label, step, `dt`, `step * dt`, BC value, and node/cohesive counts with each snapshot. Do not assemble internal forces solely to observe a state.

The intrinsic run only proves that getters and array shapes exist after a separate 8-element intrinsic initialization. Its zero `sigma_c_eff`, `insertion_stress`, and `normal`, and default nominal `delta_c_eff` do not prove any values from the historical extrinsic insertion event or the 651 traction law.

## Frozen provenance and reproduction

Inputs: MSH SHA-256 `b98402c5a7b7af36b8e5e8a2a90c7db689cb1af86a089c8e1b9a0dbf07ba2638`; material SHA-256 `c300c3091ac323eba756f2d7bfbf55e888c68d6b10f5fcb9763c0e6b6194305f`; manifest SHA-256 `21b563069550e2091474fe94877bd3ceaae094ac4c92aed3900ecdd94702dd68`; prior probe driver SHA-256 `65fb7e4a4590d1cbb0a015a4ff35a0b5261528926394c5a2cac692cf76e42be6`. The actual runtime wheel is `akantu-5.0.7.post1`, SHA-256 `f618c31bb4cdf23d86e74a10f7e94dac3dae7b0ae2cfd057d352c1b156c0654b`. The inspected local source cache is tagged v5.0.7 and hashed separately in the machine audit; it is not represented as the wheel's exact source. The cached `material_cohesive_linear_inline_impl.hh` is absent.

The new script SHA-256 is `097acee662d3fd890543e17db4f3fe955aa3007d3e2c68e359f812ea68c8a4a3`. Run it using the exact WSL commands in the paired experiment [README](../../../experiments/continuation_20261016/extrinsic_capture_preflight/README.md). The two machine results are `extrinsic_capture.json` and `intrinsic_capability_capture.json`; the copied arrays and stdout logs are alongside them. Each process was limited to 60 seconds; the complete paired data folder is about 0.12 MiB, below the 16 MiB total output cap.
