# Moving-front provider: one bounded source-level check plan

**Status: planned only.** The provider, driver, and both front states have not been imported or executed. No MLS batch, provider sampling, operator assembly, or solve has been run. `moving_front_check_driver.py` refuses execution until the matching `moving_front_check_spec.json` is changed to `READY`; after that it claims a one-run marker before doing work.

## Frozen inputs and expected behavior

Use only the exact pair in `moving_front_cases.json`: the finite front at `x=0` and the same front shifted to `x=0.0005 m`, within the same closed box and on the same 100 material-rest source coordinates. Load `local_front_support_safe_v2_20261006/support_topology.py` under the hash-specific private module name already specified in the source. Reclassify every source for each front and build a fresh BASE/H/shifted-F1 registry. Do not use the old manifest's role or certificate fields. Do not use old coefficients.

Before constructing either provider, record the current roles and counts; require 100 sources, no UNKNOWN source roles, one BASE per source, at most one enrichment per source, and no more than 375 scalar DOFs. Compute the current core radius from the current safe-v2 roles and current finite segment:

`r_core = min(distance(source center, current finite F) for each non-F1 source) - h`.

For each state use exactly `r1=min(0.00025 m, 0.25*r_core)` and `r2=min(0.00050 m, 0.50*r_core)`, with four diagonal x-y directions at each radius. Both radii are strictly inside a positive core by construction. If the measured `r_core` is nonpositive, or the two computed radii are not distinct positive machine values, block the core check without searching for replacement radii or changing the algorithm.

## Single future batch per state

Each future provider batch contains exactly 25 unique point/trace identities:

| Rows | Points | Check purpose |
|---:|---|---|
| 8 | Two radii, four diagonal directions in the current-front x-y plane, `z=0` | Core field and gradient against the current-front analytic manufactured field |
| 4 | `x=front_x-0.008` and `-0.004 m`, each with explicit PLUS and MINUS | Exact finite-Γ side values; no epsilon-offset trace approximation |
| 2 | `(front_x-0.002, ±0.0005, 0) m` | Record collar error only; no exactness criterion |
| 3 | `(front_x+0.004, y, 0) m`, with `y=-1e-8, 0, +1e-8 m` | Check both open sides and the ahead-of-front plane point; verify all H contributions vanish and the BASE/F1 rows approach continuously |
| 8 | Current front interior; both finite-front endpoints; finite-Γ x-start edge; finite-Γ z-edge; Γ interior with missing side; a point at `0.5*tolerance` from Γ; and one point outside the closed box | Each must fail closed as UNKNOWN for its stated front/edge/tolerance/domain reason |

Use the existing explicit PLUS/MINUS API for the four exact traces. Check that the basis rows are finite and have the correct analytic one-sided H/F1 branch values; do not claim that the complete field is exact at these face locations. The two collar rows are observations: report their measured values/errors without pass/fail against exact full-field reproduction. At each current-front analytic sample, use the exact analytic reference and its gradient. For this manufactured-only check, coefficients are constructed afresh from the current basis: `c=A*sqrt(R0)`, each current F1 source gets `BASE_i=c*Fbar_i` and `FRONT_F1_SHIFTED_i=c`, all other BASE and H coefficients are zero. Check only points strictly inside that state's measured core. This verifies the local shifted basis algebra and MLS partition there; it says nothing about the original fixed load or a mechanical solution.

## Identity and topology controls

Across the two states, require different surface-geometry hashes and different basis-generation hashes, and confirm every sample carries its own case, front, geometry hash, basis hash, source-coordinate digest, classifier digest, and canonical-ID digest. Construct a provenance record from the first state and require the second state to reject it. Do not remap that coefficient record; production coefficients need mechanical re-projection or re-solution after the basis change.

For each state, add two classifier-only synthetic controls that do not enter the source set or the MLS batch: a source center exactly `h` from the current finite front, and a center `0.5*tolerance` beyond that tangent distance. Safe-v2 must classify both as UNKNOWN. These four controls cover exact and near tangency at each front without moving a material source.

The provider batches total 50 point/trace rows. Reclassifying the 100 fixed sources twice contributes 200 topology records; the four classifier-only tangent controls bring the planned total to 254 records, within the 256-record/query ceiling. Keep each batch at 25 rows and stop the entire check at 30 seconds. No broad parameter sweep, operator, stiffness assembly, or solver is part of this check.

## Acceptance and limitations

The frozen tolerances are absolute affine displacement `1e-10 m`, affine gradient `1e-9` dimensionless, current-front core displacement `1e-12 m`, core gradient `1e-9` dimensionless, trace basis value `1e-12` dimensionless, and trace basis gradient `1e-9 m^-1`. Accept the check only if the affine BASE patch and current-front core manufactured field meet their respective limits; exact-side rows match the explicit one-sided basis formulas; all three ahead rows have zero H values/gradients and their finite-sample BASE/F1 variation is recorded; every edge/tolerance/outside negative is UNKNOWN; both tangency controls per state are UNKNOWN; semantic IDs remain stable under a permutation of already-classified roles; all sample tags match their state; and cross-front provenance is rejected. The ahead variation is a finite-sample diagnostic, not a global continuity proof. Collar errors are observations without a threshold. This source-only plan makes no measured error or pass claim and does not qualify fracture growth.

The one-shot driver caps the accounting at 254 records: 200 source-role classifications, four topology-only tangency controls, and 50 provider point/trace identities. It uses one provider/MLS batch per state and no extra MLS call for the affine, core, trace, tag, or permutation checks. Its 30-second limit includes imports; its 2 MiB limit covers diagnostic-owned arrays, a fixed 512 KiB transient reserve, and the serialized result (not the Python interpreter/import baseline). The result is limited to 256 KiB and is created exclusively once. No outcome is claimed in this source-only plan.

No outcome is claimed in this source-only plan. In particular, support coverage and a nonzero local core do not establish global blending consistency, crack-growth accuracy, or a solved mechanical state.
