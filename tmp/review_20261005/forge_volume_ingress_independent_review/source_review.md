# Forge volume ingress independent source review

Status: **snapshot review only; not cleared for an actual run**. The author was still changing the worker during this review. Findings below are pinned to the observed snapshot hashes and are not claims about a later final freeze.

## Snapshot reviewed

- `pgfrag/experiments/forge_volume_ingress_20261005/worker.py` — SHA-256 `9819b66d7659dc61b29a34b038b2e5c8f4863d5635662faeea97effb0fc06034` (observed mtime 2026-10-05 05:30:24 UTC)
- `pgfrag/experiments/forge_volume_ingress_20261005/CONTRACT.md` — SHA-256 `845bd7957b34549e5aef1af963bce8f995dc6b72160187d8800efe75c601ce23` (observed mtime 2026-10-05 05:28:54 UTC)
- `pgfrag/experiments/forge_volume_ingress_20261005/tmp/plan.json` — SHA-256 `901ba1538d424cfb3879e59ff69bf5ee52f28e95423517d2bd8606a6d91d0a75` (observed mtime 2026-10-05 05:22:34 UTC)

The worker/contract hashes are the reviewed snapshot, not a final freeze pin. No repeated live-source hash polling was done after root directed that the author was still editing.

## Findings

1. **Run blocker — failed-stage retention and receipts are contradicted by the worker.** In the reviewed worker snapshot, `run_bounded` deletes the staging directory through `_remove_stage` on timeout, missing child receipt, and child failure/blocked status (worker.py:990–1034). The contract says failed raw stages are retained and lists `parent_execution_receipt.json` and `worker_failure.json` (CONTRACT.md:51), but the reviewed worker has no code that writes either file. The 64 KiB parent-receipt reserve only reduces the child write budget; it does not create a receipt. Thus the cell partial ledger written during work is discarded on these failure paths. This conflicts with the requested preservation of failed-stage/raw-cell evidence and must be reconciled in the final pinned source before an actual run.

2. **Plan drift — the plan JSON is stale relative to the root-selected configuration.** The reviewed worker and contract freeze h=0.5 mm and ROI x=[7.5,10.5], y=[-1,1], z=[4,6] mm, with no clipping and exactly 6×4×4=96 cells (worker.py:35–36, 397–408; CONTRACT.md:12–16). The plan JSON still describes the old full-bounds budget and says there is no owner-selected h or ROI (plan.json:49, 111–114). Update the plan artifact or otherwise make the approved configuration provenance consistent before treating this set of files as one reviewable package.

3. **Runtime uncertainty — empty Boolean common behavior is not established by source inspection.** The worker inspects each cell common by calling `Volume()`, `Solids()`, `Faces()`, and `isValid()`; it records `ZERO_INTERSECTION` only for a zero-volume, no-solid, valid result (worker.py:517–575, 726–734). The source path is conservative if those calls succeed, but no CAD execution was authorized, so whether OpenCascade returns that usable shape state for void-only cells is unverified. An actual small run after final freeze must establish it.

## Static checks that passed

- The reviewed source does not clip the explicit ROI to source bounds. It freezes the requested root ROI and h, enforces exact axis divisibility, and requires 96 cells. The ROI intentionally extends past the x=10 mm source bound, leaving exterior cells in the ledger.
- The shared `vesselforge.brep_identity.brep_bounds_mm` helper explicitly calls `BRepBndLib.AddOptimal_s(shape, box, False, False)`, rejects void/whole/nonfinite bounds, and emits the matching policy record. This satisfies the stated geometry-only bounds policy in the reviewed snapshot.
- The worker uses exact B-rep source-solid/cell commons, checks component and ROI volumes against the declared absolute/relative tolerance, sorts cells k/j/i and components by bounds/centroid, and preserves raw per-cell data in its live partial ledger. Tie behavior for indistinguishable component sort keys was not runtime-tested.
- The path records source hashes before and after computation and refuses changed frozen inputs. Source positions are in mm; the importer call relies on CadQuery's default STEP unit of MM. Position conversion is ×1e-3 to m, B-rep volume is ×1e-9 to m³, and mass is declared density times converted B-rep volume, with no h³ approximation. The exact installed CadQuery/OCP runtime API behavior remains unverified in this review.
- Candidate centroids are obtained through `component.Center()` and classified against the source solid with the declared strict IN-only gate. The source also keeps surface semantics/parent mapping unknown and labels mechanics as unvalidated/not connected.

## Scope and next gate

This was a source-only audit. No CAD worker, tests, native solver, full suite, or production integration ran. No source, contract, plan, Git state, or asset file was changed by this audit. Do not treat this review as clearance to run. Re-review only the final worker/contract/plan pins after the writer freezes them, then root can decide whether to authorize a bounded runtime probe.