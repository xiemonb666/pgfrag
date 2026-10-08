# Final source review — forge volume ingress

Review result: **source is ready for one root-authorized, local 96-cell CAD geometry probe**. This is a source-only recommendation; no CAD command or test was run here. It does not approve full-mother discretization, calibrated mechanics, or native integration.

## Frozen pins verified

All four expected hashes match the frozen receipt:

- worker pgfrag/experiments/forge_volume_ingress_20261005/worker.py: BE02C2914B3459B727B43D9969B95EFEB07B014768AC68D9AB39BBBA8DA602D5
- contract pgfrag/experiments/forge_volume_ingress_20261005/CONTRACT.md: 7151DF72D082B0D1D728F86AA118C3F396A4D6A3E545FBF5DC76F37A7118C6F7
- plan pgfrag/experiments/forge_volume_ingress_20261005/tmp/plan.json: E7124CEB99489B858BB2D116FED5DD9C88960209DEE632F7D49D75E2FFEFC73A
- freeze receipt pgfrag/experiments/forge_volume_ingress_20261005/tmp/source_freeze.json: 3B29580E4470CFA3016968EB784CD5B1AB9C02E832A2CE1BF5DEE57C2A2A5538

The previous snapshot findings are resolved in this pin: the worker now has failure receipts and retains failed stages; the plan records the root-selected ROI and failure-evidence policy.

## Failure evidence and bounds

Worker exceptions are caught in _child_main. It marks the partial ledger and current-cell progress BLOCKED, then writes worker_failure.json with frozen input pins, stage inventory, active-cell diagnostic, and the partial-ledger hash. The child result records if this receipt could not be written.

The parent handles timeout, launch/orchestration failure, missing child result, and any non-PASS result without deleting the stage. It marks or creates a BLOCKED partial ledger, marks cell_in_progress.json, writes parent_execution_receipt.json with stop/exit details and pre-receipt inventory, and returns the retained-stage path and inventory. A compact fallback parent receipt is attempted if the full receipt does not fit. I found no stage deletion call in the worker.

The fixed pilot is h=0.5 mm over x=[7.5,10.5], y=[-1,1], z=[4,6] mm. The ROI is not clipped to source bounds. Exact-divisibility checks yield 6×4×4=96 cells and the worker rejects any other ROI or h. Candidate cells, solids, points, faces, output bytes, and elapsed time all have declared caps.

The worker reserves 64 KiB below the 20 MiB stage cap for parent evidence. Normal JSON writes check the current non-temporary files plus the new payload; the parent receipt can use the reserved space, and the stage inventory records byte counts and hashes. One narrow hard-cap caveat remains: _write_stage_json excludes .tmp files from its budget preflight, while _stage_inventory includes them. If a child is killed during an atomic temporary-file write, the retained stage may contain a partial .tmp file, and the parent does not reject a final inventory above 20 MiB. The fixed 96-cell/500-point caps make this unlikely for this local probe, but strict 20 MiB enforcement under interruption is not proven by source alone.

On timeout, the parent tries terminate, then kill, and records whether the child remained alive. If both stop attempts fail, it proceeds to mark and inventory the stage despite child_alive_after_stop=true; that rare path could race with a still-running child. This does not block a single bounded attempt, but the receipt must be checked for the stopped-child flag and stop errors.

## Geometry gates

The grid computes exact source-solid/cell B-rep commons and a separate source-solid/ROI reference common. It checks valid finite volumes, connected-solid/component caps, strict positive component volumes, topology, centroid classification, and ROI-versus-cell volume reconciliation. Exterior cells remain in the grid and must produce an inspectable zero intersection or block the run.

component.Center() supplies each B-rep component centroid. The classifier accepts only TopAbs_IN at 1e-7 mm; ON, OUT, and UNKNOWN block that candidate. The source contains no centroid projection/rescue path. A candidate is emitted only after the unmodified centroid classifies IN. Surface-parent projection remains unknown.

## Runtime boundary and recommendation

The unresolved empty-common behavior is appropriately left to the first bounded CAD execution: source review cannot establish how this installed OpenCascade runtime reports void-only cell commons. That is the purpose of the one 96-cell geometry probe, not a reason to require full-mother coverage or calibrated material parameters beforehand.

Proceed with one isolated run only if root authorizes it. Review the resulting retained stage/receipt, exact 96-cell count, zero-cell records, strict centroid states, and local volume reconciliation before interpreting any pointset. A pass would establish only this local B-rep ingress mechanism; it would not establish physical calibration, whole-mother conservation, convergence, fracture behavior, or production/native integration.

## Execution record

Read-only source and hash review only. No CAD kernel, pytest, solver, full suite, production file, source asset, or Git command was used. The earlier source_review.md and source_review.json were left intact; this final review is an additional artifact.