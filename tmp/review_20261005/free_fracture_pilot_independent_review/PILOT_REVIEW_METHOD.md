# Static source and adapter review (pre-output)

Status: **prepared, output audit not run**. The final receipt and summary hash were checked, and a read-only input inventory streamed SHA-256 hashes and recorded byte counts for the ten flat result files. No result CSV/TOML values or binary arrays have been parsed or summarized by the independent arithmetic adapter; no solver was run.

## Frozen source bindings

- `run_free_fracture_pilot.jl`: SHA-256 `9922d52c13468df5152d88c1fdb1b7c9fc6ce473a1a3a9e56507aaf37479f50b`
- `fixed_case.toml`: SHA-256 `45b9f8bb5d10522fa2e6fa91ca8d944e2d6af27396034b489e27a026783bbedd`
- `run_once.ps1`: SHA-256 `97434a0f5b17dd93ae7a0064f98515b4651aad583b0a039b7907c646c89ba69d`
- `source_freeze.toml`: SHA-256 `1f4c3e9b10a607d303bbc3435bf22ae8b591ce8f6b5651aacb6cba3c3bfd7217`; schema v1, source-only freeze, 84 declared pins.

The fixed case is the 4,500-point, 20×15×15 NoCorrection OSB pilot with 96 steps / 97 captured states at 3.75e-8 s, complete-horizon free interior, end grips, strict native CriticalStretch failure, and bounded one-chunk output. Static inspection found no new blocking inconsistency in the driver, fixed case, launcher, or pinned native recorder path. This is not a claim that the completed numeric output is correct.

## Independent adapter

`pilot_independent_review.py` is confined to this review directory. Its result path is fixed to the pilot's `run_attempt_01/results` and it does not launch Julia or modify the experiment. It now requires both an explicit `--root-authorized` switch and the author-provided final summary SHA-256. The summary hash must match the completed launcher receipt before any result files are opened; the adapter then checks the exact flat file inventory, byte lengths, SHA-256 values, typed-binary sidecars, array order/shapes/units, and source-freeze pins.

The independent reconstruction covers the strict bond-failure transition rule from point coordinates, canonical roster irreversibility, point active-neighbor counts/damage, union-find connectivity after every distinct active-roster transition (including isolated points), minimum-point-ID labels, component mass/bounds, graph-history rows, 97-state secant-grip path, free VelocityVerlet updates, energy and force-density integration, BC half-step impulse/work, and the full-step momentum identity. Repeated states with byte-equivalent active roster reuse the already established partition; their binary labels are still compared point by point, and selected-state bounds are refreshed from that state's coordinates. This bounds the graph work without changing the graph result. Raw full-step KE is separated into free and grip-masked contributions; physical constrained-grip KE remains UNKNOWN. Energy closure and fracture dissipation are diagnostic/NOT COMPUTED.

The spatial/grip association analysis reports incident-pair denominators and event counts for left grip, right grip, both endpoints in the right grip, and right-grip boundary-crossing pairs. It also aggregates right-grip incident pairs and first breaks by their exact reference x-layer pairs, with first/last event steps, and reports right-grip graph membership and mass share per selected state. This specifically enables an independent check of the reported 6,810-pair/right-grip observation without treating those pairs, their disconnection, or a graph component as an internal crack or physical fragment. No crack-plane fit, fracture surface, or fabricated surface rendering is implemented.

The author-frozen receipt has an exact bookkeeping discrepancy: its `final_attempt_bytes` is 88,924,215, while independently summing the completed attempt tree gives 88,924,268 (+53 bytes). The receipt and original experiment are preserved unchanged. The adapter records that as a failed zero-tolerance check and continues arithmetic only when the attempt tree is within the fixed case's `max_attempt_bytes` (209,715,200), the results tree is within both the fixed case's `max_output_bytes` and source freeze's `max_experiment_output_bytes` (199,229,440), and all other byte/source/exit checks pass. It then reports `AUDITED_WITH_RECEIPT_SIZE_DISCREPANCY` only if that is the sole failed check. Native arithmetic status is a separate field; any other failed check yields overall `FAIL`. The special audited-with-discrepancy state is a completed diagnostic (process exit 0), while retaining the failed receipt check in JSON.

Root source review found and directed a narrow correction before any adapter execution: half-step kinetic-work arrays are indexed `(point, component, state)`, so the independent left/right grip masks now select the point axis (`work_by_dof[left_mask, :, :]` / `right_mask`) and sum point/component axes. A prior expression mistakenly applied the 4,500-point mask to the 3-component axis. The adapter has not been run after this correction.

When authorized and given the final summary hash, the adapter produces actual-data equal-world-mm geometry/damage with deterministically sampled first-break bond segments, plus a four-panel event/graph/mass-share/raw-energy timeline. Graph count and fractions occupy separate panels with independent units, so a value near 0.95 is visible rather than sharing an axis with cumulative event counts. No displacement magnification is used and the captions keep graph components separate from physical fragments.

`run_bounded_review.ps1` is a separate parent-authorized launcher for the adapter. It requires `-RootAuthorized` and a 64-hex final summary SHA-256, starts the bundled Python adapter without a visible window, and kills it at 110 seconds. Thus the read-only arithmetic pass has a hard 110-second child-process bound, under the 120-second review budget. The wrapper has only been parser-checked; it has not been executed.

## Validation performed here

- Compiled the adapter with the VesselForge Python 3.12 runtime's `py_compile`; syntax check passed.
- Checked the five named frozen-source hashes above.
- The metadata/hash-only input inventory binds the completed receipt (`0bae8302c765effdd2614f72a2e1802d4779aa9ea4910829a1281f063c4eb557`) and summary (`1e3396d30aa124a3371b47bcd6cc40eea4118312a18d685bbfbe3fdaa60a9c8c`); its result tree is 88,921,168 bytes across ten files and manifest SHA-256 `922353f92501004041c89f6377e2636df106b0883e5220bdb2106967b855e063`. Inventory JSON SHA-256 is `480696292e9dc58cb368ebdb35102b1d727f716a42a215a1759cf6b678fce87d`.
- Independently summed the completed attempt tree to 88,924,268 bytes, confirming the receipt's recorded count differs by 53 bytes.
- Confirmed the recorded trees fit their separate exact limits: results 88,921,168 / 199,229,440 bytes, attempt 88,924,268 / 209,715,200 bytes.
- PowerShell parser check passed for the 110-second wrapper; neither wrapper nor audit adapter was run.
- Did not invoke the arithmetic adapter, parse result values, run tests or Julia, or write to experiment/source/vendor/test/fixture/oldtmp paths.

Adapter SHA-256: `e2d2f8cf2f9d71826c2d71bf9058546ea3ea6249ac4228a4d5f73149815c6187` (104,987 bytes).
Bounded runner SHA-256: `ab74c6afbdea7ee9b938c0aea44a81243a099a529388f5dd04fe915dd1f4bc7d`.
