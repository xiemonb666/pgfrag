# Retained pilot adapter review

This adapter will read only the frozen retained stage and the frozen supervisor summary listed in PILOT_AUDIT_INPUT_FREEZE.json. The summary is only a provenance/finalization/promotion cross-check; expected geometry is computed by analytic_roi_reference.py from the fixed R10/R8 cylinder and the root ROI.

The retained stage schema is:
- material_points.json: point_id, position_m, cell_volume_m3, mass_kg, and source.grid_cell_index/grid_ijk/component indices/AABBs.
- cell_ledger.json: 96 cells with raw_common_volume_mm3 and component_records. Each positive component carries raw_volume_mm3, centroid_mm, centroid_source_classification, candidate_point_id, and candidate_point_status.
- parent_execution_receipt.json and worker_failure.json: parent and worker failure evidence.
- cell_ledger.partial.json and cell_in_progress.json: retained blocked ledger/progress.

Adapter checks:
1. Verify every frozen file byte count and SHA-256 before parsing; independently compare its file inventory with the stage list in the summary for provenance only.
2. Generate the 96 analytic cell references and a separate direct full-ROI integral using the reviewed pure-stdlib formula module. Compare analytic direct ROI to analytic sum, actual reference ROI B-rep volume, actual sum of raw cell commons, sum of component volumes, and point volume sum. Report absolute deltas and the frozen tolerance; do not use producer summary volumes as expected values.
3. Cross-map each material point to the raw component via point ID and cell/component keys. Check SI volume conversion by multiplying cell_volume_m3 by 1e9 and comparing to raw mm3; check mass_kg against rho×m3 with rho=2100 kg/m3. Check position_m×1e3 against the raw component centroid_mm and its cell/component bounds.
4. Count raw component classifier states from the ledger and require every emitted point to map to exactly one positive component marked IN and CANDIDATE_CENTROID_IN. This producer classifier remains a separate recorded check; it is not used as the independent geometry oracle. Missing classifier or source-map evidence remains UNKNOWN.
5. For every positive analytic cell, independently compare the delivered volume-weighted centroid to the analytic cell centroid on x, y, and z. Require all coordinates and deltas to be finite and gate the maximum absolute axis delta at the fixed `1e-6 mm` tolerance. Report the observed maximum delta overall and per axis; do not adjust the tolerance based on results.
6. Independently check each delivered point/centroid against the analytic source annulus using strict `8 mm < sqrt(x²+y²) < 10 mm`, the source z slab `[0, 10] mm`, and its mapped grid-cell z slab. Preserve this boolean in the point CSV and record offending point IDs separately from producer classifier outcomes.
7. Build a cell audit CSV and point audit CSV, then two plots: equal-aspect XY analytic fill plus actual centroids in four z slabs; and analytic/actual cell volume, rhoV mass, and signed volume delta. Add titles that state the report-finalization failure and nonpromotion.

Do not infer successful production handoff from the material_points.json candidate file. Finalization summary and parent receipt both indicate failure/nonpromotion; the geometry ledger remains useful only as a local diagnostic recovery. No run, calibration, full-mother, native integration, or physics claim is made by this adapter.
