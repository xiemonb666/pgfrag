# SOURCE_ONLY: one frozen VF solid to a bounded material-point candidate

Status: design/source preparation only. This worker is not connected to pgfrag jobs, backend.py, the Julia runtime, or the public package API. It has not been executed; no CAD boolean, point generation, pytest, or solver run was performed.

## Frozen input and existing contract

The only allowed source asset is the copied VF fixture:

vesselforge/tmp/review_20261005/geometry_resume/migration_copy/cylinder_asset_run01

Pinned source manifest SHA-256: df781dc6c818512af8d66b164288aa0d5991097ea0166b85556119de2c9fba41

Pinned canonical STEP SHA-256: 5c68bbc29b2b0bc684cc6c551c2769872c4c490d588b9a3ea2cbb98e7f74135e

The existing PG package contract remains forge.asset.v1, as checked in pgfrag/src/pgfrag/assets.py. The existing forge_asset_surface job kind remains unchanged and currently returns UNSUPPORTED before launching a solver. This isolated candidate data format is not a new job kind and does not alter that behavior.

The input manifest declares mm, kg, seconds, right_handed_z_up, one canonical STEP, surface arrays, field arrays, mechanics placeholders, and geometry_only scope. Its mechanics remain UNVALIDATED_PROTOTYPE. Face semantics are mother-only regions. The separate migration receipt proves two finite wall-reference stations on the copied mother; it does not establish global minimum wall thickness, fragment mapping, or physics.

## Sampling rule

The one authorized candidate configuration is fixed in source: h=0.5 mm and explicit diagnostic ROI x=[7.5,10.5], y=[-1,1], z=[4,6] mm. The ROI intentionally extends beyond the outer x bound of 10 mm so the ledger contains exterior/void-only cells as well as material intersections. It is not clipped to the source bounds. Its 3 x 2 x 2 mm AABB divides exactly into 6 x 4 x 4 = 96 cells. No CLI or worker fallback may alter h or the ROI. This run is a local mechanism probe, not a whole-mother discretization.

Enumerate a deterministic Cartesian grid over that region. For every candidate cell, compute the OpenCascade B-rep solid-cell intersection common of the canonical material solid and the cell box. Do not use triangle center tests, inside-at-center plus h cubed weights, or a surface projection fallback. For each positive-volume connected solid in the common:

1. Retain the raw B-rep common volume, connected-component volume, and B-rep center of mass.
2. Require the center of mass to classify exactly IN against the canonical source solid at the declared 1e-7 mm classifier tolerance. ON, OUT, exception, or UNKNOWN blocks the whole candidate.
3. Emit one point per connected common solid with positive volume. Preserve its source solid index, grid-cell index, cell and component AABBs, exact mm coordinates/volume, and SI conversions.
4. Leave source semantic region, mother triangle, and surface-parent mapping UNKNOWN_NOT_PROJECTED. Interior volume points are not assigned to a single surface face by assumption.

The reference volume is an independent second B-rep solid-cell intersection common of the source solid with the full ROI. The sum of all per-cell connected-component volumes must agree with that ROI B-rep intersection volume within max(1e-9 mm^3, 1e-8 times the reference volume). Retain both the ROI common's component count/volumes and every cell's zero/positive volume, components, centroid result, and point provenance. If the ROI common has more than 16 connected components, block and preserve the bounded receipt with the over-cap count; if it has positive volume but no connected solids, or a cell result cannot be reconciled, block. This is a local B-rep volume reconciliation only, not whole-mother conservation, physical conservation, or discretization convergence.

Position conversion is mm to m (1e-3); volume conversion is mm^3 to m^3 (1e-9); point mass is declared rho_kg_m3 times the converted B-rep intersection component volume. No h^3 mass is assigned. Output ordering is deterministic by k, j, i, then component bounds/centroid; point IDs are one-based in that order.

## Bounds

Hard caps are at most 256 candidate grid cells, 500 emitted points, 16 connected common solids per cell, 512 common faces per cell, and 20 MiB total output. One isolated worker process is terminated after 180 seconds. The source worker rechecks elapsed time between CAD operations; the outer process timeout is the hard guard around an uninterruptible OCC call.

The migrated cylinder has geometry-only B-rep bounds [-10,10] x [-10,10] x [0,10] mm and a declared nominal wall parameter of 2 mm. The existing shared bounds helper uses OCP.BRepBndLib.AddOptimal_s with use_triangulation=False and use_shape_tolerance=False. These metadata are not a global thickness certificate. The previous full-mother budget was rejected: h=2.5 mm yields 256 bbox cells with h/t=1.25; h=3 mm yields 196 bbox cells with h/t=1.5. Passing the count cap by using cells at least as large as the nominal wall does not meet the intended geometry-resolution gate and is not the pilot.

The approved diagnostic ROI's explicit h/t ratio is 0.25 against the nominal wall parameter, or four cells across that nominal thickness. The ROI material volume estimated analytically from the cylindrical radii is about 8.0168 mm^3; this is only a planning estimate, not the CAD oracle. The exact independent ROI common and 96 cell intersections remain to be computed by a later owner-authorized run.

An analytic estimate from the frozen STEP face inventory is a straight annular cylinder with outer radius 10 mm, inner radius 8 mm, and height 10 mm, giving 360π = 1,130.973355 mm^3 (about 28.27% of the 4,000 mm^3 bounding box). For the approved partial ROI, 2 times the integral from -1 to 1 of sqrt(10^2-y^2)-sqrt(8^2-y^2) gives about 8.0168 mm^3. With declared density 2,100 kg/m^3 this would correspond to about 0.00001684 kg in the ROI. These are formula-based planning estimates, not B-rep common volumes or physical volume measurements.

## Candidate output

The worker writes only under this experiment's tmp directory:

- ingress_report.json: input hashes, region, explicit h, volume reconciliation, runtime/version, limits, and UNKNOWN/qualification fields.
- material_points.json: emitted only if every common component centroid is IN and local volume reconciliation passes; ordered positions in m, positive intersection volumes in m^3, masses in kg, source cell/component identifiers, AABBs, and explicit unknown surface mapping.
- cell_ledger.partial.json and cell_ledger.json: zero/positive/blocked cell status, every available raw common/component volume, centroid classification, AABBs, and linked candidate point provenance.
- cell_in_progress.json: atomically refreshed before and during each cell/component operation, with source pins, cell/component AABBs, raw volumes, centroid/classifier state, and the hash of the preceding partial ledger. A timeout therefore retains the last completed ledger plus a bounded record of the active cell.
- parent_execution_receipt.json and worker_failure.json: parent child-exit/timeout/retained-byte receipt and worker failure detail where relevant. Failed, timed-out, and not-started runs retain their stage; the code does not delete a failed stage. On any failed outcome, partial-ledger status is changed to BLOCKED and records the parent/worker failure summary.

The pointset is an ingress diagnostic candidate only. A failed gate retains its bounded raw stage and is never promoted to material_points.json. It is not imported by backend.py; candidate bond-pair topology, native graph initialization, loading, integration, fracture events, energy closure, fragments, and labels remain separate gates. Do not run this worker until root reviews this source revision.
