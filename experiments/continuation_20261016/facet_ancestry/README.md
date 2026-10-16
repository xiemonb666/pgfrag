# Generic source-facet ancestry

`pgfrag.facet_ancestry` is a reusable, geometry-independent-of-plane lineage
helper for native T4 meshes with optional cohesive insertion. It is prepared
for later bulk-driven workflows; it does not change the existing prescribed
planar weak-seam path.

The source table is built with `build_source_facet_table(source_nodes,
source_tets, source_node_ids, source_cell_ids, source_anchor_hash, ...)`.
Coordinates must be finite, every source T4 must have positive rest
orientation and nonzero volume, node and cell IDs must be unique, and the
caller must supply the frozen 64-hex SHA-256 anchor. Every triangular face is
classified from exact source-node-ID incidence: one parent is exterior, two
parents are internal, and more than two parents fail. Internal parents must
have opposed outward windings. By default all internal facets are eligible;
pass an explicit subset of source-ID triplets or stable facet identities when
only part of the source interface is eligible. Exterior and unknown facets
cannot be marked eligible.

A stable facet identity is the frozen source-anchor digest followed by the
sorted original source node IDs: `sha256|id0,id1,id2`. It never depends on
coordinates, current node ordinals, cohesive row ordering, or a fracture
plane. Different source IDs at the same coordinates therefore remain
different facets.

`map_native_cohesive_ancestry` requires the full pre-insertion native bulk
connectivity snapshot to exactly equal source T4 connectivity, including row
order and local corner order. Current bulk rows must preserve that lineage.
Current native nodes are mapped to source node IDs only from their cell row
and local corner ordinal, and every current bulk corner must have an exact
rest-coordinate match. The helper never performs nearest-coordinate matching.
Each current `(N, 6)` cohesive row is paired by the source IDs on its two
three-node sides, then checked against the unique eligible internal source
facet, its two original parent cells, and one current native bulk incidence
per side. An empty `(0, 6)` cohesive array is valid and reports all eligible
facets as missing from that snapshot.

Snapshot reports separate `active_facet_identities` from
`missing_eligible_facet_identities`; “active” describes cohesive connectivity
present in the supplied snapshot. `event_delta(previous_ids, current_ids)`
reports newly inserted and persistent identities and fails if an identity is
deleted, duplicated/reused within a snapshot, or mixed across source anchors.
The output's `lineage_scope` and `lineage_assumptions` make clear that
post-insertion native row/local-corner identity depends on the native backend
preserving that correspondence. Exact coordinate checks corroborate the
mapping but do not independently prove persistent cell IDs for
indistinguishable coincident cells. Stable facet IDs do not establish native
quadrature-point ordinal stability across insertion. This helper makes no
material-state, tangent, fracture, source-geometry, or solver-behavior claim
and does not stitch or weld gaps.

## Reproduction

From `pgfrag/` in PowerShell, with the project environment that provides
NumPy and pytest:

```powershell
$env:PYTHONPATH = "src"
python -m pytest tests/test_facet_ancestry.py -q
python experiments/continuation_20261016/facet_ancestry/demo_constructed_case.py
python experiments/continuation_20261016/facet_ancestry/run_saved_96t4_case.py
```

The demo builds two positively oriented T4s sharing one internal triangle,
checks an empty partial-insertion snapshot, duplicates the shared nodes in one
current cell, and reports the resulting stable facet event. Its receipt is
`constructed_case_receipt.json`. This is a constructed algorithm-validation
case only: no solver is run and no physical fracture is claimed.

`run_saved_96t4_case.py` is a separate read-only check against the historical
96-T4 source MSH, `input_native_initial.npz`, and saved step 650 B/C/D dense
states in `relicweave-synlab/tmp/review_20261004/implementation/cohesive_event_probe_20261004/coarse654_20261004_01/`.
It requires the saved initial native connectivity and rest coordinates to
match the parsed source exactly, selects the eight declared interface Tri3
facets by original MSH node tags, and validates the saved C/D parent pairs.
The C and D saved cohesive rows have their recorded order compared; the
runner also reverses D's rows in memory to exercise identity invariance while
leaving all historical files untouched. It writes only
`saved_96t4_run_receipt.json` and `input_hashes.json` here. This is connectivity
ancestry evidence only, with no new solver run or physical-fracture claim.
