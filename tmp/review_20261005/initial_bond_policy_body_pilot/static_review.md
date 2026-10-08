# Static review — initial bond policy native graph pilot

**Review state:** source prepared, static only, awaiting root review. The Julia driver and Python policy generator have not been executed. This review is limited to the newly isolated experiment files and source-order inspection of the already pinned fork/depot sources; the original project, manifest, vendor tree, depot, overlay, and `fork_run_03` were not changed.

## Native hook and horizon finding

The selected seam is the real `BondSystem` path in `initial_bond_policy_native_20261005/overlay/src/discretization/bond_system.jl`:

1. `find_bonds` creates `GridNeighborhoodSearch(search_radius=maximum_horizon(body))` and fills the raw directed per-point bond roster (lines 117–126).
2. With policy enabled, `_filter_initial_policy_candidates!` runs on the complete one-chunk raw ledger (lines 128–130).
3. `filter_bonds!` applies the ordinary physical precrack filters afterward (line 132).
4. `get_bond_data` only then derives halo points, builds `ChunkHandler`, localizes neighbors, and constructs `bond_ids` (lines 83–86).
5. `BondSystem` creates correction and kernels only after `get_bond_data` returns (lines 56–66).

This is before precrack filtering and before halo/localizer/bond-ID/correction/kernel stages. For the fixed fixture the precrack list is empty, so the observed admitted graph can only differ from the raw graph due to the initial geometry mask. Geometry exclusions remain admission decisions, not damage/failure events.

The pinned native PointNeighbors source says “smaller than” in its API prose, but the selected `GridNeighborhoodSearch` implementation uses `distance2 <= search_radius^2` at `third_party/pgfrag-runtime/depot/packages/PointNeighbors/STFG7/src/nhs_grid.jl:549`. `find_bonds!` passes the per-point material horizon as that `search_radius` in the pinned overlay. The pilot therefore uses inclusive `<=` with no horizon padding/tolerance. No fixture pair is exactly at 1.8 mm; the cross-gap distances adjacent to the threshold are `sqrt(3.25) mm ≈ 1.8028 mm` and are outside.

## Fixture and expected graph

The two boxes are `{x=0,1}×{y=0,1}×{z=0,1}` and `{x=2.5,3.5}×{y=0,1}×{z=0,1}` in millimetres. At 1.8 mm each box is a complete K8 subgraph (28 pairs, max distance `sqrt(3) mm`). Exactly four cross-gap same-`y,z` candidates lie at 1.5 mm. All other cross-gap distances exceed 1.8 mm. The expected raw roster is 60 canonical pairs / 120 directed candidates. The normal toy mask includes 56 within-box pairs and excludes four cross-gap pairs, leaving 112 directed edges and degree 7 for each point.

The external IDs are deliberately permuted relative to native CSV row indices. Every native roster comparison maps row indices through the explicit native↔external identity table and compares complete directed sets. It is not a count-only or only-listed-pairs check. The four excluded pairs are checked absent in both directions; included pairs are checked present in both directions; native `bond_ids` are read and verified as a unique contiguous ledger.

The negative bundle deliberately forces the seven otherwise included pairs incident to external ID 103 to `OUT_OR_VOID`, on top of the four actual toy gap exclusions. It has 49 included pairs and exactly one zero-degree point. This is a synthetic invalid-policy input whose only intended use is to confirm `BLOCKED_ZERO_DEGREE` fails before a native `BondSystem` can be admitted. Its classification is not a geometry claim.

## Weighted-volume scope

The pinned OSB source computes weighted volume on demand in `physics/ordinary_state_based.jl:207–219`; it is not a persistent pre-solver `weighted_volume` field. The pilot prepares an all-active `OSBStorage` scratch value and calls the actual native `calc_weighted_volume` method directly on the admitted `BondSystem`. The independent comparison reconstructs the expected included adjacency from the fixture policy and sums the OSB terms using native positions, uniform horizon, linear kernel, volume, and `NoCorrection` (`surface_correction_factor(::NoCorrection, _) == 1`). It does not call force density, create a time solver, construct `BodyChunk`/`ThreadsBodyDataHandler`, or enter physical fracture initialization. Thus the claim is specifically a native weighted-volume formula readback on a native post-mask graph, not full OSB simulation initialization.

The decision not to construct `ThreadsBodyDataHandler` or `BodyChunk` preserves the no-solver/no-force scope. `PointDecomposition(body,1)` plus direct native `BondSystem` construction exercises the exact one-body/one-chunk graph seam, but does not establish full Threads handler assembly. If root requires that later, it needs a separately approved source/run step because the pinned handler constructs `BodyChunk`s using an `AbstractTimeSolver` and calls chunk initialization.

## Evidence and limitations

- Generator produces fixture-only JSON/TOML bundles and distinct JSON/TOML expected hash values; the Julia constructor requires those values separately from the CLI and re-reads exact bytes.
- The generated audit includes a CAD-named field only because the pinned transport schema requires it. Its file is visibly named and documented as a non-CAD marker; status `READY` and accuracy `0.0` are schema placeholders. The bridge validator verifies file-byte identity but does not establish scientific classifier execution. The fixture compiler result and interval metadata are not pytest or CAD evidence.
- The application is always tagged `GEOMETRY_FILTER_FIXTURE_ONLY`. `native_graph_readback` in the static policy remains `NOT_RUN`; a future driver result is a separate fixture-run receipt.
- Output bounds are guarded in source for 125 points, 20,000 directed candidates, 60 seconds, and 20 MiB of result artifact size. Total Julia process RSS is not enforced by source; if 20 MiB refers to full process memory, a root-approved external RSS constraint or budget clarification is needed before any run.
- The future driver run must use the verified disposable `fork_run_03/project` environment without changing its manifest, source tree, depot, package pins, or original repository state. This review does not execute it.

## Status

Static source review: prepared. Julia import, policy generation, raw native search, native BondSystem graph readback, weighted-volume call, blocked zero-degree probe, CAD, solver, force, and physical events: `NOT_RUN`.
