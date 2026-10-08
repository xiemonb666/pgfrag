# Planned assertions (not executed)

The only positive case is a 16-point `GEOMETRY_FILTER_FIXTURE_ONLY` two-box cloud. The zero-degree variant is intentionally invalid and exists solely to test the gate.

| ID | Planned assertion | Evidence source |
|---|---|---|
| A01 | Fixture is exactly 16 rows and within 125-point cap; directed candidate roster stays below 20,000. | Julia source guards; fixed CSV |
| A02 | Point CSV bytes, IDs, exact row order, mm coordinates, identity frame, and mm→m body positions agree. | `InitialBondPolicy` constructor and body-scope validation; Julia source check with `atol=1e-12 m`, `rtol=0` |
| A03 | Every native point horizon is exactly `0.0018 m`; no point-specific horizon variation is accepted. | `get_point_param(body, :δ, i) == model.horizon_m == HORIZON_M` |
| A04 | Native raw `find_bonds` roster has exactly 120 directed candidates and exactly equals the independent all-ordered-pairs Euclidean predicate `distance2 <= δ^2`. | Actual native raw search vs independently built directed `(external_source, external_target)` set and source-degree vector |
| A05 | Transport state ledger equals an independent fixture-only two-component pair classification; canonical external IDs remain distinct from native row ordering. | Parsed `decision_by_pair` vs independent candidate/state map |
| A06 | Geometry filter admission receipt reproduces raw candidates and filtered bond list; raw=120 directed/60 canonical, included=112 directed/56 canonical, excluded=4 canonical, unknown=0. | `BondSystem.initial_bond_admission` |
| A07 | Filter lifecycle is raw candidate search → reciprocal geometry filter → physical precrack filter → halo/ChunkHandler/localization/bond IDs → correction/kernels. The fixture has no precracks. | Pinned overlay `bond_system.jl`; graph is constructed from actual `BondSystem` |
| A08 | Every included canonical pair appears in both native directed orientations; every excluded candidate appears in neither orientation. | Readback from `system.bonds` grouped by `system.bond_ids` |
| A09 | Actual native post-mask graph contains 112 directed edges, degree 7 for each of 16 points, identity single-chunk localization, and a unique contiguous native bond-ID ledger. | `BondSystem` fields / ChunkHandler mapping |
| A10 | Geometry exclusions remain `geometry_exclusions_are_physics_events == false`; static transport remains `native_graph_readback == NOT_RUN`. | Native admission receipt; no physics/event code is called |
| A11 | Pinned native OSB `calc_weighted_volume` read-only method returns finite positive values equal within `rtol=1e-12`, `atol=0` to an independent sum over expected included edges with default `NoCorrection` and all-active scratch state. | Actual `calc_weighted_volume` vs direct expected-mask formula |
| A12 | A second transport that forces all included pairs incident to external ID 103 out of the mask is statically `BLOCKED_ZERO_DEGREE`; constructing its native `BondSystem` throws before graph admission. | Negative transport gate and `_validate_initial_policy_scope` call order |
| A13 | Driver wall clock is at most 60 seconds and result artifact is at most 20 MiB. | Driver guards; future external timeout still required |

No assertion here proves a full ThreadsBodyDataHandler/BodyChunk initialization, production CAD adapter, classifier execution, native physics initialization, force result, solver result, crack, or physical failure. No Julia source/test was run while preparing this plan.
