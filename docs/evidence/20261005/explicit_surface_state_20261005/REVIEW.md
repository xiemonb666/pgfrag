# Explicit surface state source review and freeze

Date: 2026-10-05

## Scope

Reviewed the source-only candidate at `pgfrag/experiments/explicit_surface_state_20261005/`. No production package, vendor solver, job/result schema, observed-feature code, or other agent's work was changed.

## Review result

- One authoritative reference-space triangle geometry and one oriented boundary-edge array are stored per patch. The two side records share that geometry by construction, carry independent immutable displacement arrays, and export reverse triangle winding on the negative side.
- Inputs must be exact `numpy.ndarray` objects. Dtype, dimensionality, shape, row count, and input-byte caps are checked before dtype conversion or contiguous copying. Duck-typed `__array__` objects and ndarray subclasses are rejected. Counterexample coverage includes oversized zero-stride views with the copy path guarded.
- Duplicate faces, degenerate faces, inconsistent winding, nonmanifold edges, invalid IDs/steps, invalid lineage, and over-budget inputs fail closed. The supported boundary status means only that directed boundary edges compose one oriented cycle; branched/touching/multiple-front topology stays `UNKNOWN`. State-wide surface intersections stay `UNKNOWN_NOT_EVALUATED`.
- Patch assignment names one exact patch ID and one exact fixture body/fragment pair. A relation entry is assigned fixture metadata using the stored patch area; it does not prove two fragment geometries share a physical surface and is not geometric/model R2. A same-body internal surface cannot carry a relation, including through direct state construction. Current leaf relation entries use exact patch assignments, with no descendant Cartesian product.
- Parent IDs establish metadata lineage only. Child patch containment, patch partition, global intersection, and geometric lineage of leaf patches remain `UNKNOWN` and are not verified.
- Historical fixture relations are separate from current visibility, so wear/occlusion does not delete their history.
- All outputs are fixed to `fixture_proxy.explicit_surface_state.v1` and relation namespace `fixture_proxy_relation.explicit_surface_state.v1`. The schema and serializer fix `mechanics_qualification=NOT_COUPLED`, `mechanically_active=false`, and `production_promotion=BLOCKED_NOT_COUPLED`. No external solver proof verifier or production-label path exists.
- Caps: 100,000 vertices and 200,000 triangles per patch; 4,096 patches and 256 MiB of stored patch arrays per state; 100,000 assignment/history rows. These are storage guards, not resolution claims.

## Static checks

- Python AST parsing of the implementation and the candidate test source: passed.
- JSON parsing of `state.schema.json`: passed.
- `pytest`: not run, as requested.
- Solver/native execution: not run, as requested.

## Limits retained as explicit unknowns

This is not a mechanics-coupled implementation. It does not check surface-surface intersections, patch containment across parent/child geometry, volume partition, crack-growth criteria, support visibility, displacement-jump consumption, energy balance, contact, or solver event authority. It does not fit surfaces from OSB damage and is not a reproduction of Ai & Augarde. No `SIMULATED_GEOMETRIC_GT` or production/model R2 claim is supported.

## Frozen files

SHA-256 entries for the four candidate source artifacts are recorded in `freeze_manifest.json`. The manifest intentionally excludes itself and this review text.
