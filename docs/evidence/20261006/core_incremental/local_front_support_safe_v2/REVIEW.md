# Safe-v2 finite-support topology patch review

Status: `SOURCE_ONLY / FOCUSED_GEOMETRY_CHECK_PASS / NOT_COUPLED`. The implementation is isolated under `pgfrag/experiments/local_front_support_safe_v2_20261006/`. No frozen provider, old classifier, old fixed case, or old manifest was modified. No provider manifest was generated or changed.

## Source and set convention

`support_topology.py` is copied from `pgfrag/experiments/local_front_support_contract_20261006/support_topology.py`, whose source SHA-256 before copying was `ea20e51d2b774c7e6b8b48403359b56bab1ba6f5405dba41f34813b6aec4578e`. The new module header records that provenance.

The safe-v2 module states the topology convention explicitly: the physical domain Ω is a closed box; each support is Ω intersected with a closed Euclidean ball; the finite crack is `[x_min, front) × {y_plane} × [z_min, z_max]` (closed at its start and z ends, open at the tip); and the finite front is the closed z segment at `x=front, y=y_plane`. H requires the whole box-clipped plane section to lie in Γ, two nonempty y-side pieces, and positive clearance from F. Exact equality at an included crack endpoint that is also clipped by the box is accepted. Any actual overrun and any nonzero endpoint gap no larger than the tolerance returns UNKNOWN. No epsilon shift or coordinate snap is used.

## Source diff

The relevant changes in the copied classifier are:

```diff
- if crack_x_min < box_x_min - tol or front_x > box_x_max + tol: reject
+ if crack_x_min < box_x_min or front_x > box_x_max: reject
- if crack_z_min < box_z_min - tol or crack_z_max > box_z_max + tol: reject
+ if crack_z_min < box_z_min or crack_z_max > box_z_max: reject
- if not box.contains(source, tol): source is outside-domain UNKNOWN
+ if not box.contains(source, 0.0): source is outside-domain UNKNOWN
- reject an H section only when x/z overrun exceeds tolerance
+ reject every actual overrun; reject nonzero interior margins <= tolerance;
+ preserve exact equality at closed crack start/z endpoints
```

The finite front remains excluded from H by its existing tolerance band and positive-front-clearance predicate. The new endpoint check changes only the copied classifier; it does not alter the original 30/10/60 manifest or provider wiring.

## Focused receipt

The one focused run was:

```text
& 'C:\Users\xiemo\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' 'pgfrag/experiments/local_front_support_safe_v2_20261006/focused_geometry_check.py'
{"all_source_y_nonzero": true, "classification_queries": 106, "elapsed_limit_seconds": 10, "exact_box_clipped_endpoint": "HEAVISIDE", "finite_front_outside_box_within_tolerance": "REJECTED", "frozen_role_counts": {"BASE": 60, "FRONT": 30, "HEAVISIDE_ELIGIBLE": 10, "UNKNOWN": 0}, "frozen_scalar_dofs": 140, "interior_gap_above_tolerance": "HEAVISIDE", "nonzero_interior_gap_below_tolerance": "UNKNOWN", "overhang_below_and_above_tolerance": "UNKNOWN", "result": "PASS", "source_outside_closed_box": "UNKNOWN"}
```

The final invocation completed in about 0.4 seconds. Its 106 classifier calls comprise six explicit source classifications and the frozen 100-point lattice; the separate invalid-geometry construction check is rejected before classification. Covered cases are:

* The original sub-tolerance H overhang and a 1.5-tolerance overhang both return UNKNOWN.
* A nonzero 0.5-tolerance gap inside a closed crack end returns UNKNOWN; a 2-tolerance interior gap remains HEAVISIDE.
* The exact z endpoint clipped by Ω and included in Γ remains HEAVISIDE.
* The original front segment wholly above Ω by 0.1–0.9 tolerance is rejected as invalid geometry.
* A source center outside the closed box by 0.5 tolerance returns UNKNOWN.
* The 100 frozen sources remain 30 FRONT, 10 H, 60 BASE, 140 scalar DOFs, and all source y values remain nonzero.

There was one launch-only attempt before the successful run: the script initially walked one parent too far when resolving the repository root and stopped with `FileNotFoundError` before reading the case or classifying any geometry. The path was corrected; the receipt above is from the sole run that performed classifications.

## Limits

This check exercises only the copied pure-geometry classifier and the fixed source coordinates. It does not integrate the safe-v2 copy into the frozen provider, generate a role manifest, run MLS queries, assemble an operator, solve, or qualify the global collar. Exact endpoint acceptance relies on the shared floating-point representation of the box-clipped bound and crack endpoint; nonzero differences in the tolerance band fail closed as UNKNOWN. Broader fuzzing and cross-platform floating-point analysis were not performed.
