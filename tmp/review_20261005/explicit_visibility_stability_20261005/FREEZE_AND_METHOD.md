# Source freeze and exact-certificate method

Status: source implementation and test source prepared for algebra review. No pytest, solver, renderer, or runtime check has been run in this review phase. Existing `explicit_visibility_20261005/visibility3d.py` and its 17 source tests remain read-only.

## Scope and contract

- New source ownership: `pgfrag/experiments/explicit_visibility_stability_20261005/`.
- Review ownership: `pgfrag/tmp/review_20261005/explicit_visibility_stability_20261005/`.
- Geometry only. No nucleation/growth rule, crack identity, tetra adjacency/direction cut selection, material/ceramic physics, solver integration, or render/real-fracture claim.
- Mechanics coupling remains `NOT_COUPLED`.
- Coordinates are SI metres. The fixed visibility reference is `visibility3d.py`, SHA-256 `9d5097ca33ce899cf788ad43ef0740a34296216a8fd7aa490868c27ffe1ebb03`; it is dynamically loaded only after checking that exact hash.
- Input arrays must be exact `numpy.ndarray` objects, float64, C-contiguous, finite, and bounded before private copying: source/query row caps 500 each, triangle cap 256, 16,384 query/source pairs, and 16,384 source/query/facet work units. Above either pair/work cap the result is immutable UNKNOWN, with no Fraction conversion, float-classifier call, or per-facet result allocation.
- Every input float is interpreted as its exact dyadic rational. Input numerator/denominator bit length is at most 256. Before each Fraction add/subtract/multiply/divide, a conservative upper bound for the resulting numerator/denominator bit lengths is checked against 2,048; violations return per-facet UNKNOWN. The arithmetic budget counts bounded exact predicates, not physical length tolerances.
- Returned pair/facet arrays are copied into immutable bytes-backed storage. `stable_visible_mask` is true only after strict exact classification of every facet and agreement with the pinned float result. There is no code path from UNKNOWN to visible.

## Predicate derivation

For facet vertices `A,B,C`, define the exact normal `N=(B-A)×(C-A)`. `N=0` is degenerate and UNKNOWN. For fixed source `S` and query `Q`, use signed plane numerators `dS=N·(S-A)` and `dQ=N·(Q-A)`; these retain their signs without normalizing the plane or comparing a squared area to a length tolerance.

1. If both signs are strict and equal, the segment cannot cross the facet plane. This relation persists for all sufficiently small query perturbations, so the finite-facet result is stable visible.
2. If both are strict and opposite, `t=dS/(dS-dQ)` and `P=S+t(Q-S)` are exact. Compute barycentric coordinates of `P` in `(A,B,C)` by the exact Gram system. All three strictly positive means a transverse relative-interior crossing and remains so in a sufficiently small open query neighbourhood: stable blocked. Any strictly negative coordinate means the plane intersection lies strictly outside the finite closed triangle and remains so nearby: stable visible. No negative and at least one zero is an edge/front contact: UNKNOWN.
3. If exactly one endpoint is on the plane, test that endpoint's exact barycentric coordinates against the finite triangle. Any strict negative coordinate proves it is separated from the finite facet and supports a local visible certificate. An inside/front/edge endpoint is UNKNOWN. If both are on the plane, the segment is coplanar and UNKNOWN.

For multiple facets, aggregate only if every facet is strict. A strict blocker plus any ambiguous facet therefore remains aggregate UNKNOWN; both the blocked facet and ambiguous facet are preserved in the facet arrays. This conservative policy prevents a certificate from depending on the visibility primitive's different blocker-dominates-ambiguity convention.

The candidate pair class is then compared with the frozen float classifier evaluated on a private snapshot of the same float64 arrays. Any disagreement is UNKNOWN. This check is an implementation-consistency gate, not part of the exact continuity proof.

## Certificate limits

The proof is an ideal-real existence proof around exact dyadic input coordinates, with source and triangles fixed and query varying in 3D. No radius is computed. It is not an IEEE arbitrary-perturbation/roundoff guarantee. A rounded rotation or unit conversion changes the dyadic inputs; tests can compare decisive classifications away from boundaries but cannot claim exact predicate invariance between pre- and post-transform coordinates. The frozen float classifier adds configured SI tolerance and a 32-epsilon absolute-world-coordinate precision guard; exact certificates are withheld when its status differs. This does not establish physical correctness, crack-surface ground truth, or mechanics coupling.

## Source-only test inventory

`test_stability3d.py` contains hand cases for: strict interior crossing and strict plane miss; same-side visibility; source/query on-plane finite-facet misses; source/query front points, edge crossing, coplanarity, and degeneracy; definite blocker plus ambiguous facet under all-facets-strict aggregation; exact-vs-float tolerance mismatch; no-facet visible identity; a proper non-axis-aligned Rodrigues quarter-turn on a `1e-3 m` fixture and reversed winding; input-ratio, intermediate arithmetic, and work/pair budgets; bytes-backed immutability; array-subclass conversion hooks and shape preflight; nonfinite and non-float64 input rejection. The test file is not authorized to run before this algebra review.

## Freeze boundary

The source/test/README contents listed in `SHA256SUMS.txt` define the review candidate. No implementation edits should be made after algebra review begins; any requested change requires a new freeze and new hashes before the single focused test run.
