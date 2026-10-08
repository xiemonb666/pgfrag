# Independent finite-support topology review

Status: `FIXED_CASE_GEOMETRY_OK / GENERIC_TOLERANCE_EDGE_FINDINGS`. This review concerns only `support_topology.py`, `fixed_case.json`, `role_manifest.json`, and the existing `REVIEW.md`. It makes no claim about an MLS provider, operator, solver, or the unproved global collar.

## Fixed-case result

The frozen case is geometrically consistent with the stated finite-front contract. An independent 100-point calculation over the listed source lattice gives 30 FRONT, 10 HEAVISIDE_ELIGIBLE, and 60 BASE roles, with no UNKNOWN points. The manifest reports the same counts and 140 scalar DOFs. All source y coordinates are nonzero.

For the ten H candidates, the plane section is clipped by the box before its bounds are compared with the finite crack. At z=±8 mm the clipped section ends exactly at the corresponding box/crack endpoint; its other z bound is about ±2.255 mm. Its x upper bound is about −2.255 mm, the two-side margin is 3 mm, and the finite-front clearance is about 1.944 mm. The FRONT group has minimum front-intersection margin about 1.343 mm. These margins are many orders larger than the 1e−12 m tolerance. The finite-ahead source at (+8,+4,0) mm is BASE because its nearest finite crack/front distance is √80 mm² ≈ 8.944 mm, while its support radius is 7 mm. Thus the inspected fixed lattice is not affected by the tolerance-edge findings below.

## Findings

### [P2] H containment accepts a plane section that protrudes beyond a finite crack end

At `support_topology.py:304-307`, the lower x and z checks and the upper z check allow up to `tolerance_m` of actual protrusion. The certificate later sets `plane_section_inside_finite_crack` from the H role alone (`support_topology.py:215-216`), so the evidence says “inside” even when the bounds are outside. This conflicts with the fixed-case policy that tolerance-band cases are UNKNOWN.

Concrete geometry: keep the fixed box, h=7 mm, and tolerance 1e−12 m, but set crack z bounds to [−8 mm, 8 mm−0.5e−12 m]. For source (−8 mm, +4 mm, +4 mm), the plane-section radius is √33 mm² ≈ 5.744563 mm. Its clipped section has x upper bound −2.255437 mm and z upper bound +8 mm. The source-to-front distance is √80 mm² ≈ 8.944272 mm, its crack distance is 4 mm, and its two-side margin is 3 mm. Every other H condition passes; the z overrun is only 0.5 tolerance, so the code certifies H and records two components. In exact geometry, however, the section contains a positive-width strip beyond Γ. The support can pass through that strip, so removing Γ does not split this support into two components.

### [P2] Crack/front geometry just outside the box is accepted and measured un-clipped

`FiniteCrackGeometry.__post_init__` (`support_topology.py:76-79`) permits crack z bounds up to one tolerance outside the box. Distances in `_distance_to_front` and `_distance_to_finite_crack` are then measured to that un-clipped geometry (`support_topology.py:125-135`), even though the actual support is Ω intersected with the ball.

Concrete geometry: retain the fixed box, h=7 mm, and tolerance 1e−12 m; set the crack z interval to [8 mm+0.1e−12 m, 8 mm+0.9e−12 m]. The constructor accepts it because both endpoints are within one tolerance of the box. For source (+4 mm,+4 mm,+8 mm), the un-clipped front distance is about 5.657 mm, so classification reaches FRONT and gives a unique shifted F1 source value. But the finite front and crack both lie entirely above Ω, and therefore Ω intersect the support ball meets neither feature; the exact clipped-support role is BASE. This is a tolerance-edge false positive in the generic geometry API, not in the frozen case whose endpoints lie exactly on the box.

## Boundary checks and confidence

The fixed-case front tangent `(h,0,0)` and crack tangent `(−8 mm,h,0)` fall into the respective UNKNOWN branches. A finite-ahead support that overlaps the infinite plane but stays wholly ahead of the finite tip is not promoted to H. Box-clipped H endpoints remain valid for the fixed geometry because Γ reaches the same z bounds as Ω. These checks support the finite-geometry design and the fixed counts.

The two findings have high confidence as statements about the current generic classifier: both follow directly from its accepted geometry predicates and branch conditions, and the examples have strict margins from the unrelated front/crack thresholds. They are low practical risk for this frozen lattice because its geometric margins are millimetres while tolerance is 1e−12 m. The findings matter if this API is reused for another geometry or claims fail-closed topology certification at tolerance edges.

## Suggested correction

Either require crack/front extents to lie inside Ω exactly, or normalize the geometry to Ω before any distance or containment calculation and bind that normalized geometry into the case. For H, certify only exact section containment; if an endpoint mismatch is within tolerance, return UNKNOWN unless a documented snap rule is applied consistently before classification. Derive the `plane_section_inside_finite_crack` certificate from the actual bounds check rather than from the assigned role.

## Scope note

The requested prior review was found at `pgfrag/tmp/review_20261006/local_front_support_contract_20261006/REVIEW.md`; the top-level `tmp/review_20261006/...` path is absent in this checkout. No author files were changed. The independent lattice enumeration and two constructed tolerance-edge examples used 102 geometry evaluations total, under the stated 256-query and 30-second limits; no MLS or solver code was invoked.
