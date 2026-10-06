# Convex-segment interval upper bound (v2 source proposal)

This source-only v2 proposal tightens the interval upper bound for the same supplied piecewise-linear 3-D polylines. It reuses the pinned v1 module's unit/input validation, rational point-to-segment squared distance, dyadic square-root enclosure, sample cap, and result type. It does not modify v1 or call a solver.

For one source interval with endpoints `a,b` and any fixed target segment `S`, the function `x -> distance(x,S)` is convex: distance to a convex set is convex, and a line segment is an affine parameterization. Thus for every `x` between `a` and `b`, `distance(x,S) <= max(distance(a,S), distance(b,S))`. The target polyline is a union of segments, so `distance(x,target) = min_S distance(x,S)` and its entire interval is bounded by

`Uconvex = min_S max(distance_upper(a,S), distance_upper(b,S))`.

Each endpoint distance is the v1 exact projection followed by the v1 dyadic square-root enclosure. The interval upper bound is `min(Ulipschitz, Uconvex)`, where `Ulipschitz` is the existing endpoint-plus-length bound. Taking the minimum preserves validity because both are certified upper bounds. Per sampled point, all target-segment distance enclosures are cached; a split computes the midpoint sample once and shares it between the new intervals. The two directed passes retain the common v1 cap of 128 samples.

The five planned v2 checks received one ROOT-approved focused run in `pgfrag/tmp/review_20261006/continuous_polyline_bias_bounds_20261006/focused_run02/`. All five passed. The captured records show an exact `[1/2,1/2] m` parallel bound and `[0,0] m` for the shared-segment metre/millimetre case, both from four samples; the same-vertex-set fixture brackets the analytic `H^2=1/2 m^2`; the non-dyadic tiny-tolerance case reaches 128 samples and correctly returns `UNKNOWN_WITH_BOUNDS`; proper poses and reversal return the exact `[1/4,1/4] m` bound. The run receipt and all exact rational bounds, sample counts, and statuses are preserved in that folder.

The certificate covers only the supplied polylines in a shared frame. It does not establish a bound for an underlying crack surface, solver lineage, physical admission, or absence of tetrahedral orientation bias. Results remain `GIVEN_PIECEWISE_LINEAR_POLYLINES_ONLY` and `NOT_ADMITTED`.
