# Continuous polyline Hausdorff bounds

This experiment estimates the symmetric Hausdorff distance between two **given, ordered, piecewise-linear 3-D curves**. It addresses a measurement gap in vertex-only crack-path comparisons: the vertices can agree as an unordered point set while the line segments between them differ. The implementation does not call a solver, infer missing geometry, or claim to measure a physical crack surface.

## Contract

Call `symmetric_polyline_hausdorff_bounds(polyline_a, polyline_b, units_a=..., units_b=..., atol_m=...)` from `continuous_polyline_bias_bounds.py`. Each curve must contain 2–32 3-D vertices. Adjacent duplicates, empty curves, booleans, nonfinite values, unknown units, coordinates above 1000 m in magnitude, and scalar exact fractions exceeding 256 numerator or denominator bits are rejected. Accepted length units are `m`, `mm`, and `um`; conversion to metres uses rational factors. Floats are interpreted through their decimal spelling, then all geometry operations use exact `Fraction` arithmetic.

The tolerance is a positive value in metres. The algorithm starts with every input vertex in both directed comparisons, then bisects the interval with the largest certified upper bound. It stops when the symmetric upper-minus-lower gap is no greater than the requested tolerance. There are at most 128 point-distance samples in total. If that cap is reached first, the result is `UNKNOWN_WITH_BOUNDS`; callers must retain both bounds and must not convert that status into a pass.

`HausdorffBounds.as_record()` exposes metre bounds as exact `numerator/denominator` strings. Every result carries `geometry_scope=GIVEN_PIECEWISE_LINEAR_POLYLINES_ONLY` and `physical_admission=NOT_ADMITTED`. A physical claim still requires real solver lineage, mesh and extraction provenance, and review outside this metric.

## Bound argument

For a query point and target segment with direction `v`, the exact projection parameter is `(w·v)/(v·v)`, clamped to `[0,1]`; squared distance to the projected point is rational. Taking the minimum over target segments gives the exact squared distance to the target polyline. Its square root is enclosed on a dyadic grid: for `N=40`, integer square root of `floor(squared_distance * 2^(2N))` gives a lower endpoint, and the upper endpoint is one grid unit higher unless equality is exact. Thus each sampled point distance has a rigorous enclosure no wider than `2^-40 m`.

Distance to a fixed set is 1-Lipschitz. On a source interval of length `L`, with endpoint distances `g0` and `g1`, every interior distance is bounded above by `min(g0+x, g1+L-x)`, whose maximum is `(g0+g1+L)/2`. Substituting endpoint-distance upper bounds and a rigorously rounded-up interval length therefore bounds every point on that whole source interval. The maximum lower enclosure over sampled points bounds the directed maximum from below; the maximum interval bound over both directed traversals bounds the symmetric Hausdorff distance from above. Rational comparisons make the stopping decision exact.

These are bounds for the supplied polygonal lines in a shared coordinate frame. They are not bounds for a continuous or smooth crack surface, and they cannot establish that mesh-orientation bias is absent.

## Planned checks (source only)

`tests/test_continuous_polyline_bias_bounds.py` records planned unit and invariance checks. The source is included for review; this SOURCE ONLY pass did not execute it or evaluate the metric. `fixtures/sparse_vertex_set_counterexample.json` is a qualitative static case showing that unordered vertex-set Hausdorff can be zero while the two continuous polylines differ.
