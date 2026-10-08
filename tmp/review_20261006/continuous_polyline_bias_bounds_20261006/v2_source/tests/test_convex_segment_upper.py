"""Planned v2 checks; source-only pass does not execute these tests."""

import json
from fractions import Fraction
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from convex_segment_upper import symmetric_polyline_hausdorff_bounds_convex  # noqa: E402


class ConvexSegmentUpperTests(unittest.TestCase):
    def test_parallel_offset_is_exact_without_adaptive_splits(self):
        result = symmetric_polyline_hausdorff_bounds_convex(
            [[0, 0, 0], [1, 0, 0]], [[0, 0, Fraction(1, 2)], [1, 0, Fraction(1, 2)]],
            units_a="m", units_b="m", atol_m=Fraction(1, 100),
        )
        self.assertEqual((result.lower_m, result.upper_m), (Fraction(1, 2), Fraction(1, 2)))
        self.assertEqual(result.samples, 4)

    def test_same_segment_and_millimetre_conversion_have_exact_zero(self):
        result = symmetric_polyline_hausdorff_bounds_convex(
            [[0, 0, 0], [1, 0, 0]], [[0, 0, 0], [1000, 0, 0]],
            units_a="m", units_b="mm", atol_m=Fraction(1, 100),
        )
        self.assertEqual((result.lower_m, result.upper_m), (Fraction(0), Fraction(0)))

    def test_fixture_brackets_analytic_same_vertex_set_distance(self):
        case_path = Path(__file__).resolve().parents[1] / "fixtures" / "sparse_vertex_set_counterexample.json"
        case = json.loads(case_path.read_text(encoding="utf-8"))
        self.assertEqual(set(map(tuple, case["polyline_a"])), set(map(tuple, case["polyline_b"])))
        result = symmetric_polyline_hausdorff_bounds_convex(
            case["polyline_a"], case["polyline_b"],
            units_a=case["coordinate_units"], units_b=case["coordinate_units"],
            atol_m=Fraction(1, 100),
        )
        self.assertGreater(result.lower_m, 0)
        self.assertLessEqual(result.lower_m ** 2, Fraction(1, 2))
        self.assertLessEqual(Fraction(1, 2), result.upper_m ** 2)

    def test_non_dyadic_offset_and_tiny_tolerance_still_hit_cap(self):
        result = symmetric_polyline_hausdorff_bounds_convex(
            [[0, 0, 0], [1, 0, 0]], [[0, Fraction(1, 3), 0], [1, Fraction(1, 3), 0]],
            units_a="m", units_b="m", atol_m=Fraction(1, 10**30),
        )
        self.assertEqual(result.status, "UNKNOWN_WITH_BOUNDS")
        self.assertEqual(result.samples, 128)

    def test_proper_pose_and_reversed_order_preserve_tight_bounds(self):
        a = [[0, 0, 0], [1, 0, 0], [1, 1, 0]]
        b = [[0, 0, Fraction(1, 4)], [1, 0, Fraction(1, 4)], [1, 1, Fraction(1, 4)]]
        moved_a = [[5, -2, 7], [5, -1, 7], [4, -1, 7]]
        moved_b = [[5, -2, Fraction(29, 4)], [5, -1, Fraction(29, 4)], [4, -1, Fraction(29, 4)]]
        original = symmetric_polyline_hausdorff_bounds_convex(a, b, units_a="m", units_b="m", atol_m=Fraction(1, 100))
        posed = symmetric_polyline_hausdorff_bounds_convex(moved_a, moved_b, units_a="m", units_b="m", atol_m=Fraction(1, 100))
        reversed_order = symmetric_polyline_hausdorff_bounds_convex(a[::-1], b[::-1], units_a="m", units_b="m", atol_m=Fraction(1, 100))
        self.assertEqual((original.lower_m, original.upper_m), (Fraction(1, 4), Fraction(1, 4)))
        self.assertEqual((original.lower_m, original.upper_m), (posed.lower_m, posed.upper_m))
        self.assertEqual((original.lower_m, original.upper_m), (reversed_order.lower_m, reversed_order.upper_m))


if __name__ == "__main__":
    unittest.main()
