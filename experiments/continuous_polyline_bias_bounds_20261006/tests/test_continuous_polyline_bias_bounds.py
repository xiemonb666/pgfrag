"""Planned checks only; this file was not executed in the SOURCE ONLY pass."""

from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from continuous_polyline_bias_bounds import (  # noqa: E402
    InputError,
    symmetric_polyline_hausdorff_bounds,
)


class ContinuousPolylineBoundsTests(unittest.TestCase):
    def test_parallel_segments_bracket_known_offset(self):
        result = symmetric_polyline_hausdorff_bounds(
            [[0, 0, 0], [1, 0, 0]], [[0, Fraction(1, 2), 0], [1, Fraction(1, 2), 0]],
            units_a="m", units_b="m", atol_m=Fraction(1, 100),
        )
        self.assertLessEqual(result.lower_m, Fraction(1, 2))
        self.assertGreaterEqual(result.upper_m, Fraction(1, 2))

    def test_crossing_polyline_repartition_has_zero_lower_bound(self):
        result = symmetric_polyline_hausdorff_bounds(
            [[-1, 0, 0], [0, 0, 0], [1, 0, 0]],
            [[0, 0, 0], [-1, 0, 0], [1, 0, 0]],
            units_a="m", units_b="m", atol_m=2,
        )
        self.assertEqual(result.lower_m, 0)

    def test_three_dimensional_offset(self):
        result = symmetric_polyline_hausdorff_bounds(
            [[0, 0, 0], [1, 0, 0]], [[0, 0, Fraction(1, 2)], [1, 0, Fraction(1, 2)]],
            units_a="m", units_b="m", atol_m=Fraction(1, 100),
        )
        self.assertLessEqual(result.lower_m, Fraction(1, 2))
        self.assertGreaterEqual(result.upper_m, Fraction(1, 2))

    def test_millimetres_convert_to_the_shared_metre_frame(self):
        result = symmetric_polyline_hausdorff_bounds(
            [[0, 0, 0], [1, 0, 0]], [[0, 0, 0], [1000, 0, 0]],
            units_a="m", units_b="mm", atol_m=2,
        )
        self.assertEqual(result.lower_m, 0)
        self.assertEqual(result.upper_m, Fraction(1, 2))

    def test_same_vertex_set_fixture_brackets_analytic_hausdorff(self):
        case_path = Path(__file__).resolve().parents[1] / "fixtures" / "sparse_vertex_set_counterexample.json"
        case = json.loads(case_path.read_text(encoding="utf-8"))
        self.assertEqual(set(map(tuple, case["polyline_a"])), set(map(tuple, case["polyline_b"])))
        result = symmetric_polyline_hausdorff_bounds(
            case["polyline_a"], case["polyline_b"],
            units_a=case["coordinate_units"], units_b=case["coordinate_units"],
            atol_m=Fraction(1, 100),
        )
        self.assertGreater(result.lower_m, 0)
        self.assertLessEqual(result.lower_m ** 2, Fraction(1, 2))
        self.assertLessEqual(Fraction(1, 2), result.upper_m ** 2)

    def test_rejects_zero_length_segment_bool_nonfinite_and_unknown_units(self):
        common = {"units_a": "m", "units_b": "m", "atol_m": 1}
        with self.assertRaises(InputError):
            symmetric_polyline_hausdorff_bounds([[0, 0, 0], [0, 0, 0]], [[0, 0, 0], [1, 0, 0]], **common)
        with self.assertRaises(InputError):
            symmetric_polyline_hausdorff_bounds([[False, 0, 0], [1, 0, 0]], [[0, 0, 0], [1, 0, 0]], **common)
        with self.assertRaises(InputError):
            symmetric_polyline_hausdorff_bounds([[0, 0, 0], [float("inf"), 0, 0]], [[0, 0, 0], [1, 0, 0]], **common)
        with self.assertRaises(InputError):
            symmetric_polyline_hausdorff_bounds([[0, 0, 0], [1, 0, 0]], [[0, 0, 0], [1, 0, 0]], units_a="inch", units_b="m", atol_m=1)

    def test_cap_returns_unknown_with_certified_bounds(self):
        result = symmetric_polyline_hausdorff_bounds(
            [[0, 0, 0], [1, 0, 0]], [[0, Fraction(1, 3), 0], [1, Fraction(1, 3), 0]],
            units_a="m", units_b="m", atol_m=Fraction(1, 10**30),
        )
        self.assertEqual(result.status, "UNKNOWN_WITH_BOUNDS")
        self.assertLessEqual(result.samples, 128)
        self.assertLessEqual(result.lower_m, Fraction(1, 3))
        self.assertGreaterEqual(result.upper_m, Fraction(1, 3))

    def test_proper_rigid_pose_and_reversal_preserve_bounds(self):
        a = [[0, 0, 0], [1, 0, 0], [1, 1, 0]]
        b = [[0, 0, Fraction(1, 4)], [1, 0, Fraction(1, 4)], [1, 1, Fraction(1, 4)]]
        moved_a = [[5, -2, 7], [5, -1, 7], [4, -1, 7]]
        moved_b = [[5, -2, Fraction(29, 4)], [5, -1, Fraction(29, 4)], [4, -1, Fraction(29, 4)]]
        original = symmetric_polyline_hausdorff_bounds(a, b, units_a="m", units_b="m", atol_m=Fraction(1, 100))
        posed = symmetric_polyline_hausdorff_bounds(moved_a, moved_b, units_a="m", units_b="m", atol_m=Fraction(1, 100))
        reversed_order = symmetric_polyline_hausdorff_bounds(a[::-1], b[::-1], units_a="m", units_b="m", atol_m=Fraction(1, 100))
        self.assertEqual((original.lower_m, original.upper_m), (posed.lower_m, posed.upper_m))
        self.assertEqual((original.lower_m, original.upper_m), (reversed_order.lower_m, reversed_order.upper_m))

if __name__ == "__main__":
    unittest.main()
