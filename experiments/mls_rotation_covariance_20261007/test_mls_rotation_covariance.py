"""Future focused numeric check; do not run until ROOT authorizes it."""
import json
import sys
import unittest
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
REVIEW = ROOT / "pgfrag/tmp/review_20261007/mls_rotation_covariance_20261007/focused_run01"
sys.path.insert(0, str(HERE))
from mls_rotation_covariance import evaluate  # noqa: E402


class RotationCovarianceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = evaluate()
        REVIEW.mkdir(parents=True, exist_ok=True)
        target = REVIEW / "result.json"
        if target.exists():
            raise FileExistsError(target)
        target.write_text(json.dumps(cls.data, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        cls.rows = cls.data["results"]

    def test_fixture_is_proper_non_axis_rotation_and_budgeted(self):
        r = np.asarray(self.data["rotation"])
        self.assertLess(np.max(np.abs(r.T @ r - np.eye(3))), 1e-14)
        self.assertAlmostEqual(float(np.linalg.det(r)), 1.0, places=14)
        self.assertGreater(np.count_nonzero(np.abs(r) > 1e-12), 3)
        self.assertEqual(self.data["new_provider_queries"], 64)

    def test_all_real_provider_rows_are_accepted_at_frozen_threshold(self):
        for row in self.rows.values():
            self.assertEqual(row["statuses"], ["OK"] * 16)
            self.assertEqual(row["support_counts"], [27] * 16)
            self.assertTrue(np.isfinite(row["weights"]).all())
            self.assertTrue(np.isfinite(row["gradients_per_m"]).all())
            self.assertTrue(all(np.isfinite(row["condition_numbers"])))
            self.assertLessEqual(max(row["condition_numbers"]), 1.0e12)

    def test_rigid_rotation_preserves_weights_and_rotates_physical_gradients(self):
        c = self.data["comparisons"]
        self.assertLessEqual(c["rotation_weights_max_abs"], 1.0e-10)
        self.assertLessEqual(c["rotation_gradient_normalized_max_abs"], 1.0e-10)

    def test_uniform_mm_conversion_preserves_weights_and_metre_gradients(self):
        c = self.data["comparisons"]
        self.assertLessEqual(c["millimetre_weights_max_abs"], 1.0e-10)
        self.assertLessEqual(c["millimetre_gradient_normalized_max_abs"], 1.0e-10)

    def test_bad_anisotropic_control_is_detected_by_actual_shape_weights(self):
        self.assertGreater(self.data["comparisons"]["bad_anisotropic_weights_max_abs"], 1.0e-8)


if __name__ == "__main__":
    unittest.main()
