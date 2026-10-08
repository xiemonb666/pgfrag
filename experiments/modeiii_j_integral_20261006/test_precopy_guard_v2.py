"""Focused tests for v2's pre-allocation sample and row-shape guards."""

from __future__ import annotations

import json
from pathlib import Path
import time
import unittest
from unittest.mock import patch

import numpy as np

import j_integral_precopy_guard_v2 as guarded


RECEIPT_PATH = (
    Path(__file__).resolve().parents[3]
    / "pgfrag"
    / "tmp"
    / "review_20261006"
    / "modeiii_j_integral_20261006"
    / "precopy_guard_receipt.json"
)


def _valid_data() -> dict[str, object]:
    theta = np.arange(8, dtype=np.float64) * (2.0 * np.pi / 8.0)
    radius = 0.01
    return {
        "contour_coordinates_m": np.column_stack(
            (radius * np.cos(theta), radius * np.sin(theta), np.zeros_like(theta))
        ),
        "normals": np.column_stack((np.cos(theta), np.sin(theta), np.zeros_like(theta))),
        "weights_m": np.full(8, 2.0 * np.pi * radius / 8.0),
        "grad_u_z_m_per_m": np.zeros((8, 3), dtype=np.float64),
        "shear_modulus_Pa": 1.0e6,
        "units": dict(guarded.INPUT_UNITS),
        "field_provenance": "shape-guard-only placeholder; not a field validation",
    }


class PrecopyBudgetTests(unittest.TestCase):
    receipts: list[dict[str, object]] = []

    def test_257_rows_rejected_before_numpy_conversion(self) -> None:
        data = _valid_data()
        data["contour_coordinates_m"] = [(0.01, 0.0, 0.0)] * 257
        with patch.object(guarded.np, "asarray", side_effect=AssertionError("conversion occurred")):
            result = guarded.integrate_modeiii_j(**data)
        self.assertEqual(result.status, "UNKNOWN")
        self.assertEqual(result.reason, "sample_count_out_of_range")
        self.receipts.append(
            {
                "case": "257_rows_rejected_before_numpy_conversion",
                "status": "PASS",
                "reason": result.reason,
                "numpy_conversion_attempted": False,
            }
        )

    def test_long_row_shape_rejected_before_numpy_conversion(self) -> None:
        data = _valid_data()
        data["contour_coordinates_m"] = np.zeros((8, 4), dtype=np.float64)
        with patch.object(guarded.np, "asarray", side_effect=AssertionError("conversion occurred")):
            result = guarded.integrate_modeiii_j(**data)
        self.assertEqual(result.status, "UNKNOWN")
        self.assertEqual(result.reason, "invalid_contour_coordinates_shape")
        self.receipts.append(
            {
                "case": "wrong_trailing_shape_rejected_before_numpy_conversion",
                "status": "PASS",
                "reason": result.reason,
                "input_shape": [8, 4],
                "numpy_conversion_attempted": False,
            }
        )

    def test_small_correct_contour_shape_accepted(self) -> None:
        result = guarded.integrate_modeiii_j(**_valid_data())
        self.assertEqual(result.status, "OK")
        self.assertEqual(result.details["sample_count"], 8)
        self.receipts.append(
            {
                "case": "small_correct_contour_shape_accepted",
                "status": "PASS",
                "sample_count": result.details["sample_count"],
                "J_J_per_m2": result.J_J_per_m2,
            }
        )


def main() -> None:
    start = time.perf_counter()
    PrecopyBudgetTests.receipts = []
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(PrecopyBudgetTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    elapsed = time.perf_counter() - start
    receipt = {
        "artifact_id": "modeiii_j_integral_precopy_guard_v2",
        "status": "PASS" if result.wasSuccessful() and elapsed <= 5.0 else "FAIL",
        "elapsed_seconds": elapsed,
        "maximum_seconds": 5.0,
        "tests_run": result.testsRun,
        "failures": len(result.failures),
        "errors": len(result.errors),
        "cases": PrecopyBudgetTests.receipts,
        "field_reference_called": False,
        "analytic_224_sample_validation_rerun": False,
        "actual_static_j_computed": False,
    }
    RECEIPT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with RECEIPT_PATH.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(receipt, handle, indent=2, ensure_ascii=False, allow_nan=False)
        handle.write("\n")
    print(json.dumps(receipt, indent=2, ensure_ascii=False))
    if receipt["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
