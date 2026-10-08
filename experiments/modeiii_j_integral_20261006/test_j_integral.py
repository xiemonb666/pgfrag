"""Small contract tests for the external-data J integrator."""

from __future__ import annotations

import unittest

import numpy as np

from j_integral import INPUT_UNITS, integrate_modeiii_j


def valid_inputs() -> dict[str, object]:
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
        "units": dict(INPUT_UNITS),
        "field_provenance": "contract-test placeholder field; not a solver result",
    }


class JIntegralContractTests(unittest.TestCase):
    def test_missing_data_is_unknown(self) -> None:
        result = integrate_modeiii_j()
        self.assertEqual(result.status, "UNKNOWN")
        self.assertIsNone(result.J_J_per_m2)
        self.assertEqual(result.reason, "missing_units")

    def test_missing_field_provenance_is_unknown(self) -> None:
        data = valid_inputs()
        data["field_provenance"] = None
        result = integrate_modeiii_j(**data)
        self.assertEqual(result.status, "UNKNOWN")
        self.assertEqual(result.reason, "missing_field_provenance")

    def test_invalid_weights_are_unknown(self) -> None:
        for invalid_weight in (0.0, -1.0):
            with self.subTest(invalid_weight=invalid_weight):
                data = valid_inputs()
                weights = np.asarray(data["weights_m"]).copy()
                weights[2] = invalid_weight
                data["weights_m"] = weights
                result = integrate_modeiii_j(**data)
                self.assertEqual(result.status, "UNKNOWN")
                self.assertEqual(result.reason, "line_weights_must_be_positive")

    def test_inward_normals_are_unknown(self) -> None:
        data = valid_inputs()
        data["normals"] = -np.asarray(data["normals"])
        result = integrate_modeiii_j(**data)
        self.assertEqual(result.status, "UNKNOWN")
        self.assertEqual(result.reason, "normals_must_point_outward_around_front")

    def test_nonfinite_samples_are_unknown(self) -> None:
        data = valid_inputs()
        coordinates = np.asarray(data["contour_coordinates_m"]).copy()
        coordinates[3, 1] = np.nan
        data["contour_coordinates_m"] = coordinates
        result = integrate_modeiii_j(**data)
        self.assertEqual(result.status, "UNKNOWN")
        self.assertEqual(result.reason, "nonfinite_contour_coordinates")

    def test_wrong_units_are_unknown(self) -> None:
        data = valid_inputs()
        data["units"] = {**INPUT_UNITS, "line_weight": "mm"}
        result = integrate_modeiii_j(**data)
        self.assertEqual(result.status, "UNKNOWN")
        self.assertEqual(result.reason, "invalid_or_incomplete_units")


if __name__ == "__main__":
    unittest.main()
