"""Independent analytic checks for the source-frozen Mode III field.

This file is prepared for a single focused test authorization; it has not been
executed. Central differences stay away from the front and branch face.
"""

from __future__ import annotations

from dataclasses import FrozenInstanceError
from pathlib import Path
import sys
import unittest

import numpy as np


sys.path.insert(0, str(Path(__file__).resolve().parent))
import field_reference as reference  # noqa: E402


class ModeIIIFieldReferenceTests(unittest.TestCase):
    def setUp(self) -> None:
        self.case = reference.DEFAULT_CASE

    def test_fixed_case_is_immutable_and_has_the_frozen_si_values(self) -> None:
        self.assertEqual(self.case.radius_m, 1.0e-2)
        self.assertEqual(self.case.length_m, 1.0e-2)
        self.assertEqual(self.case.K_III_Pa_sqrt_m, 100.0)
        self.assertEqual(self.case.young_modulus_Pa, 1.0e6)
        self.assertEqual(self.case.poisson_ratio, 0.3)
        self.assertAlmostEqual(self.case.shear_modulus_Pa, 1.0e6 / 2.6)
        with self.assertRaises(FrozenInstanceError):
            self.case.radius_m = 2.0e-2

    def test_analytic_gradient_matches_central_difference_away_from_cut(self) -> None:
        point = np.array([[2.0e-3, 3.0e-3, -1.0e-3]], dtype=np.float64)
        sample = reference.evaluate_field(point)
        self.assertIs(sample.status_at(0), reference.FieldStatus.OK)
        step = 1.0e-7
        numerical = np.zeros(3, dtype=np.float64)
        for axis in range(3):
            plus = point.copy()
            minus = point.copy()
            plus[0, axis] += step
            minus[0, axis] -= step
            f_plus = reference.evaluate_field(plus).f1_m_sqrt[0]
            f_minus = reference.evaluate_field(minus).f1_m_sqrt[0]
            numerical[axis] = (f_plus - f_minus) / (2.0 * step)
        np.testing.assert_allclose(
            sample.grad_f1_per_sqrt_m[0], numerical, rtol=2.0e-7, atol=2.0e-10
        )
        np.testing.assert_allclose(numerical[2], 0.0, atol=2.0e-10)
        np.testing.assert_allclose(
            sample.displacement_gradient[0, 2],
            self.case.amplitude_sqrt_m * sample.grad_f1_per_sqrt_m[0],
            rtol=0.0,
            atol=0.0,
        )
        mu = self.case.shear_modulus_Pa
        np.testing.assert_allclose(
            sample.stress_Pa[0, 0, 2],
            mu * sample.displacement_gradient[0, 2, 0],
            rtol=2.0e-15,
        )
        np.testing.assert_allclose(
            sample.stress_Pa[0, 1, 2],
            mu * sample.displacement_gradient[0, 2, 1],
            rtol=2.0e-15,
        )

    def test_front_and_exact_branch_face_are_unknown_without_regularization(self) -> None:
        points = np.array(
            [[0.0, 0.0, 0.0], [-4.0e-3, 0.0, 0.0], [-4.0e-3, -0.0, 0.0]],
            dtype=np.float64,
        )
        result = reference.evaluate_field(points)
        self.assertEqual(
            result.statuses.tolist(),
            [
                reference.FieldStatus.UNKNOWN_FRONT.value,
                reference.FieldStatus.UNKNOWN_BRANCH_FACE.value,
                reference.FieldStatus.UNKNOWN_BRANCH_FACE.value,
            ],
        )

    def test_crack_face_traces_have_the_sqrt_r_mode_iii_jump(self) -> None:
        radius = 4.0e-3
        epsilons = (4.0e-5, 4.0e-6, 4.0e-7)
        relative_errors = []
        for epsilon in epsilons:
            points = np.array(
                [[-radius, epsilon, 0.0], [-radius, -epsilon, 0.0]],
                dtype=np.float64,
            )
            result = reference.evaluate_field(points)
            self.assertTrue(np.all(result.statuses == reference.FieldStatus.OK.value))
            jump = result.displacement_m[0, 2] - result.displacement_m[1, 2]
            exact_jump = 2.0 * self.case.amplitude_sqrt_m * np.sqrt(radius)
            relative_errors.append(abs(jump - exact_jump) / abs(exact_jump))
        self.assertGreater(relative_errors[0], relative_errors[1])
        self.assertGreater(relative_errors[1], relative_errors[2])
        self.assertLess(relative_errors[2], 1.0e-8)

    def test_positive_x_trace_and_off_crack_shadow_path_are_continuous(self) -> None:
        radius = 4.0e-3
        epsilons = (4.0e-5, 4.0e-6, 4.0e-7)
        positive_x_jumps = []
        for epsilon in epsilons:
            points = np.array(
                [[radius, epsilon, 0.0], [radius, -epsilon, 0.0]],
                dtype=np.float64,
            )
            result = reference.evaluate_field(points)
            self.assertTrue(np.all(result.statuses == reference.FieldStatus.OK.value))
            positive_x_jumps.append(
                abs(result.f1_m_sqrt[0] - result.f1_m_sqrt[1])
            )
        self.assertGreater(positive_x_jumps[0], positive_x_jumps[1])
        self.assertGreater(positive_x_jumps[1], positive_x_jumps[2])

        # x=0, y>0 is a smooth path away from both the real cut and front.
        epsilon = 1.0e-7
        points = np.array(
            [[-epsilon, 3.0e-3, 0.0], [0.0, 3.0e-3, 0.0], [epsilon, 3.0e-3, 0.0]],
            dtype=np.float64,
        )
        result = reference.evaluate_field(points)
        self.assertTrue(np.all(result.statuses == reference.FieldStatus.OK.value))
        self.assertLess(
            abs(result.f1_m_sqrt[0] - result.f1_m_sqrt[2]),
            2.0e-4,
        )
        self.assertLess(
            abs(result.f1_m_sqrt[1] - 0.5 * (result.f1_m_sqrt[0] + result.f1_m_sqrt[2])),
            2.0e-8,
        )

        # A negative-x branch plane beyond either end cap is outside the
        # manufactured domain, so it cannot create an admissible ghost cut.
        beyond_cap = np.array([[ -3.0e-3, 0.0, 0.5 * self.case.length_m + 1.0e-6]])
        outside = reference.evaluate_field(beyond_cap)
        self.assertIs(outside.status_at(0), reference.FieldStatus.OUTSIDE_DOMAIN)

    def test_crack_face_traction_limit_vanishes_on_both_sides(self) -> None:
        radius = 4.0e-3
        magnitudes = []
        for epsilon in (4.0e-5, 4.0e-6, 4.0e-7):
            points = np.array(
                [[-radius, epsilon, 0.0], [-radius, -epsilon, 0.0]],
                dtype=np.float64,
            )
            normals = np.array([[0.0, 1.0, 0.0], [0.0, -1.0, 0.0]], dtype=np.float64)
            result = reference.evaluate_field(points, normals=normals)
            self.assertTrue(np.all(result.statuses == reference.FieldStatus.OK.value))
            magnitudes.append(float(np.max(np.linalg.norm(result.traction_Pa, axis=1))))
        self.assertGreater(magnitudes[0], magnitudes[1])
        self.assertGreater(magnitudes[1], magnitudes[2])

    def test_exact_traction_has_side_wall_and_both_end_cap_contributions(self) -> None:
        radius = self.case.radius_m
        length = self.case.length_m
        diagonal = radius / (2.0 * np.sqrt(2.0))
        cap_points = np.array(
            [[diagonal, diagonal, -0.5 * length], [diagonal, diagonal, 0.5 * length]],
            dtype=np.float64,
        )
        cap_normals = np.array([[0.0, 0.0, -1.0], [0.0, 0.0, 1.0]])
        caps = reference.evaluate_field(cap_points, normals=cap_normals)
        self.assertTrue(np.all(caps.statuses == reference.FieldStatus.OK.value))
        self.assertGreater(np.linalg.norm(caps.traction_Pa[0]), 0.0)
        self.assertGreater(np.linalg.norm(caps.traction_Pa[1]), 0.0)
        self.assertGreater(abs(caps.traction_Pa[1, 0]), 0.0)
        self.assertGreater(abs(caps.traction_Pa[1, 1]), 0.0)
        np.testing.assert_allclose(caps.traction_Pa[0], -caps.traction_Pa[1], rtol=2.0e-15)

        side_point = np.array([[0.0, radius, 0.0]])
        side_normal = np.array([[0.0, 1.0, 0.0]])
        side = reference.evaluate_field(side_point, normals=side_normal)
        self.assertIs(side.status_at(0), reference.FieldStatus.OK)
        self.assertGreater(abs(side.traction_Pa[0, 2]), 0.0)

    def test_polar_quadrature_matches_finite_total_strain_energy(self) -> None:
        gauss, weights = np.polynomial.legendre.leggauss(4)
        radial = 0.5 * self.case.radius_m * (gauss + 1.0)
        radial_weights = 0.5 * self.case.radius_m * weights
        theta = math_pi_times(gauss + 1.0)
        theta_weights = math_pi_times(weights)
        axial = 0.5 * self.case.length_m * gauss
        axial_weights = 0.5 * self.case.length_m * weights
        rr, tt, zz = np.meshgrid(radial, theta, axial, indexing="ij")
        wr, wt, wz = np.meshgrid(
            radial_weights, theta_weights, axial_weights, indexing="ij"
        )
        points = np.column_stack(
            (rr.ravel() * np.cos(tt.ravel()), rr.ravel() * np.sin(tt.ravel()), zz.ravel())
        )
        result = reference.evaluate_field(points)
        self.assertEqual(points.shape[0], 64)
        self.assertTrue(np.all(result.statuses == reference.FieldStatus.OK.value))
        polar_integral = float(
            np.sum(
                result.energy_density_J_per_m3
                * rr.ravel()
                * wr.ravel()
                * wt.ravel()
                * wz.ravel()
            )
        )
        exact = self.case.exact_total_strain_energy_J
        self.assertGreater(exact, 0.0)
        self.assertLessEqual(abs(polar_integral - exact) / exact, 2.0e-13)

    def test_results_are_plain_immutable_arrays_and_query_budget_is_enforced(self) -> None:
        point = np.array([[1.0e-3, 2.0e-3, 0.0]], dtype=np.float64)
        result = reference.evaluate_field(point)
        self.assertIs(type(result.stress_Pa), np.ndarray)
        self.assertFalse(result.stress_Pa.flags.writeable)
        with self.assertRaises(ValueError):
            result.stress_Pa[0, 0, 0] = 1.0
        too_many = np.zeros((reference.MAX_QUERIES + 1, 3), dtype=np.float64)
        with self.assertRaisesRegex(reference.FieldContractError, "budgetRejected"):
            reference.evaluate_field(too_many)


def math_pi_times(values: np.ndarray) -> np.ndarray:
    return np.pi * values


if __name__ == "__main__":
    unittest.main()
