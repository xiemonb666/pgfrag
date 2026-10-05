"""Source-only acceptance tests. Deliberately not run before root review."""

from __future__ import annotations

import unittest

import numpy as np

from mls3d_mls import (
    MLS3DConfig,
    MLSStatus,
    evaluate_mls3d,
    green_lagrange_strain,
)


def grid_points(z_values: np.ndarray | None = None) -> np.ndarray:
    axis = np.linspace(-1.0, 1.0, 5)
    z_axis = axis if z_values is None else z_values
    xx, yy, zz = np.meshgrid(axis, axis, z_axis, indexing="ij")
    return np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))


def config(radius: float = 2.5) -> MLS3DConfig:
    return MLS3DConfig((0.0, 0.0, 0.0), 1.0, radius)


def proper_rotation(seed: int = 20261005) -> np.ndarray:
    q, _ = np.linalg.qr(np.random.default_rng(seed).normal(size=(3, 3)))
    if np.linalg.det(q) < 0.0:
        q[:, 0] *= -1.0
    return q


class HostileArrayLike:
    def __init__(self) -> None:
        self.array_calls = 0
        self.tolist_calls = 0

    def __array__(self, *args: object, **kwargs: object) -> np.ndarray:
        self.array_calls += 1
        raise AssertionError("__array__ must not be invoked")

    def tolist(self) -> list[float]:
        self.tolist_calls += 1
        raise AssertionError("tolist must not be invoked")


class LinearMLS3DTests(unittest.TestCase):
    def test_constant_and_affine_reproduction_with_gradient(self) -> None:
        points = grid_points()
        query = np.array([[0.13, -0.17, 0.21]])
        matrix = np.array([[1.2, -0.4, 0.7], [0.1, 2.0, -0.3], [-0.8, 0.2, 0.5]])
        offset = np.array([0.4, -1.0, 2.5])
        result = evaluate_mls3d(points, points @ matrix.T + offset, query, config())
        self.assertEqual(result.statuses, (MLSStatus.OK,))
        np.testing.assert_allclose(result.values[0], matrix @ query[0] + offset, atol=2.0e-12)
        np.testing.assert_allclose(result.gradients[0], matrix, atol=2.0e-12)

        constant = np.array([3.0, -2.0, 0.25])
        result_constant = evaluate_mls3d(
            points, np.broadcast_to(constant, (len(points), 3)), query, config()
        )
        self.assertEqual(result_constant.statuses, (MLSStatus.OK,))
        np.testing.assert_allclose(result_constant.values[0], constant, atol=2.0e-12)
        np.testing.assert_allclose(result_constant.gradients[0], 0.0, atol=2.0e-12)

    def test_rigid_translation_and_finite_rotation_have_zero_strain(self) -> None:
        points = grid_points()
        query = np.array([[0.13, -0.17, 0.21]])
        translation = np.array([0.2, -0.5, 1.1])
        translated = evaluate_mls3d(
            points,
            np.broadcast_to(translation, (len(points), 3)),
            query,
            config(),
        )
        self.assertEqual(translated.statuses, (MLSStatus.OK,))
        np.testing.assert_allclose(translated.gradients[0], 0.0, atol=2.0e-12)
        np.testing.assert_allclose(green_lagrange_strain(translated.gradients[0]), 0.0, atol=2.0e-12)

        rotation = proper_rotation(71)
        rigid_displacement = points @ rotation.T - points
        rotated = evaluate_mls3d(points, rigid_displacement, query, config())
        self.assertEqual(rotated.statuses, (MLSStatus.OK,))
        np.testing.assert_allclose(rotated.gradients[0], rotation - np.eye(3), atol=2.0e-12)
        np.testing.assert_allclose(green_lagrange_strain(rotated.gradients[0]), 0.0, atol=3.0e-12)

    def test_covariance_under_generic_proper_rotation(self) -> None:
        rng = np.random.default_rng(318)
        points = rng.uniform(-1.0, 1.0, size=(90, 3))
        values = rng.normal(size=(90, 3))
        query = np.array([[0.11, -0.07, 0.23]])
        origin = np.array([0.17, -0.12, 0.08])
        rotation = proper_rotation(2026)
        cfg = MLS3DConfig(tuple(origin), 1.7, 4.0)
        before = evaluate_mls3d(points, values, query, cfg)
        after = evaluate_mls3d(
            points @ rotation.T,
            values @ rotation.T,
            query @ rotation.T,
            MLS3DConfig(tuple(rotation @ origin), 1.7, 4.0),
        )
        self.assertEqual(before.statuses, (MLSStatus.OK,))
        self.assertEqual(after.statuses, (MLSStatus.OK,))
        np.testing.assert_allclose(after.values[0], rotation @ before.values[0], atol=2.0e-12)
        expected_gradient = rotation @ before.gradients[0] @ rotation.T
        np.testing.assert_allclose(after.gradients[0], expected_gradient, atol=3.0e-12)

    def test_analytic_gradient_matches_independent_finite_difference(self) -> None:
        points = grid_points()
        values = np.column_stack(
            (
                0.3 + points[:, 0] ** 2 - 0.2 * points[:, 1],
                points[:, 0] * points[:, 2] + 0.5 * points[:, 1],
                np.sin(points[:, 0]) + points[:, 1] * points[:, 2],
            )
        )
        query = np.array([[0.13, -0.17, 0.21]])
        cfg = config()
        analytic = evaluate_mls3d(points, values, query, cfg)
        self.assertEqual(analytic.statuses, (MLSStatus.OK,))
        h = 2.0e-6
        finite_difference = np.empty((3, 3))
        for axis in range(3):
            step = np.zeros((1, 3))
            step[0, axis] = h
            plus = evaluate_mls3d(points, values, query + step, cfg)
            minus = evaluate_mls3d(points, values, query - step, cfg)
            self.assertEqual(plus.statuses, (MLSStatus.OK,))
            self.assertEqual(minus.statuses, (MLSStatus.OK,))
            finite_difference[:, axis] = (plus.values[0] - minus.values[0]) / (2.0 * h)
        np.testing.assert_allclose(analytic.gradients[0], finite_difference, rtol=2.0e-7, atol=2.0e-8)

    def test_poor_support_rank_and_duplicates_return_unknown(self) -> None:
        query = np.zeros((1, 3))
        sparse = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]])
        poor = evaluate_mls3d(sparse, np.ones((4, 1)), query, config(radius=0.5))
        self.assertEqual(poor.statuses, (MLSStatus.UNKNOWN_POOR_SUPPORT,))
        self.assertTrue(np.isnan(poor.values[0, 0]))

        planar = grid_points(np.array([0.0]))
        rank = evaluate_mls3d(planar, np.ones((len(planar), 1)), query, config())
        self.assertEqual(rank.statuses, (MLSStatus.UNKNOWN_RANK_DEFICIENT,))

        duplicate_positions = np.array(
            [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5]]
        )
        duplicate = evaluate_mls3d(
            duplicate_positions, np.ones((5, 1)), query, config(radius=2.0)
        )
        self.assertEqual(duplicate.statuses, (MLSStatus.UNKNOWN_DUPLICATE_SOURCES,))

        switch_sources = np.array(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.2, 0.0], [0.0, 0.0, 0.2]]
        )
        at_support_boundary = evaluate_mls3d(
            switch_sources,
            np.ones((4, 1)),
            query,
            MLS3DConfig((0.0, 0.0, 0.0), 1.0, 1.0),
        )
        self.assertEqual(at_support_boundary.statuses, (MLSStatus.UNKNOWN_SUPPORT_BOUNDARY,))

        masked = evaluate_mls3d(
            grid_points(),
            np.ones((125, 1)),
            query,
            config(),
            visibility_mask=np.ones(125, dtype=bool),
            mask_locally_constant=False,
        )
        self.assertEqual(masked.statuses, (MLSStatus.UNKNOWN_VISIBILITY_SWITCH,))
        unspecified_mask = evaluate_mls3d(
            grid_points(), np.ones((125, 1)), query, config(), visibility_mask=np.ones(125, dtype=bool)
        )
        self.assertEqual(unspecified_mask.statuses, (MLSStatus.UNKNOWN_VISIBILITY_SWITCH,))

    def test_input_budgets_reject_before_conversion_or_output_allocation(self) -> None:
        query = np.zeros((1, 3))
        points = np.array(
            [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]]
        )

        with self.assertRaisesRegex(ValueError, "budgetRejected.*source_positions"):
            evaluate_mls3d(np.zeros((501, 3)), np.zeros((501, 1)), query, config())
        with self.assertRaisesRegex(ValueError, "budgetRejected.*query_positions"):
            evaluate_mls3d(points, np.zeros((4, 1)), np.zeros((501, 3)), config())
        with self.assertRaisesRegex(ValueError, "budgetRejected.*components"):
            evaluate_mls3d(points, np.zeros((4, 4)), query, config())

        # This is a zero-stride view with a huge logical component dimension.
        # Shape inspection must reject it before a finite scan or dtype copy.
        huge_view = np.broadcast_to(np.zeros((4, 1)), (4, 10**9))
        self.assertEqual(huge_view.strides[1], 0)
        with self.assertRaisesRegex(ValueError, "budgetRejected.*components"):
            evaluate_mls3d(points, huge_view, query, config())

        hostile = HostileArrayLike()
        with self.assertRaisesRegex(ValueError, "plain numpy.ndarray"):
            evaluate_mls3d(points, hostile, query, config())  # type: ignore[arg-type]
        self.assertEqual(hostile.array_calls, 0)
        self.assertEqual(hostile.tolist_calls, 0)


if __name__ == "__main__":
    unittest.main()
