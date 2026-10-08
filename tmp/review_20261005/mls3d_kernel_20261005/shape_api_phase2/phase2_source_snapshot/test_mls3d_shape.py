"""Source-only acceptance checks for the phase2 direct MLS shape API."""

from __future__ import annotations

import numpy as np
import pytest

from mls3d_mls import MLS3DConfig, MLSStatus, evaluate_mls3d
from mls3d_shape import shape_functions_3d


def grid_points(z_values: np.ndarray | None = None) -> np.ndarray:
    axis = np.linspace(-1.0, 1.0, 5)
    z_axis = axis if z_values is None else z_values
    xx, yy, zz = np.meshgrid(axis, axis, z_axis, indexing="ij")
    return np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))


def proper_rotation(seed: int = 20261005) -> np.ndarray:
    q, _ = np.linalg.qr(np.random.default_rng(seed).normal(size=(3, 3)))
    if np.linalg.det(q) < 0.0:
        q[:, 0] *= -1.0
    return q


def config(
    radius: float = 2.5,
    origin: tuple[float, float, float] = (0.0, 0.0, 0.0),
    scale: float = 1.0,
    max_condition_number: float = 1.0e12,
) -> MLS3DConfig:
    return MLS3DConfig(origin, scale, radius, max_condition_number=max_condition_number)


def test_partition_unity_and_affine_coordinate_reproduction() -> None:
    points = grid_points()
    queries = np.array([[0.13, -0.17, 0.21], [-0.23, 0.08, -0.11]])
    result = shape_functions_3d(points, queries, config())

    assert result.statuses == (MLSStatus.OK, MLSStatus.OK)
    assert result.weights.shape == (2, len(points))
    assert result.gradients.shape == (2, len(points), 3)
    assert result.support_counts.shape == (2,)
    assert result.reasons == (None, None)
    np.testing.assert_allclose(result.weights.sum(axis=1), 1.0, atol=2.0e-12)
    np.testing.assert_allclose(result.gradients.sum(axis=1), 0.0, atol=3.0e-12)
    np.testing.assert_allclose(result.weights @ points, queries, atol=3.0e-12)
    coordinate_gradient = np.einsum("qnk,ni->qik", result.gradients, points)
    np.testing.assert_allclose(coordinate_gradient, np.broadcast_to(np.eye(3), (2, 3, 3)), atol=4.0e-12)
    frozen_weights = result.weights.copy()
    for output in (result.weights, result.gradients, result.support_counts, result.condition_numbers, result.ok):
        assert not output.flags.writeable
        with pytest.raises(ValueError):
            output.setflags(write=True)
    points[:] = 99.0
    queries[:] = -99.0
    np.testing.assert_array_equal(result.weights, frozen_weights)


def test_direct_shapes_reproduce_scalar_vector_mls_values_and_gradients() -> None:
    points = grid_points()
    queries = np.array([[0.13, -0.17, 0.21], [-0.23, 0.08, -0.11]])
    scalar = 0.7 + np.sin(points[:, 0]) + points[:, 1] * points[:, 2] - 0.2 * points[:, 0] ** 2
    vector = np.column_stack(
        (
            np.sin(points[:, 0]),
            points[:, 1] * points[:, 2],
            points[:, 0] ** 2 - 0.2 * points[:, 1],
        )
    )
    cfg = config()
    shapes = shape_functions_3d(points, queries, cfg)
    scalar_reference = evaluate_mls3d(points, scalar, queries, cfg)
    vector_reference = evaluate_mls3d(points, vector, queries, cfg)

    assert shapes.statuses == scalar_reference.statuses == vector_reference.statuses == (MLSStatus.OK, MLSStatus.OK)
    np.testing.assert_allclose(shapes.weights @ scalar, scalar_reference.values[:, 0], atol=3.0e-12)
    scalar_gradient = np.einsum("qnk,n->qk", shapes.gradients, scalar)
    np.testing.assert_allclose(scalar_gradient, scalar_reference.gradients[:, 0, :], atol=4.0e-12)
    np.testing.assert_allclose(shapes.weights @ vector, vector_reference.values, atol=3.0e-12)
    vector_gradient = np.einsum("qnk,ni->qik", shapes.gradients, vector)
    np.testing.assert_allclose(vector_gradient, vector_reference.gradients, atol=5.0e-12)


def test_shape_weights_and_gradients_are_rotation_covariant() -> None:
    rng = np.random.default_rng(4321)
    points = rng.uniform(-1.0, 1.0, size=(100, 3))
    queries = np.array([[0.11, -0.07, 0.23]])
    origin = np.array([0.17, -0.12, 0.08])
    rotation = proper_rotation(2026)
    before = shape_functions_3d(points, queries, config(4.0, tuple(origin), 1.7))
    after = shape_functions_3d(
        points @ rotation.T,
        queries @ rotation.T,
        config(4.0, tuple(rotation @ origin), 1.7),
    )

    assert before.statuses == after.statuses == (MLSStatus.OK,)
    np.testing.assert_allclose(after.weights, before.weights, atol=4.0e-12)
    # With X'=R X, the coordinate-gradient column transforms as g'=R g.
    expected = np.einsum("ji,qni->qnj", rotation, before.gradients)
    np.testing.assert_allclose(after.gradients, expected, atol=5.0e-12)


def test_shape_gradient_matches_independent_central_difference() -> None:
    points = grid_points()
    query = np.array([[0.13, -0.17, 0.21]])
    cfg = config()
    analytic = shape_functions_3d(points, query, cfg)
    assert analytic.statuses == (MLSStatus.OK,)

    h = 2.0e-6
    finite_difference = np.empty((len(points), 3))
    for axis in range(3):
        step = np.zeros((1, 3))
        step[0, axis] = h
        plus = shape_functions_3d(points, query + step, cfg)
        minus = shape_functions_3d(points, query - step, cfg)
        assert plus.statuses == minus.statuses == (MLSStatus.OK,)
        finite_difference[:, axis] = (plus.weights[0] - minus.weights[0]) / (2.0 * h)
    np.testing.assert_allclose(analytic.gradients[0], finite_difference, rtol=3.0e-7, atol=3.0e-8)


def test_poor_rank_duplicate_condition_and_switch_guards_return_unknown() -> None:
    query = np.zeros((1, 3))
    sparse = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]])
    poor = shape_functions_3d(sparse, query, config(0.5))
    assert poor.statuses == (MLSStatus.UNKNOWN_POOR_SUPPORT,)
    assert poor.reasons == (MLSStatus.UNKNOWN_POOR_SUPPORT.value,)
    assert np.all(np.isnan(poor.weights[0]))

    planar = shape_functions_3d(grid_points(np.array([0.0])), query, config())
    assert planar.statuses == (MLSStatus.UNKNOWN_RANK_DEFICIENT,)

    duplicate_points = np.array(
        [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5]]
    )
    duplicate = shape_functions_3d(duplicate_points, query, config(2.0))
    assert duplicate.statuses == (MLSStatus.UNKNOWN_DUPLICATE_SOURCES,)

    ill_conditioned = shape_functions_3d(grid_points(), query, config(max_condition_number=1.0000001))
    assert ill_conditioned.statuses == (MLSStatus.UNKNOWN_ILL_CONDITIONED,)

    boundary_points = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.2, 0.0], [0.0, 0.0, 0.2]]
    )
    boundary = shape_functions_3d(boundary_points, query, config(1.0))
    assert boundary.statuses == (MLSStatus.UNKNOWN_SUPPORT_BOUNDARY,)

    masked = shape_functions_3d(
        grid_points(), query, config(), visibility_mask=np.ones(125, dtype=bool), mask_locally_constant=False
    )
    assert masked.statuses == (MLSStatus.UNKNOWN_VISIBILITY_SWITCH,)
    unspecified = shape_functions_3d(grid_points(), query, config(), visibility_mask=np.ones(125, dtype=bool))
    assert unspecified.statuses == (MLSStatus.UNKNOWN_VISIBILITY_SWITCH,)


def test_fixed_visibility_mask_keeps_hidden_source_rows_zero() -> None:
    points = grid_points()
    query = np.array([[0.13, -0.17, 0.21]])
    visibility = np.ones(len(points), dtype=bool)
    visibility[0] = False
    original_visibility = visibility.copy()

    result = shape_functions_3d(
        points,
        query,
        config(),
        visibility_mask=visibility,
        mask_locally_constant=True,
    )
    assert result.statuses == (MLSStatus.OK,)
    np.testing.assert_array_equal(visibility, original_visibility)
    np.testing.assert_array_equal(result.weights[:, ~visibility], 0.0)
    np.testing.assert_array_equal(result.gradients[:, ~visibility, :], 0.0)


def test_source_and_query_budgets_reject_before_shape_output_allocation() -> None:
    query = np.zeros((1, 3))
    points = np.array([[0.0, 0.0, 0.0]])
    many_sources = np.broadcast_to(points, (501, 3))
    with pytest.raises(ValueError, match="budgetRejected.*source_positions"):
        shape_functions_3d(many_sources, query, config())

    many_queries = np.broadcast_to(np.zeros((1, 3)), (501, 3))
    four_points = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    with pytest.raises(ValueError, match="budgetRejected.*query_positions"):
        shape_functions_3d(four_points, many_queries, config(3.0))

    huge_sources = np.broadcast_to(points, (10**9, 3))
    with pytest.raises(ValueError, match="budgetRejected.*source_positions"):
        shape_functions_3d(huge_sources, query, config())

