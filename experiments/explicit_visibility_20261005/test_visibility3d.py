"""Counterexample-focused source tests; intentionally not executed in this pass."""

from __future__ import annotations

import inspect
from pathlib import Path
import sys

import numpy as np
import pytest


sys.path.insert(0, str(Path(__file__).resolve().parent))
import visibility3d as visibility  # noqa: E402


def _triangle() -> np.ndarray:
    return np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]], dtype=np.float64)


def _classify_one(source: list[float], query: list[float], triangles: np.ndarray | None = None) -> visibility.VisibilityStatus:
    sources = np.array([source], dtype=np.float64)
    queries = np.array([query], dtype=np.float64)
    facets = _triangle() if triangles is None else triangles
    return visibility.classify_visibility(sources, queries, facets).status_at(0, 0)


def test_crossing_interior_is_blocked_and_plane_miss_is_visible() -> None:
    sources = np.array([[0.2, 0.2, 1.0]])
    queries = np.array([[0.2, 0.2, -1.0], [2.0, 2.0, -1.0]])
    result = visibility.classify_visibility(sources, queries, _triangle())

    assert result.statuses.tolist() == [["BLOCKED"], ["VISIBLE"]]


def test_query_major_pairing_distinguishes_each_source_for_each_query() -> None:
    sources = np.array([[0.2, 0.2, 1.0], [2.0, 2.0, 1.0]])
    queries = np.array([[0.2, 0.2, -1.0], [2.0, 2.0, -1.0]])
    result = visibility.classify_visibility(sources, queries, _triangle())

    assert result.statuses.tolist() == [["BLOCKED", "VISIBLE"], ["VISIBLE", "VISIBLE"]]
    assert result.visible_mask.tolist() == [[False, True], [True, True]]
    assert result.blocked_mask.tolist() == [[True, False], [False, False]]
    assert not bool(np.any(result.unknown_mask))


def test_segment_crossing_triangle_plane_outside_finite_facet_is_visible() -> None:
    assert _classify_one([0.8, 0.8, 1.0], [0.8, 0.8, -1.0]) is visibility.VisibilityStatus.VISIBLE


def test_rotated_and_translated_world_layout_preserves_decisive_classifications() -> None:
    rotation = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]])
    translation = np.array([3.0, -4.0, 8.0])
    source = np.array([[0.2, 0.2, 1.0], [2.0, 2.0, 1.0]])
    query = np.array([[0.2, 0.2, -1.0], [2.0, 2.0, -1.0]])
    triangles = _triangle()

    before = visibility.classify_visibility(source, query, triangles)
    after = visibility.classify_visibility(
        source @ rotation.T + translation,
        query @ rotation.T + translation,
        triangles @ rotation.T + translation,
    )

    assert after.statuses.tolist() == before.statuses.tolist()


def test_segment_direction_does_not_select_a_cut_side() -> None:
    triangles = _triangle()
    forward = _classify_one([0.2, 0.2, 1.0], [0.2, 0.2, -1.0], triangles)
    reverse = _classify_one([0.2, 0.2, -1.0], [0.2, 0.2, 1.0], triangles)

    assert forward is visibility.VisibilityStatus.BLOCKED
    assert reverse is forward
    assert "identity" not in inspect.signature(visibility.classify_visibility).parameters
    assert "tetrahedra" not in inspect.signature(visibility.classify_visibility).parameters


def test_endpoint_on_or_near_finite_front_is_unknown() -> None:
    assert _classify_one([0.2, 0.2, 0.0], [0.2, 0.2, 1.0]) is visibility.VisibilityStatus.UNKNOWN
    assert _classify_one([0.2, 0.2, 5.0e-13], [0.2, 0.2, 1.0]) is visibility.VisibilityStatus.UNKNOWN


def test_endpoint_in_triangle_plane_but_far_outside_facet_is_not_automatically_unknown() -> None:
    assert _classify_one([10.0, 10.0, 0.0], [10.0, 10.0, 0.0]) is visibility.VisibilityStatus.VISIBLE


def test_edge_grazing_is_unknown() -> None:
    assert _classify_one([0.5, 0.0, 1.0], [0.5, 0.0, -1.0]) is visibility.VisibilityStatus.UNKNOWN
    assert _classify_one([0.5, 1.0e-11, 1.0], [0.5, 1.0e-11, -1.0]) is visibility.VisibilityStatus.UNKNOWN


def test_coplanar_overlap_is_unknown_but_coplanar_far_miss_is_visible() -> None:
    assert _classify_one([-1.0, 0.2, 0.0], [1.0, 0.2, 0.0]) is visibility.VisibilityStatus.UNKNOWN
    assert _classify_one([-1.0, 10.0, 0.0], [1.0, 10.0, 0.0]) is visibility.VisibilityStatus.VISIBLE


def test_degenerate_facet_fails_closed_as_unknown() -> None:
    degenerate = np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]])
    assert _classify_one([0.2, 0.2, 1.0], [0.2, 0.2, -1.0], degenerate) is visibility.VisibilityStatus.UNKNOWN


def test_proven_blocker_dominates_a_different_ambiguous_facet() -> None:
    triangles = np.concatenate(
        (_triangle(), np.array([[[0.0, 0.0, 0.5], [1.0, 0.0, 0.5], [2.0, 0.0, 0.5]]])),
        axis=0,
    )
    assert _classify_one([0.2, 0.2, 1.0], [0.2, 0.2, -1.0], triangles) is visibility.VisibilityStatus.BLOCKED


def test_result_is_immutable_and_inputs_are_unchanged() -> None:
    sources = np.array([[0.2, 0.2, 1.0]])
    queries = np.array([[0.2, 0.2, -1.0]])
    triangles = _triangle()
    snapshots = (sources.copy(), queries.copy(), triangles.copy())

    result = visibility.classify_visibility(sources, queries, triangles)

    np.testing.assert_array_equal(sources, snapshots[0])
    np.testing.assert_array_equal(queries, snapshots[1])
    np.testing.assert_array_equal(triangles, snapshots[2])
    assert result.spatial_stability == "NOT_ESTABLISHED"
    for output in (result.statuses, result.visible_mask, result.blocked_mask, result.unknown_mask):
        with pytest.raises(ValueError):
            output.setflags(write=True)


def test_fake_array_hook_subclass_and_noncontiguous_view_are_rejected() -> None:
    class ArrayBomb(np.ndarray):
        def __array__(self, *args: object, **kwargs: object) -> np.ndarray:
            raise AssertionError("array conversion hook must not run")

    subclass = np.zeros((1, 3), dtype=np.float64).view(ArrayBomb)
    queries = np.zeros((1, 3), dtype=np.float64)
    with pytest.raises(ValueError, match="exact numpy.ndarray"):
        visibility.classify_visibility(subclass, queries, _triangle())

    noncontiguous = np.zeros((1, 6), dtype=np.float64)[:, ::2]
    with pytest.raises(ValueError, match="C-contiguous"):
        visibility.classify_visibility(noncontiguous, queries, _triangle())


def test_oversized_zero_stride_view_is_rejected_from_shape_before_copy() -> None:
    backing = np.zeros(1, dtype=np.float64)
    oversized_view = np.lib.stride_tricks.as_strided(backing, shape=(visibility.MAX_SOURCES + 1, 3), strides=(0, 0))
    queries = np.zeros((1, 3), dtype=np.float64)

    with pytest.raises(ValueError, match="budgetRejected"):
        visibility.classify_visibility(oversized_view, queries, _triangle())


def test_work_pair_budget_is_checked_before_private_input_copies(monkeypatch: pytest.MonkeyPatch) -> None:
    sources = np.zeros((visibility.MAX_SOURCES, 3), dtype=np.float64)
    queries = np.zeros((visibility.MAX_QUERIES, 3), dtype=np.float64)
    triangles = np.zeros((5, 3, 3), dtype=np.float64)

    def forbidden_copy(*args: object, **kwargs: object) -> None:
        pytest.fail("private copies must not be made before the work-pair budget check")

    monkeypatch.setattr(visibility, "_copy_inputs", forbidden_copy)
    with pytest.raises(ValueError, match="budgetRejected"):
        visibility.classify_visibility(sources, queries, triangles)


def test_no_facets_is_an_identity_visibility_map() -> None:
    sources = np.array([[0.2, 0.2, 1.0], [2.0, 2.0, -1.0]])
    queries = np.array([[0.2, 0.2, -1.0], [2.0, 2.0, 1.0]])
    no_facets = np.empty((0, 3, 3), dtype=np.float64)

    result = visibility.classify_visibility(sources, queries, no_facets)

    assert result.statuses.tolist() == [["VISIBLE", "VISIBLE"], ["VISIBLE", "VISIBLE"]]
    assert bool(np.all(result.visible_mask))
    assert not bool(np.any(result.blocked_mask))
    assert not bool(np.any(result.unknown_mask))


def test_generic_proper_nonaxis_quarter_turn_winding_and_si_scale_covariance() -> None:
    axis = np.array([1.0, 2.0, 3.0])
    axis /= np.linalg.norm(axis)
    skew = np.array(
        [[0.0, -axis[2], axis[1]], [axis[2], 0.0, -axis[0]], [-axis[1], axis[0], 0.0]],
        dtype=np.float64,
    )
    rotation = np.eye(3) + skew + skew @ skew  # Rodrigues, 90 degrees about a non-axis-aligned unit vector.
    np.testing.assert_allclose(rotation.T @ rotation, np.eye(3), atol=2.0e-15)
    assert np.linalg.det(rotation) > 0.0
    expected = [["BLOCKED", "VISIBLE"], ["VISIBLE", "VISIBLE"]]

    for scale in (1.0, 1.0e-3):
        sources = np.array([[0.2, 0.2, 1.0], [2.0, 2.0, 1.0]], dtype=np.float64) * scale
        queries = np.array([[0.2, 0.2, -1.0], [2.0, 2.0, -1.0]], dtype=np.float64) * scale
        triangles = _triangle() * scale
        rewound = triangles[:, [0, 2, 1], :].copy()

        baseline = visibility.classify_visibility(sources, queries, triangles).statuses.tolist()
        rotated = visibility.classify_visibility(
            sources @ rotation.T,
            queries @ rotation.T,
            triangles @ rotation.T,
        ).statuses.tolist()
        reversed_winding = visibility.classify_visibility(sources, queries, rewound).statuses.tolist()
        rotated_reversed = visibility.classify_visibility(
            sources @ rotation.T,
            queries @ rotation.T,
            rewound @ rotation.T,
        ).statuses.tolist()

        assert baseline == expected
        assert rotated == baseline
        assert reversed_winding == baseline
        assert rotated_reversed == baseline
