"""Source tests for exact finite-facet visibility stability; not run before review."""

from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import stability3d  # noqa: E402


def _fixture() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    sources = np.asarray(((0.0, 0.0, 1.0),), dtype=np.float64)
    queries = np.asarray(((0.0, 0.0, -1.0),), dtype=np.float64)
    triangles = np.asarray(
        (((-1.0, -1.0, 0.0), (1.0, -1.0, 0.0), (0.0, 1.0, 0.0)),),
        dtype=np.float64,
    )
    return sources, queries, triangles


def _classify(source: tuple[float, float, float], query: tuple[float, float, float], triangle: np.ndarray):
    return stability3d.certify_local_visibility_stability(
        np.asarray((source,), dtype=np.float64),
        np.asarray((query,), dtype=np.float64),
        np.asarray((triangle,), dtype=np.float64),
    )


def test_strict_interior_crossing_and_strict_plane_miss_are_certified() -> None:
    sources, queries, triangles = _fixture()
    crossing = stability3d.certify_local_visibility_stability(sources, queries, triangles)
    assert crossing.pair_statuses[0, 0] == "BLOCKED"
    assert crossing.stable_blocked_mask[0, 0]
    assert crossing.facet_statuses[0, 0, 0] == "STABLE_BLOCKED"
    assert crossing.facet_reasons[0, 0, 0] == "STABLE_STRICT_INTERIOR_CROSSING"

    miss = _classify((3.0, 0.0, 1.0), (3.0, 0.0, -1.0), triangles[0])
    assert miss.pair_statuses[0, 0] == "VISIBLE"
    assert miss.stable_visible_mask[0, 0]
    assert miss.facet_reasons[0, 0, 0] == "STABLE_STRICT_FINITE_MISS"


def test_same_strict_side_is_certified_visible() -> None:
    _, _, triangles = _fixture()
    result = _classify((0.0, 0.0, 1.0), (0.0, 0.0, 2.0), triangles[0])
    assert result.pair_statuses[0, 0] == "VISIBLE"
    assert result.stable_visible_mask[0, 0]
    assert result.facet_reasons[0, 0, 0] == "STABLE_SAME_STRICT_SIDE"


def test_plane_endpoint_far_outside_finite_facet_is_a_strict_miss() -> None:
    _, _, triangles = _fixture()
    source_front = _classify((3.0, 0.0, 0.0), (3.0, 0.0, 1.0), triangles[0])
    query_front = _classify((3.0, 0.0, 1.0), (3.0, 0.0, 0.0), triangles[0])
    assert source_front.pair_statuses[0, 0] == "VISIBLE"
    assert source_front.facet_reasons[0, 0, 0] == "STABLE_SOURCE_PLANE_FINITE_MISS"
    assert query_front.pair_statuses[0, 0] == "VISIBLE"
    assert query_front.facet_reasons[0, 0, 0] == "STABLE_QUERY_PLANE_FINITE_MISS"


def test_plane_front_inside_facet_edge_and_coplanar_cases_are_unknown() -> None:
    _, _, triangles = _fixture()
    triangle = triangles[0]
    source_front = _classify((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), triangle)
    query_front = _classify((0.0, 0.0, 1.0), (0.0, 0.0, 0.0), triangle)
    edge = _classify((0.0, -1.0, 1.0), (0.0, -1.0, -1.0), triangle)
    coplanar = _classify((3.0, 0.0, 0.0), (4.0, 0.0, 0.0), triangle)
    for result in (source_front, query_front, edge, coplanar):
        assert result.pair_statuses[0, 0] == "UNKNOWN"
        assert not result.stable_visible_mask[0, 0]
        assert result.unknown_mask[0, 0]
    assert source_front.facet_reasons[0, 0, 0] == "UNKNOWN_SOURCE_FRONT_OR_EDGE"
    assert query_front.facet_reasons[0, 0, 0] == "UNKNOWN_QUERY_FRONT_OR_EDGE"
    assert edge.facet_reasons[0, 0, 0] == "UNKNOWN_TRIANGLE_EDGE_OR_FRONT"
    assert coplanar.facet_reasons[0, 0, 0] == "UNKNOWN_COPLANAR_SEGMENT"


def test_degenerate_facet_is_unknown() -> None:
    degenerate = np.asarray(
        (((-1.0, 0.0, 0.0), (0.0, 0.0, 0.0), (1.0, 0.0, 0.0)),),
        dtype=np.float64,
    )
    result = _classify((0.0, 0.0, 1.0), (0.0, 0.0, -1.0), degenerate[0])
    assert result.pair_statuses[0, 0] == "UNKNOWN"
    assert result.facet_reasons[0, 0, 0] == "UNKNOWN_DEGENERATE_TRIANGLE"


def test_ambiguous_facet_keeps_pair_unknown_even_with_a_definite_blocker() -> None:
    _, _, triangle = _fixture()
    degenerate = np.asarray(((-1.0, 0.0, 0.0), (0.0, 0.0, 0.0), (1.0, 0.0, 0.0)))
    facets = np.stack((triangle[0], degenerate), axis=0)
    result = stability3d.certify_local_visibility_stability(
        np.asarray(((0.0, 0.0, 1.0),)),
        np.asarray(((0.0, 0.0, -1.0),)),
        facets,
    )
    assert result.facet_statuses[0, 0, 0] == "STABLE_BLOCKED"
    assert result.facet_statuses[0, 0, 1] == "UNKNOWN"
    assert result.pair_statuses[0, 0] == "UNKNOWN"
    assert result.pair_reasons[0, 0] == "UNKNOWN_DEGENERATE_TRIANGLE"
    assert not result.stable_blocked_mask[0, 0]


def test_exact_certificate_that_disagrees_with_float_tolerance_stays_unknown() -> None:
    _, _, triangles = _fixture()
    result = _classify((0.0, 0.0, 1.0e-14), (0.0, 0.0, -1.0e-14), triangles[0])
    assert result.facet_statuses[0, 0, 0] == "STABLE_BLOCKED"
    assert result.frozen_float_statuses[0, 0] == "UNKNOWN"
    assert result.pair_statuses[0, 0] == "UNKNOWN"
    assert result.pair_reasons[0, 0] == "UNKNOWN_FLOAT_STATUS_MISMATCH"


def test_no_facets_are_vacuously_visible_and_identity() -> None:
    sources = np.asarray(((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)), dtype=np.float64)
    queries = np.asarray(((2.0, 0.0, 0.0),), dtype=np.float64)
    facets = np.empty((0, 3, 3), dtype=np.float64)
    result = stability3d.certify_local_visibility_stability(sources, queries, facets)
    assert result.pair_statuses.shape == (1, 2)
    assert np.all(result.pair_statuses == "VISIBLE")
    assert np.all(result.stable_visible_mask)
    assert result.facet_statuses.shape == (1, 2, 0)


def _generic_proper_90_degree_rotation() -> np.ndarray:
    axis = np.asarray((1.0, 2.0, 3.0), dtype=np.float64)
    axis /= np.linalg.norm(axis)
    x, y, z = axis
    skew = np.asarray(((0.0, -z, y), (z, 0.0, -x), (-y, x, 0.0)), dtype=np.float64)
    # Rodrigues with theta = pi/2: a proper rotation about a non-axis-aligned axis.
    return np.outer(axis, axis) + skew


def test_generic_rotation_winding_reversal_and_si_scale_preserve_decisive_classes() -> None:
    sources, queries, triangles = _fixture()
    metre_scale = stability3d.certify_local_visibility_stability(sources, queries, triangles)
    scale = 1.0e-3
    base_source = sources * scale
    base_query = queries * scale
    base_triangles = triangles * scale
    baseline = stability3d.certify_local_visibility_stability(base_source, base_query, base_triangles)
    assert metre_scale.pair_statuses[0, 0] == baseline.pair_statuses[0, 0]
    assert baseline.pair_statuses[0, 0] == "BLOCKED"

    rotation = _generic_proper_90_degree_rotation()
    rotated_source = base_source @ rotation.T
    rotated_query = base_query @ rotation.T
    rotated_triangles = np.matmul(base_triangles, rotation.T)
    rotated = stability3d.certify_local_visibility_stability(
        rotated_source, rotated_query, rotated_triangles
    )
    reversed_winding = stability3d.certify_local_visibility_stability(
        rotated_source, rotated_query, rotated_triangles[:, ::-1, :].copy()
    )
    assert rotated.pair_statuses[0, 0] == baseline.pair_statuses[0, 0]
    assert reversed_winding.pair_statuses[0, 0] == baseline.pair_statuses[0, 0]
    assert rotated.stable_blocked_mask[0, 0]
    assert reversed_winding.stable_blocked_mask[0, 0]


def test_intermediate_fraction_budget_rejects_before_oversized_product() -> None:
    arithmetic = stability3d._BudgetArithmetic()
    left = Fraction(1, 1 << 1100)
    right = Fraction(1, 1 << 1100)
    with pytest.raises(stability3d._IntermediateBudgetExceeded):
        arithmetic.mul(left, right)


def test_input_fraction_bit_budget_fails_closed() -> None:
    sources, queries, triangles = _fixture()
    sources[0, 0] = np.nextafter(np.float64(0.0), np.float64(1.0))
    result = stability3d.certify_local_visibility_stability(sources, queries, triangles)
    assert result.fraction_input_max_bits > stability3d.MAX_INPUT_FRACTION_BITS
    assert result.pair_statuses[0, 0] == "UNKNOWN"
    assert result.pair_reasons[0, 0] == "UNKNOWN_INPUT_BIT_BUDGET"
    assert not result.stable_visible_mask[0, 0]


def test_segment_facet_work_budget_fails_closed_without_facet_array() -> None:
    sources = np.zeros((129, 3), dtype=np.float64)
    queries = np.ones((127, 3), dtype=np.float64)
    triangles = np.zeros((2, 3, 3), dtype=np.float64)
    result = stability3d.certify_local_visibility_stability(sources, queries, triangles)
    assert result.pair_statuses.shape == (127, 129)
    assert np.all(result.pair_statuses == "UNKNOWN")
    assert np.all(result.pair_reasons == "UNKNOWN_WORK_BUDGET")
    assert result.facet_statuses is None
    assert result.facet_reasons is None
    assert np.all(result.frozen_float_statuses == "NOT_EVALUATED")


def test_query_source_pair_budget_fails_closed_without_facet_array() -> None:
    sources = np.zeros((129, 3), dtype=np.float64)
    queries = np.ones((128, 3), dtype=np.float64)
    triangles = np.empty((0, 3, 3), dtype=np.float64)
    result = stability3d.certify_local_visibility_stability(sources, queries, triangles)
    assert result.pair_statuses.shape == (128, 129)
    assert np.all(result.pair_reasons == "UNKNOWN_PAIR_BUDGET")
    assert result.facet_statuses is None


def test_result_arrays_are_bytes_backed_immutable() -> None:
    sources, queries, triangles = _fixture()
    result = stability3d.certify_local_visibility_stability(sources, queries, triangles)
    for array in (
        result.pair_statuses,
        result.pair_reasons,
        result.stable_visible_mask,
        result.stable_blocked_mask,
        result.unknown_mask,
        result.frozen_float_statuses,
        result.facet_statuses,
        result.facet_reasons,
    ):
        assert array is not None
        assert not array.flags.writeable
        with pytest.raises(ValueError):
            array.setflags(write=True)


def test_preflight_rejects_subclasses_and_oversized_rows_before_conversion() -> None:
    class HookArray(np.ndarray):
        def __array__(self, *args, **kwargs):
            raise AssertionError("array conversion hook must not run")

    sources, queries, triangles = _fixture()
    subclassed = sources.view(HookArray)
    with pytest.raises(ValueError, match="exact numpy.ndarray"):
        stability3d.certify_local_visibility_stability(subclassed, queries, triangles)

    too_many_sources = np.zeros((stability3d.MAX_SOURCES + 1, 3), dtype=np.float64)
    with pytest.raises(ValueError, match="budgetRejected"):
        stability3d.certify_local_visibility_stability(too_many_sources, queries, triangles)


def test_nonfinite_and_non_float64_inputs_are_rejected() -> None:
    sources, queries, triangles = _fixture()
    sources[0, 0] = np.nan
    with pytest.raises(ValueError, match="finite"):
        stability3d.certify_local_visibility_stability(sources, queries, triangles)
    with pytest.raises(ValueError, match="float64"):
        stability3d.certify_local_visibility_stability(
            np.asarray(((0, 0, 1),), dtype=np.float32), queries, triangles
        )

