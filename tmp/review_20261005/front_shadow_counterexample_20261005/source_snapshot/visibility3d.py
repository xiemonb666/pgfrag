"""Bounded, source-only visibility checks for explicit triangular crack facets.

This module answers one geometric question: does a finite source-to-query
segment cross the relative interior of any supplied triangle? It has no crack
identity, tetrahedral adjacency, material model, solver hook, nucleation or
growth rule, or ceramic physics. It is not mechanics-coupled. Spatial
derivative/local-stability evidence is explicitly NOT_ESTABLISHED; a fixed
visibility pattern by itself is not evidence that visibility is locally stable.

Coordinates and the absolute tolerance use SI meters. Public inputs must be
exact, C-contiguous NumPy ndarrays and are bounded before any value scan or
copy. Tolerance is the larger of ``absolute_tolerance_m + relative_tolerance
* local_length`` and a float64 coordinate-precision guard. The latter grows
with the absolute world-coordinate offset because subtracting nearby points
at a large offset loses precision.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import math

import numpy as np


MAX_SOURCES = 500
MAX_QUERIES = 500
MAX_TRIANGLES = 256
MAX_SEGMENT_TRIANGLE_PAIRS = 1_000_000
SEGMENT_CHUNK_SIZE = 128
MAX_ABS_COORDINATE_M = 1.0e150
_FLOAT64_EPS = np.finfo(np.float64).eps


class VisibilityStatus(str, Enum):
    """Fail-closed classification for one query/source segment."""

    VISIBLE = "VISIBLE"
    BLOCKED = "BLOCKED"
    UNKNOWN = "UNKNOWN"


@dataclass(frozen=True, slots=True)
class VisibilityConfig:
    """Fixed dimensional and relative tolerances for one coordinate layout."""

    absolute_tolerance_m: float = 1.0e-12
    relative_tolerance: float = 1.0e-10

    def __post_init__(self) -> None:
        absolute = _finite_scalar("absolute_tolerance_m", self.absolute_tolerance_m)
        relative = _finite_scalar("relative_tolerance", self.relative_tolerance)
        if absolute <= 0.0:
            raise ValueError("absolute_tolerance_m must be positive")
        if absolute > MAX_ABS_COORDINATE_M:
            raise ValueError(f"absolute_tolerance_m must be <= {MAX_ABS_COORDINATE_M:g} m")
        if not 0.0 < relative < 0.5:
            raise ValueError("relative_tolerance must be in (0, 0.5)")
        object.__setattr__(self, "absolute_tolerance_m", absolute)
        object.__setattr__(self, "relative_tolerance", relative)


@dataclass(frozen=True, slots=True)
class VisibilityResult:
    """Immutable query-major status and visibility arrays.

    Array shape is ``(n_query, n_source)``. Only `visible_mask == True`
    authorizes use as query-local support. UNKNOWN remains distinct from
    BLOCKED and can never be inferred as visible from these masks.
    """

    statuses: np.ndarray
    visible_mask: np.ndarray
    blocked_mask: np.ndarray
    unknown_mask: np.ndarray
    spatial_stability: str = "NOT_ESTABLISHED"

    def status_at(self, query_index: int, source_index: int) -> VisibilityStatus:
        """Return the enum value at one query/source pair."""
        return VisibilityStatus(str(self.statuses[query_index, source_index]))


def _finite_scalar(name: str, value: object) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(
        value, (int, float, np.integer, np.floating)
    ):
        raise ValueError(f"{name} must be a finite real scalar")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be a finite real scalar")
    return result


def _preflight_array(
    name: str,
    value: object,
    *,
    tail_shape: tuple[int, ...],
    max_rows: int,
    max_bytes: int,
) -> np.ndarray:
    """Inspect only safe ndarray metadata and enforce bounds before scanning."""
    if type(value) is not np.ndarray:
        raise ValueError(f"{name} must be an exact numpy.ndarray")
    if value.ndim != len(tail_shape) + 1 or tuple(value.shape[1:]) != tail_shape:
        tail = ", ".join(str(part) for part in tail_shape)
        raise ValueError(f"{name} must have shape (n, {tail})")
    rows = int(value.shape[0])
    if rows > max_rows:
        raise ValueError(f"budgetRejected: {name} has {rows} rows; maximum is {max_rows}")
    if value.dtype.kind not in "fiu" or value.dtype.itemsize > 8:
        raise ValueError(f"{name} must use a real numeric dtype of at most 8 bytes")
    if value.nbytes > max_bytes:
        raise ValueError(f"budgetRejected: {name} is {value.nbytes} bytes; maximum is {max_bytes}")
    if not value.flags.c_contiguous:
        raise ValueError(f"{name} must be C-contiguous")
    return value


def _validate_finite_values(name: str, value: np.ndarray) -> None:
    if not bool(np.isfinite(value).all()):
        raise ValueError(f"{name} must contain only finite coordinates")
    if value.dtype.kind in "iu":
        # Avoid silently collapsing distinct integer coordinates during the
        # later float64 copy. Inputs are small, so this scan stays bounded.
        exact_limit = 1 << 53
        if value.dtype.kind == "u":
            if bool(np.any(value > exact_limit)):
                raise ValueError(f"{name} integer coordinates must be exactly representable in float64")
        elif bool(np.any(value > exact_limit)) or bool(np.any(value < -exact_limit)):
            raise ValueError(f"{name} integer coordinates must be exactly representable in float64")


def _copy_inputs(
    source_positions: np.ndarray,
    query_positions: np.ndarray,
    triangle_facets: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Make bounded private float64 copies after all public preflight checks."""
    return (
        source_positions.astype(np.float64, copy=True),
        query_positions.astype(np.float64, copy=True),
        triangle_facets.astype(np.float64, copy=True),
    )


def _norm2(vectors: np.ndarray) -> np.ndarray:
    return np.hypot(vectors[..., 0], vectors[..., 1])


def _norm3(vectors: np.ndarray) -> np.ndarray:
    return np.hypot(np.hypot(vectors[..., 0], vectors[..., 1]), vectors[..., 2])


def _cross2(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    return left[..., 0] * right[..., 1] - left[..., 1] * right[..., 0]


def _edge_signed_distances(
    points_xy: np.ndarray,
    edge_starts_xy: np.ndarray,
    edge_vectors_xy: np.ndarray,
    edge_lengths: np.ndarray,
) -> np.ndarray:
    relative = points_xy[:, None, :] - edge_starts_xy[None, :, :]
    signed_area = _cross2(edge_vectors_xy[None, :, :], relative)
    return signed_area / edge_lengths[None, :]


def _point_to_edges_2d(
    points_xy: np.ndarray,
    edge_starts_xy: np.ndarray,
    edge_vectors_xy: np.ndarray,
    edge_length_squared: np.ndarray,
) -> np.ndarray:
    relative = points_xy[:, None, :] - edge_starts_xy[None, :, :]
    parameter = np.sum(relative * edge_vectors_xy[None, :, :], axis=2)
    parameter = np.clip(parameter / edge_length_squared[None, :], 0.0, 1.0)
    closest = edge_starts_xy[None, :, :] + parameter[:, :, None] * edge_vectors_xy[None, :, :]
    return _norm2(points_xy[:, None, :] - closest)


def _point_to_segment_2d(
    points_xy: np.ndarray,
    starts_xy: np.ndarray,
    vectors_xy: np.ndarray,
) -> np.ndarray:
    length_squared = np.sum(vectors_xy * vectors_xy, axis=1)
    safe_length_squared = np.where(length_squared > 0.0, length_squared, 1.0)
    relative = points_xy - starts_xy
    parameter = np.sum(relative * vectors_xy, axis=1) / safe_length_squared
    parameter = np.clip(parameter, 0.0, 1.0)
    closest = starts_xy + parameter[:, None] * vectors_xy
    return _norm2(points_xy - closest)


def _projected_segment_near_triangle(
    start_xy: np.ndarray,
    end_xy: np.ndarray,
    edge_starts_xy: np.ndarray,
    edge_vectors_xy: np.ndarray,
    edge_lengths: np.ndarray,
    edge_length_squared: np.ndarray,
    tolerance: np.ndarray,
) -> np.ndarray:
    """Conservative finite-triangle proximity test in the triangle plane."""
    segment_vector = end_xy - start_xy
    segment_length = _norm2(segment_vector)
    start_signed = _edge_signed_distances(start_xy, edge_starts_xy, edge_vectors_xy, edge_lengths)
    end_signed = _edge_signed_distances(end_xy, edge_starts_xy, edge_vectors_xy, edge_lengths)
    start_inside = np.all(start_signed >= 0.0, axis=1)
    end_inside = np.all(end_signed >= 0.0, axis=1)
    near = start_inside | end_inside

    segment_min = np.minimum(start_xy, end_xy)
    segment_max = np.maximum(start_xy, end_xy)
    for edge_index in range(3):
        edge_start = edge_starts_xy[edge_index]
        edge_vector = edge_vectors_xy[edge_index]
        edge_end = edge_start + edge_vector
        edge_start_rows = np.broadcast_to(edge_start, start_xy.shape)
        edge_end_rows = np.broadcast_to(edge_end, start_xy.shape)
        segment_start_dist = _point_to_segment_2d(start_xy, edge_start_rows, np.broadcast_to(edge_vector, start_xy.shape))
        segment_end_dist = _point_to_segment_2d(end_xy, edge_start_rows, np.broadcast_to(edge_vector, start_xy.shape))
        edge_start_dist = _point_to_segment_2d(edge_start_rows, start_xy, segment_vector)
        edge_end_dist = _point_to_segment_2d(edge_end_rows, start_xy, segment_vector)
        edge_distance = np.minimum.reduce((segment_start_dist, segment_end_dist, edge_start_dist, edge_end_dist))

        denominator = _cross2(segment_vector, np.broadcast_to(edge_vector, segment_vector.shape))
        relative_edge_start = edge_start_rows - start_xy
        numerator_segment = _cross2(relative_edge_start, np.broadcast_to(edge_vector, segment_vector.shape))
        numerator_edge = _cross2(relative_edge_start, segment_vector)
        edge_length = edge_lengths[edge_index]
        determinant_tolerance = tolerance * np.maximum(segment_length, edge_length)
        near_parallel = np.abs(denominator) <= determinant_tolerance
        stable_intersection = ~near_parallel
        safe_denominator = np.where(stable_intersection, denominator, 1.0)
        along_segment = numerator_segment / safe_denominator
        along_edge = numerator_edge / safe_denominator
        intersects = stable_intersection & (along_segment >= 0.0) & (along_segment <= 1.0) & (along_edge >= 0.0) & (along_edge <= 1.0)
        line_separation = np.full(segment_length.shape, np.inf, dtype=np.float64)
        nonzero_segment = segment_length > 0.0
        line_separation[nonzero_segment] = (
            np.abs(_cross2(relative_edge_start[nonzero_segment], segment_vector[nonzero_segment]))
            / segment_length[nonzero_segment]
        )
        edge_min = np.minimum(edge_start, edge_end)
        edge_max = np.maximum(edge_start, edge_end)
        boxes_overlap = np.all(
            (segment_min <= edge_max[None, :] + tolerance[:, None])
            & (edge_min[None, :] <= segment_max + tolerance[:, None]),
            axis=1,
        )
        uncertain_parallel = near_parallel & boxes_overlap & (line_separation <= tolerance)
        near |= (edge_distance <= tolerance) | intersects | uncertain_parallel
    return near


def _triangle_frame(
    triangle: np.ndarray,
) -> tuple:
    a, b, c = triangle
    edge_ab = b - a
    edge_ac = c - a
    normal_cross = np.cross(edge_ab, edge_ac)
    area_twice = float(_norm3(normal_cross))
    edge_bc = c - b
    edge_lengths_3d = _norm3(np.stack((edge_ab, edge_bc, -edge_ac), axis=0))
    edge_max = float(np.max(edge_lengths_3d))
    normal = normal_cross / area_twice
    edge_ab_length = float(edge_lengths_3d[0])
    axis_u = edge_ab / edge_ab_length
    axis_v = np.cross(normal, axis_u)
    edge_starts_xy = np.array(((0.0, 0.0), (edge_ab_length, 0.0), (float(np.dot(edge_ac, axis_u)), float(np.dot(edge_ac, axis_v)))), dtype=np.float64)
    edge_vectors_xy = np.roll(edge_starts_xy, -1, axis=0) - edge_starts_xy
    edge_lengths_xy = _norm2(edge_vectors_xy)
    edge_length_squared = np.sum(edge_vectors_xy * edge_vectors_xy, axis=1)
    return (
        a,
        normal,
        axis_u,
        axis_v,
        edge_starts_xy,
        edge_vectors_xy,
        edge_lengths_xy,
        edge_max,
        area_twice,
        edge_length_squared,
    )


def _classify_against_triangle(
    starts: np.ndarray,
    ends: np.ndarray,
    frame: tuple,
    config: VisibilityConfig,
    coordinate_scale: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    (a, normal, axis_u, axis_v, edge_starts_xy, edge_vectors_xy, edge_lengths_xy,
     edge_max, _area_twice, edge_length_squared) = frame
    start_local = starts - a
    end_local = ends - a
    d0 = start_local @ normal
    d1 = end_local @ normal
    start_xy = np.column_stack((start_local @ axis_u, start_local @ axis_v))
    end_xy = np.column_stack((end_local @ axis_u, end_local @ axis_v))
    segment_length = _norm3(ends - starts)
    local_length = np.maximum(segment_length, edge_max)
    tolerance = np.maximum(
        config.absolute_tolerance_m + config.relative_tolerance * local_length,
        (32.0 * _FLOAT64_EPS) * coordinate_scale,
    )

    start_edge_distances = _point_to_edges_2d(start_xy, edge_starts_xy, edge_vectors_xy, edge_length_squared)
    end_edge_distances = _point_to_edges_2d(end_xy, edge_starts_xy, edge_vectors_xy, edge_length_squared)
    start_inside = np.all(_edge_signed_distances(start_xy, edge_starts_xy, edge_vectors_xy, edge_lengths_xy) >= 0.0, axis=1)
    end_inside = np.all(_edge_signed_distances(end_xy, edge_starts_xy, edge_vectors_xy, edge_lengths_xy) >= 0.0, axis=1)
    start_distance = np.where(
        start_inside,
        np.abs(d0),
        np.hypot(np.abs(d0), np.min(start_edge_distances, axis=1)),
    )
    end_distance = np.where(
        end_inside,
        np.abs(d1),
        np.hypot(np.abs(d1), np.min(end_edge_distances, axis=1)),
    )
    endpoint_near = (start_distance <= tolerance) | (end_distance <= tolerance)

    crosses_plane = ((d0 < 0.0) & (d1 > 0.0)) | ((d0 > 0.0) & (d1 < 0.0))
    denominator = d0 - d1
    safe_denominator = np.where(crosses_plane, denominator, 1.0)
    along_segment = np.where(crosses_plane, d0 / safe_denominator, 0.0)
    hit_xy = start_xy + along_segment[:, None] * (end_xy - start_xy)
    hit_signed = _edge_signed_distances(hit_xy, edge_starts_xy, edge_vectors_xy, edge_lengths_xy)
    hit_inside = np.all(hit_signed >= 0.0, axis=1)
    hit_edge_distances = _point_to_edges_2d(hit_xy, edge_starts_xy, edge_vectors_xy, edge_length_squared)
    hit_near_edge = np.min(hit_edge_distances, axis=1) <= tolerance

    near_coplanar = (np.abs(d0) <= tolerance) & (np.abs(d1) <= tolerance)
    definitely_blocked = crosses_plane & ~near_coplanar & ~endpoint_near & hit_inside & ~hit_near_edge

    plane_near_at_endpoint = (np.abs(d0) <= tolerance) | (np.abs(d1) <= tolerance)
    projected_near = np.zeros(starts.shape[0], dtype=bool)
    candidates = plane_near_at_endpoint & ~definitely_blocked & ~endpoint_near
    if bool(np.any(candidates)):
        candidate_indices = np.flatnonzero(candidates)
        candidate_d0 = d0[candidate_indices]
        candidate_slope = d1[candidate_indices] - candidate_d0
        candidate_tolerance = tolerance[candidate_indices]
        varying = candidate_slope != 0.0
        safe_slope = np.where(varying, candidate_slope, 1.0)
        interval_a = (-candidate_tolerance - candidate_d0) / safe_slope
        interval_b = (candidate_tolerance - candidate_d0) / safe_slope
        interval_start = np.maximum(0.0, np.minimum(interval_a, interval_b))
        interval_end = np.minimum(1.0, np.maximum(interval_a, interval_b))
        interval_valid = varying & (interval_start <= interval_end)
        constant_near = ~varying & (np.abs(candidate_d0) <= candidate_tolerance)
        interval_start[constant_near] = 0.0
        interval_end[constant_near] = 1.0
        interval_valid |= constant_near
        valid_indices = candidate_indices[interval_valid]
        if valid_indices.size:
            local_start = interval_start[interval_valid]
            local_end = interval_end[interval_valid]
            projected_vector = end_xy[valid_indices] - start_xy[valid_indices]
            near_start_xy = start_xy[valid_indices] + local_start[:, None] * projected_vector
            near_end_xy = start_xy[valid_indices] + local_end[:, None] * projected_vector
            projected_near[valid_indices] = _projected_segment_near_triangle(
                near_start_xy,
                near_end_xy,
                edge_starts_xy,
                edge_vectors_xy,
                edge_lengths_xy,
                edge_length_squared,
                tolerance[valid_indices],
            )

    uncertain = endpoint_near | (crosses_plane & hit_near_edge) | (candidates & projected_near)
    uncertain &= ~definitely_blocked
    return definitely_blocked, uncertain


def _freeze_array(array: np.ndarray) -> np.ndarray:
    frozen_bytes = array.tobytes(order="C")
    return np.frombuffer(frozen_bytes, dtype=array.dtype).reshape(array.shape)


def classify_visibility(
    source_positions: np.ndarray,
    query_positions: np.ndarray,
    triangle_facets: np.ndarray,
    config: VisibilityConfig = VisibilityConfig(),
) -> VisibilityResult:
    """Classify all source-to-query segments against explicit crack triangles.

    `statuses`, `visible_mask`, `blocked_mask`, and `unknown_mask` are
    query-major arrays with shape ``(n_query, n_source)``. A support consumer
    may use a source only where `visible_mask` is true. BLOCKED is returned
    when any strict relative-interior crossing is established, even if a
    different facet is ambiguous; otherwise any ambiguous facet yields
    UNKNOWN. Invalid or degenerate facets are never repaired or ignored.

    The function uses only finite triangle geometry. It accepts no crack IDs,
    tetrahedral adjacency, or direction/side labels. It certifies visibility
    at the exact supplied coordinates only; spatial derivative/local stability
    remains NOT_ESTABLISHED.
    """
    if type(config) is not VisibilityConfig:
        raise ValueError("config must be an exact VisibilityConfig")
    source_input = _preflight_array(
        "source_positions",
        source_positions,
        tail_shape=(3,),
        max_rows=MAX_SOURCES,
        max_bytes=MAX_SOURCES * 3 * 8,
    )
    query_input = _preflight_array(
        "query_positions",
        query_positions,
        tail_shape=(3,),
        max_rows=MAX_QUERIES,
        max_bytes=MAX_QUERIES * 3 * 8,
    )
    triangle_input = _preflight_array(
        "triangle_facets",
        triangle_facets,
        tail_shape=(3, 3),
        max_rows=MAX_TRIANGLES,
        max_bytes=MAX_TRIANGLES * 3 * 3 * 8,
    )

    n_source = int(source_input.shape[0])
    n_query = int(query_input.shape[0])
    n_triangle = int(triangle_input.shape[0])
    work_pairs = n_source * n_query * n_triangle
    if work_pairs > MAX_SEGMENT_TRIANGLE_PAIRS:
        raise ValueError(
            f"budgetRejected: segment-triangle work is {work_pairs}; maximum is {MAX_SEGMENT_TRIANGLE_PAIRS}"
        )

    _validate_finite_values("source_positions", source_input)
    _validate_finite_values("query_positions", query_input)
    _validate_finite_values("triangle_facets", triangle_input)
    maximum_coordinate = 0.0
    for array in (source_input, query_input, triangle_input):
        if array.size:
            maximum_coordinate = max(maximum_coordinate, float(np.max(np.abs(array))))
    if maximum_coordinate > MAX_ABS_COORDINATE_M:
        raise ValueError(f"coordinate magnitude exceeds safe float64 range {MAX_ABS_COORDINATE_M:g} m")

    sources, queries, triangles = _copy_inputs(source_input, query_input, triangle_input)
    n_segment = n_source * n_query
    blocked_any = np.zeros(n_segment, dtype=bool)
    unknown_any = np.zeros(n_segment, dtype=bool)

    if n_segment and n_triangle:
        triangle_frames: list[tuple | None] = []
        triangle_degenerate: list[bool] = []
        triangle_coordinate_scales: list[float] = []
        for triangle in triangles:
            coordinate_scale = float(np.max(np.abs(triangle)))
            edge_ab = triangle[1] - triangle[0]
            edge_ac = triangle[2] - triangle[0]
            normal_cross = np.cross(edge_ab, edge_ac)
            area_twice = float(_norm3(normal_cross))
            edge_bc = triangle[2] - triangle[1]
            edge_lengths = _norm3(np.stack((edge_ab, edge_bc, -edge_ac), axis=0))
            edge_max = float(np.max(edge_lengths))
            triangle_tolerance = max(
                config.absolute_tolerance_m + config.relative_tolerance * edge_max,
                32.0 * _FLOAT64_EPS * coordinate_scale,
            )
            degenerate = (
                edge_max <= triangle_tolerance
                or bool(np.any(edge_lengths <= triangle_tolerance))
                or area_twice <= triangle_tolerance * edge_max
            )
            triangle_degenerate.append(degenerate)
            triangle_coordinate_scales.append(coordinate_scale)
            triangle_frames.append(None if degenerate else _triangle_frame(triangle))

        if any(triangle_degenerate):
            # A malformed facet makes the complete cut set ambiguous. Other
            # established blockers may still reject a segment below.
            unknown_any.fill(True)

        for chunk_start in range(0, n_segment, SEGMENT_CHUNK_SIZE):
            chunk_end = min(chunk_start + SEGMENT_CHUNK_SIZE, n_segment)
            flat_indices = np.arange(chunk_start, chunk_end, dtype=np.int64)
            query_indices = flat_indices // n_source
            source_indices = flat_indices % n_source
            segment_starts = sources[source_indices]
            segment_ends = queries[query_indices]
            segment_coordinate_scale = np.maximum(
                np.max(np.abs(segment_starts), axis=1),
                np.max(np.abs(segment_ends), axis=1),
            )

            for triangle_index, triangle in enumerate(triangles):
                frame = triangle_frames[triangle_index]
                if frame is None:
                    continue
                coordinate_scale = np.maximum(
                    segment_coordinate_scale,
                    triangle_coordinate_scales[triangle_index],
                )
                chunk_blocked, chunk_unknown = _classify_against_triangle(
                    segment_starts,
                    segment_ends,
                    frame,
                    config,
                    coordinate_scale,
                )
                blocked_any[chunk_start:chunk_end] |= chunk_blocked
                unknown_any[chunk_start:chunk_end] |= chunk_unknown

    status_flat = np.full(n_segment, VisibilityStatus.VISIBLE.value, dtype="<U7")
    status_flat[unknown_any] = VisibilityStatus.UNKNOWN.value
    status_flat[blocked_any] = VisibilityStatus.BLOCKED.value
    statuses = status_flat.reshape((n_query, n_source))
    visible_mask = statuses == VisibilityStatus.VISIBLE.value
    blocked_mask = statuses == VisibilityStatus.BLOCKED.value
    unknown_mask = statuses == VisibilityStatus.UNKNOWN.value
    return VisibilityResult(
        statuses=_freeze_array(statuses),
        visible_mask=_freeze_array(visible_mask),
        blocked_mask=_freeze_array(blocked_mask),
        unknown_mask=_freeze_array(unknown_mask),
    )
