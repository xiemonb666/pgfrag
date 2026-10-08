"""Pure-math audit of volume-tetra resolution and directional geometry.

Coordinates and all lengths accepted by this module are in SI metres. Signed
tetrahedron volumes are reported in m^3. The module does not load CAD, call a
solver, alter a mesh, or infer physical adequacy from geometry alone.

Wall thickness, local wall normals, process-zone length, and PD horizon are
independent evidence fields. Caller-provided definitions and source labels
remain declarations; the module does not verify external source bytes. Global
wall qualification therefore remains unknown. A candidate scale screen needs
an explicit named design policy. A result below its limit means only that the
policy did not flag that ratio; it is not a universal sufficiency claim.

Surface tessellation tolerance, surface triangle edge length, and volume-tet
spacing are distinct quantities. This input contains only volume tetrahedra;
no surface tolerance or triangle statistic is inferred from it.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np


SCHEMA = "pgfrag.mesh_resolution_audit.v1"
_EDGE_PAIRS = tuple(itertools.combinations(range(4), 2))
_FACE_TRIPLES = tuple(itertools.combinations(range(4), 3))
_QUANTILES = (0.01, 0.05, 0.50, 0.95, 0.99)
_NUMERICAL_DEGENERACY_RATIO = 1.0e-12


def _safe_number(value: Any) -> float | None:
    """Convert a scalar to a finite Python float; booleans are not numbers here."""
    if isinstance(value, (bool, np.bool_)):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError):
        return None
    return number if math.isfinite(number) else None


def _positive_number(value: Any) -> float | None:
    number = _safe_number(value)
    return number if number is not None and number > 0.0 else None


def _text(value: Any) -> str | None:
    return value.strip() if isinstance(value, str) and value.strip() else None


def _stats(values: Sequence[float]) -> dict[str, Any]:
    clean = np.asarray([float(x) for x in values], dtype=np.float64)
    if not np.isfinite(clean).all():
        raise ValueError("_stats accepts finite values only; callers must report numeric-range failures with their cell indices")
    if clean.size == 0:
        return {
            "count": 0,
            "min": None,
            "p01": None,
            "p05": None,
            "median": None,
            "p95": None,
            "p99": None,
            "max": None,
        }
    scale = float(np.max(np.abs(clean)))
    q = np.zeros(len(_QUANTILES), dtype=np.float64) if scale == 0.0 else np.quantile(clean / scale, _QUANTILES) * scale
    return {
        "count": int(clean.size),
        "min": float(clean.min()),
        "p01": float(q[0]),
        "p05": float(q[1]),
        "median": float(q[2]),
        "p95": float(q[3]),
        "p99": float(q[4]),
        "max": float(clean.max()),
    }


def _index_matrix(cells: Any, vertex_count: int) -> tuple[np.ndarray, list[str | None]]:
    raw = np.asarray(cells, dtype=object)
    if raw.ndim != 2 or raw.shape[1] != 4:
        raise ValueError("cells must be a two-dimensional Mx4 array of zero-based vertex indices")
    indices = np.zeros(raw.shape, dtype=np.int64)
    reasons: list[str | None] = []
    for cell_index, row in enumerate(raw):
        reason = None
        for local_index, item in enumerate(row):
            if isinstance(item, (bool, np.bool_)):
                reason = "INVALID_INDEX"
                continue
            number = _safe_number(item)
            if number is None or not number.is_integer():
                reason = "INVALID_INDEX"
                continue
            index = int(number)
            if index < 0 or index >= vertex_count:
                reason = "INVALID_INDEX"
                continue
            indices[cell_index, local_index] = index
        reasons.append(reason)
    return indices, reasons


def _parse_interval(wall: Mapping[str, Any] | None) -> tuple[tuple[float, float] | None, str]:
    if not isinstance(wall, Mapping):
        return None, "UNKNOWN"
    if wall.get("scope") != "global_interval":
        return None, "UNKNOWN"
    raw = wall.get("thickness_interval_m")
    if not isinstance(raw, (list, tuple)) or len(raw) != 2:
        return None, "INCOMPLETE_CALLER_DECLARATION"
    low = _positive_number(raw[0])
    high = _positive_number(raw[1])
    if low is None or high is None or low > high:
        return None, "INVALID_INTERVAL"
    if _text(wall.get("definition")) is None or _text(wall.get("source")) is None:
        return (low, high), "INCOMPLETE_CALLER_DECLARATION"
    return (low, high), "CALLER_DECLARED_GLOBAL_INTERVAL"


def _parse_policy(policy: Mapping[str, Any] | None) -> tuple[bool, str | None, str | None]:
    if not isinstance(policy, Mapping):
        return False, None, None
    policy_id = _text(policy.get("id"))
    basis = _text(policy.get("basis"))
    return policy_id is not None and basis is not None, policy_id, basis


def _policy_limit(policy: Mapping[str, Any] | None, key: str, policy_valid: bool) -> float | None:
    if not policy_valid or not isinstance(policy, Mapping):
        return None
    return _positive_number(policy.get(key))


def _assessment(ratio: float | None, limit: float | None, *, policy_valid: bool) -> dict[str, Any]:
    if ratio is None:
        return {"status": "UNKNOWN", "ratio": None, "maximum_ratio": limit, "reason": "MEASUREMENT_OR_COVERAGE_MISSING"}
    if not policy_valid:
        return {"status": "UNKNOWN", "ratio": ratio, "maximum_ratio": limit, "reason": "NAMED_DESIGN_POLICY_MISSING"}
    if limit is None:
        return {"status": "UNKNOWN", "ratio": ratio, "maximum_ratio": None, "reason": "POLICY_LIMIT_MISSING_OR_INVALID"}
    if ratio > limit:
        return {"status": "RESOLUTION_INSUFFICIENT", "ratio": ratio, "maximum_ratio": limit, "reason": "RATIO_EXCEEDS_EXPLICIT_POLICY_LIMIT"}
    return {
        "status": "NOT_FLAGGED_BY_DESIGN_POLICY",
        "ratio": ratio,
        "maximum_ratio": limit,
        "reason": "BELOW_POLICY_LIMIT_ONLY_NOT_A_PHYSICAL_SUFFICIENCY_CLAIM",
    }


def _safe_divide(numerator: float, denominator: float) -> tuple[float | None, str | None]:
    """Divide finite values without allowing JSON inf/nan or silent zero underflow."""
    if not math.isfinite(numerator) or not math.isfinite(denominator) or numerator < 0.0 or denominator <= 0.0:
        return None, "INVALID_RATIO_INPUT"
    try:
        value = numerator / denominator
    except (OverflowError, ZeroDivisionError):
        return None, "OVERFLOW"
    if not math.isfinite(value):
        return None, "OVERFLOW"
    if numerator > 0.0 and value == 0.0:
        return None, "UNDERFLOW"
    return value, None


def _cellwise_ratios(
    metric: str,
    numerators_and_denominators: Sequence[tuple[int, float, float]],
    issues: list[dict[str, Any]],
) -> tuple[list[float], list[int]]:
    values: list[float] = []
    failed: list[tuple[int, str]] = []
    for cell_index, numerator, denominator in numerators_and_denominators:
        value, failure = _safe_divide(float(numerator), float(denominator))
        if failure is None and value is not None:
            values.append(value)
        else:
            failed.append((int(cell_index), failure or "UNKNOWN"))
    if failed:
        issues.append({
            "metric": metric,
            "count": int(len(failed)),
            "cell_indices": [cell_index for cell_index, _ in failed],
            "range_statuses": [status for _, status in failed],
        })
    return values, [cell_index for cell_index, _ in failed]


def _scaled_signed_volume(determinant: float, scale: float) -> tuple[float | None, str]:
    """Return signed SI volume, detecting range loss from dimensionless geometry."""
    if determinant == 0.0:
        return 0.0, "FINITE"
    log_magnitude = math.log(abs(determinant)) + 3.0 * math.log(scale) - math.log(6.0)
    if log_magnitude > math.log(float.fromhex("0x1.fffffffffffffp+1023")):
        return None, "OVERFLOW"
    try:
        magnitude = math.exp(log_magnitude)
    except OverflowError:
        return None, "OVERFLOW"
    if not math.isfinite(magnitude):
        return None, "OVERFLOW"
    if magnitude == 0.0:
        return None, "UNDERFLOW"
    return math.copysign(magnitude, determinant), "FINITE"


def _signed_volume_log10_abs(determinant: float, scale: float) -> float | None:
    """Keep volume magnitude available when SI m^3 is outside float range."""
    if determinant == 0.0 or scale <= 0.0:
        return None
    log_magnitude = math.log(abs(determinant)) + 3.0 * math.log(scale) - math.log(6.0)
    return log_magnitude / math.log(10.0)


def _checked_rows(values: Any, expected: int, width: int) -> tuple[list[np.ndarray | None], int, str | None]:
    if values is None:
        return [None] * expected, 0, "MISSING"
    try:
        raw = np.asarray(values, dtype=object)
    except (TypeError, ValueError):
        return [None] * expected, 0, "INVALID_SHAPE"
    if raw.shape != (expected, width):
        return [None] * expected, 0, "INVALID_SHAPE"
    rows: list[np.ndarray | None] = []
    good = 0
    for row in raw:
        parsed = [_safe_number(item) for item in row]
        if any(item is None for item in parsed):
            rows.append(None)
            continue
        vector = np.asarray(parsed, dtype=np.float64)
        norm = math.hypot(float(vector[0]), float(vector[1]), float(vector[2])) if width == 3 else None
        if width == 3:
            if norm is None or norm <= 0.0 or not math.isfinite(norm):
                rows.append(None)
                continue
            vector = vector / norm
        rows.append(vector)
        good += 1
    return rows, good, None


def audit_tetra_mesh(
    vertices_m: Any,
    cells: Any,
    *,
    wall_measurement: Mapping[str, Any] | None = None,
    wall_normal_map: Mapping[str, Any] | None = None,
    process_zone: Mapping[str, Any] | None = None,
    pd_horizon_m: Any = None,
    resolution_policy: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Audit tetra coordinates/connectivity and independently documented scales.

    ``vertices_m`` and ``cells`` are SI coordinates and zero-based tetra
    connectivity. ``wall_measurement`` may describe a caller-declared global
    thickness interval or a local cell map. ``wall_normal_map`` is a local
    cell-ordered vector map; without complete provenance it remains a numeric
    diagnostic and cannot support a wall-resolution classification.

    A policy may contain these independent, dimensionless limits:
    ``max_h_to_wall``, ``max_normal_projected_edge_to_wall``,
    ``max_h_to_process_zone``, and ``max_h_to_pd_horizon``. Every classification
    also requires non-empty policy ``id`` and ``basis``. There are no default
    physical sufficiency limits.
    """
    try:
        vertices = np.asarray(vertices_m, dtype=np.float64)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("vertices_m must be a numeric Nx3 array") from exc
    if vertices.ndim != 2 or vertices.shape[1] != 3:
        raise ValueError("vertices_m must be a numeric Nx3 array")

    indices, index_reasons = _index_matrix(cells, int(vertices.shape[0]))
    cell_count = int(indices.shape[0])
    finite_vertices = np.isfinite(vertices).all(axis=1)
    per_cell: list[dict[str, Any]] = []
    internal: list[dict[str, Any] | None] = []
    edge_lengths_all: list[float] = []
    cell_max_edges: list[float] = []
    cell_min_heights: list[float] = []
    signed_volumes: list[float] = []
    status_counts: Counter[str] = Counter()
    orientation_counts: Counter[str] = Counter()
    near_degenerate_count = 0
    zero_volume_count = 0
    numeric_range_issues: list[dict[str, Any]] = []
    invalid_index_count = sum(reason == "INVALID_INDEX" for reason in index_reasons)

    for cell_index in range(cell_count):
        index_reason = index_reasons[cell_index]
        if index_reason is not None:
            row = {"cell_index": cell_index, "status": index_reason, "orientation": "UNKNOWN", "signed_volume_m3": None, "minimum_height_m": None, "maximum_edge_m": None}
            per_cell.append(row)
            internal.append(None)
            status_counts[index_reason] += 1
            orientation_counts["unknown"] += 1
            continue

        tet = vertices[indices[cell_index]]
        if not np.isfinite(tet).all():
            row = {"cell_index": cell_index, "status": "NONFINITE_COORDINATE", "orientation": "UNKNOWN", "signed_volume_m3": None, "minimum_height_m": None, "maximum_edge_m": None}
            per_cell.append(row)
            internal.append(None)
            status_counts["NONFINITE_COORDINATE"] += 1
            orientation_counts["unknown"] += 1
            continue

        p0 = tet[0]
        edge_vectors = [tet[j] - tet[i] for i, j in _EDGE_PAIRS]
        raw_edge_lengths = [math.hypot(float(v[0]), float(v[1]), float(v[2])) for v in edge_vectors]
        if not all(math.isfinite(length) for length in raw_edge_lengths):
            row = {"cell_index": cell_index, "status": "NUMERIC_UNSTABLE", "orientation": "UNKNOWN", "signed_volume_m3": None, "minimum_height_m": None, "maximum_edge_m": None}
            per_cell.append(row)
            internal.append(None)
            status_counts["NUMERIC_UNSTABLE"] += 1
            orientation_counts["unknown"] += 1
            continue
        scale = max(raw_edge_lengths)
        if scale <= 0.0:
            row = {"cell_index": cell_index, "status": "ZERO_SIZE", "orientation": "ZERO", "signed_volume_m3": 0.0, "minimum_height_m": 0.0, "maximum_edge_m": 0.0}
            per_cell.append(row)
            internal.append({"edge_vectors": edge_vectors, "edge_lengths": raw_edge_lengths, "max_edge": 0.0, "min_height": 0.0, "volume": 0.0, "orientation": "ZERO"})
            status_counts["ZERO_SIZE"] += 1
            orientation_counts["zero"] += 1
            zero_volume_count += 1
            signed_volumes.append(0.0)
            continue

        normalized = [(tet[i] - p0) / scale for i in range(4)]
        try:
            determinant = float(np.dot(np.cross(normalized[1], normalized[2]), normalized[3]))
        except (FloatingPointError, ValueError):
            determinant = math.nan
        if not math.isfinite(determinant):
            row = {"cell_index": cell_index, "status": "NUMERIC_UNSTABLE", "orientation": "UNKNOWN", "signed_volume_m3": None, "minimum_height_m": None, "maximum_edge_m": float(scale)}
            per_cell.append(row)
            internal.append(None)
            status_counts["NUMERIC_UNSTABLE"] += 1
            orientation_counts["unknown"] += 1
            continue

        if determinant > 0.0:
            orientation = "POSITIVE"
            orientation_counts["positive"] += 1
        elif determinant < 0.0:
            orientation = "NEGATIVE"
            orientation_counts["negative"] += 1
        else:
            orientation = "ZERO"
            orientation_counts["zero"] += 1

        normalized_signed_volume = determinant / 6.0
        signed_volume, volume_range_status = _scaled_signed_volume(determinant, scale)
        signed_volume_log10_abs = _signed_volume_log10_abs(determinant, scale)

        if signed_volume is not None:
            signed_volumes.append(signed_volume)

        heights: list[float] = []
        for face in _FACE_TRIPLES:
            u = normalized[face[1]] - normalized[face[0]]
            v = normalized[face[2]] - normalized[face[0]]
            twice_area_scaled = math.hypot(*map(float, np.cross(u, v)))
            height = 0.0 if twice_area_scaled == 0.0 else abs(determinant) * scale / twice_area_scaled
            heights.append(height if math.isfinite(height) else math.nan)
        minimum_height = min((x for x in heights if math.isfinite(x)), default=math.nan)
        if math.isfinite(minimum_height):
            cell_min_heights.append(minimum_height)
        max_edge = float(scale)
        edge_lengths_all.extend(raw_edge_lengths)
        cell_max_edges.append(max_edge)
        if determinant == 0.0:
            zero_volume_count += 1
            status = "ZERO_VOLUME"
        elif volume_range_status != "FINITE":
            status = f"UNKNOWN_NUMERIC_RANGE/{volume_range_status}"
            numeric_range_issues.append({
                "cell_index": int(cell_index),
                "status": volume_range_status,
                "maximum_edge_scale_m": float(scale),
                "normalized_signed_volume_over_h_cubed": float(normalized_signed_volume),
                "signed_volume_log10_abs_m3": float(signed_volume_log10_abs) if signed_volume_log10_abs is not None else None,
                "signed_volume_m3": None,
            })
        elif minimum_height / scale <= _NUMERICAL_DEGENERACY_RATIO:
            near_degenerate_count += 1
            status = "NEAR_DEGENERATE"
        elif orientation == "NEGATIVE":
            status = "NEGATIVE_ORIENTATION"
        else:
            status = "OK"
        status_counts[status] += 1
        row = {
            "cell_index": cell_index,
            "status": status,
            "orientation": orientation,
            "signed_volume_m3": signed_volume,
            "maximum_edge_scale_m": float(scale),
            "normalized_signed_volume_over_h_cubed": float(normalized_signed_volume),
            "signed_volume_log10_abs_m3": float(signed_volume_log10_abs) if signed_volume_log10_abs is not None else None,
            "volume_computation_status": volume_range_status,
            "minimum_height_m": float(minimum_height) if math.isfinite(minimum_height) else None,
            "minimum_height_over_h": float(minimum_height / scale) if math.isfinite(minimum_height) else None,
            "maximum_edge_m": max_edge,
        }
        per_cell.append(row)
        internal.append({
            "edge_vectors": edge_vectors,
            "edge_lengths": raw_edge_lengths,
            "max_edge": max_edge,
            "min_height": minimum_height,
            "volume": signed_volume,
            "orientation": orientation,
        })

    h_p95 = _stats(cell_max_edges)["p95"]
    h_median = _stats(cell_max_edges)["median"]
    valid_cell_lengths = [(i, float(data["max_edge"])) for i, data in enumerate(internal) if data is not None]
    orientation_summary = {
        "positive": int(orientation_counts["positive"]),
        "negative": int(orientation_counts["negative"]),
        "zero": int(orientation_counts["zero"]),
        "unknown": int(orientation_counts["unknown"]),
    }
    signed_volume_stats = _stats(signed_volumes)
    signed_volume_stats["unknown_numeric_range_count"] = int(len(numeric_range_issues))
    signed_volume_stats["unknown_numeric_range_cell_indices"] = [int(item["cell_index"]) for item in numeric_range_issues]

    # The unoriented second moment describes mesh edge directions only. It has
    # no crack-path, fracture-law, or mesh-independence interpretation.
    direction_vectors: list[np.ndarray] = []
    for data in internal:
        if data is None:
            continue
        for vector in data["edge_vectors"]:
            length = math.hypot(float(vector[0]), float(vector[1]), float(vector[2]))
            if length > 0.0 and math.isfinite(length):
                direction_vectors.append(np.asarray(vector, dtype=np.float64) / length)
    if direction_vectors:
        directions = np.asarray(direction_vectors, dtype=np.float64)
        tensor = np.einsum("ni,nj->ij", directions, directions) / float(len(directions))
        eigenvalues = np.linalg.eigvalsh(tensor)
        direction_tensor: dict[str, Any] = {
            "status": "DIAGNOSTIC_ONLY",
            "edge_occurrence_count": int(len(directions)),
            "unoriented_unit_edge_second_moment": tensor.tolist(),
            "eigenvalues_ascending": eigenvalues.tolist(),
            "eigenvalue_max_over_min": float(eigenvalues[-1] / max(float(eigenvalues[0]), 1.0e-15)),
        }
    else:
        direction_tensor = {
            "status": "UNKNOWN",
            "edge_occurrence_count": 0,
            "unoriented_unit_edge_second_moment": None,
            "eigenvalues_ascending": None,
            "eigenvalue_max_over_min": None,
        }

    policy_valid, policy_id, policy_basis = _parse_policy(resolution_policy)
    wall_interval, wall_status = _parse_interval(wall_measurement)
    wall_caller_declared = wall_status == "CALLER_DECLARED_GLOBAL_INTERVAL"
    wall_source_verification_status = "UNKNOWN_NOT_PERFORMED_BY_PURE_MATH_AUDIT"
    valid_geometry_count = sum(data is not None for data in internal)
    ratio_issues: list[dict[str, Any]] = []
    local_thickness: list[float | None] = [None] * cell_count
    local_thickness_status = "MISSING"
    if isinstance(wall_measurement, Mapping) and wall_measurement.get("scope") == "local_map":
        raw_thickness = wall_measurement.get("thickness_by_cell_m")
        if raw_thickness is None:
            local_thickness_status = "MISSING"
        elif not isinstance(raw_thickness, (list, tuple)) or len(raw_thickness) != cell_count:
            local_thickness_status = "INVALID_SHAPE"
        else:
            local_thickness = [_positive_number(item) for item in raw_thickness]
            if _text(wall_measurement.get("definition")) is not None and _text(wall_measurement.get("source")) is not None:
                local_thickness_status = "CALLER_DECLARED_LOCAL_MAP"
            else:
                local_thickness_status = "INCOMPLETE_CALLER_DECLARATION"

    normal_rows, normal_count, normal_shape_status = _checked_rows(
        wall_normal_map.get("normals_by_cell") if isinstance(wall_normal_map, Mapping) else None,
        cell_count,
        3,
    )
    normal_metadata_complete = (
        isinstance(wall_normal_map, Mapping)
        and _text(wall_normal_map.get("definition")) is not None
        and _text(wall_normal_map.get("source")) is not None
    )
    mapped_valid_cell_ids = [i for i, data in enumerate(internal) if data is not None and normal_rows[i] is not None]
    normal_coverage = (len(mapped_valid_cell_ids) / valid_geometry_count) if valid_geometry_count else 0.0
    normal_evidence_complete = bool(
        normal_metadata_complete
        and valid_geometry_count > 0
        and len(mapped_valid_cell_ids) == valid_geometry_count
    )

    normal_projection_by_cell: list[float | None] = [None] * cell_count
    normal_ratio_candidates: list[tuple[int, float, float]] = []
    for cell_index in mapped_valid_cell_ids:
        n = normal_rows[cell_index]
        data = internal[cell_index]
        assert n is not None and data is not None
        projections = [abs(float(np.dot(vector, n))) for vector in data["edge_vectors"]]
        projection = max(projections, default=0.0)
        if math.isfinite(projection):
            normal_projection_by_cell[cell_index] = projection
        thickness = local_thickness[cell_index] if local_thickness[cell_index] is not None else (wall_interval[0] if wall_interval else None)
        if thickness is not None and thickness > 0.0 and math.isfinite(projection):
            normal_ratio_candidates.append((cell_index, projection, thickness))

    projected_values = [x for x in normal_projection_by_cell if x is not None]
    normal_ratio_values, normal_ratio_failed_ids = _cellwise_ratios(
        "normal_projected_edge_to_wall",
        normal_ratio_candidates,
        ratio_issues,
    )
    normal_projection_status = (
        "CALLER_DECLARED_LOCAL_NORMAL_MAP" if normal_evidence_complete
        else "INCOMPLETE_CALLER_DECLARATION" if normal_count or normal_shape_status == "INVALID_SHAPE"
        else "UNKNOWN"
    )
    normal_ratio_covered = (
        normal_evidence_complete
        and len(normal_ratio_candidates) == valid_geometry_count
        and (wall_interval is not None or local_thickness_status == "CALLER_DECLARED_LOCAL_MAP")
    )
    normal_ratio_p95 = _stats(normal_ratio_values)["p95"] if normal_ratio_covered and not normal_ratio_failed_ids else None
    normal_limit = _policy_limit(resolution_policy, "max_normal_projected_edge_to_wall", policy_valid)
    normal_assessment = _assessment(normal_ratio_p95, normal_limit, policy_valid=policy_valid)
    if not normal_evidence_complete:
        normal_assessment = {
            "status": "UNKNOWN",
            "ratio": normal_ratio_p95,
            "maximum_ratio": normal_limit,
            "reason": "COMPLETE_CALLER_DECLARED_LOCAL_WALL_NORMAL_MAPPING_REQUIRED",
        }
    elif not normal_ratio_covered or normal_ratio_failed_ids:
        normal_assessment = {
            "status": "UNKNOWN",
            "ratio": None,
            "maximum_ratio": normal_limit,
            "reason": "LOCAL_THICKNESS_OR_GLOBAL_INTERVAL_COVERAGE_REQUIRED" if not normal_ratio_failed_ids else "NONFINITE_RATIO_NUMERIC_RANGE",
        }

    global_wall_ratio_candidates = [
        (cell_index, edge_length, wall_interval[0])
        for cell_index, edge_length in valid_cell_lengths
    ] if wall_interval is not None else []
    global_wall_ratio_values, global_wall_ratio_failed_ids = _cellwise_ratios(
        "h_to_wall", global_wall_ratio_candidates, ratio_issues
    )
    h_to_wall = (
        _stats(global_wall_ratio_values)["p95"]
        if wall_interval is not None
        and len(global_wall_ratio_candidates) == valid_geometry_count
        and not global_wall_ratio_failed_ids
        else None
    )
    wall_limit = _policy_limit(resolution_policy, "max_h_to_wall", policy_valid)
    wall_candidate_assessment = _assessment(h_to_wall, wall_limit, policy_valid=policy_valid)
    wall_assessment = {
        "status": "UNKNOWN",
        "ratio": h_to_wall,
        "maximum_ratio": wall_limit,
        "reason": "ACTUAL_GLOBAL_WALL_SOURCE_VERIFICATION_NOT_PERFORMED_BY_THIS_MODULE",
    }

    local_h_wall_candidates = [
        (i, float(internal[i]["max_edge"]), float(local_thickness[i]))
        for i in range(cell_count)
        if internal[i] is not None and local_thickness[i] is not None and local_thickness_status == "CALLER_DECLARED_LOCAL_MAP"
    ]
    local_h_wall_ratios, local_h_wall_failed_ids = _cellwise_ratios(
        "local_cell_h_to_wall", local_h_wall_candidates, ratio_issues
    )

    process_length = _positive_number(process_zone.get("length_m")) if isinstance(process_zone, Mapping) else None
    process_caller_declared = bool(
        process_length is not None
        and isinstance(process_zone, Mapping)
        and _text(process_zone.get("definition")) is not None
        and _text(process_zone.get("source")) is not None
    )
    process_ratio_candidates = [(i, length, process_length) for i, length in valid_cell_lengths] if process_length is not None else []
    process_ratio_values, process_ratio_failed_ids = _cellwise_ratios(
        "h_to_process_zone", process_ratio_candidates, ratio_issues
    )
    h_to_process = (
        _stats(process_ratio_values)["p95"]
        if process_length is not None
        and len(process_ratio_candidates) == valid_geometry_count
        and not process_ratio_failed_ids
        else None
    )
    process_limit = _policy_limit(resolution_policy, "max_h_to_process_zone", policy_valid)
    process_assessment = _assessment(h_to_process, process_limit, policy_valid=policy_valid)
    if not process_caller_declared:
        process_assessment = {"status": "UNKNOWN", "ratio": h_to_process, "maximum_ratio": process_limit, "reason": "PROCESS_ZONE_LENGTH_DEFINITION_AND_SOURCE_REQUIRED"}
    elif process_ratio_failed_ids:
        process_assessment = {"status": "UNKNOWN", "ratio": None, "maximum_ratio": process_limit, "reason": "NONFINITE_RATIO_NUMERIC_RANGE"}

    horizon = _positive_number(pd_horizon_m)
    h_delta_ratio_candidates = [(i, length, horizon) for i, length in valid_cell_lengths] if horizon is not None else []
    h_delta_ratio_values, h_delta_ratio_failed_ids = _cellwise_ratios("h_over_pd_horizon", h_delta_ratio_candidates, ratio_issues)
    h_over_delta = (
        _stats(h_delta_ratio_values)["p95"]
        if horizon is not None
        and len(h_delta_ratio_candidates) == valid_geometry_count
        and not h_delta_ratio_failed_ids
        else None
    )
    delta_over_h: float | None = None
    delta_over_h_failure: str | None = None
    if horizon is not None and h_p95 is not None and h_p95 > 0.0:
        delta_over_h, delta_over_h_failure = _safe_divide(horizon, float(h_p95))
    elif horizon is not None and h_p95 == 0.0:
        delta_over_h_failure = "ZERO_DENOMINATOR"
    if delta_over_h_failure is not None:
        ratio_issues.append({"metric": "delta_over_h", "count": 1, "cell_indices": [i for i, _ in valid_cell_lengths], "range_statuses": [delta_over_h_failure], "aggregation": "scalar_h_is_p95_per_cell_maximum_edge"})
    delta_over_wall: float | None = None
    delta_over_wall_failure: str | None = None
    if horizon is not None and wall_interval is not None:
        delta_over_wall, delta_over_wall_failure = _safe_divide(horizon, wall_interval[0])
    if delta_over_wall_failure is not None:
        ratio_issues.append({"metric": "delta_over_wall_lower_bound", "count": 1, "cell_indices": [], "range_statuses": [delta_over_wall_failure], "aggregation": "scalar_source_declared_wall_interval"})
    horizon_limit = _policy_limit(resolution_policy, "max_h_to_pd_horizon", policy_valid)
    horizon_assessment = _assessment(h_over_delta, horizon_limit, policy_valid=policy_valid)
    if horizon is None:
        horizon_assessment = {"status": "UNKNOWN", "ratio": None, "maximum_ratio": horizon_limit, "reason": "PD_HORIZON_MISSING_OR_INVALID"}
    elif h_delta_ratio_failed_ids:
        horizon_assessment = {"status": "UNKNOWN", "ratio": None, "maximum_ratio": horizon_limit, "reason": "NONFINITE_RATIO_NUMERIC_RANGE"}

    wall_report: dict[str, Any] = {
        "evidence_status": wall_status,
        "scope": "global_interval" if isinstance(wall_measurement, Mapping) and wall_measurement.get("scope") == "global_interval" else "local_or_missing",
        "thickness_interval_m": [float(wall_interval[0]), float(wall_interval[1])] if wall_interval is not None else None,
        "caller_declaration": {
            "definition": _text(wall_measurement.get("definition")) if isinstance(wall_measurement, Mapping) else None,
            "source": _text(wall_measurement.get("source")) if isinstance(wall_measurement, Mapping) else None,
            "status": wall_status,
        },
        "actual_source_verification": {
            "status": wall_source_verification_status,
            "verified": False,
            "reason": "The pure-math audit does not open or independently verify source bytes or global coverage.",
        },
        "global_wall_status": "UNKNOWN",
        "global_wall_qualification_status": "UNKNOWN",
        "h_to_wall": {
            "definition": "p95 of per-cell maximum edge divided by the lower bound of the caller-declared global interval; candidate ratio only",
            "value": h_to_wall,
            "value_status": (
                "UNKNOWN_NUMERIC_RANGE" if global_wall_ratio_failed_ids
                else "CALLER_DECLARED_CANDIDATE" if h_to_wall is not None
                else "UNKNOWN"
            ),
            "policy_assessment": wall_assessment,
            "caller_declared_candidate_policy_screen": wall_candidate_assessment if wall_caller_declared else {
                "status": "UNKNOWN",
                "ratio": h_to_wall,
                "maximum_ratio": wall_limit,
                "reason": "COMPLETE_CALLER_DECLARATION_REQUIRED_FOR_CANDIDATE_SCREEN",
            },
            "numeric_range_failure_count": len(global_wall_ratio_failed_ids),
            "numeric_range_failure_cell_indices": global_wall_ratio_failed_ids,
        },
    }

    if local_thickness_status == "CALLER_DECLARED_LOCAL_MAP":
        local_h_wall_p95 = _stats(local_h_wall_ratios)["p95"]
        local_wall_limit = _policy_limit(resolution_policy, "max_h_to_wall", policy_valid)
        wall_report["local_cell_h_to_wall"] = {
            "definition": "per-cell maximum edge divided by caller-declared local thickness; candidate diagnostic only",
            "coverage_count": int(len(local_h_wall_ratios)),
            "coverage_fraction_of_valid_cells": (len(local_h_wall_ratios) / valid_geometry_count) if valid_geometry_count else 0.0,
            "quantiles": _stats(local_h_wall_ratios),
            "policy_assessment": _assessment(
                local_h_wall_p95 if len(local_h_wall_ratios) == valid_geometry_count and valid_geometry_count > 0 and not local_h_wall_failed_ids else None,
                local_wall_limit,
                policy_valid=policy_valid,
            ),
            "numeric_range_failure_count": len(local_h_wall_failed_ids),
            "numeric_range_failure_cell_indices": local_h_wall_failed_ids,
        }
    wall_report["local_normal_projection"] = {
        "evidence_status": normal_projection_status,
        "definition": "per-cell maximum absolute projection of the six tetra edge vectors onto the supplied unit wall normal; no maximum-edge substitute is used",
        "normal_map_coverage_fraction_of_valid_cells": float(normal_coverage),
        "normal_map_coverage_count": int(len(mapped_valid_cell_ids)),
        "projected_edge_m_quantiles": _stats(projected_values),
        "projected_edge_over_local_wall_quantiles": _stats(normal_ratio_values),
        "policy_assessment": normal_assessment,
        "numeric_range_failure_count": len(normal_ratio_failed_ids),
        "numeric_range_failure_cell_indices": normal_ratio_failed_ids,
        "definition_source": _text(wall_normal_map.get("definition")) if isinstance(wall_normal_map, Mapping) else None,
        "mapping_source": _text(wall_normal_map.get("source")) if isinstance(wall_normal_map, Mapping) else None,
    }

    process_report = {
        "evidence_status": "CALLER_DECLARED" if process_caller_declared else "UNKNOWN",
        "length_m": process_length if process_caller_declared else None,
        "definition": _text(process_zone.get("definition")) if isinstance(process_zone, Mapping) else None,
        "source": _text(process_zone.get("source")) if isinstance(process_zone, Mapping) else None,
        "actual_source_verification_status": "UNKNOWN_NOT_PERFORMED_BY_PURE_MATH_AUDIT",
        "h_to_process_zone": {
            "definition": "p95 of per-cell maximum edge divided by the separately caller-declared process-zone length",
            "value": h_to_process,
            "policy_assessment": process_assessment,
            "numeric_range_failure_count": len(process_ratio_failed_ids),
            "numeric_range_failure_cell_indices": process_ratio_failed_ids,
        },
    }

    pd_report = {
        "horizon_m": horizon,
        "status": "INPUT_REPORTED" if horizon is not None else "UNKNOWN",
        "delta_over_h": delta_over_h,
        "delta_over_h_status": "UNKNOWN_NUMERIC_RANGE" if delta_over_h_failure is not None else "COMPUTED" if delta_over_h is not None else "UNKNOWN",
        "delta_over_h_numeric_range_failure_count": int(delta_over_h_failure is not None),
        "delta_over_h_numeric_range_failure_cell_indices": [i for i, _ in valid_cell_lengths] if delta_over_h_failure is not None else [],
        "h_over_delta": h_over_delta,
        "delta_over_wall_lower_bound": delta_over_wall,
        "delta_over_wall_status": (
            "UNKNOWN_NUMERIC_RANGE" if delta_over_wall_failure is not None
            else "CALLER_DECLARED_CANDIDATE" if delta_over_wall is not None
            else "UNKNOWN"
        ),
        "delta_over_wall_numeric_range_failure_count": int(delta_over_wall_failure is not None),
        "h_over_delta_policy_assessment": horizon_assessment,
        "h_over_delta_numeric_range_failure_count": len(h_delta_ratio_failed_ids),
        "h_over_delta_numeric_range_failure_cell_indices": h_delta_ratio_failed_ids,
        "independence_note": "PD horizon delta is reported only from pd_horizon_m; it is never used as a process-zone length.",
    }

    if orientation_counts["positive"] + orientation_counts["negative"] + orientation_counts["zero"] == 0:
        all_positive: bool | None = None
    else:
        all_positive = bool(
            orientation_counts["positive"] == cell_count
            and orientation_counts["negative"] == 0
            and orientation_counts["zero"] == 0
            and orientation_counts["unknown"] == 0
        )

    return {
        "schema": SCHEMA,
        "scope": "pure mathematical tetra geometry and declared scale audit; no CAD, solver, fracture, or physical validation",
        "units": {"coordinates": "m", "lengths": "m", "signed_tetrahedron_volume": "m^3", "ratios": "dimensionless"},
        "input_integrity": {
            "vertex_count": int(vertices.shape[0]),
            "tetrahedron_count": cell_count,
            "nonfinite_vertex_count": int(np.count_nonzero(~finite_vertices)),
            "invalid_index_cell_count": int(invalid_index_count),
            "nonfinite_coordinate_cell_count": int(status_counts["NONFINITE_COORDINATE"]),
            "numeric_unstable_cell_count": int(status_counts["NUMERIC_UNSTABLE"]),
            "unknown_numeric_range_cell_count": int(len(numeric_range_issues)),
            "unknown_numeric_range_cell_indices": [int(item["cell_index"]) for item in numeric_range_issues],
            "has_invalid_indices": bool(invalid_index_count > 0),
            "has_nonfinite_vertices": bool(np.any(~finite_vertices)),
        },
        "tetra_geometry": {
            "signed_volume_m3": signed_volume_stats,
            "signed_orientation_counts": orientation_summary,
            "all_signed_orientations_positive": all_positive,
            "negative_orientation_is_reported_separately_from_degeneracy": True,
            "zero_volume_count": int(zero_volume_count),
            "near_degenerate_count": int(near_degenerate_count),
            "unknown_numeric_range_count": int(len(numeric_range_issues)),
            "unknown_numeric_range_issues": numeric_range_issues,
            "all_nonzero_signed_volumes_representable": bool(
                len(numeric_range_issues) == 0 and orientation_counts["unknown"] == 0
            ),
            "numerical_degeneracy_ratio_threshold": _NUMERICAL_DEGENERACY_RATIO,
            "pooled_tetra_edge_length_m": _stats(edge_lengths_all),
            "per_tetra_maximum_edge_m": _stats(cell_max_edges),
            "minimum_tetra_height_m": _stats(cell_min_heights),
            "mesh_characteristic_h_m": {
                "value": float(h_p95) if h_p95 is not None else None,
                "definition": "p95 of per-tetra maximum edge length; used only as a declared geometric comparison scale",
                "median_per_tetra_maximum_edge_m": float(h_median) if h_median is not None else None,
            },
            "cell_status_counts": {key: int(value) for key, value in sorted(status_counts.items())},
            "per_cell": per_cell,
        },
        "surface_vs_volume_discretization": {
            "surface_tessellation_tolerance_m": None,
            "surface_tessellation_tolerance_status": "NOT_PROVIDED",
            "surface_triangle_edge_length_m": _stats([]),
            "surface_triangle_edge_length_status": "NOT_PROVIDED",
            "volume_tetra_spacing_m": float(h_p95) if h_p95 is not None else None,
            "volume_tetra_spacing_definition": "p95 of per-tetra maximum edge lengths from supplied 3D tetrahedra",
            "interpretation_limit": "A surface tessellation tolerance such as 0.5 mm, or a surface triangle edge size, is not volume-tet spacing and cannot establish volume-mesh resolution.",
        },
        "resolution_policy": {
            "status": "VALID_NAMED_POLICY" if policy_valid else "UNKNOWN",
            "id": policy_id,
            "basis": policy_basis,
            "limits_are_design_inputs_not_universal_physical_constants": True,
        },
        "wall_resolution": wall_report,
        "process_zone_resolution": process_report,
        "pd_horizon": pd_report,
        "numeric_ratio_issues": ratio_issues,
        "directional_diagnostics": {
            "tetra_edge_direction_tensor": direction_tensor,
            "single_mesh_crack_direction_bias_status": "UNKNOWN",
            "independent_layout_and_native_case_status": "UNKNOWN",
            "interpretation_limit": "A single mesh edge tensor or edge-direction distribution cannot establish unbiased crack direction. Independent layouts and native-case fracture results require separate validation.",
        },
    }


def _load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as stream:
        return json.load(stream)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-json", type=Path, required=True, help="JSON containing vertices_m and cells")
    parser.add_argument("--output-json", type=Path, required=True, help="path for the audit report")
    args = parser.parse_args(argv)
    try:
        source = _load_json(args.input_json)
        report = audit_tetra_mesh(
            source["vertices_m"],
            source["cells"],
            wall_measurement=source.get("wall_measurement"),
            wall_normal_map=source.get("wall_normal_map"),
            process_zone=source.get("process_zone"),
            pd_horizon_m=source.get("pd_horizon_m"),
            resolution_policy=source.get("resolution_policy"),
        )
        rendered = json.dumps(report, indent=2, allow_nan=False)
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(rendered + "\n", encoding="utf-8")
        print(json.dumps({
            "status": "MATH_AUDIT_COMPLETE",
            "output_json": str(args.output_json),
            "tetrahedron_count": report["tetra_geometry"]["signed_orientation_counts"],
            "wall_status": report["wall_resolution"]["evidence_status"],
            "h_to_wall": report["wall_resolution"]["h_to_wall"]["value"],
            "normal_status": report["wall_resolution"]["local_normal_projection"]["evidence_status"],
            "process_zone_status": report["process_zone_resolution"]["evidence_status"],
            "direction_bias_status": report["directional_diagnostics"]["single_mesh_crack_direction_bias_status"],
        }, indent=2))
        return 0
    except (KeyError, OSError, TypeError, ValueError, json.JSONDecodeError) as exc:
        print(json.dumps({"status": "INPUT_ERROR", "error": str(exc)}), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
