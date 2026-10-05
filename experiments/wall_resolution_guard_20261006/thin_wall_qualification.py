"""Local thin-wall resolution qualification for analytic annular-cell fixtures.

This isolated review utility has no CAD or solver imports.  It accepts unit-tagged
wall samples, element geometry/max-edge data, centroid volume weights, and segment
provenance.  Its decision is deliberately local: a mesh-quality precondition is
not a physical qualification or a full-domain guarantee.
"""

from __future__ import annotations

import math
from collections import defaultdict
from typing import Any, Iterable


SCHEMA = "pgfrag.thin_wall_resolution_input.v1"
MAX_CELLS = 96
MAX_SEGMENTS = 128
LENGTH_FACTORS_M = {"m": 1.0, "mm": 1.0e-3}
VOLUME_FACTORS_M3 = {"m^3": 1.0, "m3": 1.0, "mm^3": 1.0e-9, "mm3": 1.0e-9}


class QualificationInputError(ValueError):
    """Raised for malformed, dimensionally invalid, or over-budget input."""


def _finite_number(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise QualificationInputError(f"{label} must be a finite number")
    result = float(value)
    if not math.isfinite(result):
        raise QualificationInputError(f"{label} must be a finite number")
    return result


def _quantity(value: Any, factors: dict[str, float], label: str) -> float:
    if not isinstance(value, dict) or not {"value", "unit"} <= set(value):
        raise QualificationInputError(f"{label} must contain value and unit")
    unit = value["unit"]
    if not isinstance(unit, str) or unit not in factors:
        raise QualificationInputError(f"unsupported {label} unit: {unit!r}")
    return _finite_number(value["value"], label) * factors[unit]


def _length_m(value: Any, label: str) -> float:
    return _quantity(value, LENGTH_FACTORS_M, label)


def _volume_m3(value: Any, label: str) -> float:
    return _quantity(value, VOLUME_FACTORS_M3, label)


def _vector_m(value: Any, dimensions: int, label: str) -> tuple[float, ...]:
    if not isinstance(value, dict) or not isinstance(value.get("value"), list):
        raise QualificationInputError(f"{label} must be a unit-tagged coordinate vector")
    if len(value["value"]) != dimensions:
        raise QualificationInputError(f"{label} must have {dimensions} coordinates")
    unit = value.get("unit")
    if unit not in LENGTH_FACTORS_M:
        raise QualificationInputError(f"unsupported {label} unit: {unit!r}")
    factor = LENGTH_FACTORS_M[unit]
    return tuple(_finite_number(x, label) * factor for x in value["value"])


def _source_ref(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise QualificationInputError(f"{label} requires a non-empty source_ref")
    return value


def _vec_sub(a: tuple[float, ...], b: tuple[float, ...]) -> tuple[float, ...]:
    return tuple(x - y for x, y in zip(a, b))


def _dot(a: tuple[float, ...], b: tuple[float, ...]) -> float:
    return sum(x * y for x, y in zip(a, b))


def _norm(a: tuple[float, ...]) -> float:
    return math.sqrt(_dot(a, a))


def _quantile(values: list[float], q: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    index = (len(ordered) - 1) * q
    low = math.floor(index)
    high = math.ceil(index)
    if low == high:
        return ordered[low]
    return ordered[low] + (ordered[high] - ordered[low]) * (index - low)


def _strict_positive(value: float, label: str) -> float:
    if not math.isfinite(value) or value <= 0.0:
        raise QualificationInputError(f"{label} must be positive")
    return value


def _angle_unit_value(record: Any, label: str) -> float:
    if not isinstance(record, dict) or record.get("unit") != "rad":
        raise QualificationInputError(f"{label} must be tagged in rad")
    return _finite_number(record.get("value"), label)


def _parse_geometry(payload: dict[str, Any]) -> dict[str, Any]:
    geometry = payload.get("geometry")
    if not isinstance(geometry, dict) or geometry.get("type") != "hollow_cylinder_patch":
        raise QualificationInputError("geometry.type must be hollow_cylinder_patch")
    if geometry.get("angle_unit") != "rad":
        raise QualificationInputError("geometry.angle_unit must be rad")
    center = _vector_m(geometry.get("axis_center"), 2, "geometry.axis_center")
    rin = _strict_positive(_length_m(geometry.get("inner_radius"), "inner_radius"), "inner_radius")
    rout = _strict_positive(_length_m(geometry.get("outer_radius"), "outer_radius"), "outer_radius")
    if rout <= rin:
        raise QualificationInputError("outer_radius must exceed inner_radius")
    zrange = geometry.get("z_interval")
    if not isinstance(zrange, list) or len(zrange) != 2:
        raise QualificationInputError("geometry.z_interval must contain two tagged lengths")
    zlo = _length_m(zrange[0], "z_interval[0]")
    zhi = _length_m(zrange[1], "z_interval[1]")
    if zhi <= zlo:
        raise QualificationInputError("geometry.z_interval must have positive length")
    patch = geometry.get("theta_patch")
    if not isinstance(patch, list) or len(patch) != 2:
        raise QualificationInputError("geometry.theta_patch must contain two rad values")
    theta_lo = _angle_unit_value(patch[0], "theta_patch[0]")
    theta_hi = _angle_unit_value(patch[1], "theta_patch[1]")
    if theta_hi <= theta_lo or theta_hi - theta_lo >= math.tau:
        raise QualificationInputError("geometry.theta_patch must be a positive patch shorter than 2*pi")
    ref_volume = 0.5 * (rout * rout - rin * rin) * (theta_hi - theta_lo) * (zhi - zlo)
    return {
        "axis_center_m": center,
        "inner_radius_m": rin,
        "outer_radius_m": rout,
        "z_min_m": zlo,
        "z_max_m": zhi,
        "theta_min_rad": theta_lo,
        "theta_max_rad": theta_hi,
        "reference_volume_m3": ref_volume,
    }


def _cell_geometry(cell: dict[str, Any], geometry: dict[str, Any]) -> dict[str, Any]:
    cell_id = cell.get("cell_id")
    if not isinstance(cell_id, str) or not cell_id:
        raise QualificationInputError("each cell requires a non-empty cell_id")
    source_ref = _source_ref(cell.get("source_ref"), f"cell {cell_id}")
    cg = cell.get("geometry")
    if not isinstance(cg, dict) or cg.get("type") != "annular_wedge_hex":
        raise QualificationInputError(f"cell {cell_id} geometry must be annular_wedge_hex")
    if cg.get("angle_unit") != "rad":
        raise QualificationInputError(f"cell {cell_id} angle_unit must be rad")
    ri = _length_m(cg.get("inner_radius"), f"cell {cell_id} inner_radius")
    ro = _length_m(cg.get("outer_radius"), f"cell {cell_id} outer_radius")
    zlo = _length_m(cg.get("z_min"), f"cell {cell_id} z_min")
    zhi = _length_m(cg.get("z_max"), f"cell {cell_id} z_max")
    patch = cg.get("theta_interval")
    if not isinstance(patch, list) or len(patch) != 2:
        raise QualificationInputError(f"cell {cell_id} theta_interval must have two rad values")
    tlo = _angle_unit_value(patch[0], f"cell {cell_id} theta_min")
    thi = _angle_unit_value(patch[1], f"cell {cell_id} theta_max")
    if not (geometry["inner_radius_m"] <= ri < ro <= geometry["outer_radius_m"]):
        raise QualificationInputError(f"cell {cell_id} radial interval is outside the declared wall")
    if not (geometry["z_min_m"] <= zlo < zhi <= geometry["z_max_m"]):
        raise QualificationInputError(f"cell {cell_id} z interval is outside the declared patch")
    if not (geometry["theta_min_rad"] <= tlo < thi <= geometry["theta_max_rad"]):
        raise QualificationInputError(f"cell {cell_id} angular interval is outside the declared patch")
    dtheta = thi - tlo
    analytic_volume = 0.5 * (ro * ro - ri * ri) * dtheta * (zhi - zlo)
    geometric_max_edge = max(ro - ri, ro * dtheta, zhi - zlo)
    declared_max_edge = _strict_positive(
        _length_m(cell.get("max_edge"), f"cell {cell_id} max_edge"),
        f"cell {cell_id} max_edge",
    )
    if abs(declared_max_edge - geometric_max_edge) > max(1e-12, geometric_max_edge * 1e-9):
        raise QualificationInputError(
            f"cell {cell_id} max_edge does not match its annular-wedge boundary edges"
        )
    return {
        "cell_id": cell_id,
        "source_ref": source_ref,
        "r_inner_m": ri,
        "r_outer_m": ro,
        "z_min_m": zlo,
        "z_max_m": zhi,
        "theta_min_rad": tlo,
        "theta_max_rad": thi,
        "volume_m3": analytic_volume,
        "max_edge_m": geometric_max_edge,
    }


def _check_nonoverlap(cells: list[dict[str, Any]], tolerance: float = 1e-14) -> list[dict[str, Any]]:
    overlaps: list[dict[str, Any]] = []
    for index, left in enumerate(cells):
        for right in cells[index + 1 :]:
            radial = min(left["r_outer_m"], right["r_outer_m"]) - max(left["r_inner_m"], right["r_inner_m"])
            angular = min(left["theta_max_rad"], right["theta_max_rad"]) - max(left["theta_min_rad"], right["theta_min_rad"])
            axial = min(left["z_max_m"], right["z_max_m"]) - max(left["z_min_m"], right["z_min_m"])
            if radial > tolerance and angular > tolerance and axial > tolerance:
                overlaps.append({"cell_ids": [left["cell_id"], right["cell_id"]], "status": "OVERLAP"})
    return overlaps


def _project_to_sample(point: tuple[float, float, float], inner: tuple[float, float, float], normal: tuple[float, float, float]) -> tuple[float, float]:
    delta = _vec_sub(point, inner)
    along = _dot(delta, normal)
    perpendicular = _norm(tuple(delta[i] - along * normal[i] for i in range(3)))
    return along, perpendicular


def _cell_intersection_evidence(
    raw_span: dict[str, Any],
    cell_by_id: dict[str, dict[str, Any]],
    inner: tuple[float, float, float],
    outer: tuple[float, float, float],
    normal: tuple[float, float, float],
    ray_z: float,
    geometry: dict[str, Any],
    tolerance: float,
) -> dict[str, Any]:
    cell_id = raw_span.get("cell_id")
    if cell_id not in cell_by_id:
        raise QualificationInputError(f"wall span refers to unknown cell {cell_id!r}")
    source_ref = _source_ref(raw_span.get("source_ref"), f"wall span for {cell_id}")
    start = _vector_m(raw_span.get("start"), 3, f"wall span {cell_id} start")
    end = _vector_m(raw_span.get("end"), 3, f"wall span {cell_id} end")
    s0, e0 = _project_to_sample(start, inner, normal)
    s1, e1 = _project_to_sample(end, inner, normal)
    if e0 > tolerance or e1 > tolerance:
        raise QualificationInputError(f"wall span {cell_id} is not on the sampled normal line")
    if s1 < s0:
        s0, s1 = s1, s0
    if s0 < -tolerance or s1 > _norm(_vec_sub(outer, inner)) + tolerance or s1 - s0 <= tolerance:
        raise QualificationInputError(f"wall span {cell_id} has invalid projected thickness bounds")
    cell = cell_by_id[cell_id]
    center = geometry["axis_center_m"]
    dx = inner[0] - center[0]
    dy = inner[1] - center[1]
    angle = math.atan2(normal[1], normal[0])
    while angle < geometry["theta_min_rad"]:
        angle += math.tau
    while angle >= geometry["theta_min_rad"] + math.tau:
        angle -= math.tau
    if not (
        cell["theta_min_rad"] - tolerance <= angle <= cell["theta_max_rad"] + tolerance
        and cell["z_min_m"] - tolerance <= ray_z <= cell["z_max_m"] + tolerance
    ):
        raise QualificationInputError(f"wall span {cell_id} cell does not contain the sample ray")
    expected_lo = cell["r_inner_m"] - geometry["inner_radius_m"]
    expected_hi = cell["r_outer_m"] - geometry["inner_radius_m"]
    if abs(s0 - expected_lo) > tolerance or abs(s1 - expected_hi) > tolerance:
        raise QualificationInputError(f"wall span {cell_id} disagrees with exact analytic cell geometry")
    return {
        "cell_id": cell_id,
        "source_ref": source_ref,
        "start_depth_m": s0,
        "end_depth_m": s1,
        "length_m": s1 - s0,
        "max_edge_m": cell["max_edge_m"],
    }


def _wall_sample_report(
    sample: dict[str, Any],
    sample_index: int,
    geometry: dict[str, Any],
    cell_by_id: dict[str, dict[str, Any]],
    points_by_cell: dict[str, list[dict[str, Any]]],
    tolerance: float,
) -> dict[str, Any]:
    sample_id = sample.get("sample_id")
    if not isinstance(sample_id, str) or not sample_id:
        raise QualificationInputError(f"wall_samples[{sample_index}] requires sample_id")
    sample_source = _source_ref(sample.get("source_ref"), f"wall sample {sample_id}")
    inner = _vector_m(sample.get("inner_point"), 3, f"wall sample {sample_id} inner_point")
    outer = _vector_m(sample.get("outer_point"), 3, f"wall sample {sample_id} outer_point")
    normal_raw = sample.get("normal")
    if not isinstance(normal_raw, list) or len(normal_raw) != 3:
        raise QualificationInputError(f"wall sample {sample_id} normal must have three components")
    normal = tuple(_finite_number(x, f"wall sample {sample_id} normal") for x in normal_raw)
    nlen = _norm(normal)
    if abs(nlen - 1.0) > 1e-10:
        raise QualificationInputError(f"wall sample {sample_id} normal must be unit length")
    length = _norm(_vec_sub(outer, inner))
    if length <= 0.0:
        raise QualificationInputError(f"wall sample {sample_id} has zero thickness")
    direction = tuple(x / length for x in _vec_sub(outer, inner))
    if _norm(_vec_sub(direction, normal)) > 1e-9:
        raise QualificationInputError(f"wall sample {sample_id} normal does not point inner-to-outer")
    center = geometry["axis_center_m"]
    rin_point = math.hypot(inner[0] - center[0], inner[1] - center[1])
    rout_point = math.hypot(outer[0] - center[0], outer[1] - center[1])
    if abs(rin_point - geometry["inner_radius_m"]) > tolerance or abs(rout_point - geometry["outer_radius_m"]) > tolerance:
        raise QualificationInputError(f"wall sample {sample_id} endpoints disagree with declared cylinder radii")
    if abs(inner[2] - outer[2]) > tolerance or not (geometry["z_min_m"] <= inner[2] <= geometry["z_max_m"]):
        raise QualificationInputError(f"wall sample {sample_id} is outside the declared axial patch")
    raw_intersections = sample.get("cell_intersections")
    if not isinstance(raw_intersections, list):
        raise QualificationInputError(f"wall sample {sample_id} cell_intersections must be a list")
    spans = [
        _cell_intersection_evidence(item, cell_by_id, inner, outer, normal, inner[2], geometry, tolerance)
        for item in raw_intersections
    ]
    if len({row["cell_id"] for row in spans}) != len(spans):
        raise QualificationInputError(f"wall sample {sample_id} repeats a cell intersection")
    spans.sort(key=lambda row: (row["start_depth_m"], row["end_depth_m"], row["cell_id"]))
    cursor = 0.0
    gap_count = 0
    overlap_count = 0
    for span in spans:
        if span["start_depth_m"] > cursor + tolerance:
            gap_count += 1
        if span["start_depth_m"] < cursor - tolerance:
            overlap_count += 1
        cursor = max(cursor, span["end_depth_m"])
    if cursor < length - tolerance:
        gap_count += 1
    if cursor > length + tolerance:
        overlap_count += 1

    traversed = {span["cell_id"] for span in spans}
    centroid_depths: list[float] = []
    centroid_layer_evidence: list[dict[str, Any]] = []
    for cell_id in traversed:
        point_rows = points_by_cell.get(cell_id, [])
        if not point_rows:
            continue
        total_weight = sum(point["weight_m3"] for point in point_rows)
        centroid = tuple(
            sum(point["position_m"][axis] * point["weight_m3"] for point in point_rows) / total_weight
            for axis in range(3)
        )
        along, perpendicular = _project_to_sample(centroid, inner, normal)
        if tolerance < along < length - tolerance:
            centroid_depths.append(along)
            centroid_layer_evidence.append({
                "cell_id": cell_id,
                "projected_depth_m": along,
                "lateral_offset_m": perpendicular,
                "combined_pointweight_m3": total_weight,
                "centroid_count_combined": len(point_rows),
            })
    centroid_depths.sort()
    distinct_centroid_depths: list[float] = []
    for depth in centroid_depths:
        if not distinct_centroid_depths or abs(depth - distinct_centroid_depths[-1]) > max(tolerance, length * 1e-8):
            distinct_centroid_depths.append(depth)

    ratios = [span["max_edge_m"] / length for span in spans]
    ratio_by_cell = [
        {"cell_id": span["cell_id"], "h_over_t": span["max_edge_m"] / length}
        for span in spans
    ]
    contiguous = bool(spans) and gap_count == 0 and overlap_count == 0 and abs(cursor - length) <= tolerance
    return {
        "sample_id": sample_id,
        "source_ref": sample_source,
        "wall_thickness_m": length,
        "wall_thickness_mm": length * 1e3,
        "crossed_cell_intervals": len(spans),
        "crossed_cell_ids": [span["cell_id"] for span in spans],
        "actual_interval_evidence": spans,
        "actual_thickness_centroid_layer_count": len(distinct_centroid_depths),
        "centroid_layer_depths_m": distinct_centroid_depths,
        "centroid_layer_evidence": sorted(centroid_layer_evidence, key=lambda row: (row["projected_depth_m"], row["cell_id"])),
        "actual_cell_h_over_t_values": ratio_by_cell,
        "h_over_t": {"min": min(ratios) if ratios else None, "median": _quantile(ratios, 0.5), "max": max(ratios) if ratios else None},
        "interval_partition_status": "PASS_CONTIGUOUS_LOCAL_WALL_INTERVAL" if contiguous else "FAIL_GAP_OR_OVERLAP",
        "interval_gap_count": gap_count,
        "interval_overlap_count": overlap_count,
    }


def _circle_roots(
    p0: tuple[float, float, float],
    delta: tuple[float, float, float],
    center: tuple[float, float],
    radius: float,
) -> list[float]:
    x = p0[0] - center[0]
    y = p0[1] - center[1]
    dx, dy = delta[0], delta[1]
    a = dx * dx + dy * dy
    b = 2.0 * (x * dx + y * dy)
    c = x * x + y * y - radius * radius
    if a <= 1e-30:
        return []
    discriminant = b * b - 4.0 * a * c
    if discriminant < 0.0:
        return []
    root = math.sqrt(max(0.0, discriminant))
    return [(-b - root) / (2.0 * a), (-b + root) / (2.0 * a)]


def _classify_material_point(
    point: tuple[float, float, float], geometry: dict[str, Any], tolerance: float
) -> bool:
    center = geometry["axis_center_m"]
    radius = math.hypot(point[0] - center[0], point[1] - center[1])
    return (
        geometry["inner_radius_m"] - tolerance <= radius <= geometry["outer_radius_m"] + tolerance
        and geometry["z_min_m"] - tolerance <= point[2] <= geometry["z_max_m"] + tolerance
    )


def _segment_crossing_check(segment: dict[str, Any], geometry: dict[str, Any], tolerance: float) -> dict[str, Any]:
    segment_id = segment.get("segment_id")
    if not isinstance(segment_id, str) or not segment_id:
        raise QualificationInputError("each segment requires segment_id")
    source_ref = _source_ref(segment.get("source_ref"), f"segment {segment_id}")
    start = _vector_m(segment.get("start"), 3, f"segment {segment_id} start")
    end = _vector_m(segment.get("end"), 3, f"segment {segment_id} end")
    delta = _vec_sub(end, start)
    length = _norm(delta)
    if length <= 0.0:
        raise QualificationInputError(f"segment {segment_id} has zero length")
    endpoints_material = _classify_material_point(start, geometry, tolerance) and _classify_material_point(end, geometry, tolerance)

    breaks = [0.0, 1.0]
    for radius in (geometry["inner_radius_m"], geometry["outer_radius_m"]):
        breaks.extend(root for root in _circle_roots(start, delta, geometry["axis_center_m"], radius) if 0.0 < root < 1.0)
    if abs(delta[2]) > 1e-30:
        for zbound in (geometry["z_min_m"], geometry["z_max_m"]):
            value = (zbound - start[2]) / delta[2]
            if 0.0 < value < 1.0:
                breaks.append(value)
    breaks.sort()
    unique_breaks: list[float] = []
    for value in breaks:
        if not unique_breaks or abs(value - unique_breaks[-1]) > 1e-12:
            unique_breaks.append(value)
    bad_intervals: list[dict[str, Any]] = []
    for lo, hi in zip(unique_breaks, unique_breaks[1:]):
        if hi - lo <= 1e-14:
            continue
        mid = 0.5 * (lo + hi)
        point = tuple(start[i] + mid * delta[i] for i in range(3))
        if not _classify_material_point(point, geometry, tolerance):
            radial = math.hypot(point[0] - geometry["axis_center_m"][0], point[1] - geometry["axis_center_m"][1])
            if radial < geometry["inner_radius_m"] - tolerance:
                reason = "CROSSES_HOLLOW_CAVITY"
            elif radial > geometry["outer_radius_m"] + tolerance:
                reason = "LEAVES_OUTER_MATERIAL_BOUNDARY"
            else:
                reason = "LEAVES_AXIAL_MATERIAL_BOUNDS"
            bad_intervals.append({"parameter_interval": [lo, hi], "reason": reason})
    passed = endpoints_material and not bad_intervals
    return {
        "segment_id": segment_id,
        "source_ref": source_ref,
        "length_m": length,
        "endpoint_classification": "BOTH_IN_MATERIAL" if endpoints_material else "ONE_OR_MORE_ENDPOINTS_OUTSIDE_MATERIAL",
        "checked_parameter_intervals": max(0, len(unique_breaks) - 1),
        "violating_intervals": bad_intervals,
        "status": "PASS_ENTIRE_SEGMENT_IN_MATERIAL" if passed else "FAIL_MATERIAL_CROSSING",
    }


def _declared_sample_design_matches(
    samples: list[dict[str, Any]], cells: list[dict[str, Any]], coverage: dict[str, Any], expected_count: int
) -> bool:
    """Cross-check the declared fixture sampling layout against actual ray coordinates."""
    if coverage.get("sample_design_complete") is not True:
        return False
    expected_design = "one exact wall-normal line through every local angular-sector/axial-slab center"
    if coverage.get("sample_design") != expected_design:
        return False
    theta_centers = sorted({round(0.5 * (cell["theta_min_rad"] + cell["theta_max_rad"]), 12) for cell in cells})
    z_centers = sorted({round(0.5 * (cell["z_min_m"] + cell["z_max_m"]), 15) for cell in cells})
    expected = {(theta, z) for theta in theta_centers for z in z_centers}
    if len(expected) != expected_count or len(samples) != expected_count:
        return False
    actual: set[tuple[float, float]] = set()
    # Exact sample endpoints, radii, and directions are checked after this
    # layout check. Here we reject duplicate rays that could masquerade as a
    # complete declared sample design.
    for sample in samples:
        if not isinstance(sample, dict):
            return False
        normal = sample.get("normal")
        inner = sample.get("inner_point")
        if not isinstance(normal, list) or len(normal) != 3:
            return False
        if not isinstance(inner, dict) or not isinstance(inner.get("value"), list) or len(inner["value"]) != 3 or inner.get("unit") not in LENGTH_FACTORS_M:
            return False
        nx = _finite_number(normal[0], "sample normal")
        ny = _finite_number(normal[1], "sample normal")
        theta = math.atan2(ny, nx)
        z = _finite_number(inner["value"][2], "sample z") * LENGTH_FACTORS_M[inner["unit"]]
        actual.add((round(theta, 12), round(z, 15)))
    return len(actual) == expected_count and actual == expected


def check_segments(
    segments: Iterable[dict[str, Any]], geometry_payload: dict[str, Any], *, tolerance_m: float = 1e-12
) -> dict[str, Any]:
    """Adapter entry point for bond/pair segments from an actual mesh producer.

    Each record must provide unit-tagged start/end 3D coordinates and a source_ref.
    The checker intersects the entire segment with exact inner/outer cylinder
    quadratics and z planes; material endpoints alone never imply a pass.
    """
    geometry = _parse_geometry({"geometry": geometry_payload})
    rows = list(segments)
    if len(rows) > MAX_SEGMENTS:
        raise QualificationInputError(f"segment adapter exceeds {MAX_SEGMENTS} records")
    reports = [_segment_crossing_check(row, geometry, tolerance_m) for row in rows]
    failures = [row["segment_id"] for row in reports if row["status"] != "PASS_ENTIRE_SEGMENT_IN_MATERIAL"]
    return {
        "adapter": "check_segments(unit_tagged_endpoints, exact_hollow_cylinder_geometry)",
        "segment_count": len(reports),
        "status": "UNKNOWN_NO_SEGMENTS" if not reports else ("PASS_NO_CROSSINGS_IN_CHECKED_SEGMENTS" if not failures else "FAIL_CROSSINGS_DETECTED"),
        "failure_segment_ids": failures,
        "segments": reports,
        "scope": "CHECKED_SEGMENTS_ONLY; DOES_NOT_ESTABLISH_UNCHECKED_GRAPH_COVERAGE",
    }


def qualify_wall_input(payload: dict[str, Any]) -> dict[str, Any]:
    """Validate units, geometry-derived point weights, local layers, and segments."""
    if not isinstance(payload, dict) or payload.get("schema") != SCHEMA:
        raise QualificationInputError(f"schema must be {SCHEMA}")
    declared_units = payload.get("units")
    if not isinstance(declared_units, dict):
        raise QualificationInputError("units must declare length, volume, and mass")
    if declared_units.get("length") not in LENGTH_FACTORS_M:
        raise QualificationInputError(f"unsupported declared length unit: {declared_units.get('length')!r}")
    if declared_units.get("volume") not in VOLUME_FACTORS_M3:
        raise QualificationInputError(f"unsupported declared volume unit: {declared_units.get('volume')!r}")
    if declared_units.get("mass") != "kg":
        raise QualificationInputError("only kg mass inputs are supported")
    case_id = payload.get("case_id")
    if not isinstance(case_id, str) or not case_id:
        raise QualificationInputError("case_id is required")
    density = _finite_number(payload.get("density_kg_m3"), "density_kg_m3")
    _strict_positive(density, "density_kg_m3")
    geometry = _parse_geometry(payload)
    source_ref = _source_ref(payload.get("source_ref"), "input")
    cells_raw = payload.get("cells")
    if not isinstance(cells_raw, list) or not cells_raw:
        raise QualificationInputError("cells must be a non-empty list")
    if len(cells_raw) > MAX_CELLS:
        raise QualificationInputError(f"cell input exceeds {MAX_CELLS} cells")
    cells = [_cell_geometry(cell, geometry) for cell in cells_raw]
    cell_ids = [cell["cell_id"] for cell in cells]
    if len(set(cell_ids)) != len(cell_ids):
        raise QualificationInputError("cell_id values must be unique")
    cell_by_id = {cell["cell_id"]: cell for cell in cells}
    overlaps = _check_nonoverlap(cells)

    point_rows = payload.get("centroid_pointweights")
    if not isinstance(point_rows, list) or not point_rows:
        raise QualificationInputError("centroid_pointweights must be a non-empty list")
    points_by_cell: dict[str, list[dict[str, Any]]] = defaultdict(list)
    point_reports: list[dict[str, Any]] = []
    point_volume_sums: dict[str, float] = defaultdict(float)
    point_mass_sum = 0.0
    point_volume_sum = 0.0
    tolerance = max(1e-12, min(1e-9, (geometry["outer_radius_m"] - geometry["inner_radius_m"]) * 1e-8))
    for index, point in enumerate(point_rows):
        point_id = point.get("point_id")
        if not isinstance(point_id, str) or not point_id:
            raise QualificationInputError(f"centroid_pointweights[{index}] requires point_id")
        psource = _source_ref(point.get("source_ref"), f"point {point_id}")
        cell_id = point.get("cell_id")
        if cell_id not in cell_by_id:
            raise QualificationInputError(f"point {point_id} refers to unknown cell {cell_id!r}")
        position = _vector_m(point.get("position"), 3, f"point {point_id} position")
        weight = _strict_positive(_volume_m3(point.get("point_weight"), f"point {point_id} point_weight"), f"point {point_id} point_weight")
        if not _classify_material_point(position, geometry, tolerance):
            raise QualificationInputError(f"point {point_id} centroid is outside the hollow-cylinder material")
        cg = cell_by_id[cell_id]
        angle = math.atan2(position[1] - geometry["axis_center_m"][1], position[0] - geometry["axis_center_m"][0])
        while angle < geometry["theta_min_rad"]:
            angle += math.tau
        while angle >= geometry["theta_min_rad"] + math.tau:
            angle -= math.tau
        if not (
            cg["r_inner_m"] - tolerance <= math.hypot(position[0] - geometry["axis_center_m"][0], position[1] - geometry["axis_center_m"][1]) <= cg["r_outer_m"] + tolerance
            and cg["theta_min_rad"] - tolerance <= angle <= cg["theta_max_rad"] + tolerance
            and cg["z_min_m"] - tolerance <= position[2] <= cg["z_max_m"] + tolerance
        ):
            raise QualificationInputError(f"point {point_id} centroid is outside its source cell")
        expected_mass = density * weight
        mass_record = point.get("mass")
        if mass_record is None:
            raise QualificationInputError(f"point {point_id} requires a unit-tagged mass")
        mass_kg = _quantity(mass_record, {"kg": 1.0}, f"point {point_id} mass")
        points_by_cell[cell_id].append({"point_id": point_id, "source_ref": psource, "position_m": position, "weight_m3": weight})
        point_volume_sums[cell_id] += weight
        point_volume_sum += weight
        point_mass_sum += mass_kg
        point_reports.append({
            "point_id": point_id,
            "cell_id": cell_id,
            "source_ref": psource,
            "point_weight_m3": weight,
            "mass_kg": mass_kg,
            "rho_times_weight_kg": expected_mass,
            "mass_delta_kg": mass_kg - expected_mass,
        })

    cell_volume_sum = sum(cell["volume_m3"] for cell in cells)
    reference_volume = geometry["reference_volume_m3"]
    point_volume_tolerance = max(1e-16, reference_volume * 1e-10)
    per_cell_errors = []
    for cell in cells:
        actual = point_volume_sums.get(cell["cell_id"], 0.0)
        delta = actual - cell["volume_m3"]
        if abs(delta) > max(1e-16, cell["volume_m3"] * 1e-10):
            per_cell_errors.append({"cell_id": cell["cell_id"], "geometry_volume_m3": cell["volume_m3"], "point_weight_sum_m3": actual, "delta_m3": delta})
    mass_expected = density * reference_volume
    mass_delta = point_mass_sum - mass_expected
    volume_consistency = {
        "analytic_patch_reference_volume_m3": reference_volume,
        "sum_exact_cell_geometry_volumes_m3": cell_volume_sum,
        "sum_centroid_pointweights_m3": point_volume_sum,
        "volume_tolerance_m3": point_volume_tolerance,
        "cell_geometry_partition_volume_delta_m3": cell_volume_sum - reference_volume,
        "pointweight_volume_delta_m3": point_volume_sum - reference_volume,
        "point_mass_sum_kg": point_mass_sum,
        "rho_times_reference_volume_kg": mass_expected,
        "mass_tolerance_kg": max(1e-15, abs(mass_expected) * 1e-10),
        "mass_delta_kg": mass_delta,
        "per_cell_pointweight_mismatches": per_cell_errors,
        "cell_geometry_nonoverlap_count": len(overlaps),
        "volume_status": "PASS_LOCAL_GEOMETRY_POINTWEIGHT_RECONCILIATION"
        if not per_cell_errors and not overlaps and abs(cell_volume_sum - reference_volume) <= point_volume_tolerance and abs(point_volume_sum - reference_volume) <= point_volume_tolerance
        else "FAIL_VOLUME_OR_WEIGHT_RECONCILIATION",
        "mass_status": "PASS_DENSITY_TIMES_GEOMETRIC_VOLUME"
        if not per_cell_errors and abs(mass_delta) <= max(1e-15, abs(mass_expected) * 1e-10) and all(abs(row["mass_delta_kg"]) <= max(1e-15, abs(row["rho_times_weight_kg"]) * 1e-10) for row in point_reports)
        else "FAIL_MASS_RECONCILIATION",
        "overlapping_cell_pairs": overlaps,
    }

    sampling_policy = payload.get("sampling_policy")
    if not isinstance(sampling_policy, dict):
        raise QualificationInputError("sampling_policy is required")
    minimum_intervals = sampling_policy.get("min_thickness_sampling_intervals")
    minimum_centroids = sampling_policy.get("min_distinct_centroid_layers")
    for name, value in (("min_thickness_sampling_intervals", minimum_intervals), ("min_distinct_centroid_layers", minimum_centroids)):
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise QualificationInputError(f"sampling_policy.{name} must be a positive integer")
    process_zone_status = sampling_policy.get("process_zone_status")
    if process_zone_status != "REQUIRES_SEPARATE_VALIDATION":
        raise QualificationInputError("process-zone qualification must remain explicitly separate")

    coverage = payload.get("coverage")
    coverage_present = isinstance(coverage, dict)
    wall_reports: list[dict[str, Any]] = []
    if coverage_present:
        if coverage.get("scope") != "LOCAL_WALL_PATCH_ONLY":
            raise QualificationInputError("coverage.scope must be LOCAL_WALL_PATCH_ONLY")
        expected_sample_count = coverage.get("expected_wall_sample_count")
        if isinstance(expected_sample_count, bool) or not isinstance(expected_sample_count, int) or expected_sample_count < 1:
            raise QualificationInputError("coverage.expected_wall_sample_count must be positive")
        samples = payload.get("wall_samples")
        if not isinstance(samples, list):
            raise QualificationInputError("wall_samples must be a list when coverage is declared")
        ids = [sample.get("sample_id") for sample in samples if isinstance(sample, dict)]
        if len(ids) != len(samples) or len(set(ids)) != len(ids):
            raise QualificationInputError("wall sample ids must be present and unique")
        sample_design_complete = _declared_sample_design_matches(
            samples, cells, coverage, expected_sample_count
        )
        if sample_design_complete:
            wall_reports = [
                _wall_sample_report(sample, i, geometry, cell_by_id, points_by_cell, tolerance)
                for i, sample in enumerate(samples)
            ]
        else:
            wall_reports = []
    else:
        expected_sample_count = None
        samples = payload.get("wall_samples") if isinstance(payload.get("wall_samples"), list) else []
        sample_design_complete = False

    raw_segments = payload.get("segments", [])
    if not isinstance(raw_segments, list):
        raise QualificationInputError("segments must be a list")
    if len(raw_segments) > MAX_SEGMENTS:
        raise QualificationInputError(f"segment input exceeds {MAX_SEGMENTS} segments")
    segment_result = check_segments(raw_segments, {
        "type": "hollow_cylinder_patch",
        "angle_unit": "rad",
        "axis_center": {"value": [geometry["axis_center_m"][0] * 1e3, geometry["axis_center_m"][1] * 1e3], "unit": "mm"},
        "inner_radius": {"value": geometry["inner_radius_m"] * 1e3, "unit": "mm"},
        "outer_radius": {"value": geometry["outer_radius_m"] * 1e3, "unit": "mm"},
        "z_interval": [
            {"value": geometry["z_min_m"] * 1e3, "unit": "mm"},
            {"value": geometry["z_max_m"] * 1e3, "unit": "mm"},
        ],
        "theta_patch": [
            {"value": geometry["theta_min_rad"], "unit": "rad"},
            {"value": geometry["theta_max_rad"], "unit": "rad"},
        ],
    }, tolerance_m=tolerance)

    local_domain_complete = (
        coverage_present
        and coverage.get("local_cell_geometry_complete") is True
        and sample_design_complete
        and not overlaps
        and abs(cell_volume_sum - reference_volume) <= point_volume_tolerance
    )
    local_samples_complete = bool(wall_reports) and all(
        row["interval_partition_status"] == "PASS_CONTIGUOUS_LOCAL_WALL_INTERVAL" for row in wall_reports
    )
    min_observed_intervals = min((row["crossed_cell_intervals"] for row in wall_reports), default=None)
    min_observed_centroids = min((row["actual_thickness_centroid_layer_count"] for row in wall_reports), default=None)
    max_h_over_t = max(
        (row["h_over_t"]["max"] for row in wall_reports if row["h_over_t"]["max"] is not None),
        default=None,
    )
    if volume_consistency["volume_status"].startswith("FAIL") or volume_consistency["mass_status"].startswith("FAIL"):
        local_status = "BLOCKED_INCONSISTENT_VOLUME_OR_MASS"
    elif segment_result["status"] == "FAIL_CROSSINGS_DETECTED":
        local_status = "BLOCKED_CAVITY_OR_BOUNDARY_CROSSING_SEGMENT"
    elif segment_result["status"] == "UNKNOWN_NO_SEGMENTS":
        local_status = "QUALIFICATION_UNKNOWN_NO_SEGMENT_EVIDENCE"
    elif not coverage_present or not local_domain_complete or not local_samples_complete:
        local_status = "QUALIFICATION_UNKNOWN_MISSING_OR_INCOMPLETE_LOCAL_COVERAGE"
    elif (
        min_observed_intervals is not None
        and min_observed_intervals < minimum_intervals
    ) or (
        min_observed_centroids is not None
        and min_observed_centroids < minimum_centroids
    ) or (max_h_over_t is not None and max_h_over_t >= 1.0):
        local_status = "UNDERRESOLVED"
    else:
        local_status = "LOCAL_RESOLUTION_PRECONDITION_MET"

    flattened_ratios = [
        entry["h_over_t"]
        for row in wall_reports
        for entry in row["actual_cell_h_over_t_values"]
    ]
    return {
        "schema": "pgfrag.thin_wall_resolution_report.v1",
        "case_id": case_id,
        "source_ref": source_ref,
        "inputs": {
            "declared_units": declared_units,
            "cell_count": len(cells),
            "centroid_pointweight_count": len(point_rows),
            "wall_sample_count": len(samples),
            "segment_count": len(raw_segments),
            "fixture_cell_cap": MAX_CELLS,
            "fixture_segment_cap": MAX_SEGMENTS,
        },
        "local_scope": {
            "coverage_scope": coverage.get("scope") if coverage_present else "MISSING",
            "local_cell_geometry_complete": local_domain_complete,
            "wall_sample_design_complete": sample_design_complete,
            "full_domain_coverage_declared_by_input": bool(coverage_present and coverage.get("full_domain_coverage") is True),
            "full_domain_coverage_verified": False,
            "global_status": "UNKNOWN_NO_FULL_DOMAIN_COVERAGE",
            "global_pass": False,
            "scope_note": "Only the supplied local patch and listed sample rays were examined; no local result is promoted to a whole-body result.",
        },
        "resolution": {
            "local_status": local_status,
            "engineering_precondition": {
                "min_actual_thickness_sampling_intervals_per_ray": minimum_intervals,
                "min_distinct_centroid_layers_per_ray": minimum_centroids,
                "rationale": sampling_policy.get("rationale"),
                "universality": "DECLARED_ENGINEERING_PRECONDITION_ONLY; NOT A UNIVERSAL CRACK-RESOLUTION RULE",
                "process_zone_status": process_zone_status,
            },
            "min_observed_actual_thickness_intervals": min_observed_intervals,
            "min_observed_distinct_centroid_layers": min_observed_centroids,
            "actual_h_over_t_distribution": {
                "sampled_interval_count": len(flattened_ratios),
                "min": min(flattened_ratios) if flattened_ratios else None,
                "p50": _quantile(flattened_ratios, 0.5),
                "max": max(flattened_ratios) if flattened_ratios else None,
                "interpretation": "Descriptive max-edge / sampled-wall-thickness ratio; it does not replace measured interval or centroid-layer evidence.",
            },
            "wall_samples": wall_reports,
            "physical_status": "NOT_QUALIFIED",
            "physical_pass": False,
        },
        "volume_mass_consistency": volume_consistency,
        "cavity_crossing_segments": segment_result,
        "decision_limits": [
            "A local precondition result is not a physicalPASS.",
            "Three actual intervals and three centroid layers are only the declared fixture precondition; process-zone adequacy requires separate evidence.",
            "An attractive or smooth section image is not used by any gate.",
        ],
    }


def make_analytic_hollow_cylinder_fixture(
    h_mm: float,
    *,
    coverage_present: bool = True,
    coverage_complete: bool = True,
    source_ref: str | None = None,
) -> dict[str, Any]:
    """Build a local exact annular-wedge fixture with positive geometric weights."""
    h_mm = _finite_number(h_mm, "fixture h_mm")
    if h_mm not in (0.5, 2.0):
        raise QualificationInputError("fixtures are restricted to h=0.5 mm or h=2.0 mm")
    rin_mm, rout_mm = 8.0, 10.0
    zlo_mm, zhi_mm = 4.0, 6.0
    rho = 2100.0
    nradial = int(round((rout_mm - rin_mm) / h_mm))
    ntheta = 4
    nz = 4
    theta_step = h_mm / rout_mm if h_mm == 0.5 else 0.05
    theta_lo = -0.5 * ntheta * theta_step
    theta_hi = theta_lo + ntheta * theta_step
    z_edges = [zlo_mm + (zhi_mm - zlo_mm) * i / nz for i in range(nz + 1)]
    radial_edges = [rin_mm + (rout_mm - rin_mm) * i / nradial for i in range(nradial + 1)]
    theta_edges = [theta_lo + (theta_hi - theta_lo) * i / ntheta for i in range(ntheta + 1)]
    source_ref = source_ref or f"fixture:hollow-cylinder-R8-R10mm:h{h_mm:g}mm"

    def lq(value: float) -> dict[str, Any]:
        return {"value": value, "unit": "mm"}

    def vq(value: float) -> dict[str, Any]:
        return {"value": value, "unit": "mm^3"}

    def xyzq(value: tuple[float, float, float]) -> dict[str, Any]:
        return {"value": list(value), "unit": "mm"}

    cells: list[dict[str, Any]] = []
    pointweights: list[dict[str, Any]] = []
    cell_id_at: dict[tuple[int, int, int], str] = {}
    for ir in range(nradial):
        ri, ro = radial_edges[ir], radial_edges[ir + 1]
        for ia in range(ntheta):
            tlo, thi = theta_edges[ia], theta_edges[ia + 1]
            dtheta = thi - tlo
            radial_bar = (2.0 / 3.0) * (ro**3 - ri**3) / (ro**2 - ri**2)
            angular_factor = math.sin(dtheta / 2.0) / (dtheta / 2.0)
            rcentroid = radial_bar * angular_factor
            tcentroid = 0.5 * (tlo + thi)
            for iz in range(nz):
                z0, z1 = z_edges[iz], z_edges[iz + 1]
                cell_id = f"c-r{ir:02d}-a{ia:02d}-z{iz:02d}"
                cell_id_at[(ir, ia, iz)] = cell_id
                cell_volume_mm3 = 0.5 * (ro**2 - ri**2) * dtheta * (z1 - z0)
                edge_mm = max(ro - ri, ro * dtheta, z1 - z0)
                centroid_mm = (
                    0.0 + rcentroid * math.cos(tcentroid),
                    0.0 + rcentroid * math.sin(tcentroid),
                    0.5 * (z0 + z1),
                )
                cells.append({
                    "cell_id": cell_id,
                    "source_ref": f"{source_ref}:cell:{cell_id}:analytic-boundary-edges",
                    "geometry": {
                        "type": "annular_wedge_hex",
                        "inner_radius": lq(ri),
                        "outer_radius": lq(ro),
                        "theta_interval": [
                            {"value": tlo, "unit": "rad"},
                            {"value": thi, "unit": "rad"},
                        ],
                        "angle_unit": "rad",
                        "z_min": lq(z0),
                        "z_max": lq(z1),
                    },
                    "max_edge": lq(edge_mm),
                })
                pointweights.append({
                    "point_id": f"p-{cell_id}",
                    "cell_id": cell_id,
                    "source_ref": f"{source_ref}:centroid:{cell_id}:exact-cell-COM",
                    "position": xyzq(centroid_mm),
                    "point_weight": vq(cell_volume_mm3),
                    "mass": {"value": rho * cell_volume_mm3 * 1e-9, "unit": "kg"},
                })

    wall_samples: list[dict[str, Any]] = []
    for ia in range(ntheta):
        tmid = 0.5 * (theta_edges[ia] + theta_edges[ia + 1])
        normal = [math.cos(tmid), math.sin(tmid), 0.0]
        for iz in range(nz):
            zmid = 0.5 * (z_edges[iz] + z_edges[iz + 1])
            inner = (rin_mm * normal[0], rin_mm * normal[1], zmid)
            outer = (rout_mm * normal[0], rout_mm * normal[1], zmid)
            intersections = []
            for ir in range(nradial):
                r0, r1 = radial_edges[ir], radial_edges[ir + 1]
                cell_id = cell_id_at[(ir, ia, iz)]
                intersections.append({
                    "cell_id": cell_id,
                    "source_ref": f"{source_ref}:ray-hit:{cell_id}:exact-radial-cell-span",
                    "start": xyzq((r0 * normal[0], r0 * normal[1], zmid)),
                    "end": xyzq((r1 * normal[0], r1 * normal[1], zmid)),
                })
            wall_samples.append({
                "sample_id": f"ws-a{ia:02d}-z{iz:02d}",
                "source_ref": f"{source_ref}:wall-normal-sample:theta={tmid:.12g}:z={zmid:.12g}mm",
                "normal": normal,
                "inner_point": xyzq(inner),
                "outer_point": xyzq(outer),
                "cell_intersections": intersections,
            })

    segments: list[dict[str, Any]] = []
    for ir in range(nradial):
        for ia in range(ntheta - 1):
            for iz in range(nz):
                left_id = cell_id_at[(ir, ia, iz)]
                right_id = cell_id_at[(ir, ia + 1, iz)]
                left_point = next(row for row in pointweights if row["cell_id"] == left_id)
                right_point = next(row for row in pointweights if row["cell_id"] == right_id)
                segments.append({
                    "segment_id": f"seg:{left_id}>{right_id}",
                    "source_ref": f"{source_ref}:fixture-neighbor-segment:{left_id}>{right_id}",
                    "start": left_point["position"],
                    "end": right_point["position"],
                })

    result: dict[str, Any] = {
        "schema": SCHEMA,
        "case_id": f"hollow-cylinder-t2mm-h{h_mm:g}mm",
        "source_ref": source_ref,
        "units": {"length": "mm", "volume": "mm^3", "mass": "kg"},
        "density_kg_m3": rho,
        "geometry": {
            "type": "hollow_cylinder_patch",
            "axis_center": {"value": [0.0, 0.0], "unit": "mm"},
            "inner_radius": lq(rin_mm),
            "outer_radius": lq(rout_mm),
            "theta_patch": [
                {"value": theta_lo, "unit": "rad"},
                {"value": theta_hi, "unit": "rad"},
            ],
            "angle_unit": "rad",
            "z_interval": [lq(zlo_mm), lq(zhi_mm)],
        },
        "sampling_policy": {
            "min_thickness_sampling_intervals": 3,
            "min_distinct_centroid_layers": 3,
            "process_zone_status": "REQUIRES_SEPARATE_VALIDATION",
            "rationale": "This review declared a three-interval local mesh-precondition for the fixture; crack-process-zone adequacy is material- and model-dependent and is not inferred.",
        },
        "coverage": {
            "scope": "LOCAL_WALL_PATCH_ONLY",
            "local_cell_geometry_complete": coverage_complete,
            "sample_design_complete": coverage_complete,
            "expected_wall_sample_count": len(wall_samples),
            "full_domain_coverage": False,
            "sample_design": "one exact wall-normal line through every local angular-sector/axial-slab center",
        } if coverage_present else None,
        "cells": cells,
        "centroid_pointweights": pointweights,
        "wall_samples": wall_samples,
        "segments": segments,
        "source_evidence": [
            "pgfrag/experiments/forge_volume_ingress_20261005/CONTRACT.md",
            "pgfrag/docs/evidence/20261005/forge_volume_ingress_independent_audit_20261005/audit/PILOT_ACTUAL_AUDIT.json",
            "pgfrag/docs/evidence/20261005/forge_volume_ingress_roi_20261005/failed_run/pilot_results_summary.md",
        ],
        "fixture_disclaimer": "This exact annular-wedge fixture exercises local unit/weight/layer/cavity gates; it is not the canonical CAD body and proves no mother-wide coverage.",
    }
    if not coverage_complete and result["coverage"] is not None:
        result["coverage"]["sample_design_complete"] = False
    return result


def historical_h8_assessment(*, nominal_wall_mm: float = 2.0) -> dict[str, Any]:
    """Explicit disposition for the old coarse h=8 mm proposal."""
    h_m = _length_m({"value": 8.0, "unit": "mm"}, "historical h")
    t_m = _length_m({"value": nominal_wall_mm, "unit": "mm"}, "historical nominal wall")
    return {
        "case_id": "historical-h8mm-nominal-wall-comparison",
        "h_m": h_m,
        "nominal_wall_m": t_m,
        "h_over_t": h_m / t_m,
        "status": "UNDERRESOLVED",
        "basis": "Known coarse edge is at least the nominal wall thickness; no sampled wall-ray or actual through-thickness layer evidence upgrades it.",
        "scope": "Historical nominal-parameter comparison only; it is not a recovered mesh measurement or global wall certificate.",
    }
