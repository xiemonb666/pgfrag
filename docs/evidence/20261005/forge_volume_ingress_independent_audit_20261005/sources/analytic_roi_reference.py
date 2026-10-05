"""Pure-stdlib analytic oracle for the frozen annular-cylinder ROI.

This module deliberately does not import CadQuery/OpenCascade and does not execute
any producer code. It only computes disk-minus-disk rectangle moments and the
root-selected Cartesian cell reference values.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Iterable


OUTER_RADIUS_MM = 10.0
INNER_RADIUS_MM = 8.0
ROI_MM = (7.5, 10.5, -1.0, 1.0, 4.0, 6.0)
CELL_MM = 0.5


@dataclass(frozen=True)
class PlanarMoments:
    area_mm2: float
    first_x_mm3: float
    first_y_mm3: float


@dataclass(frozen=True)
class CellReference:
    grid_cell_index: int
    grid_ijk: tuple[int, int, int]
    cell_aabb_mm: tuple[float, float, float, float, float, float]
    expected_area_mm2: float
    expected_volume_mm3: float
    expected_centroid_mm: tuple[float, float, float] | None


def _circle_s(radius: float, x: float) -> float:
    return math.sqrt(max(0.0, radius * radius - x * x))


def _primitive_s(radius: float, x: float) -> float:
    s = _circle_s(radius, x)
    ratio = max(-1.0, min(1.0, x / radius))
    return 0.5 * (x * s + radius * radius * math.asin(ratio))


def _primitive_xs(radius: float, x: float) -> float:
    return -_circle_s(radius, x) ** 3 / 3.0


def _primitive_half_s2(radius: float, x: float) -> float:
    return 0.5 * (radius * radius * x - x ** 3 / 3.0)


def _boundary_integrals(
    kind: str, value: float, radius: float, x0: float, x1: float
) -> tuple[float, float, float]:
    """Return integral(f), integral(x*f), integral(f^2/2) for one y boundary."""
    if kind == "constant":
        return (
            value * (x1 - x0),
            value * (x1 * x1 - x0 * x0) / 2.0,
            value * value * (x1 - x0) / 2.0,
        )
    if kind == "circle_plus":
        return (
            _primitive_s(radius, x1) - _primitive_s(radius, x0),
            _primitive_xs(radius, x1) - _primitive_xs(radius, x0),
            _primitive_half_s2(radius, x1)
            - _primitive_half_s2(radius, x0),
        )
    if kind == "circle_minus":
        return (
            -(_primitive_s(radius, x1) - _primitive_s(radius, x0)),
            -(_primitive_xs(radius, x1) - _primitive_xs(radius, x0)),
            _primitive_half_s2(radius, x1)
            - _primitive_half_s2(radius, x0),
        )
    raise ValueError(f"unsupported boundary kind: {kind}")


def disk_rectangle_moments(
    radius_mm: float,
    x_bounds_mm: tuple[float, float],
    y_bounds_mm: tuple[float, float],
) -> PlanarMoments:
    """Integrate a radius-r disk clipped to one axis-aligned rectangle.

    The y-slice is clipped analytically. Segment endpoints include every x where
    the circle arc crosses a rectangle y edge, so each segment has fixed active
    top/bottom boundaries and uses closed-form circle-arc primitives.
    """
    x0, x1 = map(float, x_bounds_mm)
    y0, y1 = map(float, y_bounds_mm)
    radius = float(radius_mm)
    if not math.isfinite(radius) or radius <= 0:
        raise ValueError("radius must be positive and finite")
    if not all(math.isfinite(value) for value in (x0, x1, y0, y1)):
        raise ValueError("rectangle bounds must be finite")
    if not (x0 < x1 and y0 < y1):
        raise ValueError("rectangle must have positive area")

    lo = max(x0, -radius)
    hi = min(x1, radius)
    if lo >= hi:
        return PlanarMoments(0.0, 0.0, 0.0)

    breaks = {lo, hi}
    for y_edge in (y0, y1):
        if abs(y_edge) < radius:
            crossing = math.sqrt(max(0.0, radius * radius - y_edge * y_edge))
            for x_cross in (-crossing, crossing):
                if lo < x_cross < hi:
                    breaks.add(x_cross)
    ordered = sorted(breaks)

    area = first_x = first_y = 0.0
    for left, right in zip(ordered, ordered[1:]):
        if not right > left:
            continue
        mid = (left + right) / 2.0
        arc = _circle_s(radius, mid)
        top_mid = min(y1, arc)
        bottom_mid = max(y0, -arc)
        if not top_mid > bottom_mid:
            continue

        top_kind = "circle_plus" if arc < y1 else "constant"
        top_value = 0.0 if top_kind == "circle_plus" else y1
        bottom_kind = "circle_minus" if -arc > y0 else "constant"
        bottom_value = 0.0 if bottom_kind == "circle_minus" else y0

        top = _boundary_integrals(top_kind, top_value, radius, left, right)
        bottom = _boundary_integrals(
            bottom_kind, bottom_value, radius, left, right
        )
        area += top[0] - bottom[0]
        first_x += top[1] - bottom[1]
        first_y += top[2] - bottom[2]

    if area < 0.0 and abs(area) < 1e-12:
        return PlanarMoments(0.0, 0.0, 0.0)
    if area < 0.0:
        raise ArithmeticError(f"negative disk/rectangle area: {area}")
    return PlanarMoments(area, first_x, first_y)


def annulus_rectangle_moments(
    x_bounds_mm: tuple[float, float],
    y_bounds_mm: tuple[float, float],
    *,
    outer_radius_mm: float = OUTER_RADIUS_MM,
    inner_radius_mm: float = INNER_RADIUS_MM,
) -> PlanarMoments:
    """Compute annulus/rectangle moments as outer-disk minus inner-disk."""
    outer = disk_rectangle_moments(
        outer_radius_mm, x_bounds_mm, y_bounds_mm
    )
    inner = disk_rectangle_moments(
        inner_radius_mm, x_bounds_mm, y_bounds_mm
    )
    result = PlanarMoments(
        outer.area_mm2 - inner.area_mm2,
        outer.first_x_mm3 - inner.first_x_mm3,
        outer.first_y_mm3 - inner.first_y_mm3,
    )
    if result.area_mm2 < -1e-10:
        raise ArithmeticError(f"negative annulus/rectangle area: {result.area_mm2}")
    return result


def _axis_cells(lo: float, hi: float, step: float) -> list[tuple[float, float]]:
    count = round((hi - lo) / step)
    if count < 1 or abs((hi - lo) / step - count) > 1e-12:
        raise ValueError("ROI axis is not exactly divisible by cell size")
    return [(lo + i * step, lo + (i + 1) * step) for i in range(count)]


def expected_cells() -> list[CellReference]:
    """Return 96 analytic cells in worker ordering: k, then j, then i."""
    xmin, xmax, ymin, ymax, zmin, zmax = ROI_MM
    xs = _axis_cells(xmin, xmax, CELL_MM)
    ys = _axis_cells(ymin, ymax, CELL_MM)
    zs = _axis_cells(zmin, zmax, CELL_MM)
    rows: list[CellReference] = []
    for k, (z0, z1) in enumerate(zs):
        for j, (y0, y1) in enumerate(ys):
            for i, (x0, x1) in enumerate(xs):
                moments = annulus_rectangle_moments((x0, x1), (y0, y1))
                dz = z1 - z0
                volume = moments.area_mm2 * dz
                centroid = (
                    (
                        moments.first_x_mm3 / moments.area_mm2,
                        moments.first_y_mm3 / moments.area_mm2,
                        (z0 + z1) / 2.0,
                    )
                    if moments.area_mm2 > 0.0
                    else None
                )
                rows.append(
                    CellReference(
                        grid_cell_index=len(rows),
                        grid_ijk=(i, j, k),
                        cell_aabb_mm=(x0, x1, y0, y1, z0, z1),
                        expected_area_mm2=moments.area_mm2,
                        expected_volume_mm3=volume,
                        expected_centroid_mm=centroid,
                    )
                )
    if len(rows) != 96:
        raise ArithmeticError(f"expected exactly 96 cells, got {len(rows)}")
    return rows


def expected_roi() -> dict[str, object]:
    """Return independent direct-ROI and summed-grid analytic checks."""
    xmin, xmax, ymin, ymax, zmin, zmax = ROI_MM
    direct = annulus_rectangle_moments(
        (xmin, xmax), (ymin, ymax)
    )
    depth = zmax - zmin
    direct_volume = direct.area_mm2 * depth
    direct_centroid = (
        (
            direct.first_x_mm3 / direct.area_mm2,
            direct.first_y_mm3 / direct.area_mm2,
            (zmin + zmax) / 2.0,
        )
        if direct.area_mm2 > 0.0
        else None
    )
    cells = expected_cells()
    summed_volume = math.fsum(row.expected_volume_mm3 for row in cells)
    summed_first_x = math.fsum(
        (row.expected_centroid_mm[0] * row.expected_volume_mm3)
        for row in cells
        if row.expected_centroid_mm is not None
    )
    summed_first_y = math.fsum(
        (row.expected_centroid_mm[1] * row.expected_volume_mm3)
        for row in cells
        if row.expected_centroid_mm is not None
    )
    summed_first_z = math.fsum(
        (row.expected_centroid_mm[2] * row.expected_volume_mm3)
        for row in cells
        if row.expected_centroid_mm is not None
    )
    summed_centroid = (
        (summed_first_x / summed_volume,
         summed_first_y / summed_volume,
         summed_first_z / summed_volume)
        if summed_volume > 0.0
        else None
    )
    return {
        "roi_aabb_mm": ROI_MM,
        "outer_radius_mm": OUTER_RADIUS_MM,
        "inner_radius_mm": INNER_RADIUS_MM,
        "cell_size_mm": CELL_MM,
        "cell_count": len(cells),
        "direct_planar_area_mm2": direct.area_mm2,
        "direct_volume_mm3": direct_volume,
        "direct_centroid_mm": direct_centroid,
        "sum_cell_volume_mm3": summed_volume,
        "sum_cell_centroid_mm": summed_centroid,
        "direct_minus_grid_volume_mm3": direct_volume - summed_volume,
    }


def cell_rows_as_dicts(rows: Iterable[CellReference] | None = None) -> list[dict[str, object]]:
    """Serialize expected cells for a later reviewed read-only adapter."""
    cells = list(expected_cells() if rows is None else rows)
    return [asdict(row) for row in cells]