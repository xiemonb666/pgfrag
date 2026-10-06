"""Bounded symmetric Hausdorff distance for supplied 3-D polylines only."""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
from math import isqrt


MAX_VERTICES = 32
MIN_VERTICES = 2
MAX_COORDINATE_M = Fraction(1000)
MAX_INPUT_BITS = 256
SQRT_BITS = 40
SAMPLE_CAP = 128
UNIT_SCALE = {"m": Fraction(1), "mm": Fraction(1, 1000), "um": Fraction(1, 1_000_000)}
Point = tuple[Fraction, Fraction, Fraction]
Polyline = tuple[Point, ...]


class InputError(ValueError):
    """Raised when an input cannot be interpreted under this experiment contract."""


def _number(value: object, label: str) -> Fraction:
    if isinstance(value, bool):
        raise InputError(f"{label} must be a finite number, not bool")
    if isinstance(value, Fraction):
        result = value
    elif isinstance(value, Decimal):
        if not value.is_finite():
            raise InputError(f"{label} must be finite")
        result = Fraction(value)
    elif isinstance(value, int):
        result = Fraction(value)
    elif isinstance(value, float):
        if not math.isfinite(value):
            raise InputError(f"{label} must be finite")
        result = Fraction(Decimal(str(value)))
    else:
        raise InputError(f"{label} must be int, float, Decimal, or Fraction")
    if result.numerator.bit_length() > MAX_INPUT_BITS or result.denominator.bit_length() > MAX_INPUT_BITS:
        raise InputError(f"{label} exceeds the {MAX_INPUT_BITS}-bit exact-fraction guard")
    return result


def _coordinate(value: object, factor: Fraction, label: str) -> Fraction:
    result = _number(value, label) * factor
    if result.numerator.bit_length() > MAX_INPUT_BITS or result.denominator.bit_length() > MAX_INPUT_BITS:
        raise InputError(f"{label} exceeds the {MAX_INPUT_BITS}-bit exact-fraction guard after unit conversion")
    if abs(result) > MAX_COORDINATE_M:
        raise InputError(f"{label} exceeds 1000 m in magnitude")
    return result


def _polyline(raw: object, units: object, label: str) -> Polyline:
    if not isinstance(units, str) or units not in UNIT_SCALE:
        raise InputError(f"{label} units must be one of m, mm, um")
    if isinstance(raw, (str, bytes)) or not isinstance(raw, Sequence):
        raise InputError(f"{label} must be a sequence of 3-D vertices")
    if not MIN_VERTICES <= len(raw) <= MAX_VERTICES:
        raise InputError(f"{label} must contain 2..{MAX_VERTICES} vertices")
    factor = UNIT_SCALE[units]
    points: list[Point] = []
    for i, raw_point in enumerate(raw):
        if isinstance(raw_point, (str, bytes)) or not isinstance(raw_point, Sequence) or len(raw_point) != 3:
            raise InputError(f"{label}[{i}] must have exactly three coordinates")
        point = tuple(_coordinate(v, factor, f"{label}[{i}][{j}]") for j, v in enumerate(raw_point))
        if points and point == points[-1]:
            raise InputError(f"{label} contains a zero-length segment at vertex {i}")
        points.append(point)  # type: ignore[arg-type]
    return tuple(points)


def _dot(a: Point, b: Point) -> Fraction:
    return sum((x * y for x, y in zip(a, b)), Fraction(0))


def _squared_distance(a: Point, b: Point) -> Fraction:
    delta = tuple(x - y for x, y in zip(a, b))
    return _dot(delta, delta)  # type: ignore[arg-type]


def _point_segment_squared_distance(p: Point, a: Point, b: Point) -> Fraction:
    vector = tuple(y - x for x, y in zip(a, b))
    offset = tuple(x - y for x, y in zip(p, a))
    denominator = _dot(vector, vector)  # type: ignore[arg-type]
    if denominator == 0:
        raise ArithmeticError("validated target unexpectedly contains a zero-length segment")
    projection = _dot(offset, vector) / denominator  # type: ignore[arg-type]
    if projection <= 0:
        return _squared_distance(p, a)
    if projection >= 1:
        return _squared_distance(p, b)
    closest = tuple(a[i] + projection * vector[i] for i in range(3))  # type: ignore[index]
    return _squared_distance(p, closest)  # type: ignore[arg-type]


def _point_polyline_squared_distance(point: Point, target: Polyline) -> Fraction:
    return min(_point_segment_squared_distance(point, a, b) for a, b in zip(target, target[1:]))


def _sqrt_enclosure(squared: Fraction) -> tuple[Fraction, Fraction]:
    """Return dyadic lower/upper bounds of sqrt(squared), with width <= 2^-40."""
    if squared < 0:
        raise ArithmeticError("squared distance became negative")
    scaled_floor = (squared.numerator << (2 * SQRT_BITS)) // squared.denominator
    root_floor = isqrt(scaled_floor)
    lower = Fraction(root_floor, 1 << SQRT_BITS)
    if lower * lower == squared:
        return lower, lower
    return lower, lower + Fraction(1, 1 << SQRT_BITS)


def _distance_enclosure(point: Point, target: Polyline) -> tuple[Fraction, Fraction]:
    return _sqrt_enclosure(_point_polyline_squared_distance(point, target))


@dataclass(frozen=True)
class _Interval:
    left: Point
    right: Point
    left_distance_upper: Fraction
    right_distance_upper: Fraction


@dataclass
class _DirectedState:
    target: Polyline
    intervals: list[_Interval]
    lower: Fraction


def _interval_upper(interval: _Interval) -> Fraction:
    _, length_upper = _sqrt_enclosure(_squared_distance(interval.left, interval.right))
    return (interval.left_distance_upper + interval.right_distance_upper + length_upper) / 2


def _directed_state(source: Polyline, target: Polyline) -> _DirectedState:
    bounds = [_distance_enclosure(point, target) for point in source]
    intervals = [
        _Interval(source[i], source[i + 1], bounds[i][1], bounds[i + 1][1])
        for i in range(len(source) - 1)
    ]
    return _DirectedState(target, intervals, max(lower for lower, _ in bounds))


def _midpoint(a: Point, b: Point) -> Point:
    return tuple((x + y) / 2 for x, y in zip(a, b))  # type: ignore[return-value]


def _fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


@dataclass(frozen=True)
class HausdorffBounds:
    lower_m: Fraction
    upper_m: Fraction
    requested_atol_m: Fraction
    status: str
    samples: int
    sample_cap: int = SAMPLE_CAP
    geometry_scope: str = "GIVEN_PIECEWISE_LINEAR_POLYLINES_ONLY"
    physical_admission: str = "NOT_ADMITTED"

    def as_record(self) -> dict[str, object]:
        return {
            "lower_m": _fraction_text(self.lower_m),
            "upper_m": _fraction_text(self.upper_m),
            "requested_atol_m": _fraction_text(self.requested_atol_m),
            "status": self.status,
            "samples": self.samples,
            "sample_cap": self.sample_cap,
            "geometry_scope": self.geometry_scope,
            "physical_admission": self.physical_admission,
        }


def symmetric_polyline_hausdorff_bounds(
    polyline_a: object,
    polyline_b: object,
    *,
    units_a: object,
    units_b: object,
    atol_m: object,
) -> HausdorffBounds:
    """Bound Hausdorff distance between two supplied polylines in one shared frame."""
    a = _polyline(polyline_a, units_a, "polyline_a")
    b = _polyline(polyline_b, units_b, "polyline_b")
    tolerance = _number(atol_m, "atol_m")
    if tolerance <= 0 or tolerance > MAX_COORDINATE_M:
        raise InputError("atol_m must be positive and at most 1000 m")

    states = [_directed_state(a, b), _directed_state(b, a)]
    samples = sum(state.samples for state in states)
    while True:
        lower = max(state.lower for state in states)
        best: tuple[Fraction, int, int] | None = None
        for state_index, state in enumerate(states):
            for interval_index, interval in enumerate(state.intervals):
                candidate = (_interval_upper(interval), state_index, interval_index)
                if best is None or candidate[0] > best[0]:
                    best = candidate
        if best is None:
            raise ArithmeticError("validated polylines produced no source intervals")
        upper = best[0]
        gap = upper - lower
        if gap < 0:
            raise ArithmeticError("global Hausdorff bound gap became negative")
        if gap <= tolerance:
            status = "BOUNDS_WITHIN_REQUESTED_ATOL"
            break
        if samples >= SAMPLE_CAP:
            status = "UNKNOWN_WITH_BOUNDS"
            break

        _, state_index, interval_index = best
        state = states[state_index]
        interval = state.intervals.pop(interval_index)
        middle = _midpoint(interval.left, interval.right)
        middle_lower, middle_upper = _distance_enclosure(middle, state.target)
        state.lower = max(state.lower, middle_lower)
        state.intervals.extend(
            (
                _Interval(interval.left, middle, interval.left_distance_upper, middle_upper),
                _Interval(middle, interval.right, middle_upper, interval.right_distance_upper),
            )
        )
        samples += 1

    return HausdorffBounds(lower, upper, tolerance, status, samples)
