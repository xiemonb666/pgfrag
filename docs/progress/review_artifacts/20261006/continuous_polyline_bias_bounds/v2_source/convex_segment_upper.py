"""Tighter certified interval upper bounds using convex target segments."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from continuous_polyline_bias_bounds import (
    MAX_COORDINATE_M,
    SAMPLE_CAP,
    HausdorffBounds,
    InputError,
    Point,
    Polyline,
    _midpoint,
    _number,
    _point_segment_squared_distance,
    _polyline,
    _sqrt_enclosure,
    _squared_distance,
)


@dataclass(frozen=True)
class _Sample:
    point: Point
    segment_bounds: tuple[tuple[Fraction, Fraction], ...]
    distance_lower: Fraction
    distance_upper: Fraction


@dataclass(frozen=True)
class _Interval:
    left: _Sample
    right: _Sample


@dataclass
class _DirectedState:
    target: Polyline
    intervals: list[_Interval]
    lower: Fraction


def _sample(point: Point, target: Polyline) -> _Sample:
    bounds = tuple(
        _sqrt_enclosure(_point_segment_squared_distance(point, a, b))
        for a, b in zip(target, target[1:])
    )
    return _Sample(
        point,
        bounds,
        min(lower for lower, _ in bounds),
        min(upper for _, upper in bounds),
    )


def _interval_upper(interval: _Interval) -> Fraction:
    left, right = interval.left, interval.right
    _, length_upper = _sqrt_enclosure(_squared_distance(left.point, right.point))
    lipschitz_upper = (left.distance_upper + right.distance_upper + length_upper) / 2
    convex_upper = min(
        max(left.segment_bounds[index][1], right.segment_bounds[index][1])
        for index in range(len(left.segment_bounds))
    )
    return min(lipschitz_upper, convex_upper)


def _directed_state(source: Polyline, target: Polyline) -> _DirectedState:
    samples = [_sample(point, target) for point in source]
    intervals = [_Interval(samples[i], samples[i + 1]) for i in range(len(samples) - 1)]
    return _DirectedState(target, intervals, max(sample.distance_lower for sample in samples))


def symmetric_polyline_hausdorff_bounds_convex(
    polyline_a: object,
    polyline_b: object,
    *,
    units_a: object,
    units_b: object,
    atol_m: object,
) -> HausdorffBounds:
    """Bound symmetric Hausdorff distance for supplied 3-D polylines only."""
    a = _polyline(polyline_a, units_a, "polyline_a")
    b = _polyline(polyline_b, units_b, "polyline_b")
    tolerance = _number(atol_m, "atol_m")
    if tolerance <= 0 or tolerance > MAX_COORDINATE_M:
        raise InputError("atol_m must be positive and at most 1000 m")

    states = [_directed_state(a, b), _directed_state(b, a)]
    samples = len(a) + len(b)
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
        middle = _sample(_midpoint(interval.left.point, interval.right.point), state.target)
        state.lower = max(state.lower, middle.distance_lower)
        state.intervals.extend((_Interval(interval.left, middle), _Interval(middle, interval.right)))
        samples += 1

    return HausdorffBounds(lower, upper, tolerance, status, samples)
