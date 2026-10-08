"""Exact, bounded local-stability certificates for finite triangular facets.

This is a source-only geometry primitive. For fixed float64 source and facet
coordinates, it uses exact rational predicates on the dyadic values represented
by those floats to certify whether a source-to-query segment stays visible or
blocked in some open neighbourhood of the query. It does not provide a
neighbourhood radius, an IEEE-roundoff guarantee, a curved-surface guarantee,
or physical/solver truth. A certificate is accepted only when every facet is
strictly decidable and its pair status agrees with the pinned float visibility
primitive. UNKNOWN never authorizes query-local support.

Coordinates are SI metres. Public inputs must be exact, C-contiguous float64
NumPy ndarrays and are shape/byte/work bounded before coordinate conversion or
Fraction construction. Fractions are limited to 256 input bits and 2048
intermediate bits. The ideal-real open-neighbourhood statement is about the
exact dyadic coordinates supplied to this function; a rotation, unit
conversion, or other float transformation generally creates different dyadic
inputs and is not covered by an exact-invariance claim.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import hashlib
import importlib.util
import math
from math import gcd
from pathlib import Path
import sys
from types import ModuleType

import numpy as np


MAX_SOURCES = 500
MAX_QUERIES = 500
MAX_TRIANGLES = 256
MAX_QUERY_SOURCE_PAIRS = 16_384
MAX_EXACT_SEGMENT_TRIANGLE_PAIRS = 16_384
MAX_INPUT_FRACTION_BITS = 256
MAX_INTERMEDIATE_FRACTION_BITS = 2_048
MAX_ABS_COORDINATE_M = 1.0e150

_FROZEN_CLASSIFIER_SHA256 = "9d5097ca33ce899cf788ad43ef0740a34296216a8fd7aa490868c27ffe1ebb03"
_FROZEN_CLASSIFIER_MAX_BYTES = 256_000
_FROZEN_REFERENCE_NOT_EVALUATED = "NOT_EVALUATED"

_VISIBLE = "VISIBLE"
_BLOCKED = "BLOCKED"
_UNKNOWN = "UNKNOWN"
_STABLE_VISIBLE = "STABLE_VISIBLE"
_STABLE_BLOCKED = "STABLE_BLOCKED"
_FACET_UNKNOWN = "UNKNOWN"


@dataclass(frozen=True, slots=True)
class StabilityResult:
    """Immutable query-major stability results and audit detail.

    Pair arrays have shape ``(n_query, n_source)``. Facet arrays, when
    present, have shape ``(n_query, n_source, n_triangle)``. A true
    ``stable_visible_mask`` is the only result that may authorize visible
    query-local support. True means an ideal-real open-neighbourhood
    existence certificate for the exact supplied dyadic inputs, and agreement
    with the pinned float classifier; it does not provide a radius.
    """

    pair_statuses: np.ndarray
    pair_reasons: np.ndarray
    stable_visible_mask: np.ndarray
    stable_blocked_mask: np.ndarray
    unknown_mask: np.ndarray
    frozen_float_statuses: np.ndarray
    facet_statuses: np.ndarray | None
    facet_reasons: np.ndarray | None
    fraction_input_max_bits: int
    fraction_intermediate_max_bits: int
    certificate_scope: str = "IDEAL_REAL_QUERY_OPEN_NEIGHBORHOOD_EXISTENCE_ONLY_NO_RADIUS"
    ieee_robustness: str = "NOT_CLAIMED"
    spatial_stability: str = "CERTIFIED_ONLY_WHERE_STABLE_MASK_IS_TRUE"
    mechanics_coupling: str = "NOT_COUPLED"


class _IntermediateBudgetExceeded(Exception):
    pass


class _BudgetArithmetic:
    """Fraction operators that reject oversized products before constructing them."""

    def __init__(self) -> None:
        self.maximum_bits = 0

    @staticmethod
    def _bits(value: Fraction) -> int:
        return max(abs(value.numerator).bit_length(), value.denominator.bit_length())

    def _check(self, *values: Fraction, predicted: int | None = None) -> None:
        observed = max((self._bits(value) for value in values), default=0)
        self.maximum_bits = max(self.maximum_bits, observed)
        if observed > MAX_INTERMEDIATE_FRACTION_BITS:
            raise _IntermediateBudgetExceeded
        if predicted is not None and predicted > MAX_INTERMEDIATE_FRACTION_BITS:
            raise _IntermediateBudgetExceeded

    def add(self, left: Fraction, right: Fraction) -> Fraction:
        common = gcd(left.denominator, right.denominator)
        left_scale = right.denominator // common
        right_scale = left.denominator // common
        predicted_numerator = max(
            abs(left.numerator).bit_length() + left_scale.bit_length(),
            abs(right.numerator).bit_length() + right_scale.bit_length(),
        ) + 1
        predicted = max(
            predicted_numerator,
            left.denominator.bit_length() + left_scale.bit_length(),
        )
        self._check(left, right, predicted=predicted)
        numerator = left.numerator * left_scale + right.numerator * right_scale
        denominator = left.denominator * left_scale
        result = Fraction(numerator, denominator)
        self._check(result)
        return result

    def sub(self, left: Fraction, right: Fraction) -> Fraction:
        common = gcd(left.denominator, right.denominator)
        left_scale = right.denominator // common
        right_scale = left.denominator // common
        predicted_numerator = max(
            abs(left.numerator).bit_length() + left_scale.bit_length(),
            abs(right.numerator).bit_length() + right_scale.bit_length(),
        ) + 1
        predicted = max(
            predicted_numerator,
            left.denominator.bit_length() + left_scale.bit_length(),
        )
        self._check(left, right, predicted=predicted)
        numerator = left.numerator * left_scale - right.numerator * right_scale
        denominator = left.denominator * left_scale
        result = Fraction(numerator, denominator)
        self._check(result)
        return result

    def mul(self, left: Fraction, right: Fraction) -> Fraction:
        cancel_left = gcd(abs(left.numerator), right.denominator)
        cancel_right = gcd(abs(right.numerator), left.denominator)
        left_numerator = left.numerator // cancel_left
        right_denominator = right.denominator // cancel_left
        right_numerator = right.numerator // cancel_right
        left_denominator = left.denominator // cancel_right
        predicted = max(
            abs(left_numerator).bit_length() + abs(right_numerator).bit_length(),
            left_denominator.bit_length() + right_denominator.bit_length(),
        )
        self._check(left, right, predicted=predicted)
        result = Fraction(left_numerator * right_numerator, left_denominator * right_denominator)
        self._check(result)
        return result

    def div(self, left: Fraction, right: Fraction) -> Fraction:
        if right == 0:
            raise ZeroDivisionError
        cancel_numerator = gcd(abs(left.numerator), abs(right.numerator))
        cancel_denominator = gcd(left.denominator, right.denominator)
        left_numerator = left.numerator // cancel_numerator
        right_numerator = right.numerator // cancel_numerator
        left_denominator = left.denominator // cancel_denominator
        right_denominator = right.denominator // cancel_denominator
        predicted = max(
            abs(left_numerator).bit_length() + right_denominator.bit_length(),
            left_denominator.bit_length() + abs(right_numerator).bit_length(),
        )
        self._check(left, right, predicted=predicted)
        result = Fraction(left_numerator * right_denominator, left_denominator * right_numerator)
        self._check(result)
        return result


def _freeze_array(array: np.ndarray) -> np.ndarray:
    """Return an immutable ndarray whose storage is an immutable bytes object."""
    frozen_bytes = array.tobytes(order="C")
    return np.frombuffer(frozen_bytes, dtype=array.dtype).reshape(array.shape)


def _preflight_array(
    name: str,
    value: object,
    *,
    tail_shape: tuple[int, ...],
    max_rows: int,
    max_bytes: int,
) -> np.ndarray:
    """Check ndarray metadata and budgets without invoking user array hooks."""
    if type(value) is not np.ndarray:
        raise ValueError(f"{name} must be an exact numpy.ndarray")
    if value.ndim != len(tail_shape) + 1 or tuple(value.shape[1:]) != tail_shape:
        tail = ", ".join(str(part) for part in tail_shape)
        raise ValueError(f"{name} must have shape (n, {tail})")
    rows = int(value.shape[0])
    if rows > max_rows:
        raise ValueError(f"budgetRejected: {name} has {rows} rows; maximum is {max_rows}")
    if value.dtype != np.dtype(np.float64):
        raise ValueError(f"{name} must use float64 coordinates")
    if value.nbytes > max_bytes:
        raise ValueError(f"budgetRejected: {name} is {value.nbytes} bytes; maximum is {max_bytes}")
    if not value.flags.c_contiguous:
        raise ValueError(f"{name} must be C-contiguous")
    return value


def _validated_ratios(
    named_arrays: tuple[tuple[str, np.ndarray], ...],
) -> tuple[int, float]:
    """Validate finite SI coordinates and report max exact-ratio bit length."""
    maximum_bits = 0
    maximum_coordinate = 0.0
    for name, array in named_arrays:
        for scalar in array.flat:
            value = float(scalar)
            if not math.isfinite(value):
                raise ValueError(f"{name} must contain only finite coordinates")
            absolute = abs(value)
            maximum_coordinate = max(maximum_coordinate, absolute)
            numerator, denominator = value.as_integer_ratio()
            maximum_bits = max(maximum_bits, abs(numerator).bit_length(), denominator.bit_length())
    if maximum_coordinate > MAX_ABS_COORDINATE_M:
        raise ValueError(f"coordinate magnitude exceeds safe float64 range {MAX_ABS_COORDINATE_M:g} m")
    return maximum_bits, maximum_coordinate


def _as_fraction(value: np.float64) -> Fraction:
    numerator, denominator = float(value).as_integer_ratio()
    return Fraction(numerator, denominator)


def _fraction_points(array: np.ndarray) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
    return tuple(tuple(_as_fraction(value) for value in row) for row in array)  # type: ignore[return-value]


def _vec_sub(
    arithmetic: _BudgetArithmetic,
    left: tuple[Fraction, Fraction, Fraction],
    right: tuple[Fraction, Fraction, Fraction],
) -> tuple[Fraction, Fraction, Fraction]:
    return tuple(arithmetic.sub(a, b) for a, b in zip(left, right))  # type: ignore[return-value]


def _dot(
    arithmetic: _BudgetArithmetic,
    left: tuple[Fraction, Fraction, Fraction],
    right: tuple[Fraction, Fraction, Fraction],
) -> Fraction:
    products = [arithmetic.mul(a, b) for a, b in zip(left, right)]
    return arithmetic.add(arithmetic.add(products[0], products[1]), products[2])


def _cross(
    arithmetic: _BudgetArithmetic,
    left: tuple[Fraction, Fraction, Fraction],
    right: tuple[Fraction, Fraction, Fraction],
) -> tuple[Fraction, Fraction, Fraction]:
    return (
        arithmetic.sub(arithmetic.mul(left[1], right[2]), arithmetic.mul(left[2], right[1])),
        arithmetic.sub(arithmetic.mul(left[2], right[0]), arithmetic.mul(left[0], right[2])),
        arithmetic.sub(arithmetic.mul(left[0], right[1]), arithmetic.mul(left[1], right[0])),
    )


def _exact_barycentric(
    arithmetic: _BudgetArithmetic,
    point: tuple[Fraction, Fraction, Fraction],
    triangle: tuple[
        tuple[Fraction, Fraction, Fraction],
        tuple[Fraction, Fraction, Fraction],
        tuple[Fraction, Fraction, Fraction],
    ],
) -> tuple[Fraction, Fraction, Fraction] | None:
    """Return exact barycentric coordinates, or None for a degenerate basis."""
    a, b, c = triangle
    ab = _vec_sub(arithmetic, b, a)
    ac = _vec_sub(arithmetic, c, a)
    ap = _vec_sub(arithmetic, point, a)
    dot00 = _dot(arithmetic, ab, ab)
    dot01 = _dot(arithmetic, ab, ac)
    dot02 = _dot(arithmetic, ab, ap)
    dot11 = _dot(arithmetic, ac, ac)
    dot12 = _dot(arithmetic, ac, ap)
    denominator = arithmetic.sub(arithmetic.mul(dot00, dot11), arithmetic.mul(dot01, dot01))
    if denominator <= 0:
        return None
    u_numerator = arithmetic.sub(arithmetic.mul(dot11, dot02), arithmetic.mul(dot01, dot12))
    v_numerator = arithmetic.sub(arithmetic.mul(dot00, dot12), arithmetic.mul(dot01, dot02))
    u = arithmetic.div(u_numerator, denominator)
    v = arithmetic.div(v_numerator, denominator)
    w = arithmetic.sub(arithmetic.sub(Fraction(1), u), v)
    return w, u, v


def _facet_certificate(
    source: tuple[Fraction, Fraction, Fraction],
    query: tuple[Fraction, Fraction, Fraction],
    triangle: tuple[
        tuple[Fraction, Fraction, Fraction],
        tuple[Fraction, Fraction, Fraction],
        tuple[Fraction, Fraction, Fraction],
    ],
) -> tuple[str, str, int]:
    """Return exact local facet status, reason, and maximum intermediate bits."""
    arithmetic = _BudgetArithmetic()
    try:
        a, b, c = triangle
        ab = _vec_sub(arithmetic, b, a)
        ac = _vec_sub(arithmetic, c, a)
        normal = _cross(arithmetic, ab, ac)
        if normal == (0, 0, 0):
            return _FACET_UNKNOWN, "UNKNOWN_DEGENERATE_TRIANGLE", arithmetic.maximum_bits

        source_distance = _dot(arithmetic, normal, _vec_sub(arithmetic, source, a))
        query_distance = _dot(arithmetic, normal, _vec_sub(arithmetic, query, a))

        if source_distance == 0 and query_distance == 0:
            return _FACET_UNKNOWN, "UNKNOWN_COPLANAR_SEGMENT", arithmetic.maximum_bits

        if source_distance == 0:
            barycentric = _exact_barycentric(arithmetic, source, triangle)
            if barycentric is None:
                return _FACET_UNKNOWN, "UNKNOWN_DEGENERATE_TRIANGLE", arithmetic.maximum_bits
            if any(value < 0 for value in barycentric):
                return _STABLE_VISIBLE, "STABLE_SOURCE_PLANE_FINITE_MISS", arithmetic.maximum_bits
            return _FACET_UNKNOWN, "UNKNOWN_SOURCE_FRONT_OR_EDGE", arithmetic.maximum_bits

        if query_distance == 0:
            barycentric = _exact_barycentric(arithmetic, query, triangle)
            if barycentric is None:
                return _FACET_UNKNOWN, "UNKNOWN_DEGENERATE_TRIANGLE", arithmetic.maximum_bits
            if any(value < 0 for value in barycentric):
                return _STABLE_VISIBLE, "STABLE_QUERY_PLANE_FINITE_MISS", arithmetic.maximum_bits
            return _FACET_UNKNOWN, "UNKNOWN_QUERY_FRONT_OR_EDGE", arithmetic.maximum_bits

        if (source_distance > 0 and query_distance > 0) or (
            source_distance < 0 and query_distance < 0
        ):
            return _STABLE_VISIBLE, "STABLE_SAME_STRICT_SIDE", arithmetic.maximum_bits

        parameter = arithmetic.div(source_distance, arithmetic.sub(source_distance, query_distance))
        direction = _vec_sub(arithmetic, query, source)
        intersection = tuple(
            arithmetic.add(source[index], arithmetic.mul(parameter, direction[index]))
            for index in range(3)
        )
        barycentric = _exact_barycentric(arithmetic, intersection, triangle)
        if barycentric is None:
            return _FACET_UNKNOWN, "UNKNOWN_DEGENERATE_TRIANGLE", arithmetic.maximum_bits
        if any(value < 0 for value in barycentric):
            return _STABLE_VISIBLE, "STABLE_STRICT_FINITE_MISS", arithmetic.maximum_bits
        if all(value > 0 for value in barycentric):
            return _STABLE_BLOCKED, "STABLE_STRICT_INTERIOR_CROSSING", arithmetic.maximum_bits
        return _FACET_UNKNOWN, "UNKNOWN_TRIANGLE_EDGE_OR_FRONT", arithmetic.maximum_bits
    except _IntermediateBudgetExceeded:
        return _FACET_UNKNOWN, "UNKNOWN_INTERMEDIATE_BIT_BUDGET", arithmetic.maximum_bits
    except ZeroDivisionError:
        return _FACET_UNKNOWN, "UNKNOWN_EXACT_PREDICATE_UNDECIDED", arithmetic.maximum_bits


def _load_frozen_classifier() -> tuple[ModuleType | None, str]:
    """Load the hash-pinned visibility primitive, failing closed on drift."""
    source_path = Path(__file__).resolve().parent.parent / "explicit_visibility_20261005" / "visibility3d.py"
    try:
        if not source_path.is_file() or source_path.stat().st_size > _FROZEN_CLASSIFIER_MAX_BYTES:
            return None, "UNKNOWN_REFERENCE_SOURCE_UNAVAILABLE"
        source_bytes = source_path.read_bytes()
        if hashlib.sha256(source_bytes).hexdigest() != _FROZEN_CLASSIFIER_SHA256:
            return None, "UNKNOWN_REFERENCE_HASH_MISMATCH"
        module_name = "_explicit_visibility_20261005_frozen_for_stability"
        specification = importlib.util.spec_from_file_location(module_name, source_path)
        if specification is None or specification.loader is None:
            return None, "UNKNOWN_REFERENCE_SOURCE_UNAVAILABLE"
        module = importlib.util.module_from_spec(specification)
        sys.modules[module_name] = module
        specification.loader.exec_module(module)
        return module, ""
    except Exception:
        return None, "UNKNOWN_REFERENCE_LOAD_ERROR"


def _frozen_float_statuses(
    reference: ModuleType | None,
    sources: np.ndarray,
    queries: np.ndarray,
    triangles: np.ndarray,
    failure_reason: str,
) -> tuple[np.ndarray, str]:
    if reference is None:
        return np.full((queries.shape[0], sources.shape[0]), _FROZEN_REFERENCE_NOT_EVALUATED, dtype="<U13"), failure_reason
    try:
        result = reference.classify_visibility(sources, queries, triangles)
        return np.asarray(result.statuses, dtype="<U13"), ""
    except Exception:
        return np.full((queries.shape[0], sources.shape[0]), _FROZEN_REFERENCE_NOT_EVALUATED, dtype="<U13"), "UNKNOWN_REFERENCE_CLASSIFICATION_ERROR"


def certify_local_visibility_stability(
    source_positions: np.ndarray,
    query_positions: np.ndarray,
    triangle_facets: np.ndarray,
) -> StabilityResult:
    """Certify strict query-local stability of finite crack-facet visibility.

    For one fixed source and one nondegenerate facet, exact signed plane
    distances with the same strict sign prove a local finite-facet miss. Strict
    opposite signs produce an exact plane intersection: strictly positive
    barycentric coordinates prove stable interior blocking; any strictly
    negative coordinate proves a finite-facet miss; an edge/front point is
    UNKNOWN. A source or query exactly on the plane is UNKNOWN when it lies in
    or on the finite triangle, but a strict barycentric miss proves that the
    endpoint is separated from that facet and is stable visible. Coplanar
    segments and degenerate facets are UNKNOWN.

    Aggregation is deliberately conservative: every facet must have a strict
    certificate. Any ambiguous facet keeps the pair UNKNOWN even if another
    facet is a definite blocker; per-facet states preserve both facts. The
    aggregate candidate must also match a fresh result from the hash-pinned
    float classifier before a stable mask is published.
    """
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
    pair_count = n_source * n_query
    exact_work = pair_count * n_triangle
    if pair_count > MAX_QUERY_SOURCE_PAIRS or exact_work > MAX_EXACT_SEGMENT_TRIANGLE_PAIRS:
        pair_statuses = np.full((n_query, n_source), _UNKNOWN, dtype="<U7")
        budget_reason = "UNKNOWN_PAIR_BUDGET" if pair_count > MAX_QUERY_SOURCE_PAIRS else "UNKNOWN_WORK_BUDGET"
        pair_reasons = np.full((n_query, n_source), budget_reason, dtype="<U48")
        frozen = np.full((n_query, n_source), _FROZEN_REFERENCE_NOT_EVALUATED, dtype="<U13")
        return StabilityResult(
            pair_statuses=_freeze_array(pair_statuses),
            pair_reasons=_freeze_array(pair_reasons),
            stable_visible_mask=_freeze_array(np.zeros((n_query, n_source), dtype=bool)),
            stable_blocked_mask=_freeze_array(np.zeros((n_query, n_source), dtype=bool)),
            unknown_mask=_freeze_array(np.ones((n_query, n_source), dtype=bool)),
            frozen_float_statuses=_freeze_array(frozen),
            facet_statuses=None,
            facet_reasons=None,
            fraction_input_max_bits=0,
            fraction_intermediate_max_bits=0,
        )

    # All metadata budgets have passed. Take bounded private snapshots before
    # scanning values so the reference classifier and exact predicates see the
    # same coordinates even if a caller later mutates its input arrays.
    sources_float = source_input.astype(np.float64, copy=True)
    queries_float = query_input.astype(np.float64, copy=True)
    triangles_float = triangle_input.astype(np.float64, copy=True)
    fraction_input_max_bits, _maximum_coordinate = _validated_ratios(
        (
            ("source_positions", sources_float),
            ("query_positions", queries_float),
            ("triangle_facets", triangles_float),
        )
    )
    reference, reference_failure = _load_frozen_classifier()
    frozen_statuses, classification_failure = _frozen_float_statuses(
        reference,
        sources_float,
        queries_float,
        triangles_float,
        reference_failure,
    )
    if classification_failure:
        reference = None

    pair_statuses = np.full((n_query, n_source), _UNKNOWN, dtype="<U7")
    pair_reasons = np.full((n_query, n_source), "UNKNOWN_FACET_AMBIGUITY", dtype="<U48")
    facet_statuses = np.full((n_query, n_source, n_triangle), _FACET_UNKNOWN, dtype="<U20")
    facet_reasons = np.full((n_query, n_source, n_triangle), "UNKNOWN_INPUT_BIT_BUDGET", dtype="<U48")
    maximum_intermediate_bits = 0

    if fraction_input_max_bits <= MAX_INPUT_FRACTION_BITS:
        sources = _fraction_points(sources_float)
        queries = _fraction_points(queries_float)
        triangles = _fraction_points(triangles_float)
        for query_index, query in enumerate(queries):
            for source_index, source in enumerate(sources):
                for triangle_index, triangle in enumerate(triangles):
                    status, reason, observed_bits = _facet_certificate(source, query, triangle)
                    facet_statuses[query_index, source_index, triangle_index] = status
                    facet_reasons[query_index, source_index, triangle_index] = reason
                    maximum_intermediate_bits = max(maximum_intermediate_bits, observed_bits)
    elif n_triangle == 0:
        facet_reasons = np.empty((n_query, n_source, 0), dtype="<U48")
    else:
        facet_reasons.fill("UNKNOWN_INPUT_BIT_BUDGET")

    for query_index in range(n_query):
        for source_index in range(n_source):
            reasons = facet_reasons[query_index, source_index]
            statuses = facet_statuses[query_index, source_index]
            if fraction_input_max_bits > MAX_INPUT_FRACTION_BITS and n_triangle:
                candidate = _UNKNOWN
                candidate_reason = "UNKNOWN_INPUT_BIT_BUDGET"
            elif any(status == _FACET_UNKNOWN for status in statuses):
                candidate = _UNKNOWN
                candidate_reason = "UNKNOWN_FACET_AMBIGUITY"
            elif any(status == _STABLE_BLOCKED for status in statuses):
                candidate = _BLOCKED
                candidate_reason = "ALL_FACETS_STRICT_BLOCKED_OR_MISS"
            else:
                candidate = _VISIBLE
                candidate_reason = "ALL_FACETS_STRICT_VISIBLE"

            if candidate == _UNKNOWN:
                pair_statuses[query_index, source_index] = _UNKNOWN
                if candidate_reason == "UNKNOWN_FACET_AMBIGUITY" and len(reasons):
                    # Retain a useful first facet reason while all details remain
                    # available in the complete facet arrays.
                    unknown_indices = [
                        index for index, status in enumerate(statuses) if status == _FACET_UNKNOWN
                    ]
                    if unknown_indices:
                        pair_reasons[query_index, source_index] = str(reasons[unknown_indices[0]])
                    else:
                        pair_reasons[query_index, source_index] = candidate_reason
                else:
                    pair_reasons[query_index, source_index] = candidate_reason
                continue

            if reference is None:
                pair_statuses[query_index, source_index] = _UNKNOWN
                pair_reasons[query_index, source_index] = classification_failure or reference_failure or "UNKNOWN_REFERENCE_UNAVAILABLE"
                continue

            if str(frozen_statuses[query_index, source_index]) != candidate:
                pair_statuses[query_index, source_index] = _UNKNOWN
                pair_reasons[query_index, source_index] = "UNKNOWN_FLOAT_STATUS_MISMATCH"
                continue

            pair_statuses[query_index, source_index] = candidate
            pair_reasons[query_index, source_index] = candidate_reason

    stable_visible = pair_statuses == _VISIBLE
    stable_blocked = pair_statuses == _BLOCKED
    unknown = pair_statuses == _UNKNOWN
    return StabilityResult(
        pair_statuses=_freeze_array(pair_statuses),
        pair_reasons=_freeze_array(pair_reasons),
        stable_visible_mask=_freeze_array(stable_visible),
        stable_blocked_mask=_freeze_array(stable_blocked),
        unknown_mask=_freeze_array(unknown),
        frozen_float_statuses=_freeze_array(frozen_statuses),
        facet_statuses=_freeze_array(facet_statuses),
        facet_reasons=_freeze_array(facet_reasons),
        fraction_input_max_bits=fraction_input_max_bits,
        fraction_intermediate_max_bits=maximum_intermediate_bits,
    )

