"""Bounded 3-D linear moving-least-squares interpolation primitive.

This module is a mathematical primitive for an explicit-surface mechanics
branch. It is not a fracture law, crack tracker, support-visibility model, or
CPM solver. Wendland C2 is the branch's selected engineering weight; it is not
asserted to be prescribed by Ai & Augarde.

The linear basis is evaluated in a fixed reference frame,
    p(X) = [1, (X - reference_origin) / reference_scale].
The scale is a single positive scalar so the primitive is covariant under
proper rotations. For a fixed source set and fixed visibility mask, all query
derivatives include the derivatives of the weighted moment matrix and its
linear solve. Visibility-mask changes are discrete; callers must mark queries
at/on a switching boundary as UNKNOWN.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

import numpy as np


MAX_SOURCES = 500
MAX_QUERIES = 500


class MLSStatus(str, Enum):
    OK = "OK"
    UNKNOWN_VISIBILITY_SWITCH = "UNKNOWN_VISIBILITY_SWITCH"
    UNKNOWN_SUPPORT_BOUNDARY = "UNKNOWN_SUPPORT_BOUNDARY"
    UNKNOWN_POOR_SUPPORT = "UNKNOWN_POOR_SUPPORT"
    UNKNOWN_DUPLICATE_SOURCES = "UNKNOWN_DUPLICATE_SOURCES"
    UNKNOWN_RANK_DEFICIENT = "UNKNOWN_RANK_DEFICIENT"
    UNKNOWN_ILL_CONDITIONED = "UNKNOWN_ILL_CONDITIONED"
    UNKNOWN_LINEAR_SOLVE_FAILED = "UNKNOWN_LINEAR_SOLVE_FAILED"


@dataclass(frozen=True)
class MLS3DConfig:
    """Fixed numerical scales for one material-reference coordinate system."""

    reference_origin: tuple[float, float, float]
    reference_scale: float
    support_radius: float
    max_condition_number: float = 1.0e12
    rank_relative_tolerance: float = 1.0e-12
    switch_tolerance: float = 1.0e-10
    duplicate_tolerance: float = 1.0e-12

    def __post_init__(self) -> None:
        origin = np.asarray(self.reference_origin, dtype=np.float64)
        if origin.shape != (3,) or not np.all(np.isfinite(origin)):
            raise ValueError("reference_origin must be three finite coordinates")
        object.__setattr__(self, "reference_origin", tuple(float(v) for v in origin))
        if not np.isfinite(self.reference_scale) or self.reference_scale <= 0.0:
            raise ValueError("reference_scale must be finite and positive")
        if not np.isfinite(self.support_radius) or self.support_radius <= 0.0:
            raise ValueError("support_radius must be finite and positive")
        if not np.isfinite(self.max_condition_number) or self.max_condition_number <= 1.0:
            raise ValueError("max_condition_number must be finite and greater than one")
        if not 0.0 < self.rank_relative_tolerance < 1.0:
            raise ValueError("rank_relative_tolerance must be in (0, 1)")
        if not 0.0 <= self.switch_tolerance < 0.5:
            raise ValueError("switch_tolerance must be in [0, 0.5)")
        if not 0.0 <= self.duplicate_tolerance < 1.0:
            raise ValueError("duplicate_tolerance must be in [0, 1)")


@dataclass(frozen=True)
class MLS3DResult:
    """Per-query MLS outputs; failed queries carry NaNs and an UNKNOWN status."""

    values: np.ndarray  # (n_query, n_field_components)
    gradients: np.ndarray  # (n_query, n_field_components, 3), d(value)/dX
    statuses: tuple[MLSStatus, ...]
    support_counts: np.ndarray
    condition_numbers: np.ndarray

    @property
    def ok(self) -> np.ndarray:
        """Read-only boolean success flags, one per query."""
        flags = np.fromiter((s is MLSStatus.OK for s in self.statuses), dtype=bool)
        flags.setflags(write=False)
        return flags


def wendland_c2(q: np.ndarray | float) -> tuple[np.ndarray, np.ndarray]:
    """Return Wendland C2 compact radial weight and d(weight)/d(q).

    w(q) = (1-q)^4 (4q+1) for 0 <= q < 1; it and its first two radial
    derivatives vanish continuously at q=1. Values outside support are zero.
    """
    q_arr = np.asarray(q, dtype=np.float64)
    if np.any(q_arr < 0.0) or not np.all(np.isfinite(q_arr)):
        raise ValueError("q must contain finite nonnegative values")
    inside = q_arr < 1.0
    one_minus = np.maximum(1.0 - q_arr, 0.0)
    weight = np.where(inside, one_minus**4 * (4.0 * q_arr + 1.0), 0.0)
    derivative = np.where(inside, -20.0 * q_arr * one_minus**3, 0.0)
    return weight, derivative


def _require_plain_ndarray(name: str, value: object) -> np.ndarray:
    # Exact type check avoids invoking user-defined __array__, __iter__, or
    # tolist hooks before bounded shape metadata has been inspected.
    if type(value) is not np.ndarray:
        raise ValueError(f"{name} must be a plain numpy.ndarray")
    return value


def _validate_real_finite(name: str, array: np.ndarray) -> np.ndarray:
    if array.dtype.kind not in "fiu":
        raise ValueError(f"{name} must use a real numeric dtype")
    # Callers validate dimensions and caps before scanning this array.
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite")
    return array


def _validate_points(name: str, points: np.ndarray, max_rows: int) -> np.ndarray:
    a = _require_plain_ndarray(name, points)
    if a.ndim != 2 or a.shape[1] != 3:
        raise ValueError(f"{name} must have shape (n, 3)")
    if a.shape[0] > max_rows:
        raise ValueError(f"budgetRejected: {name} has {a.shape[0]} rows; maximum is {max_rows}")
    return _validate_real_finite(name, a)


def _validate_values(values: np.ndarray, n_source: int) -> np.ndarray:
    y = _require_plain_ndarray("source_values", values)
    if y.ndim == 1:
        n_component = 1
    elif y.ndim == 2:
        n_component = y.shape[1]
    else:
        raise ValueError("source_values must have shape (n_source,) or (n_source, n_component)")
    if y.shape[0] > MAX_SOURCES:
        raise ValueError(f"budgetRejected: source_values has {y.shape[0]} rows; maximum is {MAX_SOURCES}")
    if not 1 <= n_component <= 3:
        raise ValueError(f"budgetRejected: source_values components must be between 1 and 3; got {n_component}")
    if y.shape[0] != n_source:
        raise ValueError("source_values row count must match source_positions")
    return _validate_real_finite("source_values", y)


def _finish_result(
    values: np.ndarray,
    gradients: np.ndarray,
    statuses: list[MLSStatus],
    support_counts: np.ndarray,
    condition_numbers: np.ndarray,
) -> MLS3DResult:
    for a in (values, gradients, support_counts, condition_numbers):
        a.setflags(write=False)
    return MLS3DResult(
        values=values,
        gradients=gradients,
        statuses=tuple(statuses),
        support_counts=support_counts,
        condition_numbers=condition_numbers,
    )


def evaluate_mls3d(
    source_positions: np.ndarray,
    source_values: np.ndarray,
    query_positions: np.ndarray,
    config: MLS3DConfig,
    *,
    visibility_mask: np.ndarray | None = None,
    mask_locally_constant: bool | np.ndarray | None = None,
) -> MLS3DResult:
    """Evaluate linear MLS and analytic spatial gradients at up to 500 queries.

    The boolean visibility mask is copied and treated as fixed for this call.
    When it is supplied, ``mask_locally_constant`` must explicitly be True for
    each query whose neighborhood is known to remain on the same side of any
    mask boundary. False/unspecified entries return UNKNOWN. This is a
    conditional derivative for an externally chosen mask; the mask is not a
    fracture or crack representation.

    A support-boundary query (q within ``config.switch_tolerance`` of 1 for an
    active source) is also returned UNKNOWN. Although this selected C2 kernel
    has a continuous first derivative there, this conservative rule avoids
    silently treating a discrete neighborhood switch as validated.

    Singular or ill-conditioned weighted moment matrices return UNKNOWN;
    no pseudo-inverse or regularization hides missing support. Inputs must be
    plain NumPy ndarrays. Invalid shapes, non-finite data, unsupported dtypes,
    and invalid configuration raise ValueError. Exceeding a source, query, or
    field-component cap raises ValueError with a ``budgetRejected`` marker
    before conversion or output allocation.
    """
    x_input = _validate_points("source_positions", source_positions, MAX_SOURCES)
    query_input = _validate_points("query_positions", query_positions, MAX_QUERIES)
    y_input = _validate_values(source_values, x_input.shape[0])

    n_source = x_input.shape[0]
    n_query = query_input.shape[0]
    if visibility_mask is None:
        if mask_locally_constant is not None:
            raise ValueError("mask_locally_constant requires visibility_mask")
        visible_input = None
        mask_stable_input = None
    else:
        visible_input = _require_plain_ndarray("visibility_mask", visibility_mask)
        if visible_input.shape != (n_source,) or visible_input.dtype.kind != "b":
            raise ValueError("visibility_mask must be a boolean vector with one entry per source")
        if mask_locally_constant is None:
            mask_stable_input = None
        elif isinstance(mask_locally_constant, (bool, np.bool_)):
            mask_stable_input = bool(mask_locally_constant)
        else:
            stable_input = _require_plain_ndarray("mask_locally_constant", mask_locally_constant)
            if stable_input.shape != (n_query,) or stable_input.dtype.kind != "b":
                raise ValueError("mask_locally_constant must be bool or a boolean vector with one entry per query")
            mask_stable_input = stable_input

    # Only after every input shape, hard cap, dtype, and finite value has been
    # validated do conversions or result-sized allocations occur.
    x = x_input.astype(np.float64, copy=False)
    queries = query_input.astype(np.float64, copy=False)
    y = y_input.astype(np.float64, copy=False)
    if y.ndim == 1:
        y = y[:, None]
    n_component = y.shape[1]
    values_out = np.full((n_query, n_component), np.nan, dtype=np.float64)
    gradients_out = np.full((n_query, n_component, 3), np.nan, dtype=np.float64)
    support_counts = np.zeros(n_query, dtype=np.int64)
    condition_numbers = np.full(n_query, np.nan, dtype=np.float64)
    if visible_input is None:
        visible = np.ones(n_source, dtype=bool)
        mask_stable = np.ones(n_query, dtype=bool)
    else:
        visible = np.array(visible_input, dtype=bool, copy=True)
        visible.setflags(write=False)
        if mask_stable_input is None:
            mask_stable = np.zeros(n_query, dtype=bool)
        elif isinstance(mask_stable_input, bool):
            mask_stable = np.full(n_query, mask_stable_input, dtype=bool)
        else:
            mask_stable = np.array(mask_stable_input, dtype=bool, copy=True)

    statuses = [MLSStatus.UNKNOWN_POOR_SUPPORT] * n_query

    origin = np.asarray(config.reference_origin, dtype=np.float64)
    scale = float(config.reference_scale)
    radius = float(config.support_radius)
    source_basis = np.empty((n_source, 4), dtype=np.float64)
    source_basis[:, 0] = 1.0
    source_basis[:, 1:] = (x - origin) / scale

    for qi, query in enumerate(queries):
        if not mask_stable[qi]:
            statuses[qi] = MLSStatus.UNKNOWN_VISIBILITY_SWITCH
            continue

        delta = x - query
        distance = np.linalg.norm(delta, axis=1)
        q = distance / radius
        support_counts[qi] = int(np.count_nonzero(visible & (q < 1.0)))
        boundary_indices = np.flatnonzero(visible & (np.abs(q - 1.0) <= config.switch_tolerance))
        if boundary_indices.size:
            statuses[qi] = MLSStatus.UNKNOWN_SUPPORT_BOUNDARY
            continue
        active_indices = np.flatnonzero(visible & (q < 1.0))
        if active_indices.size < 4:
            statuses[qi] = MLSStatus.UNKNOWN_POOR_SUPPORT
            continue

        local_positions = x[active_indices]
        duplicate_cutoff = config.duplicate_tolerance * radius
        local_delta = local_positions[:, None, :] - local_positions[None, :, :]
        squared = np.einsum("ijk,ijk->ij", local_delta, local_delta)
        np.fill_diagonal(squared, np.inf)
        if np.any(squared <= duplicate_cutoff * duplicate_cutoff):
            statuses[qi] = MLSStatus.UNKNOWN_DUPLICATE_SOURCES
            continue

        local_q = q[active_indices]
        weight, dweight_dq = wendland_c2(local_q)
        local_distance = distance[active_indices]
        local_delta = delta[active_indices]
        radial_derivative = np.zeros((active_indices.size, 3), dtype=np.float64)
        nonzero_distance = local_distance > 0.0
        # d q_i / d x = (x - X_i) / (radius * ||x-X_i||).
        radial_derivative[nonzero_distance] = (
            -local_delta[nonzero_distance]
            / (radius * local_distance[nonzero_distance, None])
        )
        dweight_dx = dweight_dq[:, None] * radial_derivative

        p = source_basis[active_indices]
        local_y = y[active_indices]
        moment = np.einsum("n,ni,nj->ij", weight, p, p, optimize=True)
        rhs = p.T @ (weight[:, None] * local_y)
        dmoment = np.einsum("nk,ni,nj->kij", dweight_dx, p, p, optimize=True)
        drhs = np.einsum("nk,ni,nc->kic", dweight_dx, p, local_y, optimize=True)

        try:
            singular_values = np.linalg.svd(moment, compute_uv=False)
        except np.linalg.LinAlgError:
            statuses[qi] = MLSStatus.UNKNOWN_LINEAR_SOLVE_FAILED
            continue
        largest = float(singular_values[0]) if singular_values.size else 0.0
        smallest = float(singular_values[-1]) if singular_values.size else 0.0
        if largest <= 0.0 or not np.isfinite(largest):
            statuses[qi] = MLSStatus.UNKNOWN_RANK_DEFICIENT
            continue
        rank = int(np.count_nonzero(singular_values > largest * config.rank_relative_tolerance))
        if rank < 4 or smallest <= 0.0:
            statuses[qi] = MLSStatus.UNKNOWN_RANK_DEFICIENT
            continue
        cond = largest / smallest
        condition_numbers[qi] = cond
        if not np.isfinite(cond) or cond > config.max_condition_number:
            statuses[qi] = MLSStatus.UNKNOWN_ILL_CONDITIONED
            continue

        try:
            coefficients = np.linalg.solve(moment, rhs)
            dcoefficients = np.empty((3, 4, n_component), dtype=np.float64)
            for axis in range(3):
                # Differentiate A a = b: A da/dx = db/dx - (dA/dx) a.
                dcoefficients[axis] = np.linalg.solve(
                    moment, drhs[axis] - dmoment[axis] @ coefficients
                )
        except np.linalg.LinAlgError:
            statuses[qi] = MLSStatus.UNKNOWN_LINEAR_SOLVE_FAILED
            continue

        query_basis = np.empty(4, dtype=np.float64)
        query_basis[0] = 1.0
        query_basis[1:] = (query - origin) / scale
        values_out[qi] = query_basis @ coefficients
        # The first term differentiates p(query); the second differentiates
        # coefficients through both the right-hand side and moment matrix.
        for axis in range(3):
            gradients_out[qi, :, axis] = (
                coefficients[axis + 1, :] / scale
                + query_basis @ dcoefficients[axis]
            )
        if not np.all(np.isfinite(values_out[qi])) or not np.all(np.isfinite(gradients_out[qi])):
            values_out[qi] = np.nan
            gradients_out[qi] = np.nan
            statuses[qi] = MLSStatus.UNKNOWN_LINEAR_SOLVE_FAILED
            continue
        statuses[qi] = MLSStatus.OK

    return _finish_result(values_out, gradients_out, statuses, support_counts, condition_numbers)


def green_lagrange_strain(displacement_gradient: np.ndarray) -> np.ndarray:
    """Return finite Green-Lagrange strain for a displacement gradient.

    ``displacement_gradient[..., i, j]`` is d(u_i)/d(X_j). This measure gives
    zero strain for any proper rigid rotation represented by the reproduced
    affine displacement field, unlike the infinitesimal strain at finite angle.
    """
    grad = np.asarray(displacement_gradient, dtype=np.float64)
    if grad.shape[-2:] != (3, 3) or not np.all(np.isfinite(grad)):
        raise ValueError("displacement_gradient must be finite with trailing shape (3, 3)")
    identity = np.eye(3, dtype=np.float64)
    deformation_gradient = grad + identity
    right_cauchy_green = np.einsum("...ki,...kj->...ij", deformation_gradient, deformation_gradient)
    return 0.5 * (right_cauchy_green - identity)

