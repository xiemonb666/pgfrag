"""Bounded 3-D MLS shape functions for EFGM-style assembly.

This source-only API returns the MLS nodal shape row directly, without
constructing an N-by-N identity field and routing it through the scalar MLS
interpolator. It reuses the reviewed phase1 configuration, kernel, status, and
bounded position/mask validation contract. This module remains
MATHEMATICAL_PRIMITIVE_ONLY / NOT_COUPLED; it contains no mechanics solver or
fracture state.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from mls3d_mls import (
    MAX_QUERIES,
    MAX_SOURCES,
    MLS3DConfig,
    MLSStatus,
    _require_plain_ndarray,
    _validate_points,
    wendland_c2,
)


@dataclass(frozen=True)
class ShapeFunctions3DResult:
    """Batched MLS shape data with standard query/node/axis orientation."""

    weights: np.ndarray  # (n_query, n_source), Phi[query, source]
    gradients: np.ndarray  # (n_query, n_source, 3), dPhi[query, source]/dX[axis]
    statuses: tuple[MLSStatus, ...]
    support_counts: np.ndarray  # (n_query,)
    condition_numbers: np.ndarray  # (n_query,)

    @property
    def reasons(self) -> tuple[str | None, ...]:
        """A stable textual reason for UNKNOWN rows; None for successful rows."""
        return tuple(None if status is MLSStatus.OK else status.value for status in self.statuses)

    @property
    def ok(self) -> np.ndarray:
        flags = np.fromiter((status is MLSStatus.OK for status in self.statuses), dtype=bool)
        return _immutable_array(flags)


def _fixed_visibility_inputs(
    n_source: int,
    n_query: int,
    visibility_mask: np.ndarray | None,
    mask_locally_constant: bool | np.ndarray | None,
) -> tuple[np.ndarray | None, bool | np.ndarray | None]:
    if visibility_mask is None:
        if mask_locally_constant is not None:
            raise ValueError("mask_locally_constant requires visibility_mask")
        return None, None

    mask = _require_plain_ndarray("visibility_mask", visibility_mask)
    if mask.shape != (n_source,) or mask.dtype.kind != "b":
        raise ValueError("visibility_mask must be a boolean vector with one entry per source")
    if mask_locally_constant is None:
        stable: bool | np.ndarray | None = None
    elif isinstance(mask_locally_constant, (bool, np.bool_)):
        stable = bool(mask_locally_constant)
    else:
        stable_array = _require_plain_ndarray("mask_locally_constant", mask_locally_constant)
        if stable_array.shape != (n_query,) or stable_array.dtype.kind != "b":
            raise ValueError("mask_locally_constant must be bool or a boolean vector with one entry per query")
        stable = stable_array
    return mask, stable


def _immutable_array(array: np.ndarray) -> np.ndarray:
    """Copy into immutable bytes-backed storage that cannot be made writable."""
    return np.frombuffer(array.tobytes(order="C"), dtype=array.dtype).reshape(array.shape)


def _immutable_result(
    weights: np.ndarray,
    gradients: np.ndarray,
    statuses: list[MLSStatus],
    support_counts: np.ndarray,
    condition_numbers: np.ndarray,
) -> ShapeFunctions3DResult:
    return ShapeFunctions3DResult(
        weights=_immutable_array(weights),
        gradients=_immutable_array(gradients),
        statuses=tuple(statuses),
        support_counts=_immutable_array(support_counts),
        condition_numbers=_immutable_array(condition_numbers),
    )


def shape_functions_3d(
    source_positions: np.ndarray,
    query_positions: np.ndarray,
    config: MLS3DConfig,
    *,
    visibility_mask: np.ndarray | None = None,
    mask_locally_constant: bool | np.ndarray | None = None,
) -> ShapeFunctions3DResult:
    """Return compact-support linear MLS shape values and analytic gradients.

    Result orientations are `weights[query, source]` and
    `gradients[query, source, reference_axis] = d(Phi_source)/dX_axis`.
    Successful rows contain zeros for masked or out-of-support sources;
    UNKNOWN rows contain NaNs for all source entries. The source and query
    limits are 500 each. Inputs must be plain NumPy arrays, and all shape,
    finite-value, mask, and budget checks happen before conversion or output
    allocation. Cap violations raise `ValueError` with `budgetRejected`.

    The dual solve is `A.T c=p(x)`, with `Phi_i=w_i p_i.T c`. Differentiation
    uses `A.T dc/dx_j=dp/dx_j-(dA.T/dx_j)c`, then differentiates both `w_i`
    and `c`; it does not form an identity matrix or perform N independent MLS
    interpolations. Rank, condition, support, radial-switch, and fixed-mask
    guards match the phase1 scalar/vector primitive.
    """
    x_input = _validate_points("source_positions", source_positions, MAX_SOURCES)
    query_input = _validate_points("query_positions", query_positions, MAX_QUERIES)
    n_source = x_input.shape[0]
    n_query = query_input.shape[0]
    mask_input, stable_input = _fixed_visibility_inputs(
        n_source, n_query, visibility_mask, mask_locally_constant
    )

    # All input validation and hard caps precede conversion and Q-by-N output
    # allocation. The caps bound maximum output storage to 500x500x3 values.
    x = x_input.astype(np.float64, copy=True)
    queries = query_input.astype(np.float64, copy=True)
    x.setflags(write=False)
    queries.setflags(write=False)
    weights_out = np.full((n_query, n_source), np.nan, dtype=np.float64)
    gradients_out = np.full((n_query, n_source, 3), np.nan, dtype=np.float64)
    support_counts = np.zeros(n_query, dtype=np.int64)
    condition_numbers = np.full(n_query, np.nan, dtype=np.float64)
    if mask_input is None:
        visible = np.ones(n_source, dtype=bool)
        mask_stable = np.ones(n_query, dtype=bool)
    else:
        visible = np.array(mask_input, dtype=bool, copy=True)
        visible.setflags(write=False)
        if stable_input is None:
            mask_stable = np.zeros(n_query, dtype=bool)
        elif isinstance(stable_input, bool):
            mask_stable = np.full(n_query, stable_input, dtype=bool)
        else:
            mask_stable = np.array(stable_input, dtype=bool, copy=True)

    statuses = [MLSStatus.UNKNOWN_POOR_SUPPORT] * n_query
    origin = np.asarray(config.reference_origin, dtype=np.float64)
    scale = float(config.reference_scale)
    radius = float(config.support_radius)
    source_basis = np.empty((n_source, 4), dtype=np.float64)
    source_basis[:, 0] = 1.0
    source_basis[:, 1:] = (x - origin) / scale
    if not np.all(np.isfinite(source_basis)):
        return _immutable_result(
            weights_out,
            gradients_out,
            [MLSStatus.UNKNOWN_ILL_CONDITIONED] * n_query,
            support_counts,
            condition_numbers,
        )

    for qi, query in enumerate(queries):
        if not mask_stable[qi]:
            statuses[qi] = MLSStatus.UNKNOWN_VISIBILITY_SWITCH
            continue

        delta = x - query
        distance = np.linalg.norm(delta, axis=1)
        q = distance / radius
        support_counts[qi] = int(np.count_nonzero(visible & (q < 1.0)))
        if np.any(visible & (np.abs(q - 1.0) <= config.switch_tolerance)):
            statuses[qi] = MLSStatus.UNKNOWN_SUPPORT_BOUNDARY
            continue
        active = np.flatnonzero(visible & (q < 1.0))
        if active.size < 4:
            statuses[qi] = MLSStatus.UNKNOWN_POOR_SUPPORT
            continue

        local_positions = x[active]
        duplicate_cutoff = config.duplicate_tolerance * radius
        duplicate_found = False
        for local_index in range(active.size - 1):
            later_delta = local_positions[local_index + 1 :] - local_positions[local_index]
            if np.any(np.einsum("ij,ij->i", later_delta, later_delta) <= duplicate_cutoff * duplicate_cutoff):
                duplicate_found = True
                break
        if duplicate_found:
            statuses[qi] = MLSStatus.UNKNOWN_DUPLICATE_SOURCES
            continue

        local_q = q[active]
        weight, dweight_dq = wendland_c2(local_q)
        local_distance = distance[active]
        local_delta = delta[active]
        dq_dx = np.zeros((active.size, 3), dtype=np.float64)
        nonzero_distance = local_distance > 0.0
        dq_dx[nonzero_distance] = (
            -local_delta[nonzero_distance]
            / (radius * local_distance[nonzero_distance, None])
        )
        dweight_dx = dweight_dq[:, None] * dq_dx

        p = source_basis[active]
        moment = np.einsum("n,ni,nj->ij", weight, p, p, optimize=True)
        dmoment = np.einsum("nk,ni,nj->kij", dweight_dx, p, p, optimize=True)
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

        query_basis = np.empty(4, dtype=np.float64)
        query_basis[0] = 1.0
        query_basis[1:] = (query - origin) / scale
        try:
            # Dual solve avoids constructing the N-by-N identity in the old
            # interpolation interface.
            dual_coefficients = np.linalg.solve(moment.T, query_basis)
            ddual_coefficients = np.empty((3, 4), dtype=np.float64)
            for axis in range(3):
                dp_dx = np.zeros(4, dtype=np.float64)
                dp_dx[axis + 1] = 1.0 / scale
                ddual_coefficients[axis] = np.linalg.solve(
                    moment.T, dp_dx - dmoment[axis].T @ dual_coefficients
                )
        except np.linalg.LinAlgError:
            statuses[qi] = MLSStatus.UNKNOWN_LINEAR_SOLVE_FAILED
            continue

        local_weights = weight * (p @ dual_coefficients)
        local_gradients = np.empty((active.size, 3), dtype=np.float64)
        p_dot_c = p @ dual_coefficients
        for axis in range(3):
            local_gradients[:, axis] = (
                dweight_dx[:, axis] * p_dot_c
                + weight * (p @ ddual_coefficients[axis])
            )
        if not np.all(np.isfinite(local_weights)) or not np.all(np.isfinite(local_gradients)):
            statuses[qi] = MLSStatus.UNKNOWN_LINEAR_SOLVE_FAILED
            continue

        weights_out[qi] = 0.0
        gradients_out[qi] = 0.0
        weights_out[qi, active] = local_weights
        gradients_out[qi, active, :] = local_gradients
        statuses[qi] = MLSStatus.OK

    return _immutable_result(weights_out, gradients_out, statuses, support_counts, condition_numbers)
