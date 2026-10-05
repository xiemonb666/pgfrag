"""Private adapter from frozen phase-2 MLS shapes to the EFGM operator.

This adapter uses the real ``shape_functions_3d`` API. It binds MLS source row
order to operator coefficient ID order and batches one immutable result for
the exact volume-quadrature query sequence. The only visibility mode here is
identity/no-crack: no mask is passed and ``mask_locally_constant=True`` is
never asserted. UNKNOWN MLS rows remain UNKNOWN and cause operator assembly to
block without a partial stiffness matrix.
"""

from __future__ import annotations

import numpy as np

from efgm_elastic_operator import (
    MAX_NODES,
    MAX_QUADRATURE_POINTS,
    ShapeSample,
    ShapeStatus,
)
from mls3d_mls import MLS3DConfig, MLSStatus
from mls3d_shape import ShapeFunctions3DResult, shape_functions_3d


class MLSAdapterContractError(ValueError):
    """Raised when the cached MLS batch is used with a different point order."""


def _immutable_array(array: np.ndarray) -> np.ndarray:
    """Copy bounded data onto immutable bytes-backed storage."""
    return np.frombuffer(array.tobytes(order="C"), dtype=array.dtype).reshape(array.shape)


def _checked_points(name: str, value: object, maximum_rows: int) -> np.ndarray:
    if type(value) is not np.ndarray:
        raise ValueError(f"{name} must be a plain numpy.ndarray")
    if value.ndim != 2 or value.shape[1] != 3:
        raise ValueError(f"{name} must have shape (n, 3)")
    row_count = int(value.shape[0])
    if row_count == 0:
        raise ValueError(f"{name} must contain at least one row")
    if row_count > maximum_rows:
        raise ValueError(
            f"budgetRejected: {name} has {row_count} rows; maximum is {maximum_rows}"
        )
    if value.dtype.kind not in "fiu" or not np.all(np.isfinite(value)):
        raise ValueError(f"{name} must contain finite real numeric coordinates")
    copied = np.array(value, dtype=np.float64, copy=True)
    return _immutable_array(copied)


class MLSShapeProvider:
    """Adapt one batched MLS evaluation to point-wise operator samples.

    The exact source input order is the global coefficient order, so successful
    rows are exposed as dense `0..N-1` values and gradients. This class accepts
    no visibility mask argument. A later explicit-crack adapter needs a
    surface-backed visibility computation and query-local stability evidence;
    it must not use this identity-visibility class for cracked supports.
    """

    def __init__(
        self,
        source_positions: np.ndarray,
        query_positions: np.ndarray,
        config: MLS3DConfig,
    ) -> None:
        if type(config) is not MLS3DConfig:
            raise ValueError("config must be the frozen MLS3DConfig")
        self.source_positions = _checked_points(
            "source_positions", source_positions, MAX_NODES
        )
        self.query_positions = _checked_points(
            "query_positions", query_positions, MAX_QUADRATURE_POINTS
        )
        support_ids = np.arange(self.source_positions.shape[0], dtype=np.int64)
        self.support_ids = _immutable_array(support_ids)
        # Deliberately omit both visibility_mask and mask_locally_constant.
        self.result: ShapeFunctions3DResult = shape_functions_3d(
            self.source_positions,
            self.query_positions,
            config,
        )
        self._next_query = 0

    @property
    def query_count(self) -> int:
        return int(self.query_positions.shape[0])

    @property
    def consumed_queries(self) -> int:
        return self._next_query

    def evaluate(self, material_position: np.ndarray) -> ShapeSample:
        if self._next_query >= self.query_count:
            raise MLSAdapterContractError("operator requested more queries than the MLS batch")
        if type(material_position) is not np.ndarray:
            raise MLSAdapterContractError("operator query must be a plain numpy.ndarray")
        if material_position.shape != (3,) or material_position.dtype.kind not in "fiu":
            raise MLSAdapterContractError("operator query must be a real coordinate vector of shape (3,)")
        expected = self.query_positions[self._next_query]
        if not np.array_equal(material_position, expected):
            raise MLSAdapterContractError(
                f"operator query order diverged from MLS batch at index {self._next_query}"
            )

        query_id = self._next_query
        self._next_query += 1
        status = self.result.statuses[query_id]
        if status is MLSStatus.OK:
            return ShapeSample(
                weights=self.result.weights[query_id],
                gradients=self.result.gradients[query_id],
                support_ids=self.support_ids,
                status=ShapeStatus.OK,
            )

        reason = self.result.reasons[query_id] or status.value
        support_count = int(self.result.support_counts[query_id])
        condition_number = float(self.result.condition_numbers[query_id])
        if not np.isfinite(condition_number):
            condition_text = "unknown"
        else:
            condition_text = f"{condition_number:.17g}"
        diagnostic = (
            f"MLSStatus.{reason}; support_count={support_count}; "
            f"condition_number={condition_text}"
        )
        shape_status = (
            ShapeStatus.UNKNOWN_POOR_SUPPORT
            if status is MLSStatus.UNKNOWN_POOR_SUPPORT
            else ShapeStatus.UNKNOWN
        )
        # The result's UNKNOWN arrays contain NaNs; the operator sees only the
        # explicit status and blocks before validating or copying those arrays.
        return ShapeSample(
            weights=self.result.weights[query_id],
            gradients=self.result.gradients[query_id],
            support_ids=self.support_ids,
            status=shape_status,
            diagnostic=diagnostic,
        )
