"""Continuous 3-D MLS rows for all 162 semantic enriched coefficients.

The underlying elastic assembler consumes one scalar shape row per vector
basis site. Here each source's BASE and normalized FRONT_F1 scalar shape is
computed once, then the full-vector assembler maps that same scalar row to UX,
UY, and UZ through the semantic registry. The crack is represented only by
the global unshifted Mode III enrichment in this fixed bridge; MLS remains
continuous and unmasked.
"""

from __future__ import annotations

import json
import math

import numpy as np

from enriched_dof_contract import (
    BasisKind,
    DofKey,
    DofRegistry,
    FieldComponent,
    RegistryState,
    validate_problem_budget,
)
from efgm_elastic_operator import ShapeSample, ShapeStatus
from field_reference import FieldStatus, ModeIIIFixedCase, evaluate_field
from mls3d_mls import MLS3DConfig, MLSStatus
from mls3d_shape import ShapeFunctions3DResult, shape_functions_3d

from full3d_cracked_elastic_operator import FULL3D_BASIS_KINDS, basis_row_ids


MAX_SOURCE_COUNT = 27
MAX_QUERY_COUNT = 256


class ProviderOrderError(ValueError):
    """Raised when a caller requests a coordinate outside the frozen batch order."""


def _immutable_array(array: np.ndarray) -> np.ndarray:
    return np.frombuffer(array.tobytes(order="C"), dtype=array.dtype).reshape(array.shape)


def _checked_points(name: str, value: object, maximum_rows: int) -> np.ndarray:
    if type(value) is not np.ndarray:
        raise ValueError(f"{name} must be an exact numpy.ndarray")
    if value.ndim != 2 or value.shape[1] != 3:
        raise ValueError(f"{name} must have shape (n, 3)")
    rows = int(value.shape[0])
    if not 1 <= rows <= maximum_rows:
        raise ValueError(f"budgetRejected: {name} row count must be in [1, {maximum_rows}]")
    if value.dtype.kind not in "fiu" or value.dtype.itemsize > 8:
        raise ValueError(f"{name} must use real numeric values with at most 8-byte elements")
    if not np.all(np.isfinite(value)):
        raise ValueError(f"{name} must contain only finite coordinates")
    if value.dtype.kind in "iu" and (
        np.any(value > (1 << 53)) or np.any(value < -(1 << 53))
    ):
        raise ValueError(f"{name} integer values must be exactly representable in float64")
    return _immutable_array(np.asarray(value, dtype=np.float64, order="C"))


class ContinuousMLSFull3DFrontProvider:
    """Provide cached BASE/F1 scalar rows for UX, UY, and UZ assembly.

    For ``Fbar=F1/sqrt(R)``, the scalar basis rows are ``phi_i`` and
    ``phi_i*Fbar``. Their gradients include the full product rule. All 27
    sources receive the global front term and there is no Heaviside term or
    source-side value of F1. Source coefficients remain physical
    displacements in metres.
    """

    def __init__(
        self,
        source_positions: np.ndarray,
        source_ids: tuple[str, ...],
        query_positions: np.ndarray,
        registry: DofRegistry,
        config: MLS3DConfig,
        field_case: ModeIIIFixedCase,
    ) -> None:
        if type(registry) is not DofRegistry or registry.state is not RegistryState.FROZEN:
            raise ValueError("registry must be the frozen semantic full-vector registry")
        if type(source_positions) is not np.ndarray or source_positions.ndim != 2 or source_positions.shape[1] != 3:
            raise ValueError("source_positions must be a plain ndarray of shape (n,3)")
        if type(query_positions) is not np.ndarray or query_positions.ndim != 2 or query_positions.shape[1] != 3:
            raise ValueError("query_positions must be a plain ndarray of shape (q,3)")
        source_count = int(source_positions.shape[0])
        query_count = int(query_positions.shape[0])
        if source_count != MAX_SOURCE_COUNT:
            raise ValueError("the fixed full-3D bridge requires exactly 27 source sites")
        validate_problem_budget(source_count, registry.dof_count, query_count)
        if query_count > MAX_QUERY_COUNT:
            raise ValueError(f"budgetRejected: query count exceeds {MAX_QUERY_COUNT}")
        if type(config) is not MLS3DConfig:
            raise ValueError("config must be the frozen MLS3DConfig")
        if type(field_case) is not ModeIIIFixedCase:
            raise ValueError("field_case must be the frozen ModeIIIFixedCase")
        if type(source_ids) is not tuple or len(source_ids) != source_count:
            raise ValueError("source_ids must be a tuple with one ID per source position")
        if any(type(source_id) is not str or not source_id for source_id in source_ids):
            raise ValueError("source IDs must be nonempty strings")
        if len(set(source_ids)) != source_count:
            raise ValueError("source IDs must be unique")
        if (
            config.reference_origin != (0.0, 0.0, 0.0)
            or float(config.reference_scale) != 0.01
            or float(config.support_radius) != 0.025
        ):
            raise ValueError("MLS config must retain origin 0, reference scale R=0.01 m, and h=0.025 m")
        expected_keys = {
            DofKey(source_id, component, basis, 0 if basis is BasisKind.BASE else 1)
            for source_id in source_ids
            for component in (FieldComponent.UX, FieldComponent.UY, FieldComponent.UZ)
            for basis in FULL3D_BASIS_KINDS
        }
        if registry.source_count != source_count or set(registry.ordered_keys) != expected_keys:
            raise ValueError("registry must contain UX/UY/UZ BASE and FRONT_F1 for all 27 sources")

        self.source_positions = _checked_points("source_positions", source_positions, MAX_SOURCE_COUNT)
        self.query_positions = _checked_points("query_positions", query_positions, MAX_QUERY_COUNT)
        self.source_ids = source_ids
        self.registry = registry
        self.config = config
        self.field_case = field_case
        self.basis_row_ids = basis_row_ids(source_ids)
        self.coefficient_ids = registry.canonical_ids
        self._next_query = 0

        # This instance evaluates one fixed batch of at most 256 points with
        # continuous visibility. The prepared case uses 256 volume points and
        # a separate 128-point boundary batch. No crack mask is applied.
        self.mls_result: ShapeFunctions3DResult = shape_functions_3d(
            self.source_positions,
            self.query_positions,
            config,
        )
        self.field_result = evaluate_field(self.query_positions, case=field_case)
        self.unknown_rows = self._collect_unknown_rows()
        self.values: np.ndarray | None = None
        self.gradients: np.ndarray | None = None
        if self.unknown_rows:
            return
        if not np.all(np.isfinite(self.mls_result.weights)) or not np.all(
            np.isfinite(self.mls_result.gradients)
        ):
            raise FloatingPointError("MLS returned a non-finite successful row")

        query_count = int(self.query_positions.shape[0])
        basis_count = len(self.basis_row_ids)
        values = np.empty((query_count, basis_count), dtype=np.float64)
        gradients = np.empty((query_count, basis_count, 3), dtype=np.float64)
        source_row = {source_id: row for row, source_id in enumerate(source_ids)}
        sqrt_radius = math.sqrt(field_case.radius_m)
        fbar = self.field_result.f1_m_sqrt / sqrt_radius
        grad_fbar = self.field_result.grad_f1_per_sqrt_m / sqrt_radius
        for basis_row, (source_id, basis) in enumerate(self.basis_row_ids):
            source_index = source_row[source_id]
            phi = self.mls_result.weights[:, source_index]
            grad_phi = self.mls_result.gradients[:, source_index, :]
            if basis is BasisKind.BASE:
                values[:, basis_row] = phi
                gradients[:, basis_row, :] = grad_phi
            elif basis is BasisKind.FRONT_F1_UNSHIFTED:
                values[:, basis_row] = phi * fbar
                gradients[:, basis_row, :] = (
                    fbar[:, None] * grad_phi + phi[:, None] * grad_fbar
                )
            else:
                raise ValueError(f"unsupported basis kind {basis.value}")
        if not np.all(np.isfinite(values)) or not np.all(np.isfinite(gradients)):
            raise FloatingPointError("enriched full-vector shape rows contain non-finite values")
        self.values = _immutable_array(values)
        self.gradients = _immutable_array(gradients)
        self._support_ids = _immutable_array(np.arange(basis_count, dtype=np.int64))

    def _collect_unknown_rows(self) -> list[dict[str, object]]:
        rows: list[dict[str, object]] = []
        for query_id, status in enumerate(self.mls_result.statuses):
            if status is not MLSStatus.OK:
                rows.append(
                    {
                        "query_index": query_id,
                        "source": "MLS",
                        "status": status.value,
                        "support_count": int(self.mls_result.support_counts[query_id]),
                        "condition_number": (
                            float(self.mls_result.condition_numbers[query_id])
                            if np.isfinite(self.mls_result.condition_numbers[query_id])
                            else None
                        ),
                    }
                )
        for query_id in range(len(self.field_result.statuses)):
            status = self.field_result.status_at(query_id)
            if status is not FieldStatus.OK:
                rows.append(
                    {"query_index": query_id, "source": "MODEIII_REFERENCE", "status": status.value}
                )
        return rows

    @property
    def query_count(self) -> int:
        return int(self.query_positions.shape[0])

    @property
    def consumed_queries(self) -> int:
        return self._next_query

    def evaluate(self, material_position: np.ndarray) -> ShapeSample:
        if self._next_query >= self.query_count:
            raise ProviderOrderError("provider was asked for more queries than its frozen batch")
        if type(material_position) is not np.ndarray or material_position.shape != (3,):
            raise ProviderOrderError("material_position must be a plain Vector3 ndarray")
        if material_position.dtype.kind not in "fiu" or not np.all(np.isfinite(material_position)):
            raise ProviderOrderError("material_position must contain finite real coordinates")
        query_id = self._next_query
        if not np.array_equal(np.asarray(material_position, dtype=np.float64), self.query_positions[query_id]):
            raise ProviderOrderError(f"query order diverged from the frozen batch at row {query_id}")
        self._next_query += 1

        if self.unknown_rows:
            mls_status = next(
                (
                    row["status"]
                    for row in self.unknown_rows
                    if int(row["query_index"]) == query_id and row["source"] == "MLS"
                ),
                None,
            )
            status = (
                ShapeStatus.UNKNOWN_POOR_SUPPORT
                if mls_status == MLSStatus.UNKNOWN_POOR_SUPPORT.value
                else ShapeStatus.UNKNOWN
            )
            diagnostic = json.dumps(
                {"batch_preflight_blocked": True, "unknown_rows": self.unknown_rows},
                sort_keys=True,
                separators=(",", ":"),
                allow_nan=False,
            )
            empty_i = np.empty(0, dtype=np.int64)
            return ShapeSample(
                weights=np.empty(0, dtype=np.float64),
                gradients=np.empty((0, 3), dtype=np.float64),
                support_ids=empty_i,
                status=status,
                diagnostic=diagnostic,
            )
        if self.values is None or self.gradients is None:
            raise RuntimeError("successful provider has no cached shape values")
        return ShapeSample(
            weights=self.values[query_id],
            gradients=self.gradients[query_id],
            support_ids=self._support_ids,
            status=ShapeStatus.OK,
        )
