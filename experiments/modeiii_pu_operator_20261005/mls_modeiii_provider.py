"""Private continuous-MLS / normalized-front enrichment provider.

This source-only adapter maps the real 3-D MLS nodal shapes into the frozen
scalar antiplane coefficient registry. It defines no Heaviside enrichment,
visibility mask, solver, topology qualification, or crack-growth operation.
"""

from __future__ import annotations

import math

import numpy as np

from enriched_dof_contract import (
    AntiplaneShapeSample,
    BasisKind,
    DofKey,
    DofRegistry,
    FieldComponent,
    ProviderContractError,
    ProviderStatus,
    RegistryState,
    ScalarShapeContribution,
    validate_antiplane_sample,
    validate_problem_budget,
)
from mls3d_mls import MLS3DConfig, MLSStatus
from mls3d_shape import ShapeFunctions3DResult, shape_functions_3d
from field_reference import FieldResult, FieldStatus, ModeIIIFixedCase, evaluate_field


class ProviderOrderError(ValueError):
    """Raised when a caller requests a different query than the frozen batch."""


def _immutable_array(array: np.ndarray) -> np.ndarray:
    """Make a bounded, bytes-backed snapshot that cannot be made writable."""
    return np.frombuffer(array.tobytes(order="C"), dtype=array.dtype).reshape(array.shape)


def _checked_points(name: str, value: object, maximum_rows: int) -> np.ndarray:
    if type(value) is not np.ndarray:
        raise ValueError(f"{name} must be an exact numpy.ndarray")
    if value.ndim != 2 or value.shape[1] != 3:
        raise ValueError(f"{name} must have shape (n, 3)")
    rows = int(value.shape[0])
    if rows < 1 or rows > maximum_rows:
        raise ValueError(f"budgetRejected: {name} row count must be in [1, {maximum_rows}]")
    if value.dtype.kind not in "fiu" or value.dtype.itemsize > 8:
        raise ValueError(f"{name} must use real numeric values with at most 8-byte elements")
    if not np.all(np.isfinite(value)):
        raise ValueError(f"{name} must contain only finite coordinates")
    return _immutable_array(np.asarray(value, dtype=np.float64, order="C"))


class ContinuousMLSFrontProvider:
    """Provide BASE and unshifted-F1 rows over the frozen canonical registry.

    Let ``Fbar=F1/sqrt(R)``. The row for a source i is
    ``phi_i`` for BASE and ``phi_i*Fbar`` for FRONT_F1_UNSHIFTED. Gradients
    include the full product rule. Coefficients in both families therefore
    have displacement units (m). Source values of F1 are neither requested
    nor inferred, so sources on the front and crack cut are admissible for
    this unshifted basis.

    The only visibility mode is continuous, unmasked MLS. All MLS and exact
    field queries are pre-evaluated in the fixed batch. Any non-OK MLS or
    field row is retained as an UNKNOWN condition and yields no partial row.
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
        if type(registry) is not DofRegistry:
            raise ValueError("registry must be the frozen enriched_dof_contract DofRegistry")
        if type(source_positions) is not np.ndarray or source_positions.ndim != 2 or source_positions.shape[1] != 3:
            raise ValueError("source_positions must be a plain ndarray of shape (n,3)")
        if type(query_positions) is not np.ndarray or query_positions.ndim != 2 or query_positions.shape[1] != 3:
            raise ValueError("query_positions must be a plain ndarray of shape (q,3)")
        source_count = int(source_positions.shape[0])
        query_count = int(query_positions.shape[0])
        validate_problem_budget(source_count, registry.dof_count, query_count)
        if type(config) is not MLS3DConfig:
            raise ValueError("config must be the frozen MLS3DConfig")
        if type(field_case) is not ModeIIIFixedCase:
            raise ValueError("field_case must be the frozen ModeIIIFixedCase")
        if registry.dof_count != 2 * source_count:
            raise ProviderContractError("the fixed BASE plus FRONT basis requires exactly two DOFs per source")
        if type(source_ids) is not tuple or len(source_ids) != source_count:
            raise ValueError("source_ids must be a tuple with one ID per source position")
        if any(type(source_id) is not str or not source_id for source_id in source_ids):
            raise ValueError("source IDs must be nonempty strings")
        if len(set(source_ids)) != source_count:
            raise ValueError("source IDs must be unique")

        self.source_positions = _checked_points("source_positions", source_positions, 125)
        self.query_positions = _checked_points("query_positions", query_positions, 256)
        if (
            config.reference_origin != (0.0, 0.0, 0.0)
            or config.reference_scale != 0.01
            or config.support_radius != 0.025
        ):
            raise ValueError("MLS config must retain origin 0, reference scale R=0.01 m, and h=0.025 m support")
        if registry.state is not RegistryState.FROZEN:
            raise ProviderContractError("freeze the coefficient registry before constructing the provider")
        expected_keys = {
            DofKey(source_id, FieldComponent.UZ, BasisKind.BASE, 0)
            for source_id in source_ids
        } | {
            DofKey(source_id, FieldComponent.UZ, BasisKind.FRONT_F1_UNSHIFTED, 1)
            for source_id in source_ids
        }
        if set(registry.ordered_keys) != expected_keys:
            raise ProviderContractError(
                "registry must contain exactly BASE and FRONT_F1_UNSHIFTED UZ for every source"
            )

        self.source_ids = source_ids
        self.registry = registry
        self._next_query = 0
        self.mls_result: ShapeFunctions3DResult = shape_functions_3d(
            self.source_positions,
            self.query_positions,
            config,
        )
        # Both visibility keywords are intentionally omitted. A crack mask is
        # not part of this continuous-MLS PU construction.
        self.field_result: FieldResult = evaluate_field(self.query_positions, case=field_case)
        self.unknown_rows = self._collect_unknown_rows()
        self.values: np.ndarray | None = None
        self.gradients: np.ndarray | None = None
        self.coefficient_ids = registry.canonical_ids

        if self.unknown_rows:
            return
        if not np.all(np.isfinite(self.mls_result.weights)) or not np.all(
            np.isfinite(self.mls_result.gradients)
        ):
            raise ProviderContractError("MLS returned a non-finite OK row")

        source_row = {source_id: row for row, source_id in enumerate(source_ids)}
        values = np.empty((query_count, registry.dof_count), dtype=np.float64)
        gradients = np.empty((query_count, registry.dof_count, 3), dtype=np.float64)
        sqrt_radius = math.sqrt(field_case.radius_m)
        fbar = self.field_result.f1_m_sqrt / sqrt_radius
        grad_fbar = self.field_result.grad_f1_per_sqrt_m / sqrt_radius
        for slot, key in enumerate(registry.ordered_keys):
            source_row_id = source_row.get(key.source_id)
            if source_row_id is None:
                raise ProviderContractError("registry contains an ID without a source coordinate")
            if key.field_component is not FieldComponent.UZ:
                raise ProviderContractError("the restricted scalar operator only supports UZ keys")
            phi = self.mls_result.weights[:, source_row_id]
            grad_phi = self.mls_result.gradients[:, source_row_id, :]
            if key.basis_kind is BasisKind.BASE and key.basis_index == 0:
                values[:, slot] = phi
                gradients[:, slot, :] = grad_phi
            elif key.basis_kind is BasisKind.FRONT_F1_UNSHIFTED and key.basis_index == 1:
                values[:, slot] = phi * fbar
                gradients[:, slot, :] = (
                    fbar[:, None] * grad_phi
                    + phi[:, None] * grad_fbar
                )
            else:
                raise ProviderContractError(
                    f"unsupported basis key in this fixed provider: {key.basis_kind.value}"
                )
        if not np.all(np.isfinite(values)) or not np.all(np.isfinite(gradients)):
            raise ProviderContractError("enriched scalar rows contain non-finite values")
        self.values = _immutable_array(values)
        self.gradients = _immutable_array(gradients)

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
                    {
                        "query_index": query_id,
                        "source": "MODEIII_REFERENCE",
                        "status": status.value,
                    }
                )
        return rows

    @property
    def query_count(self) -> int:
        return int(self.query_positions.shape[0])

    @property
    def consumed_queries(self) -> int:
        return self._next_query

    def evaluate(self, material_position: tuple[float, float, float]) -> AntiplaneShapeSample:
        if self._next_query >= self.query_count:
            raise ProviderOrderError("provider was asked for more queries than the frozen batch")
        if type(material_position) is not tuple or len(material_position) != 3:
            raise ProviderOrderError("material_position must be an exact Vector3 tuple")
        if any(
            isinstance(value, bool)
            or not isinstance(value, (int, float, np.integer, np.floating))
            or not np.isfinite(float(value))
            for value in material_position
        ):
            raise ProviderOrderError("material_position must contain three finite real coordinates")
        query_vector = np.asarray(tuple(float(value) for value in material_position), dtype=np.float64)
        query_id = self._next_query
        if not np.array_equal(query_vector, self.query_positions[query_id]):
            raise ProviderOrderError(f"query order diverged from the frozen batch at row {query_id}")
        self._next_query += 1

        if self.unknown_rows:
            diagnostic = json_diagnostic_for_query(self.unknown_rows)
            return AntiplaneShapeSample((), (), ProviderStatus.UNKNOWN, diagnostic)
        if self.values is None or self.gradients is None:
            raise ProviderContractError("provider row storage is unavailable")
        contributions = tuple(
            ScalarShapeContribution(
                float(self.values[query_id, slot]),
                tuple(float(value) for value in self.gradients[query_id, slot, :]),
            )
            for slot in range(self.registry.dof_count)
        )
        sample = AntiplaneShapeSample(
            coefficient_ids=self.coefficient_ids,
            contributions=contributions,
            status=ProviderStatus.OK,
        )
        return validate_antiplane_sample(sample, self.registry)


def json_diagnostic_for_query(rows: list[dict[str, object]]) -> str:
    """Return all exact batch failure records without numerical repair."""
    import json

    return json.dumps(
        {"batch_preflight_blocked": True, "unknown_rows": rows},
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )
