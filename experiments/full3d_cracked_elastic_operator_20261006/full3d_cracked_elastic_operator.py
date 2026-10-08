"""Complete 3-D elastic stiffness for the fixed enriched straight-crack case.

This module reuses the reviewed full-vector small-strain EFGM operator. The
54 scalar basis functions (BASE and normalized FRONT_F1 for each of 27 source
sites) are exposed to that operator as 54 vector basis sites. The resulting
node-major matrix is permuted into the frozen semantic registry order, which
contains UX/UY/UZ for each source and basis family. No solver, load
application, visibility mask, regularization, or crack growth is implemented.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np

from enriched_dof_contract import (
    BasisKind,
    DofKey,
    DofRegistry,
    FieldComponent,
    RegistryState,
    validate_problem_budget,
)
from efgm_elastic_operator import (
    AssemblyStatus,
    IsotropicMaterial,
    ShapeStatus,
    assemble_elastic_operator,
)


FULL3D_BASIS_KINDS = (BasisKind.BASE, BasisKind.FRONT_F1_UNSHIFTED)
FULL3D_COMPONENTS = (FieldComponent.UX, FieldComponent.UY, FieldComponent.UZ)
COMPONENT_INDEX = {
    FieldComponent.UX: 0,
    FieldComponent.UY: 1,
    FieldComponent.UZ: 2,
}
MAX_SOURCE_COUNT = 125
MAX_QUADRATURE_POINTS = 256


class Full3DShapeProvider(Protocol):
    """Provider of 54 scalar shape functions lifted to vector elasticity."""

    basis_row_ids: tuple[tuple[str, BasisKind], ...]

    def evaluate(self, material_position: np.ndarray) -> object:
        """Return one efgm_elastic_operator.ShapeSample row."""


@dataclass(frozen=True)
class Full3DAssemblyResult:
    """Full semantic-order stiffness or a fail-closed support status."""

    status: AssemblyStatus
    stiffness: np.ndarray | None
    dof_count: int
    quadrature_count: int
    blocked_quadrature_index: int | None = None
    shape_status: ShapeStatus | None = None
    diagnostic: str | None = None


def _immutable_array(array: np.ndarray) -> np.ndarray:
    return np.frombuffer(array.tobytes(order="C"), dtype=array.dtype).reshape(array.shape)


def _check_source_inputs(
    source_positions: object,
    source_ids: object,
) -> tuple[np.ndarray, tuple[str, ...]]:
    if type(source_positions) is not np.ndarray:
        raise ValueError("source_positions must be an exact numpy.ndarray")
    if source_positions.ndim != 2 or source_positions.shape[1] != 3:
        raise ValueError("source_positions must have shape (n, 3)")
    source_count = int(source_positions.shape[0])
    if not 1 <= source_count <= MAX_SOURCE_COUNT:
        raise ValueError(f"budgetRejected: source count must be in [1, {MAX_SOURCE_COUNT}]")
    if (
        source_positions.dtype.kind not in "fiu"
        or source_positions.dtype.itemsize > 8
        or not np.all(np.isfinite(source_positions))
    ):
        raise ValueError("source_positions must contain finite real coordinates")
    if source_positions.dtype.kind in "iu" and (
        np.any(source_positions > (1 << 53)) or np.any(source_positions < -(1 << 53))
    ):
        raise ValueError("integer source coordinates must be exactly representable in float64")
    if type(source_ids) is not tuple or len(source_ids) != source_count:
        raise ValueError("source_ids must be a tuple with one ID per source position")
    if any(type(source_id) is not str or not source_id for source_id in source_ids):
        raise ValueError("source IDs must be nonempty strings")
    if len(set(source_ids)) != source_count:
        raise ValueError("source IDs must be unique")
    return np.asarray(source_positions, dtype=np.float64, order="C"), source_ids


def basis_row_ids(source_ids: tuple[str, ...]) -> tuple[tuple[str, BasisKind], ...]:
    """Return deterministic scalar-basis rows, source-major then basis-major."""
    if type(source_ids) is not tuple or not source_ids:
        raise ValueError("source_ids must be a nonempty tuple")
    return tuple((source_id, basis) for source_id in source_ids for basis in FULL3D_BASIS_KINDS)


def build_full3d_registry(source_ids: tuple[str, ...]) -> DofRegistry:
    """Freeze exactly three physical components times two basis terms per source."""
    if type(source_ids) is not tuple or not source_ids:
        raise ValueError("source_ids must be a nonempty tuple")
    registry = DofRegistry()
    for source_id in source_ids:
        for component in FULL3D_COMPONENTS:
            registry.register(DofKey(source_id, component, BasisKind.BASE, 0))
            registry.register(
                DofKey(source_id, component, BasisKind.FRONT_F1_UNSHIFTED, 1)
            )
    registry.freeze()
    return registry


def basis_site_positions(
    source_positions: np.ndarray,
    source_ids: tuple[str, ...],
) -> np.ndarray:
    """Duplicate each source coordinate for its BASE and FRONT_F1 basis rows."""
    positions, ids = _check_source_inputs(source_positions, source_ids)
    source_row = {source_id: row for row, source_id in enumerate(ids)}
    rows = basis_row_ids(ids)
    result = np.ascontiguousarray(
        np.asarray([positions[source_row[source_id]] for source_id, _ in rows], dtype=np.float64)
    )
    return result


def _natural_dof_for_canonical_slot(
    source_ids: tuple[str, ...],
    registry: DofRegistry,
) -> np.ndarray:
    if type(registry) is not DofRegistry or registry.state is not RegistryState.FROZEN:
        raise ValueError("registry must be a frozen DofRegistry")
    rows = basis_row_ids(source_ids)
    expected = {
        DofKey(source_id, component, basis, 0 if basis is BasisKind.BASE else 1)
        for source_id in source_ids
        for component in FULL3D_COMPONENTS
        for basis in FULL3D_BASIS_KINDS
    }
    if registry.source_count != len(source_ids) or set(registry.ordered_keys) != expected:
        raise ValueError("registry must contain exactly UX/UY/UZ BASE and FRONT_F1 for every source")
    validate_problem_budget(len(source_ids), registry.dof_count, 0)
    if registry.dof_count != 6 * len(source_ids):
        raise ValueError("full 3-D registry must contain exactly six coefficients per source")
    basis_row = {row_id: row for row, row_id in enumerate(rows)}
    natural = np.empty(registry.dof_count, dtype=np.int64)
    for key in registry.ordered_keys:
        row = basis_row[(key.source_id, key.basis_kind)]
        natural[registry.slot_of(key)] = 3 * row + COMPONENT_INDEX[key.field_component]
    if not np.array_equal(np.sort(natural), np.arange(registry.dof_count, dtype=np.int64)):
        raise ValueError("semantic registry-to-vector-basis mapping is not a permutation")
    return natural


def assemble_full3d_enriched_operator(
    source_positions: np.ndarray,
    source_ids: tuple[str, ...],
    quadrature_positions: np.ndarray,
    quadrature_volume_weights: np.ndarray,
    material: IsotropicMaterial,
    shape_provider: Full3DShapeProvider,
    registry: DofRegistry,
) -> Full3DAssemblyResult:
    """Assemble and semantically order ``K = integral B.T D B dV``.

    Each scalar basis function supplies one vector basis site. The reused
    operator constructs the complete 3-D engineering-strain matrix with order
    ``[eps_xx, eps_yy, eps_zz, gamma_xy, gamma_yz, gamma_xz]`` and the full
    isotropic constitutive matrix, including lambda coupling and all three
    shear terms. A failed provider row returns no partial stiffness.
    """
    positions, ids = _check_source_inputs(source_positions, source_ids)
    if type(registry) is not DofRegistry or registry.state is not RegistryState.FROZEN:
        raise ValueError("registry must be frozen before assembly")
    if type(material) is not IsotropicMaterial:
        raise ValueError("material must be an IsotropicMaterial")
    expected_rows = basis_row_ids(ids)
    if getattr(shape_provider, "basis_row_ids", None) != expected_rows:
        raise ValueError("shape provider basis rows differ from the frozen source/basis order")
    natural = _natural_dof_for_canonical_slot(ids, registry)
    dof_count = registry.dof_count
    if type(quadrature_positions) is not np.ndarray or quadrature_positions.ndim != 2:
        raise ValueError("quadrature_positions must be a plain rank-2 numpy.ndarray")
    quadrature_count = int(quadrature_positions.shape[0])
    if quadrature_count > MAX_QUADRATURE_POINTS:
        raise ValueError(f"budgetRejected: quadrature count exceeds {MAX_QUADRATURE_POINTS}")

    sites = basis_site_positions(positions, ids)
    raw = assemble_elastic_operator(
        coefficient_positions=sites,
        quadrature_positions=quadrature_positions,
        quadrature_volume_weights=quadrature_volume_weights,
        material=material,
        shape_provider=shape_provider,
    )
    if raw.status is not AssemblyStatus.OK or raw.stiffness is None:
        return Full3DAssemblyResult(
            status=raw.status,
            stiffness=None,
            dof_count=dof_count,
            quadrature_count=raw.quadrature_count,
            blocked_quadrature_index=raw.blocked_quadrature_index,
            shape_status=raw.shape_status,
            diagnostic=raw.diagnostic,
        )
    stiffness = raw.stiffness[np.ix_(natural, natural)]
    if stiffness.shape != (dof_count, dof_count) or not np.all(np.isfinite(stiffness)):
        raise FloatingPointError("semantic-order stiffness is malformed or non-finite")
    return Full3DAssemblyResult(
        status=AssemblyStatus.OK,
        stiffness=_immutable_array(stiffness),
        dof_count=dof_count,
        quadrature_count=raw.quadrature_count,
    )


def rigid_mode_coefficient_vectors(
    source_positions: np.ndarray,
    source_ids: tuple[str, ...],
    registry: DofRegistry,
) -> np.ndarray:
    """Return 3 translations and 3 rotations with all enrichment entries zero.

    For BASE coefficients, translations are constant component values and
    rotations are ``omega x (x - source_centroid)`` at each source coordinate.
    This uses linear completeness of the frozen MLS shape basis. FRONT_F1
    coefficients remain exactly zero in all six diagnostic mode vectors.
    """
    positions, ids = _check_source_inputs(source_positions, source_ids)
    _natural_dof_for_canonical_slot(ids, registry)
    centered = positions - np.mean(positions, axis=0)
    source_row = {source_id: row for row, source_id in enumerate(ids)}
    modes = np.zeros((registry.dof_count, 6), dtype=np.float64)
    axes = np.eye(3, dtype=np.float64)
    for slot, key in enumerate(registry.ordered_keys):
        if key.basis_kind is not BasisKind.BASE:
            continue
        component = COMPONENT_INDEX[key.field_component]
        modes[slot, component] = 1.0
        position = centered[source_row[key.source_id]]
        for rotation_axis in range(3):
            modes[slot, 3 + rotation_axis] = np.cross(axes[rotation_axis], position)[component]
    return _immutable_array(modes)
