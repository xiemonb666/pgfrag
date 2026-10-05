"""SOURCE-ONLY registry and provider protocol for PU-enriched coefficients.

This module is a bounded contract primitive. It does not build an MLS basis,
qualify crack/front supports, assemble a stiffness matrix, or run a solver.
The antiplane provider protocol is deliberately limited to a scalar ``u_z``
field and returns one physical scalar shape value and its analytic spatial
gradient for every registered coefficient.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import json
import math
from numbers import Real
from typing import Protocol, TypeAlias


MAX_TOTAL_DOF = 375
MAX_SOURCE_NODES = 125
MAX_QUADRATURE_POINTS = 256

Vector3: TypeAlias = tuple[float, float, float]
EngineeringStrain6: TypeAlias = tuple[float, float, float, float, float, float]


class FieldComponent(str, Enum):
    """Physical displacement component attached to a coefficient."""

    UX = "ux"
    UY = "uy"
    UZ = "uz"


class BasisKind(str, Enum):
    """Basis family. Crack-surface and shifted/unshifted front terms have distinct IDs."""

    BASE = "base"
    HEAVISIDE = "heaviside"
    FRONT_F1_SHIFTED = "front_f1_shifted"
    FRONT_F1_UNSHIFTED = "front_f1_unshifted"


@dataclass(frozen=True)
class DofKey:
    """Stable semantic identity for one scalar coefficient.

    ``basis_index`` is zero for BASE and one for the current non-BASE kinds
    (the first crack-surface or front function). ``source_id`` is a stable
    caller-owned source identity; it is not a row number and is never
    synthesized from a duplicate physical node.
    """

    source_id: str
    field_component: FieldComponent
    basis_kind: BasisKind
    basis_index: int

    def __post_init__(self) -> None:
        if type(self.source_id) is not str or not self.source_id or self.source_id.strip() != self.source_id:
            raise ValueError("source_id must be a nonempty canonical string without edge whitespace")
        if type(self.field_component) is not FieldComponent:
            raise ValueError("field_component must be a FieldComponent")
        if type(self.basis_kind) is not BasisKind:
            raise ValueError("basis_kind must be a BasisKind")
        if type(self.basis_index) is not int:
            raise ValueError("basis_index must be an integer (bool is not accepted)")
        if self.basis_kind is BasisKind.BASE:
            if self.basis_index != 0:
                raise ValueError("BASE coefficients require basis_index=0")
        elif self.basis_index != 1:
            raise ValueError("this prototype requires basis_index=1 for non-BASE basis kinds")

    @property
    def canonical_id(self) -> str:
        """Return a deterministic, collision-safe ID independent of insertion order."""
        return json.dumps(
            [self.source_id, self.field_component.value, self.basis_kind.value, self.basis_index],
            ensure_ascii=True,
            separators=(",", ":"),
        )

    @property
    def sort_key(self) -> tuple[str, str, str, int]:
        return (
            self.source_id,
            self.field_component.value,
            self.basis_kind.value,
            self.basis_index,
        )


class RegistryState(str, Enum):
    OPEN = "OPEN"
    FROZEN = "FROZEN"


class DofRegistry:
    """Unique source/component/basis registry with deterministic frozen slots.

    A source that owns BASE and enrichment coefficients remains one source
    node while contributing multiple DOFs. Re-adding the same key is
    idempotent. Dense integer slots exist only after ``freeze`` and are sorted
    by the full semantic key, not by registration order.
    """

    __slots__ = ("_keys", "_source_ids", "_state", "_ordered_keys", "_slot_by_id")

    def __init__(self) -> None:
        self._keys: set[DofKey] = set()
        self._source_ids: set[str] = set()
        self._state = RegistryState.OPEN
        self._ordered_keys: tuple[DofKey, ...] | None = None
        self._slot_by_id: dict[str, int] | None = None

    @property
    def state(self) -> RegistryState:
        return self._state

    @property
    def dof_count(self) -> int:
        return len(self._keys)

    @property
    def source_count(self) -> int:
        return len(self._source_ids)

    def register(self, key: DofKey) -> str:
        """Register one key or fail before mutating the registry on a cap error."""
        if self._state is not RegistryState.OPEN:
            raise RuntimeError("registry is frozen")
        if type(key) is not DofKey:
            raise TypeError("registry accepts DofKey instances only")
        if key in self._keys:
            return key.canonical_id

        # Check both hard caps before inserting either the key or its source.
        next_dof_count = len(self._keys) + 1
        adds_source = key.source_id not in self._source_ids
        next_source_count = len(self._source_ids) + int(adds_source)
        validate_problem_budget(next_source_count, next_dof_count, 0)

        self._keys.add(key)
        if adds_source:
            self._source_ids.add(key.source_id)
        return key.canonical_id

    def freeze(self) -> tuple[DofKey, ...]:
        """Freeze registration and return keys in the stable canonical order."""
        if self._state is RegistryState.OPEN:
            if not self._keys:
                raise ValueError("cannot freeze an empty DOF registry")
            ordered = tuple(sorted(self._keys, key=lambda key: key.sort_key))
            self._ordered_keys = ordered
            self._slot_by_id = {key.canonical_id: slot for slot, key in enumerate(ordered)}
            self._state = RegistryState.FROZEN
        return self.ordered_keys

    @property
    def ordered_keys(self) -> tuple[DofKey, ...]:
        if self._state is not RegistryState.FROZEN or self._ordered_keys is None:
            raise RuntimeError("freeze the registry before requesting ordered keys")
        return self._ordered_keys

    @property
    def canonical_ids(self) -> tuple[str, ...]:
        if self._state is not RegistryState.FROZEN:
            raise RuntimeError("freeze the registry before requesting canonical_ids")
        return tuple(key.canonical_id for key in self._ordered_keys)

    def slot_of(self, key_or_id: DofKey | str) -> int:
        if self._state is not RegistryState.FROZEN or self._slot_by_id is None:
            raise RuntimeError("freeze the registry before requesting coefficient slots")
        canonical_id = key_or_id.canonical_id if type(key_or_id) is DofKey else key_or_id
        if type(canonical_id) is not str or canonical_id not in self._slot_by_id:
            raise KeyError("coefficient is not registered")
        return self._slot_by_id[canonical_id]


class ProviderStatus(str, Enum):
    OK = "OK"
    UNKNOWN = "UNKNOWN"
    BLOCKED = "BLOCKED"


@dataclass(frozen=True)
class ScalarShapeContribution:
    """One coefficient's physical scalar shape and analytic gradient."""

    value: float
    gradient: Vector3


@dataclass(frozen=True)
class AntiplaneShapeSample:
    """A full dense row over registered ``u_z`` coefficients in canonical order."""

    coefficient_ids: tuple[str, ...]
    contributions: tuple[ScalarShapeContribution, ...]
    status: ProviderStatus
    diagnostic: str | None = None


class AntiplaneShapeProvider(Protocol):
    """Provider for pure antiplane ``u_z`` field values and analytic gradients.

    This scalar contract is limited to the straight-front Mode III gate. It is
    not a general 3-D elasticity provider and does not assert support topology.
    """

    def evaluate(self, material_position: Vector3) -> AntiplaneShapeSample:
        """Return the physical ``u_z`` shape row at one material point."""


class ProviderContractError(ValueError):
    """Raised when an OK provider row does not match this registry contract."""


def _is_finite_real(value: object) -> bool:
    if isinstance(value, bool) or not isinstance(value, Real):
        return False
    try:
        return math.isfinite(float(value))
    except (OverflowError, ValueError):
        return False


def validate_problem_budget(source_count: int, dof_count: int, quadrature_count: int) -> None:
    """Validate hard limits using scalar counts before caller copies/allocations."""
    counts = (source_count, dof_count, quadrature_count)
    if any(type(count) is not int or count < 0 for count in counts):
        raise ValueError("source, DOF, and quadrature counts must be nonnegative integers")
    if source_count > MAX_SOURCE_NODES:
        raise ValueError(f"budgetRejected: source_count exceeds {MAX_SOURCE_NODES}")
    if dof_count > MAX_TOTAL_DOF:
        raise ValueError(f"budgetRejected: dof_count exceeds {MAX_TOTAL_DOF}")
    if quadrature_count > MAX_QUADRATURE_POINTS:
        raise ValueError(f"budgetRejected: quadrature_count exceeds {MAX_QUADRATURE_POINTS}")


def validate_antiplane_sample(sample: object, registry: DofRegistry) -> AntiplaneShapeSample:
    """Validate a full row and its registry alignment without allocating arrays.

    Any UNKNOWN/BLOCKED status is returned as a fail-closed row and must stop
    assembly. It may omit all numerical entries; partial non-OK rows are never
    consumed. This validator checks data shape and identity only, not topology,
    enrichment admissibility, approximation rank, or convergence.
    """
    if type(registry) is not DofRegistry or registry.state is not RegistryState.FROZEN:
        raise ProviderContractError("provider validation requires a frozen DofRegistry")
    if type(sample) is not AntiplaneShapeSample:
        raise ProviderContractError("provider must return AntiplaneShapeSample")
    if type(sample.status) is not ProviderStatus:
        raise ProviderContractError("sample.status must be a ProviderStatus")
    if sample.diagnostic is not None and type(sample.diagnostic) is not str:
        raise ProviderContractError("sample.diagnostic must be a string or None")

    if sample.status is not ProviderStatus.OK:
        if not sample.diagnostic:
            raise ProviderContractError("non-OK sample must preserve a diagnostic reason")
        if (
            type(sample.coefficient_ids) is not tuple
            or type(sample.contributions) is not tuple
            or len(sample.coefficient_ids) != 0
            or len(sample.contributions) != 0
        ):
            raise ProviderContractError("non-OK sample must not expose partial row data")
        return sample

    expected_ids = registry.canonical_ids
    if type(sample.coefficient_ids) is not tuple:
        raise ProviderContractError("coefficient_ids must be a tuple")
    if any(type(identifier) is not str for identifier in sample.coefficient_ids):
        raise ProviderContractError("coefficient_ids must contain canonical strings")
    if sample.coefficient_ids != expected_ids:
        raise ProviderContractError("coefficient_ids must exactly match frozen canonical registry order")
    if type(sample.contributions) is not tuple or len(sample.contributions) != registry.dof_count:
        raise ProviderContractError("contributions must contain one entry per registered DOF")
    for slot, contribution in enumerate(sample.contributions):
        if type(contribution) is not ScalarShapeContribution:
            raise ProviderContractError(f"contributions[{slot}] must be ScalarShapeContribution")
        if not _is_finite_real(contribution.value):
            raise ProviderContractError(f"contributions[{slot}].value must be finite real")
        if type(contribution.gradient) is not tuple or len(contribution.gradient) != 3:
            raise ProviderContractError(f"contributions[{slot}].gradient must be a 3-tuple")
        if any(not _is_finite_real(value) for value in contribution.gradient):
            raise ProviderContractError(f"contributions[{slot}].gradient must be finite real")
        if registry.ordered_keys[slot].field_component is not FieldComponent.UZ:
            raise ProviderContractError("antiplane provider only supports registered field_component=UZ")
    return sample


def antiplane_engineering_strain(gradient: Vector3) -> EngineeringStrain6:
    """Map ``grad(u_z)`` to [exx, eyy, ezz, gxy, gyz, gxz] convention.

    The three-vector is ``(d uz/dx, d uz/dy, d uz/dz)``. The return order and
    shear convention match the frozen 3-D elastic operator. In the ideal
    z-invariant Mode III field the third gradient entry vanishes analytically.
    """
    if type(gradient) is not tuple or len(gradient) != 3:
        raise ValueError("gradient must be a 3-tuple")
    if any(not _is_finite_real(value) for value in gradient):
        raise ValueError("gradient must contain finite real values")
    gx, gy, gz = (float(value) for value in gradient)
    return (0.0, 0.0, gz, 0.0, gy, gx)

