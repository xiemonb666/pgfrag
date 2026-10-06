"""Bounded enriched-face Mode-I trial primitive; source-only experimental interface.

All generalized coefficients are displacement amplitudes in metres and both
face operators are dimensionless. A caller must supply six semantic rigid
candidate vectors and their expected displacements at every shared reference
point. This verifies only those declared candidates, not the actual 93-DOF
enriched spectrum. The primitive is NOT_COUPLED to a crack front or solver.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from enum import Enum
import hashlib
import importlib.util
import math
from numbers import Real
from pathlib import Path
import sys


MAX_COEFFICIENTS, MAX_QUADRATURE_POINTS = 93, 16
GEOMETRY_ABS_TOL_M, GEOMETRY_REL_TOL = 1.0e-12, 1.0e-10
NORMAL_TOL, MODE_ABS_TOL_M, MODE_REL_TOL = 1.0e-12, 1.0e-12, 1.0e-10
SCALAR_SOURCE_PIN = "78433638608ad4216979607412fe0d4966f1423706c47f734f50dd2a8121fa17"
DERIVATIVE_SOURCE_PIN = "d2bf0dd79507e5ab26efe69b5489c3284b7c9118656ab49fffe942303ba11854"
RIGID_MODE_NAMES = ("Tx", "Ty", "Tz", "Rx", "Ry", "Rz")

_HERE = Path(__file__).resolve().parent
_EXPERIMENTS = _HERE.parent
_SCALAR_PATH = _EXPERIMENTS / "cohesive_scalar_law_20261006" / "cohesive_scalar_law.py"
_DERIVATIVE_PATH = _EXPERIMENTS / "cohesive_interface_assembly_20261006" / "cohesive_interface_assembly.py"
_DERIVATIVE_MODULE = "_pgfrag_cohesive_interface_assembly_20261006"


def _load_pinned_derivative():
    scalar_hash = hashlib.sha256(_SCALAR_PATH.read_bytes()).hexdigest()
    derivative_hash = hashlib.sha256(_DERIVATIVE_PATH.read_bytes()).hexdigest()
    if scalar_hash != SCALAR_SOURCE_PIN or derivative_hash != DERIVATIVE_SOURCE_PIN:
        raise ImportError("frozen scalar law or reviewed nodal derivative pin changed")
    spec = importlib.util.spec_from_file_location(_DERIVATIVE_MODULE, _DERIVATIVE_PATH)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load pinned derivative source: {_DERIVATIVE_PATH}")
    module = sys.modules.get(_DERIVATIVE_MODULE)
    if module is None:
        module = importlib.util.module_from_spec(spec)
        sys.modules[_DERIVATIVE_MODULE] = module
        spec.loader.exec_module(module)
    return module


_DERIVATIVE = _load_pinned_derivative()
BilinearCohesiveMaterial = _DERIVATIVE.BilinearCohesiveMaterial
CohesiveState = _DERIVATIVE.CohesiveState
TangentStatus = _DERIVATIVE.TangentStatus

PointId = tuple[str, str]
Vector3 = tuple[float, float, float]
Matrix3 = tuple[tuple[float, ...], tuple[float, ...], tuple[float, ...]]


class UnsupportedEnrichedModeError(NotImplementedError):
    """Compression and tangential/mixed-mode enriched responses are UNKNOWN."""


class TraceFamily(str, Enum):
    BASE = "BASE"
    GLOBALFRONT = "GLOBALFRONT"


def _real(name: str, value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a plain real value in its named SI unit")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _positive(name: str, value: object) -> float:
    result = _real(name, value)
    if result <= 0.0:
        raise ValueError(f"{name} must be positive")
    return result


def _ids(value: object) -> tuple[str, ...]:
    if not isinstance(value, tuple) or not 1 <= len(value) <= MAX_COEFFICIENTS:
        raise ValueError("coefficient_ids must be a frozen tuple of 1..93 semantic IDs")
    if any(not isinstance(item, str) or not item.strip() for item in value) or len(set(value)) != len(value):
        raise ValueError("semantic coefficient IDs must be unique nonempty strings")
    return value


def _v3(name: str, values: Sequence[object]) -> Vector3:
    if not isinstance(values, (tuple, list)) or len(values) != 3:
        raise ValueError(f"{name} must contain three SI values")
    return tuple(_real(name, item) for item in values)  # type: ignore[return-value]


def _m3(name: str, values: Sequence[Sequence[object]], ndof: int) -> Matrix3:
    if not isinstance(values, (tuple, list)) or len(values) != 3:
        raise ValueError(f"{name} must have shape 3 x ndof")
    rows = tuple(tuple(_real(name, item) for item in row) for row in values)
    if any(len(row) != ndof for row in rows):
        raise ValueError(f"{name} columns must match the frozen semantic coefficient order")
    return rows  # type: ignore[return-value]


def _dot(left: Sequence[float], right: Sequence[float]) -> float:
    return sum(a * b for a, b in zip(left, right))


def _apply(matrix: Matrix3, vector: Sequence[float]) -> Vector3:
    return tuple(_dot(row, vector) for row in matrix)  # type: ignore[return-value]


def _geom_close_m(actual: float, expected: float) -> bool:
    return abs(actual-expected) <= GEOMETRY_ABS_TOL_M + GEOMETRY_REL_TOL*max(abs(actual), abs(expected))


def _vclose_m(actual: Vector3, expected: Vector3) -> bool:
    return all(_geom_close_m(a, b) for a, b in zip(actual, expected))


@dataclass(frozen=True)
class RigidTraceCandidates:
    """Caller-declared Tx,Ty,Tz,Rx,Ry,Rz coefficient vectors; amplitudes are m."""

    coefficient_ids: tuple[str, ...]
    candidate_vectors_m: tuple[tuple[float, ...], ...]

    def __post_init__(self) -> None:
        ids = _ids(self.coefficient_ids)
        if not isinstance(self.candidate_vectors_m, tuple) or len(self.candidate_vectors_m) != 6:
            raise ValueError("exactly six ordered rigid candidate vectors are required")
        vectors = tuple(tuple(_real("candidate_vectors_m", x) for x in row) for row in self.candidate_vectors_m)
        if any(len(row) != len(ids) for row in vectors):
            raise ValueError("every rigid candidate vector must match coefficient_ids")
        object.__setattr__(self, "candidate_vectors_m", vectors)


@dataclass(frozen=True)
class EnrichedFacePoint:
    """One QP with paired face traces evaluated at one common world point."""

    interface_id: str
    quadrature_id: str
    basis_family: TraceFamily | str
    coefficient_ids: tuple[str, ...]
    positive_face_operator: Matrix3
    negative_face_operator: Matrix3
    normal: Vector3
    positive_reference_point_m: Vector3
    negative_reference_point_m: Vector3
    expected_rigid_displacements_m: tuple[Vector3, ...]
    one_side_area_weight_m2: float

    def __post_init__(self) -> None:
        if any(not isinstance(x, str) or not x.strip() for x in (self.interface_id, self.quadrature_id)):
            raise ValueError("interface_id and quadrature_id must be nonempty strings")
        ids = _ids(self.coefficient_ids)
        try:
            family = TraceFamily(self.basis_family)
        except (TypeError, ValueError) as exc:
            raise ValueError("basis_family must be BASE or GLOBALFRONT") from exc
        plus = _m3("positive_face_operator", self.positive_face_operator, len(ids))
        minus = _m3("negative_face_operator", self.negative_face_operator, len(ids))
        normal = _v3("normal", self.normal)
        if abs(math.sqrt(_dot(normal, normal))-1.0) > NORMAL_TOL:
            raise ValueError("normal must be unit length; it is never silently normalized")
        plus_point = _v3("positive_reference_point_m", self.positive_reference_point_m)
        minus_point = _v3("negative_reference_point_m", self.negative_reference_point_m)
        if not _vclose_m(plus_point, minus_point):
            raise ValueError("H+ and H- must refer to the same world reference point")
        if not isinstance(self.expected_rigid_displacements_m, tuple) or len(self.expected_rigid_displacements_m) != 6:
            raise ValueError("six expected point displacements in Tx..Rz order are required")
        rigid = tuple(_v3("expected_rigid_displacements_m", row) for row in self.expected_rigid_displacements_m)
        object.__setattr__(self, "basis_family", family)
        object.__setattr__(self, "positive_face_operator", plus)
        object.__setattr__(self, "negative_face_operator", minus)
        object.__setattr__(self, "normal", normal)
        object.__setattr__(self, "positive_reference_point_m", plus_point)
        object.__setattr__(self, "negative_reference_point_m", minus_point)
        object.__setattr__(self, "expected_rigid_displacements_m", rigid)
        object.__setattr__(self, "one_side_area_weight_m2", _positive("one_side_area_weight_m2", self.one_side_area_weight_m2))

    @property
    def point_id(self) -> PointId:
        return (self.interface_id, self.quadrature_id)

    @property
    def jump_operator(self) -> Matrix3:
        return tuple(tuple(a-b for a, b in zip(plus, minus)) for plus, minus in zip(self.positive_face_operator, self.negative_face_operator))  # type: ignore[return-value]


def _validate_rigid_trace(point: EnrichedFacePoint, candidates: RigidTraceCandidates) -> None:
    if candidates.coefficient_ids != point.coefficient_ids:
        raise ValueError("rigid candidates must use the exact frozen coefficient order")
    jump = point.jump_operator
    for name, candidate, expected in zip(RIGID_MODE_NAMES, candidates.candidate_vectors_m, point.expected_rigid_displacements_m):
        if not _vclose_m(_apply(point.positive_face_operator, candidate), expected):
            raise ValueError(f"H+ fails caller-supplied rigid candidate {name} at the reference point")
        if not _vclose_m(_apply(point.negative_face_operator, candidate), expected):
            raise ValueError(f"H- fails caller-supplied rigid candidate {name} at the reference point")
        if not _vclose_m(_apply(jump, candidate), (0.0, 0.0, 0.0)):
            raise ValueError(f"J does not annihilate caller-supplied rigid candidate {name}")


def _allocate_tangent(ndof: int) -> list[list[float]]:
    return [[0.0]*ndof for _ in range(ndof)]


@dataclass(frozen=True)
class PointTrialResult:
    point_id: PointId
    opening_m: float
    traction_pa: float
    k_alg_pa_per_m: float
    tangent_side: str
    tangent_status: TangentStatus
    nondifferentiable_reasons: tuple[str, ...]
    trial_state: CohesiveState
    generalized_force_n: tuple[float, ...]
    positive_face_force_at_reference_n: Vector3
    negative_face_force_at_reference_n: Vector3


@dataclass(frozen=True)
class EnrichedTrialResult:
    generalized_internal_force_n: tuple[float, ...]
    tangent_n_per_m: tuple[tuple[float, ...], ...]
    point_results: tuple[PointTrialResult, ...]
    stored_energy_j: float
    dissipated_energy_j: float
    external_work_j: float
    tangent_status: TangentStatus
    rigid_evidence_status: str = "all six caller-declared modes reproduced; not a spectral proof"
    spd_guarantee: str = "none; cohesive softening may make the tangent indefinite"


def assemble_enriched_trial(
    material: BilinearCohesiveMaterial,
    points: Sequence[EnrichedFacePoint],
    trial_coefficients_m: Sequence[object],
    committed_states: Mapping[PointId, CohesiveState],
    rigid_candidates: RigidTraceCandidates,
) -> EnrichedTrialResult:
    """Read-only Mode-I trial: f=A J.T n t; K=A (J.T n) k_alg (J.T n).T."""
    if not isinstance(material, BilinearCohesiveMaterial):
        raise TypeError("material must be the frozen scalar-law BilinearCohesiveMaterial")
    if not isinstance(points, (tuple, list)) or not points:
        raise ValueError("points must be a nonempty sequence")
    if len(points) > MAX_QUADRATURE_POINTS:
        raise ValueError("quadrature-point cap is 16")
    points = tuple(points)
    if any(not isinstance(point, EnrichedFacePoint) for point in points):
        raise TypeError("all points must be EnrichedFacePoint")
    point_ids = tuple(point.point_id for point in points)
    if len(set(point_ids)) != len(point_ids):
        raise ValueError("duplicate point IDs would double-count area and history")
    coefficient_ids = points[0].coefficient_ids
    ndof = len(coefficient_ids)
    if ndof > MAX_COEFFICIENTS:
        raise ValueError("enriched coefficient cap is 93")
    if not isinstance(rigid_candidates, RigidTraceCandidates):
        raise ValueError("six rigid candidate vectors are required; rigid qualification cannot be skipped")
    if rigid_candidates.coefficient_ids != coefficient_ids or any(p.coefficient_ids != coefficient_ids for p in points):
        raise ValueError("all points and rigid evidence must share one frozen semantic coefficient order")
    if any(p.basis_family != points[0].basis_family for p in points):
        raise ValueError("one assembly may contain only one declared BASE or GLOBALFRONT family")
    if not isinstance(committed_states, Mapping) or set(committed_states) != set(point_ids):
        raise ValueError("committed history IDs must exactly match unique interface/QP IDs")
    if any(not isinstance(state, CohesiveState) for state in committed_states.values()):
        raise TypeError("committed states must use the pinned immutable scalar state type")
    if not isinstance(trial_coefficients_m, (tuple, list)) or len(trial_coefficients_m) != ndof:
        raise ValueError("trial_coefficients_m must match the frozen semantic coefficient order")
    trial = tuple(_real("trial_coefficients_m", x) for x in trial_coefficients_m)
    for point in points:
        _validate_rigid_trace(point, rigid_candidates)

    force, tangent = [0.0]*ndof, _allocate_tangent(ndof)
    results, stored, dissipated, work = [], 0.0, 0.0, 0.0
    for point in points:
        jump_op = point.jump_operator
        jump = _apply(jump_op, trial)
        opening = _dot(point.normal, jump)
        if opening < 0.0:
            raise UnsupportedEnrichedModeError("compression response is UNKNOWN and unsupported")
        tangential = tuple(jump[i]-opening*point.normal[i] for i in range(3))
        jump_norm = math.sqrt(_dot(jump, jump))
        tangent_norm = math.sqrt(_dot(tangential, tangential))
        if tangent_norm > MODE_ABS_TOL_M+MODE_REL_TOL*jump_norm:
            raise UnsupportedEnrichedModeError("tangential slip/mixed-mode response is UNKNOWN and unsupported")
        committed = committed_states[point.point_id]
        state = _DERIVATIVE._scalar_update(material, committed, opening, mode=_DERIVATIVE._MODE_I)
        traction = _DERIVATIVE._scalar_traction(material, state, mode=_DERIVATIVE._MODE_I)
        kalg, status, side, reasons = _DERIVATIVE._algorithmic_tangent(material, committed, opening)
        gradient = tuple(sum(jump_op[c][i]*point.normal[c] for c in range(3)) for i in range(ndof))
        area = point.one_side_area_weight_m2
        point_force = tuple(area*gradient[i]*traction for i in range(ndof))
        for i in range(ndof):
            force[i] += point_force[i]
            for j in range(ndof):
                tangent[i][j] += area*gradient[i]*kalg*gradient[j]
        face_force = tuple(area*traction*x for x in point.normal)
        results.append(PointTrialResult(point.point_id, opening, traction, kalg, side, status, reasons, state, point_force, face_force, tuple(-x for x in face_force)))
        stored += area*state.stored_energy_j_per_m2
        dissipated += area*state.dissipated_energy_j_per_m2
        work += area*state.external_work_j_per_m2
    if not all(math.isfinite(x) for x in (*force, stored, dissipated, work)) or any(not math.isfinite(x) for row in tangent for x in row):
        raise FloatingPointError("enriched force, tangent, or energy overflowed")
    overall = TangentStatus.NONDIFFERENTIABLE if any(p.tangent_status is TangentStatus.NONDIFFERENTIABLE for p in results) else TangentStatus.DIFFERENTIABLE
    return EnrichedTrialResult(tuple(force), tuple(tuple(row) for row in tangent), tuple(results), stored, dissipated, work, overall)
