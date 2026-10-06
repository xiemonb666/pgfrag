"""Mode-I cohesive interface residual/tangent assembly; source-only primitive.

Each quadrature point stores positive/negative face interpolation operators
H+ and H- (3 x ndof); its explicit jump operator is J = H+ - H-. With nodal
trial displacement d in metres and unit normal n, jump=Jd and delta_n=n.Jd.
Each point contributes its one-side area weight once. Units are: J,n dimensionless,
area m^2, traction Pa, residual N, tangent N/m, and integrated energies J.

Every call evaluates the trial opening from the same caller-owned immutable
committed history. No Newton trial commits or mutates material history. Scalar
damage/energy behavior comes from the existing cohesive scalar module. Negative
softening tangents are valid; no SPD guarantee is made. Compression, tangential
slip and mixed mode are UNKNOWN and rejected. Kink tangents are labeled and
report their selected one-sided branch. This is not a solver or front coupling.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from enum import Enum
import importlib.util
import math
from numbers import Real
from pathlib import Path
import sys


_SCALAR_NAME = "_pgfrag_cohesive_scalar_law_20261006"
_SCALAR_PATH = Path(__file__).resolve().parents[1] / "cohesive_scalar_law_20261006" / "cohesive_scalar_law.py"
_SPEC = importlib.util.spec_from_file_location(_SCALAR_NAME, _SCALAR_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise ImportError(f"cannot load existing scalar law: {_SCALAR_PATH}")
_SCALAR = sys.modules.get(_SCALAR_NAME) or importlib.util.module_from_spec(_SPEC)
if _SCALAR_NAME not in sys.modules:
    sys.modules[_SCALAR_NAME] = _SCALAR
    _SPEC.loader.exec_module(_SCALAR)
BilinearCohesiveMaterial = _SCALAR.BilinearCohesiveMaterial
CohesiveState = _SCALAR.CohesiveState
_scalar_update, _scalar_traction, _MODE_I = _SCALAR.update_state, _SCALAR.traction_pa, _SCALAR.MODE_I_OPENING

PointId = tuple[str, str]
Vector3 = tuple[float, float, float]
Matrix3 = tuple[tuple[float, ...], tuple[float, ...], tuple[float, ...]]
MAX_NODES, MAX_DOF, MAX_QUADRATURE_POINTS = 12, 36, 16
_GEOM_ABS_TOL_M, _GEOM_REL_TOL = 1.0e-12, 1.0e-10
_FORCE_ABS_TOL_N, _FORCE_REL_TOL = 1.0e-12, 1.0e-10
_KINK_ABS, _KINK_REL = 1.0e-15, 1.0e-12


class UnsupportedInterfaceModeError(NotImplementedError):
    """Compression or any tangential/mixed-mode opening is UNKNOWN here."""


class TangentStatus(str, Enum):
    DIFFERENTIABLE = "differentiable"
    NONDIFFERENTIABLE = "nondifferentiable"


def _real(name: str, value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a plain real scalar in SI units")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _v3(name: str, values: Sequence[object]) -> Vector3:
    if not isinstance(values, (tuple, list)) or len(values) != 3:
        raise ValueError(f"{name} must contain three SI values")
    return tuple(_real(name, value) for value in values)  # type: ignore[return-value]


def _m3(name: str, values: Sequence[Sequence[object]]) -> Matrix3:
    if not isinstance(values, (tuple, list)) or len(values) != 3:
        raise ValueError(f"{name} must be 3 x ndof")
    rows = tuple(tuple(_real(name, value) for value in row) for row in values)
    if not rows[0] or any(len(row) != len(rows[0]) for row in rows):
        raise ValueError(f"{name} rows need the same positive ndof")
    return rows  # type: ignore[return-value]


def _dot(a: Sequence[float], b: Sequence[float]) -> float: return sum(x * y for x, y in zip(a, b))


def _apply(a: Matrix3, x: Sequence[float]) -> Vector3: return tuple(_dot(row, x) for row in a)  # type: ignore[return-value]


def _cross(a: Vector3, b: Vector3) -> Vector3:
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def _close(a: float, b: float) -> bool: return abs(a-b) <= _GEOM_ABS_TOL_M + _GEOM_REL_TOL * max(abs(a), abs(b))


def _force_close_n(a: float, b: float) -> bool:
    return abs(a-b) <= _FORCE_ABS_TOL_N + _FORCE_REL_TOL * max(abs(a), abs(b))


@dataclass(frozen=True)
class CohesiveQuadraturePoint:
    """One physical point; both face operators use one global DOF ordering."""

    interface_id: str
    quadrature_id: str
    positive_face_operator: Matrix3
    negative_face_operator: Matrix3
    normal: Vector3
    positive_face_point_m: Vector3
    negative_face_point_m: Vector3
    one_side_area_weight_m2: float

    def __post_init__(self) -> None:
        if not isinstance(self.interface_id, str) or not self.interface_id.strip():
            raise ValueError("interface_id must be nonempty")
        if not isinstance(self.quadrature_id, str) or not self.quadrature_id.strip():
            raise ValueError("quadrature_id must be nonempty")
        plus, minus = _m3("positive_face_operator", self.positive_face_operator), _m3("negative_face_operator", self.negative_face_operator)
        if len(plus[0]) != len(minus[0]):
            raise ValueError("face operators must share ndof")
        normal = _v3("normal", self.normal)
        if not math.isclose(math.sqrt(_dot(normal, normal)), 1.0, rel_tol=0.0, abs_tol=1.0e-12):
            raise ValueError("normal must be unit length; normalization is not implicit")
        area = _real("one_side_area_weight_m2", self.one_side_area_weight_m2)
        if area <= 0.0:
            raise ValueError("one_side_area_weight_m2 must be positive")
        object.__setattr__(self, "positive_face_operator", plus)
        object.__setattr__(self, "negative_face_operator", minus)
        object.__setattr__(self, "normal", normal)
        object.__setattr__(self, "positive_face_point_m", _v3("positive_face_point_m", self.positive_face_point_m))
        object.__setattr__(self, "negative_face_point_m", _v3("negative_face_point_m", self.negative_face_point_m))
        object.__setattr__(self, "one_side_area_weight_m2", area)

    @property
    def point_id(self) -> PointId:
        return self.interface_id, self.quadrature_id

    @property
    def ndof(self) -> int:
        return len(self.positive_face_operator[0])

    @property
    def jump_operator(self) -> Matrix3:
        return tuple(tuple(p-m for p, m in zip(a, b)) for a, b in zip(self.positive_face_operator, self.negative_face_operator))  # type: ignore[return-value]


@dataclass(frozen=True)
class PointTrialResult:
    point_id: PointId
    opening_m: float
    traction_pa: float
    k_alg_pa_per_m: float
    tangent_status: TangentStatus
    tangent_side: str
    nondifferentiable_reasons: tuple[str, ...]
    trial_state: CohesiveState
    positive_face_force_n: tuple[float, ...]
    negative_face_force_n: tuple[float, ...]
    positive_face_resultant_n: Vector3
    negative_face_resultant_n: Vector3


@dataclass(frozen=True)
class InterfaceTrialResult:
    internal_force_n: tuple[float, ...]
    tangent_n_per_m: tuple[tuple[float, ...], ...]
    point_results: tuple[PointTrialResult, ...]
    stored_energy_j: float
    dissipated_energy_j: float
    external_work_j: float
    tangent_status: TangentStatus
    spd_guarantee: str = "none; cohesive softening may make tangent indefinite"


def _validate_geometry(point: CohesiveQuadraturePoint, positions: tuple[Vector3, ...]) -> None:
    ndof = 3 * len(positions)
    if point.ndof != ndof:
        raise ValueError(f"{point.point_id}: J must have 3*node_count columns")
    plus, minus, jump = point.positive_face_operator, point.negative_face_operator, point.jump_operator
    flat = tuple(x for p in positions for x in p)
    for axis in range(3):
        translation = tuple(float(i % 3 == axis) for i in range(ndof))
        expected = tuple(float(i == axis) for i in range(3))
        if any(not _close(x, y) for op in (plus, minus) for x, y in zip(_apply(op, translation), expected)):
            raise ValueError(f"{point.point_id}: each face operator must reproduce rigid translations")
        if any(not _close(x, 0.0) for x in _apply(jump, translation)):
            raise ValueError(f"{point.point_id}: J must annihilate rigid translations")
    for op, declared in ((plus, point.positive_face_point_m), (minus, point.negative_face_point_m)):
        if any(not _close(x, y) for x, y in zip(_apply(op, flat), declared)):
            raise ValueError(f"{point.point_id}: face operator disagrees with its geometric point")
    if any(not _close(a, b) for a, b in zip(point.positive_face_point_m, point.negative_face_point_m)):
        raise ValueError(f"{point.point_id}: paired face geometric points must coincide")
    origin = tuple(sum(p[i] for p in positions)/len(positions) for i in range(3))
    for axis in ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)):
        rotation = tuple(x for p in positions for x in _cross(axis, tuple(p[i]-origin[i] for i in range(3))))
        scale = max(1.0e-12, *(abs(x) for x in rotation))
        if any(abs(x) > 1.0e-15+1.0e-10*scale for x in _apply(jump, rotation)):
            raise ValueError(f"{point.point_id}: J must annihilate rigid rotations at supplied node points")


def _near(a: float, b: float) -> bool:
    return abs(a-b) <= _KINK_ABS + _KINK_REL*max(abs(a), abs(b))


def _algorithmic_tangent(material, committed: CohesiveState, opening: float):
    d0, df, history = material.delta0_m, material.deltaf_m, committed.kappa_m
    reasons = tuple(name for name, corner in (("history maximum", history), ("delta0", d0), ("deltaf", df)) if _near(opening, corner))
    if opening <= history:
        tangent, side = (1.0-committed.damage)*material.k0_pa_per_m, "history-fixed/unloading-side"
    elif opening <= d0:
        tangent, side = material.k0_pa_per_m, "elastic-left-side"
    elif opening < df:
        tangent, side = -material.tpeak_pa/(df-d0), "softening-loading-side"
    else:
        tangent, side = 0.0, "fully-broken-right-side"
    status = TangentStatus.NONDIFFERENTIABLE if reasons else TangentStatus.DIFFERENTIABLE
    return tangent, status, side, reasons


def _allocate_tangent(ndof: int) -> list[list[float]]:
    return [[0.0]*ndof for _ in range(ndof)]


def assemble_trial(
    material: BilinearCohesiveMaterial,
    points: Sequence[CohesiveQuadraturePoint],
    nodal_positions_m: Sequence[Sequence[object]],
    trial_displacement_m: Sequence[object],
    committed_states: Mapping[PointId, CohesiveState],
) -> InterfaceTrialResult:
    """Read-only Newton trial: f=sum(w J.T n t), K=sum(w g k_alg g.T), g=J.T n."""
    if not isinstance(material, BilinearCohesiveMaterial):
        raise TypeError("material must be the existing BilinearCohesiveMaterial")
    if not isinstance(points, (tuple, list)) or not points:
        raise ValueError("points must be a nonempty sequence")
    if len(points) > MAX_QUADRATURE_POINTS:
        raise ValueError(f"quadrature point cap is {MAX_QUADRATURE_POINTS}")
    points = tuple(points)
    if any(not isinstance(p, CohesiveQuadraturePoint) for p in points):
        raise TypeError("all points must be CohesiveQuadraturePoint")
    ids = tuple(p.point_id for p in points)
    if len(set(ids)) != len(ids):
        raise ValueError("duplicate point ID would double-count interface area/history")
    if not isinstance(committed_states, Mapping) or set(committed_states) != set(ids):
        raise ValueError("committed history IDs must exactly match interface/quadrature IDs")
    if not isinstance(nodal_positions_m, (tuple, list)) or not nodal_positions_m:
        raise ValueError("nodal_positions_m must be a nonempty sequence")
    node_count = len(nodal_positions_m)
    ndof = 3*node_count
    if node_count > MAX_NODES or ndof > MAX_DOF:
        raise ValueError(f"nodal primitive cap is {MAX_NODES} nodes / {MAX_DOF} DOFs")
    states = dict(committed_states)
    positions = tuple(_v3("nodal_positions_m", row) for row in nodal_positions_m)
    if any(p.ndof != ndof for p in points):
        raise ValueError("all J operators need 3*node_count columns")
    if not isinstance(trial_displacement_m, (tuple, list)) or len(trial_displacement_m) != ndof:
        raise ValueError("trial_displacement_m must have 3*node_count entries")
    displacement = tuple(_real("trial_displacement_m", x) for x in trial_displacement_m)
    for point in points:
        _validate_geometry(point, positions)

    force, tangent = [0.0]*ndof, _allocate_tangent(ndof)
    results, stored, dissipated, work = [], 0.0, 0.0, 0.0
    for point in points:
        committed = states[point.point_id]
        jump = _apply(point.jump_operator, displacement)
        opening = _dot(point.normal, jump)
        if opening < 0.0:
            raise UnsupportedInterfaceModeError(f"{point.point_id}: compression is UNKNOWN")
        tangential = tuple(jump[i]-opening*point.normal[i] for i in range(3))
        if math.sqrt(_dot(tangential, tangential)) > 1.0e-12 + 1.0e-10*math.sqrt(_dot(jump, jump)):
            raise UnsupportedInterfaceModeError(f"{point.point_id}: tangential/mixed mode is UNKNOWN")
        trial = _scalar_update(material, committed, opening, mode=_MODE_I)
        traction = _scalar_traction(material, trial, mode=_MODE_I)
        kalg, status, side, reasons = _algorithmic_tangent(material, committed, opening)
        area, plus_force, minus_force = point.one_side_area_weight_m2, [0.0]*ndof, [0.0]*ndof
        gradient = tuple(sum(point.jump_operator[c][i]*point.normal[c] for c in range(3)) for i in range(ndof))
        for i in range(ndof):
            for c in range(3):
                traction_weight = point.normal[c]*traction*area
                plus_force[i] += point.positive_face_operator[c][i]*traction_weight
                minus_force[i] -= point.negative_face_operator[c][i]*traction_weight
            force[i] += plus_force[i]+minus_force[i]
        plus_resultant = tuple(sum(plus_force[3*n+c] for n in range(ndof//3)) for c in range(3))
        minus_resultant = tuple(sum(minus_force[3*n+c] for n in range(ndof//3)) for c in range(3))
        if any(not _force_close_n(a, -b) for a, b in zip(plus_resultant, minus_resultant)):
            raise ValueError(f"{point.point_id}: face forces do not balance action and reaction")
        for i in range(ndof):
            for j in range(ndof):
                tangent[i][j] += area*gradient[i]*kalg*gradient[j]
        stored += area*trial.stored_energy_j_per_m2
        dissipated += area*trial.dissipated_energy_j_per_m2
        work += area*trial.external_work_j_per_m2
        results.append(PointTrialResult(point.point_id, opening, traction, kalg, status, side, reasons, trial, tuple(plus_force), tuple(minus_force), plus_resultant, minus_resultant))
    if not all(math.isfinite(x) for x in (*force, stored, dissipated, work)) or any(not math.isfinite(x) for row in tangent for x in row):
        raise FloatingPointError("interface force, tangent, or energy overflowed")
    status = TangentStatus.NONDIFFERENTIABLE if any(r.tangent_status is TangentStatus.NONDIFFERENTIABLE for r in results) else TangentStatus.DIFFERENTIABLE
    return InterfaceTrialResult(tuple(force), tuple(tuple(row) for row in tangent), tuple(results), stored, dissipated, work, status)
