"""Source-only full-vector base/branch shape and 6-Voigt B algebra; no solver."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import json
import math
from numbers import Real
from typing import TypeAlias


Vector3: TypeAlias = tuple[float, float, float]
Matrix3: TypeAlias = tuple[Vector3, Vector3, Vector3]
BColumn6: TypeAlias = tuple[float, float, float, float, float, float]
MAX_SOURCE_COUNT = 27
DOFS_PER_SOURCE = 15
MAX_MIXED_DOFS = MAX_SOURCE_COUNT * DOFS_PER_SOURCE
FRAME_TOLERANCE = 1.0e-10
MIN_LENGTH_M = 1.0e-9
MAX_LENGTH_M = 1.0e3
COMPONENTS = ("ux", "uy", "uz")


class RowStatus(str, Enum):
    OK = "OK"
    UNKNOWN_FRONT = "UNKNOWN_FRONT"
    UNKNOWN_BRANCH_FACE = "UNKNOWN_BRANCH_FACE"


class EnrichmentContractError(ValueError):
    """Invalid source rows, SI inputs, crack frame, or branch selector."""


def _real(name: str, value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise EnrichmentContractError(f"{name} must be a finite real scalar")
    try:
        result = float(value)
    except (OverflowError, ValueError) as exc:
        raise EnrichmentContractError(f"{name} must be a finite real scalar") from exc
    if not math.isfinite(result):
        raise EnrichmentContractError(f"{name} must be a finite real scalar")
    return result


def _vector(name: str, value: object) -> Vector3:
    if type(value) not in (tuple, list) or len(value) != 3:
        raise EnrichmentContractError(f"{name} must contain exactly three scalars")
    return tuple(_real(f"{name}[{i}]", v) for i, v in enumerate(value))  # type: ignore[return-value]


def _matrix(name: str, value: object) -> Matrix3:
    if type(value) not in (tuple, list) or len(value) != 3:
        raise EnrichmentContractError(f"{name} must be 3 by 3")
    return tuple(_vector(f"{name}[{i}]", row) for i, row in enumerate(value))  # type: ignore[return-value]


def _transpose(a: Matrix3) -> Matrix3:
    return tuple(tuple(a[j][i] for j in range(3)) for i in range(3))  # type: ignore[return-value]


def _matvec(a: Matrix3, x: Vector3) -> Vector3:
    return tuple(sum(a[i][j] * x[j] for j in range(3)) for i in range(3))  # type: ignore[return-value]


def _det(a: Matrix3) -> float:
    return (a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
            - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
            + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]))


def _finite_tree(name: str, values: object) -> None:
    if isinstance(values, (tuple, list)):
        for value in values:
            _finite_tree(name, value)
    elif not math.isfinite(float(values)):
        raise FloatingPointError(f"non-finite {name}")


def mixed_dof_count(source_count: int) -> int:
    """Return 15 scalar coefficients per source: (BASE + 4 branches) x 3."""
    if type(source_count) is not int or not 1 <= source_count <= MAX_SOURCE_COUNT:
        raise EnrichmentContractError(f"source_count must be in [1, {MAX_SOURCE_COUNT}]")
    return DOFS_PER_SOURCE * source_count


def _proper_rotation(value: object) -> Matrix3:
    q = _matrix("axes_world_from_local", value)
    qtq = tuple(tuple(sum(q[k][i] * q[k][j] for k in range(3))
                      for j in range(3)) for i in range(3))
    if any(abs(qtq[i][j] - (1.0 if i == j else 0.0)) > FRAME_TOLERANCE
           for i in range(3) for j in range(3)):
        raise EnrichmentContractError("crack axes must be orthonormal")
    if abs(_det(q) - 1.0) > FRAME_TOLERANCE:
        raise EnrichmentContractError("crack axes must form a proper right-handed rotation")
    return q


def _branch_functions(x: float, y: float, reference_length_m: float,
                      branch_trace: str | None) -> tuple[float, tuple[float, ...], tuple[Vector3, ...], float] | None:
    r = math.hypot(x, y)
    if r == 0.0:
        return None
    on_cut = x < 0.0 and y == 0.0
    if on_cut and branch_trace is None:
        raise LookupError("UNKNOWN_BRANCH_FACE")
    if branch_trace is not None and not on_cut:
        raise EnrichmentContractError("branch_trace is valid only for local x < 0, y = 0")
    theta = (math.pi if branch_trace == "upper" else -math.pi) if on_cut else math.atan2(y, x)
    sh, ch = math.sin(theta / 2.0), math.cos(theta / 2.0)
    st, ct = math.sin(theta), math.cos(theta)
    h = (sh, ch, sh * st, ch * st)
    hp = (0.5 * ch, -0.5 * sh, 0.5 * ch * st + sh * ct,
          -0.5 * sh * st + ch * ct)
    q = math.sqrt(r) / math.sqrt(reference_length_m)
    scale = q / r
    values = tuple(q * item for item in h)
    gradients_local = tuple(
        (scale * (0.5 * ct * h[i] - st * hp[i]),
         scale * (0.5 * st * h[i] + ct * hp[i]), 0.0)
        for i in range(4)
    )
    _finite_tree("branch values/gradients", (values, gradients_local))
    return theta, values, gradients_local, r


def _semantic_id(source_id: str, component: str, basis_kind: str, basis_index: int) -> str:
    return json.dumps([source_id, component, basis_kind, basis_index],
                      ensure_ascii=True, separators=(",", ":"))


def _b_column(component: str, gradient: Vector3) -> BColumn6:
    gx, gy, gz = gradient
    if component == "ux":
        return (gx, 0.0, 0.0, gy, 0.0, gz)
    if component == "uy":
        return (0.0, gy, 0.0, gx, gz, 0.0)
    return (0.0, 0.0, gz, 0.0, gy, gx)


@dataclass(frozen=True, slots=True)
class ShapeContribution:
    semantic_id: str
    source_id: str
    component: str
    basis_kind: str
    basis_index: int
    shape_value: float
    shape_gradient_world_per_m: Vector3
    B_column_per_m: BColumn6


@dataclass(frozen=True, slots=True)
class EnrichedShapeRow:
    status: RowStatus
    source_ids: tuple[str, ...]
    coefficient_ids: tuple[str, ...]
    contributions: tuple[ShapeContribution, ...]
    branch_values: tuple[float, float, float, float] | None = None
    branch_gradients_world_per_m: tuple[Vector3, Vector3, Vector3, Vector3] | None = None
    theta_rad: float | None = None
    diagnostic: str | None = None


def build_mixed_enriched_row(
    point_world_m: object | None = None,
    *,
    source_ids: object,
    base_values: object,
    base_gradients_world_per_m: object,
    tip_world_m: object | None = None,
    axes_world_from_local: object,
    reference_length_m: float,
    branch_trace: str | None = None,
    local_point_from_tip_m: object | None = None,
) -> EnrichedShapeRow:
    """Build one row sorted by source, base/branch, then ux/uy/uz.

    Enrichment is ``phi=N F`` with ``grad(phi)=F grad(N)+N grad(F)``.
    Supply world point plus tip, or exact local coordinates relative to tip;
    axes are Q columns and the local form preserves exact crack-face traces.
    """
    if type(source_ids) is not tuple or not 1 <= len(source_ids) <= MAX_SOURCE_COUNT:
        raise EnrichmentContractError(f"source_ids must be a tuple of 1..{MAX_SOURCE_COUNT} IDs")
    if any(type(s) is not str or not s or s.strip() != s for s in source_ids):
        raise EnrichmentContractError("source IDs must be nonempty canonical strings")
    if len(set(source_ids)) != len(source_ids):
        raise EnrichmentContractError("source IDs must be unique")
    if type(base_values) not in (tuple, list) or len(base_values) != len(source_ids):
        raise EnrichmentContractError("base_values must have one row per source ID")
    if (type(base_gradients_world_per_m) not in (tuple, list)
            or len(base_gradients_world_per_m) != len(source_ids)):
        raise EnrichmentContractError("base gradients must have one row per source ID")
    qmat = _proper_rotation(axes_world_from_local)
    ref = _real("reference_length_m", reference_length_m)
    if not MIN_LENGTH_M <= ref <= MAX_LENGTH_M:
        raise EnrichmentContractError("reference_length_m is outside the finite SI contract")
    if branch_trace not in (None, "upper", "lower"):
        raise EnrichmentContractError("branch_trace must be None, 'upper', or 'lower'")
    if local_point_from_tip_m is not None:
        if point_world_m is not None or tip_world_m is not None:
            raise EnrichmentContractError("supply local coordinates or world point plus tip, not both")
        local = _vector("local_point_from_tip_m", local_point_from_tip_m)
    else:
        if point_world_m is None or tip_world_m is None:
            raise EnrichmentContractError("world queries require both point_world_m and tip_world_m")
        point = _vector("point_world_m", point_world_m)
        tip = _vector("tip_world_m", tip_world_m)
        delta = tuple(point[i] - tip[i] for i in range(3))
        local = _matvec(_transpose(qmat), delta)
    values_in = tuple(_real(f"base_values[{i}]", value) for i, value in enumerate(base_values))
    gradients_in = tuple(_vector(f"base_gradients_world_per_m[{i}]", row)
                         for i, row in enumerate(base_gradients_world_per_m))
    order = tuple(sorted(range(len(source_ids)), key=lambda i: source_ids[i]))
    ordered_ids = tuple(source_ids[i] for i in order)
    _finite_tree("local coordinate", local)
    try:
        branch_result = _branch_functions(local[0], local[1], ref, branch_trace)
    except LookupError:
        return EnrichedShapeRow(RowStatus.UNKNOWN_BRANCH_FACE, ordered_ids, (), (),
                                diagnostic="explicit upper/lower branch trace required")
    if branch_result is None:
        return EnrichedShapeRow(RowStatus.UNKNOWN_FRONT, ordered_ids, (), (),
                                diagnostic="branch gradients are singular at r=0")
    theta, branch_values, gradients_local, _ = branch_result
    gradients_world = tuple(_matvec(qmat, grad) for grad in gradients_local)
    _finite_tree("world branch gradients", gradients_world)
    rows: list[ShapeContribution] = []
    for source_index in order:
        source_id = source_ids[source_index]
        base = values_in[source_index]
        base_grad = gradients_in[source_index]
        scalar_rows = [("base", 0, base, base_grad)]
        for branch_index in range(4):
            f = branch_values[branch_index]
            gf = gradients_world[branch_index]
            phi = base * f
            grad_phi = tuple(f * base_grad[j] + base * gf[j] for j in range(3))
            scalar_rows.append(("front_branch", branch_index + 1, phi, grad_phi))
        for basis_kind, basis_index, shape_value, shape_gradient in scalar_rows:
            _finite_tree("enriched shape", (shape_value, shape_gradient))
            for component in COMPONENTS:
                rows.append(ShapeContribution(
                    semantic_id=_semantic_id(source_id, component, basis_kind, basis_index),
                    source_id=source_id, component=component, basis_kind=basis_kind,
                    basis_index=basis_index, shape_value=shape_value,
                    shape_gradient_world_per_m=shape_gradient,
                    B_column_per_m=_b_column(component, shape_gradient),
                ))
    coefficient_ids = tuple(row.semantic_id for row in rows)
    if len(rows) != mixed_dof_count(len(source_ids)) or len(set(coefficient_ids)) != len(rows):
        raise EnrichmentContractError("semantic mixed DOF row is incomplete or non-unique")
    return EnrichedShapeRow(
        RowStatus.OK, ordered_ids, coefficient_ids, tuple(rows),
        branch_values=branch_values, branch_gradients_world_per_m=gradients_world,
        theta_rad=theta,
    )
