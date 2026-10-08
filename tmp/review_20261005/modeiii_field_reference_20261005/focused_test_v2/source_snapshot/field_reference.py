"""Bounded exact Mode III reference field for a straight crack front.

This is an analytic manufactured-field primitive, not a solver, material
calibration, or finite-front crack solution. Coordinates and lengths are SI.
The only branch cut inside the admissible cylinder is the prescribed crack
surface x < 0, y = 0, -L/2 <= z <= L/2. Exact front and branch-face queries
are UNKNOWN because the field is singular or has two traces there.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import math

import numpy as np


MAX_QUERIES = 256
MAX_OUTPUT_BYTES = 2 * 1024 * 1024
MIN_LENGTH_M = 1.0e-9
MAX_LENGTH_M = 1.0e3
MIN_YOUNG_MODULUS_PA = 1.0e-3
MAX_YOUNG_MODULUS_PA = 1.0e15
MAX_ABS_KIII_PA_SQRT_M = 1.0e12
NORMAL_UNIT_TOLERANCE = 1.0e-12
_SQRT_TWO_OVER_PI = math.sqrt(2.0 / math.pi)


class FieldStatus(str, Enum):
    """Status for one point; placeholder zeros are unusable unless status is OK."""

    OK = "OK"
    OUTSIDE_DOMAIN = "OUTSIDE_DOMAIN"
    UNKNOWN_FRONT = "UNKNOWN_FRONT"
    UNKNOWN_BRANCH_FACE = "UNKNOWN_BRANCH_FACE"


class FieldContractError(ValueError):
    """Raised when a public input violates the finite SI field contract."""


class FieldNumericsError(FloatingPointError):
    """Raised when a mathematically finite query is not representable in float64."""


def _finite_scalar(name: str, value: object) -> float:
    numeric_types = (int, float, np.integer, np.floating)
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, numeric_types):
        raise FieldContractError(f"{name} must be a finite real scalar")
    result = float(value)
    if not math.isfinite(result):
        raise FieldContractError(f"{name} must be a finite real scalar")
    return result


@dataclass(frozen=True, slots=True)
class ModeIIIFixedCase:
    """Straight-front cylindrical manufactured case, in SI units."""

    radius_m: float = 1.0e-2
    length_m: float = 1.0e-2
    K_III_Pa_sqrt_m: float = 100.0
    young_modulus_Pa: float = 1.0e6
    poisson_ratio: float = 0.3

    def __post_init__(self) -> None:
        radius = _finite_scalar("radius_m", self.radius_m)
        length = _finite_scalar("length_m", self.length_m)
        k_iii = _finite_scalar("K_III_Pa_sqrt_m", self.K_III_Pa_sqrt_m)
        young = _finite_scalar("young_modulus_Pa", self.young_modulus_Pa)
        poisson = _finite_scalar("poisson_ratio", self.poisson_ratio)
        if not MIN_LENGTH_M <= radius <= MAX_LENGTH_M:
            raise FieldContractError(
                f"radius_m must be in [{MIN_LENGTH_M:g}, {MAX_LENGTH_M:g}] m"
            )
        if not MIN_LENGTH_M <= length <= MAX_LENGTH_M:
            raise FieldContractError(
                f"length_m must be in [{MIN_LENGTH_M:g}, {MAX_LENGTH_M:g}] m"
            )
        if abs(k_iii) > MAX_ABS_KIII_PA_SQRT_M:
            raise FieldContractError(
                f"abs(K_III_Pa_sqrt_m) must be <= {MAX_ABS_KIII_PA_SQRT_M:g}"
            )
        if not MIN_YOUNG_MODULUS_PA <= young <= MAX_YOUNG_MODULUS_PA:
            raise FieldContractError(
                "young_modulus_Pa must be in "
                f"[{MIN_YOUNG_MODULUS_PA:g}, {MAX_YOUNG_MODULUS_PA:g}] Pa"
            )
        if not -0.99 <= poisson <= 0.49:
            raise FieldContractError("poisson_ratio must be in [-0.99, 0.49]")
        object.__setattr__(self, "radius_m", radius)
        object.__setattr__(self, "length_m", length)
        object.__setattr__(self, "K_III_Pa_sqrt_m", k_iii)
        object.__setattr__(self, "young_modulus_Pa", young)
        object.__setattr__(self, "poisson_ratio", poisson)
        if not math.isfinite(self.shear_modulus_Pa):
            raise FieldContractError("derived shear modulus must be finite")
        if not math.isfinite(self.amplitude_sqrt_m):
            raise FieldContractError("derived Mode III amplitude must be finite")
        if not math.isfinite(self.exact_total_strain_energy_J):
            raise FieldContractError("derived exact strain energy must be finite")

    @property
    def shear_modulus_Pa(self) -> float:
        return self.young_modulus_Pa / (2.0 * (1.0 + self.poisson_ratio))

    @property
    def amplitude_sqrt_m(self) -> float:
        """A = K_III / mu * sqrt(2/pi), so u_z = A F1."""
        return self.K_III_Pa_sqrt_m / self.shear_modulus_Pa * _SQRT_TWO_OVER_PI

    @property
    def exact_total_strain_energy_J(self) -> float:
        """U = K_III^2 R L / (2 mu), integrated over the cylinder."""
        return (
            self.K_III_Pa_sqrt_m**2
            * self.radius_m
            * self.length_m
            / (2.0 * self.shear_modulus_Pa)
        )


DEFAULT_CASE = ModeIIIFixedCase()


@dataclass(frozen=True, slots=True)
class FieldResult:
    """Immutable query-major values; use a row only when its status is OK.

    f1_m_sqrt has units sqrt(m); grad_f1_per_sqrt_m has units 1/sqrt(m).
    displacement_m is a global xyz vector, displacement_gradient is
    du_i/dx_j, stress_Pa is a symmetric tensor, and optional traction_Pa is
    sigma @ n for each supplied unit normal.
    """

    statuses: np.ndarray
    f1_m_sqrt: np.ndarray
    grad_f1_per_sqrt_m: np.ndarray
    displacement_m: np.ndarray
    displacement_gradient: np.ndarray
    stress_Pa: np.ndarray
    energy_density_J_per_m3: np.ndarray
    traction_Pa: np.ndarray | None

    def status_at(self, query_index: int) -> FieldStatus:
        return FieldStatus(str(self.statuses[query_index]))


def _preflight_array(name: str, value: object) -> np.ndarray:
    if type(value) is not np.ndarray:
        raise FieldContractError(f"{name} must be an exact numpy.ndarray")
    if value.ndim != 2 or value.shape[1] != 3:
        raise FieldContractError(f"{name} must have shape (n, 3)")
    row_count = int(value.shape[0])
    if not 1 <= row_count <= MAX_QUERIES:
        raise FieldContractError(
            f"budgetRejected: {name} must have 1..{MAX_QUERIES} rows"
        )
    if value.dtype.kind not in "fiu" or value.dtype.itemsize > 8:
        raise FieldContractError(f"{name} must use a real numeric dtype of at most 8 bytes")
    if not value.flags.c_contiguous:
        raise FieldContractError(f"{name} must be C-contiguous")
    if not bool(np.isfinite(value).all()):
        raise FieldContractError(f"{name} must contain only finite values")
    if value.dtype.kind in "iu":
        exact_limit = 1 << 53
        if value.dtype.kind == "u":
            if bool(np.any(value > exact_limit)):
                raise FieldContractError(
                    f"{name} integer values must be exactly representable in float64"
                )
        elif bool(np.any(value > exact_limit)) or bool(np.any(value < -exact_limit)):
            raise FieldContractError(
                f"{name} integer values must be exactly representable in float64"
            )
    if value.nbytes > MAX_QUERIES * 3 * 8:
        raise FieldContractError(f"budgetRejected: {name} exceeds its byte bound")
    return value.astype(np.float64, copy=True)


def _validate_normals(normals: object, query_count: int) -> np.ndarray:
    result = _preflight_array("normals", normals)
    if result.shape[0] != query_count:
        raise FieldContractError("normals must contain one row per query")
    lengths = np.hypot(np.hypot(result[:, 0], result[:, 1]), result[:, 2])
    if bool(np.any(np.abs(lengths - 1.0) > NORMAL_UNIT_TOLERANCE)):
        raise FieldContractError("each supplied physical normal must have unit length")
    return result


def _immutable_bytes_array(array: np.ndarray | None) -> np.ndarray | None:
    """Back an output view with immutable bytes so WRITEABLE cannot be reenabled."""
    if array is None:
        return None
    return np.frombuffer(array.tobytes(order="C"), dtype=array.dtype).reshape(array.shape)


def evaluate_field(
    points_m: np.ndarray,
    *,
    case: ModeIIIFixedCase = DEFAULT_CASE,
    normals: np.ndarray | None = None,
) -> FieldResult:
    """Evaluate the exact field, analytic gradient, stress, and optional sigma@n.

    Query coordinates use global xyz in meters. Points outside the closed
    cylinder are OUTSIDE_DOMAIN. Exact front points (x=y=0) and exact branch
    face points (x<0, y=0) are UNKNOWN even when signed zero could select a
    one-sided atan2 value. No epsilon shift or singularity regularization is
    applied.
    """

    if type(case) is not ModeIIIFixedCase:
        raise FieldContractError("case must be a ModeIIIFixedCase")
    points = _preflight_array("points_m", points_m)
    query_count = int(points.shape[0])
    normal_rows = None if normals is None else _validate_normals(normals, query_count)

    statuses = np.full(query_count, FieldStatus.OUTSIDE_DOMAIN.value, dtype="<U24")
    f1 = np.zeros(query_count, dtype=np.float64)
    grad_f1 = np.zeros((query_count, 3), dtype=np.float64)
    displacement = np.zeros((query_count, 3), dtype=np.float64)
    displacement_gradient = np.zeros((query_count, 3, 3), dtype=np.float64)
    stress = np.zeros((query_count, 3, 3), dtype=np.float64)
    energy_density = np.zeros(query_count, dtype=np.float64)
    traction = None if normal_rows is None else np.zeros((query_count, 3), dtype=np.float64)

    radius = case.radius_m
    half_length = 0.5 * case.length_m
    x, y, z = points.T
    inside_box = (
        (np.abs(x) <= radius)
        & (np.abs(y) <= radius)
        & (np.abs(z) <= half_length)
    )
    radial_distance = np.zeros(query_count, dtype=np.float64)
    selected = np.flatnonzero(inside_box)
    radial_distance[selected] = np.hypot(x[selected], y[selected])
    inside_cylinder = inside_box & (radial_distance <= radius)
    statuses[~inside_cylinder] = FieldStatus.OUTSIDE_DOMAIN.value

    front = inside_cylinder & (x == 0.0) & (y == 0.0)
    statuses[front] = FieldStatus.UNKNOWN_FRONT.value
    branch_face = inside_cylinder & (x < 0.0) & (y == 0.0)
    statuses[branch_face] = FieldStatus.UNKNOWN_BRANCH_FACE.value
    valid = inside_cylinder & ~front & ~branch_face
    statuses[valid] = FieldStatus.OK.value

    indices = np.flatnonzero(valid)
    if indices.size:
        r = radial_distance[indices]
        with np.errstate(over="raise", divide="raise", invalid="raise"):
            try:
                theta = np.arctan2(y[indices], x[indices])
                root_r = np.sqrt(r)
                sin_half = np.sin(0.5 * theta)
                cos_half = np.cos(0.5 * theta)
                values = root_r * sin_half
                gradients = np.zeros((indices.size, 3), dtype=np.float64)
                gradients[:, 0] = -sin_half / (2.0 * root_r)
                gradients[:, 1] = cos_half / (2.0 * root_r)

                amplitude = case.amplitude_sqrt_m
                mu = case.shear_modulus_Pa
                scaled_gradients = amplitude * gradients
                f1[indices] = values
                grad_f1[indices] = gradients
                displacement[indices, 2] = amplitude * values
                displacement_gradient[indices, 2, :] = scaled_gradients
                tau_xz = mu * scaled_gradients[:, 0]
                tau_yz = mu * scaled_gradients[:, 1]
                stress[indices, 0, 2] = tau_xz
                stress[indices, 2, 0] = tau_xz
                stress[indices, 1, 2] = tau_yz
                stress[indices, 2, 1] = tau_yz
                energy_density[indices] = mu * amplitude**2 / (8.0 * r)
                if traction is not None and normal_rows is not None:
                    traction[indices] = np.einsum(
                        "nij,nj->ni", stress[indices], normal_rows[indices]
                    )
            except FloatingPointError as exc:
                raise FieldNumericsError(
                    "float64 overflow while evaluating the exact Mode III field"
                ) from exc

    arrays = (
        statuses,
        f1,
        grad_f1,
        displacement,
        displacement_gradient,
        stress,
        energy_density,
        traction,
    )
    output_bytes = sum(array.nbytes for array in arrays if array is not None)
    if output_bytes > MAX_OUTPUT_BYTES:
        raise FieldContractError(
            f"budgetRejected: result requires {output_bytes} bytes; "
            f"maximum is {MAX_OUTPUT_BYTES}"
        )
    immutable = tuple(_immutable_bytes_array(array) for array in arrays)
    return FieldResult(
        statuses=immutable[0],
        f1_m_sqrt=immutable[1],
        grad_f1_per_sqrt_m=immutable[2],
        displacement_m=immutable[3],
        displacement_gradient=immutable[4],
        stress_Pa=immutable[5],
        energy_density_J_per_m3=immutable[6],
        traction_Pa=immutable[7],
    )
