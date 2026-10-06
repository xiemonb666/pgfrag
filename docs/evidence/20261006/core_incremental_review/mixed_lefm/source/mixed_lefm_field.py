"""Source-ready straight-front mixed-mode LEFM field and branch basis.

Pure standard-library analytic reference. This is a manufactured elastic
field, not a solver, crack-growth rule, or material calibration.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import math
from typing import TypeAlias


Vector3: TypeAlias = tuple[float, float, float]
Matrix3: TypeAlias = tuple[Vector3, Vector3, Vector3]
Voigt6: TypeAlias = tuple[float, float, float, float, float, float]
MIN_LENGTH_M = 1.0e-9
MAX_LENGTH_M = 1.0e3
MIN_YOUNG_MODULUS_PA = 1.0e-3
MAX_YOUNG_MODULUS_PA = 1.0e15
MAX_ABS_K_PA_SQRT_M = 1.0e12
FRAME_TOLERANCE = 1.0e-10


class FieldStatus(str, Enum):
    OK = "OK"
    OUTSIDE_DOMAIN = "OUTSIDE_DOMAIN"
    UNKNOWN_FRONT = "UNKNOWN_FRONT"
    UNKNOWN_BRANCH_FACE = "UNKNOWN_BRANCH_FACE"


class FieldContractError(ValueError):
    """A malformed parameter, frame, point, or ambiguous requested trace."""


def _real(name: str, value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise FieldContractError(f"{name} must be a finite real scalar")
    result = float(value)
    if not math.isfinite(result):
        raise FieldContractError(f"{name} must be a finite real scalar")
    return result


def _vec3(name: str, value: object) -> Vector3:
    if type(value) not in (tuple, list) or len(value) != 3:
        raise FieldContractError(f"{name} must contain exactly three scalars")
    return tuple(_real(f"{name}[{i}]", item) for i, item in enumerate(value))  # type: ignore[return-value]


def _matrix3(name: str, value: object) -> Matrix3:
    if type(value) not in (tuple, list) or len(value) != 3:
        raise FieldContractError(f"{name} must be a 3 by 3 matrix")
    rows = tuple(_vec3(f"{name}[{i}]", row) for i, row in enumerate(value))
    return rows  # type: ignore[return-value]


def _transpose(a: Matrix3) -> Matrix3:
    return tuple(tuple(a[j][i] for j in range(3)) for i in range(3))  # type: ignore[return-value]


def _mat_vec(a: Matrix3, x: Vector3) -> Vector3:
    return tuple(sum(a[i][j] * x[j] for j in range(3)) for i in range(3))  # type: ignore[return-value]


def _mat_mul(a: Matrix3, b: Matrix3) -> Matrix3:
    bt = _transpose(b)
    return tuple(tuple(sum(a[i][k] * bt[j][k] for k in range(3)) for j in range(3)) for i in range(3))  # type: ignore[return-value]


def _det(a: Matrix3) -> float:
    return (
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
    )


def _finite_tree(name: str, values: object) -> None:
    if isinstance(values, (tuple, list)):
        for value in values:
            _finite_tree(name, value)
    elif not math.isfinite(float(values)):
        raise FloatingPointError(f"non-finite {name} from analytic field evaluation")


@dataclass(frozen=True, slots=True)
class CrackFrame:
    """Proper rotation from local crack xyz to world xyz, plus world origin.

    Matrix rows are world components and columns are local crack axes. The
    local crack occupies x < 0, y = 0; +x points toward the tip and +z follows
    the straight front, so local x cross local y equals local z.
    """

    origin_world_m: Vector3 = (0.0, 0.0, 0.0)
    axes_world_from_local: Matrix3 = (
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
    )

    def __post_init__(self) -> None:
        origin = _vec3("origin_world_m", self.origin_world_m)
        axes = _matrix3("axes_world_from_local", self.axes_world_from_local)
        gram = _mat_mul(_transpose(axes), axes)
        if any(abs(gram[i][j] - (1.0 if i == j else 0.0)) > FRAME_TOLERANCE
               for i in range(3) for j in range(3)):
            raise FieldContractError("crack frame must be orthonormal")
        if abs(_det(axes) - 1.0) > FRAME_TOLERANCE:
            raise FieldContractError("crack frame must be a proper right-handed rotation")
        object.__setattr__(self, "origin_world_m", origin)
        object.__setattr__(self, "axes_world_from_local", axes)


@dataclass(frozen=True, slots=True)
class MixedLEFMCase:
    """Finite cylindrical query domain and isotropic plane-strain LEFM data.

    Lengths use meters, E uses Pa, and each K uses Pa sqrt(m). The exact
    Williams leading field is independent of z inside the finite extrusion.
    """

    radius_m: float
    length_m: float
    K_I_Pa_sqrt_m: float
    K_II_Pa_sqrt_m: float
    K_III_Pa_sqrt_m: float
    young_modulus_Pa: float
    poisson_ratio: float
    frame: CrackFrame = CrackFrame()

    def __post_init__(self) -> None:
        for name in ("radius_m", "length_m", "K_I_Pa_sqrt_m", "K_II_Pa_sqrt_m",
                     "K_III_Pa_sqrt_m", "young_modulus_Pa", "poisson_ratio"):
            object.__setattr__(self, name, _real(name, getattr(self, name)))
        if not MIN_LENGTH_M <= self.radius_m <= MAX_LENGTH_M:
            raise FieldContractError(f"radius_m must be in [{MIN_LENGTH_M:g}, {MAX_LENGTH_M:g}] m")
        if not MIN_LENGTH_M <= self.length_m <= MAX_LENGTH_M:
            raise FieldContractError(f"length_m must be in [{MIN_LENGTH_M:g}, {MAX_LENGTH_M:g}] m")
        for name in ("K_I_Pa_sqrt_m", "K_II_Pa_sqrt_m", "K_III_Pa_sqrt_m"):
            if abs(getattr(self, name)) > MAX_ABS_K_PA_SQRT_M:
                raise FieldContractError(f"abs({name}) must be <= {MAX_ABS_K_PA_SQRT_M:g}")
        if not MIN_YOUNG_MODULUS_PA <= self.young_modulus_Pa <= MAX_YOUNG_MODULUS_PA:
            raise FieldContractError("young_modulus_Pa is outside the finite SI contract")
        if not -0.99 <= self.poisson_ratio <= 0.49:
            raise FieldContractError("poisson_ratio must be in [-0.99, 0.49]")
        if type(self.frame) is not CrackFrame:
            raise FieldContractError("frame must be a CrackFrame")
        _finite_tree("derived material constants", (self.shear_modulus_Pa, self.lame_lambda_Pa))

    @property
    def shear_modulus_Pa(self) -> float:
        return self.young_modulus_Pa / (2.0 * (1.0 + self.poisson_ratio))

    @property
    def lame_lambda_Pa(self) -> float:
        e, nu = self.young_modulus_Pa, self.poisson_ratio
        return e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))

    @property
    def plane_strain_E_prime_Pa(self) -> float:
        return self.young_modulus_Pa / (1.0 - self.poisson_ratio**2)


@dataclass(frozen=True, slots=True)
class BranchBasis:
    """Four dimensionless sqrt(r/R) functions and local gradients in 1/m."""

    values: tuple[float, float, float, float]
    gradients_per_m: tuple[Vector3, Vector3, Vector3, Vector3]
    theta_rad: float


@dataclass(frozen=True, slots=True)
class LEFMField:
    """One query. Tensor indices are displacement component, then coordinate."""

    status: FieldStatus
    theta_rad: float | None = None
    branch_trace: str | None = None
    displacement_world_m: Vector3 | None = None
    displacement_gradient_world: Matrix3 | None = None
    strain_tensor_world: Matrix3 | None = None
    strain_voigt_engineering: Voigt6 | None = None
    stress_tensor_world_Pa: Matrix3 | None = None
    stress_voigt_Pa: Voigt6 | None = None


def _trace_theta(x: float, y: float, branch_trace: str | None) -> tuple[float, str | None]:
    on_cut = x < 0.0 and y == 0.0
    if on_cut:
        if branch_trace not in ("upper", "lower"):
            raise FieldContractError("branch-face evaluation requires explicit upper or lower trace")
        return (math.pi if branch_trace == "upper" else -math.pi), branch_trace
    if branch_trace is not None:
        raise FieldContractError("branch_trace is valid only at x < 0, y = 0")
    return math.atan2(y, x), None


def branch_functions_local(
    x_m: float,
    y_m: float,
    *,
    reference_length_m: float,
    branch_trace: str | None = None,
) -> BranchBasis:
    """Evaluate the four Moes-Dolbow-Belytschko branch functions locally.

    The functions are normalized as sqrt(r/R), hence dimensionless; their
    analytic local x,y gradients have units 1/m. At r=0 they are rejected
    because the gradients are singular. Crack-face values require a side.
    """

    x, y = _real("x_m", x_m), _real("y_m", y_m)
    ref = _real("reference_length_m", reference_length_m)
    if not MIN_LENGTH_M <= ref <= MAX_LENGTH_M:
        raise FieldContractError("reference_length_m is outside the finite SI contract")
    r = math.hypot(x, y)
    if r == 0.0:
        raise FieldContractError("branch basis gradients are undefined at r = 0")
    theta, _ = _trace_theta(x, y, branch_trace)
    s, c = math.sin(theta / 2.0), math.cos(theta / 2.0)
    st, ct = math.sin(theta), math.cos(theta)
    h = (s, c, s * st, c * st)
    hp = (0.5 * c, -0.5 * s, 0.5 * c * st + s * ct,
          -0.5 * s * st + c * ct)
    q = math.sqrt(r) / math.sqrt(ref)
    scale = q / r
    values = tuple(q * value for value in h)
    gradients = tuple(
        (scale * (0.5 * ct * h[i] - st * hp[i]),
         scale * (0.5 * st * h[i] + ct * hp[i]), 0.0)
        for i in range(4)
    )
    _finite_tree("branch basis", (values, gradients))
    return BranchBasis(values, gradients, theta)  # type: ignore[arg-type]


def _add_angular_mode(
    displacement: list[float], gradient: list[list[float]], component: int,
    amplitude_per_sqrt_m: float, root_r: float, theta: float, h: float, hp: float,
) -> None:
    c, s = math.cos(theta), math.sin(theta)
    displacement[component] += amplitude_per_sqrt_m * root_r * h
    factor = amplitude_per_sqrt_m / root_r
    gradient[component][0] += factor * (0.5 * c * h - s * hp)
    gradient[component][1] += factor * (0.5 * s * h + c * hp)


def _strain_and_stress(gradient: list[list[float]], case: MixedLEFMCase) -> tuple[Matrix3, Voigt6, Matrix3, Voigt6]:
    strain = tuple(tuple(0.5 * (gradient[i][j] + gradient[j][i]) for j in range(3))
                   for i in range(3))
    mu, lam = case.shear_modulus_Pa, case.lame_lambda_Pa
    trace = sum(strain[i][i] for i in range(3))
    stress = tuple(tuple(2.0 * mu * strain[i][j] + (lam * trace if i == j else 0.0)
                         for j in range(3)) for i in range(3))
    strain_voigt = (strain[0][0], strain[1][1], strain[2][2],
                    2.0 * strain[0][1], 2.0 * strain[1][2], 2.0 * strain[0][2])
    stress_voigt = (stress[0][0], stress[1][1], stress[2][2],
                    stress[0][1], stress[1][2], stress[0][2])
    _finite_tree("strain/stress", (strain, stress, strain_voigt, stress_voigt))
    return strain, strain_voigt, stress, stress_voigt  # type: ignore[return-value]


def evaluate_local_point(local_xyz_m: object, *, case: MixedLEFMCase,
                         branch_trace: str | None = None) -> LEFMField:
    """Evaluate a local crack-frame point and return world-frame field tensors.

    The query cylinder is r <= radius_m and |z| <= length_m/2. This entry
    preserves exact local y=0 for explicit upper/lower crack-face traces.
    Exact front points return UNKNOWN_FRONT; an unqualified crack face returns
    UNKNOWN_BRANCH_FACE. No r cutoff, finite difference, or signed-zero side
    inference is used. The case frame rotates valid outputs into world xyz.
    """
    if type(case) is not MixedLEFMCase:
        raise FieldContractError("case must be a MixedLEFMCase")
    if branch_trace not in (None, "upper", "lower"):
        raise FieldContractError("branch_trace must be None, 'upper', or 'lower'")
    local = _vec3("local_xyz_m", local_xyz_m)
    frame = case.frame
    x, y, z = local
    on_cut = x < 0.0 and y == 0.0
    if branch_trace is not None and not on_cut:
        raise FieldContractError("branch_trace is valid only at x < 0, y = 0")
    r = math.hypot(x, y)
    if r > case.radius_m or abs(z) > 0.5 * case.length_m:
        return LEFMField(FieldStatus.OUTSIDE_DOMAIN)
    if r == 0.0:
        if branch_trace is not None:
            raise FieldContractError("branch_trace cannot select the crack front")
        return LEFMField(FieldStatus.UNKNOWN_FRONT)
    if x < 0.0 and y == 0.0 and branch_trace is None:
        return LEFMField(FieldStatus.UNKNOWN_BRANCH_FACE)
    theta, used_trace = _trace_theta(x, y, branch_trace)
    root_r = math.sqrt(r)
    sh, ch = math.sin(theta / 2.0), math.cos(theta / 2.0)
    st, ct = math.sin(theta), math.cos(theta)
    kappa = 3.0 - 4.0 * case.poisson_ratio
    mu = case.shear_modulus_Pa
    displacement = [0.0, 0.0, 0.0]
    gradient = [[0.0, 0.0, 0.0] for _ in range(3)]
    norm = 1.0 / math.sqrt(2.0 * math.pi)
    a_i = case.K_I_Pa_sqrt_m * norm / (2.0 * mu)
    a_ii = case.K_II_Pa_sqrt_m * norm / (2.0 * mu)
    _add_angular_mode(displacement, gradient, 0, a_i, root_r, theta,
                      ch * (kappa - ct), -0.5 * sh * (kappa - ct) + ch * st)
    _add_angular_mode(displacement, gradient, 1, a_i, root_r, theta,
                      sh * (kappa - ct), 0.5 * ch * (kappa - ct) + sh * st)
    _add_angular_mode(displacement, gradient, 0, a_ii, root_r, theta,
                      sh * (kappa + 2.0 + ct), 0.5 * ch * (kappa + 2.0 + ct) - sh * st)
    _add_angular_mode(displacement, gradient, 1, a_ii, root_r, theta,
                      -ch * (kappa - 2.0 + ct),
                      0.5 * sh * (kappa - 2.0 + ct) + ch * st)
    a_iii = case.K_III_Pa_sqrt_m / mu * math.sqrt(2.0 / math.pi)
    _add_angular_mode(displacement, gradient, 2, a_iii, root_r, theta, sh, 0.5 * ch)
    _finite_tree("local displacement/gradient", (displacement, gradient))
    strain, strain_voigt, stress, stress_voigt = _strain_and_stress(gradient, case)
    q = frame.axes_world_from_local
    u_world = _mat_vec(q, tuple(displacement))
    grad_world = _mat_mul(_mat_mul(q, tuple(tuple(row) for row in gradient)), _transpose(q))
    strain_world = _mat_mul(_mat_mul(q, strain), _transpose(q))
    stress_world = _mat_mul(_mat_mul(q, stress), _transpose(q))
    strain_voigt_world = (strain_world[0][0], strain_world[1][1], strain_world[2][2],
                          2.0 * strain_world[0][1], 2.0 * strain_world[1][2],
                          2.0 * strain_world[0][2])
    stress_voigt_world = (stress_world[0][0], stress_world[1][1], stress_world[2][2],
                          stress_world[0][1], stress_world[1][2], stress_world[0][2])
    _finite_tree("world field", (u_world, grad_world, strain_world, stress_world,
                                 strain_voigt_world, stress_voigt_world))
    return LEFMField(
        status=FieldStatus.OK, theta_rad=theta, branch_trace=used_trace,
        displacement_world_m=u_world, displacement_gradient_world=grad_world,
        strain_tensor_world=strain_world, strain_voigt_engineering=strain_voigt_world,
        stress_tensor_world_Pa=stress_world, stress_voigt_Pa=stress_voigt_world,
    )


def evaluate_point(point_world_m: object, *, case: MixedLEFMCase,
                   branch_trace: str | None = None) -> LEFMField:
    """Transform a world point to local coordinates and evaluate it fail-closed.

    For exact crack-face traces in a rotated frame, call evaluate_local_point
    with an exactly specified local y=0 point instead of relying on a world to
    local floating-point round trip.
    """
    if type(case) is not MixedLEFMCase:
        raise FieldContractError("case must be a MixedLEFMCase")
    if branch_trace not in (None, "upper", "lower"):
        raise FieldContractError("branch_trace must be None, 'upper', or 'lower'")
    point = _vec3("point_world_m", point_world_m)
    frame = case.frame
    delta = tuple(point[i] - frame.origin_world_m[i] for i in range(3))
    local = _mat_vec(_transpose(frame.axes_world_from_local), delta)
    _finite_tree("local coordinates", local)
    return evaluate_local_point(local, case=case, branch_trace=branch_trace)


def energy_release_rate_J_per_m2(case: MixedLEFMCase) -> float:
    """Plane-strain isotropic G=(KI²+KII²)/E' + KIII²/(2 mu).

    This is an elastic LEFM benchmark identity, not a toughness calibration
    or autonomous fracture criterion.
    """

    if type(case) is not MixedLEFMCase:
        raise FieldContractError("case must be a MixedLEFMCase")
    value = ((case.K_I_Pa_sqrt_m**2 + case.K_II_Pa_sqrt_m**2)
             / case.plane_strain_E_prime_Pa
             + case.K_III_Pa_sqrt_m**2 / (2.0 * case.shear_modulus_Pa))
    _finite_tree("energy release rate", value)
    return value
