"""Bounded small-strain 3-D elastic EFGM Galerkin operator.

This source-only primitive assembles

    K = sum_q B(q).T @ D @ B(q) * volume_weight(q)

from a caller-supplied shape provider. It does not implement crack visibility,
fracture growth, loading, a linear solver, or a benchmark. In particular, the
provider owns the MLS kernel and must return physical shape values as well as
analytic gradients; this module never treats MLS coefficients as physical
Dirichlet values.

Engineering-shear order is [eps_xx, eps_yy, eps_zz, gamma_xy, gamma_yz,
gamma_xz], where gamma_ij = du_i/dX_j + du_j/dX_i. The conjugate stress order
is [sig_xx, sig_yy, sig_zz, tau_xy, tau_yz, tau_xz]. Thus each shear diagonal
in D is mu, since tau_ij = mu * gamma_ij.

All coordinates and quadrature volume weights use one consistent length unit.
E and stress use one consistent pressure unit. K has units pressure*length
for displacement coefficients measured in the coordinate length unit; its
quadratic form has units of energy. The dimensionless fixture in the companion
contract uses E=1, side length 2, and sigma=1e-3.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Protocol

import numpy as np


MAX_DOF = 375
MAX_NODES = MAX_DOF // 3
MAX_QUADRATURE_POINTS = 256


class ShapeStatus(str, Enum):
    """Status returned by an injected shape provider at one physical point."""

    OK = "OK"
    UNKNOWN_POOR_SUPPORT = "UNKNOWN_POOR_SUPPORT"
    UNKNOWN = "UNKNOWN"
    BLOCKED = "BLOCKED"


class AssemblyStatus(str, Enum):
    OK = "OK"
    BLOCKED_UNKNOWN_SUPPORT = "BLOCKED_UNKNOWN_SUPPORT"
    BLOCKED_SHAPE_STATUS = "BLOCKED_SHAPE_STATUS"


@dataclass(frozen=True)
class ShapeSample:
    """One physical field evaluation expressed in global coefficient IDs.

    ``weights[k]`` and ``gradients[k, :]`` belong to ``support_ids[k]``.
    ``gradients[k, j]`` is d(N_k)/dX_j in the material/reference frame. The
    provider returns dense all-node rows in canonical order, so support_ids is
    exactly 0..N-1 and hidden/out-of-support contributions are zero rows.
    For a non-OK status the arrays may be empty; assembly blocks before using
    them. The provider must return exact plain ndarrays, not array subclasses.
    """

    weights: object
    gradients: object
    support_ids: object
    status: ShapeStatus
    diagnostic: str | None = None


class ShapeProvider(Protocol):
    """Adapter contract for MLS or another partition-of-unity shape source."""

    def evaluate(self, material_position: np.ndarray) -> ShapeSample:
        """Return values, analytic gradients, global support IDs, and status."""


class ShapeProviderContractError(ValueError):
    """Raised when an OK shape sample violates the bounded provider contract."""


class ElasticOperatorNumericsError(FloatingPointError):
    """Raised when finite input data overflows a constitutive or K value."""


def _provider_ndarray(name: str, value: object) -> np.ndarray:
    try:
        return _plain_ndarray(name, value)
    except ValueError as exc:
        raise ShapeProviderContractError(str(exc)) from exc


@dataclass(frozen=True)
class IsotropicMaterial:
    """Stable isotropic linear elastic material in a consistent unit system."""

    young_modulus: float
    poisson_ratio: float

    def __post_init__(self) -> None:
        e = self.young_modulus
        nu = self.poisson_ratio
        numeric_types = (int, float, np.integer, np.floating)
        if isinstance(e, (bool, np.bool_)) or not isinstance(e, numeric_types):
            raise ValueError("young_modulus must be a real scalar")
        if isinstance(nu, (bool, np.bool_)) or not isinstance(nu, numeric_types):
            raise ValueError("poisson_ratio must be a real scalar")
        e = float(e)
        nu = float(nu)
        if not np.isfinite(e) or e <= 0.0:
            raise ValueError("young_modulus must be finite and positive")
        if not np.isfinite(nu) or not -1.0 < nu < 0.5:
            raise ValueError("poisson_ratio must satisfy -1 < nu < 0.5")
        lam = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
        mu = e / (2.0 * (1.0 + nu))
        if not np.isfinite(lam) or not np.isfinite(mu) or mu <= 0.0:
            raise ValueError("derived Lamé parameters must be finite and stable")
        object.__setattr__(self, "young_modulus", e)
        object.__setattr__(self, "poisson_ratio", nu)

    @property
    def lame_lambda(self) -> float:
        e, nu = self.young_modulus, self.poisson_ratio
        return e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))

    @property
    def shear_modulus(self) -> float:
        return self.young_modulus / (2.0 * (1.0 + self.poisson_ratio))

    @property
    def constitutive_matrix(self) -> np.ndarray:
        """Return 6x6 D for the documented engineering-shear convention."""
        lam = self.lame_lambda
        mu = self.shear_modulus
        d = np.zeros((6, 6), dtype=np.float64)
        with np.errstate(over="ignore", invalid="ignore"):
            d[:3, :3] = lam
            d[0, 0] += 2.0 * mu
            d[1, 1] += 2.0 * mu
            d[2, 2] += 2.0 * mu
            d[3, 3] = mu
            d[4, 4] = mu
            d[5, 5] = mu
        if not np.all(np.isfinite(d)):
            raise ElasticOperatorNumericsError("constitutive matrix D overflowed")
        d.setflags(write=False)
        return d


@dataclass(frozen=True)
class ElasticAssemblyResult:
    """Complete operator or a fail-closed support status with no partial K."""

    status: AssemblyStatus
    stiffness: np.ndarray | None
    dof_count: int
    quadrature_count: int
    blocked_quadrature_index: int | None = None
    shape_status: ShapeStatus | None = None
    diagnostic: str | None = None


def _plain_ndarray(name: str, value: object) -> np.ndarray:
    if type(value) is not np.ndarray:
        raise ValueError(f"{name} must be a plain numpy.ndarray")
    return value


def _check_real_finite(name: str, a: np.ndarray) -> None:
    if a.dtype.kind not in "fiu":
        raise ValueError(f"{name} must use a real numeric dtype")
    if not np.all(np.isfinite(a)):
        raise ValueError(f"{name} must contain only finite values")


def _validate_coordinates(name: str, value: object, max_rows: int) -> tuple[np.ndarray, int]:
    a = _plain_ndarray(name, value)
    if a.ndim != 2 or a.shape[1] != 3:
        raise ValueError(f"{name} must have shape (n, 3)")
    count = int(a.shape[0])
    if count == 0:
        raise ValueError(f"{name} must contain at least one point")
    if count > max_rows:
        raise ValueError(f"budgetRejected: {name} has {count} rows; maximum is {max_rows}")
    _check_real_finite(name, a)
    return a, count


def _validated_sample(sample: object, node_count: int) -> ShapeSample:
    if type(sample) is not ShapeSample:
        raise ShapeProviderContractError("shape provider must return ShapeSample")
    if type(sample.status) is not ShapeStatus:
        raise ShapeProviderContractError("ShapeSample.status must be a ShapeStatus")
    if sample.diagnostic is not None and type(sample.diagnostic) is not str:
        raise ShapeProviderContractError("ShapeSample.diagnostic must be a string or None")
    if sample.status is not ShapeStatus.OK:
        if sample.status is ShapeStatus.UNKNOWN and not sample.diagnostic:
            raise ShapeProviderContractError(
                "generic UNKNOWN status must preserve its provider reason in diagnostic"
            )
        return sample

    weights = _provider_ndarray("shape weights", sample.weights)
    gradients = _provider_ndarray("shape gradients", sample.gradients)
    support_ids = _provider_ndarray("shape support_ids", sample.support_ids)
    if weights.ndim != 1:
        raise ShapeProviderContractError("shape weights must have shape (n_support,)")
    support_count = int(weights.shape[0])
    if support_count != node_count:
        raise ShapeProviderContractError(
            f"dense all-node shape rows must match node_count ({node_count})"
        )
    if gradients.shape != (support_count, 3):
        raise ShapeProviderContractError("shape gradients must have shape (n_support, 3)")
    if support_ids.shape != (support_count,) or support_ids.dtype.kind not in "iu":
        raise ShapeProviderContractError("shape support_ids must be an integer vector matching weights")
    try:
        _check_real_finite("shape weights", weights)
        _check_real_finite("shape gradients", gradients)
    except ValueError as exc:
        raise ShapeProviderContractError(str(exc)) from exc
    expected_ids = np.arange(node_count, dtype=support_ids.dtype)
    if not np.array_equal(support_ids, expected_ids):
        raise ShapeProviderContractError(
            "dense shape rows require canonical global support_ids 0..node_count-1"
        )

    # Copy only after checking every shape and support-size cap.
    copied_weights = np.array(weights, dtype=np.float64, copy=True)
    copied_gradients = np.array(gradients, dtype=np.float64, copy=True)
    copied_ids = np.array(support_ids, dtype=np.int64, copy=True)
    for a in (copied_weights, copied_gradients, copied_ids):
        a.setflags(write=False)
    return ShapeSample(copied_weights, copied_gradients, copied_ids, ShapeStatus.OK, sample.diagnostic)


def _strain_displacement_matrix(gradients: np.ndarray) -> np.ndarray:
    """Build B with engineering shear rows and node-major vector DOFs."""
    support_count = gradients.shape[0]
    b = np.zeros((6, 3 * support_count), dtype=np.float64)
    for local_id, (gx, gy, gz) in enumerate(gradients):
        col = 3 * local_id
        b[0, col] = gx
        b[1, col + 1] = gy
        b[2, col + 2] = gz
        b[3, col] = gy
        b[3, col + 1] = gx
        b[4, col + 1] = gz
        b[4, col + 2] = gy
        b[5, col] = gz
        b[5, col + 2] = gx
    return b


def assemble_elastic_operator(
    coefficient_positions: np.ndarray,
    quadrature_positions: np.ndarray,
    quadrature_volume_weights: np.ndarray,
    material: IsotropicMaterial,
    shape_provider: ShapeProvider,
) -> ElasticAssemblyResult:
    """Assemble bounded dense K; no loads, constraints, solve, or crack edits.

    Exact ndarray types, all input dimensions, the 375-DOF cap, and the
    256-point quadrature cap are checked before coordinate copies or K
    allocation. A non-OK provider status returns a blocked result with
    ``stiffness=None``. Invalid provider arrays raise
    ShapeProviderContractError. Non-finite D, local stiffness, or accumulated
    global stiffness raises ElasticOperatorNumericsError. No failure returns
    a partial matrix or replaces overflow with a clipped value.
    """
    nodes_input, node_count = _validate_coordinates(
        "coefficient_positions", coefficient_positions, MAX_NODES
    )
    dof_count = 3 * node_count
    if dof_count > MAX_DOF:
        raise ValueError(f"budgetRejected: DOF count {dof_count} exceeds {MAX_DOF}")
    qpoints_input, quadrature_count = _validate_coordinates(
        "quadrature_positions", quadrature_positions, MAX_QUADRATURE_POINTS
    )
    volume_input = _plain_ndarray("quadrature_volume_weights", quadrature_volume_weights)
    if volume_input.ndim != 1 or volume_input.shape[0] != quadrature_count:
        raise ValueError("quadrature_volume_weights must have one entry per quadrature position")
    _check_real_finite("quadrature_volume_weights", volume_input)
    if np.any(volume_input <= 0.0):
        raise ValueError("quadrature_volume_weights must be positive")
    if type(material) is not IsotropicMaterial:
        raise ValueError("material must be an IsotropicMaterial")
    evaluate = getattr(shape_provider, "evaluate", None)
    if not callable(evaluate):
        raise ValueError("shape_provider must implement evaluate(material_position)")

    # All public input caps are now known and accepted; only bounded copies follow.
    qpoints = np.array(qpoints_input, dtype=np.float64, copy=True)
    volumes = np.array(volume_input, dtype=np.float64, copy=True)
    qpoints.setflags(write=False)

    # Validate and freeze all shape rows before allocating the global matrix.
    samples: list[ShapeSample] = []
    for qid in range(quadrature_count):
        sample = _validated_sample(evaluate(qpoints[qid]), node_count)
        if sample.status is not ShapeStatus.OK:
            if sample.status is ShapeStatus.UNKNOWN_POOR_SUPPORT or sample.status is ShapeStatus.UNKNOWN:
                assembly_status = AssemblyStatus.BLOCKED_UNKNOWN_SUPPORT
            else:
                assembly_status = AssemblyStatus.BLOCKED_SHAPE_STATUS
            return ElasticAssemblyResult(
                status=assembly_status,
                stiffness=None,
                dof_count=dof_count,
                quadrature_count=quadrature_count,
                blocked_quadrature_index=qid,
                shape_status=sample.status,
                diagnostic=sample.diagnostic,
            )
        samples.append(sample)

    constitutive = material.constitutive_matrix
    stiffness = np.zeros((dof_count, dof_count), dtype=np.float64)
    dof_offsets = np.arange(3, dtype=np.int64)
    for qid, sample in enumerate(samples):
        gradients = sample.gradients
        support_ids = sample.support_ids
        b = _strain_displacement_matrix(gradients)
        with np.errstate(over="ignore", invalid="ignore"):
            local_stiffness = (b.T @ constitutive @ b) * volumes[qid]
        if not np.all(np.isfinite(local_stiffness)):
            raise ElasticOperatorNumericsError(
                f"local stiffness overflowed at quadrature index {qid}"
            )
        local_dofs = (3 * support_ids[:, None] + dof_offsets[None, :]).reshape(-1)
        with np.errstate(over="ignore", invalid="ignore"):
            stiffness[np.ix_(local_dofs, local_dofs)] += local_stiffness
        if not np.all(np.isfinite(stiffness)):
            raise ElasticOperatorNumericsError(
                f"global stiffness overflowed after quadrature index {qid}"
            )

    stiffness.setflags(write=False)
    return ElasticAssemblyResult(
        status=AssemblyStatus.OK,
        stiffness=stiffness,
        dof_count=dof_count,
        quadrature_count=quadrature_count,
    )


def physical_displacement_row(sample: ShapeSample, node_count: int) -> np.ndarray:
    """Return H(X), shape (3, 3*node_count), mapping coefficients to u(X).

    Use ``H @ d = prescribed_u`` to represent a physical Dirichlet condition.
    This function deliberately does not clamp coefficient entries: MLS
    coefficients are not generally physical displacements at their own sites.
    """
    if not isinstance(node_count, int) or isinstance(node_count, bool) or not 1 <= node_count <= MAX_NODES:
        raise ValueError(f"node_count must be in [1, {MAX_NODES}]")
    validated = _validated_sample(sample, node_count)
    if validated.status is not ShapeStatus.OK:
        raise ValueError(f"physical displacement row is unavailable for status {validated.status.value}")
    row = np.zeros((3, 3 * node_count), dtype=np.float64)
    for value, support_id in zip(validated.weights, validated.support_ids):
        offset = 3 * int(support_id)
        row[0, offset] = value
        row[1, offset + 1] = value
        row[2, offset + 2] = value
    row.setflags(write=False)
    return row


def rigid_mode_coefficient_vectors(coefficient_positions: np.ndarray) -> np.ndarray:
    """Return the six coefficient vectors for translations and rotations.

    The result has shape (3*n_nodes, 6), with node-major vector DOFs. Its
    columns represent three unit translations and three unit infinitesimal
    rotations about the centroid of the supplied material coordinates. For a
    linearly complete shape provider, K times these vectors should vanish up to
    quadrature and floating-point error. These are diagnostic mode vectors,
    not the integral gauge constraints required by a static solve.
    """
    nodes_input, node_count = _validate_coordinates(
        "coefficient_positions", coefficient_positions, MAX_NODES
    )
    dof_count = 3 * node_count
    if dof_count > MAX_DOF:
        raise ValueError(f"budgetRejected: DOF count {dof_count} exceeds {MAX_DOF}")
    nodes = np.array(nodes_input, dtype=np.float64, copy=True)
    centered = nodes - np.mean(nodes, axis=0)
    modes = np.zeros((dof_count, 6), dtype=np.float64)
    axes = np.eye(3, dtype=np.float64)
    for node_id, position in enumerate(centered):
        block = slice(3 * node_id, 3 * node_id + 3)
        modes[block, :3] = axes
        for rotation_axis in range(3):
            modes[block, 3 + rotation_axis] = np.cross(axes[rotation_axis], position)
    modes.setflags(write=False)
    return modes

