"""Mass-normalized stiffness bounds for small-strain T4/Tri3 cohesive meshes.

The bound is assembled element by element, so memory use is linear in the number
of bulk and cohesive elements.  It applies to the initial linear bulk operator
and the declared cohesive material tangent for the supplied rest-state face
normals.  The consistent Tri3 area matrix is a positive-semidefinite upper
envelope over a native one-point ``A/9 * ones`` facet quadrature; it is not an
assertion about the native cohesive quadrature itself.

This is not a full-trajectory guarantee.  Current-normal geometric derivatives,
finite-deformation terms, contact, and the evolving softening tangent are out of
scope.  A negative softening eigenvalue is a physical instability; the returned
central-difference limit addresses high-frequency numerical instability of the
positive initial/fixed-normal stiffness operator only.  A practical step should
be at most half the returned critical step, and still needs separate trajectory
and energy checks.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


_DENSE_VALIDATION_MAX_DOFS = 2_000
_NATIVE_MASS_RELATIVE_TOLERANCE = 1.0e-8


@dataclass(frozen=True, slots=True)
class CohesiveStabilityBound:
    """Result and audit quantities for a row-sum stability upper bound.

    ``lambda_upper_bound_per_s2`` bounds the largest eigenvalue of the free-DOF
    mass-normalized initial stiffness matrix by Gershgorin's theorem.  The
    optional dense eigenvalue is a small-problem diagnostic, not used to form
    the bound.  Per-node/per-DOF arrays are returned as read-only float/bool
    arrays so callers can serialize or inspect the calculation safely.
    """

    lambda_upper_bound_per_s2: float
    critical_time_step_s: float
    recommended_time_step_s: float
    bulk_lumped_node_mass_kg: np.ndarray
    mass_used_node_kg: np.ndarray
    bulk_total_lumped_mass_kg: float
    cohesive_face_areas_m2: np.ndarray
    cohesive_total_area_m2: float
    bulk_row_sum_contributions_per_s2: np.ndarray
    cohesive_row_sum_contributions_per_s2: np.ndarray
    free_dof_count: int
    blocked_dof_count: int
    total_dof_count: int
    mass_source: str
    dense_max_eigenvalue_per_s2: float | None
    dense_min_eigenvalue_per_s2: float | None
    cohesive_area_matrix_policy: str = (
        "CONSISTENT_TRI3_PSD_UPPER_ENVELOPE_OVER_NATIVE_SINGLE_POINT_A_OVER_9_ONES"
    )
    stability_scope: str = (
        "INITIAL_LINEAR_BULK_PLUS_DECLARED_COHESIVE_TANGENT_AT_REST_FACE_NORMALS; "
        "NO_CURRENT_NORMAL_GEOMETRIC_DERIVATIVE; NO_FULL_SOFTENING_OR_CONTACT_CLAIM"
    )
    interpretation: str = (
        "CENTRAL_DIFFERENCE_HIGH_FREQUENCY_LIMIT_FOR_POSITIVE_STIFFNESS; "
        "NEGATIVE_SOFTENING_EIGENVALUES_ARE_PHYSICAL_INSTABILITY"
    )


def _readonly(array: np.ndarray) -> np.ndarray:
    result = np.asarray(array).copy()
    result.setflags(write=False)
    return result


def _finite_scalar(name: str, value: float, *, positive: bool = False) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a finite scalar") from exc
    if not math.isfinite(result) or (positive and result <= 0.0):
        qualifier = "finite and positive" if positive else "finite"
        raise ValueError(f"{name} must be {qualifier}")
    return result


def _connectivity(name: str, values: np.ndarray, shape_tail: tuple[int, ...], node_count: int) -> np.ndarray:
    raw = np.asarray(values)
    if raw.ndim != 1 + len(shape_tail) or raw.shape[1:] != shape_tail:
        raise ValueError(f"{name} must have shape (count, {', '.join(map(str, shape_tail))})")
    if raw.dtype.kind not in "iu" or raw.dtype.kind == "b":
        raise ValueError(f"{name} must contain integer node indices")
    if raw.size == 0 and raw.shape[0] != 0:
        raise ValueError(f"{name} has an invalid empty shape")
    if raw.size and (np.any(raw < 0) or np.any(raw >= node_count)):
        raise ValueError(f"{name} contains a node index outside [0, {node_count})")
    result = raw.astype(np.int64, copy=False)
    if result.shape[0]:
        if np.any(np.sort(result, axis=-1)[..., 1:] == np.sort(result, axis=-1)[..., :-1]):
            raise ValueError(f"each {name} element must contain distinct nodes")
    return result


def _blocked_mask(values: np.ndarray, node_count: int) -> np.ndarray:
    raw = np.asarray(values)
    if raw.shape != (node_count, 3):
        raise ValueError(f"blocked_dofs must have shape ({node_count}, 3)")
    if raw.dtype.kind == "b":
        return raw.astype(bool, copy=True)
    if raw.dtype.kind not in "iuf" or not np.all(np.isfinite(raw)) or np.any((raw != 0) & (raw != 1)):
        raise ValueError("blocked_dofs must contain only booleans or exact 0/1 values")
    return raw.astype(bool)


def _tetrahedron_stiffness(
    coordinates_m: np.ndarray,
    young_modulus_pa: float,
    poisson_ratio: float,
) -> tuple[np.ndarray, float]:
    """Return the small-strain isotropic T4 stiffness and positive volume."""
    x0 = coordinates_m[0]
    jacobian = np.column_stack((coordinates_m[1] - x0, coordinates_m[2] - x0, coordinates_m[3] - x0))
    determinant = float(np.linalg.det(jacobian))
    if not math.isfinite(determinant) or determinant == 0.0:
        raise ValueError("tetrahedra must have finite, nondegenerate geometry")
    volume = abs(determinant) / 6.0
    if not math.isfinite(volume) or volume <= 0.0:
        raise ValueError("tetrahedra must have finite, positive volume")

    # xi = inverse(J) @ (x-x0), so each physical shape gradient is a ROW
    # of inverse(J). A transpose would fail on skew tetrahedra.
    inverse_jacobian = np.linalg.inv(jacobian)
    gradients = np.empty((4, 3), dtype=np.float64)
    gradients[1:] = inverse_jacobian
    gradients[0] = -gradients[1:].sum(axis=0)

    # Engineering strain ordering: [xx, yy, zz, yz, xz, xy].
    b = np.zeros((6, 12), dtype=np.float64)
    for node, (gx, gy, gz) in enumerate(gradients):
        column = 3 * node
        b[0, column] = gx
        b[1, column + 1] = gy
        b[2, column + 2] = gz
        b[3, column + 1] = gz
        b[3, column + 2] = gy
        b[4, column] = gz
        b[4, column + 2] = gx
        b[5, column] = gy
        b[5, column + 1] = gx

    lame_lambda = young_modulus_pa * poisson_ratio / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio))
    shear_mu = young_modulus_pa / (2.0 * (1.0 + poisson_ratio))
    constitutive = np.zeros((6, 6), dtype=np.float64)
    constitutive[:3, :3] = lame_lambda
    np.fill_diagonal(constitutive[:3, :3], lame_lambda + 2.0 * shear_mu)
    constitutive[3:, 3:] = np.eye(3) * shear_mu
    stiffness = volume * (b.T @ constitutive @ b)
    if not np.all(np.isfinite(stiffness)):
        raise ValueError("tetrahedron stiffness overflowed; check geometry and material units")
    # Symmetrize only roundoff from the matrix products; this does not regularize
    # or alter any zero/negative mass or geometric singularity.
    stiffness = 0.5 * (stiffness + stiffness.T)
    return stiffness, volume


def tri3_cohesive_stiffness_matrix(
    face_coordinates_m: np.ndarray,
    normal_stiffness_Pa_per_m: float,
    shear_stiffness_Pa_per_m: float,
) -> tuple[np.ndarray, float, np.ndarray]:
    """Build one consistent-area Tri3 cohesive stiffness matrix.

    ``face_coordinates_m`` is a (2, 3, 3) array; the two sides must have exact
    matching rest coordinates in source-ancestor node order.  DOF order in the
    returned 18x18 matrix is side 0's three xyz nodes, then side 1's three xyz
    nodes.  This helper makes the energy convention auditable and is also used
    by the global bound assembler.
    """
    coordinates = np.asarray(face_coordinates_m, dtype=np.float64)
    if coordinates.shape != (2, 3, 3) or not np.all(np.isfinite(coordinates)):
        raise ValueError("face_coordinates_m must be a finite array with shape (2, 3, 3)")
    if not np.array_equal(coordinates[0], coordinates[1]):
        raise ValueError("paired cohesive face coordinates must match exactly in ancestor order")
    kn = _finite_scalar("normal_stiffness_Pa_per_m", normal_stiffness_Pa_per_m, positive=True)
    kt = _finite_scalar("shear_stiffness_Pa_per_m", shear_stiffness_Pa_per_m, positive=True)

    cross = np.cross(coordinates[0, 1] - coordinates[0, 0], coordinates[0, 2] - coordinates[0, 0])
    double_area = float(np.linalg.norm(cross))
    if not math.isfinite(double_area) or double_area <= 0.0:
        raise ValueError("cohesive Tri3 faces must have finite, positive area")
    area = 0.5 * double_area
    normal = cross / double_area
    tangent = np.eye(3) - np.outer(normal, normal)
    traction_stiffness = kn * np.outer(normal, normal) + kt * tangent
    area_mass = (area / 12.0) * np.array(
        [[2.0, 1.0, 1.0], [1.0, 2.0, 1.0], [1.0, 1.0, 2.0]],
        dtype=np.float64,
    )
    signed_area_mass = np.block([[area_mass, -area_mass], [-area_mass, area_mass]])
    stiffness = np.kron(signed_area_mass, traction_stiffness)
    if not np.all(np.isfinite(stiffness)):
        raise ValueError("cohesive stiffness overflowed; check face geometry and tangent units")
    stiffness = 0.5 * (stiffness + stiffness.T)
    return stiffness, area, normal


def _native_mass_vector(native_point_mass_kg: np.ndarray, node_count: int) -> np.ndarray:
    raw = np.asarray(native_point_mass_kg, dtype=np.float64)
    if raw.shape == (node_count, 3):
        if not np.all(np.isfinite(raw)):
            raise ValueError("native_point_mass_kg must contain only finite values")
        scale = np.max(np.abs(raw), axis=1)
        if np.any(np.max(np.abs(raw - raw[:, :1]), axis=1) > _NATIVE_MASS_RELATIVE_TOLERANCE * scale):
            raise ValueError("native_point_mass_kg xyz components must carry the same point mass")
        raw = raw[:, 0]
    elif raw.shape != (node_count,):
        raise ValueError(f"native_point_mass_kg must have shape ({node_count},) or ({node_count}, 3)")
    if not np.all(np.isfinite(raw)) or np.any(raw < 0.0):
        raise ValueError("native_point_mass_kg must be finite and nonnegative")
    return np.asarray(raw, dtype=np.float64).copy()


def mass_normalized_stiffness_bound(
    nodes_m: np.ndarray,
    tets: np.ndarray,
    paired_faces: np.ndarray,
    rho_kg_m3: float,
    young_modulus_Pa: float,
    poisson_ratio: float,
    normal_stiffness_Pa_per_m: float,
    shear_stiffness_Pa_per_m: float,
    blocked_dofs: np.ndarray,
    native_point_mass_kg: np.ndarray | None = None,
    validate_dense_eigenvalue: bool = False,
) -> CohesiveStabilityBound:
    """Return a conservative free-DOF mass-normalized stiffness row-sum bound.

    Parameters
    ----------
    nodes_m:
        Rest coordinates in metres, shape (N, 3).
    tets:
        P1 tetrahedron node indices, shape (T, 4).
    paired_faces:
        Source-ancestor-matched Tri3 pairs, shape (F, 2, 3).  Axis 1 selects
        the two sides and the final axis preserves corresponding vertex order.
        Their rest coordinates must match exactly.
    blocked_dofs:
        Boolean (N, 3) mask where True means prescribed/constrained.
    native_point_mass_kg:
        Optional native point masses, either (N,) or (N, 3).  For a 3-vector
        per node, xyz components are checked for equality and reduced to one
        scalar mass; the vector is never summed across xyz.  Supplied values
        must match the T4-derived lumped bulk masses within relative 1e-8.
    validate_dense_eigenvalue:
        If True, assemble a dense diagnostic matrix and return its minimum and
        maximum eigenvalues.  Limited to at most 2,000 total DOFs.  This does
        not participate in the scalable row-sum bound.

    The local absolute-row contributions are accumulated before global
    cancellation.  Therefore they upper-bound the absolute row sums of the
    globally assembled matrix by the triangle inequality, including cases
    where adjacent element stiffness entries cancel.
    """
    nodes = np.asarray(nodes_m, dtype=np.float64)
    if nodes.ndim != 2 or nodes.shape[1] != 3 or nodes.shape[0] < 4:
        raise ValueError("nodes_m must have shape (N, 3) with at least four nodes")
    if not np.all(np.isfinite(nodes)):
        raise ValueError("nodes_m must contain only finite coordinates")
    node_count = nodes.shape[0]
    elements = _connectivity("tets", tets, (4,), node_count)
    if elements.shape[0] == 0:
        raise ValueError("at least one T4 element is required to define bulk point masses")
    faces = np.asarray(paired_faces)
    if faces.ndim != 3 or faces.shape[1:] != (2, 3):
        raise ValueError("paired_faces must have shape (F, 2, 3)")
    face_connectivity = _connectivity("paired_faces", faces.reshape((-1, 3)), (3,), node_count).reshape((-1, 2, 3))
    for pair in face_connectivity:
        if np.intersect1d(pair[0], pair[1]).size:
            raise ValueError("the two sides of each cohesive face must use distinct node indices")
        if not np.array_equal(nodes[pair[0]], nodes[pair[1]]):
            raise ValueError("paired cohesive face coordinates must match exactly in ancestor order")
    mask = _blocked_mask(blocked_dofs, node_count)

    rho = _finite_scalar("rho_kg_m3", rho_kg_m3, positive=True)
    young = _finite_scalar("young_modulus_Pa", young_modulus_Pa, positive=True)
    nu = _finite_scalar("poisson_ratio", poisson_ratio)
    if not (-1.0 < nu < 0.5):
        raise ValueError("poisson_ratio must be strictly between -1 and 0.5 for stable isotropic elasticity")
    kn = _finite_scalar("normal_stiffness_Pa_per_m", normal_stiffness_Pa_per_m, positive=True)
    kt = _finite_scalar("shear_stiffness_Pa_per_m", shear_stiffness_Pa_per_m, positive=True)
    validate_dense = bool(validate_dense_eigenvalue)
    total_dof_count = 3 * node_count
    if validate_dense and total_dof_count > _DENSE_VALIDATION_MAX_DOFS:
        raise ValueError(
            f"dense eigenvalue validation is limited to {_DENSE_VALIDATION_MAX_DOFS} total DOFs"
        )

    bulk_mass = np.zeros(node_count, dtype=np.float64)
    volumes = np.empty(elements.shape[0], dtype=np.float64)
    # Validate all tetrahedra and derive the lumped mass once before assembly.
    tetra_stiffnesses: list[np.ndarray] = []
    for index, nodes_in_tet in enumerate(elements):
        stiffness, volume = _tetrahedron_stiffness(nodes[nodes_in_tet], young, nu)
        tetra_stiffnesses.append(stiffness)
        volumes[index] = volume
        bulk_mass[nodes_in_tet] += (rho * volume) / 4.0
    if not np.all(np.isfinite(bulk_mass)):
        raise ValueError("bulk lumped mass overflowed; check geometry and material units")

    mass_source = "T4_DERIVED_LUMPED_BULK_MASS"
    if native_point_mass_kg is None:
        mass = bulk_mass.copy()
    else:
        native_mass = _native_mass_vector(native_point_mass_kg, node_count)
        scale = np.maximum(np.abs(native_mass), np.abs(bulk_mass))
        mismatch = np.abs(native_mass - bulk_mass)
        if np.any(mismatch > _NATIVE_MASS_RELATIVE_TOLERANCE * scale):
            bad = int(np.argmax(mismatch / np.maximum(scale, np.finfo(np.float64).tiny)))
            raise ValueError(
                "native_point_mass_kg does not match T4-derived lumped bulk mass "
                f"at node {bad}: native={native_mass[bad]:.17g} kg, derived={bulk_mass[bad]:.17g} kg"
            )
        mass = native_mass
        mass_source = "NATIVE_POINT_MASS_VERIFIED_AGAINST_T4_DERIVED"

    free_mask = ~mask.reshape(-1)
    free_dofs = np.flatnonzero(free_mask)
    if free_dofs.size == 0:
        raise ValueError("at least one free DOF is required")
    free_nodes = free_dofs // 3
    if np.any(mass[free_nodes] <= 0.0):
        bad_node = int(free_nodes[np.flatnonzero(mass[free_nodes] <= 0.0)[0]])
        raise ValueError(f"free DOF at node {bad_node} has zero or negative point mass")

    bulk_rows = np.zeros(total_dof_count, dtype=np.float64)
    cohesive_rows = np.zeros(total_dof_count, dtype=np.float64)
    dense_stiffness = np.zeros((total_dof_count, total_dof_count), dtype=np.float64) if validate_dense else None

    def add_local(
        local_stiffness: np.ndarray,
        local_nodes: np.ndarray,
        row_accumulator: np.ndarray,
    ) -> None:
        local_dofs = (3 * local_nodes[:, None] + np.arange(3, dtype=np.int64)).reshape(-1)
        local_free = free_mask[local_dofs]
        local_indices = np.flatnonzero(local_free)
        if local_indices.size:
            global_free_dofs = local_dofs[local_indices]
            local_free_nodes = global_free_dofs // 3
            root_mass = np.sqrt(mass[local_free_nodes])
            reduced_abs = np.abs(local_stiffness[np.ix_(local_indices, local_indices)])
            reduced_abs = reduced_abs / root_mass[:, None] / root_mass[None, :]
            if not np.all(np.isfinite(reduced_abs)):
                raise ValueError("mass-normalized local stiffness overflowed")
            row_accumulator[global_free_dofs] += reduced_abs.sum(axis=1)
        if dense_stiffness is not None:
            dense_stiffness[np.ix_(local_dofs, local_dofs)] += local_stiffness

    for nodes_in_tet, stiffness in zip(elements, tetra_stiffnesses, strict=True):
        add_local(stiffness, nodes_in_tet, bulk_rows)

    areas = np.empty(face_connectivity.shape[0], dtype=np.float64)
    for index, pair in enumerate(face_connectivity):
        coordinates = nodes[pair]
        stiffness, area, _normal = tri3_cohesive_stiffness_matrix(coordinates, kn, kt)
        areas[index] = area
        add_local(stiffness, pair.reshape(-1), cohesive_rows)
    if not np.all(np.isfinite(bulk_rows)) or not np.all(np.isfinite(cohesive_rows)):
        raise ValueError("row-sum stiffness bound overflowed")

    total_rows = bulk_rows + cohesive_rows
    lambda_bound = float(np.max(total_rows[free_dofs]))
    if not math.isfinite(lambda_bound) or lambda_bound <= 0.0:
        raise ValueError("computed stiffness upper bound must be finite and positive")
    critical_dt = 2.0 / math.sqrt(lambda_bound)
    dense_max: float | None = None
    dense_min: float | None = None
    if dense_stiffness is not None:
        sqrt_mass = np.sqrt(mass[free_nodes])
        free_stiffness = dense_stiffness[np.ix_(free_dofs, free_dofs)]
        normalized = free_stiffness / sqrt_mass[:, None] / sqrt_mass[None, :]
        normalized = 0.5 * (normalized + normalized.T)
        eigenvalues = np.linalg.eigvalsh(normalized)
        dense_min = float(eigenvalues[0])
        dense_max = float(eigenvalues[-1])
        if not math.isfinite(dense_min) or not math.isfinite(dense_max):
            raise ValueError("dense eigenvalue diagnostic returned non-finite values")

    return CohesiveStabilityBound(
        lambda_upper_bound_per_s2=lambda_bound,
        critical_time_step_s=critical_dt,
        recommended_time_step_s=0.5 * critical_dt,
        bulk_lumped_node_mass_kg=_readonly(bulk_mass),
        mass_used_node_kg=_readonly(mass),
        bulk_total_lumped_mass_kg=float(bulk_mass.sum()),
        cohesive_face_areas_m2=_readonly(areas),
        cohesive_total_area_m2=float(areas.sum()),
        bulk_row_sum_contributions_per_s2=_readonly(bulk_rows.reshape((node_count, 3))),
        cohesive_row_sum_contributions_per_s2=_readonly(cohesive_rows.reshape((node_count, 3))),
        free_dof_count=int(free_dofs.size),
        blocked_dof_count=int(total_dof_count - free_dofs.size),
        total_dof_count=total_dof_count,
        mass_source=mass_source,
        dense_max_eigenvalue_per_s2=dense_max,
        dense_min_eigenvalue_per_s2=dense_min,
    )
