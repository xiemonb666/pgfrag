"""External-field Mode III J-integral diagnostic for a straight z-front.

This module consumes already sampled contour data.  It does not solve an
elasticity problem, reconstruct a displacement field, or infer a stress
intensity factor.  See README.md for the input contract and derivation.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping

import numpy as np


MAX_SAMPLES = 256
NORMAL_TOLERANCE = 1.0e-8
PLANE_TOLERANCE_M = 1.0e-12
Z_GRADIENT_TOLERANCE = 1.0e-12
NORMAL_CLOSURE_REL_TOLERANCE = 5.0e-2

INPUT_UNITS = {
    "coordinates": "m",
    "normal": "1",
    "line_weight": "m",
    "grad_u_z": "m/m",
    "shear_modulus": "Pa",
}

OUTPUT_UNITS = {
    "energy_density_W": "J/m^3",
    "traction_t_z": "Pa",
    "integrand_density": "Pa",
    "line_weight": "m",
    "J": "J/m^2",
}


@dataclass(frozen=True, slots=True)
class JIntegralResult:
    """A status-bearing result. UNKNOWN results never carry a numeric J."""

    status: str
    J_J_per_m2: float | None
    reason: str | None
    details: dict[str, Any]

    def to_dict(self) -> dict[str, Any]:
        return {
            "status": self.status,
            "J_J_per_m2": self.J_J_per_m2,
            "reason": self.reason,
            **self.details,
        }


def _unknown(reason: str, **details: Any) -> JIntegralResult:
    return JIntegralResult(
        status="UNKNOWN",
        J_J_per_m2=None,
        reason=reason,
        details={"input_units_expected": dict(INPUT_UNITS), **details},
    )


def _preflight_row_shape(name: str, value: object, shape_tail: tuple[int, ...]) -> int:
    """Check sample count and trailing row shape without materializing an array."""
    if value is None:
        raise ValueError(f"missing_{name}")
    if isinstance(value, np.ndarray):
        shape = tuple(value.shape)
    elif hasattr(value, "shape"):
        try:
            shape = tuple(value.shape)  # type: ignore[attr-defined]
        except (TypeError, ValueError) as exc:
            raise ValueError(f"invalid_{name}_shape") from exc
    elif hasattr(value, "__len__") and hasattr(value, "__getitem__"):
        row_count = len(value)  # type: ignore[arg-type]
        if not 1 <= row_count <= MAX_SAMPLES:
            raise ValueError("sample_count_out_of_range")
        for row_index in range(row_count):
            row = value[row_index]  # type: ignore[index]
            if shape_tail:
                if isinstance(row, np.ndarray):
                    row_shape = tuple(row.shape)
                elif isinstance(row, (list, tuple)):
                    if any(isinstance(item, (list, tuple, np.ndarray)) for item in row):
                        row_shape = ()
                    else:
                        row_shape = (len(row),)
                else:
                    row_shape = ()
                if row_shape != shape_tail:
                    raise ValueError(f"invalid_{name}_shape")
            elif isinstance(row, (list, tuple)) or (isinstance(row, np.ndarray) and row.ndim != 0):
                raise ValueError(f"invalid_{name}_shape")
        shape = (row_count, *shape_tail)
    else:
        raise ValueError(f"invalid_{name}_type")

    if len(shape) != 1 + len(shape_tail) or tuple(shape[1:]) != shape_tail:
        raise ValueError(f"invalid_{name}_shape")
    row_count = int(shape[0])
    if not 1 <= row_count <= MAX_SAMPLES:
        raise ValueError("sample_count_out_of_range")
    return row_count


def _numeric_array(name: str, value: object, shape_tail: tuple[int, ...]) -> np.ndarray | None:
    _preflight_row_shape(name, value, shape_tail)
    try:
        array = np.asarray(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"invalid_{name}_array") from exc
    if array.dtype.kind not in "fiu" or array.dtype.itemsize > 8:
        raise ValueError(f"invalid_{name}_dtype")
    if array.ndim != 1 + len(shape_tail) or array.shape[1:] != shape_tail:
        raise ValueError(f"invalid_{name}_shape")
    if array.dtype.kind in "iu":
        exact_limit = 1 << 53
        if bool(np.any(array > exact_limit)) or (array.dtype.kind == "i" and bool(np.any(array < -exact_limit))):
            raise ValueError(f"inexact_integer_{name}")
    try:
        converted = array.astype(np.float64, copy=True)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"invalid_{name}_numeric_conversion") from exc
    if not bool(np.isfinite(converted).all()):
        raise ValueError(f"nonfinite_{name}")
    return converted


def integrate_modeiii_j(
    contour_coordinates_m: object = None,
    normals: object = None,
    weights_m: object = None,
    grad_u_z_m_per_m: object = None,
    *,
    shear_modulus_Pa: object = None,
    units: Mapping[str, str] | None = None,
    field_provenance: str | None = None,
    one_sided_crack_face: str | None = None,
) -> JIntegralResult:
    """Integrate J = ∮[W n_x - t_z u_z,x] ds from external field samples.

    Arrays are query-major and ordered along the contour, with shapes coordinates=(n,3), normals=(n,3),
    weights=(n,), and grad_u_z=(n,3).  Coordinates lie in one x-y plane
    perpendicular to the straight z-axis front.  Normals are outward unit
    normals; positive weights are physical ds values in meters.  A point on
    the exact branch face x<0,y=0 is accepted only when the caller labels the
    supplied trace explicitly as ``upper`` or ``lower``.  No coordinate is
    shifted or regularized here.

    Missing, malformed, nonfinite, or convention-inconsistent inputs produce
    status UNKNOWN and no numeric J.  The declared field provenance is
    recorded but cannot be independently verified by this integrator.
    """

    if units is None:
        return _unknown("missing_units")
    if not isinstance(units, Mapping):
        return _unknown("invalid_units_declaration")
    if dict(units) != INPUT_UNITS:
        return _unknown("invalid_or_incomplete_units", declared_units=dict(units))
    if not isinstance(field_provenance, str) or not field_provenance.strip():
        return _unknown("missing_field_provenance")
    if one_sided_crack_face not in (None, "upper", "lower"):
        return _unknown("invalid_one_sided_crack_face", one_sided_crack_face=one_sided_crack_face)

    try:
        # Check all four shapes/budgets before np.asarray or astype(copy=True)
        # can materialize any caller-provided sequence.
        preflight_counts = {}
        for name, value, tail in (
            ("contour_coordinates", contour_coordinates_m, (3,)),
            ("normals", normals, (3,)),
            ("weights", weights_m, ()),
            ("grad_u_z", grad_u_z_m_per_m, (3,)),
        ):
            preflight_counts[name] = _preflight_row_shape(name, value, tail)
        if len(set(preflight_counts.values())) != 1:
            return _unknown("sample_count_mismatch", sample_counts=preflight_counts)
        xyz = _numeric_array("contour_coordinates", contour_coordinates_m, (3,))
        normal = _numeric_array("normals", normals, (3,))
        weight = _numeric_array("weights", weights_m, ())
        gradient = _numeric_array("grad_u_z", grad_u_z_m_per_m, (3,))
    except ValueError as exc:
        return _unknown(str(exc))
    assert xyz is not None and normal is not None and weight is not None and gradient is not None
    n = int(xyz.shape[0])
    if not 1 <= n <= MAX_SAMPLES:
        return _unknown("sample_count_out_of_range", sample_count=n, maximum_samples=MAX_SAMPLES)
    if normal.shape[0] != n or weight.shape[0] != n or gradient.shape[0] != n:
        return _unknown(
            "sample_count_mismatch",
            sample_counts={
                "coordinates": n,
                "normals": int(normal.shape[0]),
                "weights": int(weight.shape[0]),
                "grad_u_z": int(gradient.shape[0]),
            },
        )
    if shear_modulus_Pa is None:
        return _unknown("missing_or_nonfinite_shear_modulus")
    if isinstance(shear_modulus_Pa, (bool, np.bool_)):
        return _unknown("invalid_shear_modulus")
    try:
        mu = float(shear_modulus_Pa)
    except (TypeError, ValueError, OverflowError):
        return _unknown("invalid_shear_modulus")
    if not math.isfinite(mu) or mu <= 0.0:
        return _unknown("shear_modulus_must_be_positive_finite")
    if bool(np.any(weight <= 0.0)):
        return _unknown("line_weights_must_be_positive", sample_count=n)

    normal_lengths = np.linalg.norm(normal, axis=1)
    if bool(np.any(np.abs(normal_lengths - 1.0) > NORMAL_TOLERANCE)):
        return _unknown("normals_must_be_unit_length")
    if bool(np.any(np.abs(normal[:, 2]) > NORMAL_TOLERANCE)):
        return _unknown("normals_must_lie_in_xy_plane")
    coordinate_scale = max(1.0, float(np.max(np.abs(xyz))))
    plane_tol = PLANE_TOLERANCE_M * coordinate_scale
    if float(np.ptp(xyz[:, 2])) > plane_tol:
        return _unknown("contour_must_lie_in_one_xy_plane")
    gradient_scale = max(1.0, float(np.max(np.abs(gradient))))
    if bool(np.any(np.abs(gradient[:, 2]) > Z_GRADIENT_TOLERANCE * gradient_scale)):
        return _unknown("field_must_be_z_invariant")

    front_points = (xyz[:, 0] == 0.0) & (xyz[:, 1] == 0.0)
    if bool(np.any(front_points)):
        return _unknown("exact_front_samples_are_singular", count=int(np.count_nonzero(front_points)))
    branch_points = (xyz[:, 0] < 0.0) & (xyz[:, 1] == 0.0)
    if bool(np.any(branch_points)) and one_sided_crack_face is None:
        return _unknown(
            "branch_face_requires_explicit_one_sided_trace",
            count=int(np.count_nonzero(branch_points)),
        )

    # These two checks catch the common inward-normal and open/inconsistent
    # contour errors without trying to reconstruct a contour from quadrature
    # points.  They are diagnostics, not a proof of contour topology.
    perimeter = float(np.sum(weight))
    normal_closure = np.sum(normal * weight[:, None], axis=0)
    closure_relative = float(np.linalg.norm(normal_closure) / perimeter)
    if closure_relative > NORMAL_CLOSURE_REL_TOLERANCE:
        return _unknown(
            "weighted_normals_do_not_close_contour",
            normal_closure_m=[float(v) for v in normal_closure],
            normal_closure_relative=closure_relative,
        )
    signed_area = 0.5 * float(np.sum((xyz[:, 0] * normal[:, 0] + xyz[:, 1] * normal[:, 1]) * weight))
    if not math.isfinite(signed_area) or signed_area <= 0.0:
        return _unknown("normals_must_point_outward_around_front", signed_area_m2=signed_area)
    xy = xyz[:, :2]
    next_xy = np.roll(xy, -1, axis=0)
    cross = xy[:, 0] * next_xy[:, 1] - xy[:, 1] * next_xy[:, 0]
    dot = np.sum(xy * next_xy, axis=1)
    winding_angle = float(np.sum(np.arctan2(cross, dot)))
    if abs(winding_angle) < math.pi:
        return _unknown(
            "ordered_contour_samples_do_not_encircle_front",
            winding_angle_rad=winding_angle,
        )

    grad_x = gradient[:, 0]
    grad_y = gradient[:, 1]
    energy_density = 0.5 * mu * (grad_x * grad_x + grad_y * grad_y)
    traction_z = mu * (gradient[:, 0] * normal[:, 0] + gradient[:, 1] * normal[:, 1])
    energy_term = energy_density * normal[:, 0]
    traction_term = -traction_z * grad_x
    integrand_density = energy_term + traction_term
    contribution = integrand_density * weight
    J = float(np.sum(contribution))
    if not all(bool(np.isfinite(a).all()) for a in (energy_density, traction_z, integrand_density, contribution)) or not math.isfinite(J):
        return _unknown("nonfinite_integral_result")

    samples = []
    for i in range(n):
        samples.append(
            {
                "coordinate_m": [float(v) for v in xyz[i]],
                "normal_outward": [float(v) for v in normal[i]],
                "weight_ds_m": float(weight[i]),
                "grad_u_z_m_per_m": [float(v) for v in gradient[i]],
                "W_J_per_m3": float(energy_density[i]),
                "t_z_Pa": float(traction_z[i]),
                "W_n_x_Pa": float(energy_term[i]),
                "minus_t_z_u_z_x_Pa": float(traction_term[i]),
                "integrand_Pa": float(integrand_density[i]),
                "weighted_contribution_J_per_m2": float(contribution[i]),
            }
        )

    details = {
        "formula": "J = integral_[Gamma] (W*n_x - t_z*du_z/dx) ds",
        "assumptions": [
            "straight z-axis crack front",
            "z-invariant antiplane Mode III field, u_x=u_y=0",
            "linear isotropic elasticity with supplied positive shear modulus",
            "closed contour around the front with outward normals",
        ],
        "normal_convention": "outward from the enclosed front region; positive signed area required",
        "field_provenance": field_provenance.strip(),
        "field_provenance_verified": False,
        "one_sided_crack_face": one_sided_crack_face,
        "sample_count": n,
        "shear_modulus_Pa": mu,
        "units": {**INPUT_UNITS, **OUTPUT_UNITS},
        "geometry_diagnostics": {
            "perimeter_weight_sum_m": perimeter,
            "weighted_normal_closure_m": [float(v) for v in normal_closure],
            "weighted_normal_closure_relative": closure_relative,
            "signed_area_from_normals_m2": signed_area,
            "winding_angle_around_front_rad": winding_angle,
            "contour_z_m": float(np.mean(xyz[:, 2])),
        },
        "energy_term_integral_J_per_m2": float(np.sum(energy_term * weight)),
        "traction_term_integral_J_per_m2": float(np.sum(traction_term * weight)),
        "samples": samples,
    }
    return JIntegralResult(status="OK", J_J_per_m2=J, reason=None, details=details)

