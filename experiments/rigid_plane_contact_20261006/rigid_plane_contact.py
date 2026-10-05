"""Small SI-unit API for frictionless unilateral point/rigid-plane contact.

This module evaluates one fixed plane's penalty force and potential. It does not
perform body/fragment broad phase, surface collision, rotation dynamics, or
multi-contact resolution.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np


_UNIT_NORMAL_TOL = 1.0e-10


def _si_scalar(name: str, value: Any, *, positive: bool = False) -> float:
    """Validate a bare finite scalar expressed in the SI unit in ``name``."""
    if hasattr(value, "units"):
        raise TypeError(f"{name} must be a bare SI scalar; convert quantities before calling")
    if isinstance(value, (str, bytes, bool, np.bool_)) or not np.isscalar(value):
        raise TypeError(f"{name} must be a bare numeric SI scalar")
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise TypeError(f"{name} must be a bare numeric SI scalar") from exc
    if not np.isfinite(result):
        raise ValueError(f"{name} must be finite")
    if positive and result <= 0.0:
        raise ValueError(f"{name} must be positive")
    return result


def _vector3(name: str, value: Any) -> np.ndarray:
    """Return a finite three-vector; vector values are expressed in SI units."""
    if hasattr(value, "units"):
        raise TypeError(f"{name} must be a bare SI vector; convert quantities before calling")
    try:
        result = np.asarray(value, dtype=float)
    except (TypeError, ValueError, OverflowError) as exc:
        raise TypeError(f"{name} must be a numeric three-vector in SI units") from exc
    if result.shape != (3,):
        raise ValueError(f"{name} must have shape (3,)")
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain only finite values")
    return result.copy()


def _bytes_backed_readonly_vector3(value: np.ndarray) -> np.ndarray:
    """Create an immutable vector whose read-only flag cannot be re-enabled."""
    immutable_bytes = value.astype(np.float64, copy=False).tobytes(order="C")
    return np.frombuffer(immutable_bytes, dtype=np.float64).reshape((3,))


@dataclass(frozen=True, slots=True)
class RigidPlanePointContact:
    """Frictionless elastic penalty contact with one fixed infinite plane.

    ``unit_normal`` points from the plane into the permitted (free) half-space.
    All dimensional arguments are bare SI values: metres, newtons per metre.
    The normal is validated as unit length and is not silently renormalized.
    """

    plane_point_m: Any
    unit_normal: Any
    stiffness_N_per_m: Any

    def __post_init__(self) -> None:
        point = _vector3("plane_point_m", self.plane_point_m)
        normal = _vector3("unit_normal", self.unit_normal)
        norm = float(np.linalg.norm(normal))
        if not np.isfinite(norm) or norm == 0.0:
            raise ValueError("unit_normal must be finite and nonzero")
        if abs(norm - 1.0) > _UNIT_NORMAL_TOL:
            raise ValueError("unit_normal must have length 1 within 1e-10")
        stiffness = _si_scalar("stiffness_N_per_m", self.stiffness_N_per_m, positive=True)
        point = _bytes_backed_readonly_vector3(point)
        normal = _bytes_backed_readonly_vector3(normal)
        object.__setattr__(self, "plane_point_m", point)
        object.__setattr__(self, "unit_normal", normal)
        object.__setattr__(self, "stiffness_N_per_m", stiffness)

    def signed_gap_m(self, position_m: Any) -> float:
        """Signed plane distance; positive is free space and negative penetrates."""
        position = _vector3("position_m", position_m)
        try:
            with np.errstate(over="raise", invalid="raise"):
                gap = float(np.dot(position - self.plane_point_m, self.unit_normal))
        except FloatingPointError as exc:
            raise FloatingPointError("signed gap overflowed for the supplied finite coordinates") from exc
        if not np.isfinite(gap):
            raise FloatingPointError("signed gap is not finite")
        return gap

    def force_N(self, position_m: Any) -> np.ndarray:
        """Return the contact force on the point mass, in newtons."""
        gap = self.signed_gap_m(position_m)
        if gap >= 0.0:
            return np.zeros(3, dtype=float)
        try:
            with np.errstate(over="raise", invalid="raise"):
                force = (-self.stiffness_N_per_m * gap) * self.unit_normal
        except FloatingPointError as exc:
            raise FloatingPointError("contact force overflowed for the supplied finite inputs") from exc
        if not np.all(np.isfinite(force)):
            raise FloatingPointError("contact force is not finite")
        return force

    def potential_J(self, position_m: Any) -> float:
        """Return stored elastic penalty energy, in joules."""
        gap = self.signed_gap_m(position_m)
        if gap >= 0.0:
            return 0.0
        try:
            with np.errstate(over="raise", invalid="raise"):
                energy = 0.5 * self.stiffness_N_per_m * gap * gap
        except FloatingPointError as exc:
            raise FloatingPointError("contact potential overflowed for the supplied finite inputs") from exc
        if not np.isfinite(energy):
            raise FloatingPointError("contact potential is not finite")
        return float(energy)


def signed_gap_m(position_m: Any, plane_point_m: Any, unit_normal: Any) -> float:
    """Convenience geometry API, with the same finite SI and unit-normal checks."""
    return RigidPlanePointContact(plane_point_m, unit_normal, 1.0).signed_gap_m(position_m)


def contact_force_N(
    position_m: Any, plane_point_m: Any, unit_normal: Any, stiffness_N_per_m: Any
) -> np.ndarray:
    """Convenience force API for one point and one fixed plane."""
    return RigidPlanePointContact(plane_point_m, unit_normal, stiffness_N_per_m).force_N(position_m)


def contact_potential_J(
    position_m: Any, plane_point_m: Any, unit_normal: Any, stiffness_N_per_m: Any
) -> float:
    """Convenience potential API for one point and one fixed plane."""
    return RigidPlanePointContact(plane_point_m, unit_normal, stiffness_N_per_m).potential_J(position_m)

