"""Seeded, continuous heterogeneity bound to material rest coordinates.

This experiment module is deliberately independent of PGFrag's solver and
production material contracts.  The field is a finite random Fourier
realization whose coefficients are generated once from ``seed`` and then
evaluated only from physical coordinates in metres.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
import math
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np


SCHEMA = "pgfrag.rest_material_field.manifest.v1"
FIELD_STATUS = "UNCALIBRATED"
REST_FRAME_ID = "material-rest-cartesian-xyz"
_E_FOLD_SHELL_Q = 2.199123071161498
_DEFAULT_MATERIAL_UNITS = (
    ("heterogeneity", "dimensionless"),
    ("positive_parameter_ratio", "dimensionless"),
)


class MaterialFieldValidationError(ValueError):
    """Raised when field parameters or evaluation coordinates are invalid."""


def _finite_float(value: Any, name: str) -> float:
    if isinstance(value, (bool, np.bool_)):
        raise MaterialFieldValidationError(f"{name} must be a finite real number")
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise MaterialFieldValidationError(f"{name} must be a finite real number") from exc
    if not math.isfinite(result):
        raise MaterialFieldValidationError(f"{name} must be finite")
    return result


def _point3(point_m: Sequence[float], name: str = "point_m") -> tuple[float, float, float]:
    try:
        array = np.asarray(point_m, dtype=np.float64)
    except (TypeError, ValueError, OverflowError) as exc:
        raise MaterialFieldValidationError(f"{name} must contain three finite coordinates in metres") from exc
    if array.shape != (3,) or not np.isfinite(array).all():
        raise MaterialFieldValidationError(f"{name} must contain three finite coordinates in metres")
    return (float(array[0]), float(array[1]), float(array[2]))


def _points3(points_m: Sequence[Sequence[float]], name: str = "points_m") -> tuple[tuple[float, float, float], ...]:
    try:
        array = np.asarray(points_m, dtype=np.float64)
    except (TypeError, ValueError, OverflowError) as exc:
        raise MaterialFieldValidationError(f"{name} must be an N-by-3 array of finite coordinates in metres") from exc
    if array.ndim != 2 or array.shape[1] != 3 or not np.isfinite(array).all():
        raise MaterialFieldValidationError(f"{name} must be an N-by-3 array of finite coordinates in metres")
    return tuple((float(row[0]), float(row[1]), float(row[2])) for row in array)


def _rotation3(value: Sequence[Sequence[float]]) -> tuple[tuple[float, float, float], ...]:
    try:
        matrix = np.asarray(value, dtype=np.float64)
    except (TypeError, ValueError, OverflowError) as exc:
        raise MaterialFieldValidationError("rotation must be a finite proper 3-by-3 matrix") from exc
    if matrix.shape != (3, 3) or not np.isfinite(matrix).all():
        raise MaterialFieldValidationError("rotation must be a finite proper 3-by-3 matrix")
    if not np.allclose(matrix.T @ matrix, np.eye(3), rtol=0.0, atol=1e-10):
        raise MaterialFieldValidationError("rotation must be orthonormal; scale and shear are not allowed")
    determinant = float(np.linalg.det(matrix))
    if not math.isclose(determinant, 1.0, rel_tol=0.0, abs_tol=1e-10):
        raise MaterialFieldValidationError("rotation must be right-handed with determinant +1")
    return tuple(tuple(float(x) for x in row) for row in matrix)


@dataclass(frozen=True)
class MaterialFieldParameters:
    """Immutable parameters for a three-dimensional rest-frame field.

    ``amplitude`` is a maximum relative deviation, not an RMS value.  Keeping
    it in [0, 1) makes the derived multiplicative map strictly positive.
    """

    seed: int
    correlation_length_m: float
    amplitude: float = 0.30
    mode_count: int = 24
    coordinate_units: str = "m"
    correlation_length_units: str = "m"
    reference_frame_id: str = REST_FRAME_ID
    material_units: tuple[tuple[str, str], ...] | Mapping[str, str] = _DEFAULT_MATERIAL_UNITS

    def __post_init__(self) -> None:
        if isinstance(self.seed, (bool, np.bool_)) or not isinstance(self.seed, (int, np.integer)):
            raise MaterialFieldValidationError("seed must be an integer")
        if int(self.seed) < 0 or int(self.seed) > (2**63 - 1):
            raise MaterialFieldValidationError("seed must be in [0, 2^63-1]")
        object.__setattr__(self, "seed", int(self.seed))

        ell = _finite_float(self.correlation_length_m, "correlation_length_m")
        if ell <= 0.0:
            raise MaterialFieldValidationError("correlation_length_m must be positive")
        object.__setattr__(self, "correlation_length_m", ell)

        amplitude = _finite_float(self.amplitude, "amplitude")
        if not 0.0 <= amplitude < 1.0:
            raise MaterialFieldValidationError("amplitude must be in [0, 1) to preserve positive mapped parameters")
        object.__setattr__(self, "amplitude", amplitude)

        if isinstance(self.mode_count, (bool, np.bool_)) or not isinstance(self.mode_count, (int, np.integer)):
            raise MaterialFieldValidationError("mode_count must be an integer in [4, 128]")
        if not 4 <= int(self.mode_count) <= 128:
            raise MaterialFieldValidationError("mode_count must be an integer in [4, 128]")
        object.__setattr__(self, "mode_count", int(self.mode_count))

        if self.coordinate_units != "m":
            raise MaterialFieldValidationError("coordinate_units must be exactly 'm'; convert coordinates explicitly")
        if self.correlation_length_units != "m":
            raise MaterialFieldValidationError("correlation_length_units must be exactly 'm'")
        if not isinstance(self.reference_frame_id, str) or not self.reference_frame_id.strip():
            raise MaterialFieldValidationError("reference_frame_id must be a non-empty string")

        units = self.material_units
        if isinstance(units, Mapping):
            pairs = tuple(sorted((str(key), str(unit)) for key, unit in units.items()))
        else:
            try:
                pairs = tuple(sorted((str(key), str(unit)) for key, unit in units))
            except (TypeError, ValueError) as exc:
                raise MaterialFieldValidationError("material_units must be a mapping or pairs of strings") from exc
        if pairs != _DEFAULT_MATERIAL_UNITS:
            raise MaterialFieldValidationError(
                "material_units must declare dimensionless heterogeneity and positive_parameter_ratio"
            )
        object.__setattr__(self, "material_units", pairs)


@dataclass(frozen=True)
class RigidFrameTransform:
    """Rigid transform from canonical material rest coordinates to query frame.

    ``point_frame = rotation_rest_to_frame @ point_rest + translation_rest_to_frame_m``.
    Evaluation in the query frame applies the inverse transform before looking
    up the field, so a rigid motion cannot move a material defect through the
    material.
    """

    rotation_rest_to_frame: tuple[tuple[float, float, float], ...] | Sequence[Sequence[float]]
    translation_rest_to_frame_m: tuple[float, float, float] | Sequence[float]
    frame_id: str = "query-frame"
    rest_frame_id: str = REST_FRAME_ID

    def __post_init__(self) -> None:
        rotation = _rotation3(self.rotation_rest_to_frame)
        translation = _point3(self.translation_rest_to_frame_m, "translation_rest_to_frame_m")
        if not isinstance(self.frame_id, str) or not self.frame_id.strip():
            raise MaterialFieldValidationError("frame_id must be a non-empty string")
        if not isinstance(self.rest_frame_id, str) or not self.rest_frame_id.strip():
            raise MaterialFieldValidationError("rest_frame_id must be a non-empty string")
        object.__setattr__(self, "rotation_rest_to_frame", rotation)
        object.__setattr__(self, "translation_rest_to_frame_m", translation)

    def to_rest(self, point_frame_m: Sequence[float]) -> tuple[float, float, float]:
        point = np.asarray(_point3(point_frame_m, "point_frame_m"), dtype=np.float64)
        rotation = np.asarray(self.rotation_rest_to_frame, dtype=np.float64)
        translation = np.asarray(self.translation_rest_to_frame_m, dtype=np.float64)
        rest = rotation.T @ (point - translation)
        return (float(rest[0]), float(rest[1]), float(rest[2]))

    def to_rest_many(self, points_frame_m: Sequence[Sequence[float]]) -> tuple[tuple[float, float, float], ...]:
        points = np.asarray(_points3(points_frame_m, "points_frame_m"), dtype=np.float64)
        if len(points) == 0:
            return ()
        rotation = np.asarray(self.rotation_rest_to_frame, dtype=np.float64)
        translation = np.asarray(self.translation_rest_to_frame_m, dtype=np.float64)
        rest = (points - translation) @ rotation
        return tuple((float(row[0]), float(row[1]), float(row[2])) for row in rest)


@dataclass(frozen=True)
class FourierMode:
    wavevector_rad_per_m: tuple[float, float, float]
    phase_rad: float


@dataclass(frozen=True)
class RestFrameMaterialField:
    """A frozen, bounded scalar field tied to immutable reference coordinates."""

    parameters: MaterialFieldParameters
    modes: tuple[FourierMode, ...] = field(init=False, repr=False)
    field_id: str = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.parameters, MaterialFieldParameters):
            raise MaterialFieldValidationError("parameters must be MaterialFieldParameters")

        rng = np.random.Generator(np.random.PCG64(self.parameters.seed))
        directions = rng.standard_normal((self.parameters.mode_count, 3))
        norms = np.linalg.norm(directions, axis=1)
        # A Gaussian direction has zero norm with probability zero; the guard
        # keeps the contract explicit even under a substituted RNG.
        if not np.isfinite(norms).all() or np.any(norms == 0.0):
            raise MaterialFieldValidationError("seed produced an invalid Fourier direction")
        directions = directions / norms[:, None]
        wave_number = _E_FOLD_SHELL_Q / self.parameters.correlation_length_m
        phases = rng.uniform(0.0, 2.0 * math.pi, self.parameters.mode_count)
        modes = tuple(
            FourierMode(
                wavevector_rad_per_m=(float(wave_number * direction[0]), float(wave_number * direction[1]), float(wave_number * direction[2])),
                phase_rad=float(phase),
            )
            for direction, phase in zip(directions, phases, strict=True)
        )
        object.__setattr__(self, "modes", modes)
        payload = self._manifest_payload()
        canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
        object.__setattr__(self, "field_id", hashlib.sha256(canonical.encode("utf-8")).hexdigest())

    def _normalized_value(self, point_m: tuple[float, float, float]) -> float:
        terms = []
        for mode in self.modes:
            phase = math.fsum(k * x for k, x in zip(mode.wavevector_rad_per_m, point_m, strict=True))
            terms.append(math.cos(phase + mode.phase_rad))
        return max(-1.0, min(1.0, math.fsum(terms) / len(self.modes)))

    def evaluate(self, point_m: Sequence[float]) -> float:
        """Return signed dimensionless heterogeneity at a rest-frame point.

        The input coordinates are in metres in ``parameters.reference_frame_id``.
        This method never derives randomness from a point, node id, or sample order.
        """

        point = _point3(point_m)
        if self.parameters.amplitude == 0.0:
            return 0.0
        return self.parameters.amplitude * self._normalized_value(point)

    def evaluate_many(self, points_m: Sequence[Sequence[float]]) -> tuple[float, ...]:
        points = _points3(points_m)
        if self.parameters.amplitude == 0.0:
            return tuple(0.0 for _ in points)
        return tuple(self.parameters.amplitude * self._normalized_value(point) for point in points)

    def evaluate_in_frame(self, point_frame_m: Sequence[float], transform: RigidFrameTransform) -> float:
        if not isinstance(transform, RigidFrameTransform):
            raise MaterialFieldValidationError("transform must be RigidFrameTransform")
        if transform.rest_frame_id != self.parameters.reference_frame_id:
            raise MaterialFieldValidationError("transform rest_frame_id does not match the field reference frame")
        return self.evaluate(transform.to_rest(point_frame_m))

    def evaluate_many_in_frame(
        self,
        points_frame_m: Sequence[Sequence[float]],
        transform: RigidFrameTransform,
    ) -> tuple[float, ...]:
        if not isinstance(transform, RigidFrameTransform):
            raise MaterialFieldValidationError("transform must be RigidFrameTransform")
        if transform.rest_frame_id != self.parameters.reference_frame_id:
            raise MaterialFieldValidationError("transform rest_frame_id does not match the field reference frame")
        return self.evaluate_many(transform.to_rest_many(points_frame_m))

    def positive_parameter_bounds(self, nominal_value: float) -> tuple[float, float]:
        """Return finite bounds for a caller-supplied, positive nominal value."""

        nominal = _finite_float(nominal_value, "nominal_value")
        if nominal <= 0.0:
            raise MaterialFieldValidationError("nominal_value must be positive")
        lower = nominal * (1.0 - self.parameters.amplitude)
        upper = nominal * (1.0 + self.parameters.amplitude)
        if not math.isfinite(lower) or not math.isfinite(upper):
            raise MaterialFieldValidationError("nominal_value is too large for finite mapped parameter bounds")
        return lower, upper

    def evaluate_positive_parameter(self, nominal_value: float, point_m: Sequence[float]) -> float:
        """Map a caller-supplied nominal positive parameter through a bounded ratio.

        ``p(x) = p0 * (1 + h(x))`` with ``h in [-amplitude, amplitude]``.
        The field does not provide nominal E, Gc, strength, or any other
        physically calibrated parameter.
        """

        nominal = _finite_float(nominal_value, "nominal_value")
        if nominal <= 0.0:
            raise MaterialFieldValidationError("nominal_value must be positive")
        point = _point3(point_m)
        value = nominal * (1.0 + self.evaluate(point))
        if not math.isfinite(value) or value <= 0.0:
            raise MaterialFieldValidationError("mapped parameter must remain finite and positive")
        return value

    def _manifest_payload(self) -> dict[str, Any]:
        ell = self.parameters.correlation_length_m
        wave_number = _E_FOLD_SHELL_Q / ell
        return {
            "schema": SCHEMA,
            "field_kind": "seeded_finite_random_fourier_mode_scalar",
            "coordinate_system": {
                "frame_id": self.parameters.reference_frame_id,
                "binding": "material_rest_coordinates",
                "coordinate_units": self.parameters.coordinate_units,
                "dimension": 3,
            },
            "parameters": {
                "seed": self.parameters.seed,
                "correlation_length_m": ell,
                "correlation_length_units": self.parameters.correlation_length_units,
                "amplitude_max_relative_deviation": self.parameters.amplitude,
                "mode_count": self.parameters.mode_count,
                "material_units": dict(self.parameters.material_units),
            },
            "field_contract": {
                "normalized_range": [-1.0, 1.0],
                "heterogeneity_range": [-self.parameters.amplitude, self.parameters.amplitude],
                "positive_parameter_map": "p(x) = p0 * (1 + h(x)); caller-supplied p0 must be positive",
                "positive_parameter_ratio_bounds": [1.0 - self.parameters.amplitude, 1.0 + self.parameters.amplitude],
                "continuous": True,
                "point_randomness": False,
                "mesh_index_randomness": False,
            },
            "spectrum_and_correlation": {
                "dimension": 3,
                "spectrum": "equal-amplitude isotropic random directions on one wave-number shell",
                "wave_number_magnitude_rad_per_m": wave_number,
                "ensemble_covariance": "C(r) = sin(k*r)/(k*r), for independent uniform phases and isotropic directions",
                "correlation_length_definition": "first positive radius where the isotropic ensemble covariance equals exp(-1)",
                "shell_root_q": _E_FOLD_SHELL_Q,
                "configured_ensemble_e_fold_length_m": ell,
                "finite_realization_caveat": "The seeded finite set of directions has discrete spectral lines and an anisotropic, oscillatory empirical correlation; correlation_length_m is the target isotropic ensemble e-fold length, not a fitted or guaranteed local realization statistic.",
            },
            "random_modes": [
                {
                    "wavevector_rad_per_m": list(mode.wavevector_rad_per_m),
                    "phase_rad": mode.phase_rad,
                }
                for mode in self.modes
            ],
            "status": FIELD_STATUS,
            "lineage": {
                "source_kind": "synthetic_seeded_field",
                "generator": "finite random Fourier modes with fixed seed",
                "upstream_artifact_sha256": None,
                "material_measurement_source": None,
                "calibration_status": FIELD_STATUS,
                "relicweave_lineage_compatibility": "metadata-only: source_kind, status, and artifact hash slot are exposed; no RelicWeave material record or calibration is claimed",
            },
        }

    def to_manifest(self) -> dict[str, Any]:
        """Return a detached JSON-compatible manifest snapshot."""

        payload = self._manifest_payload()
        canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
        payload["field_id_sha256"] = hashlib.sha256(canonical.encode("utf-8")).hexdigest()
        return payload

    def manifest_json(self, *, indent: int = 2) -> str:
        return json.dumps(self.to_manifest(), sort_keys=True, indent=indent, ensure_ascii=False) + "\n"

