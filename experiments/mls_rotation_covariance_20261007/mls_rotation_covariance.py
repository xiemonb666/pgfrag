"""Coordinate-covariance check for the accepted 27-site MLS provider only."""
from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
CASE = ROOT / "pgfrag/experiments/mixed93_actual_cohesive_bridge_20261006/case.json"
MLS_DIR = ROOT / "pgfrag/experiments/mls3d_kernel_20261005"
PINS = HERE / "source_pins.json"
TOL = 1.0e-10
MAX_QUERIES = 64


def verify_pins() -> None:
    pins = json.loads(PINS.read_text(encoding="utf-8"))
    for item in [pins["cohort_case"], *pins["files"]]:
        digest = hashlib.sha256((ROOT / item["path"]).read_bytes()).hexdigest()
        if digest != item["sha256"]:
            raise RuntimeError(f"source pin mismatch: {item['path']}")


def _provider():
    if str(MLS_DIR) not in sys.path:
        sys.path.insert(0, str(MLS_DIR))
    from mls3d_mls import MLS3DConfig
    from mls3d_shape import shape_functions_3d
    return MLS3DConfig, shape_functions_3d


def _inputs():
    case = json.loads(CASE.read_text(encoding="utf-8"))
    source = np.asarray(case["source_binding"]["source_positions_xyz_m"], dtype=np.float64)
    geo, rule = case["geometry"], case["interface_quadrature"]
    query = np.asarray([
        (geo["radius_m"] / 2.0 * (x - 1.0), 0.0, geo["length_m"] / 2.0 * z)
        for x in rule["gl4_nodes"] for z in rule["gl4_nodes"]
    ], dtype=np.float64)
    mls = case["mls"]
    config_data = {
        "reference_origin": tuple(float(v) for v in mls["reference_origin_m"]),
        "reference_scale": float(mls["reference_scale_m"]),
        "support_radius": float(mls["support_radius_m"]),
        "max_condition_number": float(mls["max_condition_number"]),
        "rank_relative_tolerance": float(mls["rank_relative_tolerance"]),
        "switch_tolerance": float(mls["switch_tolerance"]),
        "duplicate_tolerance": float(mls["duplicate_tolerance"]),
    }
    if source.shape != (27, 3) or query.shape != (16, 3):
        raise RuntimeError("fixture must retain 27 frozen sources and 16 interface queries")
    if case["mls"]["api"] != "mls3d_shape.shape_functions_3d":
        raise RuntimeError("unexpected MLS provider contract")
    return source, query, config_data, json.loads((HERE / "case.json").read_text(encoding="utf-8"))["transforms"]


def _rotation(angles) -> np.ndarray:
    a, b = float(angles["angle_z_rad"]), float(angles["angle_x_rad"])
    rz = np.array([[math.cos(a), -math.sin(a), 0.0],
                   [math.sin(a), math.cos(a), 0.0], [0.0, 0.0, 1.0]])
    rx = np.array([[1.0, 0.0, 0.0],
                   [0.0, math.cos(b), -math.sin(b)],
                   [0.0, math.sin(b), math.cos(b)]])
    return rz @ rx


def _record(result, source, query, grad_to_per_m=1.0):
    return {
        "source_xyz": source.tolist(), "query_xyz": query.tolist(),
        "weights": result.weights.tolist(),
        "gradients_input_unit_inverse": result.gradients.tolist(),
        "gradients_per_m": (result.gradients * grad_to_per_m).tolist(),
        "condition_numbers": result.condition_numbers.tolist(),
        "support_counts": result.support_counts.tolist(),
        "statuses": [getattr(s, "value", str(s)) for s in result.statuses],
    }


def _length_config(c, unit):
    return {
        "coordinate_unit": unit,
        "reference_origin": list(c["reference_origin"]),
        "reference_scale": c["reference_scale"],
        "support_radius": c["support_radius"],
    }


def evaluate():
    """Make exactly four 16-query calls to the real provider and return evidence."""
    verify_pins()
    Config, shape = _provider()
    source, query, c, transforms = _inputs()
    base_config = Config(**c)
    base = shape(source, query, base_config, visibility_mask=None)

    rotation = _rotation(transforms)
    shift = np.asarray(transforms["translation_m"], dtype=np.float64)
    origin = np.asarray(c["reference_origin"], dtype=np.float64)
    source_r, query_r = source @ rotation.T + shift, query @ rotation.T + shift
    c_r = dict(c, reference_origin=tuple(rotation @ origin + shift))
    rotated = shape(source_r, query_r, Config(**c_r), visibility_mask=None)

    unit_factor = float(transforms["uniform_unit_factor"])
    source_mm, query_mm = source * unit_factor, query * unit_factor
    c_mm = dict(c, reference_origin=tuple(unit_factor * origin),
                reference_scale=unit_factor * c["reference_scale"],
                support_radius=unit_factor * c["support_radius"])
    millimetres = shape(source_mm, query_mm, Config(**c_mm), visibility_mask=None)

    stretch = np.diag(transforms["bad_anisotropic_diagonal_scale"])
    source_bad = origin + (source - origin) @ stretch.T
    query_bad = origin + (query - origin) @ stretch.T
    anisotropic = shape(source_bad, query_bad, base_config, visibility_mask=None)
    calls = 4 * len(query)
    if calls > MAX_QUERIES:
        raise RuntimeError("new MLS query budget exceeded")
    grad_rot_expected = np.einsum("ij,qnj->qni", rotation, base.gradients)
    comparisons = {
        "rotation_weights_max_abs": float(np.max(np.abs(rotated.weights - base.weights))),
        "rotation_gradient_normalized_max_abs": float(np.max(np.abs(rotated.gradients - grad_rot_expected)) * c["reference_scale"]),
        "millimetre_weights_max_abs": float(np.max(np.abs(millimetres.weights - base.weights))),
        "millimetre_gradient_normalized_max_abs": float(np.max(np.abs(millimetres.gradients * unit_factor - base.gradients)) * c["reference_scale"]),
        "bad_anisotropic_weights_max_abs": float(np.max(np.abs(anisotropic.weights - base.weights))),
    }
    return {
        "classification": "GIVEN_SOURCE_GRID_AND_INTERFACE_QUADRATURE_ONLY_NOT_ADMITTED",
        "fixture_case_sha256": hashlib.sha256(CASE.read_bytes()).hexdigest(),
        "source_pins_sha256": hashlib.sha256(PINS.read_bytes()).hexdigest(),
        "new_provider_queries": calls,
        "tolerance": TOL,
        "rotation": rotation.tolist(), "translation_m": shift.tolist(),
        "bad_anisotropic_scale": stretch.tolist(), "millimetre_length_factor": unit_factor,
        "length_configs": {
            "identity_m": _length_config(c, "m"),
            "rigid_rotation_translation_m": _length_config(c_r, "m"),
            "uniformly_scaled_mm": _length_config(c_mm, "mm"),
            "bad_anisotropic_negative_m": _length_config(c, "m"),
        },
        "material_parameters_changed": False,
        "comparisons": comparisons,
        "results": {
            "identity_m": _record(base, source, query),
            "rigid_rotation_translation_m": _record(rotated, source_r, query_r),
            "uniformly_scaled_mm": _record(millimetres, source_mm, query_mm, 1000.0),
            "bad_anisotropic_negative_m": _record(anisotropic, source_bad, query_bad),
        },
    }
