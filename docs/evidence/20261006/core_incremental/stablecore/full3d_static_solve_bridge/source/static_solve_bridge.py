"""Read one frozen full-3D operator result and solve its six-moment gauge system.

The module imports no numerical package unless called with ``--execute``.
It never builds quadrature or calls the MLS/provider/operator/field code.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import time
from typing import Any


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CASE_PATH = HERE / "static_solve_case.json"
DEFAULT_OUTPUT = (
    ROOT
    / "pgfrag"
    / "tmp"
    / "review_20261006"
    / "full3d_static_solve_bridge_20261006"
    / "run02"
    / "static_solve_bridge_result.json"
)
CASE_ID = "modeiii_pu_full3d_static_solve_bridge_20261006"
OPERATOR_CASE_ID = "modeiii_pu_full3d_elastic_operator_20261006"
SOURCE_COUNT = 27
DOF_COUNT = 162
BASIS_COUNT = 54
VOLUME_COUNT = 256
BOUNDARY_COUNT = 128
QUERY_COUNT = VOLUME_COUNT + BOUNDARY_COUNT
REPLAY_REGION_SLICES = {
    "volume": [0, VOLUME_COUNT],
    "lateral_side": [VOLUME_COUNT, VOLUME_COUNT + 64],
    "top_cap": [VOLUME_COUNT + 64, VOLUME_COUNT + 96],
    "bottom_cap": [VOLUME_COUNT + 96, QUERY_COUNT],
}
REPLAY_QUERY_ORDER = ["volume", "lateral_side", "top_cap", "bottom_cap"]
REPLAY_RIGID_MODE_ORDER = ["tx", "ty", "tz", "rx", "ry", "rz"]
REPLAY_ARRAY_UNITS = {
    "operator.K_N_per_m": "N/m",
    "operator.boundary_load_N": "N",
    "operator.manufactured_coefficients_m": "m",
    "sources.positions_xyz_m": "m",
    "quadrature.positions_xyz_m": "m",
    "quadrature.volume_weights_m3": "m^3",
    "quadrature.boundary_area_weights_m2": "m^2",
    "basis.values": "dimensionless",
    "basis.gradients_per_m": "1/m",
    "boundary.traction_Pa": "Pa",
    "exact_field.displacement_m": "m",
    "exact_field.displacement_gradient_dimensionless": "m/m (dimensionless)",
}
GAUGE_COUNT = 6
SADDLE_COUNT = DOF_COUNT + GAUGE_COUNT
MAX_WALL_SECONDS = 60.0
MAX_NUMERIC_ARRAY_BYTES = 8 * 1024 * 1024
MAX_RESULT_BYTES = 8 * 1024 * 1024


class BridgeBlocked(RuntimeError):
    """A fail-closed input, support, or spectrum condition."""


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _read_json(path: Path, byte_limit: int) -> tuple[dict[str, Any], bytes]:
    try:
        if path.stat().st_size > byte_limit:
            raise BridgeBlocked(f"input exceeds {byte_limit} bytes")
        data = path.read_bytes()
    except OSError as exc:
        raise BridgeBlocked(f"cannot read input file: {type(exc).__name__}") from exc
    if len(data) > byte_limit:
        raise BridgeBlocked(f"input exceeds {byte_limit} bytes")
    try:
        value = json.loads(data.decode("utf-8"), parse_constant=_reject_json_constant)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise BridgeBlocked(f"invalid UTF-8 JSON input: {type(exc).__name__}") from exc
    if not isinstance(value, dict):
        raise BridgeBlocked("JSON root must be an object")
    return value, data


def _reject_json_constant(value: str) -> None:
    raise ValueError(f"non-finite JSON number {value} is forbidden")


def _case() -> dict[str, Any]:
    value, _ = _read_json(CASE_PATH, MAX_RESULT_BYTES)
    if value.get("case_id") != CASE_ID:
        raise BridgeBlocked("static solve case ID mismatch")
    return value


def _resolve_saved_result(value: str | None, case: dict[str, Any]) -> Path:
    saved = case.get("saved_result")
    if not isinstance(saved, dict):
        raise BridgeBlocked("saved_result pin is missing from the frozen case")
    raw_path = value or saved.get("path")
    if not isinstance(raw_path, str) or not raw_path.strip():
        raise BridgeBlocked("saved-result path is not frozen")
    path = Path(raw_path)
    return path.resolve() if path.is_absolute() else (ROOT / path).resolve()


def _load_pinned_input(
    path: Path, case: dict[str, Any]
) -> tuple[dict[str, Any], bytes]:
    saved = case.get("saved_result")
    expected_hash = saved.get("sha256") if isinstance(saved, dict) else None
    if not isinstance(expected_hash, str) or len(expected_hash) != 64:
        raise BridgeBlocked("saved-result SHA-256 pin is not frozen")
    result, raw = _read_json(path, MAX_RESULT_BYTES)
    actual_hash = _sha256_bytes(raw)
    if actual_hash.lower() != expected_hash.lower():
        raise BridgeBlocked("saved-result SHA-256 does not match the frozen pin")
    return result, raw


def _status_ok(value: Any) -> bool:
    if isinstance(value, str):
        normalized = value.strip().lower()
        return normalized in {"ok", "success"} or normalized.startswith("ok_")
    if isinstance(value, dict):
        return _status_ok(value.get("status"))
    return False


def _fixed_case_matches(root: dict[str, Any], case: dict[str, Any]) -> None:
    fixed = root.get("fixed_case")
    if not isinstance(fixed, dict) or fixed.get("case_id") != OPERATOR_CASE_ID:
        raise BridgeBlocked("saved result does not identify the frozen full-3D operator case")
    geometry = fixed.get("geometry")
    material = fixed.get("material_and_manufactured_field")
    if not isinstance(geometry, dict) or not isinstance(material, dict):
        raise BridgeBlocked("saved result is missing fixed geometry/material metadata")
    expected = case["manufactured_case"]
    if (
        geometry.get("radius_m") != expected["radius_m"]
        or geometry.get("length_m") != expected["length_m"]
        or material.get("young_modulus_Pa") != expected["young_modulus_Pa"]
        or material.get("poisson_ratio") != expected["poisson_ratio"]
        or material.get("K_III_Pa_sqrt_m") != expected["K_III_Pa_sqrt_m"]
    ):
        raise BridgeBlocked("saved result differs from the frozen manufactured benchmark")


def _as_f64(np: Any, value: Any, shape: tuple[int, ...], name: str) -> Any:
    try:
        array = np.asarray(value, dtype=np.float64)
    except (TypeError, ValueError, OverflowError) as exc:
        raise BridgeBlocked(f"{name} is not a numeric array") from exc
    if array.shape != shape or not bool(np.isfinite(array).all()):
        raise BridgeBlocked(f"{name} has the wrong shape or contains non-finite values")
    return np.ascontiguousarray(array)


def _engineering_strain(np: Any, gradient: Any) -> Any:
    result = np.empty((gradient.shape[0], 6), dtype=np.float64)
    result[:, 0] = gradient[:, 0, 0]
    result[:, 1] = gradient[:, 1, 1]
    result[:, 2] = gradient[:, 2, 2]
    result[:, 3] = gradient[:, 0, 1] + gradient[:, 1, 0]
    result[:, 4] = gradient[:, 1, 2] + gradient[:, 2, 1]
    result[:, 5] = gradient[:, 0, 2] + gradient[:, 2, 0]
    return result


def _reconstruct(
    np: Any,
    coefficients: Any,
    keys: list[dict[str, Any]],
    basis_rows: dict[tuple[str, str], int],
    values: Any,
    gradients: Any,
    row_slice: slice,
) -> tuple[Any, Any]:
    row_count = values[row_slice].shape[0]
    displacement = np.zeros((row_count, 3), dtype=np.float64)
    displacement_gradient = np.zeros((row_count, 3, 3), dtype=np.float64)
    component_index = {"ux": 0, "uy": 1, "uz": 2}
    for slot, key in enumerate(keys):
        row = basis_rows[(key["source_id"], key["basis_kind"])]
        component = component_index[key["component"]]
        coefficient = float(coefficients[slot])
        displacement[:, component] += coefficient * values[row_slice, row]
        displacement_gradient[:, component, :] += (
            coefficient * gradients[row_slice, row, :]
        )
    if not bool(np.isfinite(displacement).all()) or not bool(
        np.isfinite(displacement_gradient).all()
    ):
        raise BridgeBlocked("reconstructed displacement/gradient is non-finite")
    return displacement, displacement_gradient


def _rigid_modes_at_volume_centroid(
    np: Any,
    raw_modes: Any,
    saved_center: Any,
    volume_centroid: Any,
    radius: float,
) -> tuple[Any, Any]:
    translations = raw_modes[:, :3]
    result = np.array(raw_modes, dtype=np.float64, copy=True)
    axes = np.eye(3, dtype=np.float64)
    center_shift = saved_center - volume_centroid
    for axis in range(3):
        result[:, 3 + axis] = (
            raw_modes[:, 3 + axis]
            + translations @ np.cross(axes[axis], center_shift)
        ) / radius
    return result, center_shift


def _spectrum_gate(
    np: Any,
    stiffness: Any,
    scale: Any,
    raw_rigid_modes: Any,
    condition_max: float,
) -> dict[str, Any]:
    diagonal = np.diag(stiffness)
    if not bool(np.isfinite(diagonal).all()) or bool(np.any(diagonal <= 0.0)):
        raise BridgeBlocked("raw K has a non-positive diagonal; coefficient scaling is unavailable")
    scale[:] = 1.0 / np.sqrt(diagonal)
    scaled_raw = scale[:, None] * stiffness * scale[None, :]
    scaled_symmetric = 0.5 * (scaled_raw + scaled_raw.T)
    if not bool(np.isfinite(scaled_symmetric).all()):
        raise BridgeBlocked("scaled raw K is non-finite")
    eigenvalues = np.linalg.eigvalsh(scaled_symmetric)
    spectral_scale = max(float(np.max(np.abs(eigenvalues))), np.finfo(np.float64).tiny)
    tolerance = (
        64.0
        * np.finfo(np.float64).eps
        * max(1.0, condition_max)
        * spectral_scale
    )
    near_zero = np.flatnonzero(np.abs(eigenvalues) <= tolerance)
    positive = np.flatnonzero(eigenvalues > tolerance)
    significant_negative = np.flatnonzero(eigenvalues < -tolerance)

    scaled_modes = raw_rigid_modes / scale[:, None]
    rigid_residuals = np.linalg.norm(scaled_symmetric @ scaled_modes, axis=0) / np.maximum(
        spectral_scale * np.linalg.norm(scaled_modes, axis=0),
        np.finfo(np.float64).tiny,
    )
    rigid_tolerance = 128.0 * np.finfo(np.float64).eps * max(1.0, condition_max)
    if (
        near_zero.size != GAUGE_COUNT
        or positive.size != DOF_COUNT - GAUGE_COUNT
        or significant_negative.size
        or not bool(np.isfinite(rigid_residuals).all())
        or bool(np.any(rigid_residuals > rigid_tolerance))
    ):
        raise BridgeBlocked(
            "scaled K does not resolve exactly six expected rigid null modes and 156 positive modes"
        )
    return {
        "eigenvalues_scaled_symmetric_K": eigenvalues,
        "zero_tolerance": tolerance,
        "near_zero_count": int(near_zero.size),
        "positive_count": int(positive.size),
        "significant_negative_count": int(significant_negative.size),
        "strict_negative_count_including_roundoff": int(np.count_nonzero(eigenvalues < 0.0)),
        "rigid_scaled_residuals": rigid_residuals,
        "rigid_scaled_residual_tolerance": rigid_tolerance,
        "scaled_K_raw_relative_asymmetry": float(
            np.linalg.norm(scaled_raw - scaled_raw.T, ord="fro")
            / max(np.linalg.norm(scaled_raw, ord="fro"), np.finfo(np.float64).tiny)
        ),
        "scaled_K": scaled_raw,
    }


def _budget_estimate() -> dict[str, int]:
    # Float64 data read by this bridge plus K, scaled K, 168x168 saddle,
    # solve workspace, gauge/mode arrays, and saved volume field comparisons.
    saved = (
        DOF_COUNT * DOF_COUNT
        + DOF_COUNT
        + DOF_COUNT * GAUGE_COUNT
        + DOF_COUNT * GAUGE_COUNT
        + SOURCE_COUNT * 3
        + QUERY_COUNT * 3
        + VOLUME_COUNT
        + BOUNDARY_COUNT * 3
        + 2 * QUERY_COUNT
        + BOUNDARY_COUNT
        + QUERY_COUNT * BASIS_COUNT
        + QUERY_COUNT * BASIS_COUNT * 3
        + QUERY_COUNT * 3
        + QUERY_COUNT * 9
        + BOUNDARY_COUNT * 3
    ) * 8
    work = (
        3 * DOF_COUNT * DOF_COUNT
        + 3 * SADDLE_COUNT * SADDLE_COUNT
        + 4 * DOF_COUNT * GAUGE_COUNT
        + VOLUME_COUNT * (3 + 9 + 6)
    ) * 8
    return {
        "saved_float64_arrays_bytes": int(saved),
        "solve_and_field_work_arrays_bytes": int(work),
        "estimated_explicit_numeric_array_bytes": int(saved + work),
        "cap_bytes": MAX_NUMERIC_ARRAY_BYTES,
    }


def _array_budget_guard() -> dict[str, int]:
    estimate = _budget_estimate()
    if estimate["estimated_explicit_numeric_array_bytes"] > MAX_NUMERIC_ARRAY_BYTES:
        raise BridgeBlocked("explicit float64 array estimate exceeds 8 MiB")
    return estimate


def _run_numeric(saved_result: dict[str, Any], input_hash: str, case: dict[str, Any]) -> dict[str, Any]:
    import numpy as np

    started = time.monotonic()
    replay = saved_result.get("static_solve_replay")
    if not isinstance(replay, dict) or replay.get("schema") != "static_solve_replay.v1":
        raise BridgeBlocked("saved result lacks the frozen static_solve_replay.v1 data block")
    if (
        saved_result.get("case_id") != OPERATOR_CASE_ID
        or replay.get("case_id") != OPERATOR_CASE_ID
        or not _status_ok(saved_result.get("status"))
        or saved_result.get("assembly_status") != "OK"
        or saved_result.get("mechanical_eligibility") != "DIAGNOSTIC_PASS_UNQUALIFIED"
    ):
        raise BridgeBlocked("pinned operator result status/case does not authorize a diagnostic solve")
    if (
        replay.get("query_order") != REPLAY_QUERY_ORDER
        or replay.get("rigid_mode_order") != REPLAY_RIGID_MODE_ORDER
        or replay.get("array_units") != REPLAY_ARRAY_UNITS
    ):
        raise BridgeBlocked("saved replay query order, rigid-mode order, or unit contract is inconsistent")
    spectrum_input = saved_result.get("spectrum")
    if (
        not isinstance(spectrum_input, dict)
        or spectrum_input.get("gate_passed") is not True
        or spectrum_input.get("observed_near_zero_count") != GAUGE_COUNT
        or spectrum_input.get("observed_positive_count") != DOF_COUNT - GAUGE_COUNT
        or spectrum_input.get("significant_negative_count") != 0
    ):
        raise BridgeBlocked("upstream full-3D spectrum gate did not pass")
    field_reference = saved_result.get("field_reference")
    if not isinstance(field_reference, dict) or field_reference.get("all_statuses_ok") is not True:
        raise BridgeBlocked("upstream exact-field status summary is not all OK")
    _fixed_case_matches(saved_result, case)
    budget = _array_budget_guard()

    operator = replay.get("operator")
    registry = replay.get("registry")
    sources = replay.get("sources")
    quadrature = replay.get("quadrature")
    basis = replay.get("basis")
    exact = replay.get("exact_field")
    boundary = replay.get("boundary")
    if not all(isinstance(item, dict) for item in (operator, registry, sources, quadrature, basis, exact, boundary)):
        raise BridgeBlocked("saved replay block is missing an operator, registry, source, quadrature, basis, field, or boundary object")

    stiffness = _as_f64(np, operator.get("K_N_per_m"), (DOF_COUNT, DOF_COUNT), "K_N_per_m")
    force = _as_f64(np, operator.get("boundary_load_N"), (DOF_COUNT,), "boundary_load_N")
    raw_rigid = _as_f64(
        np,
        operator.get("rigid_mode_coefficients"),
        (DOF_COUNT, GAUGE_COUNT),
        "rigid_mode_coefficients",
    )
    raw_rigid_by_mode = _as_f64(
        np,
        operator.get("rigid_mode_coefficients_by_mode"),
        (GAUGE_COUNT, DOF_COUNT),
        "rigid_mode_coefficients_by_mode",
    )
    if not np.array_equal(raw_rigid_by_mode.T, raw_rigid):
        raise BridgeBlocked("saved rigid-mode coefficient orientations disagree")
    source_positions = _as_f64(
        np, sources.get("positions_xyz_m"), (SOURCE_COUNT, 3), "source positions"
    )
    positions = _as_f64(
        np, quadrature.get("positions_xyz_m"), (QUERY_COUNT, 3), "384 query positions"
    )
    volume_weights = _as_f64(
        np, quadrature.get("volume_weights_m3"), (VOLUME_COUNT,), "volume weights"
    )
    boundary_normals = _as_f64(
        np,
        quadrature.get("boundary_normals"),
        (BOUNDARY_COUNT, 3),
        "boundary normals",
    )
    boundary_weights = _as_f64(
        np,
        quadrature.get("boundary_area_weights_m2"),
        (BOUNDARY_COUNT,),
        "boundary area weights",
    )
    values = _as_f64(
        np, basis.get("values"), (QUERY_COUNT, BASIS_COUNT), "basis values"
    )
    gradients = _as_f64(
        np,
        basis.get("gradients_per_m"),
        (QUERY_COUNT, BASIS_COUNT, 3),
        "basis gradients",
    )
    exact_displacement = _as_f64(
        np,
        exact.get("displacement_m"),
        (QUERY_COUNT, 3),
        "exact displacement",
    )
    exact_gradient = _as_f64(
        np,
        exact.get("displacement_gradient_dimensionless"),
        (QUERY_COUNT, 3, 3),
        "exact displacement gradient",
    )
    traction = _as_f64(
        np,
        boundary.get("traction_Pa"),
        (BOUNDARY_COUNT, 3),
        "full-vector boundary traction",
    )

    if quadrature.get("region_slices_half_open") != REPLAY_REGION_SLICES:
        raise BridgeBlocked("saved query regions are not volume 256, side 64, and two caps of 32")
    if not bool(np.all(volume_weights > 0.0)) or not bool(np.all(boundary_weights > 0.0)):
        raise BridgeBlocked("quadrature weights must be positive")
    if np.unique(source_positions, axis=0).shape[0] != SOURCE_COUNT:
        raise BridgeBlocked("saved source positions are not unique")
    normal_lengths = np.linalg.norm(boundary_normals, axis=1)
    if not np.allclose(normal_lengths, 1.0, rtol=0.0, atol=1e-12):
        raise BridgeBlocked("saved boundary normals are not unit vectors")
    if np.unique(positions, axis=0).shape[0] != QUERY_COUNT:
        raise BridgeBlocked("the saved 384 query locations are not unique")
    volume = float(np.sum(volume_weights))
    if not math.isfinite(volume) or volume <= 0.0:
        raise BridgeBlocked("saved volume quadrature measure is non-positive or non-finite")
    volume_positions = positions[:VOLUME_COUNT]
    boundary_positions = positions[VOLUME_COUNT:]
    volume_centroid = np.einsum("q,qc->c", volume_weights, volume_positions) / volume
    source_centroid = np.mean(source_positions, axis=0)
    rigid_center_value = replay.get("rigid_mode_center_xyz_m")
    if rigid_center_value is None:
        rigid_center_value = sources.get("rigid_mode_center_xyz_m")
    saved_rigid_center = _as_f64(np, rigid_center_value, (3,), "rigid-mode center")
    if not np.allclose(saved_rigid_center, source_centroid, rtol=0.0, atol=1e-12 * 0.01):
        raise BridgeBlocked("saved rigid-mode center differs from the mean source position")
    if sources.get("rigid_mode_rotation_definition") != (
        "omega cross (source_position - source_centroid), with omega in radians and no division by R"
    ):
        raise BridgeBlocked("saved rigid-mode rotation definition is not the pinned producer convention")

    if not isinstance(exact.get("statuses"), list) or exact.get("statuses") != ["OK"] * QUERY_COUNT:
        raise BridgeBlocked("saved exact-field query status rows are not all OK")
    provider_evidence = replay.get("provider_evidence")
    if not isinstance(provider_evidence, dict):
        raise BridgeBlocked("saved replay provider evidence is missing")
    provider_condition_max = 0.0
    for batch_name, expected_count in (("volume", VOLUME_COUNT), ("boundary", BOUNDARY_COUNT)):
        batch = provider_evidence.get(batch_name)
        if (
            not isinstance(batch, dict)
            or batch.get("query_count") != expected_count
            or batch.get("consumed_queries") != expected_count
            or batch.get("mls_statuses") != ["OK"] * expected_count
            or batch.get("field_statuses") != ["OK"] * expected_count
        ):
            raise BridgeBlocked(f"saved {batch_name} provider/field status rows are incomplete or not OK")
        support_counts = _as_f64(
            np, batch.get("support_counts"), (expected_count,), f"{batch_name} support counts"
        )
        condition_numbers = _as_f64(
            np,
            batch.get("condition_numbers"),
            (expected_count,),
            f"{batch_name} condition numbers",
        )
        if (
            not bool(np.all(support_counts > 0.0))
            or not bool(np.all(support_counts == np.floor(support_counts)))
            or not bool(np.all(condition_numbers >= 1.0))
        ):
            raise BridgeBlocked(f"saved {batch_name} MLS support/condition evidence is invalid")
        provider_condition_max = max(provider_condition_max, float(np.max(condition_numbers)))

    source_ids = registry.get("source_ids")
    ordered_ids = registry.get("ordered_canonical_ids")
    raw_keys = registry.get("ordered_keys")
    # The frozen producer serializes basis_row_ids beside registry, directly
    # under static_solve_replay (the ordered DofKey fields remain in registry).
    raw_basis_ids = replay.get("basis_row_ids")
    if (
        not isinstance(source_ids, list)
        or len(source_ids) != SOURCE_COUNT
        or any(type(item) is not str for item in source_ids)
        or len(set(source_ids)) != SOURCE_COUNT
        or not isinstance(ordered_ids, list)
        or len(ordered_ids) != DOF_COUNT
        or any(type(item) is not str for item in ordered_ids)
        or len(set(ordered_ids)) != DOF_COUNT
        or not isinstance(raw_keys, list)
        or len(raw_keys) != DOF_COUNT
        or not isinstance(raw_basis_ids, list)
        or len(raw_basis_ids) != BASIS_COUNT
    ):
        raise BridgeBlocked("saved semantic registry dimensions or uniqueness are invalid")
    keys: list[dict[str, Any]] = []
    for slot, raw_key in enumerate(raw_keys):
        if not isinstance(raw_key, dict):
            raise BridgeBlocked("ordered registry key must be an object")
        component = raw_key.get("component", raw_key.get("field_component"))
        key = {
            "slot": raw_key.get("slot", slot),
            "canonical_id": raw_key.get("canonical_id"),
            "source_id": raw_key.get("source_id"),
            "component": component,
            "basis_kind": raw_key.get("basis_kind"),
            "basis_index": raw_key.get("basis_index"),
        }
        if (
            key["slot"] != slot
            or key["source_id"] not in source_ids
            or key["component"] not in {"ux", "uy", "uz"}
            or key["basis_kind"] not in {"base", "front_f1_unshifted"}
            or key["basis_index"] != (0 if key["basis_kind"] == "base" else 1)
        ):
            raise BridgeBlocked(f"malformed semantic registry entry at slot {slot}")
        canonical = json.dumps(
            [key["source_id"], key["component"], key["basis_kind"], key["basis_index"]],
            ensure_ascii=True,
            separators=(",", ":"),
        )
        if ordered_ids[slot] != canonical:
            raise BridgeBlocked(f"canonical ID and registry key differ at slot {slot}")
        if key["canonical_id"] not in (None, canonical):
            raise BridgeBlocked(f"embedded canonical ID and registry key differ at slot {slot}")
        keys.append(key)

    basis_rows: dict[tuple[str, str], int] = {}
    for row, raw_row in enumerate(raw_basis_ids):
        if isinstance(raw_row, dict):
            source_id, basis_kind = raw_row.get("source_id"), raw_row.get("basis_kind")
        elif isinstance(raw_row, (list, tuple)) and len(raw_row) == 2:
            source_id, basis_kind = raw_row
        else:
            raise BridgeBlocked(f"malformed scalar basis row ID at row {row}")
        pair = (source_id, basis_kind)
        if source_id not in source_ids or basis_kind not in {"base", "front_f1_unshifted"} or pair in basis_rows:
            raise BridgeBlocked(f"duplicate or invalid scalar basis row ID at row {row}")
        basis_rows[pair] = row
    if set(basis_rows) != {
        (source_id, basis_kind)
        for source_id in source_ids
        for basis_kind in ("base", "front_f1_unshifted")
    }:
        raise BridgeBlocked("scalar basis rows do not cover BASE/F1 for all saved sources")

    component_index = {"ux": 0, "uy": 1, "uz": 2}
    source_index = {source_id: row for row, source_id in enumerate(source_ids)}
    expected_rigid = np.zeros((DOF_COUNT, GAUGE_COUNT), dtype=np.float64)
    axes = np.eye(3, dtype=np.float64)
    for slot, key in enumerate(keys):
        if key["basis_kind"] != "base":
            continue
        component = component_index[key["component"]]
        expected_rigid[slot, :3] = axes[:, component]
        source_offset = source_positions[source_index[key["source_id"]]] - saved_rigid_center
        expected_rigid[slot, 3:] = np.cross(axes, source_offset)[:, component]
    if not np.allclose(raw_rigid, expected_rigid, rtol=0.0, atol=1e-13):
        raise BridgeBlocked("saved rigid coefficients are not the six BASE-only physical modes")

    rebuilt_force = np.zeros(DOF_COUNT, dtype=np.float64)
    for slot, key in enumerate(keys):
        basis_row = basis_rows[(key["source_id"], key["basis_kind"])]
        component = component_index[key["component"]]
        rebuilt_force[slot] = np.dot(
            boundary_weights * traction[:, component],
            values[VOLUME_COUNT:, basis_row],
        )
    load_reconstruction_error = rebuilt_force - force
    load_reconstruction_relative = float(
        np.linalg.norm(load_reconstruction_error)
        / max(float(np.linalg.norm(force)), np.finfo(np.float64).tiny)
    )
    if not math.isfinite(load_reconstruction_relative) or load_reconstruction_relative > 1e-10:
        raise BridgeBlocked("saved full-vector load does not match its boundary rows and shape values")

    material = case["manufactured_case"]
    radius = float(material["radius_m"])
    young = float(material["young_modulus_Pa"])
    poisson = float(material["poisson_ratio"])
    k_iii = float(material["K_III_Pa_sqrt_m"])
    shear = young / (2.0 * (1.0 + poisson))
    amplitude = k_iii / shear * math.sqrt(2.0 / math.pi)
    exact_coefficients = np.zeros(DOF_COUNT, dtype=np.float64)
    for slot, key in enumerate(keys):
        if key["basis_kind"] == "front_f1_unshifted" and key["component"] == "uz":
            exact_coefficients[slot] = amplitude * math.sqrt(radius)
    if np.count_nonzero(exact_coefficients) != SOURCE_COUNT:
        raise BridgeBlocked("registry does not contain the 27 exact UZ FRONT_F1 coefficients")
    saved_exact_coefficients = _as_f64(
        np,
        operator.get("manufactured_coefficients_m"),
        (DOF_COUNT,),
        "saved manufactured coefficients",
    )
    exact_coefficient_scale = max(
        float(np.max(np.abs(exact_coefficients))), np.finfo(np.float64).tiny
    )
    if not np.allclose(
        saved_exact_coefficients,
        exact_coefficients,
        rtol=0.0,
        atol=8.0 * np.finfo(np.float64).eps * exact_coefficient_scale,
    ):
        raise BridgeBlocked("saved manufactured coefficient vector differs from the registry reconstruction")

    reconstructed_exact, reconstructed_exact_gradient = _reconstruct(
        np,
        exact_coefficients,
        keys,
        basis_rows,
        values,
        gradients,
        slice(0, VOLUME_COUNT),
    )
    exact_displacement_star_error = reconstructed_exact - exact_displacement[:VOLUME_COUNT]
    exact_gradient_star_error = (
        reconstructed_exact_gradient - exact_gradient[:VOLUME_COUNT]
    )

    gauge_modes = np.zeros((GAUGE_COUNT, VOLUME_COUNT, 3), dtype=np.float64)
    gauge_modes[0, :, 0] = 1.0
    gauge_modes[1, :, 1] = 1.0
    gauge_modes[2, :, 2] = 1.0
    axes = np.eye(3, dtype=np.float64)
    centered_volume = volume_positions - volume_centroid
    for axis in range(3):
        gauge_modes[3 + axis] = np.cross(axes[axis], centered_volume) / radius

    gauge_matrix = np.zeros((GAUGE_COUNT, DOF_COUNT), dtype=np.float64)
    for slot, key in enumerate(keys):
        basis_row = basis_rows[(key["source_id"], key["basis_kind"])]
        component = component_index[key["component"]]
        gauge_matrix[:, slot] = (
            (gauge_modes[:, :, component] * volume_weights[None, :])
            @ values[:VOLUME_COUNT, basis_row]
        ) / volume
    gauge_target = np.einsum(
        "q,lqc,qc->l",
        volume_weights,
        gauge_modes,
        exact_displacement[:VOLUME_COUNT],
    ) / volume
    if not bool(np.isfinite(gauge_matrix).all()) or not bool(np.isfinite(gauge_target).all()):
        raise BridgeBlocked("physical volume gauge contains non-finite entries")
    if not gauge_target[3] > 0.0:
        raise BridgeBlocked("expected positive r_rx manufactured gauge target is missing")

    physical_rigid_modes, rigid_center_shift = _rigid_modes_at_volume_centroid(
        np, raw_rigid, saved_rigid_center, volume_centroid, radius
    )
    gauge_on_rigid = gauge_matrix @ physical_rigid_modes
    singular_values = np.linalg.svd(gauge_on_rigid, compute_uv=False)
    rank_tolerance = 64.0 * np.finfo(np.float64).eps * float(singular_values[0])
    if not bool(np.isfinite(singular_values).all()) or singular_values[-1] <= rank_tolerance:
        raise BridgeBlocked("the six volume moments do not control all six physical rigid modes")

    mls = saved_result.get("mls")
    if (
        not isinstance(mls, dict)
        or mls.get("batch_count") != 2
        or mls.get("batch_query_counts") != [VOLUME_COUNT, BOUNDARY_COUNT]
        or mls.get("provider_queries_consumed_total") != QUERY_COUNT
        or mls.get("status_counts") != {"OK": QUERY_COUNT}
    ):
        raise BridgeBlocked("saved MLS batch/status summary does not match the fixed 384-row contract")
    condition_max = float(mls.get("condition_number_max", 0.0))
    if (
        not math.isfinite(condition_max)
        or condition_max < 1.0
        or not math.isclose(condition_max, provider_condition_max, rel_tol=1e-12, abs_tol=1e-12)
    ):
        raise BridgeBlocked("saved maximum MLS condition number disagrees with provider evidence")
    scale = np.empty(DOF_COUNT, dtype=np.float64)
    spectrum = _spectrum_gate(np, stiffness, scale, raw_rigid, condition_max)

    constraint_norm = np.linalg.norm(gauge_matrix * scale[None, :], axis=1)
    if not bool(np.isfinite(constraint_norm).all()) or bool(np.any(constraint_norm <= 0.0)):
        raise BridgeBlocked("a scaled physical moment constraint has zero or invalid norm")
    constraint_scale = 1.0 / constraint_norm
    scaled_constraint = constraint_scale[:, None] * gauge_matrix * scale[None, :]
    scaled_stiffness = spectrum["scaled_K"]
    saddle = np.zeros((SADDLE_COUNT, SADDLE_COUNT), dtype=np.float64)
    saddle[:DOF_COUNT, :DOF_COUNT] = scaled_stiffness
    saddle[:DOF_COUNT, DOF_COUNT:] = scaled_constraint.T
    saddle[DOF_COUNT:, :DOF_COUNT] = scaled_constraint
    rhs = np.concatenate((scale * force, constraint_scale * gauge_target))
    if not bool(np.isfinite(saddle).all()) or not bool(np.isfinite(rhs).all()):
        raise BridgeBlocked("scaled direct saddle system is non-finite")
    if time.monotonic() - started > MAX_WALL_SECONDS:
        raise BridgeBlocked("60-second wall limit reached before the single direct solve")
    try:
        scaled_solution = np.linalg.solve(saddle, rhs)
    except np.linalg.LinAlgError as exc:
        raise BridgeBlocked(f"direct 168x168 saddle solve failed: {type(exc).__name__}") from exc
    if not bool(np.isfinite(scaled_solution).all()):
        raise BridgeBlocked("direct saddle solve returned non-finite unknowns")
    coefficients = scale * scaled_solution[:DOF_COUNT]
    multipliers = constraint_scale * scaled_solution[DOF_COUNT:]

    equilibrium_residual = stiffness @ coefficients + gauge_matrix.T @ multipliers - force
    constraint_residual = gauge_matrix @ coefficients - gauge_target
    scaled_equilibrium = scale * equilibrium_residual
    scaled_kd = scale * (stiffness @ coefficients)
    scaled_clambda = scale * (gauge_matrix.T @ multipliers)
    scaled_force = scale * force
    weak_denominator = (
        np.linalg.norm(scaled_kd)
        + np.linalg.norm(scaled_clambda)
        + np.linalg.norm(scaled_force)
    )
    weak_residual_relative = float(
        np.linalg.norm(scaled_equilibrium)
        / max(float(weak_denominator), np.finfo(np.float64).tiny)
    )
    constraint_relative = float(
        np.linalg.norm(constraint_residual)
        / max(
            float(np.linalg.norm(gauge_target)),
            float(np.linalg.norm(gauge_matrix @ coefficients)),
            np.finfo(np.float64).tiny,
        )
    )
    strain_energy = 0.5 * float(coefficients @ stiffness @ coefficients)
    external_work = float(coefficients @ force)
    multiplier_work = float(gauge_target @ multipliers)
    energy_identity_defect = 2.0 * strain_energy - external_work + multiplier_work
    saddle_backward_residual = float(
        np.linalg.norm(saddle @ scaled_solution - rhs)
        / max(float(np.linalg.norm(rhs)), np.finfo(np.float64).tiny)
    )
    scalars = (
        weak_residual_relative,
        constraint_relative,
        strain_energy,
        external_work,
        multiplier_work,
        energy_identity_defect,
        saddle_backward_residual,
    )
    if not all(math.isfinite(item) for item in scalars):
        raise BridgeBlocked("solution energy or residual evidence is non-finite")

    solved_displacement, solved_gradient = _reconstruct(
        np,
        coefficients,
        keys,
        basis_rows,
        values,
        gradients,
        slice(0, QUERY_COUNT),
    )
    solved_volume_strain = _engineering_strain(np, solved_gradient[:VOLUME_COUNT])
    exact_volume_strain = _engineering_strain(np, exact_gradient[:VOLUME_COUNT])
    displacement_error = solved_displacement[:VOLUME_COUNT] - exact_displacement[:VOLUME_COUNT]
    strain_error = solved_volume_strain - exact_volume_strain
    displacement_scale = max(
        float(np.max(np.abs(exact_displacement[:VOLUME_COUNT]))), np.finfo(np.float64).tiny
    )
    strain_scale = max(float(np.max(np.abs(exact_volume_strain))), np.finfo(np.float64).tiny)

    boundary_forces = np.einsum("q,qc->c", boundary_weights, traction)
    boundary_moments = np.einsum(
        "q,qc->c", boundary_weights, np.cross(boundary_positions, traction)
    )
    center_moment = np.einsum(
        "q,qc->c",
        boundary_weights,
        np.cross(boundary_positions - saved_rigid_center, traction),
    )
    rigid_work = raw_rigid.T @ force
    wrench_by_region = {}
    for name, section in (
        ("side", slice(0, 64)),
        ("top", slice(64, 96)),
        ("bottom", slice(96, 128)),
    ):
        wrench_by_region[name] = {
            "force_N": np.einsum(
                "q,qc->c", boundary_weights[section], traction[section]
            ),
            "moment_about_origin_Nm": np.einsum(
                "q,qc->c",
                boundary_weights[section],
                np.cross(boundary_positions[section], traction[section]),
            ),
        }
    exact_boundary_work = float(
        np.einsum(
            "qc,qc,q->",
            exact_displacement[VOLUME_COUNT:],
            traction,
            boundary_weights,
        )
    )
    reconstructed_boundary_work = float(
        np.einsum(
            "qc,qc,q->",
            solved_displacement[VOLUME_COUNT:],
            traction,
            boundary_weights,
        )
    )
    top_cap_traction_max_abs = np.max(np.abs(traction[64:96]), axis=0)
    bottom_cap_traction_max_abs = np.max(np.abs(traction[96:128]), axis=0)
    if not (
        np.max(top_cap_traction_max_abs[:2]) > 0.0
        and np.max(bottom_cap_traction_max_abs[:2]) > 0.0
    ):
        raise BridgeBlocked("saved cap traction rows omit the nonzero full-vector in-plane traction")

    lame_lambda = young * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson))
    constitutive = np.zeros((6, 6), dtype=np.float64)
    constitutive[:3, :3] = lame_lambda
    constitutive[0, 0] += 2.0 * shear
    constitutive[1, 1] += 2.0 * shear
    constitutive[2, 2] += 2.0 * shear
    constitutive[3, 3] = shear
    constitutive[4, 4] = shear
    constitutive[5, 5] = shear
    exact_saved_quadrature_energy = 0.5 * float(
        np.dot(
            volume_weights,
            np.einsum(
                "qi,ij,qj->q",
                exact_volume_strain,
                constitutive,
                exact_volume_strain,
            ),
        )
    )
    solved_saved_quadrature_energy = 0.5 * float(
        np.dot(
            volume_weights,
            np.einsum(
                "qi,ij,qj->q",
                solved_volume_strain,
                constitutive,
                solved_volume_strain,
            ),
        )
    )
    exact_coefficient_energy = 0.5 * float(
        exact_coefficients @ stiffness @ exact_coefficients
    )
    exact_coefficient_weak_defect = stiffness @ exact_coefficients - force
    exact_coefficient_gauge_residual = gauge_matrix @ exact_coefficients - gauge_target

    if time.monotonic() - started > MAX_WALL_SECONDS:
        raise BridgeBlocked("60-second single-attempt wall limit exceeded")
    result = {
        "schema": "modeiii_pu_full3d_static_solve_bridge_result.v1",
        "case_id": CASE_ID,
        "status": "SOLVED_DIAGNOSTIC_ONLY_UNQUALIFIED",
        "classification": "SYNTHETIC_MANUFACTURED_BENCHMARK_ONLY_NOT_CERAMIC",
        "source_result": {
            "sha256": input_hash,
            "case_id": OPERATOR_CASE_ID,
            "schema": saved_result.get("schema"),
            "path": case["saved_result"].get("path"),
        },
        "dimensions": {
            "source_count": SOURCE_COUNT,
            "dof_count": DOF_COUNT,
            "saved_unique_location_count": QUERY_COUNT,
            "provider_batch_sizes": [VOLUME_COUNT, BOUNDARY_COUNT],
            "additional_queries": 0,
            "saddle_dimension": SADDLE_COUNT,
            "basis_value_rows": [QUERY_COUNT, BASIS_COUNT],
            "basis_gradient_rows": [QUERY_COUNT, BASIS_COUNT, 3],
        },
        "units": {
            "coordinates": "m",
            "basis_values": "1",
            "basis_gradients": "1/m",
            "K": "N/m",
            "C": "1",
            "g_and_d": "m",
            "f_and_lambda": "N",
            "energy_and_work": "J",
            "boundary_traction": "Pa",
        },
        "gauge": {
            "definition": "six physical-volume moments from saved 256-point volume quadrature; translations and rotations omega x (x-centroid)/R",
            "volume_measure_m3": volume,
            "volume_centroid_xyz_m": volume_centroid,
            "saved_rigid_center_xyz_m": saved_rigid_center,
            "rigid_center_shift_to_volume_centroid_m": rigid_center_shift,
            "ordered_canonical_ids": ordered_ids,
            "C_dimensionless": gauge_matrix,
            "g_target_m": gauge_target,
            "expected_nonzero_target_r_rx_positive_m": float(gauge_target[3]),
            "gauge_on_normalized_rigid_modes": gauge_on_rigid,
            "gauge_on_rigid_singular_values": singular_values,
            "gauge_rank_tolerance": rank_tolerance,
            "coefficient_pinning": False,
        },
        "spectrum": {
            key: value
            for key, value in spectrum.items()
            if key != "scaled_K"
        },
        "solve": {
            "method": "direct np.linalg.solve on the 168x168 row/diagonally-scaled KKT equations; raw K/C/f/g retained",
            "diagonal_scale_S": scale,
            "constraint_row_scale": constraint_scale,
            "raw_applied_load_N": force,
            "reconstructed_full_vector_load_relative_error": load_reconstruction_relative,
            "reconstructed_full_vector_load_raw_error_N": load_reconstruction_error,
            "raw_K_sha256_float64_little_endian_C_order": _sha256_bytes(
                np.asarray(stiffness, dtype="<f8", order="C").tobytes(order="C")
            ),
            "raw_stiffness_relative_asymmetry": float(
                np.linalg.norm(stiffness - stiffness.T, ord="fro")
                / max(np.linalg.norm(stiffness, ord="fro"), np.finfo(np.float64).tiny)
            ),
            "coefficients_m": coefficients,
            "lagrange_multipliers_N": multipliers,
            "weak_equilibrium_residual_N": equilibrium_residual,
            "weak_equilibrium_residual_relative": weak_residual_relative,
            "physical_moment_constraint_residual_m": constraint_residual,
            "physical_moment_constraint_residual_relative": constraint_relative,
            "saddle_backward_residual_relative": saddle_backward_residual,
            "strain_energy_J": strain_energy,
            "saved_gradient_reconstructed_strain_energy_J": solved_saved_quadrature_energy,
            "external_work_J": external_work,
            "g_dot_lambda_J": multiplier_work,
            "energy_identity_expected": "2U-W=-g^T*lambda",
            "energy_identity_defect_2U_minus_W_plus_g_dot_lambda_J": energy_identity_defect,
            "pseudoinverse_or_ridge_used": False,
            "load_projection_or_equilibrium_repair_used": False,
            "eigenvalue_truncation_used": False,
        },
        "rigid_virtual_work_and_boundary_wrench": {
            "rigid_mode_order": ["tx", "ty", "tz", "rx", "ry", "rz"],
            "rigid_work_f_dot_r_N_or_Nm": rigid_work,
            "integrated_force_N": boundary_forces,
            "integrated_moment_about_origin_Nm": boundary_moments,
            "integrated_moment_about_rigid_center_Nm": center_moment,
            "rigid_work_minus_integrated_force_and_center_moment": rigid_work
            - np.concatenate((boundary_forces, center_moment)),
            "regions": wrench_by_region,
            "top_cap_max_abs_traction_xyz_Pa": top_cap_traction_max_abs,
            "bottom_cap_max_abs_traction_xyz_Pa": bottom_cap_traction_max_abs,
            "manufactured_exact_boundary_work_J": exact_boundary_work,
            "solved_reconstructed_boundary_work_J": reconstructed_boundary_work,
        },
        "manufactured_reference_separate_from_solution": {
            "exact_coefficients_m": exact_coefficients,
            "exact_coefficient_gauge_residual_m": exact_coefficient_gauge_residual,
            "exact_coefficient_weak_defect_N": exact_coefficient_weak_defect,
            "exact_coefficient_quadratic_energy_J": exact_coefficient_energy,
            "coefficient_reconstruction_max_displacement_error_m": np.max(
                np.abs(exact_displacement_star_error), axis=0
            ),
            "coefficient_reconstruction_max_gradient_error_dimensionless": np.max(
                np.abs(exact_gradient_star_error), axis=(0, 2)
            ),
            "coefficient_reconstruction_rows": "saved 256 volume points only",
            "saved_exact_volume_strain_energy_by_quadrature_J": exact_saved_quadrature_energy,
            "analytic_total_strain_energy_J": float(
                k_iii**2 * radius * float(material["length_m"]) / (2.0 * shear)
            ),
        },
        "solved_field_at_saved_volume_rows_only": {
            "solved_displacement_m": solved_displacement[:VOLUME_COUNT],
            "exact_displacement_m": exact_displacement[:VOLUME_COUNT],
            "solved_engineering_strain_dimensionless": solved_volume_strain,
            "exact_engineering_strain_dimensionless": exact_volume_strain,
            "engineering_strain_order": [
                "eps_xx",
                "eps_yy",
                "eps_zz",
                "gamma_xy",
                "gamma_yz",
                "gamma_xz",
            ],
            "max_relative_displacement_error": float(
                np.max(np.abs(displacement_error)) / displacement_scale
            ),
            "max_relative_strain_error": float(
                np.max(np.abs(strain_error)) / strain_scale
            ),
            "offgrid_queries_added": 0,
        },
        "array_budget": budget,
        "elapsed_seconds": time.monotonic() - started,
        "limitations": [
            "one fixed quadrature rule is not a convergence study",
            "benchmark values are synthetic and do not describe ceramic material",
            "the global unshifted F1 manufactured bridge is not a local or moving-front crack method",
            "the exact strain is singular at the mathematical front; positive saved quadrature nodes do not create a radial hole",
            "a successful algebraic solve does not qualify fracture, crack propagation, or a physical field",
        ],
    }
    return result


def _json_bytes(result: dict[str, Any]) -> bytes:
    def safe(value: Any) -> Any:
        if isinstance(value, dict):
            return {str(key): safe(item) for key, item in value.items()}
        if isinstance(value, (list, tuple)):
            return [safe(item) for item in value]
        if hasattr(value, "tolist"):
            return safe(value.tolist())
        if hasattr(value, "item"):
            return safe(value.item())
        return value

    text = json.dumps(safe(result), ensure_ascii=False, separators=(",", ":"), allow_nan=False)
    encoded = (text + "\n").encode("utf-8")
    if len(encoded) > MAX_RESULT_BYTES:
        raise BridgeBlocked("result JSON exceeds the 8 MiB cap")
    return encoded


def _worker(saved_result_path: Path, output_path: Path) -> int:
    case = _case()
    source, raw = _load_pinned_input(saved_result_path, case)
    source_hash = _sha256_bytes(raw)
    result = _run_numeric(source, source_hash, case)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_bytes(_json_bytes(result))
    return 0


def _parent(saved_result_path: Path, output_path: Path, case: dict[str, Any]) -> int:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    staging = output_path.with_name(output_path.name + ".staging")
    receipt_path = output_path.parent / "run_receipt.json"
    stdout_path = output_path.parent / "child_stdout.bin"
    stderr_path = output_path.parent / "child_stderr.bin"
    if any(path.exists() for path in (output_path, staging, receipt_path, stdout_path, stderr_path)):
        raise BridgeBlocked("prior run evidence exists; refusing a retry or overwrite")
    input_result, input_raw = _load_pinned_input(saved_result_path, case)
    input_hash = _sha256_bytes(input_raw)
    del input_result
    del input_raw
    if not isinstance(input_hash, str) or len(input_hash) != 64:
        raise BridgeBlocked("pinned input digest did not produce a SHA-256 string")
    command = [
        sys.executable,
        "-B",
        str(Path(__file__).resolve()),
        "--execute",
        "--worker",
        "--saved-result",
        str(saved_result_path),
        "--output",
        str(staging),
    ]
    started = time.monotonic()
    child_return_code: int | None = None
    stdout_data = b""
    stderr_data = b""
    spawn_error: str | None = None
    try:
        completed = subprocess.run(
            command,
            cwd=str(ROOT),
            capture_output=True,
            timeout=MAX_WALL_SECONDS,
            check=False,
        )
        child_return_code, stdout_data, stderr_data = (
            completed.returncode,
            completed.stdout,
            completed.stderr,
        )
    except subprocess.TimeoutExpired as exc:
        stdout_data = exc.stdout or b""
        stderr_data = exc.stderr or b""
    except OSError as exc:
        spawn_error = type(exc).__name__
    elapsed = time.monotonic() - started
    stdout_path.write_bytes(stdout_data)
    stderr_path.write_bytes(stderr_data)
    stage_size = staging.stat().st_size if staging.is_file() else None
    stage_data = staging.read_bytes() if stage_size is not None and 0 < stage_size <= MAX_RESULT_BYTES else None
    stage_hash = _sha256_bytes(stage_data) if stage_data is not None else None
    status = "CHILD_OK" if child_return_code == 0 and stage_data is not None else (
        "CHILD_START_FAILED" if spawn_error is not None else (
            "CHILD_TIMEOUT" if child_return_code is None else "CHILD_FAILED"
        )
    )
    primary_hash = None
    if status == "CHILD_OK":
        output_path.write_bytes(stage_data)
        primary_hash = stage_hash
        staging.unlink()
    receipt = {
        "schema": "full3d_static_solve_bridge_run_receipt.v1",
        "case_id": CASE_ID,
        "status": status,
        "child_return_code": child_return_code,
        "child_start_error": spawn_error,
        "child_elapsed_seconds": elapsed,
        "saved_input_sha256": input_hash,
        "stdout_path": stdout_path.name,
        "stdout_sha256": _sha256_bytes(stdout_data),
        "stderr_path": stderr_path.name,
        "stderr_sha256": _sha256_bytes(stderr_data),
        "stage_result_path": staging.name if stage_size is not None else None,
        "stage_result_sha256": stage_hash,
        "primary_result_path": output_path.name if primary_hash is not None else None,
        "primary_result_sha256": primary_hash,
    }
    receipt_path.write_bytes(_json_bytes(receipt))
    if status != "CHILD_OK":
        raise BridgeBlocked(f"static-solve child ended with {status} (rc={child_return_code})")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execute", action="store_true", help="run the one pinned saved-result solve")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--saved-result", type=str, default=None)
    parser.add_argument("--output", type=str, default=str(DEFAULT_OUTPUT))
    args = parser.parse_args()
    if not args.execute:
        print(
            json.dumps(
                {
                    "case_id": CASE_ID,
                    "status": "SOURCE_ONLY_NOT_RUN",
                    "execution_invoked": False,
                    "numpy_imported": False,
                },
                separators=(",", ":"),
            )
        )
        return 0
    case = _case()
    saved_path = _resolve_saved_result(args.saved_result, case)
    output_path = Path(args.output).resolve()
    try:
        if args.worker:
            return _worker(saved_path, output_path)
        return _parent(saved_path, output_path, case)
    except BridgeBlocked as exc:
        print(
            json.dumps(
                {
                    "case_id": CASE_ID,
                    "status": "BLOCKED",
                    "execution_invoked": True,
                    "reason": str(exc),
                },
                ensure_ascii=False,
                separators=(",", ":"),
            ),
            file=sys.stderr,
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
