"""Read a frozen solved Mode III coefficient vector into the pinned provider and J integrator.

This source-only adapter never assembles stiffness or runs the static solver.
Execution is gated by --execute-readback and is intended to be approved once
after review.  The worker consumes one fixed result hash, evaluates one bounded
provider batch, and writes direct field samples plus status-bearing J results.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback
from typing import Any


HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[2]
PGFRAG = PROJECT_ROOT / "pgfrag"
SOLVE_DIR = PGFRAG / "experiments" / "modeiii_pu_static_solve_v3_20261006"
J_DIR = PGFRAG / "experiments" / "modeiii_j_integral_20261006"
REVIEW_ROOT = PGFRAG / "tmp" / "review_20261006"
SOLVE_REVIEW_DIR = REVIEW_ROOT / "modeiii_pu_static_solve_v3_20261006"
HERE_REVIEW_DIR = REVIEW_ROOT / "solved_field_j_coupling_20261006"
SOURCE_FREEZE_PATH = HERE_REVIEW_DIR / "source_freeze.json"
SOURCE_FREEZE_SHA_PATH = HERE_REVIEW_DIR / "source_freeze.sha256"
INPUT_RESULT = (
    SOLVE_REVIEW_DIR
    / "runs"
    / "run01"
    / "static_solve_diagnostic.json"
)
DEFAULT_OUTPUT = HERE_REVIEW_DIR / "runs" / "run01" / "solved_field_j_coupling.json"

CASE_ID = "modeiii_solved_field_j_readback_coupling_20261006_v1"
RESULT_SCHEMA = "modeiii_solved_field_j_coupling_result.v1"
SOLVE_RESULT_SCHEMA = "modeiii_pu_antiplane_static_solve_result.v3"
SOLVE_CASE_ID = "modeiii_pu_uz_static_traction_solve_20261006_v3"
SOLVE_RESULT_SHA256 = "39cb2dabe767998baf5e84b038f39510d6ef7db10bcd6c77298b480045cb24d3"
SOLVE_SOURCE_SHA256 = "7a8ac7076a788a29d7f1bcae34bcede80b1481030910c37da1af673a2a370a9f"
SOLVE_FREEZE_SHA256 = "1e57fe4e81fca80789550aec376934f3ac450b819ce7963744d30009c901c5a0"
J_INTEGRATOR_NAME = "j_integral_precopy_guard_v3.py"
J_INTEGRATOR_SHA256 = "b4d5ca3734d0e44e7b8ed4bb690257c61e72c647d499003faae6d8fee6880b7b"

MAX_PROVIDER_QUERIES = 256
CONTOUR_QUERY_COUNT = 224
CONTROL_QUERY_COUNT = 24
PLANNED_PROVIDER_QUERIES = CONTOUR_QUERY_COUNT + CONTROL_QUERY_COUNT
MAX_RESULT_JSON_BYTES = 8 * 1024 * 1024
WALL_TIMEOUT_SECONDS = 60.0
WORKER_TIMEOUT_SECONDS = 59.0

EXPECTED_DEPENDENCIES = {
    "pgfrag/experiments/modeiii_pu_static_solve_v3_20261006/static_solve.py": SOLVE_SOURCE_SHA256,
    "pgfrag/experiments/modeiii_pu_static_solve_v3_20261006/static_solve_case.json": "13ca46e05b520624d903be5d71c55e58b4d23db621423b47da4b2a0c174a38a9",
    "pgfrag/experiments/modeiii_pu_operator_20261005/fixed_case.json": "dce0a1f3009246a52198bc0ed833ed67998144549a653e7bce9bcc17babd71b9",
    "pgfrag/experiments/modeiii_pu_operator_20261005/mls_modeiii_provider.py": "6025642fdd119a95c76d693d3edbc4785d982731560e556ce7b13bad109c6034",
    "pgfrag/experiments/enriched_dof_contract_20261005/enriched_dof_contract.py": "e3ddb3b067465bb5b269576a1b8a5100e2bb92c9c924bc09c949702c19efa9d8",
    "pgfrag/experiments/modeiii_field_reference_20261005/field_reference.py": "7b07dfb30a64c50734dc18de2778b7725320de6d5a8f275c55f301c453ef18d7",
    "pgfrag/experiments/mls3d_kernel_20261005/mls3d_shape.py": "66e680616a988182875d17aeb84fadb541867aee13fc15e0ae098f7c98839ef5",
    "pgfrag/experiments/mls3d_kernel_20261005/mls3d_mls.py": "404256d9b7546361a5704c4c3cec7c022069a9cf83a137c58a16e7a64f79ec63",
    "pgfrag/experiments/efgm_elastic_operator_20261005/efgm_elastic_operator.py": "41a5ac6ef6359ed1640ee681c94a3496c18652b29cc72c679b3e37c20dd9ac0d",
    "pgfrag/experiments/modeiii_j_integral_20261006/j_integral_precopy_guard_v3.py": J_INTEGRATOR_SHA256,
    "pgfrag/tmp/review_20261006/modeiii_pu_static_solve_v3_20261006/source_freeze.json": SOLVE_FREEZE_SHA256,
    "pgfrag/tmp/review_20261006/modeiii_pu_static_solve_v3_20261006/runs/run01/static_solve_diagnostic.json": SOLVE_RESULT_SHA256,
}


class CouplingUnknown(RuntimeError):
    """Raised when replayable readback evidence cannot be constructed."""


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _assert_source_pins() -> dict[str, str]:
    observed: dict[str, str] = {}
    mismatches: dict[str, dict[str, str]] = {}
    for relative, expected in EXPECTED_DEPENDENCIES.items():
        path = PROJECT_ROOT / Path(relative)
        actual = _sha256_path(path) if path.is_file() else "MISSING"
        observed[relative] = actual
        if actual != expected:
            mismatches[relative] = {"expected": expected, "actual": actual}

    if not SOURCE_FREEZE_PATH.is_file() or not SOURCE_FREEZE_SHA_PATH.is_file():
        mismatches["coupling_source_freeze"] = {
            "expected": "source_freeze.json plus .sha256 sidecar",
            "actual": "MISSING",
        }
    else:
        freeze_hash = _sha256_path(SOURCE_FREEZE_PATH)
        sidecar = SOURCE_FREEZE_SHA_PATH.read_text(encoding="ascii").strip().split()[0].lower()
        if freeze_hash != sidecar:
            mismatches["coupling_source_freeze"] = {
                "expected": sidecar,
                "actual": freeze_hash,
            }
        observed["pgfrag/tmp/review_20261006/solved_field_j_coupling_20261006/source_freeze.json"] = freeze_hash

    own_hash = _sha256_path(Path(__file__).resolve())
    observed["pgfrag/experiments/solved_field_j_coupling_20261006/read_solved_field_j.py"] = own_hash
    try:
        freeze = json.loads(SOURCE_FREEZE_PATH.read_text(encoding="utf-8"))
        frozen_own_hash = freeze["owned_sources"][
            "pgfrag/experiments/solved_field_j_coupling_20261006/read_solved_field_j.py"
        ].lower()
        if own_hash != frozen_own_hash:
            mismatches["read_solved_field_j.py"] = {
                "expected": frozen_own_hash,
                "actual": own_hash,
            }
    except (OSError, ValueError, KeyError, TypeError) as exc:
        mismatches["coupling_source_freeze_content"] = {
            "expected": "valid owned_sources pin",
            "actual": f"{type(exc).__name__}: {exc}",
        }

    if mismatches:
        raise CouplingUnknown(f"source_pin_mismatch:{json.dumps(mismatches, sort_keys=True)}")
    return observed


def _load_static_module() -> Any:
    source_path = SOLVE_DIR / "static_solve.py"
    spec = importlib.util.spec_from_file_location("frozen_modeiii_static_solve_v3", source_path)
    if spec is None or spec.loader is None:
        raise CouplingUnknown("static_solve_module_import_spec_unavailable")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _parse_json_without_nonfinite(path: Path) -> dict[str, Any]:
    def reject_constant(value: str) -> None:
        raise ValueError(f"nonfinite_json_constant:{value}")

    value = json.loads(path.read_text(encoding="utf-8"), parse_constant=reject_constant)
    if type(value) is not dict:
        raise CouplingUnknown("solve_result_root_must_be_object")
    return value


def _cut_parameter_for_negative_x(a: float, b: float, rotation: float) -> float:
    root = math.atan2(-a * math.sin(rotation), b * math.cos(rotation))
    candidates = [root % (2.0 * math.pi), (root + math.pi) % (2.0 * math.pi)]
    for candidate in candidates:
        x = a * math.cos(candidate) * math.cos(rotation) - b * math.sin(candidate) * math.sin(rotation)
        if x < 0.0:
            return candidate
    raise CouplingUnknown("ellipse_negative_x_cut_not_found")


def _ellipse_quadrature(np: Any, a: float, b: float, rotation: float, order_per_arc: int):
    """Copy the frozen analytic runner's two smooth arcs, without analytic fields."""
    nodes, weights = np.polynomial.legendre.leggauss(order_per_arc)
    cut_t = _cut_parameter_for_negative_x(a, b, rotation)
    xyz_rows: list[Any] = []
    normal_rows: list[Any] = []
    ds_rows: list[Any] = []
    for lower, upper in ((0.0, cut_t), (cut_t, 2.0 * math.pi)):
        mid = 0.5 * (lower + upper)
        half = 0.5 * (upper - lower)
        t = mid + half * nodes
        dt_weight = half * weights
        cos_t = np.cos(t)
        sin_t = np.sin(t)
        cos_phi = math.cos(rotation)
        sin_phi = math.sin(rotation)
        x_unrotated = a * cos_t
        y_unrotated = b * sin_t
        x = cos_phi * x_unrotated - sin_phi * y_unrotated
        y = sin_phi * x_unrotated + cos_phi * y_unrotated
        tangent_speed = np.hypot(a * sin_t, b * cos_t)
        nx_unrotated = b * cos_t
        ny_unrotated = a * sin_t
        normal_length = np.hypot(nx_unrotated, ny_unrotated)
        nx = (cos_phi * nx_unrotated - sin_phi * ny_unrotated) / normal_length
        ny = (sin_phi * nx_unrotated + cos_phi * ny_unrotated) / normal_length
        xyz_rows.append(np.column_stack((x, y, np.zeros_like(x))))
        normal_rows.append(np.column_stack((nx, ny, np.zeros_like(nx))))
        ds_rows.append(dt_weight * tangent_speed)
    return (
        np.ascontiguousarray(np.concatenate(xyz_rows, axis=0), dtype=np.float64),
        np.ascontiguousarray(np.concatenate(normal_rows, axis=0), dtype=np.float64),
        np.ascontiguousarray(np.concatenate(ds_rows, axis=0), dtype=np.float64),
    )


def _build_contours(np: Any) -> tuple[list[dict[str, Any]], Any]:
    contours: list[dict[str, Any]] = []
    all_xyz: list[Any] = []
    semiaxes = ((0.0030, 0.0018, "small"), (0.0060, 0.0036, "large"))
    rotations = (0.0, 0.37)
    orders = (12, 16)
    row_start = 0
    for a, b, size_label in semiaxes:
        for rotation in rotations:
            for order in orders:
                xyz, normals, ds = _ellipse_quadrature(np, a, b, rotation, order)
                if xyz.shape != (2 * order, 3) or normals.shape != xyz.shape or ds.shape != (2 * order,):
                    raise CouplingUnknown("frozen_contour_shape_mismatch")
                if bool(np.any((xyz[:, 0] == 0.0) & (xyz[:, 1] == 0.0))):
                    raise CouplingUnknown("contour_contains_exact_front_point")
                if bool(np.any((xyz[:, 0] < 0.0) & (xyz[:, 1] == 0.0))):
                    raise CouplingUnknown("contour_contains_exact_branch_face_point")
                stop = row_start + int(xyz.shape[0])
                contours.append(
                    {
                        "contour_id": f"{size_label}_rotation_{rotation:.2f}_gl{order}_per_arc",
                        "radius_scale_label": size_label,
                        "ellipse_semiaxes_m": {"a": a, "b": b},
                        "rotation_rad": rotation,
                        "quadrature_order_per_smooth_arc": order,
                        "smooth_arc_count": 2,
                        "provider_row_slice_half_open": [row_start, stop],
                        "coordinates_m": xyz,
                        "normals_outward": normals,
                        "weights_ds_m": ds,
                    }
                )
                all_xyz.append(xyz)
                row_start = stop
    if row_start != CONTOUR_QUERY_COUNT:
        raise CouplingUnknown(f"contour_query_budget_mismatch:{row_start}")
    return contours, np.ascontiguousarray(np.concatenate(all_xyz, axis=0), dtype=np.float64)


def _build_reflection_controls(np: Any, contours: list[dict[str, Any]]) -> tuple[Any, list[dict[str, Any]]]:
    base = next(
        entry
        for entry in contours
        if entry["radius_scale_label"] == "small"
        and entry["rotation_rad"] == 0.0
        and entry["quadrature_order_per_smooth_arc"] == 12
    )
    # Select 12 finite in-plane locations spread across both smooth arcs, then
    # query a reflected pair at z=+/-2.5 mm for each location.
    xy = base["coordinates_m"][::2, :2]
    if xy.shape != (CONTROL_QUERY_COUNT // 2, 2):
        raise CouplingUnknown("reflection_control_xy_count_mismatch")
    controls: list[dict[str, Any]] = []
    control_xyz: list[tuple[float, float, float]] = []
    z_offset = 0.0025
    for location_index, point in enumerate(xy):
        pair_indices: list[int] = []
        for sign, label in ((-1.0, "minus_z"), (1.0, "plus_z")):
            pair_indices.append(len(control_xyz))
            control_xyz.append((float(point[0]), float(point[1]), sign * z_offset))
        controls.append(
            {
                "pair_id": f"small_gl12_location_{location_index:02d}",
                "xy_m": [float(point[0]), float(point[1])],
                "minus_z_coordinates_m": [float(point[0]), float(point[1]), -z_offset],
                "plus_z_coordinates_m": [float(point[0]), float(point[1]), z_offset],
                "minus_z_query_index": pair_indices[0],
                "plus_z_query_index": pair_indices[1],
                "z_offset_m": z_offset,
            }
        )
    return np.ascontiguousarray(np.asarray(control_xyz, dtype=np.float64)), controls


def _finite_array(np: Any, value: Any, expected_shape: tuple[int, ...], name: str):
    if type(value) is not np.ndarray or value.shape != expected_shape or not bool(np.isfinite(value).all()):
        raise CouplingUnknown(f"invalid_provider_{name}_shape_or_values")
    return value


def _coefficient_binding(
    np: Any,
    result: dict[str, Any],
    level: dict[str, Any],
    registry: Any,
) -> tuple[Any, dict[str, Any]]:
    canonical_ids = result.get("registry", {}).get("ordered_canonical_coefficient_ids")
    coefficient_values = level.get("solution", {}).get("coefficients_m")
    if type(canonical_ids) is not list or len(canonical_ids) != 54 or any(type(item) is not str for item in canonical_ids):
        raise CouplingUnknown("invalid_solved_canonical_coefficient_ids")
    if len(set(canonical_ids)) != 54 or tuple(canonical_ids) != tuple(registry.canonical_ids):
        raise CouplingUnknown("solved_ids_do_not_match_frozen_provider_registry_order")
    if type(coefficient_values) is not list or len(coefficient_values) != 54:
        raise CouplingUnknown("invalid_solved_coefficient_vector")
    if any(isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(float(value)) for value in coefficient_values):
        raise CouplingUnknown("nonfinite_solved_coefficient_vector")

    # Bind the solved d values to semantic IDs first, then order that mapping
    # using the frozen provider registry's canonical ID sequence.
    coefficient_by_id = {
        canonical_ids[index]: float(coefficient_values[index])
        for index in range(len(canonical_ids))
    }
    provider_ordered = np.ascontiguousarray(
        np.asarray([coefficient_by_id[key] for key in registry.canonical_ids], dtype=np.float64)
    )
    binding_rows = [[key, coefficient_by_id[key]] for key in registry.canonical_ids]
    binding_hash = _sha256_bytes(
        json.dumps(binding_rows, ensure_ascii=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
    )
    evidence = {
        "binding_method": "source-result canonical ID -> coefficient map -> frozen provider registry canonical ID order",
        "source_registry_ids_exactly_match_provider_registry": True,
        "ordered_canonical_coefficient_ids": list(registry.canonical_ids),
        "coefficients_m_in_provider_slot_order": [float(value) for value in provider_ordered],
        "semantic_binding_sha256": binding_hash,
    }
    return provider_ordered, evidence


def _comparison_record(J_value: float | None, expected_J: float) -> dict[str, Any]:
    if J_value is None:
        return {
            "status": "NOT_COMPUTABLE_WHILE_INTEGRATOR_J_UNKNOWN",
            "analytic_reference_J_J_per_m2": expected_J,
            "absolute_difference_J_per_m2": None,
            "relative_difference": None,
        }
    delta = J_value - expected_J
    return {
        "status": "COMPARISON_ONLY",
        "analytic_reference_J_J_per_m2": expected_J,
        "absolute_difference_J_per_m2": abs(delta),
        "relative_difference": delta / expected_J,
    }


def _run_readback() -> dict[str, Any]:
    start = time.monotonic()
    started_at_utc = _utc_now()
    input_hash = _sha256_path(INPUT_RESULT)
    if input_hash != SOLVE_RESULT_SHA256:
        raise CouplingUnknown(f"solve_result_hash_mismatch:{input_hash}")
    result = _parse_json_without_nonfinite(INPUT_RESULT)
    if (
        result.get("schema") != SOLVE_RESULT_SCHEMA
        or result.get("case_id") != SOLVE_CASE_ID
        or result.get("status") != "SOLVED_DIAGNOSTIC_ONLY_UNQUALIFIED"
    ):
        raise CouplingUnknown("input_result_not_the_frozen_solved_diagnostic")
    if result.get("solver_source_sha256", "").lower() != SOLVE_SOURCE_SHA256:
        raise CouplingUnknown("input_result_solver_source_hash_mismatch")

    # Source hashes are checked before importing any numerical dependencies.
    source_pins = _assert_source_pins()
    static = _load_static_module()
    observed_static_pins = static._verify_source_pins()
    if any(observed_static_pins.get(Path(name).name) != expected for name, expected in EXPECTED_DEPENDENCIES.items() if Path(name).name in observed_static_pins):
        raise CouplingUnknown("static_solve_dependency_pin_mismatch")
    static._load_numerical_dependencies()
    np = static.np

    solve_case, base_case = static._read_and_validate_cases()
    source_positions, source_ids = static._build_sources(base_case)
    registry = static._build_registry(source_ids)
    if (
        registry.dof_count != 54
        or tuple(registry.canonical_ids) != tuple(result.get("registry", {}).get("ordered_canonical_coefficient_ids", ()))
    ):
        raise CouplingUnknown("registry_semantic_id_alignment_failed")
    _, mls_config, field_case = static._build_material_cases(base_case)

    contours, contour_xyz = _build_contours(np)
    control_xyz, control_pairs = _build_reflection_controls(np, contours)
    query_positions = np.ascontiguousarray(np.vstack((contour_xyz, control_xyz)), dtype=np.float64)
    if query_positions.shape != (PLANNED_PROVIDER_QUERIES, 3):
        raise CouplingUnknown(f"provider_query_count_mismatch:{query_positions.shape}")
    if query_positions.shape[0] > MAX_PROVIDER_QUERIES:
        raise CouplingUnknown("provider_query_budget_exceeded")

    provider = static._build_provider(
        source_positions,
        source_ids,
        query_positions,
        registry,
        mls_config,
        field_case,
    )
    if provider.unknown_rows:
        raise CouplingUnknown(f"provider_unknown_rows:{json.dumps(provider.unknown_rows, sort_keys=True)}")
    if provider.query_count != PLANNED_PROVIDER_QUERIES:
        raise CouplingUnknown("provider_did_not_use_the_single_frozen_query_batch")
    raw_values = _finite_array(np, provider.values, (PLANNED_PROVIDER_QUERIES, 54), "values")
    raw_gradients = _finite_array(np, provider.gradients, (PLANNED_PROVIDER_QUERIES, 54, 3), "gradients")
    if tuple(provider.coefficient_ids) != tuple(registry.canonical_ids):
        raise CouplingUnknown("provider_rows_are_not_in_registered_semantic_id_order")
    # _build_provider already sampled this entire bounded batch and returned
    # its frozen value/gradient matrices. Contract those rows directly; do not
    # issue a second round of per-point evaluate calls.
    values = raw_values
    gradients = raw_gradients
    if type(provider.consumed_queries) is not int or provider.consumed_queries < 0:
        raise CouplingUnknown("provider_contract_evaluation_count_malformed")

    material = base_case["material_and_manufactured_field"]
    E = float(material["young_modulus_Pa"])
    nu = float(material["poisson_ratio"])
    KIII = float(material["K_III_Pa_sqrt_m"])
    mu = float(E / (2.0 * (1.0 + nu)))
    analytic_reference_J = float(KIII * KIII / (2.0 * mu))

    level_by_name = {str(level.get("name")): level for level in result.get("volume_levels", [])}
    required_level_names = ("volume_gl4x8x4", "volume_gl4x16x4")
    if set(level_by_name) != set(required_level_names):
        raise CouplingUnknown("input_result_volume_levels_changed")

    solved_by_level: dict[str, dict[str, Any]] = {}
    for name in required_level_names:
        level = level_by_name[name]
        if level.get("status") != "OK" or level.get("solution", {}).get("status") != "OK":
            raise CouplingUnknown(f"input_solution_not_ok:{name}")
        coefficient_vector, binding = _coefficient_binding(np, result, level, registry)
        solved_uz = np.ascontiguousarray(values @ coefficient_vector, dtype=np.float64)
        solved_gradient = np.ascontiguousarray(
            np.einsum("qdi,d->qi", gradients, coefficient_vector), dtype=np.float64
        )
        _finite_array(np, solved_uz, (PLANNED_PROVIDER_QUERIES,), f"{name}_uz")
        _finite_array(np, solved_gradient, (PLANNED_PROVIDER_QUERIES, 3), f"{name}_gradient")
        solved_by_level[name] = {
            "level": level,
            "coefficients": coefficient_vector,
            "coefficient_binding": binding,
            "uz": solved_uz,
            "gradient": solved_gradient,
        }

    if str(J_DIR) not in sys.path:
        sys.path.insert(0, str(J_DIR))
    from j_integral_precopy_guard_v3 import INPUT_UNITS, integrate_modeiii_j

    contour_records: list[dict[str, Any]] = []
    for contour in contours:
        begin, stop = contour["provider_row_slice_half_open"]
        xyz = contour["coordinates_m"]
        normals = contour["normals_outward"]
        weights = contour["weights_ds_m"]
        contour_record: dict[str, Any] = {
            "contour_id": contour["contour_id"],
            "radius_scale_label": contour["radius_scale_label"],
            "ellipse_semiaxes_m": contour["ellipse_semiaxes_m"],
            "rotation_rad": contour["rotation_rad"],
            "quadrature_order_per_smooth_arc": contour["quadrature_order_per_smooth_arc"],
            "smooth_arc_count": contour["smooth_arc_count"],
            "sample_count": int(stop - begin),
            "provider_query_slice_half_open": [begin, stop],
            "coordinates_m": xyz.tolist(),
            "normals_outward": normals.tolist(),
            "weights_ds_m": weights.tolist(),
            "by_solved_volume_level": {},
        }
        for name, solved in solved_by_level.items():
            uz = np.ascontiguousarray(solved["uz"][begin:stop], dtype=np.float64)
            grad = np.ascontiguousarray(solved["gradient"][begin:stop, :], dtype=np.float64)
            provenance = (
                f"direct readback from frozen solved coefficients; source result SHA256={input_hash}; "
                f"volume level={name}; provider rows bound by canonical semantic coefficient IDs"
            )
            j_result = integrate_modeiii_j(
                contour_coordinates_m=np.ascontiguousarray(xyz, dtype=np.float64),
                normals=np.ascontiguousarray(normals, dtype=np.float64),
                weights_m=np.ascontiguousarray(weights, dtype=np.float64),
                grad_u_z_m_per_m=grad,
                shear_modulus_Pa=mu,
                units=INPUT_UNITS,
                field_provenance=provenance,
            )
            j_dict = j_result.to_dict()
            J_value = j_result.J_J_per_m2 if j_result.status == "OK" else None
            if J_value is not None and not math.isfinite(float(J_value)):
                raise CouplingUnknown(f"nonfinite_integrator_J:{contour['contour_id']}:{name}")
            contour_record["by_solved_volume_level"][name] = {
                "solved_uz_m": uz.tolist(),
                "solved_grad_u_z_m_per_m_full_xyz": grad.tolist(),
                "grad_z_max_abs_m_per_m": float(np.max(np.abs(grad[:, 2]))),
                "grad_xy_max_abs_m_per_m": float(np.max(np.abs(grad[:, :2]))),
                "integrator_result": j_dict,
                "analytic_reference_comparison": _comparison_record(
                    float(J_value) if J_value is not None else None,
                    analytic_reference_J,
                ),
                "interpretation": "midsection J diagnostic only; cannot establish full-domain 2D behavior or physical G",
            }
        contour_records.append(contour_record)

    control_start = CONTOUR_QUERY_COUNT
    control_records: list[dict[str, Any]] = []
    reflection_metrics: dict[str, Any] = {}
    for name, solved in solved_by_level.items():
        pair_deltas: list[dict[str, Any]] = []
        uz_values = solved["uz"]
        grad_values = solved["gradient"]
        for pair in control_pairs:
            minus_i = control_start + int(pair["minus_z_query_index"])
            plus_i = control_start + int(pair["plus_z_query_index"])
            minus_u = float(uz_values[minus_i])
            plus_u = float(uz_values[plus_i])
            minus_g = [float(value) for value in grad_values[minus_i]]
            plus_g = [float(value) for value in grad_values[plus_i]]
            pair_deltas.append(
                {
                    **pair,
                    "minus_z_uz_m": minus_u,
                    "plus_z_uz_m": plus_u,
                    "uz_plus_minus_difference_m": plus_u - minus_u,
                    "minus_z_grad_u_z_m_per_m_full_xyz": minus_g,
                    "plus_z_grad_u_z_m_per_m_full_xyz": plus_g,
                    "grad_plus_minus_difference_m_per_m_full_xyz": [
                        plus_g[index] - minus_g[index] for index in range(3)
                    ],
                }
            )
        plus_minus_u = [abs(row["uz_plus_minus_difference_m"]) for row in pair_deltas]
        reflected_gradient_delta = [
            math.sqrt(sum(value * value for value in row["grad_plus_minus_difference_m_per_m_full_xyz"]))
            for row in pair_deltas
        ]
        all_control_gradients = grad_values[control_start:, :]
        reflection_metrics[name] = {
            "pair_count": len(pair_deltas),
            "maximum_reflected_uz_difference_m": max(plus_minus_u),
            "maximum_reflected_gradient_difference_m_per_m": max(reflected_gradient_delta),
            "maximum_abs_grad_z_m_per_m_at_controls": float(np.max(np.abs(all_control_gradients[:, 2]))),
            "maximum_abs_grad_xy_m_per_m_at_controls": float(np.max(np.abs(all_control_gradients[:, :2]))),
            "status": "CONTROL_DIAGNOSTIC_ONLY_NO_Z_INVARIANCE_PASS_ASSIGNED",
        }
        control_records.append(
            {
                "volume_level": name,
                "points": pair_deltas,
                "summary": reflection_metrics[name],
            }
        )

    path_sensitivity: dict[str, Any] = {}
    for name in required_level_names:
        by_geometry: dict[str, list[dict[str, Any]]] = {}
        for contour in contour_records:
            key = f"{contour['radius_scale_label']}_rotation_{contour['rotation_rad']:.2f}"
            j_record = contour["by_solved_volume_level"][name]["integrator_result"]
            by_geometry.setdefault(key, []).append(
                {
                    "contour_id": contour["contour_id"],
                    "quadrature_order_per_smooth_arc": contour["quadrature_order_per_smooth_arc"],
                    "J_J_per_m2": j_record.get("J_J_per_m2"),
                    "status": j_record.get("status"),
                    "reason": j_record.get("reason"),
                }
            )
        order_differences: list[dict[str, Any]] = []
        for geometry, rows in by_geometry.items():
            by_order = {int(row["quadrature_order_per_smooth_arc"]): row for row in rows}
            gl12 = by_order.get(12, {}).get("J_J_per_m2")
            gl16 = by_order.get(16, {}).get("J_J_per_m2")
            order_differences.append(
                {
                    "geometry": geometry,
                    "GL12_J_J_per_m2": gl12,
                    "GL16_J_J_per_m2": gl16,
                    "absolute_difference_J_per_m2": (
                        abs(float(gl12) - float(gl16))
                        if gl12 is not None and gl16 is not None
                        else None
                    ),
                    "status": "AVAILABLE" if gl12 is not None and gl16 is not None else "UNKNOWN_CONTOUR_RESULT",
                }
            )
        numeric_contours = [
            {
                "contour_id": contour["contour_id"],
                "radius_scale_label": contour["radius_scale_label"],
                "rotation_rad": contour["rotation_rad"],
                "quadrature_order_per_smooth_arc": contour["quadrature_order_per_smooth_arc"],
                "J_J_per_m2": contour["by_solved_volume_level"][name]["integrator_result"].get("J_J_per_m2"),
            }
            for contour in contour_records
            if contour["by_solved_volume_level"][name]["integrator_result"].get("J_J_per_m2") is not None
        ]
        same_order_geometry_changes: list[dict[str, Any]] = []
        for left_index, left in enumerate(numeric_contours):
            for right in numeric_contours[left_index + 1 :]:
                same_geometry = (
                    left["radius_scale_label"] == right["radius_scale_label"]
                    and left["rotation_rad"] == right["rotation_rad"]
                )
                if left["quadrature_order_per_smooth_arc"] == right["quadrature_order_per_smooth_arc"] and not same_geometry:
                    same_order_geometry_changes.append(
                        {
                            "quadrature_order_per_smooth_arc": left["quadrature_order_per_smooth_arc"],
                            "contour_a": left["contour_id"],
                            "contour_b": right["contour_id"],
                            "absolute_J_difference_J_per_m2": abs(
                                float(left["J_J_per_m2"]) - float(right["J_J_per_m2"])
                            ),
                        }
                    )
        geometry_range_by_order: dict[str, dict[str, Any]] = {}
        for order in (12, 16):
            same_order = [row for row in numeric_contours if row["quadrature_order_per_smooth_arc"] == order]
            same_order_values = [float(row["J_J_per_m2"]) for row in same_order]
            geometry_range_by_order[str(order)] = {
                "numeric_geometry_count": len(same_order_values),
                "minimum_J_J_per_m2": min(same_order_values) if same_order_values else None,
                "maximum_J_J_per_m2": max(same_order_values) if same_order_values else None,
                "absolute_spread_across_size_and_rotation_J_per_m2": (
                    max(same_order_values) - min(same_order_values) if same_order_values else None
                ),
                "comparison_scope": "different ellipse size and/or rotation at this fixed GL order",
            }
        path_sensitivity[name] = {
            "groups_by_size_and_rotation": by_geometry,
            "quadrature_sensitivity_same_geometry_GL12_GL16": order_differences,
            "cross_geometry_sensitivity_same_quadrature_order": same_order_geometry_changes,
            "geometry_range_by_fixed_quadrature_order": geometry_range_by_order,
            "status": "SENSITIVITY_DIAGNOSTIC_ONLY; this comparison does not establish path independence; unavailable groups remain UNKNOWN without numeric J",
        }

    query_array_bytes = int(query_positions.nbytes + values.nbytes + gradients.nbytes)
    explicit_numeric_bytes = int(
        source_positions.nbytes
        + query_positions.nbytes
        + raw_values.nbytes
        + raw_gradients.nbytes
        + values.nbytes
        + gradients.nbytes
        + control_xyz.nbytes
        + sum(
            entry["coordinates_m"].nbytes
            + entry["normals_outward"].nbytes
            + entry["weights_ds_m"].nbytes
            for entry in contours
        )
        + sum(item["uz"].nbytes + item["gradient"].nbytes + item["coefficients"].nbytes for item in solved_by_level.values())
    )
    max_grad_z = {
        name: max(float(record["by_solved_volume_level"][name]["grad_z_max_abs_m_per_m"]) for record in contour_records)
        for name in required_level_names
    }
    unknown_j_count = sum(
        1
        for record in contour_records
        for name in required_level_names
        if record["by_solved_volume_level"][name]["integrator_result"].get("status") != "OK"
    )
    elapsed = time.monotonic() - start
    if elapsed > WALL_TIMEOUT_SECONDS:
        raise CouplingUnknown(f"worker_wall_budget_exceeded:{elapsed:.6f}")

    return {
        "schema": RESULT_SCHEMA,
        "case_id": CASE_ID,
        "status": "READBACK_COMPLETE_UNQUALIFIED",
        "classification": "SOLVED_FIELD_MIDSECTION_J_COUPLING_DIAGNOSTIC_ONLY",
        "started_at_utc": started_at_utc,
        "finished_at_utc": _utc_now(),
        "source_solution": {
            "path": str(INPUT_RESULT.relative_to(PROJECT_ROOT)).replace("\\", "/"),
            "sha256": input_hash,
            "schema": result["schema"],
            "case_id": result["case_id"],
            "status": result["status"],
            "solver_source_sha256": result["solver_source_sha256"],
            "source_freeze_sha256": SOLVE_FREEZE_SHA256,
            "levels": [
                {
                    "name": name,
                    "quadrature_rule": level_by_name[name].get("quadrature_rule"),
                    "solution_status": level_by_name[name]["solution"]["status"],
                    "semantic_coefficient_binding": solved_by_level[name]["coefficient_binding"],
                }
                for name in required_level_names
            ],
        },
        "source_hashes": source_pins,
        "provider": {
            "class": "ContinuousMLSFrontProvider",
            "source_sha256": source_pins["pgfrag/experiments/modeiii_pu_operator_20261005/mls_modeiii_provider.py"],
            "provider_batch_count": 1,
            "sampled_batch_query_count": int(provider.query_count),
            "contour_query_count": CONTOUR_QUERY_COUNT,
            "z_reflection_control_query_count": CONTROL_QUERY_COUNT,
            "maximum_query_count": MAX_PROVIDER_QUERIES,
            "unknown_row_count": len(provider.unknown_rows),
            "coefficient_ids_match_registry": True,
            "provider_condition": static._condition_summary([provider]),
            "provider_stored_numeric_bytes_values_gradients": query_array_bytes,
            "contract_evaluation_count": int(provider.consumed_queries),
            "contract_evaluation_count_semantics": "direct contraction of the single sampled batch matrices; no provider.evaluate calls",
        },
        "J_integrator": {
            "module": J_INTEGRATOR_NAME,
            "source_sha256": J_INTEGRATOR_SHA256,
            "formula": "J = integral_[Gamma] (W*n_x - t_z*du_z/dx) ds",
            "units": {**INPUT_UNITS},
            "shear_modulus_Pa": mu,
            "analytic_reference_J_J_per_m2_comparison_only": analytic_reference_J,
            "unknown_contour_level_count": unknown_j_count,
            "max_abs_grad_z_on_midsection_contours_m_per_m_by_level": max_grad_z,
            "grad_z_projection_or_replacement": False,
            "contour_geometry_plane": "global x-y at z=0",
            "interpretation": "Every z=0 contour result is diagnostic only. Integrator UNKNOWN is preserved verbatim, including any z-invariance reason. This experiment does not prove full-domain 2D behavior or physical G.",
        },
        "contours": contour_records,
        "z_reflection_controls": {
            "z_locations_m": [-0.0025, 0.0025],
            "control_point_pair_count": len(control_pairs),
            "total_provider_queries": CONTROL_QUERY_COUNT,
            "points_by_solved_volume_level": control_records,
            "interpretation": "actual u_z and all three gradient components are retained; controls assign no z-invariance PASS and do not qualify J or G",
        },
        "path_sensitivity": path_sensitivity,
        "budgets": {
            "provider_queries": int(provider.query_count),
            "provider_query_cap": MAX_PROVIDER_QUERIES,
            "wall_elapsed_seconds": elapsed,
            "wall_cap_seconds": WALL_TIMEOUT_SECONDS,
            "tracked_explicit_numeric_array_bytes_estimate": explicit_numeric_bytes,
            "tracked_array_estimate_scope": "adapter-held arrays and provider values/gradients; provider MLS/exact-field internals and Python object overhead are excluded",
            "result_json_cap_bytes": MAX_RESULT_JSON_BYTES,
            "result_json_bytes": None,
            "no_stiffness_reassembly": True,
            "no_static_resolve": True,
        },
        "limitations": [
            "The coefficient values are read from the existing frozen solve result; the adapter does not assemble K or invoke the static solver.",
            "The exact field reference is used only inside the frozen basis provider and as an independent scalar J comparison; solved field outputs come from phi/dphi contracted with the solved coefficients.",
            "The full gradient, including its native z component, is passed unchanged to the frozen J integrator.",
            "A z=0 contour is a midsection diagnostic and cannot establish full-domain two-dimensional behavior, path-independent physical G, or crack-growth readiness.",
        ],
    }


def _write_json(path: Path, value: dict[str, Any]) -> int:
    encoded = json.dumps(value, ensure_ascii=True, indent=2, allow_nan=False).encode("utf-8") + b"\n"
    if len(encoded) > MAX_RESULT_JSON_BYTES:
        raise CouplingUnknown(f"result_json_budget_exceeded:{len(encoded)}")
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {path}")
    path.write_bytes(encoded)
    return len(encoded)


def _worker(output_path: Path) -> int:
    try:
        payload = _run_readback()
        # First serialization measures the full evidence size. If the declared
        # output cap is exceeded, leave a small, explicit UNKNOWN artifact.
        payload["budgets"]["result_json_bytes"] = 0
        encoded = json.dumps(payload, ensure_ascii=True, indent=2, allow_nan=False).encode("utf-8") + b"\n"
        if len(encoded) > MAX_RESULT_JSON_BYTES:
            payload = {
                "schema": RESULT_SCHEMA,
                "case_id": CASE_ID,
                "status": "UNKNOWN_RESULT_JSON_BUDGET",
                "reason": f"result JSON would be {len(encoded)} bytes, cap is {MAX_RESULT_JSON_BYTES}",
                "source_solution_sha256": SOLVE_RESULT_SHA256,
            }
        else:
            for _ in range(3):
                payload["budgets"]["result_json_bytes"] = len(encoded)
                encoded = json.dumps(payload, ensure_ascii=True, indent=2, allow_nan=False).encode("utf-8") + b"\n"
            if len(encoded) > MAX_RESULT_JSON_BYTES:
                raise CouplingUnknown(f"result_json_budget_exceeded:{len(encoded)}")
        _write_json(output_path, payload)
        return 0 if payload.get("status") == "READBACK_COMPLETE_UNQUALIFIED" else 2
    except BaseException as exc:
        failure = {
            "schema": RESULT_SCHEMA,
            "case_id": CASE_ID,
            "status": "UNKNOWN_READBACK_FAILURE",
            "reason": f"{type(exc).__name__}: {exc}",
            "traceback": traceback.format_exc(limit=6),
            "source_solution_sha256": SOLVE_RESULT_SHA256,
            "source_hashes": {
                "static_solve.py": SOLVE_SOURCE_SHA256,
                "j_integral_precopy_guard_v3.py": J_INTEGRATOR_SHA256,
            },
        }
        try:
            if not output_path.exists():
                _write_json(output_path, failure)
        except BaseException:
            pass
        return 2


def _write_parent_failure(path: Path, status: str, reason: str, stdout: str = "", stderr: str = "") -> None:
    failure = {
        "schema": RESULT_SCHEMA,
        "case_id": CASE_ID,
        "status": status,
        "reason": reason,
        "source_solution_sha256": SOLVE_RESULT_SHA256,
        "wall_cap_seconds": WALL_TIMEOUT_SECONDS,
    }
    _write_json(path, failure)
    path.parent.mkdir(parents=True, exist_ok=True)
    (path.parent / "stdout.log").write_text(stdout, encoding="utf-8")
    (path.parent / "stderr.log").write_text(stderr, encoding="utf-8")


def _is_finite_number(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value))


def _finite_tree(value: object) -> bool:
    if isinstance(value, (int, float)):
        return _is_finite_number(value)
    if type(value) is list:
        return all(_finite_tree(item) for item in value)
    return False


def _finite_json_numbers(value: object) -> bool:
    if type(value) is dict:
        return all(type(key) is str and _finite_json_numbers(item) for key, item in value.items())
    if type(value) is list:
        return all(_finite_json_numbers(item) for item in value)
    if isinstance(value, float):
        return math.isfinite(value)
    return value is None or type(value) in (str, int, bool)


def _valid_utc_timestamp(value: object) -> bool:
    if type(value) is not str:
        return False
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        return False
    return parsed.tzinfo is not None and parsed.utcoffset() is not None and parsed.utcoffset().total_seconds() == 0


def _validate_staged_result(value: dict[str, Any], byte_count: int) -> str | None:
    """Check the complete promoted-result envelope without importing NumPy."""
    if value.get("schema") != RESULT_SCHEMA or value.get("case_id") != CASE_ID:
        return "schema_or_case_id_mismatch"
    if value.get("status") != "READBACK_COMPLETE_UNQUALIFIED":
        return f"worker_status_not_complete:{value.get('status')}"
    if not _finite_json_numbers(value):
        return "nonfinite_or_non_json_value_in_result"
    if (
        value.get("classification") != "SOLVED_FIELD_MIDSECTION_J_COUPLING_DIAGNOSTIC_ONLY"
        or not _valid_utc_timestamp(value.get("started_at_utc"))
        or not _valid_utc_timestamp(value.get("finished_at_utc"))
    ):
        return "result_classification_or_timestamps_invalid"
    source_solution = value.get("source_solution")
    if (
        type(source_solution) is not dict
        or source_solution.get("path") != str(INPUT_RESULT.relative_to(PROJECT_ROOT)).replace("\\", "/")
        or source_solution.get("sha256") != SOLVE_RESULT_SHA256
        or source_solution.get("schema") != SOLVE_RESULT_SCHEMA
        or source_solution.get("case_id") != SOLVE_CASE_ID
        or source_solution.get("status") != "SOLVED_DIAGNOSTIC_ONLY_UNQUALIFIED"
        or source_solution.get("solver_source_sha256", "").lower() != SOLVE_SOURCE_SHA256
        or source_solution.get("source_freeze_sha256", "").lower() != SOLVE_FREEZE_SHA256
    ):
        return "source_solution_identity_or_hash_mismatch"
    levels = source_solution.get("levels")
    expected_levels = ("volume_gl4x8x4", "volume_gl4x16x4")
    if (
        type(levels) is not list
        or len(levels) != 2
        or any(type(row) is not dict for row in levels)
        or tuple(row.get("name") for row in levels) != expected_levels
    ):
        return "source_solution_level_records_malformed"
    previous_ids: list[str] | None = None
    for level in levels:
        binding = level.get("semantic_coefficient_binding")
        if (
            type(binding) is not dict
            or level.get("solution_status") != "OK"
            or type(level.get("quadrature_rule")) is not str
            or not level.get("quadrature_rule")
            or binding.get("binding_method") != "source-result canonical ID -> coefficient map -> frozen provider registry canonical ID order"
            or binding.get("source_registry_ids_exactly_match_provider_registry") is not True
        ):
            return "semantic_coefficient_binding_missing"
        ids = binding.get("ordered_canonical_coefficient_ids")
        coefficients = binding.get("coefficients_m_in_provider_slot_order")
        if (
            type(ids) is not list
            or len(ids) != 54
            or any(type(item) is not str or not item for item in ids)
            or len(set(ids)) != 54
            or type(coefficients) is not list
            or len(coefficients) != 54
            or not _finite_tree(coefficients)
            or type(binding.get("semantic_binding_sha256")) is not str
            or len(binding["semantic_binding_sha256"]) != 64
            or any(char not in "0123456789abcdefABCDEF" for char in binding["semantic_binding_sha256"])
        ):
            return "semantic_coefficient_binding_invalid"
        if previous_ids is not None and ids != previous_ids:
            return "semantic_coefficient_ids_differ_between_solve_levels"
        previous_ids = ids
        expected_binding_hash = _sha256_bytes(
            json.dumps(list(zip(ids, coefficients)), ensure_ascii=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
        )
        if binding["semantic_binding_sha256"].lower() != expected_binding_hash:
            return "semantic_coefficient_binding_hash_invalid"

    provider = value.get("provider")
    source_hashes = value.get("source_hashes")
    if type(source_hashes) is not dict:
        return "source_hash_records_missing"
    for relative, expected_hash in EXPECTED_DEPENDENCIES.items():
        observed_hash = source_hashes.get(relative)
        if not isinstance(observed_hash, str) or observed_hash.lower() != expected_hash:
            return f"source_dependency_hash_mismatch:{relative}"
    owned_source_relative = "pgfrag/experiments/solved_field_j_coupling_20261006/read_solved_field_j.py"
    owned_freeze_relative = "pgfrag/tmp/review_20261006/solved_field_j_coupling_20261006/source_freeze.json"
    try:
        if (
            source_hashes.get(owned_source_relative) != _sha256_path(Path(__file__).resolve())
            or source_hashes.get(owned_freeze_relative) != _sha256_path(SOURCE_FREEZE_PATH)
        ):
            return "coupling_owned_source_or_freeze_hash_mismatch"
        freeze_value = json.loads(SOURCE_FREEZE_PATH.read_text(encoding="utf-8"))
        frozen_owned_hash = freeze_value["owned_sources"][owned_source_relative].lower()
    except (OSError, ValueError, KeyError, TypeError, AttributeError):
        return "coupling_source_freeze_unreadable"
    if source_hashes.get(owned_source_relative, "").lower() != frozen_owned_hash:
        return "coupling_source_freeze_owned_source_pin_mismatch"
    if (
        type(provider) is not dict
        or provider.get("class") != "ContinuousMLSFrontProvider"
        or provider.get("source_sha256", "").lower() != EXPECTED_DEPENDENCIES["pgfrag/experiments/modeiii_pu_operator_20261005/mls_modeiii_provider.py"]
        or provider.get("provider_batch_count") != 1
        or provider.get("sampled_batch_query_count") != PLANNED_PROVIDER_QUERIES
        or type(provider.get("contract_evaluation_count")) is not int
        or provider.get("contract_evaluation_count") < 0
        or provider.get("contract_evaluation_count") != 0
        or provider.get("contract_evaluation_count_semantics") != "direct contraction of the single sampled batch matrices; no provider.evaluate calls"
        or provider.get("contour_query_count") != CONTOUR_QUERY_COUNT
        or provider.get("z_reflection_control_query_count") != CONTROL_QUERY_COUNT
        or provider.get("maximum_query_count") != MAX_PROVIDER_QUERIES
        or provider.get("unknown_row_count") != 0
        or provider.get("coefficient_ids_match_registry") is not True
        or not _is_finite_number(provider.get("provider_stored_numeric_bytes_values_gradients"))
        or provider.get("provider_stored_numeric_bytes_values_gradients") <= 0
        or type(provider.get("provider_condition")) is not dict
    ):
        return "provider_query_or_registry_counts_malformed"
    j_integrator = value.get("J_integrator")
    if (
        type(j_integrator) is not dict
        or j_integrator.get("module") != J_INTEGRATOR_NAME
        or str(j_integrator.get("source_sha256", "")).lower() != J_INTEGRATOR_SHA256
        or j_integrator.get("grad_z_projection_or_replacement") is not False
        or not _is_finite_number(j_integrator.get("shear_modulus_Pa"))
        or float(j_integrator["shear_modulus_Pa"]) <= 0.0
        or not _is_finite_number(j_integrator.get("analytic_reference_J_J_per_m2_comparison_only"))
        or type(j_integrator.get("max_abs_grad_z_on_midsection_contours_m_per_m_by_level")) is not dict
        or set(j_integrator["max_abs_grad_z_on_midsection_contours_m_per_m_by_level"]) != set(expected_levels)
        or any(not _is_finite_number(item) for item in j_integrator["max_abs_grad_z_on_midsection_contours_m_per_m_by_level"].values())
        or type(j_integrator.get("interpretation")) is not str
        or "diagnostic only" not in j_integrator["interpretation"]
    ):
        return "J_integrator_identity_or_full_gradient_contract_invalid"

    contours = value.get("contours")
    if type(contours) is not list or len(contours) != 8:
        return "contour_count_mismatch"
    seen_ids: set[str] = set()
    seen_geometry_orders: set[tuple[str, float, int]] = set()
    expected_row_start = 0
    for contour in contours:
        if type(contour) is not dict:
            return "contour_record_not_object"
        contour_id = contour.get("contour_id")
        if type(contour_id) is not str or contour_id in seen_ids:
            return "contour_id_missing_or_duplicated"
        seen_ids.add(contour_id)
        order = contour.get("quadrature_order_per_smooth_arc")
        expected_count = 2 * order if order in (12, 16) else -1
        count = contour.get("sample_count")
        if count != expected_count:
            return f"contour_sample_count_mismatch:{contour_id}"
        geometry_key = (
            str(contour.get("radius_scale_label")),
            float(contour.get("rotation_rad")) if _is_finite_number(contour.get("rotation_rad")) else math.nan,
            int(order) if order in (12, 16) else -1,
        )
        if geometry_key in seen_geometry_orders:
            return f"contour_geometry_order_duplicated:{contour_id}"
        seen_geometry_orders.add(geometry_key)
        expected_semiaxes = {
            "small": {"a": 0.003, "b": 0.0018},
            "large": {"a": 0.006, "b": 0.0036},
        }.get(contour.get("radius_scale_label"))
        if (
            expected_semiaxes is None
            or contour.get("ellipse_semiaxes_m") != expected_semiaxes
            or contour.get("rotation_rad") not in (0.0, 0.37)
            or contour.get("smooth_arc_count") != 2
            or contour_id != f"{contour.get('radius_scale_label')}_rotation_{contour.get('rotation_rad'):.2f}_gl{order}_per_arc"
        ):
            return f"contour_geometry_metadata_invalid:{contour_id}"
        row_slice = contour.get("provider_query_slice_half_open")
        if (
            type(row_slice) is not list
            or len(row_slice) != 2
            or row_slice[0] != expected_row_start
            or row_slice[1] - row_slice[0] != count
        ):
            return f"contour_provider_slice_invalid:{contour_id}"
        expected_row_start = row_slice[1]
        for key, width in (("coordinates_m", 3), ("normals_outward", 3)):
            rows = contour.get(key)
            if type(rows) is not list or len(rows) != count or any(type(row) is not list or len(row) != width or not _finite_tree(row) for row in rows):
                return f"contour_{key}_malformed:{contour_id}"
        if any(row[2] != 0.0 for row in contour["coordinates_m"]):
            return f"contour_not_on_midsection:{contour_id}"
        weights = contour.get("weights_ds_m")
        if type(weights) is not list or len(weights) != count or not _finite_tree(weights) or any(float(item) <= 0.0 for item in weights):
            return f"contour_weights_malformed:{contour_id}"
        fields_by_level = contour.get("by_solved_volume_level")
        if type(fields_by_level) is not dict or set(fields_by_level) != set(expected_levels):
            return f"contour_level_fields_missing:{contour_id}"
        for level_name in expected_levels:
            field = fields_by_level[level_name]
            if type(field) is not dict:
                return f"contour_field_not_object:{contour_id}:{level_name}"
            uz = field.get("solved_uz_m")
            gradient = field.get("solved_grad_u_z_m_per_m_full_xyz")
            if type(uz) is not list or len(uz) != count or not _finite_tree(uz):
                return f"contour_solved_uz_invalid:{contour_id}:{level_name}"
            if type(gradient) is not list or len(gradient) != count or any(type(row) is not list or len(row) != 3 or not _finite_tree(row) for row in gradient):
                return f"contour_solved_gradient_invalid:{contour_id}:{level_name}"
            j_result = field.get("integrator_result")
            if type(j_result) is not dict or j_result.get("status") not in ("OK", "UNKNOWN"):
                return f"contour_integrator_result_invalid:{contour_id}:{level_name}"
            if j_result.get("status") == "OK" and not _is_finite_number(j_result.get("J_J_per_m2")):
                return f"contour_integrator_J_invalid:{contour_id}:{level_name}"
            if j_result.get("status") == "UNKNOWN" and (
                j_result.get("J_J_per_m2") is not None or type(j_result.get("reason")) is not str or not j_result["reason"]
            ):
                return f"contour_integrator_unknown_reason_missing:{contour_id}:{level_name}"
            comparison = field.get("analytic_reference_comparison")
            if type(comparison) is not dict or not _is_finite_number(comparison.get("analytic_reference_J_J_per_m2")):
                return f"contour_comparison_reference_invalid:{contour_id}:{level_name}"
            if not _is_finite_number(field.get("grad_z_max_abs_m_per_m")) or not _is_finite_number(field.get("grad_xy_max_abs_m_per_m")):
                return f"contour_gradient_summary_invalid:{contour_id}:{level_name}"
            if j_result.get("status") == "OK" and (
                not _is_finite_number(comparison.get("absolute_difference_J_per_m2"))
                or not _is_finite_number(comparison.get("relative_difference"))
                or comparison.get("status") != "COMPARISON_ONLY"
            ):
                return f"contour_numeric_J_comparison_invalid:{contour_id}:{level_name}"
            if j_result.get("status") == "UNKNOWN" and (
                comparison.get("status") != "NOT_COMPUTABLE_WHILE_INTEGRATOR_J_UNKNOWN"
                or comparison.get("absolute_difference_J_per_m2") is not None
                or comparison.get("relative_difference") is not None
            ):
                return f"contour_unknown_J_comparison_invalid:{contour_id}:{level_name}"
            if "midsection diagnostic only" not in str(field.get("interpretation", "")):
                return f"contour_interpretation_invalid:{contour_id}:{level_name}"
    expected_geometry_orders = {
        (size, rotation, order)
        for size in ("small", "large")
        for rotation in (0.0, 0.37)
        for order in (12, 16)
    }
    if seen_geometry_orders != expected_geometry_orders or expected_row_start != CONTOUR_QUERY_COUNT:
        return "contour_geometry_order_set_mismatch"

    controls = value.get("z_reflection_controls")
    if (
        type(controls) is not dict
        or controls.get("total_provider_queries") != CONTROL_QUERY_COUNT
        or controls.get("control_point_pair_count") != CONTROL_QUERY_COUNT // 2
        or controls.get("z_locations_m") != [-0.0025, 0.0025]
        or type(controls.get("interpretation")) is not str
        or "no z-invariance PASS" not in controls["interpretation"]
    ):
        return "z_reflection_control_counts_malformed"
    control_levels = controls.get("points_by_solved_volume_level")
    if (
        type(control_levels) is not list
        or len(control_levels) != 2
        or tuple(level.get("volume_level") for level in control_levels if type(level) is dict) != expected_levels
    ):
        return "z_reflection_control_level_count_mismatch"
    for level in control_levels:
        rows = level.get("points") if type(level) is dict else None
        if type(rows) is not list or len(rows) != CONTROL_QUERY_COUNT // 2:
            return "z_reflection_control_pairs_malformed"
        for row in rows:
            if (
                type(row) is not dict
                or not _finite_tree(row.get("xy_m"))
                or len(row.get("xy_m", [])) != 2
                or not _finite_tree(row.get("minus_z_coordinates_m"))
                or not _finite_tree(row.get("plus_z_coordinates_m"))
                or row.get("minus_z_coordinates_m") != [*row.get("xy_m", []), -0.0025]
                or row.get("plus_z_coordinates_m") != [*row.get("xy_m", []), 0.0025]
                or row.get("z_offset_m") != 0.0025
            ):
                return "z_reflection_control_coordinates_invalid"
            if not _finite_tree(row.get("minus_z_grad_u_z_m_per_m_full_xyz")) or not _finite_tree(row.get("plus_z_grad_u_z_m_per_m_full_xyz")):
                return "z_reflection_control_gradient_invalid"
            if not _is_finite_number(row.get("minus_z_uz_m")) or not _is_finite_number(row.get("plus_z_uz_m")):
                return "z_reflection_control_displacement_invalid"
            if not _finite_tree(row.get("grad_plus_minus_difference_m_per_m_full_xyz")):
                return "z_reflection_control_gradient_difference_invalid"
        summary = level.get("summary")
        if (
            type(summary) is not dict
            or summary.get("status") != "CONTROL_DIAGNOSTIC_ONLY_NO_Z_INVARIANCE_PASS_ASSIGNED"
            or summary.get("pair_count") != CONTROL_QUERY_COUNT // 2
            or not _is_finite_number(summary.get("maximum_reflected_uz_difference_m"))
            or not _is_finite_number(summary.get("maximum_reflected_gradient_difference_m_per_m"))
            or not _is_finite_number(summary.get("maximum_abs_grad_z_m_per_m_at_controls"))
            or not _is_finite_number(summary.get("maximum_abs_grad_xy_m_per_m_at_controls"))
        ):
            return "z_reflection_control_summary_invalid"

    path_sensitivity = value.get("path_sensitivity")
    if type(path_sensitivity) is not dict or set(path_sensitivity) != set(expected_levels):
        return "path_sensitivity_levels_malformed"
    expected_geometries = {
        f"{size}_rotation_{rotation:.2f}"
        for size in ("small", "large")
        for rotation in (0.0, 0.37)
    }
    for level_name in expected_levels:
        sensitivity = path_sensitivity[level_name]
        groups = sensitivity.get("groups_by_size_and_rotation") if type(sensitivity) is dict else None
        if type(groups) is not dict or set(groups) != expected_geometries:
            return f"path_sensitivity_geometry_groups_invalid:{level_name}"
        for rows in groups.values():
            if type(rows) is not list or len(rows) != 2 or {row.get("quadrature_order_per_smooth_arc") for row in rows if type(row) is dict} != {12, 16}:
                return f"path_sensitivity_quadrature_groups_invalid:{level_name}"
            for row in rows:
                if (
                    type(row) is not dict
                    or row.get("status") not in ("OK", "UNKNOWN")
                    or (row.get("status") == "OK" and not _is_finite_number(row.get("J_J_per_m2")))
                    or (row.get("status") == "UNKNOWN" and (row.get("J_J_per_m2") is not None or type(row.get("reason")) is not str))
                ):
                    return f"path_sensitivity_contour_value_invalid:{level_name}"
        quadrature = sensitivity.get("quadrature_sensitivity_same_geometry_GL12_GL16")
        if type(quadrature) is not list or len(quadrature) != 4:
            return f"path_sensitivity_quadrature_summary_invalid:{level_name}"
        for row in quadrature:
            if type(row) is not dict or row.get("geometry") not in expected_geometries:
                return f"path_sensitivity_quadrature_record_invalid:{level_name}"
            vals = (row.get("GL12_J_J_per_m2"), row.get("GL16_J_J_per_m2"))
            if any(item is not None and not _is_finite_number(item) for item in vals):
                return f"path_sensitivity_quadrature_value_invalid:{level_name}"
            difference = row.get("absolute_difference_J_per_m2")
            if row.get("status") == "AVAILABLE":
                if any(item is None for item in vals) or not _is_finite_number(difference):
                    return f"path_sensitivity_quadrature_difference_invalid:{level_name}"
            elif row.get("status") != "UNKNOWN_CONTOUR_RESULT" or difference is not None:
                return f"path_sensitivity_quadrature_unknown_invalid:{level_name}"
        cross_geometry = sensitivity.get("cross_geometry_sensitivity_same_quadrature_order")
        if type(cross_geometry) is not list:
            return f"path_sensitivity_cross_geometry_invalid:{level_name}"
        for row in cross_geometry:
            if (
                type(row) is not dict
                or row.get("quadrature_order_per_smooth_arc") not in (12, 16)
                or row.get("contour_a") == row.get("contour_b")
                or not _is_finite_number(row.get("absolute_J_difference_J_per_m2"))
            ):
                return f"path_sensitivity_cross_geometry_record_invalid:{level_name}"
        ranges = sensitivity.get("geometry_range_by_fixed_quadrature_order")
        if type(ranges) is not dict or set(ranges) != {"12", "16"}:
            return f"path_sensitivity_fixed_order_ranges_invalid:{level_name}"
        for order, record in ranges.items():
            if (
                type(record) is not dict
                or record.get("numeric_geometry_count") not in (0, 1, 2, 3, 4)
                or record.get("comparison_scope") != "different ellipse size and/or rotation at this fixed GL order"
            ):
                return f"path_sensitivity_fixed_order_record_invalid:{level_name}:{order}"
            fields = ("minimum_J_J_per_m2", "maximum_J_J_per_m2", "absolute_spread_across_size_and_rotation_J_per_m2")
            if record["numeric_geometry_count"] == 0:
                if any(record.get(field) is not None for field in fields):
                    return f"path_sensitivity_empty_fixed_order_range_invalid:{level_name}:{order}"
            elif any(not _is_finite_number(record.get(field)) for field in fields):
                return f"path_sensitivity_numeric_fixed_order_range_invalid:{level_name}:{order}"
        if sensitivity.get("status") != "SENSITIVITY_DIAGNOSTIC_ONLY; this comparison does not establish path independence; unavailable groups remain UNKNOWN without numeric J":
            return f"path_sensitivity_interpretation_invalid:{level_name}"

    budgets = value.get("budgets")
    if (
        type(budgets) is not dict
        or budgets.get("provider_queries") != PLANNED_PROVIDER_QUERIES
        or budgets.get("provider_query_cap") != MAX_PROVIDER_QUERIES
        or not _is_finite_number(budgets.get("wall_elapsed_seconds"))
        or float(budgets["wall_elapsed_seconds"]) <= 0.0
        or float(budgets["wall_elapsed_seconds"]) > WALL_TIMEOUT_SECONDS
        or budgets.get("result_json_cap_bytes") != MAX_RESULT_JSON_BYTES
        or budgets.get("result_json_bytes") != byte_count
        or byte_count <= 0
        or not _is_finite_number(budgets.get("tracked_explicit_numeric_array_bytes_estimate"))
        or budgets.get("tracked_explicit_numeric_array_bytes_estimate") <= 0
        or budgets.get("no_stiffness_reassembly") is not True
        or budgets.get("no_static_resolve") is not True
    ):
        return "worker_budget_evidence_malformed"
    return None


def _parent_run(output_path: Path) -> int:
    if output_path.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    stage_path = output_path.with_name(f".{output_path.name}.{os.getpid()}.worker.json")
    if stage_path.exists():
        raise FileExistsError(f"refusing to overwrite existing worker staging file: {stage_path}")
    argv = [sys.executable, "-B", str(Path(__file__).resolve()), "--worker-output", str(stage_path)]
    started = time.monotonic()
    stdout = ""
    stderr = ""
    child_return_code: int | None = None
    failure_status: str | None = None
    failure_reason: str | None = None
    try:
        completed = subprocess.run(
            argv,
            cwd=str(PROJECT_ROOT),
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=WORKER_TIMEOUT_SECONDS,
            check=False,
        )
        child_return_code = int(completed.returncode)
        stdout = completed.stdout
        stderr = completed.stderr
        if time.monotonic() - started > WALL_TIMEOUT_SECONDS:
            failure_status = "UNKNOWN_WALL_TIMEOUT"
            failure_reason = f"outer readback exceeded {WALL_TIMEOUT_SECONDS:.1f} seconds; worker staging output retained for review"
        elif child_return_code != 0:
            failure_status = "UNKNOWN_WORKER_NONZERO_EXIT"
            failure_reason = f"worker_return_code={child_return_code}; staging output retained for review"
        elif not stage_path.is_file():
            failure_status = "UNKNOWN_WORKER_NO_OUTPUT"
            failure_reason = f"worker_return_code={child_return_code}; child produced no JSON"
        elif stage_path.stat().st_size > MAX_RESULT_JSON_BYTES:
            failure_status = "UNKNOWN_RESULT_JSON_BUDGET"
            failure_reason = f"worker output is {stage_path.stat().st_size} bytes"
        else:
            staged = _parse_json_without_nonfinite(stage_path)
            contract_error = _validate_staged_result(staged, stage_path.stat().st_size)
            if contract_error is not None:
                failure_status = "UNKNOWN_WORKER_RESULT_CONTRACT"
                failure_reason = contract_error + "; staging output retained for review"
            else:
                if output_path.exists():
                    raise FileExistsError(f"refusing to overwrite concurrent output: {output_path}")
                os.replace(stage_path, output_path)
    except subprocess.TimeoutExpired as exc:
        stdout = exc.stdout.decode("utf-8", errors="replace") if isinstance(exc.stdout, bytes) else (exc.stdout or "")
        stderr = exc.stderr.decode("utf-8", errors="replace") if isinstance(exc.stderr, bytes) else (exc.stderr or "")
        failure_status = "UNKNOWN_WALL_TIMEOUT"
        failure_reason = f"worker exceeded {WORKER_TIMEOUT_SECONDS:.1f} seconds"
    except BaseException as exc:
        failure_status = "UNKNOWN_PARENT_LAUNCH_OR_PROMOTION"
        failure_reason = f"{type(exc).__name__}: {exc}"

    elapsed = time.monotonic() - started
    if failure_status is not None:
        _write_parent_failure(output_path, failure_status, failure_reason or "unknown", stdout, stderr)
    else:
        (output_path.parent / "stdout.log").write_text(stdout, encoding="utf-8")
        (output_path.parent / "stderr.log").write_text(stderr, encoding="utf-8")

    receipt = {
        "schema": "modeiii_solved_field_j_coupling_run_receipt.v1",
        "execution_authorized_by_cli_flag": True,
        "actual_argv": argv,
        "worker_return_code": child_return_code,
        "outer_wall_seconds": elapsed,
        "outer_timeout_seconds": WALL_TIMEOUT_SECONDS,
        "worker_timeout_seconds": WORKER_TIMEOUT_SECONDS,
        "timed_out": failure_status == "UNKNOWN_WALL_TIMEOUT",
        "result_path": str(output_path.resolve()),
        "result_bytes": output_path.stat().st_size if output_path.is_file() else None,
        "result_sha256": _sha256_path(output_path) if output_path.is_file() else None,
        "result_status": (
            _parse_json_without_nonfinite(output_path).get("status")
            if output_path.is_file()
            else failure_status
        ),
        "stdout_log": "stdout.log",
        "stderr_log": "stderr.log",
        "worker_staging_path": str(stage_path.resolve()),
        "worker_staging_preserved": stage_path.is_file(),
        "worker_staging_bytes": stage_path.stat().st_size if stage_path.is_file() else None,
        "worker_staging_sha256": _sha256_path(stage_path) if stage_path.is_file() else None,
        "failure_status": failure_status,
        "failure_reason": failure_reason,
    }
    (output_path.parent / "run_receipt.json").write_text(
        json.dumps(receipt, ensure_ascii=True, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return 0 if failure_status is None else 2


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--execute-readback", action="store_true")
    parser.add_argument("--worker-output", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker_output is not None:
        return _worker(args.worker_output.resolve())
    if not args.execute_readback:
        parser.error("source-only by default; obtain review approval, then pass --execute-readback")
    return _parent_run(args.output.resolve())


if __name__ == "__main__":
    raise SystemExit(main())
