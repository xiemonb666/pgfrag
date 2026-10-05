"""One bounded source-only Mode III PU scalar-operator diagnostic.

When separately reviewed and invoked, this runner evaluates continuous MLS
shapes, enriches them with the unshifted normalized F1 field, assembles one
restricted scalar Uz stiffness matrix, and writes JSON diagnostics. It does
not solve a boundary-value problem or qualify a crack mechanics model.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback


HERE = Path(__file__).resolve().parent
PGFRAG = HERE.parents[1]
MLS_DIR = PGFRAG / "experiments" / "mls3d_kernel_20261005"
OPERATOR_DIR = PGFRAG / "experiments" / "efgm_elastic_operator_20261005"
CONTRACT_DIR = PGFRAG / "experiments" / "enriched_dof_contract_20261005"
FIELD_DIR = PGFRAG / "experiments" / "modeiii_field_reference_20261005"
CASE_PATH = HERE / "fixed_case.json"

CASE_ID = "modeiii_pu_uz_operator_bridge_20261005"
MAX_SOURCE_COUNT = 125
MAX_TOTAL_DOFS = 375
MAX_QUADRATURE_POINTS = 256
MAX_NUMERIC_ARRAY_BYTES = 32 * 1024 * 1024
WALL_TIMEOUT_SECONDS = 30.0

# Root-reviewed dependency bytes. Any mismatch blocks before importing the
# numerical modules or allocating the fixture arrays.
EXPECTED_DEPENDENCIES = {
    "enriched_dof_contract.py": "e3ddb3b067465bb5b269576a1b8a5100e2bb92c9c924bc09c949702c19efa9d8",
    "field_reference.py": "7b07dfb30a64c50734dc18de2778b7725320de6d5a8f275c55f301c453ef18d7",
    "mls3d_shape.py": "66e680616a988182875d17aeb84fadb541867aee13fc15e0ae098f7c98839ef5",
    "mls3d_mls.py": "404256d9b7546361a5704c4c3cec7c022069a9cf83a137c58a16e7a64f79ec63",
    "efgm_elastic_operator.py": "41a5ac6ef6359ed1640ee681c94a3496c18652b29cc72c679b3e37c20dd9ac0d",
}


class DiagnosticBudgetError(RuntimeError):
    """Raised when fixed numeric-array or elapsed-time budgets are exceeded."""


class DiagnosticNumericsError(FloatingPointError):
    """Raised on nonfinite scalar stiffness or diagnostic values."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _timestamp() -> str:
    return datetime.now(timezone.utc).isoformat()


def _verify_dependencies() -> dict[str, str]:
    paths = {
        "enriched_dof_contract.py": CONTRACT_DIR / "enriched_dof_contract.py",
        "field_reference.py": FIELD_DIR / "field_reference.py",
        "mls3d_shape.py": MLS_DIR / "mls3d_shape.py",
        "mls3d_mls.py": MLS_DIR / "mls3d_mls.py",
        "efgm_elastic_operator.py": OPERATOR_DIR / "efgm_elastic_operator.py",
    }
    observed = {name: _sha256(path) for name, path in paths.items()}
    mismatches = {
        name: {"expected": EXPECTED_DEPENDENCIES[name], "actual": observed[name]}
        for name in paths
        if observed[name] != EXPECTED_DEPENDENCIES[name]
    }
    if mismatches:
        raise RuntimeError(f"read-only dependency source pin mismatch: {mismatches}")
    return observed


def _load_numerical_dependencies() -> None:
    """Import numeric modules only in the wall-timed worker process."""
    global np, AntiplaneShapeSample, BasisKind, DofKey, DofRegistry
    global FieldComponent, ProviderStatus, ScalarShapeContribution
    global antiplane_engineering_strain, validate_antiplane_sample, validate_problem_budget
    global IsotropicMaterial, MLS3DConfig, MLSStatus, shape_functions_3d
    global FieldStatus, ModeIIIFixedCase, evaluate_field, ContinuousMLSFrontProvider

    for path in (str(MLS_DIR), str(OPERATOR_DIR), str(CONTRACT_DIR), str(FIELD_DIR), str(HERE)):
        if path not in sys.path:
            sys.path.insert(0, path)
    import numpy as np
    from enriched_dof_contract import (
        AntiplaneShapeSample,
        BasisKind,
        DofKey,
        DofRegistry,
        FieldComponent,
        ProviderStatus,
        ScalarShapeContribution,
        antiplane_engineering_strain,
        validate_antiplane_sample,
        validate_problem_budget,
    )
    from efgm_elastic_operator import IsotropicMaterial
    from field_reference import FieldStatus, ModeIIIFixedCase, evaluate_field
    from mls3d_mls import MLS3DConfig, MLSStatus
    from mls3d_shape import shape_functions_3d
    from mls_modeiii_provider import ContinuousMLSFrontProvider


def _require_fixed_case(case: dict[str, object]) -> None:
    if case.get("case_id") != CASE_ID:
        raise ValueError("fixed case ID does not match this diagnostic")
    geometry = case["geometry"]
    sources = case["sources"]
    material = case["material_and_manufactured_field"]
    dofs = case["dofs"]
    mls = case["mls"]
    quadrature = case["quadrature"]
    limits = case["diagnostic"]
    assert all(isinstance(value, dict) for value in (geometry, sources, material, dofs, mls, quadrature, limits))
    volume = quadrature["volume"]
    side = quadrature["lateral_boundary"]
    caps = quadrature["end_caps"]
    rules = case["gauss_legendre_rules"]
    assert all(isinstance(value, dict) for value in (volume, side, caps, rules))
    expected = (
        geometry["domain"] == "x^2 + y^2 <= R^2, -L/2 <= z <= L/2",
        geometry["radius_m"] == 0.01,
        geometry["length_m"] == 0.01,
        geometry["z_bounds_m"] == [-0.005, 0.005],
        sources["grid_shape"] == [3, 3, 3],
        sources["axis_coordinates_m"] == [-0.005, 0.0, 0.005],
        sources["grid_spacing_m"] == 0.005,
        sources["source_count"] == 27,
        sources["support_radius_m"] == 0.025,
        sources["all_front_support_source_ids"] == "all 27 sources",
        sources["heaviside_set_J_H"] == [],
        dofs["coefficient_dof_count"] == 54,
        dofs["basis_kinds"] == ["BASE", "FRONT_F1_UNSHIFTED"],
        dofs["physical_units_for_every_coefficient"] == "m",
        material["young_modulus_Pa"] == 1_000_000.0,
        material["poisson_ratio"] == 0.3,
        material["K_III_Pa_sqrt_m"] == 100.0,
        mls["reference_origin_m"] == [0.0, 0.0, 0.0],
        mls["reference_scale_m"] == 0.01,
        mls["support_radius_m"] == 0.025,
        mls["visibility_mask"] is None,
        mls["mask_locally_constant"] is None,
        volume["point_count"] == 128,
        volume["radial_legendre_order"] == 4,
        volume["theta_legendre_order"] == 8,
        volume["z_legendre_order"] == 4,
        side["point_count"] == 64,
        caps["point_count_each"] == 32,
        quadrature["total_unique_query_count"] == 256,
        limits["maximum_source_count"] == MAX_SOURCE_COUNT,
        limits["maximum_total_dofs"] == MAX_TOTAL_DOFS,
        limits["maximum_quadrature_points"] == MAX_QUADRATURE_POINTS,
        limits["maximum_estimated_numeric_array_bytes"] == MAX_NUMERIC_ARRAY_BYTES,
        limits["wall_timeout_seconds"] == WALL_TIMEOUT_SECONDS,
        limits["run_solver"] is False,
        limits["apply_load_to_solver"] is False,
        limits["visibility_mask"] is False,
        limits["use_pseudoinverse_or_regularization"] is False,
        limits["mechanics_field_qualification"] is False,
    )
    if not all(expected):
        raise ValueError("fixed case changed a frozen geometry, MLS, material, quadrature, or budget input")
    for name, order in (("order_4", 4), ("order_8", 8), ("order_16", 16)):
        rule = rules[name]
        assert isinstance(rule, dict)
        if len(rule["nodes"]) != order or len(rule["weights"]) != order:
            raise ValueError(f"Gauss-Legendre rule {name} has an unexpected order")
        if not all(-1.0 < float(value) < 1.0 for value in rule["nodes"]):
            raise ValueError(f"Gauss-Legendre rule {name} contains an endpoint node")


def _array_budget_estimate(source_count: int, dof_count: int, query_count: int, volume_count: int) -> dict[str, int]:
    f8 = 8
    matrix_bytes = dof_count * dof_count * f8
    parts = {
        "source_coordinates": source_count * 3 * f8,
        "all_volume_and_boundary_queries": query_count * 3 * f8,
        "volume_and_boundary_quadrature_weights": (volume_count + query_count - volume_count) * f8,
        "boundary_normals": (query_count - volume_count) * 3 * f8,
        "quadrature_rule_and_intermediate_arrays": 1024 * f8,
        "MLS_weights": query_count * source_count * f8,
        "MLS_gradients": query_count * source_count * 3 * f8,
        "enriched_value_rows": query_count * dof_count * f8,
        "enriched_gradient_rows": query_count * dof_count * 3 * f8,
        "field_reference_values_tensors_and_boundary_tractions": query_count * 40 * f8,
        "constitutive_B_and_local_work": (6 * 6 + 6 * dof_count) * f8,
        "stiffness_and_spectral_numeric_work": 16 * matrix_bytes,
        "coefficient_rigid_eigenvalue_vectors": (dof_count * 4 + dof_count) * f8,
    }
    return {**parts, "estimated_explicit_numeric_array_bytes": int(sum(parts.values()))}


def _gauss_rule(case: dict[str, object], order: int) -> tuple[np.ndarray, np.ndarray]:
    rules = case["gauss_legendre_rules"]
    assert isinstance(rules, dict)
    rule = rules[f"order_{order}"]
    assert isinstance(rule, dict)
    return (
        np.asarray(rule["nodes"], dtype=np.float64),
        np.asarray(rule["weights"], dtype=np.float64),
    )


def _map_rule(nodes: np.ndarray, weights: np.ndarray, lower: float, upper: float) -> tuple[np.ndarray, np.ndarray]:
    half = 0.5 * (upper - lower)
    middle = 0.5 * (upper + lower)
    return middle + half * nodes, half * weights


def _build_sources(case: dict[str, object]) -> tuple[np.ndarray, tuple[str, ...]]:
    sources = case["sources"]
    assert isinstance(sources, dict)
    axis = np.asarray(sources["axis_coordinates_m"], dtype=np.float64)
    xyz = [(float(x), float(y), float(z)) for x, y, z in itertools.product(axis, axis, axis)]
    positions = np.ascontiguousarray(np.asarray(xyz, dtype=np.float64))
    source_ids = tuple(f"node-{row:03d}" for row in range(len(xyz)))
    return positions, source_ids


def _build_quadrature(case: dict[str, object]) -> dict[str, np.ndarray | int | float]:
    geometry = case["geometry"]
    rules = case["quadrature"]
    assert isinstance(geometry, dict) and isinstance(rules, dict)
    volume_record = rules["volume"]
    side_record = rules["lateral_boundary"]
    cap_record = rules["end_caps"]
    assert isinstance(volume_record, dict) and isinstance(side_record, dict) and isinstance(cap_record, dict)
    radius = float(geometry["radius_m"])
    length = float(geometry["length_m"])
    half_length = length / 2.0

    r4, wr4 = _gauss_rule(case, 4)
    theta8, wtheta8 = _gauss_rule(case, 8)
    z4, wz4 = _gauss_rule(case, 4)
    radial, radial_weights = _map_rule(r4, wr4, 0.0, radius)
    theta, theta_weights = _map_rule(theta8, wtheta8, -math.pi, math.pi)
    axial, axial_weights = _map_rule(z4, wz4, -half_length, half_length)
    volume_points: list[tuple[float, float, float]] = []
    volume_weights: list[float] = []
    for ir, it, iz in itertools.product(range(4), range(8), range(4)):
        r = float(radial[ir])
        angle = float(theta[it])
        volume_points.append((r * math.cos(angle), r * math.sin(angle), float(axial[iz])))
        volume_weights.append(float(radial_weights[ir] * theta_weights[it] * axial_weights[iz] * r))
    volume_xyz = np.ascontiguousarray(np.asarray(volume_points, dtype=np.float64))
    volume_w = np.ascontiguousarray(np.asarray(volume_weights, dtype=np.float64))

    theta16, wtheta16 = _gauss_rule(case, 16)
    side_theta, side_theta_weights = _map_rule(theta16, wtheta16, -math.pi, math.pi)
    side_axial, side_axial_weights = _map_rule(z4, wz4, -half_length, half_length)
    side_points: list[tuple[float, float, float]] = []
    side_normals: list[tuple[float, float, float]] = []
    side_area: list[float] = []
    for it, iz in itertools.product(range(16), range(4)):
        angle = float(side_theta[it])
        cosine, sine = math.cos(angle), math.sin(angle)
        side_points.append((radius * cosine, radius * sine, float(side_axial[iz])))
        side_normals.append((cosine, sine, 0.0))
        side_area.append(float(radius * side_theta_weights[it] * side_axial_weights[iz]))
    side_xyz = np.ascontiguousarray(np.asarray(side_points, dtype=np.float64))
    side_n = np.ascontiguousarray(np.asarray(side_normals, dtype=np.float64))
    side_w = np.ascontiguousarray(np.asarray(side_area, dtype=np.float64))

    cap_radial, cap_radial_weights = radial, radial_weights
    cap_theta, cap_theta_weights = theta, theta_weights
    cap_points: list[tuple[float, float, float]] = []
    cap_area: list[float] = []
    for ir, it in itertools.product(range(4), range(8)):
        r = float(cap_radial[ir])
        angle = float(cap_theta[it])
        cap_points.append((r * math.cos(angle), r * math.sin(angle), 0.0))
        cap_area.append(float(cap_radial_weights[ir] * cap_theta_weights[it] * r))
    cap_xy = np.asarray(cap_points, dtype=np.float64)
    cap_w = np.ascontiguousarray(np.asarray(cap_area, dtype=np.float64))
    top_xyz = np.ascontiguousarray(np.column_stack((cap_xy[:, 0], cap_xy[:, 1], np.full(32, half_length))))
    bottom_xyz = np.ascontiguousarray(np.column_stack((cap_xy[:, 0], cap_xy[:, 1], np.full(32, -half_length))))
    cap_top_n = np.ascontiguousarray(np.tile(np.array([0.0, 0.0, 1.0]), (32, 1)))
    cap_bottom_n = np.ascontiguousarray(np.tile(np.array([0.0, 0.0, -1.0]), (32, 1)))

    boundary_xyz = np.ascontiguousarray(np.vstack((side_xyz, top_xyz, bottom_xyz)))
    boundary_normals = np.ascontiguousarray(np.vstack((side_n, cap_top_n, cap_bottom_n)))
    boundary_weights = np.ascontiguousarray(np.concatenate((side_w, cap_w, cap_w)))
    query_xyz = np.ascontiguousarray(np.vstack((volume_xyz, boundary_xyz)))
    return {
        "volume_points": volume_xyz,
        "volume_weights": volume_w,
        "side_points": side_xyz,
        "side_normals": side_n,
        "side_area_weights": side_w,
        "top_points": top_xyz,
        "bottom_points": bottom_xyz,
        "cap_area_weights": cap_w,
        "boundary_points": boundary_xyz,
        "boundary_normals": boundary_normals,
        "boundary_area_weights": boundary_weights,
        "query_points": query_xyz,
        "r_min_sample_m": float(np.min(radial)),
        "r_max_sample_m": float(np.max(radial)),
        "volume_point_count": int(volume_xyz.shape[0]),
        "side_point_count": int(side_xyz.shape[0]),
        "cap_point_count_each": int(cap_w.shape[0]),
        "query_point_count": int(query_xyz.shape[0]),
    }


def _build_registry(source_ids: tuple[str, ...]) -> DofRegistry:
    registry = DofRegistry()
    for source_id in source_ids:
        registry.register(DofKey(source_id, FieldComponent.UZ, BasisKind.BASE, 0))
        registry.register(
            DofKey(source_id, FieldComponent.UZ, BasisKind.FRONT_F1_UNSHIFTED, 1)
        )
    registry.freeze()
    return registry


def _status_counts(statuses: tuple[object, ...]) -> dict[str, int]:
    all_statuses = sorted({status.value for status in statuses})
    return {name: sum(status.value == name for status in statuses) for name in all_statuses}


def _json_finite(value: float) -> float | None:
    return float(value) if math.isfinite(float(value)) else None


def _roundoff_references(volume_count: int, dof_count: int, maximum_condition: float) -> dict[str, float | int | str]:
    epsilon = float(np.finfo(np.float64).eps)
    unit_roundoff = epsilon / 2.0
    operation_count = 12 + volume_count + 2 * dof_count + dof_count * dof_count
    product = operation_count * unit_roundoff
    gamma = product / (1.0 - product)
    condition_factor = max(1.0, float(maximum_condition))
    return {
        "binary64_epsilon": epsilon,
        "unit_roundoff": unit_roundoff,
        "conservative_operation_count": operation_count,
        "gamma_n": gamma,
        "symmetry_relative_reference": 4.0 * gamma,
        "rigid_relative_reference_times_observed_condition": 8.0 * gamma * condition_factor,
        "field_relative_reference_times_observed_condition": 8.0 * gamma * condition_factor,
        "eigenvalue_zero_absolute_reference_multiplier": 4.0 * gamma,
        "interpretation": "roundoff-scale diagnostic references only; no qualification or pass/fail gate",
    }


def _blocked_result(
    case: dict[str, object],
    status: str,
    provider: object | None,
    budget: dict[str, int],
    started_at: str,
    started: float,
    dependency_hashes: dict[str, str],
    reason: str,
) -> dict[str, object]:
    mls_info: dict[str, object] | None = None
    if provider is not None:
        counts = _status_counts(provider.mls_result.statuses)
        finite = provider.mls_result.condition_numbers[np.isfinite(provider.mls_result.condition_numbers)]
        mls_info = {
            "status_counts": counts,
            "support_count_min": int(np.min(provider.mls_result.support_counts)),
            "support_count_max": int(np.max(provider.mls_result.support_counts)),
            "condition_min": float(np.min(finite)) if finite.size else None,
            "condition_max": float(np.max(finite)) if finite.size else None,
            "unknown_rows": provider.unknown_rows,
            "field_reference_status_counts": _status_counts(
                tuple(provider.field_result.status_at(i) for i in range(len(provider.field_result.statuses)))
            ),
        }
    return {
        "schema": "modeiii_pu_operator_diagnostic.v1",
        "case_id": CASE_ID,
        "status": status,
        "classification": "DIAGNOSTIC_ONLY_UNQUALIFIED",
        "started_at_utc": started_at,
        "fixed_case": case,
        "blocked_reason": reason,
        "mls": mls_info,
        "array_budget": {
            "parts_bytes": budget,
            "estimated_explicit_numeric_array_bytes": budget.get(
                "estimated_explicit_numeric_array_bytes", 0
            ),
            "cap_bytes": MAX_NUMERIC_ARRAY_BYTES,
            "python_object_overhead_and_blas_workspace_excluded": True,
        },
        "no_partial_stiffness_returned": True,
        "source_sha256": dependency_hashes,
        "elapsed_seconds": time.monotonic() - started,
        "finished_at_utc": _timestamp(),
        "limitations": case["limitations"],
    }


def _run_diagnostic() -> dict[str, object]:
    started = time.monotonic()
    started_at = _timestamp()
    dependency_hashes = _verify_dependencies()
    case = json.loads(CASE_PATH.read_text(encoding="utf-8"))
    _require_fixed_case(case)
    sources = case["sources"]
    dofs = case["dofs"]
    quadrature_case = case["quadrature"]
    material_case = case["material_and_manufactured_field"]
    mls_case = case["mls"]
    assert isinstance(sources, dict) and isinstance(dofs, dict)
    assert isinstance(quadrature_case, dict) and isinstance(material_case, dict)
    assert isinstance(mls_case, dict)
    source_count = int(sources["source_count"])
    dof_count = int(dofs["coefficient_dof_count"])
    volume_count = int(quadrature_case["volume"]["point_count"])
    query_count = int(quadrature_case["total_unique_query_count"])
    validate_problem_budget(source_count, dof_count, query_count)
    budget = _array_budget_estimate(source_count, dof_count, query_count, volume_count)
    if int(budget["estimated_explicit_numeric_array_bytes"]) > MAX_NUMERIC_ARRAY_BYTES:
        raise DiagnosticBudgetError("explicit numeric array estimate exceeds the 32 MiB cap")

    # The preflight above precedes all fixture, field, MLS, and operator arrays.
    source_positions, source_ids = _build_sources(case)
    quadrature = _build_quadrature(case)
    query_positions = quadrature["query_points"]
    if not isinstance(query_positions, np.ndarray):
        raise RuntimeError("fixed quadrature builder did not return a query array")
    if (
        source_positions.shape != (source_count, 3)
        or query_positions.shape != (query_count, 3)
        or int(quadrature["volume_point_count"]) != volume_count
        or int(quadrature["side_point_count"]) != 64
        or int(quadrature["cap_point_count_each"]) != 32
    ):
        raise ValueError("fixed source or quadrature builder returned unexpected dimensions")

    registry = _build_registry(source_ids)
    if registry.dof_count != dof_count or registry.source_count != source_count:
        raise RuntimeError("frozen registry dimensions differ from the fixed case")
    field_case = ModeIIIFixedCase(
        radius_m=float(case["geometry"]["radius_m"]),
        length_m=float(case["geometry"]["length_m"]),
        K_III_Pa_sqrt_m=float(material_case["K_III_Pa_sqrt_m"]),
        young_modulus_Pa=float(material_case["young_modulus_Pa"]),
        poisson_ratio=float(material_case["poisson_ratio"]),
    )
    mls_config = MLS3DConfig(
        reference_origin=tuple(float(value) for value in mls_case["reference_origin_m"]),
        reference_scale=float(mls_case["reference_scale_m"]),
        support_radius=float(mls_case["support_radius_m"]),
        max_condition_number=float(mls_case["max_condition_number"]),
        rank_relative_tolerance=float(mls_case["rank_relative_tolerance"]),
        switch_tolerance=float(mls_case["switch_tolerance"]),
        duplicate_tolerance=float(mls_case["duplicate_tolerance"]),
    )
    provider = ContinuousMLSFrontProvider(
        source_positions,
        source_ids,
        query_positions,
        registry,
        mls_config,
        field_case,
    )
    if time.monotonic() - started > WALL_TIMEOUT_SECONDS:
        raise DiagnosticBudgetError("30-second limit exceeded during MLS/reference evaluation")

    if provider.unknown_rows:
        result = _blocked_result(
            case,
            "BLOCKED_UNKNOWN_SUPPORT_OR_REFERENCE",
            provider,
            budget,
            started_at,
            started,
            dependency_hashes,
            "at least one MLS or analytic-reference query is UNKNOWN; no stiffness assembled",
        )
        result["dimensions"] = {
            "source_count": source_count,
            "total_scalar_dof_count": dof_count,
            "volume_quadrature_count": volume_count,
            "all_unique_query_count": query_count,
        }
        return result

    material = IsotropicMaterial(
        young_modulus=float(material_case["young_modulus_Pa"]),
        poisson_ratio=float(material_case["poisson_ratio"]),
    )
    constitutive = material.constitutive_matrix
    stiffness = np.zeros((dof_count, dof_count), dtype=np.float64)
    volume_weights = quadrature["volume_weights"]
    assert isinstance(volume_weights, np.ndarray)
    for query_id in range(volume_count):
        sample = provider.evaluate(tuple(float(value) for value in query_positions[query_id]))
        if sample.status is not ProviderStatus.OK:
            result = _blocked_result(
                case,
                "BLOCKED_PROVIDER_STATUS",
                provider,
                budget,
                started_at,
                started,
                dependency_hashes,
                sample.diagnostic or sample.status.value,
            )
            result["blocked_quadrature_index"] = query_id
            result["no_partial_stiffness_returned"] = True
            return result
        b = np.zeros((6, dof_count), dtype=np.float64)
        for coefficient_id, contribution in enumerate(sample.contributions):
            b[:, coefficient_id] = np.asarray(
                antiplane_engineering_strain(contribution.gradient), dtype=np.float64
            )
        with np.errstate(over="ignore", invalid="ignore"):
            local = (b.T @ constitutive @ b) * float(volume_weights[query_id])
            stiffness += local
        if not np.all(np.isfinite(local)) or not np.all(np.isfinite(stiffness)):
            raise DiagnosticNumericsError(
                f"scalar stiffness became nonfinite at volume quadrature row {query_id}"
            )
    if provider.consumed_queries != volume_count:
        raise RuntimeError("volume assembly did not consume its exact query prefix")
    stiffness.setflags(write=False)

    coefficients = np.zeros(dof_count, dtype=np.float64)
    front_coefficient = field_case.amplitude_sqrt_m * math.sqrt(field_case.radius_m)
    for slot, key in enumerate(registry.ordered_keys):
        if key.basis_kind is BasisKind.BASE:
            coefficients[slot] = 0.0
        elif key.basis_kind is BasisKind.FRONT_F1_UNSHIFTED:
            coefficients[slot] = front_coefficient
        else:
            raise RuntimeError(f"unexpected registered basis {key.basis_kind.value}")

    # Finish the fixed-order boundary batch through the same validated provider.
    boundary_values = np.empty(query_count - volume_count, dtype=np.float64)
    for boundary_row, query_id in enumerate(range(volume_count, query_count)):
        sample = provider.evaluate(tuple(float(value) for value in query_positions[query_id]))
        if sample.status is not ProviderStatus.OK:
            result = _blocked_result(
                case,
                "BLOCKED_PROVIDER_STATUS",
                provider,
                budget,
                started_at,
                started,
                dependency_hashes,
                sample.diagnostic or sample.status.value,
            )
            result["blocked_query_index"] = query_id
            result["stiffness_discarded"] = True
            result["no_partial_stiffness_returned"] = True
            return result
        boundary_values[boundary_row] = sum(
            coefficients[index] * contribution.value
            for index, contribution in enumerate(sample.contributions)
        )
    if provider.consumed_queries != query_count:
        raise RuntimeError("provider did not consume the complete frozen volume/boundary batch")

    exact_field = provider.field_result
    values = provider.values
    gradients = provider.gradients
    if values is None or gradients is None:
        raise RuntimeError("provider reported success without numerical rows")
    displacement_h = values @ coefficients
    gradient_h = np.einsum("qdi,d->qi", gradients, coefficients)
    exact_uz = exact_field.displacement_m[:, 2]
    exact_grad_uz = exact_field.displacement_gradient[:, 2, :]
    displacement_error = np.abs(displacement_h - exact_uz)
    gradient_error = np.abs(gradient_h - exact_grad_uz)
    displacement_scale = max(float(np.max(np.abs(exact_uz))), np.finfo(np.float64).tiny)
    gradient_scale = max(float(np.max(np.abs(exact_grad_uz))), np.finfo(np.float64).tiny)

    volume_points = quadrature["volume_points"]
    assert isinstance(volume_points, np.ndarray)
    _ = volume_points  # Fixed volume prefix is represented in query_positions.
    exact_volume_energy = float(
        np.dot(exact_field.energy_density_J_per_m3[:volume_count], volume_weights)
    )
    analytic_energy = float(field_case.exact_total_strain_energy_J)
    coefficient_energy = 0.5 * float(coefficients @ stiffness @ coefficients)

    boundary_points = quadrature["boundary_points"]
    boundary_normals = quadrature["boundary_normals"]
    boundary_area_weights = quadrature["boundary_area_weights"]
    assert isinstance(boundary_points, np.ndarray)
    assert isinstance(boundary_normals, np.ndarray)
    assert isinstance(boundary_area_weights, np.ndarray)
    boundary_field = evaluate_field(boundary_points, case=field_case, normals=boundary_normals)
    boundary_statuses = [boundary_field.status_at(i).value for i in range(boundary_points.shape[0])]
    if any(status != FieldStatus.OK.value for status in boundary_statuses):
        raise RuntimeError("fixed boundary quadrature unexpectedly hit a field UNKNOWN/outside status")
    traction = boundary_field.traction_Pa
    if traction is None:
        raise RuntimeError("boundary field reference omitted requested traction values")
    exact_boundary_uz = boundary_field.displacement_m[:, 2]
    exact_boundary_work = float(
        np.dot(exact_boundary_uz * traction[:, 2], boundary_area_weights)
    )
    discrete_boundary_work = float(
        np.dot(boundary_values * traction[:, 2], boundary_area_weights)
    )

    side_count = int(quadrature["side_point_count"])
    cap_count = int(quadrature["cap_point_count_each"])
    side_slice = slice(0, side_count)
    top_slice = slice(side_count, side_count + cap_count)
    bottom_slice = slice(side_count + cap_count, side_count + 2 * cap_count)
    side_work = float(
        np.dot(
            boundary_values[side_slice] * traction[side_slice, 2],
            boundary_area_weights[side_slice],
        )
    )
    top_work = float(
        np.dot(
            boundary_values[top_slice] * traction[top_slice, 2],
            boundary_area_weights[top_slice],
        )
    )
    bottom_work = float(
        np.dot(
            boundary_values[bottom_slice] * traction[bottom_slice, 2],
            boundary_area_weights[bottom_slice],
        )
    )
    cap_records: dict[str, object] = {}
    for name, region_slice in (("top", top_slice), ("bottom", bottom_slice)):
        region_traction = traction[region_slice]
        region_weights = boundary_area_weights[region_slice]
        cap_records[name] = {
            "max_abs_full_traction_by_xyz_Pa": [
                float(np.max(np.abs(region_traction[:, component]))) for component in range(3)
            ],
            "integrated_full_traction_resultant_N": [
                float(value) for value in np.sum(region_traction * region_weights[:, None], axis=0)
            ],
            "max_abs_scalar_tz_Pa": float(np.max(np.abs(region_traction[:, 2]))),
            "uz_virtual_work_J": float(
                np.dot(
                    boundary_values[region_slice] * region_traction[:, 2],
                    region_weights,
                )
            ),
            "x_y_full_traction_recorded_but_not_in_scalar_uz_virtual_work": True,
        }

    symmetry_relative = float(
        np.linalg.norm(stiffness - stiffness.T, ord="fro")
        / max(np.linalg.norm(stiffness, ord="fro"), np.finfo(np.float64).tiny)
    )
    stiffness_norm = float(np.linalg.norm(stiffness, ord="fro"))
    translation = np.zeros(dof_count, dtype=np.float64)
    for slot, key in enumerate(registry.ordered_keys):
        if key.basis_kind is BasisKind.BASE:
            translation[slot] = 1.0
    rigid_residual = float(
        np.linalg.norm(stiffness @ translation)
        / max(stiffness_norm * np.linalg.norm(translation), np.finfo(np.float64).tiny)
    )

    symmetric_part = 0.5 * (stiffness + stiffness.T)
    eigenvalues = np.linalg.eigvalsh(symmetric_part)
    if not np.all(np.isfinite(eigenvalues)):
        raise DiagnosticNumericsError("symmetric-part spectrum contains nonfinite eigenvalues")
    condition_numbers = provider.mls_result.condition_numbers
    finite_conditions = condition_numbers[np.isfinite(condition_numbers)]
    maximum_condition = float(np.max(finite_conditions)) if finite_conditions.size else 1.0
    minimum_condition = float(np.min(finite_conditions)) if finite_conditions.size else 1.0
    roundoff = _roundoff_references(volume_count, dof_count, maximum_condition)
    eigenvalue_scale = max(float(np.max(np.abs(eigenvalues))), np.finfo(np.float64).tiny)
    eigen_tolerance = float(roundoff["eigenvalue_zero_absolute_reference_multiplier"]) * eigenvalue_scale
    near_zero_count = int(np.count_nonzero(np.abs(eigenvalues) <= eigen_tolerance))
    negative_count = int(np.count_nonzero(eigenvalues < 0.0))
    negative_beyond_roundoff = int(np.count_nonzero(eigenvalues < -eigen_tolerance))
    max_source_to_front_distance = max(
        math.sqrt(
            float(source_positions[row, 0] ** 2 + source_positions[row, 1] ** 2)
            + max(
                abs(float(source_positions[row, 2]) - float(case["geometry"]["z_bounds_m"][0])),
                abs(float(source_positions[row, 2]) - float(case["geometry"]["z_bounds_m"][1])),
            ) ** 2
        )
        for row in range(source_count)
    )

    result: dict[str, object] = {
        "schema": "modeiii_pu_operator_diagnostic.v1",
        "case_id": CASE_ID,
        "status": "ASSEMBLED_DIAGNOSTIC_ONLY_UNQUALIFIED",
        "classification": "ENGINEERING_MANUFACTURED_FIELD_OPERATOR_DIAGNOSTIC",
        "started_at_utc": started_at,
        "fixed_case": case,
        "dimensions": {
            "source_count": source_count,
            "total_scalar_dof_count": dof_count,
            "volume_quadrature_count": volume_count,
            "lateral_boundary_point_count": side_count,
            "endcap_point_count_each": cap_count,
            "unique_MLS_query_count": query_count,
            "K_shape": list(stiffness.shape),
            "K_units": "N/m",
        },
        "registry": {
            "state": registry.state.value,
            "source_count": registry.source_count,
            "dof_count": registry.dof_count,
            "basis_key_counts": {
                kind.value: sum(key.basis_kind is kind for key in registry.ordered_keys)
                for kind in BasisKind
            },
            "all_field_components_UZ": all(
                key.field_component is FieldComponent.UZ for key in registry.ordered_keys
            ),
            "canonical_ids": list(registry.canonical_ids),
            "front_support_set_T_F": list(source_ids),
            "heaviside_support_set_J_H": [],
            "topology_note": "all 27 source-centered h-balls contain the straight front; this is a global algebraic bridge only, not local enrichment qualification",
            "maximum_source_to_front_segment_distance_m": max_source_to_front_distance,
        },
        "basis_and_coefficients": {
            "base_formula": "phi_i",
            "front_formula": "psi_i=phi_i*Fbar, Fbar=F1/sqrt(R)",
            "front_gradient_formula": "Fbar*grad(phi_i)+phi_i*grad(Fbar)",
            "coefficient_units": "m for all BASE and FRONT_F1_UNSHIFTED keys",
            "front_factor_dimensionless": True,
            "front_factor_source_value_used": False,
            "base_coefficients_m": [0.0] * source_count,
            "front_coefficient_each_m": front_coefficient,
            "front_coefficients_m": [front_coefficient] * source_count,
            "manufactured_coefficient_vector": coefficients.tolist(),
        },
        "constitutive_reduction": {
            "strain_order": "[eps_xx, eps_yy, eps_zz, gamma_xy, gamma_yz, gamma_xz]",
            "uz_strain_column": "[0, 0, dphi/dz, 0, dphi/dy, dphi/dx]",
            "scalar_energy_density": "0.5*((lambda+2*mu)*(duz/dz)^2 + mu*(duz/dx)^2 + mu*(duz/dy)^2)",
            "longitudinal_term_included": True,
            "constitutive_lambda_Pa": float(material.lame_lambda),
            "constitutive_mu_Pa": float(material.shear_modulus),
            "D_source": "frozen IsotropicMaterial.constitutive_matrix; no change to shared operator",
        },
        "mls": {
            "api": "shape_functions_3d",
            "visibility": "continuous/no mask; visibility kwargs omitted",
            "row_order": "source row i maps to source_id node-iii and registry BASE/FRONT keys for that site",
            "status_counts": _status_counts(provider.mls_result.statuses),
            "support_count_min": int(np.min(provider.mls_result.support_counts)),
            "support_count_max": int(np.max(provider.mls_result.support_counts)),
            "condition_number_min": minimum_condition,
            "condition_number_max": maximum_condition,
            "condition_numbers": [_json_finite(float(value)) for value in condition_numbers],
            "partition_sum_max_abs_error": float(
                np.max(np.abs(np.sum(provider.mls_result.weights, axis=1) - 1.0))
            ),
            "partition_gradient_sum_max_norm_per_m": float(
                np.max(np.linalg.norm(np.sum(provider.mls_result.gradients, axis=1), axis=1))
            ),
            "unknown_rows": [],
        },
        "quadrature": {
            "query_order": ["volume", "lateral", "top_cap", "bottom_cap"],
            "volume_measure_sum_m3": float(np.sum(volume_weights)),
            "analytic_cylinder_volume_m3": math.pi * float(case["geometry"]["radius_m"]) ** 2 * float(case["geometry"]["length_m"]),
            "radial_domain_m": [0.0, float(case["geometry"]["radius_m"])],
            "r_min_positive_volume_sample_m": float(quadrature["r_min_sample_m"]),
            "r_max_volume_sample_m": float(quadrature["r_max_sample_m"]),
            "radial_hole_or_clipping": False,
            "radial_rule_integrates_full_0_to_R_interval": True,
            "angular_open_branch_interval_rad": [-math.pi, math.pi],
            "branch_cut_or_front_sampled": False,
            "lateral_area_sum_m2": float(np.sum(quadrature["side_area_weights"])) ,
            "top_cap_area_sum_m2": float(np.sum(quadrature["cap_area_weights"])) ,
            "bottom_cap_area_sum_m2": float(np.sum(quadrature["cap_area_weights"])) ,
        },
        "manufactured_modeiii_field": {
            "status": "DIAGNOSTIC_ONLY",
            "displacement_max_absolute_error_m": float(np.max(displacement_error)),
            "displacement_max_relative_error": float(np.max(displacement_error)) / displacement_scale,
            "gradient_max_absolute_error_per_m": float(np.max(gradient_error)),
            "gradient_max_relative_error": float(np.max(gradient_error)) / gradient_scale,
            "reference_is_acceptance_gate": False,
            "interpretation": "global partition reproduction on T_F=all sources is algebraic evidence only",
        },
        "energy_and_work": {
            "exact_closed_form_strain_energy_J": analytic_energy,
            "independent_exact_tensor_volume_energy_quadrature_J": exact_volume_energy,
            "manufactured_K_quadratic_energy_J": coefficient_energy,
            "K_energy_vs_exact_tensor_quadrature_relative_difference": abs(coefficient_energy - exact_volume_energy) / max(abs(exact_volume_energy), np.finfo(np.float64).tiny),
            "tensor_quadrature_vs_closed_form_relative_difference": abs(exact_volume_energy - analytic_energy) / max(abs(analytic_energy), np.finfo(np.float64).tiny),
            "exact_uz_external_traction_work_J": exact_boundary_work,
            "discrete_uz_external_traction_work_J": discrete_boundary_work,
            "half_discrete_boundary_work_J": 0.5 * discrete_boundary_work,
            "side_scalar_uz_virtual_work_J": side_work,
            "top_cap_scalar_uz_virtual_work_J": top_work,
            "bottom_cap_scalar_uz_virtual_work_J": bottom_work,
            "top_and_bottom_full_traction_records": cap_records,
            "work_is_applied_to_solver": False,
            "all_energy_and_work_comparisons_are_qualification_gates": False,
        },
        "operator_diagnostics": {
            "matrix_frobenius_norm_N_per_m": float(np.linalg.norm(stiffness, ord="fro")),
            "stiffness_sha256_float64_c_order": hashlib.sha256(
                np.asarray(stiffness, dtype="<f8", order="C").tobytes(order="C")
            ).hexdigest(),
            "stiffness_matrix_N_per_m": stiffness.tolist(),
            "relative_frobenius_asymmetry": symmetry_relative,
            "spectrum_matrix": "symmetric_part=(K+K.T)/2 for eigvalsh only; raw assembled K remains unchanged",
            "eigenvalues_N_per_m": [float(value) for value in eigenvalues],
            "negative_eigenvalue_count": negative_count,
            "negative_beyond_roundoff_count": negative_beyond_roundoff,
            "near_zero_count": near_zero_count,
            "expected_scalar_rigid_mode_count": 1,
            "additional_near_zero_after_one_scalar_translation": max(near_zero_count - 1, 0),
            "near_zero_absolute_threshold_N_per_m": eigen_tolerance,
            "scalar_rigid_mode": "uniform uz translation: BASE coefficients=1, FRONT coefficients=0",
            "scalar_rigid_mode_scaled_residual": rigid_residual,
            "full_six_vector_rigid_modes_representable": False,
            "interpretation": "diagnostic spectrum only; no rank truncation, regularization, positive-definiteness enforcement, or stability claim",
        },
        "array_budget": {
            "parts_bytes": budget,
            "estimated_explicit_numeric_array_bytes": int(budget["estimated_explicit_numeric_array_bytes"]),
            "cap_bytes": MAX_NUMERIC_ARRAY_BYTES,
            "within_estimate_cap": int(budget["estimated_explicit_numeric_array_bytes"]) <= MAX_NUMERIC_ARRAY_BYTES,
            "python_object_overhead_and_blas_workspace_excluded": True,
        },
        "limitations": case["limitations"],
        "source_sha256": dependency_hashes,
        "elapsed_seconds": time.monotonic() - started,
        "wall_timeout_seconds": WALL_TIMEOUT_SECONDS,
        "within_wall_timeout": time.monotonic() - started <= WALL_TIMEOUT_SECONDS,
        "finished_at_utc": _timestamp(),
    }
    if not np.all(np.isfinite(stiffness)):
        raise DiagnosticNumericsError("assembled stiffness contains a nonfinite entry")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=HERE / "results" / "modeiii_pu_operator_diagnostic.json",
        help="JSON output path; does not modify the frozen case",
    )
    parser.add_argument("--_worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    output_path = args.output
    if args._worker:
        try:
            dependency_hashes = _verify_dependencies()
            _load_numerical_dependencies()
            result = _run_diagnostic()
            result["worker_dependency_sha256"] = dependency_hashes
            exit_code = 0 if result["status"] == "ASSEMBLED_DIAGNOSTIC_ONLY_UNQUALIFIED" else 2
        except Exception as exc:
            result = {
                "schema": "modeiii_pu_operator_diagnostic.v1",
                "case_id": CASE_ID,
                "status": "DIAGNOSTIC_ERROR",
                "classification": "DIAGNOSTIC_ONLY_UNQUALIFIED",
                "error_type": type(exc).__name__,
                "error": str(exc),
                "traceback": traceback.format_exc(limit=8),
                "no_retry": True,
                "finished_at_utc": _timestamp(),
            }
            exit_code = 2
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(
            json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        print(json.dumps(result, ensure_ascii=False, allow_nan=False))
        return exit_code

    # The parent enforces a true wall deadline even if a native numeric call
    # blocks. Numerical modules load only in the child and are inside timeout.
    output_path.parent.mkdir(parents=True, exist_ok=True)
    staging_path = output_path.with_name(f".{output_path.name}.{os.getpid()}.worker.json")
    staging_path.unlink(missing_ok=True)
    parent_started = time.monotonic()
    parent_started_at = _timestamp()
    try:
        completed = subprocess.run(
            [
                sys.executable,
                "-B",
                str(HERE / "modeiii_operator_diagnostic.py"),
                "--_worker",
                "--output",
                str(staging_path),
            ],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=WALL_TIMEOUT_SECONDS,
            check=False,
        )
        try:
            result = json.loads(staging_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            result = {
                "schema": "modeiii_pu_operator_diagnostic.v1",
                "case_id": CASE_ID,
                "status": "WORKER_FAILED",
                "worker_return_code": completed.returncode,
                "worker_stdout_tail": completed.stdout[-2000:],
                "worker_stderr_tail": completed.stderr[-2000:],
                "no_retry": True,
                "finished_at_utc": _timestamp(),
            }
            exit_code = 2
        else:
            exit_code = completed.returncode
    except subprocess.TimeoutExpired:
        result = {
            "schema": "modeiii_pu_operator_diagnostic.v1",
            "case_id": CASE_ID,
            "status": "WALL_TIMEOUT",
            "wall_timeout_seconds": WALL_TIMEOUT_SECONDS,
            "elapsed_seconds": time.monotonic() - parent_started,
            "no_retry": True,
            "no_result_published": True,
            "finished_at_utc": _timestamp(),
        }
        exit_code = 2
    finally:
        staging_path.unlink(missing_ok=True)
    parent_elapsed = time.monotonic() - parent_started
    result["parent_started_at_utc"] = parent_started_at
    result["parent_elapsed_seconds"] = parent_elapsed
    result["parent_wall_timeout_seconds"] = WALL_TIMEOUT_SECONDS
    result["parent_within_wall_timeout"] = parent_elapsed <= WALL_TIMEOUT_SECONDS
    output_path.write_text(
        json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(result, ensure_ascii=False, allow_nan=False))
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
