"""Bounded source-only antiplane traction solve for the straight-front Mode III fixture.

The runner reuses the pinned continuous MLS/front-PU provider and canonical
DOF registry. It assembles only the scalar UZ Galerkin K, side traction f, and
normalized physical mean constraint. It does not implement a new MLS or
enrichment primitive, a full-vector solver, crack growth, or plotting.
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
from typing import Any


HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[2]
PGFRAG = PROJECT_ROOT / "pgfrag"
MLS_DIR = PGFRAG / "experiments" / "mls3d_kernel_20261005"
OPERATOR_DIR = PGFRAG / "experiments" / "efgm_elastic_operator_20261005"
CONTRACT_DIR = PGFRAG / "experiments" / "enriched_dof_contract_20261005"
FIELD_DIR = PGFRAG / "experiments" / "modeiii_field_reference_20261005"
PU_OPERATOR_DIR = PGFRAG / "experiments" / "modeiii_pu_operator_20261005"
CASE_PATH = HERE / "static_solve_case.json"
BASE_CASE_PATH = PU_OPERATOR_DIR / "fixed_case.json"

CASE_ID = "modeiii_pu_uz_static_traction_solve_20261006_v3"
MAX_SOURCE_COUNT = 125
MAX_TOTAL_DOFS = 375
MAX_PROVIDER_QUERIES = 256
MAX_NUMERIC_ARRAY_BYTES = 32 * 1024 * 1024
MAX_RESULT_JSON_BYTES = 4 * 1024 * 1024
WALL_TIMEOUT_SECONDS = 30.0

# Direct byte pins. These are checked in the worker before importing NumPy or
# any numerical module. Source pins also appear in the sibling review freeze.
EXPECTED_SHA256 = {
    "static_solve_case.json": "13ca46e05b520624d903be5d71c55e58b4d23db621423b47da4b2a0c174a38a9",
    "base_operator_case.json": "dce0a1f3009246a52198bc0ed833ed67998144549a653e7bce9bcc17babd71b9",
    "mls_modeiii_provider.py": "6025642fdd119a95c76d693d3edbc4785d982731560e556ce7b13bad109c6034",
    "enriched_dof_contract.py": "e3ddb3b067465bb5b269576a1b8a5100e2bb92c9c924bc09c949702c19efa9d8",
    "field_reference.py": "7b07dfb30a64c50734dc18de2778b7725320de6d5a8f275c55f301c453ef18d7",
    "mls3d_shape.py": "66e680616a988182875d17aeb84fadb541867aee13fc15e0ae098f7c98839ef5",
    "mls3d_mls.py": "404256d9b7546361a5704c4c3cec7c022069a9cf83a137c58a16e7a64f79ec63",
    "efgm_elastic_operator.py": "41a5ac6ef6359ed1640ee681c94a3496c18652b29cc72c679b3e37c20dd9ac0d",
}

DEPENDENCY_PATHS = {
    "static_solve_case.json": CASE_PATH,
    "base_operator_case.json": BASE_CASE_PATH,
    "mls_modeiii_provider.py": PU_OPERATOR_DIR / "mls_modeiii_provider.py",
    "enriched_dof_contract.py": CONTRACT_DIR / "enriched_dof_contract.py",
    "field_reference.py": FIELD_DIR / "field_reference.py",
    "mls3d_shape.py": MLS_DIR / "mls3d_shape.py",
    "mls3d_mls.py": MLS_DIR / "mls3d_mls.py",
    "efgm_elastic_operator.py": OPERATOR_DIR / "efgm_elastic_operator.py",
}

np: Any = None
AntiplaneShapeSample: Any = None
BasisKind: Any = None
DofKey: Any = None
DofRegistry: Any = None
FieldComponent: Any = None
ProviderStatus: Any = None
antiplane_engineering_strain: Any = None
validate_problem_budget: Any = None
IsotropicMaterial: Any = None
MLS3DConfig: Any = None
ModeIIIFixedCase: Any = None
FieldStatus: Any = None
evaluate_field: Any = None
ContinuousMLSFrontProvider: Any = None


class StaticSolveUnknown(RuntimeError):
    """Raised only after an UNKNOWN result JSON has been written."""


def _timestamp() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _verify_source_pins() -> dict[str, str]:
    observed: dict[str, str] = {}
    mismatches: dict[str, dict[str, str]] = {}
    for name, path in DEPENDENCY_PATHS.items():
        if not path.is_file():
            observed[name] = "MISSING"
            mismatches[name] = {"expected": EXPECTED_SHA256[name], "actual": "MISSING"}
            continue
        actual = _sha256(path)
        observed[name] = actual
        if actual != EXPECTED_SHA256[name]:
            mismatches[name] = {"expected": EXPECTED_SHA256[name], "actual": actual}
    if mismatches:
        raise ValueError(f"source pin mismatch; imports are blocked: {mismatches}")
    return observed


def _read_and_validate_cases() -> tuple[dict[str, object], dict[str, object]]:
    solve_case = json.loads(CASE_PATH.read_text(encoding="utf-8"))
    base_case = json.loads(BASE_CASE_PATH.read_text(encoding="utf-8"))
    if solve_case.get("schema") != "modeiii_pu_antiplane_static_solve_case.v3":
        raise ValueError("static solve case schema changed")
    if solve_case.get("case_id") != CASE_ID:
        raise ValueError("static solve case ID changed")
    base_link = solve_case.get("base_operator_case")
    if not isinstance(base_link, dict) or base_link.get("case_id") != base_case.get("case_id"):
        raise ValueError("static solve no longer refers to the pinned base operator case")
    if base_case.get("case_id") != "modeiii_pu_uz_operator_bridge_20261005":
        raise ValueError("base operator case ID changed")

    geometry = base_case["geometry"]
    sources = base_case["sources"]
    material = base_case["material_and_manufactured_field"]
    mls = base_case["mls"]
    dofs = base_case["dofs"]
    if not all(isinstance(value, dict) for value in (geometry, sources, material, mls, dofs)):
        raise ValueError("pinned base case is missing a required object")
    expected_base = (
        geometry["radius_m"] == 0.01,
        geometry["length_m"] == 0.01,
        geometry["z_bounds_m"] == [-0.005, 0.005],
        sources["grid_shape"] == [3, 3, 3],
        sources["axis_coordinates_m"] == [-0.005, 0.0, 0.005],
        sources["grid_spacing_m"] == 0.005,
        sources["source_count"] == 27,
        sources["support_radius_m"] == 0.025,
        sources["front_set_T_F"] == "all 27 sources; this is only a global algebraic bridge because every fixed support contains the straight front",
        sources["heaviside_set_J_H"] == [],
        material["young_modulus_Pa"] == 1_000_000.0,
        material["poisson_ratio"] == 0.3,
        material["K_III_Pa_sqrt_m"] == 100.0,
        dofs["coefficient_dof_count"] == 54,
        dofs["basis_kinds"] == ["BASE", "FRONT_F1_UNSHIFTED"],
        mls["reference_origin_m"] == [0.0, 0.0, 0.0],
        mls["reference_scale_m"] == 0.01,
        mls["support_radius_m"] == 0.025,
        mls["visibility_mask"] is None,
        mls["mask_locally_constant"] is None,
        base_case["diagnostic"]["maximum_source_count"] == MAX_SOURCE_COUNT,
        base_case["diagnostic"]["maximum_total_dofs"] == MAX_TOTAL_DOFS,
        base_case["diagnostic"]["maximum_quadrature_points"] == MAX_PROVIDER_QUERIES,
        base_case["diagnostic"]["maximum_estimated_numeric_array_bytes"] == MAX_NUMERIC_ARRAY_BYTES,
    )
    if not all(expected_base):
        raise ValueError("pinned geometry, MLS scales, material, registry, or budget differs from the reviewed case")

    quadrature = solve_case["quadrature"]
    levels = quadrature["volume_levels"]
    budget = solve_case["budget"]
    if not isinstance(levels, list) or len(levels) != 2:
        raise ValueError("exactly two frozen volume quadrature levels are required")
    expected_levels = [(4, 8, 4, 128), (4, 16, 4, 256)]
    for level, expected in zip(levels, expected_levels):
        if not isinstance(level, dict):
            raise ValueError("each quadrature level must be an object")
        actual = (
            level.get("radial_GL"),
            level.get("theta_GL"),
            level.get("z_GL"),
            level.get("volume_point_count"),
        )
        if actual != expected:
            raise ValueError(f"volume quadrature changed from the frozen rule: {actual}")
    comparison = solve_case["quadrature_comparison"]
    if (
        comparison["refined_dimension"] != "theta angular quadrature"
        or comparison["level_128"] != "GL4 radial x GL8 theta x GL4 z"
        or comparison["level_256"] != "GL4 radial x GL16 theta x GL4 z"
        or comparison["fixed_dimensions"] != ["radial GL order 4", "axial GL order 4"]
    ):
        raise ValueError("the frozen volume-level comparison must refine theta only")
    if (
        quadrature["lateral_load_rule"]["theta_GL"] != 16
        or quadrature["lateral_load_rule"]["z_GL"] != 4
        or quadrature["lateral_load_rule"]["point_count"] != 64
        or quadrature["endcap_audit_rule"]["radial_GL"] != 4
        or quadrature["endcap_audit_rule"]["theta_GL"] != 8
        or quadrature["endcap_audit_rule"]["point_count_each"] != 32
    ):
        raise ValueError("boundary quadrature changed from the frozen rule")
    if (
        budget["wall_timeout_seconds"] != WALL_TIMEOUT_SECONDS
        or budget["maximum_source_count"] != MAX_SOURCE_COUNT
        or budget["maximum_registered_dofs"] != MAX_TOTAL_DOFS
        or budget["maximum_volume_quadrature_points_per_assembly"] != MAX_PROVIDER_QUERIES
        or budget["maximum_explicit_numeric_array_bytes"] != MAX_NUMERIC_ARRAY_BYTES
        or budget["maximum_result_json_bytes"] != MAX_RESULT_JSON_BYTES
    ):
        raise ValueError("static solve budgets changed from their frozen limits")
    if solve_case["basis_and_gauge"]["coefficient_count"] != 54:
        raise ValueError("static solve coefficient count changed")
    return solve_case, base_case


def _preflight_counts(
    solve_case: dict[str, object], base_case: dict[str, object]
) -> tuple[dict[str, int], int]:
    source_count = int(base_case["sources"]["source_count"])
    dof_count = int(base_case["dofs"]["coefficient_dof_count"])
    volume_max = max(int(level["volume_point_count"]) for level in solve_case["quadrature"]["volume_levels"])
    boundary_count = (
        int(solve_case["quadrature"]["lateral_load_rule"]["point_count"])
        + 2 * int(solve_case["quadrature"]["endcap_audit_rule"]["point_count_each"])
    )
    offgrid_count = int(solve_case["quadrature"]["offgrid_reference"]["point_count"])
    if source_count > MAX_SOURCE_COUNT or dof_count > MAX_TOTAL_DOFS:
        raise ValueError("budgetRejected: frozen source or coefficient count exceeds the hard cap")
    if volume_max > MAX_PROVIDER_QUERIES or boundary_count > MAX_PROVIDER_QUERIES or offgrid_count > MAX_PROVIDER_QUERIES:
        raise ValueError("budgetRejected: a provider/reference batch exceeds 256 query rows")

    # Conservative numeric-only estimate. It covers all fixed point batches,
    # provider values/gradients, MLS rows, stresses/tractions, two K levels,
    # spectral work, and the 55x55 scaled saddle workspace. Python objects and
    # vendor BLAS workspace are excluded and reported separately.
    all_query_count = volume_max + boundary_count + offgrid_count
    f8 = 8
    matrix_bytes = dof_count * dof_count * f8
    parts = {
        "source_positions": source_count * 3 * f8,
        "volume_boundary_offgrid_positions": all_query_count * 3 * f8,
        "boundary_normals": boundary_count * 3 * f8,
        "volume_and_boundary_weights": (volume_max + boundary_count) * f8,
        "base_mls_weights": all_query_count * source_count * f8,
        "base_mls_gradients": all_query_count * source_count * 3 * f8,
        "enriched_values": all_query_count * dof_count * f8,
        "enriched_gradients": all_query_count * dof_count * 3 * f8,
        "reference_field_stress_traction_and_scalar_outputs": all_query_count * 24 * f8,
        "two_stiffness_levels_and_saddle_spectral_work": 64 * matrix_bytes,
        "load_constraint_solution_and_offgrid_work": (all_query_count * dof_count + 8 * dof_count) * f8,
        "small_B_local_and_quadrature_rule_work": (6 * dof_count + 4096) * f8,
        "saved_result_evidence_arrays": (
            2 * matrix_bytes
            + 3 * dof_count * f8
            + offgrid_count * (3 + 4 + 2 * 4) * f8
        ),
    }
    estimate = int(sum(parts.values()))
    parts["estimated_explicit_numeric_array_bytes"] = estimate
    return parts, estimate


def _load_numerical_dependencies() -> None:
    """Import dependencies only after source pins and scalar budgets pass."""
    global np, AntiplaneShapeSample, BasisKind, DofKey, DofRegistry
    global FieldComponent, ProviderStatus, antiplane_engineering_strain
    global validate_problem_budget, IsotropicMaterial, MLS3DConfig
    global ModeIIIFixedCase, FieldStatus, evaluate_field, ContinuousMLSFrontProvider

    for directory in (MLS_DIR, OPERATOR_DIR, CONTRACT_DIR, FIELD_DIR, PU_OPERATOR_DIR):
        value = str(directory)
        if value not in sys.path:
            sys.path.insert(0, value)
    import numpy as np
    from enriched_dof_contract import (
        AntiplaneShapeSample,
        BasisKind,
        DofKey,
        DofRegistry,
        FieldComponent,
        ProviderStatus,
        antiplane_engineering_strain,
        validate_antiplane_sample,
        validate_problem_budget,
    )
    from efgm_elastic_operator import IsotropicMaterial
    from field_reference import FieldStatus, ModeIIIFixedCase, evaluate_field
    from mls3d_mls import MLS3DConfig
    from mls_modeiii_provider import ContinuousMLSFrontProvider

    # Keep the validator import visible and pinned; provider.evaluate also
    # validates every returned row before it leaves the provider boundary.
    _ = AntiplaneShapeSample, validate_antiplane_sample


def _gauss_rule(base_case: dict[str, object], order: int) -> tuple[Any, Any]:
    rules = base_case["gauss_legendre_rules"]
    rule = rules[f"order_{order}"]
    return np.asarray(rule["nodes"], dtype=np.float64), np.asarray(rule["weights"], dtype=np.float64)


def _map_rule(nodes: Any, weights: Any, lower: float, upper: float) -> tuple[Any, Any]:
    half = 0.5 * (upper - lower)
    middle = 0.5 * (upper + lower)
    return middle + half * nodes, half * weights


def _build_sources(base_case: dict[str, object]) -> tuple[Any, tuple[str, ...]]:
    axis = np.asarray(base_case["sources"]["axis_coordinates_m"], dtype=np.float64)
    rows = list(itertools.product(axis, axis, axis))
    positions = np.ascontiguousarray(np.asarray(rows, dtype=np.float64))
    source_ids = tuple(f"node-{row:03d}" for row in range(len(rows)))
    return positions, source_ids


def _build_registry(source_ids: tuple[str, ...]) -> Any:
    registry = DofRegistry()
    for source_id in source_ids:
        registry.register(DofKey(source_id, FieldComponent.UZ, BasisKind.BASE, 0))
        registry.register(
            DofKey(source_id, FieldComponent.UZ, BasisKind.FRONT_F1_UNSHIFTED, 1)
        )
    registry.freeze()
    return registry


def _build_volume_rule(base_case: dict[str, object], theta_order: int) -> tuple[Any, Any]:
    radius = float(base_case["geometry"]["radius_m"])
    length = float(base_case["geometry"]["length_m"])
    radial_nodes, radial_weights = _gauss_rule(base_case, 4)
    theta_nodes, theta_weights = _gauss_rule(base_case, theta_order)
    axial_nodes, axial_weights = _gauss_rule(base_case, 4)
    radial, wr = _map_rule(radial_nodes, radial_weights, 0.0, radius)
    theta, wt = _map_rule(theta_nodes, theta_weights, -math.pi, math.pi)
    axial, wz = _map_rule(axial_nodes, axial_weights, -0.5 * length, 0.5 * length)
    points: list[tuple[float, float, float]] = []
    volumes: list[float] = []
    for ir, it, iz in itertools.product(range(4), range(theta_order), range(4)):
        r = float(radial[ir])
        angle = float(theta[it])
        points.append((r * math.cos(angle), r * math.sin(angle), float(axial[iz])))
        volumes.append(float(wr[ir] * wt[it] * wz[iz] * r))
    xyz = np.ascontiguousarray(np.asarray(points, dtype=np.float64))
    weights = np.ascontiguousarray(np.asarray(volumes, dtype=np.float64))
    return xyz, weights


def _build_boundary_rule(base_case: dict[str, object]) -> tuple[Any, Any, Any, dict[str, int]]:
    radius = float(base_case["geometry"]["radius_m"])
    length = float(base_case["geometry"]["length_m"])
    theta_nodes, theta_weights = _gauss_rule(base_case, 16)
    axial_nodes, axial_weights = _gauss_rule(base_case, 4)
    radial_nodes, radial_weights = _gauss_rule(base_case, 4)
    cap_theta_nodes, cap_theta_weights = _gauss_rule(base_case, 8)
    theta, wt = _map_rule(theta_nodes, theta_weights, -math.pi, math.pi)
    axial, wz = _map_rule(axial_nodes, axial_weights, -0.5 * length, 0.5 * length)
    radial, wr = _map_rule(radial_nodes, radial_weights, 0.0, radius)
    cap_theta, cwt = _map_rule(cap_theta_nodes, cap_theta_weights, -math.pi, math.pi)

    side_points: list[tuple[float, float, float]] = []
    side_normals: list[tuple[float, float, float]] = []
    side_weights: list[float] = []
    for it, iz in itertools.product(range(16), range(4)):
        angle = float(theta[it])
        cx, sy = math.cos(angle), math.sin(angle)
        side_points.append((radius * cx, radius * sy, float(axial[iz])))
        side_normals.append((cx, sy, 0.0))
        side_weights.append(float(radius * wt[it] * wz[iz]))

    cap_xy: list[tuple[float, float]] = []
    cap_weights: list[float] = []
    for ir, it in itertools.product(range(4), range(8)):
        r = float(radial[ir])
        angle = float(cap_theta[it])
        cap_xy.append((r * math.cos(angle), r * math.sin(angle)))
        cap_weights.append(float(wr[ir] * cwt[it] * r))
    half_length = 0.5 * length
    top_points = [(x, y, half_length) for x, y in cap_xy]
    bottom_points = [(x, y, -half_length) for x, y in cap_xy]
    top_normals = [(0.0, 0.0, 1.0)] * len(top_points)
    bottom_normals = [(0.0, 0.0, -1.0)] * len(bottom_points)

    points = np.ascontiguousarray(
        np.asarray(side_points + top_points + bottom_points, dtype=np.float64)
    )
    normals = np.ascontiguousarray(
        np.asarray(side_normals + top_normals + bottom_normals, dtype=np.float64)
    )
    areas = np.ascontiguousarray(
        np.asarray(side_weights + cap_weights + cap_weights, dtype=np.float64)
    )
    slices = {"side_count": 64, "top_count": 32, "bottom_count": 32}
    return points, normals, areas, slices


def _build_offgrid(base_case: dict[str, object], solve_case: dict[str, object]) -> Any:
    radius = float(base_case["geometry"]["radius_m"])
    length = float(base_case["geometry"]["length_m"])
    spec = solve_case["quadrature"]["offgrid_reference"]
    theta_indices = range(16)
    points: list[tuple[float, float, float]] = []
    for r_fraction, k, z_fraction in itertools.product(
        spec["r_over_R"], theta_indices, spec["z_over_L"]
    ):
        r = float(r_fraction) * radius
        angle = -math.pi + 2.0 * math.pi * (int(k) + 0.37) / 16.0
        z = float(z_fraction) * length
        points.append((r * math.cos(angle), r * math.sin(angle), z))
    return np.ascontiguousarray(np.asarray(points, dtype=np.float64))


def _build_material_cases(base_case: dict[str, object]) -> tuple[Any, Any, Any]:
    material = base_case["material_and_manufactured_field"]
    geometry = base_case["geometry"]
    mls_case = base_case["mls"]
    material_object = IsotropicMaterial(
        young_modulus=float(material["young_modulus_Pa"]),
        poisson_ratio=float(material["poisson_ratio"]),
    )
    mls_config = MLS3DConfig(
        reference_origin=tuple(float(v) for v in mls_case["reference_origin_m"]),
        reference_scale=float(mls_case["reference_scale_m"]),
        support_radius=float(mls_case["support_radius_m"]),
        max_condition_number=float(mls_case["max_condition_number"]),
        rank_relative_tolerance=float(mls_case["rank_relative_tolerance"]),
        switch_tolerance=float(mls_case["switch_tolerance"]),
        duplicate_tolerance=float(mls_case["duplicate_tolerance"]),
    )
    field_case = ModeIIIFixedCase(
        radius_m=float(geometry["radius_m"]),
        length_m=float(geometry["length_m"]),
        K_III_Pa_sqrt_m=float(material["K_III_Pa_sqrt_m"]),
        young_modulus_Pa=float(material["young_modulus_Pa"]),
        poisson_ratio=float(material["poisson_ratio"]),
    )
    return material_object, mls_config, field_case


def _build_provider(
    source_positions: Any,
    source_ids: tuple[str, ...],
    query_positions: Any,
    registry: Any,
    mls_config: Any,
    field_case: Any,
) -> Any:
    return ContinuousMLSFrontProvider(
        source_positions,
        source_ids,
        query_positions,
        registry,
        mls_config,
        field_case,
    )


def _condition_summary(providers: list[Any]) -> dict[str, float | int | None]:
    rows: list[Any] = []
    for provider in providers:
        condition_numbers = provider.mls_result.condition_numbers
        rows.extend(float(value) for value in condition_numbers if np.isfinite(value))
    if not rows:
        return {"finite_row_count": 0, "minimum": None, "maximum": None}
    return {
        "finite_row_count": len(rows),
        "minimum": float(min(rows)),
        "maximum": float(max(rows)),
    }


def _status_counts(statuses: Any) -> dict[str, int]:
    names = sorted({status.value for status in statuses})
    return {name: sum(status.value == name for status in statuses) for name in names}


def _assemble_boundary_load(
    source_positions: Any,
    source_ids: tuple[str, ...],
    query_positions: Any,
    boundary_normals: Any,
    boundary_weights: Any,
    boundary_slices: dict[str, int],
    registry: Any,
    mls_config: Any,
    field_case: Any,
) -> dict[str, object]:
    provider = _build_provider(
        source_positions, source_ids, query_positions, registry, mls_config, field_case
    )
    if provider.unknown_rows:
        return {
            "status": "UNKNOWN_BOUNDARY_PROVIDER_ROW",
            "reason": "boundary MLS or exact-field batch contains UNKNOWN rows",
            "unknown_rows": provider.unknown_rows,
            "provider": provider,
        }
    field = evaluate_field(query_positions, case=field_case, normals=boundary_normals)
    statuses = [field.status_at(index) for index in range(query_positions.shape[0])]
    if any(status is not FieldStatus.OK for status in statuses):
        return {
            "status": "UNKNOWN_BOUNDARY_FIELD_REFERENCE",
            "reason": "a fixed boundary point is not an OK analytic reference query",
            "field_status_counts": _status_counts(statuses),
            "provider": provider,
        }
    traction = field.traction_Pa
    if traction is None:
        return {
            "status": "UNKNOWN_MISSING_ANALYTIC_TRACTION",
            "reason": "field reference returned no sigma*n tensor",
            "provider": provider,
        }

    side_count = int(boundary_slices["side_count"])
    top_count = int(boundary_slices["top_count"])
    bottom_count = int(boundary_slices["bottom_count"])
    if side_count + top_count + bottom_count != int(query_positions.shape[0]):
        return {
            "status": "UNKNOWN_BOUNDARY_DIMENSION_MISMATCH",
            "reason": "boundary point slices do not cover the fixed query batch",
            "provider": provider,
        }
    if np.any(traction[side_count:, 2] != 0.0):
        return {
            "status": "UNKNOWN_NONZERO_SCALAR_CAP_TRACTION",
            "reason": "the exact reference has nonzero t_z on an end cap; scalar UZ scope cannot silently discard it",
            "provider": provider,
        }

    dof_count = registry.dof_count
    force = np.zeros(dof_count, dtype=np.float64)
    side_force = 0.0
    side_abs_force = 0.0
    side_exact_work = 0.0
    exact_uz = field.displacement_m[:, 2]
    for query_id in range(query_positions.shape[0]):
        sample = provider.evaluate(tuple(float(value) for value in query_positions[query_id]))
        if sample.status is not ProviderStatus.OK:
            return {
                "status": "UNKNOWN_BOUNDARY_PROVIDER_STATUS",
                "reason": sample.diagnostic or sample.status.value,
                "blocked_query_index": query_id,
                "provider": provider,
            }
        if query_id < side_count:
            values = np.fromiter(
                (contribution.value for contribution in sample.contributions),
                dtype=np.float64,
                count=dof_count,
            )
            traction_z = float(traction[query_id, 2])
            area = float(boundary_weights[query_id])
            force += values * traction_z * area
            side_force += traction_z * area
            side_abs_force += abs(traction_z) * area
            side_exact_work += float(exact_uz[query_id]) * traction_z * area
    if (
        not np.all(np.isfinite(force))
        or not math.isfinite(side_force)
        or not math.isfinite(side_abs_force)
        or not math.isfinite(side_exact_work)
    ):
        return {
            "status": "UNKNOWN_NONFINITE_BOUNDARY_LOAD",
            "reason": "side traction integration produced a nonfinite scalar load or resultant",
            "provider": provider,
        }
    if provider.consumed_queries != int(query_positions.shape[0]):
        return {
            "status": "UNKNOWN_BOUNDARY_PROVIDER_ORDER",
            "reason": "provider did not consume the exact ordered boundary batch",
            "provider": provider,
        }

    cap_records: dict[str, object] = {}
    cursor = side_count
    for name, count in (("top", top_count), ("bottom", bottom_count)):
        region = slice(cursor, cursor + count)
        region_traction = traction[region]
        region_weights = boundary_weights[region]
        cap_records[name] = {
            "max_abs_full_traction_by_xyz_Pa": [
                float(np.max(np.abs(region_traction[:, axis]))) for axis in range(3)
            ],
            "integrated_full_traction_resultant_N": [
                float(value)
                for value in np.sum(region_traction * region_weights[:, None], axis=0)
            ],
            "max_abs_scalar_tz_Pa": float(np.max(np.abs(region_traction[:, 2]))),
            "uz_virtual_work_J": float(
                np.dot(exact_uz[region] * region_traction[:, 2], region_weights)
            ),
            "full_xy_traction_recorded_but_excluded_from_scalar_uz_rhs": True,
        }
        cursor += count

    return {
        "status": "OK",
        "force": force,
        "side_force_resultant_N": side_force,
        "side_integral_abs_tz_N": side_abs_force,
        "side_force_balance_ratio": abs(side_force) / max(side_abs_force, np.finfo(np.float64).tiny),
        "exact_boundary_work_J": side_exact_work,
        "cap_records": cap_records,
        "field_reference_status_counts": _status_counts(statuses),
        "provider_condition": _condition_summary([provider]),
        "provider": provider,
    }


def _assemble_volume_level(
    level: dict[str, object],
    base_case: dict[str, object],
    source_positions: Any,
    source_ids: tuple[str, ...],
    registry: Any,
    mls_config: Any,
    field_case: Any,
    constitutive: Any,
) -> dict[str, object]:
    theta_order = int(level["theta_GL"])
    points, volume_weights = _build_volume_rule(base_case, theta_order)
    expected_count = int(level["volume_point_count"])
    if points.shape != (expected_count, 3) or volume_weights.shape != (expected_count,):
        return {
            "status": "UNKNOWN_VOLUME_RULE_DIMENSION",
            "name": level["name"],
            "reason": "generated volume point count differs from the fixed case",
        }
    provider = _build_provider(
        source_positions, source_ids, points, registry, mls_config, field_case
    )
    if provider.unknown_rows:
        return {
            "status": "UNKNOWN_VOLUME_PROVIDER_ROW",
            "name": level["name"],
            "reason": "volume MLS or exact-field batch contains UNKNOWN rows",
            "unknown_rows": provider.unknown_rows,
        }

    dof_count = registry.dof_count
    stiffness = np.zeros((dof_count, dof_count), dtype=np.float64)
    physical_volume = math.pi * field_case.radius_m**2 * field_case.length_m
    mean_row = np.zeros(dof_count, dtype=np.float64)
    for query_id in range(expected_count):
        sample = provider.evaluate(tuple(float(value) for value in points[query_id]))
        if sample.status is not ProviderStatus.OK:
            return {
                "status": "UNKNOWN_VOLUME_PROVIDER_STATUS",
                "name": level["name"],
                "reason": sample.diagnostic or sample.status.value,
                "blocked_query_index": query_id,
            }
        b_matrix = np.zeros((6, dof_count), dtype=np.float64)
        value_row = np.empty(dof_count, dtype=np.float64)
        for coefficient_id, contribution in enumerate(sample.contributions):
            b_matrix[:, coefficient_id] = antiplane_engineering_strain(
                contribution.gradient
            )
            value_row[coefficient_id] = contribution.value
        weight = float(volume_weights[query_id])
        with np.errstate(over="ignore", invalid="ignore"):
            local = (b_matrix.T @ constitutive @ b_matrix) * weight
            stiffness += local
            mean_row += value_row * (weight / physical_volume)
        if not np.all(np.isfinite(local)) or not np.all(np.isfinite(stiffness)):
            return {
                "status": "UNKNOWN_NONFINITE_STIFFNESS",
                "name": level["name"],
                "reason": f"scalar stiffness became nonfinite at volume row {query_id}",
            }
        if not np.all(np.isfinite(mean_row)):
            return {
                "status": "UNKNOWN_NONFINITE_MEAN_ROW",
                "name": level["name"],
                "reason": f"normalized physical mean row became nonfinite at volume row {query_id}",
            }
    if provider.consumed_queries != expected_count:
        return {
            "status": "UNKNOWN_VOLUME_PROVIDER_ORDER",
            "name": level["name"],
            "reason": "provider did not consume the exact ordered volume batch",
        }

    exact_energy = float(
        np.dot(provider.field_result.energy_density_J_per_m3, volume_weights)
    )
    if not math.isfinite(exact_energy) or not np.all(np.isfinite(stiffness)):
        return {
            "status": "UNKNOWN_NONFINITE_VOLUME_DIAGNOSTIC",
            "name": level["name"],
            "reason": "volume energy or stiffness result is nonfinite",
        }
    condition = _condition_summary([provider])
    return {
        "status": "OK",
        "name": level["name"],
        "theta_order": theta_order,
        "volume_points": points,
        "volume_weights": volume_weights,
        "stiffness": stiffness,
        "mean_row": mean_row,
        "exact_field_volume_energy_J": exact_energy,
        "volume_measure_sum_m3": float(np.sum(volume_weights)),
        "r_min_sample_m": float(np.min(np.hypot(points[:, 0], points[:, 1]))),
        "mls_condition": condition,
        "mls_status_counts": _status_counts(provider.mls_result.statuses),
        "max_support_count": int(np.max(provider.mls_result.support_counts)),
        "provider": provider,
    }


def _manufactured_coefficients(registry: Any, field_case: Any) -> Any:
    values = np.zeros(registry.dof_count, dtype=np.float64)
    front_value = field_case.amplitude_sqrt_m * math.sqrt(field_case.radius_m)
    for slot, key in enumerate(registry.ordered_keys):
        if key.field_component is not FieldComponent.UZ:
            raise ValueError("static solve registry contains a non-UZ coefficient")
        if key.basis_kind is BasisKind.BASE and key.basis_index == 0:
            values[slot] = 0.0
        elif key.basis_kind is BasisKind.FRONT_F1_UNSHIFTED and key.basis_index == 1:
            values[slot] = front_value
        else:
            raise ValueError(f"unsupported fixed solve basis key {key.basis_kind.value}")
    return values


def _scaled_spectrum(
    stiffness: Any,
    mean_row: Any,
    registry: Any,
    volume_count: int,
    mls_condition: float,
) -> dict[str, object]:
    dof_count = int(stiffness.shape[0])
    if not np.all(np.isfinite(stiffness)):
        return {"status": "UNKNOWN_NONFINITE_STIFFNESS", "reason": "K contains a nonfinite entry"}
    diagonal = np.diag(stiffness)
    if not np.all(np.isfinite(diagonal)) or np.any(diagonal <= 0.0):
        return {"status": "UNKNOWN_NONPOSITIVE_STIFFNESS_DIAGONAL", "reason": "diagonal equilibration requires every K_aa > 0"}
    scale = 1.0 / np.sqrt(diagonal)
    scaled_raw = scale[:, None] * stiffness * scale[None, :]
    symmetric_for_spectrum = 0.5 * (scaled_raw + scaled_raw.T)
    eigenvalues = np.linalg.eigvalsh(symmetric_for_spectrum)
    if (
        not np.all(np.isfinite(eigenvalues))
        or not np.all(np.isfinite(scaled_raw))
        or not np.all(np.isfinite(symmetric_for_spectrum))
    ):
        return {"status": "UNKNOWN_NONFINITE_SPECTRUM", "reason": "scaled stiffness spectrum is nonfinite"}
    spectral_scale = max(float(np.max(np.abs(eigenvalues))), np.finfo(np.float64).tiny)
    epsilon = float(np.finfo(np.float64).eps)
    unit_roundoff = epsilon / 2.0
    operation_count = 12 + volume_count + 2 * dof_count + dof_count * dof_count
    product = operation_count * unit_roundoff
    gamma = product / (1.0 - product)
    zero_tolerance = 4.0 * gamma * spectral_scale
    near_zero = np.flatnonzero(np.abs(eigenvalues) <= zero_tolerance)
    negative_beyond = np.flatnonzero(eigenvalues < -zero_tolerance)

    translation = np.zeros(dof_count, dtype=np.float64)
    expected_translation_norm = 0
    for slot, key in enumerate(registry.ordered_keys):
        if key.basis_kind is BasisKind.BASE:
            translation[slot] = 1.0
            expected_translation_norm += 1
    scaled_translation = translation / scale
    translation_residual = float(
        np.linalg.norm(symmetric_for_spectrum @ scaled_translation)
        / max(
            np.linalg.norm(symmetric_for_spectrum, ord=2) * np.linalg.norm(scaled_translation),
            np.finfo(np.float64).tiny,
        )
    )
    translation_mean = float(np.dot(mean_row, translation))
    translation_mean_error = abs(translation_mean - 1.0)
    translation_residual_tolerance = 8.0 * gamma
    translation_mean_tolerance = 8.0 * gamma

    if not math.isfinite(translation_residual) or not math.isfinite(translation_mean_error):
        status = "UNKNOWN_NONFINITE_TRANSLATION_DIAGNOSTIC"
    elif near_zero.size == 0:
        status = "UNKNOWN_EXPECTED_TRANSLATION_NOT_RESOLVED"
    elif near_zero.size > 1:
        status = "UNKNOWN_EXTRA_NULL_MODE"
    elif negative_beyond.size:
        status = "UNKNOWN_NEGATIVE_SCALED_SPECTRUM"
    elif expected_translation_norm != registry.source_count:
        status = "UNKNOWN_TRANSLATION_BASIS_MISSING"
    elif translation_residual > translation_residual_tolerance:
        status = "UNKNOWN_EXPECTED_TRANSLATION_NOT_NULL"
    elif translation_mean_error > translation_mean_tolerance:
        status = "UNKNOWN_MEAN_GAUGE_TRANSLATION_MISMATCH"
    else:
        status = "OK"

    quotient_condition: float | None = None
    if near_zero.size == 1 and eigenvalues.size > 1 and eigenvalues[1] > zero_tolerance:
        quotient_condition = float(eigenvalues[-1] / eigenvalues[1])
    mean_scaled = scale * mean_row
    mean_scaled_norm = float(np.linalg.norm(mean_scaled))
    gauge_ok = math.isfinite(mean_scaled_norm) and mean_scaled_norm > 0.0
    normalized_mean = mean_scaled / mean_scaled_norm if gauge_ok else None
    saddle = None
    saddle_condition: float | None = None
    if gauge_ok:
        saddle = np.zeros((dof_count + 1, dof_count + 1), dtype=np.float64)
        saddle[:dof_count, :dof_count] = symmetric_for_spectrum
        saddle[:dof_count, dof_count] = normalized_mean
        saddle[dof_count, :dof_count] = normalized_mean
        saddle_condition = float(np.linalg.cond(saddle, p=2))
    if gauge_ok and (
        saddle_condition is None
        or not np.all(np.isfinite(saddle))
        or not math.isfinite(saddle_condition)
    ):
        status = "UNKNOWN_NONFINITE_OR_SINGULAR_SCALED_SADDLE"
    elif not gauge_ok and status == "OK":
        status = "UNKNOWN_SCALED_GAUGE"
    return {
        "status": status,
        "scaled_stiffness": symmetric_for_spectrum,
        "scale_vector": scale,
        "scaled_mean_unit_row": normalized_mean,
        "scaled_mean_norm": mean_scaled_norm,
        "scaled_saddle": saddle,
        "eigenvalues": eigenvalues,
        "zero_tolerance": zero_tolerance,
        "near_zero_count": int(near_zero.size),
        "additional_near_zero_count": max(int(near_zero.size) - 1, 0),
        "negative_eigenvalue_count": int(np.count_nonzero(eigenvalues < 0.0)),
        "negative_beyond_roundoff_count": int(negative_beyond.size),
        "positive_quotient_condition_number": quotient_condition,
        "scaled_saddle_condition_number_2": saddle_condition,
        "scaled_gauge_available": gauge_ok,
        "scaled_gauge_reason": None if gauge_ok else "normalized physical mean row vanishes after coefficient scaling",
        "scaled_stiffness_raw_relative_frobenius_asymmetry": float(
            np.linalg.norm(scaled_raw - scaled_raw.T, ord="fro")
            / max(np.linalg.norm(scaled_raw, ord="fro"), np.finfo(np.float64).tiny)
        ),
        "expected_translation_scaled_residual": translation_residual,
        "expected_translation_scaled_residual_tolerance": translation_residual_tolerance,
        "physical_mean_of_unit_translation": translation_mean,
        "physical_mean_translation_error_from_one": translation_mean_error,
        "physical_mean_translation_roundoff_tolerance": translation_mean_tolerance,
        "diagonal_scaling": "S=diag(K_aa^(-1/2)); Khat=S*K*S",
        "spectrum_matrix": "sym(Khat) for spectrum and symmetric Lagrange solve; raw K remains available for defects and physical residual reporting",
        "roundoff_reference": {
            "binary64_epsilon": epsilon,
            "unit_roundoff": unit_roundoff,
            "conservative_operation_count": operation_count,
            "gamma_n": gamma,
            "eigenvalue_zero_tolerance": zero_tolerance,
            "observed_maximum_mls_condition_number": float(mls_condition),
            "interpretation": "diagnostic label only, not a qualification or convergence gate",
        },
    }


def _solve_scaled_saddle(
    spectrum: dict[str, object], stiffness: Any, mean_row: Any, force: Any
) -> dict[str, object]:
    if spectrum.get("status") != "OK":
        return {"status": str(spectrum.get("status")), "reason": "scaled spectrum did not resolve exactly one expected translation"}
    scale = spectrum["scale_vector"]
    scaled_mean = spectrum["scaled_mean_unit_row"]
    saddle = spectrum["scaled_saddle"]
    dof_count = int(stiffness.shape[0])
    rhs = np.zeros(dof_count + 1, dtype=np.float64)
    rhs[:dof_count] = scale * force
    try:
        solution = np.linalg.solve(saddle, rhs)
    except np.linalg.LinAlgError as exc:
        return {"status": "UNKNOWN_FAILED_DIRECT_SOLVE", "reason": f"direct Lagrange solve failed: {type(exc).__name__}: {exc}"}
    if not np.all(np.isfinite(solution)):
        return {"status": "UNKNOWN_NONFINITE_DIRECT_SOLUTION", "reason": "direct Lagrange solve returned nonfinite unknowns"}
    scaled_coefficients = solution[:dof_count]
    scaled_multiplier = float(solution[dof_count])
    coefficients = scale * scaled_coefficients
    mean_norm = float(spectrum["scaled_mean_norm"])
    multiplier = scaled_multiplier / mean_norm
    residual = stiffness @ coefficients + mean_row * multiplier - force
    scaled_residual = scale * residual
    scaled_kd = scale * (stiffness @ coefficients)
    scaled_clambda = scale * mean_row * multiplier
    scaled_force = scale * force
    residual_denominator = (
        np.linalg.norm(scaled_kd)
        + np.linalg.norm(scaled_clambda)
        + np.linalg.norm(scaled_force)
    )
    weak_residual = float(
        np.linalg.norm(scaled_residual)
        / max(float(residual_denominator), np.finfo(np.float64).tiny)
    )
    return {
        "status": "OK",
        "coefficients_m": coefficients,
        "lagrange_multiplier_N": float(multiplier),
        "scaled_unknowns": solution,
        "weak_residual_relative": weak_residual,
        "physical_mean_m": float(np.dot(mean_row, coefficients)),
        "saddle_backward_residual_relative": float(
            np.linalg.norm(saddle @ solution - rhs)
            / max(np.linalg.norm(rhs), np.finfo(np.float64).tiny)
        ),
    }


def _offgrid_diagnostics(
    source_positions: Any,
    source_ids: tuple[str, ...],
    points: Any,
    registry: Any,
    mls_config: Any,
    field_case: Any,
    solved_levels: list[dict[str, object]],
) -> dict[str, object]:
    provider = _build_provider(
        source_positions, source_ids, points, registry, mls_config, field_case
    )
    if provider.unknown_rows:
        return {
            "status": "UNKNOWN_OFFGRID_PROVIDER_ROW",
            "unknown_rows": provider.unknown_rows,
        }
    if provider.values is None or provider.gradients is None:
        return {"status": "UNKNOWN_MISSING_OFFGRID_ROWS"}
    exact_uz = provider.field_result.displacement_m[:, 2]
    exact_gradient = provider.field_result.displacement_gradient[:, 2, :]
    if (
        not np.all(np.isfinite(points))
        or not np.all(np.isfinite(exact_uz))
        or not np.all(np.isfinite(exact_gradient))
    ):
        return {"status": "UNKNOWN_NONFINITE_OFFGRID_REFERENCE"}
    displacement_scale = max(
        abs(float(field_case.amplitude_sqrt_m)) * math.sqrt(field_case.radius_m),
        np.finfo(np.float64).tiny,
    )
    gradient_scale = max(
        abs(float(field_case.amplitude_sqrt_m)) / math.sqrt(field_case.radius_m),
        np.finfo(np.float64).tiny,
    )
    records: dict[str, object] = {}
    for level in solved_levels:
        coefficients = level["solution"]["coefficients_m"]
        uz = provider.values @ coefficients
        gradient = np.einsum("qdi,d->qi", provider.gradients, coefficients)
        if not np.all(np.isfinite(uz)) or not np.all(np.isfinite(gradient)):
            return {
                "status": "UNKNOWN_NONFINITE_OFFGRID_SOLUTION",
                "blocked_volume_level": str(level["name"]),
            }
        displacement_error = uz - exact_uz
        gradient_error = gradient - exact_gradient
        gradient_error_norm = np.linalg.norm(gradient_error, axis=1)
        if (
            not np.all(np.isfinite(displacement_error))
            or not np.all(np.isfinite(gradient_error))
            or not np.all(np.isfinite(gradient_error_norm))
        ):
            return {
                "status": "UNKNOWN_NONFINITE_OFFGRID_ERROR",
                "blocked_volume_level": str(level["name"]),
            }
        records[str(level["name"])] = {
            "displacement_max_abs_m": float(np.max(np.abs(displacement_error))),
            "displacement_max_relative_to_reference_scale": float(
                np.max(np.abs(displacement_error)) / displacement_scale
            ),
            "displacement_rms_relative_to_reference_scale": float(
                np.sqrt(np.mean(displacement_error**2)) / displacement_scale
            ),
            "gradient_max_abs_dimensionless": float(np.max(gradient_error_norm)),
            "gradient_max_relative_to_reference_scale": float(
                np.max(gradient_error_norm) / gradient_scale
            ),
            "gradient_rms_relative_to_reference_scale": float(
                np.sqrt(np.mean(gradient_error_norm**2)) / gradient_scale
            ),
            "solved_uz_m": uz.tolist(),
            "solved_grad_uz_dimensionless": gradient.tolist(),
            "displacement_reference_scale_A_sqrt_R_m": displacement_scale,
            "gradient_reference_scale_A_over_sqrt_R_dimensionless": gradient_scale,
        }
    return {
        "status": "OK",
        "point_count": int(points.shape[0]),
        "sampling": "fixed nonquadrature 4x16x4 offgrid set from static_solve_case.json",
        "branch_or_front_queries": False,
        "query_coordinates_m": points.tolist(),
        "analytic_uz_m": exact_uz.tolist(),
        "analytic_grad_uz_dimensionless": exact_gradient.tolist(),
        "gradient_field_units": "dimensionless strain-like grad(u_z); provider basis gradients are 1/m",
        "results_by_volume_level": records,
        "provider_condition": _condition_summary([provider]),
    }


def _json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if np is not None and isinstance(value, np.ndarray):
        return _json_safe(value.tolist())
    if np is not None and isinstance(value, np.generic):
        return _json_safe(value.item())
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("nonfinite numeric value cannot be serialized as evidence")
        return float(value)
    return value


def _public_level_record(
    level: dict[str, object],
    d_star: Any,
    force: Any,
    analytic_energy: float,
    roundoff_level: dict[str, object],
    solution: dict[str, object] | None,
) -> dict[str, object]:
    stiffness = level.get("stiffness")
    mean_row = level.get("mean_row")
    volume_weights = level.get("volume_weights")
    if stiffness is None or mean_row is None or volume_weights is None:
        return {key: value for key, value in level.items() if key not in {"provider", "stiffness", "mean_row", "volume_points", "volume_weights"}}
    if (
        not np.all(np.isfinite(stiffness))
        or not np.all(np.isfinite(mean_row))
        or not np.all(np.isfinite(force))
        or not np.all(np.isfinite(d_star))
    ):
        raise ValueError("UNKNOWN_NONFINITE_REPLAYABLE_OPERATOR_EVIDENCE")
    star_defect = stiffness @ d_star - force
    if not np.all(np.isfinite(star_defect)):
        raise ValueError("UNKNOWN_NONFINITE_ANALYTIC_COEFFICIENT_WEAK_DEFECT")
    spectrum = roundoff_level
    scale_vector = spectrum.get("scale_vector")
    if scale_vector is not None:
        scaled_weak_defect_relative: float | None = float(
            np.linalg.norm(scale_vector * star_defect)
            / max(
                np.linalg.norm(scale_vector * (stiffness @ d_star))
                + np.linalg.norm(scale_vector * force),
                np.finfo(np.float64).tiny,
            )
        )
    else:
        scaled_weak_defect_relative = None
    if scaled_weak_defect_relative is not None and not math.isfinite(scaled_weak_defect_relative):
        raise ValueError("UNKNOWN_NONFINITE_SCALED_ANALYTIC_WEAK_DEFECT")
    record: dict[str, object] = {
        "status": level["status"],
        "name": level["name"],
        "theta_order": level["theta_order"],
        "volume_point_count": int(level["volume_points"].shape[0]),
        "volume_measure_sum_m3": level["volume_measure_sum_m3"],
        "analytic_field_volume_energy_quadrature_J": level["exact_field_volume_energy_J"],
        "analytic_field_volume_energy_error_relative": abs(float(level["exact_field_volume_energy_J"]) - analytic_energy) / max(abs(analytic_energy), np.finfo(np.float64).tiny),
        "r_min_positive_sample_m": level["r_min_sample_m"],
        "mls_status_counts": level["mls_status_counts"],
        "mls_condition": level["mls_condition"],
        "max_support_count": level["max_support_count"],
        "relative_frobenius_asymmetry": float(
            np.linalg.norm(stiffness - stiffness.T, ord="fro")
            / max(np.linalg.norm(stiffness, ord="fro"), np.finfo(np.float64).tiny)
        ),
        "stiffness_sha256_float64_c_order": hashlib.sha256(
            np.asarray(stiffness, dtype="<f8", order="C").tobytes(order="C")
        ).hexdigest(),
        "stiffness_matrix_N_per_m": stiffness.tolist(),
        "physical_mean_constraint_row_dimensionless": mean_row.tolist(),
        "analytic_coefficient_vector_m": d_star.tolist(),
        "analytic_coefficient_weak_defect_raw_N": star_defect.tolist(),
        "analytic_coefficient_weak_defect_scaled_relative": scaled_weak_defect_relative,
        "scaled_spectrum": {
            key: value
            for key, value in spectrum.items()
            if key not in {"scaled_stiffness", "scale_vector", "scaled_mean_unit_row", "scaled_saddle"}
        },
    }
    if solution is not None:
        if solution.get("status") != "OK":
            record["solution_attempt"] = {
                "status": solution.get("status"),
                "reason": solution.get("reason"),
                "direct_solver_called": True,
                "partial_solution_returned": False,
            }
            return record
        coefficients = solution["coefficients_m"]
        strain_energy = 0.5 * float(coefficients @ stiffness @ coefficients)
        external_work = float(coefficients @ force)
        energy_error = abs(strain_energy - analytic_energy) / max(abs(analytic_energy), np.finfo(np.float64).tiny)
        work_error = abs(external_work - 2.0 * analytic_energy) / max(
            2.0 * abs(analytic_energy), np.finfo(np.float64).tiny
        )
        normalized_physical_mean = abs(float(solution["physical_mean_m"])) / max(
            float(np.max(np.abs(d_star))) if d_star.size else 0.0,
            np.finfo(np.float64).tiny,
        )
        solution_scalars = (
            strain_energy,
            external_work,
            energy_error,
            work_error,
            float(solution["lagrange_multiplier_N"]),
            float(solution["weak_residual_relative"]),
            float(solution["saddle_backward_residual_relative"]),
            float(solution["physical_mean_m"]),
            normalized_physical_mean,
        )
        if not all(math.isfinite(value) for value in solution_scalars) or not np.all(
            np.isfinite(coefficients)
        ):
            raise ValueError("UNKNOWN_NONFINITE_SOLUTION_ENERGY_WORK_OR_RESIDUAL")
        record["solution"] = {
            "status": solution["status"],
            "coefficients_m": coefficients.tolist(),
            "lagrange_multiplier_N": solution["lagrange_multiplier_N"],
            "weak_residual_relative": solution["weak_residual_relative"],
            "saddle_backward_residual_relative": solution["saddle_backward_residual_relative"],
            "physical_mean_m": solution["physical_mean_m"],
            "physical_mean_normalized_by_A_sqrt_R": normalized_physical_mean,
            "strain_energy_J": strain_energy,
            "energy_error_relative_to_closed_form": energy_error,
            "external_work_J": external_work,
            "external_work_error_relative_to_2Uexact": work_error,
        }
    return record


def _run_static_solve(
    solve_case: dict[str, object],
    base_case: dict[str, object],
    input_hashes: dict[str, str],
    budget_parts: dict[str, int],
    started_at: str,
    started: float,
) -> dict[str, object]:
    sources, source_ids = _build_sources(base_case)
    validate_problem_budget(len(source_ids), 54, 256)
    registry = _build_registry(source_ids)
    if registry.source_count != len(source_ids) or registry.dof_count != 54:
        return {
            "schema": "modeiii_pu_antiplane_static_solve_result.v3",
            "case_id": CASE_ID,
            "status": "UNKNOWN_REGISTRY_DIMENSION",
            "classification": "DIAGNOSTIC_ONLY_UNQUALIFIED",
            "reason": "frozen canonical registry has unexpected dimensions",
            "source_sha256": input_hashes,
            "finished_at_utc": _timestamp(),
        }
    material, mls_config, field_case = _build_material_cases(base_case)
    d_star = _manufactured_coefficients(registry, field_case)
    analytic_energy = float(field_case.exact_total_strain_energy_J)
    physical_volume = math.pi * field_case.radius_m**2 * field_case.length_m
    constitutive = material.constitutive_matrix

    boundary_points, boundary_normals, boundary_weights, boundary_slices = _build_boundary_rule(base_case)
    boundary_result = _assemble_boundary_load(
        sources,
        source_ids,
        boundary_points,
        boundary_normals,
        boundary_weights,
        boundary_slices,
        registry,
        mls_config,
        field_case,
    )
    if boundary_result["status"] != "OK":
        boundary_result.pop("provider", None)
        return {
            "schema": "modeiii_pu_antiplane_static_solve_result.v3",
            "case_id": CASE_ID,
            "status": str(boundary_result["status"]),
            "classification": "DIAGNOSTIC_ONLY_UNQUALIFIED",
            "started_at_utc": started_at,
            "reason": boundary_result.get("reason"),
            "boundary_diagnostic": boundary_result,
            "source_sha256": input_hashes,
            "array_budget": {"parts_bytes": budget_parts, "cap_bytes": MAX_NUMERIC_ARRAY_BYTES},
            "no_partial_stiffness_returned": True,
            "finished_at_utc": _timestamp(),
        }
    force = boundary_result["force"]
    boundary_condition = boundary_result["provider_condition"]
    boundary_public = {
        key: value
        for key, value in boundary_result.items()
        if key not in {"force", "provider"}
    }
    boundary_result.pop("provider", None)

    assembled: list[dict[str, object]] = []
    for level in solve_case["quadrature"]["volume_levels"]:
        if time.monotonic() - started > WALL_TIMEOUT_SECONDS:
            return {
                "schema": "modeiii_pu_antiplane_static_solve_result.v3",
                "case_id": CASE_ID,
                "status": "UNKNOWN_WALL_BUDGET",
                "classification": "DIAGNOSTIC_ONLY_UNQUALIFIED",
                "started_at_utc": started_at,
                "reason": "30-second work budget was reached during volume assembly",
                "completed_volume_levels": [item["name"] for item in assembled],
                "source_sha256": input_hashes,
                "array_budget": {"parts_bytes": budget_parts, "cap_bytes": MAX_NUMERIC_ARRAY_BYTES},
                "no_solver_attempted": True,
                "finished_at_utc": _timestamp(),
            }
        level_result = _assemble_volume_level(
            level,
            base_case,
            sources,
            source_ids,
            registry,
            mls_config,
            field_case,
            constitutive,
        )
        if level_result["status"] != "OK":
            return {
                "schema": "modeiii_pu_antiplane_static_solve_result.v3",
                "case_id": CASE_ID,
                "status": str(level_result["status"]),
                "classification": "DIAGNOSTIC_ONLY_UNQUALIFIED",
                "started_at_utc": started_at,
                "reason": level_result.get("reason"),
                "blocked_volume_level": level_result.get("name"),
                "completed_volume_levels": [item["name"] for item in assembled],
                "unknown_rows": level_result.get("unknown_rows", []),
                "source_sha256": input_hashes,
                "array_budget": {"parts_bytes": budget_parts, "cap_bytes": MAX_NUMERIC_ARRAY_BYTES},
                "no_partial_stiffness_returned": True,
                "finished_at_utc": _timestamp(),
            }
        assembled.append(level_result)

    spectra: list[dict[str, object]] = []
    for level in assembled:
        spectrum = _scaled_spectrum(
            level["stiffness"],
            level["mean_row"],
            registry,
            int(level["volume_points"].shape[0]),
            float(level["mls_condition"]["maximum"] or 1.0),
        )
        spectra.append(spectrum)
    blocked_spectra = [index for index, spectrum in enumerate(spectra) if spectrum.get("status") != "OK"]
    if blocked_spectra:
        level_records = [
            _public_level_record(level, d_star, force, analytic_energy, spectra[index], None)
            for index, level in enumerate(assembled)
        ]
        return {
            "schema": "modeiii_pu_antiplane_static_solve_result.v3",
            "case_id": CASE_ID,
            "status": "UNKNOWN_MATRIX_SPECTRUM",
            "classification": "DIAGNOSTIC_ONLY_UNQUALIFIED",
            "started_at_utc": started_at,
            "reason": "one or more scaled stiffness spectra do not resolve exactly one scalar translation",
            "blocked_level_names": [assembled[index]["name"] for index in blocked_spectra],
            "boundary": boundary_public,
            "applied_force_coefficients_N": force.tolist(),
            "ordered_canonical_coefficient_ids": list(registry.canonical_ids),
            "levels": level_records,
            "expected_scalar_nullity_before_gauge": 1,
            "additional_null_modes_suppressed": False,
            "no_solver_attempted": True,
            "source_sha256": input_hashes,
            "array_budget": {
                "parts_bytes": budget_parts,
                "estimated_explicit_numeric_array_bytes": budget_parts["estimated_explicit_numeric_array_bytes"],
                "cap_bytes": MAX_NUMERIC_ARRAY_BYTES,
                "within_estimate_cap": True,
                "python_object_overhead_and_blas_workspace_excluded": True,
            },
            "finished_at_utc": _timestamp(),
        }

    solved_levels: list[dict[str, object]] = []
    failed_levels: list[str] = []
    for level, spectrum in zip(assembled, spectra):
        solution = _solve_scaled_saddle(
            spectrum, level["stiffness"], level["mean_row"], force
        )
        if solution["status"] != "OK":
            failed_levels.append(str(level["name"]))
        solved_levels.append({**level, "spectrum": spectrum, "solution": solution})
    if failed_levels:
        level_records = [
            _public_level_record(level, d_star, force, analytic_energy, level["spectrum"], level["solution"])
            for level in solved_levels
        ]
        return {
            "schema": "modeiii_pu_antiplane_static_solve_result.v3",
            "case_id": CASE_ID,
            "status": "UNKNOWN_DIRECT_SOLVE_FAILURE",
            "classification": "DIAGNOSTIC_ONLY_UNQUALIFIED",
            "started_at_utc": started_at,
            "reason": "at least one fixed quadrature Lagrange system failed direct factorization or returned nonfinite values",
            "failed_volume_levels": failed_levels,
            "boundary": boundary_public,
            "applied_force_coefficients_N": force.tolist(),
            "ordered_canonical_coefficient_ids": list(registry.canonical_ids),
            "levels": level_records,
            "no_pseudoinverse_or_regularization": True,
            "source_sha256": input_hashes,
            "array_budget": {"parts_bytes": budget_parts, "cap_bytes": MAX_NUMERIC_ARRAY_BYTES},
            "finished_at_utc": _timestamp(),
        }

    offgrid_points = _build_offgrid(base_case, solve_case)
    offgrid = _offgrid_diagnostics(
        sources,
        source_ids,
        offgrid_points,
        registry,
        mls_config,
        field_case,
        solved_levels,
    )
    if offgrid["status"] != "OK":
        return {
            "schema": "modeiii_pu_antiplane_static_solve_result.v3",
            "case_id": CASE_ID,
            "status": str(offgrid["status"]),
            "classification": "DIAGNOSTIC_ONLY_UNQUALIFIED",
            "started_at_utc": started_at,
            "reason": "independent off-grid field reference or MLS batch is UNKNOWN",
            "boundary": boundary_public,
            "applied_force_coefficients_N": force.tolist(),
            "ordered_canonical_coefficient_ids": list(registry.canonical_ids),
            "levels": [
                _public_level_record(level, d_star, force, analytic_energy, level["spectrum"], level["solution"])
                for level in solved_levels
            ],
            "offgrid": offgrid,
            "source_sha256": input_hashes,
            "array_budget": {"parts_bytes": budget_parts, "cap_bytes": MAX_NUMERIC_ARRAY_BYTES},
            "finished_at_utc": _timestamp(),
        }

    level_records = [
        _public_level_record(level, d_star, force, analytic_energy, level["spectrum"], level["solution"])
        for level in solved_levels
    ]
    return {
        "schema": "modeiii_pu_antiplane_static_solve_result.v3",
        "case_id": CASE_ID,
        "status": "SOLVED_DIAGNOSTIC_ONLY_UNQUALIFIED",
        "classification": "ENGINEERING_MANUFACTURED_ANTIPLANE_DIAGNOSTIC_ONLY",
        "started_at_utc": started_at,
        "fixed_case": {"solve_case": solve_case, "base_operator_case": base_case},
        "dimensions": {
            "source_count": len(source_ids),
            "scalar_dof_count": registry.dof_count,
            "lateral_load_points": 64,
            "endcap_audit_points_each": 32,
            "offgrid_points": int(offgrid_points.shape[0]),
        },
        "registry": {
            "state": registry.state.value,
            "source_count": registry.source_count,
            "dof_count": registry.dof_count,
            "ordered_canonical_coefficient_ids": list(registry.canonical_ids),
            "matrix_axes_and_coefficient_vectors_follow_this_order": True,
            "basis_counts": {
                "BASE": sum(key.basis_kind is BasisKind.BASE for key in registry.ordered_keys),
                "FRONT_F1_UNSHIFTED": sum(key.basis_kind is BasisKind.FRONT_F1_UNSHIFTED for key in registry.ordered_keys),
            },
            "front_support_set_T_F": list(source_ids),
            "heaviside_support_set_J_H": [],
            "interpretation": "all-site front enrichment is an algebraic bridge only, not local convergence evidence",
        },
        "boundary": boundary_public,
        "applied_force_coefficients_N": force.tolist(),
        "volume_levels": level_records,
        "offgrid": offgrid,
        "expected_closed_form_strain_energy_J": analytic_energy,
        "physical_volume_m3": physical_volume,
        "solver": {
            "direct_lagrange_system": True,
            "coefficient_pinning": False,
            "pseudoinverse_or_regularization": False,
            "expected_removed_physical_modes": ["uniform uz translation"],
            "additional_null_modes_hidden": False,
        },
        "array_budget": {
            "parts_bytes": budget_parts,
            "estimated_explicit_numeric_array_bytes": budget_parts["estimated_explicit_numeric_array_bytes"],
            "cap_bytes": MAX_NUMERIC_ARRAY_BYTES,
            "maximum_result_json_bytes": MAX_RESULT_JSON_BYTES,
            "within_estimate_cap": True,
            "python_object_overhead_and_blas_workspace_excluded": True,
        },
        "source_sha256": input_hashes,
        "elapsed_seconds": time.monotonic() - started,
        "wall_timeout_seconds": WALL_TIMEOUT_SECONDS,
        "within_wall_timeout": time.monotonic() - started <= WALL_TIMEOUT_SECONDS,
        "error_values_are_diagnostics_only": True,
        "limitations": solve_case["field_scope"]["qualification"],
        "finished_at_utc": _timestamp(),
    }


def _write_result(path: Path, result: dict[str, object]) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(_json_safe(result), indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    encoded = payload.encode("utf-8")
    status = str(result.get("status", "UNKNOWN_RESULT_STATUS"))
    if len(encoded) > MAX_RESULT_JSON_BYTES:
        overflow = {
            "schema": "modeiii_pu_antiplane_static_solve_result.v3",
            "case_id": CASE_ID,
            "status": "UNKNOWN_RESULT_JSON_BUDGET",
            "classification": "DIAGNOSTIC_ONLY_UNQUALIFIED",
            "candidate_result_json_bytes": len(encoded),
            "maximum_result_json_bytes": MAX_RESULT_JSON_BYTES,
            "large_numeric_evidence_promoted": False,
            "solver_source_sha256": result.get("solver_source_sha256"),
            "finished_at_utc": _timestamp(),
        }
        payload = json.dumps(overflow, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
        encoded = payload.encode("utf-8")
        status = str(overflow["status"])
    path.write_bytes(encoded)
    return status


def _worker(output_path: Path) -> int:
    started = time.monotonic()
    started_at = _timestamp()
    result: dict[str, object]
    input_hashes: dict[str, str] = {}
    solve_case: dict[str, object] | None = None
    try:
        input_hashes = _verify_source_pins()
        solve_case, base_case = _read_and_validate_cases()
        budget_parts, estimate = _preflight_counts(solve_case, base_case)
        if estimate > MAX_NUMERIC_ARRAY_BYTES:
            result = {
                "schema": "modeiii_pu_antiplane_static_solve_result.v3",
                "case_id": CASE_ID,
                "status": "UNKNOWN_ARRAY_BUDGET",
                "classification": "DIAGNOSTIC_ONLY_UNQUALIFIED",
                "started_at_utc": started_at,
                "reason": "preflight numeric-array estimate exceeds 32 MiB",
                "array_budget": {"parts_bytes": budget_parts, "cap_bytes": MAX_NUMERIC_ARRAY_BYTES},
                "source_sha256": input_hashes,
                "numeric_imports_or_allocations_started": False,
                "finished_at_utc": _timestamp(),
            }
        else:
            _load_numerical_dependencies()
            result = _run_static_solve(
                solve_case,
                base_case,
                input_hashes,
                budget_parts,
                started_at,
                started,
            )
            if time.monotonic() - started > WALL_TIMEOUT_SECONDS:
                result["status"] = "UNKNOWN_WALL_BUDGET"
                result["reason"] = "30-second worker budget exceeded before output serialization"
                result["solver_result_is_not_promoted"] = True
    except Exception as exc:
        result = {
            "schema": "modeiii_pu_antiplane_static_solve_result.v3",
            "case_id": CASE_ID,
            "status": "UNKNOWN_SOURCE_OR_NUMERICS_ERROR",
            "classification": "DIAGNOSTIC_ONLY_UNQUALIFIED",
            "started_at_utc": started_at,
            "error_type": type(exc).__name__,
            "error": str(exc),
            "traceback": traceback.format_exc(limit=8),
            "no_retry": True,
            "source_sha256": input_hashes,
            "finished_at_utc": _timestamp(),
        }

    result["solver_source_sha256"] = _sha256(Path(__file__).resolve())
    if solve_case is not None:
        result["quadrature_comparison"] = solve_case["quadrature_comparison"]
    # The required UNKNOWN record reaches disk before the failure code is
    # returned. There are no plotting or repair paths.
    try:
        written_status = _write_result(output_path, result)
    except (TypeError, ValueError, OverflowError) as exc:
        serialization_failure = {
            "schema": "modeiii_pu_antiplane_static_solve_result.v3",
            "case_id": CASE_ID,
            "status": "UNKNOWN_RESULT_SERIALIZATION",
            "classification": "DIAGNOSTIC_ONLY_UNQUALIFIED",
            "reason": "result contained unsupported or nonfinite numeric evidence; original values were not converted to null",
            "error_type": type(exc).__name__,
            "error": str(exc),
            "solver_source_sha256": result.get("solver_source_sha256"),
            "finished_at_utc": _timestamp(),
        }
        written_status = _write_result(output_path, serialization_failure)
    if written_status != "SOLVED_DIAGNOSTIC_ONLY_UNQUALIFIED":
        raise StaticSolveUnknown(written_status)
    return 0


def _stage_metadata(staging: Path) -> dict[str, object]:
    metadata: dict[str, object] = {"worker_staging_path": str(staging.resolve())}
    try:
        metadata["worker_staging_json_exists"] = staging.exists()
        metadata["worker_staging_is_file"] = staging.is_file()
        if staging.is_file():
            size = staging.stat().st_size
            metadata["worker_staging_bytes"] = size
            if size <= MAX_RESULT_JSON_BYTES:
                metadata["worker_staging_sha256"] = _sha256(staging)
            else:
                metadata["worker_staging_sha256_not_computed"] = "staging exceeds the 4 MiB result cap"
    except OSError as exc:
        metadata["worker_staging_read_error"] = f"{type(exc).__name__}: {exc}"
    return metadata


def _finite_number(value: object) -> bool:
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        return False
    try:
        return math.isfinite(float(value))
    except (OverflowError, TypeError, ValueError):
        return False


def _finite_vector(value: object, size: int) -> bool:
    return (
        isinstance(value, list)
        and len(value) == size
        and all(_finite_number(item) for item in value)
    )


def _finite_matrix(value: object, rows: int, columns: int) -> bool:
    return (
        isinstance(value, list)
        and len(value) == rows
        and all(_finite_vector(row, columns) for row in value)
    )


def _reject_json_constant(value: str) -> None:
    raise ValueError(f"nonstandard JSON numeric constant {value!r}")


def _validate_worker_success_payload(result: object) -> str | None:
    """Validate only the frozen success envelope and replayable evidence shape."""
    if not isinstance(result, dict):
        return "root must be a JSON object"
    expected_top_level = {
        "schema", "case_id", "status", "classification", "started_at_utc",
        "finished_at_utc", "solver_source_sha256", "source_sha256", "dimensions",
        "registry", "applied_force_coefficients_N", "volume_levels", "offgrid",
        "solver", "array_budget", "fixed_case", "boundary",
        "expected_closed_form_strain_energy_J", "physical_volume_m3",
        "quadrature_comparison", "elapsed_seconds", "wall_timeout_seconds",
        "within_wall_timeout", "error_values_are_diagnostics_only", "limitations",
    }
    missing = sorted(expected_top_level.difference(result))
    if missing:
        return f"missing required top-level fields: {missing}"
    if result["schema"] != "modeiii_pu_antiplane_static_solve_result.v3":
        return "unexpected result schema"
    if result["case_id"] != CASE_ID:
        return "unexpected case ID"
    if result["status"] != "SOLVED_DIAGNOSTIC_ONLY_UNQUALIFIED":
        return "worker result is not a solved diagnostic"
    if result["classification"] != "ENGINEERING_MANUFACTURED_ANTIPLANE_DIAGNOSTIC_ONLY":
        return "unexpected result classification"
    if not isinstance(result["started_at_utc"], str) or not isinstance(result["finished_at_utc"], str):
        return "missing timestamp strings"
    if result["solver_source_sha256"] != _sha256(Path(__file__).resolve()):
        return "worker source hash does not match the frozen v3 source"
    source_hashes = result["source_sha256"]
    if not isinstance(source_hashes, dict):
        return "source_sha256 must be an object"
    if source_hashes != EXPECTED_SHA256:
        return "worker dependency hashes differ from the frozen source pins"
    if not all(_finite_number(result[name]) for name in ("expected_closed_form_strain_energy_J", "physical_volume_m3", "elapsed_seconds", "wall_timeout_seconds")):
        return "nonfinite or invalid required scalar evidence"
    if result["within_wall_timeout"] is not True or result["error_values_are_diagnostics_only"] is not True:
        return "required success flags are false"

    dimensions = result["dimensions"]
    if not isinstance(dimensions, dict) or any(dimensions.get(key) != expected for key, expected in {
        "source_count": 27, "scalar_dof_count": 54, "lateral_load_points": 64,
        "endcap_audit_points_each": 32, "offgrid_points": 256,
    }.items()):
        return "dimensions do not match the frozen case"
    registry = result["registry"]
    if not isinstance(registry, dict):
        return "registry must be an object"
    canonical_ids = registry.get("ordered_canonical_coefficient_ids")
    if (
        not isinstance(canonical_ids, list)
        or len(canonical_ids) != 54
        or not all(isinstance(item, str) and item for item in canonical_ids)
        or len(set(canonical_ids)) != 54
        or registry.get("dof_count") != 54
        or result.get("ordered_canonical_coefficient_ids", canonical_ids) != canonical_ids
    ):
        return "registry canonical coefficient IDs are missing or malformed"
    if not _finite_vector(result["applied_force_coefficients_N"], 54):
        return "applied force evidence must contain 54 finite coefficients"

    levels = result["volume_levels"]
    if not isinstance(levels, list) or len(levels) != 2:
        return "exactly two volume levels are required"
    required_level_names = {"volume_gl4x8x4", "volume_gl4x16x4"}
    seen_names: set[str] = set()
    for level in levels:
        if not isinstance(level, dict):
            return "each volume level must be an object"
        name = level.get("name")
        if not isinstance(name, str) or name not in required_level_names or name in seen_names:
            return "volume level names do not match the frozen two-level comparison"
        seen_names.add(name)
        if level.get("status") != "OK":
            return f"volume level {name} is not OK"
        if not _finite_matrix(level.get("stiffness_matrix_N_per_m"), 54, 54):
            return f"volume level {name} has invalid stiffness matrix evidence"
        if not _finite_vector(level.get("physical_mean_constraint_row_dimensionless"), 54):
            return f"volume level {name} has invalid physical mean row evidence"
        if not _finite_vector(level.get("analytic_coefficient_vector_m"), 54):
            return f"volume level {name} has invalid analytic coefficient evidence"
        if not _finite_vector(level.get("analytic_coefficient_weak_defect_raw_N"), 54):
            return f"volume level {name} has invalid analytic weak-defect evidence"
        expected_points = 128 if name == "volume_gl4x8x4" else 256
        expected_theta_order = 8 if name == "volume_gl4x8x4" else 16
        if level.get("volume_point_count") != expected_points or level.get("theta_order") != expected_theta_order:
            return f"volume level {name} metadata does not match the frozen quadrature"
        solution = level.get("solution")
        if not isinstance(solution, dict) or solution.get("status") != "OK":
            return f"volume level {name} lacks a successful solution record"
        if not _finite_vector(solution.get("coefficients_m"), 54):
            return f"volume level {name} has invalid solution coefficients"
        for key in (
            "lagrange_multiplier_N", "weak_residual_relative", "saddle_backward_residual_relative",
            "physical_mean_m", "physical_mean_normalized_by_A_sqrt_R", "strain_energy_J",
            "energy_error_relative_to_closed_form", "external_work_J", "external_work_error_relative_to_2Uexact",
        ):
            if not _finite_number(solution.get(key)):
                return f"volume level {name} has invalid scalar solution evidence: {key}"

    offgrid = result["offgrid"]
    if not isinstance(offgrid, dict) or offgrid.get("status") != "OK" or offgrid.get("point_count") != 256:
        return "off-grid evidence is missing or not OK"
    if (
        not _finite_matrix(offgrid.get("query_coordinates_m"), 256, 3)
        or not _finite_vector(offgrid.get("analytic_uz_m"), 256)
        or not _finite_matrix(offgrid.get("analytic_grad_uz_dimensionless"), 256, 3)
    ):
        return "off-grid reference arrays have invalid dimensions or values"
    solved = offgrid.get("results_by_volume_level")
    if not isinstance(solved, dict) or set(solved) != required_level_names:
        return "off-grid solved fields do not match the frozen volume levels"
    for name in required_level_names:
        record = solved[name]
        if not isinstance(record, dict):
            return f"off-grid result for {name} must be an object"
        if (
            not _finite_vector(record.get("solved_uz_m"), 256)
            or not _finite_matrix(record.get("solved_grad_uz_dimensionless"), 256, 3)
        ):
            return f"off-grid solved fields for {name} have invalid dimensions or values"
    return None


def _parent_run(output_path: Path) -> int:
    started = time.monotonic()
    deadline = started + WALL_TIMEOUT_SECONDS
    if output_path.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    staging = output_path.with_name(f".{output_path.name}.{os.getpid()}.worker.json")
    if staging.exists():
        raise FileExistsError(f"refusing to overwrite existing staging result: {staging}")
    remaining = max(0.001, deadline - time.monotonic())
    command = [sys.executable, "-B", str(Path(__file__).resolve()), "--_worker", "--output", str(staging)]
    try:
        completed = subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=remaining,
            check=False,
        )
    except subprocess.TimeoutExpired:
        timeout_result = {
            "schema": "modeiii_pu_antiplane_static_solve_result.v3",
            "case_id": CASE_ID,
            "status": "UNKNOWN_WALL_TIMEOUT",
            "classification": "DIAGNOSTIC_ONLY_UNQUALIFIED",
            "reason": "hard 30-second parent deadline expired; no retry was attempted",
            **_stage_metadata(staging),
            "staging_preserved": True,
            "no_solver_result_promoted": True,
            "finished_at_utc": _timestamp(),
        }
        _write_result(output_path, timeout_result)
        return 2

    if int(completed.returncode) != 0:
        failure_result = {
            "schema": "modeiii_pu_antiplane_static_solve_result.v3",
            "case_id": CASE_ID,
            "status": "UNKNOWN_WORKER_NONZERO_EXIT",
            "classification": "DIAGNOSTIC_ONLY_UNQUALIFIED",
            "reason": "worker exit code was nonzero; any staging file is preserved without promotion",
            "worker_return_code": int(completed.returncode),
            **_stage_metadata(staging),
            "staging_preserved": True,
            "no_solver_result_promoted": True,
            "worker_stdout": completed.stdout[-4000:],
            "worker_stderr": completed.stderr[-4000:],
            "no_retry": True,
            "finished_at_utc": _timestamp(),
        }
        _write_result(output_path, failure_result)
        return 2

    if not staging.is_file():
        failure_result = {
            "schema": "modeiii_pu_antiplane_static_solve_result.v3",
            "case_id": CASE_ID,
            "status": "UNKNOWN_WORKER_DID_NOT_WRITE_JSON",
            "classification": "DIAGNOSTIC_ONLY_UNQUALIFIED",
            "reason": "worker exited successfully but did not produce a staging JSON file",
            "worker_return_code": int(completed.returncode),
            **_stage_metadata(staging),
            "staging_preserved": True,
            "no_solver_result_promoted": True,
            "worker_stdout": completed.stdout[-4000:],
            "worker_stderr": completed.stderr[-4000:],
            "no_retry": True,
            "finished_at_utc": _timestamp(),
        }
        _write_result(output_path, failure_result)
        return 2

    try:
        with staging.open("rb") as stream:
            staged_bytes = stream.read(MAX_RESULT_JSON_BYTES + 1)
        if len(staged_bytes) > MAX_RESULT_JSON_BYTES:
            raise ValueError("staging JSON exceeds the frozen 4 MiB result cap")
        staged_sha256 = hashlib.sha256(staged_bytes).hexdigest()
        staged_result = json.loads(
            staged_bytes.decode("utf-8"), parse_constant=_reject_json_constant
        )
    except (OSError, UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        failure_result = {
            "schema": "modeiii_pu_antiplane_static_solve_result.v3",
            "case_id": CASE_ID,
            "status": "UNKNOWN_WORKER_INVALID_JSON",
            "classification": "DIAGNOSTIC_ONLY_UNQUALIFIED",
            "reason": "worker staging bytes could not be parsed as strict UTF-8 JSON",
            "parse_error_type": type(exc).__name__,
            "parse_error": str(exc),
            **_stage_metadata(staging),
            "staging_preserved": True,
            "no_solver_result_promoted": True,
            "worker_stdout": completed.stdout[-4000:],
            "worker_stderr": completed.stderr[-4000:],
            "no_retry": True,
            "finished_at_utc": _timestamp(),
        }
        _write_result(output_path, failure_result)
        return 2

    try:
        contract_error = _validate_worker_success_payload(staged_result)
    except Exception as exc:
        contract_error = f"success-envelope validation failed closed: {type(exc).__name__}: {exc}"
    if contract_error is not None:
        failure_result = {
            "schema": "modeiii_pu_antiplane_static_solve_result.v3",
            "case_id": CASE_ID,
            "status": "UNKNOWN_WORKER_RESULT_CONTRACT",
            "classification": "DIAGNOSTIC_ONLY_UNQUALIFIED",
            "reason": contract_error,
            "worker_status": staged_result.get("status") if isinstance(staged_result, dict) else None,
            "worker_return_code": int(completed.returncode),
            **_stage_metadata(staging),
            "staging_preserved": True,
            "no_solver_result_promoted": True,
            "no_retry": True,
            "finished_at_utc": _timestamp(),
        }
        _write_result(output_path, failure_result)
        return 2

    try:
        if output_path.exists():
            raise FileExistsError(f"refusing to replace concurrently created output: {output_path}")
        if _sha256(staging) != staged_sha256:
            raise OSError("staging file changed after validation")
        os.replace(staging, output_path)
    except OSError as exc:
        failure_result = {
            "schema": "modeiii_pu_antiplane_static_solve_result.v3",
            "case_id": CASE_ID,
            "status": "UNKNOWN_WORKER_PROMOTION_FAILED",
            "classification": "DIAGNOSTIC_ONLY_UNQUALIFIED",
            "reason": f"validated staging JSON was not promoted: {type(exc).__name__}: {exc}",
            **_stage_metadata(staging),
            "staging_preserved": staging.exists(),
            "no_solver_result_promoted": True,
            "no_retry": True,
            "finished_at_utc": _timestamp(),
        }
        if not output_path.exists():
            _write_result(output_path, failure_result)
        return 2
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=HERE / "results" / "static_solve_diagnostic.json",
        help="new JSON output path; existing files are not overwritten",
    )
    parser.add_argument("--_worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args._worker:
        try:
            return _worker(args.output)
        except StaticSolveUnknown as exc:
            print(f"UNKNOWN: {exc}", file=sys.stderr)
            return 2
    return _parent_run(args.output)


if __name__ == "__main__":
    raise SystemExit(main())
