"""Bounded source-only assembly of the pinned local-front scalar Poisson K.

This is an antiplane scalar bridge over the fixed finite-crack box. It is not
the full 3-D isotropic elasticity operator and does not solve a boundary-value
problem. Numerical imports and all query/operator work are confined to the
hard-timed worker process.
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
import struct
import sys
import time
import traceback


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
CASE_PATH = HERE / "fixed_case.json"
PROVIDER_DIR = REPO / "pgfrag" / "experiments" / "local_front_provider_20261006"
SUPPORT_DIR = REPO / "pgfrag" / "experiments" / "local_front_support_contract_20261006"
REGISTRY_DIR = REPO / "pgfrag" / "experiments" / "enriched_dof_contract_20261005"
MLS_DIR = REPO / "pgfrag" / "experiments" / "mls3d_kernel_20261005"
MATERIAL_DIR = REPO / "pgfrag" / "experiments" / "modeiii_pu_operator_20261005"

CASE_ID = "local_front_scalar_antiplane_poisson_operator_20261006"
MAX_SOURCE_COUNT = 125
MAX_TOTAL_DOFS = 375
MAX_QUADRATURE_POINTS = 256
MAX_EXPLICIT_NUMERIC_ARRAY_BYTES = 8 * 1024 * 1024
WALL_TIMEOUT_SECONDS = 60.0

# These are read-only roots already pinned by the local-front provider; the
# fixed case and provider identities are pinned here as well.
EXPECTED_DEPENDENCIES = {
    "local_front_operator_fixed_case.json": "bee0638bdbafdcd9dd1af02152faf3d87384503b6a029a70831f71b5369c162b",
    "local_front_provider.py": "c586ec99dee8dc61f7a74904e1c0ce196d8ca562cbc2bb6f61411a8262c96421",
    "support_fixed_case.json": "0bdbb20032e105851d2909ff54850ebc54d3f2041c25ab5e42e9102244d90ca5",
    "support_topology.py": "ea20e51d2b774c7e6b8b48403359b56bab1ba6f5405dba41f34813b6aec4578e",
    "role_manifest.json": "e5e318f80b5dfedbc8ba122f82df3cca1585e604272399d21bc35907d1d28f01",
    "enriched_dof_contract.py": "e3ddb3b067465bb5b269576a1b8a5100e2bb92c9c924bc09c949702c19efa9d8",
    "mls3d_mls.py": "404256d9b7546361a5704c4c3cec7c022069a9cf83a137c58a16e7a64f79ec63",
    "mls3d_shape.py": "66e680616a988182875d17aeb84fadb541867aee13fc15e0ae098f7c98839ef5",
    "modeiii_material_case.json": "dce0a1f3009246a52198bc0ed833ed67998144549a653e7bce9bcc17babd71b9",
}


class DiagnosticBudgetError(RuntimeError):
    """Raised when the frozen source, DOF, query, byte, or time cap is exceeded."""


class DiagnosticNumericsError(FloatingPointError):
    """Raised when an operator or diagnostic quantity is nonfinite."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _verify_dependencies() -> dict[str, str]:
    paths = {
        "local_front_operator_fixed_case.json": CASE_PATH,
        "local_front_provider.py": PROVIDER_DIR / "local_front_provider.py",
        "support_fixed_case.json": SUPPORT_DIR / "fixed_case.json",
        "support_topology.py": SUPPORT_DIR / "support_topology.py",
        "role_manifest.json": SUPPORT_DIR / "role_manifest.json",
        "enriched_dof_contract.py": REGISTRY_DIR / "enriched_dof_contract.py",
        "mls3d_mls.py": MLS_DIR / "mls3d_mls.py",
        "mls3d_shape.py": MLS_DIR / "mls3d_shape.py",
        "modeiii_material_case.json": MATERIAL_DIR / "fixed_case.json",
    }
    observed = {name: _sha256(path) for name, path in paths.items()}
    mismatches = {
        name: {"expected": EXPECTED_DEPENDENCIES[name], "actual": observed[name]}
        for name in paths
        if observed[name] != EXPECTED_DEPENDENCIES[name]
    }
    if mismatches:
        raise RuntimeError(f"read-only local-front source pins mismatch: {mismatches}")
    return observed


def _load_numerical_dependencies() -> None:
    """Load MLS/provider/registry only inside the externally timed worker."""
    global np, BasisKind, ProviderStatus, TraceSide, LocalFrontProvider
    global load_pinned_inputs, validate_problem_budget

    for path in (str(PROVIDER_DIR), str(SUPPORT_DIR), str(REGISTRY_DIR), str(MLS_DIR)):
        if path not in sys.path:
            sys.path.insert(0, path)
    import numpy as np
    from enriched_dof_contract import BasisKind, ProviderStatus, validate_problem_budget
    from local_front_provider import LocalFrontProvider, TraceSide, load_pinned_inputs


def _require_fixed_case(case: dict[str, object]) -> None:
    if case.get("case_id") != CASE_ID:
        raise ValueError("fixed case ID changed")
    geometry = case["geometry"]
    provider = case["provider"]
    roles = case["source_roles_and_dofs"]
    material = case["material_and_scalar_model"]
    quadrature = case["quadrature"]
    coverage = case["support_coverage"]
    diagnostic = case["diagnostic_budget"]
    assert all(isinstance(value, dict) for value in (geometry, provider, roles, material, quadrature, coverage, diagnostic))
    expected = (
        geometry["domain_box_m"] == [[-0.01, 0.01], [-0.01, 0.01], [-0.008, 0.008]],
        geometry["crack_surface"] == "y=0, -0.01<=x<0, -0.008<=z<=0.008",
        geometry["front_segment"] == "x=0, y=0, -0.008<=z<=0.008",
        provider["module"] == "pgfrag/experiments/local_front_provider_20261006/local_front_provider.py",
        provider["provider_class"] == "LocalFrontProvider",
        provider["input_loader"] == "load_pinned_inputs",
        provider["visibility_mask"] is False,
        roles["source_count"] == 100,
        roles["total_scalar_dof_count"] == 140,
        roles["geometric_role_counts"] == {"FRONT": 30, "HEAVISIDE_ELIGIBLE": 10, "BASE": 60},
        roles["semantic_basis_counts"] == {"BASE": 100, "FRONT_F1_SHIFTED": 30, "HEAVISIDE": 10},
        roles["modeiii_analytic_coefficients_used_as_global_exact_solution"] is False,
        material["young_modulus_Pa"] == 1_000_000.0,
        material["poisson_ratio"] == 0.3,
        material["full_3d_elasticity_claim"] is False,
        quadrature["points_per_axis_per_quadrant"] == [4, 4, 4],
        quadrature["points_per_quadrant"] == 64,
        quadrature["total_volume_points"] == 256,
        quadrature["x_intervals_m"] == [[-0.01, 0.0], [0.0, 0.01]],
        quadrature["y_intervals_m"] == [[-0.01, 0.0], [0.0, 0.01]],
        quadrature["z_interval_m"] == [-0.008, 0.008],
        quadrature["quadrant_order"] == [
            "x_negative_y_negative",
            "x_negative_y_positive",
            "x_positive_y_negative",
            "x_positive_y_positive",
        ],
        quadrature["expected_full_box_volume_m3"] == 0.0000064,
        quadrature["expected_volume_weight_sum_m3"] == 0.0000064,
        quadrature["all_mapped_gauss_weights_positive"] is True,
        coverage["fixed_mls_support_radius_m"] == 0.007,
        coverage["source_lattice_axes_m"] == {
            "x": [-0.008, -0.004, 0.0, 0.004, 0.008],
            "y": [-0.008, -0.004, 0.004, 0.008],
            "z": [-0.008, -0.004, 0.0, 0.004, 0.008],
        },
        coverage["analytical_max_nearest_source_axis_offsets_m"] == [0.002, 0.004, 0.002],
        math.isclose(
            float(coverage["analytical_max_nearest_source_distance_m"]),
            0.004898979485566356,
            rel_tol=0.0,
            abs_tol=1.0e-18,
        ),
        math.isclose(
            float(coverage["analytical_min_geometric_support_coverage_margin_m"]),
            0.002101020514433644,
            rel_tol=0.0,
            abs_tol=1.0e-18,
        ),
        diagnostic["maximum_sources"] == MAX_SOURCE_COUNT,
        diagnostic["maximum_scalar_dofs"] == MAX_TOTAL_DOFS,
        diagnostic["maximum_volume_queries"] == MAX_QUADRATURE_POINTS,
        diagnostic["maximum_estimated_explicit_numeric_array_bytes"] == MAX_EXPLICIT_NUMERIC_ARRAY_BYTES,
        diagnostic["wall_timeout_seconds"] == WALL_TIMEOUT_SECONDS,
        diagnostic["run_solver"] is False,
        diagnostic["apply_load"] is False,
        diagnostic["regularize_or_truncate"] is False,
    )
    if not all(expected):
        raise ValueError("fixed local-front operator case changed a frozen input, count, or scope boundary")
    nodes = quadrature["gauss_legendre_order_4_nodes"]
    weights = quadrature["gauss_legendre_order_4_weights"]
    if len(nodes) != 4 or len(weights) != 4:
        raise ValueError("fixed quadrature rule must be GL4")
    if nodes != [-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526]:
        raise ValueError("fixed GL4 abscissae changed")
    if weights != [0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538]:
        raise ValueError("fixed GL4 weights changed")
    if not all(-1.0 < float(node) < 1.0 for node in nodes):
        raise ValueError("GL4 nodes must lie strictly inside their mapped intervals")
    if not all(float(weight) > 0.0 for weight in weights):
        raise ValueError("GL4 weights must be positive")
    if not math.isclose(sum(float(weight) for weight in weights), 2.0, rel_tol=0.0, abs_tol=2.0e-15):
        raise ValueError("GL4 weights no longer integrate a constant on [-1,1]")


def _array_budget_estimate(source_count: int, dof_count: int, query_count: int) -> dict[str, int]:
    f8 = 8
    matrix_bytes = dof_count * dof_count * f8
    parts = {
        "source_positions": source_count * 3 * f8,
        "volume_query_positions": query_count * 3 * f8,
        "volume_weights": query_count * f8,
        "mls_weights": query_count * source_count * f8,
        "mls_gradients": query_count * source_count * 3 * f8,
        "provider_values": query_count * dof_count * f8,
        "provider_gradients": query_count * dof_count * 3 * f8,
        "registry_and_probe_vectors": dof_count * 12 * f8,
        "volume_mean_rows_and_probe_fields": (dof_count * 2 + query_count * 4) * f8,
        "stiffness_local_symmetric_and_spectral_work_arrays": 16 * matrix_bytes,
        "quadrature_and_scalar_work_arrays": 8192,
    }
    return {**parts, "estimated_explicit_numeric_array_bytes": int(sum(parts.values()))}


def _maximum_nearest_axis_offset(lower: float, upper: float, source_axis: tuple[float, ...]) -> float:
    if not source_axis or tuple(sorted(set(source_axis))) != source_axis:
        raise ValueError("source lattice axes must be nonempty, sorted, and unique")
    if source_axis[0] < lower or source_axis[-1] > upper:
        raise ValueError("source lattice lies outside the fixed domain box")
    candidates = [lower, upper]
    candidates.extend(
        0.5 * (left + right)
        for left, right in zip(source_axis[:-1], source_axis[1:], strict=True)
        if lower <= 0.5 * (left + right) <= upper
    )
    return max(min(abs(point - source) for source in source_axis) for point in candidates)


def _verify_geometric_support_coverage(bundle: object, case: dict[str, object]) -> dict[str, object]:
    """Check only the nearest-source distance bound, without invoking MLS."""
    geometry = case["geometry"]
    coverage_case = case["support_coverage"]
    assert isinstance(geometry, dict) and isinstance(coverage_case, dict)
    bounds = geometry["domain_box_m"]
    axes_case = coverage_case["source_lattice_axes_m"]
    assert isinstance(axes_case, dict)
    source_positions = tuple(
        tuple(float(coordinate) for coordinate in source.position_m)
        for source in bundle.sources
    )
    if not source_positions or any(len(position) != 3 for position in source_positions):
        raise ValueError("pinned source positions are empty or malformed")
    if any(not math.isfinite(value) for position in source_positions for value in position):
        raise ValueError("pinned source lattice contains a nonfinite coordinate")
    actual_axes = tuple(tuple(sorted({position[axis] for position in source_positions})) for axis in range(3))
    declared_axes = tuple(tuple(float(value) for value in axes_case[axis]) for axis in ("x", "y", "z"))
    if actual_axes != declared_axes:
        raise ValueError("pinned source coordinate axes differ from the fixed coverage calculation")
    expected_lattice = set(itertools.product(*declared_axes))
    if len(source_positions) != len(expected_lattice) or set(source_positions) != expected_lattice:
        raise ValueError("pinned source positions do not form the declared complete tensor lattice")

    max_axis_offsets = tuple(
        _maximum_nearest_axis_offset(float(bounds[axis][0]), float(bounds[axis][1]), declared_axes[axis])
        for axis in range(3)
    )
    max_nearest_distance = math.sqrt(sum(value * value for value in max_axis_offsets))
    support_radius = float(coverage_case["fixed_mls_support_radius_m"])
    provider_support_radius = float(bundle.geometry.support_radius_m)
    if provider_support_radius != support_radius:
        raise ValueError("pinned provider support radius differs from the fixed geometric coverage radius")
    margin = support_radius - max_nearest_distance
    expected_offsets = tuple(float(value) for value in coverage_case["analytical_max_nearest_source_axis_offsets_m"])
    expected_distance = float(coverage_case["analytical_max_nearest_source_distance_m"])
    expected_margin = float(coverage_case["analytical_min_geometric_support_coverage_margin_m"])
    if (
        len(expected_offsets) != 3
        or any(
            not math.isclose(actual, expected, rel_tol=0.0, abs_tol=1.0e-18)
            for actual, expected in zip(max_axis_offsets, expected_offsets, strict=True)
        )
        or not math.isclose(max_nearest_distance, expected_distance, rel_tol=0.0, abs_tol=1.0e-18)
        or not math.isclose(margin, expected_margin, rel_tol=0.0, abs_tol=1.0e-18)
    ):
        raise ValueError("computed nearest-source geometry differs from the fixed coverage calculation")
    if margin <= 0.0:
        raise ValueError("nearest-source geometry does not establish strict support-radius coverage")
    return {
        "source_position_count": len(source_positions),
        "complete_tensor_lattice_verified": True,
        "axis_maximum_nearest_offsets_m": list(max_axis_offsets),
        "maximum_nearest_source_distance_m": max_nearest_distance,
        "support_radius_m": support_radius,
        "provider_support_radius_m": provider_support_radius,
        "minimum_geometric_support_margin_m": margin,
        "all_domain_points_strictly_inside_at_least_one_geometric_support_ball": True,
        "mls_support_cardinality_rank_conditioning_and_status_established": False,
    }


def _mapped_rule(nodes: object, weights: object, lower: float, upper: float) -> tuple[object, object]:
    half = 0.5 * (upper - lower)
    middle = 0.5 * (upper + lower)
    return middle + half * nodes, half * weights


def _build_volume_quadrature(case: dict[str, object]) -> dict[str, object]:
    quadrature = case["quadrature"]
    geometry = case["geometry"]
    assert isinstance(quadrature, dict) and isinstance(geometry, dict)
    nodes = np.asarray(quadrature["gauss_legendre_order_4_nodes"], dtype=np.float64)
    weights = np.asarray(quadrature["gauss_legendre_order_4_weights"], dtype=np.float64)
    x_bounds, y_bounds, z_bounds = geometry["domain_box_m"]
    x_intervals = quadrature["x_intervals_m"]
    y_intervals = quadrature["y_intervals_m"]
    x_parts = [_mapped_rule(nodes, weights, float(a), float(b)) for a, b in x_intervals]
    y_parts = [_mapped_rule(nodes, weights, float(a), float(b)) for a, b in y_intervals]
    z_coords, z_weights = _mapped_rule(nodes, weights, float(z_bounds[0]), float(z_bounds[1]))

    point_rows: list[tuple[float, float, float]] = []
    volume_weights: list[float] = []
    for ix_part, iy_part in itertools.product(range(2), range(2)):
        x_coords, x_weights = x_parts[ix_part]
        y_coords, y_weights = y_parts[iy_part]
        for ix, iy, iz in itertools.product(range(4), range(4), range(4)):
            x = float(x_coords[ix])
            y = float(y_coords[iy])
            z = float(z_coords[iz])
            weight = float(x_weights[ix] * y_weights[iy] * z_weights[iz])
            point_rows.append((x, y, z))
            volume_weights.append(weight)

    points = np.ascontiguousarray(np.asarray(point_rows, dtype=np.float64))
    weights_array = np.ascontiguousarray(np.asarray(volume_weights, dtype=np.float64))
    if points.shape != (256, 3) or weights_array.shape != (256,):
        raise ValueError("quadrature builder did not return the fixed 256 point batch")
    if not np.all(np.isfinite(points)) or not np.all(np.isfinite(weights_array)):
        raise ValueError("quadrature contains a nonfinite point or weight")
    if not np.all(weights_array > 0.0):
        raise ValueError("quadrature must retain strictly positive volume weights")
    if np.any(points[:, 0] == 0.0) or np.any(points[:, 1] == 0.0):
        raise ValueError("split open GL quadrants unexpectedly sampled x=0 or the y=0 crack plane")
    if np.any(points[:, 0] <= float(x_bounds[0])) or np.any(points[:, 0] >= float(x_bounds[1])):
        raise ValueError("quadrature point reached or left the open x box boundary")
    if np.any(points[:, 1] <= float(y_bounds[0])) or np.any(points[:, 1] >= float(y_bounds[1])):
        raise ValueError("quadrature point reached or left the open y box boundary")
    if np.any(points[:, 2] <= float(z_bounds[0])) or np.any(points[:, 2] >= float(z_bounds[1])):
        raise ValueError("quadrature point reached or left the open z box boundary")
    min_front_distance = float(np.min(np.hypot(points[:, 0], points[:, 1])))
    if not min_front_distance > 0.0:
        raise ValueError("quadrature sampled the exact front")

    expected_volume = float(quadrature["expected_full_box_volume_m3"])
    volume_sum = float(np.sum(weights_array))
    summation_tolerance = 64.0 * float(np.finfo(np.float64).eps) * expected_volume
    if abs(volume_sum - expected_volume) > summation_tolerance:
        raise ValueError("GL4 tensor rule does not recover the fixed box volume to roundoff")
    if len({tuple(float(value) for value in row) for row in points}) != len(points):
        raise ValueError("quadrant quadrature contains duplicate physical points")
    return {
        "points": points,
        "weights": weights_array,
        "volume_sum_m3": volume_sum,
        "expected_volume_m3": expected_volume,
        "minimum_xy_distance_to_front_m": min_front_distance,
        "point_count": int(points.shape[0]),
    }


def _status_counts(statuses: object) -> dict[str, int]:
    names = sorted({str(status.value) for status in statuses})
    return {name: sum(status.value == name for status in statuses) for name in names}


def _json_finite(value: float) -> float | None:
    return float(value) if math.isfinite(float(value)) else None


def _semantic_order(registry: object) -> list[dict[str, object]]:
    ids = registry.canonical_ids
    rows: list[dict[str, object]] = []
    for slot, (key, canonical_id) in enumerate(zip(registry.ordered_keys, ids, strict=True)):
        rows.append(
            {
                "slot": slot,
                "canonical_id": canonical_id,
                "source_id": key.source_id,
                "field_component": key.field_component.value,
                "basis_kind": key.basis_kind.value,
                "basis_index": key.basis_index,
            }
        )
    return rows


def _blocked_result(
    *,
    status: str,
    reason: str,
    case: dict[str, object],
    quadrature: dict[str, object],
    provider: object | None,
    budget: dict[str, int],
    dependency_hashes: dict[str, str],
    geometric_support_coverage: dict[str, object],
    started_at: str,
    started: float,
) -> dict[str, object]:
    unknown_rows: list[dict[str, object]] = []
    statuses: dict[str, int] = {}
    if provider is not None:
        statuses = _status_counts(provider.row_statuses)
        for index, row_status in enumerate(provider.row_statuses):
            if row_status is not ProviderStatus.OK:
                unknown_rows.append(
                    {
                        "query_index": index,
                        "coordinate_xyz_m": [float(value) for value in provider.query_points_m[index]],
                        "provider_status": row_status.value,
                        "original_reason": provider.row_diagnostics[index],
                        "geometry_classification_reason": provider.classifications[index].reason,
                    }
                )
    return {
        "schema": "local_front_scalar_operator_diagnostic.v1",
        "case_id": CASE_ID,
        "status": status,
        "classification": "SOURCE_ONLY_SCALAR_ANTIPLANE_POISSON_DIAGNOSTIC",
        "started_at_utc": started_at,
        "fixed_case": case,
        "blocked_reason": reason,
        "dimensions": {"source_count": 100, "scalar_dof_count": 140, "volume_query_count": 256},
        "array_budget": {
            "parts_bytes": budget,
            "estimated_explicit_numeric_array_bytes": budget.get("estimated_explicit_numeric_array_bytes", 0),
            "cap_bytes": MAX_EXPLICIT_NUMERIC_ARRAY_BYTES,
            "within_estimate_cap": budget.get("estimated_explicit_numeric_array_bytes", 0) <= MAX_EXPLICIT_NUMERIC_ARRAY_BYTES,
            "python_object_and_vendor_workspace_excluded": True,
        },
        "provider_row_status_counts": statuses,
        "unknown_rows_with_unmodified_coordinates_and_reason": unknown_rows,
        "geometric_support_coverage": geometric_support_coverage,
        "quadrature_coordinates_xyz_m": quadrature.get("points").tolist() if quadrature.get("points") is not None else [],
        "no_partial_stiffness_returned": True,
        "dependency_sha256": dependency_hashes,
        "elapsed_seconds": time.monotonic() - started,
        "wall_timeout_seconds": WALL_TIMEOUT_SECONDS,
        "finished_at_utc": _utc_now(),
        "limitations": case["limitations"],
    }


def _run_diagnostic() -> dict[str, object]:
    started = time.monotonic()
    started_at = _utc_now()
    dependency_hashes = _verify_dependencies()
    case = json.loads(CASE_PATH.read_text(encoding="utf-8"))
    _require_fixed_case(case)
    roles = case["source_roles_and_dofs"]
    quadrature_case = case["quadrature"]
    material_case = case["material_and_scalar_model"]
    assert isinstance(roles, dict) and isinstance(quadrature_case, dict) and isinstance(material_case, dict)
    source_count = int(roles["source_count"])
    dof_count = int(roles["total_scalar_dof_count"])
    query_count = int(quadrature_case["total_volume_points"])
    if not (1 <= source_count <= MAX_SOURCE_COUNT):
        raise DiagnosticBudgetError("source count exceeds the hard source cap")
    if not (1 <= dof_count <= MAX_TOTAL_DOFS):
        raise DiagnosticBudgetError("scalar DOF count exceeds the hard DOF cap")
    if not (1 <= query_count <= MAX_QUADRATURE_POINTS):
        raise DiagnosticBudgetError("volume query count exceeds the hard quadrature cap")
    budget = _array_budget_estimate(source_count, dof_count, query_count)
    if budget["estimated_explicit_numeric_array_bytes"] > MAX_EXPLICIT_NUMERIC_ARRAY_BYTES:
        raise DiagnosticBudgetError("estimated explicit numeric arrays exceed the fixed 8 MiB cap")

    # Only after the scalar count and byte preflight do we import/load the
    # pinned provider and allocate query, MLS, and operator arrays.
    _load_numerical_dependencies()
    validate_problem_budget(source_count, dof_count, query_count)
    bundle = load_pinned_inputs()
    if len(bundle.sources) != source_count:
        raise ValueError("pinned provider source count does not match this case")
    geometric_support_coverage = _verify_geometric_support_coverage(bundle, case)
    quadrature = _build_volume_quadrature(case)
    points_array = quadrature["points"]
    weights_array = quadrature["weights"]
    assert isinstance(points_array, np.ndarray) and isinstance(weights_array, np.ndarray)
    if int(points_array.shape[0]) != query_count:
        raise ValueError("fixed quadrature count differs from the preflight count")

    query_points = tuple(tuple(float(value) for value in row) for row in points_array)
    trace_sides = (None,) * query_count
    provider = LocalFrontProvider(bundle, query_points, trace_sides)
    if time.monotonic() - started > WALL_TIMEOUT_SECONDS:
        raise DiagnosticBudgetError("60-second worker deadline exceeded during local-front MLS evaluation")
    if provider.registry.dof_count != dof_count or provider.registry.source_count != source_count:
        raise ValueError("provider's frozen semantic registry differs from the fixed case")

    classifications = tuple(item.status for item in provider.classifications)
    if any(status is not ProviderStatus.OK for status in classifications) or provider.unknown_query_count:
        return _blocked_result(
            status="BLOCKED_UNKNOWN_GEOMETRY_OR_MLS_SUPPORT",
            reason="at least one quadrature row is UNKNOWN; the original query coordinate and reason are retained",
            case=case,
            quadrature=quadrature,
            provider=provider,
            budget=budget,
            dependency_hashes=dependency_hashes,
            geometric_support_coverage=geometric_support_coverage,
            started_at=started_at,
            started=started,
        )
    if provider.mls_result is None or len(provider.mls_result.statuses) != query_count:
        raise RuntimeError("valid fixed volume points did not produce one MLS result row each")
    if any(status.value != "OK" for status in provider.mls_result.statuses):
        return _blocked_result(
            status="BLOCKED_UNKNOWN_GEOMETRY_OR_MLS_SUPPORT",
            reason="MLS returned a non-OK support, rank, conditioning, or solve status",
            case=case,
            quadrature=quadrature,
            provider=provider,
            budget=budget,
            dependency_hashes=dependency_hashes,
            geometric_support_coverage=geometric_support_coverage,
            started_at=started_at,
            started=started,
        )
    if provider.ok_query_count != query_count:
        raise RuntimeError("provider row success count differs from the complete quadrature batch")

    registry = provider.registry
    canonical_ids = registry.canonical_ids
    semantic_order = _semantic_order(registry)
    if len(canonical_ids) != dof_count or len(set(canonical_ids)) != dof_count:
        raise RuntimeError("registry semantic IDs are missing or duplicated")
    mu = float(material_case["young_modulus_Pa"]) / (2.0 * (1.0 + float(material_case["poisson_ratio"])))
    if not math.isfinite(mu) or mu <= 0.0:
        raise DiagnosticNumericsError("scalar Poisson shear coefficient must be finite and positive")

    stiffness = np.zeros((dof_count, dof_count), dtype=np.float64)
    mean_integral_row = np.zeros(dof_count, dtype=np.float64)
    affine_coefficients = np.zeros(dof_count, dtype=np.float64)
    translation_coefficients = np.zeros(dof_count, dtype=np.float64)
    source_position_by_id = {source.source_id: source.position_m for source in bundle.sources}
    for slot, key in enumerate(registry.ordered_keys):
        if key.basis_kind is BasisKind.BASE:
            affine_coefficients[slot] = float(source_position_by_id[key.source_id][0])
            translation_coefficients[slot] = 1.0
        else:
            affine_coefficients[slot] = 0.0
            translation_coefficients[slot] = 0.0

    affine_value_rows = np.empty(query_count, dtype=np.float64)
    affine_gradient_rows = np.empty((query_count, 3), dtype=np.float64)
    for query_index in range(query_count):
        sample = provider.sample_at(query_index)
        if sample.status is not ProviderStatus.OK:
            return _blocked_result(
                status="BLOCKED_UNKNOWN_PROVIDER_SAMPLE",
                reason=sample.diagnostic or sample.status.value,
                case=case,
                quadrature=quadrature,
                provider=provider,
                budget=budget,
                dependency_hashes=dependency_hashes,
                geometric_support_coverage=geometric_support_coverage,
                started_at=started_at,
                started=started,
            )
        if sample.coefficient_ids != canonical_ids or len(sample.contributions) != dof_count:
            raise RuntimeError(f"provider semantic row differs from canonical registry order at q={query_index}")
        values = np.asarray([item.value for item in sample.contributions], dtype=np.float64)
        gradients = np.asarray([item.gradient for item in sample.contributions], dtype=np.float64)
        if values.shape != (dof_count,) or gradients.shape != (dof_count, 3):
            raise RuntimeError(f"provider row has a malformed value/gradient shape at q={query_index}")
        if not np.all(np.isfinite(values)) or not np.all(np.isfinite(gradients)):
            raise DiagnosticNumericsError(f"provider row is nonfinite at q={query_index}")

        weight = float(weights_array[query_index])
        with np.errstate(over="ignore", invalid="ignore"):
            local = gradients @ gradients.T
            local *= mu * weight
            stiffness += local
            mean_integral_row += values * weight
            affine_value_rows[query_index] = float(values @ affine_coefficients)
            affine_gradient_rows[query_index, :] = gradients.T @ affine_coefficients
        if not np.all(np.isfinite(local)) or not np.all(np.isfinite(stiffness)):
            raise DiagnosticNumericsError(f"scalar stiffness became nonfinite at q={query_index}")
        if not np.all(np.isfinite(mean_integral_row)):
            raise DiagnosticNumericsError(f"volume mean row became nonfinite at q={query_index}")
        if not np.all(np.isfinite(affine_gradient_rows[query_index])) or not math.isfinite(
            float(affine_value_rows[query_index])
        ):
            raise DiagnosticNumericsError(f"BASE affine patch became nonfinite at q={query_index}")

    volume_sum = float(quadrature["volume_sum_m3"])
    mean_row = mean_integral_row / volume_sum
    if not np.all(np.isfinite(mean_row)):
        raise DiagnosticNumericsError("normalized volume mean row contains a nonfinite entry")
    stiffness.setflags(write=False)

    constant_residual = stiffness @ translation_coefficients
    stiffness_frobenius = float(np.linalg.norm(stiffness, ord="fro"))
    translation_residual = float(
        np.linalg.norm(constant_residual)
        / max(stiffness_frobenius * np.linalg.norm(translation_coefficients), np.finfo(np.float64).tiny)
    )
    symmetry_relative = float(
        np.linalg.norm(stiffness - stiffness.T, ord="fro")
        / max(stiffness_frobenius, np.finfo(np.float64).tiny)
    )
    symmetric_part = 0.5 * (stiffness + stiffness.T)
    eigenvalues = np.linalg.eigvalsh(symmetric_part)
    if not np.all(np.isfinite(eigenvalues)):
        raise DiagnosticNumericsError("symmetric-part eigenvalue spectrum is nonfinite")
    if not all(
        math.isfinite(value)
        for value in (stiffness_frobenius, translation_residual, symmetry_relative)
    ):
        raise DiagnosticNumericsError("operator norm or residual diagnostic is nonfinite")

    epsilon = float(np.finfo(np.float64).eps)
    unit_roundoff = epsilon / 2.0
    operation_count = 12 + query_count + 2 * dof_count + dof_count * dof_count
    product = operation_count * unit_roundoff
    gamma_n = product / (1.0 - product)
    eigen_scale = max(float(np.max(np.abs(eigenvalues))), np.finfo(np.float64).tiny)
    zero_tolerance = 4.0 * gamma_n * eigen_scale
    near_zero_count = int(np.count_nonzero(np.abs(eigenvalues) <= zero_tolerance))
    negative_count = int(np.count_nonzero(eigenvalues < 0.0))
    negative_beyond_roundoff = int(np.count_nonzero(eigenvalues < -zero_tolerance))
    spectrum_invalid_reasons = []
    if negative_beyond_roundoff:
        spectrum_invalid_reasons.append("significant_negative_eigenvalues")
    if near_zero_count != 1:
        spectrum_invalid_reasons.append("near_zero_count_differs_from_single_scalar_translation")
    mechanical_eligibility = (
        "BLOCKED_INVALID_SPECTRUM" if spectrum_invalid_reasons else "SPECTRUM_CHECK_PASS_ONLY_UNQUALIFIED"
    )
    solve_allowed = mechanical_eligibility == "SPECTRUM_CHECK_PASS_ONLY_UNQUALIFIED"
    assembly_status = "OK"
    result_status = (
        "ASSEMBLED_INVALID_SPECTRUM_DIAGNOSTIC_ONLY"
        if spectrum_invalid_reasons
        else "ASSEMBLED_DIAGNOSTIC_ONLY_UNQUALIFIED"
    )

    x_reference = points_array[:, 0]
    exact_affine_gradient = np.asarray((1.0, 0.0, 0.0), dtype=np.float64)
    affine_displacement_error = np.abs(affine_value_rows - x_reference)
    affine_gradient_error = np.linalg.norm(affine_gradient_rows - exact_affine_gradient[None, :], axis=1)
    affine_discrete_energy = 0.5 * float(affine_coefficients @ stiffness @ affine_coefficients)
    exact_tensor_energy = 0.5 * mu * volume_sum
    exact_box_energy = 0.5 * mu * float(quadrature["expected_volume_m3"])
    if not np.all(np.isfinite((affine_discrete_energy, exact_tensor_energy, exact_box_energy))):
        raise DiagnosticNumericsError("affine probe energy is nonfinite")
    normalized_mean_row = mean_row
    translated_mean = float(normalized_mean_row @ translation_coefficients)
    stiffness_bytes = np.asarray(stiffness, dtype="<f8", order="C").tobytes(order="C")

    mls_result = provider.mls_result
    support_counts = np.asarray(mls_result.support_counts)
    condition_numbers = np.asarray(mls_result.condition_numbers, dtype=np.float64)
    finite_conditions = condition_numbers[np.isfinite(condition_numbers)]
    result: dict[str, object] = {
        "schema": "local_front_scalar_operator_diagnostic.v1",
        "case_id": CASE_ID,
        "status": result_status,
        "assembly_status": assembly_status,
        "mechanical_eligibility": mechanical_eligibility,
        "mechanical_eligibility_reasons": spectrum_invalid_reasons,
        "solve_allowed": solve_allowed,
        "solver_invoked": False,
        "classification": "SOURCE_ONLY_SCALAR_ANTIPLANE_POISSON_DIAGNOSTIC",
        "started_at_utc": started_at,
        "fixed_case": case,
        "dimensions": {
            "source_count": source_count,
            "semantic_scalar_dof_count": dof_count,
            "quadrature_point_count": query_count,
            "K_shape": list(stiffness.shape),
            "K_units": "N/m",
        },
        "registry": {
            "state": registry.state.value,
            "source_count": registry.source_count,
            "dof_count": registry.dof_count,
            "semantic_order_source": "registry.ordered_keys and registry.canonical_ids",
            "canonical_ids": list(canonical_ids),
            "ordered_semantic_keys": semantic_order,
            "basis_key_counts": {
                kind.value: sum(key.basis_kind is kind for key in registry.ordered_keys)
                for kind in BasisKind
            },
        },
        "material_and_model": {
            "young_modulus_Pa": float(material_case["young_modulus_Pa"]),
            "poisson_ratio": float(material_case["poisson_ratio"]),
            "mu_Pa": mu,
            "bilinear_form": "integral mu * grad(N_a) dot grad(N_b) dV",
            "scalar_antiplane_poisson_only": True,
            "full_3d_isotropic_elasticity_claim": False,
            "full_3d_epsilon_zz_lambda_plus_2mu_term_used": False,
        },
        "support_and_mls": {
            "status_counts": _status_counts(provider.row_statuses),
            "mls_status_counts": _status_counts(mls_result.statuses),
            "support_count_min": int(np.min(support_counts)),
            "support_count_max": int(np.max(support_counts)),
            "support_counts_by_query": [int(value) for value in support_counts],
            "condition_number_min": float(np.min(finite_conditions)) if finite_conditions.size else None,
            "condition_number_max": float(np.max(finite_conditions)) if finite_conditions.size else None,
            "condition_numbers_by_query": [_json_finite(float(value)) for value in condition_numbers],
            "unknown_rows": [],
            "geometric_coverage_does_not_prove_mls_rank_or_conditioning": True,
        },
        "geometric_support_coverage": geometric_support_coverage,
        "quadrature": {
            "rule": case["quadrature"]["rule"],
            "query_order": "four xy quadrants, GL4 x GL4 x GL4, 64 rows per quadrant",
            "coordinates_xyz_m": points_array.tolist(),
            "all_volume_weights_positive": bool(np.all(weights_array > 0.0)),
            "volume_weights_m3": weights_array.tolist(),
            "volume_weight_sum_m3": volume_sum,
            "expected_box_volume_m3": float(quadrature["expected_volume_m3"]),
            "volume_sum_abs_error_m3": abs(volume_sum - float(quadrature["expected_volume_m3"])),
            "minimum_xy_distance_to_front_m": float(quadrature["minimum_xy_distance_to_front_m"]),
            "query_on_crack_plane_or_front": False,
            "endpoint_or_epsilon_shift_used": False,
            "integration_accuracy_claim": "single fixed diagnostic GL4 split rule; no convergence or singular-integration qualification",
        },
        "volume_mean_row": {
            "integrated_shape_row_m3": mean_integral_row.tolist(),
            "normalized_mean_shape_row": normalized_mean_row.tolist(),
            "units_of_normalized_row": "dimensionless",
            "semantic_order_matches_registry": True,
            "mean_of_constant_translation_m": translated_mean,
        },
        "constant_translation": {
            "coefficient_vector_m": translation_coefficients.tolist(),
            "base_entries_are_one_meter": True,
            "all_enrichment_entries_are_zero": True,
            "expected_scalar_nullity_at_least": 1,
            "scaled_Kt_residual": translation_residual,
            "interpretation": "one scalar antiplane translation; no six-vector rigid-body claim",
        },
        "affine_base_patch": {
            "field": "Uz=x",
            "coefficient_vector_m": affine_coefficients.tolist(),
            "base_coefficients_equal_source_x_m": True,
            "enrichment_coefficients_zero": True,
            "max_displacement_error_m": float(np.max(affine_displacement_error)),
            "max_gradient_error_per_m": float(np.max(affine_gradient_error)),
            "discrete_K_energy_J": affine_discrete_energy,
            "independent_exact_scalar_tensor_energy_J": exact_tensor_energy,
            "independent_closed_box_energy_J": exact_box_energy,
            "K_energy_relative_difference_vs_tensor_rule": abs(affine_discrete_energy - exact_tensor_energy) / max(abs(exact_tensor_energy), np.finfo(np.float64).tiny),
            "tensor_rule_relative_difference_vs_box_formula": abs(exact_tensor_energy - exact_box_energy) / max(abs(exact_box_energy), np.finfo(np.float64).tiny),
            "is_ModeIII_analytic_coefficient_vector": False,
            "acceptance_gate": False,
        },
        "operator_diagnostics": {
            "stiffness_sha256_float64_c_order": hashlib.sha256(stiffness_bytes).hexdigest(),
            "stiffness_matrix_N_per_m": stiffness.tolist(),
            "frobenius_norm_N_per_m": stiffness_frobenius,
            "relative_frobenius_asymmetry": symmetry_relative,
            "eigenvalue_matrix": "symmetric_part=(K+K.T)/2 for eigenvalues only; raw assembled K remains unchanged",
            "eigenvalues_N_per_m": [float(value) for value in eigenvalues],
            "raw_negative_eigenvalue_count": negative_count,
            "negative_beyond_roundoff_count": negative_beyond_roundoff,
            "near_zero_count": near_zero_count,
            "expected_scalar_translation_mode_count": 1,
            "additional_near_zero_after_translation": max(near_zero_count - 1, 0),
            "absolute_near_zero_threshold_N_per_m": zero_tolerance,
            "roundoff_reference": {
                "binary64_epsilon": epsilon,
                "unit_roundoff": unit_roundoff,
                "conservative_operation_count": operation_count,
                "gamma_n": gamma_n,
                "zero_threshold_multiplier": 4.0,
                "acceptance_gate": False,
            },
            "no_regularization_or_rank_truncation": True,
            "stability_claim": False,
            "mechanical_eligibility": mechanical_eligibility,
            "mechanical_eligibility_reasons": spectrum_invalid_reasons,
            "solve_allowed": solve_allowed,
            "solver_invoked": False,
        },
        "array_budget": {
            "parts_bytes": budget,
            "estimated_explicit_numeric_array_bytes": int(budget["estimated_explicit_numeric_array_bytes"]),
            "cap_bytes": MAX_EXPLICIT_NUMERIC_ARRAY_BYTES,
            "within_estimate_cap": int(budget["estimated_explicit_numeric_array_bytes"]) <= MAX_EXPLICIT_NUMERIC_ARRAY_BYTES,
            "python_object_and_vendor_workspace_excluded": True,
        },
        "source_sha256": dependency_hashes,
        "elapsed_seconds": time.monotonic() - started,
        "wall_timeout_seconds": WALL_TIMEOUT_SECONDS,
        "within_worker_wall_timeout": time.monotonic() - started <= WALL_TIMEOUT_SECONDS,
        "limitations": case["limitations"],
        "finished_at_utc": _utc_now(),
    }
    if not np.all(np.isfinite(stiffness)):
        raise DiagnosticNumericsError("assembled K contains a nonfinite value")
    return result


def _require_finite_json_tree(value: object, location: str = "$") -> None:
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError(f"nonfinite numeric value in worker JSON at {location}")
    elif isinstance(value, int) and not isinstance(value, bool):
        if value < -sys.float_info.max or value > sys.float_info.max:
            raise ValueError(f"out-of-range integer in worker JSON at {location}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _require_finite_json_tree(child, f"{location}[{index}]")
    elif isinstance(value, dict):
        for key, child in value.items():
            _require_finite_json_tree(child, f"{location}.{key}")


def _require_finite_vector(value: object, length: int, location: str) -> list[float]:
    if not isinstance(value, list) or len(value) != length:
        raise ValueError(f"worker result {location} must be a list of length {length}")
    converted: list[float] = []
    for index, item in enumerate(value):
        if isinstance(item, bool) or not isinstance(item, (int, float)) or not math.isfinite(float(item)):
            raise ValueError(f"worker result {location}[{index}] is not a finite number")
        converted.append(float(item))
    return converted


def _validate_complete_worker_result(result: object) -> dict[str, object]:
    if not isinstance(result, dict):
        raise ValueError("worker result root must be a JSON object")
    _require_finite_json_tree(result)
    if result.get("schema") != "local_front_scalar_operator_diagnostic.v1":
        raise ValueError("worker result has an unexpected schema")
    if result.get("case_id") != CASE_ID:
        raise ValueError("worker result case ID differs from the fixed case")
    expected_case = json.loads(CASE_PATH.read_text(encoding="utf-8"))
    if result.get("fixed_case") != expected_case:
        raise ValueError("worker result fixed case differs from the current pinned case file")
    if (
        result.get("assembly_status") != "OK"
        or result.get("solve_allowed") is not True
        or result.get("solver_invoked") is not False
    ):
        raise ValueError("worker result does not describe a complete assembly-only diagnostic")
    status = result.get("status")
    eligibility = result.get("mechanical_eligibility")
    reasons = result.get("mechanical_eligibility_reasons")
    if not isinstance(reasons, list):
        raise ValueError("worker result mechanical eligibility reasons are malformed")
    if status != "ASSEMBLED_DIAGNOSTIC_ONLY_UNQUALIFIED" or eligibility != "SPECTRUM_CHECK_PASS_ONLY_UNQUALIFIED" or reasons:
        raise ValueError("worker result spectrum did not pass the pre-promotion eligibility check")

    source_hashes = result.get("source_sha256")
    worker_hashes = result.get("worker_dependency_sha256")
    if source_hashes != EXPECTED_DEPENDENCIES or worker_hashes != EXPECTED_DEPENDENCIES:
        raise ValueError("worker result dependency hashes do not match the frozen source inputs")
    dimensions = result.get("dimensions")
    if not isinstance(dimensions, dict) or dimensions != {
        "source_count": 100,
        "semantic_scalar_dof_count": 140,
        "quadrature_point_count": 256,
        "K_shape": [140, 140],
        "K_units": "N/m",
    }:
        raise ValueError("worker result source, DOF, query, or stiffness dimensions differ from the fixed case")

    registry = result.get("registry")
    if not isinstance(registry, dict):
        raise ValueError("worker result semantic registry section is missing")
    canonical_ids = registry.get("canonical_ids")
    ordered_keys = registry.get("ordered_semantic_keys")
    if (
        not isinstance(canonical_ids, list)
        or len(canonical_ids) != 140
        or len(set(canonical_ids)) != 140
        or not isinstance(ordered_keys, list)
        or len(ordered_keys) != 140
    ):
        raise ValueError("worker result semantic registry does not contain 140 unique ordered IDs")
    if registry.get("basis_key_counts") != {"BASE": 100, "FRONT_F1_SHIFTED": 30, "HEAVISIDE": 10}:
        raise ValueError("worker result semantic basis counts differ from the fixed 100/30/10 layout")

    quadrature = result.get("quadrature")
    if not isinstance(quadrature, dict):
        raise ValueError("worker result quadrature section is missing")
    coordinates = quadrature.get("coordinates_xyz_m")
    if not isinstance(coordinates, list) or len(coordinates) != 256:
        raise ValueError("worker result quadrature coordinate count is incomplete")
    for index, point in enumerate(coordinates):
        _require_finite_vector(point, 3, f"quadrature.coordinates_xyz_m[{index}]")
    volume_weights = _require_finite_vector(quadrature.get("volume_weights_m3"), 256, "quadrature.volume_weights_m3")
    if any(weight <= 0.0 for weight in volume_weights):
        raise ValueError("worker result contains a nonpositive volume quadrature weight")
    if not math.isclose(sum(volume_weights), 6.4e-6, rel_tol=0.0, abs_tol=64.0 * sys.float_info.epsilon * 6.4e-6):
        raise ValueError("worker result quadrature weights do not recover the fixed box volume")

    mean = result.get("volume_mean_row")
    if not isinstance(mean, dict):
        raise ValueError("worker result volume mean rows are missing")
    _require_finite_vector(mean.get("integrated_shape_row_m3"), 140, "volume_mean_row.integrated_shape_row_m3")
    _require_finite_vector(mean.get("normalized_mean_shape_row"), 140, "volume_mean_row.normalized_mean_shape_row")

    operator = result.get("operator_diagnostics")
    if not isinstance(operator, dict):
        raise ValueError("worker result operator diagnostics are missing")
    matrix = operator.get("stiffness_matrix_N_per_m")
    if not isinstance(matrix, list) or len(matrix) != 140:
        raise ValueError("worker result stiffness matrix row count is incomplete")
    matrix_bytes = bytearray()
    for row_index, row in enumerate(matrix):
        numeric_row = _require_finite_vector(row, 140, f"operator_diagnostics.stiffness_matrix_N_per_m[{row_index}]")
        matrix_bytes.extend(struct.pack("<140d", *numeric_row))
    if hashlib.sha256(matrix_bytes).hexdigest() != operator.get("stiffness_sha256_float64_c_order"):
        raise ValueError("worker result stiffness matrix hash does not match its serialized values")
    _require_finite_vector(operator.get("eigenvalues_N_per_m"), 140, "operator_diagnostics.eigenvalues_N_per_m")
    if (
        operator.get("solve_allowed") is not True
        or operator.get("solver_invoked") is not False
        or operator.get("mechanical_eligibility") != eligibility
    ):
        raise ValueError("operator diagnostic solve gate does not match the top-level eligibility")
    negative_beyond_roundoff = operator.get("negative_beyond_roundoff_count")
    near_zero_count = operator.get("near_zero_count")
    if type(negative_beyond_roundoff) is not int or type(near_zero_count) is not int:
        raise ValueError("worker result spectrum counts are malformed")
    if (
        negative_beyond_roundoff != 0
        or near_zero_count != 1
        or operator.get("expected_scalar_translation_mode_count") != 1
    ):
        raise ValueError("worker result failed the single-null-mode/no-significant-negative pre-promotion gate")
    if operator.get("additional_near_zero_after_translation") != 0:
        raise ValueError("worker result extra-null count is inconsistent with the full spectrum")

    support = result.get("support_and_mls")
    if not isinstance(support, dict) or support.get("unknown_rows") != []:
        raise ValueError("worker result is missing a complete known MLS support batch")
    counts = support.get("support_counts_by_query")
    conditions = support.get("condition_numbers_by_query")
    if not isinstance(counts, list) or len(counts) != 256 or not isinstance(conditions, list) or len(conditions) != 256:
        raise ValueError("worker result support diagnostics do not cover all 256 query rows")
    if any(type(value) is not int or value <= 0 or value > 100 for value in counts):
        raise ValueError("worker result contains an invalid support count")
    if any(
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(float(value))
        or float(value) < 1.0
        for value in conditions
    ):
        raise ValueError("worker result contains an invalid MLS condition number")
    if support.get("status_counts") != {"OK": 256} or support.get("mls_status_counts") != {"OK": 256}:
        raise ValueError("worker result does not contain 256 known provider and MLS rows")

    budget = result.get("array_budget")
    expected_budget = _array_budget_estimate(100, 140, 256)
    if (
        not isinstance(budget, dict)
        or budget.get("within_estimate_cap") is not True
        or budget.get("estimated_explicit_numeric_array_bytes") != expected_budget["estimated_explicit_numeric_array_bytes"]
        or budget.get("parts_bytes") != expected_budget
        or budget.get("cap_bytes") != MAX_EXPLICIT_NUMERIC_ARRAY_BYTES
    ):
        raise ValueError("worker result explicit numeric-array estimate exceeded its fixed cap")
    coverage = result.get("geometric_support_coverage")
    if (
        not isinstance(coverage, dict)
        or coverage.get("source_position_count") != 100
        or coverage.get("complete_tensor_lattice_verified") is not True
        or coverage.get("all_domain_points_strictly_inside_at_least_one_geometric_support_ball") is not True
        or coverage.get("mls_support_cardinality_rank_conditioning_and_status_established") is not False
    ):
        raise ValueError("worker result geometric coverage record is missing or overclaims MLS status")
    for field, expected_value in (
        ("maximum_nearest_source_distance_m", 0.004898979485566356),
        ("support_radius_m", 0.007),
        ("provider_support_radius_m", 0.007),
        ("minimum_geometric_support_margin_m", 0.002101020514433644),
    ):
        actual_value = coverage.get(field)
        if not isinstance(actual_value, (int, float)) or not math.isclose(
            float(actual_value), expected_value, rel_tol=0.0, abs_tol=1.0e-18
        ):
            raise ValueError(f"worker result geometric coverage field {field} differs from the frozen case")
    return result


def _write_unknown_receipt(path: Path, receipt: dict[str, object]) -> None:
    with path.open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(receipt, indent=2, ensure_ascii=False, allow_nan=False) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=HERE / "results" / "local_front_scalar_operator.json",
        help="primary JSON output path",
    )
    parser.add_argument("--_worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    output_path = args.output.resolve()

    if args._worker:
        if output_path.exists():
            print(json.dumps({"status": "OUTPUT_EXISTS_NO_OVERWRITE", "path": str(output_path)}))
            return 3
        try:
            dependency_hashes = _verify_dependencies()
            result = _run_diagnostic()
            result["worker_dependency_sha256"] = dependency_hashes
            if result.get("assembly_status") != "OK":
                exit_code = 2
            elif result.get("mechanical_eligibility") == "BLOCKED_INVALID_SPECTRUM":
                exit_code = 3
            else:
                exit_code = 0
        except Exception as exc:
            result = {
                "schema": "local_front_scalar_operator_diagnostic.v1",
                "case_id": CASE_ID,
                "status": "DIAGNOSTIC_ERROR",
                "classification": "SOURCE_ONLY_SCALAR_ANTIPLANE_POISSON_DIAGNOSTIC",
                "error_type": type(exc).__name__,
                "error": str(exc),
                "traceback": traceback.format_exc(limit=8),
                "no_partial_stiffness_returned": True,
                "no_retry": True,
                "finished_at_utc": _utc_now(),
            }
            exit_code = 2
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
        print(
            json.dumps(
                {
                    "case_id": result.get("case_id"),
                    "status": result.get("status"),
                    "dimensions": result.get("dimensions"),
                    "blocked_reason": result.get("blocked_reason"),
                    "error_type": result.get("error_type"),
                    "error": result.get("error"),
                    "finished_at_utc": result.get("finished_at_utc"),
                },
                ensure_ascii=False,
                allow_nan=False,
            )
        )
        return exit_code

    # This parent timeout encloses import, MLS evaluation, assembly, spectrum,
    # and output generation in one child process. No retry is performed.
    receipt_path = output_path.with_name(f"{output_path.stem}.unknown_receipt.json")
    staging = output_path.with_name(f".{output_path.name}.{os.getpid()}.worker.json")
    promotion_path = output_path.with_name(f".{output_path.name}.{os.getpid()}.promotion.json")
    if output_path.exists() or receipt_path.exists() or staging.exists() or promotion_path.exists():
        print(
            json.dumps(
                {
                    "status": "OUTPUT_OR_STAGING_EXISTS_NO_OVERWRITE",
                    "primary_output": str(output_path),
                    "unknown_receipt": str(receipt_path),
                    "staging": str(staging),
                },
                ensure_ascii=False,
            )
        )
        return 3
    output_path.parent.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    started_at = _utc_now()
    argv = [sys.executable, "-B", str(HERE / "operator_diagnostic.py"), "--_worker", "--output", str(staging)]
    stdout = ""
    stderr = ""
    return_code = 125
    timed_out = False
    worker_result: object | None = None
    failure_reason: str | None = None
    try:
        completed = subprocess.run(
            argv,
            cwd=str(REPO),
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=WALL_TIMEOUT_SECONDS,
            check=False,
        )
        stdout, stderr = completed.stdout, completed.stderr
        return_code = int(completed.returncode)
        elapsed_before_validation = time.monotonic() - started
        if elapsed_before_validation > WALL_TIMEOUT_SECONDS:
            timed_out = True
            return_code = 124
            failure_reason = "parent wall deadline exceeded before promotion"
        elif return_code != 0:
            failure_reason = f"worker returned nonzero status {return_code}"
        else:
            try:
                raw_worker_result = json.loads(staging.read_text(encoding="utf-8"))
                worker_result = _validate_complete_worker_result(raw_worker_result)
                if time.monotonic() - started > WALL_TIMEOUT_SECONDS:
                    timed_out = True
                    failure_reason = "parent wall deadline exceeded during result validation"
                    return_code = 124
                    worker_result = None
            except (OSError, json.JSONDecodeError, ValueError, TypeError, OverflowError, struct.error) as exc:
                failure_reason = f"worker result validation failed: {type(exc).__name__}: {exc}"
                return_code = 125
    except subprocess.TimeoutExpired as exc:
        timed_out = True
        stdout_value = exc.stdout or ""
        stderr_value = exc.stderr or ""
        stdout = stdout_value.decode("utf-8", "replace") if isinstance(stdout_value, bytes) else stdout_value
        stderr = stderr_value.decode("utf-8", "replace") if isinstance(stderr_value, bytes) else stderr_value
        return_code = 124
        failure_reason = "worker exceeded the 60-second parent timeout"
    except OSError as exc:
        failure_reason = f"worker launch or file operation failed: {type(exc).__name__}: {exc}"
        return_code = 125

    elapsed = time.monotonic() - started
    if failure_reason is None and worker_result is not None:
        worker_result["parent_argv"] = argv
        worker_result["parent_started_at_utc"] = started_at
        worker_result["parent_elapsed_seconds"] = elapsed
        worker_result["parent_wall_timeout_seconds"] = WALL_TIMEOUT_SECONDS
        worker_result["parent_within_wall_timeout"] = elapsed <= WALL_TIMEOUT_SECONDS
        worker_result["parent_timed_out"] = False
        worker_result["parent_return_code"] = return_code
        worker_result["parent_stdout"] = stdout
        worker_result["parent_stderr"] = stderr
        try:
            with promotion_path.open("x", encoding="utf-8", newline="\n") as stream:
                stream.write(json.dumps(worker_result, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
            promotion_path.rename(output_path)
            staging.unlink(missing_ok=True)
        except OSError as exc:
            failure_reason = f"primary result promotion failed: {type(exc).__name__}: {exc}"
            return_code = 125

    if failure_reason is not None:
        receipt = {
            "schema": "local_front_scalar_operator_unknown_receipt.v1",
            "case_id": CASE_ID,
            "status": "UNKNOWN_NO_PRIMARY_RESULT",
            "reason": failure_reason,
            "argv": argv,
            "started_at_utc": started_at,
            "elapsed_seconds": elapsed,
            "wall_timeout_seconds": WALL_TIMEOUT_SECONDS,
            "timed_out": timed_out,
            "return_code": return_code,
            "stdout": stdout,
            "stderr": stderr,
            "staging_path": str(staging),
            "staging_preserved": staging.exists(),
            "staging_sha256": _sha256(staging) if staging.is_file() else None,
            "primary_output_created": output_path.exists(),
            "no_retry": True,
            "finished_at_utc": _utc_now(),
        }
        try:
            _write_unknown_receipt(receipt_path, receipt)
        except FileExistsError:
            return_code = 125

    print(
        json.dumps(
            {
                "case_id": CASE_ID,
                "status": worker_result.get("status") if worker_result is not None and failure_reason is None else "UNKNOWN_NO_PRIMARY_RESULT",
                "primary_output": str(output_path) if failure_reason is None else None,
                "unknown_receipt": str(receipt_path) if failure_reason is not None else None,
                "parent_elapsed_seconds": elapsed,
                "parent_return_code": return_code,
                "timed_out": timed_out,
            },
            ensure_ascii=False,
            allow_nan=False,
        )
    )
    return return_code


if __name__ == "__main__":
    raise SystemExit(main())
