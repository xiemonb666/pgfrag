"""One bounded source-only 162-DOF full 3-D elastic operator diagnostic.

This entry point is prepared for a later, separately authorized single run.
It assembles one complete isotropic K, forms the full-vector exact boundary
load diagnostic, compares the manufactured field and energy, and reports the
162-value spectrum. It never solves a boundary-value problem or changes a
failed matrix. Do not invoke before source review and authorization.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import itertools
import json
import math
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

CASE_ID = "modeiii_pu_full3d_elastic_operator_20261006"
SOURCE_COUNT = 27
DOF_COUNT = 162
VOLUME_COUNT = 256
SIDE_COUNT = 64
CAP_COUNT_EACH = 32
BOUNDARY_COUNT = SIDE_COUNT + 2 * CAP_COUNT_EACH
QUERY_COUNT = VOLUME_COUNT + BOUNDARY_COUNT
MAX_PROVIDER_BATCH = 256
EXPECTED_NEAR_ZERO_MODES = 6
EXPECTED_POSITIVE_MODES = 156
MAX_NUMERIC_ARRAY_BYTES = 8 * 1024 * 1024
MAX_RESULT_JSON_BYTES = 8 * 1024 * 1024
WALL_TIMEOUT_SECONDS = 60.0
STATIC_REPLAY_SCHEMA = "static_solve_replay.v1"
STATIC_REPLAY_ENCODING = "JSON nested arrays; row-major in declared dimensions"
STATIC_REPLAY_QUERY_ORDER = ["volume", "lateral_side", "top_cap", "bottom_cap"]
STATIC_REPLAY_RIGID_MODE_ORDER = ["tx", "ty", "tz", "rx", "ry", "rz"]
STATIC_REPLAY_REGION_SLICES = {
    "volume": [0, VOLUME_COUNT],
    "lateral_side": [VOLUME_COUNT, VOLUME_COUNT + SIDE_COUNT],
    "top_cap": [VOLUME_COUNT + SIDE_COUNT, VOLUME_COUNT + SIDE_COUNT + CAP_COUNT_EACH],
    "bottom_cap": [VOLUME_COUNT + SIDE_COUNT + CAP_COUNT_EACH, QUERY_COUNT],
}
STATIC_REPLAY_ARRAY_UNITS = {
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

# The first five entries are frozen dependencies from the reviewed 2025-10-05
# scalar bridge. The remaining pins freeze this reviewed operator/provider and
# the one fixed case before any numerical module is imported.
EXPECTED_DEPENDENCIES = {
    "enriched_dof_contract.py": "e3ddb3b067465bb5b269576a1b8a5100e2bb92c9c924bc09c949702c19efa9d8",
    "field_reference.py": "7b07dfb30a64c50734dc18de2778b7725320de6d5a8f275c55f301c453ef18d7",
    "mls3d_shape.py": "66e680616a988182875d17aeb84fadb541867aee13fc15e0ae098f7c98839ef5",
    "mls3d_mls.py": "404256d9b7546361a5704c4c3cec7c022069a9cf83a137c58a16e7a64f79ec63",
    "efgm_elastic_operator.py": "41a5ac6ef6359ed1640ee681c94a3496c18652b29cc72c679b3e37c20dd9ac0d",
    "full3d_cracked_elastic_operator.py": "f8853c24e6d3fceabe5835c90ca24cb47cdc186f718e6b705b0c92310f285022",
    "mls_full3d_enriched_provider.py": "53eb2e998e2d239f941e338375af2f7fd8b7a0d90d7d252a6d0b5f359ea093ce",
    "fixed_case.json": "475e63bb059a8464de20f8f530ecf90d3fa3a94163f04342b8fa75b6ab230ed6",
}


class DiagnosticBudgetError(RuntimeError):
    """Raised when fixed array or wall-time budgets are exceeded."""


class DiagnosticNumericsError(FloatingPointError):
    """Raised when stiffness, energy, traction, or spectrum is non-finite."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _timestamp() -> str:
    return datetime.now(timezone.utc).isoformat()


def _dependency_paths() -> dict[str, Path]:
    return {
        "enriched_dof_contract.py": CONTRACT_DIR / "enriched_dof_contract.py",
        "field_reference.py": FIELD_DIR / "field_reference.py",
        "mls3d_shape.py": MLS_DIR / "mls3d_shape.py",
        "mls3d_mls.py": MLS_DIR / "mls3d_mls.py",
        "efgm_elastic_operator.py": OPERATOR_DIR / "efgm_elastic_operator.py",
        "full3d_cracked_elastic_operator.py": HERE / "full3d_cracked_elastic_operator.py",
        "mls_full3d_enriched_provider.py": HERE / "mls_full3d_enriched_provider.py",
        "fixed_case.json": CASE_PATH,
    }


def _verify_dependencies() -> dict[str, str]:
    paths = _dependency_paths()
    observed = {name: _sha256(path) for name, path in paths.items()}
    mismatches = {
        name: {"expected": EXPECTED_DEPENDENCIES[name], "actual": observed[name]}
        for name in paths
        if observed[name] != EXPECTED_DEPENDENCIES[name]
    }
    if mismatches:
        raise RuntimeError(f"read-only source pin mismatch; no numerical imports were attempted: {mismatches}")
    return observed


def _load_numerical_dependencies() -> None:
    """Import numeric modules only after source hashes pass in the timed worker."""
    for path in (str(MLS_DIR), str(OPERATOR_DIR), str(CONTRACT_DIR), str(FIELD_DIR), str(HERE)):
        if path not in sys.path:
            sys.path.insert(0, path)
    global np, BasisKind, DofKey, FieldComponent, DofRegistry
    global MLS3DConfig, ModeIIIFixedCase, IsotropicMaterial
    global ContinuousMLSFull3DFrontProvider
    global assemble_full3d_enriched_operator, build_full3d_registry
    global basis_row_ids, rigid_mode_coefficient_vectors
    import numpy as np
    from enriched_dof_contract import BasisKind, DofKey, DofRegistry, FieldComponent
    from efgm_elastic_operator import IsotropicMaterial
    from field_reference import ModeIIIFixedCase
    from mls3d_mls import MLS3DConfig
    from full3d_cracked_elastic_operator import (
        assemble_full3d_enriched_operator,
        basis_row_ids,
        build_full3d_registry,
        rigid_mode_coefficient_vectors,
    )
    from mls_full3d_enriched_provider import ContinuousMLSFull3DFrontProvider


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
    boundary = case["boundary_conditions"]
    assert all(isinstance(value, dict) for value in (geometry, sources, material, dofs, mls, quadrature, limits, boundary))
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
        geometry["crack_surface"] == "x < 0, y = 0, -L/2 <= z <= L/2",
        sources["grid_shape"] == [3, 3, 3],
        sources["axis_coordinates_m"] == [-0.005, 0.0, 0.005],
        sources["grid_spacing_m"] == 0.005,
        sources["source_count"] == SOURCE_COUNT,
        sources["support_radius_m"] == 0.025,
        sources["all_front_support_source_ids"] == "all 27 sources",
        sources["heaviside_set_J_H"] == [],
        dofs["coefficient_dof_count"] == DOF_COUNT,
        dofs["field_component"] == "UX, UY, and UZ for each source and each basis family",
        dofs["basis_kinds"] == ["BASE", "FRONT_F1_UNSHIFTED"],
        dofs["coefficients_per_source"] == 6,
        material["young_modulus_Pa"] == 1_000_000.0,
        material["poisson_ratio"] == 0.3,
        material["K_III_Pa_sqrt_m"] == 100.0,
        material["expected_exact_total_strain_energy_J"] == 1.3e-6,
        material["exact_stress"] == "sigma_xz=sigma_zx=mu*duz/dx; sigma_yz=sigma_zy=mu*duz/dy; all other stress components are zero",
        material["base_coefficients"] == "all 81 UX/UY/UZ BASE coefficients equal 0 m",
        material["front_coefficient_by_component"] == "UX=UY=0 m; UZ=A*sqrt(R) m for every source",
        material["front_coefficients"] == "all UX and UY FRONT_F1 coefficients equal 0 m; all 27 UZ FRONT_F1 coefficients equal A*sqrt(R) m",
        mls["reference_origin_m"] == [0.0, 0.0, 0.0],
        mls["reference_scale_m"] == 0.01,
        mls["support_radius_m"] == 0.025,
        mls["visibility_mask"] is None,
        volume["point_count"] == VOLUME_COUNT,
        volume["radial_legendre_order"] == 4,
        volume["theta_legendre_order"] == 16,
        volume["z_legendre_order"] == 4,
        side["point_count"] == SIDE_COUNT,
        caps["point_count_each"] == CAP_COUNT_EACH,
        caps["full_vector_traction"] == "top: (mu*duz/dx, mu*duz/dy, 0); bottom: (-mu*duz/dx, -mu*duz/dy, 0)",
        caps["radial_parameterization"] == "r=R*s^2, s in (0,1), using GL4 in s",
        caps["area_weight_formula"] == "dA=r dr dtheta=2*R^2*s^3 ds dtheta",
        quadrature["total_unique_query_count"] == QUERY_COUNT,
        quadrature["maximum_queries_per_provider_batch"] == MAX_PROVIDER_BATCH,
        quadrature["provider_batch_sizes"] == [VOLUME_COUNT, BOUNDARY_COUNT],
        limits["maximum_predeclared_total_query_count"] == QUERY_COUNT,
        limits["maximum_queries_per_provider_batch"] == MAX_PROVIDER_BATCH,
        limits["provider_batch_count"] == 2,
        limits["maximum_source_count"] == SOURCE_COUNT,
        limits["maximum_total_dofs"] == DOF_COUNT,
        limits["maximum_quadrature_points"] == QUERY_COUNT,
        limits["maximum_estimated_numeric_array_bytes"] == MAX_NUMERIC_ARRAY_BYTES,
        limits["maximum_result_json_bytes"] == MAX_RESULT_JSON_BYTES,
        limits["wall_timeout_seconds"] == WALL_TIMEOUT_SECONDS,
        limits["stiffness_dimension"] == [DOF_COUNT, DOF_COUNT],
        limits["assemble_one_full3d_stiffness_matrix"] is True,
        limits["assemble_one_scalar_stiffness_matrix"] is False,
        limits["expected_rigid_modes"] == EXPECTED_NEAR_ZERO_MODES,
        limits["expected_positive_eigenvalues"] == EXPECTED_POSITIVE_MODES,
        limits["fail_on_extra_near_zero_modes"] is True,
        limits["check_six_rigid_virtual_work"] is True,
        limits["report_total_boundary_force_and_moment"] is True,
        limits["no_load_projection_or_equilibrium_repair"] is True,
        limits["run_solver"] is False,
        limits["solve_allowed"] is False,
        limits["apply_load_to_solver"] is False,
        limits["apply_full_vector_boundary_traction_to_solver"] is False,
        limits["use_pseudoinverse_or_regularization"] is False,
        limits["mechanics_field_qualification"] is False,
        boundary["loading"]["top_cap"] == "n=(0,0,1), t=(mu*duz/dx, mu*duz/dy, 0)",
        boundary["loading"]["bottom_cap"] == "n=(0,0,-1), t=(-mu*duz/dx, -mu*duz/dy, 0)",
        boundary["loading"]["crack_faces"].startswith("zero traction on both faces"),
        boundary["analytic_wrench"]["scale"] == "Q=mu*A*R^(3/2)",
        boundary["analytic_wrench"]["lateral_side_moment_about_origin_Nm"][0] == "(4/3)*Q*L",
        boundary["analytic_wrench"]["top_cap_moment_about_origin_Nm"][2] == "(4/5)*Q*R",
        boundary["analytic_wrench"]["bottom_cap_moment_about_origin_Nm"][2] == "-(4/5)*Q*R",
        boundary["analytic_wrench"]["total_force_N"] == [0, 0, 0],
        boundary["analytic_wrench"]["total_moment_about_origin_Nm"] == [0, 0, 0],
        boundary["analytic_wrench"]["interpretation"].endswith("f dot r"),
    )
    if not all(expected):
        raise ValueError("fixed case changed a frozen 3-D geometry, material benchmark, basis, load, or budget input")
    for name, order in (("order_4", 4), ("order_8", 8), ("order_16", 16)):
        rule = rules[name]
        assert isinstance(rule, dict)
        if len(rule["nodes"]) != order or len(rule["weights"]) != order:
            raise ValueError(f"Gauss-Legendre rule {name} has an unexpected order")
        if not all(-1.0 < float(value) < 1.0 for value in rule["nodes"]):
            raise ValueError(f"Gauss-Legendre rule {name} contains an endpoint node")


def _schema_fixture_check() -> None:
    """Exercise only fixed-case and JSON wire-shape validation; no numeric imports."""
    _verify_dependencies()
    case = json.loads(CASE_PATH.read_text(encoding="utf-8"))
    _require_fixed_case(case)
    source_ids = [f"node-{row:03d}" for row in range(SOURCE_COUNT)]
    ordered_keys = []
    for source_id in source_ids:
        for component in ("ux", "uy", "uz"):
            for basis_kind, basis_index in (("base", 0), ("front_f1_unshifted", 1)):
                identity = [source_id, component, basis_kind, basis_index]
                ordered_keys.append(
                    {
                        "slot": len(ordered_keys),
                        "canonical_id": json.dumps(identity, ensure_ascii=True, separators=(",", ":")),
                        "source_id": source_id,
                        "component": component,
                        "basis_kind": basis_kind,
                        "basis_index": basis_index,
                    }
                )
    canonical_ids = [record["canonical_id"] for record in ordered_keys]
    basis_row_ids_saved = [
        {"source_id": source_id, "basis_kind": basis_kind}
        for source_id in source_ids
        for basis_kind in ("base", "front_f1_unshifted")
    ]
    wire_error = _semantic_registry_wire_error(source_ids, canonical_ids, ordered_keys, basis_row_ids_saved)
    if wire_error is not None:
        raise ValueError(f"semantic wire fixture failed: {wire_error}")
    malformed_rows = list(basis_row_ids_saved)
    malformed_rows[0] = {"source_id": source_ids[0], "basis_kind": "BASE"}
    if _semantic_registry_wire_error(source_ids, canonical_ids, ordered_keys, malformed_rows) is None:
        raise ValueError("semantic wire fixture accepted a non-producer enum spelling")
    if _json_shape([[0.0] * 3 for _ in range(QUERY_COUNT)]) != (QUERY_COUNT, 3):
        raise ValueError("nested-array shape fixture failed")

    def zeros(shape: tuple[int, ...]) -> list[object]:
        if len(shape) == 1:
            return [0.0] * shape[0]
        return [zeros(shape[1:]) for _ in range(shape[0])]

    provider_evidence = {
        "volume": {
            "query_count": VOLUME_COUNT,
            "consumed_queries": VOLUME_COUNT,
            "mls_statuses": ["OK"] * VOLUME_COUNT,
            "field_statuses": ["OK"] * VOLUME_COUNT,
            "support_counts": [SOURCE_COUNT] * VOLUME_COUNT,
            "condition_numbers": [1.0] * VOLUME_COUNT,
        },
        "boundary": {
            "query_count": BOUNDARY_COUNT,
            "consumed_queries": BOUNDARY_COUNT,
            "mls_statuses": ["OK"] * BOUNDARY_COUNT,
            "field_statuses": ["OK"] * BOUNDARY_COUNT,
            "support_counts": [SOURCE_COUNT] * BOUNDARY_COUNT,
            "condition_numbers": [1.0] * BOUNDARY_COUNT,
        },
    }
    replay_fixture = {
        "schema": STATIC_REPLAY_SCHEMA,
        "case_id": CASE_ID,
        "upstream_source_sha256": EXPECTED_DEPENDENCIES,
        "array_encoding": STATIC_REPLAY_ENCODING,
        "query_order": STATIC_REPLAY_QUERY_ORDER,
        "rigid_mode_order": STATIC_REPLAY_RIGID_MODE_ORDER,
        "array_units": STATIC_REPLAY_ARRAY_UNITS,
        "operator": {
            "K_N_per_m": zeros((DOF_COUNT, DOF_COUNT)),
            "boundary_load_N": zeros((DOF_COUNT,)),
            "manufactured_coefficients_m": zeros((DOF_COUNT,)),
            "rigid_mode_coefficients": zeros((DOF_COUNT, 6)),
            "rigid_mode_coefficients_by_mode": zeros((6, DOF_COUNT)),
        },
        "registry": {
            "source_ids": source_ids,
            "ordered_canonical_ids": canonical_ids,
            "ordered_keys": ordered_keys,
        },
        "sources": {
            "positions_xyz_m": zeros((SOURCE_COUNT, 3)),
            "rigid_mode_center_xyz_m": zeros((3,)),
            "rigid_mode_rotation_definition": "omega cross (source_position - source_centroid), with omega in radians and no division by R",
        },
        "basis_row_ids": basis_row_ids_saved,
        "quadrature": {
            "positions_xyz_m": zeros((QUERY_COUNT, 3)),
            "volume_weights_m3": zeros((VOLUME_COUNT,)),
            "boundary_normals": zeros((BOUNDARY_COUNT, 3)),
            "boundary_area_weights_m2": zeros((BOUNDARY_COUNT,)),
            "region_slices_half_open": STATIC_REPLAY_REGION_SLICES,
        },
        "basis": {
            "values": zeros((QUERY_COUNT, 54)),
            "gradients_per_m": zeros((QUERY_COUNT, 54, 3)),
        },
        "boundary": {"traction_Pa": zeros((BOUNDARY_COUNT, 3))},
        "exact_field": {
            "displacement_m": zeros((QUERY_COUNT, 3)),
            "displacement_gradient_dimensionless": zeros((QUERY_COUNT, 3, 3)),
            "stress_Pa": zeros((QUERY_COUNT, 3, 3)),
            "energy_density_J_per_m3": zeros((QUERY_COUNT,)),
            "statuses": ["OK"] * QUERY_COUNT,
        },
        "provider_evidence": provider_evidence,
    }
    success_fixture = {
        "status": "OK_DIAGNOSTIC_ONLY_UNQUALIFIED",
        "assembly_status": "OK",
        "mechanical_eligibility": "DIAGNOSTIC_PASS_UNQUALIFIED",
        "solve_allowed": False,
        "source_sha256": EXPECTED_DEPENDENCIES,
        "dimensions": {"semantic_dof_count": DOF_COUNT},
        "spectrum": {
            "symmetric_eigenvalues": [0.0] * EXPECTED_NEAR_ZERO_MODES + [1.0] * EXPECTED_POSITIVE_MODES,
            "zero_tolerance": 1e-12,
            "symmetric_part_eigenvalues_count": DOF_COUNT,
            "expected_near_zero_rigid_modes": EXPECTED_NEAR_ZERO_MODES,
            "observed_near_zero_count": EXPECTED_NEAR_ZERO_MODES,
            "expected_positive_modes": EXPECTED_POSITIVE_MODES,
            "observed_positive_count": EXPECTED_POSITIVE_MODES,
            "significant_negative_count": 0,
            "gate_passed": True,
        },
        "static_solve_replay": replay_fixture,
        "mls": {
            "provider_queries_consumed_total": QUERY_COUNT,
            "status_counts": {"OK": QUERY_COUNT},
        },
    }
    promotable, reason = _promotable_success_result(success_fixture, 0)
    if not promotable:
        raise ValueError(f"complete promotion-gate schema fixture failed: {reason}")
    if _promotable_success_result(success_fixture, 3)[0]:
        raise ValueError("complete promotion gate accepted a nonzero worker return code")


def _array_budget_estimate() -> dict[str, int]:
    f8 = 8
    matrix_bytes = DOF_COUNT * DOF_COUNT * f8
    parts = {
        "source_coordinates": SOURCE_COUNT * 3 * f8,
        "54_duplicated_basis_site_coordinates": 54 * 3 * f8,
        "all_384_predeclared_volume_and_boundary_queries": QUERY_COUNT * 3 * f8,
        "volume_and_boundary_quadrature_weights": QUERY_COUNT * f8,
        "boundary_unit_normals": (SIDE_COUNT + 2 * CAP_COUNT_EACH) * 3 * f8,
        "fixed_gauss_rule_and_map_arrays": 2048 * f8,
        "MLS_weights": QUERY_COUNT * SOURCE_COUNT * f8,
        "MLS_gradients": QUERY_COUNT * SOURCE_COUNT * 3 * f8,
        "BASE_and_F1_scalar_shape_values": QUERY_COUNT * 54 * f8,
        "BASE_and_F1_scalar_shape_gradients": QUERY_COUNT * 54 * 3 * f8,
        "concatenated_two_batch_basis_values_and_gradients": QUERY_COUNT * 54 * 4 * f8,
        "boundary_shape_value_snapshot": BOUNDARY_COUNT * 54 * f8,
        "field_reference_all_outputs_and_statuses": QUERY_COUNT * 72 * f8,
        "reconstructed_and_exact_vector_fields": 2 * QUERY_COUNT * 12 * f8,
        "validated_volume_sample_copies": VOLUME_COUNT * 54 * 5 * f8,
        "full_vector_B_local_products_and_reordered_K": 6 * DOF_COUNT * f8 + 4 * matrix_bytes,
        "stiffness_symmetry_spectrum_work_and_eigensolver_workspace": 16 * matrix_bytes,
        "six_rigid_vectors_coefficients_and_load_vector": (DOF_COUNT * 8 + DOF_COUNT * 6) * f8,
        "boundary_field_traction_load_and_work_rows": (SIDE_COUNT + 2 * CAP_COUNT_EACH) * 24 * f8,
    }
    return {**parts, "estimated_explicit_numeric_array_bytes": int(sum(parts.values()))}


def _gauss_rule(case: dict[str, object], order: int) -> tuple[object, object]:
    rules = case["gauss_legendre_rules"]
    return np.asarray(rules[f"order_{order}"]["nodes"], dtype=np.float64), np.asarray(
        rules[f"order_{order}"]["weights"], dtype=np.float64
    )


def _map_rule(nodes: object, weights: object, lower: float, upper: float) -> tuple[object, object]:
    half = 0.5 * (upper - lower)
    middle = 0.5 * (upper + lower)
    return middle + half * nodes, half * weights


def _build_sources(case: dict[str, object]) -> tuple[object, tuple[str, ...]]:
    axis = np.asarray(case["sources"]["axis_coordinates_m"], dtype=np.float64)
    xyz = [(float(x), float(y), float(z)) for x, y, z in itertools.product(axis, axis, axis)]
    positions = np.ascontiguousarray(np.asarray(xyz, dtype=np.float64))
    return positions, tuple(f"node-{row:03d}" for row in range(len(xyz)))


def _build_quadrature(case: dict[str, object]) -> dict[str, object]:
    radius = float(case["geometry"]["radius_m"])
    length = float(case["geometry"]["length_m"])
    half_length = 0.5 * length
    r4, wr4 = _gauss_rule(case, 4)
    theta8, wt8 = _gauss_rule(case, 8)
    theta16, wt16 = _gauss_rule(case, 16)
    z4, wz4 = _gauss_rule(case, 4)
    radial, radial_weights = _map_rule(r4, wr4, 0.0, radius)
    theta, theta_weights = _map_rule(theta16, wt16, -math.pi, math.pi)
    axial, axial_weights = _map_rule(z4, wz4, -half_length, half_length)

    volume_points: list[tuple[float, float, float]] = []
    volume_weights: list[float] = []
    for ir, it, iz in itertools.product(range(4), range(16), range(4)):
        r, angle = float(radial[ir]), float(theta[it])
        volume_points.append((r * math.cos(angle), r * math.sin(angle), float(axial[iz])))
        volume_weights.append(float(radial_weights[ir] * theta_weights[it] * axial_weights[iz] * r))
    volume_xyz = np.ascontiguousarray(np.asarray(volume_points, dtype=np.float64))
    volume_w = np.ascontiguousarray(np.asarray(volume_weights, dtype=np.float64))

    side_theta, side_theta_weights = _map_rule(theta16, wt16, -math.pi, math.pi)
    side_points: list[tuple[float, float, float]] = []
    side_normals: list[tuple[float, float, float]] = []
    side_weights: list[float] = []
    for it, iz in itertools.product(range(16), range(4)):
        angle = float(side_theta[it])
        cosine, sine = math.cos(angle), math.sin(angle)
        side_points.append((radius * cosine, radius * sine, float(axial[iz])))
        side_normals.append((cosine, sine, 0.0))
        side_weights.append(float(radius * side_theta_weights[it] * axial_weights[iz]))
    side_xyz = np.ascontiguousarray(np.asarray(side_points, dtype=np.float64))
    side_n = np.ascontiguousarray(np.asarray(side_normals, dtype=np.float64))
    side_w = np.ascontiguousarray(np.asarray(side_weights, dtype=np.float64))

    cap_parameter, cap_parameter_weights = _map_rule(r4, wr4, 0.0, 1.0)
    cap_theta, cap_theta_weights = _map_rule(theta8, wt8, -math.pi, math.pi)
    cap_xy: list[tuple[float, float]] = []
    cap_w_values: list[float] = []
    for ir, it in itertools.product(range(4), range(8)):
        s, angle = float(cap_parameter[ir]), float(cap_theta[it])
        r = radius * s * s
        cap_xy.append((r * math.cos(angle), r * math.sin(angle)))
        cap_w_values.append(
            float(2.0 * radius * radius * s**3 * cap_parameter_weights[ir] * cap_theta_weights[it])
        )
    cap_array = np.asarray(cap_xy, dtype=np.float64)
    cap_w = np.ascontiguousarray(np.asarray(cap_w_values, dtype=np.float64))
    top_xyz = np.ascontiguousarray(np.column_stack((cap_array[:, 0], cap_array[:, 1], np.full(CAP_COUNT_EACH, half_length))))
    bottom_xyz = np.ascontiguousarray(np.column_stack((cap_array[:, 0], cap_array[:, 1], np.full(CAP_COUNT_EACH, -half_length))))
    top_n = np.ascontiguousarray(np.tile(np.array([0.0, 0.0, 1.0]), (CAP_COUNT_EACH, 1)))
    bottom_n = np.ascontiguousarray(np.tile(np.array([0.0, 0.0, -1.0]), (CAP_COUNT_EACH, 1)))

    boundary_xyz = np.ascontiguousarray(np.vstack((side_xyz, top_xyz, bottom_xyz)))
    boundary_normals = np.ascontiguousarray(np.vstack((side_n, top_n, bottom_n)))
    boundary_weights = np.ascontiguousarray(np.concatenate((side_w, cap_w, cap_w)))
    query_xyz = np.ascontiguousarray(np.vstack((volume_xyz, boundary_xyz)))
    return {
        "volume_points": volume_xyz,
        "volume_weights": volume_w,
        "side_points": side_xyz,
        "side_normals": side_n,
        "side_weights": side_w,
        "top_points": top_xyz,
        "bottom_points": bottom_xyz,
        "cap_area_weights": cap_w,
        "boundary_points": boundary_xyz,
        "boundary_normals": boundary_normals,
        "boundary_area_weights": boundary_weights,
        "query_points": query_xyz,
        "volume_count": int(volume_xyz.shape[0]),
        "side_count": int(side_xyz.shape[0]),
        "cap_count_each": int(cap_w.shape[0]),
        "query_count": int(query_xyz.shape[0]),
        "r_min_sample_m": float(np.min(radial)),
    }


def _field_case(case: dict[str, object]) -> object:
    geometry = case["geometry"]
    material = case["material_and_manufactured_field"]
    return ModeIIIFixedCase(
        radius_m=float(geometry["radius_m"]),
        length_m=float(geometry["length_m"]),
        K_III_Pa_sqrt_m=float(material["K_III_Pa_sqrt_m"]),
        young_modulus_Pa=float(material["young_modulus_Pa"]),
        poisson_ratio=float(material["poisson_ratio"]),
    )


def _mls_config(case: dict[str, object]) -> object:
    mls = case["mls"]
    return MLS3DConfig(
        reference_origin=tuple(float(value) for value in mls["reference_origin_m"]),
        reference_scale=float(mls["reference_scale_m"]),
        support_radius=float(mls["support_radius_m"]),
        max_condition_number=float(mls["max_condition_number"]),
        rank_relative_tolerance=float(mls["rank_relative_tolerance"]),
        switch_tolerance=float(mls["switch_tolerance"]),
        duplicate_tolerance=float(mls["duplicate_tolerance"]),
    )


def _status_counts(statuses: object) -> dict[str, int]:
    return {
        name.value: sum(status is name for status in statuses)
        for name in sorted(set(statuses), key=lambda status: status.value)
    }


def _text_status_counts(statuses: object) -> dict[str, int]:
    values = [str(status) for status in statuses]
    return {name: values.count(name) for name in sorted(set(values))}


def _json_shape(value: object) -> tuple[int, ...] | None:
    if not isinstance(value, list):
        return ()
    if not value:
        return (0,)
    child_shapes = [_json_shape(item) for item in value]
    if any(shape is None for shape in child_shapes) or any(shape != child_shapes[0] for shape in child_shapes):
        return None
    return (len(value),) + child_shapes[0]


def _semantic_registry_wire_error(
    source_ids: object,
    canonical_ids: object,
    ordered_keys: object,
    basis_row_ids_saved: object,
) -> str | None:
    """Validate producer JSON against the registry's actual enum wire values."""
    if (
        not isinstance(source_ids, list)
        or len(source_ids) != SOURCE_COUNT
        or any(type(source_id) is not str for source_id in source_ids)
        or len(set(source_ids)) != SOURCE_COUNT
    ):
        return "source ID list is not 27 unique entries"
    if not isinstance(canonical_ids, list) or len(canonical_ids) != DOF_COUNT or any(type(item) is not str for item in canonical_ids):
        return "canonical ID list is not 162 entries"
    if not isinstance(ordered_keys, list) or len(ordered_keys) != DOF_COUNT:
        return "ordered key list is not 162 entries"
    if not isinstance(basis_row_ids_saved, list) or len(basis_row_ids_saved) != 54:
        return "basis row ID list is not 54 entries"
    observed_keys: set[tuple[object, object, object, object]] = set()
    for index, record in enumerate(ordered_keys):
        if not isinstance(record, dict):
            return "ordered key record is not an object"
        if record.get("slot") != index or record.get("canonical_id") != canonical_ids[index]:
            return "ordered keys and canonical IDs are not slot aligned"
        identity = [record.get("source_id"), record.get("component"), record.get("basis_kind"), record.get("basis_index")]
        expected_id = json.dumps(identity, ensure_ascii=True, separators=(",", ":"))
        if record.get("canonical_id") != expected_id:
            return "canonical ID does not encode source/component/basis/index"
        source_id, component, basis_kind, basis_index = identity
        if (
            type(source_id) is not str
            or type(component) is not str
            or type(basis_kind) is not str
            or type(basis_index) is not int
            or source_id not in source_ids
            or component not in {"ux", "uy", "uz"}
        ):
            return "ordered key uses an unknown source or component wire value"
        if basis_kind not in {"base", "front_f1_unshifted"}:
            return "ordered key uses an unknown basis enum wire value"
        if basis_index != (0 if basis_kind == "base" else 1):
            return "basis index disagrees with the actual basis enum value"
        observed_keys.add(tuple(identity))
    expected_keys = {
        (source_id, component, basis_kind, 0 if basis_kind == "base" else 1)
        for source_id in source_ids
        for component in ("ux", "uy", "uz")
        for basis_kind in ("base", "front_f1_unshifted")
    }
    if observed_keys != expected_keys or len(set(canonical_ids)) != DOF_COUNT:
        return "ordered semantic keys are not the exact six-per-source registry"
    expected_basis_rows = [
        {"source_id": source_id, "basis_kind": basis_kind}
        for source_id in source_ids
        for basis_kind in ("base", "front_f1_unshifted")
    ]
    if basis_row_ids_saved != expected_basis_rows:
        return "basis row IDs do not use the producer's source-major lowercase enum values"
    return None


def _finite_json_tree(value: object) -> bool:
    """Return false for any non-finite numeric leaf in a JSON-compatible tree."""
    if isinstance(value, float):
        return math.isfinite(value)
    if isinstance(value, dict):
        return all(_finite_json_tree(item) for item in value.values())
    if isinstance(value, (list, tuple)):
        return all(_finite_json_tree(item) for item in value)
    return True


def _finite_numeric_json_tree(value: object) -> bool:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        if isinstance(value, list):
            return all(_finite_numeric_json_tree(item) for item in value)
        return False
    return math.isfinite(float(value))


def _promotable_success_result(result: object, worker_return_code: int) -> tuple[bool, str]:
    """Require every pass gate before a worker result can become the primary receipt."""
    if worker_return_code != 0 or not isinstance(result, dict):
        return False, "worker return code was nonzero or result was not an object"
    if result.get("status") != "OK_DIAGNOSTIC_ONLY_UNQUALIFIED":
        return False, "worker status did not report a diagnostic pass"
    if result.get("assembly_status") != "OK":
        return False, "full 3-D stiffness assembly status was not OK"
    if result.get("mechanical_eligibility") != "DIAGNOSTIC_PASS_UNQUALIFIED":
        return False, "mechanical spectrum eligibility gate did not pass"
    if result.get("solve_allowed") is not False:
        return False, "this diagnostic must never authorize a solver"
    if not _finite_json_tree(result):
        return False, "worker result contains a non-finite numeric value"
    if result.get("source_sha256") != EXPECTED_DEPENDENCIES:
        return False, "worker source hashes do not exactly match the frozen pins"
    dimensions = result.get("dimensions")
    spectrum = result.get("spectrum")
    if not isinstance(dimensions, dict) or not isinstance(spectrum, dict):
        return False, "dimension or spectrum record is missing"
    if dimensions.get("semantic_dof_count") != DOF_COUNT:
        return False, "semantic DOF count differs from the frozen case"
    eigenvalues = spectrum.get("symmetric_eigenvalues")
    if (
        not isinstance(eigenvalues, list)
        or len(eigenvalues) != DOF_COUNT
        or not _finite_numeric_json_tree(eigenvalues)
    ):
        return False, "the complete 162-value spectrum is missing"
    zero_tolerance = spectrum.get("zero_tolerance")
    if isinstance(zero_tolerance, bool) or not isinstance(zero_tolerance, (int, float)):
        return False, "spectrum zero tolerance is missing"
    if not math.isfinite(float(zero_tolerance)) or zero_tolerance < 0.0:
        return False, "spectrum zero tolerance is invalid"
    if spectrum.get("symmetric_part_eigenvalues_count") != DOF_COUNT:
        return False, "spectrum dimension differs from the frozen case"
    if spectrum.get("expected_near_zero_rigid_modes") != EXPECTED_NEAR_ZERO_MODES:
        return False, "the expected six rigid modes were not recorded"
    if spectrum.get("observed_near_zero_count") != EXPECTED_NEAR_ZERO_MODES:
        return False, "the spectrum does not contain exactly six near-zero modes"
    if spectrum.get("expected_positive_modes") != EXPECTED_POSITIVE_MODES:
        return False, "the expected 156 positive modes were not recorded"
    if spectrum.get("observed_positive_count") != EXPECTED_POSITIVE_MODES:
        return False, "the spectrum does not contain 156 positive modes"
    if spectrum.get("significant_negative_count") != 0 or spectrum.get("gate_passed") is not True:
        return False, "significant negative or extra zero modes block mechanics"
    actual_near_zero = sum(abs(value) <= zero_tolerance for value in eigenvalues)
    actual_positive = sum(value > zero_tolerance for value in eigenvalues)
    actual_significant_negative = sum(value < -zero_tolerance for value in eigenvalues)
    if (actual_near_zero, actual_positive, actual_significant_negative) != (
        EXPECTED_NEAR_ZERO_MODES,
        EXPECTED_POSITIVE_MODES,
        0,
    ):
        return False, "the complete spectrum values do not match the claimed rigid/positive counts"
    replay = result.get("static_solve_replay")
    if not isinstance(replay, dict) or replay.get("schema") != STATIC_REPLAY_SCHEMA:
        return False, "static solve replay block is missing or has the wrong schema"
    if replay.get("case_id") != CASE_ID or replay.get("upstream_source_sha256") != EXPECTED_DEPENDENCIES:
        return False, "static solve replay case or upstream pins differ from the frozen source"
    if replay.get("array_encoding") != STATIC_REPLAY_ENCODING:
        return False, "static solve replay encoding is not the frozen JSON row layout"
    if replay.get("query_order") != STATIC_REPLAY_QUERY_ORDER:
        return False, "static solve replay query order differs from the frozen case"
    if replay.get("rigid_mode_order") != STATIC_REPLAY_RIGID_MODE_ORDER:
        return False, "static solve replay rigid-mode order is missing"
    if replay.get("array_units") != STATIC_REPLAY_ARRAY_UNITS:
        return False, "static solve replay unit contract is missing or inconsistent"
    operator = replay.get("operator")
    registry = replay.get("registry")
    sources = replay.get("sources")
    quadrature = replay.get("quadrature")
    basis = replay.get("basis")
    boundary = replay.get("boundary")
    exact = replay.get("exact_field")
    if not all(isinstance(item, dict) for item in (operator, registry, sources, quadrature, basis, boundary, exact)):
        return False, "static solve replay is missing a required data block"
    required_shapes = [
        (operator, "K_N_per_m", (DOF_COUNT, DOF_COUNT)),
        (operator, "boundary_load_N", (DOF_COUNT,)),
        (operator, "manufactured_coefficients_m", (DOF_COUNT,)),
        (operator, "rigid_mode_coefficients", (DOF_COUNT, 6)),
        (operator, "rigid_mode_coefficients_by_mode", (6, DOF_COUNT)),
        (sources, "positions_xyz_m", (SOURCE_COUNT, 3)),
        (sources, "rigid_mode_center_xyz_m", (3,)),
        (quadrature, "positions_xyz_m", (QUERY_COUNT, 3)),
        (quadrature, "volume_weights_m3", (VOLUME_COUNT,)),
        (quadrature, "boundary_normals", (BOUNDARY_COUNT, 3)),
        (quadrature, "boundary_area_weights_m2", (BOUNDARY_COUNT,)),
        (basis, "values", (QUERY_COUNT, 54)),
        (basis, "gradients_per_m", (QUERY_COUNT, 54, 3)),
        (boundary, "traction_Pa", (BOUNDARY_COUNT, 3)),
        (exact, "displacement_m", (QUERY_COUNT, 3)),
        (exact, "displacement_gradient_dimensionless", (QUERY_COUNT, 3, 3)),
    ]
    for block, key, expected_shape in required_shapes:
        if _json_shape(block.get(key)) != expected_shape or not _finite_numeric_json_tree(block.get(key)):
            return False, f"static replay array {key} is missing or has the wrong shape"
    if _json_shape(exact.get("stress_Pa")) != (QUERY_COUNT, 3, 3) or not _finite_numeric_json_tree(exact.get("stress_Pa")):
        return False, "cached exact stress tensor is missing or has the wrong shape"
    if _json_shape(exact.get("energy_density_J_per_m3")) != (QUERY_COUNT,) or not _finite_numeric_json_tree(exact.get("energy_density_J_per_m3")):
        return False, "cached exact energy density is missing or has the wrong shape"
    source_ids = registry.get("source_ids")
    canonical_ids = registry.get("ordered_canonical_ids")
    ordered_keys = registry.get("ordered_keys")
    basis_row_ids_saved = replay.get("basis_row_ids")
    wire_error = _semantic_registry_wire_error(source_ids, canonical_ids, ordered_keys, basis_row_ids_saved)
    if wire_error is not None:
        return False, wire_error
    if not isinstance(exact.get("statuses"), list) or len(exact["statuses"]) != QUERY_COUNT:
        return False, "static replay field status rows are incomplete"
    if any(status != "OK" for status in exact["statuses"]):
        return False, "a cached exact-field query status is not OK"
    provider_evidence = replay.get("provider_evidence")
    if not isinstance(provider_evidence, dict):
        return False, "static replay provider evidence is missing"
    for batch_name, expected_count in (("volume", VOLUME_COUNT), ("boundary", BOUNDARY_COUNT)):
        batch = provider_evidence.get(batch_name)
        if not isinstance(batch, dict) or batch.get("query_count") != expected_count:
            return False, f"static replay {batch_name} provider count is invalid"
        if batch.get("consumed_queries") != expected_count:
            return False, f"static replay {batch_name} consumed count is invalid"
        if batch.get("mls_statuses") != ["OK"] * expected_count:
            return False, f"static replay {batch_name} MLS statuses are not all OK"
        if batch.get("field_statuses") != ["OK"] * expected_count:
            return False, f"static replay {batch_name} field statuses are not all OK"
        for evidence_key in ("support_counts", "condition_numbers"):
            values = batch.get(evidence_key)
            if not isinstance(values, list) or len(values) != expected_count or not _finite_numeric_json_tree(values):
                return False, f"static replay {batch_name} {evidence_key} are incomplete or non-finite"
    if quadrature.get("region_slices_half_open") != STATIC_REPLAY_REGION_SLICES:
        return False, "static replay query region slices differ from the frozen quadrature order"
    if sources.get("rigid_mode_rotation_definition") != "omega cross (source_position - source_centroid), with omega in radians and no division by R":
        return False, "static replay rigid rotation center/normalization is ambiguous"
    mls = result.get("mls")
    if not isinstance(mls, dict) or mls.get("provider_queries_consumed_total") != QUERY_COUNT:
        return False, "the full predeclared provider query set was not consumed"
    if mls.get("status_counts") != {"OK": QUERY_COUNT}:
        return False, "not every predeclared MLS/reference row reported OK"
    return True, "all source, finiteness, assembly, and six-rigid-mode gates passed"


def _blocked_result(
    case: dict[str, object],
    status: str,
    dependency_hashes: dict[str, str],
    budget: dict[str, int],
    started_at: str,
    started: float,
    reason: str,
    provider: object | tuple[object, ...] | None = None,
    assembly_status: str = "NOT_ATTEMPTED",
) -> dict[str, object]:
    mls = None
    if provider is not None:
        providers = provider if type(provider) is tuple else (provider,)
        batch_records = []
        all_statuses = []
        all_support_counts = []
        all_conditions = []
        all_unknown_rows = []
        for batch_index, batch_provider in enumerate(providers):
            finite_conditions = batch_provider.mls_result.condition_numbers[
                np.isfinite(batch_provider.mls_result.condition_numbers)
            ]
            all_statuses.extend(batch_provider.mls_result.statuses)
            all_support_counts.extend(int(value) for value in batch_provider.mls_result.support_counts)
            all_conditions.extend(float(value) for value in finite_conditions)
            all_unknown_rows.extend(
                {"batch_index": batch_index, **row} for row in batch_provider.unknown_rows
            )
            batch_records.append(
                {
                    "batch_index": batch_index,
                    "query_count": batch_provider.query_count,
                    "status_counts": _status_counts(batch_provider.mls_result.statuses),
                    "consumed_queries": batch_provider.consumed_queries,
                }
            )
        mls = {
            "batches": batch_records,
            "total_predeclared_query_count": sum(batch["query_count"] for batch in batch_records),
            "status_counts": _status_counts(tuple(all_statuses)),
            "support_count_min": min(all_support_counts) if all_support_counts else None,
            "support_count_max": max(all_support_counts) if all_support_counts else None,
            "condition_max": max(all_conditions) if all_conditions else None,
            "unknown_rows": all_unknown_rows,
        }
    return {
        "schema": "modeiii_pu_full3d_operator_diagnostic.v1",
        "case_id": CASE_ID,
        "status": status,
        "classification": "DIAGNOSTIC_ONLY_UNQUALIFIED",
        "assembly_status": assembly_status,
        "mechanical_eligibility": "BLOCKED",
        "solve_allowed": False,
        "started_at_utc": started_at,
        "fixed_case": case,
        "blocked_reason": reason,
        "mls": mls,
        "array_budget": {
            "parts_bytes": budget,
            "estimated_explicit_numeric_array_bytes": budget.get("estimated_explicit_numeric_array_bytes", 0),
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
    budget = _array_budget_estimate()
    if budget["estimated_explicit_numeric_array_bytes"] > MAX_NUMERIC_ARRAY_BYTES:
        raise DiagnosticBudgetError("explicit numeric array estimate exceeds the 8 MiB cap")

    # Public counts and source pins pass before any fixture, MLS, or matrix
    # arrays are allocated. All 384 locations are predeclared and split into
    # two provider batches of 256 volume rows and 128 boundary rows.
    _load_numerical_dependencies()
    source_positions, source_ids = _build_sources(case)
    quadrature = _build_quadrature(case)
    query_positions = quadrature["query_points"]
    if (
        source_positions.shape != (SOURCE_COUNT, 3)
        or query_positions.shape != (QUERY_COUNT, 3)
        or quadrature["volume_count"] != VOLUME_COUNT
        or quadrature["side_count"] != SIDE_COUNT
        or quadrature["cap_count_each"] != CAP_COUNT_EACH
    ):
        raise ValueError("fixed source or quadrature builder returned unexpected dimensions")

    volume_points = quadrature["volume_points"]
    boundary_points = quadrature["boundary_points"]
    if volume_points.shape[0] > MAX_PROVIDER_BATCH or boundary_points.shape[0] > MAX_PROVIDER_BATCH:
        raise DiagnosticBudgetError("a frozen provider batch exceeds 256 queries")

    registry = build_full3d_registry(source_ids)
    if registry.dof_count != DOF_COUNT or registry.source_count != SOURCE_COUNT:
        raise RuntimeError("semantic registry does not contain exactly 162 coefficients")
    field_case = _field_case(case)
    volume_provider = ContinuousMLSFull3DFrontProvider(
        source_positions,
        source_ids,
        volume_points,
        registry,
        _mls_config(case),
        field_case,
    )
    boundary_provider = ContinuousMLSFull3DFrontProvider(
        source_positions,
        source_ids,
        boundary_points,
        registry,
        _mls_config(case),
        field_case,
    )
    providers = (volume_provider, boundary_provider)
    if time.monotonic() - started > WALL_TIMEOUT_SECONDS:
        raise DiagnosticBudgetError("60-second wall limit exceeded during MLS/reference evaluation")
    if any(provider.unknown_rows for provider in providers):
        return _blocked_result(
            case,
            "BLOCKED_UNKNOWN_SUPPORT_OR_REFERENCE",
            dependency_hashes,
            budget,
            started_at,
            started,
            "at least one of the 384 predeclared MLS or exact-field rows is UNKNOWN; no stiffness was assembled",
            providers,
        )

    material_data = case["material_and_manufactured_field"]
    material = IsotropicMaterial(
        young_modulus=float(material_data["young_modulus_Pa"]),
        poisson_ratio=float(material_data["poisson_ratio"]),
    )
    volume_weights = quadrature["volume_weights"]
    assembly = assemble_full3d_enriched_operator(
        source_positions,
        source_ids,
        volume_points,
        volume_weights,
        material,
        volume_provider,
        registry,
    )
    if assembly.status.value != "OK" or assembly.stiffness is None:
        return _blocked_result(
            case,
            "BLOCKED_PROVIDER_STATUS",
            dependency_hashes,
            budget,
            started_at,
            started,
            assembly.diagnostic or assembly.status.value,
            providers,
            "BLOCKED",
        )
    stiffness = assembly.stiffness
    if not np.all(np.isfinite(stiffness)):
        raise DiagnosticNumericsError("assembled full 3-D stiffness contains non-finite entries")
    if volume_provider.consumed_queries != VOLUME_COUNT:
        raise RuntimeError("volume assembly did not consume the exact 256-query volume batch")

    coefficients = np.zeros(DOF_COUNT, dtype=np.float64)
    front_value = field_case.amplitude_sqrt_m * math.sqrt(field_case.radius_m)
    for slot, key in enumerate(registry.ordered_keys):
        if key.basis_kind is BasisKind.FRONT_F1_UNSHIFTED and key.field_component is FieldComponent.UZ:
            coefficients[slot] = front_value

    # Complete boundary rows in the frozen sequence, retaining all vector
    # components of sigma*n including cap tx and ty.
    boundary_count = BOUNDARY_COUNT
    boundary_values = np.empty((boundary_count, 54), dtype=np.float64)
    for boundary_row in range(boundary_count):
        sample = boundary_provider.evaluate(boundary_points[boundary_row])
        if sample.status.value != "OK":
            return _blocked_result(
                case,
                "BLOCKED_PROVIDER_STATUS",
                dependency_hashes,
                budget,
                started_at,
                started,
                sample.diagnostic or sample.status.value,
                providers,
                "OK",
            )
        boundary_values[boundary_row] = sample.weights
    if boundary_provider.consumed_queries != boundary_count:
        raise RuntimeError("boundary provider did not consume its exact 128-query batch")

    basis_rows = {row_id: row for row, row_id in enumerate(basis_row_ids(source_ids))}
    if (
        volume_provider.values is None
        or volume_provider.gradients is None
        or boundary_provider.values is None
        or boundary_provider.gradients is None
    ):
        raise RuntimeError("provider batches reported success without cached scalar basis rows")
    field_values = np.vstack((volume_provider.values, boundary_provider.values))
    field_gradients = np.vstack((volume_provider.gradients, boundary_provider.gradients))
    exact_displacement = np.vstack(
        (volume_provider.field_result.displacement_m, boundary_provider.field_result.displacement_m)
    )
    exact_gradient = np.vstack(
        (volume_provider.field_result.displacement_gradient, boundary_provider.field_result.displacement_gradient)
    )
    discrete_displacement = np.zeros((QUERY_COUNT, 3), dtype=np.float64)
    discrete_gradient = np.zeros((QUERY_COUNT, 3, 3), dtype=np.float64)
    for slot, key in enumerate(registry.ordered_keys):
        row = basis_rows[(key.source_id, key.basis_kind)]
        component = (FieldComponent.UX, FieldComponent.UY, FieldComponent.UZ).index(key.field_component)
        discrete_displacement[:, component] += coefficients[slot] * field_values[:, row]
        discrete_gradient[:, component, :] += coefficients[slot] * field_gradients[:, row, :]
    displacement_abs = np.abs(discrete_displacement - exact_displacement)
    gradient_abs = np.abs(discrete_gradient - exact_gradient)
    displacement_scale = max(float(np.max(np.abs(exact_displacement))), np.finfo(np.float64).tiny)
    gradient_scale = max(float(np.max(np.abs(exact_gradient))), np.finfo(np.float64).tiny)

    exact_quadrature_energy = float(
        np.dot(volume_provider.field_result.energy_density_J_per_m3, volume_weights)
    )
    analytic_energy = float(field_case.exact_total_strain_energy_J)
    coefficient_energy = 0.5 * float(coefficients @ stiffness @ coefficients)

    boundary_normals = quadrature["boundary_normals"]
    boundary_weights = quadrature["boundary_area_weights"]
    boundary_field = boundary_provider.field_result
    traction = np.einsum("qij,qj->qi", boundary_field.stress_Pa, boundary_normals)

    load_vector = np.zeros(DOF_COUNT, dtype=np.float64)
    boundary_displacement = discrete_displacement[VOLUME_COUNT:, :]
    for boundary_row in range(boundary_count):
        for basis_row in range(54):
            source_id, basis = boundary_provider.basis_row_ids[basis_row]
            for component in (FieldComponent.UX, FieldComponent.UY, FieldComponent.UZ):
                key = DofKey(source_id, component, basis, 0 if basis is BasisKind.BASE else 1)
                slot = registry.slot_of(key)
                component_index = (FieldComponent.UX, FieldComponent.UY, FieldComponent.UZ).index(component)
                load_vector[slot] += (
                    boundary_values[boundary_row, basis_row]
                    * traction[boundary_row, component_index]
                    * boundary_weights[boundary_row]
                )
    exact_boundary_work = float(
        np.sum(boundary_field.displacement_m * traction * boundary_weights[:, None])
    )
    coefficient_boundary_work = float(coefficients @ load_vector)
    discrete_boundary_work = float(
        np.sum(boundary_displacement * traction * boundary_weights[:, None])
    )

    side_slice = slice(0, SIDE_COUNT)
    top_slice = slice(SIDE_COUNT, SIDE_COUNT + CAP_COUNT_EACH)
    bottom_slice = slice(SIDE_COUNT + CAP_COUNT_EACH, boundary_count)
    def _wrench(region: slice) -> dict[str, list[float]]:
        region_traction = traction[region]
        region_points = boundary_points[region]
        region_weights = boundary_weights[region]
        force = np.sum(region_traction * region_weights[:, None], axis=0)
        moment = np.sum(np.cross(region_points, region_traction) * region_weights[:, None], axis=0)
        return {
            "force_N": [float(value) for value in force],
            "moment_about_global_origin_Nm": [float(value) for value in moment],
        }

    side_wrench = _wrench(side_slice)
    top_wrench = _wrench(top_slice)
    bottom_wrench = _wrench(bottom_slice)
    net_force = np.sum(traction * boundary_weights[:, None], axis=0)
    net_moment = np.sum(
        np.cross(boundary_points, traction) * boundary_weights[:, None], axis=0
    )
    analytic_scale = material.shear_modulus * field_case.amplitude_sqrt_m * field_case.radius_m**1.5
    analytic_side_mx = (4.0 / 3.0) * analytic_scale * field_case.length_m
    analytic_cap_fy = (4.0 / 3.0) * analytic_scale
    analytic_cap_mx = -(2.0 / 3.0) * analytic_scale * field_case.length_m
    analytic_cap_mz = (4.0 / 5.0) * analytic_scale * field_case.radius_m
    cap_traction_records = {}
    for name, region in (("top", top_slice), ("bottom", bottom_slice)):
        cap_t = traction[region]
        cap_w = boundary_weights[region]
        cap_traction_records[name] = {
            "max_abs_traction_xyz_Pa": [float(np.max(np.abs(cap_t[:, c]))) for c in range(3)],
            "integrated_resultant_N": [float(value) for value in np.sum(cap_t * cap_w[:, None], axis=0)],
            "max_abs_tx_and_ty_Pa": [float(np.max(np.abs(cap_t[:, c]))) for c in (0, 1)],
            "max_abs_tz_Pa": float(np.max(np.abs(cap_t[:, 2]))),
        }

    stiffness_norm = float(np.linalg.norm(stiffness, ord="fro"))
    symmetry_relative = float(
        np.linalg.norm(stiffness - stiffness.T, ord="fro")
        / max(stiffness_norm, np.finfo(np.float64).tiny)
    )
    rigid_vectors = rigid_mode_coefficient_vectors(source_positions, source_ids, registry)
    rigid_center = np.mean(source_positions, axis=0)
    rigid_residuals = []
    for mode_index in range(6):
        mode = rigid_vectors[:, mode_index]
        residual = np.linalg.norm(stiffness @ mode)
        denominator = max(stiffness_norm * np.linalg.norm(mode), np.finfo(np.float64).tiny)
        rigid_residuals.append(float(residual / denominator))
    rigid_load_work = rigid_vectors.T @ load_vector
    integrated_wrench = np.concatenate((net_force, net_moment))
    centroid_moment = net_moment - np.cross(rigid_center, net_force)
    rigid_wrench_about_mode_center = np.concatenate((net_force, centroid_moment))
    rigid_work_minus_wrench = rigid_load_work - rigid_wrench_about_mode_center
    rigid_work_minus_analytic_wrench = rigid_load_work.copy()
    integrated_wrench_minus_analytic = integrated_wrench.copy()

    symmetric_part = 0.5 * (stiffness + stiffness.T)
    eigenvalues = np.linalg.eigvalsh(symmetric_part)
    if not np.all(np.isfinite(eigenvalues)):
        raise DiagnosticNumericsError("symmetric stiffness spectrum contains non-finite eigenvalues")
    all_conditions = np.concatenate(
        (volume_provider.mls_result.condition_numbers, boundary_provider.mls_result.condition_numbers)
    )
    finite_condition = all_conditions[np.isfinite(all_conditions)]
    all_mls_statuses = volume_provider.mls_result.statuses + boundary_provider.mls_result.statuses
    all_support_counts = np.concatenate(
        (volume_provider.mls_result.support_counts, boundary_provider.mls_result.support_counts)
    )
    max_condition = float(np.max(finite_condition)) if finite_condition.size else 1.0
    spectral_scale = max(float(np.max(np.abs(eigenvalues))), np.finfo(np.float64).tiny)
    zero_tolerance = float(64.0 * np.finfo(np.float64).eps * max(1.0, max_condition) * spectral_scale)
    near_zero_count = int(np.count_nonzero(np.abs(eigenvalues) <= zero_tolerance))
    positive_count = int(np.count_nonzero(eigenvalues > zero_tolerance))
    significant_negative_count = int(np.count_nonzero(eigenvalues < -zero_tolerance))
    strict_negative_count = int(np.count_nonzero(eigenvalues < 0.0))
    spectrum_gate = (
        near_zero_count == EXPECTED_NEAR_ZERO_MODES
        and positive_count == EXPECTED_POSITIVE_MODES
        and significant_negative_count == 0
    )
    status = "OK_DIAGNOSTIC_ONLY_UNQUALIFIED" if spectrum_gate else "DIAGNOSTIC_FAILED_SPECTRUM"
    mechanical_eligibility = "DIAGNOSTIC_PASS_UNQUALIFIED" if spectrum_gate else "BLOCKED"
    slot_records = [
        {
            "slot": slot,
            "canonical_id": key.canonical_id,
            "source_id": key.source_id,
            "component": key.field_component.value,
            "basis_kind": key.basis_kind.value,
            "basis_index": key.basis_index,
        }
        for slot, key in enumerate(registry.ordered_keys)
    ]
    basis_id_records = [
        {"source_id": source_id, "basis_kind": basis.value}
        for source_id, basis in volume_provider.basis_row_ids
    ]
    saved_numerical_data = {
        "array_encoding": "JSON nested numeric arrays; row-major in listed dimension order",
        "canonical_coefficient_order": "DofRegistry.ordered_keys; canonical_ids and slot_records align by slot",
        "query_order": ["volume", "lateral_side", "top_cap", "bottom_cap"],
        "query_region_slices_half_open": {
            "volume": [0, VOLUME_COUNT],
            "lateral_side": [VOLUME_COUNT, VOLUME_COUNT + SIDE_COUNT],
            "top_cap": [VOLUME_COUNT + SIDE_COUNT, VOLUME_COUNT + SIDE_COUNT + CAP_COUNT_EACH],
            "bottom_cap": [VOLUME_COUNT + SIDE_COUNT + CAP_COUNT_EACH, QUERY_COUNT],
        },
        "array_units": {
            "raw_semantic_stiffness_K_N_per_m": "N/m",
            "boundary_load_vector_f_N": "N",
            "source_positions_m": "m",
            "query_positions_m": "m",
            "volume_weights_m3": "m^3",
            "boundary_area_weights_m2": "m^2",
            "basis_values_384x54": "dimensionless",
            "basis_gradients_per_m_384x54x3": "1/m",
            "boundary_traction_Pa_128x3": "Pa",
            "exact_displacement_m_384x3": "m",
            "exact_displacement_gradient_384x3x3": "dimensionless",
            "exact_stress_Pa_384x3x3": "Pa",
            "exact_energy_density_J_per_m3_384": "J/m^3",
        },
        "raw_semantic_stiffness_K_N_per_m": stiffness.tolist(),
        "boundary_load_vector_f_N": load_vector.tolist(),
        "rigid_mode_names": ["tx", "ty", "tz", "rx", "ry", "rz"],
        "rigid_mode_coefficients_by_mode": rigid_vectors.T.tolist(),
        "rigid_mode_coefficients_by_dof": rigid_vectors.tolist(),
        "registry_canonical_ids": list(registry.canonical_ids),
        "registry_slot_records": slot_records,
        "source_positions_m": source_positions.tolist(),
        "basis_row_ids": basis_id_records,
        "query_positions_m": quadrature["query_points"].tolist(),
        "volume_weights_m3": volume_weights.tolist(),
        "boundary_normals": boundary_normals.tolist(),
        "boundary_area_weights_m2": boundary_weights.tolist(),
        "basis_values_384x54": field_values.tolist(),
        "basis_gradients_per_m_384x54x3": field_gradients.tolist(),
        "boundary_traction_Pa_128x3": traction.tolist(),
        "exact_displacement_m_384x3": exact_displacement.tolist(),
        "exact_displacement_gradient_384x3x3": exact_gradient.tolist(),
        "exact_stress_Pa_384x3x3": np.vstack(
            (volume_provider.field_result.stress_Pa, boundary_field.stress_Pa)
        ).tolist(),
        "exact_energy_density_J_per_m3_384": np.concatenate(
            (
                volume_provider.field_result.energy_density_J_per_m3,
                boundary_field.energy_density_J_per_m3,
            )
        ).tolist(),
        "provider_field_evidence": {
            "volume": {
                "query_count": volume_provider.query_count,
                "consumed_queries": volume_provider.consumed_queries,
                "mls_statuses": [value.value for value in volume_provider.mls_result.statuses],
                "field_statuses": [str(value) for value in volume_provider.field_result.statuses],
                "support_counts": volume_provider.mls_result.support_counts.tolist(),
                "condition_numbers": volume_provider.mls_result.condition_numbers.tolist(),
            },
            "boundary": {
                "query_count": boundary_provider.query_count,
                "consumed_queries": boundary_provider.consumed_queries,
                "mls_statuses": [value.value for value in boundary_provider.mls_result.statuses],
                "field_statuses": [str(value) for value in boundary_field.statuses],
                "support_counts": boundary_provider.mls_result.support_counts.tolist(),
                "condition_numbers": boundary_provider.mls_result.condition_numbers.tolist(),
            },
        },
    }
    static_solve_replay = {
        "schema": STATIC_REPLAY_SCHEMA,
        "case_id": CASE_ID,
        "upstream_source_sha256": dependency_hashes,
        "array_encoding": STATIC_REPLAY_ENCODING,
        "coefficient_order": "registry.ordered_keys; canonical IDs and ordered keys align by slot",
        "rigid_mode_order": STATIC_REPLAY_RIGID_MODE_ORDER,
        "query_order": STATIC_REPLAY_QUERY_ORDER,
        "array_units": STATIC_REPLAY_ARRAY_UNITS,
        "operator": {
            "K_N_per_m": saved_numerical_data["raw_semantic_stiffness_K_N_per_m"],
            "boundary_load_N": saved_numerical_data["boundary_load_vector_f_N"],
            "manufactured_coefficients_m": coefficients.tolist(),
            "rigid_mode_coefficients": saved_numerical_data["rigid_mode_coefficients_by_dof"],
            "rigid_mode_coefficients_by_mode": saved_numerical_data["rigid_mode_coefficients_by_mode"],
        },
        "registry": {
            "source_ids": list(source_ids),
            "ordered_canonical_ids": saved_numerical_data["registry_canonical_ids"],
            "ordered_keys": saved_numerical_data["registry_slot_records"],
        },
        "sources": {
            "positions_xyz_m": saved_numerical_data["source_positions_m"],
            "rigid_mode_center_xyz_m": rigid_center.tolist(),
            "rigid_mode_rotation_definition": "omega cross (source_position - source_centroid), with omega in radians and no division by R",
        },
        "basis_row_ids": saved_numerical_data["basis_row_ids"],
        "quadrature": {
            "positions_xyz_m": saved_numerical_data["query_positions_m"],
            "volume_weights_m3": saved_numerical_data["volume_weights_m3"],
            "boundary_normals": saved_numerical_data["boundary_normals"],
            "boundary_area_weights_m2": saved_numerical_data["boundary_area_weights_m2"],
            "region_slices_half_open": STATIC_REPLAY_REGION_SLICES,
        },
        "basis": {
            "values": saved_numerical_data["basis_values_384x54"],
            "gradients_per_m": saved_numerical_data["basis_gradients_per_m_384x54x3"],
        },
        "exact_field": {
            "displacement_m": saved_numerical_data["exact_displacement_m_384x3"],
            "displacement_gradient_dimensionless": saved_numerical_data["exact_displacement_gradient_384x3x3"],
            "stress_Pa": saved_numerical_data["exact_stress_Pa_384x3x3"],
            "energy_density_J_per_m3": saved_numerical_data["exact_energy_density_J_per_m3_384"],
            "statuses": [
                str(value)
                for value in list(volume_provider.field_result.statuses)
                + list(boundary_field.statuses)
            ],
        },
        "boundary": {
            "traction_Pa": saved_numerical_data["boundary_traction_Pa_128x3"],
        },
        "provider_evidence": saved_numerical_data["provider_field_evidence"],
    }
    result = {
        "schema": "modeiii_pu_full3d_operator_diagnostic.v1",
        "case_id": CASE_ID,
        "status": status,
        "classification": "DIAGNOSTIC_ONLY_UNQUALIFIED",
        "assembly_status": "OK",
        "mechanical_eligibility": mechanical_eligibility,
        "solve_allowed": False,
        "started_at_utc": started_at,
        "fixed_case": case,
        "dimensions": {
            "source_count": SOURCE_COUNT,
            "semantic_dof_count": DOF_COUNT,
            "scalar_basis_rows": 54,
            "volume_quadrature_count": VOLUME_COUNT,
            "side_quadrature_count": SIDE_COUNT,
            "cap_quadrature_count_each": CAP_COUNT_EACH,
            "total_predeclared_query_count": QUERY_COUNT,
            "provider_batch_query_counts": [volume_provider.query_count, boundary_provider.query_count],
            "maximum_provider_batch_query_count": MAX_PROVIDER_BATCH,
            "basis_partition_of_unity_mapping": "scalar BASE/F1 rows lifted to UX/UY/UZ vector basis then reordered by DofKey canonical IDs",
        },
        "constitutive_and_operator": {
            "lambda_Pa": material.lame_lambda,
            "mu_Pa": material.shear_modulus,
            "D_order": ["sigma_xx", "sigma_yy", "sigma_zz", "tau_xy", "tau_yz", "tau_xz"],
            "engineering_strain_order": ["eps_xx", "eps_yy", "eps_zz", "gamma_xy", "gamma_yz", "gamma_xz"],
            "formula": "K = sum_q B.T @ D @ B * dV; all 6 engineering strains and full isotropic lambda/mu coupling",
            "matrix_shape": [DOF_COUNT, DOF_COUNT],
            "symmetry_relative_frobenius": symmetry_relative,
        },
        "field_and_energy": {
            "max_abs_displacement_error_by_xyz_m": [float(np.max(displacement_abs[:, c])) for c in range(3)],
            "max_relative_displacement_error": float(np.max(displacement_abs)) / displacement_scale,
            "max_abs_displacement_gradient_error_by_component_and_axis_dimensionless": [
                [float(np.max(gradient_abs[:, c, axis])) for axis in range(3)] for c in range(3)
            ],
            "max_relative_gradient_error": float(np.max(gradient_abs)) / gradient_scale,
            "analytic_exact_energy_J": analytic_energy,
            "fixed_volume_quadrature_exact_field_energy_J": exact_quadrature_energy,
            "coefficient_quadratic_energy_J": coefficient_energy,
            "one_rule_is_convergence_evidence": False,
            "benchmark_energy_is_ceramic_material_data": False,
        },
        "full_vector_boundary_load": {
            "traction_source": "exact sigma*n on cylinder side and both end caps; zero crack-face traction is analytic and not sampled across the branch cut",
            "cap_traction_xyz_max_abs_Pa": cap_traction_records,
            "boundary_load_vector_dof_count": int(load_vector.size),
            "exact_field_full_boundary_work_J": exact_boundary_work,
            "manufactured_coefficient_boundary_work_J": coefficient_boundary_work,
            "reconstructed_boundary_work_J": discrete_boundary_work,
            "traction_work_identity_abs_error_J": abs(coefficient_boundary_work - discrete_boundary_work),
            "cap_in_plane_tractions_in_load_vector": True,
            "load_applied_to_solver": False,
        },
        "integrated_boundary_wrench_and_rigid_virtual_work": {
            "work_columns_order": ["tx", "ty", "tz", "rx", "ry", "rz"],
            "rigid_work_equals_force_for_unit_translation_and_moment_for_unit_rotation": True,
            "rigid_rotation_center_xyz_m": [float(value) for value in rigid_center],
            "rotation_moment_comparison": "M_about_source_centroid = M_about_global_origin - centroid cross total_force",
            "rigid_vectors_transpose_load_vector": [float(value) for value in rigid_load_work],
            "integrated_force_N": [float(value) for value in net_force],
            "integrated_moment_about_global_origin_Nm": [float(value) for value in net_moment],
            "integrated_wrench_about_rigid_center": [float(value) for value in rigid_wrench_about_mode_center],
            "rigid_work_minus_integrated_wrench": [float(value) for value in rigid_work_minus_wrench],
            "integrated_wrench_minus_analytic_net_zero": [float(value) for value in integrated_wrench_minus_analytic],
            "rigid_work_minus_analytic_net_zero": [float(value) for value in rigid_work_minus_analytic_wrench],
            "all_six_rigid_virtual_works_checked": True,
            "side_integrated_wrench": side_wrench,
            "top_integrated_wrench": top_wrench,
            "bottom_integrated_wrench": bottom_wrench,
            "analytic_wrench_terms": {
                "scale": "Q=mu*A*R^(3/2)",
                "side_force_N": [0.0, 0.0, 0.0],
                "side_moment_about_origin_Nm": [float(analytic_side_mx), 0.0, 0.0],
                "top_cap_force_N": [0.0, float(analytic_cap_fy), 0.0],
                "top_cap_moment_about_origin_Nm": [float(analytic_cap_mx), 0.0, float(analytic_cap_mz)],
                "bottom_cap_force_N": [0.0, float(-analytic_cap_fy), 0.0],
                "bottom_cap_moment_about_origin_Nm": [float(analytic_cap_mx), 0.0, float(-analytic_cap_mz)],
                "net_force_N": [0.0, 0.0, 0.0],
                "net_moment_about_origin_Nm": [0.0, 0.0, 0.0],
            },
            "cap_radial_parameterization": "r=R*s^2, s in (0,1), dA=2*R^2*s^3 ds dtheta; GL4 x GL8",
            "no_load_projection_or_equilibrium_repair": True,
        },
        "rigid_modes": {
            "semantic_mode_columns": ["tx", "ty", "tz", "rx", "ry", "rz"],
            "BASE_coefficients_are_source_coordinate_shapes": True,
            "all_FRONT_F1_entries_are_zero": True,
            "relative_residuals": rigid_residuals,
        },
        "spectrum": {
            "symmetric_part_eigenvalues_count": int(eigenvalues.size),
            "symmetric_eigenvalues": [float(value) for value in eigenvalues],
            "expected_near_zero_rigid_modes": EXPECTED_NEAR_ZERO_MODES,
            "observed_near_zero_count": near_zero_count,
            "expected_positive_modes": EXPECTED_POSITIVE_MODES,
            "observed_positive_count": positive_count,
            "significant_negative_count": significant_negative_count,
            "strict_negative_count_including_roundoff": strict_negative_count,
            "zero_tolerance": zero_tolerance,
            "minimum_eigenvalue": float(np.min(eigenvalues)),
            "maximum_eigenvalue": float(np.max(eigenvalues)),
            "gate_passed": spectrum_gate,
            "extra_zero_modes_fail_without_matrix_repair": True,
            "pseudoinverse_or_ridge_used": False,
        },
        "mls": {
            "batch_count": 2,
            "batch_query_counts": [volume_provider.query_count, boundary_provider.query_count],
            "status_counts": _status_counts(all_mls_statuses),
            "support_count_min": int(np.min(all_support_counts)),
            "support_count_max": int(np.max(all_support_counts)),
            "condition_number_max": max_condition,
            "visibility_mask": None,
            "provider_queries_consumed_total": volume_provider.consumed_queries + boundary_provider.consumed_queries,
        },
        "field_reference": {
            "volume_status_counts": _text_status_counts(volume_provider.field_result.statuses),
            "boundary_status_counts": _text_status_counts(boundary_field.statuses),
            "all_statuses_ok": all(
                status == "OK"
                for status in list(volume_provider.field_result.statuses) + list(boundary_field.statuses)
            ),
            "traction_reconstructed_from_cached_symmetric_stress_and_normal": True,
        },
        "static_solve_replay": static_solve_replay,
        "array_budget": {
            "parts_bytes": budget,
            "estimated_explicit_numeric_array_bytes": budget["estimated_explicit_numeric_array_bytes"],
            "cap_bytes": MAX_NUMERIC_ARRAY_BYTES,
            "maximum_result_json_bytes": MAX_RESULT_JSON_BYTES,
            "python_object_overhead_and_blas_workspace_excluded": True,
        },
        "source_sha256": dependency_hashes,
        "elapsed_seconds": time.monotonic() - started,
        "finished_at_utc": _timestamp(),
        "limitations": case["limitations"],
    }
    if not _finite_json_tree(result):
        raise DiagnosticNumericsError("diagnostic result contains a non-finite numeric value")
    return result


def _worker(output_path: Path) -> int:
    started = time.monotonic()
    started_at = _timestamp()
    try:
        result = _run_diagnostic()
        exit_code = 0 if result["status"] == "OK_DIAGNOSTIC_ONLY_UNQUALIFIED" else 3
    except Exception as exc:
        result = {
            "schema": "modeiii_pu_full3d_operator_diagnostic.v1",
            "case_id": CASE_ID,
            "status": "WORKER_ERROR",
            "classification": "DIAGNOSTIC_ONLY_UNQUALIFIED",
            "assembly_status": "UNKNOWN",
            "mechanical_eligibility": "BLOCKED",
            "solve_allowed": False,
            "started_at_utc": started_at,
            "error_type": type(exc).__name__,
            "error": str(exc),
            "traceback_tail": traceback.format_exc()[-5000:],
            "no_partial_stiffness_returned": True,
            "elapsed_seconds": time.monotonic() - started,
            "finished_at_utc": _timestamp(),
        }
        exit_code = 2
    output_path.parent.mkdir(parents=True, exist_ok=True)
    serialized = json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    output_path.write_text(serialized, encoding="utf-8")
    print(
        json.dumps(
            {"status": result.get("status"), "assembly_status": result.get("assembly_status"), "stage_path": str(output_path)},
            ensure_ascii=False,
        )
    )
    return exit_code


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE.parents[2] / "tmp" / "review_20261006" / "full3d_cracked_elastic_operator_20261006" / "run01" / "result.json")
    parser.add_argument("--_worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--_schema_check", action="store_true", help=argparse.SUPPRESS)
    arguments = parser.parse_args()
    output_path = arguments.output.resolve()
    if arguments._schema_check:
        _schema_fixture_check()
        print("source schema fixture passed; no numerical module was imported")
        return 0
    if arguments._worker:
        return _worker(output_path)

    # The parent grants exactly one worker attempt. Timeout/error is recorded
    # and is never retried. This outer process does not import numerical code.
    staging_path = output_path.with_name(f".{output_path.name}.worker.json")
    if staging_path.exists() or output_path.exists():
        raise FileExistsError("single-attempt output or worker stage already exists; refusing to overwrite it")
    parent_started = time.monotonic()
    parent_started_at = _timestamp()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        completed = subprocess.run(
            [sys.executable, "-B", str(HERE / "full3d_operator_diagnostic.py"), "--_worker", "--output", str(staging_path)],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=WALL_TIMEOUT_SECONDS,
            check=False,
        )
        try:
            stage_result = json.loads(staging_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            stage_result = None
        stage_exists = staging_path.is_file()
        stage_hash = _sha256(staging_path) if stage_exists else None
        stage_bytes = staging_path.stat().st_size if stage_exists else None
        if completed.returncode == 0 and stage_bytes is not None and stage_bytes <= MAX_RESULT_JSON_BYTES:
            promotable, promotion_reason = _promotable_success_result(stage_result, completed.returncode)
        elif completed.returncode == 0:
            promotable, promotion_reason = False, "complete replay JSON exceeded the 8 MiB result cap"
        else:
            promotable, promotion_reason = False, "worker returned nonzero; only its stage diagnostic may be retained"

        if promotable:
            staging_path.replace(output_path)
            print(json.dumps({"status": stage_result["status"], "primary_path": str(output_path)}, ensure_ascii=False))
            return 0

        stage_status = stage_result.get("status") if isinstance(stage_result, dict) else None
        if stage_status == "DIAGNOSTIC_FAILED_SPECTRUM":
            receipt_status = "MECHANICAL_BLOCKED"
            failure_reason = "stiffness assembled, but spectrum has extra zero or significant negative modes"
        elif completed.returncode == 0:
            receipt_status = "PROMOTION_BLOCKED"
            failure_reason = promotion_reason
        elif stage_status:
            receipt_status = "WORKER_DIAGNOSTIC_BLOCKED"
            failure_reason = str(stage_result.get("blocked_reason", stage_result.get("error", stage_status)))
        else:
            receipt_status = "WORKER_FAILED"
            failure_reason = "worker produced no parseable stage diagnostic"
        result = {
            "schema": "modeiii_pu_full3d_operator_failure_receipt.v1",
            "case_id": CASE_ID,
            "status": receipt_status,
            "assembly_status": stage_result.get("assembly_status", "UNKNOWN") if isinstance(stage_result, dict) else "UNKNOWN",
            "mechanical_eligibility": "BLOCKED",
            "solve_allowed": False,
            "worker_return_code": completed.returncode,
            "failure_reason": failure_reason,
            "worker_stage_path": str(staging_path) if stage_exists else None,
            "worker_stage_sha256": stage_hash,
            "worker_stage_bytes": stage_bytes,
            "no_retry": True,
            "parent_started_at_utc": parent_started_at,
            "parent_elapsed_seconds": time.monotonic() - parent_started,
            "finished_at_utc": _timestamp(),
        }
        exit_code = 3 if receipt_status in {"MECHANICAL_BLOCKED", "PROMOTION_BLOCKED"} else max(completed.returncode, 2)
    except subprocess.TimeoutExpired:
        result = {
            "schema": "modeiii_pu_full3d_operator_failure_receipt.v1",
            "case_id": CASE_ID,
            "status": "WALL_TIMEOUT",
            "assembly_status": "UNKNOWN",
            "mechanical_eligibility": "BLOCKED",
            "solve_allowed": False,
            "wall_timeout_seconds": WALL_TIMEOUT_SECONDS,
            "elapsed_seconds": time.monotonic() - parent_started,
            "worker_stage_path": str(staging_path) if staging_path.is_file() else None,
            "worker_stage_sha256": _sha256(staging_path) if staging_path.is_file() else None,
            "no_retry": True,
            "parent_started_at_utc": parent_started_at,
            "finished_at_utc": _timestamp(),
        }
        exit_code = 2
    output_path.write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, allow_nan=False))
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
