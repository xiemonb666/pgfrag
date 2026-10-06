"""Source-only driver for the bounded full-3D tip-aligned quadrature study.

Only ``--source-check`` is intended during source review. Numerical imports,
query generation, MLS/reference evaluation, assembly, and solve occur only in
the future ``--execute`` path after explicit authorization. Execute under one
external 60-second process timeout; this driver never retries a failed case.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import math
import sys
import time
import traceback
from pathlib import Path


HERE = Path(__file__).resolve().parent
PGFRAG = HERE.parents[1]
ROOT = HERE.parents[2]
CASE_PATH = HERE / "study_case.json"
FROZEN_DIR = PGFRAG / "experiments" / "full3d_cracked_elastic_operator_20261006"
FIELD_DIR = PGFRAG / "experiments" / "modeiii_field_reference_20261005"
MLS_DIR = PGFRAG / "experiments" / "mls3d_kernel_20261005"
OPERATOR_DIR = PGFRAG / "experiments" / "efgm_elastic_operator_20261005"
CONTRACT_DIR = PGFRAG / "experiments" / "enriched_dof_contract_20261005"
FIXED_CASE_PATH = FROZEN_DIR / "fixed_case.json"
STUDY_OUTPUT_ROOT = PGFRAG / "tmp" / "review_20261006" / "tip_aligned_quadrature_study_20261006"
DEFAULT_OUTPUT_DIR = STUDY_OUTPUT_ROOT / "run01"

CASE_ID = "tip_aligned_full3d_quadrature_study_20261006"
SOURCE_COUNT = 27
BASIS_COUNT = 54
DOF_COUNT = 162
MAX_QUERY_BATCH = 256
MAX_WALL_SECONDS = 60.0
MAX_ARRAY_BYTES = 32 * 1024 * 1024
MAX_OUTPUT_BYTES = 20 * 1024 * 1024
RADIUS_M = 0.01
LENGTH_M = 0.01
CAP_SPLIT_M = 0.0005
TRIAL_OFFSETS_M = (0.0, 0.0005)
VOLUME_THETA_ORDERS = (32, 64)
CAP_THETA_ORDERS = (8, 16)
RIGID_ORDER = ("tx", "ty", "tz", "rx", "ry", "rz")

# Dependencies shared with the frozen operator. They are read-only runtime
# inputs; the study's own source and case hashes are recorded in its result.
EXPECTED_INPUTS = {
    "enriched_dof_contract.py": "e3ddb3b067465bb5b269576a1b8a5100e2bb92c9c924bc09c949702c19efa9d8",
    "field_reference.py": "7b07dfb30a64c50734dc18de2778b7725320de6d5a8f275c55f301c453ef18d7",
    "mls3d_shape.py": "66e680616a988182875d17aeb84fadb541867aee13fc15e0ae098f7c98839ef5",
    "mls3d_mls.py": "404256d9b7546361a5704c4c3cec7c022069a9cf83a137c58a16e7a64f79ec63",
    "efgm_elastic_operator.py": "41a5ac6ef6359ed1640ee681c94a3496c18652b29cc72c679b3e37c20dd9ac0d",
    "full3d_cracked_elastic_operator.py": "f8853c24e6d3fceabe5835c90ca24cb47cdc186f718e6b705b0c92310f285022",
    "mls_full3d_enriched_provider.py": "53eb2e998e2d239f941e338375af2f7fd8b7a0d90d7d252a6d0b5f359ea093ce",
    "fixed_case.json": "475e63bb059a8464de20f8f530ecf90d3fa3a94163f04342b8fa75b6ab230ed6",
    "prior_result.json": "3b31e7d40c1d8e4b5d829f0b093fee4edadc966bde4086a4b3f4ba48fdf89cd2",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _input_paths() -> dict[str, Path]:
    return {
        "enriched_dof_contract.py": CONTRACT_DIR / "enriched_dof_contract.py",
        "field_reference.py": FIELD_DIR / "field_reference.py",
        "mls3d_shape.py": MLS_DIR / "mls3d_shape.py",
        "mls3d_mls.py": MLS_DIR / "mls3d_mls.py",
        "efgm_elastic_operator.py": OPERATOR_DIR / "efgm_elastic_operator.py",
        "full3d_cracked_elastic_operator.py": FROZEN_DIR / "full3d_cracked_elastic_operator.py",
        "mls_full3d_enriched_provider.py": FROZEN_DIR / "mls_full3d_enriched_provider.py",
        "fixed_case.json": FIXED_CASE_PATH,
        "prior_result.json": ROOT / "tmp" / "review_20261006" / "full3d_cracked_elastic_operator_20261006" / "run01" / "result.json",
    }


def _source_check() -> dict[str, object]:
    """Pure standard-library contract fixture; intentionally builds no points."""
    case = json.loads(CASE_PATH.read_text(encoding="utf-8"))
    if case.get("case_id") != CASE_ID or case.get("status") != "SOURCE_ONLY_AWAITING_ROOT_REVIEW":
        raise ValueError("study case identity/status mismatch")
    if case.get("execution_authorization") != {
        "authorized": False, "requires_explicit_root_authorization": True, "authorized_at_utc": None
    }:
        raise ValueError("source-check fixture expects execution to remain unauthorized")
    volume = 2 * 4 * 4 * (32 + 64)
    side = 64 * 4
    cap_both_caps_per_rule = {
        order: 2 * 4 * 2 * order * 2  # radial panels x GL4 x half-circles x theta x caps
        for order in CAP_THETA_ORDERS
    }
    caps_both_a = 2 * sum(cap_both_caps_per_rule.values())
    total_rows = volume + side + caps_both_a + 128
    raw_batches = 2 * (2 + 4) + 1 + 2 * (1 + 2) + 1
    assert volume == 3072
    assert side == 256
    assert cap_both_caps_per_rule == {8: 256, 16: 512}
    assert caps_both_a == 1536
    assert total_rows == 4992 and raw_batches == 20
    assert case["query_budget"]["requested_rows"]["total_maximum_rows"] == total_rows
    assert case["query_budget"]["maximum_provider_batches_before_global_exact_dedup"] == raw_batches
    assert case["volume_rule"]["volume_points_per_case"] == {"theta32": 512, "theta64": 1024}
    assert case["boundary_rules"]["side"]["point_count"] == side
    caps = case["boundary_rules"]["end_caps"]
    assert caps["fixed_radial_split_m"] == CAP_SPLIT_M
    assert caps["fixed_split_applies_to_both_trial_tips"] is True
    assert caps["area_jacobian"] == "dA=r*(dr/ds)*ds*dtheta on each radial segment"
    assert caps["points_per_cap_by_theta_order"] == {"8_per_semicircle": 128, "16_per_semicircle": 256}
    assert caps["points_both_caps_by_theta_order"] == {"8_per_semicircle": 256, "16_per_semicircle": 512}
    assert case["fixed_physical_gauge"]["g_rx_exact_m"] == 0.000005031114771170333
    assert case["requested_diagnostics"]["deferred"].startswith("No J contour")
    assert case["numerical_trust_gates"] == {
        "weak_equilibrium_residual_relative_max": 1e-9,
        "gauge_residual_relative_max": 1e-9,
        "saddle_backward_residual_relative_max": 1e-9,
        "energy_identity_relative_max": 1e-9,
        "ktt_energy_relative_defect_max": 1e-9,
        "interpretation": "fixed-front algebraic consistency gates, not physical accuracy or fracture acceptance criteria",
    }
    assert case["mechanical_spectrum_gates"] == {
        "dof_count": DOF_COUNT,
        "zero_negative_eigenvalue_eps_factor": 64,
        "rigid_mode_relative_residual_eps_factor": 128,
        "scaled_asymmetry_relative_max": 1e-12,
        "scaled_asymmetry_relative_formula": "||Kscaled-Kscaled^T||_F/max(||Kscaled||_F,tiny)",
        "zero_negative_eigenvalue_threshold_formula": "64*eps_float64*DOF_COUNT*spectral_scale",
        "rigid_mode_relative_residual_threshold_formula": "128*eps_float64*DOF_COUNT",
        "mls_condition_number_policy": "record and enforce only through the provider's own qualification; never relax mechanical tolerances",
        "interpretation": "engineering numerical gates only; the scaled eigenspectrum check is not a strict eigenvalue proof",
    }

    # Ensure numeric code is not imported at module load and no future
    # execution function has accidentally been invoked by this fixture.
    if "numpy" in sys.modules:
        raise RuntimeError("source-check must run before numerical imports")
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    top_numeric_imports = []
    for node in tree.body:
        if isinstance(node, ast.Import):
            top_numeric_imports.extend(alias.name for alias in node.names if alias.name == "numpy")
        elif isinstance(node, ast.ImportFrom) and node.module == "numpy":
            top_numeric_imports.append(node.module)
    if top_numeric_imports:
        raise RuntimeError("numpy must remain a lazy import inside future execution")
    provider_class = next(node for node in tree.body if isinstance(node, ast.ClassDef) and node.name == "_RowsProvider")
    provider_init = next(node for node in provider_class.body if isinstance(node, ast.FunctionDef) and node.name == "__init__")
    init_args = [arg.arg for arg in provider_init.args.args]
    if "basis_ids" not in init_args:
        raise RuntimeError("cached-row adapter must accept the frozen scalar basis row IDs")
    has_basis_row_assignment = any(
        isinstance(node, ast.Assign)
        and any(isinstance(target, ast.Attribute) and target.attr == "basis_row_ids" for target in node.targets)
        and isinstance(node.value, ast.Call)
        and isinstance(node.value.func, ast.Name) and node.value.func.id == "tuple"
        and len(node.value.args) == 1
        and isinstance(node.value.args[0], ast.Name) and node.value.args[0].id == "basis_ids"
        for node in ast.walk(provider_init)
    )
    execute_fn = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "_execute")
    adapter_calls = [node for node in ast.walk(execute_fn)
                     if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                     and node.func.id == "_RowsProvider"]
    if not has_basis_row_assignment or not any(
        call.args and isinstance(call.args[-1], ast.Name) and call.args[-1].id == "basis_ids"
        for call in adapter_calls
    ):
        raise RuntimeError("cached-row adapter must expose the exact basis row tuple passed to the frozen assembler")
    frozen_tree = ast.parse((FROZEN_DIR / "full3d_cracked_elastic_operator.py").read_text(encoding="utf-8"))
    assembler_fn = next(node for node in frozen_tree.body
                        if isinstance(node, ast.FunctionDef) and node.name == "assemble_full3d_enriched_operator")
    assembler_calls = [node for node in ast.walk(assembler_fn)
                       if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                       and node.func.id == "basis_row_ids"]
    assembler_contract = any(
        isinstance(node, ast.Compare) and isinstance(node.left, ast.Call)
        and isinstance(node.left.func, ast.Name) and node.left.func.id == "getattr"
        and len(node.left.args) >= 2
        and isinstance(node.left.args[0], ast.Name) and node.left.args[0].id == "shape_provider"
        and isinstance(node.left.args[1], ast.Constant) and node.left.args[1].value == "basis_row_ids"
        for node in ast.walk(assembler_fn)
    )
    if not assembler_calls or not assembler_contract:
        raise RuntimeError("frozen assembler no longer exposes the expected basis_row_ids provider contract")
    spectrum_fn = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "_spectrum")
    spectrum_names = {node.id for node in ast.walk(spectrum_fn) if isinstance(node, ast.Name)}
    if "condition_max" in spectrum_names or not {"DOF_COUNT", "zero_factor", "rigid_factor", "spectral_scale"}.issubset(spectrum_names):
        raise RuntimeError("mechanical spectrum thresholds must use only fixed DOF-scaled machine tolerances")
    spectrum_assignments = {}
    for node in ast.walk(spectrum_fn):
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name):
                    spectrum_assignments[target.id] = ast.dump(node.value, include_attributes=False)
    expected_thresholds = {
        "tolerance": "zero_factor * np.finfo(np.float64).eps * DOF_COUNT * spectral_scale",
        "rigid_tolerance": "rigid_factor * np.finfo(np.float64).eps * DOF_COUNT",
    }
    for target, expression in expected_thresholds.items():
        expected_tree = ast.dump(ast.parse(expression, mode="eval").body, include_attributes=False)
        if spectrum_assignments.get(target) != expected_tree:
            raise RuntimeError(f"mechanical spectrum threshold formula changed: {target}")
    finite_arguments = {
        node.args[0].id for node in ast.walk(spectrum_fn)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
        and node.func.attr == "isfinite" and node.args and isinstance(node.args[0], ast.Name)
    }
    has_asymmetry_gate = any(
        isinstance(node, ast.Compare) and isinstance(node.left, ast.Name)
        and node.left.id == "scaled_asymmetry" and node.ops and isinstance(node.ops[0], ast.Gt)
        and node.comparators and isinstance(node.comparators[0], ast.Name)
        and node.comparators[0].id == "asymmetry_limit"
        for node in ast.walk(spectrum_fn)
    )
    if not {"scaled_raw", "eigenvalues"}.issubset(finite_arguments) or not has_asymmetry_gate:
        raise RuntimeError("mechanical spectrum finite-value or asymmetry gate is incomplete")
    return {
        "status": "SOURCE_CHECK_OK_NO_NUMERICAL_IMPORTS",
        "requested_rows_max": total_rows,
        "raw_batch_upper_bound": raw_batches,
        "fixed_shared_cap_split_m": CAP_SPLIT_M,
        "cached_row_basis_contract": "AST_VERIFIED",
        "mechanical_spectrum_contract": "AST_VERIFIED",
        "future_execution": "NOT_RUN",
    }


def _verify_inputs() -> dict[str, str]:
    observed = {name: _sha256(path) for name, path in _input_paths().items()}
    mismatch = {
        name: {"expected": EXPECTED_INPUTS[name], "actual": observed[name]}
        for name in observed
        if observed[name] != EXPECTED_INPUTS[name]
    }
    if mismatch:
        raise RuntimeError(f"read-only frozen input hash mismatch: {mismatch}")
    return observed


def _load_numeric_dependencies() -> dict[str, object]:
    """Import frozen APIs only after all source/result pins passed."""
    for path in (FIELD_DIR, MLS_DIR, OPERATOR_DIR, CONTRACT_DIR, FROZEN_DIR):
        text_path = str(path)
        if text_path not in sys.path:
            sys.path.insert(0, text_path)
    import numpy as np
    from enriched_dof_contract import BasisKind, DofKey, FieldComponent
    from efgm_elastic_operator import IsotropicMaterial, ShapeSample, ShapeStatus
    from field_reference import ModeIIIFixedCase
    from mls3d_mls import MLS3DConfig
    from full3d_cracked_elastic_operator import (
        assemble_full3d_enriched_operator,
        build_full3d_registry,
        rigid_mode_coefficient_vectors,
    )
    from mls_full3d_enriched_provider import ContinuousMLSFull3DFrontProvider

    return locals()


def _legendre(np: object, fixed_case: dict[str, object], order: int) -> tuple[object, object]:
    if order <= 16:
        rule = fixed_case["gauss_legendre_rules"][f"order_{order}"]
        return (
            np.asarray(rule["nodes"], dtype=np.float64),
            np.asarray(rule["weights"], dtype=np.float64),
        )
    nodes, weights = np.polynomial.legendre.leggauss(order)
    return np.asarray(nodes, dtype=np.float64), np.asarray(weights, dtype=np.float64)


def _map_rule(nodes: object, weights: object, lo: float, hi: float) -> tuple[object, object]:
    half = 0.5 * (hi - lo)
    mid = 0.5 * (hi + lo)
    return mid + half * nodes, half * weights


def _build_volume(np: object, fixed: dict[str, object], a: float, theta_order: int) -> dict[str, object]:
    radius, length = RADIUS_M, LENGTH_M
    sr, wr = _legendre(np, fixed, 4)
    sz, wz = _legendre(np, fixed, 4)
    st, wt = _legendre(np, fixed, theta_order)
    s, ws = _map_rule(sr, wr, 0.0, 1.0)
    z, wz = _map_rule(sz, wz, -0.5 * length, 0.5 * length)
    theta, wt = _map_rule(st, wt, -math.pi, math.pi)
    points, weights = [], []
    for it in range(theta_order):
        angle = float(theta[it])
        rho_max = -a * math.cos(angle) + math.sqrt(radius * radius - a * a * math.sin(angle) ** 2)
        for iz in range(4):
            for radial_i in range(4):
                rho = rho_max * float(s[radial_i]) ** 2
                points.append((a + rho * math.cos(angle), rho * math.sin(angle), float(z[iz])))
                weights.append(
                    2.0 * rho_max**2 * float(s[radial_i]) ** 3
                    * float(ws[radial_i]) * float(wt[it]) * float(wz[iz])
                )
    return {"xyz": np.ascontiguousarray(points, dtype=np.float64), "weights": np.asarray(weights, dtype=np.float64)}


def _build_side(np: object, fixed: dict[str, object]) -> dict[str, object]:
    st, wt = _legendre(np, fixed, 64)
    sz, wz = _legendre(np, fixed, 4)
    theta, wt = _map_rule(st, wt, -math.pi, math.pi)
    z, wz = _map_rule(sz, wz, -0.5 * LENGTH_M, 0.5 * LENGTH_M)
    xyz, normals, weights = [], [], []
    max_roundoff_correction = 0.0
    for it in range(64):
        angle = float(theta[it])
        nx, ny = math.cos(angle), math.sin(angle)
        for iz in range(4):
            x, y = RADIUS_M * nx, RADIUS_M * ny
            radial = math.hypot(x, y)
            if radial > RADIUS_M:
                factor = math.nextafter(RADIUS_M / radial, 0.0)
                x, y = factor * x, factor * y
                max_roundoff_correction = max(max_roundoff_correction, radial - math.hypot(x, y))
            xyz.append((x, y, float(z[iz])))
            normals.append((nx, ny, 0.0))
            weights.append(RADIUS_M * float(wt[it]) * float(wz[iz]))
    return {
        "xyz": np.ascontiguousarray(xyz, dtype=np.float64),
        "normals": np.ascontiguousarray(normals, dtype=np.float64),
        "weights": np.asarray(weights, dtype=np.float64),
        "max_inward_roundoff_correction_m": max_roundoff_correction,
    }


def _build_caps(np: object, fixed: dict[str, object], theta_order: int) -> dict[str, object]:
    sr, wr = _legendre(np, fixed, 4)
    st, wt = _legendre(np, fixed, theta_order)
    s, ws = _map_rule(sr, wr, 0.0, 1.0)
    xyz, normals, weights = [], [], []
    radial_maps = (
        ("inner", lambda q: CAP_SPLIT_M * math.sin(0.5 * math.pi * q) ** 2,
         lambda q: 0.5 * math.pi * CAP_SPLIT_M * math.sin(math.pi * q)),
        ("outer", lambda q: CAP_SPLIT_M + (RADIUS_M - CAP_SPLIT_M) * q**2,
         lambda q: 2.0 * (RADIUS_M - CAP_SPLIT_M) * q),
    )
    halves = ((-math.pi, 0.0), (0.0, math.pi))
    for cap_sign in (-1.0, 1.0):
        z = cap_sign * 0.5 * LENGTH_M
        normal = (0.0, 0.0, cap_sign)
        for _segment_name, radius_fn, derivative_fn in radial_maps:
            for lo, hi in halves:
                theta, wtheta = _map_rule(st, wt, lo, hi)
                for it in range(theta_order):
                    angle = float(theta[it])
                    ca, sa = math.cos(angle), math.sin(angle)
                    for ir in range(4):
                        sv = float(s[ir])
                        radius = radius_fn(sv)
                        drds = derivative_fn(sv)
                        xyz.append((radius * ca, radius * sa, z))
                        normals.append(normal)
                        weights.append(radius * drds * float(ws[ir]) * float(wtheta[it]))
    points_per_cap = 16 * theta_order
    return {
        "xyz": np.ascontiguousarray(xyz, dtype=np.float64),
        "normals": np.ascontiguousarray(normals, dtype=np.float64),
        "weights": np.asarray(weights, dtype=np.float64),
        "region_slices": {"bottom": [0, points_per_cap], "top": [points_per_cap, 2 * points_per_cap]},
    }


def _build_offgrid(np: object) -> dict[str, object]:
    radii = (0.03, 0.18, 0.57, 0.89)
    angles = (math.pi / 8.0, -math.pi / 8.0, 3 * math.pi / 8.0, -3 * math.pi / 8.0,
              5 * math.pi / 8.0, -5 * math.pi / 8.0, 7 * math.pi / 8.0, -7 * math.pi / 8.0)
    z_fracs = (-3.0 / 8.0, -1.0 / 8.0, 1.0 / 8.0, 3.0 / 8.0)
    xyz = [
        (RADIUS_M * r * math.cos(t), RADIUS_M * r * math.sin(t), LENGTH_M * zf)
        for r in radii for t in angles for zf in z_fracs
    ]
    return {"xyz": np.ascontiguousarray(xyz, dtype=np.float64)}


def _declared_groups(np: object, fixed: dict[str, object]) -> dict[str, dict[str, object]]:
    groups: dict[str, dict[str, object]] = {}
    shared_cap8 = _build_caps(np, fixed, 8)
    shared_cap16 = _build_caps(np, fixed, 16)
    for a in TRIAL_OFFSETS_M:
        for theta_order in VOLUME_THETA_ORDERS:
            groups[f"a{a:g}.volume{theta_order}"] = _build_volume(np, fixed, a, theta_order)
        # Repeated raw request rows are intentional: one physical cap rule is
        # specialized with each trial basis. Global xyz dedup keeps the old
        # MLS/reference APIs within the declared bound.
        groups[f"a{a:g}.cap8"] = shared_cap8
        groups[f"a{a:g}.cap16"] = shared_cap16
    groups["side"] = _build_side(np, fixed)
    groups["offgrid"] = _build_offgrid(np)
    raw_rows = sum(int(group["xyz"].shape[0]) for group in groups.values())
    if raw_rows != 4992:
        raise RuntimeError(f"predeclared raw point count changed: {raw_rows} != 4992")
    expected = {
        "a0.volume32": 512, "a0.volume64": 1024,
        "a0.0005.volume32": 512, "a0.0005.volume64": 1024,
        "a0.cap8": 256, "a0.cap16": 512,
        "a0.0005.cap8": 256, "a0.0005.cap16": 512,
        "side": 256, "offgrid": 128,
    }
    actual = {name: int(group["xyz"].shape[0]) for name, group in groups.items()}
    if actual != expected:
        raise RuntimeError(f"individual predeclared point counts changed: {actual}")
    return groups


def _exact_unique_map(np: object, groups: dict[str, dict[str, object]]) -> tuple[object, dict[str, object]]:
    # Python tuple keys make dedup exact in float64 and keep stable first-seen
    # ordering. Coordinates are only constructed here, before any query call.
    unique: list[tuple[float, float, float]] = []
    row_by_xyz: dict[tuple[float, float, float], int] = {}
    index_by_group: dict[str, object] = {}
    for name, group in groups.items():
        mapped = []
        for row in group["xyz"]:
            key = (float(row[0]), float(row[1]), float(row[2]))
            index = row_by_xyz.get(key)
            if index is None:
                index = len(unique)
                row_by_xyz[key] = index
                unique.append(key)
            mapped.append(index)
        index_by_group[name] = np.asarray(mapped, dtype=np.int64)
    return np.ascontiguousarray(unique, dtype=np.float64), index_by_group


class _RowsProvider:
    """Small adapter from already cached shifted basis rows to frozen assembler."""

    def __init__(self, np: object, ShapeSample: object, ShapeStatus: object,
                 values: object, gradients: object, expected_xyz: object,
                 basis_ids: tuple[tuple[str, object], ...]) -> None:
        self.np, self.ShapeSample, self.ShapeStatus = np, ShapeSample, ShapeStatus
        self.values, self.gradients = values, gradients
        self.expected_xyz = expected_xyz
        self.basis_row_ids = tuple(basis_ids)
        self.support = np.arange(BASIS_COUNT, dtype=np.int64)
        self.cursor = 0

    def evaluate(self, material_position: object) -> object:
        if self.cursor >= self.values.shape[0]:
            raise RuntimeError("assembler consumed more rows than its declared chunk")
        row = self.cursor
        if not self.np.array_equal(self.np.asarray(material_position, dtype=self.np.float64), self.expected_xyz[row]):
            raise RuntimeError(f"assembler query order diverged at local chunk row {row}")
        self.cursor += 1
        return self.ShapeSample(
            weights=self.values[row], gradients=self.gradients[row],
            support_ids=self.support, status=self.ShapeStatus.OK,
        )


def _shifted_basis(np: object, base_values: object, base_gradients: object, xyz: object, a: float) -> tuple[object, object]:
    dx = xyz[:, 0] - a
    y = xyz[:, 1]
    rho = np.hypot(dx, y)
    if not np.isfinite(rho).all() or np.any(rho <= 0.0):
        raise ValueError("shifted FRONT_F1 query lies at or numerically beyond the declared tip singularity")
    negative_cut = (y == 0.0) & (dx < 0.0)
    if np.any(negative_cut):
        raise ValueError("shifted FRONT_F1 query lies on the negative-x branch cut with an undeclared trace")
    theta = np.arctan2(y, dx)
    root = np.sqrt(rho)
    fbar = root * np.sin(0.5 * theta) / math.sqrt(RADIUS_M)
    grad_fbar = np.empty((xyz.shape[0], 3), dtype=np.float64)
    grad_fbar[:, 0] = -np.sin(0.5 * theta) / (2.0 * root * math.sqrt(RADIUS_M))
    grad_fbar[:, 1] = np.cos(0.5 * theta) / (2.0 * root * math.sqrt(RADIUS_M))
    grad_fbar[:, 2] = 0.0
    values = np.empty((xyz.shape[0], BASIS_COUNT), dtype=np.float64)
    gradients = np.empty((xyz.shape[0], BASIS_COUNT, 3), dtype=np.float64)
    values[:, 0::2] = base_values
    gradients[:, 0::2, :] = base_gradients
    values[:, 1::2] = base_values * fbar[:, None]
    gradients[:, 1::2, :] = (
        base_gradients * fbar[:, None, None]
        + base_values[:, :, None] * grad_fbar[:, None, :]
    )
    if not np.isfinite(values).all() or not np.isfinite(gradients).all():
        raise FloatingPointError("shifted FRONT_F1 values/gradients are non-finite")
    return values, gradients


def _strain(np: object, gradient: object) -> object:
    result = np.empty((gradient.shape[0], 6), dtype=np.float64)
    result[:, 0] = gradient[:, 0, 0]
    result[:, 1] = gradient[:, 1, 1]
    result[:, 2] = gradient[:, 2, 2]
    result[:, 3] = gradient[:, 0, 1] + gradient[:, 1, 0]
    result[:, 4] = gradient[:, 1, 2] + gradient[:, 2, 1]
    result[:, 5] = gradient[:, 0, 2] + gradient[:, 2, 0]
    return result


def _stress_tensor(np: object, stress_vector: object) -> object:
    result = np.zeros((stress_vector.shape[0], 3, 3), dtype=np.float64)
    result[:, 0, 0], result[:, 1, 1], result[:, 2, 2] = stress_vector[:, 0], stress_vector[:, 1], stress_vector[:, 2]
    result[:, 0, 1] = result[:, 1, 0] = stress_vector[:, 3]
    result[:, 1, 2] = result[:, 2, 1] = stress_vector[:, 4]
    result[:, 0, 2] = result[:, 2, 0] = stress_vector[:, 5]
    return result


def _reconstruct(np: object, keys: object, basis_rows: dict[tuple[str, object], int], values: object,
                 gradients: object, coefficients: object, component_index: dict[object, int]) -> tuple[object, object]:
    u = np.zeros((values.shape[0], 3), dtype=np.float64)
    grad = np.zeros((values.shape[0], 3, 3), dtype=np.float64)
    for slot, key in enumerate(keys):
        row = basis_rows[(key.source_id, key.basis_kind)]
        comp = component_index[key.field_component]
        u[:, comp] += coefficients[slot] * values[:, row]
        grad[:, comp, :] += coefficients[slot] * gradients[:, row, :]
    return u, grad


def _weighted_relative(np: object, delta: object, reference: object, weights: object) -> float:
    axes = tuple(range(1, delta.ndim))
    numerator = float(np.dot(weights, np.sum(delta * delta, axis=axes)))
    denominator = float(np.dot(weights, np.sum(reference * reference, axis=axes)))
    return math.sqrt(numerator / max(denominator, np.finfo(np.float64).tiny))


def _wrench(np: object, xyz: object, traction: object, weights: object) -> dict[str, object]:
    force = np.sum(weights[:, None] * traction, axis=0)
    moment = np.sum(weights[:, None] * np.cross(xyz, traction), axis=0)
    return {"force_N": force, "moment_about_origin_Nm": moment}


def _jsonable(np: object, value: object) -> object:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(key): _jsonable(np, item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(np, item) for item in value]
    return value


def _volume_moments(np: object, values: object, weights: object, xyz: object,
                    keys: object, basis_rows: dict[tuple[str, object], int], registry: object,
                    source_positions: object, source_ids: object, rigid_raw: object,
                    component_index: dict[object, int], DofKey: object, BasisKind: object,
                    components: object) -> tuple[object, object, object, object]:
    measure = float(np.sum(weights))
    centroid = np.einsum("q,qc->c", weights, xyz) / measure
    modes = np.zeros((6, xyz.shape[0], 3), dtype=np.float64)
    modes[0, :, 0] = modes[1, :, 1] = modes[2, :, 2] = 1.0
    axes = np.eye(3, dtype=np.float64)
    centered = xyz - centroid
    for axis in range(3):
        modes[3 + axis] = np.cross(axes[axis], centered) / RADIUS_M
    constraint = np.zeros((6, DOF_COUNT), dtype=np.float64)
    for slot, key in enumerate(keys):
        row = basis_rows[(key.source_id, key.basis_kind)]
        comp = component_index[key.field_component]
        constraint[:, slot] = np.einsum("q,lq->l", weights * values[:, row], modes[:, :, comp]) / measure
    source_center = np.mean(source_positions, axis=0)
    center_shift = source_center - centroid
    physical_rigid = np.array(rigid_raw, dtype=np.float64, copy=True)
    for axis in range(3):
        physical_rigid[:, 3 + axis] = (
            rigid_raw[:, 3 + axis] + rigid_raw[:, :3] @ np.cross(axes[axis], center_shift)
        ) / RADIUS_M
    if not np.isfinite(constraint).all() or not np.isfinite(physical_rigid).all():
        raise FloatingPointError("gauge moments or rigid vectors are non-finite")
    return constraint, centroid, physical_rigid, measure


def _spectrum(np: object, stiffness: object, rigid_raw: object,
              spectrum_gates: dict[str, object]) -> dict[str, object]:
    """Apply fixed-dimension engineering gates, not a strict eigenvalue proof."""
    if int(spectrum_gates["dof_count"]) != DOF_COUNT:
        raise ValueError("mechanical spectrum gate dimension differs from the frozen 162-DOF operator")
    zero_factor = float(spectrum_gates["zero_negative_eigenvalue_eps_factor"])
    rigid_factor = float(spectrum_gates["rigid_mode_relative_residual_eps_factor"])
    asymmetry_limit = float(spectrum_gates["scaled_asymmetry_relative_max"])
    diagonal = np.diag(stiffness)
    if not np.isfinite(diagonal).all() or np.any(diagonal <= 0.0):
        raise RuntimeError("raw stiffness has a non-positive diagonal; spectrum gate blocked")
    scale = 1.0 / np.sqrt(diagonal)
    if not np.isfinite(scale).all():
        raise FloatingPointError("mechanical spectrum scaling is non-finite")
    scaled_raw = scale[:, None] * stiffness * scale[None, :]
    if not np.isfinite(scaled_raw).all():
        raise FloatingPointError("scaled raw stiffness contains non-finite entries")
    scaled_asymmetry = float(np.linalg.norm(scaled_raw - scaled_raw.T) /
        max(np.linalg.norm(scaled_raw), np.finfo(np.float64).tiny))
    if not math.isfinite(scaled_asymmetry) or scaled_asymmetry > asymmetry_limit:
        raise RuntimeError(
            f"mechanical BLOCKED: scaled stiffness asymmetry {scaled_asymmetry} exceeds {asymmetry_limit}"
        )
    symmetric = 0.5 * scaled_raw + 0.5 * scaled_raw.T
    if not np.isfinite(symmetric).all():
        raise FloatingPointError("symmetrized scaled stiffness contains non-finite entries")
    eigenvalues = np.linalg.eigvalsh(symmetric)
    if not np.isfinite(eigenvalues).all():
        raise FloatingPointError("scaled stiffness eigenvalues contain non-finite entries")
    spectral_scale = max(float(np.max(np.abs(eigenvalues))), np.finfo(np.float64).tiny)
    tolerance = zero_factor * np.finfo(np.float64).eps * DOF_COUNT * spectral_scale
    near = np.abs(eigenvalues) <= tolerance
    positive = eigenvalues > tolerance
    negative = eigenvalues < -tolerance
    rigid = rigid_raw / scale[:, None]
    residuals = np.linalg.norm(symmetric @ rigid, axis=0) / np.maximum(
        spectral_scale * np.linalg.norm(rigid, axis=0), np.finfo(np.float64).tiny
    )
    rigid_tolerance = rigid_factor * np.finfo(np.float64).eps * DOF_COUNT
    passed = (np.count_nonzero(near) == 6 and np.count_nonzero(positive) == 156
              and not np.any(negative) and np.isfinite(residuals).all()
              and np.all(residuals <= rigid_tolerance))
    if not passed:
        raise RuntimeError(
            "mechanical BLOCKED: scaled K must have exactly 6 expected rigid zeros, "
            "156 positive modes, no significant negative modes, and valid rigid residuals"
        )
    return {
        "scale": scale,
        "eigenvalues": eigenvalues,
        "zero_tolerance": tolerance,
        "zero_tolerance_formula": f"{zero_factor:g}*eps_float64*DOF_COUNT*spectral_scale",
        "near_zero_count": int(np.count_nonzero(near)),
        "positive_count": int(np.count_nonzero(positive)),
        "significant_negative_count": int(np.count_nonzero(negative)),
        "strict_negative_count": int(np.count_nonzero(eigenvalues < 0.0)),
        "rigid_scaled_residuals": residuals,
        "rigid_scaled_residual_tolerance": rigid_tolerance,
        "rigid_residual_tolerance_formula": f"{rigid_factor:g}*eps_float64*DOF_COUNT",
        "scaled_asymmetry_relative": scaled_asymmetry,
        "scaled_asymmetry_relative_tolerance": asymmetry_limit,
    }


def _solve(np: object, stiffness: object, force: object, constraint: object, target: object,
           spectrum: dict[str, object], trust_gates: dict[str, object]) -> tuple[object, object, dict[str, float]]:
    scale = spectrum["scale"]
    scaled_k = scale[:, None] * stiffness * scale[None, :]
    row_norm = np.linalg.norm(constraint * scale[None, :], axis=1)
    if not np.isfinite(row_norm).all() or np.any(row_norm <= 0.0):
        raise RuntimeError("a physical volume-moment constraint is singular")
    cscale = 1.0 / row_norm
    scaled_c = cscale[:, None] * constraint * scale[None, :]
    singular_values = np.linalg.svd(scaled_c, compute_uv=False)
    # Rank is checked on the six moment equations themselves, not by projecting
    # any force or coefficient vector.
    if singular_values[-1] <= 64.0 * np.finfo(np.float64).eps * singular_values[0]:
        raise RuntimeError("physical volume-moment rows do not have full rank")
    saddle = np.zeros((DOF_COUNT + 6, DOF_COUNT + 6), dtype=np.float64)
    saddle[:DOF_COUNT, :DOF_COUNT] = scaled_k
    saddle[:DOF_COUNT, DOF_COUNT:] = scaled_c.T
    saddle[DOF_COUNT:, :DOF_COUNT] = scaled_c
    rhs = np.concatenate((scale * force, cscale * target))
    solution = np.linalg.solve(saddle, rhs)
    if not np.isfinite(solution).all():
        raise RuntimeError("direct 168x168 saddle solve returned non-finite values")
    coeff = scale * solution[:DOF_COUNT]
    multipliers = cscale * solution[DOF_COUNT:]
    equilibrium = stiffness @ coeff + constraint.T @ multipliers - force
    constraint_residual = constraint @ coeff - target
    saddle_backward = float(np.linalg.norm(saddle @ solution - rhs) /
        max(np.linalg.norm(rhs), np.finfo(np.float64).tiny))
    weak = float(np.linalg.norm(scale * equilibrium) / max(
        np.linalg.norm(scale * (stiffness @ coeff)) + np.linalg.norm(scale * (constraint.T @ multipliers))
        + np.linalg.norm(scale * force), np.finfo(np.float64).tiny))
    gauge_res = float(np.linalg.norm(constraint_residual) / max(
        np.linalg.norm(target), np.linalg.norm(constraint @ coeff), np.finfo(np.float64).tiny))
    strain_energy = 0.5 * float(coeff @ stiffness @ coeff)
    external_work = float(coeff @ force)
    multiplier_work = float(target @ multipliers)
    identity = 2.0 * strain_energy - external_work + multiplier_work
    identity_relative = abs(identity) / max(abs(2.0 * strain_energy), abs(external_work),
                                             abs(multiplier_work), np.finfo(np.float64).tiny)
    metrics = {"weak_equilibrium_residual_relative": weak,
               "gauge_residual_relative": gauge_res,
               "saddle_backward_residual_relative": saddle_backward,
               "strain_energy_J": strain_energy,
               "external_work_J": external_work,
               "gauge_multiplier_work_J": multiplier_work,
               "potential_U_minus_W_J": strain_energy - external_work,
               "energy_identity_defect_J": identity,
               "energy_identity_relative": identity_relative}
    if not all(math.isfinite(v) for v in metrics.values()):
        raise FloatingPointError("solve residual/energy metrics are non-finite")
    gate_values = {
        "weak_equilibrium_residual_relative_max": weak,
        "gauge_residual_relative_max": gauge_res,
        "saddle_backward_residual_relative_max": saddle_backward,
        "energy_identity_relative_max": identity_relative,
    }
    exceeded = {name: {"value": value, "limit": float(trust_gates[name])}
                for name, value in gate_values.items() if value > float(trust_gates[name])}
    if exceeded:
        raise RuntimeError(f"fixed-front algebraic solve trust gate failed: metrics={metrics}; exceeded={exceeded}")
    return coeff, multipliers, metrics


def _execute(output_dir: Path) -> dict[str, object]:
    output_dir = output_dir.resolve()
    allowed_root = STUDY_OUTPUT_ROOT.resolve()
    try:
        output_dir.relative_to(allowed_root)
    except ValueError as exc:
        raise ValueError(f"study output path must stay within its own review directory: {allowed_root}") from exc
    if output_dir.exists():
        if not output_dir.is_dir():
            raise FileExistsError("study output path exists and is not a directory")
        if any(output_dir.iterdir()):
            raise FileExistsError("study output directory is nonempty; refusing overwrite before imports/queries")
    if (output_dir / "arrays.npz").exists() or (output_dir / "result.json").exists():
        raise FileExistsError("single-attempt study outputs already exist; refusing overwrite")
    started = time.monotonic()
    case = json.loads(CASE_PATH.read_text(encoding="utf-8"))
    if case.get("case_id") != CASE_ID:
        raise ValueError("study case ID mismatch")
    if case.get("execution_authorization", {}).get("authorized") is not True:
        raise PermissionError("numerical execution is blocked until root explicitly authorizes the frozen source")
    input_hashes = _verify_inputs()
    deps = _load_numeric_dependencies()
    np = deps["np"]
    fixed = json.loads(FIXED_CASE_PATH.read_text(encoding="utf-8"))
    source_axis = np.asarray([-0.005, 0.0, 0.005], dtype=np.float64)
    source_xyz = np.ascontiguousarray(
        [(x, y, z) for x in source_axis for y in source_axis for z in source_axis], dtype=np.float64
    )
    source_ids = tuple(f"node-{row:03d}" for row in range(SOURCE_COUNT))
    registry = deps["build_full3d_registry"](source_ids)
    keys = registry.ordered_keys
    if registry.dof_count != DOF_COUNT:
        raise RuntimeError("semantic registry did not freeze exactly 162 vector DOFs")
    basis_ids = tuple((sid, basis) for sid in source_ids for basis in
                      (deps["BasisKind"].BASE, deps["BasisKind"].FRONT_F1_UNSHIFTED))
    basis_rows = {value: row for row, value in enumerate(basis_ids)}
    component_index = {deps["FieldComponent"].UX: 0, deps["FieldComponent"].UY: 1,
                       deps["FieldComponent"].UZ: 2}
    slot_lookup = {(key.source_id, key.field_component, key.basis_kind): registry.slot_of(key) for key in keys}
    groups = _declared_groups(np, fixed)
    unique_xyz, group_indices = _exact_unique_map(np, groups)
    raw_rows = sum(int(group["xyz"].shape[0]) for group in groups.values())
    if raw_rows != 4992 or unique_xyz.shape[0] > raw_rows:
        raise RuntimeError("global point registry violates the 4,992-row predeclaration")
    query_batches = []
    base_values_parts, base_gradient_parts = [], []
    exact_u_parts, exact_grad_parts, exact_stress_parts, exact_energy_parts = [], [], [], []
    mls_condition_parts, mls_support_parts, mls_status_parts, field_status_parts = [], [], [], []
    condition_max = 0.0
    field_case = deps["ModeIIIFixedCase"](radius_m=RADIUS_M, length_m=LENGTH_M,
        K_III_Pa_sqrt_m=100.0, young_modulus_Pa=1_000_000.0, poisson_ratio=0.3)
    config = deps["MLS3DConfig"](reference_origin=(0.0, 0.0, 0.0), reference_scale=RADIUS_M,
        support_radius=0.025, max_condition_number=1.0e12, rank_relative_tolerance=1.0e-12,
        switch_tolerance=1.0e-10, duplicate_tolerance=1.0e-12)
    for start in range(0, unique_xyz.shape[0], MAX_QUERY_BATCH):
        stop = min(start + MAX_QUERY_BATCH, unique_xyz.shape[0])
        batch = np.ascontiguousarray(unique_xyz[start:stop])
        raw = deps["ContinuousMLSFull3DFrontProvider"](
            source_xyz, source_ids, batch, registry, config, field_case
        )
        if raw.unknown_rows or raw.values is None or raw.gradients is None:
            raise RuntimeError(f"MLS/reference batch has blocked rows: {raw.unknown_rows[:4]}")
        base_values_parts.append(np.asarray(raw.mls_result.weights, dtype=np.float64))
        base_gradient_parts.append(np.asarray(raw.mls_result.gradients, dtype=np.float64))
        ref = raw.field_result
        exact_u_parts.append(np.asarray(ref.displacement_m, dtype=np.float64))
        exact_grad_parts.append(np.asarray(ref.displacement_gradient, dtype=np.float64))
        exact_stress_parts.append(np.asarray(ref.stress_Pa, dtype=np.float64))
        exact_energy_parts.append(np.asarray(ref.energy_density_J_per_m3, dtype=np.float64))
        mls_condition_parts.append(np.asarray(raw.mls_result.condition_numbers, dtype=np.float64))
        mls_support_parts.append(np.asarray(raw.mls_result.support_counts, dtype=np.int64))
        mls_status_parts.extend(status.value for status in raw.mls_result.statuses)
        field_status_parts.extend(str(status) for status in ref.statuses)
        condition_max = max(condition_max, float(np.max(raw.mls_result.condition_numbers)))
        query_batches.append({"start": start, "stop": stop, "count": stop - start,
                              "mls_status": "OK", "field_status": "OK"})
        if time.monotonic() - started > MAX_WALL_SECONDS:
            raise TimeoutError("60-second study wall limit reached during MLS/reference batches")
    base_values = np.vstack(base_values_parts)
    base_gradients = np.vstack(base_gradient_parts)
    exact_u = np.vstack(exact_u_parts)
    exact_gradient = np.vstack(exact_grad_parts)
    exact_stress = np.vstack(exact_stress_parts)
    exact_energy_density = np.concatenate(exact_energy_parts)
    mls_condition_numbers = np.concatenate(mls_condition_parts)
    mls_support_counts = np.concatenate(mls_support_parts)
    if unique_xyz.shape != (len(base_values), 3) or not all(np.isfinite(a).all() for a in
            (base_values, base_gradients, exact_u, exact_gradient, exact_stress, exact_energy_density)):
        raise RuntimeError("query replay arrays have invalid shapes or non-finite values")

    material = deps["IsotropicMaterial"](young_modulus=field_case.young_modulus_Pa,
        poisson_ratio=field_case.poisson_ratio)
    D = material.constitutive_matrix
    raw_rigid = deps["rigid_mode_coefficient_vectors"](source_xyz, source_ids, registry)
    amp = field_case.amplitude_sqrt_m
    gauge_target = np.asarray((0.0, 0.0, 0.0, amp * math.sqrt(RADIUS_M) * 16.0 / (21.0 * math.pi), 0.0, 0.0))
    if not math.isclose(float(gauge_target[3]), 5.031114771170333e-6, rel_tol=0.0, abs_tol=2e-20):
        raise RuntimeError("closed-form physical g_rx differs from the frozen source fixture")

    results: dict[str, dict[str, object]] = {}
    npz_arrays: dict[str, object] = {
        "source_positions_xyz_m": source_xyz,
        "unique_query_positions_xyz_m": unique_xyz,
        "base_mls_values": base_values,
        "base_mls_gradients_per_m": base_gradients,
        "exact_displacement_m": exact_u,
        "exact_displacement_gradient_m_per_m": exact_gradient,
        "exact_stress_Pa": exact_stress,
        "exact_energy_density_J_per_m3": exact_energy_density,
        "rigid_mode_coefficient_vectors": raw_rigid,
        "gauge_target_m": gauge_target,
        "mls_condition_numbers": mls_condition_numbers,
        "mls_support_counts": mls_support_counts,
        "mls_status_utf8": np.asarray(mls_status_parts),
        "field_status_utf8": np.asarray(field_status_parts),
        "registry_canonical_ids_utf8": np.asarray(registry.canonical_ids),
        "registry_ordered_keys_json": np.asarray([
            json.dumps([key.source_id, key.field_component.value, key.basis_kind.value, key.basis_index], separators=(",", ":"))
            for key in keys
        ]),
        "basis_row_ids_json": np.asarray([
            json.dumps([sid, basis.value], separators=(",", ":")) for sid, basis in basis_ids
        ]),
    }
    for name, group in groups.items():
        safe = name.replace(".", "_")
        npz_arrays[f"cloud_{safe}_xyz_m"] = group["xyz"]
        if "weights" in group:
            npz_arrays[f"cloud_{safe}_weights"] = group["weights"]
        if "normals" in group:
            npz_arrays[f"cloud_{safe}_normals"] = group["normals"]
        npz_arrays[f"cloud_{safe}_unique_row_ids"] = group_indices[name]

    for a in TRIAL_OFFSETS_M:
        for theta_order in VOLUME_THETA_ORDERS:
            case_name = f"a{a:g}_theta{theta_order}"
            volume_key = f"a{a:g}.volume{theta_order}"
            volume_group = groups[volume_key]
            vi = group_indices[volume_key]
            vxyz = unique_xyz[vi]
            vw = volume_group["weights"]
            vol_values, vol_grads = _shifted_basis(np, base_values[vi], base_gradients[vi], vxyz, a)
            # Assembly calls are deliberately chunked to the unchanged frozen
            # 256-point API. Each chunk's complete raw K is accumulated once.
            K = np.zeros((DOF_COUNT, DOF_COUNT), dtype=np.float64)
            chunk_records = []
            for start in range(0, len(vi), MAX_QUERY_BATCH):
                stop = min(start + MAX_QUERY_BATCH, len(vi))
                adapter = _RowsProvider(np, deps["ShapeSample"], deps["ShapeStatus"],
                    vol_values[start:stop], vol_grads[start:stop], vxyz[start:stop], basis_ids)
                assembled = deps["assemble_full3d_enriched_operator"](
                    source_xyz, source_ids, np.ascontiguousarray(vxyz[start:stop]),
                    np.ascontiguousarray(vw[start:stop]), material, adapter, registry,
                )
                if assembled.status.value != "OK" or assembled.stiffness is None or adapter.cursor != stop - start:
                    raise RuntimeError(f"chunked full-3D assembly blocked at {case_name} chunk {start}:{stop}")
                K += assembled.stiffness
                chunk_records.append({"start": start, "stop": stop, "count": stop - start, "status": "OK"})
            if not np.isfinite(K).all():
                raise FloatingPointError("summed raw semantic stiffness is non-finite")

            side_idx = group_indices["side"]
            c8_idx = group_indices[f"a{a:g}.cap8"]
            c16_idx = group_indices[f"a{a:g}.cap16"]
            off_idx = group_indices["offgrid"]
            side_group, c8_group, c16_group, off_group = groups["side"], groups[f"a{a:g}.cap8"], groups[f"a{a:g}.cap16"], groups["offgrid"]
            side_values, side_grads = _shifted_basis(np, base_values[side_idx], base_gradients[side_idx], unique_xyz[side_idx], a)
            c8_values, c8_grads = _shifted_basis(np, base_values[c8_idx], base_gradients[c8_idx], unique_xyz[c8_idx], a)
            c16_values, c16_grads = _shifted_basis(np, base_values[c16_idx], base_gradients[c16_idx], unique_xyz[c16_idx], a)
            off_values, off_grads = _shifted_basis(np, base_values[off_idx], base_gradients[off_idx], unique_xyz[off_idx], a)

            def load_for(values: object, idx: object, group: dict[str, object]) -> tuple[object, object, dict[str, object]]:
                normals = group["normals"]
                traction = np.einsum("qij,qj->qi", exact_stress[idx], normals)
                f = np.zeros(DOF_COUNT, dtype=np.float64)
                for row, (sid, basis) in enumerate(basis_ids):
                    for component_enum, component in component_index.items():
                        slot = slot_lookup[(sid, component_enum, basis)]
                        f[slot] = np.dot(group["weights"] * traction[:, component], values[:, row])
                if "region_slices" in group:
                    region_wrenches = {}
                    for region_name, limits in group["region_slices"].items():
                        region = slice(limits[0], limits[1])
                        region_wrenches[region_name] = _wrench(
                            np, group["xyz"][region], traction[region], group["weights"][region]
                        )
                    wrench = {
                        "force_N": region_wrenches["bottom"]["force_N"] + region_wrenches["top"]["force_N"],
                        "moment_about_origin_Nm": region_wrenches["bottom"]["moment_about_origin_Nm"] + region_wrenches["top"]["moment_about_origin_Nm"],
                        "by_cap": region_wrenches,
                    }
                else:
                    wrench = _wrench(np, group["xyz"], traction, group["weights"])
                return f, traction, wrench

            fside, side_t, side_wrench = load_for(side_values, side_idx, side_group)
            fc8, c8_t, c8_wrench = load_for(c8_values, c8_idx, c8_group)
            fc16, c16_t, c16_wrench = load_for(c16_values, c16_idx, c16_group)
            force = fside + fc16
            constraint, centroid, physical_rigid, measure = _volume_moments(
                np, vol_values, vw, vxyz, keys, basis_rows, registry, source_xyz, source_ids,
                raw_rigid, component_index, deps["DofKey"], deps["BasisKind"], tuple(component_index),
            )
            gauge_on_rigid = constraint @ physical_rigid
            gauge_singular = np.linalg.svd(gauge_on_rigid, compute_uv=False)
            if gauge_singular[-1] <= 64.0 * np.finfo(np.float64).eps * gauge_singular[0]:
                raise RuntimeError("six volume moments do not control the physical rigid modes")
            spectrum = _spectrum(np, K, raw_rigid, case["mechanical_spectrum_gates"])
            coeff, multipliers, solve_metrics = _solve(np, K, force, constraint, gauge_target, spectrum,
                case["numerical_trust_gates"])
            solved_u, solved_grad = _reconstruct(np, keys, basis_rows, vol_values, vol_grads, coeff, component_index)
            solved_strain = _strain(np, solved_grad)
            solved_stress = _stress_tensor(np, solved_strain @ D.T)
            reconstructed_energy = 0.5 * float(np.dot(vw, np.einsum("qi,qi->q", solved_strain @ D, solved_strain)))
            direct_energy = 0.5 * float(coeff @ K @ coeff)
            energy_relative = abs(reconstructed_energy - direct_energy) / max(abs(direct_energy), np.finfo(np.float64).tiny)
            solve_metrics["reconstructed_strain_energy_J"] = reconstructed_energy
            solve_metrics["KTT_energy_relative_defect"] = energy_relative
            ktt_tolerance = float(case["numerical_trust_gates"]["ktt_energy_relative_defect_max"])
            if not math.isfinite(energy_relative) or energy_relative > ktt_tolerance:
                raise RuntimeError(
                    f"fixed-front KTT energy trust gate failed: metrics={solve_metrics}; limit={ktt_tolerance}"
                )

            exact_coefficient = np.zeros(DOF_COUNT, dtype=np.float64)
            if a == 0.0:
                for slot, key in enumerate(keys):
                    if key.basis_kind is deps["BasisKind"].FRONT_F1_UNSHIFTED and key.field_component is deps["FieldComponent"].UZ:
                        exact_coefficient[slot] = amp * math.sqrt(RADIUS_M)
            exact_quad_energy = float(np.dot(vw, exact_energy_density[vi]))
            fine_by_cap = c16_wrench["by_cap"]
            coarse_by_cap = c8_wrench["by_cap"]
            total_force = side_wrench["force_N"] + c16_wrench["force_N"]
            total_moment = side_wrench["moment_about_origin_Nm"] + c16_wrench["moment_about_origin_Nm"]
            q_load = material.shear_modulus * amp * RADIUS_M**1.5
            analytic_wrench = {
                "side_force_N": [0.0, 0.0, 0.0],
                "side_moment_about_origin_Nm": [4.0 / 3.0 * q_load * LENGTH_M, 0.0, 0.0],
                "top_force_N": [0.0, 4.0 / 3.0 * q_load, 0.0],
                "top_moment_about_origin_Nm": [-2.0 / 3.0 * q_load * LENGTH_M, 0.0, 4.0 / 5.0 * q_load * RADIUS_M],
                "bottom_force_N": [0.0, -4.0 / 3.0 * q_load, 0.0],
                "bottom_moment_about_origin_Nm": [-2.0 / 3.0 * q_load * LENGTH_M, 0.0, -4.0 / 5.0 * q_load * RADIUS_M],
                "total_force_N": [0.0, 0.0, 0.0], "total_moment_about_origin_Nm": [0.0, 0.0, 0.0],
            }
            case_diag: dict[str, object] = {
                "status": "OK_DIAGNOSTIC_ONLY_UNQUALIFIED",
                "a_m": a, "volume_theta_order": theta_order,
                "volume_rows": len(vi), "volume_measure_m3": measure,
                "volume_centroid_xyz_m": centroid,
                "assembly_chunks": chunk_records,
                "assembly_raw_K_symmetry_relative": float(np.linalg.norm(K - K.T) / max(np.linalg.norm(K), np.finfo(np.float64).tiny)),
                "spectrum": {k: v for k, v in spectrum.items() if k not in {"scale", "eigenvalues", "rigid_scaled_residuals"}},
                "six_rigid_scaled_residuals": spectrum["rigid_scaled_residuals"],
                "eigenvalue_min_scaled": float(np.min(spectrum["eigenvalues"])),
                "eigenvalue_max_scaled": float(np.max(spectrum["eigenvalues"])),
                "gauge_target_m": gauge_target,
                "gauge_constraint_singular_values": gauge_singular,
                "solve": solve_metrics,
                "six_rigid_virtual_work_N_or_Nm": raw_rigid.T @ force,
                "boundary_wrench": {"side": side_wrench, "top_coarse": coarse_by_cap["top"],
                    "bottom_coarse": coarse_by_cap["bottom"], "top_fine": fine_by_cap["top"],
                    "bottom_fine": fine_by_cap["bottom"], "caps_coarse": c8_wrench, "caps_fine": c16_wrench,
                    "total_fine_force_N": total_force, "total_fine_moment_about_origin_Nm": total_moment,
                    "analytic_full_traction_reference": analytic_wrench,
                    "cap_force_delta_N": c16_wrench["force_N"] - c8_wrench["force_N"],
                    "cap_moment_delta_Nm": c16_wrench["moment_about_origin_Nm"] - c8_wrench["moment_about_origin_Nm"]},
                "cap_generalized_load_delta_norm_N": float(np.linalg.norm(fc16 - fc8)),
                "cap_generalized_load_delta_relative": float(np.linalg.norm(fc16 - fc8) / max(np.linalg.norm(fc16), np.finfo(np.float64).tiny)),
                "cap_external_work_coarse_J": float(coeff @ fc8),
                "cap_external_work_fine_J": float(coeff @ fc16),
                "exact_old_field_energy_quadrature_J": exact_quad_energy,
                "manufactured_coefficients_energy_J": 0.5 * float(exact_coefficient @ K @ exact_coefficient) if a == 0.0 else None,
            }
            if a == 0.0:
                ref_u, ref_grad, ref_stress = exact_u[vi], exact_gradient[vi], exact_stress[vi]
                case_diag["old_tip_weighted_relative_errors"] = {
                    "displacement": _weighted_relative(np, solved_u - ref_u, ref_u, vw),
                    "displacement_gradient": _weighted_relative(np, solved_grad - ref_grad, ref_grad, vw),
                    "full_stress": _weighted_relative(np, solved_stress - ref_stress, ref_stress, vw),
                }
            normal_sq = solved_stress[:, (0, 1, 2), (0, 1, 2)]
            shear_sq = solved_stress[:, (1, 2), (2, 0)]
            normal_and_xy = np.column_stack((normal_sq, solved_stress[:, 0, 1]))
            case_diag["solved_stress_diagnostics"] = {
                "normal_to_mode_III_shear_L2_ratio": math.sqrt(float(np.dot(vw, np.sum(normal_and_xy**2, axis=1))) /
                    max(float(np.dot(vw, np.sum(shear_sq**2, axis=1))), np.finfo(np.float64).tiny)),
                "max_abs_sigma_yy_Pa": float(np.max(np.abs(solved_stress[:, 1, 1]))),
                "exact_new_tip_error_claimed": False,
            }
            off_u, off_grad = _reconstruct(np, keys, basis_rows, off_values, off_grads, coeff, component_index)
            off_strain = _strain(np, off_grad)
            off_stress = _stress_tensor(np, off_strain @ D.T)
            if a == 0.0:
                case_diag["offgrid_old_tip_weighted_relative_errors"] = {
                    "displacement": _weighted_relative(np, off_u - exact_u[off_idx], exact_u[off_idx], np.ones(len(off_idx))),
                    "displacement_gradient": _weighted_relative(np, off_grad - exact_gradient[off_idx], exact_gradient[off_idx], np.ones(len(off_idx))),
                    "full_stress": _weighted_relative(np, off_stress - exact_stress[off_idx], exact_stress[off_idx], np.ones(len(off_idx))),
                }

            results[case_name] = case_diag
            case_prefix = case_name + "__"
            npz_arrays.update({
                case_prefix + "volume_xyz_m": vxyz, case_prefix + "volume_weights_m3": vw,
                case_prefix + "volume_basis_values": vol_values, case_prefix + "volume_basis_gradients_per_m": vol_grads,
                case_prefix + "K_N_per_m": K, case_prefix + "boundary_load_N": force,
                case_prefix + "constraint_matrix": constraint, case_prefix + "rigid_mode_coefficients": raw_rigid,
                case_prefix + "coefficients_m": coeff, case_prefix + "multipliers": multipliers,
                case_prefix + "eigenvalues_scaled": spectrum["eigenvalues"],
                case_prefix + "solved_volume_displacement_m": solved_u,
                case_prefix + "solved_volume_gradient_m_per_m": solved_grad,
                case_prefix + "solved_volume_stress_Pa": solved_stress,
                case_prefix + "offgrid_xyz_m": groups["offgrid"]["xyz"],
                case_prefix + "offgrid_basis_values": off_values,
                case_prefix + "offgrid_basis_gradients_per_m": off_grads,
                case_prefix + "solved_offgrid_displacement_m": off_u,
                case_prefix + "solved_offgrid_gradient_m_per_m": off_grad,
                case_prefix + "solved_offgrid_stress_Pa": off_stress,
                case_prefix + "side_basis_values": side_values,
                case_prefix + "side_basis_gradients_per_m": side_grads,
                case_prefix + "cap8_basis_values": c8_values,
                case_prefix + "cap8_basis_gradients_per_m": c8_grads,
                case_prefix + "cap16_basis_values": c16_values,
                case_prefix + "cap16_basis_gradients_per_m": c16_grads,
                case_prefix + "side_load_N": fside,
                case_prefix + "cap8_load_N": fc8,
                case_prefix + "cap16_load_N": fc16,
                case_prefix + "side_traction_Pa": side_t,
                case_prefix + "cap8_traction_Pa": c8_t,
                case_prefix + "cap16_traction_Pa": c16_t,
            })
            if time.monotonic() - started > MAX_WALL_SECONDS:
                raise TimeoutError("60-second study wall limit reached during case assembly/solve")

    # Finite delta-a energy differences are computed from solved U-W only;
    # no contour or additional field queries are introduced.
    cross_case = {"exact_old_tip_local_asymptotic_G_J_per_m2": 0.013,
                  "prior_finite_offset_G_context_J_per_m2": 0.005466,
                  "finite_delta_a_is_not_percent_compared_to_local_G": True}
    for theta in VOLUME_THETA_ORDERS:
        name0, name5 = f"a0_theta{theta}", f"a0.0005_theta{theta}"
        base = results[name0]["solve"]["potential_U_minus_W_J"]
        shifted = results[name5]["solve"]["potential_U_minus_W_J"]
        v0 = npz_arrays[name0 + "__solved_offgrid_displacement_m"]
        v5 = npz_arrays[name5 + "__solved_offgrid_displacement_m"]
        grad0 = npz_arrays[name0 + "__solved_offgrid_gradient_m_per_m"]
        grad5 = npz_arrays[name5 + "__solved_offgrid_gradient_m_per_m"]
        stress0 = npz_arrays[name0 + "__solved_offgrid_stress_Pa"]
        stress5 = npz_arrays[name5 + "__solved_offgrid_stress_Pa"]
        ones = np.ones(v0.shape[0], dtype=np.float64)
        cross_case[f"theta{theta}"] = {
            "G_delta_a_J_per_m2": -(shifted - base) / (LENGTH_M * 0.0005),
            "Pi_a0_J": base, "Pi_a0.5mm_J": shifted,
            "fixed_physical_offgrid_solution_change_a0_to_a0.5mm": {
                "displacement_relative_to_a0": _weighted_relative(np, v5 - v0, v0, ones),
                "gradient_relative_to_a0": _weighted_relative(np, grad5 - grad0, grad0, ones),
                "stress_relative_to_a0": _weighted_relative(np, stress5 - stress0, stress0, ones),
            },
        }
    theta_changes = {}
    for a in TRIAL_OFFSETS_M:
        name32, name64 = f"a{a:g}_theta32", f"a{a:g}_theta64"
        k32, k64 = npz_arrays[name32 + "__K_N_per_m"], npz_arrays[name64 + "__K_N_per_m"]
        f32, f64 = npz_arrays[name32 + "__boundary_load_N"], npz_arrays[name64 + "__boundary_load_N"]
        u32, u64 = npz_arrays[name32 + "__solved_offgrid_displacement_m"], npz_arrays[name64 + "__solved_offgrid_displacement_m"]
        theta_changes[f"a{a:g}"] = {
            "raw_K_relative_change": float(np.linalg.norm(k64 - k32) / max(np.linalg.norm(k64), np.finfo(np.float64).tiny)),
            "fine_boundary_load_relative_change": float(np.linalg.norm(f64 - f32) / max(np.linalg.norm(f64), np.finfo(np.float64).tiny)),
            "Pi_theta64_minus_theta32_J": results[name64]["solve"]["potential_U_minus_W_J"] - results[name32]["solve"]["potential_U_minus_W_J"],
            "offgrid_displacement_relative_change": _weighted_relative(np, u64 - u32, u32, np.ones(u32.shape[0])),
        }
    cross_case["theta32_to_theta64_at_fixed_a_and_fine_boundary"] = theta_changes
    result = {
        "schema": "tip_aligned_full3d_quadrature_study_result.v1",
        "case_id": CASE_ID,
        "status": "OK_DIAGNOSTIC_ONLY_UNQUALIFIED",
        "classification": "SYNTHETIC_MODE_III_ENGINEERING_BENCHMARK_NOT_CERAMIC_MATERIAL",
        "assembly_status": "OK",
        "mechanical_eligibility": "SPECTRUM_GATE_PASSED",
        "solve_allowed": True,
        "study_case_sha256": _sha256(CASE_PATH),
        "source_sha256": {"driver": _sha256(Path(__file__)), **input_hashes},
        "registry": {
            "source_ids": list(source_ids),
            "ordered_canonical_ids": list(registry.canonical_ids),
            "ordered_keys": [[key.source_id, key.field_component.value, key.basis_kind.value, key.basis_index] for key in keys],
            "basis_row_ids": [[sid, basis.value] for sid, basis in basis_ids],
            "component_order": ["ux", "uy", "uz"],
            "basis_order": ["base", "front_f1_unshifted"],
        },
        "region_slices": {
            "semantics": "half-open row slices local to each named side/cap cloud",
            "side": [0, 256],
            "cap8_bottom": groups["a0.cap8"]["region_slices"]["bottom"],
            "cap8_top": groups["a0.cap8"]["region_slices"]["top"],
            "cap16_bottom": groups["a0.cap16"]["region_slices"]["bottom"],
            "cap16_top": groups["a0.cap16"]["region_slices"]["top"],
        },
        "query_budget": {"raw_requested_rows": raw_rows, "unique_xyz_rows": int(len(unique_xyz)),
            "exact_duplicate_rows_elided_before_provider": raw_rows - int(len(unique_xyz)),
            "provider_max_batch_rows": MAX_QUERY_BATCH, "provider_batches": query_batches,
            "provider_batch_count": len(query_batches),
            "side_max_inward_roundoff_correction_m": groups["side"]["max_inward_roundoff_correction_m"],
            "raw_group_row_counts": {name: int(group["xyz"].shape[0]) for name, group in groups.items()}},
        "provider_reference_status": {
            "mls_status_counts": {"OK": len(mls_status_parts)},
            "old_tip_reference_status_counts": {"OK": len(field_status_parts)},
            "mls_condition_number_min": float(np.min(mls_condition_numbers)),
            "mls_condition_number_max": condition_max,
            "mls_support_count_min": int(np.min(mls_support_counts)),
            "mls_support_count_max": int(np.max(mls_support_counts)),
        },
        "analytic_gauge_rx_m": float(gauge_target[3]),
        "units": {
            "coordinates": "m", "volume_weights": "m^3", "area_weights": "m^2",
            "shape_values": "dimensionless", "basis_gradients": "1/m",
            "displacement": "m", "displacement_gradient": "m/m (dimensionless)",
            "stress_and_traction": "Pa", "stiffness": "N/m", "boundary_load": "N",
            "rigid_mode_vector_columns": {"translations": "unit displacement coefficient per unit translation",
                                           "rotations": "m displacement coefficient per radian"},
        },
        "cases": results,
        "cross_case_finite_delta_a_energy": cross_case,
        "array_budget_bytes": {
            "persisted_numeric_and_index_arrays": int(sum(array.nbytes for array in npz_arrays.values())),
            "conservative_transient_assembly_and_solve_workspace_allowance": 8 * 1024 * 1024,
            "conservative_peak_estimate": int(sum(array.nbytes for array in npz_arrays.values())) + 8 * 1024 * 1024,
        },
        "limitations": case["limitations"],
        "elapsed_seconds": time.monotonic() - started,
    }
    array_bytes = result["array_budget_bytes"]["conservative_peak_estimate"]
    if array_bytes > MAX_ARRAY_BYTES:
        raise MemoryError(f"explicit ndarray budget exceeded: {array_bytes} > {MAX_ARRAY_BYTES}")
    output_dir.mkdir(parents=True, exist_ok=True)
    npz_path = output_dir / "arrays.npz"
    json_path = output_dir / "result.json"
    if npz_path.exists() or json_path.exists():
        raise FileExistsError("single-attempt output already exists; refusing overwrite/retry")
    np.savez_compressed(npz_path, **npz_arrays)
    result["arrays_npz"] = {"path": npz_path.name, "sha256": _sha256(npz_path), "bytes": npz_path.stat().st_size}
    json_bytes = json.dumps(_jsonable(np, result), indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    json_path.write_text(json_bytes, encoding="utf-8")
    combined_bytes = npz_path.stat().st_size + json_path.stat().st_size
    if combined_bytes > MAX_OUTPUT_BYTES:
        raise MemoryError(f"combined NPZ/JSON output exceeded {MAX_OUTPUT_BYTES} bytes")
    result["measured_output_bytes"] = {"npz": npz_path.stat().st_size, "json": json_path.stat().st_size,
                                      "combined": combined_bytes}
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-check", action="store_true", help="run the standard-library source/count fixture only")
    parser.add_argument("--execute", action="store_true", help="future single bounded numerical study; requires separate authorization")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    if args.source_check:
        print(json.dumps(_source_check(), sort_keys=True))
        return 0
    if not args.execute:
        parser.error("select --source-check, or use --execute only after separate authorization")
    try:
        result = _execute(args.output_dir.resolve())
    except Exception as exc:
        # The parent process is responsible for preserving the original return
        # code/stdout/stderr under the one external 60-second timeout.
        print(json.dumps({"status": "STUDY_FAILED_NO_RETRY", "error_type": type(exc).__name__,
                          "error": str(exc), "traceback": traceback.format_exc()[-4000:]}, ensure_ascii=False))
        return 2
    print(json.dumps({"status": result["status"], "result_path": str(args.output_dir.resolve() / "result.json"),
                      "array_sha256": result["arrays_npz"]["sha256"], "elapsed_seconds": result["elapsed_seconds"],
                      "combined_output_bytes": result["measured_output_bytes"]["combined"]}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
