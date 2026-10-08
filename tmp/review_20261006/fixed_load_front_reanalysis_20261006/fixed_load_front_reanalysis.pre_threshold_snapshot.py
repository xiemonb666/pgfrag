"""One pinned fixed-load, finite-front-extension reanalysis; source-only by default."""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import struct
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CASE_PATH = HERE / "fixed_load_front_reanalysis_case.json"
CASE_ID = "modeiii_pu_full3d_fixed_load_front_reanalysis_20261006"
OPERATOR_CASE_ID = "modeiii_pu_full3d_elastic_operator_20261006"
BASELINE_CASE_ID = "modeiii_pu_full3d_static_solve_bridge_20261006"
PINS_PATH = ROOT / "tmp/review_20261006/full3d_cracked_elastic_operator_20261006/SOURCE_PINS.json"
DEFAULT_OUTPUT = ROOT / "pgfrag/tmp/review_20261006/fixed_load_front_reanalysis_20261006/run01/result.json"
SOURCE_COUNT, DOF_COUNT, BASIS_COUNT = 27, 162, 54
VOLUME_COUNT, BOUNDARY_COUNT, QUERY_COUNT = 256, 128, 384
MAX_BYTES = 8 * 1024 * 1024
REGIONS = {"volume": [0, 256], "lateral_side": [256, 320], "top_cap": [320, 352], "bottom_cap": [352, 384]}

class Blocked(RuntimeError):
    pass

def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()

def _read(path: Path) -> tuple[dict, bytes]:
    try:
        if path.stat().st_size > MAX_BYTES:
            raise Blocked("pinned JSON exceeds 8 MiB")
        raw = path.read_bytes()
        value = json.loads(raw.decode("utf-8"), parse_constant=lambda x: (_ for _ in ()).throw(ValueError(x)))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise Blocked(f"cannot read valid pinned JSON ({type(exc).__name__})") from exc
    if len(raw) > MAX_BYTES or not isinstance(value, dict):
        raise Blocked("pinned JSON exceeds its cap or has a non-object root")
    return value, raw

def _case() -> dict:
    case, _ = _read(CASE_PATH)
    if case.get("case_id") != CASE_ID:
        raise Blocked("case ID mismatch")
    return case

def _resolve(value: str, label: str) -> Path:
    path = Path(value)
    return path.resolve() if path.is_absolute() else (ROOT / path).resolve()

def _output_target(value: str) -> tuple[Path, Path]:
    target = Path(value).resolve()
    allowed = (ROOT / "pgfrag/tmp/review_20261006/fixed_load_front_reanalysis_20261006").resolve()
    staging = target.with_name(target.name + ".staging")
    if allowed not in target.parents:
        raise Blocked("output path must remain in the prepared review directory")
    if target.exists() or staging.exists():
        raise Blocked("refusing to solve over a prior result or partial staging file")
    return target, staging

def _pinned_json(path_text: str, expected: object, label: str) -> tuple[dict, str]:
    if not isinstance(expected, str) or len(expected) != 64:
        raise Blocked(f"{label} SHA-256 pin is not frozen")
    value, raw = _read(_resolve(path_text, label))
    actual = _sha(raw)
    if actual.lower() != expected.lower():
        raise Blocked(f"{label} SHA-256 differs from its frozen pin")
    return value, actual

def _verify_sources(case: dict) -> dict[str, Path]:
    pins, raw = _read(PINS_PATH)
    frozen = case["source_pins_sha256"]
    if _sha(raw).lower() != frozen.lower():
        raise Blocked("frozen upstream SOURCE_PINS.json changed")
    entries = {item["name"]: item for item in pins.get("dependencies", []) if isinstance(item, dict)}
    required = ("enriched_dof_contract.py", "efgm_elastic_operator.py", "full3d_cracked_elastic_operator.py")
    paths: dict[str, Path] = {}
    for name in required:
        item = entries.get(name)
        if not isinstance(item, dict) or not isinstance(item.get("path"), str):
            raise Blocked(f"upstream pin is missing {name}")
        path = ROOT / item["path"]
        if not path.is_file() or _sha(path.read_bytes()).lower() != str(item.get("sha256", "")).lower():
            raise Blocked(f"upstream source differs from frozen pin: {name}")
        paths[name] = path
    return paths

def _numpy_modules(paths: dict[str, Path]):
    import numpy as np
    for path in paths.values():
        if str(path.parent) not in sys.path:
            sys.path.insert(0, str(path.parent))
    from enriched_dof_contract import BasisKind, FieldComponent
    from efgm_elastic_operator import IsotropicMaterial, ShapeSample, ShapeStatus
    from full3d_cracked_elastic_operator import (
        assemble_full3d_enriched_operator,
        basis_row_ids,
        build_full3d_registry,
        rigid_mode_coefficient_vectors,
    )
    return np, BasisKind, FieldComponent, IsotropicMaterial, ShapeSample, ShapeStatus, assemble_full3d_enriched_operator, basis_row_ids, build_full3d_registry, rigid_mode_coefficient_vectors

def _array(np, value: object, shape: tuple[int, ...], name: str):
    try:
        array = np.asarray(value, dtype=np.float64)
    except (TypeError, ValueError, OverflowError) as exc:
        raise Blocked(f"{name} is not a numeric array") from exc
    if array.shape != shape or not bool(np.isfinite(array).all()):
        raise Blocked(f"{name} has wrong shape or non-finite values")
    return np.ascontiguousarray(array)

def _basis_map(raw_rows: object, source_ids: list[str]) -> dict[tuple[str, str], int]:
    if not isinstance(raw_rows, list) or len(raw_rows) != BASIS_COUNT:
        raise Blocked("saved scalar basis row map is malformed")
    result = {}
    for row, item in enumerate(raw_rows):
        if isinstance(item, dict):
            pair = (item.get("source_id"), item.get("basis_kind"))
        elif isinstance(item, (list, tuple)) and len(item) == 2:
            pair = tuple(item)
        else:
            raise Blocked(f"malformed basis row ID at {row}")
        if pair[0] not in source_ids or pair[1] not in {"base", "front_f1_unshifted"} or pair in result:
            raise Blocked(f"invalid or duplicate basis row ID at {row}")
        result[pair] = row
    expected = {(source, family) for source in source_ids for family in ("base", "front_f1_unshifted")}
    if set(result) != expected:
        raise Blocked("saved basis rows do not map all source/family pairs")
    return result

def _refront(np, positions, old_values, old_gradients, old_rows, new_rows, radius, offset):
    values = np.empty((QUERY_COUNT, BASIS_COUNT), dtype=np.float64)
    gradients = np.empty((QUERY_COUNT, BASIS_COUNT, 3), dtype=np.float64)
    for q, (x, y, _z) in enumerate(positions):
        dx = float(x) - offset
        dy = float(y)
        r = math.hypot(dx, dy)
        if r == 0.0 or (dy == 0.0 and dx < 0.0):
            raise Blocked(f"UNKNOWN shifted-front singular/branch query at frozen row {q}")
        theta = math.atan2(dy, dx)
        root = math.sqrt(r / radius)
        fbar = root * math.sin(0.5 * theta)
        scale = 0.5 / math.sqrt(radius * r)
        grad_f = np.array((-scale * math.sin(0.5 * theta), scale * math.cos(0.5 * theta), 0.0))
        for (source, family), new_row in new_rows.items():
            base_row = old_rows[(source, "base")]
            phi = old_values[q, base_row]
            grad_phi = old_gradients[q, base_row]
            if family == "base":
                values[q, new_row] = phi
                gradients[q, new_row] = grad_phi
            else:
                values[q, new_row] = phi * fbar
                gradients[q, new_row] = fbar * grad_phi + phi * grad_f
    if not bool(np.isfinite(values).all()) or not bool(np.isfinite(gradients).all()):
        raise Blocked("shifted basis contains a non-finite value")
    return values, gradients

class _CachedProvider:
    def __init__(self, np, ShapeSample, ShapeStatus, row_ids, points, values, gradients):
        self.np, self.ShapeSample, self.ShapeStatus = np, ShapeSample, ShapeStatus
        self.basis_row_ids, self.points = row_ids, points
        self.values, self.gradients = values, gradients
        self.index = 0
        self.support = np.arange(BASIS_COUNT, dtype=np.int64)
    def evaluate(self, point):
        q = self.index
        if q >= len(self.points) or not self.np.array_equal(point, self.points[q]):
            raise Blocked(f"cached assembly query order diverged at row {q}")
        self.index += 1
        return self.ShapeSample(self.values[q], self.gradients[q], self.support, self.ShapeStatus.OK)

def _run(operator_result: dict, operator_sha: str, baseline: dict, case: dict, modules: tuple):
    np, BasisKind, FieldComponent, IsotropicMaterial, ShapeSample, ShapeStatus, assemble, row_ids_fn, build_registry, rigid_fn = modules
    started = time.monotonic()
    replay = operator_result.get("static_solve_replay")
    if not isinstance(replay, dict):
        raise Blocked("pinned operator result has no replay object")
    if (operator_result.get("case_id") != OPERATOR_CASE_ID or replay.get("case_id") != OPERATOR_CASE_ID or operator_result.get("assembly_status") != "OK"
            or operator_result.get("mechanical_eligibility") != "DIAGNOSTIC_PASS_UNQUALIFIED"
            or replay.get("schema") != "static_solve_replay.v1"):
        raise Blocked("pinned operator result is not an eligible static replay")
    spectrum = operator_result.get("spectrum", {})
    if spectrum.get("gate_passed") is not True or spectrum.get("observed_near_zero_count") != 6 or spectrum.get("observed_positive_count") != 156 or spectrum.get("significant_negative_count") != 0:
        raise Blocked("baseline operator spectrum gate did not pass")
    fixed = operator_result.get("fixed_case", {})
    geom, matcase = fixed.get("geometry", {}), fixed.get("material_and_manufactured_field", {})
    geom_expected, material_expected = case["geometry"], case["material"]
    if geom.get("radius_m") != geom_expected["radius_m"] or geom.get("length_m") != geom_expected["length_m"]:
        raise Blocked("pinned cylinder geometry differs from the frozen case")
    if matcase.get("young_modulus_Pa") != material_expected["young_modulus_Pa"] or matcase.get("poisson_ratio") != material_expected["poisson_ratio"]:
        raise Blocked("pinned elastic material differs from the frozen case")
    if replay.get("query_order") != ["volume", "lateral_side", "top_cap", "bottom_cap"]:
        raise Blocked("saved query order differs from the frozen volume/side/cap order")
    reg, src = replay.get("registry", {}), replay.get("sources", {})
    quad = replay.get("quadrature", {})
    op, basis, boundary, exact = replay.get("operator", {}), replay.get("basis", {}), replay.get("boundary", {}), replay.get("exact_field", {})
    ids = reg.get("source_ids")
    if not isinstance(ids, list) or ids != [f"node-{i:03d}" for i in range(SOURCE_COUNT)]:
        raise Blocked("source IDs differ from the frozen 27-site order")
    source_positions = _array(np, src.get("positions_xyz_m"), (SOURCE_COUNT, 3), "source positions")
    axis = (-0.005, 0.0, 0.005)
    import itertools
    expected_sources = np.asarray(list(itertools.product(axis, axis, axis)), dtype=np.float64)
    if not np.array_equal(source_positions, expected_sources):
        raise Blocked("rest source coordinates differ from the frozen 3x3x3 grid")
    positions = _array(np, quad.get("positions_xyz_m"), (QUERY_COUNT, 3), "saved 384 positions")
    volume_weights = _array(np, quad.get("volume_weights_m3"), (VOLUME_COUNT,), "volume weights")
    normals = _array(np, quad.get("boundary_normals"), (BOUNDARY_COUNT, 3), "boundary normals")
    area_weights = _array(np, quad.get("boundary_area_weights_m2"), (BOUNDARY_COUNT,), "boundary area weights")
    old_values = _array(np, basis.get("values"), (QUERY_COUNT, BASIS_COUNT), "saved basis values")
    old_gradients = _array(np, basis.get("gradients_per_m"), (QUERY_COUNT, BASIS_COUNT, 3), "saved basis gradients")
    traction = _array(np, boundary.get("traction_Pa"), (BOUNDARY_COUNT, 3), "frozen traction")
    _array(np, op.get("K_N_per_m"), (DOF_COUNT, DOF_COUNT), "baseline K")
    old_force = _array(np, op.get("boundary_load_N"), (DOF_COUNT,), "baseline load")
    old_rigid = _array(np, op.get("rigid_mode_coefficients"), (DOF_COUNT, 6), "baseline rigid modes")
    if quad.get("region_slices_half_open") != REGIONS or not bool(np.all(volume_weights > 0.0)) or not bool(np.all(area_weights > 0.0)):
        raise Blocked("saved quadrature partition or positive weights differ from the frozen 384-row data")
    if exact.get("statuses") != ["OK"] * QUERY_COUNT:
        raise Blocked("a saved field row is UNKNOWN")
    baseline_solve = baseline.get("solve", {})
    old_energy, old_work = float(baseline_solve.get("strain_energy_J", math.nan)), float(baseline_solve.get("external_work_J", math.nan))
    if not math.isfinite(old_energy) or not math.isfinite(old_work):
        raise Blocked("pinned baseline energy/work input is missing or non-finite")
    base_gauge = baseline.get("gauge")
    if not isinstance(base_gauge, dict): raise Blocked("pinned baseline result has no gauge object")
    target = _array(np, base_gauge.get("g_target_m"), (6,), "fixed baseline gauge target")
    if not target[3] > 0.0: raise Blocked("pinned baseline bridge has no positive Mode III r_rx target")
    evidence = replay.get("provider_evidence", {})
    condition_max = 1.0
    for name, count in (("volume", VOLUME_COUNT), ("boundary", BOUNDARY_COUNT)):
        batch = evidence.get(name, {})
        if batch.get("query_count") != count or batch.get("consumed_queries") != count or batch.get("mls_statuses") != ["OK"] * count or batch.get("field_statuses") != ["OK"] * count:
            raise Blocked(f"saved {name} provider evidence has an UNKNOWN row")
        cond = _array(np, batch.get("condition_numbers"), (count,), f"{name} MLS condition numbers")
        if not bool(np.all(cond >= 1.0)):
            raise Blocked(f"saved {name} MLS condition evidence is invalid")
        condition_max = max(condition_max, float(np.max(cond)))
    old_rows = _basis_map(replay.get("basis_row_ids"), ids)
    keys = reg.get("ordered_keys")
    ordered_ids = reg.get("ordered_canonical_ids")
    registry = build_registry(tuple(ids))
    if not isinstance(keys, list) or len(keys) != DOF_COUNT or ordered_ids != list(registry.canonical_ids):
        raise Blocked("saved semantic registry differs from the frozen 162-slot registry")
    new_row_ids = row_ids_fn(tuple(ids))
    new_rows = {(sid, family.value): row for row, (sid, family) in enumerate(new_row_ids)}
    for slot, key in enumerate(keys):
        expected_key = registry.ordered_keys[slot]
        if (not isinstance(key, dict) or key.get("slot") != slot or key.get("canonical_id") not in (None, ordered_ids[slot])
                or key.get("source_id") != expected_key.source_id or key.get("component") != expected_key.field_component.value
                or key.get("basis_kind") != expected_key.basis_kind.value or key.get("basis_index") != expected_key.basis_index):
            raise Blocked(f"saved canonical registry row is malformed at slot {slot}")
    offset = float(case["front_extension_m"])
    radius, length = float(geom_expected["radius_m"]), float(geom_expected["length_m"])
    values, gradients = _refront(np, positions, old_values, old_gradients, old_rows, new_rows, radius, offset)
    array_bytes = sum(a.nbytes for a in (positions, volume_weights, normals, area_weights, old_values, old_gradients, traction, values, gradients)) + DOF_COUNT * DOF_COUNT * 8 * 5 + 168 * 168 * 8 * 3
    if array_bytes > MAX_BYTES:
        raise Blocked("explicit numeric-array estimate exceeds 8 MiB")
    generation = hashlib.sha256()
    generation.update(json.dumps({"algorithm": "cached_phi_times_shifted_Fbar_product_rule.v1", "front_extension_m": offset, "source_result_sha256": operator_sha, "canonical_ids": ordered_ids}, sort_keys=True, separators=(",", ":")).encode())
    generation.update(np.asarray(values, dtype="<f8", order="C").tobytes())
    generation.update(np.asarray(gradients, dtype="<f8", order="C").tobytes())
    provider = _CachedProvider(np, ShapeSample, ShapeStatus, new_row_ids, positions[:VOLUME_COUNT], values[:VOLUME_COUNT], gradients[:VOLUME_COUNT])
    material = IsotropicMaterial(float(material_expected["young_modulus_Pa"]), float(material_expected["poisson_ratio"]))
    assembled = assemble(source_positions, tuple(ids), positions[:VOLUME_COUNT], volume_weights, material, provider, registry)
    if assembled.stiffness is None or assembled.status.value != "OK" or provider.index != VOLUME_COUNT:
        raise Blocked(f"shifted full-3D assembly blocked: {assembled.status.value}")
    stiffness = np.asarray(assembled.stiffness, dtype=np.float64)
    if not bool(np.isfinite(stiffness).all()):
        raise Blocked("shifted stiffness contains non-finite entries")
    force = np.zeros(DOF_COUNT, dtype=np.float64)
    component_index = {"ux": 0, "uy": 1, "uz": 2}
    for slot, key in enumerate(keys):
        component = component_index[key.get("component", key.get("field_component"))]
        row = new_rows[(key["source_id"], key["basis_kind"])]
        force[slot] = np.dot(area_weights * traction[:, component], values[VOLUME_COUNT:, row])
    if not bool(np.isfinite(force).all()):
        raise Blocked("reintegrated fixed-load vector is non-finite")
    old_rebuilt = np.zeros(DOF_COUNT, dtype=np.float64)
    for slot, key in enumerate(keys):
        component = component_index[key.get("component", key.get("field_component"))]
        row = old_rows[(key["source_id"], key["basis_kind"])]
        old_rebuilt[slot] = np.dot(area_weights * traction[:, component], old_values[VOLUME_COUNT:, row])
    load_reconstruction_error = float(np.linalg.norm(old_rebuilt - old_force) / max(np.linalg.norm(old_force), np.finfo(np.float64).tiny))
    if not math.isfinite(load_reconstruction_error) or load_reconstruction_error > 1e-10:
        raise Blocked("pinned original load does not reconstruct from its unchanged traction rows")
    volume = float(np.sum(volume_weights))
    volume_positions, boundary_positions = positions[:VOLUME_COUNT], positions[VOLUME_COUNT:]
    centroid = np.einsum("q,qc->c", volume_weights, volume_positions) / volume
    gauge_modes = np.zeros((6, VOLUME_COUNT, 3), dtype=np.float64)
    gauge_modes[0, :, 0], gauge_modes[1, :, 1], gauge_modes[2, :, 2] = 1.0, 1.0, 1.0
    axes = np.eye(3)
    for axis_id in range(3):
        gauge_modes[3 + axis_id] = np.cross(axes[axis_id], volume_positions - centroid) / radius
    gauge = np.zeros((6, DOF_COUNT), dtype=np.float64)
    for slot, key in enumerate(keys):
        component = component_index[key.get("component", key.get("field_component"))]
        row = new_rows[(key["source_id"], key["basis_kind"])]
        gauge[:, slot] = (gauge_modes[:, :, component] * volume_weights) @ values[:VOLUME_COUNT, row] / volume
    if not bool(np.isfinite(gauge).all()):
        raise Blocked("new-basis physical-volume gauge is non-finite")
    rigid = rigid_fn(source_positions, tuple(ids), registry)
    if not np.allclose(rigid, old_rigid, rtol=0.0, atol=1e-13):
        raise Blocked("saved physical rigid modes differ from the frozen source registry")
    saved_center = _array(np, src.get("rigid_mode_center_xyz_m"), (3,), "saved rigid center")
    if not np.allclose(saved_center, np.mean(source_positions, axis=0), rtol=0.0, atol=1e-13):
        raise Blocked("saved rigid center differs from the unchanged source centroid")
    rigid_at_centroid = np.array(rigid, copy=True)
    center_shift = saved_center - centroid
    for axis_id in range(3):
        rigid_at_centroid[:, 3 + axis_id] = (rigid[:, 3 + axis_id] + rigid[:, :3] @ np.cross(axes[axis_id], center_shift)) / radius
    singular = np.linalg.svd(gauge @ rigid_at_centroid, compute_uv=False)
    if not bool(np.isfinite(singular).all()) or singular[-1] <= 64.0 * np.finfo(np.float64).eps * singular[0]:
        raise Blocked("six physical-volume moments do not control all rigid modes")
    diagonal = np.diag(stiffness)
    if not bool(np.isfinite(diagonal).all()) or bool(np.any(diagonal <= 0.0)):
        raise Blocked("shifted K has a non-positive diagonal")
    scale = 1.0 / np.sqrt(diagonal)
    if not bool(np.isfinite(scale).all()):
        raise Blocked("coefficient scaling is non-finite")
    scaled_raw = scale[:, None] * stiffness * scale[None, :]
    scaled_sym = 0.5 * (scaled_raw + scaled_raw.T)
    if not bool(np.isfinite(scaled_raw).all()) or not bool(np.isfinite(scaled_sym).all()):
        raise Blocked("scaled stiffness contains non-finite entries")
    eigenvalues = np.linalg.eigvalsh(scaled_sym)
    if not bool(np.isfinite(eigenvalues).all()):
        raise Blocked("scaled stiffness eigenvalues are non-finite")
    spectral_scale = max(float(np.max(np.abs(eigenvalues))), np.finfo(np.float64).tiny)
    zero_tol = 64.0 * np.finfo(np.float64).eps * max(1.0, condition_max) * spectral_scale
    near = int(np.count_nonzero(np.abs(eigenvalues) <= zero_tol))
    positive = int(np.count_nonzero(eigenvalues > zero_tol))
    negative = int(np.count_nonzero(eigenvalues < -zero_tol))
    scaled_rigid = rigid / scale[:, None]
    rigid_residuals = np.linalg.norm(scaled_sym @ scaled_rigid, axis=0) / np.maximum(spectral_scale * np.linalg.norm(scaled_rigid, axis=0), np.finfo(np.float64).tiny)
    if not bool(np.isfinite(rigid_residuals).all()):
        raise Blocked("scaled rigid-mode residuals are non-finite")
    if near != 6 or positive != 156 or negative != 0 or np.any(rigid_residuals > 128.0 * np.finfo(np.float64).eps * max(1.0, condition_max)):
        raise Blocked("shifted K does not resolve exactly six rigid and 156 positive modes")
    row_norm = np.linalg.norm(gauge * scale[None, :], axis=1)
    if np.any(row_norm <= 0.0) or not bool(np.isfinite(row_norm).all()):
        raise Blocked("scaled gauge has a zero or invalid row")
    constraint_scale = 1.0 / row_norm
    scaled_gauge = constraint_scale[:, None] * gauge * scale[None, :]
    saddle = np.zeros((168, 168), dtype=np.float64)
    saddle[:DOF_COUNT, :DOF_COUNT] = scaled_raw
    saddle[:DOF_COUNT, DOF_COUNT:] = scaled_gauge.T
    saddle[DOF_COUNT:, :DOF_COUNT] = scaled_gauge
    rhs = np.concatenate((scale * force, constraint_scale * target))
    if not bool(np.isfinite(saddle).all()) or not bool(np.isfinite(rhs).all()):
        raise Blocked("scaled KKT matrix or right-hand side is non-finite")
    try:
        solution = np.linalg.solve(saddle, rhs)
    except np.linalg.LinAlgError as exc:
        raise Blocked(f"direct KKT solve failed ({type(exc).__name__})") from exc
    if not bool(np.isfinite(solution).all()):
        raise Blocked("direct KKT solution contains non-finite unknowns")
    coefficients = scale * solution[:DOF_COUNT]
    multipliers = constraint_scale * solution[DOF_COUNT:]
    if not bool(np.isfinite(coefficients).all()) or not bool(np.isfinite(multipliers).all()):
        raise Blocked("solution coefficients or Lagrange multipliers are non-finite")
    equilibrium = stiffness @ coefficients + gauge.T @ multipliers - force
    gauge_residual = gauge @ coefficients - target
    strain_energy = 0.5 * float(coefficients @ stiffness @ coefficients)
    external_work = float(coefficients @ force)
    multiplier_work = float(target @ multipliers)
    energy_identity_defect = 2.0 * strain_energy - external_work + multiplier_work
    scaled_kd = scale * (stiffness @ coefficients)
    scaled_clambda = scale * (gauge.T @ multipliers)
    scaled_force = scale * force
    weak_relative = float(np.linalg.norm(scale * equilibrium) / max(float(np.linalg.norm(scaled_kd) + np.linalg.norm(scaled_clambda) + np.linalg.norm(scaled_force)), np.finfo(np.float64).tiny))
    constraint_relative = float(np.linalg.norm(gauge_residual) / max(float(np.linalg.norm(target)), float(np.linalg.norm(gauge @ coefficients)), np.finfo(np.float64).tiny))
    saddle_relative = float(np.linalg.norm(saddle @ solution - rhs) / max(float(np.linalg.norm(rhs)), np.finfo(np.float64).tiny))
    residual_gates = {"weak_equilibrium_vector_finite": bool(np.isfinite(equilibrium).all()), "weak_equilibrium_relative_finite": math.isfinite(weak_relative), "physical_gauge_vector_finite": bool(np.isfinite(gauge_residual).all()), "physical_gauge_relative_finite": math.isfinite(constraint_relative), "energy_identity_defect_finite": math.isfinite(energy_identity_defect), "saddle_backward_relative_finite": math.isfinite(saddle_relative)}
    if not all(residual_gates.values()):
        raise Blocked("weak equilibrium, gauge, or energy identity residual failed the baseline finite-value gate")
    old_pi = old_energy - old_work
    new_pi = strain_energy - external_work
    delta_a = length * offset
    release = -(new_pi - old_pi) / delta_a
    boundary_force = np.einsum("q,qc->c", area_weights, traction)
    boundary_moment = np.einsum("q,qc->c", area_weights, np.cross(boundary_positions, traction))
    center_moment = np.einsum("q,qc->c", area_weights, np.cross(boundary_positions - saved_center, traction))
    rigid_work = rigid.T @ force
    wrench_values = np.concatenate((boundary_force, boundary_moment, center_moment, rigid_work))
    if not all(math.isfinite(v) for v in (old_pi, new_pi, release, strain_energy, external_work, multiplier_work, energy_identity_defect, weak_relative, constraint_relative, saddle_relative)) or not bool(np.isfinite(wrench_values).all()):
        raise Blocked("energy or finite-change diagnostic is non-finite")
    traction_hash = _sha(np.asarray(traction, dtype="<f8", order="C").tobytes())
    quadrature_hash = _sha(np.concatenate((positions.ravel(), volume_weights, normals.ravel(), area_weights)).astype("<f8").tobytes())
    return {
        "schema": "modeiii_pu_full3d_fixed_load_front_reanalysis_result.v1", "case_id": CASE_ID,
        "status": "FINITE_CHANGE_RELEASE_DIAGNOSTIC_ONLY_NOT_COMMITTED",
        "classification": "SYNTHETIC_ENGINEERING_DIAGNOSTIC_NOT_CERAMIC_NOT_AUTONOMOUS",
        "source_results": {"operator_sha256": operator_sha, "baseline_static_solve_sha256": case["baseline_static_solve"]["sha256"]},
        "front": {"old_front_x_m": 0.0, "new_front_x_m": offset, "delta_a_m": offset, "basis_generation_hash_sha256": generation.hexdigest(), "new_tip_or_branch_sampled": False},
        "fixed_data": {"source_sites": 27, "semantic_dofs": 162, "query_rows_reused": 384, "additional_queries": 0, "quadrature_and_boundary_hash_sha256": quadrature_hash, "traction_values_little_endian_f64_sha256": traction_hash, "traction_recomputed_or_rescaled": False, "base_load_reconstruction_relative_error": load_reconstruction_error},
        "assembly": {"operator": "full 3D isotropic lambda+mu B.T D B", "K_N_per_m": stiffness, "f_N": force, "K_sha256_little_endian_f64_C_order": _sha(np.asarray(stiffness, dtype="<f8", order="C").tobytes())},
        "numeric_array_bytes_estimate": array_bytes,
        "spectrum": {"near_zero": near, "positive": positive, "significant_negative": negative, "scaled_eigenvalues": eigenvalues, "rigid_scaled_residuals": rigid_residuals},
        "gauge": {"definition": "same six physical volume moments as static bridge; rotations about saved volume centroid normalized by R", "C_dimensionless": gauge, "g_fixed_from_baseline_m": target, "coefficient_pinning": False},
        "solve": {"method": "one direct diagonally scaled 168x168 KKT solve", "coefficients_m": coefficients, "lagrange_multipliers_N": multipliers, "weak_equilibrium_residual_N": equilibrium, "weak_equilibrium_residual_relative": weak_relative, "physical_gauge_residual_m": gauge_residual, "physical_gauge_residual_relative": constraint_relative, "saddle_backward_residual_relative": saddle_relative, "energy_identity_defect_2U_minus_W_plus_g_dot_lambda_J": energy_identity_defect, "residual_gate_policy": "finite-only, matching static bridge; no magnitude cutoff is frozen", "residual_gates": residual_gates, "post_solve_residual_gate_passed": all(residual_gates.values()), "strain_energy_J": strain_energy, "external_work_J": external_work, "g_dot_lambda_J": multiplier_work, "pseudoinverse_or_ridge_used": False, "load_projection_or_repair_used": False},
        "rigid_load_check": {"f_dot_r_N_or_Nm": rigid_work, "integrated_boundary_force_N": boundary_force, "integrated_moment_about_origin_Nm": boundary_moment, "integrated_moment_about_rigid_center_Nm": center_moment, "wrench_difference": rigid_work - np.concatenate((boundary_force, center_moment))},
        "potential_energy": {"definition": "Pi=U-f.T d", "baseline_Pi_J": old_pi, "extended_Pi_J": new_pi, "delta_Pi_J": new_pi - old_pi, "one_sided_crack_extension_area_m2": delta_a, "two_physical_face_area_total_m2": 2.0 * delta_a, "finite_change_release_diagnostic_J_per_m2": release, "Gc_compared": False, "growth_committed": False},
        "limitations": ["one finite fixed-quadrature comparison is not a quadrature, MLS support, or extension-size convergence study", "the finite-change release value does not establish a crack-growth direction or threshold", "no ceramic calibration, autonomous growth, or material claim"],
        "elapsed_seconds": time.monotonic() - started,
    }

def _json_bytes(value: dict) -> bytes:
    def safe(item):
        if isinstance(item, dict): return {str(k): safe(v) for k, v in item.items()}
        if isinstance(item, (list, tuple)): return [safe(v) for v in item]
        if hasattr(item, "tolist"): return safe(item.tolist())
        if hasattr(item, "item"): return safe(item.item())
        return item
    data = (json.dumps(safe(value), ensure_ascii=False, separators=(",", ":"), allow_nan=False) + "\n").encode()
    if len(data) > MAX_BYTES: raise Blocked("result exceeds 8 MiB")
    return data

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execute", action="store_true", help="perform the single pinned reanalysis")
    parser.add_argument("--output", default=str(DEFAULT_OUTPUT))
    args = parser.parse_args()
    try:
        case = _case()
        if not args.execute:
            print(json.dumps({"case_id": CASE_ID, "status": "SOURCE_ONLY_PINNED_AWAITING_DIAGNOSTIC_AUTHORIZATION", "execution_invoked": False, "numpy_imported": False}))
            return 0
        target, staging = _output_target(args.output)
        sources = _verify_sources(case)
        operator, op_sha = _pinned_json(case["operator_result"]["path"], case["operator_result"]["sha256"], "operator result")
        baseline, base_sha = _pinned_json(case["baseline_static_solve"]["path"], case["baseline_static_solve"]["sha256"], "baseline static solve")
        if baseline.get("case_id") != BASELINE_CASE_ID or baseline.get("status") != "SOLVED_DIAGNOSTIC_ONLY_UNQUALIFIED" or baseline.get("source_result", {}).get("sha256") != op_sha:
            raise Blocked("baseline static solve is not tied to the exact pinned operator result")
        modules = _numpy_modules(sources)
        result = _run(operator, op_sha, baseline, case, modules)
        result["source_results"]["baseline_static_solve_sha256"] = base_sha
        if target.exists() or staging.exists(): raise Blocked("refusing to overwrite a prior result or partial staging file")
        target.parent.mkdir(parents=True, exist_ok=True)
        staging.write_bytes(_json_bytes(result))
        os.replace(staging, target)
        print(json.dumps({"case_id": CASE_ID, "status": result["status"], "result_path": str(target), "operator_sha256": op_sha, "baseline_static_solve_sha256": base_sha}))
        return 0
    except Blocked as exc:
        print(json.dumps({"case_id": CASE_ID, "status": "BLOCKED", "reason": str(exc)}), file=sys.stderr)
        return 2

if __name__ == "__main__":
    raise SystemExit(main())
