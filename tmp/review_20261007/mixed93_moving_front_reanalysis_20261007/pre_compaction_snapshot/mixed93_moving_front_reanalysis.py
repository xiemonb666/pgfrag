"""Source-gated, two-front 93-DOF fixed-old-load LEFM reanalysis diagnostic.

This is a numerical experiment harness, not a crack-growth or material model.
It deliberately has no ghost masking, force projection, matrix repair, or retry.
The case remains execution_authorization=false until separately reviewed.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import sys
from pathlib import Path

SCHEMA = "mixed93_moving_front_reanalysis.v1"
DOF, SCALAR, MAX_BATCH = 93, 31, 256
RIGID_COUNT = 6
SPECTRUM_FACTOR = 64.0
EPS = 2.220446049250313e-16
ASYMMETRY_TOL = 1.0e-12
RIGID_RELATIVE_TOL = 128.0 * EPS * DOF
WRENCH_TOL = 1.0e-9
KKT_TOL = 1.0e-9
SIDE_CROSSCHECK_TOL = 1.0e-12
OUTPUT_LIMIT_BYTES = 12 * 1024 * 1024
MEMORY_ESTIMATE_LIMIT_BYTES = 64 * 1024 * 1024
NEW_CAP_ROWS = 1344
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
REVIEW = ROOT / "pgfrag" / "tmp" / "review_20261007" / "mixed93_moving_front_reanalysis_20261007"
CASE_PATH = HERE / "case.json"
PINS_PATH = REVIEW / "SOURCE_PINS.json"
SOURCE_PATHS = {
    "semantic_operator": "pgfrag/experiments/mixed_enriched_operator_20261006/mixed_enriched_operator.py",
    "elastic_operator": "pgfrag/experiments/efgm_elastic_operator_20261005/efgm_elastic_operator.py",
    "mls_kernel": "pgfrag/experiments/mls3d_kernel_20261005/mls3d_mls.py",
    "mls_shape": "pgfrag/experiments/mls3d_kernel_20261005/mls3d_shape.py",
    "analytic_field": "pgfrag/experiments/mixed_lefm_benchmark_20261006/mixed_lefm_field.py",
    "tip_cache_result": "pgfrag/tmp/review_20261006/tip_aligned_quadrature_study_20261006/run01/result.json",
    "tip_cache_arrays": "pgfrag/tmp/review_20261006/tip_aligned_quadrature_study_20261006/run01/arrays.npz",
    "failed_arrays": "pgfrag/tmp/review_20261006/mixed93_cached_lefm_solve_20261006/run01/arrays.npz",
    "independent_audit": "pgfrag/tmp/review_20261006/mixed93_saved_independent/audit.json",
}


class StudyError(RuntimeError):
    pass


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise StudyError(f"cannot import pinned module {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _finite(np, name, value) -> None:
    if not bool(np.isfinite(value).all()):
        raise StudyError(f"non-finite {name}")


def _pin_inputs(case: dict) -> dict[str, Path]:
    """Validate the current false/authorized case and every immutable input."""
    if not PINS_PATH.is_file():
        raise StudyError(f"missing source pins: {PINS_PATH}")
    pins = json.loads(PINS_PATH.read_text(encoding="utf-8"))
    if pins.get("schema") != "mixed93_moving_front_reanalysis_source_pins.v1":
        raise StudyError("unexpected source pin schema")
    if json.loads(CASE_PATH.read_text(encoding="utf-8")) != case:
        raise StudyError("case changed after it was read")
    case_sha = _sha(CASE_PATH)
    if pins.get("current_case_sha256") != case_sha:
        raise StudyError(f"current case SHA differs from reviewed pin: {case_sha}")
    authorized = case.get("execution_authorization", {}).get("authorized")
    if type(authorized) is not bool or pins.get("execution_authorized") is not authorized:
        raise StudyError("case and source-pin authorization states differ")
    if pins.get("entry_sha256") != _sha(Path(__file__).resolve()):
        raise StudyError("entry source differs from its reviewed pin")
    references = pins.get("references", {})
    paths = {key: ROOT / relative for key, relative in SOURCE_PATHS.items()}
    for key, path in paths.items():
        record = references.get(key)
        relative = SOURCE_PATHS[key]
        if not isinstance(record, dict) or record.get("path") != relative:
            raise StudyError(f"missing or misordered input pin: {key}")
        actual = _sha(path)
        if record.get("sha256") != actual:
            raise StudyError(f"pinned input SHA mismatch for {key}: {actual}")
    for key, record in case.get("dependency_pins", {}).items():
        path = ROOT / record["path"]
        if not path.is_file() or _sha(path) != record["sha256"]:
            raise StudyError(f"case dependency pin mismatch: {key}")
    if case.get("cached_arrays_sha256") != _sha(paths["tip_cache_arrays"]):
        raise StudyError("tip cache arrays differ from case binding")
    if case.get("cached_result_sha256") != _sha(paths["tip_cache_result"]):
        raise StudyError("tip cache result differs from case binding")
    return paths


def _validate_case(case: dict) -> None:
    if case.get("schema") != "mixed93_moving_front_reanalysis_case.v1":
        raise StudyError("unexpected case schema")
    if case.get("status") != "SOURCE_READY_EXECUTION_GATED":
        raise StudyError("case status does not describe a source-only artifact")
    if type(case.get("execution_authorization", {}).get("authorized")) is not bool:
        raise StudyError("case execution authorization must be an explicit Boolean")
    if case.get("classification") != "TOY_LEFM_FIXED_OLD_TIP_LOAD_TWO_FRONT_REANALYSIS_NOT_CZM_NOT_GROWTH_QUALIFICATION":
        raise StudyError("case classification changed")
    quad = case["quadrature"]
    caps = quad["end_caps"]
    gates = case["predeclared_gates"]
    expected = {
        "scaled_raw_asymmetry_max": ASYMMETRY_TOL,
        "spectrum_zero_negative_tolerance": "64*eps*93*max(abs(eigenvalues))",
        "spectrum_zero_eigenvalues": 6,
        "spectrum_positive_eigenvalues": 87,
        "rigid_null_relative_max": RIGID_RELATIVE_TOL,
        "rigid_virtual_work_and_wrench_relative_max": WRENCH_TOL,
        "side_traction_scaled_crosscheck_relative_max": SIDE_CROSSCHECK_TOL,
        "weak_gauge_saddle_energy_relative_max": KKT_TOL,
        "provider_new_rows_max": NEW_CAP_ROWS,
        "provider_batch_max": MAX_BATCH,
        "working_array_ceiling_bytes": MEMORY_ESTIMATE_LIMIT_BYTES,
        "combined_output_max_bytes": OUTPUT_LIMIT_BYTES,
        "single_external_timeout_seconds": 64,
    }
    for key, value in expected.items():
        if gates.get(key) != value:
            raise StudyError(f"predeclared gate {key} differs from source constant {value!r}")
    if (caps["both_caps_rows"]["coarse"], caps["both_caps_rows"]["fine"]) != (448, 896):
        raise StudyError("dual-center cap row counts differ from the source contract")
    if (caps["coarse_full_boundary_rows_including_side"], caps["fine_full_boundary_rows_including_side"]) != (704, 1152):
        raise StudyError("boundary row counts differ from the source contract")
    if quad["volume"]["per_state_rows"] != 512 or quad["side"]["rows"] != 256:
        raise StudyError("body/side quadrature contract changed")


def _cloud(np, archive, name, with_normals=False):
    xyz = np.asarray(archive[f"cloud_{name}_xyz_m"], dtype=np.float64)
    row_ids = np.asarray(archive[f"cloud_{name}_unique_row_ids"], dtype=np.int64)
    weights = np.asarray(archive[f"cloud_{name}_weights"], dtype=np.float64)
    normals = None
    if with_normals:
        normals = np.asarray(archive[f"cloud_{name}_normals"], dtype=np.float64)
        if normals.shape != xyz.shape:
            raise StudyError(f"cached normals have a bad shape: {name}")
    return xyz, row_ids, weights, normals


def _cap_cloud(np, split_a, radius, length, theta_order):
    """Shared physical-disk chart: old-tip and shifted-tip centers, seven sectors."""
    gl_s, gl_ws = np.polynomial.legendre.leggauss(4)
    gl_t, gl_wt = np.polynomial.legendre.leggauss(theta_order)
    s_nodes, s_weights = (gl_s + 1.0) * 0.5, gl_ws * 0.5
    theta0 = math.acos(split_a / (2.0 * radius))
    theta1 = math.pi - theta0
    sectors = {
        "left": ((-math.pi, -theta0), (-theta0, 0.0), (0.0, theta0), (theta0, math.pi)),
        "right": ((-math.pi, -theta1), (-theta1, theta1), (theta1, math.pi)),
    }
    centers = {"left": 0.0, "right": split_a}
    rows_xyz, rows_normals, rows_weights, region = [], [], [], []
    for zsign in (-1.0, 1.0):
        z = zsign * length * 0.5
        for chart in ("left", "right"):
            center = centers[chart]
            for sector_index, (lo, hi) in enumerate(sectors[chart]):
                theta = (lo + hi) * 0.5 + gl_t * (hi - lo) * 0.5
                theta_weights = gl_wt * (hi - lo) * 0.5
                for ti, angle in enumerate(theta):
                    cth, sth = math.cos(float(angle)), math.sin(float(angle))
                    disk = -center * cth + math.sqrt(radius * radius - center * center * sth * sth)
                    if chart == "left":
                        cut = split_a / (2.0 * cth) if cth > 0.0 else math.inf
                    else:
                        cut = -split_a / (2.0 * cth) if cth < 0.0 else math.inf
                    rmax = min(disk, cut)
                    if not math.isfinite(rmax) or rmax <= 0.0:
                        raise StudyError("cap chart has a nonpositive radial limit")
                    for ri, sval in enumerate(s_nodes):
                        rho = rmax * float(sval) ** 2
                        area = 2.0 * rmax * rmax * float(sval) ** 3 * float(s_weights[ri]) * float(theta_weights[ti])
                        rows_xyz.append((center + rho * cth, rho * sth, z))
                        rows_normals.append((0.0, 0.0, zsign))
                        rows_weights.append(area)
                        region.append((0 if zsign < 0 else 1, 0 if chart == "left" else 1, sector_index, ti, ri))
    xyz = np.asarray(rows_xyz, dtype=np.float64)
    normals = np.asarray(rows_normals, dtype=np.float64)
    weights = np.asarray(rows_weights, dtype=np.float64)
    if xyz.shape != (4 * theta_order * 7 * 2, 3):
        raise StudyError("cap point count differs from 4×theta×7 sectors×2 caps")
    _finite(np, "cap positions", xyz)
    _finite(np, "cap weights", weights)
    if bool(np.any(weights <= 0.0)) or float(np.max(np.linalg.norm(xyz[:, :2], axis=1))) > radius * (1.0 + 1.0e-12):
        raise StudyError("cap points/weights leave the physical disk")
    return xyz, normals, weights, np.asarray(region, dtype=np.int64)


def _query_cap_mls(np, mls, source_xyz, queries, config):
    values, gradients, statuses, conditions, batches = [], [], [], [], []
    for start in range(0, len(queries), MAX_BATCH):
        stop = min(start + MAX_BATCH, len(queries))
        batch = np.ascontiguousarray(queries[start:stop])
        result = mls.shape_functions_3d(source_xyz, batch, config)
        batches.append(stop - start)
        if any(getattr(status, "value", status) != "OK" for status in result.statuses):
            bad = [(i + start, getattr(status, "value", status))
                   for i, status in enumerate(result.statuses)
                   if getattr(status, "value", status) != "OK"]
            raise StudyError(f"cap MLS query blocked: {bad[:4]}")
        values.append(np.asarray(result.weights, dtype=np.float64))
        gradients.append(np.asarray(result.gradients, dtype=np.float64))
        conditions.append(np.asarray(result.condition_numbers, dtype=np.float64))
        statuses.extend(getattr(status, "value", status) for status in result.statuses)
    all_values = np.concatenate(values, axis=0)
    all_gradients = np.concatenate(gradients, axis=0)
    all_conditions = np.concatenate(conditions)
    if all_values.shape != (len(queries), 27) or all_gradients.shape != (len(queries), 27, 3):
        raise StudyError("cap MLS returned an unexpected basis shape")
    _finite(np, "cap MLS values", all_values)
    _finite(np, "cap MLS gradients", all_gradients)
    if not bool(np.isfinite(all_conditions).all()) or max(batches) > MAX_BATCH:
        raise StudyError("cap MLS condition/batch diagnostics are invalid")
    return all_values, all_gradients, statuses, all_conditions, batches


def _basis(np, xyz, base_values, base_gradients, front_x, field_module, radius):
    values = np.empty((len(xyz), SCALAR), dtype=np.float64)
    gradients = np.empty((len(xyz), SCALAR, 3), dtype=np.float64)
    values[:, :27] = base_values
    gradients[:, :27, :] = base_gradients
    for row, point in enumerate(xyz):
        x, y = float(point[0] - front_x), float(point[1])
        rho = math.hypot(x, y)
        if rho == 0.0 or (y == 0.0 and x < 0.0):
            raise StudyError(f"basis query lies at a tip or on an undeclared branch trace: row={row}")
        branch = field_module.branch_functions_local(x, y, reference_length_m=radius)
        values[row, 27:] = branch.values
        gradients[row, 27:, :] = branch.gradients_per_m
    _finite(np, "enriched basis values", values)
    _finite(np, "enriched basis gradients", gradients)
    return values, gradients


def _physical_moment_matrix(np, values, xyz, weights, radius):
    volume = float(np.sum(weights))
    if not math.isfinite(volume) or volume <= 0.0:
        raise StudyError("invalid physical volume measure")
    c = np.zeros((6, DOF), dtype=np.float64)
    for basis_index in range(SCALAR):
        phi_w = weights * values[:, basis_index] / volume
        for component in range(3):
            col = 3 * basis_index + component
            c[component, col] = float(np.sum(phi_w))
            if component == 0:
                c[4, col] = float(np.dot(phi_w, xyz[:, 2] / radius))
                c[5, col] = float(np.dot(phi_w, -xyz[:, 1] / radius))
            elif component == 1:
                c[3, col] = float(np.dot(phi_w, -xyz[:, 2] / radius))
                c[5, col] = float(np.dot(phi_w, xyz[:, 0] / radius))
            else:
                c[3, col] = float(np.dot(phi_w, xyz[:, 1] / radius))
                c[4, col] = float(np.dot(phi_w, -xyz[:, 0] / radius))
    _finite(np, "physical moment matrix", c)
    return c, volume


class _ArrayShapeProvider:
    """Strict row-aligned adapter for the existing bounded full-3D assembler."""
    def __init__(self, np, operator, xyz, values, gradients):
        self.np, self.operator = np, operator
        self.xyz, self.values, self.gradients = xyz, values, gradients
        self.index = 0
        self.support_ids = np.arange(SCALAR, dtype=np.int64)

    def evaluate(self, material_position):
        if self.index >= len(self.xyz) or not self.np.array_equal(material_position, self.xyz[self.index]):
            raise StudyError("elastic assembler query order differs from the pinned chunk")
        j = self.index
        self.index += 1
        return self.operator.ShapeSample(self.values[j], self.gradients[j], self.support_ids,
                                         self.operator.ShapeStatus.OK)


def _assemble_k(np, elastic, source_xyz, front_x, xyz, weights, values, gradients, material):
    # Four repeated tip coordinates are bookkeeping centers for the 31 scalar
    # basis rows, never additional material particles.
    tip = np.asarray((front_x, 0.0, 0.0), dtype=np.float64)
    coefficient_positions = np.ascontiguousarray(np.vstack((source_xyz, np.repeat(tip[None, :], 4, axis=0))))
    raw = np.zeros((DOF, DOF), dtype=np.float64)
    batches = []
    for start in range(0, len(xyz), MAX_BATCH):
        stop = min(start + MAX_BATCH, len(xyz))
        qxyz = np.ascontiguousarray(xyz[start:stop])
        qweights = np.ascontiguousarray(weights[start:stop])
        provider = _ArrayShapeProvider(np, elastic, qxyz, values[start:stop], gradients[start:stop])
        result = elastic.assemble_elastic_operator(coefficient_positions, qxyz, qweights, material, provider)
        if getattr(result.status, "value", result.status) != "OK" or result.stiffness is None:
            raise StudyError(f"full 3-D K assembly blocked: {result.diagnostic}")
        if provider.index != stop - start:
            raise StudyError("full 3-D assembler did not consume its complete chunk")
        raw += result.stiffness
        batches.append(stop - start)
    _finite(np, "raw full 3-D stiffness", raw)
    return raw, batches


def _spectrum_gate(np, matrix):
    diagonal = np.diag(matrix)
    _finite(np, "raw stiffness diagonal", diagonal)
    if bool(np.any(diagonal <= 0.0)):
        raise StudyError("raw stiffness has nonpositive diagonal")
    scale_vector = 1.0 / np.sqrt(diagonal)
    scaled = scale_vector[:, None] * matrix * scale_vector[None, :]
    _finite(np, "diagonally scaled raw stiffness", scaled)
    asymmetry = float(np.linalg.norm(scaled - scaled.T) / max(float(np.linalg.norm(scaled)), np.finfo(float).tiny))
    if asymmetry > ASYMMETRY_TOL:
        raise StudyError(f"scaled raw stiffness asymmetry {asymmetry:.6g} exceeds {ASYMMETRY_TOL:g}")
    eigenvalues = np.linalg.eigvalsh(scaled)
    _finite(np, "scaled stiffness eigenvalues", eigenvalues)
    spectral_scale = float(np.max(np.abs(eigenvalues)))
    tolerance = SPECTRUM_FACTOR * np.finfo(float).eps * DOF * spectral_scale
    zeros = np.abs(eigenvalues) <= tolerance
    positives = eigenvalues > tolerance
    negatives = eigenvalues < -tolerance
    metrics = {"scaled_raw_relative_asymmetry": asymmetry, "asymmetry_limit": ASYMMETRY_TOL,
               "spectral_scale": spectral_scale, "zero_negative_tolerance": tolerance,
               "zero_count": int(zeros.sum()), "positive_count": int(positives.sum()),
               "significant_negative_count": int(negatives.sum()),
               "minimum_eigenvalue": float(eigenvalues[0]), "maximum_eigenvalue": float(eigenvalues[-1]),
               "passed": int(zeros.sum()) == 6 and int(positives.sum()) == 87 and not bool(negatives.any())}
    return scaled, scale_vector, eigenvalues, metrics


def _mode_matrices(np, operator, source_ids, source_xyz, semantic_ids, radius):
    modes = operator.rigid_body_mode_coefficients(source_ids, tuple(map(tuple, source_xyz)))
    if modes.matrix_shape != (6, DOF) or modes.coefficient_ids != semantic_ids:
        raise StudyError("BASE-only six-rigid semantic layout mismatch")
    rigid = np.frombuffer(modes.values_le_f64, dtype="<f8").copy().reshape(6, DOF)
    # Normalize rotations to (omega × x)/R, matching the moment functional.
    rigid[3:, :] /= radius
    _finite(np, "normalized rigid coefficient vectors", rigid)
    return rigid, tuple(modes.mode_ids)


def _rigid_spectrum_checks(np, kscaled, scale_vector, rigid, cphys):
    scaled_modes = rigid.T / scale_vector[:, None]
    residuals = []
    for j in range(6):
        relative = float(np.linalg.norm(kscaled @ scaled_modes[:, j], ord=np.inf) / max(
            float(np.linalg.norm(kscaled, ord=np.inf) * np.linalg.norm(scaled_modes[:, j], ord=np.inf)), np.finfo(float).tiny))
        residuals.append(relative)
    coupling = cphys @ rigid.T
    singular = np.linalg.svd(coupling, compute_uv=False)
    rank = int(np.count_nonzero(singular > RIGID_RELATIVE_TOL * float(singular[0]))) if singular[0] > 0 else 0
    passed = max(residuals) <= RIGID_RELATIVE_TOL and rank == 6
    return coupling, singular, residuals, rank, passed


def _stress_tensor(np, voigt):
    tensor = np.zeros((len(voigt), 3, 3), dtype=np.float64)
    tensor[:, 0, 0], tensor[:, 1, 1], tensor[:, 2, 2] = voigt[:, 0], voigt[:, 1], voigt[:, 2]
    tensor[:, 0, 1] = tensor[:, 1, 0] = voigt[:, 3]
    tensor[:, 1, 2] = tensor[:, 2, 1] = voigt[:, 4]
    tensor[:, 0, 2] = tensor[:, 2, 0] = voigt[:, 5]
    return tensor


def _analytic_traction(np, field_module, xyz, normals, material):
    # Slightly enlarged query bounds avoid floating-point r=R/z=L/2 rejection;
    # they do not alter any LEFM stress formula.
    case = field_module.MixedLEFMCase(
        radius_m=material["R_m"] * (1.0 + 1.0e-12),
        length_m=material["L_m"] * (1.0 + 1.0e-12),
        K_I_Pa_sqrt_m=1.0, K_II_Pa_sqrt_m=0.0, K_III_Pa_sqrt_m=0.0,
        young_modulus_Pa=material["E_Pa"], poisson_ratio=material["nu"])
    stress = np.empty((len(xyz), 6), dtype=np.float64)
    for row, point in enumerate(xyz):
        result = field_module.evaluate_local_point(tuple(float(v) for v in point), case=case)
        status = getattr(result.status, "value", result.status)
        if status != "OK" or result.stress_voigt_Pa is None:
            raise StudyError(f"old-tip analytic traction query status {status} at row {row}")
        stress[row] = result.stress_voigt_Pa
    traction = np.einsum("qij,qj->qi", _stress_tensor(np, stress), normals)
    _finite(np, "old-tip analytic traction", traction)
    return traction, stress


def _b_column(component, gradient):
    gx, gy, gz = (float(v) for v in gradient)
    if component == "ux":
        return (gx, 0.0, 0.0, gy, 0.0, gz)
    if component == "uy":
        return (0.0, gy, 0.0, gx, gz, 0.0)
    return (0.0, 0.0, gz, 0.0, gy, gx)


def _shape_row(operator, source_ids, semantic_ids, values, gradients):
    components = ("ux", "uy", "uz")
    contributions = []
    for basis_index in range(SCALAR):
        scope = "SOURCE_BASE" if basis_index < 27 else "GLOBAL_FRONT"
        owner = source_ids[basis_index] if basis_index < 27 else "global"
        semantic_basis_index = 0 if basis_index < 27 else basis_index - 26
        for ci, component in enumerate(components):
            slot = 3 * basis_index + ci
            contributions.append(operator.SemanticContribution(
                semantic_ids[slot], scope, owner, component, semantic_basis_index,
                float(values[basis_index]), tuple(float(v) for v in gradients[basis_index]),
                _b_column(component, gradients[basis_index])))
    row_ids = tuple(item.semantic_id for item in contributions)
    if row_ids != semantic_ids:
        raise StudyError("generated force shape row differs from pinned semantic registry")
    return operator.SemanticShapeRow("global_front", source_ids, semantic_ids, tuple(contributions), "OK")


def _assemble_force(np, operator, source_ids, semantic_ids, xyz, normals, weights, traction, values, gradients):
    force = np.zeros(DOF, dtype=np.float64)
    batches = []
    for start in range(0, len(xyz), MAX_BATCH):
        stop = min(start + MAX_BATCH, len(xyz))
        points = []
        for row in range(start, stop):
            shape = _shape_row(operator, source_ids, semantic_ids, values[row], gradients[row])
            points.append(operator.SurfaceTractionPoint(
                float(weights[row]), tuple(float(v) for v in traction[row]), shape))
        raw = operator.assemble_raw_force(source_ids, (tuple(points),), layout=operator.GLOBAL_FRONT_LAYOUT)
        if raw.dof_count != DOF or raw.quadrature_count != stop - start:
            raise StudyError("public semantic force assembler returned an unexpected layout/count")
        force += np.frombuffer(raw.values_le_f64, dtype="<f8")
        batches.append(stop - start)
    _finite(np, "raw generalized boundary force", force)
    return force, batches


def _wrench(np, xyz, normals, weights, traction, rigid, force, radius):
    total_force = np.einsum("q,qi->i", weights, traction)
    moment = np.einsum("q,qi->i", weights, np.cross(xyz, traction))
    resultant = float(np.sum(weights * np.linalg.norm(traction, axis=1)))
    virtual = np.zeros(6, dtype=np.float64) if force is None else rigid @ force
    reference = np.concatenate((total_force, moment / radius))
    defect = None if force is None else float(np.max(np.abs(virtual - reference)))
    numerator = max(float(np.linalg.norm(total_force)), float(np.linalg.norm(moment)) / radius,
                    0.0 if force is None else float(np.linalg.norm(virtual, ord=np.inf)))
    relative = numerator / max(resultant, np.finfo(float).tiny)
    return {"total_force_N": total_force.tolist(), "moment_about_origin_Nm": moment.tolist(),
            "rigid_virtual_work_N": virtual.tolist(), "virtual_work_vs_wrench_max_abs_N": defect,
            "integral_traction_norm_N": resultant, "normalized_wrench_virtual_work": relative,
            "passed": relative <= WRENCH_TOL}


def _solve_kkt(np, kraw, kscaled, cphys, force, scale_vector):
    cbar = cphys * scale_vector[None, :]
    row_norm = np.linalg.norm(cbar, axis=1)
    if bool(np.any(~np.isfinite(row_norm))) or bool(np.any(row_norm <= np.finfo(float).tiny)):
        raise StudyError("physical gauge has an empty or nonfinite scaled row")
    d = 1.0 / row_norm
    chat = d[:, None] * cbar
    saddle = np.zeros((DOF + 6, DOF + 6), dtype=np.float64)
    saddle[:DOF, :DOF], saddle[:DOF, DOF:], saddle[DOF:, :DOF] = kscaled, chat.T, chat
    rhs = np.concatenate((scale_vector * force, np.zeros(6, dtype=np.float64)))
    solution = np.linalg.solve(saddle, rhs)
    _finite(np, "scaled KKT solution", solution)
    y, lhat = solution[:DOF], solution[DOF:]
    q, lam = scale_vector * y, d * lhat
    _finite(np, "physical solution", q)
    _finite(np, "physical multipliers", lam)
    scaled_residual = saddle @ solution - rhs
    saddle_rel = float(np.linalg.norm(scaled_residual, ord=np.inf) / max(
        float(np.linalg.norm(saddle, ord=np.inf) * np.linalg.norm(solution, ord=np.inf)
              + np.linalg.norm(rhs, ord=np.inf)), np.finfo(float).tiny))
    weak = kraw @ q + cphys.T @ lam - force
    weak_rel = float(np.linalg.norm(weak, ord=np.inf) / max(
        float(np.linalg.norm(kraw, ord=np.inf) * np.linalg.norm(q, ord=np.inf)
              + np.linalg.norm(cphys.T @ lam, ord=np.inf) + np.linalg.norm(force, ord=np.inf)),
        np.finfo(float).tiny))
    gauge = cphys @ q
    gauge_rel = float(np.linalg.norm(gauge, ord=np.inf) / max(
        float(np.linalg.norm(cphys, ord=np.inf) * np.linalg.norm(q, ord=np.inf)), np.finfo(float).tiny))
    twice_u = float(q @ (kraw @ q))
    work = float(q @ force)
    gauge_work = 0.0
    energy_rel = abs(twice_u - work + gauge_work) / max(abs(twice_u), abs(work), abs(gauge_work), np.finfo(float).tiny)
    metrics = {"scaled_saddle_backward_error": saddle_rel, "physical_weak_relative_error": weak_rel,
               "physical_gauge_relative_error": gauge_rel, "energy_identity_relative_error": energy_rel,
               "twice_strain_energy_J": twice_u, "strain_energy_J": 0.5 * twice_u,
               "external_work_J": work, "gauge_multiplier_work_J": gauge_work,
               "potential_U_minus_W_J": 0.5 * twice_u - work,
               "passed": max(saddle_rel, weak_rel, gauge_rel, energy_rel) <= KKT_TOL}
    arrays = {"Cbar_scaled": cbar, "gauge_row_scale": d, "Chat_scaled": chat,
              "saddle_scaled": saddle, "rhs_scaled": rhs, "solution_scaled": solution,
              "lambda_scaled": lhat}
    return q, lam, metrics, arrays


def _run(case: dict, paths: dict[str, Path]):
    import numpy as np

    semantic = _module("moving_front_semantic_operator", ROOT / SOURCE_PATHS["semantic_operator"])
    elastic = _module("moving_front_elastic_operator", ROOT / SOURCE_PATHS["elastic_operator"])
    _module("mls3d_mls", ROOT / SOURCE_PATHS["mls_kernel"])
    mls = _module("mls3d_shape", ROOT / SOURCE_PATHS["mls_shape"])
    field = _module("moving_front_analytic_field", ROOT / SOURCE_PATHS["analytic_field"])
    result = json.loads(paths["tip_cache_result"].read_text(encoding="utf-8"))
    if (result.get("schema"), result.get("status"), result.get("assembly_status"),
            result.get("mechanical_eligibility"), result.get("solve_allowed")) != (
            "tip_aligned_full3d_quadrature_study_result.v1", "OK_DIAGNOSTIC_ONLY_UNQUALIFIED",
            "OK", "SPECTRUM_GATE_PASSED", True):
        raise StudyError("tip-aligned arrays do not have the pinned diagnostic cache status")
    source_ids = tuple(result["registry"]["source_ids"])
    if len(source_ids) != 27 or source_ids != tuple(sorted(source_ids)):
        raise StudyError("cached source ID order is not canonical")
    semantic_ids = semantic.semantic_coefficient_ids(source_ids, semantic.GLOBAL_FRONT_LAYOUT)
    if len(semantic_ids) != DOF:
        raise StudyError("public semantic registry is not 93 DOFs")
    geometry, material0 = case["geometry"], case["material"]
    radius, length = float(geometry["radius_m"]), float(geometry["length_m"])
    material = {**material0, "R_m": radius, "L_m": length}
    source_sha = {}
    for key, path in paths.items():
        source_sha[key] = _sha(path)

    with np.load(paths["tip_cache_arrays"], allow_pickle=False) as cache, np.load(paths["failed_arrays"], allow_pickle=False) as failed:
        source_xyz = np.asarray(cache["source_positions_xyz_m"], dtype=np.float64)
        unique_xyz = np.asarray(cache["unique_query_positions_xyz_m"], dtype=np.float64)
        base_values_all = np.asarray(cache["base_mls_values"], dtype=np.float64)
        base_gradients_all = np.asarray(cache["base_mls_gradients_per_m"], dtype=np.float64)
        if (source_xyz.shape, unique_xyz.shape, base_values_all.shape, base_gradients_all.shape) != (
                (27, 3), (4224, 3), (4224, 27), (4224, 27, 3)):
            raise StudyError("pinned cache table shapes changed")
        _finite(np, "source positions", source_xyz)
        _finite(np, "cache query positions", unique_xyz)
        _finite(np, "cached MLS values", base_values_all)
        _finite(np, "cached MLS gradients", base_gradients_all)
        clouds = {
            "old_volume": _cloud(np, cache, "a0_volume32"),
            "new_volume": _cloud(np, cache, "a0_0005_volume32"),
            "side": _cloud(np, cache, "side", True),
        }
        expected = {"old_volume": 512, "new_volume": 512, "side": 256}
        base_cloud = {}
        for name, (xyz, row_ids, weights, normals) in clouds.items():
            if xyz.shape != (expected[name], 3) or row_ids.shape != (len(xyz),) or weights.shape != (len(xyz),):
                raise StudyError(f"pinned cache cloud count/shape changed: {name}")
            if int(row_ids.min()) < 0 or int(row_ids.max()) >= len(unique_xyz) or not np.array_equal(xyz, unique_xyz[row_ids]):
                raise StudyError(f"cached cloud row IDs do not resolve to the pinned MLS table: {name}")
            _finite(np, f"{name} weights", weights)
            if bool(np.any(weights <= 0.0)):
                raise StudyError(f"nonpositive {name} quadrature weight")
            if normals is not None and (float(np.max(np.abs(np.linalg.norm(normals, axis=1) - 1.0))) > 1.0e-10):
                raise StudyError("cached side normals are not unit vectors")
            base_cloud[name] = (base_values_all[row_ids], base_gradients_all[row_ids])
        if not np.array_equal(source_xyz, np.asarray(failed["source_positions_xyz_m"], dtype=np.float64)):
            raise StudyError("failed-cache source positions differ from the frozen tip cache")

        split_a = float(case["front_states"][1]["front_x_m"])
        caps = {}
        for label, order in (("coarse", 8), ("fine", 16)):
            xyz, normals, weights, region = _cap_cloud(np, split_a, radius, length, order)
            caps[label] = {"xyz": xyz, "normals": normals, "weights": weights, "region": region}
        coarse_count, fine_count = len(caps["coarse"]["xyz"]), len(caps["fine"]["xyz"])
        cap_xyz = np.concatenate((caps["coarse"]["xyz"], caps["fine"]["xyz"]))
        if (coarse_count, fine_count, len(cap_xyz)) != (448, 896, NEW_CAP_ROWS):
            raise StudyError("new cap query count differs from the reviewed six-batch budget")
        if len(np.unique(cap_xyz, axis=0)) != NEW_CAP_ROWS:
            raise StudyError("cap query positions are not unique across coarse/fine rules")
        provider_cfg = case["provider"]
        config = mls.MLS3DConfig(reference_origin=tuple(provider_cfg["reference_origin_m"]),
            reference_scale=float(provider_cfg["reference_scale_m"]), support_radius=float(provider_cfg["support_radius_m"]),
            max_condition_number=float(provider_cfg["max_condition_number"]),
            rank_relative_tolerance=float(provider_cfg["rank_relative_tolerance"]),
            switch_tolerance=float(provider_cfg["switch_tolerance"]),
            duplicate_tolerance=float(provider_cfg["duplicate_tolerance"]))
        cap_base_values, cap_base_gradients, cap_status, cap_conditions, query_batches = _query_cap_mls(
            np, mls, source_xyz, cap_xyz, config)
        if query_batches != [256, 192, 256, 256, 256, 128]:
            raise StudyError(f"actual cap query batch order differs from declaration: {query_batches}")

        split = {"coarse": slice(0, coarse_count), "fine": slice(coarse_count, NEW_CAP_ROWS)}
        cap_base = {name: (cap_base_values[s], cap_base_gradients[s]) for name, s in split.items()}
        side_xyz, side_ids, side_weights, side_normals = clouds["side"]
        failed_xyz = np.asarray(failed["boundary_xyz_m"], dtype=np.float64)
        failed_normals = np.asarray(failed["boundary_normals"], dtype=np.float64)
        failed_weights = np.asarray(failed["boundary_weights_m2"], dtype=np.float64)
        saved_traction_100 = np.asarray(failed["pure_I_boundary_traction_Pa"], dtype=np.float64)
        if (failed_xyz.shape, failed_normals.shape, failed_weights.shape, saved_traction_100.shape) != (
                (768, 3), (768, 3), (768,), (768, 3)):
            raise StudyError("failed arrays lack the pinned side-load recovery rows")
        if not (np.array_equal(failed_xyz[:256], side_xyz) and np.array_equal(failed_normals[:256], side_normals)
                and np.array_equal(failed_weights[:256], side_weights)):
            raise StudyError("saved side load point/normal/weight binding differs from frozen side cloud")
        side_saved = saved_traction_100[:256] * float(case["frozen_load"]["scale"])
        side_analytic, _ = _analytic_traction(np, field, side_xyz, side_normals, material)
        side_check = float(np.linalg.norm(side_saved - side_analytic) / max(
            float(np.linalg.norm(side_analytic)), np.finfo(float).tiny))
        if side_check > float(case["frozen_load"]["saved_side_traction_scaled_relative_check_max"]):
            raise StudyError(f"scaled saved side traction disagrees with the old-tip KI=1 field: {side_check:.6g}")
        cap_tractions = {}
        for label, group in caps.items():
            cap_tractions[label], _ = _analytic_traction(np, field, group["xyz"], group["normals"], material)

        states = []
        for state_name, front in (("old_front", 0.0), ("candidate_front", split_a)):
            volume_key = "old_volume" if front == 0.0 else "new_volume"
            vol_xyz, vol_ids, vol_weights, _ = clouds[volume_key]
            side_values, side_gradients = _basis(np, side_xyz, *base_cloud["side"], front, field, radius)
            vol_values, vol_gradients = _basis(np, vol_xyz, *base_cloud[volume_key], front, field, radius)
            cap_basis = {}
            for label, cap in caps.items():
                cap_basis[label] = _basis(np, cap["xyz"], *cap_base[label], front, field, radius)
            cphys, volume_measure = _physical_moment_matrix(np, vol_values, vol_xyz, vol_weights, radius)
            if abs(volume_measure - math.pi * radius * radius * length) > 1.0e-10 * math.pi * radius * radius * length:
                raise StudyError(f"{state_name} front-aligned volume rule misses cylinder measure")
            kraw, kbatches = _assemble_k(np, elastic, source_xyz, front, vol_xyz, vol_weights,
                                          vol_values, vol_gradients,
                                          elastic.IsotropicMaterial(material["E_Pa"], material["nu"]))
            kscaled, scale_vector, eigenvalues, spectrum = _spectrum_gate(np, kraw)
            rigid, rigid_ids = _mode_matrices(np, semantic, source_ids, source_xyz, semantic_ids, radius)
            coupling, coupling_singular, null_residuals, coupling_rank, rigid_ok = _rigid_spectrum_checks(
                np, kscaled, scale_vector, rigid, cphys)
            states.append({"name": state_name, "front_x_m": front, "xyz": vol_xyz, "weights": vol_weights,
                "values": vol_values, "gradients": vol_gradients, "C": cphys, "volume_m3": volume_measure,
                "K": kraw, "Kscaled": kscaled, "scale": scale_vector, "eigenvalues": eigenvalues,
                "spectrum": spectrum, "rigid": rigid, "rigid_ids": rigid_ids,
                "coupling": coupling, "coupling_singular": coupling_singular,
                "null_residuals": null_residuals, "coupling_rank": coupling_rank,
                "mechanical_gate": bool(spectrum["passed"] and rigid_ok), "K_batches": kbatches,
                "side_values": side_values, "side_gradients": side_gradients, "cap_basis": cap_basis})

        if not all(state["mechanical_gate"] for state in states):
            raise StudyError("one or both raw scaled spectrum/rigid gates failed; no KKT solve attempted")

        audit = json.loads(paths["independent_audit"].read_text(encoding="utf-8"))
        if audit.get("schema") != "mixed93_saved_independent_audit.v1" or audit.get("status") != "RECOVERED_CACHE_INDEPENDENTLY_CONSISTENT_PHYSICS_UNQUALIFIED":
            raise StudyError("failed-producer cache audit no longer matches its declared source status")
        if audit.get("consistency_gates") != {"matrix_and_layout": True, "spectrum_and_rigid_modes": True,
                                                "four_case_algebra_and_fields": True}:
            raise StudyError("independent failed-cache consistency gates changed")
        audit_path = Path(audit["input_paths"]["failed_arrays_npz"]).resolve()
        if audit_path != paths["failed_arrays"].resolve():
            raise StudyError("independent audit does not bind the actual failed arrays path")
        if audit.get("input_bytes", {}).get("failed_arrays_npz") != paths["failed_arrays"].stat().st_size:
            raise StudyError("independent audit failed-arrays byte count changed")
        if audit.get("pins_pre", {}).get("failed_arrays") != _sha(paths["failed_arrays"]):
            raise StudyError("independent audit failed-arrays pre-pin changed")
        if audit.get("pins_post", {}).get("failed_arrays") != _sha(paths["failed_arrays"]):
            raise StudyError("independent audit failed-arrays post-pin changed")
        if audit.get("producer_run", {}).get("rc") != 1 or audit.get("producer_run", {}).get("result_json_present") is not False:
            raise StudyError("prior producer history no longer records rc=1 and no result JSON")

        ids_array = np.asarray(semantic_ids)
        rigid = states[0]["rigid"]
        boundary_groups = {}
        for label, cap in caps.items():
            boundary_groups[label] = {
                "xyz": np.concatenate((side_xyz, cap["xyz"])),
                "normals": np.concatenate((side_normals, cap["normals"])),
                "weights": np.concatenate((side_weights, cap["weights"])),
                "traction": np.concatenate((side_saved, cap_tractions[label])),
                "side_rows": len(side_xyz),
            }
        wrench_metrics = {label: _wrench(np, g["xyz"], g["normals"], g["weights"], g["traction"],
                                          rigid, None, radius)
                          for label, g in boundary_groups.items()}
        if not wrench_metrics["fine"]["passed"]:
            raise StudyError(f"fine fixed-load rigid wrench/virtual-work gate failed: {wrench_metrics['fine']}")

        arrays = {
            "source_positions_xyz_m": source_xyz, "source_ids": np.asarray(source_ids), "semantic_ids": ids_array,
            "rigid_mode_coefficients_normalized": rigid,
            "side_xyz_m": side_xyz, "side_normals": side_normals, "side_area_weights_m2": side_weights,
            "side_saved_traction_KI100_Pa": saved_traction_100[:256],
            "side_frozen_traction_KI1_Pa": side_saved,
            "cap_query_xyz_m": cap_xyz, "cap_query_base_mls_values": cap_base_values,
            "cap_query_base_mls_gradients_per_m": cap_base_gradients,
            "cap_query_condition_numbers": cap_conditions,
            "cap_query_region_index": np.concatenate((caps["coarse"]["region"], caps["fine"]["region"])),
            "constitutive_matrix_Pa": np.asarray(semantic.isotropic_constitutive_matrix_pa(material["E_Pa"], material["nu"]), dtype=np.float64),
        }
        metrics = {
            "schema": SCHEMA, "case_id": case["case_id"], "status": "ASSEMBLY_DIAGNOSTIC_PENDING_GATES",
            "classification": case["classification"], "execution_authorized": True,
            "physics_qualification": "UNKNOWN_NOT_CRACK_GROWTH_OR_MATERIAL_QUALIFICATION",
            "source_sha256": {**source_sha, "entry": _sha(Path(__file__).resolve()),
                "case": _sha(CASE_PATH), "source_pins": _sha(PINS_PATH)},
            "source_ids": list(source_ids), "semantic_ids": list(semantic_ids), "semantic_dof_count": DOF,
            "cached_rows_reused": 1280, "new_cap_mls_rows": NEW_CAP_ROWS,
            "new_cap_query_batches": query_batches, "new_cap_mls_statuses": sorted(set(cap_status)),
            "new_cap_mls_condition_max": float(np.max(cap_conditions)),
            "side_saved_traction_relative_crosscheck": side_check,
            "side_saved_traction_scale_applied_once": float(case["frozen_load"]["scale"]),
            "cap_rule_rows": {"coarse": coarse_count, "fine": fine_count},
            "boundary_rows": {k: len(v["xyz"]) for k, v in boundary_groups.items()},
            "fixed_physical_load_wrench": wrench_metrics,
            "states": [], "energy_release": {"definition": "-(Pi_candidate-Pi_old)/(L*delta_a)",
                "classification": "finite fixed-load two-state diagnostic only; no Gc/CZM/growth inference"},
            "limitations": ["finite fixed quadrature; no convergence proof", "a>0 candidate solution is unknown under the frozen old-tip load",
                "near-tip LEFM is a small-strain idealization", "failed source arrays have rc=1 and are independently recovered, not a successful producer receipt"],
        }
        for state in states:
            name, front = state["name"], state["front_x_m"]
            arrays.update({f"{name}_volume_xyz_m": state["xyz"], f"{name}_volume_weights_m3": state["weights"],
                f"{name}_volume_values": state["values"], f"{name}_volume_gradients_per_m": state["gradients"],
                f"{name}_K_raw_N_per_m": state["K"], f"{name}_K_scaled_raw": state["Kscaled"],
                f"{name}_K_scale_sqrt_m_per_sqrt_N": state["scale"], f"{name}_eigenvalues_scaled": state["eigenvalues"],
                f"{name}_C_physical": state["C"], f"{name}_rigid_gauge_coupling": state["coupling"],
                f"{name}_rigid_coupling_singular_values": state["coupling_singular"]})
            state_metrics = {"front_x_m": front, "raw_K_units": "N/m", "gauge_units": "m",
                "volume_measure_m3": state["volume_m3"], "K_assembly_batches": state["K_batches"],
                "spectrum": state["spectrum"], "rigid_null_relative_residuals": state["null_residuals"],
                "rigid_null_limit": RIGID_RELATIVE_TOL, "physical_gauge_rigid_rank": state["coupling_rank"],
                "rigid_mode_ids": list(state["rigid_ids"]),
                "mechanical_gate_passed": state["mechanical_gate"], "loads": {}}
            for label, group in boundary_groups.items():
                cap = caps[label]
                cap_values, cap_gradients = state["cap_basis"][label]
                values = np.concatenate((state["side_values"], cap_values))
                gradients = np.concatenate((state["side_gradients"], cap_gradients))
                force, force_batches = _assemble_force(np, semantic, source_ids, semantic_ids,
                    group["xyz"], group["normals"], group["weights"], group["traction"], values, gradients)
                label_wrench = _wrench(np, group["xyz"], group["normals"], group["weights"],
                                       group["traction"], state["rigid"], force, radius)
                state_metrics["loads"][label] = {"generalized_load_N": force.tolist(),
                    "force_assembly_batches": force_batches, "rigid_wrench_virtual_work": label_wrench}
                arrays.update({f"{name}_{label}_boundary_xyz_m": group["xyz"],
                    f"{name}_{label}_boundary_normals": group["normals"],
                    f"{name}_{label}_boundary_area_weights_m2": group["weights"],
                    f"{name}_{label}_boundary_traction_Pa": group["traction"],
                    f"{name}_side_values": state["side_values"], f"{name}_side_gradients_per_m": state["side_gradients"],
                    f"{name}_{label}_cap_values": cap_values, f"{name}_{label}_cap_gradients_per_m": cap_gradients,
                    f"{name}_{label}_boundary_load_N": force})
                if label == "fine":
                    if not label_wrench["passed"]:
                        raise StudyError(f"{name} fine generalized load violates rigid wrench gate")
                    q, lam, solve_metrics, solve_arrays = _solve_kkt(
                        np, state["K"], state["Kscaled"], state["C"], force, state["scale"])
                    if not state["mechanical_gate"] or not solve_metrics["passed"]:
                        raise StudyError(f"{name} solve gate blocked: {solve_metrics}")
                    state_metrics["solve"] = {**solve_metrics, "q_solution": q.tolist(),
                        "lagrange_multipliers_physical": lam.tolist(), "constraint_target_m": [0.0] * 6}
                    for key, value in solve_arrays.items():
                        arrays[f"{name}_{key}"] = value
                    arrays[f"{name}_q_solution"] = q
                    arrays[f"{name}_lambda_physical"] = lam
                    u = np.zeros((len(state["xyz"]), 3), dtype=np.float64)
                    for j in range(SCALAR):
                        u += state["values"][:, j, None] * q[3*j:3*j+3][None, :]
                    bmat = _strain_displacement(state["gradients"])
                    strain = np.einsum("qai,i->qa", bmat, q)
                    dmat = np.asarray(semantic.isotropic_constitutive_matrix_pa(material["E_Pa"], material["nu"]), dtype=np.float64)
                    stress = strain @ dmat.T
                    arrays.update({f"{name}_volume_solution_displacement_m": u,
                        f"{name}_volume_solution_strain_m_per_m": strain,
                        f"{name}_volume_solution_stress_Pa": stress})
            metrics["states"].append(state_metrics)
        old_pi = metrics["states"][0]["solve"]["potential_U_minus_W_J"]
        new_pi = metrics["states"][1]["solve"]["potential_U_minus_W_J"]
        delta_a = split_a
        metrics["energy_release"].update({"Pi_old_J": old_pi, "Pi_candidate_J": new_pi,
            "G_fd_J_per_m2": -(new_pi - old_pi) / (length * delta_a),
            "one_sided_projected_area_m2": length * delta_a,
            "two_face_area_m2": 2.0 * length * delta_a,
            "path_and_quadrature_sensitivity": "UNKNOWN"})
        metrics["status"] = "OK_FIXED_LOAD_REANALYSIS_DIAGNOSTIC_UNQUALIFIED"
        metrics["assembly_status"] = "OK"
        metrics["mechanical_eligibility"] = "TWO_STATE_SPECTRUM_WRENCH_GAUGE_AND_SOLVE_GATES_PASSED"
        metrics["solve_allowed"] = True
        estimated = sum(np.asarray(value).nbytes for value in arrays.values()) + 2 * MAX_BATCH * 6 * DOF * 8
        metrics["explicit_array_memory_estimate_bytes"] = int(estimated)
        metrics["explicit_array_memory_estimate_is_rss"] = False
        if estimated > MEMORY_ESTIMATE_LIMIT_BYTES:
            raise StudyError(f"explicit working array estimate exceeds {MEMORY_ESTIMATE_LIMIT_BYTES} bytes")
        return metrics, arrays, 0


def _strain_displacement(gradients):
    import numpy as np
    b = np.zeros((len(gradients), 6, DOF), dtype=np.float64)
    for basis_index in range(SCALAR):
        gx, gy, gz = (gradients[:, basis_index, axis] for axis in range(3))
        for component in range(3):
            col = 3 * basis_index + component
            if component == 0:
                b[:, 0, col], b[:, 3, col], b[:, 5, col] = gx, gy, gz
            elif component == 1:
                b[:, 1, col], b[:, 3, col], b[:, 4, col] = gy, gx, gz
            else:
                b[:, 2, col], b[:, 4, col], b[:, 5, col] = gz, gy, gx
    return b


def _write_result(output: Path, metrics: dict, arrays: dict, return_code: int) -> None:
    import numpy as np
    def json_default(value):
        if isinstance(value, np.generic):
            return value.item()
        raise TypeError(f"unsupported JSON value: {type(value).__name__}")
    output.parent.mkdir(parents=True, exist_ok=False)
    array_path = output.with_name("arrays.npz")
    np.savez_compressed(array_path, **arrays)
    metrics["arrays_npz"] = {"path": array_path.name, "sha256": _sha(array_path), "bytes": array_path.stat().st_size}
    metrics["return_code"] = return_code
    output.write_text(json.dumps(metrics, indent=2, sort_keys=True, allow_nan=False, default=json_default) + "\n", encoding="utf-8")
    total = array_path.stat().st_size + output.stat().st_size
    if total > OUTPUT_LIMIT_BYTES:
        raise StudyError(f"combined output exceeds {OUTPUT_LIMIT_BYTES} bytes")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--execute", action="store_true", help="requires separate case and source-pin authorization")
    parser.add_argument("--output", type=Path, default=REVIEW / "run01" / "result.json")
    args = parser.parse_args(argv)
    if not args.execute:
        print("SOURCE_READY_EXECUTION_DISABLED")
        return 0
    case = json.loads(CASE_PATH.read_text(encoding="utf-8"))
    _validate_case(case)
    output = args.output.resolve()
    try:
        output.relative_to(REVIEW.resolve())
    except ValueError as exc:
        raise StudyError("output must remain within this experiment's review directory") from exc
    if output.exists() or output.parent.exists():
        raise StudyError("refusing an existing output or run directory")
    paths = _pin_inputs(case)
    if case["execution_authorization"].get("authorized") is not True:
        raise StudyError("case has no explicit numerical authorization")
    metrics, arrays, code = _run(case, paths)
    _write_result(output, metrics, arrays, code)
    print(json.dumps({"status": metrics["status"], "result": str(output), "return_code": code}, sort_keys=True))
    return code


if __name__ == "__main__":
    raise SystemExit(main())
