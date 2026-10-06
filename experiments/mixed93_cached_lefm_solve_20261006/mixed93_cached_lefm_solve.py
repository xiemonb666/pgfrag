"""Cached 93-DOF mixed-mode LEFM benchmark; execution is explicitly gated."""
from __future__ import annotations
import argparse, hashlib, importlib.util, json, math, sys
from pathlib import Path
SCHEMA = "mixed93_cached_lefm_solve.v1"
DOF, SCALAR, MAX_BATCH = 93, 31, 256
RIGID_COUNT, SPECTRUM_FACTOR = 6, 64.0
KKT_TOL, ENERGY_TOL, ASYMMETRY_TOL = 1.0e-9, 1.0e-9, 1.0e-12
FIELD_RELATIVE_LIMIT, OUTPUT_LIMIT_BYTES = 0.02, 12 * 1024 * 1024
MEMORY_ESTIMATE_LIMIT_BYTES = 64 * 1024 * 1024
RIGID_RELATIVE_TOL = 128.0 * 2.220446049250313e-16 * DOF
REPRESENTATION_RELATIVE_TOL, REPRESENTATION_EPS_FACTOR = 1.0e-11, 64.0
HERE = Path(__file__).resolve().parent
ROOT, REVIEW = HERE.parents[2], HERE.parents[2] / "pgfrag" / "tmp" / "review_20261006" / "mixed93_cached_lefm_solve_20261006"
CASE_PATH, PINS_PATH = HERE / "case.json", REVIEW / "SOURCE_PINS.json"
SOURCE_PATHS = {
    "mixed_lefm_field": "pgfrag/experiments/mixed_lefm_benchmark_20261006/mixed_lefm_field.py",
    "mixed_enrichment_algebra": "pgfrag/experiments/mixed_lefm_enrichment_algebra_20261006/mixed_enrichment_algebra.py",
    "global_semantic_operator": "pgfrag/experiments/mixed_enriched_operator_20261006/mixed_enriched_operator.py",
    "cached_result": "pgfrag/tmp/review_20261006/tip_aligned_quadrature_study_20261006/run01/result.json",
    "cached_arrays": "pgfrag/tmp/review_20261006/tip_aligned_quadrature_study_20261006/run01/arrays.npz",
}
class StudyError(RuntimeError): pass
def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""): digest.update(block)
    return digest.hexdigest()
def _module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise StudyError(f"cannot import pinned source {path}")
    module = importlib.util.module_from_spec(spec); sys.modules[name] = module
    spec.loader.exec_module(module)
    return module
def _finite(np, name: str, array) -> None:
    if not np.isfinite(array).all():
        raise StudyError(f"non-finite values in {name}")

def _pin_inputs(case: dict) -> tuple[Path, Path]:
    if not PINS_PATH.is_file(): raise StudyError(f"missing source pins: {PINS_PATH}")
    pins = json.loads(PINS_PATH.read_text(encoding="utf-8"))
    if pins.get("schema") != "mixed93_source_pins.v1": raise StudyError("unexpected SOURCE_PINS schema")
    if json.loads(CASE_PATH.read_text(encoding="utf-8")) != case: raise StudyError("case object changed after it was read")
    case_sha = _sha(CASE_PATH)
    if pins.get("current_case_sha256") != case_sha:
        raise StudyError(f"current case SHA-256 differs from reviewed pin: {case_sha}")
    if pins.get("execution_authorized") is not case.get("execution_authorization", {}).get("authorized"):
        raise StudyError("SOURCE_PINS authorization state differs from the pinned case")
    _validate_case_contract(case)
    expected = pins.get("references", {})
    for key, relative in SOURCE_PATHS.items():
        record = expected.get(key)
        path = ROOT / Path(relative)
        if not isinstance(record, dict) or record.get("path") != relative:
            raise StudyError(f"missing or misordered pin for {key}")
        actual = _sha(path)
        if record.get("sha256") != actual: raise StudyError(f"pinned SHA-256 mismatch for {key}: {actual}")
    cache_result, cache_arrays = (ROOT / SOURCE_PATHS[k] for k in ("cached_result", "cached_arrays"))
    if case.get("cached_result_sha256") != _sha(cache_result): raise StudyError("case cache-result SHA mismatch")
    if case.get("cached_arrays_sha256") != _sha(cache_arrays): raise StudyError("case arrays SHA mismatch")
    if pins.get("entry_sha256") != _sha(Path(__file__).resolve()): raise StudyError("entry source differs from pin")
    return cache_result, cache_arrays

def _validate_case_contract(case: dict) -> None:
    gates = case.get("predeclared_gates", {})
    expected = {
        "scaled_raw_asymmetry_max": ASYMMETRY_TOL, "rigid_null_relative_error_max": RIGID_RELATIVE_TOL,
        "kkt_backward_error_max": KKT_TOL, "gauge_relative_error_max": KKT_TOL,
        "energy_identity_relative_error_max": ENERGY_TOL, "engineering_field_relative_error_max": FIELD_RELATIVE_LIMIT,
        "output_combined_limit_bytes": OUTPUT_LIMIT_BYTES, "explicit_working_memory_estimate_limit_bytes": MEMORY_ESTIMATE_LIMIT_BYTES,
        "exact_representation_relative_error_max": REPRESENTATION_RELATIVE_TOL,
        "exact_representation_absolute_epsilon_factor": REPRESENTATION_EPS_FACTOR,
        "physical_moment_relative_error_max": REPRESENTATION_RELATIVE_TOL,
    }
    for key, value in expected.items():
        if gates.get(key) != value:
            raise StudyError(f"case gate {key} differs from source constant {value!r}")
    if (gates.get("spectrum_zero_eigenvalues"), gates.get("spectrum_positive_eigenvalues")) != (RIGID_COUNT, DOF - RIGID_COUNT):
        raise StudyError("case spectrum counts differ from source DOF contract")
    if (gates.get("significant_negative_eigenvalues"), gates.get("spectrum_zero_negative_tolerance")) != (0, "64*eps*93*max(abs(eigenvalues))"):
        raise StudyError("case spectrum tolerance contract differs from source")
    if gates.get("future_single_run_external_timeout_seconds") != 60:
        raise StudyError("case runtime budget differs from source contract")
def _slot(scalar_index: int, component_index: int) -> int:
    return 3 * scalar_index + component_index if scalar_index < 27 else 81 + 3 * (scalar_index - 27) + component_index
def _branch_values(np, xyz, row_ids, base_values, base_gradients, branch_fn, radius):
    count = len(xyz)
    values = np.empty((count, SCALAR), dtype=np.float64)
    gradients = np.empty((count, SCALAR, 3), dtype=np.float64)
    sizes = []
    for start in range(0, count, MAX_BATCH):
        stop = min(start + MAX_BATCH, count)
        sizes.append(stop - start); ids = row_ids[start:stop]
        values[start:stop, :27] = base_values[ids]
        gradients[start:stop, :27, :] = base_gradients[ids]
        for local_index, (x, y, _z) in enumerate(xyz[start:stop]):
            branch = branch_fn(float(x), float(y), reference_length_m=radius)
            row = start + local_index; values[row, 27:] = branch.values
            gradients[row, 27:, :] = branch.gradients_per_m
    _finite(np, "enriched basis values", values)
    _finite(np, "enriched basis gradients", gradients)
    return values, gradients, sizes
def _b_matrices(np, gradients):
    batch = gradients.shape[0]
    bmat = np.zeros((batch, 6, DOF), dtype=np.float64)
    for scalar_index in range(SCALAR):
        gx, gy, gz = (gradients[:, scalar_index, axis] for axis in range(3))
        for component in range(3):
            column = _slot(scalar_index, component)
            if component == 0:
                bmat[:, 0, column], bmat[:, 3, column], bmat[:, 5, column] = gx, gy, gz
            elif component == 1:
                bmat[:, 1, column], bmat[:, 3, column], bmat[:, 4, column] = gy, gx, gz
            else:
                bmat[:, 2, column], bmat[:, 4, column], bmat[:, 5, column] = gz, gy, gx
    return bmat
def _assemble_raw_k(np, gradients, weights, constitutive):
    matrix = np.zeros((DOF, DOF), dtype=np.float64)
    batch_sizes = []
    for start in range(0, len(weights), MAX_BATCH):
        stop = min(start + MAX_BATCH, len(weights)); batch_sizes.append(stop - start)
        bmat = _b_matrices(np, gradients[start:stop])
        matrix += np.einsum("qai,ab,qbj,q->ij", bmat, constitutive, bmat, weights[start:stop], optimize=True)
    _finite(np, "raw K", matrix)
    asym = float(np.linalg.norm(matrix - matrix.T) / max(np.linalg.norm(matrix), np.finfo(float).tiny))
    return matrix, asym, batch_sizes
def _make_exact_coefficients(np, ids, material, K):
    mu = material["E_Pa"] / (2.0 * (1.0 + material["nu"]))
    kappa = 3.0 - 4.0 * material["nu"]
    root = math.sqrt(material["R_m"] / (2.0 * math.pi)) / (2.0 * mu)
    ai, aii = K[0] * root, K[1] * root
    aiii = K[2] / mu * math.sqrt(2.0 * material["R_m"] / math.pi)
    branches = {
        1: (aii * (kappa + 1.0), ai * (kappa + 1.0), aiii),
        2: (ai * (kappa - 1.0), aii * (1.0 - kappa), 0.0),
        3: (ai, aii, 0.0),
        4: (aii, -ai, 0.0),
    }
    q = np.zeros(DOF, dtype=np.float64)
    for slot, identifier in enumerate(ids):
        scope, owner, component, basis_index = json.loads(identifier)
        if scope == "GLOBAL_FRONT":
            if owner != "global" or basis_index not in branches:
                raise StudyError("unexpected global-front semantic ID")
            q[slot] = branches[basis_index][{"ux": 0, "uy": 1, "uz": 2}[component]]
        elif scope != "SOURCE_BASE" or basis_index != 0: raise StudyError("unexpected base semantic ID")
    _finite(np, "manufactured q_exact", q)
    return q
def _reference(np, field_module, material, K, xyz):
    case = field_module.MixedLEFMCase(radius_m=material["R_m"], length_m=material["L_m"],
        K_I_Pa_sqrt_m=K[0], K_II_Pa_sqrt_m=K[1], K_III_Pa_sqrt_m=K[2],
        young_modulus_Pa=material["E_Pa"], poisson_ratio=material["nu"])
    displacement, strain, stress = (np.empty((len(xyz), n), dtype=np.float64) for n in (3, 6, 6))
    for row, point in enumerate(xyz):
        result = field_module.evaluate_local_point(tuple(float(v) for v in point), case=case)
        status = getattr(result.status, "value", result.status)
        if status != "OK" or result.displacement_world_m is None or result.strain_voigt_engineering is None or result.stress_voigt_Pa is None:
            raise StudyError(f"analytic reference status {status} at cached row {row}")
        displacement[row] = result.displacement_world_m
        strain[row] = result.strain_voigt_engineering
        stress[row] = result.stress_voigt_Pa
    _finite(np, "exact displacement", displacement)
    _finite(np, "exact strain", strain)
    _finite(np, "exact stress", stress)
    return displacement, strain, stress
def _force_from_traction(np, values, weights, normals, exact_stress):
    tensor = np.zeros((len(weights), 3, 3), dtype=np.float64)
    tensor[:, 0, 0] = exact_stress[:, 0]
    tensor[:, 1, 1] = exact_stress[:, 1]
    tensor[:, 2, 2] = exact_stress[:, 2]
    tensor[:, 0, 1] = tensor[:, 1, 0] = exact_stress[:, 3]; tensor[:, 1, 2] = tensor[:, 2, 1] = exact_stress[:, 4]
    tensor[:, 0, 2] = tensor[:, 2, 0] = exact_stress[:, 5]
    traction = np.einsum("qij,qj->qi", tensor, normals)
    force = np.zeros(DOF, dtype=np.float64)
    for scalar_index in range(SCALAR):
        for component in range(3):
            force[_slot(scalar_index, component)] = np.dot(weights, values[:, scalar_index] * traction[:, component])
    _finite(np, "boundary traction", traction)
    _finite(np, "raw boundary load", force)
    return force, traction
def _relative_l2(np, error, reference, weights, component_weights=None):
    sq, ref_sq = np.square(error), np.square(reference)
    if component_weights is not None:
        sq = sq * component_weights
        ref_sq = ref_sq * component_weights
    denominator = float(np.sum(weights[:, None] * ref_sq))
    numerator = float(np.sum(weights[:, None] * sq))
    if denominator <= np.finfo(float).tiny:
        return None
    return math.sqrt(max(0.0, numerator / denominator))
def _reconstruct_field(np, values, gradients, q, constitutive):
    count = len(values)
    displacement = np.zeros((count, 3), dtype=np.float64)
    strain = np.empty((count, 6), dtype=np.float64)
    stress = np.empty((count, 6), dtype=np.float64)
    batches = []
    for start in range(0, count, MAX_BATCH):
        stop = min(start + MAX_BATCH, count)
        batches.append(stop - start)
        for scalar_index in range(SCALAR):
            for component in range(3):
                displacement[start:stop, component] += values[start:stop, scalar_index] * q[_slot(scalar_index, component)]
        bmat = _b_matrices(np, gradients[start:stop])
        strain[start:stop] = np.einsum("qaj,j->qa", bmat, q)
        stress[start:stop] = strain[start:stop] @ constitutive.T
    for name, field in (("displacement", displacement), ("strain", strain), ("stress", stress)):
        _finite(np, f"reconstructed {name}", field)
    return displacement, strain, stress, batches

def _relative_vector_gate(np, name, actual, reference):
    _finite(np, name, actual); _finite(np, f"{name} reference", reference)
    error, reference_norm = float(np.linalg.norm(actual - reference)), float(np.linalg.norm(reference))
    absolute_scale = float(np.max(np.abs(reference))) if np.size(reference) else 0.0
    tolerance = REPRESENTATION_RELATIVE_TOL * reference_norm + (
        REPRESENTATION_EPS_FACTOR * np.finfo(float).eps * absolute_scale
        * math.sqrt(max(1, np.size(reference)))
    )
    relative = error / max(reference_norm, np.finfo(float).tiny)
    return {"error_norm": error, "reference_norm": reference_norm, "relative_error": relative,
            "absolute_tolerance": tolerance, "passed": error <= tolerance}

def _weighted_field_gate(np, name, actual, reference, weights):
    _finite(np, name, actual)
    _finite(np, f"{name} reference", reference)
    if actual.shape != reference.shape or len(weights) != len(actual): raise StudyError(f"{name} shape mismatch in representation audit")
    delta_norm = math.sqrt(float(np.sum(weights[:, None] * np.square(actual - reference))))
    reference_norm = math.sqrt(float(np.sum(weights[:, None] * np.square(reference))))
    max_abs = float(np.max(np.abs(reference))) if np.size(reference) else 0.0
    absolute_tolerance = (REPRESENTATION_EPS_FACTOR * np.finfo(float).eps * max_abs
                          * math.sqrt(max(np.finfo(float).tiny, float(np.sum(weights)) * actual.shape[1])))
    tolerance = REPRESENTATION_RELATIVE_TOL * reference_norm + absolute_tolerance
    return {"error_norm": delta_norm, "reference_norm": reference_norm,
            "relative_error": delta_norm / max(reference_norm, np.finfo(float).tiny),
            "absolute_tolerance": absolute_tolerance, "tolerance": tolerance, "passed": delta_norm <= tolerance}

def _physical_moment_matrix(np, values, xyz, weights, radius):
    volume = float(np.sum(weights))
    if not math.isfinite(volume) or volume <= 0.0:
        raise StudyError("physical moment gauge has invalid volume measure")
    center, relative = np.zeros(3, dtype=np.float64), xyz
    C = np.zeros((RIGID_COUNT, DOF), dtype=np.float64)
    for scalar_index in range(SCALAR):
        phi_w = weights * values[:, scalar_index] / volume
        for component in range(3):
            column = _slot(scalar_index, component)
            C[component, column] = float(np.sum(phi_w))
            if component == 0:
                C[4, column], C[5, column] = float(np.dot(phi_w, relative[:, 2] / radius)), float(np.dot(phi_w, -relative[:, 1] / radius))
            elif component == 1:
                C[3, column], C[5, column] = float(np.dot(phi_w, -relative[:, 2] / radius)), float(np.dot(phi_w, relative[:, 0] / radius))
            else:
                C[3, column], C[4, column] = float(np.dot(phi_w, relative[:, 1] / radius)), float(np.dot(phi_w, -relative[:, 0] / radius))
    _finite(np, "physical volume-moment matrix", C)
    return C, volume, center

def _direct_physical_moments(np, displacement, xyz, weights, radius, volume):
    moments = np.concatenate((np.einsum("q,qi->i", weights, displacement) / volume,
        np.einsum("q,qi->i", weights, np.cross(xyz, displacement)) / (volume * radius)))
    _finite(np, "direct analytic physical moments", moments)
    return moments
def _scaled_spectrum_gate(np, matrix):
    diagonal = np.diag(matrix)
    _finite(np, "raw stiffness diagonal", diagonal)
    if bool(np.any(diagonal <= 0.0)): raise StudyError("raw stiffness has a nonpositive diagonal; diagonal scaling undefined")
    scale_vector = 1.0 / np.sqrt(diagonal); scaled = scale_vector[:, None] * matrix * scale_vector[None, :]
    _finite(np, "diagonally scaled raw stiffness", scaled)
    asymmetry = float(np.linalg.norm(scaled - scaled.T) / max(np.linalg.norm(scaled), np.finfo(float).tiny))
    if asymmetry > ASYMMETRY_TOL: raise StudyError(f"scaled raw K asymmetry {asymmetry:.6g} exceeds {ASYMMETRY_TOL:g}")
    eigenvalues = np.linalg.eigvalsh(scaled)
    _finite(np, "eigenvalues", eigenvalues)
    scale = float(np.max(np.abs(eigenvalues)))
    tol = SPECTRUM_FACTOR * np.finfo(float).eps * DOF * scale
    zero, negative, positive = np.abs(eigenvalues) <= tol, eigenvalues < -tol, eigenvalues > tol
    passed = int(zero.sum()) == RIGID_COUNT and int(positive.sum()) == DOF - RIGID_COUNT and not bool(negative.any())
    return scaled, scale_vector, eigenvalues, {"scaled_relative_asymmetry": asymmetry,
                         "scale": scale, "tolerance": tol,
                         "zero_count": int(zero.sum()), "positive_count": int(positive.sum()),
                         "significant_negative_count": int(negative.sum()), "passed": passed}
def _solve_case(np, Kraw, Kscaled, C, force, scale_vector, g):
    Cbar = C * scale_vector[None, :]; row_norms = np.linalg.norm(Cbar, axis=1)
    if bool(np.any(~np.isfinite(row_norms))) or bool(np.any(row_norms <= np.finfo(float).tiny)):
        raise StudyError("physical gauge has an empty or invalid scaled row")
    row_scale = 1.0 / row_norms; Chat = row_scale[:, None] * Cbar; ghat = row_scale * g
    saddle = np.zeros((DOF + RIGID_COUNT, DOF + RIGID_COUNT), dtype=np.float64)
    saddle[:DOF, :DOF] = Kscaled
    saddle[:DOF, DOF:] = Chat.T
    saddle[DOF:, :DOF] = Chat
    rhs = np.concatenate((scale_vector * force, ghat)); solution = np.linalg.solve(saddle, rhs)
    _finite(np, "KKT solution", solution)
    y, scaled_multipliers = solution[:DOF], solution[DOF:]
    q = scale_vector * y
    multipliers = row_scale * scaled_multipliers
    _finite(np, "restored physical solution", q)
    _finite(np, "restored physical multipliers", multipliers)
    residual = saddle @ solution - rhs
    kkt_rel = float(np.linalg.norm(residual, ord=np.inf) / max(
        np.linalg.norm(saddle, ord=np.inf) * np.linalg.norm(solution, ord=np.inf)
        + np.linalg.norm(rhs, ord=np.inf), np.finfo(float).tiny))
    weak_residual = Kraw @ q + C.T @ multipliers - force
    weak_rel = float(np.linalg.norm(weak_residual, ord=np.inf) / max(np.linalg.norm(Kraw, ord=np.inf) * np.linalg.norm(q, ord=np.inf) + np.linalg.norm(C.T @ multipliers, ord=np.inf) + np.linalg.norm(force, ord=np.inf), np.finfo(float).tiny))
    gauge_vector = C @ q - g
    gauge_rel = float(np.linalg.norm(gauge_vector, ord=np.inf) / max(np.linalg.norm(C, ord=np.inf) * np.linalg.norm(q, ord=np.inf) + np.linalg.norm(g, ord=np.inf), np.finfo(float).tiny))
    strain_energy_twice = float(q @ (Kraw @ q))
    load_work = float(q @ force)
    gauge_work = float(g @ multipliers)
    energy_rel = abs(strain_energy_twice - load_work + gauge_work) / max(
        abs(strain_energy_twice), abs(load_work), abs(gauge_work), np.finfo(float).tiny)
    if max(kkt_rel, weak_rel, gauge_rel, energy_rel) > KKT_TOL:
        raise StudyError(f"KKT/energy gate failed: scaled_saddle={kkt_rel:.6g}, physical_weak={weak_rel:.6g}, gauge={gauge_rel:.6g}, energy={energy_rel:.6g}")
    metrics = {"kkt_backward_error": kkt_rel, "physical_weak_relative_error": weak_rel,
                            "gauge_relative_error": gauge_rel,
                            "twice_strain_energy_J": strain_energy_twice,
                            "load_work_J": load_work, "gauge_work_J": gauge_work,
                            "energy_identity_relative_error": energy_rel}
    arrays = {"Cbar": Cbar, "row_scale": row_scale, "Chat": Chat,
              "ghat": ghat, "saddle_scaled": saddle, "rhs_scaled": rhs,
              "solution_scaled": solution, "scaled_multipliers": scaled_multipliers}
    return q, multipliers, metrics, arrays
def _load_cloud(np, archive, name, with_boundary):
    xyz = np.asarray(archive[f"cloud_{name}_xyz_m"], dtype=np.float64)
    row_ids = np.asarray(archive[f"cloud_{name}_unique_row_ids"], dtype=np.int64)
    if with_boundary:
        weights = np.asarray(archive[f"cloud_{name}_weights"], dtype=np.float64)
        normals = np.asarray(archive[f"cloud_{name}_normals"], dtype=np.float64)
        if normals.shape != xyz.shape:
            raise StudyError(f"bad cached normals in cloud {name}")
        return xyz, row_ids, weights, normals
    weights = np.asarray(archive[f"cloud_{name}_weights"], dtype=np.float64)
    return xyz, row_ids, weights, None

def _run(case: dict, cache_result_path: Path, cache_arrays_path: Path):
    import numpy as np
    material = dict(case["material"])
    material["R_m"], material["L_m"] = float(case["geometry"]["radius_m"]), float(case["geometry"]["length_m"])
    field_module = _module("mixed93_pinned_field", ROOT / SOURCE_PATHS["mixed_lefm_field"])
    operator = _module("mixed93_pinned_semantic_operator", ROOT / SOURCE_PATHS["global_semantic_operator"])
    cached_result = json.loads(cache_result_path.read_text(encoding="utf-8"))
    if (cached_result.get("schema"), cached_result.get("status"), cached_result.get("assembly_status"), cached_result.get("mechanical_eligibility"), cached_result.get("solve_allowed")) != ("tip_aligned_full3d_quadrature_study_result.v1", "OK_DIAGNOSTIC_ONLY_UNQUALIFIED", "OK", "SPECTRUM_GATE_PASSED", True): raise StudyError("frozen cache result is not an eligible diagnostic source")
    source_ids = tuple(cached_result["registry"]["source_ids"])
    if len(source_ids) != 27 or source_ids != tuple(sorted(source_ids)): raise StudyError("cached source ID order is not canonical 27-source order")
    with np.load(cache_arrays_path, allow_pickle=False) as archive:
        source_positions = np.asarray(archive["source_positions_xyz_m"], dtype=np.float64); unique_xyz = np.asarray(archive["unique_query_positions_xyz_m"], dtype=np.float64)
        base_values = np.asarray(archive["base_mls_values"], dtype=np.float64); base_gradients = np.asarray(archive["base_mls_gradients_per_m"], dtype=np.float64)
        if (source_positions.shape, unique_xyz.shape, base_values.shape, base_gradients.shape) != ((27,3),(4224,3),(4224,27),(4224,27,3)): raise StudyError("frozen cached array shape contract changed")
        for name, array in (("source positions",source_positions),("unique positions",unique_xyz),("base MLS values",base_values),("base MLS gradients",base_gradients)): _finite(np, name, array)
        vol, side, caps = (_load_cloud(np,archive,n,b) for n,b in (("a0_volume32",False),("side",True),("a0_cap16",True)))
        off_xyz = np.asarray(archive["cloud_offgrid_xyz_m"],dtype=np.float64); off_ids = np.asarray(archive["cloud_offgrid_unique_row_ids"],dtype=np.int64)
    groups = {"volume":(vol[0],vol[1],vol[2],None), "boundary":tuple(np.concatenate((side[i],caps[i])) for i in range(4)),
        "offgrid":(off_xyz,off_ids,np.ones(len(off_xyz),dtype=np.float64),None)}
    expected_counts = {"volume": 512, "boundary": 768, "offgrid": 128}
    basis, query_batches = {}, {}
    for name, (xyz, row_ids, weights, normals) in groups.items():
        if len(xyz) != expected_counts[name] or len(row_ids) != len(xyz) or len(weights) != len(xyz): raise StudyError(f"cached {name} cloud count mismatch")
        if row_ids.min() < 0 or row_ids.max() >= len(unique_xyz) or not np.array_equal(xyz, unique_xyz[row_ids]): raise StudyError(f"cached {name} points disagree with cached MLS row references")
        _finite(np, f"{name} quadrature weights", weights)
        if name != "offgrid" and bool(np.any(weights <= 0.0)): raise StudyError(f"nonpositive {name} quadrature weight")
        if normals is not None:
            _finite(np, f"{name} normals", normals)
            if float(np.max(np.abs(np.linalg.norm(normals, axis=1) - 1.0))) > 1.0e-10: raise StudyError(f"{name} normals are not unit vectors")
        basis[name] = _branch_values(np,xyz,row_ids,base_values,base_gradients,field_module.branch_functions_local,material["R_m"])
        query_batches[name] = basis[name][2]
    ids = operator.semantic_coefficient_ids(source_ids)
    if len(ids) != DOF: raise StudyError(f"semantic registry has {len(ids)} coefficients, expected {DOF}")
    D = np.asarray(operator.isotropic_constitutive_matrix_pa(material["E_Pa"], material["nu"]), dtype=np.float64)
    vol_values, vol_gradients = basis["volume"][:2]
    Kraw, raw_asymmetry, k_batches = _assemble_raw_k(np, vol_gradients, vol[2], D)
    Kscaled, scale_vector, eigenvalues, spectrum = _scaled_spectrum_gate(np, Kraw)
    rigid = operator.rigid_body_mode_coefficients(source_ids, tuple(map(tuple, source_positions)))
    if rigid.matrix_shape != (RIGID_COUNT, DOF) or rigid.coefficient_ids != ids: raise StudyError("frozen six-rigid coefficient layout mismatch")
    Ccoef = np.frombuffer(rigid.values_le_f64, dtype="<f8").copy().reshape(RIGID_COUNT, DOF)
    _finite(np, "rigid coefficient modes", Ccoef)
    rigid_singular = np.linalg.svd(Ccoef, compute_uv=False); rigid_rank_tol = RIGID_RELATIVE_TOL * float(rigid_singular[0])
    rigid_rank = int(np.count_nonzero(rigid_singular > rigid_rank_tol))
    rigid_y = Ccoef.T / scale_vector[:, None]
    rigid_null_error = float(np.linalg.norm(Kscaled @ rigid_y, ord=np.inf) / max(
        np.linalg.norm(Kscaled, ord=np.inf) * np.linalg.norm(rigid_y, ord=np.inf), np.finfo(float).tiny))
    Cphys, volume_measure, gauge_center = _physical_moment_matrix(np,vol_values,vol[0],vol[2],material["R_m"])
    rigid_coupling = Cphys @ Ccoef.T; coupling_singular = np.linalg.svd(rigid_coupling, compute_uv=False)
    coupling_tol = RIGID_RELATIVE_TOL * float(coupling_singular[0])
    physical_coupling_rank = int(np.count_nonzero(coupling_singular > coupling_tol))
    arrays = {
        "K_N_per_m": Kraw, "K_diagonally_scaled": Kscaled, "K_diagonal_scale_sqrt_m_per_sqrt_N": scale_vector,
        "constitutive_matrix_Pa": D, "eigenvalues_scaled": eigenvalues,
        "source_positions_xyz_m": source_positions, "source_ids": np.asarray(source_ids), "semantic_ids": np.asarray(ids),
        "rigid_coefficient_modes": Ccoef, "physical_volume_moment_matrix": Cphys,
        "physical_gauge_rigid_coupling": rigid_coupling,
        "volume_xyz_m": vol[0], "volume_unique_row_ids": vol[1], "volume_weights_m3": vol[2],
        "boundary_xyz_m": groups["boundary"][0], "boundary_unique_row_ids": groups["boundary"][1],
        "boundary_normals": groups["boundary"][3], "boundary_weights_m2": groups["boundary"][2],
        "offgrid_xyz_m": off_xyz, "offgrid_unique_row_ids": off_ids,
    }
    for name, (values, gradients, _sizes) in basis.items():
        arrays[f"{name}_basis_values"] = values; arrays[f"{name}_basis_gradients_per_m"] = gradients
    rigid_passed = rigid_rank == RIGID_COUNT and rigid_null_error <= RIGID_RELATIVE_TOL and physical_coupling_rank == RIGID_COUNT
    mechanics_passed = bool(spectrum["passed"] and rigid_passed)
    metrics = {
        "schema": SCHEMA, "case_id": case["case_id"], "classification": case["classification"],
        "status": "ASSEMBLY_READY_FOR_DIAGNOSTIC" if mechanics_passed else "BLOCKED_MECHANICAL_GATE",
        "assembly_status": "OK", "mechanical_eligibility": "SCALED_SPECTRUM_AND_SIX_RIGID_GATES_PASSED" if mechanics_passed else "BLOCKED",
        "solve_allowed": bool(mechanics_passed), "spectrum": spectrum,
        "rigid_coefficient_mode_rank": rigid_rank, "rigid_coefficient_mode_rank_tolerance": rigid_rank_tol,
        "rigid_scaled_null_relative_error": rigid_null_error, "rigid_scaled_null_relative_tolerance": RIGID_RELATIVE_TOL,
        "physical_gauge_rigid_coupling_rank": physical_coupling_rank, "physical_gauge_rigid_coupling_rank_tolerance": coupling_tol,
        "raw_K_relative_asymmetry_diagnostic_only": raw_asymmetry,
        "assembly_volume_batches": k_batches, "cached_cloud_batch_sizes": query_batches,
        "reconstruction_batch_limit_rows": MAX_BATCH,
        "raw_cached_rows": 1408, "unique_cached_mls_rows_touched": int(np.unique(np.concatenate((vol[1], groups["boundary"][1], off_ids))).size),
        "new_mls_queries": 0, "source_count": len(source_ids), "source_ids": list(source_ids), "dof_count": DOF,
        "rigid_mode_ids": list(rigid.mode_ids), "physical_gauge": {"definition": "volume means of [u_x,u_y,u_z, ((x-center) cross u)/R]",
            "center_xyz_m": gauge_center.tolist(), "reference_radius_m": material["R_m"], "volume_measure_m3": volume_measure,
            "constraint_targets_nonzero_manufactured": True},
        "region_row_counts": expected_counts,
        "boundary_region_slices_half_open": {"side": [0, 256], "cap_bottom": [256, 512], "cap_top": [512, 768]},
        "case_metrics": [], "limitations": ["toy isotropic LEFM benchmark only", "not crack-growth or material qualification"],
    }
    if not mechanics_passed:
        return metrics, arrays, 3

    # Complete every q_exact representation audit before the first linear solve.
    prepared = []
    exact_audits_passed = True
    strain_weights = np.asarray((1.0, 1.0, 1.0, 0.5, 0.5, 0.5))
    stress_weights = np.asarray((1.0, 1.0, 1.0, 2.0, 2.0, 2.0))
    for case_def in case["load_cases"]:
        Ktuple = tuple(float(v) for v in case_def["K_Pa_sqrt_m"]); qexact = _make_exact_coefficients(np, ids, material, Ktuple)
        exact_fields, audit_rows = {}, {}
        for name, (xyz, _row_ids, weights, _normals) in groups.items():
            values, gradients = basis[name][:2]
            ref_u, ref_strain, ref_stress = _reference(np, field_module, material, Ktuple, xyz)
            rec_u, rec_strain, rec_stress, rec_batches = _reconstruct_field(np, values, gradients, qexact, D)
            field_gates = {field:_weighted_field_gate(np,f"{name} q_exact {field}",actual,reference,weights)
                for field,actual,reference in (("displacement",rec_u,ref_u),("strain",rec_strain,ref_strain),("stress",rec_stress,ref_stress))}
            exact_fields[name] = {"u": ref_u, "strain": ref_strain, "stress": ref_stress}
            audit_rows[name] = {**field_gates, "reconstruction_batches": rec_batches}
            for key,value in exact_fields[name].items(): arrays[f"{case_def['case_id']}_{name}_exact_{key}"] = value
        g = Cphys @ qexact; direct_g = _direct_physical_moments(np,exact_fields["volume"]["u"],vol[0],vol[2],material["R_m"],volume_measure)
        moment_check = _relative_vector_gate(np,"physical moment reproduction",g,direct_g)
        audit_ok = moment_check["passed"] and all(item[f]["passed"] for item in audit_rows.values() for f in ("displacement","strain","stress"))
        exact_audits_passed = exact_audits_passed and audit_ok
        force, traction = _force_from_traction(np,basis["boundary"][0],groups["boundary"][2],groups["boundary"][3],exact_fields["boundary"]["stress"])
        prefix = case_def["case_id"]
        arrays.update({f"{prefix}_q_exact":qexact,f"{prefix}_physical_gauge_target_m":g,f"{prefix}_analytic_volume_moments_m":direct_g,
            f"{prefix}_boundary_force_N":force,f"{prefix}_boundary_traction_Pa":traction})
        prepared.append({"definition": case_def, "qexact": qexact, "g": g, "force": force, "fields": exact_fields})
        metrics["case_metrics"].append({"case_id":prefix,"K_Pa_sqrt_m":list(Ktuple),"exact_representation":audit_rows,
            "physical_gauge_target_m":g.tolist(),"analytic_volume_moments_m":direct_g.tolist(),"physical_moment_reproduction":moment_check,"exact_representation_gate_passed":audit_ok})
    if not exact_audits_passed:
        metrics.update(status="BLOCKED_EXACT_REPRESENTATION_GATE", mechanical_eligibility="SPECTRUM_AND_RIGID_GATES_PASSED_BUT_EXACT_REPRESENTATION_BLOCKED", solve_allowed=False)
        return metrics, arrays, 3

    solve_arrays = {}
    for item, metrics_case in zip(prepared, metrics["case_metrics"]):
        q, multipliers, solve_metrics, kkt_arrays = _solve_case(np,Kraw,Kscaled,Cphys,item["force"],scale_vector,item["g"])
        prefix = item["definition"]["case_id"]
        metrics_case.update({**solve_metrics,"q_solution":q.tolist(),"lagrange_multipliers_physical":multipliers.tolist(),
            "q_exact":item["qexact"].tolist(),"physical_gauge_target_m":item["g"].tolist(),"rigid_coefficient_mode_force_work":(Ccoef@item["force"]).tolist()})
        arrays.update({f"{prefix}_q_solution":q,f"{prefix}_lambda_physical":multipliers,f"{prefix}_lambda_scaled":kkt_arrays["scaled_multipliers"],
            f"{prefix}_kkt_rhs_scaled":kkt_arrays["rhs_scaled"],f"{prefix}_kkt_solution_scaled":kkt_arrays["solution_scaled"],f"{prefix}_gauge_rhs_scaled":kkt_arrays["ghat"]})
        solve_arrays = kkt_arrays
        computed = {}
        for name, (xyz, _row_ids, weights, _normals) in groups.items():
            values, gradients = basis[name][:2]; u, strain, stress, batches = _reconstruct_field(np,values,gradients,q,D)
            exact = item["fields"][name]
            computed[name] = (u, strain, stress)
            metrics_case[f"{name}_solution_reconstruction_batches"] = batches
            arrays.update({f"{prefix}_{name}_solution_u_m":u,f"{prefix}_{name}_solution_strain_m_per_m":strain,f"{prefix}_{name}_solution_stress_Pa":stress})
            for field, actual, reference, comp_weights in (("displacement",u,exact["u"],None),("strain_tensor",strain,exact["strain"],strain_weights),("stress_tensor",stress,exact["stress"],stress_weights)):
                metrics_case[f"{name}_solution_{field}_relative_L2"] = _relative_l2(np,actual-reference,reference,weights,comp_weights)
        vol_u, vol_strain, vol_stress = computed["volume"]
        _off_u, _off_strain, off_stress = computed["offgrid"]; vol_exact = item["fields"]["volume"]
        metrics_case["volume_sigma_component_max_abs_error_Pa"] = np.max(np.abs(vol_stress - vol_exact["stress"]), axis=0).tolist()
        metrics_case["volume_sigma_component_rms_Pa"] = np.sqrt(np.average(np.square(vol_stress - vol_exact["stress"]), axis=0, weights=vol[2])).tolist()
        metrics_case["volume_strain_component_max_abs_error_m_per_m"] = np.max(np.abs(vol_strain - vol_exact["strain"]), axis=0).tolist()
        metrics_case["volume_strain_component_rms_m_per_m"] = np.sqrt(np.average(np.square(vol_strain - vol_exact["strain"]), axis=0, weights=vol[2])).tolist()
        errors = [metrics_case[f"{group}_solution_{field}_relative_L2"] for group,field in (("volume","displacement"),("offgrid","displacement"),("volume","strain_tensor"),("volume","stress_tensor"),("offgrid","stress_tensor"))]
        metrics_case["engineering_error_gate_passed"] = all(value is not None and value <= FIELD_RELATIVE_LIMIT for value in errors)
    arrays.update({"Cbar_scaled":solve_arrays["Cbar"],"physical_gauge_row_scale":solve_arrays["row_scale"],"Chat_scaled":solve_arrays["Chat"],"scaled_KKT_matrix":solve_arrays["saddle_scaled"]})
    memory_estimate = 2*sum(np.asarray(v).nbytes for v in arrays.values()) + base_values.nbytes + base_gradients.nbytes + unique_xyz.nbytes + 2*MAX_BATCH*DOF*6*8 + MAX_BATCH*6*DOF*8
    if memory_estimate > MEMORY_ESTIMATE_LIMIT_BYTES: raise StudyError(f"explicit working-memory estimate {memory_estimate} exceeds cap")
    metrics["explicit_working_memory_estimate_bytes"] = int(memory_estimate)
    gates = [item["engineering_error_gate_passed"] for item in metrics["case_metrics"]]
    metrics.update(engineering_error_gate_passed=all(gates),status="OK_DIAGNOSTIC_ONLY" if all(gates) else "FAILED_ENGINEERING_ERROR_GATE",
        mechanical_eligibility="SCALED_SPECTRUM_RIGID_PHYSICAL_GAUGE_AND_SOLVE_GATES_PASSED",solve_allowed=True)
    return metrics, arrays, 0 if all(gates) else 3
def _write_result(output_json: Path, metrics: dict, arrays: dict, return_code: int) -> None:
    import numpy as np
    def json_default(value):
        if isinstance(value, np.generic):
            return value.item()
        raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")
    output_json.parent.mkdir(parents=True, exist_ok=False)
    array_path = output_json.with_name("arrays.npz")
    np.savez_compressed(array_path, **arrays)
    if array_path.stat().st_size > OUTPUT_LIMIT_BYTES:
        raise StudyError(f"compressed array output exceeds {OUTPUT_LIMIT_BYTES} bytes")
    metrics["arrays_npz"] = {"path": array_path.name, "sha256": _sha(array_path), "bytes": array_path.stat().st_size}
    metrics["output_budget_bytes"] = OUTPUT_LIMIT_BYTES
    metrics["explicit_memory_estimate_limit_bytes"] = MEMORY_ESTIMATE_LIMIT_BYTES
    metrics["return_code"] = return_code
    metrics["input_sha256"] = {
        "cached_result": _sha(ROOT / SOURCE_PATHS["cached_result"]),
        "cached_arrays": _sha(ROOT / SOURCE_PATHS["cached_arrays"]),
        "driver": _sha(Path(__file__).resolve()), "case": _sha(CASE_PATH), "source_pins": _sha(PINS_PATH),
    }
    output_json.write_text(json.dumps(metrics, indent=2, sort_keys=True, allow_nan=False, default=json_default) + "\n", encoding="utf-8")
    total = array_path.stat().st_size + output_json.stat().st_size
    if total > OUTPUT_LIMIT_BYTES:
        raise StudyError(f"combined output exceeds {OUTPUT_LIMIT_BYTES} bytes")

def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--execute", action="store_true", help="run only after case authorization is true")
    parser.add_argument("--output", type=Path, default=REVIEW / "run01" / "result.json")
    args = parser.parse_args(argv)
    if not args.execute:
        print("SOURCE_READY_EXECUTION_DISABLED")
        return 0
    case = json.loads(CASE_PATH.read_text(encoding="utf-8"))
    if case.get("schema") != SCHEMA or case.get("execution_authorization", {}).get("authorized") is not True:
        raise StudyError("case is not authorized for numerical execution")
    output = args.output.resolve()
    try:
        output.relative_to(REVIEW.resolve())
    except ValueError as exc:
        raise StudyError("output must remain inside this study's review directory") from exc
    if output.exists() or output.parent.exists():
        raise StudyError("refusing an existing output or run directory")
    cache_result, cache_arrays = _pin_inputs(case)
    metrics, arrays, result_rc = _run(case, cache_result, cache_arrays)
    _write_result(output, metrics, arrays, result_rc)
    print(json.dumps({"status": metrics["status"], "result": str(output), "return_code": result_rc}, sort_keys=True))
    return result_rc
if __name__ == "__main__":
    raise SystemExit(main())
