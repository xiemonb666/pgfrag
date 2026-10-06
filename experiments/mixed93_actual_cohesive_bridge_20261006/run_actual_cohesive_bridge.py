"""Numerical worker for the frozen 93-DOF cohesive bridge; root-go gated."""
from __future__ import annotations
import argparse, hashlib, importlib, importlib.util, json, math, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
REVIEW = ROOT / "pgfrag" / "tmp" / "review_20261006" / "mixed93_actual_cohesive_bridge_20261006"
CASE_PATH, PINS_PATH = HERE / "case.json", REVIEW / "NUMERIC_SOURCE_PINS.json"
SPECTRUM_FACTOR, KKT_TOL = 64.0, 1.0e-9
ASYMMETRY_TOL, OUTPUT_LIMIT, MEMORY_LIMIT = 1.0e-12, 12 * 1024 * 1024, 64 * 1024 * 1024
class GateFailure(RuntimeError):
    def __init__(self, gate: str, detail: object): super().__init__(f"{gate}: {detail}"); self.gate, self.detail = gate, detail
def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()
def record_gate(gates: dict, name: str, passed: bool, evidence: object) -> None:
    gates[name] = {"status": "PASS" if passed else "FAIL", "passed": bool(passed), "evidence": evidence}
    if not passed:
        raise GateFailure(name, evidence)
def require_execution(pins: dict, token: str | None) -> None:
    go = pins.get("root_go")
    if (pins.get("execution_authorized") is not True or not isinstance(go, dict)
            or go.get("authorized") is not True or not isinstance(go.get("token"), str)
            or not go["token"] or token != go["token"]):
        raise PermissionError("numeric execution is disabled until NUMERIC_SOURCE_PINS carries an exact root-go token")
def pinned_inputs(pins: dict) -> dict[str, str]:
    found = {}
    for name, ref in pins.get("references", {}).items():
        if not isinstance(ref, dict) or not isinstance(ref.get("path"), str):
            raise ValueError(f"malformed numeric source reference: {name}")
        path = ROOT / Path(ref["path"])
        if not path.is_file():
            raise FileNotFoundError(f"pinned source missing: {ref['path']}")
        actual = sha256(path)
        if actual != ref.get("sha256"):
            raise ValueError(f"pinned source hash mismatch: {name} ({actual})")
        found[name] = actual
    cache = pins.get("accepted_cache_binding", {})
    cache_path = ROOT / Path(cache.get("arrays_path", ""))
    if not cache_path.is_file():
        raise FileNotFoundError("accepted arrays.npz is missing")
    actual = sha256(cache_path)
    if actual != cache.get("arrays_sha256"):
        raise ValueError(f"accepted arrays.npz SHA-256 mismatch: {actual}")
    found["accepted_arrays_npz"] = actual
    return found
def load_registered(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load pinned module {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module
def finite(np, name: str, value) -> None:
    if not bool(np.isfinite(value).all()):
        raise GateFailure("finite_inputs", f"non-finite {name}")
def relative_error(np, actual, expected) -> float:
    return float(np.linalg.norm(actual - expected) / max(float(np.linalg.norm(expected)), np.finfo(float).tiny))
def point_ledger(assembly, points: list, coefficients: tuple) -> list[dict]:
    rows = []
    for point, item in zip(points, assembly.point_results):
        state = item.trial_state
        jump = tuple(sum(a * b for a, b in zip(row, coefficients)) for row in point.jump_operator)
        opening = sum(a * b for a, b in zip(point.normal, jump))
        shear = tuple(jump[i] - opening * point.normal[i] for i in range(3))
        rows.append({"point_id": list(item.point_id), "opening_m": item.opening_m, "traction_pa": item.traction_pa,
            "jump_vector_m": list(jump), "normal_opening_recomputed_m": opening, "tangential_jump_vector_m": list(shear),
            "k_alg_pa_per_m": item.k_alg_pa_per_m, "tangent_side": item.tangent_side,
            "tangent_status": getattr(item.tangent_status, "value", str(item.tangent_status)),
            "nondifferentiable_reasons": list(item.nondifferentiable_reasons),
            "state": {"opening_m": state.opening_m, "kappa_m": state.kappa_m, "damage": state.damage,
                "stored_energy_j_per_m2": state.stored_energy_j_per_m2,
                "dissipated_energy_j_per_m2": state.dissipated_energy_j_per_m2,
                "external_work_j_per_m2": state.external_work_j_per_m2, "update_count": state.update_count},
            "generalized_internal_force_n": list(item.generalized_force_n),
            "primitive_positive_face_internal_force_n": list(item.positive_face_force_at_reference_n),
            "physical_upper_face_traction_n_per_m2": [-item.traction_pa * n for n in point.normal],
            "physical_lower_face_traction_n_per_m2": [item.traction_pa * n for n in point.normal]})
    return rows

def spectrum_gate(np, matrix, name: str) -> dict:
    diagonal = np.diag(matrix)
    if bool(np.any(diagonal <= 0.0)):
        raise GateFailure("spectrum", f"{name} has nonpositive diagonal; no repair is allowed")
    scale_vector = 1.0 / np.sqrt(diagonal)
    scaled = scale_vector[:, None] * matrix * scale_vector[None, :]
    asymmetry = float(np.linalg.norm(scaled - scaled.T) / max(float(np.linalg.norm(scaled)), np.finfo(float).tiny))
    if asymmetry > ASYMMETRY_TOL:
        raise GateFailure("spectrum", f"{name} scaled asymmetry {asymmetry} exceeds {ASYMMETRY_TOL}")
    eigenvalues = np.linalg.eigvalsh(scaled)
    finite(np, f"{name} eigenvalues", eigenvalues)
    tolerance = SPECTRUM_FACTOR * np.finfo(float).eps * matrix.shape[0] * float(np.max(np.abs(eigenvalues)))
    zero = int(np.count_nonzero(np.abs(eigenvalues) <= tolerance))
    positive = int(np.count_nonzero(eigenvalues > tolerance))
    negative = int(np.count_nonzero(eigenvalues < -tolerance))
    return {"scaled_relative_asymmetry": asymmetry, "zero_tolerance": tolerance,
            "zero_count": zero, "positive_count": positive, "negative_count": negative,
            "passed": zero == 6 and positive == 87 and negative == 0}


def execute(case: dict, pins: dict, output_path: Path) -> tuple[dict, int]:
    gates = {name: {"status": "UNKNOWN", "passed": None, "evidence": None} for name in (
        "cache_identity", "scaled_cached_kkt", "mls", "exact_traces", "rigid_traces",
        "interface_quadrature", "primitive_baseline", "zero_tangent_control", "raw_spectrum",
        "coupled_spectrum", "coupled_kkt", "primitive_final", "mechanical_change", "energy", "resource_budget")}
    result = {"schema": "mixed93_actual_cohesive_bridge_numeric_result.v1",
              "classification": case["classification"], "status": "NOT_ACCEPTED",
              "operator_layout": "MIXED93_BASE81_GLOBAL12", "trace_family": "GLOBALFRONT",
              "numeric_modules_imported": False, "npz_opened": False,
              "counts": {"new_mls_queries": 0, "mls_queries_requested": 0, "mls_query_attempts": 0,
                         "interface_qps": 0, "primitive_assemblies": 0, "primitive_successful_returns": 0,
                         "coupled_kkt_solves": 0, "coupled_kkt_successful_returns": 0},
              "evidence": {}, "gates": gates}
    try:
        # Recheck every source pin and hash the cache before importing NumPy or opening NPZ.
        hashes = pinned_inputs(pins)
        record_gate(gates, "cache_identity", True, {"source_pin_count": len(hashes),
                    "npz_sha256": hashes["accepted_arrays_npz"], "npz_content_opened_before_hash": False})
        import numpy as np
        result["numeric_modules_imported"] = True
        experiments = ROOT / "pgfrag" / "experiments"
        mls_dir = experiments / "mls3d_kernel_20261005"
        sys.path.insert(0, str(mls_dir))
        mls_module = importlib.import_module("mls3d_mls")
        shape_module = importlib.import_module("mls3d_shape")
        field_module = load_registered("mixed93_bridge_field", experiments / "mixed_lefm_benchmark_20261006" / "mixed_lefm_field.py")
        operator = load_registered("mixed93_bridge_operator", experiments / "mixed_enriched_operator_20261006" / "mixed_enriched_operator.py")
        trace = load_registered("mixed93_bridge_trace", experiments / "enriched_cohesive_trace_20261006" / "enriched_cohesive_trace.py")
        # The trace primitive preloads its pinned derivative/scalar classes; reuse those exact class objects.
        scalar = trace._DERIVATIVE._SCALAR
        result["npz_opened"] = True
        cache = case["load_and_cache_binding"]
        arrays_path = ROOT / Path(cache["failed_serialization_arrays_path"])
        with np.load(arrays_path, allow_pickle=False) as archive:
            keys = cache["actual_array_keys"]
            source_positions = np.asarray(archive["source_positions_xyz_m"], dtype=np.float64)
            source_ids = tuple(str(v.decode("utf-8") if isinstance(v, bytes) else v) for v in archive["source_ids"].tolist())
            semantic_ids = tuple(str(v.decode("utf-8") if isinstance(v, bytes) else v) for v in archive["semantic_ids"].tolist())
            K = np.asarray(archive[keys["K"]], dtype=np.float64)
            C = np.asarray(archive[keys["physical_constraint_matrix"]], dtype=np.float64)
            rigid_saved = np.asarray(archive[keys["rigid_vectors"]], dtype=np.float64)
            q_solution = np.asarray(archive[keys["q_solution"]], dtype=np.float64)
            q_exact = np.asarray(archive[keys["q_exact"]], dtype=np.float64)
            force_saved = np.asarray(archive[keys["external_force"]], dtype=np.float64)
            gauge_saved = np.asarray(archive[keys["physical_gauge_target"]], dtype=np.float64)
            lambda_saved = np.asarray(archive["pure_I_lambda_physical"], dtype=np.float64)
        source = case["source_binding"]
        expected_ids, expected_xyz = tuple(source["source_ids"]), np.asarray(source["source_positions_xyz_m"], dtype=np.float64)
        ids_from_operator = operator.semantic_coefficient_ids(expected_ids)
        if (source_ids != expected_ids or source_positions.shape != (27, 3)
                or not np.array_equal(source_positions, expected_xyz)
                or semantic_ids != tuple(case["registry"]["semantic_coefficient_ids_in_order"])
                or semantic_ids != ids_from_operator):
            raise GateFailure("cache_identity", "cached source coordinates/IDs/semantic IDs differ from the exact frozen 93 registry")
        shapes = {"K": K.shape, "C": C.shape, "rigid": rigid_saved.shape, "q_solution": q_solution.shape,
                  "q_exact": q_exact.shape, "force": force_saved.shape, "gauge": gauge_saved.shape, "lambda": lambda_saved.shape}
        if shapes != {"K": (93, 93), "C": (6, 93), "rigid": (6, 93), "q_solution": (93,),
                      "q_exact": (93,), "force": (93,), "gauge": (6,), "lambda": (6,)}:
            raise GateFailure("cache_identity", f"cached array shapes changed: {shapes}")
        for name, array in (("K", K), ("C", C), ("rigid", rigid_saved), ("q_solution", q_solution),
                            ("q_exact", q_exact), ("force", force_saved), ("gauge", gauge_saved), ("lambda", lambda_saved)):
            finite(np, name, array)
        result["evidence"]["cache"] = {"arrays_sha256": hashes["accepted_arrays_npz"],
            "independent_audit_sha256": pins["accepted_cache_binding"]["audit_sha256"],
            "recovery_receipt_sha256": pins["accepted_cache_binding"]["receipt_sha256"],
            "actual_array_keys": cache["actual_array_keys"], "source_ids": list(source_ids),
            "source_positions_xyz_m": source_positions.tolist(), "semantic_ids": list(semantic_ids),
            "K_raw_N_per_m": K.tolist(), "C_phys": C.tolist(), "q_solution_KI100_m": q_solution.tolist(),
            "q_exact_KI100_m_diagnostic_only": q_exact.tolist(), "f_ext_KI100_N": force_saved.tolist(),
            "g_phys_KI100_m": gauge_saved.tolist(), "lambda_phys_KI100": lambda_saved.tolist(),
            "rigid_coefficient_modes": rigid_saved.tolist(),
            "q0_KI1_m": (q_solution * float(cache["linear_scale_from_cached_KI100"])).tolist(),
            "q_exact_KI1_m_diagnostic_only": (q_exact * float(cache["linear_scale_from_cached_KI100"])).tolist(),
            "f_ext_KI1_N": (force_saved * float(cache["linear_scale_from_cached_KI100"])).tolist(),
            "g_phys_KI1_m": (gauge_saved * float(cache["linear_scale_from_cached_KI100"])).tolist(),
            "lambda_phys_KI1": (lambda_saved * float(cache["linear_scale_from_cached_KI100"])).tolist()}
        loaded_bytes = sum(int(a.nbytes) for a in (source_positions, K, C, rigid_saved, q_solution,
            q_exact, force_saved, gauge_saved, lambda_saved))
        working_estimate = 2 * loaded_bytes + 2 * 16 * 27 * 4 * 8 + 2 * 16 * 3 * 93 * 8 + 3 * 93 * 93 * 8 + 8 * 1024 * 1024
        record_gate(gates, "resource_budget", working_estimate <= MEMORY_LIMIT,
                    {"conservative_working_memory_estimate_bytes": working_estimate,
                     "limit_bytes": MEMORY_LIMIT, "output_limit_bytes": OUTPUT_LIMIT})
        scale = float(cache["linear_scale_from_cached_KI100"])
        q0, force, gauge, multipliers = q_solution * scale, force_saved * scale, gauge_saved * scale, lambda_saved * scale
        weak0 = K @ q0 + C.T @ multipliers - force
        gauge0 = C @ q0 - gauge
        weak0_rel = float(np.linalg.norm(weak0, ord=np.inf) / max(float(np.linalg.norm(K, ord=np.inf) * np.linalg.norm(q0, ord=np.inf)
            + np.linalg.norm(C.T @ multipliers, ord=np.inf) + np.linalg.norm(force, ord=np.inf)), np.finfo(float).tiny))
        gauge0_rel = float(np.linalg.norm(gauge0, ord=np.inf) / max(float(np.linalg.norm(C, ord=np.inf) * np.linalg.norm(q0, ord=np.inf)
            + np.linalg.norm(gauge, ord=np.inf)), np.finfo(float).tiny))
        record_gate(gates, "scaled_cached_kkt", max(weak0_rel, gauge0_rel) <= KKT_TOL,
                    {"weak_relative_error": weak0_rel, "gauge_relative_error": gauge0_rel,
                     "q_source": "actual pure_I_q_solution scaled by 0.01", "q_exact_used_as_solution": False,
                     "q_exact_diagnostic_norm_m": float(np.linalg.norm(q_exact * scale))})
        config_data = case["mls"]
        config = mls_module.MLS3DConfig(reference_origin=tuple(config_data["reference_origin_m"]),
            reference_scale=float(config_data["reference_scale_m"]), support_radius=float(config_data["support_radius_m"]),
            max_condition_number=float(config_data["max_condition_number"]),
            rank_relative_tolerance=float(config_data["rank_relative_tolerance"]),
            switch_tolerance=float(config_data["switch_tolerance"]), duplicate_tolerance=float(config_data["duplicate_tolerance"]))
        qr = case["interface_quadrature"]; radius, length = case["geometry"]["radius_m"], case["geometry"]["length_m"]
        queries, areas, coordinates = [], [], []
        for ix, node_x in enumerate(qr["gl4_nodes"]):
            for iz, node_z in enumerate(qr["gl4_nodes"]):
                x, z = radius / 2.0 * (node_x - 1.0), length / 2.0 * node_z
                queries.append((x, 0.0, z)); coordinates.append((x, z))
                areas.append(radius * length / 4.0 * qr["gl4_weights"][ix] * qr["gl4_weights"][iz])
        queries = np.asarray(queries, dtype=np.float64); areas = np.asarray(areas, dtype=np.float64)
        result["counts"]["mls_queries_requested"] = len(queries); result["counts"]["mls_query_attempts"] += 1
        result["counts"]["new_mls_queries"] = int(len(queries))
        shape = shape_module.shape_functions_3d(source_positions, queries, config, visibility_mask=None)
        base = shape.weights
        finite(np, "MLS weights", base); finite(np, "MLS gradients", shape.gradients)
        value_moment = base @ source_positions
        grad_scaled = shape.gradients * config.reference_scale
        source_scaled = (source_positions - np.asarray(config.reference_origin)) / config.reference_scale
        grad_sum = np.sum(grad_scaled, axis=1)
        grad_moment = np.einsum("ni,qnj->qij", source_scaled, grad_scaled)
        value_error = float(np.max(np.abs(value_moment - queries)) / config.reference_scale)
        partition_error = float(np.max(np.abs(np.sum(base, axis=1) - 1.0)))
        grad_sum_error = float(np.max(np.abs(grad_sum)))
        identity_error = float(np.max(np.abs(grad_moment - np.eye(3))))
        statuses = [getattr(s, "value", str(s)) for s in shape.statuses]
        mls_pass = (len(statuses) == 16 and all(s == "OK" for s in statuses)
                    and bool(np.all(shape.support_counts == 27))
                    and bool(np.all(shape.condition_numbers <= config.max_condition_number))
                    and max(value_error, partition_error, grad_sum_error, identity_error) <= 1.0e-10)
        record_gate(gates, "mls", mls_pass, {"statuses": statuses, "support_counts": shape.support_counts.tolist(),
                    "condition_numbers": shape.condition_numbers.tolist(), "normalized_value_moment_error": value_error,
                    "partition_error": partition_error, "normalized_gradient_sum_error": grad_sum_error,
                    "normalized_gradient_moment_error": identity_error, "visibility_mask": None})
        result["evidence"]["mls"] = {"query_xyz_m": queries.tolist(), "weights_query_source": base.tolist(),
            "gradients_query_source_xyz_per_m": shape.gradients.tolist(), "statuses": statuses,
            "support_counts": shape.support_counts.tolist(), "condition_numbers": shape.condition_numbers.tolist(),
            "normalized_value_moment_error": value_error, "partition_error": partition_error,
            "normalized_gradient_sum_error": grad_sum_error, "normalized_gradient_moment_error": identity_error}
        slots = {identifier: i for i, identifier in enumerate(semantic_ids)}
        candidates = operator.rigid_body_mode_coefficients(source_ids, tuple(map(tuple, source_positions)))
        rigid_modes = np.frombuffer(candidates.values_le_f64, dtype="<f8").copy().reshape(6, 93)
        points, Hplus, Hminus, jump_ops = [], [], [], []
        exact_trace_error, other_branch_max = 0.0, 0.0
        for qi, (x, z) in enumerate(coordinates):
            upper = field_module.branch_functions_local(x, 0.0, reference_length_m=radius, branch_trace="upper").values
            lower = field_module.branch_functions_local(x, 0.0, reference_length_m=radius, branch_trace="lower").values
            analytic = math.sqrt(-x / radius)
            exact_trace_error = max(exact_trace_error, abs(upper[0] - analytic), abs(lower[0] + analytic))
            other_branch_max = max(other_branch_max, *(abs(v) for v in (*upper[1:], *lower[1:])))
            if exact_trace_error > 1.0e-14 or other_branch_max > 1.0e-14:
                raise GateFailure("exact_traces", {"F1_error": exact_trace_error, "F2_F4_roundoff": other_branch_max})
            up, lo = (analytic, 0.0, 0.0, 0.0), (-analytic, 0.0, 0.0, 0.0)
            plus, minus = np.zeros((3, 93)), np.zeros((3, 93))
            for ni, sid in enumerate(source_ids):
                for ci, component in enumerate(("ux", "uy", "uz")):
                    col = slots[json.dumps(["SOURCE_BASE", sid, component, 0], separators=(",", ":"))]
                    plus[ci, col] = minus[ci, col] = base[qi, ni]
            for branch in range(1, 5):
                for ci, component in enumerate(("ux", "uy", "uz")):
                    col = slots[json.dumps(["GLOBAL_FRONT", "global", component, branch], separators=(",", ":"))]
                    plus[ci, col], minus[ci, col] = up[branch - 1], lo[branch - 1]
            point_xyz = (float(x), 0.0, float(z))
            center = np.mean(source_positions, axis=0)
            expected = []
            for mode in range(6):
                if mode < 3:
                    expected.append(tuple(1.0 if component == mode else 0.0 for component in range(3)))
                else:
                    omega = np.eye(3)[mode - 3]
                    expected.append(tuple(np.cross(omega, np.asarray(point_xyz) - center)))
            ep = trace.EnrichedFacePoint("finite_sheet_mode_i_20261006", f"qp-{qi:02d}", trace.TraceFamily.GLOBALFRONT,
                semantic_ids, tuple(map(tuple, plus)), tuple(map(tuple, minus)), (0.0, 1.0, 0.0), point_xyz,
                point_xyz, tuple(expected), float(areas[qi]))
            points.append(ep); Hplus.append(plus); Hminus.append(minus); jump_ops.append(plus - minus)
        result["counts"]["interface_qps"] = len(points)
        result["evidence"]["interface"] = [{"quadrature_id": point.quadrature_id,
            "point_xyz_m": list(point.positive_reference_point_m), "normal_unit_lower_to_upper": list(point.normal),
            "one_side_area_weight_m2": point.one_side_area_weight_m2, "Hplus": Hplus[i].tolist(),
            "Hminus": Hminus[i].tolist(), "J": jump_ops[i].tolist(),
            "expected_rigid_displacements_m": [list(v) for v in point.expected_rigid_displacements_m]}
            for i, point in enumerate(points)]
        record_gate(gates, "exact_traces", exact_trace_error <= 1.0e-14 and other_branch_max <= 1.0e-14,
                    {"query_count": len(points), "F1_upper_lower_formula_max_error": exact_trace_error,
                     "F2_F4_raw_roundoff_max_before_exact_limit": other_branch_max,
                     "base_trace_identical": True, "only_F1_branch_jump": True, "front_avoided_by_open_GL_nodes": True})
        candidates = operator.rigid_body_mode_coefficients(source_ids, tuple(map(tuple, source_positions)))
        rigid_modes = np.frombuffer(candidates.values_le_f64, dtype="<f8").copy().reshape(6, 93)
        rigid_singular = np.linalg.svd(rigid_modes, compute_uv=False)
        rigid_rank = int(np.count_nonzero(rigid_singular > 128.0 * np.finfo(float).eps * 93 * rigid_singular[0]))
        rigid_errors = []
        for p, hp, hm, jump in zip(points, Hplus, Hminus, jump_ops):
            expected = np.asarray(p.expected_rigid_displacements_m, dtype=np.float64).T
            rigid_errors.extend((float(np.max(np.abs(hp @ rigid_modes.T - expected))),
                                 float(np.max(np.abs(hm @ rigid_modes.T - expected))),
                                 float(np.max(np.abs(jump @ rigid_modes.T)))))
        record_gate(gates, "rigid_traces", rigid_rank == 6 and rigid_saved.shape == (6, 93)
                    and bool(np.allclose(rigid_saved, rigid_modes, rtol=0.0, atol=1.0e-15))
                    and max(rigid_errors) <= 1.0e-10,
                    {"mode_order": list(candidates.mode_ids), "rank": rigid_rank,
                     "saved_mode_max_difference": float(np.max(np.abs(rigid_saved - rigid_modes))),
                     "max_Hplus_Hminus_or_jump_error_m": max(rigid_errors),
                     "translation_coefficients_are_unit_m": True, "rotation_coefficients_are_omega_cross_X_m": True})
        record_gate(gates, "interface_quadrature", len(points) == 16
                    and abs(float(np.sum(areas) - radius * length)) <= 1.0e-12 * radius * length,
                    {"qp_count": len(points), "one_side_area_sum_m2": float(np.sum(areas)), "area_counted_once": True})
        mat_data = case["cohesive_law"]["material"]
        material = trace.BilinearCohesiveMaterial(float(mat_data["k0_pa_per_m"]), float(mat_data["tpeak_pa"]), float(mat_data["gc_j_per_m2"]))
        virgin = {point.point_id: scalar.initial_state(material) for point in points}
        rigid_candidates = trace.RigidTraceCandidates(semantic_ids, tuple(map(tuple, rigid_modes)))
        result["counts"]["primitive_assemblies"] += 1
        baseline = trace.assemble_enriched_trial(material, points, tuple(q0), virgin, rigid_candidates)
        result["counts"]["primitive_successful_returns"] += 1
        G = np.zeros((93, 93), dtype=np.float64); G[82, 82] = 2.0 * material.k0_pa_per_m * radius * length
        fbase = np.asarray(baseline.generalized_internal_force_n); Kbase = np.asarray(baseline.tangent_n_per_m)
        primitive_G_error = relative_error(np, Kbase, G)
        primitive_f_error = relative_error(np, fbase, G @ q0)
        openings = np.asarray([item.opening_m for item in baseline.point_results])
        baseline_pass = (abs(G[82, 82] - 20000.0) <= 1.0e-12 * 20000.0 and np.count_nonzero(G) == 1 and primitive_G_error <= 1.0e-12
                         and primitive_f_error <= 1.0e-12 and bool(np.all(openings > 0.0))
                         and bool(np.all(openings < material.delta0_m))
                         and baseline.tangent_status.value == "differentiable")
        if not baseline_pass:
            raise GateFailure("primitive_baseline", {"KG_error": primitive_G_error, "force_error": primitive_f_error,
                "opening_min_m": float(np.min(openings)), "opening_max_m": float(np.max(openings)), "delta0_m": material.delta0_m})
        physical_upper_traction = [-float(item.traction_pa) * np.asarray(points[i].normal) for i, item in enumerate(baseline.point_results)]
        internal_face = [np.asarray(item.positive_face_force_at_reference_n) for item in baseline.point_results]
        sign_positive = all(float(np.dot(points[i].normal, internal_face[i])) > 0.0
                            and float(np.dot(points[i].normal, physical_upper_traction[i])) < 0.0
                            and float(np.linalg.norm(internal_face[i] - points[i].one_side_area_weight_m2
                                * baseline.point_results[i].traction_pa * np.asarray(points[i].normal))) <= 1.0e-12
                            and float(np.linalg.norm(points[i].one_side_area_weight_m2 * physical_upper_traction[i]
                                + internal_face[i])) <= 1.0e-12
                            for i in range(16))
        qf = float(q0 @ fbase); sum_area_open_t = float(sum(points[i].one_side_area_weight_m2 * baseline.point_results[i].opening_m * baseline.point_results[i].traction_pa for i in range(16)))
        baseline_energy = 0.5 * float(q0 @ G @ q0)
        record_gate(gates, "primitive_baseline", sign_positive and qf > 0.0 and abs(qf - sum_area_open_t) <= 1.0e-12 * max(abs(qf), 1e-30)
                    and abs(float(baseline.stored_energy_j) - baseline_energy) <= 1.0e-12 * max(abs(baseline_energy), 1e-30)
                    and float(baseline.dissipated_energy_j) == 0.0,
                    {"KG_only_nonzero_index": [82, 82], "KG_N_per_m": float(G[82, 82]), "tangent_relative_error": primitive_G_error,
                     "force_relative_error": primitive_f_error, "opening_min_m": float(np.min(openings)),
                     "opening_max_m": float(np.max(openings)), "q_dot_internal_force_j": qf,
                     "sum_area_opening_traction_j": sum_area_open_t, "stored_energy_j": float(baseline.stored_energy_j),
                     "analytic_stored_energy_j": baseline_energy, "upper_physical_traction_is_minus_tn": True,
                     "primitive_positive_face_force_is_internal_plus_Atn": True, "action_reaction_alone_is_not_gate": True})
        result["evidence"]["baseline"] = {"generalized_internal_force_N": fbase.tolist(), "tangent_N_per_m": Kbase.tolist(),
            "stored_energy_j": float(baseline.stored_energy_j), "dissipated_energy_j": float(baseline.dissipated_energy_j),
            "primitive_external_work_path_integral_j": float(baseline.external_work_j), "work_interpretation": "not_incremental_work",
            "points": point_ledger(baseline, points, tuple(q0))}
        full_state = scalar.update_state(material, scalar.initial_state(material), material.deltaf_m)
        zero_history = {point.point_id: full_state for point in points}
        result["counts"]["primitive_assemblies"] += 1
        zero = trace.assemble_enriched_trial(material, points, tuple(q0), zero_history, rigid_candidates)
        result["counts"]["primitive_successful_returns"] += 1
        zero_force, zero_tangent = np.asarray(zero.generalized_internal_force_n), np.asarray(zero.tangent_n_per_m)
        record_gate(gates, "zero_tangent_control", float(np.linalg.norm(zero_force)) == 0.0 and float(np.linalg.norm(zero_tangent)) == 0.0,
                    {"force_norm_N": float(np.linalg.norm(zero_force)), "tangent_norm_N_per_m": float(np.linalg.norm(zero_tangent)),
                     "synthetic_fully_damaged_history": True, "state_discarded": True, "physical_solve_input": False})
        raw_spec = spectrum_gate(np, K, "K_raw")
        record_gate(gates, "raw_spectrum", raw_spec["passed"], raw_spec)
        A = K + G
        coupled_spec = spectrum_gate(np, A, "K_raw_plus_KG")
        record_gate(gates, "coupled_spectrum", coupled_spec["passed"], coupled_spec)
        d = 1.0 / np.sqrt(np.diag(A)); As = d[:, None] * A * d[None, :]
        Cbar = C * d[None, :]; row_norms = np.linalg.norm(Cbar, axis=1)
        if bool(np.any(row_norms <= np.finfo(float).tiny)):
            raise GateFailure("coupled_kkt", "physical gauge contains an empty scaled row")
        row_scale = 1.0 / row_norms; Chat = row_scale[:, None] * Cbar; ghat = row_scale * gauge
        saddle = np.zeros((99, 99), dtype=np.float64); saddle[:93, :93] = As
        saddle[:93, 93:] = Chat.T; saddle[93:, :93] = Chat
        rhs = np.concatenate((d * force, ghat))
        result["counts"]["coupled_kkt_solves"] += 1
        solution = np.linalg.solve(saddle, rhs)
        result["counts"]["coupled_kkt_successful_returns"] += 1
        y, lambda_scaled = solution[:93], solution[93:]
        qc, lambda_c = d * y, row_scale * lambda_scaled
        kkt_rel = float(np.linalg.norm(saddle @ solution - rhs, ord=np.inf) / max(float(np.linalg.norm(saddle, ord=np.inf) * np.linalg.norm(solution, ord=np.inf)
            + np.linalg.norm(rhs, ord=np.inf)), np.finfo(float).tiny))
        weak = A @ qc + C.T @ lambda_c - force
        weak_rel = float(np.linalg.norm(weak, ord=np.inf) / max(float(np.linalg.norm(A, ord=np.inf) * np.linalg.norm(qc, ord=np.inf)
            + np.linalg.norm(C.T @ lambda_c, ord=np.inf) + np.linalg.norm(force, ord=np.inf)), np.finfo(float).tiny))
        gauge_rel = float(np.linalg.norm(C @ qc - gauge, ord=np.inf) / max(float(np.linalg.norm(C, ord=np.inf) * np.linalg.norm(qc, ord=np.inf)
            + np.linalg.norm(gauge, ord=np.inf)), np.finfo(float).tiny))
        ubulk, ugamma = 0.5 * float(qc @ K @ qc), 0.5 * float(qc @ G @ qc)
        force_pairing, gauge_pairing = float(qc @ force), float(gauge @ lambda_c)
        energy_rel = abs(2.0 * (ubulk + ugamma) - force_pairing + gauge_pairing) / max(abs(2.0 * (ubulk + ugamma)), abs(force_pairing), abs(gauge_pairing), np.finfo(float).tiny)
        record_gate(gates, "coupled_kkt", max(kkt_rel, weak_rel, gauge_rel, energy_rel) <= KKT_TOL,
                    {"order": 99, "backward_error": kkt_rel, "weak_relative_error": weak_rel,
                     "gauge_relative_error": gauge_rel, "energy_identity_relative_error": energy_rel})
        result["counts"]["primitive_assemblies"] += 1
        final = trace.assemble_enriched_trial(material, points, tuple(qc), virgin, rigid_candidates)
        result["counts"]["primitive_successful_returns"] += 1
        final_force, final_tangent = np.asarray(final.generalized_internal_force_n), np.asarray(final.tangent_n_per_m)
        final_openings = np.asarray([item.opening_m for item in final.point_results])
        final_pass = (relative_error(np, final_tangent, G) <= 1.0e-12 and relative_error(np, final_force, G @ qc) <= 1.0e-12
                      and bool(np.all(final_openings > 0.0)) and bool(np.all(final_openings < material.delta0_m))
                      and final.tangent_status.value == "differentiable")
        record_gate(gates, "primitive_final", final_pass, {"tangent_relative_error": relative_error(np, final_tangent, G),
                    "force_relative_error": relative_error(np, final_force, G @ qc), "opening_min_m": float(np.min(final_openings)),
                    "opening_max_m": float(np.max(final_openings)), "dissipated_energy_j": float(final.dissipated_energy_j)})
        delta_q = float(np.linalg.norm(qc - q0)); cohesive_force = float(np.linalg.norm(final_force))
        q_floor = 128.0 * np.finfo(float).eps * float(np.linalg.norm(q0))
        f_floor = 128.0 * np.finfo(float).eps * float(np.linalg.norm(force))
        record_gate(gates, "mechanical_change", delta_q > q_floor and cohesive_force > f_floor and float(np.linalg.norm(G @ qc)) > f_floor,
                    {"q_change_norm_m": delta_q, "q_roundoff_floor_m": q_floor, "cohesive_force_norm_n": cohesive_force,
                     "force_roundoff_floor_n": f_floor})
        record_gate(gates, "energy", abs(float(final.stored_energy_j) - ugamma) <= 1.0e-12 * max(abs(ugamma), 1e-30)
                    and float(final.dissipated_energy_j) == 0.0
                    and abs(float(final.external_work_j) - ugamma) <= 1.0e-12 * max(abs(ugamma), 1e-30),
                    {"bulk_stored_j": ubulk, "cohesive_stored_j": ugamma, "primitive_stored_j": float(final.stored_energy_j),
                     "primitive_dissipated_j": float(final.dissipated_energy_j), "primitive_external_work_j": float(final.external_work_j),
                     "primitive_external_work_is_path_integral_equal_to_stored_on_virgin_elastic_ramp": True,
                     "global_energy_identity_relative_error": energy_rel})
        result["evidence"]["analytic_KG_N_per_m"] = G.tolist()
        result["evidence"]["coupled"] = {"K_plus_KG_N_per_m": A.tolist(), "f_ext_N": force.tolist(),
            "g_phys_m": gauge.tolist(), "q_c_m": qc.tolist(), "lambda_c": lambda_c.tolist(),
            "generalized_internal_force_N": final_force.tolist(), "tangent_N_per_m": final_tangent.tolist(),
            "U_bulk_j": ubulk, "U_gamma_j": ugamma,
            "external_force_displacement_pairing_j": force_pairing, "gauge_multiplier_pairing_j": gauge_pairing,
            "energy_identity_interpretation": "algebraic_stationarity_pairing; not incremental_work_or_loading_path_assessment",
            "primitive_external_work_path_integral_j": float(final.external_work_j),
            "primitive_external_work_interpretation": "virgin cohesive loading-path integral; not global f_dot_q pairing",
            "points": point_ledger(final, points, tuple(qc))}
        if any(gate["passed"] is not True for gate in gates.values()): raise GateFailure("all_gates", "one or more numerical gates remain false or UNKNOWN")
        result.update(status="ACCEPTED_DIAGNOSTIC_ONLY", numeric_modules_imported=True, npz_opened=True,
            inputs={"npz_sha256": hashes["accepted_arrays_npz"], "source_ids": list(source_ids),
                    "q_solution_scale": scale, "q_exact_kept_diagnostic_only": True},
            solution={"q0_m": q0.tolist(), "q_c_m": qc.tolist(), "lambda_c": lambda_c.tolist(),
                      "q_change_norm_m": delta_q, "internal_force_norm_n": cohesive_force},
            energy={"U_bulk_j": ubulk, "U_gamma_j": ugamma,
                    "external_force_displacement_pairing_j": force_pairing, "gauge_multiplier_pairing_j": gauge_pairing,
                    "identity_is_not_incremental_work_or_loading_path_assessment": True},
            limitations=["toy fixed-crack elastic cohesive coupling only", "no material calibration, process-zone, mixed-mode, or front-growth qualification"])
        result["gates"] = gates
        return result, 0
    except GateFailure as exc:
        gates.setdefault(exc.gate, {"status": "FAIL", "passed": False, "evidence": str(exc.detail)})
        gates[exc.gate] = {"status": "FAIL", "passed": False, "evidence": exc.detail}
        result["failure"] = {"gate": exc.gate, "detail": exc.detail}
    except Exception as exc:
        result["failure"] = {"gate": "execution", "type": type(exc).__name__, "detail": str(exc)}
    result["gates"] = gates
    result["status"] = "NOT_ACCEPTED"
    return result, 3


def json_default(value):
    if type(value).__module__.split(".", 1)[0] == "numpy" and hasattr(value, "item"):
        return value.item()
    raise TypeError(f"not JSON serializable: {type(value).__name__}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="root-go gated numerical 93-DOF cohesive bridge worker")
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--root-go")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    if not args.execute or args.output is None:
        print(json.dumps({"status": "EXECUTION_DISABLED", "numeric_modules_imported": False, "npz_opened": False}, allow_nan=False))
        return 2
    try:
        case, pins = json.loads(CASE_PATH.read_text(encoding="utf-8")), json.loads(PINS_PATH.read_text(encoding="utf-8"))
        require_execution(pins, args.root_go)
    except Exception as exc:
        print(json.dumps({"status": "NOT_ACCEPTED", "failure": str(exc), "numeric_modules_imported": False, "npz_opened": False}, allow_nan=False))
        return 2
    result, code = execute(case, pins, args.output)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(result, indent=2, sort_keys=True, allow_nan=False, default=json_default) + "\n")
    print(json.dumps({"status": result["status"], "result": args.output.name,
                      "counts": result.get("counts", {}), "failure": result.get("failure")}, allow_nan=False))
    return code


if __name__ == "__main__":
    raise SystemExit(main())
