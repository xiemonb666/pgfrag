"""Root-go-gated fixed-face Mode-I opening-history diagnostic for the saved 93-DOF case."""
from __future__ import annotations

import argparse
import base64
import hashlib
import importlib.util
import json
import math
import sys
import time
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
REVIEW = ROOT / "pgfrag" / "tmp" / "review_20261007" / "mixed93_opening_control_20261007"
CASE_PATH = HERE / "case.json"
PINS_PATH = REVIEW / "SOURCE_PINS.json"
SCHEMA = "mixed93_opening_control_result.v1"
KKT_TOL = 1.0e-9
TRACE_ABS_TOL_M = 1.0e-12
TRACE_REL_TOL = 1.0e-10
ANALYTIC_NORMAL_ABS_TOL = 1.0e-14
FORCE_ABS_TOL_N = 1.0e-12
FORCE_REL_TOL = 1.0e-10
UNIT_FORCE_ABS_TOL_N_PER_M = 1.0e-12
ENERGY_REL_TOL = 2.0e-10
WALL_LIMIT_S = 20.0
OUTPUT_LIMIT_BYTES = 8 * 1024 * 1024
CONTROL_INDEX = 82
NDof = 93


class GateFailure(RuntimeError):
    def __init__(self, gate: str, detail: object):
        super().__init__(f"{gate}: {detail}")
        self.gate, self.detail = gate, detail


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_execution(pins: dict, token: str | None) -> None:
    go = pins.get("root_go")
    if (pins.get("execution_authorized") is not True or not isinstance(go, dict)
            or go.get("authorized") is not True or not isinstance(go.get("token"), str)
            or not go["token"] or token != go["token"]):
        raise PermissionError("execution is disabled until SOURCE_PINS.json has an exact root-go token")


def check_pins(case: dict, pins: dict) -> dict[str, str]:
    frozen = pins["frozen_contract"]
    if sha256(CASE_PATH) != frozen["case_sha256"]:
        raise GateFailure("source_pins", "case.json changed after source freeze")
    if sha256(Path(__file__)) != frozen["driver_sha256"]:
        raise GateFailure("source_pins", "driver changed after source freeze")
    hashes = {}
    for name, item in pins["references"].items():
        path = ROOT / Path(item["path"])
        actual = sha256(path)
        if actual != item["sha256"]:
            raise GateFailure("source_pins", {"reference": name, "actual_sha256": actual})
        hashes[name] = actual
    source = case["saved_result"]
    if hashes.get("saved_elastic_result") != source["sha256"]:
        raise GateFailure("source_pins", "case and pins disagree on saved result identity")
    return hashes


def load_registered(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot import pinned primitive {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def finite(np, name: str, value) -> None:
    if not bool(np.isfinite(value).all()):
        raise GateFailure("finite_values", f"non-finite {name}")


def close_enough(actual: float, expected: float, atol: float, rtol: float) -> bool:
    return abs(actual - expected) <= atol + rtol * max(abs(actual), abs(expected))


def envelope_traction(material, opening_m: float) -> float:
    if opening_m <= material.delta0_m:
        return material.k0_pa_per_m * opening_m
    if opening_m < material.deltaf_m:
        return material.tpeak_pa * (material.deltaf_m - opening_m) / (material.deltaf_m - material.delta0_m)
    return 0.0


def envelope_work(material, start_m: float, end_m: float) -> float:
    """Independent exact integral of the declared piecewise-linear virgin envelope."""
    if end_m <= start_m:
        return 0.0
    cuts = [start_m]
    cuts.extend(c for c in (material.delta0_m, material.deltaf_m) if start_m < c < end_m)
    cuts.append(end_m)
    return sum(0.5 * (envelope_traction(material, a) + envelope_traction(material, b)) * (b - a)
               for a, b in zip(cuts, cuts[1:]))


def law_interval_work(material, old_state, new_opening_m: float) -> float:
    """Integrate the bilinear law along this step using only committed history and opening."""
    stiffness = material.k0_pa_per_m * (1.0 - old_state.damage)
    if new_opening_m <= old_state.kappa_m:
        return 0.5 * stiffness * (new_opening_m**2 - old_state.opening_m**2)
    reload = 0.5 * stiffness * (old_state.kappa_m**2 - old_state.opening_m**2)
    return reload + envelope_work(material, old_state.kappa_m, new_opening_m)


def relative_error(np, actual, expected) -> float:
    return float(np.linalg.norm(actual - expected, ord=np.inf) /
                 max(float(np.linalg.norm(actual, ord=np.inf)), float(np.linalg.norm(expected, ord=np.inf)), 1.0e-300))


def execute(case: dict, pins: dict) -> dict:
    elapsed_start = time.monotonic()
    hashes = check_pins(case, pins)
    import numpy as np

    cache_path = ROOT / Path(case["saved_result"]["path"])
    saved = json.loads(cache_path.read_text(encoding="utf-8"))
    if saved.get("schema") != "mixed93_actual_cohesive_bridge_numeric_result.v1" or saved.get("status") != "ACCEPTED_DIAGNOSTIC_ONLY":
        raise GateFailure("saved_result", "saved result is not the pinned accepted elastic diagnostic")
    if saved.get("npz_opened") is not True:
        raise GateFailure("saved_result", "saved elastic run did not contain the previously accepted matrices")
    evidence = saved["evidence"]
    cached = evidence["cache"]
    K = np.asarray(cached["K_raw_N_per_m"], dtype=np.float64)
    C = np.asarray(cached["C_phys"], dtype=np.float64)
    rigid = np.asarray(cached["rigid_coefficient_modes"], dtype=np.float64)
    source_positions = np.asarray(cached["source_positions_xyz_m"], dtype=np.float64)
    coefficient_ids = tuple(cached["semantic_ids"])
    if K.shape != (NDof, NDof) or C.shape != (6, NDof) or rigid.shape != (6, NDof):
        raise GateFailure("saved_result", "K, C, or rigid-mode array shape changed")
    if len(coefficient_ids) != NDof or source_positions.shape != (27, 3):
        raise GateFailure("saved_result", "semantic coefficient or source-node layout changed")
    if coefficient_ids[CONTROL_INDEX] != case["control"]["semantic_id"]:
        raise GateFailure("control_identity", coefficient_ids[CONTROL_INDEX])
    for name, array in (("K_raw", K), ("C_phys", C), ("rigid_modes", rigid), ("source_positions", source_positions)):
        finite(np, name, array)

    # Reconstruct J_i from saved H+_i/H-_i and check its exact normal row analytically.
    trace = load_registered("_opening_control_pinned_trace", ROOT / Path(pins["references"]["trace_primitive"]["path"]))
    scalar = trace._DERIVATIVE._SCALAR
    mat_data = case["cohesive_law"]
    material = trace.BilinearCohesiveMaterial(float(mat_data["k0_pa_per_m"]),
        float(mat_data["tpeak_pa"]), float(mat_data["gc_j_per_m2"]))
    normal = np.asarray(case["geometry"]["normal_lower_to_upper"], dtype=np.float64)
    if not np.array_equal(normal, np.asarray((0.0, 1.0, 0.0))):
        raise GateFailure("geometry", "frozen interface normal changed")
    interfaces, points, alpha = [], [], []
    for row in evidence["interface"]:
        hp = np.asarray(row["Hplus"], dtype=np.float64)
        hm = np.asarray(row["Hminus"], dtype=np.float64)
        saved_jump = np.asarray(row["J"], dtype=np.float64)
        if hp.shape != (3, NDof) or hm.shape != (3, NDof) or saved_jump.shape != (3, NDof):
            raise GateFailure("interface_reconstruction", row.get("quadrature_id"))
        jump = hp - hm
        x = float(row["point_xyz_m"][0])
        if not x < 0.0:
            raise GateFailure("interface_reconstruction", "quadrature point touched or crossed the front")
        normal_factor = 2.0 * math.sqrt(-x / float(case["geometry"]["radius_m"]))
        normal_row = normal @ jump
        expected_row = np.zeros(NDof, dtype=np.float64)
        expected_row[CONTROL_INDEX] = normal_factor
        if not np.array_equal(jump, saved_jump) or not np.allclose(normal_row, expected_row, rtol=0.0, atol=ANALYTIC_NORMAL_ABS_TOL):
            raise GateFailure("analytic_normal_trace", {"quadrature_id": row["quadrature_id"],
                "normal_row_max_error": float(np.max(np.abs(normal_row - expected_row)))})
        point_xyz = tuple(float(v) for v in row["point_xyz_m"])
        expected_rigid = tuple(tuple(float(v) for v in v3) for v3 in row["expected_rigid_displacements_m"])
        point = trace.EnrichedFacePoint("finite_sheet_mode_i_20261006", row["quadrature_id"],
            trace.TraceFamily.GLOBALFRONT, coefficient_ids, tuple(map(tuple, hp)), tuple(map(tuple, hm)),
            tuple(float(v) for v in normal), point_xyz, point_xyz, expected_rigid,
            float(row["one_side_area_weight_m2"]))
        points.append(point)
        alpha.append(normal_factor)
        interfaces.append({"quadrature_id": row["quadrature_id"], "point_xyz_m": list(point_xyz),
            "one_side_area_weight_m2": point.one_side_area_weight_m2, "normal_lower_to_upper": list(normal),
            "Hplus": hp.tolist(), "Hminus": hm.tolist(), "J_reconstructed": jump.tolist(),
            "normal_row": normal_row.tolist(), "analytic_alpha": normal_factor,
            "expected_rigid_displacements_m": [list(v) for v in expected_rigid]})
    if len(points) != 16 or not close_enough(sum(p.one_side_area_weight_m2 for p in points),
            float(case["geometry"]["radius_m"]) * float(case["geometry"]["length_m"]), 1.0e-15, 1.0e-12):
        raise GateFailure("interface_quadrature", "expected exactly 16 one-sided weights summing to R*L")
    rigid_modes = tuple(tuple(float(v) for v in row) for row in rigid)
    rigid_candidates = trace.RigidTraceCandidates(coefficient_ids, rigid_modes)

    # One diagonally equilibrated, 100-by-100 bulk/gauge/opening KKT solve.
    diagonal = np.diag(K)
    if bool(np.any(diagonal <= 0.0)):
        raise GateFailure("unit_control_kkt", "K has a nonpositive diagonal")
    d = 1.0 / np.sqrt(diagonal)
    Kbar = d[:, None] * K * d[None, :]
    Cbar_raw = C * d[None, :]
    gauge_norms = np.linalg.norm(Cbar_raw, axis=1)
    if bool(np.any(gauge_norms <= np.finfo(float).tiny)):
        raise GateFailure("unit_control_kkt", "physical gauge contains an empty scaled row")
    gauge_scale = 1.0 / gauge_norms
    Cbar = gauge_scale[:, None] * Cbar_raw
    d_control = float(d[CONTROL_INDEX])
    unit_e = np.zeros(NDof, dtype=np.float64)
    unit_e[CONTROL_INDEX] = 1.0
    saddle = np.zeros((100, 100), dtype=np.float64)
    saddle[:NDof, :NDof] = Kbar
    saddle[:NDof, NDof:NDof+6] = Cbar.T
    saddle[NDof:NDof+6, :NDof] = Cbar
    saddle[:NDof, -1] = unit_e
    saddle[-1, :NDof] = unit_e
    rhs = np.zeros(100, dtype=np.float64)
    rhs[-1] = 1.0 / d_control
    solution = np.linalg.solve(saddle, rhs)
    y = solution[:NDof]
    unit_q = d * y
    unit_lambda_phys = gauge_scale * solution[NDof:NDof+6]
    unit_saddle_control_multiplier = solution[-1] / d_control
    kkt_error = float(np.linalg.norm(saddle @ solution - rhs, ord=np.inf) /
                      max(float(np.linalg.norm(saddle, ord=np.inf) * np.linalg.norm(solution, ord=np.inf)
                          + np.linalg.norm(rhs, ord=np.inf)), np.finfo(float).tiny))
    bulk_unit_residual = K @ unit_q + C.T @ unit_lambda_phys
    noncontrol_unit = bulk_unit_residual.copy()
    noncontrol_unit[CONTROL_INDEX] = 0.0
    if (kkt_error > KKT_TOL or not close_enough(float(unit_q[CONTROL_INDEX]), 1.0, 1.0e-12, 1.0e-10)
            or float(np.linalg.norm(C @ unit_q, ord=np.inf)) > 1.0e-9
            or float(np.linalg.norm(noncontrol_unit, ord=np.inf)) > UNIT_FORCE_ABS_TOL_N_PER_M
                + FORCE_REL_TOL * max(abs(float(bulk_unit_residual[CONTROL_INDEX])), 1.0)):
        raise GateFailure("unit_control_kkt", {"backward_error": kkt_error,
            "control_value": float(unit_q[CONTROL_INDEX]), "gauge_residual_inf": float(np.linalg.norm(C @ unit_q, ord=np.inf)),
            "noncontrol_bulk_residual_inf": float(np.linalg.norm(noncontrol_unit, ord=np.inf))})
    for point, factor in zip(points, alpha):
        jump = np.asarray(point.jump_operator, dtype=np.float64) @ unit_q
        opening = float(np.dot(normal, jump))
        tangent = jump - opening * normal
        if (not close_enough(opening, factor, TRACE_ABS_TOL_M, TRACE_REL_TOL)
                or float(np.linalg.norm(tangent)) > TRACE_ABS_TOL_M + TRACE_REL_TOL * float(np.linalg.norm(jump))):
            raise GateFailure("unit_control_trace", {"quadrature_id": point.quadrature_id,
                "opening_per_Q": opening, "analytic_alpha": factor, "tangential_jump_per_Q": tangent.tolist(),
                "decision": "FAIL; no projection or trace repair is permitted"})

    q_star = material.deltaf_m / min(alpha)
    fractions = tuple(float(v) for v in case["trajectory_fraction_Qstar"])
    if len(fractions) != 12 or fractions[0] != 0.0 or max(fractions) < 1.0:
        raise GateFailure("trajectory", "frozen 12-step opening path changed")
    committed = {p.point_id: scalar.initial_state(material) for p in points}
    trajectory = []
    q_history, damage_history, kappa_history, opening_history = [], [], [], []
    previous_q = 0.0
    cumulative_control_work = 0.0
    for step_index, fraction in enumerate(fractions):
        Q = fraction * q_star
        q = Q * unit_q
        finite(np, f"q step {step_index}", q)
        bulk_internal = K @ q
        gauge_internal = C.T @ (Q * unit_lambda_phys)
        bulk_residual = bulk_internal + gauge_internal
        assembly = trace.assemble_enriched_trial(material, points, tuple(float(v) for v in q), committed, rigid_candidates)
        gamma_force = np.asarray(assembly.generalized_internal_force_n, dtype=np.float64)
        finite(np, f"cohesive force step {step_index}", gamma_force)
        other = bulk_residual + gamma_force
        other_control = other.copy()
        other_control[CONTROL_INDEX] = 0.0
        force_scale = max(abs(float(other[CONTROL_INDEX])), float(np.linalg.norm(other, ord=np.inf)), 1.0)
        if float(np.linalg.norm(other_control, ord=np.inf)) > FORCE_ABS_TOL_N + FORCE_REL_TOL * force_scale:
            raise GateFailure("generalized_force_support", {"step": step_index,
                "noncontrol_generalized_force_inf_N": float(np.linalg.norm(other_control, ord=np.inf))})
        physical_control = float(other[CONTROL_INDEX])
        point_rows, primitive_cohesive_increment, analytic_cohesive_increment = [], 0.0, 0.0
        next_states = {}
        step_damage, step_kappa, step_opening = [], [], []
        for point, factor, item in zip(points, alpha, assembly.point_results):
            old = committed[point.point_id]
            state = item.trial_state
            jump = np.asarray(point.jump_operator, dtype=np.float64) @ q
            opening = float(np.dot(normal, jump))
            tangent = jump - opening * normal
            tangential_norm = float(np.linalg.norm(tangent))
            if opening < 0.0:
                raise GateFailure("opening_sign", {"step": step_index, "quadrature_id": point.quadrature_id})
            if (not close_enough(opening, factor * Q, TRACE_ABS_TOL_M, TRACE_REL_TOL)
                    or tangential_norm > TRACE_ABS_TOL_M + TRACE_REL_TOL * float(np.linalg.norm(jump))):
                raise GateFailure("mode_i_trace", {"step": step_index, "quadrature_id": point.quadrature_id,
                    "opening_m": opening, "expected_alpha_Q_m": factor * Q,
                    "tangential_jump_m": tangent.tolist(), "decision": "FAIL; no projection"})
            area = point.one_side_area_weight_m2
            primitive_work = area * (state.external_work_j_per_m2 - old.external_work_j_per_m2)
            law_work = area * law_interval_work(material, old, opening)
            primitive_cohesive_increment += primitive_work
            analytic_cohesive_increment += law_work
            state_row = {"opening_m": state.opening_m, "kappa_m": state.kappa_m, "damage": state.damage,
                "stored_energy_j_per_m2": state.stored_energy_j_per_m2,
                "dissipated_energy_j_per_m2": state.dissipated_energy_j_per_m2,
                "cumulative_external_work_j_per_m2": state.external_work_j_per_m2,
                "update_count": state.update_count}
            point_rows.append({"point_id": list(item.point_id), "point_xyz_m": list(point.positive_reference_point_m),
                "area_weight_m2": area, "analytic_alpha": factor, "opening_m": opening,
                "normal_opening_recomputed_m": float(np.dot(normal, jump)),
                "jump_vector_m": jump.tolist(), "tangential_jump_vector_m": tangent.tolist(),
                "tangential_norm_m": tangential_norm, "traction_pa": item.traction_pa,
                "k_alg_pa_per_m": item.k_alg_pa_per_m, "tangent_side": item.tangent_side,
                "tangent_status": getattr(item.tangent_status, "value", str(item.tangent_status)),
                "nondifferentiable_reasons": list(item.nondifferentiable_reasons), "state": state_row,
                "generalized_internal_force_N": list(item.generalized_force_n),
                "positive_face_internal_force_N": list(item.positive_face_force_at_reference_n),
                "negative_face_internal_force_N": list(item.negative_face_force_at_reference_n),
                "upper_physical_traction_N_per_m2": [-item.traction_pa * float(v) for v in normal],
                "lower_physical_traction_N_per_m2": [item.traction_pa * float(v) for v in normal],
                "cohesive_work_increment_primitive_j": primitive_work,
                "cohesive_work_increment_independent_law_j": law_work})
            next_states[point.point_id] = state
            step_damage.append(state.damage)
            step_kappa.append(state.kappa_m)
            step_opening.append(opening)
        work_scale = max(abs(primitive_cohesive_increment), abs(analytic_cohesive_increment),
                         material.gc_j_per_m2 * float(case["geometry"]["radius_m"]) * float(case["geometry"]["length_m"]), 1.0e-300)
        work_error = abs(primitive_cohesive_increment - analytic_cohesive_increment) / work_scale
        if work_error > ENERGY_REL_TOL:
            raise GateFailure("cohesive_work", {"step": step_index, "primitive_increment_j": primitive_cohesive_increment,
                "independent_law_integral_j": analytic_cohesive_increment, "relative_error": work_error})
        bulk_energy = 0.5 * float(q @ K @ q)
        cohesive_stored = float(assembly.stored_energy_j)
        cohesive_dissipated = float(assembly.dissipated_energy_j)
        bulk_work = 0.5 * float(unit_q @ K @ unit_q) * (Q*Q - previous_q*previous_q)
        total_work_increment = bulk_work + analytic_cohesive_increment
        old_stored = sum(p.one_side_area_weight_m2 * committed[p.point_id].stored_energy_j_per_m2 for p in points)
        old_diss = sum(p.one_side_area_weight_m2 * committed[p.point_id].dissipated_energy_j_per_m2 for p in points)
        energy_change = bulk_energy + cohesive_stored + cohesive_dissipated - (0.5 * float(unit_q @ K @ unit_q) * previous_q**2 + old_stored + old_diss)
        energy_error = abs(total_work_increment - energy_change) / max(abs(total_work_increment), abs(energy_change), 1.0e-300)
        if energy_error > ENERGY_REL_TOL:
            raise GateFailure("global_energy_balance", {"step": step_index, "control_work_increment_j": total_work_increment,
                "stored_plus_dissipated_change_j": energy_change, "relative_error": energy_error})
        cumulative_control_work += total_work_increment
        if any(state.damage + 1.0e-14 < committed[pid].damage or state.kappa_m + 1.0e-15 < committed[pid].kappa_m
               for pid, state in next_states.items()):
            raise GateFailure("irreversible_history", {"step": step_index})
        if step_index == 0:
            zero_corner = {"assembly_tangent_status": getattr(assembly.tangent_status, "value", str(assembly.tangent_status)),
                "point_statuses": sorted(set(row["tangent_status"] for row in point_rows)),
                "reasons": sorted(set(reason for row in point_rows for reason in row["nondifferentiable_reasons"])),
                "interpretation": "recorded primitive one-sided/corner status at Q=0; no tangent projection or solve uses it"}
        trajectory.append({"step_index": step_index, "Q_fraction_of_Qstar": fraction, "Q_m": Q,
            "q_physics_m": q.tolist(), "bulk_internal_generalized_force_N": bulk_internal.tolist(),
            "physical_gauge_generalized_force_N": gauge_internal.tolist(),
            "bulk_plus_gauge_residual_N": bulk_residual.tolist(), "cohesive_generalized_internal_force_N": gamma_force.tolist(),
            "control_applied_force_positive_opening_N": physical_control,
            "noncontrol_force_inf_N": float(np.linalg.norm(other_control, ord=np.inf)),
            "primitive_aggregate_cumulative_external_work_j": float(assembly.external_work_j),
            "cohesive_work_increment_primitive_j": primitive_cohesive_increment,
            "cohesive_work_increment_independent_law_j": analytic_cohesive_increment,
            "bulk_work_increment_exact_j": bulk_work, "control_work_increment_exact_j": total_work_increment,
            "control_work_cumulative_exact_j": cumulative_control_work, "bulk_stored_energy_j": bulk_energy,
            "cohesive_stored_energy_j": cohesive_stored, "cohesive_dissipated_energy_j": cohesive_dissipated,
            "energy_balance_relative_error": energy_error, "primitive_vs_independent_law_relative_error": work_error,
            "maximum_damage": max(step_damage), "minimum_damage": min(step_damage),
            "partial_damage_point_count": sum(0.0 < value < 1.0 for value in step_damage),
            "all_points_fully_damaged": all(value == 1.0 for value in step_damage),
            "points": point_rows})
        q_history.append(q.copy())
        damage_history.append(step_damage)
        kappa_history.append(step_kappa)
        opening_history.append(step_opening)
        committed = next_states
        previous_q = Q
        if time.monotonic() - elapsed_start > WALL_LIMIT_S:
            raise GateFailure("wall_budget", {"limit_seconds": WALL_LIMIT_S, "step_reached": step_index})

    if not trajectory[3]["partial_damage_point_count"] or not trajectory[5]["partial_damage_point_count"]:
        raise GateFailure("partial_damage_cycle", "partial damage was not represented at the initial and reloaded 0.30 Qstar states")
    if not trajectory[8]["all_points_fully_damaged"] or not trajectory[11]["all_points_fully_damaged"]:
        raise GateFailure("full_damage_cycle", "the declared peak did not damage every point")
    q_bytes = np.asarray(q_history, dtype="<f8").tobytes()
    damage_bytes = np.asarray(damage_history, dtype="<f8").tobytes()
    kappa_bytes = np.asarray(kappa_history, dtype="<f8").tobytes()
    opening_bytes = np.asarray(opening_history, dtype="<f8").tobytes()
    result = {"schema": SCHEMA, "status": "DIAGNOSTIC_ONLY", "classification": case["classification"],
        "execution_authorized_by_root_go": True, "source_hashes": hashes,
        "counts": {"new_mls_queries": 0, "provider_calls": 0, "npz_opened": False,
            "bulk_unit_control_kkt_order": 100, "bulk_unit_control_kkt_solves": 1,
            "primitive_assemblies": len(trajectory), "primitive_successful_returns": len(trajectory)},
        "saved_input_binding": {"saved_result_path": case["saved_result"]["path"],
            "saved_result_sha256": case["saved_result"]["sha256"], "coefficient_ids": list(coefficient_ids),
            "source_positions_xyz_m": source_positions.tolist(), "K_raw_N_per_m": K.tolist(),
            "C_phys": C.tolist(), "rigid_coefficient_modes": rigid.tolist(),
            "physical_external_force_N": [0.0] * NDof, "physical_gauge_target_m": [0.0] * 6,
            "control_semantic_id": coefficient_ids[CONTROL_INDEX], "control_index_zero_based": CONTROL_INDEX,
            "interfaces_reconstructed_from_saved_Hplus_Hminus": interfaces},
        "unit_control_response": {"q_per_Q": unit_q.tolist(), "physical_gauge_lambda_per_Q": unit_lambda_phys.tolist(),
            "saddle_control_multiplier_per_Q": float(unit_saddle_control_multiplier),
            "bulk_applied_control_force_per_Q_N_per_m": float(bulk_unit_residual[CONTROL_INDEX]),
            "KKT_backward_error": kkt_error, "physical_gauge_residual_inf": float(np.linalg.norm(C @ unit_q, ord=np.inf)),
            "unit_control_value": float(unit_q[CONTROL_INDEX]), "noncontrol_bulk_residual_inf": float(np.linalg.norm(noncontrol_unit, ord=np.inf)),
            "alpha_i": alpha, "tangential_trace_per_Q_verified_zero": True,
            "K_raw_N_per_m": K.tolist(), "C_phys": C.tolist()},
        "control_definition": "positive applied opening resistance = (Kq + C.T*lambda_phys + fGamma)[82]; body and prescribed boundary loads are zero",
        "Qstar_m": q_star, "trajectory_fraction_Qstar": list(fractions), "zero_opening_corner": zero_corner,
        "trajectory": trajectory,
        "trajectory_binary_le_f64": {"format": "IEEE-754 little-endian float64, base64",
            "q_physics_shape": [len(q_history), NDof], "q_physics": base64.b64encode(q_bytes).decode("ascii"),
            "damage_shape": [len(damage_history), len(points)], "damage": base64.b64encode(damage_bytes).decode("ascii"),
            "kappa_m": base64.b64encode(kappa_bytes).decode("ascii"),
            "opening_m": base64.b64encode(opening_bytes).decode("ascii")},
        "energy_interpretation": {"bulk_increment": "exact 0.5*k_bulk*(Q_new^2-Q_old^2)",
            "cohesive_increment": "independent exact piecewise-linear law integral, checked against primitive cumulative-work differences",
            "total_balance": "control path work checked against Delta(U_bulk + cohesive stored + cohesive dissipated)",
            "physical_control_force_is_recorded_at_each_step_but_no_trapezoid_is_used" : True},
        "limitations": ["fixed-face global Mode-I opening-control toy history",
            "no autonomous growth or crack advance", "no real material qualification", "no roughness or multibody claim",
            "not a fracture package"]}
    finite(np, "q history", np.asarray(q_history))
    finite(np, "damage history", np.asarray(damage_history))
    finite(np, "kappa history", np.asarray(kappa_history))
    if time.monotonic() - elapsed_start > WALL_LIMIT_S:
        raise GateFailure("wall_budget", {"limit_seconds": WALL_LIMIT_S, "phase": "serialization preparation"})
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="root-go-gated fixed-face opening-controlled cohesive history")
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--root-go")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    if not args.execute or args.output is None:
        print(json.dumps({"status": "EXECUTION_DISABLED", "numeric_modules_imported": False, "npz_opened": False}))
        return 2
    try:
        case = json.loads(CASE_PATH.read_text(encoding="utf-8"))
        pins = json.loads(PINS_PATH.read_text(encoding="utf-8"))
        require_execution(pins, args.root_go)
        allowed_output = (ROOT / Path(case["output_path"])).resolve()
        requested_output = args.output.resolve()
        if requested_output != allowed_output:
            raise GateFailure("output_path", "output must match the frozen review artifact path")
        result = execute(case, pins)
        code = 0
    except Exception as exc:
        result = {"schema": SCHEMA, "status": "NOT_ACCEPTED", "failure": {
            "gate": getattr(exc, "gate", "execution"), "type": type(exc).__name__, "detail": getattr(exc, "detail", str(exc))},
            "numeric_modules_imported": type(exc).__name__ != "PermissionError", "npz_opened": False}
        code = 3
    args.output.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if len(payload.encode("utf-8")) > OUTPUT_LIMIT_BYTES:
        result = {"schema": SCHEMA, "status": "NOT_ACCEPTED", "failure": {
            "gate": "output_budget", "detail": "result exceeded the frozen 8 MiB output cap"}}
        payload = json.dumps(result, allow_nan=False) + "\n"
        code = 3
    with args.output.open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(payload)
    print(json.dumps({"status": result["status"], "result": args.output.name,
        "primitive_assemblies": result.get("counts", {}).get("primitive_assemblies", 0),
        "failure": result.get("failure")}, allow_nan=False))
    return code


if __name__ == "__main__":
    raise SystemExit(main())
