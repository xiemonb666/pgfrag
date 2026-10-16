"""Bounded native, perfectly bonded bulk control for the frozen vessel mesh.

This is a low-load compliance diagnostic on the exact shared-node T4 mesh used
by weak_seam_v4. It intentionally creates no cohesive elements and carries no
physical-fracture or material-calibration qualification.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time
import traceback

import numpy as np


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path: Path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def inventory(root: Path):
    return {str(p.relative_to(root)).replace("\\", "/"):
            {"bytes": int(p.stat().st_size), "sha256": sha256(p)}
            for p in sorted(root.rglob("*")) if p.is_file() and p.name != "run_receipt.json"}


def load_factor(time_s: float, duration_s: float, peak_n: float) -> float:
    if not 0.0 < time_s < duration_s:
        return 0.0
    return float(peak_n * np.sin(np.pi * time_s / duration_s) ** 2)


def expected_lumped_mass(nodes: np.ndarray, tets: np.ndarray, density: float) -> np.ndarray:
    x = np.asarray(nodes, dtype=np.float64)[np.asarray(tets, dtype=np.int64)]
    jac = x[:, 1:, :] - x[:, :1, :]
    volumes = np.abs(np.linalg.det(jac)) / 6.0
    if not np.isfinite(volumes).all() or np.any(volumes <= 0.0):
        raise ValueError("frozen T4 mesh contains a non-positive or nonfinite cell volume")
    expected = np.zeros(int(len(nodes)), dtype=np.float64)
    np.add.at(expected, np.asarray(tets, dtype=np.int64).reshape(-1),
              np.repeat(float(density) * volumes / 4.0, 4))
    return expected


def compare_reference_load(generated: dict, reference_path: Path):
    reference = json.loads(reference_path.read_text(encoding="utf-8"))
    if len(reference.get("patches", [])) != 2 or len(generated.get("patches", [])) != 2:
        raise ValueError("expected two reference and two generated loading patches")
    patch_results = []
    for ref, got in zip(reference["patches"], generated["patches"]):
        if int(ref["sign"]) != int(got["sign"]):
            raise ValueError("cohesive low-load and bulk patch signs differ")
        ref_ids = np.asarray(ref["native_node_ids"], dtype=np.int64)
        got_ids = np.asarray(got["native_node_ids"], dtype=np.int64)
        ref_weights = np.asarray(ref["weights"], dtype=np.float64)
        got_weights = np.asarray(got["weights"], dtype=np.float64)
        same_ids = bool(np.array_equal(ref_ids, got_ids))
        max_abs = float(np.max(np.abs(ref_weights - got_weights))) if same_ids else None
        patch_results.append({"sign": int(got["sign"]), "same_native_patch_node_ids": same_ids,
                              "max_abs_weight_difference": max_abs,
                              "weight_count_reference": int(len(ref_weights)),
                              "weight_count_bulk": int(len(got_weights))})
    ref_com = np.asarray(reference["target_reference_com_m"], dtype=np.float64)
    got_com = np.asarray(generated["target_reference_com_m"], dtype=np.float64)
    com_error = float(np.max(np.abs(ref_com - got_com)))
    weights_match = all(item["same_native_patch_node_ids"] and
                        item["max_abs_weight_difference"] is not None and
                        item["max_abs_weight_difference"] <= 2.0e-12
                        for item in patch_results)
    return {"reference_artifact": str(reference_path.name),
            "reference_sha256": sha256(reference_path), "patches": patch_results,
            "max_abs_target_com_difference_m": com_error,
            "same_weights_within_2e-12": bool(weights_match),
            "target_com_equal_within_1e-12_m": bool(com_error <= 1.0e-12)}


def run(input_dir: Path, out: Path):
    import akantu as aka

    input_dir = input_dir.resolve()
    frozen = json.loads((input_dir / "input_freeze.json").read_text(encoding="utf-8"))
    for name, expected in frozen["files"].items():
        candidate = input_dir / name
        if not candidate.is_file() or sha256(candidate) != expected["sha256"]:
            raise ValueError(f"frozen input hash mismatch: {name}")

    # These helpers are imported only from the byte-pinned weak_seam_v4 source
    # snapshot; the input bundle is selected before importing it.
    sys.path.insert(0, str(input_dir))
    from cohesive_mesh import read_gmsh22
    from cohesive_vessel import make_balanced_patch_load

    control = json.loads((input_dir / "control.json").read_text(encoding="utf-8"))
    source_case = json.loads((input_dir / "weak_seam_case.json").read_text(encoding="utf-8"))
    cohesive_reference_case = json.loads((input_dir / "cohesive_lowload_case_reference.json").read_text(encoding="utf-8"))
    if source_case["case_name"] != "weak_seam_v4":
        raise ValueError("source case is not the frozen weak_seam_v4")
    for name, key in (("mesh.msh", "mesh.msh"), ("mesh_result.json", "mesh_result.json"),
                      ("forge_manifest.json", "forge_manifest.json"),
                      ("weak_seam_material_source.dat", "material.dat"),
                      ("cohesive_vessel.py", "cohesive_vessel.py"),
                      ("cohesive_mesh.py", "cohesive_mesh.py"),
                      ("cohesive_stability.py", "cohesive_stability.py")):
        if source_case["input_hashes"][key] != sha256(input_dir / name):
            raise ValueError(f"weak_seam_v4 embedded source pin differs: {name}")
    if control["material"] != source_case["bulk_material"]:
        raise ValueError("control bulk material differs from weak_seam_v4")
    if cohesive_reference_case["case_name"] != "bonded_cohesive_control_v1":
        raise ValueError("loading reference is not the frozen native bonded cohesive low-load control")
    if (cohesive_reference_case["source_forge_manifest_sha256"] != source_case["source_forge_manifest_sha256"] or
        cohesive_reference_case["source_mesh_result_sha256"] != source_case["source_mesh_result_sha256"] or
        cohesive_reference_case["bulk_material"]["rho_kg_m3"] != control["material"]["rho_kg_m3"] or
        cohesive_reference_case["bulk_material"]["E_pa"] != control["material"]["E_pa"] or
        cohesive_reference_case["bulk_material"]["nu"] != control["material"]["nu"] or
        cohesive_reference_case["loading"]["peak_force_N"] != control["loading"]["peak_force_N"] or
        cohesive_reference_case["loading"]["pulse_duration_s"] != control["loading"]["pulse_duration_s"] or
        cohesive_reference_case["loading"]["patch_x_threshold_m"] != control["loading"]["patch_x_threshold_m"] or
        cohesive_reference_case["time"]["total_s"] != control["time"]["total_s"]):
        raise ValueError("frozen cohesive low-load reference does not share this bulk/mesh/load contract")
    reference_result = json.loads((input_dir / "cohesive_lowload_result_reference.json").read_text(encoding="utf-8"))
    reference_bound = json.loads((input_dir / "cohesive_lowload_initial_bound_reference.json").read_text(encoding="utf-8"))
    if (reference_result.get("status") != "NATIVE_RUN_COMPLETE_DIAGNOSTIC" or
        int(reference_result.get("solver_steps", -1)) != int(control["time"]["steps"]) or
        int(reference_bound.get("steps", -1)) != int(control["time"]["steps"]) or
        not bool(reference_bound.get("all_dofs_free")) or
        abs(float(reference_bound.get("dt_s", 0.0)) - float(control["time"]["total_s"]) / int(control["time"]["steps"])) > 1.0e-18):
        raise ValueError("cohesive low-load reference run did not complete the same free-DOF time contract")
    if control["time"]["steps"] != 3244 or control["time"]["total_s"] != 50.0e-6:
        raise ValueError("reviewed fixed T/steps contract changed")
    if control["loading"]["peak_force_N"] != 100.0 or control["loading"]["pulse_duration_s"] != 25.0e-6:
        raise ValueError("reviewed 100 N / 25 us low-load contract changed")

    source = read_gmsh22(input_dir / "mesh.msh")
    aka.parseInput(str(input_dir / "material.dat"))
    mesh = aka.Mesh(3)
    mesh.read(str(input_dir / "mesh.msh"))
    native_nodes_before = np.asarray(mesh.getNodes()).copy()
    native_tets_before = np.asarray(mesh.getConnectivity(aka._tetrahedron_4)).copy()
    source_tets = np.asarray([row["nodes"] for row in source["tetrahedra"]], dtype=np.int64)
    if not np.array_equal(native_nodes_before, source["nodes"]):
        raise ValueError("native node ordering/coordinates differ from source Gmsh snapshot")
    if not np.array_equal(native_tets_before, source_tets):
        raise ValueError("native T4 ordering/connectivity differs from source Gmsh snapshot")
    model = aka.SolidMechanicsModel(mesh)
    model.initFull(_analysis_method=aka._explicit_lumped_mass)
    model.assembleMassLumped()
    rest = np.asarray(mesh.getNodes()).copy()
    tets = np.asarray(mesh.getConnectivity(aka._tetrahedron_4)).copy()
    if not np.array_equal(rest, native_nodes_before) or not np.array_equal(tets, native_tets_before):
        raise ValueError("perfectly bonded model altered initial node/tet source ordering")
    mass3 = np.asarray(model.getMass()).copy()
    if mass3.shape != (int(len(rest)), 3) or not np.isfinite(mass3).all():
        raise ValueError(f"unexpected native mass array shape/content: {mass3.shape}")
    blocked_view = np.asarray(model.getBlockedDOFs())
    blocked_view[:] = False
    blocked = np.asarray(model.getBlockedDOFs()).copy()
    if blocked.shape != rest.shape or blocked.any():
        raise ValueError("the control requires every native DOF to remain free")

    density = float(control["material"]["rho_kg_m3"])
    independent_mass = expected_lumped_mass(rest, tets, density)
    native_mass = mass3[:, 0]
    mass_diff = mass3 - independent_mass[:, None]
    max_mass_abs = float(np.max(np.abs(mass_diff)))
    rel_denominator = np.maximum(independent_mass[:, None], max(float(independent_mass.max()) * 1.0e-15, 1.0e-30))
    max_mass_rel = float(np.max(np.abs(mass_diff) / rel_denominator))
    mass_pass = bool(max_mass_abs <= 1.0e-13 and max_mass_rel <= 1.0e-9)

    ancestry = np.arange(int(len(rest)), dtype=np.int64)
    unit_load, load_report = make_balanced_patch_load(
        source, rest, ancestry, native_mass, float(control["loading"]["patch_x_threshold_m"]))
    load_compare = compare_reference_load(load_report, input_dir / "cohesive_lowload_loading_reference.json")

    total_time = float(control["time"]["total_s"])
    steps = int(control["time"]["steps"])
    if steps > int(control["limits"]["max_steps"]):
        raise ValueError("fixed control step count exceeds its declared 4000-step cap")
    dt = total_time / int(steps)
    bulk_dt = float(model.getStableTimeStep())
    dt_limit = float(control["time"]["stable_dt_fraction"]) * bulk_dt
    dt_pass = bool(dt <= dt_limit * (1.0 + 1.0e-12))
    threshold = float(control["loading"]["patch_x_threshold_m"])
    load_ids = np.flatnonzero(unit_load[:, 0])
    negative_patch = next(p for p in load_report["patches"] if int(p["sign"]) == -1)
    positive_patch = next(p for p in load_report["patches"] if int(p["sign"]) == 1)
    neg_ids = np.asarray(negative_patch["native_node_ids"], dtype=np.int64)
    pos_ids = np.asarray(positive_patch["native_node_ids"], dtype=np.int64)
    neg_weights = np.asarray(negative_patch["weights"], dtype=np.float64)
    pos_weights = np.asarray(positive_patch["weights"], dtype=np.float64)

    saved_steps = np.asarray(sorted(set([int(k) for k in range(0, int(steps) + 1, 50)] + [int(steps)])), dtype=np.int64)
    forecast_bytes = (int(steps + 1) * 26 * 8 + int(len(saved_steps)) * int(len(rest)) * 4 * 3 * 8 +
                      sum(int(p.stat().st_size) for p in input_dir.iterdir() if p.is_file()) + 1_000_000)
    if forecast_bytes > int(control["output"]["output_cap_bytes"]):
        raise ValueError(f"bounded output forecast {forecast_bytes} exceeds the 16 MiB cap")
    node_states = np.zeros((int(len(saved_steps)), int(len(rest)), 4, 3), dtype=np.float64)
    state_step_to_row = {int(step): int(row) for row, step in enumerate(saved_steps)}
    columns = ["time_s", "load_factor_N", "patch_left_mean_Ux_m", "patch_right_mean_Ux_m",
               "right_minus_left_response_m", "kinetic_energy_J", "potential_energy_J",
               "external_work_J", "external_work_minus_kinetic_and_potential_J",
               "momentum_x_kg_mps", "momentum_y_kg_mps", "momentum_z_kg_mps",
               "angular_momentum_x_kg_m2ps", "angular_momentum_y_kg_m2ps", "angular_momentum_z_kg_m2ps",
               "external_torque_x_Nm", "external_torque_y_Nm", "external_torque_z_Nm",
               "internal_torque_x_Nm", "internal_torque_y_Nm", "internal_torque_z_Nm",
               "max_abs_force_balance_residual_N", "max_abs_displacement_m", "max_abs_velocity_mps",
               "max_abs_acceleration_mps2", "max_abs_internal_force_N"]
    history = np.zeros((int(steps) + 1, int(len(columns))), dtype=np.float64)
    states_path = out / "initial_state.npz"
    np.savez_compressed(states_path, rest_nodes_m=rest.astype(np.float64, copy=False),
                        tets=tets.astype(np.int64, copy=False),
                        source_node_tags=np.asarray(source["node_tags"], dtype=np.int64),
                        source_tet_ids=np.asarray([e["element_id"] for e in source["tetrahedra"]], dtype=np.int64),
                        native_to_source=ancestry, mass_kg=mass3, blocked=blocked,
                        unit_external_load_N_per_N=unit_load)
    write_json(out / "lineage.json", {
        "native_initial_node_count": int(len(rest)), "native_initial_t4_count": int(len(tets)),
        "cohesive_element_count": 0, "native_initial_nodes_exact_source_order": True,
        "native_initial_t4_exact_source_order": True, "native_to_source_zero_based": "arange(node_count)",
        "source_zero_ancestry_exact": bool(np.array_equal(ancestry, np.arange(int(len(rest)), dtype=np.int64))),
        "all_native_dofs_free": bool(not blocked.any()),
        "mesh_sha256": sha256(input_dir / "mesh.msh")})
    write_json(out / "mass_audit.json", {
        "independent_method": "rho * abs(det(T4 edge matrix)) / 6, equally lumped to all four nodes",
        "nodes_checked": int(len(rest)), "native_components_checked": int(mass3.size),
        "native_total_mass_kg_component0": float(native_mass.sum()),
        "independent_total_mass_kg": float(independent_mass.sum()),
        "native_mass_com_xyz_m": (native_mass[:, None] * rest).sum(axis=0).tolist(),
        "max_abs_per_node_component_mass_difference_kg": max_mass_abs,
        "max_relative_per_node_component_mass_difference": max_mass_rel,
        "all_node_mass_pass": mass_pass})
    write_json(out / "loading.json", load_report)
    write_json(out / "load_weight_comparison.json", load_compare)
    bound_report = {
        "steps": int(steps), "dt_s": float(dt), "native_bulk_stable_dt_s": float(bulk_dt),
        "dt_fraction_limit": float(control["time"]["stable_dt_fraction"]),
        "dt_fraction_limit_s": float(dt_limit), "dt_within_native_fraction_limit": dt_pass,
        "mass_all_node_audit_pass": mass_pass,
        "same_lowload_cohesive_patch_node_ids_and_weights": load_compare["same_weights_within_2e-12"],
        "all_native_dofs_free": bool(not blocked.any()),
        "mesh_counts": {"nodes": int(len(rest)), "tets": int(len(tets)), "cohesive_elements": 0},
        "mass_kg": float(native_mass.sum()), "output_bytes_forecast_including_inputs_and_margin": int(forecast_bytes),
        "scope": "native fully bonded bulk stability check only"}
    write_json(out / "initial_bound.json", bound_report)

    initial_u = np.asarray(model.getDisplacement()).copy()
    initial_v = np.asarray(model.getVelocity()).copy()
    initial_a = np.asarray(model.getAcceleration()).copy()
    initial_internal = np.asarray(model.getInternalForce()).copy()
    if any(not np.isfinite(a).all() for a in (initial_u, initial_v, initial_a, initial_internal)):
        raise ValueError("nonfinite zero-load initial native state")
    initial_energies = np.asarray([float(model.getEnergy("kinetic")), float(model.getEnergy("potential"))], dtype=np.float64)
    if not np.isfinite(initial_energies).all() or np.max(np.abs(initial_energies)) > 1.0e-15:
        raise ValueError(f"the zero-load source state must have zero native energy: {initial_energies.tolist()}")
    if not mass_pass or not load_compare["same_weights_within_2e-12"] or not load_compare["target_com_equal_within_1e-12_m"]:
        raise ValueError("mass or loading-reference audit failed before time integration")
    if not dt_pass:
        raise ValueError("the fixed baseline dt does not pass the native 0.4*bulk-stable-dt limit")

    model.setTimeStep(float(dt))
    previous_u = initial_u.copy()
    previous_v = initial_v.copy()
    previous_a = initial_a.copy()
    previous_force = np.zeros_like(rest)
    previous_internal = initial_internal.copy()
    external_work = 0.0
    saved_index = 0
    peak_index = 0
    started = time.monotonic()
    completed = 0
    diagnostics_path = out / "history.npz"
    try:
        for step in range(int(steps) + 1):
            time_s = float(step) * float(dt)
            factor = load_factor(time_s, float(control["loading"]["pulse_duration_s"]),
                                 float(control["loading"]["peak_force_N"]))
            applied = factor * unit_load
            if step:
                np.asarray(model.getExternalForce())[:] = applied
                model.solveStep()
                completed = int(step)
            u = np.asarray(model.getDisplacement()).copy()
            v = np.asarray(model.getVelocity()).copy()
            a = np.asarray(model.getAcceleration()).copy()
            internal = np.asarray(model.getInternalForce()).copy()
            external = np.asarray(model.getExternalForce()).copy()
            energies = np.asarray([float(model.getEnergy("kinetic")), float(model.getEnergy("potential"))], dtype=np.float64)
            if not all(np.isfinite(q).all() for q in (u, v, a, internal, external, energies)):
                raise ValueError(f"nonfinite native state at step {step}")
            if blocked.any() or not np.array_equal(external, applied):
                raise ValueError(f"boundary constraints or force alteration detected at step {step}")
            if not np.array_equal(np.asarray(mesh.getNodes()), rest) or not np.array_equal(
                    np.asarray(mesh.getConnectivity(aka._tetrahedron_4)), tets):
                raise ValueError(f"native mesh topology/rest coordinates changed at step {step}")
            neg_mean = float(neg_weights @ u[neg_ids, 0])
            pos_mean = float(pos_weights @ u[pos_ids, 0])
            response = pos_mean - neg_mean
            if step:
                external_work += float(np.sum(0.5 * (previous_force + external) * (u - previous_u)))
            momentum = np.sum(mass3[:, :1] * v, axis=0)
            current = rest + u
            angular = np.cross(current, mass3[:, :1] * v).sum(axis=0)
            external_torque = np.cross(current, external).sum(axis=0)
            internal_torque = np.cross(current, internal).sum(axis=0)
            force_residual = mass3 * a - internal - external
            row = [time_s, factor, neg_mean, pos_mean, response, energies[0], energies[1], external_work,
                   external_work - float(energies.sum()), *momentum.tolist(), *angular.tolist(),
                   *external_torque.tolist(), *internal_torque.tolist(),
                   float(np.abs(force_residual).max()), float(np.abs(u).max()), float(np.abs(v).max()),
                   float(np.abs(a).max()), float(np.abs(internal).max())]
            history[int(step)] = np.asarray(row, dtype=np.float64)
            if step in state_step_to_row:
                node_states[state_step_to_row[int(step)]] = np.stack((u, v, a, internal), axis=1)
                saved_index += 1
            if factor > history[peak_index, 1]:
                peak_index = int(step)
            if step % 250 == 0 or int(step) == int(steps):
                np.savez_compressed(diagnostics_path, diagnostics=history[:int(step) + 1],
                                    diagnostic_columns=np.asarray(columns), node_states=node_states[:int(saved_index)],
                                    node_saved_steps=saved_steps[:int(saved_index)])
                write_json(out / "checkpoint.json", {
                    "completed_steps": int(step), "node_state_frames_committed": int(saved_index),
                    "elapsed_s": float(time.monotonic() - started), "dt_s": float(dt)})
                with (out / "progress.jsonl").open("a", encoding="utf-8") as log:
                    log.write(json.dumps({"step": int(step), "of": int(steps),
                                          "elapsed_s": float(time.monotonic() - started)}) + "\n")
            if time.monotonic() - started > float(control["limits"]["wall_time_s"]):
                raise TimeoutError("declared 120-second wall-time cap reached")
            previous_u, previous_v, previous_a = u, v, a
            previous_force, previous_internal = external, internal
    finally:
        np.savez_compressed(diagnostics_path, diagnostics=history[:int(completed) + 1],
                            diagnostic_columns=np.asarray(columns), node_states=node_states[:int(saved_index)],
                            node_saved_steps=saved_steps[:int(saved_index)])

    max_factor_step = int(peak_index)
    peak_factor = float(history[max_factor_step, 1])
    peak_response = float(history[max_factor_step, 4])
    compliance = peak_response / peak_factor if peak_factor > 0.0 else math.nan
    max_response_index = int(np.argmax(np.abs(history[:completed + 1, 4])))
    elapsed = float(time.monotonic() - started)
    diagnostic = {
        "status": "NATIVE_BONDED_BULK_CONTROL_COMPLETE_DIAGNOSTIC",
        "solver_steps_completed": int(completed), "mesh_counts": [int(len(rest)), int(len(tets)), 0],
        "dt_s": float(dt), "native_bulk_stable_dt_s": float(bulk_dt), "elapsed_s": elapsed,
        "same_initial_node_and_t4_source_order": True, "native_to_source_zero_based_arange": True,
        "all_node_mass_pass": mass_pass, "load_weight_reference_match": load_compare,
        "all_dofs_free": bool(not blocked.any()),
        "compliance_diagnostic": {
            "definition": "right(+x) minus left(-x) force-weighted patch mean Ux / actual applied resultant force",
            "classification": "transient secant compliance of the fully bonded finite bulk under this pulse; not interface compliance",
            "peak_load_step": int(max_factor_step), "peak_load_time_s": float(history[max_factor_step, 0]),
            "actual_peak_applied_force_N": peak_factor, "response_at_peak_load_m": peak_response,
            "response_divided_by_applied_force_m_per_N": float(compliance),
            "inverse_secant_stiffness_N_per_m": float(1.0 / compliance) if compliance != 0.0 else None,
            "max_abs_response_step": int(max_response_index),
            "max_abs_response_time_s": float(history[max_response_index, 0]),
            "max_abs_response_m": float(abs(history[max_response_index, 4]))},
        "max_force_balance_residual_N": float(np.max(history[:completed + 1, 21])),
        "max_abs_external_torque_Nm": np.max(np.abs(history[:completed + 1, 15:18]), axis=0).tolist(),
        "max_abs_energy_balance_residual_J": float(np.max(np.abs(history[:completed + 1, 8]))),
        "qualification": "ACTUAL_NATIVE_BONDED_BULK_COMPLIANCE_DIAGNOSTIC_ONLY; NO_COHESIVE_INTERFACE; NO_PHYSICAL_OR_MATERIAL_CALIBRATION_QUALIFICATION"}
    write_json(out / "diagnostic.json", diagnostic)
    return diagnostic


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve()
    if not out.is_dir():
        raise FileNotFoundError(f"run directory must be pre-created and non-overwriting: {out}")
    if any((out / name).exists() for name in ("started.json", "result.json", "run_receipt.json")):
        raise FileExistsError("run outputs already exist; refusing overwrite")
    start = time.time()
    source_inventory = inventory(out)
    write_json(out / "started.json", {"argv": sys.argv, "started_unix_s": start,
                "python": sys.version, "executable": sys.executable,
                "input_directory": str(args.input_dir.resolve()),
                "preexisting_run_tree_inventory": source_inventory,
                "output_cap_bytes": 16 * 1024 * 1024})
    try:
        result = run(args.input_dir, out)
    except Exception as exc:
        result = {"status": "FAILED_PRESERVED_NO_RETRY", "error": repr(exc),
                  "traceback": traceback.format_exc(),
                  "qualification": "NO_PHYSICAL_OR_MATERIAL_CALIBRATION_QUALIFICATION"}
    write_json(out / "result.json", result)
    files = inventory(out)
    total_bytes = int(sum(entry["bytes"] for entry in files.values()))
    receipt = {"argv": sys.argv, "started_unix_s": start, "ended_unix_s": time.time(),
               "status": result["status"], "persistent_bytes_before_receipt": total_bytes,
               "output_cap_bytes": 16 * 1024 * 1024, "output_cap_pass": bool(total_bytes <= 16 * 1024 * 1024),
               "files": files}
    write_json(out / "run_receipt.json", receipt)
    print(json.dumps({"result": result, "receipt_bytes": total_bytes,
                      "output_cap_pass": receipt["output_cap_pass"]}, allow_nan=False), flush=True)
    return 0 if result["status"] == "NATIVE_BONDED_BULK_CONTROL_COMPLETE_DIAGNOSTIC" and receipt["output_cap_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
