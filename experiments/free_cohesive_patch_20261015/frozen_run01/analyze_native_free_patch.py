from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any

import numpy as np

from free_patch_contract import (
    build_native_lineage,
    read_gmsh22,
    sha256,
    strict_json_dump,
)


EXP = Path(__file__).resolve().parent
CASE_PATH = EXP / "case.json"
MESH_PATH = EXP / "block_2x2x4mm.msh"
MATERIAL_PATH = EXP / "material.dat"
ENERGY_KEYS = ("kinetic", "potential", "dissipated", "reversible", "cohesive contact")
NATIVE_QP_FIELDS = ("opening", "tractions", "damage", "delta_max")


def finite_all(arrays: dict[str, np.ndarray]) -> bool:
    return all(np.isfinite(value).all() for value in arrays.values())


def trace_rows(path: Path) -> list[dict[str, Any]]:
    rows = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        try:
            rows.append(json.loads(line))
        except json.JSONDecodeError as error:
            raise ValueError(f"invalid trace JSON at line {line_number}: {error}") from error
    return rows


def load_history(directory: Path, frame_count: int) -> dict[str, np.ndarray]:
    names = (
        "opening", "tractions", "damage", "delta_max", "traction_dot_opening_Pa_m_raw",
        "interface_work_increment_J_by_qp_raw", "energies_J_raw", "displacement_m",
        "internal_force_N_raw", "acceleration_m_s2_raw", "external_force_N_raw",
        "mass_kg_raw", "boundary_readback_error_m", "native_work_cumulative_J_raw",
        "native_energy_total_J_raw", "native_work_vs_total_residual_J_raw",
        "actuator_reaction_estimate_N_raw", "top_internal_force_z_sum_N_raw",
        "top_inertial_force_z_sum_N_raw", "top_external_force_z_sum_N_raw",
        "actuator_work_cumulative_J_raw", "energy_balance_residual_J_raw",
        "commanded_top_displacement_m", "topology_counts",
    )
    arrays = {}
    for name in names:
        path = directory / f"{name}.npy"
        if not path.is_file():
            raise FileNotFoundError(f"missing raw history array: {path.name}")
        arr = np.load(path, mmap_mode="r", allow_pickle=False)
        if arr.ndim == 0 or len(arr) < frame_count:
            raise ValueError(f"raw history {name} has {len(arr)} frames, expected at least {frame_count}")
        arrays[name] = arr[:frame_count]
    return arrays


def write_summary(path: Path, summary: dict[str, Any]) -> None:
    strict_json_dump(path, summary)


def analyze(run_dir: Path) -> dict[str, Any]:
    run_dir = run_dir.resolve()
    case = json.loads(CASE_PATH.read_text(encoding="utf-8"))
    status_path = run_dir / "run_status.json"
    started_path = run_dir / "run_started.json"
    if not status_path.is_file() or not started_path.is_file():
        raise FileNotFoundError("run_started.json and run_status.json are required")
    status = json.loads(status_path.read_text(encoding="utf-8"))
    started = json.loads(started_path.read_text(encoding="utf-8"))
    frame_count = int(status.get("history_frames_saved", 0))
    if frame_count < 1:
        raise ValueError("run did not save an initialized native history frame")
    arrays = load_history(run_dir / "raw_history", frame_count)
    trace = trace_rows(run_dir / "trace.jsonl")
    if len(trace) != frame_count:
        raise ValueError(f"trace rows {len(trace)} != raw history frames {frame_count}")

    with np.load(run_dir / "initial_state.npz", allow_pickle=False) as archive:
        initial = {key: archive[key].copy() for key in archive.files}
    mesh_data = read_gmsh22(MESH_PATH)
    rebuilt_lineage = build_native_lineage(
        mesh_data,
        initial["native_rest_nodes"],
        initial["native_bulk_tets"],
        initial["native_cohesive_connectivity"],
        expected_area_m2=float(case["geometry"]["central_interface"]["rest_area_m2"]),
    )
    lineage_path = run_dir / "native_lineage.json"
    saved_lineage = json.loads(lineage_path.read_text(encoding="utf-8"))
    lineage_match = (
        saved_lineage.get("input_face_count") == rebuilt_lineage["input_face_count"]
        and saved_lineage.get("native_cohesive_element_count") == rebuilt_lineage["native_cohesive_element_count"]
        and saved_lineage.get("native_node_to_input_node_zero_based") == rebuilt_lineage["native_node_to_input_node_zero_based"]
        and saved_lineage.get("parent_facets") == rebuilt_lineage["parent_facets"]
    )

    opening = arrays["opening"]
    traction = arrays["tractions"]
    damage = arrays["damage"][:, :, 0]
    delta_max = arrays["delta_max"][:, :, 0]
    energies = arrays["energies_J_raw"]
    native_work_increment = arrays["interface_work_increment_J_by_qp_raw"]
    independent_work_cumulative = np.cumsum(native_work_increment.sum(axis=1))
    reported_work_cumulative = arrays["native_work_cumulative_J_raw"]
    native_total_energy = arrays["native_energy_total_J_raw"]
    closure_residual = independent_work_cumulative - native_total_energy
    energy_work_tolerance = float(case["energy_contract"].get("interface_work_identity_tolerance_J", 1e-12))
    negative_tolerance = float(case["energy_contract"].get("negative_native_dissipation_absolute_tolerance_J", 1e-12))
    closure_error_max = float(np.max(np.abs(closure_residual), initial=0.0))
    cumulative_match_error = float(np.max(np.abs(independent_work_cumulative - reported_work_cumulative), initial=0.0))

    expected_frames = int(case["loading"]["total_steps"]) + 1
    expected_steps = int(case["loading"]["total_steps"])
    expected_counts = np.asarray((54, 96, 8), dtype=np.int64)
    counts = arrays["topology_counts"]
    topology_stable = bool(np.array_equal(counts, np.repeat(expected_counts[None, :], frame_count, axis=0)))
    boundary_errors = arrays["boundary_readback_error_m"]
    boundary_ok = bool(np.isfinite(boundary_errors).all() and np.max(boundary_errors, initial=0.0) <= 1e-12)
    field_shape_ok = bool(
        opening.shape == (frame_count, 24, 3)
        and traction.shape == (frame_count, 24, 3)
        and damage.shape == (frame_count, 24)
        and delta_max.shape == (frame_count, 24)
    )
    native_finite = finite_all(arrays)
    damage_ok = bool(
        np.isfinite(damage).all()
        and (damage >= -1e-12).all()
        and (damage <= 1.0 + 1e-12).all()
        and np.all(np.diff(damage, axis=0) >= -1e-12)
    )
    delta_ok = bool(np.isfinite(delta_max).all() and np.all(np.diff(delta_max, axis=0) >= -1e-12))
    negative_dissipation_count = int(np.count_nonzero(energies[:, 2] < -negative_tolerance))
    negative_dissipation_raw_count = int(np.count_nonzero(energies[:, 2] < 0.0))
    negative_dissipation_ok = negative_dissipation_count == 0
    work_identity_ok = bool(closure_error_max <= energy_work_tolerance and cumulative_match_error <= energy_work_tolerance)
    history_complete = bool(
        status.get("status") == "RUN_COMPLETE_RAW_HISTORY_CLOSED"
        and int(status.get("solver_steps_completed", -1)) == expected_steps
        and frame_count == expected_frames
        and [int(row["history_index"]) for row in trace] == list(range(expected_frames))
        and [int(row["solver_step"]) for row in trace] == list(range(expected_frames))
    )

    expected_blocked = int(case["boundary_conditions"]["expected_blocked_components_after_insertion"])
    initial_blocked = initial["native_blocked_dofs"]
    interface_nodes = initial["interface_native_node_ids"]
    interior_nodes = initial["interior_native_node_ids"]
    blocked_component_count = int(np.count_nonzero(initial_blocked))
    free_interface_components = int(interface_nodes.size * 3 - np.count_nonzero(initial_blocked[interface_nodes]))
    free_interior_components = int(interior_nodes.size * 3 - np.count_nonzero(initial_blocked[interior_nodes]))
    dofs_ok = bool(
        blocked_component_count == expected_blocked
        and initial_blocked.size - blocked_component_count == int(case["boundary_conditions"]["expected_free_dof_components"])
        and free_interface_components == int(case["boundary_conditions"]["expected_free_cohesive_dof_components"])
        and free_interior_components == int(case["boundary_conditions"]["expected_free_geometric_interior_bulk_dof_components"])
        and not np.any(initial_blocked[interface_nodes])
    )

    commanded = arrays["commanded_top_displacement_m"]
    expected_commanded = np.asarray([
        0.0 if i == 0 else min(1.0, i / int(case["loading"]["ramp_steps"]))
        for i in range(frame_count)
    ])
    expected_commanded = 0.5 * float(case["loading"]["top_target_displacement_m"]) * (
        1.0 - np.cos(np.pi * expected_commanded)
    )
    expected_commanded[np.arange(frame_count) >= int(case["loading"]["ramp_steps"])] = float(case["loading"]["top_target_displacement_m"])
    loading_ok = bool(np.max(np.abs(commanded - expected_commanded), initial=0.0) <= 1e-15)

    stability = status.get("stability_contract", {})
    solver_dt = float(case["loading"]["time_step_s"])
    observed_bulk_dt = stability.get("bulk_get_stable_time_step_s")
    bulk_fraction = float(case["loading"]["bulk_stable_step_fraction_limit"])
    bulk_stability_ok = bool(
        observed_bulk_dt is not None
        and np.isfinite(float(observed_bulk_dt))
        and float(observed_bulk_dt) > 0.0
        and solver_dt <= bulk_fraction * float(observed_bulk_dt)
        and stability.get("bulk_check_is_total_cohesive_cfl") is False
    )
    lambda_upper = case["stability"].get("lambda_max_upper_bound_per_s2")
    total_stability_status = case["stability"].get("total_stability_status")
    total_stable_dt = None
    if lambda_upper is not None and np.isfinite(float(lambda_upper)) and float(lambda_upper) > 0.0:
        total_stable_dt = 2.0 / float(np.sqrt(float(lambda_upper)))
    total_fraction = float(case["stability"]["fraction_limit"])
    observed_lambda = stability.get("lambda_max_upper_bound_per_s2")
    observed_total_dt = stability.get("derived_total_stable_time_step_s")
    total_stability_ok = bool(
        total_stability_status == "ESTABLISHED_PARENT_REVIEWED_BOUND"
        and total_stable_dt is not None
        and observed_lambda is not None
        and np.isclose(float(observed_lambda), float(lambda_upper), rtol=0.0, atol=0.0)
        and observed_total_dt is not None
        and np.isclose(float(observed_total_dt), total_stable_dt, rtol=1e-12, atol=0.0)
        and stability.get("total_stability_status") == total_stability_status
        and stability.get("total_gate_passed_before_step_1") is True
        and solver_dt <= total_fraction * total_stable_dt
    )

    # Distinguish a traction sign error from a native history/API mismatch.
    # The source formula uses +t·delta(opening), and the independent integral uses
    # the same raw signed work. A negative D with a closed identity is therefore
    # in the native path/history, not something to fix by flipping traction.
    dissipated = energies[:, 2]
    if negative_dissipation_count:
        if work_identity_ok:
            dissipation_diagnosis = (
                "NEGATIVE_NATIVE_D_MATCHES_SIGNED_HISTORY_IDENTITY: traction orientation is not the cause; "
                "the signed native work minus current reversible energy is negative. The state/history mechanism "
                "still needs a constitutive or timestep investigation, so this is a diagnostic FAIL."
            )
        else:
            dissipation_diagnosis = (
                "NEGATIVE_NATIVE_D_WITH_WORK_HISTORY_MISMATCH: traction sign remains unaltered; the source-formula "
                "identity does not close, leaving material update timing, wheel/source correspondence, or quadrature "
                "mapping unresolved. This is a diagnostic FAIL."
            )
    elif work_identity_ok:
        dissipation_diagnosis = (
            "NO_MATERIAL_NEGATIVE_D_IN_THIS_INTRINSIC_RUN_AND_SIGNED_HISTORY_IDENTITY_CLOSED; the former "
            "extrinsic post-insertion negative-D episode remains unexplained and is not claimed resolved."
        )
    else:
        dissipation_diagnosis = (
            "NO_LARGE_NEGATIVE_D_BUT_WORK_HISTORY_IDENTITY_FAILED; no sign correction is applied and the "
            "API/history discrepancy blocks a sign-history closure PASS."
        )

    full_debond = None
    hold_steps = int(case["loading"]["hold_steps"])
    sigma = float(case["material"]["cohesive"]["sigma_c_Pa"])
    full_damage_frames = np.flatnonzero(np.all(np.abs(damage - 1.0) <= 1e-12, axis=1))
    if full_damage_frames.size:
        for index in full_damage_frames:
            hold_end = min(frame_count, int(index) + hold_steps + 1)
            if hold_end <= index or hold_end - index < 2:
                continue
            damage_hold = np.all(np.abs(damage[index:hold_end] - 1.0) <= 1e-12)
            traction_hold = np.linalg.norm(traction[index:hold_end], axis=2).max(initial=0.0) <= 0.01 * sigma
            if damage_hold and traction_hold:
                full_debond = {"first_history_index": int(index), "hold_end_exclusive": int(hold_end)}
                break

    source_integrity = {
        "case_sha256_matches_started": started.get("case_sha256") == sha256(CASE_PATH),
        "mesh_sha256_matches_case": sha256(MESH_PATH) == case["geometry"]["mesh_sha256"],
        "material_sha256_matches_case": sha256(MATERIAL_PATH) == case["material"]["material_sha256"],
        "runner_sha256_matches_started": started.get("runner_sha256") == sha256(EXP / "run_native_free_patch.py"),
        "contract_sha256_matches_started": started.get("contract_sha256") == sha256(EXP / "free_patch_contract.py"),
        "analyzer_sha256_matches_started": started.get("analyzer_sha256") == sha256(Path(__file__)),
        "test_sha256_matches_started": started.get("test_sha256") == sha256(EXP / "test_source_contract.py"),
        "launcher_sha256_matches_started": started.get("launcher_sha256") == sha256(EXP / "launch_free_patch.sh"),
        "stability_status_matches_case": started.get("total_stability_status") == total_stability_status,
        "stability_lambda_matches_case": started.get("lambda_max_upper_bound_per_s2") == lambda_upper,
    }
    source_integrity_ok = all(source_integrity.values())

    summary = {
        "analysis_status": "ANALYSIS_COMPLETE_DIAGNOSTIC; NO_PHYSICAL_OR_FRAGMENT_QUALIFICATION",
        "case_id": case["case_id"],
        "classification": case["classification"],
        "native_run_status": status.get("status"),
        "solver_steps_completed": status.get("solver_steps_completed"),
        "history_frames_saved": frame_count,
        "expected_history_frames": expected_frames,
        "run_directory": str(run_dir),
        "source_integrity": source_integrity,
        "gates": {
            "closed_initialized_plus_every_postsolve_history": history_complete,
            "frozen_input_and_source_bytes_verified": source_integrity_ok,
            "native_54_96_8_topology_unchanged": topology_stable,
            "exact_source_cell_ordinal_and_positive_parentfacet_lineage": lineage_match,
            "all_24_quadrature_values_saved_without_cell_averaging": field_shape_ok,
            "native_fields_and_raw_energy_arrays_finite": native_finite,
            "free_bulk_and_all_cohesive_side_dofs_contract": dofs_ok,
            "pre_and_post_solve_boundary_readback_le_1e-12m": boundary_ok,
            "prescribed_half_cosine_and_hold_contract": loading_ok,
            "bulk_only_stable_step_sanity_check": bulk_stability_ok,
            "parent_reviewed_total_cohesive_plus_bulk_stability_bound": total_stability_ok,
            "damage_finite_bounded_and_monotone_at_each_qp": damage_ok,
            "delta_max_finite_and_monotone_at_each_qp": delta_ok,
            "raw_native_dissipation_not_below_roundoff_floor": negative_dissipation_ok,
            "independent_signed_interface_work_matches_native_total_history": work_identity_ok,
        },
        "energy_and_sign_diagnostics": {
            "native_energy_keys": list(ENERGY_KEYS),
            "native_traction_convention": "energy-conjugate to native opening; raw t dot opening retained; no sign flip",
            "source_formula": "D = cumulative sum(0.5 * delta_opening dot (t + t_previous)) - 0.5 * t dot opening",
            "native_dissipation_root_cause_diagnosis": dissipation_diagnosis,
            "negative_dissipation_absolute_tolerance_J": negative_tolerance,
            "raw_dissipated_min_J": float(np.min(dissipated)),
            "raw_dissipated_min_history_index": int(np.argmin(dissipated)),
            "raw_negative_dissipation_sample_count": negative_dissipation_raw_count,
            "below_tolerance_dissipation_sample_count": negative_dissipation_count,
            "independent_signed_interface_work_final_J": float(independent_work_cumulative[-1]),
            "native_dissipated_plus_reversible_final_J_raw": float(native_total_energy[-1]),
            "max_abs_independent_work_vs_native_total_residual_J": closure_error_max,
            "max_abs_cumulative_work_array_consistency_error_J": cumulative_match_error,
            "history_identity_tolerance_J": energy_work_tolerance,
            "traction_dot_opening_Pa_m_raw_min": float(np.min(arrays["traction_dot_opening_Pa_m_raw"])),
            "traction_dot_opening_Pa_m_raw_max": float(np.max(arrays["traction_dot_opening_Pa_m_raw"])),
            "interface_work_increment_J_raw_min": float(np.min(native_work_increment)),
            "interface_work_increment_J_raw_max": float(np.max(native_work_increment)),
            "raw_energy_first_last_J": [energies[0].tolist(), energies[-1].tolist()],
            "actuator_energy_balance_residual_first_last_J_raw": [
                float(arrays["energy_balance_residual_J_raw"][0]),
                float(arrays["energy_balance_residual_J_raw"][-1]),
            ],
            "actuator_balance_is_diagnostic_only": True,
            "full_debond_state_by_declared_damage_and_hold_rule": full_debond,
        },
        "lineage": {
            "input_interface_faces": rebuilt_lineage["input_face_count"],
            "native_cohesive_elements": rebuilt_lineage["native_cohesive_element_count"],
            "positive_rest_area_sum_m2": rebuilt_lineage["positive_rest_area_sum_m2"],
            "parent_facets": rebuilt_lineage["parent_facets"],
            "native_cohesive_element_id_note": "local row index from native getConnectivity; not a persistent global ID",
        },
        "free_dof_contract": {
            "blocked_components": blocked_component_count,
            "free_components": int(initial_blocked.size - blocked_component_count),
            "free_cohesive_side_components": free_interface_components,
            "free_geometric_interior_bulk_components": free_interior_components,
        },
        "scope": case["scope_limits"],
        "interpretation": "A completed run is a bounded native coupling observation only. No fracture calibration, crack path, process-zone resolution, fragment, R2, or pretraining eligibility claim is made.",
    }
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_directory", nargs="?", default="/home/xiemo/free_cohesive_patch_20261015/free_native_intrinsic_modeI_patch_20261015_v1")
    args = parser.parse_args()
    run_dir = Path(args.run_directory).resolve()
    analysis_dir = run_dir / "analysis"
    if analysis_dir.exists() and any(analysis_dir.iterdir()):
        raise SystemExit(f"analysis output exists; refusing overwrite: {analysis_dir}")
    analysis_dir.mkdir(parents=True, exist_ok=True)
    summary = analyze(run_dir)
    summary["analysis_source_sha256"] = sha256(Path(__file__))
    summary["analysis_argv"] = sys.argv
    out = analysis_dir / "analysis.json"
    write_summary(out, summary)
    print(json.dumps({"analysis_status": summary["analysis_status"], "gates": summary["gates"],
                      "dissipation_diagnosis": summary["energy_and_sign_diagnostics"]["native_dissipation_root_cause_diagnosis"],
                      "output": str(out)}, allow_nan=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
