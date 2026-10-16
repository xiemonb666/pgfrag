"""Independent readback for the one authorized zero-load native benchmark."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[4]
PG = ROOT / "pgfrag"
EXPERIMENT = PG / "experiments/continuation_20261016/reference_vessel_benchmark"
DATA = PG / "data/continuation_20261016/reference_vessel_benchmark"
INPUTS = EXPERIMENT / "inputs"
RUN = EXPERIMENT / "native/run_attempt_01"
CASE_PATH = INPUTS / "case.json"
REPORT = DATA / "native_evolution_audit.json"
SUMMARY = DATA / "README.md"
ALL_OUTPUT_BUDGET = 256 * 1024 * 1024

LEDGER_COLUMNS = [
    "time_s", "load_factor_N", "interface_work_J", "external_work_J", "algorithmic_kinetic_J",
    "momentum_x", "momentum_y", "momentum_z", "angular_momentum_x", "angular_momentum_y", "angular_momentum_z",
    "external_torque_x_Nm", "external_torque_y_Nm", "external_torque_z_Nm",
    "ma_minus_f_internal_external_max_N", "position_recurrence_max_m", "velocity_recurrence_max_mps",
    "external_work_minus_total_alg_J", "max_damage", "min_damage", "internal_work_cumulative_J",
    "kinetic_alg_increment_identity_residual_J", "internal_work_plus_material_energy_J",
    "native_minus_nodal_kinetic_J", "internal_current_torque_x_Nm", "internal_current_torque_y_Nm",
    "internal_current_torque_z_Nm", "max_contact_traction_Pa",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def json_read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def json_write(path: Path, obj: dict) -> None:
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def output_bytes() -> int:
    roots = (EXPERIMENT, DATA)
    return sum(p.stat().st_size for root in roots for p in root.rglob("*") if p.is_file())


def main() -> int:
    checks: dict[str, bool] = {}
    details: dict[str, object] = {}

    def check(name: str, passed: bool, detail: object = None) -> None:
        checks[name] = bool(passed)
        if detail is not None:
            details[name] = detail

    report = {
        "schema": "pgfrag.native_cohesive_vessel_benchmark_audit.v1",
        "status": "AUDIT_RUNNING",
        "classification": "ZERO_LOAD_NATIVE_COST_BENCHMARK_ONLY; PLACEHOLDER_UNCALIBRATED; NO_PHYSICAL_OR_FRACTURE_QUALIFICATION",
        "checks": checks,
        "details": details,
    }
    try:
        if not RUN.is_dir():
            raise FileNotFoundError(f"native attempt directory is absent: {RUN}")
        case = json_read(CASE_PATH)
        budget = json_read(DATA / "preflight_budget.json")
        execution = json_read(EXPERIMENT / "execution_record.json")
        result = json_read(RUN / "result.json")
        started = json_read(RUN / "started.json")
        receipt = json_read(RUN / "receipt.json")
        completion = json_read(RUN / "completion.json") if (RUN / "completion.json").is_file() else {}
        checkpoint = json_read(RUN / "checkpoint.json") if (RUN / "checkpoint.json").is_file() else {}
        bound = json_read(RUN / "initial_bound.json") if (RUN / "initial_bound.json").is_file() else {}
        loading = json_read(RUN / "loading.json") if (RUN / "loading.json").is_file() else {}

        input_hashes = {p.name: sha256(p) for p in sorted(INPUTS.iterdir()) if p.is_file()}
        expected_input_hashes = {name: row["sha256"] for name, row in budget["input_freeze"]["frozen_files"].items()}
        check("frozen_input_hashes_match_preflight", input_hashes == expected_input_hashes,
              {"actual": input_hashes, "expected": expected_input_hashes})
        case_inputs_expected = {name: digest for name, digest in input_hashes.items() if name != "case.json"}
        check("wrapper_case_input_hashes_match", case["input_hashes"] == case_inputs_expected,
              {"case_input_names": sorted(case["input_hashes"]),
               "expected_noncase_input_names": sorted(case_inputs_expected)})
        check("execution_input_hashes_unchanged", execution.get("frozen_inputs_unchanged") is True
              and execution.get("frozen_input_hashes_before") == execution.get("frozen_input_hashes_after") == input_hashes,
              {"before_equals_after": execution.get("frozen_inputs_unchanged"), "returncode": execution.get("returncode"),
               "outer_elapsed_s": execution.get("outer_elapsed_s"), "timed_out": execution.get("timed_out")})
        check("wrapper_started_inventory_matches_inputs", started.get("source_input_inventory") == input_hashes,
              {"inventory_names": sorted(started.get("source_input_inventory", {})), "frozen_input_names": sorted(input_hashes)})
        check("no_forge_manifest_used", not (INPUTS / "forge_manifest.json").exists())
        check("source_geometry_ancestry_still_unbound",
              case["source_geometry"]["source_face_ancestry"]["status"] == "PRESENT_BUT_UNBOUND"
              and case["source_geometry"]["source_face_ancestry"]["exterior_triangle_counts"]["unknown"] == 35784)
        check("zero_load_case_parameters_preserved", float(case["loading"]["peak_force_N"]) == 0.0
              and case["time"]["total_s"] == 1.0e-7 and case["output"]["node_stride"] == 1
              and case["output"]["flush_stride"] == 10 and case["limits"]["max_steps"] == 256
              and case["limits"]["wall_time_s"] == 120 and case["limits"]["output_cap_bytes"] == 128 * 1024 * 1024)
        check("execution_process_completed_once", execution.get("status") == "PROCESS_EXITED_PRESERVED"
              and execution.get("returncode") == 0 and execution.get("timed_out") is False
              and execution.get("no_retry_or_resume") is True,
              {"status": execution.get("status"), "returncode": execution.get("returncode"),
               "outer_elapsed_s": execution.get("outer_elapsed_s"), "timed_out": execution.get("timed_out")})
        check("result_reports_native_completion", result.get("status") == "NATIVE_RUN_COMPLETE_DIAGNOSTIC"
              and result.get("source_inputs_unchanged") is True,
              {"status": result.get("status"), "error": result.get("error"), "traceback": result.get("traceback")})

        expected_counts = [18634, 56139, 541]
        steps = int(budget["time_forecast"]["solver_steps"])
        frames = steps + 1
        dt = float(budget["time_forecast"]["actual_dt_s"])
        result_counts = result.get("counts")
        check("native_mesh_counts_match_fullsize_v5", result_counts == expected_counts,
              {"counts": result_counts, "expected": expected_counts})
        check("time_steps_match_preflight_forecast", result.get("solver_steps") == steps
              and int(bound.get("steps", -1)) == steps and math.isclose(float(result.get("dt_s", math.nan)), dt, rel_tol=0, abs_tol=0)
              and math.isclose(float(bound.get("dt_s", math.nan)), dt, rel_tol=0, abs_tol=0)
              and math.isclose(steps * dt, float(case["time"]["total_s"]), rel_tol=2e-15, abs_tol=1e-22),
              {"steps": result.get("solver_steps"), "dt_s": result.get("dt_s"), "total_s": steps * dt,
               "forecast_steps": steps, "forecast_dt_s": dt, "initial_bound_steps": bound.get("steps")})
        check("completion_and_final_checkpoint_are_consistent",
              completion.get("steps_completed") == steps and completion.get("qp_frames_written") == frames
              and completion.get("node_frames_written") == frames
              and checkpoint.get("completed_steps") == steps and checkpoint.get("qp_frames_committed") == frames
              and checkpoint.get("node_frames_committed") == frames,
              {"completion": completion, "checkpoint": checkpoint})
        check("reported_history_forecast_matches_budget",
              result.get("raw_history_bytes_forecast") == budget["exact_history_forecast"]["wrapper_forecast_bytes_including_256_each_for_12_arrays"],
              {"result_forecast_bytes": result.get("raw_history_bytes_forecast"),
               "budget_forecast_bytes": budget["exact_history_forecast"]["wrapper_forecast_bytes_including_256_each_for_12_arrays"]})

        # The wrapper's receipt intentionally excludes itself from its own file
        # inventory. Recompute every other file byte count and digest here.
        actual_files = {
            str(p.relative_to(RUN)).replace("\\", "/"):
                {"bytes": p.stat().st_size, "sha256": sha256(p)}
            for p in sorted(RUN.rglob("*")) if p.is_file() and p.name != "receipt.json"
        }
        check("native_receipt_inventory_matches_files", receipt.get("files") == actual_files,
              {"inventory_entries": len(receipt.get("files", {})), "actual_entries": len(actual_files),
               "persistent_bytes_before_receipt": receipt.get("persistent_bytes_before_receipt"),
               "actual_persistent_bytes": sum(row["bytes"] for row in actual_files.values())})
        receipt_digest_ok = all(
            execution.get(field, {}).get("sha256") == sha256(EXPERIMENT / filename)
            and execution.get(field, {}).get("bytes") == (EXPERIMENT / filename).stat().st_size
            for field, filename in (("stdout", "native_stdout.log"), ("stderr", "native_stderr.log"))
        )
        check("captured_native_stdout_stderr_match_execution_record", receipt_digest_ok,
              {"stdout": execution.get("stdout"), "stderr": execution.get("stderr")})
        check("wrapper_receipt_status_matches_result", receipt.get("status") == result.get("status"))

        # Independent readback of all full-QP histories and every saved node
        # frame, using memory maps so verification does not duplicate 94 MB.
        hist_dir = RUN / "history"
        float_shapes = {
            "opening_m": (frames, 1623, 3), "traction_Pa": (frames, 1623, 3),
            "damage": (frames, 1623), "delta_max_m": (frames, 1623),
            "contact_traction_Pa": (frames, 1623, 3), "contact_opening_m": (frames, 1623, 3),
            "contact_energy_check": (frames, 3), "interface_work_check": (frames, 4),
            "energy_J": (frames, 5), "ledger": (frames, 28),
            "patch_displacement_x_m": (frames, int(budget["native_counts"]["balanced_patch_nodes"])),
            "node_state": (frames, 18634, 4, 3),
        }
        history: dict[str, np.ndarray] = {}
        array_report: dict[str, dict] = {}
        finite_all = True
        shape_dtype_all = True
        history_zero_metrics: dict[str, float] = {}
        history_files = {p.name for p in hist_dir.glob("*.npy")} if hist_dir.is_dir() else set()
        expected_history_files = {f"{name}.npy" for name in float_shapes} | {"node_saved_steps.npy"}
        shape_dtype_all = history_files == expected_history_files
        for name, shape in float_shapes.items():
            path = hist_dir / f"{name}.npy"
            if not path.is_file():
                array_report[name] = {"present": False}
                finite_all = shape_dtype_all = False
                continue
            arr = np.load(path, mmap_mode="r", allow_pickle=False)
            history[name] = arr
            matches = arr.shape == shape and arr.dtype == np.dtype(np.float64)
            shape_dtype_all &= bool(matches)
            finite = bool(np.isfinite(arr).all()) if matches else False
            finite_all &= finite
            array_report[name] = {"shape": list(arr.shape), "dtype": arr.dtype.str, "bytes": path.stat().st_size,
                                  "shape_dtype_match": bool(matches), "all_finite": finite}
            if matches and finite:
                history_zero_metrics[name] = float(np.max(np.abs(arr)))
        node_steps_path = hist_dir / "node_saved_steps.npy"
        if node_steps_path.is_file():
            node_steps = np.load(node_steps_path, mmap_mode="r", allow_pickle=False)
            expected_node_steps = np.arange(frames, dtype=np.int64)
            saved_steps_ok = node_steps.dtype == np.dtype(np.int64) and np.array_equal(node_steps, expected_node_steps)
            shape_dtype_all &= bool(saved_steps_ok)
            array_report["node_saved_steps"] = {"shape": list(node_steps.shape), "dtype": node_steps.dtype.str,
                                                "match_all_solver_steps": bool(saved_steps_ok)}
        else:
            saved_steps_ok = False
            shape_dtype_all = False
            array_report["node_saved_steps"] = {"present": False}
        layout = json_read(hist_dir / "layout.json") if (hist_dir / "layout.json").is_file() else {}
        check("all_history_arrays_have_expected_shape_dtype_and_step_indices", shape_dtype_all and saved_steps_ok,
              {"arrays": array_report, "expected_shapes": {k: list(v) for k, v in float_shapes.items()},
               "expected_step_indices": [0, steps]})
        check("all_history_arrays_are_finite", finite_all, {"array_finite": {k: v.get("all_finite") for k, v in array_report.items()}})
        check("history_layout_preserves_native_qp_order_and_full_qps",
              layout.get("qp_order") == "native element-major, three quadrature points per cohesive_3d_6; no averaging"
              and layout.get("all_qp_states_retained") is True and layout.get("node_stride") == 1
              and layout.get("ledger_columns") == LEDGER_COLUMNS,
              {"qp_order": layout.get("qp_order"), "all_qp_states_retained": layout.get("all_qp_states_retained"),
               "node_stride": layout.get("node_stride"), "ledger_column_count": len(layout.get("ledger_columns", []))})

        failure_steps_path = RUN / "qp_first_full_damage_step.npy"
        if failure_steps_path.is_file():
            qp_first_failure = np.load(failure_steps_path, mmap_mode="r", allow_pickle=False)
            no_full_damage = qp_first_failure.dtype == np.dtype(np.int64) and qp_first_failure.shape == (1623,) and bool(np.all(qp_first_failure == -1))
        else:
            qp_first_failure = None
            no_full_damage = False
        check("all_quadrature_points_remain_unfailed", no_full_damage,
              {"present": qp_first_failure is not None,
               "shape": list(qp_first_failure.shape) if qp_first_failure is not None else None,
               "dtype": qp_first_failure.dtype.str if qp_first_failure is not None else None,
               "unique_first_full_damage_steps": np.unique(qp_first_failure).tolist() if qp_first_failure is not None else None})

        zero_tolerances = {
            "opening_m": 1e-12, "traction_Pa": 1e-2, "damage": 1e-15,
            "contact_traction_Pa": 1e-2,
            "contact_opening_m": 1e-12, "contact_energy_check": 1e-15,
            "interface_work_check": 1e-15, "energy_J": 1e-15,
            "patch_displacement_x_m": 1e-12,
        }
        zero_checks = {name: (name in history_zero_metrics and history_zero_metrics[name] <= tol)
                       for name, tol in zero_tolerances.items()}
        # Node state concatenates displacement, velocity, acceleration, and
        # internal force components; use units appropriate to each channel.
        node = history.get("node_state")
        node_state_limits = {"displacement_m": 1e-12, "velocity_mps": 1e-9,
                             "acceleration_mps2": 1e-4, "internal_force_N": 1e-6}
        node_state_maxima = {}
        node_states_zero = node is not None and node.shape == float_shapes["node_state"]
        if node_states_zero:
            for index, (name, tolerance) in enumerate(node_state_limits.items()):
                maximum = float(np.max(np.abs(node[:, :, index, :])))
                node_state_maxima[name] = maximum
                node_states_zero &= maximum <= tolerance
        zero_checks["node_state"] = bool(node_states_zero)
        zero_tolerances["node_state"] = node_state_limits
        delta0 = float(case["interface"]["delta0_m"])
        delta_max_difference = (float(np.max(np.abs(history["delta_max_m"] - delta0)))
                               if "delta_max_m" in history else math.inf)
        delta_max_at_initial = math.isfinite(delta_max_difference) and delta_max_difference <= 1e-15
        zero_checks["delta_max_m_is_unchanged_initial_history_threshold"] = delta_max_at_initial
        ledger = history.get("ledger")
        non_time_ledger_max = float(np.max(np.abs(ledger[:, 1:]))) if ledger is not None else math.inf
        ledger_zero = math.isfinite(non_time_ledger_max) and non_time_ledger_max <= 1e-15
        zero_checks["non_time_ledger_columns_zero"] = ledger_zero
        check("zero_load_qp_energetic_and_node_states_within_explicit_tolerances", all(zero_checks.values()),
              {"maximum_absolute_by_array": history_zero_metrics, "tolerances_by_array": zero_tolerances,
               "node_state_channel_maxima": node_state_maxima,
               "delta_max_m_initial_threshold": delta0,
               "delta_max_m_max_abs_difference_from_initial_threshold": delta_max_difference,
               "ledger_non_time_columns_max_abs": non_time_ledger_max,
               "ledger_time_column_is_checked_separately": True,
               "per_array_pass": zero_checks})

        state = history.get("node_state")
        initial_state_path = RUN / "initial_state.npz"
        if initial_state_path.is_file() and ledger is not None and state is not None and state.shape == float_shapes["node_state"]:
            with np.load(initial_state_path, allow_pickle=False) as z:
                rest = np.asarray(z["rest_nodes_m"], dtype=np.float64)
                mass = np.asarray(z["mass_kg"], dtype=np.float64)
                external_unit = np.asarray(z["unit_external_load"], dtype=np.float64)
                blocked = np.asarray(z["blocked"], dtype=bool)
            shapes_ok = rest.shape == (18634, 3) and mass.shape == (18634, 3) and external_unit.shape == (18634, 3) and not blocked.any()
            force_error = []
            momentum_error = []
            angular_error = []
            for frame in range(frames):
                disp = np.asarray(state[frame, :, 0, :])
                vel = np.asarray(state[frame, :, 1, :])
                acc = np.asarray(state[frame, :, 2, :])
                internal = np.asarray(state[frame, :, 3, :])
                # peak_force_N is zero, so the actual nodal external force is exactly zero.
                residual = mass * acc - internal - (float(case["loading"]["peak_force_N"]) * external_unit)
                force_error.append(float(np.max(np.abs(residual))))
                momentum = np.sum(mass * vel, axis=0)
                angular = np.sum(np.cross(rest + disp, mass * vel), axis=0)
                momentum_error.append(float(np.max(np.abs(momentum - ledger[frame, 5:8]))))
                angular_error.append(float(np.max(np.abs(angular - ledger[frame, 8:11]))))
            force_error = np.asarray(force_error)
            momentum_error = np.asarray(momentum_error)
            angular_error = np.asarray(angular_error)
            independent_force_ok = bool(np.isfinite(force_error).all() and force_error.max(initial=0.0) <= 1e-6)
            momentum_ok = bool(np.isfinite(momentum_error).all() and momentum_error.max(initial=0.0) <= 1e-12)
            angular_ok = bool(np.isfinite(angular_error).all() and angular_error.max(initial=0.0) <= 1e-12)
            ledger_force_matches = bool(np.allclose(ledger[:, 14], force_error, rtol=2e-12, atol=1e-15))
            times_ok = bool(np.allclose(ledger[:, 0], np.arange(frames) * dt, rtol=0, atol=1e-22)
                            and np.all(ledger[:, 1] == 0.0))
            force_report = {"mass_state_shapes_ok": bool(shapes_ok),
                            "max_abs_nodal_force_balance_residual_N": float(force_error.max(initial=0.0)),
                            "max_abs_momentum_vs_ledger_kg_m_s": float(momentum_error.max(initial=0.0)),
                            "max_abs_angular_momentum_vs_ledger_kg_m2_s": float(angular_error.max(initial=0.0)),
                            "ledger_force_residual_matches_independent": ledger_force_matches,
                            "momentum_limit_kg_m_s": 1e-12, "angular_momentum_limit_kg_m2_s": 1e-12,
                            "force_balance_limit_N": 1e-6,
                            "reference_unit_load_sum_N": loading.get("reference_unit_force_sum"),
                            "reference_unit_torque_sum_m": loading.get("reference_unit_torque_sum_m"),
                            "actual_external_load_is_zero": float(case["loading"]["peak_force_N"]) == 0.0,
                            "times_and_load_factor_match": times_ok}
            check("independent_nodal_force_momentum_and_angular_momentum_audit_passes",
                  shapes_ok and independent_force_ok and momentum_ok and angular_ok and ledger_force_matches and times_ok,
                  force_report)
        else:
            check("independent_nodal_force_momentum_and_angular_momentum_audit_passes", False,
                  {"reason": "required full node history, initial state, or ledger is missing"})

        # Cost and completion values are observations, not a full-run ETA claim.
        actual_native_bytes = sum(row["bytes"] for row in actual_files.values())
        experiment_bytes_before_audit_report = output_bytes()
        audit_summary = {
            "native_result_elapsed_s": result.get("elapsed_s"),
            "outer_process_elapsed_s": execution.get("outer_elapsed_s"),
            "native_steps": result.get("solver_steps"),
            "native_steps_per_second_observed": (result.get("solver_steps", 0) / result["elapsed_s"]
                                                  if result.get("elapsed_s", 0) else None),
            "counts": result.get("counts"), "dt_s": result.get("dt_s"),
            "actual_native_run_bytes_ex_receipt": actual_native_bytes,
            "receipt_bytes": (RUN / "receipt.json").stat().st_size if (RUN / "receipt.json").exists() else None,
            "output_bytes_across_benchmark_experiment_and_data": experiment_bytes_before_audit_report,
            "declared_all_output_budget_bytes": ALL_OUTPUT_BUDGET,
            "whole_benchmark_budget_remaining_before_audit_report": ALL_OUTPUT_BUDGET - experiment_bytes_before_audit_report,
            "no_full_run_eta_claim": True,
        }
        details["observed_native_cost"] = audit_summary
        report["timings_and_cost"] = audit_summary
        report["array_readback"] = array_report
        report["checks"] = checks
        report["details"] = details
        report["status"] = "ZERO_LOAD_NATIVE_EVOLUTION_NUMERICAL_PASS" if all(checks.values()) else "ZERO_LOAD_NATIVE_EVOLUTION_NUMERICAL_FAIL_PRESERVED"
    except Exception as exc:
        report["status"] = "ZERO_LOAD_NATIVE_EVOLUTION_AUDIT_ERROR_PRESERVED"
        report["audit_error"] = repr(exc)
        report["audit_traceback"] = __import__("traceback").format_exc()
        report["checks"] = checks
        report["details"] = details

    json_write(REPORT, report)
    all_bytes = output_bytes()
    report["output_bytes_across_benchmark_experiment_and_data_including_audit"] = all_bytes
    report["declared_all_output_budget_bytes"] = ALL_OUTPUT_BUDGET
    report["all_output_budget_passes"] = all_bytes <= ALL_OUTPUT_BUDGET
    if all_bytes > ALL_OUTPUT_BUDGET:
        report["status"] = "ZERO_LOAD_NATIVE_EVOLUTION_AUDIT_FAIL_OUTPUT_BUDGET_EXCEEDED"
    json_write(REPORT, report)
    # The report changed only by its self-inclusive total field; recalculate the
    # exact final total after the second serialization and update once more.
    all_bytes = output_bytes()
    report["output_bytes_across_benchmark_experiment_and_data_including_audit"] = all_bytes
    report["all_output_budget_passes"] = all_bytes <= ALL_OUTPUT_BUDGET
    json_write(REPORT, report)

    md = [
        "# Reference vessel zero-load native benchmark audit",
        "",
        f"Status: **{report['status']}**.",
        "",
        "This is a bounded full-size v5 diagnostic-mesh cost/numerical-control benchmark with zero applied force. The bulk/interface values are uncalibrated placeholders; exterior ancestry remains unbound. It is not a physical qualification, fracture result, training result, or full-run cost guarantee.",
        "",
        f"The native wrapper reported {audit_summary.get('native_steps')} steps at {audit_summary.get('dt_s')} s, with {audit_summary.get('native_result_elapsed_s')} s inside the wrapper and {audit_summary.get('outer_process_elapsed_s')} s wall time for the outer process. These are observations from this short zero-load run, not a projected duration for 50 microseconds.",
        "",
        f"The final byte count across the two benchmark folders is {all_bytes:,} / {ALL_OUTPUT_BUDGET:,} bytes. Every check is recorded in [`native_evolution_audit.json`](native_evolution_audit.json).",
        "",
        "The attempt is preserved as run_attempt_01. No retry, resume, loading change, or tuning was performed.",
    ]
    # Iterate the self-inclusive byte total until both report and summary show
    # the exact on-disk total. Their numeric fields have stable width at this
    # benchmark size, so this normally converges on the first correction.
    all_bytes = output_bytes()
    for _ in range(8):
        report["output_bytes_across_benchmark_experiment_and_data_including_audit"] = all_bytes
        report["declared_all_output_budget_bytes"] = ALL_OUTPUT_BUDGET
        report["all_output_budget_passes"] = all_bytes <= ALL_OUTPUT_BUDGET
        if all_bytes > ALL_OUTPUT_BUDGET:
            report["status"] = "ZERO_LOAD_NATIVE_EVOLUTION_AUDIT_FAIL_OUTPUT_BUDGET_EXCEEDED"
        json_write(REPORT, report)
        md[8] = f"The final byte count across the two benchmark folders is {all_bytes:,} / {ALL_OUTPUT_BUDGET:,} bytes. Every check is recorded in [`native_evolution_audit.json`](native_evolution_audit.json)."
        SUMMARY.write_text("\n".join(md) + "\n", encoding="utf-8")
        actual = output_bytes()
        if actual == all_bytes:
            break
        all_bytes = actual
    all_bytes = output_bytes()
    report["output_bytes_across_benchmark_experiment_and_data_including_audit"] = all_bytes
    report["all_output_budget_passes"] = all_bytes <= ALL_OUTPUT_BUDGET
    if all_bytes > ALL_OUTPUT_BUDGET:
        report["status"] = "ZERO_LOAD_NATIVE_EVOLUTION_AUDIT_FAIL_OUTPUT_BUDGET_EXCEEDED"
    json_write(REPORT, report)
    print(json.dumps({"status": report["status"], "checks_passed": sum(checks.values()),
                      "checks_total": len(checks), "output_bytes": all_bytes,
                      "report": str(REPORT)}, separators=(",", ":")))
    return 0 if report["status"] == "ZERO_LOAD_NATIVE_EVOLUTION_NUMERICAL_PASS" and report["all_output_budget_passes"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
