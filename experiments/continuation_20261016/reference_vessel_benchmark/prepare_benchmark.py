"""Freeze and budget one zero-load native vessel evolution benchmark.

This preparation step copies only immutable executed inputs into this new
benchmark package, changes the requested time/resource controls, and checks
the exact NPY history forecast before the native driver is allowed to run.
It does not import Akantu or execute solver evolution.
"""
from __future__ import annotations

import hashlib
import io
import json
import math
from pathlib import Path
import shutil
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[4]
PG = ROOT / "pgfrag"
EXPERIMENT = PG / "experiments/continuation_20261016/reference_vessel_benchmark"
DATA = PG / "data/continuation_20261016/reference_vessel_benchmark"
INPUTS = EXPERIMENT / "inputs"
CASE_PATH = INPUTS / "case.json"
PREVIOUS = PG / "experiments/continuation_20261016/reference_vessel_preflight/run_attempt_02"
PREVIOUS_INPUTS = PREVIOUS / "inputs"
PREVIOUS_SOURCE_ARCHIVE = PREVIOUS / "source_archive"
PREVIOUS_RUN = PREVIOUS / "native_preflight/run01"
REPORT_PATH = DATA / "preflight_budget.json"
README_PATH = DATA / "README.md"
ALL_OUTPUT_BUDGET = 256 * 1024 * 1024
HISTORY_CAP = 128 * 1024 * 1024
NON_HISTORY_OUTPUT_ALLOWANCE = 4_000_000


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def npy_header_bytes(shape: tuple[int, ...], dtype: np.dtype) -> int:
    stream = io.BytesIO()
    np.lib.format.write_array_header_1_0(stream, {
        "descr": np.dtype(dtype).str, "fortran_order": False, "shape": shape,
    })
    return stream.tell()


def expected_layout_bytes(node_stride: int) -> int:
    ledger_names = [
        "time_s", "load_factor_N", "interface_work_J", "external_work_J", "algorithmic_kinetic_J",
        "momentum_x", "momentum_y", "momentum_z", "angular_momentum_x", "angular_momentum_y", "angular_momentum_z",
        "external_torque_x_Nm", "external_torque_y_Nm", "external_torque_z_Nm",
        "ma_minus_f_internal_external_max_N", "position_recurrence_max_m", "velocity_recurrence_max_mps",
        "external_work_minus_total_alg_J", "max_damage", "min_damage", "internal_work_cumulative_J",
        "kinetic_alg_increment_identity_residual_J", "internal_work_plus_material_energy_J",
        "native_minus_nodal_kinetic_J", "internal_current_torque_x_Nm", "internal_current_torque_y_Nm",
        "internal_current_torque_z_Nm", "max_contact_traction_Pa",
    ]
    layout = {
        "ledger_columns": ledger_names,
        "energy_columns": ["kinetic", "potential", "dissipated", "reversible", "cohesive contact"],
        "node_state_columns": ["displacement_m", "velocity_mps", "acceleration_mps2", "internal_force_N"],
        "contact_energy_check_columns": ["independent_quadrature_contact_energy_J", "native_contact_energy_J", "independent_minus_native_J"],
        "interface_work_check_columns": ["full_interface_work_J", "cohesive_work_J", "contact_work_J", "cross_work_J"],
        "qp_order": "native element-major, three quadrature points per cohesive_3d_6; no averaging",
        "node_stride": node_stride, "all_qp_states_retained": True,
    }
    return len((json.dumps(layout, indent=2, allow_nan=False) + "\n").encode("utf-8"))


def dir_bytes(path: Path) -> int:
    return sum(p.stat().st_size for p in path.rglob("*") if p.is_file())


def main() -> None:
    if INPUTS.exists() or (EXPERIMENT / "native/run_attempt_01").exists() or REPORT_PATH.exists():
        raise FileExistsError("benchmark attempt or prepared inputs already exist; this task permits one attempt only")
    DATA.mkdir(parents=True, exist_ok=True)
    INPUTS.mkdir(parents=True)

    # Preserve the exact current source bytes used by the earlier successful
    # zero-step preflight, plus the exact physical-name alias mesh/material.
    mappings = {
        "cohesive_vessel.py": PG / "src/pgfrag/cohesive_vessel.py",
        "cohesive_mesh.py": PG / "src/pgfrag/cohesive_mesh.py",
        "cohesive_stability.py": PG / "src/pgfrag/cohesive_stability.py",
        "mesh.msh": PREVIOUS_INPUTS / "mesh.msh",
        "material.dat": PREVIOUS_INPUTS / "material.dat",
        "mesh_alias_receipt.json": PREVIOUS_INPUTS / "mesh_alias_receipt.json",
        "mesh.original.msh": PREVIOUS_SOURCE_ARCHIVE / "mesh.msh.original",
    }
    for name, source in mappings.items():
        shutil.copyfile(source, INPUTS / name)
    for module in ("cohesive_vessel.py", "cohesive_mesh.py", "cohesive_stability.py"):
        previous_copy = PREVIOUS_INPUTS / module
        if sha256(INPUTS / module) != sha256(previous_copy):
            raise ValueError(f"current source does not match the prior frozen executed source: {module}")

    alias_receipt = json.loads((INPUTS / "mesh_alias_receipt.json").read_text(encoding="utf-8"))
    if alias_receipt.get("schema") != "pgfrag.physical_name_alias_receipt.v1":
        raise ValueError("unexpected physical-name alias receipt schema")
    if not all(alias_receipt.get(key) is True for key in (
        "nodes_block_byte_equal", "elements_block_byte_equal", "all_bytes_outside_physical_names_byte_equal",
        "physical_integer_tags_and_order_equal",
    )):
        raise ValueError("the alias receipt does not establish an unchanged mesh body and tags")
    if sha256(INPUTS / "mesh.msh") != alias_receipt["executed_mesh"]["sha256"]:
        raise ValueError("copied alias mesh does not match its existing receipt")
    if sha256(INPUTS / "mesh.original.msh") != alias_receipt["source_mesh"]["sha256"]:
        raise ValueError("copied original source mesh does not match the alias receipt")

    old_case = json.loads((PREVIOUS_INPUTS / "case.json").read_text(encoding="utf-8"))
    if old_case.get("schema") != "pgfrag.native_cohesive_vessel_case.v1":
        raise ValueError("unexpected preflight case schema")
    case = json.loads(json.dumps(old_case))
    case["case_name"] = "reference_vessel_v5_zero_load_100ns_native_benchmark"
    case["time"]["total_s"] = 1.0e-7
    case["limits"].update({
        "max_steps": 256,
        "wall_time_s": 120,
        "output_cap_bytes": HISTORY_CAP,
    })
    case["output"].update({"node_stride": 1, "flush_stride": 10, "all_native_qp_steps": True})
    if case["loading"]["peak_force_N"] != 0.0:
        raise ValueError("benchmark must preserve the declared zero load")
    case["limitations"] = [
        "v5 geometry/resource benchmark only; source is a diagnostic v5 candidate",
        "exterior source-face ancestry is PRESENT_BUT_UNBOUND (35784 unknown exterior triangles); no fragment or surface-paint qualification",
        "bulk/interface values remain uncalibrated placeholders; no specimen physical data",
        "100 ns zero-load evolution measures initialization and short-run native cost only; no fracture or physical qualification",
        "initial/fixed-normal stability estimate is not a full-trajectory guarantee",
    ]
    case["input_hashes"] = {name: sha256(INPUTS / name) for name in mappings}
    CASE_PATH.write_text(json.dumps(case, indent=2, allow_nan=False) + "\n", encoding="utf-8")

    total_s = float(case["time"]["total_s"])
    initial_fraction = float(case["time"]["initial_bound_fraction"])
    with (PREVIOUS_RUN / "initial_bound.json").open(encoding="utf-8") as stream:
        old_bound = json.load(stream)
    maximum_dt = min(initial_fraction * float(old_bound["critical_time_step_s"]), 0.4 * float(old_bound["bulk_only_native_dt_s"]))
    steps = math.ceil(total_s / maximum_dt) * int(case["time"].get("refinement_factor", 1))
    if steps > int(case["limits"]["max_steps"]):
        raise ValueError("declared 256-step cap is insufficient; refusing to loosen it automatically")
    dt = total_s / steps
    frames = steps + 1
    node_count, qp_count = 18634, 1623
    patch_nodes = 3708
    saved_steps = np.unique(np.r_[np.arange(0, steps + 1, case["output"]["node_stride"]), steps]).astype(np.int64)
    float_specs = {
        "opening_m": (frames, qp_count, 3), "traction_Pa": (frames, qp_count, 3),
        "damage": (frames, qp_count), "delta_max_m": (frames, qp_count),
        "contact_traction_Pa": (frames, qp_count, 3), "contact_opening_m": (frames, qp_count, 3),
        "contact_energy_check": (frames, 3), "interface_work_check": (frames, 4),
        "energy_J": (frames, 5), "ledger": (frames, 28),
        "patch_displacement_x_m": (frames, patch_nodes),
        "node_state": (len(saved_steps), node_count, 4, 3),
    }
    raw_history_bytes = sum(int(np.prod(shape, dtype=np.int64)) * np.dtype(np.float64).itemsize for shape in float_specs.values())
    wrapper_forecast = sum(int(np.prod(shape, dtype=np.int64)) * np.dtype(np.float64).itemsize + 256 for shape in float_specs.values())
    step_stream = io.BytesIO()
    np.save(step_stream, saved_steps, allow_pickle=False)
    saved_steps_bytes = len(step_stream.getvalue())
    layout_bytes = expected_layout_bytes(int(case["output"]["node_stride"]))
    history_folder_bytes = wrapper_forecast + saved_steps_bytes + layout_bytes
    wrapper_cap_check_bytes = wrapper_forecast + 2_000_000
    if wrapper_cap_check_bytes > int(case["limits"]["output_cap_bytes"]):
        raise ValueError(f"exact full-history forecast {wrapper_cap_check_bytes} exceeds the 128 MiB case cap; no run attempted")

    previous_receipt = json.loads((PREVIOUS_RUN / "receipt.json").read_text(encoding="utf-8"))
    previous_static_output = int(previous_receipt["persistent_bytes_before_receipt"])
    failure_summary_shape = (qp_count,)
    failure_summary_bytes = int(np.prod(failure_summary_shape)) * np.dtype(np.int64).itemsize + npy_header_bytes(failure_summary_shape, np.dtype(np.int64))
    prepared_bytes = dir_bytes(EXPERIMENT) + dir_bytes(DATA)
    estimated_native_outputs = previous_static_output + history_folder_bytes + failure_summary_bytes
    predicted_all_bytes = prepared_bytes + estimated_native_outputs + NON_HISTORY_OUTPUT_ALLOWANCE
    if predicted_all_bytes > ALL_OUTPUT_BUDGET:
        raise ValueError("all-output forecast including inputs exceeds 256 MiB; no run attempted")

    report = {
        "schema": "pgfrag.native_cohesive_vessel_benchmark_preflight_budget.v1",
        "status": "EXACT_HISTORY_AND_ALL_OUTPUT_BUDGET_PASS; RUN_NOT_YET_STARTED",
        "solver_evolution_performed": False,
        "classification": "PLACEHOLDER_COST_ESTIMATE_UNCALIBRATED; NO SPECIMEN PHYSICAL DATA",
        "input_freeze": {
            "frozen_files": {name: {"sha256": sha256(INPUTS / name), "bytes": (INPUTS / name).stat().st_size} for name in sorted([*mappings.keys(), "case.json"])},
            "current_modules_match_previous_executed_bundle": True,
            "alias_mesh_sha256": sha256(INPUTS / "mesh.msh"),
            "original_source_mesh_sha256": sha256(INPUTS / "mesh.original.msh"),
            "alias_receipt_sha256": sha256(INPUTS / "mesh_alias_receipt.json"),
            "physical_name_adapter_only": alias_receipt["operation"] == "REPLACE_TWO_PHYSICAL_NAME_STRINGS_ONLY" and alias_receipt["all_bytes_outside_physical_names_byte_equal"],
            "no_forge_manifest_copied": not (INPUTS / "forge_manifest.json").exists(),
            "source_geometry_candidate": case["source_geometry"]["candidate"],
            "source_v5_exterior_ancestry_status": case["source_geometry"]["source_face_ancestry"]["status"],
            "unknown_exterior_triangles": case["source_geometry"]["source_face_ancestry"]["exterior_triangle_counts"]["unknown"],
        },
        "case_changes_from_prior_zero_step_preflight": {
            "time_total_s": total_s,
            "history_node_stride": int(case["output"]["node_stride"]),
            "flush_stride": int(case["output"]["flush_stride"]),
            "max_steps": int(case["limits"]["max_steps"]),
            "wrapper_wall_time_s": int(case["limits"]["wall_time_s"]),
            "history_output_cap_bytes": int(case["limits"]["output_cap_bytes"]),
            "peak_force_N": float(case["loading"]["peak_force_N"]),
            "unchanged_physics_material_load_patch_and_stability_fraction": True,
        },
        "native_counts": {"nodes_after_seam_insertion": node_count, "bulk_T4": 56139, "cohesive_Tri3": 541,
                          "interface_QP": qp_count, "balanced_patch_nodes": patch_nodes},
        "time_forecast": {"maximum_dt_s": maximum_dt, "actual_dt_s": dt, "solver_steps": steps,
                          "frames_including_initial": frames, "node_state_frames": len(saved_steps)},
        "exact_history_forecast": {
            "float64_fields_and_shapes": {name: list(shape) for name, shape in float_specs.items()},
            "node_saved_steps_shape": [len(saved_steps)], "node_saved_steps_dtype": np.dtype(np.int64).str,
            "raw_numeric_array_bytes": raw_history_bytes,
            "wrapper_forecast_bytes_including_256_each_for_12_arrays": wrapper_forecast,
            "node_saved_steps_npy_bytes": saved_steps_bytes,
            "layout_json_bytes": layout_bytes,
            "history_folder_estimate_bytes": history_folder_bytes,
            "wrapper_cap_check_forecast_plus_2mb_bytes": wrapper_cap_check_bytes,
            "history_cap_bytes": int(case["limits"]["output_cap_bytes"]),
            "history_cap_check_passes": True,
        },
        "all_output_budget": {
            "budget_bytes": ALL_OUTPUT_BUDGET,
            "frozen_inputs_and_all_benchmark_scripts_bytes_before_budget_reports": prepared_bytes,
            "static_native_outputs_estimated_from_same_geometry_preflight_bytes": previous_static_output,
            "qp_first_full_damage_summary_bytes": failure_summary_bytes,
            "history_folder_estimate_bytes": history_folder_bytes,
            "additional_runtime_receipt_log_and_small_json_allowance_bytes": NON_HISTORY_OUTPUT_ALLOWANCE,
            "estimated_all_bytes": predicted_all_bytes,
            "remaining_bytes": ALL_OUTPUT_BUDGET - predicted_all_bytes,
            "check_passes": True,
            "scope_note": "The exact prepared byte count covers frozen input files and all benchmark scripts present before this budget report is written. The 4,000,000-byte allowance covers this budget report and README, runtime JSON/receipts/logs, and small audit report. It also includes the prior same-geometry static artifact size, all forecast histories, and first-damage summary. No wall-time guarantee is implied.",
        },
        "source_geometry_limitations": {
            "geometry_source": case["source_geometry"]["candidate"],
            "exterior_ancestry": case["source_geometry"]["source_face_ancestry"]["status"],
            "unknown_exterior_triangle_count": case["source_geometry"]["source_face_ancestry"]["exterior_triangle_counts"]["unknown"],
            "fragments_or_surface_paint_qualified": False,
        },
        "attempt_policy": "Exactly one native run_attempt_01; no restart, retuning, force change, or automatic retry.",
    }
    REPORT_PATH.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    md = f"""# Reference vessel zero-load native benchmark preflight

This frozen case is prepared for one short Akantu evolution run. The real v5 alias mesh contains 56,139 T4 and 541 cohesive Tri3 facets (1,623 QPs), with 18,634 nodes after seam insertion. It preserves the prior placeholder material/interface values and zero peak force. No specimen data or fracture claim is implied.

The original source-v5 metadata still reports exterior ancestry `PRESENT_BUT_UNBOUND` with 35,784 unknown exterior triangles. The copied mesh is the exact previously audited physical-name alias; the alias receipt and copied original mesh bind its bytes. No `forge_manifest.json` is used.

The exact selected step count is **{steps}** at **{dt:.17g} s**, with {frames} frames including the initial state. All-step float64 histories forecast **{raw_history_bytes:,} numeric bytes** and **{history_folder_bytes:,} bytes** for the existing NPY layout including arrays, node step indices, and layout. The wrapper forecast plus its 2 MB guard is **{wrapper_cap_check_bytes:,} / {HISTORY_CAP:,} bytes**, so the case cap passes before any solver launch.

The whole-benchmark budget forecast is **{predicted_all_bytes:,} / {ALL_OUTPUT_BUDGET:,} bytes**, including the frozen inputs/scripts/report, same-geometry static-output estimate, histories, first-damage summary, and a 4 MB allowance for runtime metadata and logs. The native run is still pending; this preflight is only a storage check.

The case is immutable at `{CASE_PATH.relative_to(ROOT).as_posix()}`. Run once with `run_benchmark.py`; `audit_benchmark.py` independently verifies the complete result and zero-load invariants. A failed attempt is preserved and must not be retried under this authorization.
"""
    README_PATH = DATA / "README.md"
    README_PATH.write_text(md, encoding="utf-8")
    print(json.dumps({"status": report["status"], "steps": steps, "dt_s": dt,
                      "history_bytes": history_folder_bytes, "history_cap_bytes": HISTORY_CAP,
                      "predicted_all_bytes": predicted_all_bytes, "all_output_budget_bytes": ALL_OUTPUT_BUDGET,
                      "frozen_input_dir": str(INPUTS), "report": str(REPORT_PATH)}, separators=(",", ":")))


if __name__ == "__main__":
    main()
