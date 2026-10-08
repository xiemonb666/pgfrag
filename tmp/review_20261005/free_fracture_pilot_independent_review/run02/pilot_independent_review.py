#!/usr/bin/env python3
"""Independent, read-only audit for the one-shot free-fracture pilot.

This script never starts Julia or writes outside its own review directory. It
first requires a completed, final launcher receipt, then independently derives
bond failures, per-point damage/activity, active-graph components, mechanics,
and diagnostic energy histories from the frozen point tensor and reference
roster. Producer summary fields are not used as numerical oracles.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sys
import tomllib
import traceback
from collections import Counter
from array import array
from pathlib import Path
from typing import Any, Iterator

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[4]
EXPERIMENT_DIR = REPO_ROOT / "pgfrag" / "experiments" / "free_fracture_pilot_20261005"
ATTEMPT_DIR = EXPERIMENT_DIR / "run_attempt_01"
RESULTS_DIR = ATTEMPT_DIR / "results"
OUTPUT_DIR = SCRIPT_DIR / "audit_output"

EXPECTED_SOURCE_HASHES = {
    "pgfrag/experiments/free_fracture_pilot_20261005/run_free_fracture_pilot.jl":
        "9922d52c13468df5152d88c1fdb1b7c9fc6ce473a1a3a9e56507aaf37479f50b",
    "pgfrag/experiments/free_fracture_pilot_20261005/fixed_case.toml":
        "45b9f8bb5d10522fa2e6fa91ca8d944e2d6af27396034b489e27a026783bbedd",
    "pgfrag/experiments/free_fracture_pilot_20261005/run_once.ps1":
        "97434a0f5b17dd93ae7a0064f98515b4651aad583b0a039b7907c646c89ba69d",
}
EXPECTED_FREEZE_SHA256 = "1f4c3e9b10a607d303bbc3435bf22ae8b591ce8f6b5651aacb6cba3c3bfd7217"
EXPECTED_FREEZE_PIN_COUNT = 84
EXPECTED_POINTS = 4500
EXPECTED_GRID = (20, 15, 15)
EXPECTED_DT_S = 3.75e-8
EXPECTED_STEPS = 96

STATE_FIELDS = [
    "ux_m", "uy_m", "uz_m",
    "vx_m_s", "vy_m_s", "vz_m_s",
    "vhalf_x_m_s", "vhalf_y_m_s", "vhalf_z_m_s",
    "accel_x_m_s2", "accel_y_m_s2", "accel_z_m_s2",
    "b_int_x_N_m3", "b_int_y_N_m3", "b_int_z_N_m3",
    "b_ext_x_N_m3", "b_ext_y_N_m3", "b_ext_z_N_m3",
    "osb_strain_energy_density_J_m3", "damage", "n_active_bonds_exact_float64",
]
STATE_UNITS = [
    "m", "m", "m", "m/s", "m/s", "m/s", "m/s", "m/s", "m/s",
    "m/s^2", "m/s^2", "m/s^2", "N/m^3", "N/m^3", "N/m^3",
    "N/m^3", "N/m^3", "N/m^3", "J/m^3", "1", "count",
]
POINT_REF_HEADER = [
    "point_id", "x_ref_m", "y_ref_m", "z_ref_m", "volume_m3", "mass_kg",
    "complete_spherical_horizon", "minimum_face_clearance_minus_horizon_m",
    "left_grip", "right_grip", "reference_active_bonds",
]
BOND_HEADER = [
    "point_i", "point_j", "reference_length_m", "critical_stretch",
    "fail_permit", "active_at_reference",
]
EVENT_HEADER = [
    "step", "time_s", "point_i", "point_j", "reference_length_m",
    "current_length_m", "stretch", "critical_stretch", "stretch_over_critical",
    "fail_permit", "active_before", "active_after",
]
GRAPH_HEADER = [
    "step", "time_s", "directed_bond_count", "active_directed_bonds",
    "canonical_pair_count", "active_canonical_pairs", "broken_canonical_pairs",
    "new_break_events", "cumulative_break_events", "graph_component_count",
    "largest_component_mass_kg", "largest_component_mass_fraction",
    "total_material_mass_kg", "max_bond_stretch", "max_stretch_over_critical",
    "max_point_damage",
]
ENERGY_HEADER = [
    "step", "time_s", "bc_impulse_interval_start_s", "bc_impulse_interval_end_s",
    "mass_kg", "momentum_x_kg_m_s", "momentum_y_kg_m_s", "momentum_z_kg_m_s",
    "native_net_force_x_N", "native_net_force_y_N", "native_net_force_z_N",
    "osb_strain_energy_J", "kinetic_energy_J", "free_dof_kinetic_energy_J",
    "constrained_dof_raw_full_step_kinetic_energy_J", "left_grip_u_x_mean_m",
    "right_grip_u_x_mean_m", "left_bc_impulse_x_N_s", "left_bc_impulse_y_N_s",
    "left_bc_impulse_z_N_s", "left_bc_force_average_x_N", "left_bc_force_average_y_N",
    "left_bc_force_average_z_N", "right_bc_impulse_x_N_s", "right_bc_impulse_y_N_s",
    "right_bc_impulse_z_N_s", "right_bc_force_average_x_N", "right_bc_force_average_y_N",
    "right_bc_force_average_z_N", "left_bc_discrete_kinetic_work_J",
    "right_bc_discrete_kinetic_work_J", "total_bc_discrete_kinetic_work_step_J",
    "actuator_work_trapezoid_surrogate_step_J",
    "actuator_work_trapezoid_surrogate_cumulative_J",
    "max_free_dof_impulse_residual_N_s", "momentum_update_residual_x_N_s",
    "momentum_update_residual_y_N_s", "momentum_update_residual_z_N_s",
    "energy_balance_residual_diagnostic_J",
]


class AuditHalt(RuntimeError):
    pass


class AuditState:
    def __init__(self) -> None:
        self.checks: dict[str, dict[str, Any]] = {}
        self.errors: list[str] = []

    def check(self, name: str, passed: bool, evidence: Any = None) -> bool:
        value = bool(passed)
        self.checks[name] = {"pass": value, "evidence": evidence}
        if not value:
            self.errors.append(name)
        return value

    def require(self, name: str, passed: bool, evidence: Any = None) -> None:
        if not self.check(name, passed, evidence):
            raise AuditHalt(f"{name}: {evidence}")


RECEIPT_SIZE_DISCREPANCY_CHECK = "receipt_final_attempt_bytes_match"
NATIVE_ARITHMETIC_CHECKS = (
    "point_binary_exact_element_count",
    "all_point_tensor_fields_finite",
    "graph_label_binary_exact_element_count",
    "point_reference_row_count",
    "point_reference_ids_are_complete_and_ordered",
    "reference_coordinates_form_complete_analytic_lattice",
    "point_volumes_and_masses_match_uniform_cell_case",
    "reference_horizon_columns_recompute",
    "grip_masks_match_analytic_end_layers",
    "canonical_roster_nonempty_and_bounded",
    "canonical_roster_initial_state_and_permissions",
    "canonical_reference_lengths_match_point_geometry",
    "canonical_roster_degrees_match_point_reference",
    "canonical_critical_stretch_matches_fixed_osb_formula",
    "event_rows_have_valid_reciprocal_transition_fields",
    "binary_active_bond_count_encoding_is_exact_integer",
    "first_break_events_recomputed_from_deformed_coordinates",
    "active_counts_damage_and_graph_labels_recomputed_each_state",
    "graph_history_csv_matches_independent_topology",
    "96_step_physical_time_secant_grip_path",
    "captured_half_velocity_is_the_native_grip_secant",
    "raw_full_velocity_kinetic_energy_partition",
    "energy_reaction_history_recomputed_from_binary",
    "velocity_verlet_state_time_levels_and_free_dof_updates",
    "free_dof_impulse_is_zero_and_momentum_identity_recomputed",
)


def overall_review_status(audit: AuditState) -> str:
    if not audit.errors:
        return "PASS"
    if audit.errors == [RECEIPT_SIZE_DISCREPANCY_CHECK]:
        return "AUDITED_WITH_RECEIPT_SIZE_DISCREPANCY"
    return "FAIL"


def native_arithmetic_validation(audit: AuditState) -> dict[str, Any]:
    selected = {name: audit.checks[name] for name in NATIVE_ARITHMETIC_CHECKS
                if name in audit.checks}
    missing = [name for name in NATIVE_ARITHMETIC_CHECKS if name not in selected]
    passed = (not missing and all(check["pass"] for check in selected.values()))
    return {
        "status": "PASS" if passed else ("FAIL" if selected else "NOT_COMPLETED"),
        "meaning": "Numerical agreement of frozen native states, event/graph reconstructions, mechanics and energy-reaction columns only; no physical qualification.",
        "check_count": len(selected),
        "required_check_count": len(NATIVE_ARITHMETIC_CHECKS),
        "missing_checks": missing,
        "checks": selected,
    }


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_binding(path: Path) -> dict[str, Any]:
    return {"path": str(path.resolve()), "bytes": path.stat().st_size,
            "sha256": sha256_file(path)}


def require_owned_output() -> None:
    output_resolved = OUTPUT_DIR.resolve()
    own_resolved = SCRIPT_DIR.resolve()
    if output_resolved == own_resolved or own_resolved not in output_resolved.parents:
        raise AuditHalt("audit output path is not a child of the assigned review directory")


def normalize_path(path: str | Path) -> str:
    return os.path.normcase(os.path.normpath(str(Path(path).resolve())))


def verify_source_bindings(audit: AuditState) -> dict[str, Any]:
    freeze_path = EXPERIMENT_DIR / "source_freeze.toml"
    freeze_raw = freeze_path.read_bytes()
    freeze_sha = hashlib.sha256(freeze_raw).hexdigest()
    audit.require("source_freeze_sha256_is_the_reviewed_freeze",
                  freeze_sha == EXPECTED_FREEZE_SHA256,
                  {"expected": EXPECTED_FREEZE_SHA256, "actual": freeze_sha})
    freeze = tomllib.loads(freeze_raw.decode("utf-8"))
    pins = freeze.get("pins", [])
    audit.require("frozen_source_pin_count",
                  freeze.get("pin_count") == EXPECTED_FREEZE_PIN_COUNT == len(pins),
                  {"declared": freeze.get("pin_count"), "observed": len(pins)})
    pin_paths = [str(pin["path"]) for pin in pins]
    audit.require("frozen_pin_paths_unique", len(pin_paths) == len(set(pin_paths)),
                  {"pin_count": len(pin_paths), "unique": len(set(pin_paths))})
    pin_results: list[dict[str, Any]] = []
    for pin in pins:
        relative = str(pin["path"])
        candidate = (REPO_ROOT / Path(relative.replace("/", os.sep))).resolve(strict=True)
        if REPO_ROOT.resolve() not in candidate.parents:
            raise AuditHalt(f"source-freeze path escaped repository root: {relative}")
        observed = file_binding(candidate)
        expected_bytes = int(pin["bytes"])
        expected_sha = str(pin["sha256"]).lower()
        passed = observed["bytes"] == expected_bytes and observed["sha256"] == expected_sha
        pin_results.append({"relative_path": relative, "expected_bytes": expected_bytes,
                            "expected_sha256": expected_sha,
                            "observed_bytes": observed["bytes"],
                            "observed_sha256": observed["sha256"], "pass": passed})
        if not passed:
            raise AuditHalt(f"frozen source pin changed: {relative}")
    audit.require("all_frozen_source_pins_match", all(p["pass"] for p in pin_results),
                  {"verified": len(pin_results), "failed": 0})
    named_bindings: dict[str, Any] = {}
    for relative, expected_sha in EXPECTED_SOURCE_HASHES.items():
        candidate = (REPO_ROOT / Path(relative.replace("/", os.sep))).resolve(strict=True)
        observed = file_binding(candidate)
        pin = next((p for p in pins if str(p["path"]) == relative), None)
        passed = (observed["sha256"] == expected_sha and pin is not None
                  and observed["sha256"] == str(pin["sha256"]).lower()
                  and observed["bytes"] == int(pin["bytes"]))
        audit.require(f"explicit_source_pin_{Path(relative).name}", passed,
                      {"observed": observed, "expected_sha256": expected_sha})
        named_bindings[relative] = observed
    return {"freeze": {"path": str(freeze_path.resolve()), "bytes": len(freeze_raw),
                       "sha256": freeze_sha, "schema": freeze.get("schema"),
                       "status": freeze.get("freeze_status")},
            "pins": pin_results, "named_files": named_bindings,
            "parsed_freeze": freeze}


def parse_bool(text: str) -> bool:
    value = text.strip().lower()
    if value == "true":
        return True
    if value == "false":
        return False
    raise ValueError(f"not a CSV boolean: {text!r}")


def read_csv_records(path: Path, required_header: list[str]) -> list[dict[str, str]]:
    return list(iter_csv_records(path, required_header))


def iter_csv_records(path: Path, required_header: list[str]) -> Iterator[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != required_header:
            raise AuditHalt(f"CSV header mismatch in {path.name}: {reader.fieldnames}")
        yield from reader


def verify_receipt_before_results(audit: AuditState,
                                  expected_final_summary_sha256: str) -> dict[str, Any]:
    receipt_path = ATTEMPT_DIR / "launcher_receipt.json"
    if not receipt_path.is_file():
        raise AuditHalt("final launcher receipt does not exist; refusing to inspect partial results")
    receipt_raw = receipt_path.read_bytes()
    receipt = json.loads(receipt_raw.decode("utf-8"))
    binding = {"path": str(receipt_path.resolve()), "bytes": len(receipt_raw),
               "sha256": hashlib.sha256(receipt_raw).hexdigest(),
               "schema": receipt.get("schema"), "status": receipt.get("status")}
    audit.require("receipt_schema", receipt.get("schema") == "pgfrag.free_fracture_pilot.launch_receipt.v1",
                  binding)
    fields = {
        "completed_exit_zero": receipt.get("status") == "COMPLETED_EXIT_0"
            and receipt.get("native_process_exit_code") == 0,
        "process_started_and_not_killed": receipt.get("process_started") is True
            and receipt.get("native_process_killed") is False,
        "no_timeout_or_cap_exceeded": receipt.get("timed_out") is False
            and receipt.get("disk_cap_exceeded") is False
            and receipt.get("output_cap_exceeded") is False
            and receipt.get("cap_exceeded") is False,
        "final_budget_authority": receipt.get("final_budget_ok") is True
            and receipt.get("result_authoritative") is True,
        "summary_status_provenance": receipt.get("summary_status") == "ENGINEERING_DIAGNOSTIC_COMPLETED",
        "output_path_matches_fixed_attempt": normalize_path(receipt.get("output_directory", ""))
            == normalize_path(RESULTS_DIR),
    }
    for key, passed in fields.items():
        audit.require(f"receipt_{key}", passed, {"value": receipt.get(key),
                      "output_directory": receipt.get("output_directory")})
    observed_summary_sha = str(receipt.get("summary_sha256", "")).lower()
    audit.require("authorized_final_summary_hash_matches_receipt",
                  observed_summary_sha == expected_final_summary_sha256.lower(),
                  {"authorized_expected_sha256": expected_final_summary_sha256.lower(),
                   "receipt_summary_sha256": observed_summary_sha})
    return {"path": receipt_path, "receipt": receipt, "binding": binding,
            "checks": fields}


def verify_receipt_provenance(audit: AuditState, receipt_info: dict[str, Any],
                              source_info: dict[str, Any], result_files: dict[str, Any],
                              result_bytes: int, max_attempt_bytes: int,
                              max_result_bytes: int) -> dict[str, Any]:
    receipt = receipt_info["receipt"]
    required_files = {
        "driver_sha256": "pgfrag/experiments/free_fracture_pilot_20261005/run_free_fracture_pilot.jl",
        "fixed_case_sha256": "pgfrag/experiments/free_fracture_pilot_20261005/fixed_case.toml",
        "source_freeze_sha256": "pgfrag/experiments/free_fracture_pilot_20261005/source_freeze.toml",
    }
    comparisons: dict[str, Any] = {}
    for receipt_key, relative in required_files.items():
        if receipt_key == "source_freeze_sha256":
            observed = source_info["freeze"]["sha256"]
        else:
            observed = source_info["named_files"][relative]["sha256"]
        claimed = str(receipt.get(receipt_key, "")).lower()
        passed = claimed == observed
        audit.require(f"receipt_{receipt_key}_binds", passed,
                      {"receipt": claimed, "observed": observed})
        comparisons[receipt_key] = {"receipt": claimed, "observed": observed, "pass": passed}
    pins_text = str(receipt.get("required_source_pins", ""))
    receipt_pins: dict[str, str] = {}
    for item in pins_text.split(";"):
        if not item:
            continue
        if "=" not in item:
            raise AuditHalt("malformed required_source_pins receipt field")
        name, value = item.split("=", 1)
        if name in receipt_pins:
            raise AuditHalt(f"duplicate receipt source pin: {name}")
        receipt_pins[name] = value.lower()
    expected_receipt_paths = {
        **EXPECTED_SOURCE_HASHES,
        "third_party/pgfrag-runtime/extract/julia-1.10.12/bin/julia.exe": None,
    }
    freeze_pins = {str(p["path"]): str(p["sha256"]).lower()
                   for p in source_info["parsed_freeze"].get("pins", [])}
    required_pin_checks = {}
    for relative, expected in expected_receipt_paths.items():
        if expected is None:
            expected = freeze_pins.get(relative)
        present = relative in receipt_pins
        passed = present and expected is not None and receipt_pins.get(relative) == expected
        required_pin_checks[relative] = {"expected": expected,
                                         "receipt": receipt_pins.get(relative),
                                         "pass": passed}
        audit.require("receipt_pin_" + Path(relative).name, passed,
                      required_pin_checks[relative])
    attempt_bytes = sum(p.stat().st_size for p in ATTEMPT_DIR.rglob("*") if p.is_file())
    audit.require("receipt_experiment_output_bytes_match",
                  int(receipt.get("experiment_output_bytes", -1)) == result_bytes,
                  {"receipt": receipt.get("experiment_output_bytes"), "observed": result_bytes})
    recorded_attempt_bytes = int(receipt.get("final_attempt_bytes", -1))
    audit.check(RECEIPT_SIZE_DISCREPANCY_CHECK,
                recorded_attempt_bytes == attempt_bytes,
                {"receipt_recorded_bytes": recorded_attempt_bytes,
                 "independently_observed_final_attempt_tree_bytes": attempt_bytes,
                 "signed_delta_observed_minus_receipt": attempt_bytes - recorded_attempt_bytes,
                 "tolerance_bytes": 0})
    frozen_output_cap = int(source_info["parsed_freeze"]["max_experiment_output_bytes"])
    audit.require("actual_final_attempt_tree_within_case_attempt_cap",
                  attempt_bytes <= max_attempt_bytes,
                  {"actual_final_attempt_tree_bytes": attempt_bytes,
                   "case_max_attempt_bytes": max_attempt_bytes})
    audit.require("actual_result_tree_within_case_output_cap",
                  result_bytes <= max_result_bytes,
                  {"actual_result_tree_bytes": result_bytes,
                   "case_max_output_bytes": max_result_bytes})
    audit.require("actual_result_tree_within_frozen_output_cap",
                  result_bytes <= frozen_output_cap,
                  {"actual_result_tree_bytes": result_bytes,
                   "frozen_max_experiment_output_bytes": frozen_output_cap})
    audit.require("receipt_summary_sha256_binds",
                  receipt.get("summary_sha256") == result_files["summary.toml"]["sha256"],
                  {"receipt": receipt.get("summary_sha256"),
                   "observed": result_files["summary.toml"]["sha256"]})
    log_checks = {}
    for prefix, filename in (("stdout", "native_stdout.log"), ("stderr", "native_stderr.log")):
        path = ATTEMPT_DIR / filename
        observed = file_binding(path)
        expected_path = normalize_path(ATTEMPT_DIR / filename)
        claimed_path = normalize_path(receipt.get(f"{prefix}_path", ""))
        passed = (observed["bytes"] == int(receipt.get(f"{prefix}_bytes", -1))
                  and observed["sha256"] == str(receipt.get(f"{prefix}_sha256", "")).lower()
                  and expected_path == claimed_path)
        audit.require(f"receipt_{prefix}_log_binds", passed,
                      {**observed, "receipt_bytes": receipt.get(f"{prefix}_bytes"),
                       "receipt_sha256": receipt.get(f"{prefix}_sha256"),
                       "receipt_path": receipt.get(f"{prefix}_path")})
        log_checks[prefix] = {**observed, "pass": passed}
    return {"source_hash_comparisons": comparisons,
            "required_source_pins": required_pin_checks,
            "logs": log_checks, "experiment_output_bytes": result_bytes,
            "attempt_bytes": attempt_bytes,
            "case_max_attempt_bytes": max_attempt_bytes,
            "case_max_output_bytes": max_result_bytes,
            "frozen_max_experiment_output_bytes": frozen_output_cap}


def verify_result_inventory(audit: AuditState) -> tuple[dict[str, Any], int]:
    expected_names = {
        "point_state_history.f64le", "point_state_layout.toml",
        "graph_component_labels.i32le", "graph_component_labels_layout.toml",
        "point_reference.csv", "canonical_bonds.csv", "native_failure_events.csv",
        "graph_state_history.csv", "energy_reaction_history.csv", "summary.toml",
    }
    observed_names = {p.name for p in RESULTS_DIR.iterdir() if p.is_file()}
    all_nested = [p for p in RESULTS_DIR.rglob("*") if p.is_file()]
    nested_ok = all(p.parent == RESULTS_DIR for p in all_nested)
    audit.require("result_inventory_exact_and_flat", observed_names == expected_names and nested_ok,
                  {"expected": sorted(expected_names), "observed": sorted(observed_names),
                   "nested_files": [str(p) for p in all_nested if p.parent != RESULTS_DIR]})
    files = {name: file_binding(RESULTS_DIR / name) for name in sorted(expected_names)}
    total = sum(info["bytes"] for info in files.values())
    audit.require("result_directory_under_case_output_cap",
                  total <= int(tomllib.loads((EXPERIMENT_DIR / "fixed_case.toml").read_text("utf-8"))
                               ["limits"]["max_output_bytes"]),
                  {"bytes": total})
    return files, total


def verify_binary_sidecar(audit: AuditState, binary_name: str, sidecar_name: str,
                          expected_schema: str) -> tuple[dict[str, Any], dict[str, Any]]:
    binary_path = RESULTS_DIR / binary_name
    sidecar_path = RESULTS_DIR / sidecar_name
    sidecar_raw = sidecar_path.read_bytes()
    sidecar = tomllib.loads(sidecar_raw.decode("utf-8"))
    binary_info = file_binding(binary_path)
    metadata = {
        "sidecar": {"path": str(sidecar_path.resolve()), "bytes": len(sidecar_raw),
                    "sha256": hashlib.sha256(sidecar_raw).hexdigest(),
                    "schema": sidecar.get("schema")},
        "binary": binary_info,
    }
    audit.require(f"{binary_name}_sidecar_schema", sidecar.get("schema") == expected_schema,
                  metadata["sidecar"])
    audit.require(f"{binary_name}_byte_and_digest_binding",
                  binary_info["bytes"] == int(sidecar.get("byte_count", -1))
                  and binary_info["sha256"] == str(sidecar.get("sha256", "")).lower(),
                  {**binary_info, "sidecar_bytes": sidecar.get("byte_count"),
                   "sidecar_sha256": sidecar.get("sha256")})
    audit.require(f"{binary_name}_declared_data_file", sidecar.get("data_file") == binary_name,
                  {"declared": sidecar.get("data_file"), "observed": binary_name})
    return sidecar, metadata


def to_float(row: dict[str, str], key: str) -> float:
    value = float(row[key])
    if not math.isfinite(value):
        raise AuditHalt(f"non-finite CSV value {key}: {row[key]}")
    return value


def to_int(row: dict[str, str], key: str) -> int:
    value = int(row[key])
    return value


def close_scalar(a: float, b: float, *, rtol: float = 2e-10,
                 atol: float = 1e-18) -> bool:
    return abs(a - b) <= atol + rtol * max(abs(a), abs(b))


def target_displacement(t: float, ramp: float, hold: float,
                        unload: float, amplitude: float) -> float:
    if t < 0:
        return 0.0
    if t <= ramp:
        x = t / ramp
        return amplitude * (3.0 * x * x - 2.0 * x * x * x)
    if t <= ramp + hold:
        return amplitude
    if t <= ramp + hold + unload:
        x = (t - ramp - hold) / unload
        return amplitude * (1.0 - 3.0 * x * x + 2.0 * x * x * x)
    return 0.0


def graph_components(active: np.ndarray, pairs: np.ndarray, masses: np.ndarray,
                     coords_m: np.ndarray, left_mask: np.ndarray,
                     right_mask: np.ndarray) -> tuple[np.ndarray, list[dict[str, Any]]]:
    n_points = len(masses)
    parent = np.arange(n_points + 1, dtype=np.int32)
    rank = np.zeros(n_points + 1, dtype=np.uint8)

    def find(a: int) -> int:
        root = a
        while parent[root] != root:
            parent[root] = parent[parent[root]]
            root = int(parent[root])
        return root

    for edge in np.flatnonzero(active):
        a = find(int(pairs[edge, 0]))
        b = find(int(pairs[edge, 1]))
        if a == b:
            continue
        if rank[a] < rank[b]:
            parent[a] = b
        elif rank[a] > rank[b]:
            parent[b] = a
        else:
            parent[b] = a
            rank[a] += 1

    root_min: dict[int, int] = {}
    roots = np.empty(n_points, dtype=np.int32)
    for point in range(1, n_points + 1):
        root = find(point)
        roots[point - 1] = root
        root_min[root] = min(root_min.get(root, point), point)
    labels = np.fromiter((root_min[int(root)] for root in roots),
                         dtype=np.int32, count=n_points)
    unique, inverse, point_counts = np.unique(labels, return_inverse=True,
                                               return_counts=True)
    comp_mass = np.bincount(inverse, weights=masses, minlength=len(unique))
    left_counts = np.bincount(inverse, weights=left_mask.astype(np.float64),
                              minlength=len(unique)).astype(np.int64)
    right_counts = np.bincount(inverse, weights=right_mask.astype(np.float64),
                               minlength=len(unique)).astype(np.int64)
    xyz_mm = coords_m * 1000.0
    bbox_min = np.full((len(unique), 3), np.inf, dtype=np.float64)
    bbox_max = np.full((len(unique), 3), -np.inf, dtype=np.float64)
    np.minimum.at(bbox_min, inverse, xyz_mm)
    np.maximum.at(bbox_max, inverse, xyz_mm)
    comp_records: list[dict[str, Any]] = []
    for idx, label in enumerate(unique):
        comp_records.append({
            "label_min_point_id": int(label),
            "point_count": int(point_counts[idx]),
            "mass_kg": float(comp_mass[idx]),
            "mass_fraction": float(comp_mass[idx] / np.sum(masses)),
            "left_grip_point_count": int(left_counts[idx]),
            "right_grip_point_count": int(right_counts[idx]),
            "bbox_mm": [[float(bbox_min[idx, dim]), float(bbox_max[idx, dim])]
                        for dim in range(3)],
        })
    comp_records.sort(key=lambda r: (-r["mass_kg"], r["label_min_point_id"]))
    return labels, comp_records


def component_bboxes(labels: np.ndarray, coords_m: np.ndarray) -> dict[int, list[list[float]]]:
    """Recompute only component bounds for a selected state with unchanged topology."""
    unique, inverse = np.unique(labels, return_inverse=True)
    xyz_mm = coords_m * 1000.0
    bbox_min = np.full((len(unique), 3), np.inf, dtype=np.float64)
    bbox_max = np.full((len(unique), 3), -np.inf, dtype=np.float64)
    np.minimum.at(bbox_min, inverse, xyz_mm)
    np.maximum.at(bbox_max, inverse, xyz_mm)
    return {
        int(label): [[float(bbox_min[idx, dim]), float(bbox_max[idx, dim])]
                     for dim in range(3)]
        for idx, label in enumerate(unique)
    }


def _expected_source_steps(case: dict[str, Any]) -> dict[str, Any]:
    loading = case["loading"]
    solver = case["solver"]
    durations = [float(loading[name]) for name in (
        "ramp_duration_s", "hold_duration_s", "unload_duration_s",
        "post_unload_dwell_duration_s")]
    dt = float(solver["expected_time_step_s"])
    total = sum(durations)
    return {
        "dt_s": dt,
        "total_duration_s": total,
        "n_steps": math.ceil(total / dt),
        "phase_intervals_ceiling": [math.ceil(d / dt) for d in durations],
        "phase_durations_s": durations,
    }


def audit_pilot(audit: AuditState, receipt_info: dict[str, Any]) -> dict[str, Any]:
    source_info = verify_source_bindings(audit)
    freeze = source_info["parsed_freeze"]
    case_path = EXPERIMENT_DIR / "fixed_case.toml"
    case_raw = case_path.read_bytes()
    case = tomllib.loads(case_raw.decode("utf-8"))
    audit.require("fixed_case_schema_and_source_only_status",
                  case.get("schema") == "pgfrag.free_fracture_pilot.case.v1"
                  and case.get("run_status") == "SOURCE_ONLY_NOT_RUN",
                  {"schema": case.get("schema"), "run_status": case.get("run_status")})
    audit.require("case_id_matches_frozen_source",
                  case.get("case_id") == freeze.get("case_id"),
                  {"case": case.get("case_id"), "freeze": freeze.get("case_id")})
    source_binding = file_binding(case_path)
    audit.require("fixed_case_sha_matches_freeze_and_receipt",
                  source_binding["sha256"] == EXPECTED_SOURCE_HASHES[
                      "pgfrag/experiments/free_fracture_pilot_20261005/fixed_case.toml"]
                  and source_binding["sha256"] == source_info["named_files"][
                      "pgfrag/experiments/free_fracture_pilot_20261005/fixed_case.toml"]["sha256"],
                  source_binding)
    result_files, result_bytes = verify_result_inventory(audit)
    provenance = verify_receipt_provenance(
        audit, receipt_info, source_info, result_files, result_bytes,
        int(case["limits"]["max_attempt_bytes"]),
        int(case["limits"]["max_output_bytes"]))

    point_sidecar, point_binding = verify_binary_sidecar(
        audit, "point_state_history.f64le", "point_state_layout.toml",
        "pgfrag.free_fracture_pilot.point_state_binary.v1")
    graph_sidecar, graph_binding = verify_binary_sidecar(
        audit, "graph_component_labels.i32le", "graph_component_labels_layout.toml",
        "pgfrag.free_fracture_pilot.graph_labels_binary.v1")
    geometry = case["geometry"]
    material = case["material"]
    loading = case["loading"]
    solver = case["solver"]
    limits = case["limits"]
    dt = float(solver["expected_time_step_s"])
    n_points = int(geometry["expected_point_count"])
    steps = _expected_source_steps(case)
    n_steps = int(steps["n_steps"])
    n_states = n_steps + 1
    audit.require("bounded_case_uses_expected_point_and_step_counts",
                  n_points == EXPECTED_POINTS and tuple(geometry["expected_grid_shape"]) == EXPECTED_GRID
                  and n_steps == EXPECTED_STEPS and n_steps <= int(limits["max_steps"]),
                  {"points": n_points, "grid": geometry["expected_grid_shape"],
                   "steps": n_steps, "states": n_states, "phase_interval_ceilings": steps["phase_intervals_ceiling"]})
    audit.require("frozen_dt_is_exact_expected_value",
                  dt == EXPECTED_DT_S and dt == float(point_sidecar.get("time_step_s", -1.0))
                  and dt == float(graph_sidecar.get("time_step_s", -1.0)),
                  {"case_dt": dt, "point_layout_dt": point_sidecar.get("time_step_s"),
                   "graph_layout_dt": graph_sidecar.get("time_step_s")})
    expected_shape = [21, n_points, n_states]
    audit.require("point_binary_layout_is_frozen_shape_and_names",
                  point_sidecar.get("array_shape") == expected_shape
                  and point_sidecar.get("field_names") == STATE_FIELDS
                  and point_sidecar.get("field_units") == STATE_UNITS
                  and point_sidecar.get("storage_order") ==
                      "Julia column-major tensor [field_index, point_id, state_index]; field index varies fastest, then ascending point_id, then ascending state_index"
                  and point_sidecar.get("format") == "IEEE-754 Float64, little-endian",
                  {"shape": point_sidecar.get("array_shape"), "field_count": point_sidecar.get("field_count"),
                   "state_count": point_sidecar.get("state_count")})
    audit.require("graph_binary_layout_is_frozen_shape_and_labels",
                  graph_sidecar.get("array_shape") == [n_points, n_states]
                  and graph_sidecar.get("format") == "signed Int32, little-endian"
                  and graph_sidecar.get("label_units") == "point ID; zero is invalid and must not occur",
                  {"shape": graph_sidecar.get("array_shape"), "format": graph_sidecar.get("format")})

    point_path = RESULTS_DIR / "point_state_history.f64le"
    raw_state = np.fromfile(point_path, dtype="<f8")
    expected_items = 21 * n_points * n_states
    audit.require("point_binary_exact_element_count",
                  raw_state.size == expected_items and point_path.stat().st_size == expected_items * 8,
                  {"elements": int(raw_state.size), "expected": expected_items,
                   "bytes": point_path.stat().st_size})
    tensor = raw_state.reshape((n_states, n_points, 21)).transpose(2, 1, 0)
    audit.require("all_point_tensor_fields_finite", bool(np.isfinite(tensor).all()),
                  {"shape": list(tensor.shape)})
    label_path = RESULTS_DIR / "graph_component_labels.i32le"
    raw_labels = np.fromfile(label_path, dtype="<i4")
    audit.require("graph_label_binary_exact_element_count",
                  raw_labels.size == n_points * n_states,
                  {"elements": int(raw_labels.size), "expected": n_points * n_states})
    labels_binary = raw_labels.reshape((n_states, n_points)).T

    point_rows = read_csv_records(RESULTS_DIR / "point_reference.csv", POINT_REF_HEADER)
    audit.require("point_reference_row_count", len(point_rows) == n_points,
                  {"rows": len(point_rows), "expected": n_points})
    point_ids = np.array([to_int(r, "point_id") for r in point_rows], dtype=np.int32)
    xyz = np.array([[to_float(r, key) for key in ("x_ref_m", "y_ref_m", "z_ref_m")]
                    for r in point_rows], dtype=np.float64)
    volume = np.array([to_float(r, "volume_m3") for r in point_rows], dtype=np.float64)
    mass = np.array([to_float(r, "mass_kg") for r in point_rows], dtype=np.float64)
    point_complete = np.array([parse_bool(r["complete_spherical_horizon"]) for r in point_rows])
    horizon_margin = np.array([to_float(r, "minimum_face_clearance_minus_horizon_m")
                               for r in point_rows], dtype=np.float64)
    left_mask = np.array([parse_bool(r["left_grip"]) for r in point_rows], dtype=bool)
    right_mask = np.array([parse_bool(r["right_grip"]) for r in point_rows], dtype=bool)
    ref_degrees_csv = np.array([to_int(r, "reference_active_bonds") for r in point_rows],
                               dtype=np.int32)
    audit.require("point_reference_ids_are_complete_and_ordered",
                  np.array_equal(point_ids, np.arange(1, n_points + 1, dtype=np.int32)),
                  {"first": int(point_ids[0]), "last": int(point_ids[-1]),
                   "unique": int(np.unique(point_ids).size)})
    lengths_m = np.array([float(geometry[k]) for k in
                          ("length_x_m", "length_y_m", "length_z_m")], dtype=np.float64)
    h = float(geometry["point_spacing_h_m"])
    delta_horizon = float(geometry["horizon_delta_m"])
    expected_axes = [(np.arange(n, dtype=np.float64) + 0.5) * h - length / 2
                     for n, length in zip(EXPECTED_GRID, lengths_m)]
    grid_indexes = []
    axis_errors = []
    for dim, axis in enumerate(expected_axes):
        nearest = np.abs(xyz[:, dim, None] - axis[None, :]).argmin(axis=1)
        grid_indexes.append(nearest)
        axis_errors.append(float(np.max(np.abs(xyz[:, dim] - axis[nearest]))))
    combos = np.stack(grid_indexes, axis=1)
    audit.require("reference_coordinates_form_complete_analytic_lattice",
                  xyz.shape == (n_points, 3) and len(np.unique(combos, axis=0)) == n_points
                  and max(axis_errors) <= max(1e-14, h * 1e-9),
                  {"unique_lattice_points": int(len(np.unique(combos, axis=0))),
                   "axis_max_error_m": axis_errors, "expected_grid": list(EXPECTED_GRID)})
    audit.require("point_volumes_and_masses_match_uniform_cell_case",
                  bool(np.allclose(volume, h ** 3, rtol=1e-12, atol=0.0))
                  and bool(np.allclose(mass, float(material["rho_kg_m3"]) * volume,
                                       rtol=1e-12, atol=0.0)),
                  {"volume_m3": float(volume[0]), "mass_kg": float(mass[0]),
                   "total_mass_kg": float(mass.sum()),
                   "expected_total_mass_kg": float(material["rho_kg_m3"] * np.prod(lengths_m))})
    clearance = np.min(lengths_m[None, :] / 2 - np.abs(xyz), axis=1)
    recomputed_margin = clearance - delta_horizon
    expected_complete = clearance + 64.0 * np.finfo(np.float64).eps * delta_horizon >= delta_horizon
    audit.require("reference_horizon_columns_recompute",
                  np.array_equal(point_complete, expected_complete)
                  and np.allclose(horizon_margin, recomputed_margin, rtol=1e-10, atol=1e-15),
                  {"complete_points": int(point_complete.sum()),
                   "expected_complete_points": int(expected_complete.sum()),
                   "max_margin_delta_m": float(np.max(np.abs(horizon_margin - recomputed_margin)))})
    xmin, xmax = float(xyz[:, 0].min()), float(xyz[:, 0].max())
    grip_tolerance = h * 1e-8
    expected_left = np.abs(xyz[:, 0] - xmin) <= grip_tolerance
    expected_right = np.abs(xyz[:, 0] - xmax) <= grip_tolerance
    audit.require("grip_masks_match_analytic_end_layers",
                  np.array_equal(left_mask, expected_left)
                  and np.array_equal(right_mask, expected_right)
                  and not np.any(left_mask & right_mask),
                  {"left_count": int(left_mask.sum()), "right_count": int(right_mask.sum()),
                   "overlap": int(np.count_nonzero(left_mask & right_mask))})

    bond_i = array("i")
    bond_j = array("i")
    bond_ref = array("d")
    bond_critical = array("d")
    bond_permit_bytes = bytearray()
    bond_active_bytes = bytearray()
    n_bonds = 0
    audit.require("canonical_roster_nonempty_and_bounded",
                  array("i").itemsize == 4 and array("d").itemsize == 8,
                  {"native_int_bytes": array("i").itemsize,
                   "native_double_bytes": array("d").itemsize})
    pair_index: dict[tuple[int, int], int] = {}
    previous_pair = (0, 0)
    for row in iter_csv_records(RESULTS_DIR / "canonical_bonds.csv", BOND_HEADER):
        k = n_bonds
        i = to_int(row, "point_i")
        j = to_int(row, "point_j")
        pair = (i, j)
        if not (1 <= i < j <= n_points) or pair <= previous_pair:
            raise AuditHalt(f"canonical roster endpoints/order invalid at row {k + 2}: {pair}")
        previous_pair = pair
        bond_i.append(i)
        bond_j.append(j)
        bond_ref.append(to_float(row, "reference_length_m"))
        bond_critical.append(to_float(row, "critical_stretch"))
        bond_permit_bytes.append(int(parse_bool(row["fail_permit"])))
        bond_active_bytes.append(int(parse_bool(row["active_at_reference"])))
        pair_index[pair] = k
        n_bonds += 1
    audit.require("canonical_roster_nonempty_and_bounded",
                  0 < n_bonds <= int(freeze["max_experiment_output_bytes"])
                  and n_bonds <= 2_000_000,
                  {"canonical_bonds": n_bonds})
    pairs = np.column_stack((np.frombuffer(bond_i, dtype="<i4"),
                             np.frombuffer(bond_j, dtype="<i4"))).copy()
    ref_length = np.frombuffer(bond_ref, dtype="<f8").copy()
    critical = np.frombuffer(bond_critical, dtype="<f8").copy()
    fail_permit = np.frombuffer(bond_permit_bytes, dtype=np.uint8).astype(bool)
    active_ref_flag = np.frombuffer(bond_active_bytes, dtype=np.uint8).astype(bool)
    audit.require("canonical_roster_initial_state_and_permissions",
                  bool(np.all(active_ref_flag)) and bool(np.all(fail_permit))
                  and bool(np.all(ref_length > 0.0)) and bool(np.all(critical > 0.0)),
                  {"pairs": n_bonds, "all_active_at_reference": bool(np.all(active_ref_flag)),
                   "all_fail_permitted": bool(np.all(fail_permit))})
    reference_pair_length = np.linalg.norm(xyz[pairs[:, 0] - 1] - xyz[pairs[:, 1] - 1], axis=1)
    max_ref_length_delta = float(np.max(np.abs(reference_pair_length - ref_length)))
    audit.require("canonical_reference_lengths_match_point_geometry",
                  bool(np.allclose(reference_pair_length, ref_length, rtol=1e-10, atol=1e-15)),
                  {"max_abs_delta_m": max_ref_length_delta})
    reference_degrees = np.bincount(pairs.ravel(), minlength=n_points + 1)[1:].astype(np.int32)
    audit.require("canonical_roster_degrees_match_point_reference",
                  np.array_equal(reference_degrees, ref_degrees_csv),
                  {"roster_degree_sum": int(reference_degrees.sum()),
                   "csv_degree_sum": int(ref_degrees_csv.sum()),
                   "max_degree_delta": int(np.max(np.abs(reference_degrees - ref_degrees_csv)))})
    rho = float(material["rho_kg_m3"])
    young = float(material["E_pa"])
    nu = float(material["nu"])
    gc = float(material["Gc_j_m2"])
    bulk = young / (3.0 * (1.0 - 2.0 * nu))
    formula_critical = math.sqrt(5.0 * gc / (9.0 * bulk * delta_horizon))
    audit.require("canonical_critical_stretch_matches_fixed_osb_formula",
                  bool(np.allclose(critical, formula_critical, rtol=1e-10, atol=1e-15)),
                  {"roster_first": float(critical[0]), "formula": formula_critical,
                   "max_abs_delta": float(np.max(np.abs(critical - formula_critical)))})

    events_rows = read_csv_records(RESULTS_DIR / "native_failure_events.csv", EVENT_HEADER)
    events_by_step: list[list[dict[str, Any]]] = [[] for _ in range(n_states)]
    event_step_by_edge = np.full(n_bonds, -1, dtype=np.int32)
    event_records: list[dict[str, Any]] = []
    last_sort_key = (-1, -1, -1)
    for row_number, row in enumerate(events_rows, start=2):
        step = to_int(row, "step")
        i, j = to_int(row, "point_i"), to_int(row, "point_j")
        pair = (i, j)
        edge = pair_index.get(pair)
        if edge is None:
            raise AuditHalt(f"failure event row {row_number} is absent from canonical roster: {pair}")
        event = {
            "step": step, "time_s": to_float(row, "time_s"), "point_i": i, "point_j": j,
            "edge_index": edge,
            "reference_length_m": to_float(row, "reference_length_m"),
            "current_length_m": to_float(row, "current_length_m"),
            "stretch": to_float(row, "stretch"),
            "critical_stretch": to_float(row, "critical_stretch"),
            "stretch_over_critical": to_float(row, "stretch_over_critical"),
            "fail_permit": parse_bool(row["fail_permit"]),
            "active_before": parse_bool(row["active_before"]),
            "active_after": parse_bool(row["active_after"]),
        }
        sort_key = (step, i, j)
        if sort_key <= last_sort_key:
            raise AuditHalt(f"failure events are not unique/sorted at row {row_number}: {sort_key}")
        last_sort_key = sort_key
        if not (1 <= step < n_states) or event_step_by_edge[edge] != -1:
            raise AuditHalt(f"failure event step or unique-first-event condition invalid at row {row_number}")
        event_step_by_edge[edge] = step
        events_by_step[step].append(event)
        event_records.append(event)
    audit.require("event_rows_have_valid_reciprocal_transition_fields",
                  all(e["fail_permit"] and e["active_before"] and not e["active_after"]
                      for e in event_records),
                  {"event_count": len(event_records)})

    graph_rows = read_csv_records(RESULTS_DIR / "graph_state_history.csv", GRAPH_HEADER)
    energy_rows = read_csv_records(RESULTS_DIR / "energy_reaction_history.csv", ENERGY_HEADER)
    audit.require("graph_and_energy_csv_state_counts",
                  len(graph_rows) == n_states and len(energy_rows) == n_states,
                  {"graph_rows": len(graph_rows), "energy_rows": len(energy_rows),
                   "expected_states": n_states})

    n_active_tensor = tensor[20]
    integer_active_ok = bool(np.isfinite(n_active_tensor).all()
                             and np.equal(n_active_tensor, np.rint(n_active_tensor)).all()
                             and np.all(n_active_tensor >= 0.0))
    audit.require("binary_active_bond_count_encoding_is_exact_integer",
                  integer_active_ok,
                  {"min": float(n_active_tensor.min()), "max": float(n_active_tensor.max())})

    active = np.ones(n_bonds, dtype=bool)
    cumulative = 0
    active_pair_counts = np.zeros(n_states, dtype=np.int64)
    broken_pair_counts = np.zeros(n_states, dtype=np.int64)
    new_event_counts = np.zeros(n_states, dtype=np.int64)
    graph_component_counts = np.zeros(n_states, dtype=np.int64)
    graph_largest_mass = np.zeros(n_states, dtype=np.float64)
    graph_largest_fraction = np.zeros(n_states, dtype=np.float64)
    graph_labels_calc = np.zeros((n_points, n_states), dtype=np.int32)
    graph_series: list[dict[str, Any]] = []
    selected_component_states: list[dict[str, Any]] = []
    previous_right_labels: tuple[int, ...] | None = None
    point_damage = tensor[19]
    actual_coordinates = xyz[None, :, :] + tensor[:3].transpose(2, 1, 0)
    max_stretch_state = np.zeros(n_states, dtype=np.float64)
    max_ratio_state = np.zeros(n_states, dtype=np.float64)
    threshold_ambiguous_count = 0
    first_break_witnesses: list[dict[str, Any]] = []
    geometry_active_count_mismatches = 0
    edge_i0 = pairs[:, 0] - 1
    edge_j0 = pairs[:, 1] - 1
    stretch_roundoff_tol = 128.0 * np.finfo(np.float64).eps
    cached_graph_labels: np.ndarray | None = None
    cached_components: list[dict[str, Any]] | None = None

    for state in range(n_states):
        step = state
        state_events = events_by_step[state]
        new_event_counts[state] = len(state_events)
        event_edges = {int(e["edge_index"]) for e in state_events}
        if state == 0 and event_edges:
            raise AuditHalt("native event CSV contains an event in reference state zero")
        if state > 0:
            pos = actual_coordinates[state]
            differences = pos[edge_j0] - pos[edge_i0]
            current_length = np.sqrt(np.einsum("ij,ij->i", differences, differences))
            if not np.isfinite(current_length).all() or np.any(current_length <= 0.0):
                raise AuditHalt(f"non-finite/nonpositive deformed canonical length at step {step}")
            stretches = (current_length - ref_length) / ref_length
            ratios = stretches / critical
            max_stretch_state[state] = float(np.max(stretches))
            max_ratio_state[state] = float(np.max(ratios))
            active_before = active.copy()
            boundary = stretch_roundoff_tol * np.maximum(1.0, np.maximum(np.abs(stretches), critical))
            difference_to_threshold = stretches - critical
            robust_trigger = active_before & fail_permit & (difference_to_threshold > boundary)
            borderline = active_before & fail_permit & (np.abs(difference_to_threshold) <= boundary)
            allowed_event_edges = set(np.flatnonzero(robust_trigger).tolist())
            allowed_event_edges.update(edge for edge in event_edges if borderline[edge])
            missing = allowed_event_edges - event_edges
            unexpected = event_edges - allowed_event_edges
            if missing or unexpected:
                witness_edge = int(next(iter(missing or unexpected)))
                raise AuditHalt(
                    f"step {step}: coordinate-derived first-break transitions disagree with event CSV; "
                    f"missing={len(missing)}, unexpected={len(unexpected)}, witness={pairs[witness_edge].tolist()}, "
                    f"stretch={stretches[witness_edge]:.17g}, epsilon_c={critical[witness_edge]:.17g}")
            threshold_ambiguous_count += int(np.count_nonzero(borderline))
            actual_breaks = np.zeros(n_bonds, dtype=bool)
            if event_edges:
                actual_breaks[np.fromiter(event_edges, dtype=np.int64)] = True
            if np.any(actual_breaks & ~active_before):
                raise AuditHalt(f"step {step}: reactivation or duplicate break event in canonical roster")
            if np.any(robust_trigger & ~actual_breaks):
                raise AuditHalt(f"step {step}: threshold-triggered canonical transition is missing from event CSV")
            if np.any(actual_breaks & ~(robust_trigger | borderline)):
                raise AuditHalt(f"step {step}: event is below independently recomputed strict threshold")
            active[actual_breaks] = False
            for event in state_events:
                edge = int(event["edge_index"])
                event_tol_length = max(1e-15, 4e-12 * max(ref_length[edge], current_length[edge]))
                event_tol_stretch = max(1e-13, 4e-10 * max(abs(stretches[edge]), critical[edge]))
                if (abs(event["time_s"] - step * dt) > 1e-19
                    or abs(event["reference_length_m"] - ref_length[edge]) > event_tol_length
                    or abs(event["current_length_m"] - current_length[edge]) > event_tol_length
                    or abs(event["stretch"] - stretches[edge]) > event_tol_stretch
                    or abs(event["critical_stretch"] - critical[edge]) > 1e-13
                    or abs(event["stretch_over_critical"] - ratios[edge]) > 1e-9 * max(1.0, abs(ratios[edge]))):
                    raise AuditHalt(f"step {step}: event row not supported by deformed binary coordinates: {pairs[edge].tolist()}")
                event["computed_current_length_m"] = float(current_length[edge])
                event["computed_stretch"] = float(stretches[edge])
                event["computed_stretch_over_critical"] = float(ratios[edge])
                if len(first_break_witnesses) < 50:
                    first_break_witnesses.append({k: event[k] for k in (
                        "step", "time_s", "point_i", "point_j", "reference_length_m",
                        "computed_current_length_m", "computed_stretch", "critical_stretch",
                        "computed_stretch_over_critical")})
        else:
            pos = actual_coordinates[0]
            max_stretch_state[0] = float(np.max((np.linalg.norm(
                pos[edge_j0] - pos[edge_i0], axis=1) - ref_length) / ref_length))
            max_ratio_state[0] = float(np.max(
                ((np.linalg.norm(pos[edge_j0] - pos[edge_i0], axis=1) - ref_length) / ref_length)
                / critical))
        degree = np.bincount(pairs[active].ravel(), minlength=n_points + 1)[1:].astype(np.int32)
        active_pair_counts[state] = int(active.sum())
        broken_pair_counts[state] = n_bonds - int(active.sum())
        cumulative += len(state_events)
        if cumulative != int(broken_pair_counts[state]):
            raise AuditHalt(f"step {step}: event cumulative count disagrees with active canonical roster")
        new_event_counts[state] = len(state_events)
        state_count_field = n_active_tensor[:, state]
        state_damage_field = point_damage[:, state]
        if not np.array_equal(state_count_field.astype(np.int64), degree):
            geometry_active_count_mismatches += 1
            point_witness = int(np.flatnonzero(state_count_field.astype(np.int64) != degree)[0] + 1)
            raise AuditHalt(f"step {step}: binary active-bond count differs from independently reconstructed graph "
                            f"at point {point_witness}: binary={state_count_field[point_witness-1]}, "
                            f"recomputed={degree[point_witness-1]}")
        expected_damage = 1.0 - degree / reference_degrees
        damage_delta = np.abs(state_damage_field - expected_damage)
        if float(np.max(damage_delta)) > 2e-12:
            p = int(np.argmax(damage_delta)) + 1
            raise AuditHalt(f"step {step}: native damage not reproducible from active/reference degrees "
                            f"at point {p}; delta={damage_delta[p-1]:.6g}")
        if state > 0 and np.any(state_damage_field < point_damage[:, state - 1] - 2e-14):
            raise AuditHalt(f"step {step}: damage decreased despite irreversible native bond activity")
        # A union-find pass is needed only when the active canonical roster changes.
        # For repeated states with an identical roster, the prior partition is the
        # exact graph partition; selected-state bounds are refreshed below.
        if state == 0 or state_events:
            calc_labels, components = graph_components(active, pairs, mass,
                                                       actual_coordinates[state], left_mask, right_mask)
            cached_graph_labels = calc_labels
            cached_components = components
        else:
            if cached_graph_labels is None or cached_components is None:
                raise AuditHalt("graph cache was not initialized before unchanged topology state")
            calc_labels, components = cached_graph_labels, cached_components
        graph_labels_calc[:, state] = calc_labels
        if not np.array_equal(labels_binary[:, state], calc_labels):
            bad = int(np.flatnonzero(labels_binary[:, state] != calc_labels)[0]) + 1
            raise AuditHalt(f"step {step}: graph-label binary mismatch at point {bad}: "
                            f"binary={labels_binary[bad-1,state]}, recomputed={calc_labels[bad-1]}")
        graph_component_counts[state] = len(components)
        graph_largest_mass[state] = components[0]["mass_kg"]
        graph_largest_fraction[state] = components[0]["mass_fraction"]
        if state > 0 and graph_component_counts[state] < graph_component_counts[state - 1]:
            raise AuditHalt(f"step {step}: active-graph component count decreased under irreversible bond loss")
        if state > 0 and graph_largest_mass[state] > graph_largest_mass[state - 1] + 1e-14 * mass.sum():
            raise AuditHalt(f"step {step}: largest connected-component mass increased under bond loss")
        expected_active_from_first_events = active_ref_flag & (
            (event_step_by_edge < 0) | (event_step_by_edge > step))
        if not np.array_equal(active, expected_active_from_first_events):
            raise AuditHalt(f"step {step}: reconstructed active roster permits reactivation or disagrees with first-break steps")
        labels_on_right = tuple(int(x) for x in np.unique(calc_labels[right_mask]))
        labels_on_left = tuple(int(x) for x in np.unique(calc_labels[left_mask]))
        largest = components[0]
        largest_label = int(largest["label_min_point_id"])
        right_in_largest = int(np.count_nonzero(calc_labels[right_mask] == largest_label))
        right_all_largest = right_in_largest == int(right_mask.sum())
        right_wholly_outside_largest = right_in_largest == 0
        grip_membership_changed = (previous_right_labels is not None
                                   and labels_on_right != previous_right_labels)
        graph_count_changed = (state > 0
                               and graph_component_counts[state] != graph_component_counts[state - 1])
        graph_series.append({
            "step": state, "time_s": float(state * dt),
            "graph_component_count": int(graph_component_counts[state]),
            "active_canonical_pairs": int(active_pair_counts[state]),
            "broken_canonical_pairs": int(broken_pair_counts[state]),
            "new_break_events": int(new_event_counts[state]),
            "cumulative_break_events": int(cumulative),
            "largest_component_label_min_point_id": largest_label,
            "largest_component_mass_kg": float(largest["mass_kg"]),
            "largest_component_mass_fraction": float(largest["mass_fraction"]),
            "right_grip_component_labels": list(labels_on_right),
            "left_grip_component_labels": list(labels_on_left),
            "right_grip_points_in_largest": right_in_largest,
            "right_grip_point_count": int(right_mask.sum()),
            "right_grip_fraction_in_largest": float(right_in_largest / right_mask.sum()),
            "right_grip_all_in_largest": right_all_largest,
            "right_grip_wholly_outside_largest": right_wholly_outside_largest,
            "all_right_grip_points_share_a_graph_component": len(labels_on_right) == 1,
        })
        if state == 0 or state == n_steps or state_events or graph_count_changed or grip_membership_changed:
            current_bboxes = component_bboxes(calc_labels, actual_coordinates[state])
            selected = [c for c in components
                        if c["label_min_point_id"] == largest_label
                        or c["right_grip_point_count"] > 0
                        or c["left_grip_point_count"] > 0]
            selected.sort(key=lambda c: (-c["mass_kg"], c["label_min_point_id"]))
            selected_with_current_bounds = [
                {**c, "bbox_mm": current_bboxes[int(c["label_min_point_id"])]}
                for c in selected
            ]
            selected_component_states.append({
                "step": state, "time_s": float(state * dt),
                "graph_component_count": len(components),
                "largest_component_label_min_point_id": largest_label,
                "largest_component_mass_fraction": float(largest["mass_fraction"]),
                "right_grip_component_labels": list(labels_on_right),
                "left_grip_component_labels": list(labels_on_left),
                "grip_touching_component_detail_count": len(selected),
                "grip_touching_components_top_20": selected_with_current_bounds[:20],
                "component_detail_truncated": len(selected) > 20,
            })
        previous_right_labels = labels_on_right

    audit.require("first_break_events_recomputed_from_deformed_coordinates",
                  True, {"event_count": len(event_records),
                         "threshold_ambiguous_edge_state_checks": threshold_ambiguous_count,
                         "witnesses_recorded": len(first_break_witnesses)})
    audit.require("active_counts_damage_and_graph_labels_recomputed_each_state",
                  geometry_active_count_mismatches == 0
                  and graph_labels_calc.shape == labels_binary.shape,
                  {"states": n_states, "points": n_points,
                   "active_count_mismatch_states": geometry_active_count_mismatches,
                   "max_component_count": int(graph_component_counts.max()),
                   "final_component_count": int(graph_component_counts[-1])})

    graph_numeric_deltas: dict[str, float] = {}
    graph_int_keys = {
        "step": np.arange(n_states),
        "directed_bond_count": np.full(n_states, 2 * n_bonds),
        "active_directed_bonds": 2 * active_pair_counts,
        "canonical_pair_count": np.full(n_states, n_bonds),
        "active_canonical_pairs": active_pair_counts,
        "broken_canonical_pairs": broken_pair_counts,
        "new_break_events": new_event_counts,
        "cumulative_break_events": np.cumsum(new_event_counts),
        "graph_component_count": graph_component_counts,
    }
    for name, expected in graph_int_keys.items():
        actual = np.array([to_int(r, name) for r in graph_rows], dtype=np.int64)
        if not np.array_equal(actual, expected):
            bad = int(np.flatnonzero(actual != expected)[0])
            raise AuditHalt(f"graph history {name} disagrees at state {bad}: {actual[bad]} vs {expected[bad]}")
        graph_numeric_deltas[name] = 0.0
    graph_float_expected = {
        "time_s": np.arange(n_states, dtype=np.float64) * dt,
        "largest_component_mass_kg": graph_largest_mass,
        "largest_component_mass_fraction": graph_largest_fraction,
        "total_material_mass_kg": np.full(n_states, float(mass.sum())),
        "max_bond_stretch": max_stretch_state,
        "max_stretch_over_critical": max_ratio_state,
        "max_point_damage": np.max(point_damage, axis=0),
    }
    for name, expected in graph_float_expected.items():
        actual = np.array([to_float(r, name) for r in graph_rows], dtype=np.float64)
        delta = np.abs(actual - expected)
        graph_numeric_deltas[name] = float(delta.max(initial=0.0))
        if not np.allclose(actual, expected, rtol=2e-9, atol=2e-14):
            bad = int(np.flatnonzero(~np.isclose(actual, expected, rtol=2e-9, atol=2e-14))[0])
            raise AuditHalt(f"graph history {name} disagrees at step {bad}: "
                            f"CSV={actual[bad]:.17g}, recomputed={expected[bad]:.17g}")
    audit.require("graph_history_csv_matches_independent_topology",
                  True, {"columns": len(GRAPH_HEADER), "max_abs_numeric_delta": graph_numeric_deltas})

    coords_m = actual_coordinates
    damage_max_by_state = np.max(point_damage, axis=0)
    disp = tensor[:3]
    full_v = tensor[3:6].transpose(1, 0, 2)
    vhalf = tensor[6:9].transpose(1, 0, 2)
    accel = tensor[9:12].transpose(1, 0, 2)
    bint = tensor[12:15].transpose(1, 0, 2)
    bext = tensor[15:18].transpose(1, 0, 2)
    psi = tensor[18]
    times = np.arange(n_states, dtype=np.float64) * dt
    target_critical = float(loading["target_epsilon_c_multiplier"]) * formula_critical
    grip_span = xmax - xmin
    target_u_amplitude = target_critical * grip_span
    target_u = np.array([target_displacement(float(t), *steps["phase_durations_s"][:3],
                                              target_u_amplitude) for t in times])
    right_mean_ux = np.mean(disp[0, right_mask, :], axis=0)
    left_zero_max = float(np.max(np.abs(disp[:, left_mask, :])))
    right_yz_max = float(np.max(np.abs(disp[1:3, right_mask, :])))
    right_spread = float(np.max(np.abs(disp[0, right_mask, :] - right_mean_ux[None, :])))
    grip_path_err = float(np.max(np.abs(right_mean_ux - target_u)))
    path_tol = max(256.0 * np.finfo(float).eps * target_u_amplitude, 1e-14)
    audit.require("96_step_physical_time_secant_grip_path",
                  n_steps == 96 and abs(dt - EXPECTED_DT_S) == 0.0
                  and grip_path_err <= path_tol and right_spread <= path_tol
                  and left_zero_max <= path_tol and right_yz_max <= path_tol
                  and abs(float(right_mean_ux[-1])) <= path_tol,
                  {"n_steps": n_steps, "dt_s": dt, "total_time_s": float(times[-1]),
                   "target_displacement_m": target_u_amplitude,
                   "max_right_grip_path_error_m": grip_path_err,
                   "max_right_grip_point_spread_m": right_spread,
                   "max_left_grip_displacement_m": left_zero_max,
                   "max_right_yz_displacement_m": right_yz_max,
                   "tolerance_m": path_tol})
    expected_vhalf_rx = np.zeros(n_states, dtype=np.float64)
    for s in range(1, n_states):
        expected_vhalf_rx[s] = (target_displacement(float(times[s]), *steps["phase_durations_s"][:3],
                                                      target_u_amplitude)
                                - target_displacement(float(times[s] - dt), *steps["phase_durations_s"][:3],
                                                      target_u_amplitude)) / dt
    vhalf_right_err = float(np.max(np.abs(vhalf[right_mask, 0, :] - expected_vhalf_rx[None, :])))
    vhalf_fixed_err = max(float(np.max(np.abs(vhalf[left_mask, :, :]))),
                          float(np.max(np.abs(vhalf[right_mask, 1:3, :]))))
    vhalf_tol = max(1e-12, 2e-10 * float(np.max(np.abs(expected_vhalf_rx))))
    audit.require("captured_half_velocity_is_the_native_grip_secant",
                  vhalf_right_err <= vhalf_tol and vhalf_fixed_err <= vhalf_tol,
                  {"right_grip_secant_max_abs_error_m_s": vhalf_right_err,
                   "fixed_components_max_abs_vhalf_m_s": vhalf_fixed_err,
                   "tolerance_m_s": vhalf_tol})

    constrained = np.repeat((left_mask | right_mask)[:, None], 3, axis=1)
    free_mask = ~constrained
    kinetic_all = 0.5 * np.sum(mass[:, None, None] * full_v * full_v, axis=(0, 1))
    kinetic_free = 0.5 * np.sum(mass[:, None, None] * full_v * full_v
                                * free_mask[:, :, None], axis=(0, 1))
    kinetic_grip_raw = 0.5 * np.sum(mass[:, None, None] * full_v * full_v
                                    * constrained[:, :, None], axis=(0, 1))
    strain_energy = np.sum(psi * volume[:, None], axis=0)
    momentum = np.sum(mass[:, None, None] * full_v, axis=0)
    force_net = np.sum((bint + bext) * volume[:, None, None], axis=0)
    kinetic_partition_error = kinetic_all - kinetic_free - kinetic_grip_raw
    audit.require("raw_full_velocity_kinetic_energy_partition",
                  float(np.max(np.abs(kinetic_partition_error))) <= 2e-15,
                  {"max_abs_partition_delta_J": float(np.max(np.abs(kinetic_partition_error))),
                   "status": "raw solver diagnostic; constrained-grip physical KE UNKNOWN"})

    left_impulse = np.zeros((3, n_states), dtype=np.float64)
    right_impulse = np.zeros((3, n_states), dtype=np.float64)
    left_work = np.zeros(n_states, dtype=np.float64)
    right_work = np.zeros(n_states, dtype=np.float64)
    free_impulse_max = np.zeros(n_states, dtype=np.float64)
    momentum_residual = np.zeros((3, n_states), dtype=np.float64)
    if n_states > 1:
        v_before = full_v[:, :, :-1] + accel[:, :, :-1] * (dt / 2.0)
        v_after = vhalf[:, :, 1:]
        impulse_by_dof = mass[:, None, None] * (v_after - v_before)
        work_by_dof = 0.5 * mass[:, None, None] * (v_after * v_after - v_before * v_before)
        left_impulse[:, 1:] = np.sum(impulse_by_dof[left_mask, :, :], axis=0)
        right_impulse[:, 1:] = np.sum(impulse_by_dof[right_mask, :, :], axis=0)
        left_work[1:] = np.sum(work_by_dof[left_mask, :, :], axis=(0, 1))
        right_work[1:] = np.sum(work_by_dof[right_mask, :, :], axis=(0, 1))
        free_impulse_max[1:] = np.max(np.abs(impulse_by_dof[free_mask, :]), axis=0)
        force_impulse = (dt / 2.0) * (force_net[:, :-1] + force_net[:, 1:])
        total_bc_impulse = left_impulse[:, 1:] + right_impulse[:, 1:]
        momentum_residual[:, 1:] = (momentum[:, 1:] - momentum[:, :-1]
                                     - force_impulse - total_bc_impulse)
    bc_work_step = left_work + right_work
    bc_work_cumulative = np.cumsum(bc_work_step)
    energy_residual = strain_energy + kinetic_all - (strain_energy[0] + kinetic_all[0]) - bc_work_cumulative
    actual_left_ux = np.mean(disp[0, left_mask, :], axis=0)
    actuator_step = np.zeros(n_states, dtype=np.float64)
    if n_states > 1:
        right_force = right_impulse[0] / dt
        actuator_step[1:] = 0.5 * (right_force[:-1] + right_force[1:]) * np.diff(right_mean_ux)
    actuator_cumulative = np.cumsum(actuator_step)

    energy_expected: dict[str, np.ndarray] = {
        "time_s": times,
        "bc_impulse_interval_start_s": np.maximum(np.arange(n_states) - 1, 0) * dt,
        "bc_impulse_interval_end_s": times,
        "mass_kg": np.full(n_states, mass.sum()),
        "momentum_x_kg_m_s": momentum[0], "momentum_y_kg_m_s": momentum[1],
        "momentum_z_kg_m_s": momentum[2],
        "native_net_force_x_N": force_net[0], "native_net_force_y_N": force_net[1],
        "native_net_force_z_N": force_net[2],
        "osb_strain_energy_J": strain_energy, "kinetic_energy_J": kinetic_all,
        "free_dof_kinetic_energy_J": kinetic_free,
        "constrained_dof_raw_full_step_kinetic_energy_J": kinetic_grip_raw,
        "left_grip_u_x_mean_m": actual_left_ux,
        "right_grip_u_x_mean_m": right_mean_ux,
        "left_bc_impulse_x_N_s": left_impulse[0], "left_bc_impulse_y_N_s": left_impulse[1],
        "left_bc_impulse_z_N_s": left_impulse[2],
        "left_bc_force_average_x_N": left_impulse[0] / dt,
        "left_bc_force_average_y_N": left_impulse[1] / dt,
        "left_bc_force_average_z_N": left_impulse[2] / dt,
        "right_bc_impulse_x_N_s": right_impulse[0], "right_bc_impulse_y_N_s": right_impulse[1],
        "right_bc_impulse_z_N_s": right_impulse[2],
        "right_bc_force_average_x_N": right_impulse[0] / dt,
        "right_bc_force_average_y_N": right_impulse[1] / dt,
        "right_bc_force_average_z_N": right_impulse[2] / dt,
        "left_bc_discrete_kinetic_work_J": left_work,
        "right_bc_discrete_kinetic_work_J": right_work,
        "total_bc_discrete_kinetic_work_step_J": bc_work_step,
        "actuator_work_trapezoid_surrogate_step_J": actuator_step,
        "actuator_work_trapezoid_surrogate_cumulative_J": actuator_cumulative,
        "max_free_dof_impulse_residual_N_s": free_impulse_max,
        "momentum_update_residual_x_N_s": momentum_residual[0],
        "momentum_update_residual_y_N_s": momentum_residual[1],
        "momentum_update_residual_z_N_s": momentum_residual[2],
        "energy_balance_residual_diagnostic_J": energy_residual,
    }
    energy_deltas: dict[str, float] = {}
    for name in ENERGY_HEADER:
        if name == "step":
            actual = np.array([to_int(r, name) for r in energy_rows], dtype=np.int64)
            expected = np.arange(n_states, dtype=np.int64)
            good = np.array_equal(actual, expected)
            max_delta = 0.0 if good else float(np.max(np.abs(actual - expected)))
        else:
            actual = np.array([to_float(r, name) for r in energy_rows], dtype=np.float64)
            expected = energy_expected[name]
            max_delta = float(np.max(np.abs(actual - expected)))
            if name.endswith("_s"):
                rtol, atol = 1e-10, 1e-19
            elif name.endswith("_m"):
                rtol, atol = 2e-9, 1e-16
            elif name.endswith("_J") or name.endswith("_N_s") or name.endswith("_kg_m_s"):
                rtol, atol = 2e-9, 2e-16
            elif name.endswith("_N"):
                rtol, atol = 2e-9, 2e-14
            elif name == "mass_kg":
                rtol, atol = 1e-12, 1e-16
            else:
                rtol, atol = 2e-9, 2e-16
            good = bool(np.allclose(actual, expected, rtol=rtol, atol=atol))
        energy_deltas[name] = max_delta
        if not good:
            raise AuditHalt(f"energy/reaction CSV column {name} differs from independent binary recomputation; "
                            f"max_abs_delta={max_delta:.6g}")
    audit.require("energy_reaction_history_recomputed_from_binary",
                  True, {"rows": n_states, "numeric_columns": len(ENERGY_HEADER),
                         "max_abs_column_deltas": energy_deltas})

    # Independent time-level checks, including free updates and native velocity BC impulses.
    if n_states > 1:
        u_increment = tensor[:3, :, 1:] - tensor[:3, :, :-1]
        u_from_half = vhalf[:, :, 1:].transpose(1, 0, 2) * dt
        u_error = float(np.max(np.abs(u_increment - u_from_half)))
        free_vhalf_expected = full_v[:, :, :-1] + accel[:, :, :-1] * (dt / 2.0)
        free_vhalf_error = float(np.max(np.abs((vhalf[:, :, 1:] - free_vhalf_expected)[free_mask, :])))
        full_v_error = float(np.max(np.abs(full_v[:, :, 1:] -
                                            (vhalf[:, :, 1:] + accel[:, :, 1:] * (dt / 2.0)))))
        accel_expected = (bint + bext) / rho
        accel_error = float(np.max(np.abs(accel - accel_expected)))
    else:
        u_error = free_vhalf_error = full_v_error = accel_error = 0.0
        u_increment = np.zeros((3, n_points, 0), dtype=np.float64)
        accel_expected = np.zeros_like(accel)
    u_tol = max(2e-15, 2e-10 * float(np.max(np.abs(u_increment))) if u_increment.size else 0.0)
    velocity_scale = max(float(np.max(np.abs(full_v))), float(np.max(np.abs(vhalf))), 1e-30)
    velocity_tol = max(2e-12, 2e-10 * velocity_scale)
    acceleration_scale = max(float(np.max(np.abs(accel))), 1e-30)
    acceleration_tol = max(2e-8, 2e-10 * acceleration_scale)
    audit.require("velocity_verlet_state_time_levels_and_free_dof_updates",
                  u_error <= u_tol and free_vhalf_error <= velocity_tol
                  and full_v_error <= velocity_tol and accel_error <= acceleration_tol,
                  {"max_position_drift_error_m": u_error,
                   "max_free_vhalf_update_error_m_s": free_vhalf_error,
                   "max_full_velocity_reconstruction_error_m_s": full_v_error,
                   "max_acceleration_force_density_error_m_s2": accel_error,
                   "position_tolerance_m": u_tol,
                   "velocity_tolerance_m_s": velocity_tol,
                   "acceleration_tolerance_m_s2": acceleration_tol})
    audit.require("free_dof_impulse_is_zero_and_momentum_identity_recomputed",
                  float(np.max(free_impulse_max)) <= 2e-14
                  and float(np.max(np.abs(momentum_residual))) <= 2e-13,
                  {"max_free_dof_impulse_N_s": float(np.max(free_impulse_max)),
                   "max_abs_momentum_update_residual_N_s": float(np.max(np.abs(momentum_residual)))})

    # Event spatial summary uses the measured reference geometry and the actual state at each event.
    x_edges = np.linspace(xmin, xmax, 7)
    event_x_bins = np.zeros(6, dtype=np.int64)
    event_region_counts = {"touches_left_grip": 0, "touches_right_grip": 0,
                           "near_left_within_horizon": 0, "near_right_within_horizon": 0,
                           "interior_other": 0}
    event_spatial_witnesses = []
    for event in event_records:
        i0, j0 = event["point_i"] - 1, event["point_j"] - 1
        s = event["step"]
        mid_ref_x = float((xyz[i0, 0] + xyz[j0, 0]) / 2)
        b = min(5, max(0, int(np.searchsorted(x_edges, mid_ref_x, side="right") - 1)))
        event_x_bins[b] += 1
        touches_left = bool(left_mask[i0] or left_mask[j0])
        touches_right = bool(right_mask[i0] or right_mask[j0])
        near_left = bool(min(xyz[i0, 0], xyz[j0, 0]) - xmin <= delta_horizon + 1e-12)
        near_right = bool(xmax - max(xyz[i0, 0], xyz[j0, 0]) <= delta_horizon + 1e-12)
        event_region_counts["touches_left_grip"] += int(touches_left)
        event_region_counts["touches_right_grip"] += int(touches_right)
        event_region_counts["near_left_within_horizon"] += int(near_left)
        event_region_counts["near_right_within_horizon"] += int(near_right)
        event_region_counts["interior_other"] += int(not touches_left and not touches_right
                                                     and not near_left and not near_right)
        if len(event_spatial_witnesses) < 50:
            coords_at_break = actual_coordinates[s]
            p_i, p_j = coords_at_break[i0], coords_at_break[j0]
            midpoint = (p_i + p_j) / 2
            event_spatial_witnesses.append({
                "step": s, "time_s": float(times[s]), "point_i": event["point_i"],
                "point_j": event["point_j"], "reference_midpoint_x_mm": mid_ref_x * 1000.0,
                "event_midpoint_xyz_mm": (midpoint * 1000.0).tolist(),
                "touches_right_grip": touches_right,
                "touches_left_grip": touches_left,
                "within_one_horizon_of_right_edge": near_right,
                "stretch_over_critical": event["computed_stretch_over_critical"],
            })
    event_total = len(event_records)
    if event_total:
        event_fraction_by_x_bin = (event_x_bins / event_total).tolist()
    else:
        event_fraction_by_x_bin = [0.0] * 6
    right_detached_states = [int(s) for s, row in enumerate(graph_series)
                             if row["right_grip_wholly_outside_largest"]]
    right_attached_states = [int(s) for s, row in enumerate(graph_series)
                             if row["right_grip_all_in_largest"]]
    spatial_analysis = {
        "interpretation_boundary": "All labels here describe connected components of the active canonical-bond graph only; they do not identify physical fragments or a crack surface.",
        "first_break_event_count": event_total,
        "reference_x_bins_mm": (x_edges * 1000.0).tolist(),
        "first_break_event_counts_by_reference_midpoint_x_bin": event_x_bins.tolist(),
        "first_break_event_fractions_by_reference_midpoint_x_bin": event_fraction_by_x_bin,
        "first_break_event_location_counts": event_region_counts,
        "right_grip_wholly_outside_largest_graph_component_steps": right_detached_states,
        "right_grip_all_in_largest_graph_component_steps": right_attached_states,
        "graph_states_with_right_grip_membership": [r for r in graph_series
            if r["step"] in set([0, n_steps, *right_detached_states,
                                 *[s for s in range(1, n_states)
                                   if graph_component_counts[s] != graph_component_counts[s - 1]]])],
        "event_spatial_witnesses_first_50": event_spatial_witnesses,
        "event_subset_plot_policy": "If more than 500 first-break events exist, plot 500 deterministically evenly spaced event-order samples including first and last; machine audit retains all event rows and checks.",
    }

    # Explicitly quantify whether first breaks involve either end grip and which
    # reference x-layers the right-grip incident edges connect.  These are
    # pair/graph associations only, never crack-surface or fragment evidence.
    edge_left = left_mask[pairs[:, 0] - 1] | left_mask[pairs[:, 1] - 1]
    edge_right = right_mask[pairs[:, 0] - 1] | right_mask[pairs[:, 1] - 1]
    edge_both_right = right_mask[pairs[:, 0] - 1] & right_mask[pairs[:, 1] - 1]
    edge_cross_right = edge_right & ~edge_both_right
    pair_classes = {
        "all_canonical_pairs": np.ones(n_bonds, dtype=bool),
        "incident_to_left_grip": edge_left,
        "incident_to_right_grip": edge_right,
        "both_endpoints_in_right_grip": edge_both_right,
        "crosses_right_grip_boundary": edge_cross_right,
        "incident_to_either_end_grip": edge_left | edge_right,
    }
    first_break_edge_mask = event_step_by_edge >= 0
    pair_event_association = {}
    for class_name, class_mask in pair_classes.items():
        pair_count = int(np.count_nonzero(class_mask))
        event_count = int(np.count_nonzero(class_mask & first_break_edge_mask))
        pair_event_association[class_name] = {
            "canonical_pair_count": pair_count,
            "first_break_event_count": event_count,
            "fraction_of_class_pairs_broken": (event_count / pair_count if pair_count else None),
            "fraction_of_all_first_break_events": (event_count / event_total if event_total else None),
        }
    right_grip_mass_fraction = float(mass[right_mask].sum() / mass.sum())
    right_x_layer = int(np.max(grid_indexes[0]))
    right_layer_to_layer_events: Counter[tuple[int, int]] = Counter()
    right_layer_to_layer_pairs: Counter[tuple[int, int]] = Counter()
    right_layer_to_layer_steps: dict[tuple[int, int], list[int]] = {}
    for edge, (i0, j0) in enumerate(pairs):
        xi, xj = int(grid_indexes[0][i0]), int(grid_indexes[0][j0])
        key = (min(xi, xj), max(xi, xj))
        if xi == right_x_layer or xj == right_x_layer:
            right_layer_to_layer_pairs[key] += 1
        first_step = int(event_step_by_edge[edge])
        if first_step >= 0 and (xi == right_x_layer or xj == right_x_layer):
            right_layer_to_layer_events[key] += 1
            right_layer_to_layer_steps.setdefault(key, []).append(first_step)
    right_layer_pair_rows = []
    for key in sorted(set(right_layer_to_layer_pairs) | set(right_layer_to_layer_events)):
        steps_for_pair = right_layer_to_layer_steps.get(key, [])
        right_layer_pair_rows.append({
            "reference_x_layers_mm": [float(expected_axes[0][key[0]] * 1000.0),
                                      float(expected_axes[0][key[1]] * 1000.0)],
            "canonical_pair_count": int(right_layer_to_layer_pairs.get(key, 0)),
            "first_break_event_count": int(right_layer_to_layer_events.get(key, 0)),
            "fraction_of_pair_class_broken": (
                right_layer_to_layer_events[key] / right_layer_to_layer_pairs[key]
                if right_layer_to_layer_pairs.get(key, 0) else None),
            "first_event_step": min(steps_for_pair) if steps_for_pair else None,
            "last_event_step": max(steps_for_pair) if steps_for_pair else None,
            "distinct_event_steps": sorted(set(steps_for_pair)),
        })
    spatial_analysis["end_grip_pair_event_association"] = pair_event_association
    spatial_analysis["right_grip_mass_fraction"] = right_grip_mass_fraction
    spatial_analysis["right_grip_reference_x_mm"] = float(expected_axes[0][right_x_layer] * 1000.0)
    spatial_analysis["right_grip_incident_reference_layer_pairs"] = right_layer_pair_rows
    spatial_analysis["right_grip_incident_reference_layer_pair_method"] = (
        "Canonical pairs are counted once with sorted reference x-layer indices; event counts come from independently coordinate-recomputed first breaks."
    )

    check_indices = [0, n_steps]
    check_indices.extend(int(s) for s in np.flatnonzero(new_event_counts))
    check_indices.extend(int(s) for s in np.flatnonzero(np.diff(graph_component_counts) != 0) + 1)
    check_indices = sorted(set(check_indices))
    state_metrics = {
        "time_s": times.tolist(),
        "target_grip_displacement_m": target_u.tolist(),
        "actual_right_grip_displacement_m": right_mean_ux.tolist(),
        "active_canonical_pair_count": active_pair_counts.tolist(),
        "broken_canonical_pair_count": broken_pair_counts.tolist(),
        "new_break_events": new_event_counts.tolist(),
        "cumulative_break_events": np.cumsum(new_event_counts).tolist(),
        "graph_component_count": graph_component_counts.tolist(),
        "largest_component_mass_fraction": graph_largest_fraction.tolist(),
        "max_stretch_over_critical": max_ratio_state.tolist(),
        "max_point_damage": damage_max_by_state.tolist(),
        "kinetic_energy_raw_full_velocity_all_dof_J": kinetic_all.tolist(),
        "kinetic_energy_raw_full_velocity_free_dof_J": kinetic_free.tolist(),
        "kinetic_energy_raw_full_velocity_grip_masked_dof_J": kinetic_grip_raw.tolist(),
        "osb_strain_energy_J": strain_energy.tolist(),
        "cumulative_discrete_bc_kinetic_work_J": bc_work_cumulative.tolist(),
        "energy_balance_residual_diagnostic_J": energy_residual.tolist(),
    }
    report = {
        "schema": "pgfrag.free_fracture_pilot.independent_review.v1",
        "status": overall_review_status(audit),
        "meaning_of_pass": "Independent bytes/source bindings and numeric reconstruction passed; this says nothing about physical fragment validity, convergence, calibration, surface response, grid bias, or fracture dissipation.",
        "scope": "single completed free-fracture pilot attempt; no solver invocation",
        "source_bindings": {k: v for k, v in source_info.items() if k != "parsed_freeze"},
        "case_binding": {**source_binding, "schema": case.get("schema"), "case_id": case.get("case_id")},
        "launcher_provenance": {"receipt": receipt_info["binding"], "checks": receipt_info["checks"],
                                **provenance},
        "result_inventory": {"directory": str(RESULTS_DIR.resolve()),
                             "total_bytes": result_bytes, "files": result_files},
        "binary_layouts": {"point_state": {"metadata": point_binding,
                                             "array_shape_field_point_state": expected_shape,
                                             "field_names": STATE_FIELDS, "field_units": STATE_UNITS},
                           "graph_labels": {"metadata": graph_binding,
                                            "array_shape_point_state": [n_points, n_states]}},
        "geometry_and_roster": {
            "grid_shape": list(EXPECTED_GRID), "point_count": n_points,
            "canonical_bond_count": n_bonds, "directed_bond_count": 2 * n_bonds,
            "point_reference_digest": result_files["point_reference.csv"],
            "canonical_roster_digest": result_files["canonical_bonds.csv"],
            "event_csv_digest": result_files["native_failure_events.csv"],
            "mass_total_kg": float(mass.sum()), "box_shape_m": lengths_m.tolist(),
            "point_spacing_m": h, "horizon_m": delta_horizon,
            "critical_stretch_formula_value": formula_critical,
            "critical_stretch_roster_first": float(critical[0]),
            "left_grip_point_count": int(left_mask.sum()),
            "right_grip_point_count": int(right_mask.sum()),
            "reference_degree_min_max": [int(reference_degrees.min()), int(reference_degrees.max())],
            "reference_degree_sum": int(reference_degrees.sum()),
            "all_initial_canonical_pairs_active_and_permitted": bool(np.all(active_ref_flag & fail_permit)),
        },
        "reconstructed_history": {
            "dt_s": dt, "step_count": n_steps, "state_count": n_states,
            "phase_durations_s": steps["phase_durations_s"],
            "phase_interval_ceilings": steps["phase_intervals_ceiling"],
            "physical_end_time_s": float(times[-1]),
            "grip_target_displacement_m": target_u_amplitude,
            "event_count": event_total,
            "first_event_step": int(event_records[0]["step"]) if event_records else None,
            "first_event_time_s": float(event_records[0]["time_s"]) if event_records else None,
            "final_active_canonical_pairs": int(active_pair_counts[-1]),
            "final_broken_canonical_pairs": int(broken_pair_counts[-1]),
            "final_graph_component_count": int(graph_component_counts[-1]),
            "final_largest_graph_component_mass_fraction": float(graph_largest_fraction[-1]),
            "threshold_ambiguous_edge_state_checks_within_128eps": threshold_ambiguous_count,
            "first_break_witnesses_first_50": first_break_witnesses,
            "checked_states": check_indices,
            "state_series": state_metrics,
        "selected_graph_component_states": selected_component_states,
            "spatial_and_grip_association": spatial_analysis,
        },
        "independent_recomputation": {
            "recomputed_each_state_from_all_canonical_edges_and_binary_coordinates": True,
            "damage_definition": "1 - active_bond_count_at_state / reference_neighbor_count, recomputed from the irreversible active canonical roster",
            "strict_failure_rule": "fail_permit && stretch > critical_stretch; for threshold gaps within 128*epsilon(Float64) roundoff band, the event record is accepted only as an explicitly counted ambiguous comparison and its coordinates/length/stretch/threshold are still cross-checked",
            "graph_definition": "connected components of currently active canonical bonds; all point IDs are included, including graph-isolated points",
            "kinetic_energy_semantics": "raw 0.5 * mass * |stored full-step V|^2; free and grip-masked contributions are reported separately. Native BC constrains vhalf, while stored full-step V includes post-force acceleration half-step, so physical KE at constrained grips is UNKNOWN.",
            "energy_balance_status": "diagnostic only; no fracture-dissipation or energy-closure qualification",
            "energy_csv_max_abs_deltas": energy_deltas,
            "graph_csv_max_abs_deltas": graph_numeric_deltas,
            "grip_path_max_error_m": grip_path_err,
            "grip_secant_vhalf_max_error_m_s": vhalf_right_err,
            "max_free_dof_vhalf_update_error_m_s": free_vhalf_error,
            "max_position_drift_reconstruction_error_m": u_error,
            "max_full_velocity_reconstruction_error_m_s": full_v_error,
            "max_acceleration_force_density_error_m_s2": accel_error,
            "max_free_dof_bc_impulse_residual_N_s": float(np.max(free_impulse_max)),
            "max_momentum_identity_residual_N_s": float(np.max(np.abs(momentum_residual))),
            "final_raw_ke_all_dof_J": float(kinetic_all[-1]),
            "final_raw_ke_free_dof_J": float(kinetic_free[-1]),
            "final_raw_ke_grip_masked_dof_J": float(kinetic_grip_raw[-1]),
            "final_integrated_osb_strain_energy_J": float(strain_energy[-1]),
            "final_cumulative_discrete_bc_kinetic_work_J": float(bc_work_cumulative[-1]),
            "final_energy_balance_residual_diagnostic_J": float(energy_residual[-1]),
            "physical_grip_kinetic_energy_status": "UNKNOWN",
        },
        "scientific_status": {
            "physical_material_validity": "UNKNOWN_UNCALIBRATED",
            "physical_loading_validity": "UNKNOWN_UNCALIBRATED",
            "surface_response_and_NoCorrection_bias": "UNKNOWN_NOT_QUALIFIED",
            "grid_or_horizon_bias": "UNKNOWN_NOT_ASSESSED",
            "fracture_dissipation": "UNKNOWN_NOT_COMPUTED",
            "physical_crack_surface": "UNKNOWN_NOT_RECONSTRUCTED",
            "physical_fragment_claim": "NOT_ESTABLISHED",
            "active_graph_components_are_physical_fragments": False,
        },
        "checks": audit.checks,
        "native_arithmetic_validation": native_arithmetic_validation(audit),
    }
    make_figures(report, xyz, actual_coordinates, point_damage, right_mean_ux,
                 target_u, times, event_records, graph_series, kinetic_all,
                 kinetic_free, kinetic_grip_raw, strain_energy, bc_work_cumulative,
                 energy_residual, new_event_counts, cumulative_breaks=np.cumsum(new_event_counts),
                 max_ratio=max_ratio_state)
    return report


def make_figures(report: dict[str, Any], xyz: np.ndarray, actual_coordinates: np.ndarray,
                 damage: np.ndarray, right_mean_ux: np.ndarray, target_u: np.ndarray,
                 times: np.ndarray, events: list[dict[str, Any]], graph_series: list[dict[str, Any]],
                 kinetic_all: np.ndarray, kinetic_free: np.ndarray,
                 kinetic_grip_raw: np.ndarray, strain_energy: np.ndarray,
                 bc_work_cumulative: np.ndarray, energy_residual: np.ndarray,
                 new_event_counts: np.ndarray, cumulative_breaks: np.ndarray,
                 max_ratio: np.ndarray) -> None:
    require_owned_output()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    mpl_info = {"matplotlib": matplotlib.__version__, "numpy": np.__version__}
    report["visualization_runtime"] = mpl_info
    # Actual final coordinates, in millimetres, with no displacement magnification.
    xyz_final_mm = actual_coordinates[-1] * 1000.0
    damage_final = damage[:, -1]
    fig = plt.figure(figsize=(12, 9), constrained_layout=True)
    ax = fig.add_subplot(111, projection="3d")
    vmax = max(float(np.max(damage_final)), 1e-12)
    points = ax.scatter(xyz_final_mm[:, 0], xyz_final_mm[:, 1], xyz_final_mm[:, 2],
                        c=damage_final, cmap="magma", vmin=0.0, vmax=vmax,
                        s=7, alpha=0.9, depthshade=False, linewidths=0)
    plotted_indices: list[int] = []
    if len(events) <= 500:
        plotted_indices = list(range(len(events)))
    elif events:
        plotted_indices = sorted(set(np.linspace(0, len(events) - 1, 500, dtype=int).tolist()))
    for event_index in plotted_indices:
        event = events[event_index]
        state = int(event["step"])
        p_i = actual_coordinates[state, int(event["point_i"]) - 1] * 1000.0
        p_j = actual_coordinates[state, int(event["point_j"]) - 1] * 1000.0
        ax.plot([p_i[0], p_j[0]], [p_i[1], p_j[1]], [p_i[2], p_j[2]],
                color="#00b7ff", linewidth=0.75, alpha=0.72)
    cb = fig.colorbar(points, ax=ax, pad=0.10, shrink=0.72)
    cb.set_label("Native point damage (fraction of broken reference neighbors)")
    mins = xyz_final_mm.min(axis=0)
    maxs = xyz_final_mm.max(axis=0)
    spans = np.maximum(maxs - mins, 1e-9)
    padding = 0.03 * spans
    ax.set_xlim(mins[0] - padding[0], maxs[0] + padding[0])
    ax.set_ylim(mins[1] - padding[1], maxs[1] + padding[1])
    ax.set_zlim(mins[2] - padding[2], maxs[2] + padding[2])
    ax.set_box_aspect(tuple(spans.tolist()))
    ax.set_xlabel("x actual position (mm)")
    ax.set_ylabel("y actual position (mm)")
    ax.set_zlabel("z actual position (mm)")
    ax.view_init(elev=20, azim=-57)
    ax.set_title("Free-fracture pilot: terminal point damage and first-break bond locations\n"
                 "Point cloud is the actual terminal configuration at equal-mm scale; cyan bonds are sampled at their own first-break states.")
    if not events:
        ax.text2D(0.03, 0.03, "No native threshold events recorded in this completed output.",
                  transform=ax.transAxes)
    elif len(events) > len(plotted_indices):
        ax.text2D(0.03, 0.03,
                  f"Showing {len(plotted_indices)} of {len(events)} events: evenly spaced event-order samples incl. endpoints; no surface reconstructed.",
                  transform=ax.transAxes)
    fig_path = OUTPUT_DIR / "pilot_geometry_damage_first_breaks.png"
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)

    # A compact diagnostic timeline; no energy closure/fracture dissipation claim.
    t_us = times * 1e6
    fig, axes = plt.subplots(4, 1, figsize=(13, 14), constrained_layout=True)
    ax = axes[0]
    ax.plot(t_us, right_mean_ux * 1e3, color="#1f77b4", label="measured right-grip mean Ux")
    ax.plot(t_us, target_u * 1e3, color="black", linestyle="--", linewidth=1,
            label="continuous-time target Ux")
    ax.set_ylabel("grip displacement (mm)")
    ax.set_xlabel("physical time (µs)")
    ax.grid(alpha=0.25)
    ax2 = ax.twinx()
    ax2.bar(t_us, new_event_counts, width=float(times[1] * 1e6) * 0.85 if len(times) > 1 else 0.1,
            color="#d62728", alpha=0.35, label="new first-break events / step")
    ax2.set_ylabel("new first-break events per step")
    lines, labels = ax.get_legend_handles_labels()
    bars, bar_labels = ax2.get_legend_handles_labels()
    ax.legend(lines + bars, labels + bar_labels, loc="upper left", fontsize=8)
    ax.set_title("Native grip path and first-break event timeline")

    ax = axes[1]
    component_count = np.array([r["graph_component_count"] for r in graph_series])
    largest_fraction = np.array([r["largest_component_mass_fraction"] for r in graph_series])
    ax.step(t_us, component_count, where="post", color="#9467bd", label="active-graph component count")
    ax.set_ylabel("active graph components")
    ax.set_xlabel("physical time (µs)")
    ax.grid(alpha=0.25)
    ax.legend(loc="best", fontsize=8)
    ax.set_title("Active-graph connectivity count (not physical fragment count)")

    ax = axes[2]
    right_grip_largest_fraction = np.array(
        [r["right_grip_fraction_in_largest"] for r in graph_series], dtype=np.float64)
    ax.plot(t_us, largest_fraction, color="#2ca02c", label="largest component material-mass fraction")
    ax.plot(t_us, right_grip_largest_fraction, color="#ff7f0e", linestyle="--",
            label="right-grip point fraction in largest component")
    ax.set_ylim(-0.03, 1.03)
    ax.set_ylabel("fraction (0–1)")
    ax.set_xlabel("physical time (µs)")
    ax.grid(alpha=0.25)
    ax.legend(loc="best", fontsize=8)
    ax.set_title("Largest-component mass share and right-grip membership")

    ax = axes[3]
    to_microjoule = 1e6
    ax.plot(t_us, kinetic_all * to_microjoule, label="raw full-step KE, all DOFs")
    ax.plot(t_us, kinetic_free * to_microjoule, label="raw full-step KE, free DOFs")
    ax.plot(t_us, kinetic_grip_raw * to_microjoule, label="raw KE, grip-masked DOFs")
    ax.plot(t_us, strain_energy * to_microjoule, label="integrated native OSB strain energy")
    ax.plot(t_us, bc_work_cumulative * to_microjoule, label="cumulative discrete BC kinetic work")
    ax.plot(t_us, energy_residual * to_microjoule, linestyle="--",
            label="energy-balance diagnostic residual")
    ax.set_xlabel("physical time (µs)")
    ax.set_ylabel("raw diagnostic energy (µJ)")
    ax.grid(alpha=0.25)
    ax.legend(loc="best", fontsize=8, ncol=2)
    ax.set_title("Raw solver energies and BC work; constrained-grip physical KE and fracture dissipation remain UNKNOWN")
    fig_path_timeline = OUTPUT_DIR / "pilot_timeline_graph_energy.png"
    fig.savefig(fig_path_timeline, dpi=170)
    plt.close(fig)
    report["visualizations"] = {
        fig_path.name: file_binding(fig_path),
        fig_path_timeline.name: file_binding(fig_path_timeline),
        "geometry_event_sampling": {
            "total_first_break_events": len(events),
            "plotted_events": len(plotted_indices),
            "method": "all events if <=500; otherwise 500 evenly spaced event-order samples including first and last",
            "deformation_scale": "actual coordinates in millimetres; no magnification",
            "surface_reconstruction": False,
        },
    }


def make_markdown(report: dict[str, Any]) -> str:
    status = report.get("status", "BLOCKED")
    native_arithmetic = report.get("native_arithmetic_validation", {})
    rec = report.get("reconstructed_history", {})
    independ = report.get("independent_recomputation", {})
    spatial = rec.get("spatial_and_grip_association", {})
    receipt_check = report.get("checks", {}).get(RECEIPT_SIZE_DISCREPANCY_CHECK, {})
    receipt_evidence = receipt_check.get("evidence", {})
    if receipt_check and not receipt_check.get("pass", False):
        receipt_note = (
            "\n\nReceipt bookkeeping discrepancy (exact, zero-byte tolerance): "
            f"the receipt records {receipt_evidence.get('receipt_recorded_bytes')} attempt bytes, "
            f"while the independently counted completed attempt tree has "
            f"{receipt_evidence.get('independently_observed_final_attempt_tree_bytes')} bytes "
            f"(observed − recorded = {receipt_evidence.get('signed_delta_observed_minus_receipt')} bytes). "
            "The original receipt remains unchanged; numeric review status is shown separately."
        )
    else:
        receipt_note = ""
    return f"""# Independent free-fracture pilot review

Independent audit status: **{status}**. This status covers source/result byte binding and numerical reconstruction only. It does not establish physical fragments, a crack surface, convergence, calibrated material response, surface validity, grid/horizon independence, or fracture dissipation.

Native arithmetic validation status: **{native_arithmetic.get('status', 'NOT_COMPLETED')}** ({native_arithmetic.get('check_count', 0)}/{native_arithmetic.get('required_check_count', 0)} required independent numeric checks completed). This is reported separately from the execution-receipt bookkeeping status.
{receipt_note}

The audit uses the completed one-shot receipt as an execution/provenance gate, verifies frozen source pins and output hashes, then reconstructs activity directly from the little-endian state tensor, canonical bond roster, point coordinates, and event threshold rule. Producer summary numbers are hash-bound for provenance but are not used as numerical oracles.

The state and graph tensors are independently decoded in Julia column-major order and checked against the sidecar field names, units, counts, shapes, byte totals, and SHA-256 digests. At every state the reviewer recomputes each active pair's stretch from the recorded deformed point coordinates; first breaks must be permitted and strictly beyond the CriticalStretch threshold, except that events within a 128-ULP comparison band are explicitly counted as roundoff-ambiguous. It then rebuilds per-point active degree and damage, irreversible transitions, union-find components, minimum-point-ID graph labels, component mass, event totals, and graph CSV metrics.

The native grip path is checked at all 97 state times against the fixed physical-time cubic ramp/hold/unload/dwell profile. The right grip's stored half-step velocity is checked against the secant of that path; free displacement drift, free half-step velocity, stored full velocity, acceleration/force-density relation, discrete grip impulses and work, momentum update, and the energy/reaction CSV are independently recomputed. Raw full-step kinetic energy is split into free and grip-masked DOFs. Since the native boundary condition constrains vhalf and the stored full-step velocity includes the post-force acceleration half-step, physical kinetic energy at constrained grip DOFs remains **UNKNOWN**.

Graph labels describe only connectivity in the active canonical-bond graph. A component count, even one isolated from an end grip, is not evidence of a physical fragment or an internal crack. Event locations and component bounding boxes are shown only as actual geometry and actual sampled first-break bond segments; the reviewer fits no crack plane and invents no surface. Material calibration, physical loading validity, NoCorrection free-surface response, grid/horizon bias, and fracture dissipation remain **UNKNOWN** or **NOT COMPUTED**.

The geometry figure is {report.get('visualizations', {}).get('pilot_geometry_damage_first_breaks.png', {}).get('path', 'pending')} and uses equal world-unit millimetres with no displacement exaggeration. The timeline/graph/raw-energy figure is {report.get('visualizations', {}).get('pilot_timeline_graph_energy.png', {}).get('path', 'pending')}.
"""


def write_outputs(audit: AuditState, report: dict[str, Any] | None,
                 fatal: str | None = None) -> None:
    require_owned_output()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    payload: dict[str, Any]
    if report is None:
        payload = {
            "schema": "pgfrag.free_fracture_pilot.independent_review.v1",
            "status": "FAIL" if audit.errors else "BLOCKED",
            "fatal_error": fatal,
            "meaning_of_pass": "No output state was independently accepted; see fatal_error and checks.",
            "checks": audit.checks,
            "scientific_status": {
                "physical_fragment_claim": "NOT_ESTABLISHED",
                "surface_response": "UNKNOWN",
                "grid_bias": "UNKNOWN",
                "fracture_dissipation": "NOT_COMPUTED",
            },
        }
    else:
        payload = report
        payload["checks"] = audit.checks
        payload["status"] = overall_review_status(audit)
        payload["native_arithmetic_validation"] = native_arithmetic_validation(audit)
    json_path = OUTPUT_DIR / "pilot_independent_review.json"
    with json_path.open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(payload, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    md_path = OUTPUT_DIR / "PILOT_REVIEW_REPORT.md"
    md_path.write_text(make_markdown(payload), encoding="utf-8", newline="\n")


def main() -> int:
    parser = argparse.ArgumentParser(description="Read-only PG free-fracture pilot independent audit")
    parser.add_argument("--root-authorized", action="store_true",
                        help="explicitly supplied only after the parent authorizes the final read-only audit")
    parser.add_argument("--expected-final-summary-sha256", default="",
                        help="author-provided final summary SHA-256; checked against receipt before result files are opened")
    args = parser.parse_args()
    if not args.root_authorized:
        print("AUDIT_NOT_STARTED: explicit parent authorization flag is required", file=sys.stderr)
        return 3
    expected_sha = args.expected_final_summary_sha256.strip().lower()
    if len(expected_sha) != 64 or any(char not in "0123456789abcdef" for char in expected_sha):
        print("AUDIT_NOT_STARTED: a 64-character final summary SHA-256 is required", file=sys.stderr)
        return 3
    require_owned_output()
    audit = AuditState()
    try:
        # Receipt first: a partial/failed one-shot never causes result arrays to be opened.
        receipt_info = verify_receipt_before_results(audit, expected_sha)
        report = audit_pilot(audit, receipt_info)
        write_outputs(audit, report)
        print(f"INDEPENDENT_PILOT_AUDIT_{report['status']} {OUTPUT_DIR}")
        completed_statuses = {"PASS", "AUDITED_WITH_RECEIPT_SIZE_DISCREPANCY"}
        return 0 if report["status"] in completed_statuses else 2
    except Exception as exc:
        detail = f"{type(exc).__name__}: {exc}"
        if not isinstance(exc, AuditHalt):
            detail += "\n" + traceback.format_exc()
        write_outputs(audit, None, detail)
        print(f"INDEPENDENT_PILOT_AUDIT_BLOCKED {detail}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
