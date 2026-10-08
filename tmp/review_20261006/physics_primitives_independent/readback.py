"""One bounded, implementation-independent CSV and pin readback."""
import csv
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(r"C:\Users\xiemo\Desktop\ShardForge")
OUT = ROOT / "pgfrag/tmp/review_20261006/physics_primitives_independent"
COH = ROOT / "pgfrag/tmp/review_20261006/cohesive_scalar_law_20261006"
FIELD = ROOT / "pgfrag/tmp/review_20261006/rest_material_field_20261006"
CONTACT = ROOT / "pgfrag/tmp/review_20261006/rigid_plane_contact_20261006"

def read_csv(path):
    with path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    return rows

def max_abs(values):
    return max((abs(x) for x in values), default=0.0)

# Cohesive ledger: recompute running maximum, closed-form damage, stored and
# dissipated energies, and traction from CSV state columns and receipt constants.
coh_rows = read_csv(COH / "cohesive_history.csv")
coh_receipt = json.loads((COH / "run_receipt.json").read_text(encoding="utf-8"))
p = coh_receipt["parameters"]
k0, tp, gc = p["k0_pa_per_m"], p["tpeak_pa"], p["gc_j_per_m2"]
d0, df = tp / k0, 2.0 * gc / tp
running_max = 0.0
kappa_err, damage_err, stored_err, diss_err, traction_err, ledger_err = [], [], [], [], [], []
damage_values, diss_values = [], []
for row in coh_rows:
    opening = float(row["opening_m"])
    running_max = max(running_max, opening)
    kappa = running_max
    damage = 0.0 if kappa <= d0 else (1.0 if kappa >= df else 1.0 - tp * (df-kappa) / (df-d0) / (k0*kappa))
    stored = 0.5 * k0 * (1.0-damage) * opening * opening
    diss = 0.0 if kappa <= d0 else (gc if kappa >= df else gc*(kappa-d0)/(df-d0))
    traction = k0 * (1.0-damage) * opening
    kappa_err.append(float(row["kappa_m"])-kappa)
    damage_err.append(float(row["damage"])-damage)
    stored_err.append(float(row["stored_energy_j_per_m2"])-stored)
    diss_err.append(float(row["dissipated_energy_j_per_m2"])-diss)
    traction_err.append(float(row["traction_pa"])-traction)
    ledger_err.append(float(row["external_work_j_per_m2"])-float(row["stored_energy_j_per_m2"])-float(row["dissipated_energy_j_per_m2"]))
    damage_values.append(float(row["damage"]))
    diss_values.append(float(row["dissipated_energy_j_per_m2"]))
coh_summary = {
    "rows": len(coh_rows), "parameters": {"delta0_m": d0, "deltaf_m": df},
    "max_abs_csv_minus_independent_kappa_m": max_abs(kappa_err),
    "max_abs_csv_minus_independent_damage": max_abs(damage_err),
    "max_abs_csv_minus_independent_stored_J_m2": max_abs(stored_err),
    "max_abs_csv_minus_independent_dissipation_J_m2": max_abs(diss_err),
    "max_abs_csv_minus_independent_traction_Pa": max_abs(traction_err),
    "max_abs_psi_plus_D_minus_work_J_m2": max_abs(ledger_err),
    "minimum_damage_increment": min((b-a for a,b in zip(damage_values, damage_values[1:])), default=0.0),
    "minimum_dissipation_increment_J_m2": min((b-a for a,b in zip(diss_values, diss_values[1:])), default=0.0),
    "final_kappa_m": float(coh_rows[-1]["kappa_m"]), "final_damage": float(coh_rows[-1]["damage"]),
    "final_dissipation_J_m2": float(coh_rows[-1]["dissipated_energy_j_per_m2"]),
}

# Field CSV: independently evaluate the stated finite Fourier sum from detached
# manifest coefficients, without importing or invoking the field implementation.
field_rows = read_csv(FIELD / "field_samples.csv")
manifest = json.loads((FIELD / "field_manifest.json").read_text(encoding="utf-8"))
amp = manifest["parameters"]["amplitude_max_relative_deviation"]
modes = manifest["random_modes"]
field_errors, field_values = [], []
for row in field_rows:
    xyz = (float(row["x_m"]), float(row["y_m"]), float(row["z_m"]))
    value = amp * sum(math.cos(sum(k*x for k,x in zip(mode["wavevector_rad_per_m"], xyz)) + mode["phase_rad"]) for mode in modes) / len(modes)
    field_values.append(float(row["heterogeneity_dimensionless"]))
    field_errors.append(float(row["heterogeneity_dimensionless"])-value)
spec = manifest["spectrum_and_correlation"]
q = spec["shell_root_q"]
k = spec["wave_number_magnitude_rad_per_m"]
ell = spec["configured_ensemble_e_fold_length_m"]
field_summary = {
    "rows": len(field_rows), "mode_count": len(modes),
    "max_abs_csv_minus_manifest_fourier_sum": max_abs(field_errors),
    "sample_min": min(field_values), "sample_max": max(field_values),
    "grid_unique_x": len({r["x_m"] for r in field_rows}), "grid_unique_y": len({r["y_m"] for r in field_rows}),
    "all_z_zero": all(float(r["z_m"]) == 0.0 for r in field_rows),
    "wave_number_times_e_fold_length": k*ell,
    "shell_root_q": q,
    "sin_q_over_q": math.sin(q)/q,
    "exp_minus_one": math.exp(-1.0),
    "finite_realization_covariance_claim": spec["finite_realization_caveat"],
}

# Contact saved trajectories: recompute signed plane spring law, energies,
# trapezoidal impulse, and momentum balance from positions/velocities only.
contact_rows = read_csv(CONTACT / "actual_trajectory_force_energy_impulse.csv")
contact_groups = {}
for row in contact_rows:
    contact_groups.setdefault(row["case_id"], []).append(row)
contact_case_results = []
for case_id, rows in sorted(contact_groups.items()):
    rows.sort(key=lambda row: int(row["step"]))
    forces, impulse = [], 0.0
    force_err, potential_err, kinetic_err, total_err, gap_err, saved_impulse_err = [], [], [], [], [], []
    previous = None
    for row in rows:
        y, vy, t = float(row["y_m"]), float(row["vy_m_per_s"]), float(row["time_s"])
        gap = y
        force = -100.0 * min(gap, 0.0)
        potential = 0.5 * 100.0 * min(gap, 0.0)**2
        kinetic = 0.5 * vy*vy
        if previous is not None:
            dt = t - float(previous["time_s"])
            impulse += 0.5 * (forces[-1] + force) * dt
        gap_err.append(float(row["signed_gap_m"])-gap)
        force_err.append(float(row["contact_force_y_N"])-force)
        potential_err.append(float(row["contact_potential_J"])-potential)
        kinetic_err.append(float(row["kinetic_J"])-kinetic)
        total_err.append(float(row["total_energy_J"])-(kinetic+potential))
        saved_impulse_err.append(float(row["cumulative_contact_impulse_Ns"])-impulse)
        forces.append(force)
        previous = row
    delta_p = float(rows[-1]["vy_m_per_s"]) - float(rows[0]["vy_m_per_s"])
    contact_case_results.append({
        "case_id": case_id, "rows": len(rows), "dt_s": float(rows[1]["time_s"])-float(rows[0]["time_s"]),
        "max_abs_gap_identity_m": max_abs(gap_err), "max_abs_force_error_N": max_abs(force_err),
        "max_abs_potential_error_J": max_abs(potential_err), "max_abs_kinetic_error_J": max_abs(kinetic_err),
        "max_abs_total_energy_identity_error_J": max_abs(total_err),
        "max_abs_saved_vs_independent_impulse_Ns": max_abs(saved_impulse_err),
        "independent_net_impulse_Ns": impulse, "point_delta_momentum_Ns": delta_p,
        "momentum_impulse_balance_error_Ns": delta_p-impulse,
        "max_abs_total_energy_residual_vs_2J": max(abs(float(r["total_energy_J"])-2.0) for r in rows),
    })

# Check that the focused contact receipt's source pins describe both the frozen
# source snapshot and present source, and that existing trajectory bytes remain pinned.
pin_doc = json.loads((CONTACT / "focused_readonly_check/focused_source_and_prior_artifact_pins.json").read_text(encoding="utf-8"))
snapshot = CONTACT / "focused_readonly_check/source_snapshot"
current_sources = {
    "rigid_plane_contact.py": ROOT / "pgfrag/experiments/rigid_plane_contact_20261006/rigid_plane_contact.py",
    "run_toy_diagnostic.py": ROOT / "pgfrag/experiments/rigid_plane_contact_20261006/run_toy_diagnostic.py",
    "run_focused_checks.py": ROOT / "pgfrag/experiments/rigid_plane_contact_20261006/run_focused_checks.py",
}
source_pin_results = {}
for name, path in current_sources.items():
    expected = pin_doc["focused_sources_sha256"][name]
    source_pin_results[name] = {
        "snapshot_matches_pin": hashlib.sha256((snapshot/name).read_bytes()).hexdigest() == expected,
        "current_matches_focused_snapshot": hashlib.sha256(path.read_bytes()).hexdigest() == expected,
    }
artifact_pin_results = {}
for name, details in pin_doc["prior_run_artifacts_current_sha256"].items():
    path = CONTACT/name
    artifact_pin_results[name] = {
        "current_matches_focused_pin": hashlib.sha256(path.read_bytes()).hexdigest() == details["sha256_at_focused_check"],
        "history_note": details["history_note"],
    }

csv_rows = len(coh_rows) + len(field_rows) + len(contact_rows)
summary = {
    "scope": "read-only CSV/manifests/source-pin readback; no implementation import/call, solver, CAD, random fixture, or trajectory rerun",
    "csv_data_rows_read": csv_rows,
    "csv_row_limit": 2000,
    "cohesive": coh_summary,
    "rest_material_field": field_summary,
    "rigid_plane_contact": {
        "total_trajectory_rows": len(contact_rows), "cases": contact_case_results,
        "source_pin_checks": source_pin_results, "prior_artifact_pin_checks": artifact_pin_results,
        "focused_receipt_passed_groups": json.loads((CONTACT/"focused_readonly_check/focused_check_receipt.json").read_text(encoding="utf-8"))["passed_check_count"],
        "original_run_source_snapshot_available": False,
        "comparison_errors_pre_correction_bytes_available": False,
    },
}
assert csv_rows <= 2000
(OUT / "independent_readback.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(json.dumps(summary, ensure_ascii=False, indent=2))
