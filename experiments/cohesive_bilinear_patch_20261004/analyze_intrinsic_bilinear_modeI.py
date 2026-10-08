from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

OUT = Path("/mnt/c/Users/xiemo/Desktop/ShardForge/relicweave-synlab/tmp/review_20261004/implementation/cohesive_bilinear_patch_20261004/controlled_intrinsic_bilinear_modeI_201_run01")
SIGMA = 70.0e6
GC = 7.6
DELTA_0 = 2.1714285714285714e-8
ENERGY_KEYS = ("kinetic", "potential", "dissipated", "reversible", "cohesive contact")


def analytic_traction(s):
    dc = 2.0 * GC / SIGMA
    if s <= DELTA_0:
        return SIGMA * s / DELTA_0
    if s <= dc:
        return SIGMA * (dc - s) / (dc - DELTA_0)
    return 0.0


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def analyze():
    status = json.loads((OUT / "run_status.json").read_text(encoding="utf-8"))
    hist_path = OUT / "native_history_partial_or_complete.npz"
    if not hist_path.is_file():
        summary = {"run_status": status.get("status"), "solver_steps_completed": status.get("solver_steps_completed", 0),
                   "analysis_status": "NO_NATIVE_HISTORY", "plot_status": "NOT_POSSIBLE"}
        write_json(OUT / "analysis_201.json", summary)
        return summary
    with np.load(hist_path, allow_pickle=False) as z:
        s = z["commanded_opening_m"]
        opening, traction = z["opening"], z["tractions"]
        damage, delta_max = z["damage"], z["delta_max"]
        energies = z["energies_J"]
        counts, mass_sums = z["counts"], z["mass_component_sums_kg"]
        force_z = z["upper_internal_force_z_sum_raw_N"]
        boundary_error = z["post_solve_boundary_readback_error_m"]
        rest = z["rest_nodes"][0]
        coh = z["cohesive_connectivity"][0]
        pidx, solver_step = z["point_index"], z["solver_step"]
        disp = z["displacement"]
    n = len(s)
    area_cell = np.asarray([
        0.5 * np.linalg.norm(np.cross(rest[row[1]] - rest[row[0]], rest[row[2]] - rest[row[0]]))
        for row in coh
    ])
    area = float(area_cell.sum())
    gc_area = GC * area
    opening_norm = np.linalg.norm(opening, axis=2)
    traction_norm = np.linalg.norm(traction, axis=2)
    expected = np.asarray([analytic_traction(float(x)) for x in s])
    opening_err = np.abs(opening_norm - s[:, None])
    traction_err = np.abs(traction_norm - expected[:, None])
    d = energies[:, 2]
    r = energies[:, 3]
    contact = energies[:, 4]
    kinetic = energies[:, 0]
    potential = energies[:, 1]
    finite_damage = bool(np.isfinite(damage).all())
    damage_bounds = finite_damage and bool(np.all((damage >= 0.0) & (damage <= 1.0)))
    damage_monotonic = finite_damage and bool(np.all(np.diff(damage, axis=0) >= -1e-12))
    finite_delta = bool(np.isfinite(delta_max).all())
    delta_monotonic = finite_delta and bool(np.all(np.diff(delta_max, axis=0) >= -1e-12))
    damage_means = damage.mean(axis=1)
    delta_means = delta_max.mean(axis=1)
    damage_peak_i = int(np.argmax(damage_means))
    delta_peak_i = int(np.argmax(delta_means))
    damage_cell_spreads = np.ptp(damage, axis=1)
    delta_cell_spreads = np.ptp(delta_max, axis=1)
    opening_spread = np.ptp(opening_norm, axis=1)
    traction_spread = np.ptp(traction_norm, axis=1)
    dot = np.einsum("tci,tci->tc", opening, traction)
    signed_work = 0.0
    for i in range(1, n):
        delta_vec = opening[i] - opening[i - 1]
        mean_t = 0.5 * (traction[i] + traction[i - 1])
        signed_work += float(np.sum(area_cell * np.einsum("ci,ci->c", mean_t, delta_vec)))
    initial = np.load(OUT / "initial_state.npz", allow_pickle=False)
    initial_rest = initial["rest_nodes"]
    initial_coh = initial["cohesive_connectivity"]
    faces_a = (initial_rest[initial_coh[:, :3], :2] * 1000.0).tolist()
    faces_b = (initial_rest[initial_coh[:, 3:], :2] * 1000.0).tolist()
    plot_data = {
        "opening_nm": (opening_norm.mean(axis=1) * 1e9).tolist(),
        "commanded_nm": (s * 1e9).tolist(),
        "traction_mean_MPa": (traction_norm.mean(axis=1) / 1e6).tolist(),
        "traction_min_MPa": (traction_norm.min(axis=1) / 1e6).tolist(),
        "traction_max_MPa": (traction_norm.max(axis=1) / 1e6).tolist(),
        "analytic_traction_MPa": (expected / 1e6).tolist(),
        "dissipated_uJ": (d * 1e6).tolist(),
        "reversible_uJ": (r * 1e6).tolist(),
        "GcA_uJ": float(gc_area * 1e6),
        "interface_side_A_xy_mm": faces_a,
        "interface_side_B_xy_mm": faces_b,
        "actual_plot_points": n,
    }
    write_json(OUT / "plot_data.json", plot_data)
    with (OUT / "damage_delta_history_201.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["point_index", "solver_step", "command_nm", "damage_mean", "damage_min", "damage_max",
                         "damage_cell_spread", "delta_max_mean_nm", "delta_max_min_nm", "delta_max_max_nm",
                         "delta_max_cell_spread_nm"] +
                        [f"damage_cell_{i}" for i in range(damage.shape[1])] +
                        [f"delta_max_cell_{i}_nm" for i in range(delta_max.shape[1])])
        for i in range(n):
            writer.writerow([int(pidx[i]), int(solver_step[i]), float(s[i] * 1e9), float(damage_means[i]),
                             float(damage[i].min()), float(damage[i].max()), float(damage_cell_spreads[i]),
                             float(delta_means[i] * 1e9), float(delta_max[i].min() * 1e9),
                             float(delta_max[i].max() * 1e9), float(delta_cell_spreads[i] * 1e9)] +
                            damage[i].tolist() + (delta_max[i] * 1e9).tolist())
    with (OUT / "native_vector_sign_history_201.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["point_index", "solver_step", "command_nm"] +
                        [name for c in range(opening.shape[1])
                         for name in (f"opening_z_cell_{c}_m", f"traction_z_cell_{c}_Pa",
                                      f"traction_dot_opening_cell_{c}_Pa_m")])
        for i in range(n):
            writer.writerow([int(pidx[i]), int(solver_step[i]), float(s[i] * 1e9)] +
                            [value for c in range(opening.shape[1])
                             for value in (float(opening[i, c, 2]), float(traction[i, c, 2]),
                                           float(dot[i, c]))])
    interfaces = np.array_equal(initial_rest[initial_coh[:, :3]], initial_rest[initial_coh[:, 3:]])
    topology_stable = bool(np.all(counts == np.array([54, 96, 8])) and
                           np.all(counts == counts[0]) and np.array_equal(coh, initial_coh) and
                           np.array_equal(initial_rest, rest))
    mass_gate = bool(np.all(mass_sums > 0) and np.max(np.abs(mass_sums - 4.0e-5)) <= 1e-12)
    damage_final_near_one = n == 201 and abs(float(damage[-1].mean()) - 1.0) <= 0.01
    delta_final_near_max_command = n == 201 and abs(float(delta_max[-1].max()) - float(s.max())) <= 1e-12
    checks = {
        "exact_201_points_and_solver_steps": n == 201 and status.get("solver_steps_completed") == 201 and
                                             np.array_equal(pidx, np.arange(201)) and np.array_equal(solver_step, np.arange(1, 202)),
        "distinct_201_native_dump_records": n == 201 and status.get("dump_frames_and_parsed_records") == 201,
        "stable_54_96_8_topology_and_rest": topology_stable and interfaces,
        "positive_mass_each_dof_4e-5kg": mass_gate,
        "opening_matches_command_abs_1e-12m": n == 201 and float(opening_err.max()) <= 1e-12,
        "prescribed_boundary_readback_abs_1e-12m": n == 201 and float(boundary_error.max()) <= 1e-12,
        "traction_profile_within_1pct_sigma": n == 201 and float(traction_err.max()) <= 0.01 * SIGMA,
        "initial_and_final_traction_near_zero": n == 201 and max(float(traction_norm[0].max()), float(traction_norm[-1].max())) <= 0.01 * SIGMA,
        "dissipated_energy_not_below_minus_1e-12J": bool(np.isfinite(d).all() and np.min(d) >= -1e-12),
        "final_dissipation_matches_GcA_1pct": n == 201 and gc_area > 0 and abs(float(d[-1]) - gc_area) / gc_area <= 0.01,
        "final_reversible_energy_abs_1e-12J": n == 201 and abs(float(r[-1])) <= 1e-12,
        "damage_finite_0_to_1": damage_bounds,
        "damage_non_decreasing_loading_tol_1e-12": damage_monotonic,
        "damage_final_within_0.01_of_one": damage_final_near_one,
        "delta_max_finite_and_loading_non_decreasing_tol_1e-12m": finite_delta and delta_monotonic,
        "delta_max_final_within_1e-12m_of_max_command": delta_final_near_max_command,
        "eight_cell_opening_uniformity_1e-12m": float(opening_spread.max()) <= 1e-12,
        "eight_cell_traction_norm_uniformity_1pct_sigma": float(traction_spread.max()) <= 0.01 * SIGMA,
        "eight_cell_damage_uniformity_1e-12": float(damage_cell_spreads.max()) <= 1e-12,
        "eight_cell_delta_max_uniformity_1e-12m": float(delta_cell_spreads.max()) <= 1e-12,
    }
    full = n == 201 and status.get("status") == "RUN_201_POINTS_COMPLETE"
    result = {
        "case": "controlled_intrinsic_bilinear_modeI",
        "run_status": status.get("status"),
        "solver_steps_completed": status.get("solver_steps_completed"),
        "native_dump_records_parsed": n,
        "analysis_status": "COMPLETE" if full else "PARTIAL",
        "gates": checks,
        "all_gates_pass": bool(full and all(checks.values())),
        "scope": "UNVALIDATED toy constitutive patch only; full prescribed DOFs, not free dynamics. No real-material calibration, physical fracture/fragment claim, mesh-bias result, R2, or PRETRAIN_ELIGIBLE.",
        "representation": "Native cohesive dump contains one average tuple per cohesive cell (8 tuples/frame), not 24 quadrature-point samples.",
        "units_note": "VTK does not encode physical units; interpreted by field semantics.",
        "thresholds": {"opening_abs_m": 1e-12, "traction_abs_Pa": 0.01 * SIGMA,
                       "damage_monotonic_absolute": 1e-12, "delta_max_monotonic_m": 1e-12,
                       "damage_bounds": [0.0, 1.0], "final_damage_abs_from_one": 0.01,
                       "dissipated_lower_J": -1e-12, "final_D_vs_GcA_relative": 0.01,
                       "final_reversible_abs_J": 1e-12},
        "areas_m2": {"per_native_cell": area_cell.tolist(), "total": area, "GcA_J": gc_area},
        "native_energies_J": {
            "kinetic_initial": float(kinetic[0]), "kinetic_final": float(kinetic[-1]),
            "potential_initial": float(potential[0]), "potential_final": float(potential[-1]),
            "dissipated_initial": float(d[0]), "dissipated_min": float(d.min()),
            "dissipated_final": float(d[-1]), "reversible_initial": float(r[0]),
            "reversible_final": float(r[-1]), "cohesive_contact_initial": float(contact[0]),
            "cohesive_contact_final": float(contact[-1]), "GcA_J": gc_area},
        "native_damage": {
            "initial_cell_values": damage[0].tolist(), "initial_min_mean_max": [float(damage[0].min()), float(damage[0].mean()), float(damage[0].max())],
            "peak_point_index": damage_peak_i, "peak_command_nm": float(s[damage_peak_i] * 1e9),
            "peak_cell_values": damage[damage_peak_i].tolist(),
            "peak_min_mean_max": [float(damage[damage_peak_i].min()), float(damage_means[damage_peak_i]), float(damage[damage_peak_i].max())],
            "final_cell_values": damage[-1].tolist(), "final_min_mean_max": [float(damage[-1].min()), float(damage[-1].mean()), float(damage[-1].max())],
            "max_cell_spread": float(damage_cell_spreads.max()), "bounds_pass": damage_bounds,
            "monotonicity_tolerance": -1e-12, "full_history_csv": "damage_delta_history_201.csv"},
        "native_delta_max_m": {
            "initial_cell_values": delta_max[0].tolist(), "initial_min_mean_max": [float(delta_max[0].min()), float(delta_max[0].mean()), float(delta_max[0].max())],
            "peak_point_index": delta_peak_i, "peak_command_nm": float(s[delta_peak_i] * 1e9),
            "peak_cell_values": delta_max[delta_peak_i].tolist(),
            "final_cell_values": delta_max[-1].tolist(), "final_min_mean_max": [float(delta_max[-1].min()), float(delta_max[-1].mean()), float(delta_max[-1].max())],
            "max_cell_spread_m": float(delta_cell_spreads.max()), "final_vs_max_command_abs_m": abs(float(delta_max[-1].max()) - float(s.max())),
            "monotonicity_tolerance_m": -1e-12, "full_history_csv": "damage_delta_history_201.csv"},
        "opening_traction_diagnostics": {
            "max_opening_norm_command_error_m": float(opening_err.max()),
            "max_traction_norm_profile_error_Pa": float(traction_err.max()),
            "max_traction_norm_profile_error_fraction_sigma": float(traction_err.max() / SIGMA),
            "first_cell_traction_norms_Pa": traction_norm[0].tolist(),
            "peak_expected_traction_Pa": float(expected.max()),
            "final_cell_traction_norms_Pa": traction_norm[-1].tolist(),
            "first_opening_z_m": opening[0, :, 2].tolist(),
            "first_traction_z_Pa": traction[0, :, 2].tolist(),
            "final_opening_z_m": opening[-1, :, 2].tolist(),
            "final_traction_z_Pa": traction[-1, :, 2].tolist(),
            "cell_opening_dot_traction_initial_Pa_m": dot[0].tolist(),
            "cell_opening_dot_traction_final_Pa_m": dot[-1].tolist(),
            "max_cell_opening_spread_m": float(opening_spread.max()),
            "max_cell_traction_norm_spread_Pa": float(traction_spread.max())},
        "area_traction_opening_work_signed_J": signed_work,
        "upper_internal_force_z_sum_raw_N": {"initial": float(force_z[0]), "final": float(force_z[-1])},
        "full_history_files": {"npz": "native_history_partial_or_complete.npz",
                               "per_point_npz_directory": "samples", "trace": "trace.jsonl",
                               "native_vtk_directory": "paraview",
                               "vector_sign_history_csv": "native_vector_sign_history_201.csv"},
        "energy_balance": "NOT_APPLICABLE as a pass criterion: all DOFs are prescribed. Kinetic and work values are diagnostics only.",
        "plot_png": "bilinear_patch_review.png" if full else "partial_bilinear_patch_review.png"}
    write_json(OUT / "analysis_201.json", result)
    print(json.dumps({"analysis_status": result["analysis_status"], "all_gates_pass": result["all_gates_pass"],
                      "steps": n, "failed_gates": [k for k, v in checks.items() if not v]}, allow_nan=False))
    return result


if __name__ == "__main__":
    analyze()
