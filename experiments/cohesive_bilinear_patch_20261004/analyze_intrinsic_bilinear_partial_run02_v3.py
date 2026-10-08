from __future__ import annotations

import csv
import json
import re
import hashlib
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

RAW = Path(r"C:\Users\xiemo\Desktop\ShardForge\relicweave-synlab\tmp\review_20261004\implementation\cohesive_bilinear_patch_20261004\controlled_intrinsic_bilinear_modeI_201_run02")
EXP = Path(r"C:\Users\xiemo\Desktop\ShardForge\pgfrag\experiments\cohesive_bilinear_patch_20261004")
MANIFEST = EXP / "manifest_bilinear_modeI_201_v2.json"
DERIVED = RAW / "derived_partial_v3"
ENERGY_KEYS = ("kinetic", "potential", "dissipated", "reversible", "cohesive contact")
REQUIRED_SAMPLE_KEYS = {
    "point_index", "solver_step", "commanded_opening_m", "opening", "tractions",
    "damage", "delta_max", "displacement", "internal_force", "mass", "rest_nodes",
    "bulk_tets", "cohesive_connectivity", "energies_J", "mass_component_sums_kg",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def analytic_traction(s: float, sigma: float, gc: float, delta0: float) -> float:
    deltac = 2.0 * gc / sigma
    if s <= delta0:
        return sigma * s / delta0
    if s <= deltac:
        return sigma * (deltac - s) / (deltac - delta0)
    return 0.0


def index_from_name(path: Path):
    match = re.search(r"_(\d{4})\.(?:pvtu|vtu)$", path.name)
    return int(match.group(1)) if match else None


def main():
    if DERIVED.exists() and any(DERIVED.iterdir()):
        raise SystemExit(f"derived output already exists; refusing overwrite: {DERIVED}")
    DERIVED.mkdir(parents=True, exist_ok=True)
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    sigma = float(manifest["material"]["cohesive"]["sigma_c_Pa"])
    gc = float(manifest["material"]["cohesive"]["G_c_J_m2"])
    delta0 = float(manifest["material"]["cohesive"]["delta_0_m"])
    deltac = float(manifest["material"]["cohesive"]["delta_c_theory_m"])
    run_started = json.loads((RAW / "run_started.json").read_text(encoding="utf-8"))
    external = json.loads((RAW / "external_execution_receipt_v2.json").read_text(encoding="utf-8"))
    process_check = json.loads((RAW / "post_timeout_process_check_v2.json").read_text(encoding="utf-8"))
    run_status_present = (RAW / "run_status.json").is_file()

    trace_path = RAW / "trace.jsonl"
    trace_rows = [json.loads(line) for line in trace_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    trace_ids = [int(row["point_index"]) for row in trace_rows]
    trace_by_id = {int(row["point_index"]): row for row in trace_rows}
    trace_duplicate_ids = sorted({i for i in trace_ids if trace_ids.count(i) > 1})

    sample_dir = RAW / "samples"
    sample_files = sorted(sample_dir.glob("point_*.npz"))
    sample_ids_by_file = {}
    for path in sample_files:
        match = re.fullmatch(r"point_(\d{3})\.npz", path.name)
        if match:
            sample_ids_by_file[int(match.group(1))] = path
    sample_file_ids = sorted(sample_ids_by_file)

    pvtu_files = sorted((RAW / "paraview").rglob("*.pvtu"))
    vtu_files = sorted((RAW / "paraview").rglob("*.vtu"))
    pvtu_ids = sorted(i for p in pvtu_files if (i := index_from_name(p)) is not None)
    vtu_ids = sorted(i for p in vtu_files if (i := index_from_name(p)) is not None)

    initial_path = RAW / "initial_state.npz"
    initial_data = np.load(initial_path, allow_pickle=False)
    initial = {k: initial_data[k].copy() for k in initial_data.files}
    initial_data.close()
    reference_rest = initial["rest_nodes"]
    reference_tets = initial["bulk_tets"]
    reference_coh = initial["cohesive_connectivity"]
    initial_blocked = initial.get("blocked_dofs")
    missing_initial = sorted({"rest_nodes", "bulk_tets", "cohesive_connectivity"} - set(initial))

    required_pairs = set(trace_ids) & set(sample_file_ids) & set(pvtu_ids) & set(vtu_ids)
    usable = []
    invalid_samples = []
    for index in sorted(required_pairs):
        path = sample_ids_by_file[index]
        try:
            with np.load(path, allow_pickle=False) as z:
                keys = set(z.files)
                absent = sorted(REQUIRED_SAMPLE_KEYS - keys)
                if absent:
                    raise ValueError(f"missing keys {absent}")
                row = {key: z[key].copy() for key in z.files}
            if int(row["point_index"]) != index:
                raise ValueError("sample point_index differs from filename")
            if int(row["solver_step"]) != index + 1:
                raise ValueError("solver_step is not point_index+1")
            tr = trace_by_id[index]
            if int(tr["solver_step"]) != int(row["solver_step"]):
                raise ValueError("trace/sample solver_step mismatch")
            if int(tr["dump_index"]) != index or int(tr["native_dump"]["step_index"]) != index:
                raise ValueError("trace native dump index mismatch")
            if float(tr["commanded_opening_m"]) != float(row["commanded_opening_m"]):
                raise ValueError("trace/sample command mismatch")
            for key in ("opening", "tractions", "damage", "delta_max", "displacement", "internal_force",
                        "mass", "rest_nodes", "bulk_tets", "cohesive_connectivity", "energies_J",
                        "mass_component_sums_kg"):
                if not np.isfinite(row[key]).all():
                    raise ValueError(f"nonfinite sample array {key}")
            usable.append((index, row, tr))
        except Exception as exc:
            invalid_samples.append({"point_index": index, "path": str(path), "error": repr(exc)})

    usable_ids = [item[0] for item in usable]
    n = len(usable)
    if n == 0:
        raise RuntimeError("no closed sample NPZs could be paired to trace and native VTK frames")
    if missing_initial:
        raise RuntimeError(f"initial state lacks arrays required for exact topology comparison: {missing_initial}")

    pidx = np.asarray(usable_ids, dtype=np.int64)
    solver_steps = np.asarray([int(row["solver_step"]) for _, row, _ in usable], dtype=np.int64)
    command = np.asarray([float(row["commanded_opening_m"]) for _, row, _ in usable], dtype=np.float64)
    opening = np.stack([row["opening"] for _, row, _ in usable])
    traction = np.stack([row["tractions"] for _, row, _ in usable])
    damage_raw = np.stack([row["damage"] for _, row, _ in usable])
    delta_raw = np.stack([row["delta_max"] for _, row, _ in usable])
    if damage_raw.shape != (n, 8, 1) or delta_raw.shape != (n, 8, 1):
        raise RuntimeError(f"native cell-average scalar shapes differ from (n,8,1): damage={damage_raw.shape}, delta_max={delta_raw.shape}")
    damage = damage_raw[:, :, 0].copy()
    delta_max = delta_raw[:, :, 0].copy()
    energies = np.stack([row["energies_J"] for _, row, _ in usable])
    mass_sums = np.stack([row["mass_component_sums_kg"] for _, row, _ in usable])
    counts = np.asarray([[int(tr["counts"][k]) for k in ("nodes", "bulk_tets", "cohesive_3d6")]
                         for _, _, tr in usable], dtype=np.int64)
    post_error = np.asarray([float(tr["post_solve_boundary_readback_error_m"]) for _, _, tr in usable])
    pre_error = np.asarray([float(tr["prescribed_readback_error_m"]) for _, _, tr in usable])
    blocked_readback = np.asarray([tr.get("native_blocked_dofs_all_true_pre_solve") is True
                                   for _, _, tr in usable], dtype=bool)
    force_z = np.asarray([float(tr["upper_internal_force_z_sum_raw_N"]) for _, _, tr in usable])
    rest_hist = np.stack([row["rest_nodes"] for _, row, _ in usable])
    tet_hist = np.stack([row["bulk_tets"] for _, row, _ in usable])
    coh_hist = np.stack([row["cohesive_connectivity"] for _, row, _ in usable])
    mass_hist = np.stack([row["mass"] for _, row, _ in usable])
    disp_hist = np.stack([row["displacement"] for _, row, _ in usable])
    force_hist = np.stack([row["internal_force"] for _, row, _ in usable])

    sample_index_contiguous = np.array_equal(pidx, np.arange(n))
    trace_exact_ids = trace_ids == usable_ids
    sample_vtk_exact_ids = sample_file_ids == pvtu_ids == vtu_ids == usable_ids
    counts_stable = bool(np.all(counts == np.asarray([54, 96, 8])) and np.all(counts == counts[0]))
    exact_rest = bool(np.all([np.array_equal(item, reference_rest) for item in rest_hist]))
    exact_tets = bool(np.all([np.array_equal(item, reference_tets) for item in tet_hist]))
    exact_coh = bool(np.all([np.array_equal(item, reference_coh) for item in coh_hist]))
    duplicate_face_pairing_exact = bool(np.array_equal(
        reference_rest[reference_coh[:, :3]], reference_rest[reference_coh[:, 3:]]))
    per_dof_mass_positive = bool(np.all(np.isfinite(mass_hist)) and np.all(mass_hist > 0.0))
    mass_sum_pass = bool(np.all(mass_sums > 0.0) and np.max(np.abs(mass_sums - 4.0e-5)) <= 1.0e-12)

    opening_norm = np.linalg.norm(opening, axis=2)
    traction_norm = np.linalg.norm(traction, axis=2)
    expected_norm = np.asarray([analytic_traction(float(s), sigma, gc, delta0) for s in command])
    opening_abs_error = np.abs(opening_norm - command[:, None])
    traction_profile_abs_error = np.abs(traction_norm - expected_norm[:, None])
    dot_signed = np.einsum("tci,tci->tc", opening, traction)
    opening_z = opening[:, :, 2]
    traction_z = traction[:, :, 2]
    min_damage = float(damage.min())
    max_damage = float(damage.max())
    damage_finite = bool(np.isfinite(damage_raw).all())
    damage_bounds = damage_finite and bool(np.all((damage >= 0.0) & (damage <= 1.0)))
    damage_monotone = damage_finite and bool(np.all(np.diff(damage, axis=0) >= -1.0e-12))
    delta_finite = bool(np.isfinite(delta_raw).all())
    delta_monotone = delta_finite and bool(np.all(np.diff(delta_max, axis=0) >= -1.0e-12))
    damage_spread = np.ptp(damage, axis=1)
    delta_spread = np.ptp(delta_max, axis=1)
    opening_norm_spread = np.ptp(opening_norm, axis=1)
    traction_norm_spread = np.ptp(traction_norm, axis=1)
    energies_finite = bool(np.isfinite(energies).all())
    dissipated = energies[:, 2]
    reversible = energies[:, 3]
    contact = energies[:, 4]
    negative_indices = [int(pidx[i]) for i in range(n) if dissipated[i] < 0.0]
    min_d_i = int(np.argmin(dissipated))
    area_cell = np.asarray([
        0.5 * np.linalg.norm(np.cross(reference_rest[row[1]] - reference_rest[row[0]],
                                       reference_rest[row[2]] - reference_rest[row[0]]))
        for row in reference_coh
    ])
    total_area = float(area_cell.sum())
    gc_area = gc * total_area
    signed_work = 0.0
    for i in range(1, n):
        delta_vec = opening[i] - opening[i - 1]
        mean_traction = 0.5 * (traction[i] + traction[i - 1])
        signed_work += float(np.sum(area_cell * np.einsum("ci,ci->c", mean_traction, delta_vec)))

    full_native_run = bool(run_status_present and n == 201 and
                           external.get("status") == "PYTHON_EXIT_0" and
                           not external.get("external_timeout_triggered", False))
    reached_deltac = bool(float(command.max()) >= deltac)
    gates = {
        "168_available_points_contiguous_from_zero": {
            "status": "PASS" if n == 168 and sample_index_contiguous and trace_exact_ids and sample_vtk_exact_ids else "FAIL",
            "observed_points": n, "indices_first_last": [int(pidx[0]), int(pidx[-1])]},
        "full_201_native_points_and_solver_steps": {
            "status": "PASS" if full_native_run and n == 201 and np.array_equal(pidx, np.arange(201)) and np.array_equal(solver_steps, np.arange(1, 202)) else "FAIL_INCOMPLETE",
            "observed_native_samples": n, "required": 201},
        "stable_exact_54_nodes_96_tets_8_cohesive_counts": {
            "status": "PASS_PARTIAL" if counts_stable else "FAIL",
            "counts_unique": np.unique(counts, axis=0).tolist()},
        "rest_nodes_tets_cohesive_connectivity_exact_to_initial": {
            "status": "PASS_PARTIAL" if exact_rest and exact_tets and exact_coh else "FAIL",
            "exact_rest": exact_rest, "exact_tets": exact_tets, "exact_cohesive_connectivity": exact_coh},
        "initial_cohesive_face_pairing_exact_rest_coordinates": {
            "status": "PASS" if duplicate_face_pairing_exact else "FAIL"},
        "positive_nodal_mass_and_each_component_sum_4e-5kg": {
            "status": "PASS_PARTIAL" if per_dof_mass_positive and mass_sum_pass else "FAIL",
            "all_nodal_mass_positive": per_dof_mass_positive,
            "component_sums_kg_first_last": [mass_sums[0].tolist(), mass_sums[-1].tolist()]},
        "native_blocked_dofs_all_true_and_prescribed_boundary_readback": {
            "status": "PASS_PARTIAL" if bool(blocked_readback.all()) and float(max(pre_error.max(), post_error.max())) <= 1.0e-12 else "FAIL",
            "all_step_fresh_blocked_readbacks_true": bool(blocked_readback.all()),
            "max_pre_error_m": float(pre_error.max()), "max_post_error_m": float(post_error.max())},
        "native_opening_norm_matches_command_abs_1e-12m_partial": {
            "status": "PASS_PARTIAL" if float(opening_abs_error.max()) <= 1.0e-12 else "FAIL",
            "max_abs_error_m": float(opening_abs_error.max())},
        "native_traction_norm_vs_analytic_profile_within_1pct_sigma_partial": {
            "status": "PASS_PARTIAL" if float(traction_profile_abs_error.max()) <= 0.01 * sigma else "FAIL",
            "max_abs_error_Pa": float(traction_profile_abs_error.max()),
            "threshold_Pa": 0.01 * sigma},
        "dissipated_native_energy_not_below_minus_1e-12J_partial": {
            "status": "PASS_PARTIAL" if energies_finite and float(dissipated.min()) >= -1.0e-12 else "FAIL",
            "min_raw_dissipated_J": float(dissipated.min()),
            "min_point_index": int(pidx[min_d_i]),
            "negative_value_point_indices_raw_sign_retained": negative_indices},
        "native_damage_finite_0_to_1_and_loading_monotone_partial": {
            "status": "PASS_PARTIAL" if damage_bounds and damage_monotone else "FAIL",
            "finite": damage_finite, "minimum": min_damage, "maximum": max_damage,
            "monotone_tolerance": -1.0e-12},
        "native_delta_max_finite_and_monotone_partial": {
            "status": "PASS_PARTIAL" if delta_finite and delta_monotone else "FAIL",
            "finite": delta_finite, "monotone_tolerance_m": -1.0e-12},
        "final_dissipation_matches_GcA_1pct": {
            "status": "NOT_EVALUATED_INCOMPLETE" if not full_native_run or not reached_deltac else "CHECK_REQUIRED"},
        "final_reversible_energy_abs_1e-12J": {
            "status": "NOT_EVALUATED_INCOMPLETE" if not full_native_run else "CHECK_REQUIRED"},
        "final_damage_near_one_and_delta_max_near_command": {
            "status": "NOT_EVALUATED_INCOMPLETE" if not full_native_run else "CHECK_REQUIRED"},
        "initial_and_final_traction_near_zero_full_curve": {
            "status": "NOT_EVALUATED_INCOMPLETE" if not full_native_run else "CHECK_REQUIRED"},
    }

    with (DERIVED / "partial_native_energy_history_v2.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["point_index", "solver_step", "commanded_opening_m", "kinetic_J", "potential_J",
                         "dissipated_J_raw_signed", "reversible_J_raw_signed", "cohesive_contact_J_raw_signed",
                         "upper_internal_force_z_sum_raw_N"])
        for i, (_, _, tr) in enumerate(usable):
            writer.writerow([int(pidx[i]), int(solver_steps[i]), float(command[i]),
                             *[float(x) for x in energies[i]], float(force_z[i])])

    with (DERIVED / "partial_native_vector_sign_history_v2.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["point_index", "solver_step", "cell", "commanded_opening_m",
                         "opening_x_m_raw", "opening_y_m_raw", "opening_z_m_raw",
                         "traction_x_Pa_raw", "traction_y_Pa_raw", "traction_z_Pa_raw",
                         "traction_dot_opening_Pa_m_raw_signed"])
        for i in range(n):
            for cell in range(8):
                writer.writerow([int(pidx[i]), int(solver_steps[i]), cell, float(command[i]),
                                 *opening[i, cell].tolist(), *traction[i, cell].tolist(), float(dot_signed[i, cell])])

    with (DERIVED / "partial_native_damage_delta_history_v2.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["point_index", "solver_step", "commanded_opening_m", "cell",
                         "damage_native_cell_average", "delta_max_m_native_cell_average"])
        for i in range(n):
            for cell in range(8):
                writer.writerow([int(pidx[i]), int(solver_steps[i]), float(command[i]), cell,
                                 float(damage[i, cell]), float(delta_max[i, cell])])

    derived_npz = DERIVED / "native_partial_history_derived_from_closed_samples_v2.npz"
    np.savez_compressed(
        derived_npz,
        derived_artifact_kind=np.asarray("DERIVED_PARTIAL_FROM_CLOSED_POINT_NPZ_NOT_NATIVE_RUN_OUTPUT"),
        source_run_id=np.asarray("controlled_intrinsic_bilinear_modeI_201_run02"),
        point_index=pidx, solver_step=solver_steps, commanded_opening_m=command,
        opening=opening, tractions=traction, opening_z_m=opening_z, traction_z_Pa=traction_z,
        traction_dot_opening_Pa_m_signed=dot_signed, damage_cell_average=damage,
        delta_max_cell_average_m=delta_max, energies_J=energies, energy_keys=np.asarray(ENERGY_KEYS),
        mass_component_sums_kg=mass_sums, counts=counts, rest_nodes=rest_hist,
        bulk_tets=tet_hist, cohesive_connectivity=coh_hist, mass=mass_hist,
        displacement=disp_hist, internal_force=force_hist)

    # Geometry is shown from unchanged rest coordinates only; no displacement amplification.
    face_a = reference_rest[reference_coh[:, :3], :2] * 1000.0
    face_b = reference_rest[reference_coh[:, 3:], :2] * 1000.0
    fig, axs = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    ax = axs[0, 0]
    for cell in range(8):
        ax.plot(opening_z[:, cell] * 1.0e9, traction_z[:, cell] / 1.0e6, color="steelblue", alpha=0.45, linewidth=0.9)
    ax.plot(opening_z.mean(axis=1) * 1.0e9, traction_z.mean(axis=1) / 1.0e6,
            color="navy", marker=".", linewidth=1.8, label="native signed z (cell mean)")
    xref_nm = np.linspace(0.0, max(float(command.max() * 1.0e9), deltac * 1.0e9), 300)
    yref = np.asarray([analytic_traction(float(x * 1.0e-9), sigma, gc, delta0) for x in xref_nm]) / 1.0e6
    ax.plot(xref_nm, yref, color="darkorange", linestyle="--", linewidth=1.7, label="analytic scalar reference")
    ax.set_xlabel("native opening_z (nm)")
    ax.set_ylabel("native traction_z (MPa)")
    ax.set_title("Signed native z response; no sign flips")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    ax.text(0.02, 0.02, f"observed through {command[-1]*1e9:.3f} nm; delta_c={deltac*1e9:.3f} nm",
            transform=ax.transAxes, fontsize=8)

    ax = axs[0, 1]
    x = command * 1.0e9
    ax.plot(x, energies[:, 2] * 1.0e6, label="dissipated raw signed")
    ax.plot(x, energies[:, 3] * 1.0e6, label="reversible raw signed")
    ax.plot(x, energies[:, 4] * 1.0e6, label="cohesive contact raw")
    ax.plot(x, energies[:, 1] * 1.0e6, label="bulk potential raw")
    ax.plot(x, energies[:, 0] * 1.0e6, label="kinetic raw")
    ax.axhline(gc_area * 1.0e6, color="black", linestyle=":", label="Gc × rest area reference")
    ax.set_xlabel("commanded opening (nm)")
    ax.set_ylabel("native energy (µJ)")
    ax.set_title("Native energies; incomplete loading path")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=7)

    for ax, faces, color, side in ((axs[1, 0], face_a, "seagreen", "A"),
                                   (axs[1, 1], face_b, "darkorange", "B")):
        polys = [face.tolist() for face in faces]
        collection = PolyCollection(polys, facecolors=color, edgecolors="black", alpha=0.55, linewidths=0.7)
        ax.add_collection(collection)
        all_xy = faces.reshape(-1, 2)
        ax.scatter(all_xy[:, 0], all_xy[:, 1], s=5, color="black")
        ax.set_xlim(float(all_xy[:, 0].min()) - 0.1, float(all_xy[:, 0].max()) + 0.1)
        ax.set_ylim(float(all_xy[:, 1].min()) - 0.1, float(all_xy[:, 1].max()) + 0.1)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        ax.set_title(f"Native rest interface side {side}; actual scale")
        ax.grid(True, alpha=0.2)
    fig.suptitle(f"run02 PARTIAL native evidence: {n}/201 points; external timeout, no gap exaggeration")
    png_path = DERIVED / "partial_bilinear_patch_review_v2.png"
    fig.savefig(png_path, dpi=160)
    plt.close(fig)

    first_idx = 0
    last_idx = n - 1
    peak_traction_index = int(np.argmax(traction_norm.mean(axis=1)))
    first_negative = next((int(pidx[i]) for i in range(n) if dissipated[i] < 0.0), None)
    sentinel_indices = sorted(set([0, min(20, n - 1), n - 1]))
    sentinel = {}
    for index in sentinel_indices:
        i = usable_ids.index(index)
        sentinel[str(index)] = {
            "commanded_opening_m": float(command[i]),
            "opening_z_m_raw_per_cell": opening_z[i].tolist(),
            "traction_z_Pa_raw_per_cell": traction_z[i].tolist(),
            "traction_dot_opening_Pa_m_raw_signed_per_cell": dot_signed[i].tolist(),
            "damage_raw_cell_average_per_cell": damage[i].tolist(),
            "delta_max_m_raw_cell_average_per_cell": delta_max[i].tolist(),
            "native_energies_J_raw_signed": {key: float(energies[i, j]) for j, key in enumerate(ENERGY_KEYS)},
        }

    summary = {
        "analysis_status": "PARTIAL_DERIVED_ONLY",
        "candidate_constitutive_check": "NOT_EVALUATED_INCOMPLETE",
        "native_run": {
            "run_id": manifest["run_id"],
            "runner_status_file_present": run_status_present,
            "run_started_solver_attempts": run_started.get("solver_attempts"),
            "external_status": external.get("status"),
            "external_timeout_triggered": external.get("external_timeout_triggered"),
            "external_elapsed_wall_s": external.get("elapsed_wall_s"),
            "python_exit_code": external.get("python_exit_code"),
            "wsl_host_exit_code": external.get("wsl_host_exit_code"),
            "powershell_command_exit_code": external.get("powershell_command_exit_code"),
            "target_process_check": process_check.get("status"),
            "trace_solver_steps_evidenced": int(len(trace_rows)),
            "closed_sample_npz_count": int(len(sample_files)),
            "valid_closed_sample_npz_count": n,
            "closed_pvtu_count": int(len(pvtu_files)),
            "closed_vtu_count": int(len(vtu_files)),
            "trace_indices": trace_ids,
            "sample_indices": sample_file_ids,
            "pvtu_indices": pvtu_ids,
            "vtu_indices": vtu_ids,
            "sample_invalid": invalid_samples,
            "trace_duplicate_indices": trace_duplicate_ids,
            "trace_sample_pvtu_vtu_indices_exactly_equal": bool(trace_exact_ids and sample_vtk_exact_ids),
            "solver_steps_completed_from_fsynced_trace_samples": int(solver_steps[-1]),
            "run_complete": full_native_run,
            "last_commanded_opening_m": float(command[-1]),
            "last_commanded_opening_nm": float(command[-1] * 1.0e9),
            "delta_c_theory_m": deltac,
            "delta_c_theory_nm": float(deltac * 1.0e9),
            "delta_c_reached": reached_deltac,
            "per_step_solver_wall_time": "UNKNOWN; trace records do not carry per-step wall time",
            "cause_of_300s_timeout": "UNKNOWN; elapsed wall limit is observed, no evidence attributes it specifically to filesystem I/O",
        },
        "raw_validation": {
            "indices_contiguous_from_zero": sample_index_contiguous,
            "all_trace_sample_native_frame_indices_match": bool(trace_exact_ids and sample_vtk_exact_ids),
            "all_sample_npz_closed_and_readable": len(invalid_samples) == 0 and len(sample_files) == n,
            "raw_sample_fields_finite": bool(all(np.isfinite(a).all() for a in
                (opening, traction, damage_raw, delta_raw, energies, mass_hist, rest_hist, tet_hist, coh_hist,
                 disp_hist, force_hist))),
            "native_scalar_representation": "damage/delta max are 8 per-cohesive-cell averages with one scalar tuple/cell, not QP data",
            "vtk_units": "Native VTK has no physical-unit metadata; m/Pa/J interpretation is semantic inference from field names, manifest and runner contract.",
        },
        "gates": gates,
        "counts_topology_mass_bc": {
            "counts_per_frame_unique": np.unique(counts, axis=0).tolist(),
            "rest_nodes_exact_to_initial_all_frames": exact_rest,
            "bulk_tets_exact_to_initial_all_frames": exact_tets,
            "cohesive_connectivity_exact_to_initial_all_frames": exact_coh,
            "initial_interface_face_pair_coordinates_exact": duplicate_face_pairing_exact,
            "initial_blocked_dofs_shape": list(initial_blocked.shape) if initial_blocked is not None else None,
            "initial_blocked_dofs_all_true": bool(np.all(initial_blocked)) if initial_blocked is not None else None,
            "per_step_native_blocked_readback_all_true": bool(blocked_readback.all()),
            "mass_component_sums_kg_first": mass_sums[0].tolist(),
            "mass_component_sums_kg_last": mass_sums[-1].tolist(),
            "all_nodal_mass_positive": per_dof_mass_positive,
            "maximum_pre_solve_boundary_readback_error_m": float(pre_error.max()),
            "maximum_post_solve_boundary_readback_error_m": float(post_error.max()),
        },
        "native_energy_partial": {
            "energy_key_order": list(ENERGY_KEYS),
            "initial_raw_signed_J": {key: float(energies[0, j]) for j, key in enumerate(ENERGY_KEYS)},
            "last_raw_signed_J": {key: float(energies[-1, j]) for j, key in enumerate(ENERGY_KEYS)},
            "dissipated_min_raw_signed_J": float(dissipated.min()),
            "dissipated_min_point_index": int(pidx[min_d_i]),
            "first_negative_dissipated_point_index": first_negative,
            "dissipated_negative_points_raw_sign_retained": negative_indices,
            "Gc_times_rest_area_reference_J": gc_area,
            "final_D_vs_GcA_gate": "NOT_EVALUATED_INCOMPLETE; delta_c not reached",
            "final_reversible_gate": "NOT_EVALUATED_INCOMPLETE; 201-point endpoint not reached",
        },
        "native_damage_delta_partial": {
            "damage_initial_per_cell": damage[0].tolist(),
            "damage_partial_max_per_cell": damage.max(axis=0).tolist(),
            "damage_last_per_cell": damage[-1].tolist(),
            "damage_min_max": [min_damage, max_damage],
            "damage_bounds_0_to_1_partial": damage_bounds,
            "damage_monotone_partial_tolerance": -1.0e-12,
            "damage_monotone_partial": damage_monotone,
            "delta_max_initial_m_per_cell": delta_max[0].tolist(),
            "delta_max_last_m_per_cell": delta_max[-1].tolist(),
            "delta_max_last_nm_per_cell": (delta_max[-1] * 1.0e9).tolist(),
            "delta_max_monotone_partial_tolerance_m": -1.0e-12,
            "delta_max_monotone_partial": delta_monotone,
            "final_damage_near_one_gate": "NOT_EVALUATED_INCOMPLETE; delta_c/post-peak endpoint not reached",
        },
        "opening_traction_signed_partial": {
            "max_opening_norm_vs_command_error_m": float(opening_abs_error.max()),
            "max_traction_norm_vs_analytic_profile_error_Pa": float(traction_profile_abs_error.max()),
            "partial_traction_norm_profile_relative_to_sigma": float(traction_profile_abs_error.max() / sigma),
            "peak_observed_traction_norm_point_index": int(pidx[peak_traction_index]),
            "peak_observed_traction_norm_MPa": float(traction_norm[peak_traction_index].mean() / 1.0e6),
            "first_opening_z_m_raw_per_cell": opening_z[0].tolist(),
            "first_traction_z_Pa_raw_per_cell": traction_z[0].tolist(),
            "last_opening_z_m_raw_per_cell": opening_z[-1].tolist(),
            "last_traction_z_Pa_raw_per_cell": traction_z[-1].tolist(),
            "last_opening_z_nm_raw_per_cell": (opening_z[-1] * 1.0e9).tolist(),
            "last_traction_z_MPa_raw_per_cell": (traction_z[-1] / 1.0e6).tolist(),
            "first_T_dot_opening_Pa_m_raw_signed_per_cell": dot_signed[0].tolist(),
            "last_T_dot_opening_Pa_m_raw_signed_per_cell": dot_signed[-1].tolist(),
            "area_traction_opening_work_signed_diagnostic_J": signed_work,
            "work_note": "Independent signed trapezoidal area×traction·opening diagnostic only; not the native dissipated/reversible energy.",
            "sampled_sentinels": sentinel,
        },
        "derived_outputs": {
            "partial_npz": str(derived_npz),
            "energy_csv": str(DERIVED / "partial_native_energy_history_v2.csv"),
            "signed_vector_csv": str(DERIVED / "partial_native_vector_sign_history_v2.csv"),
            "damage_delta_csv": str(DERIVED / "partial_native_damage_delta_history_v2.csv"),
            "png": str(png_path),
            "source_analyzer_sha256": sha256(Path(__file__)),
        },
        "scope": "UNVALIDATED toy intrinsic bilinear patch, controlled all-DOF prescribed loading. No free fracture, real material calibration, mesh-bias conclusion, R2 or PRETRAIN_ELIGIBLE.",
        "energy_balance": "NOT_APPLICABLE as pass criterion because all DOFs are prescribed; no dynamic balance claim.",
        "interface_area_m2": {"per_native_cell": area_cell.tolist(), "total": total_area, "GcA_J": gc_area},
        "fsync_note": "Runner JSONL trace fsyncs each line; per-point NPZ files are closed per save but runner did not fsync NPZ data.",
    }
    write_json(DERIVED / "analysis_partial_v2.json", summary)
    print(json.dumps({
        "analysis_status": summary["analysis_status"],
        "candidate_constitutive_check": summary["candidate_constitutive_check"],
        "trace_sample_pvtu_vtu_count": [len(trace_rows), len(sample_files), len(pvtu_files), len(vtu_files)],
        "indices_first_last": [int(pidx[0]), int(pidx[-1])],
        "last_commanded_opening_nm": float(command[-1] * 1.0e9),
        "delta_c_nm": float(deltac * 1.0e9),
        "dissipated_min_raw_J": float(dissipated.min()),
        "damage_last_mean": float(damage[-1].mean()),
        "png": str(png_path),
    }, allow_nan=False, sort_keys=True))
    return summary


if __name__ == "__main__":
    main()