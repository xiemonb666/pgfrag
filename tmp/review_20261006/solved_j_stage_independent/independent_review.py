from __future__ import annotations

import hashlib
import json
import math
import os
import time
from pathlib import Path

OUT = Path(__file__).resolve().parent
OUT.mkdir(parents=True, exist_ok=True)
os.environ["MPLCONFIGDIR"] = str(OUT / ".mplconfig")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(r"C:\Users\xiemo\Desktop\ShardForge")
STATIC_PATH = ROOT / "pgfrag/tmp/review_20261006/modeiii_pu_static_solve_v3_20261006/runs/run01/static_solve_diagnostic.json"
COUPLING_RUN = ROOT / "pgfrag/tmp/review_20261006/solved_field_j_coupling_20261006/runs/run01"
OUTER_RECEIPT = COUPLING_RUN / "run_receipt.json"
EXPECTED_STATIC_SHA = "39cb2dabe767998baf5e84b038f39510d6ef7db10bcd6c77298b480045cb24d3"
EXPECTED_WORKER_SHA = "1377e428a12f53bbc9dc9537c4dbf5d0c3f67a937bcdeceee6ddc6c98dfdfac9"
EXPECTED_IDS = 54
MAX_INPUT_BYTES = 3 * 1024 * 1024


def read_json(path: Path) -> tuple[dict, bytes, str]:
    raw = path.read_bytes()
    return json.loads(raw.decode("utf-8")), raw, hashlib.sha256(raw).hexdigest()


def vec(text: str, count: int) -> np.ndarray:
    out = np.fromstring(text, sep=" ", dtype=np.float64)
    if out.shape != (count,) or not np.all(np.isfinite(out)):
        raise ValueError(f"saved integrator vector malformed: {text!r}")
    return out


def fsum(values) -> float:
    return math.fsum(float(value) for value in values)


def json_safe(value):
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    return value


started = time.perf_counter()
primary, static_raw, static_sha = read_json(STATIC_PATH)
outer_receipt, receipt_raw, receipt_sha = read_json(OUTER_RECEIPT)
worker_path = Path(outer_receipt["worker_staging_path"])
worker, worker_raw, worker_sha = read_json(worker_path)
input_bytes = len(static_raw) + len(receipt_raw) + len(worker_raw)
if input_bytes > MAX_INPUT_BYTES:
    raise ValueError(f"saved input budget exceeded: {input_bytes} bytes")
if static_sha != EXPECTED_STATIC_SHA or worker_sha != EXPECTED_WORKER_SHA:
    raise ValueError("frozen saved input SHA256 mismatch")
if outer_receipt.get("worker_staging_sha256", "").lower() != worker_sha:
    raise ValueError("outer receipt does not bind the staging JSON SHA256")
if worker.get("source_solution", {}).get("sha256", "").lower() != static_sha:
    raise ValueError("worker source_solution hash does not bind primary static JSON")

provider = worker["provider"]
contours = worker["contours"]
mu = float(worker["J_integrator"]["shear_modulus_Pa"])
j_exact = float(worker["J_integrator"]["analytic_reference_J_J_per_m2_comparison_only"])
if not math.isfinite(mu) or mu <= 0 or not math.isfinite(j_exact):
    raise ValueError("invalid saved material or analytic reference")

primary_levels = {row["name"]: row for row in primary["volume_levels"]}
primary_ids = primary["registry"]["ordered_canonical_coefficient_ids"]
worker_levels = {row["name"]: row for row in worker["source_solution"]["levels"]}
if len(primary_ids) != EXPECTED_IDS or len(set(primary_ids)) != EXPECTED_IDS:
    raise ValueError("primary canonical coefficient ID list is not 54 unique IDs")
if set(primary_levels) != set(worker_levels) or len(primary_levels) != 2:
    raise ValueError("saved static/coupling solve-level sets differ")

binding_checks = []
for name in sorted(primary_levels):
    static_level = primary_levels[name]
    staged_level = worker_levels[name]
    binding = staged_level["semantic_coefficient_binding"]
    ids = binding["ordered_canonical_coefficient_ids"]
    coeffs = [float(item) for item in binding["coefficients_m_in_provider_slot_order"]]
    solved = [float(item) for item in static_level["solution"]["coefficients_m"]]
    digest_payload = json.dumps(
        list(zip(ids, coeffs)), ensure_ascii=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    binding_sha = hashlib.sha256(digest_payload).hexdigest()
    binding_checks.append(
        {
            "volume_level": name,
            "static_theta_order": static_level["theta_order"],
            "static_volume_point_count": static_level["volume_point_count"],
            "static_solution_status": static_level["solution"]["status"],
            "stage_solution_status": staged_level["solution_status"],
            "canonical_id_count": len(ids),
            "ids_equal_primary_registry_order": ids == primary_ids,
            "coefficients_equal_primary_solution_exactly": coeffs == solved,
            "maximum_absolute_coefficient_difference_m": max(abs(a - b) for a, b in zip(coeffs, solved)),
            "binding_hash_equal_saved": binding_sha.lower() == binding["semantic_binding_sha256"].lower(),
            "binding_sha256_recomputed": binding_sha,
        }
    )

if not all(
    row["ids_equal_primary_registry_order"]
    and row["coefficients_equal_primary_solution_exactly"]
    and row["binding_hash_equal_saved"]
    and row["canonical_id_count"] == EXPECTED_IDS
    for row in binding_checks
):
    raise ValueError("semantic coefficient binding does not exactly match the primary solution")

contour_records = []
flat_rows = []
unknowns = []
saved_gradient_max_abs_difference = 0.0
saved_coordinate_max_abs_difference = 0.0
saved_normal_max_abs_difference = 0.0
saved_weight_max_abs_difference = 0.0
saved_weighted_term_max_abs_difference = 0.0
slice_cursor = 0

for contour in contours:
    cid = contour["contour_id"]
    xyz = np.asarray(contour["coordinates_m"], dtype=np.float64)
    normals = np.asarray(contour["normals_outward"], dtype=np.float64)
    weights = np.asarray(contour["weights_ds_m"], dtype=np.float64)
    begin, stop = [int(item) for item in contour["provider_query_slice_half_open"]]
    n = int(contour["sample_count"])
    if xyz.shape != (n, 3) or normals.shape != (n, 3) or weights.shape != (n,):
        raise ValueError(f"saved contour shape mismatch: {cid}")
    if begin != slice_cursor or stop - begin != n:
        raise ValueError(f"saved provider slice gap or overlap: {cid}")
    if not (np.all(np.isfinite(xyz)) and np.all(np.isfinite(normals)) and np.all(np.isfinite(weights))):
        raise ValueError(f"nonfinite saved contour geometry: {cid}")
    if np.any(weights <= 0.0):
        raise ValueError(f"nonpositive saved ds weight: {cid}")

    axes = np.asarray(contour["ellipse_semiaxes_m"], dtype=np.float64)
    if axes.shape != (2,) or np.any(axes <= 0):
        raise ValueError(f"invalid saved ellipse semiaxes: {cid}")
    angle = float(contour["rotation_rad"])
    cosine, sine = math.cos(angle), math.sin(angle)
    local_x = xyz[:, 0] * cosine + xyz[:, 1] * sine
    local_y = -xyz[:, 0] * sine + xyz[:, 1] * cosine
    ellipse_residual = (local_x / axes[0]) ** 2 + (local_y / axes[1]) ** 2 - 1.0
    expected_local_normal = np.column_stack((local_x / axes[0] ** 2, local_y / axes[1] ** 2))
    expected_normal = np.column_stack(
        (
            expected_local_normal[:, 0] * cosine - expected_local_normal[:, 1] * sine,
            expected_local_normal[:, 0] * sine + expected_local_normal[:, 1] * cosine,
        )
    )
    expected_normal /= np.linalg.norm(expected_normal, axis=1)[:, None]
    normal_unit_error = float(np.max(np.abs(np.linalg.norm(normals, axis=1) - 1.0)))
    normal_geometry_error = float(np.max(np.linalg.norm(normals[:, :2] - expected_normal, axis=1)))

    record_by_level = {}
    for name in sorted(primary_levels):
        level_record = contour["by_solved_volume_level"][name]
        grad = np.asarray(level_record["solved_grad_u_z_m_per_m_full_xyz"], dtype=np.float64)
        uz = np.asarray(level_record["solved_uz_m"], dtype=np.float64)
        if grad.shape != (n, 3) or uz.shape != (n,) or not (np.all(np.isfinite(grad)) and np.all(np.isfinite(uz))):
            raise ValueError(f"saved contour solved-field shape/nonfinite mismatch: {cid}/{name}")

        ux, uy = grad[:, 0], grad[:, 1]
        nx, ny = normals[:, 0], normals[:, 1]
        energy_density_2d = (mu / 2.0) * (ux * ux + uy * uy)
        traction_times_ux = mu * (ux * nx + uy * ny) * ux
        density = energy_density_2d * nx - traction_times_ux
        contribution = density * weights
        independent_j = fsum(contribution)
        energy_term = fsum(energy_density_2d * nx * weights)
        traction_term = fsum(traction_times_ux * weights)

        integrator = level_record["integrator_result"]
        status = integrator.get("status")
        if status != "OK":
            unknowns.append(
                {"contour_id": cid, "volume_level": name, "status": status, "reason": integrator.get("reason")}
            )
        stored_j = integrator.get("J_J_per_m2")
        stored_j_delta = None if stored_j is None else independent_j - float(stored_j)

        samples = integrator.get("samples") or []
        if len(samples) != n or int(integrator.get("sample_count", -1)) != n:
            raise ValueError(f"integrator saved sample count differs from contour rows: {cid}/{name}")
        local_gradient_error = 0.0
        local_coordinate_error = 0.0
        local_normal_error = 0.0
        local_weight_error = 0.0
        local_contribution_error = 0.0
        sample_saved_contributions = []
        for index, sample in enumerate(samples):
            sample_xyz = vec(sample["coordinate_m"], 3)
            sample_n = vec(sample["normal_outward"], 3)
            sample_g = vec(sample["grad_u_z_m_per_m"], 3)
            sample_w = float(sample["weight_ds_m"])
            local_coordinate_error = max(local_coordinate_error, float(np.max(np.abs(sample_xyz - xyz[index]))))
            local_normal_error = max(local_normal_error, float(np.max(np.abs(sample_n - normals[index]))))
            local_gradient_error = max(local_gradient_error, float(np.max(np.abs(sample_g - grad[index]))))
            local_weight_error = max(local_weight_error, abs(sample_w - float(weights[index])))
            if "weighted_contribution_J_per_m2" in sample:
                sample_saved_contributions.append(float(sample["weighted_contribution_J_per_m2"]))
                local_contribution_error = max(
                    local_contribution_error,
                    abs(float(sample["weighted_contribution_J_per_m2"]) - float(contribution[index])),
                )
        saved_gradient_max_abs_difference = max(saved_gradient_max_abs_difference, local_gradient_error)
        saved_coordinate_max_abs_difference = max(saved_coordinate_max_abs_difference, local_coordinate_error)
        saved_normal_max_abs_difference = max(saved_normal_max_abs_difference, local_normal_error)
        saved_weight_max_abs_difference = max(saved_weight_max_abs_difference, local_weight_error)
        saved_weighted_term_max_abs_difference = max(saved_weighted_term_max_abs_difference, local_contribution_error)

        record_by_level[name] = {
            "saved_integrator_status": status,
            "saved_integrator_reason": integrator.get("reason"),
            "sample_count": n,
            "independent_J_J_per_m2": independent_j,
            "saved_integrator_J_J_per_m2": stored_j,
            "independent_minus_saved_J_J_per_m2": stored_j_delta,
            "J_exact_comparison_only_J_per_m2": j_exact,
            "absolute_error_to_J_exact_J_per_m2": abs(independent_j - j_exact),
            "relative_error_to_J_exact": abs(independent_j - j_exact) / abs(j_exact),
            "energy_term_J_per_m2": energy_term,
            "traction_term_J_per_m2": traction_term,
            "stored_integrator_energy_term_J_per_m2": integrator.get("energy_term_integral_J_per_m2"),
            "stored_integrator_traction_term_J_per_m2": integrator.get("traction_term_integral_J_per_m2"),
            "max_abs_grad_z_m_per_m": float(np.max(np.abs(grad[:, 2]))),
            "sample_matches": {
                "max_coordinate_abs_difference_m": local_coordinate_error,
                "max_normal_abs_difference": local_normal_error,
                "max_full_xyz_gradient_abs_difference": local_gradient_error,
                "max_ds_abs_difference_m": local_weight_error,
                "max_weighted_contribution_abs_difference_J_per_m2": local_contribution_error,
                "fsum_saved_sample_contributions_J_per_m2": fsum(sample_saved_contributions),
            },
        }
        flat_rows.append(
            {
                "contour_id": cid,
                "radius_scale_label": contour["radius_scale_label"],
                "rotation_rad": float(contour["rotation_rad"]),
                "quadrature_order_per_smooth_arc": int(contour["quadrature_order_per_smooth_arc"]),
                "provider_query_slice_half_open": [begin, stop],
                "volume_level": name,
                "independent_J_J_per_m2": independent_j,
                "saved_integrator_status": status,
                "saved_integrator_J_J_per_m2": stored_j,
                "absolute_error_to_J_exact_J_per_m2": abs(independent_j - j_exact),
                "relative_error_to_J_exact": abs(independent_j - j_exact) / abs(j_exact),
            }
        )

    contour_records.append(
        {
            "contour_id": cid,
            "radius_scale_label": contour["radius_scale_label"],
            "ellipse_semiaxes_m": axes.tolist(),
            "rotation_rad": angle,
            "quadrature_order_per_smooth_arc": int(contour["quadrature_order_per_smooth_arc"]),
            "smooth_arc_count": int(contour["smooth_arc_count"]),
            "sample_count": n,
            "provider_query_slice_half_open": [begin, stop],
            "geometry_checks": {
                "maximum_ellipse_equation_absolute_residual": float(np.max(np.abs(ellipse_residual))),
                "maximum_normal_unit_error": normal_unit_error,
                "maximum_normal_vs_ellipse_outward_normal_error": normal_geometry_error,
                "maximum_abs_contour_z_m": float(np.max(np.abs(xyz[:, 2]))),
                "ds_weight_sum_m": fsum(weights),
                "minimum_ds_weight_m": float(np.min(weights)),
            },
            "by_solved_volume_level": record_by_level,
        }
    )
    slice_cursor = stop

if slice_cursor != int(provider["contour_query_count"]):
    raise ValueError("saved contour slices do not exactly fill the provider contour rows")

controls_root = worker["z_reflection_controls"]
controls_by_level = {row["volume_level"]: row for row in controls_root["points_by_solved_volume_level"]}
control_summary = {}
for name, control_row in controls_by_level.items():
    pairs = control_row["points"]
    max_du = 0.0
    max_dgrad = 0.0
    max_abs_grad_z = 0.0
    max_pair_coordinate_error = 0.0
    pair_rows = []
    for pair in pairs:
        minus_xyz = np.asarray(pair["minus_z_coordinates_m"], dtype=np.float64)
        plus_xyz = np.asarray(pair["plus_z_coordinates_m"], dtype=np.float64)
        minus_g = np.asarray(pair["minus_z_grad_u_z_m_per_m_full_xyz"], dtype=np.float64)
        plus_g = np.asarray(pair["plus_z_grad_u_z_m_per_m_full_xyz"], dtype=np.float64)
        du = float(pair["uz_plus_minus_difference_m"])
        dg_saved = np.asarray(pair["grad_plus_minus_difference_m_per_m_full_xyz"], dtype=np.float64)
        dg = plus_g - minus_g
        coordinate_error = max(
            float(np.max(np.abs(minus_xyz[:2] - plus_xyz[:2]))),
            abs(float(minus_xyz[2]) + 0.0025),
            abs(float(plus_xyz[2]) - 0.0025),
            abs(float(pair["z_offset_m"]) - 0.0025),
        )
        max_pair_coordinate_error = max(max_pair_coordinate_error, coordinate_error)
        max_du = max(max_du, abs(du))
        max_dgrad = max(max_dgrad, float(np.linalg.norm(dg)))
        max_abs_grad_z = max(max_abs_grad_z, abs(float(minus_g[2])), abs(float(plus_g[2])))
        pair_rows.append(
            {
                "pair_id": pair["pair_id"],
                "z_offset_m": float(pair["z_offset_m"]),
                "minus_z_uz_m": float(pair["minus_z_uz_m"]),
                "plus_z_uz_m": float(pair["plus_z_uz_m"]),
                "uz_plus_minus_difference_m": du,
                "minus_z_grad_u_z_m_per_m_full_xyz": minus_g.tolist(),
                "plus_z_grad_u_z_m_per_m_full_xyz": plus_g.tolist(),
                "maximum_abs_grad_z_m_per_m_for_pair": max(abs(float(minus_g[2])), abs(float(plus_g[2]))),
                "gradient_difference_matches_saved_max_abs": float(np.max(np.abs(dg - dg_saved))),
                "xy_same_and_z_is_plus_minus_2p5mm_max_abs_error_m": coordinate_error,
            }
        )
    control_summary[name] = {
        "saved_control_status": control_row["summary"]["status"],
        "pair_count": len(pairs),
        "maximum_reflected_uz_difference_m": max_du,
        "maximum_reflected_gradient_difference_m_per_m": max_dgrad,
        "maximum_abs_grad_z_m_per_m_at_controls": max_abs_grad_z,
        "maximum_control_coordinate_error_m": max_pair_coordinate_error,
        "pairs": pair_rows,
    }

control_query_count = int(controls_root["total_provider_queries"])
if slice_cursor + control_query_count != int(provider["sampled_batch_query_count"]):
    raise ValueError("contour plus control query counts differ from saved provider batch total")

path_groups = {}
for row in flat_rows:
    geom_key = (row["radius_scale_label"], row["rotation_rad"])
    path_groups.setdefault(geom_key, {})[(row["volume_level"], row["quadrature_order_per_smooth_arc"])] = row

path_sensitivity = []
for (size_label, rotation), values in sorted(path_groups.items()):
    entry = {"radius_scale_label": size_label, "rotation_rad": rotation, "by_solved_volume_level": {}}
    for name in sorted(primary_levels):
        twelve = values[(name, 12)]
        sixteen = values[(name, 16)]
        entry["by_solved_volume_level"][name] = {
            "GL12_J_J_per_m2": twelve["independent_J_J_per_m2"],
            "GL16_J_J_per_m2": sixteen["independent_J_J_per_m2"],
            "GL12_abs_error_to_J_exact_J_per_m2": twelve["absolute_error_to_J_exact_J_per_m2"],
            "GL16_abs_error_to_J_exact_J_per_m2": sixteen["absolute_error_to_J_exact_J_per_m2"],
            "GL12_vs_GL16_absolute_difference_J_per_m2": abs(
                twelve["independent_J_J_per_m2"] - sixteen["independent_J_J_per_m2"]
            ),
            "saved_integrator_status_GL12": twelve["saved_integrator_status"],
            "saved_integrator_status_GL16": sixteen["saved_integrator_status"],
        }
    path_sensitivity.append(entry)

resolution_sensitivity = {}
for name in sorted(primary_levels):
    by_order = {}
    for order in (12, 16):
        rows = [r for r in flat_rows if r["volume_level"] == name and r["quadrature_order_per_smooth_arc"] == order]
        vals = [float(r["independent_J_J_per_m2"]) for r in rows]
        by_order[str(order)] = {
            "path_count": len(vals),
            "minimum_J_J_per_m2_across_size_and_rotation": min(vals),
            "maximum_J_J_per_m2_across_size_and_rotation": max(vals),
            "range_J_J_per_m2_across_size_and_rotation": max(vals) - min(vals),
            "maximum_absolute_error_to_J_exact_J_per_m2": max(abs(value - j_exact) for value in vals),
        }
    all_vals = [float(r["independent_J_J_per_m2"]) for r in flat_rows if r["volume_level"] == name]
    resolution_sensitivity[name] = {
        "contour_path_count": len(all_vals),
        "minimum_J_J_per_m2": min(all_vals),
        "maximum_J_J_per_m2": max(all_vals),
        "range_J_J_per_m2": max(all_vals) - min(all_vals),
        "maximum_absolute_error_to_J_exact_J_per_m2": max(abs(value - j_exact) for value in all_vals),
        "by_quadrature_order_across_size_and_rotation": by_order,
    }

unknown_count = len(unknowns)
numeric = {
    "schema": "modeiii_solved_j_saved_array_independent_review.v1",
    "interpretation": "independent evaluation of the explicitly requested midsection contour formula from saved coordinates, normals, ds weights and stored full XYZ gradients; no provider, quadrature, solver or J-integrator imports/calls",
    "formula": "fsum((mu/2*(ux^2+uy^2)*nx - mu*(ux*nx+uy*ny)*ux)*ds)",
    "inputs": {
        "static_json_path": str(STATIC_PATH),
        "static_json_sha256": static_sha,
        "static_json_bytes": len(static_raw),
        "coupling_run_receipt_path": str(OUTER_RECEIPT),
        "coupling_run_receipt_sha256": receipt_sha,
        "worker_staging_path": str(worker_path),
        "worker_staging_sha256": worker_sha,
        "worker_staging_bytes": len(worker_raw),
        "aggregate_saved_json_input_bytes": input_bytes,
        "saved_input_budget_bytes": MAX_INPUT_BYTES,
        "primary_status": primary.get("status"),
        "worker_stage_status": worker.get("status"),
        "outer_worker_return_code": outer_receipt.get("worker_return_code"),
        "outer_result_status": outer_receipt.get("result_status"),
    },
    "field_and_integrator": {
        "shear_modulus_Pa": mu,
        "J_exact_comparison_only_J_per_m2": j_exact,
        "provider_sampled_rows": int(provider["sampled_batch_query_count"]),
        "contour_rows": slice_cursor,
        "z_control_rows": control_query_count,
        "contour_count": len(contours),
        "expected_query_rows_sum_matches_248": slice_cursor + control_query_count == 248,
        "stored_gradients_keep_xyz_components": True,
        "max_saved_sample_vs_stored_gradient_abs_difference_m_per_m": saved_gradient_max_abs_difference,
        "max_saved_sample_vs_contour_coordinate_abs_difference_m": saved_coordinate_max_abs_difference,
        "max_saved_sample_vs_normal_abs_difference": saved_normal_max_abs_difference,
        "max_saved_sample_vs_ds_weight_abs_difference_m": saved_weight_max_abs_difference,
        "max_saved_sample_vs_weighted_formula_term_abs_difference_J_per_m2": saved_weighted_term_max_abs_difference,
        "semantic_54_coefficient_bindings": binding_checks,
        "controls_by_solved_volume_level": control_summary,
        "unknown_contour_integrator_records_preserved": unknowns,
        "unknown_contour_integrator_record_count": unknown_count,
    },
    "contours": contour_records,
    "path_sensitivity_separated_by_size_rotation_and_order": path_sensitivity,
    "resolution_sensitivity": resolution_sensitivity,
    "contour_level_rows": flat_rows,
    "elapsed_analysis_seconds_in_script": time.perf_counter() - started,
}

numeric_path = OUT / "numeric_results.json"
numeric_path.write_text(json.dumps(json_safe(numeric), ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")

fine_name = "volume_gl4x16x4"
fig, axes = plt.subplots(1, 3, figsize=(16.8, 5.8), gridspec_kw={"width_ratios": [1.08, 1.42, 0.9]})
ax = axes[0]
max_color = 0.0
for contour in contours:
    saved_grad = np.asarray(contour["by_solved_volume_level"][fine_name]["solved_grad_u_z_m_per_m_full_xyz"], dtype=np.float64)
    saved_normal = np.asarray(contour["normals_outward"], dtype=np.float64)
    saved_weights = np.asarray(contour["weights_ds_m"], dtype=np.float64)
    sx, sy = saved_grad[:, 0], saved_grad[:, 1]
    snx, sny = saved_normal[:, 0], saved_normal[:, 1]
    weighted_contribution = ((mu / 2.0) * (sx * sx + sy * sy) * snx - mu * (sx * snx + sy * sny) * sx) * saved_weights
    max_color = max(max_color, float(np.max(np.abs(weighted_contribution))))
max_color = max(max_color, 1e-15)
scatter_for_colorbar = None
marker_by_order = {12: "o", 16: "^"}
for contour in contours:
    xyz = np.asarray(contour["coordinates_m"], dtype=np.float64)
    normals = np.asarray(contour["normals_outward"], dtype=np.float64)
    weights = np.asarray(contour["weights_ds_m"], dtype=np.float64)
    grad = np.asarray(contour["by_solved_volume_level"][fine_name]["solved_grad_u_z_m_per_m_full_xyz"], dtype=np.float64)
    ux, uy = grad[:, 0], grad[:, 1]
    nx, ny = normals[:, 0], normals[:, 1]
    weighted = ((mu / 2.0) * (ux * ux + uy * uy) * nx - mu * (ux * nx + uy * ny) * ux) * weights
    scatter_for_colorbar = ax.scatter(
        xyz[:, 0] * 1e3,
        xyz[:, 1] * 1e3,
        c=weighted,
        cmap="coolwarm",
        vmin=-max_color,
        vmax=max_color,
        s=12 + 20 * weights / max(float(np.max(weights)), 1e-15),
        marker=marker_by_order[int(contour["quadrature_order_per_smooth_arc"])],
        alpha=0.78,
        linewidths=0,
    )
ax.set_aspect("equal", adjustable="box")
ax.set_xlabel("x (mm)")
ax.set_ylabel("y (mm)")
ax.set_title("Saved z=0 contour rows\ncolor = discrete weighted J contribution")
ax.grid(True, alpha=0.25)
ax.legend(
    handles=[
        plt.Line2D([], [], marker="o", linestyle="", color="#555555", label="GL12 / arc"),
        plt.Line2D([], [], marker="^", linestyle="", color="#555555", label="GL16 / arc"),
    ],
    loc="lower left",
    fontsize=8,
)
fig.colorbar(scatter_for_colorbar, ax=ax, fraction=0.047, pad=0.02, label="weighted contribution (J/m²)")

ax = axes[1]
x = np.arange(len(contours), dtype=np.float64)
labels = [f"{c['radius_scale_label'][0].upper()} {float(c['rotation_rad']):.2f}\nGL{c['quadrature_order_per_smooth_arc']}" for c in contours]
for name, marker, color, label in (
    ("volume_gl4x8x4", "o", "#2878b5", "coarse GL4×8×4"),
    ("volume_gl4x16x4", "s", "#d4572b", "fine GL4×16×4"),
):
    values = [r["by_solved_volume_level"][name]["independent_J_J_per_m2"] for r in contour_records]
    y = 100.0 * (np.asarray(values, dtype=np.float64) / j_exact - 1.0)
    styles = ["filled" if r["by_solved_volume_level"][name]["saved_integrator_status"] == "OK" else "open" for r in contour_records]
    ax.plot(x, y, color=color, linewidth=0.9, alpha=0.7, label=label)
    for idx, state in enumerate(styles):
        ax.scatter(
            x[idx], y[idx], marker=marker, s=34, color=color if state == "filled" else "white",
            edgecolors=color, linewidths=1.1, zorder=3,
        )
ax.axhline(0.0, color="#444444", linewidth=0.9, linestyle="--", label="Jexact = 0.013")
ax.set_xticks(x, labels, rotation=0, fontsize=7)
ax.set_ylabel("100 × (J / Jexact − 1) (%)")
ax.set_title("Each saved contour path\ncoarse and fine volume solves")
ax.grid(True, axis="y", alpha=0.25)
ax.legend(fontsize=8, loc="best")

ax = axes[2]
control_names = sorted(control_summary)
pair_count = max((len(control_summary[name]["pairs"]) for name in control_names), default=0)
pair_x = np.arange(pair_count)
for name, marker, color, label in (
    ("volume_gl4x8x4", "o", "#2878b5", "coarse"),
    ("volume_gl4x16x4", "s", "#d4572b", "fine"),
):
    pairs = control_summary[name]["pairs"]
    vals = [row["maximum_abs_grad_z_m_per_m_for_pair"] for row in pairs]
    ax.plot(pair_x[:len(vals)], vals, marker=marker, color=color, linewidth=0.9, markersize=3.5, label=label)
ax.set_xticks(pair_x, [f"{i + 1}" for i in range(pair_count)], fontsize=7)
ax.set_xlabel("saved ±z control pair")
ax.set_ylabel("max |∂u_z/∂z| (m/m)")
ax.set_title("Controls at z = ±2.5 mm\nall XYZ gradients retained")
ax.grid(True, axis="y", alpha=0.25)
ax.legend(fontsize=8)
max_delta_u = max((item["maximum_reflected_uz_difference_m"] for item in control_summary.values()), default=0.0)
max_delta_grad = max((item["maximum_reflected_gradient_difference_m_per_m"] for item in control_summary.values()), default=0.0)
ax.text(
    0.03, 0.98,
    f"max |Δu| = {max_delta_u:.3g} m\nmax |Δ∇u| = {max_delta_grad:.3g} m/m",
    transform=ax.transAxes, va="top", ha="left", fontsize=7,
    bbox={"facecolor": "white", "edgecolor": "#bbbbbb", "alpha": 0.88},
)

fig.suptitle("Independent readback of saved solved-field J samples", fontsize=13)
fig.tight_layout(rect=(0, 0, 1, 0.94))
figure_path = OUT / "saved_contour_j_review.png"
fig.savefig(figure_path, dpi=170, bbox_inches="tight")
plt.close(fig)

coarse = resolution_sensitivity["volume_gl4x8x4"]
fine = resolution_sensitivity["volume_gl4x16x4"]
review_text = (
    "# 解场 J 的独立只读复核\n\n"
    f"读取 static JSON（SHA256 `{static_sha}`）与 receipt 指定的 worker stage（SHA256 `{worker_sha}`）；总输入 {input_bytes:,} B。未导入 provider、J 积分器或 solver；未发出新查询。\n\n"
    f"两档 54 个 canonical ID 与 primary registry 顺序完全相同，stage 系数与 primary solution 逐项完全相等；248 行由 {slice_cursor} 个轮廓点和 {control_query_count} 个 z 控制点组成，8 个 slice 连续且无重叠。积分器样点的坐标、法向、ds 和完整 XYZ 梯度均与保存数组核对；梯度最大差 `{saved_gradient_max_abs_difference:.3g}`。\n\n"
    f"按指定式 `fsum((μ/2·(ux²+uy²)·nx − μ·(ux·nx+uy·ny)·ux)·ds)` 重算。Jexact={j_exact:.6g} 只作对照。粗解轮廓 J 范围 [{coarse['minimum_J_J_per_m2']:.8g}, {coarse['maximum_J_J_per_m2']:.8g}]，最大绝对误差 {coarse['maximum_absolute_error_to_J_exact_J_per_m2']:.3g}；细解范围 [{fine['minimum_J_J_per_m2']:.8g}, {fine['maximum_J_J_per_m2']:.8g}]，最大绝对误差 {fine['maximum_absolute_error_to_J_exact_J_per_m2']:.3g} J/m²。GL12/16、尺寸与旋转已分组见 numeric JSON。\n\n"
    f"±2.5 mm 控制的最大 |Δu|={max_delta_u:.3g} m，最大 |Δ∇u|={max_delta_grad:.3g} m/m；这保留了非严格二维证据。J 仅是 z=0 截面的指定泛函，不能据此认定全域二维或物理 G。保存的 UNKNOWN 数={unknown_count}；本审阅不覆盖父层 validator 的 UNKNOWN，也不设置 PASS。\n\n"
    "图仅绘制实际保存采样点，无跨裂面插值或虚构几何。执行返回码、耗时、hash 与 touched-file 清单见 `execution_receipt.json`。\n"
)
if len(review_text) > 1200:
    raise ValueError(f"REVIEW.md exceeds 1200 characters: {len(review_text)}")
(OUT / "REVIEW.md").write_text(review_text, encoding="utf-8")
print(json.dumps({"status": "OK", "analysis_elapsed_seconds": time.perf_counter() - started, "input_bytes": input_bytes, "unknown_count": unknown_count, "outputs": [str(numeric_path), str(figure_path), str(OUT / "REVIEW.md")]}, ensure_ascii=False))
