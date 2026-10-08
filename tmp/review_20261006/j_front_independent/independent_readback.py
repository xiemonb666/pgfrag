"""Independent bounded readback of saved Mode III J sample rows.

This deliberately imports neither the J-integral implementation nor the
analytic field provider or any solver. It recomputes each sample term from the
saved global coordinates, normals, gradients, and physical ds weights.
"""

from __future__ import annotations

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SAVED_JSON = (
    ROOT
    / "pgfrag"
    / "tmp"
    / "review_20261006"
    / "modeiii_j_integral_20261006"
    / "modeiii_j_integral_analytic_validation.json"
)
FIXED_CASE_JSON = (
    ROOT
    / "pgfrag"
    / "experiments"
    / "modeiii_field_reference_20261005"
    / "fixed_case.json"
)
OUTPUT_JSON = Path(__file__).with_name("independent_readback.json")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def near(actual: float, expected: float, tolerance: float, label: str) -> float:
    error = abs(actual - expected)
    scale = max(abs(actual), abs(expected), 1.0e-30)
    relative = error / scale
    require(relative <= tolerance, f"{label}: relative difference {relative:.6g}")
    return relative


def main() -> None:
    payload = json.loads(SAVED_JSON.read_text(encoding="utf-8"))
    fixed_case = json.loads(FIXED_CASE_JSON.read_text(encoding="utf-8"))
    require(
        payload["classification"] == "ANALYTIC_INTEGRATOR_VALIDATION_NOT_SOLVED_FIELD",
        "saved artifact classification was changed",
    )
    fixed_load = fixed_case["loading_and_material"]
    case_data = payload["fixed_material_source"]
    mu = float(case_data["mu_Pa"])
    k_iii = float(case_data["KIII_Pa_sqrt_m"])
    expected_j = k_iii * k_iii / (2.0 * mu)
    amplitude = k_iii / mu * math.sqrt(2.0 / math.pi)
    require(mu > 0.0 and math.isfinite(mu), "invalid saved shear modulus")
    require(expected_j > 0.0 and math.isfinite(expected_j), "invalid reference J")
    for saved_key, fixed_key in (
        ("E_Pa", "young_modulus_Pa"),
        ("nu", "poisson_ratio"),
        ("KIII_Pa_sqrt_m", "K_III_Pa_sqrt_m"),
    ):
        near(
            float(case_data[saved_key]),
            float(fixed_load[fixed_key]),
            2.0e-15,
            f"fixed material {saved_key}",
        )
    radius_limit = float(fixed_case["geometry"]["radius_m"])
    length_limit = float(fixed_case["geometry"]["length_m"])

    per_case = []
    max_gradient_abs_error = 0.0
    max_normal_abs_error = 0.0
    total_samples = 0
    branch_upper = 0
    branch_lower = 0
    for record in payload["cases"]:
        result = record["integrator_result"]
        require(result["status"] == "OK", f"{record['radius_scale_label']} integral was not OK")
        require(result["field_provenance_verified"] is False, "provenance was presented as verified")
        samples = result["samples"]
        require(len(samples) == int(record["sample_count"]), "saved sample count mismatch")
        total_samples += len(samples)
        axes = record["ellipse_semiaxes_m"]
        a, b = float(axes["a"]), float(axes["b"])
        rotation = float(record["rotation_rad"])
        c, s = math.cos(rotation), math.sin(rotation)
        contributions = []
        energy_contributions = []
        traction_contributions = []
        weighted_normal_x = []
        weighted_normal_y = []
        perimeter_terms = []
        area_terms = []
        winding_terms = []
        gradient_abs_errors = []
        normal_abs_errors = []
        branch_above = 0
        branch_below = 0

        for index, sample in enumerate(samples):
            x, y, z = (float(v) for v in sample["coordinate_m"])
            nx, ny, nz = (float(v) for v in sample["normal_outward"])
            gx, gy, gz = (float(v) for v in sample["grad_u_z_m_per_m"])
            ds = float(sample["weight_ds_m"])
            require(all(math.isfinite(v) for v in (x, y, z, nx, ny, nz, gx, gy, gz, ds)),
                    "nonfinite saved sample")
            require(ds > 0.0, "nonpositive saved ds weight")
            require(abs(z) <= 1.0e-15, "sample left the global xy contour plane")
            require(abs(nz) <= 1.0e-12, "normal left the global xy frame")
            require(abs(gz) <= 1.0e-12, "saved field gradient is not z-invariant")
            require(math.hypot(x, y) < radius_limit, "contour left the manufactured cylinder")
            require(abs(math.hypot(nx, ny) - 1.0) <= 1.0e-10, "saved normal is not unit")
            require(not (x == 0.0 and y == 0.0), "saved contour hit the singular front")
            require(not (x < 0.0 and y == 0.0), "saved contour hit the exact branch face")
            if x < 0.0 and y > 0.0:
                branch_above += 1
            if x < 0.0 and y < 0.0:
                branch_below += 1

            # Reconstruct the outward ellipse normal in its local frame, then
            # rotate it into the saved global xyz frame.
            x_local = c * x + s * y
            y_local = -s * x + c * y
            local_nx = x_local / (a * a)
            local_ny = y_local / (b * b)
            local_norm = math.hypot(local_nx, local_ny)
            expected_nx = (c * local_nx - s * local_ny) / local_norm
            expected_ny = (s * local_nx + c * local_ny) / local_norm
            normal_error = math.hypot(nx - expected_nx, ny - expected_ny)
            normal_abs_errors.append(normal_error)

            # Recompute the global analytic Mode III gradient from each saved
            # point. The two signs on the negative-x branch must remain on the
            # physical upper/lower crack faces, with no new rotated cut.
            r = math.hypot(x, y)
            theta = math.atan2(y, x)
            expected_gx = amplitude * (-math.sin(0.5 * theta) / (2.0 * math.sqrt(r)))
            expected_gy = amplitude * (math.cos(0.5 * theta) / (2.0 * math.sqrt(r)))
            gradient_error = math.hypot(gx - expected_gx, gy - expected_gy)
            gradient_abs_errors.append(gradient_error)
            gradient_scale = max(abs(expected_gx), abs(expected_gy), 1.0e-3)
            require(
                gradient_error <= 2.0e-12 * gradient_scale,
                f"global analytic gradient mismatch at sample {index}",
            )
            if x < 0.0 and y > 0.0:
                require(gx < 0.0, "upper crack-face trace has the wrong global gx sign")
            if x < 0.0 and y < 0.0:
                require(gx > 0.0, "lower crack-face trace has the wrong global gx sign")
            require(
                normal_error <= 2.0e-12,
                f"rotated physical normal mismatch at sample {index}",
            )

            # Independently apply J = integral (W*n_x - t_z*u_z,x) ds.
            w = 0.5 * mu * (gx * gx + gy * gy)
            tz = mu * (gx * nx + gy * ny)
            wnx = w * nx
            minus_tz_gx = -tz * gx
            integrand = wnx + minus_tz_gx
            contribution = integrand * ds
            contributions.append(contribution)
            energy_contributions.append(wnx * ds)
            traction_contributions.append(minus_tz_gx * ds)
            weighted_normal_x.append(nx * ds)
            weighted_normal_y.append(ny * ds)
            perimeter_terms.append(ds)
            area_terms.append(0.5 * (x * nx + y * ny) * ds)
            saved_outputs = (
                ("W_J_per_m3", w),
                ("t_z_Pa", tz),
                ("W_n_x_Pa", wnx),
                ("minus_t_z_u_z_x_Pa", minus_tz_gx),
                ("integrand_Pa", integrand),
                ("weighted_contribution_J_per_m2", contribution),
            )
            for field, independently_computed in saved_outputs:
                near(float(sample[field]), independently_computed, 3.0e-13, f"saved {field}")

            next_sample = samples[(index + 1) % len(samples)]
            x_next, y_next = (float(v) for v in next_sample["coordinate_m"][:2])
            cross = x * y_next - y * x_next
            dot = x * x_next + y * y_next
            winding_terms.append(math.atan2(cross, dot))

        require(branch_above > 0 and branch_below > 0, "both physical crack-face traces were not sampled")
        branch_upper += branch_above
        branch_lower += branch_below
        j_independent = math.fsum(contributions)
        stored_j = float(result["J_J_per_m2"])
        j_vs_saved = near(j_independent, stored_j, 2.0e-12, "independent J vs saved result")
        relative_error = (j_independent - expected_j) / expected_j
        require(
            abs(float(record["relative_error"]) - relative_error) <= 2.0e-12,
            "stored relative error differs from independent recomputation",
        )
        require(abs(relative_error) <= 1.0e-4, "independent J missed the declared tolerance")

        perimeter = math.fsum(perimeter_terms)
        closure = math.hypot(math.fsum(weighted_normal_x), math.fsum(weighted_normal_y))
        closure_relative = closure / perimeter
        require(closure_relative <= 2.0e-12, "independent weighted normal closure failed")
        signed_area = math.fsum(area_terms)
        expected_area = math.pi * a * b
        area_relative = near(signed_area, expected_area, 2.0e-12, "outward normal area")
        winding = math.fsum(winding_terms)
        near(winding, 2.0 * math.pi, 2.0e-12, "ordered global contour winding")
        energy_term = math.fsum(energy_contributions)
        traction_term = math.fsum(traction_contributions)
        require(
            abs(float(result["energy_term_integral_J_per_m2"]) - energy_term)
            <= 2.0e-12 * expected_j,
            "stored energy term integral differs from independent recomputation",
        )
        require(
            abs(float(result["traction_term_integral_J_per_m2"]) - traction_term)
            <= 2.0e-12 * expected_j,
            "stored traction term integral differs from independent recomputation",
        )
        max_gradient_abs_error = max(max_gradient_abs_error, max(gradient_abs_errors))
        max_normal_abs_error = max(max_normal_abs_error, max(normal_abs_errors))
        per_case.append(
            {
                "radius_scale_label": record["radius_scale_label"],
                "ellipse_semiaxes_m": {"a": a, "b": b},
                "rotation_rad": rotation,
                "quadrature_order_per_arc": int(record["quadrature_order_per_arc"]),
                "sample_count": len(samples),
                "branch_upper_samples": branch_above,
                "branch_lower_samples": branch_below,
                "independent_J_J_per_m2": j_independent,
                "saved_J_J_per_m2": stored_j,
                "relative_error_vs_KIII2_over_2mu": relative_error,
                "relative_difference_vs_saved_J": j_vs_saved,
                "normal_closure_relative": closure_relative,
                "signed_area_m2": signed_area,
                "signed_area_relative_error": area_relative,
                "winding_angle_rad": winding,
                "energy_term_J_per_m2": energy_term,
                "traction_term_J_per_m2": traction_term,
                "max_gradient_abs_difference": max(gradient_abs_errors),
                "max_rotated_normal_abs_difference": max(normal_abs_errors),
            }
        )

    require(total_samples == int(payload["sample_budget"]["actual_total_samples"]),
            "total sample count does not match saved budget")
    require(total_samples <= int(payload["sample_budget"]["maximum_total_samples"]),
            "saved sample budget exceeded")
    for label in {case["radius_scale_label"] for case in per_case}:
        for rotation in {case["rotation_rad"] for case in per_case}:
            pair = [
                case
                for case in per_case
                if case["radius_scale_label"] == label and case["rotation_rad"] == rotation
            ]
            pair.sort(key=lambda case: case["quadrature_order_per_arc"])
            require(
                abs(pair[1]["relative_error_vs_KIII2_over_2mu"])
                < abs(pair[0]["relative_error_vs_KIII2_over_2mu"]),
                f"quadrature did not improve for {label}, rotation={rotation}",
            )
    max_relative_error = max(abs(case["relative_error_vs_KIII2_over_2mu"]) for case in per_case)
    summary = {
        "status": "INDEPENDENT_READBACK_OK",
        "source_classification": payload["classification"],
        "input": str(SAVED_JSON),
        "method": "plain-Python fsum of saved samples; no integrator, field provider, or solver import",
        "sample_count": total_samples,
        "maximum_samples": int(payload["sample_budget"]["maximum_total_samples"]),
        "KIII_Pa_sqrt_m": k_iii,
        "mu_Pa": mu,
        "analytic_J_J_per_m2": expected_j,
        "maximum_absolute_relative_error": max_relative_error,
        "branch_upper_sample_count": branch_upper,
        "branch_lower_sample_count": branch_lower,
        "max_gradient_abs_difference": max_gradient_abs_error,
        "max_rotated_normal_abs_difference": max_normal_abs_error,
        "cases": per_case,
    }
    OUTPUT_JSON.write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "status": summary["status"],
                "sample_count": total_samples,
                "analytic_J_J_per_m2": expected_j,
                "maximum_absolute_relative_error": max_relative_error,
                "branch_upper_sample_count": branch_upper,
                "branch_lower_sample_count": branch_lower,
                "max_gradient_abs_difference": max_gradient_abs_error,
                "max_rotated_normal_abs_difference": max_normal_abs_error,
                "output": str(OUTPUT_JSON),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
