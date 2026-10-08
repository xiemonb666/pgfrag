"""Build bounded analytic contour data, run the external-data integrator, plot it."""

from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path
import sys
from typing import Any

import numpy as np
from PIL import Image, ImageDraw, ImageFont

from j_integral import INPUT_UNITS, integrate_modeiii_j


sys.dont_write_bytecode = True
REPO_ROOT = Path(__file__).resolve().parents[3]
EXPERIMENT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = REPO_ROOT / "pgfrag" / "tmp" / "review_20261006" / "modeiii_j_integral_20261006"
FIELD_REFERENCE_DIR = REPO_ROOT / "pgfrag" / "experiments" / "modeiii_field_reference_20261005"
OPERATOR_DIR = REPO_ROOT / "pgfrag" / "experiments" / "modeiii_pu_operator_20261005"


def _read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _load_field_reference():
    source = FIELD_REFERENCE_DIR / "field_reference.py"
    spec = importlib.util.spec_from_file_location("modeiii_field_reference_readonly", source)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load read-only field reference: {source}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _cut_parameter_for_negative_x(a: float, b: float, rotation: float) -> float:
    # Solve global y(t)=0 for the two ellipse intersections, then choose x<0.
    root = math.atan2(-a * math.sin(rotation), b * math.cos(rotation))
    candidates = [root % (2.0 * math.pi), (root + math.pi) % (2.0 * math.pi)]
    for candidate in candidates:
        x = a * math.cos(candidate) * math.cos(rotation) - b * math.sin(candidate) * math.sin(rotation)
        if x < 0.0:
            return candidate
    raise RuntimeError("failed to identify the negative-x branch-face crossing")


def _ellipse_quadrature(a: float, b: float, rotation: float, order_per_arc: int):
    """Gauss-Legendre samples on two smooth arcs split at the exact cut crossing."""
    nodes, weights = np.polynomial.legendre.leggauss(order_per_arc)
    cut_t = _cut_parameter_for_negative_x(a, b, rotation)
    intervals = [(0.0, cut_t), (cut_t, 2.0 * math.pi)]
    xyz_rows: list[np.ndarray] = []
    normal_rows: list[np.ndarray] = []
    ds_rows: list[np.ndarray] = []
    for lower, upper in intervals:
        mid = 0.5 * (lower + upper)
        half = 0.5 * (upper - lower)
        t = mid + half * nodes
        dt_weight = half * weights
        cos_t = np.cos(t)
        sin_t = np.sin(t)
        cos_phi = math.cos(rotation)
        sin_phi = math.sin(rotation)
        x_unrotated = a * cos_t
        y_unrotated = b * sin_t
        x = cos_phi * x_unrotated - sin_phi * y_unrotated
        y = sin_phi * x_unrotated + cos_phi * y_unrotated
        tangent_speed = np.hypot(a * sin_t, b * cos_t)
        nx_unrotated = b * cos_t
        ny_unrotated = a * sin_t
        normal_length = np.hypot(nx_unrotated, ny_unrotated)
        nx = (cos_phi * nx_unrotated - sin_phi * ny_unrotated) / normal_length
        ny = (sin_phi * nx_unrotated + cos_phi * ny_unrotated) / normal_length
        xyz_rows.append(np.column_stack((x, y, np.zeros_like(x))))
        normal_rows.append(np.column_stack((nx, ny, np.zeros_like(nx))))
        ds_rows.append(dt_weight * tangent_speed)
    return (
        np.concatenate(xyz_rows, axis=0),
        np.concatenate(normal_rows, axis=0),
        np.concatenate(ds_rows, axis=0),
    )


def _draw_plot(records: list[dict[str, Any]], expected_j: float, output_path: Path) -> None:
    width, height = 1260, 760
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)

    def font(size: int):
        try:
            return ImageFont.truetype(r"C:\Windows\Fonts\arial.ttf", size=size)
        except OSError:
            return ImageFont.load_default(size=size)

    title_font, body_font, small_font = font(28), font(18), font(15)
    draw.text((60, 24), "Mode III J from analytic field contours", fill=(15, 25, 40), font=title_font)
    draw.text((60, 61), "Relative error vs KIII^2/(2 mu) [ppm]  |  manufactured field; no solver", fill=(40, 50, 60), font=body_font)
    left, top, right, bottom = 125, 135, 1210, 590
    values = [r["relative_error_ppm"] for r in records]
    span = max(1.0, max(abs(v) for v in values) * 1.25)
    y_min, y_max = -span, span

    group_keys = list(dict.fromkeys((r["radius_scale_label"], r["rotation_rad"]) for r in records))
    group_spacing = (right - left) / len(group_keys)

    def group_center(group_index: int) -> float:
        return left + (group_index + 0.5) * group_spacing

    def py(value: float) -> float:
        return bottom - (value - y_min) * (bottom - top) / (y_max - y_min)

    draw.line((left, top, left, bottom, right, bottom), fill=(35, 45, 55), width=2)
    for tick in range(-2, 3):
        val = span * tick / 2.0
        y = py(val)
        draw.line((left, y, right, y), fill=(225, 230, 235), width=1)
        draw.text((20, y - 10), f"{val:.3g}", fill=(50, 60, 70), font=small_font)
    zero_y = py(0.0)
    draw.line((left, zero_y, right, zero_y), fill=(185, 45, 45), width=2)
    draw.text((right - 205, zero_y - 24), "exact reference: 0 ppm", fill=(155, 35, 35), font=small_font)

    colors = {12: (30, 105, 190), 16: (230, 125, 25)}
    for group_index, group in enumerate(group_keys):
        center = group_center(group_index)
        for record in [r for r in records if (r["radius_scale_label"], r["rotation_rad"]) == group]:
            offset = -12 if int(record["quadrature_order_per_arc"]) == 12 else 12
            x = center + offset
            y = py(record["relative_error_ppm"])
            color = colors[int(record["quadrature_order_per_arc"])]
            draw.ellipse((x - 7, y - 7, x + 7, y + 7), fill=color, outline=(20, 20, 20))
            draw.text((x - 14, y - 28), f"{record['quadrature_order_per_arc']}", fill=color, font=small_font)
        label = f"{group[0]}, rotation={group[1]:.2f} rad"
        draw.text((center - 63, bottom + 16), label, fill=(45, 55, 65), font=small_font)

    legend_y = 101
    for x, order in ((left, 12), (left + 255, 16)):
        color = colors[order]
        draw.ellipse((x, legend_y, x + 14, legend_y + 14), fill=color, outline=(20, 20, 20))
        draw.text((x + 22, legend_y - 3), f"Gauss order {order} per arc", fill=(45, 55, 65), font=small_font)
    draw.line((left + 540, legend_y + 7, left + 575, legend_y + 7), fill=(185, 45, 45), width=2)
    draw.text((left + 585, legend_y - 3), "analytic reference", fill=(45, 55, 65), font=small_font)

    draw.text((left, bottom + 52), f"Exact analytic value: {expected_j:.12g} J/m^2; 224 field samples recorded in JSON", fill=(45, 55, 65), font=body_font)
    draw.text((left, bottom + 80), "Rotated ellipses split at the exact branch crossing; no epsilon-shifted points", fill=(45, 55, 65), font=small_font)
    image.save(output_path, format="PNG", optimize=True)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    field_case_data = _read_json(FIELD_REFERENCE_DIR / "fixed_case.json")
    operator_case_data = _read_json(OPERATOR_DIR / "fixed_case.json")
    material = operator_case_data["material_and_manufactured_field"]
    ref_load = field_case_data["loading_and_material"]
    for key in ("young_modulus_Pa", "poisson_ratio", "K_III_Pa_sqrt_m"):
        if float(ref_load[key]) != float(material[key]):
            raise RuntimeError(f"fixed material mismatch between references: {key}")

    field_reference = _load_field_reference()
    case = field_reference.ModeIIIFixedCase(
        radius_m=float(field_case_data["geometry"]["radius_m"]),
        length_m=float(field_case_data["geometry"]["length_m"]),
        K_III_Pa_sqrt_m=float(material["K_III_Pa_sqrt_m"]),
        young_modulus_Pa=float(material["young_modulus_Pa"]),
        poisson_ratio=float(material["poisson_ratio"]),
    )
    expected_j = case.K_III_Pa_sqrt_m**2 / (2.0 * case.shear_modulus_Pa)
    records: list[dict[str, Any]] = []
    total_samples = 0
    semiaxes = [(0.0030, 0.0018, "small"), (0.0060, 0.0036, "large")]
    rotations = [0.0, 0.37]
    orders = [12, 16]

    for a, b, size_label in semiaxes:
        for rotation in rotations:
            for order in orders:
                xyz, normals, ds = _ellipse_quadrature(a, b, rotation, order)
                branch_points = (xyz[:, 0] < 0.0) & (xyz[:, 1] == 0.0)
                front_points = (xyz[:, 0] == 0.0) & (xyz[:, 1] == 0.0)
                if bool(np.any(branch_points)) or bool(np.any(front_points)):
                    raise RuntimeError("quadrature sampled an exact crack-face/front point; no epsilon shift is permitted")
                field = field_reference.evaluate_field(xyz, case=case)
                if not bool(np.all(field.statuses == "OK")):
                    raise RuntimeError(f"analytic field returned non-OK status: {field.statuses}")
                gradients = np.asarray(field.displacement_gradient[:, 2, :], dtype=np.float64)
                result = integrate_modeiii_j(
                    contour_coordinates_m=xyz,
                    normals=normals,
                    weights_m=ds,
                    grad_u_z_m_per_m=gradients,
                    shear_modulus_Pa=case.shear_modulus_Pa,
                    units=INPUT_UNITS,
                    field_provenance=(
                        "analytic exact Mode III field from modeiii_field_reference_20261005/field_reference.py; "
                        "manufactured verification input, not solved field"
                    ),
                )
                if result.status != "OK" or result.J_J_per_m2 is None:
                    raise RuntimeError(f"J integration failed: {result.to_dict()}")
                total_samples += len(xyz)
                relative_error = (result.J_J_per_m2 - expected_j) / expected_j
                records.append(
                    {
                        "radius_scale_label": size_label,
                        "ellipse_semiaxes_m": {"a": a, "b": b},
                        "rotation_rad": rotation,
                        "quadrature_order_per_arc": order,
                        "smooth_arc_count": 2,
                        "sample_count": len(xyz),
                        "analytic_reference_J_J_per_m2": expected_j,
                        "relative_error": relative_error,
                        "relative_error_ppm": relative_error * 1.0e6,
                        "integrator_result": result.to_dict(),
                    }
                )

    if total_samples > 256:
        raise RuntimeError(f"validation budget exceeded: {total_samples} > 256")
    max_relative_error = max(abs(record["relative_error"]) for record in records)
    if max_relative_error > 1.0e-4:
        raise RuntimeError(f"analytic contour error exceeds the declared 1e-4 tolerance: {max_relative_error:.6g}")
    for size_label in {record["radius_scale_label"] for record in records}:
        for rotation in rotations:
            pair = [record for record in records if record["radius_scale_label"] == size_label and record["rotation_rad"] == rotation]
            lower = next(record for record in pair if record["quadrature_order_per_arc"] == min(orders))
            higher = next(record for record in pair if record["quadrature_order_per_arc"] == max(orders))
            if abs(higher["relative_error"]) >= abs(lower["relative_error"]):
                raise RuntimeError(f"higher angular order did not reduce error for {size_label}, rotation={rotation}")

    payload = {
        "artifact_id": "modeiii_j_integral_20261006",
        "classification": "ANALYTIC_INTEGRATOR_VALIDATION_NOT_SOLVED_FIELD",
        "interpretation": "This verifies formula implementation on a manufactured analytic field. It does not establish a solved field, a measured/validated G, crack growth, or material calibration.",
        "formula": "J = integral_[Gamma] [W*n_x - t_z*(du_z/dx)] ds",
        "mode_iii_reduction": {
            "strain_energy_density": "W = mu/2 * ((du_z/dx)^2 + (du_z/dy)^2)",
            "traction_component": "t_z = mu * (grad(u_z) dot n)",
            "expected_straight_ahead_value": "J = G = KIII^2/(2*mu)",
        },
        "coordinate_and_units": {
            "front_axis": "global +z; straight front",
            "contour_plane": "global x-y at fixed z=0",
            "crack_surface": "x<0, y=0",
            "coordinates": "m",
            "outward_normals": "dimensionless unit vectors in x-y plane",
            "line_weights": "ds in m",
            "grad_u_z": "m/m (dimensionless displacement gradient)",
            "shear_modulus": "Pa",
            "energy_density_and_integrand_density": "J/m^3 and Pa, respectively",
            "J": "J/m^2",
        },
        "fixed_material_source": {
            "path": "pgfrag/experiments/modeiii_pu_operator_20261005/fixed_case.json",
            "E_Pa": case.young_modulus_Pa,
            "nu": case.poisson_ratio,
            "mu_Pa": case.shear_modulus_Pa,
            "KIII_Pa_sqrt_m": case.K_III_Pa_sqrt_m,
            "analytic_reference_J_J_per_m2": expected_j,
        },
        "field_reference_source": {
            "path": "pgfrag/experiments/modeiii_field_reference_20261005/field_reference.py",
            "fixed_case_path": "pgfrag/experiments/modeiii_field_reference_20261005/fixed_case.json",
            "role": "analytic manufactured field used solely to supply external displacement gradients",
        },
        "contour_construction": {
            "type": "rotated ellipses with outward normals and physical ds weights",
            "smooth_arc_quadrature": "Gauss-Legendre over two arcs split at the exact negative-x branch crossing; no epsilon-shifted coordinates",
            "semiaxes_m": [{"a": a, "b": b} for a, b, _ in semiaxes],
            "rotations_rad": rotations,
            "quadrature_orders_per_arc": orders,
        },
        "sample_budget": {"maximum_total_samples": 256, "actual_total_samples": total_samples},
        "maximum_absolute_relative_error": max_relative_error,
        "sources": [
            {
                "title": "J. R. Rice, A Path Independent Integral and the Approximate Analysis of Strain Concentration by Notches and Cracks (1968)",
                "doi_url": "https://doi.org/10.1115/1.3601206",
                "primary_full_text_url": "https://esag.harvard.edu/rice/015_Rice_PathIndepInt_JAM68.pdf",
                "relevant_equations": "Eq. (8a) gives the contour integrand; Eq. (10) gives J=(1+nu)KIII^2/E for Mode III, equal to KIII^2/(2mu).",
            }
        ],
        "cases": records,
    }
    json_path = OUTPUT_DIR / "modeiii_j_integral_analytic_validation.json"
    with json_path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False, allow_nan=False)
        handle.write("\n")
    png_path = OUTPUT_DIR / "modeiii_j_integral_analytic_validation.png"
    _draw_plot(records, expected_j, png_path)
    print(json.dumps({"status": "OK", "samples": total_samples, "expected_J_J_per_m2": expected_j, "max_abs_relative_error": max_relative_error, "json": str(json_path), "png": str(png_path)}, indent=2))


if __name__ == "__main__":
    main()
