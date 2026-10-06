"""Independent readback of one saved local-front K; never assembles or solves."""
from __future__ import annotations

import hashlib
import json
import math
import re
from collections import Counter
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw


ROOT = Path(r"C:\Users\xiemo\Desktop\ShardForge")
INPUT = ROOT / "pgfrag/tmp/review_20261006/local_front_operator_20261006/runs/run01/.local_front_scalar_operator.json.149216.worker.json"
OUT = ROOT / "pgfrag/tmp/review_20261006/full3d_and_local_operator_independent"
EXPECTED_INPUT_SHA = "ebdea44857ab2b1a358fc679d0b1e4996943905a5996ba83f8f8f1ba68340762"


def check(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def native(value):
    """Reject NumPy scalar leakage and non-finite values before JSON encoding."""
    if isinstance(value, dict):
        return {str(k): native(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [native(v) for v in value]
    if isinstance(value, bool) or value is None or isinstance(value, (str, int)):
        return value
    if isinstance(value, float):
        check(math.isfinite(value), "non-finite native float")
        return value
    raise TypeError(f"non-native report scalar: {type(value).__name__}")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    input_bytes = INPUT.read_bytes()
    input_sha = hashlib.sha256(input_bytes).hexdigest()
    check(input_sha == EXPECTED_INPUT_SHA, "saved staging SHA changed")
    d = json.loads(input_bytes)
    op = d["operator_diagnostics"]
    K = np.asarray(op["stiffness_matrix_N_per_m"], dtype=np.float64)
    check(K.shape == (140, 140), f"unexpected K shape {K.shape}")
    check(bool(np.isfinite(K).all()), "K contains NaN or infinity")
    ksha = hashlib.sha256(K.tobytes(order="C")).hexdigest()
    check(ksha == op["stiffness_sha256_float64_c_order"], "K byte hash differs from worker receipt")

    keys = d["registry"]["ordered_semantic_keys"]
    canonical = d["registry"]["canonical_ids"]
    check(len(keys) == len(canonical) == 140, "registry semantic vector length mismatch")
    check([k["slot"] for k in keys] == list(range(140)), "registry slots not canonical 0..139")
    check([k["canonical_id"] for k in keys] == canonical, "canonical id order mismatch")
    check(len(set(canonical)) == 140, "duplicate canonical semantic ids")
    basis_counts = Counter(k["basis_kind"] for k in keys)
    expected_counts = {"base": 100, "front_f1_shifted": 30, "heaviside": 10}
    check(dict(basis_counts) == expected_counts, f"semantic basis counts {dict(basis_counts)}")

    asym = float(np.linalg.norm(K - K.T, ord="fro") / np.linalg.norm(K, ord="fro"))
    eig = np.linalg.eigvalsh(K)
    tol = float(op["absolute_near_zero_threshold_N_per_m"])
    raw_negative = int(np.count_nonzero(eig < 0.0))
    significant_negative = int(np.count_nonzero(eig < -tol))
    near_zero = int(np.count_nonzero(np.abs(eig) <= tol))
    check(raw_negative == int(op["raw_negative_eigenvalue_count"]), "raw negative count differs")
    check(significant_negative == int(op["negative_beyond_roundoff_count"]), "significant negative count differs")
    check(near_zero == int(op["near_zero_count"]), "near-zero count differs")
    eig_saved = np.asarray(op["eigenvalues_N_per_m"], dtype=np.float64)
    check(eig_saved.shape == (140,), "saved eigenvalue vector shape mismatch")
    eig_delta = float(np.max(np.abs(eig - eig_saved)))

    # Reconstruct the scalar translation and affine Uz=x coefficients directly
    # from the saved canonical semantic IDs and fixed lattice axes.
    t = np.asarray([1.0 if k["basis_kind"] == "base" else 0.0 for k in keys], dtype=np.float64)
    scaled_kt = float(np.linalg.norm(K @ t) / (np.linalg.norm(K, ord="fro") * np.linalg.norm(t)))
    saved_t = np.asarray(d["constant_translation"]["coefficient_vector_m"], dtype=np.float64)
    check(np.array_equal(t, saved_t), "saved translation differs from semantic reconstruction")
    mean_row = np.asarray(d["volume_mean_row"]["normalized_mean_shape_row"], dtype=np.float64)
    check(mean_row.shape == (140,), "normalized mean row shape mismatch")
    mean_translation = float(mean_row @ t)

    axes = d["fixed_case"]["support_coverage"]["source_lattice_axes_m"]
    ux = []
    for k in keys:
        if k["basis_kind"] == "base":
            m = re.fullmatch(r"lf-x(\d+)-y(\d+)-z(\d+)", k["source_id"])
            check(m is not None, f"unexpected source id {k['source_id']}")
            ix = int(m.group(1))
            ux.append(float(axes["x"][ix]))
        else:
            ux.append(0.0)
    affine = np.asarray(ux, dtype=np.float64)
    affine_saved = np.asarray(d["affine_base_patch"]["coefficient_vector_m"], dtype=np.float64)
    check(np.array_equal(affine, affine_saved), "saved affine coefficients differ from semantic reconstruction")
    affine_quadratic_form = float(affine @ K @ affine)
    affine_energy = float(0.5 * affine_quadratic_form)
    mu = float(d["material_and_model"]["mu_Pa"])
    weights = np.asarray(d["quadrature"]["volume_weights_m3"], dtype=np.float64)
    coords = np.asarray(d["quadrature"]["coordinates_xyz_m"], dtype=np.float64)
    expected_volume = float(d["quadrature"]["expected_box_volume_m3"])
    volume_sum = float(np.sum(weights, dtype=np.float64))
    exact_affine_energy = float(0.5 * mu * expected_volume)
    check(weights.shape == (256,) and coords.shape == (256, 3), "quadrature saved-array shapes mismatch")
    check(bool(np.isfinite(weights).all() and np.isfinite(coords).all()), "quadrature arrays non-finite")
    check(bool(np.all(weights > 0.0)), "non-positive quadrature weight")
    check(abs(volume_sum - expected_volume) <= 1e-18, "quadrature volume sum mismatch")

    receipt = {
        "review_kind": "INDEPENDENT_SAVED_ARRAY_READBACK_ONLY",
        "promotion_result": "OUTER125_FAILED_NOT_PROMOTED",
        "input_path": str(INPUT.relative_to(ROOT)).replace("\\", "/"),
        "input_sha256": input_sha,
        "worker_status": str(d["status"]),
        "worker_assembly_status": str(d["assembly_status"]),
        "worker_solver_invoked": bool(d["solver_invoked"]),
        "K_shape": [int(x) for x in K.shape],
        "K_sha256_recomputed": ksha,
        "K_finite": bool(np.isfinite(K).all()),
        "relative_frobenius_asymmetry_recomputed": asym,
        "eigenvalue_count": int(eig.size),
        "eigenvalue_saved_vector_max_abs_delta_N_per_m": eig_delta,
        "eigenvalue_min_N_per_m": float(eig[0]),
        "eigenvalue_max_N_per_m": float(eig[-1]),
        "raw_negative_count_recomputed": raw_negative,
        "negative_beyond_roundoff_count_recomputed": significant_negative,
        "near_zero_count_recomputed": near_zero,
        "extra_near_zero_after_expected_scalar_translation": int(near_zero - 1),
        "near_zero_threshold_N_per_m": tol,
        "no_regularization_or_rank_truncation_claim_in_worker": bool(op["no_regularization_or_rank_truncation"]),
        "semantic_counts_reconstructed": {str(k): int(v) for k, v in sorted(basis_counts.items())},
        "registry_order_and_uniqueness_pass": True,
        "normalized_mean_row_sum": float(np.sum(mean_row, dtype=np.float64)),
        "translation_K_residual_scaled_recomputed": scaled_kt,
        "normalized_mean_of_translation_m": mean_translation,
        "quadrature_point_count": int(weights.size),
        "quadrature_all_weights_positive": bool(np.all(weights > 0.0)),
        "quadrature_volume_sum_m3": volume_sum,
        "expected_box_volume_m3": expected_volume,
        "quadrature_coords_finite": bool(np.isfinite(coords).all()),
        "affine_coefficients_reconstructed_from_semantic_lattice": True,
        "affine_quadratic_form_dTKd_Nm": affine_quadratic_form,
        "affine_discrete_energy_J_recomputed": affine_energy,
        "affine_exact_scalar_energy_J_recomputed": exact_affine_energy,
        "affine_energy_relative_delta": float(abs(affine_energy - exact_affine_energy) / exact_affine_energy),
        "full_3d_elasticity_claim": bool(d["material_and_model"]["full_3d_isotropic_elasticity_claim"]),
        "arrays_are_saved_source_data_only_no_reassembly_or_solver": True,
        "limitations": [
            "This independent readback does not change the failed OUTER125 promotion receipt.",
            "The artifact is a scalar antiplane Poisson operator diagnostic, not a full 3-D elasticity solve.",
            "This review independently checks persisted arrays and algebraic probes only; it does not rerun MLS or assemble K.",
        ],
    }
    receipt = native(receipt)
    out_json = OUT / "local_saved_array_audit.json"
    out_json.write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n", encoding="utf-8")

    # One small, dependency-light PNG: sorted eigenvalue magnitudes on a log scale.
    width, height = 1000, 620
    im = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(im)
    left, right, top, bottom = 96, 40, 62, 86
    x0, x1, y0, y1 = left, width - right, top, height - bottom
    draw.text((left, 18), "Saved local-front K spectrum (independent eigvalsh)", fill=(20, 20, 20))
    draw.text((left, height - 45), "eigenvalue rank (sorted ascending)", fill=(20, 20, 20))
    draw.text((14, top + 4), "log10(|lambda|/max|lambda|)", fill=(20, 20, 20))
    for yval in (0, -4, -8, -12, -16):
        y = int(y0 + (0 - yval) / 16.0 * (y1 - y0))
        draw.line((x0, y, x1, y), fill=(220, 220, 220), width=1)
        draw.text((42, y - 6), str(yval), fill=(70, 70, 70))
    draw.line((x0, y0, x0, y1), fill=(25, 25, 25), width=2)
    draw.line((x0, y1, x1, y1), fill=(25, 25, 25), width=2)
    max_abs = max(float(np.max(np.abs(eig))), np.finfo(np.float64).tiny)
    plotted = np.log10(np.maximum(np.abs(eig) / max_abs, 1e-16))
    plotted = np.clip(plotted, -16.0, 0.0)
    pts = []
    for i, yval in enumerate(plotted):
        x = int(x0 + i / max(1, len(eig) - 1) * (x1 - x0))
        y = int(y0 + (0.0 - float(yval)) / 16.0 * (y1 - y0))
        pts.append((x, y))
    draw.line(pts, fill=(35, 89, 145), width=2)
    for i, lam in enumerate(eig):
        if abs(float(lam)) <= tol or float(lam) < -tol:
            x, y = pts[i]
            color = (205, 48, 48) if float(lam) < -tol else (230, 120, 20)
            draw.ellipse((x - 4, y - 4, x + 4, y + 4), fill=color)
    draw.text((x0, y1 + 14), "0", fill=(30, 30, 30))
    draw.text((x1 - 24, y1 + 14), "139", fill=(30, 30, 30))
    draw.text((x0 + 210, y1 + 14), "orange: expected scalar translation", fill=(160, 84, 18))
    png = OUT / "local_saved_K_spectrum.png"
    im.save(png, format="PNG", optimize=True)

    # Enforce requested report size bound and confirm strict JSON readability.
    check(out_json.stat().st_size < 2 * 1024 * 1024, "audit JSON exceeds 2 MiB")
    check(png.stat().st_size < 2 * 1024 * 1024, "PNG exceeds 2 MiB")
    check(json.loads(out_json.read_text(encoding="utf-8"))["input_sha256"] == EXPECTED_INPUT_SHA, "saved JSON reread mismatch")
    print(json.dumps({
        "input_sha256": input_sha,
        "K_sha256": ksha,
        "eigen_min": float(eig[0]),
        "eigen_max": float(eig[-1]),
        "raw_negative": raw_negative,
        "significant_negative": significant_negative,
        "near_zero": near_zero,
        "translation_scaled_residual": scaled_kt,
        "affine_energy_J": affine_energy,
        "exact_energy_J": exact_affine_energy,
        "audit_json_bytes": out_json.stat().st_size,
        "spectrum_png_bytes": png.stat().st_size,
    }, allow_nan=False))


if __name__ == "__main__":
    main()
