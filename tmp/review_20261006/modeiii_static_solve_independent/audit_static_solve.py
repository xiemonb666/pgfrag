#!/usr/bin/env python3
"""One-pass read-only audit and visualization of saved Mode III static-solve fields."""
from __future__ import annotations
import hashlib
import json
import math
import re
import time
from pathlib import Path
from typing import Any
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

ROOT = Path(r"C:\Users\xiemo\Desktop\ShardForge")
RUN = ROOT / "pgfrag/tmp/review_20261006/modeiii_pu_static_solve_v3_20261006/runs/run01"
PRIMARY = RUN / "static_solve_diagnostic.json"
INPUT_RECEIPT = RUN / "run_receipt.json"
FREEZE = ROOT / "pgfrag/tmp/review_20261006/modeiii_pu_static_solve_v3_20261006/source_freeze.json"
CASE = ROOT / "pgfrag/experiments/modeiii_pu_static_solve_v3_20261006/static_solve_case.json"
OUT = ROOT / "pgfrag/tmp/review_20261006/modeiii_static_solve_independent"
OUT_JSON = OUT / "static_solve_independent_audit.json"
OUT_PNG = OUT / "static_solve_fields_q256_zslice.png"
PRIMARY_SHA = "39cb2dabe767998baf5e84b038f39510d6ef7db10bcd6c77298b480045cb24d3"
PRIMARY_BYTES = 369684

def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()

def norm(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())

def get(obj: Any, dotted: str) -> Any:
    value = obj
    for part in dotted.split("."):
        if not isinstance(value, dict) or part not in value:
            raise KeyError(dotted)
        value = value[part]
    return value

def walk(obj: Any, path: str = "$"):
    yield path, obj
    if isinstance(obj, dict):
        for key, value in obj.items():
            yield from walk(value, f"{path}.{key}")
    elif isinstance(obj, list):
        for index, value in enumerate(obj):
            yield from walk(value, f"{path}[{index}]")

def array(value: Any, label: str, shape: tuple[int, ...] | None = None) -> np.ndarray:
    result = np.asarray(value, dtype=np.float64)
    if shape is not None and result.shape != shape:
        raise ValueError(f"{label} shape {result.shape}, expected {shape}")
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{label} contains non-finite values")
    return result

def ids_family(item: Any) -> str:
    text = item if isinstance(item, str) else json.dumps(item, ensure_ascii=False)
    text = text.upper()
    base = bool(re.search(r"(?<![A-Z0-9])BASE(?![A-Z0-9])", text))
    front = bool(re.search(r"(?<![A-Z0-9])FRONT(?![A-Z0-9])", text))
    if base == front:
        raise ValueError(f"Ambiguous BASE/FRONT canonical ID: {item!r}")
    return "BASE" if base else "FRONT"

def levels_of(value: Any) -> list[tuple[str, dict[str, Any]]]:
    if isinstance(value, list):
        return [(str(x.get("name", f"level_{i}")), x) for i, x in enumerate(value)]
    if isinstance(value, dict):
        return [(str(k), v) for k, v in value.items() if isinstance(v, dict)]
    raise TypeError("volume_levels must be an array or object")

def collect_pins(obj: Any) -> dict[str, str]:
    """Read path/hash entries from the saved freeze manifest only."""
    pins: dict[str, str] = {}
    pkeys = ("path", "filepath", "relativepath", "sourcepath", "file")
    hkeys = ("sha256", "expectedsha256", "digest", "hash")
    def visit(value: Any):
        if isinstance(value, dict):
            fields = {norm(k): v for k, v in value.items()}
            path = next((fields[k] for k in pkeys if k in fields and isinstance(fields[k], str)), None)
            digest = next((fields[k] for k in hkeys if k in fields and isinstance(fields[k], str)
                           and re.fullmatch(r"[0-9a-fA-F]{64}", fields[k])), None)
            if path and digest:
                pins[str(path)] = str(digest).lower()
            for key, child in value.items():
                if isinstance(child, str) and re.fullmatch(r"[0-9a-fA-F]{64}", child) and re.search(r"\.(py|json|md|toml)$", str(key), re.I):
                    pins[str(key)] = child.lower()
                if isinstance(child, dict) and re.search(r"\.(py|json|md|toml)$", str(key), re.I):
                    child_fields = {norm(k): v for k, v in child.items()}
                    item_hash = next((child_fields[k] for k in hkeys if k in child_fields
                                      and isinstance(child_fields[k], str)
                                      and re.fullmatch(r"[0-9a-fA-F]{64}", child_fields[k])), None)
                    if item_hash:
                        pins[str(key)] = str(item_hash).lower()
                visit(child)
        elif isinstance(value, list):
            for child in value:
                visit(child)
    visit(obj)
    return pins

def rooted_path(text: str) -> Path | None:
    path = Path(text)
    if not path.is_absolute():
        path = ROOT / path
    try:
        path.resolve().relative_to(ROOT.resolve())
    except ValueError:
        return None
    return path

def hashes(paths: set[Path]) -> dict[str, str]:
    return {str(p): sha(p.read_bytes()) for p in sorted(paths, key=lambda x: str(x).lower()) if p.is_file()}

def receipt_pins(obj: Any) -> list[dict[str, str]]:
    found = []
    names = {"path", "filepath", "relativepath", "sourcepath", "file"}
    digest_names = ("sha256", "expected_sha256", "pre_sha256", "post_sha256", "digest", "hash")
    for _, value in walk(obj):
        if not isinstance(value, dict):
            continue
        fields = {norm(k): v for k, v in value.items()}
        path = next((fields[k] for k in names if k in fields and isinstance(fields[k], str)), None)
        if not path:
            continue
        for key in digest_names:
            digest = fields.get(norm(key))
            if isinstance(digest, str) and re.fullmatch(r"[0-9a-fA-F]{64}", digest):
                found.append({"path": path, "stage": key, "sha256": digest.lower()})
    return found

def main() -> int:
    t0 = time.monotonic()
    OUT.mkdir(parents=True, exist_ok=True)
    if OUT_JSON.exists() or OUT_PNG.exists():
        raise FileExistsError("Independent result or figure already exists; refusing overwrite")
    primary_raw = PRIMARY.read_bytes()
    primary_sha = sha(primary_raw)
    if len(primary_raw) != PRIMARY_BYTES or primary_sha != PRIMARY_SHA:
        raise ValueError(f"Primary pin mismatch before JSON parse: bytes={len(primary_raw)}, sha={primary_sha}")
    receipt_raw, freeze_raw, case_raw = INPUT_RECEIPT.read_bytes(), FREEZE.read_bytes(), CASE.read_bytes()
    doc, receipt, freeze, case = map(json.loads, (primary_raw, receipt_raw, freeze_raw, case_raw))
    receipt_pin_entries = receipt_pins(receipt)
    base_info = case["base_operator_case"]
    base_case_path = ROOT / base_info["path"]
    base_case_raw = base_case_path.read_bytes()
    if sha(base_case_raw) != str(base_info["sha256"]).lower():
        raise ValueError(f"Base fixed-case source hash mismatch: {base_case_path}")
    base_case = json.loads(base_case_raw)

    # Baseline hashes are captured before numerical work/rendering, then checked again at the end.
    pins = collect_pins(freeze)
    monitored = {PRIMARY, INPUT_RECEIPT, FREEZE, CASE, base_case_path}
    for raw_path in pins:
        path = rooted_path(raw_path)
        if path is not None and path.is_file():
            monitored.add(path)
    for entry in receipt_pin_entries:
        path = rooted_path(entry["path"])
        if path is not None and path.is_file():
            monitored.add(path)
    before = hashes(monitored - {PRIMARY, INPUT_RECEIPT, FREEZE, CASE, base_case_path})
    before.update({str(PRIMARY): primary_sha, str(INPUT_RECEIPT): sha(receipt_raw),
                   str(FREEZE): sha(freeze_raw), str(CASE): sha(case_raw), str(base_case_path): sha(base_case_raw)})
    if not pins:
        raise ValueError("Source freeze manifest contains no recognizable path/hash pins")
    source_pin_prechecks = []
    for raw_path, expected in pins.items():
        path = rooted_path(raw_path)
        if path is not None and str(path) in before:
            source_pin_prechecks.append({"path": str(path), "expected_sha256": expected,
                                         "current_sha256": before[str(path)], "matches": before[str(path)] == expected})
    if not source_pin_prechecks or not all(x["matches"] for x in source_pin_prechecks):
        raise ValueError(f"Source freeze pin mismatch before audit: {source_pin_prechecks}")

    E_path, nu_path = "material_and_manufactured_field.young_modulus_Pa", "material_and_manufactured_field.poisson_ratio"
    k_path = "material_and_manufactured_field.K_III_Pa_sqrt_m"
    r_path, l_path = "geometry.radius_m", "geometry.length_m"
    E, nu = float(get(base_case, E_path)), float(get(base_case, nu_path))
    KIII, R, L = float(get(base_case, k_path)), float(get(base_case, r_path)), float(get(base_case, l_path))
    z_bounds = array(get(base_case, "geometry.z_bounds_m"), "fixed-case z bounds", (2,))
    if L != float(z_bounds[1] - z_bounds[0]):
        raise ValueError(f"Fixed-case length {L} != z-bounds span {z_bounds[1]-z_bounds[0]}")
    if E <= 0 or R <= 0 or L <= 0 or not (-1 < nu < 0.5):
        raise ValueError(f"Bad fixed-case material/geometry: E={E}, nu={nu}, KIII={KIII}, R={R}, L={L}")
    mu = E / (2 * (1 + nu))
    U_exact = KIII * KIII * R * L / (2 * mu)
    exact_parameters = {"E_Pa": {"path": E_path, "value": E}, "nu": {"path": nu_path, "value": nu},
                        "KIII_Pa_sqrt_m": {"path": k_path, "value": KIII}, "R_m": {"path": r_path, "value": R},
                        "L_m": {"path": l_path, "value": L}}

    canonical = get(doc, "registry.ordered_canonical_coefficient_ids")
    if not isinstance(canonical, list) or len(canonical) != 54:
        raise ValueError("Expected 54 registry.ordered_canonical_coefficient_ids")
    families = [ids_family(x) for x in canonical]
    counts = {key: families.count(key) for key in ("BASE", "FRONT")}
    if counts != {"BASE": 27, "FRONT": 27}:
        raise ValueError(f"Expected 27 BASE and 27 FRONT entries, got {counts}")
    trans = np.asarray([1.0 if x == "BASE" else 0.0 for x in families])
    f = array(get(doc, "applied_force_coefficients_N"), "applied force", (54,))
    levels = levels_of(get(doc, "volume_levels"))
    if len(levels) != 2:
        raise ValueError(f"Expected two saved quadrature levels, got {len(levels)}")
    point_counts = {str(x["name"]): int(x["volume_point_count"]) for x in case["quadrature"]["volume_levels"]}
    if any(name not in point_counts for name, _ in levels) or sorted(point_counts[name] for name, _ in levels) != [128, 256]:
        raise ValueError(f"Saved quadrature metadata do not identify 128/256 levels: {point_counts}")
    off = get(doc, "offgrid")
    coords = array(off["query_coordinates_m"], "off-grid coordinates")
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError(f"Expected query coordinates Nx3, got {coords.shape}")
    u_ref = array(off["analytic_uz_m"], "analytic off-grid uz", (len(coords),))
    grad_ref = array(off["analytic_grad_uz_dimensionless"], "analytic off-grid gradient")
    if grad_ref.ndim != 2 or grad_ref.shape[0] != len(coords):
        raise ValueError(f"Expected analytic gradient NxD, got {grad_ref.shape}")
    level_results, saved_fields = [], {}
    by_name = off["results_by_volume_level"]
    for name, level in sorted(levels, key=lambda x: point_counts[x[0]]):
        K = array(level["stiffness_matrix_N_per_m"], f"{name} K", (54, 54))
        c = array(level["physical_mean_constraint_row_dimensionless"], f"{name} c", (54,))
        a = array(level["analytic_coefficient_vector_m"], f"{name} analytic coefficients", (54,))
        solution = level["solution"]
        d = array(solution["coefficients_m"], f"{name} d", (54,))
        lam = float(solution["lagrange_multiplier_N"])
        if not math.isfinite(lam):
            raise ValueError(f"{name} lambda non-finite")
        Kh = sha(np.asarray(K, dtype="<f8", order="C").tobytes(order="C"))
        eig = np.linalg.eigvalsh(0.5 * (K + K.T))
        residual = K @ d + c * lam - f
        constraint = float(c @ d)
        energy = float(0.5 * d @ K @ d)
        work = float(d @ f)
        tf, tc = float(trans @ f), float(trans @ c)
        compatibility = tf - lam * tc
        balance = 2 * energy + lam * constraint - work
        weak = array(level["analytic_coefficient_weak_defect_raw_N"], f"{name} analytic weak defect")
        field = by_name[name]
        u = array(field["solved_uz_m"], f"{name} solved uz", u_ref.shape)
        grad = array(field["solved_grad_uz_dimensionless"], f"{name} solved gradient", grad_ref.shape)
        du, dg = u - u_ref, grad - grad_ref
        dg_norm = np.linalg.norm(dg, axis=1)
        gr_norm = np.linalg.norm(grad_ref, axis=1)
        metrics = {
            "finite_sample_count": int(len(coords)),
            "uz_error_max_abs_m": float(np.max(np.abs(du))),
            "uz_error_rms_m": float(np.sqrt(np.mean(du * du))),
            "uz_reference_scaled_max_error": float(np.max(np.abs(du)) / np.max(np.abs(u_ref))) if np.max(np.abs(u_ref)) else None,
            "gradient_error_max_euclidean_dimensionless": float(np.max(dg_norm)),
            "gradient_error_rms_euclidean_dimensionless": float(np.sqrt(np.mean(dg_norm * dg_norm))),
            "gradient_reference_max_euclidean_dimensionless": float(np.max(gr_norm)),
            "gradient_error_max_reference_scaled": float(np.max(dg_norm) / np.max(gr_norm)) if np.max(gr_norm) else None,
        }
        residual_scale = np.linalg.norm(K @ d) + np.linalg.norm(c * lam) + np.linalg.norm(f)
        level_results.append({
            "name": name, "volume_quadrature_points": point_counts[name], "status": level.get("status"),
            "matrix": {"sha256_float64_c_order": Kh, "max_abs_asymmetry": float(np.max(np.abs(K-K.T))),
                       "eigenvalues_symmetric_part_sorted_N_per_m": eig.tolist(),
                       "raw_negative_eigenvalues": eig[eig < 0].tolist(), "minimum": float(eig[0]),
                       "second": float(eig[1]), "maximum": float(eig[-1])},
            "saddle": {"residual_vector_N": residual.tolist(), "residual_l2_N": float(np.linalg.norm(residual)),
                       "residual_max_abs_N": float(np.max(np.abs(residual))),
                       "residual_relative_to_term_norms": float(np.linalg.norm(residual)/residual_scale) if residual_scale else None,
                       "c_dot_d_m": constraint, "lambda_N": lam, "translation_dot_f_N": tf,
                       "translation_dot_c": tc, "translation_load_compatibility_N": compatibility,
                       "half_dKd_J": energy, "analytic_reference_energy_J": U_exact,
                       "energy_relative_error_percent": 100*(energy-U_exact)/U_exact,
                       "d_dot_f_J": work, "saddle_energy_work_defect_J": balance,
                       "analytic_coefficient_weak_defect_l2_N": float(np.linalg.norm(weak)),
                       "analytic_coefficient_weak_defect_max_abs_N": float(np.max(np.abs(weak)))},
            "offgrid_error": metrics})
        saved_fields[name] = (u, grad)

    high_name = max(levels, key=lambda item: point_counts[item[0]])[0]
    z0 = float(np.unique(coords[:, 2])[np.argmin(np.abs(np.unique(coords[:, 2])))])
    ix = np.flatnonzero(coords[:, 2] == z0)
    x, y = coords[ix, 0], coords[ix, 1]
    u0, _ = saved_fields[high_name]
    exact, solved = u_ref[ix], u0[ix]
    error = solved - exact
    good = np.isfinite(x) & np.isfinite(y) & np.isfinite(exact) & np.isfinite(solved)
    x, y, exact, solved, error = x[good], y[good], exact[good], solved[good], error[good]
    if not len(x):
        raise ValueError("No finite saved points on selected z layer")
    xlim = [float(coords[:, 0].min()), float(coords[:, 0].max())]
    ylim = [float(coords[:, 1].min()), float(coords[:, 1].max())]
    if xlim[0] == xlim[1]: xlim = [xlim[0]-1e-6, xlim[1]+1e-6]
    if ylim[0] == ylim[1]: ylim = [ylim[0]-1e-6, ylim[1]+1e-6]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.8), constrained_layout=True)
    vmax = max(float(np.max(np.abs(error))), np.finfo(float).tiny)
    panels = [(exact, "Analytic $u_z$", "$u_z$ (m)", "viridis", None),
              (solved, f"Saved solved $u_z$ ({high_name})", "$u_z$ (m)", "viridis", None),
              (error, "Solved − analytic", "$\\Delta u_z$ (m)", "coolwarm", Normalize(-vmax, vmax))]
    for ax, (values, title, label, cmap, norm_obj) in zip(axes, panels):
        sc = ax.scatter(x, y, c=values, s=29, cmap=cmap, norm=norm_obj, edgecolors="black", linewidths=0.25)
        if xlim[0] < 0: ax.plot([xlim[0], min(0.0, xlim[1])], [0, 0], "k--", linewidth=0.8)
        ax.set(xlim=xlim, ylim=ylim, xlabel="x (m)", ylabel="y (m)", title=f"{title}\nz={z0:.6g} m; n={len(values)}")
        ax.set_aspect("equal", adjustable="box")
        colorbar = fig.colorbar(sc, ax=ax, shrink=0.82)
        colorbar.set_label(label)
    fig.suptitle("Saved off-grid fields; scatter only, no interpolation across the crack branch")
    fig.savefig(OUT_PNG, dpi=170)
    plt.close(fig)

    after = hashes(monitored)
    checks = [{"path": p, "pre_sha256": h, "post_sha256": after.get(p), "unchanged": h == after.get(p)}
              for p, h in before.items()]
    pin_checks = []
    for raw_path, expected in pins.items():
        path = rooted_path(raw_path)
        if path is not None and str(path) in after:
            pin_checks.append({"path": str(path), "expected_sha256": expected,
                               "current_sha256": after[str(path)], "matches": after[str(path)] == expected})
    receipt_pin_checks = []
    for entry in receipt_pin_entries:
        path = rooted_path(entry["path"])
        if path is not None and str(path) in after:
            receipt_pin_checks.append({**entry, "current_sha256": after[str(path)],
                                       "matches": after[str(path)] == entry["sha256"]})
    receipt_claims = [{"path": p, "value": v} for p, v in walk(receipt)
                      if any(k in p.lower() for k in ("authorized_once", "observed", "promotion", "exit_code", "exitcode"))
                      and isinstance(v, (str, int, float, bool))]
    caveat = ("The primary run receipt's post-command authorized_once assignment had a missing-$ PowerShell error. "
              "Its observed flag is false; rc 0 is inferred from success/promotion acceptance, not directly measured by that snippet.")
    result = {"audit": {
        "inputs": {"primary_path": str(PRIMARY), "primary_bytes": len(primary_raw), "primary_sha256": primary_sha,
                    "primary_pin_matches": True, "input_receipt_path": str(INPUT_RECEIPT),
                    "input_receipt_bytes": len(receipt_raw), "input_receipt_sha256": sha(receipt_raw),
                    "source_freeze_path": str(FREEZE), "source_freeze_sha256": sha(freeze_raw),
                    "fixed_case_path": str(CASE), "fixed_case_sha256": sha(case_raw),
                    "base_fixed_case_path": str(base_case_path), "base_fixed_case_sha256": sha(base_case_raw),
                    "source_freeze_pin_count": len(pins), "source_pin_checks": pin_checks,
                    "all_source_pins_match": all(x["matches"] for x in pin_checks),
                    "run_receipt_pin_checks": receipt_pin_checks,
                    "all_resolvable_run_receipt_pins_match": all(x["matches"] for x in receipt_pin_checks),
                    "source_primary_receipt_prepost": checks,
                    "all_source_primary_receipt_hashes_unchanged": all(x["unchanged"] for x in checks)},
        "canonical_dofs": {"path": "registry.ordered_canonical_coefficient_ids", "count": 54,
                           "family_counts": counts, "translation_definition": "BASE=1, FRONT=0",
                           "translation_vector": trans.tolist()},
        "fixed_case_analytic_reference": {"parameters": {"E_Pa": E, "nu": nu, "KIII_Pa_sqrt_m": KIII, "R_m": R, "L_m": L},
                                            "paths": {"E": E_path, "nu": nu_path, "KIII": k_path, "R": r_path, "L": l_path},
                                            "mu_Pa": mu, "U_J": U_exact},
        "volume_levels": level_results,
        "plot": {"path": str(OUT_PNG), "level": high_name, "z_layer_m": z0, "point_count": len(x),
                 "physical_xy_limits_m": {"x": xlim, "y": ylim}, "interpolation": "none; saved samples plotted with scatter"},
        "input_run_receipt_scope": {"claims": receipt_claims, "caveat": caveat},
        "scope": ["Saved arrays only; no provider, assembly, integrals, or solver.",
                  "Small saddle residual is an algebraic check, not a physical field or fracture-growth qualification.",
                  "Quadrature-dependent errors are reported as sensitivities without a pass/fail threshold."],
        "elapsed_seconds": time.monotonic()-t0}}
    OUT_JSON.write_text(json.dumps(result, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    print(json.dumps({"exit_code": 0, "elapsed_seconds": result["audit"]["elapsed_seconds"],
                      "primary_sha256": primary_sha,
                      "levels": [{"name": q["name"], "K_sha256": q["matrix"]["sha256_float64_c_order"],
                                  "residual_l2_N": q["saddle"]["residual_l2_N"], "constraint_m": q["saddle"]["c_dot_d_m"],
                                  "lambda_N": q["saddle"]["lambda_N"], "translation_compatibility_N": q["saddle"]["translation_load_compatibility_N"],
                                  "energy_J": q["saddle"]["half_dKd_J"], "energy_error_percent": q["saddle"]["energy_relative_error_percent"],
                                  "grad_error_max_reference_scaled": q["offgrid_error"]["gradient_error_max_reference_scaled"]}
                                 for q in level_results],
                      "source_hashes_unchanged": result["audit"]["inputs"]["all_source_primary_receipt_hashes_unchanged"],
                      "source_pins_match": result["audit"]["inputs"]["all_source_pins_match"],
                      "json": str(OUT_JSON), "png": str(OUT_PNG)}, separators=(",", ":")))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
