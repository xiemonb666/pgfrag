#!/usr/bin/env python3
"""Single-pass read-only numerical audit of the saved Mode III operator."""

from __future__ import annotations

import hashlib
import json
import math
import re
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(r"C:\Users\xiemo\Desktop\ShardForge")
RUN_DIR = ROOT / "pgfrag/tmp/review_20261005/modeiii_pu_operator_20261005/runs/run01"
DIAG_PATH = RUN_DIR / "operator_diagnostic.json"
RECEIPT_PATH = RUN_DIR / "run_record.json"
SOURCE_DIR = ROOT / "pgfrag/experiments/modeiii_pu_operator_20261005"
OUT_DIR = ROOT / "pgfrag/tmp/review_20261005/modeiii_pu_operator_independent_20261005"
OUT_JSON = OUT_DIR / "audit_result.json"
OUT_MD = OUT_DIR / "audit_result.md"
EXPECTED_SOURCE_FREEZE = "65197acaead31faacfad967ffda94ebac9991ff7a433c86fc133e699ad591780"
EXPECTED_DIAGNOSTIC_SHA256 = "8a62d61384156921fecc585b2336f161dd572988eadae8bed46afadf105c22c2"


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def walk(value: Any, path: str = "$"):
    yield path, value
    if isinstance(value, dict):
        for k, v in value.items():
            yield from walk(v, f"{path}.{k}")
    elif isinstance(value, list):
        for i, v in enumerate(value):
            yield from walk(v, f"{path}[{i}]")


def norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", s.lower())


def list_numeric_arrays(obj: Any):
    for path, val in walk(obj):
        if not isinstance(val, list) or not val:
            continue
        try:
            arr = np.asarray(val, dtype=np.float64)
        except (TypeError, ValueError):
            continue
        if arr.size:
            yield path, val, arr


def matrix_candidates(obj: Any):
    for path, val, arr in list_numeric_arrays(obj):
        if arr.ndim == 2 and arr.shape == (54, 54):
            yield path, val, arr


def get_named_scalar(obj: Any, aliases: tuple[str, ...]):
    found = []
    for path, val in walk(obj):
        key = norm(path.rsplit(".", 1)[-1].split("[")[0])
        if key not in aliases:
            continue
        if isinstance(val, (int, float)) and not isinstance(val, bool):
            found.append((path, float(val)))
    return found


def find_prepost_pairs(obj: Any, path: str = "$"):
    pairs = []
    if isinstance(obj, dict):
        normalized = {norm(str(k)): v for k, v in obj.items()}
        before_keys = ("before", "pre", "sourcebefore", "sourcepre", "beforehash", "prehash", "sha256before", "sha256pre")
        after_keys = ("after", "post", "sourceafter", "sourcepost", "afterhash", "posthash", "sha256after", "sha256post")
        before = next((normalized[k] for k in before_keys if k in normalized), None)
        after = next((normalized[k] for k in after_keys if k in normalized), None)
        if before is not None or after is not None:
            pairs.append({"path": path, "before": before, "after": after, "equal": before == after})
        for k, v in obj.items():
            pairs.extend(find_prepost_pairs(v, f"{path}.{k}"))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            pairs.extend(find_prepost_pairs(v, f"{path}[{i}]"))
    return pairs


def main() -> int:
    started = time.monotonic()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    diag_raw = DIAG_PATH.read_bytes()
    receipt_raw = RECEIPT_PATH.read_bytes()
    diag = json.loads(diag_raw)
    receipt = json.loads(receipt_raw)

    mats = list(matrix_candidates(diag))
    mat_rows = [{"path": p, "shape": list(a.shape)} for p, _, a in mats]
    selected = None
    if mats:
        # Prefer a field explicitly named as an operator, stiffness, or K.
        selected = sorted(
            mats,
            key=lambda x: (
                not any(t in norm(x[0]) for t in ("operator", "stiffness", "matrix", "k")),
                len(x[0]),
            ),
        )[0]
    if selected is None:
        raise RuntimeError("No saved numeric 54x54 matrix found in the diagnostic JSON")

    matrix_path, _, raw_matrix = selected
    K = np.array(raw_matrix, dtype=np.float64, order="C", copy=True)
    canonical_bytes = np.asarray(K, dtype="<f8", order="C").tobytes(order="C")
    matrix_hash = sha256(canonical_bytes)
    sym_abs = float(np.max(np.abs(K - K.T)))
    sym_fro = float(np.linalg.norm(K - K.T, ord="fro"))
    norm_fro = float(np.linalg.norm(K, ord="fro"))
    eigenvalues = np.linalg.eigvalsh(K)
    positive = eigenvalues[eigenvalues > 0.0]
    quotient_cond = float(positive[-1] / positive[0]) if positive.size else None
    t = np.ones(K.shape[0], dtype=np.float64)
    kt = K @ t
    trans_residual_l2 = float(np.linalg.norm(kt))
    trans_scale = float(np.linalg.norm(K, ord=2) * np.linalg.norm(t))
    trans_rel = trans_residual_l2 / trans_scale if trans_scale else None

    vector_candidates = []
    semantic_candidates = []
    all_array_shapes = []
    for path, val, arr in list_numeric_arrays(diag):
        all_array_shapes.append({"path": path, "shape": list(arr.shape)})
        leaf = norm(path.rsplit(".", 1)[-1].split("[")[0])
        if arr.size == 54 and any(s in leaf for s in ("dstar", "displacement", "prescribed", "target", "dofvalues")):
            vector_candidates.append((path, arr.reshape(-1)))
    for path, val in walk(diag):
        leaf = norm(path.rsplit(".", 1)[-1].split("[")[0])
        if "semantic" in leaf and isinstance(val, list):
            semantic_candidates.append({"path": path, "count": len(val), "sample": val[:6]})

    # Prefer the candidate whose field name contains d-star.
    dstar = None
    if vector_candidates:
        vector_candidates.sort(key=lambda x: ("dstar" not in norm(x[0]), len(x[0])))
        dstar_path, dstar = vector_candidates[0]
    else:
        dstar_path = None
    quadratic_energy = None
    if dstar is not None:
        quadratic_energy = float(0.5 * dstar @ K @ dstar)

    aliases = {
        "E": ("e", "epa", "youngsmodulus", "elasticmodulus"),
        "nu": ("nu", "poissonratio"),
        "KIII": ("kiii", "kiiipa sqrtm", "stressintensity", "stressintensityfactor"),
        "R": ("r", "radius", "cylinderradius", "cylinderr"),
        "L": ("l", "length", "frontlength", "cylinderlength"),
        "mu": ("mu", "mupa", "shearmodulus"),
        "U": ("u", "uj", "energy", "analyticenergy", "continuumenergy", "energyexact"),
    }
    scalars = {key: get_named_scalar(diag, value) for key, value in aliases.items()}
    E, nu, KIII, R, L = 1.0e6, 0.3, 100.0, 0.01, 0.01
    mu_exact = E / (2.0 * (1.0 + nu))
    U_exact = KIII * KIII * R * L / (2.0 * mu_exact)
    algebraic_metadata_check = {
        "prescribed_values": {"E": E, "nu": nu, "KIII": KIII, "R": R, "L": L},
        "exact_mu": mu_exact,
        "exact_U": U_exact,
        "observations": {},
    }
    for name in ("E", "nu", "KIII", "R", "L"):
        algebraic_metadata_check["observations"][name] = [
            {"path": p, "value": v, "difference_from_prescribed": v - {"E": E, "nu": nu, "KIII": KIII, "R": R, "L": L}[name]}
            for p, v in scalars[name]
        ]
    algebraic_metadata_check["observations"]["mu"] = [
        {"path": p, "value": v, "difference_from_exact_formula": v - mu_exact}
        for p, v in scalars["mu"]
    ]
    algebraic_metadata_check["observations"]["U"] = [
        {"path": p, "value": v, "difference_from_exact_formula": v - U_exact}
        for p, v in scalars["U"]
    ]

    source_hashes = {}
    for f in sorted(SOURCE_DIR.rglob("*")):
        if f.is_file():
            rel = f.relative_to(SOURCE_DIR).as_posix()
            source_hashes[rel] = sha256(f.read_bytes())

    # Surface any receipt/diagnostic digest and source provenance fields for direct comparison.
    provenance_fields = []
    for label, document in (("diagnostic", diag), ("receipt", receipt)):
        for path, val in walk(document):
            path_l = path.lower()
            if any(t in path_l for t in ("source", "pre", "post", "freeze", "receipt", "diagnostic", "hash", "sha")):
                if isinstance(val, (str, int, float, bool)) or (isinstance(val, list) and len(val) <= 20):
                    provenance_fields.append({"document": label, "path": path, "value": val})

    result = {
        "audit": {
            "inputs": {
                "diagnostic_path": str(DIAG_PATH),
                "diagnostic_bytes": len(diag_raw),
                "diagnostic_sha256": sha256(diag_raw),
            "diagnostic_sha256_matches_supplied_pin": sha256(diag_raw) == EXPECTED_DIAGNOSTIC_SHA256,
            "supplied_diagnostic_sha256": EXPECTED_DIAGNOSTIC_SHA256,
                "receipt_path": str(RECEIPT_PATH),
                "receipt_bytes": len(receipt_raw),
                "receipt_sha256": sha256(receipt_raw),
                "expected_source_freeze_sha256": EXPECTED_SOURCE_FREEZE,
            "current_source_file_count": len(source_hashes),
                "source_file_hashes_current": source_hashes,
            },
            "matrix": {
                "selected_path": matrix_path,
                "shape": list(K.shape),
                "dtype": "float64 little-endian C-order",
                "sha256_float64_c": matrix_hash,
                "candidate_paths": mat_rows,
                "max_abs_asymmetry": sym_abs,
                "fro_asymmetry": sym_fro,
                "fro_norm": norm_fro,
                "eigenvalues_full_sorted": eigenvalues.tolist(),
                "positive_quotient_condition_number": quotient_cond,
                "uniform_translation_residual_l2": trans_residual_l2,
                "uniform_translation_relative_to_spectral_scale": trans_rel,
            },
            "prescribed_data": {
                "semantic_id_candidates": semantic_candidates,
                "dstar_path": dstar_path,
                "dstar_count": int(dstar.size) if dstar is not None else None,
                "dstar_values": dstar.tolist() if dstar is not None else None,
                "quadratic_energy_half_dstar_K_dstar": quadratic_energy,
                "numeric_array_shapes": all_array_shapes,
            },
            "analytic_reference": {
                "E_Pa": E,
                "nu": nu,
                "KIII_Pa_sqrt_m": KIII,
                "cylinder_R_m": R,
                "front_length_L_m": L,
                "mu_exact_Pa": mu_exact,
                "U_exact_J": U_exact,
                "formulas": {
                    "mu": "E/(2*(1+nu))",
                    "U": "KIII^2*R*L/(2*mu), integrated Mode III singular-field energy in cylindrical radius-R domain",
                },
                "metadata_numeric_candidates": scalars,
                "metadata_algebraic_consistency": algebraic_metadata_check,
            },
            "provenance_fields": provenance_fields,
            "source_prepost_pairs": {
                "diagnostic": find_prepost_pairs(diag),
                "receipt": find_prepost_pairs(receipt),
            },
            "elapsed_seconds": time.monotonic() - started,
            "limitations": [
                "Read-only audit of the saved matrix and metadata; no operator reassembly, solver, experimental provider, or field rendering.",
                "Numerical residuals and eigenvalues are reported as measurements; no rounding tolerance or PSD repair is applied.",
                "The audit scope is the fixed continuous MLS plus normalized unshifted F1 PU scalar-Uz operator with T=all 27 and J empty; it does not establish propagation, local qualification, or full 3D behavior.",
            ],
        }
    }
    OUT_JSON.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    m = result["audit"]["matrix"]
    lines = [
        "# Independent Mode III operator audit",
        "",
        f"- Diagnostic SHA-256: `{result['audit']['inputs']['diagnostic_sha256']}` ({result['audit']['inputs']['diagnostic_bytes']} bytes).",
        f"- Receipt SHA-256: `{result['audit']['inputs']['receipt_sha256']}` ({result['audit']['inputs']['receipt_bytes']} bytes).",
        f"- Matrix: `{matrix_path}`, shape 54×54, float64 C-order SHA-256 `{matrix_hash}`.",
        f"- Max absolute asymmetry: `{sym_abs:.17g}`; Frobenius asymmetry: `{sym_fro:.17g}`.",
        f"- Uniform translation residual: `||K·1||₂={trans_residual_l2:.17g}`; relative to spectral scale `{trans_rel}`.",
        f"- Positive-quotient condition number: `{quotient_cond}`.",
        f"- Full sorted eigenvalues: `{json.dumps(eigenvalues.tolist(), separators=(',', ':'))}`.",
        f"- d* source: `{dstar_path}`; count `{int(dstar.size) if dstar is not None else 'unavailable'}`; ½d*ᵀKd* `{quadratic_energy}`.",
        f"- Analytic references: μ = `{mu_exact:.17g} Pa`; U = `{U_exact:.17g} J` using U = KIII²RL/(2μ).",
        f"- Elapsed: `{result['audit']['elapsed_seconds']:.3f} s`.",
        "",
        "See `audit_result.json` for semantic-ID candidates, metadata values, provenance fields, and current frozen-source file hashes.",
        "",
        "Scope: saved fixed-operator readback only; no operator reassembly, solver, experimental provider, or field rendering. No numerical acceptance threshold or PSD repair was applied. This does not establish propagation, local qualification, or full 3D behavior.",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps({
        "exit_code": 0,
        "elapsed_seconds": result["audit"]["elapsed_seconds"],
        "diagnostic_sha256": result["audit"]["inputs"]["diagnostic_sha256"],
        "receipt_sha256": result["audit"]["inputs"]["receipt_sha256"],
        "diagnostic_sha256_matches_supplied_pin": result["audit"]["inputs"]["diagnostic_sha256_matches_supplied_pin"],
        "current_source_file_count": result["audit"]["inputs"]["current_source_file_count"],
        "matrix_sha256_float64_c": matrix_hash,
        "matrix_path": matrix_path,
        "matrix_shape": list(K.shape),
        "max_abs_asymmetry": sym_abs,
        "eigenvalues": eigenvalues.tolist(),
        "positive_quotient_condition_number": quotient_cond,
        "translation_residual_l2": trans_residual_l2,
        "dstar_path": dstar_path,
        "quadratic_energy": quadratic_energy,
        "mu_exact_Pa": mu_exact,
        "U_exact_J": U_exact,
        "output_json": str(OUT_JSON),
        "output_md": str(OUT_MD),
    }, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(json.dumps({"exit_code": 1, "error": repr(exc)}, separators=(",", ":")), file=sys.stderr)
        raise
