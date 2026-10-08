#!/usr/bin/env python3
"""Semantic-coordinate audit of the already-saved Mode III operator.

This script reads only the saved diagnostic and its run receipt. It does not
assemble an operator, call a solver/provider, modify frozen inputs, or render a
field. It writes its report beside this script when explicitly run.
"""

from __future__ import annotations

import hashlib
import json
import re
import time
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(r"C:\Users\xiemo\Desktop\ShardForge")
RUN_DIR = ROOT / "pgfrag/tmp/review_20261005/modeiii_pu_operator_20261005/runs/run01"
DIAG_PATH = RUN_DIR / "operator_diagnostic.json"
RECEIPT_PATH = RUN_DIR / "run_record.json"
OUT_DIR = ROOT / "pgfrag/tmp/review_20261005/modeiii_pu_operator_independent_20261005"
OUT_JSON = OUT_DIR / "audit_result_v2.json"
EXPECTED_DIAGNOSTIC_SHA256 = "8a62d61384156921fecc585b2336f161dd572988eadae8bed46afadf105c22c2"
EXPECTED_MATRIX_SHA256 = "bb229305198d99b7aa654736eec0aec7d4197504be3f59fbecbaf12fbb368729"

# Preserve the first pass exactly. Its scientific output remains in audit_result.json;
# this records the runtime variation instead of rewriting that receipt.
V1_RUNTIME_RECEIPT = {
    "initial_alias_launch": {
        "argv": ["python", "-B", "pgfrag/tmp/review_20261005/modeiii_pu_operator_independent_20261005/audit_modeiii.py"],
        "exit_code": 1,
        "wall_seconds": 0.3866498,
        "result": "Windows Store Python alias; no inputs were opened",
    },
    "successful_v1_launch": {
        "argv": [
            "C:\\Users\\xiemo\\.cache\\codex-runtimes\\codex-primary-runtime\\dependencies\\python\\python.exe",
            "-B",
            "pgfrag/tmp/review_20261005/modeiii_pu_operator_independent_20261005/audit_modeiii.py",
        ],
        "exit_code": 0,
        "wall_seconds": 0.4366498,
        "script_measured_seconds": 0.0940000000409782,
        "runtime_variation": "app-bundled Python was used; the authorized VesselForge venv was not used",
        "v1_result_path": "pgfrag/tmp/review_20261005/modeiii_pu_operator_independent_20261005/audit_result.json",
        "v1_receipt_sha256": "101f948ea40c40b63845c60df417efad84343d9fd41ad2ebf875a83cbab85f2a",
    },
}


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def norm(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def get_path(document: dict[str, Any], dotted: str) -> Any:
    value: Any = document
    for part in dotted.split("."):
        if not isinstance(value, dict) or part not in value:
            raise KeyError(dotted)
        value = value[part]
    return value


def first_available(documents: list[tuple[str, dict[str, Any]]], dotted_paths: list[str]):
    for document_name, document in documents:
        for dotted in dotted_paths:
            try:
                return document_name, dotted, get_path(document, dotted)
            except KeyError:
                continue
    raise KeyError(f"None of the paths exists: {dotted_paths}")


def flatten_canonical_ids(value: Any) -> list[Any]:
    """Preserve the JSON canonical order and carry map keys as semantic labels."""
    if isinstance(value, list):
        return value
    if isinstance(value, dict):
        flattened = []
        for key, child in value.items():
            if isinstance(child, list):
                flattened.extend({"semantic_key": key, "canonical_id": item} for item in child)
            else:
                flattened.append({"semantic_key": key, "canonical_id": child})
        return flattened
    raise TypeError(f"registry.canonical_ids must be a JSON list or object, got {type(value).__name__}")


def semantic_text(item: Any) -> str:
    if isinstance(item, dict):
        preferred = []
        fallback = []
        semantic_names = {
            "semantickey", "semanticid", "family", "basis", "role", "coefficientkind",
            "canonicalid", "id", "label", "name", "dof",
        }
        for key, value in item.items():
            if isinstance(value, (str, int, float)):
                target = preferred if norm(key) in semantic_names else fallback
                target.append(f"{key}={value}")
        return " | ".join(preferred + fallback + [json.dumps(item, sort_keys=True, ensure_ascii=False)])
    return str(item)


def family_for(item: Any) -> str:
    text = semantic_text(item).upper()
    has_base = re.search(r"(?<![A-Z0-9])BASE(?![A-Z0-9])", text) is not None
    has_front = re.search(r"(?<![A-Z0-9])FRONT(?![A-Z0-9])", text) is not None
    # Some canonical key schemes use an explicit B/F prefix separated from the
    # source identifier; accept those only as a final fallback.
    if not has_base and not has_front:
        has_base = re.search(r"(?:^|[|:=_\- ])B(?:$|[|:=_\- ])", text) is not None
        has_front = re.search(r"(?:^|[|:=_\- ])F(?:$|[|:=_\- ])", text) is not None
    if has_base == has_front:
        raise ValueError(f"Cannot unambiguously classify canonical ID by BASE/FRONT semantic key: {item!r}")
    return "BASE" if has_base else "FRONT"


def main() -> int:
    started = time.monotonic()
    if OUT_JSON.exists():
        raise FileExistsError(f"Refusing to overwrite existing v2 audit result: {OUT_JSON}")
    diagnostic_raw = DIAG_PATH.read_bytes()
    diagnostic_sha256 = sha256(diagnostic_raw)
    if diagnostic_sha256 != EXPECTED_DIAGNOSTIC_SHA256:
        raise ValueError(
            f"Diagnostic SHA-256 mismatch before matrix parse: expected {EXPECTED_DIAGNOSTIC_SHA256}, got {diagnostic_sha256}"
        )
    diagnostic = json.loads(diagnostic_raw)
    receipt_raw = RECEIPT_PATH.read_bytes()
    receipt = json.loads(receipt_raw)

    # These exact frozen metadata paths define this fixed case. Do not fall
    # back to guessed aliases or prompt defaults if an input is missing.
    exact_parameter_paths = {
        "E_Pa": "fixed_case.material_and_manufactured_field.young_modulus_Pa",
        "nu": "fixed_case.material_and_manufactured_field.poisson_ratio",
        "KIII_Pa_sqrt_m": "fixed_case.material_and_manufactured_field.K_III_Pa_sqrt_m",
        "R_m": "fixed_case.geometry.radius_m",
        "L_m": "fixed_case.geometry.length_m",
    }
    exact_parameters = {}
    for name, path in exact_parameter_paths.items():
        value = float(get_path(diagnostic, path))
        if not np.isfinite(value):
            raise ValueError(f"Non-finite fixed-case value at $.{path}: {value}")
        exact_parameters[name] = value
    z_bounds = np.asarray(
        get_path(diagnostic, "fixed_case.geometry.z_bounds_m"), dtype=np.float64
    ).reshape(-1)
    if z_bounds.shape != (2,) or not np.all(np.isfinite(z_bounds)) or z_bounds[1] <= z_bounds[0]:
        raise ValueError(f"Bad fixed-case z_bounds_m: {z_bounds.tolist()}")
    z_span = float(z_bounds[1] - z_bounds[0])
    if exact_parameters["L_m"] != z_span:
        raise ValueError(
            f"geometry.length_m disagrees with z_bounds_m span: {exact_parameters['L_m']} != {z_span}"
        )
    if exact_parameters["E_Pa"] <= 0 or not (-1.0 < exact_parameters["nu"] < 0.5):
        raise ValueError(f"Bad isotropic material parameters: {exact_parameters}")
    if exact_parameters["R_m"] <= 0 or exact_parameters["L_m"] <= 0:
        raise ValueError(f"Bad cylinder dimensions: {exact_parameters}")
    expected_parameters = {
        "E_Pa": 1.0e6,
        "nu": 0.3,
        "KIII_Pa_sqrt_m": 100.0,
        "R_m": 0.01,
        "L_m": 0.01,
    }
    for name, expected in expected_parameters.items():
        if exact_parameters[name] != expected:
            raise ValueError(f"Unexpected fixed-case {name}: {exact_parameters[name]} != {expected}")

    K = np.asarray(get_path(diagnostic, "operator_diagnostics.stiffness_matrix_N_per_m"), dtype=np.float64)
    if K.shape != (54, 54):
        raise ValueError(f"Expected saved 54x54 K, got {K.shape}")
    K = np.array(K, dtype=np.float64, order="C", copy=True)
    K_hash = sha256(np.asarray(K, dtype="<f8", order="C").tobytes(order="C"))

    canonical_raw = get_path(diagnostic, "registry.canonical_ids")
    canonical_ids = flatten_canonical_ids(canonical_raw)
    if len(canonical_ids) != 54:
        raise ValueError(f"Expected 54 registry.canonical_ids in matrix order, found {len(canonical_ids)}")
    families = [family_for(item) for item in canonical_ids]
    counts = {family: families.count(family) for family in ("BASE", "FRONT")}
    if counts != {"BASE": 27, "FRONT": 27}:
        raise ValueError(f"Expected 27 BASE and 27 FRONT coefficients, got {counts}")
    translation = np.asarray([1.0 if family == "BASE" else 0.0 for family in families], dtype=np.float64)

    dstar = np.asarray(
        get_path(diagnostic, "basis_and_coefficients.manufactured_coefficient_vector"),
        dtype=np.float64,
    ).reshape(-1)
    if dstar.shape != (54,):
        raise ValueError(f"Expected manufactured coefficient vector with 54 entries, got {dstar.shape}")
    quadratic_energy = float(0.5 * dstar @ K @ dstar)

    # The recorded spectrum states that eigvalsh used the symmetric part only;
    # preserve the assembled K above without modifying or repairing it.
    K_symmetric_for_spectrum = 0.5 * (K + K.T)
    eigenvalues = np.linalg.eigvalsh(K_symmetric_for_spectrum)
    spectrum_owner, spectrum_path, stored_eigenvalues_raw = first_available(
        [("diagnostic", diagnostic), ("receipt", receipt)],
        [
            "operator_diagnostics.eigenvalues_N_per_m",
            "primary_summary.operator_diagnostics.eigenvalues_N_per_m",
        ],
    )
    stored_eigenvalues = np.asarray(stored_eigenvalues_raw, dtype=np.float64).reshape(-1)
    if stored_eigenvalues.shape != eigenvalues.shape:
        raise ValueError(f"Stored spectrum shape {stored_eigenvalues.shape} != recomputed {eigenvalues.shape}")
    eig_diff = eigenvalues - stored_eigenvalues

    threshold_owner, threshold_path, threshold_raw = first_available(
        [("diagnostic", diagnostic), ("receipt", receipt)],
        [
            "operator_diagnostics.near_zero_absolute_threshold_N_per_m",
            "primary_summary.operator_diagnostics.near_zero_absolute_threshold_N_per_m",
        ],
    )
    threshold = float(threshold_raw)
    near_zero = np.abs(eigenvalues) <= threshold
    positive_quotient = eigenvalues[eigenvalues > threshold]
    quotient_condition = (
        float(positive_quotient[-1] / positive_quotient[0]) if positive_quotient.size else None
    )
    raw_negative = eigenvalues[eigenvalues < 0.0]
    mode_count_owner, mode_count_path, mode_count_raw = first_available(
        [("diagnostic", diagnostic), ("receipt", receipt)],
        [
            "operator_diagnostics.expected_scalar_rigid_mode_count",
            "primary_summary.operator_diagnostics.expected_scalar_rigid_mode_count",
        ],
    )
    documented_modes = int(mode_count_raw)

    # Recompute analytic references directly from the exact extracted inputs.
    E = exact_parameters["E_Pa"]
    nu = exact_parameters["nu"]
    KIII = exact_parameters["KIII_Pa_sqrt_m"]
    R = exact_parameters["R_m"]
    L = exact_parameters["L_m"]
    mu_exact = E / (2.0 * (1.0 + nu))
    U_exact = KIII * KIII * R * L / (2.0 * mu_exact)
    metadata_observed = {
        "expected_fixed_case_values": expected_parameters,
        "formula_inputs_extracted_by_exact_path": exact_parameters,
        "exact_formula_values": {"mu_Pa": mu_exact, "U_J": U_exact},
        "z_bounds_m": z_bounds.tolist(),
        "z_span_m": z_span,
        "length_matches_z_span_exactly": exact_parameters["L_m"] == z_span,
        "formulae": {
            "mu_Pa": "E/(2*(1+nu))",
            "U_J": "KIII^2*R*L/(2*mu)",
        },
    }

    translation_residual = K @ translation
    translation_l2 = float(np.linalg.norm(translation_residual))
    translation_scale = float(np.linalg.norm(K, ord=2) * np.linalg.norm(translation))
    translation_relative = translation_l2 / translation_scale if translation_scale else None
    asymmetry_abs = float(np.max(np.abs(K - K.T)))

    result = {
        "audit": {
            "v2_runtime": {
                "requested_interpreter": "C:/Users/xiemo/Desktop/ShardForge/third_party/vesselforge-runtime/venv/Scripts/python.exe",
                "argv_suffix": ["-B", str((OUT_DIR / "audit_modeiii_semantic_v2.py").resolve())],
                "must_complete_within_seconds": 30,
                "v1_runtime_receipt": V1_RUNTIME_RECEIPT,
            },
            "inputs": {
                "diagnostic_path": str(DIAG_PATH),
                "diagnostic_bytes": len(diagnostic_raw),
                "diagnostic_sha256": sha256(diagnostic_raw),
                "diagnostic_sha256_matches_supplied_pin": sha256(diagnostic_raw) == EXPECTED_DIAGNOSTIC_SHA256,
                "receipt_path": str(RECEIPT_PATH),
                "receipt_bytes": len(receipt_raw),
                "receipt_sha256": sha256(receipt_raw),
                "receipt_primary_summary_exists": isinstance(receipt.get("primary_summary"), dict),
            },
            "matrix": {
                "shape": list(K.shape),
                "matrix_path": "$.operator_diagnostics.stiffness_matrix_N_per_m",
                "sha256_float64_c_order": K_hash,
                "sha256_matches_v1": K_hash == EXPECTED_MATRIX_SHA256,
                "max_abs_asymmetry": asymmetry_abs,
                "eigenvalue_procedure": "eigvalsh((K+K.T)/2), as specified by saved diagnostic metadata; raw K unchanged",
                "stored_spectrum_source": {"document": spectrum_owner, "path": spectrum_path},
                "recomputed_eigenvalues_sorted": eigenvalues.tolist(),
                "stored_eigenvalues_sorted": stored_eigenvalues.tolist(),
                "eigenvalue_difference_recomputed_minus_stored": eig_diff.tolist(),
                "max_abs_eigenvalue_difference": float(np.max(np.abs(eig_diff))),
                "raw_negative_eigenvalues": raw_negative.tolist(),
                "raw_negative_eigenvalue_count": int(raw_negative.size),
                "near_zero_absolute_threshold_N_per_m": threshold,
                "near_zero_threshold_source": {"document": threshold_owner, "path": threshold_path},
                "near_zero_classification_abs_lambda_le_threshold": near_zero.tolist(),
                "near_zero_count_at_reported_threshold": int(np.count_nonzero(near_zero)),
                "reported_scalar_rigid_mode_count": documented_modes,
                "reported_scalar_rigid_mode_count_source": {"document": mode_count_owner, "path": mode_count_path},
                "positive_quotient_condition_number_excluding_reported_near_zero_modes": quotient_condition,
                "positive_quotient_min_eigenvalue_N_per_m": float(positive_quotient[0]) if positive_quotient.size else None,
                "largest_eigenvalue_N_per_m": float(eigenvalues[-1]),
            },
            "semantic_dofs": {
                "source_path": "$.registry.canonical_ids",
                "canonical_ids_in_matrix_order": [semantic_text(item) for item in canonical_ids],
                "basis_families_in_matrix_order": families,
                "family_counts": counts,
                "uniform_physical_scalar_translation_definition": "BASE coefficients=1, FRONT coefficients=0",
                "translation_vector": translation.tolist(),
                "translation_residual_Kt_l2": translation_l2,
                "translation_residual_relative_to_spectral_scale": translation_relative,
            },
            "manufactured_vector_energy": {
                "vector_path": "$.basis_and_coefficients.manufactured_coefficient_vector",
                "vector_count": int(dstar.size),
                "vector": dstar.tolist(),
                "quadratic_energy_half_dstar_K_dstar_J": quadratic_energy,
            },
            "analytic_reference_and_metadata": metadata_observed,
            "scope": [
                "Saved-matrix algebraic audit only; no operator reassembly, solver, provider import, or field rendering.",
                "Raw eigenvalue signs and the saved near-zero threshold are reported without treating the threshold as scientific acceptance; no PSD repair was applied.",
                "This case is the prescribed straight-crack all-source F1 PU scalar-Uz algebraic diagnostic, not an uncracked model or fracture-growth/local-qualification/full-3D result.",
            ],
            "elapsed_seconds": time.monotonic() - started,
        }
    }
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "exit_code": 0,
        "elapsed_seconds": result["audit"]["elapsed_seconds"],
        "diagnostic_sha256_matches_supplied_pin": result["audit"]["inputs"]["diagnostic_sha256_matches_supplied_pin"],
        "matrix_sha256": K_hash,
        "eigenvalue_count": int(eigenvalues.size),
        "raw_negative_eigenvalues": raw_negative.tolist(),
        "near_zero_count": int(np.count_nonzero(near_zero)),
        "quotient_condition_number": quotient_condition,
        "translation_residual_l2": translation_l2,
        "quadratic_energy_J": quadratic_energy,
        "mu_exact_Pa": mu_exact,
        "U_exact_J": U_exact,
        "output_json": str(OUT_JSON),
    }, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
