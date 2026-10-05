"""One bounded MLS-to-elastic-operator diagnostic; not a solver or benchmark.

When reviewed and intentionally invoked, this script evaluates actual MLS
shape values/gradients on the fixed SI cube, calls the elastic operator once,
and writes one JSON record. It never loads a crack surface, applies a load,
solves a system, tunes support, or repairs a matrix. Diagnostic thresholds
are roundoff references only; no numerical metric is promoted to an acceptance
gate.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
PGFRAG = HERE.parents[1]
MLS_DIR = PGFRAG / "experiments" / "mls3d_kernel_20261005"
OPERATOR_DIR = PGFRAG / "experiments" / "efgm_elastic_operator_20261005"


CASE_PATH = HERE / "fixed_case.json"
CASE_ID = "efgm_mls_patch_si_cube_20261005"
MAX_ARRAY_BUDGET = 32 * 1024 * 1024
WALL_BUDGET_SECONDS = 30.0


class DiagnosticBudgetError(RuntimeError):
    """Raised when a preflight array estimate or wall-time budget is exceeded."""


def _load_numerical_dependencies() -> None:
    """Load numerical modules only in the timed worker process."""
    global np, AssemblyStatus, ElasticOperatorNumericsError, IsotropicMaterial
    global rigid_mode_coefficient_vectors, assemble_elastic_operator
    global MLS3DConfig, MLSStatus, MLSShapeProvider

    for path in (str(MLS_DIR), str(OPERATOR_DIR), str(HERE)):
        if path not in sys.path:
            sys.path.insert(0, path)
    import numpy as np
    from efgm_elastic_operator import (
        AssemblyStatus,
        ElasticOperatorNumericsError,
        IsotropicMaterial,
        rigid_mode_coefficient_vectors,
        assemble_elastic_operator,
    )
    from mls3d_mls import MLS3DConfig, MLSStatus
    from mls_elastic_adapter import MLSShapeProvider


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _timestamp() -> str:
    return datetime.now(timezone.utc).isoformat()


def _require_fixed_case(case: dict[str, object]) -> None:
    if case.get("case_id") != CASE_ID:
        raise ValueError("fixed_case.json case_id does not match this source")
    domain = case["domain"]
    material = case["material"]
    mls = case["mls"]
    quadrature = case["volume_quadrature"]
    diagnostic = case["diagnostic"]
    assert isinstance(domain, dict) and isinstance(material, dict)
    assert isinstance(mls, dict) and isinstance(quadrature, dict)
    assert isinstance(diagnostic, dict)
    expected = (
        domain["bounds_m"] == [[-1.0, 1.0], [-1.0, 1.0], [-1.0, 1.0]],
        domain["coefficient_grid_shape"] == [3, 3, 3],
        domain["coefficient_axis_coordinates_m"] == [-1.0, 0.0, 1.0],
        domain["coefficient_count"] == 27,
        domain["displacement_dof_count"] == 81,
        quadrature["points_per_axis"] == 4,
        quadrature["point_count"] == 64,
        len(quadrature["gauss_legendre_nodes_on_minus1_plus1"]) == 4,
        len(quadrature["gauss_legendre_weights"]) == 4,
        quadrature["total_volume_m3"] == 8.0,
        material["young_modulus_Pa"] == 1_000_000.0,
        material["poisson_ratio"] == 0.3,
        mls["support_radius_m"] == 2.5,
        mls["reference_scale_m"] == 1.0,
        mls["reference_origin_m"] == [0.0, 0.0, 0.0],
        mls["visibility_mask"] is None,
        mls["mask_locally_constant"] is None,
        diagnostic["max_estimated_numeric_array_bytes"] == MAX_ARRAY_BUDGET,
        diagnostic["wall_timeout_seconds"] == WALL_BUDGET_SECONDS,
        diagnostic["max_displacement_dofs"] == 375,
        diagnostic["max_quadrature_points"] == 256,
    )
    if not all(expected):
        raise ValueError("fixed case changed a frozen geometry/material/MLS/budget value")


def _build_fixture(case: dict[str, object]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    domain = case["domain"]
    quadrature = case["volume_quadrature"]
    assert isinstance(domain, dict) and isinstance(quadrature, dict)
    axis = np.asarray(domain["coefficient_axis_coordinates_m"], dtype=np.float64)
    xx, yy, zz = np.meshgrid(axis, axis, axis, indexing="ij")
    source_positions = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))

    nodes_1d = np.asarray(quadrature["gauss_legendre_nodes_on_minus1_plus1"], dtype=np.float64)
    weights_1d = np.asarray(quadrature["gauss_legendre_weights"], dtype=np.float64)
    bounds = np.asarray(domain["bounds_m"], dtype=np.float64)
    midpoint = np.mean(bounds, axis=1)
    half_extent = 0.5 * (bounds[:, 1] - bounds[:, 0])
    qpoints_list: list[list[float]] = []
    volume_weights_list: list[float] = []
    for index in itertools.product(range(4), repeat=3):
        reference = nodes_1d[np.asarray(index)]
        qpoints_list.append((midpoint + half_extent * reference).tolist())
        volume_weights_list.append(float(np.prod(weights_1d[np.asarray(index)] * half_extent)))
    query_positions = np.asarray(qpoints_list, dtype=np.float64)
    volume_weights = np.asarray(volume_weights_list, dtype=np.float64)
    return source_positions, query_positions, volume_weights


def _array_budget_estimate(n_source: int, n_query: int) -> dict[str, int]:
    f8 = np.dtype(np.float64).itemsize
    i8 = np.dtype(np.int64).itemsize
    dof_count = 3 * n_source
    matrix_bytes = dof_count * dof_count * f8
    explicit = {
        "source_positions": n_source * 3 * f8,
        "quadrature_positions": n_query * 3 * f8,
        "quadrature_volume_weights": n_query * f8,
        "coefficient_displacements": n_source * 3 * f8,
        "mls_weights": n_query * n_source * f8,
        "mls_gradients": n_query * n_source * 3 * f8,
        "mls_ids_and_per_query_operator_copies": n_query * n_source * (i8 + f8 + 3 * f8 + i8),
        "operator_global_local_and_spectral_work_arrays": 8 * matrix_bytes,
        "rigid_modes_and_interpolated_fields": dof_count * 6 * f8 + n_query * (3 + 9) * f8,
    }
    return {**explicit, "estimated_explicit_array_bytes": int(sum(explicit.values()))}


def _roundoff_references(
    n_query: int, n_dof: int, maximum_condition: float
) -> dict[str, float | int | str]:
    epsilon = float(np.finfo(np.float64).eps)
    unit_roundoff = epsilon / 2.0
    # Two length-six contractions form each local B.T D B entry, n_query
    # additions accumulate K, and two dense DOF-scale products cover the
    # rigid/energy diagnostics. This intentionally errs high for a scale.
    operation_count = 2 * 6 + n_query + 2 * n_dof
    product = operation_count * unit_roundoff
    gamma = product / (1.0 - product)
    condition_factor = max(1.0, float(maximum_condition))
    return {
        "binary64_epsilon": epsilon,
        "unit_roundoff": unit_roundoff,
        "conservative_operation_count": operation_count,
        "gamma_n": gamma,
        "symmetry_relative_reference": 4.0 * gamma,
        "rigid_residual_relative_reference_times_observed_condition": 8.0 * gamma * condition_factor,
        "affine_relative_reference_times_observed_condition": 8.0 * gamma * condition_factor,
        "eigenvalue_zero_absolute_reference_multiplier": 4.0 * gamma,
        "interpretation": "roundoff-scale diagnostic references only; no pass/fail gate",
    }


def _json_float(value: float) -> float | None:
    return float(value) if math.isfinite(float(value)) else None


def _status_counts(statuses: tuple[MLSStatus, ...]) -> dict[str, int]:
    return {status.value: sum(s is status for s in statuses) for status in MLSStatus}


def _run_diagnostic() -> dict[str, object]:
    started = time.monotonic()
    started_at_utc = _timestamp()
    case = json.loads(CASE_PATH.read_text(encoding="utf-8"))
    _require_fixed_case(case)
    domain = case["domain"]
    quadrature = case["volume_quadrature"]
    assert isinstance(domain, dict) and isinstance(quadrature, dict)
    n_source = int(domain["coefficient_count"])
    n_query = int(quadrature["point_count"])
    dof_count = 3 * n_source
    if (n_source, n_query, dof_count) != (27, 64, 81):
        raise ValueError("fixed case produced unexpected MLS/operator dimensions")
    budget_parts = _array_budget_estimate(n_source, n_query)
    budget_bytes = int(budget_parts["estimated_explicit_array_bytes"])
    if budget_bytes > MAX_ARRAY_BUDGET:
        raise DiagnosticBudgetError(
            f"explicit numeric array estimate {budget_bytes} exceeds {MAX_ARRAY_BUDGET} bytes"
        )
    source_positions, query_positions, volume_weights = _build_fixture(case)
    if source_positions.shape != (n_source, 3) or query_positions.shape != (n_query, 3):
        raise ValueError("fixed fixture builder returned unexpected coordinate dimensions")

    mls_case = case["mls"]
    material_case = case["material"]
    assert isinstance(mls_case, dict) and isinstance(material_case, dict)
    mls_config = MLS3DConfig(
        reference_origin=tuple(float(v) for v in mls_case["reference_origin_m"]),
        reference_scale=float(mls_case["reference_scale_m"]),
        support_radius=float(mls_case["support_radius_m"]),
        max_condition_number=float(mls_case["max_condition_number"]),
        rank_relative_tolerance=float(mls_case["rank_relative_tolerance"]),
        switch_tolerance=float(mls_case["switch_tolerance"]),
        duplicate_tolerance=float(mls_case["duplicate_tolerance"]),
    )
    material = IsotropicMaterial(
        young_modulus=float(material_case["young_modulus_Pa"]),
        poisson_ratio=float(material_case["poisson_ratio"]),
    )

    provider = MLSShapeProvider(source_positions, query_positions, mls_config)
    mls_result = provider.result
    status_counts = _status_counts(mls_result.statuses)
    finite_conditions = mls_result.condition_numbers[np.isfinite(mls_result.condition_numbers)]
    observed_max_condition = float(np.max(finite_conditions)) if finite_conditions.size else 1.0
    roundoff = _roundoff_references(n_query, dof_count, observed_max_condition)
    if time.monotonic() - started > WALL_BUDGET_SECONDS:
        raise DiagnosticBudgetError("30-second budget exceeded after the MLS batch")

    # The affine field is a diagnostic input, not a load or a calibrated material test.
    affine_gradient = np.array(
        [[1.0e-4, -2.0e-5, 3.0e-5], [4.0e-5, -5.0e-5, 2.0e-5], [-1.0e-5, 3.0e-5, 7.0e-5]],
        dtype=np.float64,
    )
    affine_offset = np.array([0.01, -0.02, 0.005], dtype=np.float64)
    coefficient_displacements = source_positions @ affine_gradient.T + affine_offset
    exact_displacements = query_positions @ affine_gradient.T + affine_offset
    interpolation_diagnostic: dict[str, object] = {
        "status": "NOT_EVALUATED_UNKNOWN_MLS_ROW",
        "reference": "binary64 roundoff times observed maximum MLS condition number; diagnostic only",
    }
    if all(status is MLSStatus.OK for status in mls_result.statuses):
        interpolated_displacements = mls_result.weights @ coefficient_displacements
        interpolated_gradients = np.einsum(
            "qnj,nc->qcj", mls_result.gradients, coefficient_displacements
        )
        displacement_abs = np.abs(interpolated_displacements - exact_displacements)
        gradient_abs = np.abs(interpolated_gradients - affine_gradient[None, :, :])
        displacement_scale = max(float(np.max(np.abs(exact_displacements))), np.finfo(float).tiny)
        gradient_scale = max(float(np.max(np.abs(affine_gradient))), np.finfo(float).tiny)
        interpolation_diagnostic = {
            "status": "DIAGNOSTIC_ONLY",
            "displacement_max_absolute_error_m": float(np.max(displacement_abs)),
            "displacement_max_relative_error": float(np.max(displacement_abs)) / displacement_scale,
            "gradient_max_absolute_error_per_m": float(np.max(gradient_abs)),
            "gradient_max_relative_error": float(np.max(gradient_abs)) / gradient_scale,
            "reference_relative_scale": roundoff["affine_relative_reference_times_observed_condition"],
            "reference_is_acceptance_gate": False,
        }

    # Exactly one call into the frozen operator. Any MLS UNKNOWN is returned as
    # an explicit blocked result; there is no retry or matrix repair.
    assembly = assemble_elastic_operator(
        source_positions,
        query_positions,
        volume_weights,
        material,
        provider,
    )
    if provider.consumed_queries != n_query and assembly.status is AssemblyStatus.OK:
        raise RuntimeError("successful assembly did not consume the complete MLS query batch")
    if time.monotonic() - started > WALL_BUDGET_SECONDS:
        raise DiagnosticBudgetError("30-second budget exceeded after the operator assembly")

    result: dict[str, object] = {
        "schema": "efgm_mls_patch_operator_diagnostic.v1",
        "case_id": CASE_ID,
        "status": assembly.status.value,
        "started_at_utc": started_at_utc,
        "fixed_case": case,
        "dimensions": {
            "source_count": n_source,
            "quadrature_count": n_query,
            "displacement_dof_count": assembly.dof_count,
            "volume_sum_m3": float(np.sum(volume_weights)),
        },
        "identity_visibility": {
            "mode": "all candidate sources; no crack surface or visibility mask",
            "mask_locally_constant_asserted": False,
            "crack_surface_loaded": False,
        },
        "mls": {
            "api": "shape_functions_3d",
            "row_order_to_global_ids": "source row i maps to coefficient ID i, dense IDs 0..N-1",
            "shape_status_counts": status_counts,
            "support_count_min": int(np.min(mls_result.support_counts)),
            "support_count_max": int(np.max(mls_result.support_counts)),
            "maximum_finite_condition_number": observed_max_condition,
            "unknown_rows": [
                {
                    "query_index": i,
                    "status": status.value,
                    "reason": mls_result.reasons[i],
                }
                for i, status in enumerate(mls_result.statuses)
                if status is not MLSStatus.OK
            ],
            "condition_numbers": [_json_float(float(v)) for v in mls_result.condition_numbers],
        },
        "array_budget": {
            "explicit_array_estimate_bytes": budget_parts,
            "estimated_explicit_array_bytes": budget_bytes,
            "cap_bytes": MAX_ARRAY_BUDGET,
            "excludes_python_runtime_and_vendor_blas_workspace": True,
            "within_estimate_cap": budget_bytes <= MAX_ARRAY_BUDGET,
        },
        "assembly": {
            "operator_calls": 1,
            "status": assembly.status.value,
            "blocked_quadrature_index": assembly.blocked_quadrature_index,
            "shape_status": assembly.shape_status.value if assembly.shape_status else None,
            "diagnostic": assembly.diagnostic,
        },
        "affine_interpolation": interpolation_diagnostic,
        "roundoff_references": roundoff,
        "limitations": [
            "pure operator diagnostic; no solver, prescribed displacement, traction, or reaction",
            "no crack geometry, split support, visibility query, crack loading, SIF, or propagation",
            "support radius is a fixed numerical MLS support input, not a material parameter",
            "material values are an engineering fixture and are not ceramic calibration",
            "six rigid-mode residual checks establish presence only",
            "spectral extra-null count is diagnostic and does not prove EFGM stability",
        ],
    }

    if assembly.status is AssemblyStatus.OK and assembly.stiffness is not None:
        stiffness = assembly.stiffness
        finite_matrix = bool(np.all(np.isfinite(stiffness)))
        if not finite_matrix:
            raise ElasticOperatorNumericsError("operator returned non-finite K")
        symmetry_relative = float(
            np.linalg.norm(stiffness - stiffness.T, ord="fro")
            / max(np.linalg.norm(stiffness, ord="fro"), np.finfo(float).tiny)
        )
        stiffness_frobenius_norm = float(np.linalg.norm(stiffness, ord="fro"))
        rigid_modes = rigid_mode_coefficient_vectors(source_positions)
        rigid_action = stiffness @ rigid_modes
        rigid_residual_by_mode = [
            float(np.linalg.norm(rigid_action[:, mode_id]))
            / max(
                stiffness_frobenius_norm * float(np.linalg.norm(rigid_modes[:, mode_id])),
                np.finfo(float).tiny,
            )
            for mode_id in range(6)
        ]
        rigid_residual = float(
            np.linalg.norm(rigid_action, ord="fro")
            / max(
                stiffness_frobenius_norm * np.linalg.norm(rigid_modes, ord="fro"),
                np.finfo(float).tiny,
            )
        )
        # Independent tensor energy uses Lamé scalars from E,nu, not operator D.
        e_value = float(material_case["young_modulus_Pa"])
        nu_value = float(material_case["poisson_ratio"])
        lame_lambda = e_value * nu_value / ((1.0 + nu_value) * (1.0 - 2.0 * nu_value))
        mu = e_value / (2.0 * (1.0 + nu_value))
        strain = 0.5 * (affine_gradient + affine_gradient.T)
        stress = 2.0 * mu * strain + lame_lambda * float(np.trace(strain)) * np.eye(3)
        independent_tensor_energy = 0.5 * float(case["volume_quadrature"]["total_volume_m3"]) * float(
            np.sum(strain * stress)
        )
        operator_affine_energy = 0.5 * float(
            coefficient_displacements.reshape(-1)
            @ stiffness
            @ coefficient_displacements.reshape(-1)
        )
        energy_relative_difference = abs(operator_affine_energy - independent_tensor_energy) / max(
            abs(independent_tensor_energy), np.finfo(float).tiny
        )

        symmetric_part = 0.5 * (stiffness + stiffness.T)
        eigenvalues = np.linalg.eigvalsh(symmetric_part)
        spectral_scale = max(float(np.max(np.abs(eigenvalues))), np.finfo(float).tiny)
        eigenvalue_zero_tolerance = float(roundoff["eigenvalue_zero_absolute_reference_multiplier"]) * spectral_scale
        near_zero_count = int(np.count_nonzero(np.abs(eigenvalues) <= eigenvalue_zero_tolerance))
        negative_significant = eigenvalues[eigenvalues < -eigenvalue_zero_tolerance]
        raw_negative = eigenvalues[eigenvalues < 0.0]
        result["operator_diagnostics"] = {
            "K_shape": list(stiffness.shape),
            "K_units": "N/m",
            "K_finite": finite_matrix,
            "K_frobenius_norm_N_per_m": stiffness_frobenius_norm,
            "symmetry": {
                "relative_frobenius_asymmetry": symmetry_relative,
                "roundoff_reference": roundoff["symmetry_relative_reference"],
                "reference_is_acceptance_gate": False,
            },
            "rigid_modes": {
                "expected_modes_checked": 6,
                "coefficient_mode_matrix_rank": int(np.linalg.matrix_rank(rigid_modes)),
                "scaled_residual": rigid_residual,
                "mode_order": ["translation_x", "translation_y", "translation_z", "rotation_x", "rotation_y", "rotation_z"],
                "scaled_residual_by_mode": rigid_residual_by_mode,
                "roundoff_reference_times_observed_condition": roundoff[
                    "rigid_residual_relative_reference_times_observed_condition"
                ],
                "interpretation": "presence of six physical rigid modes only; extra null modes remain possible",
                "reference_is_acceptance_gate": False,
            },
            "affine_patch_energy": {
                "operator_energy_J": operator_affine_energy,
                "independent_tensor_energy_J": independent_tensor_energy,
                "relative_difference": energy_relative_difference,
                "roundoff_reference_times_observed_condition": roundoff[
                    "affine_relative_reference_times_observed_condition"
                ],
                "reference_is_acceptance_gate": False,
            },
            "spectrum": {
                "matrix_used": "symmetric_part=(K+K.T)/2 for eigvalsh only; raw K is unchanged",
                "eigenvalues_N_per_m": [float(value) for value in eigenvalues],
                "raw_negative_eigenvalue_count": int(raw_negative.size),
                "negative_beyond_roundoff_count": int(negative_significant.size),
                "near_zero_count": near_zero_count,
                "additional_near_zero_after_six_rigid_count": max(near_zero_count - 6, 0),
                "absolute_near_zero_threshold_N_per_m": eigenvalue_zero_tolerance,
                "interpretation": "diagnostic spectrum only; no positive-definiteness enforcement",
                "reference_is_acceptance_gate": False,
            },
        }
    else:
        result["operator_diagnostics"] = None

    elapsed = time.monotonic() - started
    result["timing"] = {
        "elapsed_seconds": elapsed,
        "wall_budget_seconds": WALL_BUDGET_SECONDS,
        "within_wall_budget": elapsed <= WALL_BUDGET_SECONDS,
    }
    dependencies = {
        "mls3d_shape.py": MLS_DIR / "mls3d_shape.py",
        "mls3d_mls.py": MLS_DIR / "mls3d_mls.py",
        "efgm_elastic_operator.py": OPERATOR_DIR / "efgm_elastic_operator.py",
        "mls_elastic_adapter.py": HERE / "mls_elastic_adapter.py",
        "operator_diagnostic.py": HERE / "operator_diagnostic.py",
        "fixed_case.json": CASE_PATH,
    }
    result["sha256"] = {name: _sha256(path) for name, path in dependencies.items()}
    result["finished_at_utc"] = _timestamp()
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=HERE / "results" / "operator_diagnostic.json",
        help="JSON output path; does not alter the frozen case",
    )
    parser.add_argument(
        "--_worker",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    args = parser.parse_args()
    output_path = args.output
    if args._worker:
        _load_numerical_dependencies()
        try:
            result = _run_diagnostic()
            exit_code = 0 if result["status"] == AssemblyStatus.OK.value else 2
        except Exception as exc:
            result = {
                "schema": "efgm_mls_patch_operator_diagnostic.v1",
                "case_id": CASE_ID,
                "status": "DIAGNOSTIC_ERROR",
                "error_type": type(exc).__name__,
                "error": str(exc),
                "traceback": traceback.format_exc(limit=8),
                "no_retry": True,
                "finished_at_utc": _timestamp(),
            }
            exit_code = 2
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(
            json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        print(json.dumps(result, ensure_ascii=False, allow_nan=False))
        return exit_code

    # Run the one MLS/assembly diagnostic in a child process so the stated
    # wall budget can actually terminate a stuck native numerical call. The
    # worker writes only to a unique staging file; the parent publishes one
    # final JSON record after the worker exits.
    output_path.parent.mkdir(parents=True, exist_ok=True)
    staging_path = output_path.with_name(
        f".{output_path.name}.{os.getpid()}.worker.json"
    )
    staging_path.unlink(missing_ok=True)
    parent_started = time.monotonic()
    parent_started_at = _timestamp()
    try:
        completed = subprocess.run(
            [
                sys.executable,
                str(HERE / "operator_diagnostic.py"),
                "--_worker",
                "--output",
                str(staging_path),
            ],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=WALL_BUDGET_SECONDS,
            check=False,
        )
        try:
            result = json.loads(staging_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            result = {
                "schema": "efgm_mls_patch_operator_diagnostic.v1",
                "case_id": CASE_ID,
                "status": "WORKER_FAILED",
                "worker_return_code": completed.returncode,
                "worker_stdout_tail": completed.stdout[-2000:],
                "worker_stderr_tail": completed.stderr[-2000:],
                "no_retry": True,
                "finished_at_utc": _timestamp(),
            }
            exit_code = 2
        else:
            exit_code = completed.returncode
    except subprocess.TimeoutExpired:
        result = {
            "schema": "efgm_mls_patch_operator_diagnostic.v1",
            "case_id": CASE_ID,
            "status": "WALL_TIMEOUT",
            "wall_budget_seconds": WALL_BUDGET_SECONDS,
            "elapsed_seconds": time.monotonic() - parent_started,
            "no_retry": True,
            "no_result_published": True,
            "finished_at_utc": _timestamp(),
        }
        exit_code = 2
    finally:
        staging_path.unlink(missing_ok=True)

    parent_elapsed = time.monotonic() - parent_started
    result["parent_started_at_utc"] = parent_started_at
    result["parent_elapsed_seconds"] = parent_elapsed
    result["parent_wall_budget_seconds"] = WALL_BUDGET_SECONDS
    result["parent_within_wall_budget"] = parent_elapsed <= WALL_BUDGET_SECONDS
    output_path.write_text(
        json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(result, ensure_ascii=False, allow_nan=False))
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
