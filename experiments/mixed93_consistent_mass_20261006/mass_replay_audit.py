"""Saved-array audit, blocked until the case has root pins and authorization."""

import hashlib
import json
import math
from pathlib import Path
import time
from mixed93_consistent_mass import assemble_consistent_mass
from recovery_metadata_gate import verify_recovery_contract

CASE, ROOT = Path(__file__).with_name("case.json"), Path(__file__).resolve().parents[3]


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def scaled_error(left, right, np):
    scale = np.sqrt(np.maximum(np.diag(right), np.finfo(float).tiny))
    return float(np.max(np.abs((left - right) / scale[:, None] / scale[None, :])))


def main():
    started = time.perf_counter(); case_bytes = CASE.read_bytes()
    case_pre = hashlib.sha256(case_bytes).hexdigest(); case = json.loads(case_bytes)
    pins = case.get("inputs", {})
    if case.get("execution", {}).get("authorized") is not True:
        raise RuntimeError("BLOCKED: execution.authorized must be explicitly true")
    keys = ("result_path", "result_sha256", "arrays_path", "arrays_sha256")
    if any(not pins.get(key) for key in keys):
        raise RuntimeError("BLOCKED: failed arrays path and SHA-256 pins are required")
    resolve = lambda p: Path(p) if Path(p).is_absolute() else ROOT / p
    arrays_path = resolve(pins["arrays_path"]).resolve()
    arrays_sha = sha256(arrays_path)
    if arrays_sha.lower() != pins["arrays_sha256"].lower():
        raise RuntimeError("BLOCKED: arrays SHA-256 differs from the fixed case pin")
    recovery_metadata = verify_recovery_contract(case, ROOT, arrays_path, arrays_sha)
    output = resolve(case["planned_output_dir"])
    if output.exists():
        raise RuntimeError(f"BLOCKED: refusing to overwrite {output}")
    core = Path(__file__).with_name("mixed93_consistent_mass.py")
    recovery_source = Path(__file__).with_name("recovery_metadata_gate.py")
    core_pre, audit_pre, recovery_source_pre = sha256(core), sha256(Path(__file__)), sha256(recovery_source)
    import numpy as np
    with np.load(arrays_path, allow_pickle=False) as data:
        required = {"semantic_ids", "source_ids", "source_positions_xyz_m", "volume_basis_values", "volume_weights_m3", "volume_xyz_m", "rigid_coefficient_modes"}
        missing = sorted(required.difference(data.files))
        if missing:
            raise ValueError(f"cached arrays missing: {missing}")
        ids = tuple(str(v) for v in data["semantic_ids"])
        sources = tuple(str(v) for v in data["source_ids"])
        B, w = data["volume_basis_values"], data["volume_weights_m3"]
        xyz, positions = np.asarray(data["volume_xyz_m"]), np.asarray(data["source_positions_xyz_m"])
        C = np.asarray(data["rigid_coefficient_modes"])
        if len(sources) != 27 or len(set(sources)) != 27 or C.shape != (6, 93): raise ValueError("frozen source IDs or rigid mode shape mismatch")
        if xyz.shape != (512, 3) or positions.shape != (27, 3): raise ValueError("cached position array shape mismatch")
        bases, components, owners = [], [], []
        for slot, raw_id in enumerate(ids):
            scope, owner, component, branch = json.loads(raw_id)
            if scope == "SOURCE_BASE":
                bases.append(sources.index(owner))
                if slot < 81 and slot % 3 == 0:
                    owners.append(owner)
            else:
                bases.append(27 + branch - 1)
            if scope not in ("SOURCE_BASE", "GLOBAL_FRONT"): raise ValueError(f"unsupported semantic scope: {scope}")
            components.append(("ux", "uy", "uz").index(component))
        if tuple(owners) != sources:
            raise ValueError("semantic source order differs from saved source_ids")
        raw = assemble_consistent_mass(ids, B, w, case["mass"]["density_kg_m3"])
        M = np.frombuffer(raw.values_le_f64, dtype="<f8").reshape(93, 93)
        center = np.asarray(case["geometry"]["center_xyz_m"], dtype=np.float64)
        if not np.isfinite(xyz).all() or not np.isfinite(positions).all() or not np.isfinite(C).all(): raise ValueError("cached positions and rigid modes must be finite")
        expected_modes = np.zeros_like(C); source_center = np.mean(positions, axis=0)
        for i, position in enumerate(positions):
            expected_modes[:3, 3*i:3*i+3] = np.eye(3); expected_modes[3:, 3*i:3*i+3] = np.cross(np.eye(3), position - source_center)
        rigid_coeff_error = float(np.max(np.abs(C - expected_modes)))
        diagonal = np.diag(M)
        if np.any(diagonal <= 0.0): raise ValueError("raw mass has nonpositive diagonal")
        scale = 1.0 / np.sqrt(diagonal)
        scaled = scale[:, None] * M * scale[None, :]
        asym = float(np.linalg.norm(M - M.T, ord=np.inf) / np.linalg.norm(M, ord=np.inf))
        eigen = np.linalg.eigvalsh(scaled)
        eig_tol = float(64.0 * np.finfo(float).eps * 93 * max(float(np.max(np.abs(eigen))), np.finfo(float).tiny))
        npos, nzero, nneg = (int(np.count_nonzero(test)) for test in (eigen > eig_tol, np.abs(eigen) <= eig_tol, eigen < -eig_tol))
        cholesky = False
        if asym <= 1.0e-12 and (npos, nzero, nneg) == (93, 0, 0):
            try:
                np.linalg.cholesky(scaled)
                cholesky = True
            except np.linalg.LinAlgError:
                pass
        blocks = [[i for i, c in enumerate(components) if c == k] for k in range(3)]
        block_delta = max(float(np.max(np.abs(M[np.ix_(blocks[k], blocks[k])] - M[np.ix_(blocks[0], blocks[0])]))) for k in (1, 2))
        cross = max(float(np.max(np.abs(M[np.ix_(blocks[i], blocks[j])]))) for i in range(3) for j in range(3) if i != j)
        fitted = np.zeros((512, 3, 6))
        for slot, basis in enumerate(bases):
            fitted[:, components[slot], :] += B[:, basis, None] * C[:, slot][None, :]
        physical = np.zeros_like(fitted)
        relative = xyz - center
        for axis in range(3):
            physical[:, axis, axis] = 1.0
            physical[:, :, 3 + axis] = np.cross(np.eye(3)[axis], relative)
        mass_weights = raw.density_kg_m3 * w; km = C @ M @ C.T
        kb = np.einsum("q,qcm,qcn->mn", mass_weights, fitted, fitted); kp = np.einsum("q,qcm,qcn->mn", mass_weights, physical, physical)
        radius, length = case["geometry"]["radius_m"], case["geometry"]["length_m"]; analytic_volume = math.pi * radius**2 * length
        mass_value = raw.density_kg_m3 * analytic_volume; ix = mass_value * (radius**2 / 4.0 + length**2 / 12.0); iz = mass_value * radius**2 / 2.0
        analytic = np.diag([mass_value] * 3 + [ix, ix, iz])
        e_modes, e_quadrature, e_analytic = scaled_error(km, kb, np), scaled_error(kb, kp, np), scaled_error(kp, analytic, np)
        volume_error = abs(float(np.sum(w)) - analytic_volume) / analytic_volume
        rotation_scales = [math.sqrt(radius**2 + (length / 2.0)**2)] * 2 + [radius]
        mode_scales = np.asarray([1.0, 1.0, 1.0, *rotation_scales], dtype=np.float64)
        pointwise_by_mode = np.max(np.abs(fitted - physical), axis=(0, 1)) / mode_scales
        pointwise_error = float(np.max(pointwise_by_mode))
        audit_arrays = (M, scaled, eigen, km, kb, kp, analytic, fitted, physical, pointwise_by_mode)
        if any(not np.isfinite(array).all() for array in audit_arrays): raise ValueError("mass or rigid audit produced nonfinite values")
        passed = bool((npos, nzero, nneg) == (93, 0, 0) and cholesky and asym <= 1.0e-12 and cross == 0.0
                      and block_delta == 0.0 and max(e_modes, e_quadrature, e_analytic) <= 1.0e-10 and volume_error <= 1.0e-10 and rigid_coeff_error <= 1.0e-12)
        passed = bool(passed and pointwise_error <= case["future_audit"]["rigid_pointwise_tolerance"])
        metrics = {
            "schema": "pgfrag.mixed93_consistent_mass_audit.v2", "case_id": case["case_id"],
            "status": "MASS_DIAGNOSTIC_GATES_PASSED_PRODUCER_FAILURE_RETAINED" if passed else "BLOCKED_MASS_AUDIT_GATE_PRODUCER_FAILURE_RETAINED",
            "mass_audit_gate_status": "PASS" if passed else "BLOCKED", "matrix_shape": [93, 93], "matrix_units": "kg",
            "producer_provenance": recovery_metadata,
            "density_kg_m3": raw.density_kg_m3, "density_status": "NOT_MATERIAL_CALIBRATED",
            "volume_m3": raw.volume_m3, "trace_kg": float(np.trace(M)), "relative_raw_asymmetry": asym,
            "cross_component_max_abs_kg": cross, "component_block_difference_kg": block_delta,
            "scaled_spectrum": {"positive": npos, "zero": nzero, "significant_negative": nneg,
                                "tolerance": eig_tol, "minimum": float(eigen[0]), "maximum": float(eigen[-1]),
                                "cholesky_succeeded": cholesky},
            "rigid_kinetic_block_units": {"translation_translation": "kg", "translation_rotation": "kg m", "rotation_rotation": "kg m²"},
            "rigid_kinetic_raw_M": km.tolist(), "rigid_kinetic_direct_quadrature": kp.tolist(),
            "rigid_kinetic_centered_cylinder_analytic": analytic.tolist(),
            "rigid_pointwise_error_max_normalized": pointwise_error,
            "rigid_pointwise_error_by_mode_normalized": pointwise_by_mode.tolist(),
            "rigid_pointwise_mode_scales": {"translation_dimensionless": [1.0, 1.0, 1.0], "rotation_scales_m": rotation_scales},
            "rigid_kinetic_errors": {"coefficient_projection": e_modes, "direct_quadrature": e_quadrature, "analytic_inertia": e_analytic,
                                     "volume": volume_error, "coefficient_source_reconstruction_max_abs": rigid_coeff_error},
        }
        raw_bytes = raw.values_le_f64
    recovery_path = Path(recovery_metadata["recovery_json_path"])
    if sha256(arrays_path).lower() != arrays_sha.lower() or sha256(recovery_path).lower() != recovery_metadata["recovery_json_sha256"].lower():
        raise RuntimeError("BLOCKED: a pinned input changed during mass audit")
    if Path(recovery_metadata["producer_result_json_path"]).exists():
        raise RuntimeError("BLOCKED: producer result JSON appeared during mass audit")
    case_post, core_post, audit_post, recovery_source_post = sha256(CASE), sha256(core), sha256(Path(__file__)), sha256(recovery_source)
    if (case_post, core_post, audit_post, recovery_source_post) != (case_pre, core_pre, audit_pre, recovery_source_pre):
        raise RuntimeError("BLOCKED: case or source changed during mass audit")
    metrics["elapsed_seconds_inside_checker"] = time.perf_counter() - started
    metrics["input_pins"] = {"producer_result_path": recovery_metadata["producer_result_json_path"], "producer_result_sha256": None,
                             "case_result_alias_path": recovery_metadata["case_result_alias_path"],
                             "case_result_alias_sha256": recovery_metadata["case_result_alias_sha256"],
                             "arrays_path": str(arrays_path), "arrays_sha256": arrays_sha,
                             "recovery_json_path": recovery_metadata["recovery_json_path"],
                             "recovery_json_sha256": recovery_metadata["recovery_json_sha256"],
                             "mass_source_sha256_before": core_pre, "mass_source_sha256_after": core_post,
                             "recovery_gate_source_sha256_before": recovery_source_pre,
                             "recovery_gate_source_sha256_after": recovery_source_post,
                             "audit_source_sha256_before": audit_pre, "audit_source_sha256_after": audit_post,
                             "case_sha256_before": case_pre, "case_sha256_after": case_post}
    output.mkdir(parents=True, exist_ok=False)
    (output / "mass_raw_le_f64.bin").write_bytes(raw_bytes)
    np.savetxt(output / "mass_raw_kg.csv", np.frombuffer(raw_bytes, dtype="<f8").reshape(93, 93), delimiter=",", fmt="%.17g")
    (output / "result.json").write_text(json.dumps(metrics, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    return 0 if passed else 1
if __name__ == "__main__":
    raise SystemExit(main())
