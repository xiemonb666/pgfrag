from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(r"C:\Users\xiemo\Desktop\ShardForge")
REVIEW = ROOT / "pgfrag/tmp/review_20261007/mixed93_moving_front_independent"
EXPERIMENT = ROOT / "pgfrag/experiments/mixed93_moving_front_reanalysis_20261007"
RUN = ROOT / "pgfrag/tmp/review_20261007/mixed93_moving_front_reanalysis_20261007/run01"
CAPTURE = ROOT / "pgfrag/tmp/review_20261007/mixed93_moving_front_reanalysis_20261007/capture_once01/capture_receipt.json"
ADDENDUM = ROOT / "pgfrag/tmp/review_20261007/mixed93_moving_front_reanalysis_20261007/root_go_outer_addendum.json"
PINS = ROOT / "pgfrag/tmp/review_20261007/mixed93_moving_front_reanalysis_20261007/SOURCE_PINS.json"
LIMIT_SECONDS = 10.0
LIMIT_OUTPUT_BYTES = 256 * 1024


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def pinned_paths(capture: dict) -> dict[str, Path]:
    paths = {key: ROOT / Path(key.replace("/", os.sep).replace("\\", os.sep))
             for key in capture.get("inputs_pre", {})}
    case = json.loads((EXPERIMENT / "case.json").read_text(encoding="utf-8"))
    paths["closeout.entry"] = EXPERIMENT / "mixed93_moving_front_reanalysis.py"
    paths["closeout.case"] = EXPERIMENT / "case.json"
    paths["closeout.source_pins"] = PINS
    paths["closeout.failure_json"] = RUN / "failure.json"
    paths["closeout.failure_arrays_npz"] = RUN / "failure_arrays.npz"
    paths["closeout.capture_receipt"] = CAPTURE
    paths["closeout.outer_addendum"] = ADDENDUM
    paths["closeout.audit"] = ROOT / Path(case["independent_audit_path"].replace("/", os.sep))
    return paths


def hash_set(paths: dict[str, Path]) -> dict[str, dict]:
    return {key: {"path": str(path.relative_to(ROOT)).replace("\\", "/"),
                  "sha256": sha256(path), "bytes": path.stat().st_size}
            for key, path in sorted(paths.items())}


def child_readback() -> int:
    try:
        import numpy as np

        failure = json.loads((RUN / "failure.json").read_text(encoding="utf-8"))
        archive_path = RUN / failure["replay_arrays_npz"]["path"]
        with np.load(archive_path, allow_pickle=False) as archive:
            keys = list(archive.files)
            arrays = {key: archive[key] for key in keys}

        array_headers = {key: {"shape": list(value.shape), "dtype": str(value.dtype),
                               "nbytes": int(value.nbytes)}
                         for key, value in arrays.items()}
        finite_counts = {}
        for key, value in arrays.items():
            if value.dtype.kind in "fc":
                finite_counts[key] = {"finite": int(np.isfinite(value).sum()),
                                      "total": int(value.size)}

        batches = []
        for index, record in enumerate(failure.get("provider_batches", [])):
            prefix = f"cap_batch_{index}_"
            batch = {"start": record.get("start"), "stop": record.get("stop"),
                     "requested": record.get("requested"),
                     "returned_statuses": record.get("returned_statuses"), "arrays": {}}
            for suffix in ("xyz_m", "weights", "gradients", "support_counts", "condition_numbers"):
                name = prefix + suffix
                if name in arrays:
                    value = arrays[name]
                    batch["arrays"][suffix] = array_headers[name]
            batches.append(batch)
        expected_sizes = [256, 192, 256, 256, 256, 128]
        batch_sizes = [int(item.get("requested", -1)) for item in failure.get("provider_batches", [])]
        provider_ok = (batch_sizes == expected_sizes
                       and failure.get("provider_rows_returned") == 1344
                       and len(failure.get("provider_statuses", [])) == 1344
                       and all(value == "OK" for value in failure.get("provider_statuses", [])))
        for index, count in enumerate(expected_sizes):
            for suffix, shape in (("xyz_m", [count, 3]), ("weights", [count, 27]),
                                  ("gradients", [count, 27, 3]),
                                  ("support_counts", [count]), ("condition_numbers", [count])):
                name = f"cap_batch_{index}_{suffix}"
                if name not in arrays or list(arrays[name].shape) != shape:
                    provider_ok = False
                elif arrays[name].dtype.kind in "fc" and not bool(np.isfinite(arrays[name]).all()):
                    provider_ok = False

        front_checks = {}
        eps = np.finfo(np.float64).eps
        for front in ("old_front", "candidate_front"):
            names = {suffix: f"{front}_{suffix}" for suffix in (
                "K_raw_N_per_m", "K_scaled_raw", "eigenvalues_scaled",
                "rigid_mode_coefficients", "C_physical", "rigid_gauge_coupling")}
            if any(name not in arrays for name in names.values()):
                front_checks[front] = {"status": "UNKNOWN", "missing": [
                    name for name in names.values() if name not in arrays]}
                continue
            kraw = arrays[names["K_raw_N_per_m"]]
            kstored = arrays[names["K_scaled_raw"]]
            stored_eigs = arrays[names["eigenvalues_scaled"]]
            rigid = arrays[names["rigid_mode_coefficients"]]
            cphys = arrays[names["C_physical"]]
            coupling = arrays[names["rigid_gauge_coupling"]]
            check = {"K_shape": list(kraw.shape), "K_scaled_shape": list(kstored.shape),
                     "eigenvalue_count": int(stored_eigs.size),
                     "rigid_shape": list(rigid.shape), "C_shape": list(cphys.shape),
                     "coupling_shape": list(coupling.shape)}
            shapes_ok = (kraw.shape == (93, 93) and kstored.shape == (93, 93)
                         and stored_eigs.shape == (93,) and rigid.shape == (6, 93)
                         and cphys.shape == (6, 93) and coupling.shape == (6, 6))
            check["finite"] = all(bool(np.isfinite(x).all()) for x in
                                  (kraw, kstored, stored_eigs, rigid, cphys, coupling))
            raw_norm = max(float(np.linalg.norm(kraw)), np.finfo(float).tiny)
            scaled_norm = max(float(np.linalg.norm(kstored)), np.finfo(float).tiny)
            check["raw_relative_asymmetry"] = float(np.linalg.norm(kraw-kraw.T) / raw_norm)
            check["stored_scaled_relative_asymmetry"] = float(np.linalg.norm(kstored-kstored.T) / scaled_norm)
            diag_ok = bool(np.all(np.diag(kraw) > 0.0))
            check["positive_raw_diagonal"] = diag_ok
            if diag_ok:
                scale = 1.0 / np.sqrt(np.diag(kraw))
                reconstructed = scale[:, None] * kraw * scale[None, :]
                check["scaled_matrix_reconstruction_relative_error"] = float(
                    np.linalg.norm(reconstructed-kstored) / scaled_norm)
            else:
                check["scaled_matrix_reconstruction_relative_error"] = None
            scale_vector = 1.0 / np.sqrt(np.diag(kraw)) if diag_ok else None
            if scale_vector is not None:
                modes = rigid.T / scale_vector[:, None]
                check["rigid_null_relative_residuals"] = [float(
                    np.linalg.norm(kstored @ modes[:, j], ord=np.inf) /
                    max(float(np.linalg.norm(kstored, ord=np.inf) *
                              np.linalg.norm(modes[:, j], ord=np.inf)), np.finfo(float).tiny))
                    for j in range(6)]
            else:
                check["rigid_null_relative_residuals"] = None
            expected_coupling = cphys @ rigid.T
            check["gauge_coupling_max_abs_error"] = float(np.max(np.abs(expected_coupling-coupling)))
            try:
                sv = np.linalg.svd(coupling, compute_uv=False)
                check["gauge_coupling_singular_values"] = [float(v) for v in sv]
                check["gauge_coupling_rank_default"] = int(np.linalg.matrix_rank(coupling))
            except Exception as exc:
                check["gauge_rank_status"] = f"UNKNOWN: {type(exc).__name__}: {exc}"
            spectral_scale = float(np.max(np.abs(stored_eigs))) if stored_eigs.size else 0.0
            tolerance = 64.0 * eps * 93 * spectral_scale
            check["saved_eigenvalue_counts"] = {
                "tolerance": tolerance,
                "zero": int(np.count_nonzero(np.abs(stored_eigs) <= tolerance)),
                "positive": int(np.count_nonzero(stored_eigs > tolerance)),
                "significant_negative": int(np.count_nonzero(stored_eigs < -tolerance)),
            }
            try:
                independent_eigs = np.linalg.eigvalsh(kstored)
                check["eigvalsh_recheck_max_abs_delta"] = float(
                    np.max(np.abs(independent_eigs-stored_eigs)))
                check["eigvalsh_recheck_finite"] = bool(np.isfinite(independent_eigs).all())
            except Exception as exc:
                check["eigvalsh_recheck_status"] = f"UNKNOWN: {type(exc).__name__}: {exc}"
            check["layout_and_finite"] = bool(shapes_ok and check["finite"])
            front_checks[front] = check

        uncompressed = sum(int(value.nbytes) for value in arrays.values())
        output = {
            "schema": "mixed93_saved_failure_readback.v1",
            "status": "SAVED_ARRAY_DIAGNOSTIC_ONLY",
            "failure_status": failure.get("status"),
            "last_completed_phase": failure.get("last_completed_phase"),
            "failure_error": failure.get("error"),
            "archive_header_keys": keys,
            "archive_array_headers": array_headers,
            "archive_uncompressed_bytes": uncompressed,
            "failure_json_uncompressed_bytes_match":
                uncompressed == failure.get("checkpoint_arrays_uncompressed_bytes"),
            "floating_array_finite_counts": finite_counts,
            "provider": {"batch_sizes": batch_sizes, "returned_rows": failure.get("provider_rows_returned"),
                         "status_count": len(failure.get("provider_statuses", [])),
                         "all_statuses_ok": provider_ok, "batches": batches},
            "front_matrix_checks": front_checks,
            "nonperformed": ["provider", "field", "assembly", "KKT solve", "new queries", "front update"],
            "physical_or_growth_qualification": "UNKNOWN_NOT_EVALUATED",
        }
        print(json.dumps(output, indent=2, sort_keys=True, allow_nan=False))
        return 0
    except Exception as exc:
        print(json.dumps({"schema": "mixed93_saved_failure_readback.v1",
                          "status": "UNKNOWN_READBACK_ERROR",
                          "error_type": type(exc).__name__, "error": str(exc),
                          "nonperformed": ["provider", "field", "assembly", "KKT solve", "new queries"]},
                         indent=2, sort_keys=True, allow_nan=False))
        return 0


def main() -> int:
    if sys.argv[1:] == ["--child"]:
        return child_readback()
    out_dir = REVIEW / "saved_failure_readback01"
    out_dir.mkdir(parents=True, exist_ok=False)
    capture = json.loads(CAPTURE.read_text(encoding="utf-8"))
    paths = pinned_paths(capture)
    before = hash_set(paths)
    child_argv = [sys.executable, "-B", str(Path(__file__).resolve()), "--child"]
    started = time.perf_counter()
    process = subprocess.Popen(child_argv, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    timed_out = False
    try:
        stdout, stderr = process.communicate(timeout=LIMIT_SECONDS)
    except subprocess.TimeoutExpired:
        timed_out = True
        process.kill()
        stdout, stderr = process.communicate()
    elapsed = time.perf_counter() - started
    (out_dir / "stdout.bin").write_bytes(stdout)
    (out_dir / "stderr.bin").write_bytes(stderr)
    after = hash_set(paths)
    capture_inputs = capture.get("inputs_post", {})
    pin_mismatches = []
    for key, record in capture_inputs.items():
        current = after.get(key, {}).get("sha256")
        if current != record.get("sha256"):
            pin_mismatches.append(key)
    receipt = {
        "schema": "mixed93_saved_failure_readback_capture.v1",
        "status": "CAPTURED_READBACK_DIAGNOSTIC",
        "parent_argv": [sys.executable, "-B", str(Path(__file__).resolve())],
        "child_argv": child_argv,
        "child_return_code": process.returncode,
        "child_timed_out": timed_out,
        "child_timeout_seconds": LIMIT_SECONDS,
        "child_elapsed_seconds": elapsed,
        "stdout": {"path": "stdout.bin", "bytes": len(stdout),
                   "sha256": hashlib.sha256(stdout).hexdigest()},
        "stderr": {"path": "stderr.bin", "bytes": len(stderr),
                   "sha256": hashlib.sha256(stderr).hexdigest()},
        "output_limit_bytes": LIMIT_OUTPUT_BYTES,
        "input_pins_pre": before,
        "input_pins_post": after,
        "capture_receipt_input_pin_mismatches": pin_mismatches,
        "pins_stable_during_readback": before == after,
        "failure_arrays_npz_read_by_one_child": True,
        "provider_or_solver_calls": 0,
        "receipt_included_in_output_tally": True,
    }
    receipt_path = out_dir / "receipt.json"
    for _ in range(16):
        encoded = json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False).encode("utf-8") + b"\n"
        total = len(stdout) + len(stderr) + len(encoded)
        if receipt.get("combined_saved_output_bytes") == total:
            break
        receipt["combined_saved_output_bytes"] = total
    encoded = json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False).encode("utf-8") + b"\n"
    total = len(stdout) + len(stderr) + len(encoded)
    if receipt.get("combined_saved_output_bytes") != total:
        raise RuntimeError("receipt byte tally did not reach a fixed point")
    receipt_path.write_bytes(encoded)
    print(json.dumps({"status": receipt["status"], "child_rc": process.returncode,
                      "elapsed_seconds": elapsed, "combined_saved_output_bytes": total,
                      "pins_stable": before == after, "pin_mismatches": pin_mismatches},
                     sort_keys=True))
    return 0 if (not timed_out and elapsed <= LIMIT_SECONDS and total <= LIMIT_OUTPUT_BYTES
                 and before == after and not pin_mismatches and process.returncode == 0) else 2


if __name__ == "__main__":
    raise SystemExit(main())
