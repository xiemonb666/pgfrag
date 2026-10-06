"""One-shot stdlib process capture for the authorized local study."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path.cwd().resolve()
CAPTURE = Path(__file__).resolve().parent
REVIEW = CAPTURE.parent
CASE = ROOT / "pgfrag/experiments/tip_aligned_quadrature_study_20261006/study_case.json"
ENTRY = ROOT / "pgfrag/experiments/tip_aligned_quadrature_study_20261006/tip_aligned_quadrature_study.py"
OUTPUT = ROOT / "pgfrag/tmp/review_20261006/tip_aligned_quadrature_study_20261006/run01"
PINS = REVIEW / "SOURCE_PINS.json"
PRE_CASE = REVIEW / "case_pre_authorization.json"
PRE_PINS = REVIEW / "SOURCE_PINS_pre_authorization.json"
AUTHORIZED_AT = "2026-10-06T08:49:55Z"
PYTHON = Path(r"C:\Users\xiemo\Desktop\ShardForge\third_party\vesselforge-runtime\venv\Scripts\python.exe")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def hash_map(paths: dict[str, str]) -> dict[str, str]:
    return {name: sha256(ROOT / name) for name in paths}


def artifacts() -> list[dict[str, object]]:
    if not OUTPUT.exists():
        return []
    return [
        {"path": str(path.relative_to(ROOT)), "sha256": sha256(path), "bytes": path.stat().st_size}
        for path in sorted(OUTPUT.rglob("*")) if path.is_file()
    ]


def main() -> int:
    if OUTPUT.exists():
        raise RuntimeError(f"run01 appeared after preflight; refusing rerun/overwrite: {OUTPUT}")
    if not all(path.is_file() for path in (CASE, ENTRY, PINS, PRE_CASE, PRE_PINS, PYTHON)):
        raise FileNotFoundError("study input, authorization snapshots, or required runtime is missing")
    if PRE_PINS.read_bytes() != PINS.read_bytes():
        raise RuntimeError("pre-authorization 13-pin manifest snapshot differs from current pins")

    manifest = json.loads(PINS.read_text(encoding="utf-8"))
    pre_case = json.loads(PRE_CASE.read_text(encoding="utf-8"))
    case = json.loads(CASE.read_text(encoding="utf-8"))
    if manifest.get("execution_authorized") is not False:
        raise RuntimeError("preserved source pins must retain the pre-authorization false state")
    if pre_case["execution_authorization"]["authorized"] is not False:
        raise RuntimeError("preserved pre-authorization case is not the original false authorization")
    authorization = case["execution_authorization"]
    if authorization.get("authorized") is not True or authorization.get("authorized_at_utc") != AUTHORIZED_AT:
        raise RuntimeError("case authorization does not match the recorded one-time root authorization")
    if case.get("status") != "SOURCE_ONLY_AWAITING_ROOT_REVIEW":
        raise RuntimeError("unexpected source status mutation before execution")
    pre_case["execution_authorization"]["authorized"] = True
    pre_case["execution_authorization"]["authorized_at_utc"] = AUTHORIZED_AT
    if case != pre_case:
        raise RuntimeError("case contains changes beyond the authorized and timestamp fields")

    study_pins = manifest["study_sources"]
    immutable_pins = manifest["immutable_inputs"]
    all_pins = {**study_pins, **immutable_pins}
    if sha256(PRE_CASE) != study_pins["pgfrag/experiments/tip_aligned_quadrature_study_20261006/study_case.json"]:
        raise RuntimeError("preserved pre-authorization case does not match its frozen source pin")
    before_sources = hash_map(study_pins)
    before_inputs = hash_map(immutable_pins)
    for name, expected in immutable_pins.items():
        if before_inputs[name] != expected:
            raise RuntimeError(f"immutable input pin mismatch before run: {name}")
    driver_name = "pgfrag/experiments/tip_aligned_quadrature_study_20261006/tip_aligned_quadrature_study.py"
    readme_name = "pgfrag/experiments/tip_aligned_quadrature_study_20261006/README.md"
    for name in (driver_name, readme_name):
        if before_sources[name] != study_pins[name]:
            raise RuntimeError(f"study source pin mismatch before run: {name}")
    if before_sources["pgfrag/experiments/tip_aligned_quadrature_study_20261006/study_case.json"] == study_pins[
        "pgfrag/experiments/tip_aligned_quadrature_study_20261006/study_case.json"
    ]:
        raise RuntimeError("authorized case unexpectedly still has its pre-authorization hash")

    command = [str(PYTHON), "-B", str(ENTRY), "--execute", "--output-dir", str(OUTPUT)]
    started = time.monotonic()
    returncode: int | None = None
    timed_out = False
    launch_error: str | None = None
    stdout = b""
    stderr = b""
    try:
        completed = subprocess.run(command, cwd=ROOT, capture_output=True, timeout=60, check=False)
        returncode = completed.returncode
        stdout, stderr = completed.stdout, completed.stderr
    except subprocess.TimeoutExpired as error:
        timed_out = True
        stdout = error.stdout or b""
        stderr = error.stderr or b""
    except OSError as error:
        launch_error = f"{type(error).__name__}: {error}"
        stderr = launch_error.encode("utf-8", errors="replace")
    elapsed = time.monotonic() - started

    stdout_path = CAPTURE / "stdout.bin"
    stderr_path = CAPTURE / "stderr.bin"
    stdout_path.write_bytes(stdout)
    stderr_path.write_bytes(stderr)
    after_sources = hash_map(study_pins)
    after_inputs = hash_map(immutable_pins)
    case_after_sha = sha256(CASE)
    output_artifacts = artifacts()
    receipt = {
        "schema": "tip_aligned_study_single_attempt_capture.v1",
        "authorization": {
            "authorized_at_utc": AUTHORIZED_AT,
            "case_status_at_run": case["status"],
            "status_meaning": "source-stage status retained as history; root execution authorization is separately true",
            "pre_authorization_case_snapshot": str(PRE_CASE.relative_to(ROOT)),
            "pre_authorization_case_sha256": sha256(PRE_CASE),
            "pre_authorization_pins_snapshot": str(PRE_PINS.relative_to(ROOT)),
            "pre_authorization_pins_sha256": sha256(PRE_PINS),
            "case_sha256_at_process_start": sha256(CASE),
            "case_sha256_after_process": case_after_sha,
        },
        "process": {
            "command": command,
            "cwd": str(ROOT),
            "timeout_seconds": 60,
            "timed_out": timed_out,
            "returncode": returncode,
            "launch_error": launch_error,
            "elapsed_seconds": elapsed,
            "stdout_path": str(stdout_path.relative_to(ROOT)),
            "stdout_sha256": sha256(stdout_path),
            "stdout_bytes": stdout_path.stat().st_size,
            "stderr_path": str(stderr_path.relative_to(ROOT)),
            "stderr_sha256": sha256(stderr_path),
            "stderr_bytes": stderr_path.stat().st_size,
        },
        "source_sha256_before_process": before_sources,
        "source_sha256_after_process": after_sources,
        "immutable_input_sha256_before_process": before_inputs,
        "immutable_input_sha256_after_process": after_inputs,
        "source_pin_manifest_sha256_before_process": sha256(PINS),
        "source_pin_manifest_sha256_after_process": sha256(PINS),
        "output_path": str(OUTPUT.relative_to(ROOT)),
        "output_artifacts": output_artifacts,
        "output_bytes": sum(int(item["bytes"]) for item in output_artifacts),
    }
    receipt_path = CAPTURE / "execution_receipt.json"
    receipt_path.write_text(json.dumps(receipt, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "receipt": str(receipt_path.relative_to(ROOT)),
        "returncode": returncode,
        "timed_out": timed_out,
        "launch_error": launch_error,
        "elapsed_seconds": elapsed,
        "stdout_path": str(stdout_path.relative_to(ROOT)),
        "stdout_sha256": receipt["process"]["stdout_sha256"],
        "stderr_path": str(stderr_path.relative_to(ROOT)),
        "stderr_sha256": receipt["process"]["stderr_sha256"],
        "output_artifacts": output_artifacts,
    }, ensure_ascii=False))
    if timed_out:
        return 124
    if launch_error is not None:
        return 127
    return returncode if returncode is not None else 1


if __name__ == "__main__":
    raise SystemExit(main())
