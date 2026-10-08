"""One-shot bounded launcher for the root-authorized read-only audit."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time
from datetime import datetime, timezone


WORKSPACE = Path.cwd().resolve()
REVIEW = WORKSPACE / "pgfrag/tmp/review_20261005/forge_volume_ingress_independent_review"
FREEZE_PATH = REVIEW / "PILOT_AUDIT_INPUT_FREEZE.json"
ADAPTER = REVIEW / "pilot_actual_audit.py"
OUTPUT_DIR = REVIEW / "audit_retry_hashcase_plotcompat_20261005"
ANALYTIC = REVIEW / "analytic_roi_reference.py"
METHOD = REVIEW / "PILOT_AUDIT_ADAPTER_REVIEW.md"
OLD_FREEZE = REVIEW / "PILOT_AUDIT_INPUT_FREEZE_before_centroid_gate.json"
FAILED_FREEZE = REVIEW / "PILOT_AUDIT_INPUT_FREEZE_before_partial_pin_correction.json"
FAILED_RECEIPT = REVIEW / "READONLY_AUDIT_EXECUTION_RECEIPT.json"
CASE_RECEIPT = REVIEW / "READONLY_AUDIT_PREFLIGHT_MISMATCH_RECEIPT.json"
PRECANON_FREEZE = REVIEW / "PILOT_AUDIT_INPUT_FREEZE_before_sha_lowercase_canonicalization.json"
OLD_ADAPTER = REVIEW / "pilot_actual_audit_before_hashcase_fix_06F1A1D3E64EDB3F950C260F673999048FFA2584025EE77FD922AF99A09F92B5.py"
FAILED_AUDIT_RECEIPT = REVIEW / "READONLY_AUDIT_EXECUTION_RECEIPT_corrected_freeze.json"
FAILED_AUDIT_STDOUT = REVIEW / "READONLY_AUDIT_STDOUT_corrected_freeze.txt"
FAILED_AUDIT_STDERR = REVIEW / "READONLY_AUDIT_STDERR_corrected_freeze.txt"
VALIDATION_REV2 = REVIEW / "PILOT_AUDIT_FREEZE_VALIDATION_rev2.json"
VALIDATION_REV3 = REVIEW / "PILOT_AUDIT_FREEZE_VALIDATION_rev3.json"
VALIDATION_REV4 = REVIEW / "PILOT_AUDIT_FREEZE_VALIDATION_rev4.json"
CANONICALIZATION_RECEIPT = REVIEW / "PILOT_AUDIT_FREEZE_CANONICALIZATION_rev3.json"
RENDERFIX_FREEZE_RECEIPT = REVIEW / "PILOT_AUDIT_FREEZE_RENDERFIX_rev4.json"
VALIDATION = REVIEW / "PILOT_AUDIT_FREEZE_VALIDATION_rev5.json"
PREFLIGHT_RECEIPT = REVIEW / "PILOT_AUDIT_LAUNCHER_PREFLIGHT_rev5.json"
PREVIOUS_CELL_CSV = REVIEW / "PILOT_CELL_AUDIT.csv"
PREVIOUS_POINT_CSV = REVIEW / "PILOT_POINT_AUDIT.csv"
PREVIOUS_RETRY_RECEIPT = REVIEW / "READONLY_AUDIT_EXECUTION_RECEIPT_canonicalized_freeze.json"
PREVIOUS_RETRY_STDOUT = REVIEW / "READONLY_AUDIT_STDOUT_canonicalized_freeze.txt"
PREVIOUS_RETRY_STDERR = REVIEW / "READONLY_AUDIT_STDERR_canonicalized_freeze.txt"
RUNTIME = Path(
    "C:/Users/xiemo/Desktop/ShardForge/third_party/"
    "vesselforge-runtime/venv/Scripts/python.exe"
)
RECEIPT = REVIEW / "READONLY_AUDIT_EXECUTION_RECEIPT_renderreport_retry.json"
STDOUT_PATH = REVIEW / "READONLY_AUDIT_STDOUT_renderreport_retry.txt"
STDERR_PATH = REVIEW / "READONLY_AUDIT_STDERR_renderreport_retry.txt"
OUTPUT_NAMES = (
    "PILOT_ACTUAL_AUDIT.json",
    "PILOT_CELL_AUDIT.csv",
    "PILOT_POINT_AUDIT.csv",
    "PILOT_XY_CELL_MAP.png",
    "PILOT_CELL_VOLUME_MASS_DELTA.png",
)
TIMEOUT_SECONDS = 55.0


def sha256(path: Path) -> str | None:
    if not path.is_file():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, value: dict) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


freeze = json.loads(FREEZE_PATH.read_text(encoding="utf-8"))
stage = WORKSPACE / freeze["input_scope"]["retained_stage"]
summary = WORKSPACE / freeze["input_scope"]["supervisor_summary"]


def collect_hashes() -> dict:
    stage_rows = {}
    if stage.is_dir():
        for path in sorted(stage.iterdir(), key=lambda item: item.name):
            if path.is_file():
                stage_rows[path.name] = {
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
    return {
        "launcher": {"bytes": Path(__file__).stat().st_size, "sha256": sha256(Path(__file__).resolve())},
        "adapter": {"bytes": ADAPTER.stat().st_size, "sha256": sha256(ADAPTER)},
        "analytic_source": {
            "bytes": ANALYTIC.stat().st_size,
            "sha256": sha256(ANALYTIC),
        },
        "method_review": {"bytes": METHOD.stat().st_size, "sha256": sha256(METHOD)},
        "freeze": {
            "bytes": FREEZE_PATH.stat().st_size,
            "sha256": sha256(FREEZE_PATH),
        },
        "prior_freeze": {
            "bytes": OLD_FREEZE.stat().st_size,
            "sha256": sha256(OLD_FREEZE),
        },
        "preserved_preflight_artifacts": {
            "malformed_freeze": {"bytes": FAILED_FREEZE.stat().st_size, "sha256": sha256(FAILED_FREEZE)},
            "precanonical_freeze": {"bytes": PRECANON_FREEZE.stat().st_size, "sha256": sha256(PRECANON_FREEZE)},
            "preplotcompat_freeze": {"bytes": (REVIEW / "PILOT_AUDIT_INPUT_FREEZE_before_plotcompat_20261005.json").stat().st_size, "sha256": sha256(REVIEW / "PILOT_AUDIT_INPUT_FREEZE_before_plotcompat_20261005.json")},
            "old_adapter": {"bytes": OLD_ADAPTER.stat().st_size, "sha256": sha256(OLD_ADAPTER)},
            "preretry_adapter": {"bytes": (REVIEW / "pilot_actual_audit_pre_retry_refactor_827F429C141403A5E85C7C9A67376D38B0628F9904BF59A11E8DCFE883E2EC32.py").stat().st_size, "sha256": sha256(REVIEW / "pilot_actual_audit_pre_retry_refactor_827F429C141403A5E85C7C9A67376D38B0628F9904BF59A11E8DCFE883E2EC32.py")},
            "failed_receipt": {"bytes": FAILED_RECEIPT.stat().st_size, "sha256": sha256(FAILED_RECEIPT)},
            "case_receipt": {"bytes": CASE_RECEIPT.stat().st_size, "sha256": sha256(CASE_RECEIPT)},
            "failed_audit_receipt": {"bytes": FAILED_AUDIT_RECEIPT.stat().st_size, "sha256": sha256(FAILED_AUDIT_RECEIPT)},
            "failed_audit_stdout": {"bytes": FAILED_AUDIT_STDOUT.stat().st_size, "sha256": sha256(FAILED_AUDIT_STDOUT)},
            "failed_audit_stderr": {"bytes": FAILED_AUDIT_STDERR.stat().st_size, "sha256": sha256(FAILED_AUDIT_STDERR)},
            "freeze_validation_rev2": {"bytes": VALIDATION_REV2.stat().st_size, "sha256": sha256(VALIDATION_REV2)},
            "freeze_validation_rev3": {"bytes": VALIDATION_REV3.stat().st_size, "sha256": sha256(VALIDATION_REV3)},
            "freeze_validation_rev4": {"bytes": VALIDATION_REV4.stat().st_size, "sha256": sha256(VALIDATION_REV4)},
            "freeze_validation_rev5": {"bytes": VALIDATION.stat().st_size, "sha256": sha256(VALIDATION)},
            "canonicalization_receipt": {"bytes": CANONICALIZATION_RECEIPT.stat().st_size, "sha256": sha256(CANONICALIZATION_RECEIPT)},
            "renderfix_freeze_receipt": {"bytes": RENDERFIX_FREEZE_RECEIPT.stat().st_size, "sha256": sha256(RENDERFIX_FREEZE_RECEIPT)},
            "previous_renderretry_receipt": {"bytes": PREVIOUS_RETRY_RECEIPT.stat().st_size, "sha256": sha256(PREVIOUS_RETRY_RECEIPT)},
            "previous_renderretry_stdout": {"bytes": PREVIOUS_RETRY_STDOUT.stat().st_size, "sha256": sha256(PREVIOUS_RETRY_STDOUT)},
            "previous_renderretry_stderr": {"bytes": PREVIOUS_RETRY_STDERR.stat().st_size, "sha256": sha256(PREVIOUS_RETRY_STDERR)},
            "previous_partial_cell_csv": {"bytes": PREVIOUS_CELL_CSV.stat().st_size, "sha256": sha256(PREVIOUS_CELL_CSV)},
            "previous_partial_point_csv": {"bytes": PREVIOUS_POINT_CSV.stat().st_size, "sha256": sha256(PREVIOUS_POINT_CSV)},
        },
        "summary": {
            "bytes": summary.stat().st_size if summary.is_file() else None,
            "sha256": sha256(summary),
        },
        "stage_files": stage_rows,
    }


def collect_sha_pins(value: object, path: str = "$") -> list[tuple[str, object]]:
    found: list[tuple[str, object]] = []
    if isinstance(value, dict):
        for key, item in value.items():
            child = f"{path}.{key}"
            if key == "sha256":
                found.append((child, item))
            else:
                found.extend(collect_sha_pins(item, child))
    elif isinstance(value, list):
        for index, item in enumerate(value):
            found.extend(collect_sha_pins(item, f"{path}[{index}]"))
    return found


pre_hashes = collect_hashes()
expected_stage_names = sorted(freeze["retained_stage_files"])
pin_errors = []
for pin_path, pin_value in collect_sha_pins(freeze):
    normalized = pin_value.lower() if isinstance(pin_value, str) else ""
    if re.fullmatch(r"[0-9a-f]{64}", normalized) is None:
        pin_errors.append(f"invalid 64-hex SHA-256 pin: {pin_path}")
if sorted(pre_hashes["stage_files"]) != expected_stage_names:
    pin_errors.append("retained stage file inventory differs from freeze")
for name, expected in freeze["retained_stage_files"].items():
    actual = pre_hashes["stage_files"].get(name)
    if (
        actual is None
        or actual["bytes"] != expected["bytes"]
        or actual["sha256"].upper() != expected["sha256"].upper()
    ):
        pin_errors.append(f"retained stage pin mismatch: {name}")
for actual_key, expected_key in (
    ("summary", "supervisor_summary"),
    ("analytic_source", "analytic_source"),
    ("adapter", "audit_adapter"),
    ("method_review", "method_review"),
    ("prior_freeze", "prior_source_snapshot"),
):
    actual = pre_hashes[actual_key]
    expected = freeze[expected_key]
    if actual["bytes"] != expected["bytes"] if "bytes" in expected else False:
        pin_errors.append(f"size pin mismatch: {actual_key}")
    if actual["sha256"].upper() != expected["sha256"].upper():
        pin_errors.append(f"SHA-256 pin mismatch: {actual_key}")

existing_outputs = [name for name in OUTPUT_NAMES if (OUTPUT_DIR / name).exists()]
if existing_outputs:
    pin_errors.append("refusing to overwrite pre-existing audit outputs: "
                       + ", ".join(existing_outputs))
existing_receipts = [
    path.name for path in (RECEIPT, STDOUT_PATH, STDERR_PATH) if path.exists()
]
if existing_receipts:
    pin_errors.append("refusing to overwrite pre-existing execution records: "
                       + ", ".join(existing_receipts))

launcher_argv = [str(RUNTIME), "-B", str(Path(__file__).resolve())]
audit_argv = [str(RUNTIME), "-B", str(ADAPTER)]
if "--preflight-only" in sys.argv:
    preflight_record = {
        "schema": "pgfrag.forge_volume_ingress_launcher_preflight.v1",
        "status": "PREFLIGHT_PASS_NOT_EXECUTED" if not pin_errors else "PREFLIGHT_FAIL_NOT_EXECUTED",
        "audit_invoked": False,
        "launcher_argv": launcher_argv + ["--preflight-only"],
        "launcher_sha256": sha256(Path(__file__).resolve()),
        "freeze_sha256": pre_hashes["freeze"]["sha256"],
        "source_and_input_hashes_pre": pre_hashes,
        "preflight_errors": pin_errors,
        "new_actual_execution_receipt_exists": RECEIPT.exists(),
        "audit_output_dir": str(OUTPUT_DIR),
        "audit_outputs_exist": {name: (OUTPUT_DIR / name).exists() for name in OUTPUT_NAMES},
    }
    write_json(PREFLIGHT_RECEIPT, preflight_record)
    print(json.dumps({
        "status": preflight_record["status"],
        "audit_invoked": False,
        "freeze_sha256": pre_hashes["freeze"]["sha256"],
        "preflight_errors": pin_errors,
        "receipt": PREFLIGHT_RECEIPT.name,
    }, sort_keys=True))
    raise SystemExit(0 if not pin_errors else 3)

started_utc = utc_now()
started = time.perf_counter()
stdout = ""
stderr = ""
returncode = None
timed_out = False
run_status = "NOT_RUN_PREFLIGHT_FAILED" if pin_errors else "PENDING"

if not pin_errors:
    child_env = os.environ.copy()
    child_env["MPLCONFIGDIR"] = str(OUTPUT_DIR / ".mplconfig")
    try:
        child = subprocess.run(
            audit_argv,
            cwd=str(WORKSPACE),
            env=child_env,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=TIMEOUT_SECONDS,
            check=False,
        )
        stdout = child.stdout
        stderr = child.stderr
        returncode = child.returncode
        run_status = "COMPLETED"
    except subprocess.TimeoutExpired as exc:
        timed_out = True
        run_status = "TIMED_OUT"
        stdout = exc.stdout or ""
        stderr = exc.stderr or ""
        if isinstance(stdout, bytes):
            stdout = stdout.decode("utf-8", errors="replace")
        if isinstance(stderr, bytes):
            stderr = stderr.decode("utf-8", errors="replace")

elapsed_seconds = time.perf_counter() - started
finished_utc = utc_now()
post_hashes = collect_hashes()

STDOUT_PATH.write_text(stdout, encoding="utf-8", newline="")
STDERR_PATH.write_text(stderr, encoding="utf-8", newline="")
output_hashes = {}
for name in OUTPUT_NAMES:
    path = OUTPUT_DIR / name
    output_hashes[name] = {
        "exists": path.is_file(),
        "bytes": path.stat().st_size if path.is_file() else None,
        "sha256": sha256(path),
    }

unchanged_keys = (
    "launcher", "adapter", "analytic_source", "method_review", "freeze", "prior_freeze",
    "preserved_preflight_artifacts", "summary", "stage_files",
)
inputs_unchanged = all(
    pre_hashes[key] == post_hashes[key] for key in unchanged_keys
)
receipt = {
    "schema": "pgfrag.forge_volume_ingress_readonly_audit_execution_receipt.v1",
    "status": run_status,
    "started_utc": started_utc,
    "finished_utc": finished_utc,
    "elapsed_seconds": elapsed_seconds,
    "timeout_seconds": TIMEOUT_SECONDS,
    "launcher_argv": launcher_argv,
    "launcher_sha256": sha256(Path(__file__).resolve()),
    "audit_argv": audit_argv,
    "audit_command_line": subprocess.list2cmdline(audit_argv),
    "working_directory": str(WORKSPACE),
    "audit_output_directory": str(OUTPUT_DIR),
    "runtime": str(RUNTIME),
    "return_code": returncode,
    "timed_out": timed_out,
    "preflight_errors": pin_errors,
    "source_and_input_hashes_pre": pre_hashes,
    "source_and_input_hashes_post": post_hashes,
    "source_and_inputs_unchanged": inputs_unchanged,
    "stdout_file": {
        "path": STDOUT_PATH.name,
        "bytes": STDOUT_PATH.stat().st_size,
        "sha256": sha256(STDOUT_PATH),
    },
    "stderr_file": {
        "path": STDERR_PATH.name,
        "bytes": STDERR_PATH.stat().st_size,
        "sha256": sha256(STDERR_PATH),
    },
    "outputs": output_hashes,
    "producer_or_cad_or_solver_invoked_by_launcher": False,
}
write_json(RECEIPT, receipt)
print(json.dumps({
    "status": run_status,
    "return_code": returncode,
    "timed_out": timed_out,
    "elapsed_seconds": elapsed_seconds,
    "preflight_errors": pin_errors,
    "audit_output_directory": str(OUTPUT_DIR),
    "source_and_inputs_unchanged": inputs_unchanged,
    "receipt": RECEIPT.name,
    "outputs": output_hashes,
}, sort_keys=True))

if pin_errors:
    raise SystemExit(3)
if timed_out:
    raise SystemExit(124)
if returncode != 0:
    raise SystemExit(returncode if returncode is not None else 1)
