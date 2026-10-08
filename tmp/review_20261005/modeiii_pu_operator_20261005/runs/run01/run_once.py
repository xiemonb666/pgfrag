"""Run the root-authorized frozen Mode III operator diagnostic exactly once."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import time


REPO = Path(r"C:\Users\xiemo\Desktop\ShardForge")
RUN_DIR = Path(__file__).resolve().parent
REVIEW_DIR = REPO / "pgfrag" / "tmp" / "review_20261005" / "modeiii_pu_operator_20261005"
SOURCE_DIR = REPO / "pgfrag" / "experiments" / "modeiii_pu_operator_20261005"
FREEZE_PATH = REVIEW_DIR / "source_freeze.json"
DIAGNOSTIC_PATH = SOURCE_DIR / "modeiii_operator_diagnostic.py"
CASE_PATH = SOURCE_DIR / "fixed_case.json"
PRIMARY_PATH = RUN_DIR / "operator_diagnostic.json"
APPROVED_FREEZE_SHA256 = "65197acaead31faacfad967ffda94ebac9991ff7a433c86fc133e699ad591780"
WALL_BUDGET_SECONDS = 30.0
RUNNER_TIMEOUT_SECONDS = 25.0
POST_RUN_RESERVE_SECONDS = 4.0
MAX_ARTIFACT_BYTES = 4 * 1024 * 1024
EXPECTED_RUNTIME = (REPO / "third_party" / "vesselforge-runtime" / "venv" / "Scripts" / "python.exe").resolve()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def json_write(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def frozen_paths(manifest: dict[str, object]) -> list[tuple[str, str, str]]:
    entries: list[tuple[str, str, str]] = []
    for section, role in (("owned_sources", "owned_source"), ("read_only_dependencies", "read_only_dependency")):
        values = manifest[section]
        if not isinstance(values, dict):
            raise ValueError(f"freeze manifest section {section} is malformed")
        entries.extend((role, str(path), str(expected).lower()) for path, expected in values.items())
    entries.append(("freeze_manifest", str(FREEZE_PATH.relative_to(REPO)).replace("\\", "/"), APPROVED_FREEZE_SHA256))
    return entries


def take_pins(entries: list[tuple[str, str, str]], phase: str) -> list[dict[str, object]]:
    pins: list[dict[str, object]] = []
    for role, relative_path, expected in entries:
        path = REPO / Path(relative_path)
        observed = sha256(path)
        pins.append(
            {
                "role": role,
                "path": relative_path,
                "sha256": observed,
                "expected_sha256": expected,
                "matches_approved_pin": observed == expected,
                "phase": phase,
            }
        )
    return pins


def decode_stream(value: str | bytes | None) -> str:
    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode("utf-8", "replace")
    return value


def quadrature_coordinates_for_indices(indices: set[int]) -> dict[str, object]:
    """Reconstruct only reported fixed-case query coordinates, without repair."""
    import numpy as np

    case = json.loads(CASE_PATH.read_text(encoding="utf-8"))
    rules = case["gauss_legendre_rules"]
    geometry = case["geometry"]
    radius = float(geometry["radius_m"])
    length = float(geometry["length_m"])
    half_length = length / 2.0

    def rule(order: int) -> tuple[object, object]:
        record = rules[f"order_{order}"]
        return np.asarray(record["nodes"], dtype=np.float64), np.asarray(record["weights"], dtype=np.float64)

    def mapped(nodes: object, weights: object, lower: float, upper: float) -> tuple[object, object]:
        half = 0.5 * (upper - lower)
        middle = 0.5 * (upper + lower)
        return middle + half * nodes, half * weights

    r4, wr4 = rule(4)
    theta8, wtheta8 = rule(8)
    z4, wz4 = rule(4)
    radial, radial_weights = mapped(r4, wr4, 0.0, radius)
    theta, _ = mapped(theta8, wtheta8, -math.pi, math.pi)
    axial, _ = mapped(z4, wz4, -half_length, half_length)
    points: list[tuple[float, float, float]] = []
    for ir in range(4):
        for it in range(8):
            for iz in range(4):
                r = float(radial[ir])
                angle = float(theta[it])
                points.append((r * math.cos(angle), r * math.sin(angle), float(axial[iz])))

    theta16, wtheta16 = rule(16)
    side_theta, _ = mapped(theta16, wtheta16, -math.pi, math.pi)
    side_axial, _ = mapped(z4, wz4, -half_length, half_length)
    for it in range(16):
        for iz in range(4):
            angle = float(side_theta[it])
            cosine, sine = math.cos(angle), math.sin(angle)
            points.append((radius * cosine, radius * sine, float(side_axial[iz])))

    for z_value in (half_length, -half_length):
        for ir in range(4):
            for it in range(8):
                r = float(radial[ir])
                angle = float(theta[it])
                points.append((r * math.cos(angle), r * math.sin(angle), z_value))

    rows = []
    for index in sorted(indices):
        if index < 0 or index >= len(points):
            rows.append({"query_index": index, "coordinate_status": "INDEX_OUT_OF_RANGE"})
        else:
            rows.append(
                {
                    "query_index": index,
                    "coordinate_xyz_m": [float(value) for value in points[index]],
                    "coordinate_source": "reconstructed from fixed GL tables and original fixed-case mapping; no clipping or nudging",
                }
            )
    return {
        "schema": "modeiii_pu_operator_unknown_query_coordinates.v1",
        "source_freeze_sha256": APPROVED_FREEZE_SHA256,
        "query_order": ["volume", "lateral", "top_cap", "bottom_cap"],
        "unknown_queries": rows,
        "note": "The primary diagnostic JSON and its query indices remain unchanged; this is a coordinate receipt only.",
    }


def main() -> int:
    outer_started = time.monotonic()
    start_utc = utc_now()
    protected_names = (
        "pre_run_pins.json",
        "post_run_pins.json",
        "stdout.txt",
        "stderr.txt",
        "operator_diagnostic.json",
        "run_record.json",
        "unknown_query_coordinates.json",
    )
    existing = [name for name in protected_names if (RUN_DIR / name).exists()]
    if existing:
        raise SystemExit(f"refusing to overwrite pre-existing run artifacts: {existing}")

    freeze_sha = sha256(FREEZE_PATH)
    manifest = json.loads(FREEZE_PATH.read_text(encoding="utf-8"))
    entries = frozen_paths(manifest)
    pre_pins = take_pins(entries, "pre_run")
    pin_match = freeze_sha == APPROVED_FREEZE_SHA256 and all(pin["matches_approved_pin"] for pin in pre_pins)
    json_write(
        RUN_DIR / "pre_run_pins.json",
        {
            "started_at_utc": start_utc,
            "source_freeze_manifest_sha256": freeze_sha,
            "approved_source_freeze_manifest_sha256": APPROVED_FREEZE_SHA256,
            "all_pins_match": pin_match,
            "pins": pre_pins,
        },
    )

    argv = [str(Path(sys.executable).resolve()), "-B", str(DIAGNOSTIC_PATH), "--output", str(PRIMARY_PATH)]
    runtime_match = Path(sys.executable).resolve().as_posix().casefold() == EXPECTED_RUNTIME.as_posix().casefold()
    stdout_text = ""
    stderr_text = ""
    return_code: int | None = None
    timed_out = False
    timeout_taskkill: dict[str, object] | None = None
    launched = False
    diagnostic_elapsed = 0.0
    if pin_match and runtime_match:
        elapsed_before_launch = time.monotonic() - outer_started
        child_timeout = min(
            RUNNER_TIMEOUT_SECONDS,
            max(0.0, WALL_BUDGET_SECONDS - elapsed_before_launch - POST_RUN_RESERVE_SECONDS),
        )
        if child_timeout <= 0:
            stderr_text = "No diagnostic process launched: insufficient remaining wall budget.\n"
            return_code = 124
            timed_out = True
        else:
            launch_time = time.monotonic()
            launched = True
            process = subprocess.Popen(
                argv,
                cwd=str(REPO),
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                encoding="utf-8",
                errors="replace",
                creationflags=getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0),
            )
            try:
                stdout_text, stderr_text = process.communicate(timeout=child_timeout)
                return_code = int(process.returncode)
            except subprocess.TimeoutExpired as exc:
                timed_out = True
                stdout_text = decode_stream(exc.stdout)
                stderr_text = decode_stream(exc.stderr)
                taskkill_argv = ["taskkill.exe", "/PID", str(process.pid), "/T", "/F"]
                try:
                    taskkill = subprocess.run(
                        taskkill_argv,
                        capture_output=True,
                        text=True,
                        encoding="utf-8",
                        errors="replace",
                        timeout=1.0,
                        check=False,
                    )
                    timeout_taskkill = {
                        "argv": taskkill_argv,
                        "return_code": taskkill.returncode,
                        "stdout": taskkill.stdout,
                        "stderr": taskkill.stderr,
                    }
                except subprocess.TimeoutExpired as taskkill_exc:
                    timeout_taskkill = {
                        "argv": taskkill_argv,
                        "return_code": 124,
                        "stdout": decode_stream(taskkill_exc.stdout),
                        "stderr": decode_stream(taskkill_exc.stderr),
                    }
                try:
                    stdout_text, stderr_text = process.communicate(timeout=0.5)
                    return_code = int(process.returncode) if process.returncode is not None else 124
                except subprocess.TimeoutExpired:
                    return_code = 124
                    stderr_text += "\nOuter diagnostic deadline reached; process tree termination was requested.\n"
            diagnostic_elapsed = time.monotonic() - launch_time
    elif not pin_match:
        stderr_text = "Approved freeze/source pin mismatch; diagnostic was not started.\n"
        return_code = 125
    else:
        stderr_text = f"Expected VF runtime {EXPECTED_RUNTIME}, got {Path(sys.executable).resolve()}; diagnostic was not started.\n"
        return_code = 126

    (RUN_DIR / "stdout.txt").write_text(stdout_text, encoding="utf-8")
    (RUN_DIR / "stderr.txt").write_text(stderr_text, encoding="utf-8")

    primary_summary: dict[str, object] | None = None
    primary_sha: str | None = None
    primary_bytes: int | None = None
    if PRIMARY_PATH.exists():
        primary_sha = sha256(PRIMARY_PATH)
        primary_bytes = PRIMARY_PATH.stat().st_size
        try:
            result = json.loads(PRIMARY_PATH.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            primary_summary = {"parse_error": f"{type(exc).__name__}: {exc}"}
        else:
            primary_summary = {
                "status": result.get("status"),
                "case_id": result.get("case_id"),
                "dimensions": result.get("dimensions"),
                "mls": result.get("mls"),
                "manufactured_modeiii_field": result.get("manufactured_modeiii_field"),
                "energy_and_work": result.get("energy_and_work"),
                "operator_diagnostics": result.get("operator_diagnostics"),
                "parent_elapsed_seconds": result.get("parent_elapsed_seconds"),
                "parent_within_wall_timeout": result.get("parent_within_wall_timeout"),
                "worker_dependency_sha256": result.get("worker_dependency_sha256"),
            }
            unknown_rows = (result.get("mls") or {}).get("unknown_rows", [])
            unknown_ids = {
                int(row["query_index"])
                for row in unknown_rows
                if isinstance(row, dict) and isinstance(row.get("query_index"), int)
            }
            if unknown_ids:
                try:
                    coordinate_receipt = quadrature_coordinates_for_indices(unknown_ids)
                    json_write(RUN_DIR / "unknown_query_coordinates.json", coordinate_receipt)
                    primary_summary["unknown_coordinate_receipt"] = "unknown_query_coordinates.json"
                except Exception as exc:
                    primary_summary["unknown_coordinate_receipt_error"] = f"{type(exc).__name__}: {exc}"

    post_pins = take_pins(entries, "post_run")
    post_unchanged = all(
        before["sha256"] == after["sha256"]
        for before, after in zip(pre_pins, post_pins, strict=True)
    )
    json_write(
        RUN_DIR / "post_run_pins.json",
        {
            "finished_at_utc": utc_now(),
            "all_pins_unchanged": post_unchanged,
            "pins": [
                {
                    "role": after["role"],
                    "path": after["path"],
                    "pre_sha256": before["sha256"],
                    "post_sha256": after["sha256"],
                    "matches_approved_pin": after["matches_approved_pin"],
                    "unchanged_during_run": before["sha256"] == after["sha256"],
                }
                for before, after in zip(pre_pins, post_pins, strict=True)
            ],
        },
    )

    final_elapsed = time.monotonic() - outer_started
    if not pin_match:
        status = "PIN_MISMATCH_NO_RUN"
    elif timed_out:
        status = "TIMED_OUT_NO_RETRY"
    elif primary_summary is None:
        status = "COMPLETED_NO_PRIMARY_JSON"
    else:
        diagnostic_status = primary_summary.get("status")
        status = "COMPLETED_ONCE_" + str(diagnostic_status)
    record = {
        "schema": "modeiii_pu_operator_single_run_receipt.v1",
        "run_id": "modeiii_pu_operator_run01_20261005",
        "status": status,
        "authorization": "One fixed-case Mode III PU operator diagnostic explicitly authorized by root after source freeze review.",
        "source_freeze_manifest_sha256": freeze_sha,
        "approved_source_freeze_manifest_sha256": APPROVED_FREEZE_SHA256,
        "runtime": str(Path(sys.executable).resolve()),
        "runtime_matches_authorized_vf": runtime_match,
        "python_no_bytecode_flag": "-B",
        "working_directory": str(REPO),
        "argv": argv,
        "raw_process_arguments": subprocess.list2cmdline(argv[1:]),
        "launched": launched,
        "started_at_utc": start_utc,
        "finished_at_utc": utc_now(),
        "diagnostic_process_elapsed_seconds": diagnostic_elapsed,
        "outer_elapsed_seconds_including_pin_checks_and_receipt": final_elapsed,
        "outer_hard_wall_budget_seconds": WALL_BUDGET_SECONDS,
        "runner_timeout_seconds": RUNNER_TIMEOUT_SECONDS,
        "post_run_reserve_seconds": POST_RUN_RESERVE_SECONDS,
        "within_outer_wall_budget": final_elapsed <= WALL_BUDGET_SECONDS,
        "timed_out": timed_out,
        "no_retry": True,
        "real_return_code": return_code,
        "outer_timeout_process_tree_action": timeout_taskkill,
        "primary_json_path": str(PRIMARY_PATH),
        "primary_json_sha256": primary_sha,
        "primary_json_bytes": primary_bytes,
        "primary_summary": primary_summary,
        "stdout_path": "stdout.txt",
        "stdout_sha256": sha256(RUN_DIR / "stdout.txt"),
        "stdout_bytes": (RUN_DIR / "stdout.txt").stat().st_size,
        "stderr_path": "stderr.txt",
        "stderr_sha256": sha256(RUN_DIR / "stderr.txt"),
        "stderr_bytes": (RUN_DIR / "stderr.txt").stat().st_size,
        "pre_run_pins_path": "pre_run_pins.json",
        "post_run_pins_path": "post_run_pins.json",
        "pre_run_all_pins_match": pin_match,
        "post_run_all_pins_unchanged": post_unchanged,
        "pin_count": len(entries),
        "artifact_bytes_including_this_receipt": 0,
        "artifact_cap_bytes": MAX_ARTIFACT_BYTES,
        "within_artifact_cap": True,
    }
    for _ in range(3):
        json_write(RUN_DIR / "run_record.json", record)
        artifact_bytes = sum(path.stat().st_size for path in RUN_DIR.iterdir() if path.is_file())
        if record["artifact_bytes_including_this_receipt"] == artifact_bytes:
            break
        record["artifact_bytes_including_this_receipt"] = artifact_bytes
        record["within_artifact_cap"] = artifact_bytes <= MAX_ARTIFACT_BYTES
    json_write(RUN_DIR / "run_record.json", record)
    print(json.dumps(record, ensure_ascii=False, allow_nan=False))
    return int(return_code or 0)


if __name__ == "__main__":
    raise SystemExit(main())
