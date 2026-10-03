"""Run the predeclared fixed-ROI wall-resolution study once per target h."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[3]
RUNNER = Path(__file__).with_name("run_roi_trial.py")
OUTPUT = ROOT / "relicweave-synlab/tmp/review_20261004/implementation/wall_resolution"
STEP = ROOT / "relicweave-synlab/tmp/review_20261003/implementation/vf/artifacts/four_bands_brep_semantics1_final/geometry/solid.step"
EXPECTED_STEP_SHA256 = "4623bd373aeed0acfacfbf6cc0dda65952339561779b43429c59d61b8b5a1c30"
MAX_TRIAL_SECONDS = 180.0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def run_one(label: str, arguments: list[str], timeout_s: float) -> dict[str, Any]:
    trial_dir = OUTPUT if label == "preflight" else OUTPUT / label
    command = [sys.executable, str(RUNNER), *arguments, "--output-dir", str(trial_dir)]
    environment = os.environ.copy()
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        environment[name] = "1"
    started = time.perf_counter()
    result: dict[str, Any] = {
        "label": label,
        "argv": command,
        "started_at_utc": utc_now(),
        "timeout_seconds": timeout_s,
        "thread_environment": {name: environment[name] for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")},
    }
    try:
        completed = subprocess.run(
            command, cwd=ROOT, env=environment, capture_output=True, text=True,
            encoding="utf-8", errors="replace", timeout=timeout_s, check=False,
        )
        result.update({
            "exit_code": int(completed.returncode),
            "timed_out": False,
            "stdout_tail": completed.stdout[-12000:],
            "stderr_tail": completed.stderr[-12000:],
        })
    except subprocess.TimeoutExpired as error:
        stdout = error.stdout or ""
        stderr = error.stderr or ""
        if isinstance(stdout, bytes):
            stdout = stdout.decode("utf-8", errors="replace")
        if isinstance(stderr, bytes):
            stderr = stderr.decode("utf-8", errors="replace")
        result.update({
            "exit_code": None,
            "timed_out": True,
            "stdout_tail": str(stdout)[-12000:],
            "stderr_tail": str(stderr)[-12000:],
        })
    result["elapsed_wall_seconds"] = float(time.perf_counter() - started)
    result["finished_at_utc"] = utc_now()
    report_path = trial_dir / ("preflight.json" if label == "preflight" else "report.json")
    result["report_path"] = str(report_path)
    result["report_sha256"] = sha256(report_path) if report_path.is_file() else None
    result["report_status"] = None
    if report_path.is_file():
        try:
            report = json.loads(report_path.read_text(encoding="utf-8"))
            result["report_status"] = report.get("status")
        except (OSError, json.JSONDecodeError):
            result["report_status"] = "UNREADABLE_REPORT"
    return result


def write_comparison_figure(trial_records: list[dict[str, Any]]) -> str | None:
    from PIL import Image, ImageDraw, ImageFont

    paths = []
    for record in trial_records:
        report_path = Path(str(record.get("report_path", "")))
        png_path = report_path.parent / f"actual_section_y0_h{record['mesh_size_mm']:g}mm.png"
        if not png_path.is_file():
            return None
        paths.append(png_path)
    images = [Image.open(path).convert("RGB") for path in paths]
    height = max(image.height for image in images)
    normalized = []
    for image in images:
        width = round(image.width * height / image.height)
        normalized.append(image.resize((width, height), Image.Resampling.LANCZOS))
    canvas = Image.new("RGB", (sum(image.width for image in normalized), height), "white")
    x = 0
    for image in normalized:
        canvas.paste(image, (x, 0))
        x += image.width
    output = OUTPUT / "actual_sections_h1mm_vs_h0p5mm.png"
    canvas.save(output, format="PNG", optimize=True)
    return str(output)


def main() -> int:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    started_at = utc_now()
    step_hash = sha256(STEP)
    if step_hash.lower() != EXPECTED_STEP_SHA256:
        raise SystemExit(f"canonical STEP SHA mismatch: expected {EXPECTED_STEP_SHA256}, got {step_hash}")

    preflight = run_one("preflight", ["--preflight-only"], MAX_TRIAL_SECONDS)
    run_records: list[dict[str, Any]] = []
    if not preflight["timed_out"] and preflight["exit_code"] == 0 and preflight["report_status"] == "PREFLIGHT_WITHIN_BOTH_FIXED_TET_BUDGETS":
        for size in (1.0, 0.5):
            label = f"h{size:g}mm"
            record = run_one(label, ["--mesh-size-mm", f"{size:g}"], MAX_TRIAL_SECONDS)
            record["mesh_size_mm"] = size
            run_records.append(record)
            if record["report_status"] and str(record["report_status"]).startswith("STOPPED_PREFLIGHT"):
                break

    comparison_png = write_comparison_figure(run_records) if len(run_records) == 2 else None
    artifacts: list[dict[str, Any]] = []
    for path in sorted(OUTPUT.rglob("*")):
        if path.is_file() and path.name != "run_index.json":
            artifacts.append({"path": str(path), "sha256": sha256(path), "bytes": path.stat().st_size})
    complete_trials = len(run_records) == 2
    all_success = complete_trials and all(
        record["exit_code"] == 0
        and not record["timed_out"]
        and record["report_status"] == "PASS_LOCAL_GEOMETRY_AND_MESH_CHECKS"
        for record in run_records
    )
    if preflight["report_status"] == "STOPPED_PREFLIGHT_EXACT_VOLUME_FORECAST_AT_OR_OVER_TET_LIMIT":
        overall_status = "STOPPED_PREFLIGHT_TET_BUDGET"
    elif all_success:
        overall_status = "BOTH_FIXED_TRIALS_COMPLETED_WITH_LOCAL_CHECKS_PASSING"
    elif complete_trials:
        overall_status = "BOTH_FIXED_TRIALS_COMPLETED_WITH_ONE_OR_MORE_CHECKS_FAILED"
    elif preflight["timed_out"]:
        overall_status = "PREFLIGHT_TIMEOUT_NO_MESH_ATTEMPT"
    elif preflight["exit_code"] != 0:
        overall_status = "PREFLIGHT_FAILED_NO_MESH_ATTEMPT"
    else:
        overall_status = "FIXED_TRIAL_SEQUENCE_STOPPED_EARLY"

    index = {
        "experiment_id": "pg.wall-resolution-local-native-cad-static-mesh/0.1.0",
        "overall_status": overall_status,
        "started_at_utc": started_at,
        "finished_at_utc": utc_now(),
        "python_executable": sys.executable,
        "argv": sys.argv,
        "runner_script": str(RUNNER),
        "runner_script_sha256": sha256(RUNNER),
        "source_step": str(STEP),
        "source_step_sha256": step_hash,
        "expected_source_step_sha256": EXPECTED_STEP_SHA256,
        "fixed_scope": {
            "roi_box_mm": {"xmin": -23.0, "xmax": -14.0, "ymin": -1.4, "ymax": 1.4, "zmin": 23.1, "zmax": 25.9},
            "target_h_mm": [1.0, 0.5],
            "maximum_volume_tets_per_trial": 10000,
            "maximum_wall_seconds_per_trial": MAX_TRIAL_SECONDS,
            "single_thread": True,
            "stop_when_either_exact_v_forecast_reaches_budget": True,
            "no_roi_h_budget_or_threshold_adaptation": True,
        },
        "preflight": preflight,
        "trials": run_records,
        "comparison_section_png": comparison_png,
        "artifacts_hashed_before_run_index_write": artifacts,
    }
    index_path = OUTPUT / "run_index.json"
    index_path.write_text(json.dumps(index, indent=2), encoding="utf-8")
    print(json.dumps({
        "overall_status": overall_status,
        "run_index": str(index_path),
        "run_index_sha256": sha256(index_path),
        "preflight": {"exit_code": preflight["exit_code"], "status": preflight["report_status"], "elapsed_wall_seconds": preflight["elapsed_wall_seconds"]},
        "trials": [{"h_mm": record["mesh_size_mm"], "exit_code": record["exit_code"], "timed_out": record["timed_out"], "status": record["report_status"], "elapsed_wall_seconds": record["elapsed_wall_seconds"], "report_sha256": record["report_sha256"]} for record in run_records],
        "comparison_section_png": comparison_png,
    }, indent=2))
    return 0 if all_success else 2 if run_records else 3


if __name__ == "__main__":
    raise SystemExit(main())
