"""Run exactly one bounded native zero-load vessel benchmark under WSL."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parents[4]
PG = ROOT / "pgfrag"
EXPERIMENT = PG / "experiments/continuation_20261016/reference_vessel_benchmark"
DATA = PG / "data/continuation_20261016/reference_vessel_benchmark"
INPUTS = EXPERIMENT / "inputs"
CASE = INPUTS / "case.json"
RUN = EXPERIMENT / "native/run_attempt_01"
STDOUT = EXPERIMENT / "native_stdout.log"
STDERR = EXPERIMENT / "native_stderr.log"
RECORD = EXPERIMENT / "execution_record.json"
PYTHON = Path("/mnt/c/Users/xiemo/Desktop/ShardForge/third_party/akantu-runtime/venv/bin/python")
OUTER_WALL_SECONDS = 150


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def input_hashes() -> dict[str, str]:
    return {path.name: sha256(path) for path in sorted(INPUTS.iterdir()) if path.is_file()}


def matching_live_processes() -> list[dict]:
    """Inspect /proc command lines for a live wrapper using this exact case."""
    matches: list[dict] = []
    target = os.fsencode(str(CASE))
    for entry in Path("/proc").iterdir():
        if not entry.name.isdigit():
            continue
        try:
            raw = (entry / "cmdline").read_bytes()
            stat = (entry / "stat").read_text(encoding="ascii")
        except (OSError, UnicodeError):
            continue
        args = raw.split(b"\0")
        if target in args and any(b"cohesive_vessel.py" in item for item in args):
            try:
                state = stat.rsplit(")", 1)[1].strip().split()[0]
            except (IndexError, ValueError):
                state = "?"
            matches.append({"pid": int(entry.name), "state": state,
                            "argv": [arg.decode("utf-8", errors="replace") for arg in args if arg]})
    return matches


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def main() -> int:
    if not CASE.is_file() or not (DATA / "preflight_budget.json").is_file():
        raise FileNotFoundError("run preparation and exact budget gate must pass first")
    if RUN.exists() or RECORD.exists() or STDOUT.exists() or STDERR.exists():
        raise FileExistsError("one attempt has already started or has preserved output; refusing retry")
    budget = json.loads((DATA / "preflight_budget.json").read_text(encoding="utf-8"))
    case = json.loads(CASE.read_text(encoding="utf-8"))
    if budget.get("status") != "EXACT_HISTORY_AND_ALL_OUTPUT_BUDGET_PASS; RUN_NOT_YET_STARTED":
        raise ValueError("preflight forecast did not pass its caps")
    if not budget["exact_history_forecast"]["history_cap_check_passes"] or not budget["all_output_budget"]["check_passes"]:
        raise ValueError("history or total-output budget rejects this attempt")
    if float(case["loading"]["peak_force_N"]) != 0.0:
        raise ValueError("this authorized benchmark is zero-load only")

    active = matching_live_processes()
    if active:
        raise RuntimeError(f"a process already uses this exact benchmark case; no new process launched: {active}")

    env = os.environ.copy()
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLAS_NUM_THREADS",
                 "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env[name] = "1"
    env["OMP_DYNAMIC"] = "FALSE"
    argv = [str(PYTHON), str(INPUTS / "cohesive_vessel.py"), "--case", str(CASE), "--out", str(RUN)]
    before = input_hashes()
    if before != {name: row["sha256"] for name, row in budget["input_freeze"]["frozen_files"].items()}:
        raise ValueError("frozen input hash verification failed before run")

    started_wall = time.time()
    started_monotonic = time.monotonic()
    timed_out = False
    returncode: int | None = None
    with STDOUT.open("wb") as out, STDERR.open("wb") as err:
        process = subprocess.Popen(
            argv, cwd=EXPERIMENT, env=env, stdout=out, stderr=err,
            start_new_session=True, close_fds=True,
        )
        try:
            returncode = process.wait(timeout=OUTER_WALL_SECONDS)
        except subprocess.TimeoutExpired:
            timed_out = True
            os.killpg(process.pid, signal.SIGTERM)
            try:
                returncode = process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                returncode = process.wait()
    elapsed = time.monotonic() - started_monotonic
    after = input_hashes()
    record = {
        "schema": "pgfrag.native_cohesive_vessel_benchmark_execution.v1",
        "status": "PROCESS_TIMED_OUT_PRESERVED" if timed_out else "PROCESS_EXITED_PRESERVED",
        "attempt": "run_attempt_01",
        "outer_driver_argv": sys.argv,
        "native_argv": argv,
        "cwd": str(EXPERIMENT),
        "python": str(PYTHON),
        "thread_environment": {name: env[name] for name in (
            "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLAS_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "OMP_DYNAMIC",
        )},
        "started_unix_s": started_wall,
        "ended_unix_s": time.time(),
        "outer_elapsed_s": elapsed,
        "outer_wall_cap_s": OUTER_WALL_SECONDS,
        "timed_out": timed_out,
        "returncode": returncode,
        "live_processes_before_launch": active,
        "frozen_input_hashes_before": before,
        "frozen_input_hashes_after": after,
        "frozen_inputs_unchanged": before == after,
        "stdout": {"path": STDOUT.name, "bytes": STDOUT.stat().st_size, "sha256": sha256(STDOUT)},
        "stderr": {"path": STDERR.name, "bytes": STDERR.stat().st_size, "sha256": sha256(STDERR)},
        "run_output_exists": RUN.exists(),
        "no_retry_or_resume": True,
    }
    write_json(RECORD, record)
    print(json.dumps(record, allow_nan=False, separators=(",", ":")))
    return 0 if returncode == 0 and not timed_out and before == after else 1


if __name__ == "__main__":
    raise SystemExit(main())
