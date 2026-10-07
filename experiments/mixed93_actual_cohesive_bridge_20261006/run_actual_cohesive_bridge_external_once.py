"""Independent 60-second Windows watchdog and raw capture for the pinned inner supervisor."""
from __future__ import annotations
import argparse, hashlib, json, os, subprocess, sys, time
from pathlib import Path
STARTED = time.monotonic()
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
REVIEW = ROOT / "pgfrag" / "tmp" / "review_20261006" / "mixed93_actual_cohesive_bridge_20261006"
PINS = REVIEW / "NUMERIC_SOURCE_PINS.json"
SUPERVISOR = HERE / "run_actual_cohesive_bridge_once.py"
RUNS = REVIEW / "numeric_runs"
CAPTURE = REVIEW / "ownedexternal_attempt_01"
PYTHON = Path(r"C:\Users\xiemo\Desktop\ShardForge\third_party\vesselforge-runtime\venv\Scripts\python.exe")
WALL, BYTE_LIMIT = 60.0, 12 * 1024 * 1024
def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""): h.update(block)
    return h.hexdigest()
def fingerprints(pins: dict) -> dict:
    found = {}
    for name, ref in pins["references"].items():
        path = (ROOT / ref["path"]).resolve()
        path.relative_to(ROOT)
        digest = sha(path)
        if digest != ref["sha256"]: raise ValueError(f"source hash mismatch: {name}={digest}")
        found[name] = digest
    cache = pins["accepted_cache_binding"]
    digest = sha((ROOT / cache["arrays_path"]).resolve())
    if digest != cache["arrays_sha256"]: raise ValueError(f"accepted cache raw-byte SHA mismatch: {digest}")
    found["accepted_arrays_npz_raw_sha256"] = digest
    return found
def authorize(pins: dict, token: str) -> None:
    go = pins.get("root_go")
    if (pins.get("schema") != "mixed93_actual_cohesive_bridge_numeric_source_pins.v1"
        or pins.get("execution_authorized") is not True or not isinstance(go, dict)
        or go.get("authorized") is not True or not isinstance(go.get("token"), str)
        or not go["token"] or token != go["token"]): raise PermissionError("exact root-go token is required")
    frozen = pins["frozen_contract"]
    if frozen["case_sha256"] != "fcc4b46ed101b0f6067a9f3ec10816909c9fa931edd5a8cdde4c7f46cb3f3a71": raise ValueError("frozen case pin changed")
    if frozen["original_source_pins_sha256"] != "9cb1fee478dadc080d22a872f29af80e91de8549902334be79fa518cac4554f3": raise ValueError("original static pins hash changed")
    if frozen["original_source_pins_execution_authorized"] is not False or frozen["case_execution_authorized"] is not False: raise ValueError("old case/source pins must remain disabled")
def tree_bytes(root: Path) -> int:
    return sum(p.stat().st_size for p in root.rglob("*") if p.is_file()) if root.exists() else 0
def receipt_write(path: Path, receipt: dict, fixed: int, candidate: bool) -> tuple[int, int, bool]:
    for _ in range(32):
        receipt["outer_elapsed_seconds"] = round(time.monotonic() - STARTED, 6)
        blob = (json.dumps(receipt, sort_keys=True, indent=2, allow_nan=False) + "\n").encode("utf-8")
        size, total = len(blob), fixed + len(blob)
        budget = total <= BYTE_LIMIT
        accepted = candidate and budget and time.monotonic() - STARTED <= WALL
        desired = {"receipt_bytes": size, "combined_output_bytes": total, "budget_passed": budget,
                   "status": "ACCEPTED_DIAGNOSTIC_ONLY" if accepted else "NOT_ACCEPTED", "wrapper_return_code": 0 if accepted else 3}
        if all(receipt.get(k) == v for k, v in desired.items()):
            path.write_bytes(blob)
            if path.stat().st_size != size: raise OSError("external receipt size changed after write")
            if time.monotonic() - STARTED > WALL: candidate = False; continue
            return size, total, accepted
        receipt.update(desired)
    receipt.update(status="NOT_ACCEPTED", wrapper_return_code=3, budget_passed=False, accounting_failure="receipt fixed point failed; fail closed")
    blob = (json.dumps(receipt, sort_keys=True, indent=2, allow_nan=False) + "\n").encode("utf-8")
    path.write_bytes(blob)
    return len(blob), fixed + len(blob), False
def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="one external wall-clock watchdog for the frozen bridge supervisor")
    parser.add_argument("--root-go", required=True); args = parser.parse_args(argv)
    try:
        pins = json.loads(PINS.read_text(encoding="utf-8")); authorize(pins, args.root_go)
    except Exception: return 2
    CAPTURE.mkdir(exist_ok=False)
    lock = CAPTURE / "external_once.lock"; fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    os.write(fd, b"one external attempt; retained to prevent retry\n"); os.close(fd)
    out, err = CAPTURE / "supervisor_stdout.bin", CAPTURE / "supervisor_stderr.bin"
    kill_out, kill_err = CAPTURE / "tree_kill_stdout.bin", CAPTURE / "tree_kill_stderr.bin"
    receipt_path = CAPTURE / "receipt.json"
    pre, post, pre_error, post_error = None, None, None, None
    child_rc, launch_error, timed_out, kill_rc = None, None, False, None
    child_seconds, command = 0.0, [str(PYTHON), "-B", str(SUPERVISOR), "--execute", "--root-go", args.root_go]
    pins = json.loads(PINS.read_text(encoding="utf-8"))
    try:
        if RUNS.exists(): raise FileExistsError("numeric_runs exists; this one-attempt capture refuses retry")
        if not PYTHON.is_file(): raise FileNotFoundError(str(PYTHON))
        pre = {"pins_json_sha256": sha(PINS), "references": fingerprints(pins), "external_wrapper_sha256": sha(Path(__file__).resolve())}
    except Exception as exc: pre_error = f"{type(exc).__name__}: {exc}"
    out.touch(exist_ok=False); err.touch(exist_ok=False)
    if pre_error is None:
        remaining = max(0.0, WALL - (time.monotonic() - STARTED)); began = time.monotonic()
        if remaining <= 0.0: timed_out = True
        else:
            try:
                with out.open("wb") as stdout, err.open("wb") as stderr:
                    proc = subprocess.Popen(command, cwd=ROOT, stdout=stdout, stderr=stderr, creationflags=getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0))
                    try: child_rc = proc.wait(timeout=remaining)
                    except subprocess.TimeoutExpired:
                        timed_out = True
                        try:
                            killed = subprocess.run(["taskkill.exe", "/PID", str(proc.pid), "/T", "/F"], capture_output=True, timeout=2.0, check=False)
                            kill_rc = killed.returncode; kill_out.write_bytes(killed.stdout); kill_err.write_bytes(killed.stderr)
                        except Exception as exc: launch_error = f"tree kill: {type(exc).__name__}: {exc}"
                        if proc.poll() is None: proc.kill()
                        try: proc.wait(timeout=2.0)
                        except subprocess.TimeoutExpired: launch_error = (launch_error or "") + "; child did not exit after kill"
                        child_rc = proc.returncode
            except Exception as exc: launch_error = f"{type(exc).__name__}: {exc}"
        child_seconds = time.monotonic() - began
    result, inner, parse_error = None, None, None
    for path, name in ((RUNS / "run01" / "result.json", "result"), (RUNS / "run01" / "receipt.json", "receipt")):
        if path.is_file():
            try:
                data = json.loads(path.read_text(encoding="utf-8"))
                if name == "result": result = data
                else: inner = data
            except Exception as exc: parse_error = f"{name}: {type(exc).__name__}: {exc}"
    try:
        pins_after = json.loads(PINS.read_text(encoding="utf-8"))
        post = {"pins_json_sha256": sha(PINS), "references": fingerprints(pins_after), "external_wrapper_sha256": sha(Path(__file__).resolve())}
    except Exception as exc: post_error = f"{type(exc).__name__}: {exc}"
    equal = pre is not None and post is not None and pre == post
    elapsed = time.monotonic() - STARTED
    inner_ok = inner is not None and inner.get("status") == "ACCEPTED_DIAGNOSTIC_ONLY" and inner.get("budget_passed") is True
    candidate = (pre_error is None and post_error is None and equal and not timed_out and launch_error is None
        and child_rc == 0 and elapsed <= WALL and result is not None and result.get("status") == "ACCEPTED_DIAGNOSTIC_ONLY"
        and inner_ok and (kill_rc is None))
    artifacts = {str(p.relative_to(ROOT)): {"bytes": p.stat().st_size, "sha256": sha(p)}
        for base in (RUNS, CAPTURE) if base.exists() for p in base.rglob("*") if p.is_file() and p != receipt_path}
    fixed = tree_bytes(RUNS) + tree_bytes(CAPTURE) - (receipt_path.stat().st_size if receipt_path.exists() else 0)
    receipt = {"schema": "mixed93_external_watchdog_receipt.v1", "status": "NOT_ACCEPTED",
        "classification": "external timing/capture only; no broader physical qualification",
        "command": command, "child_return_code": child_rc, "launch_error": launch_error, "timed_out": timed_out,
        "tree_kill_return_code": kill_rc, "child_seconds": child_seconds, "outer_elapsed_seconds": elapsed,
        "wall_limit_seconds": WALL, "source_pre": pre, "source_pre_error": pre_error,
        "source_post": post, "source_post_error": post_error, "source_pins_equal": equal,
        "worker_status": result.get("status") if result else None, "inner_supervisor_status": inner.get("status") if inner else None,
        "inner_parse_error": parse_error, "artifacts_before_external_receipt": artifacts,
        "lock_bytes": lock.stat().st_size, "raw_supervisor_stdout_stderr": True,
        "combined_output_limit_bytes": BYTE_LIMIT, "receipt_bytes": 0, "combined_output_bytes": 0,
        "budget_passed": False, "wrapper_return_code": 3}
    _, _, accepted = receipt_write(receipt_path, receipt, fixed, candidate)
    return 0 if accepted else 3
if __name__ == "__main__": raise SystemExit(main())
