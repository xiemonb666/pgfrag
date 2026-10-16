"""One-shot Windows capture wrapper; no numerical work is performed by this file."""
from __future__ import annotations
import hashlib, json, subprocess, time
from datetime import datetime, timezone
from pathlib import Path
WORKSPACE = Path(r"C:\Users\xiemo\Desktop\ShardForge")
ROOT = WORKSPACE.resolve(strict=True)
HERE = ROOT / "pgfrag" / "tmp" / "review_20261007" / "mixed93_opening_control_20261007"
CASE_PATH, PINS_PATH = HERE / "capture_CASE.json", HERE / "capture_PINS.json"
DEADLINE_S, OUTPUT_CAP = 20.0, 8 * 1024 * 1024
def stamp(): return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def rooted(relative): p = (ROOT / relative).resolve(strict=True); p.relative_to(ROOT); return p
def audit(case, pins):
    if ROOT != WORKSPACE.resolve(strict=True) or case["workspace_root_literal"] != str(WORKSPACE):
        raise RuntimeError("literal workspace root changed")
    frozen = pins["frozen_contract"]
    if digest(CASE_PATH) != frozen["capture_case_sha256"] or digest(Path(__file__).resolve(strict=True)) != frozen["capture_driver_sha256"]:
        raise RuntimeError("capture case or wrapper hash mismatch")
    source_case, driver = rooted(case["opening_control_case"]), rooted(case["opening_control_driver"])
    source_pin_path = rooted(case["opening_control_source_pins"])
    active = json.loads(source_pin_path.read_text(encoding="utf-8"))
    refs = active.get("references", {})
    if len(refs) != case["expected_source_reference_count"] or refs != pins["source_references"]:
        raise RuntimeError("opening source pin reference set changed")
    if active.get("execution_authorized") is not True or active.get("root_go", {}).get("authorized") is not True:
        raise PermissionError("opening source pins are not root-go authorized")
    token = active["root_go"].get("token")
    capture_go = pins.get("root_go", {})
    if (pins.get("execution_authorized") is not True or capture_go.get("authorized") is not True
            or not isinstance(token, str) or not token or capture_go.get("token") != token):
        raise PermissionError("capture pins do not carry the matching exact root-go token")
    if source_case != rooted(active["frozen_contract"]["case_path"]):
        raise RuntimeError("opening case path differs from source pins")
    if driver != rooted(active["frozen_contract"]["driver_path"]):
        raise RuntimeError("opening driver path differs from source pins")
    snapshot = rooted(case["disabled_source_pins_snapshot"])
    if digest(snapshot) != frozen["disabled_source_pins_snapshot_sha256"]:
        raise RuntimeError("preserved disabled source-pin snapshot changed")
    original_driver = rooted(case["pre_fix_driver_snapshot"])
    if digest(original_driver) != frozen["pre_fix_driver_snapshot_sha256"]:
        raise RuntimeError("preserved pre-fix driver bytes changed")
    original = json.loads(snapshot.read_text(encoding="utf-8"))
    if original.get("execution_authorized") is not False or original.get("root_go", {}).get("authorized") is not False:
        raise RuntimeError("saved source-pin bytes are not the original disabled version")
    paths = {"capture_case": CASE_PATH, "capture_driver": Path(__file__).resolve(strict=True),
             "capture_pins": PINS_PATH, "opening_case": source_case, "opening_driver": driver,
             "opening_source_pins": source_pin_path, "disabled_snapshot": snapshot,
             "pre_fix_driver_snapshot": original_driver}
    for name, item in refs.items():
        paths["reference:" + name] = rooted(item["path"])
    for name in case["trace_dependency_names"]:
        paths["trace_dependency:" + name] = rooted(refs[name]["path"])
    runtime = Path(case["python_executable_literal"]).resolve(strict=True)
    if str(runtime) != case["python_executable_literal"]:
        raise RuntimeError("pinned VF Python runtime path changed")
    return token, paths, runtime
def hashes(paths): return {name: {"path": str(p), "sha256": digest(p) if p.is_file() else None, "exists": p.is_file()} for name, p in paths.items()}
def write_receipt(path, receipt, sizes):
    receipt["file_tally_bytes"] = {**sizes, path.name: 0}
    for _ in range(8):
        raw = (json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n").encode("utf-8")
        total = sum(sizes.values()) + len(raw)
        receipt["file_tally_bytes"][path.name] = len(raw)
        receipt["combined_bytes_including_receipt"] = total
        receipt["combined_file_limit_bytes"] = OUTPUT_CAP
        receipt["combined_file_limit_pass"] = total <= OUTPUT_CAP
        if total > OUTPUT_CAP:
            receipt["status"] = "NOT_ACCEPTED"
            receipt["failure_reason"] = "combined run01 files exceeded 8 MiB including receipt"
        updated = (json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n").encode("utf-8")
        if len(updated) == len(raw) and sum(sizes.values()) + len(updated) == total:
            raw = updated
            break
    else:
        raise RuntimeError("receipt size accounting did not converge")
    with path.open("xb") as f:
        f.write(raw)
    return sum(sizes.values()) + path.stat().st_size
def main():
    started = time.monotonic()
    deadline = started + DEADLINE_S
    pre_start, pre_start_utc = time.monotonic(), stamp(); case, pins = json.loads(CASE_PATH.read_text(encoding="utf-8")), json.loads(PINS_PATH.read_text(encoding="utf-8"))
    token, paths, runtime = audit(case, pins)
    pre = hashes(paths)
    pre_end, pre_end_utc = time.monotonic(), stamp()
    if pre_end >= deadline:
        raise TimeoutError("20-second deadline expired during preflight hashing")
    run_dir = (ROOT / case["run_directory"]).resolve()
    run_dir.relative_to(ROOT)
    run_dir.mkdir(exist_ok=False)
    out, err = run_dir / case["stdout_name"], run_dir / case["stderr_name"]
    result, receipt_path = run_dir / case["result_name"], run_dir / case["receipt_name"]
    child_argv = [str(runtime), "-B", str(paths["opening_control_driver"]), "--execute", "--root-go", token, "--output", str(result)]
    child_start, child_start_utc = time.monotonic(), stamp()
    child_end, child_end_utc, proc, timed_out, tree = None, None, None, False, "not_needed"
    with out.open("xb") as stdout, err.open("xb") as stderr:
        token_forwarded = False
        try:
            proc = subprocess.Popen(child_argv, cwd=str(ROOT), stdin=subprocess.DEVNULL,
                stdout=stdout, stderr=stderr, creationflags=subprocess.CREATE_NEW_PROCESS_GROUP)
            token_forwarded = True
            try:
                proc.wait(timeout=max(0.0, deadline - time.monotonic()))
            except subprocess.TimeoutExpired:
                timed_out = True
                try:
                    killed = subprocess.run(["taskkill.exe", "/PID", str(proc.pid), "/T", "/F"],
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=1.0, check=False)
                    tree = "taskkill_tree_returncode_" + str(killed.returncode)
                except Exception as exc:
                    tree = "taskkill_failed_" + type(exc).__name__
                if proc.poll() is None:
                    proc.kill()
                try:
                    proc.wait(timeout=1.0)
                except subprocess.TimeoutExpired:
                    tree += "_direct_child_wait_timeout"
        except Exception as exc:
            tree = "spawn_failed_" + type(exc).__name__
        child_end, child_end_utc = time.monotonic(), stamp()
    post_start, post_start_utc = time.monotonic(), stamp()
    post = hashes(paths)
    post_hash_end, post_hash_end_utc = time.monotonic(), stamp()
    same = pre == post
    rc = None if proc is None else proc.returncode
    result_exists = result.is_file()
    result_info = {"exists": result_exists, "bytes": result.stat().st_size if result_exists else 0,
                   "sha256": digest(result) if result_exists else None}
    sizes = {p.name: p.stat().st_size for p in (out, err)}
    if result_exists:
        sizes[result.name] = result.stat().st_size
    post_end, post_end_utc = time.monotonic(), stamp()
    receipt = {"schema": "mixed93_opening_control_capture_receipt.v1", "status": "CAPTURED",
        "actual_returncode": rc, "timed_out": timed_out, "timeout_tree_termination": tree,
        "child_argv_root_go_forwarded_exactly": token_forwarded, "child_started": proc is not None,
        "child_stdout_raw_file": out.name, "child_stderr_raw_file": err.name,
        "result": result_info, "source_hashes_pre": pre, "source_hashes_post": post,
        "pre_post_hashes_match": same, "checker_run": False,
        "times_utc": {"preflight_started": pre_start_utc, "preflight_finished": pre_end_utc,
            "child_started": child_start_utc, "child_finished": child_end_utc,
            "post_accounting_started": post_start_utc, "post_hashing_finished": post_hash_end_utc,
            "post_accounting_finished": post_end_utc},
        "elapsed_seconds": {"preflight_including_config_and_hashes": pre_end-started,
            "run_directory_reservation": child_start-pre_end, "child_wall": child_end-child_start,
            "post_hashing": post_hash_end-post_start, "post_accounting_including_result_digest_and_sizes": post_end-post_start,
            "through_post_accounting": post_end-started,
            "deadline_seconds": DEADLINE_S}}
    if rc != 0 or timed_out or not same or post_end >= deadline or not result_exists:
        receipt["status"] = "NOT_ACCEPTED"
    total = write_receipt(receipt_path, receipt, sizes)
    if total != receipt["combined_bytes_including_receipt"]:
        raise RuntimeError("receipt byte tally mismatch")
    if time.monotonic() >= deadline:
        return 3
    return 0 if receipt["status"] == "CAPTURED" and total <= OUTPUT_CAP else 3
if __name__ == "__main__":
    raise SystemExit(main())
