"""One-shot Windows capture wrapper; no numerical work is performed here."""
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

def rooted(relative):
    path = (ROOT / relative).resolve(strict=True); path.relative_to(ROOT)
    if not path.is_file(): raise RuntimeError("pinned source is not a file: " + relative)
    return path

def audit(case, pins):
    if ROOT != WORKSPACE.resolve(strict=True) or case["workspace_root_literal"] != str(WORKSPACE): raise RuntimeError("literal workspace root changed")
    frozen = pins["frozen_contract"]
    source_case, driver = rooted(case["opening_control_case"]), rooted(case["opening_control_driver"])
    source_pin_path = rooted(case["opening_control_source_pins"])
    active = json.loads(source_pin_path.read_text(encoding="utf-8"))
    refs = active.get("references", {})
    if len(refs) != case["expected_source_reference_count"] or refs != pins["source_references"]: raise RuntimeError("opening source pin reference set changed")
    if source_case != rooted(active["frozen_contract"]["case_path"]): raise RuntimeError("opening case path differs from source pins")
    if driver != rooted(active["frozen_contract"]["driver_path"]): raise RuntimeError("opening driver path differs from source pins")
    snapshot = rooted(case["disabled_source_pins_snapshot"])
    final_snapshot = rooted(case["final_disabled_source_pins_snapshot"])
    old_driver = rooted(case["pre_fix_driver_snapshot"])
    prior_case = rooted(case["pre_vf_capture_case_snapshot"])
    prior_wrapper = rooted(case["pre_vf_capture_driver_snapshot"])
    if case["pre_vf_capture_pins_existed"] is not False or pins["pre_vf_capture_pins_existed"] is not False: raise RuntimeError("pre-VF capture-pin absence record changed")
    for saved in (snapshot, final_snapshot):
        value = json.loads(saved.read_text(encoding="utf-8"))
        if value.get("execution_authorized") is not False or value.get("root_go", {}).get("authorized") is not False: raise RuntimeError("saved source-pin bytes are not disabled")
    if active.get("execution_authorized") is not True or active.get("root_go", {}).get("authorized") is not True: raise PermissionError("opening source pins are not root-go authorized")
    token = active["root_go"].get("token")
    capture_go = pins.get("root_go", {})
    if (pins.get("execution_authorized") is not True or capture_go.get("authorized") is not True
            or not isinstance(token, str) or not token or capture_go.get("token") != token):
        raise PermissionError("capture pins do not carry the matching exact root-go token")
    if active.get("frozen_contract", {}).get("execution_disabled_when_frozen") is not True: raise RuntimeError("opening source pins frozen contract changed")
    runtime = Path(case["python_executable_literal"]).resolve(strict=True)
    if str(runtime) != case["python_executable_literal"]:
        raise RuntimeError("pinned VF Python runtime path changed")
    paths = {"capture_case": CASE_PATH, "capture_driver": Path(__file__).resolve(strict=True),
        "capture_pins": PINS_PATH, "opening_case": source_case, "opening_driver": driver,
        "opening_source_pins": source_pin_path, "disabled_snapshot": snapshot,
        "final_disabled_snapshot": final_snapshot, "pre_fix_driver_snapshot": old_driver,
        "pre_vf_capture_case_snapshot": prior_case, "pre_vf_capture_driver_snapshot": prior_wrapper}
    for name, item in refs.items():
        paths["reference:" + name] = rooted(item["path"])
    for name in case["trace_dependency_names"]:
        paths["trace_dependency:" + name] = rooted(refs[name]["path"])
    return token, paths, runtime, active, refs

def hashes(paths): return {name: {"path": str(path), "sha256": digest(path) if path.is_file() else None, "exists": path.is_file()} for name, path in paths.items()}

def verify_pre(case, pins, active, refs, pre):
    frozen, source = pins["frozen_contract"], active["frozen_contract"]
    expected = {"capture_case": frozen["capture_case_sha256"],
        "capture_driver": frozen["capture_driver_sha256"],
        "opening_case": source["case_sha256"], "opening_driver": source["driver_sha256"],
        "opening_source_pins": frozen["opening_source_pins_sha256"],
        "disabled_snapshot": frozen["disabled_source_pins_snapshot_sha256"],
        "final_disabled_snapshot": frozen["final_disabled_source_pins_snapshot_sha256"],
        "pre_fix_driver_snapshot": frozen["pre_fix_driver_snapshot_sha256"],
        "pre_vf_capture_case_snapshot": frozen["pre_vf_capture_case_snapshot_sha256"],
        "pre_vf_capture_driver_snapshot": frozen["pre_vf_capture_driver_snapshot_sha256"]}
    if source["case_sha256"] != frozen["opening_case_sha256"] or source["driver_sha256"] != frozen["opening_driver_sha256"]: raise RuntimeError("opening case/driver hash contract changed")
    for name, item in refs.items():
        expected["reference:" + name] = item["sha256"]
    for name, wanted in expected.items():
        if not pre[name]["exists"] or pre[name]["sha256"] != wanted:
            raise RuntimeError("pre-child SHA-256 mismatch: " + name)

def file_sizes(directory): return {path.relative_to(directory).as_posix(): path.stat().st_size for path in directory.rglob("*") if path.is_file()}

def write_receipt(path, receipt, sizes):
    receipt["file_tally_bytes"] = {**sizes, path.name: 0}
    for _ in range(12):
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
    with path.open("xb") as stream:
        stream.write(raw)
    return sum(sizes.values()) + path.stat().st_size

def main():
    started = time.monotonic()
    deadline = started + DEADLINE_S
    pre_start, pre_start_utc = time.monotonic(), stamp()
    case, pins = json.loads(CASE_PATH.read_text(encoding="utf-8")), json.loads(PINS_PATH.read_text(encoding="utf-8"))
    token, paths, runtime, active, refs = audit(case, pins)
    pre = hashes(paths)
    verify_pre(case, pins, active, refs, pre)
    pre_end, pre_end_utc = time.monotonic(), stamp()
    if pre_end >= deadline:
        raise TimeoutError("20-second deadline expired during preflight hashing")
    run_dir = (ROOT / case["run_directory"]).resolve()
    run_dir.relative_to(ROOT)
    run_dir.mkdir(exist_ok=False)
    out, err = run_dir / case["stdout_name"], run_dir / case["stderr_name"]
    result, receipt_path = run_dir / case["result_name"], run_dir / case["receipt_name"]
    argv = [str(runtime), "-B", str(paths["opening_driver"]), "--execute", "--root-go", token,
        "--output", str(result)]
    window_start, window_start_utc = time.monotonic(), stamp()
    process_start_utc, proc = None, None
    timed_out, token_forwarded, tree = False, False, "not_needed"
    with out.open("xb") as stdout, err.open("xb") as stderr:
        try:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                timed_out, tree = True, "deadline_expired_before_spawn"
            else:
                proc = subprocess.Popen(argv, cwd=str(ROOT), stdin=subprocess.DEVNULL,
                    stdout=stdout, stderr=stderr, creationflags=subprocess.CREATE_NEW_PROCESS_GROUP)
                process_start_utc, token_forwarded = stamp(), True
                try:
                    proc.wait(timeout=max(0.0, deadline - time.monotonic()))
                except subprocess.TimeoutExpired:
                    timed_out = True
                    try:
                        kill_wait = min(1.0, max(0.05, deadline - time.monotonic()))
                        killed = subprocess.run(["taskkill.exe", "/PID", str(proc.pid), "/T", "/F"],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=kill_wait, check=False)
                        tree = "taskkill_tree_returncode_" + str(killed.returncode)
                    except Exception as exc:
                        tree = "taskkill_failed_" + type(exc).__name__
                    if proc.poll() is None:
                        proc.kill()
                    try:
                        proc.wait(timeout=0.25)
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
    sizes = file_sizes(run_dir)
    post_end, post_end_utc = time.monotonic(), stamp()
    failure = [message for condition, message in ((rc != 0, "child return code was not zero"),
        (timed_out, "20-second deadline expired"), (not same, "pinned source hashes changed during child run"),
        (post_end >= deadline, "deadline expired during post-accounting"), (not result_exists, "result.json was not created")) if condition]
    receipt = {"schema": "mixed93_opening_control_capture_receipt.v1",
        "status": "NOT_ACCEPTED" if failure else "CAPTURED", "failure_reasons": failure,
        "actual_returncode": rc, "timed_out": timed_out, "timeout_tree_termination": tree,
        "child_root_go_token_forwarded_exactly": token_forwarded,
        "child_root_go_token_sha256": hashlib.sha256(token.encode("utf-8")).hexdigest(),
        "child_started": proc is not None, "child_stdout_raw_file": out.name,
        "child_stderr_raw_file": err.name, "result": result_info,
        "source_hashes_pre": pre, "source_hashes_post": post,
        "pre_post_hashes_match": same, "checker_run": False,
        "times_utc": {"preflight_started": pre_start_utc, "preflight_finished": pre_end_utc,
            "child_window_started": window_start_utc, "child_process_started": process_start_utc,
            "child_finished": child_end_utc, "post_accounting_started": post_start_utc,
            "post_hashing_finished": post_hash_end_utc, "post_accounting_finished": post_end_utc},
        "elapsed_seconds": {"preflight_including_config_and_hashes": pre_end-started,
            "run_directory_reservation": window_start-pre_end, "child_wall": child_end-window_start,
            "post_hashing": post_hash_end-post_start,
            "post_accounting_including_result_digest_and_sizes": post_end-post_start,
            "through_post_accounting": post_end-started,
            "remaining_before_receipt_seconds": max(0.0, deadline-post_end),
            "deadline_seconds": DEADLINE_S}}
    receipt_start = time.monotonic()
    total = write_receipt(receipt_path, receipt, sizes)
    actual_sizes = file_sizes(run_dir)
    actual_total = sum(actual_sizes.values())
    if total != actual_total or actual_sizes != receipt["file_tally_bytes"]:
        return 4
    if time.monotonic() >= deadline:
        return 3
    return 0 if receipt["status"] == "CAPTURED" and actual_total <= OUTPUT_CAP else 3

if __name__ == "__main__":
    raise SystemExit(main())
