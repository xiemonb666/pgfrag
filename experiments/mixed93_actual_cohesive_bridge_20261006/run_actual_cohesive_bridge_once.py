"""Single-attempt supervisor: SOURCE_READY by default, one process only after root-go."""
from __future__ import annotations
import argparse, ast, hashlib, json, os, subprocess, sys, time
from pathlib import Path
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
REVIEW = ROOT / "pgfrag" / "tmp" / "review_20261006" / "mixed93_actual_cohesive_bridge_20261006"
PINS = REVIEW / "NUMERIC_SOURCE_PINS.json"
STATIC, WORKER = HERE / "mixed93_actual_cohesive_bridge.py", HERE / "run_actual_cohesive_bridge.py"
RUNS = REVIEW / "numeric_runs"
OUTPUT_LIMIT, WALL_LIMIT = 12 * 1024 * 1024, 60.0
ALLOWED = {"__future__", "argparse", "ast", "hashlib", "importlib", "json", "math", "os", "pathlib", "subprocess", "sys", "time"}
def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()
def fingerprints(pins: dict) -> dict:
    hashes = {}
    for name, ref in pins["references"].items():
        path = ROOT / Path(ref["path"])
        if not path.is_file(): raise FileNotFoundError(ref["path"])
        got = sha(path)
        if got != ref["sha256"]: raise ValueError(f"numeric source pin mismatch for {name}: {got}")
        hashes[name] = got
    cache = pins["accepted_cache_binding"]
    got = sha(ROOT / Path(cache["arrays_path"]))
    if got != cache["arrays_sha256"]: raise ValueError(f"accepted NPZ byte hash mismatch: {got}")
    hashes["accepted_arrays_npz"] = got
    return hashes
def validate_ast(path: Path, cap: int, worker: bool = False) -> None:
    source = path.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(path))
    parent = {child: node for node in ast.walk(tree) for child in ast.iter_child_nodes(node)}
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            roots = [alias.name.split(".", 1)[0] for alias in node.names]
        elif isinstance(node, ast.ImportFrom):
            roots = [(node.module or "").split(".", 1)[0]]
        else:
            continue
        for root in roots:
            cur, in_execute = node, False
            while cur in parent:
                cur = parent[cur]
                in_execute |= isinstance(cur, (ast.FunctionDef, ast.AsyncFunctionDef)) and cur.name == "execute"
            if root not in ALLOWED and not (worker and root == "numpy" and in_execute):
                raise ValueError(f"unapproved import {root} in {path.name}")
    if len(source.splitlines()) > cap: raise ValueError(f"{path.name} exceeds {cap} source lines")
def source_ready(pins: dict, require_disabled: bool = True) -> dict:
    if pins.get("schema") != "mixed93_actual_cohesive_bridge_numeric_source_pins.v1":
        raise ValueError("unexpected numeric source pin schema")
    if require_disabled and (pins.get("execution_authorized") is not False or pins.get("root_go", {}).get("authorized") is not False):
        raise ValueError("numeric pins must remain execution-disabled in SOURCE_READY")
    validate_ast(WORKER, 480, worker=True); validate_ast(Path(__file__).resolve(), 160)
    checker = subprocess.run([sys.executable, "-B", str(STATIC)], cwd=ROOT, capture_output=True, timeout=20)
    if checker.returncode != 0: raise ValueError(f"frozen static checker failed rc={checker.returncode}: {checker.stderr.decode(errors='replace')[-2000:]}")
    static_result = json.loads(checker.stdout.decode("utf-8"))
    if static_result.get("status") != "SOURCE_READY" or static_result.get("numeric_modules_imported") is not False: raise ValueError("frozen checker did not return SOURCE_READY")
    hashes = fingerprints(pins)
    return {"status": "SOURCE_READY", "scope": "AST and pinned metadata/hash checks only",
            "numeric_modules_imported": False, "npz_hashed": True, "npz_content_opened": False,
            "provider_called": False, "solver_called": False, "execution_authorized": False,
            "numeric_references_verified": len(hashes), "numeric_reference_hashes": hashes, "static_checker": static_result}
def serialize_receipt(path: Path, receipt: dict, fixed_bytes: int) -> tuple[int, int]:
    receipt["receipt_bytes"], receipt["combined_output_bytes"] = 0, 0
    for _ in range(32):
        blob = (json.dumps(receipt, sort_keys=True, indent=2, allow_nan=False) + "\n").encode("utf-8")
        size, total = len(blob), fixed_bytes + len(blob)
        if receipt["receipt_bytes"] == size and receipt["combined_output_bytes"] == total:
            path.write_bytes(blob)
            if path.stat().st_size != size: raise OSError("receipt byte count changed after write")
            return size, total
        receipt["receipt_bytes"], receipt["combined_output_bytes"] = size, total
    raise RuntimeError("receipt fixed-point byte accounting did not converge")
def main(argv: list[str] | None = None) -> int:
    all_started = time.monotonic()
    parser = argparse.ArgumentParser(description="one bounded actual-93 cohesive bridge attempt")
    parser.add_argument("--execute", action="store_true"); parser.add_argument("--root-go")
    args = parser.parse_args(argv)
    pins = json.loads(PINS.read_text(encoding="utf-8"))
    if not args.execute: print(json.dumps(source_ready(pins), sort_keys=True, allow_nan=False)); return 0
    go = pins.get("root_go")
    if (pins.get("execution_authorized") is not True or not isinstance(go, dict) or go.get("authorized") is not True
            or not isinstance(go.get("token"), str) or not go["token"] or args.root_go != go["token"]):
        print(json.dumps({"status": "EXECUTION_DISABLED", "numeric_modules_imported": False, "npz_opened": False})); return 2
    ready = source_ready(pins, require_disabled=False)
    pre = {"references": ready["numeric_reference_hashes"], "pins_json_sha256": sha(PINS)}
    RUNS.mkdir(parents=True, exist_ok=True)
    lock = RUNS / ".numeric_bridge.lock"
    fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    os.write(fd, json.dumps({"pid": os.getpid(), "started_unix": time.time()}).encode("utf-8")); os.close(fd)
    run = RUNS / "run01"; run.mkdir(exist_ok=False)
    result_path = run / "result.json"
    command = [sys.executable, "-B", str(WORKER), "--execute", "--root-go", args.root_go, "--output", str(result_path)]
    remaining = max(0.0, WALL_LIMIT - (time.monotonic() - all_started))
    timed_out, launch_error, outer_rc = False, None, None
    child_started = time.monotonic(); stdout = stderr = b""
    try:
        proc = subprocess.Popen(command, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        try:
            stdout, stderr = proc.communicate(timeout=remaining)
        except subprocess.TimeoutExpired:
            timed_out = True; proc.kill(); stdout, stderr = proc.communicate()
        outer_rc = proc.returncode
    except Exception as exc:
        launch_error = f"{type(exc).__name__}: {exc}"
    elapsed, total_elapsed = time.monotonic() - child_started, time.monotonic() - all_started
    (run / "stdout.bin").write_bytes(stdout); (run / "stderr.bin").write_bytes(stderr)
    result, parse_error = None, None
    if result_path.is_file():
        try:
            result = json.loads(result_path.read_text(encoding="utf-8"))
        except Exception as exc: parse_error = f"{type(exc).__name__}: {exc}"
    post = {"references": None, "pins_json_sha256": None, "error": None}
    try:
        post = {"references": fingerprints(pins), "pins_json_sha256": sha(PINS), "error": None}
    except Exception as exc:
        post["error"] = f"{type(exc).__name__}: {exc}"
    pins_equal = pre["references"] == post["references"] and pre["pins_json_sha256"] == post["pins_json_sha256"]
    result_bytes = result_path.stat().st_size if result_path.is_file() else 0
    log_bytes = len(stdout) + len(stderr); lock_bytes = lock.stat().st_size
    fixed_bytes = log_bytes + result_bytes + lock_bytes
    accepted = outer_rc == 0 and not timed_out and launch_error is None and total_elapsed <= WALL_LIMIT and pins_equal and post["error"] is None and result is not None and result.get("status") == "ACCEPTED_DIAGNOSTIC_ONLY"
    receipt = {"schema": "mixed93_actual_cohesive_bridge_outer_receipt.v1",
        "status": "ACCEPTED_DIAGNOSTIC_ONLY" if accepted else "NOT_ACCEPTED",
        "classification": "fixed-crack toy elastic cohesive bridge; no material/growth qualification",
        "run_directory": str(run), "attempt": 1, "lock": {"path": str(lock), "bytes": lock_bytes,
            "sha256": sha(lock), "acquired_exclusively": True, "retained_to_prevent_retry": True},
        "argv": command, "outer_return_code": outer_rc, "launch_error": launch_error, "timed_out": timed_out,
        "child_seconds": elapsed, "total_seconds_from_preflight": total_elapsed, "wall_limit_seconds": WALL_LIMIT,
        "source_pre_pins": pre, "source_post_pins": post, "source_pins_equal": pins_equal,
        "files": {name: ({"bytes": (run / name).stat().st_size, "sha256": sha(run / name)} if (run / name).is_file() else None)
                  for name in ("stdout.bin", "stderr.bin", "result.json")},
        "worker_result_parse_error": parse_error, "numeric_result_status": result.get("status") if result else None,
        "counts": result.get("counts", {}) if result else {}, "log_bytes": log_bytes, "result_bytes": result_bytes,
        "lock_bytes_included": True, "receipt_bytes_included": True, "combined_output_limit_bytes": OUTPUT_LIMIT,
        "budget_passed": None, "supervisor_failure": {"launch_error": launch_error, "parse_error": parse_error,
            "post_pin_error": post["error"], "source_pins_equal": pins_equal, "timed_out": timed_out}}
    receipt_path = run / "receipt.json"
    process_ok = accepted
    for _ in range(8):
        receipt_bytes, total_bytes = serialize_receipt(receipt_path, receipt, fixed_bytes)
        budget_passed = total_bytes <= OUTPUT_LIMIT
        status = "ACCEPTED_DIAGNOSTIC_ONLY" if process_ok and budget_passed else "NOT_ACCEPTED"
        if receipt.get("budget_passed") is budget_passed and receipt["status"] == status: break
        receipt["budget_passed"], receipt["status"] = budget_passed, status
    else:
        receipt["budget_passed"], receipt["status"] = False, "NOT_ACCEPTED"
        receipt["supervisor_failure"]["receipt_accounting"] = "serialized byte fixed point failed; fail closed"
        receipt_bytes, total_bytes = serialize_receipt(receipt_path, receipt, fixed_bytes)
    accepted = receipt["status"] == "ACCEPTED_DIAGNOSTIC_ONLY" and receipt["budget_passed"] is True
    print(json.dumps({"status": receipt["status"], "run_directory": str(run), "outer_return_code": outer_rc,
                      "combined_output_bytes": total_bytes, "receipt_bytes": receipt_bytes}, sort_keys=True, allow_nan=False))
    return 0 if accepted else 3


if __name__ == "__main__":
    raise SystemExit(main())
