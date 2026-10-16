"""Exclusive one-shot capture for the saved-history checker; execution requires root GO."""
from __future__ import annotations
import argparse, hashlib, json, subprocess, sys, time
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = Path(r"C:\Users\xiemo\Desktop\ShardForge")
CASE_PATH, PINS_PATH = HERE / "case.json", HERE / "SOURCE_PINS.json"; RUNTIME = ROOT / "third_party/vesselforge-runtime/venv/Scripts/python.exe"
CAP_BYTES, CHILD_SECONDS = 128 * 1024, 10.0
def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""): h.update(block)
    return h.hexdigest()

def pinned_inputs(case, pins):
    refs = {}
    def add(label, ref):
        refs[label] = (ROOT / ref["path"], ref["sha256"])
    for group in ("source_pins", "producer_support_files"):
        for name, ref in case[group].items(): add(group + ":" + name, ref)
    add("elastic_input", case["elastic_input"])
    for name, value in case["producer"].items():
        if name.endswith("_path"):
            key = name[:-5] + "_sha256"
            if key in case["producer"]: add("producer:" + name, {"path": value, "sha256": case["producer"][key]})
    for name, ref in pins["source_files"].items(): add("own:" + name, ref)
    b = case["pending_root_binding"]
    refs["opening_result"] = (ROOT / b["opening_result_path"], b["opening_result_sha256"])
    refs["opening_receipt"] = (ROOT / b["opening_receipt_path"], b["opening_receipt_sha256"])
    return refs

def snapshot(refs, enforce=True):
    out = {}
    for name, (path, expected) in refs.items():
        path = path.resolve(); path.relative_to(ROOT); actual = digest(path) if path.is_file() else None
        if enforce and actual != expected: raise RuntimeError("pinned input changed: " + name)
        out[name] = {"path": str(path), "exists": actual is not None, "sha256": actual}
    out["own:SOURCE_PINS"] = {"path": str(PINS_PATH.resolve()), "sha256": digest(PINS_PATH)}
    return out

def write_exclusive(path, raw):
    with Path(path).open("xb") as f: f.write(raw)

def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--execute", action="store_true"); p.add_argument("--root-go", required=True)
    p.add_argument("--result", type=Path, required=True); p.add_argument("--receipt", type=Path, required=True)
    a = p.parse_args(argv)
    try:
        if not a.execute: raise PermissionError("--execute is required")
        if str(ROOT) != str(HERE.parents[3]): raise RuntimeError("workspace root changed")
        if Path(sys.executable).resolve(strict=True) != RUNTIME.resolve(strict=True): raise RuntimeError("wrong bundled runtime")
        case = json.loads(CASE_PATH.read_text(encoding="utf-8")); pins = json.loads(PINS_PATH.read_text(encoding="utf-8"))
        go = pins.get("root_go", {})
        if pins.get("execution_authorized") is not True or go.get("authorized") is not True or not go.get("token") or a.root_go != go["token"]:
            raise PermissionError("exact independent root GO is required")
        contract = case["capture_contract"]; binding = case["pending_root_binding"]
        if digest(__file__) != contract["wrapper_sha256"] or digest(__file__) != pins["source_files"]["capture_wrapper"]["sha256"]:
            raise RuntimeError("capture wrapper pin mismatch")
        for name in ("checker", "case", "readme", "design"):
            ref = pins["source_files"][name]
            if digest(ROOT / ref["path"]) != ref["sha256"]: raise RuntimeError("own source pin mismatch: " + name)
        result, receipt = ((ROOT / binding["opening_result_path"]).resolve(strict=True),
                           (ROOT / binding["opening_receipt_path"]).resolve(strict=True))
        if a.result.resolve(strict=True) != result or a.receipt.resolve(strict=True) != receipt: raise RuntimeError("input paths differ from frozen binding")
        out = (ROOT / binding["independent_output_path"]).resolve(); out.relative_to(ROOT); run_dir = (ROOT / contract["run_directory"]).resolve()
        if out.parent != run_dir or out.name != "result.json" or run_dir.exists(): raise FileExistsError("configured exclusive run directory is unavailable")
        refs = pinned_inputs(case, pins); pre = snapshot(refs)
        if digest(result) != binding["opening_result_sha256"] or digest(receipt) != binding["opening_receipt_sha256"]: raise RuntimeError("saved opening input mismatch")
        run_dir.mkdir(exist_ok=False)
        stdout_path, stderr_path = run_dir / "stdout.bin", run_dir / "stderr.bin"
        token_hash = hashlib.sha256(a.root_go.encode("utf-8")).hexdigest()
        argv_child = [str(RUNTIME), "-B", str(HERE / "check_saved_opening_history.py"), "--execute", "--root-go", a.root_go,
                      "--result", str(result), "--receipt", str(receipt)]
        started = time.monotonic(); timed_out = False; rc = None; stdout = stderr = b""
        try:
            proc = subprocess.run(argv_child, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=CHILD_SECONDS, check=False)
            rc, stdout, stderr = proc.returncode, proc.stdout, proc.stderr
        except subprocess.TimeoutExpired as exc:
            timed_out, stdout, stderr = True, exc.stdout or b"", exc.stderr or b""
        child_seconds = time.monotonic() - started; write_exclusive(stdout_path, stdout); write_exclusive(stderr_path, stderr)
        post = snapshot(refs, enforce=False); same = pre == post
        out_info = {"exists":out.is_file(),"bytes":out.stat().st_size if out.is_file() else 0,"sha256":digest(out) if out.is_file() else None}
        accepted = rc == 0 and not timed_out and same and out_info["exists"]
        rec = {"schema":"mixed93_opening_history_independent_capture.v1","status":"CAPTURED" if accepted else "NOT_ACCEPTED",
            "actual_returncode":rc,"timed_out":timed_out,"child_started":True,"child_deadline_seconds":CHILD_SECONDS,
            "child_elapsed_seconds":child_seconds,"child_root_go_token_sha256":token_hash,"child_arguments_token_redacted":True,
            "source_hashes_pre":pre,"source_hashes_post":post,"pre_post_hashes_match":same,"result":out_info,
            "stdout":{"file":"stdout.bin","bytes":len(stdout),"sha256":hashlib.sha256(stdout).hexdigest()},
            "stderr":{"file":"stderr.bin","bytes":len(stderr),"sha256":hashlib.sha256(stderr).hexdigest()},
            "output_limit_bytes":CAP_BYTES,"checker_run":True,"retry_count":0}
        base = {"result.json":out_info["bytes"]} if out_info["exists"] else {}
        base.update({"stdout.bin":len(stdout),"stderr.bin":len(stderr)}); guess = 0
        for _ in range(20):
            tally = dict(base); tally["capture_receipt.json"] = guess
            total = sum(tally.values()); rec["file_tally_bytes"] = tally; rec["combined_bytes_including_receipt"] = total
            rec["combined_file_limit_pass"] = total <= CAP_BYTES
            if total > CAP_BYTES: rec["status"] = "NOT_ACCEPTED"
            raw = (json.dumps(rec,indent=2,sort_keys=True,allow_nan=False)+"\n").encode("utf-8")
            if len(raw) == guess: break
            guess = len(raw)
        else: raise RuntimeError("receipt accounting did not stabilize")
        write_exclusive(run_dir / "capture_receipt.json", raw)
        actual = {q.name:q.stat().st_size for q in run_dir.iterdir() if q.is_file()}
        if actual != rec["file_tally_bytes"] or sum(actual.values()) != rec["combined_bytes_including_receipt"]: raise RuntimeError("final output tally mismatch")
        print(json.dumps({"status":rec["status"],"actual_returncode":rc,"timed_out":timed_out,"output":str(out),"receipt":str(run_dir/"capture_receipt.json"),"bytes":sum(actual.values())},allow_nan=False))
        return 0 if rec["status"] == "CAPTURED" and rec["combined_file_limit_pass"] else 3
    except Exception as exc:
        print(json.dumps({"status":"NOT_ACCEPTED","failure_type":type(exc).__name__,"failure":str(exc),"retry_count":0},allow_nan=False))
        return 3
if __name__ == "__main__": raise SystemExit(main())
