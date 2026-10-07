"""One-shot, root-GO-gated external capture for the independent saved audit."""
from __future__ import annotations
import argparse, hashlib, json, subprocess, sys, time
from pathlib import Path

WORKSPACE = Path(r"C:\Users\xiemo\Desktop\ShardForge")
HERE = Path(__file__).resolve().parent
CAP_BYTES = 4 * 1024 * 1024

def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()

def snap(path: Path) -> dict:
    if not path.is_file():
        return {"exists": False, "bytes": None, "sha256": None}
    return {"exists": True, "bytes": path.stat().st_size, "sha256": sha(path)}

def main() -> int:
    parser = argparse.ArgumentParser(description="Capture exactly one authorized checker run; never calls a solver.")
    parser.add_argument("--root-go", required=True)
    args = parser.parse_args()
    if HERE != WORKSPACE / "pgfrag/tmp/review_20261007/mixed93_actual_cohesive_saved_independent":
        raise SystemExit("capture package is outside its pinned owner directory")
    case = json.loads((HERE / "case.json").read_text(encoding="utf-8"))
    pins = json.loads((HERE / "pins.json").read_text(encoding="utf-8"))
    gate = pins.get("root_go", {})
    if not (pins.get("execution_authorized") is True and gate.get("authorized") is True and args.root_go == gate.get("token")):
        raise SystemExit("capture disabled until the exact second root GO is pinned")
    run = HERE / "run02"
    if run.exists():
        raise SystemExit("run02 already exists; refusing another invocation")
    tracked = {name: WORKSPACE / ref["path"] for name, ref in case["source_pins"].items()}
    tracked.update({
        "independent_saved_checker": HERE / "check_saved_result.py",
        "independent_source_contract": HERE / "check_source_contract.py",
        "independent_case": HERE / "case.json",
        "independent_pins": HERE / "pins.json",
        "capture_source": HERE / "capture_once.py",
        "producer_result": WORKSPACE / pins["result"]["path"],
        "producer_receipt": WORKSPACE / pins["receipt"]["path"],
        "producer_npz": WORKSPACE / case["producer_npz"],
    })
    missing = [name for name, path in tracked.items() if not path.is_file()]
    if missing:
        raise SystemExit("prelaunch inventory missing pinned file(s): " + ", ".join(missing))
    expected = {name: ref["sha256"] for name, ref in case["source_pins"].items()}
    expected.update({
        "independent_saved_checker": pins["saved_check_sha256"],
        "independent_source_contract": pins["source_check_sha256"],
        "independent_case": pins["case_sha256"],
        "producer_result": pins["result"]["sha256"],
        "producer_receipt": pins["receipt"]["sha256"],
        "producer_npz": case["expected_hashes"]["producer_npz"],
    })
    wrong = [name for name, digest in expected.items() if sha(tracked[name]) != digest]
    if wrong:
        raise SystemExit("prelaunch pin mismatch: " + ", ".join(wrong))
    start = time.perf_counter()
    deadline = start + 15.0
    pre = {name: snap(path) for name, path in tracked.items()}
    output = run / "result.json"
    pre["independent_output"] = snap(output)
    run.mkdir()
    argv = [sys.executable, "-B", str(HERE / "check_saved_result.py"), "--execute", "--root-go", args.root_go, "--output", str(output)]
    code = None; timed_out = False; launch_error = None; stdout = b""; stderr = b""
    try:
        remaining = deadline - time.perf_counter()
        if remaining <= 0:
            raise TimeoutError("15-second whole-run deadline exhausted before checker launch")
        proc = subprocess.run(argv, cwd=str(HERE), stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=remaining, check=False, shell=False)
        code, stdout, stderr = proc.returncode, proc.stdout, proc.stderr
    except subprocess.TimeoutExpired as exc:
        timed_out = True; stdout = exc.stdout or b""; stderr = exc.stderr or b""
    except Exception as exc:
        launch_error = f"{type(exc).__name__}: {exc}"
    (run / "stdout.bin").write_bytes(stdout)
    (run / "stderr.bin").write_bytes(stderr)
    post = {name: snap(path) for name, path in tracked.items()}
    post["independent_output"] = snap(output)
    receipt = {
        "schema": "mixed93_actual_cohesive_saved_independent_external_receipt.v1",
        "argv": argv, "cwd": str(HERE), "return_code": code, "timed_out": timed_out,
        "launch_error": launch_error, "elapsed_seconds": round(time.perf_counter() - start, 6),
        "deadline_seconds": 15.0, "pre": pre, "post": post,
        "pre_post_equal": {name: pre[name] == post[name] for name in tracked},
        "stdout": {"path": "stdout.bin", "bytes": len(stdout), "sha256": hashlib.sha256(stdout).hexdigest()},
        "stderr": {"path": "stderr.bin", "bytes": len(stderr), "sha256": hashlib.sha256(stderr).hexdigest()},
        "tree_cap_bytes": CAP_BYTES,
    }
    receipt_path = run / "external_receipt.json"
    for _ in range(2):
        receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        total = sum(path.stat().st_size for path in run.rglob("*") if path.is_file())
        receipt["output_tree_bytes"] = total
        receipt["output_tree_within_cap"] = total <= CAP_BYTES
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    final_bytes = sum(path.stat().st_size for path in run.rglob("*") if path.is_file())
    print(json.dumps({"return_code": code, "timed_out": timed_out, "elapsed_seconds": receipt["elapsed_seconds"],
                      "output_tree_bytes": final_bytes, "output_tree_within_cap": final_bytes <= CAP_BYTES,
                      "stdout": receipt["stdout"], "stderr": receipt["stderr"], "result": post["independent_output"],
                      "receipt": str(receipt_path)}, sort_keys=True))
    return code if code is not None else (124 if timed_out else 125)

if __name__ == "__main__":
    raise SystemExit(main())
