"""Bounded scalar-only test of the driver's actual JSON default function."""
import datetime
import hashlib
import json
import pathlib
import subprocess
import sys
import time

ROOT = pathlib.Path.cwd().resolve()
REVIEW = ROOT / "pgfrag/tmp/review_20261006/mixed93_cached_lefm_solve_20261006"
OUTPUT = REVIEW / "serializer_only_check01"
ENTRY = ROOT / "pgfrag/experiments/mixed93_cached_lefm_solve_20261006/mixed93_cached_lefm_solve.py"
PINS = REVIEW / "SOURCE_PINS.json"
FREEZE = REVIEW / "failure_source_freeze"
ARTIFACTS = {
    "current_entry": ENTRY,
    "current_case": ENTRY.with_name("case.json"),
    "current_pins": PINS,
    "failed_entry": FREEZE / "mixed93_cached_lefm_solve.py",
    "failed_case": FREEZE / "case.authorized.json",
    "failed_pins": FREEZE / "SOURCE_PINS.authorized.json",
    "failure_manifest": FREEZE / "failure_source_freeze.json",
    "execution_receipt": REVIEW / "execution_capture/execution_receipt.json",
    "captured_stdout": REVIEW / "execution_capture/stdout.bin",
    "captured_stderr": REVIEW / "execution_capture/stderr.bin",
    "failed_npz": REVIEW / "run01/arrays.npz",
}
PROBE = r'''import ast, json, pathlib, numpy as np
source_path = pathlib.Path(r"SOURCE_PATH")
tree = ast.parse(source_path.read_text(encoding="utf-8"), filename=str(source_path))
writer = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "_write_result")
default_node = next(node for node in writer.body if isinstance(node, ast.FunctionDef) and node.name == "json_default")
namespace = {"np": np}
exec(compile(ast.Module(body=[default_node], type_ignores=[]), str(source_path), "exec"), namespace)
default = namespace["json_default"]
payload = {"bool": np.bool_(False), "int": np.int64(3), "float": np.float64(0.25)}
encoded = json.dumps(payload, default=default, allow_nan=False, sort_keys=True)
decoded = json.loads(encoded)
assert type(decoded["bool"]) is bool and decoded["bool"] is False
assert type(decoded["int"]) is int and decoded["int"] == 3
assert type(decoded["float"]) is float and decoded["float"] == 0.25
try:
    json.dumps({"nan": np.float64(float("nan"))}, default=default, allow_nan=False)
except ValueError:
    nan_rejected = True
else:
    raise AssertionError("NaN should still fail with allow_nan=False")
try:
    json.dumps({"unsupported": object()}, default=default, allow_nan=False)
except TypeError:
    unsupported_rejected = True
else:
    raise AssertionError("unsupported object should still raise TypeError")
print(json.dumps({"status":"PASS","decoded":decoded,"decoded_types":{k:type(v).__name__ for k,v in decoded.items()},
                  "nan_rejected":nan_rejected,"unsupported_rejected":unsupported_rejected},sort_keys=True))
'''


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def hashes():
    return {name: sha256(path) for name, path in ARTIFACTS.items()}


def main():
    if OUTPUT.exists():
        raise RuntimeError(f"refusing existing serialization check directory: {OUTPUT}")
    pins = json.loads(PINS.read_text(encoding="utf-8"))
    if pins["entry_sha256"] != sha256(ENTRY) or pins["failed_run_entry_sha256"] != sha256(ARTIFACTS["failed_entry"]):
        raise RuntimeError("current/failed entry SHA evidence mismatch")
    before = hashes()
    OUTPUT.mkdir(exist_ok=False)
    probe_path = OUTPUT / "probe_actual_json_default.py"
    probe_path.write_text(PROBE.replace("SOURCE_PATH", str(ENTRY)), encoding="utf-8")
    command = [sys.executable, "-B", str(probe_path)]
    started = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
    start = time.perf_counter()
    timed_out, returncode, stdout, stderr = False, None, b"", b""
    try:
        completed = subprocess.run(command, cwd=ROOT, capture_output=True, timeout=10, check=False)
        returncode, stdout, stderr = completed.returncode, completed.stdout, completed.stderr
    except subprocess.TimeoutExpired as exc:
        timed_out, stdout, stderr = True, exc.stdout or b"", exc.stderr or b""
    elapsed = time.perf_counter() - start
    (OUTPUT / "stdout.bin").write_bytes(stdout)
    (OUTPUT / "stderr.bin").write_bytes(stderr)
    after = hashes()
    if stdout:
        try:
            probe_result = json.loads(stdout.decode("utf-8"))
        except (UnicodeDecodeError, json.JSONDecodeError):
            probe_result = None
    else:
        probe_result = None
    receipt = {
        "schema": "mixed93_serializer_only_check.v1",
        "started_utc": started,
        "timeout_seconds": 10,
        "elapsed_seconds": elapsed,
        "timed_out": timed_out,
        "actual_process_return_code": returncode,
        "command": command,
        "probe_result": probe_result,
        "stdout_bytes": len(stdout),
        "stderr_bytes": len(stderr),
        "pre_hashes": before,
        "post_hashes": after,
        "artifacts_unchanged": before == after,
        "probe_sha256": sha256(probe_path),
    }
    (OUTPUT / "result.json").write_text(json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    size = sum(path.stat().st_size for path in OUTPUT.iterdir() if path.is_file())
    if size > 1024 * 1024:
        raise RuntimeError(f"serialization check artifacts exceed 1 MiB: {size}")
    print(json.dumps({"actual_rc": returncode, "timed_out": timed_out, "elapsed_seconds": elapsed,
                      "artifacts_unchanged": before == after, "result": str(OUTPUT / "result.json")}, sort_keys=True))


if __name__ == "__main__":
    main()
