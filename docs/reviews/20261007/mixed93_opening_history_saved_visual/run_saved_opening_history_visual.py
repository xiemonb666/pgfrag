"""Standard-library source gate and one-shot, bounded saved-result renderer."""
import argparse
import ast
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
CASE, PINS = HERE / "case.json", HERE / "SOURCE_PINS.json"
RENDERER = HERE / "render_saved_opening_history.py"
MAX_OUTPUT = 2 * 1024 * 1024

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def resolve(path):
    value = Path(path)
    return value if value.is_absolute() else ROOT / value

def check_sources(case, pins):
    if pins.get("schema") != "mixed93_opening_history_saved_visual_pins.v1":
        raise RuntimeError("unexpected saved-visual pin schema")
    local = (("case", CASE), ("renderer", RENDERER), ("runner", Path(__file__).resolve()), ("readme", HERE / "README.md"))
    for name, path in local:
        source = path.read_text(encoding="utf-8")
        if len(source.splitlines()) > (180 if name == "renderer" else 140):
            raise RuntimeError(f"{name} exceeds the frozen source-line limit")
        if pins["references"][name].get("sha256") != sha(path):
            raise RuntimeError(f"source SHA-256 mismatch: {name}")
        if path.suffix == ".py":
            for node in ast.parse(source, filename=str(path)).body:
                if isinstance(node, ast.Import) and any(a.name.split(".")[0] in {"numpy", "matplotlib"} for a in node.names):
                    raise RuntimeError(f"numeric/plot import at module load: {name}")
                if isinstance(node, ast.ImportFrom) and (node.module or "").split(".")[0] in {"numpy", "matplotlib"}:
                    raise RuntimeError(f"numeric/plot import at module load: {name}")
    producer = case["producer_sources"]
    for key, ref_name in (("driver", "producer_driver"), ("case", "producer_case")):
        item = producer[key]
        ref = pins["references"][ref_name]
        path = resolve(item["path"])
        if ref.get("path") != item["path"] or ref.get("sha256") != item["sha256"] or sha(path) != item["sha256"]:
            raise RuntimeError(f"producer {key} source pin mismatch")
    planned = case["input_result"]
    if planned.get("path") != pins["input_result"].get("path") or planned.get("sha256") != pins["input_result"].get("sha256"):
        raise RuntimeError("case and pin result identities differ")
    for name, item in pins["references"].items():
        path = resolve(item["path"])
        if not path.is_file() or sha(path) != item["sha256"]:
            raise RuntimeError(f"pinned source changed: {name}")

def ready():
    case, pins = json.loads(CASE.read_text(encoding="utf-8")), json.loads(PINS.read_text(encoding="utf-8"))
    check_sources(case, pins)
    return case, pins

def pinned_paths(case, pins, result_path):
    paths = {"case": CASE, "pins": PINS, "renderer": RENDERER, "runner": Path(__file__).resolve(), "result": result_path}
    paths.update({"source:" + name: resolve(item["path"]) for name, item in pins["references"].items()})
    return paths

def digest_paths(paths):
    return {name: sha(path) for name, path in paths.items()}


def encode_capture(capture, files):
    fixed = sum(x["bytes"] for x in files)
    for _ in range(8):
        blob = (json.dumps(capture, indent=2, sort_keys=True, allow_nan=False) + "\n").encode("utf-8")
        total = fixed + len(blob)
        if capture.get("combined_output_bytes") == total: return blob, total
        capture["combined_output_bytes"] = total
    raise RuntimeError("capture byte accounting did not stabilize")

def run(case, pins, token, task_started):
    go, auth, ref = pins.get("root_go", {}), case.get("execution_authorization", {}), case["input_result"]
    if (auth.get("authorized") is not True or auth.get("authorized_by") != "ROOT" or pins.get("execution_authorized") is not True
            or go.get("authorized") is not True or not token or token != go.get("token")):
        raise RuntimeError("BLOCKED: exact root render GO is required")
    if not ref.get("sha256"):
        raise RuntimeError("BLOCKED: actual opening-history result SHA-256 is pending")
    result_path = resolve(ref["path"])
    if not result_path.is_file() or sha(result_path) != ref["sha256"]:
        raise RuntimeError("BLOCKED: result is absent or differs from its pin")
    paths, out = pinned_paths(case, pins, result_path), resolve(case["limits"]["planned_output_dir"])
    out.mkdir(parents=True, exist_ok=False)
    png = out / "opening_history_saved_diagnostic.png"
    cmd = [sys.executable, "-B", str(RENDERER), "--case", str(CASE), "--pins", str(PINS), "--output", str(png), "--root-go", token]
    before = digest_paths(paths)
    started, timed_out, rc, stdout, stderr = time.perf_counter(), False, None, b"", b""
    try:
        remaining = max(0.0, case["limits"]["max_wall_seconds"] - (started - task_started))
        proc = subprocess.run(cmd, cwd=ROOT, capture_output=True, timeout=remaining)
        rc, stdout, stderr = proc.returncode, proc.stdout, proc.stderr
    except subprocess.TimeoutExpired as exc:
        timed_out, stdout, stderr = True, exc.stdout or b"", exc.stderr or b""
    elapsed = time.perf_counter() - started
    (out / "stdout.bin").write_bytes(stdout); (out / "stderr.bin").write_bytes(stderr)
    after = digest_paths(paths)
    files = [{"path": p.name, "bytes": p.stat().st_size, "sha256": sha(p)}
             for p in sorted(out.iterdir()) if p.is_file() and p.name != "capture.json"]
    capture = {"schema": "mixed93_opening_history_saved_visual_capture.v1",
        "status": "RENDERED" if rc == 0 and not timed_out and before == after else "NOT_ACCEPTED",
        "return_code": rc, "timed_out": timed_out, "wall_seconds": time.perf_counter() - task_started,
        "subprocess_seconds": elapsed, "input_path": ref["path"], "input_sha256": ref["sha256"],
        "pins_pre": before, "pins_post": after, "outputs_excluding_capture": files,
        "byte_limit": MAX_OUTPUT, "root_go_authorized_by": auth.get("authorized_by"), "no_retry": True}
    blob, total = encode_capture(capture, files)
    if total > MAX_OUTPUT:
        capture["status"] = "BLOCKED_OUTPUT_BYTE_LIMIT"
        blob, total = encode_capture(capture, files)
    (out / "capture.json").write_bytes(blob)
    actual_total = sum(p.stat().st_size for p in out.iterdir() if p.is_file())
    if actual_total != total or time.perf_counter() - task_started > case["limits"]["max_wall_seconds"] or actual_total > MAX_OUTPUT:
        raise RuntimeError("render capture exceeded its fixed budget; no retry permitted")
    print(json.dumps({"status": capture["status"], "result_sha256": ref["sha256"], "combined_output_bytes": actual_total}, sort_keys=True))
    return 0 if capture["status"] == "RENDERED" else 1

def main():
    task_started = time.perf_counter()
    parser = argparse.ArgumentParser(); parser.add_argument("--render", action="store_true"); parser.add_argument("--root-go")
    args = parser.parse_args()
    case, pins = ready()
    if not args.render:
        state = "SOURCE_READY_RESULT_PENDING" if case["input_result"]["sha256"] is None else "SOURCE_READY_RESULT_PINNED_NO_RENDER_REQUESTED"
        print(json.dumps({"status": state, "result_sha256": case["input_result"]["sha256"], "matplotlib_imported": False, "rendered": False}, sort_keys=True))
        return 0
    return run(case, pins, args.root_go, task_started)


if __name__ == "__main__":
    raise SystemExit(main())
