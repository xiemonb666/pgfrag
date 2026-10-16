"""SOURCE_READY metadata gate and one-shot, bounded saved-result renderer."""
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
CASE, PINS, RENDERER = HERE / "case.json", HERE / "SOURCE_PINS.json", HERE / "render_saved_bridge.py"
MAX_LINES, MAX_OUTPUT = 130, 2 * 1024 * 1024


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def resolve(path):
    value = Path(path)
    return value if value.is_absolute() else ROOT / value

def check_sources(case, pins):
    if pins.get("schema") != "mixed93_actual_cohesive_saved_visual_pins.v1":
        raise RuntimeError("unexpected saved-visual pin schema")
    for name, path in (("case", CASE), ("runner", Path(__file__).resolve()), ("renderer", RENDERER), ("readme", HERE / "README.md")):
        record = pins.get("references", {}).get(name, {})
        source = path.read_text(encoding="utf-8")
        if len(source.splitlines()) > (160 if name == "renderer" else MAX_LINES):
            raise RuntimeError(f"{name} exceeds its frozen source-line limit")
        if record.get("sha256") != sha(path):
            raise RuntimeError(f"source SHA-256 mismatch: {name}")
        if path.suffix == ".py":
            tree = ast.parse(source, filename=str(path))
            for node in tree.body:
                if isinstance(node, ast.Import) and any(a.name.split(".")[0] in {"numpy", "matplotlib"} for a in node.names):
                    raise RuntimeError(f"numeric/plot import must not occur at module load: {name}")
                if isinstance(node, ast.ImportFrom) and (node.module or "").split(".")[0] in {"numpy", "matplotlib"}:
                    raise RuntimeError(f"numeric/plot import must not occur at module load: {name}")
    if len(Path(__file__).read_text(encoding="utf-8").splitlines()) > MAX_LINES:
        raise RuntimeError("runner exceeds source-line limit")
    if case.get("schema") != "mixed93_actual_cohesive_saved_visual_case.v1":
        raise RuntimeError("unexpected case schema")


def ready():
    case, pins = json.loads(CASE.read_text(encoding="utf-8")), json.loads(PINS.read_text(encoding="utf-8"))
    check_sources(case, pins)
    worker = case["producer_worker_source"]
    if not (pins.get("producer_worker_source") == worker and resolve(worker["path"]).is_file() and sha(resolve(worker["path"])) == worker["sha256"]): raise RuntimeError("numeric worker source SHA-256 mismatch")
    planned = case["input_result"]
    if planned.get("path") != pins.get("input_result", {}).get("path"):
        raise RuntimeError("case and pin result paths differ")
    if planned.get("sha256") != pins.get("input_result", {}).get("sha256"):
        raise RuntimeError("case and pin result hashes differ")
    return case, pins


def run(case, pins, token, task_started):
    go, auth, ref = pins.get("root_go", {}), case.get("execution_authorization", {}), case["input_result"]
    if (auth.get("authorized") is not True or auth.get("authorized_by") != "ROOT"
            or pins.get("execution_authorized") is not True or go.get("authorized") is not True
            or not token or token != go.get("token")):
        raise RuntimeError("BLOCKED: exact root render GO is required")
    if not ref.get("sha256") or not ref["sha256"]:
        raise RuntimeError("BLOCKED: actual result SHA-256 is still pending")
    result_path = resolve(ref["path"])
    if not result_path.is_file() or sha(result_path) != ref["sha256"]:
        raise RuntimeError("BLOCKED: result file is absent or differs from its pin")
    result = json.loads(result_path.read_text(encoding="utf-8"))
    if result.get("schema") != ref["expected_schema"] or result.get("status") != ref["expected_status"]:
        raise RuntimeError("BLOCKED: result is not the accepted saved diagnostic")
    if not result.get("gates") or any(gate.get("passed") is not True for gate in result.get("gates", {}).values()):
        raise RuntimeError("BLOCKED: one or more saved numeric gates did not pass")
    out = resolve(case["limits"]["planned_output_dir"])
    out.mkdir(parents=True, exist_ok=False)
    png = out / "bridge_saved_diagnostic.png"
    cmd = [sys.executable, "-B", str(RENDERER), "--case", str(CASE), "--pins", str(PINS),
           "--output", str(png), "--root-go", token]
    before = {"case": sha(CASE), "pins": sha(PINS), "renderer": sha(RENDERER),
              "runner": sha(Path(__file__).resolve()), "result": sha(result_path)}
    started, timed_out, rc, stdout, stderr = time.perf_counter(), False, None, b"", b""
    try:
        remaining = max(0.0, case["limits"]["max_wall_seconds"] - (started - task_started))
        proc = subprocess.run(cmd, cwd=ROOT, capture_output=True, timeout=remaining)
        rc, stdout, stderr = proc.returncode, proc.stdout, proc.stderr
    except subprocess.TimeoutExpired as exc:
        timed_out, stdout, stderr = True, exc.stdout or b"", exc.stderr or b""
    elapsed = time.perf_counter() - started
    (out / "stdout.bin").write_bytes(stdout); (out / "stderr.bin").write_bytes(stderr)
    after = {"case": sha(CASE), "pins": sha(PINS), "renderer": sha(RENDERER),
             "runner": sha(Path(__file__).resolve()), "result": sha(result_path)}
    files = []
    for path in sorted(p for p in out.iterdir() if p.is_file() and p.name != "capture.json"):
        files.append({"path": path.name, "bytes": path.stat().st_size, "sha256": sha(path)})
    capture = {"schema": "mixed93_saved_visual_render_capture.v1", "status": "RENDERED" if rc == 0 and not timed_out and before == after else "NOT_ACCEPTED",
               "return_code": rc, "timed_out": timed_out, "wall_seconds": time.perf_counter() - task_started,
               "subprocess_seconds": elapsed, "input_path": ref["path"],
               "input_sha256": ref["sha256"], "pins_pre": before, "pins_post": after,
               "outputs_excluding_capture": files, "byte_limit": MAX_OUTPUT,
               "root_go_authorized_by": auth.get("authorized_by"), "no_retry": True}
    blob = (json.dumps(capture, indent=2, sort_keys=True, allow_nan=False) + "\n").encode("utf-8")
    if sum(item["bytes"] for item in files) + len(blob) > MAX_OUTPUT:
        capture["status"] = "BLOCKED_OUTPUT_BYTE_LIMIT"
        blob = (json.dumps(capture, indent=2, sort_keys=True, allow_nan=False) + "\n").encode("utf-8")
    (out / "capture.json").write_bytes(blob)
    total = sum(p.stat().st_size for p in out.iterdir() if p.is_file())
    if time.perf_counter() - task_started > case["limits"]["max_wall_seconds"] or total > MAX_OUTPUT:
        raise RuntimeError("render capture exceeded its fixed wall/output budget; no retry permitted")
    return 0 if capture["status"] == "RENDERED" else 1


def main():
    task_started = time.perf_counter()
    parser = argparse.ArgumentParser()
    parser.add_argument("--render", action="store_true")
    parser.add_argument("--root-go")
    args = parser.parse_args()
    case, pins = ready()
    if not args.render:
        state = "SOURCE_READY_ACTUAL_RESULT_PENDING" if case["input_result"]["sha256"] is None else "SOURCE_READY_RESULT_PINNED_NO_RENDER_REQUESTED"
        print(json.dumps({"status": state, "result_sha256": case["input_result"]["sha256"],
                          "matplotlib_imported": False, "rendered": False}, sort_keys=True))
        return 0
    return run(case, pins, args.root_go, task_started)


if __name__ == "__main__":
    raise SystemExit(main())
