"""One-shot, stdlib-only process capture for the preapproved mixed93 run."""
import datetime
import hashlib
import json
import pathlib
import subprocess
import time

ROOT = pathlib.Path.cwd().resolve()
REVIEW = ROOT / "pgfrag/tmp/review_20261006/mixed93_cached_lefm_solve_20261006"
ENTRY = ROOT / "pgfrag/experiments/mixed93_cached_lefm_solve_20261006/mixed93_cached_lefm_solve.py"
CASE = ENTRY.with_name("case.json")
PINS = REVIEW / "SOURCE_PINS.json"
SNAPSHOT = REVIEW / "source_freeze_snapshot"
RUN = REVIEW / "run01"
CAPTURE = REVIEW / "execution_capture"
PYTHON = ROOT / "third_party/vesselforge-runtime/venv/Scripts/python.exe"
OUTPUT = RUN / "result.json"


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def input_hashes():
    pins = json.loads(PINS.read_text(encoding="utf-8"))
    paths = {"entry": ENTRY, "case": CASE, "source_pins_manifest": PINS}
    paths.update({key: ROOT / record["path"] for key, record in pins["references"].items()})
    return {key: sha256(path) for key, path in paths.items()}


def main():
    if RUN.exists() or OUTPUT.parent.exists():
        raise RuntimeError(f"refusing existing output directory: {RUN}")
    if not PYTHON.is_file() or not ENTRY.is_file() or not CASE.is_file() or not PINS.is_file():
        raise RuntimeError("missing approved runtime or source input")
    case = json.loads(CASE.read_text(encoding="utf-8"))
    pins = json.loads(PINS.read_text(encoding="utf-8"))
    frozen_case = SNAPSHOT / "case.unauthorized.json"
    frozen_pins = SNAPSHOT / "SOURCE_PINS.unauthorized.json"
    if case["execution_authorization"]["authorized"] is not True or pins["execution_authorized"] is not True:
        raise RuntimeError("case and source pins are not both authorized")
    if pins["entry_sha256"] != sha256(ENTRY) or pins["current_case_sha256"] != sha256(CASE):
        raise RuntimeError("entry or current authorized case SHA differs from reviewed pin")
    if pins["case_sha256_at_source_freeze"] != sha256(frozen_case):
        raise RuntimeError("historical unauthorized case snapshot does not match freeze pin")
    for key, record in pins["references"].items():
        if record["sha256"] != sha256(ROOT / record["path"]):
            raise RuntimeError(f"source reference SHA mismatch: {key}")

    CAPTURE.mkdir(exist_ok=False)
    argv = [str(PYTHON), "-B", str(ENTRY), "--execute", "--output", str(OUTPUT)]
    pre_hashes = input_hashes()
    started_utc = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
    start = time.perf_counter()
    timed_out, returncode, stdout, stderr = False, None, b"", b""
    try:
        completed = subprocess.run(argv, cwd=ROOT, capture_output=True, timeout=60, check=False)
        returncode, stdout, stderr = completed.returncode, completed.stdout, completed.stderr
    except subprocess.TimeoutExpired as exc:
        timed_out = True
        stdout, stderr = exc.stdout or b"", exc.stderr or b""
    elapsed = time.perf_counter() - start
    (CAPTURE / "stdout.bin").write_bytes(stdout)
    (CAPTURE / "stderr.bin").write_bytes(stderr)
    post_hashes = input_hashes()
    output_files = {}
    for label, path in (("result_json", OUTPUT), ("arrays_npz", RUN / "arrays.npz")):
        output_files[label] = {"exists": path.is_file()}
        if path.is_file():
            output_files[label].update(sha256=sha256(path), bytes=path.stat().st_size)
    receipt = {
        "schema": "mixed93_one_shot_execution_capture.v1",
        "started_utc": started_utc,
        "elapsed_seconds": elapsed,
        "timeout_seconds": 60,
        "timed_out": timed_out,
        "actual_process_return_code": returncode,
        "argv": argv,
        "stdout_bytes": len(stdout),
        "stderr_bytes": len(stderr),
        "stdout_text_utf8": stdout.decode("utf-8", errors="replace"),
        "stderr_text_utf8": stderr.decode("utf-8", errors="replace"),
        "stdout_file": "stdout.bin",
        "stderr_file": "stderr.bin",
        "case_authorization_delta": [
            "execution_authorization.authorized",
            "execution_authorization.authorized_by",
            "execution_authorization.authorized_utc",
        ],
        "historical_unauthorized_case_sha256": sha256(frozen_case),
        "historical_unauthorized_pins_sha256": sha256(frozen_pins),
        "pre_run_eight_pin_hashes": pre_hashes,
        "post_run_eight_pin_hashes": post_hashes,
        "input_hashes_unchanged_during_run": pre_hashes == post_hashes,
        "output_files": output_files,
    }
    (CAPTURE / "execution_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"actual_rc": returncode, "timed_out": timed_out,
                      "elapsed_seconds": elapsed, "capture": str(CAPTURE)}, sort_keys=True))


if __name__ == "__main__":
    main()
