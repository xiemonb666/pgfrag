"""Run one pinned operator-focused pytest command and preserve all receipts."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import xml.etree.ElementTree as ET


ROOT = Path(__file__).resolve().parent
SNAPSHOT = ROOT / "source_snapshot"
MANIFEST = json.loads((ROOT / "source_freeze.json").read_text(encoding="utf-8"))
LIVE_SOURCE = Path(__file__).resolve().parents[4] / "experiments" / "efgm_elastic_operator_20261005" / "efgm_elastic_operator.py"
LIVE_TESTS = LIVE_SOURCE.with_name("test_efgm_elastic_operator.py")
SNAPSHOT_SOURCE = SNAPSHOT / "efgm_elastic_operator.py"
SNAPSHOT_TESTS = SNAPSHOT / "test_efgm_elastic_operator.py"
JUNIT = ROOT / "junit.xml"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def iso_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def xml_counts(path: Path) -> dict[str, int] | None:
    if not path.exists():
        return None
    root = ET.parse(path).getroot()
    suite = root if root.tag == "testsuite" else root.find("testsuite")
    if suite is None:
        return None
    return {
        key: int(suite.attrib.get(key, "0"))
        for key in ("tests", "failures", "errors", "skipped")
    }


def write_json(name: str, value: object) -> None:
    (ROOT / name).write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


pre_hashes = {
    "live_source": sha256(LIVE_SOURCE),
    "live_test": sha256(LIVE_TESTS),
    "snapshot_source": sha256(SNAPSHOT_SOURCE),
    "snapshot_test": sha256(SNAPSHOT_TESTS),
}
expected = {
    "live_source": MANIFEST["sha256"],
    "live_test": MANIFEST["test_sha256"],
    "snapshot_source": MANIFEST["sha256"],
    "snapshot_test": MANIFEST["test_sha256"],
}
pins: dict[str, object] = {
    "started_utc": iso_now(),
    "pre_run_sha256": pre_hashes,
    "expected_sha256": expected,
}
write_json("pre_run_pins.json", pins)

stdout_text = ""
stderr_text = ""
return_code = 125
timed_out = False
elapsed_seconds = 0.0
command = [
    sys.executable,
    "-m",
    "pytest",
    "-q",
    "-p",
    "no:cacheprovider",
    str(SNAPSHOT_TESTS),
    "--junitxml",
    str(JUNIT),
]

if pre_hashes != expected:
    stderr_text = "Pinned live/snapshot source hashes do not match source_freeze.json; pytest was not started.\n"
else:
    env = os.environ.copy()
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    started = time.monotonic()
    pins["pytest_started_utc"] = iso_now()
    try:
        proc = subprocess.run(
            command,
            cwd=SNAPSHOT,
            env=env,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=30.0,
            check=False,
        )
        return_code = int(proc.returncode)
        stdout_text = proc.stdout
        stderr_text = proc.stderr
    except subprocess.TimeoutExpired as exc:
        timed_out = True
        return_code = 124
        stdout_value = exc.stdout or ""
        stderr_value = exc.stderr or ""
        stdout_text = stdout_value.decode("utf-8", "replace") if isinstance(stdout_value, bytes) else stdout_value
        stderr_text = stderr_value.decode("utf-8", "replace") if isinstance(stderr_value, bytes) else stderr_value
        stderr_text += "\nFocused pytest exceeded its 30-second timeout and was terminated.\n"
    elapsed_seconds = time.monotonic() - started

(ROOT / "stdout.log").write_text(stdout_text, encoding="utf-8")
(ROOT / "stderr.log").write_text(stderr_text, encoding="utf-8")
(ROOT / "return_code.txt").write_text(f"{return_code}\n", encoding="ascii")

post_hashes = {
    "live_source": sha256(LIVE_SOURCE),
    "live_test": sha256(LIVE_TESTS),
    "snapshot_source": sha256(SNAPSHOT_SOURCE),
    "snapshot_test": sha256(SNAPSHOT_TESTS),
}
write_json(
    "post_run_pins.json",
    {
        "finished_utc": iso_now(),
        "post_run_sha256": post_hashes,
        "snapshot_unchanged": {
            "source": pre_hashes["snapshot_source"] == post_hashes["snapshot_source"],
            "test": pre_hashes["snapshot_test"] == post_hashes["snapshot_test"],
        },
    },
)
write_json(
    "run_record.json",
    {
        "status": (
            "TIMED_OUT"
            if timed_out
            else (
                "COMPLETED_PASS"
                if pre_hashes == expected and return_code == 0
                else (
                    "COMPLETED_FAIL"
                    if pre_hashes == expected
                    else "PIN_MISMATCH_NO_RUN"
                )
            )
        ),
        "command": command,
        "cwd": str(SNAPSHOT),
        "timeout_seconds": 30,
        "elapsed_seconds": elapsed_seconds,
        "return_code": return_code,
        "junit": str(JUNIT) if JUNIT.exists() else None,
        "junit_counts": xml_counts(JUNIT),
        "pytest_invocations": 1 if pre_hashes == expected else 0,
    },
)
sys.exit(return_code)
