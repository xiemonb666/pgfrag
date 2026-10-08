"""Run the root-authorized, bounded v3 bond-visibility pytest selection once."""

from __future__ import annotations

from contextlib import redirect_stderr, redirect_stdout
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import sys
import xml.etree.ElementTree as ET


REVIEW_DIR = Path(__file__).resolve().parent
WORKSPACE_ROOT = REVIEW_DIR.parents[3]
PGFRAG_ROOT = WORKSPACE_ROOT / "pgfrag"
SOURCE_PATH = PGFRAG_ROOT / "src" / "pgfrag" / "bond_visibility.py"
TEST_PATH = PGFRAG_ROOT / "tests" / "test_bond_visibility.py"
RUN_DIR = REVIEW_DIR / "v3_run"
JUNIT_PATH = RUN_DIR / "junit.xml"
STDOUT_PATH = RUN_DIR / "stdout.txt"
STDERR_PATH = RUN_DIR / "stderr.txt"
SIDECAR_PATH = RUN_DIR / "geometry_classification_sidecar.json"
RECEIPT_PATH = RUN_DIR / "execution_receipt.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_pin(path: Path) -> dict[str, object]:
    return {"path": str(path.resolve()), "sha256": sha256(path), "byte_length": path.stat().st_size}


def matching_pyc_files(directory: Path, stem: str) -> list[str]:
    if not directory.exists():
        return []
    return sorted(str(path.resolve()) for path in directory.glob(f"{stem}*.pyc"))


def junit_summary(path: Path) -> dict[str, object]:
    if not path.exists():
        return {"exists": False}
    root = ET.parse(path).getroot()
    suites = [root] if root.tag == "testsuite" else list(root.findall("testsuite"))
    return {
        "exists": True,
        "tests": sum(int(suite.attrib.get("tests", "0")) for suite in suites),
        "failures": sum(int(suite.attrib.get("failures", "0")) for suite in suites),
        "errors": sum(int(suite.attrib.get("errors", "0")) for suite in suites),
        "skipped": sum(int(suite.attrib.get("skipped", "0")) for suite in suites),
        "suites": [dict(suite.attrib) for suite in suites],
    }


def artifact_pin(path: Path) -> dict[str, object] | None:
    if not path.exists():
        return None
    return file_pin(path)


def main() -> int:
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
    os.environ["PYTHONPATH"] = os.pathsep.join(
        (str(PGFRAG_ROOT / "src"), str(REVIEW_DIR))
    )
    os.environ["VISIBILITY_RESUME_SIDECAR"] = str(SIDECAR_PATH.resolve())
    sys.dont_write_bytecode = True
    sys.path.insert(0, str(REVIEW_DIR))
    sys.path.insert(0, str(PGFRAG_ROOT / "src"))

    import OCP
    import cadquery
    import pytest
    import pgfrag.bond_visibility as bond_visibility

    source_pyc_dir = SOURCE_PATH.parent / "__pycache__"
    test_pyc_dir = TEST_PATH.parent / "__pycache__"
    pre_pins = {"source": file_pin(SOURCE_PATH), "test": file_pin(TEST_PATH)}
    pre_pyc = {
        "bond_visibility": matching_pyc_files(source_pyc_dir, "bond_visibility"),
        "test_bond_visibility": matching_pyc_files(test_pyc_dir, "test_bond_visibility"),
    }
    started_utc = datetime.now(timezone.utc).isoformat()
    pytest_args = [
        "-p",
        "no:cacheprovider",
        "-p",
        "capture_visibility",
        f"--junitxml={JUNIT_PATH}",
        "pgfrag/tests/test_bond_visibility.py::BondVisibilityPureMathRegressionTests",
        "pgfrag/tests/test_bond_visibility.py::BondVisibilityGeometryTests::test_hollow_shell_endpoints_inside_but_segment_crosses_void",
        "pgfrag/tests/test_bond_visibility.py::BondVisibilityGeometryTests::test_hollow_shell_preserves_all_outer_and_inner_hits",
        "pgfrag/tests/test_bond_visibility.py::BondVisibilityGeometryTests::test_same_wall_and_solid_coupon_remain_inside_material",
    ]

    with STDOUT_PATH.open("w", encoding="utf-8", newline="\n") as stdout_stream:
        with STDERR_PATH.open("w", encoding="utf-8", newline="\n") as stderr_stream:
            with redirect_stdout(stdout_stream), redirect_stderr(stderr_stream):
                print(
                    json.dumps(
                        {
                            "python_executable": sys.executable,
                            "python_version": sys.version,
                            "platform": platform.platform(),
                            "pytest_version": pytest.__version__,
                            "cadquery_version": getattr(cadquery, "__version__", "unknown"),
                            "ocp_version": getattr(OCP, "__version__", "unknown"),
                            "module_file": str(Path(bond_visibility.__file__).resolve()),
                            "algorithm": bond_visibility.BOND_VISIBILITY_ALGORITHM,
                            "argv": pytest_args,
                            "environment": {
                                "PYTHONPATH": os.environ["PYTHONPATH"],
                                "PYTHONDONTWRITEBYTECODE": os.environ["PYTHONDONTWRITEBYTECODE"],
                                "pytest_cache_provider_disabled": True,
                            },
                        },
                        ensure_ascii=False,
                    )
                )
                exit_code = int(pytest.main(pytest_args))

    finished_utc = datetime.now(timezone.utc).isoformat()
    post_pins = {"source": file_pin(SOURCE_PATH), "test": file_pin(TEST_PATH)}
    post_pyc = {
        "bond_visibility": matching_pyc_files(source_pyc_dir, "bond_visibility"),
        "test_bond_visibility": matching_pyc_files(test_pyc_dir, "test_bond_visibility"),
    }
    try:
        junit = junit_summary(JUNIT_PATH)
    except Exception as exc:
        junit = {"exists": JUNIT_PATH.exists(), "parse_error": type(exc).__name__}
    sidecar: dict[str, object] = {"exists": False}
    if SIDECAR_PATH.exists():
        sidecar = file_pin(SIDECAR_PATH)
        try:
            sidecar_data = json.loads(SIDECAR_PATH.read_text(encoding="utf-8"))
            sidecar["captured_call_count"] = sidecar_data.get("captured_call_count")
            sidecar["pytest_exit_code"] = sidecar_data.get("pytest_exit_code")
        except Exception as exc:
            sidecar["parse_error"] = type(exc).__name__
    receipt = {
        "run": "bond_visibility_v3_bounded_pytest",
        "started_utc": started_utc,
        "finished_utc": finished_utc,
        "exit_code": exit_code,
        "retry_count": 0,
        "runtime": {
            "python_executable": sys.executable,
            "python_version": sys.version,
            "platform": platform.platform(),
            "pytest_version": pytest.__version__,
            "cadquery_version": getattr(cadquery, "__version__", "unknown"),
            "ocp_version": getattr(OCP, "__version__", "unknown"),
            "module_file": str(Path(bond_visibility.__file__).resolve()),
            "algorithm": bond_visibility.BOND_VISIBILITY_ALGORITHM,
            "argv": pytest_args,
            "environment": {
                "PYTHONPATH": os.environ["PYTHONPATH"],
                "PYTHONDONTWRITEBYTECODE": os.environ["PYTHONDONTWRITEBYTECODE"],
                "pytest_cache_provider_disabled": True,
                "python_dont_write_bytecode": sys.dont_write_bytecode,
            },
        },
        "source_test_pins_before": pre_pins,
        "source_test_pins_after": post_pins,
        "target_pyc_files_before": pre_pyc,
        "target_pyc_files_after": post_pyc,
        "junit": junit,
        "stdout": artifact_pin(STDOUT_PATH),
        "stderr": artifact_pin(STDERR_PATH),
        "geometry_classification_sidecar": sidecar,
        "scope": {
            "step_export_or_import": False,
            "native_solver": False,
            "test_selection": pytest_args[5:],
        },
    }
    RECEIPT_PATH.write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
