"""Focused serialization tests; importing these must never invoke CAD/OCP."""

from __future__ import annotations

import builtins
import importlib.util
from pathlib import Path
from types import ModuleType

import pytest


WORKER_PATH = Path(__file__).resolve().parents[1] / "worker.py"
SPEC = importlib.util.spec_from_file_location("forge_volume_ingress_worker", WORKER_PATH)
assert SPEC is not None and SPEC.loader is not None
worker = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(worker)


def test_total_byte_fixed_point_crosses_decimal_digit_boundary() -> None:
    report: dict[str, object] = {"schema": "test.report.v1", "output": {"total_bytes": 0}}
    for padding_size in range(12_000):
        report["padding"] = "x" * padding_size
        if len(worker._json_bytes(report)) == 9_999:
            break
    else:
        pytest.fail("could not construct the 9,999-byte pre-count boundary fixture")

    payload = worker._serialize_report_with_total(report, non_report_bytes=0)
    final_total = report["output"]["total_bytes"]
    assert 10_000 <= final_total < 100_000
    assert len(payload) == final_total
    assert payload == worker._json_bytes(report)


def test_report_writer_blocks_when_total_exceeds_output_budget(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(worker, "MAX_OUTPUT_BYTES", 1_024)
    monkeypatch.setattr(worker, "PARENT_RECEIPT_RESERVE_BYTES", 128)
    (tmp_path / "payload.bin").write_bytes(b"p" * 256)
    report = {"schema": "test.report.v1", "output": {"total_bytes": 0}, "padding": "x" * 2_000}

    with pytest.raises(worker.IngressBlocked, match="exceed"):
        worker._write_final_report(tmp_path, report)

    assert not (tmp_path / "ingress_report.json").exists()


def test_report_serialization_does_not_import_cad_or_ocp(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    (tmp_path / "payload.bin").write_bytes(b"raw")
    report = {"schema": "test.report.v1", "output": {"total_bytes": 0}}
    original_import = builtins.__import__

    def deny_geometry_import(name: str, *args: object, **kwargs: object) -> ModuleType:
        if name == "cadquery" or name == "OCP" or name.startswith("OCP."):
            raise AssertionError(f"unexpected geometry import during report serialization: {name}")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", deny_geometry_import)
    worker._write_final_report(tmp_path, report)

    saved = worker.json.loads((tmp_path / "ingress_report.json").read_text(encoding="utf-8"))
    actual_total = sum(path.stat().st_size for path in tmp_path.iterdir() if path.is_file())
    assert saved["output"]["total_bytes"] == actual_total
