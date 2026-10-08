"""Capture one run's selected in-memory bond-visibility geometry results."""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

import pgfrag.bond_visibility as bond_visibility


_TARGET_TESTS = {
    "test_hollow_shell_endpoints_inside_but_segment_crosses_void",
    "test_hollow_shell_preserves_all_outer_and_inner_hits",
    "test_same_wall_and_solid_coupon_remain_inside_material",
}
_active_test = ""
_capture_index = 0
_captured: list[dict[str, Any]] = []
_original_classify = bond_visibility.classify_candidate_pairs


def _capture_classification(*args: Any, **kwargs: Any) -> dict[str, Any]:
    global _capture_index
    result = _original_classify(*args, **kwargs)
    test_name = _active_test.rsplit("::", 1)[-1]
    if test_name in _TARGET_TESTS:
        # Validate and retain the exact returned ledger while the fixture is
        # already being classified; this does not construct or classify it again.
        json.dumps(result, allow_nan=False)
        _captured.append(
            {
                "test_nodeid": _active_test,
                "call_index": _capture_index,
                "result": result,
            }
        )
        _capture_index += 1
    return result


def pytest_configure(config: Any) -> None:
    bond_visibility.classify_candidate_pairs = _capture_classification


def pytest_runtest_setup(item: Any) -> None:
    global _active_test, _capture_index
    _active_test = item.nodeid
    _capture_index = 0


def pytest_sessionfinish(session: Any, exitstatus: int) -> None:
    sidecar_path = os.environ.get("VISIBILITY_RESUME_SIDECAR")
    if not sidecar_path:
        return
    payload = {
        "schema": "pgfrag.bond-visibility-geometry-sidecar.v1",
        "module_file": str(Path(bond_visibility.__file__).resolve()),
        "algorithm": bond_visibility.BOND_VISIBILITY_ALGORITHM,
        "pytest_exit_code": int(exitstatus),
        "captured_call_count": len(_captured),
        "results": _captured,
    }
    Path(sidecar_path).write_text(
        json.dumps(payload, allow_nan=False, indent=2) + "\n",
        encoding="utf-8",
    )
