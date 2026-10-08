#!/usr/bin/env python3
"""Static guard for the separately copied run02 reviewer adapter."""
from __future__ import annotations

import ast
import hashlib
import json
import tomllib
from pathlib import Path


RUN02_DIR = Path(__file__).resolve().parent
BASE_DIR = RUN02_DIR.parent
REPO_ROOT = RUN02_DIR.parents[4]
EXPECTED_BASE_ADAPTER_SHA256 = "e2d2f8cf2f9d71826c2d71bf9058546ea3ea6249ac4228a4d5f73149815c6187"
EXPECTED_FREEZE_SHA256 = "a6eaa3f0f4339b644a44d949e00fb220e85c4da1266dd149e94db88d29d437e3"
EXPECTED_CASE_SHA256 = "45b9f8bb5d10522fa2e6fa91ca8d944e2d6af27396034b489e27a026783bbedd"
EXPECTED_DRIVER_SHA256 = "9922d52c13468df5152d88c1fdb1b7c9fc6ce473a1a3a9e56507aaf37479f50b"
EXPECTED_LAUNCHER_SHA256 = "97434a0f5b17dd93ae7a0064f98515b4651aad583b0a039b7907c646c89ba69d"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def bind(path: Path) -> dict[str, object]:
    return {"path": str(path.resolve()), "bytes": path.stat().st_size, "sha256": sha(path)}


def require(condition: bool, message: str, errors: list[str]) -> None:
    if not condition:
        errors.append(message)


def main() -> int:
    errors: list[str] = []
    base_adapter = BASE_DIR / "pilot_independent_review.py"
    adapter = RUN02_DIR / "pilot_independent_review.py"
    original = base_adapter.read_text(encoding="utf-8")
    updated = adapter.read_text(encoding="utf-8")
    root_old = "REPO_ROOT = SCRIPT_DIR.parents[3]"
    root_new = "REPO_ROOT = SCRIPT_DIR.parents[4]"
    api_old = "clearance = np.minimum(lengths_m[None, :] / 2 - np.abs(xyz), axis=1)"
    api_new = "clearance = np.min(lengths_m[None, :] / 2 - np.abs(xyz), axis=1)"
    require(original.count(root_old) == 1, "baseline root-depth anchor count is not one", errors)
    require(original.count(api_old) == 1, "baseline invalid NumPy reduction anchor count is not one", errors)
    expected_updated = original.replace(root_old, root_new, 1).replace(api_old, api_new, 1)
    require(updated == expected_updated,
            "run02 adapter contains changes beyond root-depth and np.minimum-to-np.min corrections", errors)
    require(sha(base_adapter) == EXPECTED_BASE_ADAPTER_SHA256,
            "failed run01 adapter was changed from its frozen SHA", errors)
    require("OUTPUT_DIR = SCRIPT_DIR / \"audit_output\"" in updated,
            "run02 output is not isolated under run02/audit_output", errors)

    tree = ast.parse(updated)
    bad_ufunc_axis_calls = []
    expected_reduction_seen = False
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Attribute):
            continue
        func = node.func
        owner = func.value
        has_axis = any(keyword.arg == "axis" for keyword in node.keywords)
        if (isinstance(owner, ast.Name) and owner.id == "np"
                and func.attr in {"minimum", "maximum"} and has_axis):
            bad_ufunc_axis_calls.append({"function": f"np.{func.attr}", "line": node.lineno})
        if (isinstance(owner, ast.Name) and owner.id == "np" and func.attr == "min"
                and has_axis and any(keyword.arg == "axis" and isinstance(keyword.value, ast.Constant)
                                     and keyword.value.value == 1 for keyword in node.keywords)):
            if node.lineno == next((n.lineno for n in ast.walk(tree)
                                    if isinstance(n, ast.Assign)
                                    and any(isinstance(t, ast.Name) and t.id == "clearance"
                                            for t in n.targets)), -1):
                expected_reduction_seen = True
    require(not bad_ufunc_axis_calls,
            f"ufunc minimum/maximum still receives axis: {bad_ufunc_axis_calls}", errors)
    require(expected_reduction_seen, "correct np.min axis=1 clearance reduction was not found", errors)

    freeze = BASE_DIR / "PRE_AUDIT_FREEZE.json"
    require(sha(freeze) == EXPECTED_FREEZE_SHA256,
            "run01 pre-audit freeze manifest changed", errors)
    freeze_doc = json.loads(freeze.read_text(encoding="utf-8"))
    require(freeze_doc["review_artifacts"]["pilot_independent_review.py"]["sha256"]
            == EXPECTED_BASE_ADAPTER_SHA256, "run01 adapter freeze record mismatch", errors)

    runner_base = (BASE_DIR / "run_bounded_review.ps1").read_text(encoding="utf-8-sig")
    runner_new = (RUN02_DIR / "run_bounded_review.ps1").read_text(encoding="utf-8-sig")
    runner_expected = runner_base.replace(
        "Join-Path $reviewDir '..\\..\\..\\..'",
        "Join-Path $reviewDir '..\\..\\..\\..\\..'", 1)
    require(runner_base.count("Join-Path $reviewDir '..\\..\\..\\..'") == 1,
            "baseline runner root-depth anchor count is not one", errors)
    require(runner_new == runner_expected,
            "run02 runner contains unexpected edits", errors)

    snapshot_base = (BASE_DIR / "snapshot_bindings.py").read_text(encoding="utf-8")
    snapshot_new = (RUN02_DIR / "snapshot_bindings.py").read_text(encoding="utf-8")
    snapshot_expected = snapshot_base.replace(
        "REPO_ROOT = REVIEW_DIR.parents[3]", "REPO_ROOT = REVIEW_DIR.parents[4]", 1)
    require(snapshot_base.count("REPO_ROOT = REVIEW_DIR.parents[3]") == 1,
            "baseline snapshot root-depth anchor count is not one", errors)
    require(snapshot_new == snapshot_expected,
            "run02 snapshot helper contains unexpected edits", errors)

    inventory_base = BASE_DIR / "pilot_input_hash_inventory.json"
    inventory_copy = RUN02_DIR / "pilot_input_hash_inventory.json"
    require(inventory_base.read_bytes() == inventory_copy.read_bytes(),
            "run02 input inventory copy differs from run01 frozen inventory", errors)
    source_bindings = {
        "driver": bind(REPO_ROOT / "pgfrag/experiments/free_fracture_pilot_20261005/run_free_fracture_pilot.jl"),
        "fixed_case": bind(REPO_ROOT / "pgfrag/experiments/free_fracture_pilot_20261005/fixed_case.toml"),
        "launcher": bind(REPO_ROOT / "pgfrag/experiments/free_fracture_pilot_20261005/run_once.ps1"),
    }
    require(source_bindings["driver"]["sha256"] == EXPECTED_DRIVER_SHA256,
            "native driver source pin changed", errors)
    require(source_bindings["fixed_case"]["sha256"] == EXPECTED_CASE_SHA256,
            "fixed case source pin changed", errors)
    require(source_bindings["launcher"]["sha256"] == EXPECTED_LAUNCHER_SHA256,
            "launcher source pin changed", errors)
    source_freeze_path = REPO_ROOT / "pgfrag/experiments/free_fracture_pilot_20261005/source_freeze.toml"
    source_freeze = tomllib.loads(source_freeze_path.read_text(encoding="utf-8"))
    require(source_freeze.get("pin_count") == len(source_freeze.get("pins", [])) == 84,
            "native source freeze no longer declares 84 matching pins", errors)

    payload = {
        "schema": "pgfrag.free_fracture_pilot.run02_static_preflight.v1",
        "status": "PASS" if not errors else "FAIL",
        "run01_artifacts_preserved": True,
        "declared_run02_changes": [
            "REPO_ROOT path-depth adjustment for copy under run02/",
            "clearance reduction API np.minimum(..., axis=1) corrected to np.min(..., axis=1)",
            "run02-specific wrapper/snapshot-helper root-depth adjustments",
        ],
        "invalid_numpy_ufunc_axis_calls": bad_ufunc_axis_calls,
        "source_bindings": source_bindings,
        "native_pin_count": len(source_freeze.get("pins", [])),
        "artifacts": {name: bind(RUN02_DIR / name) for name in (
            "pilot_independent_review.py", "run_bounded_review.ps1",
            "snapshot_bindings.py", "pilot_input_hash_inventory.json")},
        "errors": errors,
    }
    (RUN02_DIR / "RUN02_STATIC_PREFLIGHT.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8")
    print(f"RUN02_STATIC_PREFLIGHT_{payload['status']} "
          f"invalid_ufunc_axis_calls={len(bad_ufunc_axis_calls)} pins={len(source_freeze.get('pins', []))}")
    return 0 if not errors else 2


if __name__ == "__main__":
    raise SystemExit(main())
