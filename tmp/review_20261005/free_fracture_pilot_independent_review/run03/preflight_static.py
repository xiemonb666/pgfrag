#!/usr/bin/env python3
"""Static preflight for isolated run03 after one-based event-pair indexing fix."""
from __future__ import annotations
import ast, hashlib, json, tomllib
from pathlib import Path

RUN_DIR = Path(__file__).resolve().parent
BASE_DIR = RUN_DIR.parent
REPO_ROOT = RUN_DIR.parents[4]
EXPECTED_RUN02_ADAPTER_SHA256 = "14179191572b25c7acdf50699c040a9493461f112d9a0b9305df0f6a2d162826"
EXPECTED_RUN02_AUDIT_SHA256 = "6e7be91780f6f4d92256e3bda6ad44880cc703dde1762a8996430f93454ca9a4"
EXPECTED_SUMMARY_SHA256 = "1e3396d30aa124a3371b47bcd6cc40eea4118312a18d685bbfbe3fdaa60a9c8c"
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

def main() -> int:
    errors: list[str] = []
    run02 = BASE_DIR / "run02"
    run02_adapter = run02 / "pilot_independent_review.py"
    adapter_path = RUN_DIR / "pilot_independent_review.py"
    old_adapter = run02_adapter.read_text(encoding="utf-8")
    adapter = adapter_path.read_text(encoding="utf-8")
    old = "    for edge, (i0, j0) in enumerate(pairs):\n        xi, xj = int(grid_indexes[0][i0]), int(grid_indexes[0][j0])"
    new = ("    for edge, (point_i, point_j) in enumerate(pairs):\n"
           "        point_i0, point_j0 = int(point_i) - 1, int(point_j) - 1\n"
           "        xi, xj = int(grid_indexes[0][point_i0]), int(grid_indexes[0][point_j0])")
    require_old_count = old_adapter.count(old) == 1
    require_new_count = adapter.count(new) == 1
    errors += [] if require_old_count else ["run02 defect anchor count is not exactly one"]
    errors += [] if require_new_count else ["run03 one-based-to-zero-based fix missing or duplicated"]
    expected = old_adapter.replace(old, new, 1)
    errors += [] if adapter == expected else ["run03 adapter differs from run02 by more than exact spatial index correction"]
    errors += [] if sha(run02_adapter) == EXPECTED_RUN02_ADAPTER_SHA256 else ["run02 adapter preservation SHA mismatch"]
    run02_audit_path = run02 / "audit_output/pilot_independent_review.json"
    errors += [] if sha(run02_audit_path) == EXPECTED_RUN02_AUDIT_SHA256 else ["run02 failed audit JSON preservation SHA mismatch"]
    run02_audit_doc = json.loads(run02_audit_path.read_text(encoding="utf-8"))
    errors += [] if run02_audit_doc.get("status") == "FAIL" else ["run02 error evidence is no longer the preserved failed run"]
    errors += [] if 'OUTPUT_DIR = SCRIPT_DIR / "audit_output"' in adapter else ["run03 output path is not isolated"]
    tree = ast.parse(adapter)
    bad_ufunc_axis = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute):
            if isinstance(node.func.value, ast.Name) and node.func.value.id == "np" and node.func.attr in {"minimum", "maximum"} and any(k.arg == "axis" for k in node.keywords):
                bad_ufunc_axis.append({"function": f"np.{node.func.attr}", "line": node.lineno})
    errors += [] if not bad_ufunc_axis else [f"NumPy ufunc minimum/maximum still receives axis: {bad_ufunc_axis}"]
    index_routes = [
        {"line": n, "text": line.strip()}
        for n, line in enumerate(adapter.splitlines(), 1)
        if any(token in line for token in ("pairs[", "enumerate(pairs)", "grid_indexes[0][", "edge_i0 =", "edge_j0 ="))
    ]
    spatial_block = adapter[adapter.index("right_layer_to_layer_steps:"):adapter.index("right_layer_pair_rows =")]
    errors += [] if "int(point_i) - 1, int(point_j) - 1" in spatial_block else ["spatial loop does not normalize one-based pair IDs exactly once"]
    errors += [] if "grid_indexes[0][point_i0]" in spatial_block and "grid_indexes[0][point_j0]" in spatial_block else ["spatial loop does not index grid arrays via normalized zero-based IDs"]
    errors += [] if "grid_indexes[0][i0]" not in spatial_block and "grid_indexes[0][j0]" not in spatial_block else ["spatial loop still has unnormalized direct pair indexing"]
    source_freeze_path = REPO_ROOT / "pgfrag/experiments/free_fracture_pilot_20261005/source_freeze.toml"
    source_freeze_raw = source_freeze_path.read_bytes()
    source_freeze = tomllib.loads(source_freeze_raw.decode("utf-8"))
    errors += [] if len(source_freeze.get("pins", [])) == source_freeze.get("pin_count") == 84 else ["source freeze pin count is not 84"]
    experiment = REPO_ROOT / "pgfrag/experiments/free_fracture_pilot_20261005"
    driver = bind(experiment / "run_free_fracture_pilot.jl")
    case = bind(experiment / "fixed_case.toml")
    errors += [] if case["sha256"] == EXPECTED_CASE_SHA256 else ["fixed case SHA mismatch"]
    launcher = bind(experiment / "run_once.ps1")
    errors += [] if driver["sha256"] == EXPECTED_DRIVER_SHA256 else ["native driver pin changed"]
    errors += [] if case["sha256"] == EXPECTED_CASE_SHA256 else ["fixed case pin changed"]
    errors += [] if launcher["sha256"] == EXPECTED_LAUNCHER_SHA256 else ["launcher pin changed"]
    inventory_base = BASE_DIR / "run02/pilot_input_hash_inventory.json"
    inventory_copy = RUN_DIR / "pilot_input_hash_inventory.json"
    errors += [] if inventory_base.read_bytes() == inventory_copy.read_bytes() else ["frozen native output inventory copy changed"]
    inventory_doc = json.loads(inventory_copy.read_text(encoding="utf-8"))
    summary_row = next((f for f in inventory_doc["files"] if f["name"] == "summary.toml"), None)
    errors += [] if summary_row and summary_row["sha256"].lower() == EXPECTED_SUMMARY_SHA256 else ["frozen final summary input SHA mismatch"]
    report = {
        "schema": "pgfrag.free_fracture_pilot.run03_static_preflight.v1",
        "status": "PASS" if not errors else "FAIL",
        "run02_adapter_preserved": bind(run02_adapter),
        "run02_audit_json": bind(run02 / "audit_output/pilot_independent_review.json"),
        "run02_failure_witness": "point ID 4500 reached NumPy array index 4500 in a length-4500 grid index array; native outputs remained unmodified.",
        "declared_run03_change": "In right-grip x-layer association, each 1-based canonical pair point ID is converted once with int(point_id)-1 and only normalized indices address grid_indexes.",
        "pairs_numpy_index_routes_reviewed": index_routes,
        "invalid_numpy_ufunc_axis_calls": bad_ufunc_axis,
        "native_driver": driver,
        "fixed_case": case,
        "launcher": launcher,
        "source_freeze_pin_count": len(source_freeze.get("pins", [])),
        "adapter": bind(adapter_path),
        "runner": bind(RUN_DIR / "run_bounded_review.ps1"),
        "snapshot_helper": bind(RUN_DIR / "snapshot_bindings.py"),
        "input_inventory": bind(inventory_copy),
        "errors": errors,
    }
    target = RUN_DIR / "RUN03_STATIC_PREFLIGHT.json"
    target.write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"RUN03_STATIC_PREFLIGHT_{report['status']} pair_routes={len(index_routes)} ufunc_axis_calls={len(bad_ufunc_axis)} pins={report['source_freeze_pin_count']}")
    return 0 if not errors else 2

if __name__ == "__main__":
    raise SystemExit(main())