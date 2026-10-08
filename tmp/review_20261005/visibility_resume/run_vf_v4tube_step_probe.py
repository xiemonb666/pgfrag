"""Run one read-only canonical STEP visibility probe for the existing VFv4 asset."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import sys


REVIEW_DIR = Path(__file__).resolve().parent
WORKSPACE_ROOT = REVIEW_DIR.parents[3]
PGFRAG_ROOT = WORKSPACE_ROOT / "pgfrag"
ASSET_ROOT = (
    WORKSPACE_ROOT
    / "vesselforge"
    / "tmp"
    / "review_20261004"
    / "wall_mesh_binding"
    / "v4_run"
    / "cylinder_asset_run01"
)
STEP_PATH = ASSET_ROOT / "geometry" / "solid.step"
EXPECTED_STEP_SHA256 = "5c68bbc29b2b0bc684cc6c551c2769872c4c490d588b9a3ea2cbb98e7f74135e"
ASSET_MANIFEST_PATH = ASSET_ROOT / "manifest.json"
WALL_REFERENCE_PATH = ASSET_ROOT / "metadata" / "wall_reference.json"
RUN_DIR = REVIEW_DIR / "vf_v4tube_step_probe"
RESULT_PATH = RUN_DIR / "probe_result.json"
RECEIPT_PATH = RUN_DIR / "execution_receipt.json"
SOURCE_PATH = PGFRAG_ROOT / "src" / "pgfrag" / "bond_visibility.py"
TEST_PATH = PGFRAG_ROOT / "tests" / "test_bond_visibility.py"
FRAME = "right_handed_z_up"
SOURCE_ACCURACY_MM = 1.0e-7
SOURCE_ACCURACY_BASIS = (
    "declared modeler CAD test resolution 1e-7 mm for this loader probe only; "
    "not a measured physical material accuracy"
)
DISTANCE_TOLERANCE_MM = 1.0e-4


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def pin(path: Path) -> dict[str, object]:
    return {
        "path": str(path.resolve()),
        "sha256": sha256(path),
        "byte_length": path.stat().st_size,
    }


def main() -> int:
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
    sys.dont_write_bytecode = True
    sys.path.insert(0, str(PGFRAG_ROOT / "src"))

    import OCP
    import cadquery
    import pgfrag.bond_visibility as bond_visibility

    started_utc = datetime.now(timezone.utc).isoformat()
    source_test_before = {"source": pin(SOURCE_PATH), "test": pin(TEST_PATH)}
    asset_before = pin(STEP_PATH) if STEP_PATH.exists() else None
    asset_manifest = pin(ASSET_MANIFEST_PATH)
    wall_reference = pin(WALL_REFERENCE_PATH)
    asset_manifest_data = json.loads(ASSET_MANIFEST_PATH.read_text(encoding="utf-8"))
    manifest_units = asset_manifest_data.get("units", {})
    manifest_step = asset_manifest_data.get("files", {}).get("geometry/solid.step", {})
    result_data: dict[str, object] = {
        "probe": "vf_v4tube_canonical_step_two_pair_visibility",
        "started_utc": started_utc,
        "asset": {
            "path": str(STEP_PATH.resolve()),
            "expected_sha256": EXPECTED_STEP_SHA256,
            "pre_import_pin": asset_before,
            "manifest_pin": asset_manifest,
            "wall_reference_pin": wall_reference,
            "manifest_length_unit": manifest_units.get("length"),
            "manifest_frame": asset_manifest_data.get("frame"),
            "manifest_step_sha256": manifest_step.get("sha256"),
        },
        "declarations": {
            "source_length_unit": "mm",
            "frame": FRAME,
            "source_accuracy_mm": SOURCE_ACCURACY_MM,
            "source_accuracy_basis": SOURCE_ACCURACY_BASIS,
            "distance_tolerance_mm": DISTANCE_TOLERANCE_MM,
            "material_domain_id": "solid-0",
        },
        "runtime": {
            "python_executable": sys.executable,
            "python_version": sys.version,
            "platform": platform.platform(),
            "cadquery_version": getattr(cadquery, "__version__", "unknown"),
            "ocp_version": getattr(OCP, "__version__", "unknown"),
            "module_file": str(Path(bond_visibility.__file__).resolve()),
            "algorithm": bond_visibility.BOND_VISIBILITY_ALGORITHM,
            "runner_argv": list(sys.orig_argv),
            "python_dont_write_bytecode": sys.dont_write_bytecode,
        },
        "source_test_pins_before": source_test_before,
        "load_canonical_step_call_count": 0,
        "classify_candidate_pairs_call_count": 0,
        "probe_status": "NOT_STARTED",
        "source": None,
        "result": None,
        "error": None,
    }

    exit_code = 1
    stage = "asset_hash_preflight"
    try:
        if asset_before is None or asset_before["sha256"] != EXPECTED_STEP_SHA256:
            result_data["probe_status"] = "STEP_SHA256_MISMATCH"
        elif (
            manifest_step.get("sha256") != EXPECTED_STEP_SHA256
            or manifest_units.get("length") != "mm"
            or asset_manifest_data.get("frame") != FRAME
        ):
            result_data["probe_status"] = "ASSET_MANIFEST_DECLARATION_MISMATCH"
        else:
            stage = "load_canonical_step"
            result_data["load_canonical_step_call_count"] = 1
            source = bond_visibility.load_canonical_step(
                STEP_PATH,
                expected_sha256=EXPECTED_STEP_SHA256,
                source_length_unit="mm",
                frame=FRAME,
                source_accuracy_mm=SOURCE_ACCURACY_MM,
                source_accuracy_basis=SOURCE_ACCURACY_BASIS,
                material_domain_ids=("solid-0",),
            )
            result_data["source"] = {
                "status": source.status,
                "source_origin": source.source_origin,
                "source_geometry_sha256": source.source_geometry_sha256,
                "expected_sha256": source.expected_sha256,
                "source_length_unit": source.source_length_unit,
                "frame": source.frame,
                "source_accuracy_mm": source.source_accuracy_mm,
                "source_accuracy_basis": source.source_accuracy_basis,
                "solid_count": len(source.solids),
                "solid_domain_ids": [solid.material_domain_id for solid in source.solids],
                "unknown_reasons": list(source.unknown_reasons),
                "backend_versions": dict(source.backend_versions),
            }
            if source.status != "READY" or len(source.solids) != 1:
                result_data["probe_status"] = "SOURCE_NOT_READY_OR_NOT_ONE_SOLID"
            else:
                stage = "classify_candidate_pairs"
                result_data["classify_candidate_pairs_call_count"] = 1
                result = bond_visibility.classify_candidate_pairs(
                    source,
                    reference_point_ids=[1, 2, 3, 4],
                    reference_positions_mm=[
                        (9.0, 0.0, 2.0),
                        (9.0, 0.0, 8.0),
                        (9.0, 0.0, 5.0),
                        (-9.0, 0.0, 5.0),
                    ],
                    point_material_domain_ids=["solid-0"] * 4,
                    candidate_pairs=[(1, 2), (3, 4)],
                    reference_frame=FRAME,
                    distance_tolerance_mm=DISTANCE_TOLERANCE_MM,
                )
                result_data["result"] = result
                stage = "result_json_validation"
                pairs = result["pairs"]
                expected_pairs = {
                    (1, 2): "IN_MATERIAL",
                    (3, 4): "OUT_OR_VOID",
                }
                matches = len(pairs) == 2
                for pair in pairs:
                    key = tuple(pair["pair"])
                    matches = matches and pair["classification"] == expected_pairs.get(key)
                    matches = matches and pair["reciprocal"]["same_classification"] is True
                    matches = matches and all(
                        row["state"] == pair["classification"]
                        for row in pair["reciprocal"]["directed_records"]
                    )
                json.dumps(result_data, ensure_ascii=False, allow_nan=False)
                result_data["probe_status"] = "EXPECTED_PAIR_CLASSIFICATIONS" if matches else "UNEXPECTED_PAIR_CLASSIFICATIONS"
                exit_code = 0 if matches else 1
    except Exception as exc:
        result_data["probe_status"] = f"{stage.upper()}_EXCEPTION"
        result_data["result"] = None
        result_data["error"] = {"stage": stage, "type": type(exc).__name__, "message": str(exc)}

    finished_utc = datetime.now(timezone.utc).isoformat()
    result_data["finished_utc"] = finished_utc
    source_test_after = {"source": pin(SOURCE_PATH), "test": pin(TEST_PATH)}
    asset_after = pin(STEP_PATH) if STEP_PATH.exists() else None
    result_data["source_test_pins_after"] = source_test_after
    result_data["asset"]["post_probe_pin"] = asset_after
    result_data["scope"] = {
        "step_import_count": result_data["load_canonical_step_call_count"],
        "step_export_count": 0,
        "mesh_tessellation_count": 0,
        "native_solver_or_graph_mutation": False,
        "pair_count_requested": 2,
        "retry_count": 0,
    }
    RESULT_PATH.write_text(
        json.dumps(result_data, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    receipt = {
        "probe": result_data["probe"],
        "started_utc": started_utc,
        "finished_utc": finished_utc,
        "exit_code": exit_code,
        "retry_count": 0,
        "load_canonical_step_call_count": result_data["load_canonical_step_call_count"],
        "classify_candidate_pairs_call_count": result_data["classify_candidate_pairs_call_count"],
        "probe_status": result_data["probe_status"],
        "source_test_pins_before": source_test_before,
        "source_test_pins_after": source_test_after,
        "asset_pre_pin": asset_before,
        "asset_post_pin": asset_after,
        "result": pin(RESULT_PATH),
        "runtime": result_data["runtime"],
        "runner_argv": list(sys.orig_argv),
        "declarations": result_data["declarations"],
        "scope": result_data["scope"],
        "error": result_data["error"],
    }
    RECEIPT_PATH.write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
