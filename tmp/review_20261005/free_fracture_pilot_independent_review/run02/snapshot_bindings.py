#!/usr/bin/env python3
"""Capture read-only source and completed-attempt byte/hash pins."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import tomllib
from pathlib import Path


REVIEW_DIR = Path(__file__).resolve().parent
REPO_ROOT = REVIEW_DIR.parents[4]
EXPERIMENT_DIR = REPO_ROOT / "pgfrag" / "experiments" / "free_fracture_pilot_20261005"
ATTEMPT_DIR = EXPERIMENT_DIR / "run_attempt_01"
RESULTS_DIR = ATTEMPT_DIR / "results"
EXPECTED_FREEZE_SHA256 = "1f4c3e9b10a607d303bbc3435bf22ae8b591ce8f6b5651aacb6cba3c3bfd7217"
EXPECTED_RECEIPT_SHA256 = "0bae8302c765effdd2614f72a2e1802d4779aa9ea4910829a1281f063c4eb557"
EXPECTED_SUMMARY_SHA256 = "1e3396d30aa124a3371b47bcd6cc40eea4118312a18d685bbfbe3fdaa60a9c8c"
EXPECTED_CASE_SHA256 = "45b9f8bb5d10522fa2e6fa91ca8d944e2d6af27396034b489e27a026783bbedd"
EXPECTED_CASE_OUTPUT_CAP = 199_229_440
EXPECTED_CASE_ATTEMPT_CAP = 209_715_200
EXPECTED_RESULT_NAMES = {
    "point_state_history.f64le", "point_state_layout.toml",
    "graph_component_labels.i32le", "graph_component_labels_layout.toml",
    "point_reference.csv", "canonical_bonds.csv", "native_failure_events.csv",
    "graph_state_history.csv", "energy_reaction_history.csv", "summary.toml",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def bind_file(path: Path, relative: str | None = None) -> dict[str, object]:
    return {
        "path": relative if relative is not None else str(path.resolve()),
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("pre", "post"))
    args = parser.parse_args()

    freeze_path = EXPERIMENT_DIR / "source_freeze.toml"
    freeze_raw = freeze_path.read_bytes()
    freeze_sha = hashlib.sha256(freeze_raw).hexdigest()
    if freeze_sha != EXPECTED_FREEZE_SHA256:
        raise SystemExit(f"source freeze changed: {freeze_sha}")
    freeze = tomllib.loads(freeze_raw.decode("utf-8"))
    pins = freeze.get("pins", [])
    if len(pins) != 84 or freeze.get("pin_count") != 84:
        raise SystemExit(f"source pin count mismatch: {len(pins)}")
    source_pins = []
    for pin in pins:
        rel = str(pin["path"])
        path = (REPO_ROOT / Path(rel.replace("/", os.sep))).resolve(strict=True)
        if REPO_ROOT.resolve() not in path.parents:
            raise SystemExit(f"source pin escapes repository root: {rel}")
        observed = bind_file(path, rel)
        observed["matches_freeze"] = (
            observed["bytes"] == int(pin["bytes"])
            and observed["sha256"] == str(pin["sha256"]).lower()
        )
        source_pins.append(observed)
    mismatches = [p for p in source_pins if not p["matches_freeze"]]

    case_path = EXPERIMENT_DIR / "fixed_case.toml"
    case_raw = case_path.read_bytes()
    case_sha = hashlib.sha256(case_raw).hexdigest()
    case = tomllib.loads(case_raw.decode("utf-8"))
    limits = case["limits"]
    case_limits_match = (
        case_sha == EXPECTED_CASE_SHA256
        and int(limits["max_output_bytes"]) == EXPECTED_CASE_OUTPUT_CAP
        and int(limits["max_attempt_bytes"]) == EXPECTED_CASE_ATTEMPT_CAP
    )

    receipt_path = ATTEMPT_DIR / "launcher_receipt.json"
    receipt_sha = sha256_file(receipt_path)
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    if receipt_sha != EXPECTED_RECEIPT_SHA256:
        raise SystemExit(f"final receipt changed: {receipt_sha}")
    if receipt.get("status") != "COMPLETED_EXIT_0":
        raise SystemExit(f"completed final receipt required: {receipt.get('status')}")
    if str(receipt.get("summary_sha256", "")).lower() != EXPECTED_SUMMARY_SHA256:
        raise SystemExit("final receipt summary SHA does not match the authorized final output")
    summary_path = RESULTS_DIR / "summary.toml"
    if sha256_file(summary_path) != EXPECTED_SUMMARY_SHA256:
        raise SystemExit("final summary bytes changed")

    result_paths = sorted((p for p in RESULTS_DIR.iterdir() if p.is_file()), key=lambda p: p.name)
    observed_names = {p.name for p in result_paths}
    if observed_names != EXPECTED_RESULT_NAMES:
        raise SystemExit(f"result file inventory changed: {sorted(observed_names)}")
    inventory_path = REVIEW_DIR / "pilot_input_hash_inventory.json"
    expected_inventory = json.loads(inventory_path.read_text(encoding="utf-8"))
    expected_file_map = {str(f["name"]): f for f in expected_inventory["files"]}
    result_files = [bind_file(p, p.name) for p in result_paths]
    for item in result_files:
        expected = expected_file_map[item["path"]]
        item["matches_frozen_input_inventory"] = (
            item["bytes"] == int(expected["bytes"])
            and item["sha256"] == str(expected["sha256"]).lower()
        )
    result_mismatches = [f for f in result_files if not f["matches_frozen_input_inventory"]]
    lines = [f"{f['path']}|{f['bytes']}|{f['sha256']}" for f in result_files]
    result_manifest_sha = hashlib.sha256(("\n".join(lines) + "\n").encode("utf-8")).hexdigest()

    attempt_files = []
    for path in sorted((p for p in ATTEMPT_DIR.rglob("*") if p.is_file()),
                       key=lambda p: str(p.relative_to(ATTEMPT_DIR)).replace("\\", "/")):
        rel = str(path.relative_to(ATTEMPT_DIR)).replace("\\", "/")
        attempt_files.append(bind_file(path, rel))
    attempt_bytes = sum(int(p["bytes"]) for p in attempt_files)
    results_bytes = sum(int(f["bytes"]) for f in result_files)
    cap_checks = {
        "case_limit_fields_match_frozen_values": case_limits_match,
        "result_tree_within_case_output_cap": results_bytes <= EXPECTED_CASE_OUTPUT_CAP,
        "result_tree_within_freeze_output_cap": results_bytes <= int(freeze["max_experiment_output_bytes"]),
        "attempt_tree_within_case_attempt_cap": attempt_bytes <= EXPECTED_CASE_ATTEMPT_CAP,
    }
    cap_mismatches = [name for name, passed in cap_checks.items() if not passed]

    payload = {
        "schema": "pgfrag.free_fracture_pilot.source_output_prepost_snapshot.v1",
        "stage": args.stage,
        "scope": "SHA-256 and byte counts only; no result values decoded or interpreted",
        "snapshot_utility": bind_file(Path(__file__).resolve()),
        "source_freeze": {**bind_file(freeze_path), "pin_count": len(source_pins),
                          "all_pins_match": not mismatches, "pins": source_pins},
        "fixed_case": {**bind_file(case_path), "max_output_bytes": limits["max_output_bytes"],
                        "max_attempt_bytes": limits["max_attempt_bytes"],
                        "matches_frozen_case_and_limits": case_limits_match},
        "final_receipt": {**bind_file(receipt_path),
                          "status": receipt.get("status"),
                          "recorded_final_attempt_bytes": receipt.get("final_attempt_bytes"),
                          "summary_sha256": receipt.get("summary_sha256")},
        "result_directory": str(RESULTS_DIR.resolve()),
        "result_file_count": len(result_files),
        "result_tree_bytes": results_bytes,
        "result_tree_manifest_sha256": result_manifest_sha,
        "all_result_files_match_prepared_inventory": not result_mismatches,
        "results": result_files,
        "attempt_tree_bytes": attempt_bytes,
        "attempt_tree_files": attempt_files,
        "case_caps": cap_checks,
        "cap_mismatches": cap_mismatches,
        "mismatch_count": len(mismatches) + len(result_mismatches) + len(cap_mismatches),
        "source_mismatches": mismatches,
        "result_mismatches": result_mismatches,
    }
    target = REVIEW_DIR / f"{args.stage}_execution_snapshot.json"
    target.write_text(json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                      encoding="utf-8")
    print(f"SNAPSHOT_{args.stage.upper()} pins={len(source_pins)} results={len(result_files)} "
          f"mismatches={payload['mismatch_count']} attempt_bytes={attempt_bytes} file={target}")
    return 0 if payload["mismatch_count"] == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
