"""Hash-only preflight for the root-authorized corrected input freeze."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
import sys


WORKSPACE = Path.cwd().resolve()
REVIEW = WORKSPACE / "pgfrag/tmp/review_20261005/forge_volume_ingress_independent_review"
FREEZE_PATH = REVIEW / "PILOT_AUDIT_INPUT_FREEZE.json"
OLD_FREEZE = REVIEW / "PILOT_AUDIT_INPUT_FREEZE_before_partial_pin_correction.json"
PRIOR_RECEIPT = REVIEW / "READONLY_AUDIT_EXECUTION_RECEIPT.json"
PRIOR_CASE_RECEIPT = REVIEW / "READONLY_AUDIT_PREFLIGHT_MISMATCH_RECEIPT.json"
OUTPUT = REVIEW / "PILOT_AUDIT_FREEZE_VALIDATION_rev5.json"
EXPECTED_ADAPTER = "827f429c141403a5e85c7c9a67376d38b0628f9904bf59a11e8dcfe883e2ec32"
EXPECTED_ANALYTIC = "f2e8a7208e2522be69236de95c5b3de824feedab98ba4e840ca89d223349cc9d"
EXPECTED_METHOD = "948f3be5df0f07c8aa0c2a9fca8cf543dc46db3b39bda53eadf7427e014cfae4"
EXPECTED_PREVIOUS_FREEZE = "670afabf71f90d5bd9858ab3624dff8afa1666aa8507c3693126803e7e6bf16e"
EXPECTED_FAILED_FREEZE = "a566327eb68fbcba13e1dcec5f3ffbb0defb17b67f01f705765091571e540f12"
EXPECTED_PRIOR_RECEIPT = "8680b247968a37f01f57f3f3487f1fffdb021f0c62a9262ff7c8d928d38dffd2"
EXPECTED_CASE_RECEIPT = "3089d76acee64e3f64f23891922034cfc8a394e60606cd58037c7c73ad5e607b"
EXPECTED_ADAPTER_SNAPSHOT = "06f1a1d3e64edb3f950c260f673999048ffa2584025ee77fd922af99a09f92b5"
EXPECTED_PRECANON_FREEZE = "ef7a781d1fc0851f531c1cf61286a5c2786f908b997a9e7b31f10b2c5d581e1b"
EXPECTED_FAILED_AUDIT_RECEIPT = "bc4f59a87836128b47d210ddb35ac78d259d35a9c04f80a222560f071891149d"
EXPECTED_FAILED_AUDIT_STDERR = "4d41cc4b37ad49bbc2a872022acff8bbe4d7e3c814f665f44c038321a931eef5"
EXPECTED_REV2_VALIDATION = "1c9c0326d44c1dd0cf1f5b763ba00601ca699e57b8e14199d641dc8a91a058b4"
EXPECTED_CANONICALIZATION_RECEIPT = "13cee0096e73234c87cb0af9886647ba4ac0fec4b1cce001bd3074febc28b29f"
EXPECTED_FAILED_REV3_VALIDATION = "62756fa1a6d6f373ec4ffea8e56dfb4b2ef39a046900f4016b46d1dbb656a1b8"
EXPECTED_CURRENT_FREEZE = "90a3ec268983453477f7a18a8d12b5e3323ead007bd4b0831bf60ad0d4635e50"
EXPECTED_RENDERFIX_FREEZE_RECEIPT = "e4d2957d27584c7f41f6d66653ad0552f68918bc6aaf3b0059a2eee2205f64f6"
EXPECTED_REV4_VALIDATION = "e6f2be2c80dcbeeffd09fd3b7df520ea643b77a8dc4bcbc8089a882dd4acea59"
EXPECTED_FRESH_FAILED_RECEIPT = "589f4a7d0b35d434cae5ef3b43244b60aabd1b1f1d174f9ed29a9cca01a3202c"
EXPECTED_FRESH_FAILED_STDOUT = "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
EXPECTED_FRESH_FAILED_STDERR = "94c436a6fed7957daaf7d640b22556a3f52f1b2f7db16bf153c8518aad943ae3"
EXPECTED_PARTIAL_CELL_CSV = "8517032040d22b1a8d22b3f331431360b14c9c630675dbd0f5a8dbe8a9a469c5"
EXPECTED_PARTIAL_POINT_CSV = "0c76db6e93a4c1f31449aecae7cf99f8ca32a894fafdc9101b11fb87899d3d74"

sys.path.insert(0, str(REVIEW))
import pilot_actual_audit as adapter


def digest(path: Path) -> str | None:
    if not path.is_file():
        return None
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def sha_fields(value: object, path: str = "$") -> list[dict]:
    found: list[dict] = []
    if isinstance(value, dict):
        for key, item in value.items():
            child_path = f"{path}.{key}"
            if key == "sha256":
                normalized = item.lower() if isinstance(item, str) else None
                found.append({
                    "json_path": child_path,
                    "raw_value": item,
                    "normalized_lowercase": normalized,
                    "valid_64_hex_after_normalization": (
                        isinstance(normalized, str)
                        and re.fullmatch(r"[0-9a-f]{64}", normalized) is not None
                    ),
                })
            else:
                found.extend(sha_fields(item, child_path))
    elif isinstance(value, list):
        for index, item in enumerate(value):
            found.extend(sha_fields(item, f"{path}[{index}]"))
    return found


freeze = json.loads(FREEZE_PATH.read_text(encoding="utf-8"))
stage = WORKSPACE / freeze["input_scope"]["retained_stage"]
summary = WORKSPACE / freeze["input_scope"]["supervisor_summary"]
stage_expected = freeze["retained_stage_files"]
stage_actual = {}
if stage.is_dir():
    for path in sorted(stage.iterdir(), key=lambda item: item.name):
        if path.is_file():
            stage_actual[path.name] = {
                "bytes": path.stat().st_size,
                "sha256": digest(path),
            }

pin_rows = sha_fields(freeze)
all_pins_valid = bool(pin_rows) and all(
    row["valid_64_hex_after_normalization"] for row in pin_rows
)
all_pins_stored_lowercase = bool(pin_rows) and all(
    row["raw_value"] == row["normalized_lowercase"] for row in pin_rows
)
stage_checks = []
for name in sorted(stage_expected):
    expected = stage_expected[name]
    actual = stage_actual.get(name)
    stage_checks.append({
        "name": name,
        "expected_bytes": expected["bytes"],
        "actual_bytes": actual["bytes"] if actual else None,
        "expected_sha256_normalized": str(expected["sha256"]).lower(),
        "actual_sha256_normalized": str(actual["sha256"]).lower() if actual else None,
        "matches": bool(
            actual is not None
            and actual["bytes"] == expected["bytes"]
            and actual["sha256"].lower() == str(expected["sha256"]).lower()
        ),
    })

source_rows = {
    "adapter": {"path": REVIEW / "pilot_actual_audit.py", "expected": EXPECTED_ADAPTER},
    "analytic_source": {"path": REVIEW / "analytic_roi_reference.py", "expected": EXPECTED_ANALYTIC},
    "method_review": {"path": REVIEW / "PILOT_AUDIT_ADAPTER_REVIEW.md", "expected": EXPECTED_METHOD},
    "previous_source_freeze": {"path": REVIEW / "PILOT_AUDIT_INPUT_FREEZE_before_centroid_gate.json", "expected": EXPECTED_PREVIOUS_FREEZE},
    "preserved_failed_freeze": {"path": OLD_FREEZE, "expected": EXPECTED_FAILED_FREEZE},
    "preserved_failed_receipt": {"path": PRIOR_RECEIPT, "expected": EXPECTED_PRIOR_RECEIPT},
    "preserved_case_normalization_receipt": {"path": PRIOR_CASE_RECEIPT, "expected": EXPECTED_CASE_RECEIPT},
    "preserved_adapter_snapshot": {
        "path": REVIEW / "pilot_actual_audit_before_hashcase_fix_06F1A1D3E64EDB3F950C260F673999048FFA2584025EE77FD922AF99A09F92B5.py",
        "expected": EXPECTED_ADAPTER_SNAPSHOT,
    },
    "preserved_precanonicalization_freeze": {
        "path": REVIEW / "PILOT_AUDIT_INPUT_FREEZE_before_sha_lowercase_canonicalization.json",
        "expected": EXPECTED_PRECANON_FREEZE,
    },
    "preserved_failed_audit_receipt": {
        "path": REVIEW / "READONLY_AUDIT_EXECUTION_RECEIPT_corrected_freeze.json",
        "expected": EXPECTED_FAILED_AUDIT_RECEIPT,
    },
    "preserved_failed_audit_stderr": {
        "path": REVIEW / "READONLY_AUDIT_STDERR_corrected_freeze.txt",
        "expected": EXPECTED_FAILED_AUDIT_STDERR,
    },
    "preserved_rev2_validation": {
        "path": REVIEW / "PILOT_AUDIT_FREEZE_VALIDATION_rev2.json",
        "expected": EXPECTED_REV2_VALIDATION,
    },
    "canonicalization_receipt": {
        "path": REVIEW / "PILOT_AUDIT_FREEZE_CANONICALIZATION_rev3.json",
        "expected": EXPECTED_CANONICALIZATION_RECEIPT,
    },
    "preserved_failed_rev3_validation": {
        "path": REVIEW / "PILOT_AUDIT_FREEZE_VALIDATION_rev3.json",
        "expected": EXPECTED_FAILED_REV3_VALIDATION,
    },
    "preserved_previous_adapter_snapshot": {
        "path": REVIEW / "pilot_actual_audit_pre_retry_refactor_827F429C141403A5E85C7C9A67376D38B0628F9904BF59A11E8DCFE883E2EC32.py",
        "expected": EXPECTED_ADAPTER,
    },
    "preserved_pre_render_freeze": {
        "path": REVIEW / "PILOT_AUDIT_INPUT_FREEZE_before_plotcompat_20261005.json",
        "expected": "52ed7e2e29f3024af461f2c938e164cb04dec4878d1d796537ef965248df3d60",
    },
    "preserved_renderfix_freeze_receipt": {
        "path": REVIEW / "PILOT_AUDIT_FREEZE_RENDERFIX_rev4.json",
        "expected": EXPECTED_RENDERFIX_FREEZE_RECEIPT,
    },
    "preserved_rev4_validation": {
        "path": REVIEW / "PILOT_AUDIT_FREEZE_VALIDATION_rev4.json",
        "expected": EXPECTED_REV4_VALIDATION,
    },
    "preserved_previous_failed_audit_receipt": {
        "path": REVIEW / "READONLY_AUDIT_EXECUTION_RECEIPT_canonicalized_freeze.json",
        "expected": EXPECTED_FRESH_FAILED_RECEIPT,
    },
    "preserved_previous_failed_audit_stdout": {
        "path": REVIEW / "READONLY_AUDIT_STDOUT_canonicalized_freeze.txt",
        "expected": EXPECTED_FRESH_FAILED_STDOUT,
    },
    "preserved_previous_failed_audit_stderr": {
        "path": REVIEW / "READONLY_AUDIT_STDERR_canonicalized_freeze.txt",
        "expected": EXPECTED_FRESH_FAILED_STDERR,
    },
    "preserved_previous_partial_cell_csv": {
        "path": REVIEW / "PILOT_CELL_AUDIT.csv",
        "expected": EXPECTED_PARTIAL_CELL_CSV,
    },
    "preserved_previous_partial_point_csv": {
        "path": REVIEW / "PILOT_POINT_AUDIT.csv",
        "expected": EXPECTED_PARTIAL_POINT_CSV,
    },
}
source_checks = {}
for name, item in source_rows.items():
    actual_hash = digest(item["path"])
    source_checks[name] = {
        "path": str(item["path"]),
        "expected_sha256": item["expected"],
        "actual_sha256": actual_hash,
        "matches": actual_hash is not None and actual_hash.lower() == item["expected"],
    }

summary_pin = freeze["supervisor_summary"]
summary_check = {
    "bytes": summary.stat().st_size if summary.is_file() else None,
    "sha256": digest(summary),
    "expected_bytes": summary_pin["bytes"],
    "expected_sha256_normalized": str(summary_pin["sha256"]).lower(),
}
summary_check["matches"] = bool(
    summary_check["bytes"] == summary_check["expected_bytes"]
    and summary_check["sha256"] is not None
    and summary_check["sha256"].lower() == summary_check["expected_sha256_normalized"]
)

actual_names_match = sorted(stage_actual) == sorted(stage_expected)
all_stage_match = actual_names_match and all(row["matches"] for row in stage_checks)
source_match = all(row["matches"] for row in source_checks.values())
freeze_hash = digest(FREEZE_PATH)
corrected_pin = str(stage_expected["cell_ledger.partial.json"]["sha256"])
summary_stage_rows = json.loads(summary.read_text(encoding="utf-8")).get(
    "retained_stage", {}
).get("files", []) if summary.is_file() else []
summary_stage_hashes = {
    str(row.get("path")): str(row.get("sha256", "")).lower()
    for row in summary_stage_rows
}
freeze_stage_hashes = {
    str(name): str(row["sha256"]).lower()
    for name, row in stage_expected.items()
}
summary_inventory_matches = summary_stage_hashes == freeze_stage_hashes
adapter.validate_sha256_pins(freeze)
adapter_input_verification = adapter.verify_inputs(freeze)
adapter_normalizes_uppercase = (
    adapter.normalize_sha256_pin("A" * 64, "preflight uppercase probe")
    == "a" * 64
)
adapter_hash_path_passes = (
    adapter_input_verification["all_frozen_inputs_match"]
    and adapter_normalizes_uppercase
)
passed = all(
    [all_pins_valid, all_pins_stored_lowercase, actual_names_match, all_stage_match,
     summary_check["matches"], source_match, summary_inventory_matches,
     adapter_hash_path_passes]
)
current_freeze_matches_expected = freeze_hash == EXPECTED_CURRENT_FREEZE
passed = passed and current_freeze_matches_expected
receipt = {
    "schema": "pgfrag.forge_volume_ingress_corrected_freeze_validation.v1",
    "status": "PREFLIGHT_PASS_NOT_EXECUTED" if passed else "PREFLIGHT_FAIL_NOT_EXECUTED",
    "audit_invoked": False,
    "freeze_path": str(FREEZE_PATH),
    "freeze_sha256_before_execution": freeze_hash,
    "freeze_sha256_matches_expected": current_freeze_matches_expected,
    "corrected_cell_ledger_partial_sha256": corrected_pin,
    "all_sha256_pins_64_hex_after_lowercase_normalization": all_pins_valid,
    "all_sha256_pins_stored_lowercase": all_pins_stored_lowercase,
    "sha256_pin_count": len(pin_rows),
    "sha256_pins": pin_rows,
    "stage_directory_inventory_matches": actual_names_match,
    "stage_file_checks": stage_checks,
    "summary_check": summary_check,
    "summary_stage_inventory_matches_freeze": summary_inventory_matches,
    "adapter_verify_inputs_case_normalization_path": {
        "all_frozen_inputs_match": adapter_input_verification["all_frozen_inputs_match"],
        "stage_file_names_match_freeze": adapter_input_verification["stage_file_names_match_freeze"],
        "summary_matches_freeze": adapter_input_verification["summary_matches_freeze"],
        "uppercase_sha_probe_normalizes": adapter_normalizes_uppercase,
    },
    "source_and_preserved_artifact_checks": source_checks,
    "all_preflight_checks_pass": passed,
}
with OUTPUT.open("w", encoding="utf-8", newline="\n") as stream:
    json.dump(receipt, stream, indent=2, sort_keys=True, allow_nan=False)
    stream.write("\n")
print(json.dumps({
    "status": receipt["status"],
    "freeze_sha256": freeze_hash,
    "pin_count": len(pin_rows),
    "all_pins_valid": all_pins_valid,
    "all_pins_stored_lowercase": all_pins_stored_lowercase,
    "stage_inventory_match": actual_names_match,
    "all_stage_match": all_stage_match,
    "summary_match": summary_check["matches"],
    "summary_inventory_match": summary_inventory_matches,
    "source_match": source_match,
    "adapter_hash_path_passes": adapter_hash_path_passes,
    "corrected_partial_hash": corrected_pin,
    "receipt": OUTPUT.name,
}, sort_keys=True))
if not passed:
    raise SystemExit(2)
