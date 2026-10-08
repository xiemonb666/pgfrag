"""Canonicalize frozen SHA-256 pin strings after source-only hash-case fix."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re


WORKSPACE = Path.cwd().resolve()
REVIEW = WORKSPACE / "pgfrag/tmp/review_20261005/forge_volume_ingress_independent_review"
FREEZE = REVIEW / "PILOT_AUDIT_INPUT_FREEZE.json"
ADAPTER = REVIEW / "pilot_actual_audit.py"
PRIOR_FREEZE = REVIEW / "PILOT_AUDIT_INPUT_FREEZE_before_sha_lowercase_canonicalization.json"
RECEIPT = REVIEW / "PILOT_AUDIT_FREEZE_CANONICALIZATION_rev3.json"
EXPECTED_PRIOR_FREEZE = "ef7a781d1fc0851f531c1cf61286a5c2786f908b997a9e7b31f10b2c5d581e1b"


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def collect_pins(value: object, path: str = "$") -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    if isinstance(value, dict):
        for key, item in value.items():
            child = f"{path}.{key}"
            if key == "sha256":
                if not isinstance(item, str) or re.fullmatch(r"[0-9A-Fa-f]{64}", item) is None:
                    raise ValueError(f"invalid source SHA-256 pin before normalization: {child}")
                rows.append((child, item))
            else:
                rows.extend(collect_pins(item, child))
    elif isinstance(value, list):
        for index, item in enumerate(value):
            rows.extend(collect_pins(item, f"{path}[{index}]"))
    return rows


if RECEIPT.exists():
    raise FileExistsError(f"refusing to overwrite {RECEIPT.name}")
prior_bytes = FREEZE.read_bytes()
prior_hash = hashlib.sha256(prior_bytes).hexdigest()
if prior_hash != EXPECTED_PRIOR_FREEZE:
    raise RuntimeError(f"input freeze changed unexpectedly: {prior_hash}")
if digest(PRIOR_FREEZE) != EXPECTED_PRIOR_FREEZE:
    raise RuntimeError("preserved pre-canonicalization freeze no longer matches its pin")

freeze = json.loads(prior_bytes.decode("utf-8"))
pins_before = collect_pins(freeze)
new_adapter_hash = digest(ADAPTER)
freeze["audit_adapter"]["sha256"] = new_adapter_hash
normalized_paths = []


def normalize(value: object, path: str = "$") -> None:
    if isinstance(value, dict):
        for key, item in value.items():
            child = f"{path}.{key}"
            if key == "sha256":
                if not isinstance(item, str) or re.fullmatch(r"[0-9A-Fa-f]{64}", item) is None:
                    raise ValueError(f"invalid source SHA-256 pin: {child}")
                lowered = item.lower()
                if lowered != item:
                    normalized_paths.append(child)
                value[key] = lowered
            else:
                normalize(item, child)
    elif isinstance(value, list):
        for index, item in enumerate(value):
            normalize(item, f"{path}[{index}]")


normalize(freeze)
pins_after = collect_pins(freeze)
all_lowercase = all(value == value.lower() for _, value in pins_after)
all_valid_64_hex = all(re.fullmatch(r"[0-9a-f]{64}", value) for _, value in pins_after)
pin_count_stable = len(pins_before) == len(pins_after)
if not (all_lowercase and all_valid_64_hex and pin_count_stable):
    raise RuntimeError("SHA-256 pin canonicalization validation failed")

with FREEZE.open("w", encoding="utf-8", newline="\n") as stream:
    json.dump(freeze, stream, indent=2, sort_keys=False, allow_nan=False)
    stream.write("\n")

new_freeze_hash = digest(FREEZE)
partial_hash = freeze["retained_stage_files"]["cell_ledger.partial.json"]["sha256"]
receipt = {
    "schema": "pgfrag.forge_volume_ingress_freeze_pin_canonicalization.v1",
    "status": "PIN_CANONICALIZATION_COMPLETE_NOT_EXECUTED",
    "audit_invoked": False,
    "prior_freeze_sha256": prior_hash,
    "new_freeze_sha256": new_freeze_hash,
    "new_adapter_sha256": new_adapter_hash,
    "normalized_pin_count": len(pins_after),
    "all_pins_64_hex_lowercase": all_lowercase and all_valid_64_hex,
    "normalized_paths": normalized_paths,
    "corrected_partial_ledger_pin": partial_hash,
    "prior_freeze_copy_sha256": digest(PRIOR_FREEZE),
}
with RECEIPT.open("w", encoding="utf-8", newline="\n") as stream:
    json.dump(receipt, stream, indent=2, sort_keys=True, allow_nan=False)
    stream.write("\n")
print(json.dumps(receipt, sort_keys=True))
