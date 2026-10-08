"""Pin the reviewed render-safe adapter into a new lowercase freeze revision."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re


WORKSPACE = Path.cwd().resolve()
REVIEW = WORKSPACE / "pgfrag/tmp/review_20261005/forge_volume_ingress_independent_review"
FREEZE = REVIEW / "PILOT_AUDIT_INPUT_FREEZE.json"
PRIOR_FREEZE = REVIEW / "PILOT_AUDIT_INPUT_FREEZE_before_plotcompat_20261005.json"
ADAPTER = REVIEW / "pilot_actual_audit.py"
SOURCE_SNAPSHOT = REVIEW / "pilot_actual_audit_pre_retry_refactor_827F429C141403A5E85C7C9A67376D38B0628F9904BF59A11E8DCFE883E2EC32.py"
RECEIPT = REVIEW / "PILOT_AUDIT_FREEZE_RENDERFIX_rev4.json"
EXPECTED_PRIOR_FREEZE = "52ed7e2e29f3024af461f2c938e164cb04dec4878d1d796537ef965248df3d60"
EXPECTED_ADAPTER = "827f429c141403a5e85c7c9a67376d38b0628f9904bf59a11e8dcfe883e2ec32"


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def pins(value: object, path: str = "$") -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    if isinstance(value, dict):
        for key, item in value.items():
            child = f"{path}.{key}"
            if key == "sha256":
                if not isinstance(item, str) or re.fullmatch(r"[0-9a-f]{64}", item) is None:
                    raise ValueError(f"noncanonical frozen SHA-256 pin: {child}")
                rows.append((child, item))
            else:
                rows.extend(pins(item, child))
    elif isinstance(value, list):
        for index, item in enumerate(value):
            rows.extend(pins(item, f"{path}[{index}]"))
    return rows


if RECEIPT.exists():
    raise FileExistsError(f"refusing to overwrite {RECEIPT.name}")
if digest(FREEZE) != EXPECTED_PRIOR_FREEZE:
    raise RuntimeError("current freeze no longer matches the preserved pre-render revision")
if digest(PRIOR_FREEZE) != EXPECTED_PRIOR_FREEZE:
    raise RuntimeError("preserved pre-render freeze hash mismatch")
if digest(SOURCE_SNAPSHOT) != EXPECTED_ADAPTER:
    raise RuntimeError("preserved adapter snapshot hash mismatch")
new_adapter_hash = digest(ADAPTER)
freeze = json.loads(FREEZE.read_text(encoding="utf-8"))
freeze["audit_adapter"]["sha256"] = new_adapter_hash
all_pins = pins(freeze)
if not all(value == value.lower() for _, value in all_pins):
    raise RuntimeError("freeze SHA-256 pins are not all lowercase")

with FREEZE.open("w", encoding="utf-8", newline="\n") as stream:
    json.dump(freeze, stream, indent=2, sort_keys=False, allow_nan=False)
    stream.write("\n")
new_freeze_hash = digest(FREEZE)
receipt = {
    "schema": "pgfrag.forge_volume_ingress_renderretry_freeze_update.v1",
    "status": "NEW_FREEZE_PINNED_NOT_EXECUTED",
    "audit_invoked": False,
    "prior_freeze_sha256": EXPECTED_PRIOR_FREEZE,
    "new_freeze_sha256": new_freeze_hash,
    "adapter_sha256": new_adapter_hash,
    "adapter_sha256_matches_preserved_source_snapshot": new_adapter_hash == EXPECTED_ADAPTER,
    "all_sha256_pins_count": len(all_pins),
    "all_sha256_pins_64_hex_lowercase": all(
        re.fullmatch(r"[0-9a-f]{64}", value) is not None for _, value in all_pins
    ),
    "freeze_sha256_pins": dict(all_pins),
    "preserved_pre_render_freeze_sha256": digest(PRIOR_FREEZE),
    "preserved_adapter_snapshot_sha256": digest(SOURCE_SNAPSHOT),
}
with RECEIPT.open("w", encoding="utf-8", newline="\n") as stream:
    json.dump(receipt, stream, indent=2, sort_keys=True, allow_nan=False)
    stream.write("\n")
print(json.dumps(receipt, sort_keys=True))
