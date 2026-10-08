"""Read-only mapping check for the one preserved local-front staging artifact."""

from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
STAGE = ROOT / "pgfrag/tmp/review_20261006/local_front_operator_20261006/runs/run01/.local_front_scalar_operator.json.149216.worker.json"
STAGE_SHA256 = "ebdea44857ab2b1a358fc679d0b1e4996943905a5996ba83f8f8f1ba68340762"
MAX_STAGE_BYTES = 1_000_000
EXPECTED = {"base": 100, "front_f1_shifted": 30, "heaviside": 10, "front_f1_unshifted": 0}


def _skip_space(text: str, offset: int) -> int:
    while offset < len(text) and text[offset].isspace():
        offset += 1
    return offset


def _read_top_level_registry(text: str) -> dict[str, object]:
    decoder = json.JSONDecoder()
    offset = _skip_space(text, 0)
    if offset >= len(text) or text[offset] != "{":
        raise ValueError("stage JSON root is not an object")
    offset += 1
    while True:
        offset = _skip_space(text, offset)
        if offset >= len(text) or text[offset] == "}":
            break
        key, offset = decoder.raw_decode(text, offset)
        offset = _skip_space(text, offset)
        if not isinstance(key, str) or offset >= len(text) or text[offset] != ":":
            raise ValueError("malformed top-level JSON key/value boundary")
        offset = _skip_space(text, offset + 1)
        if key == "registry":
            registry, _ = decoder.raw_decode(text, offset)
            if not isinstance(registry, dict):
                raise ValueError("registry entry is not an object")
            return registry
        _, offset = decoder.raw_decode(text, offset)
        offset = _skip_space(text, offset)
        if offset < len(text) and text[offset] == ",":
            offset += 1
        elif offset >= len(text) or text[offset] != "}":
            raise ValueError("malformed top-level JSON separator")
    raise ValueError("stage has no top-level registry object")


def inspect_saved_stage() -> dict[str, object]:
    if not STAGE.is_file() or STAGE.stat().st_size > MAX_STAGE_BYTES:
        raise ValueError("pinned staging artifact is missing or exceeds the 1 MB read cap")
    raw = STAGE.read_bytes()
    observed_hash = hashlib.sha256(raw).hexdigest()
    if observed_hash != STAGE_SHA256:
        raise ValueError(f"staging SHA mismatch: {observed_hash}")
    registry = _read_top_level_registry(raw.decode("utf-8"))
    ids, keys = registry.get("canonical_ids"), registry.get("ordered_semantic_keys")
    if not isinstance(ids, list) or not isinstance(keys, list) or len(ids) != 140 or len(keys) != 140:
        raise ValueError("canonical IDs and ordered semantic keys must each have 140 rows")
    if any(not isinstance(value, str) for value in ids) or len(set(ids)) != 140:
        raise ValueError("canonical IDs are malformed or duplicated")
    derived = Counter()
    for slot, row in enumerate(keys):
        if not isinstance(row, dict) or row.get("slot") != slot or row.get("canonical_id") != ids[slot]:
            raise ValueError(f"canonical semantic ordering mismatch at slot {slot}")
        kind = row.get("basis_kind")
        if kind not in EXPECTED:
            raise ValueError(f"unknown basis_kind wire value at slot {slot}: {kind!r}")
        derived[kind] += 1
    counts = {kind: derived[kind] for kind in EXPECTED}
    serialized = registry.get("basis_key_counts")
    if counts != EXPECTED or serialized != EXPECTED:
        raise ValueError(f"basis-kind counts differ from the explicit lowercase wire map: derived={counts}, serialized={serialized}")
    return {
        "schema": "local_front_basis_count_mapping_check.v1",
        "status": "PASS_MAPPING_ONLY",
        "stage_sha256": observed_hash,
        "canonical_dof_count": len(keys),
        "derived_basis_kind_counts": counts,
        "serialized_basis_key_counts": serialized,
        "accepted_unused_zero_kind": "front_f1_unshifted",
        "unknown_or_extra_basis_kinds_rejected": True,
        "provider_or_operator_invoked": False,
    }


if __name__ == "__main__":
    print(json.dumps(inspect_saved_stage(), indent=2, ensure_ascii=False, allow_nan=False))
