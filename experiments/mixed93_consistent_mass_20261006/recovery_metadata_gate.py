"""Standard-library-only gate for the failed producer's independent recovery JSON."""

import hashlib
import json
import math
from pathlib import Path


def _sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _at(record, dotted_path):
    value = record
    for part in dotted_path.split("."):
        value = value[part]
    return value


def _strings(value):
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for child in value.values():
            yield from _strings(child)
    elif isinstance(value, list):
        for child in value:
            yield from _strings(child)


def _resolved(path_text, root):
    path = Path(path_text)
    return (path if path.is_absolute() else root / path).resolve()


def verify_recovery_contract(case, workspace_root, arrays_path, arrays_sha256):
    """Read and validate only JSON metadata; call before importing NumPy."""
    binding = case["independent_recovery_contract"]
    inputs = case["inputs"]
    expected = binding
    if not binding.get("json_path") or not binding.get("json_sha256"):
        raise RuntimeError("BLOCKED: independent recovery JSON pins are not root-bound")
    if not inputs.get("arrays_path") or not inputs.get("arrays_sha256"):
        raise RuntimeError("BLOCKED: fixed failed arrays path and SHA-256 are not root-bound")
    if not inputs.get("result_path") or not inputs.get("result_sha256"):
        raise RuntimeError("BLOCKED: independent recovery JSON alias path and SHA-256 are not root-bound")
    required_values = ("expected_schema", "expected_status", "expected_saved_cache_state",
                       "expected_physical_fracture_status", "expected_dof_count", "expected_kkt_order",
                       "expected_source_count", "expected_global_branch_count", "expected_producer_return_code",
                       "expected_result_json_present", "expected_failed_arrays_bytes")
    if any(expected.get(key) is None for key in required_values):
        raise RuntimeError("BLOCKED: recovery contract values still require root binding")

    recovery_path = _resolved(binding["json_path"], workspace_root)
    if (_resolved(inputs["result_path"], workspace_root) != recovery_path
            or inputs["result_sha256"].lower() != binding["json_sha256"].lower()):
        raise RuntimeError("BLOCKED: case result alias must identify the pinned independent recovery JSON")
    if not recovery_path.is_file() or _sha256(recovery_path).lower() != binding["json_sha256"].lower():
        raise RuntimeError("BLOCKED: independent recovery JSON is absent or its SHA-256 changed")
    recovery = json.loads(recovery_path.read_text(encoding="utf-8"))
    fields = binding["field_paths"]
    get = lambda name: _at(recovery, fields[name])
    exact = (("schema", "expected_schema"), ("status", "expected_status"),
             ("saved_cache_state", "expected_saved_cache_state"),
             ("physical_fracture_status", "expected_physical_fracture_status"),
             ("dof_count", "expected_dof_count"), ("kkt_order", "expected_kkt_order"),
             ("source_count", "expected_source_count"),
             ("global_branch_count", "expected_global_branch_count"),
             ("producer_return_code", "expected_producer_return_code"),
             ("producer_result_json_present", "expected_result_json_present"))
    for field_name, expected_name in exact:
        value, target = get(field_name), expected[expected_name]
        if type(value) is not type(target) or value != target:
            raise RuntimeError(f"BLOCKED: recovery field {field_name} does not match the fixed contract")

    declared_arrays = _resolved(get("failed_arrays_path"), workspace_root)
    if declared_arrays != arrays_path.resolve() or declared_arrays != _resolved(inputs["arrays_path"], workspace_root):
        raise RuntimeError("BLOCKED: recovery NPZ path does not match the fixed mass input")
    array_pins = (get("failed_arrays_sha256_pre"), get("failed_arrays_sha256_post"))
    if any(not isinstance(pin, str) or pin.lower() != arrays_sha256.lower() for pin in array_pins):
        raise RuntimeError("BLOCKED: recovery NPZ pre/post hashes do not match the fixed mass input")
    byte_count = get("failed_arrays_bytes")
    if type(byte_count) is not int or byte_count != expected["expected_failed_arrays_bytes"] or byte_count != arrays_path.stat().st_size:
        raise RuntimeError("BLOCKED: recovery NPZ byte count does not match the pinned file")

    result_path = _resolved(get("producer_result_json_path"), workspace_root)
    if result_path != arrays_path.with_name("result.json").resolve() or result_path.exists():
        raise RuntimeError("BLOCKED: failed producer result JSON path is not the absent sibling result")
    producer = recovery["producer_run"]
    elapsed = producer.get("elapsed_seconds")
    if isinstance(elapsed, bool) or not isinstance(elapsed, (float, int)) or not math.isfinite(float(elapsed)) or elapsed <= 0:
        raise RuntimeError("BLOCKED: producer elapsed time is missing or invalid")
    failure_text = " ".join(_strings(producer)).lower()
    if (any(term.lower() not in failure_text for term in binding["producer_failure_text_terms"])
            or "serial" not in failure_text):
        raise RuntimeError("BLOCKED: recovery metadata does not retain the producer JSON-serialization failure")

    return {
        "recovery_json_path": str(recovery_path), "recovery_json_sha256": _sha256(recovery_path),
        "case_result_alias_path": str(recovery_path), "case_result_alias_sha256": binding["json_sha256"].lower(),
        "saved_cache_recovery_status": get("status"), "saved_cache_state": get("saved_cache_state"),
        "physical_fracture_status": get("physical_fracture_status"),
        "producer_return_code": get("producer_return_code"),
        "producer_elapsed_seconds": float(elapsed), "producer_result_json_present": get("producer_result_json_present"),
        "producer_result_json_path": str(result_path), "producer_failure_details": failure_text,
        "failed_arrays_path": str(declared_arrays), "failed_arrays_sha256": arrays_sha256.lower(),
        "failed_arrays_bytes": byte_count, "dof_count": get("dof_count"), "kkt_order": get("kkt_order"),
        "source_count": get("source_count"), "global_branch_count": get("global_branch_count"),
    }
