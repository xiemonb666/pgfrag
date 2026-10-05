"""Compile the strict JSON audit mask into a bounded, typed TOML transport.

This is a transport-only helper. It never calls Peridynamics, changes native state,
or authorizes initialization. The JSON audit contract remains the source artifact;
the generated TOML includes the SHA-256 of those exact JSON bytes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import stat
import tempfile
import tomllib
import unicodedata
from pathlib import Path
from typing import Any, Mapping, Sequence


AUDIT_SCHEMA = "pgfrag.material-visibility-native-mask.v1"
TRANSPORT_SCHEMA = "pgfrag.native-material-visibility-toml.v1"
STATES = ("IN_MATERIAL", "OUT_OR_VOID", "UNKNOWN")
ELIGIBILITIES = {
    "BLOCKED_UNKNOWN",
    "BLOCKED_EMPTY_CANDIDATE_SET",
    "BLOCKED_EMPTY_MASK",
    "BLOCKED_ZERO_DEGREE",
    "CANDIDATE_MASK_POSITIVE_DEGREE_NATIVE_UNVERIFIED",
}
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_KEY_RE = re.compile(r"^[A-Za-z0-9_-]+$")
MAX_AUDIT_JSON_BYTES = 64 * 1024 * 1024
MAX_NATIVE_TOML_BYTES = 128 * 1024 * 1024
MAX_POINTS = 5000
MAX_CANDIDATE_PAIRS = 250_000
MAX_EXTERNAL_ID = (1 << 63) - 1
MAX_REASON_UTF8_BYTES = 4096
MAX_PATH_UTF8_BYTES = 16_384


class TransportError(ValueError):
    """The audit contract cannot be safely represented by the native TOML transport."""


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _metadata(info: os.stat_result) -> tuple[int, int, int, int, int]:
    return (
        int(info.st_size),
        int(info.st_mtime_ns),
        int(getattr(info, "st_dev", 0)),
        int(getattr(info, "st_ino", 0)),
        int(stat.S_IFMT(info.st_mode)),
    )


def _read_stable_bounded(path: Path, label: str, limit: int) -> bytes:
    """Stat before opening, read twice with caps, and reject path/handle changes."""
    try:
        resolved = Path(path).resolve(strict=True)
        before = _metadata(resolved.stat())
    except OSError as exc:
        raise TransportError(f"cannot resolve/stat {label}: {path}") from exc
    if not stat.S_ISREG(resolved.stat().st_mode):
        raise TransportError(f"{label} must be a regular file")
    if before[0] < 0 or before[0] > limit:
        raise TransportError(f"{label} exceeds the {limit}-byte LocalTiny limit")

    def one_read() -> bytes:
        try:
            with resolved.open("rb") as stream:
                handle_before = _metadata(os.fstat(stream.fileno()))
                path_before = _metadata(resolved.stat())
                if handle_before != before or path_before != before:
                    raise TransportError(f"{label} changed after its size preflight")
                data = stream.read(before[0] + 1)
                if len(data) != before[0]:
                    raise TransportError(f"{label} changed size during bounded read")
                handle_after = _metadata(os.fstat(stream.fileno()))
            path_after = _metadata(resolved.stat())
        except OSError as exc:
            raise TransportError(f"cannot read {label}: {resolved}") from exc
        if handle_after != before or path_after != before:
            raise TransportError(f"{label} changed during bounded read")
        return data

    first = one_read()
    second = one_read()
    if first != second:
        raise TransportError(f"{label} changed between bounded reads")
    return first


def _reject_json_constant(value: str) -> None:
    raise TransportError(f"audit JSON contains non-standard number {value}")


def _reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise TransportError(f"audit JSON contains duplicate key {key!r}")
        result[key] = value
    return result


def _load_strict_json(data: bytes) -> Any:
    try:
        return json.loads(
            data.decode("utf-8", errors="strict"),
            parse_constant=_reject_json_constant,
            object_pairs_hook=_reject_duplicate_keys,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise TransportError("audit contract is not strict UTF-8 JSON") from exc


def _expect_object(value: Any, keys: set[str], label: str) -> Mapping[str, Any]:
    if not isinstance(value, dict) or set(value) != keys:
        raise TransportError(f"{label} must be an object with exactly these keys: {sorted(keys)}")
    return value


def _expect_string(value: Any, label: str, *, max_bytes: int = MAX_PATH_UTF8_BYTES) -> str:
    if not isinstance(value, str) or not value:
        raise TransportError(f"{label} must be a non-empty string")
    try:
        encoded = value.encode("utf-8", errors="strict")
    except UnicodeEncodeError as exc:
        raise TransportError(f"{label} contains an invalid Unicode scalar") from exc
    if len(encoded) > max_bytes:
        raise TransportError(f"{label} exceeds its {max_bytes}-byte string limit")
    return value


def _expect_sha(value: Any, label: str) -> str:
    if not isinstance(value, str) or _SHA256_RE.fullmatch(value) is None:
        raise TransportError(f"{label} must be a lowercase SHA-256 digest")
    return value


def _expect_int(value: Any, label: str, *, minimum: int = 0, maximum: int = MAX_EXTERNAL_ID) -> int:
    if type(value) is not int or value < minimum or value > maximum:
        raise TransportError(f"{label} must be an integer in [{minimum}, {maximum}]")
    return value


def _expect_float(value: Any, label: str, *, minimum: float | None = None, strict_minimum: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TransportError(f"{label} must be a finite number")
    try:
        converted = float(value)
    except (OverflowError, ValueError) as exc:
        raise TransportError(f"{label} is outside the supported numeric range") from exc
    if not math.isfinite(converted):
        raise TransportError(f"{label} must be finite")
    if minimum is not None and (converted <= minimum if strict_minimum else converted < minimum):
        operator = ">" if strict_minimum else ">="
        raise TransportError(f"{label} must be {operator} {minimum}")
    return converted


def _hash_point_order(point_ids: Sequence[int]) -> str:
    payload = bytearray(b"pgfrag.point-ids.v1\n")
    for point_id in point_ids:
        payload.extend(f"{point_id}\n".encode("ascii"))
    return sha256_bytes(bytes(payload))


def _hash_pair_map(rows: Sequence[Mapping[str, Any]]) -> str:
    payload = bytearray(b"pgfrag.native-pair-map.v1\n")
    for row in rows:
        left, right = row["pair"]
        payload.extend(f"{left},{right},{row['classification']}\n".encode("ascii"))
    return sha256_bytes(bytes(payload))


def _hash_included_mask(pairs: Sequence[Sequence[int]]) -> str:
    payload = bytearray(b"pgfrag.native-included-mask.v1\n")
    for left, right in pairs:
        payload.extend(f"{left},{right}\n".encode("ascii"))
    return sha256_bytes(bytes(payload))


def _validate_transform(value: Any) -> list[list[float]]:
    if not isinstance(value, list) or len(value) != 4:
        raise TransportError("native frame transform must be a 4x4 array")
    matrix: list[list[float]] = []
    for row_index, row in enumerate(value):
        if not isinstance(row, list) or len(row) != 4:
            raise TransportError("native frame transform must be a 4x4 array")
        matrix.append([_expect_float(item, f"transform[{row_index}]") for item in row])
    if matrix[3] != [0.0, 0.0, 0.0, 1.0]:
        raise TransportError("native frame transform must use homogeneous affine form")
    rotation = [row[:3] for row in matrix[:3]]
    for i in range(3):
        for j in range(3):
            dot = sum(rotation[k][i] * rotation[k][j] for k in range(3))
            expected = 1.0 if i == j else 0.0
            if not math.isclose(dot, expected, rel_tol=0.0, abs_tol=1e-12):
                raise TransportError("native frame transform rotation must be orthonormal")
    determinant = (
        rotation[0][0] * (rotation[1][1] * rotation[2][2] - rotation[1][2] * rotation[2][1])
        - rotation[0][1] * (rotation[1][0] * rotation[2][2] - rotation[1][2] * rotation[2][0])
        + rotation[0][2] * (rotation[1][0] * rotation[2][1] - rotation[1][1] * rotation[2][0])
    )
    if not math.isclose(determinant, 1.0, rel_tol=0.0, abs_tol=1e-12):
        raise TransportError("native frame transform must preserve orientation")
    return matrix


def _canonical_pair(value: Any, label: str, point_ids: set[int]) -> tuple[int, int]:
    if not isinstance(value, list) or len(value) != 2:
        raise TransportError(f"{label} must be a two-ID array")
    left = _expect_int(value[0], f"{label}[0]")
    right = _expect_int(value[1], f"{label}[1]")
    if left >= right or left not in point_ids or right not in point_ids:
        raise TransportError(f"{label} must be canonical i<j and reference known point IDs")
    return left, right


def _validate_contract(contract: Any, actual_json_sha256: str, audit_path: Path) -> dict[str, Any]:
    top_keys = {
        "schema", "static_contract_status", "native_mask_eligibility", "native_gate",
        "source", "coordinates", "points", "candidate_generation",
        "raw_euclidean_candidate_pairs", "counts", "pair_map_sha256", "mask",
        "evidence", "initial_physics_events",
    }
    root = _expect_object(contract, top_keys, "audit contract")
    if root["schema"] != AUDIT_SCHEMA or root["static_contract_status"] != "VALID":
        raise TransportError("unsupported or statically invalid audit contract")
    eligibility = root["native_mask_eligibility"]
    if not isinstance(eligibility, str) or eligibility not in ELIGIBILITIES:
        raise TransportError("unsupported native mask eligibility state")

    gate = _expect_object(root["native_gate"], {
        "status", "raw_candidate_set_reproduced", "native_graph_readback",
        "initialization_authorized", "candidate_pair_count_nonzero",
        "included_pair_count_nonzero", "zero_degree_point_ids",
    }, "native_gate")
    if (gate["status"] != "UNVERIFIED" or gate["raw_candidate_set_reproduced"] is not False
            or gate["native_graph_readback"] != "NOT_RUN"
            or gate["initialization_authorized"] is not False):
        raise TransportError("audit contract attempts to claim native verification or authorization")
    if type(gate["candidate_pair_count_nonzero"]) is not bool or type(gate["included_pair_count_nonzero"]) is not bool:
        raise TransportError("native_gate count flags must be booleans")
    if not isinstance(gate["zero_degree_point_ids"], list) or len(gate["zero_degree_point_ids"]) > MAX_POINTS:
        raise TransportError("native_gate.zero_degree_point_ids exceeds the point cap")
    reported_zero_degree_ids = [
        _expect_int(value, f"native_gate.zero_degree_point_ids[{index}]")
        for index, value in enumerate(gate["zero_degree_point_ids"])
    ]
    if len(reported_zero_degree_ids) != len(set(reported_zero_degree_ids)):
        raise TransportError("native_gate.zero_degree_point_ids must be unique")

    source = _expect_object(root["source"], {
        "source_cad_path", "source_cad_sha256", "classifier_source_sha256",
        "source_hash_status", "source_status", "source_accuracy_mm", "source_accuracy_basis",
    }, "source")
    source_cad_path = _expect_string(source["source_cad_path"], "source.source_cad_path")
    source_cad_sha = _expect_sha(source["source_cad_sha256"], "source.source_cad_sha256")
    classifier_source_sha = _expect_sha(source["classifier_source_sha256"], "source.classifier_source_sha256")
    if source["source_hash_status"] != "MATCHED_ACTUAL_FILE_BYTES" or source["source_status"] != "READY":
        raise TransportError("source must be SHA-bound to a READY canonical CAD file")
    source_accuracy = _expect_float(source["source_accuracy_mm"], "source.source_accuracy_mm", minimum=0.0)
    source_accuracy_basis = _expect_string(source["source_accuracy_basis"], "source.source_accuracy_basis", max_bytes=4096)

    coordinates = _expect_object(root["coordinates"], {
        "input_unit", "native_unit", "m_per_input_unit", "reference_frame", "native_frame",
        "native_frame_from_reference_frame_m",
    }, "coordinates")
    if (coordinates["input_unit"] != "mm" or coordinates["native_unit"] != "m"
            or coordinates["m_per_input_unit"] != 0.001):
        raise TransportError("coordinate contract must explicitly preserve mm-to-m conversion")
    reference_frame = _expect_string(coordinates["reference_frame"], "coordinates.reference_frame", max_bytes=512)
    native_frame = _expect_string(coordinates["native_frame"], "coordinates.native_frame", max_bytes=512)
    transform = _validate_transform(coordinates["native_frame_from_reference_frame_m"])

    points = _expect_object(root["points"], {
        "payload_path", "payload_sha256", "order", "point_ids_order_sha256", "rows",
    }, "points")
    point_path = _expect_string(points["payload_path"], "points.payload_path")
    point_payload_sha = _expect_sha(points["payload_sha256"], "points.payload_sha256")
    point_order_sha = _expect_sha(points["point_ids_order_sha256"], "points.point_ids_order_sha256")
    if points["order"] != "exact_payload_row_order":
        raise TransportError("point order must bind to exact payload row order")
    point_rows = points["rows"]
    if not isinstance(point_rows, list) or not point_rows or len(point_rows) > MAX_POINTS:
        raise TransportError(f"points.rows must contain 1..{MAX_POINTS} rows")
    validated_points: list[dict[str, Any]] = []
    point_ids: list[int] = []
    native_index_by_id: dict[int, int] = {}
    for row_number, raw_row in enumerate(point_rows, start=1):
        row = _expect_object(raw_row, {"point_id", "native_index", "position_mm"}, f"points.rows[{row_number - 1}]")
        point_id = _expect_int(row["point_id"], f"points.rows[{row_number - 1}].point_id")
        native_index = _expect_int(
            row["native_index"], f"points.rows[{row_number - 1}].native_index", minimum=1, maximum=MAX_POINTS
        )
        if native_index != row_number or point_id in native_index_by_id:
            raise TransportError("point native_index must be exactly 1..N in row order and IDs must be unique")
        position = row["position_mm"]
        if not isinstance(position, list) or len(position) != 3:
            raise TransportError(f"points.rows[{row_number - 1}].position_mm must have three coordinates")
        position_mm = [_expect_float(value, f"points.rows[{row_number - 1}].position_mm") for value in position]
        native_index_by_id[point_id] = native_index
        point_ids.append(point_id)
        validated_points.append({"point_id": point_id, "native_index": native_index, "position_mm": position_mm})
    if _hash_point_order(point_ids) != point_order_sha:
        raise TransportError("point_ids_order_sha256 does not match explicit row order")

    generation = _expect_object(root["candidate_generation"], {
        "metric", "horizon_m", "candidate_rule", "scope",
    }, "candidate_generation")
    if (generation["metric"] != "euclidean_reference_distance"
            or generation["candidate_rule"] != "exact_pinned_Peridynamics_find_bonds_convention"
            or generation["scope"] != "one_body_one_uniform_horizon"):
        raise TransportError("unsupported candidate-generation convention")
    horizon_m = _expect_float(generation["horizon_m"], "candidate_generation.horizon_m", minimum=0.0, strict_minimum=True)

    raw_candidates = root["raw_euclidean_candidate_pairs"]
    if not isinstance(raw_candidates, list) or len(raw_candidates) > MAX_CANDIDATE_PAIRS:
        raise TransportError(f"candidate ledger exceeds {MAX_CANDIDATE_PAIRS} pairs")
    point_id_set = set(point_ids)
    validated_candidates: list[dict[str, Any]] = []
    decisions_for_hash: list[dict[str, Any]] = []
    seen_pairs: set[tuple[int, int]] = set()
    previous_pair: tuple[int, int] | None = None
    for index, raw_row in enumerate(raw_candidates):
        row = _expect_object(raw_row, {"pair", "classification", "reason"}, f"raw candidate {index}")
        pair = _canonical_pair(row["pair"], f"raw candidate {index}.pair", point_id_set)
        if pair in seen_pairs or (previous_pair is not None and pair <= previous_pair):
            raise TransportError("raw candidate pairs must be unique and sorted lexicographically")
        seen_pairs.add(pair)
        previous_pair = pair
        state = row["classification"]
        if state not in STATES:
            raise TransportError(f"raw candidate {index} has unsupported classification")
        reason = _expect_string(row["reason"], f"raw candidate {index}.reason", max_bytes=MAX_REASON_UTF8_BYTES)
        decisions_for_hash.append({"pair": [pair[0], pair[1]], "classification": state})
        validated_candidates.append({
            "point_i": pair[0],
            "point_j": pair[1],
            "native_i": native_index_by_id[pair[0]],
            "native_j": native_index_by_id[pair[1]],
            "classification": state,
            "reason": reason,
        })
    if _hash_pair_map(decisions_for_hash) != _expect_sha(root["pair_map_sha256"], "pair_map_sha256"):
        raise TransportError("pair_map_sha256 does not match full sorted candidate ledger")

    counts = _expect_object(root["counts"], {
        "raw_euclidean_candidate_pairs", "initial_material_bond_pairs", "geometry_excluded_pairs", "unknown_pairs",
    }, "counts")
    state_counts = {state: sum(row["classification"] == state for row in decisions_for_hash) for state in STATES}
    expected_counts = {
        "raw_euclidean_candidate_pairs": len(decisions_for_hash),
        "initial_material_bond_pairs": state_counts["IN_MATERIAL"],
        "geometry_excluded_pairs": state_counts["OUT_OR_VOID"],
        "unknown_pairs": state_counts["UNKNOWN"],
    }
    for name, expected in expected_counts.items():
        if _expect_int(counts[name], f"counts.{name}", maximum=MAX_CANDIDATE_PAIRS) != expected:
            raise TransportError(f"counts.{name} disagrees with candidate ledger")

    mask = _expect_object(root["mask"], {
        "included_pairs", "geometry_excluded_pairs", "unknown_pairs", "mask_sha256",
    }, "mask")
    mask_groups = {
        "included_pairs": "IN_MATERIAL",
        "geometry_excluded_pairs": "OUT_OR_VOID",
        "unknown_pairs": "UNKNOWN",
    }
    native_mask_rows: dict[str, list[dict[str, int]]] = {}
    for field, state in mask_groups.items():
        raw_group = mask[field]
        if not isinstance(raw_group, list) or len(raw_group) > MAX_CANDIDATE_PAIRS:
            raise TransportError(f"mask.{field} exceeds the pair cap")
        expected_pairs = [tuple(row["pair"]) for row in decisions_for_hash if row["classification"] == state]
        parsed_pairs = [_canonical_pair(pair, f"mask.{field}[{i}]", point_id_set) for i, pair in enumerate(raw_group)]
        if parsed_pairs != expected_pairs:
            raise TransportError(f"mask.{field} disagrees with the full candidate state partition")
        native_mask_rows[field] = [
            {"point_i": left, "point_j": right,
             "native_i": native_index_by_id[left], "native_j": native_index_by_id[right]}
            for left, right in parsed_pairs
        ]
    included_pairs = [tuple(row["pair"]) for row in decisions_for_hash if row["classification"] == "IN_MATERIAL"]
    mask_sha = _expect_sha(mask["mask_sha256"], "mask.mask_sha256")
    if _hash_included_mask(included_pairs) != mask_sha:
        raise TransportError("mask.mask_sha256 does not match included canonical pairs")
    degrees = {point_id: 0 for point_id in point_ids}
    for left, right in included_pairs:
        degrees[left] += 1
        degrees[right] += 1
    zero_degree_ids = [point_id for point_id in point_ids if degrees[point_id] == 0]
    if (reported_zero_degree_ids != zero_degree_ids
            or gate["candidate_pair_count_nonzero"] is not bool(decisions_for_hash)
            or gate["included_pair_count_nonzero"] is not bool(included_pairs)):
        raise TransportError("native_gate static count flags disagree with contract arrays")
    if state_counts["UNKNOWN"]:
        expected_eligibility = "BLOCKED_UNKNOWN"
    elif not decisions_for_hash:
        expected_eligibility = "BLOCKED_EMPTY_CANDIDATE_SET"
    elif not included_pairs:
        expected_eligibility = "BLOCKED_EMPTY_MASK"
    elif zero_degree_ids:
        expected_eligibility = "BLOCKED_ZERO_DEGREE"
    else:
        expected_eligibility = "CANDIDATE_MASK_POSITIVE_DEGREE_NATIVE_UNVERIFIED"
    if eligibility != expected_eligibility:
        raise TransportError("native_mask_eligibility disagrees with unknown/empty/degree conditions")

    evidence = _expect_object(root["evidence"], {
        "classifier_result_path", "classifier_result_sha256", "classifier_source_path",
        "classifier_source_sha256", "raw_hit_interval_sidecar_path",
        "raw_hit_interval_sidecar_sha256", "sidecar_selection", "visibility_algorithm",
    }, "evidence")
    classifier_result_path = _expect_string(evidence["classifier_result_path"], "evidence.classifier_result_path")
    classifier_result_sha = _expect_sha(evidence["classifier_result_sha256"], "evidence.classifier_result_sha256")
    classifier_source_path = _expect_string(evidence["classifier_source_path"], "evidence.classifier_source_path")
    evidence_classifier_source_sha = _expect_sha(evidence["classifier_source_sha256"], "evidence.classifier_source_sha256")
    sidecar_path = _expect_string(evidence["raw_hit_interval_sidecar_path"], "evidence.raw_hit_interval_sidecar_path")
    sidecar_sha = _expect_sha(evidence["raw_hit_interval_sidecar_sha256"], "evidence.raw_hit_interval_sidecar_sha256")
    if evidence_classifier_source_sha != classifier_source_sha:
        raise TransportError("classifier source SHA differs between source and evidence")
    selection = _expect_object(evidence["sidecar_selection"], {"test_nodeid", "call_index"}, "sidecar_selection")
    test_nodeid = _expect_string(selection["test_nodeid"], "sidecar_selection.test_nodeid", max_bytes=4096)
    call_index = _expect_int(selection["call_index"], "sidecar_selection.call_index", maximum=(1 << 31) - 1)
    algorithm = _expect_string(evidence["visibility_algorithm"], "evidence.visibility_algorithm", max_bytes=512)

    physics = _expect_object(root["initial_physics_events"], {
        "evidence_status", "planned_physical_broken_bond_pairs", "geometry_exclusions_are_physics_events",
    }, "initial_physics_events")
    if (physics["evidence_status"] != "PLANNED_NOT_OBSERVED_NATIVE_NOT_RUN"
            or type(physics["planned_physical_broken_bond_pairs"]) is not int
            or physics["planned_physical_broken_bond_pairs"] != 0
            or physics["geometry_exclusions_are_physics_events"] is not False):
        raise TransportError("initial physics fields must remain planned, not native observations")

    return {
        "transport": {
            "schema": TRANSPORT_SCHEMA,
            "audit_contract_json_path": str(audit_path),
            "audit_contract_json_sha256": actual_json_sha256,
            "static_contract_status": "VALID",
            "native_mask_eligibility": eligibility,
            "qualification_status": "TRANSPORT_ONLY_NATIVE_UNVERIFIED",
            "initialization_authorized": False,
        },
        "source": {
            "source_cad_path": source_cad_path,
            "source_cad_sha256": source_cad_sha,
            "classifier_source_path": classifier_source_path,
            "classifier_source_sha256": classifier_source_sha,
            "source_hash_status": "MATCHED_ACTUAL_FILE_BYTES",
            "source_status": "READY",
            "source_accuracy_mm": source_accuracy,
            "source_accuracy_basis": source_accuracy_basis,
        },
        "coordinates": {
            "input_unit": "mm",
            "native_unit": "m",
            "m_per_input_unit": 0.001,
            "reference_frame": reference_frame,
            "native_frame": native_frame,
            "native_frame_from_reference_frame_m": transform,
        },
        "points": {
            "payload_path": point_path,
            "payload_sha256": point_payload_sha,
            "order": "exact_payload_row_order",
            "point_ids_order_sha256": point_order_sha,
            "count": len(validated_points),
            "rows": validated_points,
        },
        "candidate_generation": {
            "metric": "euclidean_reference_distance",
            "horizon_m": horizon_m,
            "candidate_rule": "exact_pinned_Peridynamics_find_bonds_convention",
            "scope": "one_body_one_uniform_horizon",
        },
        "raw_euclidean_candidate_pairs": validated_candidates,
        "counts": expected_counts,
        "hashes": {
            "pair_map_sha256": root["pair_map_sha256"],
            "mask_sha256": mask_sha,
        },
        "mask": {
            "included_pair_count": len(native_mask_rows["included_pairs"]),
            "included_pairs": native_mask_rows["included_pairs"],
            "geometry_excluded_pair_count": len(native_mask_rows["geometry_excluded_pairs"]),
            "geometry_excluded_pairs": native_mask_rows["geometry_excluded_pairs"],
            "unknown_pair_count": len(native_mask_rows["unknown_pairs"]),
            "unknown_pairs": native_mask_rows["unknown_pairs"],
        },
        "evidence": {
            "classifier_result_path": classifier_result_path,
            "classifier_result_sha256": classifier_result_sha,
            "raw_hit_interval_sidecar_path": sidecar_path,
            "raw_hit_interval_sidecar_sha256": sidecar_sha,
            "sidecar_test_nodeid": test_nodeid,
            "sidecar_call_index": call_index,
            "visibility_algorithm": algorithm,
        },
        "native_gate": {
            "status": "UNVERIFIED",
            "raw_candidate_set_reproduced": False,
            "native_graph_readback": "NOT_RUN",
            "initialization_authorized": False,
            "zero_degree_external_point_ids": zero_degree_ids,
        },
        "initial_physics_events": {
            "evidence_status": "PLANNED_NOT_OBSERVED_NATIVE_NOT_RUN",
            "planned_physical_broken_bond_pairs": 0,
            "geometry_exclusions_are_physics_events": False,
        },
    }


def _toml_string(value: str) -> str:
    value.encode("utf-8", errors="strict")
    result: list[str] = ['"']
    for char in value:
        code = ord(char)
        if char == '"':
            result.append('\\"')
        elif char == "\\":
            result.append("\\\\")
        elif char == "\b":
            result.append("\\b")
        elif char == "\t":
            result.append("\\t")
        elif char == "\n":
            result.append("\\n")
        elif char == "\f":
            result.append("\\f")
        elif char == "\r":
            result.append("\\r")
        else:
            if unicodedata.category(char) == "Cc":
                raise TransportError("TOML strings cannot contain unhandled control characters")
            result.append(char)
    result.append('"')
    return "".join(result)


def _toml_value(value: Any) -> str:
    if type(value) is bool:
        return "true" if value else "false"
    if type(value) is int:
        if value < -(1 << 63) or value > (1 << 63) - 1:
            raise TransportError("integer value is outside TOML's signed 64-bit range")
        return str(value)
    if type(value) is float:
        if not math.isfinite(value):
            raise TransportError("TOML transport cannot contain NaN or infinity")
        return repr(value)
    if isinstance(value, str):
        return _toml_string(value)
    if isinstance(value, list):
        if any(isinstance(item, (dict,)) for item in value):
            raise TransportError("nested TOML table arrays must be emitted as array-of-table headers")
        return "[" + ", ".join(_toml_value(item) for item in value) + "]"
    raise TransportError(f"unsupported TOML transport value: {type(value).__name__}")


def _toml_key(key: str) -> str:
    if _KEY_RE.fullmatch(key) is None:
        raise TransportError(f"unsupported TOML bare key {key!r}")
    return key


def _emit_table_body(table: Mapping[str, Any], path: tuple[str, ...], lines: list[str]) -> None:
    scalar_fields: list[tuple[str, Any]] = []
    sub_tables: list[tuple[str, Mapping[str, Any]]] = []
    table_arrays: list[tuple[str, Sequence[Mapping[str, Any]]]] = []
    for key, value in table.items():
        safe_key = _toml_key(key)
        if isinstance(value, dict):
            sub_tables.append((safe_key, value))
        elif isinstance(value, list) and value and all(isinstance(item, dict) for item in value):
            table_arrays.append((safe_key, value))
        else:
            scalar_fields.append((safe_key, value))
    for key, value in scalar_fields:
        lines.append(f"{key} = {_toml_value(value)}")
    for key, child in sub_tables:
        child_path = path + (key,)
        lines.append("")
        lines.append("[" + ".".join(_toml_key(part) for part in child_path) + "]")
        _emit_table_body(child, child_path, lines)
    for key, rows in table_arrays:
        array_path = path + (key,)
        for row in rows:
            lines.append("")
            lines.append("[[" + ".".join(_toml_key(part) for part in array_path) + "]]")
            _emit_table_body(row, array_path, lines)


def _serialize_toml(model: Mapping[str, Any]) -> bytes:
    lines: list[str] = []
    _emit_table_body(model, (), lines)
    encoded = ("\n".join(lines).rstrip() + "\n").encode("utf-8", errors="strict")
    if len(encoded) > MAX_NATIVE_TOML_BYTES:
        raise TransportError(f"native TOML exceeds {MAX_NATIVE_TOML_BYTES}-byte output limit")
    return encoded


def _typed_equal(actual: Any, expected: Any, path: str = "root") -> None:
    if type(actual) is not type(expected):
        raise TransportError(f"TOML readback type mismatch at {path}")
    if isinstance(expected, dict):
        if set(actual) != set(expected):
            raise TransportError(f"TOML readback keys mismatch at {path}")
        for key in expected:
            _typed_equal(actual[key], expected[key], f"{path}.{key}")
        return
    if isinstance(expected, list):
        if len(actual) != len(expected):
            raise TransportError(f"TOML readback array length mismatch at {path}")
        for index, (actual_item, expected_item) in enumerate(zip(actual, expected)):
            _typed_equal(actual_item, expected_item, f"{path}[{index}]")
        return
    if isinstance(expected, float):
        if (not math.isfinite(actual) or actual != expected
                or math.copysign(1.0, actual) != math.copysign(1.0, expected)):
            raise TransportError(f"TOML readback numeric mismatch at {path}")
    elif actual != expected:
        raise TransportError(f"TOML readback value mismatch at {path}")


def _parse_and_compare_toml(data: bytes, expected: Mapping[str, Any], label: str) -> None:
    if len(data) > MAX_NATIVE_TOML_BYTES:
        raise TransportError(f"{label} exceeds the {MAX_NATIVE_TOML_BYTES}-byte limit")
    try:
        parsed = tomllib.loads(data.decode("utf-8", errors="strict"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        raise TransportError(f"{label} is not valid UTF-8 TOML") from exc
    _typed_equal(parsed, expected)


def _atomic_write_and_readback(target: Path, payload: bytes, expected: Mapping[str, Any]) -> bytes:
    if target.exists() or target.is_symlink():
        raise TransportError(f"refusing to overwrite existing TOML output: {target}")
    if not target.parent.is_dir():
        raise TransportError("output directory must already exist")
    temp_path: Path | None = None
    try:
        descriptor, name = tempfile.mkstemp(prefix=f".{target.name}.", suffix=".tmp", dir=target.parent)
        temp_path = Path(name)
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        temp_bytes = _read_stable_bounded(temp_path, "temporary TOML output", MAX_NATIVE_TOML_BYTES)
        if temp_bytes != payload:
            raise TransportError("temporary TOML bytes differ from generated bytes")
        _parse_and_compare_toml(temp_bytes, expected, "temporary TOML output")
        if target.exists() or target.is_symlink():
            raise TransportError(f"refusing to overwrite existing TOML output: {target}")
        os.replace(temp_path, target)
        temp_path = None
        final_bytes = _read_stable_bounded(target, "native TOML output", MAX_NATIVE_TOML_BYTES)
        if final_bytes != payload:
            raise TransportError("atomically written TOML bytes differ from generated bytes")
        _parse_and_compare_toml(final_bytes, expected, "native TOML output")
        return final_bytes
    except OSError as exc:
        raise TransportError(f"atomic TOML output failed: {target}") from exc
    finally:
        if temp_path is not None:
            try:
                temp_path.unlink()
            except FileNotFoundError:
                pass
            except OSError:
                pass


def compile_native_mask_transport(audit_contract_path: Path, output_path: Path) -> dict[str, Any]:
    """Compile one actual JSON audit file to TOML and atomically verify the written bytes."""
    try:
        audit_path = Path(audit_contract_path).resolve(strict=True)
        requested_target = Path(output_path).expanduser()
        target = requested_target.resolve(strict=False)
    except OSError as exc:
        raise TransportError("cannot resolve audit or output path") from exc
    if audit_path == target:
        raise TransportError("TOML output path must not overwrite the JSON audit contract")
    if target.suffix.lower() != ".toml":
        raise TransportError("native transport output path must end in .toml")
    if (requested_target.exists() or requested_target.is_symlink()
            or target.exists() or target.is_symlink()):
        raise TransportError(f"refusing to overwrite existing TOML output: {target}")
    raw = _read_stable_bounded(audit_path, "JSON audit contract", MAX_AUDIT_JSON_BYTES)
    actual_json_sha = sha256_bytes(raw)
    contract = _load_strict_json(raw)
    model = _validate_contract(contract, actual_json_sha, audit_path)
    toml_bytes = _serialize_toml(model)
    _parse_and_compare_toml(toml_bytes, model, "generated native TOML")
    final_bytes = _atomic_write_and_readback(target, toml_bytes, model)
    return {
        "status": "TRANSPORT_ONLY_NATIVE_UNVERIFIED",
        "audit_contract_json_path": str(audit_path),
        "audit_contract_json_sha256": actual_json_sha,
        "native_toml_path": str(target),
        "native_toml_sha256": sha256_bytes(final_bytes),
        "native_gate_status": "UNVERIFIED",
        "initialization_authorized": False,
        "point_count": len(model["points"]["rows"]),
        "raw_candidate_pair_count": model["counts"]["raw_euclidean_candidate_pairs"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("audit_contract_json", type=Path)
    parser.add_argument("native_transport_toml", type=Path)
    args = parser.parse_args()
    try:
        receipt = compile_native_mask_transport(args.audit_contract_json, args.native_transport_toml)
    except TransportError as exc:
        parser.error(str(exc))
    print(json.dumps(receipt, sort_keys=True, separators=(",", ":"), allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
