"""Source-only contract packager candidate for a future native bond-graph gate.

This module deliberately stops before Peridynamics.jl. It binds files and validates
the static classifier ledger; it is not wired into the production runner and has not
been executed as part of the native-bridge design task.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import os
import re
from pathlib import Path
from typing import Any, Mapping, Sequence


CONTRACT_SCHEMA = "pgfrag.material-visibility-native-mask.v1"
VISIBILITY_SCHEMA = "pgfrag.bond-visibility.v1"
STATES = frozenset(("IN_MATERIAL", "OUT_OR_VOID", "UNKNOWN"))
_POINT_IDS_PREFIX = b"pgfrag.point-ids.v1\n"
_PAIR_MAP_PREFIX = b"pgfrag.native-pair-map.v1\n"
_INCLUDED_MASK_PREFIX = b"pgfrag.native-included-mask.v1\n"
_POINT_TABLE_HEADER = ("point_id", "x_mm", "y_mm", "z_mm")
MAX_TOTAL_INPUT_BYTES = 64 * 1024 * 1024
MAX_CLASSIFIER_RESULT_BYTES = 8 * 1024 * 1024
MAX_POINT_PAYLOAD_BYTES = 2 * 1024 * 1024
MAX_SOURCE_CAD_BYTES = 16 * 1024 * 1024
MAX_INTERVAL_SIDECAR_BYTES = 32 * 1024 * 1024
MAX_CLASSIFIER_SOURCE_BYTES = 2 * 1024 * 1024
MAX_POINTS = 5000
MAX_PAIRS = 250_000
MAX_SIDECAR_RESULTS = 16


class ContractError(ValueError):
    """The static visibility inputs cannot safely gate native initialization."""


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _finite_float(value: Any) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    try:
        converted = float(value)
    except (OverflowError, ValueError):
        return None
    return converted if math.isfinite(converted) else None


def _file_metadata(stat_result: Any) -> tuple[int, int, int, int]:
    return (
        int(stat_result.st_size),
        int(stat_result.st_mtime_ns),
        int(getattr(stat_result, "st_dev", 0)),
        int(getattr(stat_result, "st_ino", 0)),
    )


def _resolve_input_path(raw_path: Path, label: str) -> Path:
    try:
        return Path(raw_path).resolve(strict=True)
    except OSError as exc:
        raise ContractError(f"cannot resolve {label} path: {raw_path}") from exc


def _read_once_bounded(path: Path, label: str, expected: tuple[int, int, int, int], limit: int) -> bytes:
    try:
        with path.open("rb") as stream:
            before_handle = _file_metadata(os.fstat(stream.fileno()))
            before_path = _file_metadata(path.stat())
            if before_handle != expected or before_path != expected or expected[0] > limit:
                raise ContractError(f"{label} changed after size preflight")
            data = stream.read(expected[0] + 1)
            if len(data) != expected[0]:
                raise ContractError(f"{label} size changed during bounded read")
            after_handle = _file_metadata(os.fstat(stream.fileno()))
        after_path = _file_metadata(path.stat())
    except OSError as exc:
        raise ContractError(f"cannot read {label}: {path}") from exc
    if after_handle != expected or after_path != expected:
        raise ContractError(f"{label} changed during bounded read")
    return data


def _read_stable_bundle(specs: Sequence[tuple[str, Path, int]]) -> dict[str, bytes]:
    """Stat every input first, enforce per-file/aggregate limits, then read twice bounded."""
    prepared: dict[str, tuple[Path, tuple[int, int, int, int], int]] = {}
    seen_paths: set[Path] = set()
    total_bytes = 0
    if not specs:
        raise ContractError("contract input bundle cannot be empty")
    for label, raw_path, limit in specs:
        path = Path(raw_path)
        try:
            resolved = _resolve_input_path(path, label)
            metadata = _file_metadata(resolved.stat())
        except OSError as exc:
            raise ContractError(f"cannot stat {label}: {path}") from exc
        if resolved in seen_paths:
            raise ContractError("each contract input role must name a distinct file")
        seen_paths.add(resolved)
        if metadata[0] < 0 or metadata[0] > limit:
            raise ContractError(f"{label} exceeds its {limit}-byte input limit")
        total_bytes += metadata[0]
        if total_bytes > MAX_TOTAL_INPUT_BYTES:
            raise ContractError(f"contract inputs exceed the {MAX_TOTAL_INPUT_BYTES}-byte aggregate limit")
        prepared[label] = (resolved, metadata, limit)

    captured: dict[str, bytes] = {}
    for label, (path, metadata, limit) in prepared.items():
        first = _read_once_bounded(path, label, metadata, limit)
        second = _read_once_bounded(path, label, metadata, limit)
        if first != second:
            raise ContractError(f"{label} changed between bounded reads")
        captured[label] = first
    return captured


def _load_strict_json(data: bytes, label: str) -> Any:
    def reject_constant(value: str) -> None:
        raise ContractError(f"{label} contains non-standard JSON number {value}")

    def reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise ContractError(f"{label} contains duplicate object key {key!r}")
            result[key] = value
        return result

    try:
        return json.loads(
            data.decode("utf-8"),
            parse_constant=reject_constant,
            object_pairs_hook=reject_duplicate_keys,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ContractError(f"{label} is not strict UTF-8 JSON") from exc


def _point_rows(payload: bytes) -> list[dict[str, Any]]:
    try:
        decoded = payload.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ContractError("point payload must be UTF-8 CSV") from exc
    try:
        reader = csv.reader(io.StringIO(decoded, newline=""), strict=True)
        header = next(reader, None)
    except csv.Error as exc:
        raise ContractError("point payload is malformed CSV") from exc
    if header is None or tuple(header) != _POINT_TABLE_HEADER:
        raise ContractError("point payload header must be point_id,x_mm,y_mm,z_mm")
    parsed: list[dict[str, Any]] = []
    seen: set[int] = set()
    try:
        rows = enumerate(reader, start=2)
        for line_number, row in rows:
            if line_number - 1 > MAX_POINTS:
                raise ContractError(f"point payload exceeds the {MAX_POINTS}-point LocalTiny limit")
            if len(row) != 4:
                raise ContractError(f"point payload row {line_number} must have four columns")
            if re.fullmatch(r"(?:0|[1-9][0-9]*)", row[0]) is None:
                raise ContractError(f"point payload row {line_number} ID is not canonical non-negative decimal")
            try:
                point_id = int(row[0], 10)
                xyz = [float(value) for value in row[1:]]
            except ValueError as exc:
                raise ContractError(f"point payload row {line_number} has an invalid value") from exc
            if point_id in seen:
                raise ContractError(f"point payload row {line_number} has a duplicate ID")
            if any(not math.isfinite(value) for value in xyz):
                raise ContractError(f"point payload row {line_number} has a non-finite coordinate")
            seen.add(point_id)
            parsed.append({"point_id": point_id, "native_index": len(parsed) + 1, "position_mm": xyz})
    except csv.Error as exc:
        raise ContractError("point payload is malformed CSV") from exc
    if not parsed:
        raise ContractError("point payload must contain at least one point")
    return parsed


def _point_id_order_sha256(rows: Sequence[Mapping[str, Any]]) -> str:
    payload = bytearray(_POINT_IDS_PREFIX)
    for row in rows:
        payload.extend(f"{row['point_id']}\n".encode("ascii"))
    return sha256_bytes(bytes(payload))


def _canonical_pair(value: Any, label: str) -> tuple[int, int]:
    if not isinstance(value, list) or len(value) != 2:
        raise ContractError(f"{label} must be a two-ID array")
    left, right = value
    if (isinstance(left, bool) or not isinstance(left, int) or left < 0
            or isinstance(right, bool) or not isinstance(right, int) or right < 0
            or left >= right):
        raise ContractError(f"{label} must be canonical non-negative IDs with i < j")
    return left, right


def _pair_map_sha256(decisions: Sequence[Mapping[str, Any]]) -> str:
    payload = bytearray(_PAIR_MAP_PREFIX)
    for row in decisions:
        i, j = row["pair"]
        payload.extend(f"{i},{j},{row['classification']}\n".encode("ascii"))
    return sha256_bytes(bytes(payload))


def _included_mask_sha256(included_pairs: Sequence[Sequence[int]]) -> str:
    payload = bytearray(_INCLUDED_MASK_PREFIX)
    for i, j in included_pairs:
        payload.extend(f"{i},{j}\n".encode("ascii"))
    return sha256_bytes(bytes(payload))


def _classifier_pairs(
    result: Mapping[str, Any], point_rows: Sequence[Mapping[str, Any]]
) -> list[dict[str, Any]]:
    raw_pairs = result.get("pairs")
    if not isinstance(raw_pairs, list):
        raise ContractError("classifier result must contain a pairs array")
    if len(raw_pairs) > MAX_PAIRS:
        raise ContractError(f"classifier result exceeds the {MAX_PAIRS}-pair LocalTiny limit")
    positions = {row["point_id"]: row["position_mm"] for row in point_rows}
    decisions: list[dict[str, Any]] = []
    seen: set[tuple[int, int]] = set()
    for index, record in enumerate(raw_pairs):
        if not isinstance(record, Mapping):
            raise ContractError(f"classifier pair record {index} is not an object")
        pair = _canonical_pair(record.get("pair"), f"classifier pair record {index}.pair")
        if pair in seen:
            raise ContractError(f"duplicate classifier pair {pair}")
        seen.add(pair)
        state = record.get("classification")
        if not isinstance(state, str) or state not in STATES:
            raise ContractError(f"classifier pair {pair} has invalid classification")
        if record.get("admissible_static_candidate") is not (state == "IN_MATERIAL"):
            raise ContractError(f"classifier pair {pair} has inconsistent admissibility flag")
        if not isinstance(record.get("raw_hits"), list) or not isinstance(record.get("open_intervals"), list):
            raise ContractError(f"classifier pair {pair} is missing raw hit or interval evidence")
        reciprocal = record.get("reciprocal")
        if (not isinstance(reciprocal, Mapping)
                or reciprocal.get("same_classification") is not True
                or not isinstance(reciprocal.get("directed_records"), list)
                or len(reciprocal["directed_records"]) != 2):
            raise ContractError(f"classifier pair {pair} lacks reciprocal direction evidence")
        directions = reciprocal["directed_records"]
        expected_directions = {(pair[0], pair[1]), (pair[1], pair[0])}
        observed_directions: set[tuple[int, int]] = set()
        for direction in directions:
            if not isinstance(direction, Mapping):
                raise ContractError(f"classifier pair {pair} has malformed directed evidence")
            source_id, target_id, directed_state = (
                direction.get("source_point_id"),
                direction.get("target_point_id"),
                direction.get("state"),
            )
            if (type(source_id) is not int or type(target_id) is not int
                    or directed_state != state):
                raise ContractError(f"classifier pair {pair} has inconsistent directed evidence")
            observed_directions.add((source_id, target_id))
        if observed_directions != expected_directions:
            raise ContractError(f"classifier pair {pair} is missing a reciprocal direction")
        if pair[0] not in positions or pair[1] not in positions:
            raise ContractError(f"classifier pair {pair} references an unknown point ID")
        endpoints = record.get("endpoints")
        if not isinstance(endpoints, list) or len(endpoints) != 2:
            raise ContractError(f"classifier pair {pair} has no complete endpoint rows")
        for endpoint, expected_id in zip(endpoints, pair):
            endpoint_id = endpoint.get("point_id") if isinstance(endpoint, Mapping) else None
            if type(endpoint_id) is not int or endpoint_id != expected_id:
                raise ContractError(f"classifier pair {pair} endpoint ID/order mismatch")
            xyz = endpoint.get("position_mm")
            if not isinstance(xyz, list) or len(xyz) != 3:
                raise ContractError(f"classifier pair {pair} endpoint coordinates are missing")
            endpoint_xyz = [_finite_float(v) for v in xyz]
            if any(value is None for value in endpoint_xyz):
                raise ContractError(f"classifier pair {pair} has a non-finite endpoint")
            if endpoint_xyz != positions[expected_id]:
                raise ContractError(f"classifier pair {pair} endpoint does not match point payload")
        reason = record.get("reason")
        if not isinstance(reason, str) or not reason:
            raise ContractError(f"classifier pair {pair} has no reason")
        decisions.append({
            "pair": [pair[0], pair[1]],
            "classification": state,
            "reason": reason,
        })
    decisions.sort(key=lambda row: tuple(row["pair"]))
    return decisions


def _canonical_json_bytes(value: Any, label: str) -> bytes:
    try:
        return json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
            allow_nan=False,
        ).encode("utf-8")
    except (TypeError, ValueError, UnicodeEncodeError) as exc:
        raise ContractError(f"{label} is not strict canonicalizable JSON data") from exc


def _validate_geometry_sidecar(
    *,
    sidecar_bytes: bytes,
    classifier_result: Mapping[str, Any],
    classifier_source_path: Path,
    selection: Mapping[str, Any],
) -> dict[str, Any]:
    sidecar = _load_strict_json(sidecar_bytes, "raw-hit interval sidecar")
    if (not isinstance(sidecar, Mapping)
            or sidecar.get("schema") != "pgfrag.bond-visibility-geometry-sidecar.v1"):
        raise ContractError("unsupported raw-hit interval sidecar schema")
    if type(sidecar.get("pytest_exit_code")) is not int or sidecar["pytest_exit_code"] != 0:
        raise ContractError("raw-hit interval sidecar does not record a successful captured run")
    records = sidecar.get("results")
    if not isinstance(records, list) or not records or len(records) > MAX_SIDECAR_RESULTS:
        raise ContractError(f"raw-hit interval sidecar must contain 1..{MAX_SIDECAR_RESULTS} captured results")
    if type(sidecar.get("captured_call_count")) is not int or sidecar["captured_call_count"] != len(records):
        raise ContractError("sidecar captured_call_count does not match its results array")
    try:
        module_path = Path(sidecar.get("module_file", "")).resolve(strict=True)
        expected_module_path = Path(classifier_source_path).resolve(strict=True)
    except (OSError, TypeError) as exc:
        raise ContractError("sidecar module_file cannot be resolved") from exc
    if module_path != expected_module_path:
        raise ContractError("sidecar module_file does not name the hashed classifier source")
    algorithm = classifier_result.get("algorithm")
    if not isinstance(algorithm, str) or sidecar.get("algorithm") != algorithm:
        raise ContractError("sidecar algorithm does not match classifier result")
    if not isinstance(selection, Mapping):
        raise ContractError("sidecar selection must identify one test_nodeid and call_index")
    test_nodeid = selection.get("test_nodeid")
    call_index = selection.get("call_index")
    if not isinstance(test_nodeid, str) or not test_nodeid or type(call_index) is not int or call_index < 0:
        raise ContractError("sidecar selection must identify one test_nodeid and non-negative call_index")
    seen_selectors: set[tuple[str, int]] = set()
    matches: list[Mapping[str, Any]] = []
    for record in records:
        if not isinstance(record, Mapping):
            raise ContractError("sidecar result entry is not an object")
        nodeid = record.get("test_nodeid")
        index = record.get("call_index")
        if not isinstance(nodeid, str) or type(index) is not int or index < 0:
            raise ContractError("sidecar result entry has invalid test_nodeid/call_index")
        selector = (nodeid, index)
        if selector in seen_selectors:
            raise ContractError("sidecar contains duplicate test_nodeid/call_index entries")
        seen_selectors.add(selector)
        if selector == (test_nodeid, call_index):
            captured = record.get("result")
            if not isinstance(captured, Mapping):
                raise ContractError("selected sidecar entry has no classifier result")
            matches.append(captured)
    if len(matches) != 1:
        raise ContractError("sidecar selection did not resolve to exactly one result")
    if _canonical_json_bytes(matches[0], "selected sidecar result") != _canonical_json_bytes(
            classifier_result, "classifier result"):
        raise ContractError("selected sidecar result differs from classifier result, including pair hit/interval evidence")
    return {"test_nodeid": test_nodeid, "call_index": call_index}


def build_contract(
    *,
    classifier_result_path: Path,
    point_payload_path: Path,
    source_cad_path: Path,
    raw_hit_interval_sidecar_path: Path,
    classifier_source_path: Path,
    sidecar_selection: Mapping[str, Any],
    horizon_m: float,
    reference_frame: str,
    native_frame: str,
    native_frame_from_reference_frame_m: Sequence[Sequence[float]],
) -> dict[str, Any]:
    """Build a pinned mask contract; UNKNOWN remains present and blocks native eligibility.

    The caller supplies the exact inputs associated with the completed static
    classifier result. This function only packages and validates those inputs; a
    future native constructor must still independently reproduce the horizon
    candidate set and read the contract back before initialization.
    """
    resolved_paths = {
        "classifier result": _resolve_input_path(classifier_result_path, "classifier result"),
        "ordered point payload": _resolve_input_path(point_payload_path, "ordered point payload"),
        "source CAD": _resolve_input_path(source_cad_path, "source CAD"),
        "raw-hit interval sidecar": _resolve_input_path(raw_hit_interval_sidecar_path, "raw-hit interval sidecar"),
        "classifier source": _resolve_input_path(classifier_source_path, "classifier source"),
    }
    resolved_source_path = resolved_paths["classifier source"]
    captured = _read_stable_bundle((
        ("classifier result", resolved_paths["classifier result"], MAX_CLASSIFIER_RESULT_BYTES),
        ("ordered point payload", resolved_paths["ordered point payload"], MAX_POINT_PAYLOAD_BYTES),
        ("source CAD", resolved_paths["source CAD"], MAX_SOURCE_CAD_BYTES),
        ("raw-hit interval sidecar", resolved_paths["raw-hit interval sidecar"], MAX_INTERVAL_SIDECAR_BYTES),
        ("classifier source", resolved_source_path, MAX_CLASSIFIER_SOURCE_BYTES),
    ))
    classifier_bytes = captured["classifier result"]
    point_bytes = captured["ordered point payload"]
    cad_bytes = captured["source CAD"]
    interval_bytes = captured["raw-hit interval sidecar"]
    classifier_source_bytes = captured["classifier source"]
    result = _load_strict_json(classifier_bytes, "classifier result")
    if (not isinstance(result, Mapping)
            or result.get("schema") != VISIBILITY_SCHEMA
            or result.get("status") != "STATIC_CLASSIFICATION_COMPLETE"):
        raise ContractError("unexpected or incomplete classifier result")
    source = result.get("source")
    reference = result.get("reference_coordinates")
    topology = result.get("initial_topology_ledger")
    if not isinstance(source, Mapping) or not isinstance(reference, Mapping) or not isinstance(topology, Mapping):
        raise ContractError("classifier result is missing source, reference, or topology metadata")

    cad_sha = sha256_bytes(cad_bytes)
    if (source.get("source_origin") != "CANONICAL_STEP_FILE"
            or source.get("status") != "READY"
            or source.get("sha256_status") != "MATCHED_ACTUAL_FILE_BYTES"
            or source.get("source_geometry_sha256") != cad_sha
            or source.get("expected_sha256") != cad_sha):
        raise ContractError("source CAD must be READY and SHA-bound to the exact input file")
    if source.get("source_length_unit") != "mm" or reference.get("coordinate_unit") != "mm":
        raise ContractError("classifier and native bridge input coordinates must both be mm")
    if (not isinstance(reference_frame, str) or not reference_frame
            or not isinstance(native_frame, str) or not native_frame
            or reference.get("frame") != reference_frame or source.get("frame") != reference_frame):
        raise ContractError("reference/native frame identifiers are missing or reference differs from classifier")
    accuracy = _finite_float(source.get("source_accuracy_mm"))
    accuracy_basis = source.get("source_accuracy_basis")
    if (accuracy is None or accuracy < 0
            or not isinstance(accuracy_basis, str) or not accuracy_basis.strip()):
        raise ContractError("source accuracy and its declared basis must be finite and explicit")
    horizon_value = _finite_float(horizon_m)
    if horizon_value is None or horizon_value <= 0:
        raise ContractError("horizon_m must be finite and positive")
    if (not isinstance(native_frame_from_reference_frame_m, (list, tuple))
            or any(not isinstance(row, (list, tuple)) for row in native_frame_from_reference_frame_m)):
        raise ContractError("native frame transform must be a finite 4x4 matrix")
    transform = [list(row) for row in native_frame_from_reference_frame_m]
    finite_transform = [[_finite_float(value) for value in row] for row in transform]
    if (len(transform) != 4 or any(len(row) != 4 for row in transform)
            or any(value is None for row in finite_transform for value in row)):
        raise ContractError("native frame transform must be a finite 4x4 matrix")
    matrix = finite_transform
    if matrix[3] != [0.0, 0.0, 0.0, 1.0]:
        raise ContractError("native frame transform must be homogeneous affine form")
    rotation = [row[:3] for row in matrix[:3]]
    for i in range(3):
        for j in range(3):
            dot = sum(rotation[k][i] * rotation[k][j] for k in range(3))
            expected = 1.0 if i == j else 0.0
            if not math.isclose(dot, expected, rel_tol=0.0, abs_tol=1e-12):
                raise ContractError("native frame transform must have an orthonormal rotation")
    determinant = (
        rotation[0][0] * (rotation[1][1] * rotation[2][2] - rotation[1][2] * rotation[2][1])
        - rotation[0][1] * (rotation[1][0] * rotation[2][2] - rotation[1][2] * rotation[2][0])
        + rotation[0][2] * (rotation[1][0] * rotation[2][1] - rotation[1][1] * rotation[2][0])
    )
    if not math.isclose(determinant, 1.0, rel_tol=0.0, abs_tol=1e-12):
        raise ContractError("native frame transform must preserve orientation")

    point_rows = _point_rows(point_bytes)
    point_ids = [row["point_id"] for row in point_rows]
    if len(point_ids) != len(set(point_ids)) or any(type(value) is not int or value < 0 for value in point_ids):
        raise ContractError("point IDs must be unique non-negative integers in exact Body row order")
    point_count = reference.get("point_count")
    if type(point_count) is not int or point_count != len(point_rows):
        raise ContractError("point payload count differs from classifier reference count")
    decisions = _classifier_pairs(result, point_rows)
    counts = {state: sum(row["classification"] == state for row in decisions) for state in STATES}
    reported = topology.get("classification_counts")
    if (not isinstance(reported, Mapping)
            or any(type(reported.get(state)) is not int or reported.get(state) != counts[state]
                   for state in STATES)):
        raise ContractError("classifier topology counts disagree with pair decisions")
    if (type(topology.get("candidate_undirected_pair_count")) is not int
            or topology.get("candidate_undirected_pair_count") != len(decisions)):
        raise ContractError("classifier candidate count disagrees with pair decisions")
    if any(type(topology.get(key)) is not int for key in (
            "initial_material_bond_candidate_count", "geometry_excluded_pair_count", "unknown_pair_count")):
        raise ContractError("classifier topology summary counts must be integers")
    if (topology["initial_material_bond_candidate_count"] != counts["IN_MATERIAL"]
            or topology["geometry_excluded_pair_count"] != counts["OUT_OR_VOID"]
            or topology["unknown_pair_count"] != counts["UNKNOWN"]):
        raise ContractError("classifier topology summary counts disagree with pair decisions")

    algorithm = result.get("algorithm")
    if not isinstance(algorithm, str) or not algorithm:
        raise ContractError("classifier algorithm identifier is missing")
    sidecar_selector = _validate_geometry_sidecar(
        sidecar_bytes=interval_bytes,
        classifier_result=result,
        classifier_source_path=resolved_source_path,
        selection=sidecar_selection,
    )
    classifier_source_sha256 = sha256_bytes(classifier_source_bytes)
    evidence = {
        "classifier_result_sha256": sha256_bytes(classifier_bytes),
        "classifier_result_path": str(resolved_paths["classifier result"]),
        "classifier_source_path": str(resolved_source_path),
        "classifier_source_sha256": classifier_source_sha256,
        "raw_hit_interval_sidecar_path": str(resolved_paths["raw-hit interval sidecar"]),
        "raw_hit_interval_sidecar_sha256": sha256_bytes(interval_bytes),
        "sidecar_selection": sidecar_selector,
        "visibility_algorithm": algorithm,
    }
    included_pairs = [row["pair"] for row in decisions if row["classification"] == "IN_MATERIAL"]
    excluded_pairs = [row["pair"] for row in decisions if row["classification"] == "OUT_OR_VOID"]
    unknown_pairs = [row["pair"] for row in decisions if row["classification"] == "UNKNOWN"]
    included_degrees = {point_id: 0 for point_id in point_ids}
    for left, right in included_pairs:
        included_degrees[left] += 1
        included_degrees[right] += 1
    zero_degree_point_ids = [point_id for point_id in point_ids if included_degrees[point_id] == 0]
    if counts["UNKNOWN"]:
        native_mask_eligibility = "BLOCKED_UNKNOWN"
    elif not decisions:
        native_mask_eligibility = "BLOCKED_EMPTY_CANDIDATE_SET"
    elif not included_pairs:
        native_mask_eligibility = "BLOCKED_EMPTY_MASK"
    elif zero_degree_point_ids:
        native_mask_eligibility = "BLOCKED_ZERO_DEGREE"
    else:
        native_mask_eligibility = "CANDIDATE_MASK_POSITIVE_DEGREE_NATIVE_UNVERIFIED"
    return {
        "schema": CONTRACT_SCHEMA,
        "static_contract_status": "VALID",
        "native_mask_eligibility": native_mask_eligibility,
        "native_gate": {
            "status": "UNVERIFIED",
            "raw_candidate_set_reproduced": False,
            "native_graph_readback": "NOT_RUN",
            "initialization_authorized": False,
            "candidate_pair_count_nonzero": bool(decisions),
            "included_pair_count_nonzero": bool(included_pairs),
            "zero_degree_point_ids": zero_degree_point_ids,
        },
        "source": {
            "source_cad_path": str(resolved_paths["source CAD"]),
            "source_cad_sha256": cad_sha,
            "classifier_source_sha256": classifier_source_sha256,
            "source_hash_status": source["sha256_status"],
            "source_status": source["status"],
            "source_accuracy_mm": accuracy,
            "source_accuracy_basis": accuracy_basis,
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
            "payload_path": str(resolved_paths["ordered point payload"]),
            "payload_sha256": sha256_bytes(point_bytes),
            "order": "exact_payload_row_order",
            "point_ids_order_sha256": _point_id_order_sha256(point_rows),
            "rows": point_rows,
        },
        "candidate_generation": {
            "metric": "euclidean_reference_distance",
            "horizon_m": horizon_value,
            "candidate_rule": "exact_pinned_Peridynamics_find_bonds_convention",
            "scope": "one_body_one_uniform_horizon",
        },
        "raw_euclidean_candidate_pairs": decisions,
        "counts": {
            "raw_euclidean_candidate_pairs": len(decisions),
            "initial_material_bond_pairs": counts["IN_MATERIAL"],
            "geometry_excluded_pairs": counts["OUT_OR_VOID"],
            "unknown_pairs": counts["UNKNOWN"],
        },
        "pair_map_sha256": _pair_map_sha256(decisions),
        "mask": {
            "included_pairs": included_pairs,
            "geometry_excluded_pairs": excluded_pairs,
            "unknown_pairs": unknown_pairs,
            "mask_sha256": _included_mask_sha256(included_pairs),
        },
        "evidence": evidence,
        "initial_physics_events": {
            "evidence_status": "PLANNED_NOT_OBSERVED_NATIVE_NOT_RUN",
            "planned_physical_broken_bond_pairs": 0,
            "geometry_exclusions_are_physics_events": False,
        },
    }


def strict_contract_json_bytes(contract: Mapping[str, Any]) -> bytes:
    """Serialize deterministically; reject NaN/Infinity before a contract is pinned."""
    try:
        value = json.dumps(
            contract,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
            allow_nan=False,
        )
    except (TypeError, ValueError) as exc:
        raise ContractError("contract cannot be serialized as strict JSON") from exc
    return (value + "\n").encode("utf-8")
