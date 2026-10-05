"""Bounded, fixture-only binding of the saved planar split meshes and pair ledger.

This deliberately validates one exported two-piece OBJ fixture.  It does not
establish a general partition proof, physical qualification, or solver coupling.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import sys
from pathlib import Path
from typing import Any


MAX_FILE_BYTES = 256 * 1024
MAX_FRAGMENTS = 2
MAX_VERTICES_TOTAL = 64
MAX_TRIANGLES_TOTAL = 64
COORD_TOL = 1.0e-12
METRIC_REL_TOL = 1.0e-10
METRIC_ABS_TOL = 1.0e-12
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


class BindingValidationError(ValueError):
    """A bounded input is incomplete, inconsistent, or fails geometry checks."""


def _fail(message: str) -> None:
    raise BindingValidationError(message)


def _read_bounded(path: Path, *, role: str, max_bytes: int = MAX_FILE_BYTES) -> bytes:
    """Preflight size, then read at most cap+1 bytes before parsing/copying."""
    cap = min(MAX_FILE_BYTES, max_bytes)
    try:
        stat = path.stat()
    except FileNotFoundError:
        _fail(f"missing {role}: {path.name}")
    if not path.is_file():
        _fail(f"{role} is not a regular file: {path.name}")
    if stat.st_size > cap:
        _fail(f"{role} exceeds {cap} byte limit: {path.name}")
    try:
        with path.open("rb") as stream:
            raw = stream.read(cap + 1)
    except OSError as exc:
        _fail(f"cannot read {role} {path.name}: {exc}")
    if len(raw) > cap:
        _fail(f"{role} grew beyond {cap} byte limit while reading: {path.name}")
    return raw


def _reject_constant(token: str) -> None:
    raise BindingValidationError(f"non-finite JSON token is forbidden: {token}")


def _unique_object_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            _fail(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def _load_json_bytes(raw: bytes, *, role: str) -> Any:
    try:
        decoded = raw.decode("utf-8", errors="strict")
        value = json.loads(
            decoded,
            object_pairs_hook=_unique_object_pairs,
            parse_constant=_reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        _fail(f"invalid UTF-8 JSON in {role}: {exc}")
    _assert_finite_tree(value, role=role)
    return value


def _assert_finite_tree(value: Any, *, role: str) -> None:
    if isinstance(value, float) and not math.isfinite(value):
        _fail(f"non-finite number in {role}")
    if isinstance(value, dict):
        for child in value.values():
            _assert_finite_tree(child, role=role)
    elif isinstance(value, list):
        for child in value:
            _assert_finite_tree(child, role=role)


def _mapping(value: Any, *, field: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        _fail(f"{field} must be an object")
    return value


def _list(value: Any, *, field: str) -> list[Any]:
    if not isinstance(value, list):
        _fail(f"{field} must be an array")
    return value


def _text(value: Any, *, field: str) -> str:
    if not isinstance(value, str) or not value:
        _fail(f"{field} must be a non-empty string")
    return value


def _finite_number(value: Any, *, field: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        _fail(f"{field} must be numeric")
    result = float(value)
    if not math.isfinite(result):
        _fail(f"{field} must be finite")
    return result


def _finite_obj_number(value: str, *, field: str) -> float:
    try:
        result = float(value)
    except ValueError:
        _fail(f"{field} must be a decimal number")
    if not math.isfinite(result):
        _fail(f"{field} must be finite")
    return result


def _vec3(value: Any, *, field: str) -> tuple[float, float, float]:
    items = _list(value, field=field)
    if len(items) != 3:
        _fail(f"{field} must contain exactly three coordinates")
    return tuple(_finite_number(item, field=field) for item in items)  # type: ignore[return-value]


def _sha256(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _check_hash(raw: bytes, expected: Any, *, role: str) -> str:
    if not isinstance(expected, str) or not _SHA256_RE.fullmatch(expected):
        _fail(f"{role} expected_sha256 must be 64 lowercase hexadecimal characters")
    actual = _sha256(raw)
    if actual != expected:
        _fail(f"{role} SHA-256 mismatch: expected {expected}, got {actual}")
    return actual


def _source_path(root: Path, name: Any, *, role: str) -> Path:
    filename = _text(name, field=f"{role}.file")
    if Path(filename).name != filename or filename in {".", ".."}:
        _fail(f"{role}.file must be a basename within the fixture directory")
    path = (root / filename).resolve()
    try:
        path.relative_to(root.resolve())
    except ValueError:
        _fail(f"{role}.file resolves outside the fixture directory")
    return path


def _near(left: float, right: float, *, rel: float = METRIC_REL_TOL, abs_: float = METRIC_ABS_TOL) -> bool:
    return math.isclose(left, right, rel_tol=rel, abs_tol=abs_)


def _same_vec(left: tuple[float, float, float], right: tuple[float, float, float]) -> bool:
    return all(a == b for a, b in zip(left, right))


def _add(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def _sub(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _cross(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _dot(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _norm(a: tuple[float, float, float]) -> float:
    return math.sqrt(_dot(a, a))


def _unit(a: tuple[float, float, float], *, field: str) -> tuple[float, float, float]:
    magnitude = _norm(a)
    if not math.isfinite(magnitude) or magnitude <= 0.0:
        _fail(f"{field} has zero or invalid magnitude")
    return (a[0] / magnitude, a[1] / magnitude, a[2] / magnitude)


def _cyclic_equal(actual: tuple[int, int, int], expected: tuple[int, int, int]) -> bool:
    return actual in (
        expected,
        (expected[1], expected[2], expected[0]),
        (expected[2], expected[0], expected[1]),
    )


def _cyclic_reverse_equal(actual: tuple[int, int, int], expected: tuple[int, int, int]) -> bool:
    reversed_expected = (expected[0], expected[2], expected[1])
    return _cyclic_equal(actual, reversed_expected)


def _parse_obj(
    raw: bytes,
    *,
    role: str,
    expected_object: str,
    max_vertices: int,
    max_triangles: int,
) -> dict[str, Any]:
    try:
        text = raw.decode("utf-8", errors="strict")
    except UnicodeDecodeError as exc:
        _fail(f"invalid UTF-8 OBJ for {role}: {exc}")
    vertices: list[tuple[float, float, float]] = []
    faces: list[tuple[int, int, int]] = []
    groups: list[str] = []
    active_group = ""
    object_name: str | None = None
    for line_no, line in enumerate(text.splitlines(), 1):
        parts = line.split()
        if not parts or parts[0].startswith("#"):
            continue
        kind = parts[0]
        if kind == "o":
            if len(parts) != 2 or object_name is not None:
                _fail(f"{role} OBJ line {line_no}: expected one simple object name")
            object_name = parts[1]
        elif kind == "g":
            if len(parts) != 2:
                _fail(f"{role} OBJ line {line_no}: group must have one name")
            active_group = parts[1]
        elif kind == "v":
            if len(parts) != 4:
                _fail(f"{role} OBJ line {line_no}: only three-coordinate vertices are allowed")
            if len(vertices) >= max_vertices:
                _fail(f"{role} OBJ exceeds the {max_vertices} vertex cap")
            vertices.append(tuple(_finite_obj_number(part, field=f"{role} OBJ line {line_no}") for part in parts[1:]))  # type: ignore[arg-type]
        elif kind == "f":
            if len(parts) != 4:
                _fail(f"{role} OBJ line {line_no}: only triangular faces are allowed")
            if len(faces) >= max_triangles:
                _fail(f"{role} OBJ exceeds the {max_triangles} triangle cap")
            indices: list[int] = []
            for token in parts[1:]:
                if not token.isdecimal() or token.startswith("0") or len(token) > 2:
                    _fail(f"{role} OBJ line {line_no}: face indices must be positive absolute integers")
                index = int(token) - 1
                if index < 0 or index >= len(vertices):
                    _fail(f"{role} OBJ line {line_no}: face index is out of range")
                indices.append(index)
            if len(set(indices)) != 3:
                _fail(f"{role} OBJ line {line_no}: degenerate repeated face index")
            faces.append(tuple(indices))  # type: ignore[arg-type]
            groups.append(active_group)
        else:
            _fail(f"{role} OBJ line {line_no}: unsupported statement {kind!r}")
    if object_name != expected_object:
        _fail(f"{role} OBJ object name mismatch: expected {expected_object!r}, got {object_name!r}")
    if not vertices or not faces:
        _fail(f"{role} OBJ must contain vertices and triangles")
    return {"vertices": vertices, "faces": faces, "groups": groups, "object": object_name}


def _mesh_topology(mesh: dict[str, Any], *, role: str) -> dict[str, Any]:
    faces: list[tuple[int, int, int]] = mesh["faces"]
    edge_directions: dict[tuple[int, int], list[tuple[int, int]]] = {}
    for face in faces:
        for start, end in ((face[0], face[1]), (face[1], face[2]), (face[2], face[0])):
            key = (min(start, end), max(start, end))
            edge_directions.setdefault(key, []).append((start, end))
    if not edge_directions:
        _fail(f"{role} has no edges")
    for edge, directions in edge_directions.items():
        if len(directions) != 2:
            _fail(f"{role} is not closed manifold at edge {edge}: incidence {len(directions)}")
        if directions[0] != (directions[1][1], directions[1][0]):
            _fail(f"{role} has inconsistent directed edge winding at edge {edge}")

    vertices: list[tuple[float, float, float]] = mesh["vertices"]
    signed_volume = 0.0
    triangle_areas: list[float] = []
    triangle_normals: list[tuple[float, float, float]] = []
    for face in faces:
        a, b, c = (vertices[index] for index in face)
        signed_volume += _dot(a, _cross(b, c)) / 6.0
        cross = _cross(_sub(b, a), _sub(c, a))
        twice_area = _norm(cross)
        if not math.isfinite(twice_area) or twice_area <= 0.0:
            _fail(f"{role} has zero-area or invalid triangle")
        triangle_areas.append(twice_area / 2.0)
        triangle_normals.append(_unit(cross, field=f"{role} triangle normal"))
    if not math.isfinite(signed_volume) or signed_volume <= 0.0:
        _fail(f"{role} is closed but inward-oriented or has non-positive signed volume ({signed_volume})")
    return {
        "edge_count": len(edge_directions),
        "closed_oriented": True,
        "signed_volume_m3": signed_volume,
        "triangle_areas_m2": triangle_areas,
        "triangle_normals": triangle_normals,
    }


def _load_source_json(
    root: Path, ref: Any, *, role: str, max_bytes: int = MAX_FILE_BYTES
) -> tuple[dict[str, Any], str]:
    ref_obj = _mapping(ref, field=role)
    path = _source_path(root, ref_obj.get("file"), role=role)
    raw = _read_bounded(path, role=role, max_bytes=max_bytes)
    actual_hash = _check_hash(raw, ref_obj.get("sha256"), role=role)
    value = _mapping(_load_json_bytes(raw, role=role), field=role)
    return value, actual_hash


def validate_binding(binding_manifest_path: str | os.PathLike[str]) -> dict[str, Any]:
    """Verify one bounded, hash-pinned, two-piece planar OBJ fixture export."""
    manifest_path = Path(binding_manifest_path).resolve()
    manifest_raw = _read_bounded(manifest_path, role="binding manifest")
    manifest = _mapping(_load_json_bytes(manifest_raw, role="binding manifest"), field="binding manifest")
    if manifest.get("schema_version") != "pgfrag.planar-mesh-binding.v1":
        _fail("unsupported binding manifest schema_version")
    if manifest.get("qualification") != "FIXTURE_GEOMETRY_ONLY":
        _fail("this validator accepts FIXTURE_GEOMETRY_ONLY inputs only")
    if manifest.get("mechanics_qualification") != "NOT_COUPLED":
        _fail("fixture binding must explicitly declare mechanics_qualification=NOT_COUPLED")
    if type(manifest.get("physical_pack_count")) is not int or manifest.get("physical_pack_count") != 0:
        _fail("fixture binding must declare physical_pack_count=0")
    if manifest.get("mesh_coordinate_unit") != "m":
        _fail("mesh_coordinate_unit must be exactly 'm' for this metre-valued fixture")

    limits = _mapping(manifest.get("limits"), field="limits")
    limit_values = {
        "max_fragments": MAX_FRAGMENTS,
        "max_vertices_total": MAX_VERTICES_TOTAL,
        "max_triangles_total": MAX_TRIANGLES_TOTAL,
        "max_bytes_each": MAX_FILE_BYTES,
    }
    for field, cap in limit_values.items():
        value = limits.get(field)
        if isinstance(value, bool) or not isinstance(value, int) or value <= 0 or value > cap:
            _fail(f"limits.{field} must be an integer in [1, {cap}]")

    fixture_dir_value = _text(manifest.get("source_fixture_dir"), field="source_fixture_dir")
    fixture_dir_path = Path(fixture_dir_value)
    if fixture_dir_path.is_absolute():
        _fail("source_fixture_dir must be relative to the binding manifest")
    fixture_root = (manifest_path.parent / fixture_dir_path).resolve()
    # The allowed input fixture lives under pgfrag; reject references elsewhere.
    package_root = manifest_path.parent.parent
    pgfrag_root = package_root.parent.parent.resolve()
    try:
        fixture_root.relative_to(pgfrag_root)
    except ValueError:
        _fail("source_fixture_dir must resolve within the pgfrag tree")
    if not fixture_root.is_dir():
        _fail("source_fixture_dir does not exist")

    output, output_hash = _load_source_json(
        fixture_root, manifest.get("output_file"), role="fixture output", max_bytes=limits["max_bytes_each"]
    )
    ledger, ledger_hash = _load_source_json(
        fixture_root, manifest.get("pair_ledger_file"), role="pair ledger", max_bytes=limits["max_bytes_each"]
    )

    if output.get("schema_version") != "pgfrag.explicit-fragment-fixture-output.v1":
        _fail("unsupported fixture output schema")
    if output.get("qualification") != "FIXTURE_GEOMETRY_ONLY":
        _fail("fixture output does not declare FIXTURE_GEOMETRY_ONLY")
    if output.get("mechanics_qualification") != "NOT_COUPLED":
        _fail("fixture output does not declare NOT_COUPLED")
    if output.get("physical_fracture_claim") != "NOT_CLAIMED":
        _fail("fixture output makes an unexpected physical fracture claim")
    if ledger.get("schema_version") != "pgfrag.explicit-fragment-interface-pair-ledger.v1":
        _fail("unsupported pair ledger schema")
    if ledger.get("qualification") != "FIXTURE_GEOMETRY_ONLY":
        _fail("pair ledger does not declare FIXTURE_GEOMETRY_ONLY")

    embedded_ledger = _mapping(output.get("interface_pair_ledger"), field="output.interface_pair_ledger")
    if embedded_ledger != ledger:
        _fail("embedded output pair ledger differs from the hash-pinned pair ledger")

    fragments_raw = _list(manifest.get("fragments"), field="fragments")
    if len(fragments_raw) != 2 or len(fragments_raw) > limits["max_fragments"]:
        _fail("binding must describe exactly two fragments and respect max_fragments")
    fragments_by_side: dict[str, dict[str, Any]] = {}
    for fragment_raw in fragments_raw:
        fragment = _mapping(fragment_raw, field="fragment")
        side = _text(fragment.get("side"), field="fragment.side")
        if side not in {"negative", "positive"} or side in fragments_by_side:
            _fail("fragments must have unique negative and positive sides")
        fragments_by_side[side] = fragment
    if set(fragments_by_side) != {"negative", "positive"}:
        _fail("fragments must include negative and positive sides")

    pieces = _list(output.get("pieces"), field="output.pieces")
    if len(pieces) != 2:
        _fail("fixture output must contain exactly two pieces")
    pieces_by_side: dict[str, dict[str, Any]] = {}
    meshes_by_side: dict[str, dict[str, Any]] = {}
    mesh_hashes: dict[str, str] = {}
    total_vertices = 0
    total_triangles = 0
    for piece_raw in pieces:
        piece = _mapping(piece_raw, field="output piece")
        side = _text(piece.get("side"), field="piece.side")
        if side in pieces_by_side or side not in fragments_by_side:
            _fail("fixture output has an unknown or duplicate piece side")
        fragment = fragments_by_side[side]
        piece_id = _text(piece.get("piece_id"), field=f"{side}.piece_id")
        if fragment.get("piece_id") != piece_id:
            _fail(f"{side} piece_id differs between binding and output")

        vertex_rows = _list(piece.get("vertices_m"), field=f"{side}.vertices_m")
        triangle_rows = _list(piece.get("triangles"), field=f"{side}.triangles")
        vertex_ids = _list(piece.get("vertex_ids"), field=f"{side}.vertex_ids")
        if not vertex_rows or len(vertex_rows) != len(vertex_ids):
            _fail(f"{side} vertex IDs and coordinates must be non-empty and aligned")
        if len(vertex_rows) > MAX_VERTICES_TOTAL or len(triangle_rows) > MAX_TRIANGLES_TOTAL:
            _fail(f"{side} embedded mesh exceeds geometry count caps")
        expected_object = _text(fragment.get("obj_object"), field=f"{side}.obj_object")
        obj_path = _source_path(fixture_root, fragment.get("mesh_file"), role=f"{side} mesh")
        obj_raw = _read_bounded(obj_path, role=f"{side} mesh", max_bytes=limits["max_bytes_each"])
        mesh_hash = _check_hash(obj_raw, fragment.get("mesh_sha256"), role=f"{side} mesh")
        mesh = _parse_obj(
            obj_raw,
            role=side,
            expected_object=expected_object,
            max_vertices=min(MAX_VERTICES_TOTAL, limits["max_vertices_total"] - total_vertices),
            max_triangles=min(MAX_TRIANGLES_TOTAL, limits["max_triangles_total"] - total_triangles),
        )
        if len(mesh["vertices"]) != len(vertex_rows) or len(mesh["faces"]) != len(triangle_rows):
            _fail(f"{side} OBJ counts differ from the embedded output piece")

        embedded_vertices = [_vec3(row, field=f"{side}.vertices_m[{index}]") for index, row in enumerate(vertex_rows)]
        for index, (actual, embedded) in enumerate(zip(mesh["vertices"], embedded_vertices)):
            if not _same_vec(actual, embedded):
                _fail(f"{side} OBJ vertex {index} differs from output vertices_m")

        embedded_faces: list[tuple[int, int, int]] = []
        for face_index, row in enumerate(triangle_rows):
            face_values = _list(row, field=f"{side}.triangles[{face_index}]")
            if len(face_values) != 3 or any(isinstance(index, bool) or not isinstance(index, int) for index in face_values):
                _fail(f"{side} triangle {face_index} must contain three integer vertex indices")
            face = tuple(face_values)
            if any(index < 0 or index >= len(mesh["vertices"]) for index in face) or len(set(face)) != 3:
                _fail(f"{side} triangle {face_index} has invalid vertex indices")
            embedded_faces.append(face)  # type: ignore[arg-type]
        if mesh["faces"] != embedded_faces:
            _fail(f"{side} OBJ face indices disagree with output triangles")

        for vertex_index, vertex_id in enumerate(vertex_ids):
            if vertex_id != f"{side}:v{vertex_index}":
                _fail(f"{side} output vertex ID {vertex_index} is not canonical")
        topology = _mesh_topology(mesh, role=f"{side} mesh")
        piece_metrics = _mapping(piece.get("metrics"), field=f"{side}.metrics")
        declared_volume = _finite_number(piece_metrics.get("volume_m3"), field=f"{side}.volume_m3")
        if not _near(topology["signed_volume_m3"], declared_volume):
            _fail(f"{side} actual signed volume differs from the embedded volume_m3")

        pieces_by_side[side] = piece
        meshes_by_side[side] = {**mesh, "topology": topology}
        mesh_hashes[side] = mesh_hash
        total_vertices += len(mesh["vertices"])
        total_triangles += len(mesh["faces"])

    if set(pieces_by_side) != {"negative", "positive"}:
        _fail("fixture output must provide negative and positive pieces")
    if total_vertices > min(MAX_VERTICES_TOTAL, limits["max_vertices_total"]):
        _fail("combined OBJ vertex count exceeds cap")
    if total_triangles > min(MAX_TRIANGLES_TOTAL, limits["max_triangles_total"]):
        _fail("combined OBJ triangle count exceeds cap")

    pair_surface_id = _text(ledger.get("surface_id"), field="pair_ledger.surface_id")
    cut_surface = _mapping(output.get("cut_surface"), field="output.cut_surface")
    if cut_surface.get("surface_id") != pair_surface_id:
        _fail("cut surface ID differs from the pair ledger")
    cut_vertices_rows = _list(cut_surface.get("vertices_m"), field="cut_surface.vertices_m")
    cut_vertices = [_vec3(row, field=f"cut_surface.vertices_m[{index}]") for index, row in enumerate(cut_vertices_rows)]
    cut_triangles_rows = _list(cut_surface.get("triangles"), field="cut_surface.triangles")
    cut_triangles: list[tuple[int, int, int]] = []
    for index, row in enumerate(cut_triangles_rows):
        items = _list(row, field=f"cut_surface.triangles[{index}]")
        if len(items) != 3 or any(isinstance(item, bool) or not isinstance(item, int) for item in items):
            _fail(f"cut surface triangle {index} is malformed")
        triangle = tuple(items)
        if any(vertex_index < 0 or vertex_index >= len(cut_vertices) for vertex_index in triangle):
            _fail(f"cut surface triangle {index} references a missing mother vertex")
        cut_triangles.append(triangle)  # type: ignore[arg-type]

    vertex_pairs = _list(ledger.get("vertex_pairs"), field="pair_ledger.vertex_pairs")
    input_vertex_map: dict[int, dict[str, Any]] = {}
    for pair_index, pair_raw in enumerate(vertex_pairs):
        pair = _mapping(pair_raw, field=f"vertex_pairs[{pair_index}]")
        input_index = pair.get("input_surface_vertex_index")
        if isinstance(input_index, bool) or not isinstance(input_index, int) or input_index in input_vertex_map:
            _fail("vertex pair input surface indices must be unique integers")
        if input_index < 0 or input_index >= len(cut_vertices):
            _fail("vertex pair references a missing input surface vertex")
        if not _same_vec(_vec3(pair.get("mother_coordinate_m"), field="vertex pair mother_coordinate_m"), cut_vertices[input_index]):
            _fail("vertex pair canonical mother coordinate differs from cut surface")
        input_vertex_map[input_index] = pair
        for side in ("negative", "positive"):
            local_index = pair.get(f"{side}_local_vertex_index")
            if isinstance(local_index, bool) or not isinstance(local_index, int):
                _fail(f"vertex pair {side} local index must be an integer")
            piece = pieces_by_side[side]
            mesh = meshes_by_side[side]
            vertex_id = pair.get(f"{side}_vertex_id")
            if local_index < 0 or local_index >= len(mesh["vertices"]):
                _fail(f"vertex pair {side} local index is out of range")
            if vertex_id != piece["vertex_ids"][local_index]:
                _fail(f"vertex pair {side} ID does not bind to actual mesh vertex index")
            if not _same_vec(mesh["vertices"][local_index], cut_vertices[input_index]):
                _fail(f"vertex pair {side} actual mesh coordinate differs from mother coordinate")
    if set(input_vertex_map) != set(range(len(cut_vertices))):
        _fail("pair ledger must bind every input surface vertex exactly once")

    reciprocal_pairs = _list(ledger.get("reciprocal_vertex_pairs"), field="pair_ledger.reciprocal_vertex_pairs")
    if len(reciprocal_pairs) != len(vertex_pairs):
        _fail("reciprocal vertex pair count differs from forward map")
    reciprocal_seen: set[tuple[str, str]] = set()
    for reciprocal_index, reciprocal_raw in enumerate(reciprocal_pairs):
        reciprocal = _mapping(reciprocal_raw, field=f"reciprocal_vertex_pairs[{reciprocal_index}]")
        ids = (reciprocal.get("negative_vertex_id"), reciprocal.get("positive_vertex_id"))
        reciprocal_seen.add((str(ids[0]), str(ids[1])))
        coordinate = _vec3(reciprocal.get("mother_coordinate_m"), field="reciprocal mother_coordinate_m")
        found = False
        for pair in vertex_pairs:
            if ids == (pair.get("positive_vertex_id"), pair.get("negative_vertex_id")):
                found = _same_vec(coordinate, _vec3(pair.get("mother_coordinate_m"), field="pair mother_coordinate_m"))
                break
        if not found:
            _fail("reciprocal vertex pair is not the reversed forward map at the same mother coordinate")
    if len(reciprocal_seen) != len(vertex_pairs):
        _fail("reciprocal vertex pairs contain duplicate mappings")

    triangle_pairs = _list(ledger.get("triangle_pairs"), field="pair_ledger.triangle_pairs")
    if len(triangle_pairs) != len(cut_triangles):
        _fail("triangle pair count differs from explicit cut surface triangles")
    pair_reports: list[dict[str, Any]] = []
    paired_area = 0.0
    for pair_index, pair_raw in enumerate(triangle_pairs):
        pair = _mapping(pair_raw, field=f"triangle_pairs[{pair_index}]")
        input_triangle_index = pair.get("input_surface_triangle_index")
        if isinstance(input_triangle_index, bool) or not isinstance(input_triangle_index, int) or input_triangle_index != pair_index:
            _fail("triangle pair indices must bind each cut surface triangle in order")
        if input_triangle_index < 0 or input_triangle_index >= len(cut_triangles):
            _fail("triangle pair references a missing input surface triangle")
        expected_mother_indices = cut_triangles[input_triangle_index]
        expected_mother_points = tuple(cut_vertices[index] for index in expected_mother_indices)
        ledger_mother_points = tuple(
            _vec3(row, field=f"triangle pair {pair_index} mother coordinate")
            for row in _list(pair.get("mother_coordinates_m"), field="triangle_pair.mother_coordinates_m")
        )
        if len(ledger_mother_points) != 3 or ledger_mother_points != expected_mother_points:
            _fail("triangle pair mother coordinates differ from the explicit surface triangle")

        expected_local: dict[str, tuple[int, int, int]] = {}
        for side in ("negative", "positive"):
            indices = tuple(input_vertex_map[index][f"{side}_local_vertex_index"] for index in expected_mother_indices)
            expected_local[side] = indices  # type: ignore[assignment]
            face_index = pair.get(f"{side}_face_index")
            if isinstance(face_index, bool) or not isinstance(face_index, int):
                _fail(f"triangle pair {side}_face_index must be an integer")
            faces: list[tuple[int, int, int]] = meshes_by_side[side]["faces"]
            if face_index < 0 or face_index >= len(faces):
                _fail(f"triangle pair references missing {side} face {face_index}")
            expected_face_id = f"{pieces_by_side[side]['piece_id']}:f{face_index}"
            if pair.get(f"{side}_face_id") != expected_face_id:
                _fail(f"triangle pair {side} face ID does not match the referenced output face index")
            if faces[face_index] != tuple(_list(pieces_by_side[side]["triangles"][face_index], field="piece triangle")):
                _fail(f"triangle pair {side} face is not bound to actual OBJ face index")
            lineage = _list(pieces_by_side[side].get("face_lineage"), field=f"{side}.face_lineage")
            if face_index >= len(lineage):
                _fail(f"triangle pair {side} face has no lineage entry")
            lineage_entry = _mapping(lineage[face_index], field=f"{side}.face_lineage[{face_index}]")
            if (
                lineage_entry.get("kind") != "explicit_fracture_face_fixture"
                or lineage_entry.get("surface_id") != pair_surface_id
                or lineage_entry.get("pair_id") != pair.get("pair_id")
            ):
                _fail(f"triangle pair {side} face lineage does not bind to this pair ID")
            expected_group = _text(fragments_by_side[side].get("fracture_face_group"), field=f"{side}.fracture_face_group")
            if meshes_by_side[side]["groups"][face_index] != expected_group:
                _fail(f"triangle pair {side} face OBJ group does not bind to the fracture surface")

        negative_face = meshes_by_side["negative"]["faces"][pair["negative_face_index"]]
        positive_face = meshes_by_side["positive"]["faces"][pair["positive_face_index"]]
        if negative_face != expected_local["negative"]:
            _fail("negative paired face vertices do not follow canonical mother-coordinate order")
        if not _cyclic_reverse_equal(positive_face, expected_local["positive"]):
            _fail("positive paired face does not have reverse winding against the canonical mother order")

        actual_points: dict[str, tuple[tuple[float, float, float], ...]] = {}
        for side, face in (("negative", negative_face), ("positive", positive_face)):
            actual_points[side] = tuple(meshes_by_side[side]["vertices"][index] for index in face)
        if actual_points["negative"] != expected_mother_points:
            _fail("negative paired face actual OBJ coordinates differ from canonical mother coordinates")
        if set(actual_points["positive"]) != set(expected_mother_points):
            _fail("positive paired face actual OBJ coordinates differ from canonical mother coordinates")

        negative_normal = meshes_by_side["negative"]["topology"]["triangle_normals"][pair["negative_face_index"]]
        positive_normal = meshes_by_side["positive"]["topology"]["triangle_normals"][pair["positive_face_index"]]
        computed_dot = _dot(negative_normal, positive_normal)
        if computed_dot >= -1.0 + 1.0e-10:
            _fail("actual paired outward normals are not opposite")
        expected_normal = _unit(
            _vec3(_mapping(cut_surface.get("plane"), field="cut_surface.plane").get("normal"), field="cut_surface.plane.normal"),
            field="cut surface plane normal",
        )
        if _dot(negative_normal, expected_normal) < 1.0 - 1.0e-10:
            _fail("negative actual paired normal disagrees with cut surface outward orientation")
        declared_negative = _unit(_vec3(pair.get("negative_outward_normal"), field="negative_outward_normal"), field="negative_outward_normal")
        declared_positive = _unit(_vec3(pair.get("positive_outward_normal"), field="positive_outward_normal"), field="positive_outward_normal")
        if _dot(negative_normal, declared_negative) < 1.0 - 1.0e-10 or _dot(positive_normal, declared_positive) < 1.0 - 1.0e-10:
            _fail("actual paired normals differ from ledger normals")

        negative_area = meshes_by_side["negative"]["topology"]["triangle_areas_m2"][pair["negative_face_index"]]
        positive_area = meshes_by_side["positive"]["topology"]["triangle_areas_m2"][pair["positive_face_index"]]
        if not _near(negative_area, positive_area):
            _fail("paired actual triangles have different areas")
        paired_area += negative_area
        pair_reports.append(
            {
                "pair_id": pair.get("pair_id"),
                "negative_face_id": pair["negative_face_id"],
                "positive_face_id": pair["positive_face_id"],
                "negative_face_index": pair["negative_face_index"],
                "positive_face_index": pair["positive_face_index"],
                "area_m2": negative_area,
                "actual_outward_normal_dot": computed_dot,
            }
        )

    declared_area = _finite_number(cut_surface.get("area_m2"), field="cut_surface.area_m2")
    if not _near(paired_area, declared_area):
        _fail("actual paired area differs from declared cut surface area_m2")

    original_body = _mapping(output.get("original_body"), field="output.original_body")
    original_volume = _finite_number(original_body.get("volume_m3"), field="original_body.volume_m3")
    fragment_volume_sum = sum(meshes_by_side[side]["topology"]["signed_volume_m3"] for side in ("negative", "positive"))
    volume_error = fragment_volume_sum - original_volume
    if not _near(fragment_volume_sum, original_volume):
        _fail("actual fragment signed volumes do not numerically conserve original volume")

    ignored_output_checks = output.get("checks", {})
    ignored_piece_topology = [piece.get("topology", {}) for piece in pieces]
    report = {
        "schema_version": "pgfrag.mesh-binding-receipt.v1",
        "validation_status": "PASS",
        "qualification": "FIXTURE_GEOMETRY_ONLY",
        "mechanics_qualification": "NOT_COUPLED",
        "physical_pack_count": 0,
        "physical_fracture_claim": "NOT_CLAIMED",
        "coordinate_unit": "m",
        "source_assets": {
            "fixture_output_file": manifest["output_file"]["file"],
            "fixture_output_sha256": output_hash,
            "pair_ledger_file": manifest["pair_ledger_file"]["file"],
            "pair_ledger_sha256": ledger_hash,
            "fragment_mesh_files": {
                side: {
                    "file": fragments_by_side[side]["mesh_file"],
                    "sha256": mesh_hashes[side],
                    "vertex_count": len(meshes_by_side[side]["vertices"]),
                    "triangle_count": len(meshes_by_side[side]["faces"]),
                }
                for side in ("negative", "positive")
            },
            "separate_fragment_meshes_json_present": False,
        },
        "geometry_checks": {
            "mesh_bytes_match_pinned_sha256": True,
            "mesh_vertices_match_embedded_output_coordinates": True,
            "mesh_faces_match_embedded_output_indices": True,
            "pair_faces_bind_to_actual_fragment_face_ids": True,
            "pair_vertices_bind_to_actual_mesh_indices_and_mother_coordinates": True,
            "paired_face_winding_is_reciprocal": True,
            "all_mesh_edges_are_closed_with_opposite_directed_incidence": True,
            "both_meshes_are_outward_oriented_by_positive_signed_volume": True,
            "paired_actual_normals_are_opposed": True,
            "paired_actual_area_matches_cut_surface_area": True,
            "fragment_signed_volumes_match_declared_values": True,
            "fragment_volume_sum_matches_original_numerically": True,
        },
        "fragments": {
            side: {
                "piece_id": pieces_by_side[side]["piece_id"],
                "edge_count": meshes_by_side[side]["topology"]["edge_count"],
                "closed_oriented": meshes_by_side[side]["topology"]["closed_oriented"],
                "signed_volume_m3": meshes_by_side[side]["topology"]["signed_volume_m3"],
            }
            for side in ("negative", "positive")
        },
        "paired_fracture_faces": pair_reports,
        "interface_area_m2": paired_area,
        "volume_conservation": {
            "original_volume_m3": original_volume,
            "fragment_volume_sum_m3": fragment_volume_sum,
            "error_m3": volume_error,
            "status": "NUMERICALLY_CONSERVED_WITHIN_TOLERANCE",
        },
        "partition_geometry_status": "VOLUME_SUM_CONSERVED_ONLY_NO_OVERLAP_OR_FULL_COVERAGE_PROOF",
        "no_overlap_proof": "NOT_ESTABLISHED",
        "full_coverage_proof": "NOT_ESTABLISHED",
        "intersection_test": "NOT_RUN",
        "producer_boolean_claims_ignored": {
            "output_checks_object_ignored": bool(ignored_output_checks),
            "piece_topology_objects_ignored": bool(ignored_piece_topology),
            "reason": "This receipt is based on parsed OBJ bytes and numeric fields; producer booleans cannot establish topology or physical qualification.",
        },
        "a12_use": {
            "observation_only": True,
            "mother_geometry_or_clean_gt_available_to_inference": False,
            "pretraining_qualified": False,
        },
    }
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, help="hash-pinned planar fixture binding manifest")
    parser.add_argument("--output", help="write the passing receipt JSON to this path")
    args = parser.parse_args(argv)
    try:
        report = validate_binding(args.manifest)
    except BindingValidationError as exc:
        print(json.dumps({"validation_status": "REJECTED", "reason": str(exc)}, indent=2), file=sys.stderr)
        return 2
    serialized = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output:
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized, encoding="utf-8")
    else:
        print(serialized, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
