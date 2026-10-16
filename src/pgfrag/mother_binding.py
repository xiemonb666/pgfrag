"""Read-only binding between a VF source asset and a verified tetra partition.

The binding carries stable exterior source-face IDs through mesh serialization.
It is geometry/data plumbing only: it does not create a fracture, observations,
physical qualification, or training eligibility. Mesh rows are joined by
original Gmsh node/cell tags and ordered tetrahedron corners, never by row
position or coordinate proximity. Coordinate comparison is only a serialization
consistency check.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import zipfile
from types import MappingProxyType
from typing import Any, Mapping

import numpy as np

from .cohesive_mesh import read_gmsh22


SCHEMA = "pgfrag.native_mother_binding.v1"
MAX_JSON_BYTES = 16 * 1024 * 1024
MAX_NPZ_BYTES = 64 * 1024 * 1024
MAX_NPZ_UNCOMPRESSED_BYTES = 256 * 1024 * 1024
MAX_MSH_BYTES = 64 * 1024 * 1024
MAX_TETRAHEDRA = 1_000_000
MAX_NODES = 2_000_000
MAX_TRIANGLES = 4_000_000
MAX_OUTPUT_BYTES = 32 * 1024 * 1024
COORDINATE_ROUNDTRIP_RTOL = 8e-16
COORDINATE_ROUNDTRIP_ATOL_M = 1e-30

_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_POSITIVE_CLAIMS = frozenset({
    "physical_pass", "production_physical_pass", "pretrain_eligible",
    "training_eligible", "deployment_eligible", "solver_qualified",
    "fracture_qualified", "observation_training_eligible",
})
_QUALIFICATION_STATUS_KEYS = frozenset({
    "physical_status", "physical_material_status", "material_qualification_status",
    "calibration_status", "solver_qualification", "fracture_or_fragment_result",
    "training_status", "observation_status",
})
_POSITIVE_STATUS_PREFIXES = ("PASS", "QUALIFIED", "CALIBRATED", "VALIDATED",
                             "ELIGIBLE", "APPROVED", "READY")


class MotherBindingError(ValueError):
    """Raised when source identity, ancestry, or qualification cannot be bound."""


@dataclass(frozen=True, slots=True)
class MotherBinding:
    """Immutable source-identity summary and tag-based mapping arrays.

    Arrays use original Gmsh tags rather than native row numbers. The only row
    references retained are explicit NPZ storage locators such as `unused_node_rows`.
    """

    binding_id: str
    mother_id: str
    summary: Mapping[str, Any]
    source_face_ids: tuple[str, ...]
    source_face_roles: tuple[str, ...]
    node_tags: np.ndarray
    unused_node_rows: np.ndarray
    cell_tags: np.ndarray
    cell_corner_node_tags: np.ndarray
    exterior_triangle_node_tags: np.ndarray
    exterior_source_face_index: np.ndarray
    interface_triangle_node_tags: np.ndarray
    interface_parent_cell_tags: np.ndarray


def _fail(message: str, *, field: str | None = None) -> None:
    suffix = f" [{field}]" if field else ""
    raise MotherBindingError(message + suffix)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _bounded_file(path: Path, limit: int, field: str) -> Path:
    value = Path(path)
    if not value.is_file():
        _fail("required input file is missing", field=field)
    size = value.stat().st_size
    if size < 0 or size > limit:
        _fail("input file exceeds its byte limit", field=field)
    return value


def _read_json(path: Path, field: str) -> dict[str, Any]:
    value = _bounded_file(path, MAX_JSON_BYTES, field)
    try:
        raw = value.read_text(encoding="utf-8")
        result = json.loads(raw, parse_constant=lambda token: _fail(
            f"non-finite JSON constant {token}", field=field))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        _fail(f"cannot read valid UTF-8 JSON: {exc}", field=field)
    if not isinstance(result, dict):
        _fail("JSON root must be an object", field=field)
    return result


def _is_digest(value: Any) -> bool:
    return isinstance(value, str) and _SHA256_RE.fullmatch(value) is not None


def _safe_relative(root: Path, relative: str, field: str) -> Path:
    pure = PurePosixPath(str(relative).replace("\\", "/"))
    if pure.is_absolute() or not pure.parts or any(part in ("", ".", "..") for part in pure.parts):
        _fail("artifact path must be a safe relative path", field=field)
    candidate = root.joinpath(*pure.parts).resolve()
    try:
        candidate.relative_to(root.resolve())
    except ValueError:
        _fail("artifact path escapes its source root", field=field)
    return candidate


def _canonical_bytes(value: Any) -> bytes:
    try:
        return json.dumps(value, sort_keys=True, separators=(",", ":"),
                          ensure_ascii=False, allow_nan=False).encode("utf-8")
    except (TypeError, ValueError) as exc:
        _fail(f"value is not canonical JSON: {exc}")


def _deep_freeze(value: Any) -> Any:
    if isinstance(value, Mapping):
        return MappingProxyType({str(key): _deep_freeze(item) for key, item in value.items()})
    if isinstance(value, list):
        return tuple(_deep_freeze(item) for item in value)
    return value


def _deep_thaw(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _deep_thaw(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_deep_thaw(item) for item in value]
    return value


def _reject_promoted_claims(value: Any, field: str = "input") -> None:
    """Reject explicit positive physical/solver/fracture/training claims.

    Structural CAD `PASS` and candidate status labels are intentionally outside
    this list. They are preserved as separate source fields in the output.
    """
    if isinstance(value, Mapping):
        for key, item in value.items():
            normalized = str(key).strip().lower()
            if normalized in _POSITIVE_CLAIMS and item is True:
                _fail("geometry-only binding cannot accept a positive physical or eligibility claim",
                      field=f"{field}.{key}")
            if normalized in _QUALIFICATION_STATUS_KEYS and isinstance(item, str):
                status = item.strip().upper()
                if status.startswith(_POSITIVE_STATUS_PREFIXES):
                    _fail("geometry-only binding cannot accept promoted qualification status",
                          field=f"{field}.{key}")
            _reject_promoted_claims(item, f"{field}.{key}")
    elif isinstance(value, list):
        for index, item in enumerate(value):
            _reject_promoted_claims(item, f"{field}[{index}]")


def _verify_receipt_files(result_path: Path, result: Mapping[str, Any],
                          audit_path: Path, audit: Mapping[str, Any]
                          ) -> tuple[Path, Path, dict[str, dict[str, Any]]]:
    if result.get("schema") != "vesselforge.generic_step_tetra_partition.v1":
        _fail("partition result schema is unsupported", field="result.schema")
    if result.get("status") != "GEOMETRY_DIAGNOSTIC_ONLY":
        _fail("partition result is not marked geometry-diagnostic-only", field="result.status")
    if audit.get("schema") != "vesselforge.reference_v6_witness_tetra_run.v1" \
            or audit.get("status") != "PASS_WITNESS_ANCESTRY_TETRA_MESH":
        _fail("audit receipt does not contain the expected completed ancestry audit",
              field="audit_receipt.schema/status")
    if audit.get("result_sha256") != _sha256(result_path) \
            or audit.get("result_json_bytes") != result_path.stat().st_size:
        _fail("partition result does not match its trusted audit receipt",
              field="audit_receipt.result_sha256/result_json_bytes")
    files = result.get("files")
    if not isinstance(files, Mapping) or not {"mesh.msh", "mesh.npz"}.issubset(files):
        _fail("partition result must receipt-bind mesh.msh and mesh.npz", field="result.files")
    verified: dict[str, dict[str, Any]] = {}
    paths: dict[str, Path] = {}
    for name in ("mesh.msh", "mesh.npz"):
        record = files[name]
        if not isinstance(record, Mapping):
            _fail("partition file receipt row must be an object", field=f"result.files.{name}")
        digest, size = record.get("sha256"), record.get("bytes")
        if not _is_digest(digest) or isinstance(size, bool) or not isinstance(size, int) or size < 0:
            _fail("partition file receipt requires SHA-256 and byte count",
                  field=f"result.files.{name}")
        path = _bounded_file(result_path.parent / name,
                             MAX_MSH_BYTES if name.endswith(".msh") else MAX_NPZ_BYTES,
                             f"partition.{name}")
        if path.stat().st_size != size or _sha256(path) != digest:
            _fail("partition file differs from its result receipt", field=f"result.files.{name}")
        verified[name] = {"sha256": digest, "bytes": size}
        paths[name] = path
    return paths["mesh.msh"], paths["mesh.npz"], verified


def _verify_source_hash_chain(
    source_root: Path,
    canonical_step_path: Path,
    native_face_identity_path: Path,
    witness_manifest_path: Path,
    v15_receipt_path: Path,
    label_import_receipt_path: Path,
    canonical_step_input_key: str,
    audit: Mapping[str, Any],
) -> tuple[dict[str, Path], dict[str, dict[str, Any]], dict[str, Any],
           dict[str, Any], dict[str, Any], dict[str, Any], dict[str, Any], dict[str, Any]]:
    source_paths = {
        "solid.step": _safe_relative(source_root, "geometry/solid.step", "source_asset"),
        "section.json": _safe_relative(source_root, "geometry/section.json", "source_asset"),
        "face_lineage.json": _safe_relative(source_root, "geometry/face_lineage.json", "source_asset"),
        "geometry_report.json": _safe_relative(source_root, "geometry_report.json", "source_asset"),
        "native_face_identity.json": Path(native_face_identity_path),
        "witness_manifest.json": Path(witness_manifest_path),
        "v15_receipt.json": Path(v15_receipt_path),
        "label_import_receipt.json": Path(label_import_receipt_path),
        str(canonical_step_input_key): Path(canonical_step_path),
    }
    if len(source_paths) != 9 or len(set(source_paths)) != len(source_paths):
        _fail("canonical step input key collides with another source artifact", field="source paths")
    input_hashes = audit.get("input_hashes")
    if not isinstance(input_hashes, Mapping):
        _fail("audit receipt is missing input hashes", field="audit_receipt.input_hashes")
    for name, path in source_paths.items():
        expected = input_hashes.get(name)
        bounded = _bounded_file(path, MAX_JSON_BYTES if name.endswith(".json") else MAX_NPZ_BYTES,
                                f"source.{name}")
        actual = _sha256(bounded)
        if not _is_digest(expected) or actual != expected:
            _fail("source file hash does not match the trusted audit receipt",
                  field=f"audit_receipt.input_hashes.{name}")
    # Run receipts may list additional build inputs, but each listed source
    # artifact must also agree with the independently nested lineage record.
    nested = audit.get("source_geometry_and_lineage")
    nested_hashes = nested.get("asset_hashes") if isinstance(nested, Mapping) else None
    if not isinstance(nested_hashes, Mapping):
        _fail("audit receipt is missing nested source geometry hashes",
              field="audit_receipt.source_geometry_and_lineage.asset_hashes")
    for name in ("solid.step", "section.json", "face_lineage.json", "geometry_report.json",
                 "native_face_identity.json", "witness_manifest.json", "v15_receipt.json",
                 "label_import_receipt.json", canonical_step_input_key):
        if nested_hashes.get(name) != input_hashes.get(name):
            _fail("nested source hash disagrees with audit input hash",
                  field=f"source_geometry_and_lineage.asset_hashes.{name}")

    source_values = {
        "geometry_report": _read_json(source_paths["geometry_report.json"], "geometry_report.json"),
        "section": _read_json(source_paths["section.json"], "section.json"),
        "lineage": _read_json(source_paths["face_lineage.json"], "face_lineage.json"),
        "identity": _read_json(source_paths["native_face_identity.json"], "native_face_identity.json"),
        "witness": _read_json(source_paths["witness_manifest.json"], "witness_manifest.json"),
        "v15": _read_json(source_paths["v15_receipt.json"], "v15_receipt.json"),
        "label": _read_json(source_paths["label_import_receipt.json"], "label_import_receipt.json"),
    }
    report = source_values["geometry_report"]
    artifact_hashes = report.get("measurements", {}).get("artifact_sha256", {})
    if not isinstance(artifact_hashes, Mapping) or not artifact_hashes:
        _fail("source geometry report has no artifact hash table",
              field="geometry_report.measurements.artifact_sha256")
    optional_missing_artifacts: list[str] = []
    required_geometry_artifacts = {"solid.step", "section.json", "face_lineage.json", "surface.npz"}
    for relative, expected in artifact_hashes.items():
        if not _is_digest(expected):
            _fail("source geometry artifact hash is invalid", field=f"geometry_report.{relative}")
        path = _safe_relative(source_root / "geometry", str(relative),
                              "geometry_report.measurements.artifact_sha256")
        if not path.is_file() and str(relative) not in required_geometry_artifacts:
            optional_missing_artifacts.append(str(relative))
            continue
        bounded = _bounded_file(path, MAX_NPZ_BYTES, f"source_asset.{relative}")
        if _sha256(bounded) != expected:
            _fail("source asset file differs from its geometry report", field=f"source_asset.{relative}")
    surface_path = _safe_relative(source_root, "geometry/surface.npz", "source_asset")
    if not surface_path.is_file():
        _fail("source asset surface mesh is required", field="source_asset.geometry.surface.npz")
    source_paths["surface.npz"] = surface_path
    provenance_path = _safe_relative(source_root, "provenance.json", "source_asset")
    source_paths["provenance.json"] = provenance_path
    source_values["optional_missing_artifacts"] = optional_missing_artifacts

    expected_source_hashes = {
        "solid.step": _sha256(source_paths["solid.step"]),
        "section.json": _sha256(source_paths["section.json"]),
        "face_lineage.json": _sha256(source_paths["face_lineage.json"]),
        "geometry_report.json": _sha256(source_paths["geometry_report.json"]),
    }
    identity = source_values["identity"]
    witness = source_values["witness"]
    v15 = source_values["v15"]
    label = source_values["label"]
    if identity.get("schema") != "vesselforge.daily_native_face_identity.v1" \
            or identity.get("source_hashes") != expected_source_hashes:
        _fail("native face identity is not bound to this source asset", field="native_face_identity")
    if v15.get("schema") != "vesselforge.daily_native_face_witness_probe.v1" \
            or v15.get("status") != "PASS_CAD_PARENT_MAP_WITHIN_DECLARED_OCC_TOLERANCES" \
            or v15.get("source_hashes") != expected_source_hashes:
        _fail("v15 receipt does not verify the source CAD ancestry chain", field="v15_receipt")
    if label.get("schema") != "vesselforge.daily_named_step_label_probe.v1" \
            or label.get("source_hashes") != expected_source_hashes:
        _fail("label import receipt does not bind this source asset", field="label_import_receipt")
    canonical_hash = _sha256(source_paths[canonical_step_input_key])
    if witness.get("schema") != "vesselforge.native_face_witness_manifest.v1" \
            or witness.get("mapping_status") != "EXACT_NATIVE_FACE_WITNESSES" \
            or witness.get("source_step_sha256") != expected_source_hashes["solid.step"] \
            or witness.get("source_face_lineage_sha256") != expected_source_hashes["face_lineage.json"] \
            or witness.get("native_body_step_sha256") != canonical_hash:
        _fail("witness manifest is not bound to the source and canonical STEP",
              field="witness_manifest")
    if v15.get("native_face_identity_sha256") != _sha256(source_paths["native_face_identity.json"]) \
            or v15.get("native_face_witness_manifest_sha256") != _sha256(source_paths["witness_manifest.json"]) \
            or v15.get("native_body_step_sha256") != canonical_hash:
        _fail("v15 receipt hash bindings disagree", field="v15_receipt.hash_bindings")
    if label.get("native_body_step_sha256") != canonical_hash:
        _fail("label import receipt is bound to another canonical STEP",
              field="label_import_receipt.native_body_step_sha256")
    return source_paths, source_values, report, identity, witness, v15, label, expected_source_hashes


def _face_table(lineage: Mapping[str, Any], identity: Mapping[str, Any],
                witness: Mapping[str, Any], witness_directory: Path
                ) -> tuple[list[dict[str, Any]], dict[str, int]]:
    source_faces = lineage.get("faces")
    identities = identity.get("faces")
    witnesses = witness.get("witnesses")
    if not isinstance(source_faces, Mapping) or not isinstance(identities, list) \
            or not isinstance(witnesses, list) or not source_faces:
        _fail("source face lineage, native identity, and witness rows are required")
    if len(source_faces) != len(identities) or len(source_faces) != len(witnesses) \
            or len(source_faces) > 4096:
        _fail("source face identity tables have different or excessive row counts")
    source_by_index: dict[int, dict[str, Any]] = {}
    for key, value in source_faces.items():
        try:
            index = int(key)
        except (TypeError, ValueError):
            _fail("source face lineage keys must be integer indices", field="face_lineage.faces")
        if index < 0 or str(index) != str(key) or index in source_by_index or not isinstance(value, Mapping):
            _fail("source face lineage indices are invalid or duplicated", field="face_lineage.faces")
        role = value.get("role")
        segments = value.get("source_segment_indices")
        if not isinstance(role, str) or not role.strip() or not isinstance(segments, list) \
                or not segments or any(isinstance(v, bool) or not isinstance(v, int) or v < 0
                                       for v in segments) or len(set(segments)) != len(segments):
            _fail("source face role and segment lineage are invalid", field=f"face_lineage.faces.{key}")
        source_by_index[index] = {"role": role, "source_segment_indices": sorted(segments)}

    by_fid: dict[str, dict[str, Any]] = {}
    by_native_index: dict[int, str] = {}
    for row in identities:
        if not isinstance(row, Mapping):
            _fail("native face identity row must be an object", field="native_face_identity.faces")
        fid = row.get("source_original_brep_face_id")
        index = row.get("source_original_brep_face_index")
        native_index = row.get("native_face_index")
        if not isinstance(fid, str) or not fid.strip() or isinstance(index, bool) \
                or not isinstance(index, int) or isinstance(native_index, bool) \
                or not isinstance(native_index, int) or index not in source_by_index \
                or fid in by_fid or native_index in by_native_index:
            _fail("native face identity contains missing or duplicate IDs", field="native_face_identity.faces")
        role = row.get("role")
        segments = row.get("source_segment_indices")
        expected = source_by_index[index]
        if role != expected["role"] or not isinstance(segments, list) \
                or sorted(segments) != expected["source_segment_indices"]:
            _fail("native identity role/segments conflict with source lineage",
                  field=f"native_face_identity.faces.{fid}")
        by_fid[fid] = {"source_face_id": fid, "source_face_index": index,
                        "native_face_index": native_index, "role": role,
                        "source_segment_indices": expected["source_segment_indices"]}
        by_native_index[native_index] = fid
    if set(row["source_face_index"] for row in by_fid.values()) != set(source_by_index):
        _fail("native identity table does not cover every source face index")

    seen_witness: set[str] = set()
    witness_root = Path(witness_directory).resolve()
    if not witness_root.is_dir():
        _fail("witness directory is missing", field="witness_directory")
    for row in witnesses:
        if not isinstance(row, Mapping):
            _fail("witness row must be an object", field="witness_manifest.witnesses")
        fid = row.get("source_face_id")
        record = by_fid.get(str(fid))
        digest, size = row.get("witness_sha256"), row.get("witness_bytes")
        path_value = row.get("witness_path")
        if record is None or fid in seen_witness or not _is_digest(digest) \
                or isinstance(size, bool) or not isinstance(size, int) or size < 1 \
                or row.get("native_face_child_verified") is not True \
                or row.get("witness_face_entity_count") != 1 \
                or row.get("role") != record["role"] or not isinstance(path_value, str):
            _fail("witness row is incomplete, conflicting, or ambiguous",
                  field=f"witness_manifest.witnesses.{fid}")
        witness_name = Path(path_value).name
        if witness_name in ("", ".", "..") or Path(witness_name).name != witness_name:
            _fail("witness filename is invalid", field=f"witness_manifest.witnesses.{fid}")
        path = _bounded_file(witness_root / witness_name, MAX_NPZ_BYTES,
                             f"witness.{fid}")
        if path.stat().st_size != size or _sha256(path) != digest:
            _fail("witness file does not match its manifest", field=f"witness.{fid}")
        record["witness_sha256"] = digest
        record["witness_bytes"] = size
        seen_witness.add(str(fid))
    if seen_witness != set(by_fid):
        _fail("witness manifest does not cover exactly the source face identity table")
    ordered = sorted(by_fid.values(), key=lambda row: row["source_face_id"])
    return ordered, {row["source_face_id"]: i for i, row in enumerate(ordered)}


def _load_npz(path: Path) -> dict[str, np.ndarray]:
    bounded = _bounded_file(path, MAX_NPZ_BYTES, "mesh.npz")
    try:
        with zipfile.ZipFile(bounded) as archive:
            entries = archive.infolist()
            if len(entries) > 64 or sum(entry.file_size for entry in entries) > MAX_NPZ_UNCOMPRESSED_BYTES:
                _fail("NPZ has excessive members or uncompressed size", field="mesh.npz")
            names = [entry.filename for entry in entries]
            if len(names) != len(set(names)) or any(
                    PurePosixPath(name).is_absolute() or len(PurePosixPath(name).parts) != 1
                    or not name.endswith(".npy") for name in names):
                _fail("NPZ member names are unsafe or duplicated", field="mesh.npz")
        with np.load(bounded, allow_pickle=False) as data:
            arrays = {name: np.array(data[name], copy=True) for name in data.files}
    except (OSError, ValueError, EOFError, zipfile.BadZipFile) as exc:
        _fail(f"cannot load safe NPZ: {exc}", field="mesh.npz")
    if any(value.dtype.hasobject for value in arrays.values()):
        _fail("object arrays are forbidden", field="mesh.npz")
    required = {
        "nodes_m", "tets", "original_node_tags", "original_cell_tags",
        "source_partition_volume", "volume_m3", "interface_triangles",
        "interface_parents", "interface_parent_partition", "exterior_triangles",
        "exterior_surface_tags", "exterior_source_face_ids", "exterior_roles",
        "exterior_ancestry_exact", "partition_entity_tags",
    }
    if not required.issubset(arrays):
        _fail("partition NPZ is missing required ancestry arrays", field="mesh.npz")
    return arrays


def _int_array(arrays: Mapping[str, np.ndarray], key: str, shape: tuple[int, ...] | None = None
               ) -> np.ndarray:
    value = np.asarray(arrays[key])
    if value.dtype.kind not in "iu" or value.dtype.kind == "b" \
            or (shape is not None and value.shape != shape):
        _fail("array dtype or shape violates the integer ancestry contract", field=f"mesh.npz.{key}")
    return value.astype(np.int64, copy=False)


def _validate_mesh_arrays(arrays: Mapping[str, np.ndarray], mesh_path: Path,
                          face_rows: list[dict[str, Any]], fid_to_index: Mapping[str, int],
                          exterior_physical_name: str,
                          interface_physical_name: str) -> dict[str, np.ndarray | int]:
    nodes = np.asarray(arrays["nodes_m"])
    if nodes.dtype.kind != "f" or nodes.dtype.itemsize != 8 or nodes.ndim != 2 \
            or nodes.shape[1] != 3 or not 1 <= len(nodes) <= MAX_NODES \
            or not np.isfinite(nodes).all():
        _fail("nodes_m must be a bounded finite float64 (N,3) array", field="mesh.npz.nodes_m")
    tets = _int_array(arrays, "tets")
    n_tets = len(tets)
    if tets.ndim != 2 or tets.shape[1] != 4 or n_tets < 1 or n_tets > MAX_TETRAHEDRA \
            or np.any(tets < 0) or np.any(tets >= len(nodes)):
        _fail("tets must be bounded valid T4 row connectivity", field="mesh.npz.tets")
    node_tags = _int_array(arrays, "original_node_tags", (len(nodes),))
    cell_tags = _int_array(arrays, "original_cell_tags", (n_tets,))
    if np.any(node_tags <= 0) or len(np.unique(node_tags)) != len(node_tags):
        _fail("original node tags must be positive and unique", field="mesh.npz.original_node_tags")
    if np.any(cell_tags <= 0) or len(np.unique(cell_tags)) != len(cell_tags):
        _fail("original cell tags must be positive and unique", field="mesh.npz.original_cell_tags")
    if np.any(np.sort(tets, axis=1)[:, 1:] == np.sort(tets, axis=1)[:, :-1]):
        _fail("T4 rows must contain four distinct nodes", field="mesh.npz.tets")
    partitions = _int_array(arrays, "source_partition_volume", (n_tets,))
    volumes = np.asarray(arrays["volume_m3"])
    if volumes.dtype.kind != "f" or volumes.shape != (n_tets,) or not np.isfinite(volumes).all() \
            or np.any(volumes <= 0.0):
        _fail("volume_m3 must contain positive finite rows", field="mesh.npz.volume_m3")

    face_incidence: dict[tuple[int, int, int], list[int]] = {}
    local_faces = ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))
    for tet_row, tet in enumerate(tets):
        tags = node_tags[tet]
        for corners in local_faces:
            key = tuple(sorted(int(tags[i]) for i in corners))
            face_incidence.setdefault(key, []).append(tet_row)
    if any(len(parents) > 2 for parents in face_incidence.values()):
        _fail("T4 incidence contains a nonmanifold facet", field="mesh.npz.tets")

    exterior = _int_array(arrays, "exterior_triangles")
    interface = _int_array(arrays, "interface_triangles")
    if exterior.ndim != 2 or exterior.shape[1] != 3 or not 1 <= len(exterior) <= MAX_TRIANGLES \
            or interface.ndim != 2 or interface.shape[1] != 3 or len(interface) > MAX_TRIANGLES \
            or np.any(exterior < 0) or np.any(exterior >= len(nodes)) \
            or np.any(interface < 0) or np.any(interface >= len(nodes)):
        _fail("exterior and interface triangles must be bounded valid Tri3 row arrays")
    ext_tags = node_tags[exterior]
    int_tags = node_tags[interface]
    ext_keys = [tuple(sorted(map(int, row))) for row in ext_tags]
    int_keys = [tuple(sorted(map(int, row))) for row in int_tags]
    if len(set(ext_keys)) != len(ext_keys) or len(set(int_keys)) != len(int_keys) \
            or set(ext_keys) & set(int_keys):
        _fail("exterior and interface triangles must be unique and disjoint")
    if any(len(face_incidence.get(key, ())) != 1 for key in ext_keys):
        _fail("each exterior triangle must have exactly one T4 parent")
    if any(len(face_incidence.get(key, ())) != 2 for key in int_keys):
        _fail("each partition interface triangle must have exactly two T4 parents")
    cross_partition_internal = {
        key for key, parents_for_face in face_incidence.items()
        if len(parents_for_face) == 2
        and partitions[parents_for_face[0]] != partitions[parents_for_face[1]]
    }
    if set(int_keys) != cross_partition_internal:
        _fail("interface triangles do not cover every and only cross-partition internal facet")
    ext_fids = np.asarray(arrays["exterior_source_face_ids"])
    ext_roles = np.asarray(arrays["exterior_roles"])
    exact = np.asarray(arrays["exterior_ancestry_exact"])
    if ext_fids.shape != (len(exterior),) or ext_fids.dtype.kind not in "US" \
            or ext_roles.shape != (len(exterior),) or ext_roles.dtype.kind not in "US" \
            or exact.shape != (len(exterior),) or exact.dtype.kind != "b" or not exact.all():
        _fail("every exterior triangle needs a unique verified FID/role sidecar row")
    ext_face_indices = np.empty(len(exterior), dtype=np.int32)
    for index, (fid_raw, role_raw) in enumerate(zip(ext_fids, ext_roles, strict=True)):
        fid, role = str(fid_raw), str(role_raw)
        table_index = fid_to_index.get(fid)
        if table_index is None or not fid or role != face_rows[table_index]["role"]:
            _fail("exterior FID is unknown or conflicts with source role",
                  field=f"mesh.npz.exterior_source_face_ids[{index}]")
        ext_face_indices[index] = table_index
    if set(map(str, ext_fids.tolist())) != set(fid_to_index):
        _fail("not every verified source face identity appears on the exterior mesh",
              field="mesh.npz.exterior_source_face_ids")

    parents = _int_array(arrays, "interface_parents")
    parent_partition = _int_array(arrays, "interface_parent_partition")
    if parents.shape != (len(interface), 2) or parent_partition.shape != (len(interface), 2) \
            or np.any(parents < 0) or np.any(parents >= n_tets):
        _fail("interface parent rows/partitions do not align with interface triangles")
    interface_parent_tags = np.empty_like(parents)
    for index, (key, pair) in enumerate(zip(int_keys, parents, strict=True)):
        expected_parents = face_incidence[key]
        if len(set(map(int, pair))) != 2 or set(map(int, pair)) != set(expected_parents):
            _fail("interface parent rows do not match exact T4 incidence",
                  field=f"mesh.npz.interface_parents[{index}]")
        pair_parts = partitions[pair]
        if pair_parts[0] == pair_parts[1] or set(map(int, pair_parts)) != set(
                map(int, parent_partition[index])):
            _fail("interface parents must belong to distinct declared partitions",
                  field=f"mesh.npz.interface_parent_partition[{index}]")
        interface_parent_tags[index] = cell_tags[pair]

    # Independently parse and join the MSH by its original node and cell tags.
    msh = read_gmsh22(_bounded_file(mesh_path, MAX_MSH_BYTES, "mesh.msh"))
    msh_nodes = {int(tag): row for row, tag in enumerate(msh["node_tags"])}
    npz_node_row = {int(tag): row for row, tag in enumerate(node_tags)}
    used_npz_rows = np.unique(tets.reshape(-1))
    used_node_tags = {int(node_tags[row]) for row in used_npz_rows}
    if set(msh_nodes) != used_node_tags:
        _fail("MSH node tags must be exactly the T4-referenced NPZ tags; unused NPZ rows remain explicit",
              field="MSH.Nodes")
    if len(msh_nodes) != len(msh["node_tags"]):
        _fail("MSH node tags are not unique", field="MSH.Nodes")
    for tag, msh_row in msh_nodes.items():
        npz_row = npz_node_row.get(tag)
        if npz_row is None:
            _fail("MSH node tag is absent from the source NPZ", field="MSH.Nodes")
        if not np.allclose(msh["nodes"][msh_row], nodes[npz_row],
                           rtol=COORDINATE_ROUNDTRIP_RTOL,
                           atol=COORDINATE_ROUNDTRIP_ATOL_M):
            _fail("MSH/NPZ coordinates disagree beyond serialization tolerance",
                  field=f"MSH.Nodes.{tag}")

    msh_tets = msh["tetrahedra"]
    npz_cell_row = {int(tag): row for row, tag in enumerate(cell_tags)}
    if len(msh_tets) != n_tets or len({int(row["element_id"]) for row in msh_tets}) != n_tets \
            or {int(row["element_id"]) for row in msh_tets} != set(npz_cell_row):
        _fail("MSH T4 tags do not cover the NPZ original cell tags exactly", field="MSH.Elements.T4")
    for row in msh_tets:
        cell_tag = int(row["element_id"])
        npz_row = npz_cell_row[cell_tag]
        expected_corner_tags = tuple(int(v) for v in node_tags[tets[npz_row]])
        if tuple(map(int, row["node_tags"])) != expected_corner_tags:
            _fail("ordered MSH T4 corner tags differ from frozen NPZ connectivity",
                  field=f"MSH.T4.{cell_tag}")

    msh_triangles: dict[tuple[int, int, int], str | None] = {}
    for row in msh["triangles"]:
        key = tuple(sorted(map(int, row["node_tags"])))
        if key in msh_triangles:
            _fail("MSH repeats an explicit triangle", field="MSH.Elements.Tri3")
        msh_triangles[key] = row["physical_name"]
    if set(msh_triangles) != set(ext_keys) | set(int_keys):
        _fail("MSH explicit Tri3 faces do not equal the sidecar exterior/interface sets",
              field="MSH.Elements.Tri3")
    if any(msh_triangles[key] != exterior_physical_name for key in ext_keys):
        _fail("exterior Tri3 physical name conflicts with its FID sidecar",
              field="MSH.Elements.source_exterior")
    if any(msh_triangles[key] != interface_physical_name for key in int_keys):
        _fail("two-parent interface Tri3 must remain separately named and FID-null",
              field="MSH.Elements.cut_interface")

    unused_rows = np.setdiff1d(np.arange(len(nodes), dtype=np.int64), used_npz_rows,
                               assume_unique=True)
    return {
        "node_tags": node_tags,
        "unused_node_rows": unused_rows,
        "cell_tags": cell_tags,
        "cell_corner_node_tags": node_tags[tets],
        "exterior_triangle_node_tags": ext_tags,
        "exterior_source_face_index": ext_face_indices,
        "interface_triangle_node_tags": int_tags,
        "interface_parent_cell_tags": interface_parent_tags,
        "npz_node_count": len(nodes),
        "msh_node_count": len(msh_nodes),
        "tetrahedron_count": n_tets,
        "exterior_triangle_count": len(exterior),
        "interface_triangle_count": len(interface),
        "unused_node_count": len(unused_rows),
    }


def _freeze_array(value: np.ndarray) -> np.ndarray:
    result = np.array(value, copy=True)
    result.setflags(write=False)
    return result


def build_mother_binding(
    *,
    source_asset_dir: str | os.PathLike[str],
    partition_result_path: str | os.PathLike[str],
    audit_receipt_path: str | os.PathLike[str],
    trusted_audit_receipt_sha256: str,
    canonical_step_path: str | os.PathLike[str],
    canonical_step_input_key: str,
    native_face_identity_path: str | os.PathLike[str],
    witness_manifest_path: str | os.PathLike[str],
    v15_receipt_path: str | os.PathLike[str],
    label_import_receipt_path: str | os.PathLike[str],
    witness_directory: str | os.PathLike[str],
    exterior_physical_name: str = "source_exterior",
    interface_physical_name: str = "cut_interface",
) -> MotherBinding:
    """Build a read-only content binding for a source asset and partition output.

    ``source_asset_dir`` follows the VF reference candidate layout: its geometry
    subdirectory contains `solid.step`, `section.json`, `face_lineage.json`, and
    `surface.npz`; the root contains `geometry_report.json` and `provenance.json`.
    Other identity, witness, canonical STEP, result, audit, and witness directory
    paths are explicit. ``canonical_step_input_key`` names that file in the
    partition audit's input hash table; no asset basename is used as identity.

    The caller must supply the audit SHA-256 as a review anchor. This verifies
    that the passed receipt is the one the caller trusted; it is not a signature
    or a statement of physical review. The result intentionally records no
    fracture or observations and always sets training eligibility false.
    """
    if not _is_digest(trusted_audit_receipt_sha256):
        _fail("caller-supplied trusted audit SHA-256 is required",
              field="trusted_audit_receipt_sha256")
    root = Path(source_asset_dir).resolve()
    result_path = _bounded_file(Path(partition_result_path), MAX_JSON_BYTES, "partition_result")
    audit_path = _bounded_file(Path(audit_receipt_path), MAX_JSON_BYTES, "audit_receipt")
    audit_sha = _sha256(audit_path)
    if audit_sha != trusted_audit_receipt_sha256.lower():
        _fail("audit receipt does not match the caller-supplied trusted hash",
              field="trusted_audit_receipt_sha256")
    result = _read_json(result_path, "partition_result")
    audit = _read_json(audit_path, "audit_receipt")
    _reject_promoted_claims(result, "partition_result")
    _reject_promoted_claims(audit, "audit_receipt")
    msh_path, npz_path, file_hashes = _verify_receipt_files(result_path, result, audit_path, audit)

    result_qualification = result.get("qualification")
    if not isinstance(result_qualification, Mapping) \
            or result_qualification.get("fracture_or_fragment_result") != "NOT_CREATED" \
            or result_qualification.get("solver_qualification") != "NOT_ESTABLISHED":
        _fail("geometry binding requires explicit no-fracture/unqualified-solver statuses",
              field="partition_result.qualification")
    independent = audit.get("independent_mesh_audit")
    if not isinstance(independent, Mapping) or independent.get("status") != "PASS" \
            or independent.get("all_exterior_triangles_have_known_manifest_role") is not True \
            or independent.get("every_interface_triangle_has_two_distinct_partition_parents") is not True \
            or independent.get("interface_source_face_id_count") != 0:
        _fail("independent mesh audit does not verify exterior/interface separation",
              field="audit_receipt.independent_mesh_audit")
    nested_audit = independent.get("independent_audit")
    audit_checks = nested_audit.get("checks") if isinstance(nested_audit, Mapping) else None
    if not isinstance(audit_checks, Mapping) or not audit_checks \
            or any(value is not True for value in audit_checks.values()):
        _fail("independent tetra mesh checks are incomplete or failed",
              field="audit_receipt.independent_mesh_audit.independent_audit.checks")

    paths, source_values, report, identity, witness, _v15, _label, base_hashes = \
        _verify_source_hash_chain(root, Path(canonical_step_path), Path(native_face_identity_path),
                                  Path(witness_manifest_path), Path(v15_receipt_path),
                                  Path(label_import_receipt_path), canonical_step_input_key, audit)
    result_ancestry = result.get("source_face_ancestry")
    if not isinstance(result_ancestry, Mapping) \
            or result_ancestry.get("status") != "VERIFIED_OCC_WITNESS_ANCESTRY" \
            or result_ancestry.get("manifest_sha256") != _sha256(Path(witness_manifest_path)) \
            or result_ancestry.get("canonical_step_sha256") != _sha256(Path(canonical_step_path)) \
            or result_ancestry.get("upstream_source_step_sha256") != base_hashes["solid.step"] \
            or result_ancestry.get("source_face_lineage_sha256") != base_hashes["face_lineage.json"]:
        _fail("partition result does not bind the verified source-face ancestry chain",
              field="result.source_face_ancestry")
    provenance_path = root / "provenance.json"
    provenance = _read_json(provenance_path, "provenance.json")
    _reject_promoted_claims(provenance, "source_provenance")
    if report.get("candidate_geometry_status") != "PASS_STRUCTURAL_CAD_ONLY":
        _fail("source asset has not passed its structural geometry gate",
              field="geometry_report.candidate_geometry_status")
    report_checks = report.get("checks")
    if not isinstance(report_checks, Mapping) or not report_checks \
            or any(value is not True for value in report_checks.values()):
        _fail("source structural geometry report contains a missing or failed check",
              field="geometry_report.checks")
    if not isinstance(provenance.get("candidate_status"), str) \
            or not isinstance(provenance.get("source_profile_fidelity"), str):
        _fail("source provenance must preserve candidate and profile-fidelity statuses",
              field="provenance")

    lineage = source_values["lineage"]
    face_rows, fid_to_index = _face_table(lineage, identity, witness, Path(witness_directory))
    arrays = _load_npz(npz_path)
    mesh_binding = _validate_mesh_arrays(arrays, msh_path, face_rows, fid_to_index,
                                         exterior_physical_name, interface_physical_name)
    if int(result.get("counts", {}).get("tetrahedra", -1)) != mesh_binding["tetrahedron_count"] \
            or int(result.get("counts", {}).get("exterior_triangles", -1)) != mesh_binding["exterior_triangle_count"] \
            or int(result.get("counts", {}).get("interface_triangles", -1)) != mesh_binding["interface_triangle_count"]:
        _fail("partition result counts disagree with validated mesh arrays", field="result.counts")

    source_hashes = {
        "source_step": base_hashes["solid.step"],
        "canonical_step": _sha256(Path(canonical_step_path)),
        "section": base_hashes["section.json"],
        "face_lineage": base_hashes["face_lineage.json"],
        "source_surface": _sha256(paths["surface.npz"]),
        "source_provenance": _sha256(provenance_path),
        "geometry_report": base_hashes["geometry_report.json"],
        "native_face_identity": _sha256(Path(native_face_identity_path)),
        "witness_manifest": _sha256(Path(witness_manifest_path)),
        "v15_receipt": _sha256(Path(v15_receipt_path)),
        "label_import_receipt": _sha256(Path(label_import_receipt_path)),
        "partition_result": _sha256(result_path),
        "partition_msh": file_hashes["mesh.msh"]["sha256"],
        "partition_npz": file_hashes["mesh.npz"]["sha256"],
        "trusted_audit_receipt": audit_sha,
    }
    table_for_id = [{key: row[key] for key in (
        "source_face_id", "source_face_index", "native_face_index", "role",
        "source_segment_indices", "witness_sha256", "witness_bytes")}
        for row in face_rows]
    binding_id = hashlib.sha256(_canonical_bytes({
        "schema": SCHEMA,
        "source_hashes": source_hashes,
        "source_face_table": table_for_id,
    })).hexdigest()
    source_record = provenance.get("record")
    mother_id = "UNSPECIFIED_SOURCE_ASSET"
    if isinstance(source_record, Mapping):
        identifier = source_record.get("object_number")
        if isinstance(identifier, str) and identifier.strip():
            mother_id = identifier.strip()
    qualification = {
        "geometry_ancestry_status": "VERIFIED_OCC_WITNESS_ANCESTRY_WITHIN_DECLARED_OCC_TOLERANCES",
        "structural_geometry_status": report.get("candidate_geometry_status"),
        "candidate_status": provenance.get("candidate_status"),
        "source_profile_fidelity": provenance.get("source_profile_fidelity"),
        "physical_material_status": "UNKNOWN",
        "physical_material_source_statement": provenance.get("physical_material_properties", "UNKNOWN"),
        "solver_status": result["qualification"]["solver_qualification"],
        "fracture_status": result["qualification"]["fracture_or_fragment_result"],
        "observation_status": "NOT_CREATED_BY_GEOMETRY_BINDING",
        "physical_pass": False,
        "pretrain_eligible": False,
        "deployment_eligible": False,
    }
    summary: dict[str, Any] = {
        "schema": SCHEMA,
        "binding_id": binding_id,
        "mother_id": mother_id,
        "source_hashes": source_hashes,
        "source_face_table": table_for_id,
        "source_artifact_status": {
            "verified_geometry_files": ["solid.step", "section.json", "face_lineage.json", "surface.npz"],
            "missing_nonbinding_artifacts": source_values["optional_missing_artifacts"],
        },
        "counts": {key: int(value) for key, value in mesh_binding.items()
                   if key.endswith("_count")},
        "unused_node_rows_sha256": hashlib.sha256(
            np.asarray(mesh_binding["unused_node_rows"], dtype="<i8").tobytes()).hexdigest(),
        "qualification": qualification,
        "lineage_scope": {
            "mesh_identity": "original Gmsh node/cell tags and exact ordered T4 corner tags",
            "post_insertion_bulk_row_local_corner_identity": (
                "CONDITIONAL_ON_NATIVE_BACKEND_PRESERVING_FROZEN_INITIAL_ROWS"),
            "coordinate_check": (
                "MSH/NPZ serialization consistency only; tag/connectivity identity is primary; "
                "no nearest-coordinate inference"),
            "coordinate_roundtrip_rtol": COORDINATE_ROUNDTRIP_RTOL,
            "coordinate_roundtrip_atol_m": COORDINATE_ROUNDTRIP_ATOL_M,
            "cohesive_qp_ordinal_stability": "NOT_CLAIMED",
            "interface_geometry": "TWO_PARENT_PARTITION_INTERFACE; SOURCE_EXTERIOR_FID_NULL",
            "source_face_identity": (
                "CONDITIONAL_CAD_ANCESTRY_WITHIN_DECLARED_OCC_TOLERANCES; "
                "NOT_LITERAL_TOPODS_IDENTITY_OR_EXACT_METROLOGY"),
        },
        "arrays_file": "binding_arrays.npz",
        "array_semantics": {
            "node_tags": "all original NPZ node tags, including rows unused by serialized T4s",
            "unused_node_rows": "NPZ storage row indices not referenced by any T4; not node identity",
            "cell_corner_node_tags": "exact ordered original node tags for each original cell tag",
            "exterior_source_face_index": "index into source_face_table; never a numeric role code",
            "interface_parent_cell_tags": "two distinct original cell tags for each FID-null interface triangle",
        },
    }
    return MotherBinding(
        binding_id=binding_id,
        mother_id=mother_id,
        summary=_deep_freeze(summary),
        source_face_ids=tuple(row["source_face_id"] for row in face_rows),
        source_face_roles=tuple(row["role"] for row in face_rows),
        node_tags=_freeze_array(mesh_binding["node_tags"]),
        unused_node_rows=_freeze_array(mesh_binding["unused_node_rows"]),
        cell_tags=_freeze_array(mesh_binding["cell_tags"]),
        cell_corner_node_tags=_freeze_array(mesh_binding["cell_corner_node_tags"]),
        exterior_triangle_node_tags=_freeze_array(mesh_binding["exterior_triangle_node_tags"]),
        exterior_source_face_index=_freeze_array(mesh_binding["exterior_source_face_index"]),
        interface_triangle_node_tags=_freeze_array(mesh_binding["interface_triangle_node_tags"]),
        interface_parent_cell_tags=_freeze_array(mesh_binding["interface_parent_cell_tags"]),
    )


def write_mother_binding(binding: MotherBinding,
                         output_dir: str | os.PathLike[str], *,
                         max_output_bytes: int = MAX_OUTPUT_BYTES) -> dict[str, Any]:
    """Write bounded mapping arrays and JSON receipt without copying source arrays."""
    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    arrays_path = root / "binding_arrays.npz"
    receipt_path = root / "binding_receipt.json"
    if arrays_path.exists() or receipt_path.exists():
        _fail("binding output files already exist", field="output_dir")
    try:
        with arrays_path.open("xb") as stream:
            np.savez_compressed(
                stream,
                node_tags=binding.node_tags,
                unused_node_rows=binding.unused_node_rows,
                cell_tags=binding.cell_tags,
                cell_corner_node_tags=binding.cell_corner_node_tags,
                exterior_triangle_node_tags=binding.exterior_triangle_node_tags,
                exterior_source_face_index=binding.exterior_source_face_index,
                interface_triangle_node_tags=binding.interface_triangle_node_tags,
                interface_parent_cell_tags=binding.interface_parent_cell_tags,
            )
        array_size = arrays_path.stat().st_size
        if array_size > MAX_NPZ_BYTES:
            _fail("binding arrays exceed NPZ output cap", field="binding_arrays.npz")
        receipt = _deep_thaw(binding.summary)
        receipt["files"] = {"binding_arrays.npz": {
            "bytes": array_size, "sha256": _sha256(arrays_path)}}
        raw = json.dumps(receipt, indent=2, sort_keys=True, ensure_ascii=False,
                         allow_nan=False).encode("utf-8") + b"\n"
        if len(raw) > MAX_JSON_BYTES or array_size + len(raw) > max_output_bytes:
            _fail("binding output exceeds configured total output cap", field="output_dir")
        with receipt_path.open("xb") as stream:
            stream.write(raw)
    except Exception:
        for path in (arrays_path, receipt_path):
            try:
                path.unlink(missing_ok=True)
            except OSError:
                pass
        raise
    return receipt


__all__ = [
    "SCHEMA", "MAX_JSON_BYTES", "MAX_NPZ_BYTES", "MAX_MSH_BYTES",
    "MAX_TETRAHEDRA", "MAX_OUTPUT_BYTES", "MotherBindingError",
    "MotherBinding", "build_mother_binding", "write_mother_binding",
]
