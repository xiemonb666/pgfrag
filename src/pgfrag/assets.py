from __future__ import annotations

import json
import math
from pathlib import Path, PurePosixPath, PureWindowsPath
from typing import Any
from zipfile import BadZipFile, ZipFile

import numpy as np

from .errors import ValidationError
from .hashing import sha256_file


def _package_file(root: Path, relative: Any) -> Path:
    if not isinstance(relative, str) or not relative or "\\" in relative or ":" in relative or "\x00" in relative:
        raise ValidationError(f"invalid package-relative path: {relative!r}")
    pure = PurePosixPath(relative)
    if pure.is_absolute() or any(part in {"", ".", ".."} for part in pure.parts):
        raise ValidationError(f"unsafe package-relative path: {relative!r}")
    root_real = root.resolve(strict=True)
    cursor = root
    for part in pure.parts:
        cursor = cursor / part
        if cursor.is_symlink():
            raise ValidationError(f"symlink is not permitted in asset packages: {relative!r}")
    try:
        resolved = cursor.resolve(strict=True)
        resolved.relative_to(root_real)
    except (OSError, ValueError) as exc:
        raise ValidationError(f"asset path is missing or escapes package root: {relative!r}") from exc
    if not resolved.is_file():
        raise ValidationError(f"asset entry is not a regular file: {relative!r}")
    return resolved


def _check_npz(path: Path, *, max_declared_uncompressed_bytes: int) -> int:
    try:
        with ZipFile(path) as archive:
            total = 0
            for member in archive.infolist():
                name = PurePosixPath(member.filename)
                if name.is_absolute() or any(part in {"", ".", ".."} for part in name.parts):
                    raise ValidationError(f"unsafe ZIP member in NPZ file {path.name!r}")
                if member.file_size < 0:
                    raise ValidationError(f"invalid declared member size in NPZ file {path.name!r}")
                total += member.file_size
            if total > max_declared_uncompressed_bytes:
                raise ValidationError(f"NPZ declares {total} uncompressed bytes, above the configured cap")
            return total
    except (OSError, BadZipFile) as exc:
        raise ValidationError(f"invalid NPZ archive {path.name!r}: {exc}") from exc


def _load_npz(path: Path, *, max_declared_uncompressed_bytes: int) -> tuple[dict[str, np.ndarray], int]:
    declared_bytes = _check_npz(path, max_declared_uncompressed_bytes=max_declared_uncompressed_bytes)
    try:
        with np.load(path, allow_pickle=False) as archive:
            arrays = {key: np.asarray(archive[key]) for key in archive.files}
    except (KeyError, ValueError, OSError) as exc:
        raise ValidationError(f"surface or fields NPZ is incomplete/invalid: {exc}") from exc
    return arrays, declared_bytes


def _require_array(arrays: dict[str, np.ndarray], key: str) -> np.ndarray:
    if key not in arrays:
        raise ValidationError(f"missing array {key!r}")
    return arrays[key]


def _floating_n3(array: np.ndarray, name: str, rows: int | None = None) -> None:
    if array.dtype.kind != "f" or array.ndim != 2 or array.shape[1] != 3 or (rows is not None and array.shape != (rows, 3)):
        raise ValidationError(f"{name} must be a floating-point N×3 array")
    if len(array) < 4 or not np.isfinite(array).all():
        raise ValidationError(f"{name} must have at least four finite rows")


def validate_asset(package_dir: str | Path, *, max_npz_declared_uncompressed_bytes: int = 512 * 1024 * 1024) -> dict[str, Any]:
    """Fail-closed forge.asset.v1 package-contract validator.

    This mirrors the current VesselForge contract locally without importing
    VesselForge. Geometry report status is checked against its hashed report,
    but PG-Frag does not rerun the producer's geometric audit.
    """
    root = Path(package_dir).resolve(strict=True)
    if not root.is_dir():
        raise ValidationError(f"asset package directory does not exist: {root}")
    manifest_path = _package_file(root, "manifest.json")
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValidationError(f"manifest.json is not valid UTF-8 JSON: {exc}") from exc
    if not isinstance(manifest, dict) or manifest.get("schema") != "forge.asset.v1":
        raise ValidationError("schema must be forge.asset.v1")
    if manifest.get("units") != {"length": "mm", "time": "s", "mass": "kg"}:
        raise ValidationError("units do not match forge.asset.v1 canonical units")
    if manifest.get("frame") != "right_handed_z_up":
        raise ValidationError("frame must be right_handed_z_up")
    for key in ("geometry", "fields", "material", "ornament", "files", "provenance", "validation"):
        if not isinstance(manifest.get(key), dict):
            raise ValidationError(f"manifest.{key} must be an object")

    files = manifest["files"]
    if "manifest.json" in files or not files:
        raise ValidationError("files must enumerate package payloads and exclude manifest.json")
    verified_files: list[dict[str, Any]] = []
    for relative, expected in sorted(files.items()):
        if not isinstance(expected, dict) or not {"sha256", "bytes"}.issubset(expected):
            raise ValidationError(f"files[{relative!r}] must provide sha256 and bytes")
        digest = expected["sha256"]
        if not isinstance(digest, str) or len(digest) != 64 or any(c not in "0123456789abcdefABCDEF" for c in digest):
            raise ValidationError(f"files[{relative!r}].sha256 must be hexadecimal SHA-256")
        path = _package_file(root, relative)
        if isinstance(expected["bytes"], bool) or path.stat().st_size != expected["bytes"]:
            raise ValidationError(f"byte count mismatch for {relative!r}")
        actual = sha256_file(path)
        if actual.lower() != digest.lower():
            raise ValidationError(f"sha256 mismatch for {relative!r}")
        verified_files.append({"path": relative, "bytes": path.stat().st_size, "sha256": actual})

    geometry, fields = manifest["geometry"], manifest["fields"]
    surface_name = geometry.get("surface_npz")
    fields_name = fields.get("surface_fields_npz")
    report_name = manifest["validation"].get("report")
    if geometry.get("solid_step") is not None and geometry["solid_step"] not in files:
        raise ValidationError("geometry.solid_step must reference a hash-listed file")
    for relative in (surface_name, fields_name, report_name):
        if not isinstance(relative, str) or relative == "manifest.json" or relative not in files:
            raise ValidationError(f"manifest references a payload not enumerated in files: {relative!r}")
    if manifest["validation"].get("scope") != "geometry_only":
        raise ValidationError("prototype validation scope must be geometry_only")
    declared_status = manifest["validation"].get("status")
    if declared_status not in {"PASS", "FAIL", "PENDING"}:
        raise ValidationError("validation status is invalid")
    report_path = _package_file(root, report_name)
    try:
        report = json.loads(report_path.read_text(encoding="utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValidationError(f"geometry validation report is not valid UTF-8 JSON: {exc}") from exc
    if not isinstance(report, dict) or report.get("scope") != "geometry_only":
        raise ValidationError("geometry validation report must declare geometry_only scope")
    if report.get("status") != declared_status:
        raise ValidationError("manifest validation status disagrees with the referenced report")

    surface_path = _package_file(root, surface_name)
    fields_path = _package_file(root, fields_name)
    surface, surface_payload_bytes = _load_npz(surface_path, max_declared_uncompressed_bytes=max_npz_declared_uncompressed_bytes)
    field_data, fields_payload_bytes = _load_npz(fields_path, max_declared_uncompressed_bytes=max_npz_declared_uncompressed_bytes)
    vertices = _require_array(surface, "vertices_mm")
    faces = _require_array(surface, "faces")
    _floating_n3(vertices, "vertices_mm")
    if faces.dtype.kind not in "iu" or faces.ndim != 2 or faces.shape[1] != 3 or len(faces) < 4:
        raise ValidationError("faces must be an integer F×3 array")
    if int(faces.min()) < 0 or int(faces.max()) >= len(vertices):
        raise ValidationError("faces contains an out-of-range vertex index")
    if np.any((faces[:, 0] == faces[:, 1]) | (faces[:, 1] == faces[:, 2]) | (faces[:, 2] == faces[:, 0])):
        raise ValidationError("faces may not repeat a vertex index")

    n_vertices, n_faces = len(vertices), len(faces)
    region = _require_array(field_data, "face_region")
    albedo = _require_array(field_data, "albedo_linear_rgb")
    roughness = _require_array(field_data, "roughness")
    material_xyz = _require_array(field_data, "vertices_material_xyz_mm")
    if region.shape != (n_faces,) or region.dtype.kind not in "iu":
        raise ValidationError("face_region must be an integer array with one value per face")
    if albedo.shape != (n_faces, 3) or albedo.dtype.kind != "f" or not np.isfinite(albedo).all():
        raise ValidationError("albedo_linear_rgb must be a finite floating F×3 array")
    if np.any(albedo < 0) or np.any(albedo > 1):
        raise ValidationError("albedo_linear_rgb values must lie in [0,1]")
    if roughness.shape != (n_faces,) or roughness.dtype.kind != "f" or not np.isfinite(roughness).all():
        raise ValidationError("roughness must be a finite floating F array")
    if np.any(roughness < 0) or np.any(roughness > 1):
        raise ValidationError("roughness values must lie in [0,1]")
    _floating_n3(material_xyz, "vertices_material_xyz_mm", rows=n_vertices)
    if not np.array_equal(material_xyz, vertices):
        raise ValidationError("initial asset material coordinates must exactly equal canonical vertices")

    material = manifest["material"]
    material_id = material.get("id")
    if not isinstance(material_id, str) or not material_id.strip():
        raise ValidationError("material.id is required")
    mechanics = material.get("mechanics")
    if not isinstance(mechanics, dict) or not {"rho_kg_m3", "E_pa", "nu", "Gc_j_m2"}.issubset(mechanics):
        raise ValidationError("material.mechanics must use canonical SI-named values")
    for name in ("rho_kg_m3", "E_pa", "Gc_j_m2"):
        value = mechanics.get(name)
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
            raise ValidationError(f"material.mechanics.{name} must be finite and positive")
    poisson = mechanics.get("nu")
    if isinstance(poisson, bool) or not isinstance(poisson, (int, float)) or not math.isfinite(poisson) or not (-1.0 < poisson < 0.5):
        raise ValidationError("material.mechanics.nu must be finite and in (-1, 0.5)")
    if material.get("calibration_status") != "UNVALIDATED_PROTOTYPE":
        raise ValidationError("prototype mechanics must remain UNVALIDATED_PROTOTYPE")

    ornament = manifest["ornament"]
    recipe = ornament.get("recipe")
    if ornament.get("kind") != "painted" or not isinstance(recipe, dict):
        raise ValidationError("prototype ornament must be a painted rule recipe")
    if recipe.get("material_binding_id") != material_id:
        raise ValidationError("ornament material binding does not match material.id")
    if recipe.get("geometry_relief_mm") != 0.0:
        raise ValidationError("painted-only prototype must declare zero geometric relief")
    if recipe.get("historical_style_validation") != "NOT_VALIDATED":
        raise ValidationError("historical style status must remain NOT_VALIDATED")
    allowed_regions = {1, 2, 3, 4, 5, 255}
    if not set(np.unique(region).astype(int)).issubset(allowed_regions):
        raise ValidationError("face_region contains an undeclared semantic code")
    if set(recipe.get("target_region_codes", [])) - allowed_regions:
        raise ValidationError("ornament targets contain an unknown semantic region")

    return {
        "schema": "forge.asset.v1",
        "contract_status": "PASS",
        "manifest_sha256": sha256_file(manifest_path),
        "manifest_path": str(manifest_path),
        "surface": {"vertices": n_vertices, "faces": n_faces, "path": surface_name},
        "fields_path": fields_name,
        "verified_files": verified_files,
        "material_id": material_id,
        "ornament_material_binding_id": recipe["material_binding_id"],
        "calibration_status": material["calibration_status"],
        "declared_geometry_validation": {
            "status": declared_status, "scope": "geometry_only", "report": report_name,
            "report_status_matches_manifest": True,
            "meaning": "report status was checked and payload hash verified; producer geometry audit was not rerun",
        },
        "vertices_material_equal_rest_vertices": True,
        "npz_declared_uncompressed_bytes": {surface_name: surface_payload_bytes, fields_name: fields_payload_bytes},
        "npz_payload_limit_note": "ZIP-declared uncompressed sizes were checked before NumPy loading; this is not a hard process memory limit.",
        "mechanical_solve_supported_from_surface_asset": False,
        "unsupported_reason": "surface-to-volume point discretization and void-topology classification are not implemented in this prototype",
        "reserved_region_semantics": "region code 255 is carried as an asset semantic code; PG-Frag does not infer or label it as dynamic fracture",
    }
