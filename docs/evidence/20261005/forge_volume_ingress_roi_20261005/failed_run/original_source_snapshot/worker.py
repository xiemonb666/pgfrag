"""Isolated, bounded VF STEP-to-material-point ingress candidate.

SOURCE_ONLY: this file has not been run. It is intentionally outside pgfrag's
production package and is not called by backend.py or any solver path.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import multiprocessing as mp
import os
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Iterable


EXPERIMENT_DIR = Path(__file__).resolve().parent
REPO_ROOT = EXPERIMENT_DIR.parents[2]
TMP_ROOT = (EXPERIMENT_DIR / "tmp").resolve()
FROZEN_ASSET = (
    REPO_ROOT
    / "vesselforge/tmp/review_20261005/geometry_resume/migration_copy/cylinder_asset_run01"
).resolve()
EXPECTED_MANIFEST_SHA256 = "df781dc6c818512af8d66b164288aa0d5991097ea0166b85556119de2c9fba41"
EXPECTED_STEP_SHA256 = "5c68bbc29b2b0bc684cc6c551c2769872c4c490d588b9a3ea2cbb98e7f74135e"
EXPECTED_VF_PYTHON = (
    REPO_ROOT / "third_party/vesselforge-runtime/venv/Scripts/python.exe"
).resolve()
DIAGNOSTIC_ROI_BOUNDS_MM = (7.5, 10.5, -1.0, 1.0, 4.0, 6.0)
DIAGNOSTIC_CELL_SIZE_MM = 0.5

MAX_CANDIDATE_CELLS = 256
MAX_MATERIAL_POINTS = 500
MAX_CONNECTED_SOLIDS_PER_CELL = 16
MAX_REFERENCE_ROI_COMPONENTS = 16
MAX_COMMON_FACES_PER_CELL = 512
MAX_OUTPUT_BYTES = 20 * 1024 * 1024
PARENT_RECEIPT_RESERVE_BYTES = 64 * 1024
WORKER_TIMEOUT_SECONDS = 180
INSIDE_TOLERANCE_MM = 1e-7
VOLUME_RELATIVE_TOLERANCE = 1e-8
VOLUME_ABSOLUTE_TOLERANCE_MM3 = 1e-9


class IngressBlocked(RuntimeError):
    """A fail-closed input, geometry, budget, or conservation condition."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while True:
            chunk = stream.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def _json_bytes(value: Any) -> bytes:
    return (
        json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False)
        + "\n"
    ).encode("utf-8")


def _package_payload_path(asset_dir: Path, relative: str) -> Path:
    if not isinstance(relative, str) or not relative:
        raise IngressBlocked("manifest contains an empty payload path")
    normalized = relative.replace("\\", "/")
    rel_path = Path(normalized)
    if rel_path.is_absolute() or any(part in {"", ".", ".."} for part in normalized.split("/")):
        raise IngressBlocked(f"unsafe package-relative payload path: {relative!r}")
    candidate = asset_dir / Path(*normalized.split("/"))
    if candidate.is_symlink():
        raise IngressBlocked(f"symlinks are not accepted in the frozen package: {relative!r}")
    path = candidate.resolve(strict=True)
    if not path.is_relative_to(asset_dir):
        raise IngressBlocked(f"payload resolves outside the asset package: {relative!r}")
    if not path.is_file():
        raise IngressBlocked(f"payload is not a regular package file: {relative!r}")
    return path


def _write_stage_json(
    stage_dir: Path,
    name: str,
    value: Any,
    *,
    reserve_parent_receipt: bool = True,
) -> None:
    path = stage_dir / name
    payload = _json_bytes(value)
    existing_bytes = sum(
        item.stat().st_size
        for item in stage_dir.iterdir()
        if item.is_file() and item != path and not item.name.endswith(".tmp")
    )
    budget = MAX_OUTPUT_BYTES - (
        PARENT_RECEIPT_RESERVE_BYTES if reserve_parent_receipt else 0
    )
    if existing_bytes + len(payload) > budget:
        raise IngressBlocked(
            f"bounded stage would exceed {budget} bytes while writing {name}"
        )
    temporary = stage_dir / f".{name}.tmp"
    temporary.write_bytes(payload)
    os.replace(temporary, path)


def _stage_inventory(stage_dir: Path) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    total = 0
    for path in sorted(stage_dir.iterdir(), key=lambda item: item.name):
        if not path.is_file():
            continue
        size = path.stat().st_size
        total += size
        rows.append({"path": path.name, "bytes": size, "sha256": _sha256(path)})
    return {"total_bytes": total, "files": rows}


def _read_json_record(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
        return value if isinstance(value, dict) else {"read_error": "JSON root is not an object"}
    except Exception as exc:
        return {"read_error": f"{type(exc).__name__}: {exc}"}


def _load_frozen_package() -> tuple[dict[str, Any], dict[str, str], dict[str, Any]]:
    asset_dir = FROZEN_ASSET.resolve(strict=True)
    if not asset_dir.is_dir():
        raise IngressBlocked("frozen VF asset path is not a directory")
    manifest_path = (asset_dir / "manifest.json").resolve(strict=True)
    manifest_sha = _sha256(manifest_path)
    if manifest_sha != EXPECTED_MANIFEST_SHA256:
        raise IngressBlocked(
            f"frozen manifest SHA mismatch: expected {EXPECTED_MANIFEST_SHA256}, got {manifest_sha}"
        )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != "forge.asset.v1":
        raise IngressBlocked("input is not forge.asset.v1")
    if manifest.get("units") != {"length": "mm", "time": "s", "mass": "kg"}:
        raise IngressBlocked("forge.asset.v1 units do not match the ingress contract")
    if manifest.get("frame") != "right_handed_z_up":
        raise IngressBlocked("unsupported coordinate frame")
    if manifest.get("validation", {}).get("scope") != "geometry_only":
        raise IngressBlocked("input validation scope must remain geometry_only")

    files = manifest.get("files")
    if not isinstance(files, dict) or not files:
        raise IngressBlocked("manifest.files must list the frozen payloads")
    input_hashes: dict[str, str] = {"manifest.json": manifest_sha}
    for relative, record in sorted(files.items()):
        if not isinstance(record, dict):
            raise IngressBlocked(f"invalid manifest file record: {relative!r}")
        path = _package_payload_path(asset_dir, relative)
        actual = _sha256(path)
        expected = record.get("sha256")
        expected_bytes = record.get("bytes")
        if actual.lower() != str(expected).lower():
            raise IngressBlocked(f"payload SHA mismatch for {relative!r}")
        if isinstance(expected_bytes, bool) or path.stat().st_size != expected_bytes:
            raise IngressBlocked(f"payload byte length mismatch for {relative!r}")
        input_hashes[relative] = actual

    step_relative = manifest.get("geometry", {}).get("solid_step")
    if not isinstance(step_relative, str) or step_relative not in files:
        raise IngressBlocked("manifest.geometry.solid_step must be a hash-listed payload")
    if input_hashes[step_relative] != EXPECTED_STEP_SHA256:
        raise IngressBlocked(
            f"frozen STEP SHA mismatch: expected {EXPECTED_STEP_SHA256}, got {input_hashes[step_relative]}"
        )

    semantics_path = _package_payload_path(asset_dir, "metadata/face_semantics.json")
    semantics = json.loads(semantics_path.read_text(encoding="utf-8"))
    if semantics.get("status") != "PASS":
        raise IngressBlocked("canonical mother face semantics are not PASS")
    semantic_faces = semantics.get("faces")
    if not isinstance(semantic_faces, list) or not semantic_faces:
        raise IngressBlocked("canonical mother face semantic inventory is missing")
    for face in semantic_faces:
        if face.get("status") != "CLASSIFIED":
            raise IngressBlocked("mother face semantic inventory contains an unresolved face")
        if face.get("canonical_step_sha256") != EXPECTED_STEP_SHA256:
            raise IngressBlocked("face semantics do not identify the frozen canonical STEP")

    wall_reference_path = _package_payload_path(asset_dir, "metadata/wall_reference.json")
    wall_reference = json.loads(wall_reference_path.read_text(encoding="utf-8"))
    if wall_reference.get("schema") != "vf.wall-reference.v1":
        raise IngressBlocked("wall reference metadata schema is unknown")
    if wall_reference.get("qualification", {}).get("global_minimum_wall_status") != "NOT_ESTABLISHED":
        raise IngressBlocked("unexpected wall qualification status; stop for review")
    for station in wall_reference.get("stations", []):
        source_evidence = station.get("source_classification_evidence", {})
        if source_evidence.get("canonical_step_sha256") != EXPECTED_STEP_SHA256:
            raise IngressBlocked("wall station does not identify the frozen canonical STEP")

    mechanics = manifest.get("material", {}).get("mechanics", {})
    rho = mechanics.get("rho_kg_m3")
    if isinstance(rho, bool) or not isinstance(rho, (int, float)) or not math.isfinite(rho) or rho <= 0:
        raise IngressBlocked("material rho_kg_m3 must be finite and positive")
    if manifest.get("material", {}).get("calibration_status") != "UNVALIDATED_PROTOTYPE":
        raise IngressBlocked("mechanics calibration status must remain UNVALIDATED_PROTOTYPE")
    declared_wall = manifest.get("provenance", {}).get("spec", {}).get("parameters", {}).get("wall_mm")
    if isinstance(declared_wall, bool) or not isinstance(declared_wall, (int, float)) or declared_wall <= 0:
        declared_wall = None

    context = {
        "asset_dir": asset_dir,
        "step_path": _package_payload_path(asset_dir, step_relative),
        "step_relative": step_relative,
        "rho_kg_m3": float(rho),
        "declared_wall_parameter_mm": float(declared_wall) if declared_wall is not None else None,
        "semantics": semantics,
        "wall_reference": wall_reference,
    }
    return manifest, input_hashes, context


def _bbox_values(shape: Any, brep_bounds_mm: Any) -> tuple[float, float, float, float, float, float]:
    result = tuple(float(value) for value in brep_bounds_mm(shape))
    if not all(math.isfinite(value) for value in result):
        raise IngressBlocked("B-rep returned non-finite bounds")
    if not (result[0] < result[1] and result[2] < result[3] and result[4] < result[5]):
        raise IngressBlocked("B-rep bounds are degenerate")
    return result


def _make_box(cq: Any, bounds: tuple[float, float, float, float, float, float]) -> Any:
    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    return cq.Solid.makeBox(
        xmax - xmin,
        ymax - ymin,
        zmax - zmin,
        cq.Vector(xmin, ymin, zmin),
    )


def _grid_axis(lo: float, hi: float, step: float) -> list[tuple[float, float]]:
    span = hi - lo
    quotient = span / step
    count = int(round(quotient))
    if count < 1 or abs(quotient - count) > 1e-12 * max(1.0, abs(quotient)):
        raise IngressBlocked(
            f"ROI axis span {span:g} mm is not exactly divisible by h={step:g} mm"
        )
    cells: list[tuple[float, float]] = []
    for index in range(count):
        start = lo + index * step
        end = lo + (index + 1) * step
        if not start < end:
            raise IngressBlocked("grid axis produced a non-positive cell")
        cells.append((start, end))
    if abs(cells[0][0] - lo) > 1e-12 or abs(cells[-1][1] - hi) > 1e-12 * max(1.0, abs(hi)):
        raise IngressBlocked("grid did not exactly cover the declared ROI")
    return cells


def _volume_tolerance(reference_mm3: float) -> float:
    return max(
        VOLUME_ABSOLUTE_TOLERANCE_MM3,
        VOLUME_RELATIVE_TOLERANCE * abs(reference_mm3),
    )


def _component_sort_key(solid: Any, brep_bounds_mm: Any) -> tuple[float, ...]:
    bounds = _bbox_values(solid, brep_bounds_mm)
    center = tuple(float(value) for value in solid.Center().toTuple())
    return bounds + center


def _classify_point_in_source(source_solid: Any, point_mm: tuple[float, float, float]) -> str:
    from OCP.BRepClass3d import BRepClass3d_SolidClassifier
    from OCP.TopAbs import TopAbs_IN, TopAbs_ON, TopAbs_OUT
    from OCP.gp import gp_Pnt

    try:
        classifier = BRepClass3d_SolidClassifier()
        classifier.Load(source_solid.wrapped)
        classifier.Perform(gp_Pnt(*point_mm), INSIDE_TOLERANCE_MM)
        state = classifier.State()
    except Exception:
        return "UNKNOWN"
    if state == TopAbs_IN:
        return "IN"
    if state == TopAbs_ON:
        return "ON"
    if state == TopAbs_OUT:
        return "OUT"
    return "UNKNOWN"


def _solid_topology_report(solid: Any) -> dict[str, Any]:
    from OCP.BRepCheck import BRepCheck_Analyzer
    from OCP.TopAbs import TopAbs_SHELL
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS

    analyzer = BRepCheck_Analyzer(solid.wrapped)
    shell_explorer = TopExp_Explorer(solid.wrapped, TopAbs_SHELL)
    shell_count = 0
    open_shell_count = 0
    while shell_explorer.More():
        shell_count += 1
        shell = TopoDS.Shell_s(shell_explorer.Current())
        if not shell.Closed():
            open_shell_count += 1
        shell_explorer.Next()
    valid = bool(analyzer.IsValid())
    closed = shell_count == 1 and open_shell_count == 0
    return {
        "brep_check_valid": valid,
        "shell_count": shell_count,
        "open_shell_count": open_shell_count,
        "one_closed_shell": closed,
        "status": "PASS" if valid and closed else "BLOCKED",
    }


def _resolve_roi(
    requested_min: list[float] | None,
    requested_max: list[float] | None,
    source_bounds: tuple[float, float, float, float, float, float],
) -> tuple[tuple[float, float, float, float, float, float], str]:
    if (requested_min is None) != (requested_max is None):
        raise IngressBlocked("ROI minimum and maximum must be supplied together")
    if requested_min is None:
        return source_bounds, "FULL_SOURCE_BREP_BOUNDS"
    values = tuple(float(value) for value in requested_min + requested_max)
    if len(values) != 6 or not all(math.isfinite(value) for value in values):
        raise IngressBlocked("ROI bounds must contain six finite values")
    roi = (values[0], values[3], values[1], values[4], values[2], values[5])
    if not (roi[0] < roi[1] and roi[2] < roi[3] and roi[4] < roi[5]):
        raise IngressBlocked("ROI bounds are degenerate")
    return roi, "EXPLICIT_DIAGNOSTIC_ROI_PARTIAL_SOURCE_INTERSECTION"


def _check_time(started: float) -> None:
    if time.monotonic() - started > WORKER_TIMEOUT_SECONDS:
        raise IngressBlocked("worker elapsed-time budget exceeded")


def _compute_candidate(request: dict[str, Any], stage_dir: Path) -> dict[str, Any]:
    import cadquery as cq

    vf_source = str(REPO_ROOT / "vesselforge/src")
    if vf_source not in sys.path:
        sys.path.insert(0, vf_source)
    from vesselforge.brep_identity import brep_bounds_mm, brep_face_bounds_policy_record

    started = time.monotonic()
    manifest, input_hashes_before, context = _load_frozen_package()
    source_shape = cq.importers.importStep(str(context["step_path"])).val()
    _check_time(started)
    if not source_shape.isValid():
        raise IngressBlocked("canonical STEP import is not a valid shape")
    source_solids = list(source_shape.Solids())
    if len(source_solids) != 1:
        raise IngressBlocked("frozen input must contain exactly one connected material solid")
    source_solid = source_solids[0]
    source_topology = _solid_topology_report(source_solid)
    if source_topology["status"] != "PASS":
        raise IngressBlocked(f"source material solid is not one valid closed shell: {source_topology}")
    source_bounds = _bbox_values(source_solid, brep_bounds_mm)
    roi_min = request["roi_min_mm"]
    roi_max = request["roi_max_mm"]
    roi_bounds, roi_mode = _resolve_roi(roi_min, roi_max, source_bounds)
    cell_size = float(request["cell_size_mm"])
    if not math.isfinite(cell_size) or cell_size != DIAGNOSTIC_CELL_SIZE_MM:
        raise IngressBlocked("worker cell size is frozen at the root-selected 0.5 mm")
    if tuple(roi_bounds) != DIAGNOSTIC_ROI_BOUNDS_MM:
        raise IngressBlocked("worker ROI is frozen to the root-selected diagnostic ROI")

    x_cells = _grid_axis(roi_bounds[0], roi_bounds[1], cell_size)
    y_cells = _grid_axis(roi_bounds[2], roi_bounds[3], cell_size)
    z_cells = _grid_axis(roi_bounds[4], roi_bounds[5], cell_size)
    candidate_count = len(x_cells) * len(y_cells) * len(z_cells)
    if candidate_count != 96 or candidate_count > MAX_CANDIDATE_CELLS:
        raise IngressBlocked(f"frozen diagnostic grid must contain exactly 96 cells; got {candidate_count}")

    roi_common = source_solid.intersect(_make_box(cq, roi_bounds))
    _check_time(started)
    roi_shape_valid = bool(roi_common.isValid())
    roi_common_volume_mm3 = float(roi_common.Volume())
    roi_common_solids = list(roi_common.Solids())
    reference_component_records: list[dict[str, Any]] = []
    reference_blockers: list[str] = []
    if not roi_shape_valid:
        reference_blockers.append("ROI_COMMON_INVALID")
    if not math.isfinite(roi_common_volume_mm3) or roi_common_volume_mm3 < 0:
        reference_blockers.append("ROI_COMMON_VOLUME_INVALID")
    if len(roi_common_solids) > MAX_REFERENCE_ROI_COMPONENTS:
        reference_blockers.append("ROI_REFERENCE_COMPONENT_COUNT_EXCEEDS_CAP")
    reference_components_truncated = len(roi_common_solids) > MAX_REFERENCE_ROI_COMPONENTS
    for component_index, component in enumerate(roi_common_solids[:MAX_REFERENCE_ROI_COMPONENTS]):
        component_volume = float(component.Volume())
        reference_component_records.append(
            {
                "component_index": component_index,
                "volume_mm3": component_volume if math.isfinite(component_volume) else None,
                "aabb_mm": _bbox_values(component, brep_bounds_mm),
                "topology": _solid_topology_report(component),
            }
        )
        if not math.isfinite(component_volume) or component_volume <= 0:
            reference_blockers.append(f"ROI_REFERENCE_COMPONENT_{component_index}_NONPOSITIVE_VOLUME")
        if reference_component_records[-1]["topology"]["status"] != "PASS":
            reference_blockers.append(f"ROI_REFERENCE_COMPONENT_{component_index}_TOPOLOGY_BLOCKED")
    reference_components_volume_mm3 = sum(
        row["volume_mm3"] for row in reference_component_records if row["volume_mm3"] is not None
    )
    if len(roi_common_solids) == 0 and roi_common_volume_mm3 > 0:
        reference_blockers.append("ROI_COMMON_POSITIVE_VOLUME_WITHOUT_SOLID_COMPONENTS")
    reference_components_match = (
        not reference_components_truncated
        and math.isfinite(roi_common_volume_mm3)
        and abs(reference_components_volume_mm3 - roi_common_volume_mm3)
        <= _volume_tolerance(roi_common_volume_mm3)
    )
    if not reference_components_match:
        reference_blockers.append("ROI_COMMON_AND_COMPONENT_VOLUMES_DISAGREE_OR_UNKNOWN")
    reference_volume_mm3 = reference_components_volume_mm3

    cells: list[dict[str, Any]] = []
    points: list[dict[str, Any]] = []
    blocking_reasons = list(reference_blockers)
    sum_cell_common_volume_mm3 = 0.0
    sum_component_volume_mm3 = 0.0
    grid_volume_complete = True
    point_id = 1

    def write_partial(status: str) -> None:
        _write_stage_json(
            stage_dir,
            "cell_ledger.partial.json",
            {
                "schema": "pgfrag.volume_ingress.cell_ledger.v1",
                "status": status,
                "region_mode": roi_mode,
                "roi_aabb_mm": list(roi_bounds),
                "source_brep_aabb_mm": list(source_bounds),
                "cell_size_mm": cell_size,
                "candidate_cell_count": candidate_count,
                "processed_cell_count": len(cells),
                "reference_roi_common": {
                    "raw_volume_mm3": roi_common_volume_mm3,
                    "raw_connected_solid_count": len(roi_common_solids),
                    "components": reference_component_records,
                    "components_truncated_by_cap": reference_components_truncated,
                },
                "frozen_input_payload_hashes_before": input_hashes_before,
                "cells": cells,
                "blocking_reasons_so_far": list(blocking_reasons),
            },
        )

    def write_cell_progress(
        status: str,
        cell_record: dict[str, Any],
        current_component: dict[str, Any] | None = None,
    ) -> None:
        previous_ledger = stage_dir / "cell_ledger.partial.json"
        _write_stage_json(
            stage_dir,
            "cell_in_progress.json",
            {
                "schema": "pgfrag.volume_ingress.cell_progress.v1",
                "status": status,
                "frozen_input_payload_hashes_before": input_hashes_before,
                "source_brep_aabb_mm": list(source_bounds),
                "roi_aabb_mm": list(roi_bounds),
                "cell_size_mm": cell_size,
                "previous_cell_ledger": {
                    "path": "cell_ledger.partial.json",
                    "sha256": _sha256(previous_ledger) if previous_ledger.is_file() else None,
                    "completed_cell_count": len(cells),
                },
                "current_cell": cell_record,
                "current_component": current_component,
                "interpretation": (
                    "Durable diagnostic snapshot of the active or most recently completed cell; "
                    "not a qualified point or native solver record."
                ),
            },
        )

    write_partial("RUNNING")
    for k, (z0, z1) in enumerate(z_cells):
        for j, (y0, y1) in enumerate(y_cells):
            for i, (x0, x1) in enumerate(x_cells):
                _check_time(started)
                cell_index = len(cells)
                cell_bounds = (x0, x1, y0, y1, z0, z1)
                cell_record: dict[str, Any] = {
                    "grid_cell_index": cell_index,
                    "grid_ijk": [i, j, k],
                    "cell_aabb_mm": list(cell_bounds),
                    "common_operation_status": "NOT_STARTED",
                    "raw_common_volume_mm3": None,
                    "raw_connected_solid_count": None,
                    "common_face_count": None,
                    "component_records": [],
                    "candidate_point_ids": [],
                    "cell_status": "UNKNOWN",
                }
                write_cell_progress("CELL_BOOLEAN_STARTED", cell_record)
                try:
                    common = source_solid.intersect(_make_box(cq, cell_bounds))
                    cell_record["common_operation_status"] = "BOOLEAN_RETURNED"
                    write_cell_progress("CELL_BOOLEAN_RETURNED", cell_record)
                except Exception as exc:
                    cell_record["common_operation_status"] = "BOOLEAN_EXCEPTION"
                    cell_record["cell_status"] = "BLOCKED_BOOLEAN_EXCEPTION"
                    cell_record["failure_reason"] = f"{type(exc).__name__}: {exc}"
                    grid_volume_complete = False
                    blocking_reasons.append(f"CELL_{cell_index}_BOOLEAN_EXCEPTION")
                    cells.append(cell_record)
                    write_partial("RUNNING_WITH_BLOCKED_CELL")
                    write_cell_progress("CELL_BOOLEAN_BLOCKED", cell_record)
                    continue

                try:
                    common_volume = float(common.Volume())
                    cell_record["raw_common_volume_mm3"] = (
                        common_volume if math.isfinite(common_volume) else None
                    )
                    write_cell_progress("COMMON_VOLUME_READ", cell_record)
                except Exception as exc:
                    cell_record["common_operation_status"] = "COMMON_INSPECTION_EXCEPTION"
                    cell_record["cell_status"] = "BLOCKED_COMMON_INSPECTION_EXCEPTION"
                    cell_record["failure_reason"] = f"{type(exc).__name__}: {exc}"
                    grid_volume_complete = False
                    blocking_reasons.append(f"CELL_{cell_index}_COMMON_INSPECTION_EXCEPTION")
                    cells.append(cell_record)
                    write_partial("RUNNING_WITH_BLOCKED_CELL")
                    write_cell_progress("COMMON_VOLUME_BLOCKED", cell_record)
                    continue
                try:
                    common_solids = list(common.Solids())
                    cell_record["raw_connected_solid_count"] = len(common_solids)
                    write_cell_progress("COMMON_SOLIDS_ENUMERATED", cell_record)
                except Exception as exc:
                    cell_record["common_operation_status"] = "COMMON_SOLID_ENUMERATION_EXCEPTION"
                    cell_record["cell_status"] = "BLOCKED_COMMON_SOLID_ENUMERATION_EXCEPTION"
                    cell_record["failure_reason"] = f"{type(exc).__name__}: {exc}"
                    grid_volume_complete = False
                    blocking_reasons.append(f"CELL_{cell_index}_COMMON_SOLID_ENUMERATION_EXCEPTION")
                    cells.append(cell_record)
                    write_partial("RUNNING_WITH_BLOCKED_CELL")
                    write_cell_progress("COMMON_SOLID_ENUMERATION_BLOCKED", cell_record)
                    continue
                try:
                    common_faces_count = len(common.Faces())
                    cell_record["common_face_count"] = common_faces_count
                    face_count_known = True
                except Exception as exc:
                    common_faces_count = None
                    face_count_known = False
                    cell_record["common_face_count_status"] = "UNKNOWN"
                    cell_record["common_face_enumeration_failure"] = (
                        f"{type(exc).__name__}: {exc}"
                    )
                    grid_volume_complete = False
                    blocking_reasons.append(f"CELL_{cell_index}_COMMON_FACE_ENUMERATION_EXCEPTION")
                try:
                    common_valid = bool(common.isValid())
                    common_valid_known = True
                except Exception as exc:
                    common_valid = False
                    common_valid_known = False
                    cell_record["common_validity_status"] = "UNKNOWN"
                    cell_record["common_validity_failure"] = f"{type(exc).__name__}: {exc}"
                    grid_volume_complete = False
                    blocking_reasons.append(f"CELL_{cell_index}_COMMON_VALIDITY_EXCEPTION")
                cell_record["common_operation_status"] = "COMPLETED"
                cell_record["raw_connected_solid_count"] = len(common_solids)
                write_cell_progress("COMMON_INSPECTION_COMPLETED", cell_record)
                _check_time(started)
                if not math.isfinite(common_volume) or common_volume < 0:
                    cell_record["cell_status"] = "BLOCKED_COMMON_VOLUME_INVALID"
                    blocking_reasons.append(f"CELL_{cell_index}_COMMON_VOLUME_INVALID")
                    grid_volume_complete = False
                elif not common_valid_known or not common_valid:
                    cell_record["cell_status"] = "BLOCKED_COMMON_INVALID"
                    blocking_reasons.append(f"CELL_{cell_index}_COMMON_INVALID")
                    grid_volume_complete = False
                elif not face_count_known:
                    cell_record["cell_status"] = "BLOCKED_COMMON_FACE_COUNT_UNKNOWN"
                    grid_volume_complete = False
                elif common_faces_count > MAX_COMMON_FACES_PER_CELL:
                    cell_record["cell_status"] = "BLOCKED_COMMON_FACE_CAP"
                    blocking_reasons.append(f"CELL_{cell_index}_COMMON_FACE_CAP")
                    grid_volume_complete = False
                else:
                    sum_cell_common_volume_mm3 += common_volume

                if len(common_solids) > MAX_CONNECTED_SOLIDS_PER_CELL:
                    cell_record["cell_status"] = "BLOCKED_COMPONENT_CAP"
                    cell_record["component_records_truncated_by_cap"] = True
                    blocking_reasons.append(f"CELL_{cell_index}_COMPONENT_CAP")
                    grid_volume_complete = False
                else:
                    cell_record["component_records_truncated_by_cap"] = False

                positive_component_sum = 0.0
                ordered_solids: list[tuple[int, Any]] = list(enumerate(common_solids))
                try:
                    ordered_solids = sorted(
                        ordered_solids,
                        key=lambda item: _component_sort_key(item[1], brep_bounds_mm),
                    )
                except Exception as exc:
                    cell_record["component_ordering_status"] = "UNKNOWN"
                    cell_record["ordering_failure"] = f"{type(exc).__name__}: {exc}"
                    blocking_reasons.append(f"CELL_{cell_index}_COMPONENT_ORDER_UNKNOWN")
                    grid_volume_complete = False
                for output_component_index, (raw_component_index, component) in enumerate(
                    ordered_solids[:MAX_CONNECTED_SOLIDS_PER_CELL]
                ):
                    component_record: dict[str, Any] = {
                        "raw_component_index": raw_component_index,
                        "intersection_component_index": output_component_index,
                        "raw_volume_mm3": None,
                        "component_aabb_mm": None,
                        "centroid_mm": None,
                        "centroid_source_classification": "NOT_CLASSIFIED",
                        "component_topology": None,
                        "candidate_point_id": None,
                        "candidate_point_status": "NOT_EMITTED",
                    }
                    write_cell_progress(
                        "COMPONENT_INSPECTION_STARTED", cell_record, component_record
                    )
                    try:
                        component_volume = float(component.Volume())
                        component_record["raw_volume_mm3"] = (
                            component_volume if math.isfinite(component_volume) else None
                        )
                        component_bounds = _bbox_values(component, brep_bounds_mm)
                        component_record["component_aabb_mm"] = list(component_bounds)
                        topology = _solid_topology_report(component)
                        component_record["component_topology"] = topology
                        write_cell_progress(
                            "COMPONENT_GEOMETRY_INSPECTED", cell_record, component_record
                        )
                    except Exception as exc:
                        component_record["component_inspection_failure"] = (
                            f"{type(exc).__name__}: {exc}"
                        )
                        component_record["candidate_point_status"] = "BLOCKED_COMPONENT_INSPECTION"
                        cell_record["component_records"].append(component_record)
                        blocking_reasons.append(
                            f"CELL_{cell_index}_COMPONENT_{raw_component_index}_INSPECTION_EXCEPTION"
                        )
                        write_cell_progress(
                            "COMPONENT_INSPECTION_BLOCKED", cell_record, component_record
                        )
                        continue

                    if not math.isfinite(component_volume) or component_volume <= 0:
                        component_record["candidate_point_status"] = "BLOCKED_NONPOSITIVE_VOLUME"
                        cell_record["component_records"].append(component_record)
                        blocking_reasons.append(
                            f"CELL_{cell_index}_COMPONENT_{raw_component_index}_NONPOSITIVE_VOLUME"
                        )
                        write_cell_progress(
                            "COMPONENT_VOLUME_BLOCKED", cell_record, component_record
                        )
                        continue
                    positive_component_sum += component_volume
                    try:
                        center = tuple(float(value) for value in component.Center().toTuple())
                    except Exception as exc:
                        component_record["centroid_source_classification"] = "UNKNOWN"
                        component_record["centroid_failure"] = f"{type(exc).__name__}: {exc}"
                        component_record["candidate_point_status"] = "BLOCKED_CENTROID_EXCEPTION"
                        cell_record["component_records"].append(component_record)
                        blocking_reasons.append(
                            f"CELL_{cell_index}_COMPONENT_{raw_component_index}_CENTROID_EXCEPTION"
                        )
                        write_cell_progress(
                            "COMPONENT_CENTROID_BLOCKED", cell_record, component_record
                        )
                        continue
                    component_record["centroid_mm"] = list(center)
                    write_cell_progress(
                        "COMPONENT_CENTROID_READ", cell_record, component_record
                    )
                    if len(center) != 3 or not all(math.isfinite(value) for value in center):
                        component_record["centroid_source_classification"] = "UNKNOWN"
                        component_record["candidate_point_status"] = "BLOCKED_CENTROID_NONFINITE"
                        cell_record["component_records"].append(component_record)
                        blocking_reasons.append(
                            f"CELL_{cell_index}_COMPONENT_{raw_component_index}_CENTROID_NONFINITE"
                        )
                        write_cell_progress(
                            "COMPONENT_CENTROID_BLOCKED", cell_record, component_record
                        )
                        continue
                    classification = _classify_point_in_source(source_solid, center)
                    component_record["centroid_source_classification"] = classification
                    write_cell_progress(
                        "COMPONENT_CENTROID_CLASSIFIED", cell_record, component_record
                    )
                    if topology["status"] != "PASS":
                        component_record["candidate_point_status"] = "BLOCKED_COMPONENT_TOPOLOGY"
                        cell_record["component_records"].append(component_record)
                        blocking_reasons.append(
                            f"CELL_{cell_index}_COMPONENT_{raw_component_index}_TOPOLOGY"
                        )
                        write_cell_progress(
                            "COMPONENT_TOPOLOGY_BLOCKED", cell_record, component_record
                        )
                        continue
                    if classification != "IN":
                        component_record["candidate_point_status"] = (
                            f"BLOCKED_CENTROID_{classification}"
                        )
                        cell_record["component_records"].append(component_record)
                        blocking_reasons.append(
                            f"CELL_{cell_index}_COMPONENT_{raw_component_index}_CENTROID_{classification}"
                        )
                        write_cell_progress(
                            "COMPONENT_CENTROID_CLASSIFICATION_BLOCKED",
                            cell_record,
                            component_record,
                        )
                        continue
                    if len(points) >= MAX_MATERIAL_POINTS:
                        component_record["candidate_point_status"] = "BLOCKED_POINT_CAP"
                        cell_record["component_records"].append(component_record)
                        blocking_reasons.append("MATERIAL_POINT_CAP_EXCEEDED")
                        write_cell_progress(
                            "COMPONENT_POINT_CAP_BLOCKED", cell_record, component_record
                        )
                        continue

                    volume_m3 = component_volume * 1e-9
                    point = {
                        "point_id": point_id,
                        "position_m": [value * 1e-3 for value in center],
                        "cell_volume_m3": volume_m3,
                        "mass_kg": context["rho_kg_m3"] * volume_m3,
                        "source": {
                            "input_brep_solid_index": 0,
                            "grid_cell_index": cell_index,
                            "grid_ijk": [i, j, k],
                            "intersection_component_index": output_component_index,
                            "raw_component_index": raw_component_index,
                            "cell_aabb_mm": list(cell_bounds),
                            "component_aabb_mm": list(component_bounds),
                            "source_semantic_region_code": None,
                            "source_semantic_region_status": "UNKNOWN_NOT_PROJECTED",
                            "surface_parent_triangle": None,
                            "surface_parent_mapping_status": "UNKNOWN_NOT_PROJECTED",
                        },
                    }
                    points.append(point)
                    cell_record["candidate_point_ids"].append(point_id)
                    component_record["candidate_point_id"] = point_id
                    component_record["candidate_point_status"] = "CANDIDATE_CENTROID_IN"
                    point_id += 1
                    cell_record["component_records"].append(component_record)

                sum_component_volume_mm3 += positive_component_sum
                if abs(positive_component_sum - common_volume) > _volume_tolerance(common_volume):
                    cell_record["cell_status"] = "BLOCKED_COMPONENT_VOLUME_MISMATCH"
                    blocking_reasons.append(f"CELL_{cell_index}_COMPONENT_VOLUME_MISMATCH")
                    grid_volume_complete = False
                elif cell_record["cell_status"].startswith("BLOCKED"):
                    pass
                elif common_volume == 0 and not common_solids:
                    cell_record["cell_status"] = "ZERO_INTERSECTION"
                elif any(
                    row["candidate_point_status"] != "CANDIDATE_CENTROID_IN"
                    for row in cell_record["component_records"]
                ):
                    cell_record["cell_status"] = "BLOCKED_COMPONENT_OR_CENTROID"
                else:
                    cell_record["cell_status"] = "POSITIVE_COMPONENTS_CENTROIDS_IN"
                cells.append(cell_record)
                write_partial(
                    "RUNNING_WITH_BLOCKED_CELL"
                    if cell_record["cell_status"].startswith("BLOCKED")
                    else "RUNNING"
                )
                write_cell_progress("CELL_COMPLETED", cell_record)

    volume_conservation_status = "UNKNOWN_INCOMPLETE_CELL_GRID"
    conservation_error_mm3: float | None = None
    conservation_tolerance_mm3 = _volume_tolerance(reference_volume_mm3)
    if grid_volume_complete and len(cells) == candidate_count and reference_components_match:
        conservation_error_mm3 = sum_component_volume_mm3 - reference_volume_mm3
        if abs(conservation_error_mm3) <= conservation_tolerance_mm3:
            volume_conservation_status = "PASS_LOCAL_ROI_BREP_VOLUME_RECONCILIATION"
        else:
            volume_conservation_status = "FAIL_LOCAL_ROI_BREP_VOLUME_RECONCILIATION"
            blocking_reasons.append("LOCAL_ROI_VOLUME_RECONCILIATION_FAILED")

    if len(cells) != candidate_count:
        blocking_reasons.append("GRID_DID_NOT_COMPLETE_ALL_96_CELLS")
    if not points:
        blocking_reasons.append("NO_CENTROID_IN_MATERIAL_POINTS")

    input_hashes_after = {
        "manifest.json": _sha256(context["asset_dir"] / "manifest.json")
    }
    for relative in sorted(manifest["files"]):
        input_hashes_after[relative] = _sha256(
            _package_payload_path(context["asset_dir"], relative)
        )
    if input_hashes_before != input_hashes_after:
        blocking_reasons.append("FROZEN_INPUT_PACKAGE_CHANGED_DURING_READ")

    unique_reasons = sorted(set(blocking_reasons))
    candidate_pass = not unique_reasons and volume_conservation_status.startswith("PASS")
    wall_mm = context["declared_wall_parameter_mm"]
    h_over_declared_wall = cell_size / wall_mm if wall_mm else None
    cell_ledger = {
        "schema": "pgfrag.volume_ingress.cell_ledger.v1",
        "status": "PASS_LOCAL_DIAGNOSTIC_RAW_LEDGER" if candidate_pass else "BLOCKED_RAW_DIAGNOSTIC_LEDGER",
        "region_mode": roi_mode,
        "roi_aabb_mm": list(roi_bounds),
        "source_brep_aabb_mm": list(source_bounds),
        "source_brep_bounds_policy": brep_face_bounds_policy_record(),
        "cell_size_mm": cell_size,
        "candidate_cell_count": candidate_count,
        "processed_cell_count": len(cells),
        "material_point_candidate_count": len(points),
        "reference_roi_common": {
            "raw_common_volume_mm3": roi_common_volume_mm3,
            "raw_connected_solid_count": len(roi_common_solids),
            "raw_component_volumes_mm3": [
                row["volume_mm3"] for row in reference_component_records
            ],
            "components_truncated_by_cap": reference_components_truncated,
            "components": reference_component_records,
        },
        "sum_cell_common_volumes_mm3": sum_cell_common_volume_mm3,
        "sum_connected_component_volumes_mm3": sum_component_volume_mm3,
        "volume_conservation_status": volume_conservation_status,
        "cells": cells,
        "blocking_reasons": unique_reasons,
    }
    write_partial("COMPLETE_RAW_LEDGER" if candidate_pass else "COMPLETE_BLOCKED_RAW_LEDGER")
    _write_stage_json(stage_dir, "cell_ledger.json", cell_ledger)

    pointset: dict[str, Any] | None = None
    if candidate_pass:
        pointset = {
            "schema": "pgfrag.volume_ingress.pointset_candidate.v1",
            "status": "PASS_LOCAL_BREP_INGRESS_DIAGNOSTIC_ONLY",
            "units": {
                "position": "m",
                "cell_volume": "m^3",
                "mass": "kg",
                "source_aabb": "mm",
            },
            "frame": "right_handed_z_up",
            "point_count": len(points),
            "points": points,
            "provenance": {
                "vf_asset_schema": manifest["schema"],
                "vf_asset_manifest_sha256": input_hashes_before["manifest.json"],
                "canonical_step_sha256": input_hashes_before[context["step_relative"]],
                "input_payload_sha256": input_hashes_before,
                "rho_kg_m3": context["rho_kg_m3"],
                "mechanical_calibration_status": manifest["material"]["calibration_status"],
                "source_coordinate_frame": manifest["frame"],
                "source_length_unit": manifest["units"]["length"],
                "canonical_geometry_provider": "OpenCascade B-rep solid-cell intersection; no h^3 approximation",
                "source_solid_count": len(source_solids),
                "source_material_solid_index": 0,
                "source_semantic_codes_available": manifest.get("provenance", {}).get("spec", {}).get(
                    "surface_semantics", {}
                ),
                "source_semantic_projection_status": "UNKNOWN_NOT_PROJECTED",
                "surface_parent_mapping_status": "UNKNOWN_NOT_PROJECTED",
                "native_pgfrag_integration_status": "NOT_CONNECTED",
                "candidate_pair_classification_status": "NOT_RUN_SEPARATE_GATE",
                "fragment_extraction_status": "NOT_IMPLEMENTED_BY_THIS_WORKER",
            },
        }
        _write_stage_json(stage_dir, "material_points.json", pointset)

    report = {
        "schema": "pgfrag.volume_ingress.report.v1",
        "status": (
            "PASS_LOCAL_BREP_INGRESS_DIAGNOSTIC_ONLY"
            if candidate_pass
            else "BLOCKED_LOCAL_BREP_INGRESS_DIAGNOSTIC"
        ),
        "interpretation": (
            "Local B-rep ingress mechanism diagnostic only; never native graph initialization, "
            "physics qualification, or whole-mother discretization."
        ),
        "input": {
            "asset_relative_path": context["asset_dir"].relative_to(REPO_ROOT).as_posix(),
            "manifest_sha256_before": input_hashes_before["manifest.json"],
            "manifest_sha256_after": input_hashes_after["manifest.json"],
            "canonical_step_sha256_before": input_hashes_before[context["step_relative"]],
            "canonical_step_sha256_after": input_hashes_after[context["step_relative"]],
            "payload_hashes_before": input_hashes_before,
            "payload_hashes_after": input_hashes_after,
            "source_solid_count": len(source_solids),
            "source_topology": source_topology,
            "source_brep_aabb_mm": list(source_bounds),
            "source_brep_bounds_policy": brep_face_bounds_policy_record(),
        },
        "discretization": {
            "region_mode": roi_mode,
            "roi_aabb_mm": list(roi_bounds),
            "roi_intentionally_extends_outside_source_xmax": True,
            "cell_size_mm": cell_size,
            "candidate_cell_count": candidate_count,
            "processed_cell_count": len(cells),
            "material_point_candidate_count": len(points),
            "declared_wall_parameter_mm": wall_mm,
            "h_over_declared_wall_parameter": h_over_declared_wall,
            "nominal_cells_per_declared_wall": wall_mm / cell_size if wall_mm else None,
            "resolution_qualification": "LOCAL_DIAGNOSTIC_ONLY_NOT_PHYSICS_READY",
            "whole_mother_coverage_status": "PARTIAL_ROI_ONLY",
        },
        "volume_conservation": {
            "scope": "explicit local ROI only",
            "status": volume_conservation_status,
            "reference_roi_brep_common_volume_mm3": roi_common_volume_mm3,
            "reference_roi_connected_component_volumes_mm3": [
                row["volume_mm3"] for row in reference_component_records
            ],
            "sum_cell_brep_common_volumes_mm3": sum_cell_common_volume_mm3,
            "sum_cell_connected_component_volumes_mm3": sum_component_volume_mm3,
            "difference_mm3": conservation_error_mm3,
            "tolerance_mm3": conservation_tolerance_mm3,
            "relative_tolerance": VOLUME_RELATIVE_TOLERANCE,
            "absolute_tolerance_mm3": VOLUME_ABSOLUTE_TOLERANCE_MM3,
            "method": "OpenCascade B-rep solid-cell intersection volumes; no h^3 approximation",
            "whole_mother_conservation": "UNKNOWN_NOT_RUN",
        },
        "centroid_gate": {
            "status": (
                "ALL_POSITIVE_COMPONENT_CENTROIDS_IN"
                if not any("CENTROID_" in reason for reason in unique_reasons)
                else "BLOCKED_WITH_RAW_CLASSIFICATIONS_RETAINED"
            ),
            "classifier": "OCP.BRepClass3d_SolidClassifier",
            "tolerance_mm": INSIDE_TOLERANCE_MM,
            "on_out_unknown_action": "BLOCK_CANDIDATE_NO_SURFACE_PROJECTION",
        },
        "qualification": {
            "mechanics": "UNVALIDATED_PROTOTYPE",
            "native_solver": "NOT_CONNECTED",
            "candidate_bond_mask": "NOT_CONNECTED",
            "surface_parent_mapping": "UNKNOWN_NOT_PROJECTED",
            "fragment_mapping": "UNKNOWN_NOT_DERIVED",
            "whole_geometry_wall": "NOT_ESTABLISHED",
            "physics_ready": False,
        },
        "blocking_reasons": unique_reasons,
        "limits": {
            "max_candidate_cells": MAX_CANDIDATE_CELLS,
            "max_material_points": MAX_MATERIAL_POINTS,
            "max_connected_solids_per_cell": MAX_CONNECTED_SOLIDS_PER_CELL,
            "max_reference_roi_components": MAX_REFERENCE_ROI_COMPONENTS,
            "max_common_faces_per_cell": MAX_COMMON_FACES_PER_CELL,
            "max_output_bytes": MAX_OUTPUT_BYTES,
            "parent_receipt_reserve_bytes": PARENT_RECEIPT_RESERVE_BYTES,
            "worker_timeout_seconds": WORKER_TIMEOUT_SECONDS,
        },
        "runtime": {
            "python_executable": sys.executable,
            "python_version": sys.version.split()[0],
            "cadquery_version": cq.__version__,
            "provider": "CadQuery/OpenCascade",
            "elapsed_seconds": time.monotonic() - started,
        },
    }
    _write_final_report(stage_dir, report)
    return {
        "status": report["status"],
        "candidate_pass": candidate_pass,
        "point_count": len(points),
        "candidate_cell_count": candidate_count,
        "processed_cell_count": len(cells),
        "reference_roi_brep_common_volume_mm3": roi_common_volume_mm3,
        "sum_cell_connected_common_volumes_mm3": sum_component_volume_mm3,
        "blocking_reasons": unique_reasons,
        "output_bytes": _stage_inventory(stage_dir)["total_bytes"],
    }


def _write_final_report(stage_dir: Path, report: dict[str, Any]) -> None:
    non_report_files = [
        path for path in stage_dir.iterdir()
        if path.is_file() and path.name != "ingress_report.json" and not path.name.endswith(".tmp")
    ]
    report["output"] = {
        "total_bytes": 0,
        "files": {
            path.name: {"bytes": path.stat().st_size, "sha256": _sha256(path)}
            for path in sorted(non_report_files, key=lambda item: item.name)
        },
    }
    non_report_bytes = sum(path.stat().st_size for path in non_report_files)
    report_payload = b""
    for _ in range(5):
        report["output"]["total_bytes"] = 0
        report_payload = _json_bytes(report)
        reported_total = non_report_bytes + len(report_payload)
        report["output"]["total_bytes"] = reported_total
        report_payload = _json_bytes(report)
        if non_report_bytes + len(report_payload) == reported_total:
            break
    else:
        raise IngressBlocked("could not stabilize report-inclusive output byte count")
    _write_stage_json(stage_dir, "ingress_report.json", report)


def _mark_stage_blocked(stage: Path, cause: dict[str, Any]) -> list[str]:
    errors: list[str] = []
    parent_status = str(cause.get("status", "BLOCKED_PARENT_EXECUTION"))
    failure_summary = {
        "status": parent_status,
        "failure_class": cause.get("failure_class"),
        "reason": cause.get("reason"),
        "timed_out": cause.get("timed_out"),
        "child_exit_code": cause.get("child_exit_code"),
        "launch_error": cause.get("launch_error"),
        "child_alive_after_stop": cause.get("child_alive_after_stop"),
        "process_stop_errors": cause.get("process_stop_errors"),
    }
    ledger_blocked_status = (
        "BLOCKED_WORKER_EXCEPTION"
        if parent_status == "BLOCKED_WORKER_EXCEPTION"
        else "BLOCKED_PARENT_EXECUTION"
    )
    targets = (
        ("cell_ledger.partial.json", ledger_blocked_status),
        ("cell_ledger.json", ledger_blocked_status),
        ("cell_in_progress.json", parent_status),
    )
    for filename, blocked_status in targets:
        path = stage / filename
        if not path.is_file():
            if filename != "cell_ledger.partial.json":
                continue
            record: dict[str, Any] = {
                "schema": "pgfrag.volume_ingress.cell_ledger.v1",
                "status": blocked_status,
                "processed_cell_count": 0,
                "candidate_cell_count": 96,
                "cells": [],
            }
        else:
            try:
                record = json.loads(path.read_text(encoding="utf-8"))
            except Exception as exc:
                errors.append(f"{filename} could not be read for blocked marking: {type(exc).__name__}: {exc}")
                continue
        old_status = record.get("status")
        if not str(old_status).startswith("BLOCKED"):
            record["status_before_parent_block"] = old_status
        record["status"] = blocked_status
        record["parent_execution_failure"] = failure_summary
        try:
            _write_stage_json(stage, filename, record, reserve_parent_receipt=True)
        except Exception as exc:
            errors.append(f"{filename} blocked marking failed: {type(exc).__name__}: {exc}")
    return errors


def _child_main(request: dict[str, Any], stage_path: str, result_queue: Any) -> None:
    stage = Path(stage_path)
    try:
        summary = _compute_candidate(request, stage)
        result_queue.put(
            {
                "status": "PASS" if summary["candidate_pass"] else "BLOCKED",
                "summary": summary,
            }
        )
    except BaseException as exc:
        blocking_marker_errors = _mark_stage_blocked(
            stage,
            {
                "status": "BLOCKED_WORKER_EXCEPTION",
                "failure_class": type(exc).__name__,
                "reason": str(exc),
            },
        )
        progress_path = stage / "cell_in_progress.json"
        try:
            active_cell_diagnostic = (
                json.loads(progress_path.read_text(encoding="utf-8"))
                if progress_path.is_file()
                else None
            )
        except Exception as progress_exc:
            active_cell_diagnostic = {
                "read_error": f"{type(progress_exc).__name__}: {progress_exc}"
            }
        failure = {
            "schema": "pgfrag.volume_ingress.worker_failure.v1",
            "status": "BLOCKED",
            "failure_class": type(exc).__name__,
            "reason": str(exc),
            "frozen_input_pins": {
                "asset_relative_path": "vesselforge/tmp/review_20261005/geometry_resume/migration_copy/cylinder_asset_run01",
                "manifest_sha256": EXPECTED_MANIFEST_SHA256,
                "canonical_step_sha256": EXPECTED_STEP_SHA256,
            },
            "roi_aabb_mm": list(DIAGNOSTIC_ROI_BOUNDS_MM),
            "cell_size_mm": DIAGNOSTIC_CELL_SIZE_MM,
            "stage_inventory_before_failure_receipt": _stage_inventory(stage),
            "active_cell_diagnostic": active_cell_diagnostic,
            "blocked_marker_errors": blocking_marker_errors,
            "raw_partial_ledger": (
                {
                    "path": "cell_ledger.partial.json",
                    "status": _read_json_record(stage / "cell_ledger.partial.json").get("status")
                    if (stage / "cell_ledger.partial.json").is_file()
                    else "UNKNOWN_NOT_WRITTEN",
                    "sha256": (
                        _sha256(stage / "cell_ledger.partial.json")
                        if (stage / "cell_ledger.partial.json").is_file()
                        else None
                    ),
                }
            ),
        }
        failure_receipt_write_error: str | None = None
        try:
            _write_stage_json(
                stage,
                "worker_failure.json",
                failure,
                reserve_parent_receipt=True,
            )
        except Exception as receipt_exc:
            failure_receipt_write_error = f"{type(receipt_exc).__name__}: {receipt_exc}"
        result_queue.put(
            {
                "status": "BLOCKED",
                "error_type": type(exc).__name__,
                "reason": str(exc),
                "stage_retained": True,
                "worker_failure_receipt_written": failure_receipt_write_error is None,
                "worker_failure_receipt_write_error": failure_receipt_write_error,
                "blocked_marker_errors": blocking_marker_errors,
            }
        )


def _safe_tmp_path(path: Path) -> Path:
    resolved = path.resolve(strict=False)
    if not resolved.is_relative_to(TMP_ROOT):
        raise IngressBlocked("output and staging paths must remain under this experiment's tmp")
    return resolved


def _stop_child_process(process: Any) -> tuple[bool, int | None, list[str]]:
    errors: list[str] = []
    for action in ("terminate", "kill"):
        try:
            if not process.is_alive():
                break
            if action == "kill" and not hasattr(process, "kill"):
                errors.append("child remains alive and Process.kill is unavailable")
                break
            getattr(process, action)()
            process.join(2)
        except Exception as exc:
            errors.append(f"{action} failed: {type(exc).__name__}: {exc}")
    try:
        still_alive = bool(process.is_alive())
    except Exception as exc:
        still_alive = True
        errors.append(f"alive-state read failed: {type(exc).__name__}: {exc}")
    try:
        exit_code = process.exitcode
    except Exception as exc:
        exit_code = None
        errors.append(f"exit-code read failed: {type(exc).__name__}: {exc}")
    return still_alive, exit_code, errors


def run_bounded(request: dict[str, Any], output_dir: Path) -> dict[str, Any]:
    TMP_ROOT.mkdir(parents=True, exist_ok=True)
    output_dir = _safe_tmp_path(output_dir)
    if output_dir.exists():
        raise IngressBlocked("output path already exists; refusing to overwrite")
    if float(request.get("cell_size_mm", math.nan)) != DIAGNOSTIC_CELL_SIZE_MM:
        raise IngressBlocked("request must use the frozen root-selected h=0.5 mm")
    if tuple(request.get("roi_bounds_mm", ())) != DIAGNOSTIC_ROI_BOUNDS_MM:
        raise IngressBlocked("request must use the frozen root-selected diagnostic ROI")
    stage = Path(tempfile.mkdtemp(prefix=".forge_ingress_stage_", dir=TMP_ROOT))
    request_record = {
        "schema": "pgfrag.volume_ingress.request.v1",
        "status": "SOURCE_ONLY_PILOT_REQUEST",
        "asset_relative_path": "vesselforge/tmp/review_20261005/geometry_resume/migration_copy/cylinder_asset_run01",
        "expected_manifest_sha256": EXPECTED_MANIFEST_SHA256,
        "expected_canonical_step_sha256": EXPECTED_STEP_SHA256,
        "cell_size_mm": DIAGNOSTIC_CELL_SIZE_MM,
        "roi_aabb_mm": list(DIAGNOSTIC_ROI_BOUNDS_MM),
        "expected_grid_shape": [6, 4, 4],
        "expected_candidate_cell_count": 96,
        "limits": {
            "max_candidate_cells": MAX_CANDIDATE_CELLS,
            "max_material_points": MAX_MATERIAL_POINTS,
            "max_connected_solids_per_cell": MAX_CONNECTED_SOLIDS_PER_CELL,
            "max_reference_roi_components": MAX_REFERENCE_ROI_COMPONENTS,
            "max_common_faces_per_cell": MAX_COMMON_FACES_PER_CELL,
            "max_output_bytes": MAX_OUTPUT_BYTES,
            "hard_timeout_seconds": WORKER_TIMEOUT_SECONDS,
        },
        "runtime_required": str(EXPECTED_VF_PYTHON),
        "native_integration": "NOT_CONNECTED",
    }
    request_write_error: str | None = None
    try:
        _write_stage_json(stage, "request.json", request_record)
    except Exception as exc:
        request_write_error = f"{type(exc).__name__}: {exc}"
    started = time.monotonic()
    timed_out = False
    launch_error = request_write_error
    process: Any | None = None
    result_queue: Any | None = None
    process_started = False
    child_result: dict[str, Any] | None = None
    exit_code: int | None = None
    child_alive_after_stop = False
    process_stop_errors: list[str] = []
    if launch_error is None:
        try:
            context = mp.get_context("spawn")
            result_queue = context.Queue(maxsize=1)
            process = context.Process(
                target=_child_main,
                args=(request, str(stage), result_queue),
                name="forge-volume-ingress-source-only",
            )
            process.start()
            process_started = True
            process.join(WORKER_TIMEOUT_SECONDS)
            if process.is_alive():
                timed_out = True
                child_alive_after_stop, exit_code, process_stop_errors = _stop_child_process(process)
            else:
                exit_code = process.exitcode
            if not timed_out:
                try:
                    child_result = result_queue.get(timeout=2)
                except Exception:
                    child_result = None
        except Exception as exc:
            launch_error = f"{type(exc).__name__}: {exc}"
            if process is not None and (
                process_started or getattr(process, "_popen", None) is not None
            ):
                child_alive_after_stop, exit_code, process_stop_errors = _stop_child_process(process)
    if result_queue is not None:
        try:
            result_queue.close()
        except Exception:
            pass
    elapsed = time.monotonic() - started
    child_status = (
        child_result.get("status")
        if child_result is not None
        else ("TIMEOUT" if timed_out else "NO_CHILD_RECEIPT")
    )
    failure_reason = (
        "CHILD_TIMEOUT"
        if timed_out
        else "PARENT_LAUNCH_OR_ORCHESTRATION_ERROR"
        if launch_error is not None
        else "CHILD_DID_NOT_COMPLETE_PASS"
    )
    marker_errors: list[str] = []
    run_failed = (
        timed_out
        or launch_error is not None
        or child_status != "PASS"
        or exit_code != 0
    )
    if run_failed:
        marker_errors = _mark_stage_blocked(
            stage,
            {
                "status": "BLOCKED_PARENT_TIMEOUT" if timed_out else "BLOCKED_PARENT_EXECUTION",
                "failure_class": failure_reason,
                "reason": launch_error or (child_result or {}).get("reason") or child_status,
                "timed_out": timed_out,
                "child_exit_code": exit_code,
                "launch_error": launch_error,
                "child_alive_after_stop": child_alive_after_stop,
                "process_stop_errors": process_stop_errors,
            },
        )
    parent_receipt = {
        "schema": "pgfrag.volume_ingress.parent_execution_receipt.v1",
        "status": (
            "CHILD_COMPLETED"
            if child_status == "PASS" and exit_code == 0
            else "BLOCKED_STAGE_RETAINED"
        ),
        "child_status": child_status,
        "child_exit_code": exit_code,
        "timed_out": timed_out,
        "launch_error": launch_error,
        "elapsed_seconds": elapsed,
        "hard_timeout_seconds": WORKER_TIMEOUT_SECONDS,
        "candidate_run_promoted": False,
        "child_alive_after_stop": child_alive_after_stop,
        "process_stop_errors": process_stop_errors,
        "stage_blocking_marker_errors": marker_errors,
        "partial_ledger_status": (
            _read_json_record(stage / "cell_ledger.partial.json").get("status")
            if (stage / "cell_ledger.partial.json").is_file()
            else "UNKNOWN_NOT_WRITTEN"
        ),
        "cell_progress_receipt": (
            {
                "path": "cell_in_progress.json",
                "sha256": _sha256(stage / "cell_in_progress.json"),
                "status": _read_json_record(stage / "cell_in_progress.json").get("status"),
            }
            if (stage / "cell_in_progress.json").is_file()
            else None
        ),
        "stage_inventory_before_parent_receipt": _stage_inventory(stage),
        "child_result": child_result,
        "request_record_sha256": hashlib.sha256(_json_bytes(request_record)).hexdigest(),
    }
    try:
        _write_stage_json(
            stage,
            "parent_execution_receipt.json",
            parent_receipt,
            reserve_parent_receipt=False,
        )
    except Exception as exc:
        run_failed = True
        marker_errors.extend(
            _mark_stage_blocked(
                stage,
                {
                    "status": "BLOCKED_PARENT_EXECUTION",
                    "failure_class": "PARENT_RECEIPT_WRITE_FAILURE",
                    "reason": f"{type(exc).__name__}: {exc}",
                    "timed_out": timed_out,
                    "child_exit_code": exit_code,
                    "launch_error": launch_error,
                    "child_alive_after_stop": child_alive_after_stop,
                    "process_stop_errors": process_stop_errors,
                },
            )
        )
        parent_receipt["receipt_write_failure"] = f"{type(exc).__name__}: {exc}"
        fallback_receipt = {
            "schema": "pgfrag.volume_ingress.parent_execution_receipt.v1",
            "status": "BLOCKED_STAGE_RETAINED",
            "child_status": child_status,
            "child_exit_code": exit_code,
            "timed_out": timed_out,
            "launch_error": launch_error,
            "elapsed_seconds": elapsed,
            "hard_timeout_seconds": WORKER_TIMEOUT_SECONDS,
            "candidate_run_promoted": False,
            "full_receipt_write_failure": f"{type(exc).__name__}: {exc}",
            "retained_stage": stage.relative_to(TMP_ROOT).as_posix(),
            "child_alive_after_stop": child_alive_after_stop,
            "process_stop_errors": process_stop_errors,
            "stage_blocking_marker_errors": marker_errors,
            "child_result": child_result,
        }
        try:
            _write_stage_json(
                stage,
                "parent_execution_receipt.json",
                fallback_receipt,
                reserve_parent_receipt=False,
            )
            parent_receipt = fallback_receipt
        except Exception as fallback_exc:
            parent_receipt["fallback_receipt_write_failure"] = (
                f"{type(fallback_exc).__name__}: {fallback_exc}"
            )
    stage_inventory_after = _stage_inventory(stage)
    if run_failed:
        return {
            "status": "BLOCKED_RAW_STAGE_RETAINED",
            "retained_stage": stage.relative_to(TMP_ROOT).as_posix(),
            "timed_out": timed_out,
            "child_exit_code": exit_code,
            "elapsed_seconds": elapsed,
            "stage_inventory": stage_inventory_after,
            "child_result": child_result,
        }
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    os.replace(stage, output_dir)
    return {
        "status": child_result["summary"]["status"],
        "output_dir": output_dir.relative_to(TMP_ROOT).as_posix(),
        **child_result["summary"],
        "parent_receipt": "parent_execution_receipt.json",
        "stage_inventory": stage_inventory_after,
    }


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="SOURCE_ONLY bounded B-rep common material-point candidate"
    )
    parser.add_argument(
        "--run-bounded-pilot",
        action="store_true",
        help="explicitly run the frozen 96-cell diagnostic after root source review",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=TMP_ROOT / "candidate_run",
    )
    args = parser.parse_args(argv)
    if not args.run_bounded_pilot:
        parser.error("SOURCE_ONLY: wait for root review before passing --run-bounded-pilot")
    actual_python = Path(sys.executable).resolve()
    if os.path.normcase(str(actual_python)) != os.path.normcase(str(EXPECTED_VF_PYTHON)):
        parser.error(f"run only from the frozen VesselForge OCP Python: {EXPECTED_VF_PYTHON}")
    if not EXPECTED_VF_PYTHON.is_file():
        parser.error(f"frozen VesselForge Python runtime is missing: {EXPECTED_VF_PYTHON}")
    request = {
        "cell_size_mm": DIAGNOSTIC_CELL_SIZE_MM,
        "roi_min_mm": [
            DIAGNOSTIC_ROI_BOUNDS_MM[0],
            DIAGNOSTIC_ROI_BOUNDS_MM[2],
            DIAGNOSTIC_ROI_BOUNDS_MM[4],
        ],
        "roi_max_mm": [
            DIAGNOSTIC_ROI_BOUNDS_MM[1],
            DIAGNOSTIC_ROI_BOUNDS_MM[3],
            DIAGNOSTIC_ROI_BOUNDS_MM[5],
        ],
        "roi_bounds_mm": list(DIAGNOSTIC_ROI_BOUNDS_MM),
    }
    try:
        result = run_bounded(request, args.output_dir)
    except Exception as exc:
        print(
            json.dumps(
                {"status": "BLOCKED", "error_type": type(exc).__name__, "reason": str(exc)},
                ensure_ascii=False,
                sort_keys=True,
            ),
            file=sys.stderr,
        )
        return 2
    print(json.dumps(result, ensure_ascii=False, sort_keys=True))
    return 0 if result["status"] == "PASS_LOCAL_BREP_INGRESS_DIAGNOSTIC_ONLY" else 2


if __name__ == "__main__":
    raise SystemExit(main())
