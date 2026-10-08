"""Strict, source-only PG fracture artifact contract for A12 integration.

This experiment defines a small export boundary. It does not run a solver, infer
fracture surfaces, qualify materials, or alter the existing PG/A12 readers.
"""

from __future__ import annotations

from dataclasses import dataclass
import base64
import hashlib
import json
import math
import os
from pathlib import Path
import re
import stat
from typing import Any, Callable, Mapping


SCHEMA = "pgfrag.fracture_artifact.v1"
OBSERVATION_SCHEMA = "a12.fracture_observation.v1"
MAX_DOCUMENT_BYTES = 32_000_000
MAX_PAYLOAD_FILE_BYTES = 8_000_000
MAX_PAYLOAD_TOTAL_BYTES = 30_000_000
MAX_ENTITY_ROWS = 100_000
MAX_SURFACE_TRIANGLES = 200_000
SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")
OBSERVATION_ID_RE = re.compile(r"obs_[0-9a-f]{16}\Z")

QUALIFICATION_CLASSES = frozenset({
    "FIXTURE_GEOMETRY_ONLY",
    "SIMULATION_UNCALIBRATED",
    "PHYSICALLY_VALIDATED",
})

_OBSERVATION_FORBIDDEN_KEYS = frozenset({
    "motherid", "mothermeshsha256", "mothercoordinatesm", "cleangroundtruth",
    "groundtruth", "truth", "eventid", "pairid", "surfaceid",
    "reciprocalsurfaceid", "lineage", "parentfragmentid", "fragmentid",
    "matefragmentid", "r2label", "cleansurfacesha256", "answerfield",
})


@dataclass(frozen=True)
class ValidationReport:
    valid: bool
    errors: tuple[dict[str, str], ...]
    qualification_class: str | None
    physical_receipt_verified: bool
    pretraining_qualified: bool
    observation_projection_safe: bool

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema": "pgfrag.fracture_artifact_validation.v1",
            "valid": self.valid,
            "errors": [dict(error) for error in self.errors],
            "qualification_class": self.qualification_class,
            "physical_receipt_verified": self.physical_receipt_verified,
            "pretraining_qualified": self.pretraining_qualified,
            "observation_projection_safe": self.observation_projection_safe,
        }


class ArtifactContractError(ValueError):
    """Raised when strict reading or export fails closed."""

    def __init__(self, report: ValidationReport | None = None, code: str = "INVALID_ARTIFACT") -> None:
        self.report = report
        self.code = code
        if report is None:
            message = code
        else:
            message = "; ".join(f"{e['code']} at {e['path']}" for e in report.errors) or code
        super().__init__(message)


def _sha256(value: Any) -> bool:
    return isinstance(value, str) and SHA256_RE.fullmatch(value) is not None


def _strict_object_keys(
    value: Any,
    required: set[str],
    optional: set[str],
    path: str,
    errors: list[dict[str, str]],
) -> bool:
    if not isinstance(value, dict):
        errors.append({"code": "TYPE_OBJECT_REQUIRED", "path": path})
        return False
    for key in sorted(required.difference(value)):
        errors.append({"code": "MISSING_FIELD", "path": f"{path}.{key}"})
    for key in sorted(set(value).difference(required | optional)):
        errors.append({"code": "UNKNOWN_FIELD", "path": f"{path}.{key}"})
    return True


def _identifier(value: Any, path: str, errors: list[dict[str, str]]) -> bool:
    if not isinstance(value, str) or not value.strip() or value.strip() != value:
        errors.append({"code": "INVALID_IDENTIFIER", "path": path})
        return False
    return True


def _known_value(
    value: Any,
    path: str,
    expected_unit: str,
    errors: list[dict[str, str]],
    *,
    positive: bool = False,
    minimum: float | None = None,
    maximum: float | None = None,
) -> tuple[str | None, float | None]:
    if not _strict_object_keys(value, {"state"}, {"value", "unit", "reason"}, path, errors):
        return None, None
    state = value.get("state")
    if state == "KNOWN":
        if not {"value", "unit"}.issubset(value) or "reason" in value:
            errors.append({"code": "KNOWN_VALUE_SHAPE", "path": path})
            return state, None
        number = value.get("value")
        if isinstance(number, bool) or not isinstance(number, (int, float)) or not math.isfinite(number):
            errors.append({"code": "NONFINITE_OR_NONNUMERIC_VALUE", "path": f"{path}.value"})
            return state, None
        number = float(number)
        if value.get("unit") != expected_unit:
            errors.append({"code": "UNIT_MISMATCH", "path": f"{path}.unit"})
        if positive and number <= 0:
            errors.append({"code": "KNOWN_VALUE_MUST_BE_POSITIVE", "path": f"{path}.value"})
        if minimum is not None and number < minimum:
            errors.append({"code": "VALUE_BELOW_MINIMUM", "path": f"{path}.value"})
        if maximum is not None and number > maximum:
            errors.append({"code": "VALUE_ABOVE_MAXIMUM", "path": f"{path}.value"})
        return state, number
    if state in {"UNKNOWN", "NOT_APPLICABLE"}:
        if "value" in value or "unit" in value or not isinstance(value.get("reason"), str) or not value["reason"].strip():
            errors.append({"code": "EXPLICIT_UNKNOWN_REASON_REQUIRED", "path": path})
        return state, None
    errors.append({"code": "INVALID_KNOWLEDGE_STATE", "path": f"{path}.state"})
    return None, None


def _known_sha(value: Any, path: str, errors: list[dict[str, str]]) -> str | None:
    if not _strict_object_keys(value, {"state"}, {"value", "reason"}, path, errors):
        return None
    state = value.get("state")
    if state == "KNOWN":
        if set(value) != {"state", "value"}:
            errors.append({"code": "KNOWN_VALUE_SHAPE", "path": path})
            return state
        if not _sha256(value.get("value")):
            errors.append({"code": "INVALID_SHA256", "path": f"{path}.value"})
    elif state in {"UNKNOWN", "NOT_APPLICABLE"}:
        if set(value) != {"state", "reason"} or not isinstance(value.get("reason"), str) or not value["reason"].strip():
            errors.append({"code": "EXPLICIT_UNKNOWN_REASON_REQUIRED", "path": path})
    else:
        errors.append({"code": "INVALID_KNOWLEDGE_STATE", "path": f"{path}.state"})
    return state


def _scan_nonfinite(value: Any, path: str, errors: list[dict[str, str]]) -> None:
    if isinstance(value, float) and not math.isfinite(value):
        errors.append({"code": "NONFINITE_NUMBER", "path": path})
    elif isinstance(value, dict):
        for key, child in value.items():
            _scan_nonfinite(child, f"{path}.{key}", errors)
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _scan_nonfinite(child, f"{path}[{index}]", errors)


def _normalized_key(key: str) -> str:
    return "".join(ch for ch in key.lower() if ch.isalnum())


def _scan_observation_leaks(value: Any, path: str, errors: list[dict[str, str]]) -> None:
    if isinstance(value, dict):
        for key, child in value.items():
            if _normalized_key(str(key)) in _OBSERVATION_FORBIDDEN_KEYS:
                errors.append({"code": "GT_LEAK_IN_OBSERVATION", "path": f"{path}.{key}"})
            _scan_observation_leaks(child, f"{path}.{key}", errors)
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _scan_observation_leaks(child, f"{path}[{index}]", errors)


def _finite_vectors(value: Any, path: str, errors: list[dict[str, str]], *, minimum_rows: int = 3) -> list[list[float]]:
    if not isinstance(value, list) or not minimum_rows <= len(value) <= MAX_ENTITY_ROWS:
        errors.append({"code": "VECTOR_ROWS_INVALID", "path": path})
        return []
    rows: list[list[float]] = []
    for row_index, row in enumerate(value):
        if not isinstance(row, list) or len(row) != 3:
            errors.append({"code": "VECTOR_SHAPE_INVALID", "path": f"{path}[{row_index}]"})
            continue
        if any(isinstance(component, bool) or not isinstance(component, (float, int)) or not math.isfinite(component) for component in row):
            errors.append({"code": "NONFINITE_OR_NONNUMERIC_VECTOR", "path": f"{path}[{row_index}]"})
            continue
        rows.append([float(component) for component in row])
    return rows


def _triangle_rows(value: Any, vertex_count: int, path: str, errors: list[dict[str, str]]) -> list[tuple[int, int, int]]:
    if not isinstance(value, list) or not 1 <= len(value) <= MAX_SURFACE_TRIANGLES:
        errors.append({"code": "SURFACE_TRIANGLE_ROWS_INVALID", "path": path})
        return []
    result: list[tuple[int, int, int]] = []
    seen: set[tuple[int, int, int]] = set()
    for index, row in enumerate(value):
        row_path = f"{path}[{index}]"
        if not isinstance(row, list) or len(row) != 3 or any(isinstance(v, bool) or not isinstance(v, int) for v in row):
            errors.append({"code": "SURFACE_TRIANGLE_SHAPE_INVALID", "path": row_path})
            continue
        a, b, c = row
        if min(a, b, c) < 0 or max(a, b, c) >= vertex_count or len({a, b, c}) != 3:
            errors.append({"code": "SURFACE_TRIANGLE_INDEX_INVALID", "path": row_path})
            continue
        key = tuple(sorted((a, b, c)))
        if key in seen:
            errors.append({"code": "DUPLICATE_SURFACE_TRIANGLE", "path": row_path})
        seen.add(key)
        result.append((a, b, c))
    return result


def _triangle_area(vertices: list[list[float]], triangle: tuple[int, int, int]) -> float:
    a, b, c = (vertices[index] for index in triangle)
    u = [b[axis] - a[axis] for axis in range(3)]
    v = [c[axis] - a[axis] for axis in range(3)]
    cross = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]
    return 0.5 * math.sqrt(sum(component * component for component in cross))


def _canonical_oriented_triangle(triangle: tuple[int, int, int]) -> tuple[int, int, int]:
    rotations = (triangle, (triangle[1], triangle[2], triangle[0]), (triangle[2], triangle[0], triangle[1]))
    return min(rotations)


def _validate_artifact(
    document: Any,
    physical_receipt_verifier: Callable[[Mapping[str, Any]], bool] | None,
) -> ValidationReport:
    errors: list[dict[str, str]] = []
    physical_verified = False
    class_name: str | None = None
    if not _strict_object_keys(document, {
        "schema", "artifact_id", "qualification", "mother", "material", "load",
        "solver", "discretization", "fragmentation_basis", "fragments",
        "fracture_surface_pairs", "events", "observation", "clean_ground_truth",
    }, set(), "$", errors):
        return ValidationReport(False, tuple(errors), None, False, False, False)

    _scan_nonfinite(document, "$", errors)
    if document.get("schema") != SCHEMA:
        errors.append({"code": "SCHEMA_VERSION_UNSUPPORTED", "path": "$.schema"})
    _identifier(document.get("artifact_id"), "$.artifact_id", errors)

    qualification = document.get("qualification")
    if _strict_object_keys(qualification, {"class", "physical_validation_receipt"}, set(), "$.qualification", errors):
        class_name = qualification.get("class")
        if class_name not in QUALIFICATION_CLASSES:
            errors.append({"code": "INVALID_QUALIFICATION_CLASS", "path": "$.qualification.class"})
        receipt = qualification.get("physical_validation_receipt")
        if class_name == "PHYSICALLY_VALIDATED":
            if receipt is None:
                errors.append({"code": "PHYSICAL_RECEIPT_UNVERIFIED", "path": "$.qualification.physical_validation_receipt"})
            if not _strict_object_keys(receipt, {
                "receipt_id", "evidence_manifest_sha256", "criteria_sha256", "verifier_id", "signature_sha256",
            }, set(), "$.qualification.physical_validation_receipt", errors):
                receipt = None
            else:
                for field in ("receipt_id", "verifier_id"):
                    _identifier(receipt.get(field), f"$.qualification.physical_validation_receipt.{field}", errors)
                for field in ("evidence_manifest_sha256", "criteria_sha256", "signature_sha256"):
                    if not _sha256(receipt.get(field)):
                        errors.append({"code": "INVALID_SHA256", "path": f"$.qualification.physical_validation_receipt.{field}"})
                if physical_receipt_verifier is None:
                    errors.append({"code": "PHYSICAL_RECEIPT_UNVERIFIED", "path": "$.qualification.physical_validation_receipt"})
                else:
                    try:
                        physical_verified = physical_receipt_verifier(receipt) is True
                    except Exception:
                        physical_verified = False
                    if not physical_verified:
                        errors.append({"code": "PHYSICAL_RECEIPT_REJECTED", "path": "$.qualification.physical_validation_receipt"})
        elif receipt is not None:
            errors.append({"code": "UNEXPECTED_PHYSICAL_RECEIPT", "path": "$.qualification.physical_validation_receipt"})

    mother = document.get("mother")
    if _strict_object_keys(mother, {"mother_id", "mesh_path", "mesh_sha256", "reference_frame"}, set(), "$.mother", errors):
        _identifier(mother.get("mother_id"), "$.mother.mother_id", errors)
        _relative_payload_path(mother.get("mesh_path"), "$.mother.mesh_path", errors)
        if not _sha256(mother.get("mesh_sha256")):
            errors.append({"code": "INVALID_SHA256", "path": "$.mother.mesh_sha256"})
        frame = mother.get("reference_frame")
        if _strict_object_keys(frame, {"name", "length_unit", "mass_unit", "volume_unit", "time_unit"}, set(), "$.mother.reference_frame", errors):
            expected_frame = {"name": "right_handed_z_up", "length_unit": "m", "mass_unit": "kg", "volume_unit": "m3", "time_unit": "s"}
            for key, expected in expected_frame.items():
                if frame.get(key) != expected:
                    errors.append({"code": "REFERENCE_FRAME_OR_UNIT_MISMATCH", "path": f"$.mother.reference_frame.{key}"})

    material_density: float | None = None
    material = document.get("material")
    if _strict_object_keys(material, {"qualification", "parameters"}, set(), "$.material", errors):
        material_qualification = material.get("qualification")
        if material_qualification not in {"UNCALIBRATED", "CALIBRATED"}:
            errors.append({"code": "INVALID_MATERIAL_QUALIFICATION", "path": "$.material.qualification"})
        params = material.get("parameters")
        if _strict_object_keys(params, {
            "density_kg_m3", "young_modulus_pa", "poisson_ratio", "fracture_toughness_j_m2",
        }, set(), "$.material.parameters", errors):
            parameter_specs = {
                "density_kg_m3": ("kg/m3", True, None, None),
                "young_modulus_pa": ("Pa", True, None, None),
                "poisson_ratio": ("1", False, -1.0, 0.5),
                "fracture_toughness_j_m2": ("J/m2", True, None, None),
            }
            parameter_states: dict[str, str | None] = {}
            for key, (unit, positive, minimum, maximum) in parameter_specs.items():
                state, _number = _known_value(params.get(key), f"$.material.parameters.{key}", unit, errors,
                                               positive=positive, minimum=minimum, maximum=maximum)
                parameter_states[key] = state
                if key == "density_kg_m3" and state == "KNOWN":
                    material_density = _number
            if class_name in {"FIXTURE_GEOMETRY_ONLY", "SIMULATION_UNCALIBRATED"} and material_qualification != "UNCALIBRATED":
                errors.append({"code": "UNCALIBRATED_CLASS_MATERIAL_MISMATCH", "path": "$.material.qualification"})
            if class_name == "SIMULATION_UNCALIBRATED" and any(state != "KNOWN" for state in parameter_states.values()):
                errors.append({"code": "SIMULATION_PARAMETER_UNKNOWN", "path": "$.material.parameters"})
            if class_name == "PHYSICALLY_VALIDATED" and material_qualification != "CALIBRATED":
                errors.append({"code": "PHYSICAL_CLASS_REQUIRES_CALIBRATED_MATERIAL", "path": "$.material.qualification"})
            if class_name == "PHYSICALLY_VALIDATED" and any(state != "KNOWN" for state in parameter_states.values()):
                errors.append({"code": "PHYSICAL_CLASS_PARAMETER_UNKNOWN", "path": "$.material.parameters"})

    load = document.get("load")
    if _strict_object_keys(load, {"kind", "source_sha256", "boundary_conditions_sha256", "peak_force_n", "time_window_s"}, set(), "$.load", errors):
        _identifier(load.get("kind"), "$.load.kind", errors)
        for field in ("source_sha256", "boundary_conditions_sha256"):
            if not _sha256(load.get(field)):
                errors.append({"code": "INVALID_SHA256", "path": f"$.load.{field}"})
        _known_value(load.get("peak_force_n"), "$.load.peak_force_n", "N", errors, minimum=0.0)
        _known_value(load.get("time_window_s"), "$.load.time_window_s", "s", errors, minimum=0.0)
        if class_name == "FIXTURE_GEOMETRY_ONLY" and load.get("kind") != "FIXTURE_NO_PHYSICAL_LOAD":
            errors.append({"code": "FIXTURE_LOAD_KIND_INVALID", "path": "$.load.kind"})

    solver = document.get("solver")
    if _strict_object_keys(solver, {
        "execution_status", "solver_id", "version", "source_sha256", "build_sha256",
        "job_sha256", "result_sha256", "run_id",
    }, set(), "$.solver", errors):
        status = solver.get("execution_status")
        if status not in {"NOT_RUN", "COMPLETED"}:
            errors.append({"code": "INVALID_SOLVER_STATUS", "path": "$.solver.execution_status"})
        for field in ("solver_id", "version", "run_id"):
            _identifier(solver.get(field), f"$.solver.{field}", errors)
        for field in ("source_sha256", "build_sha256", "job_sha256", "result_sha256"):
            state = _known_sha(solver.get(field), f"$.solver.{field}", errors)
            if status == "COMPLETED" and state != "KNOWN":
                errors.append({"code": "COMPLETED_SOLVER_PROVENANCE_UNKNOWN", "path": f"$.solver.{field}"})
        if class_name in {"SIMULATION_UNCALIBRATED", "PHYSICALLY_VALIDATED"} and status != "COMPLETED":
            errors.append({"code": "QUALIFIED_CLASS_REQUIRES_COMPLETED_SOLVER", "path": "$.solver.execution_status"})
        if class_name == "FIXTURE_GEOMETRY_ONLY" and status != "NOT_RUN":
            errors.append({"code": "FIXTURE_MUST_NOT_CLAIM_SOLVER_RUN", "path": "$.solver.execution_status"})

    discretization = document.get("discretization")
    if _strict_object_keys(discretization, {"method", "recipe_sha256", "point_roster_sha256", "spacing_m", "horizon_m"}, set(), "$.discretization", errors):
        _identifier(discretization.get("method"), "$.discretization.method", errors)
        for field in ("recipe_sha256", "point_roster_sha256"):
            if not _sha256(discretization.get(field)):
                errors.append({"code": "INVALID_SHA256", "path": f"$.discretization.{field}"})
        _known_value(discretization.get("spacing_m"), "$.discretization.spacing_m", "m", errors, positive=True)
        _known_value(discretization.get("horizon_m"), "$.discretization.horizon_m", "m", errors, positive=True)

    basis = document.get("fragmentation_basis")
    if basis not in {"CLOSED_MESH_PARTITION", "SOLVER_STATE_AND_CLOSED_MESH"}:
        errors.append({"code": "GRAPH_COMPONENTS_CANNOT_DEFINE_FRAGMENTS", "path": "$.fragmentation_basis"})

    fragment_rows = document.get("fragments")
    fragment_ids: set[str] = set()
    fragment_by_id: dict[str, Mapping[str, Any]] = {}
    fragment_measurements: dict[str, tuple[float | None, float | None]] = {}
    if not isinstance(fragment_rows, list) or not 1 <= len(fragment_rows) <= MAX_ENTITY_ROWS:
        errors.append({"code": "FRAGMENT_ROWS_INVALID", "path": "$.fragments"})
        fragment_rows = []
    for index, fragment in enumerate(fragment_rows):
        path = f"$.fragments[{index}]"
        if not _strict_object_keys(fragment, {"fragment_id", "closed_mesh", "mass_kg", "volume_m3", "lineage"}, set(), path, errors):
            continue
        fragment_id = fragment.get("fragment_id")
        _identifier(fragment_id, f"{path}.fragment_id", errors)
        if isinstance(fragment_id, str):
            if fragment_id in fragment_ids:
                errors.append({"code": "DUPLICATE_FRAGMENT_ID", "path": f"{path}.fragment_id"})
            fragment_ids.add(fragment_id)
            fragment_by_id[fragment_id] = fragment
        mesh = fragment.get("closed_mesh")
        if _strict_object_keys(mesh, {"path", "sha256", "format", "frame", "length_unit", "watertight", "closure_method"}, set(), f"{path}.closed_mesh", errors):
            _relative_payload_path(mesh.get("path"), f"{path}.closed_mesh.path", errors)
            if not _sha256(mesh.get("sha256")):
                errors.append({"code": "INVALID_SHA256", "path": f"{path}.closed_mesh.sha256"})
            if mesh.get("format") not in {"PLY", "OBJ", "VTK", "GLB", "MESH_REFERENCE"}:
                errors.append({"code": "UNSUPPORTED_FRAGMENT_MESH_FORMAT", "path": f"{path}.closed_mesh.format"})
            if mesh.get("frame") != "mother_reference" or mesh.get("length_unit") != "m":
                errors.append({"code": "FRAGMENT_MESH_FRAME_OR_UNIT_MISMATCH", "path": f"{path}.closed_mesh"})
            if mesh.get("watertight") is not True:
                errors.append({"code": "FRAGMENT_MESH_NOT_PROVEN_CLOSED", "path": f"{path}.closed_mesh.watertight"})
            _identifier(mesh.get("closure_method"), f"{path}.closed_mesh.closure_method", errors)
        mass_state, mass_value = _known_value(fragment.get("mass_kg"), f"{path}.mass_kg", "kg", errors, positive=True)
        volume_state, volume_value = _known_value(fragment.get("volume_m3"), f"{path}.volume_m3", "m3", errors, positive=True)
        if isinstance(fragment_id, str):
            fragment_measurements[fragment_id] = (mass_value if mass_state == "KNOWN" else None, volume_value if volume_state == "KNOWN" else None)
        lineage = fragment.get("lineage")
        if _strict_object_keys(lineage, {"state", "parent_fragment_ids", "creation_event_id"}, {"reason"}, f"{path}.lineage", errors):
            if lineage.get("state") == "KNOWN":
                if "reason" in lineage or not isinstance(lineage.get("parent_fragment_ids"), list):
                    errors.append({"code": "LINEAGE_SHAPE_INVALID", "path": f"{path}.lineage"})
                else:
                    for pindex, parent_id in enumerate(lineage["parent_fragment_ids"]):
                        _identifier(parent_id, f"{path}.lineage.parent_fragment_ids[{pindex}]", errors)
                if lineage.get("creation_event_id") is not None:
                    _identifier(lineage.get("creation_event_id"), f"{path}.lineage.creation_event_id", errors)
            elif lineage.get("state") == "UNKNOWN":
                if "parent_fragment_ids" in lineage or "creation_event_id" in lineage or not isinstance(lineage.get("reason"), str) or not lineage["reason"].strip():
                    errors.append({"code": "EXPLICIT_UNKNOWN_REASON_REQUIRED", "path": f"{path}.lineage"})
            else:
                errors.append({"code": "INVALID_KNOWLEDGE_STATE", "path": f"{path}.lineage.state"})
        if class_name == "PHYSICALLY_VALIDATED" and (mass_state != "KNOWN" or volume_state != "KNOWN"):
            errors.append({"code": "PHYSICAL_CLASS_REQUIRES_FRAGMENT_MASS_VOLUME", "path": path})
        if material_density is not None and mass_value is not None and volume_value is not None:
            if not math.isclose(mass_value, material_density * volume_value, rel_tol=1e-4, abs_tol=1e-12):
                errors.append({"code": "FRAGMENT_MASS_VOLUME_DENSITY_MISMATCH", "path": path})

    event_rows = document.get("events")
    events: dict[str, Mapping[str, Any]] = {}
    if not isinstance(event_rows, list) or not 1 <= len(event_rows) <= MAX_ENTITY_ROWS:
        errors.append({"code": "EVENT_ROWS_INVALID", "path": "$.events"})
        event_rows = []
    previous_step = -1
    for index, event in enumerate(event_rows):
        path = f"$.events[{index}]"
        if not _strict_object_keys(event, {"event_id", "time_s", "solver_step", "parent_event_id", "fracture_pair_ids", "kind"}, set(), path, errors):
            continue
        event_id = event.get("event_id")
        _identifier(event_id, f"{path}.event_id", errors)
        if isinstance(event_id, str) and event_id in events:
            errors.append({"code": "DUPLICATE_EVENT_ID", "path": f"{path}.event_id"})
        elif isinstance(event_id, str):
            events[event_id] = event
        time_state, event_time = _known_value(event.get("time_s"), f"{path}.time_s", "s", errors, minimum=0.0)
        step_state, step = _known_value(event.get("solver_step"), f"{path}.solver_step", "1", errors, minimum=0.0)
        if step_state == "KNOWN" and step is not None:
            if not step.is_integer():
                errors.append({"code": "SOLVER_STEP_NOT_INTEGER", "path": f"{path}.solver_step.value"})
            if step < previous_step:
                errors.append({"code": "EVENT_TIME_ORDER_INVALID", "path": path})
            previous_step = int(step)
        if class_name in {"SIMULATION_UNCALIBRATED", "PHYSICALLY_VALIDATED"} and (time_state != "KNOWN" or step_state != "KNOWN"):
            errors.append({"code": "SIMULATED_EVENT_TIME_UNKNOWN", "path": path})
        _identifier(event.get("kind"), f"{path}.kind", errors)
        if event.get("parent_event_id") is not None:
            _identifier(event.get("parent_event_id"), f"{path}.parent_event_id", errors)
        pair_ids = event.get("fracture_pair_ids")
        if not isinstance(pair_ids, list) or len(pair_ids) > MAX_ENTITY_ROWS:
            errors.append({"code": "EVENT_PAIR_IDS_INVALID", "path": f"{path}.fracture_pair_ids"})
        else:
            for pair_index, pair_id in enumerate(pair_ids):
                _identifier(pair_id, f"{path}.fracture_pair_ids[{pair_index}]", errors)

    pair_rows = document.get("fracture_surface_pairs")
    pair_ids: set[str] = set()
    surface_ids: set[str] = set()
    pairs_by_id: dict[str, Mapping[str, Any]] = {}
    if not isinstance(pair_rows, list) or not 1 <= len(pair_rows) <= MAX_ENTITY_ROWS:
        errors.append({"code": "FRACTURE_PAIR_ROWS_INVALID", "path": "$.fracture_surface_pairs"})
        pair_rows = []
    for index, pair in enumerate(pair_rows):
        path = f"$.fracture_surface_pairs[{index}]"
        if not _strict_object_keys(pair, {"pair_id", "event_id", "negative_fragment_id", "positive_fragment_id", "reference_area_m2", "negative_surface", "positive_surface"}, set(), path, errors):
            continue
        pair_id = pair.get("pair_id")
        _identifier(pair_id, f"{path}.pair_id", errors)
        if isinstance(pair_id, str):
            if pair_id in pair_ids:
                errors.append({"code": "DUPLICATE_PAIR_ID", "path": f"{path}.pair_id"})
            pair_ids.add(pair_id)
            pairs_by_id[pair_id] = pair
        event_id = pair.get("event_id")
        _identifier(event_id, f"{path}.event_id", errors)
        if event_id not in events:
            errors.append({"code": "PAIR_EVENT_NOT_FOUND", "path": f"{path}.event_id"})
        negative_fragment = pair.get("negative_fragment_id")
        positive_fragment = pair.get("positive_fragment_id")
        if negative_fragment == positive_fragment:
            errors.append({"code": "PAIR_FRAGMENTS_NOT_DISTINCT", "path": path})
        side_coordinates: dict[str, list[list[float]]] = {}
        side_normals: dict[str, list[list[float]]] = {}
        side_triangles: dict[str, list[tuple[int, int, int]]] = {}
        for side_name, fragment_id in (("negative_surface", negative_fragment), ("positive_surface", positive_fragment)):
            side = pair.get(side_name)
            side_path = f"{path}.{side_name}"
            if not _strict_object_keys(side, {"surface_id", "reciprocal_surface_id", "fragment_id", "mother_coordinates", "normals", "triangles", "parent_surface_id"}, set(), side_path, errors):
                continue
            for field in ("surface_id", "reciprocal_surface_id"):
                _identifier(side.get(field), f"{side_path}.{field}", errors)
            surface_id = side.get("surface_id")
            if isinstance(surface_id, str):
                if surface_id in surface_ids:
                    errors.append({"code": "DUPLICATE_SURFACE_ID", "path": f"{side_path}.surface_id"})
                surface_ids.add(surface_id)
            if side.get("fragment_id") != fragment_id or fragment_id not in fragment_by_id:
                errors.append({"code": "SURFACE_FRAGMENT_PAIR_MISMATCH", "path": f"{side_path}.fragment_id"})
            reciprocal = pair.get("positive_surface" if side_name == "negative_surface" else "negative_surface")
            if isinstance(reciprocal, Mapping) and side.get("reciprocal_surface_id") != reciprocal.get("surface_id"):
                errors.append({"code": "RECIPROCAL_SURFACE_ID_MISMATCH", "path": f"{side_path}.reciprocal_surface_id"})
            coords_obj = side.get("mother_coordinates")
            coords: list[list[float]] = []
            if _strict_object_keys(coords_obj, {"frame", "length_unit", "vertices"}, set(), f"{side_path}.mother_coordinates", errors):
                if coords_obj.get("frame") != "mother_reference" or coords_obj.get("length_unit") != "m":
                    errors.append({"code": "UNIT_MISMATCH", "path": f"{side_path}.mother_coordinates"})
                coords = _finite_vectors(coords_obj.get("vertices"), f"{side_path}.mother_coordinates.vertices", errors)
            normals_obj = side.get("normals")
            normals: list[list[float]] = []
            if _strict_object_keys(normals_obj, {"frame", "unit", "vectors"}, set(), f"{side_path}.normals", errors):
                if normals_obj.get("frame") != "mother_reference" or normals_obj.get("unit") != "1":
                    errors.append({"code": "UNIT_MISMATCH", "path": f"{side_path}.normals"})
                normals = _finite_vectors(normals_obj.get("vectors"), f"{side_path}.normals.vectors", errors)
                if len(normals) != len(coords):
                    errors.append({"code": "NORMAL_COORDINATE_COUNT_MISMATCH", "path": f"{side_path}.normals"})
                for nindex, normal in enumerate(normals):
                    norm = math.sqrt(sum(component * component for component in normal))
                    if not math.isclose(norm, 1.0, rel_tol=1e-5, abs_tol=1e-5):
                        errors.append({"code": "NORMAL_NOT_UNIT_LENGTH", "path": f"{side_path}.normals.vectors[{nindex}]"})
            if side.get("parent_surface_id") is not None:
                _identifier(side.get("parent_surface_id"), f"{side_path}.parent_surface_id", errors)
            side_coordinates[side_name] = coords
            side_normals[side_name] = normals
            side_triangles[side_name] = _triangle_rows(side.get("triangles"), len(coords), f"{side_path}.triangles", errors)
        negative = pair.get("negative_surface")
        positive = pair.get("positive_surface")
        if isinstance(negative, Mapping) and isinstance(positive, Mapping):
            neg_coords, pos_coords = side_coordinates.get("negative_surface", []), side_coordinates.get("positive_surface", [])
            if neg_coords and pos_coords and (len(neg_coords) != len(pos_coords) or any(
                any(not math.isclose(a, b, rel_tol=0.0, abs_tol=1e-12) for a, b in zip(na, pa))
                for na, pa in zip(neg_coords, pos_coords)
            )):
                errors.append({"code": "PAIRED_MOTHER_COORDINATES_DIFFER", "path": path})
            neg_normals, pos_normals = side_normals.get("negative_surface", []), side_normals.get("positive_surface", [])
            if len(neg_normals) == len(pos_normals):
                for nindex, (nn, pn) in enumerate(zip(neg_normals, pos_normals)):
                    dot = sum(a * b for a, b in zip(nn, pn))
                    if not math.isclose(dot, -1.0, rel_tol=1e-4, abs_tol=1e-4):
                        errors.append({"code": "PAIRED_NORMALS_NOT_OPPOSITE", "path": f"{path}.normals[{nindex}]"})
            neg_triangles = side_triangles.get("negative_surface", [])
            pos_triangles = side_triangles.get("positive_surface", [])
            neg_oriented = sorted(_canonical_oriented_triangle(t) for t in neg_triangles)
            pos_oriented_reversed = sorted(_canonical_oriented_triangle((t[0], t[2], t[1])) for t in pos_triangles)
            if neg_oriented != pos_oriented_reversed:
                errors.append({"code": "PAIRED_SURFACE_TRIANGLES_NOT_RECIPROCAL", "path": path})
            neg_area = sum(_triangle_area(neg_coords, t) for t in neg_triangles) if neg_coords else 0.0
            pos_area = sum(_triangle_area(pos_coords, t) for t in pos_triangles) if pos_coords else 0.0
            if neg_area <= 0 or not math.isclose(neg_area, pos_area, rel_tol=1e-6, abs_tol=1e-15):
                errors.append({"code": "PAIRED_SURFACE_AREA_INVALID", "path": path})
        declared_area_state, declared_area = _known_value(pair.get("reference_area_m2"), f"{path}.reference_area_m2", "m2", errors, positive=True)
        if declared_area_state == "KNOWN" and declared_area is not None and isinstance(negative, Mapping):
            neg_coords = side_coordinates.get("negative_surface", [])
            neg_triangles = side_triangles.get("negative_surface", [])
            measured_area = sum(_triangle_area(neg_coords, t) for t in neg_triangles) if neg_coords else 0.0
            if not math.isclose(measured_area, declared_area, rel_tol=1e-4, abs_tol=1e-15):
                errors.append({"code": "DECLARED_SURFACE_AREA_MISMATCH", "path": f"{path}.reference_area_m2"})

    for event_id, event in events.items():
        parent_event_id = event.get("parent_event_id")
        if parent_event_id is not None and parent_event_id not in events:
            errors.append({"code": "PARENT_EVENT_NOT_FOUND", "path": f"$.events[{event_id}].parent_event_id"})
        elif parent_event_id is not None:
            parent_time = events[parent_event_id].get("time_s", {})
            current_time = event.get("time_s", {})
            if isinstance(parent_time, Mapping) and isinstance(current_time, Mapping) and parent_time.get("state") == current_time.get("state") == "KNOWN":
                if isinstance(parent_time.get("value"), (float, int)) and isinstance(current_time.get("value"), (float, int)) and current_time["value"] < parent_time["value"]:
                    errors.append({"code": "EVENT_TIME_ORDER_INVALID", "path": f"$.events[{event_id}].time_s"})
        declared = event.get("fracture_pair_ids")
        if isinstance(declared, list):
            for pair_id in declared:
                if pair_id not in pairs_by_id or pairs_by_id.get(pair_id, {}).get("event_id") != event_id:
                    errors.append({"code": "EVENT_PAIR_LINEAGE_MISMATCH", "path": f"$.events[{event_id}].fracture_pair_ids"})
    for pair_id, pair in pairs_by_id.items():
        event = events.get(pair.get("event_id"))
        if not isinstance(event, Mapping) or pair_id not in event.get("fracture_pair_ids", []):
            errors.append({"code": "PAIR_EVENT_LINEAGE_MISMATCH", "path": f"$.fracture_surface_pairs[{pair_id}].event_id"})

    for fragment_id, fragment in fragment_by_id.items():
        lineage = fragment.get("lineage")
        if not isinstance(lineage, Mapping):
            continue
        if lineage.get("state") == "KNOWN":
            for parent_id in lineage.get("parent_fragment_ids", []):
                if parent_id not in fragment_by_id:
                    errors.append({"code": "FRAGMENT_PARENT_NOT_FOUND", "path": f"$.fragments[{fragment_id}].lineage.parent_fragment_ids"})
            creation_event = lineage.get("creation_event_id")
            if creation_event is not None and creation_event not in events:
                errors.append({"code": "FRAGMENT_CREATION_EVENT_NOT_FOUND", "path": f"$.fragments[{fragment_id}].lineage.creation_event_id"})

    observation = document.get("observation")
    _scan_observation_leaks(observation, "$.observation", errors)
    if _strict_object_keys(observation, {"capture_id", "frame", "length_unit", "processing_recipe_sha256", "fragment_observations"}, set(), "$.observation", errors):
        _identifier(observation.get("capture_id"), "$.observation.capture_id", errors)
        if observation.get("frame") != "sensor_observation" or observation.get("length_unit") != "m":
            errors.append({"code": "OBSERVATION_FRAME_OR_UNIT_INVALID", "path": "$.observation"})
        if not _sha256(observation.get("processing_recipe_sha256")):
            errors.append({"code": "INVALID_SHA256", "path": "$.observation.processing_recipe_sha256"})
        observed_rows = observation.get("fragment_observations")
        if not isinstance(observed_rows, list) or not 1 <= len(observed_rows) <= MAX_ENTITY_ROWS:
            errors.append({"code": "OBSERVATION_FRAGMENT_ROWS_INVALID", "path": "$.observation.fragment_observations"})
        else:
            observed_ids: set[str] = set()
            for index, observed in enumerate(observed_rows):
                path = f"$.observation.fragment_observations[{index}]"
                if not _strict_object_keys(observed, {"observation_fragment_id", "observed_mesh", "visibility_mask", "sample_count", "visible_area_m2"}, set(), path, errors):
                    continue
                oid = observed.get("observation_fragment_id")
                if not isinstance(oid, str) or OBSERVATION_ID_RE.fullmatch(oid) is None:
                    errors.append({"code": "OBSERVATION_ID_NOT_OPAQUE", "path": f"{path}.observation_fragment_id"})
                elif oid in observed_ids:
                    errors.append({"code": "DUPLICATE_OBSERVATION_ID", "path": f"{path}.observation_fragment_id"})
                observed_ids.add(oid) if isinstance(oid, str) else None
                for field in ("observed_mesh", "visibility_mask"):
                    payload = observed.get(field)
                    payload_path = f"{path}.{field}"
                    if _strict_object_keys(payload, {"path", "sha256", "format"}, set(), payload_path, errors):
                        _relative_payload_path(payload.get("path"), f"{payload_path}.path", errors)
                        if not _sha256(payload.get("sha256")):
                            errors.append({"code": "INVALID_SHA256", "path": f"{payload_path}.sha256"})
                        allowed_observation_formats = {"PLY"} if field == "observed_mesh" else {"JSON"}
                        if payload.get("format") not in allowed_observation_formats:
                            errors.append({"code": "UNSUPPORTED_OBSERVATION_PAYLOAD_FORMAT", "path": f"{payload_path}.format"})
                _known_value(observed.get("sample_count"), f"{path}.sample_count", "1", errors, minimum=0.0)
                _known_value(observed.get("visible_area_m2"), f"{path}.visible_area_m2", "m2", errors, minimum=0.0)

    truth = document.get("clean_ground_truth")
    if _strict_object_keys(truth, {"namespace", "labels"}, set(), "$.clean_ground_truth", errors):
        if truth.get("namespace") != "pgfrag.clean_ground_truth.v1":
            errors.append({"code": "GROUND_TRUTH_NAMESPACE_INVALID", "path": "$.clean_ground_truth.namespace"})
        labels = truth.get("labels")
        if not isinstance(labels, list) or not 1 <= len(labels) <= MAX_ENTITY_ROWS:
            errors.append({"code": "MISSING_OR_INVALID_LABELS", "path": "$.clean_ground_truth.labels"})
        else:
            label_ids: set[str] = set()
            for index, label in enumerate(labels):
                path = f"$.clean_ground_truth.labels[{index}]"
                if not _strict_object_keys(label, {"pair_id", "relation_kind", "negative_fragment_id", "positive_fragment_id", "event_id", "label_provenance_sha256"}, set(), path, errors):
                    continue
                pid = label.get("pair_id")
                if pid in label_ids:
                    errors.append({"code": "DUPLICATE_TRUTH_LABEL", "path": f"{path}.pair_id"})
                label_ids.add(pid) if isinstance(pid, str) else None
                pair = pairs_by_id.get(pid)
                if pair is None or any(label.get(key) != pair.get(key) for key in ("negative_fragment_id", "positive_fragment_id", "event_id")):
                    errors.append({"code": "TRUTH_LABEL_PAIR_MISMATCH", "path": path})
                relation = label.get("relation_kind")
                if class_name == "FIXTURE_GEOMETRY_ONLY" and relation != "FIXTURE_PAIR":
                    errors.append({"code": "FIXTURE_TRUTH_CANNOT_CLAIM_PHYSICAL_RELATION", "path": f"{path}.relation_kind"})
                if class_name == "SIMULATION_UNCALIBRATED" and relation not in {"SIMULATED_CANDIDATE", "UNKNOWN"}:
                    errors.append({"code": "UNCALIBRATED_SIMULATION_RELATION_INVALID", "path": f"{path}.relation_kind"})
                if class_name == "PHYSICALLY_VALIDATED" and relation not in {"R2", "UNKNOWN"}:
                    errors.append({"code": "PHYSICAL_TRUTH_RELATION_INVALID", "path": f"{path}.relation_kind"})
                if not _sha256(label.get("label_provenance_sha256")):
                    errors.append({"code": "INVALID_SHA256", "path": f"{path}.label_provenance_sha256"})
            if pair_ids.difference(label_ids):
                errors.append({"code": "PAIR_LABEL_MISSING", "path": "$.clean_ground_truth.labels"})

    valid = not errors
    pretraining_qualified = valid and class_name == "PHYSICALLY_VALIDATED" and physical_verified
    return ValidationReport(valid, tuple(errors), class_name, physical_verified, pretraining_qualified, valid)


def validate_artifact(
    document: Any,
    *,
    physical_receipt_verifier: Callable[[Mapping[str, Any]], bool] | None = None,
) -> ValidationReport:
    """Validate an in-memory artifact; qualification never comes from a checkbox."""
    try:
        return _validate_artifact(document, physical_receipt_verifier)
    except Exception:
        return ValidationReport(False, ({"code": "SAFE_VALIDATION_FAILURE", "path": "$"},), None, False, False, False)


def _reject_constant(token: str) -> None:
    raise ValueError(f"non-finite JSON constant: {token}")


def _relative_payload_path(value: Any, path: str, errors: list[dict[str, str]]) -> bool:
    if not isinstance(value, str) or not value or value.strip() != value or "\\" in value or ":" in value:
        errors.append({"code": "PAYLOAD_PATH_INVALID", "path": path})
        return False
    parts = value.split("/")
    if value.startswith("/") or any(part in {"", ".", ".."} for part in parts):
        errors.append({"code": "PAYLOAD_PATH_INVALID", "path": path})
        return False
    return True


def _is_linklike(path: Path) -> bool:
    try:
        if path.is_symlink():
            return True
        is_junction = getattr(path, "is_junction", None)
        return bool(is_junction is not None and is_junction())
    except OSError:
        return True


def _verified_payload(root: Path, reference: Mapping[str, Any], *, include_content: bool = False) -> dict[str, Any]:
    raw_path = reference["path"]
    parts = raw_path.split("/")
    candidate = root
    for part in parts:
        candidate = candidate / part
        if _is_linklike(candidate):
            raise ArtifactContractError(code="PAYLOAD_LINK_REJECTED")
    try:
        resolved = candidate.resolve(strict=True)
        resolved.relative_to(root)
        file_stat = resolved.stat()
        if not resolved.is_file() or file_stat.st_size > MAX_PAYLOAD_FILE_BYTES:
            raise ArtifactContractError(code="PAYLOAD_FILE_INVALID_OR_TOO_LARGE")
        with resolved.open("rb") as stream:
            opened_stat = os.fstat(stream.fileno())
            if not stat.S_ISREG(opened_stat.st_mode) or opened_stat.st_size > MAX_PAYLOAD_FILE_BYTES:
                raise ArtifactContractError(code="PAYLOAD_FILE_INVALID_OR_TOO_LARGE")
            raw = stream.read(MAX_PAYLOAD_FILE_BYTES + 1)
            if len(raw) > MAX_PAYLOAD_FILE_BYTES or len(raw) != opened_stat.st_size:
                raise ArtifactContractError(code="PAYLOAD_FILE_INVALID_OR_TOO_LARGE")
    except ArtifactContractError:
        raise
    except (OSError, ValueError) as exc:
        raise ArtifactContractError(code="PAYLOAD_FILE_UNAVAILABLE_OR_OUTSIDE_ROOT") from exc
    if hashlib.sha256(raw).hexdigest() != reference["sha256"]:
        raise ArtifactContractError(code="PAYLOAD_SHA256_MISMATCH")
    result = {"bytes_verified": True, "bytes": len(raw), "sha256_verified": True, "sha256": reference["sha256"]}
    if include_content:
        result["content_bytes"] = raw
    else:
        result["path"] = raw_path
    return result


def _verify_artifact_payloads(document: Mapping[str, Any], root: Path, *, include_observations: bool = True) -> list[dict[str, Any]]:
    references: list[Mapping[str, Any]] = [
        {"path": document["mother"]["mesh_path"], "sha256": document["mother"]["mesh_sha256"]},
        *(fragment["closed_mesh"] for fragment in document["fragments"]),
    ]
    if include_observations:
        for observation in document["observation"]["fragment_observations"]:
            references.extend((observation["observed_mesh"], observation["visibility_mask"]))
    results: list[dict[str, Any]] = []
    total = 0
    for reference in references:
        item = _verified_payload(root, reference)
        total += item["bytes"]
        if total > MAX_PAYLOAD_TOTAL_BYTES:
            raise ArtifactContractError(code="PAYLOAD_TOTAL_TOO_LARGE")
        results.append(item)
    return results


def _unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("duplicate JSON object key")
        result[key] = value
    return result


def _parse_observed_ply(raw: bytes) -> dict[str, int]:
    try:
        lines = raw.decode("ascii").splitlines()
    except UnicodeDecodeError as exc:
        raise ArtifactContractError(code="OBSERVATION_PLY_NOT_ASCII") from exc
    header_end = next((index for index, line in enumerate(lines[:16]) if line == "end_header"), None)
    if header_end is None or header_end != 8 or lines[:2] != ["ply", "format ascii 1.0"]:
        raise ArtifactContractError(code="OBSERVATION_PLY_HEADER_UNSUPPORTED")
    header = lines[:9]
    if header[2] != header[2].strip() or not header[2].startswith("element vertex "):
        raise ArtifactContractError(code="OBSERVATION_PLY_HEADER_UNSUPPORTED")
    if header[3:6] != ["property float x", "property float y", "property float z"]:
        raise ArtifactContractError(code="OBSERVATION_PLY_HEADER_UNSUPPORTED")
    if not header[6].startswith("element face ") or header[7] != "property list uchar int vertex_indices" or header[8] != "end_header":
        raise ArtifactContractError(code="OBSERVATION_PLY_HEADER_UNSUPPORTED")
    # The minimal fixture writer puts end_header after the face declaration and
    # property; keep the exact accepted header shape small and metadata-free.
    vertex_token = header[2].split(" ")
    face_token = header[6].split(" ")
    if len(vertex_token) != 3 or len(face_token) != 3:
        raise ArtifactContractError(code="OBSERVATION_PLY_HEADER_UNSUPPORTED")
    try:
        vertex_count, face_count = int(vertex_token[2]), int(face_token[2])
    except ValueError as exc:
        raise ArtifactContractError(code="OBSERVATION_PLY_HEADER_UNSUPPORTED") from exc
    if not 1 <= vertex_count <= MAX_ENTITY_ROWS or not 0 <= face_count <= MAX_SURFACE_TRIANGLES:
        raise ArtifactContractError(code="OBSERVATION_PLY_COUNT_INVALID")
    data_lines = lines[header_end + 1:]
    if len(data_lines) != vertex_count + face_count:
        raise ArtifactContractError(code="OBSERVATION_PLY_ROW_COUNT_MISMATCH")
    vertices: list[tuple[float, float, float]] = []
    for line in data_lines[:vertex_count]:
        pieces = line.split()
        if len(pieces) != 3:
            raise ArtifactContractError(code="OBSERVATION_PLY_VERTEX_INVALID")
        try:
            point = tuple(float(token) for token in pieces)
        except ValueError as exc:
            raise ArtifactContractError(code="OBSERVATION_PLY_VERTEX_INVALID") from exc
        if not all(math.isfinite(component) for component in point):
            raise ArtifactContractError(code="OBSERVATION_PLY_VERTEX_NONFINITE")
        vertices.append(point)
    for line in data_lines[vertex_count:]:
        pieces = line.split()
        if len(pieces) != 4 or pieces[0] != "3":
            raise ArtifactContractError(code="OBSERVATION_PLY_FACE_INVALID")
        try:
            indices = tuple(int(token) for token in pieces[1:])
        except ValueError as exc:
            raise ArtifactContractError(code="OBSERVATION_PLY_FACE_INVALID") from exc
        if len(set(indices)) != 3 or min(indices) < 0 or max(indices) >= vertex_count:
            raise ArtifactContractError(code="OBSERVATION_PLY_FACE_INVALID")
    return {"vertex_count": vertex_count, "face_count": face_count}


def _parse_visibility_mask(raw: bytes, vertex_count: int) -> list[int]:
    try:
        data = json.loads(raw.decode("utf-8"), parse_constant=_reject_constant, object_pairs_hook=_unique_object)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise ArtifactContractError(code="VISIBILITY_MASK_STRICT_JSON_REJECTED") from exc
    if not isinstance(data, dict) or set(data) != {"visible_vertex_indices"}:
        raise ArtifactContractError(code="VISIBILITY_MASK_SCHEMA_REJECTED")
    indices = data["visible_vertex_indices"]
    if not isinstance(indices, list) or len(indices) > MAX_ENTITY_ROWS:
        raise ArtifactContractError(code="VISIBILITY_MASK_ROWS_INVALID")
    if any(isinstance(index, bool) or not isinstance(index, int) or index < 0 or index >= vertex_count for index in indices):
        raise ArtifactContractError(code="VISIBILITY_MASK_INDEX_INVALID")
    if len(set(indices)) != len(indices):
        raise ArtifactContractError(code="VISIBILITY_MASK_DUPLICATE_INDEX")
    return indices


def _load_document(path: str | os.PathLike[str], expected_artifact_sha256: str | None = None) -> dict[str, Any]:
    source = Path(path)
    if _is_linklike(source):
        raise ArtifactContractError(code="ARTIFACT_LINK_REJECTED")
    try:
        with source.open("rb") as stream:
            file_stat = os.fstat(stream.fileno())
            if not stat.S_ISREG(file_stat.st_mode) or file_stat.st_size > MAX_DOCUMENT_BYTES:
                raise ArtifactContractError(code="DOCUMENT_TOO_LARGE")
            raw = stream.read(MAX_DOCUMENT_BYTES + 1)
            if len(raw) > MAX_DOCUMENT_BYTES or len(raw) != file_stat.st_size:
                raise ArtifactContractError(code="DOCUMENT_TOO_LARGE")
    except ArtifactContractError:
        raise
    except OSError as exc:
        raise ArtifactContractError(code="ARTIFACT_FILE_UNAVAILABLE") from exc
    if len(raw) > MAX_DOCUMENT_BYTES:
        raise ArtifactContractError(code="DOCUMENT_TOO_LARGE")
    if expected_artifact_sha256 is not None:
        if not _sha256(expected_artifact_sha256):
            raise ArtifactContractError(code="EXPECTED_ARTIFACT_SHA256_INVALID")
        if hashlib.sha256(raw).hexdigest() != expected_artifact_sha256:
            raise ArtifactContractError(code="ARTIFACT_SHA256_MISMATCH")
    try:
        document = json.loads(raw.decode("utf-8"), parse_constant=_reject_constant, object_pairs_hook=_unique_object)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise ArtifactContractError(code="STRICT_JSON_REJECTED") from exc
    if not isinstance(document, dict):
        raise ArtifactContractError(code="JSON_ROOT_OBJECT_REQUIRED")
    return document


def export_artifact(
    document: Mapping[str, Any],
    path: str | os.PathLike[str],
    *,
    physical_receipt_verifier: Callable[[Mapping[str, Any]], bool] | None = None,
) -> dict[str, Any]:
    """Validate and atomically write a strict UTF-8 JSON artifact."""
    report = validate_artifact(document, physical_receipt_verifier=physical_receipt_verifier)
    if not report.valid:
        raise ArtifactContractError(report)
    try:
        raw = (json.dumps(document, sort_keys=True, indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode("utf-8")
    except (TypeError, ValueError) as exc:
        raise ArtifactContractError(code="STRICT_JSON_SERIALIZATION_REJECTED") from exc
    if len(raw) > MAX_DOCUMENT_BYTES:
        raise ArtifactContractError(code="DOCUMENT_TOO_LARGE")
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    try:
        payload_root = target.parent.resolve(strict=True)
        _verify_artifact_payloads(document, payload_root)
    except KeyError as exc:
        raise ArtifactContractError(code="PAYLOAD_REFERENCE_MISSING") from exc
    temporary = target.with_name(target.name + ".tmp")
    try:
        temporary.write_bytes(raw)
        temporary.replace(target)
    finally:
        if temporary.exists():
            temporary.unlink()
    return {
        "path": str(target),
        "bytes": len(raw),
        "sha256": hashlib.sha256(raw).hexdigest(),
        "validation": report.as_dict(),
    }


def read_observation(
    path: str | os.PathLike[str],
    *,
    expected_artifact_sha256: str | None = None,
    physical_receipt_verifier: Callable[[Mapping[str, Any]], bool] | None = None,
) -> dict[str, Any]:
    """Return only the observation channel accepted as A12 inference input."""
    document = _load_document(path, expected_artifact_sha256)
    report = validate_artifact(document, physical_receipt_verifier=physical_receipt_verifier)
    if not report.valid:
        raise ArtifactContractError(report)
    source_root = Path(path).resolve(strict=True).parent
    observation = document["observation"]
    payloads: list[dict[str, Any]] = []
    projected_rows: list[dict[str, Any]] = []
    for observed in observation["fragment_observations"]:
        observation_fragment_id = observed["observation_fragment_id"]
        mesh_reference = observed["observed_mesh"]
        mask_reference = observed["visibility_mask"]
        mesh_payload = _verified_payload(source_root, mesh_reference, include_content=True)
        mesh_raw = mesh_payload.pop("content_bytes")
        mesh_counts = _parse_observed_ply(mesh_raw)
        mask_payload = _verified_payload(source_root, mask_reference, include_content=True)
        mask_raw = mask_payload.pop("content_bytes")
        mask_indices = _parse_visibility_mask(mask_raw, mesh_counts["vertex_count"])
        count_measure = observed["sample_count"]
        if count_measure.get("state") == "KNOWN" and float(count_measure["value"]) != len(mask_indices):
            raise ArtifactContractError(code="OBSERVATION_SAMPLE_COUNT_MISMATCH")
        mesh_payload["content_base64"] = base64.b64encode(mesh_raw).decode("ascii")
        mesh_payload.update({"role": "observed_mesh", "observation_fragment_id": observation_fragment_id, "format": mesh_reference["format"], **mesh_counts})
        mask_payload["content_base64"] = base64.b64encode(mask_raw).decode("ascii")
        mask_payload.update({"role": "visibility_mask", "observation_fragment_id": observation_fragment_id, "format": mask_reference["format"], "visible_vertex_count": len(mask_indices)})
        payloads.extend((mesh_payload, mask_payload))
        if sum(item["bytes"] for item in payloads) > MAX_PAYLOAD_TOTAL_BYTES:
            raise ArtifactContractError(code="PAYLOAD_TOTAL_TOO_LARGE")
        projected_rows.append({
            "observation_fragment_id": observation_fragment_id,
            "observed_mesh": {"sha256": observed["observed_mesh"]["sha256"], "format": observed["observed_mesh"]["format"]},
            "visibility_mask": {"sha256": observed["visibility_mask"]["sha256"], "format": observed["visibility_mask"]["format"]},
            "sample_count": observed["sample_count"],
            "visible_area_m2": observed["visible_area_m2"],
        })
    return {
        "schema": OBSERVATION_SCHEMA,
        "verified_payloads": payloads,
        "observation": {
            "capture_id": observation["capture_id"],
            "frame": observation["frame"],
            "length_unit": observation["length_unit"],
            "processing_recipe_sha256": observation["processing_recipe_sha256"],
            "fragment_observations": projected_rows,
        },
    }


def read_clean_ground_truth_for_audit(
    path: str | os.PathLike[str],
    *,
    purpose: str,
    expected_artifact_sha256: str | None = None,
    physical_receipt_verifier: Callable[[Mapping[str, Any]], bool] | None = None,
) -> dict[str, Any]:
    """Read clean GT on an explicit audit-only path, never the inference projection."""
    if purpose != "offline_audit":
        raise ArtifactContractError(code="GROUND_TRUTH_AUDIT_PURPOSE_REQUIRED")
    document = _load_document(path, expected_artifact_sha256)
    report = validate_artifact(document, physical_receipt_verifier=physical_receipt_verifier)
    if not report.valid:
        raise ArtifactContractError(report)
    source_root = Path(path).resolve(strict=True).parent
    _verify_artifact_payloads(document, source_root, include_observations=False)
    return {
        "role": "AUDIT_ONLY_CLEAN_GROUND_TRUTH",
        "qualification_class": report.qualification_class,
        "clean_ground_truth": json.loads(json.dumps(document["clean_ground_truth"], allow_nan=False)),
    }


__all__ = [
    "ArtifactContractError",
    "ValidationReport",
    "export_artifact",
    "read_clean_ground_truth_for_audit",
    "read_observation",
    "validate_artifact",
]
