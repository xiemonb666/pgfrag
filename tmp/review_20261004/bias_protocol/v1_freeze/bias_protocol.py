"""Fail-closed protocol checks for independent discretization/direction studies.

This module checks study design and evidence provenance only. It never creates a
mesh or point cloud and never invokes a mechanics solver.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from collections.abc import Mapping, Sequence
from typing import Any

SCHEMA = "pgfrag.bias-study.v1"
READINESS_SCHEMA = "pgfrag.bias-readiness.v1"
_SHA256 = re.compile(r"^[0-9a-fA-F]{64}$")
_METRIC_NAMES = (
    "energy",
    "response",
    "mass_distribution",
    "crack_direction",
    "fragment_distribution",
)
_PHYSICAL_FIELDS = (
    ("physical_object_sha256",),
    ("boundary_load_sha256",),
    ("material", "model"),
    ("material", "parameters_sha256"),
    ("material", "rest_coordinate_field_sha256"),
    ("material", "field_seed"),
    ("material", "field_sha256"),
    ("solver_contract_sha256",),
)


def canonical_sha256(value: Any) -> str:
    """Hash JSON data with stable key order and no non-finite numbers."""
    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _is_sha256(value: Any) -> bool:
    return isinstance(value, str) and _SHA256.fullmatch(value) is not None


def _positive_number(value: Any) -> bool:
    return (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and math.isfinite(float(value))
        and float(value) > 0.0
    )


def _nonempty(value: Any) -> bool:
    return isinstance(value, str) and bool(value.strip())


def _at_path(value: Any, path: tuple[str, ...]) -> Any:
    for key in path:
        if not isinstance(value, Mapping) or key not in value:
            return None
        value = value[key]
    return value


def _add_unknown(unknowns: list[dict[str, str]], path: str, reason: str) -> None:
    unknowns.append({"field": path, "status": "UNKNOWN", "reason": reason})


def _check_hash(
    value: Any, path: str, errors: list[dict[str, str]], unknowns: list[dict[str, str]]
) -> bool:
    if value is None:
        _add_unknown(unknowns, path, "required SHA-256 evidence was not supplied")
        return False
    if not _is_sha256(value):
        errors.append({"field": path, "reason": "must be a 64-character hexadecimal SHA-256"})
        return False
    return True


def _physical_case_check(
    physical_case: Any, errors: list[dict[str, str]], unknowns: list[dict[str, str]]
) -> None:
    if not isinstance(physical_case, Mapping):
        _add_unknown(unknowns, "physical_case", "physical case identity is missing")
        return
    for path in _PHYSICAL_FIELDS:
        value = _at_path(physical_case, path)
        label = "physical_case." + ".".join(path)
        if value is None:
            _add_unknown(unknowns, label, "required fixed physical-case field was not supplied")
        elif path[-1].endswith("_sha256"):
            _check_hash(value, label, errors, unknowns)
        elif path == ("material", "model") and not _nonempty(value):
            errors.append({"field": label, "reason": "must be a non-empty material model identifier"})
        elif path == ("material", "field_seed") and (
            isinstance(value, bool) or not isinstance(value, (str, int)) or not str(value).strip()
        ):
            errors.append({"field": label, "reason": "must be a non-empty string or integer seed"})


def _registration_check(
    registration: Any, errors: list[dict[str, str]], unknowns: list[dict[str, str]]
) -> bool:
    if not isinstance(registration, Mapping):
        _add_unknown(unknowns, "registration", "pre-registration record is missing")
        return False
    if registration.get("status") != "PRE_REGISTERED":
        _add_unknown(unknowns, "registration.status", "study is not declared PRE_REGISTERED")
        return False
    if not _nonempty(registration.get("frozen_at_utc")):
        _add_unknown(unknowns, "registration.frozen_at_utc", "freeze time was not supplied")
        return False
    return _check_hash(
        registration.get("protocol_sha256"),
        "registration.protocol_sha256",
        errors,
        unknowns,
    )


def _signature_record(
    discretization: Mapping[str, Any],
    run_id: str,
    errors: list[dict[str, str]],
    unknowns: list[dict[str, str]],
) -> dict[str, Any]:
    signatures = discretization.get("signatures")
    if not isinstance(signatures, Mapping):
        _add_unknown(unknowns, f"runs[{run_id}].discretization.signatures", "discrete signatures were not supplied")
        return {"artifact_sha256": None, "rotation_invariant": None}
    artifact = signatures.get("artifact_sha256")
    artifact_valid = _check_hash(
        artifact,
        f"runs[{run_id}].discretization.signatures.artifact_sha256",
        errors,
        unknowns,
    )
    invariant = signatures.get("rotation_invariant_discretization")
    invariant_digest = None
    invariant_declared = False
    if invariant is not None:
        if not isinstance(invariant, Mapping):
            errors.append(
                {
                    "field": f"runs[{run_id}].discretization.signatures.rotation_invariant_discretization",
                    "reason": "must be an object",
                }
            )
        else:
            invariant_hash_ok = _check_hash(
                invariant.get("sha256"),
                f"runs[{run_id}].discretization.signatures.rotation_invariant_discretization.sha256",
                errors,
                unknowns,
            )
            source_hash_ok = _check_hash(
                invariant.get("source_artifact_sha256"),
                f"runs[{run_id}].discretization.signatures.rotation_invariant_discretization.source_artifact_sha256",
                errors,
                unknowns,
            )
            evidence_hash_ok = _check_hash(
                invariant.get("calculation_evidence_sha256"),
                f"runs[{run_id}].discretization.signatures.rotation_invariant_discretization.calculation_evidence_sha256",
                errors,
                unknowns,
            )
            if not _nonempty(invariant.get("algorithm")):
                _add_unknown(
                    unknowns,
                    f"runs[{run_id}].discretization.signatures.rotation_invariant_discretization.algorithm",
                    "signature algorithm was not supplied",
                )
            if (
                source_hash_ok
                and artifact_valid
                and invariant.get("source_artifact_sha256", "").lower() != artifact.lower()
            ):
                errors.append(
                    {
                        "field": f"runs[{run_id}].discretization.signatures.rotation_invariant_discretization.source_artifact_sha256",
                        "reason": "must identify this run's exact mesh or point-cloud artifact",
                    }
                )
            invariant_declared = invariant_hash_ok and source_hash_ok and evidence_hash_ok and _nonempty(
                invariant.get("algorithm")
            )
            if invariant_declared:
                invariant_digest = invariant["sha256"].lower()
    return {
        "artifact_sha256": artifact.lower() if artifact_valid else None,
        "rotation_invariant": invariant_digest,
        "rotation_invariant_declared": invariant_declared,
        # The checker has no mesh/point payload reader yet. Caller-provided
        # hashes and calculation receipts remain declarations, not verification.
        "rotation_invariant_verification": "UNVERIFIED" if invariant_declared else "UNKNOWN",
    }


def _run_checks(
    runs: Any, errors: list[dict[str, str]], unknowns: list[dict[str, str]]
) -> list[dict[str, Any]]:
    if not isinstance(runs, list):
        errors.append({"field": "runs", "reason": "must be an array"})
        return []
    run_checks: list[dict[str, Any]] = []
    seen_ids: set[str] = set()
    for index, run in enumerate(runs):
        prefix = f"runs[{index}]"
        if not isinstance(run, Mapping):
            errors.append({"field": prefix, "reason": "must be an object"})
            continue
        run_id = run.get("run_id")
        if not _nonempty(run_id):
            errors.append({"field": f"{prefix}.run_id", "reason": "must be a non-empty string"})
            run_id = f"index-{index}"
        if run_id in seen_ids:
            errors.append({"field": f"{prefix}.run_id", "reason": "must be unique within the study"})
        seen_ids.add(run_id)

        direction = run.get("direction")
        if not _nonempty(direction):
            _add_unknown(unknowns, f"{prefix}.direction", "meshing/point-generation direction label is missing")

        disc = run.get("discretization")
        if not isinstance(disc, Mapping):
            errors.append({"field": f"{prefix}.discretization", "reason": "must be an object"})
            disc = {}
        kind = disc.get("kind")
        if kind not in {"cohesive_tetra_mesh", "pd_point_cloud"}:
            errors.append(
                {
                    "field": f"{prefix}.discretization.kind",
                    "reason": "must be cohesive_tetra_mesh or pd_point_cloud",
                }
            )
        for key in ("generator_id", "generator_sha256"):
            value = disc.get(key)
            label = f"{prefix}.discretization.{key}"
            if value is None:
                _add_unknown(unknowns, label, "generator provenance was not supplied")
            elif key.endswith("_sha256"):
                _check_hash(value, label, errors, unknowns)
            elif not _nonempty(value):
                errors.append({"field": label, "reason": "must be a non-empty string"})

        seed = disc.get("generation_seed")
        if seed is None:
            _add_unknown(unknowns, f"{prefix}.discretization.generation_seed", "generation seed was not supplied")
        elif isinstance(seed, bool) or not isinstance(seed, (int, str)) or not str(seed).strip():
            errors.append(
                {
                    "field": f"{prefix}.discretization.generation_seed",
                    "reason": "must be a non-empty string or integer",
                }
            )

        spacing = disc.get("spacing")
        expected_symbol = "h" if kind == "cohesive_tetra_mesh" else "dx"
        spacing_value = None
        if not isinstance(spacing, Mapping):
            _add_unknown(unknowns, f"{prefix}.discretization.spacing", "explicit h or dx was not supplied")
        else:
            if spacing.get("symbol") != expected_symbol:
                errors.append(
                    {
                        "field": f"{prefix}.discretization.spacing.symbol",
                        "reason": f"{kind} requires explicit {expected_symbol}",
                    }
                )
            if spacing.get("value_m") is None:
                _add_unknown(
                    unknowns,
                    f"{prefix}.discretization.spacing.value_m",
                    "positive explicit spacing in metres was not supplied",
                )
            elif not _positive_number(spacing.get("value_m")):
                errors.append(
                    {
                        "field": f"{prefix}.discretization.spacing.value_m",
                        "reason": "must be finite and positive",
                    }
                )
            else:
                spacing_value = float(spacing["value_m"])

        horizon = None
        process_zone_status = "NOT_APPLICABLE"
        if kind == "pd_point_cloud":
            if disc.get("horizon_m") is None:
                _add_unknown(
                    unknowns,
                    f"{prefix}.discretization.horizon_m",
                    "native PD horizon is unknown; it is not inferred from mesh edges or thickness",
                )
            elif _positive_number(disc.get("horizon_m")):
                horizon = float(disc["horizon_m"])
            else:
                errors.append(
                    {
                        "field": f"{prefix}.discretization.horizon_m",
                        "reason": "a supplied native PD horizon must be finite and positive",
                    }
                )
            process_zone = disc.get("process_zone")
            if process_zone is None:
                process_zone_status = "UNKNOWN"
                _add_unknown(
                    unknowns,
                    f"{prefix}.discretization.process_zone",
                    "process-zone length/definition/source is unknown; no wall-layer count is inferred",
                )
            elif not isinstance(process_zone, Mapping):
                errors.append(
                    {
                        "field": f"{prefix}.discretization.process_zone",
                        "reason": "must be an object when supplied",
                    }
                )
            elif process_zone.get("status") == "UNKNOWN":
                process_zone_status = "UNKNOWN"
                _add_unknown(
                    unknowns,
                    f"{prefix}.discretization.process_zone",
                    "process-zone length/definition/source is unknown; no wall-layer count is inferred",
                )
            elif process_zone.get("status") != "KNOWN":
                errors.append(
                    {
                        "field": f"{prefix}.discretization.process_zone.status",
                        "reason": "must be KNOWN or UNKNOWN",
                    }
                )
            elif not _positive_number(process_zone.get("length_m")):
                errors.append(
                    {
                        "field": f"{prefix}.discretization.process_zone.length_m",
                        "reason": "a KNOWN process-zone length must be finite and positive",
                    }
                )
            else:
                definition_ok = _check_hash(
                    process_zone.get("definition_sha256"),
                    f"{prefix}.discretization.process_zone.definition_sha256",
                    errors,
                    unknowns,
                )
                source_ok = _check_hash(
                    process_zone.get("source_sha256"),
                    f"{prefix}.discretization.process_zone.source_sha256",
                    errors,
                    unknowns,
                )
                if definition_ok and source_ok:
                    process_zone_status = "DECLARED_KNOWN_UNVERIFIED"

        signatures = _signature_record(disc, str(run_id), errors, unknowns)
        lineage = disc.get("construction_lineage")
        lineage_kind = None
        parent_signature = None
        if not isinstance(lineage, Mapping):
            _add_unknown(unknowns, f"{prefix}.discretization.construction_lineage", "construction lineage was not supplied")
        else:
            lineage_kind = lineage.get("kind")
            if lineage_kind not in {"independent_generation", "rigid_rotation", "refinement", "unknown"}:
                _add_unknown(
                    unknowns,
                    f"{prefix}.discretization.construction_lineage.kind",
                    "lineage kind is missing or unrecognized",
                )
                lineage_kind = "unknown"
            parent_signature = lineage.get("parent_artifact_sha256")
            if lineage_kind == "rigid_rotation":
                if _check_hash(
                    parent_signature,
                    f"{prefix}.discretization.construction_lineage.parent_artifact_sha256",
                    errors,
                    unknowns,
                ):
                    parent_signature = parent_signature.lower()
                transform_hash = lineage.get("transform_sha256")
                _check_hash(
                    transform_hash,
                    f"{prefix}.discretization.construction_lineage.transform_sha256",
                    errors,
                    unknowns,
                )
            elif lineage_kind == "independent_generation":
                source_geometry_hash = lineage.get("source_geometry_sha256")
                _check_hash(
                    source_geometry_hash,
                    f"{prefix}.discretization.construction_lineage.source_geometry_sha256",
                    errors,
                    unknowns,
                )

        solver = run.get("solver")
        dt = None
        end_time = None
        if not isinstance(solver, Mapping):
            _add_unknown(unknowns, f"{prefix}.solver", "solver settings were not supplied")
        else:
            if solver.get("dt_s") is None:
                _add_unknown(unknowns, f"{prefix}.solver.dt_s", "positive actual or planned dt was not supplied")
            elif _positive_number(solver.get("dt_s")):
                dt = float(solver["dt_s"])
            else:
                errors.append({"field": f"{prefix}.solver.dt_s", "reason": "must be finite and positive"})
            if solver.get("end_time_s") is None:
                _add_unknown(unknowns, f"{prefix}.solver.end_time_s", "positive comparison end time was not supplied")
            elif _positive_number(solver.get("end_time_s")):
                end_time = float(solver["end_time_s"])
            else:
                errors.append(
                    {"field": f"{prefix}.solver.end_time_s", "reason": "must be finite and positive"}
                )

        run_checks.append(
            {
                "run_id": str(run_id),
                "direction": direction if _nonempty(direction) else None,
                "kind": kind,
                "generator_id": disc.get("generator_id"),
                "generator_sha256": (
                    disc.get("generator_sha256", "").lower()
                    if _is_sha256(disc.get("generator_sha256"))
                    else None
                ),
                "generation_seed": seed,
                "spacing_m": spacing_value,
                "horizon_m": horizon,
                "process_zone_status": process_zone_status,
                "dt_s": dt,
                "end_time_s": end_time,
                "lineage_kind": lineage_kind,
                "parent_artifact_sha256": parent_signature,
                **signatures,
                "raw": run,
            }
        )
    return run_checks


def validate_manifest(manifest: Any) -> dict[str, Any]:
    """Validate a study manifest without running mechanics or trusting outcomes."""
    errors: list[dict[str, str]] = []
    unknowns: list[dict[str, str]] = []
    if not isinstance(manifest, Mapping):
        return {
            "schema": "pgfrag.bias-validation.v1",
            "valid": False,
            "errors": [{"field": "$", "reason": "manifest must be an object"}],
            "unknowns": [],
            "manifest_sha256": None,
        }
    if manifest.get("schema") != SCHEMA:
        errors.append({"field": "schema", "reason": f"must be {SCHEMA}"})
    if not _nonempty(manifest.get("study_id")):
        errors.append({"field": "study_id", "reason": "must be a non-empty string"})
    _registration_check(manifest.get("registration"), errors, unknowns)
    _physical_case_check(manifest.get("physical_case"), errors, unknowns)
    run_checks = _run_checks(manifest.get("runs"), errors, unknowns)

    physical_case = manifest.get("physical_case")
    object_hash = (
        physical_case.get("physical_object_sha256")
        if isinstance(physical_case, Mapping)
        else None
    )
    for run in run_checks:
        raw = run["raw"]
        raw_disc = raw.get("discretization")
        lineage = raw_disc.get("construction_lineage") if isinstance(raw_disc, Mapping) else None
        if (
            isinstance(lineage, Mapping)
            and lineage.get("kind") == "independent_generation"
            and _is_sha256(lineage.get("source_geometry_sha256"))
            and _is_sha256(object_hash)
            and lineage["source_geometry_sha256"].lower() != object_hash.lower()
        ):
            errors.append(
                {
                    "field": f"runs[{run['run_id']}].discretization.construction_lineage.source_geometry_sha256",
                    "reason": "must match physical_case.physical_object_sha256",
                }
            )

    try:
        digest = canonical_sha256(manifest)
    except (TypeError, ValueError):
        digest = None
        errors.append({"field": "$", "reason": "manifest is not finite, canonical JSON data"})
    return {
        "schema": "pgfrag.bias-validation.v1",
        "valid": not errors,
        "errors": errors,
        "unknowns": unknowns,
        "manifest_sha256": digest,
        "run_checks": [
            {key: value for key, value in item.items() if key != "raw"}
            for item in run_checks
        ],
    }


def _physical_signature(physical_case: Any) -> tuple[Any, ...] | None:
    if not isinstance(physical_case, Mapping):
        return None
    signature = tuple(_at_path(physical_case, path) for path in _PHYSICAL_FIELDS)
    if any(value is None for value in signature):
        return None
    return signature


def _same_physical_case(left: Any, right: Any) -> str:
    left_signature = _physical_signature(left)
    right_signature = _physical_signature(right)
    if left_signature is None or right_signature is None:
        return "UNKNOWN"
    return "MATCH" if left_signature == right_signature else "MISMATCH"


def _same_float(left: Any, right: Any) -> bool:
    return isinstance(left, (int, float)) and isinstance(right, (int, float)) and float(left) == float(right)


def _lineage_points_to(left: Mapping[str, Any], right: Mapping[str, Any]) -> bool:
    return (
        left.get("lineage_kind") == "rigid_rotation"
        and left.get("parent_artifact_sha256") is not None
        and left.get("parent_artifact_sha256") == right.get("artifact_sha256")
    ) or (
        right.get("lineage_kind") == "rigid_rotation"
        and right.get("parent_artifact_sha256") is not None
        and right.get("parent_artifact_sha256") == left.get("artifact_sha256")
    )


def _pair_audit(left: Mapping[str, Any], right: Mapping[str, Any]) -> dict[str, Any]:
    left_artifact = left.get("artifact_sha256")
    right_artifact = right.get("artifact_sha256")
    if left_artifact is not None and left_artifact == right_artifact:
        relation = "DUPLICATE_DISCRETIZATION"
        evidence = "EXACT_ARTIFACT_SHA256_MATCH"
    elif _lineage_points_to(left, right):
        relation = "SAME_SOURCE_RIGID_ROTATION"
        evidence = "DECLARED_LINEAGE_LINK_MATCHES_RECORDED_ARTIFACT"
    elif (
        left.get("rotation_invariant") is not None
        and left.get("rotation_invariant") == right.get("rotation_invariant")
    ):
        relation = "SAME_ROTATION_INVARIANT_SIGNATURE"
        evidence = "UNVERIFIED_CALLER_SIGNATURE_MATCH"
    elif left_artifact is None or right_artifact is None:
        relation = "UNKNOWN"
        evidence = "MISSING_ARTIFACT_SIGNATURE"
    elif left.get("rotation_invariant") is None or right.get("rotation_invariant") is None:
        relation = "UNKNOWN"
        evidence = "MISSING_ROTATION_INVARIANT_SIGNATURE"
    else:
        relation = "DIFFERENT_DECLARED_DISCRETIZATION_SIGNATURES"
        evidence = "UNVERIFIED_CALLER_SIGNATURES"
    return {
        "left_run_id": left.get("run_id"),
        "right_run_id": right.get("run_id"),
        "relationship": relation,
        "evidence_status": evidence,
        "independent_discretization_evidence_status": (
            "UNVERIFIED" if relation == "DIFFERENT_DECLARED_DISCRETIZATION_SIGNATURES" else "UNKNOWN"
        ),
    }


def _is_generation_pair(left: Mapping[str, Any], right: Mapping[str, Any]) -> tuple[bool, str]:
    relation = _pair_audit(left, right)["relationship"]
    if relation in {
        "DUPLICATE_DISCRETIZATION",
        "SAME_SOURCE_RIGID_ROTATION",
        "SAME_ROTATION_INVARIANT_SIGNATURE",
    }:
        return False, relation
    if left.get("lineage_kind") != "independent_generation" or right.get("lineage_kind") != "independent_generation":
        return False, "UNKNOWN_CONSTRUCTION_LINEAGE"
    if (
        not _nonempty(left.get("generator_id"))
        or left.get("generator_id") != right.get("generator_id")
        or left.get("generator_sha256") is None
        or left.get("generator_sha256") != right.get("generator_sha256")
    ):
        return False, "GENERATOR_MISMATCH_OR_UNKNOWN"
    if relation == "UNKNOWN":
        return False, "UNKNOWN_SIGNATURE_RELATION"
    left_seed = left.get("generation_seed")
    right_seed = right.get("generation_seed")
    if left_seed is None or right_seed is None or str(left_seed).strip() == str(right_seed).strip():
        return False, "SAME_GENERATION_SEED"
    return True, "UNVERIFIED_INDEPENDENT_CANDIDATE"


def _run_metrics(run: Mapping[str, Any]) -> dict[str, Any]:
    results = run.get("results")
    results = results if isinstance(results, Mapping) else {}
    actual = results.get("actual_run")
    actual = actual if isinstance(actual, Mapping) else {}
    surface = results.get("qualified_surface")
    surface = surface if isinstance(surface, Mapping) else {}
    receipt = actual.get("receipt_sha256")
    outputs = actual.get("outputs")
    output_hashes: set[str] = set()
    actual_is_real = (
        actual.get("status") == "COMPLETED"
        and _is_sha256(receipt)
        and isinstance(outputs, list)
    )
    if actual_is_real:
        for output in outputs:
            if isinstance(output, Mapping) and _is_sha256(output.get("sha256")):
                output_hashes.add(output["sha256"].lower())
        actual_is_real = bool(output_hashes)
    surface_hash = surface.get("sha256")
    surface_is_qualified = (
        surface.get("status") == "QUALIFIED"
        and _is_sha256(surface_hash)
        and _is_sha256(surface.get("qualification_report_sha256"))
        and _is_sha256(surface.get("source_run_receipt_sha256"))
        and actual_is_real
        and surface.get("source_run_receipt_sha256", "").lower() == receipt.lower()
    )
    declared_metrics = results.get("metrics")
    declared_metrics = declared_metrics if isinstance(declared_metrics, Mapping) else {}
    reports: dict[str, Any] = {}
    for name in _METRIC_NAMES:
        metric = declared_metrics.get(name)
        metric = metric if isinstance(metric, Mapping) else {}
        provenance = metric.get("provenance")
        provenance = provenance if isinstance(provenance, Mapping) else {}
        source_hash = provenance.get("source_output_sha256")
        source_run = provenance.get("actual_run_receipt_sha256")
        source_surface = provenance.get("qualified_surface_sha256")
        evidence_ok = (
            actual_is_real
            and surface_is_qualified
            and _is_sha256(source_hash)
            and source_hash.lower() in output_hashes
            and _is_sha256(source_run)
            and source_run.lower() == receipt.lower()
            and _is_sha256(source_surface)
            and source_surface.lower() == surface_hash.lower()
            and provenance.get("scope") == "physical"
            and "value" in metric
        )
        reports[name] = {
            "status": "PROVENANCE_LINKED_UNVERIFIED" if evidence_ok else "NOT_RUN",
            "verification_status": (
                "DECLARED_HASH_LINKS_UNVERIFIED" if evidence_ok else "INCOMPLETE_OR_MISMATCHED"
            ),
            "value": metric.get("value") if evidence_ok else None,
            "reason": (
                "hash references link the declared value to a completed run output and qualified surface; source bytes were not rehashed"
                if evidence_ok
                else "actual-run output and qualified-surface provenance are incomplete or mismatched"
            ),
        }
    return reports


def comparison_readiness(manifests: Any) -> dict[str, Any]:
    """Report preregistered design readiness and evidence status for one or more studies."""
    if isinstance(manifests, Mapping):
        study_manifests = [manifests]
    elif isinstance(manifests, Sequence) and not isinstance(manifests, (str, bytes)):
        study_manifests = list(manifests)
    else:
        study_manifests = []

    validations = [validate_manifest(item) for item in study_manifests]
    errors = [
        {"study_index": index, **error}
        for index, validation in enumerate(validations)
        for error in validation["errors"]
    ]
    unknowns = [
        {"study_index": index, **unknown}
        for index, validation in enumerate(validations)
        for unknown in validation["unknowns"]
    ]
    if not study_manifests:
        errors.append({"study_index": 0, "field": "$", "reason": "at least one manifest is required"})

    physical_status = (
        "UNKNOWN"
        if any(
            _physical_signature(
                item.get("physical_case") if isinstance(item, Mapping) else None
            )
            is None
            for item in study_manifests
        )
        else "MATCH"
    )
    for left_index in range(len(study_manifests)):
        for right_index in range(left_index + 1, len(study_manifests)):
            status = _same_physical_case(
                study_manifests[left_index].get("physical_case")
                if isinstance(study_manifests[left_index], Mapping)
                else None,
                study_manifests[right_index].get("physical_case")
                if isinstance(study_manifests[right_index], Mapping)
                else None,
            )
            if status == "MISMATCH":
                physical_status = "MISMATCH"
            elif status == "UNKNOWN" and physical_status != "MISMATCH":
                physical_status = "UNKNOWN"

    flattened: list[dict[str, Any]] = []
    registrations_ready = True
    pd_resolution_status = "NOT_APPLICABLE"
    for index, manifest in enumerate(study_manifests):
        validation = validations[index]
        if not isinstance(manifest, Mapping):
            registrations_ready = False
            continue
        registration = manifest.get("registration")
        if not (
            isinstance(registration, Mapping)
            and registration.get("status") == "PRE_REGISTERED"
            and _is_sha256(registration.get("protocol_sha256"))
            and _nonempty(registration.get("frozen_at_utc"))
        ):
            registrations_ready = False
        raw_runs = manifest.get("runs")
        if not isinstance(raw_runs, list):
            continue
        run_checks = {item["run_id"]: item for item in validation.get("run_checks", [])}
        for run in raw_runs:
            if not isinstance(run, Mapping):
                continue
            run_id = run.get("run_id")
            check = run_checks.get(str(run_id), {})
            entry = {**check, "physical_case": manifest.get("physical_case"), "raw": run}
            flattened.append(entry)
            if check.get("kind") == "pd_point_cloud":
                if pd_resolution_status == "NOT_APPLICABLE":
                    pd_resolution_status = "UNVERIFIED"
                if (
                    check.get("horizon_m") is None
                    or check.get("process_zone_status") != "DECLARED_KNOWN_UNVERIFIED"
                ):
                    pd_resolution_status = "UNKNOWN"

    pair_audits = [
        _pair_audit(flattened[i], flattened[j])
        for i in range(len(flattened))
        for j in range(i + 1, len(flattened))
    ]
    direction_design = "NOT_READY"
    direction_pair: dict[str, Any] | None = None
    for i, left in enumerate(flattened):
        for right in flattened[i + 1 :]:
            if (
                left.get("kind") != right.get("kind")
                or left.get("generator_id") != right.get("generator_id")
                or left.get("generator_sha256") != right.get("generator_sha256")
                or not _same_float(left.get("spacing_m"), right.get("spacing_m"))
                or not _same_float(left.get("dt_s"), right.get("dt_s"))
                or not _nonempty(left.get("direction"))
                or not _nonempty(right.get("direction"))
                or left.get("direction") == right.get("direction")
            ):
                continue
            pair_ok, reason = _is_generation_pair(left, right)
            if pair_ok:
                direction_design = "READY"
                direction_pair = {
                    "run_ids": [left.get("run_id"), right.get("run_id")],
                    "status": "UNVERIFIED_INDEPENDENT_CANDIDATE",
                    "reason": "different declared generation seeds and directions; signatures were not recomputed",
                }
                break
            if reason not in {"SAME_GENERATION_SEED", "UNKNOWN_SIGNATURE_RELATION"}:
                direction_design = "UNKNOWN" if direction_design != "READY" else direction_design
        if direction_design == "READY":
            break

    spatial_design = "NOT_READY"
    spatial_pair: dict[str, Any] | None = None
    for i, left in enumerate(flattened):
        for right in flattened[i + 1 :]:
            if (
                left.get("kind") != right.get("kind")
                or not _nonempty(left.get("generator_id"))
                or left.get("generator_id") != right.get("generator_id")
                or left.get("generator_sha256") is None
                or left.get("generator_sha256") != right.get("generator_sha256")
                or left.get("direction") != right.get("direction")
                or _same_float(left.get("spacing_m"), right.get("spacing_m"))
                or not _same_float(left.get("dt_s"), right.get("dt_s"))
            ):
                continue
            relation = _pair_audit(left, right)["relationship"]
            if relation == "DIFFERENT_DECLARED_DISCRETIZATION_SIGNATURES":
                spatial_design = "READY"
                spatial_pair = {
                    "run_ids": [left.get("run_id"), right.get("run_id")],
                    "spacing_m": [left.get("spacing_m"), right.get("spacing_m")],
                    "signature_status": "UNVERIFIED",
                }
                break
            if relation == "UNKNOWN":
                spatial_design = "UNKNOWN"
        if spatial_design == "READY":
            break

    dt_design = "NOT_READY"
    dt_pair: dict[str, Any] | None = None
    for i, left in enumerate(flattened):
        for right in flattened[i + 1 :]:
            if (
                left.get("artifact_sha256") is None
                or left.get("artifact_sha256") != right.get("artifact_sha256")
                or not _same_float(left.get("spacing_m"), right.get("spacing_m"))
                or not _same_float(left.get("end_time_s"), right.get("end_time_s"))
                or left.get("dt_s") is None
                or right.get("dt_s") is None
                or _same_float(left.get("dt_s"), right.get("dt_s"))
            ):
                continue
            left_disc = left.get("raw", {}).get("discretization", {})
            right_disc = right.get("raw", {}).get("discretization", {})
            if not isinstance(left_disc, Mapping) or not isinstance(right_disc, Mapping):
                continue
            if left_disc.get("horizon_m") != right_disc.get("horizon_m"):
                continue
            if left_disc.get("process_zone") != right_disc.get("process_zone"):
                continue
            dt_design = "READY"
            dt_pair = {
                "run_ids": [left.get("run_id"), right.get("run_id")],
                "dt_s": [left.get("dt_s"), right.get("dt_s")],
                "same_discretization_sha256": left.get("artifact_sha256"),
            }
            break
        if dt_design == "READY":
            break

    design_checks = {
        "independent_generated_direction_seed": {
            "status": direction_design,
            "candidate": direction_pair,
            "interpretation": "a common rigid rotation cannot satisfy this check",
        },
        "spatial_resolution": {"status": spatial_design, "candidate": spatial_pair},
        "time_step": {"status": dt_design, "candidate": dt_pair},
    }
    required_design_ready = all(
        design_checks[key]["status"] == "READY" for key in design_checks
    )
    has_unknown_design = any(
        design_checks[key]["status"] == "UNKNOWN" for key in design_checks
    )
    if physical_status == "MISMATCH":
        status = "NOT_COMPARABLE"
    elif errors:
        status = "INVALID_MANIFEST"
    elif physical_status == "UNKNOWN":
        status = "NOT_READY"
    elif not registrations_ready:
        status = "NOT_READY"
    elif required_design_ready:
        any_run_evidence = any(
            isinstance(run.get("results"), Mapping)
            and isinstance(run["results"].get("actual_run"), Mapping)
            and run["results"]["actual_run"].get("status") not in (None, "NOT_RUN")
            for item in flattened
            for run in [item.get("raw", {})]
            if isinstance(run, Mapping)
        )
        status = (
            "RUN_EVIDENCE_PRESENT_ANALYSIS_PENDING"
            if any_run_evidence
            else "READY_NOT_RUN"
        )
    elif has_unknown_design or unknowns:
        status = "NOT_READY"
    else:
        status = "NOT_READY"

    metrics = {
        str(item.get("run_id")): _run_metrics(item.get("raw", {}))
        for item in flattened
    }
    has_rotation_pair = any(
        item["relationship"] == "SAME_SOURCE_RIGID_ROTATION"
        for item in pair_audits
    )
    return {
        "schema": READINESS_SCHEMA,
        "status": status,
        "physical_case_comparability": physical_status,
        "design_checks": design_checks,
        "independent_discretization_evidence_status": "UNVERIFIED",
        "peridynamic_resolution_qualification": pd_resolution_status,
        "pair_audits": pair_audits,
        "metric_evidence": metrics,
        "claims": {
            "objectivity": "DESIGNED_NOT_RUN" if has_rotation_pair else "NOT_TESTED",
            "mesh_or_point_cloud_independence": "NOT_ESTABLISHED",
            "directional_bias_removed": "NOT_ESTABLISHED",
            "fragment_count_target": "NOT_REQUIRED",
            "result_artifact_hashes": "DECLARED_NOT_RECOMPUTED",
        },
        "representation_limits": {
            "cohesive_crack_path": "TETRAHEDRAL_FACET_CONSTRAINED",
            "peridynamic_point_independence": "NOT_ASSUMED",
        },
        "errors": errors,
        "unknowns": unknowns,
        "manifest_validation": validations,
    }


__all__ = ["SCHEMA", "READINESS_SCHEMA", "canonical_sha256", "validate_manifest", "comparison_readiness"]
