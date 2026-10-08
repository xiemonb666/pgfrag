from __future__ import annotations

import copy
import hashlib

from pgfrag.bias_protocol import comparison_readiness, validate_manifest


def digest(label: str) -> str:
    return hashlib.sha256(label.encode("utf-8")).hexdigest()


PHYSICAL_OBJECT = digest("physical-object-v1")


def make_run(
    run_id: str,
    *,
    direction: str,
    seed: int,
    dx_m: float,
    dt_s: float,
    artifact_key: str,
    invariant_key: str,
    process_zone: dict | None = None,
    results: dict | None = None,
) -> dict:
    artifact_sha = digest(f"artifact:{artifact_key}")
    return {
        "run_id": run_id,
        "direction": direction,
        "discretization": {
            "kind": "pd_point_cloud",
            "generator_id": "review-point-generator",
            "generator_sha256": digest("generator-source-v1"),
            "generation_seed": seed,
            "spacing": {"symbol": "dx", "value_m": dx_m},
            # Omitted horizon/process-zone values intentionally remain UNKNOWN.
            **({"process_zone": process_zone} if process_zone is not None else {}),
            "signatures": {
                "artifact_sha256": artifact_sha,
                "rotation_invariant_discretization": {
                    "sha256": digest(f"invariant:{invariant_key}"),
                    "algorithm": "declared-canonical-point-set-v1",
                    "source_artifact_sha256": artifact_sha,
                    "calculation_evidence_sha256": digest(f"calculation-record:{artifact_key}"),
                },
            },
            "construction_lineage": {
                "kind": "independent_generation",
                "source_geometry_sha256": PHYSICAL_OBJECT,
            },
        },
        "solver": {"dt_s": dt_s, "end_time_s": 1.0e-4},
        "results": results or {
            "actual_run": {"status": "NOT_RUN"},
            "qualified_surface": {"status": "NOT_AVAILABLE"},
            "metrics": {},
        },
    }


def make_manifest() -> dict:
    coarse_x = make_run(
        "coarse-x", direction="x", seed=11, dx_m=1.0e-3, dt_s=1.0e-6,
        artifact_key="coarse-x-seed11", invariant_key="coarse-x-seed11",
    )
    coarse_y = make_run(
        "coarse-y", direction="y", seed=29, dx_m=1.0e-3, dt_s=1.0e-6,
        artifact_key="coarse-y-seed29", invariant_key="coarse-y-seed29",
    )
    fine_x = make_run(
        "fine-x", direction="x", seed=47, dx_m=0.5e-3, dt_s=1.0e-6,
        artifact_key="fine-x-seed47", invariant_key="fine-x-seed47",
    )
    coarse_x_half_dt = copy.deepcopy(coarse_x)
    coarse_x_half_dt["run_id"] = "coarse-x-half-dt"
    coarse_x_half_dt["solver"]["dt_s"] = 0.5e-6
    return {
        "schema": "pgfrag.bias-study.v1",
        "study_id": "review-preregistered-design",
        "registration": {
            "status": "PRE_REGISTERED",
            "frozen_at_utc": "2026-10-04T00:00:00Z",
            "protocol_sha256": digest("bias-protocol-spec-v1"),
        },
        "physical_case": {
            "physical_object_sha256": PHYSICAL_OBJECT,
            "boundary_load_sha256": digest("fixed-boundary-and-load-v1"),
            "material": {
                "model": "PD-OSB",
                "parameters_sha256": digest("fixed-material-parameters-v1"),
                "rest_coordinate_field_sha256": digest("fixed-rest-coordinate-field-v1"),
                "field_seed": 101,
                "field_sha256": digest("fixed-material-field-v1"),
            },
            "solver_contract_sha256": digest("fixed-solver-contract-v1"),
        },
        "runs": [coarse_x, coarse_y, fine_x, coarse_x_half_dt],
    }


def test_rotation_of_same_discrete_source_does_not_count_as_independent_remesh() -> None:
    manifest = make_manifest()
    left = manifest["runs"][0]
    rotated = copy.deepcopy(left)
    rotated["run_id"] = "same-source-rotated"
    rotated["direction"] = "rotated-90-deg"
    rotated["discretization"]["generation_seed"] = 999
    rotated["discretization"]["signatures"]["artifact_sha256"] = digest("rotated-file-bytes")
    rotated_signature = rotated["discretization"]["signatures"]["rotation_invariant_discretization"]
    rotated_signature["source_artifact_sha256"] = digest("rotated-file-bytes")
    rotated_signature["calculation_evidence_sha256"] = digest("rotated-calculation-record")
    # The claimed rotation-invariant value still matches the unchanged source.
    manifest["runs"] = [left, rotated]

    report = comparison_readiness(manifest)

    matching = [
        item for item in report["pair_audits"]
        if {item["left_run_id"], item["right_run_id"]} == {"coarse-x", "same-source-rotated"}
    ]
    assert matching[0]["relationship"] == "SAME_ROTATION_INVARIANT_SIGNATURE"
    assert matching[0]["evidence_status"] == "UNVERIFIED_CALLER_SIGNATURE_MATCH"
    assert report["design_checks"]["independent_generated_direction_seed"]["status"] != "READY"
    assert report["claims"]["objectivity"] == "NOT_TESTED"

    lineage_rotation = copy.deepcopy(rotated)
    lineage_rotation["run_id"] = "lineage-linked-rotation"
    lineage_rotation["discretization"]["construction_lineage"] = {
        "kind": "rigid_rotation",
        "parent_artifact_sha256": left["discretization"]["signatures"]["artifact_sha256"],
        "transform_sha256": digest("rigid-rotation-transform"),
    }
    manifest["runs"] = [left, lineage_rotation]
    lineage_report = comparison_readiness(manifest)
    assert lineage_report["pair_audits"][0]["relationship"] == "SAME_SOURCE_RIGID_ROTATION"
    assert lineage_report["claims"]["objectivity"] == "DESIGNED_NOT_RUN"


def test_same_artifact_with_different_seed_is_still_a_duplicate() -> None:
    manifest = make_manifest()
    left = copy.deepcopy(manifest["runs"][0])
    right = copy.deepcopy(left)
    right["run_id"] = "different-seed-same-bytes"
    right["direction"] = "y"
    right["discretization"]["generation_seed"] = 303
    manifest["runs"] = [left, right]

    report = comparison_readiness(manifest)

    assert report["pair_audits"][0]["relationship"] == "DUPLICATE_DISCRETIZATION"
    assert report["design_checks"]["independent_generated_direction_seed"]["status"] != "READY"


def test_missing_rest_coordinate_field_hash_leaves_physical_comparison_unknown() -> None:
    manifest = make_manifest()
    del manifest["physical_case"]["material"]["rest_coordinate_field_sha256"]

    validation = validate_manifest(manifest)
    report = comparison_readiness(manifest)

    assert validation["valid"] is True
    assert any("rest_coordinate_field_sha256" in item["field"] for item in validation["unknowns"])
    assert report["physical_case_comparability"] == "UNKNOWN"
    assert report["status"] == "NOT_READY"


def test_cross_material_studies_are_not_comparable() -> None:
    first = make_manifest()
    second = copy.deepcopy(first)
    second["study_id"] = "different-material"
    second["physical_case"]["material"]["model"] = "cohesive-linear"

    report = comparison_readiness([first, second])

    assert report["physical_case_comparability"] == "MISMATCH"
    assert report["status"] == "NOT_COMPARABLE"


def test_geometry_only_metric_claim_is_not_physics_evidence_and_zero_is_not_pass() -> None:
    manifest = make_manifest()
    run = manifest["runs"][0]
    run["results"] = {
        "actual_run": {"status": "NOT_RUN"},
        "qualified_surface": {"status": "NOT_AVAILABLE"},
        "metrics": {
            "energy": {
                "value": 0.0,
                "provenance": {
                    "scope": "geometry_only",
                    "source_output_sha256": digest("geometry-output"),
                    "actual_run_receipt_sha256": digest("nonexistent-run"),
                    "qualified_surface_sha256": digest("unqualified-surface"),
                },
            }
        },
    }

    report = comparison_readiness(manifest)
    energy = report["metric_evidence"]["coarse-x"]["energy"]
    fragments = report["metric_evidence"]["coarse-x"]["fragment_distribution"]

    assert energy == {
        "status": "NOT_RUN",
        "verification_status": "INCOMPLETE_OR_MISMATCHED",
        "value": None,
        "reason": "actual-run output and qualified-surface provenance are incomplete or mismatched",
    }
    assert fragments["status"] == "NOT_RUN"
    assert fragments["value"] is None
    assert report["claims"]["fragment_count_target"] == "NOT_REQUIRED"

    receipt_sha = digest("completed-run-receipt")
    output_sha = digest("actual-run-output")
    surface_sha = digest("qualified-surface")
    run["results"] = {
        "actual_run": {
            "status": "COMPLETED",
            "receipt_sha256": receipt_sha,
            "outputs": [{"path": "physics/result.vtu", "sha256": output_sha}],
        },
        "qualified_surface": {
            "status": "QUALIFIED",
            "sha256": surface_sha,
            "qualification_report_sha256": digest("surface-qualification-report"),
            "source_run_receipt_sha256": receipt_sha,
        },
        "metrics": {
            "energy": {
                "value": 0.0,
                "provenance": {
                    "scope": "physical",
                    "source_output_sha256": output_sha,
                    "actual_run_receipt_sha256": receipt_sha,
                    "qualified_surface_sha256": surface_sha,
                },
            }
        },
    }
    evidenced_report = comparison_readiness(manifest)
    evidenced = evidenced_report["metric_evidence"]["coarse-x"]["energy"]
    assert evidenced_report["status"] == "RUN_EVIDENCE_PRESENT_ANALYSIS_PENDING"
    assert evidenced["status"] == "PROVENANCE_LINKED_UNVERIFIED"
    assert evidenced["verification_status"] == "DECLARED_HASH_LINKS_UNVERIFIED"
    assert evidenced["value"] == 0.0


def test_valid_preregistered_matrix_reports_ready_not_run_and_keeps_pd_scale_unknown() -> None:
    manifest = make_manifest()

    validation = validate_manifest(manifest)
    report = comparison_readiness(manifest)

    assert validation["valid"] is True
    assert report["status"] == "READY_NOT_RUN"
    assert report["design_checks"]["independent_generated_direction_seed"]["status"] == "READY"
    assert report["design_checks"]["spatial_resolution"]["status"] == "READY"
    assert report["design_checks"]["time_step"]["status"] == "READY"
    assert report["independent_discretization_evidence_status"] == "UNVERIFIED"
    assert report["peridynamic_resolution_qualification"] == "UNKNOWN"
    assert report["claims"]["mesh_or_point_cloud_independence"] == "NOT_ESTABLISHED"
    assert all(
        metric["status"] == "NOT_RUN"
        for run_metrics in report["metric_evidence"].values()
        for metric in run_metrics.values()
    )
