from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import sys
import unittest


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from fracture_artifact_contract import (  # noqa: E402
    ArtifactContractError,
    export_artifact,
    read_clean_ground_truth_for_audit,
    read_observation,
    validate_artifact,
)


FIXTURES = HERE / "fixtures"


def load_case(name: str) -> dict:
    return json.loads((FIXTURES / name).read_text(encoding="utf-8"))


def error_codes(report) -> set[str]:
    return {item["code"] for item in report.errors}


def count_fixture_rows(document: dict) -> int:
    return sum(len(document[key]) for key in ("fragments", "fracture_surface_pairs", "events")) + len(
        document["observation"]["fragment_observations"]
    ) + len(document["clean_ground_truth"]["labels"])


class FractureArtifactContractTests(unittest.TestCase):
    def test_complete_geometry_fixture_is_valid_but_never_pretraining_qualified(self) -> None:
        document = load_case("complete_fixture.json")
        report = validate_artifact(document)
        self.assertTrue(report.valid, report.as_dict())
        self.assertEqual(report.qualification_class, "FIXTURE_GEOMETRY_ONLY")
        self.assertFalse(report.pretraining_qualified)
        self.assertLessEqual(count_fixture_rows(document), 64)

    def test_missing_labels_is_distinct_from_explicit_unknown_and_genuine_zero(self) -> None:
        document = load_case("unknown_zero_distinction.json")
        report = validate_artifact(document)
        self.assertTrue(report.valid, report.as_dict())
        self.assertEqual(document["material"]["parameters"]["fracture_toughness_j_m2"]["state"], "UNKNOWN")
        self.assertEqual(document["load"]["peak_force_n"], {"state": "KNOWN", "value": 0.0, "unit": "N"})
        self.assertEqual(document["observation"]["fragment_observations"][1]["sample_count"]["value"], 0.0)

        missing = validate_artifact(load_case("missinglabels.json"))
        self.assertFalse(missing.valid)
        self.assertIn("MISSING_FIELD", error_codes(missing))
        self.assertIn("MISSING_OR_INVALID_LABELS", error_codes(missing))

        absent_value = copy.deepcopy(document)
        del absent_value["load"]["peak_force_n"]["value"]
        malformed_known = validate_artifact(absent_value)
        self.assertIn("KNOWN_VALUE_SHAPE", error_codes(malformed_known))

    def test_cross_fragment_pair_and_coordinate_unit_mismatch_fail_closed(self) -> None:
        wrong_pair = validate_artifact(load_case("crossfragmentwrongpair.json"))
        self.assertFalse(wrong_pair.valid)
        self.assertIn("SURFACE_FRAGMENT_PAIR_MISMATCH", error_codes(wrong_pair))

        wrong_units = validate_artifact(load_case("units_mismatch.json"))
        self.assertFalse(wrong_units.valid)
        self.assertIn("UNIT_MISMATCH", error_codes(wrong_units))

    def test_surface_pair_requires_reversed_triangles_and_measured_area(self) -> None:
        same_winding = load_case("complete_fixture.json")
        same_winding["fracture_surface_pairs"][0]["positive_surface"]["triangles"] = [[0, 2, 1]]
        winding_report = validate_artifact(same_winding)
        self.assertFalse(winding_report.valid)
        self.assertIn("PAIRED_SURFACE_TRIANGLES_NOT_RECIPROCAL", error_codes(winding_report))

        wrong_area = load_case("complete_fixture.json")
        wrong_area["fracture_surface_pairs"][0]["reference_area_m2"]["value"] *= 2.0
        area_report = validate_artifact(wrong_area)
        self.assertFalse(area_report.valid)
        self.assertIn("DECLARED_SURFACE_AREA_MISMATCH", error_codes(area_report))

    def test_nonfinite_json_and_gt_leak_are_rejected(self) -> None:
        raw_nonfinite = FIXTURES / "nonfinite.json"
        with self.assertRaises(ArtifactContractError) as caught:
            read_observation(raw_nonfinite)
        self.assertEqual(caught.exception.code, "STRICT_JSON_REJECTED")

        parsed_nonfinite = json.loads(raw_nonfinite.read_text(encoding="utf-8"))
        report = validate_artifact(parsed_nonfinite)
        self.assertFalse(report.valid)
        self.assertIn("NONFINITE_NUMBER", error_codes(report))

        leaked = validate_artifact(load_case("gtleak.json"))
        self.assertFalse(leaked.valid)
        self.assertIn("GT_LEAK_IN_OBSERVATION", error_codes(leaked))

        # Observation payloads use a tiny geometry-only PLY subset; metadata
        # comments are rejected before their bytes reach the inference view.
        leaking_mesh = FIXTURES / "payload" / "payload_metadata_leak.test.ply"
        leaking_artifact = FIXTURES / "payload_metadata_leak.test.json"
        leaking_mesh.write_bytes(
            b"ply\nformat ascii 1.0\ncomment mother_id hidden\nelement vertex 1\n"
            b"property float x\nproperty float y\nproperty float z\n"
            b"element face 0\nproperty list uchar int vertex_indices\nend_header\n0 0 0\n"
        )
        try:
            document = load_case("complete_fixture.json")
            reference = document["observation"]["fragment_observations"][0]["observed_mesh"]
            reference["path"] = "payload/payload_metadata_leak.test.ply"
            reference["sha256"] = hashlib.sha256(leaking_mesh.read_bytes()).hexdigest()
            leaking_artifact.write_text(json.dumps(document, allow_nan=False), encoding="utf-8")
            with self.assertRaises(ArtifactContractError) as caught:
                read_observation(leaking_artifact)
            self.assertEqual(caught.exception.code, "OBSERVATION_PLY_HEADER_UNSUPPORTED")
        finally:
            leaking_artifact.unlink(missing_ok=True)
            leaking_mesh.unlink(missing_ok=True)

    def test_graph_components_do_not_substitute_for_closed_fragments(self) -> None:
        document = load_case("complete_fixture.json")
        document["fragmentation_basis"] = "GRAPH_COMPONENTS_ONLY"
        report = validate_artifact(document)
        self.assertFalse(report.valid)
        self.assertIn("GRAPH_COMPONENTS_CANNOT_DEFINE_FRAGMENTS", error_codes(report))

        document = load_case("complete_fixture.json")
        del document["fragments"][0]["closed_mesh"]
        report = validate_artifact(document)
        self.assertFalse(report.valid)
        self.assertIn("MISSING_FIELD", error_codes(report))

    def test_physical_checkbox_cannot_promote_a_fixture(self) -> None:
        document = load_case("complete_fixture.json")
        document["qualification"]["physical_validation_passed"] = True
        document["qualification"]["class"] = "PHYSICALLY_VALIDATED"
        report = validate_artifact(document)
        self.assertFalse(report.valid)
        self.assertFalse(report.pretraining_qualified)
        self.assertIn("UNKNOWN_FIELD", error_codes(report))
        self.assertIn("PHYSICAL_RECEIPT_UNVERIFIED", error_codes(report))

    def test_inference_projection_contains_observations_only_and_gt_is_explicit_audit_path(self) -> None:
        source = FIXTURES / "complete_fixture.json"
        projected = read_observation(source)
        self.assertEqual(projected["schema"], "a12.fracture_observation.v1")
        self.assertEqual(set(projected), {"schema", "verified_payloads", "observation"})
        self.assertTrue(all(item["sha256_verified"] for item in projected["verified_payloads"]))
        self.assertTrue(all("content_base64" in item for item in projected["verified_payloads"]))
        flattened = json.dumps(projected).lower()
        for forbidden in ("mother_fixture", "clean_ground_truth", "event_fixture", "pair_fixture", "fragment_fixture"):
            self.assertNotIn(forbidden, flattened)

        with self.assertRaises(ArtifactContractError) as caught:
            read_clean_ground_truth_for_audit(source, purpose="inference")
        self.assertEqual(caught.exception.code, "GROUND_TRUTH_AUDIT_PURPOSE_REQUIRED")
        audit_truth = read_clean_ground_truth_for_audit(source, purpose="offline_audit")
        self.assertEqual(audit_truth["role"], "AUDIT_ONLY_CLEAN_GROUND_TRUTH")
        self.assertEqual(audit_truth["qualification_class"], "FIXTURE_GEOMETRY_ONLY")

    def test_export_and_read_roundtrip_uses_strict_json_and_only_writes_review_path(self) -> None:
        target = FIXTURES / "export_roundtrip.test.json"
        bad_payload_target = FIXTURES / "payload_hash_mismatch.test.json"
        try:
            result = export_artifact(load_case("complete_fixture.json"), target)
            self.assertEqual(result["validation"]["valid"], True)
            self.assertEqual(result["validation"]["pretraining_qualified"], False)
            digest = hashlib.sha256(target.read_bytes()).hexdigest()
            self.assertEqual(read_observation(target, expected_artifact_sha256=digest)["schema"], "a12.fracture_observation.v1")
            with self.assertRaises(ArtifactContractError) as caught:
                read_observation(target, expected_artifact_sha256="0" * 64)
            self.assertEqual(caught.exception.code, "ARTIFACT_SHA256_MISMATCH")

            bad_doc = load_case("complete_fixture.json")
            bad_doc["observation"]["fragment_observations"][0]["observed_mesh"]["sha256"] = "0" * 64
            bad_payload_target.write_text(json.dumps(bad_doc, allow_nan=False), encoding="utf-8")
            with self.assertRaises(ArtifactContractError) as caught:
                read_observation(bad_payload_target)
            self.assertEqual(caught.exception.code, "PAYLOAD_SHA256_MISMATCH")
        finally:
            target.unlink(missing_ok=True)
            target.with_name(target.name + ".tmp").unlink(missing_ok=True)
            bad_payload_target.unlink(missing_ok=True)


if __name__ == "__main__":
    unittest.main()
