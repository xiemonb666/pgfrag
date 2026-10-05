from __future__ import annotations

import copy
import hashlib
import json
import os
import shutil
import tempfile
import unittest
from pathlib import Path

from pgfrag.experiments.fracture_artifact_contract_20261006.mesh_binding_validator import (
    BindingValidationError,
    MAX_FILE_BYTES,
    validate_binding,
)


PACKAGE_DIR = Path(__file__).resolve().parent
FIXTURES_DIR = PACKAGE_DIR / "fixtures"
REVIEW_DIR = PACKAGE_DIR.parent.parent / "tmp" / "review_20261006" / "fracture_artifact_contract_20261006"
BASE_MANIFEST = FIXTURES_DIR / "planar_box_mesh_binding_manifest.json"
SOURCE_FIXTURE = (FIXTURES_DIR / "../../../tmp/review_20261006/explicit_fragment_fixture_20261006").resolve()


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _dump_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=False, allow_nan=False) + "\n", encoding="utf-8")


class MeshBindingValidatorTests(unittest.TestCase):
    def setUp(self) -> None:
        REVIEW_DIR.mkdir(parents=True, exist_ok=True)
        self.temp = tempfile.TemporaryDirectory(dir=REVIEW_DIR, prefix="mesh-binding-case-")
        self.case_root = Path(self.temp.name)
        self.fixture_root = self.case_root / "source_fixture"
        self.fixture_root.mkdir()
        for name in (
            "explicit_fragment_fixture_output.json",
            "interface_pair_ledger.json",
            "piece_negative.obj",
            "piece_positive.obj",
        ):
            shutil.copyfile(SOURCE_FIXTURE / name, self.fixture_root / name)
        self.manifest = copy.deepcopy(_read_json(BASE_MANIFEST))
        self.manifest["source_fixture_dir"] = Path(os.path.relpath(self.fixture_root, FIXTURES_DIR)).as_posix()
        self.manifest_path = FIXTURES_DIR / f".mesh-binding-test-{os.getpid()}-{id(self)}.json"

    def tearDown(self) -> None:
        try:
            self.manifest_path.unlink(missing_ok=True)
        finally:
            self.temp.cleanup()

    def _pin_json_inputs(self) -> None:
        for ref_key in ("output_file", "pair_ledger_file"):
            ref = self.manifest[ref_key]
            ref["sha256"] = hashlib.sha256((self.fixture_root / ref["file"]).read_bytes()).hexdigest()

    def _pin_mesh(self, side: str) -> None:
        fragment = next(item for item in self.manifest["fragments"] if item["side"] == side)
        fragment["mesh_sha256"] = hashlib.sha256((self.fixture_root / fragment["mesh_file"]).read_bytes()).hexdigest()

    def _run(self) -> dict:
        _dump_json(self.manifest_path, self.manifest)
        return validate_binding(self.manifest_path)

    def test_existing_fixture_binds_actual_obj_faces_and_metrics(self) -> None:
        report = validate_binding(BASE_MANIFEST)
        self.assertEqual(report["validation_status"], "PASS")
        self.assertEqual(report["qualification"], "FIXTURE_GEOMETRY_ONLY")
        self.assertEqual(report["mechanics_qualification"], "NOT_COUPLED")
        self.assertEqual(report["physical_pack_count"], 0)
        self.assertEqual(report["fragments"]["negative"]["signed_volume_m3"], 4.4399999999999995)
        self.assertEqual(report["fragments"]["positive"]["signed_volume_m3"], 3.56)
        self.assertAlmostEqual(report["interface_area_m2"], 4.354859354789773)
        self.assertEqual(report["partition_geometry_status"], "VOLUME_SUM_CONSERVED_ONLY_NO_OVERLAP_OR_FULL_COVERAGE_PROOF")
        self.assertTrue(report["producer_boolean_claims_ignored"]["output_checks_object_ignored"])

    def test_missing_pair_face_is_rejected(self) -> None:
        output_path = self.fixture_root / "explicit_fragment_fixture_output.json"
        ledger_path = self.fixture_root / "interface_pair_ledger.json"
        output = _read_json(output_path)
        ledger = _read_json(ledger_path)
        ledger["triangle_pairs"][0]["negative_face_index"] = 99
        ledger["triangle_pairs"][0]["negative_face_id"] = "fixture_piece_negative:f99"
        output["interface_pair_ledger"] = copy.deepcopy(ledger)
        _dump_json(ledger_path, ledger)
        _dump_json(output_path, output)
        self._pin_json_inputs()
        with self.assertRaisesRegex(BindingValidationError, "missing negative face 99"):
            self._run()

    def test_edited_bytes_fail_hash_before_geometry_parse(self) -> None:
        path = self.fixture_root / "piece_negative.obj"
        raw = path.read_bytes()
        path.write_bytes(raw.replace(b"0.69999999999999996", b"0.79999999999999996", 1))
        with self.assertRaisesRegex(BindingValidationError, "SHA-256 mismatch"):
            self._run()

    def test_coordinate_change_with_repinned_hash_fails_embedded_binding(self) -> None:
        path = self.fixture_root / "piece_negative.obj"
        raw = path.read_bytes()
        path.write_bytes(raw.replace(b"0.69999999999999996", b"0.79999999999999996", 1))
        self._pin_mesh("negative")
        with self.assertRaisesRegex(BindingValidationError, "differs from output vertices_m"):
            self._run()

    def test_missing_mesh_is_rejected(self) -> None:
        (self.fixture_root / "piece_positive.obj").unlink()
        with self.assertRaisesRegex(BindingValidationError, "missing positive mesh"):
            self._run()

    def test_closed_mesh_with_wrong_global_winding_is_rejected(self) -> None:
        path = self.fixture_root / "piece_negative.obj"
        lines = path.read_text(encoding="utf-8").splitlines()
        flipped_lines = []
        for line in lines:
            if line.startswith("f "):
                _, a, b, c = line.split()
                flipped_lines.append(f"f {a} {c} {b}")
            else:
                flipped_lines.append(line)
        path.write_text("\n".join(flipped_lines) + "\n", encoding="utf-8")
        output_path = self.fixture_root / "explicit_fragment_fixture_output.json"
        output = _read_json(output_path)
        negative = next(piece for piece in output["pieces"] if piece["side"] == "negative")
        negative["triangles"] = [[face[0], face[2], face[1]] for face in negative["triangles"]]
        _dump_json(output_path, output)
        self._pin_mesh("negative")
        self._pin_json_inputs()
        with self.assertRaisesRegex(BindingValidationError, "inward-oriented or has non-positive signed volume"):
            self._run()

    def test_unit_mismatch_is_rejected(self) -> None:
        self.manifest["mesh_coordinate_unit"] = "mm"
        with self.assertRaisesRegex(BindingValidationError, "must be exactly 'm'"):
            self._run()

    def test_mesh_byte_cap_is_checked_before_parsing(self) -> None:
        path = self.fixture_root / "piece_negative.obj"
        path.write_bytes(b"x" * (MAX_FILE_BYTES + 1))
        with self.assertRaisesRegex(BindingValidationError, "exceeds 262144 byte limit"):
            self._run()


if __name__ == "__main__":
    unittest.main()
