from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from pgfrag.assets import validate_asset
from pgfrag.errors import ValidationError
from .helpers import make_asset, temporary_directory


class ForgeAssetContractTests(unittest.TestCase):
    def test_contract_fixture_passes_and_does_not_claim_geometry_recomputed(self) -> None:
        with temporary_directory() as temporary:
            root = make_asset(Path(temporary) / "asset")
            result = validate_asset(root)
        self.assertEqual(result["contract_status"], "PASS")
        self.assertTrue(result["declared_geometry_validation"]["report_status_matches_manifest"])
        self.assertFalse(result["mechanical_solve_supported_from_surface_asset"])
        self.assertIn("not a hard process memory limit", result["npz_payload_limit_note"])

    def test_vessel_forge_four_bands_conformance_fixture(self) -> None:
        candidate = Path(__file__).resolve().parents[2] / "relicweave-synlab" / "tmp" / "review_20261003" / "implementation" / "vf" / "artifacts" / "four_bands"
        if not candidate.is_dir():
            self.skipTest("read-only VesselForge conformance sample is unavailable")
        result = validate_asset(candidate)
        self.assertEqual(result["schema"], "forge.asset.v1")
        self.assertEqual(result["surface"]["faces"], 38480)
        self.assertEqual(result["declared_geometry_validation"]["status"], "FAIL")

    def test_manifest_self_hash_is_rejected(self) -> None:
        with temporary_directory() as temporary:
            root = make_asset(Path(temporary) / "asset")
            path = root / "manifest.json"
            manifest = json.loads(path.read_text(encoding="utf-8"))
            manifest["files"]["manifest.json"] = {"bytes": 1, "sha256": "0" * 64}
            path.write_text(json.dumps(manifest), encoding="utf-8")
            with self.assertRaisesRegex(ValidationError, "exclude manifest.json"):
                validate_asset(root)

    def test_material_binding_must_match(self) -> None:
        with temporary_directory() as temporary:
            root = make_asset(Path(temporary) / "asset", binding_id="other-material")
            with self.assertRaisesRegex(ValidationError, "binding does not match"):
                validate_asset(root)

    def test_material_coordinates_are_required_and_equal_to_rest_vertices(self) -> None:
        with temporary_directory() as temporary:
            root = make_asset(Path(temporary) / "asset", material_xyz_offset=1.0)
            with self.assertRaisesRegex(ValidationError, "exactly equal"):
                validate_asset(root)

    def test_asset_poisson_domain_follows_contract_and_is_not_solver_restriction(self) -> None:
        with temporary_directory() as temporary:
            root = make_asset(Path(temporary) / "asset")
            manifest_path = root / "manifest.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["material"]["mechanics"]["nu"] = -0.2
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            self.assertEqual(validate_asset(root)["contract_status"], "PASS")


if __name__ == "__main__":
    unittest.main()
