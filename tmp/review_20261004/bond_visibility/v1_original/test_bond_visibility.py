"""Static CAD visibility classifier controls; intentionally not solver tests.

These fixtures create small in-memory OCC solids only when the test suite is
run. The current review stage records the fixture plan but does not execute
pytest, CAD operations, or STEP export.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
import tempfile
import unittest


try:
    import cadquery as cq
except ImportError:  # Keep non-CAD environments able to collect the suite.
    cq = None

from pgfrag.bond_visibility import (
    _source_from_shape,
    classify_candidate_pairs,
    load_canonical_step,
)


FRAME = "right_handed_z_up"
DISTANCE_TOLERANCE_MM = 1e-4


def _source(shape, *, domains=None, accuracy_mm=0.005):
    return _source_from_shape(
        shape,
        source_length_unit="mm",
        frame=FRAME,
        source_accuracy_mm=accuracy_mm,
        source_accuracy_basis=(
            "analytic fixture dimensions and OCC construction; test-only provenance"
            if accuracy_mm is not None
            else None
        ),
        material_domain_ids=domains,
        backend_versions={"cadquery": "test-fixture", "ocp": "test-fixture"},
    )


def _classify(source, points, domains, pairs, *, point_ids=None, tolerance=DISTANCE_TOLERANCE_MM):
    ids = list(range(len(points))) if point_ids is None else list(point_ids)
    return classify_candidate_pairs(
        source,
        reference_point_ids=ids,
        reference_positions_mm=points,
        point_material_domain_ids=domains,
        candidate_pairs=pairs,
        reference_frame=FRAME,
        distance_tolerance_mm=tolerance,
    )


@unittest.skipIf(cq is None, "CadQuery/OCP runtime is required for geometric fixtures")
class BondVisibilityGeometryTests(unittest.TestCase):
    def test_hollow_shell_endpoints_inside_but_segment_crosses_void(self) -> None:
        outer = cq.Workplane("XY").box(20.0, 20.0, 20.0)
        cavity = cq.Workplane("XY").box(12.0, 12.0, 12.0)
        shell = outer.cut(cavity).val().wrapped
        result = _classify(
            _source(shell, domains=("body",)),
            [(-8.0, 0.0, 0.0), (8.0, 0.0, 0.0)],
            ["body", "body"],
            [(0, 1)],
        )
        pair = result["pairs"][0]
        self.assertEqual(pair["endpoints"][0]["state"], "IN_MATERIAL")
        self.assertEqual(pair["endpoints"][1]["state"], "IN_MATERIAL")
        self.assertEqual(pair["classification"], "OUT_OR_VOID")
        self.assertEqual(pair["reason"], "AIR_PATH")
        self.assertEqual(
            [row["occupancy"] for row in pair["open_intervals"]],
            ["IN_MATERIAL", "OUT_OR_VOID", "IN_MATERIAL"],
        )
        self.assertEqual(len(pair["raw_hits"]), 2)

    def test_hollow_shell_preserves_all_outer_and_inner_hits(self) -> None:
        shell = (
            cq.Workplane("XY").box(20.0, 20.0, 20.0)
            .cut(cq.Workplane("XY").box(12.0, 12.0, 12.0))
            .val().wrapped
        )
        pair = _classify(
            _source(shell, domains=("body",)),
            [(-11.0, 0.0, 0.0), (11.0, 0.0, 0.0)],
            ["body", "body"],
            [(0, 1)],
        )["pairs"][0]
        self.assertEqual(pair["classification"], "OUT_OR_VOID")
        self.assertEqual(len(pair["raw_hits"]), 4)
        values = [hit["parameter_mm_from_canonical_pair_start"] for hit in pair["raw_hits"]]
        for actual, expected in zip(values, (1.0, 5.0, 17.0, 21.0)):
            self.assertAlmostEqual(actual, expected, places=5)
        self.assertEqual(
            [row["occupancy"] for row in pair["open_intervals"]],
            ["OUT_OR_VOID", "IN_MATERIAL", "OUT_OR_VOID", "IN_MATERIAL", "OUT_OR_VOID"],
        )

    def test_same_wall_and_solid_coupon_remain_inside_material(self) -> None:
        shell = (
            cq.Workplane("XY").box(20.0, 20.0, 20.0)
            .cut(cq.Workplane("XY").box(12.0, 12.0, 12.0))
            .val().wrapped
        )
        shell_result = _classify(
            _source(shell, domains=("body",)),
            [(-8.0, 0.0, 0.0), (-7.0, 0.0, 0.0)],
            ["body", "body"],
            [(0, 1)],
        )
        self.assertEqual(shell_result["pairs"][0]["classification"], "IN_MATERIAL")

        coupon = cq.Workplane("XY").box(10.0, 10.0, 10.0).val().wrapped
        coupon_result = _classify(
            _source(coupon, domains=("coupon",)),
            [(-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)],
            ["coupon", "coupon"],
            [(0, 1)],
        )
        self.assertEqual(coupon_result["pairs"][0]["classification"], "IN_MATERIAL")
        self.assertEqual(coupon_result["initial_topology_ledger"]["whole_asset_physics_status"], "NOT_ASSESSED")

    def test_separate_handle_arms_are_distinct_domains_across_gap(self) -> None:
        left = cq.Workplane("XY").box(4.0, 2.0, 2.0).translate((-3.0, 0.0, 0.0)).val()
        right = cq.Workplane("XY").box(4.0, 2.0, 2.0).translate((3.0, 0.0, 0.0)).val()
        compound = cq.Compound.makeCompound([left, right]).wrapped
        source = _source(compound, domains=("left_arm", "right_arm"))
        pair = _classify(
            source,
            [(-3.0, 0.0, 0.0), (3.0, 0.0, 0.0)],
            ["left_arm", "right_arm"],
            [(0, 1)],
        )["pairs"][0]
        self.assertEqual(pair["classification"], "OUT_OR_VOID")
        self.assertEqual(pair["reason"], "DIFFERENT_MATERIAL_DOMAIN")
        self.assertIn("OUT_OR_VOID", [row["occupancy"] for row in pair["open_intervals"]])

    def test_tangent_and_on_boundary_paths_are_unknown(self) -> None:
        sphere = cq.Workplane("XY").sphere(5.0).val().wrapped
        tangent = _classify(
            _source(sphere, domains=("sphere",)),
            [(-6.0, 5.0, 0.0), (6.0, 5.0, 0.0)],
            ["sphere", "sphere"],
            [(0, 1)],
        )["pairs"][0]
        self.assertEqual(tangent["classification"], "UNKNOWN")

        shell = (
            cq.Workplane("XY").box(20.0, 20.0, 20.0)
            .cut(cq.Workplane("XY").box(12.0, 12.0, 12.0))
            .val().wrapped
        )
        on_face = _classify(
            _source(shell, domains=("body",)),
            [(-8.0, 6.0, 0.0), (-7.0, 6.0, 0.0)],
            ["body", "body"],
            [(0, 1)],
        )["pairs"][0]
        self.assertEqual(on_face["classification"], "UNKNOWN")

    def test_thin_material_interval_within_error_band_is_unknown(self) -> None:
        thin_solid = cq.Workplane("XY").box(0.01, 4.0, 4.0).val().wrapped
        pair = _classify(
            _source(thin_solid, domains=("thin",), accuracy_mm=0.01),
            [(-1.0, 0.0, 0.0), (1.0, 0.0, 0.0)],
            ["thin", "thin"],
            [(0, 1)],
            tolerance=0.001,
        )["pairs"][0]
        self.assertEqual(pair["classification"], "UNKNOWN")
        self.assertIn(
            "THIN_INTERVAL_WITHIN_ERROR_BAND",
            [interval["reason"] for interval in pair["open_intervals"]],
        )

    def test_non_solid_or_precision_unknown_source_cannot_make_material_bonds(self) -> None:
        wire = cq.Workplane("XY").rect(4.0, 4.0).val().wrapped
        wire_source = _source(wire)
        wire_result = _classify(
            wire_source,
            [(-1.0, 0.0, 0.0), (1.0, 0.0, 0.0)],
            ["unknown", "unknown"],
            [(0, 1)],
        )
        self.assertEqual(wire_result["pairs"][0]["classification"], "UNKNOWN")
        self.assertEqual(wire_result["initial_topology_ledger"]["whole_asset_physics_status"], "NOT_ASSESSED")

        coupon = cq.Workplane("XY").box(10.0, 10.0, 10.0).val().wrapped
        precision_source = _source(coupon, accuracy_mm=None)
        precision_result = _classify(
            precision_source,
            [(-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)],
            ["solid-0", "solid-0"],
            [(0, 1)],
        )
        self.assertEqual(precision_result["pairs"][0]["classification"], "UNKNOWN")
        self.assertIsNone(precision_result["tolerances"]["combined_occupancy_uncertainty_band_mm"])

    def test_directed_candidate_input_is_canonical_and_reciprocal(self) -> None:
        coupon = cq.Workplane("XY").box(10.0, 10.0, 10.0).val().wrapped
        result = _classify(
            _source(coupon),
            [(-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)],
            ["solid-0", "solid-0"],
            [(20, 10), (10, 20)],
            point_ids=(20, 10),
        )
        self.assertEqual(result["initial_topology_ledger"]["candidate_directed_input_count"], 2)
        self.assertEqual(result["initial_topology_ledger"]["candidate_undirected_pair_count"], 1)
        pair = result["pairs"][0]
        self.assertEqual(pair["pair"], [10, 20])
        self.assertEqual(pair["classification"], "IN_MATERIAL")
        self.assertTrue(pair["reciprocal"]["same_classification"])
        self.assertEqual(len(pair["reciprocal"]["directed_records"]), 2)


class CanonicalStepBindingTests(unittest.TestCase):
    def test_actual_file_hash_mismatch_is_rejected_before_step_import(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "canonical.step"
            payload = b"not a STEP file; hash mismatch must reject before CAD import\n"
            path.write_bytes(payload)
            actual = hashlib.sha256(payload).hexdigest()
            with self.assertRaisesRegex(ValueError, "canonical STEP SHA-256 mismatch"):
                load_canonical_step(
                    path,
                    expected_sha256="0" * 64 if actual != "0" * 64 else "1" * 64,
                    source_length_unit="mm",
                    frame=FRAME,
                    source_accuracy_mm=0.01,
                    source_accuracy_basis="test manifest accuracy declaration",
                )


if __name__ == "__main__":
    unittest.main()
