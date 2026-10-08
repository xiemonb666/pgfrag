"""Static CAD visibility classifier controls; intentionally not solver tests.

These fixtures create small in-memory OCC solids only when the test suite is
run. The current review stage records the fixture plan but does not execute
pytest, CAD operations, or STEP export.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch


try:
    import cadquery as cq
except ImportError:  # Keep non-CAD environments able to collect the suite.
    cq = None

from pgfrag.bond_visibility import (
    CanonicalStepSource,
    MaterialSolid,
    _cluster_hits,
    _combined_uncertainty_band,
    _safe_midpoint,
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


class BondVisibilityPureMathRegressionTests(unittest.TestCase):
    @staticmethod
    def _hit(parameter_mm):
        return {
            "solid_index": 0,
            "material_domain_id": "solid-0",
            "raw_hit_index": 1,
            "face_index_zero_based": 4,
            "parameter_mm_from_canonical_pair_start": float(parameter_mm),
            "point_mm": None,
            "point_m": None,
            "point_status": "SYNTHETIC_MATH_CONTROL",
            "face_identity_status": "MAPPED",
            "cluster_index": None,
        }

    def test_single_link_chain_with_span_above_tolerance_is_unresolved(self) -> None:
        tolerance = 1.0
        hits = [self._hit(value) for value in (0.0, 0.75 * tolerance, 1.5 * tolerance)]
        events, reasons = _cluster_hits(hits, tolerance)
        self.assertEqual(len(events), 1)
        self.assertEqual(events[0]["distinct_parameter_count"], 3)
        self.assertEqual(events[0]["parameter_span_mm"], 1.5 * tolerance)
        self.assertTrue(events[0]["unresolved"])
        self.assertTrue(events[0]["diagnostic_cut_only"])
        self.assertIn("UNRESOLVED_NEARBY_BOUNDARY_HITS", reasons)
        self.assertEqual(
            [hit["parameter_mm_from_canonical_pair_start"] for hit in hits],
            [0.0, 0.75 * tolerance, 1.5 * tolerance],
        )
        json.dumps({"hits": hits, "events": events, "reasons": reasons}, allow_nan=False)

    def test_near_distinct_hits_retain_short_void_as_unresolved(self) -> None:
        tolerance = 1.0
        cavity_width = 0.25 * tolerance
        hits = [self._hit(10.0), self._hit(10.0 + cavity_width)]
        events, reasons = _cluster_hits(hits, tolerance)
        self.assertEqual(len(events), 1)
        self.assertEqual(events[0]["raw_hit_count"], 2)
        self.assertEqual(events[0]["distinct_parameter_count"], 2)
        self.assertEqual(events[0]["parameter_span_mm"], cavity_width)
        self.assertEqual(events[0]["resolution_status"], "UNRESOLVED_DISTINCT_PARAMETERS")
        self.assertIn("UNRESOLVED_NEARBY_BOUNDARY_HITS", reasons)
        self.assertEqual(len(hits), 2)
        json.dumps({"hits": hits, "events": events, "reasons": reasons}, allow_nan=False)

    def test_unresolved_chain_forces_unknown_even_when_diagnostic_intervals_are_inside(self) -> None:
        source = CanonicalStepSource(
            source_path=None,
            source_origin="IN_MEMORY_TEST_FIXTURE",
            source_geometry_sha256=None,
            expected_sha256=None,
            source_length_unit="mm",
            frame=FRAME,
            source_accuracy_mm=0.01,
            source_accuracy_basis="pure math state-control fixture",
            shape=object(),
            solids=(MaterialSolid(0, "solid-0", object(), (object(),)),),
            status="READY",
            unknown_reasons=(),
            backend_versions={"ocp": "not-called"},
        )
        for parameters in ((10.0, 10.75, 11.5), (10.0, 10.25)):
            with self.subTest(parameters=parameters):
                synthetic_hits = [self._hit(value) for value in parameters]
                with patch(
                    "pgfrag.bond_visibility._collect_raw_hits",
                    return_value=(synthetic_hits, []),
                ), patch(
                    "pgfrag.bond_visibility._point_classifier_state",
                    return_value="IN",
                ):
                    result = classify_candidate_pairs(
                        source,
                        reference_point_ids=[10, 20],
                        reference_positions_mm=[(0.0, 0.0, 0.0), (20.0, 0.0, 0.0)],
                        point_material_domain_ids=["solid-0", "solid-0"],
                        candidate_pairs=[(10, 20)],
                        reference_frame=FRAME,
                        distance_tolerance_mm=1.0,
                    )
                pair = result["pairs"][0]
                self.assertTrue(all(row["occupancy"] == "IN_MATERIAL" for row in pair["open_intervals"]))
                self.assertEqual(pair["classification"], "UNKNOWN")
                self.assertFalse(pair["admissible_static_candidate"])
                self.assertTrue(pair["diagnostic_intervals_only"])
                self.assertTrue(all(row["diagnostic_only"] for row in pair["open_intervals"]))
                self.assertEqual(len(pair["raw_hits"]), len(parameters))
                self.assertEqual(pair["reason"], "UNRESOLVED_NEARBY_BOUNDARY_HITS")
                json.dumps(result, allow_nan=False)

    def test_numeric_range_failures_are_unknown_and_json_finite(self) -> None:
        band, reason = _combined_uncertainty_band(1.0e308, 1.0e308)
        self.assertIsNone(band)
        self.assertEqual(reason, "NUMERIC_RANGE_OVERFLOW_UNCERTAINTY_BAND")
        doubled_band, doubled_reason = _combined_uncertainty_band(1.0e308, 1.0)
        self.assertIsNone(doubled_band)
        self.assertEqual(doubled_reason, "NUMERIC_RANGE_OVERFLOW_UNCERTAINTY_BAND")

        source = CanonicalStepSource(
            source_path=None,
            source_origin="IN_MEMORY_TEST_FIXTURE",
            source_geometry_sha256=None,
            expected_sha256=None,
            source_length_unit="mm",
            frame=FRAME,
            source_accuracy_mm=0.01,
            source_accuracy_basis="pure math numeric-range control",
            shape=None,
            solids=(),
            status="UNKNOWN",
            unknown_reasons=("SOURCE_NOT_AVAILABLE",),
            backend_versions={},
        )
        result = classify_candidate_pairs(
            source,
            reference_point_ids=[1, 2],
            reference_positions_mm=[(-1.7e308, 0.0, 0.0), (1.7e308, 0.0, 0.0)],
            point_material_domain_ids=[None, None],
            candidate_pairs=[(1, 2)],
            reference_frame=FRAME,
            distance_tolerance_mm=1.0e-4,
        )
        pair = result["pairs"][0]
        self.assertEqual(pair["classification"], "UNKNOWN")
        self.assertEqual(pair["reason"], "NUMERIC_RANGE_OVERFLOW_SEGMENT_VECTOR")
        self.assertIsNone(pair["segment_length_mm"])
        self.assertEqual(pair["endpoints"][0]["position_mm"][0], -1.7e308)
        self.assertEqual(pair["endpoints"][1]["position_mm"][0], 1.7e308)
        encoded = json.dumps(result, allow_nan=False)
        self.assertNotIn("Infinity", encoded)
        self.assertNotIn("NaN", encoded)

        self.assertTrue(math.isfinite(_safe_midpoint(1.6e308, 1.7e308)))
        large_hits = [self._hit(1.6e308), self._hit(1.7e308)]
        large_events, large_reasons = _cluster_hits(large_hits, 1.1e307)
        self.assertEqual(len(large_events), 1)
        self.assertTrue(large_events[0]["unresolved"])
        self.assertAlmostEqual(large_events[0]["parameter_mm"] / 1.0e308, 1.65)
        self.assertEqual(large_reasons, ["UNRESOLVED_NEARBY_BOUNDARY_HITS"])
        json.dumps({"hits": large_hits, "events": large_events}, allow_nan=False)

    def test_overflowing_uncertainty_band_keeps_raw_hits_without_classifier_calls(self) -> None:
        source = CanonicalStepSource(
            source_path=None,
            source_origin="IN_MEMORY_TEST_FIXTURE",
            source_geometry_sha256=None,
            expected_sha256=None,
            source_length_unit="mm",
            frame=FRAME,
            source_accuracy_mm=1.0e308,
            source_accuracy_basis="pure math numeric-range control",
            shape=object(),
            solids=(MaterialSolid(0, "solid-0", object(), (object(),)),),
            status="READY",
            unknown_reasons=(),
            backend_versions={"ocp": "not-called"},
        )
        raw_hit = self._hit(1.0)
        with patch(
            "pgfrag.bond_visibility._collect_raw_hits",
            return_value=([raw_hit], []),
        ), patch(
            "pgfrag.bond_visibility._point_classifier_state",
            return_value="IN",
        ) as classify_point:
            result = classify_candidate_pairs(
                source,
                reference_point_ids=[1, 2],
                reference_positions_mm=[(0.0, 0.0, 0.0), (2.0, 0.0, 0.0)],
                point_material_domain_ids=["solid-0", "solid-0"],
                candidate_pairs=[(1, 2)],
                reference_frame=FRAME,
                distance_tolerance_mm=1.0e308,
            )
        classify_point.assert_not_called()
        pair = result["pairs"][0]
        self.assertEqual(pair["classification"], "UNKNOWN")
        self.assertEqual(pair["reason"], "NUMERIC_RANGE_OVERFLOW_UNCERTAINTY_BAND")
        self.assertEqual(len(pair["raw_hits"]), 1)
        self.assertTrue(pair["diagnostic_intervals_only"])
        self.assertTrue(all(interval["diagnostic_only"] for interval in pair["open_intervals"]))
        self.assertIsNone(result["tolerances"]["combined_occupancy_uncertainty_band_mm"])
        self.assertEqual(
            result["tolerances"]["numeric_range_status"],
            "NUMERIC_RANGE_OVERFLOW_UNCERTAINTY_BAND",
        )
        json.dumps(result, allow_nan=False)


if __name__ == "__main__":
    unittest.main()
