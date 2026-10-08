"""Narrow unit tests for the source-only bridge packager; intentionally not run here."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


MODULE_PATH = Path(__file__).with_name("bridge_contract.py")
SPEC = importlib.util.spec_from_file_location("pgfrag_bridge_contract_candidate", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
bridge = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(bridge)


def _pair_record(left: int, right: int, state: str, positions: dict[int, list[float]]) -> dict:
    return {
        "pair": [left, right],
        "classification": state,
        "admissible_static_candidate": state == "IN_MATERIAL",
        "raw_hits": [{"t": 0.5, "kind": "fixture"}],
        "open_intervals": [{"t0": 0.0, "t1": 1.0, "state": "fixture"}],
        "reciprocal": {
            "same_classification": True,
            "directed_records": [
                {"source_point_id": left, "target_point_id": right, "state": state},
                {"source_point_id": right, "target_point_id": left, "state": state},
            ],
        },
        "endpoints": [
            {"point_id": left, "position_mm": positions[left]},
            {"point_id": right, "position_mm": positions[right]},
        ],
        "reason": "synthetic unit-test row",
    }


class BridgeContractTests(unittest.TestCase):
    def make_inputs(self, root: Path, point_ids: list[int], pair_states: list[tuple[int, int, str]]):
        positions = {point_id: [float(index), 0.0, 0.0] for index, point_id in enumerate(point_ids)}
        payload = "point_id,x_mm,y_mm,z_mm\n" + "".join(
            f"{point_id},{positions[point_id][0]},0,0\n" for point_id in point_ids
        )
        point_path = root / "points.csv"
        point_path.write_text(payload, encoding="utf-8", newline="")
        cad_path = root / "source.step"
        cad_path.write_bytes(b"synthetic STEP bytes used only as a SHA fixture\n")
        source_path = root / "bond_visibility.py"
        source_path.write_text("# synthetic classifier source pin\n", encoding="utf-8")

        pairs = [_pair_record(i, j, state, positions) for i, j, state in pair_states]
        states = ("IN_MATERIAL", "OUT_OR_VOID", "UNKNOWN")
        counts = {state: sum(row["classification"] == state for row in pairs) for state in states}
        cad_sha = hashlib.sha256(cad_path.read_bytes()).hexdigest()
        result = {
            "schema": bridge.VISIBILITY_SCHEMA,
            "status": "STATIC_CLASSIFICATION_COMPLETE",
            "algorithm": "unit-test-classifier",
            "source": {
                "source_origin": "CANONICAL_STEP_FILE",
                "status": "READY",
                "sha256_status": "MATCHED_ACTUAL_FILE_BYTES",
                "source_geometry_sha256": cad_sha,
                "expected_sha256": cad_sha,
                "source_length_unit": "mm",
                "frame": "fixture-reference",
                "source_accuracy_mm": 0.001,
                "source_accuracy_basis": "declared model tolerance; not measured physical accuracy",
            },
            "reference_coordinates": {
                "coordinate_unit": "mm",
                "frame": "fixture-reference",
                "point_count": len(point_ids),
            },
            "initial_topology_ledger": {
                "classification_counts": counts,
                "candidate_undirected_pair_count": len(pairs),
                "initial_material_bond_candidate_count": counts["IN_MATERIAL"],
                "geometry_excluded_pair_count": counts["OUT_OR_VOID"],
                "unknown_pair_count": counts["UNKNOWN"],
            },
            "pairs": pairs,
        }
        result_path = root / "classifier-result.json"
        result_path.write_text(json.dumps(result, separators=(",", ":")), encoding="utf-8")
        selection = {"test_nodeid": "test_bridge_contract.py::captured_result", "call_index": 0}
        sidecar = {
            "schema": "pgfrag.bond-visibility-geometry-sidecar.v1",
            "module_file": str(source_path.resolve()),
            "algorithm": result["algorithm"],
            "pytest_exit_code": 0,
            "captured_call_count": 1,
            "results": [{**selection, "result": result}],
        }
        sidecar_path = root / "geometry-sidecar.json"
        sidecar_path.write_text(json.dumps(sidecar, separators=(",", ":")), encoding="utf-8")
        args = {
            "classifier_result_path": result_path,
            "point_payload_path": point_path,
            "source_cad_path": cad_path,
            "raw_hit_interval_sidecar_path": sidecar_path,
            "classifier_source_path": source_path,
            "sidecar_selection": selection,
            "horizon_m": 0.01,
            "reference_frame": "fixture-reference",
            "native_frame": "fixture-native",
            "native_frame_from_reference_frame_m": [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
        }
        return args, result, sidecar, sidecar_path, source_path

    def test_source_sha_is_computed_from_the_actual_source_file(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            args, _, _, _, source_path = self.make_inputs(root, [0, 4], [(0, 4, "IN_MATERIAL")])
            contract = bridge.build_contract(**args)
            actual_sha = hashlib.sha256(source_path.read_bytes()).hexdigest()
            self.assertEqual(contract["evidence"]["classifier_source_sha256"], actual_sha)
            self.assertEqual(contract["source"]["classifier_source_sha256"], actual_sha)
            self.assertEqual(contract["evidence"]["classifier_source_path"], str(source_path.resolve()))
            self.assertNotIn("preflight_status", contract)
            self.assertEqual(
                contract["initial_physics_events"]["evidence_status"],
                "PLANNED_NOT_OBSERVED_NATIVE_NOT_RUN",
            )

    def test_sidecar_raw_hit_mismatch_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            args, _, sidecar, sidecar_path, _ = self.make_inputs(
                Path(temp), [0, 4], [(0, 4, "IN_MATERIAL")]
            )
            sidecar["results"][0]["result"]["pairs"][0]["raw_hits"] = []
            sidecar_path.write_text(json.dumps(sidecar), encoding="utf-8")
            with self.assertRaises(bridge.ContractError):
                bridge.build_contract(**args)

    def test_sidecar_open_interval_mismatch_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            args, _, sidecar, sidecar_path, _ = self.make_inputs(
                Path(temp), [0, 4], [(0, 4, "IN_MATERIAL")]
            )
            sidecar["results"][0]["result"]["pairs"][0]["open_intervals"] = []
            sidecar_path.write_text(json.dumps(sidecar), encoding="utf-8")
            with self.assertRaises(bridge.ContractError):
                bridge.build_contract(**args)

    def test_unknown_is_retained_and_blocks_the_whole_native_asset(self):
        with tempfile.TemporaryDirectory() as temp:
            args, _, _, _, _ = self.make_inputs(Path(temp), [0, 4], [(0, 4, "UNKNOWN")])
            contract = bridge.build_contract(**args)
            self.assertEqual(contract["static_contract_status"], "VALID")
            self.assertEqual(contract["native_mask_eligibility"], "BLOCKED_UNKNOWN")
            self.assertEqual(contract["mask"]["unknown_pairs"], [[0, 4]])
            self.assertFalse(contract["native_gate"]["initialization_authorized"])

    def test_empty_candidate_empty_mask_and_zero_degree_are_blocked(self):
        cases = [
            ([0, 4], [], "BLOCKED_EMPTY_CANDIDATE_SET", [0, 4]),
            ([0, 4], [(0, 4, "OUT_OR_VOID")], "BLOCKED_EMPTY_MASK", [0, 4]),
            (
                [0, 4, 9],
                [(0, 4, "IN_MATERIAL"), (0, 9, "OUT_OR_VOID"), (4, 9, "OUT_OR_VOID")],
                "BLOCKED_ZERO_DEGREE",
                [9],
            ),
        ]
        for point_ids, pair_states, expected, zero_degree in cases:
            with self.subTest(expected=expected), tempfile.TemporaryDirectory() as temp:
                args, _, _, _, _ = self.make_inputs(Path(temp), point_ids, pair_states)
                contract = bridge.build_contract(**args)
                self.assertEqual(contract["native_mask_eligibility"], expected)
                self.assertEqual(contract["native_gate"]["zero_degree_point_ids"], zero_degree)
                self.assertFalse(contract["native_gate"]["initialization_authorized"])

    def test_invalid_accuracy_basis_does_not_enter_contract(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            args, result, sidecar, _, _ = self.make_inputs(root, [0, 4], [(0, 4, "IN_MATERIAL")])
            result["source"]["source_accuracy_basis"] = ["untrusted", "non-string"]
            args["classifier_result_path"].write_text(json.dumps(result), encoding="utf-8")
            sidecar["results"][0]["result"] = result
            args["raw_hit_interval_sidecar_path"].write_text(json.dumps(sidecar), encoding="utf-8")
            with self.assertRaises(bridge.ContractError):
                bridge.build_contract(**args)

    def test_per_file_and_aggregate_budgets_reject_before_any_read(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "too-large.txt"
            path.write_bytes(b"12345")
            with patch.object(bridge, "_read_once_bounded") as bounded_read:
                with self.assertRaises(bridge.ContractError):
                    bridge._read_stable_bundle((("fixture", path, 4),))
                bounded_read.assert_not_called()

        with tempfile.TemporaryDirectory() as temp:
            first = Path(temp) / "first.txt"
            second = Path(temp) / "second.txt"
            first.write_bytes(b"123")
            second.write_bytes(b"456")
            old_total = bridge.MAX_TOTAL_INPUT_BYTES
            try:
                bridge.MAX_TOTAL_INPUT_BYTES = 5
                with patch.object(bridge, "_read_once_bounded") as bounded_read:
                    with self.assertRaises(bridge.ContractError):
                        bridge._read_stable_bundle((("first", first, 4), ("second", second, 4)))
                    bounded_read.assert_not_called()
            finally:
                bridge.MAX_TOTAL_INPUT_BYTES = old_total

    def test_point_and_pair_caps_are_enforced(self):
        original_points = bridge.MAX_POINTS
        original_pairs = bridge.MAX_PAIRS
        try:
            bridge.MAX_POINTS = 1
            with self.assertRaises(bridge.ContractError):
                bridge._point_rows(b"point_id,x_mm,y_mm,z_mm\n0,0,0,0\n1,1,0,0\n")
            bridge.MAX_PAIRS = 0
            with self.assertRaises(bridge.ContractError):
                bridge._classifier_pairs({"pairs": [{}]}, [])
        finally:
            bridge.MAX_POINTS = original_points
            bridge.MAX_PAIRS = original_pairs

    def test_sidecar_selection_and_strict_duplicate_key_checks(self):
        with tempfile.TemporaryDirectory() as temp:
            args, _, sidecar, sidecar_path, _ = self.make_inputs(
                Path(temp), [0, 4], [(0, 4, "IN_MATERIAL")]
            )
            sidecar["results"].append(dict(sidecar["results"][0]))
            sidecar["captured_call_count"] = 2
            sidecar_path.write_text(json.dumps(sidecar), encoding="utf-8")
            with self.assertRaises(bridge.ContractError):
                bridge.build_contract(**args)
            sidecar_path.write_text(
                '{"schema":"pgfrag.bond-visibility-geometry-sidecar.v1",'
                '"schema":"pgfrag.bond-visibility-geometry-sidecar.v1"}',
                encoding="utf-8",
            )
            with self.assertRaises(bridge.ContractError):
                bridge.build_contract(**args)


if __name__ == "__main__":
    unittest.main()
