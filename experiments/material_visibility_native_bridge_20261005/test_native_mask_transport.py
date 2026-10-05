"""Focused transport tests; source only in this handoff and intentionally not run."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import tempfile
import tomllib
import unittest
from pathlib import Path
from unittest.mock import patch


MODULE_PATH = Path(__file__).with_name("native_mask_transport.py")
SPEC = importlib.util.spec_from_file_location("pgfrag_native_mask_transport_candidate", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
transport = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(transport)


def _point_hash(point_ids: list[int]) -> str:
    payload = b"pgfrag.point-ids.v1\n" + b"".join(f"{point_id}\n".encode("ascii") for point_id in point_ids)
    return hashlib.sha256(payload).hexdigest()


def _pair_hash(rows: list[dict]) -> str:
    payload = bytearray(b"pgfrag.native-pair-map.v1\n")
    for row in rows:
        left, right = row["pair"]
        payload.extend(f"{left},{right},{row['classification']}\n".encode("ascii"))
    return hashlib.sha256(payload).hexdigest()


def _mask_hash(pairs: list[list[int]]) -> str:
    payload = bytearray(b"pgfrag.native-included-mask.v1\n")
    for left, right in pairs:
        payload.extend(f"{left},{right}\n".encode("ascii"))
    return hashlib.sha256(payload).hexdigest()


def _valid_contract(root: Path, point_ids: list[int] | None = None, states: list[tuple[int, int, str]] | None = None) -> dict:
    point_ids = point_ids or [0, 8]
    states = states if states is not None else [(0, 8, "IN_MATERIAL")]
    positions = {point_id: [float(index), float(index) * 0.25, 0.0] for index, point_id in enumerate(point_ids)}
    rows = [
        {"point_id": point_id, "native_index": index, "position_mm": positions[point_id]}
        for index, point_id in enumerate(point_ids, start=1)
    ]
    decisions = [
        {"pair": [left, right], "classification": state, "reason": f"fixture {left}-{right} {state}"}
        for left, right, state in states
    ]
    decisions.sort(key=lambda row: tuple(row["pair"]))
    included = [row["pair"] for row in decisions if row["classification"] == "IN_MATERIAL"]
    excluded = [row["pair"] for row in decisions if row["classification"] == "OUT_OR_VOID"]
    unknown = [row["pair"] for row in decisions if row["classification"] == "UNKNOWN"]
    point_degrees = {point_id: 0 for point_id in point_ids}
    for left, right in included:
        point_degrees[left] += 1
        point_degrees[right] += 1
    zero_degrees = [point_id for point_id in point_ids if point_degrees[point_id] == 0]
    count_by_state = {state: sum(row["classification"] == state for row in decisions) for state in transport.STATES}
    if count_by_state["UNKNOWN"]:
        eligibility = "BLOCKED_UNKNOWN"
    elif not decisions:
        eligibility = "BLOCKED_EMPTY_CANDIDATE_SET"
    elif not included:
        eligibility = "BLOCKED_EMPTY_MASK"
    elif zero_degrees:
        eligibility = "BLOCKED_ZERO_DEGREE"
    else:
        eligibility = "CANDIDATE_MASK_POSITIVE_DEGREE_NATIVE_UNVERIFIED"
    digest = lambda ch: ch * 64
    return {
        "schema": transport.AUDIT_SCHEMA,
        "static_contract_status": "VALID",
        "native_mask_eligibility": eligibility,
        "native_gate": {
            "status": "UNVERIFIED",
            "raw_candidate_set_reproduced": False,
            "native_graph_readback": "NOT_RUN",
            "initialization_authorized": False,
            "candidate_pair_count_nonzero": bool(decisions),
            "included_pair_count_nonzero": bool(included),
            "zero_degree_point_ids": zero_degrees,
        },
        "source": {
            "source_cad_path": str(root / "geometry.step"),
            "source_cad_sha256": digest("a"),
            "classifier_source_sha256": digest("b"),
            "source_hash_status": "MATCHED_ACTUAL_FILE_BYTES",
            "source_status": "READY",
            "source_accuracy_mm": 0.001,
            "source_accuracy_basis": "declared model precision; not measured physical accuracy",
        },
        "coordinates": {
            "input_unit": "mm",
            "native_unit": "m",
            "m_per_input_unit": 0.001,
            "reference_frame": "fixture-mm",
            "native_frame": "fixture-m",
            "native_frame_from_reference_frame_m": [
                [1.0, 0.0, 0.0, 0.01],
                [0.0, 1.0, 0.0, 0.02],
                [0.0, 0.0, 1.0, 0.03],
                [0.0, 0.0, 0.0, 1.0],
            ],
        },
        "points": {
            "payload_path": str(root / "points.csv"),
            "payload_sha256": digest("c"),
            "order": "exact_payload_row_order",
            "point_ids_order_sha256": _point_hash(point_ids),
            "rows": rows,
        },
        "candidate_generation": {
            "metric": "euclidean_reference_distance",
            "horizon_m": 0.01,
            "candidate_rule": "exact_pinned_Peridynamics_find_bonds_convention",
            "scope": "one_body_one_uniform_horizon",
        },
        "raw_euclidean_candidate_pairs": decisions,
        "counts": {
            "raw_euclidean_candidate_pairs": len(decisions),
            "initial_material_bond_pairs": count_by_state["IN_MATERIAL"],
            "geometry_excluded_pairs": count_by_state["OUT_OR_VOID"],
            "unknown_pairs": count_by_state["UNKNOWN"],
        },
        "pair_map_sha256": _pair_hash(decisions),
        "mask": {
            "included_pairs": included,
            "geometry_excluded_pairs": excluded,
            "unknown_pairs": unknown,
            "mask_sha256": _mask_hash(included),
        },
        "evidence": {
            "classifier_result_path": str(root / "result.json"),
            "classifier_result_sha256": digest("d"),
            "classifier_source_path": str(root / "bond_visibility.py"),
            "classifier_source_sha256": digest("b"),
            "raw_hit_interval_sidecar_path": str(root / "sidecar.json"),
            "raw_hit_interval_sidecar_sha256": digest("e"),
            "sidecar_selection": {"test_nodeid": "test_fixture::report", "call_index": 0},
            "visibility_algorithm": "fixture-algorithm",
        },
        "initial_physics_events": {
            "evidence_status": "PLANNED_NOT_OBSERVED_NATIVE_NOT_RUN",
            "planned_physical_broken_bond_pairs": 0,
            "geometry_exclusions_are_physics_events": False,
        },
    }


def _write_contract(path: Path, contract: dict) -> bytes:
    raw = (json.dumps(contract, ensure_ascii=False, allow_nan=False, separators=(",", ":")) + "\n").encode("utf-8")
    path.write_bytes(raw)
    return raw


class NativeMaskTransportTests(unittest.TestCase):
    def test_actual_json_sha_and_complete_toml_readback(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            audit_path = root / "audit.json"
            output_path = root / "native.toml"
            raw = _write_contract(audit_path, _valid_contract(root))
            receipt = transport.compile_native_mask_transport(audit_path, output_path)
            parsed = tomllib.loads(output_path.read_text(encoding="utf-8"))
            self.assertEqual(receipt["audit_contract_json_sha256"], hashlib.sha256(raw).hexdigest())
            self.assertEqual(parsed["transport"]["audit_contract_json_sha256"], hashlib.sha256(raw).hexdigest())
            self.assertEqual(receipt["native_toml_sha256"], hashlib.sha256(output_path.read_bytes()).hexdigest())
            self.assertEqual(parsed["points"]["rows"][1]["point_id"], 8)
            self.assertEqual(parsed["points"]["rows"][1]["native_index"], 2)
            self.assertEqual(parsed["raw_euclidean_candidate_pairs"][0]["native_i"], 1)
            self.assertFalse(parsed["transport"]["initialization_authorized"])
            self.assertEqual(parsed["native_gate"]["status"], "UNVERIFIED")
            self.assertEqual(receipt["status"], "TRANSPORT_ONLY_NATIVE_UNVERIFIED")

    def test_duplicate_keys_and_nonfinite_json_are_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            audit = root / "audit.json"
            output = root / "native.toml"
            audit.write_text('{"schema":"one","schema":"two"}', encoding="utf-8")
            with self.assertRaises(transport.TransportError):
                transport.compile_native_mask_transport(audit, output)
            audit.write_text('{"number":NaN}', encoding="utf-8")
            with self.assertRaises(transport.TransportError):
                transport.compile_native_mask_transport(audit, output)

    def test_pair_map_and_included_mask_checksums_are_recomputed(self):
        for field in ("pair_map_sha256", "mask.mask_sha256"):
            with self.subTest(field=field), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                contract = _valid_contract(root)
                if field == "pair_map_sha256":
                    contract["pair_map_sha256"] = "f" * 64
                else:
                    contract["mask"]["mask_sha256"] = "f" * 64
                audit = root / "audit.json"
                _write_contract(audit, contract)
                with self.assertRaises(transport.TransportError):
                    transport.compile_native_mask_transport(audit, root / "native.toml")

    def test_external_ids_must_map_to_sequential_native_indices_and_point_hash(self):
        for mutation in ("native_index", "point_ids_order_sha256"):
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                contract = _valid_contract(root)
                if mutation == "native_index":
                    contract["points"]["rows"][1]["native_index"] = 1
                else:
                    contract["points"]["point_ids_order_sha256"] = "f" * 64
                audit = root / "audit.json"
                _write_contract(audit, contract)
                with self.assertRaises(transport.TransportError):
                    transport.compile_native_mask_transport(audit, root / "native.toml")

    def test_unsorted_external_ids_keep_row_order_and_map_canonical_pairs_explicitly(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            contract = _valid_contract(
                root,
                point_ids=[8, 0, 20],
                states=[
                    (0, 8, "IN_MATERIAL"),
                    (0, 20, "IN_MATERIAL"),
                    (8, 20, "IN_MATERIAL"),
                ],
            )
            audit = root / "audit.json"
            _write_contract(audit, contract)
            output = root / "native.toml"
            transport.compile_native_mask_transport(audit, output)
            parsed = tomllib.loads(output.read_text(encoding="utf-8"))

            self.assertEqual(
                [(row["point_id"], row["native_index"]) for row in parsed["points"]["rows"]],
                [(8, 1), (0, 2), (20, 3)],
            )
            self.assertEqual(
                [(row["point_i"], row["point_j"], row["native_i"], row["native_j"])
                 for row in parsed["raw_euclidean_candidate_pairs"]],
                [(0, 8, 2, 1), (0, 20, 2, 3), (8, 20, 1, 3)],
            )
            self.assertEqual(
                [(row["native_i"], row["native_j"]) for row in parsed["mask"]["included_pairs"]],
                [(2, 1), (2, 3), (1, 3)],
            )

    def test_unknown_candidate_remains_blocked_and_never_authorizes_initialization(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            contract = _valid_contract(root, states=[(0, 8, "UNKNOWN")])
            audit = root / "audit.json"
            _write_contract(audit, contract)
            output = root / "native.toml"
            receipt = transport.compile_native_mask_transport(audit, output)
            parsed = tomllib.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(parsed["transport"]["native_mask_eligibility"], "BLOCKED_UNKNOWN")
            self.assertEqual(parsed["mask"]["unknown_pairs"][0]["point_i"], 0)
            self.assertFalse(parsed["transport"]["initialization_authorized"])
            self.assertEqual(receipt["native_gate_status"], "UNVERIFIED")

    def test_empty_candidate_mask_and_zero_degree_states_are_preserved_as_blocked(self):
        cases = [
            ([0, 8], [], "BLOCKED_EMPTY_CANDIDATE_SET"),
            ([0, 8], [(0, 8, "OUT_OR_VOID")], "BLOCKED_EMPTY_MASK"),
            ([0, 8, 20], [(0, 8, "IN_MATERIAL"), (0, 20, "OUT_OR_VOID"), (8, 20, "OUT_OR_VOID")], "BLOCKED_ZERO_DEGREE"),
        ]
        for point_ids, states, expected in cases:
            with self.subTest(expected=expected), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                contract = _valid_contract(root, point_ids=point_ids, states=states)
                audit = root / "audit.json"
                _write_contract(audit, contract)
                output = root / "native.toml"
                transport.compile_native_mask_transport(audit, output)
                parsed = tomllib.loads(output.read_text(encoding="utf-8"))
                self.assertEqual(parsed["transport"]["native_mask_eligibility"], expected)
                self.assertFalse(parsed["native_gate"]["initialization_authorized"])

    def test_unsupported_fields_and_native_pass_claims_are_rejected(self):
        for mutation in ("extra_key", "authorization"):
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                contract = _valid_contract(root)
                if mutation == "extra_key":
                    contract["unrecognized"] = "value"
                else:
                    contract["native_gate"]["initialization_authorized"] = True
                audit = root / "audit.json"
                _write_contract(audit, contract)
                with self.assertRaises(transport.TransportError):
                    transport.compile_native_mask_transport(audit, root / "native.toml")

    def test_input_size_limit_rejects_before_open(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            audit = root / "oversized.json"
            audit.write_bytes(b"{}" * 32)
            original_limit = transport.MAX_AUDIT_JSON_BYTES
            try:
                transport.MAX_AUDIT_JSON_BYTES = 8
                with patch.object(Path, "open", side_effect=AssertionError("file must not open before stat budget")):
                    with self.assertRaises(transport.TransportError):
                        transport.compile_native_mask_transport(audit, root / "native.toml")
            finally:
                transport.MAX_AUDIT_JSON_BYTES = original_limit

    def test_failed_compile_does_not_replace_existing_output(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            audit = root / "audit.json"
            output = root / "native.toml"
            output.write_bytes(b"previous transport bytes\n")
            audit.write_text('{"schema":"invalid"}', encoding="utf-8")
            with self.assertRaises(transport.TransportError):
                transport.compile_native_mask_transport(audit, output)
            self.assertEqual(output.read_bytes(), b"previous transport bytes\n")
            self.assertEqual(list(root.glob(".native.toml.*.tmp")), [])

    def test_valid_compile_refuses_to_overwrite_existing_output(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            audit = root / "audit.json"
            output = root / "native.toml"
            _write_contract(audit, _valid_contract(root))
            previous = b"reviewed transport evidence\n"
            output.write_bytes(previous)

            with self.assertRaises(transport.TransportError):
                transport.compile_native_mask_transport(audit, output)

            self.assertEqual(output.read_bytes(), previous)
            self.assertEqual(list(root.glob(".native.toml.*.tmp")), [])


if __name__ == "__main__":
    unittest.main()
