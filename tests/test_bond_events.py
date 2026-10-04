from __future__ import annotations

from dataclasses import replace
import hashlib
import json
import math
from pathlib import Path
import unittest

from pgfrag.bond_events import (
    BondEventValidationError,
    DirectedBondRecord,
    NativeBondFrame,
    NativeStateProvenance,
    export_native_bond_events,
)


FIXTURE_PATH = (
    Path(__file__).resolve().parents[1]
    / "tmp"
    / "review_20261004"
    / "bond_events"
    / "native_shaped_fixture.json"
)


def _raw_fixture() -> dict:
    return json.loads(FIXTURE_PATH.read_text(encoding="utf-8"))


def _frames(raw: dict) -> tuple[NativeBondFrame, ...]:
    return tuple(
        NativeBondFrame(
            step=frame["step"],
            time_s=frame["time_s"],
            step_dt_s=frame.get("step_dt_s"),
            directed_bonds=tuple(DirectedBondRecord(**record) for record in frame["directed_bonds"]),
        )
        for frame in raw["frames"]
    )


def _provenance() -> NativeStateProvenance:
    archive_hash = hashlib.sha256(FIXTURE_PATH.read_bytes()).hexdigest()
    return NativeStateProvenance(
        native_source="hand-authored native-shaped fixture; not solver output",
        native_interface="DirectedBondRecord contract fixture",
        native_source_sha256=None,
        frame_archive_sha256=archive_hash,
    )


def _export(raw: dict, **overrides):
    kwargs = {
        "point_ids": [int(point_id) for point_id in raw["point_ids"]],
        "rest_positions_m": {
            int(point_id): value for point_id, value in raw["rest_positions_m"].items()
        },
        "point_volumes_m3": {
            int(point_id): value for point_id, value in raw["point_volumes_m3"].items()
        },
        "density_kg_m3": raw["density_kg_m3"],
        "initial_bonds": [tuple(edge) for edge in raw["initial_bonds"]],
        "frames": _frames(raw),
        "expected_final_step": raw["expected_final_step"],
        "expected_dt_s": raw["expected_dt_s"],
        "provenance": _provenance(),
    }
    kwargs.update(overrides)
    return export_native_bond_events(**kwargs)


def _replace_pair_state(frame: NativeBondFrame, pair: tuple[int, int], active: bool) -> NativeBondFrame:
    first, second = pair
    records = tuple(
        replace(record, active=active)
        if (record.source_point_id, record.target_point_id) in ((first, second), (second, first))
        else record
        for record in frame.directed_bonds
    )
    return replace(frame, directed_bonds=records)


class NativeBondEventTests(unittest.TestCase):
    def test_first_true_to_false_event_has_step_time_and_directed_field_lineage(self) -> None:
        result = _export(_raw_fixture())
        self.assertEqual(result.timeline_status, "KNOWN")
        self.assertEqual(result.time_step_consistency, "CONSTANT_DT_VALIDATED")
        self.assertEqual(len(result.break_events), 1)
        event = result.break_events[0]
        self.assertEqual((event.pair.point_i_global, event.pair.point_j_global), (20, 30))
        self.assertEqual((event.step, event.time_s, event.previous_step, event.previous_time_s), (2, 2e-6, 1, 1e-6))
        self.assertEqual(
            tuple(ref.native_field for ref in event.broken_native_fields),
            ("storage.bond_active[3]", "storage.bond_active[4]"),
        )
        self.assertEqual(
            tuple((ref.source_point_id, ref.target_point_id) for ref in event.previously_active_native_fields),
            ((20, 30), (30, 20)),
        )

    def test_components_include_isolated_points_and_all_point_mass(self) -> None:
        result = _export(_raw_fixture())
        self.assertEqual(len(result.frame_components), 3)
        for frame in result.frame_components:
            components = frame.components
            self.assertEqual(sum(len(component.point_ids) for component in components), 4)
            self.assertAlmostEqual(sum(component.mass_kg for component in components), 8e-6)
        final_components = result.frame_components[-1].components
        self.assertEqual(
            [(component.component_id, component.point_ids) for component in final_components],
            [(10, (10, 20)), (30, (30,)), (40, (40,))],
        )
        for component, expected_mass in zip(final_components, (4e-6, 2e-6, 2e-6), strict=True):
            self.assertAlmostEqual(component.mass_kg, expected_mass)

    def test_initially_broken_pair_is_initial_state_not_a_new_event(self) -> None:
        raw = _raw_fixture()
        frames = tuple(_replace_pair_state(frame, (10, 20), False) for frame in _frames(raw))
        result = _export(raw, frames=frames)
        self.assertEqual(
            [(item.pair.point_i_global, item.pair.point_j_global) for item in result.initial_state_broken_bonds],
            [(10, 20)],
        )
        self.assertEqual(
            [(event.pair.point_i_global, event.pair.point_j_global) for event in result.break_events],
            [(20, 30)],
        )

    def test_endpoint_only_has_unknown_timeline_and_no_inferred_events(self) -> None:
        raw = _raw_fixture()
        endpoint = (_frames(raw)[-1],)
        result = _export(raw, frames=endpoint)
        self.assertEqual(result.timeline_status, "UNKNOWN")
        self.assertEqual(result.break_events, ())
        self.assertEqual(len(result.frame_components), 1)
        self.assertIn("only_one_or_zero_native_frame; endpoint chronology is unavailable", result.timeline_unknown_reasons)
        encoded = result.to_dict()
        self.assertEqual(encoded["interpretation"]["R2"], "UNKNOWN_NOT_DERIVED")
        self.assertEqual(encoded["interpretation"]["verified_rigid_fragments"], "UNKNOWN_NOT_DERIVED")
        self.assertEqual(encoded["break_events_status"], "UNKNOWN_NO_INFERENCE")
        self.assertEqual(encoded["native_state_provenance"]["qualification"], "SUPPLIED_NATIVE_STATE_CALLER_PROVIDED_UNVERIFIED")
        self.assertEqual(
            encoded["native_state_provenance"]["frame_archive_sha256"],
            hashlib.sha256(FIXTURE_PATH.read_bytes()).hexdigest(),
        )
        self.assertNotIn("fragment_count", encoded)

    def test_missing_time_keeps_timeline_unknown_and_suppresses_events(self) -> None:
        raw = _raw_fixture()
        frames = list(_frames(raw))
        frames[1] = replace(frames[1], time_s=None)
        result = _export(raw, frames=tuple(frames))
        self.assertEqual(result.timeline_status, "UNKNOWN")
        self.assertEqual(result.break_events, ())
        self.assertIn("one_or_more_native_frame_times_unavailable", result.timeline_unknown_reasons)
        self.assertIsNone(result.frame_components[1].time_s)

    def test_finite_times_without_dt_evidence_do_not_qualify_timeline(self) -> None:
        raw = _raw_fixture()
        result = _export(raw, expected_dt_s=None)
        self.assertEqual(result.timeline_status, "UNKNOWN")
        self.assertEqual(result.time_step_consistency, "UNKNOWN")
        self.assertEqual(result.break_events, ())

    def test_explicit_variable_dt_intervals_validate_time_and_preserve_event(self) -> None:
        raw = _raw_fixture()
        frames = list(_frames(raw))
        frames = [
            replace(frames[0], time_s=7e-6, step_dt_s=None),
            replace(frames[1], time_s=8e-6, step_dt_s=1e-6),
            replace(frames[2], time_s=10e-6, step_dt_s=2e-6),
        ]
        result = _export(raw, frames=tuple(frames), expected_dt_s=None)
        self.assertEqual(result.timeline_status, "KNOWN")
        self.assertEqual(result.time_step_consistency, "VARIABLE_DT_INTERVALS_VALIDATED")
        self.assertEqual((result.break_events[0].step, result.break_events[0].time_s), (2, 10e-6))

    def test_large_origin_unresolvable_clock_is_unknown_for_constant_and_variable_dt(self) -> None:
        raw = _raw_fixture()
        times = (1.0, 1.0000000000000002, 1.0000000000000004)
        for constant_mode in (True, False):
            frames = tuple(
                replace(
                    frame,
                    time_s=timestamp,
                    step_dt_s=None if index == 0 or constant_mode else 1e-16,
                )
                for index, (frame, timestamp) in enumerate(zip(_frames(raw), times, strict=True))
            )
            result = _export(
                raw,
                frames=frames,
                expected_dt_s=1e-16 if constant_mode else None,
            )
            self.assertEqual(result.timeline_status, "UNKNOWN_TIME_RESOLUTION")
            self.assertEqual(result.break_events, ())
            self.assertIn(
                "timestamp_rounding_uncertainty_exceeds_relative_interval_limit",
                result.timeline_unknown_reasons,
            )
            self.assertEqual(result.to_dict()["timestamp_rounding_relative_error_limit"], 1e-6)

    def test_small_dt_at_zero_and_resolvable_dt_at_large_origin_are_known(self) -> None:
        raw = _raw_fixture()
        base_frames = _frames(raw)
        small_at_zero = tuple(
            replace(frame, time_s=timestamp)
            for frame, timestamp in zip(base_frames, (0.0, 1e-16, 2e-16), strict=True)
        )
        small_result = _export(raw, frames=small_at_zero, expected_dt_s=1e-16)
        self.assertEqual(small_result.timeline_status, "KNOWN")
        self.assertEqual((small_result.break_events[0].step, small_result.break_events[0].time_s), (2, 2e-16))

        resolvable_at_large_origin = tuple(
            replace(frame, time_s=timestamp)
            for frame, timestamp in zip(base_frames, (1.0, 1.000000001, 1.000000002), strict=True)
        )
        large_result = _export(raw, frames=resolvable_at_large_origin, expected_dt_s=1e-9)
        self.assertEqual(large_result.timeline_status, "KNOWN")
        self.assertEqual(large_result.time_step_consistency, "CONSTANT_DT_VALIDATED")

    def test_constant_dt_checks_cumulative_time_drift(self) -> None:
        raw = _raw_fixture()
        initial_records = _frames(raw)[0].directed_bonds
        times = [1.0]
        for _ in range(4):
            times.append(times[-1] + 1e-9)
        frames = tuple(
            NativeBondFrame(step, timestamp, initial_records)
            for step, timestamp in enumerate(times)
        )
        with self.assertRaisesRegex(BondEventValidationError, "cumulative native time"):
            _export(raw, frames=frames, expected_final_step=4, expected_dt_s=1e-9)

    def test_continuous_prefix_is_unknown_but_preserves_components_mass_and_raw_lineage(self) -> None:
        raw = _raw_fixture()
        result = _export(raw, expected_final_step=3)
        self.assertEqual(result.timeline_status, "UNKNOWN")
        self.assertEqual(result.break_events, ())
        self.assertIn("continuous_frame_prefix_does_not_reach_expected_final_step", result.timeline_unknown_reasons)
        self.assertEqual(len(result.frame_components), 3)
        self.assertEqual([len(frame.native_directed_bonds) for frame in result.frame_components], [4, 4, 4])
        final_frame = result.frame_components[-1]
        self.assertEqual(final_frame.step, 2)
        self.assertEqual(len(final_frame.components), 3)
        self.assertAlmostEqual(sum(component.mass_kg for component in final_frame.components), 8e-6)
        self.assertEqual(len(final_frame.native_directed_bonds), 4)
        self.assertEqual(final_frame.native_directed_bonds[2].native_field, "storage.bond_active[3]")
        encoded = result.to_dict()
        self.assertEqual(encoded["break_events_status"], "UNKNOWN_NO_INFERENCE")
        self.assertEqual(len(encoded["frame_components"][-1]["native_directed_bond_states"]), 4)
        self.assertEqual(encoded["frame_components"][-1]["native_directed_bond_states"][2]["active"], False)
        self.assertNotIn("no_fracture", encoded)

    def test_reciprocal_mismatch_is_rejected(self) -> None:
        raw = _raw_fixture()
        frames = list(_frames(raw))
        frame = frames[1]
        records = list(frame.directed_bonds)
        records[0] = replace(records[0], active=False)
        frames[1] = replace(frame, directed_bonds=tuple(records))
        with self.assertRaisesRegex(BondEventValidationError, "reciprocal native bond states disagree"):
            _export(raw, frames=tuple(frames))

    def test_duplicate_directed_state_is_rejected(self) -> None:
        raw = _raw_fixture()
        frames = list(_frames(raw))
        frames[1] = replace(frames[1], directed_bonds=frames[1].directed_bonds + (frames[1].directed_bonds[0],))
        with self.assertRaisesRegex(BondEventValidationError, "duplicate directed bond state"):
            _export(raw, frames=tuple(frames))

    def test_missing_reciprocal_directed_state_is_rejected(self) -> None:
        raw = _raw_fixture()
        frames = list(_frames(raw))
        frames[1] = replace(frames[1], directed_bonds=frames[1].directed_bonds[:-1])
        with self.assertRaisesRegex(BondEventValidationError, "missing directed native bond state"):
            _export(raw, frames=tuple(frames))

    def test_illegal_healing_is_rejected(self) -> None:
        raw = _raw_fixture()
        frames = list(_frames(raw))
        frames[0] = _replace_pair_state(frames[0], (10, 20), False)
        with self.assertRaisesRegex(BondEventValidationError, "illegal native bond healing"):
            _export(raw, frames=tuple(frames))

    def test_missing_frame_and_reversed_frames_are_rejected(self) -> None:
        raw = _raw_fixture()
        frames = _frames(raw)
        with self.assertRaisesRegex(BondEventValidationError, "missing native solver frame"):
            _export(raw, frames=(frames[0], frames[2]))
        with self.assertRaisesRegex(BondEventValidationError, "strictly increasing"):
            _export(raw, frames=(frames[1], frames[0]))

    def test_time_that_disagrees_with_constant_dt_is_rejected(self) -> None:
        raw = _raw_fixture()
        frames = list(_frames(raw))
        frames[2] = replace(frames[2], time_s=3e-6)
        with self.assertRaisesRegex(BondEventValidationError, "disagrees with expected_dt_s"):
            _export(raw, frames=tuple(frames))

    def test_self_bond_duplicate_initial_pair_and_nonfinite_geometry_are_rejected(self) -> None:
        raw = _raw_fixture()
        frames = list(_frames(raw))
        bad_record = DirectedBondRecord(10, 10, True, "storage.bond_active[self]")
        frames[0] = replace(frames[0], directed_bonds=frames[0].directed_bonds + (bad_record,))
        with self.assertRaisesRegex(BondEventValidationError, "self-bond"):
            _export(raw, frames=tuple(frames))
        with self.assertRaisesRegex(BondEventValidationError, "duplicate canonical initial bond"):
            _export(raw, initial_bonds=raw["initial_bonds"] + [[20, 10]])
        positions = {int(key): value for key, value in raw["rest_positions_m"].items()}
        positions[10] = (math.nan, 0.0, 0.0)
        with self.assertRaisesRegex(BondEventValidationError, "must be finite"):
            _export(raw, rest_positions_m=positions)

    def test_single_material_point_is_legal_and_conserves_its_mass(self) -> None:
        raw = _raw_fixture()
        frames = (
            NativeBondFrame(0, 4.0, (), None),
            NativeBondFrame(1, 4.000001, (), None),
        )
        result = _export(
            raw,
            point_ids=[1],
            rest_positions_m={1: (0.0, 0.0, 0.0)},
            point_volumes_m3={1: 2e-9},
            density_kg_m3=3000.0,
            initial_bonds=(),
            frames=frames,
            expected_final_step=1,
            expected_dt_s=1e-6,
        )
        self.assertEqual(result.timeline_status, "KNOWN")
        self.assertEqual(result.break_events, ())
        self.assertEqual(result.frame_components[0].components[0].point_ids, (1,))
        self.assertAlmostEqual(result.frame_components[0].components[0].mass_kg, 6e-6)


if __name__ == "__main__":
    unittest.main()
