"""Targeted v2 energy-gate tests using explicit fixtures only."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import griffith_front_energy_gate_v2 as gate  # noqa: E402
import griffith_front_step as front  # noqa: E402
from test_griffith_front_step import fixture_setup  # noqa: E402


def setup_proposal(*, G: float = 15.0, Gc: float = 10.0):
    old, drive, frame, bounds = fixture_setup(G=G, Gc=Gc)
    initial = gate.EnergyGateSnapshot(old)
    proposed = gate.propose_step(
        initial, drive, frame=frame,
        requested_advance_m=0.01, step_cap_m=0.02, bounds=bounds,
    )
    candidate = proposed.controller.candidates[-1]
    reanalyzed = front.make_fixture_reanalyzed_state(proposed.active_surface, candidate, displacement_z_m=1.0e-6)
    return initial, proposed, candidate, reanalyzed


def base_evidence(proposed, candidate, state, **overrides):
    evidence = gate.make_fixture_evidence(
        candidate,
        proposed.active_surface.state_sha256,
        state,
        delta_stored_elastic_energy_J=0.0,
        external_work_increment_J=candidate.Gc_J_per_m2 * candidate.strip.area_m2,
        abs_tol_J=0.0,
        relative_tol=0.0,
        reference_energy_J=1.0,
    )
    return gate.replace(evidence, **overrides)


def last_gate(snapshot: gate.EnergyGateSnapshot) -> gate.EnergyGateEvent:
    return snapshot.energy_gate_ledger[-1]


class GriffithFrontEnergyGateV2Tests(unittest.TestCase):
    def confirm(self, proposed, candidate, state, evidence):
        return gate.confirm_candidate(
            proposed, candidate.candidate_id, reanalyzed_surface=state, evidence=evidence
        )

    def test_balanced_fixture_energy_commits_fixture_only(self):
        initial, proposed, candidate, state = setup_proposal()
        result = self.confirm(proposed, candidate, state, base_evidence(proposed, candidate, state))
        self.assertEqual(last_gate(result).status, "COMMITTED_FIXTURE_ONLY")
        self.assertEqual(result.active_surface.version_id, result.controller.committed_versions[-1].version_id)
        self.assertEqual(result.active_surface.state.mechanics_qualification, "NOT_COUPLED")
        self.assertEqual(len(result.controller.committed_versions), 2)
        self.assertEqual(result.active_surface.version_status, "COMMITTED_FIXTURE_ONLY")
        self.assertIsNone(result.pending_candidate_id)
        self.assertEqual(last_gate(result).residual_W_minus_delta_stored_minus_Gc_delta_area_J, 0.0)
        self.assertAlmostEqual(last_gate(result).derived_delta_total_potential_J, -0.01)
        self.assertEqual(result.controller.committed_versions[0].state_sha256, initial.active_surface.state_sha256)

    def test_external_reanalysis_api_is_available_but_fake_receipt_stays_experimental(self):
        _, proposed, candidate, state = setup_proposal()
        evidence = base_evidence(
            proposed, candidate, state,
            evidence_origin="EXTERNAL_RECEIPT",
            mechanical_status="REANALYZED_CURRENT",
            caller_provenance=(("test_only", "synthetic hashes; no solver invoked"),),
        )
        result = self.confirm(proposed, candidate, state, evidence)
        self.assertEqual(last_gate(result).status, "COMMITTED_EXPERIMENTAL_REANALYSIS")
        self.assertEqual(result.active_surface.version_status, "COMMITTED_EXPERIMENTAL_REANALYSIS")
        self.assertEqual(result.active_surface.qualification, "EXPERIMENTAL_UNCALIBRATED")
        self.assertEqual(result.active_surface.state.mechanics_qualification, "NOT_COUPLED")
        self.assertIn("synthetic hashes", last_gate(result).caller_provenance[0][1])

    def test_large_positive_residual_rejects_and_retains_old_surface(self):
        initial, proposed, candidate, state = setup_proposal()
        evidence = base_evidence(
            proposed, candidate, state, external_work_increment_J=1.0e6,
            abs_tol_J=1.0e-6, relative_tol=0.0,
        )
        result = self.confirm(proposed, candidate, state, evidence)
        self.assertEqual(last_gate(result).status, "NEEDS_RECOMPUTE")
        self.assertEqual(result.active_surface.version_id, initial.active_surface.version_id)
        self.assertEqual(len(result.controller.committed_versions), 1)
        self.assertEqual(result.pending_candidate_id, candidate.candidate_id)
        self.assertGreater(last_gate(result).residual_W_minus_delta_stored_minus_Gc_delta_area_J, 0.0)

    def test_large_negative_residual_rejects(self):
        initial, proposed, candidate, state = setup_proposal()
        evidence = base_evidence(
            proposed, candidate, state, delta_stored_elastic_energy_J=1.0e6,
            abs_tol_J=1.0e-6, relative_tol=0.0,
        )
        result = self.confirm(proposed, candidate, state, evidence)
        self.assertEqual(last_gate(result).status, "NEEDS_RECOMPUTE")
        self.assertEqual(result.active_surface.version_id, initial.active_surface.version_id)
        self.assertLess(last_gate(result).residual_W_minus_delta_stored_minus_Gc_delta_area_J, 0.0)

    def test_missing_tolerance_is_needs_recompute(self):
        _, proposed, candidate, state = setup_proposal()
        evidence = base_evidence(proposed, candidate, state, abs_tol_J=None)
        result = self.confirm(proposed, candidate, state, evidence)
        self.assertEqual(last_gate(result).status, "NEEDS_RECOMPUTE")
        self.assertIn("tolerance", last_gate(result).reason)

    def test_total_potential_input_is_not_bridged(self):
        initial, proposed, candidate, state = setup_proposal()
        evidence = base_evidence(
            proposed, candidate, state,
            delta_total_potential_J=-0.01,
            delta_stored_elastic_energy_J=None,
        )
        result = self.confirm(proposed, candidate, state, evidence)
        self.assertEqual(last_gate(result).status, "NEEDS_RECOMPUTE")
        self.assertIn("not_bridged", last_gate(result).reason)
        self.assertEqual(result.active_surface.version_id, initial.active_surface.version_id)
        self.assertEqual(result.to_payload()["legacy_delta_potential_bridge"], "NOT_BRIDGED")

    def test_dynamic_and_mixed_cohesive_modes_are_unknown(self):
        for mode in ("DYNAMIC", "MIXED_COHESIVE"):
            with self.subTest(mode=mode):
                initial, proposed, candidate, state = setup_proposal()
                evidence = base_evidence(proposed, candidate, state, physics_mode=mode)
                result = self.confirm(proposed, candidate, state, evidence)
                self.assertEqual(last_gate(result).status, "UNKNOWN")
                self.assertEqual(result.active_surface.version_id, initial.active_surface.version_id)

    def test_zero_tolerances_are_valid_when_residual_is_exactly_zero(self):
        _, proposed, candidate, state = setup_proposal()
        evidence = base_evidence(
            proposed, candidate, state,
            delta_stored_elastic_energy_J=0.0,
            external_work_increment_J=candidate.Gc_J_per_m2 * candidate.strip.area_m2,
            abs_tol_J=0.0, relative_tol=0.0, reference_energy_J=1.0,
        )
        result = self.confirm(proposed, candidate, state, evidence)
        self.assertEqual(last_gate(result).status, "COMMITTED_FIXTURE_ONLY")
        self.assertEqual(last_gate(result).allowed_abs_residual_J, 0.0)

    def test_negative_tolerances_and_nonpositive_reference_reject(self):
        cases = (
            {"abs_tol_J": -1.0e-9},
            {"relative_tol": -1.0e-9},
            {"reference_energy_J": 0.0},
            {"reference_energy_J": -1.0},
        )
        for changes in cases:
            with self.subTest(changes=changes):
                initial, proposed, candidate, state = setup_proposal()
                evidence = base_evidence(proposed, candidate, state, **changes)
                result = self.confirm(proposed, candidate, state, evidence)
                self.assertEqual(last_gate(result).status, "NEEDS_RECOMPUTE")
                self.assertEqual(result.active_surface.version_id, initial.active_surface.version_id)

    def test_parent_binding_mismatch_rejects(self):
        initial, proposed, candidate, state = setup_proposal()
        evidence = base_evidence(
            proposed, candidate, state,
            bound_parent_surface_state_sha256="0" * 64,
        )
        result = self.confirm(proposed, candidate, state, evidence)
        self.assertEqual(last_gate(result).status, "NEEDS_RECOMPUTE")
        self.assertIn("parent_surface_hash", last_gate(result).reason)
        self.assertEqual(result.active_surface.version_id, initial.active_surface.version_id)

    def test_units_cannot_be_replaced_by_review_like_metadata(self):
        initial, proposed, candidate, state = setup_proposal()
        evidence = base_evidence(
            proposed, candidate, state,
            energy_unit="kJ", caller_provenance=(("reviewed", "yes"),),
        )
        result = self.confirm(proposed, candidate, state, evidence)
        self.assertEqual(last_gate(result).status, "NEEDS_RECOMPUTE")
        self.assertIn("units", last_gate(result).reason)
        self.assertEqual(result.active_surface.version_id, initial.active_surface.version_id)


if __name__ == "__main__":
    unittest.main(verbosity=2)

