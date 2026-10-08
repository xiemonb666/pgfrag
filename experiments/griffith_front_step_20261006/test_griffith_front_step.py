"""Fast explicit-fixture tests; no solver or native mechanics are invoked."""

from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import griffith_front_step as front  # noqa: E402

from explicit_surface_state import CrackPatch, ExplicitSurfaceState, SurfaceSide  # noqa: E402


def fixture_seed() -> ExplicitSurfaceState:
    vertices = np.asarray(
        [[-0.10, 0.0, -0.05], [-0.10, 0.0, 0.05], [0.0, 0.0, 0.05], [0.0, 0.0, -0.05]],
        dtype=np.float64,
    )
    triangles = np.asarray([[0, 1, 2], [0, 2, 3]], dtype=np.int64)
    displacement = np.zeros_like(vertices)
    patch = CrackPatch(
        surface_id="fixture-surface",
        patch_id="fixture-root-patch",
        creation_step=0,
        parent_patch_id=None,
        reference_vertices_m=vertices,
        triangles=triangles,
        negative_side=SurfaceSide(-1, "fixture-body-minus", displacement),
        positive_side=SurfaceSide(+1, "fixture-body-plus", displacement),
    )
    return ExplicitSurfaceState((patch,))


def fixture_setup(*, G: float = 20.0, Gc: float = 10.0):
    snapshot = front.initialize(fixture_seed())
    evidence = front.make_fixture_driving_evidence(snapshot.active_surface, G=G, Gc=Gc)
    frame = front.FrontFrame(0.0, 0.0, -0.05, 0.05, "+x")
    bounds = front.DomainBounds(-0.20, 0.10, -0.01, 0.01, -0.10, 0.10)
    return snapshot, evidence, frame, bounds


def last_status(snapshot: front.ControllerSnapshot) -> str:
    return snapshot.event_ledger[-1].status


class GriffithFrontStepTests(unittest.TestCase):
    def propose(self, snapshot, evidence, frame, bounds, **kwargs):
        return front.propose_step(
            snapshot, evidence, frame=frame,
            requested_advance_m=kwargs.pop("requested_advance_m", 0.01),
            step_cap_m=kwargs.pop("step_cap_m", 0.02),
            bounds=bounds, **kwargs,
        )

    def test_arrest_when_G_below_Gc(self):
        snapshot, evidence, frame, bounds = fixture_setup(G=9.0, Gc=10.0)
        result = self.propose(snapshot, evidence, frame, bounds)
        self.assertEqual(last_status(result), "ARRESTED")
        self.assertIsNone(result.pending_candidate_id)
        self.assertEqual(result.active_surface.state_sha256, snapshot.active_surface.state_sha256)

    def test_G_above_Gc_proposes_geometry_only(self):
        snapshot, evidence, frame, bounds = fixture_setup()
        result = self.propose(snapshot, evidence, frame, bounds)
        self.assertEqual(last_status(result), "PROPOSED_WAITING_REANALYSIS")
        candidate = result.candidates[-1]
        self.assertEqual(candidate.status, "PROPOSED_WAITING_REANALYSIS")
        self.assertEqual(candidate.strip.area_m2, 0.001)
        self.assertEqual(result.active_surface.state_sha256, snapshot.active_surface.state_sha256)
        self.assertEqual(len(result.committed_versions), 1)

    def test_no_automatic_commit_and_reject_retains_old_version(self):
        snapshot, evidence, frame, bounds = fixture_setup()
        proposed = self.propose(snapshot, evidence, frame, bounds)
        candidate_id = proposed.pending_candidate_id
        self.assertEqual(len(proposed.committed_versions), 1)
        self.assertEqual(proposed.active_surface.version_id, snapshot.active_surface.version_id)
        rejected = front.reject_candidate(proposed, candidate_id, "review fixture rejection")
        self.assertEqual(last_status(rejected), "REJECTED_BY_REVIEW")
        self.assertIsNone(rejected.pending_candidate_id)
        self.assertEqual(rejected.active_surface.version_id, snapshot.active_surface.version_id)

    def test_missing_energy_rolls_back_to_old_committed_surface(self):
        snapshot, evidence, frame, bounds = fixture_setup()
        proposed = self.propose(snapshot, evidence, frame, bounds)
        candidate = proposed.candidates[-1]
        state = front.make_fixture_reanalyzed_state(proposed.active_surface, candidate)
        incomplete = front.make_fixture_reanalysis_evidence(
            candidate, state, proposed.active_surface.state_sha256, energy_complete=False
        )
        result = front.confirm_candidate(
            proposed, candidate.candidate_id, reanalyzed_surface=state, evidence=incomplete
        )
        self.assertEqual(last_status(result), "NEEDS_RECOMPUTE")
        self.assertEqual(result.active_surface.version_id, snapshot.active_surface.version_id)
        self.assertEqual(len(result.committed_versions), 1)
        self.assertEqual(result.pending_candidate_id, candidate.candidate_id)

    def test_bound_fixture_reanalysis_can_exercise_commit_branch_without_solver_claim(self):
        snapshot, evidence, frame, bounds = fixture_setup()
        proposed = self.propose(snapshot, evidence, frame, bounds)
        candidate = proposed.candidates[-1]
        state = front.make_fixture_reanalyzed_state(
            proposed.active_surface, candidate, displacement_z_m=2.0e-6
        )
        receipt = front.make_fixture_reanalysis_evidence(
            candidate, state, proposed.active_surface.state_sha256
        )
        result = front.confirm_candidate(
            proposed, candidate.candidate_id, reanalyzed_surface=state, evidence=receipt
        )
        self.assertEqual(last_status(result), "COMMITTED_FIXTURE_ONLY")
        self.assertEqual(len(result.committed_versions), 2)
        self.assertEqual(result.active_surface.version_id, result.committed_versions[-1].version_id)
        self.assertEqual(result.candidates[-1].status, "COMMITTED_FIXTURE_ONLY")
        self.assertEqual(result.candidates[-1].evidence_origin, "EXPLICIT_FIXTURE_NOT_SOLVER")
        self.assertEqual(result.active_surface.state.mechanics_qualification, "NOT_COUPLED")
        self.assertIn("energy_closure_residual_J_external_minus_delta_potential_minus_Gc_delta_area", dict(result.candidates[-1].energy_ledger))
        self.assertEqual(len({event.event_id for event in result.event_ledger}), len(result.event_ledger))

    def test_negative_Gc_is_rejected(self):
        snapshot, evidence, frame, bounds = fixture_setup(G=20.0, Gc=-1.0)
        result = self.propose(snapshot, evidence, frame, bounds)
        self.assertEqual(last_status(result), "REJECTED_INVALID_GC")
        self.assertIsNone(result.pending_candidate_id)

    def test_reversed_direction_is_rejected(self):
        snapshot, evidence, frame, bounds = fixture_setup()
        reversed_frame = front.FrontFrame(0.0, 0.0, -0.05, 0.05, "-x")
        result = self.propose(snapshot, evidence, reversed_frame, bounds)
        self.assertEqual(last_status(result), "REJECTED_REVERSED_DIRECTION")
        self.assertIsNone(result.pending_candidate_id)

    def test_out_of_domain_candidate_is_rejected(self):
        snapshot, evidence, frame, _ = fixture_setup()
        bounds = front.DomainBounds(-0.20, 0.005, -0.01, 0.01, -0.10, 0.10)
        result = self.propose(snapshot, evidence, frame, bounds, requested_advance_m=0.01)
        self.assertEqual(last_status(result), "REJECTED_OUT_OF_DOMAIN")
        self.assertIsNone(result.pending_candidate_id)

    def test_seed_frame_mismatch_is_rejected(self):
        snapshot, evidence, _, bounds = fixture_setup()
        bad_frame = front.FrontFrame(0.001, 0.0, -0.05, 0.05, "+x")
        result = self.propose(snapshot, evidence, bad_frame, bounds)
        self.assertEqual(last_status(result), "REJECTED_SEED_FRAME_MISMATCH")
        self.assertIsNone(result.pending_candidate_id)

    def test_missing_driving_evidence_is_unknown(self):
        snapshot, _, frame, bounds = fixture_setup()
        result = self.propose(snapshot, None, frame, bounds)
        self.assertEqual(last_status(result), "UNKNOWN")

    def test_unreviewed_contour_or_path_sensitivity_is_unknown(self):
        snapshot, _, frame, bounds = fixture_setup()
        evidence = front.make_fixture_driving_evidence(
            snapshot.active_surface, G=20.0, Gc=10.0, path_status="NOT_CHECKED"
        )
        result = self.propose(snapshot, evidence, frame, bounds)
        self.assertEqual(last_status(result), "UNKNOWN")
        self.assertIsNone(result.pending_candidate_id)

    def test_stale_field_needs_recompute(self):
        snapshot, _, frame, bounds = fixture_setup()
        evidence = front.make_fixture_driving_evidence(snapshot.active_surface, G=20.0, Gc=10.0)
        evidence = front.replace(evidence, bound_surface_state_sha256="0" * 64)
        result = self.propose(snapshot, evidence, frame, bounds)
        self.assertEqual(last_status(result), "NEEDS_RECOMPUTE")

    def test_uncoupled_field_needs_recompute(self):
        snapshot, _, frame, bounds = fixture_setup()
        evidence = front.make_fixture_driving_evidence(
            snapshot.active_surface, G=20.0, Gc=10.0, coupling_status="NOT_COUPLED"
        )
        result = self.propose(snapshot, evidence, frame, bounds)
        self.assertEqual(last_status(result), "NEEDS_RECOMPUTE")

    def test_requested_advance_cannot_exceed_cap(self):
        snapshot, evidence, frame, bounds = fixture_setup()
        result = self.propose(snapshot, evidence, frame, bounds, requested_advance_m=0.03, step_cap_m=0.02)
        self.assertEqual(last_status(result), "REJECTED_STEP_CAP")

    def test_equal_G_and_Gc_does_not_propose(self):
        snapshot, evidence, frame, bounds = fixture_setup(G=10.0, Gc=10.0)
        result = self.propose(snapshot, evidence, frame, bounds)
        self.assertEqual(last_status(result), "ARRESTED_AT_EQUALITY")
        self.assertIsNone(result.pending_candidate_id)


if __name__ == "__main__":
    unittest.main(verbosity=2)
