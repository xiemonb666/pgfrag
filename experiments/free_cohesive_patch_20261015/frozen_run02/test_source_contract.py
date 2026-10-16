from __future__ import annotations

import json
import unittest
from pathlib import Path

import numpy as np

from free_patch_contract import (
    input_interface_contract,
    quadrature_signed_work,
    read_gmsh22,
    sha256,
    surface_nodes,
)


HERE = Path(__file__).resolve().parent
CASE = json.loads((HERE / "case.json").read_text(encoding="utf-8"))


class FrozenSourceContractTests(unittest.TestCase):
    def test_mesh_and_material_hashes_match_frozen_case(self) -> None:
        self.assertEqual(sha256(HERE / "block_2x2x4mm.msh"), CASE["geometry"]["mesh_sha256"])
        self.assertEqual(sha256(HERE / "material.dat"), CASE["material"]["material_sha256"])

    def test_physical_interface_is_eight_two_parent_positive_area_faces(self) -> None:
        mesh = read_gmsh22(HERE / "block_2x2x4mm.msh")
        interface = input_interface_contract(mesh)
        self.assertEqual(len(mesh["nodes"]), 45)
        self.assertEqual(len(mesh["tetrahedra"]), 96)
        self.assertEqual(interface["input_face_count"], 8)
        self.assertAlmostEqual(interface["rest_area_sum_m2"], 4e-6, delta=1e-12)
        for face in interface["faces"].values():
            self.assertEqual(len(face["parent_input_tetrahedron_ids"]), 2)
            self.assertGreater(face["positive_rest_area_m2"], 0.0)

    def test_boundary_groups_and_interior_node_are_exact_input_id_contracts(self) -> None:
        mesh = read_gmsh22(HERE / "block_2x2x4mm.msh")
        self.assertEqual(len(surface_nodes(mesh, "top")), 9)
        self.assertEqual(len(surface_nodes(mesh, "bottom")), 9)
        node_id = CASE["boundary_conditions"]["geometric_interior_input_node_id_zero_based"]
        np.testing.assert_array_equal(
            mesh["nodes"][node_id],
            np.asarray(CASE["boundary_conditions"]["geometric_interior_input_node_coordinate_m"]),
        )

    def test_loading_ramp_exactly_reaches_hold_value(self) -> None:
        loading = CASE["loading"]
        ramp = int(loading["ramp_steps"])
        target = float(loading["top_target_displacement_m"])
        self.assertEqual(ramp + int(loading["hold_steps"]), int(loading["total_steps"]))
        self.assertAlmostEqual(0.5 * target * (1.0 - np.cos(np.pi * ramp / ramp)), target)
        self.assertEqual(loading["bulk_stable_step_fraction_limit"], 0.03)
        self.assertIn("not a total cohesive-plus-bulk CFL gate", loading["bulk_stable_step_gate"])

    def test_total_stability_is_explicitly_pending_and_blocks_solver_step(self) -> None:
        stability = CASE["stability"]
        self.assertEqual(stability["total_stability_status"], "NOT_ESTABLISHED_REQUIRES_ROOT_BOUND")
        self.assertIsNone(stability["lambda_max_upper_bound_per_s2"])
        self.assertIn("Before solver step 1", stability["gate"])
        source = (HERE / "run_native_free_patch.py").read_text(encoding="utf-8")
        self.assertIn("2.0 / float(np.sqrt(float(lambda_upper)))", source)
        self.assertIn("if not total_gate:", source)
        self.assertIn("Stop before solver step 1.", source)
        self.assertIn('"bulk_check_is_total_cohesive_cfl": False', source)


class SignedWorkContractTests(unittest.TestCase):
    def test_native_work_uses_energy_conjugate_plus_sign(self) -> None:
        opening0 = np.zeros((1, 3))
        traction0 = np.zeros((1, 3))
        opening1 = np.asarray([[0.0, 0.0, 2.0]])
        traction1 = np.asarray([[0.0, 0.0, 3.0]])
        work = quadrature_signed_work(opening0, traction0, opening1, traction1, np.asarray([0.5]))
        np.testing.assert_allclose(work, [1.5])
        # Reversing both face orientation vectors leaves their conjugate dot product unchanged.
        work_reversed = quadrature_signed_work(-opening0, -traction0, -opening1, -traction1, np.asarray([0.5]))
        np.testing.assert_allclose(work_reversed, work)

    def test_unloading_work_keeps_its_native_negative_sign(self) -> None:
        opening0 = np.asarray([[0.0, 0.0, 2.0]])
        traction0 = np.asarray([[0.0, 0.0, 3.0]])
        opening1 = np.asarray([[0.0, 0.0, 1.0]])
        traction1 = np.asarray([[0.0, 0.0, 2.0]])
        work = quadrature_signed_work(opening0, traction0, opening1, traction1, np.asarray([0.5]))
        self.assertLess(float(work[0]), 0.0)
        self.assertEqual(float(work[0]), -1.25)

    def test_case_does_not_allow_resisting_sign_flip_or_abs_cleanup(self) -> None:
        contract = CASE["energy_contract"]
        self.assertIn("energy-conjugate", contract["native_sign_policy"])
        self.assertIn("never apply abs()", contract["native_sign_policy"])
        self.assertIn("no negative sign is applied", contract["cohesive_traction_opening_work"])
        self.assertIn("FAIL_DIAGNOSTIC_RETAIN_RAW", CASE["acceptance"]["negative_native_dissipation"])


class SourceWiringTests(unittest.TestCase):
    def test_runner_reads_only_exact_native_energy_keys_and_no_cell_dump(self) -> None:
        source = (HERE / "run_native_free_patch.py").read_text(encoding="utf-8")
        self.assertIn('"cohesive contact"', source)
        self.assertNotIn('getEnergy("contact")', source)
        self.assertEqual(source.count("model.getEnergy(key)"), 1)
        self.assertNotIn("model.dump(", source)
        self.assertIn("model.flattenInternalReal(name, aka._ek_cohesive, aka._not_ghost)", source)
        self.assertIn("history.write(frame, row)", source)

    def test_runner_preserves_free_mask_and_has_single_attempt_no_retry(self) -> None:
        source = (HERE / "run_native_free_patch.py").read_text(encoding="utf-8")
        self.assertIn("blocked[:, :] = False", source)
        self.assertIn("native blocked-DOF mask changed", source)
        self.assertIn('"retry_policy": "none"', source)
        self.assertNotIn("np.abs(energies", source)
        self.assertIn("native_work_vs_total_residual_J_raw", source)

    def test_analyzer_treats_negative_energy_and_history_mismatch_as_fail_diagnostics(self) -> None:
        source = (HERE / "analyze_native_free_patch.py").read_text(encoding="utf-8")
        self.assertIn("NEGATIVE_NATIVE_D_MATCHES_SIGNED_HISTORY_IDENTITY", source)
        self.assertIn("NEGATIVE_NATIVE_D_WITH_WORK_HISTORY_MISMATCH", source)
        self.assertIn("raw_native_dissipation_not_below_roundoff_floor", source)
        self.assertIn("independent_signed_interface_work_matches_native_total_history", source)
        self.assertNotIn("np.abs(energies[:, 2])", source)
        self.assertIn("parent_reviewed_total_cohesive_plus_bulk_stability_bound", source)


if __name__ == "__main__":
    unittest.main(verbosity=2)
