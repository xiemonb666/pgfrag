"""Focused acceptance tests for the isolated scalar cohesive material law."""

from __future__ import annotations

import math
import unittest

from cohesive_scalar_law import (
    BilinearCohesiveMaterial,
    UnsupportedModeError,
    dissipated_energy_j_per_m2,
    initial_state,
    stored_energy_j_per_m2,
    traction_pa,
    update_state,
)


def fixture_material() -> BilinearCohesiveMaterial:
    # Illustrative SI fixture only; not a ceramic calibration.
    return BilinearCohesiveMaterial(k0_pa_per_m=1.0e12, tpeak_pa=1.0e7, gc_j_per_m2=100.0)


class BilinearCohesiveLawTests(unittest.TestCase):
    def test_monotonic_envelope_area_is_gc_and_has_expected_corners(self) -> None:
        material = fixture_material()
        self.assertAlmostEqual(material.delta0_m, 1.0e-5, places=18)
        self.assertAlmostEqual(material.deltaf_m, 2.0e-5, places=18)
        count = 200
        openings = [material.deltaf_m * i / count for i in range(count + 1)]
        state = initial_state(material)
        tractions = []
        for opening in openings:
            state = update_state(material, state, opening)
            tractions.append(traction_pa(material, state))
        area = sum(
            0.5 * (left + right) * (x1 - x0)
            for x0, x1, left, right in zip(openings, openings[1:], tractions, tractions[1:])
        )
        self.assertAlmostEqual(area, material.gc_j_per_m2, delta=1.0e-11)
        self.assertAlmostEqual(tractions[100], material.tpeak_pa, delta=1.0e-7)
        self.assertEqual(tractions[-1], 0.0)
        self.assertAlmostEqual(state.dissipated_energy_j_per_m2, material.gc_j_per_m2, delta=1.0e-12)
        self.assertAlmostEqual(state.external_work_j_per_m2, material.gc_j_per_m2, delta=1.0e-11)

    def test_load_unload_reload_preserves_energy_ledger_and_never_dissipates_on_unload(self) -> None:
        material = fixture_material()
        state = initial_state(material)
        path = [0.25, 0.50, 0.80, 0.35, 0.65, 0.80, 0.94, 0.20]
        previous_damage = 0.0
        previous_kappa = 0.0
        previous_dissipation = 0.0
        damage_at_peak = 0.0
        for fraction in path:
            old_state = state
            state = update_state(material, state, fraction * material.deltaf_m)
            self.assertGreaterEqual(state.kappa_m, previous_kappa)
            self.assertGreaterEqual(state.damage + 1.0e-15, previous_damage)
            self.assertGreaterEqual(state.dissipated_energy_j_per_m2 + 1.0e-12, previous_dissipation)
            self.assertGreaterEqual(state.dissipated_energy_j_per_m2, 0.0)
            self.assertAlmostEqual(
                state.external_work_j_per_m2,
                state.stored_energy_j_per_m2 + state.dissipated_energy_j_per_m2,
                delta=2.0e-10,
            )
            if fraction < previous_kappa / material.deltaf_m:
                self.assertEqual(state.dissipated_energy_j_per_m2, old_state.dissipated_energy_j_per_m2)
            if fraction == 0.80:
                if damage_at_peak == 0.0:
                    damage_at_peak = state.damage
                else:
                    self.assertAlmostEqual(state.damage, damage_at_peak, delta=1.0e-14)
            previous_damage = state.damage
            previous_kappa = state.kappa_m
            previous_dissipation = state.dissipated_energy_j_per_m2

        # Reaching the former peak on reload leaves its damage history unchanged.
        self.assertGreater(damage_at_peak, 0.0)
        self.assertEqual(stored_energy_j_per_m2(material, state), state.stored_energy_j_per_m2)
        self.assertEqual(dissipated_energy_j_per_m2(material, state), state.dissipated_energy_j_per_m2)

    def test_complete_failure_has_zero_traction_and_gc_dissipation(self) -> None:
        material = fixture_material()
        state = update_state(material, initial_state(material), 1.2 * material.deltaf_m)
        self.assertEqual(state.damage, 1.0)
        self.assertEqual(traction_pa(material, state), 0.0)
        self.assertEqual(state.stored_energy_j_per_m2, 0.0)
        self.assertAlmostEqual(state.dissipated_energy_j_per_m2, material.gc_j_per_m2, delta=1.0e-12)
        after_unload = update_state(material, state, 0.0)
        self.assertEqual(after_unload.dissipated_energy_j_per_m2, state.dissipated_energy_j_per_m2)
        self.assertEqual(after_unload.external_work_j_per_m2, state.external_work_j_per_m2)

    def test_illegal_parameters_units_openings_and_modes_are_rejected(self) -> None:
        bad_parameter_sets = [
            (0.0, 1.0, 1.0),
            (1.0, -1.0, 1.0),
            (1.0, 1.0, math.inf),
            (1.0, 1.0, math.nan),
            ("1 MPa/mm", 1.0, 1.0),
            (True, 1.0, 1.0),
            (1.0e12, 1.0e7, 1.0),  # deltaf <= delta0
        ]
        for args in bad_parameter_sets:
            with self.subTest(args=args), self.assertRaises(ValueError):
                BilinearCohesiveMaterial(*args)

        material = fixture_material()
        state = initial_state(material)
        with self.assertRaises(UnsupportedModeError):
            update_state(material, state, -1.0e-6)
        for mode in ("compression", "mixed_mode", "mode_ii_sliding"):
            with self.subTest(mode=mode), self.assertRaises(UnsupportedModeError):
                update_state(material, state, 0.0, mode=mode)
            with self.subTest(mode=mode), self.assertRaises(UnsupportedModeError):
                traction_pa(material, state, mode=mode)

    def test_history_cannot_be_reused_with_another_material(self) -> None:
        material = fixture_material()
        state = update_state(material, initial_state(material), 0.4 * material.deltaf_m)
        other = BilinearCohesiveMaterial(2.0e12, 1.0e7, 100.0)
        with self.assertRaises(ValueError):
            update_state(other, state, state.opening_m)


if __name__ == "__main__":
    unittest.main()
