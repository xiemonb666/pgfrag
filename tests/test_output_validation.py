from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from pgfrag.backend import OUTPUT_DATA_DICTIONARY, _validate_solver_outputs
from pgfrag.errors import ValidationError
from .helpers import base_job, temporary_directory


def write_csv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)


def make_valid_output(root: Path) -> dict:
    write_csv(root / "points.csv",
              ["point_id", "x_ref_m", "y_ref_m", "z_ref_m", "x_current_m", "y_current_m", "z_current_m",
               "vx_mps", "vy_mps", "vz_mps", "damage", "n_active_bonds", "strain_energy_density_J_m3",
               "volume_m3", "mass_kg", "intact_bond_component_id"],
              [[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, "1e-9", "1e-6", 1],
               [2, .001, 0, 0, .001, 0, 0, 0, 0, 0, 0, 1, 2, "1e-9", "1e-6", 1]])
    write_csv(root / "bond_states.csv",
              ["point_i_global", "point_j_global", "active", "reference_length_m", "current_length_m"],
              [[1, 2, True, .001, .001]])
    write_csv(root / "broken_bonds.csv",
              ["point_i_global", "point_j_global", "reference_length_m", "current_length_m"], [])
    write_csv(root / "components.csv", ["point_id", "component_id", "point_mass_kg", "damage"],
              [[1, 1, "1e-6", 0], [2, 1, "1e-6", 0]])
    write_csv(root / "isolated_points.csv", ["point_id", "x_ref_m", "y_ref_m", "z_ref_m", "mass_kg"], [])
    return {"status": "COMPLETED", "n_points": 2, "n_steps": 1, "dt_s": 1e-6, "time_s": 1e-6,
            "requested_dt_mode": "EXPLICIT", "native_stability_bound_dt_s": 2e-6,
            "accepted_dt_s": 1e-6, "native_stability_check": "PASS_EXPLICIT_DT_WITHIN_NATIVE_BOUND",
            "native_stability_scope": "single-body fixture only; not a general OSB/contact stability guarantee",
            "max_damage": 0.0, "mean_damage": 0.0, "strain_energy_j": 4e-9, "kinetic_energy_j": 0.0,
            "point_mass_kg": 2e-6, "point_volume_sum_mm3": 2.0, "wall_time_s": 0.01,
            "intact_bond_connected_component_count": 1}


class ResultValidationTests(unittest.TestCase):
    def test_output_dictionary_describes_every_emitted_csv_column(self) -> None:
        for table_name, table in OUTPUT_DATA_DICTIONARY["tables"].items():
            if not table_name.endswith(".csv"):
                continue
            for column in table["columns"].values():
                self.assertIn("dtype", column)
                self.assertIn("unit", column)
                self.assertIn("frame", column)
                self.assertIn("role", column)
                self.assertIn("missing", column)

    def test_valid_endpoint_state_passes_integrity_checks(self) -> None:
        with temporary_directory() as temporary:
            root = Path(temporary)
            summary = make_valid_output(root)
            job = base_job()
            job["material"]["rho_kg_m3"] = 1000.0
            job["solver"].update(steps=1, dt_s=1e-6)
            checks = _validate_solver_outputs(root, summary, job)
        self.assertTrue(checks["point_ids_complete_unique"])
        self.assertTrue(checks["energy_integrals_match_point_state"])
        self.assertTrue(checks["component_partition_matches_graph"])
        self.assertTrue(checks["native_time_step_bound_validated"])

    def test_explicit_dt_above_native_bound_is_rejected(self) -> None:
        with temporary_directory() as temporary:
            root = Path(temporary)
            summary = make_valid_output(root)
            summary["native_stability_bound_dt_s"] = 0.5e-6
            job = base_job()
            job["material"]["rho_kg_m3"] = 1000.0
            job["solver"].update(steps=1, dt_s=1e-6)
            with self.assertRaisesRegex(ValidationError, "exceeds the recorded upstream native stability bound"):
                _validate_solver_outputs(root, summary, job)

    def test_nan_in_current_state_is_fail_closed(self) -> None:
        with temporary_directory() as temporary:
            root = Path(temporary)
            summary = make_valid_output(root)
            path = root / "points.csv"
            with path.open(encoding="utf-8") as stream:
                rows = list(csv.reader(stream))
            rows[1][4] = "NaN"
            with path.open("w", newline="", encoding="utf-8") as stream:
                csv.writer(stream).writerows(rows)
            job = base_job()
            job["material"]["rho_kg_m3"] = 1000.0
            job["solver"].update(steps=1, dt_s=1e-6)
            with self.assertRaisesRegex(ValidationError, "non-finite"):
                _validate_solver_outputs(root, summary, job)

    def test_component_map_cannot_hide_broken_bond_graph(self) -> None:
        with temporary_directory() as temporary:
            root = Path(temporary)
            summary = make_valid_output(root)
            path = root / "components.csv"
            with path.open(encoding="utf-8") as stream:
                rows = list(csv.reader(stream))
            rows[2][1] = "2"
            with path.open("w", newline="", encoding="utf-8") as stream:
                csv.writer(stream).writerows(rows)
            job = base_job()
            job["material"]["rho_kg_m3"] = 1000.0
            job["solver"].update(steps=1, dt_s=1e-6)
            with self.assertRaisesRegex(ValidationError, "labels disagree|partition disagrees"):
                _validate_solver_outputs(root, summary, job)


if __name__ == "__main__":
    unittest.main()
