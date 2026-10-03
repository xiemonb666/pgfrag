from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from pgfrag.errors import UnsupportedCapability, ValidationError
from pgfrag.jobs import safe_output_path, validate_job
from .helpers import base_job, temporary_directory


class JobValidationTests(unittest.TestCase):
    def test_solver_contract_distinguishes_asset_nu_from_backend_supported_nu(self) -> None:
        job = base_job()
        job["material"]["nu"] = -0.2
        with self.assertRaisesRegex(UnsupportedCapability, "supports only 0 < nu < 0.5"):
            validate_job(job)
        job["material"]["nu"] = 0.24
        self.assertEqual(validate_job(job)["material"]["nu"], 0.24)

    def test_bb_does_not_silently_coerce_unsupported_nu(self) -> None:
        job = base_job()
        job["material"].update(model="BB", nu=0.24)
        with self.assertRaisesRegex(UnsupportedCapability, "requires nu=0.25"):
            validate_job(job)

    def test_velocity_verlet_rejects_native_open_upper_safety_boundary(self) -> None:
        job = base_job()
        job["solver"]["safety_factor"] = 1.0
        with self.assertRaisesRegex(ValidationError, "0 < safety_factor < 1"):
            validate_job(job)

    def test_multithread_request_is_not_advertised_without_multichunk_export(self) -> None:
        job = base_job()
        job["resources"]["threads"] = 2
        with self.assertRaisesRegex(UnsupportedCapability, "one Julia thread/chunk"):
            validate_job(job)

    def test_fragment_count_targets_are_forbidden_at_any_depth(self) -> None:
        job = base_job()
        job["loading"]["target_fragments"] = 4
        with self.assertRaisesRegex(ValidationError, "may not request or target fragment counts"):
            validate_job(job)

    def test_lattice_resource_guard(self) -> None:
        job = base_job()
        job["initial_condition"].update(dimensions_mm=[20, 20, 20], dx_mm=0.5, horizon_mm=1.5)
        with self.assertRaisesRegex(ValidationError, "exceeds resources.max_points"):
            validate_job(job)

    def test_output_path_parent_traversal_is_rejected(self) -> None:
        job = base_job()
        job["output_dir"] = "../outside"
        with self.assertRaisesRegex(ValidationError, "parent traversal"):
            validate_job(job)

    def test_safe_output_path_stays_below_job_directory(self) -> None:
        with temporary_directory() as temporary:
            root = Path(temporary)
            job_path = root / "job.json"
            job_path.write_text("{}", encoding="utf-8")
            result = safe_output_path(job_path, "nested/run")
            self.assertEqual(result, root / "nested" / "run")


if __name__ == "__main__":
    unittest.main()
