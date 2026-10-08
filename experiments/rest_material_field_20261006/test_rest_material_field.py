"""Runnable focused checks: python -m unittest test_rest_material_field.py"""

import unittest

from focused_checks import run_focused_checks


class RestMaterialFieldContractTests(unittest.TestCase):
    def test_required_material_field_contract(self) -> None:
        receipt = run_focused_checks()
        self.assertEqual(receipt["status"], "PASS")
        self.assertEqual(receipt["assertion_groups"], 5)
        self.assertEqual(receipt["point_evaluation_attempts"], 31)


if __name__ == "__main__":
    unittest.main()
