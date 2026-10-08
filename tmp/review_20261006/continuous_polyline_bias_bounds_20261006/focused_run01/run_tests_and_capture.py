"""Single approved unittest run with passive capture of returned metric records."""

from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path


SOURCE_DIR = Path(
    r"C:\Users\xiemo\Desktop\ShardForge\pgfrag\experiments\continuous_polyline_bias_bounds_20261006"
)
TEST_DIR = SOURCE_DIR / "tests"
RECORD_PATH = Path(__file__).with_name("metric_return_records.json")
sys.path.insert(0, str(SOURCE_DIR))

import continuous_polyline_bias_bounds as metric_module


original_metric = metric_module.symmetric_polyline_hausdorff_bounds
active_test = {"case_sequence": 0, "test_id": None}
result_sequence = 0
metric_returns: list[dict[str, object]] = []


def capture_return(*args: object, **kwargs: object) -> object:
    global result_sequence
    result = original_metric(*args, **kwargs)
    result_sequence += 1
    metric_returns.append(
        {
            "metric_result_sequence": result_sequence,
            "test_case_sequence": active_test["case_sequence"],
            "test_id": active_test["test_id"],
            "record": result.as_record(),
        }
    )
    return result


metric_module.symmetric_polyline_hausdorff_bounds = capture_return


class RecordingResult(unittest.TextTestResult):
    def startTest(self, test: unittest.case.TestCase) -> None:
        active_test["case_sequence"] += 1
        active_test["test_id"] = test.id()
        super().startTest(test)

    def stopTest(self, test: unittest.case.TestCase) -> None:
        super().stopTest(test)
        active_test["test_id"] = None


class RecordingRunner(unittest.TextTestRunner):
    resultclass = RecordingResult


suite = unittest.defaultTestLoader.discover(
    start_dir=str(TEST_DIR),
    pattern="test_continuous_polyline_bias_bounds.py",
)
result = RecordingRunner(stream=sys.stdout, verbosity=2).run(suite)
RECORD_PATH.write_text(
    json.dumps(
        {
            "capture_scope": "returns_from_the_single_unittest_run_only",
            "python_executable": sys.executable,
            "python_version": sys.version,
            "bytecode_writes_disabled": sys.dont_write_bytecode,
            "test_groups_run": result.testsRun,
            "run_successful": result.wasSuccessful(),
            "metric_return_count": len(metric_returns),
            "metric_returns": metric_returns,
        },
        indent=2,
        sort_keys=True,
    )
    + "\n",
    encoding="utf-8",
)
sys.exit(0 if result.wasSuccessful() else 1)
