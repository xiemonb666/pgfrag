"""Run only the small rigid-plane API checks; never starts the toy trajectory."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

from run_toy_diagnostic import FOCUSED_CHECK_DETAILS, run_self_checks


HERE = Path(__file__).resolve().parent
DEFAULT_OUT = HERE.parent / "focused_readonly_check"
SOURCE_NAMES = ("rigid_plane_contact.py", "run_toy_diagnostic.py", "run_focused_checks.py")
HISTORIC_ARTIFACTS = (
    "actual_trajectory_force_energy_impulse.csv",
    "comparison_errors.csv",
    "contact_trajectory_energy.png",
    "diagnostic_summary.json",
    "check_receipt.json",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    checked = run_self_checks()
    if "plane_point_and_normal_cannot_be_unlocked" not in checked:
        raise AssertionError("immutable plane-array check was not recorded")
    source_pins = {name: _sha256(HERE / name) for name in SOURCE_NAMES}
    artifact_root = HERE.parents[1]
    artifact_pins: dict[str, dict[str, Any]] = {}
    for name in HISTORIC_ARTIFACTS:
        path = artifact_root / name
        artifact_pins[name] = {
            "sha256_at_focused_check": _sha256(path),
            "present": path.is_file(),
            "history_note": (
                "comparison_errors.csv was post-run corrected for metric-column semantics; its pre-correction bytes were not preserved"
                if name == "comparison_errors.csv"
                else "existing bounded-run artifact; pinned in its current post-run state"
            ),
        }
    pins = {
        "schema": "rigid_plane_contact_focused_pins.v1",
        "focused_source_snapshot_dir": str(HERE),
        "focused_sources_sha256": source_pins,
        "prior_run_artifacts_current_sha256": artifact_pins,
    }
    historical = {
        "first_toy_run_source_snapshot_saved_at_run_time": False,
        "first_toy_run_source_sha256_captured_at_run_time": False,
        "first_toy_run_source_sha256": None,
        "first_toy_run_trajectory_csv_and_png_still_present": True,
        "first_toy_run_comparison_errors_original_bytes_preserved": False,
        "first_toy_run_comparison_errors_was_post_run_corrected": True,
        "no_trajectory_or_energy_csv_png_was_regenerated_by_this_focused_check": True,
    }
    receipt = {
        "schema": "rigid_plane_contact_focused_check.v1",
        "status": "PASS",
        "check_kind": "run_self_checks only; no simulate/main invocation",
        "passed_check_count": len(checked),
        "passed_checks": checked,
        "details": FOCUSED_CHECK_DETAILS,
        "historical_first_run_record": historical,
        "pin_file": "focused_source_and_prior_artifact_pins.json",
        "source_sha256": source_pins,
    }
    (outdir / "focused_source_and_prior_artifact_pins.json").write_text(
        json.dumps(pins, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    (outdir / "focused_check_receipt.json").write_text(
        json.dumps(receipt, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
