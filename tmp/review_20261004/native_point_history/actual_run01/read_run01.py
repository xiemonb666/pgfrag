from __future__ import annotations

import hashlib
import json
from pathlib import Path

from pgfrag.native_point_history import NATIVE_ORIGIN_STATUS, read_native_point_history


ROOT = Path(r"C:\Users\xiemo\Desktop\ShardForge\pgfrag")
SNAPSHOT = ROOT / "tmp" / "review_20261004" / "pd_point_event_probe" / "snapshots" / "zero_load_control_run01"
POINTS = ROOT / "experiments" / "pd_point_event_probe_20261004" / "input" / "coupon_uniform_voxel_5cube_125.points.csv"
METADATA = ROOT / "experiments" / "pd_point_event_probe_20261004" / "input" / "coupon_uniform_voxel_5cube_125.metadata.toml"
OBSERVER_SOURCE = ROOT / "experiments" / "pd_point_event_probe_20261004" / "run_pd_point_event_probe.jl"
OUT = ROOT / "tmp" / "review_20261004" / "native_point_history" / "actual_run01"
REPORT_PATH = OUT / "adapter_report.json"
MAPPING_PATH = OUT / "bond_events_input_mapping.json"


def fingerprint(path: Path) -> dict[str, object]:
    try:
        data = path.read_bytes()
    except OSError as exc:
        return {"reference": str(path), "bytes": None, "sha256": None, "read_error": str(exc)}
    return {"reference": str(path), "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def main() -> int:
    if REPORT_PATH.exists() or MAPPING_PATH.exists():
        raise FileExistsError("refusing to overwrite a previous actual run01 adapter result")
    source_paths = (
        SNAPSHOT / "summary.toml",
        SNAPSHOT / "point_state_history.csv",
        SNAPSHOT / "directed_bond_state_history.csv",
        SNAPSHOT / "canonical_bond_state_history.csv",
        POINTS,
        METADATA,
        OBSERVER_SOURCE,
    )
    before = [fingerprint(path) for path in source_paths]
    result = read_native_point_history(
        SNAPSHOT,
        POINTS,
        METADATA,
        observer_source_path=OBSERVER_SOURCE,
    )
    after = [fingerprint(path) for path in source_paths]
    input_bytes_stable = before == after
    mapping = result.bond_events_mapping
    if mapping is not None:
        MAPPING_PATH.write_text(json.dumps(mapping, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    report = {
        "schema": "pgfrag.native_point_history.actual_run01_report.v1",
        "read_only_adapter_invocation": True,
        "solver_or_bond_events_detector_invoked": False,
        "native_origin_status": NATIVE_ORIGIN_STATUS,
        "physical_fragment_status": "UNKNOWN_NOT_DERIVED",
        "fracture_surface_status": "UNKNOWN_NOT_DERIVED",
        "energy_status": "UNKNOWN_NOT_DERIVED",
        "R2_status": "UNKNOWN_NOT_DERIVED",
        "reader_result": result.to_dict(),
        "mapping_file": MAPPING_PATH.name if mapping is not None else None,
        "input_bytes_stable_before_after": input_bytes_stable,
        "input_fingerprints_before": before,
        "input_fingerprints_after": after,
    }
    REPORT_PATH.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result.status,
        "native_origin_status": result.native_origin_status,
        "point_count": len(result.point_ids),
        "validated_frame_count": len(result.frames),
        "expected_final_step": result.expected_final_step,
        "initial_canonical_bond_count": len(result.initial_bonds),
        "mapping_available": mapping is not None,
        "reasons": list(result.reasons),
        "input_bytes_stable_before_after": input_bytes_stable,
        "report_path": str(REPORT_PATH),
        "mapping_path": str(MAPPING_PATH) if mapping is not None else None,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
