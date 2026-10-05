"""Build transport-only policies for the fixed two-box native graph fixture.

This is synthetic fixture plumbing. It does not read CAD, run pytest, or represent
an actual material-visibility classification. Its normal variant labels within-box
pairs IN_MATERIAL and cross-gap pairs OUT_OR_VOID. Its blocked variant additionally
forces one point to zero degree so the native fail-closed gate can be exercised.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from types import ModuleType
from typing import Any


ROOT = Path(__file__).resolve().parent
POINTS_CSV = ROOT / "fixedcase_points.csv"
COMPONENTS_CSV = ROOT / "fixedcase_components.csv"
SOURCE_MARKER = ROOT / "fixedcase_source_marker.txt"
HORIZON_M = 0.0018
ISOLATED_EXTERNAL_ID = 103
MAX_POINTS = 125
MAX_DIRECTED_PAIRS = 20_000
MAX_OUTPUT_BYTES = 20 * 1024 * 1024


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical_json(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"),
                        ensure_ascii=False, allow_nan=False) + "\n").encode("utf-8")


def read_case() -> tuple[list[dict[str, Any]], dict[int, str]]:
    with POINTS_CSV.open("r", encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows or len(rows) > MAX_POINTS:
        raise ValueError("fixed fixture point count is outside its source budget")
    if list(rows[0]) != ["point_id", "x_mm", "y_mm", "z_mm"]:
        raise ValueError("fixed fixture CSV header changed")
    points: list[dict[str, Any]] = []
    seen: set[int] = set()
    for native_index, row in enumerate(rows, start=1):
        point_id = int(row["point_id"])
        if point_id in seen or point_id < 0:
            raise ValueError("fixed fixture point IDs must be unique non-negative integers")
        seen.add(point_id)
        xyz = [float(row[name]) for name in ("x_mm", "y_mm", "z_mm")]
        if not all(math.isfinite(value) for value in xyz):
            raise ValueError("fixed fixture contains a non-finite coordinate")
        points.append({"point_id": point_id, "native_index": native_index,
                       "position_mm": xyz})
    with COMPONENTS_CSV.open("r", encoding="utf-8", newline="") as stream:
        components = {int(row["point_id"]): row["component"]
                      for row in csv.DictReader(stream)}
    if set(components) != seen or set(components.values()) != {"left", "right"}:
        raise ValueError("fixed fixture component labels must cover the exact point IDs")
    return points, components


def hash_point_order(points: list[dict[str, Any]]) -> str:
    body = "pgfrag.point-ids.v1\n" + "".join(
        f"{row['point_id']}\n" for row in points)
    return sha256(body.encode("ascii"))


def hash_pair_map(rows: list[dict[str, Any]]) -> str:
    body = "pgfrag.native-pair-map.v1\n" + "".join(
        f"{row['pair'][0]},{row['pair'][1]},{row['classification']}\n"
        for row in rows)
    return sha256(body.encode("ascii"))


def hash_included(pairs: list[list[int]]) -> str:
    body = "pgfrag.native-included-mask.v1\n" + "".join(
        f"{pair[0]},{pair[1]}\n" for pair in pairs)
    return sha256(body.encode("ascii"))


def create_bundle(bundle_dir: Path, *, isolate_id: int | None) -> dict[str, Any]:
    bundle_dir.mkdir(parents=False, exist_ok=False)
    points, components = read_case()
    point_by_id = {row["point_id"]: row for row in points}
    if isolate_id is not None and isolate_id not in point_by_id:
        raise ValueError("zero-degree probe ID is absent from fixture points")

    decisions: list[dict[str, Any]] = []
    for source_index, source in enumerate(points):
        for target in points[source_index + 1:]:
            source_xyz = source["position_mm"]
            target_xyz = target["position_mm"]
            distance2_m = sum(((source_xyz[axis] - target_xyz[axis]) * 0.001) ** 2
                              for axis in range(3))
            if distance2_m > HORIZON_M ** 2:
                continue
            pair = sorted((source["point_id"], target["point_id"]))
            if components[pair[0]] == components[pair[1]]:
                state = "IN_MATERIAL"
                reason = "fixture component membership: same closed box"
            else:
                state = "OUT_OR_VOID"
                reason = "fixture void: segment spans the explicit empty x-gap"
            if isolate_id is not None and isolate_id in pair and state == "IN_MATERIAL":
                state = "OUT_OR_VOID"
                reason = "NEGATIVE_TEST_ONLY: forced exclusion to probe zero-degree rejection"
            decisions.append({"pair": pair, "classification": state, "reason": reason})
    decisions.sort(key=lambda row: tuple(row["pair"]))
    if len(decisions) * 2 > MAX_DIRECTED_PAIRS:
        raise ValueError("fixed fixture directed candidate count exceeds its source budget")

    groups = {
        state: [row["pair"] for row in decisions if row["classification"] == state]
        for state in ("IN_MATERIAL", "OUT_OR_VOID", "UNKNOWN")
    }
    degrees = {row["point_id"]: 0 for row in points}
    for left, right in groups["IN_MATERIAL"]:
        degrees[left] += 1
        degrees[right] += 1
    zero_degree_ids = [row["point_id"] for row in points
                       if degrees[row["point_id"]] == 0]
    if isolate_id is None:
        eligibility = "CANDIDATE_MASK_POSITIVE_DEGREE_NATIVE_UNVERIFIED"
    elif zero_degree_ids == [isolate_id]:
        eligibility = "BLOCKED_ZERO_DEGREE"
    else:
        raise ValueError("forced zero-degree fixture variant did not isolate exactly one point")

    result_path = bundle_dir / "fixture_classifier_result.json"
    sidecar_path = bundle_dir / "fixture_interval_evidence.json"
    result_bytes = canonical_json({
        "label": "GEOMETRY_FILTER_FIXTURE_ONLY",
        "classification_method": "fixed_box_component_membership",
        "candidate_pair_count": len(decisions),
        "included_pair_count": len(groups["IN_MATERIAL"]),
        "geometry_excluded_pair_count": len(groups["OUT_OR_VOID"]),
        "unknown_pair_count": len(groups["UNKNOWN"]),
        "zero_degree_external_point_ids": zero_degree_ids,
        "negative_variant_is_not_geometry_truth": isolate_id is not None,
    })
    cross_gap_rows = [row for row in decisions
                      if components[row["pair"][0]] != components[row["pair"][1]]]
    sidecar_bytes = canonical_json({
        "label": "GEOMETRY_FILTER_FIXTURE_ONLY",
        "meaning": "explicit point-set gap metadata; not CAD ray-hit evidence",
        "void_x_interval_mm": [1.0, 2.5],
        "cross_gap_candidate_pairs": [row["pair"] for row in cross_gap_rows],
        "forced_zero_degree_pair_ids": isolate_id,
    })
    result_path.write_bytes(result_bytes)
    sidecar_path.write_bytes(sidecar_bytes)

    points_bytes = POINTS_CSV.read_bytes()
    source_bytes = SOURCE_MARKER.read_bytes()
    classifier_bytes = Path(__file__).read_bytes()
    raw_candidates = [{
        "pair": row["pair"],
        "classification": row["classification"],
        "reason": row["reason"],
    } for row in decisions]
    included_pairs = groups["IN_MATERIAL"]
    audit = {
        "schema": "pgfrag.material-visibility-native-mask.v1",
        "static_contract_status": "VALID",
        "native_mask_eligibility": eligibility,
        "native_gate": {
            "status": "UNVERIFIED",
            "raw_candidate_set_reproduced": False,
            "native_graph_readback": "NOT_RUN",
            "initialization_authorized": False,
            "candidate_pair_count_nonzero": bool(decisions),
            "included_pair_count_nonzero": bool(included_pairs),
            "zero_degree_point_ids": zero_degree_ids,
        },
        "source": {
            "source_cad_path": str(SOURCE_MARKER.resolve()),
            "source_cad_sha256": sha256(source_bytes),
            "classifier_source_sha256": sha256(classifier_bytes),
            "source_hash_status": "MATCHED_ACTUAL_FILE_BYTES",
            "source_status": "READY",
            "source_accuracy_mm": 0.0,
            "source_accuracy_basis": (
                "GEOMETRY_FILTER_FIXTURE_ONLY schema placeholder; zero is not a CAD accuracy claim"
            ),
        },
        "coordinates": {
            "input_unit": "mm",
            "native_unit": "m",
            "m_per_input_unit": 0.001,
            "reference_frame": "fixedcase_reference_mm",
            "native_frame": "Peridynamics_Body_position_m",
            "native_frame_from_reference_frame_m": [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
        },
        "points": {
            "payload_path": str(POINTS_CSV.resolve()),
            "payload_sha256": sha256(points_bytes),
            "order": "exact_payload_row_order",
            "point_ids_order_sha256": hash_point_order(points),
            "rows": points,
        },
        "candidate_generation": {
            "metric": "euclidean_reference_distance",
            "horizon_m": HORIZON_M,
            "candidate_rule": "exact_pinned_Peridynamics_find_bonds_convention",
            "scope": "one_body_one_uniform_horizon",
        },
        "raw_euclidean_candidate_pairs": raw_candidates,
        "counts": {
            "raw_euclidean_candidate_pairs": len(decisions),
            "initial_material_bond_pairs": len(included_pairs),
            "geometry_excluded_pairs": len(groups["OUT_OR_VOID"]),
            "unknown_pairs": len(groups["UNKNOWN"]),
        },
        "pair_map_sha256": hash_pair_map(raw_candidates),
        "mask": {
            "included_pairs": included_pairs,
            "geometry_excluded_pairs": groups["OUT_OR_VOID"],
            "unknown_pairs": groups["UNKNOWN"],
            "mask_sha256": hash_included(included_pairs),
        },
        "evidence": {
            "classifier_result_path": str(result_path.resolve()),
            "classifier_result_sha256": sha256(result_bytes),
            "classifier_source_path": str(Path(__file__).resolve()),
            "classifier_source_sha256": sha256(classifier_bytes),
            "raw_hit_interval_sidecar_path": str(sidecar_path.resolve()),
            "raw_hit_interval_sidecar_sha256": sha256(sidecar_bytes),
            "sidecar_selection": {
                "test_nodeid": "GEOMETRY_FILTER_FIXTURE_ONLY::fixed_two_box_rule",
                "call_index": 1,
            },
            "visibility_algorithm": "fixed_box_component_membership_fixture_v1",
        },
        "initial_physics_events": {
            "evidence_status": "PLANNED_NOT_OBSERVED_NATIVE_NOT_RUN",
            "planned_physical_broken_bond_pairs": 0,
            "geometry_exclusions_are_physics_events": False,
        },
    }
    audit_path = bundle_dir / "audit_contract.json"
    audit_bytes = canonical_json(audit)
    audit_path.write_bytes(audit_bytes)

    compiler_path = ROOT.parent / "material_visibility_native_bridge_20261005" / "native_mask_transport.py"
    compiler = ModuleType("pgfrag_fixture_native_mask_transport")
    compiler.__file__ = str(compiler_path.resolve())
    compiler_source = compiler_path.read_text(encoding="utf-8")
    exec(compile(compiler_source, str(compiler_path), "exec"), compiler.__dict__)
    receipt = compiler.compile_native_mask_transport(audit_path, bundle_dir / "native_mask.toml")
    receipt["fixture_label"] = "GEOMETRY_FILTER_FIXTURE_ONLY"
    receipt["source_status_is_schema_token_only"] = True
    receipt_path = bundle_dir / "external_pins_fixture_only.json"
    receipt_path.write_bytes(canonical_json(receipt))
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_root", type=Path,
                        help="new disposable output directory; existing paths are rejected")
    args = parser.parse_args()
    output_root = args.output_root.expanduser().resolve(strict=False)
    scratch_root = (ROOT.parents[1] / "tmp" / "initial_bond_policy_body_pilot_20261005").resolve(strict=False)
    try:
        output_root.relative_to(scratch_root)
    except ValueError:
        parser.error(f"fixture output must stay under the dedicated scratch root: {scratch_root}")
    if output_root == scratch_root:
        parser.error("choose a new child directory beneath the dedicated scratch root")
    if output_root.exists() or output_root.is_symlink():
        parser.error(f"refusing to reuse existing output directory: {output_root}")
    output_root.mkdir(parents=True, exist_ok=False)
    positive = create_bundle(output_root / "positive_mask", isolate_id=None)
    zero_degree = create_bundle(output_root / "zero_degree_block", isolate_id=ISOLATED_EXTERNAL_ID)
    output_bytes = sum(path.stat().st_size for path in output_root.rglob("*") if path.is_file())
    if output_bytes > MAX_OUTPUT_BYTES:
        parser.error("fixture output exceeds its 20 MiB artifact budget")
    print(canonical_json({
        "fixture_label": "GEOMETRY_FILTER_FIXTURE_ONLY",
        "positive_mask": positive,
        "zero_degree_block": zero_degree,
        "native_body_and_graph": "NOT_RUN",
        "solver_force_cad": "NOT_RUN",
    }).decode("utf-8"), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
