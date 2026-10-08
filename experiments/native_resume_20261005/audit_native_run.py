#!/usr/bin/env python3
"""Read-only audit of the completed 2026-10-04 native affine threshold run."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import tomllib
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "pgfrag"
EXPERIMENT = BASE / "experiments" / "pd_loaded_diagnostics_20261004"
OLD_REVIEW = BASE / "tmp" / "review_20261004" / "pd_loaded_diagnostics"
RUN_ID = "affine_threshold_b_run01"
RECEIPTS = OLD_REVIEW / "receipts"
SNAPSHOT = OLD_REVIEW / "snapshots" / RUN_ID
REVIEW_OUTPUT = BASE / "tmp" / "review_20261005" / "native_resume"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def bool_csv(value: str) -> bool:
    lowered = value.strip().lower()
    require(lowered in {"true", "false"}, f"invalid bool cell: {value!r}")
    return lowered == "true"


def close(a: float, b: float, *, rtol: float = 2e-12, atol: float = 2e-15) -> bool:
    return math.isclose(a, b, rel_tol=rtol, abs_tol=atol)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    REVIEW_OUTPUT.mkdir(parents=True, exist_ok=True)
    receipt_path = RECEIPTS / f"{RUN_ID}.launch_receipt.json"
    preflight_path = RECEIPTS / f"{RUN_ID}.preflight.json"
    stdout_path = RECEIPTS / f"{RUN_ID}.stdout.log"
    stderr_path = RECEIPTS / f"{RUN_ID}.stderr.log"
    manifest_path = EXPERIMENT / "manifest.json"
    a_figure_path = BASE / "tmp" / "review_20261004" / "pd_loaded_diagnostics" / "reports" / "loaded_elastic_a8_science.png"
    b_figure_path = REVIEW_OUTPUT / "native_b_affine_threshold_audit.png"
    driver_path = EXPERIMENT / "run_affine_critical_stretch_probe.jl"
    case_path = EXPERIMENT / "case_b_affine_threshold.toml"
    common_path = EXPERIMENT / "common_loaded_case.toml"
    points_path = BASE / "experiments" / "pd_point_event_probe_20261004" / "input" / "coupon_uniform_voxel_5cube_125.points.csv"
    metadata_path = BASE / "experiments" / "pd_point_event_probe_20261004" / "input" / "coupon_uniform_voxel_5cube_125.metadata.toml"
    directed_path = SNAPSHOT / "affine_native_bond_states.csv"
    pair_path = SNAPSHOT / "affine_canonical_pair_states.csv"
    summary_path = SNAPSHOT / "summary.toml"

    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    preflight = json.loads(preflight_path.read_text(encoding="utf-8"))
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    summary = tomllib.loads(summary_path.read_text(encoding="utf-8"))
    case = tomllib.loads(case_path.read_text(encoding="utf-8"))
    common = tomllib.loads(common_path.read_text(encoding="utf-8"))
    metadata = tomllib.loads(metadata_path.read_text(encoding="utf-8"))

    require(receipt["status"] == "COMPLETED", "B receipt status is not completed")
    require(receipt["wrapper_exit_code"] == 0 and receipt["child_exit_code"] == 0,
            "B wrapper/child exit was not zero")
    require(receipt["child_exit_observed"] and receipt["wait_completed"] and not receipt["timed_out"],
            "B child exit was not observed cleanly")
    require(receipt["launch_error"] is None, "B receipt contains a launch error")
    require(preflight["status"] == "PREFLIGHT_PASS_SINGLE_RUN_ARMED", "unexpected B preflight")
    require(preflight["driver_sha256"] == sha256(driver_path), "driver hash differs from preflight")
    require(receipt["manifest_sha256_before"] == sha256(manifest_path), "manifest hash differs from launch receipt")
    require(receipt["preflight_sha256"] == sha256(preflight_path), "preflight hash differs from launch receipt")

    pins_before = receipt["source_pins_before"]
    pins_after = receipt["source_pins_after"]
    require(len(pins_before) == len(pins_after) == 26, "expected exactly 26 source pins")
    require(receipt["source_pin_count_before"] == receipt["source_pin_count_after"] == 26,
            "receipt source pin counts disagree")
    current_pin_audit = []
    for before, after in zip(pins_before, pins_after, strict=True):
        require(before["path"] == after["path"], "source pin order/path changed between receipt snapshots")
        path = ROOT / Path(before["path"])
        actual_sha = sha256(path)
        actual_bytes = path.stat().st_size
        match = actual_sha == before["sha256"] == after["sha256"] and actual_bytes == before["bytes"] == after["bytes"]
        require(match and after.get("matches_expected") is True,
                f"source pin mismatch: {before['path']}")
        current_pin_audit.append({
            "path": before["path"],
            "bytes": actual_bytes,
            "expected_sha256": before["sha256"],
            "actual_sha256": actual_sha,
            "matches_receipt_before_and_after": match,
        })

    runtime_path = ROOT / "third_party" / "pgfrag-runtime" / "extract" / "julia-1.10.12" / "bin" / "julia.exe"
    require(sha256(runtime_path) == receipt["runtime_sha256"], "Julia executable hash differs from receipt")
    require(sha256(stdout_path) == receipt["stdout"]["sha256"] and stdout_path.stat().st_size == receipt["stdout"]["bytes"],
            "stdout hash/length differs from receipt")
    require(sha256(stderr_path) == receipt["stderr"]["sha256"] and stderr_path.stat().st_size == receipt["stderr"]["bytes"],
            "stderr hash/length differs from receipt")
    expected_artifacts = {Path(item["path"]).name: item for item in receipt["output_artifacts"]}
    actual_artifacts = {}
    for path in (directed_path, pair_path, summary_path):
        record = expected_artifacts[path.name]
        actual_sha = sha256(path)
        actual_bytes = path.stat().st_size
        require(actual_sha == record["sha256"] and actual_bytes == record["bytes"],
                f"output artifact differs from launch receipt: {path.name}")
        actual_artifacts[path.name] = {"path": str(path), "bytes": actual_bytes, "sha256": actual_sha}

    require(sha256(points_path) == metadata["point_data_sha256"], "fixture SHA differs from metadata")
    require(summary["point_csv_sha256"] == metadata["point_data_sha256"], "fixture SHA differs from B summary")
    require(summary["point_count"] == metadata["point_count"] == 125, "unexpected point count")
    xyz: dict[int, tuple[float, float, float]] = {}
    volumes = []
    with points_path.open("r", encoding="utf-8", newline="") as stream:
        for point in csv.DictReader(stream):
            point_id = int(point["point_id"])
            xyz[point_id] = (float(point["x_m"]), float(point["y_m"]), float(point["z_m"]))
            volumes.append(float(point["volume_m3"]))
    require(len(xyz) == 125 and len(volumes) == 125, "point file row count mismatch")
    axis_values = [sorted({position[axis] for position in xyz.values()}) for axis in range(3)]
    require(all(len(axis) == 5 for axis in axis_values), "fixture is not a 5x5x5 cell-center lattice")
    require(all(close(right - left, 0.001, atol=2e-15) for axis in axis_values for left, right in zip(axis, axis[1:])),
            "fixture point spacing is not 1 mm")
    require(all(close(axis[0], 0.0005, atol=2e-15) and close(axis[-1], 0.0045, atol=2e-15)
                for axis in axis_values), "fixture is not centered at 0.5..4.5 mm inside a 5 mm domain")
    require(close(sum(volumes), 1.25e-7, rtol=2e-15, atol=1e-22), "fixture volume sum mismatch")

    E = float(common["material"]["E_pa"])
    nu = float(common["material"]["nu"])
    Gc = float(common["material"]["Gc_j_m2"])
    delta = float(common["geometry"]["horizon_m"])
    K = E / (3.0 * (1.0 - 2.0 * nu))
    expected_epsilon_c = math.sqrt(5.0 * Gc / (9.0 * K * delta))
    observed_epsilon_c = float(summary["native_epsilon_c"])
    require(close(observed_epsilon_c, expected_epsilon_c, rtol=2e-14, atol=1e-16),
            "native epsilon_c does not match source formula for frozen inputs")
    expected_graph = case["expected_graph"]
    directed_expected = int(expected_graph["directed_bonds"])
    canonical_expected = int(expected_graph["canonical_bonds"])
    require(summary["native_directed_bonds"] == directed_expected == 6124, "directed graph count mismatch")
    require(summary["native_canonical_bonds"] == canonical_expected == 3062, "canonical graph count mismatch")

    directed_rows = read_csv(directed_path)
    canonical_rows = read_csv(pair_path)
    require(len(directed_rows) == 3 * directed_expected, "directed CSV total rows mismatch")
    require(len(canonical_rows) == 3 * canonical_expected, "canonical CSV total rows mismatch")
    state_specs = case["states"]
    expected_states = {str(item["state_id"]): float(item["epsilon_c_multiplier"]) for item in state_specs}
    require(expected_states == {"rest": 0.0, "below": 0.5, "above": 1.5}, "frozen state multipliers changed")
    directed_by_state_pair: dict[str, dict[tuple[int, int], list[dict[str, object]]]] = defaultdict(lambda: defaultdict(list))
    directed_keys = set()
    directed_state_metrics = {}
    for raw in directed_rows:
        state = raw["state_id"]
        require(state in expected_states, f"unexpected directed state {state}")
        source = int(raw["source_point_id"])
        neighbor = int(raw["neighbor_point_id"])
        canonical_i = int(raw["canonical_i_id"])
        canonical_j = int(raw["canonical_j_id"])
        require(source in xyz and neighbor in xyz and source != neighbor, "invalid directed endpoints")
        require((canonical_i, canonical_j) == (min(source, neighbor), max(source, neighbor)),
                "directed row canonical ids do not match endpoints")
        key = (state, source, neighbor)
        require(key not in directed_keys, f"duplicate directed row {key}")
        directed_keys.add(key)
        L = float(raw["reference_length_m"])
        l = float(raw["current_length_m"])
        stretch = float(raw["actual_stretch"])
        eps = float(raw["native_epsilon_c"])
        multiplier = float(raw["eta_over_native_epsilon_c"])
        eta = float(raw["eta"])
        before = bool_csv(raw["active_before"])
        after = bool_csv(raw["active_after"])
        fail_permit = bool_csv(raw["fail_permit"])
        transition = bool_csv(raw["transition"])
        require(all(math.isfinite(v) for v in (L, l, stretch, eps, eta)) and L > 0 and l > 0 and eps > 0,
                "non-finite or non-positive bond data")
        require(close(multiplier, expected_states[state], rtol=0, atol=0), "state multiplier mismatch")
        require(close(eta, multiplier * eps, rtol=5e-15, atol=1e-18), "eta is not multiplier times native epsilon_c")
        require(close(eps, observed_epsilon_c, rtol=1e-12, atol=1e-15), "per-bond epsilon_c varies")
        X_i = xyz[source]
        X_j = xyz[neighbor]
        dX = tuple(X_j[axis] - X_i[axis] for axis in range(3))
        expected_L = math.sqrt(sum(value * value for value in dX))
        expected_l = math.sqrt(dX[0] ** 2 + dX[1] ** 2 + ((1.0 + eta) * dX[2]) ** 2)
        expected_stretch = (expected_l - expected_L) / expected_L
        require(close(L, expected_L, rtol=2e-12, atol=1e-15), "reference bond length differs from point geometry")
        require(close(l, expected_l, rtol=2e-12, atol=1e-15), "native current bond length differs from prescribed affine geometry")
        require(close(stretch, (l - L) / L, rtol=1e-11, atol=2e-15), "stretch is inconsistent with recorded lengths")
        require(close(stretch, expected_stretch, rtol=2e-10, atol=3e-15), "stretch differs from geometric affine prediction")
        expected_after = before and not (stretch > eps and fail_permit)
        require(after == expected_after, "native active_after violates threshold/fail-permit condition")
        require(transition == (before != after), "transition flag mismatch")
        require(before, "fresh-state branch did not start with active bond")
        directed_by_state_pair[state][(canonical_i, canonical_j)].append({
            "source": source,
            "neighbor": neighbor,
            "L": L,
            "l": l,
            "stretch": stretch,
            "epsilon_c": eps,
            "fail_permit": fail_permit,
            "before": before,
            "after": after,
            "transition": transition,
        })
    for state, pairs in directed_by_state_pair.items():
        require(len(pairs) == canonical_expected, f"{state}: canonical pair count differs")
        directed_count = sum(len(directions) for directions in pairs.values())
        require(directed_count == directed_expected, f"{state}: directed row count differs")
        over_threshold = sum(1 for directions in pairs.values() for row in directions[:1]
                             if row["stretch"] > row["epsilon_c"] and row["fail_permit"])
        broken_canonical = sum(1 for directions in pairs.values() for row in directions[:1] if not row["after"])
        directed_transitions = sum(1 for directions in pairs.values() for row in directions if row["transition"])
        unique_pair_rows = [directions[0] for directions in pairs.values()]
        broken_ratios = [row["stretch"] / row["epsilon_c"] for row in unique_pair_rows if not row["after"]]
        retained_ratios = [row["stretch"] / row["epsilon_c"] for row in unique_pair_rows if row["after"]]
        fail_permit_true_directed = sum(int(bool(row["fail_permit"])) for directions in pairs.values() for row in directions)
        fail_permit_false_directed = directed_count - fail_permit_true_directed
        require(all(len(directions) == 2 for directions in pairs.values()), f"{state}: pair missing reciprocal row")
        reciprocity_differences = 0
        for pair, directions in pairs.items():
            first, second = directions
            require((first["source"], first["neighbor"]) == (second["neighbor"], second["source"]),
                    f"{state} {pair}: directions are not reciprocal")
            for field in ("L", "l", "stretch", "epsilon_c"):
                require(close(float(first[field]), float(second[field]), rtol=1e-12, atol=2e-15),
                        f"{state} {pair}: reciprocal {field} differs")
            for field in ("fail_permit", "before", "after", "transition"):
                require(first[field] == second[field], f"{state} {pair}: reciprocal {field} differs")
            reciprocity_differences += int(first["stretch"] != second["stretch"])
        require(reciprocity_differences == 0, f"{state}: reciprocal stretch records differ")
        require(over_threshold == broken_canonical, f"{state}: over-threshold pair and broken pair counts differ")
        if broken_ratios:
            require(min(broken_ratios) > 1.0, f"{state}: a broken pair does not exceed native epsilon_c")
        if retained_ratios:
            require(max(retained_ratios) <= 1.0 + 2e-10,
                    f"{state}: a fail-permitted retained pair exceeds native epsilon_c")
        directed_state_metrics[state] = {
            "directed_rows": directed_count,
            "canonical_pairs": len(pairs),
            "fail_permit_true_directed_rows": fail_permit_true_directed,
            "fail_permit_false_directed_rows": fail_permit_false_directed,
            "over_threshold_canonical_pairs": over_threshold,
            "broken_canonical_pairs": broken_canonical,
            "directed_transition_count": directed_transitions,
            "reciprocal_mismatches": reciprocity_differences,
            "min_broken_stretch_over_epsilon_c": min(broken_ratios) if broken_ratios else None,
            "max_retained_stretch_over_epsilon_c": max(retained_ratios) if retained_ratios else None,
            "min_stretch_over_epsilon_c": min(row["stretch"] / row["epsilon_c"] for directions in pairs.values() for row in directions[:1]),
            "max_stretch_over_epsilon_c": max(row["stretch"] / row["epsilon_c"] for directions in pairs.values() for row in directions[:1]),
        }

    pair_by_state = {}
    pair_keys = set()
    for raw in canonical_rows:
        state = raw["state_id"]
        require(state in expected_states, f"unexpected canonical state {state}")
        pair = (int(raw["canonical_i_id"]), int(raw["canonical_j_id"]))
        key = (state, pair)
        require(key not in pair_keys, f"duplicate canonical row {key}")
        pair_keys.add(key)
        require(pair in directed_by_state_pair[state], f"canonical row missing directed pair {key}")
        raw_directions = directed_by_state_pair[state][pair]
        first = raw_directions[0]
        for csv_name, internal_name in (("reference_length_m", "L"), ("current_length_m", "l"),
                                        ("actual_stretch", "stretch"), ("native_epsilon_c", "epsilon_c")):
            require(close(float(raw[csv_name]), float(first[internal_name]), rtol=1e-12, atol=2e-15),
                    f"canonical row value differs from directions: {key} {csv_name}")
        for csv_name, internal_name in (("fail_permit", "fail_permit"), ("active_before", "before"),
                                        ("active_after", "after"), ("transition", "transition")):
            require(bool_csv(raw[csv_name]) == first[internal_name], f"canonical row state differs from directions: {key} {csv_name}")
        pair_by_state.setdefault(state, {})[pair] = raw
    require(all(len(pair_by_state.get(state, {})) == canonical_expected for state in expected_states),
            "canonical CSV does not have one row for every pair/state")

    state_summary = {}
    for index, (state, expected_multiplier) in enumerate(expected_states.items(), start=1):
        row = summary[f"state_{index}"]
        metric = directed_state_metrics[state]
        require(row["state_id"] == state, "summary state ordering differs from frozen case")
        require(row["directed_state_rows"] == metric["directed_rows"], f"{state}: summary directed count differs")
        require(row["canonical_pairs"] == metric["canonical_pairs"], f"{state}: summary canonical count differs")
        require(row["directed_transition_count"] == metric["directed_transition_count"], f"{state}: summary directed transitions differ")
        require(row["canonical_transition_count"] == metric["broken_canonical_pairs"], f"{state}: summary canonical transitions differ")
        require(row["reciprocal_status"] == "PASS", f"{state}: summary reciprocity not PASS")
        state_summary[state] = {**metric, "summary_directed_transition_count": row["directed_transition_count"],
                                "summary_canonical_transition_count": row["canonical_transition_count"],
                                "reciprocal_status": row["reciprocal_status"],
                                "eta_over_epsilon_c": expected_multiplier}
    require(state_summary["rest"]["broken_canonical_pairs"] == 0, "rest control broke a bond")
    require(state_summary["below"]["broken_canonical_pairs"] == 0, "below-threshold control broke a bond")
    require(state_summary["above"]["broken_canonical_pairs"] > 0, "above-threshold control did not break a bond")
    require(summary["status"] == "PASS_NATIVE_BRANCH_DIAGNOSTIC_ONLY", "B summary did not pass its narrow gate")
    require(summary["native_time_integration"] is False and summary["manual_bond_activity_edits"] is False,
            "B summary does not confirm the frozen execution contract")

    h = axis_values[0][1] - axis_values[0][0]
    domain_side = max(axis_values[0]) - min(axis_values[0]) + h
    maximum_nearest_surface_clearance = domain_side / 2.0
    require(delta > maximum_nearest_surface_clearance,
            "this frozen geometry unexpectedly contains a full-horizon interior point")
    require(a_figure_path.is_file() and b_figure_path.is_file(), "expected A/B visual artifact is missing")
    a_figure_expected = manifest["A"]["figure"]
    require(sha256(a_figure_path) == a_figure_expected["sha256"], "A visual artifact differs from frozen manifest")

    audit = {
        "schema": "pgfrag.native_resume.audit.v1",
        "status": "PASS_READ_ONLY_NATIVE_THRESHOLD_AUDIT",
        "run_id": RUN_ID,
        "run_execution": {
            "receipt_path": str(receipt_path),
            "receipt_status": receipt["status"],
            "wrapper_exit_code": receipt["wrapper_exit_code"],
            "child_exit_code": receipt["child_exit_code"],
            "timed_out": receipt["timed_out"],
            "elapsed_seconds": receipt["elapsed_seconds"],
            "started_utc": receipt["started_utc"],
            "finished_utc": receipt["finished_utc"],
            "execution_disposition": "existing completed run audited; no solver rerun",
            "native_time_integration": False,
            "manual_bond_activity_edits": False,
        },
        "frozen_sources": {
            "source_pin_count": len(current_pin_audit),
            "all_26_current_hashes_match_receipt_before_and_after": all(item["matches_receipt_before_and_after"] for item in current_pin_audit),
            "pins": current_pin_audit,
            "driver_sha256": sha256(driver_path),
            "runtime_sha256": sha256(runtime_path),
            "manifest_sha256": sha256(manifest_path),
            "manifest_still_says": manifest["status"],
        },
        "raw_evidence": {
            "directed_csv": actual_artifacts[directed_path.name],
            "canonical_csv": actual_artifacts[pair_path.name],
            "summary_toml": actual_artifacts[summary_path.name],
            "a_loaded_observer_png": {"path": str(a_figure_path), "bytes": a_figure_path.stat().st_size, "sha256": sha256(a_figure_path)},
            "b_threshold_audit_png": {"path": str(b_figure_path), "bytes": b_figure_path.stat().st_size, "sha256": sha256(b_figure_path)},
            "stdout": {"path": str(stdout_path), "bytes": stdout_path.stat().st_size, "sha256": sha256(stdout_path)},
            "stderr": {"path": str(stderr_path), "bytes": stderr_path.stat().st_size, "sha256": sha256(stderr_path)},
            "fixture_points_csv": {"path": str(points_path), "bytes": points_path.stat().st_size, "sha256": sha256(points_path)},
            "fixture_metadata_toml": {"path": str(metadata_path), "bytes": metadata_path.stat().st_size, "sha256": sha256(metadata_path)},
            "launch_receipt_sha256": sha256(receipt_path),
            "preflight_receipt_sha256": sha256(preflight_path),
        },
        "fixture_and_principle_checks": {
            "point_count": len(xyz),
            "lattice_axis_point_counts": [len(axis) for axis in axis_values],
            "uniform_spacing_m": 0.001,
            "total_volume_m3": sum(volumes),
            "horizon_to_spacing_ratio": delta / h,
            "domain_side_m": domain_side,
            "horizon_to_domain_side_ratio": delta / domain_side,
            "max_point_clearance_to_nearest_external_surface_m": maximum_nearest_surface_clearance,
            "has_any_full_horizon_interior_point": False,
            "surface_correction": "NoCorrection default from OSBMaterial() in frozen native driver",
            "material": {"E_pa": E, "nu": nu, "Gc_j_m2": Gc, "rho_kg_m3": float(common["material"]["rho_kg_m3"]),
                         "horizon_m": delta, "K_pa_from_formula": K},
            "epsilon_c_formula": "sqrt(5 Gc / (9 K delta))",
            "epsilon_c_formula_result": expected_epsilon_c,
            "epsilon_c_native_reported": observed_epsilon_c,
            "directed_bonds_per_state": directed_expected,
            "canonical_bonds_per_state": canonical_expected,
            "states": state_summary,
            "per_directed_row_checks": [
                "finite positive reference/current length and finite stretch/epsilon_c",
                "reference length equals point-coordinate geometry",
                "current length equals x'=x, y'=y, z'=(1+eta)z prescribed affine geometry",
                "recorded stretch equals (current-reference)/reference and independent geometry calculation",
                "active_after equals active_before AND NOT(stretch > epsilon_c AND fail_permit)",
                "each directed row has exactly one reciprocal row with matching native measures and states",
            ],
            "per_canonical_row_checks": ["one canonical row per unique pair per state", "canonical values/states agree with the two directed rows"],
        },
        "physical_scope": {
            "qualification": summary["physical_qualification"],
            "energy": summary["energy_ledger"],
            "contact": summary["contact"],
            "free_fracture_or_fragment_claim": summary["free_fracture_or_fragment_claim"],
            "interpretation": "Native CriticalStretch state transition under prescribed all-point affine deformation only; no time integration, no free crack, no material calibration.",
        },
        "visual_review": {
            "a_image_review": "Existing 1x rest/current overlay is consistent with nanometre-scale motion on a 5 mm coupon; displacement and force traces are labeled; image states no crack, contact, or energy-closure claim.",
            "b_image_review": "Threshold histogram separates the 632 above-threshold broken canonical pairs from 2,430 retained pairs; reciprocal directed stretch points lie on y=x; title/footer state the all-point prescribed affine map and no-free-fracture scope.",
            "b_image_source": "only the existing affine_native_bond_states.csv and affine_canonical_pair_states.csv; no solver run",
        },
        "assertions_passed": True,
    }
    out_path = REVIEW_OUTPUT / "native_resume_audit.json"
    out_path.write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": audit["status"],
        "audit_path": str(out_path),
        "source_pin_count": len(current_pin_audit),
        "states": state_summary,
        "epsilon_c": observed_epsilon_c,
        "raw_csv_hashes": {
            "directed": actual_artifacts[directed_path.name]["sha256"],
            "canonical": actual_artifacts[pair_path.name]["sha256"],
        },
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
