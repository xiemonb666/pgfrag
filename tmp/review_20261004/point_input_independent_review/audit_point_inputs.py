#!/usr/bin/env python3
"""Read-only arithmetic and contract audit for the frozen PD point fixtures."""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
import pathlib
import sys
import tomllib
from collections import Counter


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_pointset(experiment: pathlib.Path, stem: str) -> dict:
    input_dir = experiment / "input"
    csv_path = input_dir / f"{stem}.points.csv"
    metadata_path = input_dir / f"{stem}.metadata.toml"
    metadata = tomllib.loads(metadata_path.read_text(encoding="utf-8"))
    with csv_path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        header = reader.fieldnames
        rows = []
        for row in reader:
            rows.append(
                {
                    "point_id": int(row["point_id"]),
                    "x": float(row["x_m"]),
                    "y": float(row["y_m"]),
                    "z": float(row["z_m"]),
                    "volume": float(row["volume_m3"]),
                }
            )
    return {
        "stem": stem,
        "csv_path": csv_path,
        "metadata_path": metadata_path,
        "csv_sha256": sha256(csv_path),
        "metadata_sha256": sha256(metadata_path),
        "header": header,
        "metadata": metadata,
        "rows": rows,
    }


def close(actual: float, expected: float, *, rel: float = 1e-12, abs_: float = 1e-24) -> bool:
    return math.isclose(actual, expected, rel_tol=rel, abs_tol=abs_)


def pairwise_distance_spectrum(rows: list[dict]) -> list[float]:
    points = [(row["x"], row["y"], row["z"]) for row in rows]
    values = []
    for first, second in itertools.combinations(points, 2):
        values.append(math.dist(first, second))
    return sorted(values)


def analyze(pointset: dict, expected_header: list[str], expected_ids: list[int], density: float) -> dict:
    rows = pointset["rows"]
    metadata = pointset["metadata"]
    ids = [row["point_id"] for row in rows]
    volumes = [row["volume"] for row in rows]
    coords = {axis: [row[axis] for row in rows] for axis in ("x", "y", "z")}
    volume_sum = math.fsum(volumes)
    com = [math.fsum(row[axis] * row["volume"] for row in rows) / volume_sum for axis in ("x", "y", "z")]
    centroid = [math.fsum(coords[axis]) / len(rows) for axis in ("x", "y", "z")]
    expected_meta_com = [float(value) for value in metadata["declared_com_m"]]
    bottom_ids = [int(value) for value in metadata["load_sets"]["bottom_point_ids"]]
    top_ids = [int(value) for value in metadata["load_sets"]["top_point_ids"]]
    by_id = {row["point_id"]: row for row in rows}
    cell_counts = Counter()
    offsets = []
    cell_size_m = 0.001
    for row in rows:
        cell = tuple(int(math.floor(row[axis] / cell_size_m)) for axis in ("x", "y", "z"))
        cell_counts[cell] += 1
        center = tuple((index + 0.5) * cell_size_m for index in cell)
        offsets.append(tuple((row[axis] - center[index]) / cell_size_m for index, axis in enumerate(("x", "y", "z"))))

    failures = []
    if pointset["header"] != expected_header:
        failures.append("CSV header does not exactly match the frozen SI schema")
    if len(rows) != 125:
        failures.append(f"expected 125 point rows, found {len(rows)}")
    if ids != expected_ids or len(ids) != len(set(ids)):
        failures.append("point IDs are not unique stable IDs 1..125 in source row order")
    if any(not math.isfinite(value) for row in rows for value in (row["x"], row["y"], row["z"], row["volume"])):
        failures.append("non-finite coordinate or volume found")
    if any(row["volume"] <= 0 for row in rows):
        failures.append("non-positive point volume found")
    if any(not (0 <= row[axis] <= 0.005) for row in rows for axis in ("x", "y", "z")):
        failures.append("a coordinate is outside the declared 5 mm cube")
    if any(not close(value, 1.0e-9, rel=1e-14, abs_=1e-24) for value in volumes):
        failures.append("per-point cell volume differs from the declared 1e-9 m^3")
    if not close(volume_sum, 125.0e-9):
        failures.append("calculated volume sum differs from 125e-9 m^3")
    if not close(volume_sum, float(metadata["declared_volume_sum_m3"])):
        failures.append("calculated volume sum disagrees with metadata")
    if any(not close(actual, expected) for actual, expected in zip(com, expected_meta_com, strict=True)):
        failures.append("calculated volume-weighted COM disagrees with metadata")
    if len(cell_counts) != 125 or set(cell_counts.values()) != {1}:
        failures.append("the fixture does not place exactly one point in every original cell")
    if sorted(bottom_ids) != list(range(1, 26)) or sorted(top_ids) != list(range(101, 126)):
        failures.append("declared axial load sets do not match source point IDs 1..25 and 101..125")
    if set(bottom_ids) & set(top_ids):
        failures.append("top and bottom point-ID sets overlap")
    if any(by_id[point_id]["z"] >= 0.0025 for point_id in bottom_ids):
        failures.append("bottom IDs do not stay in the bottom source cell layer")
    if any(by_id[point_id]["z"] <= 0.0025 for point_id in top_ids):
        failures.append("top IDs do not stay in the top source cell layer")
    if metadata.get("coordinate_units") != "m" or metadata.get("volume_units") != "m^3":
        failures.append("metadata does not declare SI metres and cubic metres")
    if metadata.get("point_count") != 125:
        failures.append("metadata point_count does not equal 125")
    if metadata.get("quadrature", {}).get("voronoi_claim") is not False or metadata.get("quadrature", {}).get("real_shell_claim") is not False:
        failures.append("metadata lacks explicit non-Voronoi/non-shell claims")

    return {
        "id": metadata["input_id"],
        "csv_sha256": pointset["csv_sha256"],
        "metadata_sha256": pointset["metadata_sha256"],
        "row_count": len(rows),
        "ids_are_1_to_125_in_row_order": ids == expected_ids and len(ids) == len(set(ids)),
        "coordinate_units": metadata.get("coordinate_units"),
        "volume_units": metadata.get("volume_units"),
        "all_coordinates_and_volumes_finite": all(math.isfinite(value) for row in rows for value in (row["x"], row["y"], row["z"], row["volume"])),
        "all_volumes_positive": all(value > 0 for value in volumes),
        "per_point_volume_m3_min": min(volumes),
        "per_point_volume_m3_max": max(volumes),
        "computed_volume_sum_m3": volume_sum,
        "metadata_volume_sum_m3": float(metadata["declared_volume_sum_m3"]),
        "computed_mass_kg_at_declared_density": density * volume_sum,
        "density_kg_m3": density,
        "computed_com_m": com,
        "metadata_com_m": expected_meta_com,
        "unweighted_centroid_m": centroid,
        "computed_bbox_m": {axis: [min(values), max(values)] for axis, values in coords.items()},
        "declared_bottom_ids": bottom_ids,
        "declared_top_ids": top_ids,
        "bottom_z_bbox_m": [min(by_id[value]["z"] for value in bottom_ids), max(by_id[value]["z"] for value in bottom_ids)],
        "top_z_bbox_m": [min(by_id[value]["z"] for value in top_ids), max(by_id[value]["z"] for value in top_ids)],
        "cell_count": len(cell_counts),
        "cell_occupancy_counts": dict(sorted(Counter(cell_counts.values()).items())),
        "quadrature_method": metadata.get("quadrature", {}).get("method"),
        "volume_semantics": metadata.get("quadrature", {}).get("volume_semantics"),
        "voronoi_claim": metadata.get("quadrature", {}).get("voronoi_claim"),
        "real_shell_claim": metadata.get("quadrature", {}).get("real_shell_claim"),
        "failures": failures,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment", required=True, type=pathlib.Path)
    parser.add_argument("--prefreeze", required=True, type=pathlib.Path)
    args = parser.parse_args()
    experiment = args.experiment.resolve()

    expected_hashes = {}
    for line in args.prefreeze.read_text(encoding="ascii").splitlines():
        if not line.strip():
            continue
        digest, relative = line.split("  ", 1)
        expected_hashes[relative.replace("\\", "/")] = digest.lower()
    observed_hashes = {}
    for path in sorted(experiment.rglob("*")):
        if path.is_file():
            observed_hashes[path.relative_to(experiment).as_posix()] = sha256(path)
    critical_paths = {
        "make_point_inputs.ps1",
        "common_case.toml",
        "case_zero_load_control.toml",
        "case_tinyload_event_diagnostic.toml",
        "input/coupon_uniform_voxel_5cube_125.points.csv",
        "input/coupon_uniform_voxel_5cube_125.metadata.toml",
        "input/coupon_cell_jitter_5cube_seed17.points.csv",
        "input/coupon_cell_jitter_5cube_seed17.metadata.toml",
    }
    hash_changes = {
        name: {"prefreeze": expected_hashes.get(name), "observed": digest}
        for name, digest in observed_hashes.items()
        if name in critical_paths and expected_hashes.get(name) != digest
    }
    missing_after_prefreeze = sorted(name for name in critical_paths if name not in observed_hashes)
    added_after_prefreeze = sorted(name for name in critical_paths if name not in expected_hashes)

    generator = experiment / "make_point_inputs.ps1"
    generator_hash = sha256(generator)
    uniform = read_pointset(experiment, "coupon_uniform_voxel_5cube_125")
    jitter = read_pointset(experiment, "coupon_cell_jitter_5cube_seed17")
    header = ["point_id", "x_m", "y_m", "z_m", "volume_m3"]
    expected_ids = list(range(1, 126))
    common = tomllib.loads((experiment / "common_case.toml").read_text(encoding="utf-8"))
    density = float(common["material"]["rho_kg_m3"])
    uniform_result = analyze(uniform, header, expected_ids, density)
    jitter_result = analyze(jitter, header, expected_ids, density)

    uniform_by_id = {row["point_id"]: row for row in uniform["rows"]}
    jitter_by_id = {row["point_id"]: row for row in jitter["rows"]}
    cell_size_m = 0.001
    offsets = []
    max_abs_offset = 0.0
    for point_id in expected_ids:
        base = uniform_by_id[point_id]
        moved = jitter_by_id[point_id]
        vector = tuple((moved[axis] - base[axis]) / cell_size_m for axis in ("x", "y", "z"))
        offsets.append(vector)
        max_abs_offset = max(max_abs_offset, *(abs(value) for value in vector))
    component_histograms = {
        axis: {format(value, ".1f"): count for value, count in sorted(Counter(round(vector[index], 1) for vector in offsets).items())}
        for index, axis in enumerate(("x", "y", "z"))
    }
    component_means = [math.fsum(vector[index] for vector in offsets) / len(offsets) for index in range(3)]
    nonzero_count = sum(any(value != 0.0 for value in vector) for vector in offsets)
    formula_match_count = 0
    for point_id, vector in zip(expected_ids, offsets, strict=True):
        source_index = point_id - 1
        i, j, k = source_index % 5, (source_index // 5) % 5, source_index // 25
        expected_offset = (
            (((i * 17 + j * 3 + k * 4) % 5) - 2) * 0.1,
            (((i * 4 + j * 17 + k * 2) % 5) - 2) * 0.1,
            (((i * 3 + j * 2 + k * 17) % 5) - 2) * 0.1,
        )
        if all(close(actual, expected, rel=0.0, abs_=2e-13) for actual, expected in zip(vector, expected_offset, strict=True)):
            formula_match_count += 1

    uniform_spectrum = pairwise_distance_spectrum(uniform["rows"])
    jitter_spectrum = pairwise_distance_spectrum(jitter["rows"])
    spectral_differences = [
        abs(first - second)
        for first, second in zip(uniform_spectrum, jitter_spectrum, strict=True)
        if round(first, 12) != round(second, 12)
    ]
    manifest = json.loads((experiment / "manifest.json").read_text(encoding="utf-8"))
    manifest_prefix = experiment.relative_to(experiment.parents[2]).as_posix() + "/"
    manifest_hashes = {
        entry["path"][len(manifest_prefix):]: entry["sha256"].lower()
        for entry in manifest.get("file_hashes", [])
        if entry["path"].startswith(manifest_prefix)
    }
    manifest_checks = {}
    for relative in (
        "make_point_inputs.ps1",
        "common_case.toml",
        "case_zero_load_control.toml",
        "case_tinyload_event_diagnostic.toml",
        "input/coupon_uniform_voxel_5cube_125.points.csv",
        "input/coupon_uniform_voxel_5cube_125.metadata.toml",
        "input/coupon_cell_jitter_5cube_seed17.points.csv",
        "input/coupon_cell_jitter_5cube_seed17.metadata.toml",
    ):
        manifest_checks[relative] = {
            "manifest_sha256": manifest_hashes.get(relative),
            "actual_sha256": observed_hashes.get(relative),
            "matches": manifest_hashes.get(relative) == observed_hashes.get(relative),
        }
    fixture_manifest = manifest.get("fixture", {})
    manifest_fixture_checks = {
        "primary_input_id_matches": fixture_manifest.get("primary_input_id") == uniform["metadata"].get("input_id"),
        "primary_csv_path_matches": fixture_manifest.get("primary_points") == f"{manifest_prefix}input/{uniform['stem']}.points.csv",
        "primary_metadata_path_matches": fixture_manifest.get("primary_metadata") == f"{manifest_prefix}input/{uniform['stem']}.metadata.toml",
        "jitter_input_id_matches": fixture_manifest.get("jitter_candidate_input_id") == jitter["metadata"].get("input_id"),
        "jitter_csv_path_matches": fixture_manifest.get("jitter_candidate_points") == f"{manifest_prefix}input/{jitter['stem']}.points.csv",
        "point_count_is_125": fixture_manifest.get("point_count") == 125,
        "extent_is_5mm_cube": fixture_manifest.get("declared_cube_extent_m") == [0.005, 0.005, 0.005],
        "volume_sum_matches_125e_minus_9": close(float(fixture_manifest.get("point_volume_sum_m3", math.nan)), 125.0e-9),
        "com_matches_cube_center": all(close(float(actual), 0.0025) for actual in fixture_manifest.get("volume_weighted_com_m", [])),
    }

    cases = {}
    for case_name in ("case_zero_load_control.toml", "case_tinyload_event_diagnostic.toml"):
        cases[case_name] = tomllib.loads((experiment / case_name).read_text(encoding="utf-8"))
    launch_text = (experiment / "launch_zero_load_probe.ps1").read_text(encoding="utf-8")
    shared_conditions = {
        "horizon_m": float(common["geometry"]["horizon_m"]),
        "material": {key: common["material"][key] for key in ("model", "rho_kg_m3", "E_pa", "nu", "Gc_j_m2", "calibration_status")},
        "solver": {key: common["solver"][key] for key in ("safety_factor", "dt_mode")},
        "case_conditions": {
            name: {
                "load": data["load"],
                "steps": data["solver"]["steps"],
            }
            for name, data in cases.items()
        },
        "launcher_hardcodes_uniform_input": "coupon_uniform_voxel_5cube_125.points.csv" in launch_text and "coupon_uniform_voxel_5cube_125.metadata.toml" in launch_text,
        "launcher_hardcodes_zero_load_case": "case_zero_load_control.toml" in launch_text,
        "launcher_has_jitter_candidate_input": "coupon_cell_jitter_5cube_seed17" in launch_text,
        "launcher_sha256_at_audit": sha256(experiment / "launch_zero_load_probe.ps1"),
        "launcher_same_conditions_ab_comparison_exists": False,
        "interpretation": "Both input contracts share the same common horizon/material and identical source-ID axial sets. The checked-in launcher binds only the uniform fixture to zero_load_control; it does not run the jitter candidate, so no executed or configured same-case A/B comparison is established by these files.",
    }
    all_failures = uniform_result["failures"] + jitter_result["failures"]
    if uniform["metadata"].get("source", {}).get("source_sha256", "").lower() != generator_hash:
        all_failures.append("uniform metadata source_sha256 does not match generator")
    if jitter["metadata"].get("source", {}).get("source_sha256", "").lower() != generator_hash:
        all_failures.append("jitter metadata source_sha256 does not match generator")
    if uniform["metadata"].get("point_data_sha256", "").lower() != uniform["csv_sha256"]:
        all_failures.append("uniform metadata point_data_sha256 does not match CSV")
    if jitter["metadata"].get("point_data_sha256", "").lower() != jitter["csv_sha256"]:
        all_failures.append("jitter metadata point_data_sha256 does not match CSV")
    if hash_changes or missing_after_prefreeze or added_after_prefreeze:
        all_failures.append("experiment source/input tree changed since prefreeze")
    if any(not record["matches"] for record in manifest_checks.values()):
        all_failures.append("manifest does not pin all reviewed source/input hashes")
    if not all(manifest_fixture_checks.values()):
        all_failures.append("manifest fixture declarations do not agree with reviewed point inputs")
    if max_abs_offset > 0.2 + 1e-12:
        all_failures.append("jitter exceeded 0.2 cell widths")
    if formula_match_count != 125:
        all_failures.append("jitter coordinates do not reproduce the declared deterministic source formula")

    output = {
        "audit": "PG-Frag prepared point inputs independent arithmetic review",
        "python_version": sys.version,
        "argv": sys.argv,
        "experiment_directory": str(experiment),
        "prefreeze_sha_file": str(args.prefreeze.resolve()),
        "prefreeze_file_count": len(expected_hashes),
        "critical_prefreeze_path_count": len(critical_paths),
        "relevant_file_sha256": {name: observed_hashes.get(name) for name in sorted(critical_paths)},
        "prefreeze_hash_changes": hash_changes,
        "files_missing_after_prefreeze": missing_after_prefreeze,
        "files_added_after_prefreeze": added_after_prefreeze,
        "generator_sha256": generator_hash,
        "manifest_sha256_at_audit": sha256(experiment / "manifest.json"),
        "manifest_status_at_audit": manifest.get("status"),
        "manifest_fixture_at_audit": manifest.get("fixture", {}),
        "manifest_fixture_checks": manifest_fixture_checks,
        "manifest_hash_checks": manifest_checks,
        "fixture_checks": {"uniform": uniform_result, "jitter": jitter_result},
        "jitter_analysis": {
            "source": "make_point_inputs.ps1 uses fixed integer pattern seed 17; offsets are deterministic modular arithmetic, not independent random samples.",
            "formula_matches_by_point_id": formula_match_count,
            "one_point_per_original_cell": len(jitter_result["cell_occupancy_counts"]) == 1 and jitter_result["cell_count"] == 125 and jitter_result["cell_occupancy_counts"].get(1) == 125,
            "max_abs_offset_cell_widths": max_abs_offset,
            "component_offset_histograms_cell_widths": component_histograms,
            "mean_offset_cell_widths": component_means,
            "points_with_any_nonzero_offset": nonzero_count,
            "points_with_all_zero_offset": 125 - nonzero_count,
            "does_not_establish_elimination_of_directional_grid_bias": True,
        },
        "between_fixture_comparison": {
            "same_point_ids_and_row_order": [row["point_id"] for row in uniform["rows"]] == [row["point_id"] for row in jitter["rows"]],
            "point_ids_with_coordinate_change": sum(any(uniform_by_id[point_id][axis] != jitter_by_id[point_id][axis] for axis in ("x", "y", "z")) for point_id in expected_ids),
            "pairwise_distance_count_each": len(uniform_spectrum),
            "pairwise_distance_spectrum_mismatches_at_1e-12_m_resolution": len(spectral_differences),
            "max_sorted_pairwise_spectrum_difference_m": max(spectral_differences, default=0.0),
            "rigid_rotation_or_point_permutation_equivalent": len(spectral_differences) == 0,
            "interpretation": "A different sorted pairwise-distance multiset rules out equivalence under any rigid rotation/translation plus permutation of points. This proves geometric difference, not improved physical representativeness.",
        },
        "shared_conditions": shared_conditions,
        "overall_static_arithmetic_status": "PASS" if not all_failures else "FAIL",
        "failures": all_failures,
        "limits": [
            "No Julia, solver, native run, pytest, or launch preflight was executed.",
            "This verifies declared fixture arithmetic and metadata only; it does not validate runtime loading or solver behavior.",
            "The weights are voxel cell quadrature assignments, not exact Voronoi volumes and not a parent specimen volume.",
            "The fixed deterministic jitter pattern does not demonstrate elimination of grid or orientation bias.",
            "Material constants are declared UNVALIDATED_PROTOTYPE values and do not establish artifact eligibility.",
        ],
    }
    print(json.dumps(output, indent=2, sort_keys=True))
    return 0 if output["overall_static_arithmetic_status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
