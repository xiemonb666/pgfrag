"""Read-only independent audit for the retained single-ROI pilot.

Reads only the frozen retained stage and its supervisor summary. Expected
geometry comes from analytic_roi_reference.py; this script never imports
CadQuery/OpenCascade or invokes the producer worker. Run only after root review
and explicit authorization.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time

from analytic_roi_reference import (
    CELL_MM, INNER_RADIUS_MM, OUTER_RADIUS_MM, ROI_MM,
    expected_cells, expected_roi,
)

AUDIT_DIR = Path(__file__).resolve().parent
WORKSPACE = Path(__file__).resolve().parents[4]
STAGE = WORKSPACE / "pgfrag/experiments/forge_volume_ingress_20261005/tmp/.forge_ingress_stage__2srs6ad"
SUMMARY = WORKSPACE / "pgfrag/experiments/forge_volume_ingress_20261005/tmp/candidate_run01_supervisor/pilot_results_summary.json"
FREEZE = AUDIT_DIR / "PILOT_AUDIT_INPUT_FREEZE.json"
OUTPUT_JSON = AUDIT_DIR / "PILOT_ACTUAL_AUDIT.json"
CELL_CSV = AUDIT_DIR / "PILOT_CELL_AUDIT.csv"
POINT_CSV = AUDIT_DIR / "PILOT_POINT_AUDIT.csv"
MAP_PNG = AUDIT_DIR / "PILOT_XY_CELL_MAP.png"
DELTA_PNG = AUDIT_DIR / "PILOT_CELL_VOLUME_MASS_DELTA.png"

RHO_KG_M3 = 2100.0
VOLUME_RTOL = 1e-8
VOLUME_ATOL_MM3 = 1e-9
POSITION_ATOL_MM = 1e-9
CELL_CENTROID_TOL_MM = 1e-6


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as stream:
        value = json.load(stream)
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path.name}")
    return value


def close(a: float, b: float, *, rtol: float, atol: float) -> bool:
    return math.isfinite(a) and math.isfinite(b) and abs(a - b) <= max(
        atol, rtol * abs(b)
    )


def check(name: str, passed: bool, details: dict | None = None) -> dict:
    return {"name": name, "passed": bool(passed), "details": details or {}}


def verify_inputs(freeze: dict) -> dict:
    expected_stage = freeze["retained_stage_files"]
    actual_names = sorted(p.name for p in STAGE.iterdir() if p.is_file())
    frozen_names = sorted(expected_stage)
    file_rows = []
    matches = actual_names == frozen_names
    for name in frozen_names:
        path = STAGE / name
        exists = path.is_file()
        size = path.stat().st_size if exists else None
        digest = sha256(path) if exists else None
        expected = expected_stage[name]
        row = {
            "path": name,
            "exists": exists,
            "bytes": size,
            "sha256": digest,
            "matches_freeze": exists
            and size == expected["bytes"]
            and digest == expected["sha256"],
        }
        matches &= row["matches_freeze"]
        file_rows.append(row)
    summary_pin = freeze["supervisor_summary"]
    summary_matches = (
        SUMMARY.stat().st_size == summary_pin["bytes"]
        and sha256(SUMMARY) == summary_pin["sha256"]
    )
    total = sum(row["bytes"] or 0 for row in file_rows)
    return {
        "stage_path": str(STAGE),
        "summary_path": str(SUMMARY),
        "stage_file_names_match_freeze": actual_names == frozen_names,
        "stage_files": file_rows,
        "stage_total_bytes": total,
        "stage_within_20_mib": total <= 20 * 1024 * 1024,
        "summary_matches_freeze": summary_matches,
        "all_frozen_inputs_match": bool(matches and summary_matches),
    }


def inside(point: tuple[float, float, float], box: tuple[float, ...]) -> bool:
    x0, x1, y0, y1, z0, z1 = box
    return (
        x0 - POSITION_ATOL_MM <= point[0] <= x1 + POSITION_ATOL_MM
        and y0 - POSITION_ATOL_MM <= point[1] <= y1 + POSITION_ATOL_MM
        and z0 - POSITION_ATOL_MM <= point[2] <= z1 + POSITION_ATOL_MM
    )


def write_csv(path: Path, columns: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def make_figures(analytic_rows: list, points: list[dict], cell_rows: list[dict]) -> dict:
    os.environ.setdefault("MPLCONFIGDIR", str(AUDIT_DIR / ".mplconfig"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import cm, colors
    from matplotlib.patches import Rectangle

    xmin, xmax, ymin, ymax, _, _ = ROI_MM
    area_by_ij = {
        (r.grid_ijk[0], r.grid_ijk[1]): r.expected_area_mm2
        for r in analytic_rows
    }
    max_fraction = max(
        (a / CELL_MM**2 for a in area_by_ij.values()), default=1.0
    )
    norm = colors.Normalize(0.0, max(max_fraction, 1e-15))
    cmap = cm.get_cmap("viridis")
    points_by_k = {k: [] for k in range(4)}
    for point in points:
        source = point["source"]
        k = int(source["grid_ijk"][2])
        xyz = [float(v) * 1000.0 for v in point["position_m"]]
        points_by_k[k].append((xyz[0], xyz[1]))

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True, sharey=True)
    angles = [2 * math.pi * n / 720 for n in range(721)]
    for k, ax in enumerate(axes.flat):
        for j in range(4):
            for i in range(6):
                x0, y0 = xmin + i * CELL_MM, ymin + j * CELL_MM
                fraction = area_by_ij[(i, j)] / CELL_MM**2
                ax.add_patch(
                    Rectangle(
                        (x0, y0), CELL_MM, CELL_MM,
                        facecolor=cmap(norm(fraction)),
                        edgecolor="white", linewidth=0.7,
                    )
                )
        for radius, style, color in (
            (OUTER_RADIUS_MM, "-", "black"),
            (INNER_RADIUS_MM, "--", "dimgray"),
        ):
            ax.plot(
                [radius * math.cos(t) for t in angles],
                [radius * math.sin(t) for t in angles],
                linestyle=style, color=color, linewidth=1.0,
            )
        slab_points = points_by_k[k]
        if slab_points:
            ax.scatter(
                [p[0] for p in slab_points], [p[1] for p in slab_points],
                s=26, facecolor="white", edgecolor="black", linewidth=0.65,
                zorder=5,
            )
        zmid = ROI_MM[4] + (k + 0.5) * CELL_MM
        ax.set_title(f"z slab {k}, centroid z={zmid:.2f} mm")
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
    mapper = cm.ScalarMappable(norm=norm, cmap=cmap)
    mapper.set_array([])
    fig.colorbar(mapper, ax=list(axes.flat), shrink=0.82,
                 label="analytic annular material area fraction per XY cell")
    fig.suptitle("Retained local geometry; final report failed, result not promoted")
    fig.tight_layout()
    fig.savefig(MAP_PNG, dpi=170)
    plt.close(fig)

    indices = [r["grid_cell_index"] for r in cell_rows]
    analytic_vol = [r["analytic_volume_mm3"] for r in cell_rows]
    raw_vol = [r["raw_common_volume_mm3"] for r in cell_rows]
    point_vol = [r["delivered_point_volume_mm3"] for r in cell_rows]
    analytic_mass = [r["analytic_mass_kg"] for r in cell_rows]
    point_mass = [r["delivered_point_mass_kg"] for r in cell_rows]
    fig, axes = plt.subplots(3, 1, figsize=(13, 9), sharex=True)
    axes[0].plot(indices, analytic_vol, label="analytic cell volume")
    axes[0].scatter(indices, raw_vol, s=12, label="raw B-rep common")
    axes[0].scatter(indices, point_vol, s=12, marker="x",
                    label="sum of delivered point volumes")
    axes[0].set_ylabel("volume (mm³)")
    axes[0].legend(fontsize=8)
    axes[1].plot(indices, analytic_mass, label="rho × analytic volume")
    axes[1].scatter(indices, point_mass, s=12, label="sum of point mass")
    axes[1].set_ylabel("mass (kg)")
    axes[1].legend(fontsize=8)
    axes[2].axhline(0.0, color="black", linewidth=0.8)
    axes[2].plot(
        indices,
        [point_vol[n] - analytic_vol[n] for n in range(len(indices))],
        label="delivered − analytic",
    )
    axes[2].plot(
        indices,
        [raw_vol[n] - analytic_vol[n] for n in range(len(indices))],
        label="raw B-rep − analytic",
    )
    axes[2].set_xlabel("stable cell index (k, then j, then i)")
    axes[2].set_ylabel("signed volume delta (mm³)")
    axes[2].legend(fontsize=8)
    for ax in axes:
        ax.grid(True, alpha=0.25)
    fig.suptitle("Per-cell checks; local diagnostic only, not promoted")
    fig.tight_layout()
    fig.savefig(DELTA_PNG, dpi=170)
    plt.close(fig)
    return {
        "xy_map_png": MAP_PNG.name,
        "cell_volume_mass_delta_png": DELTA_PNG.name,
    }


def run_audit() -> dict:
    started = time.perf_counter()
    freeze = read_json(FREEZE)
    frozen = verify_inputs(freeze)
    if not frozen["all_frozen_inputs_match"]:
        raise RuntimeError("stage or summary does not match the audit input freeze")

    summary = read_json(SUMMARY)
    pointset = read_json(STAGE / "material_points.json")
    ledger = read_json(STAGE / "cell_ledger.json")
    partial = read_json(STAGE / "cell_ledger.partial.json")
    parent = read_json(STAGE / "parent_execution_receipt.json")
    worker = read_json(STAGE / "worker_failure.json")
    analytic_rows = expected_cells()
    analytic_roi = expected_roi()
    checks: list[dict] = []

    def add(name: str, passed: bool, details: dict | None = None) -> None:
        checks.append(check(name, passed, details))

    expected = {r.grid_cell_index: r for r in analytic_rows}
    cells = ledger.get("cells", [])
    points = pointset.get("points", [])
    cell_by_id = {int(c["grid_cell_index"]): c for c in cells}
    point_by_id = {int(p["point_id"]): p for p in points}
    point_ids = [int(p["point_id"]) for p in points]
    duplicate_point_rows = sorted(
        pid for pid in set(point_ids) if point_ids.count(pid) > 1
    )
    summary_stage_rows = summary.get("retained_stage", {}).get("files", [])
    summary_stage_hashes = {
        row.get("path"): row.get("sha256") for row in summary_stage_rows
    }
    frozen_stage_hashes = {
        row["path"]: row["sha256"] for row in frozen["stage_files"]
    }
    summary_stage_hash_match = (
        summary_stage_hashes == frozen_stage_hashes
        and summary.get("retained_stage", {}).get("total_bytes_including_tmp")
        == frozen["stage_total_bytes"]
    )
    add(
        "retained_stage_inventory_matches_summary",
        summary_stage_hash_match,
        {
            "summary_inventory_matches_frozen_stage": summary_stage_hash_match,
            "summary_stage_files": len(summary_stage_hashes),
            "actual_stage_files": len(frozen_stage_hashes),
        },
    )
    add(
        "96_cells",
        len(cells) == 96
        and len(cell_by_id) == 96
        and ledger.get("candidate_cell_count") == 96
        and ledger.get("processed_cell_count") == 96,
        {
            "ledger_cells": len(cells),
            "unique_cell_ids": len(cell_by_id),
            "points": len(points),
            "duplicate_point_ids": duplicate_point_rows,
        },
    )

    bad_cell_maps: list[int] = []
    for index, cell in cell_by_id.items():
        ref = expected.get(index)
        if ref is None or tuple(cell.get("grid_ijk", ())) != ref.grid_ijk:
            bad_cell_maps.append(index)
            continue
        aabb = tuple(float(v) for v in cell.get("cell_aabb_mm", []))
        if len(aabb) != 6 or any(
            abs(aabb[i] - ref.cell_aabb_mm[i]) > 1e-12 for i in range(6)
        ):
            bad_cell_maps.append(index)
    add("cell_index_ijk_aabb_map", not bad_cell_maps,
        {"bad_indices": sorted(set(bad_cell_maps))})

    raw_cell_sum = 0.0
    component_sum = 0.0
    point_sum_m3 = 0.0
    mass_sum = 0.0
    class_counts: dict[str, int] = {}
    topology_counts: dict[str, int] = {}
    candidate_component_by_id: dict[int, tuple[int, dict]] = {}
    link_errors: list[str] = []
    point_conversion_errors: list[int] = []
    mass_errors: list[int] = []
    position_errors: list[int] = []
    bound_errors: list[int] = []
    analytic_membership_errors: list[int] = []
    centroid_error_cells: list[int] = []
    max_centroid_delta_by_axis = [0.0, 0.0, 0.0]
    max_centroid_abs_delta_mm = 0.0
    classifier_errors: list[int] = []
    duplicate_point_ids: list[int] = []
    positive_components = 0
    cell_rows: list[dict] = []
    point_rows: list[dict] = []
    raw_cell_volume_errors: list[int] = []
    zero_cell_errors: list[int] = []
    component_volume_errors: list[int] = []
    positive_status_errors: list[int] = []
    source_mapping_errors: list[int] = []

    for index in sorted(expected):
        ref = expected[index]
        cell = cell_by_id.get(index)
        if cell is None:
            continue
        raw_mm3 = float(cell.get("raw_common_volume_mm3") or 0.0)
        raw_cell_sum += raw_mm3
        components = cell.get("component_records", [])
        comp_sum = math.fsum(float(c.get("raw_volume_mm3") or 0.0)
                             for c in components)
        component_sum += comp_sum
        point_mm3_sum = 0.0
        mass_kg_sum = 0.0
        mx = my = mz = 0.0
        for comp in components:
            positive_components += 1
            state = str(comp.get("centroid_source_classification", "UNKNOWN"))
            class_counts[state] = class_counts.get(state, 0) + 1
            topo = (comp.get("component_topology") or {}).get("status", "UNKNOWN")
            topology_counts[str(topo)] = topology_counts.get(str(topo), 0) + 1
            pid_value = comp.get("candidate_point_id")
            if pid_value is None:
                link_errors.append(f"cell {index}: component has no point id")
                continue
            pid = int(pid_value)
            if pid in candidate_component_by_id:
                duplicate_point_ids.append(pid)
            candidate_component_by_id[pid] = (index, comp)
            point = point_by_id.get(pid)
            if point is None:
                link_errors.append(f"component point id {pid} absent from pointset")
                continue
            src = point.get("source", {})
            if (
                int(src.get("grid_cell_index", -1)) != index
                or tuple(src.get("grid_ijk", ())) != ref.grid_ijk
                or int(src.get("raw_component_index", -1))
                != int(comp.get("raw_component_index", -2))
                or int(src.get("intersection_component_index", -1))
                != int(comp.get("intersection_component_index", -2))
            ):
                link_errors.append(f"point {pid}: cell/component map mismatch")
            volume_m3 = float(point["cell_volume_m3"])
            volume_mm3 = volume_m3 * 1e9
            component_mm3 = float(comp.get("raw_volume_mm3") or 0.0)
            if not close(
                volume_mm3, component_mm3,
                rtol=1e-12, atol=1e-10,
            ):
                point_conversion_errors.append(pid)
            mass = float(point["mass_kg"])
            expected_mass = RHO_KG_M3 * volume_m3
            if not close(mass, expected_mass, rtol=1e-12, atol=1e-15):
                mass_errors.append(pid)
            position_m = tuple(float(v) for v in point["position_m"])
            position_mm = tuple(v * 1000.0 for v in position_m)
            centroid = tuple(float(v) for v in comp.get("centroid_mm", []))
            position_ok = len(centroid) == 3 and all(
                abs(position_mm[n] - centroid[n]) <= POSITION_ATOL_MM
                for n in range(3)
            )
            if not position_ok:
                position_errors.append(pid)
            if not inside(position_mm, ref.cell_aabb_mm):
                bound_errors.append(pid)
            radial_distance = math.hypot(position_mm[0], position_mm[1])
            z_in_source = 0.0 <= position_mm[2] <= 10.0
            z_in_cell = (
                ref.cell_aabb_mm[4] <= position_mm[2] <= ref.cell_aabb_mm[5]
            )
            analytic_member = (
                INNER_RADIUS_MM < radial_distance < OUTER_RADIUS_MM
                and z_in_source
                and z_in_cell
            )
            if not analytic_member:
                analytic_membership_errors.append(pid)
            if (
                state != "IN"
                or comp.get("candidate_point_status") != "CANDIDATE_CENTROID_IN"
            ):
                classifier_errors.append(pid)
            point_mm3_sum += volume_mm3
            mass_kg_sum += mass
            mx += position_mm[0] * volume_mm3
            my += position_mm[1] * volume_mm3
            mz += position_mm[2] * volume_mm3
            point_sum_m3 += volume_m3
            mass_sum += mass
            point_rows.append(
                {
                    "point_id": pid,
                    "grid_cell_index": index,
                    "grid_i": ref.grid_ijk[0],
                    "grid_j": ref.grid_ijk[1],
                    "grid_k": ref.grid_ijk[2],
                    "raw_component_index": src.get("raw_component_index"),
                    "intersection_component_index": src.get(
                        "intersection_component_index"
                    ),
                    "x_m": position_m[0], "y_m": position_m[1],
                    "z_m": position_m[2],
                    "x_mm_from_m": position_mm[0],
                    "y_mm_from_m": position_mm[1],
                    "z_mm_from_m": position_mm[2],
                    "ledger_centroid_x_mm": centroid[0] if len(centroid) == 3 else None,
                    "ledger_centroid_y_mm": centroid[1] if len(centroid) == 3 else None,
                    "ledger_centroid_z_mm": centroid[2] if len(centroid) == 3 else None,
                    "position_matches_raw_ledger_centroid": position_ok,
                    "analytic_annulus_and_slab_membership": analytic_member,
                    "point_volume_m3": volume_m3,
                    "point_volume_converted_mm3": volume_mm3,
                    "component_raw_volume_mm3": component_mm3,
                    "volume_delta_mm3": volume_mm3 - component_mm3,
                    "mass_kg": mass,
                    "rho_times_volume_kg": expected_mass,
                    "mass_delta_kg": mass - expected_mass,
                    "classifier": state,
                    "candidate_point_status": comp.get("candidate_point_status"),
                    "source_semantic_region_status": src.get(
                        "source_semantic_region_status"
                    ),
                    "surface_parent_mapping_status": src.get(
                        "surface_parent_mapping_status"
                    ),
                }
            )

        tol = max(VOLUME_ATOL_MM3, VOLUME_RTOL * abs(ref.expected_volume_mm3))
        if ref.expected_volume_mm3 > 0:
            if abs(raw_mm3 - ref.expected_volume_mm3) > tol:
                raw_cell_volume_errors.append(index)
            if cell.get("cell_status") != "POSITIVE_COMPONENTS_CENTROIDS_IN":
                positive_status_errors.append(index)
        elif (
            abs(raw_mm3) > VOLUME_ATOL_MM3
            or components
            or cell.get("candidate_point_ids")
            or cell.get("cell_status") != "ZERO_INTERSECTION"
        ):
            zero_cell_errors.append(index)
        if not close(comp_sum, raw_mm3, rtol=VOLUME_RTOL, atol=VOLUME_ATOL_MM3):
            component_volume_errors.append(index)
        actual_centroid = (
            (mx / point_mm3_sum, my / point_mm3_sum, mz / point_mm3_sum)
            if point_mm3_sum > 0 else (None, None, None)
        )
        expected_centroid = ref.expected_centroid_mm or (None, None, None)
        centroid_delta = (None, None, None)
        centroid_delta_max = None
        centroid_matches_analytic = ref.expected_centroid_mm is None
        if ref.expected_centroid_mm is not None:
            if (
                point_mm3_sum > 0.0
                and all(
                    value is not None and math.isfinite(float(value))
                    for value in actual_centroid
                )
                and all(math.isfinite(value) for value in ref.expected_centroid_mm)
            ):
                centroid_delta = tuple(
                    float(actual_centroid[n]) - ref.expected_centroid_mm[n]
                    for n in range(3)
                )
                centroid_delta_max = max(abs(value) for value in centroid_delta)
                max_centroid_abs_delta_mm = max(
                    max_centroid_abs_delta_mm, centroid_delta_max
                )
                for axis in range(3):
                    max_centroid_delta_by_axis[axis] = max(
                        max_centroid_delta_by_axis[axis],
                        abs(centroid_delta[axis]),
                    )
                centroid_matches_analytic = (
                    centroid_delta_max <= CELL_CENTROID_TOL_MM
                )
            else:
                centroid_matches_analytic = False
            if not centroid_matches_analytic:
                centroid_error_cells.append(index)
        cell_rows.append(
            {
                "grid_cell_index": index,
                "grid_i": ref.grid_ijk[0], "grid_j": ref.grid_ijk[1],
                "grid_k": ref.grid_ijk[2],
                "cell_aabb_mm": json.dumps(ref.cell_aabb_mm),
                "analytic_area_mm2": ref.expected_area_mm2,
                "analytic_volume_mm3": ref.expected_volume_mm3,
                "raw_common_volume_mm3": raw_mm3,
                "raw_component_volume_mm3": comp_sum,
                "delivered_point_volume_mm3": point_mm3_sum,
                "delivered_point_volume_m3": point_mm3_sum * 1e-9,
                "analytic_mass_kg": ref.expected_volume_mm3 * 1e-9 * RHO_KG_M3,
                "delivered_point_mass_kg": mass_kg_sum,
                "raw_minus_analytic_mm3": raw_mm3 - ref.expected_volume_mm3,
                "points_minus_analytic_mm3": point_mm3_sum - ref.expected_volume_mm3,
                "tolerance_mm3": tol,
                "analytic_centroid_x_mm": expected_centroid[0],
                "analytic_centroid_y_mm": expected_centroid[1],
                "analytic_centroid_z_mm": expected_centroid[2],
                "delivered_centroid_x_mm": actual_centroid[0],
                "delivered_centroid_y_mm": actual_centroid[1],
                "delivered_centroid_z_mm": actual_centroid[2],
                "centroid_delta_x_mm": centroid_delta[0],
                "centroid_delta_y_mm": centroid_delta[1],
                "centroid_delta_z_mm": centroid_delta[2],
                "centroid_delta_max_abs_mm": centroid_delta_max,
                "analytic_centroid_tolerance_mm": CELL_CENTROID_TOL_MM,
                "centroid_matches_analytic": centroid_matches_analytic,
                "component_count": len(components),
                "classifier_in_count": sum(
                    c.get("centroid_source_classification") == "IN"
                    for c in components
                ),
                "cell_status": cell.get("cell_status"),
            }
        )

    unlinked_points = sorted(set(point_by_id) - set(candidate_component_by_id))
    component_ids_missing_rows = sorted(set(candidate_component_by_id) - set(point_by_id))
    add(
        "cell_common_matches_independent_analytic",
        not raw_cell_volume_errors
        and not zero_cell_errors
        and not positive_status_errors,
        {
            "positive_cell_mismatches": raw_cell_volume_errors,
            "positive_cell_status_errors": positive_status_errors,
            "zero_cell_mismatches": zero_cell_errors,
        },
    )
    add(
        "cell_components_reconcile",
        not component_volume_errors,
        {"mismatch_indices": component_volume_errors},
    )
    add(
        "component_point_mapping",
        not link_errors
        and not unlinked_points
        and not component_ids_missing_rows
        and not duplicate_point_ids,
        {
            "errors": link_errors,
            "unlinked_point_ids": unlinked_points,
            "component_ids_missing_point_rows": component_ids_missing_rows,
            "duplicate_ids": duplicate_point_ids,
        },
    )
    add(
        "SI_volume_and_rhoV",
        not point_conversion_errors and not mass_errors,
        {
            "volume_conversion_error_ids": point_conversion_errors,
            "mass_error_ids": mass_errors,
            "rho_kg_m3": RHO_KG_M3,
        },
    )
    add(
        "strict_classifier_and_centroid_mapping",
        class_counts == {"IN": positive_components}
        and topology_counts == {"PASS": positive_components}
        and not classifier_errors
        and not position_errors
        and not bound_errors,
        {
            "classifier_counts": class_counts,
            "topology_counts": topology_counts,
            "classifier_error_ids": classifier_errors,
            "position_error_ids": position_errors,
            "outside_cell_ids": bound_errors,
        },
    )
    add(
        "volume_weighted_cell_centroid_matches_independent_analytic",
        not centroid_error_cells,
        {
            "tolerance_mm": CELL_CENTROID_TOL_MM,
            "axes_checked": ["x", "y", "z"],
            "max_absolute_delta_by_axis_mm": {
                "x": max_centroid_delta_by_axis[0],
                "y": max_centroid_delta_by_axis[1],
                "z": max_centroid_delta_by_axis[2],
            },
            "max_absolute_delta_mm": max_centroid_abs_delta_mm,
            "mismatched_positive_cell_indices": centroid_error_cells,
        },
    )
    add(
        "analytic_annulus_and_slab_membership",
        not analytic_membership_errors,
        {
            "strict_radial_interval_mm": [
                INNER_RADIUS_MM,
                OUTER_RADIUS_MM,
            ],
            "source_z_interval_mm": [0.0, 10.0],
            "point_ids_outside_analytic_material_or_cell_slab": (
                analytic_membership_errors
            ),
        },
    )
    source_mapping_errors = [
        int(point["point_id"])
        for point in points
        if point.get("source", {}).get("source_semantic_region_status")
        != "UNKNOWN_NOT_PROJECTED"
        or point.get("source", {}).get("surface_parent_mapping_status")
        != "UNKNOWN_NOT_PROJECTED"
    ]
    provenance_rho = pointset.get("provenance", {}).get("rho_kg_m3")
    add(
        "declared_density_and_unprojected_surface_mapping",
        (provenance_rho is None or float(provenance_rho) == RHO_KG_M3)
        and not source_mapping_errors,
        {
            "fixture_declared_rho_kg_m3": RHO_KG_M3,
            "pointset_provenance_rho_kg_m3": provenance_rho,
            "points_with_unexpected_surface_mapping": source_mapping_errors,
        },
    )
    add(
        "all_cells_and_points_map",
        len(cells) == 96
        and len(points) == pointset.get("point_count")
        and len(point_by_id) == len(points)
        and not duplicate_point_rows
        and not bad_cell_maps,
        {
            "cell_count": len(cells),
            "point_count": len(points),
            "pointset_declared_count": pointset.get("point_count"),
            "duplicate_point_ids": duplicate_point_rows,
            "bad_cell_map_indices": bad_cell_maps,
        },
    )

    direct = float(analytic_roi["direct_volume_mm3"])
    analytic_sum = math.fsum(r.expected_volume_mm3 for r in analytic_rows)
    ref_record = ledger.get("reference_roi_common", {})
    raw_reference = ref_record.get(
        "raw_common_volume_mm3", ref_record.get("raw_volume_mm3")
    )
    raw_reference = float(raw_reference) if raw_reference is not None else None
    raw_reference_minus_sum = (
        raw_reference - raw_cell_sum if raw_reference is not None else None
    )
    point_sum_mm3 = point_sum_m3 * 1e9
    tol = max(VOLUME_ATOL_MM3, VOLUME_RTOL * abs(direct))
    whole = {
        "independent_analytic_direct_roi_mm3": direct,
        "independent_analytic_sum_96_cells_mm3": analytic_sum,
        "analytic_direct_minus_grid_mm3": analytic_roi[
            "direct_minus_grid_volume_mm3"
        ],
        "raw_brep_reference_roi_mm3": raw_reference,
        "raw_sum_96_cell_commons_mm3": raw_cell_sum,
        "raw_reference_minus_cell_sum_mm3": raw_reference_minus_sum,
        "raw_connected_component_sum_mm3": component_sum,
        "delivered_point_volume_sum_mm3": point_sum_mm3,
        "volume_tolerance_mm3": tol,
        "reference_matches_analytic": raw_reference is not None
        and close(raw_reference, direct, rtol=VOLUME_RTOL, atol=VOLUME_ATOL_MM3),
        "cell_sum_matches_analytic": close(
            raw_cell_sum, direct, rtol=VOLUME_RTOL, atol=VOLUME_ATOL_MM3
        ),
        "reference_matches_cell_sum": raw_reference_minus_sum is not None
        and abs(raw_reference_minus_sum) <= tol,
        "point_sum_matches_analytic": close(
            point_sum_mm3, direct, rtol=VOLUME_RTOL, atol=VOLUME_ATOL_MM3
        ),
        "component_sum_matches_cell_sum": close(
            component_sum, raw_cell_sum, rtol=VOLUME_RTOL, atol=VOLUME_ATOL_MM3
        ),
        "emitted_mass_sum_kg": mass_sum,
        "rho_times_point_volume_sum_kg": RHO_KG_M3 * point_sum_m3,
    }
    whole["all_volume_reconciliations_pass"] = all(
        [
            whole["reference_matches_analytic"],
            whole["cell_sum_matches_analytic"],
            whole["reference_matches_cell_sum"],
            whole["point_sum_matches_analytic"],
            whole["component_sum_matches_cell_sum"],
        ]
    )
    add(
        "independent_reference_vs_96_cell_reconciliation",
        whole["all_volume_reconciliations_pass"],
        whole,
    )

    finalization = summary.get("finalization_failure", {})
    promoted = parent.get("candidate_run_promoted")
    no_promotion = (
        summary.get("pilot_final_status") == "BLOCKED_REPORT_FINALIZATION; local geometry ledger retained"
        and summary.get("run", {}).get("outer_exit_code") == 2
        and finalization.get("report_attempt_was_promoted") is False
        and finalization.get("ingress_report_present") is False
        and parent.get("status") == "BLOCKED_STAGE_RETAINED"
        and promoted is False
        and not (STAGE / "ingress_report.json").exists()
    )
    add(
        "finalization_failure_not_promoted",
        no_promotion,
        {
            "final_status": summary.get("pilot_final_status"),
            "outer_exit_code": summary.get("run", {}).get("outer_exit_code"),
            "finalization_reason": finalization.get("reason"),
            "parent_receipt_status": parent.get("status"),
            "candidate_run_promoted": promoted,
            "ingress_report_in_retained_stage": (STAGE / "ingress_report.json").exists(),
            "raw_candidate_file_retained": (STAGE / "material_points.json").is_file(),
            "cell_ledger_status": ledger.get("status"),
            "partial_ledger_status": partial.get("status"),
            "worker_failure_class": worker.get("failure_class"),
        },
    )

    write_csv(
        CELL_CSV,
        [
            "grid_cell_index", "grid_i", "grid_j", "grid_k", "cell_aabb_mm",
            "analytic_area_mm2", "analytic_volume_mm3", "raw_common_volume_mm3",
            "raw_component_volume_mm3", "delivered_point_volume_mm3",
            "delivered_point_volume_m3", "analytic_mass_kg",
            "delivered_point_mass_kg", "raw_minus_analytic_mm3",
            "points_minus_analytic_mm3", "tolerance_mm3",
            "analytic_centroid_x_mm", "analytic_centroid_y_mm",
            "analytic_centroid_z_mm", "delivered_centroid_x_mm",
            "delivered_centroid_y_mm", "delivered_centroid_z_mm",
            "centroid_delta_x_mm", "centroid_delta_y_mm",
            "centroid_delta_z_mm", "centroid_delta_max_abs_mm",
            "analytic_centroid_tolerance_mm", "centroid_matches_analytic",
            "component_count", "classifier_in_count", "cell_status",
        ],
        cell_rows,
    )
    write_csv(
        POINT_CSV,
        [
            "point_id", "grid_cell_index", "grid_i", "grid_j", "grid_k",
            "raw_component_index", "intersection_component_index",
            "x_m", "y_m", "z_m", "x_mm_from_m", "y_mm_from_m", "z_mm_from_m",
            "ledger_centroid_x_mm", "ledger_centroid_y_mm", "ledger_centroid_z_mm",
            "position_matches_raw_ledger_centroid",
            "analytic_annulus_and_slab_membership", "point_volume_m3",
            "point_volume_converted_mm3", "component_raw_volume_mm3",
            "volume_delta_mm3", "mass_kg", "rho_times_volume_kg",
            "mass_delta_kg", "classifier", "candidate_point_status",
            "source_semantic_region_status", "surface_parent_mapping_status",
        ],
        point_rows,
    )
    figures = make_figures(analytic_rows, points, cell_rows)
    all_checks_pass = all(row["passed"] for row in checks)
    report = {
        "schema": "pgfrag.forge_volume_ingress_independent_pilot_audit.v1",
        "status": (
            "GEOMETRY_DIAGNOSTIC_AUDITED_NOT_PROMOTED"
            if all_checks_pass and no_promotion
            else "AUDIT_DISCREPANCY_NOT_PROMOTED"
        ),
        "expected_geometry_source": (
            "independent analytic R10/R8 disk-minus-disk 1D arc integrals; "
            "producer volumes were not used as expected values"
        ),
        "analytic_roi": analytic_roi,
        "centroid_comparison": {
            "positive_cell_tolerance_mm": CELL_CENTROID_TOL_MM,
            "axes": ["x", "y", "z"],
            "max_absolute_delta_by_axis_mm": {
                "x": max_centroid_delta_by_axis[0],
                "y": max_centroid_delta_by_axis[1],
                "z": max_centroid_delta_by_axis[2],
            },
            "max_absolute_delta_mm": max_centroid_abs_delta_mm,
            "mismatched_positive_cell_indices": centroid_error_cells,
        },
        "analytic_membership": {
            "strict_radial_interval_mm": [
                INNER_RADIUS_MM,
                OUTER_RADIUS_MM,
            ],
            "source_z_interval_mm": [0.0, 10.0],
            "point_ids_outside_analytic_material_or_cell_slab": (
                analytic_membership_errors
            ),
        },
        "input_freeze_verification": frozen,
        "raw_counts": {
            "ledger_cells": len(cells),
            "point_rows": len(points),
            "positive_component_rows": positive_components,
            "classifier_counts": class_counts,
            "topology_counts": topology_counts,
            "cell_status_counts": {
                str(status): sum(c.get("cell_status") == status for c in cells)
                for status in sorted({str(c.get("cell_status")) for c in cells})
            },
        },
        "whole_roi_reconciliation": whole,
        "producer_summary_arithmetic_crosscheck_only": {
            "summary_pointset_count": summary.get("local_geometry_results", {}).get(
                "pointset_count"
            ),
            "raw_pointset_count": len(points),
            "summary_cell_count": summary.get("local_geometry_results", {}).get(
                "processed_cells"
            ),
            "raw_cell_count": len(cells),
            "summary_raw_reference_mm3": summary.get(
                "local_geometry_results", {}
            ).get("reference_roi_common_volume_mm3"),
            "independent_analytic_reference_mm3": analytic_roi[
                "direct_volume_mm3"
            ],
            "summary_raw_cell_sum_mm3": summary.get(
                "local_geometry_results", {}
            ).get("sum_96_cell_common_volumes_mm3"),
            "independent_raw_cell_sum_mm3": raw_cell_sum,
            "summary_values_are_not_used_as_expected_geometry": True,
        },
        "finalization_and_promotion": {
            "final_status": summary.get("pilot_final_status"),
            "finalization_failure": finalization,
            "parent_receipt": {
                "status": parent.get("status"),
                "candidate_run_promoted": promoted,
            },
            "ingress_report_present": (STAGE / "ingress_report.json").exists(),
            "candidate_file_retained_but_not_promoted": (
                STAGE / "material_points.json"
            ).is_file() and promoted is False,
        },
        "checks": checks,
        "all_checks_pass": all_checks_pass,
        "outputs": {
            "cell_csv": CELL_CSV.name,
            "point_csv": POINT_CSV.name,
            **figures,
        },
        "elapsed_seconds": time.perf_counter() - started,
        "qualification": {
            "local_geometry_diagnostic_only": True,
            "scientific_material_calibration": "UNKNOWN_UNVALIDATED",
            "whole_mother_coverage": "UNKNOWN_NOT_RUN",
            "native_integration": "NOT_CONNECTED",
            "solver_or_physics": "NOT_RUN",
        },
    }
    with OUTPUT_JSON.open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(report, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps({
        "status": report["status"],
        "all_checks_pass": all_checks_pass,
        "failed_checks": [c["name"] for c in checks if not c["passed"]],
        "output_json": OUTPUT_JSON.name,
        "outputs": report["outputs"],
        "elapsed_seconds": report["elapsed_seconds"],
    }, sort_keys=True))
    return report


if __name__ == "__main__":
    try:
        run_audit()
    except Exception as exc:
        print(json.dumps({
            "status": "AUDIT_NOT_COMPLETED",
            "error_type": type(exc).__name__,
            "reason": str(exc),
        }, sort_keys=True), file=sys.stderr)
        raise SystemExit(2)
