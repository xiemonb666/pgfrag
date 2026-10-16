#!/usr/bin/env python3
"""One-shot, zero-step public FEEngine layout probe for the frozen coupon.

This process initializes the accepted extrinsic model, assembles initial
lumped mass exactly once, registers automatic insertion eligibility, and then
reads public shape/count/connectivity/normal APIs only. It never applies a BC,
checks cohesive stress, solves a step, or performs an observation assembly.
"""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import sys
import time
import traceback

import numpy as np


ROOT = Path(__file__).resolve().parents[4]
CASE_ROOT = ROOT / "relicweave-synlab/tmp/review_20261004/implementation/cohesive_event_probe_20261004/coarse654_20261004_01"
INPUT = CASE_ROOT / "snapshot/input"
MANIFEST = CASE_ROOT / "snapshot/manifest.json"
WHEEL = ROOT / "third_party/akantu-runtime/downloads/akantu-5.0.7.post1-cp312-cp312-manylinux_2_17_x86_64.manylinux2014_x86_64.whl"
DATA = ROOT / "pgfrag/data/continuation_20261016/extrinsic_onset_refinement/shape_only_preflight"
MSH_SHA256 = "b98402c5a7b7af36b8e5e8a2a90c7db689cb1af86a089c8e1b9a0dbf07ba2638"
MATERIAL_SHA256 = "c300c3091ac323eba756f2d7bfbf55e888c68d6b10f5fcb9763c0e6b6194305f"
MANIFEST_SHA256 = "21b563069550e2091474fe94877bd3ceaae094ac4c92aed3900ecdd94702dd68"
WHEEL_SHA256 = "f618c31bb4cdf23d86e74a10f7e94dac3dae7b0ae2cfd057d352c1b156c0654b"
OUTPUT_CAP_BYTES = 16 * 1024 * 1024


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, value) -> None:
    with path.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


def array_copy(value) -> np.ndarray:
    return np.asarray(value).copy()


def exact_source_facet_map(mesh, aka, source_table, source_nodes: np.ndarray, check_mask: np.ndarray) -> dict[str, object]:
    facet_mesh = mesh.getMeshFacets()
    facet_nodes = array_copy(facet_mesh.getNodes())
    if facet_nodes.shape != source_nodes.shape or not np.array_equal(facet_nodes, source_nodes):
        raise RuntimeError("facet-mesh node rows are not an exact match to source/native node rows")
    connectivity = array_copy(facet_mesh.getConnectivity(aka._triangle_3)).astype(np.int64, copy=False)
    if connectivity.ndim != 2 or connectivity.shape[1] != 3:
        raise RuntimeError(f"native triangle connectivity shape is unexpected: {connectivity.shape}")
    if check_mask.ndim != 1 or len(check_mask) != len(connectivity):
        raise RuntimeError(f"candidate mask/connectivity counts differ: {check_mask.shape} / {connectivity.shape}")
    if connectivity.size and (connectivity.min() < 0 or connectivity.max() >= len(source_nodes)):
        raise RuntimeError("native facet connectivity references an out-of-range node row")

    source_id_by_row = tuple(map(int, source_table.source_node_ids))
    all_records = []
    for slot, node_rows in enumerate(connectivity):
        ids = tuple(sorted(source_id_by_row[int(row)] for row in node_rows))
        facet = source_table.facets_by_node_ids.get(ids)
        if facet is None:
            raise RuntimeError(f"native facet slot {slot} does not map to an exact source T4 facet: {ids}")
        points = source_nodes[node_rows]
        cross = np.cross(points[1] - points[0], points[2] - points[0])
        norm = float(np.linalg.norm(cross))
        if not math.isfinite(norm) or norm <= 0:
            raise RuntimeError(f"native facet slot {slot} is degenerate")
        all_records.append({
            "native_triangle_slot": int(slot),
            "source_node_ids_sorted": list(ids),
            "source_facet_identity": facet.identity,
            "source_facet_classification": facet.classification,
            "native_connectivity_node_rows_ordered": list(map(int, node_rows)),
            "area_m2_from_exact_native_coordinates": 0.5 * norm,
            "geometry_unit_normal_from_connectivity_order": (cross / norm).tolist(),
        })
    slots = np.flatnonzero(check_mask).astype(int).tolist()
    candidates = [all_records[slot] for slot in slots]
    candidate_fids = sorted(str(item["source_facet_identity"]) for item in candidates)
    expected_fids = sorted(source_table.eligible_internal_facets)
    if candidate_fids != expected_fids:
        raise RuntimeError(f"native eligible facets differ from frozen FIDs: {candidate_fids} vs {expected_fids}")
    if len(candidate_fids) != len(set(candidate_fids)) or any(item["source_facet_classification"] != "internal" for item in candidates):
        raise RuntimeError("native eligible FID mapping is duplicated or not internal")
    return {
        "facet_mesh": facet_mesh,
        "connectivity": connectivity,
        "slot_records": all_records,
        "candidate_slots": slots,
        "candidate_records": candidates,
        "expected_fids": expected_fids,
    }


def normalize_shape_layout(raw: np.ndarray, nfacets: int, nqp: int, nnodes: int, selected_slots: list[int]) -> dict[str, object]:
    """Normalize only a supported FEEngine storage layout; always retain raw separately."""
    raw = np.asarray(raw)
    result: dict[str, object] = {
        "raw_shape": list(raw.shape),
        "raw_dtype": str(raw.dtype),
        "raw_c_contiguous": bool(raw.flags.c_contiguous),
        "raw_order_values_preserved": True,
    }
    normalized: np.ndarray | None = None
    layout = "UNKNOWN_RAW_ONLY"
    if raw.ndim == 3 and raw.shape == (nfacets, nqp, nnodes):
        normalized = raw.copy()
        layout = "PER_FACET_QP_NODE"
    elif raw.ndim == 2 and raw.shape == (nqp, nnodes):
        normalized = np.broadcast_to(raw[None, :, :], (nfacets, nqp, nnodes)).copy()
        layout = "SHARED_TYPE_LEVEL_QP_NODE"
    elif raw.ndim == 2 and raw.shape == (nfacets * nqp, nnodes):
        grouped = raw.reshape(nfacets, nqp, nnodes)
        if np.allclose(grouped, grouped[0:1], rtol=0.0, atol=1e-14):
            normalized = grouped.copy()
            layout = "FLAT_ELEMENT_MAJOR_REPEATED_QP_NODE"
        else:
            result["normalization_reason"] = "flattened storage does not repeat one exact type-level rule per facet slot"
    else:
        result["normalization_reason"] = "raw layout does not match a documented count/slot-compatible public shape layout"
    result["layout"] = layout
    if normalized is not None:
        result["normalized_shape"] = list(normalized.shape)
        result["selected_candidate_shape_functions_in_native_qp_order"] = normalized[selected_slots].tolist()
        result["qp_slot_pairing"] = "one row of getShapes per QP in API-returned order; candidate facet slots remain ascending ElementTypeMap order"
        result["normalized_array"] = normalized[selected_slots].copy()
    return result


def check_symmetric_triangle_rule(candidate_shapes: np.ndarray | None) -> dict[str, object]:
    if candidate_shapes is None:
        return {"status": "UNKNOWN_RAW_ONLY", "reason": "shape layout could not be safely normalized"}
    shapes = np.asarray(candidate_shapes, dtype=np.float64)
    expected = np.asarray([
        [2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0],
        [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0],
        [1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0],
    ])
    if shapes.ndim != 3 or shapes.shape[1:] != (3, 3):
        return {"status": "UNKNOWN_RAW_ONLY", "reason": "normalized candidate shape is not (facet,3,3)", "shape": list(shapes.shape)}
    failures = []
    evidence = []
    for candidate_index, local in enumerate(shapes):
        partition_ok = bool(np.allclose(local.sum(axis=1), 1.0, rtol=0.0, atol=1e-12))
        each_point_ok = all(
            any(np.allclose(np.sort(row), np.sort(point), rtol=0.0, atol=1e-12) for point in expected)
            for row in local
        )
        each_permutation_once = all(
            sum(np.allclose(row, point, rtol=0.0, atol=1e-12) for row in local) == 1
            for point in expected
        )
        ok = bool(partition_ok and each_point_ok and each_permutation_once)
        if not ok:
            failures.append(candidate_index)
        evidence.append({
            "candidate_index_in_native_slot_order": candidate_index,
            "shape_values_in_native_qp_order": local.tolist(),
            "partition_of_unity_ok": partition_ok,
            "symmetric_three_point_rule_ok": bool(each_point_ok and each_permutation_once),
        })
    return {
        "status": "SYMMETRIC_THREE_POINT_TRI3_SHAPE_RULE_VERIFIED" if not failures else "UNKNOWN_RAW_ONLY",
        "failed_candidate_indices": failures,
        "candidate_evidence": evidence,
    }


def check_normals(candidate_normals: np.ndarray, candidate_records: list[dict[str, object]]) -> dict[str, object]:
    normals = np.asarray(candidate_normals, dtype=np.float64)
    if normals.ndim != 3 or normals.shape != (len(candidate_records), 3, 3):
        return {"status": "UNKNOWN_RAW_ONLY", "reason": "native normal array does not normalize to (candidate facet,3,xyz)", "shape": list(normals.shape)}
    records = []
    failures = []
    for index, item in enumerate(candidate_records):
        local = normals[index]
        geom = np.asarray(item["geometry_unit_normal_from_connectivity_order"], dtype=np.float64)
        lengths = np.linalg.norm(local, axis=1)
        dots = local @ geom
        unit_ok = bool(np.allclose(lengths, 1.0, rtol=0.0, atol=1e-8))
        parallel_ok = bool(np.allclose(np.abs(dots), 1.0, rtol=0.0, atol=1e-8))
        orientation_ok = bool(np.all(np.sign(dots) == np.sign(dots[0])))
        ok = bool(unit_ok and parallel_ok and orientation_ok)
        if not ok:
            failures.append(index)
        records.append({
            "native_triangle_slot": item["native_triangle_slot"],
            "source_facet_identity": item["source_facet_identity"],
            "normals_native_qp_order": local.tolist(),
            "normal_lengths": lengths.tolist(),
            "dot_with_exact_connectivity_normal": dots.tolist(),
            "orientation_sign_against_connectivity_order": int(np.sign(dots[0])),
            "unit_ok": unit_ok,
            "parallel_ok": parallel_ok,
            "one_orientation_across_qps_ok": orientation_ok,
        })
    return {
        "status": "NATIVE_NORMALS_MATCH_EXACT_TRIANGLE_GEOMETRY" if not failures else "UNKNOWN_RAW_ONLY",
        "failed_candidate_indices": failures,
        "candidate_evidence": records,
    }


def run() -> int:
    DATA.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    status_path = DATA / "preflight_status.json"
    write_json(status_path, {
        "status": "RUNNING_ZERO_STEP_SHAPE_ONLY",
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "solver_steps": 0,
        "cohesive_checks": 0,
        "observation_assemblies": 0,
    })
    try:
        for name, path, expected in (
            ("source_msh", INPUT / "block_2x2x4mm.msh", MSH_SHA256),
            ("material_dat", INPUT / "material.dat", MATERIAL_SHA256),
            ("manifest", MANIFEST, MANIFEST_SHA256),
            ("runtime_wheel", WHEEL, WHEEL_SHA256),
        ):
            if sha256_file(path) != expected:
                raise RuntimeError(f"frozen {name} hash mismatch")

        import akantu as aka

        if getattr(aka, "__version__", None) != "5.0.7-post1":
            raise RuntimeError(f"unexpected runtime version {getattr(aka, '__version__', None)!r}")
        sys.path.insert(0, str(ROOT / "pgfrag/src"))
        from pgfrag.cohesive_mesh import read_gmsh22
        from pgfrag.facet_ancestry import build_source_facet_table

        parsed = read_gmsh22(INPUT / "block_2x2x4mm.msh")
        source_nodes = np.asarray(parsed["nodes"], dtype=np.float64)
        source_tets = np.asarray([row["nodes"] for row in parsed["tetrahedra"]], dtype=np.int64)
        source_node_ids = list(map(int, parsed["node_tags"]))
        source_cell_ids = [int(row["element_id"]) for row in parsed["tetrahedra"]]
        source_interface = [
            tuple(map(int, row["node_tags"]))
            for row in parsed["triangles"] if row["physical_name"] == "interface"
        ]
        source_table = build_source_facet_table(
            source_nodes, source_tets, source_node_ids, source_cell_ids,
            MSH_SHA256, eligible_internal_facets=source_interface,
        )
        if len(source_tets) != 96 or len(source_nodes) != 45 or len(source_interface) != 8:
            raise RuntimeError("the frozen source mesh did not match 96 T4, 45 nodes, and eight interface Tri3 facets")

        aka.parseInput(str(INPUT / "material.dat"))
        mesh = aka.Mesh(3)
        mesh.read(str(INPUT / "block_2x2x4mm.msh"))
        model = aka.SolidMechanicsModelCohesive(mesh)
        model.initFull(_analysis_method=aka._explicit_lumped_mass, _is_extrinsic=True)
        model.assembleMassLumped()
        inserter = model.getElementInserter()
        inserter.addPhysicalSurface("interface")
        model.updateAutomaticInsertion()

        native_nodes = array_copy(mesh.getNodes())
        native_tets = array_copy(mesh.getConnectivity(aka._tetrahedron_4)).astype(np.int64, copy=False)
        check_mask = array_copy(inserter.getCheckFacets()(aka._triangle_3)).reshape(-1).astype(bool)
        inserted_mask = array_copy(inserter.getInsertionFacets()(aka._triangle_3)).reshape(-1).astype(bool)
        if not np.array_equal(native_nodes, source_nodes) or not np.array_equal(native_tets, source_tets):
            raise RuntimeError("native node/T4 row and corner order differs from frozen source mesh")
        if int(check_mask.sum()) != 8 or int(inserted_mask.sum()) != 0:
            raise RuntimeError(f"native eligibility counts are not 8/0: {int(check_mask.sum())}/{int(inserted_mask.sum())}")

        facet_map = exact_source_facet_map(mesh, aka, source_table, native_nodes, check_mask)
        facet_engine = model.getFEEngine("FacetsFEEngine")
        nqp = int(facet_engine.getNbIntegrationPoints(aka._triangle_3, aka._not_ghost))
        raw_shapes = array_copy(facet_engine.getShapes(aka._triangle_3, aka._not_ghost, 0))
        shape_layout = normalize_shape_layout(raw_shapes, len(facet_map["connectivity"]), nqp, 3, facet_map["candidate_slots"])
        candidate_shapes = shape_layout.get("normalized_array")
        shape_rule = check_symmetric_triangle_rule(candidate_shapes)

        raw_normals = array_copy(facet_engine.getNormalsOnIntegrationPoints(aka._triangle_3, aka._not_ghost))
        grouped_normals: np.ndarray | None = None
        normal_layout: dict[str, object]
        if raw_normals.ndim == 2 and raw_normals.shape == (len(facet_map["connectivity"]) * nqp, 3):
            grouped_normals = raw_normals.reshape(len(facet_map["connectivity"]), nqp, 3)
            normal_layout = {
                "status": "ELEMENT_MAJOR_QP_GROUPING_FROM_NQP_AND_EXACT_TRIANGLE_SLOT_COUNT",
                "raw_shape": list(raw_normals.shape),
                "normalized_shape": list(grouped_normals.shape),
                "grouping": "contiguous nqp rows per ascending triangle_3 ElementTypeMap slot",
            }
        elif raw_normals.ndim == 3 and raw_normals.shape == (len(facet_map["connectivity"]), nqp, 3):
            grouped_normals = raw_normals.copy()
            normal_layout = {
                "status": "PER_FACET_QP_XYZ",
                "raw_shape": list(raw_normals.shape),
                "normalized_shape": list(grouped_normals.shape),
            }
        else:
            normal_layout = {"status": "UNKNOWN_RAW_ONLY", "raw_shape": list(raw_normals.shape), "nqp": nqp}
        if grouped_normals is not None and nqp == 3:
            selected_normals = grouped_normals[check_mask].copy()
            normals_check = check_normals(selected_normals, facet_map["candidate_records"])
        else:
            selected_normals = None
            normals_check = {"status": "UNKNOWN_RAW_ONLY", "reason": "native normal grouping or QP count not supported by the observed layout"}

        candidate_area_checks = []
        for item in facet_map["candidate_records"]:
            area = float(item["area_m2_from_exact_native_coordinates"])
            candidate_area_checks.append({
                "native_triangle_slot": item["native_triangle_slot"],
                "source_facet_identity": item["source_facet_identity"],
                "area_m2": area,
                "matches_actual_frozen_interface_area_0_5e_minus6": math.isclose(area, 0.5e-6, rel_tol=0.0, abs_tol=1e-18),
            })
        all_area_checks = all(item["matches_actual_frozen_interface_area_0_5e_minus6"] for item in candidate_area_checks)
        measure_verified = (
            all_area_checks
            and shape_rule.get("status") == "SYMMETRIC_THREE_POINT_TRI3_SHAPE_RULE_VERIFIED"
            and normals_check.get("status") == "NATIVE_NORMALS_MATCH_EXACT_TRIANGLE_GEOMETRY"
        )

        raw_arrays = {
            "source_nodes": source_nodes,
            "native_nodes": native_nodes,
            "source_tetra_connectivity": source_tets,
            "native_tetra_connectivity": native_tets,
            "facet_triangle_connectivity_all_slots": facet_map["connectivity"],
            "eligible_candidate_check_mask": check_mask,
            "insertion_mask_initial": inserted_mask,
            "facet_triangle_getShapes_raw": raw_shapes,
            "facet_triangle_normals_raw": raw_normals,
        }
        if candidate_shapes is not None:
            raw_arrays["candidate_shapes_normalized_native_qp_order"] = candidate_shapes
        if selected_normals is not None:
            raw_arrays["candidate_normals_normalized_native_qp_order"] = selected_normals
        with (DATA / "zero_step_native_arrays.npz").open("xb") as stream:
            np.savez_compressed(stream, **raw_arrays)

        summary = {
            "schema": "pgfrag.extrinsic_onset_refinement.shape_only_preflight.v1",
            "status": "VERIFIED_COUPON_AFFINE_TRI3" if measure_verified else "CAPTURED_ZERO_STEP_LAYOUT_WITH_UNKNOWN_MEASURE",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "elapsed_wall_seconds": time.monotonic() - started,
            "runtime": {
                "akantu_distribution": getattr(aka, "__version__", "UNKNOWN"),
                "python_executable": sys.executable,
                "python_version": platform.python_version(),
                "platform": platform.platform(),
                "thread_environment": {name: os.environ.get(name) for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")},
            },
            "frozen_input_hashes": {
                "source_msh_sha256": sha256_file(INPUT / "block_2x2x4mm.msh"),
                "material_dat_sha256": sha256_file(INPUT / "material.dat"),
                "manifest_sha256": sha256_file(MANIFEST),
                "wheel_sha256": sha256_file(WHEEL),
            },
            "native_initialization": {
                "mesh_elements": int(mesh.getNbElement(type=aka._tetrahedron_4)),
                "mesh_nodes": int(mesh.getNbNodes()),
                "cohesive_elements": int(mesh.getNbElement(type=aka._cohesive_3d_6)),
                "eligible_candidate_facets": int(check_mask.sum()),
                "inserted_facets": int(inserted_mask.sum()),
                "analysis_method": "_explicit_lumped_mass",
                "insertion_mode": "extrinsic",
                "initial_mass_assemblies": 1,
                "solver_steps": 0,
                "cohesive_checks": 0,
                "BC_calls": 0,
                "observation_assemblies": 0,
            },
            "public_api_layout": {
                "FEEngine.getNbIntegrationPoints": {"type": "_triangle_3", "value": nqp},
                "FEEngine.getShapes": {
                    "method": "public getShapes(_triangle_3,_not_ghost,0)",
                    "layout": {key: value for key, value in shape_layout.items() if key != "normalized_array"},
                    "symmetric_rule_check": shape_rule,
                },
                "FEEngine.getNormalsOnIntegrationPoints": {
                    "method": "public getNormalsOnIntegrationPoints(_triangle_3,_not_ghost)",
                    "layout": normal_layout,
                    "candidate_orientation_check": normals_check,
                },
                "Mesh.getMeshFacets_and_getConnectivity": {
                    "status": "EXACT_NODE_ROW_AND_STABLE_FID_MAP",
                    "native_triangle_slot_count": int(len(facet_map["connectivity"])),
                    "candidate_slots_ascending": facet_map["candidate_slots"],
                    "candidate_records_in_native_slot_order": facet_map["candidate_records"],
                    "all_slot_source_fids": [item["source_facet_identity"] for item in facet_map["slot_records"]],
                    "candidate_fids_equal_frozen_interface_set": True,
                    "connectivity_array_shape": list(facet_map["connectivity"].shape),
                    "facet_nodes_exactly_equal_source_nodes": True,
                    "native_tetra_connectivity_exactly_equal_source": True,
                },
            },
            "coupon_measure": {
                "status": "VERIFIED_COUPON_AFFINE_TRI3" if measure_verified else "UNKNOWN_RAW_ONLY",
                "scope": "only the eight planar affine Tri3 facets in this exact frozen coupon",
                "derivation": "exact native facet connectivity and coordinates; exact stable source FIDs; supported native three-point symmetric shape table; public native normals unit and parallel to each exact facet geometry normal; area/3 per point",
                "weights_source": "coupon-specific geometric measure; no per-QP weight/Jacobian getter is claimed",
                "candidate_area_checks": candidate_area_checks,
                "qp_weight_records": [
                    {
                        "native_triangle_slot": item["native_triangle_slot"],
                        "source_facet_identity": item["source_facet_identity"],
                        "area_m2": item["area_m2_from_exact_native_coordinates"],
                        "qp_weights_m2_native_order": [item["area_m2_from_exact_native_coordinates"] / 3.0] * 3 if measure_verified else None,
                    }
                    for item in facet_map["candidate_records"]
                ],
                "normal_orientation": "native sign is recorded separately per FID; no MSH winding orientation is assumed",
            },
            "raw_arrays_path": "pgfrag/data/continuation_20261016/extrinsic_onset_refinement/shape_only_preflight/zero_step_native_arrays.npz",
            "limits": {
                "output_cap_bytes": OUTPUT_CAP_BYTES,
                "wall_cap_seconds": 60,
                "native_solver_evolution": "none",
                "no_retries": True,
            },
        }
        write_json(DATA / "shape_only_preflight.json", summary)
        total_bytes = sum(path.stat().st_size for path in DATA.iterdir() if path.is_file())
        if total_bytes > OUTPUT_CAP_BYTES:
            raise RuntimeError(f"preflight output {total_bytes} bytes exceeds 16 MiB cap")
        write_json(DATA / "preflight_status_final.json", {
            "status": "COMPLETED_ZERO_STEP_SHAPE_ONLY",
            "elapsed_wall_seconds": summary["elapsed_wall_seconds"],
            "output_bytes": total_bytes,
            "output_cap_bytes": OUTPUT_CAP_BYTES,
            "solver_steps": 0,
            "cohesive_checks": 0,
            "observation_assemblies": 0,
            "measure_status": summary["coupon_measure"]["status"],
        })
        print(json.dumps({"status": summary["status"], "measure_status": summary["coupon_measure"]["status"], "output_bytes": total_bytes, "elapsed_wall_seconds": summary["elapsed_wall_seconds"]}, sort_keys=True))
        return 0
    except BaseException as error:
        failure = {
            "status": "FAILED_ZERO_STEP_SHAPE_ONLY",
            "elapsed_wall_seconds": time.monotonic() - started,
            "error": repr(error),
            "traceback": traceback.format_exc(),
            "solver_steps": 0,
            "cohesive_checks": 0,
            "observation_assemblies": 0,
        }
        write_json(DATA / "preflight_failure.json", failure)
        raise


if __name__ == "__main__":
    raise SystemExit(run())
