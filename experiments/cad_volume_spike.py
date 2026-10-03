"""Bounded CAD-to-volume quadrature candidate for the existing VF STEP asset.

This experiment imports the manifest-selected canonical STEP into Gmsh's
OpenCASCADE kernel, creates at most two deliberately coarse first-order
tetrahedral meshes, and emits one centroid/positive-volume quadrature candidate.
It is an offline geometry spike; it does not create a PG job or claim fragment
or fracture support.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import os
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any

import cadquery as cq
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.patches import Rectangle
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from OCP.BRepBndLib import BRepBndLib
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.IntCurvesFace import IntCurvesFace_ShapeIntersector
from OCP.TopAbs import TopAbs_IN, TopAbs_ON, TopAbs_OUT
from OCP.Bnd import Bnd_Box
from OCP.gp import gp_Dir, gp_Lin, gp_Pnt


EXPERIMENT_ID = "pg.cad-volume-centroid-quadrature-candidate/0.1.0"
MAX_VOLUME_CELLS = 2500
MESH_SIZE_TRIALS_MM = (8.0, 6.0)
PROPOSED_VOLUME_RELATIVE_ERROR_LIMIT = 0.01
CLASSIFIER_TOLERANCE_MM = 1e-7
INTERIOR_SAMPLE_COUNT = 512
INTERIOR_SAMPLE_SEED = 20261003
RHO_STATUS = "UNVALIDATED_PROTOTYPE"

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_ASSET_DIR = (
    ROOT
    / "relicweave-synlab"
    / "tmp"
    / "review_20261003"
    / "implementation"
    / "vf"
    / "artifacts"
    / "four_bands_brep_semantics1_final"
)
DEFAULT_OUTPUT_DIR = (
    ROOT
    / "relicweave-synlab"
    / "tmp"
    / "review_20261003"
    / "implementation"
    / "pg"
    / "volume_candidate"
)
GMSH_SDK_LIB = (
    ROOT
    / "third_party"
    / "pgfrag-geometry-runtime"
    / "sdk"
    / "gmsh-4.15.2-Windows64-sdk"
    / "lib"
)
GMSH_SDK_ARCHIVE = (
    ROOT / "third_party" / "pgfrag-geometry-runtime" / "gmsh-4.15.2-Windows64-sdk.zip"
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_gmsh() -> Any:
    if not GMSH_SDK_LIB.is_dir():
        raise FileNotFoundError(
            f"isolated Gmsh SDK library directory not found: {GMSH_SDK_LIB}"
        )
    # The Gmsh Python module and DLL live together in the official SDK. Keep the
    # DLL search path scoped to this process; no system PATH or other runtime is
    # modified.
    if hasattr(os, "add_dll_directory"):
        os.add_dll_directory(str(GMSH_SDK_LIB))
    sys.path.insert(0, str(GMSH_SDK_LIB))
    import gmsh  # type: ignore[import-not-found]

    return gmsh


def _manifest_source(asset_dir: Path) -> tuple[dict[str, Any], Path, str]:
    manifest_path = asset_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    relative_step = str(manifest["geometry"]["solid_step"])
    step_path = (asset_dir / relative_step).resolve(strict=True)
    expected = manifest["files"][relative_step]
    actual_hash = _sha256(step_path)
    actual_bytes = step_path.stat().st_size
    if actual_hash.lower() != str(expected["sha256"]).lower():
        raise ValueError(
            f"canonical STEP SHA-256 mismatch: manifest={expected['sha256']} actual={actual_hash}"
        )
    if actual_bytes != int(expected["bytes"]):
        raise ValueError(
            f"canonical STEP byte-count mismatch: manifest={expected['bytes']} actual={actual_bytes}"
        )
    return manifest, step_path, actual_hash


def _bbox_mm(shape: cq.Shape) -> dict[str, float]:
    box = Bnd_Box()
    BRepBndLib.AddOptimal_s(shape.wrapped, box, False, False)
    xmin, ymin, zmin, xmax, ymax, zmax = box.Get()
    return {
        "xmin_mm": float(xmin),
        "ymin_mm": float(ymin),
        "zmin_mm": float(zmin),
        "xmax_mm": float(xmax),
        "ymax_mm": float(ymax),
        "zmax_mm": float(zmax),
        "xlen_mm": float(xmax - xmin),
        "ylen_mm": float(ymax - ymin),
        "zlen_mm": float(zmax - zmin),
    }


def _classifier(solid: cq.Shape) -> BRepClass3d_SolidClassifier:
    classifier = BRepClass3d_SolidClassifier()
    classifier.Load(solid.wrapped)
    return classifier


def _classify_one(
    classifier: BRepClass3d_SolidClassifier,
    point: np.ndarray | tuple[float, float, float],
) -> str:
    xyz = np.asarray(point, dtype=np.float64)
    if xyz.shape != (3,) or not np.isfinite(xyz).all():
        return "UNKNOWN"
    try:
        classifier.Perform(gp_Pnt(float(xyz[0]), float(xyz[1]), float(xyz[2])), CLASSIFIER_TOLERANCE_MM)
        state = classifier.State()
    except Exception:
        return "UNKNOWN"
    if state == TopAbs_IN:
        return "IN"
    if state == TopAbs_OUT:
        return "OUT"
    if state == TopAbs_ON:
        return "UNKNOWN"
    return "UNKNOWN"


def _classify_points(solid: cq.Shape, points: np.ndarray) -> tuple[list[str], dict[str, int]]:
    checker = _classifier(solid)
    states = [_classify_one(checker, point) for point in points]
    counts = Counter(states)
    return states, {key: int(counts.get(key, 0)) for key in ("IN", "OUT", "UNKNOWN")}


def _tetra_faces(connectivity: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    local_faces = ((0, 1, 2), (0, 3, 1), (0, 2, 3), (1, 3, 2))
    counts: Counter[tuple[int, int, int]] = Counter()
    oriented: dict[tuple[int, int, int], tuple[int, int, int]] = {}
    for tet in connectivity:
        for face in local_faces:
            ids = tuple(int(tet[i]) for i in face)
            key = tuple(sorted(ids))
            counts[key] += 1
            oriented.setdefault(key, ids)
    boundary = np.asarray(
        [oriented[key] for key, count in counts.items() if count == 1], dtype=np.int64
    ).reshape((-1, 3))
    nonmanifold_face_count = sum(1 for count in counts.values() if count > 2)
    return boundary, np.asarray([nonmanifold_face_count], dtype=np.int64)


def _extract_tetra_mesh(gmsh: Any, volume_tag: int) -> dict[str, Any]:
    element_types, element_tags_by_type, node_tags_by_type = gmsh.model.mesh.getElements(3, volume_tag)
    tetra_type = int(gmsh.model.mesh.getElementType("tetrahedron", 1))
    if not element_types:
        raise ValueError("Gmsh produced no volume elements")
    unsupported = [int(kind) for kind in element_types if int(kind) != tetra_type]
    if unsupported:
        raise ValueError(f"expected linear tetrahedra only; received volume element types {unsupported}")

    source_cell_ids = np.concatenate(
        [np.asarray(tags, dtype=np.int64).reshape(-1) for tags in element_tags_by_type]
    )
    raw_node_tags = np.concatenate(
        [np.asarray(tags, dtype=np.int64).reshape(-1, 4) for tags in node_tags_by_type], axis=0
    )
    if len(source_cell_ids) != len(raw_node_tags):
        raise ValueError("Gmsh element ID and tetra connectivity counts disagree")
    source_node_ids, node_coordinates, _ = gmsh.model.mesh.getNodes()
    source_node_ids = np.asarray(source_node_ids, dtype=np.int64).reshape(-1)
    node_coordinates = np.asarray(node_coordinates, dtype=np.float64).reshape((-1, 3))
    node_id_unique = np.unique(source_node_ids).size == source_node_ids.size
    if not node_id_unique:
        raise ValueError("Gmsh returned duplicate source node IDs")
    node_row_by_id = {int(node_id): row for row, node_id in enumerate(source_node_ids)}
    try:
        node_rows = np.asarray(
            [[node_row_by_id[int(tag)] for tag in tet] for tet in raw_node_tags], dtype=np.int64
        )
    except KeyError as error:
        raise ValueError(f"tetra connectivity references an unknown Gmsh node ID {error}") from error
    tetra_xyz = node_coordinates[node_rows]

    edge1 = tetra_xyz[:, 1, :] - tetra_xyz[:, 0, :]
    edge2 = tetra_xyz[:, 2, :] - tetra_xyz[:, 0, :]
    edge3 = tetra_xyz[:, 3, :] - tetra_xyz[:, 0, :]
    signed_volume_mm3 = np.einsum("ij,ij->i", np.cross(edge1, edge2), edge3) / 6.0
    cell_volume_mm3 = np.abs(signed_volume_mm3)
    points_mm = np.mean(tetra_xyz, axis=1)
    finite_cells = np.isfinite(tetra_xyz).all(axis=(1, 2)) & np.isfinite(cell_volume_mm3)
    positive_cells = cell_volume_mm3 > 0.0
    centroid_unique_count = int(np.unique(np.round(points_mm, decimals=10), axis=0).shape[0])
    boundary_faces, nonmanifold_count = _tetra_faces(raw_node_tags)

    return {
        "sourcecellIDs": source_cell_ids,
        "source_node_ids": source_node_ids,
        "tetra_node_ids": raw_node_tags,
        "node_coordinates_mm": node_coordinates,
        "points_mm": points_mm,
        "cell_volume_mm3": cell_volume_mm3,
        "signed_volume_mm3": signed_volume_mm3,
        "boundary_node_ids": boundary_faces,
        "node_id_unique": bool(node_id_unique),
        "cell_id_unique": bool(np.unique(source_cell_ids).size == source_cell_ids.size),
        "connectivity_references_known_nodes": bool(np.isin(raw_node_tags, source_node_ids).all()),
        "finite_cell_count": int(np.count_nonzero(finite_cells)),
        "positive_cell_volume_count": int(np.count_nonzero(positive_cells)),
        "zero_or_negative_unsigned_volume_count": int(np.count_nonzero(~positive_cells)),
        "negative_signed_orientation_count": int(np.count_nonzero(signed_volume_mm3 < 0.0)),
        "centroid_unique_count": centroid_unique_count,
        "boundary_face_count": int(boundary_faces.shape[0]),
        "nonmanifold_tetra_face_count": int(nonmanifold_count[0]),
        "cell_count": int(source_cell_ids.size),
    }


def _sample_interior_points(mesh: dict[str, Any]) -> np.ndarray:
    points = mesh["points_mm"]
    count = min(INTERIOR_SAMPLE_COUNT, int(points.shape[0]))
    selected = np.linspace(0, points.shape[0] - 1, count, dtype=np.int64)
    tetra_ids = mesh["tetra_node_ids"][selected]
    node_row_by_id = {
        int(node_id): row for row, node_id in enumerate(mesh["source_node_ids"])
    }
    tetra_rows = np.asarray(
        [[node_row_by_id[int(node_id)] for node_id in tet] for tet in tetra_ids], dtype=np.int64
    )
    tetra_xyz = mesh["node_coordinates_mm"][tetra_rows]
    rng = np.random.default_rng(INTERIOR_SAMPLE_SEED)
    weights = rng.dirichlet(np.ones(4, dtype=np.float64), size=count)
    return np.einsum("ij,ijk->ik", weights, tetra_xyz)


def _cavity_probe(solid: cq.Shape, mesh: dict[str, Any]) -> dict[str, Any]:
    # This bounded central cylinder lies clear of the vessel wall and handle for
    # the manifest's stated z range. It probes an empty cavity subset, not the
    # entire non-material complement of the solid.
    z_values = (8.0, 16.0, 24.0, 32.0, 40.0, 48.0, 54.0)
    radial_values = (0.0, 4.0, 8.0)
    angles = np.linspace(0.0, 2.0 * math.pi, 8, endpoint=False)
    samples: list[tuple[float, float, float]] = []
    for z in z_values:
        samples.append((0.0, 0.0, z))
        for radius in radial_values[1:]:
            samples.extend(
                (float(radius * math.cos(theta)), float(radius * math.sin(theta)), z)
                for theta in angles
            )
    states, counts = _classify_points(solid, np.asarray(samples, dtype=np.float64))
    centroids = mesh["points_mm"]
    radii = np.hypot(centroids[:, 0], centroids[:, 1])
    in_probe_region = (
        (radii <= 8.0)
        & (centroids[:, 2] >= min(z_values))
        & (centroids[:, 2] <= max(z_values))
    )
    return {
        "method": "native BRepClass3d_SolidClassifier over 119 deterministic finite points in r<=8 mm, z=8..54 mm; this local cylinder is a sampled cavity subset, not a full-complement proof",
        "sample_count": int(len(samples)),
        "sample_state_counts": counts,
        "sample_points_first_12_mm": [list(point) for point in samples[:12]],
        "sample_states_first_12": states[:12],
        "mesh_centroid_count_in_sampled_cavity_region": int(np.count_nonzero(in_probe_region)),
        "sampled_cavity_region_mm": {"radial_limit_mm": 8.0, "z_min_mm": min(z_values), "z_max_mm": max(z_values)},
        "candidate_has_no_centroid_in_sampled_cavity_region": bool(not np.any(in_probe_region)),
    }


def _segment_profile(solid: cq.Shape, first: np.ndarray, second: np.ndarray) -> dict[str, Any]:
    first_state, _ = _classify_points(solid, np.asarray([first], dtype=np.float64))
    second_state, _ = _classify_points(solid, np.asarray([second], dtype=np.float64))
    direction_vector = np.asarray(second, dtype=np.float64) - np.asarray(first, dtype=np.float64)
    length = float(np.linalg.norm(direction_vector))
    result: dict[str, Any] = {
        "endpoint_state_a": first_state[0],
        "endpoint_state_b": second_state[0],
        "segment_length_mm": length,
        "partition_method": "IntCurvesFace_ShapeIntersector trims the finite segment at native OCC face hits; BRepClass3d_SolidClassifier classifies each resulting open-interval midpoint",
        "intervals": [],
    }
    if not math.isfinite(length) or length < 1e-8:
        result.update({"status": "UNKNOWN", "reason": "degenerate_segment"})
        return result
    if first_state[0] != "IN" or second_state[0] != "IN":
        result.update({"status": "UNKNOWN", "reason": "both endpoints must classify IN for this audit"})
        return result

    direction = direction_vector / length
    line = gp_Lin(
        gp_Pnt(float(first[0]), float(first[1]), float(first[2])),
        gp_Dir(float(direction[0]), float(direction[1]), float(direction[2])),
    )
    intersector = IntCurvesFace_ShapeIntersector()
    intersector.Load(solid.wrapped, CLASSIFIER_TOLERANCE_MM)
    try:
        intersector.Perform(line, 0.0, length)
        if not intersector.IsDone():
            result.update({"status": "UNKNOWN", "reason": "OCC segment/shape intersector did not complete"})
            return result
        events = [0.0, length]
        for hit_index in range(1, int(intersector.NbPnt()) + 1):
            parameter = float(intersector.WParameter(hit_index))
            if math.isfinite(parameter) and 1e-7 < parameter < length - 1e-7:
                events.append(parameter)
    except Exception as error:
        result.update({"status": "UNKNOWN", "reason": f"OCC segment/shape intersection error: {error}"})
        return result

    events.sort()
    unique_events: list[float] = []
    for event in events:
        if not unique_events or abs(event - unique_events[-1]) > 1e-6:
            unique_events.append(event)
    classifier = _classifier(solid)
    intervals: list[dict[str, Any]] = []
    for start, end in zip(unique_events[:-1], unique_events[1:]):
        if end - start <= 1e-7:
            continue
        midpoint_parameter = 0.5 * (start + end)
        midpoint = np.asarray(first, dtype=np.float64) + midpoint_parameter * direction
        state = _classify_one(classifier, midpoint)
        radius = float(math.hypot(midpoint[0], midpoint[1]))
        is_sampled_cavity = bool(radius <= 8.0 and 8.0 <= midpoint[2] <= 54.0)
        intervals.append(
            {
                "start_parameter_mm": float(start),
                "end_parameter_mm": float(end),
                "midpoint_xyz_mm": midpoint.tolist(),
                "state": state,
                "midpoint_in_sampled_cavity_subset": is_sampled_cavity,
            }
        )
    result["intervals"] = intervals
    states = [record["state"] for record in intervals]
    crosses_cavity = any(
        record["state"] == "OUT" and record["midpoint_in_sampled_cavity_subset"]
        for record in intervals
    )
    if "UNKNOWN" in states or not states:
        status = "UNKNOWN"
    elif crosses_cavity:
        status = "CROSSES_SAMPLED_CAVITY"
    elif all(state == "IN" for state in states):
        status = "FULL_SEGMENT_CAD_INSIDE"
    elif "OUT" in states:
        status = "OUTSIDE_SOLID_OTHER_REGION"
    else:
        status = "UNKNOWN"
    result["status"] = status
    result["intersection_hit_count"] = int(intersector.NbPnt())
    result["unique_partition_parameter_count"] = int(len(unique_events))
    return result


def _bond_segment_audit(solid: cq.Shape, mesh: dict[str, Any], centroid_states: list[str]) -> dict[str, Any]:
    points = mesh["points_mm"]
    radial = np.hypot(points[:, 0], points[:, 1])
    z = points[:, 2]
    inside_ids = np.asarray([i for i, state in enumerate(centroid_states) if state == "IN"], dtype=np.int64)
    central_body = inside_ids[
        (radial[inside_ids] >= 11.0)
        & (radial[inside_ids] <= 23.0)
        & (np.abs(points[inside_ids, 0]) <= 8.0)
        & (z[inside_ids] >= 14.0)
        & (z[inside_ids] <= 48.0)
    ]
    left = central_body[points[central_body, 1] <= -10.0]
    right = central_body[points[central_body, 1] >= 10.0]
    crosswall_records: list[dict[str, Any]] = []
    candidate_pairs: list[tuple[float, int, int]] = []
    for i in left:
        candidates = right[
            (np.abs(points[right, 0] - points[i, 0]) <= 4.0)
            & (np.abs(points[right, 2] - points[i, 2]) <= 4.0)
        ]
        for j in candidates:
            delta = points[int(j)] - points[int(i)]
            score = (
                abs(float(delta[2]))
                + 0.5 * abs(float(delta[0]))
                + 0.05 * abs(float(delta[1]))
            )
            candidate_pairs.append((score, int(i), int(j)))
    candidate_pairs.sort()
    seen: set[tuple[int, int]] = set()
    for _, i, j in candidate_pairs[:120]:
        key = tuple(sorted((i, j)))
        if key in seen:
            continue
        seen.add(key)
        profile = _segment_profile(solid, points[i], points[j])
        if profile["status"] == "CROSSES_SAMPLED_CAVITY":
            crosswall_records.append(
                {
                    "sourcecellID_a": int(mesh["sourcecellIDs"][i]),
                    "sourcecellID_b": int(mesh["sourcecellIDs"][j]),
                    "point_a_mm": points[i].tolist(),
                    "point_b_mm": points[j].tolist(),
                    **profile,
                }
            )
            if len(crosswall_records) >= 3:
                break

    side = central_body[np.abs(points[central_body, 1]) >= 10.0]
    samewall_records: list[dict[str, Any]] = []
    if side.size:
        stride = max(1, side.size // 300)
        for i in side[::stride][:300]:
            distances = np.linalg.norm(points[side] - points[int(i)], axis=1)
            order = np.argsort(distances)
            for candidate_position in order[1:10]:
                j = int(side[candidate_position])
                if not (0.25 < distances[candidate_position] <= 8.0):
                    continue
                profile = _segment_profile(solid, points[int(i)], points[j])
                if profile["status"] == "FULL_SEGMENT_CAD_INSIDE":
                    samewall_records.append(
                        {
                            "sourcecellID_a": int(mesh["sourcecellIDs"][int(i)]),
                            "sourcecellID_b": int(mesh["sourcecellIDs"][j]),
                            "point_a_mm": points[int(i)].tolist(),
                            "point_b_mm": points[j].tolist(),
                            **profile,
                        }
                    )
                    if len(samewall_records) >= 3:
                        break
            if len(samewall_records) >= 3:
                break

    return {
        "status": "PROBE_ONLY_NOT_NATIVE_BOND_GRAPH",
        "native_bond_graph_created": False,
        "horizon_mm": None,
        "endpoint_policy": "both pair endpoints must independently classify IN; endpoint inclusion alone does not classify the intervening segment",
        "segment_classifier": "IntCurvesFace_ShapeIntersector partitions each bounded centroid-to-centroid line at OCC face intersections; BRepClass3d_SolidClassifier labels each interval midpoint IN/OUT/UNKNOWN",
        "crosswall_pair_search_scope": "up to 120 deterministic opposite-y central-body centroid pair candidates; |dx|<=4 mm, |dz|<=4 mm, |x|<=8 mm, 14<=z<=48 mm, both endpoints OCC IN",
        "samewall_pair_search_scope": "deterministic nearest centroid pairs on one central-body y-side with endpoint distance<=8 mm; both endpoints OCC IN",
        "cavity_subset_for_pair_label": "radial<=8 mm and 8<=z<=54 mm, which is also probed separately as an empty CAD subset",
        "crosswall_candidates": crosswall_records,
        "samewall_controls": samewall_records,
        "interpretation": "These finite geometric line probes only demonstrate cavity obstruction for selected endpoint pairs. They do not create, infer, or validate a native Peridynamics bond graph and are not tied to a chosen horizon.",
    }


def _run_mesh_trial(
    gmsh: Any,
    step_path: Path,
    output_dir: Path,
    mesh_size_mm: float,
    index: int,
    cad_shape: cq.Shape,
    cad_volume_mm3: float,
    density_kg_m3: float,
    max_volume_cells: int,
    policy_id: str,
) -> dict[str, Any]:
    trial: dict[str, Any] = {
        "mesh_size_target_mm": mesh_size_mm,
        "max_volume_cells": max_volume_cells,
        "policy_id": policy_id,
        "status": "NOT_STARTED",
        "mesh_file": None,
        "candidate_arrays": None,
    }
    try:
        gmsh.clear()
        gmsh.model.add(f"cad_volume_trial_{index}_h{mesh_size_mm:g}mm")
        gmsh.option.setNumber("General.Terminal", 1)
        gmsh.option.setNumber("Geometry.OCCImportLabels", 0)
        gmsh.option.setNumber("Geometry.OCCScaling", 1.0)
        gmsh.option.setNumber("Mesh.MeshSizeMin", mesh_size_mm)
        gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_size_mm)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)
        gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
        gmsh.option.setNumber("Mesh.Binary", 0)

        imported = gmsh.model.occ.importShapes(str(step_path), True, "step")
        gmsh.model.occ.synchronize()
        volumes = gmsh.model.getEntities(3)
        if len(volumes) != 1 or len(imported) != 1 or int(imported[0][0]) != 3:
            trial.update(
                {
                    "status": "FAIL_IMPORTED_VOLUME_TOPOLOGY",
                    "imported_entities": [[int(dim), int(tag)] for dim, tag in imported],
                    "volume_entities": [[int(dim), int(tag)] for dim, tag in volumes],
                }
            )
            return trial
        volume_tag = int(volumes[0][1])
        gmsh_occ_volume_mm3 = float(gmsh.model.occ.getMass(3, volume_tag))
        gmsh_bbox = gmsh.model.occ.getBoundingBox(3, volume_tag)
        gmsh_bbox_mm = {
            "xmin_mm": float(gmsh_bbox[0]),
            "ymin_mm": float(gmsh_bbox[1]),
            "zmin_mm": float(gmsh_bbox[2]),
            "xmax_mm": float(gmsh_bbox[3]),
            "ymax_mm": float(gmsh_bbox[4]),
            "zmax_mm": float(gmsh_bbox[5]),
        }
        gmsh.model.mesh.generate(3)
        mesh = _extract_tetra_mesh(gmsh, volume_tag)
        cell_count = int(mesh["cell_count"])
        trial.update(
            {
                "gmsh_occ_volume_mm3": gmsh_occ_volume_mm3,
                "gmsh_occ_volume_relative_delta_vs_cadquery": float(
                    (gmsh_occ_volume_mm3 - cad_volume_mm3) / cad_volume_mm3
                ),
                "gmsh_occ_bbox_mm": gmsh_bbox_mm,
                "volume_cell_count": cell_count,
                "source_node_count": int(mesh["source_node_ids"].size),
            }
        )
        if cell_count > max_volume_cells:
            trial.update(
                {
                    "status": "FAIL_CELL_BUDGET",
                    "structured_failure": {
                        "code": "volume_cell_count_exceeds_hard_limit",
                        "observed": cell_count,
                        "limit": max_volume_cells,
                    },
                }
            )
            return trial

        mesh_volume_mm3 = float(np.sum(mesh["cell_volume_mm3"], dtype=np.float64))
        relative_delta = float((mesh_volume_mm3 - cad_volume_mm3) / cad_volume_mm3)
        volume_error = abs(relative_delta)
        node_row_by_id = {
            int(node_id): row for row, node_id in enumerate(mesh["source_node_ids"])
        }
        tetra_rows = np.asarray(
            [[node_row_by_id[int(node_id)] for node_id in tet] for tet in mesh["tetra_node_ids"]],
            dtype=np.int64,
        )
        tetra_xyz = mesh["node_coordinates_mm"][tetra_rows]
        centroid_states, centroid_state_counts = _classify_points(cad_shape, mesh["points_mm"])
        interior_samples = _sample_interior_points(mesh)
        interior_states, interior_state_counts = _classify_points(cad_shape, interior_samples)
        cavity = _cavity_probe(cad_shape, mesh)
        all_cell_mass_kg = mesh["cell_volume_mm3"] * (1e-9 * float(density_kg_m3))
        total_mass_kg = float(np.sum(all_cell_mass_kg, dtype=np.float64))
        mesh_volume_mass_kg = float(mesh_volume_mm3 * 1e-9 * float(density_kg_m3))
        cad_volume_mass_kg = float(cad_volume_mm3 * 1e-9 * float(density_kg_m3))
        mass_balance_relative_error = abs(total_mass_kg - mesh_volume_mass_kg) / max(mesh_volume_mass_kg, 1e-30)
        mass_delta_vs_cad_relative = (total_mass_kg - cad_volume_mass_kg) / cad_volume_mass_kg

        finite_ok = mesh["finite_cell_count"] == cell_count and np.isfinite(mesh["points_mm"]).all()
        positive_ok = mesh["positive_cell_volume_count"] == cell_count
        unique_ok = mesh["centroid_unique_count"] == cell_count
        centroid_cad_ok = centroid_state_counts == {"IN": cell_count, "OUT": 0, "UNKNOWN": 0}
        interior_cad_ok = interior_state_counts.get("IN", 0) == len(interior_states) and interior_state_counts.get("UNKNOWN", 0) == 0
        cavity_ok = (
            cavity["sample_state_counts"] == {"IN": 0, "OUT": cavity["sample_count"], "UNKNOWN": 0}
            and cavity["candidate_has_no_centroid_in_sampled_cavity_region"]
        )
        mass_ok = mass_balance_relative_error <= 1e-12
        volume_ok = volume_error <= PROPOSED_VOLUME_RELATIVE_ERROR_LIMIT
        identity_ok = (
            mesh["node_id_unique"]
            and mesh["cell_id_unique"]
            and mesh["connectivity_references_known_nodes"]
        )
        checks = {
            "finite_coordinates_and_volumes": bool(finite_ok),
            "all_cell_volumes_positive": bool(positive_ok),
            "unique_centroids_at_1e-10_mm_rounding": bool(unique_ok),
            "node_and_cell_ID_integrity": bool(identity_ok),
            "all_centroids_native_OCC_IN": bool(centroid_cad_ok),
            "stratified_tet_interior_samples_native_OCC_IN": bool(interior_cad_ok),
            "sampled_cavity_is_empty_and_has_no_centroid_points": bool(cavity_ok),
            "density_times_volume_mass_balance": bool(mass_ok),
            "proposed_CAD_volume_error_at_most_1_percent": bool(volume_ok),
        }
        passed = all(checks.values())
        mesh_file = output_dir / f"gmsh_mesh_h{mesh_size_mm:g}mm.msh"
        gmsh.write(str(mesh_file))
        trial.update(
            {
                "status": "PASS_CANDIDATE" if passed else "FAIL_CANDIDATE_ACCEPTANCE",
                "mesh_file": mesh_file.name,
                "total_tetra_volume_mm3": mesh_volume_mm3,
                "signed_volume_delta_mm3": float(mesh_volume_mm3 - cad_volume_mm3),
                "relative_volume_delta": relative_delta,
                "absolute_relative_volume_error": volume_error,
                "proposed_relative_volume_error_limit": PROPOSED_VOLUME_RELATIVE_ERROR_LIMIT,
                "checks": checks,
                "mesh_integrity": {
                    "finite_cell_count": mesh["finite_cell_count"],
                    "positive_cell_volume_count": mesh["positive_cell_volume_count"],
                    "zero_or_negative_unsigned_volume_count": mesh["zero_or_negative_unsigned_volume_count"],
                    "negative_signed_orientation_count": mesh["negative_signed_orientation_count"],
                    "centroid_unique_count": mesh["centroid_unique_count"],
                    "node_id_unique": mesh["node_id_unique"],
                    "cell_id_unique": mesh["cell_id_unique"],
                    "connectivity_references_known_nodes": mesh["connectivity_references_known_nodes"],
                    "boundary_face_count": mesh["boundary_face_count"],
                    "nonmanifold_tetra_face_count": mesh["nonmanifold_tetra_face_count"],
                },
                "native_cad_centroid_classification": {
                    "classifier": "OpenCascade BRepClass3d_SolidClassifier on canonical STEP-imported CadQuery/OCC solid",
                    "tolerance_mm": CLASSIFIER_TOLERANCE_MM,
                    "scope": "all linear-tetrahedron centroids",
                    "counts": centroid_state_counts,
                    "unknown_first_20_cell_ids": [
                        int(mesh["sourcecellIDs"][i])
                        for i, state in enumerate(centroid_states)
                        if state == "UNKNOWN"
                    ][:20],
                    "outside_first_20_cell_ids": [
                        int(mesh["sourcecellIDs"][i])
                        for i, state in enumerate(centroid_states)
                        if state == "OUT"
                    ][:20],
                },
                "native_cad_finite_interior_point_sampling": {
                    "classifier": "OpenCascade BRepClass3d_SolidClassifier",
                    "count": len(interior_states),
                    "scope": "one deterministic Dirichlet barycentric sample in each of up to 512 evenly source-ID-spaced tetrahedra; samples lie strictly inside the linear tetrahedra, not on a uniform whole-solid grid",
                    "seed": INTERIOR_SAMPLE_SEED,
                    "state_counts": interior_state_counts,
                    "unknown_count": int(interior_state_counts.get("UNKNOWN", 0)),
                    "outside_count": int(interior_state_counts.get("OUT", 0)),
                },
                "cavity_sample": cavity,
                "mass_quadrature": {
                    "density_kg_m3": float(density_kg_m3),
                    "density_status": RHO_STATUS,
                    "cell_mass_formula": "cell_volume_mm3 * 1e-9 m3/mm3 * rho_kg_m3",
                    "candidate_total_mass_kg": total_mass_kg,
                    "density_times_mesh_volume_kg": mesh_volume_mass_kg,
                    "mass_balance_relative_error_within_mesh_quadrature": float(mass_balance_relative_error),
                    "cad_reference_mass_kg": cad_volume_mass_kg,
                    "relative_mass_delta_vs_exact_CAD_solid": float(mass_delta_vs_cad_relative),
                    "interpretation": "Point masses conserve the candidate's tetrahedral volume sum by construction; density and mechanics remain unvalidated prototype placeholders.",
                },
                "cell_extent_stats_mm": {
                    "definition": "maximum of the six node-to-node distances in each linear tetrahedron",
                    "min": float(np.min(np.max(np.linalg.norm(tetra_xyz[:, :, None, :] - tetra_xyz[:, None, :, :], axis=-1), axis=(1, 2)))),
                    "median": float(np.median(np.max(np.linalg.norm(tetra_xyz[:, :, None, :] - tetra_xyz[:, None, :, :], axis=-1), axis=(1, 2)))),
                    "p95": float(np.percentile(np.max(np.linalg.norm(tetra_xyz[:, :, None, :] - tetra_xyz[:, None, :, :], axis=-1), axis=(1, 2)), 95)),
                    "max": float(np.max(np.max(np.linalg.norm(tetra_xyz[:, :, None, :] - tetra_xyz[:, None, :, :], axis=-1), axis=(1, 2)))),
                },
                "candidate_arrays": mesh,
                "centroid_states": centroid_states,
            }
        )
        return trial
    except Exception as error:
        trial.update({"status": "FAIL_MESH_TRIAL_EXCEPTION", "error": f"{type(error).__name__}: {error}"})
        return trial


def _save_candidate_npz(
    output_dir: Path,
    trial: dict[str, Any],
    step_hash: str,
    asset_dir: Path,
    manifest: dict[str, Any],
    gmsh_version: str,
) -> Path:
    mesh = trial["candidate_arrays"]
    candidate_path = output_dir / "cad_volume_candidate.npz"
    provenance = {
        "experiment_id": EXPERIMENT_ID,
        "source_step": str((asset_dir / manifest["geometry"]["solid_step"]).as_posix()),
        "source_step_sha256": step_hash,
        "gmsh_version": gmsh_version,
        "gmsh_kernel": "OpenCASCADE STEP import; 3D first-order linear tetrahedra",
        "mesh_size_target_mm": trial["mesh_size_target_mm"],
        "quadrature": "one cell centroid with positive tetrahedral volume weight per linear tetrahedron",
        "volume_error_relative": trial["relative_volume_delta"],
        "candidate_status": trial["status"],
        "budget_policy_id": trial["policy_id"],
        "max_volume_cells": trial["max_volume_cells"],
        "frame": manifest["frame"],
        "length_unit": manifest["units"]["length"],
        "cell_volume_unit": "mm3",
        "density_kg_m3": manifest["material"]["mechanics"]["rho_kg_m3"],
        "density_status": manifest["material"]["calibration_status"],
        "not_a_solver_job": True,
    }
    np.savez_compressed(
        candidate_path,
        points_mm=mesh["points_mm"].astype(np.float64, copy=False),
        cell_volume_mm3=mesh["cell_volume_mm3"].astype(np.float64, copy=False),
        sourcecellIDs=mesh["sourcecellIDs"].astype(np.int64, copy=False),
        source_node_ids=mesh["source_node_ids"].astype(np.int64, copy=False),
        tetra_node_ids=mesh["tetra_node_ids"].astype(np.int64, copy=False),
        source_step_sha256=np.asarray(step_hash),
        frame=np.asarray(str(manifest["frame"])),
        units_length=np.asarray(str(manifest["units"]["length"])),
        units_cell_volume=np.asarray("mm3"),
        provenance_json=np.asarray(json.dumps(provenance, sort_keys=True)),
    )
    return candidate_path


def _render_candidate_png(
    output_dir: Path, trial: dict[str, Any], bbox: dict[str, float]
) -> Path:
    mesh = trial["candidate_arrays"]
    points = mesh["points_mm"]
    node_ids = mesh["source_node_ids"]
    coordinates = mesh["node_coordinates_mm"]
    node_row_by_id = {int(node_id): row for row, node_id in enumerate(node_ids)}
    connectivity = mesh["tetra_node_ids"]
    boundary_faces = mesh["boundary_node_ids"]
    boundary_rows = np.asarray(
        [[node_row_by_id[int(node_id)] for node_id in face] for face in boundary_faces], dtype=np.int64
    )
    boundary_xyz = coordinates[boundary_rows]
    tetra_rows = np.asarray(
        [[node_row_by_id[int(node_id)] for node_id in tet] for tet in connectivity], dtype=np.int64
    )
    tetra_xyz = coordinates[tetra_rows]
    extents = np.max(
        np.linalg.norm(tetra_xyz[:, :, None, :] - tetra_xyz[:, None, :, :], axis=-1), axis=(1, 2)
    )
    fig = plt.figure(figsize=(13.2, 6.8), constrained_layout=True)
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    shell = Poly3DCollection(
        boundary_xyz,
        facecolor="#c9c7c1",
        edgecolor="#66645f",
        linewidth=0.12,
        alpha=0.56,
    )
    ax.add_collection3d(shell)
    # Plot a sparse, deterministic wireframe sample so cell extent is visible
    # independently of the much smaller centroid point symbols.
    edge_pairs = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
    sample_cell_ids = np.linspace(0, len(tetra_xyz) - 1, min(90, len(tetra_xyz)), dtype=np.int64)
    for cell_id in sample_cell_ids:
        tet = tetra_xyz[int(cell_id)]
        for first, second in edge_pairs:
            ax.plot(
                tet[[first, second], 0],
                tet[[first, second], 1],
                tet[[first, second], 2],
                color="#77746d",
                linewidth=0.28,
                alpha=0.22,
            )
    ax.scatter(points[:, 0], points[:, 1], points[:, 2], s=2.0, c="#242321", alpha=0.55, depthshade=False)
    ax.set_xlim(bbox["xmin_mm"] - 2, bbox["xmax_mm"] + 2)
    ax.set_ylim(bbox["ymin_mm"] - 2, bbox["ymax_mm"] + 2)
    ax.set_zlim(bbox["zmin_mm"] - 2, bbox["zmax_mm"] + 2)
    ax.set_box_aspect((bbox["xlen_mm"], bbox["ylen_mm"], bbox["zlen_mm"]))
    ax.view_init(elev=27, azim=-62)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_zlabel("z (mm)")
    ax.set_title("Canonical CAD boundary + sampled tetra extents")
    ax.text2D(
        0.02,
        0.97,
        "small dots = centroids; thin lines = 90 tetra cell extents",
        transform=ax.transAxes,
        va="top",
        fontsize=8,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 2.0},
    )

    ax2 = fig.add_subplot(1, 2, 2)
    section_polygons: list[np.ndarray] = []
    edge_pairs_for_section = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
    section_tolerance_mm = 1e-9
    for tet in tetra_xyz:
        intersections: list[np.ndarray] = []
        for vertex in tet:
            if abs(float(vertex[1])) <= section_tolerance_mm:
                intersections.append(vertex[[0, 2]])
        for i, j in edge_pairs_for_section:
            yi = float(tet[i, 1])
            yj = float(tet[j, 1])
            if yi * yj < 0.0:
                fraction = -yi / (yj - yi)
                point = tet[i] + fraction * (tet[j] - tet[i])
                intersections.append(point[[0, 2]])
        unique: list[np.ndarray] = []
        for point in intersections:
            if all(float(np.linalg.norm(point - seen)) > 1e-8 for seen in unique):
                unique.append(np.asarray(point, dtype=np.float64))
        if len(unique) < 3:
            continue
        polygon = np.asarray(unique, dtype=np.float64)
        center = np.mean(polygon, axis=0)
        angle = np.arctan2(polygon[:, 1] - center[1], polygon[:, 0] - center[0])
        section_polygons.append(polygon[np.argsort(angle)])
    section_collection = PolyCollection(
        section_polygons,
        facecolor="#d1cfca",
        edgecolor="#5d5b56",
        linewidth=0.3,
        alpha=0.88,
    )
    ax2.add_collection(section_collection)
    centroid_slab_halfwidth_mm = 3.0
    central_slab = np.abs(points[:, 1]) <= centroid_slab_halfwidth_mm
    ax2.scatter(
        points[central_slab, 0],
        points[central_slab, 2],
        s=5.0,
        c="#1f1e1b",
        alpha=0.8,
        linewidths=0,
        zorder=3,
    )
    ax2.set_xlim(bbox["xmin_mm"] - 1.0, bbox["xmax_mm"] + 1.0)
    ax2.set_ylim(bbox["zmin_mm"] - 1.0, bbox["zmax_mm"] + 1.0)
    ax2.set_aspect("equal", adjustable="box")
    ax2.set_xlabel("x (mm), y=0 section")
    ax2.set_ylabel("z (mm)")
    ax2.set_title("Actual y=0 tetra section + central-slab centroids")
    ax2.add_patch(
        Rectangle(
            (-8.0, 8.0),
            16.0,
            46.0,
            fill=False,
            edgecolor="#6e6b65",
            linewidth=0.7,
            linestyle="--",
        )
    )
    ax2.annotate(
        "empty cavity probe subset\n119 OCC samples: all OUT",
        xy=(0.0, 31.0),
        xytext=(8.0, 28.0),
        fontsize=7,
        color="#4d4b46",
        arrowprops={"arrowstyle": "->", "color": "#6e6b65", "lw": 0.6},
    )
    ax2.grid(True, linewidth=0.35, alpha=0.25)
    fig.suptitle(
        f"Neutral geometry-only candidate | h={trial['mesh_size_target_mm']:g} mm | "
        f"{mesh['cell_count']} tetra cells | median max edge {np.median(extents):.2f} mm",
        fontsize=12,
    )
    path = output_dir / "cad_volume_candidate.png"
    fig.savefig(path, dpi=180, facecolor="white")
    plt.close(fig)
    return path


def _write_report(output_dir: Path, report: dict[str, Any]) -> None:
    (output_dir / "report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8"
    )
    trial_rows = []
    for trial in report["mesh_trials"]:
        volume_error = trial.get("absolute_relative_volume_error")
        trial_rows.append(
            "| {h:g} | {status} | {cells} | {vol} | {cents} |".format(
                h=trial["mesh_size_target_mm"],
                status=trial["status"],
                cells=trial.get("volume_cell_count", "—"),
                vol=(f"{100.0 * volume_error:.5f}%" if volume_error is not None else "—"),
                cents=(trial.get("native_cad_centroid_classification", {}).get("counts", "—")),
            )
        )
    selected = report.get("selected_candidate")
    cad = report["canonical_cad"]
    cavity = selected.get("cavity_sample", {}) if selected else {}
    audit = report.get("bond_segment_audit", {})
    policy = report.get("run_policy", {})
    lines = [
        "# Bounded CAD-volume candidate spike",
        "",
        f"**Status:** `{report['status']}`. This is a geometry representation candidate only; it is not a PG solver job or a fracture/fragment result.",
        f"- Run policy: `{policy.get('policy_id', 'unknown')}`; cell ceiling {policy.get('max_volume_cells', 'unknown')}; sizes {policy.get('preselected_mesh_sizes_mm', 'unknown')} mm.",
        f"- Policy note: {policy.get('policy_note', 'unavailable')}",
        "",
        "## Canonical input",
        "",
        f"- STEP: `{report['source']['step_relative_path']}`",
        f"- SHA-256: `{report['source']['step_sha256']}` (matches manifest; {report['source']['step_bytes']} bytes)",
        f"- Coordinate frame / length unit: `{report['coordinate_convention']['frame']}` / `{report['coordinate_convention']['length_unit']}`",
        f"- B-rep: {cad['solid_count']} valid solid, {cad['shell_count']} shell; volume {cad['volume_mm3']:.6f} mm³; bbox x=[{cad['bbox_mm']['xmin_mm']:.3f}, {cad['bbox_mm']['xmax_mm']:.3f}], y=[{cad['bbox_mm']['ymin_mm']:.3f}, {cad['bbox_mm']['ymax_mm']:.3f}], z=[{cad['bbox_mm']['zmin_mm']:.3f}, {cad['bbox_mm']['zmax_mm']:.3f}] mm",
        "",
        "## Mesh trials",
        "",
        "| Target size (mm) | Status | Volume cells | Absolute volume error | Native OCC centroid states |",
        "|---:|---|---:|---:|---|",
        *trial_rows,
        "",
        f"The proposed CAD volume error gate is {100.0 * PROPOSED_VOLUME_RELATIVE_ERROR_LIMIT:.2f}%; it was not adjusted to force acceptance.",
        "",
    ]
    if selected:
        lines.extend(
            [
                "## Saved quadrature candidate",
                "",
                f"- NPZ: `{selected['npz_file']}`",
                f"- PNG: `{selected['png_file']}`",
                f"- Selected target size / cells: {selected['mesh_size_target_mm']:g} mm / {selected['volume_cell_count']}",
                f"- Tetra volume sum: {selected['total_tetra_volume_mm3']:.6f} mm³; signed relative error {selected['relative_volume_delta']:+.6%}",
                f"- Cell-volume quadrature mass: {selected['mass_quadrature']['candidate_total_mass_kg']:.9f} kg at rho={selected['mass_quadrature']['density_kg_m3']:g} kg/m³ (`{selected['mass_quadrature']['density_status']}`)",
                f"- Point-mass conservation error against rho×tetra-volume sum: {selected['mass_quadrature']['mass_balance_relative_error_within_mesh_quadrature']:.3e}",
                f"- Centroid OCC states: `{selected['native_cad_centroid_classification']['counts']}`; finite interior sample states: `{selected['native_cad_finite_interior_point_sampling']['state_counts']}`",
                f"- Empty cavity subset: `{cavity.get('sample_state_counts')}` for {cavity.get('sample_count')} finite OCC-classified probes; candidate centroids in that subset: {cavity.get('mesh_centroid_count_in_sampled_cavity_region')}",
                "",
            ]
        )
    lines.extend(
        [
            "## Segment obstruction probe",
            "",
            f"Status: `{audit.get('status', 'NOT_RUN')}`; cross-wall examples: {len(audit.get('crosswall_candidates', []))}; same-wall controls: {len(audit.get('samewall_controls', []))}.",
            "Each tested segment requires both endpoints to classify OCC IN. OCC face intersections partition the finite line, then native solid classification labels interval midpoints. This remains a small geometry probe with no horizon and no native PG bond graph.",
            "",
            "## Scope and limitations",
            "",
            "- Linear tetrahedral boundary facets approximate curved CAD faces. The cell-centroid checks and finite samples are evidence for this mesh, not a global boundary proof.",
            "- The centroid set is an irregular discretization; directional bias, free-surface correction, and horizon selection remain pending.",
            "- The reported 4 mm wall is the upstream analytic/profile check. This experiment does not prove a global minimum wall thickness over the fused solid.",
            "- rho=2100 kg/m³ and the other mechanics values are unvalidated placeholders. Mass conservation here is only density times the candidate volume weights.",
            "- No crack surfaces, R2 fragment geometry, supported PG asset/job schema, or solver output are created.",
            "",
            "## Reproduction",
            "",
            "```powershell",
            report["reproduction"]["powershell_command"],
            "```",
            "",
            "The command uses the existing VesselForge Python environment read-only plus the isolated official Gmsh 4.15.2 SDK under `third_party/pgfrag-geometry-runtime`; it installs no packages and writes only this spike's output directory.",
            "Official references: [Gmsh documentation](https://gmsh.info/doc/texinfo/) (OCC `importShapes`, STEP tutorial t20, 3D tetrahedral mesh API) and [CadQuery STEP import documentation](https://cadquery.readthedocs.io/en/latest/importexport.html).",
            "",
        ]
    )
    (output_dir / "report.md").write_text("\n".join(lines), encoding="utf-8")


def run(
    asset_dir: Path,
    output_dir: Path,
    *,
    max_volume_cells: int = MAX_VOLUME_CELLS,
    mesh_sizes_mm: tuple[float, ...] | None = None,
    policy_id: str = "bounded-preselected-2500cells",
    policy_note: str = "Original bounded spike policy: at most two preselected coarse sizes, hard ceiling 2500 volume tetrahedra.",
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    mesh_sizes = tuple(float(size) for size in (mesh_sizes_mm or MESH_SIZE_TRIALS_MM))
    if not 1 <= len(mesh_sizes) <= 2:
        raise ValueError("the spike accepts one or two preselected mesh sizes")
    if any(not math.isfinite(size) or size <= 0.0 for size in mesh_sizes):
        raise ValueError("mesh sizes must be positive finite millimeter values")
    if int(max_volume_cells) <= 0:
        raise ValueError("volume-cell budget must be a positive integer")
    max_volume_cells = int(max_volume_cells)
    manifest, step_path, step_hash = _manifest_source(asset_dir)
    # The canonical validation artifact sits at the manifest-selected validation
    # path. Load it as independent upstream evidence; it remains read-only.
    validation_path = asset_dir / manifest["validation"]["report"]
    validation = json.loads(validation_path.read_text(encoding="utf-8"))
    if manifest.get("frame") != "right_handed_z_up" or manifest.get("units", {}).get("length") != "mm":
        raise ValueError("canonical asset is not in the expected right-handed Z-up millimeter convention")
    if validation.get("status") != "PASS":
        raise ValueError("canonical source geometry report is not PASS")
    if int(validation["measurements"]["solid_count"]) != 1 or int(validation["measurements"]["shell_count"]) != 1:
        raise ValueError("canonical source geometry does not report one fused solid and one shell")

    expected_rho_from_manifest = float(manifest["material"]["mechanics"]["rho_kg_m3"])
    imported = cq.importers.importStep(str(step_path))
    cad_shape = imported.val()
    solids = cad_shape.Solids()
    shells = cad_shape.Shells()
    if len(solids) != 1 or len(shells) != 1 or not cad_shape.isValid():
        raise ValueError(
            f"CadQuery/OCC canonical STEP check failed: solids={len(solids)} shells={len(shells)} valid={cad_shape.isValid()}"
        )
    cad_solid = solids[0]
    cad_volume_mm3 = float(cad_solid.Volume())
    cad_bbox_mm = _bbox_mm(cad_solid)
    expected_cad_volume_mm3 = float(validation["measurements"]["solid_volume_mm3"])
    if abs(cad_volume_mm3 - expected_cad_volume_mm3) / expected_cad_volume_mm3 > 1e-8:
        raise ValueError("re-imported STEP solid volume disagrees with the manifest validation report")

    gmsh = _load_gmsh()
    gmsh.initialize()
    gmsh_version = str(gmsh.__version__)
    gmsh_sdk_sha256 = _sha256(GMSH_SDK_ARCHIVE) if GMSH_SDK_ARCHIVE.is_file() else None
    trials: list[dict[str, Any]] = []
    try:
        for index, mesh_size in enumerate(mesh_sizes):
            if index == 1 and trials and trials[0].get("volume_cell_count") is not None:
                first_count = int(trials[0]["volume_cell_count"])
                predicted_second = math.ceil(
                    first_count
                    * (mesh_sizes[0] / mesh_size) ** 3
                    * 1.20
                )
                if predicted_second > max_volume_cells:
                    trials.append(
                        {
                            "mesh_size_target_mm": mesh_size,
                            "max_volume_cells": max_volume_cells,
                            "policy_id": policy_id,
                            "status": "SKIPPED_PREDICTED_CELL_BUDGET",
                            "structured_failure": {
                                "code": "predicted_finer_trial_exceeds_cell_budget",
                                "prediction_basis": "first-trial cell count scaled by target-size ratio cubed and a 20% packing margin; guard only, not a measured cell count",
                                "predicted_volume_cells": predicted_second,
                                "limit": max_volume_cells,
                            },
                        }
                    )
                    break
            trial = _run_mesh_trial(
                gmsh,
                step_path,
                output_dir,
                mesh_size,
                index,
                cad_solid,
                cad_volume_mm3,
                expected_rho_from_manifest,
                max_volume_cells,
                policy_id,
            )
            trials.append(trial)
            if trial.get("volume_cell_count", 0) > max_volume_cells:
                break
    finally:
        gmsh.finalize()

    valid_array_trials = [trial for trial in trials if trial.get("candidate_arrays") is not None]
    accepted_trials = [trial for trial in valid_array_trials if trial["status"] == "PASS_CANDIDATE"]
    if accepted_trials:
        selected = min(accepted_trials, key=lambda trial: trial["absolute_relative_volume_error"])
        overall_status = "PASS_CANDIDATE"
    elif valid_array_trials:
        selected = min(
            valid_array_trials,
            key=lambda trial: trial.get("absolute_relative_volume_error", float("inf")),
        )
        overall_status = "FAIL_CANDIDATE_ACCEPTANCE"
    else:
        selected = None
        overall_status = "FAIL_NO_BOUNDED_VOLUME_CANDIDATE"

    candidate_path: Path | None = None
    png_path: Path | None = None
    bond_audit: dict[str, Any] = {"status": "NOT_RUN_NO_MESH"}
    if selected is not None:
        candidate_path = _save_candidate_npz(
            output_dir, selected, step_hash, asset_dir, manifest, gmsh_version
        )
        png_path = _render_candidate_png(output_dir, selected, cad_bbox_mm)
        bond_audit = _bond_segment_audit(
            cad_solid, selected["candidate_arrays"], selected["centroid_states"]
        )

    source_relative = step_path.relative_to(ROOT).as_posix()
    selected_summary = None
    if selected is not None:
        keys_to_keep = [
            "mesh_size_target_mm",
            "status",
            "volume_cell_count",
            "source_node_count",
            "mesh_file",
            "total_tetra_volume_mm3",
            "signed_volume_delta_mm3",
            "relative_volume_delta",
            "absolute_relative_volume_error",
            "proposed_relative_volume_error_limit",
            "checks",
            "mesh_integrity",
            "native_cad_centroid_classification",
            "native_cad_finite_interior_point_sampling",
            "cavity_sample",
            "mass_quadrature",
            "cell_extent_stats_mm",
        ]
        selected_summary = {key: selected[key] for key in keys_to_keep if key in selected}
        selected_summary["npz_file"] = candidate_path.name if candidate_path else None
        selected_summary["png_file"] = png_path.name if png_path else None

    report = {
        "experiment_id": EXPERIMENT_ID,
        "status": overall_status,
        "decision": "A native CAD-to-volume point-quadrature representation was exercised on the existing canonical STEP; only PASS_CANDIDATE meets every declared check. This is not an assertion that PG supports this input or that it produces realistic fracture fragments.",
        "source": {
            "manifest_relative_path": (asset_dir / "manifest.json").relative_to(ROOT).as_posix(),
            "step_relative_path": source_relative,
            "step_sha256": step_hash,
            "step_bytes": step_path.stat().st_size,
            "manifest_hash_check": "PASS",
            "asset_producer": manifest["provenance"]["producer"],
            "asset_backend": manifest["provenance"]["backend"],
            "asset_backend_version": manifest["provenance"]["backend_version"],
        },
        "coordinate_convention": {
            "frame": manifest["frame"],
            "length_unit": manifest["units"]["length"],
            "quadrature_point": "linear tetrahedron centroid in canonical asset-local coordinates",
            "volume_weight_unit": "mm3",
        },
        "canonical_cad": {
            "source_validation_status": validation["status"],
            "solid_count": len(solids),
            "shell_count": len(shells),
            "valid": bool(cad_solid.isValid()),
            "volume_mm3": cad_volume_mm3,
            "validation_report_volume_mm3": expected_cad_volume_mm3,
            "bbox_mm": cad_bbox_mm,
            "nominal_prototype_parameters": manifest["provenance"]["spec"]["parameters"],
        },
        "mesher": {
            "name": "Gmsh",
            "version": gmsh_version,
            "official_sdk_archive": GMSH_SDK_ARCHIVE.relative_to(ROOT).as_posix(),
            "official_sdk_archive_sha256": gmsh_sdk_sha256,
            "cad_kernel": "OpenCASCADE",
            "step_import_method": "gmsh.model.occ.importShapes(step, highestDimOnly=True, format='step')",
            "mesh_type": "first-order linear tetrahedron",
            "quadrature_rule": "one cell centroid weighted by positive tetrahedral cell volume",
            "mesh_size_trials_mm_preselected": list(mesh_sizes),
            "hard_volume_cell_limit": max_volume_cells,
            "mesh_size_from_curvature": "disabled; a separate OCC-vs-tetra volume and native solid-classifier audit measures resulting approximation",
        },
        "acceptance_policy": {
            "proposed_cad_volume_relative_error_limit": PROPOSED_VOLUME_RELATIVE_ERROR_LIMIT,
            "threshold_adjusted_after_results": False,
            "pass_requires": [
                f"at most {max_volume_cells} linear tetrahedron cells",
                "finite node coordinates, centroid coordinates and cell volumes",
                "strictly positive cell volumes",
                "unique centroid coordinates at 1e-10 mm rounding",
                "unique cell/node IDs and valid node references",
                "every centroid and every finite stratified interior sample classifies OCC IN",
                "all bounded central cavity probes classify OCC OUT and no point centroid is in that sampled cavity subset",
                "rho times cell-volume quadrature balances its own tetra volume sum",
                "absolute tetra total volume error against native CAD <= 1%",
            ],
        },
        "mesh_trials": [
            {key: value for key, value in trial.items() if key not in ("candidate_arrays", "centroid_states")}
            for trial in trials
        ],
        "selected_candidate": selected_summary,
        "run_policy": {
            "policy_id": policy_id,
            "max_volume_cells": max_volume_cells,
            "preselected_mesh_sizes_mm": list(mesh_sizes),
            "policy_note": policy_note,
            "original_2500_cell_h8_report_preserved": policy_id != "bounded-preselected-2500cells",
            "prior_attempt_report": "../report.json" if policy_id != "bounded-preselected-2500cells" else None,
        },
        "bond_segment_audit": bond_audit,
        "material": {
            "rho_kg_m3": expected_rho_from_manifest,
            "status": manifest["material"]["calibration_status"],
            "mechanics_provenance": manifest["material"]["provenance"]["mechanics_source"],
        },
        "limitations": [
            "Linear tetrahedral boundary facets approximate curved CAD faces; native centroid and finite interior samples are finite evidence, not a global boundary proof.",
            "The point set is an irregular centroid discretization; directional bias, free-surface correction, and horizon selection are pending.",
            "The mesh does not establish a global minimum wall thickness over the fused CAD solid.",
            "The segment obstruction probe is not tied to a horizon and does not create or validate a Peridynamics bond graph.",
            "Density and mechanics are unvalidated prototype placeholders; rho-times-volume mass conservation is only candidate quadrature bookkeeping.",
            "No solver job, crack surface, realistic R2 fragment geometry, or actual PG fragment result was produced.",
        ],
        "reproduction": {
            "powershell_command": (
                "& 'third_party/vesselforge-runtime/venv/Scripts/python.exe' -u "
                "'pgfrag/experiments/cad_volume_spike.py' "
                "--asset-dir 'relicweave-synlab/tmp/review_20261003/implementation/vf/artifacts/four_bands_brep_semantics1_final' "
                f"--output-dir '{output_dir.relative_to(ROOT).as_posix()}' "
                f"--mesh-size-mm {' '.join(f'{size:g}' for size in mesh_sizes)} "
                f"--max-volume-cells {max_volume_cells} --policy-id '{policy_id}' "
                f"--policy-note '{policy_note}' "
                f"2>&1 | Tee-Object -FilePath '{(output_dir / 'run.log').relative_to(ROOT).as_posix()}'"
            ),
            "runtime": "Existing third_party/vesselforge-runtime/venv (read-only) with official Gmsh 4.15.2 Windows64 SDK from third_party/pgfrag-geometry-runtime/sdk; no package installation performed by the spike.",
            "mesh_sizes_mm": list(mesh_sizes),
            "run_date_local": time.strftime("%Y-%m-%d %H:%M:%S %Z"),
        },
        "official_references": [
            {
                "title": "Gmsh documentation: OCC STEP import, t20, and tetrahedral mesh API",
                "url": "https://gmsh.info/doc/texinfo/",
            },
            {
                "title": "CadQuery documentation: Importing STEP files",
                "url": "https://cadquery.readthedocs.io/en/latest/importexport.html",
            },
        ],
    }
    _write_report(output_dir, report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--asset-dir", type=Path, default=DEFAULT_ASSET_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--mesh-size-mm",
        type=float,
        nargs="+",
        default=None,
        help="one or two explicitly preselected mesh sizes; defaults to the original (8, 6) mm pair",
    )
    parser.add_argument("--max-volume-cells", type=int, default=MAX_VOLUME_CELLS)
    parser.add_argument("--policy-id", default="bounded-preselected-2500cells")
    parser.add_argument(
        "--policy-note",
        default="Original bounded spike policy: at most two preselected coarse sizes, hard ceiling 2500 volume tetrahedra.",
    )
    args = parser.parse_args()
    try:
        report = run(
            args.asset_dir.resolve(),
            args.output_dir.resolve(),
            max_volume_cells=args.max_volume_cells,
            mesh_sizes_mm=tuple(args.mesh_size_mm) if args.mesh_size_mm else None,
            policy_id=args.policy_id,
            policy_note=args.policy_note,
        )
    except Exception as error:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        failure = {
            "experiment_id": EXPERIMENT_ID,
            "status": "FAIL_SETUP_OR_CANONICAL_INPUT",
            "error": f"{type(error).__name__}: {error}",
            "hard_volume_cell_limit": args.max_volume_cells,
            "mesh_size_trials_mm_preselected": args.mesh_size_mm or list(MESH_SIZE_TRIALS_MM),
            "policy_id": args.policy_id,
            "threshold_adjusted_after_results": False,
        }
        (args.output_dir.resolve() / "report.json").write_text(
            json.dumps(failure, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        (args.output_dir.resolve() / "report.md").write_text(
            "# Bounded CAD-volume candidate spike\n\n"
            "**Status:** `FAIL_SETUP_OR_CANONICAL_INPUT`\n\n"
            f"Error: `{failure['error']}`\n",
            encoding="utf-8",
        )
        print(json.dumps(failure, indent=2))
        return 1
    print(
        json.dumps(
            {
                "status": report["status"],
                "source_step_sha256": report["source"]["step_sha256"],
                "mesh_trials": [
                    {
                        "h_mm": trial["mesh_size_target_mm"],
                        "status": trial["status"],
                        "cells": trial.get("volume_cell_count"),
                        "abs_volume_error": trial.get("absolute_relative_volume_error"),
                    }
                    for trial in report["mesh_trials"]
                ],
                "selected_candidate": report["selected_candidate"]["npz_file"]
                if report["selected_candidate"]
                else None,
                "report_json": str((args.output_dir.resolve() / "report.json")),
            },
            indent=2,
        )
    )
    return 0 if report["status"] == "PASS_CANDIDATE" else 1


if __name__ == "__main__":
    raise SystemExit(main())
