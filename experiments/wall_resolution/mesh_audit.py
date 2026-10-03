"""Bounded quality and directional audit for existing whole-STEP tetra meshes.

This read-only audit does not create solver input or modify the source spike.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import cadquery as cq
import numpy as np
from OCP.IntCurvesFace import IntCurvesFace_ShapeIntersector
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.TopAbs import TopAbs_IN, TopAbs_ON, TopAbs_OUT
from OCP.gp import gp_Dir, gp_Lin, gp_Pnt
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[3]
SDK = ROOT / "third_party/pgfrag-geometry-runtime/sdk/gmsh-4.15.2-Windows64-sdk/lib"
ASSET = ROOT / "relicweave-synlab/tmp/review_20261003/implementation/vf/artifacts/four_bands_brep_semantics1_final"
EXPECTED_SHA = "4623bd373aeed0acfacfbf6cc0dda65952339561779b43429c59d61b8b5a1c30"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def qstats(x: np.ndarray) -> dict[str, float | int | None]:
    x = np.asarray(x, dtype=np.float64)
    x = x[np.isfinite(x)]
    if not len(x):
        return {"count": 0, "min": None, "p01": None, "p05": None, "median": None, "p95": None, "p99": None, "max": None}
    q = np.quantile(x, [0.01, 0.05, 0.5, 0.95, 0.99])
    return {"count": int(len(x)), "min": float(x.min()), "p01": float(q[0]), "p05": float(q[1]), "median": float(q[2]), "p95": float(q[3]), "p99": float(q[4]), "max": float(x.max())}


def mesh_from_msh(path: Path) -> dict[str, Any]:
    if not SDK.is_dir():
        raise FileNotFoundError(f"Gmsh SDK missing: {SDK}")
    if hasattr(os, "add_dll_directory"):
        os.add_dll_directory(str(SDK))
    sys.path.insert(0, str(SDK))
    import gmsh  # type: ignore

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try:
        gmsh.open(str(path))
        types, tags_by_type, nodes_by_type = gmsh.model.mesh.getElements(3)
        tetra_type = int(gmsh.model.mesh.getElementType("tetrahedron", 1))
        if not types or any(int(t) != tetra_type for t in types):
            raise ValueError(f"expected only first-order tets in {path}; got {list(map(int, types))}")
        tags = np.concatenate([np.asarray(x, dtype=np.int64) for x in tags_by_type])
        conn_ids = np.concatenate([np.asarray(x, dtype=np.int64).reshape(-1, 4) for x in nodes_by_type])
        node_ids, coords, _ = gmsh.model.mesh.getNodes()
        node_ids = np.asarray(node_ids, dtype=np.int64)
        xyz = np.asarray(coords, dtype=np.float64).reshape(-1, 3)
        row = {int(k): i for i, k in enumerate(node_ids)}
        conn = np.asarray([[row[int(k)] for k in tet] for tet in conn_ids], dtype=np.int64)
        tet_xyz = xyz[conn]
        p = tet_xyz
        det = np.einsum("ij,ij->i", np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), p[:, 3] - p[:, 0])
        signed_v = det / 6.0
        volume = np.abs(signed_v)
        cent = p.mean(axis=1)
        edge_ij = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
        edges = np.stack([np.linalg.norm(p[:, i] - p[:, j], axis=1) for i, j in edge_ij], axis=1)
        max_edge = edges.max(axis=1)
        min_edge = edges.min(axis=1)
        # Inradius and circumradius form a scale-free element shape indicator;
        # R/r is 3 for a regular tetrahedron and grows without bound for slivers.
        face_area = np.zeros(len(p), dtype=np.float64)
        face_normals = []
        for omit in range(4):
            ids = [i for i in range(4) if i != omit]
            u = p[:, ids[1]] - p[:, ids[0]]
            v = p[:, ids[2]] - p[:, ids[0]]
            n = np.cross(u, v)
            area2 = np.linalg.norm(n, axis=1)
            face_area += 0.5 * area2
            # Orient the normal away from the opposite tetra vertex.
            inside = np.einsum("ij,ij->i", n, p[:, omit] - p[:, ids[0]]) > 0
            n[inside] *= -1
            norm = np.maximum(np.linalg.norm(n, axis=1), 1e-300)
            face_normals.append(n / norm[:, None])
        inradius = 3.0 * volume / np.maximum(face_area, 1e-300)
        matrices = 2.0 * (p[:, 1:] - p[:, :1])
        rhs = np.sum((p[:, 1:] - p[:, :1]) ** 2, axis=2)
        circum = np.full(len(p), np.nan)
        good = np.isfinite(matrices).all(axis=(1, 2)) & (np.abs(signed_v) > 0)
        circum_center = np.linalg.solve(matrices[good], rhs[good][..., None])[..., 0]
        circum[good] = np.linalg.norm(circum_center, axis=1)
        # Matrix above solves coordinates relative to p0.
        radius_ratio = circum / np.maximum(inradius, 1e-300)
        radius_edge = circum / np.maximum(min_edge, 1e-300)
        dihedrals = []
        for i in range(4):
            for j in range(i + 1, 4):
                cosine = np.einsum("ij,ij->i", face_normals[i], face_normals[j]).clip(-1, 1)
                dihedrals.append(np.degrees(np.arccos(-cosine)))
        min_dihedral = np.min(np.stack(dihedrals, axis=1), axis=1)
        # Face adjacency and cell-to-cell centroid spacing are independent of
        # the global target size and expose abrupt local size changes.
        faces: dict[tuple[int, int, int], list[int]] = defaultdict(list)
        for ci, t in enumerate(conn_ids):
            for tri in ((t[0], t[1], t[2]), (t[0], t[3], t[1]), (t[0], t[2], t[3]), (t[1], t[3], t[2])):
                faces[tuple(sorted(map(int, tri)))].append(ci)
        pairs = [(cs[0], cs[1]) for cs in faces.values() if len(cs) == 2]
        pair = np.asarray(pairs, dtype=np.int64).reshape(-1, 2)
        neighbor_spacing = np.linalg.norm(cent[pair[:, 0]] - cent[pair[:, 1]], axis=1) if len(pair) else np.array([])
        neighbor_vol_ratio = np.maximum(volume[pair[:, 0]], volume[pair[:, 1]]) / np.maximum(np.minimum(volume[pair[:, 0]], volume[pair[:, 1]]), 1e-300) if len(pair) else np.array([])
        return {
            "mesh_path": str(path), "mesh_sha256": sha256(path), "node_ids": node_ids,
            "xyz": xyz, "connectivity_ids": conn_ids, "connectivity": conn,
            "tet_xyz": tet_xyz, "cell_ids": tags, "centroids": cent,
            "signed_volume": signed_v, "volume": volume, "edges": edges,
            "max_edge": max_edge, "min_edge": min_edge, "radius_ratio": radius_ratio,
            "radius_edge": radius_edge, "min_dihedral": min_dihedral,
            "neighbor_spacing": neighbor_spacing, "neighbor_volume_ratio": neighbor_vol_ratio,
        }
    finally:
        gmsh.finalize()


def radial_hits(solid: cq.Shape, z: float, theta: float) -> list[float]:
    direction = (math.cos(theta), math.sin(theta), 0.0)
    line = gp_Lin(gp_Pnt(0.0, 0.0, float(z)), gp_Dir(*map(float, direction)))
    intersector = IntCurvesFace_ShapeIntersector()
    intersector.Load(solid.wrapped, 1e-7)
    intersector.Perform(line, 0.0, 40.0)
    if not intersector.IsDone():
        return []
    values = sorted(float(intersector.WParameter(i)) for i in range(1, int(intersector.NbPnt()) + 1))
    unique: list[float] = []
    for x in values:
        if not unique or abs(x - unique[-1]) > 1e-6:
            unique.append(x)
    return unique


def wall_coordinate(solid: cq.Shape, point: np.ndarray) -> float | None:
    theta = math.atan2(float(point[1]), float(point[0]))
    hits = radial_hits(solid, float(point[2]), theta)
    r = float(np.hypot(point[0], point[1]))
    if len(hits) < 2:
        return None
    # Body-only negative-x sector has ordered radial hits [inner, outer]; for
    # general rays with a handle, keep the material interval containing r.
    for a, b in zip(hits[::2], hits[1::2]):
        if a - 1e-3 <= r <= b + 1e-3:
            return (b - r) / (b - a) if b > a else None
    return None


def tensor_from_vectors(vec: np.ndarray, weights: np.ndarray) -> np.ndarray | None:
    lengths = np.linalg.norm(vec, axis=1)
    ok = np.isfinite(vec).all(axis=1) & np.isfinite(lengths) & (lengths > 1e-12) & (weights > 0)
    if not np.any(ok):
        return None
    u = vec[ok] / lengths[ok, None]
    w = weights[ok] / weights[ok].sum()
    return np.einsum("i,ij,ik->jk", w, u, u)


def classify_points(solid: cq.Shape, points: np.ndarray) -> list[str]:
    checker = BRepClass3d_SolidClassifier()
    checker.Load(solid.wrapped)
    labels: list[str] = []
    for point in np.asarray(points, dtype=np.float64):
        try:
            checker.Perform(gp_Pnt(*map(float, point)), 1e-7)
            state = checker.State()
        except Exception:
            state = None
        labels.append("IN" if state == TopAbs_IN else "OUT" if state == TopAbs_OUT else "UNKNOWN" if state == TopAbs_ON else "UNKNOWN")
    return labels


def geometry_clipped_fabric_reference(
    mesh: dict[str, Any],
    solid: cq.Shape,
    center_ids: np.ndarray,
    candidate_deltas_mm: tuple[float, ...],
    *,
    monte_carlo_samples: int = 5000,
    seed: int = 20261003,
) -> dict[str, Any]:
    """Compare centroid-neighbor fabric with the exact native-solid clipped ball.

    Uniform-in-volume Monte Carlo directions are classified against the native
    source CAD solid. This retains true wall/free-surface truncation as the
    physical reference instead of demanding isotropy at boundaries.
    """
    points = mesh["centroids"]
    cell_volume = mesh["volume"]
    labels = classify_points(solid, points)
    in_ids = np.asarray([i for i, state in enumerate(labels) if state == "IN"], dtype=np.int64)
    inside_points = points[in_ids]
    inside_volume = cell_volume[in_ids]
    tree = cKDTree(inside_points)
    checker = BRepClass3d_SolidClassifier()
    checker.Load(solid.wrapped)
    center_count_before_native_classification = int(len(center_ids))
    center_ids = np.asarray([i for i in center_ids if labels[int(i)] == "IN"], dtype=np.int64)
    summaries: dict[str, Any] = {}
    for delta_index, delta in enumerate(candidate_deltas_mm):
        records: list[dict[str, Any]] = []
        rng = np.random.default_rng(seed + delta_index)
        for ordinal, center_id in enumerate(center_ids):
            center = points[int(center_id)]
            candidate_local_ids = tree.query_ball_point(center, float(delta))
            neighbor_ids = in_ids[np.asarray(candidate_local_ids, dtype=np.int64)]
            neighbor_ids = neighbor_ids[neighbor_ids != int(center_id)]
            discrete = tensor_from_vectors(points[neighbor_ids] - center, cell_volume[neighbor_ids]) if len(neighbor_ids) else None
            if discrete is None:
                continue
            gaussian = rng.normal(size=(monte_carlo_samples, 3))
            gaussian /= np.maximum(np.linalg.norm(gaussian, axis=1, keepdims=True), 1e-300)
            radius = float(delta) * rng.random(monte_carlo_samples) ** (1.0 / 3.0)
            sample_points = center + gaussian * radius[:, None]
            in_mask = np.zeros(monte_carlo_samples, dtype=bool)
            for sample_index, sample in enumerate(sample_points):
                try:
                    checker.Perform(gp_Pnt(*map(float, sample)), 1e-7)
                    in_mask[sample_index] = checker.State() == TopAbs_IN
                except Exception:
                    in_mask[sample_index] = False
            native = tensor_from_vectors(gaussian[in_mask], np.ones(int(np.count_nonzero(in_mask)), dtype=np.float64)) if np.any(in_mask) else None
            if native is None:
                continue
            native_eig = np.linalg.eigvalsh(native).clip(0, 1)
            discrete_eig = np.linalg.eigvalsh(discrete).clip(0, 1)
            records.append({
                "center_id": int(center_id),
                "center_xyz_mm": center.tolist(),
                "neighbor_count": int(len(neighbor_ids)),
                "native_clipped_mc_inside_count": int(np.count_nonzero(in_mask)),
                "native_clipped_fabric_eigenvalues": native_eig.tolist(),
                "mesh_fabric_eigenvalues": discrete_eig.tolist(),
                "native_clipped_eigenvalue_max_over_min": float(native_eig[-1] / max(native_eig[0], 1e-12)),
                "mesh_eigenvalue_max_over_min": float(discrete_eig[-1] / max(discrete_eig[0], 1e-12)),
                "fabric_frobenius_error_vs_native_clipped_reference": float(np.linalg.norm(discrete - native, ord="fro")),
            })
        summaries[f"delta_{delta:g}mm_candidate_only"] = {
            "definition": "mesh tensor is cell-volume-weighted centroid-neighbor unit-direction second moment; physical reference is uniform-in-volume Monte Carlo inside the same-radius ball classified IN against the exact source CAD solid",
            "candidate_delta_is_not_a_selected_or_validated_horizon": True,
            "center_count": int(len(center_ids)),
            "centers_with_a_discrete_and_native_tensor": int(len(records)),
            "monte_carlo_samples_per_center": int(monte_carlo_samples),
            "per_center_neighbor_count": qstats(np.asarray([r["neighbor_count"] for r in records], dtype=np.float64)),
            "per_center_native_physical_fabric_anisotropy_max_over_min": qstats(np.asarray([r["native_clipped_eigenvalue_max_over_min"] for r in records], dtype=np.float64)),
            "per_center_discrete_mesh_fabric_anisotropy_max_over_min": qstats(np.asarray([r["mesh_eigenvalue_max_over_min"] for r in records], dtype=np.float64)),
            "per_center_fabric_frobenius_error_vs_native": qstats(np.asarray([r["fabric_frobenius_error_vs_native_clipped_reference"] for r in records], dtype=np.float64)),
            "records": records,
        }
    return {
        "wall_center_selection_count_before_native_classification": center_count_before_native_classification,
        "wall_center_native_IN_count": int(len(center_ids)),
        "native_centroid_states": {name: int(labels.count(name)) for name in ("IN", "OUT", "UNKNOWN")},
        "native_source_inside_centroid_count": int(len(in_ids)),
        "candidate_deltas": summaries,
    }


def summarize_mesh(
    mesh: dict[str, Any],
    size_mm: float,
    solid: cq.Shape,
    *,
    candidate_deltas_mm: tuple[float, ...] = (4.0, 6.0),
    wall_center_mask: np.ndarray | None = None,
    monte_carlo_samples: int = 5000,
    report_wall_distribution: bool = True,
) -> dict[str, Any]:
    vol = mesh["volume"]
    normalized_wall = []
    # Straight-ish opposite-handle band, away from the rim/base and profile
    # knots. This audits actual native radial intersections; box coordinates
    # are not treated as wall thickness.
    p = mesh["centroids"]
    r = np.hypot(p[:, 0], p[:, 1])
    theta = np.arctan2(p[:, 1], p[:, 0])
    roi = ((p[:, 0] < -14.5) & (np.abs(p[:, 1]) <= 3.0) & (p[:, 2] >= 22.0) & (p[:, 2] <= 27.0) & (np.abs(np.abs(theta) - math.pi) <= math.radians(12)))
    if report_wall_distribution:
        for point in p[roi]:
            val = wall_coordinate(solid, point)
            if val is not None:
                normalized_wall.append(val)

    # Neighbor-direction fabric at explicitly predeclared diagnostic radii.
    # These deltas are exploratory scales, never a selected material horizon.
    labels = classify_points(solid, p)
    in_indices = np.asarray([i for i, state in enumerate(labels) if state == "IN"], dtype=np.int64)
    inside_points = p[in_indices]
    inside_volumes = vol[in_indices]
    tree = cKDTree(inside_points)
    direction: dict[str, Any] = {}
    for delta in candidate_deltas_mm:
        neighborhoods = tree.query_ball_point(inside_points, delta)
        counts = []
        anis = []
        zeros = 0
        for i, ids in enumerate(neighborhoods):
            ids = [j for j in ids if j != i]
            counts.append(len(ids))
            if len(ids) < 3:
                zeros += 1
                continue
            fabric = tensor_from_vectors(inside_points[ids] - inside_points[i], inside_volumes[np.asarray(ids)])
            if fabric is None:
                zeros += 1
                continue
            eig = np.linalg.eigvalsh(fabric).clip(0, 1)
            anis.append(float(eig[-1] / max(eig[0], 1e-12)))
        direction[f"delta_{delta:g}mm_candidate_only"] = {
            "point_count_native_IN_only": int(len(inside_points)), "neighbor_count_excluding_self": qstats(np.asarray(counts)),
            "points_with_fewer_than_3_neighbors": int(zeros),
            "fabric_eigenvalue_max_over_min_for_centers_with_at_least_3_neighbors": qstats(np.asarray(anis)),
            "definition": "cell-volume-weighted second moment of unit centroid-neighbor directions; diagnostic radius only, not a PG horizon selection",
        }

    face_counts: dict[tuple[int, int, int], int] = defaultdict(int)
    for t in mesh["connectivity_ids"]:
        for tri in ((t[0], t[1], t[2]), (t[0], t[3], t[1]), (t[0], t[2], t[3]), (t[1], t[3], t[2])):
            face_counts[tuple(sorted(map(int, tri)))] += 1
    center_mask = roi if wall_center_mask is None else np.asarray(wall_center_mask, dtype=bool)
    fabric_reference = geometry_clipped_fabric_reference(
        mesh,
        solid,
        np.flatnonzero(center_mask),
        candidate_deltas_mm,
        monte_carlo_samples=monte_carlo_samples,
    )
    return {
        "mesh_path": mesh["mesh_path"], "mesh_sha256": mesh["mesh_sha256"],
        "target_size_mm_from_report": size_mm,
        "tetra_count": int(len(vol)), "node_count": int(len(mesh["node_ids"])),
        "signed_orientation": {"positive": int(np.count_nonzero(mesh["signed_volume"] > 0)), "negative": int(np.count_nonzero(mesh["signed_volume"] < 0)), "zero": int(np.count_nonzero(mesh["signed_volume"] == 0))},
        "tetra_volume_mm3": {"sum": float(vol.sum()), "mean": float(vol.mean()), "cv": float(vol.std() / vol.mean()), "quantiles": qstats(vol)},
        "edge_length_mm": qstats(mesh["edges"].reshape(-1)),
        "tetra_max_edge_mm": qstats(mesh["max_edge"]),
        "max_edge_over_target_size": qstats(mesh["max_edge"] / size_mm),
        "tetra_min_edge_mm": qstats(mesh["min_edge"]),
        "radius_ratio_circumradius_over_inradius_regular_equals_3": qstats(mesh["radius_ratio"]),
        "circumradius_over_min_edge_regular_0_612": qstats(mesh["radius_edge"]),
        "minimum_dihedral_angle_deg": qstats(mesh["min_dihedral"]),
        "sliver_counts": {"min_dihedral_lt_5deg": int(np.count_nonzero(mesh["min_dihedral"] < 5)), "min_dihedral_lt_10deg": int(np.count_nonzero(mesh["min_dihedral"] < 10)), "radius_ratio_gt_10": int(np.count_nonzero(mesh["radius_ratio"] > 10)), "radius_ratio_gt_100": int(np.count_nonzero(mesh["radius_ratio"] > 100))},
        "face_adjacency": {"internal_face_count": int(sum(v == 2 for v in face_counts.values())), "boundary_face_count": int(sum(v == 1 for v in face_counts.values())), "nonmanifold_face_count": int(sum(v > 2 for v in face_counts.values()))},
        "neighbor_spacing_mm_face_adjacent_centroid_to_centroid": qstats(mesh["neighbor_spacing"]),
        "neighbor_cell_volume_ratio": qstats(mesh["neighbor_volume_ratio"]),
        "wall_roi_definition": "x < -14.5 mm, |y| <= 3 mm, z=22..27 mm, angular distance from -X <=12 deg; actual native radial solid intersections at each centroid determine through-wall position",
        "wall_roi_centroid_count": int(np.count_nonzero(roi)),
        "wall_coordinate_outer_to_inner_0_to_1": qstats(np.asarray(normalized_wall)) if report_wall_distribution else None,
        "wall_coordinate_approximate_distinct_layers_10pct_bins": int(len(set(np.clip(np.floor(np.asarray(normalized_wall) * 10), 0, 9).astype(int)))) if normalized_wall else None,
        "directional_fabric": direction,
        "native_clipped_directional_reference": fabric_reference,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh-h8", type=Path, required=True)
    parser.add_argument("--mesh-h6", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    args = parser.parse_args()
    step = ASSET / "geometry/solid.step"
    actual_hash = sha256(step)
    if actual_hash != EXPECTED_SHA:
        raise ValueError(f"canonical STEP hash mismatch: {actual_hash}")
    solid = cq.importers.importStep(str(step)).val().Solids()[0]
    reports = []
    for size, file in ((8.0, args.mesh_h8), (6.0, args.mesh_h6)):
        reports.append(summarize_mesh(mesh_from_msh(file), size, solid))
    result = {
        "experiment_id": "pg.wall-resolution-existing-whole-mesh-audit/0.1.0",
        "source_step": str(step), "source_step_sha256": actual_hash,
        "scope": "read-only diagnosis of existing whole-CAD h8/h6 first-order tetra meshes; no solver, no re-meshing",
        "candidate_delta_mm": [4.0, 6.0],
        "candidate_delta_note": "exploratory direction-support diagnostics only; neither value is a selected/validated material horizon and delta is not mesh edge length; per-center native-solid-clipped Monte Carlo tensors separate physical CAD truncation from mesh-cloud direction error",
        "meshes": reports,
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps({"output_json": str(args.output_json), "source_step_sha256": actual_hash, "mesh_summaries": [{k: v for k, v in m.items() if k in ("target_size_mm_from_report", "tetra_count", "tetra_volume_mm3", "edge_length_mm", "radius_ratio_circumradius_over_inradius_regular_equals_3", "minimum_dihedral_angle_deg", "sliver_counts", "wall_roi_centroid_count", "wall_coordinate_outer_to_inner_0_to_1")} for m in reports]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
