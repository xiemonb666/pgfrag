"""One-shot, bounded native-CAD clipped thin-wall mesh trial.

This produces geometry/mesh diagnostics only. It never writes a PG solver job.
"""

from __future__ import annotations

import argparse
import hashlib
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
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.IntCurvesFace import IntCurvesFace_ShapeIntersector
from OCP.TopAbs import TopAbs_IN, TopAbs_ON, TopAbs_OUT
from OCP.gp import gp_Dir, gp_Lin, gp_Pnt

import mesh_audit as ma

ROOT = Path(__file__).resolve().parents[3]
ASSET = ROOT / "relicweave-synlab/tmp/review_20261003/implementation/vf/artifacts/four_bands_brep_semantics1_final"
STEP = ASSET / "geometry/solid.step"
OUTPUT_ROOT = ROOT / "relicweave-synlab/tmp/review_20261004/implementation/wall_resolution"
SDK = ROOT / "third_party/pgfrag-geometry-runtime/sdk/gmsh-4.15.2-Windows64-sdk/lib"
EXPECTED_SHA = "4623bd373aeed0acfacfbf6cc0dda65952339561779b43429c59d61b8b5a1c30"
ROI_BOX_MM = {"xmin": -23.0, "xmax": -14.0, "ymin": -1.4, "ymax": 1.4, "zmin": 23.1, "zmax": 25.9}
MAX_TETS = 10_000
MAX_WALL_SECONDS_PER_TRIAL = 180.0
FORECAST_MARGIN = 1.15
PREDECLARED_CANDIDATE_DELTA_MM = 1.0
INTERIOR_SAMPLE_COUNT = 512
INTERIOR_SAMPLE_SEED = 20261003


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def classify_one(checker: BRepClass3d_SolidClassifier, point: np.ndarray) -> str:
    try:
        checker.Perform(gp_Pnt(*map(float, point)), 1e-7)
        state = checker.State()
    except Exception:
        return "UNKNOWN"
    return "IN" if state == TopAbs_IN else "OUT" if state == TopAbs_OUT else "UNKNOWN" if state == TopAbs_ON else "UNKNOWN"


def finite_interior_samples(mesh: dict[str, Any]) -> np.ndarray:
    count = min(INTERIOR_SAMPLE_COUNT, int(mesh["tet_xyz"].shape[0]))
    selected = np.linspace(0, mesh["tet_xyz"].shape[0] - 1, count, dtype=np.int64)
    rng = np.random.default_rng(INTERIOR_SAMPLE_SEED)
    weights = rng.dirichlet(np.ones(4, dtype=np.float64), size=count)
    return np.einsum("ij,ijk->ik", weights, mesh["tet_xyz"][selected])


def direct_native_normal_thickness(solid: cq.Shape) -> dict[str, Any]:
    """Measure outer-to-opposing-inner span along native outer-face normals."""
    faces = solid.Faces()
    cone_faces = [(i, face) for i, face in enumerate(faces) if face.geomType() == "CONE"]
    z_values = np.linspace(23.12, 25.88, 7)
    theta_values = np.radians(np.linspace(-6.0, 6.0, 7)) + math.pi
    records: list[dict[str, Any]] = []
    for z in z_values:
        for theta in theta_values:
            direction = (math.cos(float(theta)), math.sin(float(theta)), 0.0)
            radial_line = gp_Lin(gp_Pnt(0.0, 0.0, float(z)), gp_Dir(*direction))
            radial_intersector = IntCurvesFace_ShapeIntersector()
            radial_intersector.Load(solid.wrapped, 1e-7)
            radial_intersector.Perform(radial_line, 0.0, 40.0)
            if not radial_intersector.IsDone():
                records.append({"z_mm": float(z), "theta_rad": float(theta), "status": "radial_intersection_failed"})
                continue
            hits = sorted(float(radial_intersector.WParameter(i)) for i in range(1, int(radial_intersector.NbPnt()) + 1))
            unique: list[float] = []
            for value in hits:
                if not unique or abs(value - unique[-1]) > 1e-6:
                    unique.append(value)
            if len(unique) < 2:
                records.append({"z_mm": float(z), "theta_rad": float(theta), "status": "fewer_than_two_native_wall_hits", "radial_hits_mm": unique})
                continue
            rin, rout = unique[-2], unique[-1]
            point = np.array([rout * direction[0], rout * direction[1], z], dtype=np.float64)
            vertex = cq.Vertex.makeVertex(*map(float, point))
            face_index, outer_face = min(cone_faces, key=lambda item: vertex.distance(item[1]))
            normal = np.asarray(outer_face.normalAt(tuple(map(float, point))).toTuple(), dtype=np.float64)
            if np.dot(normal[:2], direction[:2]) < 0:
                normal *= -1.0
            normal /= np.linalg.norm(normal)
            eps = 0.01
            shifted = point - eps * normal
            inward = -normal
            line = gp_Lin(gp_Pnt(*map(float, shifted)), gp_Dir(*map(float, inward)))
            intersector = IntCurvesFace_ShapeIntersector()
            intersector.Load(solid.wrapped, 1e-7)
            intersector.Perform(line, 0.0, 8.0)
            if not intersector.IsDone():
                records.append({"z_mm": float(z), "theta_rad": float(theta), "status": "normal_intersection_failed", "outer_face_index": face_index})
                continue
            candidates = sorted(float(intersector.WParameter(i)) for i in range(1, int(intersector.NbPnt()) + 1) if float(intersector.WParameter(i)) > 0.1)
            if not candidates:
                records.append({"z_mm": float(z), "theta_rad": float(theta), "status": "no_opposing_normal_hit", "outer_face_index": face_index})
                continue
            thickness = eps + candidates[0]
            records.append({
                "z_mm": float(z), "theta_rad": float(theta), "status": "measured_native_normal_hit",
                "outer_face_index": int(face_index), "outer_point_xyz_mm": point.tolist(),
                "outer_unit_normal": normal.tolist(), "opposing_inner_point_xyz_mm": (shifted + candidates[0] * inward).tolist(),
                "normal_thickness_mm": float(thickness), "same_z_radial_thickness_mm": float(rout - rin),
                "normal_hit_distance_mm_from_0p01_inside_start": float(candidates[0]),
            })
    valid = [record for record in records if record.get("status") == "measured_native_normal_hit"]
    thicknesses = np.asarray([record["normal_thickness_mm"] for record in valid], dtype=np.float64)
    radial_values = np.asarray([record["same_z_radial_thickness_mm"] for record in valid], dtype=np.float64)
    return {
        "method": "native STEP B-rep radial ray locates the outer surface; native cone face normal starts 0.01mm inside and IntCurvesFace_ShapeIntersector finds the opposing inner B-rep hit; no ROI-box dimension is used",
        "measurement_count": int(len(records)), "successful_native_normal_measurements": int(len(valid)),
        "normal_thickness_mm": ma.qstats(thicknesses), "same_z_radial_thickness_mm": ma.qstats(radial_values),
        "samples": records,
    }


def _forecast_tets(exact_volume_mm3: float, mesh_size_mm: float) -> float:
    existing_density_per_target_cube = max(
        (1313.0 / 25419.097419858103) * 8.0**3,
        (2030.0 / 25528.120818621566) * 6.0**3,
    )
    return float(exact_volume_mm3 * existing_density_per_target_cube / mesh_size_mm**3 * FORECAST_MARGIN)


def _load_gmsh() -> Any:
    if not SDK.is_dir():
        raise FileNotFoundError(f"isolated Gmsh SDK not found: {SDK}")
    if hasattr(os, "add_dll_directory"):
        os.add_dll_directory(str(SDK))
    if str(SDK) not in sys.path:
        sys.path.insert(0, str(SDK))
    import gmsh  # type: ignore

    return gmsh


def native_roi_preflight(output_dir: Path) -> dict[str, Any]:
    """Measure the fixed native CAD common and forecast both fixed mesh sizes before meshing."""
    output_dir.mkdir(parents=True, exist_ok=True)
    step_hash = sha256(STEP)
    if step_hash != EXPECTED_SHA:
        raise ValueError(f"canonical STEP SHA mismatch: expected={EXPECTED_SHA} actual={step_hash}")
    manifest = json.loads((ASSET / "manifest.json").read_text(encoding="utf-8"))
    expected_manifest_hash = str(manifest["files"]["geometry/solid.step"]["sha256"])
    if expected_manifest_hash.lower() != step_hash.lower():
        raise ValueError(f"manifest STEP SHA mismatch: expected={expected_manifest_hash} actual={step_hash}")

    gmsh = _load_gmsh()
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("Geometry.OCCScaling", 1.0)
        gmsh.model.add("wall_resolution_fixed_roi_exact_volume_preflight")
        imported = gmsh.model.occ.importShapes(str(STEP), True, "step")
        gmsh.model.occ.synchronize()
        source_volumes = gmsh.model.getEntities(3)
        if len(source_volumes) != 1 or len(imported) != 1:
            raise RuntimeError(f"canonical STEP must import as one native solid: imported={imported} volumes={source_volumes}")
        box_tag = gmsh.model.occ.addBox(
            ROI_BOX_MM["xmin"], ROI_BOX_MM["ymin"], ROI_BOX_MM["zmin"],
            ROI_BOX_MM["xmax"] - ROI_BOX_MM["xmin"],
            ROI_BOX_MM["ymax"] - ROI_BOX_MM["ymin"],
            ROI_BOX_MM["zmax"] - ROI_BOX_MM["zmin"],
        )
        common, _ = gmsh.model.occ.intersect(
            [(3, int(source_volumes[0][1]))], [(3, int(box_tag))],
            removeObject=True, removeTool=True,
        )
        gmsh.model.occ.synchronize()
        roi_volumes = gmsh.model.getEntities(3)
        if len(roi_volumes) != 1 or len(common) != 1 or int(common[0][0]) != 3:
            raise RuntimeError(f"fixed CAD∩box ROI must produce one connected native volume: common={common} volumes={roi_volumes}")
        volume_tag = int(roi_volumes[0][1])
        exact_volume = float(gmsh.model.occ.getMass(3, volume_tag))
        roi_bbox = list(map(float, gmsh.model.occ.getBoundingBox(3, volume_tag)))
        if not math.isfinite(exact_volume) or exact_volume <= 0:
            raise RuntimeError(f"invalid native OCC common volume: {exact_volume}")
    finally:
        gmsh.finalize()

    forecasts = {
        str(size): _forecast_tets(exact_volume, size)
        for size in (1.0, 0.5)
    }
    within = all(value < MAX_TETS for value in forecasts.values())
    report = {
        "experiment_id": "pg.wall-resolution-local-native-cad-static-mesh/0.1.0",
        "status": "PREFLIGHT_WITHIN_BOTH_FIXED_TET_BUDGETS" if within else "STOPPED_PREFLIGHT_EXACT_VOLUME_FORECAST_AT_OR_OVER_TET_LIMIT",
        "source_step": str(STEP), "source_step_sha256": step_hash,
        "python_executable": sys.executable, "argv": sys.argv,
        "roi_box_mm": ROI_BOX_MM,
        "roi_construction": "native OpenCASCADE CAD common of the canonical STEP solid and the fixed box; no triangulated clipping or changed ROI",
        "roi_exact_occ_volume_mm3": exact_volume, "roi_occ_bbox_mm": roi_bbox,
        "max_tetra_cells_per_trial": MAX_TETS,
        "forecast_margin_fraction": FORECAST_MARGIN,
        "preflight_forecast_tetra_count_by_fixed_h_mm": forecasts,
        "preflight_gate": "both forecasts must be strictly below 10000; otherwise neither 3D mesh trial is launched",
        "meshing_started": False,
    }
    (output_dir / "preflight.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    return report


def boundary_quality(gmsh: Any, volume_tag: int) -> dict[str, Any]:
    boundary = gmsh.model.getBoundary([(3, int(volume_tag))], combined=False, oriented=False, recursive=False)
    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    node_tags = np.asarray(node_tags, dtype=np.int64)
    xyz = np.asarray(coords, dtype=np.float64).reshape(-1, 3)
    rows = {int(tag): i for i, tag in enumerate(node_tags)}
    records = []
    for dim, tag in boundary:
        if int(dim) != 2:
            continue
        tag = int(tag)
        geom_type = str(gmsh.model.getType(2, tag))
        box = np.asarray(gmsh.model.getBoundingBox(2, tag), dtype=np.float64)
        plane_label = None
        if geom_type == "Plane":
            for axis, low_key, high_key, label_low, label_high in (
                (0, "xmin", "xmax", "x cut", "x cut"),
                (1, "ymin", "ymax", "negative-y cut", "positive-y cut"),
                (2, "zmin", "zmax", "lower-z cut", "upper-z cut"),
            ):
                if abs(box[axis] - ROI_BOX_MM[low_key]) <= 2e-5 and abs(box[axis + 3] - ROI_BOX_MM[low_key]) <= 2e-5:
                    plane_label = label_low
                if abs(box[axis] - ROI_BOX_MM[high_key]) <= 2e-5 and abs(box[axis + 3] - ROI_BOX_MM[high_key]) <= 2e-5:
                    plane_label = label_high
        group = "artificial_ROI_box_cut" if plane_label else "native_CAD_boundary"
        types, tags_by_type, nodes_by_type = gmsh.model.mesh.getElements(2, tag)
        group_rec = {"surface_tag": tag, "gmsh_geometric_type": geom_type, "classification": group, "cut_label": plane_label, "bbox_mm": box.tolist(), "triangle_count": 0, "other_surface_element_types": [], "surface_area_mm2": 0.0, "triangle_min_angle_deg": []}
        for etype, etags, conn_flat in zip(types, tags_by_type, nodes_by_type):
            props = gmsh.model.mesh.getElementProperties(int(etype))
            num_nodes = int(props[3])
            name = str(props[0])
            conn_ids = np.asarray(conn_flat, dtype=np.int64).reshape(-1, num_nodes)
            if name not in ("Triangle 3", "Triangle") or num_nodes != 3:
                group_rec["other_surface_element_types"].append({"type": int(etype), "name": name, "count": int(len(conn_ids))})
                continue
            xyz_tri = np.asarray([[xyz[rows[int(node)]] for node in tri] for tri in conn_ids], dtype=np.float64)
            edge = np.stack([np.linalg.norm(xyz_tri[:, 1] - xyz_tri[:, 0], axis=1), np.linalg.norm(xyz_tri[:, 2] - xyz_tri[:, 1], axis=1), np.linalg.norm(xyz_tri[:, 0] - xyz_tri[:, 2], axis=1)], axis=1)
            area = 0.5 * np.linalg.norm(np.cross(xyz_tri[:, 1] - xyz_tri[:, 0], xyz_tri[:, 2] - xyz_tri[:, 0]), axis=1)
            angles = []
            for corner in range(3):
                a = xyz_tri[:, (corner + 1) % 3] - xyz_tri[:, corner]
                b = xyz_tri[:, (corner + 2) % 3] - xyz_tri[:, corner]
                cosine = np.einsum("ij,ij->i", a, b) / np.maximum(np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1), 1e-300)
                angles.append(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))
            group_rec["triangle_count"] += int(len(conn_ids))
            group_rec["surface_area_mm2"] += float(area.sum())
            group_rec["triangle_min_angle_deg"].extend(np.min(np.stack(angles, axis=1), axis=1).tolist())
        records.append(group_rec)
    summarized: dict[str, Any] = {}
    for group in ("native_CAD_boundary", "artificial_ROI_box_cut"):
        selected = [record for record in records if record["classification"] == group]
        summarized[group] = {
            "surface_count": int(len(selected)), "triangle_count": int(sum(r["triangle_count"] for r in selected)),
            "surface_area_mm2": float(sum(r["surface_area_mm2"] for r in selected)),
            "triangle_min_angle_deg": ma.qstats(np.asarray([a for r in selected for a in r["triangle_min_angle_deg"]], dtype=np.float64)),
        }
    return {"method": "boundary surface triangles are grouped as native curved CAD faces or planar ROI box cuts by native OCC entity type and exact ROI-plane coordinate; cut surfaces are explicitly artificial patch boundaries", "groups": summarized, "surfaces": records}


def normal_wall_distribution(
    mesh: dict[str, Any],
    solid: cq.Shape,
    centroid_labels: list[str],
    native_normal_report: dict[str, Any],
    mesh_size_mm: float,
) -> tuple[dict[str, Any], np.ndarray, np.ndarray]:
    """Map mesh centroids to native outer-normal depth across the CAD wall."""
    points = mesh["centroids"]
    faces = solid.Faces()
    measured = [record for record in native_normal_report["samples"] if record.get("status") == "measured_native_normal_hit"]
    if not measured:
        raise RuntimeError("native normal wall distribution requires at least one direct opposing-face measurement")
    depth_by_cell = np.full(len(points), np.nan, dtype=np.float64)
    fraction_by_cell = np.full(len(points), np.nan, dtype=np.float64)
    candidate_mask = np.asarray([
        float(point[0]) < -14.5 and abs(float(point[1])) <= 1.45 and 23.05 <= float(point[2]) <= 25.95
        for point in points
    ], dtype=bool)
    failed = Counter()
    records: list[dict[str, Any]] = []
    for cell_id in np.flatnonzero(candidate_mask):
        point = points[int(cell_id)]
        if centroid_labels[int(cell_id)] != "IN":
            failed["centroid_not_native_IN"] += 1
            continue
        theta = math.atan2(float(point[1]), float(point[0]))
        hits = ma.radial_hits(solid, float(point[2]), theta)
        radius = float(np.hypot(point[0], point[1]))
        interval = None
        for inner, outer in zip(hits[::2], hits[1::2]):
            if inner - 1e-3 <= radius <= outer + 1e-3:
                interval = (float(inner), float(outer))
                break
        if interval is None:
            failed["no_native_radial_material_interval"] += 1
            continue
        inner_radius, outer_radius = interval
        radial_direction = np.asarray([math.cos(theta), math.sin(theta), 0.0], dtype=np.float64)
        outer_point = np.asarray([outer_radius * radial_direction[0], outer_radius * radial_direction[1], float(point[2])], dtype=np.float64)
        vertex = cq.Vertex.makeVertex(*map(float, outer_point))
        try:
            face_index, outer_face = min(
                ((index, face) for index, face in enumerate(faces) if face.geomType() == "CONE"),
                key=lambda item: vertex.distance(item[1]),
            )
            normal = np.asarray(outer_face.normalAt(tuple(map(float, outer_point))).toTuple(), dtype=np.float64)
        except Exception:
            failed["outer_cone_normal_unavailable"] += 1
            continue
        radial_cosine = abs(float(np.dot(normal[:2], radial_direction[:2]) / max(np.linalg.norm(normal), 1e-300)))
        if not math.isfinite(radial_cosine) or radial_cosine <= 0.0:
            failed["invalid_radial_normal_projection"] += 1
            continue
        depth = max(0.0, (outer_radius - radius) * radial_cosine)
        angular_diffs = [
            math.atan2(math.sin(theta - float(record["theta_rad"])), math.cos(theta - float(record["theta_rad"])))
            for record in measured
        ]
        nearest_record = min(
            zip(measured, angular_diffs),
            key=lambda item: ((float(point[2]) - float(item[0]["z_mm"])) / 2.76) ** 2 + (item[1] / math.radians(12.0)) ** 2,
        )[0]
        normal_span = float(nearest_record["normal_thickness_mm"])
        fraction = depth / normal_span if normal_span > 0.0 else math.nan
        depth_by_cell[int(cell_id)] = depth
        fraction_by_cell[int(cell_id)] = fraction
        records.append({
            "cell_id": int(mesh["cell_ids"][int(cell_id)]),
            "centroid_xyz_mm": point.tolist(),
            "native_radial_interval_mm": [inner_radius, outer_radius],
            "outer_cone_face_index": int(face_index),
            "outer_normal_radial_projection_cosine": radial_cosine,
            "normal_depth_from_outer_mm": depth,
            "local_direct_native_normal_span_mm": normal_span,
            "outer_to_inner_normal_fraction": fraction,
            "normal_span_sample_z_mm": float(nearest_record["z_mm"]),
            "normal_span_sample_theta_rad": float(nearest_record["theta_rad"]),
        })

    valid_depths = depth_by_cell[np.isfinite(depth_by_cell)]
    valid_fractions = fraction_by_cell[np.isfinite(fraction_by_cell)]
    native_in_candidate_count = int(sum(centroid_labels[int(index)] == "IN" for index in np.flatnonzero(candidate_mask)))
    bins = np.linspace(0.0, 1.0, 11)
    counts, _ = np.histogram(np.clip(valid_fractions, 0.0, 1.0), bins=bins) if len(valid_fractions) else (np.zeros(10, dtype=np.int64), bins)
    thickness = native_normal_report["normal_thickness_mm"]
    wall_distribution = {
        "method": "per-cell centroid native radial interval; the outer CAD cone normal at each cell's outer radial hit projects outer-to-centroid radial depth onto the local normal; each depth fraction is normalized by the nearest direct native outer-to-opposing-inner normal measurement on the canonical STEP",
        "candidate_cell_count": int(np.count_nonzero(candidate_mask)),
        "native_IN_cell_count": native_in_candidate_count,
        "measured_cell_count": int(len(records)),
        "failed_cell_counts_by_reason": dict(failed),
        "direct_native_normal_span_mm": thickness,
        "direct_native_normal_sample_count": int(len(measured)),
        "native_normal_wall_thickness_over_mesh_target_h_candidate_only": float(thickness["median"] / mesh_size_mm),
        "outer_to_inner_normal_depth_mm": ma.qstats(valid_depths),
        "outer_to_inner_normal_fraction_0_to_1": ma.qstats(valid_fractions),
        "normal_fraction_10pct_band_counts": [
            {"band": f"{bins[i]:.1f}-{bins[i + 1]:.1f}", "centroid_count": int(counts[i])}
            for i in range(10)
        ],
        "normal_fraction_10pct_bands_occupied": [int(index) for index in np.flatnonzero(counts)],
        "normal_fraction_outside_0_to_1_count_before_histogram_clipping": int(np.count_nonzero((valid_fractions < 0.0) | (valid_fractions > 1.0))),
        "per_cell_records": records,
    }
    return wall_distribution, depth_by_cell, fraction_by_cell


def _section_polygon(tet_xyz: np.ndarray, plane_y_mm: float = 0.0) -> np.ndarray | None:
    intersections: list[np.ndarray] = []
    edges = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
    tol = 1e-10
    for left, right in edges:
        p0, p1 = tet_xyz[left], tet_xyz[right]
        d0, d1 = float(p0[1] - plane_y_mm), float(p1[1] - plane_y_mm)
        if abs(d0) <= tol:
            intersections.append(p0[[0, 2]])
        if abs(d1) <= tol:
            intersections.append(p1[[0, 2]])
        if d0 * d1 < -(tol * tol):
            weight = -d0 / (d1 - d0)
            point = p0 + weight * (p1 - p0)
            intersections.append(point[[0, 2]])
    unique: list[np.ndarray] = []
    for point in intersections:
        if not any(float(np.linalg.norm(point - known)) <= 1e-8 for known in unique):
            unique.append(point)
    if len(unique) < 3:
        return None
    polygon = np.asarray(unique, dtype=np.float64)
    center = polygon.mean(axis=0)
    angles = np.arctan2(polygon[:, 1] - center[1], polygon[:, 0] - center[0])
    return polygon[np.argsort(angles)]


def render_section_png(
    mesh: dict[str, Any],
    solid: cq.Shape,
    normal_fraction_by_cell: np.ndarray,
    mesh_size_mm: float,
    output_path: Path,
) -> dict[str, Any]:
    """Render the actual tetrahedral mesh intersection with the interior y=0 plane."""
    polygons: list[np.ndarray] = []
    values: list[float] = []
    for cell_index, tet in enumerate(mesh["tet_xyz"]):
        polygon = _section_polygon(tet, 0.0)
        if polygon is not None:
            polygons.append(polygon)
            values.append(float(normal_fraction_by_cell[cell_index]))

    fig, ax = plt.subplots(figsize=(9, 6), constrained_layout=True)
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#c9cdd2")
    collection = PolyCollection(
        polygons,
        array=np.ma.masked_invalid(np.asarray(values, dtype=np.float64)),
        cmap=cmap, clim=(0.0, 1.0), edgecolors=(0.08, 0.10, 0.13, 0.22), linewidths=0.12,
    )
    ax.add_collection(collection)

    z_samples = np.linspace(ROI_BOX_MM["zmin"] + 0.01, ROI_BOX_MM["zmax"] - 0.01, 64)
    outer_x: list[float] = []
    inner_x: list[float] = []
    surface_z: list[float] = []
    for z in z_samples:
        hits = ma.radial_hits(solid, float(z), math.pi)
        if len(hits) >= 2:
            inner_x.append(-float(hits[-2]))
            outer_x.append(-float(hits[-1]))
            surface_z.append(float(z))
    if surface_z:
        ax.plot(outer_x, surface_z, color="#101820", linewidth=1.2, label="native STEP outer wall")
        ax.plot(inner_x, surface_z, color="#101820", linewidth=1.2, linestyle="--", label="native STEP inner wall")
    ax.axvline(ROI_BOX_MM["xmin"], color="#b34040", linestyle=":", linewidth=1.1, label="artificial ROI cut plane")
    ax.axvline(ROI_BOX_MM["xmax"], color="#b34040", linestyle=":", linewidth=1.1)
    ax.axhline(ROI_BOX_MM["zmin"], color="#b34040", linestyle=":", linewidth=1.1)
    ax.axhline(ROI_BOX_MM["zmax"], color="#b34040", linestyle=":", linewidth=1.1)
    ax.set_xlim(ROI_BOX_MM["xmin"] - 0.3, ROI_BOX_MM["xmax"] + 0.3)
    ax.set_ylim(ROI_BOX_MM["zmin"] - 0.15, ROI_BOX_MM["zmax"] + 0.15)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("z (mm)")
    ax.set_title(f"Actual tetra mesh section at y = 0 mm, target h = {mesh_size_mm:g} mm")
    colorbar = fig.colorbar(collection, ax=ax, shrink=0.9)
    colorbar.set_label("outer-to-inner native normal depth fraction (0 to 1)")
    ax.legend(loc="lower right", fontsize=8)
    fig.text(0.02, 0.005, "The red boundaries are artificial CAD∩box patch cuts. This interior section is a mesh diagnostic, not a crack surface.", fontsize=8)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    return {
        "path": str(output_path),
        "sha256": sha256(output_path),
        "section_plane": "y=0 mm, intersected with the actual first-order tetrahedral mesh",
        "cross_section_polygon_count": int(len(polygons)),
        "native_step_outer_inner_boundary_sample_count": int(len(surface_z)),
        "artificial_cut_surfaces_are_cracks": False,
    }


def run_trial(mesh_size_mm: float, output_dir: Path) -> dict[str, Any]:
    started = time.perf_counter()
    output_dir.mkdir(parents=True, exist_ok=True)
    step_hash = sha256(STEP)
    if step_hash != EXPECTED_SHA:
        raise ValueError(f"canonical STEP SHA mismatch: expected={EXPECTED_SHA} actual={step_hash}")
    manifest = json.loads((ASSET / "manifest.json").read_text(encoding="utf-8"))
    cad_solid = cq.importers.importStep(str(STEP)).val().Solids()[0]
    native_normal_report = direct_native_normal_thickness(cad_solid)
    status = "NOT_STARTED"
    gmsh_mesh_file = output_dir / f"roi_gmsh_h{mesh_size_mm:g}mm.msh"
    gmsh_mesh_sha = None
    tet_count = 0
    exact_roi_volume = None
    tet_volumes: np.ndarray | None = None
    elapsed_mesh_s = 0.0
    boundary = None
    mesher = None
    if not (0.0 < mesh_size_mm <= 1.0) or mesh_size_mm not in (0.5, 1.0):
        raise ValueError("only the two predeclared ROI target sizes 1.0 mm and 0.5 mm are allowed")
    if not native_normal_report["successful_native_normal_measurements"]:
        raise RuntimeError("native opposing-face normal thickness could not be measured; stop before meshing")
    gmsh = _load_gmsh()
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 1)
        gmsh.option.setNumber("General.NumThreads", 1)
        gmsh.option.setNumber("Mesh.MaxNumThreads1D", 1)
        gmsh.option.setNumber("Mesh.MaxNumThreads2D", 1)
        gmsh.option.setNumber("Mesh.MaxNumThreads3D", 1)
        gmsh.option.setNumber("Geometry.OCCImportLabels", 0)
        gmsh.option.setNumber("Geometry.OCCScaling", 1.0)
        gmsh.option.setNumber("Mesh.MeshSizeMin", mesh_size_mm)
        gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_size_mm)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)
        gmsh.option.setNumber("Mesh.ElementOrder", 1)
        gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
        gmsh.option.setNumber("Mesh.Binary", 0)
        gmsh.model.add(f"wall_resolution_roi_h{mesh_size_mm:g}mm")
        imported = gmsh.model.occ.importShapes(str(STEP), True, "step")
        gmsh.model.occ.synchronize()
        source_volumes = gmsh.model.getEntities(3)
        if len(source_volumes) != 1 or len(imported) != 1:
            raise RuntimeError(f"native source STEP volume topology unexpected: imported={imported} volumes={source_volumes}")
        box_tag = gmsh.model.occ.addBox(
            ROI_BOX_MM["xmin"], ROI_BOX_MM["ymin"], ROI_BOX_MM["zmin"],
            ROI_BOX_MM["xmax"] - ROI_BOX_MM["xmin"],
            ROI_BOX_MM["ymax"] - ROI_BOX_MM["ymin"],
            ROI_BOX_MM["zmax"] - ROI_BOX_MM["zmin"],
        )
        common, _ = gmsh.model.occ.intersect([(3, int(source_volumes[0][1]))], [(3, int(box_tag))], removeObject=True, removeTool=True)
        gmsh.model.occ.synchronize()
        roi_volumes = gmsh.model.getEntities(3)
        if len(roi_volumes) != 1 or len(common) != 1 or int(common[0][0]) != 3:
            raise RuntimeError(f"ROI common must produce one connected native CAD volume; common={common} volumes={roi_volumes}")
        volume_tag = int(roi_volumes[0][1])
        exact_roi_volume = float(gmsh.model.occ.getMass(3, volume_tag))
        roi_bbox = list(map(float, gmsh.model.occ.getBoundingBox(3, volume_tag)))
        if not math.isfinite(exact_roi_volume) or exact_roi_volume <= 0:
            raise RuntimeError(f"invalid native OCC common volume: {exact_roi_volume}")
        # Conservative preflight from the preserved existing whole-vessel runs.
        # It never alters the fixed ROI, h, cell limit, or acceptance threshold.
        forecast_tets = _forecast_tets(exact_roi_volume, mesh_size_mm)
        if forecast_tets >= MAX_TETS:
            status = "STOPPED_PREFLIGHT_FORECAST_AT_OR_OVER_CELL_LIMIT"
            report = {
                "experiment_id": "pg.wall-resolution-local-native-cad-static-mesh/0.1.0",
                "status": status,
                "python_executable": sys.executable,
                "argv": sys.argv,
                "mesh_size_mm": mesh_size_mm,
                "source_step": str(STEP), "source_step_sha256": step_hash,
                "roi_box_mm": ROI_BOX_MM, "roi_exact_occ_volume_mm3": exact_roi_volume,
                "roi_occ_bbox_mm": roi_bbox, "preflight_forecast_tetra_count_with_15pct_margin": forecast_tets,
                "max_tetra_cells_per_trial": MAX_TETS,
                "meshing_started": False,
                "note": "No Gmsh 3D mesh was generated because the fixed ROI and unchanged target size forecast at/over the budget; ROI/h was not changed and no retry was attempted.",
                "native_normal_thickness": native_normal_report,
            }
            (output_dir / "report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
            return report
        print(f"WALL_RESOLUTION_PREFLIGHT roi_volume_mm3={exact_roi_volume:.12f} mesh_size_mm={mesh_size_mm:g} forecast_tets_with_margin={forecast_tets:.1f} cap={MAX_TETS}")
        if mesh_size_mm == 0.5:
            gmsh.option.setNumber("General.Terminal", 1)
        tic = time.perf_counter()
        gmsh.model.mesh.generate(3)
        elapsed_mesh_s = time.perf_counter() - tic
        types, tags_by_type, _ = gmsh.model.mesh.getElements(3, volume_tag)
        tet_type = int(gmsh.model.mesh.getElementType("tetrahedron", 1))
        if any(int(element_type) != tet_type for element_type in types):
            raise RuntimeError(f"ROI mesh produced non-linear-tetra volume elements: {list(map(int, types))}")
        tet_count = int(sum(len(tags) for tags in tags_by_type))
        gmsh.write(str(gmsh_mesh_file))
        gmsh_mesh_sha = sha256(gmsh_mesh_file)
        boundary = boundary_quality(gmsh, volume_tag)
        if tet_count > MAX_TETS:
            status = "FAIL_CELL_LIMIT_EXCEEDED_FIRST_ATTEMPT_PRESERVED"
        else:
            status = "MESH_GENERATED_WITHIN_CELL_LIMIT"
        mesher = {
            "name": "Gmsh", "version": str(gmsh.option.getString("General.Version")) if hasattr(gmsh.option, "getString") else "4.15.2",
            "kernel": "native OpenCASCADE STEP import and CAD∩box common", "sdk": str(SDK),
            "mesh_algorithm_3d": 1, "element_order": 1,
            "general_num_threads": int(gmsh.option.getNumber("General.NumThreads")),
            "mesh_max_num_threads_1d": int(gmsh.option.getNumber("Mesh.MaxNumThreads1D")),
            "mesh_max_num_threads_2d": int(gmsh.option.getNumber("Mesh.MaxNumThreads2D")),
            "mesh_max_num_threads_3d": int(gmsh.option.getNumber("Mesh.MaxNumThreads3D")),
            "target_size_mm": float(mesh_size_mm), "mesh_size_min_mm": float(gmsh.option.getNumber("Mesh.MeshSizeMin")),
            "mesh_size_max_mm": float(gmsh.option.getNumber("Mesh.MeshSizeMax")),
            "mesh_size_from_curvature": 0, "mesh_size_from_points": 0, "mesh_size_extend_from_boundary": 0,
            "mesh_generation_wall_seconds": float(elapsed_mesh_s), "max_wall_seconds_per_trial": MAX_WALL_SECONDS_PER_TRIAL,
            "tetra_count": tet_count, "node_count": int(len(gmsh.model.mesh.getNodes()[0])),
            "mesh_file": str(gmsh_mesh_file), "mesh_sha256": gmsh_mesh_sha,
            "only_one_initial_mesh_attempt": True,
        }
    finally:
        gmsh.finalize()

    mesh = ma.mesh_from_msh(gmsh_mesh_file)
    mesh_volume = float(np.sum(mesh["volume"], dtype=np.float64))
    relative_volume_error = float((mesh_volume - float(exact_roi_volume)) / float(exact_roi_volume))
    source_classifier = BRepClass3d_SolidClassifier()
    source_classifier.Load(cad_solid.wrapped)
    centroid_labels = [classify_one(source_classifier, point) for point in mesh["centroids"]]
    samples = finite_interior_samples(mesh)
    sample_labels = [classify_one(source_classifier, point) for point in samples]
    points = mesh["centroids"]
    cut_plane_clearance_mm = np.min(np.stack((
        points[:, 0] - ROI_BOX_MM["xmin"], ROI_BOX_MM["xmax"] - points[:, 0],
        points[:, 1] - ROI_BOX_MM["ymin"], ROI_BOX_MM["ymax"] - points[:, 1],
        points[:, 2] - ROI_BOX_MM["zmin"], ROI_BOX_MM["zmax"] - points[:, 2],
    ), axis=1), axis=1)
    core_mask_before_cut_clearance = np.asarray([
        float(point[0]) < -14.5 and abs(float(point[1])) <= 0.3 and 24.2 <= float(point[2]) <= 24.8
        for point in points
    ], dtype=bool)
    wall_mask = core_mask_before_cut_clearance & (cut_plane_clearance_mm >= 1.1)
    wall_dist, normal_depth_by_cell, normal_fraction_by_cell = normal_wall_distribution(
        mesh, cad_solid, centroid_labels, native_normal_report, mesh_size_mm,
    )
    core_clearance = cut_plane_clearance_mm[wall_mask]
    wall_dist["delta_1mm_diagnostic_core"] = {
        "definition": "centroid x < -14.5 mm, |y| <= 0.3 mm, z=24.2..24.8 mm, and distance >=1.1 mm from every one of the six artificial ROI cut planes",
        "diagnostic_delta_mm_candidate_only": PREDECLARED_CANDIDATE_DELTA_MM,
        "diagnostic_center_count_before_all_six_cut_plane_clearance": int(np.count_nonzero(core_mask_before_cut_clearance)),
        "eligible_centers_after_all_six_cut_plane_clearance": int(np.count_nonzero(wall_mask)),
        "all_eligible_center_clearances_at_least_1p1mm": bool(not len(core_clearance) or np.all(core_clearance >= 1.1)),
        "eligible_center_clearance_to_nearest_artificial_cut_plane_mm": ma.qstats(core_clearance),
        "role": "local directional-fabric diagnostic core only; delta is neither a selected material horizon nor a calibrated length",
    }
    section_png = render_section_png(
        mesh, cad_solid, normal_fraction_by_cell, mesh_size_mm,
        output_dir / f"actual_section_y0_h{mesh_size_mm:g}mm.png",
    )
    local_candidate_deltas = (PREDECLARED_CANDIDATE_DELTA_MM,)
    quality = ma.summarize_mesh(
        mesh, mesh_size_mm, cad_solid,
        candidate_deltas_mm=local_candidate_deltas,
        wall_center_mask=wall_mask,
        monte_carlo_samples=5000,
        report_wall_distribution=False,
    )
    state_counts = Counter(centroid_labels)
    sample_counts = Counter(sample_labels)
    checks = {
        "within_10000_tetra_budget": bool(tet_count <= MAX_TETS),
        "all_tetra_signed_orientations_positive": bool(np.all(mesh["signed_volume"] > 0)),
        "all_cell_volumes_finite_positive": bool(np.isfinite(mesh["volume"]).all() and np.all(mesh["volume"] > 0)),
        "absolute_mesh_volume_error_at_most_1pct_unchanged_gate": bool(abs(relative_volume_error) <= 0.01),
        "all_centroids_native_source_IN": bool(state_counts["IN"] == tet_count),
        "all_finite_stratified_samples_native_source_IN": bool(sample_counts["IN"] == len(samples)),
        "all_delta_core_centers_at_least_1p1mm_from_each_roi_cut_plane": bool(wall_dist["delta_1mm_diagnostic_core"]["all_eligible_center_clearances_at_least_1p1mm"]),
        "actual_y0_mesh_section_png_written": bool(Path(section_png["path"]).is_file()),
        "no_solver_or_fracture_surface_created": True,
    }
    if status != "FAIL_CELL_LIMIT_EXCEEDED_FIRST_ATTEMPT_PRESERVED":
        status = "PASS_LOCAL_GEOMETRY_AND_MESH_CHECKS" if all(checks.values()) else "FAIL_LOCAL_GEOMETRY_AND_MESH_CHECKS"
    npz = output_dir / "roi_mesh_arrays.npz"
    provenance = {
        "experiment_id": "pg.wall-resolution-local-native-cad-static-mesh/0.1.0",
        "source_step_sha256": step_hash, "roi_box_mm": ROI_BOX_MM,
        "mesh_size_target_mm": mesh_size_mm, "max_volume_tets": MAX_TETS,
        "diagnostic_delta_candidate_mm": PREDECLARED_CANDIDATE_DELTA_MM,
        "diagnostic_delta_is_selected_material_horizon": False,
        "roi_clip_surfaces_are_fracture_surfaces": False,
        "solver_job_created": False,
    }
    np.savez_compressed(
        npz,
        points_mm=mesh["centroids"], cell_volume_mm3=mesh["volume"],
        sourcecellIDs=mesh["cell_ids"], source_node_ids=mesh["node_ids"],
        tetra_node_ids=mesh["connectivity_ids"], node_coordinates_mm=mesh["xyz"],
        normal_depth_from_outer_mm=normal_depth_by_cell,
        outer_to_inner_normal_fraction=normal_fraction_by_cell,
        nearest_roi_cut_plane_clearance_mm=cut_plane_clearance_mm,
        source_step_sha256=np.asarray(step_hash), mesh_sha256=np.asarray(gmsh_mesh_sha),
        mesh_size_target_mm=np.asarray(mesh_size_mm), roi_box_mm=np.asarray([ROI_BOX_MM[k] for k in ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax")]),
        provenance_json=np.asarray(json.dumps(provenance, sort_keys=True)),
    )
    report = {
        "experiment_id": "pg.wall-resolution-local-native-cad-static-mesh/0.1.0",
        "status": status,
        "python_executable": sys.executable,
        "argv": sys.argv,
        "source": {"step_path": str(STEP), "step_sha256": step_hash, "manifest_expected_sha256": EXPECTED_SHA, "hash_check": "PASS", "length_unit": "mm", "frame": manifest["frame"]},
        "roi": {
            "definition": "one native CAD sidewall patch opposite the +X handle, on the z20–28 mm conical meridian face; Gmsh OpenCASCADE common(source STEP solid, native box)",
            "roi_box_mm": ROI_BOX_MM, "actual_material_roi_bbox_mm": roi_bbox,
            "exact_native_occ_common_volume_mm3": float(exact_roi_volume),
            "physical_wall_thickness_mm": native_normal_report["normal_thickness_mm"],
            "direct_native_normal_measurements": native_normal_report,
            "artificial_cut_surfaces": "planar faces created wherever the fixed box clips the canonical CAD solid; these are finite-patch boundaries for this local mesh, never fracture/R2 surfaces",
            "normal_wall_measurement_is_not_box_width": True,
        },
        "trial": {
            "mesh_size_target_mm": float(mesh_size_mm), "tetra_count": int(tet_count), "node_count": int(mesher["node_count"] if mesher else len(mesh["node_ids"])),
            "mesh_generation_wall_seconds": float(elapsed_mesh_s), "max_wall_seconds_per_trial": MAX_WALL_SECONDS_PER_TRIAL,
            "cell_budget_max": MAX_TETS, "conservative_preflight_tet_forecast": float(forecast_tets),
            "tetra_volume_sum_mm3": mesh_volume, "signed_relative_volume_error": relative_volume_error,
            "absolute_volume_error_pct": abs(relative_volume_error) * 100.0,
            "mesh_msh_path": str(gmsh_mesh_file), "mesh_msh_sha256": gmsh_mesh_sha,
            "mesh_npz_path": str(npz), "mesh_npz_sha256": sha256(npz),
            "mesher": mesher,
            "native_centroid_state_counts": {key: int(state_counts.get(key, 0)) for key in ("IN", "OUT", "UNKNOWN")},
            "native_finite_stratified_interior_sample": {"count": int(len(samples)), "seed": INTERIOR_SAMPLE_SEED, "state_counts": {key: int(sample_counts.get(key, 0)) for key in ("IN", "OUT", "UNKNOWN")}, "method": "512 deterministic tet-stratified cells; Dirichlet(1,1,1,1) barycentric point in each selected tetrahedron; each point classified against the untouched canonical source B-rep"},
            "through_wall_normal_sampling": wall_dist,
            "actual_y0_section_png": section_png,
            "volume_mesh_quality_and_direction": quality,
            "boundary_quality": boundary,
            "checks": checks,
        },
        "interpretation": {
            "scope": "bounded local static geometry/mesh study only",
            "delta_note": "δ=1.0mm is a predeclared local neighborhood-fabric diagnostic only. Native material horizon, surface correction, and mechanics are not selected or validated.",
            "not_conflated": "tet edge length and the candidate diagnostic δ are reported separately; neither establishes the other",
            "no_fracture_claim": "No bonds, solver, crack, fracture plane, or R2 fragments were generated. Box cut surfaces are artificial ROI boundaries.",
            "physics_limit": "0.5mm vs 1.0mm geometry sampling is not spatial or temporal physical convergence. The local native normal thickness divided by h is a candidate sampling ratio only; t/h≈8 at h=0.5mm is not a universal physical sufficiency condition. Follow-up must hold material parameters and nonlocal horizon δ fixed while refining spatial h, then separately hold mesh fixed and halve dt for time convergence; also assess δ/h neighborhood support, surface correction, thin-wall bending response, and mesh-direction/remeshing dependence.",
            "continuum_fabric_reference": "the δ=1mm candidate fabric comparison uses uniform-in-volume Monte Carlo samples in the same-radius ball classified inside the untouched canonical STEP; this accounts for true physical wall truncation. The 1mm value is still exploratory, not a selected or calibrated material length.",
            "rotation_limit": "Rigidly rotating an unchanged mesh and loads tests coordinate covariance only; it does not test remeshing-direction dependence. To test direction sensitivity, keep the CAD and rest-coordinate material field/load fixed, change meshing/frame orientation, rebuild the point cloud and bonds, then compare crack path, strain energy and fragment mass distribution.",
            "orientation_claim_limit": "A crack constrained to tetra faces is not evidence of physical fracture. This static mesh study does not establish a tet-face fracture law or compare meshfree error.",
        },
        "timing": {"total_wall_seconds": float(time.perf_counter() - started), "mesh_generation_wall_seconds": float(elapsed_mesh_s), "max_trial_wall_seconds": MAX_WALL_SECONDS_PER_TRIAL},
    }
    (output_dir / "report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({"status": status, "mesh_size_mm": mesh_size_mm, "tetra_count": tet_count, "volume_error_pct": abs(relative_volume_error) * 100.0, "native_centroid_states": report["trial"]["native_centroid_state_counts"], "sample_states": report["trial"]["native_finite_stratified_interior_sample"]["state_counts"], "wall_normal_thickness": native_normal_report["normal_thickness_mm"], "mesh_seconds": elapsed_mesh_s, "report_json": str(output_dir / "report.json")}, indent=2))
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh-size-mm", type=float, choices=(0.5, 1.0))
    parser.add_argument("--preflight-only", action="store_true", help="measure exact fixed native CAD∩box volume and forecast both fixed h values without meshing")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.preflight_only:
        if args.mesh_size_mm is not None:
            parser.error("--preflight-only cannot be combined with --mesh-size-mm")
        report = native_roi_preflight(args.output_dir.resolve())
        return 0 if report["status"] == "PREFLIGHT_WITHIN_BOTH_FIXED_TET_BUDGETS" else 3
    if args.mesh_size_mm is None:
        parser.error("--mesh-size-mm is required unless --preflight-only is selected")
    report = run_trial(float(args.mesh_size_mm), args.output_dir.resolve())
    if report["status"].startswith("STOPPED_PREFLIGHT"):
        return 3
    return 0 if report["status"] not in ("FAIL_CELL_LIMIT_EXCEEDED_FIRST_ATTEMPT_PRESERVED", "FAIL_LOCAL_GEOMETRY_AND_MESH_CHECKS") else 2


if __name__ == "__main__":
    raise SystemExit(main())
