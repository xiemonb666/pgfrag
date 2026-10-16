"""Render the frozen native cohesive-fragment packet with VTK.

This creates a two-panel scientific preview from the saved t=50 us native
pose. The assembly view preserves the packet coordinates. The two individual
fragment views use separate cameras to expose each fracture surface; actors
receive no display translation, rotation, or scale.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_PACKET = ROOT / "pgfrag/data/continuation_20261015/native_cohesive_fragments_v4"
VISUALIZATION_ROOT = ROOT / "pgfrag/experiments/cohesive_vessel_20261015"
DEFAULT_OUTPUT = VISUALIZATION_ROOT / "visualization_v4/native_fragments_v4_preview_01"
WIDTH, HEIGHT = 2000, 1100
BACKGROUND_TOP = (0.075, 0.088, 0.102)
BACKGROUND_BOTTOM = (0.020, 0.026, 0.034)
LABEL_WHITE = (0.93, 0.93, 0.90)
LABEL_MUTED = (0.67, 0.72, 0.75)
LABEL_ACCENT = (0.98, 0.69, 0.36)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def strict_json_dump(path: Path, value: Any) -> None:
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def verify_packet(packet: Path) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any], dict[str, Any]]:
    manifest_path = packet / "manifest.json"
    manifest = read_json(manifest_path)
    if manifest.get("schema") != "pgfrag.native_cohesive_fragments.v1":
        raise ValueError(f"unexpected native fragment packet schema: {manifest.get('schema')!r}")
    verified: dict[str, dict[str, Any]] = {}
    for name, expected in manifest["files"].items():
        path = packet / name
        actual_bytes = path.stat().st_size
        actual_sha = sha256(path)
        if actual_bytes != int(expected["bytes"]) or actual_sha != expected["sha256"]:
            raise ValueError(f"native packet hash/size mismatch for {name}")
        verified[name] = {"bytes": actual_bytes, "sha256": actual_sha}

    fragments = read_json(packet / "fragments.json")
    state_file = np.load(packet / "state.npz", allow_pickle=False)
    state = {key: state_file[key] for key in state_file.files}
    forge = read_json(packet / "forge_manifest.json")
    if forge.get("material", {}).get("calibration_status") != "UNVALIDATED_PROTOTYPE":
        raise ValueError("the frozen appearance/material status is not the expected uncalibrated prototype")
    return manifest, fragments, state, forge | {"_verified_files": verified, "_manifest_sha256": sha256(manifest_path)}


def verify_counts_and_coverage(
    manifest: dict[str, Any], fragments: dict[str, Any], state: dict[str, np.ndarray]
) -> dict[str, Any]:
    counts = manifest["counts"]
    if len(fragments.get("components", [])) != 2 or counts.get("fragments") != 2:
        raise ValueError("native packet must contain exactly two mechanically connected fragments")
    nodes = state["current_nodes_m"]
    tets = state["tets"]
    triangles = state["triangles_native"]
    if nodes.shape != (counts["nodes"], 3):
        raise ValueError("saved native current-node count does not match packet manifest")
    if tets.shape != (counts["tetrahedra"], 4):
        raise ValueError("saved native T4 count does not match packet manifest")
    if triangles.shape != (len(fragments["boundary_faces"]), 3):
        raise ValueError("saved native surface triangles do not match fragment boundary face rows")
    if len(triangles) != 2208:
        raise ValueError(f"expected all 2,208 native boundary triangles, got {len(triangles)}")

    component_tet_rows = [int(row) for comp in fragments["components"] for row in comp["tet_rows"]]
    if len(component_tet_rows) != len(tets) or sorted(component_tet_rows) != list(range(len(tets))):
        raise ValueError("fragment components do not cover every native T4 row exactly once")
    component_face_counts = [len(comp["oriented_native_triangles"]) for comp in fragments["components"]]
    if sum(component_face_counts) != len(triangles):
        raise ValueError("component native triangles do not cover the saved boundary triangles")
    if not fragments.get("complete_tet_volume_coverage"):
        raise ValueError("fragment result does not assert complete T4 volume coverage")
    if not all(comp.get("closed_solid") for comp in fragments["components"]):
        raise ValueError("a native component boundary is not closed")
    if len(fragments.get("failed_interface_face_pairs", [])) != counts.get("paired_failed_facets"):
        raise ValueError("failed cohesive facet pair count does not match the packet manifest")

    regions, region_counts = np.unique(state["face_region_code"], return_counts=True)
    kinds, kind_counts = np.unique(state["face_kind"], return_counts=True)
    component_ids, face_component_counts = np.unique(state["face_component"], return_counts=True)
    return {
        "node_count": int(len(nodes)),
        "tetrahedron_count": int(len(tets)),
        "component_tetrahedron_counts": [int(len(comp["tet_rows"])) for comp in fragments["components"]],
        "tetrahedron_rows_covered_once": int(len(component_tet_rows)),
        "boundary_triangle_count": int(len(triangles)),
        "component_boundary_triangle_counts": [int(n) for n in component_face_counts],
        "face_region_counts": {str(int(key)): int(value) for key, value in zip(regions, region_counts, strict=True)},
        "face_kind_counts": {str(int(key)): int(value) for key, value in zip(kinds, kind_counts, strict=True)},
        "component_boundary_face_counts": {str(int(key)): int(value) for key, value in zip(component_ids, face_component_counts, strict=True)},
        "failed_cohesive_facet_count": int(len(fragments["failed_interface_face_pairs"])),
        "all_native_tets_covered_exactly_once": True,
        "all_native_boundary_triangles_rendered": True,
    }


def verify_area_pairing_winding_invariance() -> dict[str, Any]:
    """Check exact rest-area pairing under a different cohesive face winding."""
    sys.path.insert(0, str(ROOT / "pgfrag/src"))
    from pgfrag.cohesive_fragments import extract_cohesive_fragments

    rest = np.asarray([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0],
    ])
    tets = np.asarray([[0, 1, 2, 3], [4, 6, 5, 7]], dtype=np.int64)
    # Side B has a cyclically different winding/order; source ancestry resolves
    # correspondence and paired area is evaluated in that correspondence order.
    cohesive = np.asarray([[0, 1, 2, 5, 6, 4]], dtype=np.int64)
    source = np.asarray([0, 1, 2, 3, 0, 1, 2, 4], dtype=np.int64)
    result = extract_cohesive_fragments(
        rest,
        tets,
        cohesive,
        facet_damage_qp=[[1.0]],
        tractions_qp=[[[0.0, 0.0, 0.0]]],
        damage_tolerance=0.01,
        traction_tolerance=0.0,
        native_node_to_source=source,
    )
    pair = result["failed_interface_face_pairs"][0]
    if pair["side_a_rest_area"] != pair["side_b_rest_area"]:
        raise AssertionError("exact rest-area check changed under cohesive side winding")
    return {
        "status": "PASS",
        "method": "synthetic native T4 pair with exact duplicate rest nodes and cyclic side-B ordering",
        "source_node_ancestry_used": True,
        "side_a_rest_area_m2": pair["side_a_rest_area"],
        "side_b_rest_area_m2": pair["side_b_rest_area"],
        "areas_exactly_equal": True,
        "extractor_source_path": str(ROOT / "pgfrag/src/pgfrag/cohesive_fragments.py"),
        "extractor_source_sha256": sha256(ROOT / "pgfrag/src/pgfrag/cohesive_fragments.py"),
    }


def linear_to_srgb(rgb: np.ndarray) -> np.ndarray:
    values = np.asarray(rgb, dtype=np.float64)
    return np.where(values <= 0.0031308, 12.92 * values, 1.055 * np.power(values, 1.0 / 2.4) - 0.055)


def make_texture(forge: dict[str, Any], rest_z_m: np.ndarray, pixel_count: int = 1024) -> tuple[vtk.vtkTexture, dict[str, Any]]:
    sys.path.insert(0, str(ROOT / "vesselforge/src"))
    from vesselforge.spherical_bowl import evaluate_bowl_paint_at_xyz

    recipe = forge["ornament"]["recipe"]
    z_min_m = float(np.min(rest_z_m))
    z_max_m = float(np.max(rest_z_m))
    samples_u = (np.arange(pixel_count, dtype=np.float64) + 0.5) / pixel_count
    sample_z_m = z_min_m + samples_u * (z_max_m - z_min_m)
    xyz_mm = np.zeros((pixel_count, 3), dtype=np.float64)
    xyz_mm[:, 2] = sample_z_m * 1000.0
    colors_linear = evaluate_bowl_paint_at_xyz(xyz_mm, recipe)
    colors_srgb = np.clip(linear_to_srgb(colors_linear), 0.0, 1.0)
    pixels = np.rint(colors_srgb * 255.0).astype(np.uint8)

    image = vtk.vtkImageData()
    image.SetDimensions(pixel_count, 1, 1)
    image.AllocateScalars(vtk.VTK_UNSIGNED_CHAR, 3)
    image.GetPointData().SetScalars(numpy_to_vtk(pixels, deep=True, array_type=vtk.VTK_UNSIGNED_CHAR))
    texture = vtk.vtkTexture()
    texture.SetInputData(image)
    texture.SetInterpolate(False)
    texture.RepeatOff()
    texture.EdgeClampOn()

    return texture, {
        "method": "continuous recipe evaluated at 1,024 rest-Z texel centers; linear albedo converted to sRGB bytes for display",
        "texture_dimensions": [pixel_count, 1],
        "rest_z_source_range_m": [z_min_m, z_max_m],
        "region_codes_textured": [1],
        "recipe": recipe,
    }


def make_polydata(points: np.ndarray, rest_z_m: np.ndarray, triangles: np.ndarray) -> vtk.vtkPolyData:
    polydata = vtk.vtkPolyData()
    vtk_points = vtk.vtkPoints()
    vtk_points.SetData(numpy_to_vtk(np.asarray(points, dtype=np.float64), deep=True))
    polydata.SetPoints(vtk_points)

    flat = np.empty((len(triangles), 4), dtype=np.int64)
    flat[:, 0] = 3
    flat[:, 1:] = np.asarray(triangles, dtype=np.int64)
    cells = vtk.vtkCellArray()
    cells.SetCells(len(triangles), numpy_to_vtkIdTypeArray(flat.reshape(-1), deep=True))
    polydata.SetPolys(cells)

    zmin, zmax = float(np.min(rest_z_m)), float(np.max(rest_z_m))
    u = (np.asarray(rest_z_m, dtype=np.float64) - zmin) / (zmax - zmin)
    tcoords_np = np.column_stack((u, np.full(len(u), 0.5, dtype=np.float64)))
    tcoords = numpy_to_vtk(tcoords_np, deep=True)
    tcoords.SetNumberOfComponents(2)
    tcoords.SetName("rest_z_texture_coordinates")
    polydata.GetPointData().SetTCoords(tcoords)
    return polydata


def make_surface_actor(
    points_m: np.ndarray,
    rest_z_m: np.ndarray,
    triangles: np.ndarray,
    linear_color: list[float] | tuple[float, ...],
    texture: vtk.vtkTexture | None = None,
) -> vtk.vtkActor:
    polydata = make_polydata(points_m, rest_z_m, triangles)
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputData(polydata)
    normals.ComputePointNormalsOn()
    normals.ComputeCellNormalsOff()
    normals.ConsistencyOn()
    normals.AutoOrientNormalsOff()
    normals.SplittingOff()
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(normals.GetOutputPort())
    mapper.ScalarVisibilityOff()
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    if texture is not None:
        actor.SetTexture(texture)
        actor.GetProperty().SetColor(1.0, 1.0, 1.0)
    else:
        display_rgb = linear_to_srgb(np.asarray(linear_color, dtype=np.float64))
        actor.GetProperty().SetColor(*map(float, display_rgb))
    prop = actor.GetProperty()
    prop.SetInterpolationToPhong()
    prop.SetAmbient(0.22)
    prop.SetDiffuse(0.78)
    prop.SetSpecular(0.10)
    prop.SetSpecularPower(22.0)
    prop.SetEdgeVisibility(False)
    return actor


def add_label(
    renderer: vtk.vtkRenderer,
    value: str,
    x: float,
    y: float,
    size: int,
    color: tuple[float, float, float] = LABEL_WHITE,
    *,
    bold: bool = False,
) -> vtk.vtkTextActor:
    actor = vtk.vtkTextActor()
    actor.SetInput(value)
    position = actor.GetPositionCoordinate()
    position.SetCoordinateSystemToNormalizedViewport()
    actor.SetPosition(float(x), float(y))
    text = actor.GetTextProperty()
    text.SetFontFamilyToArial()
    text.SetFontSize(int(size))
    text.SetColor(*color)
    text.SetBold(bool(bold))
    text.SetShadow(False)
    renderer.AddActor2D(actor)
    return actor


def add_camera(
    renderer: vtk.vtkRenderer,
    focal_point_m: np.ndarray,
    direction: np.ndarray,
    view_up: np.ndarray,
    parallel_scale_m: float,
    *,
    distance_m: float = 0.08,
) -> dict[str, Any]:
    direction = np.asarray(direction, dtype=np.float64)
    direction /= np.linalg.norm(direction)
    position = np.asarray(focal_point_m, dtype=np.float64) + direction * float(distance_m)
    camera = vtk.vtkCamera()
    camera.ParallelProjectionOn()
    camera.SetFocalPoint(*map(float, focal_point_m))
    camera.SetPosition(*map(float, position))
    camera.SetViewUp(*map(float, view_up))
    camera.SetParallelScale(float(parallel_scale_m))
    renderer.SetActiveCamera(camera)
    renderer.ResetCameraClippingRange()
    return {
        "projection": "parallel",
        "focal_point_m": np.asarray(focal_point_m, dtype=np.float64).tolist(),
        "position_m": position.tolist(),
        "view_direction_unit": direction.tolist(),
        "view_up": np.asarray(view_up, dtype=np.float64).tolist(),
        "parallel_scale_m": float(parallel_scale_m),
        "camera_distance_m": float(distance_m),
    }


def make_renderer() -> vtk.vtkRenderer:
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(*BACKGROUND_BOTTOM)
    renderer.SetBackground2(*BACKGROUND_TOP)
    renderer.GradientBackgroundOn()
    renderer.SetAmbient(0.26, 0.26, 0.26)
    renderer.SetAutomaticLightCreation(True)
    renderer.SetInteractive(False)
    return renderer


def add_component_actors(
    renderer: vtk.vtkRenderer,
    component_index: int,
    state: dict[str, np.ndarray],
    texture: vtk.vtkTexture,
    base_linear_rgb: list[float],
) -> None:
    component_mask = state["face_component"] == component_index
    regions = state["face_region_code"]
    triangles = state["triangles_native"]
    for region_code in (1, 2, 3, 255):
        selected = triangles[component_mask & (regions == region_code)]
        if len(selected) == 0:
            continue
        actor = make_surface_actor(
            state["current_nodes_m"],
            state["rest_nodes_m"][:, 2],
            selected,
            base_linear_rgb,
            texture if region_code == 1 else None,
        )
        renderer.AddActor(actor)


def render_preview(
    packet: Path,
    output: Path,
    manifest: dict[str, Any],
    fragments: dict[str, Any],
    state: dict[str, np.ndarray],
    forge: dict[str, Any],
    source_counts: dict[str, Any],
    winding_check: dict[str, Any],
) -> dict[str, Any]:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing preview output: {output}")
    if not output.resolve().is_relative_to(VISUALIZATION_ROOT.resolve()):
        raise ValueError(f"preview output must stay under {VISUALIZATION_ROOT}")

    current = state["current_nodes_m"]
    component_centers_x = [
        float(current[np.asarray(comp["native_node_ids"], dtype=np.int64), 0].mean())
        for comp in fragments["components"]
    ]
    left_comp = min(range(2), key=lambda i: component_centers_x[i])
    right_comp = max(range(2), key=lambda i: component_centers_x[i])
    left_nodes = np.asarray(fragments["components"][left_comp]["native_node_ids"], dtype=np.int64)
    right_nodes = np.asarray(fragments["components"][right_comp]["native_node_ids"], dtype=np.int64)
    measured_gap_m = float(current[right_nodes, 0].min() - current[left_nodes, 0].max())
    declared_gap_m = float(read_json(packet / "separation.json")["final_global_x_gap_m"])
    if not np.isclose(measured_gap_m, declared_gap_m, rtol=0.0, atol=1.0e-12):
        raise ValueError(f"current component gap disagrees with packet separation receipt: {measured_gap_m} vs {declared_gap_m}")

    output.parent.mkdir(parents=True, exist_ok=True)
    output.mkdir(parents=False)

    paint_texture, texture_receipt = make_texture(forge, state["rest_nodes_m"][:, 2])
    base_linear_rgb = list(forge["material"]["appearance"]["base_albedo_linear_rgb"])
    render_window = vtk.vtkRenderWindow()
    render_window.SetSize(WIDTH, HEIGHT)
    render_window.SetOffScreenRendering(1)
    render_window.SetMultiSamples(8)
    render_window.SetNumberOfLayers(2)
    render_window.SetWindowName("Native cohesive fragment preview")

    assembly_renderer = make_renderer()
    assembly_renderer.SetViewport(0.025, 0.035, 0.555, 0.885)
    for component_index in range(2):
        add_component_actors(assembly_renderer, component_index, state, paint_texture, base_linear_rgb)
    render_window.AddRenderer(assembly_renderer)
    world_bounds = np.asarray(state["current_nodes_m"], dtype=np.float64)
    assembly_focal = 0.5 * (world_bounds.min(axis=0) + world_bounds.max(axis=0))
    assembly_camera = add_camera(
        assembly_renderer,
        assembly_focal,
        np.asarray([0.0, -1.0, 0.10]),
        np.asarray([0.0, 0.0, 1.0]),
        0.0145,
        distance_m=0.08,
    )
    add_label(assembly_renderer, "ASSEMBLY  |  NATIVE SAVED POSE", 0.035, 0.955, 20, LABEL_ACCENT, bold=True)
    add_label(assembly_renderer, f"{measured_gap_m * 1.0e6:.3f} µm x-gap  ·  coordinates unchanged", 0.035, 0.915, 15, LABEL_WHITE)

    comp_renderers: list[vtk.vtkRenderer] = []
    fragment_cameras: list[dict[str, Any]] = []
    viewport_rows = [
        (0.590, 0.475, 0.980, 0.875),
        (0.590, 0.055, 0.980, 0.455),
    ]
    centers: list[np.ndarray] = []
    for comp_index, comp in enumerate(fragments["components"]):
        ids = np.asarray(comp["native_node_ids"], dtype=np.int64)
        center = state["current_nodes_m"][ids].mean(axis=0)
        centers.append(center)
        renderer = make_renderer()
        renderer.SetViewport(*viewport_rows[comp_index])
        add_component_actors(renderer, comp_index, state, paint_texture, base_linear_rgb)
        crack_side_sign = -1.0 if float(center[0]) > 0.0 else 1.0
        direction = np.asarray([crack_side_sign, -0.38, 0.20], dtype=np.float64)
        camera = add_camera(
            renderer,
            center,
            direction,
            np.asarray([0.0, 0.0, 1.0]),
            0.0138,
            distance_m=0.065,
        )
        comp_renderers.append(renderer)
        fragment_cameras.append(camera)
        n_tets = len(comp["tet_rows"])
        n_faces = len(comp["oriented_native_triangles"])
        add_label(renderer, f"FRAGMENT {comp_index}  ·  {n_tets:,} T4  ·  {n_faces:,} native faces", 0.045, 0.925, 16, LABEL_ACCENT, bold=True)
        add_label(renderer, "independent camera  ·  no simulated or actor motion", 0.045, 0.875, 13, LABEL_MUTED)
    for renderer in comp_renderers:
        render_window.AddRenderer(renderer)

    header = vtk.vtkRenderer()
    header.SetLayer(1)
    header.SetViewport(0.0, 0.900, 1.0, 1.0)
    header.SetBackground(0.025, 0.031, 0.040)
    header.SetBackgroundAlpha(1.0)
    header.SetInteractive(False)
    render_window.AddRenderer(header)
    title = add_label(header, "WEAK SEAM  /  UNCALIBRATED", 0.035, 0.55, 29, LABEL_WHITE, bold=True)
    add_label(header, "NATIVE t = 50 µs  ·  frozen same-frame state  ·  two independent fracture-face cameras at right", 0.035, 0.15, 16, LABEL_MUTED)

    # Actors inherit exact packet current coordinates and default identity
    # transforms. The individual panes differ only by camera and focal point.
    render_window.Render()
    capture = vtk.vtkWindowToImageFilter()
    capture.SetInput(render_window)
    capture.SetInputBufferTypeToRGB()
    capture.ReadFrontBufferOff()
    capture.Update()
    png_path = output / "native_t50us_weak_seam_two_panel.png"
    writer = vtk.vtkPNGWriter()
    writer.SetFileName(str(png_path))
    writer.SetInputConnection(capture.GetOutputPort())
    writer.Write()
    if not png_path.exists() or png_path.stat().st_size < 50_000:
        raise RuntimeError("VTK did not produce the expected full-size PNG preview")

    image_receipt = {
        "file": png_path.name,
        "sha256": sha256(png_path),
        "bytes": png_path.stat().st_size,
        "dimensions_px": [WIDTH, HEIGHT],
        "rendered_cell_count": source_counts["boundary_triangle_count"],
        "rendered_triangle_count_by_component": source_counts["component_boundary_triangle_counts"],
    }
    receipt = {
        "schema": "pgfrag.native_fragment_visualization_receipt.v1",
        "source_packet": {
            "directory": str(packet.resolve()),
            "manifest_sha256": sha256(packet / "manifest.json"),
            "manifest_schema": manifest["schema"],
            "verified_packet_files": manifest["files"],
            "all_listed_packet_files_verified": True,
            "forge_manifest_sha256": sha256(packet / "forge_manifest.json"),
            "fragments_json_sha256": sha256(packet / "fragments.json"),
            "state_npz_sha256": sha256(packet / "state.npz"),
        },
        "renderer_code": {
            "path": str(Path(__file__).resolve()),
            "sha256": sha256(Path(__file__).resolve()),
            "vtk_version": vtk.vtkVersion.GetVTKVersion(),
            "numpy_version": np.__version__,
            "python_executable": sys.executable,
        },
        "extractor_winding_area_check": winding_check,
        "native_pose": {
            "frame_label": "native packet state at case time.total_s",
            "time_s": float(read_json(packet / "case.json")["time"]["total_s"]),
            "time_us": float(read_json(packet / "case.json")["time"]["total_s"]) * 1.0e6,
            "current_nodes_source": "state.npz/current_nodes_m",
            "model_coordinate_units": "metres",
            "display_unit_scale": 1.0,
            "component_x_gap_m": measured_gap_m,
            "component_x_gap_um": measured_gap_m * 1.0e6,
        },
        "source_counts": source_counts,
        "paint": texture_receipt,
        "display_transform": {
            "actor_translation_m": [0.0, 0.0, 0.0],
            "actor_rotation_degrees": [0.0, 0.0, 0.0],
            "actor_scale": [1.0, 1.0, 1.0],
            "assembly_geometry": "unchanged current_nodes_m from the packet",
            "individual_fragment_geometry": "same unchanged current_nodes_m; each renderer filters to one component",
            "individual_fragment_camera_reorientation_only": True,
            "simulation_state_modified": False,
            "native_tet_mesh_drawn_as": "all saved native boundary triangles; no downsampling or caps",
        },
        "cameras": {
            "window_dimensions_px": [WIDTH, HEIGHT],
            "assembly_viewport_normalized": [0.025, 0.035, 0.555, 0.885],
            "assembly": assembly_camera,
            "fragment_0_viewport_normalized": list(viewport_rows[0]),
            "fragment_0_camera": fragment_cameras[0],
            "fragment_1_viewport_normalized": list(viewport_rows[1]),
            "fragment_1_camera": fragment_cameras[1],
            "background_linear_rgb_top": list(BACKGROUND_TOP),
            "background_linear_rgb_bottom": list(BACKGROUND_BOTTOM),
        },
        "labels": [
            "WEAK SEAM / UNCALIBRATED",
            "NATIVE t = 50 µs",
            "ASSEMBLY / NATIVE SAVED POSE",
            "independent camera / no simulated or actor motion",
        ],
        "image": image_receipt,
        "qualification": {
            "status": manifest["qualification"]["status"],
            "material_calibration_status": manifest["qualification"]["material_calibration_status"],
            "physical_control_status": manifest["qualification"].get("physical_control_status", "NOT_REPORTED"),
            "physical_pass": manifest["qualification"]["physical_pass"],
            "typed_R2_qualification": manifest["qualification"]["typed_R2_qualification"],
            "pretrain_eligible": manifest["qualification"]["pretrain_eligible"],
        },
    }
    strict_json_dump(output / "render_receipt.json", receipt)
    render_window.Finalize()
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--packet", type=Path, default=DEFAULT_PACKET)
    parser.add_argument("--out-dir", "--out", dest="out_dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    packet = args.packet.resolve()
    output = args.out_dir.resolve()
    manifest, fragments, state, forge_bundle = verify_packet(packet)
    verified_files = forge_bundle.pop("_verified_files")
    manifest_sha = forge_bundle.pop("_manifest_sha256")
    source_counts = verify_counts_and_coverage(manifest, fragments, state)
    winding_check = verify_area_pairing_winding_invariance()
    forge_bundle["_verified_files"] = verified_files
    forge_bundle["_manifest_sha256"] = manifest_sha
    receipt = render_preview(packet, output, manifest, fragments, state, forge_bundle, source_counts, winding_check)
    print(json.dumps({
        "output": str(output),
        "image": str(output / receipt["image"]["file"]),
        "image_sha256": receipt["image"]["sha256"],
        "receipt": str(output / "render_receipt.json"),
        "source_counts": source_counts,
        "component_x_gap_um": receipt["native_pose"]["component_x_gap_um"],
    }, indent=2))


if __name__ == "__main__":
    main()
