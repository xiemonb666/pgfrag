#!/usr/bin/env python
"""One bounded OCC-native clipped Cartesian volume-cell candidate experiment.

Static geometry only: this does not create an accepted PG job, solver input,
native bond graph, realistic fragment model, or physical fracture result.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import logging
import math
import os
import sys
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")

import cadquery as cq
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from OCP.BRepAlgoAPI import BRepAlgoAPI_Section
from OCP.BRepBndLib import BRepBndLib
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.BRepGProp import BRepGProp
from OCP.Bnd import Bnd_Box
from OCP.BOPAlgo import BOPAlgo_BOP, BOPAlgo_Operation
from OCP.GProp import GProp_GProps
from OCP.TopAbs import TopAbs_COMPOUND, TopAbs_IN, TopAbs_ON, TopAbs_OUT
from OCP.TopTools import TopTools_ListOfShape
from OCP.gp import gp_Ax3, gp_Dir, gp_Pln, gp_Pnt


ROOT = Path(__file__).resolve().parents[2]
MANIFEST = (
    ROOT / "relicweave-synlab" / "tmp" / "review_20261003" / "implementation"
    / "vf" / "artifacts" / "four_bands_brep_semantics1_final" / "manifest.json"
)
OUTPUT_ROOT = (
    ROOT / "relicweave-synlab" / "tmp" / "review_20261004" / "implementation"
    / "cad_clipped_cells"
)
OUTPUT_DIR = OUTPUT_ROOT / "fixed8mm_singlethread_case03_native_empty_compound_policy"
SPACING_MM = 8.0
MAX_BOXES = 512
BOOLEAN_WALL_LIMIT_SECONDS = 300.0
VOLUME_ERROR_LIMIT = 0.01
CLASSIFIER_TOL_MM = 1.0e-7
EXPECTED_STEP_SHA256 = "4623bd373aeed0acfacfbf6cc0dda65952339561779b43429c59d61b8b5a1c30"


class ExperimentFailure(RuntimeError):
    pass


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def configure_log(output_dir: Path) -> logging.Logger:
    logger = logging.getLogger("cad_clipped_cells")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    formatter = logging.Formatter("%(asctime)sZ %(levelname)s %(message)s")
    formatter.converter = time.gmtime
    for handler in (
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(output_dir / "run.log", mode="w", encoding="utf-8"),
    ):
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    return logger


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def exact_volume_and_centroid(shape: cq.Shape) -> tuple[float, tuple[float, float, float]]:
    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape.wrapped, props, False, False, False)
    center = props.CentreOfMass()
    return float(props.Mass()), (float(center.X()), float(center.Y()), float(center.Z()))


def exact_bbox(shape: cq.Shape) -> tuple[float, float, float, float, float, float]:
    box = Bnd_Box()
    BRepBndLib.AddOptimal_s(shape.wrapped, box, False, False)
    raw = tuple(float(v) for v in box.Get())
    if len(raw) != 6 or not all(math.isfinite(v) for v in raw):
        raise ExperimentFailure(f"invalid native BRep bbox: {raw}")
    # Bnd_Box.Get uses xmin,ymin,zmin,xmax,ymax,zmax; normalize to
    # xmin,xmax,ymin,ymax,zmin,zmax for grid and report records.
    xmin, ymin, zmin, xmax, ymax, zmax = raw
    return xmin, xmax, ymin, ymax, zmin, zmax


def axis_edges(lo: float, hi: float, step: float) -> list[float]:
    count = int(math.ceil((hi - lo) / step))
    edges = [lo + i * step for i in range(count)] + [hi]
    if len(edges) < 2 or any(b <= a for a, b in zip(edges, edges[1:])):
        raise ExperimentFailure(f"invalid grid axis [{lo}, {hi}] at {step} mm")
    edges[-1] = hi
    return edges


def make_box_records(
    bbox: tuple[float, float, float, float, float, float], spacing: float
) -> tuple[list[dict[str, Any]], tuple[list[float], list[float], list[float]]]:
    xmin, xmax, ymin, ymax, zmin, zmax = bbox
    xs, ys, zs = (
        axis_edges(xmin, xmax, spacing),
        axis_edges(ymin, ymax, spacing),
        axis_edges(zmin, zmax, spacing),
    )
    boxes: list[dict[str, Any]] = []
    for i in range(len(xs) - 1):
        for j in range(len(ys) - 1):
            for k in range(len(zs) - 1):
                bounds = (xs[i], xs[i + 1], ys[j], ys[j + 1], zs[k], zs[k + 1])
                boxes.append({
                    "cell_id": len(boxes),
                    "ijk": [i, j, k],
                    "bounds_mm": list(bounds),
                    "aabb_volume_mm3": (
                        (bounds[1] - bounds[0])
                        * (bounds[3] - bounds[2])
                        * (bounds[5] - bounds[4])
                    ),
                })
    return boxes, (xs, ys, zs)


def top_tools(items: list[cq.Shape]) -> TopTools_ListOfShape:
    result = TopTools_ListOfShape()
    for shape in items:
        result.Append(shape.wrapped)
    return result


def topology_counts(shape: cq.Shape) -> dict[str, int]:
    return {
        "solids": len(shape.Solids()),
        "shells": len(shape.Shells()),
        "faces": len(shape.Faces()),
        "wires": len(shape.Wires()),
        "edges": len(shape.Edges()),
        "vertices": len(shape.Vertices()),
    }


def native_common(
    cad_solid: cq.Solid, box: cq.Solid
) -> tuple[cq.Shape, bool, bool, bool, bool, str]:
    # BOPAlgo_BOP is the native OCC Boolean Common builder. This OCP binding
    # exposes HasErrors/GetReport here (BRepAlgoAPI_Common does not bind them).
    op = BOPAlgo_BOP()
    op.SetArguments(top_tools([cad_solid]))
    op.SetTools(top_tools([box]))
    op.SetOperation(BOPAlgo_Operation.BOPAlgo_COMMON)
    op.SetRunParallel(False)
    op.Perform()
    has_errors = bool(op.HasErrors())
    has_warnings = bool(op.HasWarnings())
    report_buffer = io.BytesIO()
    op.GetReport().DumpJson(report_buffer)
    report_json = report_buffer.getvalue().decode("utf-8", errors="replace")
    native_result = op.Shape()
    result = cq.Shape.cast(native_result)
    # BOPAlgo_BOP has no IsDone() binding. An algorithm result is usable only
    # when it reports no errors; warnings are retained in the ledger/report.
    return result, not has_errors, bool(native_result.IsNull()), has_errors, has_warnings, report_json


def json_vec3(values: tuple[float, float, float]) -> list[float | None]:
    return [float(v) if math.isfinite(v) else None for v in values]


def state_name(state: Any) -> str:
    if state == TopAbs_IN:
        return "IN"
    if state == TopAbs_OUT:
        return "OUT"
    if state == TopAbs_ON:
        return "ON"
    return f"UNKNOWN({state})"


def classify(
    classifier: BRepClass3d_SolidClassifier, point: tuple[float, float, float]
) -> str:
    classifier.Perform(gp_Pnt(*point), CLASSIFIER_TOL_MM)
    return state_name(classifier.State())


def save_shape_brep(shape: cq.Shape, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    shape.exportBrep(str(path))


def exact_section_points(shape: cq.Shape, plane_face: cq.Face) -> list[np.ndarray]:
    op = BRepAlgoAPI_Section()
    op.SetArguments(top_tools([shape]))
    op.SetTools(top_tools([plane_face]))
    op.SetRunParallel(False)
    op.Build()
    if not op.IsDone():
        raise ExperimentFailure("OCC BRepAlgoAPI_Section did not complete")
    section_shape = cq.Shape.cast(op.Shape())
    paths: list[np.ndarray] = []
    for edge in section_shape.Edges():
        # Curve discretization is for the figure only. It is not candidate data.
        samples, _parameters = edge.sample(32)
        coords = np.asarray([[p.x, p.z] for p in samples], dtype=np.float64)
        if len(coords) >= 2 and np.isfinite(coords).all():
            paths.append(coords)
    return paths


def make_figure(
    output_path: Path,
    cad_solid: cq.Solid,
    box_records: list[dict[str, Any]],
    component_records: list[dict[str, Any]],
    component_shapes: list[cq.Solid],
    points_mm: np.ndarray,
    bbox: tuple[float, float, float, float, float, float],
    logger: logging.Logger,
) -> dict[str, Any]:
    xmin, xmax, ymin, ymax, zmin, zmax = bbox
    pad = 5.0
    # Set plane-local U=world +X and V=world +Z explicitly. With a -Y normal,
    # V = normal × U = +Z, so these UV bounds are the intended world X/Z bounds.
    plane_axes = gp_Ax3(
        gp_Pnt(0.0, 0.0, 0.0), gp_Dir(0.0, -1.0, 0.0), gp_Dir(1.0, 0.0, 0.0)
    )
    plane = BRepBuilderAPI_MakeFace(
        gp_Pln(plane_axes),
        xmin - pad, xmax + pad, zmin - pad, zmax + pad,
    ).Face()
    plane_face = cq.Face(plane)
    cad_paths = exact_section_points(cad_solid, plane_face)
    fig, ax = plt.subplots(figsize=(9.2, 8.0), constrained_layout=True)
    for i, path in enumerate(cad_paths):
        ax.plot(
            path[:, 0], path[:, 1], color="#222222", linewidth=1.7,
            label="canonical CAD section at y=0" if i == 0 else None, zorder=4,
        )
    component_paths = 0
    for record, shape in zip(component_records, component_shapes):
        cell = box_records[record["cell_id"]]
        if cell["bounds_mm"][2] > 0.0 or cell["bounds_mm"][3] < 0.0:
            continue
        for path in exact_section_points(shape, plane_face):
            ax.plot(path[:, 0], path[:, 1], color="#718096", linewidth=0.55, alpha=0.6)
            component_paths += 1
    if len(points_mm):
        point_states = [c["center_state_in_canonical_cad"] for c in component_records]
        in_mask = np.asarray([state == "IN" for state in point_states], dtype=bool)
        out_mask = np.asarray([state == "OUT" for state in point_states], dtype=bool)
        if in_mask.any():
            ax.scatter(
                points_mm[in_mask, 0], points_mm[in_mask, 2], s=12, color="#176b87",
                edgecolors="white", linewidths=0.25,
                label="all native COMs projected to XZ; IN in source CAD",
                zorder=5,
            )
        if out_mask.any():
            ax.scatter(
                points_mm[out_mask, 0], points_mm[out_mask, 2], s=48, marker="x",
                color="#c53030", linewidths=1.8,
                label="native COM classified OUT (candidate failure)", zorder=6,
            )
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(xmin - pad, xmax + pad)
    ax.set_ylim(zmin - pad, zmax + pad)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("z (mm)")
    ax.set_title("OCC-clipped 8 mm cells: exact y=0 sections and volume quadrature")
    ax.grid(False)
    ax.legend(loc="upper right", fontsize=8, framealpha=0.94)
    fig.savefig(output_path, dpi=190)
    plt.close(fig)
    logger.info(
        "PNG wrote %s; source OCC curves=%d clipped OCC curves=%d",
        output_path.name, len(cad_paths), component_paths,
    )
    cad_samples = np.concatenate(cad_paths, axis=0) if cad_paths else np.empty((0, 2))
    if len(cad_samples):
        section_bounds = [
            float(np.min(cad_samples[:, 0])), float(np.max(cad_samples[:, 0])),
            float(np.min(cad_samples[:, 1])), float(np.max(cad_samples[:, 1])),
        ]
    else:
        section_bounds = [None, None, None, None]
    body_wall_bands = {}
    for z_center in (8.0, 16.0, 24.0, 32.0, 40.0, 48.0, 52.0):
        nearby = cad_samples[np.abs(cad_samples[:, 1] - z_center) <= 1.0]
        body_wall_bands[str(int(z_center))] = {
            "negative_x_wall_sampled": bool(len(nearby) and np.any(nearby[:, 0] < -1.0)),
            "positive_x_wall_sampled": bool(len(nearby) and np.any(nearby[:, 0] > 1.0)),
            "sample_count": int(len(nearby)),
        }
    return {
        "plane_world_y_mm": 0.0,
        "plane_uv_world_axes": {"u": "+X", "v": "+Z", "normal": "-Y"},
        "source_cad_section_curve_count": len(cad_paths),
        "clipped_component_section_curve_count": component_paths,
        "source_cad_section_bounds_xmin_xmax_zmin_zmax_mm": section_bounds,
        "expected_source_bbox_xmin_xmax_zmin_zmax_mm": [xmin, xmax, zmin, zmax],
        "source_section_spans_z0_to_z56": bool(
            len(cad_samples)
            and abs(section_bounds[2] - zmin) <= 0.1
            and abs(section_bounds[3] - zmax) <= 0.1
        ),
        "source_section_has_left_and_right_x": bool(
            len(cad_samples)
            and np.any(cad_samples[:, 0] < -1.0)
            and np.any(cad_samples[:, 0] > 1.0)
        ),
        "body_wall_samples_by_z_mm": body_wall_bands,
        "native_com_markers_are_xz_projections_not_y0_membership": True,
        "figure_png_sha256": sha256_file(output_path),
    }


def cavity_probes() -> list[tuple[float, float, float]]:
    probes: list[tuple[float, float, float]] = []
    for z in (8.0, 16.0, 24.0, 32.0, 40.0, 48.0, 52.0):
        probes.append((0.0, 0.0, z))
        for radius in (4.0, 8.0):
            for angle in (0.0, 90.0, 180.0, 270.0):
                theta = math.radians(angle)
                probes.append((radius * math.cos(theta), radius * math.sin(theta), z))
    return probes


def write_json(path: Path, payload: Any) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def write_report_md(path: Path, report: dict[str, Any]) -> None:
    volume_error_pct = report["geometry"].get("relative_volume_error_pct")
    volume_error_display = (
        "NOT_EVALUATED_INCOMPLETE_PARTITION"
        if volume_error_pct is None
        else f"{volume_error_pct:.6f}% (gate {report['policy']['max_relative_volume_error_pct']:.3f}%)"
    )
    mass_error_pct = report["material"].get("relative_mass_error_pct")
    mass_error_display = (
        "NOT_EVALUATED_INCOMPLETE_PARTITION"
        if mass_error_pct is None else f"{mass_error_pct:.6f}%"
    )
    lines = [
        "# OCC-native clipped Cartesian volume-cell spike",
        "",
        f"Status: {report['status']}",
        f"Run state: {report.get('run_state', 'COMPLETE')}",
        "",
        "Static geometry candidate only. This is not an accepted PG job, solver "
        "input, native bond graph, realistic fragment result, or R2 result.",
        "",
        "## Source and fixed policy",
        "",
        f"- Canonical STEP: {report['source']['step_path']}",
        f"- SHA-256: {report['source']['step_sha256']} (manifest verified)",
        f"- Frame and length unit: {report['source']['frame']}, {report['source']['length_unit']}",
        f"- Fixed grid: {report['policy']['spacing_mm']} mm; OCC Common single-thread; "
        f"box cap {report['policy']['max_candidate_boxes']}",
        f"- Processed boxes: {report['counts']['processed_boxes']} / "
        f"{report['counts']['candidate_boxes']}",
        f"- Exact null Common results: {report['counts']['exact_null_common']}",
        f"- Total 3D-material empty cells (null or legal empty compound): "
        f"{report['counts']['empty_3d_material_results']}",
        f"- Native valid empty-material compound results: "
        f"{report['counts'].get('native_empty_compound', 0)}",
        f"- Positive connected native solids retained: {report['counts']['positive_components']}",
        f"- Non-null Common without positive solids: {report['counts']['nonvolume_nonempty_commons']}",
        f"- Exceptions and geometry gate failures: {report['counts']['exceptions']}",
        f"- Boolean time: {report['timing']['boolean_seconds']:.3f} s; whole run: "
        f"{report['timing']['whole_run_seconds']:.3f} s "
        f"(both limited to {report['policy']['boolean_wall_limit_seconds']:.0f} s)",
        "",
        "## Geometry and volume quadrature",
        "",
        f"- Exact CAD B-rep volume: {report['geometry']['cad_volume_mm3']:.9f} mm3",
        f"- Processed clipped-component volume subtotal: "
        f"{report['geometry']['processed_component_volume_mm3']:.9f} mm3",
        f"- Retained subtotal / full CAD volume: "
        f"{report['geometry']['processed_component_volume_fraction_of_cad']:.6%}",
        f"- Full-partition relative volume error: {volume_error_display}",
        f"- CAD B-rep bbox in mm: {report['geometry']['bbox_mm']}",
        f"- Finite positive components: {report['checks']['finite_positive_components']} / "
        f"{report['counts']['positive_components']}",
        f"- Native COM classifications: {report['checks']['centroid_state_counts']}",
        f"- Unique COMs at 1e-9 mm: {report['checks']['unique_centroid_count']} / "
        f"{report['counts']['positive_components']}",
        "- Each retained connected solid contributes one native center of mass and its "
        "exact positive B-rep volume. Separate solids inside one box remain separate.",
        "",
        "## Cavity and mass",
        "",
        f"- Cavity probes: {report['checks']['cavity_probe_state_counts']} over "
        f"{report['checks']['cavity_probe_scope']}",
        f"- Manifest rho {report['material']['rho_kg_m3']} kg/m3; processed subtotal mass: "
        f"{report['material']['processed_partial_mass_kg']:.12f} kg",
        f"- Full-partition mass error: {mass_error_display}",
        f"- Material calibration: {report['material']['calibration_status']}",
        "",
        "## Gates and limits",
        "",
    ]
    for name, passed in report["gates"].items():
        lines.append(f"- {name}: {'PASS' if passed else 'FAIL'}")
    figure_diagnostics = report.get("figure_diagnostics") or {}
    if figure_diagnostics:
        bounds = figure_diagnostics.get(
            "source_cad_section_bounds_xmin_xmax_zmin_zmax_mm", []
        )
        lines.extend([
            "",
            "## Exact section figure diagnostic",
            "",
            f"- Corrected PNG SHA-256: {figure_diagnostics.get('figure_png_sha256')}",
            f"- Source OCC y=0 section curves: {figure_diagnostics.get('source_cad_section_curve_count')}",
            f"- Sampled section bounds xmin,xmax,zmin,zmax (mm): {bounds}",
            f"- Both body walls sampled across z bands: "
            f"{figure_diagnostics.get('body_wall_samples_by_z_mm')}",
            f"- COM markers are XZ projections, not points constrained to y=0: "
            f"{figure_diagnostics.get('native_com_markers_are_xz_projections_not_y0_membership')}",
        ])
        if figure_diagnostics.get("superseded_previous_png"):
            lines.append(
                f"- Superseded clipped figure retained at: "
                f"{figure_diagnostics['superseded_previous_png']}"
            )
    lines.extend([
        "- Exact OCC common partition addresses CAD geometric volume coverage only. "
        "An 8 mm cell spacing exceeds the nominal 4 mm wall and is not adequate "
        "evidence of physical-resolution quality.",
        "- Cell size does not guarantee global minimum-wall resolution. One centroid "
        "per connected clipped solid is irregular quadrature and may bias directional "
        "fields; free-surface correction remains pending.",
        "- Physical resolution, curvature-aware refinement, crack paths, fracture "
        "energy, bond calibration, mechanics, and R2 are NOT_VALIDATED or PENDING.",
        "- Rho is an UNVALIDATED_PROTOTYPE placeholder; the mass result demonstrates "
        "only arithmetic conservation under that placeholder.",
        "- A native COM outside the canonical source is a candidate FAIL. The exact "
        "native centroid is retained for evidence; it is not moved and no alternate "
        "interior point is substituted for the same quadrature rule.",
        "",
        "## Artifacts",
        "",
        f"- JSON report: {report['artifacts']['report_json']}",
        f"- Per-box ledger: {report['artifacts']['cell_records_json']}",
        f"- Append-only run log: {report['artifacts']['run_log']}",
        f"- Candidate arrays: {report['artifacts']['candidate_npz']}",
        f"- Every retained clipped solid: {report['artifacts']['component_brep_dir']}",
        f"- Figure: {report['artifacts']['figure_png']}",
        "",
        "Reproduction from the repository root in PowerShell:",
        "",
        "    & '.\\third_party\\vesselforge-runtime\\venv\\Scripts\\python.exe' "
        "'pgfrag\\experiments\\cad_clipped_cells_spike.py'",
        "",
    ])
    if "execution" in report:
        execution = report["execution"]
        lines.extend([
            "## Process execution",
            "",
            f"- Argv: `{json.dumps(execution['argv'], ensure_ascii=False)}`",
            f"- Exit code: {execution['exit_code']}",
            f"- Runtime: {execution['runtime_seconds']:.3f} s",
            f"- Started: {execution['started_utc']}; finished: {execution['finished_utc']}",
            "",
        ])
    path.write_text("\n".join(lines), encoding="utf-8")


def run(output_dir: Path, logger: logging.Logger) -> dict[str, Any]:
    report: dict[str, Any] = {
        "schema": "pgfrag.cad_clipped_cells_candidate_spike/1.0",
        "created_utc": utc_now(),
        "status": "FAIL_STATIC_NATIVE_GEOMETRY_CANDIDATE",
        "source": {},
        "policy": {
            "spacing_mm": SPACING_MM,
            "max_candidate_boxes": MAX_BOXES,
            "boolean_wall_limit_seconds": BOOLEAN_WALL_LIMIT_SECONDS,
            "max_relative_volume_error_pct": VOLUME_ERROR_LIMIT * 100.0,
            "occ_boolean": "OCP.BOPAlgo_BOP / BOPAlgo_COMMON",
            "run_parallel": False,
            "native_operation_errors_and_warnings_reported": True,
            "empty_3d_material_policy": (
                "native result may be omitted only when null, or when it is a valid "
                "TopAbs_COMPOUND "
                "zero-volume shape with no solids, shells, faces, wires, edges, or vertices"
            ),
            "centroid_out_policy": "FAIL immediately; no random resampling",
            "pg_job_solver_bond_graph": "NOT_CREATED",
            "training_eligibility": "INELIGIBLE",
            "classifier": "OCP.BRepClass3d_SolidClassifier",
            "volume_integral": "OCP.BRepGProp.VolumeProperties_s with UseTriangulation=False",
            "quadrature": "one exact OCC center of mass + exact native volume per connected common solid",
        },
        "geometry": {},
        "counts": {},
        "checks": {},
        "material": {},
        "timing": {},
        "gates": {},
        "failures": [],
        "artifacts": {},
    }
    run_started = time.monotonic()
    boolean_seconds = 0.0
    box_records: list[dict[str, Any]] = []
    component_records: list[dict[str, Any]] = []
    component_shapes: list[cq.Solid] = []
    exceptions: list[dict[str, Any]] = []
    operation_failures = 0
    nonvolume_nonempty = 0
    exact_empty = 0
    exact_null = 0
    native_empty_compound = 0
    processed = 0
    abort_reason: str | None = None

    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    step_rel = manifest["geometry"]["solid_step"]
    step_path = (MANIFEST.parent / step_rel).resolve()
    expected_file = manifest["files"][step_rel]
    step_hash = sha256_file(step_path)
    hash_matches = (
        step_hash == expected_file["sha256"] == EXPECTED_STEP_SHA256
        and step_path.stat().st_size == expected_file["bytes"]
    )
    report["source"] = {
        "manifest_path": str(MANIFEST),
        "step_path": str(step_path),
        "step_bytes": step_path.stat().st_size,
        "step_sha256": step_hash,
        "manifest_sha256": expected_file["sha256"],
        "manifest_hash_verified": hash_matches,
        "frame": manifest["frame"],
        "length_unit": manifest["units"]["length"],
        "mass_unit": manifest["units"]["mass"],
        "backend": manifest["provenance"]["backend"],
        "backend_version": manifest["provenance"]["backend_version"],
    }
    if not hash_matches:
        raise ExperimentFailure(f"canonical STEP hash/size mismatch: {step_hash}")
    if manifest["frame"] != "right_handed_z_up" or manifest["units"]["length"] != "mm":
        raise ExperimentFailure("unexpected source frame or length unit")
    logger.info("canonical STEP verified: %s (%d bytes)", step_hash, step_path.stat().st_size)

    imported = cq.importers.importStep(str(step_path))
    imported_shape = imported.val()
    source_solids = imported_shape.Solids()
    if len(source_solids) != 1:
        raise ExperimentFailure(f"expected one canonical solid, found {len(source_solids)}")
    cad_solid = source_solids[0]
    source_valid = bool(cad_solid.isValid())
    if not source_valid:
        raise ExperimentFailure("canonical native solid is invalid")
    cad_volume, cad_center = exact_volume_and_centroid(cad_solid)
    bbox = exact_bbox(cad_solid)
    xmin, xmax, ymin, ymax, zmin, zmax = bbox
    boxes, axes = make_box_records(bbox, SPACING_MM)
    if len(boxes) > MAX_BOXES:
        raise ExperimentFailure(
            f"fixed grid requires {len(boxes)} boxes, greater than cap {MAX_BOXES}"
        )
    aabb_sum = math.fsum(row["aabb_volume_mm3"] for row in boxes)
    bbox_volume = (xmax - xmin) * (ymax - ymin) * (zmax - zmin)
    coverage_error = abs(aabb_sum - bbox_volume) / bbox_volume
    edges_ok = all(
        math.isclose(edges[0], lo, abs_tol=1e-10)
        and math.isclose(edges[-1], hi, abs_tol=1e-10)
        and all(b > a for a, b in zip(edges, edges[1:]))
        for edges, lo, hi in (
            (axes[0], xmin, xmax), (axes[1], ymin, ymax), (axes[2], zmin, zmax)
        )
    )
    if coverage_error > 1e-12 or not edges_ok:
        raise ExperimentFailure("Cartesian cells do not exactly cover the CAD bbox")
    report["geometry"] = {
        "source_shape_valid": source_valid,
        "cad_solid_count": len(source_solids),
        "cad_volume_mm3": cad_volume,
        "manifest_expected_volume_mm3": manifest["provenance"]["spec"].get(
            "expected_volume_mm3"
        ),
        "cad_centroid_mm": list(cad_center),
        "bbox_mm": list(bbox),
        "bbox_volume_mm3": bbox_volume,
        "aabb_partition_volume_mm3": aabb_sum,
        "aabb_partition_relative_error": coverage_error,
        "axis_edges_mm": [list(axis) for axis in axes],
        "axis_edges_exact_cover": edges_ok,
    }
    report["counts"]["candidate_boxes"] = len(boxes)
    logger.info(
        "CAD volume %.9f mm3; bbox %s; grid %d x %d x %d = %d cells",
        cad_volume, bbox, len(axes[0]) - 1, len(axes[1]) - 1,
        len(axes[2]) - 1, len(boxes),
    )

    cell_jsonl = output_dir / "cell_records.jsonl"
    classifier = BRepClass3d_SolidClassifier(cad_solid.wrapped)
    with cell_jsonl.open("w", encoding="utf-8", newline="\n") as ledger:
        for index, box_record in enumerate(boxes):
            if (boolean_seconds >= BOOLEAN_WALL_LIMIT_SECONDS
                    or time.monotonic() - run_started >= BOOLEAN_WALL_LIMIT_SECONDS):
                abort_reason = (
                    f"cumulative native Common or whole-run time reached "
                    f"{BOOLEAN_WALL_LIMIT_SECONDS:.0f} s"
                )
                logger.error("STOP: %s", abort_reason)
                for remaining in boxes[index:]:
                    row = {**remaining, "status": "NOT_RUN_BOOLEAN_TIME_LIMIT",
                           "exception": abort_reason}
                    box_records.append(row)
                    ledger.write(json.dumps(row, sort_keys=True) + "\n")
                ledger.flush()
                break
            x0, x1, y0, y1, z0, z1 = box_record["bounds_mm"]
            native_box = cq.Solid.makeBox(
                x1 - x0, y1 - y0, z1 - z0, pnt=cq.Vector(x0, y0, z0)
            )
            row: dict[str, Any] = {
                **box_record,
                "status": "PENDING",
                "boolean_threads": "single",
                "common_is_done": None,
                "common_is_null": None,
                "components": [],
            }
            started = time.monotonic()
            try:
                (common_shape, is_done, is_null, has_errors, has_warnings,
                 kernel_report_json) = native_common(cad_solid, native_box)
                elapsed = time.monotonic() - started
                boolean_seconds += elapsed
                row["boolean_seconds"] = elapsed
                row["common_is_done"] = is_done
                row["common_is_null"] = is_null
                row["common_has_errors"] = has_errors
                row["common_has_warnings"] = has_warnings
                row["common_kernel_report_json"] = kernel_report_json
                if has_warnings:
                    logger.warning(
                        "cell %d BOPAlgo_BOP COMMON reported warnings: %s",
                        box_record["cell_id"], kernel_report_json,
                    )
                if not is_done:
                    operation_failures += 1
                    row["status"] = "KERNEL_COMMON_NOT_DONE"
                    row["exception"] = "BOPAlgo_BOP.HasErrors() returned true"
                    abort_reason = row["exception"]
                elif is_null:
                    exact_empty += 1
                    exact_null += 1
                    row["status"] = "EXACT_EMPTY_COMMON"
                elif not common_shape.isValid():
                    operation_failures += 1
                    row["status"] = "NON_NULL_COMMON_INVALID"
                    row["exception"] = "non-null Common result failed native IsValid()"
                    path = Path("nonvolume_nonempty") / f"cell_{box_record['cell_id']:04d}.brep"
                    save_shape_brep(common_shape, output_dir / path)
                    row["nonvolume_common_brep"] = str(path)
                    abort_reason = row["exception"]
                else:
                    solids = common_shape.Solids()
                    common_topology = topology_counts(common_shape)
                    positive_here = 0
                    for solid_index, clipped_solid in enumerate(solids):
                        component_id = len(component_records)
                        volume, center = exact_volume_and_centroid(clipped_solid)
                        finite = math.isfinite(volume) and all(math.isfinite(v) for v in center)
                        positive = finite and volume > 0.0
                        valid = bool(clipped_solid.isValid())
                        relative_path = (
                            Path("components")
                            / f"component_{component_id:05d}_cell_{box_record['cell_id']:04d}.brep"
                        )
                        save_shape_brep(clipped_solid, output_dir / relative_path)
                        center_state = classify(classifier, center) if finite else "UNKNOWN_NONFINITE"
                        component = {
                            "component_id": component_id,
                            "cell_id": box_record["cell_id"],
                            "connected_solid_index_within_common": solid_index,
                            "shape_brep": str(relative_path),
                            "volume_mm3": volume,
                            "center_of_mass_mm": list(center),
                            "center_state_in_canonical_cad": center_state,
                            "finite_volume_and_center": finite,
                            "positive_volume": positive,
                            "native_shape_valid": valid,
                            "bounds_mm": list(exact_bbox(clipped_solid)),
                        }
                        component_records.append(component)
                        component_shapes.append(clipped_solid)
                        row["components"].append(component)
                        if positive:
                            positive_here += 1
                        else:
                            operation_failures += 1
                            row["exception"] = (
                                f"connected common solid {solid_index} is not finite positive volume"
                            )
                            abort_reason = row["exception"]
                        if not valid:
                            operation_failures += 1
                            row["exception"] = f"connected common solid {solid_index} is invalid"
                            abort_reason = row["exception"]
                        if center_state != "IN":
                            operation_failures += 1
                            row.setdefault("center_classification_failures", []).append({
                                "component_id": component_id,
                                "state": center_state,
                                "center_mm": list(center),
                            })
                            abort_reason = (
                                f"native center of mass for component {component_id} "
                                f"classified {center_state}; no resampling attempted"
                            )
                    row["positive_components_in_cell"] = positive_here
                    if not solids:
                        common_props = GProp_GProps()
                        BRepGProp.VolumeProperties_s(
                            common_shape.wrapped, common_props, False, False, False
                        )
                        common_material_volume = float(common_props.Mass())
                        row["common_topology_counts"] = common_topology
                        row["common_material_volume_mm3"] = common_material_volume
                        empty_topology = all(value == 0 for value in common_topology.values())
                        if (
                            common_shape.isValid()
                            and common_shape.wrapped.ShapeType() == TopAbs_COMPOUND
                            and empty_topology
                            and math.isfinite(common_material_volume)
                            and common_material_volume == 0.0
                        ):
                            native_empty_compound += 1
                            exact_empty += 1
                            row["status"] = "EMPTY_3D_MATERIAL"
                        else:
                            nonvolume_nonempty += 1
                            operation_failures += 1
                            row["status"] = "NON_NULL_COMMON_WITHOUT_SOLID"
                            row["exception"] = (
                                "non-null Common has nonempty topology or nonzero/nonfinite "
                                "material volume but no solid"
                            )
                            path = Path("nonvolume_nonempty") / f"cell_{box_record['cell_id']:04d}.brep"
                            save_shape_brep(common_shape, output_dir / path)
                            row["nonvolume_common_brep"] = str(path)
                            abort_reason = row["exception"]
                    elif row.get("exception") or abort_reason:
                        row["status"] = "COMPONENT_GATE_FAILURE"
                    else:
                        row["status"] = "POSITIVE_COMPONENTS_RETAINED"
                if abort_reason:
                    row["run_abort_after_this_cell"] = abort_reason
                box_records.append(row)
                processed += 1
                ledger.write(json.dumps(row, sort_keys=True, allow_nan=False) + "\n")
                ledger.flush()
                if (index + 1) % 25 == 0 or index + 1 == len(boxes) or abort_reason:
                    logger.info(
                        "cell %d/%d %s; solids=%d; boolean %.3fs, cumulative %.3fs",
                        index + 1, len(boxes), row["status"], len(row["components"]),
                        elapsed, boolean_seconds,
                    )
                if abort_reason:
                    for remaining in boxes[index + 1:]:
                        missing = {
                            **remaining,
                            "status": "NOT_RUN_AFTER_KERNEL_OR_GEOMETRY_FAILURE",
                            "exception": abort_reason,
                        }
                        box_records.append(missing)
                        ledger.write(json.dumps(missing, sort_keys=True) + "\n")
                    ledger.flush()
                    break
                if (boolean_seconds > BOOLEAN_WALL_LIMIT_SECONDS
                        or time.monotonic() - run_started > BOOLEAN_WALL_LIMIT_SECONDS):
                    abort_reason = (
                        f"cumulative native Common or whole-run time exceeded "
                        f"{BOOLEAN_WALL_LIMIT_SECONDS:.0f} s"
                    )
                    logger.error("STOP: %s", abort_reason)
                    for remaining in boxes[index + 1:]:
                        missing = {**remaining, "status": "NOT_RUN_BOOLEAN_TIME_LIMIT",
                                   "exception": abort_reason}
                        box_records.append(missing)
                        ledger.write(json.dumps(missing, sort_keys=True) + "\n")
                    ledger.flush()
                    break
            except Exception as exc:
                elapsed = time.monotonic() - started
                boolean_seconds += elapsed
                exceptions.append({
                    "cell_id": box_record["cell_id"],
                    "exception_type": type(exc).__name__,
                    "message": str(exc),
                    "elapsed_seconds": elapsed,
                })
                row["status"] = "KERNEL_OR_COMPONENT_EXCEPTION"
                row["boolean_seconds"] = elapsed
                row["exception"] = f"{type(exc).__name__}: {exc}"
                box_records.append(row)
                ledger.write(json.dumps(row, sort_keys=True) + "\n")
                ledger.flush()
                abort_reason = row["exception"]
                logger.exception("OCC operation failed at cell %d", box_record["cell_id"])
                for remaining in boxes[index + 1:]:
                    missing = {
                        **remaining,
                        "status": "NOT_RUN_AFTER_KERNEL_OR_GEOMETRY_FAILURE",
                        "exception": abort_reason,
                    }
                    box_records.append(missing)
                    ledger.write(json.dumps(missing, sort_keys=True) + "\n")
                ledger.flush()
                break

    all_processed = processed == len(boxes) and len(box_records) == len(boxes)
    run_state = (
        "COMPLETE" if all_processed else
        "INCOMPLETE_COM_OUT" if abort_reason and "center of mass" in abort_reason else
        "INCOMPLETE_AFTER_KERNEL_OR_GEOMETRY_FAILURE"
    )
    points = np.asarray(
        [c["center_of_mass_mm"] for c in component_records], dtype=np.float64
    ).reshape((-1, 3))
    volumes_mm3 = np.asarray(
        [c["volume_mm3"] for c in component_records], dtype=np.float64
    )
    cell_ids = np.asarray([c["cell_id"] for c in component_records], dtype=np.int64)
    component_ids = np.asarray([c["component_id"] for c in component_records], dtype=np.int64)
    finite_positive = bool(
        len(volumes_mm3) > 0 and np.isfinite(points).all()
        and np.isfinite(volumes_mm3).all() and (volumes_mm3 > 0.0).all()
    )
    volume_sum = float(math.fsum(float(v) for v in volumes_mm3))
    partial_volume_fraction = volume_sum / cad_volume if cad_volume else math.nan
    volume_error = (
        abs(volume_sum - cad_volume) / cad_volume
        if all_processed and cad_volume else None
    )
    unique_centroids = {
        tuple(round(float(v), 9) for v in point) for point in points
    }
    centroid_states: dict[str, int] = {}
    for component in component_records:
        state = component["center_state_in_canonical_cad"]
        centroid_states[state] = centroid_states.get(state, 0) + 1

    cavity_results = [
        {"point_mm": list(point), "state": classify(classifier, point)}
        for point in cavity_probes()
    ]
    cavity_states: dict[str, int] = {}
    for probe in cavity_results:
        cavity_states[probe["state"]] = cavity_states.get(probe["state"], 0) + 1
    cavity_empty = bool(cavity_results) and all(
        probe["state"] == "OUT" for probe in cavity_results
    )

    material = manifest["material"]
    rho = float(material["mechanics"]["rho_kg_m3"])
    mass_partial = volume_sum * 1.0e-9 * rho
    mass_candidate = mass_partial if all_processed else None
    mass_cad = cad_volume * 1.0e-9 * rho
    mass_error = (
        abs(mass_partial - mass_cad) / mass_cad
        if all_processed and mass_cad else None
    )
    provenance = {
        "schema": report["schema"],
        "representation": "OCC-native CAD-clipped Cartesian volume cells",
        "static_geometry_candidate_only": all_processed,
        "partial_diagnostic_only": not all_processed,
        "partition_complete": all_processed,
        "run_state": run_state,
        "processed_boxes": processed,
        "candidate_boxes": len(boxes),
        "not_accepted_pg_job": True,
        "not_solver_input": True,
        "not_bond_graph": True,
        "training_eligible": False,
        "source_step_sha256": step_hash,
        "frame": manifest["frame"],
        "length_unit": "mm",
        "volume_unit": "mm3",
        "mass_unit": "kg",
        "grid_spacing_mm": SPACING_MM,
        "method": report["policy"],
    }
    component_bounds = np.asarray(
        [c["bounds_mm"] for c in component_records], dtype=np.float64
    ).reshape((-1, 6))
    parent_bounds = np.asarray(
        [box_records[c["cell_id"]]["bounds_mm"] for c in component_records],
        dtype=np.float64,
    ).reshape((-1, 6))
    output_npz = output_dir / "cad_clipped_cells_candidate.npz"
    points_m = points * 1.0e-3
    volumes_m3 = volumes_mm3 * 1.0e-9
    masses_kg = volumes_m3 * rho
    np.savez_compressed(
        output_npz,
        points_mm=points,
        points_m=points_m,
        cell_volume_mm3=volumes_mm3,
        cell_volume_m3=volumes_m3,
        sourcecellIDs=cell_ids,
        componentIDs=component_ids,
        rho_kg_m3=np.asarray(rho, dtype=np.float64),
        cell_mass_kg=masses_kg,
        native_center_state=np.asarray(
            [c["center_state_in_canonical_cad"] for c in component_records], dtype="U24"
        ),
        component_bounds_mm=component_bounds,
        sourcebox_bounds_mm=parent_bounds,
        source_step_sha256=np.asarray(step_hash, dtype="U64"),
        frame=np.asarray(manifest["frame"], dtype="U32"),
        length_unit=np.asarray("mm", dtype="U4"),
        volume_unit=np.asarray("mm3", dtype="U8"),
        mass_unit=np.asarray("kg", dtype="U4"),
        neutral_candidate_provenance_json=np.asarray(
            json.dumps(provenance, sort_keys=True), dtype="U4096"
        ),
    )
    write_json(output_dir / "cell_records.json", box_records)
    write_json(output_dir / "component_records.json", component_records)
    write_json(output_dir / "cavity_probe_results.json", cavity_results)

    figure_path = output_dir / "cad_clipped_cells_y0.png"
    figure_error: str | None = None
    figure_diagnostics: dict[str, Any] = {}
    try:
        figure_diagnostics = make_figure(
            figure_path, cad_solid, box_records, component_records,
            component_shapes, points, bbox, logger,
        )
    except Exception as exc:
        figure_error = f"{type(exc).__name__}: {exc}"
        exceptions.append({
            "stage": "figure",
            "exception_type": type(exc).__name__,
            "message": str(exc),
        })
        logger.exception("native-section figure failed")

    all_centroids_inside = bool(component_records) and all(
        c["center_state_in_canonical_cad"] == "IN" for c in component_records
    )
    all_shapes_valid = bool(component_records) and all(
        c["native_shape_valid"] for c in component_records
    )
    no_nonempty_nonsolid = nonvolume_nonempty == 0 and all(
        row.get("status") not in {
            "NON_NULL_COMMON_INVALID", "NON_NULL_COMMON_WITHOUT_SOLID"
        } for row in box_records
    )
    boolean_within_time = boolean_seconds <= BOOLEAN_WALL_LIMIT_SECONDS
    whole_run_seconds = time.monotonic() - run_started
    run_within_time = whole_run_seconds <= BOOLEAN_WALL_LIMIT_SECONDS
    no_exceptions = not exceptions and operation_failures == 0
    volume_gate = (
        all_processed and volume_error is not None
        and volume_error <= VOLUME_ERROR_LIMIT
    )
    units_and_weights_ok = bool(
        len(volumes_mm3) > 0
        and manifest["units"]["length"] == "mm"
        and math.isfinite(rho) and rho > 0.0
        and np.isfinite(points).all() and np.isfinite(points_m).all()
        and np.isfinite(volumes_mm3).all() and (volumes_mm3 > 0.0).all()
        and np.isfinite(volumes_m3).all() and (volumes_m3 > 0.0).all()
        and np.isfinite(masses_kg).all() and (masses_kg > 0.0).all()
        and np.allclose(points_m, points * 1.0e-3, rtol=0.0, atol=0.0)
        and np.allclose(volumes_m3, volumes_mm3 * 1.0e-9, rtol=0.0, atol=0.0)
    )
    native_kernel_clean = all(
        not row.get("common_has_errors", False)
        and not row.get("common_has_warnings", False)
        for row in box_records
        if row.get("common_has_errors") is not None
    )
    gates = {
        "canonical_manifest_and_step_hash_match": hash_matches,
        "canonical_one_valid_solid": source_valid and len(source_solids) == 1,
        "fixed_cartesian_grid_covers_exact_bbox": edges_ok and coverage_error <= 1e-12,
        "candidate_box_count_at_most_512": len(boxes) <= MAX_BOXES,
        "every_box_common_processed": all_processed,
        "native_boolean_single_thread": True,
        "boolean_wall_time_at_most_300_seconds": boolean_within_time,
        "whole_run_wall_time_at_most_300_seconds": run_within_time,
        "native_boolean_reports_no_errors_or_warnings": native_kernel_clean,
        "only_exact_or_native_valid_empty_3d_material_omitted": no_nonempty_nonsolid,
        "all_retained_components_valid_finite_positive": finite_positive and all_shapes_valid,
        "all_component_coms_natively_inside": all_centroids_inside,
        "cavity_probes_classify_out": cavity_empty,
        "native_volume_error_at_most_1_percent": volume_gate,
        "all_weights_and_units_finite_positive": units_and_weights_ok,
        "manifest_rho_mass_conservation_arithmetic": (
            all_processed and mass_error is not None
            and math.isfinite(mass_error) and mass_error <= VOLUME_ERROR_LIMIT
        ),
        "actual_occ_sections_png_created": figure_error is None and figure_path.exists(),
        "corrected_section_covers_left_right_walls_and_full_z_bbox": bool(
            figure_diagnostics
            and figure_diagnostics.get("source_section_spans_z0_to_z56")
            and figure_diagnostics.get("source_section_has_left_and_right_x")
            and all(
                band["negative_x_wall_sampled"] and band["positive_x_wall_sampled"]
                for band in figure_diagnostics.get("body_wall_samples_by_z_mm", {}).values()
            )
        ),
        "physics_and_pg_acceptance": False,
    }
    report["counts"] = {
        **report["counts"],
        "processed_boxes": processed,
        "box_records": len(box_records),
        "exact_null_common": exact_null,
        "empty_3d_material_results": exact_empty,
        "native_empty_compound": native_empty_compound,
        "positive_components": len(component_records),
        "nonvolume_nonempty_commons": nonvolume_nonempty,
        "kernel_error_results": sum(bool(row.get("common_has_errors")) for row in box_records),
        "kernel_warning_results": sum(bool(row.get("common_has_warnings")) for row in box_records),
        "exceptions": len(exceptions) + operation_failures,
        "operation_gate_failures": operation_failures,
        "component_brep_files": len(component_records),
    }
    report["geometry"].update({
        "candidate_volume_mm3": volume_sum if all_processed else None,
        "processed_component_volume_mm3": volume_sum,
        "processed_component_volume_fraction_of_cad": partial_volume_fraction,
        "volume_accuracy_status": (
            "EVALUATED_COMPLETE_PARTITION" if all_processed
            else "NOT_EVALUATED_INCOMPLETE_PARTITION"
        ),
        "run_state": run_state,
        "absolute_volume_error_mm3": (
            abs(volume_sum - cad_volume) if all_processed else None
        ),
        "relative_volume_error": volume_error,
        "relative_volume_error_pct": volume_error * 100.0 if volume_error is not None else None,
    })
    report["checks"] = {
        "finite_positive_components": sum(
            bool(c["finite_volume_and_center"] and c["positive_volume"])
            for c in component_records
        ),
        "centroid_state_counts": centroid_states,
        "all_centroids_inside": all_centroids_inside,
        "unique_centroid_count": len(unique_centroids),
        "cavity_probe_count": len(cavity_results),
        "cavity_probe_state_counts": cavity_states,
        "cavity_probe_scope": (
            "axis-centered points at radii 0, 4, 8 mm (four azimuths for nonzero "
            "radii), z=8,16,24,32,40,48,52 mm"
        ),
        "all_cavity_probes_out": cavity_empty,
        "all_component_shapes_valid": all_shapes_valid,
        "source_cad_section_bounds_xmin_xmax_zmin_zmax_mm": figure_diagnostics.get(
            "source_cad_section_bounds_xmin_xmax_zmin_zmax_mm"
        ),
        "native_com_markers_are_xz_projections_not_y0_membership": True,
    }
    report["figure_diagnostics"] = figure_diagnostics
    report["material"] = {
        "rho_kg_m3": rho,
        "calibration_status": material["calibration_status"],
        "candidate_mass_kg": mass_candidate,
        "processed_partial_mass_kg": mass_partial,
        "cad_mass_kg": mass_cad,
        "relative_mass_error": mass_error,
        "relative_mass_error_pct": mass_error * 100.0 if mass_error is not None else None,
        "mass_accuracy_status": (
            "EVALUATED_COMPLETE_PARTITION" if all_processed
            else "NOT_EVALUATED_INCOMPLETE_PARTITION"
        ),
        "density_is_unvalidated_placeholder": (
            material["calibration_status"] == "UNVALIDATED_PROTOTYPE"
        ),
    }
    report["timing"] = {
        "boolean_seconds": boolean_seconds,
        "whole_run_seconds": whole_run_seconds,
        "boolean_wall_limit_seconds": BOOLEAN_WALL_LIMIT_SECONDS,
    }
    report["gates"] = gates
    gate_reasons = {
        "canonical_manifest_and_step_hash_match": (
            f"expected={EXPECTED_STEP_SHA256}; actual={step_hash}"
        ),
        "canonical_one_valid_solid": (
            f"native solids={len(source_solids)}; valid={source_valid}"
        ),
        "fixed_cartesian_grid_covers_exact_bbox": (
            f"relative bbox partition error={coverage_error:.12g}; exact_edges={edges_ok}"
        ),
        "candidate_box_count_at_most_512": (
            f"candidate boxes={len(boxes)}; cap={MAX_BOXES}"
        ),
        "every_box_common_processed": f"processed={processed}; required={len(boxes)}",
        "boolean_wall_time_at_most_300_seconds": (
            f"boolean_seconds={boolean_seconds:.6f}; limit={BOOLEAN_WALL_LIMIT_SECONDS:.0f}"
        ),
        "whole_run_wall_time_at_most_300_seconds": (
            f"whole_run_seconds={whole_run_seconds:.6f}; limit={BOOLEAN_WALL_LIMIT_SECONDS:.0f}"
        ),
        "native_boolean_reports_no_errors_or_warnings": (
            f"error_results={report['counts'].get('kernel_error_results', 0)}; "
            f"warning_results={report['counts'].get('kernel_warning_results', 0)}"
        ),
        "only_exact_or_native_valid_empty_3d_material_omitted": (
            f"nonvolume_nonempty={nonvolume_nonempty}; native_empty_compound={native_empty_compound}"
        ),
        "all_retained_components_valid_finite_positive": (
            f"finite_positive={finite_positive}; all_valid={all_shapes_valid}; "
            f"component_count={len(component_records)}"
        ),
        "all_component_coms_natively_inside": f"native COM states={centroid_states}",
        "cavity_probes_classify_out": f"native cavity probe states={cavity_states}",
        "native_volume_error_at_most_1_percent": (
            "NOT_EVALUATED_INCOMPLETE_PARTITION"
            if volume_error is None else
            f"relative_error_pct={volume_error * 100.0:.6f}; limit=1.000000"
        ),
        "all_weights_and_units_finite_positive": (
            "validated exact mm-to-m conversion, finite positive component volumes, "
            "finite positive masses, and manifest millimeter units"
        ),
        "manifest_rho_mass_conservation_arithmetic": (
            "NOT_EVALUATED_INCOMPLETE_PARTITION"
            if mass_error is None else
            f"rho={rho}; relative_mass_error_pct={mass_error * 100.0:.6f}"
        ),
        "actual_occ_sections_png_created": (
            figure_error or f"exists={figure_path.exists()}"
        ),
        "physics_and_pg_acceptance": (
            "false by design: static geometry candidate only; no solver/job/physical qualification"
        ),
    }
    report["gate_failure_reasons"] = {
        name: gate_reasons.get(name, abort_reason or "gate did not pass")
        for name, passed in gates.items() if not passed
    }
    report["failures"] = [
        {
            "gate": name,
            "reason": report["gate_failure_reasons"][name],
        }
        for name, passed in gates.items() if not passed
    ]
    report["exceptions_detail"] = exceptions
    report["abort_reason"] = abort_reason
    report["run_state"] = run_state
    report["artifacts"] = {
        "report_json": str(output_dir / "report.json"),
        "report_md": str(output_dir / "report.md"),
        "cell_records_json": str(output_dir / "cell_records.json"),
        "cell_records_jsonl": str(cell_jsonl),
        "component_records_json": str(output_dir / "component_records.json"),
        "cavity_probes_json": str(output_dir / "cavity_probe_results.json"),
        "run_log": str(output_dir / "run.log"),
        "candidate_npz": str(output_npz),
        "component_brep_dir": str(output_dir / "components"),
        "nonvolume_common_dir": str(output_dir / "nonvolume_nonempty"),
        "figure_png": str(figure_path),
        "figure_png_sha256": sha256_file(figure_path) if figure_path.exists() else None,
        "figure_diagnostics_json": str(output_dir / "figure_diagnostics.json"),
    }
    if figure_error:
        report["figure_error"] = figure_error
    nominal_wall = float(manifest["provenance"]["spec"]["parameters"]["wall_mm"])
    report["limitations"] = {
        "nominal_source_wall_mm": nominal_wall,
        "grid_spacing_exceeds_nominal_wall": SPACING_MM > nominal_wall,
        "global_minimum_wall_guarantee": False,
        "curvature_aware_refinement": "PENDING",
        "free_surface_correction": "PENDING",
        "directional_quadrature_bias": "PENDING",
        "physical_resolution": "NOT_VALIDATED",
        "fracture_crack_paths_and_bonds": "NOT_IMPLEMENTED",
        "solver_or_pg_job_support": "NOT_IMPLEMENTED",
        "training_eligibility": "INELIGIBLE",
        "r2": "NOT_IMPLEMENTED",
    }
    geometry_gates_pass = all(
        passed for name, passed in gates.items() if name != "physics_and_pg_acceptance"
    )
    report["status"] = (
        "PASS_STATIC_NATIVE_GEOMETRY_CANDIDATE_ONLY"
        if geometry_gates_pass else "FAIL_STATIC_NATIVE_GEOMETRY_CANDIDATE"
    )
    write_json(output_dir / "report.json", report)
    write_report_md(output_dir / "report.md", report)
    volume_summary = (
        "NOT_EVALUATED_INCOMPLETE_PARTITION"
        if volume_error is None else f"{volume_error * 100.0:.6f}%"
    )
    logger.info(
        "final status=%s run-state=%s boxes=%d/%d components=%d "
        "full-partition-volume-error=%s COM=%s cavity=%s boolean=%.3fs",
        report["status"], run_state, processed, len(boxes), len(component_records),
        volume_summary, centroid_states, cavity_states, boolean_seconds,
    )
    return report


def normalize_saved_incomplete_report() -> int:
    """Mark the saved COM_OUT run incomplete without rerunning any geometry."""
    output_dir = OUTPUT_DIR
    report_path = output_dir / "report.json"
    if not report_path.exists():
        print(f"No saved candidate report to normalize: {report_path}", file=sys.stderr)
        return 2
    report = json.loads(report_path.read_text(encoding="utf-8"))
    abort_reason = str(report.get("abort_reason", ""))
    processed = int(report.get("counts", {}).get("processed_boxes", 0))
    candidate_boxes = int(report.get("counts", {}).get("candidate_boxes", 0))
    if "center of mass" not in abort_reason or processed >= candidate_boxes:
        print("Saved run is not the expected incomplete COM_OUT case", file=sys.stderr)
        return 2

    geometry = report["geometry"]
    material = report["material"]
    partial_volume = float(geometry["candidate_volume_mm3"])
    cad_volume = float(geometry["cad_volume_mm3"])
    partial_fraction = partial_volume / cad_volume
    partial_mass = float(material["candidate_mass_kg"])
    report["status"] = "FAIL_STATIC_NATIVE_GEOMETRY_CANDIDATE"
    report["run_state"] = "INCOMPLETE_COM_OUT"
    report["geometry"].update({
        "candidate_volume_mm3": None,
        "processed_component_volume_mm3": partial_volume,
        "processed_component_volume_fraction_of_cad": partial_fraction,
        "volume_accuracy_status": "NOT_EVALUATED_INCOMPLETE_PARTITION",
        "run_state": "INCOMPLETE_COM_OUT",
        "absolute_volume_error_mm3": None,
        "relative_volume_error": None,
        "relative_volume_error_pct": None,
    })
    report["material"].update({
        "candidate_mass_kg": None,
        "processed_partial_mass_kg": partial_mass,
        "relative_mass_error": None,
        "relative_mass_error_pct": None,
        "mass_accuracy_status": "NOT_EVALUATED_INCOMPLETE_PARTITION",
    })
    report.setdefault("volume_evaluation", {}).update({
        "status": "NOT_EVALUATED_INCOMPLETE_PARTITION",
        "partition_complete": False,
        "processed_boxes": processed,
        "candidate_boxes": candidate_boxes,
        "processed_component_volume_subtotal_mm3": partial_volume,
        "processed_component_volume_fraction_of_full_cad": partial_fraction,
        "full_partition_relative_volume_error_pct": None,
        "partial_subtotal_difference_is_not_a_partition_error": True,
    })
    report.setdefault("gate_failure_reasons", {}).update({
        "every_box_common_processed": f"processed={processed}/{candidate_boxes}; stopped on native COM OUT",
        "native_volume_error_at_most_1_percent": "NOT_EVALUATED_INCOMPLETE_PARTITION",
        "manifest_rho_mass_conservation_arithmetic": "NOT_EVALUATED_INCOMPLETE_PARTITION",
        "all_component_coms_natively_inside": abort_reason,
        "actual_occ_sections_png_created": (
            "previous figure failed axis-aligned section coverage review and is superseded"
        ),
        "corrected_section_covers_left_right_walls_and_full_z_bbox": (
            "awaiting saved-geometry-only diagnostic render"
        ),
    })
    report.setdefault("gates", {}).update({
        "every_box_common_processed": False,
        "native_volume_error_at_most_1_percent": False,
        "manifest_rho_mass_conservation_arithmetic": False,
        "all_component_coms_natively_inside": False,
        "actual_occ_sections_png_created": False,
        "corrected_section_covers_left_right_walls_and_full_z_bbox": False,
    })
    report["failures"] = [
        {
            "gate": name,
            "reason": report["gate_failure_reasons"].get(name, "gate did not pass"),
        }
        for name, passed in report["gates"].items() if not passed
    ]
    report.setdefault("checks", {}).update({
        "all_boxes_processed": False,
        "all_centroids_inside": False,
        "full_partition_volume_accuracy": "NOT_EVALUATED_INCOMPLETE_PARTITION",
    })
    report.setdefault("limitations", {}).update({
        "full_partition_volume_accuracy": "NOT_EVALUATED_INCOMPLETE_PARTITION",
        "training_eligibility": "INELIGIBLE",
        "solver_or_pg_job_support": "NOT_IMPLEMENTED",
    })
    started_dt = datetime.fromisoformat(report["created_utc"])
    experiment_runtime = float(report.get("timing", {}).get("whole_run_seconds", 0.0))
    report["execution"] = {
        "argv": [sys.executable, str(Path(__file__).resolve())],
        "raw_sys_argv": [str(Path(__file__).resolve())],
        "started_utc": report["created_utc"],
        "finished_utc": (started_dt + timedelta(seconds=experiment_runtime)).isoformat(
            timespec="milliseconds"
        ),
        "runtime_seconds": experiment_runtime,
        "runtime_scope": "experiment from canonical STEP import through corrected-source figure generation; excludes interpreter startup",
        "exit_code": 1,
        "python_executable": sys.executable,
        "run_state": "INCOMPLETE_COM_OUT",
    }

    npz_path = Path(report["artifacts"]["candidate_npz"])
    with np.load(npz_path, allow_pickle=False) as archive:
        arrays = {key: archive[key] for key in archive.files}
    provenance = json.loads(str(arrays["neutral_candidate_provenance_json"].item()))
    provenance.update({
        "static_geometry_candidate_only": False,
        "partial_diagnostic_only": True,
        "partition_complete": False,
        "run_state": "INCOMPLETE_COM_OUT",
        "processed_boxes": processed,
        "candidate_boxes": candidate_boxes,
        "full_partition_volume_accuracy": "NOT_EVALUATED_INCOMPLETE_PARTITION",
        "not_accepted_pg_job": True,
        "not_solver_input": True,
        "not_bond_graph": True,
        "training_eligible": False,
    })
    arrays["neutral_candidate_provenance_json"] = np.asarray(
        json.dumps(provenance, sort_keys=True), dtype="U4096"
    )
    arrays["partition_complete"] = np.asarray(False, dtype=np.bool_)
    arrays["run_state"] = np.asarray("INCOMPLETE_COM_OUT", dtype="U32")
    arrays["processed_box_count"] = np.asarray(processed, dtype=np.int64)
    arrays["candidate_box_count"] = np.asarray(candidate_boxes, dtype=np.int64)
    arrays["processed_component_volume_fraction_of_cad"] = np.asarray(
        partial_fraction, dtype=np.float64
    )
    np.savez_compressed(npz_path, **arrays)
    report["artifacts"]["candidate_npz_sha256"] = sha256_file(npz_path)
    report["normalization_execution"] = {
        "argv": [sys.executable, str(Path(__file__).resolve()), "--normalize-saved-incomplete"],
        "finished_utc": utc_now(),
        "exit_code": 0,
        "performed_common_boolean": False,
    }
    write_json(report_path, report)
    write_report_md(output_dir / "report.md", report)
    with Path(report["artifacts"]["run_log"]).open(
        "a", encoding="utf-8", newline="\n"
    ) as stream:
        stream.write(
            f"{utc_now()} AUDIT: marked {processed}/{candidate_boxes} COM_OUT run as "
            "INCOMPLETE_COM_OUT; partial retained volume is subtotal coverage only; "
            "full-partition volume and mass accuracy NOT_EVALUATED; no Common rerun\n"
        )
    return 0


def render_saved_diagnostic() -> int:
    """Re-render saved OCC geometry only; this path performs no Boolean Common."""
    output_dir = OUTPUT_DIR
    report_path = output_dir / "report.json"
    if not report_path.exists():
        print(f"No saved candidate report to render: {report_path}", file=sys.stderr)
        return 2
    report = json.loads(report_path.read_text(encoding="utf-8"))
    started = time.monotonic()
    render_started_utc = utc_now()
    old_figure = Path(report["artifacts"]["figure_png"])
    old_figure_sha256 = sha256_file(old_figure) if old_figure.exists() else None
    new_figure = output_dir / "cad_clipped_cells_y0_ax3_corrected.png"
    if new_figure.exists():
        print(f"Refusing to overwrite diagnostic PNG: {new_figure}", file=sys.stderr)
        return 2

    step_path = Path(report["source"]["step_path"])
    actual_sha = sha256_file(step_path)
    if actual_sha != EXPECTED_STEP_SHA256:
        raise ExperimentFailure(f"saved diagnostic source STEP hash mismatch: {actual_sha}")
    imported = cq.importers.importStep(str(step_path)).val()
    solids = imported.Solids()
    if len(solids) != 1 or not solids[0].isValid():
        raise ExperimentFailure("saved diagnostic source is not one valid native solid")
    cad_solid = solids[0]

    component_records = json.loads(
        Path(report["artifacts"]["component_records_json"]).read_text(encoding="utf-8")
    )
    box_records = json.loads(
        Path(report["artifacts"]["cell_records_json"]).read_text(encoding="utf-8")
    )
    component_shapes = [
        cq.importers.importBrep(str(output_dir / component["shape_brep"])).val()
        for component in component_records
    ]
    points = np.asarray(
        [component["center_of_mass_mm"] for component in component_records],
        dtype=np.float64,
    ).reshape((-1, 3))
    bbox = tuple(float(v) for v in report["geometry"]["bbox_mm"])

    logger = logging.getLogger("cad_clipped_cells_diagnostic_render")
    logger.handlers.clear()
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)sZ %(levelname)s %(message)s")
    formatter.converter = time.gmtime
    for handler in (
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(output_dir / "run.log", mode="a", encoding="utf-8"),
    ):
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    diagnostics = make_figure(
        new_figure, cad_solid, box_records, component_records,
        component_shapes, points, bbox, logger,
    )
    wall_bands = diagnostics["body_wall_samples_by_z_mm"]
    wall_bands_ok = bool(wall_bands) and all(
        entry["negative_x_wall_sampled"] and entry["positive_x_wall_sampled"]
        for entry in wall_bands.values()
    )
    expected_xmin, expected_xmax, expected_zmin, expected_zmax = (
        diagnostics["expected_source_bbox_xmin_xmax_zmin_zmax_mm"]
    )
    actual_xmin, actual_xmax, actual_zmin, actual_zmax = (
        diagnostics["source_cad_section_bounds_xmin_xmax_zmin_zmax_mm"]
    )
    source_bbox_covered = bool(
        actual_xmin is not None
        and actual_xmin <= expected_xmin + 0.1
        and actual_xmax >= expected_xmax - 0.1
        and abs(actual_zmin - expected_zmin) <= 0.1
        and abs(actual_zmax - expected_zmax) <= 0.1
    )
    section_gate = bool(
        diagnostics["source_section_spans_z0_to_z56"]
        and diagnostics["source_section_has_left_and_right_x"]
        and wall_bands_ok
        and source_bbox_covered
    )
    superseded_figure = old_figure.with_name(
        "cad_clipped_cells_y0_SUPERSEDED_wrong_plane_uv.png"
    )
    if old_figure.exists() and not superseded_figure.exists():
        old_figure.rename(superseded_figure)
    render_finished = utc_now()
    report["figure_diagnostics"] = {
        **diagnostics,
        "source_bbox_section_extent_check_passed": source_bbox_covered,
        "left_and_right_body_wall_bands_passed": wall_bands_ok,
        "corrected_section_coverage_gate": section_gate,
        "superseded_previous_png_reason": (
            "previous bounded gp_Pln face used world X/Z values as local U/V without "
            "aligning plane axes, clipping left-wall and upper source sections"
        ),
        "superseded_previous_png": str(superseded_figure),
    }
    report["artifacts"]["superseded_figure_png"] = str(superseded_figure)
    report["artifacts"]["figure_png"] = str(new_figure)
    report["artifacts"]["figure_png_sha256"] = sha256_file(new_figure)
    report["artifacts"]["superseded_figure_png_sha256"] = old_figure_sha256
    report["artifacts"]["figure_diagnostics_json"] = str(
        output_dir / "figure_diagnostics.json"
    )
    report["gates"]["actual_occ_sections_png_created"] = new_figure.exists()
    report["gates"]["corrected_section_covers_left_right_walls_and_full_z_bbox"] = section_gate
    report.get("gate_failure_reasons", {}).pop("actual_occ_sections_png_created", None)
    report.setdefault("checks", {}).update({
        "source_cad_section_bounds_xmin_xmax_zmin_zmax_mm": (
            diagnostics["source_cad_section_bounds_xmin_xmax_zmin_zmax_mm"]
        ),
        "source_section_has_left_and_right_x": diagnostics["source_section_has_left_and_right_x"],
        "source_section_spans_z0_to_z56": diagnostics["source_section_spans_z0_to_z56"],
        "native_com_markers_are_xz_projections_not_y0_membership": True,
    })
    if not section_gate:
        report.setdefault("gate_failure_reasons", {})[
            "corrected_section_covers_left_right_walls_and_full_z_bbox"
        ] = "corrected exact OCC section did not meet source bbox/wall coverage checks"
    else:
        report.get("gate_failure_reasons", {}).pop(
            "corrected_section_covers_left_right_walls_and_full_z_bbox", None
        )
    report["failures"] = [
        {"gate": name, "reason": report.get("gate_failure_reasons", {}).get(name, "gate did not pass")}
        for name, passed in report["gates"].items() if not passed
    ]
    report["figure_render_execution"] = {
        "argv": [sys.executable, str(Path(__file__).resolve()), "--render-saved-diagnostic"],
        "started_utc": render_started_utc,
        "finished_utc": render_finished,
        "runtime_seconds": time.monotonic() - started,
        "exit_code": 0 if section_gate else 1,
        "performed_common_boolean": False,
        "source_step_sha256": actual_sha,
        "reloaded_component_brep_count": len(component_shapes),
    }
    if not section_gate:
        report["figure_render_execution"]["failure_reason"] = (
            "corrected exact OCC section did not meet source bbox/wall coverage checks"
        )
    if "execution" not in report:
        started_dt = datetime.fromisoformat(report["created_utc"])
        experiment_runtime = float(report.get("timing", {}).get("whole_run_seconds", 0.0))
        report["execution"] = {
            "argv": [sys.executable, "pgfrag\\experiments\\cad_clipped_cells_spike.py"],
            "raw_sys_argv": ["pgfrag\\experiments\\cad_clipped_cells_spike.py"],
            "shell_command": (
                "& '.\\third_party\\vesselforge-runtime\\venv\\Scripts\\python.exe' "
                "'pgfrag\\experiments\\cad_clipped_cells_spike.py'"
            ),
            "working_directory": str(ROOT),
            "started_utc": report["created_utc"],
            "finished_utc": (started_dt + timedelta(seconds=experiment_runtime)).isoformat(
                timespec="milliseconds"
            ),
            "runtime_seconds": experiment_runtime,
            "runtime_scope": "reported experiment runtime through source-section preparation; excludes interpreter startup and final report serialization",
            "exit_code": 1,
            "python_executable": sys.executable,
            "run_state": report.get("run_state", "INCOMPLETE_COM_OUT"),
        }
    write_json(output_dir / "figure_diagnostics.json", report["figure_diagnostics"])
    write_json(report_path, report)
    write_report_md(output_dir / "report.md", report)
    with (output_dir / "run.log").open("a", encoding="utf-8", newline="\n") as stream:
        stream.write(
            f"{render_finished} AUDIT figure-only rerender: corrected Ax3 plane; "
            f"old_png={superseded_figure.name}; new_png={new_figure.name}; "
            f"section_gate={section_gate}; no Common performed\n"
        )
    logger.info(
        "figure-only diagnostic render complete; section gate=%s; old PNG superseded",
        section_gate,
    )
    return 0 if section_gate else 1


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--render-saved-diagnostic", action="store_true",
        help="rebuild only the saved OCC section PNG; perform no Boolean Common",
    )
    parser.add_argument(
        "--normalize-saved-incomplete", action="store_true",
        help="label the saved interrupted COM_OUT result; perform no geometry",
    )
    args = parser.parse_args()
    if args.normalize_saved_incomplete:
        return normalize_saved_incomplete_report()
    if args.render_saved_diagnostic:
        return render_saved_diagnostic()
    if OUTPUT_DIR.exists() and any(OUTPUT_DIR.iterdir()):
        print(f"Refusing to overwrite existing artifacts: {OUTPUT_DIR}", file=sys.stderr)
        return 2
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    logger = configure_log(OUTPUT_DIR)
    logger.info("one fixed case: 8 mm, <=512 boxes, native OCC Common, single-thread")
    invocation_started = time.monotonic()
    invocation_started_utc = utc_now()
    try:
        report = run(OUTPUT_DIR, logger)
    except Exception as exc:
        logger.exception("experiment aborted before normal report completion")
        report = {
            "schema": "pgfrag.cad_clipped_cells_candidate_spike/1.0",
            "created_utc": utc_now(),
            "status": "FAIL_STATIC_NATIVE_GEOMETRY_CANDIDATE",
            "fatal_exception": {"type": type(exc).__name__, "message": str(exc)},
            "artifacts": {"run_log": str(OUTPUT_DIR / "run.log")},
        }
    exit_code = 0 if report["status"] == "PASS_STATIC_NATIVE_GEOMETRY_CANDIDATE_ONLY" else 1
    report["execution"] = {
        "argv": [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
        "raw_sys_argv": list(sys.argv),
        "started_utc": invocation_started_utc,
        "finished_utc": utc_now(),
        "runtime_seconds": time.monotonic() - invocation_started,
        "runtime_scope": "whole Python process invocation",
        "exit_code": exit_code,
        "python_executable": sys.executable,
        "run_state": report.get("run_state", "COMPLETE"),
    }
    write_json(OUTPUT_DIR / "report.json", report)
    if "source" in report:
        write_report_md(OUTPUT_DIR / "report.md", report)
    logger.info(
        "process exit=%d runtime=%.3fs",
        exit_code, report["execution"]["runtime_seconds"],
    )
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
