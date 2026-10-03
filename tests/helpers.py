from __future__ import annotations

import hashlib
from contextlib import contextmanager
import json
from itertools import count
from pathlib import Path
from typing import Any

import numpy as np


_TEMP_CASES = count()


@contextmanager
def temporary_directory():
    """Use pre-created writable test scratch beneath the authorized PG audit area."""
    base = (Path(__file__).resolve().parents[2] / "relicweave-synlab" / "tmp" /
            "review_20261003" / "implementation" / "pg" / "test_cases")
    case = base / f"case_{next(_TEMP_CASES):03d}"
    if not case.is_dir():
        raise RuntimeError(f"test scratch directory was not prepared: {case}")
    yield str(case)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def make_asset(root: Path, *, status: str = "PASS", binding_id: str = "mat-1",
               material_xyz_offset: float = 0.0) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    for sub in ("geometry", "fields", "validation"):
        (root / sub).mkdir(exist_ok=True)
    vertices = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float64)
    faces = np.array([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]], dtype=np.int32)
    region = np.array([1, 2, 3, 255], dtype=np.uint8)
    np.savez(root / "geometry" / "surface.npz", vertices_mm=vertices, faces=faces)
    np.savez(root / "fields" / "surface.npz", face_region=region,
             albedo_linear_rgb=np.full((4, 3), 0.5, dtype=np.float32),
             roughness=np.full(4, 0.75, dtype=np.float32),
             vertices_material_xyz_mm=vertices + material_xyz_offset)
    report = {"scope": "geometry_only", "status": status, "tests": {"closed": status == "PASS"}}
    (root / "validation" / "geometry.json").write_text(json.dumps(report), encoding="utf-8")
    payloads = [root / "geometry" / "surface.npz", root / "fields" / "surface.npz", root / "validation" / "geometry.json"]
    files = {p.relative_to(root).as_posix(): {"bytes": p.stat().st_size, "sha256": sha256(p)} for p in payloads}
    manifest: dict[str, Any] = {
        "schema": "forge.asset.v1",
        "units": {"length": "mm", "time": "s", "mass": "kg"},
        "frame": "right_handed_z_up",
        "geometry": {"surface_npz": "geometry/surface.npz"},
        "fields": {"surface_fields_npz": "fields/surface.npz"},
        "material": {"id": "mat-1", "mechanics": {"rho_kg_m3": 2100.0, "E_pa": 5e9, "nu": 0.25, "Gc_j_m2": 45.0},
                     "calibration_status": "UNVALIDATED_PROTOTYPE"},
        "ornament": {"kind": "painted", "recipe": {"material_binding_id": binding_id, "geometry_relief_mm": 0.0,
                     "historical_style_validation": "NOT_VALIDATED", "target_region_codes": [1]}},
        "files": files,
        "provenance": {"producer": "conformance-fixture", "version": "1"},
        "validation": {"report": "validation/geometry.json", "status": status, "scope": "geometry_only"},
    }
    (root / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return root


def base_job() -> dict[str, Any]:
    return {
        "schema": "pgfrag.job.v1", "job_id": "test-job-01", "backend_id": "peridynamics_jl_headless_0_5_4",
        "initial_condition": {"kind": "uniform_box_bar_fixture", "dimensions_mm": [6, 6, 16], "dx_mm": 2, "horizon_mm": 6},
        "material": {"model": "OSB", "rho_kg_m3": 2100, "E_pa": 5e9, "nu": 0.24, "Gc_j_m2": 45,
                     "calibration_status": "UNVALIDATED_PROTOTYPE", "provenance": "unvalidated test values"},
        "loading": {"kind": "symmetric_axial_velocity", "speed_m_s": 0.0001},
        "solver": {"steps": 5, "dt_s": None, "safety_factor": 0.8},
        "resources": {"profile": "local_smoke", "max_points": 1000, "max_steps": 500, "threads": 1, "timeout_seconds": 120},
        "output_dir": "output",
    }
