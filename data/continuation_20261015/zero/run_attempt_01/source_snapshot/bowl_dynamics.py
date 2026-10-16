"""Bounded native bond-based peridynamic dynamics for an analytic spherical bowl.

This module prepares immutable run attempts and packages the native Julia
trajectory.  It intentionally does not launch a solve unless ``run_prepared``
is called (or the module CLI receives ``--execute``).
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from typing import Any, Iterable

import numpy as np


SCHEMA = "pgfrag.native_bowl_result.v1"
CASE_SCHEMA = "pgfrag.native_bowl_case.v1"
REPO_ROOT = Path(__file__).resolve().parents[2]
JULIA_PROJECT = REPO_ROOT / "src" / "pgfrag" / "resources" / "julia"
JULIA_RUNNER = JULIA_PROJECT / "run_bowl_dynamics.jl"
DATA_ROOT = REPO_ROOT / "data" / "continuation_20261015"
MAX_POINTS = 12_000
MAX_STEPS = 1_000
MAX_PAIRS = 2_000_000
MAX_RUNTIME_S = 120.0
MAX_OUTPUT_BYTES = 256 * 1024 * 1024
FROZEN_RADII_M = (12e-3, 8e-3, (0.0, 0.0, 12e-3))
FROZEN_MECHANICS = {
    "rho_kg_m3": 2300.0,
    "E_pa": 70.0e9,
    "nu": 0.25,
    "Gc_j_m2": 7.6,
}
FROZEN_SPACING_M = (1e-3, 0.8e-3)
FROZEN_HORIZON_M = 2.5e-3


@dataclass(frozen=True)
class BowlGeometry:
    outer_radius_m: float
    inner_radius_m: float
    center_m: tuple[float, float, float]

    @property
    def analytic_volume_m3(self) -> float:
        return 2.0 * math.pi / 3.0 * (self.outer_radius_m**3 - self.inner_radius_m**3)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_json(path: Path) -> tuple[dict[str, Any], bytes]:
    raw = path.read_bytes()
    value = json.loads(raw.decode("utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected a JSON object: {path}")
    return value, raw


def _manifest_path(asset: str | os.PathLike[str]) -> Path:
    candidate = Path(asset).expanduser().resolve(strict=True)
    if candidate.is_dir():
        candidate = candidate / "manifest.json"
    if candidate.name != "manifest.json" or not candidate.is_file():
        raise ValueError("asset input must be a forge.asset.v1 package or its manifest.json")
    return candidate


def _load_asset_spec(
    asset: str | os.PathLike[str],
) -> tuple[Path, dict[str, Any], BowlGeometry, str, dict[str, float]]:
    manifest_path = _manifest_path(asset)
    manifest, raw = _read_json(manifest_path)
    if manifest.get("schema") != "forge.asset.v1":
        raise ValueError("source manifest schema must be forge.asset.v1")
    if manifest.get("units") != {"length": "mm", "time": "s", "mass": "kg"}:
        raise ValueError("source asset units are not canonical")
    if manifest.get("frame") != "right_handed_z_up":
        raise ValueError("source asset frame must be right_handed_z_up")
    validation = manifest.get("validation", {})
    if validation.get("scope") != "geometry_only" or validation.get("status") != "PASS":
        raise ValueError("source bowl must carry a PASS geometry_only validation declaration")
    provenance = manifest.get("provenance", {})
    spec = provenance.get("spec")
    if not isinstance(spec, dict):
        raise ValueError("source manifest provenance.spec is required")
    shape = spec.get("shape", spec.get("prototype_form"))
    if shape not in {"concentric_lower_hemisphere", "concentric_lower_hemisphere_bowl"}:
        raise ValueError("source asset is not the declared concentric lower-hemisphere bowl")
    try:
        outer_mm = float(spec["outer_radius_mm"])
        inner_mm = float(spec["inner_radius_mm"])
        center_mm = tuple(float(x) for x in spec["center_xyz_mm"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("source manifest has incomplete spherical-bowl parameters") from exc
    if len(center_mm) != 3 or not all(math.isfinite(x) for x in (outer_mm, inner_mm, *center_mm)):
        raise ValueError("source bowl parameters must be finite")

    material_record = manifest.get("material")
    if not isinstance(material_record, dict) or not isinstance(material_record.get("mechanics"), dict):
        raise ValueError("source manifest material.mechanics is required")
    source_mechanics: dict[str, float] = {}
    for key in FROZEN_MECHANICS:
        try:
            source_mechanics[key] = float(material_record["mechanics"][key])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"source manifest material.mechanics.{key} is required") from exc
        if not math.isfinite(source_mechanics[key]) or source_mechanics[key] != FROZEN_MECHANICS[key]:
            raise ValueError(f"source material {key} does not match the frozen illustrative case")
    if material_record.get("calibration_status") != "UNVALIDATED_PROTOTYPE":
        raise ValueError("source bowl mechanics must remain UNVALIDATED_PROTOTYPE")

    # Verify every manifest-listed payload by streaming bytes. Do not load the
    # surface/field NPZ arrays into memory; the solver uses the analytic spec.
    files = manifest.get("files")
    if not isinstance(files, dict) or not files:
        raise ValueError("source manifest.files must enumerate payload hashes")
    package_root = manifest_path.parent.resolve(strict=True)
    for relative, record in sorted(files.items()):
        if not isinstance(relative, str) or "\\" in relative or Path(relative).is_absolute() or ".." in Path(relative).parts:
            raise ValueError(f"unsafe source payload path: {relative!r}")
        payload = (package_root / relative).resolve(strict=True)
        payload.relative_to(package_root)
        if not payload.is_file() or not isinstance(record, dict):
            raise ValueError(f"invalid source payload record: {relative!r}")
        if payload.stat().st_size != int(record.get("bytes", -1)) or sha256_file(payload) != str(record.get("sha256", "")).lower():
            raise ValueError(f"source payload digest/size mismatch: {relative}")
    geometry = BowlGeometry(outer_mm * 1e-3, inner_mm * 1e-3, tuple(x * 1e-3 for x in center_mm))
    return manifest_path, manifest, geometry, hashlib.sha256(raw).hexdigest(), source_mechanics


def generate_centroid_lattice(
    geometry: BowlGeometry,
    spacing_m: float = 1e-3,
    *,
    max_points: int = MAX_POINTS,
) -> np.ndarray:
    """Mirror VesselForge's deterministic cell-centroid spherical-bowl sampler."""
    if not math.isfinite(spacing_m) or spacing_m <= 0:
        raise ValueError("spacing_m must be finite and positive")
    cx, cy, cz = geometry.center_m
    r = geometry.outer_radius_m
    low = np.asarray((cx - r, cy - r, cz - r), dtype=np.float64)
    high = np.asarray((cx + r, cy + r, cz), dtype=np.float64)
    axes = []
    for lo, hi in zip(low, high):
        count = int(math.floor((hi - lo) / spacing_m + 1e-12))
        axes.append(lo + (np.arange(count, dtype=np.float64) + 0.5) * spacing_m)
    xx, yy, zz = np.meshgrid(*axes, indexing="ij")
    candidates = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
    centered = candidates - np.asarray(geometry.center_m, dtype=np.float64)
    radial = np.linalg.norm(centered, axis=1)
    points = np.asarray(candidates[(radial <= geometry.outer_radius_m) & (radial >= geometry.inner_radius_m)], dtype=np.float64)
    if len(points) == 0 or len(points) > max_points:
        raise ValueError(f"bowl lattice point count outside budget: {len(points)}")
    return points


def segment_minimum_radius_m(a: Iterable[float], b: Iterable[float], center: Iterable[float]) -> float:
    """Return the exact closest radius of a closed segment to a sphere center."""
    start = np.asarray(tuple(a), dtype=np.float64)
    end = np.asarray(tuple(b), dtype=np.float64)
    origin = np.asarray(tuple(center), dtype=np.float64)
    if start.shape != (3,) or end.shape != (3,) or origin.shape != (3,):
        raise ValueError("segment endpoints and sphere center must be three-vectors")
    direction = end - start
    denominator = float(np.dot(direction, direction))
    t = 0.0 if denominator == 0.0 else float(np.clip(-np.dot(start - origin, direction) / denominator, 0.0, 1.0))
    return float(np.linalg.norm(start + t * direction - origin))


def bond_crosses_inner_void(a: Iterable[float], b: Iterable[float], geometry: BowlGeometry) -> bool:
    return segment_minimum_radius_m(a, b, geometry.center_m) < geometry.inner_radius_m


def independent_pair_energy_joule(
    x_m: np.ndarray,
    bond_i: np.ndarray,
    bond_j: np.ndarray,
    rest_length_m: np.ndarray,
    volume_m3: np.ndarray,
    bond_constant_n_m4: float,
    active: np.ndarray | None = None,
) -> float:
    """Independent BB pair potential using the frozen discrete bond formula."""
    x = np.asarray(x_m, dtype=np.float64)
    i = np.asarray(bond_i, dtype=np.int64)
    j = np.asarray(bond_j, dtype=np.int64)
    L = np.asarray(rest_length_m, dtype=np.float64)
    volume = np.asarray(volume_m3, dtype=np.float64)
    if x.ndim != 2 or x.shape[1] != 3 or i.shape != j.shape or i.shape != L.shape:
        raise ValueError("point and bond arrays have inconsistent shapes")
    if len(i) and (np.min(i) < 0 or np.max(j) >= len(x)):
        raise ValueError("bond endpoint is outside the point array")
    stretch = (np.linalg.norm(x[j] - x[i], axis=1) - L) / L
    pair_u = 0.5 * bond_constant_n_m4 * volume[i] * volume[j] * L * stretch**2
    if active is not None:
        mask = np.asarray(active, dtype=bool)
        if mask.shape != pair_u.shape:
            raise ValueError("active bond mask shape disagrees with bond roster")
        pair_u = pair_u[mask]
    return float(np.sum(pair_u, dtype=np.float64))


def rigid_transform(points: np.ndarray, rotation: np.ndarray, translation: np.ndarray) -> np.ndarray:
    """Apply a proper rigid transform to row-vector coordinates."""
    xyz = np.asarray(points, dtype=np.float64)
    R = np.asarray(rotation, dtype=np.float64)
    t = np.asarray(translation, dtype=np.float64)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or R.shape != (3, 3) or t.shape != (3,):
        raise ValueError("invalid rigid-transform shapes")
    return xyz @ R.T + t


def _write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")


def _toml_quote(text: str) -> str:
    return json.dumps(text)


def _write_case_toml(path: Path, case: dict[str, Any]) -> None:
    g = case["geometry"]
    m = case["material"]
    load = case["loading"]
    lines = [
        f"schema = {_toml_quote(CASE_SCHEMA)}",
        f"run_id = {_toml_quote(case['run_id'])}",
        f"points_csv = {_toml_quote('points.csv')}",
        f"source_forge_manifest_sha256 = {_toml_quote(case['source_forge_manifest_sha256'])}",
        f"source_manifest_path = {_toml_quote(case['source_manifest_path'])}",
        "",
        "[geometry]",
        f"outer_radius_m = {g['outer_radius_m']:.17g}",
        f"inner_radius_m = {g['inner_radius_m']:.17g}",
        f"center_x_m = {g['center_m'][0]:.17g}",
        f"center_y_m = {g['center_m'][1]:.17g}",
        f"center_z_m = {g['center_m'][2]:.17g}",
        f"h_m = {g['h_m']:.17g}",
        f"spacing_m = {g['h_m']:.17g}",
        f"delta_m = {g['delta_m']:.17g}",
        f"n_points = {g['n_points']}",
        f"analytic_volume_m3 = {g['analytic_volume_m3']:.17g}",
        "",
        "[material]",
        f"rho_kg_m3 = {m['rho_kg_m3']:.17g}",
        f"E_pa = {m['E_pa']:.17g}",
        f"nu = {m['nu']:.17g}",
        f"Gc_j_m2 = {m['Gc_j_m2']:.17g}",
        f"calibration_status = {_toml_quote(m['calibration_status'])}",
        f"total_mass_kg = {m['total_mass_kg']:.17g}",
        "",
        "[loading]",
        f"mode = {_toml_quote(load['mode'])}",
        f"alpha_s_1 = {load['alpha_s_1']:.17g}",
        f"momentum_constraint = {_toml_quote(load['momentum_constraint'])}",
        "",
        "[solver]",
        f"steps = {case['solver']['steps']}",
        f"safety_factor = {case['solver']['safety_factor']:.17g}",
        f"time_resolution_case = {_toml_quote(case['solver']['time_resolution_case'])}",
        f"threads = {case['solver']['threads']}",
        f"max_runtime_s = {MAX_RUNTIME_S:.17g}",
        f"max_output_bytes = {MAX_OUTPUT_BYTES}",
        f"max_points = {MAX_POINTS}",
        f"max_steps = {MAX_STEPS}",
        f"max_pairs = {MAX_PAIRS}",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def _next_attempt(root: Path) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    for index in range(1, 10_000):
        candidate = root / f"run_attempt_{index:02d}"
        try:
            candidate.mkdir()
        except FileExistsError:
            continue
        return candidate
    raise RuntimeError("no unused run-attempt directory is available")


def _create_source_snapshot(run_dir: Path, manifest_path: Path) -> list[dict[str, Any]]:
    """Copy the exact producer, native runner, locked project, vendor, and asset bytes."""
    snapshot = run_dir / "source_snapshot"
    snapshot.mkdir()
    fixed = (
        (Path(__file__).resolve(), snapshot / "bowl_dynamics.py"),
        (JULIA_RUNNER, snapshot / "run_bowl_dynamics.jl"),
        (JULIA_PROJECT / "Project.toml", snapshot / "Project.toml"),
        (JULIA_PROJECT / "Manifest.toml", snapshot / "Manifest.toml"),
        (JULIA_PROJECT / "backend_provenance.json", snapshot / "backend_provenance.json"),
    )
    for source, destination in fixed:
        if not source.is_file():
            raise FileNotFoundError(f"required source snapshot input is missing: {source}")
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
    vendor_root = JULIA_PROJECT / "vendor" / "Peridynamics.jl"
    if not vendor_root.is_dir():
        raise FileNotFoundError(f"pinned Peridynamics.jl source is missing: {vendor_root}")
    shutil.copytree(vendor_root, snapshot / "vendor" / "Peridynamics.jl")
    shutil.copy2(manifest_path, snapshot / "forge_manifest.json")
    records = []
    for source in sorted(snapshot.rglob("*")):
        if source.is_file():
            records.append({
                "path": source.relative_to(run_dir).as_posix(),
                "sha256": sha256_file(source),
                "bytes": source.stat().st_size,
            })
    if not records:
        raise RuntimeError("source snapshot unexpectedly contains no files")
    return records


def _verify_source_snapshot(run_dir: Path, *, verify_canonical: bool = True) -> list[dict[str, Any]]:
    """Fail closed if source bytes used by a prepared attempt have changed."""
    prepared, _ = _read_json(run_dir / "argv.json")
    records = prepared.get("source_receipts")
    if not isinstance(records, list) or not records:
        raise ValueError("prepared attempt has no source byte receipts")
    by_path: dict[str, dict[str, Any]] = {}
    for record in records:
        if not isinstance(record, dict) or not isinstance(record.get("path"), str):
            raise ValueError("malformed source byte receipt")
        relative = record["path"]
        if relative in by_path or not relative.startswith("source_snapshot/"):
            raise ValueError(f"unsafe or duplicate snapshot receipt path: {relative}")
        path = run_dir / Path(relative)
        if not path.is_file() or path.stat().st_size != int(record.get("bytes", -1)) \
                or sha256_file(path) != str(record.get("sha256", "")).lower():
            raise ValueError(f"pinned source snapshot changed: {relative}")
        by_path[relative] = record
    manifest_path = run_dir / "source_snapshot" / "forge_manifest.json"
    if sha256_file(manifest_path) != _read_json(run_dir / "case.json")[0]["source_forge_manifest_sha256"]:
        raise ValueError("snapshotted Forge manifest no longer matches the case hash")
    if not verify_canonical:
        return records

    canonical: dict[str, Path] = {
        "source_snapshot/bowl_dynamics.py": Path(__file__).resolve(),
        "source_snapshot/run_bowl_dynamics.jl": JULIA_RUNNER,
        "source_snapshot/Project.toml": JULIA_PROJECT / "Project.toml",
        "source_snapshot/Manifest.toml": JULIA_PROJECT / "Manifest.toml",
        "source_snapshot/backend_provenance.json": JULIA_PROJECT / "backend_provenance.json",
    }
    vendor_prefix = "source_snapshot/vendor/Peridynamics.jl/"
    for relative in by_path:
        if relative.startswith(vendor_prefix):
            suffix = relative[len(vendor_prefix):]
            canonical[relative] = JULIA_PROJECT / "vendor" / "Peridynamics.jl" / Path(suffix)
    for relative, source in canonical.items():
        record = by_path.get(relative)
        if record is None or not source.is_file() or source.stat().st_size != int(record["bytes"]) \
                or sha256_file(source) != record["sha256"]:
            raise ValueError(f"canonical pinned source differs from prepared snapshot: {relative}")
    live_vendor_files = {
        path.relative_to(JULIA_PROJECT / "vendor" / "Peridynamics.jl").as_posix()
        for path in (JULIA_PROJECT / "vendor" / "Peridynamics.jl").rglob("*") if path.is_file()
    }
    snapshot_vendor_files = {
        path[len(vendor_prefix):] for path in by_path if path.startswith(vendor_prefix)
    }
    if live_vendor_files != snapshot_vendor_files:
        raise ValueError("pinned Peridynamics.jl source file inventory changed")
    return records


def prepare_run(
    asset: str | os.PathLike[str],
    *,
    runs_root: str | os.PathLike[str] = DATA_ROOT,
    alpha_s_1: float = 0.0,
    steps: int = 400,
    safety_factor: float = 0.5,
    spacing_m: float = 1e-3,
) -> Path:
    """Create a new source-pinned run attempt without launching Julia."""
    if isinstance(steps, bool) or int(steps) not in {400, 800}:
        raise ValueError("the frozen dynamics cases use 400 baseline or 800 half-dt steps")
    if not math.isfinite(float(safety_factor)):
        raise ValueError("safety_factor must be finite")
    expected_safety = 0.5 if int(steps) == 400 else 0.25
    if float(safety_factor) != expected_safety:
        raise ValueError(f"{int(steps)} steps require the frozen safety_factor={expected_safety}")
    if not math.isfinite(float(alpha_s_1)):
        raise ValueError("loading parameter must be finite")
    if alpha_s_1 not in {0.0, 100.0, 300.0}:
        raise ValueError("alpha_s_1 must be one of the frozen values 0, 100, or 300")
    if float(spacing_m) not in FROZEN_SPACING_M:
        raise ValueError(f"spacing_m must be one of the frozen spatial cases {FROZEN_SPACING_M}")
    manifest_path, manifest, geometry, manifest_hash, source_mechanics = _load_asset_spec(asset)
    if (geometry.outer_radius_m, geometry.inner_radius_m, geometry.center_m) != FROZEN_RADII_M:
        raise ValueError("this frozen run case requires Ro=12 mm, Ri=8 mm, center=(0,0,12) mm")
    h_m, delta_m = float(spacing_m), FROZEN_HORIZON_M
    points_m = generate_centroid_lattice(geometry, h_m)
    analytic_volume = geometry.analytic_volume_m3
    volume_m3 = len(points_m) * h_m**3
    relative_volume_error = abs(volume_m3 - analytic_volume) / analytic_volume
    if len(points_m) > MAX_POINTS:
        raise ValueError("lattice exceeds native point limit")

    run_dir = _next_attempt(Path(runs_root).expanduser().resolve())
    time_case = "baseline" if int(steps) == 400 else "half_dt_sensitivity"
    spatial_case = "h1mm_fixed_delta" if h_m == 1e-3 else "h0p8mm_fixed_delta"
    run_id = f"native_bowl_{spatial_case}_{int(alpha_s_1)}s-1_{steps}steps_{run_dir.name}"
    point_mass = float(source_mechanics["rho_kg_m3"] * h_m**3)
    case = {
        "schema": CASE_SCHEMA,
        "run_id": run_id,
        "source_forge_manifest_sha256": manifest_hash,
        "source_manifest_path": str(manifest_path),
        "source_forge_manifest": manifest,
        "geometry": {
            "shape": "concentric_lower_hemisphere_shell",
            "outer_radius_m": geometry.outer_radius_m,
            "inner_radius_m": geometry.inner_radius_m,
            "center_m": list(geometry.center_m),
            "h_m": h_m,
            "spacing_m": h_m,
            "delta_m": delta_m,
            "n_points": len(points_m),
            "analytic_volume_m3": analytic_volume,
            "lattice_volume_m3": volume_m3,
            "relative_volume_error": relative_volume_error,
            "volume_estimator": "centroid_count_times_h_cubed; boundary cells uncut",
            "source_point_sampler": "VesselForge spherical-bowl generate_solver_points centroid-lattice algorithm mirrored in PG-Frag",
        },
        "material": {
            **source_mechanics,
            "calibration_status": manifest["material"]["calibration_status"],
            "source_binding": "exact values copied from forge.asset.v1 material.mechanics",
            "point_mass_kg": point_mass,
            "total_mass_kg": point_mass * len(points_m),
            "point_mass_model": "rho_kg_m3 * spacing_m^3 for every selected centroid; no uniform rescaling",
        },
        "loading": {
            "mode": "zero" if alpha_s_1 == 0.0 else "imposed_strain_rate_release",
            "alpha_s_1": float(alpha_s_1),
            "rate_tensor_s_1": [[alpha_s_1, 0.0, 0.0], [0.0, -0.25 * alpha_s_1, 0.0], [0.0, 0.0, 0.0]],
            "momentum_constraint": "subtract mass-mean translation, then solve inertia tensor for zero angular momentum",
            "description": "freely evolving initial state; no grip, precrack, impact, or drop",
        },
        "solver": {
            "name": "Peridynamics.jl 0.5.4 VelocityVerlet, BBMaterial{NoCorrection}, CriticalStretch",
            "steps": int(steps),
            "safety_factor": float(safety_factor),
            "time_resolution_case": time_case,
            "spatial_sensitivity_case": spatial_case,
            "threads": 1,
            "max_runtime_s": MAX_RUNTIME_S,
            "max_output_bytes": MAX_OUTPUT_BYTES,
        },
    }
    # The exact text format is also part of the source input receipt.
    with (run_dir / "points.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(("x_m", "y_m", "z_m"))
        writer.writerows(points_m.tolist())
    _write_case_toml(run_dir / "case.toml", case)
    source_records = _create_source_snapshot(run_dir, manifest_path)
    _write_json(run_dir / "case.json", case)
    _write_json(run_dir / "argv.json", {
        "argv": ["<julia-executable>", f"--project={JULIA_PROJECT}", str(run_dir / "source_snapshot" / "run_bowl_dynamics.jl"), str(run_dir / "case.toml"), str(run_dir)],
        "cwd": str(REPO_ROOT),
        "run_id": run_id,
        "prepared_only": True,
        "argv_sha256": hashlib.sha256(json.dumps(
            ["<julia-executable>", f"--project={JULIA_PROJECT}", str(run_dir / "source_snapshot" / "run_bowl_dynamics.jl"), str(run_dir / "case.toml"), str(run_dir)],
            ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest(),
        "source_receipts": source_records,
    })
    _write_json(run_dir / "prepared.json", {
        "schema": CASE_SCHEMA,
        "status": "SOURCE_READY_NOT_RUN",
        "run_id": run_id,
        "source_forge_manifest_sha256": manifest_hash,
        "source_receipts": source_records,
        "geometry": case["geometry"],
        "material": case["material"],
        "case_sha256": sha256_file(run_dir / "case.json"),
        "limitations": ["geometry/material are prototypes; solve not launched by preparation"],
    })
    return run_dir


def _parse_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def _load_raw_states(raw_dir: Path, n_steps: int, n_points: int) -> dict[str, np.ndarray]:
    # Julia writes vec(Array(S,N,3)): state is the fastest-varying index.
    # Explicit Fortran order is part of the cross-language byte contract.
    n_states = n_steps + 1
    arrays = {
        "position_m": np.fromfile(raw_dir / "position.f64le", dtype="<f8").reshape(n_states, n_points, 3, order="F"),
        "velocity_mps": np.fromfile(raw_dir / "velocity.f64le", dtype="<f8").reshape(n_states, n_points, 3, order="F"),
        "damage": np.fromfile(raw_dir / "damage.f32le", dtype="<f4").reshape(n_states, n_points, order="F"),
        "component_label": np.fromfile(raw_dir / "component.i32le", dtype="<i4").reshape(n_states, n_points, order="F"),
        "native_sed_J_m3": np.fromfile(raw_dir / "native_sed.f32le", dtype="<f4").reshape(n_states, n_points, order="F"),
    }
    checks = _parse_csv(raw_dir / "checkpoints.csv")
    if checks:
        arrays["checkpoint_step"] = np.asarray([int(row["step"]) for row in checks], dtype=np.int32)
        arrays["checkpoint_acceleration_mps2"] = np.fromfile(
            raw_dir / "checkpoint_acceleration.f64le", dtype="<f8"
        ).reshape(len(checks), n_points, 3, order="F")
        arrays["checkpoint_internal_force_N"] = np.fromfile(
            raw_dir / "checkpoint_force.f64le", dtype="<f8"
        ).reshape(len(checks), n_points, 3, order="F")
    return arrays


def package_run(run_dir: str | os.PathLike[str]) -> dict[str, Any]:
    """Package one completed Julia attempt into its immutable result contract."""
    root = Path(run_dir).resolve(strict=True)
    _verify_source_snapshot(root)
    case, _ = _read_json(root / "case.json")
    native_status = _read_json(root / "native_status.json")[0]
    if native_status.get("status") != "NATIVE_COMPLETE":
        raise RuntimeError(f"native solve did not complete: {native_status.get('status')}")
    n_points = int(case["geometry"]["n_points"])
    n_steps = int(native_status["steps"])
    if not 0 < n_steps <= MAX_STEPS or n_points > MAX_POINTS:
        raise ValueError("native result exceeds the frozen point/step budget")
    raw_dir = root / "native_raw"
    states = _load_raw_states(raw_dir, n_steps, n_points)
    point_ids = np.arange(n_points, dtype=np.int32)
    volumes = np.full(n_points, float(case["geometry"]["h_m"]) ** 3, dtype=np.float64)
    rho = float(case["material"]["rho_kg_m3"])
    point_mass = rho * volumes
    if len(point_mass) != n_points or np.any(point_mass <= 0) or not np.isfinite(point_mass).all():
        raise ValueError("point masses do not cover every native point with positive finite values")

    bond_rows = _parse_csv(raw_dir / "bonds.csv")
    if len(bond_rows) > MAX_PAIRS:
        raise ValueError("native canonical pair count exceeds the frozen limit")
    bond_i = np.asarray([int(row["bond_i"]) for row in bond_rows], dtype=np.int32)
    bond_j = np.asarray([int(row["bond_j"]) for row in bond_rows], dtype=np.int32)
    rest_length = np.asarray([float(row["rest_length_m"]) for row in bond_rows], dtype=np.float64)
    epsilon_c = np.asarray([float(row["epsilon_c"]) for row in bond_rows], dtype=np.float64)
    n_bonds = len(bond_rows)
    if n_bonds and (np.any(bond_i >= bond_j) or np.any(bond_i < 0) or np.any(bond_j >= n_points)):
        raise ValueError("canonical bond endpoints are not sorted valid 0-based IDs")

    event_rows = _parse_csv(raw_dir / "events.csv")
    event_bond = np.asarray([int(row["bond_index"]) for row in event_rows], dtype=np.int32)
    event_step = np.asarray([int(row["step"]) for row in event_rows], dtype=np.int32)
    event_time = np.asarray([float(row["time_s"]) for row in event_rows], dtype=np.float64)
    event_stretch = np.asarray([float(row["stretch"]) for row in event_rows], dtype=np.float64)
    event_energy = np.asarray([float(row["pair_energy_J"]) for row in event_rows], dtype=np.float64)
    first_step = np.full(n_bonds, -1, dtype=np.int32)
    first_time = np.full(n_bonds, -1.0, dtype=np.float64)
    first_stretch = np.full(n_bonds, -1.0, dtype=np.float64)
    if len(event_bond):
        if np.any(event_bond < 0) or np.any(event_bond >= n_bonds) or np.any(event_step <= 0) or np.any(event_step > n_steps):
            raise ValueError("native break-event row lies outside its canonical pair/state range")
        if len(np.unique(event_bond)) != len(event_bond):
            raise ValueError("native event table contains duplicate first-break events")
        first_step[event_bond] = event_step
        first_time[event_bond] = event_time
        first_stretch[event_bond] = event_stretch
    active = np.ones((n_steps + 1, n_bonds), dtype=np.uint8)
    for state_idx, step in enumerate(range(n_steps + 1)):
        active[state_idx, first_step >= 0] = (first_step[first_step >= 0] > step).astype(np.uint8)
    active_packed = np.packbits(active, axis=1, bitorder="little")
    native_active_packed = np.fromfile(raw_dir / "bond_active_packed.u8", dtype=np.uint8).reshape(
        n_steps + 1, math.ceil(n_bonds / 8), order="F"
    )
    if not np.array_equal(native_active_packed, active_packed):
        raise ValueError("native active-bond history does not match the irreversible first-break event log")
    active_packed = native_active_packed

    ledger_rows = _parse_csv(raw_dir / "energy.csv")
    if len(ledger_rows) != n_steps + 1:
        raise ValueError("native energy ledger does not contain every state from 0 through n_steps")
    step = np.asarray([int(row["step"]) for row in ledger_rows], dtype=np.int32)
    time_s = np.asarray([float(row["time_s"]) for row in ledger_rows], dtype=np.float64)
    energy = {key: np.asarray([float(row[key]) for row in ledger_rows], dtype=np.float64) for key in (
        "kinetic_energy_J", "pair_potential_energy_J", "break_dissipation_cumulative_J",
        "native_strain_energy_J", "total_force_x_N", "total_force_y_N", "total_force_z_N",
        "momentum_x_kg_m_s", "momentum_y_kg_m_s", "momentum_z_kg_m_s",
        "angular_momentum_x_kg_m2_s", "angular_momentum_y_kg_m2_s", "angular_momentum_z_kg_m2_s",
    )}
    states.update({
        "point_id": point_ids,
        "x_ref_m": np.asarray([[float(row["x_m"]), float(row["y_m"]), float(row["z_m"])] for row in _parse_csv(root / "points.csv")], dtype=np.float64),
        "volume_m3": volumes,
        "point_mass_kg": point_mass,
        "step": step,
        "time_s": time_s,
        "bond_i": bond_i,
        "bond_j": bond_j,
        "bond_rest_length_m": rest_length,
        "bond_epsilon_c": epsilon_c,
        "bond_first_break_step": first_step,
        "bond_first_break_time_s": first_time,
        "bond_first_break_stretch": first_stretch,
        "bond_active_packed": active_packed,
        "component_label": states["component_label"].astype(np.int32, copy=False),
        "event_bond_index": event_bond,
        "event_step": event_step,
        "event_time_s": event_time,
        "event_stretch": event_stretch,
        "event_pair_energy_J": event_energy,
        **energy,
    })
    if states["x_ref_m"].shape != (n_points, 3):
        raise ValueError("point reference CSV does not cover all points")
    if states["position_m"].shape != (n_steps + 1, n_points, 3) or states["velocity_mps"].shape != (n_steps + 1, n_points, 3):
        raise ValueError("native point state arrays have an unexpected shape")
    if n_bonds and states["bond_active_packed"].shape != (n_steps + 1, math.ceil(n_bonds / 8)):
        raise ValueError("packed active state array has an unexpected shape")

    arrays_path = root / "arrays.npz"
    np.savez_compressed(arrays_path, **states)
    peak_materialized_bytes = sum(path.stat().st_size for path in root.rglob("*") if path.is_file())
    if peak_materialized_bytes > MAX_OUTPUT_BYTES:
        arrays_path.unlink(missing_ok=True)
        raise ValueError(f"raw plus packaged trajectory exceeded materialized output budget: {peak_materialized_bytes} bytes")
    transient_files = (
        "position.f64le", "velocity.f64le", "damage.f32le", "component.i32le",
        "native_sed.f32le", "bond_active_packed.u8", "checkpoint_acceleration.f64le",
        "checkpoint_force.f64le",
    )
    transient_bytes = sum((raw_dir / name).stat().st_size for name in transient_files if (raw_dir / name).is_file())
    for name in transient_files:
        candidate = raw_dir / name
        if candidate.exists():
            candidate.unlink()
    outputs = []
    output_paths = (
        root / "arrays.npz", root / "points.csv", root / "case.toml", root / "case.json",
        root / "native_status.json", root / "native_preflight.json", raw_dir / "bonds.csv", raw_dir / "excluded_bonds.csv",
        raw_dir / "events.csv", raw_dir / "energy.csv", raw_dir / "checkpoints.csv",
    )
    for candidate in output_paths:
        if candidate.is_file():
            outputs.append({"path": candidate.relative_to(root).as_posix(), "sha256": sha256_file(candidate), "bytes": candidate.stat().st_size})
    source_receipts = _read_json(root / "argv.json")[0]["source_receipts"]
    source_manifest_hash = case["source_forge_manifest_sha256"]
    total_final_bytes = sum(path.stat().st_size for path in root.rglob("*") if path.is_file())
    persistent_bytes = total_final_bytes
    if persistent_bytes > MAX_OUTPUT_BYTES:
        arrays_path.unlink(missing_ok=True)
        raise ValueError(f"run attempt exceeded persistent output budget: {persistent_bytes} bytes")

    manifest = case["source_forge_manifest"]
    results = {
        "point_count": n_points,
        "canonical_initial_pair_count": n_bonds,
        "initial_geometric_exclusion_count": int(native_status.get("initial_geometric_exclusion_count", -1)),
        "broken_pair_count": len(event_bond),
        "analytic_volume_m3": float(case["geometry"]["analytic_volume_m3"]),
        "lattice_volume_m3": float(np.sum(volumes, dtype=np.float64)),
        "relative_volume_error": abs(float(np.sum(volumes)) - float(case["geometry"]["analytic_volume_m3"])) / float(case["geometry"]["analytic_volume_m3"]),
        "reported_point_mass_sum_kg": float(np.sum(point_mass, dtype=np.float64)),
        "initial_kinetic_energy_J": float(energy["kinetic_energy_J"][0]),
        "final_kinetic_energy_J": float(energy["kinetic_energy_J"][-1]),
        "final_pair_potential_energy_J": float(energy["pair_potential_energy_J"][-1]),
        "final_break_dissipation_cumulative_J": float(energy["break_dissipation_cumulative_J"][-1]),
        "persistent_output_bytes": persistent_bytes,
        "peak_materialized_output_bytes": peak_materialized_bytes,
        "transient_state_bytes_packaged_into_arrays": transient_bytes,
    }
    result = {
        "schema": SCHEMA,
        "status": "PASS",
        "run_id": case["run_id"],
        "source_forge_manifest_sha256": source_manifest_hash,
        "source_receipts": source_receipts,
        "source_forge_asset": {
            "manifest_path_at_run": case["source_manifest_path"],
            "manifest_schema": manifest["schema"],
            "geometry_scope": manifest["validation"]["scope"],
            "declared_geometry_status": manifest["validation"]["status"],
        },
        "geometry": case["geometry"],
        "material": {
            "rho_kg_m3": rho,
            "E_pa": float(case["material"]["E_pa"]),
            "nu": float(case["material"]["nu"]),
            "Gc_j_m2": float(case["material"]["Gc_j_m2"]),
            "calibration_status": case["material"]["calibration_status"],
            "model_mass_status": "DISCRETE_MODEL_MASS_FROM_APPROXIMATE_H3_VOLUMES_NOT_MEASURED",
            "total_mass_kg": float(np.sum(point_mass, dtype=np.float64)),
        },
        "loading": case["loading"],
        "solver": {**case["solver"], **native_status},
        "results": results,
        "array_contract": {
            "format": "NumPy npz, no object arrays",
            "units": "SI; explicit units are embedded in every dimensional key",
            "bond_active_packed_bitorder": "little; bit k is bond index k within each state row; 1=active after break evaluation at that state",
            "component_label": "0-based representative point ID of each active-graph component; includes isolated points",
            "damage": "native failed neighbor fraction after the current step; zero for points with zero initial degree",
            "first_break_step": "step at which strict CriticalStretch failure was evaluated at current positions; active flag is post-removal",
            "first_break_energy": "independent pair potential immediately before deletion at that same state; discrete event overshoot is retained",
            "required_keys": ["point_id", "x_ref_m", "position_m", "velocity_mps", "damage", "step", "time_s", "volume_m3", "point_mass_kg", "bond_i", "bond_j", "bond_rest_length_m", "bond_first_break_step", "bond_first_break_time_s", "bond_first_break_stretch", "bond_active_packed", "component_label", "event_bond_index", "event_step", "event_time_s", "event_stretch", "event_pair_energy_J"],
        },
        "arrays": {"path": "arrays.npz", "sha256": sha256_file(arrays_path), "bytes": arrays_path.stat().st_size, "format": "npz"},
        "outputs": outputs,
        "limitations": [
            "centroid-selected cells use h^3 point volumes and are not cut at the curved wall; volume and mass are approximate lattice model values",
            "BBMaterial fixes nu=0.25 and the illustrative prototype parameters are not calibrated to a ceramic specimen",
            "NoCorrection is used; free surfaces and discrete spherical-shell sampling can bias stiffness and fracture onset",
            "break_dissipation_cumulative_J sums independent Uij at discrete deletion; it is not a Gc times crack-area calibration",
            "the imposed affine velocity field is a free-release initial condition, not a drop or impact simulation",
            "no contact, external grip, precrack, fragment count target, convergence claim, or experimental validation is asserted",
        ],
    }
    _write_json(root / "result.json", result)
    # Full-state native binaries are packaged in arrays.npz; native CSV ledgers
    # and the exact source snapshot remain for independent readback.
    return result


def run_prepared(run_dir: str | os.PathLike[str], julia_executable: str | os.PathLike[str], *, timeout_s: float = MAX_RUNTIME_S) -> dict[str, Any]:
    """Run one previously prepared attempt and package it; never overwrite it."""
    root = Path(run_dir).resolve(strict=True)
    if any((root / name).exists() for name in ("native_status.json", "result.json", "run_receipt.json")):
        raise FileExistsError("run attempt already contains native output; refusing to overwrite")
    prepared_argv, _ = _read_json(root / "argv.json")
    source_records = _verify_source_snapshot(root)
    executable = Path(julia_executable).expanduser().resolve(strict=True)
    snapshot_runner = root / "source_snapshot" / "run_bowl_dynamics.jl"
    argv = [str(executable), f"--project={JULIA_PROJECT}", str(snapshot_runner), str(root / "case.toml"), str(root)]
    argv_sha256 = hashlib.sha256(json.dumps(argv, ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest()
    timeout = min(float(timeout_s), MAX_RUNTIME_S)
    receipt = {
        "schema": "pgfrag.native_bowl_execution_receipt.v1",
        "run_id": prepared_argv["run_id"],
        "argv": argv,
        "argv_sha256": argv_sha256,
        "prepared_argv_sha256": prepared_argv["argv_sha256"],
        "cwd": str(REPO_ROOT),
        "timeout_s": timeout,
        "source_receipts": source_records,
        "julia_project": str(JULIA_PROJECT),
        "julia_depot": str(REPO_ROOT / "third_party" / "pgfrag-runtime" / "depot"),
        "julia_threads": 1,
        "status": "STARTED",
    }
    _write_json(root / "run_receipt.json", receipt)
    started = time.monotonic()
    env = os.environ.copy()
    env["JULIA_DEPOT_PATH"] = receipt["julia_depot"]
    env["JULIA_NUM_THREADS"] = "1"
    env["JULIA_BLAS_NUM_THREADS"] = "1"
    env["JULIA_PKG_OFFLINE"] = "true"
    try:
        _verify_source_snapshot(root)
        completed = subprocess.run(argv, cwd=REPO_ROOT, env=env, capture_output=True, text=True,
                                   timeout=timeout, check=False)
        (root / "native_stdout.log").write_text(completed.stdout, encoding="utf-8")
        (root / "native_stderr.log").write_text(completed.stderr, encoding="utf-8")
        if completed.returncode != 0:
            raise RuntimeError(f"Julia native runner exited with {completed.returncode}")
        _verify_source_snapshot(root)
        result = package_run(root)
        result["solver"]["wall_runtime_s"] = time.monotonic() - started
        _write_json(root / "result.json", result)
        receipt["status"] = "PASS"
        receipt["wall_runtime_s"] = time.monotonic() - started
        receipt["native_stdout_sha256"] = sha256_file(root / "native_stdout.log")
        receipt["native_stderr_sha256"] = sha256_file(root / "native_stderr.log")
        _write_json(root / "run_receipt.json", receipt)
        return result
    except subprocess.TimeoutExpired as exc:
        stdout = exc.stdout.decode(errors="replace") if isinstance(exc.stdout, bytes) else (exc.stdout or "")
        stderr = exc.stderr.decode(errors="replace") if isinstance(exc.stderr, bytes) else (exc.stderr or "")
        (root / "native_stdout.log").write_text(stdout, encoding="utf-8")
        (root / "native_stderr.log").write_text(stderr, encoding="utf-8")
        failed = {"schema": SCHEMA, "status": "FAILED_TIMEOUT", "run_id": _read_json(root / "case.json")[0]["run_id"],
                  "elapsed_s": time.monotonic() - started, "timeout_s": timeout,
                  "argv": argv, "source_receipts": source_records}
        _write_json(root / "result.json", failed)
        receipt["status"] = "FAILED_TIMEOUT"
        receipt["wall_runtime_s"] = time.monotonic() - started
        _write_json(root / "run_receipt.json", receipt)
        return failed
    except Exception as exc:
        failed = {"schema": SCHEMA, "status": "FAILED", "run_id": _read_json(root / "case.json")[0]["run_id"],
                  "elapsed_s": time.monotonic() - started, "error": f"{type(exc).__name__}: {exc}",
                  "argv": argv, "source_receipts": source_records}
        _write_json(root / "result.json", failed)
        receipt["status"] = "FAILED"
        receipt["wall_runtime_s"] = time.monotonic() - started
        receipt["error"] = f"{type(exc).__name__}: {exc}"
        _write_json(root / "run_receipt.json", receipt)
        return failed


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Prepare or run a bounded native spherical-bowl BB case")
    parser.add_argument("asset", help="forge.asset.v1 spherical-bowl package or manifest.json")
    parser.add_argument("--runs-root", default=str(DATA_ROOT))
    parser.add_argument("--alpha", type=float, default=0.0, choices=(0.0, 100.0, 300.0), help="frozen imposed strain-rate scale in s^-1")
    parser.add_argument("--steps", type=int, default=400)
    parser.add_argument("--safety-factor", type=float, default=0.5)
    parser.add_argument("--spacing-m", type=float, default=1e-3, choices=FROZEN_SPACING_M)
    parser.add_argument("--julia", default=os.environ.get("PGFRAG_JULIA_EXE", ""))
    parser.add_argument("--execute", action="store_true", help="launch native Julia solve after preparing the attempt")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        run_dir = prepare_run(args.asset, runs_root=args.runs_root, alpha_s_1=args.alpha,
                              steps=args.steps, safety_factor=args.safety_factor, spacing_m=args.spacing_m)
        if not args.execute:
            print(json.dumps({"status": "SOURCE_READY_NOT_RUN", "run_dir": str(run_dir)}, indent=2))
            return 0
        if not args.julia:
            raise ValueError("--execute requires --julia or PGFRAG_JULIA_EXE")
        result = run_prepared(run_dir, args.julia)
        print(json.dumps({"status": result.get("status"), "run_dir": str(run_dir), "result_path": str(run_dir / "result.json")}, indent=2))
        return 0 if result.get("status") == "PASS" else 2
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"pgfrag bowl dynamics: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

