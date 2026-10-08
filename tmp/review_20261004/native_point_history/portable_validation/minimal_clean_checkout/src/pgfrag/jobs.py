from __future__ import annotations

import json
import math
from pathlib import Path, PurePosixPath, PureWindowsPath
from typing import Any

from . import BACKEND_ID
from .errors import UnsupportedCapability, ValidationError

_FORBIDDEN_COUNT_KEYS = {"target_fragments", "target_fragment_count", "target_pieces", "piece_count", "piece_range", "count_model"}


def _has_forbidden_count_control(value: Any) -> bool:
    if isinstance(value, dict):
        return any(str(key).lower() in _FORBIDDEN_COUNT_KEYS or _has_forbidden_count_control(child)
                   for key, child in value.items())
    if isinstance(value, list):
        return any(_has_forbidden_count_control(child) for child in value)
    return False


def _number(value: Any, name: str, *, minimum: float | None = None, exclusive_minimum: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValidationError(f"{name} must be a finite number")
    result = float(value)
    if minimum is not None and (result <= minimum if exclusive_minimum else result < minimum):
        op = ">" if exclusive_minimum else ">="
        raise ValidationError(f"{name} must be {op} {minimum}")
    return result


def _relative_path(value: Any, name: str) -> str:
    if not isinstance(value, str) or not value or "\x00" in value:
        raise ValidationError(f"{name} must be a non-empty relative path")
    posix = PurePosixPath(value.replace("\\", "/"))
    win = PureWindowsPath(value)
    if posix.is_absolute() or win.is_absolute() or win.drive or any(part in ("..", "") for part in posix.parts):
        raise ValidationError(f"{name} may not be absolute or contain parent traversal")
    return value


def safe_output_path(job_path: Path, relative: str) -> Path:
    name = _relative_path(relative, "output_dir")
    root = job_path.resolve(strict=True).parent
    destination = root.joinpath(*PurePosixPath(name.replace("\\", "/")).parts)
    resolved = destination.resolve(strict=False)
    if root != resolved and root not in resolved.parents:
        raise ValidationError("output_dir escapes the job directory")
    return resolved


def validate_job(job: Any) -> dict[str, Any]:
    if not isinstance(job, dict):
        raise ValidationError("job must be a JSON object")
    if _has_forbidden_count_control(job):
        raise ValidationError("jobs may not request or target fragment counts")
    expected = {"schema", "job_id", "backend_id", "initial_condition", "material", "loading", "solver", "resources", "output_dir"}
    if set(job) != expected:
        missing = sorted(expected - set(job))
        extra = sorted(set(job) - expected)
        raise ValidationError(f"job keys mismatch; missing={missing}, extra={extra}")
    if job["schema"] != "pgfrag.job.v1":
        raise ValidationError("schema must be 'pgfrag.job.v1'")
    job_id = job["job_id"]
    if not isinstance(job_id, str) or not job_id.strip() or len(job_id) > 128:
        raise ValidationError("job_id must be a non-empty string of at most 128 characters")
    if job["backend_id"] != BACKEND_ID:
        raise UnsupportedCapability(f"backend {job['backend_id']!r} is not installed in this PG-Frag build")

    initial = job["initial_condition"]
    if not isinstance(initial, dict):
        raise ValidationError("initial_condition must be an object")
    kind = initial.get("kind")
    if kind == "uniform_box_bar_fixture":
        if set(initial) != {"kind", "dimensions_mm", "dx_mm", "horizon_mm"}:
            raise ValidationError("uniform_box_bar_fixture requires only kind, dimensions_mm, dx_mm, horizon_mm")
        dims = initial["dimensions_mm"]
        if not isinstance(dims, list) or len(dims) != 3:
            raise ValidationError("dimensions_mm must contain three values")
        dims = [_number(v, f"dimensions_mm[{i}]", minimum=0.0, exclusive_minimum=True) for i, v in enumerate(dims)]
        dx = _number(initial["dx_mm"], "dx_mm", minimum=0.0, exclusive_minimum=True)
        horizon = _number(initial["horizon_mm"], "horizon_mm", minimum=0.0, exclusive_minimum=True)
        if any(length < dx for length in dims):
            raise ValidationError("each box dimension must be at least one point spacing")
        if horizon <= dx:
            raise ValidationError("horizon_mm must exceed dx_mm")
        estimated_counts = [max(1, math.floor(length / dx + 1e-10)) for length in dims]
        estimated_points = math.prod(estimated_counts)
        initial = {"kind": kind, "dimensions_mm": dims, "dx_mm": dx, "horizon_mm": horizon,
                   "estimated_lattice_counts": estimated_counts, "estimated_points": estimated_points}
    elif kind == "forge_asset_surface":
        if set(initial) != {"kind", "package_dir", "manifest_sha256"}:
            raise ValidationError("forge_asset_surface requires kind, package_dir, manifest_sha256")
        _relative_path(initial["package_dir"], "initial_condition.package_dir")
        digest = initial["manifest_sha256"]
        if not isinstance(digest, str) or len(digest) != 64 or any(c not in "0123456789abcdefABCDEF" for c in digest):
            raise ValidationError("manifest_sha256 must be a 64-character hexadecimal SHA-256")
    else:
        raise UnsupportedCapability(f"initial condition {kind!r} is not supported")

    material = job["material"]
    if not isinstance(material, dict) or set(material) != {"model", "rho_kg_m3", "E_pa", "nu", "Gc_j_m2", "calibration_status", "provenance"}:
        raise ValidationError("material requires model, rho_kg_m3, E_pa, nu, Gc_j_m2, calibration_status, provenance")
    model = material["model"]
    if model not in {"BB", "OSB"}:
        raise UnsupportedCapability(f"material model {model!r} is not supported")
    rho = _number(material["rho_kg_m3"], "rho_kg_m3", minimum=0, exclusive_minimum=True)
    young = _number(material["E_pa"], "E_pa", minimum=0, exclusive_minimum=True)
    nu = _number(material["nu"], "nu")
    gc = _number(material["Gc_j_m2"], "Gc_j_m2", minimum=0, exclusive_minimum=True)
    if not -1.0 < nu < 0.5:
        raise ValidationError("nu must be in the forge.asset.v1 contract domain (-1, 0.5)")
    if not 0 < nu < 0.5:
        raise UnsupportedCapability("this Peridynamics.jl material adapter supports only 0 < nu < 0.5; forge.asset.v1 still records its wider (-1, 0.5) input contract")
    if model == "BB" and abs(nu - 0.25) > 1e-12:
        raise UnsupportedCapability("Peridynamics.jl 3D BBMaterial requires nu=0.25; select native OSBMaterial for other supported nu values")
    if material["calibration_status"] != "UNVALIDATED_PROTOTYPE":
        raise UnsupportedCapability("this build only accepts UNVALIDATED_PROTOTYPE parameters and cannot claim calibrated material behavior")
    if not isinstance(material["provenance"], str) or not material["provenance"].strip():
        raise ValidationError("material.provenance must explain parameter origin")

    loading = job["loading"]
    if not isinstance(loading, dict) or set(loading) != {"kind", "speed_m_s"}:
        raise ValidationError("loading requires kind and speed_m_s")
    if loading["kind"] not in {"zero_load", "symmetric_axial_velocity"}:
        raise UnsupportedCapability(f"loading kind {loading['kind']!r} is not supported")
    speed = _number(loading["speed_m_s"], "speed_m_s", minimum=0.0)
    if loading["kind"] == "zero_load" and speed != 0:
        raise ValidationError("zero_load requires speed_m_s=0")
    if loading["kind"] == "symmetric_axial_velocity" and speed <= 0:
        raise ValidationError("symmetric_axial_velocity requires speed_m_s>0")

    solver = job["solver"]
    if not isinstance(solver, dict) or set(solver) != {"steps", "dt_s", "safety_factor"}:
        raise ValidationError("solver requires steps, dt_s, safety_factor")
    steps = solver["steps"]
    if isinstance(steps, bool) or not isinstance(steps, int) or steps <= 0:
        raise ValidationError("solver.steps must be a positive integer")
    dt = solver["dt_s"]
    if dt is not None:
        dt = _number(dt, "solver.dt_s", minimum=0, exclusive_minimum=True)
    safety = _number(solver["safety_factor"], "solver.safety_factor", minimum=0, exclusive_minimum=True)
    if safety >= 1:
        raise ValidationError("safety_factor must be < 1; native Peridynamics VelocityVerlet requires 0 < safety_factor < 1")

    resources = job["resources"]
    if not isinstance(resources, dict) or set(resources) != {"profile", "max_points", "max_steps", "threads", "timeout_seconds"}:
        raise ValidationError("resources require profile, max_points, max_steps, threads, timeout_seconds")
    if resources["profile"] != "local_smoke":
        raise UnsupportedCapability("only the bounded local_smoke resource profile is implemented")
    for key in ("max_points", "max_steps", "threads", "timeout_seconds"):
        value = resources[key]
        if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
            raise ValidationError(f"resources.{key} must be a positive integer")
    if resources["threads"] != 1:
        raise UnsupportedCapability("this fixture exporter supports one Julia thread/chunk; multi-chunk output is not implemented")
    if resources["max_points"] > 1000 or resources["max_steps"] > 500 or resources["timeout_seconds"] > 120:
        raise UnsupportedCapability("local_smoke hard ceiling is 1000 points, 500 steps, 1 Julia thread, 120 seconds")
    if kind == "uniform_box_bar_fixture":
        if initial["estimated_points"] > resources["max_points"]:
            raise ValidationError(f"estimated fixture point count {initial['estimated_points']} exceeds resources.max_points")
    if steps > resources["max_steps"]:
        raise ValidationError("solver.steps exceeds resources.max_steps")

    output_dir = _relative_path(job["output_dir"], "output_dir")
    return {**job, "initial_condition": initial, "material": {**material, "rho_kg_m3": rho, "E_pa": young, "nu": nu, "Gc_j_m2": gc},
            "loading": {**loading, "speed_m_s": speed},
            "solver": {**solver, "dt_s": dt, "safety_factor": safety}, "output_dir": output_dir}


def load_job(path: str | Path) -> tuple[dict[str, Any], Path]:
    job_path = Path(path).resolve(strict=True)
    try:
        raw = json.loads(job_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValidationError(f"cannot read job JSON: {exc}") from exc
    return validate_job(raw), job_path
