"""Independent read-only reviewer for the frozen NoFailureOSB coupon output.

The adapter now follows the frozen 21-field TOML sidecar and exact CSV headers
from the static source. It does not open summary.toml or trust its PASS fields.
The script is read-only with respect to solver inputs and writes only beneath
this review directory using the approved VF Python runtime (NumPy/Matplotlib).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import re
import tomllib
from dataclasses import dataclass
from pathlib import Path
from typing import Any


REQUIRED_FIELDS = ("U", "V", "vhalf", "accel", "bint", "bext", "psi", "damage", "activecount")
EXPECTED_UNITS = {
    "U": "m",
    "V": "m/s",
    "vhalf": "m/s",
    "accel": "m/s^2",
    "bint": "N/m^3",
    "bext": "N/m^3",
    "psi": "J/m^3",
    "damage": "1",
    "activecount": "count",
}
VECTOR_FIELDS = ("U", "V", "vhalf", "accel", "bint", "bext")
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[3]
FIXED_CASE_PATH = REPO_ROOT / "pgfrag/experiments/free_evolution_coupon_20261005/fixed_case.toml"
POINT_REFERENCE_HEADER = (
    "point_id", "x_ref_m", "y_ref_m", "z_ref_m", "volume_m3", "mass_kg",
    "complete_spherical_horizon", "minimum_face_clearance_minus_horizon_m",
    "left_grip", "right_grip", "reference_active_bonds",
)
ENERGY_HISTORY_HEADER = (
    "step", "time_s", "bc_impulse_interval_start_s", "bc_impulse_interval_end_s",
    "mass_kg", "momentum_x_kg_m_s", "momentum_y_kg_m_s", "momentum_z_kg_m_s",
    "native_net_force_x_N", "native_net_force_y_N", "native_net_force_z_N",
    "osb_strain_energy_J", "kinetic_energy_J", "left_grip_u_x_mean_m",
    "right_grip_u_x_mean_m", "left_bc_impulse_x_N_s", "left_bc_impulse_y_N_s",
    "left_bc_impulse_z_N_s", "left_bc_force_average_x_N", "left_bc_force_average_y_N",
    "left_bc_force_average_z_N", "right_bc_impulse_x_N_s", "right_bc_impulse_y_N_s",
    "right_bc_impulse_z_N_s", "right_bc_force_average_x_N", "right_bc_force_average_y_N",
    "right_bc_force_average_z_N", "left_bc_discrete_kinetic_work_J",
    "right_bc_discrete_kinetic_work_J", "total_bc_discrete_kinetic_work_step_J",
    "actuator_work_trapezoid_surrogate_step_J",
    "actuator_work_trapezoid_surrogate_cumulative_J",
    "max_free_dof_impulse_residual_N_s", "momentum_update_residual_x_N_s",
    "momentum_update_residual_y_N_s", "momentum_update_residual_z_N_s",
    "energy_balance_residual_diagnostic_J",
)


@dataclass
class CanonicalCase:
    """Schema-neutral arrays produced by the frozen-sidecar adapter.

    All state arrays are aligned to ``step_ids`` and ``point_ids``. Vector
    arrays have shape (steps, points, 3); scalar point arrays have shape
    (steps, points). ``grip_target_U_m`` has shape (steps, points, 3), with
    target values defined only for grip points. ``energy_native`` contains
    native CSV columns keyed by canonical names ``kinetic_energy_J`` and
    ``strain_energy_J`` and must align to ``time_s``.
    """

    source_root: Path
    sidecar_path: Path
    source_bindings: dict[str, dict[str, Any]]
    point_ids: Any
    step_ids: Any
    time_s: Any
    dt_s: Any
    xyz_ref_m: Any
    mass_kg: Any
    volume_m3: Any
    left_grip_mask: Any
    right_grip_mask: Any
    grip_mask: Any
    grip_target_U_m: Any
    fields: dict[str, Any]
    energy_native: dict[str, Any]
    energy_table: dict[str, Any]
    energy_time_s: Any
    metadata: dict[str, Any]


def fail_closed_input_path(root: Path, relative: str) -> Path:
    """Resolve a sidecar path without allowing it to escape the input root."""
    rel = Path(relative)
    if rel.is_absolute() or ".." in rel.parts or "\\" in relative or ":" in relative:
        raise ValueError(f"sidecar path is not canonical root-relative syntax: {relative!r}")
    unresolved = root / rel
    current = root
    for part in rel.parts:
        current = current / part
        if current.is_symlink():
            raise ValueError(f"symlink input is rejected: {relative!r}")
    candidate = unresolved.resolve(strict=True)
    resolved_root = root.resolve(strict=True)
    try:
        candidate.relative_to(resolved_root)
    except ValueError as exc:
        raise ValueError(f"sidecar path escapes input root: {relative!r}") from exc
    if not candidate.is_file():
        raise ValueError(f"sidecar input is not a regular file: {relative!r}")
    return candidate


def verify_bytes(path: Path, expected_bytes: int, expected_sha256: str) -> dict[str, Any]:
    digest = hashlib.sha256()
    actual_bytes = 0
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            actual_bytes += len(chunk)
            digest.update(chunk)
    actual_sha256 = digest.hexdigest()
    return {
        "path": str(path),
        "expected_bytes": int(expected_bytes),
        "actual_bytes": actual_bytes,
        "expected_sha256": expected_sha256.lower(),
        "actual_sha256": actual_sha256,
        "pass": actual_bytes == int(expected_bytes) and actual_sha256 == expected_sha256.lower(),
    }


def read_typed_binary(path: Path, dtype_token: str, shape: tuple[int, ...]) -> Any:
    """Read an exact-length typed array after the schema adapter validates dtype."""
    import numpy as np

    allowed = {"<f8": np.dtype("<f8"), "<i8": np.dtype("<i8"), "<i4": np.dtype("<i4"), "<u8": np.dtype("<u8"), "<u4": np.dtype("<u4")}
    if dtype_token not in allowed:
        raise ValueError(f"unsupported or non-little-endian dtype token: {dtype_token!r}")
    dtype = allowed[dtype_token]
    expected = math.prod(shape) * dtype.itemsize
    if path.stat().st_size != expected:
        raise ValueError(f"binary byte count/shape mismatch for {path.name}: expected {expected}, got {path.stat().st_size}")
    values = np.fromfile(path, dtype=dtype)
    return values.reshape(shape)


def load_frozen_case(input_root: Path) -> CanonicalCase:
    """Load the frozen 21-field TOML layout and its two CSV companions.

    This function is intentionally never called before the parent authorizes a
    review run. It does not open summary.toml or trust producer PASS claims.
    """
    import numpy as np

    root = input_root.resolve(strict=True)
    if not root.is_dir():
        raise ValueError("input root must name the coupon result directory")
    sidecar_path = fail_closed_input_path(root, "point_state_layout.toml")
    binary_path_ref = fail_closed_input_path(root, "point_state_history.f64le")
    point_csv_path = fail_closed_input_path(root, "point_reference.csv")
    energy_csv_path = fail_closed_input_path(root, "energy_reaction_history.csv")

    sidecar_raw = sidecar_path.read_bytes()
    sidecar = tomllib.loads(sidecar_raw.decode("utf-8"))
    if sidecar.get("schema") != "pgfrag.free_evolution_coupon.point_state_binary.v1":
        raise ValueError("unsupported point-state layout schema")
    expected_fields = [
        "ux_m", "uy_m", "uz_m",
        "vx_m_s", "vy_m_s", "vz_m_s",
        "vhalf_x_m_s", "vhalf_y_m_s", "vhalf_z_m_s",
        "accel_x_m_s2", "accel_y_m_s2", "accel_z_m_s2",
        "b_int_x_N_m3", "b_int_y_N_m3", "b_int_z_N_m3",
        "b_ext_x_N_m3", "b_ext_y_N_m3", "b_ext_z_N_m3",
        "osb_strain_energy_density_J_m3", "damage", "n_active_bonds_exact_float64",
    ]
    expected_field_units = [
        "m", "m", "m", "m/s", "m/s", "m/s", "m/s", "m/s", "m/s",
        "m/s^2", "m/s^2", "m/s^2", "N/m^3", "N/m^3", "N/m^3",
        "N/m^3", "N/m^3", "N/m^3", "J/m^3", "1", "count",
    ]
    if sidecar.get("format") != "IEEE-754 Float64, little-endian":
        raise ValueError("point-state format is not the expected little-endian Float64 encoding")
    expected_order = "Julia column-major tensor [field_index, point_id, state_index]; field index varies fastest, then ascending point_id, then ascending state_index"
    if sidecar.get("storage_order") != expected_order:
        raise ValueError("point-state storage order differs from the frozen Julia column-major layout")
    if sidecar.get("field_names") != expected_fields or sidecar.get("field_units") != expected_field_units:
        raise ValueError("point-state field names/order/units differ from the frozen 21-field contract")
    if sidecar.get("data_file") != "point_state_history.f64le":
        raise ValueError("sidecar data_file does not name point_state_history.f64le")
    if sidecar.get("point_id_mapping") != "point_id is 1..point_count; analytic reference coordinates and volumes are in point_reference.csv":
        raise ValueError("unsupported point-ID mapping")
    if sidecar.get("state_index_mapping") != "state_index = step + 1; step 0 is the reference state, step 1..n_steps are native VelocityVerlet output states":
        raise ValueError("unsupported state-index mapping")
    expected_time_levels = {
        "displacement": "at t_n",
        "full_velocity": "at t_n",
        "velocity_half": "for n>=1 centered at t_(n-1/2), after BC evaluated at t_n; n=0 is the initialized zero placeholder",
        "acceleration": "at t_n after native force update",
        "internal_force_density": "at t_n",
        "external_force_density": "at t_n",
        "osb_strain_energy_density": "at t_n from native export_field",
        "damage_and_active_bond_count": "at t_n; failure is disabled in this case",
    }
    if sidecar.get("state_time_levels") != expected_time_levels:
        raise ValueError("sidecar time-level convention differs from the pinned native VelocityVerlet hook")
    if sidecar.get("active_bond_count_encoding") != "Float64 with exact integer representability checked before writing":
        raise ValueError("unsupported active-bond count encoding")

    point_raw = point_csv_path.read_bytes()
    point_reader = csv.DictReader(io.StringIO(point_raw.decode("utf-8"), newline=""))
    if tuple(point_reader.fieldnames or ()) != POINT_REFERENCE_HEADER:
        raise ValueError("point_reference.csv header does not match the fixed source schema")
    point_rows = list(point_reader)
    energy_raw = energy_csv_path.read_bytes()
    energy_reader = csv.DictReader(io.StringIO(energy_raw.decode("utf-8"), newline=""))
    if tuple(energy_reader.fieldnames or ()) != ENERGY_HISTORY_HEADER:
        raise ValueError("energy_reaction_history.csv header does not match the fixed source schema")
    energy_rows = list(energy_reader)

    point_count = int(sidecar["point_count"])
    state_count = int(sidecar["state_count"])
    field_count = int(sidecar["field_count"])
    shape = tuple(int(v) for v in sidecar["array_shape"])
    if field_count != 21 or shape != (21, point_count, state_count):
        raise ValueError(f"sidecar shape/field count mismatch: shape={shape}, fields={field_count}")
    if int(sidecar["byte_count"]) != 8 * field_count * point_count * state_count:
        raise ValueError("sidecar byte_count does not equal Float64 tensor size")
    if len(point_rows) != point_count or len(energy_rows) != state_count:
        raise ValueError("CSV row counts do not match sidecar point/state counts")
    fixed_case_raw = FIXED_CASE_PATH.read_bytes()
    fixed_case = tomllib.loads(fixed_case_raw.decode("utf-8"))
    if point_count != int(fixed_case["geometry"]["expected_point_count"]):
        raise ValueError("point count does not match the frozen fixed case")
    if fixed_case["material"].get("failure_enabled") is not False or fixed_case["material"].get("surface_correction") != "NoCorrection":
        raise ValueError("result is not the frozen NoFailureOSB/NoCorrection diagnostic case")

    point_ids = np.asarray([int(row["point_id"]) for row in point_rows], dtype=np.int64)
    if not np.array_equal(point_ids, np.arange(1, point_count + 1, dtype=np.int64)):
        raise ValueError("point_reference point IDs must be complete canonical 1..N")
    xyz = np.column_stack([[float(row[f"{axis}_ref_m"]) for row in point_rows] for axis in "xyz"])
    volume = np.asarray([float(row["volume_m3"]) for row in point_rows], dtype=np.float64)
    mass = np.asarray([float(row["mass_kg"]) for row in point_rows], dtype=np.float64)
    ref_active = np.asarray([int(row["reference_active_bonds"]) for row in point_rows], dtype=np.int64)
    left_mask = np.asarray([row["left_grip"] == "true" for row in point_rows], dtype=bool)
    right_mask = np.asarray([row["right_grip"] == "true" for row in point_rows], dtype=bool)
    if any(row["left_grip"] not in {"true", "false"} or row["right_grip"] not in {"true", "false"} for row in point_rows):
        raise ValueError("grip CSV columns must be lowercase true/false")
    if np.any(left_mask & right_mask) or not left_mask.any() or not right_mask.any():
        raise ValueError("left and right grip rosters must be non-empty and disjoint")
    if not np.all(np.isfinite(xyz)) or not np.all(np.isfinite(volume)) or not np.all(np.isfinite(mass)):
        raise ValueError("geometry CSV contains non-finite coordinates/volumes/masses")

    geometry = fixed_case["geometry"]
    material = fixed_case["material"]
    loading = fixed_case["loading"]
    solver = fixed_case["solver"]
    limits = fixed_case["limits"]
    h = float(geometry["point_spacing_h_m"])
    rho = float(material["rho_kg_m3"])
    lengths = np.asarray([float(geometry[f"length_{axis}_m"]) for axis in "xyz"], dtype=np.float64)
    grid_shape = np.asarray(geometry["expected_grid_shape"], dtype=np.int64)
    if len(point_rows) != math.prod(grid_shape) or not np.allclose(volume, h**3, rtol=1e-12, atol=0.0):
        raise ValueError("point geometry does not match the analytic cell-centre volume/grid")
    if not np.allclose(mass, rho * volume, rtol=1e-12, atol=1e-18):
        raise ValueError("point mass is not fixed-case density times point volume")
    if [len(np.unique(xyz[:, axis])) for axis in range(3)] != grid_shape.tolist():
        raise ValueError("reference coordinate grid shape disagrees with fixed_case.toml")
    expected_centers = [np.linspace(-length / 2 + h / 2, length / 2 - h / 2, int(count)) for length, count in zip(lengths, grid_shape)]
    if any(not np.allclose(np.unique(xyz[:, axis]), expected_centers[axis], rtol=0.0, atol=64 * np.finfo(float).eps * max(lengths)) for axis in range(3)):
        raise ValueError("reference coordinate values disagree with the fixed analytic box")
    xmin, xmax = float(xyz[:, 0].min()), float(xyz[:, 0].max())
    tol = h * 1e-8
    expected_left = np.abs(xyz[:, 0] - xmin) <= tol
    expected_right = np.abs(xyz[:, 0] - xmax) <= tol
    if not np.array_equal(left_mask, expected_left) or not np.array_equal(right_mask, expected_right):
        raise ValueError("point_reference grip flags do not match the fixed outer x cell-centre layers")

    dt_s = float(sidecar["time_step_s"])
    if not math.isfinite(dt_s) or dt_s <= 0 or state_count < 2:
        raise ValueError("sidecar needs a finite positive time step and at least two states")
    n_steps = state_count - 1
    if n_steps > int(limits["max_steps"]) or n_steps > int(solver.get("max_steps", limits["max_steps"])):
        raise ValueError("state count exceeds the fixed step budget")
    total_duration = sum(float(loading[k]) for k in ("ramp_duration_s", "hold_duration_s", "unload_duration_s", "post_unload_dwell_duration_s"))
    expected_steps = math.ceil(total_duration / dt_s)
    if n_steps != expected_steps:
        raise ValueError(f"state count does not match ceil(total physical duration/dt): {n_steps} vs {expected_steps}")
    phase_durations = [float(loading[k]) for k in ("ramp_duration_s", "hold_duration_s", "unload_duration_s", "post_unload_dwell_duration_s")]
    positive_phase_intervals = [math.ceil(duration / dt_s) for duration in phase_durations if duration > 0]
    min_steps_per_phase = int(limits["min_steps_per_phase"])
    phase_resolution_cap = min(duration for duration in phase_durations if duration > 0) / min_steps_per_phase
    if dt_s > phase_resolution_cap * (1 + 1e-12) or not positive_phase_intervals or min(positive_phase_intervals) < min_steps_per_phase:
        raise ValueError("time step does not provide the fixed minimum resolution per positive-duration phase")
    times = np.arange(state_count, dtype=np.float64) * dt_s
    if times[-1] < total_duration or times[-1] - total_duration >= dt_s * (1 + 1e-12):
        raise ValueError("solver end time does not bracket the fixed physical schedule by less than one step")

    target_relative_strain = float(loading["target_engineering_strain"])
    target_displacement = target_relative_strain * (xmax - xmin)
    ramp = phase_durations[0]
    hold = phase_durations[1]
    unload = phase_durations[2]

    def target_u(t: float) -> float:
        t = max(t, 0.0)
        if t <= ramp:
            xi = t / ramp
            return target_displacement * (3 * xi**2 - 2 * xi**3)
        if t <= ramp + hold:
            return target_displacement
        if t <= ramp + hold + unload:
            xi = (t - ramp - hold) / unload
            return target_displacement * (1 - 3 * xi**2 + 2 * xi**3)
        return 0.0

    target_u_x = np.asarray([target_u(float(t)) for t in times], dtype=np.float64)
    grip_target = np.zeros((state_count, point_count, 3), dtype=np.float64)
    grip_target[:, right_mask, 0] = target_u_x[:, None]
    n_states, n_points, n_fields = state_count, point_count, field_count
    binary_rel = sidecar["data_file"]
    binary_path = fail_closed_input_path(root, binary_rel)
    binary_binding = verify_bytes(binary_path, int(sidecar["byte_count"]), str(sidecar["sha256"]))
    if not binary_binding["pass"]:
        raise ValueError("point-state binary does not match sidecar byte count and SHA-256")
    raw = read_typed_binary(binary_path, "<f8", (n_states, n_points, n_fields))
    if not np.isfinite(raw).all():
        raise ValueError("point-state binary contains NaN/Inf before field extraction")

    fields: dict[str, Any] = {}
    ranges = {
        "U": (0, 3), "V": (3, 6), "vhalf": (6, 9), "accel": (9, 12),
        "bint": (12, 15), "bext": (15, 18),
    }
    for name, (start, stop) in ranges.items():
        fields[name] = raw[:, :, start:stop]
    fields["psi"] = raw[:, :, 18]
    fields["damage"] = raw[:, :, 19]
    active_float = raw[:, :, 20]
    if not np.all(np.isfinite(active_float)) or np.any(active_float < 0) or not np.all(active_float == np.rint(active_float)):
        raise ValueError("Float64 active-bond counts are not finite exact nonnegative integers")
    if np.any(active_float > np.iinfo(np.int64).max):
        raise ValueError("active-bond count exceeds exact integer output range")
    fields["activecount"] = active_float.astype(np.int64)
    if not np.array_equal(fields["activecount"], np.broadcast_to(ref_active[None, :], (state_count, point_count))):
        raise ValueError("native active-bond count differs from point_reference or changes over time")

    if [int(row["step"]) for row in energy_rows] != list(range(state_count)):
        raise ValueError("energy history step IDs must enumerate 0..n_steps")
    energy_table: dict[str, Any] = {}
    for name in ENERGY_HISTORY_HEADER[1:]:
        energy_table[name] = np.asarray([float(row[name]) for row in energy_rows], dtype=np.float64)
    if not all(np.isfinite(values).all() for values in energy_table.values()):
        raise ValueError("energy/reaction CSV contains NaN or Inf")
    energy_native = {
        "kinetic_energy_J": energy_table["kinetic_energy_J"],
        "strain_energy_J": energy_table["osb_strain_energy_J"],
    }
    energy_time_s = energy_table["time_s"]

    def observed_binding(path: Path) -> dict[str, Any]:
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        return {"path": str(path), "bytes": path.stat().st_size, "sha256_observed": digest, "expected_digest_in_sidecar": False}

    source_bindings = {
        "point_state_binary": binary_binding,
        "point_state_layout_toml": {"path": str(sidecar_path), "bytes": len(sidecar_raw), "sha256_observed": hashlib.sha256(sidecar_raw).hexdigest(), "pass": True},
        "point_reference_csv": observed_binding(point_csv_path),
        "energy_reaction_history_csv": observed_binding(energy_csv_path),
        "fixed_case_toml": {"path": str(FIXED_CASE_PATH), "bytes": len(fixed_case_raw), "sha256_observed": hashlib.sha256(fixed_case_raw).hexdigest(), "pass": True},
    }
    output_files = []
    output_bytes = 0
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        size = path.stat().st_size
        output_bytes += size
        output_files.append({"relative_path": path.relative_to(root).as_posix(), "bytes": size})
    max_output_bytes = int(limits["max_output_bytes"])
    if output_bytes > max_output_bytes:
        raise ValueError("actual experiment result directory exceeds the fixed output budget")
    max_binary_bytes = 8 * 21 * int(geometry["expected_point_count"]) * (int(limits["max_steps"]) + 1)
    if max_binary_bytes + 5 * 1024**2 > max_output_bytes:
        raise ValueError("fixed-case maximum point tensor plus reserve exceeds its output budget")

    metadata = {
        "case_id": fixed_case["case_id"],
        "point_count": point_count,
        "step_count": n_steps,
        "binary_array_shape_field_point_state": [field_count, point_count, state_count],
        "binary_array_dtype": "<f8",
        "binary_storage_order": sidecar["storage_order"],
        "dt_s": dt_s,
        "phase_resolution_cap_s": phase_resolution_cap,
        "phase_resolution_cap_pass": dt_s <= phase_resolution_cap * (1 + 1e-12),
        "total_physical_duration_s": total_duration,
        "solver_end_time_s": float(times[-1]),
        "solver_end_time_overrun_s": float(times[-1] - total_duration),
        "target_engineering_strain": target_relative_strain,
        "target_grip_displacement_m": target_displacement,
        "left_grip_point_count": int(left_mask.sum()),
        "right_grip_point_count": int(right_mask.sum()),
        "interior_free_point_count": int((~(left_mask | right_mask)).sum()),
        "reference_total_mass_kg": float(mass.sum()),
        "reference_active_bonds_per_point_min": int(ref_active.min()),
        "reference_active_bonds_per_point_max": int(ref_active.max()),
        "reference_active_bond_incidence_sum": int(ref_active.sum()),
        "density_kg_m3": rho,
        "point_geometry_matches_fixed_case": True,
        "steps_match_fixed_physical_duration": n_steps == expected_steps,
        "steps_within_fixed_budget": n_steps <= int(limits["max_steps"]),
        "expected_binary_bytes_at_actual_steps": 8 * 21 * point_count * state_count,
        "actual_result_directory_bytes": output_bytes,
        "max_result_directory_bytes": max_output_bytes,
        "result_directory_files": output_files,
        "grip_displacement_axis_index": 0,
        "schedule": {"ramp_duration_s": ramp, "hold_duration_s": hold, "unload_duration_s": unload, "post_unload_dwell_duration_s": phase_durations[3]},
        "energy_balance_semantics": "diagnostic only: native staggered full-step energy minus accumulated half-step BC kinetic jumps; not qualified energy closure",
        "actuator_work_semantics": "surrogate only; excluded from the driver's energy residual",
    }
    return CanonicalCase(
        source_root=root,
        sidecar_path=sidecar_path,
        source_bindings=source_bindings,
        point_ids=point_ids,
        step_ids=np.arange(state_count, dtype=np.int64),
        time_s=times,
        dt_s=dt_s,
        xyz_ref_m=xyz,
        mass_kg=mass,
        volume_m3=volume,
        left_grip_mask=left_mask,
        right_grip_mask=right_mask,
        grip_mask=left_mask | right_mask,
        grip_target_U_m=grip_target,
        fields=fields,
        energy_native=energy_native,
        energy_table=energy_table,
        energy_time_s=energy_time_s,
        metadata=metadata,
    )


def _as_step_dt(dt_s: Any, n_steps: int) -> Any:
    import numpy as np

    dt = np.asarray(dt_s, dtype=np.float64)
    if dt.ndim == 0:
        return np.full(n_steps, float(dt), dtype=np.float64)
    if dt.shape != (n_steps,):
        raise ValueError(f"dt_s must be scalar or length step_count, got {dt.shape}")
    return dt


def _all_finite(case: CanonicalCase) -> dict[str, Any]:
    import numpy as np

    evidence: dict[str, Any] = {}
    for name, array in case.fields.items():
        values = np.asarray(array)
        finite = np.isfinite(values)
        evidence[name] = {"all_finite": bool(finite.all()), "nonfinite_count": int(values.size - finite.sum())}
    for name, array in (("xyz_ref_m", case.xyz_ref_m), ("mass_kg", case.mass_kg), ("volume_m3", case.volume_m3), ("time_s", case.time_s), ("grip_target_U_m", case.grip_target_U_m)):
        values = np.asarray(array)
        finite = np.isfinite(values)
        evidence[name] = {"all_finite": bool(finite.all()), "nonfinite_count": int(values.size - finite.sum())}
    for name, array in case.energy_table.items():
        values = np.asarray(array)
        finite = np.isfinite(values)
        evidence[f"energy_csv.{name}"] = {"all_finite": bool(finite.all()), "nonfinite_count": int(values.size - finite.sum())}
    return evidence


def _close(actual: Any, expected: Any, *, rtol: float, atol: float) -> tuple[bool, float, float]:
    import numpy as np

    a = np.asarray(actual, dtype=np.float64)
    b = np.asarray(expected, dtype=np.float64)
    if a.shape != b.shape:
        return False, math.inf, math.inf
    delta = np.abs(a - b)
    scale = np.maximum(np.abs(a), np.abs(b))
    relative = np.divide(delta, scale, out=np.zeros_like(delta), where=scale > 0)
    return bool(np.allclose(a, b, rtol=rtol, atol=atol)), float(delta.max(initial=0.0)), float(relative.max(initial=0.0))


def audit_case(case: CanonicalCase) -> dict[str, Any]:
    """Recompute state/energy/momentum checks from typed arrays, not summaries."""
    import numpy as np

    fields = {name: np.asarray(case.fields[name]) for name in REQUIRED_FIELDS}
    u, v, vhalf = fields["U"], fields["V"], fields["vhalf"]
    accel, bint, bext = fields["accel"], fields["bint"], fields["bext"]
    psi, damage, activecount = fields["psi"], fields["damage"], fields["activecount"]
    times = np.asarray(case.time_s, dtype=np.float64)
    point_ids = np.asarray(case.point_ids)
    step_ids = np.asarray(case.step_ids)
    mass = np.asarray(case.mass_kg, dtype=np.float64)
    volume = np.asarray(case.volume_m3, dtype=np.float64)
    grip = np.asarray(case.grip_mask, dtype=bool)
    n_states, n_points = len(times), len(point_ids)
    n_steps = n_states - 1

    expected_shapes = {
        **{name: (n_states, n_points, 3) for name in VECTOR_FIELDS},
        **{name: (n_states, n_points) for name in ("psi", "damage", "activecount")},
    }
    shape_checks = {name: {"actual": list(fields[name].shape), "expected": list(expected_shapes[name]), "pass": fields[name].shape == expected_shapes[name]} for name in REQUIRED_FIELDS}
    if any(not item["pass"] for item in shape_checks.values()):
        raise ValueError(f"state array shape mismatch: {shape_checks}")
    if mass.shape != (n_points,) or volume.shape != (n_points,) or grip.shape != (n_points,):
        raise ValueError("geometry mass/volume/grip roster does not align with point IDs")
    if len(set(point_ids.tolist())) != n_points or len(set(step_ids.tolist())) != n_states:
        raise ValueError("point or step IDs are not unique")
    if np.any(np.diff(times) <= 0) or np.any(mass <= 0) or np.any(volume <= 0):
        raise ValueError("time must increase and point masses/volumes must be positive")
    if np.any(activecount < 0) or not np.all(activecount == np.rint(activecount)):
        raise ValueError("activecount contains negative or non-integer values")

    # Core independent aggregates, using dictionary units and actual arrays.
    kinetic_j = 0.5 * np.sum(mass[None, :] * np.sum(v * v, axis=2), axis=1)
    strain_j = np.sum(psi * volume[None, :], axis=1)
    momentum = np.sum(mass[None, :, None] * v, axis=1)
    energy_checks: dict[str, Any] = {}
    for key, recomputed in (("kinetic_energy_J", kinetic_j), ("strain_energy_J", strain_j)):
        native = np.asarray(case.energy_native[key], dtype=np.float64)
        matched, max_abs, max_rel = _close(native, recomputed, rtol=1e-10, atol=1e-12)
        energy_checks[key] = {"native_rows": int(native.size), "recomputed_rows": int(recomputed.size), "pass": matched, "max_abs_delta": max_abs, "max_relative_delta": max_rel}

    dt = _as_step_dt(case.dt_s, n_steps)
    # For row k>=1, native storage.velocity_half is the actual half velocity
    # after BC at t_k and is centered at t_(k-1/2); row 0 is a zero placeholder.
    vhalf_predicted = v[:-1] + 0.5 * accel[:-1] * dt[:, None, None]
    bc_impulse_by_point = mass[None, :, None] * (vhalf[1:] - vhalf_predicted)
    left_mask = np.asarray(case.left_grip_mask, dtype=bool)
    right_mask = np.asarray(case.right_grip_mask, dtype=bool)
    free_mask = ~grip
    left_impulse = np.zeros((n_steps, 3), dtype=np.float64)
    right_impulse = np.zeros((n_steps, 3), dtype=np.float64)
    left_work = np.zeros(n_steps, dtype=np.float64)
    right_work = np.zeros(n_steps, dtype=np.float64)
    max_free_impulse = np.zeros(n_steps, dtype=np.float64)
    free_impulse_pass = True
    max_free_impulse_ratio = 0.0
    for k in range(n_steps):
        left_impulse[k] = np.sum(bc_impulse_by_point[k, left_mask, :], axis=0)
        right_impulse[k] = np.sum(bc_impulse_by_point[k, right_mask, :], axis=0)
        before = vhalf_predicted[k]
        after = vhalf[k + 1]
        left_work[k] = np.sum(0.5 * mass[left_mask, None] * (after[left_mask] ** 2 - before[left_mask] ** 2))
        right_work[k] = np.sum(0.5 * mass[right_mask, None] * (after[right_mask] ** 2 - before[right_mask] ** 2))
        max_free_impulse[k] = float(np.max(np.abs(bc_impulse_by_point[k, free_mask, :]), initial=0.0))
        local_roundoff_scale = mass[free_mask, None] * (
            np.abs(vhalf[k + 1, free_mask])
            + np.abs(v[k, free_mask])
            + np.abs(accel[k, free_mask] * dt[k] / 2)
        )
        local_tolerance = 256 * np.finfo(np.float64).eps * local_roundoff_scale + 1e-30
        local_ratio = np.divide(
            np.abs(bc_impulse_by_point[k, free_mask]),
            local_tolerance,
            out=np.zeros_like(local_tolerance),
            where=local_tolerance > 0,
        )
        step_ratio = float(np.max(local_ratio, initial=0.0))
        max_free_impulse_ratio = max(max_free_impulse_ratio, step_ratio)
        free_impulse_pass = free_impulse_pass and step_ratio <= 1.0
    bc_impulse_by_step = left_impulse + right_impulse
    free_impulse_by_step = np.asarray([
        np.sum(bc_impulse_by_point[k, free_mask, :], axis=0) for k in range(n_steps)
    ])
    kinetic_jump_j = np.diff(kinetic_j)
    momentum_jump = np.diff(momentum, axis=0)
    force_sum = np.sum((bint + bext) * volume[None, :, None], axis=1)
    force_sum_prev = force_sum[:-1]
    force_sum_next = force_sum[1:]
    trapezoid_force_impulse = 0.5 * dt[:, None] * (force_sum_prev + force_sum_next)
    momentum_rhs = trapezoid_force_impulse + bc_impulse_by_step
    momentum_residual = momentum_jump - momentum_rhs

    acceleration_from_forces = (bint + bext) * volume[None, :, None] / mass[None, :, None]
    acceleration_ok, accel_abs, accel_rel = _close(accel, acceleration_from_forces, rtol=1e-9, atol=1e-10)
    free_accel_ok, free_accel_abs, free_accel_rel = _close(accel[:, free_mask], acceleration_from_forces[:, free_mask], rtol=1e-9, atol=1e-10)
    free_vhalf_ok, free_vhalf_abs, free_vhalf_rel = _close(vhalf[1:, free_mask], vhalf_predicted[:, free_mask], rtol=1e-9, atol=1e-12)
    full_velocity_predicted = vhalf[1:] + accel[1:] * dt[:, None, None] / 2
    velocity_update_ok, velocity_update_abs, velocity_update_rel = _close(v[1:], full_velocity_predicted, rtol=1e-9, atol=1e-12)
    displacement_predicted = u[:-1] + vhalf[1:] * dt[:, None, None]
    free_u_ok, free_u_abs, free_u_rel = _close(u[1:, free_mask], displacement_predicted[:, free_mask], rtol=1e-9, atol=1e-12)

    target = np.asarray(case.grip_target_U_m, dtype=np.float64)
    if target.shape != u.shape:
        raise ValueError(f"grip target shape {target.shape} must equal U shape {u.shape}")
    target_error = u[:, grip, :] - target[:, grip, :]
    grip_path_tolerance = max(256 * np.spacing(float(case.metadata["target_grip_displacement_m"])), 1e-14)
    grip_u_abs = float(np.max(np.abs(target_error), initial=0.0))
    grip_target_scale = max(float(np.max(np.abs(target[:, grip, :]), initial=0.0)), 1e-300)
    grip_u_rel = grip_u_abs / grip_target_scale
    grip_u_ok = grip_u_abs <= grip_path_tolerance
    right_grip_point_spread = np.max(np.abs(u[:, right_mask, 0] - np.mean(u[:, right_mask, 0], axis=1)[:, None]), axis=1)
    left_grip_zero_ok = bool(np.max(np.abs(u[:, left_mask, :]), initial=0.0) <= grip_path_tolerance)
    right_grip_yz_zero_ok = bool(np.max(np.abs(u[:, right_mask, 1:3]), initial=0.0) <= grip_path_tolerance)
    grip_spread_ok = bool(np.max(right_grip_point_spread, initial=0.0) <= grip_path_tolerance)
    no_failure_damage_ok = bool(np.allclose(damage, 0.0, rtol=0.0, atol=1e-14))
    active_constant = bool(np.all(activecount == activecount[0:1, :]))
    vhalf_initial_zero = bool(np.all(vhalf[0] == 0.0))

    energy_time = np.asarray(case.energy_time_s, dtype=np.float64)
    energy_time_ok, energy_time_abs, _ = _close(energy_time, times, rtol=1e-12, atol=1e-15)
    binding_ok = case.source_bindings["point_state_binary"].get("pass") is True
    all_finite = _all_finite(case)
    finite_ok = all(item["all_finite"] for item in all_finite.values())
    shape_ok = all(item["pass"] for item in shape_checks.values())

    state_mass = np.full(n_states, float(mass.sum()), dtype=np.float64)
    cumulative_bc_work = np.cumsum(np.concatenate(([0.0], left_work + right_work)))
    energy_balance_diag = strain_j + kinetic_j - (strain_j[0] + kinetic_j[0]) - cumulative_bc_work
    total_energy = strain_j + kinetic_j
    maximum_total_energy = float(np.max(total_energy))
    final_total_energy = float(total_energy[-1])
    final_energy_residual = float(energy_balance_diag[-1])
    final_energy_residual_fraction = abs(final_energy_residual) / max(maximum_total_energy, 1e-300)
    final_energy_residual_fraction_of_final = abs(final_energy_residual) / max(final_total_energy, 1e-300)
    schedule = case.metadata["schedule"]
    dwell_start_s = schedule["ramp_duration_s"] + schedule["hold_duration_s"] + schedule["unload_duration_s"]
    dwell_start_index = int(np.argmin(np.abs(times - dwell_start_s)))
    dwell_kinetic_change = float(kinetic_j[-1] - kinetic_j[dwell_start_index])
    dwell_strain_change = float(strain_j[-1] - strain_j[dwell_start_index])
    right_u_x = np.mean(u[:, right_mask, 0], axis=1)
    left_u_x = np.mean(u[:, left_mask, 0], axis=1)
    left_reaction_force = np.vstack((np.zeros((1, 3)), left_impulse / dt[:, None]))
    right_reaction_force = np.vstack((np.zeros((1, 3)), right_impulse / dt[:, None]))
    bc_work_step = np.concatenate(([0.0], left_work + right_work))
    left_work_state = np.concatenate(([0.0], left_work))
    right_work_state = np.concatenate(([0.0], right_work))
    actuator_step = np.zeros(n_states, dtype=np.float64)
    for k in range(1, n_states):
        actuator_step[k] = 0.5 * (right_reaction_force[k - 1, 0] + right_reaction_force[k, 0]) * (right_u_x[k] - right_u_x[k - 1])
    actuator_cumulative = np.cumsum(actuator_step)
    max_free_by_state = np.concatenate(([0.0], max_free_impulse))
    momentum_residual_state = np.vstack((np.zeros((1, 3)), momentum_residual))
    force_state = force_sum
    momentum_state = momentum

    csv_expectations: dict[str, Any] = {
        "time_s": times,
        "bc_impulse_interval_start_s": np.concatenate(([0.0], np.arange(n_steps, dtype=np.float64) * float(case.dt_s))),
        "bc_impulse_interval_end_s": times,
        "mass_kg": state_mass,
        "momentum_x_kg_m_s": momentum_state[:, 0],
        "momentum_y_kg_m_s": momentum_state[:, 1],
        "momentum_z_kg_m_s": momentum_state[:, 2],
        "native_net_force_x_N": force_state[:, 0],
        "native_net_force_y_N": force_state[:, 1],
        "native_net_force_z_N": force_state[:, 2],
        "osb_strain_energy_J": strain_j,
        "kinetic_energy_J": kinetic_j,
        "left_grip_u_x_mean_m": left_u_x,
        "right_grip_u_x_mean_m": right_u_x,
        "left_bc_impulse_x_N_s": np.vstack((np.zeros((1, 3)), left_impulse))[:, 0],
        "left_bc_impulse_y_N_s": np.vstack((np.zeros((1, 3)), left_impulse))[:, 1],
        "left_bc_impulse_z_N_s": np.vstack((np.zeros((1, 3)), left_impulse))[:, 2],
        "left_bc_force_average_x_N": left_reaction_force[:, 0],
        "left_bc_force_average_y_N": left_reaction_force[:, 1],
        "left_bc_force_average_z_N": left_reaction_force[:, 2],
        "right_bc_impulse_x_N_s": np.vstack((np.zeros((1, 3)), right_impulse))[:, 0],
        "right_bc_impulse_y_N_s": np.vstack((np.zeros((1, 3)), right_impulse))[:, 1],
        "right_bc_impulse_z_N_s": np.vstack((np.zeros((1, 3)), right_impulse))[:, 2],
        "right_bc_force_average_x_N": right_reaction_force[:, 0],
        "right_bc_force_average_y_N": right_reaction_force[:, 1],
        "right_bc_force_average_z_N": right_reaction_force[:, 2],
        "left_bc_discrete_kinetic_work_J": left_work_state,
        "right_bc_discrete_kinetic_work_J": right_work_state,
        "total_bc_discrete_kinetic_work_step_J": bc_work_step,
        "actuator_work_trapezoid_surrogate_step_J": actuator_step,
        "actuator_work_trapezoid_surrogate_cumulative_J": actuator_cumulative,
        "max_free_dof_impulse_residual_N_s": max_free_by_state,
        "momentum_update_residual_x_N_s": momentum_residual_state[:, 0],
        "momentum_update_residual_y_N_s": momentum_residual_state[:, 1],
        "momentum_update_residual_z_N_s": momentum_residual_state[:, 2],
        "energy_balance_residual_diagnostic_J": energy_balance_diag,
    }
    csv_checks: dict[str, Any] = {}
    for column, expected in csv_expectations.items():
        actual = np.asarray(case.energy_table[column], dtype=np.float64)
        matched, max_abs, max_rel = _close(actual, expected, rtol=1e-9, atol=1e-12)
        csv_checks[column] = {"pass": matched, "max_abs_delta": max_abs, "max_relative_delta": max_rel}
    csv_checks_pass = all(item["pass"] for item in csv_checks.values())
    momentum_ok, momentum_abs, momentum_rel = _close(momentum_jump, momentum_rhs, rtol=1e-9, atol=1e-18)
    velocity_update_ok = velocity_update_ok and vhalf_initial_zero

    checks = {
        "binary_file_sha256_and_byte_counts_match_toml": binding_ok,
        "binary_shapes_match_sidecar": shape_ok,
        "all_float_state_and_geometry_values_finite": finite_ok,
        "point_state_bond_activity_constant_no_failure": active_constant,
        "damage_remains_zero_no_failure": no_failure_damage_ok,
        "state_mass_from_geometry_is_constant": bool(np.all(state_mass == state_mass[0])),
        "native_acceleration_matches_force_density_times_volume_over_mass": acceleration_ok,
        "native_velocity_verlet_full_velocity_level_matches": velocity_update_ok,
        "native_kinetic_energy_matches_sum_0p5_m_v2": energy_checks["kinetic_energy_J"]["pass"],
        "native_strain_energy_matches_sum_psi_volume": energy_checks["strain_energy_J"]["pass"],
        "native_energy_time_matches_state_time": energy_time_ok,
        "fixed_secant_grip_U_matches_target": grip_u_ok,
        "left_grip_U_zero_and_right_grip_yz_zero": left_grip_zero_ok and right_grip_yz_zero_ok,
        "right_grip_points_share_target_path": grip_spread_ok,
        "free_dof_acceleration_matches_force_density_times_volume_over_mass": free_accel_ok,
        "free_dof_half_step_update_has_no_boundary_impulse": free_vhalf_ok,
        "free_dof_displacement_matches_vhalf_integral": free_u_ok,
        "free_dof_impulse_within_native_roundoff_tolerance": free_impulse_pass,
        "full_step_momentum_identity_force_and_BC_impulse": momentum_ok,
        "energy_reaction_history_csv_matches_independent_recomputation": csv_checks_pass,
    }
    return {
        "schema": "free_coupon_independent_review.report.v1",
        "status": "PASS" if all(checks.values()) else "FAIL",
        "scope": "native no-failure coupon state, energy, impulse, and momentum arithmetic only",
        "source_bindings": case.source_bindings,
        "metadata": case.metadata,
        "units_expected": EXPECTED_UNITS,
        "point_count": n_points,
        "step_count": n_steps,
        "step_ids_first_last": [int(step_ids[0]), int(step_ids[-1])],
        "time_s_first_last": [float(times[0]), float(times[-1])],
        "shape_checks": shape_checks,
        "finite_checks": all_finite,
        "energy_checks": energy_checks,
        "energy_qualification": {
            "status": "NOT_QUALIFIED",
            "final_diagnostic_residual_J": final_energy_residual,
            "final_total_energy_J": final_total_energy,
            "maximum_total_energy_J": maximum_total_energy,
            "absolute_residual_fraction_of_final_total_energy": final_energy_residual_fraction_of_final,
            "absolute_residual_fraction_of_max_total_energy": final_energy_residual_fraction,
            "diagnostic_only": True,
            "note": "The short post-unload dwell has not demonstrated stationarity; this reversible no-failure run does not qualify energy closure or physical fracture behavior.",
        },
        "post_unload_dwell": {
            "scheduled_start_time_s": float(dwell_start_s),
            "start_time_s": float(times[dwell_start_index]),
            "end_time_s": float(times[-1]),
            "kinetic_energy_change_J": dwell_kinetic_change,
            "strain_energy_change_J": dwell_strain_change,
            "final_max_free_speed_m_s": float(np.max(np.linalg.norm(v[-1, ~grip, :], axis=1), initial=0.0)),
        },
        "energy_reaction_csv_checks": csv_checks,
        "derived": {
            "state_mass_kg": state_mass.tolist(),
            "kinetic_energy_J": kinetic_j.tolist(),
            "strain_energy_J": strain_j.tolist(),
            "momentum_kg_mps": momentum.tolist(),
            "kinetic_jump_J": kinetic_jump_j.tolist(),
            "bc_impulse_per_step_Ns": bc_impulse_by_step.tolist(),
            "bc_impulse_cumulative_Ns": np.cumsum(bc_impulse_by_step, axis=0).tolist(),
            "left_bc_impulse_per_step_Ns": left_impulse.tolist(),
            "right_bc_impulse_per_step_Ns": right_impulse.tolist(),
            "left_bc_kinetic_work_jump_J": left_work.tolist(),
            "right_bc_kinetic_work_jump_J": right_work.tolist(),
            "full_step_kinetic_jump_J": kinetic_jump_j.tolist(),
            "free_dof_impulse_per_step_Ns": free_impulse_by_step.tolist(),
            "max_free_dof_impulse_per_step_Ns": max_free_impulse.tolist(),
            "momentum_jump_kg_mps": momentum_jump.tolist(),
            "momentum_identity_rhs_kg_mps": momentum_rhs.tolist(),
            "momentum_identity_residual_kg_mps": momentum_residual.tolist(),
            "native_net_force_N": force_sum.tolist(),
            "energy_balance_diagnostic_J": energy_balance_diag.tolist(),
            "actuator_work_trapezoid_surrogate_J": actuator_cumulative.tolist(),
        },
        "checks": checks,
        "numerical_residuals": {
            "grip_target_U_max_abs_m": grip_u_abs,
            "grip_target_U_max_relative": grip_u_rel,
            "grip_target_tolerance_m": grip_path_tolerance,
            "right_grip_point_spread_max_m": float(np.max(right_grip_point_spread, initial=0.0)),
            "free_acceleration_max_abs_mps2": free_accel_abs,
            "free_acceleration_max_relative": free_accel_rel,
            "free_half_step_max_abs_mps": free_vhalf_abs,
            "free_half_step_max_relative": free_vhalf_rel,
            "free_displacement_max_abs_m": free_u_abs,
            "free_displacement_max_relative": free_u_rel,
            "native_energy_time_max_abs_s": energy_time_abs,
            "native_velocity_update_max_abs_mps": velocity_update_abs,
            "native_velocity_update_max_relative": velocity_update_rel,
            "native_acceleration_update_max_abs_mps2": accel_abs,
            "native_acceleration_update_max_relative": accel_rel,
            "momentum_identity_max_abs_kg_mps": momentum_abs,
            "momentum_identity_max_relative": momentum_rel,
            "free_impulse_max_abs_N_s": float(np.max(max_free_impulse, initial=0.0)),
        "free_impulse_max_roundoff_tolerance_ratio": max_free_impulse_ratio,
        },
        "no_failure_state": {"damage_max": float(np.max(damage, initial=0.0)), "activecount_constant": active_constant, "reference_activecount_match": True},
        "interpretation_limits": [
            "NoFailureOSB is a reversible implementation diagnostic; it does not establish physical fracture qualification.",
            "Point clouds and curves show only values present in the verified native arrays; no crack surface or fracture geometry is synthesized.",
            "A passing arithmetic audit does not establish calibrated material physics or general solver validity.",
        ],
    }


def make_figures(case: CanonicalCase, report: dict[str, Any], output_dir: Path) -> None:
    """Plot actual points and recomputed curves only; no surface/crack rendering."""
    import matplotlib.pyplot as plt
    import numpy as np

    fields = {name: np.asarray(case.fields[name]) for name in REQUIRED_FIELDS}
    u = fields["U"]
    xyz = np.asarray(case.xyz_ref_m, dtype=np.float64)
    left_mask = np.asarray(case.left_grip_mask, dtype=bool)
    right_mask = np.asarray(case.right_grip_mask, dtype=bool)
    times = np.asarray(case.time_s, dtype=np.float64)
    kinetic = np.asarray(report["derived"]["kinetic_energy_J"])
    strain = np.asarray(report["derived"]["strain_energy_J"])
    target = np.asarray(case.grip_target_U_m)
    grip_ids = np.flatnonzero(right_mask)
    peak_index = int(np.argmax(np.mean(target[:, right_mask, 0], axis=1)))
    schedule = case.metadata["schedule"]
    unload_end_s = schedule["ramp_duration_s"] + schedule["hold_duration_s"] + schedule["unload_duration_s"]
    unload_end_index = int(np.argmin(np.abs(times - unload_end_s)))
    state_indices = [0, peak_index, unload_end_index, len(times) - 1]
    actual_coords = [xyz + u[k] for k in state_indices]
    labels = ["initial", "peak load", "nearest state after unload", "post-unload dwell end"]
    all_xyz = np.concatenate(actual_coords, axis=0)
    lows, highs = all_xyz.min(axis=0), all_xyz.max(axis=0)
    center = (lows + highs) / 2.0
    half = max(float(np.max(highs - lows)) / 2.0, 1e-12)

    from matplotlib.colors import Normalize

    fig = plt.figure(figsize=(14, 10))
    fig.subplots_adjust(left=0.04, right=0.89, top=0.90, bottom=0.06, wspace=0.12, hspace=0.22)
    max_displacement_mm = max(float(np.max(np.linalg.norm(u[state], axis=1) * 1000.0, initial=0.0)) for state in state_indices)
    color_norm = Normalize(vmin=0.0, vmax=max(max_displacement_mm, 1e-15))
    color_artist = None
    for panel, (coords, label, index) in enumerate(zip(actual_coords, labels, state_indices), start=1):
        ax = fig.add_subplot(2, 2, panel, projection="3d")
        displacement_mm = np.linalg.norm(u[index], axis=1) * 1000.0
        color_artist = ax.scatter(coords[:, 0] * 1000, coords[:, 1] * 1000, coords[:, 2] * 1000, c=displacement_mm, s=5, cmap="viridis", norm=color_norm)
        ax.set_title(f"{label}, t={times[index]:.6g} s")
        ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)"); ax.set_zlabel("z (mm)")
        ax.set_xlim((center[0] - half) * 1000, (center[0] + half) * 1000)
        ax.set_ylim((center[1] - half) * 1000, (center[1] + half) * 1000)
        ax.set_zlim((center[2] - half) * 1000, (center[2] + half) * 1000)
        ax.set_box_aspect((1, 1, 1))
    fig.suptitle("Native coupon points only; equal world-unit axes; no fracture surface")
    colorbar_axes = fig.add_axes([0.92, 0.20, 0.018, 0.60])
    fig.colorbar(color_artist, cax=colorbar_axes, label="|U| (mm; color only)")
    fig.savefig(output_dir / "coupon_points_3d.png", dpi=160)
    plt.close(fig)

    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    ax = axs[0, 0]
    if grip_ids.size:
        axis_index = int(case.metadata.get("grip_displacement_axis_index", 0))
        actual_grip = np.mean(u[:, grip_ids, axis_index], axis=1)
        target_grip = np.mean(target[:, grip_ids, axis_index], axis=1)
        ax.plot(times, actual_grip * 1000, label="actual grip U")
        ax.plot(times, target_grip * 1000, "--", label="secant target U")
    ax.set_title("Grip displacement"); ax.set_xlabel("time (s)"); ax.set_ylabel("U (mm)"); ax.legend()

    ax = axs[0, 1]
    ax.plot(times, kinetic, label="recomputed KE")
    ax.plot(times, strain, label="recomputed ΣψV")
    for key, label in (("kinetic_energy_J", "native KE"), ("strain_energy_J", "native strain energy")):
        ax.plot(case.energy_time_s, case.energy_native[key], linestyle="--", label=label)
    ax.set_title("Native and recomputed energy"); ax.set_xlabel("time (s)"); ax.set_ylabel("energy (J)"); ax.legend()

    ax = axs[1, 0]
    impulse = np.asarray(report["derived"]["bc_impulse_per_step_Ns"])
    impulse_time = times[1:]
    for axis, label in enumerate("xyz"):
        ax.plot(impulse_time, impulse[:, axis], label=f"J{label}")
    ax.set_title("Actual BC impulse per step"); ax.set_xlabel("time (s)"); ax.set_ylabel("impulse (N·s)"); ax.legend()

    ax = axs[1, 1]
    ax.plot(times[1:], np.diff(kinetic), label="full-step ΔKE")
    bc_work = np.asarray(report["derived"]["left_bc_kinetic_work_jump_J"]) + np.asarray(report["derived"]["right_bc_kinetic_work_jump_J"])
    ax.plot(times[1:], bc_work, label="half-step BC kinetic jump")
    ax.set_title("Kinetic energy jumps (staggered levels)"); ax.set_xlabel("time (s)"); ax.set_ylabel("energy jump (J)"); ax.legend()
    fig.tight_layout()
    fig.savefig(output_dir / "coupon_diagnostics.png", dpi=160)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, required=True, help="explicit frozen result root")
    parser.add_argument("--output-dir", type=Path, default=SCRIPT_DIR / "audit_output", help="must be a child of this review directory")
    args = parser.parse_args()
    output_dir = args.output_dir.resolve(strict=False)
    review_dir = SCRIPT_DIR.resolve(strict=True)
    input_root = args.input_root.resolve(strict=True)
    if output_dir == review_dir or not output_dir.is_relative_to(review_dir):
        raise ValueError(f"output directory must be a child of {review_dir}")
    if output_dir == input_root or input_root in output_dir.parents or output_dir in input_root.parents:
        raise ValueError("output directory must not overlap the frozen solver result directory")
    case = load_frozen_case(input_root)
    report = audit_case(case)
    output_dir.mkdir(parents=True, exist_ok=True)
    make_figures(case, report, output_dir)
    (output_dir / "independent_review.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "checks": report["checks"], "output_dir": str(output_dir)}, indent=2))


if __name__ == "__main__":
    main()
