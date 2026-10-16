"""Independent read-only audit of native cohesive-vessel runs.

This module intentionally imports only NumPy and the Python standard library.
It rebuilds mesh ancestry, lumped mass, loads, energy identities, and saved
nodal metrics from the immutable run inputs and histories.  It never modifies
the run or frozen source directories.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np


class AuditInputError(ValueError):
    """Raised when a run or frozen input cannot be trusted or read as declared."""


class AuditArithmeticError(ValueError):
    """Raised for malformed numerical layouts, distinct from a numerical mismatch."""


LEDGER_COLUMNS = {
    "time_s": 0,
    "load_factor_N": 1,
    "interface_work_J": 2,
    "external_work_J": 3,
    "algorithmic_kinetic_J": 4,
    "momentum_x": 5,
    "momentum_y": 6,
    "momentum_z": 7,
    "angular_momentum_x": 8,
    "angular_momentum_y": 9,
    "angular_momentum_z": 10,
    "external_torque_x_Nm": 11,
    "external_torque_y_Nm": 12,
    "external_torque_z_Nm": 13,
    "ma_minus_f_internal_external_max_N": 14,
    "position_recurrence_max_m": 15,
    "velocity_recurrence_max_mps": 16,
    "external_work_minus_total_alg_J": 17,
    "max_damage": 18,
    "min_damage": 19,
    "internal_work_cumulative_J": 20,
    "kinetic_alg_increment_identity_residual_J": 21,
    "internal_work_plus_material_energy_J": 22,
    "native_minus_nodal_kinetic_J": 23,
    "internal_current_torque_x_Nm": 24,
    "internal_current_torque_y_Nm": 25,
    "internal_current_torque_z_Nm": 26,
    "max_contact_traction_Pa": 27,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise AuditInputError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise AuditInputError(f"expected a JSON object in {path}")
    return value


def _relative_path(root: Path, member: str) -> Path:
    if not isinstance(member, str) or not member or Path(member).is_absolute():
        raise AuditInputError(f"invalid relative file path in inventory: {member!r}")
    candidate = (root / member).resolve()
    try:
        candidate.relative_to(root.resolve())
    except ValueError as error:
        raise AuditInputError(f"inventory path escapes its root: {member!r}") from error
    return candidate


def _verify_receipt(run: Path) -> dict[str, Any]:
    receipt = _read_json(run / "receipt.json")
    declared = receipt.get("files")
    if not isinstance(declared, dict):
        raise AuditInputError("receipt.json has no files object")
    verified: dict[str, dict[str, Any]] = {}
    errors: list[str] = []
    for name, expected in declared.items():
        if not isinstance(expected, dict) or "sha256" not in expected or "bytes" not in expected:
            errors.append(f"malformed receipt entry: {name}")
            continue
        try:
            path = _relative_path(run, name)
            actual_bytes = path.stat().st_size
            actual_hash = sha256(path)
        except (OSError, AuditInputError) as error:
            errors.append(f"cannot verify {name}: {error}")
            continue
        valid = actual_bytes == int(expected["bytes"]) and actual_hash == str(expected["sha256"])
        verified[name] = {
            "bytes": actual_bytes,
            "expected_bytes": int(expected["bytes"]),
            "sha256": actual_hash,
            "expected_sha256": str(expected["sha256"]),
            "match": bool(valid),
        }
        if not valid:
            errors.append(f"receipt mismatch: {name}")
    actual_names = {
        path.relative_to(run).as_posix()
        for path in run.rglob("*")
        if path.is_file() and path.name != "receipt.json"
    }
    declared_names = set(declared)
    if actual_names != declared_names:
        errors.append(
            "receipt file inventory differs from run tree "
            f"(missing={sorted(declared_names - actual_names)}, extra={sorted(actual_names - declared_names)})"
        )
    total_bytes = sum(item["bytes"] for item in verified.values())
    if int(receipt.get("persistent_bytes_before_receipt", -1)) != total_bytes:
        errors.append("receipt persistent byte total does not equal the verified file sizes")
    if errors:
        raise AuditInputError("; ".join(errors))
    return {
        "status": receipt.get("status"),
        "file_count": len(verified),
        "verified_bytes": total_bytes,
        "all_sha256_and_byte_counts_match": True,
        "files": verified,
    }


def _verify_frozen_inputs(run: Path, frozen: Path) -> dict[str, Any]:
    freeze = _read_json(frozen / "source_freeze.json")
    case = _read_json(frozen / "case.json")
    started = _read_json(run / "started.json")
    result = _read_json(run / "result.json")
    files = freeze.get("files")
    inventory = started.get("source_input_inventory")
    if not isinstance(files, dict) or not isinstance(inventory, dict):
        raise AuditInputError("frozen source or started.json is missing a file inventory")

    errors: list[str] = []
    frozen_verified: dict[str, Any] = {}
    for name, expected in files.items():
        path = _relative_path(frozen, name)
        if not path.is_file():
            errors.append(f"frozen input is missing: {name}")
            continue
        actual = {"sha256": sha256(path), "bytes": path.stat().st_size}
        matches = actual == {"sha256": expected.get("sha256"), "bytes": expected.get("bytes")}
        frozen_verified[name] = {**actual, "match": bool(matches)}
        if not matches:
            errors.append(f"frozen source hash or byte count mismatch: {name}")

    current_top_level = {
        path.name for path in frozen.iterdir() if path.is_file()
    }
    expected_started_names = set(inventory)
    if current_top_level != expected_started_names:
        errors.append(
            "frozen top-level inventory differs from started.json "
            f"(missing={sorted(expected_started_names-current_top_level)}, extra={sorted(current_top_level-expected_started_names)})"
        )
    started_verified: dict[str, Any] = {}
    for name, expected_hash in inventory.items():
        path = _relative_path(frozen, name)
        if not path.is_file():
            errors.append(f"started source inventory file is missing: {name}")
            continue
        actual_hash = sha256(path)
        started_verified[name] = {"sha256": actual_hash, "expected_sha256": expected_hash, "match": actual_hash == expected_hash}
        if actual_hash != expected_hash:
            errors.append(f"run-start source inventory mismatch: {name}")

    case_hashes = case.get("input_hashes")
    if not isinstance(case_hashes, dict):
        errors.append("case.json has no input_hashes object")
        case_hashes = {}
    case_inputs_verified: dict[str, Any] = {}
    for name, expected_hash in case_hashes.items():
        path = _relative_path(frozen, name)
        if not path.is_file():
            errors.append(f"case input is missing: {name}")
            continue
        actual_hash = sha256(path)
        case_inputs_verified[name] = {"sha256": actual_hash, "expected_sha256": expected_hash, "match": actual_hash == expected_hash}
        if actual_hash != expected_hash:
            errors.append(f"case input hash mismatch: {name}")
    if result.get("source_inputs_unchanged") is not True:
        errors.append("run result does not affirm source_inputs_unchanged")
    if errors:
        raise AuditInputError("; ".join(errors))
    return {
        "frozen_file_count": len(frozen_verified),
        "run_start_inventory_count": len(started_verified),
        "case_input_count": len(case_inputs_verified),
        "frozen_bytes_match_source_freeze": True,
        "run_start_inventory_matches_frozen_bytes": True,
        "case_input_hashes_match_frozen_bytes": True,
        "result_source_inputs_unchanged_claim_rechecked": True,
        "frozen_files": frozen_verified,
        "run_start_inventory": started_verified,
        "case_inputs": case_inputs_verified,
    }


def _load_history_array(path: Path) -> np.ndarray:
    try:
        array = np.load(path, mmap_mode="r", allow_pickle=False)
    except Exception as error:
        raise AuditInputError(f"history array is not readable as a standard NPY file: {path.name}: {error}") from error
    if not isinstance(array, np.ndarray):
        raise AuditInputError(f"history member is not a single NPY array: {path.name}")
    if array.dtype.hasobject:
        raise AuditInputError(f"object arrays are not accepted: {path.name}")
    return array


def _load_initial_state(path: Path) -> dict[str, np.ndarray]:
    try:
        with np.load(path, allow_pickle=False) as archive:
            result = {key: archive[key] for key in archive.files}
    except Exception as error:
        raise AuditInputError(f"initial_state.npz is unreadable: {error}") from error
    if any(value.dtype.hasobject for value in result.values()):
        raise AuditInputError("initial_state.npz contains an object array")
    return result


def _read_msh22(path: Path) -> dict[str, Any]:
    """Read only ASCII Gmsh 2.2 nodes and linear triangle/tetra connectivity."""
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as error:
        raise AuditInputError(f"cannot read frozen MSH: {error}") from error
    sections: dict[str, list[str]] = {}
    cursor = 0
    while cursor < len(lines):
        line = lines[cursor].strip()
        if line.startswith("$") and not line.startswith("$End"):
            section = line[1:]
            end_marker = "$End" + section
            cursor += 1
            rows: list[str] = []
            while cursor < len(lines) and lines[cursor].strip() != end_marker:
                rows.append(lines[cursor].strip())
                cursor += 1
            if cursor == len(lines):
                raise AuditInputError(f"unterminated MSH section {section}")
            sections[section] = rows
        cursor += 1
    mesh_format = sections.get("MeshFormat", [])
    if not mesh_format or mesh_format[0].split()[:2] != ["2.2", "0"]:
        raise AuditInputError("frozen mesh is not ASCII Gmsh MSH 2.2")

    physical_names: dict[tuple[int, int], str] = {}
    names = sections.get("PhysicalNames", [])
    if names:
        if len(names) != int(names[0]) + 1:
            raise AuditInputError("MSH physical-name count mismatch")
        for row in names[1:]:
            fields = row.split(maxsplit=2)
            if len(fields) != 3:
                raise AuditInputError("malformed MSH physical name")
            physical_names[(int(fields[0]), int(fields[1]))] = json.loads(fields[2])

    node_rows = sections.get("Nodes", [])
    if not node_rows or len(node_rows) != int(node_rows[0]) + 1:
        raise AuditInputError("MSH node count mismatch")
    node_tags: list[int] = []
    coordinates: list[list[float]] = []
    for row in node_rows[1:]:
        values = row.split()
        if len(values) != 4:
            raise AuditInputError("malformed MSH node row")
        node_tags.append(int(values[0]))
        coordinates.append([float(value) for value in values[1:]])
    if len(set(node_tags)) != len(node_tags):
        raise AuditInputError("duplicate node tags in frozen MSH")
    tag_to_row = {tag: index for index, tag in enumerate(node_tags)}

    element_rows = sections.get("Elements", [])
    if not element_rows or len(element_rows) != int(element_rows[0]) + 1:
        raise AuditInputError("MSH element count mismatch")
    tetrahedra: list[list[int]] = []
    tetrahedron_tags: list[int] = []
    triangles: dict[str, list[list[int]]] = {"interface": [], "exterior": []}
    for row in element_rows[1:]:
        values = [int(value) for value in row.split()]
        if len(values) < 3:
            raise AuditInputError("malformed MSH element row")
        element_tag, element_type, tag_count = values[:3]
        if len(values) < 3 + tag_count:
            raise AuditInputError(f"malformed MSH element tags at {element_tag}")
        element_tags = values[3:3 + tag_count]
        node_ids = values[3 + tag_count:]
        if element_type == 4:
            if len(node_ids) != 4:
                raise AuditInputError(f"nonlinear/non-T4 volume element at {element_tag}")
            tetrahedra.append([tag_to_row[tag] for tag in node_ids])
            tetrahedron_tags.append(element_tag)
        elif element_type == 2:
            if len(node_ids) != 3:
                raise AuditInputError(f"nonlinear/non-Tri3 surface element at {element_tag}")
            physical_name = physical_names.get((2, element_tags[0])) if element_tags else None
            if physical_name in triangles:
                triangles[physical_name].append([tag_to_row[tag] for tag in node_ids])

    result = {
        "nodes_m": np.asarray(coordinates, dtype=np.float64),
        "node_tags": np.asarray(node_tags, dtype=np.int64),
        "tets": np.asarray(tetrahedra, dtype=np.int64),
        "tet_ids": np.asarray(tetrahedron_tags, dtype=np.int64),
        "triangles": {name: np.asarray(rows, dtype=np.int64) for name, rows in triangles.items()},
    }
    if not np.isfinite(result["nodes_m"]).all() or not len(result["tets"]):
        raise AuditInputError("frozen MSH has non-finite nodes or no T4 elements")
    return result


def qp_trapezoid_work(opening_m: Any, traction_pa: Any, qp_area_m2: Any) -> np.ndarray:
    """Return cumulative all-QP trapezoidal traction work, including time zero."""
    opening = np.asarray(opening_m, dtype=np.float64)
    traction = np.asarray(traction_pa, dtype=np.float64)
    area = np.asarray(qp_area_m2, dtype=np.float64)
    if opening.ndim != 3 or opening.shape[-1] != 3 or traction.shape != opening.shape:
        raise ValueError("opening and traction must have matching shape (steps, QP, 3)")
    if area.shape != (opening.shape[1],) or len(opening) == 0:
        raise ValueError("one positive quadrature area is required per cohesive point")
    if np.any(area <= 0) or not np.isfinite(area).all():
        raise ValueError("quadrature areas must be finite and positive")
    increments = 0.5 * (traction[1:] + traction[:-1]) * np.diff(opening, axis=0)
    increments = np.einsum("tq,q->t", increments.sum(axis=2), area)
    return np.r_[0.0, np.cumsum(increments)]


def qp_trapezoid_work_by_point(opening_m: Any, traction_pa: Any, qp_area_m2: Any) -> np.ndarray:
    """Return cumulative trapezoidal work per quadrature point, shape (steps, QP)."""
    opening = np.asarray(opening_m, dtype=np.float64)
    traction = np.asarray(traction_pa, dtype=np.float64)
    area = np.asarray(qp_area_m2, dtype=np.float64)
    if opening.ndim != 3 or opening.shape[-1] != 3 or traction.shape != opening.shape:
        raise ValueError("opening and traction must have matching shape (steps, QP, 3)")
    if area.shape != (opening.shape[1],) or len(opening) == 0:
        raise ValueError("one positive quadrature area is required per cohesive point")
    if np.any(area <= 0) or not np.isfinite(area).all():
        raise ValueError("quadrature areas must be finite and positive")
    cumulative = np.zeros((len(opening), opening.shape[1]), dtype=np.float64)
    # Integrate one interval at a time. Contact-opening histories can be large;
    # avoid materializing another full (steps, QP, 3) increment tensor.
    for step in range(len(opening) - 1):
        mean_traction = 0.5 * (traction[step] + traction[step + 1])
        delta_opening = opening[step + 1] - opening[step]
        increment = np.einsum("qi,qi,q->q", mean_traction, delta_opening, area)
        cumulative[step + 1] = cumulative[step] + increment
    return cumulative


def quadrature_contact_energy_by_point(
    contact_opening_m: Any, contact_traction_pa: Any, qp_area_m2: Any
) -> np.ndarray:
    """Independently compute 1/2 A tc·c at every frame and quadrature point."""
    opening = np.asarray(contact_opening_m, dtype=np.float64)
    traction = np.asarray(contact_traction_pa, dtype=np.float64)
    area = np.asarray(qp_area_m2, dtype=np.float64)
    if opening.ndim != 3 or opening.shape[-1] != 3 or traction.shape != opening.shape:
        raise ValueError("contact opening and traction must have matching shape (steps, QP, 3)")
    if area.shape != (opening.shape[1],) or len(opening) == 0:
        raise ValueError("one positive quadrature area is required per contact point")
    if np.any(area <= 0) or not np.isfinite(area).all():
        raise ValueError("quadrature areas must be finite and positive")
    return 0.5 * np.einsum("sqi,sqi,q->sq", opening, traction, area)


def interface_work_decomposition(
    opening_m: Any,
    traction_pa: Any,
    contact_opening_m: Any,
    contact_traction_pa: Any,
    qp_area_m2: Any,
) -> np.ndarray:
    """Return cumulative [Wfull, Wcoh, Wc, Wcross] independently at each frame."""
    opening = np.asarray(opening_m, dtype=np.float64)
    traction = np.asarray(traction_pa, dtype=np.float64)
    contact_opening = np.asarray(contact_opening_m, dtype=np.float64)
    contact_traction = np.asarray(contact_traction_pa, dtype=np.float64)
    if not (opening.shape == traction.shape == contact_opening.shape == contact_traction.shape):
        raise ValueError("cohesive/contact openings and tractions must share shape (steps, QP, 3)")
    area = np.asarray(qp_area_m2, dtype=np.float64)
    if opening.ndim != 3 or opening.shape[-1] != 3 or area.shape != (opening.shape[1],):
        raise ValueError("opening histories must be (steps, QP, 3) with one area per QP")
    if not len(opening) or np.any(area <= 0.0) or not np.isfinite(area).all():
        raise ValueError("interface work requires nonempty histories and finite positive areas")
    cumulative = np.zeros((len(opening), 4), dtype=np.float64)
    for step in range(len(opening) - 1):
        cohesive_increment = np.einsum(
            "qi,qi,q->", 0.5 * (traction[step] + traction[step + 1]),
            opening[step + 1] - opening[step], area,
        )
        contact_increment = np.einsum(
            "qi,qi,q->", 0.5 * (contact_traction[step] + contact_traction[step + 1]),
            contact_opening[step + 1] - contact_opening[step], area,
        )
        full_increment = np.einsum(
            "qi,qi,q->",
            0.5 * (
                traction[step] + contact_traction[step]
                + traction[step + 1] + contact_traction[step + 1]
            ),
            (opening[step + 1] - opening[step])
            + (contact_opening[step + 1] - contact_opening[step]),
            area,
        )
        cumulative[step + 1, 1] = cumulative[step, 1] + cohesive_increment
        cumulative[step + 1, 2] = cumulative[step, 2] + contact_increment
        cumulative[step + 1, 0] = cumulative[step, 0] + full_increment
        cumulative[step + 1, 3] = cumulative[step + 1, 0] - cumulative[step + 1, 1] - cumulative[step + 1, 2]
    return cumulative


def qp_failure_steps(
    damage: Any,
    cohesive_traction_pa: Any,
    contact_traction_pa: Any,
    damage_tolerance: float,
    traction_zero_tolerance_pa: float,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Find first full-damage/zero-total-traction step and verify it stays failed."""
    damage_array = np.asarray(damage, dtype=np.float64)
    cohesive = np.asarray(cohesive_traction_pa, dtype=np.float64)
    contact = np.asarray(contact_traction_pa, dtype=np.float64)
    if damage_array.ndim != 2 or cohesive.shape != contact.shape or cohesive.shape != (*damage_array.shape, 3):
        raise ValueError("damage and traction histories must have shapes (steps, QP) and (steps, QP, 3)")
    if not (0.0 <= damage_tolerance < 1.0 and traction_zero_tolerance_pa >= 0.0):
        raise ValueError("invalid damage or traction failure tolerance")
    total_traction_norm = np.linalg.norm(cohesive + contact, axis=2)
    failure_condition = (damage_array >= 1.0 - damage_tolerance) & (total_traction_norm <= traction_zero_tolerance_pa)
    has_failure = np.any(failure_condition, axis=0)
    first = np.where(has_failure, np.argmax(failure_condition, axis=0), -1).astype(np.int64)
    time = np.arange(len(damage_array), dtype=np.int64)[:, None]
    after_failure = (first[None, :] >= 0) & (time >= first[None, :])
    all_qps_failed = bool(np.all(has_failure))
    max_after = float(np.max(total_traction_norm[after_failure])) if np.any(after_failure) else 0.0
    stays_zero = all_qps_failed and max_after <= traction_zero_tolerance_pa
    return first, {
        "all_qps_reached_full_damage_and_zero_total_traction": all_qps_failed,
        "post_failure_total_traction_max_norm_Pa": max_after,
        "post_failure_total_traction_stays_zero": bool(stays_zero),
    }


def reconstruct_unit_load(loading: dict[str, Any], node_count: int) -> np.ndarray:
    """Rebuild signed x-force weights from the frozen patch node IDs and weights."""
    unit = np.zeros((node_count, 3), dtype=np.float64)
    claimed: set[int] = set()
    patches = loading.get("patches")
    if not isinstance(patches, list) or len(patches) != 2:
        raise AuditInputError("loading.json must describe the two signed load patches")
    signs = set()
    for patch in patches:
        sign = int(patch.get("sign", 0))
        ids = np.asarray(patch.get("native_node_ids", []), dtype=np.int64)
        weights = np.asarray(patch.get("weights", []), dtype=np.float64)
        if sign not in (-1, 1) or sign in signs:
            raise AuditInputError("loading patch signs must be exactly -1 and +1")
        signs.add(sign)
        if ids.ndim != 1 or weights.shape != ids.shape or not len(ids):
            raise AuditInputError("loading patch node/weight lengths disagree")
        if np.any(ids < 0) or np.any(ids >= node_count) or len(np.unique(ids)) != len(ids):
            raise AuditInputError("loading patch node IDs are repeated or outside the native mesh")
        if claimed.intersection(map(int, ids)):
            raise AuditInputError("positive and negative loading patches overlap")
        claimed.update(map(int, ids))
        if not np.isfinite(weights).all() or np.any(weights <= 0) or not np.isclose(weights.sum(), 1.0, rtol=0, atol=1e-12):
            raise AuditInputError("each loading patch must have positive unit-sum weights")
        unit[ids, 0] = sign * weights
        reference_line = np.asarray(patch.get("reference_force_line_yz_m", []), dtype=np.float64)
        if reference_line.shape != (2,):
            raise AuditInputError("loading reference_force_line_yz_m must have two entries")
    if signs != {-1, 1}:
        raise AuditInputError("loading.json must contain one patch for each sign")
    return unit


def load_factor(times_s: Any, duration_s: float, peak_force_n: float) -> np.ndarray:
    times = np.asarray(times_s, dtype=np.float64)
    if duration_s <= 0 or peak_force_n <= 0 or not np.isfinite([duration_s, peak_force_n]).all():
        raise ValueError("pulse duration and peak force must be finite and positive")
    factors = np.zeros_like(times)
    active = (times > 0.0) & (times < duration_s)
    factors[active] = peak_force_n * np.sin(np.pi * times[active] / duration_s) ** 2
    return factors


def reconstruct_external_work(
    patch_displacement_x_m: Any,
    unit_external_load_x: Any,
    times_s: Any,
    duration_s: float,
    peak_force_n: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Rebuild nodal patch work from saved patch displacement and force profile."""
    displacement = np.asarray(patch_displacement_x_m, dtype=np.float64)
    unit = np.asarray(unit_external_load_x, dtype=np.float64)
    times = np.asarray(times_s, dtype=np.float64)
    if displacement.ndim != 2 or unit.shape != (displacement.shape[1],) or times.shape != (len(displacement),):
        raise ValueError("patch displacement, active nodal load, and time arrays disagree")
    factors = load_factor(times, duration_s, peak_force_n)
    force = factors[:, None] * unit[None, :]
    increments = 0.5 * (force[1:] + force[:-1]) * np.diff(displacement, axis=0)
    cumulative_work = np.r_[0.0, np.sum(increments, axis=1).cumsum()]
    return cumulative_work, factors


def tetra_lumped_mass_and_strain_energy(
    rest_nodes_m: Any,
    tets: Any,
    displacement_m: Any,
    rho_kg_m3: float,
    young_pa: float,
    poisson: float,
) -> tuple[np.ndarray, float, float]:
    """Rebuild T4 lumped masses, total volume, and linear-elastic small-strain U."""
    rest = np.asarray(rest_nodes_m, dtype=np.float64)
    cells = np.asarray(tets, dtype=np.int64)
    displacement = np.asarray(displacement_m, dtype=np.float64)
    if rest.ndim != 2 or rest.shape[1] != 3 or cells.ndim != 2 or cells.shape[1] != 4:
        raise ValueError("rest_nodes_m and tets must have shapes (N, 3) and (T, 4)")
    if displacement.shape != rest.shape:
        raise ValueError("displacement_m must match rest_nodes_m")
    if not (rho_kg_m3 > 0 and young_pa > 0 and -1.0 < poisson < 0.5):
        raise ValueError("invalid isotropic material constants")
    x = rest[cells]
    j0 = np.swapaxes(x[:, 1:] - x[:, :1], 1, 2)
    det_j0 = np.linalg.det(j0)
    if np.any(det_j0 <= 0) or not np.isfinite(det_j0).all():
        raise ValueError("T4 rest connectivity must have finite positive orientation")
    volumes = det_j0 / 6.0
    masses = np.zeros(len(rest), dtype=np.float64)
    nodal_mass = rho_kg_m3 * volumes / 4.0
    for local_node in range(4):
        np.add.at(masses, cells[:, local_node], nodal_mass)

    u = displacement[cells]
    du = np.swapaxes(u[:, 1:] - u[:, :1], 1, 2)
    deformation_gradient = du @ np.linalg.inv(j0)
    strain = 0.5 * (deformation_gradient + np.swapaxes(deformation_gradient, 1, 2))
    trace = np.trace(strain, axis1=1, axis2=2)
    mu = young_pa / (2.0 * (1.0 + poisson))
    lam = young_pa * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson))
    energy_density = 0.5 * lam * trace**2 + mu * np.einsum("tij,tij->t", strain, strain)
    return masses, float(volumes.sum()), float(np.dot(volumes, energy_density))


def central_difference_residuals(
    u_previous: Any,
    u_current: Any,
    v_previous: Any,
    v_current: Any,
    a_previous: Any,
    a_current: Any,
    dt_s: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return position and velocity recurrence residuals for one adjacent step."""
    up = np.asarray(u_previous, dtype=np.float64)
    uc = np.asarray(u_current, dtype=np.float64)
    vp = np.asarray(v_previous, dtype=np.float64)
    vc = np.asarray(v_current, dtype=np.float64)
    ap = np.asarray(a_previous, dtype=np.float64)
    ac = np.asarray(a_current, dtype=np.float64)
    if not (up.shape == uc.shape == vp.shape == vc.shape == ap.shape == ac.shape) or dt_s <= 0:
        raise ValueError("central-difference state shapes must match and dt must be positive")
    position_residual = uc - up - dt_s * vp - 0.5 * dt_s**2 * ap
    velocity_residual = vc - vp - 0.5 * dt_s * (ap + ac)
    return position_residual, velocity_residual


def _prepare_t4(rest: np.ndarray, tets: np.ndarray, rho: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    points = rest[tets]
    j0 = np.swapaxes(points[:, 1:] - points[:, :1], 1, 2)
    determinant = np.linalg.det(j0)
    if np.any(determinant <= 0.0) or not np.isfinite(determinant).all():
        raise AuditArithmeticError("native T4 rest mesh has nonpositive or non-finite signed volume")
    volumes = determinant / 6.0
    inverse_j0 = np.linalg.inv(j0)
    masses = np.zeros(len(rest), dtype=np.float64)
    element_mass = rho * volumes / 4.0
    for local_node in range(4):
        np.add.at(masses, tets[:, local_node], element_mass)
    return volumes, inverse_j0, masses, determinant


def _bulk_t4_energy(
    tets: np.ndarray,
    displacement: np.ndarray,
    volumes: np.ndarray,
    inverse_j0: np.ndarray,
    young_pa: float,
    poisson: float,
) -> float:
    element_u = displacement[tets]
    displacement_jacobian = np.swapaxes(element_u[:, 1:] - element_u[:, :1], 1, 2)
    grad_u = displacement_jacobian @ inverse_j0
    strain = 0.5 * (grad_u + np.swapaxes(grad_u, 1, 2))
    trace = np.trace(strain, axis1=1, axis2=2)
    mu = young_pa / (2.0 * (1.0 + poisson))
    lam = young_pa * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson))
    density = 0.5 * lam * trace**2 + mu * np.einsum("tij,tij->t", strain, strain)
    return float(np.dot(volumes, density))


def _deformation_metrics(
    rest: np.ndarray,
    tets: np.ndarray,
    displacement: np.ndarray,
    inverse_j0: np.ndarray,
) -> dict[str, float]:
    points = rest[tets] + displacement[tets]
    j = np.swapaxes(points[:, 1:] - points[:, :1], 1, 2)
    deformation_gradient = j @ inverse_j0
    strain = 0.5 * (deformation_gradient + np.swapaxes(deformation_gradient, 1, 2)) - np.eye(3)
    principal_strain = np.linalg.eigvalsh(strain)
    left, singular_values, right = np.linalg.svd(deformation_gradient, full_matrices=False)
    rotation = left @ right
    cosine = np.clip((np.trace(rotation, axis1=1, axis2=2) - 1.0) / 2.0, -1.0, 1.0)
    return {
        "max_abs_principal_small_strain": float(np.max(np.abs(principal_strain))),
        "min_det_F": float(np.linalg.det(deformation_gradient).min()),
        "max_polar_rotation_rad": float(np.arccos(cosine).max()),
        "min_principal_stretch": float(singular_values.min()),
        "max_principal_stretch": float(singular_values.max()),
    }


def _maxabs(value: Any) -> float:
    array = np.asarray(value, dtype=np.float64)
    return float(np.max(np.abs(array))) if array.size else 0.0


def _relative_error(error: Any, scale: float) -> float:
    return _maxabs(error) / max(abs(float(scale)), np.finfo(np.float64).tiny)


def _finite_stats(values: Any) -> dict[str, float]:
    array = np.asarray(values, dtype=np.float64)
    return {
        "min": float(np.min(array)),
        "max": float(np.max(array)),
        "max_abs": float(np.max(np.abs(array))),
    }


def _row_set(values: np.ndarray) -> set[tuple[int, ...]]:
    return {tuple(sorted(map(int, row))) for row in np.asarray(values)}


def _history_arrays(run: Path, include_contact_work: bool = False) -> tuple[dict[str, np.ndarray], np.ndarray]:
    history_dir = run / "history"
    names = (
        "opening_m",
        "traction_Pa",
        "damage",
        "delta_max_m",
        "contact_traction_Pa",
        "energy_J",
        "ledger",
        "patch_displacement_x_m",
        "node_state",
    )
    arrays = {name: _load_history_array(history_dir / f"{name}.npy") for name in names}
    if include_contact_work:
        for name in ("contact_opening_m", "contact_energy_check", "interface_work_check"):
            arrays[name] = _load_history_array(history_dir / f"{name}.npy")
    saved_steps = _load_history_array(history_dir / "node_saved_steps.npy")
    return arrays, saved_steps


def _expect_shape(name: str, array: np.ndarray, shape: tuple[int, ...]) -> None:
    if array.shape != shape:
        raise AuditInputError(f"{name} has shape {array.shape}, expected {shape}")


def _all_finite(array: np.ndarray) -> bool:
    if array.dtype.kind not in "fc":
        return True
    if array.ndim == 0:
        return bool(np.isfinite(array))
    # Scan large memmaps in time slabs to keep peak memory bounded.
    slab = max(1, min(64, array.shape[0]))
    return all(bool(np.isfinite(array[start:start + slab]).all()) for start in range(0, array.shape[0], slab))


def _audit_native_run(run: Path, frozen: Path, integrity: dict[str, Any], receipt: dict[str, Any]) -> dict[str, Any]:
    result = _read_json(run / "result.json")
    completion = _read_json(run / "completion.json")
    case = _read_json(frozen / "case.json")
    loading = _read_json(run / "loading.json")
    initial = _load_initial_state(run / "initial_state.npz")
    interface_case = case["interface"]
    try:
        contact_penalty = float(interface_case.get("penalty_Pa_m", 0.0))
    except (TypeError, ValueError, OverflowError) as error:
        raise AuditInputError("interface penalty_Pa_m must be a finite nonnegative scalar") from error
    if not np.isfinite(contact_penalty) or contact_penalty < 0.0:
        raise AuditInputError("interface penalty_Pa_m must be finite and nonnegative")
    compression_contact_enabled = contact_penalty > 0.0
    history, node_saved_steps = _history_arrays(run, include_contact_work=compression_contact_enabled)

    steps = int(result.get("solver_steps", -1))
    dt = float(result.get("dt_s", 0.0))
    if steps <= 0 or not np.isfinite(dt) or dt <= 0:
        raise AuditInputError("result.json does not declare a positive solver step count and dt")
    frame_count = steps + 1
    rest = np.asarray(initial.get("rest_nodes_m"), dtype=np.float64)
    tets = np.asarray(initial.get("tets"), dtype=np.int64)
    cohesive = np.asarray(initial.get("cohesive"), dtype=np.int64)
    ancestry = np.asarray(initial.get("native_to_source"), dtype=np.int64)
    qp_area = np.asarray(initial.get("qp_area_m2"), dtype=np.float64)
    node_count = len(rest)
    qp_count = len(qp_area)
    if rest.shape != (node_count, 3) or tets.ndim != 2 or tets.shape[1] != 4 or cohesive.ndim != 2 or cohesive.shape[1] != 6:
        raise AuditInputError("initial_state.npz has malformed node, T4, or cohesive connectivity")
    if ancestry.shape != (node_count,) or np.any(ancestry < 0):
        raise AuditInputError("native_to_source ancestry has an invalid shape or negative row")
    if result.get("counts") != [node_count, len(tets), len(cohesive)]:
        raise AuditInputError("result.json element/node counts disagree with initial_state.npz")
    if completion.get("steps_completed") != steps or completion.get("qp_frames_written") != frame_count:
        raise AuditInputError("completion.json does not confirm all requested solver/QP frames")

    expected_shapes = {
        "opening_m": (frame_count, qp_count, 3),
        "traction_Pa": (frame_count, qp_count, 3),
        "damage": (frame_count, qp_count),
        "delta_max_m": (frame_count, qp_count),
        "contact_traction_Pa": (frame_count, qp_count, 3),
        "energy_J": (frame_count, 5),
        "ledger": (frame_count, len(LEDGER_COLUMNS)),
    }
    if compression_contact_enabled:
        expected_shapes.update({
            "contact_opening_m": (frame_count, qp_count, 3),
            "contact_energy_check": (frame_count, 3),
            "interface_work_check": (frame_count, 4),
        })
    for name, shape in expected_shapes.items():
        _expect_shape(name, history[name], shape)
    active_nodes = np.flatnonzero(np.asarray(initial["unit_external_load"])[:, 0])
    _expect_shape("patch_displacement_x_m", history["patch_displacement_x_m"], (frame_count, len(active_nodes)))
    _expect_shape("node_saved_steps", node_saved_steps, (len(history["node_state"]),))
    _expect_shape("node_state", history["node_state"], (len(node_saved_steps), node_count, 4, 3))
    if len(node_saved_steps) == 0 or np.any(np.diff(node_saved_steps) <= 0) or node_saved_steps[0] != 0 or node_saved_steps[-1] != steps:
        raise AuditInputError("node_saved_steps is incomplete, duplicated, or out of order")
    if any(not _all_finite(array) for array in (*history.values(), node_saved_steps)):
        raise AuditArithmeticError("one or more complete native history arrays contain non-finite values")

    node_stride = int(case["output"]["node_stride"])
    expected_saved_steps = np.unique(np.r_[np.arange(0, steps + 1, node_stride), steps]).astype(np.int64)
    if not np.array_equal(node_saved_steps, expected_saved_steps):
        raise AuditInputError("node_saved_steps do not match the frozen node_stride and final-step contract")
    all_node_steps_observed = np.array_equal(node_saved_steps, np.arange(frame_count, dtype=np.int64))

    msh = _read_msh22(frozen / "mesh.msh")
    source_tets = msh["tets"]
    source_nodes = msh["nodes_m"]
    source_interface = msh["triangles"]["interface"]
    source_exterior = msh["triangles"]["exterior"]
    if not np.array_equal(np.asarray(initial["source_node_tags"], dtype=np.int64), msh["node_tags"]):
        raise AuditInputError("initial-state source node tags disagree with frozen MSH")
    if not np.array_equal(np.asarray(initial["source_tet_ids"], dtype=np.int64), msh["tet_ids"]):
        raise AuditInputError("initial-state source T4 IDs disagree with frozen MSH")
    if len(ancestry) != len(rest) or np.any(ancestry >= len(source_nodes)):
        raise AuditInputError("native-to-source ancestry references a missing source node")
    node_lineage_exact = np.array_equal(rest, source_nodes[ancestry])
    mapped_tets = ancestry[tets]
    tet_lineage_exact = np.array_equal(mapped_tets, source_tets)
    exterior_lineage_exact = _row_set(np.asarray(initial["source_exterior_faces"])) == _row_set(source_exterior)
    if not (node_lineage_exact and tet_lineage_exact and exterior_lineage_exact):
        raise AuditInputError("native nodes/T4s/exterior faces do not map exactly to the frozen MSH rows")
    pairs = np.asarray(initial["cohesive_pairs"], dtype=np.int64)
    if pairs.shape != (len(cohesive), 2, 3):
        raise AuditInputError("cohesive_pairs has an unexpected shape")
    pair_source_a = ancestry[pairs[:, 0]]
    pair_source_b = ancestry[pairs[:, 1]]
    cohesive_source_faces = [tuple(sorted(map(int, face))) for face in pair_source_a]
    cohesive_source_faces_b = [tuple(sorted(map(int, face))) for face in pair_source_b]
    interface_source_faces = _row_set(source_interface)
    cohesive_lineage_exact = (
        set(cohesive_source_faces) == interface_source_faces
        and len(set(cohesive_source_faces)) == len(cohesive_source_faces)
        and cohesive_source_faces == cohesive_source_faces_b
        and np.array_equal(pairs[:, 0], cohesive[:, :3])
    )
    native_triangle_area_a = 0.5 * np.linalg.norm(
        np.cross(
            rest[cohesive[:, 1]] - rest[cohesive[:, 0]],
            rest[cohesive[:, 2]] - rest[cohesive[:, 0]],
        ),
        axis=1,
    )
    native_triangle_area_b = 0.5 * np.linalg.norm(
        np.cross(
            rest[cohesive[:, 4]] - rest[cohesive[:, 3]],
            rest[cohesive[:, 5]] - rest[cohesive[:, 3]],
        ),
        axis=1,
    )
    qp_area_match = np.allclose(
        qp_area.reshape(len(cohesive), 3),
        np.repeat((native_triangle_area_a / 3.0)[:, None], 3, axis=1),
        rtol=1e-12,
        atol=1e-18,
    ) and np.allclose(native_triangle_area_a, native_triangle_area_b, rtol=1e-12, atol=1e-18)
    if not cohesive_lineage_exact or not qp_area_match:
        raise AuditInputError("cohesive side pairing or native QP area weights disagree with the frozen interface")

    bulk = case["bulk_material"]
    interface = case["interface"]
    loading_case = case["loading"]
    rho = float(bulk["rho_kg_m3"])
    young = float(bulk["E_pa"])
    poisson = float(bulk["nu"])
    gc = float(interface["Gc_J_m2"])
    volumes, inverse_j0, rebuilt_mass, determinants = _prepare_t4(rest, tets, rho)
    stored_mass3 = np.asarray(initial["mass_kg"], dtype=np.float64)
    if stored_mass3.shape != (node_count, 3):
        raise AuditInputError("native mass array must have shape (nodes, 3)")
    mass_rel_error = float(
        np.max(np.abs(stored_mass3 - rebuilt_mass[:, None]))
        / max(float(rebuilt_mass.max()), np.finfo(np.float64).tiny)
    )
    mass_matches_t4 = mass_rel_error <= 1e-12 and np.allclose(
        stored_mass3[:, 0], stored_mass3[:, 1], rtol=0.0, atol=1e-20
    ) and np.allclose(stored_mass3[:, 0], stored_mass3[:, 2], rtol=0.0, atol=1e-20)

    unit_load = reconstruct_unit_load(loading, node_count)
    stored_unit_load = np.asarray(initial["unit_external_load"], dtype=np.float64)
    if stored_unit_load.shape != unit_load.shape:
        raise AuditInputError("stored unit_external_load has an unexpected shape")
    unit_load_matches = np.allclose(unit_load, stored_unit_load, rtol=0.0, atol=2e-15)
    loading_line_errors: list[float] = []
    for patch in loading["patches"]:
        ids = np.asarray(patch["native_node_ids"], dtype=np.int64)
        weights = np.asarray(patch["weights"], dtype=np.float64)
        claimed_line = np.asarray(patch["reference_force_line_yz_m"], dtype=np.float64)
        loading_line_errors.append(float(np.max(np.abs(weights @ rest[ids, 1:] - claimed_line))))
    reference_resultant = unit_load.sum(axis=0)
    reference_torque = np.cross(rest, unit_load).sum(axis=0)
    load_factor_history = load_factor(
        np.arange(frame_count, dtype=np.float64) * dt,
        float(loading_case["pulse_duration_s"]),
        float(loading_case["peak_force_N"]),
    )
    stored_times = np.arange(frame_count, dtype=np.float64) * dt
    ledger = history["ledger"]
    energy = history["energy_J"]
    times_match = np.allclose(ledger[:, 0], stored_times, rtol=0.0, atol=max(1e-20, dt * 1e-12))
    load_factor_match = np.allclose(ledger[:, 1], load_factor_history, rtol=1e-12, atol=1e-12)

    reconstructed_external_work, _ = reconstruct_external_work(
        history["patch_displacement_x_m"],
        unit_load[active_nodes, 0],
        stored_times,
        float(loading_case["pulse_duration_s"]),
        float(loading_case["peak_force_N"]),
    )
    external_work_residual = reconstructed_external_work - ledger[:, 3]
    external_work_scale = max(_maxabs(reconstructed_external_work), _maxabs(ledger[:, 3]), 1e-30)
    external_work_relative_error = _maxabs(external_work_residual) / external_work_scale

    damage = history["damage"]
    delta_max = history["delta_max_m"]
    damage_tol = float(case["gates"]["damage_tolerance"])
    damage_steps = np.arange(frame_count, dtype=np.int64)[:, None]
    first_full = np.min(np.where(damage >= 1.0 - damage_tol, damage_steps, frame_count), axis=0)
    first_full = np.where(first_full == frame_count, -1, first_full)
    first_full_saved = _load_history_array(run / "qp_first_full_damage_step.npy")
    _expect_shape("qp_first_full_damage_step", first_full_saved, (qp_count,))
    first_full_match = np.array_equal(first_full, first_full_saved)
    damage_increment = np.diff(damage, axis=0)
    damage_irreversible = bool(np.min(damage_increment) >= -1e-12)
    damage_in_bounds = bool(float(damage.min()) >= -damage_tol and float(damage.max()) <= 1.0 + damage_tol)
    delta_increment = np.diff(delta_max, axis=0)
    delta_irreversible = bool(float(delta_max.min()) >= 0.0 and float(delta_increment.min()) >= -1e-14 * max(float(delta_max.max()), 1e-30))
    all_qp_full = bool(np.all(first_full >= 0) and np.all(damage[-1] >= 1.0 - damage_tol))
    discrete_damaged_area = np.einsum("tq,q->t", damage, qp_area)
    native_dissipated = energy[:, 2]
    d_nonnegative_tol = max(gc * float(qp_area.sum()), 1e-30) * 1e-12
    d_nonnegative = bool(float(native_dissipated.min()) >= -d_nonnegative_tol)
    d_monotone = bool(float(np.diff(native_dissipated).min()) >= -d_nonnegative_tol)
    fracture_energy_scale = max(gc * float(discrete_damaged_area[-1]), 1e-30)
    d_area_residual = native_dissipated - gc * discrete_damaged_area
    d_area_relative_error = _maxabs(d_area_residual) / fracture_energy_scale
    final_d_area_relative_error = abs(float(d_area_residual[-1])) / fracture_energy_scale

    qp_work = qp_trapezoid_work(history["opening_m"], history["traction_Pa"], qp_area)
    native_cohesive_energy = energy[:, 2] + energy[:, 3]
    qp_work_residual = qp_work - native_cohesive_energy
    qp_work_error = _maxabs(qp_work_residual)
    qp_work_relative_fracture = qp_work_error / fracture_energy_scale
    total_energy_from_ledger = ledger[:, 4] + energy[:, 1:].sum(axis=1)
    total_energy_scale = max(_maxabs(total_energy_from_ledger), fracture_energy_scale, 1e-30)
    qp_work_relative_total = qp_work_error / total_energy_scale
    tolerance = {
        "mass_relative": 1e-12,
        "external_work_vs_ledger_relative": 1e-10,
        "cohesive_work_vs_fracture_energy_relative": 1e-2,
        "dissipation_vs_Gc_times_discrete_area_relative": 1e-2,
        "state_and_force_recurrence_relative": 1e-8,
        "energy_identity_relative_to_fracture_energy": 1e-2,
        "bulk_T4_energy_relative_to_native_potential": 1e-8,
    }

    contact = history["contact_traction_Pa"]
    cohesive_traction = history["traction_Pa"]
    traction_zero_tolerance = float(case["gates"]["traction_zero_tolerance_Pa"])
    contact_norm = np.linalg.norm(contact, axis=2)
    contact_max = _maxabs(contact)
    contact_energy_max = _maxabs(energy[:, 4])
    contact_zero = (
        contact_max <= traction_zero_tolerance
        and contact_energy_max <= 1e-12 * max(total_energy_scale, 1e-30)
    )
    first_failure, qp_failure = qp_failure_steps(
        damage, cohesive_traction, contact, damage_tol, traction_zero_tolerance
    )
    qp_after_failure = (first_failure[None, :] >= 0) & (
        np.arange(frame_count, dtype=np.int64)[:, None] >= first_failure[None, :]
    )
    post_failure_contact_max = (
        float(np.max(contact_norm[qp_after_failure])) if np.any(qp_after_failure) else 0.0
    )
    post_failure_contact_zero = bool(
        qp_failure["all_qps_reached_full_damage_and_zero_total_traction"]
        and post_failure_contact_max <= traction_zero_tolerance
    )
    qp_failure.update({
        "first_failure_step_min": int(first_failure[first_failure >= 0].min())
        if np.any(first_failure >= 0) else None,
        "first_failure_step_max": int(first_failure.max()) if np.any(first_failure >= 0) else None,
        "first_failure_steps": first_failure.tolist(),
        "post_failure_contact_traction_max_norm_Pa": post_failure_contact_max,
        "post_failure_contact_traction_stays_zero": post_failure_contact_zero,
        "traction_zero_tolerance_Pa": traction_zero_tolerance,
    })

    contact_work_report: dict[str, Any] = {
        "policy": "ZERO_CONTACT_REQUIRED" if not compression_contact_enabled else "PREBREAK_COMPRESSION_PENALTY",
        "contact_traction_max_abs_Pa": contact_max,
        "native_contact_energy_max_abs_J": contact_energy_max,
        "contact_traction_and_energy_zero": bool(contact_zero) if not compression_contact_enabled else None,
    }
    compression_contact_checks: dict[str, bool | None] = {
        "zero_contact_case_traction_and_energy_are_zero": bool(contact_zero)
        if not compression_contact_enabled else None,
        "compression_contact_energy_matches_native_quadrature": None,
        "compression_contact_work_matches_delta_contact_energy": None,
        "compression_contact_check_arrays_match_independent_reconstruction": None,
        "compression_penalty_within_fixed_normal_Kn_envelope": None,
        "contact_after_breaking_is_disabled": None,
        "contact_traction_remains_zero_after_qp_failure": bool(post_failure_contact_zero),
    }
    contact_stiffness_report: dict[str, Any] = {"status": "NOT_APPLICABLE"}
    if compression_contact_enabled:
        contact_opening = history["contact_opening_m"]
        contact_energy_by_qp = quadrature_contact_energy_by_point(contact_opening, contact, qp_area)
        independent_contact_energy = contact_energy_by_qp.sum(axis=1)
        native_contact_energy = energy[:, 4]
        contact_energy_residual = independent_contact_energy - native_contact_energy
        contact_energy_relative_fracture = _maxabs(contact_energy_residual) / fracture_energy_scale
        contact_energy_relative_total = _maxabs(contact_energy_residual) / total_energy_scale
        contact_energy_nonnegative_tolerance = max(
            fracture_energy_scale, _maxabs(independent_contact_energy), 1e-30
        ) * 1e-12
        contact_energy_nonnegative = bool(
            np.min(contact_energy_by_qp) >= -contact_energy_nonnegative_tolerance
        )

        contact_work_by_qp = qp_trapezoid_work_by_point(contact_opening, contact, qp_area)
        contact_work = contact_work_by_qp.sum(axis=1)
        contact_energy_change = independent_contact_energy - independent_contact_energy[0]
        contact_work_residual = contact_work - contact_energy_change
        contact_work_residual_by_qp = contact_work_by_qp - (
            contact_energy_by_qp - contact_energy_by_qp[0]
        )
        contact_work_relative_fracture = _maxabs(contact_work_residual_by_qp) / fracture_energy_scale
        contact_work_relative_total = _maxabs(contact_work_residual_by_qp) / total_energy_scale

        interface_work = interface_work_decomposition(
            history["opening_m"], cohesive_traction, contact_opening, contact, qp_area
        )
        contact_energy_check_expected = np.column_stack((
            independent_contact_energy, native_contact_energy, contact_energy_residual
        ))
        contact_energy_check_error = _maxabs(
            history["contact_energy_check"] - contact_energy_check_expected
        )
        contact_work_check_error = _maxabs(history["interface_work_check"] - interface_work)
        contact_check_array_scale = max(
            fracture_energy_scale,
            _maxabs(contact_energy_check_expected),
            _maxabs(interface_work),
            1e-30,
        )
        contact_check_arrays_match = (
            max(contact_energy_check_error, contact_work_check_error)
            <= 1e-10 * contact_check_array_scale
        )

        contact_active = contact_norm > traction_zero_tolerance
        branch_switch_mask = contact_active[1:] != contact_active[:-1]
        branch_switch_coords = np.argwhere(branch_switch_mask)
        contact_step_residual_by_qp = np.diff(contact_work_by_qp, axis=0) - np.diff(
            contact_energy_by_qp, axis=0
        )
        if len(branch_switch_coords):
            switched_residuals = contact_step_residual_by_qp[branch_switch_mask]
            worst_switch = int(np.argmax(np.abs(switched_residuals)))
            worst_switch_step = int(branch_switch_coords[worst_switch, 0] + 1)
            worst_switch_qp = int(branch_switch_coords[worst_switch, 1])
            switch_residual_max = float(np.abs(switched_residuals[worst_switch]))
        else:
            worst_switch_step = None
            worst_switch_qp = None
            switch_residual_max = 0.0

        contact_opening_norm = np.linalg.norm(contact_opening, axis=2)
        active_gap = contact_opening_norm[contact_active]
        interface = case["interface"]
        sigma_c = float(interface["sigma_c_Pa"])
        delta0 = float(interface["delta0_m"])
        if not np.isfinite(sigma_c) or not np.isfinite(delta0) or sigma_c <= 0.0 or delta0 <= 0.0:
            raise AuditInputError("compression-contact Kn requires finite positive sigma_c and delta0")
        kn = sigma_c / delta0
        penalty_within_kn = bool(contact_penalty <= kn * (1.0 + 1e-12))
        contact_after_breaking = interface.get("contact_after_breaking", False)
        contact_after_breaking_disabled = contact_after_breaking is False
        beta_kappa_match = interface.get("beta") == 1 and interface.get("kappa") == 1
        contact_stiffness_report = {
            "status": "OBSERVED",
            "penalty_Pa_per_m": contact_penalty,
            "cohesive_normal_stiffness_Kn_Pa_per_m": kn,
            "cohesive_tangential_stiffness_kappa_Kn_Pa_per_m": float(interface["kappa"]) * kn,
            "penalty_le_Kn": penalty_within_kn,
            "penalty_to_Kn_ratio": contact_penalty / kn,
            "beta_is_one": interface.get("beta") == 1,
            "kappa_is_one": interface.get("kappa") == 1,
            "contact_after_breaking": contact_after_breaking,
            "contact_after_breaking_disabled": contact_after_breaking_disabled,
            "stability_scope": "fixed-normal tangent envelope only; not a nonlinear or full-trajectory geometry bound",
        }
        contact_work_report.update({
            "contact_energy_independent_final_J": float(independent_contact_energy[-1]),
            "native_contact_energy_final_J": float(native_contact_energy[-1]),
            "contact_energy_max_abs_difference_J": _maxabs(contact_energy_residual),
            "contact_energy_difference_relative_to_fracture_scale": contact_energy_relative_fracture,
            "contact_energy_difference_relative_to_total_scale": contact_energy_relative_total,
            "contact_energy_min_independent_J": float(np.min(contact_energy_by_qp)),
            "contact_energy_nonnegative": contact_energy_nonnegative,
            "contact_work_final_J": float(contact_work[-1]),
            "contact_work_minus_delta_independent_energy_max_abs_J": _maxabs(contact_work_residual),
            "contact_work_minus_delta_energy_max_abs_per_qp_J": _maxabs(contact_work_residual_by_qp),
            "contact_work_balance_relative_to_fracture_scale": contact_work_relative_fracture,
            "contact_work_balance_relative_to_total_scale": contact_work_relative_total,
            "full_interface_work_final_J": float(interface_work[-1, 0]),
            "cohesive_work_final_J": float(interface_work[-1, 1]),
            "contact_work_from_full_decomposition_final_J": float(interface_work[-1, 2]),
            "cross_work_final_J": float(interface_work[-1, 3]),
            "full_interface_work_history_max_abs_J": _maxabs(interface_work[:, 0]),
            "cross_work_history_max_abs_J": _maxabs(interface_work[:, 3]),
            "contact_energy_check_array_max_abs_difference_J": contact_energy_check_error,
            "interface_work_check_array_max_abs_difference_J": contact_work_check_error,
            "contact_check_array_tolerance_J": 1e-10 * contact_check_array_scale,
            "contact_branch_switch_count": int(len(branch_switch_coords)),
            "contact_branch_switch_residual_max_abs_J": switch_residual_max,
            "contact_branch_switch_residual_relative_to_fracture_scale": switch_residual_max / fracture_energy_scale,
            "worst_contact_branch_switch_step": worst_switch_step,
            "worst_contact_branch_switch_qp": worst_switch_qp,
            "contact_opening_vector_norm_max_m": float(contact_opening_norm.max()),
            "contact_opening_vector_norm_max_when_contact_active_m": float(active_gap.max()) if active_gap.size else 0.0,
            "native_contact_penetration_vector_norm_max_m": float(active_gap.max()) if active_gap.size else 0.0,
            "contact_opening_norm_interpretation": "native penalty contact-opening vector c magnitude; not the full paired-face separation",
            "contact_active_qp_frame_count": int(np.count_nonzero(contact_active)),
        })
        compression_contact_checks.update({
            "compression_contact_energy_matches_native_quadrature": (
                contact_energy_nonnegative
                and contact_energy_relative_fracture <= tolerance["energy_identity_relative_to_fracture_energy"]
            ),
            "compression_contact_work_matches_delta_contact_energy": (
                contact_work_relative_fracture <= tolerance["energy_identity_relative_to_fracture_energy"]
            ),
            "compression_contact_check_arrays_match_independent_reconstruction": contact_check_arrays_match,
            "compression_penalty_within_fixed_normal_Kn_envelope": penalty_within_kn and beta_kappa_match,
            "contact_after_breaking_is_disabled": contact_after_breaking_disabled,
            "contact_traction_remains_zero_after_qp_failure": post_failure_contact_zero,
        })
    else:
        contact_work_report.update({
            "zero_contact_required": True,
            "contact_opening_available": False,
        })

    force_report: dict[str, Any] = {"status": "NOT_OBSERVED", "observed_node_frames": int(len(node_saved_steps))}
    strict_step_report: dict[str, Any] = {"status": "NOT_OBSERVED", "reason": "saved nodal states do not include every solver step"}
    state_reconstruction: dict[str, Any] = {
        "status": "NOT_OBSERVED",
        "saved_node_frames": int(len(node_saved_steps)),
        "solver_steps": steps,
        "all_node_states_each_step": bool(all_node_steps_observed),
    }
    energy_balance: dict[str, Any] = {"status": "NOT_OBSERVED"}
    deformation_aggregate = {
        "max_abs_principal_small_strain": 0.0,
        "min_det_F": float("inf"),
        "max_polar_rotation_rad": 0.0,
        "min_principal_stretch": float("inf"),
        "max_principal_stretch": 0.0,
    }
    frame_quantities: dict[int, dict[str, Any]] = {}
    step_force_errors: list[float] = []
    step_force_scales: list[float] = []
    step_bulk_energy: list[float] = []
    step_kalg: list[float] = []
    step_momentum: list[np.ndarray] = []
    step_angular: list[np.ndarray] = []
    step_external_torque: list[np.ndarray] = []
    step_internal_torque: list[np.ndarray] = []
    ledger_comparison: dict[str, float] = {}
    sampled_pair_gap_max_m = 0.0
    sampled_full_gap_reconstruction_error_m = 0.0
    sampled_gap_max_step: int | None = None

    blocked = np.asarray(initial["blocked"], dtype=bool)
    if blocked.shape != (node_count, 3):
        raise AuditInputError("blocked DOF mask has an unexpected shape")
    all_dofs_free = not bool(blocked.any())
    u_scale = v_scale = a_scale = 0.0
    for frame_index, step in enumerate(node_saved_steps):
        state = np.asarray(history["node_state"][frame_index], dtype=np.float64)
        displacement = state[:, 0, :]
        velocity = state[:, 1, :]
        acceleration = state[:, 2, :]
        internal = state[:, 3, :]
        if compression_contact_enabled:
            current_positions = rest + displacement
            paired_face_gap = (
                current_positions[pairs[:, 0]].mean(axis=1)
                - current_positions[pairs[:, 1]].mean(axis=1)
            )
            qp_full_gap = (
                np.asarray(history["opening_m"][int(step)])
                + np.asarray(history["contact_opening_m"][int(step)])
            ).reshape(len(cohesive), 3, 3).mean(axis=1)
            gap_norms = np.linalg.norm(paired_face_gap, axis=1)
            if float(gap_norms.max()) > sampled_pair_gap_max_m:
                sampled_pair_gap_max_m = float(gap_norms.max())
                sampled_gap_max_step = int(step)
            sampled_full_gap_reconstruction_error_m = max(
                sampled_full_gap_reconstruction_error_m,
                _maxabs(qp_full_gap - paired_face_gap),
            )
        external = unit_load * load_factor_history[int(step)]
        force_residual = rebuilt_mass[:, None] * acceleration - internal - external
        step_force_errors.append(_maxabs(force_residual))
        step_force_scales.append(max(_maxabs(rebuilt_mass[:, None] * acceleration), _maxabs(internal), _maxabs(external), 1e-30))
        kinetic = 0.5 * float(np.sum(rebuilt_mass[:, None] * velocity * velocity))
        kalg = kinetic - dt**2 * float(np.sum(rebuilt_mass[:, None] * acceleration * acceleration)) / 8.0
        momentum = np.sum(rebuilt_mass[:, None] * velocity, axis=0)
        angular = np.cross(rest + displacement, rebuilt_mass[:, None] * velocity).sum(axis=0)
        external_torque = np.cross(rest + displacement, external).sum(axis=0)
        internal_torque = np.cross(rest + displacement, internal).sum(axis=0)
        bulk_u = _bulk_t4_energy(tets, displacement, volumes, inverse_j0, young, poisson)
        deformation = _deformation_metrics(rest, tets, displacement, inverse_j0)
        for key, value in deformation.items():
            if key == "min_det_F" or key == "min_principal_stretch":
                deformation_aggregate[key] = min(deformation_aggregate[key], value)
            else:
                deformation_aggregate[key] = max(deformation_aggregate[key], value)
        step_bulk_energy.append(bulk_u)
        step_kalg.append(kalg)
        step_momentum.append(momentum)
        step_angular.append(angular)
        step_external_torque.append(external_torque)
        step_internal_torque.append(internal_torque)
        frame_quantities[int(step)] = {
            "kinetic_J": kinetic,
            "algorithmic_kinetic_J": kalg,
            "bulk_T4_small_strain_energy_J": bulk_u,
            "momentum_kg_m_s": momentum,
            "angular_momentum_kg_m2_s": angular,
            "external_torque_Nm": external_torque,
            "internal_torque_Nm": internal_torque,
            "force_residual_max_N": _maxabs(force_residual),
            "force_scale_max_N": max(_maxabs(rebuilt_mass[:, None] * acceleration), _maxabs(internal), _maxabs(external), 1e-30),
            "deformation": deformation,
        }
        u_scale = max(u_scale, _maxabs(displacement))
        v_scale = max(v_scale, _maxabs(velocity))
        a_scale = max(a_scale, _maxabs(acceleration))

    if compression_contact_enabled:
        contact_work_report.update({
            "sampled_native_paired_face_centroid_gap_max_norm_m": sampled_pair_gap_max_m,
            "sampled_native_paired_face_gap_max_step": sampled_gap_max_step,
            "sampled_qp_full_opening_vs_paired_face_gap_max_abs_error_m": sampled_full_gap_reconstruction_error_m,
            "paired_face_gap_reconstruction_uses_saved_node_frames_only": True,
        })

    if frame_quantities:
        compare_indices = np.asarray(node_saved_steps, dtype=np.int64)
        frame_kalg = np.asarray(step_kalg)
        frame_bulk = np.asarray(step_bulk_energy)
        frame_momentum = np.asarray(step_momentum)
        frame_angular = np.asarray(step_angular)
        frame_ext_torque = np.asarray(step_external_torque)
        frame_int_torque = np.asarray(step_internal_torque)
        frame_ledger = ledger[compare_indices]
        frame_energy = energy[compare_indices]
        state_reconstruction = {
            "status": "PASS_SAMPLED" if not all_node_steps_observed else "PASS_ALL_STEPS",
            "saved_node_frames": int(len(node_saved_steps)),
            "solver_steps": steps,
            "all_node_states_each_step": bool(all_node_steps_observed),
            "kinetic_energy_max_abs_difference_J": _maxabs(frame_energy[:, 0] - np.asarray([q["kinetic_J"] for q in frame_quantities.values()])),
            "algorithmic_kinetic_max_abs_difference_J": _maxabs(frame_ledger[:, 4] - frame_kalg),
            "bulk_T4_energy_max_abs_difference_J": _maxabs(frame_energy[:, 1] - frame_bulk),
            "momentum_max_abs_difference_kg_m_s": _maxabs(frame_ledger[:, 5:8] - frame_momentum),
            "angular_momentum_max_abs_difference_kg_m2_s": _maxabs(frame_ledger[:, 8:11] - frame_angular),
            "external_torque_max_abs_difference_Nm": _maxabs(frame_ledger[:, 11:14] - frame_ext_torque),
            "internal_torque_max_abs_difference_Nm": _maxabs(frame_ledger[:, 24:27] - frame_int_torque),
            "force_balance_max_abs_N": max(step_force_errors),
            "force_balance_max_relative_to_force_scale": max(
                error / scale for error, scale in zip(step_force_errors, step_force_scales, strict=True)
            ),
            "native_minus_rebuilt_plain_kinetic_max_abs_J": _maxabs(
                ledger[compare_indices, 23] - (frame_energy[:, 0] - np.asarray([q["kinetic_J"] for q in frame_quantities.values()]))
            ),
        }
        sampled_potential_scale = max(_maxabs(frame_energy[:, 1]), fracture_energy_scale, 1e-30)
        state_reconstruction["bulk_T4_energy_max_relative_to_fracture_scale"] = state_reconstruction["bulk_T4_energy_max_abs_difference_J"] / fracture_energy_scale
        state_reconstruction["bulk_T4_energy_max_relative_to_native_potential"] = state_reconstruction["bulk_T4_energy_max_abs_difference_J"] / sampled_potential_scale
        state_reconstruction["deformation_metrics"] = deformation_aggregate
        # Compare the saved ledger at each observed state rather than treating
        # sparse samples as if they represented the intervening steps.
        saved_force_ledger_error = _maxabs(frame_ledger[:, 14] - np.asarray(step_force_errors))
        state_reconstruction["ledger_force_residual_max_abs_difference_N"] = saved_force_ledger_error

    strict_node_checks: dict[str, Any]
    if all_node_steps_observed:
        node_state = history["node_state"]
        pos_residuals: list[float] = []
        vel_residuals: list[float] = []
        alg_increment_residuals: list[float] = []
        internal_work_increments: list[float] = []
        for step in range(1, frame_count):
            previous = np.asarray(node_state[step - 1], dtype=np.float64)
            current = np.asarray(node_state[step], dtype=np.float64)
            u_prev, v_prev, a_prev, f_prev = (previous[:, i, :] for i in range(4))
            u_cur, v_cur, a_cur, f_cur = (current[:, i, :] for i in range(4))
            pos_res, vel_res = central_difference_residuals(
                u_prev, u_cur, v_prev, v_cur, a_prev, a_cur, dt
            )
            pos_residuals.append(_maxabs(pos_res))
            vel_residuals.append(_maxabs(vel_res))
            f_ext_prev = unit_load * load_factor_history[step - 1]
            f_ext_cur = unit_load * load_factor_history[step]
            du = u_cur - u_prev
            predicted_dkalg = 0.5 * float(np.sum((f_prev + f_cur + f_ext_prev + f_ext_cur) * du))
            actual_dkalg = frame_quantities[step]["algorithmic_kinetic_J"] - frame_quantities[step - 1]["algorithmic_kinetic_J"]
            alg_increment_residuals.append(actual_dkalg - predicted_dkalg)
            internal_work_increments.append(0.5 * float(np.sum((f_prev + f_cur) * du)))
        central_position_scale = max(u_scale, dt * v_scale, 0.5 * dt**2 * a_scale, 1e-30)
        central_velocity_scale = max(v_scale, dt * a_scale, 1e-30)
        alg_residual = np.asarray(alg_increment_residuals)
        # The saved Akantu internal-force vector is a restoring force, so its
        # signed work is negative as material energy accumulates.  Preserve the
        # force-work sign and check W_internal + Delta(U + D + reversible + Ec) = 0.
        internal_work = np.r_[0.0, np.cumsum(internal_work_increments)]
        alg_k = np.asarray([frame_quantities[step]["algorithmic_kinetic_J"] for step in range(frame_count)])
        bulk_u = np.asarray(step_bulk_energy)
        internal_energy_change = (
            (bulk_u - bulk_u[0])
            + (energy[:, 2] - energy[0, 2])
            + (energy[:, 3] - energy[0, 3])
            + (energy[:, 4] - energy[0, 4])
        )
        internal_work_residual = internal_work + internal_energy_change
        total_balance_residual = reconstructed_external_work - (
            (alg_k - alg_k[0])
            + (bulk_u - bulk_u[0])
            + (energy[:, 2] - energy[0, 2])
            + (energy[:, 3] - energy[0, 3])
            + (energy[:, 4] - energy[0, 4])
        )
        actual_total_energy = alg_k + bulk_u + energy[:, 2] + energy[:, 3] + energy[:, 4]
        total_scale = max(_maxabs(actual_total_energy), fracture_energy_scale, 1e-30)
        alg_increment_scale = max(fracture_energy_scale, _maxabs(alg_k), 1e-30)
        force_relative = max(error / scale for error, scale in zip(step_force_errors, step_force_scales, strict=True))
        ledger_comparison = {
            "position_recurrence_max_abs_m": max(pos_residuals),
            "position_recurrence_relative_to_state_scale": max(pos_residuals) / central_position_scale,
            "velocity_recurrence_max_abs_mps": max(vel_residuals),
            "velocity_recurrence_relative_to_state_scale": max(vel_residuals) / central_velocity_scale,
            "Ma_minus_internal_minus_external_max_relative_to_force_scale": force_relative,
            "algorithmic_kinetic_increment_identity_max_abs_J": _maxabs(alg_residual),
            "algorithmic_kinetic_increment_identity_relative_to_fracture_energy": _maxabs(alg_residual) / fracture_energy_scale,
            "algorithmic_kinetic_increment_identity_relative_to_total_energy": _maxabs(alg_residual) / total_scale,
            "internal_force_work_plus_delta_U_D_reversible_contact_max_abs_J": _maxabs(internal_work_residual),
            "internal_force_work_balance_relative_to_fracture_energy": _maxabs(internal_work_residual) / fracture_energy_scale,
            "internal_force_work_balance_relative_to_total_energy": _maxabs(internal_work_residual) / total_scale,
            "external_work_minus_delta_Kalg_U_D_reversible_contact_max_abs_J": _maxabs(total_balance_residual),
            "external_work_balance_relative_to_fracture_energy": _maxabs(total_balance_residual) / fracture_energy_scale,
            "external_work_balance_relative_to_total_energy": _maxabs(total_balance_residual) / total_scale,
            "native_ledger_external_energy_residual_max_abs_J": _maxabs(ledger[:, 17]),
            "native_ledger_internal_energy_residual_max_abs_J": _maxabs(ledger[:, 22]),
            "native_ledger_kalg_increment_identity_residual_max_abs_J": _maxabs(ledger[:, 21]),
        }
        strict_node_checks = {
            "status": "OBSERVED",
            "all_node_states_each_solver_step": True,
            **ledger_comparison,
        }
        energy_balance = {
            "status": "OBSERVED",
            "reconstructed_external_work_J": float(reconstructed_external_work[-1]),
            "internal_work_from_saved_nodal_forces_J": float(internal_work[-1]),
            "material_energy_increase_delta_U_D_reversible_contact_J": float(internal_energy_change[-1]),
            "bulk_T4_small_strain_energy_final_J": float(bulk_u[-1]),
            "native_contact_energy_increase_J": float(energy[-1, 4] - energy[0, 4]),
            "fracture_energy_scale_Gc_times_discrete_damaged_area_J": fracture_energy_scale,
            "total_energy_scale_J": total_scale,
            "internal_force_work_balance_relative_to_fracture_energy": ledger_comparison["internal_force_work_balance_relative_to_fracture_energy"],
            "internal_force_work_balance_relative_to_total_energy": ledger_comparison["internal_force_work_balance_relative_to_total_energy"],
            "external_work_balance_relative_to_fracture_energy": ledger_comparison["external_work_balance_relative_to_fracture_energy"],
            "external_work_balance_relative_to_total_energy": ledger_comparison["external_work_balance_relative_to_total_energy"],
            "large_kinetic_energy_does_not_set_the_fracture_normalized_acceptance_scale": True,
        }
    else:
        strict_node_checks = {
            "status": "NOT_OBSERVED",
            "all_node_states_each_solver_step": False,
            "saved_node_frames": int(len(node_saved_steps)),
            "solver_steps": steps,
            "reason": "only sampled nodal states were saved; full-step central-difference and nodal-work identities are not inferred",
        }

    # The independently rebuilt per-frame state metrics are used to check the
    # original case gates without replacing or relaxing them.
    gates = case["gates"]
    deformation_gates = {
        "max_abs_principal_small_strain": {
            "observed": float(deformation_aggregate["max_abs_principal_small_strain"]),
            "limit": float(gates["max_abs_principal_small_strain"]),
            "pass": deformation_aggregate["max_abs_principal_small_strain"] <= gates["max_abs_principal_small_strain"],
        },
        "max_polar_rotation_rad": {
            "observed": float(deformation_aggregate["max_polar_rotation_rad"]),
            "limit": float(gates["max_polar_rotation_rad"]),
            "pass": deformation_aggregate["max_polar_rotation_rad"] <= gates["max_polar_rotation_rad"],
        },
        "min_det_F": {
            "observed": float(deformation_aggregate["min_det_F"]),
            "limit": float(gates["min_det_F"]),
            "pass": deformation_aggregate["min_det_F"] >= gates["min_det_F"],
        },
    }
    free_start = int(np.searchsorted(stored_times, float(loading_case["pulse_duration_s"]), side="left"))
    free_phase_report: dict[str, Any] = {"status": "NOT_OBSERVED"}
    if all_node_steps_observed:
        momentum_series = np.asarray(step_momentum)
        angular_series = np.asarray(step_angular)
        ext_torque_series = np.asarray(step_external_torque)
        int_torque_series = np.asarray(step_internal_torque)
        if free_start < frame_count:
            momentum_drift = momentum_series[free_start:] - momentum_series[free_start]
            angular_drift = angular_series[free_start:] - angular_series[free_start]
            free_phase_report = {
                "status": "OBSERVED",
                "first_free_step": free_start,
                "first_free_time_s": float(stored_times[free_start]),
                "free_phase_external_force_max_N": _maxabs(
                    unit_load[None, :, :] * load_factor_history[free_start:, None, None]
                ),
                "linear_momentum_drift_max_kg_m_s": _maxabs(momentum_drift),
                "angular_momentum_drift_max_kg_m2_s": _maxabs(angular_drift),
                "angular_momentum_is_reported_as_observed_not_required_exact_under_small_strain_formulation": True,
                "external_torque_max_norm_Nm": float(np.linalg.norm(ext_torque_series, axis=1).max()),
                "internal_torque_max_norm_Nm": float(np.linalg.norm(int_torque_series, axis=1).max()),
                "external_torque_at_final_step_Nm": ext_torque_series[-1].tolist(),
                "internal_torque_at_final_step_Nm": int_torque_series[-1].tolist(),
            }
    else:
        free_phase_report = {
            "status": "NOT_OBSERVED",
            "reason": "free-phase momentum and current-torque histories require every nodal state",
        }

    observed_checks: dict[str, bool | None] = {
        "receipt_all_sha256_and_byte_counts_match": True,
        "frozen_input_bytes_and_run_start_inventory_unchanged": True,
        "native_initial_mesh_maps_exactly_to_frozen_MSH": bool(node_lineage_exact and tet_lineage_exact and exterior_lineage_exact),
        "T4_rebuilt_mass_matches_saved_native_mass": bool(mass_matches_t4),
        "cohesive_pairs_and_quadrature_area_match_frozen_interface": bool(cohesive_lineage_exact and qp_area_match),
        "loading_weights_and_profile_reconstructed": bool(unit_load_matches and times_match and load_factor_match and max(loading_line_errors, default=0.0) <= 1e-12),
        "all_native_QP_states_complete_and_finite": True,
        "damage_irreversible_and_bounded": bool(damage_irreversible and damage_in_bounds),
        "delta_max_irreversible": bool(delta_irreversible),
        "first_full_damage_step_matches_every_QP_history": bool(first_full_match),
        "all_QPs_reached_full_damage": bool(all_qp_full),
        "native_dissipation_nonnegative_and_monotone": bool(d_nonnegative and d_monotone),
        "native_dissipation_matches_Gc_times_actual_discrete_damaged_area": bool(d_area_relative_error <= tolerance["dissipation_vs_Gc_times_discrete_area_relative"]),
        "all_QP_trapezoidal_traction_opening_work_matches_native_D_plus_reversible": bool(qp_work_relative_fracture <= tolerance["cohesive_work_vs_fracture_energy_relative"]),
        "contact_traction_and_contact_energy_are_zero": bool(contact_zero)
        if not compression_contact_enabled else None,
        "qp_failure_requires_full_damage_and_zero_total_traction_then_stays_zero": bool(
            qp_failure["all_qps_reached_full_damage_and_zero_total_traction"]
            and qp_failure["post_failure_total_traction_stays_zero"]
        ),
        "external_work_rebuilt_from_patch_displacement_and_frozen_weights": bool(external_work_relative_error <= tolerance["external_work_vs_ledger_relative"]),
        "all_native_DOFs_free": bool(all_dofs_free),
        "strict_all_step_node_checks": bool(all_node_steps_observed) if all_node_steps_observed else None,
        **compression_contact_checks,
    }
    if all_node_steps_observed:
        observed_checks.update({
            "saved_nodal_K_P_L_and_torques_reconstructed": state_reconstruction["status"] == "PASS_ALL_STEPS",
            "Ma_equals_internal_plus_external_on_all_free_DOFs": state_reconstruction["force_balance_max_relative_to_force_scale"] <= tolerance["state_and_force_recurrence_relative"],
            "central_difference_position_and_velocity_recurrences": (
                strict_node_checks["position_recurrence_relative_to_state_scale"] <= tolerance["state_and_force_recurrence_relative"]
                and strict_node_checks["velocity_recurrence_relative_to_state_scale"] <= tolerance["state_and_force_recurrence_relative"]
            ),
            "algorithmic_kinetic_increment_identity": strict_node_checks["algorithmic_kinetic_increment_identity_relative_to_fracture_energy"] <= tolerance["energy_identity_relative_to_fracture_energy"],
            "internal_force_work_balances_delta_bulk_U_plus_D_plus_reversible": strict_node_checks["internal_force_work_balance_relative_to_fracture_energy"] <= tolerance["energy_identity_relative_to_fracture_energy"],
            "external_work_matches_delta_Kalg_plus_U_plus_D_plus_reversible": strict_node_checks["external_work_balance_relative_to_fracture_energy"] <= tolerance["energy_identity_relative_to_fracture_energy"],
            "independent_bulk_T4_small_strain_energy_matches_native_potential": state_reconstruction["bulk_T4_energy_max_relative_to_native_potential"] <= tolerance["bulk_T4_energy_relative_to_native_potential"],
            "source_deformation_applicability_gates": all(value["pass"] for value in deformation_gates.values()),
        })
    else:
        for key in (
            "saved_nodal_K_P_L_and_torques_reconstructed",
            "Ma_equals_internal_plus_external_on_all_free_DOFs",
            "central_difference_position_and_velocity_recurrences",
            "algorithmic_kinetic_increment_identity",
            "internal_force_work_balances_delta_bulk_U_plus_D_plus_reversible",
            "external_work_matches_delta_Kalg_plus_U_plus_D_plus_reversible",
            "independent_bulk_T4_small_strain_energy_matches_native_potential",
            "source_deformation_applicability_gates",
        ):
            observed_checks[key] = None
    failures = [name for name, passed in observed_checks.items() if passed is False]
    return {
        "schema": "pgfrag.independent_native_vessel_audit.v1",
        "status": "AUDIT_ARITHMETIC_FAIL" if failures else "AUDIT_ARITHMETIC_PASS",
        "qualification_status": "NOT_PHYSICAL_QUALIFIED",
        "qualification_reason": "preinserted weak seam and uncalibrated prototype material; numerical identities do not establish physical qualification",
        "coverage_status": "ALL_NATIVE_STEPS_AND_NODES" if all_node_steps_observed else "PARTIAL_SAMPLED_NODE_STATES",
        "run": {
            "path": str(run.resolve()),
            "frozen_path": str(frozen.resolve()),
            "solver_steps": steps,
            "dt_s": dt,
            "total_time_s": float(stored_times[-1]),
            "native_node_count": node_count,
            "bulk_T4_count": len(tets),
            "cohesive_Tri3_count": len(cohesive),
            "cohesive_QP_count": qp_count,
            "node_saved_frames": int(len(node_saved_steps)),
            "node_stride": node_stride,
            "all_nodes_every_solver_step_saved": bool(all_node_steps_observed),
            "runner_claimed_result_status": result.get("status"),
        },
        "integrity": integrity,
        "receipt": receipt,
        "mesh_and_mass": {
            "frozen_MSH_node_count": len(source_nodes),
            "frozen_MSH_T4_count": len(source_tets),
            "frozen_interface_Tri3_count": len(source_interface),
            "native_to_source_coordinates_exact": bool(node_lineage_exact),
            "native_T4_rows_and_local_ordinals_exact_after_ancestry_map": bool(tet_lineage_exact),
            "native_cohesive_pairs_match_all_interface_faces": bool(cohesive_lineage_exact),
            "T4_volume_m3": float(volumes.sum()),
            "T4_lumped_mass_kg": float(rebuilt_mass.sum()),
            "saved_native_mass_relative_max_error": mass_rel_error,
            "all_native_DOFs_free": bool(all_dofs_free),
            "cohesive_face_area_m2": float(native_triangle_area_a.sum()),
            "quadrature_area_sum_m2": float(qp_area.sum()),
            "QP_area_weights_match_face_area_divided_by_three": bool(qp_area_match),
        },
        "loading": {
            "unit_load_matches_saved_initial_input": bool(unit_load_matches),
            "reference_resultant_N_per_N": reference_resultant.tolist(),
            "reference_torque_Nm_per_N": reference_torque.tolist(),
            "maximum_force_line_yz_reconstruction_error_m": max(loading_line_errors, default=0.0),
            "load_factor_matches_frozen_profile": bool(load_factor_match),
            "external_work_rebuilt_J": float(reconstructed_external_work[-1]),
            "native_ledger_external_work_J": float(ledger[-1, 3]),
            "external_work_max_abs_residual_J": _maxabs(external_work_residual),
            "external_work_max_relative_residual": external_work_relative_error,
        },
        "cohesive_history": {
            "all_QP_arrays_complete_and_finite": True,
            "damage_min": float(damage.min()),
            "damage_max": float(damage.max()),
            "damage_minimum_time_increment": float(damage_increment.min()),
            "damage_irreversible": damage_irreversible,
            "delta_max_min_m": float(delta_max.min()),
            "delta_max_max_m": float(delta_max.max()),
            "delta_max_minimum_increment_m": float(delta_increment.min()),
            "delta_max_irreversible": delta_irreversible,
            "first_full_damage_step_min": int(first_full[first_full >= 0].min()) if np.any(first_full >= 0) else None,
            "first_full_damage_step_max": int(first_full.max()) if np.any(first_full >= 0) else None,
            "first_full_damage_steps_match_saved_per_QP": first_full_match,
            "all_QPs_full_at_final_step": all_qp_full,
            "actual_discrete_damaged_area_final_m2": float(discrete_damaged_area[-1]),
            "Gc_times_actual_discrete_damaged_area_final_J": fracture_energy_scale,
            "native_dissipated_energy_final_J": float(native_dissipated[-1]),
            "dissipation_vs_Gc_Ah_max_abs_error_J": _maxabs(d_area_residual),
            "dissipation_vs_Gc_Ah_max_relative_to_fracture_energy": d_area_relative_error,
            "dissipation_vs_Gc_Ah_final_relative_to_fracture_energy": final_d_area_relative_error,
            "native_dissipated_energy_min_J": float(native_dissipated.min()),
            "native_dissipated_energy_monotone": d_monotone,
            "contact_traction_max_abs_Pa": contact_max,
            "cohesive_contact_energy_max_abs_J": contact_energy_max,
            "qp_trapezoidal_work_final_J": float(qp_work[-1]),
            "native_D_plus_reversible_final_J": float(native_cohesive_energy[-1]),
            "qp_work_max_abs_residual_J": qp_work_error,
            "qp_work_max_relative_to_fracture_energy": qp_work_relative_fracture,
            "qp_work_max_relative_to_total_energy": qp_work_relative_total,
            "fracture_energy_normalization_J": fracture_energy_scale,
            "total_energy_normalization_J": total_energy_scale,
        },
        "qp_failure": qp_failure,
        "compression_contact": contact_work_report,
        "contact_stiffness_scope": contact_stiffness_report,
        "saved_node_state": state_reconstruction,
        "strict_all_step_dynamics": strict_node_checks,
        "energy_balance": energy_balance,
        "free_phase_momentum_and_current_torque": free_phase_report,
        "source_deformation_applicability_gates": deformation_gates,
        "audit_tolerances": tolerance,
        "checks": observed_checks,
        "failed_checks": failures,
        "limitations": [
            "This audit checks serialized arithmetic and the discrete linear-elastic model only; a pass is not material calibration or physical qualification.",
            "Angular momentum and current internal/external torques are reported from saved states. No all-time exact angular-momentum claim is made for the small-strain formulation.",
            "A run with sparse node_saved_steps reports strict all-step nodal identities as NOT_OBSERVED, even when sampled states and all-step QP histories are available.",
            "When compression penalty is enabled, contact opening magnitude is reported from the native contact-opening vector; paired-face gap reconstruction is limited to saved nodal frames.",
            "The compression penalty check covers only the fixed-normal tangent envelope when penalty <= Kn and beta=kappa=1; it is not a nonlinear or full-trajectory stability bound.",
        ],
    }


def audit_run(run_dir: Path, frozen_dir: Path, output_dir: Path) -> dict[str, Any]:
    run = Path(run_dir).resolve(strict=True)
    frozen = Path(frozen_dir).resolve(strict=True)
    output = Path(output_dir).resolve()
    if run == output or run in output.parents or frozen == output or frozen in output.parents:
        raise ValueError("audit output must be outside run and frozen input directories")
    if output.exists():
        raise FileExistsError(f"audit output directory already exists: {output}")
    output.mkdir(parents=True, exist_ok=False)
    base: dict[str, Any] = {
        "schema": "pgfrag.independent_native_vessel_audit.v1",
        "qualification_status": "NOT_PHYSICAL_QUALIFIED",
        "run_path": str(run),
        "frozen_path": str(frozen),
        "output_path": str(output),
    }
    try:
        receipt = _verify_receipt(run)
        integrity = _verify_frozen_inputs(run, frozen)
    except AuditInputError as error:
        base.update({
            "status": "INPUT_ERROR",
            "coverage_status": "NOT_OBSERVED",
            "input_error": str(error),
            "limitations": ["Input integrity or readability failed; no numerical qualification was attempted."],
        })
        (output / "audit.json").write_text(json.dumps(base, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        return base
    base["receipt"] = receipt
    base["integrity"] = integrity
    try:
        numerical = _audit_native_run(run, frozen, integrity, receipt)
    except AuditInputError as error:
        base.update({
            "status": "INPUT_ERROR",
            "coverage_status": "NOT_OBSERVED",
            "input_error": str(error),
            "limitations": ["Receipt and frozen source integrity were verified, but a serialized numerical input was malformed or incomplete."],
        })
    except AuditArithmeticError as error:
        base.update({
            "status": "AUDIT_ARITHMETIC_FAIL",
            "coverage_status": "NUMERICAL_INPUTS_READ",
            "audit_error": str(error),
            "limitations": ["Numerical histories were readable, but an independent arithmetic precondition failed."],
        })
    else:
        base = numerical
        base["output_path"] = str(output)
    (output / "audit.json").write_text(json.dumps(base, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return base


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Read-only independent audit of a native cohesive-vessel run")
    parser.add_argument("--run", type=Path, required=True, help="completed native run output directory")
    parser.add_argument("--frozen", type=Path, required=True, help="immutable frozen input bundle directory")
    parser.add_argument("--out", type=Path, required=True, help="new directory for audit.json")
    arguments = parser.parse_args(argv)
    try:
        report = audit_run(arguments.run, arguments.frozen, arguments.out)
    except Exception as error:
        print(f"audit failed before report creation: {error}", file=sys.stderr)
        return 2
    print(json.dumps({
        "status": report.get("status"),
        "qualification_status": report.get("qualification_status"),
        "audit_json": str(arguments.out.resolve() / "audit.json"),
        "input_error": report.get("input_error"),
        "failed_checks": report.get("failed_checks"),
    }, allow_nan=False))
    return 0 if report.get("status") in ("INPUT_ERROR", "AUDIT_ARITHMETIC_PASS") else 1


if __name__ == "__main__":
    raise SystemExit(main())

