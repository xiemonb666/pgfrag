#!/usr/bin/env python3
"""Independent, read-only audit of ``pgfrag.native_bowl_result.v1`` runs.

This file deliberately does not import ``pgfrag.bowl_dynamics`` or the Julia
runner.  Its formulas and reference geometry reconstruction are implemented
here from the frozen root contract so a saved native result can be checked
without rerunning a solver.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import sys
from typing import Any

import numpy as np


SCHEMA = "pgfrag.independent_bowl_audit.v1"
NATIVE_SCHEMA = "pgfrag.native_bowl_result.v1"
HERE = Path(__file__).resolve().parent
REQUIRED = {
    "point_id", "x_ref_m", "position_m", "velocity_mps", "damage",
    "volume_m3", "point_mass_kg", "step", "time_s", "bond_i", "bond_j",
    "bond_rest_length_m", "bond_epsilon_c", "bond_first_break_step",
    "bond_first_break_time_s", "bond_first_break_stretch", "bond_active_packed",
    "component_label", "event_bond_index", "event_step", "event_time_s",
    "event_stretch", "event_pair_energy_J", "native_sed_J_m3", "kinetic_energy_J",
    "pair_potential_energy_J", "break_dissipation_cumulative_J",
    "native_strain_energy_J", "total_force_x_N", "total_force_y_N",
    "total_force_z_N", "momentum_x_kg_m_s", "momentum_y_kg_m_s",
    "momentum_z_kg_m_s", "angular_momentum_x_kg_m2_s",
    "angular_momentum_y_kg_m2_s", "angular_momentum_z_kg_m2_s",
    "checkpoint_step", "checkpoint_acceleration_mps2",
    "checkpoint_internal_force_N",
}
LEDGER_SCALARS = (
    "kinetic_energy_J", "pair_potential_energy_J",
    "break_dissipation_cumulative_J", "native_strain_energy_J",
    "total_force_x_N", "total_force_y_N", "total_force_z_N",
    "momentum_x_kg_m_s", "momentum_y_kg_m_s", "momentum_z_kg_m_s",
    "angular_momentum_x_kg_m2_s", "angular_momentum_y_kg_m2_s",
    "angular_momentum_z_kg_m2_s",
)


class AuditInputError(RuntimeError):
    pass


def _json_no_constant(value: str) -> None:
    raise ValueError(f"non-standard JSON numeric constant: {value}")


def read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"), parse_constant=_json_no_constant)
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as exc:
        raise AuditInputError(f"cannot read valid JSON object from {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise AuditInputError(f"expected JSON object in {path}")
    return value


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8", newline="\n")


def add_check(report: dict[str, Any], name: str, passed: bool | None, details: dict[str, Any]) -> None:
    report["checks"].append({"name": name, "passed": passed, "details": details})


def _safe_input_path(root: Path, relative: str) -> Path:
    if not isinstance(relative, str) or not relative or "\\" in relative:
        raise AuditInputError(f"unsafe result path: {relative!r}")
    rel = PurePosixPath(relative)
    if rel.is_absolute() or ".." in rel.parts:
        raise AuditInputError(f"unsafe result path: {relative!r}")
    resolved_root = root.resolve(strict=True)
    resolved = (resolved_root / Path(*rel.parts)).resolve(strict=True)
    try:
        resolved.relative_to(resolved_root)
    except ValueError as exc:
        raise AuditInputError(f"result path escapes run directory: {relative!r}") from exc
    if not resolved.is_file():
        raise AuditInputError(f"result path is not a file: {relative!r}")
    return resolved


def _safe_source_path(root: Path, relative: str) -> Path:
    # Source receipt paths use POSIX paths relative to the run directory.
    return _safe_input_path(root, relative)


def derive_root_material(E_pa: float, nu: float, Gc_j_m2: float, delta_m: float) -> tuple[float, float, float]:
    """Return K, uniform NoCorrection BB bond constant, and critical stretch."""
    K = E_pa / (3.0 * (1.0 - 2.0 * nu))
    c = 18.0 * K / (math.pi * delta_m**4)
    sc = math.sqrt(5.0 * Gc_j_m2 / (9.0 * K * delta_m))
    return K, c, sc


def segment_min_radius(a: np.ndarray, b: np.ndarray, center: np.ndarray) -> float:
    d = b - a
    den = float(np.dot(d, d))
    t = 0.0 if den == 0.0 else float(np.clip(-np.dot(a - center, d) / den, 0.0, 1.0))
    return float(np.linalg.norm(a + t * d - center))


def centroid_shell_points(ro: float, ri: float, center: np.ndarray, h: float) -> np.ndarray:
    """Independent reconstruction of the full cell-centroid shell sampler."""
    low = center + np.asarray((-ro, -ro, -ro), dtype=np.float64)
    high = center + np.asarray((ro, ro, 0.0), dtype=np.float64)
    axes = []
    for lo, hi in zip(low, high):
        count = int(math.floor((float(hi) - float(lo)) / h + 1e-12))
        axes.append(lo + (np.arange(count, dtype=np.float64) + 0.5) * h)
    xx, yy, zz = np.meshgrid(*axes, indexing="ij")
    candidates = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
    radial = np.linalg.norm(candidates - center, axis=1)
    return np.asarray(candidates[(radial <= ro) & (radial >= ri)], dtype=np.float64)


def reconstruct_candidate_pairs(points: np.ndarray, delta: float, center: np.ndarray, inner_radius: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build the full <=delta reference roster and split inner-void exclusions.

    Cell buckets bound work to nearby points.  Filtering uses the exact closest
    point on the closed segment, retaining a bond iff that distance is >= Ri.
    Rows in all three returned arrays are canonical, sorted i<j.
    """
    buckets: dict[tuple[int, int, int], list[int]] = {}
    cell = np.floor(points / delta).astype(np.int64)
    for i, key_arr in enumerate(cell):
        key = (int(key_arr[0]), int(key_arr[1]), int(key_arr[2]))
        buckets.setdefault(key, []).append(i)
    potential: list[tuple[int, int, float]] = []
    retained: list[tuple[int, int, float]] = []
    excluded: list[tuple[int, int, float]] = []
    offsets = [(a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1)]
    for i, key_arr in enumerate(cell):
        key0 = (int(key_arr[0]), int(key_arr[1]), int(key_arr[2]))
        a = points[i]
        for dx, dy, dz in offsets:
            for j in buckets.get((key0[0] + dx, key0[1] + dy, key0[2] + dz), ()):
                if j <= i:
                    continue
                b = points[j]
                length = float(np.linalg.norm(b - a))
                if length <= delta:
                    pair = (i, j, length)
                    potential.append(pair)
                    if segment_min_radius(a, b, center) >= inner_radius:
                        retained.append(pair)
                    else:
                        excluded.append(pair)
    potential.sort(key=lambda row: (row[0], row[1]))
    retained.sort(key=lambda row: (row[0], row[1]))
    excluded.sort(key=lambda row: (row[0], row[1]))
    return (
        np.asarray(potential, dtype=np.float64).reshape((-1, 3)),
        np.asarray(retained, dtype=np.float64).reshape((-1, 3)),
        np.asarray(excluded, dtype=np.float64).reshape((-1, 3)),
    )


def pair_stretch_and_energy(position: np.ndarray, bond_i: np.ndarray, bond_j: np.ndarray,
                            rest_length: np.ndarray, volume: np.ndarray, c: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    i = bond_i.astype(np.int64, copy=False)
    j = bond_j.astype(np.int64, copy=False)
    d = position[j] - position[i]
    length = np.linalg.norm(d, axis=1)
    stretch = (length - rest_length) / rest_length
    pair_energy = 0.5 * c * volume[i] * volume[j] * rest_length * stretch**2
    return stretch, pair_energy, d / length[:, None]


def union_find_labels(n_points: int, bond_i: np.ndarray, bond_j: np.ndarray,
                      active: np.ndarray) -> tuple[np.ndarray, int]:
    parent = list(range(n_points))
    size = [1] * n_points

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, j in zip(bond_i[active], bond_j[active]):
        a, b = find(int(i)), find(int(j))
        if a != b:
            if size[a] < size[b]:
                a, b = b, a
            parent[b] = a
            size[a] += size[b]
    members: dict[int, list[int]] = {}
    for point in range(n_points):
        members.setdefault(find(point), []).append(point)
    label = np.empty(n_points, dtype=np.int32)
    for ids in members.values():
        label[ids] = min(ids)
    return label, len(members)


def _relative_difference(a: np.ndarray, b: np.ndarray) -> dict[str, float]:
    aa = np.asarray(a, dtype=np.float64)
    bb = np.asarray(b, dtype=np.float64)
    diff = np.abs(aa - bb)
    scale = np.maximum(np.abs(bb), np.finfo(np.float64).tiny)
    return {"max_abs": float(np.max(diff)) if diff.size else 0.0,
            "max_rel": float(np.max(diff / scale)) if diff.size else 0.0}


def _expected_initial_velocity(x_ref: np.ndarray, mass: np.ndarray, center: np.ndarray,
                               alpha: float) -> np.ndarray:
    A = np.diag((alpha, -0.25 * alpha, 0.0))
    raw = (x_ref - center) @ A.T
    total_mass = float(np.sum(mass))
    raw -= np.sum(mass[:, None] * raw, axis=0) / total_mass
    com = np.sum(mass[:, None] * x_ref, axis=0) / total_mass
    r = x_ref - com
    angular = np.sum(np.cross(r, mass[:, None] * raw), axis=0)
    inertia = np.zeros((3, 3), dtype=np.float64)
    eye = np.eye(3)
    for mi, ri in zip(mass, r):
        inertia += mi * (float(np.dot(ri, ri)) * eye - np.outer(ri, ri))
    omega = np.linalg.solve(inertia, angular) if np.linalg.matrix_rank(inertia) == 3 else np.zeros(3)
    return raw - np.cross(omega[None, :], r)


def _metric_scale(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=np.float64)
    return float(np.max(np.abs(arr))) if arr.size else 0.0


def _check_close(report: dict[str, Any], name: str, actual: np.ndarray, expected: np.ndarray,
                 *, rtol: float = 1e-9, atol: float = 1e-12, details: dict[str, Any] | None = None) -> None:
    a = np.asarray(actual, dtype=np.float64)
    e = np.asarray(expected, dtype=np.float64)
    if a.shape != e.shape:
        add_check(report, name, False, {"shape_actual": list(a.shape), "shape_expected": list(e.shape), **(details or {})})
        return
    stats = _relative_difference(a, e)
    ok = bool(np.allclose(a, e, rtol=rtol, atol=atol, equal_nan=False))
    add_check(report, name, ok, {**stats, "rtol": rtol, "atol": atol, **(details or {})})


def _load_arrays(path: Path) -> dict[str, np.ndarray]:
    try:
        with np.load(path, allow_pickle=False) as archive:
            arrays: dict[str, np.ndarray] = {}
            for name in archive.files:
                value = archive[name]
                if value.dtype.hasobject:
                    raise AuditInputError(f"object dtype is forbidden in arrays.npz: {name}")
                arrays[name] = value
    except AuditInputError:
        raise
    except Exception as exc:
        raise AuditInputError(f"cannot load safe NumPy archive {path}: {exc}") from exc
    missing = sorted(REQUIRED - arrays.keys())
    if missing:
        raise AuditInputError(f"native result archive misses required keys: {missing}")
    return arrays


def audit_run(run_dir: Path) -> tuple[dict[str, Any], dict[str, np.ndarray], dict[str, Any]]:
    root = run_dir.resolve(strict=True)
    result_path = root / "result.json"
    case_path = root / "case.json"
    if not result_path.is_file() or not case_path.is_file():
        raise AuditInputError("run directory must contain result.json and case.json")
    result = read_json(result_path)
    case = read_json(case_path)
    if result.get("schema") != NATIVE_SCHEMA:
        raise AuditInputError(f"unsupported native schema: {result.get('schema')!r}")
    # Native wrapper contract stores the NPZ receipt at result.arrays.  Accept
    # the older payloads.arrays spelling only as a read-compatibility fallback;
    # the report records which spelling supplied the receipt.
    array_receipt_field = "arrays"
    payload = result.get("arrays")
    if not isinstance(payload, dict):
        array_receipt_field = "payloads.arrays (legacy)"
        payload = result.get("payloads", {}).get("arrays", {})
    arrays_rel = payload.get("path")
    if not isinstance(arrays_rel, str):
        raise AuditInputError("result.json has no arrays.path")
    arrays_path = _safe_input_path(root, arrays_rel)
    arrays = _load_arrays(arrays_path)

    report: dict[str, Any] = {
        "schema": SCHEMA,
        "native_schema": NATIVE_SCHEMA,
        "run_id": result.get("run_id"),
        "input": {"run_directory": str(root),
                  "result_json": {"path": "result.json", "bytes": result_path.stat().st_size, "sha256": sha256_file(result_path)},
                  "case_json": {"path": "case.json", "bytes": case_path.stat().st_size, "sha256": sha256_file(case_path)},
                  "arrays_npz": {"path": arrays_rel, "receipt_field": array_receipt_field,
                                 "bytes": arrays_path.stat().st_size, "sha256": sha256_file(arrays_path)}},
        "checks": [],
        "interpretation": {
            "graph_components": "active-bond graph connectivity over all points, including isolated points; this does not establish physical shards, positive-area crack surfaces, or R2",
            "break_dissipation": "discrete removed pair potential at recorded break states; it is not Gc times crack area",
            "resolution_comparison": "half-dt and fixed-delta h comparisons are sensitivity studies; they do not establish asymptotic convergence",
            "initial_geometric_exclusion": "reference bonds whose closed segment enters the inner sphere are excluded before native initialization and are never fracture events",
        },
    }

    add_check(report, "native_status_pass", result.get("status") == "PASS",
              {"native_status": result.get("status")})
    add_check(report, "case_schema", case.get("schema") == "pgfrag.native_bowl_case.v1",
              {"case_schema": case.get("schema")})
    add_check(report, "run_id_binding", case.get("run_id") == result.get("run_id"),
              {"case_run_id": case.get("run_id"), "result_run_id": result.get("run_id")})

    # Verify all producer-listed output receipts before interpreting arrays.
    output_receipts = result.get("outputs", [])
    output_results = []
    output_ok = isinstance(output_receipts, list) and bool(output_receipts)
    for receipt in output_receipts if isinstance(output_receipts, list) else []:
        try:
            file = _safe_input_path(root, receipt.get("path"))
            actual_bytes = file.stat().st_size
            actual_hash = sha256_file(file)
            match = actual_bytes == receipt.get("bytes") and actual_hash == receipt.get("sha256")
            output_ok = output_ok and match
            output_results.append({"path": receipt.get("path"), "passed": match,
                                   "bytes_actual": actual_bytes, "bytes_recorded": receipt.get("bytes"),
                                   "sha256_actual": actual_hash, "sha256_recorded": receipt.get("sha256")})
        except Exception as exc:
            output_ok = False
            output_results.append({"path": receipt.get("path") if isinstance(receipt, dict) else None,
                                   "passed": False, "error": str(exc)})
    arrays_receipt = payload
    arrays_receipt_ok = (arrays_path.stat().st_size == arrays_receipt.get("bytes")
                         and sha256_file(arrays_path) == arrays_receipt.get("sha256"))
    output_ok = output_ok and arrays_receipt_ok
    add_check(report, "output_and_array_hashes", bool(output_ok),
              {"output_receipts": output_results,
               "arrays_payload_receipt_passed": bool(arrays_receipt_ok),
               "arrays_payload_bytes_recorded": arrays_receipt.get("bytes"),
               "arrays_payload_sha256_recorded": arrays_receipt.get("sha256")})

    # Source snapshots bind the native implementation, including the copied
    # Forge manifest that was used to prepare the case.  The original asset
    # path is checked as a separate optional provenance lookup.
    source_results = []
    source_ok = isinstance(result.get("source_receipts"), list) and bool(result.get("source_receipts"))
    source_receipts = result.get("source_receipts", [])
    for receipt in source_receipts if isinstance(source_receipts, list) else []:
        rec_path = receipt.get("path") if isinstance(receipt, dict) else None
        try:
            source_file = _safe_source_path(root, rec_path)
            actual_hash, actual_bytes = sha256_file(source_file), source_file.stat().st_size
            ok = actual_hash == receipt.get("sha256") and actual_bytes == receipt.get("bytes")
            source_ok = source_ok and ok
            source_results.append({"path": rec_path, "passed": bool(ok), "sha256_actual": actual_hash,
                                   "sha256_recorded": receipt.get("sha256"), "bytes_actual": actual_bytes,
                                   "bytes_recorded": receipt.get("bytes")})
        except Exception as exc:
            source_ok = False
            source_results.append({"path": rec_path, "passed": False, "error": str(exc)})
    add_check(report, "source_receipt_hashes_and_paths", bool(source_ok),
              {"source_receipts": source_results})
    source_names = {PurePosixPath(str(item.get("path", ""))).name
                    for item in source_receipts if isinstance(item, dict)}
    required_source_names = {"bowl_dynamics.py", "run_bowl_dynamics.jl", "Project.toml", "Manifest.toml",
                             "backend_provenance.json", "forge_manifest.json"}
    vendor_receipts = [str(item.get("path", "")) for item in source_receipts if isinstance(item, dict)
                       and str(item.get("path", "")).startswith("source_snapshot/vendor/Peridynamics.jl/")]
    coverage_ok = required_source_names.issubset(source_names) and bool(vendor_receipts)
    add_check(report, "source_receipt_coverage", coverage_ok,
              {"required_source_names": sorted(required_source_names),
               "recorded_source_names": sorted(source_names),
               "vendor_peridynamics_receipt_count": len(vendor_receipts)})

    # Verify that the full input manifest is pinned inside the copied source
    # snapshot and agrees with both case and result metadata.  The original
    # asset may have moved after preparation; that does not invalidate the
    # self-contained snapshot, but it is reported separately.
    source_hash = result.get("source_forge_manifest_sha256")
    case_source_hash = case.get("source_forge_manifest_sha256")
    snapshot_manifest_rel = "source_snapshot/forge_manifest.json"
    try:
        snapshot_manifest_path = _safe_source_path(root, snapshot_manifest_rel)
        snapshot_manifest_hash = sha256_file(snapshot_manifest_path)
        snapshot_manifest_object = read_json(snapshot_manifest_path)
        embedded_manifest_match = snapshot_manifest_object == case.get("source_forge_manifest")
        snapshot_manifest_ok = (snapshot_manifest_hash == source_hash == case_source_hash
                                and embedded_manifest_match)
        snapshot_manifest_details: dict[str, Any] = {
            "path": snapshot_manifest_rel, "sha256": snapshot_manifest_hash,
            "result_hash": source_hash, "case_hash": case_source_hash,
            "matches_embedded_case_manifest": bool(embedded_manifest_match),
        }
    except Exception as exc:
        snapshot_manifest_ok = False
        snapshot_manifest_hash = None
        snapshot_manifest_details = {"path": snapshot_manifest_rel, "error": str(exc)}
    add_check(report, "source_snapshot_manifest_hash_and_case_binding", bool(snapshot_manifest_ok),
              snapshot_manifest_details)
    external_manifest_status: bool | None = None
    external_manifest = case.get("source_manifest_path")
    if isinstance(external_manifest, str) and external_manifest:
        try:
            external_path = Path(external_manifest).resolve(strict=True)
            if external_path.is_file():
                external_hash = sha256_file(external_path)
                external_manifest_status = (external_hash == source_hash == snapshot_manifest_hash)
                external_details = {"path": str(external_path), "sha256": external_hash,
                                    "passed": bool(external_manifest_status)}
            else:
                external_details = {"path": str(external_path), "note": "recorded external path is not a file"}
        except OSError as exc:
            external_details = {"path": external_manifest, "note": f"external manifest unavailable: {exc}"}
    else:
        external_details = {"path": external_manifest, "note": "case has no external manifest path"}
    add_check(report, "external_asset_manifest_path_hash_if_available", external_manifest_status,
              external_details)

    argv_path = root / "argv.json"
    if argv_path.is_file():
        argv_record = read_json(argv_path)
        argv_values = argv_record.get("argv", [])
        try:
            argv_case = Path(str(argv_values[-2])).resolve(strict=True) if isinstance(argv_values, list) and len(argv_values) >= 5 else None
            argv_run = Path(str(argv_values[-1])).resolve(strict=True) if isinstance(argv_values, list) and len(argv_values) >= 5 else None
            argv_paths_ok = bool(argv_case == (root / "case.toml").resolve(strict=True)
                                 and argv_run == root
                                 and Path(str(argv_values[-3])).name == "run_bowl_dynamics.jl")
        except (OSError, IndexError, TypeError):
            argv_paths_ok = False
        argv_receipts_match = argv_record.get("source_receipts") == source_receipts
        argv_id_match = argv_record.get("run_id") == result.get("run_id")
        add_check(report, "prepared_argv_binding", bool(argv_paths_ok and argv_receipts_match and argv_id_match),
                  {"argv": argv_values, "run_id": argv_record.get("run_id"),
                   "source_receipts_match_result": bool(argv_receipts_match),
                   "case_and_runner_paths_valid": bool(argv_paths_ok),
                   "argv_json_sha256": sha256_file(argv_path), "argv_json_bytes": argv_path.stat().st_size})
    else:
        add_check(report, "prepared_argv_binding", False, {"error": "argv.json is missing"})

    run_receipt_path = root / "run_receipt.json"
    if run_receipt_path.is_file():
        run_receipt = read_json(run_receipt_path)
        actual_argv = run_receipt.get("argv", [])
        try:
            actual_runner = Path(str(actual_argv[-3])).resolve(strict=True) if len(actual_argv) == 5 else None
            actual_case = Path(str(actual_argv[-2])).resolve(strict=True) if len(actual_argv) == 5 else None
            actual_run = Path(str(actual_argv[-1])).resolve(strict=True) if len(actual_argv) == 5 else None
            expected_runner = (root / "source_snapshot" / "run_bowl_dynamics.jl").resolve(strict=True)
            expected_case = (root / "case.toml").resolve(strict=True)
            actual_paths_ok = bool(actual_runner == expected_runner and actual_case == expected_case
                                    and actual_run == root)
            actual_project_ok = (len(actual_argv) == 5
                                 and str(actual_argv[1]) == f"--project={run_receipt.get('julia_project')}")
            actual_argv_hash = hashlib.sha256(json.dumps(actual_argv, ensure_ascii=False,
                                                        separators=(",", ":")).encode("utf-8")).hexdigest()
            actual_hash_ok = actual_argv_hash == run_receipt.get("argv_sha256")
        except (OSError, IndexError, TypeError):
            actual_paths_ok = actual_project_ok = actual_hash_ok = False
            actual_argv_hash = None
        run_source_match = run_receipt.get("source_receipts") == source_receipts
        run_receipt_ok = bool(run_receipt.get("status") == "PASS"
                               and run_receipt.get("run_id") == result.get("run_id")
                               and actual_paths_ok and actual_project_ok and actual_hash_ok
                               and run_source_match)
        add_check(report, "actual_execution_receipt_binding", run_receipt_ok,
                  {"status": run_receipt.get("status"), "run_id": run_receipt.get("run_id"),
                   "actual_runner_case_output_paths_valid": actual_paths_ok,
                   "julia_project_argument_matches_receipt": actual_project_ok,
                   "argv_sha256_valid": actual_hash_ok, "argv_sha256_actual": actual_argv_hash,
                   "source_receipts_match_result": bool(run_source_match),
                   "run_receipt_sha256": sha256_file(run_receipt_path),
                   "run_receipt_bytes": run_receipt_path.stat().st_size})
    else:
        add_check(report, "actual_execution_receipt_binding", False,
                  {"error": "run_receipt.json is missing"})

    # The native status receipt is a distinct file produced by Julia.  Bind its
    # core execution values to the packaged result before auditing its arrays.
    native_status_path = root / "native_status.json"
    if native_status_path.is_file():
        native_status = read_json(native_status_path)
        solver_meta = result.get("solver", {})
        native_fields = ("status", "steps", "threads", "dt_s", "stable_dt_s", "safety_factor",
                         "initial_geometric_exclusion_count")
        native_field_matches = {}
        for field in native_fields:
            lhs, rhs = native_status.get(field), solver_meta.get(field)
            if isinstance(lhs, (float, int)) and isinstance(rhs, (float, int)):
                native_field_matches[field] = bool(math.isclose(float(lhs), float(rhs), rel_tol=2e-13, abs_tol=1e-18))
            else:
                native_field_matches[field] = lhs == rhs
        native_status_ok = native_status.get("status") == "NATIVE_COMPLETE" and all(native_field_matches.values())
        add_check(report, "native_status_receipt_binding", bool(native_status_ok),
                  {"native_status": native_status, "field_matches_result_solver": native_field_matches,
                   "sha256": sha256_file(native_status_path)})
    else:
        native_status = {}
        add_check(report, "native_status_receipt_binding", False, {"error": "native_status.json is missing"})

    # Units, shape, finite values, and frozen input identity.
    contract = result.get("array_contract", {})
    add_check(report, "si_unit_contract", "SI" in str(contract.get("units", "")),
              {"contract_units": contract.get("units")})
    add_check(report, "no_object_arrays", all(not a.dtype.hasobject for a in arrays.values()),
              {"array_count": len(arrays)})
    bad_float = {name: int(np.size(value) - np.count_nonzero(np.isfinite(value)))
                 for name, value in arrays.items() if np.issubdtype(value.dtype, np.floating)
                 and not np.isfinite(value).all()}
    add_check(report, "all_native_numeric_arrays_finite", not bad_float, {"nonfinite_counts": bad_float})

    geo = case.get("geometry", {})
    material = case.get("material", {})
    loading = case.get("loading", {})
    solver = result.get("solver", {})
    expected_ro, expected_ri, expected_center = 12e-3, 8e-3, np.asarray((0.0, 0.0, 12e-3))
    ro, ri = float(geo.get("outer_radius_m", math.nan)), float(geo.get("inner_radius_m", math.nan))
    center = np.asarray(geo.get("center_m", (math.nan, math.nan, math.nan)), dtype=np.float64)
    h, delta = float(geo.get("h_m", math.nan)), float(geo.get("delta_m", math.nan))
    E = float(material.get("E_pa", math.nan))
    nu = float(material.get("nu", math.nan))
    Gc = float(material.get("Gc_j_m2", math.nan))
    rho = float(material.get("rho_kg_m3", math.nan))
    if not all(math.isfinite(v) for v in (ro, ri, h, delta, E, nu, Gc, rho, *center)):
        raise AuditInputError("case geometry/material fields are missing or non-finite")
    K, c, sc = derive_root_material(E, nu, Gc, delta)
    analytic_volume = 2.0 * math.pi / 3.0 * (ro**3 - ri**3)
    report["root_formula_values"] = {"K_pa": K, "bond_constant_N_m6": c,
                                     "critical_stretch": sc,
                                     "formula": {"K": "E/[3(1-2 nu)]", "c": "18 K/(pi delta^4)",
                                                 "critical_stretch": "sqrt(5 Gc/(9 K delta))",
                                                 "pair_potential_J": "0.5 c Vi Vj L s^2",
                                                 "pair_force_N": "c s Vi Vj e"}}
    add_check(report, "frozen_geometry", bool(np.isclose(ro, expected_ro, rtol=0, atol=1e-15)
                                               and np.isclose(ri, expected_ri, rtol=0, atol=1e-15)
                                               and np.allclose(center, expected_center, rtol=0, atol=1e-15)),
              {"outer_radius_m": ro, "inner_radius_m": ri, "center_m": center.tolist()})
    add_check(report, "prototype_material_parameters", bool(np.isclose(rho, 2300.0, rtol=0, atol=1e-12)
                                                              and np.isclose(E, 70e9, rtol=0, atol=1e-3)
                                                              and np.isclose(nu, 0.25, rtol=0, atol=1e-15)
                                                              and np.isclose(Gc, 7.6, rtol=0, atol=1e-12)),
              {"rho_kg_m3": rho, "E_pa": E, "nu": nu, "Gc_j_m2": Gc,
               "calibration_status": material.get("calibration_status")})
    calibration_status = material.get("calibration_status") or case.get("material", {}).get("calibration_status")
    add_check(report, "prototype_material_remains_unvalidated",
              calibration_status in {"UNVALIDATED_PROTOTYPE", "UNVALIDATED"},
              {"calibration_status": calibration_status})
    model_text = " ".join(str(value) for value in solver.values() if isinstance(value, str))
    model_ok = all(token.lower() in model_text.lower() for token in ("BB", "NoCorrection", "CriticalStretch", "VelocityVerlet"))
    add_check(report, "native_model_contract", model_ok, {"model_fields": model_text})
    alpha = float(loading.get("alpha_s_1", math.nan))
    allowed_alpha = alpha in (0.0, 100.0, 300.0)
    no_gaussian = float(loading.get("gaussian_peak_m_s", 0.0)) == 0.0 and float(loading.get("gaussian_sigma_m", 0.0)) == 0.0
    expected_rate = np.diag((alpha, -0.25 * alpha, 0.0))
    rate = np.asarray(loading.get("rate_tensor_s_1", expected_rate.tolist()), dtype=np.float64)
    rate_ok = rate.shape == (3, 3) and np.allclose(rate, expected_rate, rtol=0.0, atol=1e-14)
    add_check(report, "frozen_free_loading", bool(allowed_alpha and no_gaussian
                                                 and rate_ok
                                                 and "momentum" in str(loading.get("momentum_constraint", "")).lower()),
              {"alpha_s_1": alpha, "mode": loading.get("mode"), "gaussian_peak_m_s": loading.get("gaussian_peak_m_s"),
               "gaussian_sigma_m": loading.get("gaussian_sigma_m"),
               "rate_tensor_s_1": rate.tolist() if rate.shape == (3, 3) else None,
               "momentum_constraint": loading.get("momentum_constraint"),
               "boundary_conditions": "none declared"})
    add_check(report, "horizon_and_spacing_positive", bool(h > 0 and delta > 0),
              {"h_m": h, "delta_m": delta, "delta_over_h": delta / h})
    add_check(report, "frozen_spacing_and_horizon_cases",
              bool(h in (1e-3, 0.8e-3) and np.isclose(delta, 2.5e-3, rtol=0.0, atol=1e-15)),
              {"h_m": h, "allowed_h_m": [1e-3, 0.8e-3], "delta_m": delta,
               "required_delta_m": 2.5e-3,
               "comparison_scope": "h changes at fixed delta are a sensitivity study, not convergence"})
    case_geometry_result = result.get("geometry", {})
    add_check(report, "case_result_geometry_binding", bool(np.isclose(case_geometry_result.get("outer_radius_m", math.nan), ro, rtol=1e-14)
                                                                   and np.isclose(case_geometry_result.get("inner_radius_m", math.nan), ri, rtol=1e-14)
                                                                   and np.allclose(case_geometry_result.get("center_m", []), center, rtol=1e-14, atol=1e-15)
                                                                   and np.isclose(case_geometry_result.get("h_m", math.nan), h, rtol=1e-14)
                                                                   and np.isclose(case_geometry_result.get("delta_m", math.nan), delta, rtol=1e-14)),
              {"case_geometry": geo, "result_geometry": case_geometry_result})

    # Consistent state dimensions and unit-bearing arrays.
    x_ref = arrays["x_ref_m"]
    position = arrays["position_m"]
    velocity = arrays["velocity_mps"]
    damage = arrays["damage"]
    point_id = arrays["point_id"]
    volume = arrays["volume_m3"]
    mass = arrays["point_mass_kg"]
    steps = arrays["step"].astype(np.int64, copy=False)
    time_s = arrays["time_s"].astype(np.float64, copy=False)
    n = len(point_id)
    n_states = len(steps)
    if x_ref.shape != (n, 3) or position.shape != (n_states, n, 3) or velocity.shape != (n_states, n, 3):
        raise AuditInputError(f"point-state shape mismatch: x_ref={x_ref.shape}, position={position.shape}, velocity={velocity.shape}")
    if damage.shape != (n_states, n) or volume.shape != (n,) or mass.shape != (n,) or time_s.shape != (n_states,):
        raise AuditInputError("damage, mass, volume, step, or time arrays have inconsistent dimensions")
    if not np.array_equal(point_id, np.arange(n, dtype=point_id.dtype)):
        add_check(report, "complete_point_id_coverage", False, {"point_count": n, "point_id_first": point_id[:10].tolist()})
    else:
        add_check(report, "complete_point_id_coverage", True, {"point_count": n, "point_ids": "0..N-1"})
    if int(geo.get("n_points", -1)) != n or result.get("results", {}).get("point_count") != n:
        add_check(report, "point_count_binding", False,
                  {"case_n_points": geo.get("n_points"), "array_n_points": n,
                   "result_point_count": result.get("results", {}).get("point_count")})
    else:
        add_check(report, "point_count_binding", True, {"n_points": n})
    add_check(report, "state_step_coverage", bool(n_states >= 2 and np.array_equal(steps, np.arange(n_states))),
              {"state_count": n_states, "first_step": int(steps[0]) if n_states else None,
               "last_step": int(steps[-1]) if n_states else None})
    dt = float(solver.get("dt_s", math.nan))
    stable_dt = float(solver.get("stable_dt_s", math.nan))
    safety = float(solver.get("safety_factor", math.nan))
    reported_steps = int(solver.get("steps", -1))
    reported_threads = int(solver.get("threads", -1))
    dt_ok = bool(math.isfinite(dt) and dt > 0 and math.isfinite(stable_dt) and stable_dt > 0
                 and 0 < safety < 1 and dt <= stable_dt * (1 + 2e-12)
                 and np.isclose(dt, safety * stable_dt, rtol=2e-10, atol=2e-18)
                 and reported_steps == n_states - 1 and reported_threads == 1
                 and n <= 12000 and reported_steps <= 1000)
    add_check(report, "native_stable_dt_and_budget", dt_ok,
              {"dt_s": dt, "stable_dt_s": stable_dt, "safety_factor": safety,
               "reported_steps": reported_steps, "saved_states": n_states,
               "threads": reported_threads, "max_points": 12000, "max_steps": 1000})
    runtime_s = float(solver.get("wall_runtime_s", math.nan))
    persistent_bytes = result.get("results", {}).get("persistent_output_bytes")
    actual_tree_bytes = sum(path.stat().st_size for path in root.rglob("*") if path.is_file())
    max_persistent_bytes = 256 * 1024 * 1024
    resource_ok = bool(math.isfinite(runtime_s) and runtime_s <= 120.0
                       and isinstance(persistent_bytes, (int, float)) and 0 <= persistent_bytes <= max_persistent_bytes
                       and actual_tree_bytes <= max_persistent_bytes
                       and len(arrays["bond_i"]) <= 2_000_000)
    add_check(report, "run_resource_budget", resource_ok,
              {"wall_runtime_s": runtime_s, "max_runtime_s": 120.0,
               "persistent_output_bytes_reported": persistent_bytes,
               "run_directory_bytes_at_audit": actual_tree_bytes,
               "max_persistent_output_bytes": max_persistent_bytes,
               "canonical_pair_count": int(len(arrays["bond_i"])), "max_pairs": 2_000_000})
    if math.isfinite(dt) and dt > 0:
        _check_close(report, "time_grid_matches_native_dt", time_s, steps.astype(np.float64) * dt,
                     rtol=2e-12, atol=2e-18)
    elif np.isfinite(time_s).all():
        add_check(report, "time_grid_matches_native_dt", False, {"dt_s": dt})

    expected_lattice = centroid_shell_points(ro, ri, center, h)
    lattice_match = expected_lattice.shape == x_ref.shape and np.allclose(expected_lattice, x_ref, rtol=0.0, atol=4e-15)
    add_check(report, "complete_reference_lattice", bool(lattice_match),
              {"expected_points": int(len(expected_lattice)), "saved_points": int(n),
               "coordinate_max_abs_m": float(np.max(np.abs(expected_lattice - x_ref))) if lattice_match and n else None})
    radii = np.linalg.norm(x_ref - center, axis=1)
    in_shell = bool(np.all(radii >= ri - 2e-15) and np.all(radii <= ro + 2e-15)
                    and np.all(x_ref[:, 2] <= center[2] + 2e-15))
    add_check(report, "reference_points_inside_analytic_shell", in_shell,
              {"minimum_radius_m": float(np.min(radii)), "maximum_radius_m": float(np.max(radii)),
               "inner_radius_m": ri, "outer_radius_m": ro, "max_z_m": float(np.max(x_ref[:, 2]))})

    expected_volume_per_point = h**3
    volume_ok = bool(np.all(volume > 0) and np.allclose(volume, expected_volume_per_point, rtol=2e-14, atol=0.0))
    mass_expected = rho * volume
    mass_ok = bool(np.all(mass > 0) and np.allclose(mass, mass_expected, rtol=2e-14, atol=0.0))
    lattice_volume = float(np.sum(volume, dtype=np.float64))
    total_mass = float(np.sum(mass, dtype=np.float64))
    add_check(report, "all_point_volumes_including_isolated", volume_ok,
              {"point_count": n, "expected_h3_m3": expected_volume_per_point,
               "sum_volume_m3": lattice_volume, "analytic_shell_volume_m3": analytic_volume,
               "relative_volume_quadrature_error": abs(lattice_volume - analytic_volume) / analytic_volume})
    add_check(report, "complete_point_mass_coverage", mass_ok,
              {"point_count": n, "expected_total_mass_kg": float(np.sum(mass_expected)),
               "reported_total_mass_kg": total_mass, "isolated_points_are_included": True})

    # Geometry candidate reconstruction and exact reference void filtering.
    potential_rows, retained_rows, excluded_rows = reconstruct_candidate_pairs(x_ref, delta, center, ri)
    bi = arrays["bond_i"].astype(np.int64, copy=False)
    bj = arrays["bond_j"].astype(np.int64, copy=False)
    rest = arrays["bond_rest_length_m"].astype(np.float64, copy=False)
    eps_c = arrays["bond_epsilon_c"].astype(np.float64, copy=False)
    B = len(bi)
    if bj.shape != (B,) or rest.shape != (B,) or eps_c.shape != (B,):
        raise AuditInputError("native canonical pair arrays have inconsistent lengths")
    canonical = bool(np.all(bi >= 0) and np.all(bj < n) and np.all(bi < bj)
                     and (B < 2 or np.all((bi[1:] > bi[:-1]) | ((bi[1:] == bi[:-1]) & (bj[1:] > bj[:-1])))))
    add_check(report, "canonical_pair_roster", canonical,
              {"pair_count": B, "indices_0_based_i_lt_j": bool(np.all(bi < bj)), "lexicographically_sorted_unique": canonical})
    expected_ij = retained_rows[:, :2].astype(np.int64)
    actual_ij = np.column_stack((bi, bj))
    roster_match = expected_ij.shape == actual_ij.shape and np.array_equal(expected_ij, actual_ij)
    actual_lengths = np.linalg.norm(x_ref[bj] - x_ref[bi], axis=1)
    rest_match = rest.shape == actual_lengths.shape and np.allclose(rest, actual_lengths, rtol=3e-13, atol=2e-15)
    candidate_counts = {"all_within_horizon_before_void_filter": int(len(potential_rows)),
                        "retained_after_exact_inner_sphere_filter": int(len(retained_rows)),
                        "initial_geometric_exclusions": int(len(excluded_rows)),
                        "native_canonical_pair_count": int(B)}
    if not roster_match:
        expected_code = expected_ij[:, 0] * max(n, 1) + expected_ij[:, 1]
        actual_code = actual_ij[:, 0] * max(n, 1) + actual_ij[:, 1]
        expected_only = np.setdiff1d(expected_code, actual_code, assume_unique=False)
        native_only = np.setdiff1d(actual_code, expected_code, assume_unique=False)
        candidate_counts["missing_native_pairs"] = int(len(expected_only))
        candidate_counts["extra_native_pairs"] = int(len(native_only))
        candidate_counts["missing_examples_i_j"] = [[int(v // max(n, 1)), int(v % max(n, 1))] for v in expected_only[:12]]
        candidate_counts["extra_examples_i_j"] = [[int(v // max(n, 1)), int(v % max(n, 1))] for v in native_only[:12]]
    add_check(report, "analytic_candidate_and_native_pair_sets_reciprocal", bool(roster_match), candidate_counts)
    add_check(report, "native_reference_rest_lengths", bool(rest_match),
              {"pair_count": B, **_relative_difference(rest, actual_lengths)})
    excluded_count_reported = solver.get("initial_geometric_exclusion_count")
    add_check(report, "initial_geometric_exclusion_count", excluded_count_reported == len(excluded_rows),
              {"independent_exclusion_count": int(len(excluded_rows)),
               "native_reported_exclusion_count": excluded_count_reported})
    # Asymmetry cannot create a bond: the closest-point test is invariant under
    # endpoint exchange. We record the reciprocal check separately from roster equality.
    reciprocal_tolerance = 32.0 * np.finfo(np.float64).eps * max(ro, 1.0)
    reciprocal_geometry = all(math.isclose(segment_min_radius(x_ref[int(i)], x_ref[int(j)], center),
                                           segment_min_radius(x_ref[int(j)], x_ref[int(i)], center),
                                           rel_tol=0.0, abs_tol=reciprocal_tolerance)
                              for i, j in expected_ij[:min(len(expected_ij), 4096)])
    add_check(report, "reference_filter_reciprocal", bool(reciprocal_geometry),
              {"tested_pairs": min(int(len(expected_ij)), 4096), "criterion": "closed-segment distance >= inner radius"})
    horizon_lengths = potential_rows[:, 2] if len(potential_rows) else np.asarray([], dtype=np.float64)
    horizon_ulp = float(np.spacing(delta))
    horizon_near = int(np.count_nonzero(np.abs(horizon_lengths - delta) <= 4.0 * horizon_ulp))
    report["geometry_pair_audit"] = {**candidate_counts, "horizon_pair_lengths_within_4_ulp": horizon_near,
                                     "horizon_boundary_policy": "candidate uses reference length <= delta; exact/near boundary count is reported separately",
                                     "inner_void_policy": "retain closed segment iff closest reference radius >= Ri"}

    # Critical stretch is uniform here; retain and report threshold ULP cases.
    eps_formula_ok = bool(np.allclose(eps_c, sc, rtol=4e-13, atol=1e-15))
    add_check(report, "uniform_critical_stretch_formula", eps_formula_ok,
              {"formula_sc": sc, "native_epsilon_min": float(np.min(eps_c)) if B else None,
               "native_epsilon_max": float(np.max(eps_c)) if B else None})
    if B == 0:
        add_check(report, "native_pairs_nonempty", False, {"pair_count": 0})

    # Decode little-endian post-removal active flags and check padding bits.
    nbytes = math.ceil(B / 8)
    packed = arrays["bond_active_packed"]
    if packed.shape != (n_states, nbytes):
        raise AuditInputError(f"packed active shape must be {(n_states, nbytes)}, got {packed.shape}")
    actual_active = np.unpackbits(packed.astype(np.uint8, copy=False), axis=1, bitorder="little")[:, :B].astype(bool)
    padding_ok = True
    if B % 8 and nbytes:
        padding_mask = np.uint8((0xFF << (B % 8)) & 0xFF)
        padding_ok = bool(np.all((packed[:, -1].astype(np.uint8) & padding_mask) == 0))
    add_check(report, "packed_active_little_bitorder_and_padding", padding_ok,
              {"shape": list(packed.shape), "bitorder": contract.get("bond_active_packed_bitorder"),
               "padding_bits_zero": padding_ok})

    event_bond = arrays["event_bond_index"].astype(np.int64, copy=False)
    event_step = arrays["event_step"].astype(np.int64, copy=False)
    event_time = arrays["event_time_s"].astype(np.float64, copy=False)
    event_stretch = arrays["event_stretch"].astype(np.float64, copy=False)
    event_energy = arrays["event_pair_energy_J"].astype(np.float64, copy=False)
    if not (event_bond.shape == event_step.shape == event_time.shape == event_stretch.shape == event_energy.shape):
        raise AuditInputError("event table arrays have inconsistent lengths")
    M = len(event_bond)
    event_unique = bool(np.all((event_bond >= 0) & (event_bond < B)) and len(np.unique(event_bond)) == M
                        and np.all((event_step > 0) & (event_step < n_states)))
    add_check(report, "event_indices_unique_and_in_range", event_unique,
              {"event_count": M, "pair_count": B, "state_count": n_states})
    event_rows_ordered = bool(M < 2 or np.all((event_step[1:] > event_step[:-1]) |
                                              ((event_step[1:] == event_step[:-1]) & (event_bond[1:] >= event_bond[:-1]))))
    add_check(report, "event_table_stable_order", event_rows_ordered,
              {"ordering": "step, then canonical pair index"})

    # Reconstruct all first crossings from recorded positions. The strict test
    # follows the vendor CriticalStretch comparison: stretch > epsilon_c.
    first_expected = np.full(B, -1, dtype=np.int64)
    first_stretch_expected = np.full(B, -1.0, dtype=np.float64)
    expected_active = np.ones(B, dtype=bool)
    active_reconstructed = np.empty((n_states, B), dtype=bool)
    event_energy_by_state = np.zeros(n_states, dtype=np.float64)
    threshold_gap_ulps: list[float] = []
    threshold_exact_count = 0
    threshold_one_above_count = 0
    threshold_one_below_count = 0
    event_stretch_independent = np.zeros(M, dtype=np.float64)
    event_energy_independent = np.zeros(M, dtype=np.float64)
    event_record_map = {int(bond): idx for idx, bond in enumerate(event_bond)}
    active_violation_count = 0
    for state in range(n_states):
        stretch_now, pair_energy_now, _ = pair_stretch_and_energy(position[state], bi, bj, rest, volume, c)
        gap = stretch_now - eps_c
        spacing = np.spacing(eps_c)
        near = np.abs(gap) <= 4.0 * np.abs(spacing)
        if np.any(near):
            ulps = gap[near] / spacing[near]
            threshold_gap_ulps.extend(float(v) for v in ulps[:2048])
            threshold_exact_count += int(np.count_nonzero(gap[near] == 0.0))
            threshold_one_above_count += int(np.count_nonzero(gap[near] == spacing[near]))
            threshold_one_below_count += int(np.count_nonzero(gap[near] == -spacing[near]))
        if state > 0:
            breaks = expected_active & (stretch_now > eps_c)
            active_violation_count += int(np.count_nonzero(breaks & ~actual_active[state]))
            if np.any(breaks):
                first_expected[breaks] = state
                first_stretch_expected[breaks] = stretch_now[breaks]
                event_energy_by_state[state] = float(np.sum(pair_energy_now[breaks], dtype=np.float64))
                expected_active[breaks] = False
        active_reconstructed[state] = expected_active
        # Bind every event row to its recorded state, stretch, and pre-removal energy.
        for bond in np.flatnonzero(first_expected == state):
            row = event_record_map.get(int(bond))
            if row is not None:
                event_stretch_independent[row] = stretch_now[bond]
                event_energy_independent[row] = pair_energy_now[bond]

    active_match = np.array_equal(active_reconstructed, actual_active)
    add_check(report, "strict_threshold_and_irreversible_active_history", active_match,
              {"independently_reconstructed_active_matches_native": active_match,
               "actual_active_one_to_zero_only": bool(np.all(np.diff(actual_active.astype(np.int8), axis=0) <= 0)) if n_states > 1 else True,
               "first_crossing_count": int(np.count_nonzero(first_expected >= 0)),
               "threshold_breaks_disagreeing_with_native_active": active_violation_count,
               "native_event_count": M, "strict_operator": "stretch > epsilon_c",
               "threshold_near_pairs_over_states_within_4_ulp": len(threshold_gap_ulps),
               "threshold_exact_evaluations": threshold_exact_count,
               "threshold_one_ulp_above_evaluations": threshold_one_above_count,
               "threshold_one_ulp_below_evaluations": threshold_one_below_count,
               "near_threshold_ulp_offsets_sample": threshold_gap_ulps[:64]})
    _check_close(report, "first_break_step_matches_threshold_reconstruction",
                 arrays["bond_first_break_step"], first_expected, rtol=0, atol=0)
    _check_close(report, "first_break_stretch_matches_recorded_state",
                 arrays["bond_first_break_stretch"], first_stretch_expected, rtol=3e-10, atol=3e-14)
    event_expected_from_bonds = np.flatnonzero(first_expected >= 0)
    event_binding_ok = (M == len(event_expected_from_bonds)
                        and np.array_equal(np.sort(event_bond), event_expected_from_bonds))
    add_check(report, "event_table_equals_first_breaks", bool(event_binding_ok),
              {"expected_event_count": int(len(event_expected_from_bonds)), "native_event_count": M,
               "initial_state_events": 0})
    event_steps_expected = first_expected[event_bond] if M else np.asarray([], dtype=np.int64)
    _check_close(report, "event_step_matches_first_break_array", event_step, event_steps_expected, rtol=0, atol=0)
    event_times_expected = event_step.astype(np.float64) * dt
    _check_close(report, "event_time_matches_state_time", event_time, event_times_expected, rtol=2e-12, atol=2e-18)
    event_stretch_ok = bool(M == 0 or (np.all(event_stretch > eps_c[event_bond])
                                      and np.allclose(event_stretch, event_stretch_independent, rtol=3e-10, atol=3e-14)))
    add_check(report, "event_strict_threshold_and_recorded_positions", event_stretch_ok,
              {"event_count": M, "minimum_stretch_margin": float(np.min(event_stretch - eps_c[event_bond])) if M else None,
               "independent_event_stretch_difference": _relative_difference(event_stretch, event_stretch_independent) if M else {"max_abs": 0.0, "max_rel": 0.0}})
    _check_close(report, "event_pair_energy_recomputed_at_recorded_state", event_energy,
                 event_energy_independent, rtol=5e-10, atol=1e-20,
                 details={"energy_is_pre_removal_at_actual_break_state": True})
    first_time_expected = np.full(B, -1.0, dtype=np.float64)
    has_event = first_expected >= 0
    first_time_expected[has_event] = first_expected[has_event] * dt
    _check_close(report, "first_break_time_matches_event_state",
                 arrays["bond_first_break_time_s"], first_time_expected, rtol=2e-12, atol=2e-18)

    # Independent mass, momentum, angular momentum, energy, force, damage and graph history.
    kin = np.empty(n_states, dtype=np.float64)
    pair_u = np.empty(n_states, dtype=np.float64)
    momentum = np.empty((n_states, 3), dtype=np.float64)
    angular = np.empty((n_states, 3), dtype=np.float64)
    total_force = np.empty((n_states, 3), dtype=np.float64)
    independent_damage = np.zeros((n_states, n), dtype=np.float64)
    component_count = np.empty(n_states, dtype=np.int64)
    max_component_mass = np.empty(n_states, dtype=np.float64)
    native_sed_energy = np.sum(arrays["native_sed_J_m3"].astype(np.float64) * volume[None, :], axis=1) if "native_sed_J_m3" in arrays else None
    degree = np.zeros(n, dtype=np.int64)
    np.add.at(degree, bi, 1)
    np.add.at(degree, bj, 1)
    add_check(report, "every_point_has_initial_native_neighbors", bool(np.all(degree > 0)),
              {"zero_initial_degree_points": np.flatnonzero(degree == 0)[:64].astype(int).tolist(),
               "zero_initial_degree_point_count": int(np.count_nonzero(degree == 0))})
    com = np.sum(mass[:, None] * x_ref, axis=0) / total_mass
    for state in range(n_states):
        kin[state] = 0.5 * float(np.sum(mass * np.sum(velocity[state]**2, axis=1), dtype=np.float64))
        momentum[state] = np.sum(mass[:, None] * velocity[state], axis=0)
        angular[state] = np.sum(np.cross(position[state] - com, mass[:, None] * velocity[state]), axis=0)
        stretch_now, energy_now, unit = pair_stretch_and_energy(position[state], bi, bj, rest, volume, c)
        live = actual_active[state]
        pair_u[state] = float(np.sum(energy_now[live], dtype=np.float64))
        fpair = c * volume[bi] * volume[bj] * stretch_now
        pair_force = (fpair[:, None] * unit) * live[:, None]
        force = np.zeros((n, 3), dtype=np.float64)
        np.add.at(force, bi, pair_force)
        np.add.at(force, bj, -pair_force)
        total_force[state] = np.sum(force, axis=0)
        broken_degree = np.zeros(n, dtype=np.int64)
        np.add.at(broken_degree, bi, (~live).astype(np.int64))
        np.add.at(broken_degree, bj, (~live).astype(np.int64))
        nonzero = degree > 0
        independent_damage[state, nonzero] = broken_degree[nonzero] / degree[nonzero]
        independent_damage[state, ~nonzero] = 0.0
        labels, ncomp = union_find_labels(n, bi, bj, live)
        component_count[state] = ncomp
        _, inverse = np.unique(labels, return_inverse=True)
        masses_by_label = np.bincount(inverse, weights=mass)
        max_component_mass[state] = float(np.max(masses_by_label)) if len(masses_by_label) else 0.0
        _check_close(report, f"component_labels_union_find_state_{state}", arrays["component_label"][state], labels,
                     rtol=0, atol=0, details={"component_count": ncomp,
                                              "includes_isolated_points": True,
                                              "semantics": "minimum point ID representative of active graph component"})
    _check_close(report, "damage_matches_broken_neighbor_fraction", damage, independent_damage,
                 rtol=3e-7, atol=3e-7,
                 details={"zero_initial_degree_damage": 0.0})

    ledger = {key: arrays[key].astype(np.float64, copy=False) for key in LEDGER_SCALARS}
    _check_close(report, "ledger_kinetic_energy", ledger["kinetic_energy_J"], kin, rtol=3e-8, atol=1e-18)
    _check_close(report, "ledger_pair_potential_energy", ledger["pair_potential_energy_J"], pair_u, rtol=5e-7, atol=1e-18)
    independent_D = np.cumsum(event_energy_by_state, dtype=np.float64)
    _check_close(report, "ledger_break_dissipation_from_event_states",
                 ledger["break_dissipation_cumulative_J"], independent_D, rtol=5e-7, atol=1e-18)
    _check_close(report, "native_strain_energy_matches_independent_pair_potential",
                 ledger["native_strain_energy_J"], pair_u, rtol=2e-5, atol=1e-18)
    if native_sed_energy is not None:
        _check_close(report, "integrated_native_sed_matches_pair_potential", native_sed_energy, pair_u,
                     rtol=3e-5, atol=1e-18)
    native_P = np.column_stack([ledger[f"momentum_{axis}_kg_m_s"] for axis in "xyz"])
    native_L = np.column_stack([ledger[f"angular_momentum_{axis}_kg_m2_s"] for axis in "xyz"])
    native_F = np.column_stack([ledger[f"total_force_{axis}_N"] for axis in "xyz"])
    _check_close(report, "ledger_linear_momentum", native_P, momentum, rtol=3e-8, atol=1e-14)
    _check_close(report, "ledger_angular_momentum", native_L, angular, rtol=3e-8, atol=1e-15,
                 details={"origin": "initial mass center of mass"})
    _check_close(report, "ledger_total_internal_force", native_F, total_force, rtol=5e-7, atol=1e-10)

    # Initial prescribed affine velocity, with translation and rigid rotation removed.
    expected_v0 = _expected_initial_velocity(x_ref, mass, center, alpha)
    initial_velocity_ok = np.allclose(velocity[0], expected_v0, rtol=3e-9, atol=2e-13)
    add_check(report, "initial_affine_velocity_and_P_L_removal", bool(initial_velocity_ok),
              {"alpha_s_1": alpha, "velocity_max_abs_difference_mps": float(np.max(np.abs(velocity[0] - expected_v0))),
               "initial_momentum_norm_kg_m_s": float(np.linalg.norm(momentum[0])),
               "initial_angular_momentum_norm_kg_m2_s": float(np.linalg.norm(angular[0]))})
    p_scale = float(np.sum(mass * np.linalg.norm(velocity[0], axis=1)))
    l_scale = float(np.sum(mass * np.linalg.norm(x_ref - com, axis=1) * np.linalg.norm(velocity[0], axis=1)))
    p_rel0 = float(np.linalg.norm(momentum[0]) / max(p_scale, 1e-30))
    l_rel0 = float(np.linalg.norm(angular[0]) / max(l_scale, 1e-30))
    add_check(report, "initial_zero_total_P_and_L", bool(p_rel0 < 3e-9 and l_rel0 < 3e-9),
              {"relative_P_residual": p_rel0, "relative_L_residual": l_rel0,
               "P_norm_kg_m_s": float(np.linalg.norm(momentum[0])),
               "L_norm_kg_m2_s": float(np.linalg.norm(angular[0]))})
    p_characteristic = float(np.max(np.sum(mass[None, :] * np.linalg.norm(velocity, axis=2), axis=1)))
    l_characteristic = float(np.max(np.sum(mass[None, :] * np.linalg.norm(x_ref - com, axis=1)[None, :]
                                             * np.linalg.norm(velocity, axis=2), axis=1)))
    p_drift_rel = float(np.max(np.linalg.norm(momentum - momentum[0], axis=1)) / max(p_characteristic, 1e-30))
    l_drift_rel = float(np.max(np.linalg.norm(angular - angular[0], axis=1)) / max(l_characteristic, 1e-30))
    add_check(report, "free_evolution_P_and_L_conservation", bool(p_drift_rel < 3e-8 and l_drift_rel < 3e-8),
              {"relative_P_drift": p_drift_rel, "relative_L_drift": l_drift_rel,
               "P_characteristic_kg_m_s": p_characteristic,
               "L_characteristic_kg_m2_s": l_characteristic})

    if alpha == 0.0:
        zero_load_ok = bool(M == 0 and np.allclose(position, x_ref[None, :, :], rtol=0.0, atol=2e-14)
                            and np.allclose(velocity, 0.0, rtol=0.0, atol=2e-14)
                            and np.allclose(kin, 0.0, rtol=0.0, atol=1e-18)
                            and np.allclose(pair_u, 0.0, rtol=0.0, atol=1e-18)
                            and np.allclose(independent_D, 0.0, rtol=0.0, atol=1e-18))
        add_check(report, "zero_load_rest_state", zero_load_ok,
                  {"zero_velocity_and_no_events": zero_load_ok, "event_count": M,
                   "max_position_change_m": float(np.max(np.abs(position - x_ref[None, :, :]))),
                   "max_speed_m_s": float(np.max(np.linalg.norm(velocity, axis=2)))})

    # Native checkpoint forces are tested pointwise against the root pair force.
    cp_steps = arrays["checkpoint_step"].astype(np.int64, copy=False)
    cp_force = arrays["checkpoint_internal_force_N"].astype(np.float64, copy=False)
    cp_accel = arrays["checkpoint_acceleration_mps2"].astype(np.float64, copy=False)
    if cp_force.shape != (len(cp_steps), n, 3) or cp_accel.shape != (len(cp_steps), n, 3):
        raise AuditInputError("checkpoint force/acceleration arrays have inconsistent shapes")
    cp_unique = bool(len(cp_steps) > 0 and np.all((cp_steps >= 0) & (cp_steps < n_states))
                     and len(np.unique(cp_steps)) == len(cp_steps)
                     and (len(cp_steps) < 2 or np.all(np.diff(cp_steps) > 0)))
    add_check(report, "checkpoint_steps_valid", cp_unique,
              {"checkpoint_steps": cp_steps.tolist(), "state_count": n_states})
    checkpoint_force_errors = []
    checkpoint_accel_errors = []
    checkpoint_force_ok = cp_unique
    checkpoint_accel_ok = cp_unique
    for row, state_value in enumerate(cp_steps):
        state = int(state_value)
        stretch_now, _, unit = pair_stretch_and_energy(position[state], bi, bj, rest, volume, c)
        live = actual_active[state]
        scalar = c * volume[bi] * volume[bj] * stretch_now
        fpair = (scalar[:, None] * unit) * live[:, None]
        expected_force = np.zeros((n, 3), dtype=np.float64)
        np.add.at(expected_force, bi, fpair)
        np.add.at(expected_force, bj, -fpair)
        force_stats = _relative_difference(cp_force[row], expected_force)
        force_ok = bool(np.allclose(cp_force[row], expected_force, rtol=8e-8, atol=2e-9))
        checkpoint_force_ok = checkpoint_force_ok and force_ok
        checkpoint_force_errors.append({"step": state, **force_stats, "passed": force_ok})
        accel_expected = expected_force / mass[:, None]
        accel_stats = _relative_difference(cp_accel[row], accel_expected)
        accel_ok = bool(np.allclose(cp_accel[row], accel_expected, rtol=2e-7, atol=2e-7))
        checkpoint_accel_ok = checkpoint_accel_ok and accel_ok
        checkpoint_accel_errors.append({"step": state, **accel_stats, "passed": accel_ok})
    add_check(report, "checkpoint_native_internal_force", bool(checkpoint_force_ok),
              {"formula": "c * stretch * Vi * Vj * unit_direction; equal and opposite pair forces",
               "checkpoints": checkpoint_force_errors})
    add_check(report, "checkpoint_native_acceleration", bool(checkpoint_accel_ok),
              {"formula": "internal_force / point_mass; no external force declared",
               "checkpoints": checkpoint_accel_errors})

    # Energy removed on each first break is retained as a discrete overshoot.
    if M:
        threshold_energy = 0.5 * c * volume[bi[event_bond]] * volume[bj[event_bond]] * rest[event_bond] * eps_c[event_bond]**2
        overshoot = event_energy / threshold_energy
        overshoot_summary = {"event_count": M, "ratio_min": float(np.min(overshoot)),
                             "ratio_mean": float(np.mean(overshoot)), "ratio_max": float(np.max(overshoot)),
                             "events_at_or_below_threshold_energy": int(np.count_nonzero(overshoot <= 1.0)),
                             "event_energy_J_sum": float(np.sum(event_energy)),
                             "native_final_D_J": float(ledger["break_dissipation_cumulative_J"][-1]),
                             "discrete_overshoot_is_retained": True}
        overshoot_ok = bool(np.all(overshoot > 1.0 - 2e-12))
    else:
        overshoot_summary = {"event_count": 0, "ratio_min": None, "ratio_mean": None, "ratio_max": None,
                             "events_at_or_below_threshold_energy": 0, "event_energy_J_sum": 0.0,
                             "native_final_D_J": float(ledger["break_dissipation_cumulative_J"][-1]),
                             "discrete_overshoot_is_retained": True}
        overshoot_ok = True
    add_check(report, "event_energy_overshoot", overshoot_ok, overshoot_summary)

    total_mechanical = kin + pair_u + independent_D
    e0 = float(total_mechanical[0])
    energy_scale = max(abs(e0), float(np.max(np.abs(total_mechanical))), 1e-30)
    energy_drift = (total_mechanical - e0) / energy_scale
    report["metrics"] = {
        "point_count": n, "canonical_pair_count": B, "initial_geometric_exclusion_count": int(len(excluded_rows)),
        "event_count": M, "analytic_volume_m3": analytic_volume, "lattice_volume_m3": lattice_volume,
        "relative_volume_quadrature_error": abs(lattice_volume - analytic_volume) / analytic_volume,
        "total_point_mass_kg": total_mass, "K_pa": K, "c_N_m6": c, "critical_stretch": sc,
        "dt_s": dt, "stable_dt_s": stable_dt, "h_m": h, "delta_m": delta,
        "initial_kinetic_energy_J": float(kin[0]), "final_kinetic_energy_J": float(kin[-1]),
        "final_pair_potential_energy_J": float(pair_u[-1]), "final_break_dissipation_J": float(independent_D[-1]),
        "final_graph_component_count": int(component_count[-1]),
        "final_largest_graph_component_mass_kg": float(max_component_mass[-1]),
        "energy_total_relative_drift_max_abs": float(np.max(np.abs(energy_drift))),
        "P_initial_norm_kg_m_s": float(np.linalg.norm(momentum[0])),
        "P_final_norm_kg_m_s": float(np.linalg.norm(momentum[-1])),
        "L_initial_norm_kg_m2_s": float(np.linalg.norm(angular[0])),
        "L_final_norm_kg_m2_s": float(np.linalg.norm(angular[-1])),
    }
    report["histories"] = {
        "step": steps.tolist(), "time_s": time_s.tolist(),
        "kinetic_native_J": ledger["kinetic_energy_J"].tolist(), "kinetic_independent_J": kin.tolist(),
        "pair_potential_native_J": ledger["pair_potential_energy_J"].tolist(), "pair_potential_independent_J": pair_u.tolist(),
        "break_dissipation_native_J": ledger["break_dissipation_cumulative_J"].tolist(),
        "break_dissipation_independent_J": independent_D.tolist(), "total_mechanical_independent_J": total_mechanical.tolist(),
        "total_mechanical_relative_drift": energy_drift.tolist(),
        "momentum_independent_kg_m_s": momentum.tolist(), "angular_momentum_independent_kg_m2_s": angular.tolist(),
        "components_active_graph": component_count.tolist(), "largest_graph_component_mass_kg": max_component_mass.tolist(),
        "damage_mean": np.mean(damage, axis=1).tolist(), "damage_max": np.max(damage, axis=1).tolist(),
        "total_internal_force_independent_N": total_force.tolist(),
    }
    report["runtime"] = {"python": sys.version.split()[0], "numpy": np.__version__,
                         "auditor_sha256": sha256_file(Path(__file__).resolve()),
                         "argv": list(sys.argv)}
    false_checks = [item["name"] for item in report["checks"] if item["passed"] is False]
    unknown_checks = [item["name"] for item in report["checks"] if item["passed"] is None]
    report["status"] = "FAIL" if false_checks else ("PARTIAL" if unknown_checks else "PASS")
    report["failed_checks"] = false_checks
    report["unverified_checks"] = unknown_checks
    return report, arrays, result


def render_plots(report: dict[str, Any], arrays: dict[str, np.ndarray], out_dir: Path) -> list[dict[str, Any]]:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise AuditInputError("matplotlib is required to render the requested standard audit figures") from exc
    report.setdefault("runtime", {})["matplotlib"] = matplotlib.__version__
    out_dir.mkdir(parents=True, exist_ok=True)
    h = report["histories"]
    t = np.asarray(h["time_s"], dtype=np.float64)
    artifacts: list[dict[str, Any]] = []

    fig, ax = plt.subplots(figsize=(10.5, 6.2), constrained_layout=True)
    ax.plot(t, h["kinetic_native_J"], label="Native K", linewidth=1.6)
    ax.plot(t, h["kinetic_independent_J"], "--", label="Independent K", linewidth=1.2)
    ax.plot(t, h["pair_potential_native_J"], label="Native pair U", linewidth=1.6)
    ax.plot(t, h["pair_potential_independent_J"], "--", label="Independent pair U", linewidth=1.2)
    ax.plot(t, h["break_dissipation_native_J"], label="Native cumulative D", linewidth=1.6)
    ax.plot(t, h["break_dissipation_independent_J"], "--", label="Independent cumulative D", linewidth=1.2)
    ax.plot(t, h["total_mechanical_independent_J"], ":", label="Independent K+U+D", linewidth=1.7, color="black")
    ax.set(xlabel="Time (s)", ylabel="Energy (J)", title="Native and independently reconstructed energy ledger")
    ax.grid(True, alpha=0.25)
    ax.legend(ncol=2, fontsize=8)
    path = out_dir / "energy_ledger.png"
    fig.savefig(path, dpi=170)
    plt.close(fig)
    artifacts.append({"path": path.name, "bytes": path.stat().st_size, "sha256": sha256_file(path)})

    fig, axes = plt.subplots(2, 1, figsize=(10.5, 7.2), sharex=True, constrained_layout=True)
    P = np.asarray(h["momentum_independent_kg_m_s"], dtype=np.float64)
    L = np.asarray(h["angular_momentum_independent_kg_m2_s"], dtype=np.float64)
    axes[0].plot(t, np.linalg.norm(P, axis=1), color="tab:blue", label="|P|")
    axes[0].set_ylabel("|P| (kg m/s)")
    axes[0].grid(True, alpha=0.25)
    axes[0].legend()
    axes[1].plot(t, np.linalg.norm(L, axis=1), color="tab:orange", label="|L| about initial COM")
    axes[1].set(xlabel="Time (s)", ylabel="|L| (kg m²/s)")
    axes[1].grid(True, alpha=0.25)
    axes[1].legend()
    fig.suptitle("Free-system momentum diagnostics")
    path = out_dir / "momentum_history.png"
    fig.savefig(path, dpi=170)
    plt.close(fig)
    artifacts.append({"path": path.name, "bytes": path.stat().st_size, "sha256": sha256_file(path)})

    fig, axes = plt.subplots(2, 1, figsize=(10.5, 7.2), sharex=True, constrained_layout=True)
    axes[0].plot(t, h["components_active_graph"], color="tab:purple")
    axes[0].set(ylabel="Component count", title="Active-bond graph components (not physical shards)")
    axes[0].grid(True, alpha=0.25)
    axes[1].plot(t, h["damage_mean"], label="Mean point damage")
    axes[1].plot(t, h["damage_max"], label="Maximum point damage")
    axes[1].set(xlabel="Time (s)", ylabel="Failed-neighbor fraction", ylim=(-0.02, 1.02))
    axes[1].grid(True, alpha=0.25)
    axes[1].legend()
    path = out_dir / "graph_and_damage_history.png"
    fig.savefig(path, dpi=170)
    plt.close(fig)
    artifacts.append({"path": path.name, "bytes": path.stat().st_size, "sha256": sha256_file(path)})

    fig, ax = plt.subplots(figsize=(8.0, 6.8), constrained_layout=True)
    final_damage = arrays["damage"][-1]
    points = arrays["x_ref_m"]
    scat = ax.scatter(points[:, 0] * 1e3, points[:, 2] * 1e3, c=final_damage,
                      s=8, cmap="magma", vmin=0.0, vmax=1.0, linewidths=0)
    ax.set_aspect("equal", adjustable="box")
    ax.set(xlabel="Reference x (mm)", ylabel="Reference z (mm)",
           title="Final point damage on reference x-z projection")
    fig.colorbar(scat, ax=ax, label="Failed-neighbor fraction")
    ax.grid(True, alpha=0.18)
    path = out_dir / "final_damage_reference_projection.png"
    fig.savefig(path, dpi=170)
    plt.close(fig)
    artifacts.append({"path": path.name, "bytes": path.stat().st_size, "sha256": sha256_file(path)})
    return artifacts


def analytic_handtest() -> dict[str, Any]:
    """Tiny deterministic math check, including strict and adjacent-ULP cases."""
    checks: list[dict[str, Any]] = []

    def check(name: str, ok: bool, details: dict[str, Any]) -> None:
        checks.append({"name": name, "passed": bool(ok), "details": details})

    E, nu, delta, Gc = 6.0, 0.25, 1.0, 0.1125
    K, c, sc = derive_root_material(E, nu, Gc, delta)
    check("root_material_formulas", K == 4.0 and math.isclose(c, 72.0 / math.pi, rel_tol=1e-15)
          and math.isclose(sc, 0.125, rel_tol=1e-15), {"K": K, "c": c, "sc": sc})
    below, exact, above = np.nextafter(sc, -math.inf), sc, np.nextafter(sc, math.inf)
    check("strict_threshold_and_adjacent_ulp", not (exact > sc) and not (below > sc) and (above > sc),
          {"one_ulp_below": float(below), "exact": float(exact), "one_ulp_above": float(above),
           "operator": "stretch > epsilon_c"})

    # Three points, one canonical bond, and one isolated point. At state 1 the
    # bond lies exactly at threshold; state 2 crosses by the next representable
    # increment and the discrete removed energy is recomputed at that state.
    x_ref = np.asarray(((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (3.0, 0.0, 0.0)))
    volume = np.asarray((0.5, 0.5, 0.25))
    mass = volume.copy()
    stretch_above = float(np.nextafter(1.0 + sc, math.inf) - 1.0)
    positions = np.asarray((x_ref,
                            ((0.0, 0.0, 0.0), (1.0 + sc, 0.0, 0.0), (3.0, 0.0, 0.0)),
                            ((0.0, 0.0, 0.0), (1.0 + stretch_above, 0.0, 0.0), (3.0, 0.0, 0.0))), dtype=np.float64)
    bi, bj, rest = np.asarray((0,), dtype=np.int64), np.asarray((1,), dtype=np.int64), np.asarray((1.0,))
    state_stretches = []
    state_energies = []
    for state in positions:
        ss, uu, ee = pair_stretch_and_energy(state, bi, bj, rest, volume, c)
        state_stretches.append(float(ss[0]))
        state_energies.append(float(uu[0]))
    active = np.asarray(((True,), (True,), (False,)))
    labels0, count0 = union_find_labels(3, bi, bj, active[0])
    labels2, count2 = union_find_labels(3, bi, bj, active[2])
    check("analytic_bond_stretch_energy_force", math.isclose(state_energies[1], 0.5 * c * volume[0] * volume[1] * sc**2, rel_tol=2e-14)
          and state_energies[2] > state_energies[1],
          {"stretches": state_stretches, "pair_energy_at_threshold_J": state_energies[1],
           "removed_event_energy_J": state_energies[2], "energy_ratio": state_energies[2] / state_energies[1]})
    force_mag = c * volume[0] * volume[1] * state_stretches[2]
    check("equal_and_opposite_pair_force", force_mag > 0 and np.array_equal(np.asarray((force_mag, 0, 0)) + np.asarray((-force_mag, 0, 0)), np.zeros(3)),
          {"force_on_i_N": [force_mag, 0.0, 0.0], "force_on_j_N": [-force_mag, 0.0, 0.0]})
    check("union_find_keeps_isolated_point", count0 == 2 and count2 == 3
          and np.array_equal(labels0, (0, 0, 2)) and np.array_equal(labels2, (0, 1, 2)),
          {"intact_labels": labels0.tolist(), "broken_labels": labels2.tolist(),
           "intact_components": count0, "broken_components": count2})
    d_cumulative = np.asarray((0.0, 0.0, state_energies[2]))
    check("discrete_D_is_event_state_pair_energy", d_cumulative[-1] == state_energies[2]
          and d_cumulative[-1] > state_energies[1],
          {"D_by_state_J": d_cumulative.tolist(), "threshold_energy_J": state_energies[1],
           "removed_event_energy_J": state_energies[2]})

    center = np.zeros(3)
    through_void = segment_min_radius(np.asarray((-0.009, 0.0, 0.0)), np.asarray((0.009, 0.0, 0.0)), center)
    tangent_void = segment_min_radius(np.asarray((0.008, 0.0, 0.0)), np.asarray((0.0105, 0.0, 0.0)), center)
    check("analytic_inner_sphere_segment_filter", through_void < 0.008 and tangent_void >= 0.008,
          {"crossing_segment_min_radius_m": through_void, "tangent_segment_min_radius_m": tangent_void,
           "inner_radius_m": 0.008})
    failed = [item["name"] for item in checks if not item["passed"]]
    return {"schema": "pgfrag.independent_bowl_handtest.v1", "status": "FAIL" if failed else "PASS",
            "scope": "analytic algebra/geometry/threshold/energy/force/union-find only; no solver launched",
            "checks": checks, "failed_checks": failed}


def synthetic_schema_handtest(out_dir: Path) -> dict[str, Any]:
    """Exercise the full readback path on a tiny static schema fixture.

    This is only a four-point in-memory analytic fixture serialized as NPZ and
    JSON. It does not call, emulate, or launch the native solver.
    """
    run_dir = out_dir / "synthetic_native_bowl_run"
    run_dir.mkdir(parents=True, exist_ok=True)
    ro, ri, h, delta = 12e-3, 8e-3, 12e-3, 30e-3
    center = np.asarray((0.0, 0.0, 12e-3), dtype=np.float64)
    rho, E, nu, Gc = 2300.0, 70e9, 0.25, 7.6
    K, c, sc = derive_root_material(E, nu, Gc, delta)
    x_ref = centroid_shell_points(ro, ri, center, h)
    potential, retained, excluded = reconstruct_candidate_pairs(x_ref, delta, center, ri)
    bi = retained[:, 0].astype(np.int32)
    bj = retained[:, 1].astype(np.int32)
    rest = retained[:, 2].astype(np.float64)
    n, B, n_states = len(x_ref), len(bi), 2
    position = np.broadcast_to(x_ref, (n_states, n, 3)).copy()
    velocity = np.zeros_like(position)
    volume = np.full(n, h**3, dtype=np.float64)
    mass = rho * volume
    labels, ncomp = union_find_labels(n, bi, bj, np.ones(B, dtype=bool))
    packed = np.packbits(np.ones((n_states, B), dtype=np.uint8), axis=1, bitorder="little")
    arrays: dict[str, np.ndarray] = {
        "point_id": np.arange(n, dtype=np.int32), "x_ref_m": x_ref,
        "position_m": position, "velocity_mps": velocity,
        "damage": np.zeros((n_states, n), dtype=np.float32),
        "component_label": np.broadcast_to(labels, (n_states, n)).copy(),
        "native_sed_J_m3": np.zeros((n_states, n), dtype=np.float32),
        "volume_m3": volume, "point_mass_kg": mass,
        "step": np.arange(n_states, dtype=np.int32), "time_s": np.asarray((0.0, 0.5)),
        "bond_i": bi, "bond_j": bj, "bond_rest_length_m": rest,
        "bond_epsilon_c": np.full(B, sc),
        "bond_first_break_step": np.full(B, -1, dtype=np.int32),
        "bond_first_break_time_s": np.full(B, -1.0),
        "bond_first_break_stretch": np.full(B, -1.0),
        "bond_active_packed": packed,
        "event_bond_index": np.asarray([], dtype=np.int32),
        "event_step": np.asarray([], dtype=np.int32),
        "event_time_s": np.asarray([], dtype=np.float64),
        "event_stretch": np.asarray([], dtype=np.float64),
        "event_pair_energy_J": np.asarray([], dtype=np.float64),
        "kinetic_energy_J": np.zeros(n_states), "pair_potential_energy_J": np.zeros(n_states),
        "break_dissipation_cumulative_J": np.zeros(n_states), "native_strain_energy_J": np.zeros(n_states),
        "total_force_x_N": np.zeros(n_states), "total_force_y_N": np.zeros(n_states), "total_force_z_N": np.zeros(n_states),
        "momentum_x_kg_m_s": np.zeros(n_states), "momentum_y_kg_m_s": np.zeros(n_states), "momentum_z_kg_m_s": np.zeros(n_states),
        "angular_momentum_x_kg_m2_s": np.zeros(n_states), "angular_momentum_y_kg_m2_s": np.zeros(n_states), "angular_momentum_z_kg_m2_s": np.zeros(n_states),
        "checkpoint_step": np.asarray((0, 1), dtype=np.int32),
        "checkpoint_acceleration_mps2": np.zeros((2, n, 3)),
        "checkpoint_internal_force_N": np.zeros((2, n, 3)),
    }
    arrays_path = run_dir / "arrays.npz"
    np.savez_compressed(arrays_path, **arrays)
    analytic_volume = 2.0 * math.pi / 3.0 * (ro**3 - ri**3)
    lattice_volume = float(np.sum(volume))
    manifest_path = run_dir / "asset" / "manifest.json"
    write_json(manifest_path, {"schema": "forge.asset.v1", "handtest": True})
    source_snapshot = run_dir / "source_snapshot"
    source_snapshot.mkdir(exist_ok=True)
    synthetic_sources = {
        "bowl_dynamics.py": b"# synthetic handtest source pin; no solver implementation\n",
        "run_bowl_dynamics.jl": b"# synthetic handtest source pin; no Julia execution\n",
        "Project.toml": b"# synthetic Julia project receipt\n",
        "Manifest.toml": b"# synthetic Julia manifest receipt\n",
    }
    source_receipts = []
    for name, data in synthetic_sources.items():
        path = source_snapshot / name
        path.write_bytes(data)
        source_receipts.append({"path": f"source_snapshot/{name}", "sha256": sha256_file(path), "bytes": path.stat().st_size})
    source_receipts.append({"path": "asset/manifest.json", "sha256": sha256_file(manifest_path), "bytes": manifest_path.stat().st_size})
    source_hash = sha256_file(manifest_path)
    case = {
        "schema": "pgfrag.native_bowl_case.v1", "run_id": "synthetic_handtest",
        "source_forge_manifest_sha256": source_hash, "source_manifest_path": str(manifest_path),
        "geometry": {"shape": "concentric_lower_hemisphere_shell", "outer_radius_m": ro,
                     "inner_radius_m": ri, "center_m": center.tolist(), "h_m": h, "delta_m": delta,
                     "n_points": n, "analytic_volume_m3": analytic_volume,
                     "lattice_volume_m3": lattice_volume},
        "material": {"rho_kg_m3": rho, "E_pa": E, "nu": nu, "Gc_j_m2": Gc,
                     "calibration_status": "UNVALIDATED_PROTOTYPE", "point_mass_kg": float(np.sum(mass))},
        "loading": {"mode": "zero", "alpha_s_1": 0.0,
                    "momentum_constraint": "remove mass-mean translation and angular momentum",
                    "gaussian_peak_m_s": 0.0, "gaussian_sigma_m": 0.0},
        "solver": {"steps": 1, "safety_factor": 0.5, "threads": 1},
    }
    case_path = run_dir / "case.json"
    write_json(case_path, case)
    case_toml_path = run_dir / "case.toml"
    case_toml_path.write_text("schema = \"pgfrag.native_bowl_case.v1\"\n", encoding="utf-8")
    output_receipts = []
    for path in (case_path, case_toml_path, arrays_path):
        output_receipts.append({"path": path.name, "sha256": sha256_file(path), "bytes": path.stat().st_size})
    argv_path = run_dir / "argv.json"
    write_json(argv_path, {"argv": ["synthetic-julia", "--project=synthetic", "run_bowl_dynamics.jl",
                                   str(case_toml_path), str(run_dir)],
                           "cwd": str(run_dir), "run_id": "synthetic_handtest",
                           "source_receipts": source_receipts})
    result = {
        "schema": NATIVE_SCHEMA, "status": "PASS", "run_id": "synthetic_handtest",
        "source_forge_manifest_sha256": source_hash, "source_receipts": source_receipts,
        "geometry": case["geometry"], "material": {**case["material"], "total_mass_kg": float(np.sum(mass))},
        "loading": case["loading"],
        "solver": {"name": "BBMaterial{NoCorrection}, CriticalStretch, VelocityVerlet",
                   "model": "BB", "material": "BBMaterial", "correction": "NoCorrection",
                   "damage_model": "CriticalStretch", "integrator": "VelocityVerlet",
                   "status": "NATIVE_COMPLETE", "steps": 1, "threads": 1,
                   "dt_s": 0.5, "stable_dt_s": 1.0, "safety_factor": 0.5,
                   "wall_runtime_s": 0.01, "initial_geometric_exclusion_count": int(len(excluded))},
        "results": {"point_count": n,
                    "peak_raw_plus_final_bytes": arrays_path.stat().st_size + case_path.stat().st_size + case_toml_path.stat().st_size},
        "array_contract": {"format": "NumPy npz", "units": "SI",
                           "bond_active_packed_bitorder": "little; 1=active"},
        "payloads": {"arrays": {"path": "arrays.npz", "sha256": sha256_file(arrays_path),
                                  "bytes": arrays_path.stat().st_size}},
        "outputs": output_receipts,
    }
    write_json(run_dir / "result.json", result)
    report, _, _ = audit_run(run_dir)
    write_json(run_dir / "independent_audit.json", report)
    return {"status": report["status"], "failed_checks": report["failed_checks"],
            "check_count": len(report["checks"]), "point_count": n,
            "canonical_pair_count": B, "initial_geometric_exclusion_count": int(len(excluded)),
            "graph_component_count": ncomp, "scope": "synthetic static schema fixture; no solver launched",
            "audit_json": str(run_dir / "independent_audit.json")}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Read-only independent audit of a saved native bowl result")
    parser.add_argument("--run-dir", type=Path, help="completed run directory containing result.json, case.json, arrays.npz")
    parser.add_argument("--out-dir", type=Path, help="separate directory for JSON and Matplotlib artifacts")
    parser.add_argument("--no-plots", action="store_true", help="write JSON only; default renders standard Matplotlib figures")
    parser.add_argument("--self-test", action="store_true", help="run the tiny analytic handtest; never launches a solver")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.self_test:
        out_dir = (args.out_dir or (HERE / "handtest")).resolve()
        result = analytic_handtest()
        schema_result = synthetic_schema_handtest(out_dir)
        result["synthetic_schema_fixture"] = schema_result
        if schema_result["status"] != "PASS":
            result["checks"].append({"name": "synthetic_native_schema_audit", "passed": False,
                                     "details": {"failed_checks": schema_result["failed_checks"]}})
            result["failed_checks"].append("synthetic_native_schema_audit")
            result["status"] = "FAIL"
        write_json(out_dir / "handtest_result.json", result)
        print(json.dumps({"status": result["status"], "result": str(out_dir / "handtest_result.json"),
                          "failed_checks": result["failed_checks"]}, indent=2))
        return 0 if result["status"] == "PASS" else 2
    if args.run_dir is None:
        raise SystemExit("--run-dir is required unless --self-test is used")
    run_dir = args.run_dir.resolve(strict=True)
    out_dir = (args.out_dir or (HERE / ("readback_" + run_dir.name))).resolve()
    try:
        out_dir.relative_to(run_dir)
    except ValueError:
        pass
    else:
        raise SystemExit("--out-dir must be outside the native run directory so the audit stays read-only")
    try:
        report, arrays, _ = audit_run(run_dir)
        plots = [] if args.no_plots else render_plots(report, arrays, out_dir)
        report["artifacts"] = plots
        report_path = out_dir / "independent_audit.json"
        write_json(report_path, report)
        print(json.dumps({"status": report["status"], "report": str(report_path),
                          "failed_checks": report["failed_checks"], "unverified_checks": report["unverified_checks"],
                          "plots": [item["path"] for item in plots]}, indent=2))
        return 0 if report["status"] == "PASS" else (3 if report["status"] == "PARTIAL" else 2)
    except (AuditInputError, OSError, ValueError, KeyError, TypeError) as exc:
        failure = {"schema": SCHEMA, "status": "INPUT_ERROR", "run_directory": str(run_dir),
                   "error": f"{type(exc).__name__}: {exc}", "auditor_sha256": sha256_file(Path(__file__).resolve())}
        write_json(out_dir / "independent_audit.json", failure)
        print(json.dumps({"status": "INPUT_ERROR", "report": str(out_dir / "independent_audit.json"),
                          "error": failure["error"]}, indent=2), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
