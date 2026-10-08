"""Read pinned static replay arrays and form the consistent 3-D mass matrix."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CASE_PATH = HERE / "fixed_case.json"
REVIEW_ROOT = ROOT / "pgfrag" / "tmp" / "review_20261006" / "full3d_consistent_mass_replay_20261006"
CASE_ID = "modeiii_pu_full3d_consistent_mass_replay_20261006"
OPERATOR_CASE_ID = "modeiii_pu_full3d_elastic_operator_20261006"
SOURCE_COUNT, BASIS_COUNT, DOF_COUNT, VOLUME_COUNT, QUERY_COUNT = 27, 54, 162, 256, 384
MAX_BYTES = 8 * 1024 * 1024
COMPONENTS = ("ux", "uy", "uz")
BASIS_KINDS = ("base", "front_f1_unshifted")
COMPONENT_INDEX = {name: i for i, name in enumerate(COMPONENTS)}


class ReplayBlocked(RuntimeError):
    """Fail-closed source, input, or numerical-contract error."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _inside(path: Path, parent: Path) -> bool:
    try:
        path.resolve().relative_to(parent.resolve())
        return True
    except ValueError:
        return False


def _shape(value: object) -> tuple[int, ...] | None:
    if not isinstance(value, list): return ()
    if not value: return (0,)
    children = [_shape(item) for item in value]
    if any(item is None for item in children) or any(item != children[0] for item in children): return None
    return (len(value),) + children[0]


def _finite(value: object) -> bool:
    if isinstance(value, bool): return False
    if isinstance(value, (int, float)): return math.isfinite(float(value))
    if isinstance(value, list): return all(_finite(item) for item in value)
    return False


def _read_json(path: Path, cap: int) -> dict[str, object]:
    if not path.is_file() or path.stat().st_size > cap: raise ReplayBlocked(f"missing or oversized JSON: {path}")
    try: value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc: raise ReplayBlocked(f"invalid JSON: {exc}") from exc
    if not isinstance(value, dict): raise ReplayBlocked("JSON root must be an object")
    return value


def _pins(root: Path, path: Path) -> tuple[dict[str, str], str]:
    pins = _read_json(path, 1024 * 1024)
    if pins.get("case_id") != OPERATOR_CASE_ID: raise ReplayBlocked("operator pin case ID differs")
    deps, entry = pins.get("dependencies"), pins.get("entrypoint")
    if not isinstance(deps, list) or not isinstance(entry, dict): raise ReplayBlocked("operator pins are incomplete")
    expected: dict[str, str] = {}
    for item in [entry, *deps]:
        if not isinstance(item, dict): raise ReplayBlocked("malformed operator source pin")
        rel, digest = item.get("path"), item.get("sha256")
        name = Path(str(rel)).name if item is entry else item.get("name")
        if type(name) is not str or type(rel) is not str or not re.fullmatch(r"[0-9a-f]{64}", str(digest)):
            raise ReplayBlocked("invalid operator source pin")
        source = (root / rel).resolve()
        if not _inside(source, root) or not source.is_file() or _sha256(source) != digest:
            raise ReplayBlocked(f"operator source pin mismatch: {rel}")
        if item is not entry:
            if name in expected: raise ReplayBlocked("duplicate source pin name")
            expected[name] = digest
    return expected, _sha256(path)


def _registry() -> tuple[list[dict[str, object]], list[str], list[dict[str, str]]]:
    sources = [f"node-{i:03d}" for i in range(SOURCE_COUNT)]
    rows = [{"source_id": source, "basis_kind": basis} for source in sources for basis in BASIS_KINDS]
    records, ids = [], []
    for source in sources:
        for component in COMPONENTS:
            for basis, index in (("base", 0), ("front_f1_unshifted", 1)):
                canonical = json.dumps([source, component, basis, index], separators=(",", ":"))
                ids.append(canonical)
                records.append({"slot": len(records), "canonical_id": canonical, "source_id": source,
                                "component": component, "basis_kind": basis, "basis_index": index})
    return records, ids, rows


def _preflight(result: dict[str, object], expected: dict[str, str]) -> dict[str, object]:
    if result.get("case_id") != OPERATOR_CASE_ID: raise ReplayBlocked("operator case ID differs")
    if (result.get("status") != "OK_DIAGNOSTIC_ONLY_UNQUALIFIED" or result.get("assembly_status") != "OK"
            or result.get("mechanical_eligibility") != "DIAGNOSTIC_PASS_UNQUALIFIED" or result.get("solve_allowed") is not False):
        raise ReplayBlocked("operator gate did not pass with solve disabled")
    if result.get("source_sha256") != expected: raise ReplayBlocked("operator result source pins differ")
    spectrum = result.get("spectrum")
    if not isinstance(spectrum, dict): raise ReplayBlocked("operator result has no spectrum")
    eig = spectrum.get("symmetric_eigenvalues")
    if (not isinstance(eig, list) or len(eig) != DOF_COUNT or not _finite(eig) or spectrum.get("gate_passed") is not True
            or spectrum.get("observed_near_zero_count") != 6 or spectrum.get("observed_positive_count") != 156
            or spectrum.get("significant_negative_count") != 0):
        raise ReplayBlocked("operator spectrum did not pass its six/156 gate")
    replay = result.get("static_solve_replay")
    if not isinstance(replay, dict) or replay.get("schema") != "static_solve_replay.v1":
        raise ReplayBlocked("static_solve_replay.v1 is missing")
    if replay.get("case_id") != OPERATOR_CASE_ID or replay.get("upstream_source_sha256") != expected:
        raise ReplayBlocked("replay case or source pins differ")
    if (replay.get("array_encoding") != "JSON nested arrays; row-major in declared dimensions"
            or replay.get("query_order") != ["volume", "lateral_side", "top_cap", "bottom_cap"]
            or replay.get("rigid_mode_order") != ["tx", "ty", "tz", "rx", "ry", "rz"]):
        raise ReplayBlocked("replay array, query, or rigid-mode order differs")
    units = replay.get("array_units")
    if not isinstance(units, dict) or units.get("basis.values") != "dimensionless" or units.get("quadrature.volume_weights_m3") != "m^3":
        raise ReplayBlocked("cached basis or quadrature units are missing")
    records, ids, rows = _registry()
    reg = replay.get("registry")
    if not isinstance(reg, dict) or reg.get("source_ids") != [f"node-{i:03d}" for i in range(SOURCE_COUNT)]:
        raise ReplayBlocked("source IDs differ from the fixed source order")
    if reg.get("ordered_canonical_ids") != ids or reg.get("ordered_keys") != records or replay.get("basis_row_ids") != rows:
        raise ReplayBlocked("registry and basis rows do not match the semantic UX/UY/UZ mapping")
    op, src, quad, basis = (replay.get(key) for key in ("operator", "sources", "quadrature", "basis"))
    if not all(isinstance(block, dict) for block in (op, src, quad, basis)):
        raise ReplayBlocked("replay data blocks are incomplete")
    required = ((src, "positions_xyz_m", (27, 3)), (src, "rigid_mode_center_xyz_m", (3,)),
                (op, "rigid_mode_coefficients", (162, 6)), (quad, "positions_xyz_m", (384, 3)),
                (quad, "volume_weights_m3", (256,)), (basis, "values", (384, 54)))
    for block, key, shape in required:
        if _shape(block.get(key)) != shape or not _finite(block.get(key)):
            raise ReplayBlocked(f"cached {key} must be finite with shape {shape}")
    if any(float(w) <= 0.0 for w in quad["volume_weights_m3"]): raise ReplayBlocked("volume weights must be positive")
    regions = {"volume": [0, 256], "lateral_side": [256, 320], "top_cap": [320, 352], "bottom_cap": [352, 384]}
    if quad.get("region_slices_half_open") != regions: raise ReplayBlocked("volume region slice differs")
    return replay


def _calculate(replay: dict[str, object], density: float, input_path: Path, input_sha: str,
               pins_sha: str, case_sha: str, started: float) -> tuple[dict[str, object], list[str]]:
    import numpy as np
    registry = replay["registry"]; records = registry["ordered_keys"]
    row_lookup = {(r["source_id"], r["basis_kind"]): i for i, r in enumerate(replay["basis_row_ids"])}
    scalar_slots = np.asarray([row_lookup[(r["source_id"], r["basis_kind"])] for r in records])
    components = np.asarray([COMPONENT_INDEX[r["component"]] for r in records])
    values = np.asarray(replay["basis"]["values"], dtype=np.float64)
    weights = np.asarray(replay["quadrature"]["volume_weights_m3"], dtype=np.float64)
    positions = np.asarray(replay["quadrature"]["positions_xyz_m"][:VOLUME_COUNT], dtype=np.float64)
    sources = np.asarray(replay["sources"]["positions_xyz_m"], dtype=np.float64)
    center = np.asarray(replay["sources"]["rigid_mode_center_xyz_m"], dtype=np.float64)
    modes = np.asarray(replay["operator"]["rigid_mode_coefficients"], dtype=np.float64)
    volume = float(np.sum(weights)); mass_scale = density * volume
    weighted = values[:VOLUME_COUNT] * np.sqrt(density * weights)[:, None]
    scalar_mass = weighted.T @ weighted
    same = components[:, None] == components[None, :]
    mass = scalar_mass[scalar_slots[:, None], scalar_slots[None, :]] * same
    norm = float(np.linalg.norm(mass, "fro")); sym = float(np.linalg.norm(mass - mass.T, "fro"))
    sym_rel = sym / max(norm, np.finfo(float).tiny)
    slots = {name: np.flatnonzero(components == i) for i, name in enumerate(COMPONENTS)}
    blocks = {name: mass[np.ix_(ix, ix)] for name, ix in slots.items()}
    traces = {name: float(np.trace(block)) for name, block in blocks.items()}
    scalar_trace, trace = float(np.trace(scalar_mass)), float(np.trace(mass))
    trace_rel = abs(trace - 3 * scalar_trace) / max(abs(3 * scalar_trace), np.finfo(float).tiny)
    block_rel = max(float(np.linalg.norm(blocks[k] - blocks["ux"], "fro")) /
                    max(float(np.linalg.norm(blocks["ux"], "fro")), np.finfo(float).tiny) for k in COMPONENTS)
    cross_max = float(np.max(np.abs(mass[~same])))
    scaled = mass / mass_scale; eig = np.linalg.eigvalsh(scaled)
    zero_tol = 1e-12 * max(float(np.max(np.abs(eig))), np.finfo(float).tiny)
    zeros = int(np.count_nonzero(np.abs(eig) <= zero_tol)); negatives = int(np.count_nonzero(eig < -zero_tol))
    positives = int(np.count_nonzero(eig > zero_tol)); spectrum_ok = zeros == negatives == 0 and positives == DOF_COUNT
    chol_ok, chol_error = False, None
    if spectrum_ok:
        try: np.linalg.cholesky(scaled); chol_ok = True
        except np.linalg.LinAlgError as exc: chol_error = str(exc)
    physical = np.zeros((VOLUME_COUNT, 3, 6))
    eye = np.eye(3)
    for axis in range(3):
        physical[:, axis, axis] = 1.0
        physical[:, :, 3 + axis] = np.cross(np.broadcast_to(eye[axis], positions.shape), positions - center)
    coef_rows = np.zeros((BASIS_COUNT, 3, 6))
    for slot, key in enumerate(records):
        coef_rows[row_lookup[(key["source_id"], key["basis_kind"])], COMPONENT_INDEX[key["component"]]] = modes[slot]
    reconstructed = np.einsum("qr,rcm->qcm", values[:VOLUME_COUNT], coef_rows, optimize=True)
    mode_scale = np.asarray([1., 1., 1., 100., 100., 100.])
    recon_err = float(np.max(np.abs((reconstructed - physical) * mode_scale[None, None, :])))
    direct = np.einsum("q,qcm,qcn->mn", density * weights, physical, physical, optimize=True)
    from_coeff = modes.T @ mass @ modes
    dscaled = direct * mode_scale[:, None] * mode_scale[None, :] / mass_scale
    cscaled = from_coeff * mode_scale[:, None] * mode_scale[None, :] / mass_scale
    rigid_err = float(np.max(np.abs(cscaled - dscaled)))
    translation = np.diag(direct)[:3]
    translation_err = float(np.max(np.abs(translation - mass_scale)) / mass_scale)
    r, length = .01, .01
    inertia_expected = np.asarray([mass_scale * (r*r/4 + length*length/12)] * 2 + [mass_scale*r*r/2])
    inertia = np.diag(direct)[3:]
    inertia_err = float(np.max(np.abs(inertia - inertia_expected) / inertia_expected))
    centroid_err = float(np.max(np.abs(np.mean(sources, axis=0) - center)))
    analytic_volume = math.pi * r*r * length
    volume_err = abs(volume - analytic_volume) / analytic_volume
    checks = {"symmetry": sym_rel <= 1e-12, "component_and_trace_consistency": block_rel <= 1e-12 and trace_rel <= 1e-12,
              "cross_component_zero": cross_max == 0.0, "spectrum_no_negative_or_zero_modes": spectrum_ok,
              "unmodified_cholesky": chol_ok, "rigid_field_reconstruction": recon_err <= 1e-8,
              "all_six_kinetic_forms": rigid_err <= 1e-8, "translation_mass_rhoV": translation_err <= 1e-10,
              "quadrature_rotation_inertia": inertia_err <= 1e-10, "volume_matches_fixed_cylinder": volume_err <= 1e-12,
              "saved_rotation_center_is_source_centroid": centroid_err <= 1e-12}
    passed = all(checks.values())
    result = {
        "schema": "full3d_consistent_mass_replay_result.v1", "case_id": CASE_ID, "operator_case_id": OPERATOR_CASE_ID,
        "status": "MASS_MATRIX_PASS_UNQUALIFIED" if passed else "MASS_MATRIX_BLOCKED",
        "classification": "NOT_MATERIAL_CALIBRATED", "solve_allowed": False,
        "input": {"path": str(input_path), "sha256": input_sha, "operator_source_pins_sha256": pins_sha, "mass_case_sha256": case_sha},
        "density_kg_per_m3": density, "density_role": "fixed toy diagnostic only; not calibrated material data",
        "dimensions": {"semantic_dofs": 162, "scalar_basis_rows": 54, "volume_points": 256},
        "integrated_volume_m3": volume, "analytic_cylinder_volume_m3": analytic_volume, "characteristic_mass_kg": mass_scale,
        "mass_matrix_unit": "kg", "mass_matrix_kg": mass.tolist(),
        "component_and_trace": {"component_order": list(COMPONENTS), "component_block_traces_kg": traces,
            "scalar_basis_trace_kg": scalar_trace, "full_vector_trace_kg": trace,
            "trace_vs_three_scalar_traces_relative_error": trace_rel, "block_relative_difference_max": block_rel,
            "cross_component_max_abs_kg": cross_max},
        "symmetry": {"frobenius_residual_kg": sym, "relative_frobenius_residual": sym_rel},
        "spectrum": {"scale_kg": mass_scale, "scale_definition": "rho * sum(saved volume_weights_m3)",
            "eigenvalues_M_over_rhoV": eig.tolist(), "eigenvalues_M_kg": (eig * mass_scale).tolist(),
            "relative_zero_tolerance": 1e-12, "applied_zero_tolerance": zero_tol, "expected_zero_modes": 0,
            "observed_near_zero_modes": zeros, "expected_positive_modes": 162, "observed_positive_modes": positives,
            "significant_negative_modes": negatives, "gate_passed": spectrum_ok,
            "ridge_pseudoinverse_or_truncation_used": False, "unmodified_scaled_cholesky_passed": chol_ok,
            "cholesky_failure": chol_error, "future_cholesky_expectation": "PASS" if chol_ok else "BLOCKED_OR_NOT_ATTEMPTED"},
        "rigid_mode_physical_kinetic_forms": {"order": ["tx", "ty", "tz", "rx", "ry", "rz"],
            "units": "translation/translation kg; translation/rotation kg*m; rotation/rotation kg*m^2",
            "direct_volume_quadrature_6x6": direct.tolist(), "cached_coefficient_form_6x6": from_coeff.tolist(),
            "translation_mass_diagonal_kg": translation.tolist(), "rhoV_kg": mass_scale,
            "rotation_inertia_diagonal_kg_m2": inertia.tolist(), "analytic_cylinder_inertia_kg_m2": inertia_expected.tolist(),
            "rotation_center_xyz_m": center.tolist(), "scaled_rigid_field_error": recon_err,
            "normalized_kinetic_form_error": rigid_err, "translation_relative_error": translation_err,
            "inertia_relative_error": inertia_err, "source_centroid_error_m": centroid_err},
        "checks": checks, "execution": {"elapsed_seconds": time.monotonic() - started,
            "cached_arrays_only": True, "provider_or_mls_called": False, "operator_reassembled": False,
            "solver_called": False, "cad_called": False},
        "limitations": ["rho=1 is a toy diagnostic, not material calibration",
            "mass consistency does not qualify stiffness, fracture, or a physical material",
            "any negative or near-zero mass mode blocks without tuning or matrix repair"]}
    return result, registry["ordered_canonical_ids"]


def _write(result: dict[str, object], ids: list[str], output: Path, csv_path: Path) -> None:
    if not _inside(output, REVIEW_ROOT) or not _inside(csv_path, REVIEW_ROOT): raise ReplayBlocked("output escapes review directory")
    if output.resolve() == csv_path.resolve() or output.exists() or csv_path.exists(): raise ReplayBlocked("refusing overwrite/retry")
    output.parent.mkdir(parents=True, exist_ok=True)
    tj, tc = output.with_name("." + output.name + ".tmp"), csv_path.with_name("." + csv_path.name + ".tmp")
    if tj.exists() or tc.exists(): raise ReplayBlocked("prior partial output exists; refusing retry")
    with tc.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream); writer.writerow(["canonical_id", *ids])
        for key, row in zip(ids, result["mass_matrix_kg"]): writer.writerow([key, *row])
    if tc.stat().st_size > MAX_BYTES:
        tc.unlink(missing_ok=True); raise ReplayBlocked("matrix CSV exceeds 8 MiB")
    data = json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    if len(data.encode("utf-8")) > MAX_BYTES:
        tc.unlink(missing_ok=True); raise ReplayBlocked("result JSON exceeds 8 MiB")
    tj.write_text(data, encoding="utf-8"); tc.replace(csv_path); tj.replace(output)


def main() -> int:
    started = time.monotonic()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expected-input-sha256", required=True)
    parser.add_argument("--expected-case-sha256", required=True)
    args = parser.parse_args()
    for value in (args.expected_input_sha256, args.expected_case_sha256):
        if not re.fullmatch(r"[0-9a-f]{64}", value): raise ReplayBlocked("SHA-256 pins must be 64 lowercase hex characters")
    case_sha = _sha256(CASE_PATH)
    if case_sha != args.expected_case_sha256: raise ReplayBlocked("fixed-case SHA-256 mismatch")
    case = _read_json(CASE_PATH, 1024 * 1024)
    if case.get("case_id") != CASE_ID or case.get("status") != "SOURCE_ONLY_WAITING_FOR_ROOT_AUTHORIZATION":
        raise ReplayBlocked("fixed-case identity/status changed")
    execution, spec = case.get("execution"), case.get("input")
    if not isinstance(execution, dict) or execution.get("authorized") is not True:
        raise ReplayBlocked("execution is disabled until root authorizes the pinned input")
    if not isinstance(spec, dict) or spec.get("sha256") != args.expected_input_sha256:
        raise ReplayBlocked("fixed-case operator result SHA-256 differs from supplied root pin")
    input_path = (ROOT / spec["path"]).resolve()
    if not _inside(input_path, ROOT) or not input_path.is_file() or input_path.stat().st_size > MAX_BYTES:
        raise ReplayBlocked("operator input is missing, outside workspace, or exceeds 8 MiB")
    input_sha = _sha256(input_path)
    if input_sha != args.expected_input_sha256: raise ReplayBlocked("operator result SHA-256 mismatch")
    expected, pins_sha = _pins(ROOT, ROOT / spec["source_pins_path"])
    replay = _preflight(_read_json(input_path, MAX_BYTES), expected)
    density = float(case["physical_model"]["density_kg_per_m3"])
    if not math.isfinite(density) or density != 1.0: raise ReplayBlocked("fixed toy density must remain rho=1 kg/m^3")
    result, ids = _calculate(replay, density, input_path, input_sha, pins_sha, case_sha, started)
    output, csv_path = (ROOT / p for p in execution["outputs"])
    _write(result, ids, output, csv_path)
    print(json.dumps({"status": result["status"], "result": str(output), "matrix_csv": str(csv_path),
                      "elapsed_seconds": result["execution"]["elapsed_seconds"]}, ensure_ascii=False))
    return 0 if result["status"] == "MASS_MATRIX_PASS_UNQUALIFIED" else 3


if __name__ == "__main__":
    try: raise SystemExit(main())
    except ReplayBlocked as exc:
        print(f"MASS_REPLAY_BLOCKED: {exc}", file=sys.stderr); raise SystemExit(2)
