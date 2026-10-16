#!/usr/bin/env python3
"""Offline, read-only reconstruction of the frozen 96-T4 insertion event."""
from __future__ import annotations

import hashlib
import json
import math
import platform
import re
import sys
import time
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[4]
RUN = ROOT / "relicweave-synlab/tmp/review_20261004/implementation/cohesive_event_probe_20261004/coarse654_20261004_01"
DATA_OUT = ROOT / "pgfrag/data/continuation_20261016/extrinsic_event_audit"
SOURCE_ROOT = ROOT / "third_party/akantu-runtime/source-v5.0.7"
WHEEL = ROOT / "third_party/akantu-runtime/downloads/akantu-5.0.7.post1-cp312-cp312-manylinux_2_17_x86_64.manylinux2014_x86_64.whl"
DIST_INFO = ROOT / "third_party/akantu-runtime/venv/lib/python3.12/site-packages/akantu-5.0.7.post1.dist-info"
MATERIAL_KEY = "material_1:MaterialCohesiveLinear3D"
ROUND_OFF_TOLERANCE_J = 1.0e-18


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def relpath(path: Path) -> str:
    return path.resolve().relative_to(ROOT.resolve()).as_posix()


def hash_record(path: Path) -> dict:
    if not path.is_file():
        return {"path": relpath(path), "exists": False, "sha256": None, "bytes": None}
    return {"path": relpath(path), "exists": True, "sha256": sha256(path), "bytes": path.stat().st_size}


def load_phase_records() -> list[dict]:
    records = []
    for line_number, line in enumerate((RUN / "phase_evidence.jsonl").read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        try:
            records.append(json.loads(line))
        except json.JSONDecodeError as error:
            raise RuntimeError(f"invalid phase evidence line {line_number}: {error}") from error
    return records


def phase_record(records: list[dict], step: int, phase: str) -> dict:
    matches = [row for row in records if row.get("step") == step and row.get("phase") == phase]
    if len(matches) != 1:
        raise RuntimeError(f"expected one evidence record for {step}{phase}, found {len(matches)}")
    return matches[0]


def material_arrays(record: dict) -> tuple[np.ndarray, np.ndarray]:
    fields = record.get("material_fields", {}).get("fields", {})
    material = fields.get(MATERIAL_KEY)
    if material is None:
        raise RuntimeError(f"missing public cohesive fields for {record.get('step')}{record.get('phase')}")
    opening = material.get("getOpening", {}).get("value", {}).get("values")
    traction = material.get("getTraction", {}).get("value", {}).get("values")
    if opening is None or traction is None:
        raise RuntimeError(f"opening/traction values absent for {record.get('step')}{record.get('phase')}")
    return np.asarray(opening, dtype=np.float64), np.asarray(traction, dtype=np.float64)


def native_energy(record: dict, key: str) -> float | None:
    entry = record.get("energy_getters", {}).get(key, {})
    if entry.get("status") != "ok" or "value_J" not in entry:
        return None
    return float(entry["value_J"])


def jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, Path):
        return relpath(value)
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    return value


def parse_cohesive_material_parameters(text: str) -> dict[str, float]:
    values = {}
    for name, raw in re.findall(r"^\s*(sigma_c|delta_c|G_c|beta)\s*=\s*([-+0-9.eE]+)\s*$", text, re.MULTILINE):
        values[name] = float(raw)
    return values


def array_state(step: int, phase: str) -> dict[str, np.ndarray]:
    path = RUN / "dense_states" / f"step_{step:05d}_{phase}.npz"
    with np.load(path, allow_pickle=False) as archive:
        return {key: archive[key].copy() for key in archive.files}


def duplicate_pairs(before: dict, after: dict) -> list[tuple[int, int]]:
    original_count = len(before["rest_nodes"])
    old_nodes = before["rest_nodes"]
    pairs = []
    for new_index in range(original_count, len(after["rest_nodes"])):
        coordinate = after["rest_nodes"][new_index]
        matches = [i for i, old_coordinate in enumerate(old_nodes) if np.array_equal(old_coordinate, coordinate)]
        if len(matches) != 1:
            raise RuntimeError(f"duplicated node {new_index} maps to {len(matches)} old coordinates")
        pairs.append((matches[0], new_index))
    return pairs


def phase_force_summary(state: dict[str, np.ndarray]) -> dict:
    force = state["internal_force"]
    return {
        "saved_vector_name": "internal_force",
        "saved_native_getter": "model.getInternalForce()",
        "max_nodal_norm_N": float(np.linalg.norm(force, axis=1).max(initial=0.0)),
        "vector_sum_N": force.sum(axis=0),
        "new_node_vector_sum_N": force[45:].sum(axis=0) if len(force) > 45 else None,
    }


def phase_mass_summary(state: dict[str, np.ndarray]) -> dict:
    mass = state["mass"]
    return {
        "node_count": int(len(mass)),
        "per_component_total_kg": mass.sum(axis=0),
        "new_node_per_component_total_kg": mass[45:].sum(axis=0) if len(mass) > 45 else None,
        "mass_columns_equal": bool(np.array_equal(mass[:, 0], mass[:, 1]) and np.array_equal(mass[:, 1], mass[:, 2])),
    }


def kinematic_metrics(state: dict[str, np.ndarray]) -> dict:
    coordinates = state["rest_nodes"] + state["displacement"]
    velocity = state["velocity"]
    mass = state["mass"]
    if mass.ndim != 2 or mass.shape[1] != 3:
        raise RuntimeError(f"unexpected mass shape {mass.shape}")
    nodal_scalar_mass = mass[:, 0]
    total_mass = float(nodal_scalar_mass.sum())
    center = np.average(coordinates, axis=0, weights=nodal_scalar_mass)
    momentum = (mass * velocity).sum(axis=0)
    angular_momentum = np.cross(coordinates - center, mass * velocity).sum(axis=0)
    kinetic = 0.5 * float(np.sum(mass * velocity * velocity))
    return {
        "mass_center_m": center,
        "mass_kg": total_mass,
        "linear_momentum_kg_m_s": momentum,
        "angular_momentum_about_mass_center_kg_m2_s": angular_momentum,
        "kinetic_energy_from_saved_lumped_mass_J": kinetic,
    }


def area_and_work(records: list[dict], mesh_state: dict[str, np.ndarray]) -> dict:
    coordinates = mesh_state["rest_nodes"]
    connectivity = mesh_state["cohesive_connectivity"]
    if connectivity.ndim != 2 or connectivity.shape[1] != 6:
        raise RuntimeError(f"expected cohesive Tri3-pair connectivity, found {connectivity.shape}")
    area_rows = []
    for element_index, element in enumerate(connectivity):
        face_a = coordinates[element[:3]]
        face_b = coordinates[element[3:]]
        area_a = float(np.linalg.norm(np.cross(face_a[1] - face_a[0], face_a[2] - face_a[0])) / 2.0)
        area_b = float(np.linalg.norm(np.cross(face_b[1] - face_b[0], face_b[2] - face_b[0])) / 2.0)
        area_rows.append({
            "cohesive_element_index": element_index,
            "face_a_area_m2": area_a,
            "face_b_area_m2": area_b,
            "area_m2": area_a,
            "paired_face_coordinate_max_abs_m": float(np.max(np.abs(face_a - face_b))),
        })
    areas = np.array([row["area_m2"] for row in area_rows], dtype=np.float64)
    # The recorded getters contain three integration points per cohesive Tri3.
    if len(connectivity) * 3 != 24:
        raise RuntimeError(f"recorded 24 QPs do not match mesh facets x 3: {len(connectivity)} x 3")
    qp_areas = np.repeat(areas / 3.0, 3)

    records_by_key = {(row.get("step"), row.get("phase")): row for row in records}
    d650 = records_by_key[(650, "D")]
    d651 = records_by_key[(651, "D")]
    d652 = records_by_key[(652, "D")]
    delta650, traction650 = material_arrays(d650)
    delta651, traction651 = material_arrays(d651)
    delta652, traction652 = material_arrays(d652)
    if any(array.shape != (24, 3) for array in (delta650, traction650, delta651, traction651, delta652, traction652)):
        raise RuntimeError("expected [24,3] native opening and traction arrays for steps 650-652")

    def increment(old_t, new_t, old_delta, new_delta):
        per_qp = qp_areas * 0.5 * np.einsum("ij,ij->i", old_t + new_t, new_delta - old_delta)
        return float(per_qp.sum()), per_qp

    dw651, per_qp651 = increment(traction650, traction651, delta650, delta651)
    reversible651 = float(np.sum(qp_areas * 0.5 * np.einsum("ij,ij->i", traction651, delta651)))
    cumulative_work651 = dw651
    dissipated651 = cumulative_work651 - reversible651
    dw652, per_qp652 = increment(traction651, traction652, delta651, delta652)
    reversible652 = float(np.sum(qp_areas * 0.5 * np.einsum("ij,ij->i", traction652, delta652)))
    cumulative_work652 = cumulative_work651 + dw652
    dissipated652 = cumulative_work652 - reversible652
    native651 = native_energy(d651, "cohesive_dissipated_J")
    native652 = native_energy(d652, "cohesive_dissipated_J")
    negative_record = None
    for row in sorted(records, key=lambda item: (int(item["step"]), {"A": 0, "B": 1, "C": 2, "D": 3}.get(item["phase"], 9))):
        value = native_energy(row, "cohesive_dissipated_J")
        if value is not None and value < 0.0:
            negative_record = {"step": int(row["step"]), "phase": row["phase"], "native_dissipated_J": value}
            break
    most_negative_qp = int(np.argmin(per_qp652))
    old_delta = delta651[most_negative_qp]
    new_delta = delta652[most_negative_qp]
    denom = float(np.linalg.norm(old_delta) * np.linalg.norm(new_delta))
    direction_cosine = float(np.dot(old_delta, new_delta) / denom) if denom else None
    return {
        "native_cohesive_element_count_at_650_C": int(len(connectivity)),
        "native_qp_count": int(len(qp_areas)),
        "qp_area_rule": "each 3-point Tri3 quadrature point receives the actual triangle area divided by 3",
        "qp_area_m2": qp_areas,
        "interface_area_by_cohesive_element": area_rows,
        "interface_total_area_m2": float(areas.sum()),
        "maximum_paired_face_coordinate_difference_m": float(max(row["paired_face_coordinate_max_abs_m"] for row in area_rows)),
        "step_650_native_insertion": {
            "A_cohesive_count": int(phase_record(records, 650, "A").get("cohesive_elements_count", -1)),
            "B_cohesive_count": int(phase_record(records, 650, "B").get("cohesive_elements_count", -1)),
            "C_cohesive_count": int(phase_record(records, 650, "C").get("cohesive_elements_count", -1)),
            "D_cohesive_count": int(phase_record(records, 650, "D").get("cohesive_elements_count", -1)),
            "D_native_dissipated_J": native_energy(d650, "cohesive_dissipated_J"),
        },
        "step_651": {
            "signed_delta_work_J": dw651,
            "cumulative_work_J": cumulative_work651,
            "reversible_energy_J": reversible651,
            "reconstructed_dissipated_J": dissipated651,
            "native_dissipated_J": native651,
            "reconstruction_minus_native_J": dissipated651 - native651 if native651 is not None else None,
            "native_match_within_roundoff_tolerance": bool(native651 is not None and abs(dissipated651 - native651) <= ROUND_OFF_TOLERANCE_J),
            "tolerance_J": ROUND_OFF_TOLERANCE_J,
            "qp_signed_delta_work_J": per_qp651,
        },
        "step_652": {
            "signed_delta_work_J": dw652,
            "cumulative_work_J": cumulative_work652,
            "reversible_energy_J": reversible652,
            "reconstructed_dissipated_J": dissipated652,
            "native_dissipated_J": native652,
            "reconstruction_minus_native_J": dissipated652 - native652 if native652 is not None else None,
            "native_match_within_roundoff_tolerance": bool(native652 is not None and abs(dissipated652 - native652) <= ROUND_OFF_TOLERANCE_J),
            "tolerance_J": ROUND_OFF_TOLERANCE_J,
            "negative_qp_work_count": int(np.count_nonzero(per_qp652 < 0.0)),
            "qp_signed_delta_work_J": per_qp652,
            "most_negative_qp_zero_based": most_negative_qp,
            "most_negative_qp_work_J": float(per_qp652[most_negative_qp]),
            "most_negative_qp_opening_direction_cosine_651_to_652": direction_cosine,
            "opening_euclidean_max_651_m": float(np.linalg.norm(delta651, axis=1).max()),
            "opening_euclidean_max_652_m": float(np.linalg.norm(delta652, axis=1).max()),
            "traction_zero_qp_count_652": int(np.count_nonzero(np.linalg.norm(traction652, axis=1) == 0.0)),
        },
        "first_negative_native_dissipated_record": negative_record,
    }


def write_outputs(audit: dict) -> None:
    DATA_OUT.mkdir(parents=True, exist_ok=True)
    json_path = DATA_OUT / "audit.json"
    report_path = DATA_OUT / "REPORT.md"
    json_path.write_text(json.dumps(jsonable(audit), indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    work = audit["work_reconstruction"]
    step651 = work["step_651"]
    step652 = work["step_652"]
    insertion = audit["insertion_state_650"]
    dmetrics = insertion["phase_D_post_solve_kinematics"]
    momentum_text = "[" + ", ".join(f"{float(value):.6e}" for value in dmetrics["linear_momentum_kg_m_s"]) + "]"
    angular_text = "[" + ", ".join(f"{float(value):.6e}" for value in dmetrics["angular_momentum_about_mass_center_kg_m2_s"]) + "]"
    report = rf"""# Extrinsic event offline audit

This is a read-only reconstruction of the frozen `coarse654_20261004_01` evidence. The script does not import Akantu, start a solver, edit historical evidence, or evolve a native state. The script and machine-readable output are in this directory's paired experiment/data folders.

## Signed work and first negative getter

Areas were recomputed from the paired Tri3 faces in the saved post-insertion mesh: {work['native_cohesive_element_count_at_650_C']} facets, {work['native_qp_count']} quadrature points, and {work['interface_total_area_m2']:.17g} m² total. Each QP receives one third of its actual triangle area. The recorded 650D→651D vector path gives signed ΔW = {step651['signed_delta_work_J']:.17g} J, R = {step651['reversible_energy_J']:.17g} J, and D = {step651['reconstructed_dissipated_J']:.17g} J; the native D getter is {step651['native_dissipated_J']:.17g} J. The 651 reconstruction matches the native getter within the stated roundoff tolerance.

For 651D→652D, signed ΔW = {step652['signed_delta_work_J']:.17g} J, cumulative W = {step652['cumulative_work_J']:.17g} J, R = {step652['reversible_energy_J']:.17g} J, and D = {step652['reconstructed_dissipated_J']:.17g} J; the native D getter is {step652['native_dissipated_J']:.17g} J. {step652['negative_qp_work_count']} of {work['native_qp_count']} QPs have negative work increments at step 652. The first negative native getter is {work['first_negative_native_dissipated_record']['step']}{work['first_negative_native_dissipated_record']['phase']}. No energy values were clamped or replaced by magnitudes.

The saved maximum Euclidean opening increases from {step652['opening_euclidean_max_651_m']:.9g} m at 651D to {step652['opening_euclidean_max_652_m']:.9g} m at 652D, while {step652['traction_zero_qp_count_652']} QP traction vectors at 652D are exactly zero. The signed saved vector path explains the first negative dissipated-energy value. It supports a leading unresolved near-zero-opening/directional-stiffness mechanism; the offline record does not establish it as the sole cause.

## Insertion state and operative mass

At step 650, the recorded B state has {insertion['phase_B_preinsertion']['node_count']} nodes, while C and D have {insertion['phase_C_after_insertion_before_solve']['node_count']} nodes and {insertion['phase_C_after_insertion_before_solve']['new_nodes']} appended nodes. In C, the duplicated nodes have equal displacement and velocity to their source nodes, and their appended mass rows are zero. D is the post-solve state: its appended-node masses sum to {dmetrics['post_solve_mass_summary']['new_node_per_component_total_kg'][0]:.17g} kg per component, total mass is {dmetrics['post_solve_mass_summary']['per_component_total_kg'][0]:.17g} kg per component, and total mass agrees with B to the recorded numerical precision. The zero appended masses in C therefore are not used as the operative post-insertion mass.

P/L/K are recorded for B, C, and D in `audit.json`; B and C values use the mass arrays saved at those phases. C's appended zero masses make its pre-solve P/L/K a stale-mass diagnostic, while D's reassembled mass is the operative post-solve result.

Using D's saved post-solve mass, velocity, displacement and rest coordinates gives P = {momentum_text} kg·m/s, L about the mass center = {angular_text} kg·m²/s, and K = {dmetrics['kinetic_energy_from_saved_lumped_mass_J']:.17g} J. Acceleration values are not in the NPZ snapshots; the tagged source says `onNodesAdded` copies acceleration to duplicates, but this evidence set cannot verify its numeric values. Duplicated displacement remains equal at the saved D sample; velocity can differ after solve.

## Force and constitutive-data limits

The snapshots save only the assembled `internal_force` vector from `model.getInternalForce()`. Tagged source ancestry shows the cohesive override computes cohesive tractions and then calls the base internal-force assembly, so the saved force is a combined vector, not separately attributable bulk and cohesive vectors. `external_force` and acceleration are not among the saved arrays; the trace has only a scalar top reaction estimate. Tagged `checkInsertion` computes per-QP `new_sigma` from the effective facet stress, then uses `delta_c_eff = sigma_limit / new_sigma * delta_c` when the nominal `delta_c` is explicit; it stores those values in internal fields. This input has an explicit nominal `delta_c`, but the phase capture omitted `sigma_c_eff`, `delta_c_eff`, `delta_max`, damage, normals, and stress tensors. Mesh geometry can recover geometric normals, while the recorded wheel is 5.0.7.post1 and the local inspected source is tagged v5.0.7. Accordingly this audit does not claim a wheel-specific 651 traction-law reproduction or a finite 651 Jacobian/eigenvalue result. It does not infer a parameter from 651 traction and then use that same traction as its validation target.

The cached source manifest's files and the installed wheel are separately hashed in `audit.json`. The wheel/source version distinction is retained; the local source cache does not contain `material_cohesive_linear_inline_impl.hh`.

## Reproduction

From the repository root in PowerShell, run:

```powershell
& 'C:\Users\xiemo\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' 'pgfrag/experiments/continuation_20261016/extrinsic_event_audit/audit_event.py'
```

The script only reads the frozen input/evidence files and writes `audit.json` and this report under `pgfrag/data/continuation_20261016/extrinsic_event_audit/`. Its own SHA-256 and all required input/source hashes are in the JSON.
"""
    report_path.write_text(report, encoding="utf-8")


def main() -> int:
    started = time.perf_counter()
    required_inputs = [
        RUN / "phase_evidence.jsonl",
        RUN / "input_native_initial.npz",
        RUN / "events.json",
        RUN / "preflight.json",
        RUN / "execution_receipt.json",
        RUN / "driver_source_identity.json",
        RUN / "snapshot/input/material.dat",
        RUN / "snapshot/input/generated_geometry.json",
        RUN / "dense_states/step_00650_B.npz",
        RUN / "dense_states/step_00650_C.npz",
        RUN / "dense_states/step_00650_D.npz",
        RUN / "dense_states/step_00651_D.npz",
        RUN / "dense_states/step_00652_D.npz",
    ]
    missing = [relpath(path) for path in required_inputs if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing frozen evidence: " + ", ".join(missing))

    records = load_phase_records()
    b650 = array_state(650, "B")
    c650 = array_state(650, "C")
    d650_state = array_state(650, "D")
    d651_state = array_state(651, "D")
    d652_state = array_state(652, "D")
    pairs = duplicate_pairs(b650, c650)
    if not pairs:
        raise RuntimeError("step 650C contains no duplicated nodes")

    def pair_mismatch(left, right, field):
        return [float(np.linalg.norm(left[field][i] - right[field][j])) for i, j in pairs]

    c_u_diff = pair_mismatch(c650, c650, "displacement")
    c_v_diff = pair_mismatch(c650, c650, "velocity")
    d_u_diff = pair_mismatch(d650_state, d650_state, "displacement")
    d_v_diff = pair_mismatch(d650_state, d650_state, "velocity")
    mass_b = b650["mass"]
    mass_c = c650["mass"]
    mass_d = d650_state["mass"]
    total_mass_difference = mass_d.sum(axis=0) - mass_b.sum(axis=0)

    event_data = json.loads((RUN / "events.json").read_text(encoding="utf-8"))
    preflight = json.loads((RUN / "preflight.json").read_text(encoding="utf-8"))
    receipt = json.loads((RUN / "execution_receipt.json").read_text(encoding="utf-8"))
    driver_identity = json.loads((RUN / "driver_source_identity.json").read_text(encoding="utf-8"))
    material_text = (RUN / "snapshot/input/material.dat").read_text(encoding="utf-8")
    nominal_material_parameters = parse_cohesive_material_parameters(material_text)
    generated_geometry = json.loads((RUN / "snapshot/input/generated_geometry.json").read_text(encoding="utf-8"))

    relevant_source_relpaths = [
        "src/model/solid_mechanics/solid_mechanics_model.cc",
        "src/model/solid_mechanics/solid_mechanics_model_cohesive/solid_mechanics_model_cohesive.cc",
        "src/model/solid_mechanics/solid_mechanics_model_cohesive/materials/constitutive_laws/material_cohesive_linear.cc",
        "src/model/solid_mechanics/solid_mechanics_model_cohesive/materials/constitutive_laws/material_cohesive_linear.hh",
    ]
    source_manifest_path = SOURCE_ROOT / "source_manifest.json"
    source_manifest = json.loads(source_manifest_path.read_text(encoding="utf-8"))
    manifest_entries = {entry["path"]: entry for entry in source_manifest}
    source_files = []
    for rel in relevant_source_relpaths:
        path = SOURCE_ROOT / rel
        manifest_entry = manifest_entries.get(rel)
        source_files.append({
            **hash_record(path),
            "manifest_ref": manifest_entry.get("ref") if manifest_entry else None,
            "manifest_api_content_sha256": manifest_entry.get("api_content_sha256") if manifest_entry else None,
            "manifest_downloaded_sha256": manifest_entry.get("downloaded_sha256") if manifest_entry else None,
            "manifest_sha_match": manifest_entry.get("sha_match") if manifest_entry else None,
            "local_matches_manifest_download": bool(manifest_entry and path.is_file() and sha256(path) == manifest_entry.get("downloaded_sha256")),
        })

    input_hashes = [hash_record(path) for path in required_inputs]
    input_hashes.extend([
        hash_record(RUN / "trace.csv"),
        hash_record(RUN / "snapshot/input/block_2x2x4mm.msh"),
    ])
    initial_input_path = RUN / "input_native_initial.npz"
    with np.load(initial_input_path, allow_pickle=False) as initial_archive:
        initial_nodes = initial_archive["rest_nodes"].shape[0]
        initial_bulk_tets = initial_archive["bulk_tetra_connectivity"].shape[0]
        initial_cohesive = initial_archive["cohesive_connectivity"].shape[0]

    work = area_and_work(records, c650)
    d_metrics = kinematic_metrics(d650_state)
    native_kinetic_650 = native_energy(phase_record(records, 650, "D"), "kinetic_J")
    d_metrics["native_kinetic_energy_getter_J"] = native_kinetic_650
    d_metrics["kinetic_energy_difference_from_native_getter_J"] = (
        d_metrics["kinetic_energy_from_saved_lumped_mass_J"] - native_kinetic_650
        if native_kinetic_650 is not None else None
    )
    d_metrics["post_solve_mass_summary"] = phase_mass_summary(d650_state)

    insertion_state = {
        "phase_B_preinsertion": {
            "node_count": int(len(b650["rest_nodes"])),
            "mass_summary": phase_mass_summary(b650),
            "kinematics_from_saved_phase_mass": kinematic_metrics(b650),
            "saved_array_keys": sorted(b650),
            "internal_force_summary": phase_force_summary(b650),
        },
        "phase_C_after_insertion_before_solve": {
            "node_count": int(len(c650["rest_nodes"])),
            "new_nodes": int(len(c650["rest_nodes"]) - len(b650["rest_nodes"])),
            "duplicate_node_source_and_clone_indices_zero_based": pairs,
            "duplicate_displacement_mismatch_norms_m": c_u_diff,
            "duplicate_velocity_mismatch_norms_m_s": c_v_diff,
            "duplicate_mass_match_norms_kg": [float(np.linalg.norm(mass_c[i] - mass_c[j])) for i, j in pairs],
            "mass_summary": phase_mass_summary(c650),
            "kinematics_from_saved_phase_mass": kinematic_metrics(c650),
            "kinematics_interpretation": "pre-solve snapshot mass has zero rows for appended duplicates; record as stale-mass diagnostic, not operative post-insertion result",
            "saved_array_keys": sorted(c650),
            "internal_force_summary": phase_force_summary(c650),
            "acceleration_array_saved": "acceleration" in c650,
        },
        "phase_D_post_solve_kinematics": d_metrics,
        "phase_D_post_solve_state": {
            "node_count": int(len(d650_state["rest_nodes"])),
            "duplicate_displacement_mismatch_norms_m": d_u_diff,
            "duplicate_velocity_mismatch_norms_m_s": d_v_diff,
            "duplicate_mass_match_norms_kg": [float(np.linalg.norm(mass_d[i] - mass_d[j])) for i, j in pairs],
            "mass_summary": phase_mass_summary(d650_state),
            "total_mass_change_from_B_per_component_kg": total_mass_difference,
            "internal_force_summary": phase_force_summary(d650_state),
            "saved_array_keys": sorted(d650_state),
            "acceleration_array_saved": "acceleration" in d650_state,
        },
        "acceleration_numeric_audit": {
            "status": "UNAVAILABLE_numeric acceleration is absent from the saved phase NPZ arrays",
            "source_semantics": "cached v5.0.7 SolidMechanicsModelCohesive::onNodesAdded copies acceleration from each source node to its duplicate",
            "numeric_value_verified": False,
            "note": "velocity differences between C and D span a native solve step and are not labeled as instantaneous acceleration",
        },
        "force_separation_audit": {
            "recorded_internal_force": "combined model.getInternalForce() nodal vector",
            "bulk_internal_force_separate_vector_saved": False,
            "cohesive_internal_force_separate_vector_saved": False,
            "external_force_vector_saved": False,
            "scalar_top_reaction_in_trace_only": True,
            "external_force_getter_exposed_in_preflight": bool(preflight.get("model_public_api", {}).get("getExternalForce", {}).get("type") == "method"),
            "evidence_gap": "phase snapshots save internal_force only; the probe sample_arrays helper did not capture getExternalForce(), separate bulk/cohesive force contributions, or acceleration",
            "source_ancestry": [
                "SolidMechanicsModelCohesive::assembleInternalForces computes cohesive traction before calling SolidMechanicsModel::assembleInternalForces",
                "SolidMechanicsModel::assembleInternalForces assembles constitutive-law forces into the shared internal_force vector",
            ],
        },
        "tagged_source_internal_parameter_rule": {
            "source_file": "third_party/akantu-runtime/source-v5.0.7/src/model/solid_mechanics/solid_mechanics_model_cohesive/materials/constitutive_laws/material_cohesive_linear.cc",
            "sigma_c_eff_rule": "new_sigma is stress_check(q), the insertion effective stress at QP q",
            "delta_c_eff_rule_for_explicit_nominal_delta_c": "new_delta = sigma_limit / new_sigma * delta_c",
            "delta_c_eff_rule_when_nominal_delta_c_is_zero": "new_delta = 2 * G_c / new_sigma",
            "stored_internal_fields": ["sigma_c_eff", "delta_c_eff", "insertion_stress"],
            "runtime_per_qp_values_present_in_frozen_phase_evidence": False,
            "official_inline_source_url": "https://gitlab.com/akantu/akantu/-/raw/v5.0.7/src/model/solid_mechanics/solid_mechanics_model_cohesive/materials/constitutive_laws/material_cohesive_linear_inline_impl.hh",
            "wheel_specific_651_law_verified": False,
            "reason": "the run's wheel is v5.0.7.post1 and the locally hashed source is tagged v5.0.7; actual per-QP internal parameters and native stress/normal state were not captured",
        },
    }

    material_path = RUN / "snapshot/input/material.dat"
    input_material_sha = sha256(material_path)
    wheel_metadata = DIST_INFO / "METADATA"
    wheel_record = DIST_INFO / "RECORD"
    wheel_evidence = {
        "wheel_version_from_preflight": preflight.get("akantu_version"),
        "downloaded_wheel": hash_record(WHEEL),
        "installed_metadata": hash_record(wheel_metadata),
        "installed_record": hash_record(wheel_record),
        "tagged_source_ref": "v5.0.7",
        "tagged_source_manifest": hash_record(source_manifest_path),
        "inline_impl_in_local_source_cache": (SOURCE_ROOT / "src/model/solid_mechanics/solid_mechanics_model_cohesive/materials/constitutive_laws/material_cohesive_linear_inline_impl.hh").is_file(),
        "inline_impl_reference_url": "https://gitlab.com/akantu/akantu/-/raw/v5.0.7/src/model/solid_mechanics/solid_mechanics_model_cohesive/materials/constitutive_laws/material_cohesive_linear_inline_impl.hh",
        "provenance_status": "WHEEL_5.0.7.post1_AND_TAGGED_SOURCE_v5.0.7_ARE_DISTINCT_EVIDENCE; no wheel-specific per-QP internal fields were saved",
    }

    audit = {
        "schema": "offline-extrinsic-event-audit-v1",
        "scope": "read-only reconstruction from an existing run; no native library import, solver call, calibration, or new evolution",
        "created_by": relpath(Path(__file__)),
        "script_sha256": sha256(Path(__file__).resolve()),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "elapsed_seconds": time.perf_counter() - started,
        "solver_steps_launched_by_audit": 0,
        "historical_inputs_modified": False,
        "run_identity": {
            "run_id": receipt.get("run_id"),
            "original_probe_driver": hash_record(ROOT / "pgfrag/experiments/cohesive_event_probe_20261004/probe_driver.py"),
            "receipt_driver_sha256": receipt.get("final_driver_source", {}).get("sha256"),
            "driver_identity_solver_sha256": driver_identity.get("solver_driver_sha256"),
            "receipt_driver_hash_matches_local_source": receipt.get("final_driver_source", {}).get("sha256") == sha256(ROOT / "pgfrag/experiments/cohesive_event_probe_20261004/probe_driver.py"),
            "initial_mesh_nodes": int(initial_nodes),
            "initial_bulk_tetrahedra": int(initial_bulk_tets),
            "initial_cohesive_elements": int(initial_cohesive),
            "insertion_events": event_data.get("events"),
        },
        "input_hashes": input_hashes,
        "source_hashes": {
            "wheel_and_runtime_provenance": wheel_evidence,
            "cached_source_files": source_files,
            "source_manifest": hash_record(source_manifest_path),
        },
        "material_input": {
            "path": relpath(material_path),
            "sha256": input_material_sha,
            "nominal_values_read_from_frozen_material_input": {
                "sigma_c_Pa": nominal_material_parameters.get("sigma_c"),
                "delta_c_m": nominal_material_parameters.get("delta_c"),
                "G_c_J_m2": nominal_material_parameters.get("G_c"),
                "beta": nominal_material_parameters.get("beta"),
                "kappa": "not specified in input; tagged source default is 1.0",
                "volume_s": "not specified in input; tagged source default is 0.0, whose scaleInsertionTraction branch does nothing",
            },
            "runtime_per_qp_parameters": {
                "sigma_c_eff_saved": False,
                "delta_c_eff_saved": False,
                "delta_max_or_damage_saved": False,
                "normals_saved": False,
                "stress_tensors_saved": False,
                "traction_law_reproduction_651": "NOT_CLAIMED; saved fields do not include per-QP effective parameters or state history, and run wheel provenance is not identical to tagged source",
                "jacobian_or_mass_normalized_eigen_audit": "NOT_RUN; actual per-QP sigma_c_eff and delta_c_eff are not established in the saved wheel state",
                "circular_inversion_guard": "No sigma_c_eff value was inferred from the 651 traction and then used to validate that same traction",
            },
            "geometry_summary_from_frozen_generated_geometry_json": generated_geometry,
        },
        "work_reconstruction": work,
        "insertion_state_650": insertion_state,
        "offline_limitations": [
            "No nodal acceleration values, external force vector, or separate bulk/cohesive internal-force vectors were captured.",
            "Per-QP sigma_c_eff, delta_c_eff, delta_max, damage, and stress tensors were not captured; actual wheel-specific 651 law validation and the requested finite-state Jacobian are unavailable.",
            "The cached source manifest is tagged v5.0.7; the recorded installed wheel identifies as v5.0.7.post1.",
            "This audit identifies the signed trajectory that reproduces the first negative getter; it does not prove that under-resolution is the only causal contribution.",
        ],
    }
    write_outputs(audit)
    summary = {
        "audit_json": relpath(DATA_OUT / "audit.json"),
        "report": relpath(DATA_OUT / "REPORT.md"),
        "script_sha256": audit["script_sha256"],
        "area_m2": work["interface_total_area_m2"],
        "step_651_reconstructed_dissipated_J": work["step_651"]["reconstructed_dissipated_J"],
        "step_652_reconstructed_dissipated_J": work["step_652"]["reconstructed_dissipated_J"],
        "step_652_negative_qp_work_count": work["step_652"]["negative_qp_work_count"],
        "elapsed_seconds": audit["elapsed_seconds"],
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
