"""Single-use, step-capped phase probe cloned from the frozen cohesive contract."""
from __future__ import annotations

import csv
from datetime import datetime, timezone
import difflib
import hashlib
import inspect
import json
import math
import os
from pathlib import Path
import platform
import re
import sys
import time
import traceback

import numpy as np


WORK = Path(__file__).resolve().parent
REPO = WORK.parents[2]
RUN_ID = "coarse654_20261004_01"
RUN_DIR = REPO / "relicweave-synlab/tmp/review_20261004/implementation/cohesive_event_probe_20261004" / RUN_ID
SNAPSHOT = RUN_DIR / "snapshot"
INPUT = SNAPSHOT / "input"
MANIFEST_PATH = SNAPSHOT / "manifest.json"
REFERENCE = RUN_DIR / "reference"
ORIGINAL_SCRIPT = WORK / "run_contract.original.py"
FROZEN_SCRIPT = REPO / "pgfrag/experiments/cohesive_capability/run_contract.py"
FROZEN_MANIFEST = REPO / "relicweave-synlab/tmp/review_20261004/implementation/cohesive_capability/manifest.json"
FROZEN_CASE = REPO / "relicweave-synlab/tmp/review_20261004/implementation/cohesive_capability/modeI_dt"
MAX_STEPS = 654
MAX_WALL_S = 300.0
FOCUS_STEPS = set(range(647, 655))
EXPECTED_RUNNER_SHA256 = "67dd1eb98acc5f5369b9cbf1a9bfbd7f6681068172e537768e788274cb10652d"
EXPECTED_MANIFEST_SHA256 = "21b563069550e2091474fe94877bd3ceaae094ac4c92aed3900ecdd94702dd68"
ENERGY_REQUESTS = {
    "kinetic_J": "kinetic",
    "bulk_potential_J": "potential",
    "cohesive_dissipated_J": "dissipated",
    "cohesive_reversible_J": "reversible",
    "cohesive_contact_native_key_J": "cohesive contact",
    "contact_alias_diagnostic_J": "contact",
}
ENERGY_FIELDS_FOR_MATCH = (
    "kinetic_J", "bulk_potential_J", "cohesive_dissipated_J", "cohesive_reversible_J"
)
MATERIAL_FIELD_RE = re.compile(r"traction|opening|separation|damage|history|maximum|internal.?variable|cohesive|stress", re.I)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        value = float(value)
    if isinstance(value, float) and not math.isfinite(value):
        return repr(value)
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    return repr(value)


def write_json(path: Path, value) -> None:
    path.write_text(json.dumps(jsonable(value), indent=2, sort_keys=True) + "\n", encoding="utf-8")


def verify_and_record_snapshots() -> dict:
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    items = []
    pairs = [(FROZEN_SCRIPT, SNAPSHOT / "run_contract.original.py"),
             (FROZEN_MANIFEST, MANIFEST_PATH)]
    for name in ("block_2x2x4mm.msh", "material.dat", "generated_geometry.json"):
        pairs.append((FROZEN_SCRIPT.parent / "input" / name, INPUT / name))
    pairs.extend([
        (FROZEN_CASE / "trace.csv", REFERENCE / "trace.csv"),
        (FROZEN_CASE / "status.json", REFERENCE / "status.json"),
    ])
    for original, snapshot in pairs:
        items.append({
            "original_path": str(original), "snapshot_path": str(snapshot),
            "original_exists": original.is_file(), "snapshot_exists": snapshot.is_file(),
            "original_bytes": original.stat().st_size if original.is_file() else None,
            "snapshot_bytes": snapshot.stat().st_size if snapshot.is_file() else None,
            "original_sha256": sha256(original) if original.is_file() else None,
            "snapshot_sha256": sha256(snapshot) if snapshot.is_file() else None,
            "identical": original.is_file() and snapshot.is_file() and sha256(original) == sha256(snapshot),
        })
    good = all(item["identical"] for item in items)
    original_sha = sha256(FROZEN_SCRIPT) if FROZEN_SCRIPT.is_file() else None
    manifest_sha = sha256(MANIFEST_PATH) if MANIFEST_PATH.is_file() else None
    integrity = {
        "run_id": RUN_ID, "captured_utc_before_solver": datetime.now(timezone.utc).isoformat(),
        "source_sha256_expected": EXPECTED_RUNNER_SHA256, "source_sha256_actual": original_sha,
        "manifest_sha256_expected": EXPECTED_MANIFEST_SHA256, "manifest_sha256_actual": manifest_sha,
        "manifest_source_expected_runner_sha256": manifest.get("source_sha256", {}).get("run_contract.py"),
        "snapshot_pairs_identical": good, "snapshot_pairs": items,
        "solver_attempts": 0, "solver_cap_steps": MAX_STEPS,
    }
    if original_sha != EXPECTED_RUNNER_SHA256 or manifest_sha != EXPECTED_MANIFEST_SHA256:
        raise RuntimeError("frozen runner/manifest SHA mismatch; solver not started")
    if manifest.get("source_sha256", {}).get("run_contract.py") != EXPECTED_RUNNER_SHA256:
        raise RuntimeError("snapshot manifest source SHA mismatch; solver not started")
    if not good:
        raise RuntimeError("source/input/manifest/reference snapshot mismatch; solver not started")
    for rel, expected in manifest["input_sha256"].items():
        candidate = INPUT.parent / rel
        if not candidate.is_file() or sha256(candidate).lower() != expected.lower():
            raise RuntimeError(f"snapshot input hash mismatch {rel}; solver not started")
    baseline_status = json.loads((REFERENCE / "status.json").read_text(encoding="utf-8"))
    if sha256(REFERENCE / "trace.csv") != baseline_status.get("trace_csv_sha256"):
        raise RuntimeError("baseline reference trace differs from frozen status; solver not started")
    write_json(RUN_DIR / "pre_run_snapshot_integrity.json", integrity)
    return integrity


def make_source_diff() -> None:
    old = ORIGINAL_SCRIPT.read_text(encoding="utf-8").splitlines(keepends=True)
    new = Path(__file__).read_text(encoding="utf-8").splitlines(keepends=True)
    diff = "".join(difflib.unified_diff(old, new, fromfile="run_contract.original.py", tofile="probe_driver.py"))
    (WORK / "run_contract_instrumented.diff").write_text(diff, encoding="utf-8")
    (RUN_DIR / "snapshot" / "run_contract_instrumented.diff").write_text(diff, encoding="utf-8")


def inspect_public_methods(obj, field_only=False) -> dict:
    result = {}
    for name in sorted(n for n in dir(obj) if not n.startswith("_")):
        if field_only and not MATERIAL_FIELD_RE.search(name):
            continue
        try:
            member = getattr(obj, name)
        except Exception as error:
            result[name] = {"status": "attribute_error", "error": repr(error)}
            continue
        if not callable(member):
            if field_only and MATERIAL_FIELD_RE.search(name):
                result[name] = {"status": "public_property", "value": jsonable(member)}
            continue
        doc = inspect.getdoc(member)
        entry = {"status": "public_method", "doc": doc[:1200] if doc else None}
        if field_only:
            try:
                signature = inspect.signature(member)
                entry["signature"] = str(signature)
                required = [p for p in signature.parameters.values()
                            if p.default is inspect.Parameter.empty
                            and p.kind in (inspect.Parameter.POSITIONAL_ONLY, inspect.Parameter.POSITIONAL_OR_KEYWORD)]
                if not required:
                    try:
                        entry["noarg_value"] = jsonable(member())
                        entry["noarg_call_status"] = "ok"
                    except Exception as error:
                        entry["noarg_call_status"] = "error"
                        entry["noarg_call_error"] = repr(error)
                else:
                    entry["noarg_call_status"] = "requires_arguments"
            except Exception as error:
                entry["signature_status"] = "unknown"
                entry["signature_error"] = repr(error)
                entry["noarg_call_status"] = "not_called_signature_unknown"
        result[name] = entry
    return result


def energy_values(model) -> dict:
    values = {}
    for label, key in ENERGY_REQUESTS.items():
        try:
            value = float(model.getEnergy(key))
            values[label] = {"status": "ok", "requested_key": key, "value_J": value}
        except Exception as error:
            values[label] = {"status": "unknown_or_unsupported", "requested_key": key, "error": repr(error)}
    return values


def material_handles(model) -> tuple[list, dict]:
    access = {"nb_material_method": "UNKNOWN", "material_getter": "UNKNOWN", "errors": []}
    count = None
    nb = getattr(model, "getNbMaterial", None)
    if callable(nb):
        access["nb_material_method"] = "getNbMaterial"
        try:
            count = int(nb())
            access["material_count"] = count
        except Exception as error:
            access["errors"].append({"method": "getNbMaterial", "error": repr(error)})
    for name in ("getMaterial", "getMaterialByName"):
        getter = getattr(model, name, None)
        if not callable(getter):
            continue
        access["material_getter"] = name
        handles = []
        if name == "getMaterialByName":
            access["errors"].append({"method": name, "status": "requires known material name"})
            continue
        for index in range(max(1, count or 1)):
            try:
                handles.append(getter(index))
            except Exception as error:
                access["errors"].append({"method": name, "index": index, "error": repr(error)})
                break
        if handles:
            return handles, access
    getter = getattr(model, "getMaterials", None)
    if callable(getter):
        access["material_getter"] = "getMaterials"
        try:
            values = getter()
            if isinstance(values, (list, tuple)):
                return list(values), access
            access["getMaterials_return_type"] = type(values).__name__
        except Exception as error:
            access["errors"].append({"method": "getMaterials", "error": repr(error)})
    return [], access


def observe_materials(model, call_fields=True) -> dict:
    handles, access = material_handles(model)
    observed = {}
    for index, handle in enumerate(handles):
        key = f"material_{index}:{type(handle).__name__}"
        observed[key] = inspect_public_methods(handle, field_only=True) if call_fields else {}
    return {"access": access, "fields": observed, "status": "observed" if handles else "UNKNOWN_no_public_material_handle"}


def build_model(aka, manifest):
    aka.parseInput(str(INPUT / "material.dat"))
    mesh = aka.Mesh(3)
    mesh.read(str(INPUT / "block_2x2x4mm.msh"))
    model = aka.SolidMechanicsModelCohesive(mesh)
    model.initFull(_analysis_method=aka._explicit_lumped_mass, _is_extrinsic=True)
    model.assembleMassLumped()
    inserter = model.getElementInserter()
    inserter.addPhysicalSurface("interface")
    model.updateAutomaticInsertion()
    check_map = inserter.getCheckFacets()
    eligible = int(np.count_nonzero(check_map(aka._triangle_3)))
    eligible_total = int(check_map(aka._triangle_3).size)
    already_inserted = int(np.count_nonzero(inserter.getInsertionFacets()(aka._triangle_3)))
    if eligible != manifest["expected_insertion_facets"] or already_inserted != 0:
        raise RuntimeError(f"pre-step insertion selection mismatch {eligible}/{eligible_total}; already={already_inserted}")
    return mesh, model, eligible, eligible_total


def preflight() -> int:
    integrity = verify_and_record_snapshots()
    make_source_diff()
    import akantu as aka
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    mesh, model, eligible, eligible_total = build_model(aka, manifest)
    material_handles_found, material_access = material_handles(model)
    material_details = []
    for index, material in enumerate(material_handles_found):
        material_details.append({"index": index, "type": type(material).__name__,
                                 "fields": inspect_public_methods(material, field_only=True)})
    module_material_classes = {}
    for name in sorted(n for n in dir(aka) if "material" in n.lower() and "cohes" in n.lower()):
        obj = getattr(aka, name)
        module_material_classes[name] = {"type": type(obj).__name__, "doc": inspect.getdoc(obj),
                                         "public_methods": inspect_public_methods(obj) if isinstance(obj, type) else None}
    model_public = {}
    for name in sorted(n for n in dir(model) if not n.startswith("_") and
                       re.search(r"material|energy|traction|opening|damage|history|stress|dump|time|displacement|velocity|force", n, re.I)):
        member = getattr(model, name, None)
        model_public[name] = {"type": type(member).__name__, "doc": inspect.getdoc(member) if callable(member) else None}
    result = {
        "preflight_status": "PASS_no_solver_steps", "solver_steps": 0,
        "run_id": RUN_ID, "snapshot_integrity": integrity,
        "python": sys.executable, "python_version": platform.python_version(),
        "akantu_version": getattr(aka, "__version__", "UNKNOWN"),
        "numpy_version": np.__version__, "eligible_interface_facets": eligible,
        "eligible_face_slots": eligible_total, "mesh_nodes": int(mesh.getNbNodes()),
        "initial_native_cohesive_count": int(mesh.getNbElement(type=aka._cohesive_3d_6)),
        "model_public_api": model_public,
        "energy_key_public_probe": energy_values(model),
        "material_access": material_access, "material_objects": material_details,
        "module_cohesive_material_classes": module_material_classes,
        "dump_field_fallback": "UNKNOWN unless dump documentation exposes a cohesive traction/opening/history field; no private API used",
    }
    write_json(RUN_DIR / "preflight.json", result)
    print(json.dumps({k: result[k] for k in ("preflight_status", "solver_steps", "python", "akantu_version",
                                             "eligible_interface_facets", "eligible_face_slots", "energy_key_public_probe",
                                             "material_access", "module_cohesive_material_classes")}, indent=2))
    return 0


def array_copy(value) -> np.ndarray:
    return np.asarray(value).copy()


def mass_components(model) -> list[float]:
    mass = array_copy(model.getMass())
    if mass.ndim == 2:
        return [float(mass[:, axis].sum()) for axis in range(mass.shape[1])]
    return [float(mass.sum())]


def current_time_value(model):
    getter = getattr(model, "getTime", None)
    if not callable(getter):
        return None
    try:
        return float(getter())
    except Exception:
        return None


def sample_arrays(mesh, model, aka) -> dict:
    # Every call re-fetches and copies after any possible native node duplication.
    ntet = int(mesh.getNbElement(type=aka._tetrahedron_4))
    ncoh = int(mesh.getNbElement(type=aka._cohesive_3d_6))
    return {
        "rest_nodes": array_copy(mesh.getNodes()),
        "displacement": array_copy(model.getDisplacement()),
        "velocity": array_copy(model.getVelocity()),
        "internal_force": array_copy(model.getInternalForce()),
        "bulk_tetra_connectivity": array_copy(mesh.getConnectivity(aka._tetrahedron_4)),
        "cohesive_connectivity": (array_copy(mesh.getConnectivity(aka._cohesive_3d_6)) if ncoh else np.empty((0, 6), dtype=np.int64)),
        "mass": array_copy(model.getMass()),
        "mesh_nodes_count": int(mesh.getNbNodes()), "bulk_tetrahedra_count": ntet, "cohesive_elements_count": ncoh,
    }


def public_material_observation(handles, preflight_materials, preflight_access) -> dict:
    result = {"access": preflight_access, "fields": {},
              "status": "observed" if handles else "UNKNOWN_no_public_material_handle"}
    for index, handle in enumerate(handles):
        key = f"material_{index}:{type(handle).__name__}"
        plan = preflight_materials[index] if index < len(preflight_materials) else None
        if plan is None or plan.get("type") != type(handle).__name__:
            result["fields"][key] = {"status": "UNKNOWN_material_type_not_seen_during_preflight"}
            continue
        values = {}
        for name, entry in plan.get("fields", {}).items():
            if entry.get("noarg_call_status") == "ok":
                try:
                    values[name] = {"status": "ok", "value": jsonable(getattr(handle, name)())}
                except Exception as error:
                    values[name] = {"status": "error", "error": repr(error)}
        result["fields"][key] = values if values else {"status": "UNKNOWN_no_public_zero_argument_state_getter"}
    return result


def compare_probe_to_baseline(probe_rows: dict, baseline_path: Path) -> dict:
    with baseline_path.open(newline="", encoding="utf-8-sig") as stream:
        baseline_rows = {int(row["step"]): row for row in csv.DictReader(stream)}
    names = {
        "kinetic_J": "kinetic_J", "bulk_potential_J": "bulk_potential_J",
        "cohesive_dissipated_J": "cohesive_dissipated_J", "cohesive_reversible_J": "cohesive_reversible_J",
        "actuator_force_estimate_N": "actuator_force_estimate_N",
        "mass_component_0_kg": "mass_component_0_kg", "mass_component_1_kg": "mass_component_1_kg",
        "mass_component_2_kg": "mass_component_2_kg",
        "cohesive_insertions_this_step": "cohesive_insertions_this_step",
        "cohesive_insertions_cumulative": "cohesive_insertions_cumulative",
        "mesh_nodes": "mesh_nodes", "bulk_tetrahedra": "bulk_tetrahedra", "cohesive_elements": "cohesive_elements",
    }
    energy_names = set(ENERGY_FIELDS_FOR_MATCH)
    rows_compared = 0
    errors = []
    maxima = {key: {"max_abs": 0.0, "max_rel": 0.0} for key in names}
    for step in range(1, MAX_STEPS + 1):
        left, right = probe_rows.get(step), baseline_rows.get(step)
        if left is None or right is None:
            errors.append({"step": step, "reason": "missing row", "probe": left is not None, "baseline": right is not None})
            continue
        rows_compared += 1
        for probe_key, old_key in names.items():
            try:
                av, bv = float(left[probe_key]), float(right[old_key])
            except Exception as error:
                errors.append({"step": step, "field": probe_key, "reason": "unavailable", "error": repr(error)})
                continue
            diff = abs(av - bv)
            rel = diff / abs(bv) if bv else (0.0 if diff == 0 else math.inf)
            maxima[probe_key]["max_abs"] = max(maxima[probe_key]["max_abs"], diff)
            maxima[probe_key]["max_rel"] = max(maxima[probe_key]["max_rel"], rel)
            if probe_key in energy_names:
                passed = diff <= 1e-12 if bv == 0.0 else (diff <= 1e-12 and rel <= 1e-10)
            elif "insertions" in probe_key or probe_key in ("mesh_nodes", "bulk_tetrahedra", "cohesive_elements"):
                passed = av == bv
            else:
                passed = diff <= 1e-12 if bv == 0.0 else rel <= 1e-10
            if not passed:
                errors.append({"step": step, "field": probe_key, "probe": av, "baseline": bv,
                               "abs_diff": diff, "relative_diff": rel})
    return {
        "baseline_trace_sha256": sha256(baseline_path), "rows_compared_1_through_654": rows_compared,
        "comparison_tolerances": {"energy": "zero baseline: abs<=1e-12 J; otherwise abs<=max(1e-12 J, 1e-10*|baseline|)",
                                   "nonenergy_float": "zero baseline: abs<=1e-12; otherwise relative<=1e-10",
                                   "counts": "exact"},
        "maxima_by_field": maxima, "passed": rows_compared == MAX_STEPS and not errors,
        "error_count": len(errors), "errors_first_100": errors[:100],
        "position_comparison": "UNKNOWN_no_per_step_baseline_native_position_frames",
    }


def save_phase(phase_path: Path, dense_dir: Path, step: int, phase: str, model_time, bc_time,
               mesh, model, aka, current_work, previous_work, insertion_delta,
               material_cache, preflight_materials, preflight_material_access, arrays=True):
    state = sample_arrays(mesh, model, aka) if arrays else None
    energies = energy_values(model)
    material = public_material_observation(material_cache, preflight_materials, preflight_material_access)
    record = {
        "step": step, "phase": phase,
        "time_before_step_s": model_time,
        "prescribed_boundary_time_s": bc_time,
        "prescribed_boundary_time_is_set": phase in ("B", "C", "D"),
        "energy_getters": energies,
        "material_fields": material,
        "material_field_last_updated_phase": "UNKNOWN_native update/commit timing not exposed by inspected public API",
        "external_work_estimate_J": current_work,
        "external_work_before_D_update_J": previous_work if phase == "D" else None,
        "external_work_update_occurs_in_phase": "D_after_solve_from_sampled_reaction_and_delta_u",
        "cohesive_insertions_this_step": insertion_delta if phase in ("C", "D") else 0,
    }
    if state is not None:
        record.update({k: v for k, v in state.items() if not isinstance(v, np.ndarray)})
        out = dense_dir / f"step_{step:05d}_{phase}.npz"
        np.savez_compressed(out, **{k: v for k, v in state.items() if isinstance(v, np.ndarray)})
        record["native_array_snapshot"] = out.name
    with phase_path.open("a", encoding="utf-8") as stream:
        stream.write(json.dumps(jsonable(record), sort_keys=True) + "\n")
        stream.flush()
        os.fsync(stream.fileno())
    return state, energies


def run_once() -> int:
    if (RUN_DIR / "run_status.json").exists() or (RUN_DIR / "trace.csv").exists():
        raise RuntimeError("single-use run output already exists; refusing to overwrite or retry")
    integrity = verify_and_record_snapshots()
    make_source_diff()
    preflight_path = RUN_DIR / "preflight.json"
    if not preflight_path.is_file():
        raise RuntimeError("required preflight.json is absent; solver not started")
    preflight_result = json.loads(preflight_path.read_text(encoding="utf-8"))
    if preflight_result.get("preflight_status") != "PASS_no_solver_steps":
        raise RuntimeError("preflight did not pass; solver not started")
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    dense_dir = RUN_DIR / "dense_states"
    dense_dir.mkdir(exist_ok=True)
    start_utc = datetime.now(timezone.utc).isoformat()
    start = time.monotonic()
    status = {
        "run_id": RUN_ID, "case": "modeI_dt", "status": "RUNNING_SINGLE_SOLVER_ATTEMPT",
        "start_utc": start_utc, "argv": sys.argv, "cwd": os.getcwd(),
        "python": sys.executable, "python_version": platform.python_version(),
        "platform": platform.platform(), "solver_attempts": 1,
        "steps_original_case": 1707, "steps_cap": MAX_STEPS, "thread_limit": "OMP_NUM_THREADS=1",
        "dt_s": float(json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))["cases"]["modeI_dt"]["dt_s"]),
        "pre_run_snapshot_integrity": integrity,
        "instrumented_source_sha256": sha256(Path(__file__)),
        "baseline_trace_sha256": sha256(REFERENCE / "trace.csv"),
    }
    write_json(RUN_DIR / "run_status.json", status)
    try:
        import akantu as aka
        manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
        case = manifest["cases"]["modeI_dt"]
        mesh, model, eligible, eligible_total = build_model(aka, manifest)
        preflight_materials = preflight_result.get("material_objects", [])
        preflight_material_access = preflight_result.get("material_access", {})
        material_cache, _run_material_access = material_handles(model)
        original_nodes = array_copy(mesh.getNodes())
        original_tets = array_copy(mesh.getConnectivity(aka._tetrahedron_4))
        np.savez_compressed(RUN_DIR / "input_native_initial.npz", rest_nodes=original_nodes,
                            bulk_tetra_connectivity=original_tets,
                            cohesive_connectivity=np.empty((0, 6), dtype=np.int64),
                            mass=array_copy(model.getMass()))
        top_node_ids = np.flatnonzero(np.abs(original_nodes[:, 2] - 4.0e-3) < 1.0e-15)
        if len(top_node_ids) != 9:
            raise RuntimeError(f"expected 9 top nodes, found {len(top_node_ids)}")
        model.applyBC(aka.FixedValue(0.0, aka._z), "bottom")

        class RampTop(aka.DirichletFunctor):
            def __init__(self, amplitude, duration):
                super().__init__(aka._z)
                self.axis, self.amplitude, self.duration, self.time = int(aka._z), amplitude, duration, 0.0
            def set_time(self, value):
                self.time = value
            def __call__(self, node, flags, disp, coord):
                flags[self.axis] = True
                disp[self.axis] = self.amplitude * min(max(self.time / self.duration, 0.0), 1.0)

        ramp = RampTop(float(manifest["loading"]["top_displacement_m"]), float(manifest["loading"]["ramp_duration_s"]))
        dt = float(case["dt_s"])
        model.setTimeStep(dt)
        initial_mass = array_copy(model.getMass())
        mass_before_components = [float(initial_mass[:, axis].sum()) for axis in range(initial_mass.shape[1])]
        phase_path = RUN_DIR / "phase_evidence.jsonl"
        trace_path = RUN_DIR / "trace.csv"
        events_path = RUN_DIR / "events.json"
        trace_fields = [
            "step", "time_s", "time_before_step_s", "imposed_top_displacement_m",
            "top_internal_reaction_raw_N", "actuator_force_estimate_N", "external_work_J",
            "kinetic_J", "bulk_potential_J", "cohesive_dissipated_J", "cohesive_reversible_J",
            "cohesive_contact_native_key_J", "contact_alias_diagnostic_J", "energy_ledger_residual_J",
            "mass_component_0_kg", "mass_component_1_kg", "mass_component_2_kg",
            "mesh_nodes", "bulk_tetrahedra", "cohesive_elements", "cohesive_insertions_this_step",
            "cohesive_insertions_cumulative",
        ]
        current_work = 0.0
        previous_force = 0.0
        previous_displacement = 0.0
        total_native_insertions = 0
        inserted_seen = int(mesh.getNbElement(type=aka._cohesive_3d_6))
        events = []
        probe_rows = {}

        def write_d_row(writer, step, time_value, time_before, insertion_delta, internal_force=None):
            nonlocal current_work, previous_force, previous_displacement, total_native_insertions
            imposed = float(manifest["loading"]["top_displacement_m"]) * min(max(time_value / float(manifest["loading"]["ramp_duration_s"]), 0.0), 1.0)
            force = array_copy(model.getInternalForce()) if internal_force is None else internal_force
            raw = float(force[top_node_ids, 2].sum())
            actuator = -raw
            before = current_work
            current_work += 0.5 * (previous_force + actuator) * (imposed - previous_displacement)
            previous_force, previous_displacement = actuator, imposed
            e = energy_values(model)
            mass = mass_components(model)
            values = {name: (entry.get("value_J") if entry.get("status") == "ok" else None) for name, entry in e.items()}
            total_native_insertions += insertion_delta
            ledger_keys = ("kinetic_J", "bulk_potential_J", "cohesive_reversible_J", "cohesive_dissipated_J", "contact_alias_diagnostic_J")
            ledger_known = all(e[k].get("status") == "ok" for k in ENERGY_REQUESTS if k in ledger_keys)
            storage = sum(values[k] for k in ledger_keys) if ledger_known else None
            row = {
                "step": step, "time_s": time_value, "time_before_step_s": time_before,
                "imposed_top_displacement_m": imposed, "top_internal_reaction_raw_N": raw,
                "actuator_force_estimate_N": actuator, "external_work_J": current_work,
                "kinetic_J": values.get("kinetic_J"), "bulk_potential_J": values.get("bulk_potential_J"),
                "cohesive_dissipated_J": values.get("cohesive_dissipated_J"),
                "cohesive_reversible_J": values.get("cohesive_reversible_J"),
                "cohesive_contact_native_key_J": values.get("cohesive_contact_native_key_J"),
                "contact_alias_diagnostic_J": values.get("contact_alias_diagnostic_J"),
                "energy_ledger_residual_J": current_work - storage if storage is not None else None,
                "mass_component_0_kg": mass[0], "mass_component_1_kg": mass[1],
                "mass_component_2_kg": mass[2] if len(mass) > 2 else "",
                "mesh_nodes": int(mesh.getNbNodes()),
                "bulk_tetrahedra": int(mesh.getNbElement(type=aka._tetrahedron_4)),
                "cohesive_elements": int(mesh.getNbElement(type=aka._cohesive_3d_6)),
                "cohesive_insertions_this_step": insertion_delta,
                "cohesive_insertions_cumulative": total_native_insertions,
            }
            writer.writerow(row)
            trace_stream.flush()
            probe_rows[step] = row
            return row, e, before

        with trace_path.open("w", newline="", encoding="utf-8") as trace_stream:
            writer = csv.DictWriter(trace_stream, fieldnames=trace_fields)
            writer.writeheader()
            writer.writerow({
                "step": 0, "time_s": 0.0, "time_before_step_s": current_time_value(model),
                "imposed_top_displacement_m": 0.0, "top_internal_reaction_raw_N": 0.0,
                "actuator_force_estimate_N": 0.0, "external_work_J": 0.0,
                **{k: (v.get("value_J") if v.get("status") == "ok" else None) for k, v in energy_values(model).items()},
                "energy_ledger_residual_J": (0.0 if all(v.get("status") == "ok" for v in energy_values(model).values()) else None),
                "mass_component_0_kg": mass_before_components[0], "mass_component_1_kg": mass_before_components[1],
                "mass_component_2_kg": mass_before_components[2] if len(mass_before_components) > 2 else "",
                "mesh_nodes": int(mesh.getNbNodes()), "bulk_tetrahedra": int(mesh.getNbElement(type=aka._tetrahedron_4)),
                "cohesive_elements": 0, "cohesive_insertions_this_step": 0, "cohesive_insertions_cumulative": 0,
            })
            trace_stream.flush()

            for step in range(1, MAX_STEPS + 1):
                time_before = current_time_value(model)
                prescribed_time = step * dt
                insertion_delta = 0
                if step in FOCUS_STEPS:
                    save_phase(phase_path, dense_dir, step, "A", time_before, prescribed_time,
                               mesh, model, aka, current_work, current_work, 0,
                               material_cache, preflight_materials, preflight_material_access)
                ramp.set_time(prescribed_time)
                model.applyBC(ramp, "top")
                if step in FOCUS_STEPS:
                    save_phase(phase_path, dense_dir, step, "B", time_before, prescribed_time,
                               mesh, model, aka, current_work, current_work, 0,
                               material_cache, preflight_materials, preflight_material_access)
                api_return = int(model.checkCohesiveStress())
                cohesive_now = int(mesh.getNbElement(type=aka._cohesive_3d_6))
                insertion_delta = cohesive_now - inserted_seen
                if insertion_delta < 0 or api_return != insertion_delta:
                    raise RuntimeError(f"cohesive insertion count disagreement at step {step}: api={api_return}, native_delta={insertion_delta}")
                if insertion_delta:
                    events.append({"step": step, "prescribed_boundary_time_s": prescribed_time,
                                   "api_return": api_return, "cohesive_count": cohesive_now,
                                   "mesh_nodes": int(mesh.getNbNodes())})
                if step in FOCUS_STEPS:
                    save_phase(phase_path, dense_dir, step, "C", time_before, prescribed_time,
                               mesh, model, aka, current_work, current_work, insertion_delta,
                               material_cache, preflight_materials, preflight_material_access)
                model.solveStep()
                # Re-fetch post-solve state; these arrays are copied after any duplication at C.
                d_state = sample_arrays(mesh, model, aka) if step in FOCUS_STEPS else None
                force_for_row = d_state["internal_force"] if d_state is not None else None
                row, d_energies, work_before = write_d_row(writer, step, prescribed_time, time_before, insertion_delta, force_for_row)
                if step in FOCUS_STEPS:
                    # Save D from the already re-fetched state and record energy values from the D trace sample.
                    dense_path = dense_dir / f"step_{step:05d}_D.npz"
                    np.savez_compressed(dense_path, **{k: v for k, v in d_state.items() if isinstance(v, np.ndarray)})
                    phase_record = {
                        "step": step, "phase": "D", "time_before_step_s": time_before,
                        "solver_time_after_s": current_time_value(model),
                        "prescribed_boundary_time_s": prescribed_time, "prescribed_boundary_time_is_set": True,
                        "energy_getters": d_energies,
                        "material_fields": public_material_observation(material_cache, preflight_materials, preflight_material_access),
                        "material_field_last_updated_phase": "UNKNOWN_native update/commit timing not exposed by inspected public API",
                        "external_work_estimate_J": current_work, "external_work_before_D_update_J": work_before,
                        "external_work_update_occurs_in_phase": "D_after_solve_from_sampled_reaction_and_delta_u",
                        "cohesive_insertions_this_step": insertion_delta,
                        "native_array_snapshot": dense_path.name,
                        "mesh_nodes_count": d_state["mesh_nodes_count"], "bulk_tetrahedra_count": d_state["bulk_tetrahedra_count"],
                        "cohesive_elements_count": d_state["cohesive_elements_count"],
                    }
                    with phase_path.open("a", encoding="utf-8") as stream:
                        stream.write(json.dumps(jsonable(phase_record), sort_keys=True) + "\n")
                        stream.flush(); os.fsync(stream.fileno())
                inserted_seen = cohesive_now
                if time.monotonic() - start > MAX_WALL_S:
                    raise TimeoutError(f"300 s cap reached at step {step}; no retry")

        write_json(events_path, {"events": events, "early_stop_after_step": MAX_STEPS,
                                 "not_a_complete_1707_step_baseline": True})
        comparison = compare_probe_to_baseline(probe_rows, REFERENCE / "trace.csv")
        write_json(RUN_DIR / "baseline_comparison.json", comparison)
        status.update({
            "status": "DIAGNOSTIC_EARLY_STOP_NOT_COMPLETE",
            "end_utc": datetime.now(timezone.utc).isoformat(), "elapsed_s": time.monotonic() - start,
            "steps_completed": MAX_STEPS, "original_case_steps": 1707, "complete_case": False,
            "eligible_interface_facets": eligible, "eligible_facet_slots_total": eligible_total,
            "insertions": events, "comparison_passed": comparison["passed"],
            "native_state_capture_steps": sorted(FOCUS_STEPS), "phase_evidence": "phase_evidence.jsonl",
            "native_material_field_update_phase": "UNKNOWN unless values are exposed by public no-argument material getters",
            "external_work_definition": "same reaction/displacement trapezoid estimate as baseline; updated at D only, after solveStep",
            "solver_retries": 0,
        })
        write_json(RUN_DIR / "run_status.json", status)
        write_conclusion(status, comparison)
        return 0
    except Exception as error:
        status.update({"status": "FAILED_NO_RETRY", "end_utc": datetime.now(timezone.utc).isoformat(),
                       "elapsed_s": time.monotonic() - start, "exception": repr(error),
                       "traceback": traceback.format_exc(), "solver_retries": 0})
        write_json(RUN_DIR / "run_status.json", status)
        write_conclusion(status, None)
        raise


def write_conclusion(status, comparison) -> None:
    result = {
        "run_id": RUN_ID, "status": status.get("status"),
        "solver_attempts": status.get("solver_attempts", 1), "solver_retries": 0,
        "steps_completed": status.get("steps_completed", 0), "complete_1707_step_case": False,
        "first_negative_dissipated_step": "UNKNOWN_until_trace_parse",
        "baseline_comparison": comparison["passed"] if comparison else "NOT_COMPLETED",
        "baseline_position_comparison": "UNKNOWN_no_per_step_native_position_frames",
        "native_traction_opening_damage": "UNKNOWN_unless_public_material_getters_returned_values; inspect phase_evidence.jsonl",
        "native_material_update_phase": "UNKNOWN_unless observed field values changed across phases",
        "contact_key_audit": "see phase_evidence.jsonl: cohesive contact is canonical; contact is legacy diagnostic only",
        "source_diff": "../experiments/cohesive_event_probe_20261004/run_contract_instrumented.diff",
        "phase_evidence": "phase_evidence.jsonl", "full_trace": "trace.csv",
        "source_sha256": status.get("pre_run_snapshot_integrity", {}).get("source_sha256_actual"),
        "manifest_sha256": status.get("pre_run_snapshot_integrity", {}).get("manifest_sha256_actual"),
        "wheel_source_identity": "UNKNOWN_installed akantu 5.0.7.post1 vs reviewed v5.0.7 tag",
    }
    trace_path = RUN_DIR / "trace.csv"
    if trace_path.is_file():
        try:
            with trace_path.open(newline="", encoding="utf-8-sig") as stream:
                rows = list(csv.DictReader(stream))
            negatives = [row for row in rows if row.get("cohesive_dissipated_J") not in (None, "") and float(row["cohesive_dissipated_J"]) < 0]
            if negatives:
                result["first_negative_dissipated_step"] = int(negatives[0]["step"])
                result["first_negative_dissipated_J"] = float(negatives[0]["cohesive_dissipated_J"])
                event_steps = [e["step"] for e in json.loads((RUN_DIR / "events.json").read_text(encoding="utf-8"))["events"]]
                result["steps_after_first_insertion"] = int(negatives[0]["step"]) - min(event_steps) if event_steps else None
        except Exception as error:
            result["trace_parse_error"] = repr(error)
    write_json(RUN_DIR / "probe_conclusion.json", result)


def main() -> int:
    if len(sys.argv) == 2 and sys.argv[1] == "--preflight":
        return preflight()
    if len(sys.argv) == 2 and sys.argv[1] == "--run":
        return run_once()
    raise SystemExit("usage: probe_driver.py --preflight | --run")


if __name__ == "__main__":
    sys.exit(main())
