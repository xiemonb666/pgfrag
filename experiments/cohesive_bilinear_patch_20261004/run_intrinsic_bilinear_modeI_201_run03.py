from __future__ import annotations

import hashlib
import json
import os
import sys
import time
import traceback
from pathlib import Path

import akantu as aka
import numpy as np

from parse_cohesive_vtu import new_pvtu_files, parse_pvtu

EXP = Path(__file__).resolve().parent
INPUT = EXP / "input"
MANIFEST_PATH = EXP / "manifest_bilinear_modeI_201_run03.json"
OUT = Path("/home/xiemo/cohesive_bilinear_patch_20261004/controlled_intrinsic_bilinear_modeI_201_run03")
DT = 1.0e-9
SIGMA = 70.0e6
GC = 7.6
DELTA_0 = 2.1714285714285714e-8
DELTA_C = 2.0 * GC / SIGMA
ENERGY_KEYS = ("kinetic", "potential", "dissipated", "reversible", "cohesive contact")
FIELD_NAMES = ("opening", "tractions", "damage", "delta max")
OBS = {}
HISTORY = []
SOLVER_STEPS_COMPLETED = 0


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def arr_copy(value) -> np.ndarray:
    return np.asarray(value).copy()


def finite_number(value):
    x = float(value)
    return x if np.isfinite(x) else None


def sha_check(path: Path, expected: str):
    actual = sha256(path)
    if actual != expected:
        raise RuntimeError(f"frozen input hash mismatch: {path.name} {actual} != {expected}")
    return actual


def write_json(path: Path, value):
    def safe(item):
        if isinstance(item, dict):
            return {str(k): safe(v) for k, v in item.items()}
        if isinstance(item, (list, tuple)):
            return [safe(v) for v in item]
        if isinstance(item, np.ndarray):
            return safe(item.tolist())
        if isinstance(item, (np.floating, float)):
            return finite_number(item)
        if isinstance(item, np.integer):
            return int(item)
        return item
    path.write_text(json.dumps(safe(value), indent=2, allow_nan=False) + "\n", encoding="utf-8")


def energy_values(model):
    return {key: float(model.getEnergy(key)) for key in ENERGY_KEYS}


def time_value(model):
    try:
        return float(model.getTime())
    except Exception:
        return None


def component_sums(a: np.ndarray):
    return [float(a[:, i].sum()) for i in range(a.shape[1])] if a.ndim == 2 else [float(a.sum())]


def grid_points():
    a = np.linspace(0.0, DELTA_0, 21)
    b = np.linspace(DELTA_0, DELTA_C, 161)[1:]
    c = np.linspace(DELTA_C, 1.2 * DELTA_C, 21)[1:]
    points = np.concatenate((a, b, c))
    if len(points) != 201 or np.any(np.diff(points) <= 0):
        raise RuntimeError(f"loading grid invalid: {len(points)} points")
    return points


def save_aggregate():
    if not HISTORY:
        return
    keys = ("opening", "tractions", "damage", "delta max", "displacement",
            "internal_force", "mass", "rest_nodes", "bulk_tets", "cohesive_connectivity")
    data = {k.replace(" ", "_"): np.stack([r[k] for r in HISTORY]) for k in keys}
    data.update({
        "point_index": np.asarray([r["point_index"] for r in HISTORY], dtype=np.int64),
        "solver_step": np.asarray([r["solver_step"] for r in HISTORY], dtype=np.int64),
        "commanded_opening_m": np.asarray([r["commanded_opening_m"] for r in HISTORY], dtype=np.float64),
        "energies_J": np.asarray([[r["energies_J"][key] for key in ENERGY_KEYS] for r in HISTORY], dtype=np.float64),
        "mass_component_sums_kg": np.asarray([r["mass_component_sums_kg"] for r in HISTORY], dtype=np.float64),
        "counts": np.asarray([[r["counts"][k] for k in ("nodes", "bulk_tets", "cohesive_3d6")] for r in HISTORY], dtype=np.int64),
        "upper_internal_force_z_sum_raw_N": np.asarray([r["upper_internal_force_z_sum_raw_N"] for r in HISTORY], dtype=np.float64),
        "post_solve_boundary_readback_error_m": np.asarray(
            [r["post_solve_boundary_readback_error_m"] for r in HISTORY], dtype=np.float64),
    })
    np.savez_compressed(OUT / "native_history_partial_or_complete.npz", **data)


def init_and_run():
    global SOLVER_STEPS_COMPLETED
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    sha_check(INPUT / "block_2x2x4mm.msh", manifest["mesh"]["sha256"])
    sha_check(INPUT / "material.dat", manifest["material"]["sha256"])
    sha_check(Path(__file__), manifest["sources"]["driver"]["sha256"])
    sha_check(Path(__file__).with_name("parse_cohesive_vtu.py"), manifest["sources"]["vtk_parser"]["sha256"])
    aka.parseInput(str(INPUT / "material.dat"))
    mesh = aka.Mesh(3)
    mesh.read(str(INPUT / "block_2x2x4mm.msh"))
    groups = [{"name": g.getName(), "elements": len(g)} for g in mesh.elementGroups()]
    n0, t0 = arr_copy(mesh.getNodes()), arr_copy(mesh.getConnectivity(aka._tetrahedron_4))
    if n0.shape != (45, 3) or t0.shape != (96, 4):
        raise RuntimeError(f"input mesh shape unexpected: nodes={n0.shape} tets={t0.shape}")
    model = aka.SolidMechanicsModelCohesive(mesh)
    model.initFull(_analysis_method=aka._explicit_lumped_mass, _is_extrinsic=False)
    model.assembleMassLumped()
    model.setTimeStep(DT)
    rest = arr_copy(mesh.getNodes())
    tets = arr_copy(mesh.getConnectivity(aka._tetrahedron_4))
    coh = arr_copy(mesh.getConnectivity(aka._cohesive_3d_6))
    mass = arr_copy(model.getMass())
    if (len(rest), len(tets), len(coh)) != (54, 96, 8):
        raise RuntimeError(f"native input counts unexpected: {(len(rest), len(tets), len(coh))}")
    centroids = rest[tets].mean(axis=1)[:, 2]
    lower = set(map(int, tets[centroids < 0.002].ravel()))
    upper = set(map(int, tets[centroids > 0.002].ravel()))
    if lower & upper or lower | upper != set(range(54)) or (len(lower), len(upper)) != (27, 27):
        raise RuntimeError("native centroid half partition must be disjoint lower27/upper27 covering all nodes")
    lower_ids, upper_ids = np.asarray(sorted(lower)), np.asarray(sorted(upper))
    displacement_view = np.asarray(model.getDisplacement())
    blocked_view = np.asarray(model.getBlockedDOFs())
    if displacement_view.shape != (54, 3) or blocked_view.shape != (54, 3):
        raise RuntimeError(f"native writable view shape unexpected: {displacement_view.shape}, {blocked_view.shape}")
    if not displacement_view.flags.writeable or not blocked_view.flags.writeable:
        raise RuntimeError("native displacement/blocked_dofs view is not writable")
    blocked_view[:, :] = True
    blocked_initial_readback = np.asarray(model.getBlockedDOFs())
    if blocked_initial_readback.shape != (54, 3) or not np.all(blocked_initial_readback):
        raise RuntimeError("native getBlockedDOFs readback failed immediately after initial prescription")
    displacement_view[lower_ids, :] = 0.0
    displacement_view[upper_ids, :] = 0.0
    for name in FIELD_NAMES:
        model.addDumpFieldToDumper("cohesive elements", name)
    energies0 = energy_values(model)
    counts0 = {"nodes": len(rest), "bulk_tets": len(tets), "cohesive_3d6": len(coh)}
    mass_sums0 = component_sums(mass)
    OBS.update({
        "physical_groups": groups,
        "counts_initial": counts0,
        "lower_node_ids": lower_ids.tolist(),
        "upper_node_ids": upper_ids.tolist(),
        "initial_mass_component_sums_kg": mass_sums0,
        "initial_energies_J": energies0,
        "solver_dt_s": DT,
        "grid_point_count": 201,
        "solver_steps_completed": 0,
    })
    np.savez_compressed(OUT / "initial_state.npz", input_nodes=n0, input_bulk_tets=t0,
                        rest_nodes=rest, bulk_tets=tets, cohesive_connectivity=coh,
                        mass=mass, lower_node_ids=lower_ids, upper_node_ids=upper_ids,
                        blocked_dofs=blocked_initial_readback.copy())
    points = grid_points()
    model.setTimeStep(DT)
    trace = OUT / "trace.jsonl"
    started = time.monotonic()
    if trace.exists():
        raise RuntimeError("trace already exists; refusing rerun")
    with trace.open("x", encoding="utf-8") as stream:
        for point_index, s in enumerate(points):
            disp_live = np.asarray(model.getDisplacement())
            blocked_live = np.asarray(model.getBlockedDOFs())
            blocked_live[:, :] = True
            disp_live[lower_ids, :] = 0.0
            disp_live[upper_ids, :] = (0.0, 0.0, float(s))
            blocked_readback = np.asarray(model.getBlockedDOFs())
            if blocked_readback.shape != (54, 3) or not np.all(blocked_readback):
                raise RuntimeError(f"native blocked-DOF pre-solve readback failed at point {point_index}")
            prescribed = arr_copy(model.getDisplacement())
            pre_error = max(float(np.max(np.abs(prescribed[lower_ids]))),
                            float(np.max(np.abs(prescribed[upper_ids] - (0.0, 0.0, float(s))))))
            if pre_error > 1e-12:
                raise RuntimeError(f"prescribed displacement readback error at point {point_index}: {pre_error}")
            time_before = time_value(model)
            model.solveStep()
            SOLVER_STEPS_COMPLETED += 1
            after_displacement = arr_copy(model.getDisplacement())
            post_error = max(float(np.max(np.abs(after_displacement[lower_ids]))),
                             float(np.max(np.abs(after_displacement[upper_ids] - (0.0, 0.0, float(s))))))
            internal_force = arr_copy(model.getInternalForce())
            mass_now = arr_copy(model.getMass())
            nodes_now = arr_copy(mesh.getNodes())
            tets_now = arr_copy(mesh.getConnectivity(aka._tetrahedron_4))
            coh_now = arr_copy(mesh.getConnectivity(aka._cohesive_3d_6))
            counts = {"nodes": len(nodes_now), "bulk_tets": len(tets_now), "cohesive_3d6": len(coh_now)}
            if counts != counts0:
                raise RuntimeError(f"native topology changed at point {point_index}: {counts}")
            before = {p.resolve() for p in OUT.rglob("*.pvtu")}
            model.dump("cohesive elements", point_index)
            created = new_pvtu_files(OUT, before)
            expected_name = f"_{point_index:04d}.pvtu"
            created = {p for p in created if p.name.endswith(expected_name)}
            if len(created) != 1:
                raise RuntimeError(f"dump index {point_index} did not create exactly one distinct PVTU: {sorted(map(str, created))}")
            pvtu = next(iter(created))
            parsed = parse_pvtu(pvtu, OUT)
            if parsed["step_index"] != point_index:
                raise RuntimeError(f"parsed dump step {parsed['step_index']} != requested {point_index}")
            fields = parsed["fields"]
            energies = energy_values(model)
            z_force = float(internal_force[upper_ids, 2].sum())
            record = {
                "point_index": point_index,
                "solver_step": SOLVER_STEPS_COMPLETED,
                "dump_index": point_index,
                "commanded_opening_m": float(s),
                "prescribed_readback_error_m": pre_error,
                "native_blocked_dofs_all_true_pre_solve": bool(np.all(blocked_readback)),
                "post_solve_boundary_readback_error_m": post_error,
                "time_before_step_s": time_before,
                "time_after_step_s": time_value(model),
                "counts": counts,
                "mass_component_sums_kg": component_sums(mass_now),
                "upper_internal_force_z_sum_raw_N": finite_number(z_force),
                "energies_J": {key: finite_number(value) for key, value in energies.items()},
                "native_dump": {"pvtu": str(pvtu.relative_to(OUT)), "vtu": str(Path(parsed["vtu"]).relative_to(OUT)),
                                "step_index": parsed["step_index"], "cells": parsed["cells"], "points": parsed["points"]},
                "fields": {name: {"shape": list(value.shape),
                                  "finite": bool(np.isfinite(value).all()),
                                  "minimum": finite_number(np.min(value)),
                                  "maximum": finite_number(np.max(value))}
                           for name, value in fields.items()},
                "opening_norm_mean_m": finite_number(np.linalg.norm(fields["opening"], axis=1).mean()),
                "traction_norm_mean_Pa": finite_number(np.linalg.norm(fields["tractions"], axis=1).mean()),
                "damage_mean": finite_number(fields["damage"].mean()),
                "delta_max_mean_m": finite_number(fields["delta max"].mean()),
            }
            sample_path = OUT / "samples" / f"point_{point_index:03d}.npz"
            np.savez_compressed(sample_path, point_index=np.int64(point_index),
                solver_step=np.int64(SOLVER_STEPS_COMPLETED), commanded_opening_m=np.float64(s),
                opening=fields["opening"], tractions=fields["tractions"], damage=fields["damage"],
                delta_max=fields["delta max"], displacement=after_displacement,
                internal_force=internal_force, mass=mass_now, rest_nodes=nodes_now,
                bulk_tets=tets_now, cohesive_connectivity=coh_now,
                energies_J=np.asarray([energies[k] for k in ENERGY_KEYS], dtype=np.float64),
                mass_component_sums_kg=np.asarray(component_sums(mass_now), dtype=np.float64))
            stream.write(json.dumps(record, allow_nan=False) + "\n")
            stream.flush()
            os.fsync(stream.fileno())
            record.update({
                "opening": fields["opening"].copy(),
                "tractions": fields["tractions"].copy(),
                "damage": fields["damage"][:, 0].copy(),
                "delta max": fields["delta max"][:, 0].copy(),
                "displacement": after_displacement.copy(),
                "internal_force": internal_force.copy(),
                "mass": mass_now.copy(),
                "rest_nodes": nodes_now.copy(),
                "bulk_tets": tets_now.copy(),
                "cohesive_connectivity": coh_now.copy(),
            })
            HISTORY.append(record)
            if time.monotonic() - started > 300.0:
                raise TimeoutError(f"300 s wall limit reached after solver step {SOLVER_STEPS_COMPLETED}")
    all_pvtu = sorted(OUT.rglob("*.pvtu"))
    indices = [int(p.stem.rsplit("_", 1)[1]) for p in all_pvtu if p.stem.rsplit("_", 1)[-1].isdigit()]
    if len(all_pvtu) != 201 or sorted(indices) != list(range(201)) or len(HISTORY) != 201:
        raise RuntimeError(f"final 201-frame integrity failed: pvtu={len(all_pvtu)}, indices={indices}, records={len(HISTORY)}")
    save_aggregate()
    OBS["solver_steps_completed"] = SOLVER_STEPS_COMPLETED
    OBS["dump_frames_and_parsed_records"] = len(HISTORY)
    OBS["status"] = "RUN_201_POINTS_COMPLETE"
    return OBS


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    if (OUT / "run_started.json").exists() or (OUT / "run_status.json").exists():
        raise SystemExit("run output already started; refusing retry")
    (OUT / "samples").mkdir(exist_ok=True)
    start = {"argv": sys.argv, "source_sha256": sha256(Path(__file__)),
             "manifest_sha256": sha256(MANIFEST_PATH), "solver_attempts": 1}
    write_json(OUT / "run_started.json", start)
    result = {}
    try:
        os.chdir(OUT)
        result = init_and_run()
    except Exception as e:
        result = {**OBS, "status": "FAIL_NO_RETRY", "error": repr(e),
                  "traceback": traceback.format_exc(), "solver_steps_completed": SOLVER_STEPS_COMPLETED}
    finally:
        try:
            save_aggregate()
        except Exception as e:
            result = {**result, "aggregate_save_error": repr(e)}
        result = {**result, "argv": sys.argv, "solver_steps_completed": SOLVER_STEPS_COMPLETED}
        result.setdefault("status", "FAIL_NO_RETRY")
        write_json(OUT / "run_status.json", result)
    print(json.dumps(result, allow_nan=False, sort_keys=True))
    return 0 if result.get("status") == "RUN_201_POINTS_COMPLETE" else 2


if __name__ == "__main__":
    raise SystemExit(main())
