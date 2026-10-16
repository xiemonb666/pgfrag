from __future__ import annotations

import json
import os
import sys
import time
import traceback
from pathlib import Path
from typing import Any

import numpy as np

from free_patch_contract import (
    build_native_lineage,
    input_interface_contract,
    quadrature_signed_work,
    read_gmsh22,
    sha256,
    strict_json_dump,
    surface_nodes,
)


EXP = Path(__file__).resolve().parent
CASE_PATH = EXP / "case.json"
MESH_PATH = EXP / "block_2x2x4mm.msh"
MATERIAL_PATH = EXP / "material.dat"
OUT = Path("/home/xiemo/free_cohesive_patch_20261015/free_native_intrinsic_modeI_patch_20261015_v1")
ENERGY_KEYS = ("kinetic", "potential", "dissipated", "reversible", "cohesive contact")
FIELD_NAMES = ("opening", "tractions", "damage", "delta max")
SOLVER_STEPS_COMPLETED = 0
HISTORY_FRAMES_SAVED = 0
OBS: dict[str, Any] = {}


def safe_float(value: Any) -> float | None:
    x = float(value)
    return x if np.isfinite(x) else None


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, np.ndarray):
        return json_safe(value.tolist())
    if isinstance(value, (np.floating, float)):
        return safe_float(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    return value


def write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(json_safe(value), indent=2, allow_nan=False) + "\n", encoding="utf-8")


class RawHistory:
    """Preallocated NPY arrays preserve every native value without JSON coercion."""

    SPEC = {
        "opening": ((24, 3), np.float64),
        "tractions": ((24, 3), np.float64),
        "damage": ((24, 1), np.float64),
        "delta_max": ((24, 1), np.float64),
        "traction_dot_opening_Pa_m_raw": ((24,), np.float64),
        "interface_work_increment_J_by_qp_raw": ((24,), np.float64),
        "energies_J_raw": ((len(ENERGY_KEYS),), np.float64),
        "displacement_m": ((54, 3), np.float64),
        "internal_force_N_raw": ((54, 3), np.float64),
        "acceleration_m_s2_raw": ((54, 3), np.float64),
        "external_force_N_raw": ((54, 3), np.float64),
        "mass_kg_raw": ((54, 3), np.float64),
        "boundary_readback_error_m": ((), np.float64),
        "native_work_cumulative_J_raw": ((), np.float64),
        "native_energy_total_J_raw": ((), np.float64),
        "native_work_vs_total_residual_J_raw": ((), np.float64),
        "actuator_reaction_estimate_N_raw": ((), np.float64),
        "top_internal_force_z_sum_N_raw": ((), np.float64),
        "top_inertial_force_z_sum_N_raw": ((), np.float64),
        "top_external_force_z_sum_N_raw": ((), np.float64),
        "actuator_work_cumulative_J_raw": ((), np.float64),
        "energy_balance_residual_J_raw": ((), np.float64),
        "commanded_top_displacement_m": ((), np.float64),
        "topology_counts": ((3,), np.int64),
    }

    def __init__(self, directory: Path, frames: int):
        directory.mkdir(parents=True, exist_ok=False)
        self.arrays = {
            name: np.lib.format.open_memmap(
                directory / f"{name}.npy", mode="w+", dtype=dtype,
                shape=(frames, *shape),
            )
            for name, (shape, dtype) in self.SPEC.items()
        }
        self.frames = frames

    def write(self, index: int, row: dict[str, Any]) -> None:
        if index < 0 or index >= self.frames:
            raise IndexError(f"history frame out of bounds: {index}")
        for name, array in self.arrays.items():
            if name not in row:
                raise KeyError(f"history row missing {name}")
            value = np.asarray(row[name], dtype=array.dtype)
            if value.shape != array.shape[1:]:
                raise ValueError(f"{name} shape {value.shape} != {array.shape[1:]}")
            array[index] = value
        for array in self.arrays.values():
            array.flush()


def commanded_opening(case: dict[str, Any], solver_step: int) -> float:
    load = case["loading"]
    ramp_steps = int(load["ramp_steps"])
    target = float(load["top_target_displacement_m"])
    if solver_step >= ramp_steps:
        return target
    return 0.5 * target * (1.0 - float(np.cos(np.pi * solver_step / ramp_steps)))


def native_fields(model: Any, aka: Any) -> dict[str, np.ndarray]:
    fields = {
        name: np.asarray(model.flattenInternalReal(name, aka._ek_cohesive, aka._not_ghost)(aka._cohesive_3d_6, aka._not_ghost)).copy()
        for name in FIELD_NAMES
    }
    expected = {"opening": (24, 3), "tractions": (24, 3), "damage": (24, 1), "delta max": (24, 1)}
    for name, shape in expected.items():
        if fields[name].shape != shape:
            raise RuntimeError(f"native {name} shape {fields[name].shape} != {shape}")
    return fields


def native_energies(model: Any) -> np.ndarray:
    # getEnergy is a read of stored constitutive state in the reviewed v5.0.7 source.
    # Read each channel once per state; keep the exact returned signs in the NPY history.
    return np.asarray([float(model.getEnergy(key)) for key in ENERGY_KEYS], dtype=np.float64)


def force_terms(model: Any) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    internal = np.asarray(model.getInternalForce()).copy()
    acceleration = np.asarray(model.getAcceleration()).copy()
    external = np.asarray(model.getExternalForce()).copy()
    mass = np.asarray(model.getMass()).copy()
    displacement = np.asarray(model.getDisplacement()).copy()
    for name, value in (("internal force", internal), ("acceleration", acceleration),
                        ("external force", external), ("mass", mass), ("displacement", displacement)):
        if value.shape != (54, 3):
            raise RuntimeError(f"native {name} shape {value.shape} != (54, 3)")
    return displacement, internal, acceleration, external, mass


def boundary_error(displacement: np.ndarray, top: np.ndarray, bottom: np.ndarray,
                   pin_x: np.ndarray, pin_y: np.ndarray, commanded: float) -> float:
    errors = [
        float(np.max(np.abs(displacement[top, 2] - commanded), initial=0.0)),
        float(np.max(np.abs(displacement[bottom, 2]), initial=0.0)),
        float(np.max(np.abs(displacement[pin_x, 0]), initial=0.0)),
        float(np.max(np.abs(displacement[pin_x, 1]), initial=0.0)),
        float(np.max(np.abs(displacement[pin_y, 1]), initial=0.0)),
    ]
    return max(errors, default=0.0)


def energy_record(energies: np.ndarray) -> dict[str, float | None]:
    return {key: safe_float(value) for key, value in zip(ENERGY_KEYS, energies, strict=True)}


def run_native() -> dict[str, Any]:
    global SOLVER_STEPS_COMPLETED, HISTORY_FRAMES_SAVED
    import akantu as aka

    case = json.loads(CASE_PATH.read_text(encoding="utf-8"))
    if case["schema"] != "akantu-free-cohesive-patch-v1":
        raise RuntimeError("unexpected case schema")
    if sha256(MESH_PATH) != case["geometry"]["mesh_sha256"]:
        raise RuntimeError("frozen mesh SHA-256 mismatch")
    if sha256(MATERIAL_PATH) != case["material"]["material_sha256"]:
        raise RuntimeError("frozen material SHA-256 mismatch")
    input_mesh = read_gmsh22(MESH_PATH)
    source_face_contract = input_interface_contract(input_mesh)
    if len(input_mesh["nodes"]) != case["geometry"]["input_nodes"]:
        raise RuntimeError("input node count differs from frozen case")
    if len(input_mesh["tetrahedra"]) != case["geometry"]["input_bulk_tetrahedra"]:
        raise RuntimeError("input tetrahedron count differs from frozen case")
    if source_face_contract["input_face_count"] != case["geometry"]["central_interface"]["input_triangles"]:
        raise RuntimeError("input interface triangle count differs from frozen case")
    if not np.isclose(source_face_contract["rest_area_sum_m2"],
                      case["geometry"]["central_interface"]["rest_area_m2"], rtol=0.0, atol=1e-12):
        raise RuntimeError("input interface area differs from frozen case")

    aka.parseInput(str(MATERIAL_PATH))
    mesh = aka.Mesh(3)
    mesh.read(str(MESH_PATH))
    input_nodes_from_api = np.asarray(mesh.getNodes()).copy()
    input_tets_from_api = np.asarray(mesh.getConnectivity(aka._tetrahedron_4)).copy()
    if input_nodes_from_api.shape != (45, 3) or input_tets_from_api.shape != (96, 4):
        raise RuntimeError(f"Akantu input mesh shapes unexpected: {input_nodes_from_api.shape}, {input_tets_from_api.shape}")
    if not np.array_equal(input_nodes_from_api, input_mesh["nodes"]):
        raise RuntimeError("Akantu node rows do not match the frozen MSH node row order")
    if not np.array_equal(input_tets_from_api, np.asarray([e["nodes"] for e in input_mesh["tetrahedra"]])):
        raise RuntimeError("Akantu tetrahedron rows do not match the frozen MSH tetrahedron row order")

    model = aka.SolidMechanicsModelCohesive(mesh)
    model.initFull(_analysis_method=aka._explicit_lumped_mass, _is_extrinsic=False)
    model.assembleMassLumped()
    bulk_stable_dt = float(model.getStableTimeStep())
    dt = float(case["loading"]["time_step_s"])
    bulk_fraction = float(case["loading"]["bulk_stable_step_fraction_limit"])
    max_bulk_dt = bulk_fraction * bulk_stable_dt
    if not (np.isfinite(bulk_stable_dt) and bulk_stable_dt > 0.0 and dt <= max_bulk_dt):
        raise RuntimeError(
            f"bulk-only stable-step check failed: requested={dt}, bulk_stable={bulk_stable_dt}, "
            f"fraction={bulk_fraction}, max={max_bulk_dt}"
        )

    stability = case["stability"]
    total_status = str(stability["total_stability_status"])
    lambda_upper = stability.get("lambda_max_upper_bound_per_s2")
    total_stable_dt = None
    if lambda_upper is not None and np.isfinite(float(lambda_upper)) and float(lambda_upper) > 0.0:
        # Standard central-difference upper limit from the reviewed whole-system
        # mass-normalized stiffness bound. The case's fraction is applied below.
        total_stable_dt = 2.0 / float(np.sqrt(float(lambda_upper)))
    total_fraction = float(stability["fraction_limit"])
    total_gate = (
        total_status == "ESTABLISHED_PARENT_REVIEWED_BOUND"
        and total_stable_dt is not None
        and np.isfinite(total_stable_dt)
        and total_stable_dt > 0.0
        and 0.0 < total_fraction <= 1.0
        and dt <= total_fraction * total_stable_dt
    )
    OBS["stability_contract"] = {
        "bulk_get_stable_time_step_s": bulk_stable_dt,
        "bulk_fraction_limit": bulk_fraction,
        "bulk_check_max_dt_s": max_bulk_dt,
        "bulk_check_is_total_cohesive_cfl": False,
        "total_stability_status": total_status,
        "lambda_max_upper_bound_per_s2": lambda_upper,
        "derived_total_stable_time_step_s": total_stable_dt,
        "total_fraction_limit": total_fraction,
        "solver_dt_s": dt,
        "total_gate_passed_before_step_1": bool(total_gate),
    }
    if not total_gate:
        raise RuntimeError(
            "total cohesive-plus-bulk stability gate is not established or failed; "
            f"status={total_status}, lambda_upper={lambda_upper}, derived_dt={total_stable_dt}, "
            f"requested_dt={dt}, fraction={total_fraction}. Stop before solver step 1."
        )
    model.setTimeStep(dt)

    native_rest = np.asarray(mesh.getNodes()).copy()
    native_tets = np.asarray(mesh.getConnectivity(aka._tetrahedron_4)).copy()
    native_cohesive = np.asarray(mesh.getConnectivity(aka._cohesive_3d_6)).copy()
    mass0 = np.asarray(model.getMass()).copy()
    # Bind the supplied initial stiffness bound to the measured native topology
    # and post-insertion mass from which the parent calculated it.
    reference_path = EXP / stability["preflight_state_file"]
    if sha256(reference_path) != stability["preflight_state_sha256"]:
        raise RuntimeError("initial-bound preflight archive hash changed")
    with np.load(reference_path, allow_pickle=False) as reference:
        for key, actual in (("native_nodes", native_rest), ("native_bulk_tets", native_tets),
                            ("native_cohesive_connectivity", native_cohesive), ("mass_kg", mass0)):
            if not np.array_equal(reference[key], actual):
                raise RuntimeError(f"native state no longer matches the initial-bound input: {key}")
    counts0 = (len(native_rest), len(native_tets), len(native_cohesive))
    if counts0 != (54, 96, 8):
        raise RuntimeError(f"native topology after insertion unexpected: {counts0}")
    lineage = build_native_lineage(
        input_mesh, native_rest, native_tets, native_cohesive,
        expected_area_m2=float(case["geometry"]["central_interface"]["rest_area_m2"]),
    )
    native_to_input = {int(k): int(v) for k, v in lineage["native_node_to_input_node_zero_based"].items()}
    input_to_native: dict[int, list[int]] = {}
    for native_id, input_id in native_to_input.items():
        input_to_native.setdefault(input_id, []).append(native_id)

    bc = case["boundary_conditions"]
    top_input = surface_nodes(input_mesh, "top")
    bottom_input = surface_nodes(input_mesh, "bottom")
    top = np.asarray(sorted(n for node in top_input for n in input_to_native[node]), dtype=np.int64)
    bottom = np.asarray(sorted(n for node in bottom_input for n in input_to_native[node]), dtype=np.int64)
    pin_native: dict[int, list[int]] = {}
    for pin in bc["minimal_in_plane_pins"]:
        source_id = int(pin["input_node_id_zero_based"])
        native_ids = sorted(input_to_native[source_id])
        if len(native_ids) != 1 or not np.array_equal(native_rest[native_ids[0]], np.asarray(pin["coordinate_m"])):
            raise RuntimeError(f"pin node identity/rest coordinate check failed: {source_id}")
        pin_native[source_id] = native_ids
    pin_x = np.asarray(pin_native[0], dtype=np.int64)
    pin_y = np.asarray(pin_native[2], dtype=np.int64)
    if len(top) != 9 or len(bottom) != 9:
        raise RuntimeError(f"physical top/bottom input groups map to unexpected native rows: {len(top)}, {len(bottom)}")

    displacement = np.asarray(model.getDisplacement())
    blocked = np.asarray(model.getBlockedDOFs())
    if displacement.shape != (54, 3) or blocked.shape != (54, 3):
        raise RuntimeError(f"native displacement/blocked shape unexpected: {displacement.shape}, {blocked.shape}")
    if not displacement.flags.writeable or not blocked.flags.writeable:
        raise RuntimeError("native displacement or blocked-DOF view is not writable")
    blocked[:, :] = False
    blocked[bottom, 2] = True
    blocked[top, 2] = True
    blocked[pin_x, 0:2] = True
    blocked[pin_y, 1] = True
    expected_blocked = int(case["boundary_conditions"]["expected_blocked_components_after_insertion"])
    if int(np.count_nonzero(blocked)) != expected_blocked:
        raise RuntimeError(f"blocked component count {np.count_nonzero(blocked)} != {expected_blocked}")
    interface_native_nodes = np.unique(native_cohesive)
    interior_source_id = int(bc["geometric_interior_input_node_id_zero_based"])
    expected_interior_xyz = np.asarray(bc["geometric_interior_input_node_coordinate_m"], dtype=np.float64)
    if not np.array_equal(input_mesh["nodes"][interior_source_id], expected_interior_xyz):
        raise RuntimeError("declared interior bulk node exact input coordinate check failed")
    interior_native = np.asarray(sorted(input_to_native[interior_source_id]), dtype=np.int64)
    expected_interface_free = int(bc["expected_free_cohesive_dof_components"])
    expected_interior_free = int(bc["expected_free_geometric_interior_bulk_dof_components"])
    if (len(interface_native_nodes) * 3 - int(np.count_nonzero(blocked[interface_native_nodes])) != expected_interface_free
            or len(interior_native) * 3 - int(np.count_nonzero(blocked[interior_native])) != expected_interior_free):
        raise RuntimeError("free cohesive or geometric-interior bulk component count failed")
    expected_free = int(bc["expected_free_dof_components"])
    if blocked.size - int(np.count_nonzero(blocked)) != expected_free:
        raise RuntimeError("total free native DOF component count failed")

    displacement[:, :] = 0.0
    initial_blocked = blocked.copy()
    initial_state_path = OUT / "initial_state.npz"
    if initial_state_path.exists():
        raise RuntimeError("initial_state.npz exists; refusing to overwrite a prior run")
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        initial_state_path,
        input_nodes=input_nodes_from_api,
        input_bulk_tets=input_tets_from_api,
        native_rest_nodes=native_rest,
        native_bulk_tets=native_tets,
        native_cohesive_connectivity=native_cohesive,
        native_mass_kg=mass0,
        native_blocked_dofs=initial_blocked,
        top_native_node_ids=top,
        bottom_native_node_ids=bottom,
        pin_x_native_node_ids=pin_x,
        pin_y_native_node_ids=pin_y,
        interface_native_node_ids=interface_native_nodes,
        interior_native_node_ids=interior_native,
        interface_positive_rest_area_per_qp_m2=np.repeat(
            np.asarray([row["positive_rest_area_m2"] / 3.0 for row in lineage["parent_facets"]]), 3
        ),
    )
    write_json(OUT / "native_lineage.json", lineage)
    (OUT / "raw_history").mkdir(exist_ok=False)
    history = RawHistory(OUT / "raw_history", int(case["loading"]["total_steps"]) + 1)
    trace_path = OUT / "trace.jsonl"
    if trace_path.exists():
        raise RuntimeError("trace.jsonl exists; refusing to append to a prior run")
    started = time.monotonic()
    work_cumulative = 0.0
    actuator_work_cumulative = 0.0
    previous_opening: np.ndarray | None = None
    previous_traction: np.ndarray | None = None
    previous_actuator_reaction: float | None = None
    previous_command = 0.0
    previous_fields = native_fields(model, aka)
    previous_energies = native_energies(model)
    if not all(np.isfinite(a).all() for a in previous_fields.values()) or not np.isfinite(previous_energies).all():
        raise RuntimeError("initial native fields or energies are non-finite")
    disp0, force0, accel0, ext0, mass_now0 = force_terms(model)
    with np.load(initial_state_path, allow_pickle=False) as initial_archive:
        qpareas = initial_archive["interface_positive_rest_area_per_qp_m2"].copy()
    dot0 = np.einsum("qi,qi->q", previous_fields["opening"], previous_fields["tractions"])
    error0 = boundary_error(disp0, top, bottom, pin_x, pin_y, 0.0)
    top_internal0 = float(force0[top, 2].sum())
    top_inertial0 = float((mass_now0[top, 2] * accel0[top, 2]).sum())
    top_external0 = float(ext0[top, 2].sum())
    actuator_reaction0 = top_inertial0 - top_internal0 - top_external0
    row0 = {
        "opening": previous_fields["opening"], "tractions": previous_fields["tractions"],
        "damage": previous_fields["damage"], "delta_max": previous_fields["delta max"],
        "traction_dot_opening_Pa_m_raw": dot0,
        "interface_work_increment_J_by_qp_raw": np.zeros(24),
        "energies_J_raw": previous_energies,
        "displacement_m": disp0, "internal_force_N_raw": force0,
        "acceleration_m_s2_raw": accel0, "external_force_N_raw": ext0,
        "mass_kg_raw": mass_now0,
        "boundary_readback_error_m": error0,
        "native_work_cumulative_J_raw": 0.0,
        "native_energy_total_J_raw": previous_energies[2] + previous_energies[3],
        "native_work_vs_total_residual_J_raw": -(previous_energies[2] + previous_energies[3]),
        "actuator_reaction_estimate_N_raw": actuator_reaction0,
        "top_internal_force_z_sum_N_raw": top_internal0,
        "top_inertial_force_z_sum_N_raw": top_inertial0,
        "top_external_force_z_sum_N_raw": top_external0,
        "actuator_work_cumulative_J_raw": 0.0,
        "energy_balance_residual_J_raw": -(float(previous_energies[0] + previous_energies[1]
                                                  + previous_energies[2] + previous_energies[3]
                                                  + previous_energies[4])),
        "commanded_top_displacement_m": 0.0,
        "topology_counts": np.asarray(counts0, dtype=np.int64),
    }
    history.write(0, row0)
    HISTORY_FRAMES_SAVED = 1
    with trace_path.open("x", encoding="utf-8") as trace:
        trace.write(json.dumps({
            "history_index": 0, "solver_step": 0, "initialized_state": True,
            "commanded_top_displacement_m": 0.0,
            "counts": {"nodes": counts0[0], "bulk_tets": counts0[1], "cohesive_3d6": counts0[2]},
            "blocked_dof_components": int(np.count_nonzero(initial_blocked)),
            "free_dof_components": int(initial_blocked.size - np.count_nonzero(initial_blocked)),
            "boundary_readback_error_m": error0,
            "energies_J_raw": energy_record(previous_energies),
            "native_traction_dot_opening_sum_Pa_m_raw": float(dot0.sum()),
            "native_interface_work_cumulative_J_raw": 0.0,
            "native_energy_total_J_raw": float(previous_energies[2] + previous_energies[3]),
            "native_work_vs_total_residual_J_raw": float(row0["native_work_vs_total_residual_J_raw"]),
            "top_internal_force_z_sum_N_raw": top_internal0,
            "top_inertial_force_z_sum_N_raw": top_inertial0,
            "top_external_force_z_sum_N_raw": top_external0,
        }, allow_nan=False) + "\n")
        trace.flush()
        os.fsync(trace.fileno())

        previous_opening = previous_fields["opening"].copy()
        previous_traction = previous_fields["tractions"].copy()
        previous_actuator_reaction = actuator_reaction0
        previous_command = 0.0
        model.setTimeStep(dt)
        for step in range(1, int(case["loading"]["total_steps"]) + 1):
            commanded = commanded_opening(case, step)
            # Set only the boundary values; all remaining bulk and cohesive DOFs evolve freely.
            displacement_live = np.asarray(model.getDisplacement())
            blocked_live = np.asarray(model.getBlockedDOFs())
            if blocked_live.shape != (54, 3) or not np.array_equal(blocked_live, initial_blocked):
                raise RuntimeError(f"native blocked-DOF mask changed before step {step}")
            displacement_live[top, 2] = commanded
            displacement_live[bottom, 2] = 0.0
            displacement_live[pin_x, 0:2] = 0.0
            displacement_live[pin_y, 1] = 0.0
            pre_disp = displacement_live.copy()
            pre_error = boundary_error(pre_disp, top, bottom, pin_x, pin_y, commanded)
            if pre_error > 1e-12:
                raise RuntimeError(f"pre-solve boundary readback failed at step {step}: {pre_error}")

            model.solveStep()
            SOLVER_STEPS_COMPLETED += 1
            fields = native_fields(model, aka)
            energies = native_energies(model)
            displacement_after, internal, acceleration, external, mass_now = force_terms(model)
            rest_now = np.asarray(mesh.getNodes()).copy()
            tets_now = np.asarray(mesh.getConnectivity(aka._tetrahedron_4)).copy()
            coh_now = np.asarray(mesh.getConnectivity(aka._cohesive_3d_6)).copy()
            counts = (len(rest_now), len(tets_now), len(coh_now))
            if counts != counts0 or not np.array_equal(rest_now, native_rest) \
                    or not np.array_equal(tets_now, native_tets) or not np.array_equal(coh_now, native_cohesive):
                raise RuntimeError(f"native topology/rest state changed at step {step}: {counts}")
            blocked_after = np.asarray(model.getBlockedDOFs()).copy()
            if not np.array_equal(blocked_after, initial_blocked):
                raise RuntimeError(f"native blocked-DOF mask changed at step {step}")
            post_error = boundary_error(displacement_after, top, bottom, pin_x, pin_y, commanded)
            if post_error > 1e-12:
                raise RuntimeError(f"post-solve boundary readback failed at step {step}: {post_error}")
            if not all(np.isfinite(value).all() for value in fields.values()) or not np.isfinite(energies).all():
                raise RuntimeError(f"non-finite native cohesive field or energy at step {step}")

            qp_work_increment = quadrature_signed_work(
                previous_opening, previous_traction,
                fields["opening"], fields["tractions"], qpareas,
            )
            work_cumulative += float(qp_work_increment.sum())
            dot_now = np.einsum("qi,qi->q", fields["opening"], fields["tractions"])
            top_internal = float(internal[top, 2].sum())
            top_inertial = float((mass_now[top, 2] * acceleration[top, 2]).sum())
            top_external = float(external[top, 2].sum())
            actuator_reaction = top_inertial - top_internal - top_external
            if previous_actuator_reaction is not None:
                actuator_work_cumulative += 0.5 * (previous_actuator_reaction + actuator_reaction) * (commanded - previous_command)
            native_total = float(energies[2] + energies[3])
            native_work_residual = work_cumulative - native_total
            balance_residual = actuator_work_cumulative - float(energies.sum())
            frame = step
            row = {
                "opening": fields["opening"], "tractions": fields["tractions"],
                "damage": fields["damage"], "delta_max": fields["delta max"],
                "traction_dot_opening_Pa_m_raw": dot_now,
                "interface_work_increment_J_by_qp_raw": qp_work_increment,
                "energies_J_raw": energies,
                "displacement_m": displacement_after,
                "internal_force_N_raw": internal,
                "acceleration_m_s2_raw": acceleration,
                "external_force_N_raw": external,
                "mass_kg_raw": mass_now,
                "boundary_readback_error_m": post_error,
                "native_work_cumulative_J_raw": work_cumulative,
                "native_energy_total_J_raw": native_total,
                "native_work_vs_total_residual_J_raw": native_work_residual,
                "actuator_reaction_estimate_N_raw": actuator_reaction,
                "top_internal_force_z_sum_N_raw": top_internal,
                "top_inertial_force_z_sum_N_raw": top_inertial,
                "top_external_force_z_sum_N_raw": top_external,
                "actuator_work_cumulative_J_raw": actuator_work_cumulative,
                "energy_balance_residual_J_raw": balance_residual,
                "commanded_top_displacement_m": commanded,
                "topology_counts": np.asarray(counts, dtype=np.int64),
            }
            history.write(frame, row)
            HISTORY_FRAMES_SAVED += 1
            trace_row = {
                "history_index": frame, "solver_step": step, "initialized_state": False,
                "commanded_top_displacement_m": commanded,
                "pre_solve_boundary_readback_error_m": pre_error,
                "post_solve_boundary_readback_error_m": post_error,
                "native_blocked_dof_mask_unchanged": True,
                "counts": {"nodes": counts[0], "bulk_tets": counts[1], "cohesive_3d6": counts[2]},
                "energies_J_raw": energy_record(energies),
                "native_traction_dot_opening_sum_Pa_m_raw": float(dot_now.sum()),
                "interface_work_increment_J_raw": float(qp_work_increment.sum()),
                "native_interface_work_cumulative_J_raw": work_cumulative,
                "native_energy_total_J_raw": native_total,
                "native_work_vs_total_residual_J_raw": native_work_residual,
                "actuator_reaction_estimate_N_raw": actuator_reaction,
                "top_internal_force_z_sum_N_raw": top_internal,
                "top_inertial_force_z_sum_N_raw": top_inertial,
                "top_external_force_z_sum_N_raw": top_external,
                "actuator_work_cumulative_J_raw": actuator_work_cumulative,
                "energy_balance_residual_J_raw": balance_residual,
                "damage_min_raw": float(fields["damage"].min()),
                "damage_max_raw": float(fields["damage"].max()),
                "delta_max_min_m_raw": float(fields["delta max"].min()),
                "delta_max_max_m_raw": float(fields["delta max"].max()),
            }
            trace.write(json.dumps(trace_row, allow_nan=False) + "\n")
            trace.flush()
            os.fsync(trace.fileno())
            previous_opening = fields["opening"].copy()
            previous_traction = fields["tractions"].copy()
            previous_actuator_reaction = actuator_reaction
            previous_command = commanded
            if step % 100 == 0:
                OBS["last_solver_step_checkpoint"] = step
            if time.monotonic() - started > float(case["loading"]["wall_time_cap_s"]):
                raise TimeoutError(f"wall-time cap reached after step {step}")

    required_frames = int(case["loading"]["total_steps"]) + 1
    if HISTORY_FRAMES_SAVED != required_frames:
        raise RuntimeError(f"history frame count {HISTORY_FRAMES_SAVED} != {required_frames}")
    return {
        "status": "RUN_COMPLETE_RAW_HISTORY_CLOSED",
        "solver_steps_completed": SOLVER_STEPS_COMPLETED,
        "history_frames_saved": HISTORY_FRAMES_SAVED,
        "expected_history_frames": required_frames,
        "counts": {"nodes": 54, "bulk_tets": 96, "cohesive_3d6": 8},
        "blocked_dof_components": expected_blocked,
        "free_dof_components": expected_free,
        "free_cohesive_dof_components": expected_interface_free,
        "free_geometric_interior_bulk_dof_components": expected_interior_free,
        "stability_contract": OBS["stability_contract"],
        "bulk_stable_step_s": bulk_stable_dt,
        "solver_dt_s": dt,
        "solver_dt_over_bulk_stable_step": dt / bulk_stable_dt,
        "solver_dt_over_derived_total_stable_step": dt / float(total_stable_dt),
        "native_lineage": "native_lineage.json",
        "initial_state": "initial_state.npz",
        "raw_history_directory": "raw_history",
        "trace": "trace.jsonl",
        "solver_attempts": 1,
        "native_dissipation_sign_status": "RAW_RETAINED; analyzer determines diagnostic gate",
        "classification": case["classification"],
    }


def main() -> int:
    if OUT.exists() and any(OUT.iterdir()):
        raise SystemExit(f"output exists and is not empty; refusing rerun: {OUT}")
    OUT.mkdir(parents=True, exist_ok=True)
    case = json.loads(CASE_PATH.read_text(encoding="utf-8"))
    start = {
        "argv": sys.argv,
        "case_sha256": sha256(CASE_PATH),
        "mesh_sha256": sha256(MESH_PATH),
        "material_sha256": sha256(MATERIAL_PATH),
        "runner_sha256": sha256(Path(__file__)),
        "contract_sha256": sha256(EXP / "free_patch_contract.py"),
        "analyzer_sha256": sha256(EXP / "analyze_native_free_patch.py"),
        "test_sha256": sha256(EXP / "test_source_contract.py"),
        "launcher_sha256": sha256(EXP / "launch_free_patch.sh"),
        "total_stability_status": case["stability"]["total_stability_status"],
        "lambda_max_upper_bound_per_s2": case["stability"]["lambda_max_upper_bound_per_s2"],
        "python": sys.version.split()[0],
        "solver_attempts": 1,
        "retry_policy": "none",
        "started_unix_time_s": time.time(),
    }
    write_json(OUT / "run_started.json", start)
    result: dict[str, Any] = {}
    try:
        os.chdir(OUT)
        result = run_native()
    except Exception as error:
        result = {
            **OBS,
            "status": "FAIL_NO_RETRY",
            "error": repr(error),
            "traceback": traceback.format_exc(),
            "solver_steps_completed": SOLVER_STEPS_COMPLETED,
            "history_frames_saved": HISTORY_FRAMES_SAVED,
        }
    result.update({"argv": sys.argv, "ended_unix_time_s": time.time(), "solver_attempts": 1})
    write_json(OUT / "run_status.json", result)
    print(json.dumps(json_safe(result), allow_nan=False, sort_keys=True))
    return 0 if result.get("status") == "RUN_COMPLETE_RAW_HISTORY_CLOSED" else 2


if __name__ == "__main__":
    raise SystemExit(main())
