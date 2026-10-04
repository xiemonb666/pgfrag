"""Run one fixed contract case against the isolated Akantu Python wheel."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import sys
import time
import traceback

import akantu as aka
import numpy as np


ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "input"
REVIEW = ROOT.parents[2] / "relicweave-synlab" / "tmp" / "review_20261004" / "implementation" / "cohesive_capability"
MANIFEST_PATH = REVIEW / "manifest.json"
MAX_WALL_S = 300.0
MAX_STEPS = 10_000


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def array_copy(value) -> np.ndarray:
    return np.asarray(value).copy()


def scalar_mass_by_dof_component(mass: np.ndarray) -> list[float]:
    if mass.ndim == 2:
        return [float(mass[:, component].sum()) for component in range(mass.shape[1])]
    return [float(mass.sum())]


class RampTop(aka.DirichletFunctor):
    def __init__(self, amplitude: float, duration: float):
        super().__init__(aka._z)
        self.axis = int(aka._z)
        self.amplitude = amplitude
        self.duration = duration
        self.time = 0.0

    def set_time(self, value: float) -> None:
        self.time = value

    def __call__(self, node, flags, disp, coord):
        flags[self.axis] = True
        disp[self.axis] = self.amplitude * min(max(self.time / self.duration, 0.0), 1.0)


def mesh_element_count(mesh, element_type) -> int:
    return int(mesh.getNbElement(type=element_type))


def mass_summary(model) -> dict:
    values = array_copy(model.getMass())
    return {
        "shape": list(values.shape),
        "raw_sum_over_all_dof_components_kg": float(values.sum()),
        "component_sums_kg": scalar_mass_by_dof_component(values),
        "array": values,
    }


def energy_summary(model) -> dict[str, float]:
    names = ("kinetic", "potential", "dissipated", "reversible", "contact")
    return {name: float(model.getEnergy(name)) for name in names}


def read_interface_triangles(mesh_path: Path) -> list[tuple[int, int, int]]:
    """Read physical group 2 directly from the frozen MSH input, without spatial matching."""
    lines = mesh_path.read_text(encoding="ascii").splitlines()
    start = lines.index("$Elements") + 2
    end = lines.index("$EndElements")
    result = []
    for line in lines[start:end]:
        fields = [int(item) for item in line.split()]
        _element_id, element_type, tag_count = fields[:3]
        if element_type != 2 or tag_count < 1 or fields[3] != 2:
            continue
        result.append(tuple(node_id - 1 for node_id in fields[5:8]))
    return result


def parent_face_index(tetrahedra: np.ndarray, face_key: tuple[int, int, int]) -> list[int]:
    matches = []
    for element_id, tet in enumerate(tetrahedra):
        if len(set(face_key).intersection(int(node) for node in tet)) == 3:
            matches.append(element_id)
    return matches


def construct_native_lineage(
    original_nodes: np.ndarray,
    original_tets: np.ndarray,
    original_interface_faces: list[tuple[int, int, int]],
    final_nodes: np.ndarray,
    final_tets: np.ndarray,
    cohesive_connectivity: np.ndarray,
) -> dict:
    """Trace material ancestry through native bulk element IDs/local vertex ordinals.

    This deliberately does not use geometric proximity. It accepts a map only when
    post-insertion bulk connectivity has the same element count/order and gives one
    conflict-free current-node -> input-material-node map under the native
    element-ID/local-ordinal correspondence.
    """
    n_original = len(original_nodes)
    if len(final_nodes) < n_original or len(final_tets) != len(original_tets):
        return {
            "status": "UNKNOWN",
            "reason": "post-insertion node or bulk-element cardinality is inconsistent with the native pre-insertion snapshot",
            "method": "native bulk element ID + local vertex ordinal only; no coordinate matching",
            "facets": [],
        }

    current_to_material: dict[int, int] = {}
    for bulk_element_id, (source_tet, post_tet) in enumerate(zip(original_tets, final_tets)):
        if len(source_tet) != 4 or len(post_tet) != 4:
            return {
                "status": "UNKNOWN",
                "reason": f"unexpected tetrahedron connectivity width at native bulk element {bulk_element_id}",
                "method": "native bulk element ID + local vertex ordinal only; no coordinate matching",
                "facets": [],
            }
        for local_ordinal, (source_node, current_node) in enumerate(zip(source_tet, post_tet)):
            source_node = int(source_node)
            current_node = int(current_node)
            if current_node < n_original and current_node != source_node:
                return {
                    "status": "UNKNOWN",
                    "reason": f"native input-node ID {current_node} changed ancestry at bulk element {bulk_element_id}, local vertex {local_ordinal}",
                    "method": "native bulk element ID + local vertex ordinal only; no coordinate matching",
                    "facets": [],
                }
            previous = current_to_material.get(current_node)
            if previous is not None and previous != source_node:
                return {
                    "status": "UNKNOWN",
                    "reason": f"current node {current_node} maps to conflicting input material IDs {previous} and {source_node}",
                    "method": "native bulk element ID + local vertex ordinal only; no coordinate matching",
                    "facets": [],
                }
            current_to_material[current_node] = source_node

    used_post_nodes = {int(node) for node in final_tets.flat}
    if not used_post_nodes.issubset(current_to_material):
        return {
            "status": "UNKNOWN",
            "reason": "some post-insertion bulk nodes have no cell-ID/local-ordinal material ancestry",
            "method": "native bulk element ID + local vertex ordinal only; no coordinate matching",
            "facets": [],
        }
    if any(current_to_material.get(node) != node for node in range(n_original)):
        return {
            "status": "UNKNOWN",
            "reason": "one or more original input-node IDs are absent or no longer self-identify in native bulk connectivity",
            "method": "native bulk element ID + local vertex ordinal only; no coordinate matching",
            "facets": [],
        }

    face_to_parent = {
        tuple(sorted(face)): parent_face_index(original_tets, tuple(sorted(face)))
        for face in original_interface_faces
    }
    mapped = []
    side_area_sums_m2 = [0.0, 0.0]
    side_normal_dot_products = []
    for cohesive_id, row in enumerate(cohesive_connectivity):
        if len(row) != 6:
            return {
                "status": "UNKNOWN",
                "reason": f"unexpected cohesive connectivity width {len(row)}",
                "method": "no coordinate matching used",
                "facets": mapped,
            }
        side_a = [int(node) for node in row[:3]]
        side_b = [int(node) for node in row[3:]]
        if any(node not in current_to_material for node in side_a + side_b):
            return {
                "status": "UNKNOWN",
                "reason": "native cohesive connectivity includes a node without cell-ID/local-ordinal material ancestry",
                "method": "native bulk element ID + local vertex ordinal only; no coordinate matching",
                "facets": mapped,
            }
        material_side_a = tuple(sorted(current_to_material[node] for node in side_a))
        material_side_b = tuple(sorted(current_to_material[node] for node in side_b))
        if material_side_a != material_side_b or material_side_a not in face_to_parent:
            return {
                "status": "UNKNOWN",
                "reason": "native cohesive sides do not map through bulk cell identity to the same declared input interface triangle",
                "method": "native bulk element ID + local vertex ordinal only; no coordinate matching",
                "facets": mapped,
            }
        source_face = material_side_a
        parents = face_to_parent[source_face]
        if len(parents) != 2:
            return {
                "status": "UNKNOWN",
                "reason": f"declared interface face has {len(parents)} input bulk parents, expected two",
                "method": "native bulk element ID + local vertex ordinal only; no coordinate matching",
                "facets": mapped,
            }
        parent_by_side = []
        for side in (side_a, side_b):
            side_parents = [
                tet_id
                for tet_id, tet in enumerate(final_tets)
                if len(set(side).intersection(int(node) for node in tet)) == 3
            ]
            parent_by_side.append(side_parents)
        if any(len(side_parents) != 1 for side_parents in parent_by_side):
            return {
                "status": "UNKNOWN",
                "reason": "native post-insertion bulk connectivity does not give one parent tetrahedron per side",
                "method": "native bulk element ID + local vertex ordinal only; no coordinate matching",
                "facets": mapped,
            }
        side_geometry = []
        for side_index, side in enumerate((side_a, side_b)):
            points = final_nodes[np.asarray(side, dtype=np.int64)]
            cross = np.cross(points[1] - points[0], points[2] - points[0])
            norm = float(np.linalg.norm(cross))
            if norm == 0.0:
                return {
                    "status": "UNKNOWN",
                    "reason": f"native cohesive element {cohesive_id} has a zero-area side",
                    "method": "native bulk element ID + local vertex ordinal only; no coordinate matching",
                    "facets": mapped,
                }
            area = 0.5 * norm
            normal = cross / norm
            side_area_sums_m2[side_index] += area
            side_geometry.append({"area_m2": area, "unit_normal": [float(v) for v in normal]})
        normal_dot = float(np.dot(side_geometry[0]["unit_normal"], side_geometry[1]["unit_normal"]))
        side_normal_dot_products.append(normal_dot)
        mapped.append(
            {
                "cohesive_element_id": cohesive_id,
                "parent_input_face_node_ids": list(source_face),
                "side_a_node_ids": side_a,
                "side_b_node_ids": side_b,
                "side_a_material_node_ids": list(current_to_material[node] for node in side_a),
                "side_b_material_node_ids": list(current_to_material[node] for node in side_b),
                "current_node_to_input_material_node_id": {
                    str(node): current_to_material[node]
                    for node in side_a + side_b
                },
                "input_parent_tetrahedron_ids": parents,
                "post_insert_parent_tetrahedron_ids_by_cohesive_side": parent_by_side,
                "rest_side_geometry": side_geometry,
                "rest_side_normal_dot_product": normal_dot,
            }
        )
    return {
        "status": "EXACT_FOR_CONTRACT_MESH_FROM_NATIVE_CONNECTIVITY",
        "reason": "post-insertion bulk tetrahedra retained the pre-insertion native element count/order and yielded a conflict-free current-node -> input-material-node map under identical native element IDs and local vertex ordinals; each cohesive side maps to the same declared input face and one bulk parent per side",
        "method": "native bulk element IDs + local vertex ordinals + cohesive connectivity only; no coordinate-nearest matching",
        "bulk_elements_before_after": [len(original_tets), len(final_tets)],
        "current_bulk_nodes_mapped": len(used_post_nodes),
        "current_nodes_with_input_material_ancestry": len(current_to_material),
        "current_node_to_input_material_node_id": {str(k): v for k, v in sorted(current_to_material.items())},
        "cohesive_rest_side_area_sums_m2": side_area_sums_m2,
        "cohesive_rest_side_normal_dot_products": side_normal_dot_products,
        "facets": mapped,
    }


def main() -> int:
    if len(sys.argv) != 2 or sys.argv[1] not in {"zero_load", "modeI_dt", "modeI_dt2"}:
        raise SystemExit("usage: run_contract.py {zero_load|modeI_dt|modeI_dt2}")
    case_name = sys.argv[1]
    if not MANIFEST_PATH.exists():
        raise RuntimeError("Frozen manifest is required before any simulation case starts")
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    case = manifest["cases"][case_name]
    output = REVIEW / case_name
    output.mkdir(parents=True, exist_ok=True)
    start_utc = datetime.now(timezone.utc).isoformat()
    start_mono = time.monotonic()
    status = {
        "case": case_name,
        "start_utc": start_utc,
        "argv": sys.argv,
        "cwd": os.getcwd(),
        "python": sys.executable,
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "manifest_sha256": sha256(MANIFEST_PATH),
        "input_hashes_valid": True,
        "status": "RUNNING",
    }
    for relative, expected in manifest["input_sha256"].items():
        actual = sha256(ROOT / relative)
        if actual.lower() != expected.lower():
            status["input_hashes_valid"] = False
            status["input_hash_mismatch"] = {"path": relative, "expected": expected, "actual": actual}
            (output / "status.json").write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
            raise RuntimeError(f"Frozen input hash mismatch: {relative}")

    try:
        aka.parseInput(str(INPUT / "material.dat"))
        mesh = aka.Mesh(3)
        mesh.read(str(INPUT / "block_2x2x4mm.msh"))
        original_nodes = array_copy(mesh.getNodes())
        original_tets = array_copy(mesh.getConnectivity(aka._tetrahedron_4))
        original_interface_faces = read_interface_triangles(INPUT / "block_2x2x4mm.msh")
        model = aka.SolidMechanicsModelCohesive(mesh)
        model.initFull(_analysis_method=aka._explicit_lumped_mass, _is_extrinsic=True)
        model.assembleMassLumped()

        inserter = model.getElementInserter()
        inserter.addPhysicalSurface("interface")
        model.updateAutomaticInsertion()
        check_map = inserter.getCheckFacets()
        eligible = int(np.count_nonzero(check_map(aka._triangle_3)))
        eligible_total = int(check_map(aka._triangle_3).size)
        insertion_map = inserter.getInsertionFacets()
        already_inserted = int(np.count_nonzero(insertion_map(aka._triangle_3)))
        if eligible != manifest["expected_insertion_facets"] or already_inserted != 0:
            raise RuntimeError(
                f"pre-step interface selection mismatch: {eligible}/{eligible_total} check facets, "
                f"{already_inserted} already inserted"
            )

        mass_before = mass_summary(model)
        initial_energies = energy_summary(model)
        if case_name == "zero_load":
            model.applyBC(aka.FixedValue(0.0, aka._z), "bottom")
            model.applyBC(aka.FixedValue(0.0, aka._z), "top")
            ramp = None
        else:
            model.applyBC(aka.FixedValue(0.0, aka._z), "bottom")
            ramp = RampTop(
                float(manifest["loading"]["top_displacement_m"]),
                float(manifest["loading"]["ramp_duration_s"]),
            )
        top_node_ids = np.flatnonzero(np.abs(original_nodes[:, 2] - 4.0e-3) < 1.0e-15)
        if len(top_node_ids) != 9:
            raise RuntimeError(f"expected 9 top nodes, found {len(top_node_ids)}")

        dt = float(case["dt_s"])
        steps = int(case["steps"])
        if steps > MAX_STEPS:
            steps = MAX_STEPS
            status["status"] = "NOT_REACHED_STEP_LIMIT"
        model.setTimeStep(dt)

        trace_path = output / "trace.csv"
        field_names = (
            "step",
            "time_s",
            "imposed_top_displacement_m",
            "top_internal_reaction_raw_N",
            "actuator_force_estimate_N",
            "external_work_J",
            "kinetic_J",
            "bulk_potential_J",
            "cohesive_dissipated_J",
            "cohesive_reversible_J",
            "cohesive_contact_energy_J",
            "energy_ledger_residual_J",
            "mass_component_0_kg",
            "mass_component_1_kg",
            "mass_component_2_kg",
            "mesh_nodes",
            "bulk_tetrahedra",
            "cohesive_elements",
            "cohesive_insertions_this_step",
            "cohesive_insertions_cumulative",
        )
        rows = []
        events = []
        current_work = 0.0
        previous_force = 0.0
        previous_displacement = 0.0
        total_native_insertions = 0
        total_api_insertions = 0
        inserted_count_seen = mesh_element_count(mesh, aka._cohesive_3d_6)
        event_first_insert_step = None
        event_first_insert_time = None
        full_debond_step = None
        full_debond_time = None
        expected_dissipation = float(manifest["expected_interface_dissipation_J"])
        mass_before_components = mass_before["component_sums_kg"]

        def append_row(step: int, time_value: float, native_inserted_count: int = 0) -> dict:
            nonlocal current_work, previous_force, previous_displacement, total_native_insertions
            displacement = float(manifest["loading"]["top_displacement_m"]) * min(
                max(time_value / float(manifest["loading"]["ramp_duration_s"]), 0.0), 1.0
            ) if ramp is not None else 0.0
            force_raw = 0.0
            if ramp is not None:
                internal = array_copy(model.getInternalForce())
                force_raw = float(internal[top_node_ids, 2].sum())
            actuator_force = -force_raw
            current_work += 0.5 * (previous_force + actuator_force) * (
                displacement - previous_displacement
            )
            previous_force = actuator_force
            previous_displacement = displacement
            energies = energy_summary(model)
            mass_now = mass_summary(model)["component_sums_kg"]
            total_native_insertions += native_inserted_count
            storage_and_loss = (
                energies["kinetic"]
                + energies["potential"]
                + energies["reversible"]
                + energies["dissipated"]
                + energies["contact"]
            )
            row = {
                "step": step,
                "time_s": time_value,
                "imposed_top_displacement_m": displacement,
                "top_internal_reaction_raw_N": force_raw,
                "actuator_force_estimate_N": actuator_force,
                "external_work_J": current_work,
                "kinetic_J": energies["kinetic"],
                "bulk_potential_J": energies["potential"],
                "cohesive_dissipated_J": energies["dissipated"],
                "cohesive_reversible_J": energies["reversible"],
                "cohesive_contact_energy_J": energies["contact"],
                "energy_ledger_residual_J": current_work - storage_and_loss,
                "mass_component_0_kg": mass_now[0],
                "mass_component_1_kg": mass_now[1],
                "mass_component_2_kg": mass_now[2] if len(mass_now) > 2 else math.nan,
                "mesh_nodes": mesh.getNbNodes(),
                "bulk_tetrahedra": mesh_element_count(mesh, aka._tetrahedron_4),
                "cohesive_elements": mesh_element_count(mesh, aka._cohesive_3d_6),
                "cohesive_insertions_this_step": native_inserted_count,
                "cohesive_insertions_cumulative": total_native_insertions,
            }
            return row

        initial_row = append_row(0, 0.0)
        rows.append(initial_row)
        for step in range(1, steps + 1):
            current_time = step * dt
            if ramp is not None:
                ramp.set_time(current_time)
                model.applyBC(ramp, "top")
            new_count = int(model.checkCohesiveStress())
            current_cohesive_count = mesh_element_count(mesh, aka._cohesive_3d_6)
            if current_cohesive_count < inserted_count_seen:
                raise RuntimeError("cohesive connectivity count decreased during the run")
            actual_added = current_cohesive_count - inserted_count_seen
            if new_count != actual_added:
                events.append(
                    {
                        "kind": "insertion_count_api_disagreement",
                        "step": step,
                        "time_s": current_time,
                        "checkCohesiveStress_return": new_count,
                        "cohesive_connectivity_count_delta": actual_added,
                    }
                )
            total_api_insertions += new_count
            if new_count > 0 or actual_added > 0:
                if event_first_insert_step is None:
                    event_first_insert_step = step
                    event_first_insert_time = current_time
                events.append(
                    {
                        "kind": "dynamic_cohesive_insertion",
                        "step": step,
                        "time_s": current_time,
                        "checkCohesiveStress_return": new_count,
                        "cohesive_elements_after_check": current_cohesive_count,
                        "mesh_nodes_after_check": mesh.getNbNodes(),
                    }
                )
            model.solveStep()
            inserted_count_seen = current_cohesive_count
            row = append_row(step, current_time, actual_added)
            rows.append(row)
            if full_debond_step is None and expected_dissipation > 0:
                if row["cohesive_dissipated_J"] >= expected_dissipation * 0.999:
                    full_debond_step = step
                    full_debond_time = current_time
            if time.monotonic() - start_mono > MAX_WALL_S:
                status["status"] = "NOT_REACHED_WALL_LIMIT"
                events.append(
                    {"kind": "wall_limit", "step": step, "time_s": current_time, "limit_s": MAX_WALL_S}
                )
                break

        with trace_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=field_names)
            writer.writeheader()
            writer.writerows(rows)

        # Re-fetch and copy native arrays after insertion; do not reuse pre-insertion views.
        final_nodes = array_copy(mesh.getNodes())
        final_tets = array_copy(mesh.getConnectivity(aka._tetrahedron_4))
        final_cohesive_count = mesh_element_count(mesh, aka._cohesive_3d_6)
        final_cohesive_connectivity = (
            array_copy(mesh.getConnectivity(aka._cohesive_3d_6))
            if final_cohesive_count
            else np.empty((0, 6), dtype=np.int64)
        )
        final_mass = mass_summary(model)
        final_energy = energy_summary(model)
        lineage = construct_native_lineage(
            original_nodes,
            original_tets,
            original_interface_faces,
            final_nodes,
            final_tets,
            final_cohesive_connectivity,
        )

        np.savez_compressed(
            output / "native_topology.npz",
            original_rest_nodes=original_nodes,
            original_bulk_tetra_connectivity=original_tets,
            original_interface_triangle_connectivity=np.asarray(original_interface_faces, dtype=np.int64),
            final_rest_nodes=final_nodes,
            final_bulk_tetra_connectivity=final_tets,
            final_cohesive_connectivity=final_cohesive_connectivity,
            final_displacement=array_copy(model.getDisplacement()),
            final_velocity=array_copy(model.getVelocity()),
            mass_before=mass_before["array"],
            mass_after=final_mass["array"],
        )
        (output / "events.json").write_text(json.dumps(events, indent=2) + "\n", encoding="utf-8")
        (output / "lineage.json").write_text(json.dumps(lineage, indent=2) + "\n", encoding="utf-8")

        fragment_summary = {"api_present": hasattr(aka, "FragmentManager"), "status": "NOT_RUN"}
        if hasattr(aka, "FragmentManager"):
            try:
                fragment_manager = aka.FragmentManager(model)
                fragment_manager.buildFragments(1.0)
                fragment_summary = {
                    "api_present": True,
                    "status": "READ",
                    "fragment_count": int(fragment_manager.getNbFragment()),
                    "fragment_group_names": [group.getName() for group in fragment_manager.elementGroups()],
                }
            except Exception as error:
                fragment_summary = {
                    "api_present": True,
                    "status": "ERROR",
                    "error": repr(error),
                }

        target_error = (
            abs(final_energy["dissipated"] - expected_dissipation) / expected_dissipation
            if expected_dissipation > 0
            else None
        )
        mass_diff = [after - before for before, after in zip(mass_before_components, final_mass["component_sums_kg"])]
        mass_conserved = all(abs(delta) <= 1e-14 * max(abs(before), 1e-300) for delta, before in zip(mass_diff, mass_before_components))
        if status["status"] == "RUNNING":
            status["status"] = "COMPLETED" if len(rows) == int(case["steps"]) + 1 else "NOT_REACHED"
        status.update(
            {
                "end_utc": datetime.now(timezone.utc).isoformat(),
                "elapsed_s": time.monotonic() - start_mono,
                "stable_time_step_s": float(model.getStableTimeStep()),
                "time_step_s": dt,
                "steps_requested": int(case["steps"]),
                "steps_completed": len(rows) - 1,
                "endtime_target_s": float(manifest["loading"]["endtime_s"]),
                "endtime_reached_s": (len(rows) - 1) * dt,
                "endtime_reached": len(rows) - 1 == int(case["steps"]),
                "inserted_interface_facets_expected": int(manifest["expected_insertion_facets"]),
                "eligible_facets_checked": eligible,
                "eligible_facet_slots_total": eligible_total,
                "dynamic_insertions_native_connectivity": total_native_insertions,
                "dynamic_insertions_api_returned": total_api_insertions,
                "cohesive_elements_final": final_cohesive_count,
                "mesh_nodes_before_after": [len(original_nodes), len(final_nodes)],
                "bulk_tetrahedra_before_after": [len(original_tets), len(final_tets)],
                "first_insertion_step": event_first_insert_step,
                "first_insertion_time_s": event_first_insert_time,
                "complete_debond_step_by_dissipation_target": full_debond_step,
                "complete_debond_time_s_by_dissipation_target": full_debond_time,
                "dissipated_energy_target_J": expected_dissipation,
                "dissipated_energy_final_J": final_energy["dissipated"],
                "dissipated_energy_relative_error": target_error,
                "mass_component_sums_before_kg": mass_before_components,
                "mass_component_sums_after_kg": final_mass["component_sums_kg"],
                "mass_component_sum_deltas_kg": mass_diff,
                "mass_conserved_by_component": mass_conserved,
                "initial_energies_J": initial_energies,
                "final_energies_J": final_energy,
                "external_work_final_J": current_work,
                "energy_ledger_residual_final_J": rows[-1]["energy_ledger_residual_J"],
                "energy_ledger_relative_residual": (
                    rows[-1]["energy_ledger_residual_J"] / current_work if current_work else None
                ),
                "lineage_status": lineage["status"],
                "native_connectivity_shape": list(final_cohesive_connectivity.shape),
                "native_parent_mapping_method": lineage["method"],
                "fragment_summary": fragment_summary,
                "contact_status": "NOT_RUN",
                "roughness_status": "NOT_VALIDATED",
                "material_status": "UNVALIDATED_TOY",
                "source_file_sha256": sha256(Path(__file__)),
                "trace_csv_sha256": sha256(trace_path),
                "native_topology_sha256": sha256(output / "native_topology.npz"),
            }
        )
        (output / "status.json").write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
        return 0
    except Exception as error:
        status.update(
            {
                "end_utc": datetime.now(timezone.utc).isoformat(),
                "elapsed_s": time.monotonic() - start_mono,
                "status": "FAILED",
                "exception": repr(error),
                "traceback": traceback.format_exc(),
            }
        )
        (output / "status.json").write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
        raise


if __name__ == "__main__":
    raise SystemExit(main())
