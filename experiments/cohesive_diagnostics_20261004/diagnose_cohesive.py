"""Read-only geometry, lineage, kinematics, and insertion-energy audit."""
import csv
import hashlib
import json
import math
import re
from collections import Counter
from itertools import combinations
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
INPUT = ROOT / "pgfrag/experiments/cohesive_capability/input"
RUNNER = ROOT / "pgfrag/experiments/cohesive_capability/run_contract.py"
DATA = ROOT / "relicweave-synlab/tmp/review_20261004/implementation/cohesive_capability"
OUT = ROOT / "relicweave-synlab/tmp/review_20261004/implementation/cohesive_diagnostics_20261004"
SOURCE = ROOT / "third_party/akantu-runtime/source-v5.0.7"
CASES = ("zero_load", "modeI_dt", "modeI_dt2")
FACE_ORDINALS = ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))
PLANE_Z = 0.002


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_msh(path):
    lines = path.read_text(encoding="ascii").splitlines()
    physical_names = {}
    nodes = {}
    elements = []
    i = 0
    while i < len(lines):
        marker = lines[i]
        if marker == "$PhysicalNames":
            count = int(lines[i + 1])
            for line in lines[i + 2 : i + 2 + count]:
                dim, tag, name = line.split(maxsplit=2)
                physical_names[(int(dim), int(tag))] = name.strip('"')
            i += count + 3
        elif marker == "$Nodes":
            count = int(lines[i + 1])
            for line in lines[i + 2 : i + 2 + count]:
                values = line.split()
                nodes[int(values[0])] = np.asarray([float(v) for v in values[1:4]], dtype=np.float64)
            i += count + 3
        elif marker == "$Elements":
            count = int(lines[i + 1])
            for line in lines[i + 2 : i + 2 + count]:
                values = [int(v) for v in line.split()]
                element_id, element_type, n_tags = values[:3]
                tags = values[3 : 3 + n_tags]
                node_ids = values[3 + n_tags :]
                elements.append({"id": element_id, "type": element_type,
                                 "physical": tags[0] if tags else 0,
                                 "nodes_1based": node_ids,
                                 "nodes_0based": [v - 1 for v in node_ids]})
            i += count + 3
        else:
            i += 1
    ordered_ids = sorted(nodes)
    node_array = np.asarray([nodes[node_id] for node_id in ordered_ids], dtype=np.float64)
    tets = np.asarray([e["nodes_0based"] for e in elements if e["type"] == 4], dtype=np.int64)
    interface_triangles = np.asarray([e["nodes_0based"] for e in elements if e["type"] == 2 and physical_names.get((2, e["physical"])) == "interface"], dtype=np.int64)
    surfaces = {
        physical_names.get((2, e["physical"]), f"physical_{e['physical']}"): e["nodes_0based"]
        for e in elements if e["type"] == 2
    }
    return {"nodes_by_id": nodes, "ordered_ids": ordered_ids, "nodes": node_array,
            "elements": elements, "tets": tets, "interface_triangles": interface_triangles,
            "physical_names": physical_names, "surface_elements": surfaces}


def signed_tet_volumes(nodes, tets):
    p = nodes[tets]
    return np.linalg.det(np.stack((p[:, 1] - p[:, 0], p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]), axis=2)) / 6.0


def volume_summary(volumes):
    volumes = np.asarray(volumes, dtype=np.float64)
    return {
        "tetrahedra": int(volumes.size),
        "negative_count": int(np.count_nonzero(volumes < 0.0)),
        "zero_count": int(np.count_nonzero(volumes == 0.0)),
        "positive_count": int(np.count_nonzero(volumes > 0.0)),
        "minimum_signed_m3": float(volumes.min()) if volumes.size else None,
        "maximum_signed_m3": float(volumes.max()) if volumes.size else None,
        "sum_signed_m3": float(volumes.sum()),
        "sum_absolute_m3_reported_separately": float(np.abs(volumes).sum()),
    }


def face_counts(tets):
    counts = Counter()
    for tet in tets:
        for local in FACE_ORDINALS:
            counts[tuple(sorted(int(tet[j]) for j in local))] += 1
    return counts


def face_summary(counts):
    return {
        "face_incidence_slots": int(sum(counts.values())),
        "unique_faces": len(counts),
        "boundary_faces_count1": sum(value == 1 for value in counts.values()),
        "interior_faces_count2": sum(value == 2 for value in counts.values()),
        "nonmanifold_faces_gt2": sum(value > 2 for value in counts.values()),
        "zero_incidence_faces": sum(value == 0 for value in counts.values()),
        "maximum_incidence": max(counts.values(), default=0),
    }


def rows_by_step(path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    return {int(row["step"]): row for row in rows}


def f(row, key):
    return float(row[key])


ENERGY_FIELDS = (
    "kinetic_J", "bulk_potential_J", "cohesive_dissipated_J",
    "cohesive_reversible_J", "cohesive_contact_energy_J",
    "external_work_J", "energy_ledger_residual_J",
)
REACTION_FIELDS = ("imposed_top_displacement_m", "top_internal_reaction_raw_N", "actuator_force_estimate_N")


def trace_event_audit(rows, events, dt):
    audits = []
    for event in events:
        if event.get("kind") != "dynamic_cohesive_insertion":
            continue
        step = int(event["step"])
        before, after = rows.get(step - 1), rows.get(step)
        following = rows.get(step + 1)
        audit = {
            "event": event,
            "dt_s": dt,
            "sampling_order": "top BC set; checkCohesiveStress may insert; solveStep; append_row samples reaction/work and native energies",
            "pre_row_step": step - 1 if before else None,
            "post_row_step": step if after else None,
            "pre_row_time_s": f(before, "time_s") if before else None,
            "post_row_time_s": f(after, "time_s") if after else None,
            "native_row_energy_delta_post_minus_pre_J": {},
            "native_row_reaction_work_samples": {},
            "following_step_dissipated_J": f(following, "cohesive_dissipated_J") if following else None,
            "following_step_dissipated_delta_J": None,
        }
        if before and after:
            audit["native_row_energy_delta_post_minus_pre_J"] = {
                key: f(after, key) - f(before, key) for key in ENERGY_FIELDS
            }
            audit["native_row_reaction_work_samples"] = {
                key: {"pre": f(before, key), "post": f(after, key)} for key in REACTION_FIELDS
            }
        if after and following:
            audit["following_step_dissipated_delta_J"] = f(following, "cohesive_dissipated_J") - f(after, "cohesive_dissipated_J")
        audits.append(audit)
    return audits


def first_negative_dissipated_window(rows, events, padding=2):
    negative_steps = sorted(step for step, row in rows.items() if f(row, "cohesive_dissipated_J") < 0.0)
    if not negative_steps:
        return {"first_negative_step": None, "rows_from_first_negative_minus_2_to_plus_2": []}
    first_step = negative_steps[0]
    event_steps = {int(event["step"]): event for event in events if event.get("kind") == "dynamic_cohesive_insertion"}
    window = []
    for step in range(first_step - padding, first_step + padding + 1):
        row = rows.get(step)
        if row is None:
            continue
        previous = rows.get(step - 1)
        window.append({
            "step": step,
            "time_s": f(row, "time_s"),
            "time_step_s": f(row, "time_s") - f(previous, "time_s") if previous else None,
            "dissipated_J": f(row, "cohesive_dissipated_J"),
            "delta_dissipated_from_previous_row_J": f(row, "cohesive_dissipated_J") - f(previous, "cohesive_dissipated_J") if previous else None,
            "kinetic_J": f(row, "kinetic_J"),
            "bulk_potential_J": f(row, "bulk_potential_J"),
            "reversible_J": f(row, "cohesive_reversible_J"),
            "external_work_J": f(row, "external_work_J"),
            "energy_ledger_residual_J": f(row, "energy_ledger_residual_J"),
            "insertion_event": event_steps.get(step),
        })
    return {
        "first_negative_step": first_step,
        "first_negative_time_s": f(rows[first_step], "time_s"),
        "first_negative_steps_after_first_insertion": first_step - min(event_steps) if event_steps else None,
        "rows_from_first_negative_minus_2_to_plus_2": window,
    }


def stats(values):
    arr = np.asarray(values, dtype=np.float64)
    return {"min": float(arr.min()) if arr.size else None,
            "mean": float(arr.mean()) if arr.size else None,
            "max": float(arr.max()) if arr.size else None,
            "count": int(arr.size)}


mesh = parse_msh(INPUT / "block_2x2x4mm.msh")
manifest = json.loads((DATA / "manifest.json").read_text(encoding="utf-8"))
generated = json.loads((INPUT / "generated_geometry.json").read_text(encoding="utf-8"))
input_tet_vol = signed_tet_volumes(mesh["nodes"], mesh["tets"])
input_faces = face_counts(mesh["tets"])
interface_sets = {tuple(sorted(int(v) for v in tri)) for tri in mesh["interface_triangles"]}
input_interface_counts = {str(tri): input_faces.get(tri, 0) for tri in sorted(interface_sets)}
input_face_info = face_summary(input_faces)

result = {
    "schema": "cohesive_diagnostics_20261004.v1",
    "readonly_inputs": {},
    "official_tag_source": {},
    "mesh": {},
    "native_cases": {},
    "boundary_conditions": {},
    "energy_timing_assessment": {},
    "api_semantics_and_limits": {},
    "assessment": {},
}

input_paths = {
    "run_contract.py": RUNNER,
    "block_2x2x4mm.msh": INPUT / "block_2x2x4mm.msh",
    "material.dat": INPUT / "material.dat",
    "generated_geometry.json": INPUT / "generated_geometry.json",
    "manifest.json": DATA / "manifest.json",
    "native_visualization_data.json": DATA / "native_visualization_data.json",
    "api_semantics_review.json": DATA / "api_semantics_review.json",
    "contract_review.json": DATA / "contract_review.json",
}
for name in CASES:
    case_dir = DATA / name
    for filename in ("status.json", "events.json", "trace.csv", "lineage.json", "native_topology.npz"):
        input_paths[f"{name}/{filename}"] = case_dir / filename
result["readonly_inputs"] = {name: {"path": str(path), "exists": path.exists(), "bytes": path.stat().st_size if path.exists() else None, "sha256": sha(path) if path.exists() else None} for name, path in input_paths.items()}

tag_sources = {
    "solid_mechanics_model.cc": SOURCE / "src/model/solid_mechanics/solid_mechanics_model.cc",
    "solid_mechanics_model_cohesive.cc": SOURCE / "src/model/solid_mechanics/solid_mechanics_model_cohesive/solid_mechanics_model_cohesive.cc",
    "material_cohesive.cc": SOURCE / "src/model/solid_mechanics/solid_mechanics_model_cohesive/materials/material_cohesive.cc",
    "material_cohesive_linear.cc": SOURCE / "src/model/solid_mechanics/solid_mechanics_model_cohesive/materials/constitutive_laws/material_cohesive_linear.cc",
    "source_manifest.json": SOURCE / "source_manifest.json",
}
result["official_tag_source"] = {name: {"path": str(path), "sha256": sha(path) if path.exists() else None, "exists": path.exists()} for name, path in tag_sources.items()}

input_coords_match_count = 0
input_tets_match_count = 0
interface_checks_by_case = {}
volumes_by_case = {}
pairings_by_case = {}
openings_by_case = {}
energy_by_case = {}
post_boundary_by_case = {}
exact_rest_by_case = {}

for case_name in CASES:
    case_dir = DATA / case_name
    status = json.loads((case_dir / "status.json").read_text(encoding="utf-8"))
    events = json.loads((case_dir / "events.json").read_text(encoding="utf-8"))
    lineage = json.loads((case_dir / "lineage.json").read_text(encoding="utf-8"))
    node_map = {int(k): int(v) for k, v in lineage["current_node_to_input_material_node_id"].items()}
    with np.load(case_dir / "native_topology.npz") as data:
        native_pre_nodes = np.array(data["original_rest_nodes"], copy=True)
        native_pre_tets = np.array(data["original_bulk_tetra_connectivity"], dtype=np.int64, copy=True)
        native_pre_interface = np.array(data["original_interface_triangle_connectivity"], dtype=np.int64, copy=True)
        native_post_rest = np.array(data["final_rest_nodes"], copy=True)
        native_post_tets = np.array(data["final_bulk_tetra_connectivity"], dtype=np.int64, copy=True)
        cohesive = np.array(data["final_cohesive_connectivity"], dtype=np.int64, copy=True)
        displacement = np.array(data["final_displacement"], copy=True)
        velocity = np.array(data["final_velocity"], copy=True)

    pre_coords_equal_input = native_pre_nodes.shape == mesh["nodes"].shape and np.array_equal(native_pre_nodes, mesh["nodes"])
    pre_tets_equal_input = native_pre_tets.shape == mesh["tets"].shape and np.array_equal(native_pre_tets, mesh["tets"])
    backmapped_post_tets = np.asarray([[node_map[int(n)] for n in tet] for tet in native_post_tets], dtype=np.int64)
    post_tets_backmap_equal_input = backmapped_post_tets.shape == mesh["tets"].shape and np.array_equal(backmapped_post_tets, mesh["tets"])
    exact_rest = []
    for current_id in range(native_post_rest.shape[0]):
        ancestor = node_map.get(current_id)
        exact_rest.append(ancestor is not None and np.array_equal(native_post_rest[current_id], mesh["nodes"][ancestor]))
    exact_rest_by_case[case_name] = {"current_nodes": int(native_post_rest.shape[0]), "mapped_nodes": len(node_map), "exact_coordinate_match_count": int(sum(exact_rest)), "exact_all": bool(len(exact_rest) == len(node_map) and all(exact_rest)), "missing_current_ids": sorted(set(range(native_post_rest.shape[0])) - set(node_map))}

    pre_v = signed_tet_volumes(native_pre_nodes, native_pre_tets)
    post_rest_v = signed_tet_volumes(native_post_rest, native_post_tets)
    post_deformed_v = signed_tet_volumes(native_post_rest + displacement, native_post_tets)
    pre_faces = face_counts(native_pre_tets)
    post_faces = face_counts(native_post_tets)
    pre_face_summary, post_face_summary = face_summary(pre_faces), face_summary(post_faces)

    post_interface_boundary_faces = {}
    for face, count in post_faces.items():
        if count != 1:
            continue
        ancestor_face = tuple(sorted(node_map[int(n)] for n in face))
        if ancestor_face in interface_sets:
            post_interface_boundary_faces.setdefault(ancestor_face, []).append(face)
    split_expected = case_name != "zero_load"
    expected_post_boundary_faces = 96 if split_expected else 80
    expected_post_interior_faces = 144 if split_expected else 152
    interface_checks = {
        "post_topology_state": "eight_interface_triangles_each_split_into_two_boundary_sides" if split_expected else "uninserted_interface_triangles_remain_shared_interior_faces",
        "post_face_count_expectation_derivation": "384 tetra face incidences; after split, 80 exterior + 16 interface boundary faces = 96 count-1, leaving 288 incidences / 2 = 144 count-2 faces; before insertion, 80 count-1 and 304 / 2 = 152 count-2 faces.",
        "input_interface_triangle_count": int(mesh["interface_triangles"].shape[0]),
        "native_pre_interface_triangle_count": int(native_pre_interface.shape[0]),
        "native_pre_interface_exact_input_set": {tuple(sorted(int(v) for v in tri)) for tri in native_pre_interface} == interface_sets,
        "each_input_interface_face_incidence_pre": {str(tri): pre_faces.get(tri, 0) for tri in sorted(interface_sets)},
        "post_boundary_face_count_per_input_interface_triangle": {str(tri): len(post_interface_boundary_faces.get(tri, [])) for tri in sorted(interface_sets)},
        "post_two_distinct_boundary_sides_each_interface_face": (all(len(post_interface_boundary_faces.get(tri, [])) == 2 and post_interface_boundary_faces[tri][0] != post_interface_boundary_faces[tri][1] for tri in interface_sets) if split_expected else None),
        "post_uninserted_interface_faces_remain_interior_each": (None if split_expected else all(post_faces.get(tri, 0) == 2 for tri in interface_sets)),
        "input_pre_and_native_pre_face_conforming": pre_face_summary["nonmanifold_faces_gt2"] == 0 and pre_face_summary["boundary_faces_count1"] == 80 and all(pre_faces.get(tri, 0) == 2 for tri in interface_sets),
        "native_post_face_conforming": post_face_summary["nonmanifold_faces_gt2"] == 0 and post_face_summary["zero_incidence_faces"] == 0 and post_face_summary["maximum_incidence"] <= 2 and post_face_summary["boundary_faces_count1"] == expected_post_boundary_faces and post_face_summary["interior_faces_count2"] == expected_post_interior_faces,
    }
    interface_checks_by_case[case_name] = interface_checks
    volumes_by_case[case_name] = {
        "input": volume_summary(input_tet_vol),
        "native_pre": volume_summary(pre_v),
        "native_post_rest": volume_summary(post_rest_v),
        "native_post_deformed": volume_summary(post_deformed_v),
        "native_pre_rest_coordinates_exact_input": pre_coords_equal_input,
        "native_pre_connectivity_exact_input": pre_tets_equal_input,
        "native_post_tetra_connectivity_backmaps_exact_input": post_tets_backmap_equal_input,
        "target_signed_reference_volume_m3": generated["bulk_volume_m3"],
    }

    parent_rows = []
    cell_centers = native_post_rest[native_post_tets].mean(axis=1)
    cohesive_details = []
    opening_records = []
    interface_node_pairs = {}
    for cohesive_id, row in enumerate(cohesive):
        side_nodes = [row[:3].astype(int), row[3:].astype(int)]
        ancestor_ids = [[node_map[int(n)] for n in side] for side in side_nodes]
        parents = []
        for side in side_nodes:
            face_set = set(int(n) for n in side)
            parent_hits = [i for i, tet in enumerate(native_post_tets) if face_set.issubset(set(int(n) for n in tet))]
            parents.append(parent_hits)
        side_center_z = [float(np.mean(cell_centers[hits, 2])) if hits else None for hits in parents]
        upper_side = 0 if side_center_z[0] is not None and side_center_z[0] > PLANE_Z else (1 if side_center_z[1] is not None and side_center_z[1] > PLANE_Z else None)
        lower_side = 1 - upper_side if upper_side is not None else None
        ordinal_equal = len(ancestor_ids[0]) == 3 and len(ancestor_ids[1]) == 3 and all(a == b for a, b in zip(ancestor_ids[0], ancestor_ids[1]))
        sorted_set_equal = sorted(ancestor_ids[0]) == sorted(ancestor_ids[1])
        rest_pair_equal = all(np.array_equal(native_post_rest[int(a)], native_post_rest[int(b)]) for a, b in zip(side_nodes[0], side_nodes[1]))
        cohesive_details.append({
            "cohesive_element_row_id": cohesive_id,
            "native_connectivity_six_ids": [int(v) for v in row],
            "side_0_current_ids": [int(v) for v in side_nodes[0]],
            "side_1_current_ids": [int(v) for v in side_nodes[1]],
            "side_0_input_ancestor_ids_zero_based": ancestor_ids[0],
            "side_1_input_ancestor_ids_zero_based": ancestor_ids[1],
            "side_0_input_ancestor_ids_gmsh_1based": [v + 1 for v in ancestor_ids[0]],
            "side_1_input_ancestor_ids_gmsh_1based": [v + 1 for v in ancestor_ids[1]],
            "pairwise_local_ordinal_ancestor_equal": ordinal_equal,
            "sorted_ancestor_sets_equal": sorted_set_equal,
            "paired_rest_coordinates_exact_equal": rest_pair_equal,
            "adjacent_native_post_tet_rows_by_side": parents,
            "adjacent_tet_centroid_z_m_by_side": side_center_z,
            "upper_side_index": upper_side,
            "lower_side_index": lower_side,
        })
        if upper_side is not None and lower_side is not None:
            upper_ids, lower_ids = side_nodes[upper_side], side_nodes[lower_side]
            for ordinal, (upper_node, lower_node) in enumerate(zip(upper_ids, lower_ids)):
                ancestor = node_map[int(upper_node)]
                upper_disp = displacement[int(upper_node)]
                lower_disp = displacement[int(lower_node)]
                jump = (native_post_rest[int(upper_node)] + upper_disp) - (native_post_rest[int(lower_node)] + lower_disp)
                opening_records.append({
                    "cohesive_element_row_id": cohesive_id,
                    "local_ordinal": ordinal,
                    "input_ancestor_id_gmsh_1based": ancestor + 1,
                    "rest_pair_exact": bool(np.array_equal(native_post_rest[int(upper_node)], native_post_rest[int(lower_node)])),
                    "upper_node_id": int(upper_node),
                    "lower_node_id": int(lower_node),
                    "upper_displacement_m": upper_disp.tolist(),
                    "lower_displacement_m": lower_disp.tolist(),
                    "opening_vector_upper_minus_lower_m": jump.tolist(),
                    "normal_opening_plus_z_m": float(jump[2]),
                    "tangential_opening_vector_xy_m": jump[:2].tolist(),
                    "tangential_opening_magnitude_m": float(np.linalg.norm(jump[:2])),
                    "input_x_y_m": mesh["nodes"][ancestor, :2].tolist(),
                })
                interface_node_pairs.setdefault(ancestor, (int(upper_node), int(lower_node)))

    pairings_by_case[case_name] = {
        "cohesive_element_count": int(cohesive.shape[0]),
        "all_eight_facets_present": int(cohesive.shape[0]) == 8 if case_name != "zero_load" else int(cohesive.shape[0]) == 0,
        "pairwise_local_ordinal_equal_count": sum(d["pairwise_local_ordinal_ancestor_equal"] for d in cohesive_details),
        "sorted_set_equal_count": sum(d["sorted_ancestor_sets_equal"] for d in cohesive_details),
        "rest_pair_exact_equal_count": sum(d["paired_rest_coordinates_exact_equal"] for d in cohesive_details),
        "details": cohesive_details,
    }
    if opening_records:
        normals = [r["normal_opening_plus_z_m"] for r in opening_records]
        tangential = [r["tangential_opening_magnitude_m"] for r in opening_records]
        upper_mean = np.mean([displacement[pair[0], :2] for pair in interface_node_pairs.values()], axis=0)
        lower_mean = np.mean([displacement[pair[1], :2] for pair in interface_node_pairs.values()], axis=0)
        reflected_normal_diffs = []
        by_xy = {(round(r["input_x_y_m"][0], 15), round(r["input_x_y_m"][1], 15)): r["normal_opening_plus_z_m"] for r in opening_records}
        for (x, y), value in list(by_xy.items()):
            for key in ((round(-x, 15), y), (x, round(-y, 15))):
                if key in by_xy:
                    reflected_normal_diffs.append(abs(value - by_xy[key]))
        openings_by_case[case_name] = {
            "pair_count": len(opening_records),
            "normal_opening_plus_z_m": stats(normals),
            "tangential_opening_magnitude_m": stats(tangential),
            "max_tangential_to_mean_normal_ratio": float(max(tangential) / abs(np.mean(normals))) if abs(np.mean(normals)) > 0 else None,
            "normal_opening_positive_pair_count": sum(v > 0.0 for v in normals),
            "normal_opening_negative_pair_count": sum(v < 0.0 for v in normals),
            "interface_upper_mean_xy_displacement_m": upper_mean.tolist(),
            "interface_lower_mean_xy_displacement_m": lower_mean.tolist(),
            "interface_upper_minus_lower_mean_xy_shift_m": (upper_mean - lower_mean).tolist(),
            "max_normal_opening_reflection_mismatch_x_or_y_m": max(reflected_normal_diffs) if reflected_normal_diffs else None,
            "pairs": opening_records,
        }
    else:
        openings_by_case[case_name] = {"pair_count": 0, "normal_opening_plus_z_m": None, "tangential_opening_magnitude_m": None, "pairs": []}

    rows = rows_by_step(case_dir / "trace.csv")
    dt = float(status["time_step_s"])
    audits = trace_event_audit(rows, events, dt)
    energy_by_case[case_name] = {
        "status": status["status"],
        "steps_completed": status["steps_completed"],
        "time_step_s": dt,
        "stable_time_step_at_end_s": status["stable_time_step_s"],
        "dt_over_end_stable_dt": dt / float(status["stable_time_step_s"]),
        "first_insertion_step": status["first_insertion_step"],
        "first_insertion_time_s": status["first_insertion_time_s"],
        "events": events,
        "event_sample_deltas": audits,
        "first_negative_dissipated_window": first_negative_dissipated_window(rows, events),
        "final_energy_fields_J": status["final_energies_J"],
        "dissipated_energy_target_J": status["dissipated_energy_target_J"],
        "dissipated_energy_relative_error_raw": status["dissipated_energy_relative_error"],
        "external_work_final_J": status["external_work_final_J"],
        "energy_ledger_residual_final_J": status["energy_ledger_residual_final_J"],
        "trace_row_count": len(rows),
        "trace_sha256_matches_status": sha(case_dir / "trace.csv") == status["trace_csv_sha256"],
    }
    post_boundary_by_case[case_name] = {"pre": pre_face_summary, "post": post_face_summary, "interface_two_boundary_sides_per_triangle": interface_checks["post_boundary_face_count_per_input_interface_triangle"]}

mesh_info = {
    "input_node_count": len(mesh["nodes"]),
    "input_tetrahedron_count": int(mesh["tets"].shape[0]),
    "input_physical_interface_triangle_count": int(mesh["interface_triangles"].shape[0]),
    "input_face_conformity": input_face_info,
    "input_interface_face_incidence_each": input_interface_counts,
    "generated_geometry_bulk_volume_m3": generated["bulk_volume_m3"],
    "input_tetra_signed_volume": volume_summary(input_tet_vol),
}
result["mesh"] = mesh_info
for name in CASES:
    result["native_cases"][name] = {
        "signed_volumes_m3": volumes_by_case[name],
        "post_face_conformity": post_boundary_by_case[name],
        "interface_boundary_checks": interface_checks_by_case[name],
        "exact_rest_coordinates_by_ancestry": exact_rest_by_case[name],
        "cohesive_pairing": pairings_by_case[name],
        "modeI_opening": openings_by_case[name],
        "energy_and_insertion_timing": energy_by_case[name],
    }

runner_text = RUNNER.read_text(encoding="utf-8")
bc_lines = [line.strip() for line in runner_text.splitlines() if ".applyBC(" in line]
bc_dof_literals = sorted(set(re.findall(r"aka\._([xyz])", "\n".join(bc_lines))))
result["boundary_conditions"] = {
    "runner_applyBC_calls": bc_lines,
    "explicitly_constrained_dofs": bc_dof_literals,
    "in_plane_x_y_constrained": "x" in bc_dof_literals or "y" in bc_dof_literals,
    "free_rigid_modes_from_declared_constraints": ["global x translation", "global y translation", "rotation about z"],
    "assessment": "Only z is constrained: zero_load fixes top and bottom z; mode-I fixes bottom z and prescribes top z. No x/y anchors or rotation pin are declared. In explicit dynamics these modes are not an algebraic singularity, but can permit lateral drift if asymmetric insertion/forces excite them.",
}

api_review = json.loads((DATA / "api_semantics_review.json").read_text(encoding="utf-8"))
result["energy_timing_assessment"] = {
    "runner_energy_sampling": "append_row calls energy_summary after model.solveStep; energy_summary samples model.getEnergy for kinetic, potential, dissipated, reversible, contact.",
    "runner_event_order": "At each step: set top BC at current_time; call checkCohesiveStress (may insert); then solveStep; then append_row and sample force/work/energy.",
    "external_work_trace_semantics": "runner external_work_J is a custom cumulative trapezoid of sampled constrained reaction times imposed-displacement increments; it is not native getEnergy('external work').",
    "native_getEnergy_external_work_source_review": api_review["findings"]["getEnergy_external_work"]["finding"],
    "native_getEnergy_dissipated_source_review": api_review["findings"]["dissipated_energy_sign"]["finding"],
    "native_energy_issue_root_cause": "UNKNOWN; trace sampling coalesces insertion and the ensuing solve into one row, and wheel/source correspondence is UNKNOWN.",
    "wheel_source_correspondence": api_review["source_review"]["exact_binary_source_correspondence"],
    "loaded_cases": {name: energy_by_case[name] for name in ("modeI_dt", "modeI_dt2")},
}
result["api_semantics_and_limits"] = {
    "wheel": api_review["installed_wheel"],
    "tag": api_review["source_review"],
    "reviewed_source_hashes": {name: value["sha256"] for name, value in result["official_tag_source"].items() if name.endswith(".cc")},
    "cohesive_damage_or_failure_export": "UNKNOWN; runner/API review did not export damage/traction/full-failure state.",
    "time_resolved_insertion_energy_or_traction_history": "UNKNOWN; available trace is per-step after solve, and native API does not expose event/lineage as the existing API probe found.",
}

all_mesh_volume_ok = all(
    volumes_by_case[name][which]["negative_count"] == 0
    and volumes_by_case[name][which]["zero_count"] == 0
    and math.isclose(volumes_by_case[name][which]["sum_signed_m3"], generated["bulk_volume_m3"], rel_tol=1e-13, abs_tol=1e-25)
    for name in CASES for which in ("native_pre", "native_post_rest")
)
all_pairing_ok = all(pairings_by_case[name]["pairwise_local_ordinal_equal_count"] == 8 and pairings_by_case[name]["rest_pair_exact_equal_count"] == 8 for name in ("modeI_dt", "modeI_dt2"))
all_rest_ok = all(exact_rest_by_case[name]["exact_all"] and volumes_by_case[name]["native_pre_rest_coordinates_exact_input"] and volumes_by_case[name]["native_pre_connectivity_exact_input"] and volumes_by_case[name]["native_post_tetra_connectivity_backmaps_exact_input"] for name in CASES)
all_face_ok = all(
    interface_checks_by_case[name]["input_pre_and_native_pre_face_conforming"]
    and interface_checks_by_case[name]["native_post_face_conforming"]
    and (interface_checks_by_case[name]["post_two_distinct_boundary_sides_each_interface_face"] is True if name != "zero_load" else interface_checks_by_case[name]["post_uninserted_interface_faces_remain_interior_each"] is True)
    for name in CASES
)
result["assessment"] = {
    "checks": {
        "runner_sha_matches_frozen_expected": sha(RUNNER) == "67dd1eb98acc5f5369b9cbf1a9bfbd7f6681068172e537768e788274cb10652d",
        "input_96_tet_signed_volume_all_positive_and_16e-9_m3": volume_summary(input_tet_vol)["negative_count"] == 0 and volume_summary(input_tet_vol)["zero_count"] == 0 and math.isclose(float(input_tet_vol.sum()), 1.6e-8, rel_tol=1e-13, abs_tol=1e-25),
        "native_pre_post_signed_reference_volumes_all_positive_and_16e-9_m3": all_mesh_volume_ok,
        "pre_post_interface_faces_conforming_and_split_as_eight_pairs": all_face_ok,
        "post_nodes_exact_rest_match_input_ancestors_no_nearest": all_rest_ok,
        "coh6_each_local_ordinal_pair_has_same_ancestor": all_pairing_ok,
        "zero_load_has_no_insertion_and_all_zero_energies": energy_by_case["zero_load"]["first_insertion_step"] is None and all(float(v) == 0.0 for v in energy_by_case["zero_load"]["final_energy_fields_J"].values()),
    },
    "most_supported_interpretation": "Signed volumes, exact rest ancestry, face incidences under the correct inserted/uninserted topology, and coh6 ordinal pairing are consistent, so the inspected evidence does not support malformed mesh or sorted-only pairing as the negative-dissipation cause. The first insertion row has D=0 and the next row D>0; D first turns negative two steps after the first insertion in both loaded traces (coarse step 652; dt/2 step 1299), with dt/2 step 1299 preceding the second insertion at 1300. This localizes the observed failure to subsequent native cohesive state evolution or its sampled history, but exact API versus integration cause remains UNKNOWN.",
    "single_minimal_next_diagnostic": "Instrument one replay only around dt/2 steps 1297-1301 (first 2-facet insertion, first negative D at 1299, and second insertion at 1300): sample native K/U/D/reversible/contact and cohesive traction/opening immediately before and after each checkCohesiveStress, then immediately after solveStep; record reaction, imposed BC, dt, wheel hash, and all insertion state. This spans the first-negative row and two following steps to separate insertion from integration/history evolution without changing mesh/material/loading or running a full horizon.",
}

OUT.mkdir(parents=True, exist_ok=True)
json_path = OUT / "cohesive_diagnostics.json"
json_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
print(json.dumps({
    "output": str(json_path),
    "runner_sha256": result["readonly_inputs"]["run_contract.py"]["sha256"],
    "mesh": result["mesh"],
    "checks": result["assessment"]["checks"],
    "openings": {name: {k: v for k, v in result["native_cases"][name]["modeI_opening"].items() if k != "pairs"} for name in ("modeI_dt", "modeI_dt2")},
    "energy": {name: {k: v for k, v in result["native_cases"][name]["energy_and_insertion_timing"].items() if k != "events"} for name in ("modeI_dt", "modeI_dt2")},
}, indent=2))
