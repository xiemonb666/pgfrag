"""Summarize fixed contract outputs and inspect native paired-face geometry."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
REVIEW = ROOT.parents[2] / "relicweave-synlab" / "tmp" / "review_20261004" / "implementation" / "cohesive_capability"


def unit_outward_normal(points: np.ndarray, opposite: np.ndarray) -> tuple[float, np.ndarray]:
    cross = np.cross(points[1] - points[0], points[2] - points[0])
    norm = float(np.linalg.norm(cross))
    if norm == 0.0:
        raise ValueError("zero-area face")
    normal = cross / norm
    centroid = points.mean(axis=0)
    if float(np.dot(normal, opposite - centroid)) > 0.0:
        normal = -normal
    return 0.5 * norm, normal


def inspect_case(name: str, target_area: float) -> dict:
    directory = REVIEW / name
    status = json.loads((directory / "status.json").read_text(encoding="utf-8"))
    lineage = json.loads((directory / "lineage.json").read_text(encoding="utf-8"))
    data = np.load(directory / "native_topology.npz")
    nodes = data["final_rest_nodes"]
    displacement = data["final_displacement"]
    deformed_nodes = nodes + displacement
    tets = data["final_bulk_tetra_connectivity"]
    original_tets = data["original_bulk_tetra_connectivity"]
    cohesive = data["final_cohesive_connectivity"]
    candidate_map: dict[int, int] = {}
    mapping_conflicts = []
    ordinal_pairs_checked = 0
    for element_id, (source_tet, post_tet) in enumerate(zip(original_tets, tets)):
        for local_ordinal, (source_node, current_node) in enumerate(zip(source_tet, post_tet)):
            ordinal_pairs_checked += 1
            old = candidate_map.get(int(current_node))
            if old is not None and old != int(source_node):
                mapping_conflicts.append(
                    {
                        "current_node_id": int(current_node),
                        "prior_input_material_node_id": old,
                        "conflicting_input_material_node_id": int(source_node),
                        "native_bulk_element_id": element_id,
                        "local_vertex_ordinal": local_ordinal,
                    }
                )
            else:
                candidate_map[int(current_node)] = int(source_node)
    row_backmap_matches = [
        bool(np.array_equal(np.asarray([candidate_map[int(n)] for n in post_tet]), source_tet))
        for source_tet, post_tet in zip(original_tets, tets)
    ]
    unchanged_rows = [bool(np.array_equal(a, b)) for a, b in zip(original_tets, tets)]
    original_node_identity_errors = [
        node for node in range(len(data["original_rest_nodes"])) if candidate_map.get(node) != node
    ]
    mapping_validation = {
        "method": "same native bulk element row ID + local vertex ordinal; no coordinates used",
        "native_bulk_rows_before_after": [int(len(original_tets)), int(len(tets))],
        "native_bulk_element_ids_and_local_ordinals_checked": int(ordinal_pairs_checked),
        "current_to_material_mapping_conflict_count": len(mapping_conflicts),
        "mapping_conflicts": mapping_conflicts,
        "input_node_identity_errors": original_node_identity_errors,
        "input_nodes_identity_checked": int(len(data["original_rest_nodes"])),
        "post_bulk_rows_backmap_exact_count": int(sum(row_backmap_matches)),
        "post_bulk_rows_backmap_checked_count": int(len(row_backmap_matches)),
        "bulk_connectivity_rows_bytewise_unchanged_count": int(sum(unchanged_rows)),
        "bulk_connectivity_rows_changed_count": int(len(unchanged_rows) - sum(unchanged_rows)),
        "current_nodes_mapped": int(len(set(map(int, tets.flat)))),
        "current_nodes_with_material_ancestry": int(len(candidate_map)),
        "pass": (
            len(original_tets) == len(tets)
            and not mapping_conflicts
            and not original_node_identity_errors
            and all(row_backmap_matches)
        ),
    }
    if lineage.get("status") == "EXACT_FOR_CONTRACT_MESH_FROM_NATIVE_CONNECTIVITY" and not mapping_validation["pass"]:
        raise ValueError(f"runner's EXACT lineage status failed independent mapping audit: {mapping_validation}")
    lineage["cell_local_ordinal_mapping_validation"] = mapping_validation
    sides_area = [0.0, 0.0]
    normal_dots = []
    paired_gap_max = 0.0
    facets = []
    visualization_facets = []
    for facet in lineage["facets"]:
        eid = int(facet["cohesive_element_id"])
        conn = cohesive[eid]
        side_nodes = [list(map(int, conn[:3])), list(map(int, conn[3:]))]
        parent_lists = facet["post_insert_parent_tetrahedron_ids_by_cohesive_side"]
        material_maps = [
            {int(k): int(v) for k, v in facet["current_node_to_input_material_node_id"].items()}
        ]
        material_map = material_maps[0]
        coords_by_side = []
        points_by_side = []
        areas = []
        normals = []
        for side_idx, (side, parent_ids) in enumerate(zip(side_nodes, parent_lists)):
            if len(parent_ids) != 1:
                raise ValueError(f"facet {eid}, side {side_idx} has {len(parent_ids)} parents")
            parent_tet = tets[int(parent_ids[0])]
            opposite_ids = [int(n) for n in parent_tet if int(n) not in side]
            if len(opposite_ids) != 1:
                raise ValueError(f"facet {eid}, side {side_idx} lacks one tetrahedron opposite node")
            area, normal = unit_outward_normal(nodes[np.asarray(side)], nodes[opposite_ids[0]])
            areas.append(area)
            normals.append(normal)
            sides_area[side_idx] += area
            by_material = {
                material_map[node]: nodes[node].copy()
                for node in side
            }
            coords_by_side.append(by_material)
            points_by_side.append(
                {
                    str(material_map[node]): {
                        "current_node_id": node,
                        "rest_xyz_m": [float(x) for x in nodes[node]],
                        "displacement_xyz_m": [float(x) for x in displacement[node]],
                        "deformed_xyz_m": [float(x) for x in deformed_nodes[node]],
                    }
                    for node in side
                }
            )
        shared_material = set(coords_by_side[0]) & set(coords_by_side[1])
        if len(shared_material) != 3:
            raise ValueError(f"facet {eid} material ancestry did not pair three nodes")
        gap = max(float(np.linalg.norm(coords_by_side[0][mid] - coords_by_side[1][mid])) for mid in shared_material)
        paired_gap_max = max(paired_gap_max, gap)
        dot = float(np.dot(normals[0], normals[1]))
        normal_dots.append(dot)
        facets.append(
            {
                "cohesive_element_id": eid,
                "parent_bulk_element_ids_by_side": [list(map(int, p)) for p in parent_lists],
                "side_area_m2": areas,
                "outward_unit_normals": [n.tolist() for n in normals],
                "outward_normal_dot": dot,
                "paired_rest_gap_m": gap,
                "shared_input_material_node_ids": sorted(shared_material),
            }
        )
        visualization_facets.append(
            {
                "cohesive_element_id": eid,
                "side_a": points_by_side[0],
                "side_b": points_by_side[1],
            }
        )
    final_energy = status["final_energies_J"]
    geometry_validation = {
        "method": "native cohesive side connectivity + native adjacent tetrahedron connectivity + rest coordinates; no nearest-coordinate matching",
        "paired_rest_gap_max_m": paired_gap_max,
        "cohesive_rest_side_area_sums_m2": sides_area,
        "target_rest_area_m2_per_side": target_area,
        "area_relative_errors": [abs(a - target_area) / target_area for a in sides_area] if target_area else None,
        "outward_normal_dot_min_max": [min(normal_dots), max(normal_dots)] if normal_dots else None,
        "facets": facets,
    }
    lineage["cohesive_rest_geometry_validation"] = geometry_validation
    (directory / "lineage.json").write_text(json.dumps(lineage, indent=2) + "\n", encoding="utf-8")
    return {
        "case": name,
        "run_status": status["status"],
        "steps_completed": status["steps_completed"],
        "endtime_reached": status["endtime_reached"],
        "first_insertion_step": status["first_insertion_step"],
        "first_insertion_time_s": status["first_insertion_time_s"],
        "native_insertions": status["dynamic_insertions_native_connectivity"],
        "api_insertions": status["dynamic_insertions_api_returned"],
        "cohesive_element_count": len(cohesive),
        "lineage_status": lineage["status"],
        "bulk_elements_before_after": lineage["bulk_elements_before_after"],
        "current_nodes_with_material_ancestry": lineage["current_nodes_with_input_material_ancestry"],
        "cell_local_ordinal_mapping_validation": mapping_validation,
        "paired_rest_gap_max_m": paired_gap_max,
        "cohesive_rest_side_area_sums_m2": sides_area,
        "target_rest_area_m2_per_side": target_area,
        "area_relative_errors": [abs(a - target_area) / target_area for a in sides_area] if target_area else None,
        "outward_normal_dot_min_max": [min(normal_dots), max(normal_dots)] if normal_dots else None,
        "facets": facets,
        "mass_component_sums_before_kg": status["mass_component_sums_before_kg"],
        "mass_component_sums_after_kg": status["mass_component_sums_after_kg"],
        "mass_conserved_by_component": status["mass_conserved_by_component"],
        "energies_J": final_energy,
        "external_work_J": status["external_work_final_J"],
        "ledger_residual_J": status["energy_ledger_residual_final_J"],
        "ledger_relative_residual": status["energy_ledger_relative_residual"],
        "cohesive_dissipation_signed_J": status["dissipated_energy_final_J"],
        "cohesive_dissipation_target_J": status["dissipated_energy_target_J"],
        "dissipation_target_gate": "NOT_REACHED" if status["complete_debond_step_by_dissipation_target"] is None else "REACHED",
        "fragment_summary": status["fragment_summary"],
        "maximum_nodal_displacement_norm_m": float(np.linalg.norm(displacement, axis=1).max()) if len(displacement) else 0.0,
        "visualization_facets": visualization_facets,
    }


def read_trace(name: str) -> list[dict]:
    with (REVIEW / name / "trace.csv").open(newline="", encoding="utf-8") as stream:
        return [{key: float(value) if key != "step" else int(value) for key, value in row.items()} for row in csv.DictReader(stream)]


manifest = json.loads((REVIEW / "manifest.json").read_text(encoding="utf-8"))
zero = inspect_case("zero_load", 0.0)
coarse = inspect_case("modeI_dt", manifest["geometry"]["interface_area_m2"])
fine = inspect_case("modeI_dt2", manifest["geometry"]["interface_area_m2"])
coarse_d = coarse["cohesive_dissipation_signed_J"]
fine_d = fine["cohesive_dissipation_signed_J"]
result = {
    "manifest_sha256": __import__("hashlib").sha256((REVIEW / "manifest.json").read_bytes()).hexdigest(),
    "frozen_endtime_s": manifest["loading"]["endtime_s"],
    "expected_full_debond_GcA_J": manifest["expected_interface_dissipation_J"],
    "zero_load": zero,
    "modeI_dt": coarse,
    "modeI_dt2": fine,
    "dt_refinement": {
        "same_endtime": coarse["endtime_reached"] and fine["endtime_reached"],
        "coarse_signed_dissipation_J": coarse_d,
        "fine_signed_dissipation_J": fine_d,
        "absolute_difference_J": abs(fine_d - coarse_d),
        "absolute_difference_over_coarse_abs": abs(fine_d - coarse_d) / abs(coarse_d) if coarse_d else None,
        "gate": "NOT_REACHED because neither case reached the frozen positive Gc*A dissipation target",
    },
    "nonclaims": manifest["nonclaims"],
}
output = REVIEW / "contract_review.json"
output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
(REVIEW / "native_visualization_data.json").write_text(
    json.dumps(
        {
            "note": "Actual rest and deformed node coordinates from native_topology.npz; no geometry reconstruction by nearest matching.",
            "modeI_dt": {
                "maximum_nodal_displacement_norm_m": coarse["maximum_nodal_displacement_norm_m"],
                "cohesive_facets": coarse["visualization_facets"],
            },
            "modeI_dt2": {
                "maximum_nodal_displacement_norm_m": fine["maximum_nodal_displacement_norm_m"],
                "cohesive_facets": fine["visualization_facets"],
            },
        },
        indent=2,
    )
    + "\n",
    encoding="utf-8",
)
print(json.dumps(result, indent=2))
