"""Pure-Python mesh, lineage, and energy contracts for the free cohesive patch."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np


TET_FACES = ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_gmsh22(path: Path) -> dict[str, Any]:
    """Read the MSH 2.2 subset used by this frozen tetrahedral contract mesh."""
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    sections: dict[str, list[str]] = {}
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("$") and not line.startswith("$End"):
            name = line[1:]
            end_marker = "$End" + name
            i += 1
            payload: list[str] = []
            while i < len(lines) and lines[i].strip() != end_marker:
                payload.append(lines[i].strip())
                i += 1
            if i == len(lines):
                raise ValueError(f"unterminated Gmsh section {name}")
            sections[name] = payload
        i += 1

    fmt = sections.get("MeshFormat", [])
    if not fmt or fmt[0].split()[0] != "2.2" or fmt[0].split()[1] != "0":
        raise ValueError("expected ASCII Gmsh MSH 2.2")

    names: dict[tuple[int, int], str] = {}
    physical = sections.get("PhysicalNames", [])
    if physical:
        for row in physical[1:]:
            fields = row.split(maxsplit=2)
            if len(fields) != 3:
                raise ValueError(f"malformed PhysicalNames row: {row!r}")
            names[(int(fields[0]), int(fields[1]))] = json.loads(fields[2])

    node_rows = sections.get("Nodes", [])
    if not node_rows:
        raise ValueError("missing Nodes section")
    node_count = int(node_rows[0])
    node_tags: list[int] = []
    node_coordinates: list[list[float]] = []
    for row in node_rows[1:]:
        values = row.split()
        if len(values) != 4:
            raise ValueError(f"malformed node row: {row!r}")
        node_tags.append(int(values[0]))
        node_coordinates.append([float(v) for v in values[1:]])
    if len(node_tags) != node_count or len(set(node_tags)) != node_count:
        raise ValueError("node count or tag uniqueness mismatch")
    node_row_by_tag = {tag: row for row, tag in enumerate(node_tags)}

    element_rows = sections.get("Elements", [])
    if not element_rows:
        raise ValueError("missing Elements section")
    element_count = int(element_rows[0])
    elements: list[dict[str, Any]] = []
    for row in element_rows[1:]:
        values = [int(v) for v in row.split()]
        if len(values) < 3:
            raise ValueError(f"malformed element row: {row!r}")
        element_id, element_type, tag_count = values[:3]
        tags = values[3 : 3 + tag_count]
        node_ids = values[3 + tag_count :]
        if len(tags) != tag_count:
            raise ValueError(f"malformed element tags: {row!r}")
        try:
            nodes = [node_row_by_tag[node_tag] for node_tag in node_ids]
        except KeyError as error:
            raise ValueError(f"element references absent node tag {error.args[0]}") from error
        elements.append({
            "element_id": element_id,
            "element_type": element_type,
            "physical_tag": tags[0] if tags else None,
            "physical_name": names.get((2 if element_type == 2 else 3, tags[0])) if tags else None,
            "node_tags": node_ids,
            "nodes": nodes,
        })
    if len(elements) != element_count:
        raise ValueError("element count mismatch")

    tets = [e for e in elements if e["element_type"] == 4]
    triangles = [e for e in elements if e["element_type"] == 2]
    if not tets or not triangles:
        raise ValueError("contract mesh requires tetrahedra and triangle groups")
    return {
        "path": str(Path(path)),
        "physical_names": names,
        "node_tags": node_tags,
        "node_row_by_tag": node_row_by_tag,
        "nodes": np.asarray(node_coordinates, dtype=np.float64),
        "elements": elements,
        "tetrahedra": tets,
        "triangles": triangles,
    }


def surface_nodes(mesh_data: dict[str, Any], physical_name: str) -> set[int]:
    rows = [e for e in mesh_data["triangles"] if e["physical_name"] == physical_name]
    if not rows:
        raise ValueError(f"physical surface group {physical_name!r} is empty")
    return {int(n) for row in rows for n in row["nodes"]}


def triangle_area(nodes: np.ndarray, face: list[int] | tuple[int, ...]) -> float:
    p = nodes[np.asarray(face, dtype=np.int64)]
    return 0.5 * float(np.linalg.norm(np.cross(p[1] - p[0], p[2] - p[0])))


def input_interface_contract(mesh_data: dict[str, Any], interface_name: str = "interface") -> dict[str, Any]:
    """Build exact face-to-parent incidence from Gmsh tags, not coordinates."""
    incidence: dict[tuple[int, int, int], list[dict[str, int]]] = {}
    for tet_row, element in enumerate(mesh_data["tetrahedra"]):
        if len(element["nodes"]) != 4:
            raise ValueError("only linear tetrahedra are supported")
        for local_face, ordinals in enumerate(TET_FACES):
            face_nodes = tuple(int(element["nodes"][ordinal]) for ordinal in ordinals)
            incidence.setdefault(tuple(sorted(face_nodes)), []).append({
                "input_tetrahedron_row": tet_row,
                "input_tetrahedron_id": int(element["element_id"]),
                "local_face_ordinal": local_face,
            })

    interface_rows = [e for e in mesh_data["triangles"] if e["physical_name"] == interface_name]
    if not interface_rows:
        raise ValueError(f"physical surface group {interface_name!r} is empty")
    faces: dict[tuple[int, int, int], dict[str, Any]] = {}
    for element in interface_rows:
        nodes = [int(n) for n in element["nodes"]]
        key = tuple(sorted(nodes))
        if key in faces:
            raise ValueError(f"duplicate declared interface face {element['element_id']}")
        parents = incidence.get(key, [])
        if len(parents) != 2:
            raise ValueError(f"interface face {element['element_id']} has {len(parents)} input parents")
        faces[key] = {
            "input_face_id": int(element["element_id"]),
            "input_face_node_ids": nodes,
            "input_face_node_key": list(key),
            "parent_input_tetrahedron_ids": [p["input_tetrahedron_id"] for p in parents],
            "parent_input_tetrahedron_rows": [p["input_tetrahedron_row"] for p in parents],
            "parent_local_face_ordinals": [p["local_face_ordinal"] for p in parents],
            "positive_rest_area_m2": triangle_area(mesh_data["nodes"], nodes),
        }
    return {
        "faces": faces,
        "input_face_count": len(faces),
        "rest_area_sum_m2": float(sum(f["positive_rest_area_m2"] for f in faces.values())),
    }


def build_native_lineage(
    mesh_data: dict[str, Any],
    native_rest_nodes: np.ndarray,
    native_bulk_tets: np.ndarray,
    native_cohesive: np.ndarray,
    interface_name: str = "interface",
    expected_area_m2: float = 4.0e-6,
    area_tolerance_m2: float = 1.0e-12,
) -> dict[str, Any]:
    """Map each native cohesive side through source tetrahedron row/local ordinal."""
    src_nodes = mesh_data["nodes"]
    src_tets = np.asarray([row["nodes"] for row in mesh_data["tetrahedra"]], dtype=np.int64)
    src_tet_ids = [int(row["element_id"]) for row in mesh_data["tetrahedra"]]
    native_rest_nodes = np.asarray(native_rest_nodes, dtype=np.float64)
    native_bulk_tets = np.asarray(native_bulk_tets, dtype=np.int64)
    native_cohesive = np.asarray(native_cohesive, dtype=np.int64)
    if native_rest_nodes.ndim != 2 or native_rest_nodes.shape[1] != 3:
        raise ValueError(f"native node array shape invalid: {native_rest_nodes.shape}")
    if native_bulk_tets.shape != src_tets.shape:
        raise ValueError(f"native/input bulk shape mismatch: {native_bulk_tets.shape} != {src_tets.shape}")
    if native_cohesive.ndim != 2 or native_cohesive.shape[1] != 6:
        raise ValueError(f"native cohesive connectivity shape invalid: {native_cohesive.shape}")

    native_to_input: dict[int, int] = {}
    mapping_rows: list[dict[str, int]] = []
    for tet_row, (input_row, native_row) in enumerate(zip(src_tets, native_bulk_tets, strict=True)):
        for local_ordinal, (source_node, native_node) in enumerate(zip(input_row, native_row, strict=True)):
            source_node, native_node = int(source_node), int(native_node)
            previous = native_to_input.setdefault(native_node, source_node)
            if previous != source_node:
                raise ValueError(f"native node {native_node} has conflicting source ancestry")
            mapping_rows.append({
                "native_bulk_element_row": tet_row,
                "input_tetrahedron_id": src_tet_ids[tet_row],
                "local_vertex_ordinal": local_ordinal,
                "input_node_id_zero_based": source_node,
                "native_node_id": native_node,
            })
            if not np.array_equal(native_rest_nodes[native_node], src_nodes[source_node]):
                raise ValueError("source-cell/local-ordinal rest coordinate consistency failed")
    if set(native_to_input) != set(range(len(native_rest_nodes))):
        raise ValueError("bulk ancestry does not cover every native node")

    source_faces = input_interface_contract(mesh_data, interface_name)
    interface_faces = source_faces["faces"]
    native_face_incidence: dict[tuple[int, int, int], list[tuple[int, int]]] = {}
    for tet_row, tet in enumerate(native_bulk_tets):
        for local_face, ordinals in enumerate(TET_FACES):
            key = tuple(sorted(int(tet[o]) for o in ordinals))
            native_face_incidence.setdefault(key, []).append((tet_row, local_face))

    seen_faces: set[tuple[int, int, int]] = set()
    records: list[dict[str, Any]] = []
    lower_z = float(np.min(src_nodes[:, 2]))
    plane_z = float(np.median(src_nodes[:, 2]))
    del lower_z
    for cohesive_row, coh in enumerate(native_cohesive):
        sides = (list(map(int, coh[:3])), list(map(int, coh[3:])))
        source_sides = [tuple(sorted(native_to_input[node] for node in side)) for side in sides]
        if source_sides[0] != source_sides[1]:
            raise ValueError(f"cohesive row {cohesive_row} sides do not pair by source node IDs")
        source_key = source_sides[0]
        if source_key not in interface_faces:
            raise ValueError(f"cohesive row {cohesive_row} is not a declared input interface face")
        if source_key in seen_faces:
            raise ValueError("an input interface face has multiple native cohesive elements")
        seen_faces.add(source_key)
        source_face = interface_faces[source_key]
        side_parents: list[dict[str, Any]] = []
        areas: list[float] = []
        for side_nodes in sides:
            native_key = tuple(sorted(side_nodes))
            adjacent = native_face_incidence.get(native_key, [])
            if len(adjacent) != 1:
                raise ValueError(f"cohesive side {cohesive_row} has {len(adjacent)} post-insertion bulk parents")
            parent_row, parent_face_ordinal = adjacent[0]
            parent_input_id = src_tet_ids[parent_row]
            parent_centroid_z = float(src_nodes[src_tets[parent_row], 2].mean())
            side_parents.append({
                "native_bulk_element_row": int(parent_row),
                "input_tetrahedron_id": int(parent_input_id),
                "local_face_ordinal": int(parent_face_ordinal),
                "side": "lower" if parent_centroid_z < plane_z else "upper",
            })
            areas.append(triangle_area(native_rest_nodes, side_nodes))
        if side_parents[0]["side"] == side_parents[1]["side"]:
            raise ValueError("cohesive face parents must lie on opposite sides of the declared plane")
        if not np.isclose(areas[0], areas[1], rtol=0.0, atol=area_tolerance_m2):
            raise ValueError("paired cohesive rest-side areas differ")
        if areas[0] <= 0.0 or areas[1] <= 0.0:
            raise ValueError("cohesive rest-side area must be positive")
        by_side = {parent["side"]: parent["input_tetrahedron_id"] for parent in side_parents}
        records.append({
            "input_face_id": source_face["input_face_id"],
            "input_face_node_ids": source_face["input_face_node_ids"],
            "parent_input_tetrahedron_ids": source_face["parent_input_tetrahedron_ids"],
            "native_cohesive_element_id": int(cohesive_row),
            "native_cohesive_element_id_semantics": "local row index returned by getConnectivity; no persistent global cohesive ID is exposed",
            "side_a_native_node_ids": sides[0],
            "side_b_native_node_ids": sides[1],
            "side_parent_tetrahedron_ids": [p["input_tetrahedron_id"] for p in side_parents],
            "side_parent_by_geometric_side": by_side,
            "post_insertion_bulk_parent_rows": [p["native_bulk_element_row"] for p in side_parents],
            "positive_rest_area_m2": areas[0],
            "side_a_positive_rest_area_m2": areas[0],
            "side_b_positive_rest_area_m2": areas[1],
        })
    if seen_faces != set(interface_faces):
        raise ValueError(f"native cohesive/source face bijection failed: {len(seen_faces)} of {len(interface_faces)}")
    total_area = float(sum(row["positive_rest_area_m2"] for row in records))
    if not np.isclose(total_area, expected_area_m2, rtol=0.0, atol=area_tolerance_m2):
        raise ValueError(f"total positive rest interface area {total_area} != {expected_area_m2}")
    return {
        "method": "input tetrahedron element ID + local vertex ordinal + native cohesive connectivity; no coordinate node matching",
        "native_node_to_input_node_zero_based": {str(k): int(v) for k, v in sorted(native_to_input.items())},
        "native_bulk_local_ordinal_lineage_rows": mapping_rows,
        "parent_facets": records,
        "input_face_count": len(interface_faces),
        "native_cohesive_element_count": len(native_cohesive),
        "input_rest_area_sum_m2": source_faces["rest_area_sum_m2"],
        "positive_rest_area_sum_m2": total_area,
    }


def quadrature_signed_work(
    opening_previous: np.ndarray,
    traction_previous: np.ndarray,
    opening_current: np.ndarray,
    traction_current: np.ndarray,
    positive_area_per_qp_m2: np.ndarray,
) -> np.ndarray:
    """Return raw signed source-formula work for each quadrature point."""
    op0 = np.asarray(opening_previous, dtype=np.float64)
    tr0 = np.asarray(traction_previous, dtype=np.float64)
    op1 = np.asarray(opening_current, dtype=np.float64)
    tr1 = np.asarray(traction_current, dtype=np.float64)
    area = np.asarray(positive_area_per_qp_m2, dtype=np.float64)
    if not (op0.shape == tr0.shape == op1.shape == tr1.shape) or op0.ndim != 2 or op0.shape[1] != 3:
        raise ValueError("opening/traction frames must share shape (n_qp, 3)")
    if area.shape != (len(op0),) or np.any(area <= 0.0):
        raise ValueError("quadrature areas must be positive and match the number of points")
    return area * 0.5 * np.einsum("qi,qi->q", op1 - op0, tr1 + tr0)


def strict_json_dump(path: Path, value: Any) -> None:
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
