"""Exact mesh and cohesive-interface contracts for linear tetrahedral meshes.

The functions here deliberately use element/node lineage and exact coordinates.
They do not infer ancestry from proximity and do not claim that a cohesive
element has fractured.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np


TET_FACES = ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))

# Gmsh MSH 2.2 element type -> (dimension, number of nodes). The reader only
# consumes linear triangles and tetrahedra, but knows enough other types to
# reject unsupported surface cells explicitly.
_GMSH_ELEMENT_TYPES: dict[int, tuple[int, int]] = {
    1: (1, 2),   # 2-node line
    2: (2, 3),   # 3-node triangle
    3: (2, 4),   # 4-node quadrilateral
    4: (3, 4),   # 4-node tetrahedron
    5: (3, 8),   # 8-node hexahedron
    6: (3, 6),   # 6-node prism
    7: (3, 5),   # 5-node pyramid
    8: (1, 3),   # 3-node second-order line
    9: (2, 6),   # 6-node second-order triangle
    10: (2, 9),  # 9-node second-order quadrilateral
    11: (3, 10), # 10-node second-order tetrahedron
    12: (3, 27), # 27-node second-order hexahedron
    13: (3, 18), # 18-node second-order prism
    14: (3, 14), # 14-node second-order pyramid
    15: (0, 1),  # point
    16: (2, 8),  # 8-node second-order quadrilateral
    17: (3, 20), # 20-node second-order hexahedron
    18: (3, 15), # 15-node second-order prism
    19: (3, 13), # 13-node second-order pyramid
    20: (2, 9),  # 9-node third-order incomplete triangle
    21: (2, 10), # 10-node third-order triangle
    22: (2, 12), # 12-node fourth-order incomplete triangle
    23: (2, 15), # 15-node fourth-order triangle
    24: (2, 15), # 15-node fifth-order incomplete triangle
    25: (2, 21), # 21-node fifth-order triangle
    26: (1, 4),
    27: (1, 5),
    28: (1, 6),
    29: (3, 20), # 20-node third-order tetrahedron
    30: (3, 35),
    31: (3, 56),
    32: (3, 22),
    33: (3, 28),
}


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of a file."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _tetrahedron_face_incidence(tetrahedra: list[dict[str, Any]]) -> dict[tuple[int, int, int], list[dict[str, int]]]:
    incidence: dict[tuple[int, int, int], list[dict[str, int]]] = {}
    for tet_row, element in enumerate(tetrahedra):
        if len(element["nodes"]) != 4:
            raise ValueError("only linear tetrahedra are supported")
        for local_face, ordinals in enumerate(TET_FACES):
            face_nodes = tuple(int(element["nodes"][ordinal]) for ordinal in ordinals)
            if len(set(face_nodes)) != 3:
                raise ValueError(f"tetrahedron {element['element_id']} has a repeated face node")
            incidence.setdefault(tuple(sorted(face_nodes)), []).append({
                "input_tetrahedron_row": tet_row,
                "input_tetrahedron_id": int(element["element_id"]),
                "local_face_ordinal": local_face,
            })
    nonmanifold = [key for key, parents in incidence.items() if len(parents) > 2]
    if nonmanifold:
        raise ValueError(f"nonmanifold source tetrahedron face has {len(incidence[nonmanifold[0]])} parents")
    return incidence


def read_gmsh22(path: Path) -> dict[str, Any]:
    """Read ASCII Gmsh MSH 2.2 linear tetrahedra and physical triangles.

    Every explicit linear triangle is retained, including physical exterior
    faces, so callers can identify source exterior face keys. Known 2D
    element types other than linear triangles are rejected explicitly.
    """
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
    if not fmt:
        raise ValueError("missing MeshFormat section")
    format_fields = fmt[0].split()
    if len(format_fields) < 3 or format_fields[0] != "2.2" or format_fields[1] != "0":
        raise ValueError("expected ASCII Gmsh MSH 2.2")

    names: dict[tuple[int, int], str] = {}
    physical = sections.get("PhysicalNames", [])
    if physical:
        try:
            count = int(physical[0])
        except ValueError as error:
            raise ValueError("malformed PhysicalNames count") from error
        if len(physical) != count + 1:
            raise ValueError("PhysicalNames count mismatch")
        for row in physical[1:]:
            fields = row.split(maxsplit=2)
            if len(fields) != 3:
                raise ValueError(f"malformed PhysicalNames row: {row!r}")
            names[(int(fields[0]), int(fields[1]))] = json.loads(fields[2])

    node_rows = sections.get("Nodes", [])
    if not node_rows:
        raise ValueError("missing Nodes section")
    try:
        node_count = int(node_rows[0])
    except ValueError as error:
        raise ValueError("malformed Nodes count") from error
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
    nodes = np.asarray(node_coordinates, dtype=np.float64)
    if not np.all(np.isfinite(nodes)):
        raise ValueError("node coordinates must be finite")
    node_row_by_tag = {tag: row for row, tag in enumerate(node_tags)}

    element_rows = sections.get("Elements", [])
    if not element_rows:
        raise ValueError("missing Elements section")
    try:
        element_count = int(element_rows[0])
    except ValueError as error:
        raise ValueError("malformed Elements count") from error
    elements: list[dict[str, Any]] = []
    for row in element_rows[1:]:
        try:
            values = [int(v) for v in row.split()]
        except ValueError as error:
            raise ValueError(f"malformed element row: {row!r}") from error
        if len(values) < 3:
            raise ValueError(f"malformed element row: {row!r}")
        element_id, element_type, tag_count = values[:3]
        if tag_count < 0 or len(values) < 3 + tag_count:
            raise ValueError(f"malformed element tags: {row!r}")
        tags = values[3 : 3 + tag_count]
        node_ids = values[3 + tag_count :]
        element_shape = _GMSH_ELEMENT_TYPES.get(element_type)
        if element_shape is None:
            raise ValueError(f"unsupported Gmsh element type {element_type} (element {element_id})")
        dimension, expected_node_count = element_shape
        if dimension == 2 and element_type != 2:
            raise ValueError(
                f"unsupported surface element type {element_type} on element {element_id}; "
                "only linear 3-node triangles (Gmsh type 2) are supported"
            )
        if dimension == 3 and element_type != 4:
            raise ValueError(
                f"unsupported volume element type {element_type} on element {element_id}; "
                "only linear 4-node tetrahedra (Gmsh type 4) are supported"
            )
        if len(node_ids) != expected_node_count:
            raise ValueError(
                f"element {element_id} type {element_type} requires {expected_node_count} nodes, got {len(node_ids)}"
            )
        try:
            nodes_in_element = [node_row_by_tag[node_tag] for node_tag in node_ids]
        except KeyError as error:
            raise ValueError(f"element references absent node tag {error.args[0]}") from error
        elements.append({
            "element_id": element_id,
            "element_type": element_type,
            "physical_tag": tags[0] if tags else None,
            "physical_name": names.get((dimension, tags[0])) if tags else None,
            "node_tags": node_ids,
            "nodes": nodes_in_element,
        })
    if len(elements) != element_count:
        raise ValueError("element count mismatch")

    tets = [element for element in elements if element["element_type"] == 4]
    triangles = [element for element in elements if element["element_type"] == 2]
    if not tets or not triangles:
        raise ValueError("mesh requires linear tetrahedra and explicit triangle groups")
    face_incidence = _tetrahedron_face_incidence(tets)
    for element in triangles:
        key = tuple(sorted(map(int, element["nodes"])))
        parent_count = len(face_incidence.get(key, []))
        if parent_count not in (1, 2):
            raise ValueError(
                f"physical triangle {element['element_id']} has {parent_count} source tetrahedron parents"
            )
    # Keep all explicit physical triangles (both exterior and interface). The
    # incidence table lets callers select only exterior source face keys.
    exterior_triangles = [
        element for element in triangles
        if len(face_incidence.get(tuple(sorted(map(int, element["nodes"]))), [])) == 1
    ]
    return {
        "path": str(Path(path)),
        "physical_names": names,
        "node_tags": node_tags,
        "node_row_by_tag": node_row_by_tag,
        "nodes": nodes,
        "elements": elements,
        "tetrahedra": tets,
        "triangles": triangles,
        "exterior_triangles": exterior_triangles,
        "tetrahedron_face_incidence": face_incidence,
    }


def surface_nodes(mesh_data: dict[str, Any], physical_name: str) -> set[int]:
    rows = [element for element in mesh_data["triangles"] if element["physical_name"] == physical_name]
    if not rows:
        raise ValueError(f"physical surface group {physical_name!r} is empty")
    return {int(node) for row in rows for node in row["nodes"]}


def triangle_area(nodes: np.ndarray, face: list[int] | tuple[int, ...]) -> float:
    coordinates = np.asarray(nodes, dtype=np.float64)
    face_ids = np.asarray(face, dtype=np.int64)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError("triangle coordinate array must have shape (n, 3)")
    if face_ids.shape != (3,) or len(set(map(int, face_ids))) != 3:
        raise ValueError("a triangle face must contain three distinct node indices")
    if np.any(face_ids < 0) or np.any(face_ids >= len(coordinates)):
        raise ValueError("triangle face references a node outside the coordinate array")
    points = coordinates[face_ids]
    if not np.all(np.isfinite(points)):
        raise ValueError("triangle coordinates must be finite")
    return 0.5 * float(np.linalg.norm(np.cross(points[1] - points[0], points[2] - points[0])))


def input_interface_contract(mesh_data: dict[str, Any], interface_name: str = "interface") -> dict[str, Any]:
    """Build exact face-to-parent incidence from the source tetrahedra."""
    incidence = mesh_data.get("tetrahedron_face_incidence")
    if incidence is None:
        incidence = _tetrahedron_face_incidence(mesh_data["tetrahedra"])
    interface_rows = [element for element in mesh_data["triangles"] if element["physical_name"] == interface_name]
    if not interface_rows:
        raise ValueError(f"physical surface group {interface_name!r} is empty")
    faces: dict[tuple[int, int, int], dict[str, Any]] = {}
    for element in interface_rows:
        nodes = [int(node) for node in element["nodes"]]
        if len(nodes) != 3 or len(set(nodes)) != 3:
            raise ValueError(f"interface face {element['element_id']} is not a three-node triangle")
        key = tuple(sorted(nodes))
        if key in faces:
            raise ValueError(f"duplicate declared interface face {element['element_id']}")
        parents = incidence.get(key, [])
        if len(parents) != 2:
            raise ValueError(f"interface face {element['element_id']} has {len(parents)} input parents")
        area = triangle_area(mesh_data["nodes"], nodes)
        if area <= 0.0:
            raise ValueError(f"interface face {element['element_id']} has nonpositive rest area")
        faces[key] = {
            "input_face_id": int(element["element_id"]),
            "input_face_node_ids": nodes,
            "input_face_node_key": list(key),
            "parent_input_tetrahedron_ids": [int(parent["input_tetrahedron_id"]) for parent in parents],
            "parent_input_tetrahedron_rows": [int(parent["input_tetrahedron_row"]) for parent in parents],
            "parent_local_face_ordinals": [int(parent["local_face_ordinal"]) for parent in parents],
            "positive_rest_area_m2": area,
        }
    return {
        "faces": faces,
        "input_face_count": len(faces),
        "rest_area_sum_m2": float(sum(face["positive_rest_area_m2"] for face in faces.values())),
    }


def build_native_lineage(
    mesh_data: dict[str, Any],
    native_rest_nodes: np.ndarray,
    native_bulk_tets: np.ndarray,
    native_cohesive: np.ndarray,
    interface_name: str = "interface",
    *,
    expected_area_m2: float,
    native_initial_bulk_tets: np.ndarray,
    plane_axis: int = 0,
    plane_coordinate_m: float = 0.0,
    area_tolerance_m2: float = 1.0e-12,
) -> dict[str, Any]:
    """Trace native cohesive faces to input source faces and parent tetrahedra.

    ``native_initial_bulk_tets`` must be the exact zero-based connectivity
    snapshot from ``getConnectivity`` before cohesive insertion. Its complete
    array must exactly equal the MSH tetrahedron rows and node ordinals. The
    post-insertion rows retain that row/local-ordinal lineage while duplicated
    node IDs are checked against exact source rest coordinates.
    """
    src_nodes = np.asarray(mesh_data["nodes"], dtype=np.float64)
    src_tets = np.asarray([row["nodes"] for row in mesh_data["tetrahedra"]], dtype=np.int64)
    src_tet_ids = [int(row["element_id"]) for row in mesh_data["tetrahedra"]]
    native_rest_nodes = np.asarray(native_rest_nodes, dtype=np.float64)
    native_bulk_tets = np.asarray(native_bulk_tets, dtype=np.int64)
    native_initial_bulk_tets = np.asarray(native_initial_bulk_tets, dtype=np.int64)
    native_cohesive = np.asarray(native_cohesive, dtype=np.int64)
    if isinstance(plane_axis, bool) or not isinstance(plane_axis, (int, np.integer)) or not 0 <= int(plane_axis) < 3:
        raise ValueError("plane_axis must be one of 0, 1, or 2")
    plane_axis = int(plane_axis)
    plane_coordinate_m = float(plane_coordinate_m)
    expected_area_m2 = float(expected_area_m2)
    area_tolerance_m2 = float(area_tolerance_m2)
    if not np.isfinite(plane_coordinate_m):
        raise ValueError("plane_coordinate_m must be finite")
    if not np.isfinite(expected_area_m2) or expected_area_m2 <= 0.0:
        raise ValueError("expected_area_m2 must be finite and positive")
    if not np.isfinite(area_tolerance_m2) or area_tolerance_m2 < 0.0:
        raise ValueError("area_tolerance_m2 must be finite and nonnegative")
    if src_nodes.ndim != 2 or src_nodes.shape[1] != 3 or not np.all(np.isfinite(src_nodes)):
        raise ValueError("source node coordinates must be a finite array with shape (n, 3)")
    if native_rest_nodes.ndim != 2 or native_rest_nodes.shape[1] != 3 or not np.all(np.isfinite(native_rest_nodes)):
        raise ValueError(f"native node array shape invalid: {native_rest_nodes.shape}")
    if src_tets.ndim != 2 or src_tets.shape[1] != 4 or np.any(src_tets < 0) or np.any(src_tets >= len(src_nodes)):
        raise ValueError("source tetrahedron connectivity must have shape (n, 4) and valid node rows")
    if native_initial_bulk_tets.shape != src_tets.shape:
        raise ValueError(
            "native initial/input bulk shape mismatch: "
            f"{native_initial_bulk_tets.shape} != {src_tets.shape}"
        )
    if not np.array_equal(native_initial_bulk_tets, src_tets):
        raise ValueError("source MSH tetrahedron order/local-vertex connectivity differs from native initial getConnectivity")
    if native_bulk_tets.shape != src_tets.shape:
        raise ValueError(f"native/input post-insertion bulk shape mismatch: {native_bulk_tets.shape} != {src_tets.shape}")
    if native_cohesive.ndim != 2 or native_cohesive.shape[1] != 6:
        raise ValueError(f"native cohesive connectivity shape invalid: {native_cohesive.shape}")
    if native_bulk_tets.size and (np.any(native_bulk_tets < 0) or np.any(native_bulk_tets >= len(native_rest_nodes))):
        raise ValueError("post-insertion bulk connectivity references an absent native node")
    if native_cohesive.size and (np.any(native_cohesive < 0) or np.any(native_cohesive >= len(native_rest_nodes))):
        raise ValueError("cohesive connectivity references an absent native node")

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
                raise ValueError("source-cell/local-ordinal exact rest coordinate consistency failed")
    if set(native_to_input) != set(range(len(native_rest_nodes))):
        raise ValueError("bulk ancestry does not cover every native node")

    source_faces = input_interface_contract(mesh_data, interface_name)
    interface_faces = source_faces["faces"]
    native_face_incidence: dict[tuple[int, int, int], list[tuple[int, int]]] = {}
    for tet_row, tet in enumerate(native_bulk_tets):
        for local_face, ordinals in enumerate(TET_FACES):
            key = tuple(sorted(int(tet[ordinal]) for ordinal in ordinals))
            native_face_incidence.setdefault(key, []).append((tet_row, local_face))
    native_nonmanifold = [key for key, parents in native_face_incidence.items() if len(parents) > 2]
    if native_nonmanifold:
        raise ValueError(
            f"nonmanifold post-insertion bulk face has {len(native_face_incidence[native_nonmanifold[0]])} parents"
        )

    seen_faces: set[tuple[int, int, int]] = set()
    records: list[dict[str, Any]] = []
    for cohesive_row, cohesive in enumerate(native_cohesive):
        sides = (list(map(int, cohesive[:3])), list(map(int, cohesive[3:])))
        if len(set(sides[0])) != 3 or len(set(sides[1])) != 3:
            raise ValueError(f"cohesive row {cohesive_row} contains a repeated side node")
        try:
            source_sides = [tuple(sorted(native_to_input[node] for node in side)) for side in sides]
        except KeyError as error:
            raise ValueError(f"cohesive row {cohesive_row} references a node without bulk ancestry: {error.args[0]}") from error
        if source_sides[0] != source_sides[1]:
            raise ValueError(f"cohesive row {cohesive_row} sides do not pair by exact source node IDs")
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
            parent_centroid_axis = float(src_nodes[src_tets[parent_row], plane_axis].mean())
            signed_distance = parent_centroid_axis - plane_coordinate_m
            if abs(signed_distance) <= 1.0e-12:
                raise ValueError(
                    f"cohesive side {cohesive_row} parent centroid is not more than 1e-12 m from the declared plane"
                )
            side_parents.append({
                "native_bulk_element_row": int(parent_row),
                "input_tetrahedron_id": int(parent_input_id),
                "local_face_ordinal": int(parent_face_ordinal),
                "signed_centroid_distance_m": signed_distance,
                "side": "negative" if signed_distance < 0.0 else "positive",
            })
            areas.append(triangle_area(native_rest_nodes, side_nodes))
        observed_parent_ids = {parent["input_tetrahedron_id"] for parent in side_parents}
        expected_parent_ids = set(source_face["parent_input_tetrahedron_ids"])
        if observed_parent_ids != expected_parent_ids or len(side_parents) != 2:
            raise ValueError("native cohesive side parents do not match the exact source face ancestors")
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
            "side_parent_tetrahedron_ids": [parent["input_tetrahedron_id"] for parent in side_parents],
            "side_parent_by_geometric_side": by_side,
            "post_insertion_bulk_parent_rows": [parent["native_bulk_element_row"] for parent in side_parents],
            "positive_rest_area_m2": areas[0],
            "side_a_positive_rest_area_m2": areas[0],
            "side_b_positive_rest_area_m2": areas[1],
        })
    if seen_faces != set(interface_faces):
        raise ValueError(f"native cohesive/source face bijection failed: {len(seen_faces)} of {len(interface_faces)}")
    total_area = float(sum(record["positive_rest_area_m2"] for record in records))
    if not np.isclose(total_area, expected_area_m2, rtol=0.0, atol=area_tolerance_m2):
        raise ValueError(f"total positive rest interface area {total_area} != {expected_area_m2}")
    return {
        "method": "exact initial MSH/native connectivity + input tetrahedron row/local vertex ordinal + native cohesive connectivity; no proximity node matching",
        "native_node_to_input_node_zero_based": {str(key): int(value) for key, value in sorted(native_to_input.items())},
        "native_bulk_local_ordinal_lineage_rows": mapping_rows,
        "parent_facets": records,
        "input_face_count": len(interface_faces),
        "native_cohesive_element_count": len(native_cohesive),
        "input_rest_area_sum_m2": source_faces["rest_area_sum_m2"],
        "positive_rest_area_sum_m2": total_area,
        "plane_axis": plane_axis,
        "plane_coordinate_m": plane_coordinate_m,
        "fracture_status_claimed": False,
    }


def quadrature_signed_work(
    opening_previous: np.ndarray,
    traction_previous: np.ndarray,
    opening_current: np.ndarray,
    traction_current: np.ndarray,
    positive_area_per_qp_m2: np.ndarray,
) -> np.ndarray:
    """Return raw signed source-formula work for each quadrature point."""
    opening_0 = np.asarray(opening_previous, dtype=np.float64)
    traction_0 = np.asarray(traction_previous, dtype=np.float64)
    opening_1 = np.asarray(opening_current, dtype=np.float64)
    traction_1 = np.asarray(traction_current, dtype=np.float64)
    area = np.asarray(positive_area_per_qp_m2, dtype=np.float64)
    if not (opening_0.shape == traction_0.shape == opening_1.shape == traction_1.shape):
        raise ValueError("opening/traction frames must share shape (n_qp, 3)")
    if opening_0.ndim != 2 or opening_0.shape[1] != 3 or not all(
        np.all(np.isfinite(values)) for values in (opening_0, traction_0, opening_1, traction_1)
    ):
        raise ValueError("opening/traction frames must be finite arrays with shape (n_qp, 3)")
    if area.shape != (len(opening_0),) or not np.all(np.isfinite(area)) or np.any(area <= 0.0):
        raise ValueError("quadrature areas must be finite, positive, and match the number of points")
    return area * 0.5 * np.einsum("qi,qi->q", opening_1 - opening_0, traction_1 + traction_0)
