"""Native tetrahedral connected components and cohesive fragment surfaces.

The extractor works from native T4 and COH3D6 connectivity. Cohesive local
ordinals define the node pairing when no source-node lineage is supplied;
when lineage is supplied, it is used to verify and resolve the pairings.
"""
from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np


# Each tuple is an outward-oriented face for a positively oriented tetrahedron
# (the omitted vertex is opposite the face).
TET_FACES: tuple[tuple[int, int, int], ...] = (
    (1, 2, 3),
    (0, 3, 2),
    (0, 1, 3),
    (0, 2, 1),
)


def _connectivity(value: Any, width: int, name: str, *, allow_empty: bool = False) -> np.ndarray:
    raw = np.asarray(value)
    if allow_empty and raw.size == 0:
        return np.empty((0, width), dtype=np.int64)
    if raw.ndim != 2 or raw.shape[1] != width:
        raise ValueError(f"{name} must have shape (n, {width}), got {raw.shape}")
    try:
        numeric = np.asarray(raw, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must contain integer node indices") from error
    if not np.all(np.isfinite(numeric)) or not np.all(numeric == np.floor(numeric)):
        raise ValueError(f"{name} must contain finite integer node indices")
    return numeric.astype(np.int64)


def _coordinates(value: Any, name: str) -> np.ndarray:
    result = np.asarray(value, dtype=np.float64)
    if result.ndim != 2 or result.shape[1] != 3:
        raise ValueError(f"{name} must have shape (n, 3), got {result.shape}")
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain only finite coordinates")
    return result


def _json_scalar(value: Any, name: str) -> Any:
    if isinstance(value, np.generic):
        value = value.item()
    if not isinstance(value, (str, int, float, bool)) and value is not None:
        raise ValueError(f"{name} values must be scalar strings or numbers")
    if isinstance(value, float) and not np.isfinite(value):
        raise ValueError(f"{name} values must be finite")
    return value


def _labels(value: Any, count: int, name: str) -> list[Any] | None:
    if value is None:
        return None
    rows = list(value)
    if len(rows) != count:
        raise ValueError(f"{name} must contain {count} labels")
    return [_json_scalar(row, name) for row in rows]


def _signed_tet_volume(nodes: np.ndarray, tet: np.ndarray) -> float:
    p0, p1, p2, p3 = nodes[tet]
    return float(np.linalg.det(np.stack((p1 - p0, p2 - p0, p3 - p0))) / 6.0)


def _triangle_area(nodes: np.ndarray, tri: Sequence[int]) -> float:
    p = nodes[np.asarray(tri, dtype=np.int64)]
    return 0.5 * float(np.linalg.norm(np.cross(p[1] - p[0], p[2] - p[0])))


class _DisjointSet:
    def __init__(self, size: int) -> None:
        self.parent = list(range(size))

    def find(self, item: int) -> int:
        root = item
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[item] != item:
            nxt = self.parent[item]
            self.parent[item] = root
            item = nxt
        return root

    def union(self, first: int, second: int) -> None:
        a, b = self.find(first), self.find(second)
        if a == b:
            return
        low, high = sorted((a, b))
        self.parent[high] = low


def _source_exterior_map(
    source_exterior_faces: Any,
) -> tuple[dict[tuple[Any, Any, Any], Any] | None, list[dict[str, Any]] | None]:
    if source_exterior_faces is None:
        return None, None
    if isinstance(source_exterior_faces, Mapping):
        rows = list(source_exterior_faces.items())
    else:
        rows = [(i, face) for i, face in enumerate(source_exterior_faces)]
    by_key: dict[tuple[Any, Any, Any], Any] = {}
    output: list[dict[str, Any]] = []
    for face_id, face_nodes in rows:
        clean_id = _json_scalar(face_id, "source exterior face ID")
        face = list(face_nodes)
        if len(face) != 3:
            raise ValueError("each source exterior face must contain exactly three source node IDs")
        clean_nodes = [_json_scalar(node, "source exterior node ID") for node in face]
        if len(set(clean_nodes)) != 3:
            raise ValueError(f"source exterior face {clean_id!r} has repeated source node IDs")
        key = tuple(sorted(clean_nodes))
        if key in by_key:
            raise ValueError(f"duplicate source exterior face ancestry for source nodes {key!r}")
        by_key[key] = clean_id
        output.append({"source_face_id": clean_id, "source_node_ids": clean_nodes, "source_node_key": list(key)})
    return by_key, output


def _edge_closure(triangles: list[list[int]], node_map: _DisjointSet) -> dict[str, Any]:
    edges: dict[tuple[int, int], list[list[int]]] = defaultdict(list)
    vertices: set[int] = set()
    degenerate_faces: list[int] = []
    for tri_index, tri in enumerate(triangles):
        mapped = [node_map.find(int(node)) for node in tri]
        vertices.update(mapped)
        if len(set(mapped)) != 3:
            degenerate_faces.append(tri_index)
        for i in range(3):
            start, end = mapped[i], mapped[(i + 1) % 3]
            key = (min(start, end), max(start, end))
            direction = 1 if start < end else -1
            edges[key].append([int(start), int(end), direction])

    edge_rows: list[dict[str, Any]] = []
    edge_ok = True
    for key in sorted(edges):
        uses = edges[key]
        reverse_pair = len(uses) == 2 and uses[0][2] == -uses[1][2]
        if not reverse_pair:
            edge_ok = False
        edge_rows.append({
            "topology_edge_node_ids": list(key),
            "incidence_count": len(uses),
            "twice_with_opposite_orientation": bool(reverse_pair),
            "oriented_uses": [[row[0], row[1]] for row in uses],
        })
    closed = edge_ok and not degenerate_faces and bool(triangles)
    euler = int(len(vertices) - len(edges) + len(triangles))
    return {
        "closed": bool(closed),
        "every_edge_twice_in_reverse": bool(closed),
        "vertex_count": len(vertices),
        "edge_count": len(edges),
        "face_count": len(triangles),
        "euler_characteristic": euler,
        "degenerate_face_indices": degenerate_faces,
        "edges": edge_rows,
    }


def extract_cohesive_fragments(
    rest_nodes: Any,
    native_tets: Any,
    native_cohesive: Any,
    facet_damage_qp: Any,
    tractions_qp: Any,
    *,
    damage_tolerance: float,
    traction_tolerance: float,
    current_nodes: Any | None = None,
    native_node_to_source: Any | None = None,
    parent_facet_ids: Any | None = None,
    parent_source_ids: Any | None = None,
    source_exterior_faces: Any | None = None,
    rho: float | None = None,
) -> dict[str, Any]:
    """Extract mechanically connected native tetrahedral fragments.

    A cohesive parent pair remains connected unless *every* quadrature point
    reaches ``damage >= 1 - damage_tolerance`` and has traction norm no greater
    than ``traction_tolerance``. Ordinary T4 faces with two incident tets join
    those tets. Components are formed from all bulk tetrahedra, so volume
    coverage is explicit and no external volume label is used.

    ``source_exterior_faces`` is optional. When supplied, it must be a mapping
    from source face IDs to three source node IDs (or an iterable of triangles,
    whose row indices become IDs), and it is checked as a complete bijection
    against non-cohesive native boundary faces. This requires
    ``native_node_to_source``.
    """
    rest = _coordinates(rest_nodes, "rest_nodes")
    tets = _connectivity(native_tets, 4, "native_tets")
    cohesive = _connectivity(native_cohesive, 6, "native_cohesive", allow_empty=True)
    if len(rest) == 0 or len(tets) == 0:
        raise ValueError("at least one node and one native tetrahedron are required")
    if tets.min() < 0 or tets.max() >= len(rest):
        raise ValueError("native_tets contains a node index outside rest_nodes")
    if cohesive.size and (cohesive.min() < 0 or cohesive.max() >= len(rest)):
        raise ValueError("native_cohesive contains a node index outside rest_nodes")
    if any(len(set(map(int, tet))) != 4 for tet in tets):
        raise ValueError("each native tetrahedron must contain four distinct node IDs")

    deformed = None if current_nodes is None else _coordinates(current_nodes, "current_nodes")
    if deformed is not None and deformed.shape != rest.shape:
        raise ValueError("current_nodes must have the same shape as rest_nodes")

    if not np.isfinite(damage_tolerance) or damage_tolerance < 0.0 or damage_tolerance > 1.0:
        raise ValueError("damage_tolerance must be finite and in [0, 1]")
    if not np.isfinite(traction_tolerance) or traction_tolerance < 0.0:
        raise ValueError("traction_tolerance must be finite and nonnegative")
    if rho is not None and (not np.isfinite(rho) or rho < 0.0):
        raise ValueError("rho must be finite and nonnegative")

    facet_damage = np.asarray(facet_damage_qp, dtype=np.float64)
    tractions = np.asarray(tractions_qp, dtype=np.float64)
    facet_count = len(cohesive)
    if facet_damage.ndim != 2 or facet_damage.shape[0] != facet_count or facet_damage.shape[1] == 0:
        raise ValueError(f"facet_damage_qp must have shape ({facet_count}, Q) with Q > 0")
    if tractions.shape != (*facet_damage.shape, 3):
        raise ValueError(f"tractions_qp must have shape {(*facet_damage.shape, 3)}")
    if not np.all(np.isfinite(facet_damage)) or not np.all(np.isfinite(tractions)):
        raise ValueError("facet damage and traction samples must be finite")
    traction_norms = np.linalg.norm(tractions, axis=2)
    qp_damage_pass = facet_damage >= (1.0 - float(damage_tolerance))
    qp_traction_pass = traction_norms <= float(traction_tolerance)
    qp_joint_pass = qp_damage_pass & qp_traction_pass
    facet_failed = np.all(qp_joint_pass, axis=1)

    source_nodes: list[Any] | None = None
    if native_node_to_source is not None:
        raw_source_nodes = np.asarray(native_node_to_source, dtype=object)
        if raw_source_nodes.shape != (len(rest),):
            raise ValueError(f"native_node_to_source must have shape ({len(rest)},)")
        source_nodes = [_json_scalar(v, "native_node_to_source") for v in raw_source_nodes]
    parent_facets = _labels(parent_facet_ids, facet_count, "parent_facet_ids")
    parent_sources = _labels(parent_source_ids, facet_count, "parent_source_ids")
    if parent_facets is not None:
        try:
            if len(set(parent_facets)) != len(parent_facets):
                raise ValueError("parent_facet_ids must be unique per native cohesive facet")
        except TypeError as error:
            raise ValueError("parent_facet_ids must be hashable scalar labels") from error

    source_exterior_by_key, source_exterior_rows = _source_exterior_map(source_exterior_faces)
    if source_exterior_by_key is not None and source_nodes is None:
        raise ValueError("source exterior face verification requires native_node_to_source")

    rest_volumes = np.asarray([_signed_tet_volume(rest, tet) for tet in tets], dtype=np.float64)
    if np.any(~np.isfinite(rest_volumes)) or np.any(rest_volumes <= 0.0):
        bad = np.flatnonzero((~np.isfinite(rest_volumes)) | (rest_volumes <= 0.0)).tolist()
        raise ValueError(f"rest tetrahedron signed volume must be positive; bad tet rows: {bad}")
    current_volumes: np.ndarray | None = None
    if deformed is not None:
        current_volumes = np.asarray([_signed_tet_volume(deformed, tet) for tet in tets], dtype=np.float64)
        if np.any(~np.isfinite(current_volumes)) or np.any(current_volumes <= 0.0):
            bad = np.flatnonzero((~np.isfinite(current_volumes)) | (current_volumes <= 0.0)).tolist()
            raise ValueError(f"current tetrahedron signed volume must be positive; bad tet rows: {bad}")

    face_incidence: dict[tuple[int, int, int], list[dict[str, Any]]] = defaultdict(list)
    for tet_row, tet in enumerate(tets):
        center = rest[tet].mean(axis=0)
        for local_face, ordinals in enumerate(TET_FACES):
            tri = [int(tet[o]) for o in ordinals]
            points = rest[np.asarray(tri, dtype=np.int64)]
            normal = np.cross(points[1] - points[0], points[2] - points[0])
            face_center = points.mean(axis=0)
            if float(np.dot(normal, face_center - center)) < 0.0:
                tri[1], tri[2] = tri[2], tri[1]
            key = tuple(sorted(tri))
            face_incidence[key].append({
                "tet_row": int(tet_row),
                "local_face_ordinal": int(local_face),
                "native_node_ids": tri,
                "native_node_key": list(key),
            })
    overfull = [(key, len(rows)) for key, rows in face_incidence.items() if len(rows) > 2]
    if overfull:
        raise ValueError(f"non-manifold native tetrahedron face incidence: {overfull[:5]}")

    cohesive_records: list[dict[str, Any]] = []
    cohesive_side_by_key: dict[tuple[int, int, int], dict[str, Any]] = {}
    for facet_row, conn in enumerate(cohesive):
        side_a = [int(v) for v in conn[:3]]
        side_b = [int(v) for v in conn[3:]]
        if len(set(side_a)) != 3 or len(set(side_b)) != 3:
            raise ValueError(f"cohesive facet row {facet_row} has a degenerate side")
        if set(side_a) & set(side_b):
            raise ValueError(f"cohesive facet row {facet_row} sides must use distinct native node IDs")
        side_keys = [tuple(sorted(side_a)), tuple(sorted(side_b))]
        if side_keys[0] == side_keys[1]:
            raise ValueError(f"cohesive facet row {facet_row} repeats the same native face on both sides")
        parents: list[dict[str, int]] = []
        for side_index, key in enumerate(side_keys):
            incident = face_incidence.get(key, [])
            if len(incident) != 1:
                raise ValueError(
                    f"cohesive facet row {facet_row} side {side_index} has {len(incident)} native T4 parents; expected exactly one"
                )
            if key in cohesive_side_by_key:
                raise ValueError(f"native cohesive side face {key} is assigned more than once")
            parents.append({
                "tet_row": int(incident[0]["tet_row"]),
                "local_face_ordinal": int(incident[0]["local_face_ordinal"]),
            })
        if parents[0]["tet_row"] == parents[1]["tet_row"]:
            raise ValueError(f"cohesive facet row {facet_row} has both sides on the same native T4")

        if source_nodes is not None:
            mapped_a = [source_nodes[node] for node in side_a]
            mapped_b = [source_nodes[node] for node in side_b]
            if len(set(mapped_a)) != 3 or len(set(mapped_b)) != 3:
                raise ValueError(f"cohesive facet row {facet_row} has repeated source nodes on a side")
            if set(mapped_a) != set(mapped_b):
                raise ValueError(f"cohesive facet row {facet_row} sides do not pair by source node ancestry")
            b_by_source = {source_nodes[node]: node for node in side_b}
            node_pairs = [[node, b_by_source[source_nodes[node]]] for node in side_a]
        else:
            # COH3D6 local node ordinals 1..3 and 4..6 are corresponding
            # vertices. The source ancestry argument above can verify this.
            node_pairs = [[side_a[i], side_b[i]] for i in range(3)]

        paired_side_a_nodes = [int(pair[0]) for pair in node_pairs]
        paired_side_b_nodes = [int(pair[1]) for pair in node_pairs]
        for node_a, node_b in node_pairs:
            if not np.array_equal(rest[int(node_a)], rest[int(node_b)]):
                raise ValueError(
                    f"cohesive facet row {facet_row} paired rest coordinates do not match exactly"
                )

        is_failed = bool(facet_failed[facet_row])
        # Use the correspondence order on each side. Exact coordinate pairing
        # then makes the two rest-area calculations directly comparable.
        areas_rest = [
            _triangle_area(rest, paired_side_a_nodes),
            _triangle_area(rest, paired_side_b_nodes),
        ]
        if any(not np.isfinite(area) or area <= 0.0 for area in areas_rest):
            raise ValueError(f"cohesive facet row {facet_row} has a nonpositive rest-side area")
        if areas_rest[0] != areas_rest[1]:
            raise ValueError(f"cohesive facet row {facet_row} paired rest-side areas do not match exactly")
        areas_current = None if deformed is None else [
            _triangle_area(deformed, side_a), _triangle_area(deformed, side_b)
        ]
        record: dict[str, Any] = {
            "facet_row": int(facet_row),
            "parent_facet_id": None if parent_facets is None else parent_facets[facet_row],
            "parent_source_id": None if parent_sources is None else parent_sources[facet_row],
            "failed": is_failed,
            "side_a_native_node_ids": side_a,
            "side_b_native_node_ids": side_b,
            "side_a_native_node_key": list(side_keys[0]),
            "side_b_native_node_key": list(side_keys[1]),
            "side_parent_tet_rows": [parents[0]["tet_row"], parents[1]["tet_row"]],
            "side_parent_local_face_ordinals": [parents[0]["local_face_ordinal"], parents[1]["local_face_ordinal"]],
            "paired_native_node_ids": node_pairs,
            "side_a_rest_area": float(areas_rest[0]),
            "side_b_rest_area": float(areas_rest[1]),
            "paired_rest_area_mean": float(sum(areas_rest) / 2.0),
            "side_a_current_area": None if areas_current is None else float(areas_current[0]),
            "side_b_current_area": None if areas_current is None else float(areas_current[1]),
            "paired_current_area_mean": None if areas_current is None else float(sum(areas_current) / 2.0),
            "damage_qp": facet_damage[facet_row].tolist(),
            "traction_norm_qp": traction_norms[facet_row].tolist(),
            "qp_damage_pass": qp_damage_pass[facet_row].tolist(),
            "qp_traction_pass": qp_traction_pass[facet_row].tolist(),
            "qp_joint_failure_pass": qp_joint_pass[facet_row].tolist(),
        }
        cohesive_records.append(record)
        for side_index, key in enumerate(side_keys):
            cohesive_side_by_key[key] = {
                "facet_row": int(facet_row),
                "side_index": int(side_index),
                "parent_tet_row": parents[side_index]["tet_row"],
                "parent_local_face_ordinal": parents[side_index]["local_face_ordinal"],
            }

    adjacency: list[set[int]] = [set() for _ in range(len(tets))]
    for key, incident in face_incidence.items():
        if len(incident) == 2:
            first, second = int(incident[0]["tet_row"]), int(incident[1]["tet_row"])
            adjacency[first].add(second)
            adjacency[second].add(first)
    for record in cohesive_records:
        if not record["failed"]:
            first, second = map(int, record["side_parent_tet_rows"])
            adjacency[first].add(second)
            adjacency[second].add(first)

    component_rows: list[list[int]] = []
    unseen = set(range(len(tets)))
    while unseen:
        seed = min(unseen)
        stack = [seed]
        unseen.remove(seed)
        found: list[int] = []
        while stack:
            row = stack.pop()
            found.append(row)
            for neighbor in sorted(adjacency[row], reverse=True):
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    stack.append(neighbor)
        component_rows.append(sorted(found))
    component_rows.sort(key=lambda rows: rows[0])
    tet_component = {tet_row: comp_index for comp_index, rows in enumerate(component_rows) for tet_row in rows}
    if sorted(tet_component) != list(range(len(tets))) or sum(map(len, component_rows)) != len(tets):
        raise RuntimeError("internal volume coverage partition error")

    node_components: dict[int, set[int]] = defaultdict(set)
    component_node_ids: list[list[int]] = []
    for comp_index, rows in enumerate(component_rows):
        nodes_in_component = sorted({int(node) for row in rows for node in tets[row]})
        component_node_ids.append(nodes_in_component)
        for node in nodes_in_component:
            node_components[node].add(comp_index)
    shared_nodes = {node: sorted(owners) for node, owners in node_components.items() if len(owners) > 1}

    # Verify original-exterior ancestry before constructing component records.
    exterior_source_id_by_face: dict[tuple[int, int, int], Any] = {}
    if source_exterior_by_key is not None:
        seen_source_exterior: dict[tuple[Any, Any, Any], int] = defaultdict(int)
        for key, incident in face_incidence.items():
            if len(incident) != 1 or key in cohesive_side_by_key:
                continue
            source_key = tuple(sorted(source_nodes[node] for node in key))  # type: ignore[index]
            source_face_id = source_exterior_by_key.get(source_key)
            if source_face_id is None:
                raise ValueError(
                    f"native exterior face {key} has no matching source exterior face ancestry (source key {source_key!r})"
                )
            seen_source_exterior[source_key] += 1
            if seen_source_exterior[source_key] > 1:
                raise ValueError(f"source exterior face {source_face_id!r} maps to multiple native boundary faces")
            exterior_source_id_by_face[key] = source_face_id
        missing = set(source_exterior_by_key) - set(seen_source_exterior)
        if missing:
            raise ValueError(f"source exterior faces have no native boundary descendants: {sorted(map(repr, missing))}")

    all_boundary_faces: list[dict[str, Any]] = []
    component_boundary_records: list[list[dict[str, Any]]] = [[] for _ in component_rows]
    for comp_index, rows in enumerate(component_rows):
        for tet_row in rows:
            for local_face, ordinals in enumerate(TET_FACES):
                tet = tets[tet_row]
                oriented = [int(tet[o]) for o in ordinals]
                points = rest[np.asarray(oriented, dtype=np.int64)]
                if float(np.dot(np.cross(points[1] - points[0], points[2] - points[0]), points.mean(axis=0) - rest[tet].mean(axis=0))) < 0.0:
                    oriented[1], oriented[2] = oriented[2], oriented[1]
                key = tuple(sorted(oriented))
                incident = face_incidence[key]
                if len(incident) == 2:
                    continue
                side_info = cohesive_side_by_key.get(key)
                if side_info is not None and not cohesive_records[side_info["facet_row"]]["failed"]:
                    continue
                if side_info is not None:
                    surface_kind = "failed_cohesive_interface"
                    source_face_id = None
                    ancestry_verified = True
                else:
                    surface_kind = "original_exterior"
                    source_face_id = exterior_source_id_by_face.get(key)
                    ancestry_verified = source_exterior_by_key is not None
                face_record = {
                    "global_face_index": len(all_boundary_faces),
                    "component_index": comp_index,
                    "tet_row": int(tet_row),
                    "local_face_ordinal": int(local_face),
                    "native_node_ids": oriented,
                    "native_node_key": list(key),
                    "surface_kind": surface_kind,
                    "facet_row": None if side_info is None else int(side_info["facet_row"]),
                    "facet_side_index": None if side_info is None else int(side_info["side_index"]),
                    "parent_facet_id": None if side_info is None or parent_facets is None else parent_facets[side_info["facet_row"]],
                    "parent_source_id": None if side_info is None or parent_sources is None else parent_sources[side_info["facet_row"]],
                    "source_exterior_face_id": source_face_id,
                    "source_exterior_ancestry_verified": bool(ancestry_verified),
                    "rest_area": _triangle_area(rest, oriented),
                    "current_area": None if deformed is None else _triangle_area(deformed, oriented),
                }
                all_boundary_faces.append(face_record)
                component_boundary_records[comp_index].append(face_record)

    failed_interface_pairs: list[dict[str, Any]] = []
    for record in cohesive_records:
        if not record["failed"]:
            continue
        face_a = next((face for face in all_boundary_faces if face["facet_row"] == record["facet_row"] and face["facet_side_index"] == 0), None)
        face_b = next((face for face in all_boundary_faces if face["facet_row"] == record["facet_row"] and face["facet_side_index"] == 1), None)
        if face_a is None or face_b is None:
            raise ValueError(f"failed cohesive facet row {record['facet_row']} does not map to two boundary faces")
        if (face_a["tet_row"], face_a["local_face_ordinal"]) != (
            record["side_parent_tet_rows"][0], record["side_parent_local_face_ordinals"][0]
        ):
            raise ValueError(f"failed cohesive facet row {record['facet_row']} side A parent incidence changed")
        if (face_b["tet_row"], face_b["local_face_ordinal"]) != (
            record["side_parent_tet_rows"][1], record["side_parent_local_face_ordinals"][1]
        ):
            raise ValueError(f"failed cohesive facet row {record['facet_row']} side B parent incidence changed")
        failed_interface_pairs.append({
            "facet_row": record["facet_row"],
            "parent_facet_id": record["parent_facet_id"],
            "parent_source_id": record["parent_source_id"],
            "side_a_boundary_face_index": face_a["global_face_index"],
            "side_b_boundary_face_index": face_b["global_face_index"],
            "side_a_component_index": face_a["component_index"],
            "side_b_component_index": face_b["component_index"],
            "side_a_native_node_ids": list(face_a["native_node_ids"]),
            "side_b_native_node_ids": list(face_b["native_node_ids"]),
            "side_a_native_node_key": list(face_a["native_node_key"]),
            "side_b_native_node_key": list(face_b["native_node_key"]),
            "side_parent_tet_rows": list(record["side_parent_tet_rows"]),
            "side_parent_local_face_ordinals": list(record["side_parent_local_face_ordinals"]),
            "side_a_rest_area": record["side_a_rest_area"],
            "side_b_rest_area": record["side_b_rest_area"],
            "paired_rest_area_mean": record["paired_rest_area_mean"],
            "side_a_current_area": record["side_a_current_area"],
            "side_b_current_area": record["side_b_current_area"],
            "paired_current_area_mean": record["paired_current_area_mean"],
        })

    components: list[dict[str, Any]] = []
    for comp_index, rows in enumerate(component_rows):
        faces = component_boundary_records[comp_index]
        triangles = [list(face["native_node_ids"]) for face in faces]
        native_nodes = _DisjointSet(len(rest))
        topology_nodes = _DisjointSet(len(rest))
        active_intact_facets: list[int] = []
        for record in cohesive_records:
            if record["failed"] or tet_component[record["side_parent_tet_rows"][0]] != comp_index:
                continue
            if tet_component[record["side_parent_tet_rows"][1]] != comp_index:
                raise RuntimeError("unfailed cohesive parent pair crosses mechanical components")
            active_intact_facets.append(int(record["facet_row"]))
            for node_a, node_b in record["paired_native_node_ids"]:
                topology_nodes.union(int(node_a), int(node_b))
        native_edge_audit = _edge_closure(triangles, native_nodes)
        welded_edge_audit = _edge_closure(triangles, topology_nodes)
        node_ids = component_node_ids[comp_index]
        cross_shared = sorted(node for node in node_ids if node in shared_nodes)
        rest_volume = float(np.sum(rest_volumes[np.asarray(rows, dtype=np.int64)]))
        current_volume = None if current_volumes is None else float(np.sum(current_volumes[np.asarray(rows, dtype=np.int64)]))
        mass = None if rho is None else float(float(rho) * rest_volume)
        reasons: list[str] = []
        if not native_edge_audit["closed"]:
            reasons.append("boundary_edge_incidence_not_closed")
        if active_intact_facets:
            reasons.append("active_cohesive_interface_not_geometrically_closed")
        if cross_shared:
            reasons.append("native_node_shared_with_another_mechanical_component")
        if source_exterior_by_key is None:
            reasons.append("source_exterior_ancestry_unverified")
        eligible = not reasons
        components.append({
            "component_index": comp_index,
            "tet_rows": rows,
            "native_node_ids": node_ids,
            "face_indices": [int(face["global_face_index"]) for face in faces],
            "boundary_faces": faces,
            "oriented_native_triangles": triangles,
            "topology_triangles_after_intact_cohesive_node_pairing": [
                [topology_nodes.find(int(node)) for node in tri] for tri in triangles
            ],
            "closed_edge_audit": native_edge_audit,
            "native_boundary_edge_audit": native_edge_audit,
            "intact_cohesive_welded_edge_audit": welded_edge_audit,
            "closed_solid": bool(native_edge_audit["closed"]),
            "active_intact_cohesive_facet_rows": active_intact_facets,
            "shared_native_node_ids_with_other_components": cross_shared,
            "rest_bulk_volume": rest_volume,
            "current_bulk_volume": current_volume,
            "rho": None if rho is None else float(rho),
            "rest_bulk_mass": mass,
            "qualified_closed_mechanical_fragment": bool(eligible),
            "qualification_reasons": reasons,
            "source_exterior_ancestry_verified": bool(source_exterior_by_key is not None),
        })

    total_rest_volume = float(np.sum(rest_volumes))
    total_current_volume = None if current_volumes is None else float(np.sum(current_volumes))
    summed_component_rest = float(sum(component["rest_bulk_volume"] for component in components))
    if not np.isclose(summed_component_rest, total_rest_volume, rtol=1e-13, atol=0.0):
        raise RuntimeError("component rest volumes do not cover the native bulk volume")
    if current_volumes is not None:
        summed_component_current = float(sum(component["current_bulk_volume"] for component in components))
        if not np.isclose(summed_component_current, total_current_volume, rtol=1e-13, atol=0.0):
            raise RuntimeError("component current volumes do not cover the native bulk volume")

    return {
        "method": "native T4 face incidence plus COH3D6 parent-side incidence and quadrature-wise failure gate",
        "connectivity_rules": {
            "ordinary_shared_native_face": "joins its two tetrahedra",
            "cohesive_parent_pair": "joins unless every QP passes both declared failure thresholds",
            "shared_native_node_without_shared_face": "does not join; disqualifies every affected component",
            "all_native_tets_included": True,
        },
        "failure_thresholds": {
            "damage_tolerance": float(damage_tolerance),
            "damage_threshold_inclusive": float(1.0 - damage_tolerance),
            "traction_norm_tolerance_inclusive": float(traction_tolerance),
            "facet_fails_only_if_all_qps_pass_both": True,
        },
        "rest_nodes": rest.tolist(),
        "current_nodes": None if deformed is None else deformed.tolist(),
        "native_tets": tets.tolist(),
        "native_cohesive": cohesive.tolist(),
        "native_node_to_source": source_nodes,
        "source_exterior_faces": source_exterior_rows,
        "source_exterior_ancestry_provided": bool(source_exterior_by_key is not None),
        "native_tet_count": len(tets),
        "native_cohesive_facet_count": len(cohesive),
        "facet_failure_evidence": cohesive_records,
        "failed_facet_rows": [int(i) for i in np.flatnonzero(facet_failed)],
        "failed_interface_face_pairs": failed_interface_pairs,
        "boundary_faces": all_boundary_faces,
        "components": components,
        "complete_tet_volume_coverage": True,
        "total_rest_bulk_volume": total_rest_volume,
        "total_current_bulk_volume": total_current_volume,
        "rho": None if rho is None else float(rho),
        "total_rest_bulk_mass": None if rho is None else float(float(rho) * total_rest_volume),
        "component_rest_volume_sum_matches_total": True,
        "component_current_volume_sum_matches_total": None if current_volumes is None else True,
        "all_components_qualified_closed_mechanical_fragments": all(
            component["qualified_closed_mechanical_fragment"] for component in components
        ),
        "qualification_scope": "topological and native-connectivity eligibility only; it does not establish material calibration or an R2 result",
    }
