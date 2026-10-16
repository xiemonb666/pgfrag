"""Freeze the v6 mesh and cost forecast for one native zero-step preflight.

The forecast is computed before the native Akantu process is allowed to run.
No solver time steps or history arrays are created by this script.
"""
from __future__ import annotations

import hashlib
import json
import math
import re
import shutil
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[4]
PG = ROOT / "pgfrag"
EXPERIMENT = PG / "experiments/continuation_20261016/reference_v6_preflight"
DATA = PG / "data/continuation_20261016/reference_v6_preflight"
INPUTS = EXPERIMENT / "inputs"
RUN = EXPERIMENT / "native_preflight/run01"
VF_RUN = ROOT / "vesselforge/experiments/continuation_20261016/reference_v6_tetra/run05/partition_output"
VF_RESULT = VF_RUN / "result.json"
VF_MSH = VF_RUN / "mesh.msh"
VF_NPZ = VF_RUN / "mesh.npz"
CANONICAL_STEP = ROOT / "vesselforge/experiments/continuation_20261015/named_step/daily_witness_probe_20261016_v15/daily_native_named.step"
PREVIOUS_INPUTS = PG / "experiments/continuation_20261016/reference_vessel_preflight/run_attempt_02/inputs"
MODULES = ("cohesive_vessel.py", "cohesive_mesh.py", "cohesive_stability.py")
ALL_OUTPUT_CAP = 256 * 1024 * 1024
HISTORY_CAP = 64 * 1024 * 1024
NON_LINEAGE_OUTPUT_ALLOWANCE = 24 * 1024 * 1024
FORECAST_REPORT_ALLOWANCE = 2 * 1024 * 1024


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def tree_bytes(path: Path) -> int:
    return sum(item.stat().st_size for item in path.rglob("*") if item.is_file())


def section_span(data: bytes, name: str) -> tuple[int, int, bytes]:
    pattern = re.compile(rb"(?m)^\$" + name.encode("ascii") + rb"\r?\n.*?^\$End" + name.encode("ascii") + rb"\r?$", re.S)
    match = pattern.search(data)
    if not match:
        raise ValueError(f"missing or malformed ${name} section")
    return match.start(), match.end(), match.group(0)


def raw_element_records(data: bytes) -> list[dict[str, object]]:
    _, _, block = section_span(data, "Elements")
    rows = block.splitlines()[1:-1]
    if not rows:
        raise ValueError("empty Elements block")
    count = int(rows[0])
    parsed: list[dict[str, object]] = []
    for row_ordinal, row in enumerate(rows[1:]):
        fields = [int(token) for token in row.split()]
        if len(fields) < 3:
            raise ValueError(f"malformed MSH element row {row_ordinal}")
        element_id, element_type, tag_count = fields[:3]
        tags = fields[3:3 + tag_count]
        node_tags = fields[3 + tag_count:]
        if tag_count < 0 or len(fields) != 3 + tag_count + len(node_tags):
            raise ValueError(f"malformed MSH tags on element {element_id}")
        parsed.append({
            "row_ordinal": row_ordinal,
            "element_id": element_id,
            "element_type": element_type,
            "tags": tags,
            "physical_tag": tags[0] if tags else None,
            "geometrical_entity_tag": tags[1] if len(tags) > 1 else None,
            "node_tags": node_tags,
        })
    if len(parsed) != count:
        raise ValueError(f"MSH element count mismatch: {len(parsed)} != {count}")
    return parsed


def mesh_alias_and_fid_map(source_msh: Path, result: dict[str, object]) -> dict[str, object]:
    original = source_msh.read_bytes()
    names_start, names_end, names_block = section_span(original, "PhysicalNames")
    replacements = {
        b'2 2 "cut_interface"': b'2 2 "interface"',
        b'2 3 "source_exterior"': b'2 3 "exterior"',
    }
    aliased_names = names_block
    replacement_counts: dict[str, int] = {}
    for old, new in replacements.items():
        count = aliased_names.count(old)
        if count != 1:
            raise ValueError(f"expected one physical-name row {old!r}; found {count}")
        aliased_names = aliased_names.replace(old, new, 1)
        replacement_counts[old.decode("ascii")] = count
    alias = original[:names_start] + aliased_names + original[names_end:]
    derived_path = INPUTS / "mesh.msh"
    derived_path.write_bytes(alias)

    src_nodes_start, src_nodes_end, src_nodes_block = section_span(original, "Nodes")
    dst_nodes_start, dst_nodes_end, dst_nodes_block = section_span(alias, "Nodes")
    src_elements_start, src_elements_end, src_elements_block = section_span(original, "Elements")
    dst_elements_start, dst_elements_end, dst_elements_block = section_span(alias, "Elements")
    if src_nodes_block != dst_nodes_block or src_elements_block != dst_elements_block:
        raise ValueError("raw MSH Nodes or Elements blocks changed during physical-name aliasing")
    if original[:names_start] != alias[:names_start] or original[names_end:] != alias[names_start + len(aliased_names):]:
        raise ValueError("bytes outside PhysicalNames changed during physical-name aliasing")

    sys.path.insert(0, str(INPUTS))
    from cohesive_mesh import input_interface_contract, read_gmsh22  # noqa: E402

    parsed_original = read_gmsh22(source_msh)
    parsed_alias = read_gmsh22(derived_path)
    if not np.array_equal(parsed_original["nodes"], parsed_alias["nodes"]):
        raise ValueError("parsed MSH node coordinates/order changed in alias")
    if parsed_original["node_tags"] != parsed_alias["node_tags"]:
        raise ValueError("parsed MSH node tag order changed in alias")

    def element_signature(rows: list[dict[str, object]]) -> list[tuple[object, ...]]:
        return [(
            int(row["element_id"]), int(row["element_type"]),
            int(row["physical_tag"]) if row["physical_tag"] is not None else None,
            tuple(int(tag) for tag in row["node_tags"]),
            tuple(int(node) for node in row["nodes"]),
        ) for row in rows]

    if element_signature(parsed_original["elements"]) != element_signature(parsed_alias["elements"]):
        raise ValueError("parsed MSH connectivity, physical tags, or element order changed in alias")
    expected_original_names = {
        (2, 2): "cut_interface", (2, 3): "source_exterior", (3, 1): "partitioned_body",
    }
    expected_alias_names = {
        (2, 2): "interface", (2, 3): "exterior", (3, 1): "partitioned_body",
    }
    if parsed_original["physical_names"] != expected_original_names:
        raise ValueError(f"unexpected source MSH physical names: {parsed_original['physical_names']}")
    if parsed_alias["physical_names"] != expected_alias_names:
        raise ValueError(f"unexpected alias MSH physical names: {parsed_alias['physical_names']}")
    if parsed_original["physical_names"][(3, 1)] != parsed_alias["physical_names"][(3, 1)]:
        raise ValueError("volume physical group changed")

    raw_original = raw_element_records(original)
    raw_alias = raw_element_records(alias)
    raw_by_id = {int(row["element_id"]): row for row in raw_original}
    if raw_original != raw_alias:
        raise ValueError("raw element IDs, tag IDs, node tags, or row order changed in alias")

    with np.load(INPUTS / "source_mesh.npz", allow_pickle=False) as source_npz:
        npz_nodes = np.asarray(source_npz["nodes_m"], dtype=np.float64)
        npz_node_tags = np.asarray(source_npz["original_node_tags"], dtype=np.int64)
        npz_tets = np.asarray(source_npz["tets"], dtype=np.int64)
        npz_cell_tags = np.asarray(source_npz["original_cell_tags"], dtype=np.int64)
        npz_interface = np.asarray(source_npz["interface_triangles"], dtype=np.int64)
        npz_exterior = np.asarray(source_npz["exterior_triangles"], dtype=np.int64)
        npz_exterior_entities = np.asarray(source_npz["exterior_surface_tags"], dtype=np.int64)
        npz_exterior_fids = np.asarray(source_npz["exterior_source_face_ids"]).astype(str)
        npz_exterior_roles = np.asarray(source_npz["exterior_roles"]).astype(str)
        npz_exterior_exact = np.asarray(source_npz["exterior_ancestry_exact"], dtype=bool)
    if npz_nodes.ndim != 2 or npz_nodes.shape[1] != 3 or npz_node_tags.shape != (len(npz_nodes),):
        raise ValueError("run05 mesh.npz node coordinates/tags have unexpected shapes")
    if (len(np.unique(npz_node_tags)) != len(npz_node_tags) or npz_tets.ndim != 2 or npz_tets.shape[1] != 4
            or npz_cell_tags.shape != (len(npz_tets),) or len(np.unique(npz_cell_tags)) != len(npz_cell_tags)):
        raise ValueError("run05 mesh.npz node tags, original T4 element tags, or connectivity are invalid")
    if np.any(npz_tets < 0) or np.any(npz_tets >= len(npz_nodes)):
        raise ValueError("run05 mesh.npz T4 connectivity references an absent node")
    npz_referenced_ordinals = np.unique(npz_tets)
    npz_unused_ordinals = np.setdiff1d(np.arange(len(npz_nodes), dtype=np.int64), npz_referenced_ordinals)
    npz_unused_tags = npz_node_tags[npz_unused_ordinals]
    msh_tags = np.asarray(parsed_original["node_tags"], dtype=np.int64)
    unused_present_in_msh = sorted(set(npz_unused_tags.tolist()) & set(msh_tags.tolist()))
    if set(msh_tags.tolist()) != set(npz_node_tags[npz_referenced_ordinals].tolist()):
        raise ValueError("MSH node tags do not equal the run05 NPZ T4-referenced node tags")
    npz_row_by_tag = {int(tag): row for row, tag in enumerate(npz_node_tags)}
    msh_row_by_tag = {int(tag): row for row, tag in enumerate(msh_tags)}
    max_coordinate_abs_error = 0.0
    max_coordinate_relative_error = 0.0
    for tag, msh_row in msh_row_by_tag.items():
        npz_row = npz_row_by_tag[tag]
        msh_xyz = np.asarray(parsed_original["nodes"][msh_row], dtype=np.float64)
        npz_xyz = np.asarray(npz_nodes[npz_row], dtype=np.float64)
        delta = np.abs(msh_xyz - npz_xyz)
        max_coordinate_abs_error = max(max_coordinate_abs_error, float(delta.max(initial=0.0)))
        reference_norm = float(np.linalg.norm(npz_xyz))
        if reference_norm > 0.0:
            max_coordinate_relative_error = max(max_coordinate_relative_error, float(np.linalg.norm(delta) / reference_norm))
        if not np.allclose(msh_xyz, npz_xyz, rtol=8e-16, atol=1e-30):
            raise ValueError(f"MSH and run05 NPZ coordinates exceed declared serialization tolerance for original node tag {tag}")
    npz_tet_row_by_cell_tag = {int(tag): row for row, tag in enumerate(npz_cell_tags)}
    msh_tet_element_ids = np.asarray([int(row["element_id"]) for row in parsed_original["tetrahedra"]], dtype=np.int64)
    if set(msh_tet_element_ids.tolist()) != set(npz_cell_tags.tolist()):
        raise ValueError("MSH T4 element IDs do not equal run05 NPZ original_cell_tags")
    msh_to_npz_tet_rows = np.asarray([npz_tet_row_by_cell_tag[int(tag)] for tag in msh_tet_element_ids], dtype=np.int64)
    npz_mapped_tet_corner_tags = npz_node_tags[npz_tets[msh_to_npz_tet_rows]]
    msh_tet_corner_tags = np.asarray([row["node_tags"] for row in parsed_original["tetrahedra"]], dtype=np.int64)
    if not np.array_equal(npz_mapped_tet_corner_tags, msh_tet_corner_tags):
        bad_row = int(np.flatnonzero(np.any(npz_mapped_tet_corner_tags != msh_tet_corner_tags, axis=1))[0])
        raise ValueError(f"MSH T4 original element ID {msh_tet_element_ids[bad_row]} has different ordered corner nodeTags in run05 NPZ")
    cell_row_map_path = INPUTS / "t4_cell_row_correspondence.npz"
    np.savez_compressed(
        cell_row_map_path,
        msh_tet_raw_row_ordinals=np.asarray([raw_by_id[int(row["element_id"])]["row_ordinal"] for row in parsed_original["tetrahedra"]], dtype=np.int64),
        msh_tet_element_ids=msh_tet_element_ids,
        npz_tet_rows_by_original_cell_tag=msh_to_npz_tet_rows,
        ordered_original_node_tags_by_corner=msh_tet_corner_tags,
    )
    vf_ancestry = result["source_face_ancestry"]
    if vf_ancestry["status"] != "VERIFIED_OCC_WITNESS_ANCESTRY":
        raise ValueError("run05 does not report verified OCC witness ancestry")
    cut = vf_ancestry["cut_witness_ancestry"]
    if cut["status"] != "VERIFIED_OCC_WITNESS_ANCESTRY" or cut["interface_has_zero_source_face_ancestry"] is not True:
        raise ValueError("run05 cut ancestry does not verify unique exterior parents and an unowned interface")
    entity_map = {int(key): value for key, value in cut["source_by_exterior_surface"].items()}
    expected_surface_tags = {int(tag) for tag in cut["exterior_surface_tags"]}
    if set(entity_map) != expected_surface_tags:
        raise ValueError("run05 exterior entity mapping does not match its declared exterior surface tags")

    fid_ids = sorted({str(value["source_face_id"]) for value in entity_map.values()})
    fid_index = {fid: index for index, fid in enumerate(fid_ids)}
    role_names = sorted({str(value["role"]) for value in entity_map.values()})
    role_index = {role: index for index, role in enumerate(role_names)}
    alias_triangles = {int(row["element_id"]): row for row in parsed_alias["triangles"]}
    exterior_rows = [row for row in parsed_alias["triangles"] if row["physical_tag"] == 3]
    if len(exterior_rows) != int(result["counts"]["exterior_triangles"]):
        raise ValueError("MSH exterior triangle count differs from run05 result")
    if len([row for row in parsed_alias["triangles"] if row["physical_tag"] == 2]) != int(result["counts"]["interface_triangles"]):
        raise ValueError("MSH interface triangle count differs from run05 result")
    if len(parsed_alias["tetrahedra"]) != int(result["counts"]["tetrahedra"]):
        raise ValueError("MSH T4 count differs from run05 result")
    npz_ext_by_key: dict[tuple[int, int, int], int] = {}
    for npz_row, triangle in enumerate(npz_exterior):
        if triangle.shape != (3,) or np.any(triangle < 0) or np.any(triangle >= len(npz_node_tags)):
            raise ValueError(f"run05 NPZ exterior triangle row {npz_row} has invalid node ordinals")
        key = tuple(sorted(int(npz_node_tags[node]) for node in triangle))
        if key in npz_ext_by_key:
            raise ValueError("run05 NPZ contains duplicate exterior sourceTri node-tag keys")
        npz_ext_by_key[key] = npz_row
    if len(npz_exterior) != len(exterior_rows) or len(npz_exterior_entities) != len(exterior_rows) or len(npz_exterior_fids) != len(exterior_rows) or len(npz_exterior_roles) != len(exterior_rows) or len(npz_exterior_exact) != len(exterior_rows):
        raise ValueError("run05 NPZ exterior identity arrays have inconsistent lengths")
    if not bool(npz_exterior_exact.all()):
        raise ValueError("run05 NPZ contains exterior triangles without exact ancestry")
    npz_interface_keys = {
        tuple(sorted(int(npz_node_tags[node]) for node in triangle))
        for triangle in npz_interface
    }
    msh_interface_keys = {
        tuple(sorted(map(int, row["node_tags"])))
        for row in parsed_alias["triangles"] if row["physical_tag"] == 2
    }
    if npz_interface_keys != msh_interface_keys:
        raise ValueError("run05 NPZ and MSH interface sourceTri node-tag key sets differ")

    raw_exterior = []
    for row in exterior_rows:
        raw = raw_by_id[int(row["element_id"])]
        entity_tag = raw["geometrical_entity_tag"]
        if entity_tag is None or int(entity_tag) not in entity_map:
            raise ValueError(f"exterior triangle {row['element_id']} has no verified source-face entity mapping")
        face = entity_map[int(entity_tag)]
        if tuple(sorted(map(int, row["node_tags"]))) != tuple(sorted(map(int, raw["node_tags"]))):
            raise ValueError(f"raw/parser node tags differ on exterior element {row['element_id']}")
        key = tuple(sorted(int(tag) for tag in raw["node_tags"]))
        npz_row = npz_ext_by_key.get(key)
        if npz_row is None:
            raise ValueError(f"MSH exterior sourceTri key has no exact run05 NPZ row: {key}")
        if (int(npz_exterior_entities[npz_row]) != int(entity_tag)
                or str(npz_exterior_fids[npz_row]) != str(face["source_face_id"])
                or str(npz_exterior_roles[npz_row]) != str(face["role"])):
            raise ValueError(f"run05 NPZ FID/role/entity binding conflicts with OCC witness ancestry for {key}")
        raw_exterior.append({
            "element_id": int(row["element_id"]),
            "row_ordinal": int(raw["row_ordinal"]),
            "physical_tag": 3,
            "geometrical_entity_tag": int(entity_tag),
            "node_tags": [int(tag) for tag in raw["node_tags"]],
            "source_tri_key_node_tags_sorted": sorted(int(tag) for tag in raw["node_tags"]),
            "source_surface_tag": int(face["source_surface_tag"]),
            "source_face_id": str(face["source_face_id"]),
            "role": str(face["role"]),
        })
    observed_entity_tags = {int(row["geometrical_entity_tag"]) for row in raw_exterior}
    if observed_entity_tags != expected_surface_tags:
        raise ValueError("not every OCC exterior entity tag is represented in the MSH exterior triangles")
    if len({tuple(row["source_tri_key_node_tags_sorted"]) for row in raw_exterior}) != len(raw_exterior):
        raise ValueError("duplicate exact sourceTri node-tag key in MSH exterior triangles")
    if {tuple(row["source_tri_key_node_tags_sorted"]) for row in raw_exterior} != set(npz_ext_by_key):
        raise ValueError("MSH and NPZ exterior sourceTri node-tag key sets differ")

    np.savez_compressed(
        INPUTS / "known_exterior_fid_roles.npz",
        source_exterior_element_ids=np.asarray([row["element_id"] for row in raw_exterior], dtype=np.int64),
        source_exterior_raw_element_row_ordinals=np.asarray([row["row_ordinal"] for row in raw_exterior], dtype=np.int64),
        source_exterior_physical_tags=np.asarray([row["physical_tag"] for row in raw_exterior], dtype=np.int64),
        source_exterior_geometrical_entity_tags=np.asarray([row["geometrical_entity_tag"] for row in raw_exterior], dtype=np.int64),
        source_exterior_node_tags=np.asarray([row["node_tags"] for row in raw_exterior], dtype=np.int64),
        source_exterior_source_tri_key_node_tags_sorted=np.asarray([row["source_tri_key_node_tags_sorted"] for row in raw_exterior], dtype=np.int64),
        source_exterior_source_surface_tags=np.asarray([row["source_surface_tag"] for row in raw_exterior], dtype=np.int64),
        source_exterior_fid_index=np.asarray([fid_index[row["source_face_id"]] for row in raw_exterior], dtype=np.int16),
        source_exterior_role_index=np.asarray([role_index[row["role"]] for row in raw_exterior], dtype=np.int8),
        source_face_ids=np.asarray(fid_ids, dtype="U64"),
        source_face_roles=np.asarray([next(str(value["role"]) for value in entity_map.values() if value["source_face_id"] == fid) for fid in fid_ids], dtype="U32"),
        role_names=np.asarray(role_names, dtype="U32"),
    )

    alias_receipt = {
        "schema": "pgfrag.physical_name_alias_receipt.v1",
        "operation": "REPLACE_TWO_PHYSICAL_NAME_STRINGS_ONLY",
        "source_mesh": {"path": "mesh.original.msh", "bytes": len(original), "sha256": hashlib.sha256(original).hexdigest()},
        "executed_mesh": {"path": "mesh.msh", "bytes": len(alias), "sha256": hashlib.sha256(alias).hexdigest()},
        "physical_name_replacements": replacement_counts,
        "physical_names_original": [{"dimension": d, "physical_tag": t, "name": n} for (d, t), n in expected_original_names.items()],
        "physical_names_derived": [{"dimension": d, "physical_tag": t, "name": n} for (d, t), n in expected_alias_names.items()],
        "nodes_block_byte_equal": src_nodes_block == dst_nodes_block,
        "elements_block_byte_equal": src_elements_block == dst_elements_block,
        "all_bytes_outside_physical_names_byte_equal": original[:names_start] == alias[:names_start] and original[names_end:] == alias[names_start + len(aliased_names):],
        "raw_element_records_and_order_equal": raw_original == raw_alias,
        "physical_integer_tags_and_order_equal": all(a["tags"] == b["tags"] and a["element_id"] == b["element_id"] and a["row_ordinal"] == b["row_ordinal"] for a, b in zip(raw_original, raw_alias, strict=True)),
        "parsed_node_coordinates_and_order_equal": True,
        "parsed_node_tag_order_equal": True,
        "parsed_connectivity_element_ids_physical_tags_and_order_equal": True,
        "volume_physical_group_unchanged": True,
        "source_node_count": len(parsed_original["nodes"]),
        "source_tetrahedron_count": len(parsed_original["tetrahedra"]),
        "source_triangle_count": len(parsed_original["triangles"]),
        "source_interface_triangle_count": len([row for row in parsed_original["triangles"] if row["physical_tag"] == 2]),
        "source_exterior_triangle_count": len(exterior_rows),
        "known_exterior_fid_role_rows": len(raw_exterior),
        "known_exterior_fid_count": len(fid_ids),
        "run05_npz_total_nodes": len(npz_nodes),
        "run05_npz_tet_referenced_nodes": len(npz_referenced_ordinals),
        "run05_npz_unused_node_count": len(npz_unused_ordinals),
        "run05_npz_unused_original_node_tags": [int(tag) for tag in npz_unused_tags],
        "run05_npz_unused_tags_present_in_msh": [int(tag) for tag in unused_present_in_msh],
        "run05_npz_to_msh_coordinate_roundtrip_max_abs_error_m": max_coordinate_abs_error,
        "run05_npz_to_msh_coordinate_roundtrip_max_relative_error": max_coordinate_relative_error,
        "run05_npz_to_msh_coordinate_roundtrip_tolerance": {"rtol": 8e-16, "atol_m": 1e-30},
        "msh_node_tags_equal_npz_tet_referenced_node_tags": True,
        "run05_npz_and_msh_tetra_row_local_order_equal": False,
        "run05_npz_and_msh_tetra_correspondence": "EXACT_ORIGINAL_CELL_TAG_AND_ORDERED_CORNER_NODE_TAGS",
        "run05_npz_tet_original_cell_tags_unique": True,
        "run05_npz_msh_tet_exact_cell_tag_set_equal": True,
        "run05_npz_tet_local_corner_order_matches_msh_after_original_cell_tag_map": True,
        "run05_npz_tet_original_cell_tag_mapping_is_identity": bool(np.array_equal(msh_to_npz_tet_rows, np.arange(len(msh_to_npz_tet_rows)))),
        "msh_t4_to_npz_row_correspondence_path": "t4_cell_row_correspondence.npz",
        "msh_t4_to_npz_row_correspondence_sha256": sha256(cell_row_map_path),
        "msh_t4_to_npz_row_correspondence_count": len(msh_tet_element_ids),
        "msh_t4_to_npz_row_correspondence_is_identity": bool(np.array_equal(msh_to_npz_tet_rows, np.arange(len(msh_to_npz_tet_rows)))),
        "run05_npz_and_msh_interface_sourceTri_key_sets_equal": True,
        "run05_npz_and_msh_exterior_sourceTri_key_sets_equal": True,
        "run05_npz_all_exterior_ancestry_exact": True,
        "known_exterior_fid_mapping_path": "known_exterior_fid_roles.npz",
        "known_exterior_fid_mapping_sha256": sha256(INPUTS / "known_exterior_fid_roles.npz"),
        "fid_entity_mapping": [
            {"geometrical_entity_tag": tag, "source_surface_tag": int(value["source_surface_tag"]),
             "source_face_id": str(value["source_face_id"]), "role": str(value["role"])}
            for tag, value in sorted(entity_map.items())
        ],
        "exact_binding_method": "MSH2 physical/geometrical tags plus exact sorted original nodeTags sourceTri keys; no coordinate or nearest matching",
    }
    write_json(INPUTS / "mesh_alias_receipt.json", alias_receipt)
    # Keep a separate row table so reviewers can inspect the exact 44,182 bindings.
    write_json(INPUTS / "known_exterior_fid_role_rows.json", {
        "schema": "pgfrag.vf_v6_exact_msh_source_exterior_fid_role_rows.v1",
        "source_mesh_sha256": hashlib.sha256(original).hexdigest(),
        "source_mesh_result_sha256": sha256(VF_RESULT),
        "source_face_lineage_sha256": vf_ancestry["source_face_lineage_sha256"],
        "canonical_step_sha256": vf_ancestry["canonical_step_sha256"],
        "upstream_source_step_sha256": vf_ancestry["upstream_source_step_sha256"],
        "rows": raw_exterior,
    })
    return {
        "source": parsed_original,
        "alias": parsed_alias,
        "alias_receipt": alias_receipt,
        "entity_map": entity_map,
        "fid_ids": fid_ids,
        "role_names": role_names,
        "exterior_rows": raw_exterior,
        "npz_audit": {
            "npz_nodes": len(npz_nodes),
            "npz_tet_referenced_nodes": len(npz_referenced_ordinals),
            "npz_unused_ordinals": [int(row) for row in npz_unused_ordinals],
            "npz_unused_original_node_tags": [int(tag) for tag in npz_unused_tags],
            "npz_unused_tags_present_in_msh": [int(tag) for tag in unused_present_in_msh],
            "max_coordinate_abs_error_m": max_coordinate_abs_error,
            "max_coordinate_relative_error": max_coordinate_relative_error,
            "coordinate_comparison": "np.allclose(rtol=8e-16, atol=1e-30); MSH ASCII serialization can differ from NPZ binary64 by rounding",
            "msh_node_tags_equal_npz_tet_referenced_tags": True,
            "msh_tets_equal_npz_rows": False,
            "msh_tet_element_ids_equal_npz_original_cell_tags": True,
            "msh_to_npz_t4_row_permutation_identity": bool(np.array_equal(msh_to_npz_tet_rows, np.arange(len(msh_to_npz_tet_rows)))),
            "msh_tet_to_npz_tet_rows_by_original_cell_tag": msh_to_npz_tet_rows.tolist(),
            "t4_cell_row_correspondence_path": "t4_cell_row_correspondence.npz",
            "t4_cell_row_correspondence_sha256": sha256(cell_row_map_path),
            "interface_and_exterior_exact_key_sets_equal": True,
        },
    }


def synthetic_lineage_forecast(mesh: dict[str, object], rho: float) -> tuple[int, dict[str, object]]:
    """Build the exact lineage JSON shape with a deterministic synthetic seam split."""
    from cohesive_mesh import build_native_lineage, input_interface_contract

    nodes = np.asarray(mesh["nodes"], dtype=np.float64)
    source_tets = np.asarray([row["nodes"] for row in mesh["tetrahedra"]], dtype=np.int64)
    source_tet_ids = [int(row["element_id"]) for row in mesh["tetrahedra"]]
    contract = input_interface_contract(mesh, interface_name="interface")
    interface_nodes = sorted({int(node) for face in contract["faces"].values() for node in face["input_face_node_ids"]})
    duplicate = {source_node: len(nodes) + index for index, source_node in enumerate(interface_nodes)}
    rest = np.vstack((nodes, nodes[np.asarray(interface_nodes, dtype=np.int64)]))
    centroids_x = nodes[source_tets, 0].mean(axis=1)
    if np.any(np.abs(centroids_x) <= 1e-12):
        raise ValueError("synthetic lineage forecast found a T4 centroid on the interface plane")
    positive_mask = centroids_x > 0.0
    native_tets = source_tets.copy()
    for tet_row in np.flatnonzero(positive_mask):
        for local_vertex in range(4):
            source_node = int(source_tets[tet_row, local_vertex])
            if source_node in duplicate:
                native_tets[tet_row, local_vertex] = duplicate[source_node]
    native_cohesive: list[list[int]] = []
    for face in contract["faces"].values():
        parent_rows = face["parent_input_tetrahedron_rows"]
        positive_parents = [row for row in parent_rows if positive_mask[row]]
        negative_parents = [row for row in parent_rows if not positive_mask[row]]
        if len(positive_parents) != 1 or len(negative_parents) != 1:
            raise ValueError("synthetic lineage forecast could not identify opposite-side interface parents")
        source_nodes = [int(node) for node in face["input_face_node_ids"]]
        side_positive = [duplicate.get(node, node) for node in source_nodes]
        native_cohesive.append(side_positive + source_nodes)
    lineage = build_native_lineage(
        mesh, rest, native_tets, np.asarray(native_cohesive, dtype=np.int64),
        native_initial_bulk_tets=source_tets,
        expected_area_m2=float(contract["rest_area_sum_m2"]),
        plane_axis=0, plane_coordinate_m=0.0,
    )
    serialized = (json.dumps(lineage, indent=2, allow_nan=False) + "\n").encode("utf-8")
    # Actual Akantu row and cohesive order can change numeric values, but not the
    # number of records or their bounded integer digit widths. Add a 10% safety
    # margin plus 1 MiB for any serializer/order-dependent length difference.
    forecast_bytes = int(math.ceil(len(serialized) * 1.10)) + 1024 * 1024
    return forecast_bytes, {
        "synthetic_exact_lineage_bytes": len(serialized),
        "lineage_forecast_bytes_with_10pct_plus_1MiB_margin": forecast_bytes,
        "synthetic_input_node_count": len(nodes),
        "synthetic_native_node_count": len(rest),
        "synthetic_bulk_T4_count": len(source_tets),
        "synthetic_cohesive_Tri3_count": len(native_cohesive),
        "synthetic_interface_node_count_duplicated": len(interface_nodes),
        "all_lineage_rows_use_input_tet_row_and_local_vertex_order": True,
        "synthetic_material_density_kg_m3": rho,
    }


def main() -> int:
    if INPUTS.exists() or RUN.exists() or (DATA / "preflight_budget.json").exists() or (DATA / "run_control.json").exists():
        raise FileExistsError("frozen inputs or a native attempt already exist; exactly one attempt is permitted")
    for path in (EXPERIMENT, DATA):
        path.mkdir(parents=True, exist_ok=True)
    INPUTS.mkdir()

    current_modules = {name: PG / "src/pgfrag" / name for name in MODULES}
    for name, source in current_modules.items():
        shutil.copyfile(source, INPUTS / name)
    shutil.copyfile(PREVIOUS_INPUTS / "material.dat", INPUTS / "material.dat")
    shutil.copyfile(VF_MSH, INPUTS / "mesh.original.msh")
    shutil.copyfile(VF_RESULT, INPUTS / "source_mesh_result.json")
    shutil.copyfile(VF_NPZ, INPUTS / "source_mesh.npz")

    result = json.loads((INPUTS / "source_mesh_result.json").read_text(encoding="utf-8"))
    if result.get("schema") != "vesselforge.generic_step_tetra_partition.v1" or result.get("status") != "GEOMETRY_DIAGNOSTIC_ONLY":
        raise ValueError("unexpected run05 source result schema or status")
    if sha256(CANONICAL_STEP) != result["source_step_sha256"]:
        raise ValueError("canonical v15 STEP does not match run05 source hash")
    if sha256(VF_MSH) != result["files"]["mesh.msh"]["sha256"] or sha256(VF_NPZ) != result["files"]["mesh.npz"]["sha256"]:
        raise ValueError("run05 source mesh hashes do not match the result manifest")
    if result["source_face_ancestry"]["source_face_lineage_sha256"] != "9e6c789989cb0f5f42ed9c37f0264bc02899a9c29f571254938a4781640bcfd2":
        raise ValueError("run05 source face lineage hash differs from the assigned v6 record")
    if result["source_face_ancestry"]["upstream_source_step_sha256"] != "1cd34bc5f5b3cf8b6b937afeb7345b982c9de01165f0d570ca11b83f3ea31c0e":
        raise ValueError("run05 upstream original v6 STEP hash differs from the assigned record")

    mesh_audit = mesh_alias_and_fid_map(INPUTS / "mesh.original.msh", result)
    source_mesh = mesh_audit["source"]
    raw_alias_hash = sha256(INPUTS / "mesh.msh")
    alias_receipt_hash = sha256(INPUTS / "mesh_alias_receipt.json")
    fid_roles_hash = sha256(INPUTS / "known_exterior_fid_roles.npz")
    fid_rows_hash = sha256(INPUTS / "known_exterior_fid_role_rows.json")

    fid_role_records = []
    cut_map = result["source_face_ancestry"]["cut_witness_ancestry"]["source_by_exterior_surface"]
    for source_face_id in mesh_audit["fid_ids"]:
        descendants = [
            int(entity_tag) for entity_tag, row in cut_map.items()
            if str(row["source_face_id"]) == source_face_id
        ]
        face_record = next(row for row in cut_map.values() if str(row["source_face_id"]) == source_face_id)
        fid_role_records.append({
            "source_face_id": source_face_id,
            "role": str(face_record["role"]),
            "source_surface_tag": int(face_record["source_surface_tag"]),
            "exterior_output_entity_tags": sorted(descendants),
        })
    provenance = {
        "schema": "pgfrag.reference_v6_native_preflight_geometry_provenance.v1",
        "classification": "GEOMETRY_DIAGNOSTIC_ONLY; REFERENCE CANDIDATE, NOT MEASURED SOURCE PROFILE",
        "candidate": "daily_aberdeen_108x53_v6_named_step_wall2mm_reference",
        "vessel_reference": "Daily Aberdeen 108x53 reference candidate; wall profile is not measured or physically qualified",
        "canonical_step": {
            "path": "vesselforge/experiments/continuation_20261015/named_step/daily_witness_probe_20261016_v15/daily_native_named.step",
            "bytes": CANONICAL_STEP.stat().st_size,
            "sha256": sha256(CANONICAL_STEP),
            "declared_length_unit": result["source_step_declared_length_unit"],
            "unit_conversion": result["step_unit_conversion"],
        },
        "upstream_original_v6_step_sha256": result["source_face_ancestry"]["upstream_source_step_sha256"],
        "source_face_lineage_sha256": result["source_face_ancestry"]["source_face_lineage_sha256"],
        "witness_manifest_sha256": result["source_face_ancestry"]["manifest_sha256"],
        "mesh_sources": {
            "mesh_msh": {"path": "mesh.original.msh", "bytes": (INPUTS / "mesh.original.msh").stat().st_size, "sha256": sha256(INPUTS / "mesh.original.msh")},
            "mesh_npz": {"path": "source_mesh.npz", "bytes": (INPUTS / "source_mesh.npz").stat().st_size, "sha256": sha256(INPUTS / "source_mesh.npz")},
            "result_json": {"path": "source_mesh_result.json", "bytes": (INPUTS / "source_mesh_result.json").stat().st_size, "sha256": sha256(INPUTS / "source_mesh_result.json")},
            "executed_alias_msh": {"path": "mesh.msh", "bytes": (INPUTS / "mesh.msh").stat().st_size, "sha256": raw_alias_hash},
            "alias_receipt": {"path": "mesh_alias_receipt.json", "bytes": (INPUTS / "mesh_alias_receipt.json").stat().st_size, "sha256": alias_receipt_hash},
        },
        "mesh_counts": {
            "source_faces": int(result["counts"]["source_faces"]),
            "source_volumes": int(result["counts"]["source_volumes"]),
            "partition_volumes": int(result["counts"]["partition_volumes"]),
            "mesh_nodes_reported_by_vf": int(result["counts"]["nodes"]),
            "msh_source_node_rows": len(source_mesh["nodes"]),
            "T4": int(result["counts"]["tetrahedra"]),
            "interface_Tri3": int(result["counts"]["interface_triangles"]),
            "exterior_Tri3": int(result["counts"]["exterior_triangles"]),
        },
        "node_inventory": {
            "mesh_npz_total_nodes": mesh_audit["npz_audit"]["npz_nodes"],
            "mesh_npz_T4_referenced_nodes": mesh_audit["npz_audit"]["npz_tet_referenced_nodes"],
            "mesh_npz_unused_node_count": len(mesh_audit["npz_audit"]["npz_unused_original_node_tags"]),
            "mesh_npz_unused_original_node_tags": mesh_audit["npz_audit"]["npz_unused_original_node_tags"],
            "mesh_msh_node_rows": len(source_mesh["nodes"]),
            "mesh_msh_tag_set_equals_npz_T4_referenced_tags": mesh_audit["npz_audit"]["msh_node_tags_equal_npz_tet_referenced_tags"],
            "unused_npz_node_tags_present_in_original_msh": mesh_audit["npz_audit"]["npz_unused_tags_present_in_msh"],
            "native_mesh_source": "preserved mesh.msh parsed by the frozen wrapper; no nodes are inserted, removed, renumbered, or filtered by this preflight package",
            "native_retention_and_mass_status": "TO_BE_REPORTED_FROM_ACTUAL_AKANTU_INITIAL_STATE",
        },
        "source_face_ancestry": {
            "status": result["source_face_ancestry"]["status"],
            "identity_scope": result["source_face_ancestry"]["identity_scope"],
            "propagation_method": result["source_face_ancestry"]["propagation_method"],
            "spatial_or_nearest_face_matching_used": result["source_face_ancestry"]["spatial_or_nearest_face_matching_used"],
            "exterior_triangle_counts": result["source_face_ancestry"]["exterior_triangle_counts"],
            "cut_surface_mapping_status": result["source_face_ancestry"]["cut_witness_ancestry"]["status"],
            "interface_has_zero_source_face_ancestry": result["source_face_ancestry"]["cut_witness_ancestry"]["interface_has_zero_source_face_ancestry"],
            "every_exterior_surface_has_exactly_one_source_face": result["source_face_ancestry"]["cut_witness_ancestry"]["every_exterior_surface_has_exactly_one_source_face"],
            "exact_msh_fid_role_arrays": {
                "path": "known_exterior_fid_roles.npz", "sha256": fid_roles_hash,
                "rows_path": "known_exterior_fid_role_rows.json", "rows_sha256": fid_rows_hash,
                "rows": len(mesh_audit["exterior_rows"]), "known_source_face_ids": fid_role_records,
                "binding_method": "MSH2 physical/geometrical tags and exact sorted original nodeTags sourceTri keys; no nearest/order guess",
            },
        },
        "mesh_geometry_diagnostics": {
            "source_cad_volume_m3": result["cad"]["source_volume_m3"],
            "tetrahedral_mesh_volume_m3": result["mesh"]["tetrahedral_volume_m3"],
            "cad_mesh_volume_relative_error": result["mesh"]["cad_mesh_volume_relative_error"],
            "interface_area_m2": result["mesh"]["triangular_interface_area_m2"],
            "occ_mesh_interface_area_relative_error": result["mesh"]["occ_mesh_interface_area_relative_error"],
            "min_tetrahedron_altitude_m": result["mesh"]["audit"]["quality"]["minimum_altitude_m"]["min"],
            "foot_fillet_local_area_errors": "approximately -9.7%; retained as a local limitation and not treated as whole-feature qualification",
            "full_feature_geometry_qualification": "NOT_ESTABLISHED",
        },
        "measured_profile_or_material_qualification": "NOT_ESTABLISHED",
        "forge_asset_manifest": "NOT_AVAILABLE_AND_NOT_REQUIRED; no substitute manifest was copied",
    }
    write_json(INPUTS / "source_geometry_provenance.json", provenance)

    old_case = json.loads((PREVIOUS_INPUTS / "case.json").read_text(encoding="utf-8"))
    if old_case.get("schema") != "pgfrag.native_cohesive_vessel_case.v1":
        raise ValueError("unexpected reference-vessel case schema")
    case = json.loads(json.dumps(old_case))
    case["classification"] = "COST_INITIALIZATION_ONLY_UNVALIDATED"
    case["source_geometry"] = {
        "status": "GEOMETRY_DIAGNOSTIC_ONLY",
        "candidate": provenance["candidate"],
        "vessel_reference": provenance["vessel_reference"],
        "source_step": provenance["canonical_step"],
        "upstream_original_v6_step_sha256": provenance["upstream_original_v6_step_sha256"],
        "source_face_lineage_sha256": provenance["source_face_lineage_sha256"],
        "witness_manifest_sha256": provenance["witness_manifest_sha256"],
        "source_step_declared_unit": result["source_step_declared_length_unit"],
        "source_step_unit_conversion": result["step_unit_conversion"],
        "mesh_msh": provenance["mesh_sources"]["mesh_msh"],
        "mesh_npz": provenance["mesh_sources"]["mesh_npz"],
        "mesh_result": provenance["mesh_sources"]["result_json"],
        "mesh_counts": provenance["mesh_counts"],
        "node_inventory": provenance["node_inventory"],
        "source_face_ancestry": provenance["source_face_ancestry"],
        "mesh_geometry_diagnostics": provenance["mesh_geometry_diagnostics"],
        "measured_profile_or_material_qualification": "NOT_ESTABLISHED",
    }
    case["bulk_material"]["calibration"] = "UNVALIDATED_COST_INITIALIZATION_ONLY_PLACEHOLDER; NO SPECIMEN DATA"
    case["interface"]["calibration"] = "UNVALIDATED_COST_INITIALIZATION_ONLY_PLACEHOLDER; NO SPECIMEN DATA"
    case["interface"]["weakplane_label"] = "PREINSERTED_WEAK_SEAM_UNCALIBRATED"
    case["limits"].update({
        "max_steps": 100000,
        "wall_time_s": 300,
        "output_cap_bytes": HISTORY_CAP,
    })
    case["output"].update({"node_stride": 50, "flush_stride": 250, "all_native_qp_steps": True})
    case["case_name"] = "reference_v6_zero_load_native_initialization_preflight"
    case["limitations"] = [
        "v6 named-STEP reference candidate geometry only; wall profile and whole-feature shape are not measured or physically qualified",
        "verified OCC witness ancestry is conditional within declared CAD tolerances and does not establish exact metrology; local foot-fillet mesh area errors are approximately -9.7%",
        "bulk and interface material values remain UNVALIDATED cost-initialization placeholders with no specimen physical data",
        "one native initialization-only preflight; no solver evolution, fracture, fragment result, training eligibility, or physical qualification",
        "native initial positive-stiffness bound is not a full-trajectory stability guarantee",
        "known exterior FID/role rows bind to exact original MSH sourceTri nodeTags and geometrical entity tags; no nearest or coordinate matching",
    ]

    input_names = [
        *MODULES, "material.dat", "mesh.msh", "mesh.original.msh", "mesh_alias_receipt.json",
        "source_mesh_result.json", "source_mesh.npz", "known_exterior_fid_roles.npz", "t4_cell_row_correspondence.npz",
        "known_exterior_fid_role_rows.json", "source_geometry_provenance.json",
    ]
    case["input_hashes"] = {name: sha256(INPUTS / name) for name in sorted(input_names)}
    case["mesh_alias_receipt"] = {
        "path": "mesh_alias_receipt.json", "sha256": alias_receipt_hash,
        "source_mesh_sha256": sha256(INPUTS / "mesh.original.msh"),
        "executed_mesh_sha256": raw_alias_hash,
        "operation": "REPLACE_TWO_PHYSICAL_NAME_STRINGS_ONLY",
        "known_exterior_fid_role_mapping_sha256": fid_roles_hash,
    }
    write_json(INPUTS / "case.json", case)

    # Record exact MSH parser equivalence before any native execution.
    parser_equivalence = {
        "schema": "pgfrag.v6_msh_alias_parser_equivalence.v1",
        "status": "PASS_RAW_AND_PARSER_EQUIVALENCE",
        "parser": "frozen inputs/cohesive_mesh.py read_gmsh22",
        "source_msh_sha256": sha256(INPUTS / "mesh.original.msh"),
        "executed_alias_msh_sha256": raw_alias_hash,
        "raw_nodes_block_byte_equal": mesh_audit["alias_receipt"]["nodes_block_byte_equal"],
        "raw_elements_block_byte_equal": mesh_audit["alias_receipt"]["elements_block_byte_equal"],
        "all_bytes_outside_physical_names_byte_equal": mesh_audit["alias_receipt"]["all_bytes_outside_physical_names_byte_equal"],
        "raw_element_records_ids_tags_connectivity_and_order_equal": mesh_audit["alias_receipt"]["raw_element_records_and_order_equal"],
        "parsed_node_coordinates_and_order_equal": mesh_audit["alias_receipt"]["parsed_node_coordinates_and_order_equal"],
        "parsed_node_tags_order_equal": mesh_audit["alias_receipt"]["parsed_node_tag_order_equal"],
        "parsed_connectivity_element_ids_physical_tags_and_order_equal": mesh_audit["alias_receipt"]["parsed_connectivity_element_ids_physical_tags_and_order_equal"],
        "run05_t4_cross_format_match": "EXACT_ELEMENT_ID_MAP_TO_ORIGINAL_CELL_TAGS_AND_ORDERED_CORNER_NODE_TAGS",
        "run05_t4_msh_row_to_npz_row_correspondence_path": "t4_cell_row_correspondence.npz",
        "run05_t4_msh_row_to_npz_row_correspondence_sha256": sha256(INPUTS / "t4_cell_row_correspondence.npz"),
        "volume_physical_group_unchanged": mesh_audit["alias_receipt"]["volume_physical_group_unchanged"],
        "physical_names_original": mesh_audit["alias_receipt"]["physical_names_original"],
        "physical_names_derived": mesh_audit["alias_receipt"]["physical_names_derived"],
        "counts": {"msh_nodes": len(source_mesh["nodes"]), "T4": len(source_mesh["tetrahedra"]),
                   "interface_Tri3": len([row for row in source_mesh["triangles"] if row["physical_tag"] == 2]),
                   "exterior_Tri3": len(mesh_audit["exterior_rows"])},
        "known_fid_role_binding_rows": len(mesh_audit["exterior_rows"]),
        "known_fid_role_binding_sha256": fid_roles_hash,
        "msh_t4_to_npz_cell_row_correspondence_sha256": sha256(INPUTS / "t4_cell_row_correspondence.npz"),
        "method": "raw MSH block byte comparison plus frozen parser connectivity/tag/order comparison; FID/role bound by exact sorted source nodeTags keys",
    }
    write_json(DATA / "parser_equivalence.json", parser_equivalence)

    source_density = float(case["bulk_material"]["rho_kg_m3"])
    lineage_forecast, lineage_detail = synthetic_lineage_forecast(mesh_audit["alias"], source_density)
    input_bytes = tree_bytes(EXPERIMENT) + tree_bytes(DATA)
    projected_total = input_bytes + lineage_forecast + NON_LINEAGE_OUTPUT_ALLOWANCE + FORECAST_REPORT_ALLOWANCE
    budget = {
        "schema": "pgfrag.reference_v6_native_initialization_preflight_budget.v1",
        "status": "WHOLE_OUTPUT_FORECAST_PASS_RUN_NOT_STARTED",
        "native_solver_execution_started": False,
        "solver_steps_forecast_as_executed": 0,
        "solver_history_arrays_allocated": False,
        "classification": "COST_INITIALIZATION_ONLY_UNVALIDATED; PREINSERTED_WEAK_SEAM_UNCALIBRATED",
        "input_freeze": {
            "case_sha256": sha256(INPUTS / "case.json"),
            "case_bytes": (INPUTS / "case.json").stat().st_size,
            "frozen_input_hashes": case["input_hashes"],
            "current_module_hashes_match_frozen_copies": all(sha256(current_modules[name]) == sha256(INPUTS / name) for name in MODULES),
            "source_mesh_sha256": sha256(INPUTS / "mesh.original.msh"),
            "alias_mesh_sha256": raw_alias_hash,
            "alias_receipt_sha256": alias_receipt_hash,
            "known_fid_role_array_sha256": fid_roles_hash,
            "t4_cell_row_correspondence_sha256": sha256(INPUTS / "t4_cell_row_correspondence.npz"),
            "runner_script_hashes": {name: sha256(EXPERIMENT / name) for name in (
                "prepare_preflight.py", "run_preflight.ps1", "run_native_preflight.sh")},
            "forge_asset_manifest_copied": False,
        },
        "native_preflight_output_forecast": {
            "synthetic_lineage": lineage_detail,
            "lineage_forecast_bytes": lineage_forecast,
            "history_cap_bytes": HISTORY_CAP,
            "history_cap_scope": "history arrays only; preflight_only returns before history allocation",
            "native_history_bytes_forecast": 0,
            "non_lineage_native_output_allowance_bytes": NON_LINEAGE_OUTPUT_ALLOWANCE,
            "forecast_report_self_allowance_bytes": FORECAST_REPORT_ALLOWANCE,
            "all_prepared_input_and_script_bytes": input_bytes,
            "whole_inputs_and_outputs_forecast_bytes": projected_total,
            "whole_inputs_and_outputs_cap_bytes": ALL_OUTPUT_CAP,
            "remaining_budget_bytes": ALL_OUTPUT_CAP - projected_total,
            "pass": projected_total <= ALL_OUTPUT_CAP,
            "no_runtime_or_solver_cost_claim": True,
        },
        "actual_v6_mesh_counts": {
            "MSH_nodes": len(source_mesh["nodes"]), "T4": len(source_mesh["tetrahedra"]),
            "interface_Tri3": len([row for row in source_mesh["triangles"] if row["physical_tag"] == 2]),
            "exterior_Tri3": len(mesh_audit["exterior_rows"]),
            "source_face_ids": len(mesh_audit["fid_ids"]),
            "partition_volumes": int(result["counts"]["partition_volumes"]),
        },
        "one_attempt_policy": "Exactly one native Akantu initialization-only invocation; preserve any failure and do not retry, retune, or change parameters.",
    }
    if projected_total > ALL_OUTPUT_CAP:
        budget["status"] = "WHOLE_OUTPUT_FORECAST_FAIL_RUN_NOT_STARTED"
        write_json(DATA / "preflight_budget.json", budget)
        raise ValueError("whole inputs+outputs forecast exceeds the declared 256 MiB cap; no native run permitted")
    write_json(DATA / "preflight_budget.json", budget)

    # The command is duplicated in the shell runner, run_control and the
    # execution record so the exact runtime and resource limits are reviewable.
    wsl_root = "/mnt/c/Users/xiemo/Desktop/ShardForge"
    wsl_experiment = f"{wsl_root}/pgfrag/experiments/continuation_20261016/reference_v6_preflight"
    wsl_runtime = f"{wsl_root}/third_party/akantu-runtime/venv/bin/python"
    argv = [
        wsl_experiment + "/inputs/cohesive_vessel.py",
        "--case", wsl_experiment + "/inputs/case.json",
        "--out", wsl_experiment + "/native_preflight/run01",
        "--preflight-only",
    ]
    run_control = {
        "schema": "pgfrag.reference_v6_native_initialization_preflight_control.v1",
        "status": "READY_FOR_SINGLE_NATIVE_INITIALIZATION_PREFLIGHT",
        "run_mode": "NATIVE_INITIALIZATION_ONLY_ZERO_SOLVER_STEPS",
        "runtime": "WSL Ubuntu; Akantu 5.0.7-post1; Python 3.12.15; NumPy 2.5.3",
        "runtime_python": wsl_runtime,
        "command_argv": argv,
        "command_shell": "bash run_native_preflight.sh; timeout --foreground --signal=TERM --kill-after=5s 300s; taskset -c 0; external PowerShell watchdog 330s",
        "external_watchdog_s": 330,
        "wrapper_wall_time_s": int(case["limits"]["wall_time_s"]),
        "output_budget_bytes": ALL_OUTPUT_CAP,
        "history_cap_bytes": HISTORY_CAP,
        "thread_environment": {
            "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE",
        },
        "cpu_affinity": "CPU 0 via taskset -c 0",
        "case_sha256": sha256(INPUTS / "case.json"),
        "input_hashes": case["input_hashes"],
        "source_geometry_hashes": {
            "canonical_step_sha256": provenance["canonical_step"]["sha256"],
            "upstream_original_v6_step_sha256": provenance["upstream_original_v6_step_sha256"],
            "source_face_lineage_sha256": provenance["source_face_lineage_sha256"],
            "run05_result_sha256": sha256(INPUTS / "source_mesh_result.json"),
            "original_msh_sha256": sha256(INPUTS / "mesh.original.msh"),
            "alias_msh_sha256": raw_alias_hash,
        },
        "runner_script_hashes": {name: sha256(EXPERIMENT / name) for name in (
            "prepare_preflight.py", "run_preflight.ps1", "run_native_preflight.sh")},
        "stdout_path": "native_preflight/stdout.log",
        "stderr_path": "native_preflight/stderr.log",
        "exit_code_path": "native_preflight/native_exit_code.txt",
        "forecast_status": budget["status"],
        "forecasted_whole_inputs_and_outputs_bytes": projected_total,
    }
    write_json(DATA / "run_control.json", run_control)

    readme = f"""# Reference vessel v6 native initialization preflight

Prepared for one native Akantu initialization-only invocation from the verified v6 tetra mesh. It uses the current `pgfrag/src/pgfrag` wrapper, mesh, and stability source frozen byte-for-byte under `inputs/`. The run command uses `--preflight-only`; it creates no solver steps or history arrays.

The MSH copy used by the wrapper changes only the two physical names needed by the wrapper (`cut_interface` → `interface`, `source_exterior` → `exterior`). `mesh.original.msh` is preserved byte-for-byte. [`parser_equivalence.json`](parser_equivalence.json) records raw Nodes/Elements block equality, raw tag/order equality, parsed connectivity equality, and the unchanged volume physical group. Exterior FID and role arrays are bound to each original MSH exterior triangle by its exact sorted original nodeTags key plus the MSH geometrical entity tag.

The source is a Daily Aberdeen 108x53, nominal 2 mm wall reference candidate, not a measured or physically qualified source profile. The run05 result reports 86,624 T4, 356 interface Tri3, and 44,182 exterior Tri3; OCC witness ancestry is conditional within declared CAD tolerances. The meshed CAD-volume error is {result['mesh']['cad_mesh_volume_relative_error']:.9%}; local foot-fillet area errors remain approximately −9.7% and do not qualify the whole feature. Material/interface values remain unvalidated cost-initialization placeholders, and the interface remains a preinserted weak seam.

The prepared forecast reserves {projected_total:,} bytes of the {ALL_OUTPUT_CAP:,}-byte whole-input/output budget. The full-history cap is {HISTORY_CAP:,} bytes; `preflight_only` returns before history allocation. No runtime duration or solver-evolution cost is implied. Canonical STEP, upstream source STEP, source face lineage, run05 result, original MSH, alias MSH, and FID/role-array hashes are frozen in `inputs/case.json`, `inputs/source_geometry_provenance.json`, and `run_control.json`.

The single attempt is run by `run_preflight.ps1` and audited by [`audit_preflight.py`](../../data/continuation_20261016/reference_v6_preflight/audit_preflight.py). The attempt is diagnostic initialization only; it does not establish real-world geometry review, material validity, fracture, training eligibility, or solver qualification.
"""
    (DATA / "README.md").write_text(readme, encoding="utf-8")
    # Update the forecast after adding the README and run-control artifacts.
    input_bytes = tree_bytes(EXPERIMENT) + tree_bytes(DATA)
    projected_total = input_bytes + lineage_forecast + NON_LINEAGE_OUTPUT_ALLOWANCE + FORECAST_REPORT_ALLOWANCE
    budget["native_preflight_output_forecast"]["all_prepared_input_and_script_bytes"] = input_bytes
    budget["native_preflight_output_forecast"]["whole_inputs_and_outputs_forecast_bytes"] = projected_total
    budget["native_preflight_output_forecast"]["remaining_budget_bytes"] = ALL_OUTPUT_CAP - projected_total
    budget["native_preflight_output_forecast"]["pass"] = projected_total <= ALL_OUTPUT_CAP
    if projected_total > ALL_OUTPUT_CAP:
        budget["status"] = "WHOLE_OUTPUT_FORECAST_FAIL_RUN_NOT_STARTED"
        run_control["status"] = "BLOCKED_WHOLE_OUTPUT_FORECAST_EXCEEDED"
        write_json(DATA / "preflight_budget.json", budget)
        write_json(DATA / "run_control.json", run_control)
        raise ValueError("updated whole inputs+outputs forecast exceeds 256 MiB; no native run permitted")
    write_json(DATA / "preflight_budget.json", budget)
    run_control["forecasted_whole_inputs_and_outputs_bytes"] = projected_total
    write_json(DATA / "run_control.json", run_control)
    print(json.dumps({
        "status": "READY_FOR_SINGLE_NATIVE_INITIALIZATION_PREFLIGHT",
        "counts": budget["actual_v6_mesh_counts"],
        "lineage_forecast_bytes": lineage_forecast,
        "whole_output_forecast_bytes": projected_total,
        "whole_output_cap_bytes": ALL_OUTPUT_CAP,
        "case_sha256": sha256(INPUTS / "case.json"),
        "input_directory": str(INPUTS),
    }, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
