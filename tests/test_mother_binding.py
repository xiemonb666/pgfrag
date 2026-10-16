from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pytest

from pgfrag.mother_binding import (
    MotherBindingError,
    build_mother_binding,
    write_mother_binding,
)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _positive_tet(nodes: np.ndarray, values: list[int]) -> list[int]:
    points = nodes[np.asarray(values, dtype=np.int64)]
    signed = float(np.dot(points[1] - points[0],
                          np.cross(points[2] - points[0], points[3] - points[0])))
    assert signed != 0.0
    return values if signed > 0.0 else [values[1], values[0], values[2], values[3]]


def _fixture(tmp_path: Path, *, reorder_storage: bool = False,
             exterior_fid: str | None = None, exterior_role: str = "foot",
             identity_role: str = "foot", interface_parent_rows: list[int] | None = None,
             provenance_extra: dict[str, object] | None = None) -> dict[str, Path | str]:
    root = tmp_path / "candidate-name-is-not-identity"
    geometry = root / "geometry"
    geometry.mkdir(parents=True)
    nodes = np.asarray([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0], [0.0, 0.0, -1.0], [2.0, 2.0, 2.0],
    ], dtype=np.float64)  # Last node is deliberately unused by the MSH/T4s.
    tets = np.asarray([
        _positive_tet(nodes, [0, 1, 2, 3]),
        _positive_tet(nodes, [0, 2, 1, 4]),
    ], dtype=np.int64)
    node_tags = np.asarray([91, 17, 804, 31, 712, 9901], dtype=np.int64)
    cell_tags = np.asarray([7009, 401], dtype=np.int64)
    partition = np.asarray([0, 1], dtype=np.int32)
    if reorder_storage:
        node_order = np.asarray([5, 2, 4, 0, 3, 1], dtype=np.int64)
        old_to_new = np.empty(len(node_order), dtype=np.int64)
        old_to_new[node_order] = np.arange(len(node_order))
        nodes = nodes[node_order]
        node_tags = node_tags[node_order]
        tets = old_to_new[tets]
        cell_order = np.asarray([1, 0], dtype=np.int64)
        tets, cell_tags, partition = tets[cell_order], cell_tags[cell_order], partition[cell_order]

    incidence: dict[tuple[int, int, int], list[int]] = {}
    local_faces = ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))
    for row, tet in enumerate(tets):
        for corners in local_faces:
            key = tuple(sorted(int(tet[index]) for index in corners))
            incidence.setdefault(key, []).append(row)
    exterior = np.asarray([key for key, parents in incidence.items() if len(parents) == 1],
                          dtype=np.int32).reshape(-1, 3)
    interface = np.asarray([key for key, parents in incidence.items() if len(parents) == 2],
                           dtype=np.int32).reshape(-1, 3)
    assert exterior.shape == (6, 3) and interface.shape == (1, 3)
    interface_key = tuple(sorted(map(int, interface[0])))
    parents = incidence[interface_key]
    parent_rows = interface_parent_rows or parents
    interface_parents = np.asarray([parent_rows], dtype=np.int32)
    interface_parent_partition = partition[interface_parents]
    if exterior_fid is None:
        exterior_fid = "FID-FOOT-SURFACE"
    exterior_fids = np.asarray([exterior_fid] * len(exterior), dtype="U64")
    exterior_roles = np.asarray([exterior_role] * len(exterior), dtype="U32")
    arrays = {
        "nodes_m": nodes,
        "tets": tets.astype(np.int32),
        "original_node_tags": node_tags,
        "original_cell_tags": cell_tags,
        "source_partition_volume": partition,
        "volume_m3": np.full(len(tets), 1.0 / 6.0, dtype=np.float64),
        "interface_triangles": interface,
        "interface_parents": interface_parents,
        "interface_parent_partition": interface_parent_partition.astype(np.int8),
        "exterior_triangles": exterior,
        "exterior_surface_tags": np.full(len(exterior), 17, dtype=np.int32),
        "exterior_source_face_ids": exterior_fids,
        "exterior_roles": exterior_roles,
        "exterior_ancestry_exact": np.ones(len(exterior), dtype=np.bool_),
        "partition_entity_tags": np.asarray([3, 4], dtype=np.int32),
    }

    # The asset source face index is 5 while its role is foot. The binding must
    # take the role from the witnessed FID table, never from a numeric face code.
    source_step = geometry / "solid.step"
    source_step.write_bytes(b"source step bytes\n")
    section = {"schema": "vf.meridian_section.v1",
               "coordinate_frame": "radius_mm_vs_z_mm; revolved about Z axis; Z-up",
               "segments": [{"role": "foot", "curve": "line",
                             "points_rz_mm": [[0.0, 0.0], [1.0, 0.0]]}]}
    _write_json(geometry / "section.json", section)
    lineage = {"face_count": 1, "faces": {"5": {
        "role": "foot", "source_segment_indices": [0],
        "lineage_basis": "test source ancestry"}}}
    _write_json(geometry / "face_lineage.json", lineage)
    surface_path = geometry / "surface.npz"
    np.savez_compressed(surface_path,
                        vertices_mm=np.asarray([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]]),
                        faces=np.asarray([[0, 1, 2]], dtype=np.int32))
    report = {
        "candidate_geometry_status": "PASS_STRUCTURAL_CAD_ONLY",
        "candidate_status": "SIZE_REFERENCE_GEOMETRY_CANDIDATE_ONLY",
        "checks": {"mesh_is_watertight": True},
        "measurements": {"artifact_sha256": {
            "solid.step": _sha(source_step),
            "section.json": _sha(geometry / "section.json"),
            "face_lineage.json": _sha(geometry / "face_lineage.json"),
            "surface.npz": _sha(surface_path),
        }},
    }
    _write_json(root / "geometry_report.json", report)
    provenance: dict[str, object] = {
        "schema": "vf.revolved_profile_candidate_provenance.v1",
        "candidate_status": "SIZE_REFERENCE_GEOMETRY_CANDIDATE_ONLY",
        "source_profile_fidelity": "NOT_ESTABLISHED; test candidate profile",
        "physical_material_properties": "UNKNOWN; no material data",
        "not_claimed": ["not exact source profile", "not calibrated"],
        "record": {"object_number": "TEST-MOTHER-001"},
    }
    if provenance_extra:
        provenance.update(provenance_extra)
    _write_json(root / "provenance.json", provenance)

    canonical_step = tmp_path / "canonical.step"
    canonical_step.write_bytes(b"canonical step bytes\n")
    witness_dir = tmp_path / "witnesses"
    witness_dir.mkdir()
    witness_path = witness_dir / "FID-FOOT-SURFACE.step"
    witness_path.write_bytes(b"one face witness\n")
    base_hashes = {
        "solid.step": _sha(source_step),
        "section.json": _sha(geometry / "section.json"),
        "face_lineage.json": _sha(geometry / "face_lineage.json"),
        "geometry_report.json": _sha(root / "geometry_report.json"),
    }
    identity = {
        "schema": "vesselforge.daily_native_face_identity.v1",
        "source_hashes": base_hashes,
        "faces": [{
            "native_face_index": 2,
            "source_original_brep_face_id": "FID-FOOT-SURFACE",
            "source_original_brep_face_index": 5,
            "role": identity_role,
            "source_segment_indices": [0],
        }],
    }
    identity_path = tmp_path / "native_face_identity.json"
    _write_json(identity_path, identity)
    witness = {
        "schema": "vesselforge.native_face_witness_manifest.v1",
        "mapping_status": "EXACT_NATIVE_FACE_WITNESSES",
        "source_step_sha256": base_hashes["solid.step"],
        "source_face_lineage_sha256": base_hashes["face_lineage.json"],
        "native_body_step_sha256": _sha(canonical_step),
        "native_solid_face_count": 1,
        "witness_count": 1,
        "identity_scope": "OCC witness ancestry within declared tolerances",
        "witnesses": [{
            "source_face_id": "FID-FOOT-SURFACE", "role": "foot",
            "native_face_child_verified": True, "witness_face_entity_count": 1,
            "witness_sha256": _sha(witness_path), "witness_bytes": witness_path.stat().st_size,
            "witness_path": str(witness_path),
        }],
    }
    witness_manifest_path = tmp_path / "witness_manifest.json"
    _write_json(witness_manifest_path, witness)
    v15 = {
        "schema": "vesselforge.daily_native_face_witness_probe.v1",
        "status": "PASS_CAD_PARENT_MAP_WITHIN_DECLARED_OCC_TOLERANCES",
        "source_hashes": base_hashes,
        "native_body_step_sha256": _sha(canonical_step),
        "native_face_identity_sha256": _sha(identity_path),
        "native_face_witness_manifest_sha256": _sha(witness_manifest_path),
    }
    v15_path = tmp_path / "v15_receipt.json"
    _write_json(v15_path, v15)
    label = {
        "schema": "vesselforge.daily_named_step_label_probe.v1",
        "source_hashes": base_hashes,
        "native_body_step_sha256": _sha(canonical_step),
    }
    label_path = tmp_path / "label_import_receipt.json"
    _write_json(label_path, label)

    # Gmsh serialization orders are independent of NPZ row order. Original
    # tags and exact T4 corner sequences carry identity across both.
    msh_path = tmp_path / "partition" / "mesh.msh"
    msh_path.parent.mkdir()
    used_rows = np.unique(tets.reshape(-1))
    msh_node_rows = list(reversed(used_rows.tolist()))
    msh_node_tags = {int(node_tags[row]): row for row in msh_node_rows}
    msh_cell_order = list(reversed(range(len(tets))))
    msh_lines = [
        "$MeshFormat", "2.2 0 8", "$EndMeshFormat",
        "$PhysicalNames", "3", '2 2 "cut_interface"',
        '2 3 "source_exterior"', '3 1 "partitioned_body"', "$EndPhysicalNames",
        "$Nodes", str(len(msh_node_rows)),
    ]
    for row in msh_node_rows:
        tag = int(node_tags[row])
        xyz = " ".join(format(float(value), ".17g") for value in nodes[row])
        msh_lines.append(f"{tag} {xyz}")
    msh_lines.extend(["$EndNodes", "$Elements", str(len(exterior) + len(interface) + len(tets))])
    element_tag = 1
    for face in reversed(exterior.tolist()):
        tags = [int(node_tags[row]) for row in face]
        msh_lines.append(f"{element_tag} 2 2 3 17 " + " ".join(map(str, tags)))
        element_tag += 1
    for face in interface.tolist():
        tags = [int(node_tags[row]) for row in face]
        msh_lines.append(f"{element_tag} 2 2 2 18 " + " ".join(map(str, tags)))
        element_tag += 1
    for row in msh_cell_order:
        tags = [int(node_tags[int(node_row)]) for node_row in tets[row]]
        entity = int(partition[row]) + 1
        msh_lines.append(f"{int(cell_tags[row])} 4 2 1 {entity} " + " ".join(map(str, tags)))
    msh_lines.extend(["$EndElements", ""])
    msh_path.write_text("\n".join(msh_lines), encoding="ascii")
    npz_path = msh_path.with_suffix(".npz")
    np.savez_compressed(npz_path, **arrays)

    partition_result_path = msh_path.parent / "result.json"
    canonical_key = "canonical_step_input.step"
    result = {
        "schema": "vesselforge.generic_step_tetra_partition.v1",
        "status": "GEOMETRY_DIAGNOSTIC_ONLY",
        "units": "SI",
        "source_step": str(canonical_step),
        "source_step_sha256": _sha(canonical_step),
        "qualification": {
            "real_world_geometry_review": "NOT_ESTABLISHED_BY_MESHING",
            "physical_material_validity": "NOT_ASSESSED",
            "fracture_or_fragment_result": "NOT_CREATED",
            "solver_qualification": "NOT_ESTABLISHED",
        },
        "counts": {"nodes": len(nodes), "tetrahedra": len(tets),
                   "exterior_triangles": len(exterior), "interface_triangles": len(interface)},
        "source_face_ancestry": {
            "status": "VERIFIED_OCC_WITNESS_ANCESTRY",
            "manifest_sha256": _sha(witness_manifest_path),
            "canonical_step_sha256": _sha(canonical_step),
            "upstream_source_step_sha256": _sha(source_step),
            "source_face_lineage_sha256": _sha(geometry / "face_lineage.json"),
            "identity_scope": "OCC witness ancestry within declared kernel tolerances",
        },
        "files": {
            name: {"sha256": _sha(path), "bytes": path.stat().st_size}
            for name, path in (("mesh.msh", msh_path), ("mesh.npz", npz_path))
        },
    }
    _write_json(partition_result_path, result)
    audit_input_paths = {
        "solid.step": source_step,
        "section.json": geometry / "section.json",
        "face_lineage.json": geometry / "face_lineage.json",
        "geometry_report.json": root / "geometry_report.json",
        "native_face_identity.json": identity_path,
        "witness_manifest.json": witness_manifest_path,
        "v15_receipt.json": v15_path,
        "label_import_receipt.json": label_path,
        canonical_key: canonical_step,
    }
    input_hashes = {name: _sha(path) for name, path in audit_input_paths.items()}
    audit = {
        "schema": "vesselforge.reference_v6_witness_tetra_run.v1",
        "status": "PASS_WITNESS_ANCESTRY_TETRA_MESH",
        "result_sha256": _sha(partition_result_path),
        "result_json_bytes": partition_result_path.stat().st_size,
        "input_hashes": input_hashes,
        "source_geometry_and_lineage": {"asset_hashes": dict(input_hashes)},
        "independent_mesh_audit": {
            "status": "PASS",
            "all_exterior_triangles_have_known_manifest_role": True,
            "every_interface_triangle_has_two_distinct_partition_parents": True,
            "interface_source_face_id_count": 0,
            "independent_audit": {"checks": {
                "positive_oriented_non_degenerate_tetrahedra": True,
                "interface_matches_all_cross_partition_faces": True,
                "two_parent_face_windings_opposed": True,
            }},
        },
    }
    audit_path = msh_path.parent / "audit_receipt.json"
    _write_json(audit_path, audit)
    return {
        "source_asset_dir": root,
        "partition_result_path": partition_result_path,
        "audit_receipt_path": audit_path,
        "trusted_audit_receipt_sha256": _sha(audit_path),
        "canonical_step_path": canonical_step,
        "canonical_step_input_key": canonical_key,
        "native_face_identity_path": identity_path,
        "witness_manifest_path": witness_manifest_path,
        "v15_receipt_path": v15_path,
        "label_import_receipt_path": label_path,
        "witness_directory": witness_dir,
        "msh_path": msh_path,
        "npz_path": npz_path,
        "partition_result": partition_result_path,
        "source_step": source_step,
        "arrays": arrays,
    }


def _build(fixture: dict[str, Path | str]):
    args = {key: fixture[key] for key in (
        "source_asset_dir", "partition_result_path", "audit_receipt_path",
        "trusted_audit_receipt_sha256", "canonical_step_path", "canonical_step_input_key",
        "native_face_identity_path", "witness_manifest_path", "v15_receipt_path",
        "label_import_receipt_path", "witness_directory")}
    return build_mother_binding(**args)  # type: ignore[arg-type]


def test_binding_joins_reordered_npz_and_msh_by_tags_and_preserves_unknown_status(tmp_path):
    fixture = _fixture(tmp_path, reorder_storage=True)
    binding = _build(fixture)

    assert binding.binding_id == binding.summary["binding_id"]
    assert binding.mother_id == "TEST-MOTHER-001"
    assert binding.source_face_roles == ("foot",)
    assert binding.summary["qualification"]["structural_geometry_status"] == "PASS_STRUCTURAL_CAD_ONLY"
    assert binding.summary["qualification"]["source_profile_fidelity"].startswith("NOT_ESTABLISHED")
    assert binding.summary["qualification"]["physical_material_status"] == "UNKNOWN"
    assert binding.summary["qualification"]["fracture_status"] == "NOT_CREATED"
    assert binding.summary["qualification"]["observation_status"] == "NOT_CREATED_BY_GEOMETRY_BINDING"
    assert binding.summary["qualification"]["pretrain_eligible"] is False
    assert binding.summary["counts"]["npz_node_count"] == 6
    assert binding.summary["counts"]["msh_node_count"] == 5
    assert binding.summary["counts"]["unused_node_count"] == 1
    assert binding.interface_parent_cell_tags.shape == (1, 2)
    assert len(set(binding.interface_parent_cell_tags[0].tolist())) == 2
    assert binding.cell_corner_node_tags.shape == (2, 4)
    assert binding.exterior_source_face_index.shape == (6,)
    assert not binding.node_tags.flags.writeable
    assert not binding.exterior_triangle_node_tags.flags.writeable
    with pytest.raises(TypeError):
        binding.summary["binding_id"] = "changed"  # type: ignore[index]


def test_numeric_face_index_does_not_supply_semantic_role(tmp_path):
    fixture = _fixture(tmp_path, exterior_fid="FID-FOOT-SURFACE", exterior_role="base")
    with pytest.raises(MotherBindingError, match="unknown or conflicts"):
        _build(fixture)


def test_unknown_exterior_fid_fails_closed(tmp_path):
    fixture = _fixture(tmp_path, exterior_fid="FID-NOT-IN-SOURCE-TABLE")
    with pytest.raises(MotherBindingError, match="unknown or conflicts"):
        _build(fixture)


def test_identity_role_conflict_with_source_lineage_fails(tmp_path):
    fixture = _fixture(tmp_path, identity_role="body_outer")
    with pytest.raises(MotherBindingError, match="role/segments conflict"):
        _build(fixture)


def test_duplicate_interface_parent_rows_fail(tmp_path):
    fixture = _fixture(tmp_path, interface_parent_rows=[0, 0])
    with pytest.raises(MotherBindingError, match="parent rows do not match exact T4 incidence"):
        _build(fixture)


def test_audit_review_anchor_and_source_hash_chain_are_required(tmp_path):
    fixture = _fixture(tmp_path)
    original_anchor = str(fixture["trusted_audit_receipt_sha256"])
    audit_path = Path(fixture["audit_receipt_path"])
    audit_path.write_text(audit_path.read_text(encoding="utf-8") + " ", encoding="utf-8")
    fixture["trusted_audit_receipt_sha256"] = original_anchor
    with pytest.raises(MotherBindingError, match="caller-supplied trusted hash"):
        _build(fixture)


def test_positive_training_claim_is_rejected_even_for_geometry_only_binding(tmp_path):
    fixture = _fixture(tmp_path, provenance_extra={"pretrain_eligible": True})
    with pytest.raises(MotherBindingError, match="positive physical or eligibility claim"):
        _build(fixture)


def test_output_is_bounded_hashed_and_does_not_duplicate_source_mesh_json(tmp_path):
    fixture = _fixture(tmp_path)
    binding = _build(fixture)
    output = tmp_path / "binding-output"
    receipt = write_mother_binding(binding, output)

    arrays_path = output / "binding_arrays.npz"
    receipt_path = output / "binding_receipt.json"
    saved = json.loads(receipt_path.read_text(encoding="utf-8"))
    assert receipt["binding_id"] == saved["binding_id"]
    assert saved["files"]["binding_arrays.npz"]["sha256"] == _sha(arrays_path)
    assert arrays_path.stat().st_size + receipt_path.stat().st_size < 32 * 1024 * 1024
    assert "cell_corner_node_tags" not in saved
    with np.load(arrays_path, allow_pickle=False) as arrays:
        assert arrays["cell_corner_node_tags"].shape == (2, 4)
        assert arrays["unused_node_rows"].shape == (1,)
        assert arrays["interface_parent_cell_tags"].shape == (1, 2)


def test_tampered_partition_sidecar_fails_before_parsing(tmp_path):
    fixture = _fixture(tmp_path)
    npz_path = Path(fixture["npz_path"])
    npz_path.write_bytes(npz_path.read_bytes() + b"tamper")
    with pytest.raises(MotherBindingError, match="differs from its result receipt"):
        _build(fixture)
