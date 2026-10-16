"""Independently audit the saved native v6 zero-step initialization artifacts."""
from __future__ import annotations

import hashlib
import json
import math
import re
import shlex
import sys
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[4]
EXPERIMENT = ROOT / "pgfrag/experiments/continuation_20261016/reference_v6_preflight"
INPUTS = EXPERIMENT / "inputs"
RUN = EXPERIMENT / "native_preflight/run01"
DATA = ROOT / "pgfrag/data/continuation_20261016/reference_v6_preflight"
ALL_OUTPUT_CAP = 256 * 1024 * 1024
HISTORY_CAP = 64 * 1024 * 1024


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


def section(data: bytes, name: str) -> tuple[int, int, bytes]:
    pattern = re.compile(rb"(?m)^\$" + name.encode("ascii") + rb"\r?\n.*?^\$End" + name.encode("ascii") + rb"\r?$", re.S)
    match = pattern.search(data)
    if not match:
        raise ValueError(f"missing MSH2 ${name} block")
    return match.start(), match.end(), match.group(0)


def parse_msh22(path: Path) -> dict[str, Any]:
    """Minimal independent ASCII MSH2 parser retaining source row ordinals."""
    raw = path.read_bytes()
    _, _, format_block = section(raw, "MeshFormat")
    fmt = format_block.splitlines()[1].split()
    if len(fmt) != 3 or fmt[0] != b"2.2" or fmt[1] != b"0":
        raise ValueError("audit parser requires ASCII MSH2")
    _, _, nodes_block = section(raw, "Nodes")
    node_lines = nodes_block.splitlines()
    node_count = int(node_lines[1])
    node_tags: list[int] = []
    nodes: list[list[float]] = []
    for row in node_lines[2:-1]:
        fields = row.split()
        if len(fields) != 4:
            raise ValueError("malformed MSH2 node row")
        node_tags.append(int(fields[0]))
        nodes.append([float(fields[1]), float(fields[2]), float(fields[3])])
    if len(nodes) != node_count or len(set(node_tags)) != node_count:
        raise ValueError("MSH2 node count or node tags are invalid")
    node_row = {tag: row for row, tag in enumerate(node_tags)}

    _, _, elements_block = section(raw, "Elements")
    element_lines = elements_block.splitlines()
    element_count = int(element_lines[1])
    elements: list[dict[str, Any]] = []
    for row_ordinal, line in enumerate(element_lines[2:-1]):
        fields = [int(token) for token in line.split()]
        if len(fields) < 3:
            raise ValueError(f"malformed MSH2 element row {row_ordinal}")
        element_id, element_type, tag_count = fields[:3]
        if tag_count < 0 or len(fields) < 3 + tag_count:
            raise ValueError(f"malformed MSH2 tags on element row {row_ordinal}")
        tags = fields[3:3 + tag_count]
        corner_tags = fields[3 + tag_count:]
        elements.append({
            "row_ordinal": row_ordinal,
            "element_id": element_id,
            "element_type": element_type,
            "tags": tags,
            "physical_tag": tags[0] if tags else None,
            "geometrical_entity_tag": tags[1] if len(tags) > 1 else None,
            "node_tags": corner_tags,
            "nodes": [node_row[tag] for tag in corner_tags],
        })
    if len(elements) != element_count or len({e["element_id"] for e in elements}) != element_count:
        raise ValueError("MSH2 element count or IDs are invalid")

    _, _, names_block = section(raw, "PhysicalNames")
    name_lines = names_block.splitlines()
    names: dict[tuple[int, int], str] = {}
    for line in name_lines[2:-1]:
        dim_text, tag_text, _ = line.split(maxsplit=2)
        dim, tag = int(dim_text), int(tag_text)
        name = shlex.split(line.decode("ascii"), posix=True)[2]
        names[(dim, tag)] = name
    return {
        "raw": raw,
        "nodes_block": section(raw, "Nodes")[2],
        "elements_block": section(raw, "Elements")[2],
        "physical_names_block": names_block,
        "node_tags": np.asarray(node_tags, dtype=np.int64),
        "nodes": np.asarray(nodes, dtype=np.float64),
        "elements": elements,
        "physical_names": names,
        "tetrahedra": [e for e in elements if e["element_type"] == 4],
        "triangles": [e for e in elements if e["element_type"] == 2],
    }


def triangle_area(nodes: np.ndarray, ids: np.ndarray) -> float:
    xyz = nodes[np.asarray(ids, dtype=np.int64)]
    return float(0.5 * np.linalg.norm(np.cross(xyz[1] - xyz[0], xyz[2] - xyz[0])))


def positive_entropy_weights(points_yz: np.ndarray, base: np.ndarray, target_yz: np.ndarray) -> np.ndarray:
    """Recompute the positive relative-entropy patch projection from inputs."""
    points = np.asarray(points_yz, dtype=np.float64)
    base = np.asarray(base, dtype=np.float64)
    target = np.asarray(target_yz, dtype=np.float64)
    scale = max(float(np.ptp(points, axis=0).max()), 1e-12)
    q = (points - target) / scale
    lam = np.zeros(2)
    logbase = np.log(base / base.sum())

    def state(value: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        logits = logbase + q @ value
        logits -= logits.max()
        weights = np.exp(logits)
        weights /= weights.sum()
        mean = weights @ q
        covariance = (q - mean).T @ (weights[:, None] * (q - mean))
        return weights, mean, covariance

    for _ in range(100):
        weights, mean, covariance = state(lam)
        if np.linalg.norm(mean, ord=np.inf) * scale <= 1e-11:
            if np.any(weights <= 0.0) or not np.isclose(weights.sum(), 1.0, atol=1e-13):
                raise ValueError("independent positive load projection is not normalized")
            return weights
        direction = np.linalg.solve(covariance, mean)
        error = float(mean @ mean)
        fraction = 1.0
        for _ in range(35):
            candidate = lam - fraction * direction
            _, next_mean, _ = state(candidate)
            if float(next_mean @ next_mean) < error:
                lam = candidate
                break
            fraction *= 0.5
        else:
            raise ValueError("independent positive load projection failed to converge")
    raise ValueError("independent positive load projection exceeded iterations")


def history_forecast(steps: int, facets: int, node_count: int, patch_node_count: int, node_stride: int) -> dict[str, Any]:
    qp_count = facets * 3
    qp_frames = steps + 1
    saved_steps = np.unique(np.r_[np.arange(0, steps + 1, node_stride), steps]).astype(np.int64)
    shapes = {
        "opening_m": (qp_frames, qp_count, 3),
        "traction_Pa": (qp_frames, qp_count, 3),
        "damage": (qp_frames, qp_count),
        "delta_max_m": (qp_frames, qp_count),
        "contact_traction_Pa": (qp_frames, qp_count, 3),
        "contact_opening_m": (qp_frames, qp_count, 3),
        "contact_energy_check": (qp_frames, 3),
        "interface_work_check": (qp_frames, 4),
        "energy_J": (qp_frames, 5),
        "ledger": (qp_frames, 28),
        "patch_displacement_x_m": (qp_frames, patch_node_count),
        "node_state": (len(saved_steps), node_count, 4, 3),
    }
    arrays = {name: {"shape": list(shape), "elements": int(np.prod(shape)), "bytes": int(np.prod(shape)) * 8,
                     "runner_forecast_bytes_including_256_byte_allowance": int(np.prod(shape)) * 8 + 256}
              for name, shape in shapes.items()}
    forecast = sum(item["runner_forecast_bytes_including_256_byte_allowance"] for item in arrays.values())
    return {
        "node_stride": node_stride,
        "qp_frames_at_stride_independent_all_qp_states": qp_frames,
        "node_state_frames": len(saved_steps),
        "saved_step_first_last": [int(saved_steps[0]), int(saved_steps[-1])],
        "array_layout": arrays,
        "array_bytes_without_allowances": sum(item["bytes"] for item in arrays.values()),
        "runner_forecast_bytes": forecast,
        "runner_guard_forecast_plus_2m_bytes": forecast + 2_000_000,
        "history_cap_bytes": HISTORY_CAP,
        "within_declared_history_cap": forecast + 2_000_000 <= HISTORY_CAP,
        "scope": "exact float64 array-shape forecast from saved native counts, selected step count, patch node count, and requested node stride; no history allocation",
    }


def main() -> int:
    checks: list[dict[str, Any]] = []
    hard_failures: list[str] = []
    limitations: list[str] = []

    def add_check(name: str, outcome: str, detail: Any = None) -> None:
        checks.append({"name": name, "outcome": outcome, "detail": detail})
        if outcome == "FAIL":
            hard_failures.append(name)
        if outcome == "LIMITATION":
            limitations.append(name)

    report: dict[str, Any] = {
        "schema": "pgfrag.reference_v6_native_initialization_preflight_audit.v1",
        "status": "AUDIT_INCOMPLETE",
        "classification": "COST_INITIALIZATION_ONLY_UNVALIDATED; PREINSERTED_WEAK_SEAM_UNCALIBRATED",
        "solver_steps_executed": 0,
        "retry_or_retune_performed": False,
        "checks": checks,
        "hard_failures": hard_failures,
        "limitations": limitations,
    }
    try:
        case_path = INPUTS / "case.json"
        case = json.loads(case_path.read_text(encoding="utf-8"))
        result = json.loads((RUN / "result.json").read_text(encoding="utf-8")) if (RUN / "result.json").exists() else {}
        execution_path = DATA / "execution_record.json"
        execution = json.loads(execution_path.read_text(encoding="utf-8")) if execution_path.exists() else {}
        watchdog_path = DATA / "process_watchdog.json"
        watchdog = json.loads(watchdog_path.read_text(encoding="utf-8")) if watchdog_path.exists() else {}
        started = json.loads((RUN / "started.json").read_text(encoding="utf-8")) if (RUN / "started.json").exists() else {}
        run_receipt = json.loads((RUN / "receipt.json").read_text(encoding="utf-8")) if (RUN / "receipt.json").exists() else {}
        bound_report = json.loads((RUN / "initial_bound.json").read_text(encoding="utf-8")) if (RUN / "initial_bound.json").exists() else {}
        load_report = json.loads((RUN / "loading.json").read_text(encoding="utf-8")) if (RUN / "loading.json").exists() else {}
        lineage = json.loads((RUN / "lineage.json").read_text(encoding="utf-8")) if (RUN / "lineage.json").exists() else {}
        provenance = json.loads((INPUTS / "source_geometry_provenance.json").read_text(encoding="utf-8"))
        source_result = json.loads((INPUTS / "source_mesh_result.json").read_text(encoding="utf-8"))

        original_msh = parse_msh22(INPUTS / "mesh.original.msh")
        alias_msh = parse_msh22(INPUTS / "mesh.msh")
        original_bytes = original_msh["raw"]
        alias_bytes = alias_msh["raw"]
        name_start, name_end, original_names_block = section(original_bytes, "PhysicalNames")
        alias_names_start, alias_names_end, alias_names_block = section(alias_bytes, "PhysicalNames")
        expected_alias_names_block = original_names_block.replace(b'2 2 "cut_interface"', b'2 2 "interface"').replace(b'2 3 "source_exterior"', b'2 3 "exterior"')
        expected_alias_bytes = original_bytes[:name_start] + expected_alias_names_block + original_bytes[name_end:]
        source_name_contract = {(2, 2): "cut_interface", (2, 3): "source_exterior", (3, 1): "partitioned_body"}
        alias_name_contract = {(2, 2): "interface", (2, 3): "exterior", (3, 1): "partitioned_body"}
        original_signature = [(e["element_id"], e["element_type"], e["tags"], e["node_tags"]) for e in original_msh["elements"]]
        alias_signature = [(e["element_id"], e["element_type"], e["tags"], e["node_tags"]) for e in alias_msh["elements"]]
        add_check("physical_name_alias_changes_only_the_two_declared_names", alias_bytes == expected_alias_bytes, {
            "expected_names": alias_name_contract, "actual_names": alias_msh["physical_names"],
            "original_and_alias_names_blocks_differ": original_names_block != alias_names_block,
        })
        add_check("raw_nodes_and_elements_blocks_byte_identical", original_msh["nodes_block"] == alias_msh["nodes_block"] and original_msh["elements_block"] == alias_msh["elements_block"], {
            "nodes_bytes": len(original_msh["nodes_block"]), "elements_bytes": len(original_msh["elements_block"]),
        })
        add_check("parsed_msh_connectivity_integer_tags_element_ids_and_order_identical", original_signature == alias_signature, {
            "elements": len(original_signature), "physical_names_original": original_msh["physical_names"],
            "physical_names_alias": alias_msh["physical_names"],
            "volume_physical_group_unchanged": original_msh["physical_names"].get((3, 1)) == alias_msh["physical_names"].get((3, 1)),
        })
        add_check("msh_original_and_alias_physical_name_contract", original_msh["physical_names"] == source_name_contract and alias_msh["physical_names"] == alias_name_contract)

        msh_tets = original_msh["tetrahedra"]
        msh_exterior = [e for e in original_msh["triangles"] if e["physical_tag"] == 3]
        msh_interface = [e for e in original_msh["triangles"] if e["physical_tag"] == 2]
        npz_path = INPUTS / "source_mesh.npz"
        with np.load(npz_path, allow_pickle=False) as npz:
            npz_nodes = np.asarray(npz["nodes_m"], dtype=np.float64)
            npz_node_tags = np.asarray(npz["original_node_tags"], dtype=np.int64)
            npz_tets = np.asarray(npz["tets"], dtype=np.int64)
            npz_cell_tags = np.asarray(npz["original_cell_tags"], dtype=np.int64)
            npz_interface = np.asarray(npz["interface_triangles"], dtype=np.int64)
            npz_exterior = np.asarray(npz["exterior_triangles"], dtype=np.int64)
            npz_exterior_entities = np.asarray(npz["exterior_surface_tags"], dtype=np.int64)
            npz_exterior_fids = np.asarray(npz["exterior_source_face_ids"]).astype(str)
            npz_exterior_roles = np.asarray(npz["exterior_roles"]).astype(str)
            npz_exterior_exact = np.asarray(npz["exterior_ancestry_exact"], dtype=bool)
        npz_referenced = np.unique(npz_tets)
        npz_unused_ordinals = np.setdiff1d(np.arange(len(npz_nodes), dtype=np.int64), npz_referenced)
        npz_unused_tags = npz_node_tags[npz_unused_ordinals]
        msh_tag_to_row = {int(tag): row for row, tag in enumerate(original_msh["node_tags"])}
        npz_tag_to_row = {int(tag): row for row, tag in enumerate(npz_node_tags)}
        msh_tag_set = set(msh_tag_to_row)
        npz_referenced_tags = set(npz_node_tags[npz_referenced].tolist())
        unused_present = sorted(set(npz_unused_tags.tolist()) & msh_tag_set)
        add_check("msh_nodes_are_exactly_npz_t4_referenced_original_node_tags", msh_tag_set == npz_referenced_tags, {
            "msh_nodes": len(msh_tag_set), "npz_total_nodes": len(npz_nodes),
            "npz_t4_referenced_nodes": len(npz_referenced), "npz_unused_node_count": len(npz_unused_tags),
            "npz_unused_original_node_tags": [int(x) for x in npz_unused_tags],
            "npz_unused_tags_present_in_msh": unused_present,
        })
        coordinate_deltas = []
        coordinate_relatives = []
        for tag in msh_tag_set:
            left = original_msh["nodes"][msh_tag_to_row[tag]]
            right = npz_nodes[npz_tag_to_row[tag]]
            delta = np.abs(left - right)
            coordinate_deltas.append(float(delta.max(initial=0.0)))
            reference_norm = float(np.linalg.norm(right))
            if reference_norm > 0.0:
                coordinate_relatives.append(float(np.linalg.norm(delta) / reference_norm))
            if not np.allclose(left, right, rtol=8e-16, atol=1e-30):
                raise ValueError(f"NPZ↔MSH coordinate tolerance exceeded for original node tag {tag}")
        max_coord_abs = max(coordinate_deltas, default=0.0)
        max_coord_rel = max(coordinate_relatives, default=0.0)
        add_check("npz_msh_coordinates_within_ascii_serialization_tolerance", max_coord_abs <= 6.94e-18 or np.allclose(
            np.vstack([original_msh["nodes"][msh_tag_to_row[int(tag)]] for tag in sorted(msh_tag_set)]),
            np.vstack([npz_nodes[npz_tag_to_row[int(tag)]] for tag in sorted(msh_tag_set)]), rtol=8e-16, atol=1e-30), {
                "rtol": 8e-16, "atol_m": 1e-30, "max_abs_error_m": max_coord_abs,
                "max_vector_relative_error": max_coord_rel,
            })

        npz_cell_row_by_id = {int(tag): row for row, tag in enumerate(npz_cell_tags)}
        msh_tet_ids = np.asarray([e["element_id"] for e in msh_tets], dtype=np.int64)
        msh_to_npz_rows = np.asarray([npz_cell_row_by_id[int(tag)] for tag in msh_tet_ids], dtype=np.int64)
        ordered_corner_tags = np.asarray([e["node_tags"] for e in msh_tets], dtype=np.int64)
        mapped_npz_corner_tags = npz_node_tags[npz_tets[msh_to_npz_rows]]
        cell_map = np.load(INPUTS / "t4_cell_row_correspondence.npz", allow_pickle=False)
        map_rows_equal = np.array_equal(cell_map["msh_tet_raw_row_ordinals"], np.asarray([e["row_ordinal"] for e in msh_tets]))
        map_ids_equal = np.array_equal(cell_map["msh_tet_element_ids"], msh_tet_ids)
        map_npz_rows_equal = np.array_equal(cell_map["npz_tet_rows_by_original_cell_tag"], msh_to_npz_rows)
        map_corner_tags_equal = np.array_equal(cell_map["ordered_original_node_tags_by_corner"], ordered_corner_tags)
        ordered_cell_match = np.array_equal(mapped_npz_corner_tags, ordered_corner_tags)
        add_check("t4_cell_correspondence_by_original_cell_tag_and_ordered_corner_node_tags", (
            len(set(msh_tet_ids.tolist())) == len(msh_tet_ids) == len(npz_cell_tags)
            and set(msh_tet_ids.tolist()) == set(npz_cell_tags.tolist())
            and ordered_cell_match and map_rows_equal and map_ids_equal and map_npz_rows_equal and map_corner_tags_equal
        ), {
            "t4_count": len(msh_tets), "msh_raw_row_ordinals": [int(e["row_ordinal"]) for e in msh_tets[:3]],
            "first_msh_t4_element_ids": msh_tet_ids[:3].tolist(),
            "first_npz_rows_by_original_cell_tag": msh_to_npz_rows[:3].tolist(),
            "row_correspondence_is_identity": bool(np.array_equal(msh_to_npz_rows, np.arange(len(msh_to_npz_rows)))),
            "all_ordered_corner_tags_match": ordered_cell_match,
            "proof_npz_sha256": sha256(INPUTS / "t4_cell_row_correspondence.npz"),
        })

        # Independently rebind every exterior MSH triangle to exact OCC FID/role ancestry.
        ancestry = source_result["source_face_ancestry"]["cut_witness_ancestry"]
        expected_entity = {int(k): v for k, v in ancestry["source_by_exterior_surface"].items()}
        rows_doc = json.loads((INPUTS / "known_exterior_fid_role_rows.json").read_text(encoding="utf-8"))
        row_records = rows_doc["rows"]
        key_to_row = {tuple(sorted(int(tag) for tag in row["node_tags"])): row for row in row_records}
        exterior_bindings_good = len(msh_exterior) == len(row_records) == len(npz_exterior) == len(npz_exterior_exact) == int(source_result["counts"]["exterior_triangles"])
        bound_ids: set[int] = set()
        for row in msh_exterior:
            raw_key = tuple(sorted(int(tag) for tag in row["node_tags"]))
            record = key_to_row.get(raw_key)
            entity_tag = int(row["geometrical_entity_tag"])
            expected = expected_entity.get(entity_tag)
            if record is None or expected is None:
                exterior_bindings_good = False
                continue
            exterior_bindings_good &= (
                int(record["element_id"]) == int(row["element_id"])
                and int(record["row_ordinal"]) == int(row["row_ordinal"])
                and int(record["physical_tag"]) == 3
                and int(record["geometrical_entity_tag"]) == entity_tag
                and record["source_face_id"] == expected["source_face_id"]
                and record["role"] == expected["role"]
                and record["node_tags"] == row["node_tags"]
            )
            bound_ids.add(int(row["element_id"]))
        array_path = INPUTS / "known_exterior_fid_roles.npz"
        with np.load(array_path, allow_pickle=False) as arrays:
            ext_ids = arrays["source_exterior_element_ids"].astype(np.int64)
            ext_rows = arrays["source_exterior_raw_element_row_ordinals"].astype(np.int64)
            ext_phys = arrays["source_exterior_physical_tags"].astype(np.int64)
            ext_geom = arrays["source_exterior_geometrical_entity_tags"].astype(np.int64)
            ext_tags = arrays["source_exterior_node_tags"].astype(np.int64)
            fid_strings = arrays["source_face_ids"].astype(str)
            role_names = arrays["role_names"].astype(str)
            ext_fids = fid_strings[arrays["source_exterior_fid_index"].astype(np.int64)]
            ext_roles = role_names[arrays["source_exterior_role_index"].astype(np.int64)]
            for i, row in enumerate(msh_exterior):
                expected_face = expected_entity.get(int(row["geometrical_entity_tag"]), {})
                exterior_bindings_good &= (
                    int(ext_ids[i]) == int(row["element_id"])
                    and int(ext_rows[i]) == int(row["row_ordinal"])
                    and int(ext_phys[i]) == int(row["physical_tag"]) == 3
                    and int(ext_geom[i]) == int(row["geometrical_entity_tag"])
                    and np.array_equal(ext_tags[i], np.asarray(row["node_tags"], dtype=np.int64))
                    and str(ext_fids[i]) == str(expected_face.get("source_face_id"))
                    and str(ext_roles[i]) == str(expected_face.get("role"))
                )
            exterior_bindings_good &= bool(npz_exterior_exact.all()) and ancestry["interface_has_zero_source_face_ancestry"] is True
        add_check("exterior_fid_role_arrays_bind_by_exact_source_node_tags_and_msh_entity_tags", exterior_bindings_good, {
            "exterior_rows": len(msh_exterior), "verified_bound_element_ids": len(bound_ids),
            "known_face_ids": len(set(row["source_face_id"] for row in row_records)),
            "exact_ancestry_rows": int(np.count_nonzero(npz_exterior_exact)),
            "interface_has_zero_source_face_ancestry": ancestry["interface_has_zero_source_face_ancestry"],
            "nearest_or_coordinate_matching_used": False,
            "known_array_sha256": sha256(array_path),
            "known_rows_sha256": sha256(INPUTS / "known_exterior_fid_role_rows.json"),
        })

        # Immutable input/source chain and exact process/runtime evidence.
        hash_mismatches = []
        for name, expected_hash in case["input_hashes"].items():
            path = INPUTS / name
            actual = sha256(path) if path.exists() else None
            if actual != expected_hash:
                hash_mismatches.append({"name": name, "expected": expected_hash, "actual": actual})
        source_hashes = {
            "canonical_step": ROOT / provenance["canonical_step"]["path"],
            "run05_msh": INPUTS / "mesh.original.msh",
            "run05_npz": INPUTS / "source_mesh.npz",
            "run05_result": INPUTS / "source_mesh_result.json",
        }
        source_hash_check = []
        for key, path in source_hashes.items():
            expected_hash = (provenance["canonical_step"]["sha256"] if key == "canonical_step"
                             else provenance["mesh_sources"][{"run05_msh": "mesh_msh", "run05_npz": "mesh_npz", "run05_result": "result_json"}[key]]["sha256"])
            actual = sha256(path) if path.exists() else None
            source_hash_check.append({"key": key, "path": str(path), "expected": expected_hash, "actual": actual, "match": actual == expected_hash})
        add_check("all_frozen_input_hashes_unchanged", not hash_mismatches, {"mismatches": hash_mismatches})
        add_check("canonical_step_and_v6_mesh_source_hashes_match_frozen_provenance", all(x["match"] for x in source_hash_check), source_hash_check)
        module_hashes = {}
        for name in ("cohesive_vessel.py", "cohesive_mesh.py", "cohesive_stability.py"):
            current = ROOT / "pgfrag/src/pgfrag" / name
            frozen = INPUTS / name
            module_hashes[name] = {
                "current_sha256": sha256(current), "frozen_sha256": sha256(frozen),
                "case_expected_sha256": case["input_hashes"][name],
                "unchanged": sha256(current) == sha256(frozen) == case["input_hashes"][name],
            }
        add_check("native_source_modules_unchanged_from_frozen_inputs", all(v["unchanged"] for v in module_hashes.values()), module_hashes)
        add_check("execution_record_and_watchdog_confirm_single_native_preflight", (
            execution.get("retry_count") == 0 and execution.get("retry_or_retune_performed") is False
            and execution.get("exit_code") == 0 and execution.get("watchdog_timed_out") is False
            and watchdog.get("external_watchdog_s") == 330 and watchdog.get("current_module_bytes_match_frozen_copies") is True
            and execution.get("runner_scripts_unchanged") is True and watchdog.get("runner_scripts_unchanged") is True
        ), {
            "runtime": execution.get("runtime"), "runtime_probe": execution.get("runtime_probe"),
            "command_argv": execution.get("command_argv"), "exit_code": execution.get("exit_code"),
            "watchdog_timed_out": execution.get("watchdog_timed_out"), "elapsed_wall_s": execution.get("elapsed_wall_s"),
            "retry_count": execution.get("retry_count"), "thread_environment": execution.get("thread_environment"),
            "cpu_affinity": execution.get("native_cpu_affinity"),
            "runner_script_hashes_before": execution.get("runner_script_hashes_before"),
            "runner_script_hashes_after": execution.get("runner_script_hashes_after"),
            "runner_scripts_unchanged": execution.get("runner_scripts_unchanged"),
        })
        add_check("result_is_zero_step_initialization_only", result.get("status") == "ZERO_STEP_PREFLIGHT_PASS" and result.get("solver_steps") == 0, {
            "native_result": result, "run_directory": str(RUN), "history_directory_exists": (RUN / "history").exists(),
        })
        add_check("no_solver_history_or_evolution_outputs_created", not (RUN / "history").exists() and not (RUN / "completion.json").exists(), {
            "history_directory_exists": (RUN / "history").exists(), "completion_exists": (RUN / "completion.json").exists(),
        })

        state_path = RUN / "initial_state.npz"
        state_data: dict[str, np.ndarray] = {}
        with np.load(state_path, allow_pickle=False) as archive:
            state_data = {key: np.asarray(archive[key]) for key in archive.files}
        rest = np.asarray(state_data["rest_nodes_m"], dtype=np.float64)
        native_tets = np.asarray(state_data["tets"], dtype=np.int64)
        cohesive = np.asarray(state_data["cohesive"], dtype=np.int64)
        mass = np.asarray(state_data["mass_kg"], dtype=np.float64)
        blocked = np.asarray(state_data["blocked"], dtype=bool)
        native_to_source = np.asarray(state_data["native_to_source"], dtype=np.int64)
        state_source_tags = np.asarray(state_data["source_node_tags"], dtype=np.int64)
        source_tet_ids = np.asarray(state_data["source_tet_ids"], dtype=np.int64)
        source_exterior_faces = np.asarray(state_data["source_exterior_faces"], dtype=np.int64)
        pairs = np.asarray(state_data["cohesive_pairs"], dtype=np.int64)
        qp_areas = np.asarray(state_data["qp_area_m2"], dtype=np.float64)
        unit_load = np.asarray(state_data["unit_external_load"], dtype=np.float64)

        input_tets = np.asarray([e["nodes"] for e in alias_msh["tetrahedra"]], dtype=np.int64)
        input_tet_ids = np.asarray([e["element_id"] for e in alias_msh["tetrahedra"]], dtype=np.int64)
        native_initial_map = lineage.get("native_node_to_input_node_zero_based", {})
        local_rows = lineage.get("native_bulk_local_ordinal_lineage_rows", [])
        lineage_row_identity = len(local_rows) == len(input_tets) * 4
        reconstructed_native_tets = np.full_like(input_tets, -1)
        reconstructed_source_tets = np.full_like(input_tets, -1)
        for row in local_rows:
            bulk_row = int(row["native_bulk_element_row"])
            local = int(row["local_vertex_ordinal"])
            if not (0 <= bulk_row < len(input_tets) and 0 <= local < 4):
                lineage_row_identity = False
                continue
            lineage_row_identity &= (
                int(row["input_tetrahedron_id"]) == int(input_tet_ids[bulk_row])
                and int(row["input_node_id_zero_based"]) == int(input_tets[bulk_row, local])
            )
            reconstructed_native_tets[bulk_row, local] = int(row["native_node_id"])
            reconstructed_source_tets[bulk_row, local] = int(row["input_node_id_zero_based"])
        state_row_identity = (
            np.array_equal(source_tet_ids, input_tet_ids)
            and np.array_equal(reconstructed_source_tets, input_tets)
            and np.array_equal(reconstructed_native_tets, native_tets)
        )
        add_check("native_initial_msh_t4_row_and_local_corner_lineage", lineage_row_identity and state_row_identity, {
            "MSH_t4_count": len(input_tets), "state_t4_count": len(native_tets),
            "source_tet_ids_exact_in_msh_order": np.array_equal(source_tet_ids, input_tet_ids),
            "native_lineage_records": len(local_rows),
            "row_local_lineage_reconstructs_native_post_insertion_tets": np.array_equal(reconstructed_native_tets, native_tets),
        })
        source_nodes = alias_msh["nodes"]
        source_tags = alias_msh["node_tags"]
        native_rows_in_bounds = native_to_source.shape == (len(rest),) and np.all(native_to_source >= 0) and np.all(native_to_source < len(source_nodes))
        ancestry_coords_exact = native_rows_in_bounds and np.array_equal(rest, source_nodes[native_to_source])
        add_check("native_rest_coordinates_and_source_node_tag_order_match_msh", (
            state_source_tags.shape == source_tags.shape and np.array_equal(state_source_tags, source_tags)
            and native_rows_in_bounds and ancestry_coords_exact
        ), {
            "MSH_source_nodes": len(source_nodes), "native_initial_nodes_after_cohesive_insertion": len(rest),
            "source_node_tags_exact_order": np.array_equal(state_source_tags, source_tags),
            "native_to_source_rows_in_bounds": native_rows_in_bounds,
            "native_rest_coordinates_exact_by_ancestry": ancestry_coords_exact,
        })
        msh_triangle_rows = [e for e in alias_msh["triangles"] if e["physical_tag"] == 3]
        source_exterior_match = np.array_equal(source_exterior_faces, np.asarray([e["nodes"] for e in msh_triangle_rows], dtype=np.int64))
        add_check("saved_source_exterior_faces_match_original_msh_row_order", source_exterior_match, {
            "saved_exterior_triangles": len(source_exterior_faces), "MSH_exterior_triangles": len(msh_triangle_rows),
        })

        # Native T4 determinant volumes and independently reconstructed density mass.
        tet_xyz = rest[native_tets]
        signed_volumes = np.linalg.det(np.transpose(tet_xyz[:, 1:] - tet_xyz[:, :1], (0, 2, 1))) / 6.0
        # np.linalg.det above uses vectors as columns; this agrees with signed tetra orientation.
        volume_sum = float(signed_volumes.sum())
        negative_or_degenerate = int(np.count_nonzero(signed_volumes <= 0.0))
        expected_mass = np.zeros(len(rest), dtype=np.float64)
        for row, volume in zip(native_tets, signed_volumes, strict=True):
            expected_mass[row] += float(case["bulk_material"]["rho_kg_m3"]) * float(volume) / 4.0
        mass_shape_ok = mass.shape == (len(rest), 3) and blocked.shape == (len(rest), 3)
        mass_component_equal = mass_shape_ok and np.allclose(mass, mass[:, :1], rtol=0.0, atol=0.0)
        mass_relative_max = float(np.max(np.abs(mass[:, 0] - expected_mass) / np.maximum(expected_mass, 1e-300))) if mass_shape_ok else None
        zero_mass_nodes = np.flatnonzero((mass[:, 0] <= 0.0) if mass_shape_ok else np.asarray([], dtype=bool))
        unblocked_zero_mass_dof_count = int(np.count_nonzero((mass <= 0.0) & (~blocked))) if mass_shape_ok else None
        expected_total_mass = float(case["bulk_material"]["rho_kg_m3"]) * volume_sum
        actual_total_mass = float(mass[:, 0].sum()) if mass_shape_ok else None
        mass_ok = bool(mass_shape_ok and mass_component_equal and np.allclose(mass[:, 0], expected_mass, rtol=1e-8, atol=1e-300)
                       and np.isclose(actual_total_mass, expected_total_mass, rtol=1e-10, atol=0.0))
        volume_expected = float(source_result["mesh"]["tetrahedral_volume_m3"])
        add_check("native_t4_signed_volumes_positive_and_match_source_mesh", negative_or_degenerate == 0 and np.isclose(volume_sum, volume_expected, rtol=1e-10, atol=1e-18), {
            "positive_t4": int(np.count_nonzero(signed_volumes > 0.0)), "nonpositive_t4": negative_or_degenerate,
            "min_signed_volume_m3": float(signed_volumes.min()), "sum_signed_volume_m3": volume_sum,
            "run05_expected_mesh_volume_m3": volume_expected,
        })
        add_check("native_lumped_mass_matches_independent_rho_volume_mass", mass_ok, {
            "density_kg_m3": case["bulk_material"]["rho_kg_m3"], "expected_total_mass_kg": expected_total_mass,
            "native_total_mass_kg": actual_total_mass, "maximum_node_relative_mass_error": mass_relative_max,
            "native_mass_components_equal": mass_component_equal, "mass_shape": list(mass.shape),
        })
        dofs_free = bool(mass_shape_ok and not blocked.any())
        add_check("native_dof_mask_is_all_free", dofs_free, {
            "native_nodes": len(rest), "native_dof_count": int(blocked.size),
            "blocked_dof_count": int(np.count_nonzero(blocked)),
        })
        add_check("no_unblocked_native_zero_mass_dofs", unblocked_zero_mass_dof_count == 0, {
            "zero_mass_node_count": len(zero_mass_nodes),
            "zero_mass_native_node_ids": zero_mass_nodes.tolist(),
            "unblocked_zero_mass_dof_count": unblocked_zero_mass_dof_count,
        })
        if unblocked_zero_mass_dof_count:
            add_check("future_evolution_mass_eligibility", "LIMITATION", "Unblocked zero-mass DOFs must block any future evolution; no automatic pruning or BC change is permitted.")

        native_unused_matches = [int(tag) for tag in npz_unused_tags if int(tag) in set(map(int, state_source_tags))]
        add_check("npz_unused_original_tags_not_silently_pruned_or_promoted_to_native_nodes", len(npz_unused_tags) == 169 and not native_unused_matches, {
            "npz_unused_count": len(npz_unused_tags), "npz_unused_original_node_tags": [int(x) for x in npz_unused_tags],
            "unused_tags_in_original_msh": unused_present,
            "unused_tags_retained_as_native_source_nodes": native_unused_matches,
            "native_mesh_node_count_observed": len(state_source_tags),
            "native_interpretation": "Akantu reads the preserved MSH node section; the 169 NPZ-only tags are absent from that section and have no native row, mass, or DOF. No mass is inferred for them.",
        })
        if len(npz_unused_tags) != 0:
            add_check("unused_npz_nodes_are_not_part_of_native_evolution_mesh", "LIMITATION", {
                "count": len(npz_unused_tags), "tags": [int(x) for x in npz_unused_tags],
                "status": "not represented in source MSH or native state; no topology mutation or renumbering performed",
            })

        # Pair rest coordinates, facet area, and three equal quadrature weights.
        pairs_shape_ok = pairs.shape == (len(cohesive), 2, 3)
        pair_side_coordinates_equal = pairs_shape_ok and np.array_equal(rest[pairs[:, 0]], rest[pairs[:, 1]])
        cohesive_side_a = cohesive[:, :3]
        cohesive_side_b = cohesive[:, 3:]
        pair_connectivity_equal = pairs_shape_ok and np.array_equal(pairs[:, 0], cohesive_side_a) and np.array_equal(pairs[:, 1], cohesive_side_b)
        facet_areas = np.asarray([triangle_area(rest, face) for face in pairs[:, 0]], dtype=np.float64) if pairs_shape_ok else np.zeros(0)
        side_b_areas = np.asarray([triangle_area(rest, face) for face in pairs[:, 1]], dtype=np.float64) if pairs_shape_ok else np.zeros(0)
        qp_expected = np.repeat(facet_areas / 3.0, 3)
        qp_match = qp_areas.shape == qp_expected.shape and np.allclose(qp_areas, qp_expected, rtol=1e-12, atol=1e-18)
        qp_area_sum = float(qp_areas.sum())
        area_report = float(bound_report.get("interface_area_m2", np.nan))
        add_check("cohesive_pair_rest_coordinates_and_tri3_quadrature_areas", (
            pair_connectivity_equal and pair_side_coordinates_equal
            and np.allclose(facet_areas, side_b_areas, rtol=0.0, atol=1e-18)
            and qp_match and np.isclose(qp_area_sum, float(facet_areas.sum()), rtol=1e-12, atol=1e-18)
        ), {
            "cohesive_tri3_count": len(cohesive), "paired_rest_coordinates_exact": pair_side_coordinates_equal,
            "saved_pairs_follow_native_cohesive_connectivity": pair_connectivity_equal,
            "qp_count": len(qp_areas), "qp_area_formula": "independent Tri3 rest area / 3 at each of three points",
            "sum_facet_area_m2": float(facet_areas.sum()), "sum_qp_area_m2": qp_area_sum,
            "bound_report_interface_area_m2": area_report,
            "max_side_area_difference_m2": float(np.max(np.abs(facet_areas - side_b_areas))) if len(facet_areas) else 0.0,
        })

        # Independently rebuild the balanced unit patch from exterior MSH triangles and native mass.
        patch_report: dict[str, Any] = {"available": False}
        if mass_shape_ok and native_rows_in_bounds:
            source_weights = np.zeros(len(source_nodes), dtype=np.float64)
            for face in msh_exterior:
                area = triangle_area(source_nodes, np.asarray(face["nodes"], dtype=np.int64))
                source_weights[np.asarray(face["nodes"], dtype=np.int64)] += area / 3.0
            target = np.sum(mass[:, 0, None] * rest, axis=0) / mass[:, 0].sum()
            threshold = float(case["loading"]["patch_x_threshold_m"])
            rebuilt = np.zeros_like(rest)
            expected_patches = []
            for sign in (-1, 1):
                ids = np.flatnonzero((sign * rest[:, 0] > threshold) & (source_weights[native_to_source] > 0.0))
                base = source_weights[native_to_source[ids]]
                weights = positive_entropy_weights(rest[ids, 1:], base, target[1:])
                rebuilt[ids, 0] = sign * weights
                expected_patches.append({"sign": sign, "native_node_ids": ids, "weights": weights,
                                         "target_yz_m": (weights @ rest[ids, 1:]).tolist(),
                                         "weight_sum": float(weights.sum()), "minimum_weight": float(weights.min())})
            force_sum = rebuilt.sum(axis=0)
            torque_sum = np.cross(rest, rebuilt).sum(axis=0)
            load_match = unit_load.shape == rebuilt.shape and np.allclose(unit_load, rebuilt, rtol=2e-11, atol=1e-14)
            json_patch_match = len(load_report.get("patches", [])) == 2
            for expected, saved in zip(expected_patches, load_report.get("patches", []), strict=False):
                json_patch_match &= (expected["sign"] == int(saved["sign"])
                                     and np.array_equal(expected["native_node_ids"], np.asarray(saved["native_node_ids"], dtype=np.int64))
                                     and np.allclose(expected["weights"], np.asarray(saved["weights"], dtype=np.float64), rtol=2e-11, atol=1e-14))
            json_patch_match &= np.allclose(target, np.asarray(load_report.get("target_reference_com_m", [])), rtol=1e-12, atol=1e-14)
            balance_ok = (np.linalg.norm(force_sum) <= 1e-12 and np.linalg.norm(torque_sum) <= 1e-11
                          and all(np.all(item["weights"] > 0) and abs(item["weight_sum"] - 1.0) <= 1e-12 for item in expected_patches)
                          and all(np.linalg.norm(np.asarray(item["target_yz_m"]) - target[1:]) <= 1e-11 for item in expected_patches))
            patch_report = {
                "available": True, "peak_force_N": float(case["loading"]["peak_force_N"]),
                "pulse_duration_s": float(case["loading"]["pulse_duration_s"]), "threshold_abs_x_m": threshold,
                "reference_mass_com_xyz_m": target.tolist(), "patches": [
                    {"sign": item["sign"], "native_node_count": len(item["native_node_ids"]),
                     "native_node_ids": item["native_node_ids"].tolist(), "weight_sum": item["weight_sum"],
                     "minimum_positive_weight": item["minimum_weight"], "weighted_reference_yz_m": item["target_yz_m"]}
                    for item in expected_patches
                ],
                "independently_rebuilt_load_matches_saved_load": bool(load_match),
                "independently_rebuilt_load_matches_loading_json": bool(json_patch_match),
                "unit_force_sum_xyz": force_sum.tolist(), "unit_torque_sum_xyz_Nm": torque_sum.tolist(),
                "positive_weights_and_each_patch_COM_yz_balanced": bool(balance_ok),
                "actual_pulse_magnitude_is_zero": float(case["loading"]["peak_force_N"]) == 0.0,
            }
            add_check("positive_balanced_com_yz_patches_and_zero_magnitude_pulse", load_match and json_patch_match and balance_ok and patch_report["actual_pulse_magnitude_is_zero"], patch_report)

        # Recompute the initial stability bound from saved native arrays; use the observed Akantu bulk dt.
        stability_path = INPUTS / "cohesive_stability.py"
        sys.path.insert(0, str(INPUTS))
        from cohesive_stability import mass_normalized_stiffness_bound  # noqa: E402

        material = case["bulk_material"]
        interface = case["interface"]
        kn = float(interface["sigma_c_Pa"]) / float(interface["delta0_m"])
        penalty = float(interface["penalty_Pa_m"])
        recomputed_bound = mass_normalized_stiffness_bound(
            rest, native_tets, pairs, float(material["rho_kg_m3"]), float(material["E_pa"]), float(material["nu"]),
            kn, kn, blocked, mass,
        )
        critical_dt = float(recomputed_bound.critical_time_step_s)
        bulk_dt = float(bound_report["bulk_only_native_dt_s"])
        max_dt = min(float(case["time"]["initial_bound_fraction"]) * critical_dt, 0.4 * bulk_dt)
        steps = math.ceil(float(case["time"]["total_s"]) / max_dt) * int(case["time"].get("refinement_factor", 1))
        selected_dt = float(case["time"]["total_s"]) / steps
        bound_match = (
            np.isclose(float(bound_report["critical_time_step_s"]), critical_dt, rtol=1e-12, atol=0.0)
            and np.isclose(float(bound_report["dt_s"]), selected_dt, rtol=1e-12, atol=0.0)
            and int(bound_report["steps"]) == steps
            and np.isclose(float(bound_report["lambda_upper_bound_per_s2"]), float(recomputed_bound.lambda_upper_bound_per_s2), rtol=1e-12, atol=0.0)
        )
        add_check("initial_positive_stiffness_bound_and_selected_timestep_recomputed", bound_match, {
            "lambda_upper_bound_per_s2": float(recomputed_bound.lambda_upper_bound_per_s2),
            "critical_time_step_s": critical_dt, "observed_akantu_bulk_dt_s": bulk_dt,
            "selected_dt_s": selected_dt, "selected_steps_for_50us": steps,
            "initial_bound_fraction": case["time"]["initial_bound_fraction"],
            "penalty_Pa_per_m": penalty, "kn_Pa_per_m": kn,
            "penalty_le_kn": penalty <= kn,
            "scope": "initial positive-stiffness bound only; not a full-trajectory guarantee",
        })
        history50 = history_forecast(steps, len(cohesive), len(rest), int(np.count_nonzero(unit_load[:, 0])), 50)
        history1 = history_forecast(steps, len(cohesive), len(rest), int(np.count_nonzero(unit_load[:, 0])), 1)
        add_check("full_history_stride_50_byte_forecast_recorded_without_allocation", "PASS" if history50["within_declared_history_cap"] else "LIMITATION", history50)
        add_check("full_history_stride_1_byte_forecast_recorded_without_allocation", "PASS" if history1["within_declared_history_cap"] else "LIMITATION", history1)

        output_bytes = tree_bytes(EXPERIMENT) + tree_bytes(DATA)
        actual_lineage_bytes = (RUN / "lineage.json").stat().st_size if (RUN / "lineage.json").exists() else None
        budget = json.loads((DATA / "preflight_budget.json").read_text(encoding="utf-8"))
        forecast = int(budget["native_preflight_output_forecast"]["lineage_forecast_bytes"])
        add_check("whole_workspace_outputs_within_256mib_and_lineage_forecast", (
            output_bytes <= ALL_OUTPUT_CAP and actual_lineage_bytes is not None and actual_lineage_bytes <= forecast
        ), {
            "measured_assigned_experiment_and_data_bytes": output_bytes, "cap_bytes": ALL_OUTPUT_CAP,
            "native_lineage_json_bytes": actual_lineage_bytes, "lineage_forecast_bytes": forecast,
            "input_output_cap_bytes": int(budget["native_preflight_output_forecast"]["whole_inputs_and_outputs_cap_bytes"]),
        })

        report.update({
            "status": "PASS_NATIVE_INITIALIZATION_ONLY" if not hard_failures else "FAIL_PRESERVED_NO_RETRY",
            "native_result_status": result.get("status"),
            "native_solver_steps_executed": result.get("solver_steps"),
            "native_initial_counts_reported": result.get("counts"),
            "native_initial_state_file": {"path": str(state_path), "sha256": sha256(state_path), "bytes": state_path.stat().st_size},
            "native_initialization_only_zero_field_gate": {
                "outcome": "PASSED_BY_FROZEN_RUNNER_GATE" if result.get("status") == "ZERO_STEP_PREFLIGHT_PASS" else "NOT_PASSED",
                "evidence": "runner checked opening, traction, damage, contact opening and contact traction max absolute value <= 1e-15 before preflight return",
                "field_arrays_serialized": False,
                "limitation": "current frozen runner does not include cohesive initial-field arrays in initial_state.npz; actual values cannot be re-read from saved artifacts",
            },
            "source_mesh_counts": {"msh_nodes": len(original_msh["nodes"]), "T4": len(msh_tets),
                                   "interface_Tri3": len(msh_interface), "exterior_Tri3": len(msh_exterior)},
            "npz_only_unused_nodes": {
                "count": len(npz_unused_tags), "original_node_tags": [int(x) for x in npz_unused_tags],
                "present_in_source_msh": unused_present, "present_as_native_source_rows": native_unused_matches,
                "retained_by_akantu": False if not native_unused_matches else True,
                "mass_and_dof_status": "NO_NATIVE_MASS_OR_DOFS; original tags are absent from the MSH Nodes block and do not appear in native source_node_tags",
            },
            "native_initial_state": {
                "native_nodes_after_cohesive_insertion": len(rest), "native_T4": len(native_tets),
                "cohesive_Tri3": len(cohesive), "free_dofs": int(blocked.size - np.count_nonzero(blocked)),
                "blocked_dofs": int(np.count_nonzero(blocked)), "unblocked_zero_mass_dofs": unblocked_zero_mass_dof_count,
                "native_mass_total_kg": actual_total_mass, "independent_rhoV_mass_kg": expected_total_mass,
                "native_mesh_volume_m3": volume_sum, "signed_positive_t4_count": int(np.count_nonzero(signed_volumes > 0.0)),
                "max_native_node_mass_relative_error": mass_relative_max,
                "cohesive_rest_area_m2": float(facet_areas.sum()), "qp_area_sum_m2": qp_area_sum,
            },
            "native_loading_audit": patch_report,
            "stability_and_history_audit": {"recomputed_bound_report_match": bound_match, "history_stride_50": history50, "history_stride_1": history1},
            "module_hashes": module_hashes,
            "source_hashes": source_hash_check,
            "execution_record": execution,
            "process_watchdog": watchdog,
            "run_receipt": run_receipt,
            "whole_assigned_output_bytes": output_bytes,
            "whole_output_cap_bytes": ALL_OUTPUT_CAP,
            "checks": checks,
            "hard_failures": hard_failures,
            "limitations": limitations,
            "physical_or_training_qualification": "NOT_ESTABLISHED; initialization cost diagnostic only",
            "retry_or_retune_performed": False,
        })
    except Exception as exc:
        report["status"] = "AUDIT_INCOMPLETE_FAIL_PRESERVED_NO_RETRY"
        report["fatal_error"] = repr(exc)
        report["checks"] = checks
        report["hard_failures"] = hard_failures + [repr(exc)]
        report["limitations"] = limitations
        report["retry_or_retune_performed"] = False
    report_path = DATA / "native_initialization_audit.json"
    write_json(report_path, report)
    print(json.dumps({"status": report["status"], "hard_failures": report.get("hard_failures", []),
                      "limitations": report.get("limitations", []), "report": str(report_path)}, separators=(",", ":")))
    return 0 if report["status"] == "PASS_NATIVE_INITIALIZATION_ONLY" else 2


if __name__ == "__main__":
    raise SystemExit(main())
