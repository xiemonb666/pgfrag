"""Read-only ancestry check on the saved 96-T4 cohesive-event-probe snapshots.

This script loads historical MSH/NPZ artifacts only.  It does not start PG,
invoke a solver, rewrite source evidence, or modify the saved input snapshots.
"""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np

from pgfrag.cohesive_mesh import read_gmsh22
from pgfrag.facet_ancestry import (
    build_source_facet_table,
    event_delta,
    map_native_cohesive_ancestry,
)


WORKSPACE_ROOT = Path(__file__).resolve().parents[4]
CASE_ROOT = (
    WORKSPACE_ROOT
    / "relicweave-synlab"
    / "tmp"
    / "review_20261004"
    / "implementation"
    / "cohesive_event_probe_20261004"
    / "coarse654_20261004_01"
)
OUTPUT_DIR = Path(__file__).resolve().parent
MSH_PATH = CASE_ROOT / "snapshot" / "input" / "block_2x2x4mm.msh"
INITIAL_PATH = CASE_ROOT / "input_native_initial.npz"
DENSE_PATHS = {
    step: CASE_ROOT / "dense_states" / f"step_00650_{step}.npz"
    for step in ("B", "C", "D")
}


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_npz(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as archive:
        return {name: archive[name].copy() for name in archive.files}


def _relative(path: Path) -> str:
    return path.resolve().relative_to(WORKSPACE_ROOT).as_posix()


def _hash_record(path: Path) -> dict[str, object]:
    return {
        "workspace_relative_path": _relative(path),
        "byte_count": path.stat().st_size,
        "sha256": _sha256_file(path),
    }


def _input_hash_inventory() -> dict[str, object]:
    return {
        "historical_inputs": {
            "source_msh": _hash_record(MSH_PATH),
            "input_native_initial_npz": _hash_record(INITIAL_PATH),
            **{f"dense_step_650_{step}": _hash_record(path) for step, path in DENSE_PATHS.items()},
        },
        "analysis_code": {
            "facet_ancestry_py": _hash_record(WORKSPACE_ROOT / "pgfrag" / "src" / "pgfrag" / "facet_ancestry.py"),
            "cohesive_mesh_reader_py": _hash_record(WORKSPACE_ROOT / "pgfrag" / "src" / "pgfrag" / "cohesive_mesh.py"),
            "read_only_runner_py": _hash_record(Path(__file__)),
        },
    }


def main() -> None:
    for path in [MSH_PATH, INITIAL_PATH, *DENSE_PATHS.values()]:
        if not path.is_file():
            raise FileNotFoundError(f"required historical input is absent: {_relative(path)}")
    input_hashes_before = _input_hash_inventory()

    mesh = read_gmsh22(MSH_PATH)
    source_tets = np.asarray([tet["nodes"] for tet in mesh["tetrahedra"]], dtype=np.int64)
    source_node_ids = [int(node_id) for node_id in mesh["node_tags"]]
    source_cell_ids = [int(tet["element_id"]) for tet in mesh["tetrahedra"]]
    source_interface_triangles = [
        tuple(int(node_id) for node_id in triangle["node_tags"])
        for triangle in mesh["triangles"]
        if triangle["physical_name"] == "interface"
    ]
    source_hash = _sha256_file(MSH_PATH)

    initial = _load_npz(INITIAL_PATH)
    required_initial = {"rest_nodes", "bulk_tetra_connectivity", "cohesive_connectivity"}
    if not required_initial.issubset(initial):
        raise ValueError(f"initial native snapshot lacks required arrays: {sorted(required_initial - set(initial))}")
    if not np.array_equal(initial["bulk_tetra_connectivity"], source_tets):
        raise ValueError(
            "saved native initial bulk connectivity does not exactly equal parsed source MSH row/local-corner connectivity; "
            "no remapping or guessing was attempted"
        )
    if not np.array_equal(initial["rest_nodes"], mesh["nodes"]):
        raise ValueError(
            "saved native initial rest coordinates do not exactly equal parsed source MSH coordinates; "
            "no coordinate-based remapping was attempted"
        )
    if initial["cohesive_connectivity"].shape != (0, 6):
        raise ValueError("saved native initial snapshot is not an empty (0, 6) cohesive state")

    source = build_source_facet_table(
        mesh["nodes"],
        source_tets,
        source_node_ids,
        source_cell_ids,
        source_hash,
        eligible_internal_facets=source_interface_triangles,
    )
    if len(source_interface_triangles) != 8 or len(source.eligible_internal_facets) != 8:
        raise ValueError(
            "the saved source MSH does not provide exactly eight eligible internal interface Tri3 facets"
        )

    snapshots: dict[str, dict[str, np.ndarray]] = {"initial": initial}
    reports: dict[str, dict[str, object]] = {}
    for step, path in DENSE_PATHS.items():
        snapshot = _load_npz(path)
        required = {"rest_nodes", "bulk_tetra_connectivity", "cohesive_connectivity"}
        if not required.issubset(snapshot):
            raise ValueError(f"step 650 {step} snapshot lacks required arrays: {sorted(required - set(snapshot))}")
        snapshots[step] = snapshot
        reports[step] = map_native_cohesive_ancestry(
            source,
            snapshot["rest_nodes"],
            snapshot["bulk_tetra_connectivity"],
            snapshot["cohesive_connectivity"],
            native_initial_bulk_tets=initial["bulk_tetra_connectivity"],
        )

    if reports["B"]["active_facet_identities"] != []:
        raise ValueError("saved step 650 B should be empty cohesive connectivity")
    if len(reports["B"]["missing_eligible_facet_identities"]) != 8:
        raise ValueError("saved step 650 B did not report all eight eligible facets as missing")
    for step in ("C", "D"):
        if len(reports[step]["active_facet_identities"]) != 8:
            raise ValueError(f"saved step 650 {step} did not map all eight active cohesive facets")
        if reports[step]["missing_eligible_facet_identities"] != []:
            raise ValueError(f"saved step 650 {step} still has eligible facets missing")
        if snapshots[step]["cohesive_connectivity"].shape != (8, 6):
            raise ValueError(f"saved step 650 {step} does not contain eight six-node cohesive rows")
        for facet_record in reports[step]["active_facets"]:
            known_parents = set(facet_record["source_parent_cell_ids"])
            side_parents = facet_record["side_parent_source_cell_ids"]
            if len(known_parents) != 2 or set(side_parents) != known_parents or len(side_parents) != 2:
                raise ValueError(f"step 650 {step} has an incomplete source parent pair")
            if not set(facet_record["source_parent_cell_rows"]).issubset(range(len(source_cell_ids))):
                raise ValueError(f"step 650 {step} has an unknown source parent row")

    c_ids = reports["C"]["active_facet_identities"]
    d_ids = reports["D"]["active_facet_identities"]
    if c_ids != d_ids:
        raise ValueError("step 650 C and D stable active source facet identity sets differ")
    delta_bc = event_delta(reports["B"]["active_facet_identities"], c_ids)
    delta_bd = event_delta(reports["B"]["active_facet_identities"], d_ids)

    # The stored C and D row orders happen to agree. Reverse D's cohesive rows
    # in memory as a targeted row-order invariance check; no saved input is
    # edited and the helper still identifies the same eight stable facets.
    permuted_d_cohesive = snapshots["D"]["cohesive_connectivity"][::-1].copy()
    permuted_d_report = map_native_cohesive_ancestry(
        source,
        snapshots["D"]["rest_nodes"],
        snapshots["D"]["bulk_tetra_connectivity"],
        permuted_d_cohesive,
        native_initial_bulk_tets=initial["bulk_tetra_connectivity"],
    )
    if permuted_d_report["active_facet_identities"] != d_ids:
        raise ValueError("reversing the saved D cohesive row order changed stable source facet identities")
    original_d_rows = {
        record["facet_identity"]: record["native_cohesive_row"]
        for record in reports["D"]["active_facets"]
    }
    permuted_d_rows = {
        record["facet_identity"]: record["native_cohesive_row"]
        for record in permuted_d_report["active_facets"]
    }
    expected_permuted_rows = {
        identity: len(original_d_rows) - 1 - row for identity, row in original_d_rows.items()
    }
    if permuted_d_rows != expected_permuted_rows:
        raise ValueError("in-memory D cohesive-row reversal did not produce the expected reversed row mapping")

    input_hashes_after = _input_hash_inventory()
    if input_hashes_after != input_hashes_before:
        raise RuntimeError("historical inputs or analysis source changed during the read-only ancestry check")
    output_input_hashes = {
        "before_read_only_analysis": input_hashes_before,
        "after_read_only_analysis": input_hashes_after,
        "unchanged_during_analysis": True,
    }

    active_summaries = {
        step: [
            {
                "facet_identity": record["facet_identity"],
                "source_node_ids": record["source_node_ids"],
                "source_parent_cell_ids": record["source_parent_cell_ids"],
                "source_parent_cell_rows": record["source_parent_cell_rows"],
                "source_parent_local_face_ordinals": record["source_parent_local_face_ordinals"],
                "side_parent_source_cell_ids": record["side_parent_source_cell_ids"],
                "side_parent_source_cell_rows": record["side_parent_source_cell_rows"],
                "native_cohesive_row": record["native_cohesive_row"],
            }
            for record in sorted(
                reports[step]["active_facets"], key=lambda item: item["native_cohesive_row"]
            )
        ]
        for step in ("C", "D")
    }
    receipt = {
        "status": "READ_ONLY_SAVED_NATIVE_ANCESTRY_DEMONSTRATION",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "case_label": "historical 96-T4 coarse coupon; step 650 B/C/D saved dense states",
        "historical_case_root": _relative(CASE_ROOT),
        "input_hashes_file": "input_hashes.json",
        "source_anchor_hash": source_hash,
        "source_mesh_counts": {
            "source_node_count": int(len(mesh["nodes"])),
            "source_tetrahedron_count": int(len(source_tets)),
            "source_cell_id_count": int(len(set(source_cell_ids))),
            "interface_tri3_count": len(source_interface_triangles),
            "eligible_internal_facet_count": len(source.eligible_internal_facets),
            "eligible_internal_facets_with_two_parents": sum(
                len(source.facets_by_identity[identity].parents) == 2
                for identity in source.eligible_internal_facets
            ),
        },
        "source_id_basis": {
            "source_node_ids": "original Gmsh MSH 2.2 node tags",
            "source_cell_ids": "original Gmsh type-4 tetrahedron element IDs",
        },
        "native_initial_contract": {
            "exact_native_initial_source_connectivity": True,
            "exact_native_initial_source_rest_coordinates": True,
            "empty_native_initial_cohesive_shape": list(initial["cohesive_connectivity"].shape),
        },
        "snapshots": {
            step: {
                "native_rest_node_count": int(len(snapshots[step]["rest_nodes"])),
                "native_bulk_tetrahedron_count": int(len(snapshots[step]["bulk_tetra_connectivity"])),
                "native_cohesive_shape": list(snapshots[step]["cohesive_connectivity"].shape),
                "active_facet_identities": reports[step]["active_facet_identities"],
                "missing_eligible_facet_identities": reports[step]["missing_eligible_facet_identities"],
                "lineage_scope": reports[step]["lineage_scope"],
                "lineage_assumptions": reports[step]["lineage_assumptions"],
                "qp_ordinal_stability_claimed": reports[step]["qp_ordinal_stability_claimed"],
            }
            for step in ("B", "C", "D")
        },
        "event_deltas_from_empty_B": {"B_to_C": delta_bc, "B_to_D": delta_bd},
        "active_facet_parent_pairs": active_summaries,
        "cohesive_row_order_check": {
            "saved_C_and_D_identity_sets_equal": True,
            "saved_C_and_D_native_cohesive_row_orders_equal": [
                record["facet_identity"] for record in active_summaries["C"]
            ]
            == [record["facet_identity"] for record in active_summaries["D"]],
            "in_memory_reversed_D_rows_preserve_all_eight_stable_identities": True,
            "in_memory_reversed_D_row_identity_mapping_verified": True,
            "in_memory_row_permutation": list(
                reversed(range(len(snapshots["D"]["cohesive_connectivity"])))
            ),
        },
        "source_parent_and_side_pair_validation": "all eight saved C/D cohesive facets map to two known original MSH parent cell IDs, one parent incidence per side",
        "solver_invoked": False,
        "saved_historical_evidence_modified": False,
        "physical_fracture_claim": False,
        "material_or_qp_state_claim": False,
    }

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUTPUT_DIR / "input_hashes.json").write_text(
        json.dumps(output_input_hashes, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (OUTPUT_DIR / "saved_96t4_run_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "status": receipt["status"],
                "source_anchor_hash": source_hash,
                "source_mesh_counts": receipt["source_mesh_counts"],
                "snapshot_active_facet_counts": {
                    step: len(reports[step]["active_facet_identities"]) for step in ("B", "C", "D")
                },
                "saved_C_D_row_orders_equal": receipt["cohesive_row_order_check"][
                    "saved_C_and_D_native_cohesive_row_orders_equal"
                ],
                "reversed_D_rows_keep_stable_ids": receipt["cohesive_row_order_check"][
                    "in_memory_reversed_D_rows_preserve_all_eight_stable_identities"
                ],
                "historical_inputs_unchanged": output_input_hashes["unchanged_during_analysis"],
                "receipt": (OUTPUT_DIR / "saved_96t4_run_receipt.json").name,
                "input_hashes": (OUTPUT_DIR / "input_hashes.json").name,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
