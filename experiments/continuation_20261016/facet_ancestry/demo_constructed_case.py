"""Constructed two-tetrahedron validation of generic facet ancestry only."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from pgfrag.facet_ancestry import build_source_facet_table, event_delta, map_native_cohesive_ancestry


def main() -> None:
    source_nodes = np.asarray(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ],
        dtype=np.float64,
    )
    source_tets = np.asarray([[0, 1, 2, 3], [1, 0, 2, 4]], dtype=np.int64)
    node_ids = [101, 17, 83, 44, 129]
    cell_ids = [9001, 9017]
    source_payload = {
        "case": "constructed two T4s sharing one triangle",
        "nodes": source_nodes.tolist(),
        "tets": source_tets.tolist(),
        "source_node_ids": node_ids,
        "source_cell_ids": cell_ids,
    }
    anchor = hashlib.sha256(
        json.dumps(source_payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    source = build_source_facet_table(source_nodes, source_tets, node_ids, cell_ids, anchor)
    if len(source.internal_facets) != 1:
        raise RuntimeError("constructed source did not produce its expected single internal facet")
    facet = source.internal_facets[0]

    # An empty native cohesive table is a valid partial-insertion snapshot.
    empty_snapshot = map_native_cohesive_ancestry(
        source,
        source_nodes,
        source_tets,
        np.empty((0, 6), dtype=np.int64),
        native_initial_bulk_tets=source_tets,
    )

    # Duplicate the shared source nodes in the second current bulk cell.  Each
    # duplicate inherits source ancestry solely through row/local-corner
    # lineage and must retain the exact corresponding source rest coordinate.
    row_by_source_id = {source_id: row for row, source_id in enumerate(node_ids)}
    native_nodes = [row.copy() for row in source_nodes]
    duplicate_by_source_id: dict[int, int] = {}
    for source_id in facet.source_node_ids:
        duplicate_by_source_id[source_id] = len(native_nodes)
        native_nodes.append(source_nodes[row_by_source_id[source_id]].copy())
    native_nodes_array = np.asarray(native_nodes, dtype=np.float64)

    native_bulk = np.array(source_tets, copy=True)
    for local_ordinal, source_node_row in enumerate(source_tets[1]):
        source_id = node_ids[int(source_node_row)]
        if source_id in duplicate_by_source_id:
            native_bulk[1, local_ordinal] = duplicate_by_source_id[source_id]
    side_a = [0, 2, 1]
    side_b = [duplicate_by_source_id[source_id] for source_id in (node_ids[0], node_ids[2], node_ids[1])]
    native_cohesive = np.asarray([side_a + side_b], dtype=np.int64)
    inserted_snapshot = map_native_cohesive_ancestry(
        source,
        native_nodes_array,
        native_bulk,
        native_cohesive,
        native_initial_bulk_tets=source_tets,
    )
    delta = event_delta(
        empty_snapshot["active_facet_identities"],
        inserted_snapshot["active_facet_identities"],
    )

    receipt = {
        "case_status": "CONSTRUCTED_ALGORITHM_VALIDATION_ONLY",
        "case_label": source_payload["case"],
        "source_anchor_hash": anchor,
        "source_tetrahedron_count": int(len(source_tets)),
        "source_facet_count": int(len(source.facets_by_identity)),
        "eligible_internal_facet_identities": list(source.eligible_facet_identities()),
        "empty_snapshot_active_facet_identities": empty_snapshot["active_facet_identities"],
        "empty_snapshot_missing_eligible_facet_identities": empty_snapshot[
            "missing_eligible_facet_identities"
        ],
        "inserted_snapshot_active_facet_identities": inserted_snapshot["active_facet_identities"],
        "inserted_snapshot_missing_eligible_facet_identities": inserted_snapshot[
            "missing_eligible_facet_identities"
        ],
        "event_delta": delta,
        "native_node_count_after_duplication": int(len(native_nodes_array)),
        "native_node_to_source_node_id": inserted_snapshot["native_node_to_source_node_id"],
        "interpretation": inserted_snapshot["interpretation"],
        "solver_run": False,
        "physical_fracture_claim": False,
    }
    output = Path(__file__).with_name("constructed_case_receipt.json")
    output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
