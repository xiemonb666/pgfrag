from __future__ import annotations

import hashlib

import numpy as np
import pytest

from pgfrag.facet_ancestry import (
    build_source_facet_table,
    event_delta,
    make_facet_identity,
    map_native_cohesive_ancestry,
)


def _anchor(label: str = "facet-ancestry-test") -> str:
    return hashlib.sha256(label.encode("utf-8")).hexdigest()


def _positive_tet(nodes: np.ndarray, tet: list[int]) -> list[int]:
    points = nodes[np.asarray(tet, dtype=np.int64)]
    volume6 = float(np.dot(points[1] - points[0], np.cross(points[2] - points[0], points[3] - points[0])))
    assert volume6 != 0.0
    return tet if volume6 > 0.0 else [tet[1], tet[0], tet[2], tet[3]]


def _two_tet_source(*, anchor: str | None = None, cell_offset: int = 0):
    nodes = np.asarray(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ],
        dtype=np.float64,
    )
    tets = np.asarray(
        [
            _positive_tet(nodes, [0, 1, 2, 3]),
            _positive_tet(nodes, [0, 1, 2, 4]),
        ],
        dtype=np.int64,
    )
    source = build_source_facet_table(
        nodes,
        tets,
        [101, 17, 83, 44, 129],
        [9001 + cell_offset, 9017 + cell_offset],
        anchor or _anchor(),
    )
    shared = next(facet for facet in source.internal_facets)
    return source, shared


def _split_facet_snapshot(source, facet):
    """Build a node-duplicated, one-cohesive-row snapshot for a source facet."""
    native_nodes = [row.copy() for row in source.source_nodes]
    source_row_by_id = {source_id: row for row, source_id in enumerate(source.source_node_ids)}
    duplicate_by_source_id: dict[int, int] = {}
    for source_id in facet.source_node_ids:
        duplicate_by_source_id[source_id] = len(native_nodes)
        native_nodes.append(source.source_nodes[source_row_by_id[source_id]].copy())

    native_bulk = np.array(source.source_tets, dtype=np.int64, copy=True)
    parent_rows = [parent.source_cell_row for parent in facet.parents]
    split_row = parent_rows[1]
    for local_ordinal, source_node_row in enumerate(source.source_tets[split_row]):
        source_id = source.source_node_ids[int(source_node_row)]
        if source_id in duplicate_by_source_id:
            native_bulk[split_row, local_ordinal] = duplicate_by_source_id[source_id]

    side_nodes: list[list[int]] = []
    for parent_index, parent_row in enumerate(parent_rows):
        side: list[int] = []
        for source_node_row in source.source_tets[parent_row]:
            source_id = source.source_node_ids[int(source_node_row)]
            if source_id in facet.source_node_ids:
                if parent_index == 1:
                    side.append(duplicate_by_source_id[source_id])
                else:
                    side.append(int(source_node_row))
        side_nodes.append(side)
    # Deliberately vary side order to exercise ID-set pairing rather than
    # assuming that native cohesive side corner order is identical.
    side_nodes[0] = list(reversed(side_nodes[0]))
    cohesive = np.asarray([side_nodes[0] + side_nodes[1]], dtype=np.int64)
    return np.asarray(native_nodes), native_bulk, cohesive


def test_one_tetrahedron_builds_four_exterior_facets_and_no_eligible_internal_facets():
    nodes = np.asarray([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    tets = np.asarray([_positive_tet(nodes, [0, 1, 2, 3])], dtype=np.int64)
    table = build_source_facet_table(nodes, tets, [7, 9, 11, 13], [41], _anchor())

    assert len(table.facets_by_identity) == 4
    assert table.internal_facets == ()
    assert table.eligible_internal_facets == frozenset()
    assert {facet.classification for facet in table.facets_by_identity.values()} == {"exterior"}


def test_two_positive_tetrahedra_have_opposed_parent_windings_and_stable_identity():
    table, shared = _two_tet_source()

    assert shared.classification == "internal"
    assert len(shared.parents) == 2
    assert {parent.source_cell_id for parent in shared.parents} == {9001, 9017}
    assert shared.source_node_ids == tuple(sorted(shared.source_node_ids))
    assert shared.identity == make_facet_identity(_anchor(), shared.source_node_ids)
    parities = []
    for parent in shared.parents:
        ids = parent.outward_ordered_source_corner_ids
        parities.append(sum(ids[a] > ids[b] for a in range(3) for b in range(a + 1, 3)) % 2)
    assert parities[0] != parities[1]
    assert table.eligible_internal_facets == frozenset({shared.identity})


def test_duplicate_coordinates_with_different_source_ids_remain_distinct_facets():
    one = np.asarray([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    nodes = np.vstack([one, one])
    tets = np.asarray([[0, 1, 2, 3], [4, 5, 6, 7]], dtype=np.int64)
    table = build_source_facet_table(
        nodes,
        tets,
        [1, 2, 3, 4, 101, 102, 103, 104],
        [12, 13],
        _anchor("coincident-nodes"),
    )

    first = make_facet_identity(_anchor("coincident-nodes"), [1, 2, 3])
    second = make_facet_identity(_anchor("coincident-nodes"), [101, 102, 103])
    assert first != second
    assert first in table.facets_by_identity
    assert second in table.facets_by_identity
    assert table.facets_by_identity[first].classification == "exterior"
    assert table.facets_by_identity[second].classification == "exterior"


def test_eligible_facets_can_be_explicitly_empty_or_selected_by_identity_or_source_ids():
    source, shared = _two_tet_source()
    empty = build_source_facet_table(
        source.source_nodes,
        source.source_tets,
        source.source_node_ids,
        source.source_cell_ids,
        source.source_anchor_hash,
        eligible_internal_facets=[],
    )
    selected_by_identity = build_source_facet_table(
        source.source_nodes,
        source.source_tets,
        source.source_node_ids,
        source.source_cell_ids,
        source.source_anchor_hash,
        eligible_internal_facets=[shared.identity],
    )
    selected_by_ids = build_source_facet_table(
        source.source_nodes,
        source.source_tets,
        source.source_node_ids,
        source.source_cell_ids,
        source.source_anchor_hash,
        eligible_internal_facets=[shared.source_node_ids],
    )

    assert empty.eligible_internal_facets == frozenset()
    assert selected_by_identity.eligible_internal_facets == frozenset({shared.identity})
    assert selected_by_ids.eligible_internal_facets == frozenset({shared.identity})


def test_exterior_or_unknown_facet_cannot_be_marked_eligible():
    source, shared = _two_tet_source()
    exterior = next(facet for facet in source.facets_by_identity.values() if facet.classification == "exterior")
    with pytest.raises(ValueError, match="exterior, not internal"):
        build_source_facet_table(
            source.source_nodes,
            source.source_tets,
            source.source_node_ids,
            source.source_cell_ids,
            source.source_anchor_hash,
            eligible_internal_facets=[exterior.source_node_ids],
        )
    with pytest.raises(ValueError, match="absent from the source mesh"):
        build_source_facet_table(
            source.source_nodes,
            source.source_tets,
            source.source_node_ids,
            source.source_cell_ids,
            source.source_anchor_hash,
            eligible_internal_facets=[(9991, 9992, 9993)],
        )
    assert shared.identity in source.eligible_internal_facets


def test_degenerate_inverted_nonmanifold_and_same_side_parents_fail_closed():
    nodes = np.asarray([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    with pytest.raises(ValueError, match="positive rest orientation"):
        build_source_facet_table(nodes, np.asarray([[0, 2, 1, 3]]), [1, 2, 3, 4], [5], _anchor())
    with pytest.raises(ValueError, match="repeated node row"):
        build_source_facet_table(nodes, np.asarray([[0, 1, 2, 2]]), [1, 2, 3, 4], [5], _anchor())

    nonmanifold_nodes = np.asarray(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, -1], [0, 0, 2]],
        dtype=float,
    )
    nonmanifold_tets = np.asarray(
        [
            _positive_tet(nonmanifold_nodes, [0, 1, 2, 3]),
            _positive_tet(nonmanifold_nodes, [0, 1, 2, 4]),
            _positive_tet(nonmanifold_nodes, [0, 1, 2, 5]),
        ],
        dtype=np.int64,
    )
    with pytest.raises(ValueError, match="nonmanifold source facet"):
        build_source_facet_table(nonmanifold_nodes, nonmanifold_tets, list(range(6)), [10, 11, 12], _anchor())

    same_side_nodes = np.asarray(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [0.25, 0.25, 2]], dtype=float
    )
    same_side_tets = np.asarray(
        [
            _positive_tet(same_side_nodes, [0, 1, 2, 3]),
            _positive_tet(same_side_nodes, [0, 1, 2, 4]),
        ],
        dtype=np.int64,
    )
    with pytest.raises(ValueError, match="opposed outward windings"):
        build_source_facet_table(same_side_nodes, same_side_tets, list(range(5)), [20, 21], _anchor())


def test_source_ids_and_source_anchor_are_validated():
    source, _ = _two_tet_source()
    with pytest.raises(ValueError, match="source node IDs must be unique"):
        build_source_facet_table(
            source.source_nodes,
            source.source_tets,
            [1, 1, 3, 4, 5],
            source.source_cell_ids,
            source.source_anchor_hash,
        )
    with pytest.raises(ValueError, match="source cell IDs must be unique"):
        build_source_facet_table(
            source.source_nodes,
            source.source_tets,
            source.source_node_ids,
            [1, 1],
            source.source_anchor_hash,
        )
    with pytest.raises(ValueError, match="64-character SHA-256"):
        build_source_facet_table(
            source.source_nodes,
            source.source_tets,
            source.source_node_ids,
            source.source_cell_ids,
            "not-frozen",
        )


def test_empty_cohesive_snapshot_is_valid_and_reports_all_eligible_facets_missing():
    source, shared = _two_tet_source()
    report = map_native_cohesive_ancestry(
        source,
        source.source_nodes,
        source.source_tets,
        np.empty((0, 6), dtype=np.int64),
        native_initial_bulk_tets=source.source_tets,
    )

    assert report["active_facet_identities"] == []
    assert report["missing_eligible_facet_identities"] == [shared.identity]
    assert report["active_cohesive_count"] == 0
    assert report["interpretation"].startswith("connectivity ancestry only")


def test_node_duplicated_native_pair_maps_to_one_eligible_source_facet():
    source, shared = _two_tet_source()
    rest, bulk, cohesive = _split_facet_snapshot(source, shared)
    report = map_native_cohesive_ancestry(
        source,
        rest,
        bulk,
        cohesive,
        native_initial_bulk_tets=source.source_tets,
    )

    assert report["active_facet_identities"] == [shared.identity]
    assert report["missing_eligible_facet_identities"] == []
    assert report["active_facets"][0]["source_parent_cell_ids"] == [9001, 9017]
    assert report["active_facets"][0]["side_parent_source_cell_ids"] == [9001, 9017]
    assert len(report["native_node_to_source_node_id"]) == len(rest)


def test_native_cohesive_row_reordering_does_not_change_stable_facet_identities():
    first_nodes = np.asarray(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, -1]], dtype=float
    )
    second_nodes = first_nodes + np.asarray([3.0, 0.0, 0.0])
    nodes = np.vstack([first_nodes, second_nodes])
    first_tets = np.asarray(
        [_positive_tet(first_nodes, [0, 1, 2, 3]), _positive_tet(first_nodes, [0, 1, 2, 4])],
        dtype=np.int64,
    )
    second_tets = np.asarray(
        [_positive_tet(second_nodes, [0, 1, 2, 3]), _positive_tet(second_nodes, [0, 1, 2, 4])],
        dtype=np.int64,
    ) + 5
    source_tets = np.vstack([first_tets, second_tets])
    source = build_source_facet_table(
        nodes,
        source_tets,
        [10, 20, 30, 40, 50, 110, 120, 130, 140, 150],
        [201, 202, 301, 302],
        _anchor("two-facets"),
    )
    facets = list(source.internal_facets)
    assert len(facets) == 2

    native_nodes = [row.copy() for row in source.source_nodes]
    native_bulk = np.array(source.source_tets, copy=True)
    cohesive_rows = []
    for facet in facets:
        duplicates: dict[int, int] = {}
        for source_id in facet.source_node_ids:
            source_node_row = source.source_node_ids.index(source_id)
            duplicates[source_id] = len(native_nodes)
            native_nodes.append(source.source_nodes[source_node_row].copy())
        row_a, row_b = [parent.source_cell_row for parent in facet.parents]
        for local_ordinal, source_node_row in enumerate(source.source_tets[row_b]):
            source_id = source.source_node_ids[int(source_node_row)]
            if source_id in duplicates:
                native_bulk[row_b, local_ordinal] = duplicates[source_id]
        sides = []
        for parent_index, parent_row in enumerate((row_a, row_b)):
            side = []
            for source_node_row in source.source_tets[parent_row]:
                source_id = source.source_node_ids[int(source_node_row)]
                if source_id in facet.source_node_ids:
                    side.append(duplicates[source_id] if parent_index else int(source_node_row))
            sides.append(side)
        cohesive_rows.append(sides[0] + sides[1])

    rest = np.asarray(native_nodes)
    cohesive = np.asarray(cohesive_rows, dtype=np.int64)
    first_report = map_native_cohesive_ancestry(
        source, rest, native_bulk, cohesive, native_initial_bulk_tets=source.source_tets
    )
    reordered_report = map_native_cohesive_ancestry(
        source, rest, native_bulk, cohesive[::-1], native_initial_bulk_tets=source.source_tets
    )

    assert first_report["active_facet_identities"] == reordered_report["active_facet_identities"]
    first_by_identity = {row["facet_identity"]: dict(row) for row in first_report["active_facets"]}
    reordered_by_identity = {row["facet_identity"]: dict(row) for row in reordered_report["active_facets"]}
    for rows in (first_by_identity, reordered_by_identity):
        for row in rows.values():
            row.pop("native_cohesive_row")
    assert first_by_identity == reordered_by_identity


def test_native_initial_and_current_bulk_lineage_must_remain_exact():
    source, shared = _two_tet_source()
    rest, bulk, cohesive = _split_facet_snapshot(source, shared)
    changed_initial = np.array(source.source_tets, copy=True)
    changed_initial[[0, 1]] = changed_initial[[1, 0]]
    with pytest.raises(ValueError, match="initial bulk connectivity differs"):
        map_native_cohesive_ancestry(
            source, rest, bulk, cohesive, native_initial_bulk_tets=changed_initial
        )

    reordered_current = np.array(bulk[[1, 0]], copy=True)
    with pytest.raises(ValueError, match="exact rest coordinate consistency"):
        map_native_cohesive_ancestry(
            source, rest, reordered_current, cohesive, native_initial_bulk_tets=source.source_tets
        )


def test_exact_rest_coordinates_are_required_for_every_current_bulk_corner():
    source, shared = _two_tet_source()
    rest, bulk, cohesive = _split_facet_snapshot(source, shared)
    rest[5, 0] = np.nextafter(rest[5, 0], np.inf)

    with pytest.raises(ValueError, match="exact rest coordinate consistency"):
        map_native_cohesive_ancestry(
            source, rest, bulk, cohesive, native_initial_bulk_tets=source.source_tets
        )


def test_unmapped_native_nodes_conflicting_corner_ancestry_and_bad_connectivity_fail():
    source, shared = _two_tet_source()
    rest, bulk, cohesive = _split_facet_snapshot(source, shared)
    with_extra = np.vstack([rest, [9.0, 9.0, 9.0]])
    with pytest.raises(ValueError, match="lack source bulk-corner ancestry"):
        map_native_cohesive_ancestry(
            source, with_extra, bulk, cohesive, native_initial_bulk_tets=source.source_tets
        )

    conflict_bulk = np.array(bulk, copy=True)
    # A shared native node cannot stand for two distinct source IDs.
    conflict_bulk[1, 0] = conflict_bulk[0, 0]
    with pytest.raises(ValueError, match="conflicting source ancestry"):
        map_native_cohesive_ancestry(
            source, rest, conflict_bulk, cohesive, native_initial_bulk_tets=source.source_tets
        )

    with pytest.raises(ValueError, match=r"shape \(n, 6\)"):
        map_native_cohesive_ancestry(
            source, rest, bulk, np.asarray([[0, 1, 2, 3, 4]], dtype=np.int64),
            native_initial_bulk_tets=source.source_tets,
        )


def test_duplicate_native_cohesive_birth_for_one_source_facet_is_rejected():
    source, shared = _two_tet_source()
    rest, bulk, cohesive = _split_facet_snapshot(source, shared)
    duplicate_rows = np.vstack([cohesive, cohesive])
    with pytest.raises(ValueError, match="duplicate native cohesive births"):
        map_native_cohesive_ancestry(
            source, rest, bulk, duplicate_rows, native_initial_bulk_tets=source.source_tets
        )


def test_cohesive_facet_must_be_explicitly_eligible_and_have_one_current_parent_per_side():
    source, shared = _two_tet_source()
    rest, bulk, cohesive = _split_facet_snapshot(source, shared)
    ineligible = build_source_facet_table(
        source.source_nodes,
        source.source_tets,
        source.source_node_ids,
        source.source_cell_ids,
        source.source_anchor_hash,
        eligible_internal_facets=[],
    )
    with pytest.raises(ValueError, match="ineligible source facet"):
        map_native_cohesive_ancestry(
            ineligible, rest, bulk, cohesive, native_initial_bulk_tets=source.source_tets
        )

    unsplit_bulk = np.array(source.source_tets, copy=True)
    with pytest.raises(ValueError, match="exactly one"):
        map_native_cohesive_ancestry(
            source,
            source.source_nodes,
            unsplit_bulk,
            np.asarray([[0, 1, 2, 0, 1, 2]], dtype=np.int64),
            native_initial_bulk_tets=source.source_tets,
        )


def test_event_delta_allows_monotone_partial_insertion_and_rejects_deletion_reuse_or_mixed_anchor():
    first_source, first_facet = _two_tet_source(anchor=_anchor("event-a"))
    second_source, second_facet = _two_tet_source(anchor=_anchor("event-b"))
    another_first_anchor_facet = make_facet_identity(first_source.source_anchor_hash, [701, 702, 703])
    empty = event_delta([], [first_facet.identity])
    partial = event_delta(
        [first_facet.identity], [first_facet.identity, another_first_anchor_facet]
    )
    persistent = event_delta([first_facet.identity], [first_facet.identity])

    assert empty["newly_active_facet_identities"] == [first_facet.identity]
    assert partial["newly_active_facet_identities"] == [another_first_anchor_facet]
    assert persistent["persistent_active_facet_identities"] == [first_facet.identity]
    with pytest.raises(ValueError, match="disappeared"):
        event_delta([first_facet.identity], [])
    with pytest.raises(ValueError, match="reused/duplicate"):
        event_delta([], [first_facet.identity, first_facet.identity])
    with pytest.raises(ValueError, match="different source-anchor hashes"):
        event_delta([first_facet.identity], [second_facet.identity])
    assert first_source.source_anchor_hash != second_source.source_anchor_hash

