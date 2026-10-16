from __future__ import annotations

import numpy as np
import pytest

from pgfrag.cohesive_fragments import extract_cohesive_fragments


SOURCE_EXTERIOR = {
    10: [1, 2, 3],
    11: [0, 2, 3],
    12: [0, 1, 3],
    13: [0, 2, 4],
    14: [0, 1, 4],
    15: [1, 2, 4],
}


def two_tet_patch() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    # Two positively oriented tetrahedra share a geometric face represented
    # with separate native node IDs, as after cohesive insertion.
    rest = np.asarray([
        [0.0, 0.0, 0.0],  # source 0, upper side
        [1.0, 0.0, 0.0],  # source 1, upper side
        [0.0, 1.0, 0.0],  # source 2, upper side
        [0.0, 0.0, 1.0],  # upper apex
        [0.0, 0.0, 0.0],  # source 0, lower side
        [1.0, 0.0, 0.0],  # source 1, lower side
        [0.0, 1.0, 0.0],  # source 2, lower side
        [0.0, 0.0, -1.0],  # lower apex
    ], dtype=np.float64)
    tets = np.asarray([[0, 1, 2, 3], [4, 6, 5, 7]], dtype=np.int64)
    cohesive = np.asarray([[0, 1, 2, 4, 5, 6]], dtype=np.int64)
    source = np.asarray([0, 1, 2, 3, 0, 1, 2, 4], dtype=np.int64)
    return rest, tets, cohesive, source


def empty_cohesive() -> np.ndarray:
    return np.empty((0, 6), dtype=np.int64)


def test_intact_and_fully_failed_cohesive_parent_pairs() -> None:
    rest, tets, cohesive, source = two_tet_patch()
    intact = extract_cohesive_fragments(
        rest, tets, cohesive,
        facet_damage_qp=[[0.2, 0.95]],
        tractions_qp=[[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]],
        damage_tolerance=0.01,
        traction_tolerance=1.0e-8,
        current_nodes=rest.copy(),
        native_node_to_source=source,
        parent_facet_ids=[81],
        parent_source_ids=[9001],
        source_exterior_faces=SOURCE_EXTERIOR,
        rho=12.0,
    )
    assert intact["failed_facet_rows"] == []
    assert len(intact["components"]) == 1
    intact_fragment = intact["components"][0]
    assert not intact_fragment["closed_solid"]
    assert intact_fragment["intact_cohesive_welded_edge_audit"]["closed"]
    assert not intact_fragment["qualified_closed_mechanical_fragment"]
    assert "active_cohesive_interface_not_geometrically_closed" in intact_fragment["qualification_reasons"]
    assert len(intact_fragment["oriented_native_triangles"]) == 6
    assert not intact_fragment["closed_edge_audit"]["every_edge_twice_in_reverse"]
    assert all(face["surface_kind"] == "original_exterior" for face in intact_fragment["boundary_faces"])
    assert all(face["source_exterior_ancestry_verified"] for face in intact_fragment["boundary_faces"])
    assert intact_fragment["rest_bulk_volume"] == pytest.approx(1.0 / 3.0)
    assert intact_fragment["current_bulk_volume"] == pytest.approx(1.0 / 3.0)
    assert intact_fragment["rest_bulk_mass"] == pytest.approx(4.0)
    assert intact["component_rest_volume_sum_matches_total"]
    assert intact["total_rest_bulk_mass"] == pytest.approx(4.0)

    failed = extract_cohesive_fragments(
        rest, tets, cohesive,
        facet_damage_qp=[[1.0, 1.0]],
        tractions_qp=[[[0.0, 0.0, 0.0], [1.0e-10, 0.0, 0.0]]],
        damage_tolerance=0.01,
        traction_tolerance=1.0e-9,
        native_node_to_source=source,
        parent_facet_ids=[81],
        parent_source_ids=[9001],
        source_exterior_faces=SOURCE_EXTERIOR,
        rho=12.0,
    )
    assert failed["failed_facet_rows"] == [0]
    assert len(failed["components"]) == 2
    assert all(fragment["closed_solid"] for fragment in failed["components"])
    assert all(fragment["qualified_closed_mechanical_fragment"] for fragment in failed["components"])
    assert all(len(fragment["oriented_native_triangles"]) == 4 for fragment in failed["components"])
    assert all(fragment["rest_bulk_volume"] == pytest.approx(1.0 / 6.0) for fragment in failed["components"])
    assert all(fragment["rest_bulk_mass"] == pytest.approx(2.0) for fragment in failed["components"])
    pair = failed["failed_interface_face_pairs"][0]
    assert pair["parent_facet_id"] == 81
    assert pair["parent_source_id"] == 9001
    assert pair["side_a_native_node_key"] == [0, 1, 2]
    assert pair["side_b_native_node_key"] == [4, 5, 6]
    assert pair["side_parent_tet_rows"] == [0, 1]
    assert pair["side_a_rest_area"] == pytest.approx(0.5)
    assert pair["side_b_rest_area"] == pytest.approx(0.5)


def test_one_qp_that_does_not_pass_prevents_facet_failure() -> None:
    rest, tets, cohesive, source = two_tet_patch()
    result = extract_cohesive_fragments(
        rest, tets, cohesive,
        facet_damage_qp=[[1.0, 0.98]],
        tractions_qp=[[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]],
        damage_tolerance=0.01,
        traction_tolerance=1.0e-8,
        native_node_to_source=source,
    )
    assert result["facet_failure_evidence"][0]["qp_joint_failure_pass"] == [True, False]
    assert result["failed_facet_rows"] == []
    assert len(result["components"]) == 1

    traction_result = extract_cohesive_fragments(
        rest, tets, cohesive,
        facet_damage_qp=[[1.0, 1.0]],
        tractions_qp=[[[0.0, 0.0, 0.0], [1.0e-3, 0.0, 0.0]]],
        damage_tolerance=0.01,
        traction_tolerance=1.0e-8,
        native_node_to_source=source,
    )
    assert traction_result["failed_facet_rows"] == []
    assert len(traction_result["components"]) == 1


def test_ordinary_shared_native_face_connects_tetrahedra() -> None:
    rest = np.asarray([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0], [0.0, 0.0, -1.0],
    ])
    tets = np.asarray([[0, 1, 2, 3], [0, 2, 1, 4]], dtype=np.int64)
    result = extract_cohesive_fragments(
        rest, tets, empty_cohesive(),
        facet_damage_qp=np.empty((0, 1)),
        tractions_qp=np.empty((0, 1, 3)),
        damage_tolerance=0.01,
        traction_tolerance=1.0e-8,
    )
    assert len(result["components"]) == 1
    assert result["components"][0]["tet_rows"] == [0, 1]
    assert result["components"][0]["closed_solid"]
    assert not result["components"][0]["qualified_closed_mechanical_fragment"]
    assert "source_exterior_ancestry_unverified" in result["components"][0]["qualification_reasons"]
    assert len(result["components"][0]["oriented_native_triangles"]) == 6


def test_paired_rest_coordinates_must_match_for_source_and_ordinal_pairing() -> None:
    rest, tets, cohesive, source = two_tet_patch()
    shifted = rest.copy()
    shifted[4, 0] = 0.01
    with pytest.raises(ValueError, match="paired rest coordinates do not match exactly"):
        extract_cohesive_fragments(
            shifted, tets, cohesive,
            facet_damage_qp=[[1.0]],
            tractions_qp=[[[0.0, 0.0, 0.0]]],
            damage_tolerance=0.01,
            traction_tolerance=1.0e-8,
            native_node_to_source=source,
        )

    permuted = cohesive.copy()
    permuted[0, 3:] = [5, 6, 4]
    with pytest.raises(ValueError, match="paired rest coordinates do not match exactly"):
        extract_cohesive_fragments(
            rest, tets, permuted,
            facet_damage_qp=[[1.0]],
            tractions_qp=[[[0.0, 0.0, 0.0]]],
            damage_tolerance=0.01,
            traction_tolerance=1.0e-8,
        )


def test_shared_native_node_across_components_denies_fragment_qualification() -> None:
    rest = np.asarray([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
        [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0],
    ])
    tets = np.asarray([[0, 1, 2, 3], [0, 5, 4, 6]], dtype=np.int64)
    result = extract_cohesive_fragments(
        rest, tets, empty_cohesive(),
        facet_damage_qp=np.empty((0, 1)),
        tractions_qp=np.empty((0, 1, 3)),
        damage_tolerance=0.01,
        traction_tolerance=1.0e-8,
    )
    assert len(result["components"]) == 2
    assert all(fragment["closed_solid"] for fragment in result["components"])
    assert all(not fragment["qualified_closed_mechanical_fragment"] for fragment in result["components"])
    assert all(fragment["shared_native_node_ids_with_other_components"] == [0] for fragment in result["components"])
    assert all("native_node_shared_with_another_mechanical_component" in fragment["qualification_reasons"]
               for fragment in result["components"])


def test_source_exterior_ancestry_must_cover_every_noncohesive_boundary_face() -> None:
    rest, tets, cohesive, source = two_tet_patch()
    source_exterior = {
        10: [1, 2, 3],
        11: [0, 2, 3],
        12: [0, 1, 3],
        13: [0, 2, 4],
        14: [0, 1, 4],
        # Intentionally omit the sixth outer face [1, 2, 4].
    }
    with pytest.raises(ValueError, match="no matching source exterior face ancestry"):
        extract_cohesive_fragments(
            rest, tets, cohesive,
            facet_damage_qp=[[0.0]],
            tractions_qp=[[[0.0, 0.0, 0.0]]],
            damage_tolerance=0.01,
            traction_tolerance=1.0e-8,
            native_node_to_source=source,
            source_exterior_faces=source_exterior,
        )


def test_wrong_cohesive_side_source_pairing_is_rejected() -> None:
    rest, tets, cohesive, source = two_tet_patch()
    wrong = cohesive.copy()
    wrong[0, 5] = 7  # A lower apex, not a source-corresponding base vertex.
    with pytest.raises(ValueError, match="side 1 has 0 native T4 parents|do not pair by source node ancestry"):
        extract_cohesive_fragments(
            rest, tets, wrong,
            facet_damage_qp=[[1.0]],
            tractions_qp=[[[0.0, 0.0, 0.0]]],
            damage_tolerance=0.01,
            traction_tolerance=1.0e-8,
            native_node_to_source=source,
        )


def test_nonpositive_rest_or_current_volume_is_rejected() -> None:
    rest, tets, cohesive, source = two_tet_patch()
    inverted = rest.copy()
    inverted[[5, 6]] = inverted[[6, 5]]
    with pytest.raises(ValueError, match="current tetrahedron signed volume must be positive"):
        extract_cohesive_fragments(
            rest, tets, cohesive,
            facet_damage_qp=[[0.0]],
            tractions_qp=[[[0.0, 0.0, 0.0]]],
            damage_tolerance=0.01,
            traction_tolerance=1.0e-8,
            current_nodes=inverted,
            native_node_to_source=source,
        )

    bad_tets = tets.copy()
    bad_tets[0, [1, 2]] = bad_tets[0, [2, 1]]
    with pytest.raises(ValueError, match="rest tetrahedron signed volume must be positive"):
        extract_cohesive_fragments(
            rest, bad_tets, cohesive,
            facet_damage_qp=[[0.0]],
            tractions_qp=[[[0.0, 0.0, 0.0]]],
            damage_tolerance=0.01,
            traction_tolerance=1.0e-8,
            native_node_to_source=source,
        )
