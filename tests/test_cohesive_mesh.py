from __future__ import annotations

import numpy as np
import pytest

from pgfrag.cohesive_mesh import (
    build_native_lineage,
    input_interface_contract,
    quadrature_signed_work,
    read_gmsh22,
    triangle_area,
)


def _mesh_data() -> dict:
    # The shared triangle is in the x=0 plane. Both parent centroids have
    # z > median(source z), so a hard-coded z/median split labels them alike.
    nodes = np.asarray([
        [0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [-1.0, 0.2, 0.2],
        [1.0, 0.4, 0.2],
    ])
    tets = [
        {"element_id": 101, "nodes": [3, 0, 1, 2]},
        {"element_id": 202, "nodes": [4, 0, 2, 1]},
    ]
    triangles = [{
        "element_id": 501,
        "element_type": 2,
        "physical_tag": 1,
        "physical_name": "interface",
        "node_tags": [1, 2, 3],
        "nodes": [0, 1, 2],
    }]
    return {"nodes": nodes, "tetrahedra": tets, "triangles": triangles}


def _lineage_inputs() -> tuple[dict, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    mesh_data = _mesh_data()
    native_nodes = np.vstack([mesh_data["nodes"], mesh_data["nodes"][[0, 1, 2]]])
    initial_tets = np.asarray([row["nodes"] for row in mesh_data["tetrahedra"]])
    post_tets = np.asarray([
        [3, 0, 1, 2],
        [4, 5, 7, 6],
    ])
    cohesive = np.asarray([[0, 1, 2, 5, 6, 7]])
    return mesh_data, native_nodes, initial_tets, post_tets, cohesive


def test_x_plane_lineage_uses_parent_centroids_and_positive_area() -> None:
    mesh_data, native_nodes, initial_tets, post_tets, cohesive = _lineage_inputs()
    result = build_native_lineage(
        mesh_data,
        native_nodes,
        post_tets,
        cohesive,
        expected_area_m2=0.5,
        native_initial_bulk_tets=initial_tets,
        plane_axis=0,
        plane_coordinate_m=0.0,
    )

    facet = result["parent_facets"][0]
    assert result["positive_rest_area_sum_m2"] == pytest.approx(0.5)
    assert facet["positive_rest_area_m2"] == pytest.approx(0.5)
    assert facet["side_a_positive_rest_area_m2"] > 0.0
    assert facet["side_b_positive_rest_area_m2"] > 0.0
    assert set(facet["side_parent_by_geometric_side"]) == {"negative", "positive"}
    assert facet["side_parent_by_geometric_side"] == {"negative": 101, "positive": 202}
    assert facet["post_insertion_bulk_parent_rows"] == [0, 1]
    assert result["fracture_status_claimed"] is False


def test_initial_connectivity_must_match_source_element_and_local_order() -> None:
    mesh_data, native_nodes, initial_tets, post_tets, cohesive = _lineage_inputs()
    wrong_order = initial_tets.copy()
    wrong_order[[0, 1]] = wrong_order[[1, 0]]
    with pytest.raises(ValueError, match="initial getConnectivity"):
        build_native_lineage(
            mesh_data,
            native_nodes,
            post_tets,
            cohesive,
            expected_area_m2=0.5,
            native_initial_bulk_tets=wrong_order,
        )

    wrong_local_ordinal = initial_tets.copy()
    wrong_local_ordinal[0, [0, 1]] = wrong_local_ordinal[0, [1, 0]]
    with pytest.raises(ValueError, match="initial getConnectivity"):
        build_native_lineage(
            mesh_data,
            native_nodes,
            post_tets,
            cohesive,
            expected_area_m2=0.5,
            native_initial_bulk_tets=wrong_local_ordinal,
        )


def test_wrong_post_insertion_ancestor_id_fails_exact_lineage() -> None:
    mesh_data, native_nodes, initial_tets, post_tets, cohesive = _lineage_inputs()
    wrong_ancestor = post_tets.copy()
    wrong_ancestor[1, 3] = 3
    with pytest.raises(ValueError, match="ancestry|coordinate"):
        build_native_lineage(
            mesh_data,
            native_nodes,
            wrong_ancestor,
            cohesive,
            expected_area_m2=0.5,
            native_initial_bulk_tets=initial_tets,
        )


def test_two_post_insertion_parents_on_a_cohesive_side_fail() -> None:
    mesh_data, native_nodes, initial_tets, post_tets, cohesive = _lineage_inputs()
    uncut_neighbor = post_tets.copy()
    uncut_neighbor[1] = [4, 0, 2, 1]
    duplicate_side = np.asarray([[0, 1, 2, 0, 1, 2]])
    with pytest.raises(ValueError, match="2 post-insertion bulk parents"):
        build_native_lineage(
            mesh_data,
            native_nodes[:5],
            uncut_neighbor,
            duplicate_side,
            expected_area_m2=0.5,
            native_initial_bulk_tets=initial_tets,
        )


def test_missing_cohesive_face_fails_bijection() -> None:
    mesh_data, native_nodes, initial_tets, post_tets, _cohesive = _lineage_inputs()
    with pytest.raises(ValueError, match="bijection"):
        build_native_lineage(
            mesh_data,
            native_nodes,
            post_tets,
            np.empty((0, 6), dtype=np.int64),
            expected_area_m2=0.5,
            native_initial_bulk_tets=initial_tets,
        )


def test_input_face_requires_two_and_only_two_source_parents() -> None:
    mesh_data = _mesh_data()
    mesh_data["tetrahedra"].append({"element_id": 303, "nodes": [5, 0, 1, 2]})
    mesh_data["nodes"] = np.vstack([mesh_data["nodes"], [[-1.0, 0.6, 0.2]]])
    with pytest.raises(ValueError, match="nonmanifold source tetrahedron face"):
        input_interface_contract(mesh_data)


def test_gmsh_reader_keeps_exterior_physical_triangles(tmp_path) -> None:
    path = tmp_path / "source.msh"
    path.write_text(
        "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n"
        "$PhysicalNames\n2\n2 1 \"interface\"\n2 2 \"exterior\"\n$EndPhysicalNames\n"
        "$Nodes\n5\n"
        "1 0 0 0\n2 0 1 0\n3 0 0 1\n4 -1 0.2 0.2\n5 1 0.4 0.2\n"
        "$EndNodes\n"
        "$Elements\n3\n"
        "10 4 0 4 1 2 3\n"
        "11 4 0 5 1 3 2\n"
        "12 2 2 1 7 1 2 3\n"
        "$EndElements\n",
        encoding="utf-8",
    )
    mesh_data = read_gmsh22(path)
    assert len(mesh_data["triangles"]) == 1
    assert mesh_data["triangles"][0]["physical_name"] == "interface"
    assert mesh_data["triangles"][0]["nodes"] == [0, 1, 2]
    # Interface triangles are deliberately not put in the exterior subset.
    assert mesh_data["exterior_triangles"] == []

    exterior_path = tmp_path / "exterior.msh"
    exterior_path.write_text(
        path.read_text(encoding="utf-8").replace(
            "$Elements\n3\n", "$Elements\n4\n"
        ).replace(
            "$EndElements\n",
            "13 2 2 2 8 1 2 4\n$EndElements\n",
        ),
        encoding="utf-8",
    )
    exterior_mesh = read_gmsh22(exterior_path)
    assert [row["physical_name"] for row in exterior_mesh["exterior_triangles"]] == ["exterior"]
    assert exterior_mesh["exterior_triangles"][0]["nodes"] == [0, 1, 3]


def test_gmsh_reader_rejects_unsupported_surface_element_type(tmp_path) -> None:
    path = tmp_path / "quadratic.msh"
    path.write_text(
        "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n"
        "$Nodes\n6\n"
        "1 0 0 0\n2 1 0 0\n3 0 1 0\n4 0.5 0 0\n5 0.5 0.5 0\n6 0 0.5 0\n"
        "$EndNodes\n"
        "$Elements\n1\n20 9 0 1 2 3 4 5 6\n$EndElements\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="unsupported surface element type 9"):
        read_gmsh22(path)


def test_triangle_area_and_quadrature_signed_work() -> None:
    mesh_data = _mesh_data()
    assert triangle_area(mesh_data["nodes"], [0, 1, 2]) == pytest.approx(0.5)
    work = quadrature_signed_work(
        np.zeros((1, 3)),
        np.asarray([[2.0, 0.0, 0.0]]),
        np.asarray([[0.25, 0.0, 0.0]]),
        np.asarray([[4.0, 0.0, 0.0]]),
        np.asarray([0.5]),
    )
    assert work == pytest.approx([0.375])
