import json
import math
import sys

import numpy as np
import pytest

from pgfrag.mesh_resolution_audit import audit_tetra_mesh


def regular_tetra(edge_m=1.0e-3):
    a = edge_m
    return np.asarray(
        [
            [0.0, 0.0, 0.0],
            [a, 0.0, 0.0],
            [a / 2.0, math.sqrt(3.0) * a / 2.0, 0.0],
            [a / 2.0, math.sqrt(3.0) * a / 6.0, math.sqrt(2.0 / 3.0) * a],
        ],
        dtype=np.float64,
    )


def explicit_policy(**overrides):
    return {
        "id": "unit-test-resolution-policy-v1",
        "basis": "test-only explicit dimensional comparison limits",
        **overrides,
    }


def test_regular_tetrahedron_reports_si_volume_edges_and_heights():
    report = audit_tetra_mesh(regular_tetra(), [[0, 1, 2, 3]])
    geometry = report["tetra_geometry"]
    expected_volume = math.sqrt(2.0) * (1.0e-3**3) / 12.0
    assert geometry["signed_volume_m3"]["median"] == pytest.approx(expected_volume)
    assert geometry["signed_orientation_counts"] == {"positive": 1, "negative": 0, "zero": 0, "unknown": 0}
    assert geometry["pooled_tetra_edge_length_m"]["median"] == pytest.approx(1.0e-3)
    assert geometry["minimum_tetra_height_m"]["median"] == pytest.approx(math.sqrt(2.0 / 3.0) * 1.0e-3)
    assert geometry["mesh_characteristic_h_m"]["value"] == pytest.approx(1.0e-3)
    assert geometry["all_signed_orientations_positive"] is True
    surface_volume = report["surface_vs_volume_discretization"]
    assert surface_volume["surface_tessellation_tolerance_status"] == "NOT_PROVIDED"
    assert surface_volume["surface_triangle_edge_length_status"] == "NOT_PROVIDED"
    assert surface_volume["volume_tetra_spacing_m"] == pytest.approx(1.0e-3)
    json.dumps(report, allow_nan=False)


def test_sliver_is_separately_reported_as_near_degenerate():
    vertices = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0e-15, 1.0e-15, 1.0e-15]]
    report = audit_tetra_mesh(vertices, [[0, 1, 2, 3]])
    cell = report["tetra_geometry"]["per_cell"][0]
    assert cell["orientation"] == "POSITIVE"
    assert cell["status"] == "NEAR_DEGENERATE"
    assert report["tetra_geometry"]["near_degenerate_count"] == 1


def test_negative_orientation_and_zero_volume_are_distinct():
    vertices = np.vstack([regular_tetra(), [[0.0, 0.0, 0.0]]])
    report = audit_tetra_mesh(vertices, [[0, 1, 2, 3], [0, 2, 1, 3], [0, 1, 2, 4]])
    geometry = report["tetra_geometry"]
    assert geometry["per_cell"][1]["status"] == "NEGATIVE_ORIENTATION"
    assert geometry["per_cell"][1]["orientation"] == "NEGATIVE"
    assert geometry["per_cell"][2]["status"] == "ZERO_VOLUME"
    assert geometry["per_cell"][2]["orientation"] == "ZERO"
    assert geometry["signed_orientation_counts"] == {"positive": 1, "negative": 1, "zero": 1, "unknown": 0}
    assert geometry["zero_volume_count"] == 1
    assert geometry["all_signed_orientations_positive"] is False


def test_caller_declared_nearwall_candidate_screen_and_normal_resolution():
    report = audit_tetra_mesh(
        regular_tetra(),
        [[0, 1, 2, 3]],
        wall_measurement={
            "scope": "global_interval",
            "thickness_interval_m": [2.0e-4, 3.0e-4],
            "definition": "source outer-to-inner wall interval",
            "source": "hand-constructed analytical fixture",
        },
        wall_normal_map={
            "normals_by_cell": [[1.0, 0.0, 0.0]],
            "definition": "local wall normal",
            "source": "hand-constructed analytical fixture",
        },
        resolution_policy=explicit_policy(
            max_h_to_wall=0.5,
            max_normal_projected_edge_to_wall=0.5,
        ),
    )
    wall = report["wall_resolution"]
    assert wall["h_to_wall"]["value"] == pytest.approx(5.0)
    assert wall["global_wall_status"] == "UNKNOWN"
    assert wall["actual_source_verification"]["verified"] is False
    assert wall["h_to_wall"]["policy_assessment"]["status"] == "UNKNOWN"
    assert wall["h_to_wall"]["caller_declared_candidate_policy_screen"]["status"] == "RESOLUTION_INSUFFICIENT"
    assert wall["local_normal_projection"]["projected_edge_m_quantiles"]["p95"] == pytest.approx(1.0e-3)
    assert wall["local_normal_projection"]["policy_assessment"]["status"] == "RESOLUTION_INSUFFICIENT"
    assert report["directional_diagnostics"]["single_mesh_crack_direction_bias_status"] == "UNKNOWN"


def test_tiny_regular_tetrahedron_keeps_scale_and_marks_volume_underflow_unknown():
    report = audit_tetra_mesh(regular_tetra(edge_m=1.0e-120), [[0, 1, 2, 3]])
    geometry = report["tetra_geometry"]
    cell = geometry["per_cell"][0]
    assert cell["orientation"] == "POSITIVE"
    assert cell["status"] == "UNKNOWN_NUMERIC_RANGE/UNDERFLOW"
    assert cell["signed_volume_m3"] is None
    assert cell["volume_computation_status"] == "UNDERFLOW"
    assert cell["maximum_edge_scale_m"] == pytest.approx(1.0e-120)
    assert cell["normalized_signed_volume_over_h_cubed"] > 0.0
    assert cell["signed_volume_log10_abs_m3"] == pytest.approx(math.log10(math.sqrt(2.0) / 12.0) - 360.0)
    issue = geometry["unknown_numeric_range_issues"][0]
    assert issue["cell_index"] == 0
    assert issue["status"] == "UNDERFLOW"
    assert issue["maximum_edge_scale_m"] == pytest.approx(1.0e-120)
    assert issue["normalized_signed_volume_over_h_cubed"] == pytest.approx(math.sqrt(2.0) / 12.0)
    assert issue["signed_volume_log10_abs_m3"] == pytest.approx(math.log10(math.sqrt(2.0) / 12.0) - 360.0)
    assert issue["signed_volume_m3"] is None
    assert geometry["all_nonzero_signed_volumes_representable"] is False
    json.dumps(report, allow_nan=False)


def test_extreme_delta_to_h_and_h_to_delta_range_failures_are_reported():
    report = audit_tetra_mesh(
        regular_tetra(edge_m=1.0e-100),
        [[0, 1, 2, 3]],
        pd_horizon_m=1.0e308,
    )
    horizon = report["pd_horizon"]
    assert horizon["h_over_delta"] is None
    assert horizon["h_over_delta_numeric_range_failure_count"] == 1
    assert horizon["h_over_delta_numeric_range_failure_cell_indices"] == [0]
    assert horizon["delta_over_h"] is None
    assert horizon["delta_over_h_status"] == "UNKNOWN_NUMERIC_RANGE"
    assert horizon["delta_over_h_numeric_range_failure_count"] == 1
    assert horizon["delta_over_h_numeric_range_failure_cell_indices"] == [0]
    json.dumps(report, allow_nan=False)


def test_overflow_has_one_terminal_status_count_per_cell():
    edge_m = 1.0e104
    expected_log_volume = math.log(math.sqrt(2.0) / 12.0) + 3.0 * math.log(edge_m)
    assert expected_log_volume > math.log(sys.float_info.max)
    report = audit_tetra_mesh(regular_tetra(edge_m=edge_m), [[0, 1, 2, 3]])
    geometry = report["tetra_geometry"]
    assert geometry["per_cell"][0]["status"] == "UNKNOWN_NUMERIC_RANGE/OVERFLOW"
    assert geometry["cell_status_counts"]["UNKNOWN_NUMERIC_RANGE/OVERFLOW"] == 1
    assert sum(geometry["cell_status_counts"].values()) == 1
    assert geometry["unknown_numeric_range_count"] == 1
    assert geometry["unknown_numeric_range_cell_indices"] == [0]


def test_nonfinite_scale_ratios_are_unknown_and_keep_cell_indices():
    a = 1.0e308
    vertices = [
        [0.0, 0.0, 0.0],
        [a, 0.0, 0.0],
        [a / 2.0, math.sqrt(3.0) * a / 2.0, 0.0],
        [a / 2.0, math.sqrt(3.0) * a / 6.0, math.sqrt(2.0 / 3.0) * a],
    ]
    report = audit_tetra_mesh(
        vertices,
        [[0, 1, 2, 3]],
        wall_measurement={
            "scope": "global_interval",
            "thickness_interval_m": [1.0e-308, 2.0e-308],
            "definition": "declared extreme interval for numeric-range test",
            "source": "synthetic test input",
        },
        process_zone={"length_m": 1.0e-308, "definition": "test length", "source": "synthetic test input"},
        pd_horizon_m=1.0e-308,
        resolution_policy={"id": "numeric-range-v1", "basis": "test", "max_h_to_wall": 0.5, "max_h_to_process_zone": 0.5, "max_h_to_pd_horizon": 0.5},
    )
    wall = report["wall_resolution"]["h_to_wall"]
    process = report["process_zone_resolution"]["h_to_process_zone"]
    horizon = report["pd_horizon"]
    assert wall["value"] is None
    assert wall["value_status"] == "UNKNOWN_NUMERIC_RANGE"
    assert wall["policy_assessment"]["status"] == "UNKNOWN"
    assert wall["numeric_range_failure_count"] == 1
    assert wall["numeric_range_failure_cell_indices"] == [0]
    assert process["value"] is None
    assert process["policy_assessment"]["status"] == "UNKNOWN"
    assert process["numeric_range_failure_cell_indices"] == [0]
    assert horizon["h_over_delta"] is None
    assert horizon["h_over_delta_policy_assessment"]["status"] == "UNKNOWN"
    assert horizon["h_over_delta_numeric_range_failure_cell_indices"] == [0]
    assert horizon["delta_over_h"] is None
    assert horizon["delta_over_h_status"] == "UNKNOWN_NUMERIC_RANGE"
    json.dumps(report, allow_nan=False)


def test_absent_normal_mapping_stays_unknown_without_edge_fallback():
    report = audit_tetra_mesh(
        regular_tetra(),
        [[0, 1, 2, 3]],
        wall_measurement={
            "scope": "global_interval",
            "thickness_interval_m": [2.0e-4, 3.0e-4],
            "definition": "source outer-to-inner wall interval",
            "source": "fixture",
        },
        resolution_policy=explicit_policy(max_h_to_wall=0.5, max_normal_projected_edge_to_wall=0.5),
    )
    normal = report["wall_resolution"]["local_normal_projection"]
    assert report["wall_resolution"]["h_to_wall"]["policy_assessment"]["status"] == "UNKNOWN"
    assert report["wall_resolution"]["h_to_wall"]["caller_declared_candidate_policy_screen"]["status"] == "RESOLUTION_INSUFFICIENT"
    assert normal["evidence_status"] == "UNKNOWN"
    assert normal["projected_edge_m_quantiles"]["count"] == 0
    assert normal["policy_assessment"]["status"] == "UNKNOWN"
    assert normal["policy_assessment"]["ratio"] is None

    no_wall_evidence = audit_tetra_mesh(
        regular_tetra(),
        [[0, 1, 2, 3]],
        resolution_policy=explicit_policy(max_h_to_wall=0.5, max_normal_projected_edge_to_wall=0.5),
    )
    assert no_wall_evidence["wall_resolution"]["global_wall_status"] == "UNKNOWN"
    assert no_wall_evidence["wall_resolution"]["h_to_wall"]["policy_assessment"]["status"] == "UNKNOWN"
    assert no_wall_evidence["wall_resolution"]["local_normal_projection"]["policy_assessment"]["status"] == "UNKNOWN"


def test_missing_process_zone_stays_unknown_and_does_not_reuse_pd_horizon():
    report = audit_tetra_mesh(
        regular_tetra(),
        [[0, 1, 2, 3]],
        pd_horizon_m=4.0e-4,
        resolution_policy=explicit_policy(max_h_to_process_zone=0.2, max_h_to_pd_horizon=0.5),
    )
    assert report["process_zone_resolution"]["evidence_status"] == "UNKNOWN"
    assert report["process_zone_resolution"]["h_to_process_zone"]["value"] is None
    assert report["process_zone_resolution"]["h_to_process_zone"]["policy_assessment"]["status"] == "UNKNOWN"
    assert report["pd_horizon"]["delta_over_h"] == pytest.approx(0.4)
    assert report["pd_horizon"]["h_over_delta"] == pytest.approx(2.5)
    assert report["pd_horizon"]["h_over_delta_policy_assessment"]["status"] == "RESOLUTION_INSUFFICIENT"


def test_invalid_indices_and_nonfinite_coordinates_are_separate_and_json_safe():
    vertices = [
        [0.0, 0.0, 0.0],
        [1.0e-3, 0.0, 0.0],
        [5.0e-4, 8.660254037844386e-4, 0.0],
        [5.0e-4, 2.886751345948129e-4, 8.164965809277261e-4],
        [float("nan"), 0.0, 0.0],
    ]
    report = audit_tetra_mesh(vertices, [[0, 1, 2, 3], [0, 1, 2, 99], [0, 1, 2, 4]])
    integrity = report["input_integrity"]
    assert integrity["invalid_index_cell_count"] == 1
    assert integrity["nonfinite_vertex_count"] == 1
    assert integrity["nonfinite_coordinate_cell_count"] == 1
    assert integrity["has_invalid_indices"] is True
    assert integrity["has_nonfinite_vertices"] is True
    statuses = [cell["status"] for cell in report["tetra_geometry"]["per_cell"]]
    assert statuses == ["OK", "INVALID_INDEX", "NONFINITE_COORDINATE"]
    json.dumps(report, allow_nan=False)
