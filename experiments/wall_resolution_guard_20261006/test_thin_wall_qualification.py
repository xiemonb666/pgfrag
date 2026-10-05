from __future__ import annotations

import sys
from pathlib import Path

import pytest


HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from thin_wall_qualification import (  # noqa: E402
    QualificationInputError,
    _length_m,
    check_segments,
    historical_h8_assessment,
    make_analytic_hollow_cylinder_fixture,
    qualify_wall_input,
)


def test_mm_converts_to_m_and_fine_fixture_proves_four_real_layers() -> None:
    assert _length_m({"value": 2.0, "unit": "mm"}, "test") == pytest.approx(0.002)
    report = qualify_wall_input(make_analytic_hollow_cylinder_fixture(0.5))
    assert report["resolution"]["local_status"] == "LOCAL_RESOLUTION_PRECONDITION_MET"
    assert report["resolution"]["min_observed_actual_thickness_intervals"] == 4
    assert report["resolution"]["min_observed_distinct_centroid_layers"] == 4
    assert report["resolution"]["wall_samples"][0]["wall_thickness_m"] == pytest.approx(0.002)
    assert report["resolution"]["physical_pass"] is False
    assert report["local_scope"]["global_pass"] is False
    assert report["resolution"]["actual_h_over_t_distribution"]["min"] == pytest.approx(0.25)
    assert report["resolution"]["actual_h_over_t_distribution"]["max"] == pytest.approx(0.25)


def test_unsupported_or_dimensionally_bad_units_fail_closed() -> None:
    payload = make_analytic_hollow_cylinder_fixture(0.5)
    payload["cells"][0]["max_edge"]["unit"] = "inch"
    with pytest.raises(QualificationInputError, match="unsupported"):
        qualify_wall_input(payload)

    payload = make_analytic_hollow_cylinder_fixture(0.5)
    payload["units"]["volume"] = "mm"
    with pytest.raises(QualificationInputError, match="unsupported declared volume unit"):
        qualify_wall_input(payload)


def test_positive_h_cubed_pointweights_do_not_pass_geometry_weight_reconciliation() -> None:
    payload = make_analytic_hollow_cylinder_fixture(0.5)
    for row in payload["centroid_pointweights"]:
        row["point_weight"] = {"value": 0.5**3, "unit": "mm^3"}
        row["mass"] = {
            "value": payload["density_kg_m3"] * 0.5**3 * 1e-9,
            "unit": "kg",
        }

    report = qualify_wall_input(payload)
    assert report["volume_mass_consistency"]["pointweight_volume_delta_m3"] != pytest.approx(0.0, abs=1e-16)
    assert report["volume_mass_consistency"]["per_cell_pointweight_mismatches"]
    assert report["volume_mass_consistency"]["volume_status"] == "FAIL_VOLUME_OR_WEIGHT_RECONCILIATION"
    assert report["resolution"]["local_status"] == "BLOCKED_INCONSISTENT_VOLUME_OR_MASS"


def test_coarse_h_equals_wall_is_underresolved_from_actual_one_interval() -> None:
    report = qualify_wall_input(make_analytic_hollow_cylinder_fixture(2.0))
    assert report["resolution"]["local_status"] == "UNDERRESOLVED"
    assert report["resolution"]["min_observed_actual_thickness_intervals"] == 1
    assert report["resolution"]["min_observed_distinct_centroid_layers"] == 1
    assert report["resolution"]["actual_h_over_t_distribution"]["max"] == pytest.approx(1.0)
    old = historical_h8_assessment()
    assert old["status"] == "UNDERRESOLVED"
    assert old["h_over_t"] == pytest.approx(4.0)


def test_missing_coverage_is_unknown_and_never_global_pass() -> None:
    payload = make_analytic_hollow_cylinder_fixture(0.5, coverage_present=False)
    report = qualify_wall_input(payload)
    assert report["resolution"]["local_status"] == "QUALIFICATION_UNKNOWN_MISSING_OR_INCOMPLETE_LOCAL_COVERAGE"
    assert report["local_scope"]["global_status"] == "UNKNOWN_NO_FULL_DOMAIN_COVERAGE"
    assert report["local_scope"]["global_pass"] is False
    assert report["resolution"]["physical_status"] == "NOT_QUALIFIED"
    assert report["resolution"]["physical_pass"] is False

    # An input claim alone cannot certify the whole body; this evaluator only
    # proves the supplied local patch and sample design.
    payload = make_analytic_hollow_cylinder_fixture(0.5)
    payload["coverage"]["full_domain_coverage"] = True
    claimed = qualify_wall_input(payload)
    assert claimed["local_scope"]["full_domain_coverage_declared_by_input"] is True
    assert claimed["local_scope"]["full_domain_coverage_verified"] is False
    assert claimed["local_scope"]["global_status"] == "UNKNOWN_NO_FULL_DOMAIN_COVERAGE"
    assert claimed["local_scope"]["global_pass"] is False


def test_exact_quadratic_segment_check_catches_material_endpoints_through_cavity() -> None:
    payload = make_analytic_hollow_cylinder_fixture(0.5)
    segments = [
        {
            "segment_id": "long-through-hole",
            "source_ref": "test:two-material-centroids-connected-across-hollow",
            "start": {"value": [-9.0, 0.0, 5.0], "unit": "mm"},
            "end": {"value": [9.0, 0.0, 5.0], "unit": "mm"},
        }
    ]
    report = check_segments(segments, payload["geometry"])
    assert report["segments"][0]["endpoint_classification"] == "BOTH_IN_MATERIAL"
    assert report["status"] == "FAIL_CROSSINGS_DETECTED"
    assert report["segments"][0]["violating_intervals"][0]["reason"] == "CROSSES_HOLLOW_CAVITY"


def test_empty_segment_adapter_input_is_unknown_not_an_empty_pass() -> None:
    payload = make_analytic_hollow_cylinder_fixture(0.5)
    report = check_segments([], payload["geometry"])
    assert report["segment_count"] == 0
    assert report["status"] == "UNKNOWN_NO_SEGMENTS"
    payload["segments"] = []
    qualification = qualify_wall_input(payload)
    assert qualification["resolution"]["local_status"] == "QUALIFICATION_UNKNOWN_NO_SEGMENT_EVIDENCE"


def test_declared_but_incomplete_wall_sample_design_is_unknown() -> None:
    payload = make_analytic_hollow_cylinder_fixture(0.5, coverage_complete=False)
    payload["wall_samples"] = payload["wall_samples"][:-1]
    report = qualify_wall_input(payload)
    assert report["resolution"]["local_status"] == "QUALIFICATION_UNKNOWN_MISSING_OR_INCOMPLETE_LOCAL_COVERAGE"
    assert report["resolution"]["min_observed_actual_thickness_intervals"] is None
    assert report["local_scope"]["global_pass"] is False


def test_duplicate_rays_cannot_masquerade_as_complete_local_coverage() -> None:
    payload = make_analytic_hollow_cylinder_fixture(0.5)
    repeated = payload["wall_samples"][0]
    for index, sample in enumerate(payload["wall_samples"]):
        sample["normal"] = list(repeated["normal"])
        sample["inner_point"] = dict(repeated["inner_point"])
        sample["outer_point"] = dict(repeated["outer_point"])
        sample["cell_intersections"] = [dict(row) for row in repeated["cell_intersections"]]
        sample["sample_id"] = f"duplicate-ray-{index}"
    report = qualify_wall_input(payload)
    assert report["local_scope"]["wall_sample_design_complete"] is False
    assert report["resolution"]["local_status"] == "QUALIFICATION_UNKNOWN_MISSING_OR_INCOMPLETE_LOCAL_COVERAGE"


def test_fine_fixture_has_local_volume_and_mass_consistency_from_exact_cell_geometry() -> None:
    report = qualify_wall_input(make_analytic_hollow_cylinder_fixture(0.5))
    consistency = report["volume_mass_consistency"]
    assert consistency["volume_status"] == "PASS_LOCAL_GEOMETRY_POINTWEIGHT_RECONCILIATION"
    assert consistency["mass_status"] == "PASS_DENSITY_TIMES_GEOMETRIC_VOLUME"
    assert consistency["analytic_patch_reference_volume_m3"] == pytest.approx(
        consistency["sum_exact_cell_geometry_volumes_m3"], rel=1e-10
    )
    assert consistency["rho_times_reference_volume_kg"] == pytest.approx(
        consistency["point_mass_sum_kg"], rel=1e-10
    )

