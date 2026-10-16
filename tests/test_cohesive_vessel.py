from __future__ import annotations

import numpy as np
import pytest
from types import SimpleNamespace

from pgfrag.cohesive_vessel import (
    contact_energy_check,
    deformation_metrics,
    fields,
    interface_work_increment,
    load_factor,
    make_balanced_patch_load,
    positive_centroid_weights,
    validate_contact_penalty,
)


def test_positive_centroid_weights_project_irregular_area_weights_to_target() -> None:
    points = np.asarray([
        [-1.2, -0.8],
        [-0.9, 1.4],
        [0.4, -1.1],
        [1.8, 0.2],
        [0.3, 2.0],
    ])
    base_area_weights = np.asarray([0.3, 1.7, 4.2, 0.8, 2.6])
    target = np.asarray([0.1, 0.25])

    weights = positive_centroid_weights(points, base_area_weights, target)

    assert np.all(weights > 0.0)
    assert weights.sum() == pytest.approx(1.0, abs=1e-13)
    assert weights @ points == pytest.approx(target, abs=1e-11)


def test_positive_centroid_weights_reject_infeasible_target() -> None:
    points = np.asarray([[-1.0, -1.0], [-1.0, 1.0], [1.0, -1.0], [1.0, 1.0]])
    with pytest.raises(ValueError):
        positive_centroid_weights(points, np.ones(4), np.asarray([2.0, 0.0]))


def _two_patch_source_and_native_nodes() -> tuple[dict, np.ndarray, np.ndarray, np.ndarray]:
    # Two square patches face outward in +/-x. Their common yz centroid is the
    # origin, so unit opposite x forces have zero net force and reference torque.
    yz = np.asarray([[-1.0, -1.0], [-1.0, 1.0], [1.0, -1.0], [1.0, 1.0]])
    native_nodes = np.vstack([
        np.column_stack((np.full(4, -1.0), yz)),
        np.column_stack((np.full(4, 1.0), yz)),
    ])
    triangles = [
        {"physical_name": "exterior", "nodes": [0, 1, 2]},
        {"physical_name": "exterior", "nodes": [1, 2, 3]},
        {"physical_name": "exterior", "nodes": [4, 5, 6]},
        {"physical_name": "exterior", "nodes": [5, 6, 7]},
    ]
    source = {"nodes": native_nodes.copy(), "triangles": triangles}
    ancestry = np.arange(len(native_nodes), dtype=np.int64)
    mass = np.ones(len(native_nodes), dtype=float)
    return source, native_nodes, ancestry, mass


def test_balanced_patch_load_has_unit_opposed_forces_and_common_reference_line() -> None:
    source, native_nodes, ancestry, mass = _two_patch_source_and_native_nodes()

    load, report = make_balanced_patch_load(source, native_nodes, ancestry, mass, threshold_m=0.5)

    assert load.shape == native_nodes.shape
    assert np.all(load[:, 1:] == 0.0)
    assert load[:, 0].min() < 0.0 < load[:, 0].max()
    assert load[:, 0].sum() == pytest.approx(0.0, abs=1e-13)
    assert np.cross(native_nodes, load).sum(axis=0) == pytest.approx(np.zeros(3), abs=1e-12)
    assert report["reference_unit_force_sum"] == pytest.approx([0.0, 0.0, 0.0], abs=1e-13)
    assert report["reference_unit_torque_sum_m"] == pytest.approx([0.0, 0.0, 0.0], abs=1e-12)
    assert report["target_reference_com_m"] == pytest.approx([0.0, 0.0, 0.0])
    assert len(report["patches"]) == 2
    assert {patch["sign"] for patch in report["patches"]} == {-1, 1}
    for patch in report["patches"]:
        assert patch["weights"] == pytest.approx([abs(load[node, 0]) for node in patch["native_node_ids"]])
        assert sum(patch["weights"]) == pytest.approx(1.0, abs=1e-13)
        assert patch["reference_force_line_yz_m"] == pytest.approx([0.0, 0.0], abs=1e-11)


def test_sin_squared_load_pulse_integrates_to_peak_times_duration_over_two() -> None:
    peak, duration = 37.0, 0.42
    time = np.linspace(0.0, duration, 100_001)
    values = np.asarray([load_factor(t, duration, peak) for t in time])
    integral = float(np.sum(0.5 * (values[:-1] + values[1:]) * np.diff(time)))

    assert load_factor(0.0, duration, peak) == 0.0
    assert load_factor(duration, duration, peak) == 0.0
    assert load_factor(-0.1, duration, peak) == 0.0
    assert load_factor(duration + 0.1, duration, peak) == 0.0
    assert integral == pytest.approx(peak * duration / 2.0, rel=1e-9)


def test_affine_skew_tet_metrics_match_strain_stretch_and_polar_rotation() -> None:
    # The reference tetrahedron is skew, so this exercises the local edge-matrix
    # inverse instead of an axis-aligned special case.
    rest = np.asarray([
        [0.2, -0.3, 0.5],
        [1.3, -0.1, 0.4],
        [0.45, 0.85, 0.7],
        [0.1, 0.15, 1.55],
    ])
    angle = 0.17
    rotation = np.asarray([
        [np.cos(angle), -np.sin(angle), 0.0],
        [np.sin(angle), np.cos(angle), 0.0],
        [0.0, 0.0, 1.0],
    ])
    stretch = np.asarray([
        [1.012, 0.004, -0.002],
        [0.004, 0.991, 0.003],
        [-0.002, 0.003, 1.006],
    ])
    deformation = rotation @ stretch
    translation = np.asarray([3.2, -1.4, 0.7])
    current = rest @ deformation.T + translation
    tets = np.asarray([[0, 1, 2, 3]], dtype=np.int64)

    metrics = deformation_metrics(rest, current, tets)

    strain = 0.5 * (deformation + deformation.T) - np.eye(3)
    expected_principal_strain = np.linalg.eigvalsh(strain)
    singular_values = np.linalg.svd(deformation, compute_uv=False)
    assert metrics["max_abs_principal_small_strain"] == pytest.approx(
        np.max(np.abs(expected_principal_strain)), rel=1e-12
    )
    assert metrics["min_det_F"] == pytest.approx(np.linalg.det(deformation), rel=1e-12)
    assert metrics["max_polar_rotation_rad"] == pytest.approx(angle, abs=1e-12)
    assert metrics["min_principal_stretch"] == pytest.approx(singular_values.min(), rel=1e-12)
    assert metrics["max_principal_stretch"] == pytest.approx(singular_values.max(), rel=1e-12)

    shifted = deformation_metrics(rest + translation, current + translation, tets)
    assert shifted == pytest.approx(metrics, rel=1e-12, abs=1e-12)


def test_fields_reads_native_contact_opening_as_three_qp_vectors() -> None:
    values = {
        "opening": np.zeros((3, 3)),
        "tractions": np.ones((3, 3)),
        "damage": np.zeros((3, 1)),
        "delta max": np.zeros((3, 1)),
        "contact_tractions": -np.ones((3, 3)),
        "contact_opening": 2.0 * np.ones((3, 3)),
    }

    class FakeModel:
        def __init__(self):
            self.requested_fields = []

        def flattenInternalReal(self, name, element_kind, ghost_kind):
            self.requested_fields.append((name, element_kind, ghost_kind))
            return lambda element_type, ghost_flag: values[name]

    model = FakeModel()
    aka = SimpleNamespace(
        _ek_cohesive="cohesive_kind",
        _not_ghost="not_ghost",
        _cohesive_3d_6="cohesive_tri3",
    )
    result = fields(model, aka, facets=1)

    assert ("contact_opening", "cohesive_kind", "not_ghost") in model.requested_fields
    assert result["contact_opening"].shape == (3, 3)
    assert np.array_equal(result["contact_opening"], values["contact_opening"])


def test_normal_compression_contact_energy_matches_quadrature_formula() -> None:
    # Three equal quadrature shares on a 0.06 m² face, each with a 1e-6 m
    # normal compression opening and matching 1e9 Pa contact traction.
    contact_opening = np.tile([0.0, 0.0, -1e-6], (3, 1))
    contact_traction = np.tile([0.0, 0.0, -1e9], (3, 1))
    qp_area = np.full(3, 0.02)

    check = contact_energy_check(contact_opening, contact_traction, qp_area, native_energy_J=29.0)

    assert check == pytest.approx([30.0, 29.0, 1.0], abs=1e-12)


def test_full_interface_work_records_cohesive_contact_and_cross_terms() -> None:
    zero = np.zeros((1, 3))
    cohesive_opening_1 = np.asarray([[0.3, 0.0, 0.0]])
    cohesive_traction_0 = np.asarray([[1.0, 0.0, 0.0]])
    cohesive_traction_1 = np.asarray([[3.0, 0.0, 0.0]])
    contact_opening_1 = np.asarray([[0.2, 0.0, 0.0]])
    contact_traction_1 = np.asarray([[4.0, 0.0, 0.0]])

    increments = interface_work_increment(
        zero,
        cohesive_traction_0,
        zero,
        zero,
        cohesive_opening_1,
        cohesive_traction_1,
        contact_opening_1,
        contact_traction_1,
        np.asarray([2.0]),
    )

    assert increments == pytest.approx([4.0, 1.2, 0.8, 2.0], abs=1e-13)
    assert increments[0] == pytest.approx(increments[1] + increments[2] + increments[3])


def test_fixed_normal_contact_penalty_must_fit_kn_envelope() -> None:
    assert validate_contact_penalty(0.0, 1e15) == 0.0
    assert validate_contact_penalty(1e15, 1e15) == 1e15
    for invalid in (-1.0, 1.01e15, np.inf, np.nan):
        with pytest.raises(ValueError, match="0 <= penalty <= Kn"):
            validate_contact_penalty(invalid, 1e15)
