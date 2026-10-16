"""Small analytic checks for the independent native-vessel audit helpers."""
from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pytest


_AUDIT_DIR = Path(__file__).resolve().parents[1] / "experiments" / "cohesive_vessel_20261015"
sys.path.insert(0, str(_AUDIT_DIR))
import audit_native_vessel as audit  # noqa: E402


def test_qp_traction_work_integrates_vector_opening_with_area() -> None:
    opening = np.array([[[0.0, 0.0, 0.0]], [[1.0, 0.0, 0.0]], [[2.0, 0.0, 0.0]]])
    traction = np.array([[[2.0, 0.0, 0.0]], [[4.0, 0.0, 0.0]], [[2.0, 0.0, 0.0]]])

    cumulative = audit.qp_trapezoid_work(opening, traction, np.array([0.5]))

    # Each interval has mean traction 3 Pa, 1 m opening, and 0.5 m^2 area.
    np.testing.assert_allclose(cumulative, [0.0, 1.5, 3.0], rtol=0.0, atol=1e-14)


def test_t4_mass_and_affine_small_strain_energy_match_closed_form() -> None:
    rest = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        dtype=float,
    )
    strain = 0.01
    translation = np.array([0.2, -0.3, 0.4])
    displacement = strain * rest + translation

    mass, volume, energy = audit.tetra_lumped_mass_and_strain_energy(
        rest,
        np.array([[0, 1, 2, 3]]),
        displacement,
        rho_kg_m3=12.0,
        young_pa=120.0,
        poisson=0.25,
    )

    assert volume == pytest.approx(1.0 / 6.0)
    np.testing.assert_allclose(mass, np.full(4, 0.5), rtol=0.0, atol=1e-14)
    # epsilon = diag(0.01, 0.01, 0.01); lambda = mu = 48 Pa.
    assert energy == pytest.approx(0.006, rel=1e-13, abs=1e-15)


def test_t4_rejects_inverted_connectivity_instead_of_taking_absolute_volume() -> None:
    rest = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        dtype=float,
    )

    with pytest.raises(ValueError, match="positive orientation"):
        audit.tetra_lumped_mass_and_strain_energy(
            rest,
            np.array([[0, 2, 1, 3]]),
            np.zeros_like(rest),
            rho_kg_m3=1.0,
            young_pa=1.0,
            poisson=0.25,
        )


def test_external_work_uses_signed_patch_weights_and_pulse_profile() -> None:
    work, force = audit.reconstruct_external_work(
        patch_displacement_x_m=np.array([[0.0], [1.0], [2.0]]),
        unit_external_load_x=np.array([1.0]),
        times_s=np.array([0.0, 0.5, 1.0]),
        duration_s=1.0,
        peak_force_n=2.0,
    )

    np.testing.assert_allclose(force, [0.0, 2.0, 0.0], rtol=0.0, atol=1e-14)
    np.testing.assert_allclose(work, [0.0, 1.0, 2.0], rtol=0.0, atol=1e-14)


def test_central_difference_recurrence_checks_both_position_and_velocity() -> None:
    position, velocity = audit.central_difference_residuals(
        u_previous=np.array([[0.0]]),
        u_current=np.array([[1.5]]),
        v_previous=np.array([[2.0]]),
        v_current=np.array([[4.5]]),
        a_previous=np.array([[4.0]]),
        a_current=np.array([[6.0]]),
        dt_s=0.5,
    )

    np.testing.assert_allclose(position, 0.0, rtol=0.0, atol=1e-14)
    np.testing.assert_allclose(velocity, 0.0, rtol=0.0, atol=1e-14)


def test_contact_penalty_energy_and_closed_cycle_work_are_reversible() -> None:
    contact_opening = np.array([[[0.0, 0.0, 0.0]], [[0.1, 0.0, 0.0]], [[0.0, 0.0, 0.0]]])
    contact_traction = 2.0 * contact_opening
    area = np.array([2.0])

    ec_by_qp = audit.quadrature_contact_energy_by_point(contact_opening, contact_traction, area)
    work = audit.interface_work_decomposition(
        opening_m=np.zeros_like(contact_opening),
        traction_pa=np.zeros_like(contact_traction),
        contact_opening_m=contact_opening,
        contact_traction_pa=contact_traction,
        qp_area_m2=area,
    )

    np.testing.assert_allclose(ec_by_qp.sum(axis=1), [0.0, 0.02, 0.0], rtol=0.0, atol=1e-14)
    np.testing.assert_allclose(work[:, 2], [0.0, 0.02, 0.0], rtol=0.0, atol=1e-14)
    np.testing.assert_allclose(work[:, 0], work[:, 2], rtol=0.0, atol=1e-14)
    np.testing.assert_allclose(work[:, 1], 0.0, rtol=0.0, atol=1e-14)
    np.testing.assert_allclose(work[:, 3], 0.0, rtol=0.0, atol=1e-14)


def test_interface_work_reports_nonzero_cross_work_separately() -> None:
    opening = np.array([[[0.0, 0.0, 0.0]], [[0.1, 0.0, 0.0]]])
    traction = np.array([[[0.0, 0.0, 0.0]], [[0.0, 0.2, 0.0]]])
    contact_opening = np.array([[[0.0, 0.0, 0.0]], [[0.0, 0.1, 0.0]]])
    contact_traction = np.array([[[0.0, 0.0, 0.0]], [[0.3, 0.0, 0.0]]])

    work = audit.interface_work_decomposition(
        opening, traction, contact_opening, contact_traction, np.array([1.0])
    )

    assert work[-1, 1] == pytest.approx(0.0, abs=1e-14)
    assert work[-1, 2] == pytest.approx(0.0, abs=1e-14)
    assert work[-1, 3] == pytest.approx(0.025, rel=1e-13)
    assert work[-1, 0] == pytest.approx(work[-1, 3], rel=1e-13)


def test_qp_failure_requires_zero_total_traction_and_remains_zero() -> None:
    damage = np.array([[0.99], [1.0], [1.0], [1.0]])
    traction = np.array([[[2.0, 0.0, 0.0]], [[1.0, 0.0, 0.0]], [[0.0, 0.0, 0.0]], [[3.0, 0.0, 0.0]]])
    contact = np.zeros_like(traction)

    first, report = audit.qp_failure_steps(damage, traction, contact, 0.01, 0.1)

    np.testing.assert_array_equal(first, [2])
    assert report["all_qps_reached_full_damage_and_zero_total_traction"]
    assert not report["post_failure_total_traction_stays_zero"]
    assert report["post_failure_total_traction_max_norm_Pa"] == pytest.approx(3.0)
