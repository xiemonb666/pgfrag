import numpy as np
import pytest

from pgfrag.cohesive_stability import (
    _tetrahedron_stiffness,
    mass_normalized_stiffness_bound,
    tri3_cohesive_stiffness_matrix,
)


def test_skew_tetra_affine_energy_matches_isotropic_continuum():
    # A skew Jacobian catches the row/column shape-gradient distinction.
    x = np.array([[0., 0., 0.], [.002, 0., 0.], [.0004, .001, 0.], [.0001, .0003, .0015]])
    young, nu = 70e9, .22
    stiffness, volume = _tetrahedron_stiffness(x, young, nu)
    strain = np.array([[1e-4, 3e-5, 0.], [3e-5, 2e-4, 0.], [0., 0., -1e-4]])
    displacement = (x @ strain.T).ravel()
    mu = young / (2 * (1 + nu))
    lam = young * nu / ((1 + nu) * (1 - 2 * nu))
    expected = volume * (.5 * lam * np.trace(strain)**2 + mu * np.sum(strain**2))
    assert .5 * displacement @ stiffness @ displacement == pytest.approx(expected, rel=1e-12)
    rigid = np.tile([.01, -.02, .03], 4)
    assert np.linalg.norm(stiffness @ rigid) < 1e-7


def test_uniform_interface_opening_energy_is_single_area_energy():
    side = np.array([[0., 0., 0.], [.002, 0., 0.], [0., .002, 0.]])
    stiffness, area, normal = tri3_cohesive_stiffness_matrix(np.stack([side, side]), 1e15, 3e15)
    delta = 1e-8
    displacement = np.concatenate([np.zeros(9), np.tile(delta * normal, 3)])
    assert .5 * displacement @ stiffness @ displacement == pytest.approx(.5 * 1e15 * area * delta**2)
    assert area == pytest.approx(2e-6)


def _tiny_mesh():
    side = np.array([[0., 0., 0.], [.002, 0., 0.], [0., .002, 0.]])
    nodes = np.concatenate([side, [[0., 0., -.002]], side, [[0., 0., .002]]])
    tets = np.array([[0, 1, 2, 3], [4, 5, 6, 7]], dtype=np.int32)
    pairs = np.array([[[0, 1, 2], [4, 5, 6]]], dtype=np.int32)
    return nodes, tets, pairs


def test_mass_bound_covers_dense_eigenvalue_and_density_scaling():
    nodes, tets, pairs = _tiny_mesh()
    args = (nodes, tets, pairs, 2500., 70e9, .22, 1e15, 9e15, np.zeros((8, 3), bool))
    result = mass_normalized_stiffness_bound(*args, validate_dense_eigenvalue=True)
    assert result.lambda_upper_bound_per_s2 >= result.dense_max_eigenvalue_per_s2
    twice = mass_normalized_stiffness_bound(nodes, tets, pairs, 5000., *args[4:])
    assert twice.lambda_upper_bound_per_s2 == pytest.approx(result.lambda_upper_bound_per_s2 / 2)
    assert twice.critical_time_step_s == pytest.approx(result.critical_time_step_s * np.sqrt(2))


def test_native_mass_and_duplicate_face_nodes_fail_closed():
    nodes, tets, pairs = _tiny_mesh()
    args = (nodes, tets, pairs, 2500., 70e9, .22, 1e15, 9e15, np.zeros((8, 3), bool))
    with pytest.raises(ValueError):
        mass_normalized_stiffness_bound(*args, native_point_mass_kg=np.ones(8))
    pairs[0, 0, 1] = pairs[0, 0, 0]
    with pytest.raises(ValueError, match="distinct"):
        mass_normalized_stiffness_bound(*args)
