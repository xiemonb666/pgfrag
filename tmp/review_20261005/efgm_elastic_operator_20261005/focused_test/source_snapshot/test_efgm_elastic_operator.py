"""Source-only acceptance tests for the bounded elastic Galerkin operator.

These tests are prepared for root review and have intentionally not been run.
The fixture uses an independent global affine-complete shape provider, not
the assembly's B/D implementation or the not-yet-integrated MLS adapter.
The provider spans only four global polynomial modes; therefore the fixture's
K may contain additional null modes. The rigid residual check verifies that
the six physical rigid modes are present, not that they are the only null modes
or that a local MLS/EFGM support is stable.
"""

from __future__ import annotations

import importlib.util
from itertools import product
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import numpy as np


MODULE_PATH = Path(__file__).with_name("efgm_elastic_operator.py")
SPEC = importlib.util.spec_from_file_location("pgfrag_efgm_elastic_operator", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
elastic = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = elastic
SPEC.loader.exec_module(elastic)


def bounded_fixture() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return 27 sites and 64 volume points on a side-2 cube."""
    axis = np.array([-1.0, 0.0, 1.0], dtype=np.float64)
    xx, yy, zz = np.meshgrid(axis, axis, axis, indexing="ij")
    nodes = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))

    gauss = (-1.0 / np.sqrt(3.0), 1.0 / np.sqrt(3.0))
    cell_centers = (-0.5, 0.5)
    quadrature = []
    for center in product(cell_centers, repeat=3):
        for local in product(gauss, repeat=3):
            quadrature.append(np.asarray(center) + 0.5 * np.asarray(local))
    qpoints = np.asarray(quadrature, dtype=np.float64)
    volume_weights = np.full(qpoints.shape[0], 0.125, dtype=np.float64)
    return nodes, qpoints, volume_weights


class GlobalAffineShapeProvider:
    """Independent dense shape rows that exactly reproduce affine fields."""

    def __init__(self, coefficient_positions: np.ndarray) -> None:
        self.node_count = coefficient_positions.shape[0]
        basis = np.column_stack(
            (np.ones(self.node_count, dtype=np.float64), coefficient_positions)
        )
        self.right_inverse = np.linalg.solve(basis.T @ basis, basis.T)
        self.support_ids = np.arange(self.node_count, dtype=np.int64)
        self.gradient_rows = self.right_inverse[1:, :].T.copy()

    def evaluate(self, material_position: np.ndarray):
        query_basis = np.concatenate(([1.0], material_position))
        weights = query_basis @ self.right_inverse
        return elastic.ShapeSample(
            weights=np.asarray(weights, dtype=np.float64),
            gradients=self.gradient_rows.copy(),
            support_ids=self.support_ids.copy(),
            status=elastic.ShapeStatus.OK,
        )


class PoorSupportProvider:
    def __init__(self, status: elastic.ShapeStatus, diagnostic: str | None = None) -> None:
        self.status = status
        self.diagnostic = diagnostic
        self.calls = 0

    def evaluate(self, material_position: np.ndarray):
        self.calls += 1
        return elastic.ShapeSample(None, None, None, self.status, self.diagnostic)


class MalformedShapeProvider:
    def __init__(self, coefficient_positions: np.ndarray, *, bad_ids: bool = False) -> None:
        self.node_count = coefficient_positions.shape[0]
        self.bad_ids = bad_ids

    def evaluate(self, material_position: np.ndarray):
        weights = np.full(self.node_count, 1.0 / self.node_count)
        gradient_columns = 3 if self.bad_ids else 2
        gradients = np.zeros((self.node_count, gradient_columns), dtype=np.float64)
        ids = np.arange(self.node_count, dtype=np.int64)
        if self.bad_ids:
            ids[-1] = 0
        return elastic.ShapeSample(weights, gradients, ids, elastic.ShapeStatus.OK)


class HugeGradientProvider:
    def __init__(self, node_count: int, gradient_value: float) -> None:
        self.node_count = node_count
        self.gradient_value = gradient_value
        self.support_ids = np.arange(node_count, dtype=np.int64)

    def evaluate(self, material_position: np.ndarray):
        gradients = np.zeros((self.node_count, 3), dtype=np.float64)
        gradients[0, 0] = self.gradient_value
        return elastic.ShapeSample(
            weights=np.zeros(self.node_count, dtype=np.float64),
            gradients=gradients,
            support_ids=self.support_ids.copy(),
            status=elastic.ShapeStatus.OK,
        )


class ExplodingProvider:
    def evaluate(self, material_position: np.ndarray):
        raise AssertionError("provider must not run after an input budget rejection")


class ElasticOperatorTests(unittest.TestCase):
    def setUp(self) -> None:
        self.nodes, self.qpoints, self.volumes = bounded_fixture()
        self.material = elastic.IsotropicMaterial(young_modulus=1.0, poisson_ratio=0.3)
        self.provider = GlobalAffineShapeProvider(self.nodes)

    def assemble(self):
        return elastic.assemble_elastic_operator(
            self.nodes,
            self.qpoints,
            self.volumes,
            self.material,
            self.provider,
        )

    def test_fixture_size_and_galerkin_symmetry_six_rigid_modes(self) -> None:
        self.assertEqual(self.nodes.shape, (27, 3))
        self.assertEqual(self.qpoints.shape, (64, 3))
        self.assertAlmostEqual(float(np.sum(self.volumes)), 8.0)
        result = self.assemble()
        self.assertIs(result.status, elastic.AssemblyStatus.OK)
        self.assertEqual(result.stiffness.shape, (81, 81))
        stiffness = result.stiffness
        self.assertIsNotNone(stiffness)

        relative_asymmetry = np.linalg.norm(stiffness - stiffness.T, ord="fro") / max(
            np.linalg.norm(stiffness, ord="fro"), 1.0
        )
        self.assertLessEqual(relative_asymmetry, 1.0e-12)

        modes = elastic.rigid_mode_coefficient_vectors(self.nodes)
        self.assertEqual(modes.shape, (81, 6))
        self.assertEqual(np.linalg.matrix_rank(modes), 6)
        scaled_residual = np.linalg.norm(stiffness @ modes, ord="fro") / max(
            np.linalg.norm(stiffness, ord="fro") * np.linalg.norm(modes, ord="fro"), 1.0
        )
        self.assertLessEqual(scaled_residual, 1.0e-10)

    def test_affine_constant_strain_energy_matches_independent_tensor_energy(self) -> None:
        result = self.assemble()
        self.assertIs(result.status, elastic.AssemblyStatus.OK)
        gradient = np.array(
            [[0.13, -0.07, 0.05], [0.11, -0.04, 0.09], [-0.06, 0.08, 0.17]],
            dtype=np.float64,
        )
        offset = np.array([0.2, -0.3, 0.1], dtype=np.float64)
        coefficients = (self.nodes @ gradient.T + offset).reshape(-1)
        computed = 0.5 * float(coefficients @ result.stiffness @ coefficients)

        strain = 0.5 * (gradient + gradient.T)
        lame_lambda = self.material.lame_lambda
        mu = self.material.shear_modulus
        stress = 2.0 * mu * strain + lame_lambda * np.trace(strain) * np.eye(3)
        independent = 0.5 * 8.0 * float(np.sum(strain * stress))
        self.assertGreater(independent, 0.0)
        self.assertLessEqual(abs(computed - independent) / independent, 1.0e-10)

        sample = self.provider.evaluate(np.array([0.23, -0.19, 0.31]))
        field_row = elastic.physical_displacement_row(sample, self.nodes.shape[0])
        physical_displacement = field_row @ coefficients
        expected_displacement = gradient @ np.array([0.23, -0.19, 0.31]) + offset
        np.testing.assert_allclose(physical_displacement, expected_displacement, atol=1.0e-12)

    def test_material_constitutive_matrix_is_positive_and_stable_range_is_enforced(self) -> None:
        d = self.material.constitutive_matrix
        np.testing.assert_allclose(d, d.T, atol=0.0)
        self.assertTrue(np.all(np.linalg.eigvalsh(d) > 0.0))
        for young_modulus, poisson_ratio in ((0.0, 0.3), (-1.0, 0.3), (1.0, -1.0), (1.0, 0.5)):
            with self.subTest(E=young_modulus, nu=poisson_ratio):
                with self.assertRaises(ValueError):
                    elastic.IsotropicMaterial(young_modulus, poisson_ratio)

    def test_unknown_support_blocks_without_a_partial_matrix(self) -> None:
        provider = PoorSupportProvider(elastic.ShapeStatus.UNKNOWN_POOR_SUPPORT)
        result = elastic.assemble_elastic_operator(
            self.nodes, self.qpoints, self.volumes, self.material, provider
        )
        self.assertIs(result.status, elastic.AssemblyStatus.BLOCKED_UNKNOWN_SUPPORT)
        self.assertIsNone(result.stiffness)
        self.assertEqual(result.blocked_quadrature_index, 0)
        self.assertEqual(provider.calls, 1)

        reason = "UNKNOWN_ILL_CONDITIONED: cond=1.2e13"
        unknown = PoorSupportProvider(elastic.ShapeStatus.UNKNOWN, reason)
        result = elastic.assemble_elastic_operator(
            self.nodes, self.qpoints, self.volumes, self.material, unknown
        )
        self.assertIsNone(result.stiffness)
        self.assertEqual(result.diagnostic, reason)

    def test_generic_unknown_requires_preserved_provider_reason(self) -> None:
        provider = PoorSupportProvider(elastic.ShapeStatus.UNKNOWN)
        with self.assertRaises(elastic.ShapeProviderContractError):
            elastic.assemble_elastic_operator(
                self.nodes, self.qpoints, self.volumes, self.material, provider
            )

    def test_malformed_provider_rows_fail_closed(self) -> None:
        bad_gradient_shape = MalformedShapeProvider(self.nodes)
        with self.assertRaisesRegex(elastic.ShapeProviderContractError, "shape gradients"):
            elastic.assemble_elastic_operator(
                self.nodes, self.qpoints, self.volumes, self.material, bad_gradient_shape
            )

        duplicate_ids = MalformedShapeProvider(self.nodes, bad_ids=True)
        with self.assertRaisesRegex(elastic.ShapeProviderContractError, "canonical global support_ids"):
            elastic.assemble_elastic_operator(
                self.nodes, self.qpoints, self.volumes, self.material, duplicate_ids
            )

    def test_dof_and_quadrature_caps_reject_before_copies_or_provider_calls(self) -> None:
        too_many_nodes = np.zeros((126, 3), dtype=np.float64)
        with patch.object(elastic.np, "array", side_effect=AssertionError("unexpected copy")):
            with self.assertRaisesRegex(ValueError, "budgetRejected"):
                elastic.assemble_elastic_operator(
                    too_many_nodes,
                    self.qpoints,
                    self.volumes,
                    self.material,
                    ExplodingProvider(),
                )

        too_many_qpoints = np.zeros((257, 3), dtype=np.float64)
        too_many_volumes = np.ones(257, dtype=np.float64)
        with patch.object(elastic.np, "array", side_effect=AssertionError("unexpected copy")):
            with self.assertRaisesRegex(ValueError, "budgetRejected"):
                elastic.assemble_elastic_operator(
                    self.nodes,
                    too_many_qpoints,
                    too_many_volumes,
                    self.material,
                    ExplodingProvider(),
                )

    def test_ndarray_subclass_is_rejected_before_array_conversion(self) -> None:
        class HostileNdarray(np.ndarray):
            def __new__(cls, base: np.ndarray):
                return base.view(cls)

            def __array__(self, *args, **kwargs):
                raise AssertionError("array conversion must not be invoked")

        subclass_nodes = self.nodes.view(HostileNdarray)
        with self.assertRaisesRegex(ValueError, "plain numpy.ndarray"):
            elastic.assemble_elastic_operator(
                subclass_nodes,
                self.qpoints,
                self.volumes,
                self.material,
                ExplodingProvider(),
            )

    def test_nonfinite_input_cannot_return_an_overflowed_local_operator(self) -> None:
        provider = HugeGradientProvider(self.nodes.shape[0], 1.0e200)
        one_point = np.array([[0.0, 0.0, 0.0]], dtype=np.float64)
        huge_volume = np.array([1.0e308], dtype=np.float64)
        with self.assertRaises(elastic.ElasticOperatorNumericsError):
            elastic.assemble_elastic_operator(
                self.nodes,
                one_point,
                huge_volume,
                self.material,
                provider,
            )

    def test_finite_local_operators_cannot_overflow_global_accumulation(self) -> None:
        provider = HugeGradientProvider(self.nodes.shape[0], 1.0e152)
        two_points = np.array([[0.0, 0.0, 0.0], [0.25, 0.0, 0.0]], dtype=np.float64)
        large_volumes = np.array([1.0e4, 1.0e4], dtype=np.float64)
        with self.assertRaises(elastic.ElasticOperatorNumericsError):
            elastic.assemble_elastic_operator(
                self.nodes,
                two_points,
                large_volumes,
                self.material,
                provider,
            )

    def test_finite_lame_parameters_that_overflow_D_are_rejected(self) -> None:
        material = elastic.IsotropicMaterial(young_modulus=9.0e307, poisson_ratio=0.4)
        with self.assertRaises(elastic.ElasticOperatorNumericsError):
            _ = material.constitutive_matrix


if __name__ == "__main__":
    unittest.main()
