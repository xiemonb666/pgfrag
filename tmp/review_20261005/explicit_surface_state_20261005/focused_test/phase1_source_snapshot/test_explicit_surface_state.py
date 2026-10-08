"""Counterexample-focused candidate tests. Intentionally not run in this handoff."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest
from unittest.mock import patch as mock_patch

import numpy as np


MODULE_PATH = Path(__file__).with_name("explicit_surface_state.py")
SPEC = importlib.util.spec_from_file_location("pgfrag_explicit_surface_candidate", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
surface = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = surface
SPEC.loader.exec_module(surface)


class ArrayLikeTrap:
    def __init__(self):
        self.conversion_calls = 0

    def __array__(self, *args, **kwargs):
        self.conversion_calls += 1
        raise AssertionError("array-like conversion must not run")


def _patch(
    patch_id="root",
    *,
    parent_patch_id=None,
    creation_step=0,
    surface_id="surface-0",
    body_negative="body-0",
    body_positive="body-0",
    negative_displacement=None,
    positive_displacement=None,
    declare_front_unknown=False,
    vertices=None,
    triangles=None,
):
    if vertices is None:
        vertices = np.asarray(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
            dtype=np.float64,
        )
    if triangles is None:
        triangles = np.asarray([[0, 1, 2], [0, 2, 3]], dtype=np.int64)
    if negative_displacement is None:
        negative_displacement = np.zeros_like(vertices, dtype=np.float64)
    if positive_displacement is None:
        positive_displacement = np.zeros_like(vertices, dtype=np.float64)
    return surface.CrackPatch(
        surface_id=surface_id,
        patch_id=patch_id,
        creation_step=creation_step,
        parent_patch_id=parent_patch_id,
        reference_vertices_m=vertices,
        triangles=triangles,
        negative_side=surface.SurfaceSide(-1, body_negative, negative_displacement),
        positive_side=surface.SurfaceSide(1, body_positive, positive_displacement),
        declare_front_unknown=declare_front_unknown,
    )


def _separated(state, patch_id, negative_fragment, positive_fragment, *, negative_body=None, positive_body=None, step=1):
    if negative_body is None:
        negative_body = f"body-{negative_fragment}"
    if positive_body is None:
        positive_body = f"body-{positive_fragment}"
    return state.assign_patch(
        surface.PatchAssignment(
            patch_id=patch_id,
            negative_body_id=negative_body,
            positive_body_id=positive_body,
            negative_fragment_id=negative_fragment,
            positive_fragment_id=positive_fragment,
            separated_step=step,
        )
    )


class ExplicitSurfaceStateTests(unittest.TestCase):
    def test_array_like_and_fake_array_protocol_are_rejected_without_conversion(self):
        fake_displacement = ArrayLikeTrap()
        with self.assertRaisesRegex(surface.SurfaceStateError, "must be a numpy.ndarray"):
            surface.SurfaceSide(-1, "body-0", fake_displacement)
        self.assertEqual(fake_displacement.conversion_calls, 0)

        fake_vertices = ArrayLikeTrap()
        zero_displacement = np.zeros((3, 3), dtype=np.float64)
        with self.assertRaisesRegex(surface.SurfaceStateError, "must be a numpy.ndarray"):
            _patch(
                vertices=fake_vertices,
                triangles=np.asarray([[0, 1, 2]], dtype=np.int64),
                negative_displacement=zero_displacement,
                positive_displacement=zero_displacement,
            )
        self.assertEqual(fake_vertices.conversion_calls, 0)

        fake_triangles = ArrayLikeTrap()
        with self.assertRaisesRegex(surface.SurfaceStateError, "must be a numpy.ndarray"):
            _patch(
                vertices=np.asarray([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
                triangles=fake_triangles,
                negative_displacement=zero_displacement,
                positive_displacement=zero_displacement,
            )
        self.assertEqual(fake_triangles.conversion_calls, 0)

    def test_oversized_zero_stride_views_are_rejected_before_contiguous_copy(self):
        vertex_source = np.asarray([0.0, 0.0, 0.0], dtype=np.float64)
        huge_displacements = np.broadcast_to(
            vertex_source, (surface.MAX_VERTICES_PER_PATCH + 1, 3)
        )
        with mock_patch.object(
            surface.np,
            "ascontiguousarray",
            side_effect=AssertionError("copy attempted before row cap"),
        ):
            with self.assertRaisesRegex(surface.SurfaceStateError, "row count exceeds cap"):
                surface.SurfaceSide(-1, "body-0", huge_displacements)

        face_source = np.asarray([[0, 1, 2]], dtype=np.int64)
        huge_faces = np.broadcast_to(face_source, (surface.MAX_TRIANGLES_PER_PATCH + 1, 3))
        with mock_patch.object(
            surface.np,
            "ascontiguousarray",
            side_effect=AssertionError("copy attempted before triangle cap"),
        ):
            with self.assertRaisesRegex(surface.SurfaceStateError, "row count exceeds cap"):
                surface._int_array(huge_faces, shape_tail=(3,), name="triangles")

    def test_arrays_are_defensively_copied_and_backed_by_immutable_bytes(self):
        vertices = np.asarray(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
            dtype=np.float64,
        )
        triangles = np.asarray([[0, 1, 2], [0, 2, 3]], dtype=np.int64)
        patch = _patch(vertices=vertices, triangles=triangles)
        vertices[:] = 999.0
        triangles[:] = 0
        self.assertEqual(patch.reference_area_m2, 1.0)
        self.assertEqual(patch.triangles.tolist(), [[0, 1, 2], [0, 2, 3]])
        with self.assertRaises(ValueError):
            patch.reference_vertices_m.setflags(write=True)

    def test_duplicate_triangle_reversed_winding_is_rejected(self):
        vertices = np.asarray(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=np.float64
        )
        with self.assertRaisesRegex(surface.SurfaceStateError, "duplicate triangles"):
            _patch(vertices=vertices, triangles=[[0, 1, 2], [2, 1, 0]])

    def test_zero_area_and_nearly_degenerate_triangles_are_rejected(self):
        collinear = np.asarray(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]], dtype=np.float64
        )
        with self.assertRaisesRegex(surface.SurfaceStateError, "degenerate"):
            _patch(vertices=collinear, triangles=[[0, 1, 2]])

    def test_three_faces_on_one_edge_are_nonmanifold_not_a_front_guess(self):
        vertices = np.asarray(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, -1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
            dtype=np.float64,
        )
        with self.assertRaisesRegex(surface.SurfaceStateError, "nonmanifold edge"):
            _patch(vertices=vertices, triangles=[[0, 1, 2], [1, 0, 3], [0, 1, 4]])

    def test_touching_front_cycles_remain_unknown_junction(self):
        # Two individually valid triangles touch only at one vertex. A naive
        # edge-count check sees no nonmanifold edge, but the combined front has
        # a multi-branch vertex and must remain explicitly unknown.
        vertices = np.asarray(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [-1.0, 0.0, 0.0],
                [0.0, -1.0, 0.0],
            ],
            dtype=np.float64,
        )
        patch = _patch(vertices=vertices, triangles=[[0, 1, 2], [0, 3, 4]])
        self.assertEqual(patch.front_status, "UNKNOWN_UNSUPPORTED_BRANCH_OR_INTERSECTION")
        state = surface.ExplicitSurfaceState((patch,))
        state = _separated(state, "root", "f-neg", "f-pos")
        self.assertEqual(state.historical_proxy_relations, ())

    def test_both_sides_share_reference_mesh_reverse_winding_but_not_displacement(self):
        negative = np.zeros((4, 3), dtype=np.float64)
        positive = np.zeros((4, 3), dtype=np.float64)
        negative[:, 2] = -0.1
        positive[:, 2] = 0.25
        patch = _patch(negative_displacement=negative, positive_displacement=positive)
        np.testing.assert_array_equal(patch.deformed_triangles(+1), patch.triangles)
        np.testing.assert_array_equal(patch.deformed_triangles(-1), patch.triangles[:, ::-1])
        np.testing.assert_array_equal(patch.deformed_vertices(+1)[:, :2], patch.reference_vertices_m[:, :2])
        np.testing.assert_array_equal(patch.deformed_vertices(-1)[:, :2], patch.reference_vertices_m[:, :2])
        self.assertFalse(np.array_equal(patch.deformed_vertices(+1), patch.deformed_vertices(-1)))
        self.assertEqual(patch.geometry_sha256, _patch().geometry_sha256)

    def test_same_body_internal_surface_never_creates_fragment_relation(self):
        state = surface.ExplicitSurfaceState((_patch(),))
        state = _separated(
            state,
            "root",
            "fragment-neg",
            "fragment-pos",
            negative_body="still-one-body",
            positive_body="still-one-body",
        )
        self.assertEqual(state.historical_proxy_relations, ())
        self.assertEqual(state.current_leaf_proxy_relations(), ())

    def test_manually_injected_relation_cannot_bypass_same_body_gate(self):
        patch = _patch()
        assignment = surface.PatchAssignment("root", "one-body", "one-body", "f-neg", "f-pos", 1)
        forged = surface.ProxyRelation(
            relation_namespace=surface.FIXTURE_RELATION_NAMESPACE,
            relation_kind="R2",
            surface_id="surface-0",
            patch_id="root",
            original_patch_id="root",
            negative_fragment_id="f-neg",
            positive_fragment_id="f-pos",
            separated_step=1,
            shared_reference_area_m2=patch.reference_area_m2,
        )
        with self.assertRaisesRegex(surface.SurfaceStateError, "same-body internal cracks"):
            surface.ExplicitSurfaceState((patch,), (assignment,), (forged,))

    def test_positive_area_distinct_body_and_fragment_pair_is_explicitly_proxy_only(self):
        state = surface.ExplicitSurfaceState((_patch(),))
        state = _separated(state, "root", "f-neg", "f-pos")
        relation = state.current_leaf_proxy_relations()[0]
        self.assertEqual(relation.relation_kind, "R2")
        self.assertEqual(relation.relation_namespace, surface.FIXTURE_RELATION_NAMESPACE)
        self.assertGreater(relation.shared_reference_area_m2, 0.0)
        self.assertEqual(relation.original_patch_id, "root")

    def test_wear_visibility_does_not_rewrite_historical_relation(self):
        state = surface.ExplicitSurfaceState((_patch(),))
        state = _separated(state, "root", "f-neg", "f-pos")
        historical_before = state.historical_proxy_relations
        hidden = state.with_current_visibility("root", False)
        self.assertEqual(hidden.historical_proxy_relations, historical_before)
        self.assertEqual(hidden.current_leaf_proxy_relations(), historical_before)
        self.assertEqual(hidden.current_visible_leaf_proxy_relations(), ())

    def test_descendant_pairs_follow_patch_lineage_without_cartesian_product(self):
        root = _patch("root")
        child_a = _patch("child-a", parent_patch_id="root", creation_step=2)
        child_b = _patch("child-b", parent_patch_id="root", creation_step=2)
        state = surface.ExplicitSurfaceState((root, child_a, child_b))
        state = _separated(state, "root", "old-neg", "old-pos", step=1)
        state = _separated(state, "child-a", "neg-a", "pos-a", step=3)
        state = _separated(state, "child-b", "neg-b", "pos-b", step=3)
        pairs = {
            (relation.patch_id, relation.negative_fragment_id, relation.positive_fragment_id)
            for relation in state.current_leaf_proxy_relations()
        }
        self.assertEqual(pairs, {("child-a", "neg-a", "pos-a"), ("child-b", "neg-b", "pos-b")})
        self.assertNotIn(("child-a", "neg-a", "pos-b"), pairs)
        self.assertNotIn(("child-b", "neg-b", "pos-a"), pairs)
        self.assertEqual(len(state.historical_proxy_relations), 3)

    def test_missing_parent_cross_surface_parent_and_nonincreasing_step_fail_closed(self):
        root = _patch("root")
        with self.assertRaisesRegex(surface.SurfaceStateError, "unknown parent"):
            surface.ExplicitSurfaceState((_patch("orphan", parent_patch_id="absent", creation_step=1),))
        child = _patch("child", parent_patch_id="root", creation_step=1, surface_id="other-surface")
        with self.assertRaisesRegex(surface.SurfaceStateError, "share one surface_id"):
            surface.ExplicitSurfaceState((root, child))
        with self.assertRaisesRegex(surface.SurfaceStateError, "must follow its parent"):
            surface.ExplicitSurfaceState((root, _patch("same-step", parent_patch_id="root", creation_step=0)))

    def test_reassignment_is_not_a_split_and_needs_a_descendant_patch(self):
        state = surface.ExplicitSurfaceState((_patch(),))
        state = state.assign_patch(
            surface.PatchAssignment("root", "b0", "b1", "f0", "f1", separated_step=1)
        )
        with self.assertRaisesRegex(surface.SurfaceStateError, "add a descendant patch"):
            state.assign_patch(surface.PatchAssignment("root", "b2", "b3", "f2", "f3", separated_step=2))

    def test_schema_and_serialization_have_no_production_label_leak(self):
        state = _separated(surface.ExplicitSurfaceState((_patch(),)), "root", "f-neg", "f-pos")
        payload = state.to_payload()
        encoded = state.to_json()
        self.assertEqual(payload["schema_version"], surface.SCHEMA_VERSION)
        self.assertEqual(payload["namespace"], surface.FIXTURE_PROXY_NAMESPACE)
        self.assertEqual(payload["mechanics_qualification"], "NOT_COUPLED")
        self.assertFalse(payload["mechanically_active"])
        self.assertEqual(payload["production_promotion"], "BLOCKED_NOT_COUPLED")
        self.assertEqual(payload["proxy_relation_namespace"], surface.FIXTURE_RELATION_NAMESPACE)
        self.assertNotIn('"R2":', encoded)
        self.assertNotIn('"model_r2":', encoded)
        self.assertNotIn("SIMULATED_GEOMETRIC_GT", encoded)
        self.assertNotIn("physical_r2", encoded.lower())
        self.assertEqual(payload["counts"]["mechanically_active_count"], 0)
        with self.assertRaises(surface.UnsupportedPromotionError):
            state.assert_production_eligible()

    def test_patch_count_budget_rejects_before_state_construction(self):
        original_cap = surface.MAX_PATCHES
        try:
            surface.MAX_PATCHES = 0
            with self.assertRaisesRegex(surface.SurfaceStateError, "patches count exceeds cap"):
                surface.ExplicitSurfaceState((_patch(),))
        finally:
            surface.MAX_PATCHES = original_cap

    def test_explicitly_unknown_front_cannot_emit_even_fixture_relation(self):
        state = surface.ExplicitSurfaceState((_patch(declare_front_unknown=True),))
        state = _separated(state, "root", "f-neg", "f-pos")
        self.assertEqual(state.historical_proxy_relations, ())


if __name__ == "__main__":
    unittest.main()
