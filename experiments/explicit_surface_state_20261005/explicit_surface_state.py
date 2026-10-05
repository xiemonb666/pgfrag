"""Source-only, fail-closed explicit crack-surface state candidate.

This module is deliberately isolated from ``pgfrag`` production code.  It stores
reference-space patch geometry and side kinematics, but it has no solver hook,
mechanical coupling proof verifier, observed-feature path, or production label
writer.  Any adjacency emitted here is namespaced fixture/proxy data.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
from typing import Any

import numpy as np


SCHEMA_VERSION = "pgfrag.explicit-crack-surface-state.v1"
FIXTURE_PROXY_NAMESPACE = "fixture_proxy.explicit_surface_state.v1"
FIXTURE_RELATION_NAMESPACE = "fixture_proxy_relation.explicit_surface_state.v1"
MECHANICS_QUALIFICATION = "NOT_COUPLED"
PRODUCTION_PROMOTION = "BLOCKED_NOT_COUPLED"
MAX_PATCHES = 4096
MAX_VERTICES_PER_PATCH = 100_000
MAX_TRIANGLES_PER_PATCH = 200_000
MAX_INPUT_NDARRAY_BYTES = 8 * 1024 * 1024
MAX_TOTAL_PATCH_ARRAY_BYTES = 256 * 1024 * 1024
MAX_ASSIGNMENTS = 100_000
MAX_HISTORICAL_PROXY_RELATIONS = 100_000


class SurfaceStateError(ValueError):
    """Raised when geometry, patch topology, or lineage is invalid."""


class UnsupportedPromotionError(RuntimeError):
    """Raised because this source-only candidate cannot produce production GT."""


def _identifier(value: str, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise SurfaceStateError(f"{name} must be a non-empty string")
    return value


def _step(value: int, name: str) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise SurfaceStateError(f"{name} must be an integer step")
    value = int(value)
    if value < 0:
        raise SurfaceStateError(f"{name} must be non-negative")
    return value


def _bounded_tuple(values: Any, *, maximum: int, name: str) -> tuple[Any, ...]:
    if isinstance(values, (tuple, list)) and len(values) > maximum:
        raise SurfaceStateError(f"{name} count exceeds cap {maximum}")
    result: list[Any] = []
    try:
        iterator = iter(values)
    except TypeError as exc:
        raise SurfaceStateError(f"{name} must be iterable") from exc
    for value in iterator:
        if len(result) >= maximum:
            raise SurfaceStateError(f"{name} count exceeds cap {maximum}")
        result.append(value)
    return tuple(result)


def _strict_ndarray_2d(
    value: Any,
    *,
    shape_tail: tuple[int, ...],
    name: str,
    max_rows: int,
    min_rows: int,
    allowed_kinds: str,
    max_itemsize: int,
) -> np.ndarray:
    """Preflight an exact ndarray before any copy or array-like conversion."""

    if type(value) is not np.ndarray:
        raise SurfaceStateError(f"{name} must be a numpy.ndarray")
    if value.ndim != 2 or tuple(value.shape[1:]) != shape_tail:
        raise SurfaceStateError(f"{name} must have shape (N, {', '.join(map(str, shape_tail))})")
    rows = int(value.shape[0])
    if rows < min_rows:
        raise SurfaceStateError(f"{name} must contain at least {min_rows} rows")
    if rows > max_rows:
        raise SurfaceStateError(f"{name} row count exceeds cap {max_rows}")
    if value.dtype.kind not in allowed_kinds or value.dtype.itemsize > max_itemsize:
        raise SurfaceStateError(f"{name} has an unsupported input dtype")
    if value.nbytes > MAX_INPUT_NDARRAY_BYTES:
        raise SurfaceStateError(f"{name} input bytes exceed cap {MAX_INPUT_NDARRAY_BYTES}")
    return value


def _frozen_array(value: Any, *, dtype: np.dtype, shape_tail: tuple[int, ...], name: str) -> np.ndarray:
    """Copy to canonical little-endian bytes, then expose a non-writeable view.

    A normal ``array.setflags(write=False)`` is reversible when the array owns
    writable memory.  An immutable ``bytes`` backing store makes the published
    ndarray read-only even if a caller tries to re-enable writes.
    """

    source = _strict_ndarray_2d(
        value,
        shape_tail=shape_tail,
        name=name,
        max_rows=MAX_VERTICES_PER_PATCH,
        min_rows=1,
        allowed_kinds="fiu",
        max_itemsize=16,
    )
    canonical_dtype = np.dtype(dtype).newbyteorder("<")
    try:
        array = np.ascontiguousarray(source, dtype=canonical_dtype)
    except (TypeError, ValueError, OverflowError) as exc:
        raise SurfaceStateError(f"{name} cannot be represented as {canonical_dtype}") from exc
    if array.size == 0:
        raise SurfaceStateError(f"{name} must not be empty")
    if canonical_dtype.kind == "f" and not np.isfinite(array).all():
        raise SurfaceStateError(f"{name} must contain only finite values")
    raw = array.tobytes(order="C")
    frozen = np.frombuffer(raw, dtype=canonical_dtype).reshape(array.shape)
    return frozen


def _int_array(value: Any, *, shape_tail: tuple[int, ...], name: str) -> np.ndarray:
    source = _strict_ndarray_2d(
        value,
        shape_tail=shape_tail,
        name=name,
        max_rows=MAX_TRIANGLES_PER_PATCH,
        min_rows=1,
        allowed_kinds="iu",
        max_itemsize=8,
    )
    if source.dtype.kind == "u" and np.any(source > np.iinfo(np.int64).max):
        raise SurfaceStateError(f"{name} indices do not fit signed 64-bit integers")
    try:
        canonical = np.ascontiguousarray(source, dtype=np.dtype("<i8"))
    except (TypeError, ValueError, OverflowError) as exc:
        raise SurfaceStateError(f"{name} indices do not fit signed 64-bit integers") from exc
    return np.frombuffer(canonical.tobytes(order="C"), dtype=np.dtype("<i8")).reshape(canonical.shape)


@dataclass(frozen=True, slots=True)
class SurfaceSide:
    """One side's independent displacement field over the common reference mesh."""

    side_sign: int
    body_id: str
    displacement_m: np.ndarray
    fragment_id: str | None = None

    def __post_init__(self) -> None:
        if (
            isinstance(self.side_sign, (bool, np.bool_))
            or not isinstance(self.side_sign, (int, np.integer))
            or self.side_sign not in (-1, 1)
        ):
            raise SurfaceStateError("side_sign must be -1 or +1")
        object.__setattr__(self, "side_sign", int(self.side_sign))
        object.__setattr__(self, "body_id", _identifier(self.body_id, "body_id"))
        if self.fragment_id is not None:
            object.__setattr__(self, "fragment_id", _identifier(self.fragment_id, "fragment_id"))
        displacement = _frozen_array(
            self.displacement_m, dtype=np.dtype("<f8"), shape_tail=(3,), name="displacement_m"
        )
        if displacement.ndim != 2:
            raise SurfaceStateError("displacement_m must have shape (vertex_count, 3)")
        object.__setattr__(self, "displacement_m", displacement)


@dataclass(frozen=True, slots=True)
class CrackPatch:
    """An authoritative reference-space patch with two views of one geometry."""

    surface_id: str
    patch_id: str
    creation_step: int
    parent_patch_id: str | None
    reference_vertices_m: np.ndarray
    triangles: np.ndarray
    negative_side: SurfaceSide
    positive_side: SurfaceSide
    declare_front_unknown: bool = False
    front_edges: np.ndarray = field(init=False, repr=False)
    front_status: str = field(init=False)
    reference_area_m2: float = field(init=False)
    geometry_sha256: str = field(init=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "surface_id", _identifier(self.surface_id, "surface_id"))
        object.__setattr__(self, "patch_id", _identifier(self.patch_id, "patch_id"))
        object.__setattr__(self, "creation_step", _step(self.creation_step, "creation_step"))
        if self.parent_patch_id is not None:
            object.__setattr__(self, "parent_patch_id", _identifier(self.parent_patch_id, "parent_patch_id"))
            if self.parent_patch_id == self.patch_id:
                raise SurfaceStateError("a patch cannot be its own parent")
        if not isinstance(self.negative_side, SurfaceSide) or self.negative_side.side_sign != -1:
            raise SurfaceStateError("negative_side must have side_sign=-1")
        if not isinstance(self.positive_side, SurfaceSide) or self.positive_side.side_sign != 1:
            raise SurfaceStateError("positive_side must have side_sign=+1")
        if not isinstance(self.declare_front_unknown, (bool, np.bool_)):
            raise SurfaceStateError("declare_front_unknown must be boolean")

        raw_vertices = _strict_ndarray_2d(
            self.reference_vertices_m,
            shape_tail=(3,),
            name="reference_vertices_m",
            max_rows=MAX_VERTICES_PER_PATCH,
            min_rows=3,
            allowed_kinds="fiu",
            max_itemsize=16,
        )
        raw_triangles = _strict_ndarray_2d(
            self.triangles,
            shape_tail=(3,),
            name="triangles",
            max_rows=MAX_TRIANGLES_PER_PATCH,
            min_rows=1,
            allowed_kinds="iu",
            max_itemsize=8,
        )
        vertices = _frozen_array(
            raw_vertices, dtype=np.dtype("<f8"), shape_tail=(3,), name="reference_vertices_m"
        )
        triangles = _int_array(raw_triangles, shape_tail=(3,), name="triangles")
        if vertices.ndim != 2 or vertices.shape[0] < 3:
            raise SurfaceStateError("reference_vertices_m must have shape (N, 3), N >= 3")
        if triangles.ndim != 2 or triangles.shape[0] < 1:
            raise SurfaceStateError("triangles must have shape (M, 3), M >= 1")
        if self.negative_side.displacement_m.shape != vertices.shape:
            raise SurfaceStateError("negative-side displacement must match the common reference vertices")
        if self.positive_side.displacement_m.shape != vertices.shape:
            raise SurfaceStateError("positive-side displacement must match the common reference vertices")
        if np.any(triangles < 0) or np.any(triangles >= vertices.shape[0]):
            raise SurfaceStateError("triangle index is outside the reference vertex array")
        if np.any(np.sort(triangles, axis=1)[:, 1:] == np.sort(triangles, axis=1)[:, :-1]):
            raise SurfaceStateError("a triangle repeats a vertex index")
        canonical_faces = np.sort(triangles, axis=1)
        if np.unique(canonical_faces, axis=0).shape[0] != triangles.shape[0]:
            raise SurfaceStateError("duplicate triangles are not allowed, including reverse-wound duplicates")
        if np.unique(vertices, axis=0).shape[0] != vertices.shape[0]:
            raise SurfaceStateError("duplicate reference vertex coordinates are ambiguous and not allowed")

        extent = float(np.max(np.ptp(vertices, axis=0)))
        scale2 = max(extent * extent, np.finfo(np.float64).tiny)
        cross_vectors = np.cross(
            vertices[triangles[:, 1]] - vertices[triangles[:, 0]],
            vertices[triangles[:, 2]] - vertices[triangles[:, 0]],
        )
        double_areas = np.linalg.norm(cross_vectors, axis=1)
        if not np.isfinite(double_areas).all() or np.any(double_areas <= 64.0 * np.finfo(np.float64).eps * scale2):
            raise SurfaceStateError("degenerate or numerically unresolved reference triangle")

        edge_uses: dict[tuple[int, int], list[tuple[int, int]]] = {}
        for a, b, c in triangles.tolist():
            for left, right in ((a, b), (b, c), (c, a)):
                key = (left, right) if left < right else (right, left)
                edge_uses.setdefault(key, []).append((left, right))
        boundary: list[tuple[int, int]] = []
        for key, uses in edge_uses.items():
            if len(uses) > 2:
                raise SurfaceStateError(f"nonmanifold edge {key} has {len(uses)} incident triangles")
            if len(uses) == 1:
                boundary.append(uses[0])
            elif uses[0] != (uses[1][1], uses[1][0]):
                raise SurfaceStateError(f"interior edge {key} has inconsistent triangle winding")
        if not boundary:
            front_status = "UNKNOWN_NO_ORIENTED_FRONT"
        else:
            incoming: dict[int, int] = {}
            outgoing: dict[int, int] = {}
            undirected: dict[int, set[int]] = {}
            for left, right in boundary:
                outgoing[left] = outgoing.get(left, 0) + 1
                incoming[right] = incoming.get(right, 0) + 1
                undirected.setdefault(left, set()).add(right)
                undirected.setdefault(right, set()).add(left)
            boundary_vertices = set(incoming) | set(outgoing)
            is_single_oriented_cycle = all(
                incoming.get(vertex, 0) == 1 and outgoing.get(vertex, 0) == 1
                for vertex in boundary_vertices
            )
            unseen = set(boundary_vertices)
            components = 0
            while unseen:
                components += 1
                stack = [unseen.pop()]
                while stack:
                    current = stack.pop()
                    for neighbor in undirected[current]:
                        if neighbor in unseen:
                            unseen.remove(neighbor)
                            stack.append(neighbor)
            front_status = (
                "SUPPORTED_SINGLE_ORIENTED_BOUNDARY_LOOP_TOPOLOGY"
                if is_single_oriented_cycle and components == 1 and not self.declare_front_unknown
                else "UNKNOWN_UNSUPPORTED_BRANCH_OR_INTERSECTION"
            )

        area = float(np.sum(double_areas) * 0.5)
        if not np.isfinite(area) or area <= 0.0:
            raise SurfaceStateError("reference patch must have finite positive area")
        geometry_sha = hashlib.sha256(
            b"pgfrag.explicit-crack-reference-geometry.v1\n"
            + vertices.tobytes(order="C")
            + b"--triangles--\n"
            + triangles.tobytes(order="C")
        ).hexdigest()
        frozen_boundary = np.asarray(boundary, dtype=np.dtype("<i8")).reshape((-1, 2))
        frozen_boundary = np.frombuffer(frozen_boundary.tobytes(), dtype=np.dtype("<i8")).reshape((-1, 2))

        object.__setattr__(self, "reference_vertices_m", vertices)
        object.__setattr__(self, "triangles", triangles)
        object.__setattr__(self, "front_edges", frozen_boundary)
        object.__setattr__(self, "front_status", front_status)
        object.__setattr__(self, "reference_area_m2", area)
        object.__setattr__(self, "geometry_sha256", geometry_sha)

    def deformed_vertices(self, side_sign: int) -> np.ndarray:
        if side_sign == -1:
            side = self.negative_side
        elif side_sign == 1:
            side = self.positive_side
        else:
            raise SurfaceStateError("side_sign must be -1 or +1")
        deformed = self.reference_vertices_m + side.displacement_m
        if not np.isfinite(deformed).all():
            raise SurfaceStateError("deformed side coordinates overflow to a non-finite value")
        raw = np.ascontiguousarray(deformed, dtype=np.dtype("<f8")).tobytes()
        return np.frombuffer(raw, dtype=np.dtype("<f8")).reshape(deformed.shape)

    def deformed_triangles(self, side_sign: int) -> np.ndarray:
        if side_sign == 1:
            triangles = self.triangles
        elif side_sign == -1:
            triangles = self.triangles[:, ::-1]
        else:
            raise SurfaceStateError("side_sign must be -1 or +1")
        return np.frombuffer(np.ascontiguousarray(triangles, dtype=np.dtype("<i8")).tobytes(), dtype=np.dtype("<i8")).reshape(triangles.shape)


@dataclass(frozen=True, slots=True)
class PatchAssignment:
    """Exact ownership of one patch by a pair of final fragment/body IDs."""

    patch_id: str
    negative_body_id: str
    positive_body_id: str
    negative_fragment_id: str
    positive_fragment_id: str
    separated_step: int

    def __post_init__(self) -> None:
        for field_name in (
            "patch_id",
            "negative_body_id",
            "positive_body_id",
            "negative_fragment_id",
            "positive_fragment_id",
        ):
            object.__setattr__(self, field_name, _identifier(getattr(self, field_name), field_name))
        object.__setattr__(self, "separated_step", _step(self.separated_step, "separated_step"))


@dataclass(frozen=True, slots=True)
class ProxyRelation:
    """Assigned fixture relation metadata; never a geometric or production/model R2 proof."""

    relation_namespace: str
    relation_kind: str
    surface_id: str
    patch_id: str
    original_patch_id: str
    negative_fragment_id: str
    positive_fragment_id: str
    separated_step: int
    shared_reference_area_m2: float

    def __post_init__(self) -> None:
        if self.relation_namespace != FIXTURE_RELATION_NAMESPACE:
            raise SurfaceStateError("proxy relations must remain in the fixture relation namespace")
        if self.relation_kind != "R2":
            raise SurfaceStateError("only the namespaced fixture relation kind R2 is supported")
        for field_name in (
            "surface_id",
            "patch_id",
            "original_patch_id",
            "negative_fragment_id",
            "positive_fragment_id",
        ):
            object.__setattr__(self, field_name, _identifier(getattr(self, field_name), field_name))
        if self.negative_fragment_id == self.positive_fragment_id:
            raise SurfaceStateError("a proxy relation requires distinct final fragments")
        object.__setattr__(self, "separated_step", _step(self.separated_step, "separated_step"))
        area = float(self.shared_reference_area_m2)
        if not np.isfinite(area) or area <= 0.0:
            raise SurfaceStateError("a proxy relation requires positive finite shared reference area")
        object.__setattr__(self, "shared_reference_area_m2", area)


@dataclass(frozen=True, slots=True)
class ExplicitSurfaceState:
    """Immutable bounded state and fixture-only historical relation ledger."""

    patches: tuple[CrackPatch, ...]
    assignments: tuple[PatchAssignment, ...] = ()
    historical_proxy_relations: tuple[ProxyRelation, ...] = ()
    current_visibility: tuple[tuple[str, bool], ...] = ()
    schema_version: str = field(default=SCHEMA_VERSION, init=False)
    namespace: str = field(default=FIXTURE_PROXY_NAMESPACE, init=False)
    mechanics_qualification: str = field(default=MECHANICS_QUALIFICATION, init=False)
    mechanically_active: bool = field(default=False, init=False)
    production_promotion: str = field(default=PRODUCTION_PROMOTION, init=False)
    surface_intersection_status: str = field(default="UNKNOWN_NOT_EVALUATED", init=False)

    def __post_init__(self) -> None:
        patches = _bounded_tuple(self.patches, maximum=MAX_PATCHES, name="patches")
        assignments = _bounded_tuple(self.assignments, maximum=MAX_ASSIGNMENTS, name="assignments")
        relations = _bounded_tuple(
            self.historical_proxy_relations,
            maximum=MAX_HISTORICAL_PROXY_RELATIONS,
            name="historical_proxy_relations",
        )
        visibility = _bounded_tuple(self.current_visibility, maximum=MAX_PATCHES, name="current_visibility")
        if not patches or any(not isinstance(patch, CrackPatch) for patch in patches):
            raise SurfaceStateError("state needs one or more validated CrackPatch values")
        if any(not isinstance(item, PatchAssignment) for item in assignments):
            raise SurfaceStateError("assignments must contain PatchAssignment values")
        if any(not isinstance(item, ProxyRelation) for item in relations):
            raise SurfaceStateError("historical_proxy_relations must contain ProxyRelation values")
        total_array_bytes = sum(
            patch.reference_vertices_m.nbytes
            + patch.triangles.nbytes
            + patch.front_edges.nbytes
            + patch.negative_side.displacement_m.nbytes
            + patch.positive_side.displacement_m.nbytes
            for patch in patches
        )
        if total_array_bytes > MAX_TOTAL_PATCH_ARRAY_BYTES:
            raise SurfaceStateError(
                f"state array bytes exceed cap {MAX_TOTAL_PATCH_ARRAY_BYTES}"
            )
        patch_by_id: dict[str, CrackPatch] = {}
        roots_by_surface: dict[str, int] = {}
        for patch in patches:
            if patch.patch_id in patch_by_id:
                raise SurfaceStateError(f"duplicate patch_id {patch.patch_id!r}")
            patch_by_id[patch.patch_id] = patch
            if patch.parent_patch_id is None:
                roots_by_surface[patch.surface_id] = roots_by_surface.get(patch.surface_id, 0) + 1
        for patch in patches:
            if patch.parent_patch_id is not None:
                parent = patch_by_id.get(patch.parent_patch_id)
                if parent is None:
                    raise SurfaceStateError(f"unknown parent patch {patch.parent_patch_id!r}")
                if parent.surface_id != patch.surface_id:
                    raise SurfaceStateError("parent and child patches must share one surface_id")
                if patch.creation_step <= parent.creation_step:
                    raise SurfaceStateError("child patch creation_step must follow its parent")
                visited = {patch.patch_id}
                ancestor = parent
                while ancestor is not None:
                    if ancestor.patch_id in visited:
                        raise SurfaceStateError("patch lineage contains a cycle")
                    visited.add(ancestor.patch_id)
                    ancestor = patch_by_id.get(ancestor.parent_patch_id) if ancestor.parent_patch_id else None
        if any(count != 1 for count in roots_by_surface.values()):
            raise SurfaceStateError("each surface_id must have exactly one root patch")
        if set(roots_by_surface) != {patch.surface_id for patch in patches}:
            raise SurfaceStateError("each surface_id must declare a root patch")

        assignment_by_patch: dict[str, PatchAssignment] = {}
        for assignment in assignments:
            if assignment.patch_id in assignment_by_patch:
                raise SurfaceStateError("a patch can have only one final ownership assignment")
            patch = patch_by_id.get(assignment.patch_id)
            if patch is None:
                raise SurfaceStateError(f"assignment references unknown patch {assignment.patch_id!r}")
            if assignment.separated_step < patch.creation_step:
                raise SurfaceStateError("separated_step cannot precede patch creation")
            assignment_by_patch[assignment.patch_id] = assignment
        related_patch_ids: set[str] = set()
        for relation in relations:
            patch = patch_by_id.get(relation.patch_id)
            assignment = assignment_by_patch.get(relation.patch_id)
            if patch is None or assignment is None:
                raise SurfaceStateError("historical proxy relation must reference an assigned patch")
            if relation.patch_id in related_patch_ids:
                raise SurfaceStateError("historical proxy relation may occur only once per patch")
            related_patch_ids.add(relation.patch_id)
            if relation.surface_id != patch.surface_id:
                raise SurfaceStateError("relation surface_id does not match its patch")
            if relation.original_patch_id != self._root_patch_id_in(patch_by_id, patch.patch_id):
                raise SurfaceStateError("relation original_patch_id does not match patch lineage")
            if relation.negative_fragment_id != assignment.negative_fragment_id or relation.positive_fragment_id != assignment.positive_fragment_id:
                raise SurfaceStateError("historical relation fragment pair does not match its patch assignment")
            if assignment.negative_body_id == assignment.positive_body_id:
                raise SurfaceStateError("same-body internal cracks cannot carry a fragment relation")
            if relation.shared_reference_area_m2 != patch.reference_area_m2:
                raise SurfaceStateError("historical relation area does not match authoritative patch geometry")
            if relation.separated_step != assignment.separated_step:
                raise SurfaceStateError("historical relation step does not match its patch assignment")
            if patch.front_status != "SUPPORTED_SINGLE_ORIENTED_BOUNDARY_LOOP_TOPOLOGY":
                raise SurfaceStateError("unknown front topology cannot carry a fragment relation")
        seen_visibility: set[str] = set()
        for patch_id, visible in visibility:
            if patch_id not in patch_by_id:
                raise SurfaceStateError(f"visibility references unknown patch {patch_id!r}")
            if patch_id in seen_visibility:
                raise SurfaceStateError("visibility may contain only one current state per patch")
            if not isinstance(visible, (bool, np.bool_)):
                raise SurfaceStateError("visibility state must be boolean")
            seen_visibility.add(patch_id)
        object.__setattr__(self, "patches", patches)
        object.__setattr__(self, "assignments", assignments)
        object.__setattr__(self, "historical_proxy_relations", relations)
        object.__setattr__(self, "current_visibility", tuple((pid, bool(value)) for pid, value in visibility))

    @property
    def patch_count(self) -> int:
        return len(self.patches)

    @property
    def triangle_count(self) -> int:
        return sum(int(patch.triangles.shape[0]) for patch in self.patches)

    @property
    def proxy_relation_count(self) -> int:
        return len(self.historical_proxy_relations)

    def _patch_map(self) -> dict[str, CrackPatch]:
        return {patch.patch_id: patch for patch in self.patches}

    @staticmethod
    def _root_patch_id_in(patch_by_id: dict[str, CrackPatch], patch_id: str) -> str:
        patch = patch_by_id[patch_id]
        while patch.parent_patch_id is not None:
            patch = patch_by_id[patch.parent_patch_id]
        return patch.patch_id

    def _leaf_patch_ids(self) -> set[str]:
        parents = {patch.parent_patch_id for patch in self.patches if patch.parent_patch_id is not None}
        return {patch.patch_id for patch in self.patches if patch.patch_id not in parents}

    def assign_patch(self, assignment: PatchAssignment) -> "ExplicitSurfaceState":
        """Record exact ownership for this patch; do not expand to descendants."""

        patch_by_id = self._patch_map()
        patch = patch_by_id.get(assignment.patch_id)
        if patch is None:
            raise SurfaceStateError(f"assignment references unknown patch {assignment.patch_id!r}")
        if any(item.patch_id == assignment.patch_id for item in self.assignments):
            raise SurfaceStateError("patch assignment is immutable; add a descendant patch for later splits")
        assignments = self.assignments + (assignment,)
        relations = self.historical_proxy_relations
        separated = (
            assignment.negative_body_id != assignment.positive_body_id
            and assignment.negative_fragment_id != assignment.positive_fragment_id
            and patch.front_status == "SUPPORTED_SINGLE_ORIENTED_BOUNDARY_LOOP_TOPOLOGY"
        )
        if separated:
            relations += (
                ProxyRelation(
                    relation_namespace=FIXTURE_RELATION_NAMESPACE,
                    relation_kind="R2",
                    surface_id=patch.surface_id,
                    patch_id=patch.patch_id,
                    original_patch_id=self.original_patch_id(patch.patch_id),
                    negative_fragment_id=assignment.negative_fragment_id,
                    positive_fragment_id=assignment.positive_fragment_id,
                    separated_step=assignment.separated_step,
                    shared_reference_area_m2=patch.reference_area_m2,
                ),
            )
        return ExplicitSurfaceState(
            patches=self.patches,
            assignments=assignments,
            historical_proxy_relations=relations,
            current_visibility=self.current_visibility,
        )

    def original_patch_id(self, patch_id: str) -> str:
        patch_by_id = self._patch_map()
        patch = patch_by_id.get(patch_id)
        if patch is None:
            raise SurfaceStateError(f"unknown patch_id {patch_id!r}")
        while patch.parent_patch_id is not None:
            patch = patch_by_id[patch.parent_patch_id]
        return patch.patch_id

    def with_current_visibility(self, patch_id: str, visible: bool) -> "ExplicitSurfaceState":
        """Update current compatibility visibility without editing history."""

        if patch_id not in self._patch_map():
            raise SurfaceStateError(f"visibility references unknown patch {patch_id!r}")
        if not isinstance(visible, (bool, np.bool_)):
            raise SurfaceStateError("visible must be boolean")
        visibility = dict(self.current_visibility)
        visibility[patch_id] = bool(visible)
        return ExplicitSurfaceState(
            patches=self.patches,
            assignments=self.assignments,
            historical_proxy_relations=self.historical_proxy_relations,
            current_visibility=tuple(sorted(visibility.items())),
        )

    def current_leaf_proxy_relations(self) -> tuple[ProxyRelation, ...]:
        """Return assigned metadata rows for leaf patch IDs, without geometry verification."""

        leaf_ids = self._leaf_patch_ids()
        historical_by_patch = {relation.patch_id: relation for relation in self.historical_proxy_relations}
        return tuple(
            historical_by_patch[patch_id]
            for patch_id in sorted(leaf_ids)
            if patch_id in historical_by_patch
        )

    def current_visible_leaf_proxy_relations(self) -> tuple[ProxyRelation, ...]:
        visible = dict(self.current_visibility)
        return tuple(
            relation
            for relation in self.current_leaf_proxy_relations()
            if visible.get(relation.patch_id, True)
        )

    def assert_production_eligible(self) -> None:
        """Fail closed: this module has no external solver proof verifier."""

        raise UnsupportedPromotionError(
            "NOT_COUPLED: no external solver proof binds this surface to displacement "
            "discontinuity, support visibility, and fracture evolution"
        )

    def to_payload(self) -> dict[str, Any]:
        """Serialize only the isolated fixture/proxy state and its own namespace."""

        patches = []
        for patch in self.patches:
            patches.append(
                {
                    "surface_id": patch.surface_id,
                    "patch_id": patch.patch_id,
                    "creation_step": patch.creation_step,
                    "parent_patch_id": patch.parent_patch_id,
                    "reference_vertices_m": patch.reference_vertices_m.tolist(),
                    "triangles": patch.triangles.tolist(),
                    "oriented_front_edges": patch.front_edges.tolist(),
                    "front_status": patch.front_status,
                    "reference_area_m2": patch.reference_area_m2,
                    "geometry_sha256": patch.geometry_sha256,
                    "negative_side": {
                        "side_sign": -1,
                        "body_id": patch.negative_side.body_id,
                        "fragment_id": patch.negative_side.fragment_id,
                        "displacement_m": patch.negative_side.displacement_m.tolist(),
                    },
                    "positive_side": {
                        "side_sign": 1,
                        "body_id": patch.positive_side.body_id,
                        "fragment_id": patch.positive_side.fragment_id,
                        "displacement_m": patch.positive_side.displacement_m.tolist(),
                    },
                }
            )
        return {
            "schema_version": self.schema_version,
            "namespace": self.namespace,
            "mechanics_qualification": self.mechanics_qualification,
            "mechanically_active": self.mechanically_active,
            "production_promotion": self.production_promotion,
            "surface_intersection_status": self.surface_intersection_status,
            "patches": patches,
            "patch_assignments": [
                {
                    "patch_id": item.patch_id,
                    "negative_body_id": item.negative_body_id,
                    "positive_body_id": item.positive_body_id,
                    "negative_fragment_id": item.negative_fragment_id,
                    "positive_fragment_id": item.positive_fragment_id,
                    "separated_step": item.separated_step,
                }
                for item in self.assignments
            ],
            "proxy_relation_namespace": FIXTURE_RELATION_NAMESPACE,
            "historical_proxy_relations": [
                {
                    "relation_kind": item.relation_kind,
                    "surface_id": item.surface_id,
                    "patch_id": item.patch_id,
                    "original_patch_id": item.original_patch_id,
                    "negative_fragment_id": item.negative_fragment_id,
                    "positive_fragment_id": item.positive_fragment_id,
                    "separated_step": item.separated_step,
                    "shared_reference_area_m2": item.shared_reference_area_m2,
                }
                for item in self.historical_proxy_relations
            ],
            "current_visibility": [
                {"patch_id": patch_id, "visible": visible}
                for patch_id, visible in self.current_visibility
            ],
            "counts": {
                "patch_count": self.patch_count,
                "triangle_count": self.triangle_count,
                "historical_proxy_relation_count": self.proxy_relation_count,
                "mechanically_active_count": 0,
            },
        }

    def to_json(self) -> str:
        return json.dumps(self.to_payload(), sort_keys=True, separators=(",", ":"), allow_nan=False)


def validate_surface_state(state: ExplicitSurfaceState) -> None:
    """Public validation entry point for callers loading composed states."""

    if not isinstance(state, ExplicitSurfaceState):
        raise SurfaceStateError("state must be ExplicitSurfaceState")
    # Construction is the validator; the frozen instance has already passed it.
    # Reconstructing also protects integrations that use object.__new__ manually.
    ExplicitSurfaceState(
        patches=state.patches,
        assignments=state.assignments,
        historical_proxy_relations=state.historical_proxy_relations,
        current_visibility=state.current_visibility,
    )

