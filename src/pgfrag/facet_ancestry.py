"""Exact source-facet ancestry for native tetrahedral meshes.

This module is deliberately combinatorial.  It binds current native T4
corners to frozen source T4 rows and local vertex ordinals, checks exact rest
coordinates, and pairs cohesive sides by original source node IDs.  It does
not infer ancestry from coordinates, inspect a fracture plane, or interpret
material state, quadrature state, or fracture physics.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import re
from types import MappingProxyType
from typing import Iterable, Mapping, Sequence

import numpy as np


# Face ordinals match pgfrag.cohesive_mesh.TET_FACES.  The order of each
# triplet is used only to name its stable local face ordinal; outward winding
# is computed from source rest coordinates below.
TET_FACES: tuple[tuple[int, int, int], ...] = (
    (0, 1, 2),
    (0, 1, 3),
    (0, 2, 3),
    (1, 2, 3),
)

FacetKey = tuple[int, int, int]
FacetIdentity = str

_SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")


@dataclass(frozen=True)
class FacetParent:
    """One source cell's outward-oriented incidence on a source facet."""

    source_cell_id: int
    source_cell_row: int
    local_face_ordinal: int
    outward_ordered_source_corner_ids: tuple[int, int, int]


@dataclass(frozen=True)
class SourceFacet:
    """A source triangle keyed by original source node IDs."""

    identity: FacetIdentity
    source_node_ids: FacetKey
    classification: str  # ``exterior`` or ``internal``
    parents: tuple[FacetParent, ...]


@dataclass(frozen=True)
class SourceFacetTable:
    """Frozen source mesh ancestry and the explicitly eligible internal set."""

    source_anchor_hash: str
    source_nodes: np.ndarray = field(repr=False, compare=False)
    source_tets: np.ndarray = field(repr=False, compare=False)
    source_node_ids: tuple[int, ...]
    source_cell_ids: tuple[int, ...]
    facets_by_node_ids: Mapping[FacetKey, SourceFacet] = field(repr=False, compare=False)
    facets_by_identity: Mapping[FacetIdentity, SourceFacet] = field(repr=False, compare=False)
    eligible_internal_facets: frozenset[FacetIdentity]

    @property
    def internal_facets(self) -> tuple[SourceFacet, ...]:
        """Return every source internal facet in stable identity order."""
        return tuple(
            self.facets_by_identity[key]
            for key in sorted(
                identity
                for identity, facet in self.facets_by_identity.items()
                if facet.classification == "internal"
            )
        )

    def eligible_facet_identities(self) -> tuple[FacetIdentity, ...]:
        """Stable, sorted eligible identities for snapshot reporting."""
        return tuple(sorted(self.eligible_internal_facets))


def make_facet_identity(source_anchor_hash: str, source_node_ids: Iterable[int]) -> FacetIdentity:
    """Build a stable identity from a frozen SHA-256 anchor and three source IDs.

    The canonical form is ``<sha256>|<id0>,<id1>,<id2>`` with the source IDs
    sorted numerically.  Coordinates and native row numbers are intentionally
    absent, so coincident but separately identified source nodes stay distinct.
    """
    anchor = _normalize_anchor_hash(source_anchor_hash)
    ids = _coerce_int_values(source_node_ids, label="facet source node IDs")
    if len(ids) != 3 or len(set(ids)) != 3:
        raise ValueError("a facet identity requires exactly three distinct source node IDs")
    return f"{anchor}|{','.join(str(value) for value in sorted(ids))}"


def build_source_facet_table(
    source_nodes: np.ndarray,
    source_tets: np.ndarray,
    source_node_ids: Sequence[int],
    source_cell_ids: Sequence[int],
    source_anchor_hash: str,
    *,
    eligible_internal_facets: Iterable[FacetIdentity | Sequence[int]] | None = None,
) -> SourceFacetTable:
    """Build source facet incidence from finite coordinates and positive T4s.

    ``source_tets`` contains zero-based rows into ``source_nodes``.  Its row
    order and each row's local corner order are the frozen lineage anchor used
    later for native mapping.  ``source_node_ids`` and ``source_cell_ids`` are
    unique externally meaningful IDs (for example, original MSH tags and cell
    tags).  ``source_anchor_hash`` must be the caller's frozen 64-hex SHA-256
    of the source artifact or source-anchor record.

    With ``eligible_internal_facets=None``, all source internal facets are
    eligible.  Otherwise the iterable may contain stable identity strings or
    three-node ID sequences.  An exterior or unknown facet is rejected rather
    than silently ignored.
    """
    anchor = _normalize_anchor_hash(source_anchor_hash)
    nodes = _finite_coordinates(source_nodes, label="source nodes")
    tets = _integer_matrix(source_tets, columns=4, label="source tetrahedra")
    if len(nodes) == 0:
        raise ValueError("source mesh must contain at least one node")
    if len(tets) == 0:
        raise ValueError("source mesh must contain at least one tetrahedron")
    if np.any(tets < 0) or np.any(tets >= len(nodes)):
        raise ValueError("source tetrahedron connectivity references an absent source node row")
    node_ids = _unique_ids(source_node_ids, expected=len(nodes), label="source node IDs")
    cell_ids = _unique_ids(source_cell_ids, expected=len(tets), label="source cell IDs")

    # Keep private immutable copies so later caller edits cannot change the
    # table's source lineage without changing its frozen hash anchor.
    nodes = np.array(nodes, dtype=np.float64, copy=True)
    tets = np.array(tets, dtype=np.int64, copy=True)
    nodes.setflags(write=False)
    tets.setflags(write=False)

    incidence: dict[FacetKey, list[FacetParent]] = {}
    for cell_row, tet in enumerate(tets):
        if len(set(map(int, tet))) != 4:
            raise ValueError(f"source tetrahedron row {cell_row} has a repeated node row")
        points = nodes[tet]
        volume6 = _signed_six_volume(points)
        if not np.isfinite(volume6) or volume6 <= 0.0:
            raise ValueError(
                f"source tetrahedron row {cell_row} must have finite positive rest orientation and nonzero volume"
            )
        for local_face_ordinal, face_ordinals in enumerate(TET_FACES):
            ordered_ordinals = list(face_ordinals)
            p0, p1, p2 = points[ordered_ordinals]
            opposite_ordinal = next(index for index in range(4) if index not in face_ordinals)
            opposite = points[opposite_ordinal]
            normal = np.cross(p1 - p0, p2 - p0)
            inward_test = float(np.dot(normal, opposite - p0))
            if not np.isfinite(inward_test) or inward_test == 0.0:
                raise ValueError(f"source tetrahedron row {cell_row} has a degenerate face")
            if inward_test > 0.0:
                ordered_ordinals[1], ordered_ordinals[2] = ordered_ordinals[2], ordered_ordinals[1]
            ordered_ids = tuple(node_ids[int(tet[ordinal])] for ordinal in ordered_ordinals)
            key = tuple(sorted(ordered_ids))
            incidence.setdefault(key, []).append(
                FacetParent(
                    source_cell_id=cell_ids[cell_row],
                    source_cell_row=cell_row,
                    local_face_ordinal=local_face_ordinal,
                    outward_ordered_source_corner_ids=ordered_ids,
                )
            )

    for key, parents in incidence.items():
        if len(parents) > 2:
            raise ValueError(
                f"nonmanifold source facet {key} has {len(parents)} source parents"
            )
        if len(parents) == 2 and _winding_parity(
            parents[0].outward_ordered_source_corner_ids
        ) == _winding_parity(parents[1].outward_ordered_source_corner_ids):
            raise ValueError(f"internal source facet {key} parents do not have opposed outward windings")

    facets_by_node_ids: dict[FacetKey, SourceFacet] = {}
    facets_by_identity: dict[FacetIdentity, SourceFacet] = {}
    for key, parents_list in incidence.items():
        parents = tuple(sorted(parents_list, key=lambda parent: parent.source_cell_row))
        identity = make_facet_identity(anchor, key)
        facet = SourceFacet(
            identity=identity,
            source_node_ids=key,
            classification="internal" if len(parents) == 2 else "exterior",
            parents=parents,
        )
        facets_by_node_ids[key] = facet
        facets_by_identity[identity] = facet

    if eligible_internal_facets is None:
        eligible = {
            facet.identity
            for facet in facets_by_identity.values()
            if facet.classification == "internal"
        }
    else:
        eligible = set()
        for requested in eligible_internal_facets:
            if isinstance(requested, str):
                _parse_facet_identity(requested)
                identity = requested.lower().split("|", 1)[0] + "|" + requested.split("|", 1)[1]
                facet = facets_by_identity.get(identity)
            else:
                ids = _coerce_int_values(requested, label="eligible facet source node IDs")
                if len(ids) != 3 or len(set(ids)) != 3:
                    raise ValueError("each eligible facet must contain exactly three distinct source node IDs")
                facet = facets_by_node_ids.get(tuple(sorted(ids)))
            if facet is None:
                raise ValueError(f"eligible source facet {requested!r} is absent from the source mesh")
            if facet.classification != "internal":
                raise ValueError(f"eligible source facet {facet.identity} is exterior, not internal")
            if facet.identity in eligible:
                raise ValueError(f"duplicate eligible source facet {facet.identity}")
            eligible.add(facet.identity)

    return SourceFacetTable(
        source_anchor_hash=anchor,
        source_nodes=nodes,
        source_tets=tets,
        source_node_ids=node_ids,
        source_cell_ids=cell_ids,
        facets_by_node_ids=MappingProxyType(facets_by_node_ids),
        facets_by_identity=MappingProxyType(facets_by_identity),
        eligible_internal_facets=frozenset(eligible),
    )


def map_native_cohesive_ancestry(
    source: SourceFacetTable,
    native_rest_nodes: np.ndarray,
    native_bulk_tets: np.ndarray,
    native_cohesive: np.ndarray,
    *,
    native_initial_bulk_tets: np.ndarray,
) -> dict[str, object]:
    """Map a current native cohesive snapshot to exact eligible source facets.

    The complete ``native_initial_bulk_tets`` snapshot must exactly equal the
    frozen source T4 connectivity array, including row and local-corner order.
    Current bulk rows must preserve that same shape/order/local-corner lineage;
    only node IDs may differ through native node duplication.  Each current
    bulk corner is mapped by its cell row and local ordinal, and its rest
    coordinate must exactly equal the corresponding source coordinate.

    Empty cohesive connectivity is valid and reports all eligible facets as
    missing from this snapshot.  Stable facet identity is independent of
    cohesive row ordering.  Output rows are sorted by stable identity, while
    ``native_cohesive_row`` is only a snapshot-local row number.  This helper
    makes no QP ordinal, material-state, fracture, source-geometry, or solver
    behavior claim.
    """
    rest = _finite_coordinates(native_rest_nodes, label="native rest nodes")
    bulk = _integer_matrix(native_bulk_tets, columns=4, label="native bulk tetrahedra")
    initial = _integer_matrix(native_initial_bulk_tets, columns=4, label="native initial bulk tetrahedra")
    cohesive = _integer_matrix(native_cohesive, columns=6, label="native cohesive connectivity")

    if rest.ndim != 2 or rest.shape[1] != 3:
        raise ValueError(f"native rest node array must have shape (n, 3), got {rest.shape}")
    if initial.shape != source.source_tets.shape:
        raise ValueError(
            "native initial/source bulk shape mismatch: "
            f"{initial.shape} != {source.source_tets.shape}"
        )
    if not np.array_equal(initial, source.source_tets):
        raise ValueError("native initial bulk connectivity differs from frozen source row/local-corner connectivity")
    if bulk.shape != source.source_tets.shape:
        raise ValueError(f"current native/source bulk shape mismatch: {bulk.shape} != {source.source_tets.shape}")
    if len(cohesive.shape) != 2 or cohesive.shape[1] != 6:
        raise ValueError(f"native cohesive connectivity must have shape (n, 6), got {cohesive.shape}")
    if np.any(bulk < 0) or np.any(bulk >= len(rest)):
        raise ValueError("current native bulk connectivity references an absent native node")
    if np.any(cohesive < 0) or np.any(cohesive >= len(rest)):
        raise ValueError("native cohesive connectivity references an absent native node")

    native_to_source: dict[int, int] = {}
    for cell_row, (source_tet, native_tet) in enumerate(zip(source.source_tets, bulk, strict=True)):
        if len(set(map(int, native_tet))) != 4:
            raise ValueError(f"current native bulk row {cell_row} has a repeated node ID")
        for local_ordinal, (source_node_row, native_node) in enumerate(zip(source_tet, native_tet, strict=True)):
            source_node_row = int(source_node_row)
            native_node = int(native_node)
            source_id = source.source_node_ids[source_node_row]
            previous = native_to_source.setdefault(native_node, source_id)
            if previous != source_id:
                raise ValueError(f"native node {native_node} has conflicting source ancestry")
            if not np.array_equal(rest[native_node], source.source_nodes[source_node_row]):
                raise ValueError(
                    "source-cell/local-corner exact rest coordinate consistency failed "
                    f"at bulk row {cell_row}, local corner {local_ordinal}"
                )

    if set(native_to_source) != set(range(len(rest))):
        unmapped = sorted(set(range(len(rest))) - set(native_to_source))
        raise ValueError(f"native rest nodes lack source bulk-corner ancestry: {unmapped[:8]}")

    native_face_incidence: dict[tuple[int, int, int], list[tuple[int, int]]] = {}
    for cell_row, tet in enumerate(bulk):
        for local_face_ordinal, face_ordinals in enumerate(TET_FACES):
            native_key = tuple(sorted(int(tet[ordinal]) for ordinal in face_ordinals))
            native_face_incidence.setdefault(native_key, []).append((cell_row, local_face_ordinal))
    for key, parents in native_face_incidence.items():
        if len(parents) > 2:
            raise ValueError(f"nonmanifold current native bulk face {key} has {len(parents)} parents")

    records_by_identity: dict[FacetIdentity, dict[str, object]] = {}
    for cohesive_row, cohesive_element in enumerate(cohesive):
        sides = (tuple(map(int, cohesive_element[:3])), tuple(map(int, cohesive_element[3:])))
        if any(len(set(side)) != 3 for side in sides):
            raise ValueError(f"native cohesive row {cohesive_row} has a repeated side node")
        try:
            source_side_keys = tuple(
                tuple(sorted(native_to_source[node] for node in side))
                for side in sides
            )
        except KeyError as error:
            raise ValueError(
                f"native cohesive row {cohesive_row} references a node without bulk ancestry: {error.args[0]}"
            ) from error
        if source_side_keys[0] != source_side_keys[1]:
            raise ValueError(f"native cohesive row {cohesive_row} sides do not pair by exact source node IDs")
        source_key = source_side_keys[0]
        facet = source.facets_by_node_ids.get(source_key)
        if facet is None:
            raise ValueError(f"native cohesive row {cohesive_row} does not identify a unique source facet {source_key}")
        if facet.classification != "internal" or len(facet.parents) != 2:
            raise ValueError(f"native cohesive row {cohesive_row} maps to a source exterior facet {facet.identity}")
        if facet.identity not in source.eligible_internal_facets:
            raise ValueError(f"native cohesive row {cohesive_row} maps to ineligible source facet {facet.identity}")
        if facet.identity in records_by_identity:
            raise ValueError(f"source facet {facet.identity} has duplicate native cohesive births")

        side_parent_rows: list[int] = []
        side_parent_local_ordinals: list[int] = []
        for side_index, side in enumerate(sides):
            native_key = tuple(sorted(side))
            adjacent = native_face_incidence.get(native_key, [])
            if len(adjacent) != 1:
                raise ValueError(
                    f"native cohesive row {cohesive_row} side {side_index} has "
                    f"{len(adjacent)} current native bulk parents; expected exactly one"
                )
            cell_row, local_face_ordinal = adjacent[0]
            side_parent_rows.append(cell_row)
            side_parent_local_ordinals.append(local_face_ordinal)

        expected_parent_pairs = {
            (parent.source_cell_row, parent.local_face_ordinal)
            for parent in facet.parents
        }
        observed_parent_pairs = set(zip(side_parent_rows, side_parent_local_ordinals, strict=True))
        if len(side_parent_rows) != 2 or observed_parent_pairs != expected_parent_pairs:
            raise ValueError(
                f"native cohesive row {cohesive_row} sides do not have exactly the two original source parent cells"
            )

        records_by_identity[facet.identity] = {
            "facet_identity": facet.identity,
            "source_node_ids": list(facet.source_node_ids),
            "source_parent_cell_ids": [parent.source_cell_id for parent in facet.parents],
            "source_parent_cell_rows": [parent.source_cell_row for parent in facet.parents],
            "source_parent_local_face_ordinals": [parent.local_face_ordinal for parent in facet.parents],
            "source_parent_outward_corner_ids": [
                list(parent.outward_ordered_source_corner_ids) for parent in facet.parents
            ],
            "native_cohesive_row": cohesive_row,
            "native_cohesive_row_semantics": "snapshot-local getConnectivity row; no persistent global cohesive ID is implied",
            "side_parent_source_cell_ids": [source.source_cell_ids[row] for row in side_parent_rows],
            "side_parent_source_cell_rows": side_parent_rows,
            "side_parent_local_face_ordinals": side_parent_local_ordinals,
        }

    active_ids = sorted(records_by_identity)
    return {
        "source_anchor_hash": source.source_anchor_hash,
        "lineage_scope": (
            "Source ancestry is assigned from frozen source cell row/local-corner ordinals and exact rest-coordinate checks. "
            "Post-insertion native bulk row/local-corner identity is conditional on the native backend preserving that correspondence."
        ),
        "lineage_assumptions": [
            "The supplied pre-insertion native bulk connectivity exactly equals the frozen source row/local-corner connectivity.",
            "The native backend preserves source bulk row and local-corner identity after cohesive insertion; this API has no independent persistent native cell IDs.",
            "Exact rest-coordinate consistency checks corroborate the row/local-corner mapping but do not independently prove persistent cell identity for indistinguishable coincident cells.",
        ],
        "native_node_to_source_node_id": {
            str(native_node): source_id for native_node, source_id in sorted(native_to_source.items())
        },
        "active_facet_identities": active_ids,
        "missing_eligible_facet_identities": sorted(source.eligible_internal_facets - set(active_ids)),
        "active_facets": [records_by_identity[identity] for identity in active_ids],
        "eligible_internal_facet_count": len(source.eligible_internal_facets),
        "active_cohesive_count": len(active_ids),
        "qp_ordinal_stability_claimed": False,
        "interpretation": "connectivity ancestry only; not a fracture, material-state, QP-order, or solver claim",
    }


def event_delta(
    previous_active_facets: Iterable[FacetIdentity],
    current_active_facets: Iterable[FacetIdentity],
) -> dict[str, object]:
    """Return monotone insertions between two snapshots, rejecting deletion/reuse.

    Identities must be canonical and unique in each snapshot.  A transition
    that drops a previously active identity is rejected, so a later transition
    cannot silently reinterpret a deleted identity as a new event.  Repeated
    identities within a snapshot are also rejected as reused/duplicated IDs.
    """
    previous = _identity_set(previous_active_facets, label="previous active facets")
    current = _identity_set(current_active_facets, label="current active facets")
    anchors = {identity.split("|", 1)[0] for identity in previous | current}
    if len(anchors) > 1:
        raise ValueError("event snapshots use different source-anchor hashes")
    deleted = sorted(previous - current)
    if deleted:
        raise ValueError(f"source facet identities disappeared from the current snapshot: {deleted}")
    return {
        "source_anchor_hash": next(iter(anchors), None),
        "newly_active_facet_identities": sorted(current - previous),
        "persistent_active_facet_identities": sorted(previous & current),
        "previous_active_count": len(previous),
        "current_active_count": len(current),
    }


def _normalize_anchor_hash(value: str) -> str:
    if not isinstance(value, str) or not _SHA256_RE.fullmatch(value):
        raise ValueError("source_anchor_hash must be a frozen 64-character SHA-256 hex digest")
    return value.lower()


def _parse_facet_identity(identity: str) -> tuple[str, FacetKey]:
    if not isinstance(identity, str) or identity.count("|") != 1:
        raise ValueError("facet identity must have canonical '<sha256>|<id0>,<id1>,<id2>' form")
    anchor_text, node_text = identity.split("|", 1)
    anchor = _normalize_anchor_hash(anchor_text)
    pieces = node_text.split(",")
    if len(pieces) != 3:
        raise ValueError("facet identity must contain exactly three source node IDs")
    try:
        ids = tuple(int(piece, 10) for piece in pieces)
    except ValueError as error:
        raise ValueError("facet identity source node IDs must be canonical decimal integers") from error
    if any(str(value) != piece for value, piece in zip(ids, pieces, strict=True)):
        raise ValueError("facet identity source node IDs must use canonical decimal integers")
    if len(set(ids)) != 3 or tuple(sorted(ids)) != ids:
        raise ValueError("facet identity source node IDs must be distinct and sorted")
    return anchor, ids  # type: ignore[return-value]


def _identity_set(values: Iterable[FacetIdentity], *, label: str) -> set[FacetIdentity]:
    result: set[FacetIdentity] = set()
    for value in values:
        anchor, node_ids = _parse_facet_identity(value)
        canonical = f"{anchor}|{','.join(str(item) for item in node_ids)}"
        if value != canonical:
            raise ValueError(f"{label} contains a noncanonical facet identity")
        if canonical in result:
            raise ValueError(f"{label} contains a reused/duplicate source facet identity: {canonical}")
        result.add(canonical)
    return result


def _coerce_int_values(values: Iterable[int], *, label: str) -> tuple[int, ...]:
    try:
        items = tuple(values)
    except TypeError as error:
        raise ValueError(f"{label} must be an iterable of integers") from error
    result: list[int] = []
    for item in items:
        if isinstance(item, (bool, np.bool_)) or not isinstance(item, (int, np.integer)):
            raise ValueError(f"{label} must contain integers only")
        result.append(int(item))
    return tuple(result)


def _unique_ids(values: Sequence[int], *, expected: int, label: str) -> tuple[int, ...]:
    ids = _coerce_int_values(values, label=label)
    if len(ids) != expected:
        raise ValueError(f"{label} count mismatch: got {len(ids)}, expected {expected}")
    if len(set(ids)) != len(ids):
        raise ValueError(f"{label} must be unique")
    return ids


def _integer_matrix(values: np.ndarray, *, columns: int, label: str) -> np.ndarray:
    array = np.asarray(values)
    if array.ndim != 2 or array.shape[1] != columns:
        raise ValueError(f"{label} must have shape (n, {columns}), got {array.shape}")
    if array.dtype.kind not in "iu":
        raise ValueError(f"{label} connectivity must contain integers only")
    if array.dtype.kind == "u" and array.size and np.max(array) > np.iinfo(np.int64).max:
        raise ValueError(f"{label} connectivity contains an integer outside int64 range")
    return np.asarray(array, dtype=np.int64)


def _finite_coordinates(values: np.ndarray, *, label: str) -> np.ndarray:
    try:
        array = np.asarray(values, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be numeric coordinates") from error
    if array.ndim != 2 or array.shape[1] != 3:
        raise ValueError(f"{label} must have shape (n, 3), got {array.shape}")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{label} coordinates must be finite")
    return array


def _signed_six_volume(points: np.ndarray) -> float:
    return float(np.dot(points[1] - points[0], np.cross(points[2] - points[0], points[3] - points[0])))


def _winding_parity(ordered_ids: Sequence[int]) -> int:
    inversions = sum(
        1
        for left in range(3)
        for right in range(left + 1, 3)
        if ordered_ids[left] > ordered_ids[right]
    )
    return inversions % 2
