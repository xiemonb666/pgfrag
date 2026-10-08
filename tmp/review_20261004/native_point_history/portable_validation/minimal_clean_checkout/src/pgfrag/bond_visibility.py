"""Fail-closed static visibility classification for PD horizon candidates.

This module is intentionally not connected to the native solver. It accepts a
canonical STEP file, hashes the bytes actually read, and classifies each
undirected reference-coordinate candidate against all B-rep intersections and
the open intervals between them. It never estimates visibility from endpoint
membership or equal-distance samples.

Only ``IN_MATERIAL`` is a geometrically eligible candidate. ``OUT_OR_VOID``
and ``UNKNOWN`` remain distinct in the returned static topology ledger.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
import hashlib
import math
from numbers import Integral, Real
from pathlib import Path
import re
from typing import Any


BOND_VISIBILITY_SCHEMA = "pgfrag.bond-visibility.v1"
BOND_VISIBILITY_ALGORITHM = "pgfrag/canonical-step-all-hits-interval-occupancy/1.0.0"
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_MM_TO_M = 1.0e-3


@dataclass(frozen=True, slots=True)
class MaterialSolid:
    """One validated, closed source solid and its assigned material-domain ID."""

    solid_index: int
    material_domain_id: str | None
    shape: Any
    faces: tuple[Any, ...]


@dataclass(frozen=True, slots=True)
class CanonicalStepSource:
    """Loaded geometry and the evidence needed to interpret its coordinates."""

    source_path: str | None
    source_origin: str
    source_geometry_sha256: str | None
    expected_sha256: str | None
    source_length_unit: str | None
    frame: str | None
    source_accuracy_mm: float | None
    source_accuracy_basis: str | None
    shape: Any | None
    solids: tuple[MaterialSolid, ...]
    status: str
    unknown_reasons: tuple[str, ...]
    backend_versions: Mapping[str, str]


def _is_sha256(value: Any) -> bool:
    return isinstance(value, str) and _SHA256_RE.fullmatch(value) is not None


def _finite_real(value: Any) -> bool:
    if isinstance(value, bool) or not isinstance(value, Real):
        return False
    try:
        return math.isfinite(float(value))
    except (OverflowError, ValueError):
        return False


def _xyz_mm(value: Any, *, label: str) -> tuple[float, float, float]:
    if isinstance(value, (str, bytes, Mapping)):
        raise ValueError(f"{label} must be a three-component sequence in mm")
    try:
        component_count = len(value)
    except (TypeError, ValueError):
        component_count = -1
    if component_count != 3:
        raise ValueError(f"{label} must be a three-component sequence in mm")
    if not all(_finite_real(component) for component in value):
        raise ValueError(f"{label} must contain finite non-boolean numbers")
    return tuple(float(component) for component in value)  # type: ignore[return-value]


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _unknown_source(
    *,
    source_path: str | None,
    source_origin: str,
    source_geometry_sha256: str | None,
    expected_sha256: str | None,
    source_length_unit: str | None,
    frame: str | None,
    source_accuracy_mm: float | None,
    source_accuracy_basis: str | None,
    shape: Any | None,
    solids: tuple[MaterialSolid, ...] = (),
    reasons: Iterable[str],
    backend_versions: Mapping[str, str] | None = None,
) -> CanonicalStepSource:
    return CanonicalStepSource(
        source_path=source_path,
        source_origin=source_origin,
        source_geometry_sha256=source_geometry_sha256,
        expected_sha256=expected_sha256,
        source_length_unit=source_length_unit,
        frame=frame,
        source_accuracy_mm=source_accuracy_mm,
        source_accuracy_basis=source_accuracy_basis,
        shape=shape,
        solids=solids,
        status="UNKNOWN",
        unknown_reasons=tuple(sorted(set(reasons))),
        backend_versions=dict(backend_versions or {}),
    )


def _occ_topology(shape: Any) -> tuple[tuple[Any, ...], tuple[str, ...]]:
    """Return direct solid leaves; flag any open or non-solid leaf."""
    from OCP.BRepCheck import BRepCheck_Analyzer
    from OCP.TopAbs import TopAbs_COMPOUND, TopAbs_COMPSOLID, TopAbs_SHELL, TopAbs_SOLID
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS, TopoDS_Iterator

    solids: list[Any] = []
    errors: list[str] = []

    def walk(current: Any) -> None:
        shape_type = current.ShapeType()
        if shape_type == TopAbs_SOLID:
            solid = TopoDS.Solid_s(current)
            analyzer = BRepCheck_Analyzer(solid)
            if not analyzer.IsValid():
                errors.append("SOURCE_SOLID_INVALID_OR_NONMANIFOLD")
            shell_explorer = TopExp_Explorer(solid, TopAbs_SHELL)
            shell_count = 0
            while shell_explorer.More():
                shell_count += 1
                shell = TopoDS.Shell_s(shell_explorer.Current())
                if not shell.Closed():
                    errors.append("SOURCE_SHELL_OPEN")
                shell_explorer.Next()
            if shell_count == 0:
                errors.append("SOURCE_SOLID_HAS_NO_SHELL")
            solids.append(solid)
            return
        if shape_type in (TopAbs_COMPOUND, TopAbs_COMPSOLID):
            iterator = TopoDS_Iterator()
            iterator.Initialize(current)
            had_child = False
            while iterator.More():
                had_child = True
                walk(iterator.Value())
                iterator.Next()
            if not had_child:
                errors.append("SOURCE_COMPOSITE_EMPTY")
            return
        errors.append("SOURCE_CONTAINS_NON_SOLID_LEAF")

    if shape.IsNull():
        return (), ("SOURCE_SHAPE_NULL",)
    walk(shape)
    if not solids:
        errors.append("SOURCE_HAS_NO_SOLIDS")
    return tuple(solids), tuple(sorted(set(errors)))


def _solid_faces(shape: Any) -> tuple[Any, ...]:
    from OCP.TopAbs import TopAbs_FACE
    from OCP.TopExp import TopExp_Explorer

    faces: list[Any] = []
    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    while explorer.More():
        faces.append(explorer.Current())
        explorer.Next()
    return tuple(faces)


def _source_from_shape(
    shape: Any,
    *,
    source_length_unit: str,
    frame: str,
    source_accuracy_mm: float | None,
    source_accuracy_basis: str | None,
    material_domain_ids: Sequence[str] | None = None,
    backend_versions: Mapping[str, str] | None = None,
) -> CanonicalStepSource:
    """Inspect an OCC shape; direct calls are intended for in-memory tests."""
    reasons: list[str] = []
    if source_length_unit != "mm":
        reasons.append("SOURCE_LENGTH_UNIT_NOT_MM")
    if not isinstance(frame, str) or not frame.strip():
        reasons.append("SOURCE_FRAME_MISSING")
    if not _finite_real(source_accuracy_mm) or float(source_accuracy_mm) < 0.0:
        reasons.append("SOURCE_PRECISION_UNKNOWN")
        accuracy = None
    else:
        accuracy = float(source_accuracy_mm)
    if not isinstance(source_accuracy_basis, str) or not source_accuracy_basis.strip():
        reasons.append("SOURCE_PRECISION_BASIS_MISSING")

    try:
        raw_solids, topology_errors = _occ_topology(shape)
    except Exception:
        raw_solids, topology_errors = (), ("SOURCE_TOPOLOGY_INSPECTION_FAILED",)
    reasons.extend(topology_errors)

    if material_domain_ids is None:
        if len(raw_solids) == 1:
            domain_ids: tuple[str | None, ...] = ("solid-0",)
        else:
            domain_ids = tuple(None for _ in raw_solids)
            if len(raw_solids) > 1:
                reasons.append("MULTI_SOLID_DOMAIN_IDS_REQUIRED")
    else:
        domain_ids = tuple(material_domain_ids)
        if len(domain_ids) != len(raw_solids):
            reasons.append("MATERIAL_DOMAIN_SOLID_COUNT_MISMATCH")
            domain_ids = tuple(None for _ in raw_solids)
        if any(not isinstance(value, str) or not value.strip() for value in domain_ids):
            reasons.append("MATERIAL_DOMAIN_ID_INVALID")

    solid_records = tuple(
        MaterialSolid(
            solid_index=index,
            material_domain_id=domain_ids[index] if index < len(domain_ids) else None,
            shape=solid,
            faces=_solid_faces(solid),
        )
        for index, solid in enumerate(raw_solids)
    )
    clean_reasons = tuple(sorted(set(reasons)))
    return CanonicalStepSource(
        source_path=None,
        source_origin="IN_MEMORY_TEST_FIXTURE",
        source_geometry_sha256=None,
        expected_sha256=None,
        source_length_unit=source_length_unit,
        frame=frame if isinstance(frame, str) else None,
        source_accuracy_mm=accuracy,
        source_accuracy_basis=source_accuracy_basis,
        shape=shape,
        solids=solid_records,
        status="READY" if not clean_reasons else "UNKNOWN",
        unknown_reasons=clean_reasons,
        backend_versions=dict(backend_versions or {}),
    )


def load_canonical_step(
    canonical_step_path: str | Path,
    *,
    expected_sha256: str,
    source_length_unit: str,
    frame: str,
    source_accuracy_mm: float | None,
    source_accuracy_basis: str | None,
    material_domain_ids: Sequence[str] | None = None,
) -> CanonicalStepSource:
    """Load canonical STEP geometry in mm after authenticating its current bytes.

    STEP input must be attested as mm by its source manifest. CadQuery's STEP
    importer is explicitly asked to normalize the imported coordinates to mm;
    the result records the conversion to SI as ``1 mm = 1e-3 m``. A mismatch
    between the caller's expected SHA-256 and the file's actual bytes raises.
    Geometry or precision defects return an UNKNOWN source for pair-level
    fail-closed classification.
    """
    if not _is_sha256(expected_sha256):
        raise ValueError("expected_sha256 must be 64 lowercase hexadecimal characters")
    if source_length_unit != "mm":
        raise ValueError("canonical STEP source_length_unit must be explicitly 'mm'")
    if not isinstance(frame, str) or not frame.strip():
        raise ValueError("frame must be an explicit non-empty reference-frame label")
    path = Path(canonical_step_path).expanduser().resolve(strict=True)
    actual_sha256 = _sha256_file(path)
    if actual_sha256 != expected_sha256:
        raise ValueError(
            "canonical STEP SHA-256 mismatch: "
            f"expected {expected_sha256}, actual {actual_sha256}"
        )

    versions: dict[str, str] = {}
    try:
        import cadquery as cq
        import OCP

        versions = {
            "cadquery": str(getattr(cq, "__version__", "unknown")),
            "ocp": str(getattr(OCP, "__version__", "unknown")),
        }
        imported = cq.importers.importStep(str(path), unit="MM")
        roots = list(imported.objects)
        if not roots:
            return _unknown_source(
                source_path=str(path), source_origin="CANONICAL_STEP_FILE",
                source_geometry_sha256=actual_sha256, expected_sha256=expected_sha256,
                source_length_unit=source_length_unit, frame=frame,
                source_accuracy_mm=(float(source_accuracy_mm) if _finite_real(source_accuracy_mm) else None),
                source_accuracy_basis=source_accuracy_basis, shape=None,
                reasons=("STEP_IMPORT_RETURNED_NO_ROOT_SHAPES",), backend_versions=versions,
            )
        if len(roots) == 1:
            imported_shape = roots[0].wrapped
        else:
            imported_shape = cq.Compound.makeCompound(roots).wrapped
    except Exception as exc:
        return _unknown_source(
            source_path=str(path), source_origin="CANONICAL_STEP_FILE",
            source_geometry_sha256=actual_sha256, expected_sha256=expected_sha256,
            source_length_unit=source_length_unit, frame=frame,
            source_accuracy_mm=(float(source_accuracy_mm) if _finite_real(source_accuracy_mm) else None),
            source_accuracy_basis=source_accuracy_basis, shape=None,
            reasons=(f"STEP_IMPORT_FAILED_{type(exc).__name__}",), backend_versions=versions,
        )

    inspected = _source_from_shape(
        imported_shape,
        source_length_unit=source_length_unit,
        frame=frame,
        source_accuracy_mm=source_accuracy_mm,
        source_accuracy_basis=source_accuracy_basis,
        material_domain_ids=material_domain_ids,
        backend_versions=versions,
    )
    # Rebind the in-memory topology inspection to the actual, freshly hashed file.
    return CanonicalStepSource(
        source_path=str(path),
        source_origin="CANONICAL_STEP_FILE",
        source_geometry_sha256=actual_sha256,
        expected_sha256=expected_sha256,
        source_length_unit=source_length_unit,
        frame=frame,
        source_accuracy_mm=inspected.source_accuracy_mm,
        source_accuracy_basis=source_accuracy_basis,
        shape=imported_shape,
        solids=inspected.solids,
        status=inspected.status,
        unknown_reasons=inspected.unknown_reasons,
        backend_versions=versions,
    )


def _point_classifier_state(solid: MaterialSolid, point_mm: Sequence[float], tolerance_mm: float) -> str:
    from OCP.BRepClass3d import BRepClass3d_SolidClassifier
    from OCP.TopAbs import TopAbs_IN, TopAbs_ON, TopAbs_OUT
    from OCP.gp import gp_Pnt

    try:
        classifier = BRepClass3d_SolidClassifier()
        classifier.Load(solid.shape)
        classifier.Perform(gp_Pnt(*map(float, point_mm)), float(tolerance_mm))
        state = classifier.State()
    except Exception:
        return "UNKNOWN"
    if state == TopAbs_IN:
        return "IN"
    if state == TopAbs_OUT:
        return "OUT"
    if state == TopAbs_ON:
        return "ON"
    return "UNKNOWN"


def _collect_raw_hits(
    source: CanonicalStepSource,
    start_mm: Sequence[float],
    unit_direction: Sequence[float],
    length_mm: float,
    distance_tolerance_mm: float,
) -> tuple[list[dict[str, Any]], list[str]]:
    from OCP.IntCurvesFace import IntCurvesFace_ShapeIntersector
    from OCP.gp import gp_Dir, gp_Lin, gp_Pnt

    hits: list[dict[str, Any]] = []
    failures: list[str] = []
    shapes: list[tuple[int | None, str | None, Any, tuple[Any, ...]]] = [
        (solid.solid_index, solid.material_domain_id, solid.shape, solid.faces)
        for solid in source.solids
    ]
    if not shapes and source.shape is not None:
        shapes.append((None, None, source.shape, _solid_faces(source.shape)))

    line = gp_Lin(gp_Pnt(*map(float, start_mm)), gp_Dir(*map(float, unit_direction)))
    for solid_index, domain_id, shape, faces in shapes:
        try:
            intersector = IntCurvesFace_ShapeIntersector()
            intersector.Load(shape, float(distance_tolerance_mm))
            intersector.Perform(line, 0.0, float(length_mm))
            if not intersector.IsDone():
                failures.append("INTERSECTOR_NOT_DONE")
                continue
            for hit_index in range(1, int(intersector.NbPnt()) + 1):
                parameter_mm = float(intersector.WParameter(hit_index))
                if not math.isfinite(parameter_mm):
                    failures.append("NONFINITE_INTERSECTION_PARAMETER")
                    continue
                hit_face = intersector.Face(hit_index)
                local_face_index = next(
                    (index for index, face in enumerate(faces) if hit_face.IsSame(face)),
                    None,
                )
                point = [
                    float(start_mm[axis]) + parameter_mm * float(unit_direction[axis])
                    for axis in range(3)
                ]
                hits.append(
                    {
                        "solid_index": solid_index,
                        "material_domain_id": domain_id,
                        "raw_hit_index": hit_index,
                        "face_index_zero_based": local_face_index,
                        "parameter_mm_from_canonical_pair_start": parameter_mm,
                        "point_mm": point,
                        "point_m": [value * _MM_TO_M for value in point],
                        "face_identity_status": "MAPPED" if local_face_index is not None else "UNKNOWN",
                        "cluster_index": None,
                    }
                )
                if local_face_index is None:
                    failures.append("INTERSECTION_FACE_ID_UNRESOLVED")
        except Exception:
            failures.append("INTERSECTOR_FAILED")
    hits.sort(
        key=lambda item: (
            item["parameter_mm_from_canonical_pair_start"],
            -1 if item["solid_index"] is None else item["solid_index"],
            -1 if item["face_index_zero_based"] is None else item["face_index_zero_based"],
            item["raw_hit_index"],
        )
    )
    return hits, sorted(set(failures))


def _cluster_hits(
    hits: list[dict[str, Any]], distance_tolerance_mm: float
) -> tuple[list[dict[str, Any]], list[str]]:
    groups: list[list[dict[str, Any]]] = []
    for hit in hits:
        parameter = float(hit["parameter_mm_from_canonical_pair_start"])
        if not groups:
            groups.append([hit])
            continue
        previous = float(groups[-1][-1]["parameter_mm_from_canonical_pair_start"])
        if parameter - previous <= distance_tolerance_mm:
            groups[-1].append(hit)
        else:
            groups.append([hit])

    events: list[dict[str, Any]] = []
    reasons: list[str] = []
    for cluster_index, group in enumerate(groups):
        parameters = [float(hit["parameter_mm_from_canonical_pair_start"]) for hit in group]
        low = min(parameters)
        high = max(parameters)
        exact_parameters = {value for value in parameters}
        solid_faces = {
            (hit["solid_index"], hit["face_index_zero_based"])
            for hit in group
        }
        if high - low > 0.0 and high - low <= distance_tolerance_mm:
            reasons.append("UNRESOLVED_NEARBY_BOUNDARY_HITS")
        if len(solid_faces) > 1:
            reasons.append("AMBIGUOUS_COLOCATED_BOUNDARY_HITS")
        for hit in group:
            hit["cluster_index"] = cluster_index
        events.append(
            {
                "cluster_index": cluster_index,
                "parameter_mm": math.fsum(parameters) / len(parameters),
                "parameter_span_mm": high - low,
                "raw_hit_count": len(group),
                "distinct_parameter_count": len(exact_parameters),
                "solid_face_refs": [
                    {"solid_index": solid_index, "face_index_zero_based": face_index}
                    for solid_index, face_index in sorted(
                        solid_faces,
                        key=lambda value: (
                            -1 if value[0] is None else value[0],
                            -1 if value[1] is None else value[1],
                        ),
                    )
                ],
            }
        )
    return events, sorted(set(reasons))


def _interval_occupancy(
    source: CanonicalStepSource,
    midpoint_mm: Sequence[float],
    tolerance_mm: float,
) -> tuple[str, str | None, list[dict[str, Any]]]:
    states = [
        {
            "solid_index": solid.solid_index,
            "material_domain_id": solid.material_domain_id,
            "state": _point_classifier_state(solid, midpoint_mm, tolerance_mm),
        }
        for solid in source.solids
    ]
    if source.status != "READY":
        return "UNKNOWN", None, states
    inside = [state for state in states if state["state"] == "IN"]
    on = [state for state in states if state["state"] in {"ON", "UNKNOWN"}]
    if len(inside) > 1:
        return "UNKNOWN", None, states
    if on:
        return "UNKNOWN", None, states
    if len(inside) == 1:
        return "IN_MATERIAL", inside[0]["material_domain_id"], states
    return "OUT_OR_VOID", None, states


def _endpoint_result(
    source: CanonicalStepSource,
    point_mm: Sequence[float],
    declared_domain_id: str | None,
    tolerance_mm: float,
) -> dict[str, Any]:
    states = [
        {
            "solid_index": solid.solid_index,
            "material_domain_id": solid.material_domain_id,
            "state": _point_classifier_state(solid, point_mm, tolerance_mm),
        }
        for solid in source.solids
    ]
    if source.status != "READY":
        return {"state": "UNKNOWN", "reason": "SOURCE_NOT_READY", "solid_states": states}
    if not isinstance(declared_domain_id, str) or not declared_domain_id.strip():
        return {"state": "UNKNOWN", "reason": "ENDPOINT_DOMAIN_UNASSIGNED", "solid_states": states}
    owned = [state for state in states if state["material_domain_id"] == declared_domain_id]
    if not owned:
        return {"state": "UNKNOWN", "reason": "ENDPOINT_DOMAIN_NOT_IN_SOURCE", "solid_states": states}
    if any(state["state"] in {"ON", "UNKNOWN"} for state in states):
        return {"state": "UNKNOWN", "reason": "ENDPOINT_ON_OR_UNCERTAIN_BOUNDARY", "solid_states": states}
    owned_inside = [state for state in owned if state["state"] == "IN"]
    all_inside = [state for state in states if state["state"] == "IN"]
    if not owned_inside:
        if all_inside:
            return {"state": "OUT_OR_VOID", "reason": "ENDPOINT_DOMAIN_MISMATCH", "solid_states": states}
        return {"state": "OUT_OR_VOID", "reason": "ENDPOINT_OUTSIDE_MATERIAL", "solid_states": states}
    if len(owned_inside) != 1 or len(all_inside) != 1:
        return {"state": "UNKNOWN", "reason": "ENDPOINT_OVERLAPPING_SOLIDS", "solid_states": states}
    return {"state": "IN_MATERIAL", "reason": "ENDPOINT_IN_DECLARED_DOMAIN", "solid_states": states}


def _pair_intervals(
    source: CanonicalStepSource,
    start_mm: Sequence[float],
    unit_direction: Sequence[float],
    length_mm: float,
    events: Sequence[Mapping[str, Any]],
    uncertainty_band_mm: float | None,
) -> tuple[list[dict[str, Any]], list[str]]:
    cuts = [0.0]
    cuts.extend(
        float(event["parameter_mm"])
        for event in events
        if 0.0 < float(event["parameter_mm"]) < length_mm
    )
    cuts.append(float(length_mm))
    cuts.sort()
    intervals: list[dict[str, Any]] = []
    reasons: list[str] = []
    for start_t, end_t in zip(cuts, cuts[1:]):
        width = end_t - start_t
        if width <= 0.0:
            reasons.append("ZERO_WIDTH_INTERVAL_FROM_INTERSECTIONS")
            continue
        midpoint_t = (start_t + end_t) * 0.5
        midpoint = [
            float(start_mm[axis]) + midpoint_t * float(unit_direction[axis])
            for axis in range(3)
        ]
        if uncertainty_band_mm is None:
            occupancy, domain_id, solid_states = "UNKNOWN", None, []
            reason = "SOURCE_PRECISION_UNKNOWN"
        elif width <= 2.0 * uncertainty_band_mm:
            # Keep the exact interval and its diagnostic classifier states, but
            # do not trust a layer at or below the combined source/error band.
            _, domain_id, solid_states = _interval_occupancy(
                source, midpoint, uncertainty_band_mm
            )
            occupancy = "UNKNOWN"
            reason = "THIN_INTERVAL_WITHIN_ERROR_BAND"
            reasons.append(reason)
        else:
            occupancy, domain_id, solid_states = _interval_occupancy(
                source, midpoint, uncertainty_band_mm
            )
            reason = {
                "IN_MATERIAL": "OPEN_INTERVAL_INSIDE_SOLID",
                "OUT_OR_VOID": "OPEN_INTERVAL_OUTSIDE_ALL_SOLIDS",
                "UNKNOWN": "OPEN_INTERVAL_OCCUPANCY_UNKNOWN",
            }[occupancy]
            if occupancy == "UNKNOWN":
                reasons.append(reason)
        intervals.append(
            {
                "t_start_mm": start_t,
                "t_end_mm": end_t,
                "length_mm": width,
                "length_m": width * _MM_TO_M,
                "midpoint_t_mm": midpoint_t,
                "midpoint_mm": midpoint,
                "midpoint_m": [value * _MM_TO_M for value in midpoint],
                "occupancy": occupancy,
                "material_domain_id": domain_id,
                "reason": reason,
                "solid_classifier_states": solid_states,
            }
        )
    return intervals, sorted(set(reasons))


def classify_candidate_pairs(
    source: CanonicalStepSource,
    *,
    reference_point_ids: Sequence[int],
    reference_positions_mm: Sequence[Sequence[Real]],
    point_material_domain_ids: Sequence[str | None],
    candidate_pairs: Iterable[tuple[int, int]],
    reference_frame: str,
    distance_tolerance_mm: float,
) -> dict[str, Any]:
    """Classify canonical undirected horizon pairs in the fixed reference frame.

    Coordinates and distances are supplied in mm. ``distance_tolerance_mm``
    is explicit and is used for OCC intersections; the occupancy error band is
    ``source_accuracy_mm + distance_tolerance_mm``. An open interval whose
    length is at most twice that band is UNKNOWN. The output is a static
    topology ledger only and does not mutate or initialize a native system.
    """
    if not _finite_real(distance_tolerance_mm) or float(distance_tolerance_mm) <= 0.0:
        raise ValueError("distance_tolerance_mm must be finite and positive")
    if reference_frame != source.frame:
        raise ValueError("reference_frame must exactly match the canonical STEP frame")
    raw_ids = list(reference_point_ids)
    positions = list(reference_positions_mm)
    domains = list(point_material_domain_ids)
    if len(raw_ids) != len(positions) or len(raw_ids) != len(domains):
        raise ValueError("point IDs, reference positions, and material domains must have equal lengths")
    if any(not isinstance(value, Integral) or isinstance(value, bool) or value < 0 for value in raw_ids):
        raise ValueError("reference_point_ids must be unique non-negative integers")
    ids = [int(value) for value in raw_ids]
    if len(set(ids)) != len(ids):
        raise ValueError("reference_point_ids must be unique non-negative integers")
    by_id = {point_id: row for row, point_id in enumerate(ids)}
    xyz_by_id = {
        point_id: _xyz_mm(positions[row], label=f"reference position for point {point_id}")
        for row, point_id in enumerate(ids)
    }

    canonical_pairs: set[tuple[int, int]] = set()
    directed_input_count = 0
    for pair in candidate_pairs:
        if isinstance(pair, (str, bytes, Mapping)):
            raise ValueError("each candidate pair must contain exactly two point IDs")
        try:
            if len(pair) != 2:
                raise ValueError("each candidate pair must contain exactly two point IDs")
            point_i, point_j = pair[0], pair[1]
        except (TypeError, IndexError):
            raise ValueError("each candidate pair must contain exactly two point IDs") from None
        if (
            not isinstance(point_i, Integral) or isinstance(point_i, bool)
            or not isinstance(point_j, Integral) or isinstance(point_j, bool)
            or point_i == point_j
        ):
            raise ValueError("candidate pairs require distinct integer point IDs")
        point_i, point_j = int(point_i), int(point_j)
        if point_i not in by_id or point_j not in by_id:
            raise ValueError("candidate pair references an unknown point ID")
        canonical_pairs.add((min(point_i, point_j), max(point_i, point_j)))
        directed_input_count += 1

    accuracy_known = _finite_real(source.source_accuracy_mm) and float(source.source_accuracy_mm) >= 0.0
    uncertainty_band_mm = (
        float(source.source_accuracy_mm) + float(distance_tolerance_mm)
        if accuracy_known
        else None
    )
    pair_records: list[dict[str, Any]] = []
    for point_i, point_j in sorted(canonical_pairs):
        row_i, row_j = by_id[point_i], by_id[point_j]
        start_mm = xyz_by_id[point_i]
        end_mm = xyz_by_id[point_j]
        vector = tuple(end_mm[axis] - start_mm[axis] for axis in range(3))
        length_mm = math.sqrt(math.fsum(component * component for component in vector))
        pair_reasons: list[str] = list(source.unknown_reasons)
        if not math.isfinite(length_mm) or length_mm <= 0.0:
            pair_records.append(
                _empty_pair_record(point_i, point_j, start_mm, end_mm, "UNKNOWN", "ZERO_LENGTH_SEGMENT", length_mm)
            )
            continue
        unit_direction = tuple(component / length_mm for component in vector)
        raw_hits, hit_failures = _collect_raw_hits(
            source, start_mm, unit_direction, length_mm, float(distance_tolerance_mm)
        )
        hit_events, hit_reasons = _cluster_hits(raw_hits, float(distance_tolerance_mm))
        pair_reasons.extend(hit_failures)
        pair_reasons.extend(hit_reasons)
        pair_domains_differ = domains[row_i] != domains[row_j]
        endpoint_i = _endpoint_result(source, start_mm, domains[row_i], uncertainty_band_mm or 0.0)
        endpoint_j = _endpoint_result(source, end_mm, domains[row_j], uncertainty_band_mm or 0.0)
        intervals, interval_reasons = _pair_intervals(
            source,
            start_mm,
            unit_direction,
            length_mm,
            hit_events,
            uncertainty_band_mm,
        )
        pair_reasons.extend(interval_reasons)

        tangent_events: list[int] = []
        for event in hit_events:
            event_t = float(event["parameter_mm"])
            before = next(
                (item for item in intervals if math.isclose(float(item["t_end_mm"]), event_t, rel_tol=0.0, abs_tol=float(distance_tolerance_mm))),
                None,
            )
            after = next(
                (item for item in intervals if math.isclose(float(item["t_start_mm"]), event_t, rel_tol=0.0, abs_tol=float(distance_tolerance_mm))),
                None,
            )
            if (
                before is not None and after is not None
                and before["occupancy"] == after["occupancy"]
                and before["material_domain_id"] == after["material_domain_id"]
                and before["occupancy"] != "UNKNOWN"
            ):
                tangent_events.append(int(event["cluster_index"]))
        if tangent_events:
            pair_reasons.append("TANGENT_OR_NON_CROSSING_HIT")

        state = "UNKNOWN"
        primary_reason = "GEOMETRY_UNKNOWN"
        if source.status != "READY":
            state, primary_reason = "UNKNOWN", "SOURCE_NOT_CLASSIFIABLE"
        elif hit_failures:
            state, primary_reason = "UNKNOWN", "INTERSECTION_RESULT_INCOMPLETE"
        elif any(item["state"] == "UNKNOWN" for item in (endpoint_i, endpoint_j)):
            state, primary_reason = "UNKNOWN", "ENDPOINT_CLASSIFICATION_UNKNOWN"
        elif hit_reasons:
            state, primary_reason = "UNKNOWN", hit_reasons[0]
        elif tangent_events:
            state, primary_reason = "UNKNOWN", "TANGENT_OR_NON_CROSSING_HIT"
        elif any(item["occupancy"] == "UNKNOWN" for item in intervals):
            state, primary_reason = "UNKNOWN", "OPEN_INTERVAL_OCCUPANCY_UNKNOWN"
        elif pair_domains_differ:
            state, primary_reason = "OUT_OR_VOID", "DIFFERENT_MATERIAL_DOMAIN"
        elif any(item["state"] == "OUT_OR_VOID" for item in (endpoint_i, endpoint_j)):
            state, primary_reason = "OUT_OR_VOID", "ENDPOINT_OUTSIDE_OR_WRONG_DOMAIN"
        elif any(item["occupancy"] == "OUT_OR_VOID" for item in intervals):
            state, primary_reason = "OUT_OR_VOID", "AIR_PATH"
        elif any(item["material_domain_id"] != domains[row_i] for item in intervals):
            state, primary_reason = "OUT_OR_VOID", "DISCONNECTED_MATERIAL_DOMAIN"
        else:
            state, primary_reason = "IN_MATERIAL", "WHOLE_OPEN_SEGMENT_IN_DECLARED_MATERIAL_DOMAIN"

        directions = [
            {"source_point_id": point_i, "target_point_id": point_j, "state": state},
            {"source_point_id": point_j, "target_point_id": point_i, "state": state},
        ]
        pair_records.append(
            {
                "pair": [point_i, point_j],
                "classification": state,
                "reason": primary_reason,
                "reason_codes": sorted(set(pair_reasons + [primary_reason])),
                "reciprocal": {"directed_records": directions, "same_classification": True},
                "material_domain_ids": [domains[row_i], domains[row_j]],
                "endpoints": [
                    {
                        "point_id": point_i,
                        "position_mm": list(start_mm),
                        "position_m": [value * _MM_TO_M for value in start_mm],
                        **endpoint_i,
                    },
                    {
                        "point_id": point_j,
                        "position_mm": list(end_mm),
                        "position_m": [value * _MM_TO_M for value in end_mm],
                        **endpoint_j,
                    },
                ],
                "segment_length_mm": length_mm,
                "segment_length_m": length_mm * _MM_TO_M,
                "raw_hits": raw_hits,
                "hit_events": hit_events,
                "open_intervals": intervals,
                "tangent_or_non_crossing_event_clusters": tangent_events,
                "admissible_static_candidate": state == "IN_MATERIAL",
            }
        )

    counts = {state: 0 for state in ("IN_MATERIAL", "OUT_OR_VOID", "UNKNOWN")}
    for record in pair_records:
        counts[record["classification"]] += 1
    return {
        "schema": BOND_VISIBILITY_SCHEMA,
        "status": "STATIC_CLASSIFICATION_COMPLETE",
        "algorithm": BOND_VISIBILITY_ALGORITHM,
        "scope": "static reference-coordinate candidate classification only; native graph is not connected",
        "source": {
            "source_origin": source.source_origin,
            "source_path": source.source_path,
            "source_geometry_sha256": source.source_geometry_sha256,
            "expected_sha256": source.expected_sha256,
            "sha256_status": (
                "MATCHED_ACTUAL_FILE_BYTES"
                if source.source_origin == "CANONICAL_STEP_FILE"
                and source.source_geometry_sha256 == source.expected_sha256
                else "UNBOUND_IN_MEMORY_FIXTURE"
                if source.source_origin == "IN_MEMORY_TEST_FIXTURE"
                else "UNKNOWN"
            ),
            "status": source.status,
            "unknown_reasons": list(source.unknown_reasons),
            "source_length_unit": source.source_length_unit,
            "si_conversion": {"m_per_mm": _MM_TO_M, "volume_m3_per_mm3": 1.0e-9},
            "frame": source.frame,
            "source_accuracy_mm": source.source_accuracy_mm,
            "source_accuracy_basis": source.source_accuracy_basis,
            "backend_versions": dict(source.backend_versions),
            "solid_count": len(source.solids),
            "material_domain_ids": [solid.material_domain_id for solid in source.solids],
        },
        "tolerances": {
            "distance_tolerance_mm": float(distance_tolerance_mm),
            "source_accuracy_mm": source.source_accuracy_mm,
            "combined_occupancy_uncertainty_band_mm": uncertainty_band_mm,
            "thin_open_interval_unknown_at_or_below_mm": (
                2.0 * uncertainty_band_mm if uncertainty_band_mm is not None else None
            ),
            "basis": "explicit OCC distance tolerance plus caller-attested source accuracy; no tolerance is inferred",
        },
        "reference_coordinates": {
            "coordinate_unit": "mm",
            "frame": reference_frame,
            "configuration": "fixed undeformed reference coordinates",
            "point_count": len(ids),
        },
        "initial_topology_ledger": {
            "ledger_kind": "static material-visibility mask over Euclidean horizon candidates",
            "candidate_directed_input_count": directed_input_count,
            "candidate_undirected_pair_count": len(pair_records),
            "classification_counts": counts,
            "initial_material_bond_candidate_count": counts["IN_MATERIAL"],
            "geometry_excluded_pair_count": counts["OUT_OR_VOID"],
            "unknown_pair_count": counts["UNKNOWN"],
            "unknown_retained_separately": True,
            "native_graph_mutated": False,
            "native_neighbor_counts_recomputed": False,
            "native_physics_initialized": False,
            "whole_asset_physics_status": "NOT_ASSESSED",
        },
        "pairs": pair_records,
    }


def _empty_pair_record(
    point_i: int,
    point_j: int,
    start_mm: Sequence[float],
    end_mm: Sequence[float],
    state: str,
    reason: str,
    length_mm: float,
) -> dict[str, Any]:
    directions = [
        {"source_point_id": point_i, "target_point_id": point_j, "state": state},
        {"source_point_id": point_j, "target_point_id": point_i, "state": state},
    ]
    return {
        "pair": [point_i, point_j],
        "classification": state,
        "reason": reason,
        "reason_codes": [reason],
        "reciprocal": {"directed_records": directions, "same_classification": True},
        "endpoints": [
            {"point_id": point_i, "position_mm": list(start_mm)},
            {"point_id": point_j, "position_mm": list(end_mm)},
        ],
        "segment_length_mm": length_mm,
        "segment_length_m": length_mm * _MM_TO_M if math.isfinite(length_mm) else None,
        "raw_hits": [],
        "hit_events": [],
        "open_intervals": [],
        "tangent_or_non_crossing_event_clusters": [],
        "admissible_static_candidate": False,
    }


__all__ = [
    "BOND_VISIBILITY_ALGORITHM",
    "BOND_VISIBILITY_SCHEMA",
    "CanonicalStepSource",
    "MaterialSolid",
    "classify_candidate_pairs",
    "load_canonical_step",
]
