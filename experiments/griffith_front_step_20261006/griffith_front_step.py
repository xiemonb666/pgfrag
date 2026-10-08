"""One-step, fail-closed Mode III front-growth controller.

This is an isolated experimental state machine.  It proposes geometry from a
declared +x finite straight front and a supplied J/Gc comparison.  A proposal
never changes the active surface.  Confirmation requires a newly solved
``ExplicitSurfaceState`` bound to the proposal plus signed energy-ledger
evidence.  Receipts are data supplied by an upstream solver/reviewer; this
module does not solve mechanics or independently verify those receipts.

All committed versions remain ``NOT_COUPLED`` under the existing
``explicit_surface_state`` contract and are labelled
``EXPERIMENTAL_UNCALIBRATED``.  A state-machine commit means only that the
required receipts were supplied and bound consistently; it is not a physical
qualification or an energy-balance pass.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import hashlib
import json
import math
import re
import sys
from pathlib import Path
from typing import Any, Mapping

import numpy as np


_CONTRACT_DIR = Path(__file__).resolve().parents[1] / "explicit_surface_state_20261005"
if str(_CONTRACT_DIR) not in sys.path:
    sys.path.insert(0, str(_CONTRACT_DIR))
from explicit_surface_state import (  # noqa: E402
    CrackPatch,
    ExplicitSurfaceState,
    SurfaceSide,
    SurfaceStateError,
    validate_surface_state,
)


EXPERIMENTAL_QUALIFICATION = "EXPERIMENTAL_UNCALIBRATED"
ENERGY_UNIT = "J/m^2"
HASH_PATTERN = re.compile(r"^[0-9a-f]{64}$")
SENSITIVITY_REVIEWED = "REVIEWED"


class FrontStepError(ValueError):
    """Invalid controller state or a transition not allowed by this kernel."""


def _finite_number(value: object) -> float | None:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, float, np.integer, np.floating)):
        return None
    try:
        number = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError, OverflowError):
        return None
    return number if math.isfinite(number) else None


def _sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _valid_hash(value: object) -> bool:
    return isinstance(value, str) and bool(HASH_PATTERN.fullmatch(value))


def _stable_id(kind: str, *parts: str) -> str:
    digest = _sha256_text("\n".join((kind, *parts)))
    return f"{kind}:{digest[:24]}"


def _json_hash(payload: Mapping[str, Any]) -> str:
    def json_default(value: object) -> Any:
        if isinstance(value, np.generic):
            return value.item()
        raise TypeError(f"unsupported hash payload value {type(value).__name__}")

    body = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False, default=json_default)
    return _sha256_text(body)


@dataclass(frozen=True, slots=True)
class FrontFrame:
    """Declared Mode III front: a finite z-line, advancing only in +x."""

    front_x_m: float
    front_y_m: float
    z_min_m: float
    z_max_m: float
    direction: str = "+x"

    def to_payload(self) -> dict[str, Any]:
        return {
            "front_x_m": self.front_x_m,
            "front_y_m": self.front_y_m,
            "z_min_m": self.z_min_m,
            "z_max_m": self.z_max_m,
            "direction": self.direction,
        }


@dataclass(frozen=True, slots=True)
class DomainBounds:
    """Axis-aligned finite domain bounds in reference meters."""

    x_min_m: float
    x_max_m: float
    y_min_m: float
    y_max_m: float
    z_min_m: float
    z_max_m: float

    def to_payload(self) -> dict[str, float]:
        return {
            "x_min_m": self.x_min_m,
            "x_max_m": self.x_max_m,
            "y_min_m": self.y_min_m,
            "y_max_m": self.y_max_m,
            "z_min_m": self.z_min_m,
            "z_max_m": self.z_max_m,
        }


@dataclass(frozen=True, slots=True)
class DrivingEvidence:
    """One current J/Gc assessment with explicit provenance and review states.

    ``REVIEWED`` says a declared sensitivity review receipt exists.  It does
    not mean the values passed a threshold or are path independent.
    """

    G_J_per_m2: float | None
    Gc_J_per_m2: float | None
    G_unit: str | None
    Gc_unit: str | None
    j_status: str | None
    solver_result_sha256: str | None
    source_sha256: str | None
    Gc_source_sha256: str | None
    bound_surface_state_sha256: str | None
    j_result_id: str | None
    coupling_status: str | None
    contour_sensitivity_status: str | None
    contour_review_sha256: str | None
    path_sensitivity_status: str | None
    path_review_sha256: str | None
    Gc_qualification: str = EXPERIMENTAL_QUALIFICATION

    def to_payload(self) -> dict[str, Any]:
        return {
            "G_J_per_m2": self.G_J_per_m2,
            "Gc_J_per_m2": self.Gc_J_per_m2,
            "G_unit": self.G_unit,
            "Gc_unit": self.Gc_unit,
            "j_status": self.j_status,
            "solver_result_sha256": self.solver_result_sha256,
            "source_sha256": self.source_sha256,
            "Gc_source_sha256": self.Gc_source_sha256,
            "bound_surface_state_sha256": self.bound_surface_state_sha256,
            "j_result_id": self.j_result_id,
            "coupling_status": self.coupling_status,
            "contour_sensitivity_status": self.contour_sensitivity_status,
            "contour_review_sha256": self.contour_review_sha256,
            "path_sensitivity_status": self.path_sensitivity_status,
            "path_review_sha256": self.path_review_sha256,
            "Gc_qualification": self.Gc_qualification,
        }


@dataclass(frozen=True, slots=True)
class SurfaceVersion:
    version_id: str
    parent_version_id: str | None
    state_sha256: str
    state: ExplicitSurfaceState
    qualification: str = EXPERIMENTAL_QUALIFICATION
    version_status: str = "COMMITTED"

    def to_payload(self, *, include_state: bool = False) -> dict[str, Any]:
        payload = {
            "version_id": self.version_id,
            "parent_version_id": self.parent_version_id,
            "state_sha256": self.state_sha256,
            "qualification": self.qualification,
            "version_status": self.version_status,
            "mechanics_qualification": self.state.mechanics_qualification,
            "mechanically_active": self.state.mechanically_active,
            "surface_geometry_sha256": [p.geometry_sha256 for p in self.state.patches],
        }
        if include_state:
            payload["explicit_surface_state"] = self.state.to_payload()
        return payload


@dataclass(frozen=True, slots=True)
class StripGeometry:
    patch_id: str
    parent_patch_id: str
    surface_id: str
    creation_step: int
    vertices_m: tuple[tuple[float, float, float], ...]
    triangles: tuple[tuple[int, int, int], ...]
    geometry_sha256: str
    area_m2: float

    def to_payload(self) -> dict[str, Any]:
        return {
            "patch_id": self.patch_id,
            "parent_patch_id": self.parent_patch_id,
            "surface_id": self.surface_id,
            "creation_step": self.creation_step,
            "reference_vertices_m": [list(row) for row in self.vertices_m],
            "triangles": [list(row) for row in self.triangles],
            "geometry_sha256": self.geometry_sha256,
            "reference_area_m2": self.area_m2,
        }


@dataclass(frozen=True, slots=True)
class Candidate:
    candidate_id: str
    parent_version_id: str
    parent_state_sha256: str
    geometry_sha256: str
    strip: StripGeometry
    frame: FrontFrame
    bounds: DomainBounds
    requested_advance_m: float
    step_cap_m: float
    G_J_per_m2: float
    Gc_J_per_m2: float
    driving_evidence_sha256: str
    status: str = "PROPOSED_WAITING_REANALYSIS"
    qualification: str = EXPERIMENTAL_QUALIFICATION
    evidence_origin: str | None = None
    energy_ledger: tuple[tuple[str, float], ...] = ()

    def to_payload(self) -> dict[str, Any]:
        return {
            "candidate_id": self.candidate_id,
            "parent_version_id": self.parent_version_id,
            "parent_state_sha256": self.parent_state_sha256,
            "candidate_geometry_sha256": self.geometry_sha256,
            "strip": self.strip.to_payload(),
            "frame": self.frame.to_payload(),
            "domain_bounds": self.bounds.to_payload(),
            "requested_advance_m": self.requested_advance_m,
            "step_cap_m": self.step_cap_m,
            "G_J_per_m2": self.G_J_per_m2,
            "Gc_J_per_m2": self.Gc_J_per_m2,
            "driving_evidence_sha256": self.driving_evidence_sha256,
            "status": self.status,
            "qualification": self.qualification,
            "evidence_origin": self.evidence_origin,
            "energy_ledger_J": {key: value for key, value in self.energy_ledger},
            "commit_semantics": "state activation after bound receipts; no physical qualification or energy-balance pass",
        }


@dataclass(frozen=True, slots=True)
class Event:
    event_id: str
    sequence: int
    kind: str
    status: str
    active_surface_version_id: str
    candidate_id: str | None
    reason: str | None
    evidence_hashes: tuple[str, ...] = ()

    def to_payload(self) -> dict[str, Any]:
        return {
            "event_id": self.event_id,
            "sequence": self.sequence,
            "kind": self.kind,
            "status": self.status,
            "active_surface_version_id": self.active_surface_version_id,
            "candidate_id": self.candidate_id,
            "reason": self.reason,
            "evidence_hashes": list(self.evidence_hashes),
        }


@dataclass(frozen=True, slots=True)
class ControllerSnapshot:
    committed_versions: tuple[SurfaceVersion, ...]
    active_surface_version_id: str
    candidates: tuple[Candidate, ...] = ()
    pending_candidate_id: str | None = None
    event_ledger: tuple[Event, ...] = ()

    @property
    def active_surface(self) -> SurfaceVersion:
        for version in self.committed_versions:
            if version.version_id == self.active_surface_version_id:
                return version
        raise FrontStepError("active surface version is absent from immutable version history")

    def to_payload(self, *, include_surfaces: bool = False) -> dict[str, Any]:
        return {
            "schema_version": "pgfrag.griffith-front-step.experimental.v1",
            "qualification": EXPERIMENTAL_QUALIFICATION,
            "active_surface_version_id": self.active_surface_version_id,
            "pending_candidate_id": self.pending_candidate_id,
            "committed_versions": [v.to_payload(include_state=include_surfaces) for v in self.committed_versions],
            "candidates": [candidate.to_payload() for candidate in self.candidates],
            "event_ledger": [event.to_payload() for event in self.event_ledger],
            "controller_claim": "EXPERIMENTAL_STATE_MACHINE_ONLY",
        }


def initialize(initial_state: ExplicitSurfaceState) -> ControllerSnapshot:
    """Freeze the caller's starting surface as the first committed version."""

    validate_surface_state(initial_state)
    state_sha = _sha256_text(initial_state.to_json())
    version_id = _stable_id("surface-version", state_sha)
    version = SurfaceVersion(version_id, None, state_sha, initial_state)
    return ControllerSnapshot((version,), version_id)


def _append_event(
    snapshot: ControllerSnapshot,
    *,
    kind: str,
    status: str,
    candidate_id: str | None = None,
    reason: str | None = None,
    evidence_hashes: tuple[str, ...] = (),
    active_surface_version_id: str | None = None,
) -> ControllerSnapshot:
    active_id = active_surface_version_id or snapshot.active_surface_version_id
    sequence = len(snapshot.event_ledger)
    eid = _stable_id(
        "front-event",
        str(sequence),
        kind,
        status,
        active_id,
        candidate_id or "",
        reason or "",
        *evidence_hashes,
    )
    event = Event(eid, sequence, kind, status, active_id, candidate_id, reason, evidence_hashes)
    return replace(snapshot, event_ledger=snapshot.event_ledger + (event,))


def _latest_leaf_patch(state: ExplicitSurfaceState) -> CrackPatch | None:
    parents = {p.parent_patch_id for p in state.patches if p.parent_patch_id is not None}
    leaves = [p for p in state.patches if p.patch_id not in parents]
    if not leaves:
        return None
    # This one-step kernel accepts only one declared straight front patch.
    return leaves[-1] if len(leaves) == 1 else None


def _frame_matches_patch(patch: CrackPatch, frame: FrontFrame) -> bool:
    values = [frame.front_x_m, frame.front_y_m, frame.z_min_m, frame.z_max_m]
    if any(_finite_number(value) is None for value in values) or frame.z_max_m <= frame.z_min_m:
        return False
    if frame.direction != "+x":
        return False
    vertices = patch.reference_vertices_m
    scale = max(1.0, float(np.max(np.abs(vertices))), *(abs(v) for v in values))
    tol = 64.0 * np.finfo(np.float64).eps * scale
    # Restrict this candidate to one planar x-z patch and locate the declared
    # finite line in the contract's oriented boundary edges.
    if np.max(np.abs(vertices[:, 1] - frame.front_y_m)) > tol:
        return False
    matches = 0
    for left, right in patch.front_edges.tolist():
        edge = vertices[[left, right]]
        if (
            np.max(np.abs(edge[:, 0] - frame.front_x_m)) <= tol
            and np.max(np.abs(edge[:, 1] - frame.front_y_m)) <= tol
            and abs(float(np.min(edge[:, 2])) - frame.z_min_m) <= tol
            and abs(float(np.max(edge[:, 2])) - frame.z_max_m) <= tol
        ):
            matches += 1
    return matches == 1


def _valid_bounds(bounds: DomainBounds) -> bool:
    values = tuple(bounds.to_payload().values())
    if any(_finite_number(v) is None for v in values):
        return False
    return (
        bounds.x_min_m < bounds.x_max_m
        and bounds.y_min_m < bounds.y_max_m
        and bounds.z_min_m < bounds.z_max_m
    )


def _make_strip(patch: CrackPatch, frame: FrontFrame, advance_m: float, step: int) -> StripGeometry:
    x0 = float(frame.front_x_m)
    x1 = x0 + advance_m
    y = float(frame.front_y_m)
    z0 = float(frame.z_min_m)
    z1 = float(frame.z_max_m)
    vertices = np.asarray(
        [[x0, y, z0], [x0, y, z1], [x1, y, z1], [x1, y, z0]],
        dtype=np.dtype("<f8"),
    )
    triangles = np.asarray([[0, 1, 2], [0, 2, 3]], dtype=np.dtype("<i8"))
    geom_sha = hashlib.sha256(
        b"pgfrag.explicit-crack-reference-geometry.v1\n"
        + vertices.tobytes(order="C")
        + b"--triangles--\n"
        + triangles.tobytes(order="C")
    ).hexdigest()
    patch_id = _stable_id("patch", patch.patch_id, str(step), geom_sha)
    return StripGeometry(
        patch_id=patch_id,
        parent_patch_id=patch.patch_id,
        surface_id=patch.surface_id,
        creation_step=step,
        vertices_m=tuple(tuple(float(v) for v in row) for row in vertices.tolist()),
        triangles=tuple(tuple(int(v) for v in row) for row in triangles.tolist()),
        geometry_sha256=geom_sha,
        area_m2=advance_m * (z1 - z0),
    )


def _candidate_geometry_sha(parent_sha: str, strip: StripGeometry) -> str:
    return _json_hash({"parent_surface_state_sha256": parent_sha, "strip": strip.to_payload()})


def _pending(snapshot: ControllerSnapshot) -> Candidate | None:
    if snapshot.pending_candidate_id is None:
        return None
    for candidate in snapshot.candidates:
        if candidate.candidate_id == snapshot.pending_candidate_id:
            return candidate
    return None


def propose_step(
    snapshot: ControllerSnapshot,
    evidence: DrivingEvidence | None,
    *,
    frame: FrontFrame,
    requested_advance_m: float | None,
    step_cap_m: float | None,
    bounds: DomainBounds,
) -> ControllerSnapshot:
    """Arrest, return UNKNOWN/NEEDS_RECOMPUTE, reject, or propose one strip."""

    if _pending(snapshot) is not None:
        return _append_event(
            snapshot, kind="PROPOSE", status="BLOCKED_PENDING_CANDIDATE",
            candidate_id=snapshot.pending_candidate_id,
            reason="a candidate must be confirmed or rejected before another proposal",
        )
    active = snapshot.active_surface
    if evidence is None:
        return _append_event(snapshot, kind="ASSESS", status="UNKNOWN", reason="missing_driving_evidence")
    e_payload = evidence.to_payload()
    try:
        e_sha = _json_hash(e_payload)
    except (TypeError, ValueError, OverflowError):
        return _append_event(
            snapshot, kind="ASSESS", status="UNKNOWN",
            reason="invalid_or_nonfinite_driving_evidence_payload",
        )

    if evidence.coupling_status != "COUPLED_CURRENT":
        status = "NEEDS_RECOMPUTE" if evidence.coupling_status in ("STALE", "NOT_COUPLED") else "UNKNOWN"
        return _append_event(
            snapshot, kind="ASSESS", status=status,
            reason=f"driving_field_{evidence.coupling_status or 'status_missing'}",
            evidence_hashes=(e_sha,),
        )
    if evidence.bound_surface_state_sha256 != active.state_sha256:
        return _append_event(
            snapshot, kind="ASSESS", status="NEEDS_RECOMPUTE",
            reason="driving_field_bound_to_stale_or_different_surface",
            evidence_hashes=(e_sha,),
        )
    if not all(_valid_hash(value) for value in (
        evidence.solver_result_sha256, evidence.source_sha256, evidence.Gc_source_sha256,
        evidence.contour_review_sha256, evidence.path_review_sha256,
    )) or not isinstance(evidence.j_result_id, str) or not evidence.j_result_id.strip():
        return _append_event(
            snapshot, kind="ASSESS", status="UNKNOWN",
            reason="missing_or_invalid_solver_source_or_sensitivity_receipt",
            evidence_hashes=(e_sha,),
        )
    if evidence.j_status != "OK":
        return _append_event(
            snapshot, kind="ASSESS", status="UNKNOWN",
            reason="J_result_status_not_OK", evidence_hashes=(e_sha,),
        )
    if evidence.G_unit != ENERGY_UNIT or evidence.Gc_unit != ENERGY_UNIT:
        return _append_event(
            snapshot, kind="ASSESS", status="UNKNOWN", reason="missing_or_invalid_energy_units",
            evidence_hashes=(e_sha,),
        )
    if (
        evidence.contour_sensitivity_status != SENSITIVITY_REVIEWED
        or evidence.path_sensitivity_status != SENSITIVITY_REVIEWED
    ):
        return _append_event(
            snapshot, kind="ASSESS", status="UNKNOWN",
            reason="contour_or_path_sensitivity_not_reviewed",
            evidence_hashes=(e_sha,),
        )
    if evidence.Gc_qualification != EXPERIMENTAL_QUALIFICATION:
        return _append_event(
            snapshot, kind="ASSESS", status="UNKNOWN",
            reason="Gc_qualification_must_remain_experimental_uncalibrated",
            evidence_hashes=(e_sha,),
        )
    G = _finite_number(evidence.G_J_per_m2)
    Gc = _finite_number(evidence.Gc_J_per_m2)
    if G is None or G < 0:
        return _append_event(
            snapshot, kind="ASSESS", status="UNKNOWN", reason="missing_or_invalid_energy_release_G",
            evidence_hashes=(e_sha,),
        )
    if Gc is None or Gc <= 0:
        return _append_event(
            snapshot, kind="ASSESS", status="REJECTED_INVALID_GC", reason="Gc_must_be_positive_finite",
            evidence_hashes=(e_sha,),
        )
    if G < Gc:
        return _append_event(
            snapshot, kind="ASSESS", status="ARRESTED", reason="G_below_Gc",
            evidence_hashes=(e_sha,),
        )
    if G == Gc:
        return _append_event(
            snapshot, kind="ASSESS", status="ARRESTED_AT_EQUALITY",
            reason="strict_proposal_rule_requires_G_greater_than_Gc",
            evidence_hashes=(e_sha,),
        )
    if not _valid_bounds(bounds):
        return _append_event(
            snapshot, kind="GEOMETRY", status="REJECTED_INVALID_DOMAIN", reason="invalid_domain_bounds",
            evidence_hashes=(e_sha,),
        )
    if frame.direction != "+x":
        return _append_event(
            snapshot, kind="GEOMETRY", status="REJECTED_REVERSED_DIRECTION",
            reason="only_the_declared_plus_x_direction_is_supported", evidence_hashes=(e_sha,),
        )
    parent_patch = _latest_leaf_patch(active.state)
    if parent_patch is None or not _frame_matches_patch(parent_patch, frame):
        return _append_event(
            snapshot, kind="GEOMETRY", status="REJECTED_SEED_FRAME_MISMATCH",
            reason="declared_finite_z_front_does_not_match_one_active_leaf_boundary_edge",
            evidence_hashes=(e_sha,),
        )
    advance = _finite_number(requested_advance_m)
    cap = _finite_number(step_cap_m)
    if advance is None or cap is None or advance <= 0 or cap <= 0:
        return _append_event(
            snapshot, kind="GEOMETRY", status="UNKNOWN", reason="missing_or_invalid_step_or_step_cap",
            evidence_hashes=(e_sha,),
        )
    if advance > cap:
        return _append_event(
            snapshot, kind="GEOMETRY", status="REJECTED_STEP_CAP",
            reason="requested_advance_exceeds_allowed_geometric_step_cap",
            evidence_hashes=(e_sha,),
        )
    x1 = frame.front_x_m + advance
    if (
        frame.front_x_m < bounds.x_min_m or x1 > bounds.x_max_m
        or frame.front_y_m < bounds.y_min_m or frame.front_y_m > bounds.y_max_m
        or frame.z_min_m < bounds.z_min_m or frame.z_max_m > bounds.z_max_m
    ):
        return _append_event(
            snapshot, kind="GEOMETRY", status="REJECTED_OUT_OF_DOMAIN",
            reason="proposed_strip_exceeds_domain_bounds", evidence_hashes=(e_sha,),
        )

    step = max(p.creation_step for p in active.state.patches) + 1
    strip = _make_strip(parent_patch, frame, advance, step)
    geometry_sha = _candidate_geometry_sha(active.state_sha256, strip)
    candidate_id = _stable_id(
        "front-candidate", active.version_id, geometry_sha, e_sha,
        _json_hash({
            "requested_advance_m": advance,
            "step_cap_m": cap,
            "frame": frame.to_payload(),
            "domain_bounds": bounds.to_payload(),
        }),
    )
    candidate = Candidate(
        candidate_id=candidate_id,
        parent_version_id=active.version_id,
        parent_state_sha256=active.state_sha256,
        geometry_sha256=geometry_sha,
        strip=strip,
        frame=frame,
        bounds=bounds,
        requested_advance_m=advance,
        step_cap_m=cap,
        G_J_per_m2=G,
        Gc_J_per_m2=Gc,
        driving_evidence_sha256=e_sha,
    )
    updated = replace(
        snapshot,
        candidates=snapshot.candidates + (candidate,),
        pending_candidate_id=candidate_id,
    )
    return _append_event(
        updated, kind="PROPOSE", status="PROPOSED_WAITING_REANALYSIS",
        candidate_id=candidate_id, reason="G_greater_than_Gc_candidate_only",
        evidence_hashes=(e_sha, active.state_sha256, geometry_sha),
    )


@dataclass(frozen=True, slots=True)
class ReanalysisEvidence:
    candidate_id: str | None
    mechanical_status: str | None
    solver_result_sha256: str | None
    source_sha256: str | None
    bound_parent_surface_state_sha256: str | None
    bound_candidate_geometry_sha256: str | None
    reanalyzed_surface_state_sha256: str | None
    energy_ledger_sha256: str | None
    energy_review_status: str | None
    energy_review_sha256: str | None
    delta_potential_J: float | None
    external_work_J: float | None
    Gc_delta_area_J: float | None
    area_increment_m2: float | None
    energy_unit: str | None = "J"
    area_unit: str | None = "m^2"
    evidence_origin: str = "EXTERNAL_RECEIPT"

    def to_payload(self) -> dict[str, Any]:
        return {
            "candidate_id": self.candidate_id,
            "mechanical_status": self.mechanical_status,
            "solver_result_sha256": self.solver_result_sha256,
            "source_sha256": self.source_sha256,
            "bound_parent_surface_state_sha256": self.bound_parent_surface_state_sha256,
            "bound_candidate_geometry_sha256": self.bound_candidate_geometry_sha256,
            "reanalyzed_surface_state_sha256": self.reanalyzed_surface_state_sha256,
            "energy_ledger_sha256": self.energy_ledger_sha256,
            "energy_review_status": self.energy_review_status,
            "energy_review_sha256": self.energy_review_sha256,
            "delta_potential_J": self.delta_potential_J,
            "external_work_J": self.external_work_J,
            "Gc_delta_area_J": self.Gc_delta_area_J,
            "area_increment_m2": self.area_increment_m2,
            "energy_unit": self.energy_unit,
            "area_unit": self.area_unit,
            "evidence_origin": self.evidence_origin,
        }


def _candidate(snapshot: ControllerSnapshot, candidate_id: str) -> Candidate | None:
    return next((c for c in snapshot.candidates if c.candidate_id == candidate_id), None)


def _state_geometry_matches_candidate(
    state: ExplicitSurfaceState, parent: SurfaceVersion, candidate: Candidate
) -> bool:
    try:
        validate_surface_state(state)
    except (SurfaceStateError, TypeError, ValueError):
        return False
    prior_by_id = {p.patch_id: p for p in parent.state.patches}
    state_by_id = {p.patch_id: p for p in state.patches}
    if set(state_by_id) != set(prior_by_id) | {candidate.strip.patch_id}:
        return False
    if state.assignments != parent.state.assignments:
        return False
    if state.historical_proxy_relations != parent.state.historical_proxy_relations:
        return False
    if state.current_visibility != parent.state.current_visibility:
        return False
    for patch_id, old in prior_by_id.items():
        new = state_by_id[patch_id]
        if (
            new.surface_id != old.surface_id
            or new.creation_step != old.creation_step
            or new.parent_patch_id != old.parent_patch_id
            or new.geometry_sha256 != old.geometry_sha256
            or new.negative_side.body_id != old.negative_side.body_id
            or new.negative_side.fragment_id != old.negative_side.fragment_id
            or new.positive_side.body_id != old.positive_side.body_id
            or new.positive_side.fragment_id != old.positive_side.fragment_id
        ):
            return False
    new_patch = state_by_id[candidate.strip.patch_id]
    parent_patch = prior_by_id[candidate.strip.parent_patch_id]
    return (
        new_patch.surface_id == candidate.strip.surface_id
        and new_patch.patch_id == candidate.strip.patch_id
        and new_patch.parent_patch_id == candidate.strip.parent_patch_id
        and new_patch.creation_step == candidate.strip.creation_step
        and new_patch.geometry_sha256 == candidate.strip.geometry_sha256
        and new_patch.negative_side.body_id == parent_patch.negative_side.body_id
        and new_patch.negative_side.fragment_id == parent_patch.negative_side.fragment_id
        and new_patch.positive_side.body_id == parent_patch.positive_side.body_id
        and new_patch.positive_side.fragment_id == parent_patch.positive_side.fragment_id
        and math.isclose(new_patch.reference_area_m2, candidate.strip.area_m2, rel_tol=1e-12, abs_tol=1e-18)
    )


def _mark_candidate(snapshot: ControllerSnapshot, candidate: Candidate, **changes: Any) -> ControllerSnapshot:
    updated = replace(candidate, **changes)
    candidates = tuple(updated if item.candidate_id == candidate.candidate_id else item for item in snapshot.candidates)
    return replace(snapshot, candidates=candidates)


def confirm_candidate(
    snapshot: ControllerSnapshot,
    candidate_id: str,
    *,
    reanalyzed_surface: ExplicitSurfaceState | None,
    evidence: ReanalysisEvidence | None,
) -> ControllerSnapshot:
    """Commit only after candidate-bound solve, state, and energy receipts exist."""

    candidate = _candidate(snapshot, candidate_id)
    if candidate is None or snapshot.pending_candidate_id != candidate_id:
        return _append_event(
            snapshot, kind="CONFIRM", status="REJECTED_UNKNOWN_CANDIDATE",
            candidate_id=candidate_id, reason="candidate_is_not_pending",
        )
    if candidate.status not in ("PROPOSED_WAITING_REANALYSIS", "NEEDS_RECOMPUTE"):
        return _append_event(
            snapshot, kind="CONFIRM", status="REJECTED_INVALID_CANDIDATE_STATE",
            candidate_id=candidate_id, reason=f"candidate_status_{candidate.status}",
        )
    active = snapshot.active_surface
    missing: list[str] = []
    if active.version_id != candidate.parent_version_id or active.state_sha256 != candidate.parent_state_sha256:
        missing.append("candidate_parent_is_no_longer_active")
    if evidence is None:
        missing.append("missing_reanalysis_and_energy_receipt")
        e_payload: dict[str, Any] = {}
        e_sha = ""
    else:
        e_payload = evidence.to_payload()
        try:
            e_sha = _json_hash(e_payload)
        except (TypeError, ValueError, OverflowError):
            e_sha = ""
            missing.append("invalid_or_nonfinite_reanalysis_evidence_payload")
        if evidence.candidate_id != candidate_id:
            missing.append("candidate_id_mismatch")
        mechanical_binding = (evidence.mechanical_status, evidence.evidence_origin)
        if mechanical_binding not in (
            ("REANALYZED_CURRENT", "EXTERNAL_RECEIPT"),
            ("EXPLICIT_FIXTURE_ONLY", "EXPLICIT_FIXTURE_NOT_SOLVER"),
        ):
            missing.append("mechanical_reanalysis_not_current")
        if evidence.bound_parent_surface_state_sha256 != candidate.parent_state_sha256:
            missing.append("parent_surface_binding_mismatch")
        if evidence.bound_candidate_geometry_sha256 != candidate.geometry_sha256:
            missing.append("candidate_geometry_binding_mismatch")
        if not all(_valid_hash(value) for value in (
            evidence.solver_result_sha256, evidence.source_sha256,
            evidence.reanalyzed_surface_state_sha256, evidence.energy_ledger_sha256,
            evidence.energy_review_sha256,
        )):
            missing.append("missing_or_invalid_reanalysis_source_or_energy_receipt_hash")
        if reanalyzed_surface is None:
            missing.append("missing_reanalyzed_surface_state")
        else:
            try:
                actual_state_sha = _sha256_text(reanalyzed_surface.to_json())
            except (SurfaceStateError, TypeError, ValueError):
                actual_state_sha = ""
            if actual_state_sha != evidence.reanalyzed_surface_state_sha256:
                missing.append("reanalyzed_surface_state_hash_mismatch")
            elif not _state_geometry_matches_candidate(reanalyzed_surface, active, candidate):
                missing.append("reanalyzed_surface_geometry_does_not_match_candidate")
        if evidence.energy_review_status != "REVIEWED" or evidence.energy_unit != "J" or evidence.area_unit != "m^2":
            missing.append("energy_or_area_review_not_available_with_declared_units")
        numeric = {
            "delta_potential_J": evidence.delta_potential_J,
            "external_work_J": evidence.external_work_J,
            "Gc_delta_area_J": evidence.Gc_delta_area_J,
            "area_increment_m2": evidence.area_increment_m2,
        }
        if any(_finite_number(value) is None for value in numeric.values()):
            missing.append("missing_or_nonfinite_energy_terms")
        else:
            dpot = float(evidence.delta_potential_J)  # type: ignore[arg-type]
            work = float(evidence.external_work_J)  # type: ignore[arg-type]
            gc_area = float(evidence.Gc_delta_area_J)  # type: ignore[arg-type]
            area = float(evidence.area_increment_m2)  # type: ignore[arg-type]
            expected_area = candidate.strip.area_m2
            expected_gc_area = candidate.Gc_J_per_m2 * expected_area
            if not math.isclose(area, expected_area, rel_tol=1e-12, abs_tol=1e-18):
                missing.append("energy_area_does_not_match_candidate_geometry")
            if not math.isclose(gc_area, expected_gc_area, rel_tol=1e-12, abs_tol=1e-18):
                missing.append("Gc_delta_area_does_not_match_declared_Gc_and_candidate_area")
    if missing:
        updated = _mark_candidate(snapshot, candidate, status="NEEDS_RECOMPUTE")
        return _append_event(
            updated, kind="CONFIRM", status="NEEDS_RECOMPUTE",
            candidate_id=candidate_id, reason=";".join(missing),
            evidence_hashes=tuple(value for value in (e_sha, candidate.geometry_sha256) if value),
        )

    assert evidence is not None and reanalyzed_surface is not None
    dpot = float(evidence.delta_potential_J)  # type: ignore[arg-type]
    work = float(evidence.external_work_J)  # type: ignore[arg-type]
    gc_area = float(evidence.Gc_delta_area_J)  # type: ignore[arg-type]
    residual = work - dpot - gc_area
    state_sha = _sha256_text(reanalyzed_surface.to_json())
    new_version_id = _stable_id("surface-version", candidate.candidate_id, state_sha)
    new_version = SurfaceVersion(
        version_id=new_version_id,
        parent_version_id=active.version_id,
        state_sha256=state_sha,
        state=reanalyzed_surface,
        version_status=(
            "FIXTURE_ONLY_TRANSITION"
            if evidence.evidence_origin == "EXPLICIT_FIXTURE_NOT_SOLVER"
            else "COMMITTED_WITH_BOUND_REANALYSIS_AND_ENERGY_EVIDENCE"
        ),
    )
    ledger = (
        ("delta_potential_J", dpot),
        ("external_work_J", work),
        ("Gc_delta_area_J", gc_area),
        ("energy_closure_residual_J_external_minus_delta_potential_minus_Gc_delta_area", residual),
    )
    committed_candidate = replace(
        candidate,
        status=(
            "COMMITTED_FIXTURE_ONLY"
            if evidence.evidence_origin == "EXPLICIT_FIXTURE_NOT_SOLVER"
            else "COMMITTED"
        ),
        evidence_origin=evidence.evidence_origin,
        energy_ledger=ledger,
    )
    candidates = tuple(
        committed_candidate if item.candidate_id == candidate_id else item
        for item in snapshot.candidates
    )
    updated = replace(
        snapshot,
        committed_versions=snapshot.committed_versions + (new_version,),
        active_surface_version_id=new_version_id,
        candidates=candidates,
        pending_candidate_id=None,
    )
    return _append_event(
        updated,
        kind="COMMIT",
        status=(
            "COMMITTED_FIXTURE_ONLY"
            if evidence.evidence_origin == "EXPLICIT_FIXTURE_NOT_SOLVER"
            else "COMMITTED_WITH_BOUND_REANALYSIS_AND_ENERGY_EVIDENCE"
        ),
        candidate_id=candidate_id,
        reason=(
            "synthetic_fixture_transition_only"
            if evidence.evidence_origin == "EXPLICIT_FIXTURE_NOT_SOLVER"
            else "evidence_complete; residual_recorded_without_balance_claim"
        ),
        evidence_hashes=(e_sha, state_sha, candidate.geometry_sha256),
        active_surface_version_id=new_version_id,
    )


def reject_candidate(snapshot: ControllerSnapshot, candidate_id: str, reason: str) -> ControllerSnapshot:
    """Reject one pending candidate while retaining all committed versions."""

    candidate = _candidate(snapshot, candidate_id)
    if candidate is None or snapshot.pending_candidate_id != candidate_id:
        return _append_event(
            snapshot, kind="REJECT", status="REJECTED_UNKNOWN_CANDIDATE",
            candidate_id=candidate_id, reason="candidate_is_not_pending",
        )
    if not isinstance(reason, str) or not reason.strip():
        raise FrontStepError("rejection requires a non-empty reason")
    updated = _mark_candidate(snapshot, candidate, status="REJECTED_BY_REVIEW")
    updated = replace(updated, pending_candidate_id=None)
    return _append_event(
        updated, kind="REJECT", status="REJECTED_BY_REVIEW",
        candidate_id=candidate_id, reason=reason.strip(),
        evidence_hashes=(candidate.geometry_sha256,),
    )


def make_fixture_reanalyzed_state(
    parent: SurfaceVersion,
    candidate: Candidate,
    *,
    displacement_z_m: float = 0.0,
) -> ExplicitSurfaceState:
    """Create a clearly synthetic explicit-state fixture for transition tests."""

    old_patch = next(p for p in parent.state.patches if p.patch_id == candidate.strip.parent_patch_id)
    values = np.asarray(candidate.strip.vertices_m, dtype=np.float64)
    triangles = np.asarray(candidate.strip.triangles, dtype=np.int64)
    negative_disp = np.zeros_like(values)
    positive_disp = np.zeros_like(values)
    positive_disp[:, 2] = displacement_z_m
    strip_patch = CrackPatch(
        surface_id=candidate.strip.surface_id,
        patch_id=candidate.strip.patch_id,
        creation_step=candidate.strip.creation_step,
        parent_patch_id=candidate.strip.parent_patch_id,
        reference_vertices_m=values,
        triangles=triangles,
        negative_side=SurfaceSide(-1, old_patch.negative_side.body_id, negative_disp, old_patch.negative_side.fragment_id),
        positive_side=SurfaceSide(+1, old_patch.positive_side.body_id, positive_disp, old_patch.positive_side.fragment_id),
    )
    return ExplicitSurfaceState(
        patches=parent.state.patches + (strip_patch,),
        assignments=parent.state.assignments,
        historical_proxy_relations=parent.state.historical_proxy_relations,
        current_visibility=parent.state.current_visibility,
    )


def make_fixture_driving_evidence(
    surface: SurfaceVersion,
    *,
    G: float,
    Gc: float,
    coupling_status: str = "COUPLED_CURRENT",
    contour_status: str = SENSITIVITY_REVIEWED,
    path_status: str = SENSITIVITY_REVIEWED,
) -> DrivingEvidence:
    """Deterministic hash-labelled fixture inputs; these are never solver data."""

    fixture = lambda name: _sha256_text("explicit-fixture-only:" + name)
    return DrivingEvidence(
        G_J_per_m2=G,
        Gc_J_per_m2=Gc,
        G_unit=ENERGY_UNIT,
        Gc_unit=ENERGY_UNIT,
        j_status="OK",
        solver_result_sha256=fixture("solver-result"),
        source_sha256=fixture("source"),
        Gc_source_sha256=fixture("Gc-source"),
        bound_surface_state_sha256=surface.state_sha256,
        j_result_id="fixture-j-result-not-solver",
        coupling_status=coupling_status,
        contour_sensitivity_status=contour_status,
        contour_review_sha256=fixture("contour-review"),
        path_sensitivity_status=path_status,
        path_review_sha256=fixture("path-review"),
    )


def make_fixture_reanalysis_evidence(
    candidate: Candidate,
    state: ExplicitSurfaceState,
    parent_state_sha256: str,
    *,
    energy_complete: bool = True,
) -> ReanalysisEvidence:
    """Deterministic synthetic receipt values used only by bounded tests."""

    fixture = lambda name: _sha256_text("explicit-fixture-only:" + name)
    return ReanalysisEvidence(
        candidate_id=candidate.candidate_id,
        mechanical_status="EXPLICIT_FIXTURE_ONLY",
        solver_result_sha256=fixture("reanalysis-solver-result"),
        source_sha256=fixture("reanalysis-source"),
        bound_parent_surface_state_sha256=parent_state_sha256,
        bound_candidate_geometry_sha256=candidate.geometry_sha256,
        reanalyzed_surface_state_sha256=_sha256_text(state.to_json()),
        energy_ledger_sha256=fixture("energy-ledger") if energy_complete else None,
        energy_review_status="REVIEWED" if energy_complete else "NOT_CHECKED",
        energy_review_sha256=fixture("energy-review") if energy_complete else None,
        delta_potential_J=0.10 if energy_complete else None,
        external_work_J=0.13 if energy_complete else None,
        Gc_delta_area_J=candidate.Gc_J_per_m2 * candidate.strip.area_m2 if energy_complete else None,
        area_increment_m2=candidate.strip.area_m2 if energy_complete else None,
        evidence_origin="EXPLICIT_FIXTURE_NOT_SOLVER",
    )

