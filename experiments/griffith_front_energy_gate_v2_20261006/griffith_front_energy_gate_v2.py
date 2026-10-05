"""Strict quasi-static outer energy gate for the Griffith front-step kernel.

This isolated v2 wrapper reuses the v1 geometry proposal and immutable surface
types, but it deliberately does not call v1 ``confirm_candidate``.  The v1
``delta_potential_J`` field has no frozen energy convention, so this module
does not bridge it.  The only supported balance is

    residual = external_work_increment - delta_stored_elastic_energy - Gc*dA

for a quasi-static elastic step with no inertia and no other dissipation.  It
records ``delta_total_potential = delta_stored_elastic - external_work`` for
clarity, but never subtracts that quantity in the residual a second time.

The gate validates numeric completeness, units, state/hash bindings, and the
caller-declared tolerance arithmetic.  It cannot validate the physical origin
or adequacy of caller receipts.  External-path commits remain explicitly
experimental and are never called physically validated.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import hashlib
import json
import math
import re
import sys
from pathlib import Path
from typing import Any

import numpy as np


_V1_DIR = Path(__file__).resolve().parents[1] / "griffith_front_step_20261006"
if str(_V1_DIR) not in sys.path:
    sys.path.insert(0, str(_V1_DIR))
import griffith_front_step as front  # noqa: E402


ENERGY_UNIT = "J"
AREA_UNIT = "m^2"
SUPPORTED_PHYSICS_MODE = "QUASISTATIC_ELASTIC_NO_INERTIA_NO_OTHER_DISSIPATION"
SUPPORTED_ENERGY_SEMANTICS = "DELTA_STORED_ELASTIC_PLUS_EXTERNAL_WORK_SEPARATE"
HASH_PATTERN = re.compile(r"^[0-9a-f]{64}$")
EXPERIMENTAL_QUALIFICATION = "EXPERIMENTAL_UNCALIBRATED"


def _hash_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _hash_payload(payload: Any) -> str:
    def convert(value: object) -> Any:
        if isinstance(value, np.generic):
            return value.item()
        raise TypeError(f"unsupported hash payload value {type(value).__name__}")

    raw = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), allow_nan=False, default=convert
    ).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _valid_hash(value: object) -> bool:
    return isinstance(value, str) and bool(HASH_PATTERN.fullmatch(value))


def _number(value: object) -> float | None:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, float, np.integer, np.floating)):
        return None
    result = float(value)
    return result if math.isfinite(result) else None


def _maybe_number(value: object) -> float | None:
    result = _number(value)
    if result is not None:
        return result
    return None


@dataclass(frozen=True, slots=True)
class QuasiStaticEnergyEvidence:
    """Candidate-bound input to the v2 gate.

    ``caller_provenance`` is serialized as metadata only.  No provenance,
    review, or status string replaces the explicit numeric terms and tolerance.
    """

    candidate_id: str | None
    mechanical_status: str | None
    evidence_origin: str = "EXTERNAL_RECEIPT"
    physics_mode: str | None = SUPPORTED_PHYSICS_MODE
    energy_terms_semantics: str | None = SUPPORTED_ENERGY_SEMANTICS
    energy_unit: str | None = ENERGY_UNIT
    area_unit: str | None = AREA_UNIT
    solver_result_sha256: str | None = None
    source_sha256: str | None = None
    bound_parent_surface_state_sha256: str | None = None
    bound_candidate_geometry_sha256: str | None = None
    reanalyzed_surface_state_sha256: str | None = None
    energy_ledger_sha256: str | None = None
    delta_stored_elastic_energy_J: float | None = None
    external_work_increment_J: float | None = None
    Gc_delta_area_J: float | None = None
    area_increment_m2: float | None = None
    abs_tol_J: float | None = None
    relative_tol: float | None = None
    reference_energy_J: float | None = None
    tolerance_basis: str | None = None
    tolerance_basis_sha256: str | None = None
    # Supplied only to detect and reject accidental passage of the old,
    # ambiguous total-potential quantity as if it were stored elastic energy.
    delta_total_potential_J: float | None = None
    caller_provenance: tuple[tuple[str, str], ...] = ()

    def to_payload(self) -> dict[str, Any]:
        return {
            "candidate_id": self.candidate_id,
            "mechanical_status": self.mechanical_status,
            "evidence_origin": self.evidence_origin,
            "physics_mode": self.physics_mode,
            "energy_terms_semantics": self.energy_terms_semantics,
            "energy_unit": self.energy_unit,
            "area_unit": self.area_unit,
            "solver_result_sha256": self.solver_result_sha256,
            "source_sha256": self.source_sha256,
            "bound_parent_surface_state_sha256": self.bound_parent_surface_state_sha256,
            "bound_candidate_geometry_sha256": self.bound_candidate_geometry_sha256,
            "reanalyzed_surface_state_sha256": self.reanalyzed_surface_state_sha256,
            "energy_ledger_sha256": self.energy_ledger_sha256,
            "delta_stored_elastic_energy_J": self.delta_stored_elastic_energy_J,
            "external_work_increment_J": self.external_work_increment_J,
            "Gc_delta_area_J": self.Gc_delta_area_J,
            "area_increment_m2": self.area_increment_m2,
            "abs_tol_J": self.abs_tol_J,
            "relative_tol": self.relative_tol,
            "reference_energy_J": self.reference_energy_J,
            "tolerance_basis": self.tolerance_basis,
            "tolerance_basis_sha256": self.tolerance_basis_sha256,
            "delta_total_potential_J": self.delta_total_potential_J,
            "caller_provenance_metadata_only": [list(row) for row in self.caller_provenance],
        }


@dataclass(frozen=True, slots=True)
class EnergyGateEvent:
    event_id: str
    sequence: int
    candidate_id: str | None
    status: str
    reason: str
    evidence_sha256: str | None
    evidence_origin: str | None
    physics_mode: str | None
    energy_terms_semantics: str | None
    energy_unit: str | None
    area_unit: str | None
    delta_stored_elastic_energy_J: float | None
    external_work_increment_J: float | None
    Gc_delta_area_J: float | None
    area_increment_m2: float | None
    residual_W_minus_delta_stored_minus_Gc_delta_area_J: float | None
    abs_tol_J: float | None
    relative_tol: float | None
    reference_energy_J: float | None
    relative_tolerance_contribution_J: float | None
    allowed_abs_residual_J: float | None
    derived_delta_total_potential_J: float | None
    tolerance_basis: str | None
    tolerance_basis_sha256: str | None
    caller_provenance: tuple[tuple[str, str], ...] = ()

    def to_payload(self) -> dict[str, Any]:
        return {
            "event_id": self.event_id,
            "sequence": self.sequence,
            "candidate_id": self.candidate_id,
            "status": self.status,
            "reason": self.reason,
            "evidence_sha256": self.evidence_sha256,
            "evidence_origin": self.evidence_origin,
            "physics_mode": self.physics_mode,
            "energy_terms_semantics": self.energy_terms_semantics,
            "energy_unit": self.energy_unit,
            "area_unit": self.area_unit,
            "delta_stored_elastic_energy_J": self.delta_stored_elastic_energy_J,
            "external_work_increment_J": self.external_work_increment_J,
            "Gc_delta_area_J": self.Gc_delta_area_J,
            "area_increment_m2": self.area_increment_m2,
            "residual_W_minus_delta_stored_minus_Gc_delta_area_J": self.residual_W_minus_delta_stored_minus_Gc_delta_area_J,
            "abs_tol_J": self.abs_tol_J,
            "relative_tol": self.relative_tol,
            "reference_energy_J": self.reference_energy_J,
            "relative_tolerance_contribution_J": self.relative_tolerance_contribution_J,
            "allowed_abs_residual_J": self.allowed_abs_residual_J,
            "derived_delta_total_potential_J": self.derived_delta_total_potential_J,
            "tolerance_basis": self.tolerance_basis,
            "tolerance_basis_sha256": self.tolerance_basis_sha256,
            "caller_provenance_metadata_only": [list(row) for row in self.caller_provenance],
        }


@dataclass(frozen=True, slots=True)
class EnergyGateSnapshot:
    controller: front.ControllerSnapshot
    energy_gate_ledger: tuple[EnergyGateEvent, ...] = ()

    @property
    def active_surface(self) -> front.SurfaceVersion:
        return self.controller.active_surface

    @property
    def pending_candidate_id(self) -> str | None:
        return self.controller.pending_candidate_id

    def to_payload(self, *, include_surfaces: bool = False) -> dict[str, Any]:
        return {
            "schema_version": "pgfrag.griffith-front-energy-gate.v2.experimental.v1",
            "controller": self.controller.to_payload(include_surfaces=include_surfaces),
            "energy_gate_ledger": [event.to_payload() for event in self.energy_gate_ledger],
            "energy_gate_claim": "QUASISTATIC_TOLERANCE_CHECK_ONLY_NOT_PHYSICAL_VALIDATION",
            "legacy_delta_potential_bridge": "NOT_BRIDGED",
        }


def initialize(initial_state: front.ExplicitSurfaceState) -> EnergyGateSnapshot:
    return EnergyGateSnapshot(front.initialize(initial_state))


def propose_step(
    snapshot: EnergyGateSnapshot,
    driving_evidence: front.DrivingEvidence | None,
    *,
    frame: front.FrontFrame,
    requested_advance_m: float | None,
    step_cap_m: float | None,
    bounds: front.DomainBounds,
) -> EnergyGateSnapshot:
    return replace(
        snapshot,
        controller=front.propose_step(
            snapshot.controller,
            driving_evidence,
            frame=frame,
            requested_advance_m=requested_advance_m,
            step_cap_m=step_cap_m,
            bounds=bounds,
        ),
    )


def _candidate(snapshot: front.ControllerSnapshot, candidate_id: str) -> front.Candidate | None:
    return next((candidate for candidate in snapshot.candidates if candidate.candidate_id == candidate_id), None)


def _new_gate_event(
    snapshot: EnergyGateSnapshot,
    candidate_id: str | None,
    status: str,
    reason: str,
    evidence: QuasiStaticEnergyEvidence | None,
    *,
    residual: float | None = None,
    rel_contribution: float | None = None,
    allowed: float | None = None,
    derived_total: float | None = None,
    evidence_sha256: str | None = None,
) -> EnergyGateEvent:
    sequence = len(snapshot.energy_gate_ledger)
    if evidence is not None and evidence_sha256 is None:
        try:
            evidence_sha256 = _hash_payload(evidence.to_payload())
        except (TypeError, ValueError, OverflowError):
            evidence_sha256 = None
    def safe(value: object) -> float | None:
        return _maybe_number(value)
    eid = front._stable_id(
        "qs-energy-event", str(sequence), candidate_id or "", status, reason,
        evidence_sha256 or "",
    )
    return EnergyGateEvent(
        event_id=eid,
        sequence=sequence,
        candidate_id=candidate_id,
        status=status,
        reason=reason,
        evidence_sha256=evidence_sha256,
        evidence_origin=evidence.evidence_origin if evidence else None,
        physics_mode=evidence.physics_mode if evidence else None,
        energy_terms_semantics=evidence.energy_terms_semantics if evidence else None,
        energy_unit=evidence.energy_unit if evidence else None,
        area_unit=evidence.area_unit if evidence else None,
        delta_stored_elastic_energy_J=safe(evidence.delta_stored_elastic_energy_J) if evidence else None,
        external_work_increment_J=safe(evidence.external_work_increment_J) if evidence else None,
        Gc_delta_area_J=safe(evidence.Gc_delta_area_J) if evidence else None,
        area_increment_m2=safe(evidence.area_increment_m2) if evidence else None,
        residual_W_minus_delta_stored_minus_Gc_delta_area_J=residual,
        abs_tol_J=safe(evidence.abs_tol_J) if evidence else None,
        relative_tol=safe(evidence.relative_tol) if evidence else None,
        reference_energy_J=safe(evidence.reference_energy_J) if evidence else None,
        relative_tolerance_contribution_J=rel_contribution,
        allowed_abs_residual_J=allowed,
        derived_delta_total_potential_J=derived_total,
        tolerance_basis=evidence.tolerance_basis if evidence else None,
        tolerance_basis_sha256=evidence.tolerance_basis_sha256 if evidence else None,
        caller_provenance=evidence.caller_provenance if evidence else (),
    )


def _fail(
    snapshot: EnergyGateSnapshot,
    candidate: front.Candidate | None,
    candidate_id: str,
    evidence: QuasiStaticEnergyEvidence | None,
    *,
    status: str,
    reason: str,
    evidence_sha256: str | None = None,
    residual: float | None = None,
    rel_contribution: float | None = None,
    allowed: float | None = None,
    derived_total: float | None = None,
) -> EnergyGateSnapshot:
    controller = snapshot.controller
    if candidate is not None and controller.pending_candidate_id == candidate_id:
        if status == "NEEDS_RECOMPUTE":
            controller = front._mark_candidate(controller, candidate, status="NEEDS_RECOMPUTE")
        controller = front._append_event(
            controller,
            kind="CONFIRM_QUASISTATIC_ENERGY_GATE",
            status=status,
            candidate_id=candidate_id,
            reason=reason,
            evidence_hashes=tuple(
                value for value in (evidence_sha256, candidate.geometry_sha256) if value
            ),
        )
    else:
        controller = front._append_event(
            controller,
            kind="CONFIRM_QUASISTATIC_ENERGY_GATE",
            status=status,
            candidate_id=candidate_id,
            reason="candidate_is_not_pending;" + reason,
        )
    event = _new_gate_event(
        snapshot, candidate_id, status, reason, evidence,
        residual=residual, rel_contribution=rel_contribution,
        allowed=allowed, derived_total=derived_total,
        evidence_sha256=evidence_sha256,
    )
    return EnergyGateSnapshot(controller, snapshot.energy_gate_ledger + (event,))


def confirm_candidate(
    snapshot: EnergyGateSnapshot,
    candidate_id: str,
    *,
    reanalyzed_surface: front.ExplicitSurfaceState | None,
    evidence: QuasiStaticEnergyEvidence | None,
) -> EnergyGateSnapshot:
    """Confirm only after numeric quasi-static energy closure is within tolerance."""

    controller = snapshot.controller
    candidate = _candidate(controller, candidate_id)
    if candidate is None or controller.pending_candidate_id != candidate_id:
        return _fail(
            snapshot, candidate, candidate_id, evidence, status="REJECTED_UNKNOWN_CANDIDATE",
            reason="candidate_is_not_pending",
        )
    if candidate.status not in ("PROPOSED_WAITING_REANALYSIS", "NEEDS_RECOMPUTE"):
        return _fail(
            snapshot, candidate, candidate_id, evidence, status="REJECTED_INVALID_CANDIDATE_STATE",
            reason="candidate_state_not_confirmable",
        )
    active = controller.active_surface
    if active.version_id != candidate.parent_version_id or active.state_sha256 != candidate.parent_state_sha256:
        return _fail(
            snapshot, candidate, candidate_id, evidence, status="NEEDS_RECOMPUTE",
            reason="candidate_parent_is_no_longer_active",
        )
    if evidence is None:
        return _fail(
            snapshot, candidate, candidate_id, None, status="NEEDS_RECOMPUTE",
            reason="missing_quasistatic_energy_and_reanalysis_evidence",
        )
    try:
        evidence_sha = _hash_payload(evidence.to_payload())
    except (TypeError, ValueError, OverflowError):
        return _fail(
            snapshot, candidate, candidate_id, evidence, status="NEEDS_RECOMPUTE",
            reason="invalid_or_nonfinite_energy_evidence_payload",
        )

    if evidence.physics_mode != SUPPORTED_PHYSICS_MODE:
        return _fail(
            snapshot, candidate, candidate_id, evidence, status="UNKNOWN",
            reason="dynamic_or_mixed_dissipative_physics_mode_unsupported",
            evidence_sha256=evidence_sha,
        )
    if evidence.delta_total_potential_J is not None:
        return _fail(
            snapshot, candidate, candidate_id, evidence, status="NEEDS_RECOMPUTE",
            reason="total_potential_change_is_not_stored_elastic_energy_and_legacy_field_is_not_bridged",
            evidence_sha256=evidence_sha,
        )
    if evidence.energy_terms_semantics != SUPPORTED_ENERGY_SEMANTICS:
        return _fail(
            snapshot, candidate, candidate_id, evidence, status="NEEDS_RECOMPUTE",
            reason="energy_term_semantics_missing_or_incompatible",
            evidence_sha256=evidence_sha,
        )

    if (evidence.mechanical_status, evidence.evidence_origin) not in (
        ("REANALYZED_CURRENT", "EXTERNAL_RECEIPT"),
        ("EXPLICIT_FIXTURE_ONLY", "EXPLICIT_FIXTURE_NOT_SOLVER"),
    ):
        return _fail(
            snapshot, candidate, candidate_id, evidence, status="NEEDS_RECOMPUTE",
            reason="missing_or_uncoupled_candidate_reanalysis",
            evidence_sha256=evidence_sha,
        )
    if evidence.candidate_id != candidate_id:
        return _fail(
            snapshot, candidate, candidate_id, evidence, status="NEEDS_RECOMPUTE",
            reason="candidate_id_binding_mismatch",
            evidence_sha256=evidence_sha,
        )
    if evidence.bound_parent_surface_state_sha256 != candidate.parent_state_sha256:
        return _fail(
            snapshot, candidate, candidate_id, evidence, status="NEEDS_RECOMPUTE",
            reason="parent_surface_hash_binding_mismatch",
            evidence_sha256=evidence_sha,
        )
    if evidence.bound_candidate_geometry_sha256 != candidate.geometry_sha256:
        return _fail(
            snapshot, candidate, candidate_id, evidence, status="NEEDS_RECOMPUTE",
            reason="candidate_geometry_hash_binding_mismatch",
            evidence_sha256=evidence_sha,
        )
    if not all(_valid_hash(value) for value in (
        evidence.solver_result_sha256,
        evidence.source_sha256,
        evidence.reanalyzed_surface_state_sha256,
        evidence.energy_ledger_sha256,
        evidence.tolerance_basis_sha256,
    )):
        return _fail(
            snapshot, candidate, candidate_id, evidence, status="NEEDS_RECOMPUTE",
            reason="missing_or_invalid_solve_source_energy_or_tolerance_basis_hash",
            evidence_sha256=evidence_sha,
        )
    if not isinstance(evidence.tolerance_basis, str) or not evidence.tolerance_basis.strip():
        return _fail(
            snapshot, candidate, candidate_id, evidence, status="NEEDS_RECOMPUTE",
            reason="missing_tolerance_basis",
            evidence_sha256=evidence_sha,
        )
    if evidence.energy_unit != ENERGY_UNIT or evidence.area_unit != AREA_UNIT:
        return _fail(
            snapshot, candidate, candidate_id, evidence, status="NEEDS_RECOMPUTE",
            reason="missing_or_invalid_energy_or_area_units",
            evidence_sha256=evidence_sha,
        )

    abs_tol = _number(evidence.abs_tol_J)
    rel_tol = _number(evidence.relative_tol)
    reference = _number(evidence.reference_energy_J)
    if abs_tol is None or rel_tol is None or reference is None:
        return _fail(
            snapshot, candidate, candidate_id, evidence, status="NEEDS_RECOMPUTE",
            reason="missing_or_nonfinite_numeric_tolerance_parameters",
            evidence_sha256=evidence_sha,
        )
    if abs_tol < 0 or rel_tol < 0 or reference <= 0:
        return _fail(
            snapshot, candidate, candidate_id, evidence, status="NEEDS_RECOMPUTE",
            reason="tolerances_must_be_nonnegative_and_reference_energy_positive",
            evidence_sha256=evidence_sha,
        )
    stored = _number(evidence.delta_stored_elastic_energy_J)
    work = _number(evidence.external_work_increment_J)
    gc_area = _number(evidence.Gc_delta_area_J)
    area = _number(evidence.area_increment_m2)
    if stored is None or work is None or gc_area is None or area is None:
        return _fail(
            snapshot, candidate, candidate_id, evidence, status="NEEDS_RECOMPUTE",
            reason="missing_or_nonfinite_energy_ledger_terms",
            evidence_sha256=evidence_sha,
        )
    if gc_area <= 0 or area <= 0:
        return _fail(
            snapshot, candidate, candidate_id, evidence, status="NEEDS_RECOMPUTE",
            reason="fracture_energy_and_area_increment_must_be_positive",
            evidence_sha256=evidence_sha,
        )
    if not math.isclose(area, candidate.strip.area_m2, rel_tol=1e-12, abs_tol=1e-18):
        return _fail(
            snapshot, candidate, candidate_id, evidence, status="NEEDS_RECOMPUTE",
            reason="area_increment_does_not_match_proposed_geometry",
            evidence_sha256=evidence_sha,
        )
    expected_gc_area = candidate.Gc_J_per_m2 * candidate.strip.area_m2
    if not math.isclose(gc_area, expected_gc_area, rel_tol=1e-12, abs_tol=1e-18):
        return _fail(
            snapshot, candidate, candidate_id, evidence, status="NEEDS_RECOMPUTE",
            reason="Gc_delta_area_does_not_match_candidate_Gc_and_area",
            evidence_sha256=evidence_sha,
        )

    relative_contribution = rel_tol * reference
    allowed = abs_tol + relative_contribution
    residual = work - stored - gc_area
    derived_total_potential = stored - work
    if not all(math.isfinite(value) for value in (relative_contribution, allowed, residual, derived_total_potential)):
        return _fail(
            snapshot, candidate, candidate_id, evidence, status="NEEDS_RECOMPUTE",
            reason="nonfinite_residual_or_tolerance_limit",
            evidence_sha256=evidence_sha,
        )
    if abs(residual) > allowed:
        return _fail(
            snapshot, candidate, candidate_id, evidence, status="NEEDS_RECOMPUTE",
            reason="quasistatic_energy_residual_exceeds_explicit_tolerance",
            evidence_sha256=evidence_sha,
            residual=residual,
            rel_contribution=relative_contribution,
            allowed=allowed,
            derived_total=derived_total_potential,
        )
    if reanalyzed_surface is None:
        return _fail(
            snapshot, candidate, candidate_id, evidence, status="NEEDS_RECOMPUTE",
            reason="missing_reanalyzed_surface_state",
            evidence_sha256=evidence_sha,
            residual=residual,
            rel_contribution=relative_contribution,
            allowed=allowed,
            derived_total=derived_total_potential,
        )
    try:
        state_sha = _hash_text(reanalyzed_surface.to_json())
    except (TypeError, ValueError, OverflowError):
        state_sha = ""
    if state_sha != evidence.reanalyzed_surface_state_sha256:
        return _fail(
            snapshot, candidate, candidate_id, evidence, status="NEEDS_RECOMPUTE",
            reason="reanalyzed_surface_state_hash_mismatch",
            evidence_sha256=evidence_sha,
            residual=residual,
            rel_contribution=relative_contribution,
            allowed=allowed,
            derived_total=derived_total_potential,
        )
    if not front._state_geometry_matches_candidate(reanalyzed_surface, active, candidate):
        return _fail(
            snapshot, candidate, candidate_id, evidence, status="NEEDS_RECOMPUTE",
            reason="reanalyzed_surface_geometry_or_ownership_mismatch",
            evidence_sha256=evidence_sha,
            residual=residual,
            rel_contribution=relative_contribution,
            allowed=allowed,
            derived_total=derived_total_potential,
        )

    fixture_only = evidence.evidence_origin == "EXPLICIT_FIXTURE_NOT_SOLVER"
    commit_status = "COMMITTED_FIXTURE_ONLY" if fixture_only else "COMMITTED_EXPERIMENTAL_REANALYSIS"
    version_id = front._stable_id("surface-version-v2", candidate.candidate_id, state_sha)
    version = front.SurfaceVersion(
        version_id=version_id,
        parent_version_id=active.version_id,
        state_sha256=state_sha,
        state=reanalyzed_surface,
        qualification=EXPERIMENTAL_QUALIFICATION,
        version_status=commit_status,
    )
    committed_candidate = replace(
        candidate,
        status=commit_status,
        evidence_origin=evidence.evidence_origin,
        # The v1 field is intentionally not bridged. Full semantics live in
        # this module's separate typed energy-gate event ledger.
        energy_ledger=(),
    )
    candidates = tuple(
        committed_candidate if item.candidate_id == candidate_id else item
        for item in controller.candidates
    )
    updated_controller = replace(
        controller,
        committed_versions=controller.committed_versions + (version,),
        active_surface_version_id=version_id,
        candidates=candidates,
        pending_candidate_id=None,
    )
    updated_controller = front._append_event(
        updated_controller,
        kind="CONFIRM_QUASISTATIC_ENERGY_GATE",
        status=commit_status,
        candidate_id=candidate_id,
        reason="explicit_numeric_energy_gate_passed; no_physical_validation_claim",
        evidence_hashes=(evidence_sha, state_sha, candidate.geometry_sha256),
        active_surface_version_id=version_id,
    )
    event = _new_gate_event(
        snapshot,
        candidate_id,
        commit_status,
        "explicit_numeric_energy_gate_passed; no_physical_validation_claim",
        evidence,
        residual=residual,
        rel_contribution=relative_contribution,
        allowed=allowed,
        derived_total=derived_total_potential,
        evidence_sha256=evidence_sha,
    )
    return EnergyGateSnapshot(updated_controller, snapshot.energy_gate_ledger + (event,))


def reject_candidate(snapshot: EnergyGateSnapshot, candidate_id: str, reason: str) -> EnergyGateSnapshot:
    return replace(snapshot, controller=front.reject_candidate(snapshot.controller, candidate_id, reason))


def make_fixture_evidence(
    candidate: front.Candidate,
    parent_state_sha256: str,
    reanalyzed_surface: front.ExplicitSurfaceState,
    *,
    origin: str = "EXPLICIT_FIXTURE_NOT_SOLVER",
    mechanical_status: str = "EXPLICIT_FIXTURE_ONLY",
    delta_stored_elastic_energy_J: float | None = 0.10,
    external_work_increment_J: float | None = 0.13,
    abs_tol_J: float | None = 0.0,
    relative_tol: float | None = 0.0,
    reference_energy_J: float | None = 1.0,
    physics_mode: str = SUPPORTED_PHYSICS_MODE,
    tolerance_basis: str | None = "explicit unit-test tolerance fixture; not a validated numerical tolerance",
    delta_total_potential_J: float | None = None,
) -> QuasiStaticEnergyEvidence:
    """Synthetic test evidence with conspicuous fixture hashes and labels."""

    digest = lambda item: _hash_text("explicit-test-fixture-not-solver:" + item)
    area = candidate.strip.area_m2
    return QuasiStaticEnergyEvidence(
        candidate_id=candidate.candidate_id,
        mechanical_status=mechanical_status,
        evidence_origin=origin,
        physics_mode=physics_mode,
        energy_terms_semantics=SUPPORTED_ENERGY_SEMANTICS,
        solver_result_sha256=digest("solver-result"),
        source_sha256=digest("source"),
        bound_parent_surface_state_sha256=parent_state_sha256,
        bound_candidate_geometry_sha256=candidate.geometry_sha256,
        reanalyzed_surface_state_sha256=_hash_text(reanalyzed_surface.to_json()),
        energy_ledger_sha256=digest("energy-ledger"),
        delta_stored_elastic_energy_J=delta_stored_elastic_energy_J,
        external_work_increment_J=external_work_increment_J,
        Gc_delta_area_J=candidate.Gc_J_per_m2 * area,
        area_increment_m2=area,
        abs_tol_J=abs_tol_J,
        relative_tol=relative_tol,
        reference_energy_J=reference_energy_J,
        tolerance_basis=tolerance_basis,
        tolerance_basis_sha256=digest("tolerance-basis"),
        delta_total_potential_J=delta_total_potential_J,
        caller_provenance=(("fixture", "explicit-only"), ("solver_invoked", "false")),
    )

