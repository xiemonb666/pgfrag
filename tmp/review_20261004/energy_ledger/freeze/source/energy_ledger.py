"""Strict, source-aware energy-history audit primitives for PG fracture work.

This module audits supplied histories only.  It does not run a solver, establish
that a caller's evidence is scientifically sound, or calibrate material data.

Native-source observations that motivate the contract (all citations are local
to this workspace):

* Akantu v5.0.7 ``MaterialCohesive::getEnergy`` exposes ``reversible``,
  ``dissipated``, and ``cohesive contact``.  Unknown names return 0 in that
  source, so a missing/unchecked ``contact`` selector cannot be treated as a
  measured zero.  See
  ``third_party/akantu-runtime/source-v5.0.7/src/model/solid_mechanics/solid_mechanics_model_cohesive/materials/material_cohesive.cc``.
* Akantu's model ``potential`` is material energy; ``external work`` is a
  per-call contribution in the inspected source.  The existing report says the
  extrinsic runner instead integrated sampled boundary reaction and imposed
  displacement.  See ``relicweave-synlab/docs/a12/COHESIVE_CAPABILITY_20261004.md``.
* Peridynamics.jl exports OSB point strain-energy density; the PG fixture sums
  that field against point volumes and separately sums kinetic energy.  The
  fixture itself labels total conservation UNKNOWN because imposed boundary
  velocities do external work and no complete ledger exists.  See
  ``pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/physics/ordinary_state_based.jl``
  and ``pgfrag/src/pgfrag/resources/julia/run_fixture.jl``.
* The native short-range PD contact implementation assembles force density;
  the inspected code does not export a contact-energy channel.  See
  ``pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/physics/short_range_force_contact.jl``.

The explicit ``Qualification`` and ``EvidenceRef`` inputs are caller claims
with provenance pointers.  Their presence is not itself proof; free-energy
closure is refused unless every required channel, the physical clock, the
boundary-work definition, and an independent measurement basis are supplied.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Mapping, Sequence


_ENERGY_TO_J = {
    "J": 1.0,
    "kJ": 1.0e3,
    "mJ": 1.0e-3,
    "uJ": 1.0e-6,
    "µJ": 1.0e-6,
    "nJ": 1.0e-9,
    "Wh": 3600.0,
    "kWh": 3.6e6,
}
_TIME_TO_S = {"s": 1.0, "ms": 1.0e-3, "us": 1.0e-6, "µs": 1.0e-6, "ns": 1.0e-9}
_ROLES = {"kinetic", "stored", "dissipated", "external_work", "diagnostic"}
_REGIMES = {"free_dynamics", "prescribed_boundary_dynamics", "all_dof_prescribed_constitutive"}
_STATUS_ONLY_METHODS = {"pass", "fail", "unknown", "verified", "not applicable", "n/a"}


@dataclass(frozen=True)
class Clock:
    """Ordered sample coordinates plus their clock definition and provenance."""

    values: Sequence[float]
    unit: str | None
    definition: str | None
    source: str | None
    identity: str | None


@dataclass(frozen=True)
class EnergyChannel:
    """One signed raw energy series with explicit semantic and source metadata.

    ``key`` in ``EnergyHistory.channels`` is the caller's canonical semantic
    name. ``source_key`` is the exact upstream getter/field name; this module
    never aliases or falls back between names.
    """

    values: Sequence[float]
    unit: str
    definition: str | None
    source: str | None
    source_key: str | None
    clock_id: str | None
    role: str | None
    cumulative: bool = False


@dataclass(frozen=True)
class EvidenceRef:
    """Structured evidence pointer; a free-form status such as ``PASS`` is not one."""

    source: str
    method: str


@dataclass(frozen=True)
class BoundaryWorkDefinition:
    """Definition for cumulative external work supplied through boundaries."""

    boundary: str
    force_or_reaction: str
    displacement_or_velocity: str
    integration_rule: str
    sign_convention: str
    source: str
    clock_id: str


@dataclass(frozen=True)
class EnergyTolerance:
    """Predeclared thresholds, all explicit and expressed in joules."""

    negative_dissipation_abs_j: float
    negative_dissipation_definition: str
    negative_dissipation_source: str
    balance_abs_j: float
    balance_rel: float
    balance_definition: str
    balance_source: str
    constitutive_relative_max: float | None = None
    constitutive_tolerance_definition: str | None = None
    constitutive_tolerance_source: str | None = None


@dataclass(frozen=True)
class EnergyHistory:
    clock: Clock
    channels: Mapping[str, EnergyChannel]
    tolerance: EnergyTolerance


@dataclass(frozen=True)
class Qualification:
    """Caller-supplied physical scope and evidence references.

    ``required_balance_channels`` is the explicit accounting contract for the
    case. It must include every term required by that model, including any
    contact channel that has not been proved absent. ``verified_channels`` is
    keyed by those exact canonical names.
    """

    regime: str
    source: str | None
    method: str | None
    required_balance_channels: Sequence[str]
    verified_channels: Mapping[str, EvidenceRef]
    independent_balance_basis: EvidenceRef | None
    boundary_work: BoundaryWorkDefinition | None
    constitutive_target_j: float | None = None
    constitutive_target_definition: str | None = None
    constitutive_target_source: str | None = None


def _nonempty(value: object) -> bool:
    return isinstance(value, str) and bool(value.strip())


def _valid_evidence(value: object) -> bool:
    if not isinstance(value, EvidenceRef):
        return False
    if not _nonempty(value.source) or not _nonempty(value.method):
        return False
    source = value.source.strip().casefold()
    method = value.method.strip().casefold()
    return source not in _STATUS_ONLY_METHODS and method not in _STATUS_ONLY_METHODS


def _as_finite_tuple(values: Sequence[float], label: str) -> tuple[float, ...]:
    if isinstance(values, (str, bytes)):
        raise ValueError(f"{label} must be a finite numeric sequence")
    try:
        result = tuple(float(value) for value in values)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be a finite numeric sequence") from exc
    if not result:
        raise ValueError(f"{label} must contain at least one sample")
    if any(not math.isfinite(value) for value in result):
        raise ValueError(f"{label} contains a non-finite sample")
    return result


def _validate_inputs(history: EnergyHistory) -> tuple[tuple[float, ...], dict[str, tuple[float, ...]]]:
    if not isinstance(history, EnergyHistory):
        raise TypeError("history must be an EnergyHistory")
    times = _as_finite_tuple(history.clock.values, "clock.values")
    if any(right <= left for left, right in zip(times, times[1:])):
        raise ValueError("clock.values must be strictly increasing")
    converted: dict[str, tuple[float, ...]] = {}
    for name, channel in history.channels.items():
        if not isinstance(name, str) or not name:
            raise ValueError("channel names must be non-empty strings")
        if not isinstance(channel, EnergyChannel):
            raise TypeError(f"channel {name!r} must be an EnergyChannel")
        if channel.unit not in _ENERGY_TO_J:
            raise ValueError(
                f"channel {name!r} unit {channel.unit!r} is not a supported energy unit convertible to J"
            )
        values = _as_finite_tuple(channel.values, f"channel {name!r}.values")
        if len(values) != len(times):
            raise ValueError(f"channel {name!r} sample count does not match the clock")
        if channel.role is not None and channel.role not in _ROLES:
            raise ValueError(f"channel {name!r} has unsupported role {channel.role!r}")
        converted[name] = tuple(value * _ENERGY_TO_J[channel.unit] for value in values)
    _validate_tolerance(history.tolerance)
    return times, converted


def _validate_tolerance(tolerance: EnergyTolerance) -> None:
    if not isinstance(tolerance, EnergyTolerance):
        raise TypeError("history.tolerance must be an EnergyTolerance")
    numeric = (
        ("negative_dissipation_abs_j", tolerance.negative_dissipation_abs_j),
        ("balance_abs_j", tolerance.balance_abs_j),
        ("balance_rel", tolerance.balance_rel),
    )
    for name, value in numeric:
        if not math.isfinite(float(value)) or float(value) < 0:
            raise ValueError(f"tolerance {name} must be finite and non-negative")
    for name, value in (
        ("negative_dissipation_definition", tolerance.negative_dissipation_definition),
        ("negative_dissipation_source", tolerance.negative_dissipation_source),
        ("balance_definition", tolerance.balance_definition),
        ("balance_source", tolerance.balance_source),
    ):
        if not _nonempty(value):
            raise ValueError(f"tolerance {name} must be explicit")
    if tolerance.constitutive_relative_max is not None:
        if not math.isfinite(float(tolerance.constitutive_relative_max)) or tolerance.constitutive_relative_max < 0:
            raise ValueError("constitutive_relative_max must be finite and non-negative")
        if not _nonempty(tolerance.constitutive_tolerance_definition):
            raise ValueError("constitutive tolerance definition is required when its threshold is supplied")
        if not _nonempty(tolerance.constitutive_tolerance_source):
            raise ValueError("constitutive tolerance source is required when its threshold is supplied")


def _channel_metadata_status(name: str, channel: EnergyChannel, clock: Clock) -> tuple[str, list[str]]:
    missing: list[str] = []
    if not _nonempty(channel.definition):
        missing.append("definition")
    if not _nonempty(channel.source):
        missing.append("source")
    if not _nonempty(channel.source_key):
        missing.append("source_key")
    if not _nonempty(clock.identity) or not _nonempty(channel.clock_id):
        missing.append("clock_identity")
    elif channel.clock_id != clock.identity:
        missing.append("clock_mismatch")
    return ("KNOWN" if not missing else "UNKNOWN"), missing


def _clock_is_physical(clock: Clock) -> bool:
    return (
        clock.unit in _TIME_TO_S
        and _nonempty(clock.definition)
        and _nonempty(clock.source)
        and _nonempty(clock.identity)
    )


def _boundary_work_is_defined(work: BoundaryWorkDefinition | None, clock: Clock) -> bool:
    if not isinstance(work, BoundaryWorkDefinition):
        return False
    return all(
        _nonempty(item)
        for item in (
            work.boundary,
            work.force_or_reaction,
            work.displacement_or_velocity,
            work.integration_rule,
            work.source,
            work.clock_id,
        )
    ) and work.sign_convention == "positive_into_system" and work.clock_id == clock.identity


def _dissipation_audit(
    history: EnergyHistory, converted: Mapping[str, tuple[float, ...]]
) -> dict[str, object]:
    tolerance = history.tolerance.negative_dissipation_abs_j
    records: list[dict[str, object]] = []
    failures: list[dict[str, object]] = []
    roundoff: list[dict[str, object]] = []
    dissipated_names = [
        name for name, channel in history.channels.items() if channel.role == "dissipated"
    ]
    if not dissipated_names:
        return {
            "status": "UNKNOWN",
            "reason": "no channel with role='dissipated' was supplied",
            "negative_tolerance_J": tolerance,
            "tolerance_definition": history.tolerance.negative_dissipation_definition,
            "tolerance_source": history.tolerance.negative_dissipation_source,
            "channels": [],
        }
    for name in dissipated_names:
        channel = history.channels[name]
        vals_j = converted[name]
        minimum = min(vals_j)
        negatives = [
            {"sample_index": index, "raw_value": channel.values[index], "value_J": value}
            for index, value in enumerate(vals_j)
            if value < 0.0
        ]
        decreases = []
        if channel.cumulative:
            for index, (before, after) in enumerate(zip(vals_j, vals_j[1:]), start=1):
                if after - before < -tolerance:
                    decreases.append(
                        {
                            "sample_index": index,
                            "previous_value_J": before,
                            "raw_value": channel.values[index],
                            "value_J": after,
                        }
                    )
        record = {
            "channel": name,
            "source_key": channel.source_key,
            "minimum_raw_value": channel.values[vals_j.index(minimum)],
            "minimum_J": minimum,
            "negative_samples": negatives,
            "cumulative_decreases_beyond_tolerance": decreases,
            "cumulative": channel.cumulative,
        }
        records.append(record)
        if minimum < -tolerance or decreases:
            failures.append(record)
        elif minimum < 0.0:
            roundoff.append(record)
    return {
        "status": "FAIL" if failures else ("PASS_WITH_ROUNDOFF" if roundoff else "PASS"),
        "negative_tolerance_J": tolerance,
        "tolerance_definition": history.tolerance.negative_dissipation_definition,
        "tolerance_source": history.tolerance.negative_dissipation_source,
        "raw_values_are_preserved": True,
        "raw_values_are_clamped_or_abs_transformed": False,
        "channels": records,
        "failing_channels": [record["channel"] for record in failures],
        "roundoff_only_channels": [record["channel"] for record in roundoff],
    }


def _constitutive_comparison(
    history: EnergyHistory,
    qualification: Qualification,
    converted: Mapping[str, tuple[float, ...]],
) -> dict[str, object]:
    target = qualification.constitutive_target_j
    if target is None:
        return {"status": "UNKNOWN", "reason": "no constitutive target was supplied"}
    if not math.isfinite(float(target)) or target <= 0:
        return {"status": "UNKNOWN", "reason": "constitutive target must be finite and positive"}
    target_metadata_ok = (
        _nonempty(qualification.constitutive_target_definition)
        and _nonempty(qualification.constitutive_target_source)
        and history.tolerance.constitutive_relative_max is not None
    )
    dissipated = [
        name for name, channel in history.channels.items() if channel.role == "dissipated"
    ]
    if len(dissipated) != 1:
        return {
            "status": "UNKNOWN",
            "reason": "constitutive comparison requires exactly one readable dissipated channel",
        }
    name = dissipated[0]
    channel_metadata_status, missing_channel_metadata = _channel_metadata_status(
        name, history.channels[name], history.clock
    )
    final_j = converted[name][-1]
    relative_error = abs(final_j - float(target)) / float(target)
    tolerance = history.tolerance.constitutive_relative_max
    channel_evidence_ok = _valid_evidence(qualification.verified_channels.get(name))
    if not target_metadata_ok or channel_metadata_status != "KNOWN" or not channel_evidence_ok:
        status = "UNKNOWN"
    else:
        status = "MATCHED_REFERENCE" if relative_error <= tolerance else "OUTSIDE_REFERENCE"
    return {
        "status": status,
        "scope": "constitutive comparison only; not free-energy closure or material calibration",
        "channel": name,
        "final_raw_value": history.channels[name].values[-1],
        "final_J": final_j,
        "target_J": float(target),
        "relative_error": relative_error,
        "relative_tolerance": tolerance,
        "target_definition": qualification.constitutive_target_definition,
        "target_source": qualification.constitutive_target_source,
        "channel_metadata_status": channel_metadata_status,
        "missing_or_mismatched_channel_metadata": missing_channel_metadata,
        "channel_semantics_evidence_status": (
            "VERIFIED_BY_CALLER_REFERENCE" if channel_evidence_ok else "UNKNOWN"
        ),
    }


def audit_energy_history(history: EnergyHistory, qualification: Qualification) -> dict[str, object]:
    """Audit finite ordered signed-energy histories without filling missing terms.

    The energy-balance residual uses the first supplied sample as the explicit
    baseline and computes, at each later sample,

    ``Δ(kinetic + stored) + Δdissipated - Δexternal_work``.

    A free-dynamics result is calculated only when every required channel is
    present, has matching explicit clock/provenance metadata and an
    ``EvidenceRef`` supplied for its exact key, with qualified independent
    measurement basis and boundary work.  Evidence references remain caller
    attestations; this function cannot validate files or experimental method.
    """

    times, converted = _validate_inputs(history)
    if not isinstance(qualification, Qualification):
        raise TypeError("qualification must be a Qualification")
    if qualification.regime not in _REGIMES:
        raise ValueError(f"unsupported qualification regime {qualification.regime!r}")

    channel_results: dict[str, dict[str, object]] = {}
    for name, channel in history.channels.items():
        metadata_status, missing_metadata = _channel_metadata_status(name, channel, history.clock)
        evidence_status = "VERIFIED_BY_CALLER_REFERENCE" if _valid_evidence(qualification.verified_channels.get(name)) else "UNKNOWN"
        channel_results[name] = {
            "status": "READ",
            "role": channel.role,
            "unit": channel.unit,
            "values_raw": list(channel.values),
            "values_J": list(converted[name]),
            "definition": channel.definition,
            "source": channel.source,
            "source_key": channel.source_key,
            "clock_id": channel.clock_id,
            "clock_metadata_status": metadata_status,
            "missing_or_mismatched_metadata": missing_metadata,
            "semantic_evidence_status": evidence_status,
        }
    missing_required = [
        key
        for key in qualification.required_balance_channels
        if key not in history.channels
    ]
    for key in missing_required:
        channel_results[key] = {
            "status": "UNKNOWN",
            "reason": "required channel was not supplied; no zero/default value was substituted",
            "role": "UNKNOWN",
            "unit": "UNKNOWN",
            "values_raw": None,
            "values_J": None,
            "definition": None,
            "source": None,
            "source_key": None,
            "clock_id": None,
            "clock_metadata_status": "UNKNOWN",
            "missing_or_mismatched_metadata": ["channel missing"],
            "semantic_evidence_status": "UNKNOWN",
        }
    unknown_required = [
        key
        for key in qualification.required_balance_channels
        if key in history.channels
        and (
            channel_results[key]["clock_metadata_status"] != "KNOWN"
            or channel_results[key]["semantic_evidence_status"] != "VERIFIED_BY_CALLER_REFERENCE"
            or history.channels[key].role not in {"kinetic", "stored", "dissipated", "external_work"}
        )
    ]
    dissipation = _dissipation_audit(history, converted)

    baseline = {
        "sample_index": 0,
        "clock_value_raw": times[0],
        "clock_unit": history.clock.unit,
        "policy": "first sample of the supplied ordered window; no hidden zeroing",
        "raw_energy_values": {name: channel.values[0] for name, channel in history.channels.items()},
    }

    if qualification.regime == "all_dof_prescribed_constitutive":
        balance: dict[str, object] = {
            "status": "NOT_APPLICABLE",
            "reason": "all-DOF-prescribed constitutive case; boundary motion is not an independent free-dynamics balance",
        }
    elif qualification.regime != "free_dynamics":
        balance = {
            "status": "UNKNOWN",
            "reason": "free-energy balance is only evaluated for explicitly qualified free_dynamics histories",
        }
    else:
        dynamics_evidence_ok = _valid_evidence(
            EvidenceRef(qualification.source, qualification.method)
            if _nonempty(qualification.source) and _nonempty(qualification.method)
            else None
        )
        independence_ok = _valid_evidence(qualification.independent_balance_basis)
        clock_ok = _clock_is_physical(history.clock)
        work_ok = _boundary_work_is_defined(qualification.boundary_work, history.clock)
        required_role_names = {
            role: [
                name
                for name in qualification.required_balance_channels
                if name in history.channels and history.channels[name].role == role
            ]
            for role in ("kinetic", "stored", "dissipated", "external_work")
        }
        no_missing_roles = all(required_role_names.values())
        problems = []
        if missing_required:
            problems.append("required channels missing: " + ", ".join(missing_required))
        if unknown_required:
            problems.append("required channels lack matching metadata or exact-key evidence: " + ", ".join(unknown_required))
        if not dynamics_evidence_ok:
            problems.append("free-dynamics qualification evidence is missing or status-only")
        if not independence_ok:
            problems.append("independent measurement-basis evidence is missing or status-only")
        if not clock_ok:
            problems.append("physical time unit, clock definition, source, or identity is missing")
        if not work_ok:
            problems.append("boundary external-work definition, sign, integration, source, or clock is missing")
        if not no_missing_roles:
            problems.append("required balance role is missing: " + ", ".join(role for role, names in required_role_names.items() if not names))
        if any(
            name in history.channels
            and history.channels[name].role in {"dissipated", "external_work"}
            and not history.channels[name].cumulative
            for name in qualification.required_balance_channels
        ):
            problems.append("dissipated and external-work balance histories must be cumulative")

        if problems:
            balance = {
                "status": "UNKNOWN",
                "reason": "; ".join(problems),
                "missing_required_channels": missing_required,
                "unknown_required_channels": unknown_required,
            }
        else:
            kinetic_names = required_role_names["kinetic"]
            stored_names = required_role_names["stored"]
            dissipated_names = required_role_names["dissipated"]
            work_names = required_role_names["external_work"]
            kinetic = [sum(converted[name][index] for name in kinetic_names) for index in range(len(times))]
            stored = [sum(converted[name][index] for name in stored_names) for index in range(len(times))]
            dissipated_values = [sum(converted[name][index] for name in dissipated_names) for index in range(len(times))]
            work = [sum(converted[name][index] for name in work_names) for index in range(len(times))]
            delta_kinetic = [value - kinetic[0] for value in kinetic]
            delta_stored = [value - stored[0] for value in stored]
            delta_dissipated = [value - dissipated_values[0] for value in dissipated_values]
            delta_work = [value - work[0] for value in work]
            residuals = [
                delta_kinetic[i] + delta_stored[i] + delta_dissipated[i] - delta_work[i]
                for i in range(len(times))
            ]
            scales = [
                abs(delta_kinetic[i]) + abs(delta_stored[i]) + abs(delta_dissipated[i]) + abs(delta_work[i])
                for i in range(len(times))
            ]
            limits = [max(history.tolerance.balance_abs_j, history.tolerance.balance_rel * scale) for scale in scales]
            failures = [index for index, residual in enumerate(residuals) if abs(residual) > limits[index]]
            worst_index = max(range(len(residuals)), key=lambda index: abs(residuals[index]))
            balance = {
                "status": "FAIL" if failures else "PASS",
                "formula": "Δ(kinetic + stored) + Δdissipated - Δexternal_work",
                "baseline": baseline,
                "kinetic_channels": kinetic_names,
                "stored_channels": stored_names,
                "dissipated_channels": dissipated_names,
                "external_work_channels": work_names,
                "sample_count": len(times),
                "residuals_J": residuals,
                "scales_J": scales,
                "allowed_residuals_J": limits,
                "failing_sample_indices": failures,
                "worst_sample_index": worst_index,
                "worst_residual_J": residuals[worst_index],
                "worst_scale_J": scales[worst_index],
                "tolerance_definition": history.tolerance.balance_definition,
                "tolerance_source": history.tolerance.balance_source,
                "qualification_is_caller_attestation": True,
            }

    constitutive = (
        _constitutive_comparison(history, qualification, converted)
        if qualification.regime == "all_dof_prescribed_constitutive"
        else {"status": "NOT_APPLICABLE", "reason": "no all-DOF constitutive-only target was requested"}
    )
    if dissipation["status"] == "FAIL" or balance["status"] == "FAIL":
        overall = "FAIL"
    elif qualification.regime == "all_dof_prescribed_constitutive":
        overall = "CONSTITUTIVE_ONLY"
    elif balance["status"] == "PASS" and dissipation["status"] in {"PASS", "PASS_WITH_ROUNDOFF"}:
        overall = "PASS_ENERGY_AUDIT_ONLY"
    else:
        overall = "UNKNOWN"

    clock_status = "KNOWN" if _clock_is_physical(history.clock) else "UNKNOWN"
    return {
        "schema_version": "pgfrag.energy-ledger.v1",
        "overall_status": overall,
        "qualification": {
            "regime": qualification.regime,
            "source": qualification.source,
            "method": qualification.method,
            "qualification_is_caller_attestation": True,
            "required_balance_channels": list(qualification.required_balance_channels),
            "independent_balance_basis": (
                {"source": qualification.independent_balance_basis.source, "method": qualification.independent_balance_basis.method}
                if _valid_evidence(qualification.independent_balance_basis)
                else "UNKNOWN"
            ),
        },
        "clock": {
            "status": clock_status,
            "values_raw": list(times),
            "unit": history.clock.unit,
            "definition": history.clock.definition,
            "source": history.clock.source,
            "identity": history.clock.identity,
            "physical_time_scale_to_s": _TIME_TO_S.get(history.clock.unit),
        },
        "baseline": baseline,
        "channels": channel_results,
        "missing_required_channels": missing_required,
        "dissipation": dissipation,
        "free_energy_balance": balance,
        "constitutive_reference_comparison": constitutive,
        "material_calibration": {
            "status": "UNKNOWN",
            "reason": "this ledger does not calibrate material parameters or validate a target material independently",
        },
        "interpretation_limit": "arithmetic audit of caller-supplied histories; not solver validation, physical qualification, or calibration",
    }


__all__ = [
    "BoundaryWorkDefinition",
    "Clock",
    "EnergyChannel",
    "EnergyHistory",
    "EnergyTolerance",
    "EvidenceRef",
    "Qualification",
    "audit_energy_history",
]
