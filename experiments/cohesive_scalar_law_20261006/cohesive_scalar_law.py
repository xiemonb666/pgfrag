"""Bounded scalar Mode-I bilinear cohesive law with an explicit energy ledger.

All inputs and outputs use SI units: opening in metres, traction and strength
in pascals, stiffness in Pa/m, and energy per new crack area in J/m^2.  The
law is scalar opening-only. Negative opening (compression) and every mode other
than ``mode_i_opening`` are rejected as unsupported/UNKNOWN.

The monotonic envelope is elastic up to ``delta0 = tpeak/K0`` and then
softens linearly to zero at ``deltaf = 2*Gc/tpeak``. A maximum-opening history
variable makes damage irreversible. On unloading/reloading the traction is
``(1-d)*K0*delta``. The material state records recoverable stored energy,
irreversible dissipation, and signed external work per unit area.

This is an isolated experiment primitive. Its numerical fixture values are
illustrative and are not a ceramic material calibration.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from numbers import Real


MODE_I_OPENING = "mode_i_opening"


class UnsupportedModeError(NotImplementedError):
    """The requested compression or non-opening mode is outside this law."""


def _positive_si(name: str, value: object) -> float:
    """Accept plain real scalars only; field names define the expected SI unit."""
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a plain real scalar in its named SI unit")
    try:
        result = float(value)
    except (OverflowError, TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be finite and positive") from exc
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return result


def _finite_nonnegative(name: str, value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a plain real scalar in metres")
    try:
        result = float(value)
    except (OverflowError, TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be finite and nonnegative") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite and nonnegative")
    if result < 0.0:
        raise UnsupportedModeError(
            "negative opening is compression; compression response is UNKNOWN and unsupported"
        )
    return result


def _check_mode(mode: object) -> None:
    if mode != MODE_I_OPENING or not isinstance(mode, str):
        raise UnsupportedModeError(
            "only mode_i_opening is implemented; compression and mixed-mode response are UNKNOWN"
        )


@dataclass(frozen=True)
class BilinearCohesiveMaterial:
    """Parameters for a single scalar opening law, in the named SI units."""

    k0_pa_per_m: float
    tpeak_pa: float
    gc_j_per_m2: float

    def __post_init__(self) -> None:
        k0 = _positive_si("k0_pa_per_m", self.k0_pa_per_m)
        tpeak = _positive_si("tpeak_pa", self.tpeak_pa)
        gc = _positive_si("gc_j_per_m2", self.gc_j_per_m2)
        delta0 = tpeak / k0
        deltaf = (2.0 * gc) / tpeak
        if not math.isfinite(delta0) or delta0 <= 0.0:
            raise ValueError("derived delta0=tpeak/K0 must be finite and positive")
        if not math.isfinite(deltaf) or deltaf <= 0.0:
            raise ValueError("derived deltaf=2*Gc/tpeak must be finite and positive")
        if deltaf <= delta0:
            raise ValueError("parameters require deltaf > delta0 (2*Gc/tpeak > tpeak/K0)")
        object.__setattr__(self, "k0_pa_per_m", k0)
        object.__setattr__(self, "tpeak_pa", tpeak)
        object.__setattr__(self, "gc_j_per_m2", gc)

    @property
    def delta0_m(self) -> float:
        return self.tpeak_pa / self.k0_pa_per_m

    @property
    def deltaf_m(self) -> float:
        return (2.0 * self.gc_j_per_m2) / self.tpeak_pa

    @property
    def signature(self) -> tuple[float, float, float]:
        return (self.k0_pa_per_m, self.tpeak_pa, self.gc_j_per_m2)


@dataclass(frozen=True)
class CohesiveState:
    """One immutable state, with all energies expressed per unit area."""

    opening_m: float
    kappa_m: float
    damage: float
    stored_energy_j_per_m2: float
    dissipated_energy_j_per_m2: float
    external_work_j_per_m2: float
    update_count: int
    material_signature: tuple[float, float, float]


def _damage_at(material: BilinearCohesiveMaterial, kappa_m: float) -> float:
    d0, df = material.delta0_m, material.deltaf_m
    if kappa_m <= d0:
        return 0.0
    if kappa_m >= df:
        return 1.0
    envelope = material.tpeak_pa * (df - kappa_m) / (df - d0)
    damage = 1.0 - envelope / (material.k0_pa_per_m * kappa_m)
    # Protect the closed-form value from a one-ulp excursion at either bound.
    return min(1.0, max(0.0, damage))


def _dissipation_at(material: BilinearCohesiveMaterial, kappa_m: float) -> float:
    """Exact cumulative damage dissipation for this bilinear law."""
    d0, df = material.delta0_m, material.deltaf_m
    if kappa_m <= d0:
        return 0.0
    if kappa_m >= df:
        return material.gc_j_per_m2
    return material.gc_j_per_m2 * ((kappa_m - d0) / (df - d0))


def _envelope_traction(material: BilinearCohesiveMaterial, opening_m: float) -> float:
    d0, df = material.delta0_m, material.deltaf_m
    if opening_m <= d0:
        return material.k0_pa_per_m * opening_m
    if opening_m < df:
        return material.tpeak_pa * ((df - opening_m) / (df - d0))
    return 0.0


def _envelope_work(material: BilinearCohesiveMaterial, start_m: float, end_m: float) -> float:
    """Integrate the piecewise-linear monotonic envelope exactly."""
    if end_m <= start_m:
        return 0.0
    cuts = [start_m]
    for corner in (material.delta0_m, material.deltaf_m):
        if start_m < corner < end_m:
            cuts.append(corner)
    cuts.append(end_m)
    total = 0.0
    for left, right in zip(cuts, cuts[1:]):
        tl = _envelope_traction(material, left)
        tr = _envelope_traction(material, right)
        total += 0.5 * (tl + tr) * (right - left)
    if not math.isfinite(total):
        raise FloatingPointError("monotonic envelope work overflowed")
    return total


def _stored_energy(material: BilinearCohesiveMaterial, opening_m: float, damage: float) -> float:
    stiffness = material.k0_pa_per_m * (1.0 - damage)
    if stiffness == 0.0:
        return 0.0
    energy = 0.5 * stiffness * opening_m * opening_m
    if not math.isfinite(energy):
        raise FloatingPointError("stored cohesive energy overflowed")
    return energy


def initial_state(material: BilinearCohesiveMaterial) -> CohesiveState:
    """Return the virgin state for exactly one material parameter set."""
    if not isinstance(material, BilinearCohesiveMaterial):
        raise TypeError("material must be BilinearCohesiveMaterial")
    return CohesiveState(
        opening_m=0.0,
        kappa_m=0.0,
        damage=0.0,
        stored_energy_j_per_m2=0.0,
        dissipated_energy_j_per_m2=0.0,
        external_work_j_per_m2=0.0,
        update_count=0,
        material_signature=material.signature,
    )


def _validate_state(material: BilinearCohesiveMaterial, state: CohesiveState) -> None:
    if not isinstance(material, BilinearCohesiveMaterial):
        raise TypeError("material must be BilinearCohesiveMaterial")
    if not isinstance(state, CohesiveState):
        raise TypeError("state must be CohesiveState")
    if state.material_signature != material.signature:
        raise ValueError("state belongs to different cohesive material parameters")
    values = (
        state.opening_m,
        state.kappa_m,
        state.damage,
        state.stored_energy_j_per_m2,
        state.dissipated_energy_j_per_m2,
        state.external_work_j_per_m2,
    )
    if not all(math.isfinite(value) for value in values):
        raise ValueError("state fields must be finite")
    if state.opening_m < 0.0 or state.kappa_m < state.opening_m:
        raise ValueError("state requires 0 <= opening_m <= kappa_m")
    if not 0.0 <= state.damage <= 1.0 or state.update_count < 0:
        raise ValueError("state damage/update_count is invalid")
    expected_damage = _damage_at(material, state.kappa_m)
    expected_stored = _stored_energy(material, state.opening_m, expected_damage)
    expected_dissipated = _dissipation_at(material, state.kappa_m)
    scale = max(1.0, material.gc_j_per_m2, abs(state.external_work_j_per_m2))
    tolerance = 2.0e-11 * scale
    if abs(state.damage - expected_damage) > 2.0e-12:
        raise ValueError("state damage does not match its maximum-opening history")
    if abs(state.stored_energy_j_per_m2 - expected_stored) > tolerance:
        raise ValueError("state stored energy is inconsistent with its opening and damage")
    if abs(state.dissipated_energy_j_per_m2 - expected_dissipated) > tolerance:
        raise ValueError("state dissipation is inconsistent with its maximum-opening history")
    if abs(state.external_work_j_per_m2 - expected_stored - expected_dissipated) > tolerance:
        raise ValueError("state energy ledger does not satisfy work = stored + dissipated")


def update_state(
    material: BilinearCohesiveMaterial,
    state: CohesiveState,
    opening_m: float,
    *,
    mode: str = MODE_I_OPENING,
) -> CohesiveState:
    """Advance opening and energy state; damage changes only when kappa grows.

    Work is integrated exactly along the assumed quasi-static piecewise path:
    fixed-damage elastic return/reload to the prior maximum, then the linear
    envelope as the maximum opening grows. Therefore unloading changes stored
    energy and signed external work but never creates negative dissipation.
    """
    _check_mode(mode)
    _validate_state(material, state)
    new_opening = _finite_nonnegative("opening_m", opening_m)
    old_kappa = state.kappa_m
    new_kappa = max(old_kappa, new_opening)
    damage = _damage_at(material, new_kappa)

    if new_kappa == old_kappa:
        stiffness = material.k0_pa_per_m * (1.0 - state.damage)
        if stiffness == 0.0:
            work_increment = 0.0
        else:
            # The opening is bounded by the prior maximum while stiffness > 0.
            work_increment = (
                0.5
                * stiffness
                * (new_opening - state.opening_m)
                * (new_opening + state.opening_m)
            )
    else:
        # If this is a reload from below the old maximum, first recover/re-store
        # the elastic energy at fixed damage, then follow the envelope.
        stiffness = material.k0_pa_per_m * (1.0 - state.damage)
        if stiffness == 0.0:
            reload_work = 0.0
        else:
            reload_work = (
                0.5
                * stiffness
                * (old_kappa - state.opening_m)
                * (old_kappa + state.opening_m)
            )
        work_increment = reload_work + _envelope_work(material, old_kappa, new_kappa)

    dissipated = _dissipation_at(material, new_kappa)
    stored = _stored_energy(material, new_opening, damage)
    work = state.external_work_j_per_m2 + work_increment
    if not all(math.isfinite(value) for value in (work_increment, dissipated, stored, work)):
        raise FloatingPointError("cohesive energy update overflowed")
    next_state = CohesiveState(
        opening_m=new_opening,
        kappa_m=new_kappa,
        damage=damage,
        stored_energy_j_per_m2=stored,
        dissipated_energy_j_per_m2=dissipated,
        external_work_j_per_m2=work,
        update_count=state.update_count + 1,
        material_signature=material.signature,
    )
    _validate_state(material, next_state)
    return next_state


def traction_pa(
    material: BilinearCohesiveMaterial,
    state: CohesiveState,
    *,
    mode: str = MODE_I_OPENING,
) -> float:
    """Return the nonnegative scalar opening traction in Pa at ``state``."""
    _check_mode(mode)
    _validate_state(material, state)
    value = material.k0_pa_per_m * (1.0 - state.damage) * state.opening_m
    if not math.isfinite(value):
        raise FloatingPointError("cohesive traction overflowed")
    return value


def stored_energy_j_per_m2(
    material: BilinearCohesiveMaterial,
    state: CohesiveState,
) -> float:
    """Return recoverable elastic energy per unit area at the current opening."""
    _validate_state(material, state)
    return state.stored_energy_j_per_m2


def dissipated_energy_j_per_m2(
    material: BilinearCohesiveMaterial,
    state: CohesiveState,
) -> float:
    """Return cumulative, nonnegative damage dissipation per unit area."""
    _validate_state(material, state)
    return state.dissipated_energy_j_per_m2
