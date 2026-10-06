"""Raw, semantic-lifted consistent mass for the frozen mixed93 basis."""

from dataclasses import dataclass
import json
import math
from typing import Sequence


COMPONENTS = ("ux", "uy", "uz")
COMPONENT_INDEX = {name: i for i, name in enumerate(COMPONENTS)}
SOURCE_COUNT = 27
BRANCH_COUNT = 4
SCALAR_BASIS_COUNT = SOURCE_COUNT + BRANCH_COUNT
DOF_COUNT = 3 * SCALAR_BASIS_COUNT
VOLUME_POINT_COUNT = 512
MAX_CHUNK_ROWS = 256
MAX_WEIGHT_M3 = 1.0e-4
MAX_TOTAL_VOLUME_M3 = 1.0e-3
MAX_DENSITY_KG_M3 = 1.0e7


@dataclass(frozen=True)
class RawConsistentMass:
    """Immutable row-major little-endian float64 serialization of raw M."""

    semantic_ids: tuple[str, ...]
    shape: tuple[int, int]
    quadrature_count: int
    density_kg_m3: float
    volume_m3: float
    values_le_f64: bytes


def _semantic_lift(ids: Sequence[str]) -> tuple[list[list[int]], list[list[int]]]:
    if len(ids) != DOF_COUNT:
        raise ValueError(f"expected exactly {DOF_COUNT} semantic IDs")
    slots_by_component = [[], [], []]
    basis_by_component = [[], [], []]
    seen_sources = set()
    source_order = []
    for source_index in range(SOURCE_COUNT):
        owner = None
        for component_index, component in enumerate(COMPONENTS):
            slot = 3 * source_index + component_index
            raw = ids[slot]
            try:
                item = json.loads(raw)
            except (TypeError, ValueError) as exc:
                raise ValueError(f"semantic ID {slot} is not canonical JSON") from exc
            if (not isinstance(item, list) or len(item) != 4
                    or item[0] != "SOURCE_BASE" or not isinstance(item[1], str)
                    or not item[1] or item[1] == "global" or item[2] != component
                    or type(item[3]) is not int or item[3] != 0):
                raise ValueError(f"semantic ID {slot} violates source-major SOURCE_BASE layout")
            if json.dumps(item, ensure_ascii=True, separators=(",", ":")) != raw:
                raise ValueError(f"semantic ID {slot} is not in frozen canonical encoding")
            if owner is None:
                owner = item[1]
                if owner in seen_sources:
                    raise ValueError("source owners must be unique and contiguous")
                seen_sources.add(owner)
                source_order.append(owner)
            elif item[1] != owner:
                raise ValueError("each source's ux,uy,uz IDs must be contiguous")
            slots_by_component[component_index].append(slot)
            basis_by_component[component_index].append(source_index)
    if len(seen_sources) != SOURCE_COUNT or source_order != sorted(source_order):
        raise ValueError("source owners must match the 27-source frozen sorted order")
    for branch in range(1, BRANCH_COUNT + 1):
        for component_index, component in enumerate(COMPONENTS):
            slot = 3 * SOURCE_COUNT + 3 * (branch - 1) + component_index
            raw = ids[slot]
            try:
                item = json.loads(raw)
            except (TypeError, ValueError) as exc:
                raise ValueError(f"semantic ID {slot} is not canonical JSON") from exc
            if (item != ["GLOBAL_FRONT", "global", component, branch]
                    or type(item[3]) is not int):
                raise ValueError(f"semantic ID {slot} violates branch-major GLOBAL_FRONT layout")
            if json.dumps(item, ensure_ascii=True, separators=(",", ":")) != raw:
                raise ValueError(f"semantic ID {slot} is not in frozen canonical encoding")
            slots_by_component[component_index].append(slot)
            basis_by_component[component_index].append(SOURCE_COUNT + branch - 1)
    return slots_by_component, basis_by_component


def assemble_consistent_mass(ids: Sequence[str], basis_values: object,
                             weights: object, rho: float) -> RawConsistentMass:
    """Assemble raw M=ρ∫HᵀH dV; no symmetrization, lumping, or repair."""
    import numpy as np

    semantic_ids = tuple(str(identifier) for identifier in ids)
    slots_by_component, basis_by_component = _semantic_lift(semantic_ids)
    values = np.asarray(basis_values, dtype=np.float64)
    volume_weights = np.asarray(weights, dtype=np.float64)
    if values.shape != (VOLUME_POINT_COUNT, SCALAR_BASIS_COUNT):
        raise ValueError("volume_basis_values must have shape (512,31)")
    if volume_weights.shape != (VOLUME_POINT_COUNT,):
        raise ValueError("volume_weights_m3 must have shape (512,)")
    if isinstance(rho, bool):
        raise ValueError("rho must be a numeric SI density, not bool")
    density = float(rho)
    if not math.isfinite(density) or not 0.0 < density <= MAX_DENSITY_KG_M3:
        raise ValueError("rho must be finite, positive, and within the SI safety cap")
    if not np.isfinite(values).all():
        raise ValueError("volume basis values must be finite")
    if (not np.isfinite(volume_weights).all() or np.any(volume_weights <= 0.0)
            or np.any(volume_weights > MAX_WEIGHT_M3)):
        raise ValueError("volume weights must be finite positive m³ within the SI safety cap")
    volume = float(np.sum(volume_weights, dtype=np.float64))
    if not math.isfinite(volume) or volume > MAX_TOTAL_VOLUME_M3:
        raise ValueError("quadrature volume exceeds the fixed SI safety cap")

    scalar_mass = np.zeros((SCALAR_BASIS_COUNT, SCALAR_BASIS_COUNT), dtype=np.float64)
    for start in range(0, VOLUME_POINT_COUNT, MAX_CHUNK_ROWS):
        stop = min(start + MAX_CHUNK_ROWS, VOLUME_POINT_COUNT)
        block = values[start:stop]
        weighted_block = block * (density * volume_weights[start:stop])[:, None]
        scalar_mass += block.T @ weighted_block
    if not np.isfinite(scalar_mass).all():
        raise ValueError("consistent scalar mass overflowed or became nonfinite")

    raw_mass = np.zeros((DOF_COUNT, DOF_COUNT), dtype=np.float64)
    for component_index in range(len(COMPONENTS)):
        slots = slots_by_component[component_index]
        basis = basis_by_component[component_index]
        raw_mass[np.ix_(slots, slots)] = scalar_mass[np.ix_(basis, basis)]
    payload = np.asarray(raw_mass, dtype="<f8").tobytes(order="C")
    return RawConsistentMass(semantic_ids, (DOF_COUNT, DOF_COUNT),
                             VOLUME_POINT_COUNT, density, volume, payload)
