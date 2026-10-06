"""Bounded raw elastic operator for a 93-DOF global-front semantic layout.
It uses per-source base shapes and four global tip branches; the 405-DOF
source-owned layout is diagnostic only. No provider, solve, truncation, or
fracture behavior is implemented here."""
from __future__ import annotations
from dataclasses import dataclass
import json
import math
from numbers import Real
import struct

MAX_SOURCE_COUNT, MAX_ROWS_PER_CHUNK, MAX_CHUNKS = 27, 256, 256
MAX_GLOBAL_SCALAR_BASIS, MAX_GLOBAL_VECTOR_DOFS = 31, 93
MAX_DIAGNOSTIC_SCALAR_BASIS, MAX_DIAGNOSTIC_VECTOR_DOFS = 135, 405
GLOBAL_FRONT_LAYOUT = "global_front"
PER_SOURCE_DIAGNOSTIC_LAYOUT = "per_source_front_diagnostic"
COMPONENTS, BRANCHES = ("ux", "uy", "uz"), (1, 2, 3, 4)

class MixedOperatorInputError(ValueError):
    """Malformed, non-finite, ambiguous, or semantically misordered input."""

@dataclass(frozen=True, slots=True)
class SemanticContribution:
    semantic_id: str
    scope: str
    owner_id: str
    component: str
    basis_index: int
    shape_value: float
    shape_gradient_world_per_m: tuple[float, float, float]
    B_column_per_m: tuple[float, float, float, float, float, float]

@dataclass(frozen=True, slots=True)
class SemanticShapeRow:
    layout: str
    source_ids: tuple[str, ...]
    coefficient_ids: tuple[str, ...]
    contributions: tuple[SemanticContribution, ...]
    status: str = "OK"

@dataclass(frozen=True, slots=True)
class VolumePoint:
    weight_m3: float
    shape_row: object

@dataclass(frozen=True, slots=True)
class SurfaceTractionPoint:
    weight_m2: float
    traction_pa: tuple[float, float, float]
    shape_row: object

@dataclass(frozen=True, slots=True)
class RawStiffness:
    layout: str
    source_ids: tuple[str, ...]
    coefficient_ids: tuple[str, ...]
    dof_count: int
    quadrature_count: int
    constitutive_matrix_pa: tuple[tuple[float, ...], ...]
    values_le_f64: bytes
    matrix_shape: tuple[int, int]
    units: str = "N/m"

@dataclass(frozen=True, slots=True)
class RawForce:
    layout: str
    source_ids: tuple[str, ...]
    coefficient_ids: tuple[str, ...]
    dof_count: int
    quadrature_count: int
    values_le_f64: bytes
    vector_shape: tuple[int, ...]
    units: str = "N"

@dataclass(frozen=True, slots=True)
class CoefficientModes:
    layout: str
    mode_ids: tuple[str, ...]
    coefficient_ids: tuple[str, ...]
    mode_count: int
    dof_count: int
    values_le_f64: bytes
    matrix_shape: tuple[int, int]

def _real(name: str, value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise MixedOperatorInputError(f"{name} must be a finite real")
    try: result = float(value)
    except (OverflowError, ValueError) as exc: raise MixedOperatorInputError(f"{name} must be finite") from exc
    if not math.isfinite(result): raise MixedOperatorInputError(f"{name} must be finite")
    return result

def _vec(name: str, value: object, n: int) -> tuple[float, ...]:
    if type(value) not in (tuple, list) or len(value) != n:
        raise MixedOperatorInputError(f"{name} must contain exactly {n} values")
    return tuple(_real(f"{name}[{i}]", v) for i, v in enumerate(value))

def _sources(value: object) -> tuple[str, ...]:
    if type(value) is not tuple or not 1 <= len(value) <= MAX_SOURCE_COUNT:
        raise MixedOperatorInputError("source_ids must be a 1..27 tuple")
    if any(type(s) is not str or not s or s.strip() != s for s in value) or len(set(value)) != len(value):
        raise MixedOperatorInputError("source IDs must be unique canonical strings")
    return tuple(sorted(value))

def _sid(scope: str, owner: str, component: str, index: int) -> str:
    return json.dumps([scope, owner, component, index], ensure_ascii=True, separators=(",", ":"))

def _specs(sources: tuple[str, ...], layout: str):
    if layout == GLOBAL_FRONT_LAYOUT:
        return ([("SOURCE_BASE", s, c, 0) for s in sources for c in COMPONENTS]
                + [("GLOBAL_FRONT", "global", c, b) for b in BRANCHES for c in COMPONENTS])
    if layout == PER_SOURCE_DIAGNOSTIC_LAYOUT:
        return ([("SOURCE_BASE", s, c, 0) for s in sources for c in COMPONENTS]
                + [("SOURCE_FRONT_DIAGNOSTIC", s, c, b) for s in sources for b in BRANCHES for c in COMPONENTS])
    raise MixedOperatorInputError("unknown semantic layout")

def semantic_coefficient_ids(source_ids: object, layout: str = GLOBAL_FRONT_LAYOUT) -> tuple[str, ...]:
    return tuple(_sid(*spec) for spec in _specs(_sources(source_ids), layout))

def vector_dof_count(source_count: int, layout: str = GLOBAL_FRONT_LAYOUT) -> int:
    if type(source_count) is not int or not 1 <= source_count <= MAX_SOURCE_COUNT:
        raise MixedOperatorInputError("source count outside 1..27")
    count = len(_specs(tuple(f"s{i}" for i in range(source_count)), layout))
    cap = MAX_GLOBAL_VECTOR_DOFS if layout == GLOBAL_FRONT_LAYOUT else MAX_DIAGNOSTIC_VECTOR_DOFS
    if count > cap: raise MixedOperatorInputError("semantic layout exceeds its explicit DOF cap")
    return count

def scalar_basis_count(source_count: int, layout: str = GLOBAL_FRONT_LAYOUT) -> int:
    return vector_dof_count(source_count, layout) // 3

def _b(component: str, g: tuple[float, ...]) -> tuple[float, ...]:
    gx, gy, gz = g
    if component == "ux": return (gx, 0.0, 0.0, gy, 0.0, gz)
    if component == "uy": return (0.0, gy, 0.0, gx, gz, 0.0)
    return (0.0, 0.0, gz, 0.0, gy, gx)

def _old_id(source: str, component: str, kind: str, index: int) -> str:
    return json.dumps([source, component, kind, index], ensure_ascii=True, separators=(",", ":"))

def adapt_mixed_shape_row(source_row: object, source_ids: object,
                          layout: str = GLOBAL_FRONT_LAYOUT) -> SemanticShapeRow:
    """Map the 405-row helper to 93 global or explicit redundant 405 semantics."""
    sources = _sources(source_ids)
    specs = _specs(sources, layout)
    status = getattr(source_row, "status", None)
    if getattr(status, "value", status) != "OK" or getattr(source_row, "source_ids", None) != sources:
        raise MixedOperatorInputError("mixed helper row is non-OK or source IDs are misordered")
    rows, ids = getattr(source_row, "contributions", None), getattr(source_row, "coefficient_ids", None)
    fv, fg = getattr(source_row, "branch_values", None), getattr(source_row, "branch_gradients_world_per_m", None)
    if type(fv) is not tuple or len(fv) != 4 or type(fg) is not tuple or len(fg) != 4:
        raise MixedOperatorInputError("mixed helper row must expose four branch values/gradients")
    fv = tuple(_real(f"branch value {i}", v) for i, v in enumerate(fv))
    fg = tuple(_vec(f"branch gradient {i}", v, 3) for i, v in enumerate(fg))
    old_specs = [(s, c, k, b) for s in sources for k, b in (("base", 0), *(('front_branch', i) for i in BRANCHES)) for c in COMPONENTS]
    if type(rows) is not tuple or type(ids) is not tuple or len(rows) != len(old_specs) or len(ids) != len(old_specs):
        raise MixedOperatorInputError("mixed helper row must be a complete five-function row")
    data = {}
    for j, (item, (s, c, kind, index)) in enumerate(zip(rows, old_specs)):
        if (ids[j] != _old_id(s, c, kind, index) or getattr(item, "semantic_id", None) != ids[j]
                or getattr(item, "source_id", None) != s or getattr(item, "component", None) != c
                or getattr(item, "basis_kind", None) != kind or type(getattr(item, "basis_index", None)) is not int
                or item.basis_index != index):
            raise MixedOperatorInputError(f"mixed helper ID/order mismatch at row {j}")
        v = _real(f"helper row {j} value", getattr(item, "shape_value", None))
        g = _vec(f"helper row {j} gradient", getattr(item, "shape_gradient_world_per_m", None), 3)
        bc = _vec(f"helper row {j} B", getattr(item, "B_column_per_m", None), 6)
        if bc != _b(c, g): raise MixedOperatorInputError(f"helper row {j} B/gradient mismatch")
        if kind == "front_branch":
            n, ng = data[(s, c, "base", 0)]
            expected_g = tuple(fv[index-1] * ng[a] + n * fg[index-1][a] for a in range(3))
            if v != n * fv[index-1] or g != expected_g:
                raise MixedOperatorInputError(f"helper row {j} does not match global product rule")
        data[(s, c, kind, index)] = (v, g)
    if ids != tuple(_old_id(*spec) for spec in old_specs):
        raise MixedOperatorInputError("mixed helper coefficient_ids are not canonical")
    out = []
    for scope, owner, component, index in specs:
        if scope == "SOURCE_BASE": v, g = data[(owner, component, "base", 0)]
        elif scope == "SOURCE_FRONT_DIAGNOSTIC": v, g = data[(owner, component, "front_branch", index)]
        else:
            v = _real(f"branch value {index}", fv[index - 1])
            g = _vec(f"branch gradient {index}", fg[index - 1], 3)
        out.append(SemanticContribution(_sid(scope, owner, component, index), scope, owner,
                                        component, index, v, g, _b(component, g)))
    out_ids = tuple(item.semantic_id for item in out)
    expected = tuple(_sid(*spec) for spec in specs)
    if out_ids != expected: raise MixedOperatorInputError("adapted semantic row order mismatch")
    if layout == GLOBAL_FRONT_LAYOUT:
        for s in sources:
            base = [data[(s, c, "base", 0)] for c in COMPONENTS]
            if any(item != base[0] for item in base[1:]):
                raise MixedOperatorInputError("base scalar rows disagree across vector components")
    return SemanticShapeRow(layout, sources, expected, tuple(out))

def isotropic_constitutive_matrix_pa(young_modulus_pa: object, poisson_ratio: object) -> tuple[tuple[float, ...], ...]:
    """Full isotropic C for [exx,eyy,ezz,gamma_xy,gamma_yz,gamma_xz]."""
    e, nu = _real("E_pa", young_modulus_pa), _real("nu", poisson_ratio)
    if e <= 0.0 or not -1.0 < nu < 0.5: raise MixedOperatorInputError("require E>0 Pa and -1<nu<0.5")
    mu, lam = e / (2.0 * (1.0 + nu)), e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    c = [[0.0] * 6 for _ in range(6)]
    for i in range(3):
        for j in range(3): c[i][j] = lam + (2.0 * mu if i == j else 0.0)
    for i in range(3, 6): c[i][i] = mu
    if any(not math.isfinite(v) for row in c for v in row): raise MixedOperatorInputError("C overflowed")
    return tuple(tuple(row) for row in c)

def _c_matrix(value: object) -> tuple[tuple[float, ...], ...]:
    if type(value) not in (tuple, list) or len(value) != 6: raise MixedOperatorInputError("C must be 6x6 Pa")
    c = tuple(_vec(f"C[{i}]", row, 6) for i, row in enumerate(value))
    if any(c[i][j] != c[j][i] for i in range(6) for j in range(6)):
        raise MixedOperatorInputError("C must be exactly symmetric")
    lower = [[0.0] * 6 for _ in range(6)]
    for i in range(6):
        for j in range(i + 1):
            residual = c[i][j] - sum(lower[i][k] * lower[j][k] for k in range(j))
            if i == j:
                if residual <= 0.0: raise MixedOperatorInputError("C must be positive definite")
                lower[i][j] = math.sqrt(residual)
            else: lower[i][j] = residual / lower[j][j]
            if not math.isfinite(lower[i][j]): raise MixedOperatorInputError("C Cholesky check overflowed")
    return c

def _validated(row: object, sources: tuple[str, ...], layout: str):
    if getattr(row, "status", None) != "OK" or getattr(row, "layout", None) != layout:
        raise MixedOperatorInputError("shape row status/layout mismatch")
    if getattr(row, "source_ids", None) != sources: raise MixedOperatorInputError("shape source order mismatch")
    specs, ids = _specs(sources, layout), semantic_coefficient_ids(sources, layout)
    items = getattr(row, "contributions", None)
    if getattr(row, "coefficient_ids", None) != ids or type(items) is not tuple or len(items) != len(ids):
        raise MixedOperatorInputError("shape row semantic IDs/contributions are incomplete")
    values, columns, scalar_fields = [], [], {}
    for i, (item, (scope, owner, component, index)) in enumerate(zip(items, specs)):
        if (getattr(item, "semantic_id", None) != ids[i] or getattr(item, "scope", None) != scope
                or getattr(item, "owner_id", None) != owner or getattr(item, "component", None) != component
                or type(getattr(item, "basis_index", None)) is not int or item.basis_index != index):
            raise MixedOperatorInputError(f"shape contribution {i} ID/order mismatch")
        v = _real(f"shape value {i}", getattr(item, "shape_value", None))
        g = _vec(f"shape gradient {i}", getattr(item, "shape_gradient_world_per_m", None), 3)
        b = _vec(f"B column {i}", getattr(item, "B_column_per_m", None), 6)
        if b != _b(component, g): raise MixedOperatorInputError(f"B/gradient mismatch in column {i}")
        key = (scope, owner, index)
        if key in scalar_fields and scalar_fields[key] != (v, g):
            raise MixedOperatorInputError(f"scalar field differs across components for {key}")
        scalar_fields[key] = (v, g)
        values.append(v); columns.append(b)
    return ids, tuple(values), tuple(columns)

def _chunks(chunks: object, row_type: type, name: str):
    if type(chunks) is not tuple or not 1 <= len(chunks) <= MAX_CHUNKS:
        raise MixedOperatorInputError(f"{name} must contain 1..256 chunks")
    if any(type(chunk) is not tuple or not 1 <= len(chunk) <= MAX_ROWS_PER_CHUNK for chunk in chunks):
        raise MixedOperatorInputError(f"{name} chunks must each contain 1..256 rows")
    if any(type(point) is not row_type for chunk in chunks for point in chunk):
        raise MixedOperatorInputError(f"{name} contains an unexpected row type")
    return sum(map(len, chunks))

def _pack(values: list[float]) -> bytes:
    if any(not math.isfinite(v) for v in values): raise FloatingPointError("raw output became non-finite")
    return struct.pack(f"<{len(values)}d", *values)

def assemble_raw_stiffness(source_ids: object, volume_chunks: object, constitutive_matrix_pa: object,
                           layout: str = GLOBAL_FRONT_LAYOUT) -> RawStiffness:
    """Return raw K=sum(w_m3*B.T*C*B); a 256-row chunk never exceeds its cap."""
    sources = _sources(source_ids)
    q_count = _chunks(volume_chunks, VolumePoint, "volume_chunks")
    c, ids = _c_matrix(constitutive_matrix_pa), semantic_coefficient_ids(sources, layout)
    n, k = len(ids), [0.0] * (len(ids) * len(ids))
    for chunk in volume_chunks:
        for point in chunk:
            w = _real("weight_m3", point.weight_m3)
            if w <= 0.0: raise MixedOperatorInputError("volume weight must be positive SI m^3")
            row_ids, _, bcols = _validated(point.shape_row, sources, layout)
            if row_ids != ids: raise MixedOperatorInputError("volume coefficient IDs changed")
            cb = tuple(tuple(sum(c[i][j] * bcols[d][j] for j in range(6))
                             for i in range(6)) for d in range(n))
            for i, bi in enumerate(bcols):
                base = i * n
                for j in range(n): k[base + j] += w * sum(bi[a] * cb[j][a] for a in range(6))
    return RawStiffness(layout, sources, ids, n, q_count, c, _pack(k), (n, n))

def assemble_raw_force(source_ids: object, surface_chunks: object,
                       layout: str = GLOBAL_FRONT_LAYOUT) -> RawForce:
    """Return raw f=sum(w_m2*Nvec.T*traction_pa) without load projection."""
    sources = _sources(source_ids)
    q_count = _chunks(surface_chunks, SurfaceTractionPoint, "surface_chunks")
    ids, specs = semantic_coefficient_ids(sources, layout), _specs(sources, layout)
    n, f = len(ids), [0.0] * len(ids)
    ci = {component: i for i, component in enumerate(COMPONENTS)}
    for chunk in surface_chunks:
        for point in chunk:
            w = _real("weight_m2", point.weight_m2)
            if w <= 0.0: raise MixedOperatorInputError("surface weight must be positive SI m^2")
            t = _vec("traction_pa", point.traction_pa, 3)
            row_ids, values, _ = _validated(point.shape_row, sources, layout)
            if row_ids != ids: raise MixedOperatorInputError("surface coefficient IDs changed")
            for i, value in enumerate(values): f[i] += w * value * t[ci[specs[i][2]]]
    return RawForce(layout, sources, ids, n, q_count, _pack(f), (n,))

def _positions(source_ids: object, positions_m: object) -> tuple[tuple[str, ...], dict[str, tuple[float, ...]]]:
    if type(source_ids) is not tuple or type(positions_m) not in (tuple, list) or len(source_ids) != len(positions_m):
        raise MixedOperatorInputError("positions must pair with source IDs")
    sources = _sources(source_ids)
    raw = {s: _vec(f"positions_m[{i}]", positions_m[i], 3) for i, s in enumerate(source_ids)}
    return sources, {s: raw[s] for s in sources}

def _modes(sources: tuple[str, ...], layout: str, names: tuple[str, ...], rows: list[list[float]]) -> CoefficientModes:
    ids = semantic_coefficient_ids(sources, layout)
    if len(rows) != len(names) or any(len(row) != len(ids) for row in rows):
        raise MixedOperatorInputError("mode coefficient matrix shape mismatch")
    return CoefficientModes(layout, names, ids, len(names), len(ids),
                            _pack([v for row in rows for v in row]), (len(names), len(ids)))

def rigid_body_mode_coefficients(source_ids: object, source_positions_world_m: object,
                                 layout: str = GLOBAL_FRONT_LAYOUT) -> CoefficientModes:
    """Six physical candidates; exact nullness needs base partition/linear reproduction."""
    sources, p = _positions(source_ids, source_positions_world_m)
    ids, slots = semantic_coefficient_ids(sources, layout), {}
    slots = {identifier: i for i, identifier in enumerate(ids)}
    rows = [[0.0] * len(ids) for _ in range(6)]
    center = tuple(sum(p[s][j] for s in sources) / len(sources) for j in range(3))
    for s in sources:
        x = tuple(p[s][j] - center[j] for j in range(3))
        for c, component in enumerate(COMPONENTS):
            slot = slots[_sid("SOURCE_BASE", s, component, 0)]
            rows[c][slot] = 1.0
            for axis in range(3):
                w = tuple(1.0 if j == axis else 0.0 for j in range(3))
                cross = (w[1]*x[2]-w[2]*x[1], w[2]*x[0]-w[0]*x[2], w[0]*x[1]-w[1]*x[0])
                rows[3 + axis][slot] = cross[c]
    names = ("translation_x", "translation_y", "translation_z", "rotation_x", "rotation_y", "rotation_z")
    return _modes(sources, layout, names, rows)

def branch_dependency_mode_coefficients(source_ids: object,
                                        source_local_positions_from_tip_m: object) -> CoefficientModes:
    """Six extra zero-field candidates for the redundant per-source 405 layout.

    Uses x*F3-y*F4+y*F1=0 and x*F4+y*F3-y*F2=0 under linear reproduction.
    """
    sources, p = _positions(source_ids, source_local_positions_from_tip_m)
    layout, ids = PER_SOURCE_DIAGNOSTIC_LAYOUT, semantic_coefficient_ids(sources, PER_SOURCE_DIAGNOSTIC_LAYOUT)
    slots, names, rows = {v: i for i, v in enumerate(ids)}, [], []
    for identity in (1, 2):
        for component in COMPONENTS:
            row = [0.0] * len(ids)
            for s in sources:
                x, y = p[s][:2]
                terms = ((1, y), (3, x), (4, -y)) if identity == 1 else ((2, -y), (3, y), (4, x))
                for b, value in terms: row[slots[_sid("SOURCE_FRONT_DIAGNOSTIC", s, component, b)]] = value
            rows.append(row); names.append(f"branch_identity_{identity}_{component}")
    return _modes(sources, layout, tuple(names), rows)
