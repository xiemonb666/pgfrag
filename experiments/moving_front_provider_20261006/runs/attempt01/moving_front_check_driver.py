"""One frozen, bounded diagnostic for the two moving-front source states.

SOURCE ONLY until root marks the matching check spec READY. This driver makes
two provider batches (25 point/trace identities each), no separate MLS calls,
no operator, no stiffness assembly, and no solver call. It has no CLI inputs.
"""

from __future__ import annotations

import hashlib
import json
import math
import time
from pathlib import Path
from typing import Any, Mapping, Sequence

_STARTED_MONOTONIC = time.monotonic()

HERE = Path(__file__).resolve().parent
SPEC_PATH = HERE / "moving_front_check_spec.json"
if __name__ == "__main__":
    _preflight_spec = json.loads(SPEC_PATH.read_text(encoding="utf-8"))
    if _preflight_spec.get("status") != "READY":
        raise SystemExit("execution blocked until root sets moving_front_check_spec.json status to READY")

import numpy as np  # noqa: E402
import moving_front_provider as mf  # noqa: E402


REPO_ROOT = HERE.parents[2]
OUTPUT_PATH = HERE / "moving_front_check_result.json"
RUN_MARKER_PATH = HERE / "moving_front_check.run-once"
PROVIDER_PATH = HERE / "moving_front_provider.py"
CASES_PATH = HERE / "moving_front_cases.json"

PROVIDER_SHA256 = "b3e88622220fae0ed63961b943ba710e5689b470d7df00361c297934c42dfac6"
CASES_SHA256 = "2766245b25bf93eefcf6836befe81536917888936b4543666c86ebe9d6b72100"
SOURCE_COORDINATE_RECORDS_SHA256 = "a8687a23353b6465aa0cc3fcda49a57c240297c38c141a017a167f3c526b23bb"

SOURCE_COUNT = 100
PROVIDER_POINTS_PER_STATE = 25
PROVIDER_POINTS_TOTAL = 50
SOURCE_RECLASSIFICATIONS_TOTAL = 200
TANGENCY_CONTROLS_TOTAL = 4
RECORDS_TOTAL = 254
RECORDS_HARD_CAP = 256
MAX_WALL_SECONDS = 30.0
MAX_DIAGNOSTIC_BYTES = 2 * 1024 * 1024
MAX_RESULT_JSON_BYTES = 256 * 1024
TRANSIENT_RESERVE_BYTES = 512 * 1024

AFFINE_VALUE_TOL_M = 1.0e-10
AFFINE_GRAD_TOL = 1.0e-9
CORE_VALUE_TOL_M = 1.0e-12
CORE_GRAD_TOL = 1.0e-9
TRACE_BASIS_VALUE_TOL = 1.0e-12
TRACE_BASIS_GRAD_TOL_PER_M = 1.0e-9
CORE_RADII_CAPS_M = (0.00025, 0.0005)
CORE_RADII_FRACTIONS = (0.25, 0.5)
FACE_X_OFFSETS_M = (-0.008, -0.004)
COLLAR_OFFSET_M = 0.002
COLLAR_Y_M = 0.0005
AHEAD_X_OFFSET_M = 0.004
AHEAD_Y_EPSILON_M = 1.0e-8

AFFINE_INTERCEPT_M = 0.2
AFFINE_GRADIENT = (0.25, -0.5, 0.75)


class CheckFailure(RuntimeError):
    """A frozen check precondition or acceptance assertion failed."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise CheckFailure(message)


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _sha256_file(path: Path) -> str:
    return _sha256_bytes(path.read_bytes())


def _deadline(stage: str) -> None:
    elapsed = time.monotonic() - _STARTED_MONOTONIC
    if elapsed > MAX_WALL_SECONDS:
        raise CheckFailure(f"30-second budget exceeded at {stage}: {elapsed:.6f} s")


def _expected_spec() -> dict[str, Any]:
    return {
        "schema": "moving-front-provider-check-spec.v1",
        "budgets": {
            "states": 2,
            "sources_per_state": SOURCE_COUNT,
            "maximum_scalar_dofs": 375,
            "provider_points_per_state": PROVIDER_POINTS_PER_STATE,
            "provider_points_total": PROVIDER_POINTS_TOTAL,
            "source_reclassifications_total": SOURCE_RECLASSIFICATIONS_TOTAL,
            "pure_topology_tangency_controls_total": TANGENCY_CONTROLS_TOTAL,
            "planned_records_total": RECORDS_TOTAL,
            "hard_records_cap": RECORDS_HARD_CAP,
            "maximum_wall_seconds": MAX_WALL_SECONDS,
            "maximum_diagnostic_owned_bytes": MAX_DIAGNOSTIC_BYTES,
            "maximum_result_json_bytes": MAX_RESULT_JSON_BYTES,
            "transient_workspace_reserve_bytes": TRANSIENT_RESERVE_BYTES,
            "mls_batches_total": 2,
        },
        "tolerances": {
            "affine_displacement_abs_m": AFFINE_VALUE_TOL_M,
            "affine_gradient_abs_dimensionless": AFFINE_GRAD_TOL,
            "core_displacement_abs_m": CORE_VALUE_TOL_M,
            "core_gradient_abs_dimensionless": CORE_GRAD_TOL,
            "trace_basis_value_abs_dimensionless": TRACE_BASIS_VALUE_TOL,
            "trace_basis_gradient_abs_per_m": TRACE_BASIS_GRAD_TOL_PER_M,
        },
        "core_radius_rule": {
            "definition": "min distance from current finite F to every currently non-F1 source center minus h",
            "radius_caps_m": list(CORE_RADII_CAPS_M),
            "radius_fractions_of_measured_positive_core": list(CORE_RADII_FRACTIONS),
            "formula": [
                "r1=min(0.00025 m, 0.25*r_core)",
                "r2=min(0.00050 m, 0.50*r_core)",
            ],
            "nonpositive_core_action": "BLOCK; no automatic radius search or algorithm adjustment",
        },
        "query_layout": {
            "core_points": 8,
            "exact_face_traces": 4,
            "collar_observations": 2,
            "ahead_finite_sample_points": 3,
            "expected_unknowns": 8,
            "face_x_offsets_m": list(FACE_X_OFFSETS_M),
            "collar_x_offset_m": COLLAR_OFFSET_M,
            "collar_abs_y_m": COLLAR_Y_M,
            "ahead_x_offset_m": AHEAD_X_OFFSET_M,
            "ahead_abs_y_epsilon_m": AHEAD_Y_EPSILON_M,
        },
        "manufactured_affine_fixture": {
            "intercept_m": AFFINE_INTERCEPT_M,
            "gradient_dimensionless": list(AFFINE_GRADIENT),
        },
        "source_coordinate_record_digest": {
            "source_count": SOURCE_COUNT,
            "coordinate_records_sha256": SOURCE_COORDINATE_RECORDS_SHA256,
            "role_fields_are_not_used": True,
        },
        "continuity_scope": "The three y=-1e-8,0,+1e-8 m ahead rows are finite-sample diagnostics only, not a global continuity proof.",
        "collar_scope": "Record errors only; no acceptance threshold and no full-field exactness claim.",
        "coefficient_scope": "Fresh analytic manufactured coefficients only; no old coefficients, transfer, remap, or mechanical solve.",
        "output_path": "moving_front_check_result.json",
        "run_marker_path": "moving_front_check.run-once",
        "output_policy": "Create once with exclusive mode; do not overwrite an existing result.",
        "memory_scope": "Diagnostic-owned NumPy arrays plus the fixed transient reserve and serialized result; interpreter/import baseline excluded.",
    }


def _load_frozen_spec() -> tuple[dict[str, Any], str, str]:
    _require(SPEC_PATH.is_file(), "frozen check spec is missing")
    raw = SPEC_PATH.read_bytes()
    spec = json.loads(raw.decode("utf-8"))
    _require(spec.get("status") == "READY", "execution blocked until root sets the frozen spec status to READY")
    expected = _expected_spec()
    for key, value in expected.items():
        _require(spec.get(key) == value, f"check spec changed at {key}")
    _require(SOURCE_COORDINATE_RECORDS_SHA256 == mf.COORDINATE_RECORDS_SHA256, "coordinate-record digest differs between driver and provider pins")

    driver_hash = _sha256_file(Path(__file__).resolve())
    source_files = spec.get("source_files")
    _require(isinstance(source_files, dict), "check spec lacks frozen source file pins")
    expected_paths = {
        "driver": "pgfrag/experiments/moving_front_provider_20261006/moving_front_check_driver.py",
        "provider": "pgfrag/experiments/moving_front_provider_20261006/moving_front_provider.py",
        "case_pair": "pgfrag/experiments/moving_front_provider_20261006/moving_front_cases.json",
        "safe_v2": mf.SAFE_CLASSIFIER_REL,
        "fixed_source_case": mf.SOURCE_CASE_REL,
        "coordinate_manifest": mf.COORDINATE_MANIFEST_REL,
        "registry": mf.REGISTRY_REL,
        "mls_core": mf.MLS_CORE_REL,
        "mls_shape": mf.MLS_SHAPE_REL,
        "material_case": mf.MATERIAL_CASE_REL,
    }
    _require(set(source_files) == set(expected_paths), "check spec source pin names changed")
    for name, relative_path in expected_paths.items():
        entry = source_files.get(name)
        _require(isinstance(entry, dict) and set(entry) == {"path", "sha256"}, f"malformed frozen pin: {name}")
        _require(entry["path"] == relative_path, f"frozen input path changed: {name}")
    _require(source_files.get("driver", {}).get("sha256") == driver_hash, "driver source hash differs from frozen spec")
    _require(source_files.get("provider", {}).get("sha256") == PROVIDER_SHA256, "provider pin in spec differs from driver")
    _require(source_files.get("case_pair", {}).get("sha256") == CASES_SHA256, "case-pair pin in spec differs from driver")

    for name, path in (
        ("provider", PROVIDER_PATH),
        ("case_pair", CASES_PATH),
        ("safe_v2", REPO_ROOT / mf.SAFE_CLASSIFIER_REL),
        ("fixed_source_case", REPO_ROOT / mf.SOURCE_CASE_REL),
        ("coordinate_manifest", REPO_ROOT / mf.COORDINATE_MANIFEST_REL),
        ("registry", REPO_ROOT / mf.REGISTRY_REL),
        ("mls_core", REPO_ROOT / mf.MLS_CORE_REL),
        ("mls_shape", REPO_ROOT / mf.MLS_SHAPE_REL),
        ("material_case", REPO_ROOT / mf.MATERIAL_CASE_REL),
    ):
        _require(path.is_file(), f"frozen input is missing: {name}")
        expected_hash = source_files.get(name, {}).get("sha256")
        _require(type(expected_hash) is str and _sha256_file(path) == expected_hash, f"frozen input hash changed: {name}")

    _require(_sha256_file(PROVIDER_PATH) == PROVIDER_SHA256, "provider source changed after review")
    _require(_sha256_file(CASES_PATH) == CASES_SHA256, "case-pair input changed after review")
    _require(Path(mf.__file__).resolve() == PROVIDER_PATH.resolve(), "provider module resolved to an unexpected path")
    _deadline("frozen input validation")
    return spec, _sha256_bytes(raw), driver_hash


def _distance_to_front(point: Sequence[float], state: mf.FrontState) -> float:
    x, y, z = (float(value) for value in point)
    z_lo, z_hi = state.geometry.crack_z_bounds_m
    z_gap = max(z_lo - z, 0.0, z - z_hi)
    return math.sqrt(
        (x - state.geometry.front_x_m) ** 2
        + (y - state.geometry.y_plane_m) ** 2
        + z_gap**2
    )


def _current_core_radius(state: mf.FrontState) -> float:
    non_f1 = [source for source in state.sources if source.basis_role != "FRONT_F1_SHIFTED"]
    _require(non_f1, "core radius is undefined because there are no non-F1 sources")
    return min(_distance_to_front(source.position_m, state) for source in non_f1) - state.geometry.support_radius_m


def _make_queries(state: mf.FrontState, r_core_m: float) -> tuple[list[tuple[float, float, float]], list[mf.TraceSide | None], list[dict[str, Any]], tuple[int, ...]]:
    if not math.isfinite(r_core_m) or r_core_m <= 0.0:
        raise CheckFailure("BLOCKED: current state has no positive measured core radius")
    radii = tuple(
        min(cap, fraction * r_core_m)
        for cap, fraction in zip(CORE_RADII_CAPS_M, CORE_RADII_FRACTIONS)
    )
    if any(not math.isfinite(value) or value <= 0.0 for value in radii) or radii[1] <= radii[0]:
        raise CheckFailure("BLOCKED: frozen core-radius formula did not yield two distinct positive radii")

    points: list[tuple[float, float, float]] = []
    sides: list[mf.TraceSide | None] = []
    specs: list[dict[str, Any]] = []
    core_indices: list[int] = []

    def add(label: str, point: Sequence[float], side: mf.TraceSide | None, expected: str, reason: str | None = None) -> int:
        index = len(points)
        point3 = tuple(float(value) for value in point)
        _require(len(point3) == 3 and all(math.isfinite(value) for value in point3), f"invalid frozen query {label}")
        points.append(point3)
        sides.append(side)
        specs.append(
            {
                "label": label,
                "point_m": list(point3),
                "trace_side": None if side is None else side.value,
                "expected": expected,
                "expected_reason": reason,
            }
        )
        return index

    for radius_index, radius in enumerate(radii):
        delta = radius / math.sqrt(2.0)
        for sx, sy, label in (
            (-1.0, -1.0, "mm"),
            (-1.0, 1.0, "mp"),
            (1.0, -1.0, "pm"),
            (1.0, 1.0, "pp"),
        ):
            core_indices.append(
                add(
                    f"core-r{radius_index}-{label}",
                    (
                        state.front_x_m + sx * delta,
                        state.geometry.y_plane_m + sy * delta,
                        0.0,
                    ),
                    None,
                    "OK",
                )
            )

    trace_indices: list[int] = []
    for face_offset in FACE_X_OFFSETS_M:
        x = state.front_x_m + face_offset
        for side in (mf.TraceSide.PLUS, mf.TraceSide.MINUS):
            trace_indices.append(add(f"trace-{face_offset:g}-{side.value.lower()}", (x, state.geometry.y_plane_m, 0.0), side, "OK"))

    collar_indices = []
    for y_sign in (-1.0, 1.0):
        collar_indices.append(
            add(
                f"collar-{int(y_sign):+d}",
                (state.front_x_m - COLLAR_OFFSET_M, state.geometry.y_plane_m + y_sign * COLLAR_Y_M, 0.0),
                None,
                "OK",
            )
        )

    ahead_indices = []
    for y in (-AHEAD_Y_EPSILON_M, 0.0, AHEAD_Y_EPSILON_M):
        ahead_indices.append(
            add(
                f"ahead-y{y:+.0e}",
                (state.front_x_m + AHEAD_X_OFFSET_M, state.geometry.y_plane_m + y, 0.0),
                None,
                "OK",
            )
        )

    crack_x_min = state.geometry.crack_x_min_m
    crack_z_min, crack_z_max = state.geometry.crack_z_bounds_m
    tolerance = state.geometry.tolerance_m
    unknowns = (
        ("front-interior", (state.front_x_m, state.geometry.y_plane_m, 0.0), None, "front_or_finite_front_endpoint"),
        ("front-z-min-endpoint", (state.front_x_m, state.geometry.y_plane_m, crack_z_min), None, "front_or_finite_front_endpoint"),
        ("front-z-max-endpoint", (state.front_x_m, state.geometry.y_plane_m, crack_z_max), None, "front_or_finite_front_endpoint"),
        ("gamma-x-start-edge", (crack_x_min, state.geometry.y_plane_m, 0.0), mf.TraceSide.PLUS, "finite_crack_end_or_edge"),
        ("gamma-z-edge", (state.front_x_m - 0.004, state.geometry.y_plane_m, crack_z_max), mf.TraceSide.MINUS, "finite_crack_end_or_edge"),
        ("gamma-missing-side", (state.front_x_m - 0.004, state.geometry.y_plane_m, 0.0), None, "explicit_plus_minus_trace_required"),
        ("gamma-near-cut-band", (state.front_x_m - 0.004, state.geometry.y_plane_m + 0.5 * tolerance, 0.0), None, "ambiguous_crack_plane_tolerance_band"),
        ("outside-box", (0.0101, 0.003, 0.0), None, "outside_closed_domain"),
    )
    unknown_indices = []
    for label, point, side, reason in unknowns:
        unknown_indices.append(add(label, point, side, "UNKNOWN", reason))

    _require(len(points) == PROVIDER_POINTS_PER_STATE, "query layout is not exactly 25 rows")
    _require(len(set((point, None if side is None else side.value) for point, side in zip(points, sides))) == len(points), "query identities are not unique")
    _require(len(core_indices) == 8 and len(trace_indices) == 4 and len(collar_indices) == 2 and len(ahead_indices) == 3 and len(unknown_indices) == 8, "frozen query-category counts changed")
    for index in core_indices:
        _require(_distance_to_front(points[index], state) < r_core_m, "core point is not strictly inside its measured core")
    return points, sides, specs, tuple(core_indices)


def _verify_position_storage(state: mf.FrontState, bundle: mf.InputBundle) -> None:
    array = state.source_positions
    _require(array.shape == (SOURCE_COUNT, 3), "source-position shape changed")
    _require(not array.flags.writeable and not array.flags.owndata, "source positions are not a non-owning read-only view")
    _require(np.array_equal(array, np.asarray([source.position_m for source in state.sources], dtype=np.float64)), "source-position bytes do not match classified material sources")
    _require(tuple(source.position_m for source in state.sources) == tuple(source.position_m for source in bundle.rest_sources), "front state changed material-rest source coordinates or order")
    _require(len({tuple(row) for row in array.tolist()}) == SOURCE_COUNT, "material-rest source coordinates are not unique")
    try:
        array.setflags(write=True)
    except ValueError:
        pass
    else:
        array.setflags(write=False)
        raise CheckFailure("source-position array allowed WRITEABLE to be re-enabled")
    _require(not array.flags.writeable, "source-position array became writeable during immutability probe")


def _verify_state_contract(state: mf.FrontState, bundle: mf.InputBundle) -> None:
    _verify_position_storage(state, bundle)
    _require(len(state.sources) == SOURCE_COUNT, "current state does not contain exactly 100 sources")
    _require(len({source.source_id for source in state.sources}) == SOURCE_COUNT, "current source IDs are not unique")
    _require(all(source.position_m[1] != 0.0 for source in state.sources), "a shifted-F1 source lies on y=0")
    _require(state.registry.source_count == SOURCE_COUNT, "registry source count differs from the material-rest set")
    _require(state.registry.dof_count <= 375, "current scalar DOF count exceeds 375")
    _require(sum(state.role_counts.values()) == SOURCE_COUNT, "current geometric role counts do not sum to 100")
    _require(set(state.role_counts) == {"FRONT", "HEAVISIDE_ELIGIBLE", "BASE"}, "current role keys differ from the frozen contract")
    for source in state.sources:
        _require(
            len(source.basis_roles) in (1, 2)
            and source.basis_roles[0] == "BASE"
            and source.basis_roles.count("BASE") == 1
            and all(role in ("HEAVISIDE", "FRONT_F1_SHIFTED") for role in source.basis_roles[1:]),
            f"invalid or redundant basis roles for {source.source_id}",
        )
    registered_base_ids = {
        key.source_id for key in state.registry.ordered_keys if key.basis_kind is mf.BasisKind.BASE
    }
    _require(registered_base_ids == {source.source_id for source in state.sources}, "registry lacks exactly one BASE key per material source")
    _require(state.registry.dof_count == sum(len(source.basis_roles) for source in state.sources), "registry key count suggests a duplicate or omitted basis role")
    permuted = mf._new_registry(tuple(reversed(state.sources)))
    _require(permuted.canonical_ids == state.registry.canonical_ids, "semantic coefficient IDs changed under already-classified source permutation")
    _require(permuted.source_count == state.registry.source_count and permuted.dof_count == state.registry.dof_count, "permuted registry counts changed")


def _sample(provider: mf.MovingFrontProvider, state: mf.FrontState, index: int) -> mf.MovingFrontSample:
    sample = provider.sample_at(index)
    tag = sample.provenance
    _require(tag.case_id == state.case_id, f"sample {index} case tag mismatch")
    _require(tag.front_x_m == state.front_x_m, f"sample {index} front tag mismatch")
    _require(tag.surface_geometry_sha256 == state.surface_geometry_sha256, f"sample {index} geometry tag mismatch")
    _require(tag.basis_generation_sha256 == state.basis_generation_sha256, f"sample {index} basis tag mismatch")
    _require(tag.source_coordinate_manifest_sha256 == state.source_coordinate_manifest_sha256, f"sample {index} coordinate tag mismatch")
    _require(tag.safe_classifier_sha256 == state.safe_classifier_sha256, f"sample {index} classifier tag mismatch")
    _require(tag.canonical_ids_sha256 == state.canonical_ids_sha256, f"sample {index} canonical-ID tag mismatch")
    if sample.row.status is mf.ProviderStatus.OK:
        _require(sample.row.coefficient_ids == state.registry.canonical_ids, f"sample {index} semantic IDs differ from registry")
    return sample


def _affine_coefficients(state: mf.FrontState) -> dict[str, float]:
    by_id = {source.source_id: source for source in state.sources}
    values: dict[str, float] = {}
    for key in state.registry.ordered_keys:
        source = by_id[key.source_id]
        if key.basis_kind is mf.BasisKind.BASE:
            x, y, z = source.position_m
            values[key.canonical_id] = AFFINE_INTERCEPT_M + sum(a * b for a, b in zip(AFFINE_GRADIENT, (x, y, z)))
        else:
            values[key.canonical_id] = 0.0
    return values


def _manufactured_coefficients(state: mf.FrontState, bundle: mf.InputBundle) -> dict[str, float]:
    material = bundle.material_case["material_and_manufactured_field"]
    young = float(material["young_modulus_Pa"])
    poisson = float(material["poisson_ratio"])
    stress_intensity = float(material["K_III_Pa_sqrt_m"])
    mu = young / (2.0 * (1.0 + poisson))
    amplitude = (stress_intensity / mu) * math.sqrt(2.0 / math.pi)
    coefficient = amplitude * math.sqrt(mf.R0_M)
    by_id = {source.source_id: source for source in state.sources}
    values: dict[str, float] = {}
    for key in state.registry.ordered_keys:
        source = by_id[key.source_id]
        if key.basis_kind is mf.BasisKind.BASE:
            values[key.canonical_id] = (
                coefficient * float(source.shifted_fbar_source)
                if source.basis_role == "FRONT_F1_SHIFTED" and source.shifted_fbar_source is not None
                else 0.0
            )
        elif key.basis_kind is mf.BasisKind.FRONT_F1_SHIFTED:
            values[key.canonical_id] = coefficient
        else:
            values[key.canonical_id] = 0.0
    return values


def _evaluate_row(sample: mf.MovingFrontSample, coefficient_by_id: Mapping[str, float]) -> tuple[float, tuple[float, float, float]]:
    _require(sample.row.status is mf.ProviderStatus.OK, f"cannot evaluate non-OK provider row: {sample.row.diagnostic}")
    value = 0.0
    gradient = [0.0, 0.0, 0.0]
    for coefficient_id, contribution in zip(sample.row.coefficient_ids, sample.row.contributions):
        coefficient = coefficient_by_id[coefficient_id]
        value += coefficient * contribution.value
        for axis in range(3):
            gradient[axis] += coefficient * contribution.gradient[axis]
    return value, (gradient[0], gradient[1], gradient[2])


def _check_affine_patch(
    batch: mf.MovingFrontProvider,
    state: mf.FrontState,
    query_specs: Sequence[Mapping[str, Any]],
) -> dict[str, float | int]:
    coefficients = _affine_coefficients(state)
    worst_value = 0.0
    worst_gradient = 0.0
    checked = 0
    for index, query in enumerate(query_specs):
        if query["expected"] != "OK":
            continue
        sample = _sample(batch, state, index)
        _require(sample.row.status is mf.ProviderStatus.OK, f"valid affine query {query['label']} returned UNKNOWN: {sample.row.diagnostic}")
        actual, grad = _evaluate_row(sample, coefficients)
        point = query["point_m"]
        expected = AFFINE_INTERCEPT_M + sum(a * b for a, b in zip(AFFINE_GRADIENT, point))
        value_error = abs(actual - expected)
        gradient_error = max(abs(a - b) for a, b in zip(grad, AFFINE_GRADIENT))
        worst_value = max(worst_value, value_error)
        worst_gradient = max(worst_gradient, gradient_error)
        _require(value_error <= AFFINE_VALUE_TOL_M, f"affine value error at {query['label']} exceeds {AFFINE_VALUE_TOL_M:g} m")
        _require(gradient_error <= AFFINE_GRAD_TOL, f"affine gradient error at {query['label']} exceeds {AFFINE_GRAD_TOL:g}")
        checked += 1
        _deadline("affine BASE patch checks")
    return {
        "valid_rows_checked": checked,
        "value_tolerance_abs_m": AFFINE_VALUE_TOL_M,
        "gradient_tolerance_abs_dimensionless": AFFINE_GRAD_TOL,
        "maximum_value_error_m": worst_value,
        "maximum_gradient_error_dimensionless": worst_gradient,
    }


def _check_core_manufactured(
    batch: mf.MovingFrontProvider,
    state: mf.FrontState,
    bundle: mf.InputBundle,
    points: Sequence[Sequence[float]],
    core_indices: Sequence[int],
    r_core: float,
) -> dict[str, float | int]:
    coefficients = _manufactured_coefficients(state, bundle)
    material = bundle.material_case["material_and_manufactured_field"]
    worst_value = 0.0
    worst_gradient = 0.0
    for index in core_indices:
        sample = _sample(batch, state, index)
        _require(sample.row.status is mf.ProviderStatus.OK, f"core query {index} returned UNKNOWN: {sample.row.diagnostic}")
        actual, gradient = _evaluate_row(sample, coefficients)
        expected, expected_gradient = mf.analytic_modeiii_reference(
            points[index],
            state.front_x_m,
            float(material["young_modulus_Pa"]),
            float(material["poisson_ratio"]),
            float(material["K_III_Pa_sqrt_m"]),
        )
        value_error = abs(actual - expected)
        gradient_error = max(abs(a - b) for a, b in zip(gradient, expected_gradient))
        worst_value = max(worst_value, value_error)
        worst_gradient = max(worst_gradient, gradient_error)
        _require(value_error <= CORE_VALUE_TOL_M, f"current-front core displacement error exceeds {CORE_VALUE_TOL_M:g} m")
        _require(gradient_error <= CORE_GRAD_TOL, f"current-front core gradient error exceeds {CORE_GRAD_TOL:g}")
        _deadline("core manufactured-field checks")
    return {
        "points_checked": len(core_indices),
        "measured_core_radius_m": r_core,
        "displacement_tolerance_abs_m": CORE_VALUE_TOL_M,
        "gradient_tolerance_abs_dimensionless": CORE_GRAD_TOL,
        "maximum_displacement_error_m": worst_value,
        "maximum_gradient_error_dimensionless": worst_gradient,
    }


def _one_sided_fbar(point: Sequence[float], state: mf.FrontState, side: mf.TraceSide) -> tuple[float, tuple[float, float, float]]:
    radius = state.front_x_m - float(point[0])
    _require(radius > 0.0, "trace fixture is not behind its current front")
    sign = 1.0 if side is mf.TraceSide.PLUS else -1.0
    value = sign * math.sqrt(radius / mf.R0_M)
    x_gradient = -sign / (2.0 * math.sqrt(mf.R0_M * radius))
    return value, (x_gradient, 0.0, 0.0)


def _check_exact_trace_basis(
    batch: mf.MovingFrontProvider,
    state: mf.FrontState,
    query_specs: Sequence[Mapping[str, Any]],
) -> dict[str, float | int]:
    source_by_id = {source.source_id: source for source in state.sources}
    worst_value = 0.0
    worst_gradient = 0.0
    checked = 0
    for index, query in enumerate(query_specs):
        if not str(query["label"]).startswith("trace-"):
            continue
        sample = _sample(batch, state, index)
        _require(sample.row.status is mf.ProviderStatus.OK, f"exact trace {query['label']} returned UNKNOWN: {sample.row.diagnostic}")
        side = mf.TraceSide(query["trace_side"])
        h_query = 1.0 if side is mf.TraceSide.PLUS else -1.0
        fbar_query, grad_fbar_query = _one_sided_fbar(query["point_m"], state, side)
        mls_row = batch._mls_row_by_query[index]
        phi_all = batch.mls_result.weights[mls_row]
        grad_all = batch.mls_result.gradients[mls_row]
        for slot, key in enumerate(state.registry.ordered_keys):
            source = source_by_id[key.source_id]
            source_row = state.source_row_by_id[key.source_id]
            phi = float(phi_all[source_row])
            grad_phi = grad_all[source_row]
            if key.basis_kind is mf.BasisKind.BASE:
                expected_value, expected_gradient = phi, grad_phi
            elif key.basis_kind is mf.BasisKind.HEAVISIDE:
                _require(source.source_heaviside is not None, "H trace source lacks its current-side value")
                factor = h_query - float(source.source_heaviside)
                expected_value, expected_gradient = phi * factor, factor * grad_phi
            elif key.basis_kind is mf.BasisKind.FRONT_F1_SHIFTED:
                _require(source.shifted_fbar_source is not None, "F1 trace source lacks its current shifted value")
                difference = fbar_query - float(source.shifted_fbar_source)
                expected_value = phi * difference
                expected_gradient = difference * grad_phi + phi * np.asarray(grad_fbar_query, dtype=np.float64)
            else:
                raise CheckFailure(f"unsupported trace basis kind: {key.basis_kind.value}")
            actual = sample.row.contributions[slot]
            value_error = abs(actual.value - expected_value)
            gradient_error = max(abs(actual.gradient[a] - float(expected_gradient[a])) for a in range(3))
            worst_value = max(worst_value, value_error)
            worst_gradient = max(worst_gradient, gradient_error)
            _require(math.isfinite(actual.value) and all(math.isfinite(v) for v in actual.gradient), "exact trace basis row is nonfinite")
            _require(value_error <= TRACE_BASIS_VALUE_TOL, "exact trace basis value differs from its analytic one-sided branch")
            _require(gradient_error <= TRACE_BASIS_GRAD_TOL_PER_M, "exact trace basis gradient differs from its analytic one-sided branch")
        checked += 1
        _deadline("exact trace basis checks")
    _require(checked == 4, "exact finite-Gamma trace row count changed")
    return {
        "trace_rows_checked": checked,
        "basis_value_tolerance_abs_dimensionless": TRACE_BASIS_VALUE_TOL,
        "basis_gradient_tolerance_abs_per_m": TRACE_BASIS_GRAD_TOL_PER_M,
        "maximum_basis_value_residual": worst_value,
        "maximum_basis_gradient_residual_per_m": worst_gradient,
        "complete_field_exactness_claimed": False,
    }


def _check_unknowns_and_tags(
    batch: mf.MovingFrontProvider,
    state: mf.FrontState,
    query_specs: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    audit: list[dict[str, Any]] = []
    for index, query in enumerate(query_specs):
        sample = _sample(batch, state, index)
        status = sample.row.status.value
        reason = batch.row_reasons[index]
        if query["expected"] == "UNKNOWN":
            _require(sample.row.status is mf.ProviderStatus.UNKNOWN, f"negative query {query['label']} was not UNKNOWN")
            _require(reason == query["expected_reason"], f"negative query {query['label']} had reason {reason!r}")
        else:
            _require(batch.classifications[index][0] is mf.ProviderStatus.OK, f"expected-valid query {query['label']} failed geometric classification")
        audit.append(
            {
                "label": query["label"],
                "point_m": query["point_m"],
                "trace_side": query["trace_side"],
                "expected": query["expected"],
                "status": status,
                "reason": reason,
                "tags_match_current_state": True,
            }
        )
    return audit


def _check_ahead_and_collar(
    batch: mf.MovingFrontProvider,
    state: mf.FrontState,
    bundle: mf.InputBundle,
    query_specs: Sequence[Mapping[str, Any]],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    coeffs = _manufactured_coefficients(state, bundle)
    ahead: list[tuple[int, mf.MovingFrontSample]] = []
    h_max = 0.0
    for index, query in enumerate(query_specs):
        if not str(query["label"]).startswith("ahead-"):
            continue
        sample = _sample(batch, state, index)
        _require(sample.row.status is mf.ProviderStatus.OK, f"ahead finite-sample query {query['label']} returned UNKNOWN")
        h_contributions = [
            contribution
            for key, contribution in zip(state.registry.ordered_keys, sample.row.contributions)
            if key.basis_kind is mf.BasisKind.HEAVISIDE
        ]
        for contribution in h_contributions:
            h_max = max(h_max, abs(contribution.value), *(abs(value) for value in contribution.gradient))
        _require(all(contribution.value == 0.0 and all(value == 0.0 for value in contribution.gradient) for contribution in h_contributions), "ahead-of-front H value/gradient is nonzero")
        ahead.append((index, sample))
    _require(len(ahead) == 3, "ahead finite-sample count changed")
    states = []
    for index, sample in ahead:
        value, gradient = _evaluate_row(sample, coeffs)
        states.append(
            {
                "label": query_specs[index]["label"],
                "field_value_m": value,
                "field_gradient_dimensionless": list(gradient),
            }
        )
    value_spread = max(item["field_value_m"] for item in states) - min(item["field_value_m"] for item in states)
    gradient_spread = max(
        max(item["field_gradient_dimensionless"][axis] for item in states)
        - min(item["field_gradient_dimensionless"][axis] for item in states)
        for axis in range(3)
    )
    continuity = {
        "sample_rows": states,
        "heaviside_max_abs_value_or_gradient": h_max,
        "manufactured_field_value_spread_m": value_spread,
        "manufactured_field_gradient_component_spread_dimensionless": gradient_spread,
        "scope": "finite-sample diagnostic only; not a global continuity proof",
    }

    collar: list[dict[str, Any]] = []
    material = bundle.material_case["material_and_manufactured_field"]
    for index, query in enumerate(query_specs):
        if not str(query["label"]).startswith("collar-"):
            continue
        sample = _sample(batch, state, index)
        _require(sample.row.status is mf.ProviderStatus.OK, f"collar observation {query['label']} returned UNKNOWN")
        actual, gradient = _evaluate_row(sample, coeffs)
        reference, reference_gradient = mf.analytic_modeiii_reference(
            query["point_m"],
            state.front_x_m,
            float(material["young_modulus_Pa"]),
            float(material["poisson_ratio"]),
            float(material["K_III_Pa_sqrt_m"]),
        )
        collar.append(
            {
                "label": query["label"],
                "point_m": query["point_m"],
                "manufactured_value_m": actual,
                "analytic_reference_m": reference,
                "absolute_value_error_m": abs(actual - reference),
                "manufactured_gradient_dimensionless": list(gradient),
                "analytic_gradient_dimensionless": list(reference_gradient),
                "maximum_gradient_error_dimensionless": max(abs(a - b) for a, b in zip(gradient, reference_gradient)),
                "acceptance_threshold": None,
            }
        )
    _require(len(collar) == 2, "collar observation count changed")
    return continuity, collar


def _array_bytes(state: mf.FrontState, batch: mf.MovingFrontProvider, previous_state: mf.FrontState | None) -> int:
    arrays = [state.source_positions]
    if previous_state is not None:
        arrays.append(previous_state.source_positions)
    arrays.extend((batch.values, batch.gradients))
    result = batch.mls_result
    if result is not None:
        arrays.extend((result.weights, result.gradients, result.support_counts, result.condition_numbers))
    return sum(int(array.nbytes) for array in arrays)


def _control_tangencies(state: mf.FrontState, safe: Any) -> list[dict[str, Any]]:
    geometry = state.geometry
    h = geometry.support_radius_m
    tolerance = geometry.tolerance_m
    probes = (
        ("exact", geometry.front_x_m + h),
        ("within-tolerance", geometry.front_x_m + h + 0.5 * tolerance),
    )
    records: list[dict[str, Any]] = []
    for label, x in probes:
        source_id = f"control-{state.case_id}-{label}"
        record = safe.classify_support(
            source_id,
            (x, geometry.y_plane_m, 0.0),
            geometry,
        )
        _require(record.get("status") == "UNKNOWN", f"{label} tangency control was not UNKNOWN")
        _require(record.get("reason_code") == "support_tangent_to_finite_front", f"{label} tangency control had unexpected reason")
        records.append(
            {
                "control_id": source_id,
                "position_m": [x, geometry.y_plane_m, 0.0],
                "status": record["status"],
                "reason": record["reason_code"],
                "enters_material_source_set": False,
                "enters_mls_batch": False,
            }
        )
    return records


def _run_state(
    bundle: mf.InputBundle,
    front_x_m: float,
    safe: Any,
    wrappers: dict[str, Any],
    previous_state: mf.FrontState | None,
) -> tuple[mf.FrontState, dict[str, Any]]:
    before_sources = wrappers["counts"]["source_classifications"]
    state = mf.build_front_state(bundle, front_x_m)
    state_classification_count = wrappers["counts"]["source_classifications"] - before_sources
    _require(state_classification_count == SOURCE_COUNT, "each front must reclassify exactly 100 material-rest sources")
    _verify_state_contract(state, bundle)
    r_core = _current_core_radius(state)
    points, sides, query_specs, core_indices = _make_queries(state, r_core)
    controls_before = wrappers["counts"]["tangency_controls"]
    controls = _control_tangencies(state, safe)
    _require(wrappers["counts"]["tangency_controls"] - controls_before == 2, "each front must have exactly two pure topology controls")

    mls_before = wrappers["counts"]["mls_batches"]
    batch = mf.MovingFrontProvider(state, points, sides)
    _require(wrappers["counts"]["mls_batches"] - mls_before == 1, "each state must make exactly one MLS batch")
    _require(len(batch.query_points_m) == PROVIDER_POINTS_PER_STATE, "provider batch did not receive exactly 25 point/trace identities")
    _require(sum(status is mf.ProviderStatus.OK for status in batch.row_statuses) == 17, "provider did not classify the 17 planned query rows as geometrically valid")
    _require(batch.mls_result is not None, "provider skipped its one planned MLS batch")
    _require(len(batch.mls_result.statuses) == 17, "MLS batch contains an unexpected number of valid rows")
    for index, query in enumerate(query_specs):
        if query["expected"] == "OK":
            _require(batch.row_statuses[index] is mf.ProviderStatus.OK, f"planned valid row failed MLS/provider: {query['label']} ({batch.row_reasons[index]})")
        else:
            _require(batch.row_statuses[index] is mf.ProviderStatus.UNKNOWN, f"planned negative was not UNKNOWN: {query['label']}")

    affine = _check_affine_patch(batch, state, query_specs)
    core = _check_core_manufactured(batch, state, bundle, points, core_indices, r_core)
    traces = _check_exact_trace_basis(batch, state, query_specs)
    query_audit = _check_unknowns_and_tags(batch, state, query_specs)
    continuity, collar = _check_ahead_and_collar(batch, state, bundle, query_specs)
    array_bytes = _array_bytes(state, batch, previous_state)
    _require(array_bytes + TRANSIENT_RESERVE_BYTES + MAX_RESULT_JSON_BYTES <= MAX_DIAGNOSTIC_BYTES, "diagnostic-owned workspace can exceed 2 MiB")
    current_canonical_ids_sha = state.canonical_ids_sha256
    current_provenance = mf.CoefficientProvenance(
        case_id=state.case_id,
        front_x_m=state.front_x_m,
        surface_geometry_sha256=state.surface_geometry_sha256,
        basis_generation_sha256=state.basis_generation_sha256,
        canonical_ids_sha256=current_canonical_ids_sha,
        source_solution_id=f"source-only-manufactured-control-{state.case_id}",
    )
    mf.require_current_coefficient_provenance(state, current_provenance)
    role_summary = {
        "geometric_role_counts": dict(state.role_counts),
        "source_count": state.registry.source_count,
        "scalar_dof_count": state.registry.dof_count,
        "canonical_ids_sha256": current_canonical_ids_sha,
        "semantic_ids_stable_under_reversed_classified_source_order": True,
        "permutation_reclassification_count": 0,
        "permutation_mls_batch_count": 0,
    }
    report = {
        "case_id": state.case_id,
        "front_x_m": state.front_x_m,
        "surface_geometry_sha256": state.surface_geometry_sha256,
        "basis_generation_sha256": state.basis_generation_sha256,
        "source_coordinate_manifest_sha256": state.source_coordinate_manifest_sha256,
        "safe_classifier_sha256": state.safe_classifier_sha256,
        "role_and_dof_summary": role_summary,
        "source_reclassifications": state_classification_count,
        "measured_core_radius_m": r_core,
        "core_probe_radii_m": [
            min(cap, fraction * r_core)
            for cap, fraction in zip(CORE_RADII_CAPS_M, CORE_RADII_FRACTIONS)
        ],
        "affine_base_patch": affine,
        "current_front_manufactured_core": core,
        "exact_trace_basis": traces,
        "ahead_finite_sample_diagnostic": continuity,
        "collar_observations": collar,
        "provider_point_count": len(points),
        "provider_batch_count": 1,
        "mls_batch_count": 1,
        "sample_provenance_tags_checked": len(points),
        "query_audit": query_audit,
        "topology_tangency_controls": controls,
        "diagnostic_owned_array_bytes": array_bytes,
        "same_state_provenance_identity_accepted": True,
        "mechanical_solution_or_coefficient_transfer": False,
    }
    _deadline(f"completed state front_x={front_x_m:g}")
    return state, report


def _install_counting_wrappers() -> dict[str, Any]:
    counts = {"source_classifications": 0, "tangency_controls": 0, "mls_batches": 0}
    safe = mf.load_safe_v2_classifier()
    original_classify = safe.classify_support
    original_shape = mf.shape_functions_3d

    def classify_counted(source_id: str, *args: Any, **kwargs: Any) -> Any:
        if source_id.startswith("control-"):
            counts["tangency_controls"] += 1
            if counts["tangency_controls"] > TANGENCY_CONTROLS_TOTAL:
                raise CheckFailure("pure topology control count exceeded its frozen cap")
        else:
            counts["source_classifications"] += 1
            if counts["source_classifications"] > SOURCE_RECLASSIFICATIONS_TOTAL:
                raise CheckFailure("material source reclassification count exceeded 200")
        return original_classify(source_id, *args, **kwargs)

    def shape_counted(*args: Any, **kwargs: Any) -> Any:
        counts["mls_batches"] += 1
        if counts["mls_batches"] > 2:
            raise CheckFailure("an extra MLS batch was requested")
        return original_shape(*args, **kwargs)

    safe.classify_support = classify_counted
    mf.shape_functions_3d = shape_counted
    return {
        "counts": counts,
        "safe": safe,
        "original_classify": original_classify,
        "original_shape": original_shape,
    }


def _restore_counting_wrappers(wrappers: Mapping[str, Any]) -> None:
    wrappers["safe"].classify_support = wrappers["original_classify"]
    mf.shape_functions_3d = wrappers["original_shape"]


def run_check(frozen: tuple[dict[str, Any], str, str] | None = None) -> dict[str, Any]:
    _require(RUN_MARKER_PATH.is_file(), "use the one-shot script entrypoint; no execution marker was claimed")
    _require(not OUTPUT_PATH.exists(), "check result already exists; refusing a repeat run")
    spec, spec_hash, driver_hash = _load_frozen_spec() if frozen is None else frozen
    bundle = mf.load_pinned_inputs()
    _require(len(bundle.rest_sources) == SOURCE_COUNT, "pinned coordinate binding does not contain 100 material-rest sources")
    _require(
        bundle.coordinate_manifest.get("source_records_sha256") == SOURCE_COORDINATE_RECORDS_SHA256,
        "coordinate-only source record digest changed",
    )
    _deadline("pinned provider inputs loaded")
    wrappers = _install_counting_wrappers()
    first_state: mf.FrontState | None = None
    baseline_provenance: mf.CoefficientProvenance | None = None
    state_reports: list[dict[str, Any]] = []
    try:
        front_values = tuple(float(row["front_x_m"]) for row in bundle.case_pair["moving_front_cases"])
        _require(front_values == mf.EXPECTED_FRONT_X_M, "frozen two-front case order or values changed")
        for ordinal, front_x in enumerate(front_values):
            previous = first_state
            state, report = _run_state(bundle, front_x, wrappers["safe"], wrappers, previous)
            state_reports.append(report)
            if ordinal == 0:
                first_state = state
                baseline_provenance = mf.CoefficientProvenance(
                    case_id=state.case_id,
                    front_x_m=state.front_x_m,
                    surface_geometry_sha256=state.surface_geometry_sha256,
                    basis_generation_sha256=state.basis_generation_sha256,
                    canonical_ids_sha256=state.canonical_ids_sha256,
                    source_solution_id="baseline-provenance-control-no-coefficients",
                )
                mf.require_current_coefficient_provenance(state, baseline_provenance)
            else:
                _require(first_state is not None and baseline_provenance is not None, "baseline state/provenance is absent")
                mf.assert_distinct_front_generations(first_state, state)
                try:
                    mf.require_current_coefficient_provenance(state, baseline_provenance)
                except mf.StaleBasisProvenanceError:
                    stale_rejected = True
                else:
                    raise CheckFailure("cross-front provenance was accepted")
                state_reports[-1]["cross_front_stale_provenance_rejected"] = stale_rejected
                state_reports[-1]["prior_surface_geometry_sha256"] = first_state.surface_geometry_sha256
                state_reports[-1]["prior_basis_generation_sha256"] = first_state.basis_generation_sha256
                state_reports[-1]["no_coefficient_data_loaded"] = True
            _deadline("between front states")
        _require(len(state_reports) == 2, "exactly two front states were not checked")
        counts = wrappers["counts"]
        _require(counts["source_classifications"] == SOURCE_RECLASSIFICATIONS_TOTAL, "total safe-v2 material source reclassifications differ from 200")
        _require(counts["tangency_controls"] == TANGENCY_CONTROLS_TOTAL, "total pure topology controls differ from four")
        _require(counts["mls_batches"] == 2, "total MLS batches differ from two")
        provider_points = sum(report["provider_point_count"] for report in state_reports)
        _require(provider_points == PROVIDER_POINTS_TOTAL, "total provider point/trace identities differ from 50")
        records = counts["source_classifications"] + counts["tangency_controls"] + provider_points
        _require(records == RECORDS_TOTAL and records <= RECORDS_HARD_CAP, "strict 254-record budget was violated")
        _require(state_reports[0]["surface_geometry_sha256"] != state_reports[1]["surface_geometry_sha256"], "surface hashes failed to differ")
        _require(state_reports[0]["basis_generation_sha256"] != state_reports[1]["basis_generation_sha256"], "basis-generation hashes failed to differ")

        result = {
            "schema": "moving-front-provider-check-result.v1",
            "status": "PASS",
            "check_spec_sha256": spec_hash,
            "driver_source_sha256": driver_hash,
            "provider_source_sha256": PROVIDER_SHA256,
            "case_pair_sha256": CASES_SHA256,
            "frozen_input_hashes": {name: value["sha256"] for name, value in spec["source_files"].items()},
            "budget": {
                "source_reclassifications": counts["source_classifications"],
                "pure_topology_controls": counts["tangency_controls"],
                "provider_point_trace_rows": provider_points,
                "mls_batches": counts["mls_batches"],
                "total_records": records,
                "record_hard_cap": RECORDS_HARD_CAP,
                "maximum_wall_seconds": MAX_WALL_SECONDS,
            },
            "states": state_reports,
            "cross_front_hashes_differ": True,
            "cross_front_old_provenance_rejected": bool(state_reports[1].get("cross_front_stale_provenance_rejected")),
            "old_coefficients_read_or_migrated": False,
            "operator_or_solver_run": False,
            "scope": "bounded source-provider diagnostic only; no loaded mechanical solution or growth claim",
        }
    finally:
        _restore_counting_wrappers(wrappers)

    elapsed = time.monotonic() - _STARTED_MONOTONIC
    _require(elapsed <= MAX_WALL_SECONDS, f"30-second wall-time budget exceeded: {elapsed:.6f} s")
    result["elapsed_wall_seconds"] = elapsed
    active_arrays = max(report["diagnostic_owned_array_bytes"] for report in state_reports)
    result["memory_budget"] = {
        "peak_owned_numpy_array_bytes": active_arrays,
        "transient_reserve_bytes": TRANSIENT_RESERVE_BYTES,
        "serialized_result_bytes": 0,
        "estimated_diagnostic_owned_bytes": 0,
        "maximum_diagnostic_owned_bytes": MAX_DIAGNOSTIC_BYTES,
        "interpreter_import_baseline_included": False,
    }
    result["elapsed_wall_seconds"] = time.monotonic() - _STARTED_MONOTONIC
    _deadline("final result serialization")
    for _ in range(3):
        encoded = json.dumps(result, sort_keys=True, indent=2, ensure_ascii=True, allow_nan=False).encode("utf-8") + b"\n"
        _require(len(encoded) <= MAX_RESULT_JSON_BYTES, "final serialized result exceeds 256 KiB")
        estimated_bytes = active_arrays + TRANSIENT_RESERVE_BYTES + len(encoded)
        _require(estimated_bytes <= MAX_DIAGNOSTIC_BYTES, "final diagnostic-owned working set exceeds 2 MiB")
        if result["memory_budget"]["serialized_result_bytes"] == len(encoded) and result["memory_budget"]["estimated_diagnostic_owned_bytes"] == estimated_bytes:
            break
        result["memory_budget"]["serialized_result_bytes"] = len(encoded)
        result["memory_budget"]["estimated_diagnostic_owned_bytes"] = estimated_bytes
    final_encoded = json.dumps(result, sort_keys=True, indent=2, ensure_ascii=True, allow_nan=False).encode("utf-8") + b"\n"
    _require(result["memory_budget"]["serialized_result_bytes"] == len(final_encoded), "serialized memory accounting did not converge")
    _require(
        result["memory_budget"]["estimated_diagnostic_owned_bytes"]
        == active_arrays + TRANSIENT_RESERVE_BYTES + len(final_encoded),
        "estimated working-set accounting did not converge",
    )
    _require(len(final_encoded) <= MAX_RESULT_JSON_BYTES, "result JSON exceeds the frozen 256 KiB cap")
    _require(active_arrays + TRANSIENT_RESERVE_BYTES + len(final_encoded) <= MAX_DIAGNOSTIC_BYTES, "final estimated working set exceeds 2 MiB")
    _deadline("final result serialization completed")
    return result


def main() -> int:
    _require(not RUN_MARKER_PATH.exists(), "this frozen diagnostic has already been attempted")
    _require(not OUTPUT_PATH.exists(), "refusing to overwrite an existing check result")
    frozen = _load_frozen_spec()
    with RUN_MARKER_PATH.open("x", encoding="utf-8") as marker:
        marker.write(
            "single approved execution claimed; failures do not clear this marker\n"
            f"check_spec_sha256={frozen[1]}\n"
            f"driver_source_sha256={frozen[2]}\n"
        )
    result = run_check(frozen)
    encoded = json.dumps(result, sort_keys=True, indent=2, ensure_ascii=True, allow_nan=False).encode("utf-8") + b"\n"
    if len(encoded) > MAX_RESULT_JSON_BYTES:
        raise CheckFailure("result grew beyond the frozen JSON output limit")
    with OUTPUT_PATH.open("xb") as output:
        output.write(encoded)
    _deadline("result file write")
    print(f"{result['status']}: {OUTPUT_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
