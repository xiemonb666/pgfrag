"""Pinned source-only 3-D MLS provider for the 100-source local-front case.

This adapter is deliberately bound to one exact fixed case, classifier, role
manifest, registry contract, and MLS implementation. It does not generalize
the finite-crack geometry API and it is not coupled to an operator or solver.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any, Mapping, Sequence

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[3]
SUPPORT_DIR = REPO_ROOT / "pgfrag/experiments/local_front_support_contract_20261006"
REGISTRY_DIR = REPO_ROOT / "pgfrag/experiments/enriched_dof_contract_20261005"
MLS_DIR = REPO_ROOT / "pgfrag/experiments/mls3d_kernel_20261005"
MATERIAL_CASE_PATH = REPO_ROOT / "pgfrag/experiments/modeiii_pu_operator_20261005/fixed_case.json"

for _path in (str(SUPPORT_DIR), str(REGISTRY_DIR), str(MLS_DIR)):
    if _path not in sys.path:
        sys.path.insert(0, _path)

from enriched_dof_contract import (  # noqa: E402
    AntiplaneShapeSample,
    BasisKind,
    DofKey,
    DofRegistry,
    FieldComponent,
    ProviderContractError,
    ProviderStatus,
    RegistryState,
    ScalarShapeContribution,
    validate_antiplane_sample,
    validate_problem_budget,
)
from mls3d_mls import MLS3DConfig, MLSStatus  # noqa: E402
from mls3d_shape import ShapeFunctions3DResult, shape_functions_3d  # noqa: E402
from support_topology import (  # noqa: E402
    FiniteCrackGeometry,
    ManifestBlockedError,
    classify_support,
    verify_manifest_binding,
)


SUPPORT_CASE_REL = "pgfrag/experiments/local_front_support_contract_20261006/fixed_case.json"
CLASSIFIER_REL = "pgfrag/experiments/local_front_support_contract_20261006/support_topology.py"
MANIFEST_REL = "pgfrag/experiments/local_front_support_contract_20261006/role_manifest.json"
REGISTRY_REL = "pgfrag/experiments/enriched_dof_contract_20261005/enriched_dof_contract.py"
MLS_CORE_REL = "pgfrag/experiments/mls3d_kernel_20261005/mls3d_mls.py"
MLS_SHAPE_REL = "pgfrag/experiments/mls3d_kernel_20261005/mls3d_shape.py"
MATERIAL_CASE_REL = "pgfrag/experiments/modeiii_pu_operator_20261005/fixed_case.json"

# Raw file pins: replacing the geometry, classifier, manifest, MLS, registry,
# or manufactured-field fixture requires a new deliberate review and provider.
PINNED_SHA256 = {
    SUPPORT_CASE_REL: "0bdbb20032e105851d2909ff54850ebc54d3f2041c25ab5e42e9102244d90ca5",
    CLASSIFIER_REL: "ea20e51d2b774c7e6b8b48403359b56bab1ba6f5405dba41f34813b6aec4578e",
    MANIFEST_REL: "e5e318f80b5dfedbc8ba122f82df3cca1585e604272399d21bc35907d1d28f01",
    REGISTRY_REL: "e3ddb3b067465bb5b269576a1b8a5100e2bb92c9c924bc09c949702c19efa9d8",
    MLS_CORE_REL: "404256d9b7546361a5704c4c3cec7c022069a9cf83a137c58a16e7a64f79ec63",
    MLS_SHAPE_REL: "66e680616a988182875d17aeb84fadb541867aee13fc15e0ae098f7c98839ef5",
    MATERIAL_CASE_REL: "dce0a1f3009246a52198bc0ed833ed67998144549a653e7bce9bcc17babd71b9",
}
PINNED_MANIFEST_PAYLOAD_SHA256 = "b98799bfd6121370b7a8569c9bfc2e7e5958e8666687a797b0588b6412d10b35"
EXPECTED_CASE_ID = "modeiii_local_front_pu_support_roles_20261006"
EXPECTED_COUNTS = {"source_count": 100, "scalar_dof_count": 140}
EXPECTED_ROLE_COUNTS = {"BASE": 60, "FRONT": 30, "HEAVISIDE_ELIGIBLE": 10}
ROLE_TO_BASIS_ROLES = {
    "BASE": ("BASE",),
    "FRONT": ("BASE", "FRONT_F1_SHIFTED"),
    "HEAVISIDE_ELIGIBLE": ("BASE", "HEAVISIDE"),
}
R0_M = 0.01
FRONT_CORE_RADIUS_M = 0.001
MAX_BATCH_QUERIES = 256


class TraceSide(str, Enum):
    """Exact one-sided trace on the open finite crack patch."""

    PLUS = "PLUS"
    MINUS = "MINUS"


@dataclass(frozen=True)
class SourceRecord:
    source_id: str
    position_m: tuple[float, float, float]
    geometric_role: str
    basis_roles: tuple[str, ...]


@dataclass(frozen=True)
class InputBundle:
    case: Mapping[str, Any]
    manifest: Mapping[str, Any]
    sources: tuple[SourceRecord, ...]
    geometry: FiniteCrackGeometry
    field_material: Mapping[str, Any]
    mls_config: MLS3DConfig


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def pinned_source_hashes() -> dict[str, str]:
    """Return current dependency hashes, failing closed on any moved input."""
    actual: dict[str, str] = {}
    for relative_path, expected in PINNED_SHA256.items():
        path = REPO_ROOT / relative_path
        if not path.is_file():
            raise ManifestBlockedError(f"pinned source is missing: {relative_path}")
        digest = sha256_file(path)
        if digest != expected:
            raise ManifestBlockedError(f"pinned source changed: {relative_path}")
        actual[relative_path] = digest
    return actual


def _json_file(path: Path) -> tuple[dict[str, Any], bytes]:
    raw = path.read_bytes()
    try:
        obj = json.loads(raw.decode("utf-8-sig"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ManifestBlockedError(f"invalid frozen JSON file {path.name}: {exc}") from exc
    if not isinstance(obj, dict):
        raise ManifestBlockedError(f"frozen JSON root must be an object: {path.name}")
    return obj, raw


def load_pinned_inputs() -> InputBundle:
    """Load and validate the sole geometry/role fixture this provider accepts."""
    source_hashes = pinned_source_hashes()
    case, case_raw = _json_file(REPO_ROOT / SUPPORT_CASE_REL)
    manifest, _manifest_raw = _json_file(REPO_ROOT / MANIFEST_REL)
    material_case, _material_raw = _json_file(REPO_ROOT / MATERIAL_CASE_REL)
    classifier_raw = (REPO_ROOT / CLASSIFIER_REL).read_bytes()
    if sha256_bytes(case_raw) != PINNED_SHA256[SUPPORT_CASE_REL]:
        raise ManifestBlockedError("local-front geometry fixture is not the pinned fixed case")
    if sha256_bytes(classifier_raw) != PINNED_SHA256[CLASSIFIER_REL]:
        raise ManifestBlockedError("support classifier is not the pinned reviewed source")
    if manifest.get("manifest_sha256") != PINNED_MANIFEST_PAYLOAD_SHA256:
        raise ManifestBlockedError("role manifest payload digest is not the pinned digest")
    verify_manifest_binding(manifest, case_raw, classifier_raw)

    if case.get("case_id") != EXPECTED_CASE_ID:
        raise ManifestBlockedError("provider is bound to one exact local-front case ID")
    geometry_data = case.get("geometry", {})
    crack = geometry_data.get("crack", {})
    expected_geometry = {
        "domain_box_m": [[-0.01, 0.01], [-0.01, 0.01], [-0.008, 0.008]],
        "crack": {
            "x_min_m": -0.01,
            "front_x_m": 0.0,
            "y_plane_m": 0.0,
            "z_bounds_m": [-0.008, 0.008],
        },
    }
    if geometry_data.get("domain_box_m") != expected_geometry["domain_box_m"]:
        raise ManifestBlockedError("provider rejects any changed domain geometry")
    for name, expected in expected_geometry["crack"].items():
        if crack.get(name) != expected:
            raise ManifestBlockedError(f"provider rejects any changed crack geometry ({name})")
    if case.get("support", {}).get("radius_m") != 0.007:
        raise ManifestBlockedError("provider requires the frozen h=0.007 m support radius")
    if case.get("topology", {}).get("tolerance_m") != 1.0e-12:
        raise ManifestBlockedError("provider requires the frozen topology tolerance")
    if manifest.get("case_id") != EXPECTED_CASE_ID:
        raise ManifestBlockedError("role manifest names a different local-front case")
    if manifest.get("source_count") != EXPECTED_COUNTS["source_count"]:
        raise ManifestBlockedError("the source count is not the frozen count of 100")
    if manifest.get("scalar_dof_count") != EXPECTED_COUNTS["scalar_dof_count"]:
        raise ManifestBlockedError("the DOF count is not the frozen count of 140")
    if manifest.get("role_counts") != EXPECTED_ROLE_COUNTS:
        raise ManifestBlockedError("role counts differ from the frozen 30/10/60 split")
    if manifest.get("maximum_queries_per_batch") != MAX_BATCH_QUERIES:
        raise ManifestBlockedError("manifest query cap differs from the provider cap")

    raw_sources = manifest.get("sources")
    if not isinstance(raw_sources, list) or len(raw_sources) != 100:
        raise ManifestBlockedError("role manifest must contain exactly 100 source records")
    sources: list[SourceRecord] = []
    seen_ids: set[str] = set()
    for row in raw_sources:
        if not isinstance(row, dict):
            raise ManifestBlockedError("malformed source record")
        source_id = row.get("source_id")
        position = row.get("position_m")
        geometric_role = row.get("geometric_role")
        basis_roles = row.get("basis_roles")
        if (
            type(source_id) is not str
            or source_id in seen_ids
            or not isinstance(position, list)
            or len(position) != 3
            or not isinstance(basis_roles, list)
        ):
            raise ManifestBlockedError("source records need unique IDs, 3-D positions, and basis roles")
        if any(isinstance(value, bool) or not isinstance(value, (float, int)) for value in position):
            raise ManifestBlockedError(f"invalid coordinates for {source_id}")
        if any(not math.isfinite(float(value)) for value in position):
            raise ManifestBlockedError(f"non-finite coordinates for {source_id}")
        expected_roles = ROLE_TO_BASIS_ROLES.get(geometric_role)
        if expected_roles is None or tuple(basis_roles) != expected_roles:
            raise ManifestBlockedError(f"unsupported role assignment for {source_id}")
        if len(expected_roles) > 2:
            raise ManifestBlockedError(f"source receives more than one enrichment: {source_id}")
        if abs(float(position[1])) <= 1.0e-12:
            raise ManifestBlockedError(f"shifted source value is undefined on the crack plane: {source_id}")
        seen_ids.add(source_id)
        sources.append(
            SourceRecord(
                source_id=source_id,
                position_m=tuple(float(value) for value in position),
                geometric_role=geometric_role,
                basis_roles=expected_roles,
            )
        )
    if len(seen_ids) != 100:
        raise ManifestBlockedError("source IDs are not unique")

    material = material_case.get("material_and_manufactured_field", {})
    if (
        material.get("young_modulus_Pa") != 1_000_000.0
        or material.get("poisson_ratio") != 0.3
        or material.get("K_III_Pa_sqrt_m") != 100.0
        or material_case.get("geometry", {}).get("radius_m") != R0_M
    ):
        raise ManifestBlockedError("manufactured Mode III constants differ from their frozen fixture")
    mls = material_case.get("mls", {})
    config = MLS3DConfig(
        reference_origin=tuple(float(value) for value in mls.get("reference_origin_m", ())),
        reference_scale=float(mls.get("reference_scale_m")),
        support_radius=float(case["support"]["radius_m"]),
        max_condition_number=float(mls.get("max_condition_number")),
        rank_relative_tolerance=float(mls.get("rank_relative_tolerance")),
        switch_tolerance=float(mls.get("switch_tolerance")),
        duplicate_tolerance=float(mls.get("duplicate_tolerance")),
    )
    if config.reference_origin != (0.0, 0.0, 0.0) or config.reference_scale != R0_M:
        raise ManifestBlockedError("provider requires the pinned origin and R0 reference scale")
    _ = source_hashes
    return InputBundle(
        case=case,
        manifest=manifest,
        sources=tuple(sources),
        geometry=FiniteCrackGeometry.from_case(case),
        field_material=material,
        mls_config=config,
    )


def build_registry(sources: Sequence[SourceRecord], *, reverse_registration: bool = False) -> DofRegistry:
    """Build only the BASE plus manifest-authorized single enrichment keys."""
    registry = DofRegistry()
    source_rows = list(sources)
    if reverse_registration:
        source_rows.reverse()
    for source in source_rows:
        for role in source.basis_roles:
            if role == "BASE":
                kind, basis_index = BasisKind.BASE, 0
            elif role == "HEAVISIDE":
                kind, basis_index = BasisKind.HEAVISIDE, 1
            elif role == "FRONT_F1_SHIFTED":
                kind, basis_index = BasisKind.FRONT_F1_SHIFTED, 1
            else:
                raise ManifestBlockedError(f"unsupported manifest basis role {role!r}")
            registry.register(DofKey(source.source_id, FieldComponent.UZ, kind, basis_index))
    registry.freeze()
    validate_problem_budget(len(sources), registry.dof_count, 0)
    if len(sources) != 100 or registry.source_count != 100 or registry.dof_count != 140:
        raise ProviderContractError("registry does not match the 100-source / 140-DOF frozen manifest")
    return registry


@dataclass(frozen=True)
class QueryClassification:
    status: ProviderStatus
    reason: str
    h_value: float | None
    is_finite_crack_trace: bool


def classify_query(
    point_m: Sequence[float],
    trace_side: TraceSide | None,
    geometry: FiniteCrackGeometry,
) -> QueryClassification:
    """Preflight a point against the exact pinned finite front and crack patch."""
    if isinstance(point_m, (str, bytes)) or len(point_m) != 3:
        return QueryClassification(ProviderStatus.UNKNOWN, "invalid_point_shape", None, False)
    if any(isinstance(value, bool) or not isinstance(value, (int, float, np.integer, np.floating)) for value in point_m):
        return QueryClassification(ProviderStatus.UNKNOWN, "non_numeric_point", None, False)
    point = tuple(float(value) for value in point_m)
    if any(not math.isfinite(value) for value in point):
        return QueryClassification(ProviderStatus.UNKNOWN, "non_finite_point", None, False)
    if trace_side is not None and type(trace_side) is not TraceSide:
        return QueryClassification(ProviderStatus.UNKNOWN, "unsupported_trace_side", None, False)
    if not geometry.box.contains(point, 0.0):
        return QueryClassification(ProviderStatus.UNKNOWN, "outside_domain", None, False)

    x, y, z = point
    z_gap = max(geometry.crack_z_bounds_m[0] - z, 0.0, z - geometry.crack_z_bounds_m[1])
    distance_to_front = math.sqrt(
        (x - geometry.front_x_m) ** 2 + (y - geometry.y_plane_m) ** 2 + z_gap**2
    )
    if distance_to_front <= geometry.tolerance_m:
        return QueryClassification(ProviderStatus.UNKNOWN, "front_or_front_endpoint", None, False)

    near_branch = abs(y - geometry.y_plane_m) <= geometry.tolerance_m and x < geometry.front_x_m
    if near_branch and y != geometry.y_plane_m:
        return QueryClassification(ProviderStatus.UNKNOWN, "ambiguous_branch_tolerance_band", None, False)
    if y == geometry.y_plane_m and x < geometry.front_x_m:
        x_interior = geometry.crack_x_min_m + geometry.tolerance_m < x < geometry.front_x_m - geometry.tolerance_m
        z_interior = geometry.crack_z_bounds_m[0] + geometry.tolerance_m < z < geometry.crack_z_bounds_m[1] - geometry.tolerance_m
        if not (x_interior and z_interior):
            reason = "finite_crack_boundary_or_end"
            return QueryClassification(ProviderStatus.UNKNOWN, reason, None, False)
        if trace_side is None:
            return QueryClassification(ProviderStatus.UNKNOWN, "explicit_plus_minus_trace_required", None, False)
        h_value = 1.0 if trace_side is TraceSide.PLUS else -1.0
        return QueryClassification(ProviderStatus.OK, "finite_gamma_exact_trace", h_value, True)
    if trace_side is not None:
        return QueryClassification(ProviderStatus.UNKNOWN, "trace_side_only_valid_on_finite_gamma", None, False)

    # The positive-x plane is ahead of the finite crack. It is not classified
    # as an infinite crack plane; H has the neutral value there and certified
    # H supports must evaluate to zero by their compact-support clearance.
    if y == geometry.y_plane_m and x >= geometry.front_x_m:
        h_value = 0.0
    else:
        h_value = 1.0 if y > geometry.y_plane_m else -1.0
    return QueryClassification(ProviderStatus.OK, "open_domain_point", h_value, False)


def normalized_f1_and_gradient(
    point_m: Sequence[float],
    trace_side: TraceSide | None,
) -> tuple[float, tuple[float, float, float]]:
    """Return F1/sqrt(R0) and its analytic one-sided gradient."""
    x, y, _z = (float(value) for value in point_m)
    radius = math.hypot(x, y)
    if radius == 0.0:
        raise ValueError("F1 is undefined on the front")
    if y == 0.0 and x < 0.0:
        if trace_side is TraceSide.PLUS:
            theta = math.pi
        elif trace_side is TraceSide.MINUS:
            theta = -math.pi
        else:
            raise ValueError("an exact finite-crack trace requires PLUS or MINUS")
    else:
        theta = math.atan2(y, x)
    sine = math.sin(0.5 * theta)
    cosine = math.cos(0.5 * theta)
    value = math.sqrt(radius / R0_M) * sine
    denominator = 2.0 * math.sqrt(R0_M * radius)
    return value, (-sine / denominator, cosine / denominator, 0.0)


def modeiii_amplitude(material: Mapping[str, Any]) -> float:
    young = float(material["young_modulus_Pa"])
    poisson = float(material["poisson_ratio"])
    k_iii = float(material["K_III_Pa_sqrt_m"])
    mu = young / (2.0 * (1.0 + poisson))
    return (k_iii / mu) * math.sqrt(2.0 / math.pi)


def exact_modeiii_uz_and_gradient(
    point_m: Sequence[float],
    trace_side: TraceSide | None,
    amplitude_sqrt_m: float,
) -> tuple[float, tuple[float, float, float]]:
    fbar, grad_fbar = normalized_f1_and_gradient(point_m, trace_side)
    factor = amplitude_sqrt_m * math.sqrt(R0_M)
    return factor * fbar, tuple(factor * value for value in grad_fbar)


class LocalFrontProvider:
    """Evaluate dense BASE/H/shifted-F1 rows over the frozen 140-ID registry."""

    def __init__(
        self,
        bundle: InputBundle,
        query_points_m: Sequence[Sequence[float]],
        trace_sides: Sequence[TraceSide | None],
    ) -> None:
        if len(query_points_m) < 1 or len(query_points_m) > MAX_BATCH_QUERIES:
            raise ValueError(f"budgetRejected: query count must be in [1, {MAX_BATCH_QUERIES}]")
        if len(query_points_m) != len(trace_sides):
            raise ValueError("one optional trace side is required per query")
        self.bundle = bundle
        self.registry = build_registry(bundle.sources)
        self.sources = bundle.sources
        self.source_ids = tuple(source.source_id for source in self.sources)
        self.source_positions = np.asarray([source.position_m for source in self.sources], dtype=np.float64)
        self._source_row = {source.source_id: row for row, source in enumerate(self.sources)}
        self._source_by_id = {source.source_id: source for source in self.sources}
        self._enrichment_values: dict[str, dict[str, float]] = {}
        self._source_h = {source.source_id: (1.0 if source.position_m[1] > 0.0 else -1.0) for source in self.sources}
        for source in self.sources:
            fbar_i, _ = normalized_f1_and_gradient(source.position_m, None)
            self._enrichment_values[source.source_id] = {"F1": fbar_i}

        self.query_points_m = tuple(tuple(float(value) for value in point) for point in query_points_m)
        self.trace_sides = tuple(trace_sides)
        if any(len(point) != 3 for point in self.query_points_m):
            raise ValueError("queries must each contain three coordinates")
        self._query_key_to_index: dict[tuple[tuple[float, float, float], str | None], int] = {}
        for index, (point, side) in enumerate(zip(self.query_points_m, self.trace_sides)):
            key = (point, None if side is None else side.value)
            if key in self._query_key_to_index:
                raise ValueError("query point/trace identities must be unique in a frozen batch")
            self._query_key_to_index[key] = index
        validate_problem_budget(len(self.sources), self.registry.dof_count, len(self.query_points_m))

        self.classifications = tuple(
            classify_query(point, side, bundle.geometry)
            for point, side in zip(self.query_points_m, self.trace_sides)
        )
        self.values = np.full((len(self.query_points_m), self.registry.dof_count), np.nan, dtype=np.float64)
        self.gradients = np.full((len(self.query_points_m), self.registry.dof_count, 3), np.nan, dtype=np.float64)
        self.row_statuses: list[ProviderStatus] = [ProviderStatus.UNKNOWN] * len(self.query_points_m)
        self.row_diagnostics: list[str | None] = [classification.reason for classification in self.classifications]

        valid_indices = [
            index for index, classification in enumerate(self.classifications)
            if classification.status is ProviderStatus.OK
        ]
        self.mls_result: ShapeFunctions3DResult | None = None
        self._mls_row_by_query: dict[int, int] = {}
        if valid_indices:
            mls_queries = np.asarray([self.query_points_m[index] for index in valid_indices], dtype=np.float64)
            self.mls_result = shape_functions_3d(self.source_positions, mls_queries, bundle.mls_config)
            self._mls_row_by_query = {query_index: row for row, query_index in enumerate(valid_indices)}
            for query_index in valid_indices:
                mls_row = self._mls_row_by_query[query_index]
                status = self.mls_result.statuses[mls_row]
                if status is not MLSStatus.OK:
                    self.row_statuses[query_index] = ProviderStatus.UNKNOWN
                    self.row_diagnostics[query_index] = f"MLS:{status.value}"
                    continue
                try:
                    self._build_row(query_index, mls_row)
                except (ValueError, FloatingPointError, ProviderContractError) as exc:
                    self.row_statuses[query_index] = ProviderStatus.UNKNOWN
                    self.row_diagnostics[query_index] = f"row_build:{type(exc).__name__}:{exc}"
                    continue
                self.row_statuses[query_index] = ProviderStatus.OK
                self.row_diagnostics[query_index] = None

        self.coefficient_ids = self.registry.canonical_ids
        self.ok_query_count = sum(status is ProviderStatus.OK for status in self.row_statuses)
        self.unknown_query_count = len(self.row_statuses) - self.ok_query_count

    def _build_row(self, query_index: int, mls_row: int) -> None:
        if self.mls_result is None:
            raise ProviderContractError("MLS row storage is unavailable")
        point = self.query_points_m[query_index]
        side = self.trace_sides[query_index]
        classification = self.classifications[query_index]
        assert classification.h_value is not None
        fbar_q, grad_fbar_q = normalized_f1_and_gradient(point, side)
        phi_all = self.mls_result.weights[mls_row]
        grad_all = self.mls_result.gradients[mls_row]
        for slot, key in enumerate(self.registry.ordered_keys):
            source_row = self._source_row[key.source_id]
            source = self._source_by_id[key.source_id]
            phi = float(phi_all[source_row])
            grad_phi = grad_all[source_row]
            if key.basis_kind is BasisKind.BASE and key.basis_index == 0:
                value = phi
                gradient = grad_phi
            elif key.basis_kind is BasisKind.HEAVISIDE and key.basis_index == 1:
                factor = classification.h_value - self._source_h[key.source_id]
                value = phi * factor
                gradient = factor * grad_phi
            elif key.basis_kind is BasisKind.FRONT_F1_SHIFTED and key.basis_index == 1:
                fbar_source = self._enrichment_values[key.source_id]["F1"]
                difference = fbar_q - fbar_source
                value = phi * difference
                gradient = difference * grad_phi + phi * np.asarray(grad_fbar_q, dtype=np.float64)
            else:
                raise ProviderContractError(f"unsupported registry basis: {key.basis_kind.value}")
            if not math.isfinite(float(value)) or not np.all(np.isfinite(gradient)):
                raise ProviderContractError("provider row contains non-finite values")
            self.values[query_index, slot] = float(value)
            self.gradients[query_index, slot, :] = np.asarray(gradient, dtype=np.float64)

    def query_index(self, point_m: Sequence[float], trace_side: TraceSide | None = None) -> int:
        key = (tuple(float(value) for value in point_m), None if trace_side is None else trace_side.value)
        if key not in self._query_key_to_index:
            raise KeyError("query point and trace side are not in the frozen diagnostic batch")
        return self._query_key_to_index[key]

    def row_arrays(self, query_index: int) -> tuple[np.ndarray, np.ndarray] | None:
        if type(query_index) is not int or not 0 <= query_index < len(self.query_points_m):
            raise IndexError("query index outside the frozen batch")
        if self.row_statuses[query_index] is not ProviderStatus.OK:
            return None
        return self.values[query_index], self.gradients[query_index]

    def sample_at(self, query_index: int) -> AntiplaneShapeSample:
        if type(query_index) is not int or not 0 <= query_index < len(self.query_points_m):
            raise IndexError("query index outside the frozen batch")
        status = self.row_statuses[query_index]
        if status is not ProviderStatus.OK:
            sample = AntiplaneShapeSample((), (), ProviderStatus.UNKNOWN, self.row_diagnostics[query_index] or "unknown query")
            return validate_antiplane_sample(sample, self.registry)
        contributions = tuple(
            ScalarShapeContribution(
                float(self.values[query_index, slot]),
                tuple(float(value) for value in self.gradients[query_index, slot, :]),
            )
            for slot in range(self.registry.dof_count)
        )
        sample = AntiplaneShapeSample(self.coefficient_ids, contributions, ProviderStatus.OK)
        return validate_antiplane_sample(sample, self.registry)

    def evaluate(
        self,
        material_position: tuple[float, float, float],
        *,
        trace_side: TraceSide | None = None,
    ) -> AntiplaneShapeSample:
        if type(material_position) is not tuple or len(material_position) != 3:
            raise ValueError("material_position must be a three-value tuple")
        return self.sample_at(self.query_index(material_position, trace_side))


def topology_negative_controls(bundle: InputBundle) -> list[dict[str, Any]]:
    """Small robustly separated tangent/end/outside checks on the pinned case."""
    h = float(bundle.case["support"]["radius_m"])
    tol = float(bundle.case["topology"]["tolerance_m"])
    controls = [
        ("support_tangent_to_finite_front", (-h, 0.0, 0.0), "support_tangent_to_finite_front"),
        ("support_near_tangent_to_finite_front", (-(h + 0.5 * tol), 0.0, 0.0), "support_tangent_to_finite_front"),
        ("support_tangent_to_crack", (-0.005, h, 0.0), "support_tangent_to_finite_crack"),
        ("support_center_outside_domain", (0.0101, 0.004, 0.0), "source_center_outside_domain"),
    ]
    result: list[dict[str, Any]] = []
    for label, point, expected_reason in controls:
        record = classify_support(f"diagnostic-{label}", point, bundle.geometry)
        reason = record.get("reason_code")
        status = record.get("topology_status")
        result.append(
            {
                "case": label,
                "source_position_m": list(point),
                "expected_status": "UNKNOWN",
                "expected_reason": expected_reason,
                "actual_status": status,
                "actual_reason": reason,
                "passed": status == "UNKNOWN" and reason == expected_reason,
            }
        )
    return result


def topology_box_clipped_endpoint_control(bundle: InputBundle) -> dict[str, Any]:
    """Recheck a manifest H source whose clipped disk reaches the closed z end."""
    candidates = [
        source for source in bundle.sources
        if source.geometric_role == "HEAVISIDE_ELIGIBLE"
        and abs(source.position_m[2] - 0.008) <= 1.0e-15
    ]
    if not candidates:
        raise ManifestBlockedError("pinned role manifest has no H source at a box-clipped z endpoint")
    source = candidates[0]
    certificate = classify_support(source.source_id, source.position_m, bundle.geometry)
    section = certificate.get("plane_section_bounds_m") or {}
    expected = (
        certificate.get("topology_status") == "CERTIFIED"
        and certificate.get("geometric_role") == "HEAVISIDE_ELIGIBLE"
        and abs(float(section.get("z_max_m", math.nan)) - 0.008) <= 1.0e-15
        and float(certificate.get("front_clearance_m", -math.inf)) > bundle.geometry.tolerance_m
    )
    return {
        "source_id": source.source_id,
        "source_position_m": list(source.position_m),
        "actual_status": certificate.get("topology_status"),
        "actual_role": certificate.get("geometric_role"),
        "reason": certificate.get("reason_code"),
        "box_clipped_plane_section_bounds_m": section,
        "expected_plane_z_max_m": 0.008,
        "passed": expected,
        "scope_note": "This checks one endpoint source from the exact pinned case; it does not certify arbitrary clipped finite-crack geometries.",
    }

