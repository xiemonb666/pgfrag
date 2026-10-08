"""Source-only continuous-MLS provider skeleton for two pinned front states.

This module reclassifies the same material-rest source coordinates for each
allowed front position. It never imports the frozen local provider or consumes
its role assignments or coefficient values. No operator, solver, or coefficient
transfer is implemented here.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import copy
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
from types import ModuleType
from typing import Any, Mapping, Sequence

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
CASES_PATH = HERE / "moving_front_cases.json"
SUPPORT_ROOT = REPO_ROOT / "pgfrag/experiments/local_front_support_safe_v2_20261006"
REGISTRY_ROOT = REPO_ROOT / "pgfrag/experiments/enriched_dof_contract_20261005"
MLS_ROOT = REPO_ROOT / "pgfrag/experiments/mls3d_kernel_20261005"
SOURCE_CASE_PATH = REPO_ROOT / "pgfrag/experiments/local_front_support_contract_20261006/fixed_case.json"
COORDINATE_MANIFEST_PATH = REPO_ROOT / "pgfrag/experiments/local_front_support_contract_20261006/role_manifest.json"
MATERIAL_CASE_PATH = REPO_ROOT / "pgfrag/experiments/modeiii_pu_operator_20261005/fixed_case.json"

for _dependency_path in (str(REGISTRY_ROOT), str(MLS_ROOT)):
    if _dependency_path not in sys.path:
        sys.path.insert(0, _dependency_path)

from enriched_dof_contract import (  # noqa: E402
    AntiplaneShapeSample,
    BasisKind,
    DofKey,
    DofRegistry,
    FieldComponent,
    ProviderContractError,
    ProviderStatus,
    ScalarShapeContribution,
    validate_antiplane_sample,
    validate_problem_budget,
)
from mls3d_mls import MLS3DConfig, MLSStatus  # noqa: E402
from mls3d_shape import ShapeFunctions3DResult, shape_functions_3d  # noqa: E402


CASES_JSON_SHA256 = "2766245b25bf93eebcf6836befe81536917888936b4543666c86ebe9d6b72100"
SAFE_CLASSIFIER_REL = "pgfrag/experiments/local_front_support_safe_v2_20261006/support_topology.py"
SAFE_CLASSIFIER_SHA256 = "919cb60a6d65ae6b55a0bd12b569ec6627066b2a72dfcca52145209336588664"
SOURCE_CASE_REL = "pgfrag/experiments/local_front_support_contract_20261006/fixed_case.json"
SOURCE_CASE_SHA256 = "0bdbb20032e105851d2909ff54850ebc54d3f2041c25ab5e42e9102244d90ca5"
COORDINATE_MANIFEST_REL = "pgfrag/experiments/local_front_support_contract_20261006/role_manifest.json"
COORDINATE_MANIFEST_SHA256 = "e5e318f80b5dfedbc8ba122f82df3cca1585e604272399d21bc35907d1d28f01"
COORDINATE_RECORDS_SHA256 = "a8687a23353b6465aa0cc3fcda49a57c240297c38c141a017a167f3c526b23bb"
REGISTRY_REL = "pgfrag/experiments/enriched_dof_contract_20261005/enriched_dof_contract.py"
REGISTRY_SHA256 = "e3ddb3b067465bb5b269576a1b8a5100e2bb92c9c924bc09c949702c19efa9d8"
MLS_CORE_REL = "pgfrag/experiments/mls3d_kernel_20261005/mls3d_mls.py"
MLS_CORE_SHA256 = "404256d9b7546361a5704c4c3cec7c022069a9cf83a137c58a16e7a64f79ec63"
MLS_SHAPE_REL = "pgfrag/experiments/mls3d_kernel_20261005/mls3d_shape.py"
MLS_SHAPE_SHA256 = "66e680616a988182875d17aeb84fadb541867aee13fc15e0ae098f7c98839ef5"
MATERIAL_CASE_REL = "pgfrag/experiments/modeiii_pu_operator_20261005/fixed_case.json"
MATERIAL_CASE_SHA256 = "dce0a1f3009246a52198bc0ed833ed67998144549a653e7bce9bcc17babd71b9"
EXPECTED_FRONT_X_M = (0.0, 0.0005)
R0_M = 0.01
H_M = 0.007
MAX_QUERY_POINTS = 256


class MovingFrontBlockedError(ValueError):
    """The current geometry, role classification, or MLS batch is unsupported."""


class StaleBasisProvenanceError(ValueError):
    """A coefficient bundle was created for another surface or basis generation."""


class TraceSide(str, Enum):
    PLUS = "PLUS"
    MINUS = "MINUS"


@dataclass(frozen=True)
class RestSource:
    source_id: str
    lattice_index: tuple[int, int, int]
    position_m: tuple[float, float, float]


@dataclass(frozen=True)
class CurrentSourceRole:
    source_id: str
    position_m: tuple[float, float, float]
    geometric_role: str
    basis_role: str
    basis_roles: tuple[str, ...]
    certificate_reason: str
    shifted_fbar_source: float | None
    source_heaviside: float | None


@dataclass(frozen=True)
class CoefficientProvenance:
    """Tags required for coefficient use; it carries no coefficient values."""

    case_id: str
    front_x_m: float
    surface_geometry_sha256: str
    basis_generation_sha256: str
    canonical_ids_sha256: str
    source_solution_id: str


@dataclass(frozen=True)
class SampleProvenance:
    case_id: str
    front_x_m: float
    surface_geometry_sha256: str
    basis_generation_sha256: str
    source_coordinate_manifest_sha256: str
    safe_classifier_sha256: str
    canonical_ids_sha256: str


@dataclass(frozen=True)
class MovingFrontSample:
    provenance: SampleProvenance
    row: AntiplaneShapeSample


@dataclass(frozen=True)
class InputBundle:
    case_pair: Mapping[str, Any]
    source_case: Mapping[str, Any]
    coordinate_manifest: Mapping[str, Any]
    rest_sources: tuple[RestSource, ...]
    material_case: Mapping[str, Any]
    raw_pins: Mapping[str, str]


@dataclass(frozen=True)
class FrontState:
    case_id: str
    front_x_m: float
    geometry: Any
    case_geometry: Mapping[str, Any]
    sources: tuple[CurrentSourceRole, ...]
    registry: DofRegistry
    source_positions: np.ndarray
    source_row_by_id: Mapping[str, int]
    surface_geometry_sha256: str
    basis_generation_sha256: str
    source_coordinate_manifest_sha256: str
    safe_classifier_sha256: str
    mls_config: MLS3DConfig
    role_counts: Mapping[str, int]

    @property
    def canonical_ids_sha256(self) -> str:
        return sha256_bytes(_canonical_json(list(self.registry.canonical_ids)))


def _canonical_json(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False).encode("utf-8")


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _sha256_path(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def _check_pin(relative_path: str, expected_sha256: str) -> str:
    path = REPO_ROOT / relative_path
    if not path.is_file():
        raise MovingFrontBlockedError(f"required pinned source is missing: {relative_path}")
    actual = _sha256_path(path)
    if actual != expected_sha256:
        raise MovingFrontBlockedError(f"pinned source changed: {relative_path}")
    return actual


def _verify_loaded_module_origin(module_name: str, relative_path: str) -> None:
    """Fail closed if Python resolved a pinned dependency through another path."""
    module = sys.modules.get(module_name)
    expected_path = (REPO_ROOT / relative_path).resolve()
    actual_name = getattr(module, "__file__", None)
    if actual_name is None or Path(actual_name).resolve() != expected_path:
        raise MovingFrontBlockedError(f"loaded module {module_name} did not resolve to {relative_path}")


def _read_json_pinned(path: Path) -> tuple[dict[str, Any], bytes]:
    raw = path.read_bytes()
    try:
        parsed = json.loads(raw.decode("utf-8-sig"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise MovingFrontBlockedError(f"invalid frozen JSON {path.name}: {exc}") from exc
    if not isinstance(parsed, dict):
        raise MovingFrontBlockedError(f"expected an object in {path.name}")
    return parsed, raw


def load_safe_v2_classifier() -> ModuleType:
    """Load safe-v2 under a hash-specific module name, never as support_topology."""
    source_path = REPO_ROOT / SAFE_CLASSIFIER_REL
    actual = _check_pin(SAFE_CLASSIFIER_REL, SAFE_CLASSIFIER_SHA256)
    module_name = f"_pgfrag_local_front_support_safe_v2_{actual[:20]}"
    cached = sys.modules.get(module_name)
    if cached is not None:
        return cached
    spec = importlib.util.spec_from_file_location(module_name, source_path)
    if spec is None or spec.loader is None:
        raise MovingFrontBlockedError("could not create an isolated safe-v2 classifier import")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    if Path(module.__file__).resolve() != source_path.resolve():
        raise MovingFrontBlockedError("isolated classifier import resolved to an unexpected path")
    return module


def _expected_lattice_rows(case: Mapping[str, Any]) -> list[dict[str, Any]]:
    axes = case["sources"]["axes_m"]
    rows: list[dict[str, Any]] = []
    for ix, x in enumerate(axes["x"]):
        for iy, y in enumerate(axes["y"]):
            for iz, z in enumerate(axes["z"]):
                rows.append(
                    {
                        "source_id": f"lf-x{ix:02d}-y{iy:02d}-z{iz:02d}",
                        "lattice_index": [ix, iy, iz],
                        "position_m": [float(x), float(y), float(z)],
                    }
                )
    return rows


def load_pinned_inputs() -> InputBundle:
    """Read the old case only for rest coordinates; discard its role/cert data."""
    raw_pins = {
        "moving_front_cases.json": _check_pin(
            "pgfrag/experiments/moving_front_provider_20261006/moving_front_cases.json",
            CASES_JSON_SHA256,
        ),
        SOURCE_CASE_REL: _check_pin(SOURCE_CASE_REL, SOURCE_CASE_SHA256),
        COORDINATE_MANIFEST_REL: _check_pin(COORDINATE_MANIFEST_REL, COORDINATE_MANIFEST_SHA256),
        SAFE_CLASSIFIER_REL: _check_pin(SAFE_CLASSIFIER_REL, SAFE_CLASSIFIER_SHA256),
        REGISTRY_REL: _check_pin(REGISTRY_REL, REGISTRY_SHA256),
        MLS_CORE_REL: _check_pin(MLS_CORE_REL, MLS_CORE_SHA256),
        MLS_SHAPE_REL: _check_pin(MLS_SHAPE_REL, MLS_SHAPE_SHA256),
        MATERIAL_CASE_REL: _check_pin(MATERIAL_CASE_REL, MATERIAL_CASE_SHA256),
    }
    for module_name, relative_path in (
        ("enriched_dof_contract", REGISTRY_REL),
        ("mls3d_mls", MLS_CORE_REL),
        ("mls3d_shape", MLS_SHAPE_REL),
    ):
        _verify_loaded_module_origin(module_name, relative_path)
    case_pair, _case_pair_raw = _read_json_pinned(CASES_PATH)
    source_case, _source_case_raw = _read_json_pinned(SOURCE_CASE_PATH)
    coordinate_manifest, _manifest_raw = _read_json_pinned(COORDINATE_MANIFEST_PATH)
    material_case, _material_raw = _read_json_pinned(MATERIAL_CASE_PATH)

    fixed = case_pair["fixed_rest_domain_m"]
    if fixed != {
        "domain_box": [[-0.01, 0.01], [-0.01, 0.01], [-0.008, 0.008]],
        "crack_x_min_m": -0.01,
        "crack_y_m": 0.0,
        "crack_z_bounds_m": [-0.008, 0.008],
        "front_z_bounds_m": [-0.008, 0.008],
        "support_radius_h_m": H_M,
        "front_reference_radius_R0_m": R0_M,
        "topology_tolerance_m": 1.0e-12,
    }:
        raise MovingFrontBlockedError("case pair differs from the frozen box, rest crack, h, R0, or tolerance")
    if source_case.get("case_id") != "modeiii_local_front_pu_support_roles_20261006":
        raise MovingFrontBlockedError("unexpected material-rest source case")
    if coordinate_manifest.get("source_records_sha256") != COORDINATE_RECORDS_SHA256:
        raise MovingFrontBlockedError("coordinate source record digest changed")
    if coordinate_manifest.get("source_count") != 100 or len(coordinate_manifest.get("sources", [])) != 100:
        raise MovingFrontBlockedError("the material-rest coordinate set must contain 100 sources")
    if coordinate_manifest.get("manifest_sha256") != "b98799bfd6121370b7a8569c9bfc2e7e5958e8666687a797b0588b6412d10b35":
        raise MovingFrontBlockedError("pinned coordinate manifest payload digest changed")
    if sha256_bytes(_canonical_json(coordinate_manifest["sources"])) != COORDINATE_RECORDS_SHA256:
        raise MovingFrontBlockedError("coordinate records do not match their frozen digest")
    if sha256_bytes(
        _canonical_json({key: value for key, value in coordinate_manifest.items() if key != "manifest_sha256"})
    ) != coordinate_manifest["manifest_sha256"]:
        raise MovingFrontBlockedError("coordinate manifest payload digest is invalid")

    expected = _expected_lattice_rows(source_case)
    actual_records = coordinate_manifest["sources"]
    rest_sources: list[RestSource] = []
    if len(expected) != len(actual_records):
        raise MovingFrontBlockedError("coordinate manifest row count differs from the frozen lattice")
    for expected_row, record in zip(expected, actual_records):
        # Deliberately access only source_id/lattice_index/position_m. Old
        # geometric_role, basis_role(s), and certificates are never consulted.
        if (
            record.get("source_id") != expected_row["source_id"]
            or record.get("lattice_index") != expected_row["lattice_index"]
            or record.get("position_m") != expected_row["position_m"]
        ):
            raise MovingFrontBlockedError("rest source ID/coordinate binding changed")
        position = tuple(float(value) for value in expected_row["position_m"])
        if position[1] == 0.0:
            raise MovingFrontBlockedError("all fixed shifted-enrichment source y coordinates must be nonzero")
        rest_sources.append(
            RestSource(
                source_id=expected_row["source_id"],
                lattice_index=tuple(expected_row["lattice_index"]),
                position_m=position,
            )
        )

    material = material_case.get("material_and_manufactured_field", {})
    if (
        material.get("young_modulus_Pa") != 1_000_000.0
        or material.get("poisson_ratio") != 0.3
        or material.get("K_III_Pa_sqrt_m") != 100.0
        or material_case.get("geometry", {}).get("radius_m") != R0_M
    ):
        raise MovingFrontBlockedError("mathematical Mode III comparison constants differ from their pinned fixture")
    return InputBundle(
        case_pair=case_pair,
        source_case=source_case,
        coordinate_manifest=coordinate_manifest,
        rest_sources=tuple(rest_sources),
        material_case=material_case,
        raw_pins=raw_pins,
    )


def _role_to_keys(geometric_role: str, basis_role: str) -> tuple[str, ...]:
    if geometric_role == "BASE" and basis_role == "BASE":
        return ("BASE",)
    if geometric_role == "FRONT" and basis_role == "FRONT_F1_SHIFTED":
        return ("BASE", "FRONT_F1_SHIFTED")
    if geometric_role == "HEAVISIDE_ELIGIBLE" and basis_role == "HEAVISIDE":
        return ("BASE", "HEAVISIDE")
    raise MovingFrontBlockedError(f"unknown or inconsistent current role {geometric_role}/{basis_role}")


def _sha_surface_geometry(case: Mapping[str, Any]) -> str:
    geometry = case["geometry"]
    crack = geometry["crack"]
    front = geometry["front"]
    payload = {
        "domain_box_m": geometry["domain_box_m"],
        "finite_Gamma": {
            "x_min_m": crack["x_min_m"],
            "front_x_m": crack["front_x_m"],
            "y_plane_m": crack["y_plane_m"],
            "z_bounds_m": crack["z_bounds_m"],
            "x_start_closed": True,
            "front_excluded_from_Gamma": True,
            "z_ends_closed": True,
        },
        "finite_front_segment": {
            "x_m": front["x_m"],
            "y_m": front["y_m"],
            "z_bounds_m": front["z_bounds_m"],
            "segment_ends_closed": True,
        },
    }
    return sha256_bytes(_canonical_json(payload))


def _front_relative_fbar(
    point_m: Sequence[float],
    front_x_m: float,
    y_plane_m: float,
    trace_side: TraceSide | None,
) -> tuple[float, tuple[float, float, float]]:
    xi = float(point_m[0]) - front_x_m
    eta = float(point_m[1]) - y_plane_m
    radius = math.hypot(xi, eta)
    if radius == 0.0:
        raise ValueError("shifted F1 is undefined on the current finite front")
    if eta == 0.0 and xi < 0.0:
        if trace_side is TraceSide.PLUS:
            theta = math.pi
        elif trace_side is TraceSide.MINUS:
            theta = -math.pi
        else:
            raise ValueError("exact current-crack trace requires an explicit PLUS/MINUS side")
    else:
        theta = math.atan2(eta, xi)
    half = 0.5 * theta
    sine = math.sin(half)
    cosine = math.cos(half)
    value = math.sqrt(radius / R0_M) * sine
    denominator = 2.0 * math.sqrt(R0_M * radius)
    return value, (-sine / denominator, cosine / denominator, 0.0)


def _basis_generation_hash(
    surface_geometry_sha256: str,
    roles: Sequence[CurrentSourceRole],
    registry: DofRegistry,
    classifier_sha256: str,
) -> str:
    payload = {
        "surface_geometry_sha256": surface_geometry_sha256,
        "current_source_roles": [
            {
                "source_id": row.source_id,
                "position_m": list(row.position_m),
                "geometric_role": row.geometric_role,
                "basis_role": row.basis_role,
                "basis_roles": list(row.basis_roles),
                "current_shifted_Fbar_source": row.shifted_fbar_source,
                "source_Heaviside": row.source_heaviside,
            }
            for row in roles
        ],
        "canonical_ids": list(registry.canonical_ids),
        "h_m": H_M,
        "R0_m": R0_M,
        "safe_classifier_sha256": classifier_sha256,
    }
    return sha256_bytes(_canonical_json(payload))


def _new_registry(current_sources: Sequence[CurrentSourceRole]) -> DofRegistry:
    registry = DofRegistry()
    for source in current_sources:
        for role in source.basis_roles:
            if role == "BASE":
                basis_kind, basis_index = BasisKind.BASE, 0
            elif role == "HEAVISIDE":
                basis_kind, basis_index = BasisKind.HEAVISIDE, 1
            elif role == "FRONT_F1_SHIFTED":
                basis_kind, basis_index = BasisKind.FRONT_F1_SHIFTED, 1
            else:
                raise MovingFrontBlockedError(f"unsupported current basis role {role!r}")
            registry.register(DofKey(source.source_id, FieldComponent.UZ, basis_kind, basis_index))
    registry.freeze()
    validate_problem_budget(len(current_sources), registry.dof_count, 0)
    if registry.source_count != len(current_sources):
        raise ProviderContractError("each material source must have one BASE key")
    if any(len(source.basis_roles) not in (1, 2) or source.basis_roles[0] != "BASE" for source in current_sources):
        raise ProviderContractError("each source may have BASE and at most one current enrichment")
    return registry


def build_front_state(bundle: InputBundle, front_x_m: float) -> FrontState:
    """Reclassify all rest sources and freeze a fresh role/DOF state for one front."""
    if type(front_x_m) not in (int, float) or float(front_x_m) not in EXPECTED_FRONT_X_M:
        raise MovingFrontBlockedError("only the two explicitly reviewed front positions are supported")
    case_record_by_x = {
        float(row["front_x_m"]): row for row in bundle.case_pair["moving_front_cases"]
    }
    front_x_m = float(front_x_m)
    if front_x_m not in case_record_by_x:
        raise MovingFrontBlockedError("front location has no pinned case record")
    case_record = case_record_by_x[front_x_m]
    safe = load_safe_v2_classifier()
    moving_case = copy.deepcopy(bundle.source_case)
    crack = moving_case["geometry"]["crack"]
    crack["front_x_m"] = front_x_m
    crack["definition"] = f"finite rectangle on y=0 from x=-0.01 m up to but excluding front x={front_x_m:g} m, with z in [-0.008,0.008] m"
    moving_case["geometry"]["front"]["x_m"] = front_x_m
    moving_case["geometry"]["front"]["definition"] = f"closed finite straight segment at x={front_x_m:g} m, y=0, z in [-0.008,0.008] m"
    expected_crack = {
        "x_interval_m": [-0.01, front_x_m],
        "x_front_is_excluded_from_Gamma": True,
        "y_m": 0.0,
        "z_interval_m": [-0.008, 0.008],
    }
    expected_front = {"x_m": front_x_m, "y_m": 0.0, "z_interval_m": [-0.008, 0.008]}
    if case_record.get("finite_crack") != expected_crack or case_record.get("front_segment") != expected_front:
        raise MovingFrontBlockedError("pinned front record does not match its finite Gamma/F geometry")
    geometry = safe.FiniteCrackGeometry.from_case(moving_case)

    current_sources: list[CurrentSourceRole] = []
    for source in bundle.rest_sources:
        certificate = safe.classify_support(source.source_id, source.position_m, geometry)
        if certificate.get("status") != "CERTIFIED":
            raise MovingFrontBlockedError(
                f"current geometry leaves source {source.source_id} UNKNOWN: {certificate.get('reason_code')}"
            )
        geometric_role = str(certificate["geometric_role"])
        basis_role = str(certificate["basis_role"])
        basis_roles = _role_to_keys(geometric_role, basis_role)
        fbar_source: float | None = None
        source_h: float | None = None
        if basis_role == "FRONT_F1_SHIFTED":
            fbar_source, _ = _front_relative_fbar(
                source.position_m, front_x_m, geometry.y_plane_m, None
            )
        elif basis_role == "HEAVISIDE":
            source_h = 1.0 if source.position_m[1] > geometry.y_plane_m else -1.0
        current_sources.append(
            CurrentSourceRole(
                source_id=source.source_id,
                position_m=source.position_m,
                geometric_role=geometric_role,
                basis_role=basis_role,
                basis_roles=basis_roles,
                certificate_reason=str(certificate["reason_code"]),
                shifted_fbar_source=fbar_source,
                source_heaviside=source_h,
            )
        )

    role_counts = {
        role: sum(source.geometric_role == role for source in current_sources)
        for role in ("FRONT", "HEAVISIDE_ELIGIBLE", "BASE")
    }
    registry = _new_registry(current_sources)
    # A read-only flag on an owning ndarray can be switched back on by a
    # caller. Back the public view with immutable bytes so WRITEABLE cannot be
    # re-enabled through setflags(write=True).
    position_bytes = np.asarray(
        [source.position_m for source in current_sources], dtype=np.float64
    ).tobytes(order="C")
    positions = np.frombuffer(position_bytes, dtype=np.float64).reshape((len(current_sources), 3))
    if positions.flags.writeable or positions.flags.owndata:
        raise ProviderContractError("source positions must be bytes-backed and immutable")
    source_row_by_id = {source.source_id: row for row, source in enumerate(current_sources)}
    surface_hash = _sha_surface_geometry(moving_case)
    basis_hash = _basis_generation_hash(surface_hash, current_sources, registry, SAFE_CLASSIFIER_SHA256)
    mls_data = bundle.material_case["mls"]
    config = MLS3DConfig(
        reference_origin=tuple(float(value) for value in mls_data["reference_origin_m"]),
        reference_scale=float(mls_data["reference_scale_m"]),
        support_radius=H_M,
        max_condition_number=float(mls_data["max_condition_number"]),
        rank_relative_tolerance=float(mls_data["rank_relative_tolerance"]),
        switch_tolerance=float(mls_data["switch_tolerance"]),
        duplicate_tolerance=float(mls_data["duplicate_tolerance"]),
    )
    if config.reference_origin != (0.0, 0.0, 0.0) or config.reference_scale != R0_M:
        raise MovingFrontBlockedError("continuous MLS reference frame changed")
    return FrontState(
        case_id=str(case_record["case_id"]),
        front_x_m=front_x_m,
        geometry=geometry,
        case_geometry=moving_case["geometry"],
        sources=tuple(current_sources),
        registry=registry,
        source_positions=positions,
        source_row_by_id=source_row_by_id,
        surface_geometry_sha256=surface_hash,
        basis_generation_sha256=basis_hash,
        source_coordinate_manifest_sha256=COORDINATE_RECORDS_SHA256,
        safe_classifier_sha256=SAFE_CLASSIFIER_SHA256,
        mls_config=config,
        role_counts=role_counts,
    )


def assert_distinct_front_generations(previous: FrontState, current: FrontState) -> None:
    """Reject identical geometry/basis hashes when front_x changes."""
    if previous.front_x_m == current.front_x_m:
        raise MovingFrontBlockedError("comparison requires two distinct front positions")
    if previous.surface_geometry_sha256 == current.surface_geometry_sha256:
        raise MovingFrontBlockedError("surface geometry hash failed to change with front geometry")
    if previous.basis_generation_sha256 == current.basis_generation_sha256:
        raise MovingFrontBlockedError("basis generation hash failed to change with front geometry")


def require_current_coefficient_provenance(state: FrontState, provenance: CoefficientProvenance) -> None:
    """Check identity only; coefficient mapping or mechanical reprojection is absent."""
    matches = (
        provenance.case_id == state.case_id
        and type(provenance.source_solution_id) is str
        and bool(provenance.source_solution_id.strip())
        and provenance.front_x_m == state.front_x_m
        and provenance.surface_geometry_sha256 == state.surface_geometry_sha256
        and provenance.basis_generation_sha256 == state.basis_generation_sha256
        and provenance.canonical_ids_sha256 == state.canonical_ids_sha256
    )
    if not matches:
        raise StaleBasisProvenanceError(
            "coefficient provenance belongs to another surface/basis generation; reproject or re-solve before use"
        )


def _classify_query(state: FrontState, point: Sequence[float], side: TraceSide | None) -> tuple[ProviderStatus, str, float | None]:
    if len(point) != 3 or any(not math.isfinite(float(value)) for value in point):
        return ProviderStatus.UNKNOWN, "invalid_or_nonfinite_query", None
    if side is not None and type(side) is not TraceSide:
        return ProviderStatus.UNKNOWN, "unsupported_trace_side", None
    geometry = state.geometry
    query = tuple(float(value) for value in point)
    if not geometry.box.contains(query, 0.0):
        return ProviderStatus.UNKNOWN, "outside_closed_domain", None
    x, y, z = query
    z_gap = max(geometry.crack_z_bounds_m[0] - z, 0.0, z - geometry.crack_z_bounds_m[1])
    distance_to_front = math.sqrt(
        (x - geometry.front_x_m) ** 2 + (y - geometry.y_plane_m) ** 2 + z_gap**2
    )
    if distance_to_front <= geometry.tolerance_m:
        return ProviderStatus.UNKNOWN, "front_or_finite_front_endpoint", None
    near_crack_plane = abs(y - geometry.y_plane_m) <= geometry.tolerance_m and x < geometry.front_x_m
    if near_crack_plane and y != geometry.y_plane_m:
        return ProviderStatus.UNKNOWN, "ambiguous_crack_plane_tolerance_band", None
    if y == geometry.y_plane_m and x < geometry.front_x_m:
        interior_x = geometry.crack_x_min_m + geometry.tolerance_m < x < geometry.front_x_m - geometry.tolerance_m
        interior_z = geometry.crack_z_bounds_m[0] + geometry.tolerance_m < z < geometry.crack_z_bounds_m[1] - geometry.tolerance_m
        if not (interior_x and interior_z):
            return ProviderStatus.UNKNOWN, "finite_crack_end_or_edge", None
        if side is None:
            return ProviderStatus.UNKNOWN, "explicit_plus_minus_trace_required", None
        return ProviderStatus.OK, "exact_finite_gamma_trace", 1.0 if side is TraceSide.PLUS else -1.0
    if side is not None:
        return ProviderStatus.UNKNOWN, "trace_side_only_valid_on_current_finite_gamma", None
    if y == geometry.y_plane_m and x >= geometry.front_x_m:
        # Ahead of the finite front is not classified as an infinite cut.
        h_value = 0.0
    else:
        h_value = 1.0 if y > geometry.y_plane_m else -1.0
    return ProviderStatus.OK, "open_domain_point", h_value


class MovingFrontProvider:
    """Dense BASE/H/shifted-F1 rows tagged with one moving surface generation."""

    def __init__(
        self,
        state: FrontState,
        query_points_m: Sequence[Sequence[float]],
        trace_sides: Sequence[TraceSide | None],
    ) -> None:
        if len(query_points_m) < 1 or len(query_points_m) > MAX_QUERY_POINTS:
            raise ValueError(f"budgetRejected: query batch must contain 1..{MAX_QUERY_POINTS} rows")
        if len(query_points_m) != len(trace_sides):
            raise ValueError("one optional trace side is required per query")
        validate_problem_budget(len(state.sources), state.registry.dof_count, len(query_points_m))
        self.state = state
        self.query_points_m = tuple(tuple(float(value) for value in point) for point in query_points_m)
        self.trace_sides = tuple(trace_sides)
        self._lookup: dict[tuple[tuple[float, float, float], str | None], int] = {}
        for index, (point, side) in enumerate(zip(self.query_points_m, self.trace_sides)):
            if len(point) != 3:
                raise ValueError("each query must contain exactly three coordinates")
            key = (point, None if side is None else side.value)
            if key in self._lookup:
                raise ValueError("query point/trace identities must be unique")
            self._lookup[key] = index
        self.classifications = tuple(
            _classify_query(state, point, side)
            for point, side in zip(self.query_points_m, self.trace_sides)
        )
        self.values = np.full((len(self.query_points_m), state.registry.dof_count), np.nan, dtype=np.float64)
        self.gradients = np.full((len(self.query_points_m), state.registry.dof_count, 3), np.nan, dtype=np.float64)
        self.row_statuses: list[ProviderStatus] = [ProviderStatus.UNKNOWN] * len(self.query_points_m)
        self.row_reasons = [classification[1] for classification in self.classifications]
        valid_indices = [index for index, result in enumerate(self.classifications) if result[0] is ProviderStatus.OK]
        self._mls_row_by_query: dict[int, int] = {}
        self.mls_result: ShapeFunctions3DResult | None = None
        if valid_indices:
            mls_queries = np.asarray([self.query_points_m[index] for index in valid_indices], dtype=np.float64)
            # Deliberately omit visibility_mask: continuous MLS, all 100 rest sources.
            self.mls_result = shape_functions_3d(state.source_positions, mls_queries, state.mls_config)
            self._mls_row_by_query = {query_index: row for row, query_index in enumerate(valid_indices)}
            for query_index in valid_indices:
                mls_row = self._mls_row_by_query[query_index]
                status = self.mls_result.statuses[mls_row]
                if status is not MLSStatus.OK:
                    self.row_reasons[query_index] = f"MLS:{status.value}"
                    continue
                try:
                    self._build_row(query_index, mls_row)
                except (ValueError, FloatingPointError, ProviderContractError) as exc:
                    self.row_reasons[query_index] = f"row_build:{type(exc).__name__}:{exc}"
                    continue
                self.row_statuses[query_index] = ProviderStatus.OK
                self.row_reasons[query_index] = "OK"
        self.provenance = SampleProvenance(
            case_id=state.case_id,
            front_x_m=state.front_x_m,
            surface_geometry_sha256=state.surface_geometry_sha256,
            basis_generation_sha256=state.basis_generation_sha256,
            source_coordinate_manifest_sha256=state.source_coordinate_manifest_sha256,
            safe_classifier_sha256=state.safe_classifier_sha256,
            canonical_ids_sha256=state.canonical_ids_sha256,
        )

    def _build_row(self, query_index: int, mls_row: int) -> None:
        if self.mls_result is None:
            raise ProviderContractError("MLS result is absent")
        point = self.query_points_m[query_index]
        side = self.trace_sides[query_index]
        h_value = self.classifications[query_index][2]
        if h_value is None:
            raise ProviderContractError("valid query has no Heaviside side value")
        fbar_q, grad_fbar_q = _front_relative_fbar(
            point, self.state.front_x_m, self.state.geometry.y_plane_m, side
        )
        phi_all = self.mls_result.weights[mls_row]
        grad_all = self.mls_result.gradients[mls_row]
        source_by_id = {source.source_id: source for source in self.state.sources}
        for slot, key in enumerate(self.state.registry.ordered_keys):
            source_row = self.state.source_row_by_id[key.source_id]
            source = source_by_id[key.source_id]
            phi = float(phi_all[source_row])
            grad_phi = grad_all[source_row]
            if key.basis_kind is BasisKind.BASE:
                value, gradient = phi, grad_phi
            elif key.basis_kind is BasisKind.HEAVISIDE:
                source_h = source.source_heaviside
                if source_h is None:
                    raise ProviderContractError("H source lacks its current shifted source-side value")
                factor = float(h_value) - source_h
                value, gradient = phi * factor, factor * grad_phi
            elif key.basis_kind is BasisKind.FRONT_F1_SHIFTED:
                source_fbar = source.shifted_fbar_source
                if source_fbar is None:
                    raise ProviderContractError("F1 source lacks its current shifted source value")
                difference = fbar_q - source_fbar
                value = phi * difference
                gradient = difference * grad_phi + phi * np.asarray(grad_fbar_q, dtype=np.float64)
            else:
                raise ProviderContractError(f"unsupported moving-front basis: {key.basis_kind.value}")
            if not math.isfinite(float(value)) or not np.all(np.isfinite(gradient)):
                raise ProviderContractError("moving-front row contains non-finite data")
            self.values[query_index, slot] = float(value)
            self.gradients[query_index, slot, :] = np.asarray(gradient, dtype=np.float64)

    def sample_at(self, query_index: int) -> MovingFrontSample:
        if type(query_index) is not int or not 0 <= query_index < len(self.query_points_m):
            raise IndexError("query index outside frozen batch")
        if self.row_statuses[query_index] is not ProviderStatus.OK:
            row = AntiplaneShapeSample((), (), ProviderStatus.UNKNOWN, self.row_reasons[query_index])
        else:
            contributions = tuple(
                ScalarShapeContribution(
                    float(self.values[query_index, slot]),
                    tuple(float(value) for value in self.gradients[query_index, slot, :]),
                )
                for slot in range(self.state.registry.dof_count)
            )
            row = AntiplaneShapeSample(
                self.state.registry.canonical_ids,
                contributions,
                ProviderStatus.OK,
            )
        return MovingFrontSample(
            provenance=self.provenance,
            row=validate_antiplane_sample(row, self.state.registry),
        )

    def evaluate(self, point_m: tuple[float, float, float], *, trace_side: TraceSide | None = None) -> MovingFrontSample:
        if type(point_m) is not tuple or len(point_m) != 3:
            raise ValueError("query point must be a three-value tuple")
        key = (tuple(float(value) for value in point_m), None if trace_side is None else trace_side.value)
        if key not in self._lookup:
            raise KeyError("query and trace identity are not part of the frozen batch")
        return self.sample_at(self._lookup[key])


def analytic_modeiii_reference(
    point_m: Sequence[float],
    front_x_m: float,
    young_modulus_Pa: float,
    poisson_ratio: float,
    K_III_Pa_sqrt_m: float,
    trace_side: TraceSide | None = None,
) -> tuple[float, tuple[float, float, float]]:
    """Mathematical current-front manufactured reference, not a solver state."""
    mu = young_modulus_Pa / (2.0 * (1.0 + poisson_ratio))
    amplitude = (K_III_Pa_sqrt_m / mu) * math.sqrt(2.0 / math.pi)
    fbar, grad_fbar = _front_relative_fbar(point_m, front_x_m, 0.0, trace_side)
    coefficient = amplitude * math.sqrt(R0_M)
    return coefficient * fbar, tuple(coefficient * value for value in grad_fbar)

