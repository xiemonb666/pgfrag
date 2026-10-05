"""Safe-v2 finite-crack support topology classifier.

This file is copied from
``pgfrag/experiments/local_front_support_contract_20261006/support_topology.py``
(SHA-256 ea20e51d2b774c7e6b8b48403359b56bab1ba6f5405dba41f34813b6aec4578e)
with the geometry/tolerance fixes documented alongside this copy. It is
source-only geometry code: it does not import MLS, allocate provider rows,
assemble an operator, or run a solver.

Set convention: Omega is a closed axis-aligned box and each support is Omega
intersected with a closed Euclidean ball. The finite crack is
``x in [x_min, front), y = y_plane, z in [z_min, z_max]``; its x start and z
ends are included, while its front tip is excluded. The finite front is the
closed segment ``x = front, y = y_plane, z in [z_min, z_max]``. H is certified
only when the entire clipped plane section is inside the finite crack and the
front has positive clearance. Exact contact with a closed box/clipped endpoint
is allowed. A nonzero endpoint gap within the topology tolerance is UNKNOWN;
no geometry is shifted by an epsilon.
"""
from __future__ import annotations

from dataclasses import dataclass
import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

SCHEMA = "local-front-support-role-manifest.v1"


class GeometryInputError(ValueError):
    """Malformed geometry or source input."""


class ManifestBlockedError(ValueError):
    """The selected source set contains UNKNOWN or exceeds a frozen count."""


@dataclass(frozen=True)
class Box3D:
    """Closed axis-aligned physical domain, bounds ordered x,y,z in metres."""

    bounds_m: tuple[tuple[float, float], tuple[float, float], tuple[float, float]]

    def __post_init__(self) -> None:
        if len(self.bounds_m) != 3:
            raise GeometryInputError("box must have x, y, and z intervals")
        for interval in self.bounds_m:
            if len(interval) != 2 or not all(_finite_real(v) for v in interval):
                raise GeometryInputError("box intervals must contain two finite metre values")
            if interval[0] >= interval[1]:
                raise GeometryInputError("box interval lower bound must be below upper bound")

    def contains(self, point_m: Sequence[float], tolerance_m: float) -> bool:
        return all(
            lower - tolerance_m <= value <= upper + tolerance_m
            for value, (lower, upper) in zip(point_m, self.bounds_m)
        )


@dataclass(frozen=True)
class FiniteCrackGeometry:
    """Finite planar crack x in [x_min,x_front), y=y_plane, z in [z_min,z_max]."""

    box: Box3D
    crack_x_min_m: float
    front_x_m: float
    y_plane_m: float
    crack_z_bounds_m: tuple[float, float]
    support_radius_m: float
    tolerance_m: float

    def __post_init__(self) -> None:
        values = (
            self.crack_x_min_m,
            self.front_x_m,
            self.y_plane_m,
            self.support_radius_m,
            self.tolerance_m,
            *self.crack_z_bounds_m,
        )
        if not all(_finite_real(value) for value in values):
            raise GeometryInputError("crack geometry values must be finite metre values")
        if self.crack_x_min_m >= self.front_x_m:
            raise GeometryInputError("crack x extent must have positive length")
        if self.crack_z_bounds_m[0] >= self.crack_z_bounds_m[1]:
            raise GeometryInputError("crack z extent must have positive length")
        if self.support_radius_m <= 0.0 or self.tolerance_m <= 0.0:
            raise GeometryInputError("support radius and topology tolerance must be positive")
        bx, by, bz = self.box.bounds_m
        if self.crack_x_min_m < bx[0] or self.front_x_m > bx[1]:
            raise GeometryInputError("finite crack x extent must lie inside the physical box")
        if not (by[0] + self.tolerance_m < self.y_plane_m < by[1] - self.tolerance_m):
            raise GeometryInputError("crack plane must be strictly inside the box in y")
        if self.crack_z_bounds_m[0] < bz[0] or self.crack_z_bounds_m[1] > bz[1]:
            raise GeometryInputError("finite crack z extent must lie inside the physical box")

    @classmethod
    def from_case(cls, case: Mapping[str, Any]) -> "FiniteCrackGeometry":
        try:
            domain = case["geometry"]["domain_box_m"]
            crack = case["geometry"]["crack"]
            h = case["support"]["radius_m"]
            tolerance = case["topology"]["tolerance_m"]
            box = Box3D(tuple(tuple(float(v) for v in interval) for interval in domain))
            return cls(
                box=box,
                crack_x_min_m=float(crack["x_min_m"]),
                front_x_m=float(crack["front_x_m"]),
                y_plane_m=float(crack["y_plane_m"]),
                crack_z_bounds_m=tuple(float(v) for v in crack["z_bounds_m"]),
                support_radius_m=float(h),
                tolerance_m=float(tolerance),
            )
        except (KeyError, TypeError, ValueError) as exc:
            raise GeometryInputError(f"invalid case geometry: {exc}") from exc


def _finite_real(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value))


def _point3(point_m: Sequence[float]) -> tuple[float, float, float]:
    if isinstance(point_m, (str, bytes)) or len(point_m) != 3:
        raise GeometryInputError("source position must contain exactly three coordinates in metres")
    if any(not _finite_real(value) for value in point_m):
        raise GeometryInputError("source coordinates must be finite real values")
    return tuple(float(value) for value in point_m)  # type: ignore[return-value]


def _interval_distance(value: float, interval: tuple[float, float]) -> float:
    return max(interval[0] - value, 0.0, value - interval[1])


def _distance_to_front(point: tuple[float, float, float], geometry: FiniteCrackGeometry) -> float:
    x, y, z = point
    z_gap = _interval_distance(z, geometry.crack_z_bounds_m)
    return math.sqrt((x - geometry.front_x_m) ** 2 + (y - geometry.y_plane_m) ** 2 + z_gap**2)


def _distance_to_finite_crack(point: tuple[float, float, float], geometry: FiniteCrackGeometry) -> float:
    x, y, z = point
    dx = _interval_distance(x, (geometry.crack_x_min_m, geometry.front_x_m))
    dz = _interval_distance(z, geometry.crack_z_bounds_m)
    return math.sqrt(dx**2 + (y - geometry.y_plane_m) ** 2 + dz**2)


def _point_on_crack(point: tuple[float, float, float], geometry: FiniteCrackGeometry) -> bool:
    x, y, z = point
    tol = geometry.tolerance_m
    return (
        abs(y - geometry.y_plane_m) <= tol
        and geometry.crack_x_min_m - tol <= x < geometry.front_x_m - tol
        and geometry.crack_z_bounds_m[0] - tol <= z <= geometry.crack_z_bounds_m[1] + tol
    )


def _point_on_front(point: tuple[float, float, float], geometry: FiniteCrackGeometry) -> bool:
    x, y, z = point
    tol = geometry.tolerance_m
    return (
        abs(x - geometry.front_x_m) <= tol
        and abs(y - geometry.y_plane_m) <= tol
        and geometry.crack_z_bounds_m[0] - tol <= z <= geometry.crack_z_bounds_m[1] + tol
    )


def _plane_section_bounds(
    point: tuple[float, float, float], geometry: FiniteCrackGeometry
) -> dict[str, float] | None:
    """Coordinate bounds of (Omega intersect ball) intersect y=y_plane.

    Since the source center is required inside the convex box, each individual
    x/z extremum is attained at the corresponding axis point of the disk before
    clipping by the box. These bounds are used only as a containment certificate
    for the finite crack rectangle, never to infer an infinite-plane cut.
    """
    x, y, z = point
    h = geometry.support_radius_m
    radicand = h * h - (y - geometry.y_plane_m) ** 2
    if radicand < 0.0:
        return None
    rho = math.sqrt(max(0.0, radicand))
    (x_box, _y_box, z_box) = geometry.box.bounds_m
    return {
        "x_min_m": max(x_box[0], x - rho),
        "x_max_m": min(x_box[1], x + rho),
        "z_min_m": max(z_box[0], z - rho),
        "z_max_m": min(z_box[1], z + rho),
        "disk_radius_m": rho,
    }


def _certificate(
    source_id: str,
    point: tuple[float, float, float],
    geometry: FiniteCrackGeometry,
    *,
    geometric_role: str,
    basis_role: str,
    reason_code: str,
    front_distance: float,
    crack_distance: float,
    plane_bounds: dict[str, float] | None = None,
    split_component_count: int | None = None,
    both_sides_margin: float | None = None,
) -> dict[str, Any]:
    h = geometry.support_radius_m
    return {
        "source_id": source_id,
        "position_m": list(point),
        "support_radius_m": h,
        "geometric_role": geometric_role,
        "basis_role": basis_role,
        "topology_status": "UNKNOWN" if geometric_role == "UNKNOWN" else "CERTIFIED",
        "basis_roles": ["BASE"] + ([] if basis_role in ("BASE", "UNKNOWN") else [basis_role]),
        "status": "UNKNOWN" if basis_role == "UNKNOWN" else "CERTIFIED",
        "reason_code": reason_code,
        "distance_to_finite_front_m": front_distance,
        "front_clearance_m": front_distance - h,
        "front_intersection_margin_m": h - front_distance,
        "distance_to_finite_crack_m": crack_distance,
        "crack_clearance_m": crack_distance - h,
        "plane_section_y_equals_crack_m": geometry.y_plane_m,
        "plane_section_bounds_m": plane_bounds,
        "plane_section_inside_finite_crack": plane_bounds is not None and geometric_role == "HEAVISIDE_ELIGIBLE",
        "positive_side_margin_m": both_sides_margin,
        "connected_components_after_crack_removal": split_component_count,
        "convexity_argument": (
            "Omega intersect closed_ball is convex; after removing the entire plane section contained in finite Gamma, its y>plane and y<plane pieces are convex and exactly two"
            if split_component_count == 2
            else None
        ),
    }


def classify_support(
    source_id: str,
    source_position_m: Sequence[float],
    geometry: FiniteCrackGeometry,
) -> dict[str, Any]:
    """Return a finite-geometry role certificate or fail-closed UNKNOWN.

    FRONT takes priority over HEAVISIDE. H is certified only if the complete
    plane section of the box-clipped support lies inside the finite crack patch,
    both support sides are nonempty, and the support has positive clearance from
    the finite front. Crack/front tangencies and partial finite-patch cuts are
    UNKNOWN. The H and F1 bases here are shifted, so sources on their branch set
    cannot be assigned a made-up source-side value.
    """
    if type(source_id) is not str or not source_id or source_id.strip() != source_id:
        raise GeometryInputError("source_id must be a nonempty canonical string")
    point = _point3(source_position_m)
    h = geometry.support_radius_m
    tol = geometry.tolerance_m
    # Source centers are domain data, not approximate geometric predicates.
    # A center outside the closed physical box is out of domain at any distance.
    if not geometry.box.contains(point, 0.0):
        df = _distance_to_front(point, geometry)
        dg = _distance_to_finite_crack(point, geometry)
        return _certificate(source_id, point, geometry, geometric_role="UNKNOWN", basis_role="UNKNOWN",
                            reason_code="source_center_outside_domain", front_distance=df, crack_distance=dg)

    d_front = _distance_to_front(point, geometry)
    d_crack = _distance_to_finite_crack(point, geometry)

    if abs(d_front - h) <= tol:
        return _certificate(source_id, point, geometry, geometric_role="UNKNOWN", basis_role="UNKNOWN",
                            reason_code="support_tangent_to_finite_front", front_distance=d_front,
                            crack_distance=d_crack)
    if d_front < h - tol:
        if _point_on_front(point, geometry):
            reason = "front_source_value_undefined_for_shifted_F1"
        elif _point_on_crack(point, geometry):
            reason = "crack_source_value_undefined_for_shifted_F1"
        else:
            reason = "support_contains_finite_front_and_shifted_F1_source_value_is_unique"
        basis_role = "UNKNOWN" if reason in {
            "front_source_value_undefined_for_shifted_F1",
            "crack_source_value_undefined_for_shifted_F1",
        } else "FRONT_F1_SHIFTED"
        return _certificate(source_id, point, geometry, geometric_role="FRONT", basis_role=basis_role,
                            reason_code=reason, front_distance=d_front, crack_distance=d_crack)

    if abs(d_crack - h) <= tol:
        return _certificate(source_id, point, geometry, geometric_role="UNKNOWN", basis_role="UNKNOWN",
                            reason_code="support_tangent_to_finite_crack", front_distance=d_front,
                            crack_distance=d_crack)
    if d_crack > h + tol:
        return _certificate(source_id, point, geometry, geometric_role="BASE", basis_role="BASE",
                            reason_code="support_disjoint_from_finite_crack_and_front", front_distance=d_front,
                            crack_distance=d_crack)

    # The support intersects the finite crack, but does not contain the front.
    # Certify a full split by checking the entire box-clipped plane section.
    plane = _plane_section_bounds(point, geometry)
    dy = abs(point[1] - geometry.y_plane_m)
    side_margin = h - dy
    if plane is None or side_margin <= tol:
        return _certificate(source_id, point, geometry, geometric_role="UNKNOWN", basis_role="UNKNOWN",
                            reason_code="crack_intersection_has_no_two_sided_plane_section", front_distance=d_front,
                            crack_distance=d_crack, plane_bounds=plane, both_sides_margin=side_margin)

    x_lo = plane["x_min_m"]
    x_hi = plane["x_max_m"]
    z_lo = plane["z_min_m"]
    z_hi = plane["z_max_m"]
    crack_z_lo, crack_z_hi = geometry.crack_z_bounds_m
    x_width = x_hi - x_lo
    z_width = z_hi - z_lo

    if abs(x_hi - geometry.front_x_m) <= tol:
        reason = "plane_section_tangent_to_finite_front"
    elif x_hi > geometry.front_x_m + tol:
        reason = "plane_section_crosses_ahead_of_finite_front"
    elif x_lo < geometry.crack_x_min_m:
        reason = "plane_section_extends_before_finite_crack_start"
    elif z_lo < crack_z_lo or z_hi > crack_z_hi:
        reason = "plane_section_extends_beyond_finite_crack_end"
    elif (
        0.0 < x_lo - geometry.crack_x_min_m <= tol
        or 0.0 < z_lo - crack_z_lo <= tol
        or 0.0 < crack_z_hi - z_hi <= tol
    ):
        reason = "plane_section_within_tolerance_of_closed_finite_crack_end"
    elif x_width <= tol or z_width <= tol:
        reason = "plane_section_tangent_to_support_or_box_boundary"
    else:
        reason = "full_finite_crack_section_splits_support_into_two_components"

    if reason != "full_finite_crack_section_splits_support_into_two_components":
        return _certificate(source_id, point, geometry, geometric_role="UNKNOWN", basis_role="UNKNOWN",
                            reason_code=reason, front_distance=d_front, crack_distance=d_crack,
                            plane_bounds=plane, both_sides_margin=side_margin)

    if _point_on_crack(point, geometry) or abs(point[1] - geometry.y_plane_m) <= tol:
        return _certificate(source_id, point, geometry, geometric_role="HEAVISIDE_ELIGIBLE",
                            basis_role="UNKNOWN", reason_code="crack_source_side_value_undefined_for_shifted_H",
                            front_distance=d_front, crack_distance=d_crack, plane_bounds=plane,
                            split_component_count=2, both_sides_margin=side_margin)

    return _certificate(source_id, point, geometry, geometric_role="HEAVISIDE_ELIGIBLE",
                        basis_role="HEAVISIDE", reason_code=reason, front_distance=d_front,
                        crack_distance=d_crack, plane_bounds=plane, split_component_count=2,
                        both_sides_margin=side_margin)


def _canonical_json(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
                      allow_nan=False).encode("utf-8")


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _lattice_source_rows(case: Mapping[str, Any]) -> list[dict[str, Any]]:
    axes = case["sources"]["axes_m"]
    xs, ys, zs = (list(map(float, axes[axis])) for axis in ("x", "y", "z"))
    rows: list[dict[str, Any]] = []
    for ix, x in enumerate(xs):
        for iy, y in enumerate(ys):
            for iz, z in enumerate(zs):
                source_id = f"lf-x{ix:02d}-y{iy:02d}-z{iz:02d}"
                rows.append({"source_id": source_id, "lattice_index": [ix, iy, iz],
                             "position_m": [x, y, z]})
    return rows


def _manifest_payload(manifest: Mapping[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in manifest.items() if key != "manifest_sha256"}


def build_role_manifest(
    case: Mapping[str, Any],
    case_file_bytes: bytes,
    classifier_file_bytes: bytes,
    source_rows: Iterable[Mapping[str, Any]] | None = None,
) -> dict[str, Any]:
    """Build deterministic manifest; optional rows enable permutation checks."""
    geometry = FiniteCrackGeometry.from_case(case)
    expected_rows = _lattice_source_rows(case)
    rows = list(expected_rows if source_rows is None else source_rows)
    expected_by_id = {row["source_id"]: row for row in expected_rows}
    seen_ids: set[str] = set()
    seen_positions: set[tuple[float, float, float]] = set()
    canonical_rows: list[dict[str, Any]] = []
    for row in rows:
        try:
            source_id = row["source_id"]
            lattice_index = list(row["lattice_index"])
            point = _point3(row["position_m"])
        except (KeyError, TypeError) as exc:
            raise GeometryInputError(f"malformed source row: {exc}") from exc
        if source_id in seen_ids or source_id not in expected_by_id:
            raise GeometryInputError("source IDs must be unique and match the frozen lattice IDs")
        expected = expected_by_id[source_id]
        if lattice_index != expected["lattice_index"] or point != tuple(expected["position_m"]):
            raise GeometryInputError("source ID, lattice index, and coordinates do not match fixed-case binding")
        if point in seen_positions:
            raise GeometryInputError("duplicate source coordinates are not allowed")
        seen_ids.add(source_id)
        seen_positions.add(point)
        canonical_rows.append({"source_id": source_id, "lattice_index": lattice_index,
                               "position_m": list(point)})
    if seen_ids != set(expected_by_id):
        raise GeometryInputError("source rows do not contain every frozen lattice ID exactly once")
    canonical_rows.sort(key=lambda row: row["source_id"])

    records: list[dict[str, Any]] = []
    for row in canonical_rows:
        cert = classify_support(row["source_id"], row["position_m"], geometry)
        if cert["status"] != "CERTIFIED":
            raise ManifestBlockedError(f"{row['source_id']} is UNKNOWN: {cert['reason_code']}")
        record = {
            **row,
            "geometric_role": cert["geometric_role"],
            "basis_role": cert["basis_role"],
            "basis_roles": cert["basis_roles"],
            "certificate": {key: value for key, value in cert.items()
                            if key not in {"source_id", "position_m", "geometric_role", "basis_role", "basis_roles", "status"}},
        }
        if len(record["basis_roles"]) not in (1, 2) or record["basis_roles"][0] != "BASE":
            raise ManifestBlockedError("every source must have BASE and at most one enrichment")
        records.append(record)

    role_counts = {
        role: sum(record["geometric_role"] == role for record in records)
        for role in ("FRONT", "HEAVISIDE_ELIGIBLE", "BASE")
    }
    dof_count = sum(len(record["basis_roles"]) for record in records)
    expected_counts = case["expected_counts"]
    if len(records) != int(case["sources"]["source_count"]):
        raise ManifestBlockedError("source count differs from the fixed case")
    if role_counts != expected_counts["geometric_roles"]:
        raise ManifestBlockedError(f"geometric role counts differ: {role_counts}")
    if dof_count != int(expected_counts["scalar_dofs"]):
        raise ManifestBlockedError(f"scalar DOF count differs: {dof_count}")
    budgets = case["budgets"]
    if len(records) > int(budgets["maximum_sources"]):
        raise ManifestBlockedError("source count exceeds the frozen cap")
    if dof_count > int(budgets["maximum_scalar_dofs"]):
        raise ManifestBlockedError("scalar DOF count exceeds the frozen cap")
    if int(budgets["maximum_queries_per_batch"]) > 256:
        raise ManifestBlockedError("query batch cap exceeds 256")

    recordset_sha256 = sha256_bytes(_canonical_json(records))
    manifest: dict[str, Any] = {
        "schema": SCHEMA,
        "case_id": case["case_id"],
        "status": "SOURCE_ONLY_CERTIFIED_GEOMETRY_NOT_COUPLED",
        "units": "m",
        "case_file_sha256": sha256_bytes(case_file_bytes),
        "classifier_file_sha256": sha256_bytes(classifier_file_bytes),
        "source_count": len(records),
        "scalar_dof_count": dof_count,
        "role_counts": role_counts,
        "derived_geometry_bounds": case["derived_geometry_bounds"],
        "maximum_queries_per_batch": int(budgets["maximum_queries_per_batch"]),
        "source_records_sha256": recordset_sha256,
        "sources": records,
    }
    manifest["manifest_sha256"] = sha256_bytes(_canonical_json(manifest))
    return manifest


def verify_manifest_binding(
    manifest: Mapping[str, Any],
    case_file_bytes: bytes,
    classifier_file_bytes: bytes,
) -> None:
    """Check file pins, source-record digest, manifest digest, and key uniqueness."""
    if manifest.get("schema") != SCHEMA:
        raise ManifestBlockedError("unsupported role manifest schema")
    if manifest.get("case_file_sha256") != sha256_bytes(case_file_bytes):
        raise ManifestBlockedError("fixed-case file pin does not match manifest")
    if manifest.get("classifier_file_sha256") != sha256_bytes(classifier_file_bytes):
        raise ManifestBlockedError("classifier source pin does not match manifest")
    records = manifest.get("sources")
    if not isinstance(records, list):
        raise ManifestBlockedError("manifest source records are missing")
    if manifest.get("source_records_sha256") != sha256_bytes(_canonical_json(records)):
        raise ManifestBlockedError("source record digest mismatch")
    if manifest.get("manifest_sha256") != sha256_bytes(_canonical_json(_manifest_payload(manifest))):
        raise ManifestBlockedError("manifest digest mismatch")
    ids = [record.get("source_id") for record in records]
    positions = [tuple(record.get("position_m", ())) for record in records]
    if len(set(ids)) != len(ids) or len(set(positions)) != len(positions):
        raise ManifestBlockedError("source IDs and coordinates must be unique")
    try:
        case = json.loads(case_file_bytes.decode("utf-8-sig"))
        expected_rows = {row["source_id"]: row for row in _lattice_source_rows(case)}
    except (UnicodeDecodeError, json.JSONDecodeError, KeyError, TypeError, ValueError) as exc:
        raise ManifestBlockedError(f"fixed-case binding cannot be read: {exc}") from exc
    if len(records) != len(expected_rows) or set(ids) != set(expected_rows):
        raise ManifestBlockedError("manifest IDs do not match the frozen source lattice")
    actual_roles = {"FRONT": 0, "HEAVISIDE_ELIGIBLE": 0, "BASE": 0}
    actual_dofs = 0
    for record in records:
        roles = record.get("basis_roles")
        if not isinstance(roles, list) or not roles or roles[0] != "BASE" or len(roles) > 2:
            raise ManifestBlockedError("each source must have BASE and at most one enrichment")
        source_id = record["source_id"]
        expected = expected_rows[source_id]
        if record.get("lattice_index") != expected["lattice_index"] or tuple(record.get("position_m", ())) != tuple(expected["position_m"]):
            raise ManifestBlockedError("source ID does not bind to its frozen lattice coordinates")
        geometric_role = record.get("geometric_role")
        role_to_keys = {
            "FRONT": ["BASE", "FRONT_F1_SHIFTED"],
            "HEAVISIDE_ELIGIBLE": ["BASE", "HEAVISIDE"],
            "BASE": ["BASE"],
        }
        if geometric_role not in role_to_keys or roles != role_to_keys[geometric_role]:
            raise ManifestBlockedError("basis roles do not match the single certified geometric role")
        if roles != ["BASE"] and abs(float(record["position_m"][1]) - float(case["geometry"]["crack"]["y_plane_m"])) <= float(case["topology"]["tolerance_m"]):
            raise ManifestBlockedError("shifted enrichment source lies on the branch set")
        actual_roles[geometric_role] += 1
        actual_dofs += len(roles)
    if actual_roles != manifest.get("role_counts") or actual_dofs != manifest.get("scalar_dof_count") or len(records) != manifest.get("source_count"):
        raise ManifestBlockedError("manifest counts do not match bound source records")


def load_case(path: str) -> tuple[dict[str, Any], bytes]:
    """Read the frozen JSON case and preserve its exact byte pin."""
    with open(path, "rb") as stream:
        raw = stream.read()
    try:
        case = json.loads(raw.decode("utf-8-sig"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise GeometryInputError(f"fixed-case JSON is invalid: {exc}") from exc
    if not isinstance(case, dict):
        raise GeometryInputError("fixed-case root must be a JSON object")
    return case, raw


def write_role_manifest(case_path: str, output_path: str) -> dict[str, Any]:
    """Materialize the frozen role records without invoking a provider or solver."""
    case, case_bytes = load_case(case_path)
    classifier_bytes = Path(__file__).read_bytes()
    manifest = build_role_manifest(case, case_bytes, classifier_bytes)
    Path(output_path).write_bytes(_canonical_json(manifest) + b"\n")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the finite-crack support role manifest")
    parser.add_argument("--case", required=True, help="fixed-case JSON path")
    parser.add_argument("--output", required=True, help="role manifest JSON path")
    args = parser.parse_args()
    manifest = write_role_manifest(args.case, args.output)
    print(json.dumps({
        "case_id": manifest["case_id"],
        "source_count": manifest["source_count"],
        "scalar_dof_count": manifest["scalar_dof_count"],
        "role_counts": manifest["role_counts"],
        "manifest_sha256": manifest["manifest_sha256"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
