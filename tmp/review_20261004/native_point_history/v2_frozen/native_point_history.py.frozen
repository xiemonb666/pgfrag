"""Read the explicit CSV/TOML archive emitted by the native PD point probe.

This module treats every archive as caller-supplied, unverified data.  It only
parses and cross-checks the frozen CSV schema; it never imports bond_events or
claims that a summary status or process exit code authenticates an archive.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
import hashlib
import io
import json
import math
from pathlib import Path
import re
import tomllib
from typing import Any, Mapping


NATIVE_ORIGIN_STATUS = "CALLER_SUPPLIED_ARCHIVE_NOT_AUTHENTICATED"
POINT_STATE_HEADER = (
    "step,time_s,dt_s,native_body_index,point_id,x_ref_m,y_ref_m,z_ref_m,"
    "x_current_m,y_current_m,z_current_m,ux_m,uy_m,uz_m,vx_mps,vy_mps,"
    "vz_mps,ax_mps2,ay_mps2,az_mps2,b_int_x_N_m3,b_int_y_N_m3,b_int_z_N_m3,"
    "b_ext_x_N_m3,b_ext_y_N_m3,b_ext_z_N_m3,damage,n_active_bonds,volume_m3,mass_kg"
).split(",")
DIRECTED_STATE_HEADER = (
    "step,time_s,dt_s,directed_source_point_id,directed_neighbor_point_id,"
    "canonical_point_i_id,canonical_point_j_id,directed_active,canonical_active"
).split(",")
CANONICAL_STATE_HEADER = "step,time_s,dt_s,point_i_id,point_j_id,active".split(",")
POINT_INPUT_HEADER = "point_id,x_m,y_m,z_m,volume_m3".split(",")
SNAPSHOT_FILES = (
    ("summary", "summary.toml"),
    ("point_state_history", "point_state_history.csv"),
    ("directed_bond_state_history", "directed_bond_state_history.csv"),
    ("canonical_bond_state_history", "canonical_bond_state_history.csv"),
)


@dataclass(frozen=True, slots=True)
class FileFingerprint:
    role: str
    reference: str
    bytes: int | None
    sha256: str | None


@dataclass(frozen=True, slots=True)
class InputPoint:
    point_id: int
    native_body_index: int
    rest_position_m: tuple[float, float, float]
    volume_m3: float


@dataclass(frozen=True, slots=True)
class NativePointStateRow:
    step: int
    time_s: float
    accepted_dt_s: float
    native_body_index: int
    point_id: int
    rest_position_m: tuple[float, float, float]
    volume_m3: float
    mass_kg: float
    n_active_bonds: int
    line_number: int


@dataclass(frozen=True, slots=True)
class NativeDirectedStateRow:
    step: int
    time_s: float
    accepted_dt_s: float
    source_point_id: int
    target_point_id: int
    canonical_point_i_id: int
    canonical_point_j_id: int
    active: bool
    canonical_active: bool
    native_field: str
    line_number: int


@dataclass(frozen=True, slots=True)
class NativeCanonicalStateRow:
    step: int
    time_s: float
    accepted_dt_s: float
    point_i_id: int
    point_j_id: int
    active: bool
    line_number: int


@dataclass(frozen=True, slots=True)
class BondEventsDirectedRecord:
    """Raw field-compatible record for the frozen bond_events v1 contract."""

    source_point_id: int
    target_point_id: int
    active: bool
    native_field: str


@dataclass(frozen=True, slots=True)
class RawNativeBondFrame:
    """Validated native frame plus the source rows that formed it.

    ``step_dt_s`` follows the frozen bond_events v1 convention and is None at
    step zero. ``accepted_dt_s`` preserves the probe's dt column at every step.
    """

    step: int
    time_s: float
    step_dt_s: float | None
    accepted_dt_s: float
    points: tuple[NativePointStateRow, ...]
    directed_bonds: tuple[BondEventsDirectedRecord, ...]
    canonical_bonds: tuple[tuple[int, int, bool], ...]


@dataclass(frozen=True, slots=True)
class NativePointHistory:
    status: str
    native_origin_status: str
    reasons: tuple[str, ...]
    snapshot_reference: str
    source_reference: str | None
    input_id: str | None
    summary_status: str | None
    expected_final_step: int | None
    expected_dt_s: float | None
    density_kg_m3: float | None
    point_ids: tuple[int, ...]
    input_points: tuple[InputPoint, ...]
    rest_positions_m: Mapping[int, tuple[float, float, float]]
    point_volumes_m3: Mapping[int, float]
    point_masses_kg: Mapping[int, float]
    initial_bonds: tuple[tuple[int, int], ...]
    frames: tuple[RawNativeBondFrame, ...]
    raw_point_rows: tuple[NativePointStateRow, ...]
    raw_directed_rows: tuple[NativeDirectedStateRow, ...]
    raw_canonical_rows: tuple[NativeCanonicalStateRow, ...]
    summary: Mapping[str, Any]
    files: tuple[FileFingerprint, ...]
    snapshot_manifest_sha256: str | None

    @property
    def bond_events_mapping(self) -> dict[str, Any] | None:
        """Return plain bond_events-shaped inputs only for a complete archive.

        The caller must instantiate its own bond_events dataclasses and
        provenance type. This adapter intentionally does not import that
        module, which is being revised independently.
        """

        if self.status != "READY" or self.snapshot_manifest_sha256 is None:
            return None
        observer = next((item for item in self.files if item.role == "observer_schema_source"), None)
        return {
            "point_ids": list(self.point_ids),
            "rest_positions_m": dict(self.rest_positions_m),
            "point_volumes_m3": dict(self.point_volumes_m3),
            "density_kg_m3": self.density_kg_m3,
            "initial_bonds": [list(pair) for pair in self.initial_bonds],
            "frames": [
                {
                    "step": frame.step,
                    "time_s": frame.time_s,
                    "step_dt_s": frame.step_dt_s,
                    "directed_bonds": [
                        {
                            "source_point_id": record.source_point_id,
                            "target_point_id": record.target_point_id,
                            "active": record.active,
                            "native_field": record.native_field,
                        }
                        for record in frame.directed_bonds
                    ],
                }
                for frame in self.frames
            ],
            "expected_final_step": self.expected_final_step,
            "expected_dt_s": self.expected_dt_s,
            "provenance": {
                "native_source": self.source_reference or "caller-supplied native point archive",
                "native_interface": (
                    "pd_point_event_probe CSV schema: directed_active and canonical active rows"
                ),
                "frame_archive_sha256": self.snapshot_manifest_sha256,
                "native_source_sha256": observer.sha256 if observer else None,
                "native_origin_status": NATIVE_ORIGIN_STATUS,
            },
        }

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-friendly diagnostic view, including raw fingerprints."""

        return {
            "schema": "pgfrag.native_point_history.v1",
            "status": self.status,
            "native_origin_status": self.native_origin_status,
            "reasons": list(self.reasons),
            "snapshot_reference": self.snapshot_reference,
            "source_reference": self.source_reference,
            "input_id": self.input_id,
            "summary_status": self.summary_status,
            "expected_final_step": self.expected_final_step,
            "expected_dt_s": self.expected_dt_s,
            "density_kg_m3": self.density_kg_m3,
            "point_ids": list(self.point_ids),
            "initial_bonds": [list(pair) for pair in self.initial_bonds],
            "raw_prefix": {
                "point_rows": len(self.raw_point_rows),
                "directed_rows": len(self.raw_directed_rows),
                "canonical_rows": len(self.raw_canonical_rows),
                "validated_frames": len(self.frames),
            },
            "snapshot_manifest_sha256": self.snapshot_manifest_sha256,
            "files": [
                {
                    "role": item.role,
                    "reference": item.reference,
                    "bytes": item.bytes,
                    "sha256": item.sha256,
                }
                for item in self.files
            ],
            "bond_events_mapping_available": self.bond_events_mapping is not None,
            "interpretation": {
                "native_origin_authentication": "NOT_PERFORMED_CALLER_SUPPLIED_ARCHIVE",
                "point_state_projection": (
                    "native point CSV numeric columns are validated, while projected records retain the requested "
                    "identity, rest geometry, volume, mass, timing, and active-bond count; complete CSV bytes are "
                    "represented by their SHA-256 and byte count"
                ),
                "physical_fragments": "UNKNOWN_NOT_DERIVED",
                "positive_area_crack_surface": "UNKNOWN_NOT_DERIVED",
                "energy": "UNKNOWN_NOT_DERIVED",
                "R2": "UNKNOWN_NOT_DERIVED",
            },
        }


@dataclass(frozen=True, slots=True)
class _CsvTable:
    rows: tuple[dict[str, str], ...]
    errors: tuple[str, ...]


def read_native_point_history(
    snapshot_dir: str | Path,
    point_csv_path: str | Path,
    point_metadata_path: str | Path,
    *,
    observer_source_path: str | Path | None = None,
) -> NativePointHistory:
    """Read and validate one explicit native point-history snapshot.

    A missing or empty snapshot is ``NOT_RUN``. Any started but structurally
    incomplete, truncated, or inconsistent snapshot is ``INCOMPLETE``. Only a
    complete prefix covering the caller-supplied summary's steps is ``READY``;
    that status is a structural result and does not authenticate provenance.
    """

    snapshot = Path(snapshot_dir).expanduser().resolve()
    points_path = Path(point_csv_path).expanduser().resolve()
    metadata_path = Path(point_metadata_path).expanduser().resolve()
    if observer_source_path is None:
        observer_path = Path(__file__).resolve().parents[2] / "experiments" / "pd_point_event_probe_20261004" / "run_pd_point_event_probe.jl"
    else:
        observer_path = Path(observer_source_path).expanduser().resolve()

    fingerprints: list[FileFingerprint] = []

    def capture(role: str, path: Path) -> bytes | None:
        try:
            data = path.read_bytes()
        except OSError:
            fingerprints.append(FileFingerprint(role, str(path), None, None))
            return None
        fingerprints.append(
            FileFingerprint(role, str(path), len(data), hashlib.sha256(data).hexdigest())
        )
        return data

    points_bytes = capture("point_input_csv", points_path)
    metadata_bytes = capture("point_input_metadata", metadata_path)
    capture("observer_schema_source", observer_path)
    snapshot_present = snapshot.is_dir()
    snapshot_data: dict[str, bytes | None] = {}
    for role, name in SNAPSHOT_FILES:
        data = capture(f"snapshot_{role}", snapshot / name)
        snapshot_data[role] = data

    reasons: list[str] = []
    metadata: dict[str, Any] = {}
    source_reference: str | None = None
    input_id: str | None = None
    input_points: tuple[InputPoint, ...] = ()
    input_errors: list[str] = []
    if points_bytes is None:
        input_errors.append(f"missing point input CSV: {points_path}")
    if metadata_bytes is None:
        input_errors.append(f"missing point input metadata: {metadata_path}")
    if points_bytes is not None and metadata_bytes is not None:
        try:
            metadata = tomllib.loads(metadata_bytes.decode("utf-8"))
        except (UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
            input_errors.append(f"point input metadata is unreadable: {exc}")
        else:
            source = metadata.get("source")
            if isinstance(source, dict):
                source_reference = _nonempty_string(source.get("reference"))
            input_id = _nonempty_string(metadata.get("input_id"))
            input_points, point_errors = _parse_input_points(points_bytes, metadata)
            input_errors.extend(point_errors)
            actual_hash = hashlib.sha256(points_bytes).hexdigest()
            declared_hash = metadata.get("point_data_sha256")
            if declared_hash != actual_hash:
                input_errors.append("point CSV SHA-256 disagrees with point_data_sha256 in metadata")
            if not source_reference:
                input_errors.append("point metadata source.reference is missing")
            if metadata.get("schema") != "pgfrag.pd.pointset.v1":
                input_errors.append("unsupported point-input metadata schema")
            if metadata.get("coordinate_units") != "m" or metadata.get("volume_units") != "m^3":
                input_errors.append("point input units must be metres and cubic metres")
            if not isinstance(source, dict):
                input_errors.append("point metadata source provenance is missing")
            elif not all(_nonempty_string(source.get(key)) for key in ("kind", "description", "reference")):
                input_errors.append("point metadata source kind, description, and reference are required")
            elif not _is_sha256(source.get("source_sha256")):
                input_errors.append("point metadata source.source_sha256 must be a SHA-256 digest")

    output_count = sum(snapshot_data[role] is not None for role, _ in SNAPSHOT_FILES)
    if output_count == 0:
        # An empty output directory can be left by a launcher that failed before
        # writing results. It is not an incomplete native frame archive.
        return _result(
            status="NOT_RUN",
            snapshot=snapshot,
            reasons=tuple(input_errors),
            source_reference=source_reference,
            input_id=input_id,
            input_points=input_points,
            fingerprints=fingerprints,
            summary={},
        )
    if not snapshot_present:
        reasons.append(f"snapshot directory is missing: {snapshot}")
    elif output_count < len(SNAPSHOT_FILES):
        missing = [name for role, name in SNAPSHOT_FILES if snapshot_data[role] is None]
        reasons.append(f"snapshot is incomplete; missing required output(s): {', '.join(missing)}")
    reasons.extend(input_errors)

    summary: dict[str, Any] = {}
    summary_bytes = snapshot_data.get("summary")
    if summary_bytes is not None:
        try:
            summary = tomllib.loads(summary_bytes.decode("utf-8"))
        except (UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
            reasons.append(f"summary.toml is unreadable: {exc}")

    table_specs = (
        ("point_state_history", POINT_STATE_HEADER),
        ("directed_bond_state_history", DIRECTED_STATE_HEADER),
        ("canonical_bond_state_history", CANONICAL_STATE_HEADER),
    )
    tables: dict[str, _CsvTable] = {}
    for role, header in table_specs:
        data = snapshot_data.get(role)
        if data is None:
            tables[role] = _CsvTable((), ())
            continue
        table = _read_csv_table(data, header, role, require_final_newline=True)
        tables[role] = table
        reasons.extend(table.errors)

    raw_point_rows, point_row_errors = _decode_point_rows(tables["point_state_history"])
    raw_directed_rows, directed_row_errors = _decode_directed_rows(tables["directed_bond_state_history"])
    raw_canonical_rows, canonical_row_errors = _decode_canonical_rows(tables["canonical_bond_state_history"])
    reasons.extend(point_row_errors)
    reasons.extend(directed_row_errors)
    reasons.extend(canonical_row_errors)

    expected_final_step = _summary_int(summary, "n_steps", reasons, minimum=1)
    expected_dt_s = _summary_float(summary, "dt_s", reasons, positive=True)
    density = _summary_float(summary, "rho_kg_m3", reasons, positive=True)
    if expected_final_step is not None and expected_dt_s is not None:
        final_time = _summary_float(summary, "time_s", reasons, positive=False)
        expected_time = _finite_product(expected_final_step, expected_dt_s)
        if expected_time is None:
            reasons.append("summary n_steps × accepted dt_s product is not finite")
        elif final_time is not None and not _close(final_time, expected_time, expected_dt_s):
            reasons.append("summary time_s disagrees with n_steps × accepted dt_s")
    if summary and "status" in summary:
        summary_status = str(summary["status"])
    else:
        summary_status = None
    if points_bytes is not None and summary.get("point_input_sha256") != hashlib.sha256(points_bytes).hexdigest():
        reasons.append("summary point_input_sha256 disagrees with the supplied point CSV")
    if input_id is not None and summary.get("input_id") != input_id:
        reasons.append("summary input_id disagrees with the supplied point metadata")

    point_ids = tuple(point.point_id for point in input_points)
    point_set = set(point_ids)
    input_by_id = {point.point_id: point for point in input_points}
    point_groups = _group_rows(raw_point_rows)
    directed_groups = _group_rows(raw_directed_rows)
    canonical_groups = _group_rows(raw_canonical_rows)
    _check_row_step_order(raw_point_rows, "point state history", reasons)
    _check_row_step_order(raw_directed_rows, "directed bond state history", reasons)
    _check_row_step_order(raw_canonical_rows, "canonical bond state history", reasons)

    if metadata and metadata.get("point_count") != len(input_points):
        reasons.append("parsed point count disagrees with metadata point_count")
    _check_summary_count(summary, "source_point_count", len(input_points), reasons)
    source_volume_sum = _finite_sum((point.volume_m3 for point in input_points), "source point volume sum", reasons)
    summary_volume_sum = _summary_float(summary, "point_volume_sum_m3", reasons, positive=True)
    if source_volume_sum is not None and summary_volume_sum is not None and not _close(
        source_volume_sum, summary_volume_sum, source_volume_sum
    ):
        reasons.append("summary point_volume_sum_m3 disagrees with the supplied point CSV")
    if density is not None and source_volume_sum is not None:
        point_mass_products = [_finite_product(density, point.volume_m3) for point in input_points]
        if any(value is None or value <= 0.0 for value in point_mass_products):
            reasons.append("density × point volume product is not positive finite for the supplied point CSV")
        else:
            point_mass_sum = _finite_sum(point_mass_products, "density × point-volume product sum", reasons)
            expected_mass_sum = _finite_product(density, source_volume_sum)
            if expected_mass_sum is None or expected_mass_sum <= 0.0:
                reasons.append("summary density × point-volume-sum product is not positive finite")
            else:
                summary_mass_sum = _summary_float(summary, "point_mass_sum_kg", reasons, positive=True)
                if point_mass_sum is not None and not _close(point_mass_sum, expected_mass_sum, expected_mass_sum):
                    reasons.append("density × point-volume products disagree with density × source volume sum")
                if summary_mass_sum is not None and not _close(summary_mass_sum, expected_mass_sum, expected_mass_sum):
                    reasons.append("summary point_mass_sum_kg disagrees with density × source volume sum")

    frames: list[RawNativeBondFrame] = []
    previous_canonical: dict[tuple[int, int], bool] | None = None
    initial_bonds: tuple[tuple[int, int], ...] = ()
    graph_pairs: set[tuple[int, int]] | None = None
    expected_steps = expected_final_step
    if expected_steps is not None:
        for label, groups in (
            ("point", point_groups),
            ("directed bond", directed_groups),
            ("canonical bond", canonical_groups),
        ):
            if any(step > expected_steps for step in groups):
                reasons.append(f"{label} history includes a step beyond summary n_steps")
    if expected_steps is not None and point_groups:
        observed_max = max(point_groups)
        if observed_max > expected_steps:
            reasons.append("point history includes a step beyond summary n_steps")
        last_to_check = min(observed_max, expected_steps)
        for step in range(last_to_check + 1):
            point_rows = point_groups.get(step)
            directed_rows = directed_groups.get(step, [])
            canonical_rows = canonical_groups.get(step, [])
            if point_rows is None:
                reasons.append(f"point history is missing native frame step {step}")
                break
            if len(point_rows) != len(point_ids):
                reasons.append(
                    f"point frame {step} has {len(point_rows)} rows; expected {len(point_ids)}"
                )
                break
            ids_here = [row.point_id for row in point_rows]
            indexes_here = [row.native_body_index for row in point_rows]
            if len(set(ids_here)) != len(ids_here):
                reasons.append(f"point frame {step} contains duplicate stable point IDs")
                break
            if set(ids_here) != point_set:
                missing = sorted(point_set - set(ids_here))
                extra = sorted(set(ids_here) - point_set)
                reasons.append(f"point frame {step} stable point set mismatch (missing={missing}, extra={extra})")
                break
            if set(indexes_here) != set(range(1, len(input_points) + 1)) or len(set(indexes_here)) != len(indexes_here):
                reasons.append(f"point frame {step} has duplicate, missing, or invalid native body indexes")
                break
            row_by_id = {row.point_id: row for row in point_rows}
            if any(row_by_id[point.point_id].native_body_index != point.native_body_index for point in input_points):
                reasons.append(f"point frame {step} native body indexes disagree with input CSV row order")
                break
            bad_point: str | None = None
            for point_id, row in row_by_id.items():
                source = input_by_id[point_id]
                if row.rest_position_m != source.rest_position_m:
                    bad_point = f"point frame {step} rest position for point {point_id} disagrees with source CSV"
                    break
                if row.volume_m3 != source.volume_m3:
                    bad_point = f"point frame {step} volume for point {point_id} disagrees with source CSV"
                    break
                if density is not None:
                    expected_mass = density * row.volume_m3
                    if not math.isfinite(expected_mass) or expected_mass <= 0.0:
                        bad_point = f"point frame {step} density × volume mass product is not positive finite for point {point_id}"
                        break
                    if not _close(row.mass_kg, expected_mass, max(abs(row.mass_kg), abs(expected_mass))):
                        bad_point = f"point frame {step} mass for point {point_id} disagrees with density × volume"
                        break
            if bad_point:
                reasons.append(bad_point)
                break

            clocks = [(row.time_s, row.accepted_dt_s) for row in point_rows]
            if not _clock_group_valid(step, clocks, expected_dt_s):
                reasons.append(f"point frame {step} has inconsistent time or accepted dt")
                break
            time_s, accepted_dt = clocks[0]
            if directed_rows and not _clock_group_valid(
                step, [(row.time_s, row.accepted_dt_s) for row in directed_rows], expected_dt_s
            ):
                reasons.append(f"directed bond frame {step} has inconsistent time or accepted dt")
                break
            if canonical_rows and not _clock_group_valid(
                step, [(row.time_s, row.accepted_dt_s) for row in canonical_rows], expected_dt_s
            ):
                reasons.append(f"canonical bond frame {step} has inconsistent time or accepted dt")
                break
            if directed_rows and any(
                row.time_s != time_s or row.accepted_dt_s != accepted_dt for row in directed_rows
            ):
                reasons.append(f"directed bond frame {step} clock disagrees with point state frame")
                break
            if canonical_rows and any(
                row.time_s != time_s or row.accepted_dt_s != accepted_dt for row in canonical_rows
            ):
                reasons.append(f"canonical bond frame {step} clock disagrees with point state frame")
                break

            directed_by_pair: dict[tuple[int, int], dict[tuple[int, int], NativeDirectedStateRow]] = {}
            invalid_bond = False
            for row in directed_rows:
                direction = (row.source_point_id, row.target_point_id)
                pair = (min(direction), max(direction))
                if row.source_point_id == row.target_point_id or not set(direction) <= point_set:
                    reasons.append(f"directed bond frame {step} contains a self-bond or unknown point ID")
                    invalid_bond = True
                    break
                if (row.canonical_point_i_id, row.canonical_point_j_id) != pair:
                    reasons.append(f"directed bond frame {step} canonical pair disagrees with directed endpoints")
                    invalid_bond = True
                    break
                by_direction = directed_by_pair.setdefault(pair, {})
                if direction in by_direction:
                    reasons.append(f"directed bond frame {step} contains duplicate direction {direction}")
                    invalid_bond = True
                    break
                by_direction[direction] = row
            if invalid_bond:
                break
            canonical_state: dict[tuple[int, int], bool] = {}
            for row in canonical_rows:
                pair = (row.point_i_id, row.point_j_id)
                if row.point_i_id >= row.point_j_id or not set(pair) <= point_set:
                    reasons.append(f"canonical bond frame {step} contains a non-canonical or unknown endpoint pair")
                    invalid_bond = True
                    break
                if pair in canonical_state:
                    reasons.append(f"canonical bond frame {step} contains duplicate pair {pair}")
                    invalid_bond = True
                    break
                canonical_state[pair] = row.active
            if invalid_bond:
                break
            directed_pairs = set(directed_by_pair)
            canonical_pairs = set(canonical_state)
            if directed_pairs != canonical_pairs:
                reasons.append(f"directed/canonical bond key sets disagree at step {step}")
                break
            if graph_pairs is None:
                graph_pairs = canonical_pairs
                initial_bonds = tuple(sorted(canonical_pairs))
            elif canonical_pairs != graph_pairs:
                reasons.append(f"native canonical bond key set changed at step {step}")
                break

            canonical_directed_records: list[BondEventsDirectedRecord] = []
            for pair in sorted(canonical_pairs):
                directions = directed_by_pair[pair]
                forward = directions.get(pair)
                reverse = directions.get((pair[1], pair[0]))
                if forward is None or reverse is None:
                    reasons.append(f"directed bond frame {step} is missing a reciprocal direction for {pair}")
                    invalid_bond = True
                    break
                if forward.active != reverse.active:
                    reasons.append(f"reciprocal native bond states disagree for {pair} at step {step}")
                    invalid_bond = True
                    break
                if forward.active != canonical_state[pair]:
                    reasons.append(f"directed/canonical native bond states disagree for {pair} at step {step}")
                    invalid_bond = True
                    break
                if forward.canonical_active != canonical_state[pair] or reverse.canonical_active != canonical_state[pair]:
                    reasons.append(f"directed canonical_active column disagrees for {pair} at step {step}")
                    invalid_bond = True
                    break
                canonical_directed_records.extend(
                    (
                        BondEventsDirectedRecord(pair[0], pair[1], forward.active, forward.native_field),
                        BondEventsDirectedRecord(pair[1], pair[0], reverse.active, reverse.native_field),
                    )
                )
            if invalid_bond:
                break
            active_counts = {point_id: 0 for point_id in point_set}
            for row in directed_rows:
                if row.active:
                    active_counts[row.source_point_id] += 1
            if any(row_by_id[point_id].n_active_bonds != active_counts[point_id] for point_id in point_set):
                reasons.append(f"point frame {step} n_active_bonds disagrees with directed native active states")
                break
            if step == 0 and any(not active for active in canonical_state.values()):
                reasons.append("initial native frame contains an inactive bond")
                break
            if previous_canonical is not None:
                healed = sorted(pair for pair, active in canonical_state.items() if not previous_canonical[pair] and active)
                if healed:
                    reasons.append(f"native bond healing detected at step {step}: {healed}")
                    break
            previous_canonical = canonical_state
            frames.append(
                RawNativeBondFrame(
                    step=step,
                    time_s=time_s,
                    step_dt_s=None if step == 0 else accepted_dt,
                    accepted_dt_s=accepted_dt,
                    points=tuple(sorted(point_rows, key=lambda row: row.native_body_index)),
                    directed_bonds=tuple(canonical_directed_records),
                    canonical_bonds=tuple((first, second, active) for (first, second), active in sorted(canonical_state.items())),
                )
            )
    elif expected_steps is not None:
        reasons.append("point history contains no complete frame rows")
    if expected_final_step is not None and frames and frames[-1].step != expected_final_step:
        reasons.append(
            f"point history ends at step {frames[-1].step}; summary expects final step {expected_final_step}"
        )
    if expected_final_step is not None:
        _check_summary_count(summary, "sample_count_including_step_zero", expected_final_step + 1, reasons)
    _check_summary_count(summary, "directed_bond_state_rows", len(raw_directed_rows), reasons)
    _check_summary_count(summary, "canonical_bond_state_rows", len(raw_canonical_rows), reasons)
    if frames:
        observed_mass_sum = _finite_sum((row.mass_kg for row in frames[0].points), "step-zero point mass sum", reasons)
        summary_mass_sum = _summary_float(summary, "point_mass_sum_kg", reasons, positive=True)
        if observed_mass_sum is not None and summary_mass_sum is not None and not _close(
            observed_mass_sum, summary_mass_sum, observed_mass_sum
        ):
            reasons.append("summary point_mass_sum_kg disagrees with the step-zero point state")

    if not input_points:
        rest_positions: dict[int, tuple[float, float, float]] = {}
        volumes: dict[int, float] = {}
        masses: dict[int, float] = {}
    else:
        rest_positions = {point.point_id: point.rest_position_m for point in input_points}
        volumes = {point.point_id: point.volume_m3 for point in input_points}
        if frames:
            masses = {row.point_id: row.mass_kg for row in frames[0].points}
        else:
            masses = {}

    snapshot_fingerprints = [item for item in fingerprints if item.role.startswith("snapshot_")]
    snapshot_manifest_sha256 = None
    if all(item.sha256 is not None for item in snapshot_fingerprints) and len(snapshot_fingerprints) == len(SNAPSHOT_FILES):
        manifest_payload = [
            {"role": item.role, "reference": Path(item.reference).name, "bytes": item.bytes, "sha256": item.sha256}
            for item in sorted(snapshot_fingerprints, key=lambda item: item.role)
        ]
        encoded = json.dumps(manifest_payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
        snapshot_manifest_sha256 = hashlib.sha256(encoded).hexdigest()

    if reasons:
        status = "INCOMPLETE"
    elif expected_final_step is None or expected_dt_s is None or density is None:
        status = "INCOMPLETE"
        reasons.append("summary timing or density fields are unavailable")
    elif len(frames) != expected_final_step + 1 or frames[0].step != 0:
        status = "INCOMPLETE"
        reasons.append("validated native frame prefix does not cover step zero through the expected final step")
    else:
        status = "READY"
    return NativePointHistory(
        status=status,
        native_origin_status=NATIVE_ORIGIN_STATUS,
        reasons=tuple(dict.fromkeys(reasons)),
        snapshot_reference=str(snapshot),
        source_reference=source_reference,
        input_id=input_id,
        summary_status=summary_status,
        expected_final_step=expected_final_step,
        expected_dt_s=expected_dt_s,
        density_kg_m3=density,
        point_ids=point_ids,
        input_points=input_points,
        rest_positions_m=rest_positions,
        point_volumes_m3=volumes,
        point_masses_kg=masses,
        initial_bonds=initial_bonds,
        frames=tuple(frames),
        raw_point_rows=raw_point_rows,
        raw_directed_rows=raw_directed_rows,
        raw_canonical_rows=raw_canonical_rows,
        summary=summary,
        files=tuple(fingerprints),
        snapshot_manifest_sha256=snapshot_manifest_sha256,
    )


def _result(
    *,
    status: str,
    snapshot: Path,
    reasons: tuple[str, ...],
    source_reference: str | None,
    input_id: str | None,
    input_points: tuple[InputPoint, ...],
    fingerprints: list[FileFingerprint],
    summary: Mapping[str, Any],
) -> NativePointHistory:
    density_raw = summary.get("rho_kg_m3")
    density = float(density_raw) if isinstance(density_raw, (int, float)) and not isinstance(density_raw, bool) else None
    expected_dt_raw = summary.get("dt_s")
    expected_dt = float(expected_dt_raw) if isinstance(expected_dt_raw, (int, float)) and not isinstance(expected_dt_raw, bool) else None
    expected_step_raw = summary.get("n_steps")
    expected_step = expected_step_raw if isinstance(expected_step_raw, int) and not isinstance(expected_step_raw, bool) else None
    points = {point.point_id: point.rest_position_m for point in input_points}
    volumes = {point.point_id: point.volume_m3 for point in input_points}
    return NativePointHistory(
        status=status,
        native_origin_status=NATIVE_ORIGIN_STATUS,
        reasons=reasons,
        snapshot_reference=str(snapshot),
        source_reference=source_reference,
        input_id=input_id,
        summary_status=_nonempty_string(summary.get("status")),
        expected_final_step=expected_step,
        expected_dt_s=expected_dt,
        density_kg_m3=density,
        point_ids=tuple(point.point_id for point in input_points),
        input_points=input_points,
        rest_positions_m=points,
        point_volumes_m3=volumes,
        point_masses_kg={},
        initial_bonds=(),
        frames=(),
        raw_point_rows=(),
        raw_directed_rows=(),
        raw_canonical_rows=(),
        summary=summary,
        files=tuple(fingerprints),
        snapshot_manifest_sha256=None,
    )


def _parse_input_points(data: bytes, metadata: Mapping[str, Any]) -> tuple[tuple[InputPoint, ...], list[str]]:
    table = _read_csv_table(data, POINT_INPUT_HEADER, "point input CSV", require_final_newline=False)
    errors = list(table.errors)
    records: list[InputPoint] = []
    seen: set[int] = set()
    for index, row in enumerate(table.rows, start=2):
        try:
            point_id = _int_field(row["point_id"], f"point input row {index} point_id", minimum=1)
            xyz = tuple(_float_field(row[key], f"point input row {index} {key}") for key in ("x_m", "y_m", "z_m"))
            volume = _float_field(row["volume_m3"], f"point input row {index} volume_m3", positive=True)
        except ValueError as exc:
            errors.append(str(exc))
            break
        if point_id in seen:
            errors.append(f"point input CSV contains duplicate point ID {point_id}")
            break
        seen.add(point_id)
        records.append(InputPoint(point_id, len(records) + 1, xyz, volume))
    if table.rows and not records and not errors:
        errors.append("point input CSV has no valid point rows")
    declared_count = metadata.get("point_count")
    if isinstance(declared_count, int) and not isinstance(declared_count, bool) and declared_count != len(records):
        errors.append("point input row count disagrees with metadata point_count")
    if records and isinstance(metadata.get("declared_volume_sum_m3"), (int, float)) and not isinstance(
        metadata.get("declared_volume_sum_m3"), bool
    ):
        actual_volume_sum = _finite_sum((point.volume_m3 for point in records), "source point volume sum", errors)
        try:
            declared_volume_sum = float(metadata["declared_volume_sum_m3"])
        except (OverflowError, TypeError, ValueError):
            declared_volume_sum = math.nan
        if not math.isfinite(declared_volume_sum):
            errors.append("metadata declared_volume_sum_m3 must be finite")
        elif actual_volume_sum is not None and not _close(actual_volume_sum, declared_volume_sum, actual_volume_sum):
            errors.append("point CSV volume sum disagrees with metadata declared_volume_sum_m3")
    declared_com = metadata.get("declared_com_m")
    if records and isinstance(declared_com, list) and len(declared_com) == 3:
        total_volume = _finite_sum((point.volume_m3 for point in records), "source point volume sum", errors)
        if total_volume is None or total_volume <= 0.0:
            return tuple(records), errors
        weighted_coordinates: list[float] = []
        for axis in range(3):
            terms = [point.rest_position_m[axis] * point.volume_m3 for point in records]
            if not all(math.isfinite(term) for term in terms):
                errors.append(f"point CSV weighted coordinate product is not finite on axis {axis}")
                return tuple(records), errors
            weighted_sum = _finite_sum(terms, f"point CSV weighted coordinate sum on axis {axis}", errors)
            if weighted_sum is None:
                return tuple(records), errors
            coordinate = weighted_sum / total_volume
            if not math.isfinite(coordinate):
                errors.append(f"point CSV volume-weighted coordinate is not finite on axis {axis}")
                return tuple(records), errors
            weighted_coordinates.append(coordinate)
        measured = tuple(weighted_coordinates)
        try:
            declared = tuple(float(value) for value in declared_com)
        except (OverflowError, TypeError, ValueError):
            errors.append("metadata declared_com_m must contain three numeric values")
        else:
            if not all(math.isfinite(value) for value in declared) or any(
                not _close(actual, expected, total_volume) for actual, expected in zip(measured, declared, strict=True)
            ):
                errors.append("point CSV volume-weighted COM disagrees with metadata declared_com_m")
    return tuple(records), errors


def _read_csv_table(data: bytes, expected_header: list[str], label: str, *, require_final_newline: bool) -> _CsvTable:
    errors: list[str] = []
    parsed_data = data
    if require_final_newline and data and not data.endswith((b"\n", b"\r")):
        errors.append(f"{label} is truncated: final CSV row has no line terminator")
        final_newline = max(data.rfind(b"\n"), data.rfind(b"\r"))
        parsed_data = data[: final_newline + 1] if final_newline >= 0 else b""
    try:
        decoded = parsed_data.decode("utf-8")
    except UnicodeDecodeError as exc:
        return _CsvTable((), tuple(errors + [f"{label} is not UTF-8: {exc}"]))
    reader = csv.reader(io.StringIO(decoded, newline=""), strict=True)
    try:
        header = next(reader)
    except (StopIteration, csv.Error) as exc:
        return _CsvTable((), tuple(errors + [f"{label} is missing a valid CSV header: {exc}"]))
    if header != expected_header:
        return _CsvTable((), tuple(errors + [f"{label} header does not match the explicit native CSV schema"]))
    rows: list[dict[str, str]] = []
    try:
        for cells in reader:
            if not cells or all(cell == "" for cell in cells):
                continue
            if len(cells) != len(header):
                errors.append(
                    f"{label} row ending at line {reader.line_num} has {len(cells)} fields; expected {len(header)}"
                )
                break
            rows.append(dict(zip(header, cells, strict=True)))
    except csv.Error as exc:
        errors.append(f"{label} CSV parsing stopped at line {reader.line_num}: {exc}")
    return _CsvTable(tuple(rows), tuple(errors))


def _decode_point_rows(table: _CsvTable) -> tuple[tuple[NativePointStateRow, ...], list[str]]:
    rows: list[NativePointStateRow] = []
    errors: list[str] = []
    for line, row in enumerate(table.rows, start=2):
        try:
            rest_position = tuple(
                _float_field(row[key], f"point row {line} {key}") for key in ("x_ref_m", "y_ref_m", "z_ref_m")
            )
            current_position = tuple(
                _float_field(row[key], f"point row {line} {key}")
                for key in ("x_current_m", "y_current_m", "z_current_m")
            )
            displacement = tuple(
                _float_field(row[key], f"point row {line} {key}") for key in ("ux_m", "uy_m", "uz_m")
            )
            for key in (
                "vx_mps", "vy_mps", "vz_mps", "ax_mps2", "ay_mps2", "az_mps2",
                "b_int_x_N_m3", "b_int_y_N_m3", "b_int_z_N_m3",
                "b_ext_x_N_m3", "b_ext_y_N_m3", "b_ext_z_N_m3",
            ):
                _float_field(row[key], f"point row {line} {key}")
            damage = _float_field(row["damage"], f"point row {line} damage")
            if not 0.0 <= damage <= 1.0:
                raise ValueError(f"point row {line} damage must be in [0, 1]")
            active_bonds = _int_field(row["n_active_bonds"], f"point row {line} n_active_bonds", minimum=0)
            for axis, (rest, current, delta) in enumerate(zip(rest_position, current_position, displacement, strict=True)):
                expected_current = rest + delta
                if not math.isfinite(expected_current):
                    raise ValueError(f"point row {line} rest position plus displacement overflows on axis {axis}")
                scale = max(abs(rest), abs(current), abs(delta), abs(expected_current), 1e-300)
                if not _close(current, expected_current, scale):
                    raise ValueError(f"point row {line} current position disagrees with rest position plus displacement")
            rows.append(
                NativePointStateRow(
                    step=_int_field(row["step"], f"point row {line} step", minimum=0),
                    time_s=_float_field(row["time_s"], f"point row {line} time_s", nonnegative=True),
                    accepted_dt_s=_float_field(row["dt_s"], f"point row {line} dt_s", positive=True),
                    native_body_index=_int_field(row["native_body_index"], f"point row {line} native_body_index", minimum=1),
                    point_id=_int_field(row["point_id"], f"point row {line} point_id", minimum=1),
                    rest_position_m=rest_position,
                    volume_m3=_float_field(row["volume_m3"], f"point row {line} volume_m3", positive=True),
                    mass_kg=_float_field(row["mass_kg"], f"point row {line} mass_kg", positive=True),
                    n_active_bonds=active_bonds,
                    line_number=line,
                )
            )
        except ValueError as exc:
            errors.append(str(exc))
            break
    return tuple(rows), errors


def _decode_directed_rows(table: _CsvTable) -> tuple[tuple[NativeDirectedStateRow, ...], list[str]]:
    rows: list[NativeDirectedStateRow] = []
    errors: list[str] = []
    for line, row in enumerate(table.rows, start=2):
        try:
            rows.append(
                NativeDirectedStateRow(
                    step=_int_field(row["step"], f"directed row {line} step", minimum=0),
                    time_s=_float_field(row["time_s"], f"directed row {line} time_s", nonnegative=True),
                    accepted_dt_s=_float_field(row["dt_s"], f"directed row {line} dt_s", positive=True),
                    source_point_id=_int_field(row["directed_source_point_id"], f"directed row {line} source ID", minimum=1),
                    target_point_id=_int_field(row["directed_neighbor_point_id"], f"directed row {line} target ID", minimum=1),
                    canonical_point_i_id=_int_field(row["canonical_point_i_id"], f"directed row {line} canonical i", minimum=1),
                    canonical_point_j_id=_int_field(row["canonical_point_j_id"], f"directed row {line} canonical j", minimum=1),
                    active=_bool_field(row["directed_active"], f"directed row {line} directed_active"),
                    canonical_active=_bool_field(row["canonical_active"], f"directed row {line} canonical_active"),
                    native_field=f"directed_bond_state_history.csv:line-{line}:directed_active",
                    line_number=line,
                )
            )
        except ValueError as exc:
            errors.append(str(exc))
            break
    return tuple(rows), errors


def _decode_canonical_rows(table: _CsvTable) -> tuple[tuple[NativeCanonicalStateRow, ...], list[str]]:
    rows: list[NativeCanonicalStateRow] = []
    errors: list[str] = []
    for line, row in enumerate(table.rows, start=2):
        try:
            rows.append(
                NativeCanonicalStateRow(
                    step=_int_field(row["step"], f"canonical row {line} step", minimum=0),
                    time_s=_float_field(row["time_s"], f"canonical row {line} time_s", nonnegative=True),
                    accepted_dt_s=_float_field(row["dt_s"], f"canonical row {line} dt_s", positive=True),
                    point_i_id=_int_field(row["point_i_id"], f"canonical row {line} point_i_id", minimum=1),
                    point_j_id=_int_field(row["point_j_id"], f"canonical row {line} point_j_id", minimum=1),
                    active=_bool_field(row["active"], f"canonical row {line} active"),
                    line_number=line,
                )
            )
        except ValueError as exc:
            errors.append(str(exc))
            break
    return tuple(rows), errors


def _group_rows(rows: tuple[Any, ...]) -> dict[int, list[Any]]:
    groups: dict[int, list[Any]] = {}
    previous_step: int | None = None
    for row in rows:
        step = row.step
        if previous_step is not None and step < previous_step:
            # Preserve grouped data for diagnosis; sequence validation below
            # will make the snapshot incomplete instead of reordering evidence.
            pass
        groups.setdefault(step, []).append(row)
        previous_step = step
    return groups


def _check_row_step_order(rows: tuple[Any, ...], label: str, reasons: list[str]) -> None:
    observed: list[int] = []
    for row in rows:
        if not observed or row.step != observed[-1]:
            observed.append(row.step)
    if not observed:
        return
    if observed[0] != 0 or observed != list(range(observed[-1] + 1)):
        reasons.append(f"{label} rows are not a contiguous step-zero prefix: {observed}")


def _clock_group_valid(step: int, clocks: list[tuple[float, float]], expected_dt: float | None) -> bool:
    if not clocks:
        return False
    first_time, first_dt = clocks[0]
    if any(time_s != first_time or dt_s != first_dt for time_s, dt_s in clocks[1:]):
        return False
    if step == 0 and first_time != 0.0:
        return False
    expected_time = _finite_product(step, first_dt)
    if expected_time is None or not _close(first_time, expected_time, first_dt):
        return False
    return expected_dt is None or _close(first_dt, expected_dt, expected_dt)


def _summary_int(summary: Mapping[str, Any], key: str, reasons: list[str], *, minimum: int) -> int | None:
    value = summary.get(key)
    if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
        reasons.append(f"summary field {key} must be an integer >= {minimum}")
        return None
    return value


def _summary_float(summary: Mapping[str, Any], key: str, reasons: list[str], *, positive: bool | None = None) -> float | None:
    value = summary.get(key)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        reasons.append(f"summary field {key} must be finite numeric data")
        return None
    try:
        result = float(value)
    except (OverflowError, TypeError, ValueError):
        reasons.append(f"summary field {key} must be finite numeric data")
        return None
    if not math.isfinite(result):
        reasons.append(f"summary field {key} must be finite numeric data")
        return None
    if positive is True and result <= 0.0:
        reasons.append(f"summary field {key} must be positive")
        return None
    if positive is False and result < 0.0:
        reasons.append(f"summary field {key} must be non-negative")
        return None
    return result


def _check_summary_count(summary: Mapping[str, Any], key: str, actual: int, reasons: list[str]) -> None:
    value = summary.get(key)
    if not isinstance(value, int) or isinstance(value, bool) or value != actual:
        reasons.append(f"summary field {key} disagrees with parsed rows (expected {actual}, declared {value!r})")


def _int_field(value: str, label: str, *, minimum: int) -> int:
    if re.fullmatch(r"[0-9]+", value) is None:
        raise ValueError(f"{label} must be an integer")
    result = int(value)
    if result < minimum:
        raise ValueError(f"{label} must be >= {minimum}")
    return result


def _float_field(value: str, label: str, *, positive: bool = False, nonnegative: bool = False) -> float:
    try:
        result = float(value)
    except ValueError as exc:
        raise ValueError(f"{label} must be numeric") from exc
    if not math.isfinite(result):
        raise ValueError(f"{label} must be finite")
    if positive and result <= 0.0:
        raise ValueError(f"{label} must be positive")
    if nonnegative and result < 0.0:
        raise ValueError(f"{label} must be non-negative")
    return result


def _bool_field(value: str, label: str) -> bool:
    if value == "true":
        return True
    if value == "false":
        return False
    raise ValueError(f"{label} must be the exact lowercase string 'true' or 'false'")


def _close(actual: float, expected: float, scale: float) -> bool:
    if not all(math.isfinite(value) for value in (actual, expected, scale)):
        return False
    magnitude = max(abs(actual), abs(expected), abs(scale), 1e-300)
    return abs(actual - expected) <= max(1e-10 * abs(scale), 32.0 * math.ulp(magnitude))


def _finite_sum(values: Any, label: str, errors: list[str]) -> float | None:
    try:
        total = math.fsum(values)
    except OverflowError:
        errors.append(f"{label} overflowed")
        return None
    if not math.isfinite(total):
        errors.append(f"{label} is not finite")
        return None
    return total


def _finite_product(left: Any, right: Any) -> float | None:
    try:
        left_float = float(left)
        right_float = float(right)
    except (OverflowError, TypeError, ValueError):
        return None
    if not math.isfinite(left_float) or not math.isfinite(right_float):
        return None
    product = left_float * right_float
    return product if math.isfinite(product) else None


def _is_sha256(value: Any) -> bool:
    return isinstance(value, str) and re.fullmatch(r"[0-9a-fA-F]{64}", value) is not None


def _nonempty_string(value: Any) -> str | None:
    return value.strip() if isinstance(value, str) and value.strip() else None
