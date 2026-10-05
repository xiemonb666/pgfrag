"""Run one bounded local-front provider check; no operator or solver."""

from __future__ import annotations

import json
import math
from pathlib import Path
import sys
import time
from typing import Any

import numpy as np
from PIL import Image, ImageDraw, ImageFont

from local_front_provider import (
    MAX_BATCH_QUERIES,
    PINNED_MANIFEST_PAYLOAD_SHA256,
    REPO_ROOT,
    TraceSide,
    LocalFrontProvider,
    MLSStatus,
    build_registry,
    exact_modeiii_uz_and_gradient,
    load_pinned_inputs,
    modeiii_amplitude,
    normalized_f1_and_gradient,
    pinned_source_hashes,
    topology_box_clipped_endpoint_control,
    topology_negative_controls,
)


EXPERIMENT_DIR = Path(__file__).resolve().parent
OUTPUT_JSON = EXPERIMENT_DIR / "provider_diagnostic.json"
OUTPUT_PNG = EXPERIMENT_DIR / "local_front_section.png"
PROVIDER_SOURCE = EXPERIMENT_DIR / "local_front_provider.py"
WALL_LIMIT_S = 30.0


def _finite(value: float) -> float | None:
    return float(value) if math.isfinite(float(value)) else None


def _add_query(
    points: list[tuple[float, float, float]],
    sides: list[TraceSide | None],
    metadata: list[dict[str, Any]],
    label: str,
    point: tuple[float, float, float],
    side: TraceSide | None = None,
    kind: str = "special",
) -> int:
    index = len(points)
    points.append(tuple(float(value) for value in point))
    sides.append(side)
    metadata.append(
        {
            "query_id": label,
            "query_index": index,
            "point_m": [float(value) for value in point],
            "trace_side": None if side is None else side.value,
            "kind": kind,
        }
    )
    return index


def _build_queries() -> tuple[
    list[tuple[float, float, float]],
    list[TraceSide | None],
    list[dict[str, Any]],
    list[float],
    list[float],
    list[list[int]],
    dict[str, int],
    dict[str, int],
    dict[str, int],
]:
    points: list[tuple[float, float, float]] = []
    sides: list[TraceSide | None] = []
    metadata: list[dict[str, Any]] = []
    x_values = np.delete(np.linspace(-0.009, 0.009, 15), 7)
    y_values = np.delete(np.linspace(-0.0096, 0.0096, 17), 8)
    grid_indices: list[list[int]] = []
    for y_index, y in enumerate(y_values):
        row: list[int] = []
        for x_index, x in enumerate(x_values):
            row.append(
                _add_query(
                    points,
                    sides,
                    metadata,
                    f"section_y{y_index:02d}_x{x_index:02d}",
                    (float(x), float(y), 0.0),
                    kind="section_grid",
                )
            )
        grid_indices.append(row)

    core_indices: dict[str, int] = {}
    for sx in (-1, 1):
        for sy in (-1, 1):
            label = f"core_{'p' if sx > 0 else 'm'}x_{'p' if sy > 0 else 'm'}y"
            core_indices[label] = _add_query(
                points,
                sides,
                metadata,
                label,
                (sx * 0.00045, sy * 0.00045, 0.0),
            )

    trace_indices: dict[str, int] = {}
    for x in (-0.008, -0.006, -0.004):
        x_label = f"{x:.3f}"
        for side in (TraceSide.PLUS, TraceSide.MINUS):
            label = f"gamma_x{x_label}_{side.value.lower()}"
            trace_indices[label] = _add_query(
                points,
                sides,
                metadata,
                label,
                (x, 0.0, 0.0),
                side=side,
                kind="finite_gamma_trace",
            )

    ahead_index = _add_query(
        points,
        sides,
        metadata,
        "ahead_of_front_plane_no_trace",
        (0.004, 0.0, 0.0),
        kind="front_ghost_check",
    )
    negative_indices = {
        "front_core": _add_query(points, sides, metadata, "front_core", (0.0, 0.0, 0.0)),
        "finite_front_endpoint_plus": _add_query(
            points, sides, metadata, "finite_front_endpoint_plus", (0.0, 0.0, 0.008)
        ),
        "crack_face_missing_trace": _add_query(
            points, sides, metadata, "crack_face_missing_trace", (-0.006, 0.0, 0.0)
        ),
        "ambiguous_cut_tolerance_band": _add_query(
            points, sides, metadata, "ambiguous_cut_tolerance_band", (-0.006, 0.5e-12, 0.0)
        ),
        "finite_crack_end": _add_query(
            points,
            sides,
            metadata,
            "finite_crack_end",
            (-0.006, 0.0, 0.008),
            side=TraceSide.PLUS,
        ),
        "out_of_domain": _add_query(
            points, sides, metadata, "out_of_domain", (0.0101, 0.002, 0.0)
        ),
    }
    # Exactly one source lies at this support-radius distance. The MLS kernel
    # switch guard must mark this query UNKNOWN instead of smoothing across it.
    negative_indices["tangent_mls_support"] = _add_query(
        points,
        sides,
        metadata,
        "tangent_mls_support",
        (-0.001, -0.008, -0.008),
    )
    if len(points) != 242 or len(points) > MAX_BATCH_QUERIES:
        raise RuntimeError(f"unexpected focused-check query count: {len(points)}")
    _ = ahead_index
    return points, sides, metadata, list(map(float, x_values)), list(map(float, y_values)), grid_indices, core_indices, trace_indices, negative_indices


def _coefficient_vector(provider: LocalFrontProvider, amplitude: float) -> np.ndarray:
    source_by_id = {source.source_id: source for source in provider.sources}
    front_coefficient = amplitude * math.sqrt(0.01)
    coefficients = np.zeros(provider.registry.dof_count, dtype=np.float64)
    for slot, key in enumerate(provider.registry.ordered_keys):
        source = source_by_id[key.source_id]
        if key.basis_kind.value == "base" and source.geometric_role == "FRONT":
            source_fbar, _ = normalized_f1_and_gradient(source.position_m, None)
            coefficients[slot] = front_coefficient * source_fbar
        elif key.basis_kind.value == "front_f1_shifted":
            coefficients[slot] = front_coefficient
        # BASE coefficients outside FRONT and all HEAVISIDE coefficients are 0.
    return coefficients


def _scalar_row(provider: LocalFrontProvider, query_index: int, coefficients: np.ndarray) -> tuple[float, tuple[float, float, float]]:
    arrays = provider.row_arrays(query_index)
    if arrays is None:
        raise ValueError(f"query {query_index} does not have an OK provider row")
    values, gradients = arrays
    value = float(values @ coefficients)
    gradient = gradients.T @ coefficients
    return value, tuple(float(component) for component in gradient)


def _build_section_data(
    provider: LocalFrontProvider,
    x_values: list[float],
    y_values: list[float],
    grid_indices: list[list[int]],
    coefficients: np.ndarray,
    amplitude: float,
) -> dict[str, Any]:
    exact_matrix: list[list[float | None]] = []
    local_matrix: list[list[float | None]] = []
    error_matrix: list[list[float | None]] = []
    unknown_cells: list[dict[str, Any]] = []
    successful_errors: list[float] = []
    for j, y in enumerate(y_values):
        exact_row: list[float | None] = []
        local_row: list[float | None] = []
        error_row: list[float | None] = []
        for i, x in enumerate(x_values):
            query_index = grid_indices[j][i]
            if provider.row_statuses[query_index].value != "OK":
                exact_row.append(None)
                local_row.append(None)
                error_row.append(None)
                unknown_cells.append(
                    {
                        "query_index": query_index,
                        "query_id": f"section_y{j:02d}_x{i:02d}",
                        "point_m": list(provider.query_points_m[query_index]),
                        "provider_status": provider.row_statuses[query_index].value,
                        "reason": provider.row_diagnostics[query_index],
                    }
                )
                continue
            point = provider.query_points_m[query_index]
            exact, _ = exact_modeiii_uz_and_gradient(point, None, amplitude)
            local, _ = _scalar_row(provider, query_index, coefficients)
            error = abs(local - exact)
            exact_row.append(exact)
            local_row.append(local)
            error_row.append(error)
            successful_errors.append(error)
        exact_matrix.append(exact_row)
        local_matrix.append(local_row)
        error_matrix.append(error_row)
    error_array = np.asarray(successful_errors, dtype=np.float64)
    summary = {
        "sampled_grid_shape_yx": [len(y_values), len(x_values)],
        "successful_grid_query_count": int(error_array.size),
        "unknown_grid_query_count": len(unknown_cells),
        "mean_absolute_error_m": float(np.mean(error_array)) if error_array.size else None,
        "rms_error_m": float(np.sqrt(np.mean(error_array**2))) if error_array.size else None,
        "p95_absolute_error_m": float(np.percentile(error_array, 95.0)) if error_array.size else None,
        "maximum_absolute_error_m": float(np.max(error_array)) if error_array.size else None,
        "interpretation": "whole-box z=0 section sample; collar error is reported data, not an exactness target outside the r<1 mm front core",
    }
    return {
        "z_m": 0.0,
        "x_m": x_values,
        "y_m": y_values,
        "exact_uz_m": exact_matrix,
        "local_reconstructed_uz_m": local_matrix,
        "absolute_error_m": error_matrix,
        "unknown_cells": unknown_cells,
        "summary": summary,
    }


def _font(size: int) -> ImageFont.ImageFont:
    try:
        return ImageFont.truetype("arial.ttf", size=size)
    except OSError:
        return ImageFont.load_default()


def _diverging_color(value: float, limit: float) -> tuple[int, int, int]:
    t = min(1.0, max(0.0, 0.5 + 0.5 * value / limit))
    if t <= 0.5:
        u = 2.0 * t
        return (int(42 + 213 * u), int(89 + 166 * u), int(172 + 83 * u))
    u = 2.0 * (t - 0.5)
    return (255, int(255 - 178 * u), int(255 - 214 * u))


def _error_color(value: float, limit: float) -> tuple[int, int, int]:
    t = min(1.0, max(0.0, value / limit))
    return (int(35 + 215 * t), int(79 + 95 * (1.0 - t)), int(145 - 95 * t))


def _edges(values: list[float], lower: float, upper: float) -> list[float]:
    result = [lower]
    result.extend(0.5 * (values[index] + values[index + 1]) for index in range(len(values) - 1))
    result.append(upper)
    return result


def _draw_panel(
    draw: ImageDraw.ImageDraw,
    origin_x: int,
    title: str,
    matrix: list[list[float | None]],
    x_values: list[float],
    y_values: list[float],
    *,
    diverging: bool,
    value_limit: float,
    font_title: ImageFont.ImageFont,
    font_small: ImageFont.ImageFont,
) -> None:
    left = origin_x + 58
    top = 92
    size = 330
    right = left + size
    bottom = top + size
    draw.text((origin_x + 10, 28), title, fill=(20, 30, 45), font=font_title)
    x_edges = _edges(x_values, -0.01, 0.01)
    y_edges = _edges(y_values, -0.01, 0.01)
    for j, y in enumerate(y_values):
        for i, x in enumerate(x_values):
            value = matrix[j][i]
            if value is None:
                color = (205, 207, 211)
            elif diverging:
                color = _diverging_color(float(value), value_limit)
            else:
                color = _error_color(float(value), value_limit)
            px0 = left + int((x_edges[i] + 0.01) / 0.02 * size)
            px1 = left + int((x_edges[i + 1] + 0.01) / 0.02 * size)
            py0 = top + int((0.01 - y_edges[j + 1]) / 0.02 * size)
            py1 = top + int((0.01 - y_edges[j]) / 0.02 * size)
            draw.rectangle((px0, py0, max(px0, px1 - 1), max(py0, py1 - 1)), fill=color)
    draw.rectangle((left, top, right, bottom), outline=(35, 42, 54), width=2)
    # Mark the branch line and front without interpolating field values onto it.
    crack_y = top + size // 2
    front_x = left + size // 2
    draw.line((left, crack_y, front_x, crack_y), fill=(30, 30, 30), width=2)
    draw.line((front_x, top, front_x, bottom), fill=(75, 75, 75), width=1)
    draw.ellipse((front_x - 4, crack_y - 4, front_x + 4, crack_y + 4), fill=(0, 0, 0))
    draw.text((left - 12, bottom + 9), "-0.01", fill=(30, 30, 30), font=font_small)
    draw.text((front_x - 4, bottom + 9), "0", fill=(30, 30, 30), font=font_small)
    draw.text((right - 35, bottom + 9), "+0.01", fill=(30, 30, 30), font=font_small)
    draw.text((left + size // 2 - 8, bottom + 28), "x (m)", fill=(30, 30, 30), font=font_small)
    draw.text((origin_x + 4, top - 2), "+y", fill=(30, 30, 30), font=font_small)
    draw.text((origin_x + 4, bottom - 10), "-y", fill=(30, 30, 30), font=font_small)

    bar_left = right + 12
    bar_right = bar_left + 18
    for py in range(top, bottom):
        fraction = 1.0 - (py - top) / max(1, size - 1)
        sample_value = (fraction * 2.0 - 1.0) * value_limit if diverging else fraction * value_limit
        color = _diverging_color(sample_value, value_limit) if diverging else _error_color(sample_value, value_limit)
        draw.line((bar_left, py, bar_right, py), fill=color)
    draw.rectangle((bar_left, top, bar_right, bottom), outline=(35, 42, 54), width=1)
    draw.text((bar_right + 4, top - 2), f"{value_limit:.1e}", fill=(30, 30, 30), font=font_small)
    if diverging:
        draw.text((bar_right + 4, top + size // 2 - 7), "0", fill=(30, 30, 30), font=font_small)
        draw.text((bar_right + 4, bottom - 9), f"-{value_limit:.1e}", fill=(30, 30, 30), font=font_small)
    else:
        draw.text((bar_right + 4, bottom - 9), "0", fill=(30, 30, 30), font=font_small)


def _write_section_png(section: dict[str, Any], output_path: Path) -> None:
    x_values = [float(value) for value in section["x_m"]]
    y_values = [float(value) for value in section["y_m"]]
    matrices = [section["exact_uz_m"], section["local_reconstructed_uz_m"], section["absolute_error_m"]]
    direct_values = [float(value) for matrix in matrices[:2] for row in matrix for value in row if value is not None]
    errors = [float(value) for row in matrices[2] for value in row if value is not None]
    uz_limit = max((abs(value) for value in direct_values), default=1.0e-12)
    error_limit = max(errors, default=1.0e-12)
    uz_limit = max(uz_limit, 1.0e-12)
    error_limit = max(error_limit, 1.0e-12)
    image = Image.new("RGB", (1480, 520), (255, 255, 255))
    draw = ImageDraw.Draw(image)
    title_font = _font(19)
    small_font = _font(13)
    _draw_panel(draw, 5, "Exact Mode III uz (m)", matrices[0], x_values, y_values, diverging=True,
                value_limit=uz_limit, font_title=title_font, font_small=small_font)
    _draw_panel(draw, 495, "Local shifted-F1 PU uz (m)", matrices[1], x_values, y_values, diverging=True,
                value_limit=uz_limit, font_title=title_font, font_small=small_font)
    _draw_panel(draw, 985, "Absolute error (m)", matrices[2], x_values, y_values, diverging=False,
                value_limit=error_limit, font_title=title_font, font_small=small_font)
    draw.text(
        (20, 498),
        "Direct z=0 MLS samples on both open sides; black line marks finite Gamma, grey line marks x=0. Trace values are checked separately.",
        fill=(45, 52, 62),
        font=small_font,
    )
    image.save(output_path, format="PNG", optimize=True)


def _round_registry_order_check(bundle: Any, provider: LocalFrontProvider) -> dict[str, Any]:
    reverse_registry = build_registry(bundle.sources, reverse_registration=True)
    same = reverse_registry.canonical_ids == provider.registry.canonical_ids
    return {
        "reverse_registration_canonical_ids_unchanged": same,
        "canonical_id_count": len(provider.registry.canonical_ids),
        "all_canonical_ids_unique": len(set(provider.registry.canonical_ids)) == provider.registry.dof_count,
        "provider_ids_equal_registry_ids": provider.coefficient_ids == provider.registry.canonical_ids,
    }


def run() -> dict[str, Any]:
    started = time.perf_counter()
    bundle = load_pinned_inputs()
    points, sides, metadata, x_values, y_values, grid_indices, core_indices, trace_indices, negative_indices = _build_queries()
    provider = LocalFrontProvider(bundle, points, sides)
    amplitude = modeiii_amplitude(bundle.field_material)
    front_coefficient = amplitude * math.sqrt(0.01)
    coefficients = _coefficient_vector(provider, amplitude)

    row_validation = []
    for query_index, sample in enumerate(provider.sample_at(index) for index in range(len(points))):
        row_validation.append(
            {
                "query_index": query_index,
                "status": sample.status.value,
                "dense_id_count": len(sample.coefficient_ids),
                "dense_value_count": len(sample.contributions),
                "diagnostic": sample.diagnostic,
            }
        )
    contract_rows_ok = all(
        (row["status"] == "OK" and row["dense_id_count"] == 140 and row["dense_value_count"] == 140)
        or (row["status"] == "UNKNOWN" and row["dense_id_count"] == 0 and row["dense_value_count"] == 0)
        for row in row_validation
    )

    affine_intercept = 0.125
    affine_gradient = np.asarray((-0.7, 1.2, 0.25), dtype=np.float64)
    affine_coefficients = np.asarray(
        [affine_intercept + float(np.dot(affine_gradient, source.position_m)) for source in provider.sources],
        dtype=np.float64,
    )
    base_slots = {
        source.source_id: provider.registry.slot_of(
            next(
                key for key in provider.registry.ordered_keys
                if key.source_id == source.source_id and key.basis_kind.value == "base"
            )
        )
        for source in provider.sources
    }
    affine_value_errors: list[float] = []
    affine_gradient_errors: list[float] = []
    for query_index, point in enumerate(provider.query_points_m):
        arrays = provider.row_arrays(query_index)
        if arrays is None:
            continue
        values, gradients = arrays
        slots = [base_slots[source.source_id] for source in provider.sources]
        approx = float(values[slots] @ affine_coefficients)
        exact = affine_intercept + float(np.dot(affine_gradient, point))
        grad_approx = gradients[slots].T @ affine_coefficients
        affine_value_errors.append(abs(approx - exact))
        affine_gradient_errors.append(float(np.max(np.abs(grad_approx - affine_gradient))))

    core_rows: list[dict[str, Any]] = []
    core_field_errors: list[float] = []
    core_gradient_errors: list[float] = []
    core_partition_sums: list[float] = []
    front_ids = {source.source_id for source in provider.sources if source.geometric_role == "FRONT"}
    front_base_slots = [base_slots[source_id] for source_id in front_ids]
    for label, query_index in core_indices.items():
        point = provider.query_points_m[query_index]
        radius = math.hypot(point[0], point[1])
        local, local_gradient = _scalar_row(provider, query_index, coefficients)
        exact, exact_gradient = exact_modeiii_uz_and_gradient(point, None, amplitude)
        field_error = abs(local - exact)
        gradient_error = max(abs(a - b) for a, b in zip(local_gradient, exact_gradient))
        values, _gradients = provider.row_arrays(query_index)  # type: ignore[misc]
        partition_sum = float(np.sum(values[front_base_slots]))
        core_field_errors.append(field_error)
        core_gradient_errors.append(gradient_error)
        core_partition_sums.append(partition_sum)
        core_rows.append(
            {
                "query_id": label,
                "point_m": list(point),
                "radius_m": radius,
                "inside_core_radius": radius < 0.001,
                "front_base_partition_sum_observed_only": partition_sum,
                "exact_uz_m": exact,
                "local_uz_m": local,
                "absolute_field_error_m": field_error,
                "exact_gradient_m_per_m": list(exact_gradient),
                "local_gradient_m_per_m": list(local_gradient),
                "maximum_gradient_error_m_per_m": gradient_error,
            }
        )

    trace_rows: list[dict[str, Any]] = []
    max_base_jump = 0.0
    max_front_jump_identity_error = 0.0
    max_heaviside_jump_identity_error = 0.0
    for x in (-0.008, -0.006, -0.004):
        plus_index = trace_indices[f"gamma_x{x:.3f}_plus"]
        minus_index = trace_indices[f"gamma_x{x:.3f}_minus"]
        plus_values, _ = provider.row_arrays(plus_index)  # type: ignore[misc]
        minus_values, _ = provider.row_arrays(minus_index)  # type: ignore[misc]
        fbar_plus, _ = normalized_f1_and_gradient((x, 0.0, 0.0), TraceSide.PLUS)
        fbar_minus, _ = normalized_f1_and_gradient((x, 0.0, 0.0), TraceSide.MINUS)
        fbar_jump = fbar_plus - fbar_minus
        h_entries: list[dict[str, Any]] = []
        f_entries: list[dict[str, Any]] = []
        h_jump = 0.0
        f_jump = 0.0
        for slot, key in enumerate(provider.registry.ordered_keys):
            difference = float(plus_values[slot] - minus_values[slot])
            if key.basis_kind.value == "base":
                max_base_jump = max(max_base_jump, abs(difference))
            elif key.basis_kind.value == "heaviside":
                phi = float(plus_values[base_slots[key.source_id]])
                expected = 2.0 * phi
                max_heaviside_jump_identity_error = max(max_heaviside_jump_identity_error, abs(difference - expected))
                h_jump += difference
                h_entries.append(
                    {"source_id": key.source_id, "phi": phi, "jump_plus_minus": difference, "expected_2phi": expected}
                )
            elif key.basis_kind.value == "front_f1_shifted":
                phi = float(plus_values[base_slots[key.source_id]])
                expected = phi * fbar_jump
                max_front_jump_identity_error = max(max_front_jump_identity_error, abs(difference - expected))
                f_jump += difference
                f_entries.append(
                    {"source_id": key.source_id, "phi": phi, "jump_plus_minus": difference, "expected_phi_times_Fbar_jump": expected}
                )
        exact_plus, _ = exact_modeiii_uz_and_gradient((x, 0.0, 0.0), TraceSide.PLUS, amplitude)
        exact_minus, _ = exact_modeiii_uz_and_gradient((x, 0.0, 0.0), TraceSide.MINUS, amplitude)
        local_plus, _ = _scalar_row(provider, plus_index, coefficients)
        local_minus, _ = _scalar_row(provider, minus_index, coefficients)
        trace_rows.append(
            {
                "x_m": x,
                "z_m": 0.0,
                "exact_gamma_plus_fbar": fbar_plus,
                "exact_gamma_minus_fbar": fbar_minus,
                "exact_uz_plus_m": exact_plus,
                "exact_uz_minus_m": exact_minus,
                "local_uz_plus_m": local_plus,
                "local_uz_minus_m": local_minus,
                "base_row_plus_minus_max_abs": max(
                    abs(float(plus_values[base_slots[source.source_id]] - minus_values[base_slots[source.source_id]]))
                    for source in provider.sources
                ),
                "heaviside_all_one_coeff_jump_m": h_jump,
                "heaviside_jump_by_source": h_entries,
                "front_all_one_coeff_jump_m": f_jump,
                "front_jump_by_source": f_entries,
            }
        )

    ahead_index = next(
        index for index, row in enumerate(metadata) if row["query_id"] == "ahead_of_front_plane_no_trace"
    )
    ahead_values, ahead_gradients = provider.row_arrays(ahead_index)  # type: ignore[misc]
    heaviside_slots = [
        slot for slot, key in enumerate(provider.registry.ordered_keys) if key.basis_kind.value == "heaviside"
    ]
    ahead_h_value_max = float(np.max(np.abs(ahead_values[heaviside_slots]))) if heaviside_slots else 0.0
    ahead_h_gradient_max = float(np.max(np.abs(ahead_gradients[heaviside_slots]))) if heaviside_slots else 0.0

    section = _build_section_data(provider, x_values, y_values, grid_indices, coefficients, amplitude)
    topology_controls = topology_negative_controls(bundle)
    box_clipped_endpoint = topology_box_clipped_endpoint_control(bundle)
    negative_query_rows = []
    expected_negative_reasons = {
        "front_core": "front_or_front_endpoint",
        "finite_front_endpoint_plus": "front_or_front_endpoint",
        "crack_face_missing_trace": "explicit_plus_minus_trace_required",
        "ambiguous_cut_tolerance_band": "ambiguous_branch_tolerance_band",
        "finite_crack_end": "finite_crack_boundary_or_end",
        "out_of_domain": "outside_domain",
        "tangent_mls_support": "MLS:UNKNOWN_SUPPORT_BOUNDARY",
    }
    for label, query_index in negative_indices.items():
        actual = provider.row_diagnostics[query_index]
        status = provider.row_statuses[query_index].value
        expected = expected_negative_reasons[label]
        negative_query_rows.append(
            {
                "query_id": label,
                "point_m": list(provider.query_points_m[query_index]),
                "provider_status": status,
                "actual_reason": actual,
                "expected_reason": expected,
                "passed": status == "UNKNOWN" and actual == expected,
            }
        )

    registry_order_check = _round_registry_order_check(bundle, provider)
    condition_values = []
    support_values = []
    if provider.mls_result is not None:
        condition_values = [
            float(value) for value in provider.mls_result.condition_numbers if math.isfinite(float(value))
        ]
        support_values = [int(value) for value in provider.mls_result.support_counts]
    unknown_rows = [
        {
            "query_index": index,
            "query_id": metadata[index]["query_id"],
            "point_m": list(provider.query_points_m[index]),
            "status": provider.row_statuses[index].value,
            "reason": provider.row_diagnostics[index],
        }
        for index in range(len(points))
        if provider.row_statuses[index].value != "OK"
    ]
    source_role_binding = {
        "source_count": len(bundle.sources),
        "unique_source_id_count": len({source.source_id for source in bundle.sources}),
        "basis_role_counts": {
            role: sum(source.geometric_role == role for source in bundle.sources)
            for role in ("BASE", "FRONT", "HEAVISIDE_ELIGIBLE")
        },
        "base_keys": sum(key.basis_kind.value == "base" for key in provider.registry.ordered_keys),
        "front_f1_shifted_keys": sum(key.basis_kind.value == "front_f1_shifted" for key in provider.registry.ordered_keys),
        "heaviside_keys": sum(key.basis_kind.value == "heaviside" for key in provider.registry.ordered_keys),
        "maximum_enrichments_per_source": max(len(source.basis_roles) - 1 for source in bundle.sources),
        "all_sources_retained": len(provider.source_positions) == 100 and provider.mls_result is not None and provider.mls_result.weights.shape[1] == 100,
    }
    all_checks_passed = (
        contract_rows_ok
        and registry_order_check["reverse_registration_canonical_ids_unchanged"]
        and registry_order_check["all_canonical_ids_unique"]
        and source_role_binding["source_count"] == 100
        and source_role_binding["base_keys"] == 100
        and source_role_binding["front_f1_shifted_keys"] == 30
        and source_role_binding["heaviside_keys"] == 10
        and source_role_binding["maximum_enrichments_per_source"] == 1
        and source_role_binding["all_sources_retained"]
        and bool(affine_value_errors)
        and max(affine_value_errors) < 1.0e-9
        and max(affine_gradient_errors) < 1.0e-7
        and len(core_rows) == 4
        and max(core_field_errors) < 1.0e-10
        and max(core_gradient_errors) < 1.0e-7
        and max_base_jump < 1.0e-14
        and max_front_jump_identity_error < 1.0e-14
        and max_heaviside_jump_identity_error < 1.0e-14
        and ahead_h_value_max == 0.0
        and ahead_h_gradient_max == 0.0
        and all(row["passed"] for row in negative_query_rows)
        and all(row["passed"] for row in topology_controls)
        and box_clipped_endpoint["passed"]
        and section["summary"]["successful_grid_query_count"] > 0
    )
    elapsed = time.perf_counter() - started
    report: dict[str, Any] = {
        "schema": "local-front-provider-diagnostic.v1",
        "status": "PASS" if all_checks_passed and elapsed <= WALL_LIMIT_S else "FAIL",
        "scope": "source-only continuous-MLS provider; no operator, stiffness assembly, solver, CAD, or fracture evolution",
        "pinned_inputs": {
            "raw_sha256": pinned_source_hashes(),
            "manifest_payload_sha256": PINNED_MANIFEST_PAYLOAD_SHA256,
            "case_id": bundle.case["case_id"],
            "units": "m",
        },
        "implementation_sha256": {
            "provider": __import__("hashlib").sha256(PROVIDER_SOURCE.read_bytes()).hexdigest(),
            "diagnostic_runner": __import__("hashlib").sha256(Path(__file__).read_bytes()).hexdigest(),
        },
        "budgets": {
            "source_count": len(bundle.sources),
            "scalar_dof_count": provider.registry.dof_count,
            "query_count_in_frozen_batch": len(points),
            "mls_query_count_after_fail_closed_geometry_preflight": len(provider._mls_row_by_query),
            "maximum_queries_per_batch": MAX_BATCH_QUERIES,
            "within_source_cap_125": len(bundle.sources) <= 125,
            "within_scalar_dof_cap_375": provider.registry.dof_count <= 375,
            "within_query_cap_256": len(points) <= 256,
            "visibility": "continuous MLS; no visibility mask and no source deletion",
        },
        "roles_and_semantic_ids": {
            **source_role_binding,
            **registry_order_check,
            "canonical_ids_first_6": list(provider.registry.canonical_ids[:6]),
            "canonical_ids_last_6": list(provider.registry.canonical_ids[-6:]),
        },
        "provider_contract_validation": {
            "rows_checked": len(row_validation),
            "ok_rows": provider.ok_query_count,
            "unknown_rows": provider.unknown_query_count,
            "all_rows_were_dense_or_empty_unknown": contract_rows_ok,
            "row_summaries": row_validation,
        },
        "mls_diagnostics": {
            "reference_origin_m": list(bundle.mls_config.reference_origin),
            "reference_scale_m": bundle.mls_config.reference_scale,
            "support_radius_h_m": bundle.mls_config.support_radius,
            "kernel": "Wendland C2",
            "basis": "linear 3-D MLS",
            "input_source_order": "role manifest sequence; source_id binds every canonical coefficient key",
            "status_counts": {
                status.value: sum(
                    item.value == status.value for item in (provider.mls_result.statuses if provider.mls_result else ())
                )
                for status in MLSStatus
            },
            "support_count_min": min(support_values) if support_values else None,
            "support_count_max": max(support_values) if support_values else None,
            "condition_number_min": min(condition_values) if condition_values else None,
            "condition_number_max": max(condition_values) if condition_values else None,
        },
        "shape_formulas": {
            "BASE": "phi_i",
            "HEAVISIDE": "phi_i*(H(q)-H(X_i)); H exists only for manifest-certified fully split supports; no source deletion",
            "FRONT_F1_SHIFTED": "phi_i*(Fbar(q)-Fbar(X_i)); Fbar=F1/sqrt(R0), R0=0.01 m",
            "front_gradient": "(Fbar(q)-Fbar(X_i))*grad(phi_i)+phi_i*grad(Fbar(q))",
            "heaviside_gradient_on_each_open_side": "(H(q)-H(X_i))*grad(phi_i); the distributional jump is represented by exact sided rows, not an epsilon offset",
            "partition_normalization": "none; signed MLS weights are not divided by any enriched partition sum",
            "finite_gamma_traces": "PLUS/MINUS evaluated at the exact y=0 coordinate; no epsilon displacement",
        },
        "affine_base_reproduction": {
            "test_field": "0.125 + (-0.7)*x + 1.2*y + 0.25*z",
            "successful_rows_checked": len(affine_value_errors),
            "maximum_value_error": max(affine_value_errors) if affine_value_errors else None,
            "maximum_gradient_component_error": max(affine_gradient_errors) if affine_gradient_errors else None,
            "passed": bool(affine_value_errors) and max(affine_value_errors) < 1.0e-9 and max(affine_gradient_errors) < 1.0e-7,
        },
        "core_shifted_f1_reproduction": {
            "front_core_radius_m": 0.001,
            "exact_target": "u=A*sqrt(r)*sin(theta/2)",
            "coefficient_choice": "BASE_i=c*Fbar(X_i) on T_F; FRONT_i=c; H_i=0; c=A*sqrt(R0)",
            "amplitude_A_sqrt_m": amplitude,
            "front_coefficient_c_m": front_coefficient,
            "no_partition_normalization": True,
            "queries": core_rows,
            "maximum_field_error_m": max(core_field_errors) if core_field_errors else None,
            "maximum_gradient_error_m_per_m": max(core_gradient_errors) if core_gradient_errors else None,
            "observed_front_base_partition_sums": core_partition_sums,
            "passed": len(core_rows) == 4 and max(core_field_errors) < 1.0e-10 and max(core_gradient_errors) < 1.0e-7,
        },
        "finite_gamma_trace_checks": {
            "traces": trace_rows,
            "maximum_base_plus_minus_jump": max_base_jump,
            "maximum_front_shifted_jump_identity_error": max_front_jump_identity_error,
            "maximum_heaviside_jump_identity_error": max_heaviside_jump_identity_error,
            "passed_trace_identities": max_base_jump < 1.0e-14 and max_front_jump_identity_error < 1.0e-14 and max_heaviside_jump_identity_error < 1.0e-14,
        },
        "no_infinite_plane_ghost_ahead_of_front": {
            "query_id": "ahead_of_front_plane_no_trace",
            "point_m": list(provider.query_points_m[ahead_index]),
            "provider_status": provider.row_statuses[ahead_index].value,
            "trace_side_required": False,
            "maximum_absolute_heaviside_value": ahead_h_value_max,
            "maximum_absolute_heaviside_gradient": ahead_h_gradient_max,
            "passed": provider.row_statuses[ahead_index].value == "OK" and ahead_h_value_max == 0.0 and ahead_h_gradient_max == 0.0,
        },
        "negative_query_cases": negative_query_rows,
        "finite_support_topology_negative_controls": topology_controls,
        "box_clipped_endpoint_topology_control": box_clipped_endpoint,
        "whole_box_z0_section": section,
        "unknown_queries": unknown_rows,
        "wall_clock_seconds": elapsed,
        "wall_clock_limit_seconds": WALL_LIMIT_S,
        "within_wall_clock_limit": elapsed <= WALL_LIMIT_S,
        "known_generalization_limit": "The provider rejects changed case/classifier/manifest hashes and does not claim the pinned classifier is a proved generic finite-geometry API. Independent review identified tolerance-band holes for arbitrary near-edge support inputs; the frozen manifest margins are far above those bands.",
    }
    _write_section_png(section, OUTPUT_PNG)
    # PNG rendering is also inside the one bounded diagnostic wall budget.
    elapsed = time.perf_counter() - started
    report["wall_clock_seconds"] = elapsed
    report["within_wall_clock_limit"] = elapsed <= WALL_LIMIT_S
    if not report["within_wall_clock_limit"]:
        report["status"] = "FAIL"
    OUTPUT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    return report


if __name__ == "__main__":
    result = run()
    print(json.dumps({
        "status": result["status"],
        "query_count": result["budgets"]["query_count_in_frozen_batch"],
        "mls_query_count": result["budgets"]["mls_query_count_after_fail_closed_geometry_preflight"],
        "provider_ok_rows": result["provider_contract_validation"]["ok_rows"],
        "provider_unknown_rows": result["provider_contract_validation"]["unknown_rows"],
        "affine_max_value_error": result["affine_base_reproduction"]["maximum_value_error"],
        "core_max_field_error_m": result["core_shifted_f1_reproduction"]["maximum_field_error_m"],
        "core_max_gradient_error": result["core_shifted_f1_reproduction"]["maximum_gradient_error_m_per_m"],
        "h_trace_jump_identity_error": result["finite_gamma_trace_checks"]["maximum_heaviside_jump_identity_error"],
        "ahead_h_max": result["no_infinite_plane_ghost_ahead_of_front"]["maximum_absolute_heaviside_value"],
        "whole_section_max_error_m": result["whole_box_z0_section"]["summary"]["maximum_absolute_error_m"],
        "wall_clock_seconds": result["wall_clock_seconds"],
        "json": str(OUTPUT_JSON),
        "png": str(OUTPUT_PNG),
    }, indent=2))
