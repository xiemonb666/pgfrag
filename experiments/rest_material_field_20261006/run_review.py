"""Run one bounded field-contract check and save review artifacts."""

from __future__ import annotations

from datetime import datetime, timezone
import csv
import hashlib
import json
from pathlib import Path
import platform
import shutil
import sys
import time
from typing import Any

import numpy as np
import PIL
from PIL import Image, ImageDraw, ImageFont

from focused_checks import CountedField, QueryLedger, run_focused_checks
from rest_material_field import MaterialFieldParameters, RestFrameMaterialField


HERE = Path(__file__).resolve().parent
PGFRAG_ROOT = HERE.parents[1]
OUTPUT = PGFRAG_ROOT / "tmp" / "review_20261006" / "rest_material_field_20261006"
SAMPLE_SIDE = 14
DOMAIN_HALF_WIDTH_M = 0.018
MAX_QUERY_ATTEMPTS = 256
MAX_RUNTIME_SECONDS = 30.0
FAILED_IMPLEMENTATION_ATTEMPTS = 26


def _write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, sort_keys=True, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _font(size: int) -> ImageFont.ImageFont:
    for candidate in (r"C:\Windows\Fonts\arial.ttf", r"C:\Windows\Fonts\segoeui.ttf"):
        if Path(candidate).is_file():
            return ImageFont.truetype(candidate, size=size)
    try:
        return ImageFont.load_default(size=size)
    except TypeError:
        return ImageFont.load_default()


def _interpolated_diverging_colormap(values: np.ndarray, limit: float) -> Image.Image:
    # Compact blue-white-red map.  Colour interpolation affects display only;
    # all underlying values are exact stored field evaluations.
    navy = np.asarray([39.0, 73.0, 158.0])
    white = np.asarray([247.0, 247.0, 244.0])
    red = np.asarray([178.0, 24.0, 43.0])
    t = np.clip((values / limit + 1.0) * 0.5, 0.0, 1.0)
    rgb = np.empty((*values.shape, 3), dtype=np.uint8)
    lower = t <= 0.5
    lower_t = (t[lower] * 2.0)[:, None]
    upper_t = ((t[~lower] - 0.5) * 2.0)[:, None]
    rgb[lower] = np.rint(navy + lower_t * (white - navy)).astype(np.uint8)
    rgb[~lower] = np.rint(white + upper_t * (red - white)).astype(np.uint8)
    return Image.fromarray(rgb, mode="RGB")


def _draw_field_plot(path: Path, values: np.ndarray, x_values_m: np.ndarray, y_values_m: np.ndarray, amplitude: float) -> dict[str, float]:
    sample_min = float(np.min(values))
    sample_max = float(np.max(values))
    plot_limit = max(abs(sample_min), abs(sample_max), 1e-12)

    canvas = Image.new("RGB", (1120, 820), (255, 255, 255))
    draw = ImageDraw.Draw(canvas)
    title_font, body_font, small_font = _font(28), _font(18), _font(14)
    draw.text((96, 24), "Seeded rest-frame material heterogeneity", fill=(22, 31, 48), font=title_font)
    draw.text(
        (98, 62),
        f"z = 0 mm | {values.size} evaluated points | display stretched to observed ±{plot_limit:.4f}",
        fill=(67, 78, 96),
        font=body_font,
    )

    left, top, side = 130, 125, 600
    image = _interpolated_diverging_colormap(np.flipud(values), plot_limit)
    image = image.resize((side, side), resample=Image.Resampling.BICUBIC)
    canvas.paste(image, (left, top))
    draw.rectangle((left, top, left + side, top + side), outline=(35, 43, 57), width=2)

    half_mm = float(max(abs(x_values_m[0]), abs(x_values_m[-1])) * 1000.0)
    for fraction, label in ((0.0, f"-{half_mm:g}"), (0.5, "0"), (1.0, f"+{half_mm:g}")):
        x = left + int(round(fraction * side))
        draw.line((x, top + side, x, top + side + 8), fill=(35, 43, 57), width=2)
        draw.text((x - 18, top + side + 12), label, fill=(35, 43, 57), font=small_font)
        y = top + int(round((1.0 - fraction) * side))
        draw.line((left - 8, y, left, y), fill=(35, 43, 57), width=2)
        draw.text((left - 52, y - 8), label, fill=(35, 43, 57), font=small_font)

    draw.text((left + 245, top + side + 40), "x (mm)", fill=(35, 43, 57), font=body_font)
    y_label = Image.new("RGBA", (120, 28), (255, 255, 255, 0))
    ImageDraw.Draw(y_label).text((0, 0), "y (mm)", fill=(35, 43, 57), font=body_font)
    canvas.paste(y_label.rotate(90, expand=True), (18, top + 240), y_label.rotate(90, expand=True))

    bar_left = 780
    bar_top = top
    bar_width = 36
    bar_height = side
    colorbar_values = np.linspace(1.0, -1.0, bar_height, dtype=np.float64).reshape(-1, 1)
    colorbar = _interpolated_diverging_colormap(colorbar_values, 1.0).resize((bar_width, bar_height))
    canvas.paste(colorbar, (bar_left, bar_top))
    draw.rectangle((bar_left, bar_top, bar_left + bar_width, bar_top + bar_height), outline=(35, 43, 57), width=1)
    draw.text((bar_left + 49, bar_top - 4), f"+{plot_limit:.3f}", fill=(35, 43, 57), font=small_font)
    draw.text((bar_left + 49, bar_top + bar_height // 2 - 8), "0", fill=(35, 43, 57), font=small_font)
    draw.text((bar_left + 49, bar_top + bar_height - 18), f"-{plot_limit:.3f}", fill=(35, 43, 57), font=small_font)
    draw.text((bar_left, bar_top + bar_height + 20), "heterogeneity", fill=(35, 43, 57), font=small_font)

    notes = [
        "Status: UNCALIBRATED",
        "Frame: material rest; coords m",
        "Ensemble e-fold target: 12 mm",
        f"Max relative deviation: ±{amplitude:.2f}",
        "Colour range: saved sample min/max",
        "No Weibull or calibration claim",
        "Image interpolation adds no queries",
    ]
    note_y = 150
    for line in notes:
        draw.text((855, note_y), line, fill=(44, 56, 73), font=small_font)
        note_y += 34
    canvas.save(path, format="PNG", optimize=True)
    return {"display_symmetric_limit": plot_limit, "sample_min": sample_min, "sample_max": sample_max}


def main() -> int:
    start = time.perf_counter()
    OUTPUT.mkdir(parents=True, exist_ok=True)

    # One corrected focused validation pass; these 31 checks include the
    # rejected non-finite point attempt.  The previous implementation attempt
    # used 26 points before its first failed assertion, which is included in
    # cumulative query accounting below.
    checks = run_focused_checks()
    _write_json(OUTPUT / "focused_checks.json", checks)

    field = RestFrameMaterialField(
        MaterialFieldParameters(seed=20261006, correlation_length_m=0.012, amplitude=0.30)
    )
    sample_ledger = QueryLedger()
    counted = CountedField(field, sample_ledger)
    x_values = np.linspace(-DOMAIN_HALF_WIDTH_M, DOMAIN_HALF_WIDTH_M, SAMPLE_SIDE)
    y_values = np.linspace(-DOMAIN_HALF_WIDTH_M, DOMAIN_HALF_WIDTH_M, SAMPLE_SIDE)
    query_points = np.asarray([(x, y, 0.0) for y in y_values for x in x_values], dtype=np.float64)
    sampled_values = np.asarray(counted.evaluate_many(query_points), dtype=np.float64).reshape(SAMPLE_SIDE, SAMPLE_SIDE)
    total_queries = int(checks["point_evaluation_attempts"] + sample_ledger.count)
    if sample_ledger.count != SAMPLE_SIDE * SAMPLE_SIDE:
        raise AssertionError("sample grid evaluation count does not match its 15-by-15 coordinates")
    cumulative_queries = FAILED_IMPLEMENTATION_ATTEMPTS + total_queries
    if cumulative_queries > MAX_QUERY_ATTEMPTS:
        raise AssertionError(f"cumulative field query budget exceeded: {cumulative_queries} > {MAX_QUERY_ATTEMPTS}")

    manifest_path = OUTPUT / "field_manifest.json"
    manifest_path.write_text(field.manifest_json(), encoding="utf-8")
    samples_path = OUTPUT / "field_samples.csv"
    with samples_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["x_m", "y_m", "z_m", "heterogeneity_dimensionless"])
        for y_index, y in enumerate(y_values):
            for x_index, x in enumerate(x_values):
                writer.writerow(
                    [
                        format(float(x), ".17g"),
                        format(float(y), ".17g"),
                        "0",
                        format(float(sampled_values[y_index, x_index]), ".17g"),
                    ]
                )

    plot_path = OUTPUT / "field_samples.png"
    plot_summary = _draw_field_plot(plot_path, sampled_values, x_values, y_values, field.parameters.amplitude)
    method_copy = OUTPUT / "METHOD_AND_LIMITS.md"
    shutil.copyfile(HERE / "METHOD_AND_LIMITS.md", method_copy)

    check_seconds = time.perf_counter() - start
    if check_seconds > MAX_RUNTIME_SECONDS:
        raise AssertionError(f"focused checks and artifact preparation exceeded {MAX_RUNTIME_SECONDS:g} seconds")

    files = [manifest_path, samples_path, plot_path, OUTPUT / "focused_checks.json", method_copy]
    artifact_hashes = {path.name: _sha256(path) for path in files}
    _write_json(OUTPUT / "artifacts_sha256.json", artifact_hashes)

    receipt = {
        "schema": "pgfrag.rest_material_field.review_receipt.v1",
        "status": "PASS",
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "field_id_sha256": field.field_id,
        "field_status": "UNCALIBRATED",
        "query_accounting": {
            "limit": MAX_QUERY_ATTEMPTS,
            "prior_failed_implementation_attempts": FAILED_IMPLEMENTATION_ATTEMPTS,
            "focused_check_attempts_including_rejected_query": int(checks["point_evaluation_attempts"]),
            "saved_sample_points": int(sample_ledger.count),
            "current_artifact_run_attempts": total_queries,
            "cumulative_attempts_including_prior_failed_run": cumulative_queries,
            "png_render_and_hash_calls": 0,
        },
        "execution": {
            "runtime_seconds": round(check_seconds, 6),
            "runtime_limit_seconds": MAX_RUNTIME_SECONDS,
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pillow": PIL.__version__,
            "solver_invoked": False,
            "production_or_solver_paths_modified": False,
            "git_commands_used": False,
        },
        "checks": checks,
        "sample_grid": {
            "shape": [SAMPLE_SIDE, SAMPLE_SIDE],
            "plane": "z=0 m in the material rest frame",
            "domain_m": [[float(x_values[0]), float(x_values[-1])], [float(y_values[0]), float(y_values[-1])]],
            "sample_statistics": {
                "minimum": float(np.min(sampled_values)),
                "maximum": float(np.max(sampled_values)),
                "mean": float(np.mean(sampled_values)),
                "rms_about_zero": float(np.sqrt(np.mean(sampled_values**2))),
            },
            "plot_display": plot_summary,
        },
        "pins": [
            {"path": path.name, "sha256": artifact_hashes[path.name]}
            for path in files
        ],
        "scope": [
            "new isolated experiment module and focused checks only",
            "no solver or production-material integration",
            "no measured-material calibration or sourced nominal E/Gc/strength values",
            "RelicWeave-compatible lineage metadata is exposed as metadata-only; no upstream material artifact is claimed",
        ],
    }
    _write_json(OUTPUT / "receipt.json", receipt)
    print(json.dumps({"status": "PASS", "output": str(OUTPUT), "queries": total_queries, "runtime_seconds": check_seconds}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
