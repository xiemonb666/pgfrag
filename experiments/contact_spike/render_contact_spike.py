"""Render the CSV state emitted by run_contact_spike.jl using the installed Pillow runtime."""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[3]
CONTACT_SOURCE = ROOT / "pgfrag" / "vendor" / "Peridynamics.jl" / "src" / "physics" / "short_range_force_contact.jl"
VERLET_SOURCE = ROOT / "pgfrag" / "vendor" / "Peridynamics.jl" / "src" / "time_solvers" / "velocity_verlet.jl"
MULTIBODY_SOURCE = ROOT / "pgfrag" / "vendor" / "Peridynamics.jl" / "src" / "discretization" / "multibody_setup.jl"
RUNNER_SOURCE = ROOT / "pgfrag" / "experiments" / "contact_spike" / "run_contact_spike.jl"

CASES = (
    ("no_contact_control", "No-contact control"),
    ("contact_baseline", "Native contact, 0.5 x BB bond bound"),
    ("contact_half_dt", "Native contact, dt / 2"),
)
COLORS = {"a": (32, 104, 190), "b": (216, 112, 40)}
FORCE_COLORS = {"a": (32, 158, 100), "b": (156, 72, 160)}


def file_hash(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def number(row, key):
    value = row.get(key, "")
    return float(value) if value not in (None, "") else math.nan


def grouped_trajectories(rows):
    result = {}
    for row in rows:
        key = (row["body"], int(row["body_local_point_id"]))
        result.setdefault(key, []).append((int(row["step"]), number(row, "x_m") * 1000.0,
                                           number(row, "z_m") * 1000.0))
    for points in result.values():
        points.sort(key=lambda item: item[0])
    return result


def final_point_rows(rows):
    final_step = max((int(row["step"]) for row in rows), default=-1)
    return [row for row in rows if int(row["step"]) == final_step], final_step


def bounds_for(rows):
    xs = [number(row, "x_m") * 1000.0 for row in rows]
    zs = [number(row, "z_m") * 1000.0 for row in rows]
    xs = [x for x in xs if math.isfinite(x)]
    zs = [z for z in zs if math.isfinite(z)]
    if not xs or not zs:
        return (-2.0, 2.0, -4.0, 4.0)
    xmin, xmax = min(xs), max(xs)
    zmin, zmax = min(zs), max(zs)
    xpad = max((xmax - xmin) * 0.12, 0.35)
    zpad = max((zmax - zmin) * 0.12, 0.35)
    return xmin - xpad, xmax + xpad, zmin - zpad, zmax + zpad


def draw_axes(draw, box, limits, x_label, z_label, font, color=(65, 65, 65)):
    left, top, right, bottom = box
    xmin, xmax, zmin, zmax = limits
    draw.rectangle(box, outline=color, width=2)
    for i in range(5):
        frac = i / 4.0
        x = left + frac * (right - left)
        z = top + frac * (bottom - top)
        draw.line((x, top, x, bottom), fill=(224, 228, 232), width=1)
        draw.line((left, z, right, z), fill=(224, 228, 232), width=1)
        xv = xmin + frac * (xmax - xmin)
        zv = zmax - frac * (zmax - zmin)
        draw.text((x - 24, bottom + 7), f"{xv:.2g}", fill=(45, 45, 45), font=font)
        draw.text((left - 60, z - 7), f"{zv:.2g}", fill=(45, 45, 45), font=font)
    draw.text(((left + right) / 2 - 22, bottom + 31), x_label, fill=(25, 25, 25), font=font)
    draw.text((left - 62, top - 23), z_label, fill=(25, 25, 25), font=font)


def map_xy(point, box, limits):
    left, top, right, bottom = box
    xmin, xmax, zmin, zmax = limits
    x = left + (point[0] - xmin) / max(xmax - xmin, 1e-12) * (right - left)
    y = bottom - (point[1] - zmin) / max(zmax - zmin, 1e-12) * (bottom - top)
    return x, y


def draw_path_panel(draw, box, rows, font):
    limits = bounds_for(rows)
    draw_axes(draw, box, limits, "x (mm)", "z (mm)", font)
    trajectories = grouped_trajectories(rows)
    for (body, _point_id), samples in trajectories.items():
        points = [map_xy((x, z), box, limits) for _step, x, z in samples]
        if len(points) > 1:
            draw.line(points, fill=COLORS[body], width=2)
        if points:
            x0, y0 = points[0]
            xf, yf = points[-1]
            draw.ellipse((x0 - 4, y0 - 4, x0 + 4, y0 + 4), outline=COLORS[body], width=2)
            draw.rectangle((xf - 4, yf - 4, xf + 4, yf + 4), fill=COLORS[body])
    draw.text((box[0] + 5, box[1] + 5), "Trajectories: circles = initial, squares = final", fill=(30, 30, 30), font=font)
    return limits


def draw_force_panel(draw, box, rows, font):
    limits = bounds_for(rows)
    draw_axes(draw, box, limits, "x (mm)", "z (mm)", font)
    finals, final_step = final_point_rows(rows)
    max_force = max((math.hypot(number(row, "native_contact_force_x_N"),
                                number(row, "native_contact_force_z_N")) for row in finals), default=0.0)
    mm_per_newton = 0.0 if not math.isfinite(max_force) or max_force <= 0 else 0.9 / max_force
    for row in finals:
        body = row["body"]
        xy = map_xy((number(row, "x_m") * 1000.0, number(row, "z_m") * 1000.0), box, limits)
        draw.ellipse((xy[0] - 4, xy[1] - 4, xy[0] + 4, xy[1] + 4), fill=COLORS[body])
        fx = number(row, "native_contact_force_x_N")
        fz = number(row, "native_contact_force_z_N")
        if mm_per_newton and math.isfinite(fx) and math.isfinite(fz):
            # Plot arrows in physical point coordinates; visual scaling is disclosed below.
            scale = mm_per_newton
            end = map_xy((number(row, "x_m") * 1000.0 + fx * scale,
                          number(row, "z_m") * 1000.0 + fz * scale), box, limits)
            draw.line((xy[0], xy[1], end[0], end[1]), fill=FORCE_COLORS[body], width=3)
            draw.ellipse((end[0] - 3, end[1] - 3, end[0] + 3, end[1] + 3), fill=FORCE_COLORS[body])
    if max_force > 0:
        scale_text = f"arrows scaled: 0.9 mm = {max_force:.4g} N (largest final point force)"
    else:
        scale_text = "final native contact force = 0 N at all points"
    draw.text((box[0] + 5, box[1] + 5), f"Final contact force; sample step {final_step}", fill=(30, 30, 30), font=font)
    draw.text((box[0] + 5, box[3] - 17), scale_text, fill=(35, 35, 35), font=font)


def final_state_by_key(rows):
    finals, step = final_point_rows(rows)
    return {(row["body"], int(row["body_local_point_id"])): row for row in finals}, step


def compare_same_time(outdir):
    baseline_path = outdir / "contact_baseline" / "point_steps.csv"
    half_path = outdir / "contact_half_dt" / "point_steps.csv"
    baseline_summary_path = outdir / "contact_baseline" / "summary.json"
    half_summary_path = outdir / "contact_half_dt" / "summary.json"
    if not all(path.exists() for path in (baseline_path, half_path, baseline_summary_path, half_summary_path)):
        return {"status": "NOT_AVAILABLE", "reason": "both contact cases must have output files"}
    base_summary = read_json(baseline_summary_path)
    half_summary = read_json(half_summary_path)
    if base_summary.get("status") != "COMPLETED" or half_summary.get("status") != "COMPLETED":
        return {"status": "NOT_AVAILABLE", "reason": "one or both contact cases did not complete"}
    base_rows = read_csv(baseline_path)
    half_rows = read_csv(half_path)
    base, base_step = final_state_by_key(base_rows)
    half, half_step = final_state_by_key(half_rows)
    keys = sorted(set(base) & set(half))
    if len(keys) != 16:
        return {"status": "NOT_AVAILABLE", "reason": f"expected 16 matching points, got {len(keys)}"}
    position_deltas = []
    velocity_deltas = []
    for key in keys:
        position_deltas.append(math.sqrt(sum((number(base[key], axis) - number(half[key], axis)) ** 2
                                             for axis in ("x_m", "y_m", "z_m"))))
        velocity_deltas.append(math.sqrt(sum((number(base[key], axis) - number(half[key], axis)) ** 2
                                             for axis in ("vx_mps", "vy_mps", "vz_mps"))))
    return {
        "status": "COMPLETED",
        "comparison_time_s": base_summary.get("end_time_s"),
        "baseline_dt_s": base_summary.get("dt_s"),
        "half_dt_s": half_summary.get("dt_s"),
        "baseline_final_step": base_step,
        "half_dt_final_step": half_step,
        "matching_points": len(keys),
        "position_difference_rms_m": math.sqrt(sum(value * value for value in position_deltas) / len(keys)),
        "position_difference_max_m": max(position_deltas),
        "velocity_difference_rms_m_s": math.sqrt(sum(value * value for value in velocity_deltas) / len(keys)),
        "velocity_difference_max_m_s": max(velocity_deltas),
        "interpretation": "same-end-time pointwise difference between native contact runs; not a convergence proof",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    outdir = Path(args.outdir).resolve()
    report_path = outdir / "report.json"
    report = read_json(report_path)
    plot_path = outdir / "contact_spike_actual_xz.png"
    canvas = Image.new("RGB", (2200, 1550), (255, 255, 255))
    draw = ImageDraw.Draw(canvas)
    font = ImageFont.load_default()
    title_font = ImageFont.load_default(size=22)
    label_font = ImageFont.load_default(size=14)
    draw.text((64, 28), "RelicWeave PG-Frag: native point-contact spike", fill=(15, 30, 50), font=title_font)
    draw.text((64, 60), "Discrete point-to-point state; no CAD surface reconstruction or fracture display", fill=(55, 55, 55), font=label_font)

    statuses = {case.get("case_id"): case.get("status") for case in report.get("cases", [])}
    top = 100
    row_height = 440
    left_box_x = 100
    right_box_x = 1170
    panel_w = 900
    panel_h = 300
    has_saved_state = False
    for index, (case_id, label) in enumerate(CASES):
        y = top + index * row_height
        draw.text((left_box_x, y), f"{label} | status={statuses.get(case_id, 'NO OUTPUT')}",
                  fill=(10, 25, 40), font=label_font)
        point_path = outdir / case_id / "point_steps.csv"
        if not point_path.exists():
            draw.text((left_box_x, y + 36), "No point-state CSV was produced for this case.", fill=(180, 20, 20), font=font)
            continue
        rows = read_csv(point_path)
        if not rows:
            draw.text((left_box_x, y + 36), "No point-state rows were recorded.", fill=(180, 20, 20), font=font)
            continue
        has_saved_state = True
        left_box = (left_box_x, y + 40, left_box_x + panel_w, y + 40 + panel_h)
        right_box = (right_box_x, y + 40, right_box_x + panel_w, y + 40 + panel_h)
        draw_path_panel(draw, left_box, rows, font)
        draw_force_panel(draw, right_box, rows, font)

    if not has_saved_state:
        raise RuntimeError("no actual point state was recorded; refusing to emit an empty state/force plot")

    hashes = {
        "contact_src_sha256": file_hash(CONTACT_SOURCE),
        "verlet_src_sha256": file_hash(VERLET_SOURCE),
        "multibody_src_sha256": file_hash(MULTIBODY_SOURCE),
        "julia_runner_sha256": file_hash(RUNNER_SOURCE),
        "plotter_sha256": file_hash(Path(__file__).resolve()),
    }
    footer_y = 1435
    draw.text((64, footer_y), "Fields: x_m,z_m,native_contact_force_x_N,native_contact_force_z_N; positions shown in mm; force arrows in N.",
              fill=(20, 20, 20), font=font)
    draw.text((64, footer_y + 18), f"SHA256 contact={hashes['contact_src_sha256']}", fill=(35, 35, 35), font=font)
    draw.text((64, footer_y + 36), f"SHA256 runner={hashes['julia_runner_sha256']} plotter={hashes['plotter_sha256']}",
              fill=(35, 35, 35), font=font)

    canvas.save(plot_path, format="PNG", optimize=True)
    report["same_end_time_contact_dt_comparison"] = compare_same_time(outdir)
    report["plot"] = {
        "path": str(plot_path),
        "sha256": file_hash(plot_path),
        "plotted_fields": ["x_m", "z_m", "native_contact_force_x_N", "native_contact_force_z_N"],
        "coordinate_display_unit": "mm",
        "force_display_unit": "N; arrow scale written into image",
        "source_hashes_sha256": hashes,
    }
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"PLOT_COMPLETED path={plot_path} sha256={report['plot']['sha256']}")
    if report.get("status") != "COMPLETED":
        print("PLOT_NOTE some native cases failed; image includes only saved point rows and status labels")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
