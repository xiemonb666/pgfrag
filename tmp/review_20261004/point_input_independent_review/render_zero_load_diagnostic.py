#!/usr/bin/env python3
"""Render one static diagnostic view of the captured native zero-load run."""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
import pathlib
import platform
import sys
import tomllib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def number(row: dict[str, str], key: str) -> float:
    value = float(row[key])
    if not math.isfinite(value):
        raise ValueError(f"non-finite {key} in captured CSV")
    return value


def fmt(value: float) -> str:
    return f"{value:.4g}"


def value_range(values: list[float]) -> tuple[float, float]:
    if not values:
        raise ValueError("cannot summarize an empty sample")
    return min(values), max(values)


def vector_magnitude(row: dict[str, str], keys: tuple[str, str, str]) -> float:
    return math.sqrt(sum(number(row, key) ** 2 for key in keys))


def point_by_id(rows: list[dict[str, str]], step: int) -> dict[int, dict[str, str]]:
    selected = [row for row in rows if int(row["step"]) == step]
    result = {int(row["point_id"]): row for row in selected}
    if len(selected) != 125 or len(result) != 125 or set(result) != set(range(1, 126)):
        raise ValueError(f"step {step} must contain exactly one row for every source point ID 1..125")
    return result


def draw_cube(ax: object, extent_mm: float = 5.0) -> None:
    corners = list(itertools.product((0.0, extent_mm), repeat=3))
    for start in corners:
        for axis in range(3):
            end = list(start)
            end[axis] = extent_mm - end[axis]
            end_tuple = tuple(end)
            if start[axis] == 0.0:
                ax.plot(
                    [start[0], end_tuple[0]],
                    [start[1], end_tuple[1]],
                    [start[2], end_tuple[2]],
                    color="#64748b",
                    linewidth=0.75,
                    alpha=0.7,
                    zorder=0,
                )


def main() -> int:
    if len(sys.argv) != 4:
        raise SystemExit("usage: render_zero_load_diagnostic.py SNAPSHOT_DIR RECEIPT_DIR OUTPUT_PNG")
    snapshot_dir = pathlib.Path(sys.argv[1]).resolve()
    receipt_dir = pathlib.Path(sys.argv[2]).resolve()
    output_png = pathlib.Path(sys.argv[3]).resolve()

    input_paths = {
        "point_state_history.csv": snapshot_dir / "point_state_history.csv",
        "canonical_bond_state_history.csv": snapshot_dir / "canonical_bond_state_history.csv",
        "canonical_bond_transitions.csv": snapshot_dir / "canonical_bond_transitions.csv",
        "directed_bond_state_history.csv": snapshot_dir / "directed_bond_state_history.csv",
        "directed_bond_transitions.csv": snapshot_dir / "directed_bond_transitions.csv",
        "summary.toml": snapshot_dir / "summary.toml",
        "launch_receipt.json": receipt_dir / "zero_load_control_run01.launch_receipt.json",
        "freeze_receipt.json": receipt_dir / "zero_load_control_run01.freeze_receipt.json",
        "stdout.log": receipt_dir / "zero_load_control_run01.stdout.log",
        "stderr.log": receipt_dir / "zero_load_control_run01.stderr.log",
    }
    missing = [name for name, path in input_paths.items() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"missing frozen run inputs: {missing}")
    input_hashes = {name: sha256(path) for name, path in input_paths.items()}

    point_rows = read_csv(input_paths["point_state_history.csv"])
    canonical_rows = read_csv(input_paths["canonical_bond_state_history.csv"])
    canonical_transitions = read_csv(input_paths["canonical_bond_transitions.csv"])
    directed_rows = read_csv(input_paths["directed_bond_state_history.csv"])
    directed_transitions = read_csv(input_paths["directed_bond_transitions.csv"])
    summary = tomllib.loads(input_paths["summary.toml"].read_text(encoding="utf-8-sig"))
    launch = json.loads(input_paths["launch_receipt.json"].read_text(encoding="utf-8-sig"))
    freeze = json.loads(input_paths["freeze_receipt.json"].read_text(encoding="utf-8-sig"))
    run_stdout = input_paths["stdout.log"].read_text(encoding="utf-8-sig")
    run_stderr = input_paths["stderr.log"].read_text(encoding="utf-8-sig")

    steps = sorted({int(row["step"]) for row in point_rows})
    if steps != [0, 1]:
        raise ValueError(f"expected captured samples at steps 0 and 1, found {steps}")
    points = {step: point_by_id(point_rows, step) for step in steps}
    times = {step: sorted({number(row, "time_s") for row in point_rows if int(row["step"]) == step}) for step in steps}
    if any(len(values) != 1 for values in times.values()):
        raise ValueError("each step must have one captured time_s")
    sample_metrics = {}
    for step in steps:
        rows = [row for row in point_rows if int(row["step"]) == step]
        displacement = [vector_magnitude(row, ("ux_m", "uy_m", "uz_m")) for row in rows]
        velocity = [vector_magnitude(row, ("vx_mps", "vy_mps", "vz_mps")) for row in rows]
        damage = [number(row, "damage") for row in rows]
        sample_metrics[step] = {
            "time_s": times[step][0],
            "dt_s": sorted({number(row, "dt_s") for row in rows})[0],
            "displacement_m_range": list(value_range(displacement)),
            "velocity_m_s_range": list(value_range(velocity)),
            "damage_range": list(value_range(damage)),
            "max_displacement_norm_m": max(displacement),
            "max_velocity_norm_m_s": max(velocity),
        }

    if launch.get("status") != "COMPLETED" or int(launch.get("julia_process_exit_code", -1)) != 0:
        raise ValueError("the captured launcher receipt does not report a completed process exit 0")
    frozen_source_hashes = {record["path"]: str(record["sha256"]).lower() for record in freeze.get("source_file_hashes", [])}
    input_point_csv = next((path for path in frozen_source_hashes if path.endswith("coupon_uniform_voxel_5cube_125.points.csv")), None)
    if freeze.get("status", "").startswith("STATIC_PREFLIGHT_PASSED") is False or input_point_csv is None:
        raise ValueError("the frozen prelaunch receipt does not identify the prepared uniform point CSV")
    if frozen_source_hashes[input_point_csv] != str(summary.get("point_input_sha256", "")).lower():
        raise ValueError("the native summary input SHA does not match its frozen prelaunch CSV pin")
    if summary.get("status") != "COMPLETED" or summary.get("load_kind") != "zero_load":
        raise ValueError("the captured run summary is not a completed zero-load case")
    if int(summary.get("source_point_count", 0)) != 125 or int(summary.get("n_steps", 0)) != 1 or int(summary.get("sample_count_including_step_zero", 0)) != 2:
        raise ValueError("the summary does not confirm 125 points, one native step, and two samples")
    if "PD_POINT_EVENT_PROBE_COMPLETED points=125 steps=1" not in run_stdout:
        raise ValueError("the captured native stdout does not confirm 125 points and one step")
    if "precompilation errored" not in run_stderr or "ArgumentError:" not in run_stderr:
        raise ValueError("the captured native stderr no longer contains the reported precompile error")
    if len(canonical_transitions) != 0 or len(directed_transitions) != 0:
        raise ValueError("the captured zero-load run contains transition rows; inspect before plotting")

    canonical_step_one = [row for row in canonical_rows if int(row["step"]) == 1]
    canonical_step_zero = [row for row in canonical_rows if int(row["step"]) == 0]
    if not canonical_step_one or not canonical_step_zero:
        raise ValueError("canonical bond states are missing step 0 or step 1")
    edges = sorted(
        {
            (int(row["point_i_id"]), int(row["point_j_id"]))
            for row in canonical_step_one
            if row["active"].strip().lower() == "true"
        }
    )
    if not edges:
        raise ValueError("captured step-one canonical graph has no active bonds")
    target_edge_count = 100
    edge_stride = max(1, math.ceil(len(edges) / target_edge_count))
    sampled_edges = edges[::edge_stride]

    reference = points[0]
    current = points[1]
    rest_mm = {
        point_id: tuple(number(row, key) * 1000.0 for key in ("x_ref_m", "y_ref_m", "z_ref_m"))
        for point_id, row in reference.items()
    }
    current_mm = {
        point_id: tuple(number(row, key) * 1000.0 for key in ("x_current_m", "y_current_m", "z_current_m"))
        for point_id, row in current.items()
    }
    max_rest_current_delta_m = max(
        math.dist(
            tuple(number(current[point_id], key) for key in ("x_current_m", "y_current_m", "z_current_m")),
            tuple(number(reference[point_id], key) for key in ("x_ref_m", "y_ref_m", "z_ref_m")),
        )
        for point_id in range(1, 126)
    )

    fig = plt.figure(figsize=(10.0, 6.4), facecolor="white")
    ax = fig.add_axes([0.025, 0.12, 0.59, 0.78], projection="3d")
    text_ax = fig.add_axes([0.64, 0.10, 0.34, 0.80])
    text_ax.axis("off")
    draw_cube(ax)

    for edge_index, (point_i, point_j) in enumerate(sampled_edges):
        first = rest_mm[point_i]
        second = rest_mm[point_j]
        ax.plot(
            [first[0], second[0]],
            [first[1], second[1]],
            [first[2], second[2]],
            color="#94a3b8",
            linewidth=0.42,
            alpha=0.22,
            zorder=1,
            label="Sparse captured-bond sample" if edge_index == 0 else None,
        )

    ax.scatter(
        [rest_mm[index][0] for index in range(1, 126)],
        [rest_mm[index][1] for index in range(1, 126)],
        [rest_mm[index][2] for index in range(1, 126)],
        s=22,
        marker="o",
        facecolors="none",
        edgecolors="#2563eb",
        linewidths=0.8,
        alpha=0.9,
        label="Rest points, step 0",
        depthshade=False,
        zorder=3,
    )
    ax.scatter(
        [current_mm[index][0] for index in range(1, 126)],
        [current_mm[index][1] for index in range(1, 126)],
        [current_mm[index][2] for index in range(1, 126)],
        s=20,
        marker="x",
        color="#dc2626",
        linewidths=0.9,
        alpha=0.9,
        label="Current points, step 1",
        depthshade=False,
        zorder=4,
    )

    ax.set_xlim(-0.2, 5.2)
    ax.set_ylim(-0.2, 5.2)
    ax.set_zlim(-0.2, 5.2)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xticks((0, 2.5, 5))
    ax.set_yticks((0, 2.5, 5))
    ax.set_zticks((0, 2.5, 5))
    ax.set_xlabel("X (mm)", labelpad=2)
    ax.set_ylabel("Y (mm)", labelpad=2)
    ax.set_zlabel("Z (mm)", labelpad=2)
    ax.tick_params(labelsize=7, pad=0)
    ax.view_init(elev=23, azim=38)
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper left", bbox_to_anchor=(-0.04, 1.02), fontsize=7, frameon=False)

    fig.suptitle("PG-Frag native zero-load geometry + callback diagnostic", fontsize=14, fontweight="bold", y=0.97)
    text_ax.text(0.0, 0.98, "NATIVE RUN", fontsize=9, fontweight="bold", color="#334155", va="top")
    text_ax.text(0.0, 0.91, "COMPLETED | Julia process exit 0", fontsize=9, color="#166534", va="top")
    text_ax.text(0.0, 0.84, "125 native points | one zero-load step", fontsize=8.5, color="#0f172a", va="top")
    text_ax.text(0.0, 0.77, "Synthetic 5 mm cube | horizon 3 mm", fontsize=8, color="#475569", va="top")
    text_ax.text(0.0, 0.73, f"input_id: {summary.get('input_id', 'UNKNOWN')}", fontsize=7.1, color="#64748b", va="top")
    text_ax.text(0.0, 0.69, "CAPTURED SAMPLES", fontsize=9, fontweight="bold", color="#334155", va="top")

    y_position = 0.62
    for step in steps:
        metrics = sample_metrics[step]
        disp = metrics["displacement_m_range"]
        vel = metrics["velocity_m_s_range"]
        damage = metrics["damage_range"]
        sample_text = (
            f"Step {step} | t = {fmt(metrics['time_s'])} s\n"
            f"||U||: {fmt(disp[0])} to {fmt(disp[1])} m\n"
            f"||v||: {fmt(vel[0])} to {fmt(vel[1])} m/s\n"
            f"damage: {fmt(damage[0])} to {fmt(damage[1])}"
        )
        text_ax.text(0.0, y_position, sample_text, fontsize=8, color="#0f172a", va="top", linespacing=1.28, family="monospace")
        y_position -= 0.22

    text_ax.text(0.0, 0.24, "OBSERVED", fontsize=9, fontweight="bold", color="#334155", va="top")
    text_ax.text(
        0.0,
        0.18,
        f"dt = {fmt(sample_metrics[1]['dt_s'])} s | directed/canonical transitions = 0/0\n"
        f"Bond sample: {len(sampled_edges)}/{len(edges)} captured intact pairs.\n"
        "Grey bonds show connectivity only, not fracture.\n"
        "Rest/current markers overlap at U = 0.\n"
        "Synthetic fixture; not a real artifact/specimen.\n"
        "Stderr records a Peridynamics precompile error.",
        fontsize=7.1,
        color="#475569",
        va="top",
        linespacing=1.35,
    )

    fig.savefig(output_png, dpi=120, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    output = {
        "status": "RENDERED",
        "argv": sys.argv,
        "python_version": platform.python_version(),
        "matplotlib_version": matplotlib.__version__,
        "input_sha256": input_hashes,
        "launch_receipt_status": launch.get("status"),
        "julia_process_exit_code": launch.get("julia_process_exit_code"),
        "summary_status": summary.get("status"),
        "summary_load_kind": summary.get("load_kind"),
        "freeze_receipt_status": freeze.get("status"),
        "frozen_point_csv_path": input_point_csv,
        "frozen_point_csv_sha256": frozen_source_hashes[input_point_csv],
        "summary_point_input_sha256": summary.get("point_input_sha256"),
        "native_stdout_completion_line": run_stdout.strip(),
        "native_stderr_contains_precompile_error": "precompilation errored" in run_stderr,
        "native_stderr_excerpt": run_stderr.strip(),
        "steps": steps,
        "sample_metrics": sample_metrics,
        "point_count_per_step": {str(step): len(points[step]) for step in steps},
        "max_rest_to_current_point_delta_m": max_rest_current_delta_m,
        "canonical_pairs_at_step_one": len(edges),
        "canonical_active_pairs_rendered": len(sampled_edges),
        "canonical_edge_stride": edge_stride,
        "canonical_transitions": len(canonical_transitions),
        "directed_transitions": len(directed_transitions),
        "directed_state_rows": len(directed_rows),
        "summary_canonical_state_rows": summary.get("canonical_bond_state_rows"),
        "summary_directed_state_rows": summary.get("directed_bond_state_rows"),
        "output_png": str(output_png),
        "output_png_bytes": output_png.stat().st_size,
        "output_png_sha256": sha256(output_png),
        "max_png_bytes_requested": 200000,
        "png_within_size_limit": output_png.stat().st_size <= 200000,
        "max_rest_to_current_point_delta_m_note": "Measured from captured CSV values; no synthetic displacement was added.",
        "bond_wireframe_note": "Sparse deterministic sample of captured active native canonical bonds; visual connectivity only and does not depict fracture.",
    }
    print(json.dumps(output, indent=2, sort_keys=True))
    return 0 if output["png_within_size_limit"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
