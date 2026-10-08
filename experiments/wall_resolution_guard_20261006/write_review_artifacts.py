from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Rectangle


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
REVIEW_DIR = REPO_ROOT / "pgfrag" / "tmp" / "review_20261006" / "wall_resolution_guard_20261006"
sys.path.insert(0, str(HERE))

from thin_wall_qualification import (  # noqa: E402
    historical_h8_assessment,
    make_analytic_hollow_cylinder_fixture,
    qualify_wall_input,
)


def _json_write(path: Path, value: object) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def _fixture_cross_section(ax: plt.Axes, payload: dict, title: str) -> tuple[list[float], list[float]]:
    cells = payload["cells"]
    point_by_cell = {row["cell_id"]: row for row in payload["centroid_pointweights"]}
    radial_cells = []
    for cell in cells:
        cg = cell["geometry"]
        tmid = 0.5 * (cg["theta_interval"][0]["value"] + cg["theta_interval"][1]["value"])
        if abs(tmid - 0.0) <= 0.05 + 1e-12:
            radial_cells.append((cg["z_min"]["value"], cg["inner_radius"]["value"], cell, point_by_cell[cell["cell_id"]], tmid))
    # Keep one actual angular cut; every row is an exact radial/z cell section.
    cut_thetas = {round(row[4], 12) for row in radial_cells}
    chosen_theta = min(cut_thetas, key=abs) if cut_thetas else 0.0
    radial_cells = [row for row in radial_cells if abs(row[4] - chosen_theta) < 1e-10]
    values_mm3 = [row[3]["point_weight"]["value"] for row in radial_cells]
    norm = Normalize(vmin=min(values_mm3), vmax=max(values_mm3) if max(values_mm3) > min(values_mm3) else min(values_mm3) + 1e-12)
    cmap = plt.get_cmap("viridis")
    centroid_x: list[float] = []
    centroid_y: list[float] = []
    for z0, r0, cell, point, theta in radial_cells:
        cg = cell["geometry"]
        r1 = cg["outer_radius"]["value"]
        z1 = cg["z_max"]["value"]
        weight = point["point_weight"]["value"]
        ax.add_patch(Rectangle((r0 - 8.0, z0), r1 - r0, z1 - z0, facecolor=cmap(norm(weight)), edgecolor="white", linewidth=1.0))
        px, py, pz = point["position"]["value"]
        depth = (px * math.cos(theta) + py * math.sin(theta)) - 8.0
        centroid_x.append(depth)
        centroid_y.append(pz)
        ax.text((r0 + r1) * 0.5 - 8.0, (z0 + z1) * 0.5, f"{weight:.3f}", ha="center", va="center", fontsize=7, color="white")
    ax.scatter(centroid_x, centroid_y, marker="o", s=22, facecolor="white", edgecolor="black", linewidth=0.65, zorder=4, label="actual centroid position")
    report = qualify_wall_input(payload)
    relevant = [row for row in report["resolution"]["wall_samples"] if row["sample_id"].startswith("ws-a01-")]
    actual_counts = sorted({row["crossed_cell_intervals"] for row in relevant})
    centroid_counts = sorted({row["actual_thickness_centroid_layer_count"] for row in relevant})
    ax.set_title(
        f"{title}\nselected ray: {actual_counts[0] if actual_counts else 0} actual spans / "
        f"{centroid_counts[0] if centroid_counts else 0} centroid layers",
        fontsize=11,
    )
    ax.set_xlim(0.0, 2.0)
    ax.set_ylim(4.0, 6.0)
    ax.set_xlabel("depth from inner wall (mm)")
    ax.set_ylabel("z (mm)")
    ax.set_xticks([0.0, 0.5, 1.0, 1.5, 2.0])
    ax.set_yticks([4.0, 4.5, 5.0, 5.5, 6.0])
    ax.grid(color="#596579", alpha=0.18, linewidth=0.5)
    ax.set_aspect("equal", adjustable="box")
    ax.legend(loc="upper right", fontsize=7)
    return values_mm3, [v for v in centroid_x]


def _write_section_figure(fine: dict, coarse: dict, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 6.8), constrained_layout=False)
    fig.subplots_adjust(left=0.075, right=0.90, bottom=0.17, top=0.79, wspace=0.20)
    fine_values, _ = _fixture_cross_section(axes[0], fine, "Resolved fixture: h=0.5 mm, t=2 mm")
    coarse_values, _ = _fixture_cross_section(axes[1], coarse, "Underresolved fixture: h=2 mm, t=2 mm")
    all_values = fine_values + coarse_values
    norm = Normalize(vmin=min(all_values), vmax=max(all_values))
    scalar = ScalarMappable(norm=norm, cmap="viridis")
    scalar.set_array([])
    colorbar = fig.colorbar(scalar, ax=axes, shrink=0.78, pad=0.025)
    colorbar.set_label("exact cell volume / centroid pointweight (mm³)")
    fig.suptitle(
        "Actual analytic section samples from hollow-cylinder cell geometry\n"
        "Rin=8 mm, Rout=10 mm; colors show exact annular-wedge weights; dots show projected exact centroids",
        fontsize=13,
    )
    fig.text(
        0.5,
        0.025,
        "Fixture-only local evidence. The declared ≥3-span/≥3-centroid precondition is not universal crack guidance; "
        "process-zone, physical, and whole-domain qualification remain unevaluated.",
        ha="center",
        va="bottom",
        fontsize=9,
    )
    fig.savefig(path, dpi=170, facecolor="white")
    plt.close(fig)


def main() -> int:
    REVIEW_DIR.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    fine = make_analytic_hollow_cylinder_fixture(0.5)
    coarse = make_analytic_hollow_cylinder_fixture(2.0)
    fine_report = qualify_wall_input(fine)
    coarse_report = qualify_wall_input(coarse)
    elapsed = time.perf_counter() - started
    if elapsed > 30.0:
        raise RuntimeError(f"analytic qualification exceeded 30 s: {elapsed:.3f}")
    historical = historical_h8_assessment(nominal_wall_mm=2.0)
    _json_write(REVIEW_DIR / "thin_wall_qualification_inputs.json", {
        "schema": "pgfrag.thin_wall_resolution_fixture_bundle.v1",
        "fixture_disclaimer": "Analytic geometry only; no CAD kernel, solver, or original source writer was invoked.",
        "cases": [fine, coarse],
        "historical_h8_comparison": historical,
    })
    _json_write(REVIEW_DIR / "thin_wall_qualification_report.json", {
        "schema": "pgfrag.thin_wall_resolution_review.v1",
        "status": "LOCAL_PRECONDITION_MET_FOR_H05_FIXTURE__H2_UNDERRESOLVED__NO_PHYSICAL_OR_GLOBAL_PASS",
        "cases": [fine_report, coarse_report],
        "historical_h8_comparison": historical,
        "execution": {
            "elapsed_seconds": elapsed,
            "time_limit_seconds": 30,
            "max_fixture_cells": 96,
            "max_fixture_segments": 128,
            "cad_invoked": False,
            "solver_invoked": False,
            "original_ingress_source_modified": False,
        },
    })
    _write_section_figure(fine, coarse, REVIEW_DIR / "actual_section_samples.png")
    _write_review_note(fine_report, coarse_report, historical, elapsed)
    print(json.dumps({
        "review_dir": str(REVIEW_DIR),
        "elapsed_seconds": elapsed,
        "h0_5": fine_report["resolution"]["local_status"],
        "h0_5_min_intervals": fine_report["resolution"]["min_observed_actual_thickness_intervals"],
        "h2": coarse_report["resolution"]["local_status"],
        "h2_min_intervals": coarse_report["resolution"]["min_observed_actual_thickness_intervals"],
        "historical_h8": historical["status"],
    }, ensure_ascii=False, indent=2))
    return 0


def _write_review_note(fine: dict, coarse: dict, historical: dict, elapsed: float) -> None:
    body = f"""# Thin-wall resolution guard review

Status: **local analytic fixture qualification only**. The utility is isolated under `pgfrag/experiments/wall_resolution_guard_20261006/`; it does not modify or import the old ingress worker, invoke CAD, or launch a solver.

The old source is the bounded 96-cell ROI ingress experiment at `pgfrag/experiments/forge_volume_ingress_20261005/`. Its recovered local ledger reconciles 8.01676271854445 mm³ across the 96 cells and 1.68352017089433e-05 kg at 2100 kg/m³. The independent audit reports exact R=8/10 mm annular-cylinder ROI integration and local point-weight checks. Those records do not establish global minimum wall thickness or thickness-direction layer count; the producer also remained blocked at report finalization and did not promote a candidate.

This checker's exact local fixture is a hollow cylinder with `t=2 mm`, `Rin=8 mm`, `Rout=10 mm`, and a 2 mm axial patch. `h=0.5 mm` uses {fine['inputs']['cell_count']} exact annular-wedge cells. Every declared local normal sample traverses {fine['resolution']['min_observed_actual_thickness_intervals']} actual cell intervals and contains {fine['resolution']['min_observed_distinct_centroid_layers']} distinct projected centroid layers. Its result is **{fine['resolution']['local_status']}**. The `h=2 mm` fixture has {coarse['inputs']['cell_count']} cells but only {coarse['resolution']['min_observed_actual_thickness_intervals']} actual interval and centroid layer, so its result is **{coarse['resolution']['local_status']}**. The weights are exact sector volumes `0.5 × (Rout²−Rin²) × Δtheta × Δz`; no `h³` substitute is accepted.

The historical `h=8 mm` against the source's nominal 2 mm wall parameter is explicitly **{historical['status']}** (`h/t={historical['h_over_t']:.1f}`). This is a coarse nominal-parameter comparison, not a recovered local mesh measurement or a global wall certificate.

The declared three-interval and three-centroid thresholds are fixture-specific engineering preconditions. They are not a universal crack-resolution rule. Process-zone adequacy is a separate required check. Even the fine fixture reports `physical_status=NOT_QUALIFIED`, `physical_pass=false`, and `global_pass=false`; missing full-domain coverage can only produce a global UNKNOWN. The line-segment adapter uses exact quadratic intersections with inner and outer cylinder radii plus axial bounds; it rejects a long segment whose two material endpoints are joined through the hollow cavity.

The qualification run took {elapsed:.4f} s (limit 30 s), with {fine['inputs']['cell_count']} / {coarse['inputs']['cell_count']} cells and {fine['inputs']['segment_count']} / {coarse['inputs']['segment_count']} segments. See `thin_wall_qualification_report.json`, `thin_wall_qualification_inputs.json`, and `actual_section_samples.png` for raw review artifacts.

Run tests with:

```powershell
third_party/vesselforge-runtime/venv/Scripts/python.exe -m pytest -q pgfrag/experiments/wall_resolution_guard_20261006/test_thin_wall_qualification.py
```
"""
    (REVIEW_DIR / "REVIEW.md").write_text(body, encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())

