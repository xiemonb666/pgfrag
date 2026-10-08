"""Render-only refresh of the retained 96-cell audit plots.

Reads the already-frozen audit JSON and CSVs. It does not load CAD, invoke the
producer, recompute geometry, or alter any prior output. New PNGs and a render
receipt go into a new child directory so the earlier figures remain intact.
"""

from __future__ import annotations

import ast
import csv
import hashlib
import json
import math
import traceback
from pathlib import Path


REVIEW_DIR = Path(__file__).resolve().parent
AUDIT_DIR = REVIEW_DIR / "audit_retry_hashcase_plotcompat_20261005"
OUTPUT_DIR = AUDIT_DIR / "render_only_revised_colorbar_20261005"
AUDIT_JSON = AUDIT_DIR / "PILOT_ACTUAL_AUDIT.json"
CELL_CSV = AUDIT_DIR / "PILOT_CELL_AUDIT.csv"
POINT_CSV = AUDIT_DIR / "PILOT_POINT_AUDIT.csv"
EXPECTED_SHA256 = {
    AUDIT_JSON.name: "de33b3028cc1ad4b8d735b1f2fff9d7632fb15e24fd6b5fc6a6f515020ae217f",
    CELL_CSV.name: "8517032040d22b1a8d22b3f331431360b14c9c630675dbd0f5a8dbe8a9a469c5",
    POINT_CSV.name: "0c76db6e93a4c1f31449aecae7cf99f8ca32a894fafdc9101b11fb87899d3d74",
}
RECEIPT_PATH = OUTPUT_DIR / "RENDER_ONLY_RECEIPT.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_receipt(receipt: dict) -> None:
    RECEIPT_PATH.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    receipt = {
        "schema": "pgfrag.forge_volume_ingress_render_only_receipt.v1",
        "status": "PENDING",
        "source_mode": "existing_audit_json_and_csv_only",
        "geometry_recomputed": False,
        "cad_or_producer_invoked": False,
        "inputs": {},
        "outputs": {},
        "exception": None,
    }
    try:
        for path in (AUDIT_JSON, CELL_CSV, POINT_CSV):
            actual = sha256(path)
            expected = EXPECTED_SHA256[path.name]
            receipt["inputs"][path.name] = {
                "path": str(path),
                "sha256": actual,
                "expected_sha256": expected,
                "matches_expected": actual == expected,
            }
            if actual != expected:
                raise RuntimeError(f"frozen render input SHA-256 mismatch: {path.name}")

        audit = json.loads(AUDIT_JSON.read_text(encoding="utf-8"))
        if audit.get("status") != "GEOMETRY_DIAGNOSTIC_AUDITED_NOT_PROMOTED":
            raise RuntimeError("retained numeric audit status is not the expected diagnostic")
        if audit.get("all_checks_pass") is not True:
            raise RuntimeError("retained numeric audit check verdict is not all-pass")

        with CELL_CSV.open("r", newline="", encoding="utf-8") as stream:
            cells = list(csv.DictReader(stream))
        with POINT_CSV.open("r", newline="", encoding="utf-8") as stream:
            points = list(csv.DictReader(stream))
        if len(cells) != 96 or len(points) != 80:
            raise RuntimeError("retained CSV row counts do not match the frozen audit")

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib import colors
        from matplotlib.patches import Rectangle
        from matplotlib.ticker import MaxNLocator

        roi = audit["analytic_roi"]["roi_aabb_mm"]
        xmin, xmax, ymin, ymax, zmin, zmax = (float(value) for value in roi)
        inner_radius = float(audit["analytic_roi"]["inner_radius_mm"])
        outer_radius = float(audit["analytic_roi"]["outer_radius_mm"])
        cell_size = float(audit["analytic_roi"]["cell_size_mm"])
        indices = sorted(cells, key=lambda row: int(row["grid_cell_index"]))

        # XY view: use a dedicated colorbar axis so it cannot overlap the map.
        fraction_by_row = []
        for row in indices:
            aabb = ast.literal_eval(row["cell_aabb_mm"])
            area = float(row["analytic_area_mm2"])
            fraction_by_row.append(area / ((float(aabb[1]) - float(aabb[0])) *
                                           (float(aabb[3]) - float(aabb[2]))))
        norm = colors.Normalize(0.0, max(max(fraction_by_row, default=0.0), 1e-15))
        cmap = matplotlib.colormaps["viridis"]
        points_by_k: dict[int, list[tuple[float, float]]] = {k: [] for k in range(4)}
        for row in points:
            k = int(row["grid_k"])
            points_by_k[k].append((float(row["x_mm_from_m"]), float(row["y_mm_from_m"])))

        fig = plt.figure(figsize=(12.5, 8.0), constrained_layout=True)
        grid = fig.add_gridspec(
            2, 3, width_ratios=(1.0, 1.0, 0.055), height_ratios=(1.0, 1.0),
            wspace=0.22, hspace=0.32,
        )
        axes = [
            fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1]),
            fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1]),
        ]
        colorbar_ax = fig.add_subplot(grid[:, 2])
        angles = [2.0 * 3.141592653589793 * n / 720 for n in range(721)]
        circle_cos = [math.cos(t) for t in angles]
        circle_sin = [math.sin(t) for t in angles]
        for k, ax in enumerate(axes):
            for row in indices:
                if int(row["grid_k"]) != k:
                    continue
                aabb = ast.literal_eval(row["cell_aabb_mm"])
                fraction = float(row["analytic_area_mm2"]) / (
                    (float(aabb[1]) - float(aabb[0])) *
                    (float(aabb[3]) - float(aabb[2]))
                )
                ax.add_patch(Rectangle(
                    (float(aabb[0]), float(aabb[2])),
                    float(aabb[1]) - float(aabb[0]),
                    float(aabb[3]) - float(aabb[2]),
                    facecolor=cmap(norm(fraction)), edgecolor="white", linewidth=0.7,
                ))
            for radius, linestyle, color in (
                (outer_radius, "-", "black"), (inner_radius, "--", "dimgray")
            ):
                ax.plot([radius * x for x in circle_cos],
                        [radius * y for y in circle_sin],
                        linestyle=linestyle, color=color, linewidth=1.0)
            slab_points = points_by_k[k]
            if slab_points:
                ax.scatter(
                    [point[0] for point in slab_points],
                    [point[1] for point in slab_points],
                    s=26, facecolor="white", edgecolor="black", linewidth=0.65, zorder=5,
                )
            zmid = zmin + (k + 0.5) * cell_size
            ax.set_title(f"z slab {k}, centroid z={zmid:.2f} mm")
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymin, ymax)
            ax.set_aspect("equal", adjustable="box")
            ax.set_xlabel("x (mm)")
            ax.set_ylabel("y (mm)")
        mapper = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
        mapper.set_array([])
        fig.colorbar(
            mapper, cax=colorbar_ax,
            label="analytic annular material area fraction per XY cell",
        )
        fig.suptitle("Retained local geometry; producer final report blocked; result not promoted")
        map_path = OUTPUT_DIR / "PILOT_XY_CELL_MAP_grid_colorbar.png"
        fig.savefig(map_path, dpi=170)
        plt.close(fig)

        # Per-cell chart is reconstructed only from already-recorded CSV values.
        cell_ids = [int(row["grid_cell_index"]) for row in indices]
        analytic_vol = [float(row["analytic_volume_mm3"]) for row in indices]
        raw_vol = [float(row["raw_common_volume_mm3"]) for row in indices]
        point_vol = [float(row["delivered_point_volume_mm3"]) for row in indices]
        analytic_mass = [float(row["analytic_mass_kg"]) for row in indices]
        point_mass = [float(row["delivered_point_mass_kg"]) for row in indices]
        point_delta = [float(row["points_minus_analytic_mm3"]) for row in indices]
        raw_delta = [float(row["raw_minus_analytic_mm3"]) for row in indices]
        fig, axes = plt.subplots(3, 1, figsize=(13, 9), sharex=True, constrained_layout=True)
        axes[0].plot(cell_ids, analytic_vol, label="analytic cell volume")
        axes[0].scatter(cell_ids, raw_vol, s=12, label="raw B-rep common")
        axes[0].scatter(cell_ids, point_vol, s=12, marker="x",
                        label="sum of delivered point volumes")
        axes[0].set_ylabel("volume (mm³)")
        axes[0].legend(fontsize=8)
        axes[1].plot(cell_ids, analytic_mass, label="rho × analytic volume")
        axes[1].scatter(cell_ids, point_mass, s=12, label="sum of point mass")
        axes[1].set_ylabel("mass (kg)")
        axes[1].legend(fontsize=8)
        axes[2].axhline(0.0, color="black", linewidth=0.8)
        axes[2].plot(cell_ids, point_delta, label="delivered − analytic")
        axes[2].plot(cell_ids, raw_delta, label="raw B-rep − analytic")
        axes[2].set_xlabel("stable cell index (k, then j, then i)")
        axes[2].set_ylabel("signed volume delta (mm³)")
        axes[2].legend(fontsize=8)
        for ax in axes:
            ax.grid(True, alpha=0.25)
            ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=12))
        fig.suptitle("Per-cell checks; local diagnostic only, not promoted")
        delta_path = OUTPUT_DIR / "PILOT_CELL_VOLUME_MASS_DELTA_render_only.png"
        fig.savefig(delta_path, dpi=170)
        plt.close(fig)

        receipt["status"] = "COMPLETE"
        receipt["outputs"] = {
            map_path.name: {"path": str(map_path), "sha256": sha256(map_path)},
            delta_path.name: {"path": str(delta_path), "sha256": sha256(delta_path)},
        }
    except Exception as exc:  # Keep render outcome separate from numeric audit JSON.
        receipt["status"] = "FAILED"
        receipt["exception"] = {
            "type": type(exc).__name__,
            "message": str(exc),
            "traceback": traceback.format_exc(),
        }
        raise
    finally:
        write_receipt(receipt)


if __name__ == "__main__":
    main()
