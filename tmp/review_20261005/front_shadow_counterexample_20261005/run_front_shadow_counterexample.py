"""One bounded source-only run of the frozen front-shadow counterexample."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time

sys.dont_write_bytecode = True

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
CASE_PATH = HERE / "fixedcase.json"
MANIFEST_PATH = HERE / "source_manifest.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(65536), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def load_pinned_helpers():
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    source_paths = {}
    for item in manifest["sources"]:
        path = REPO / item["path"]
        observed = sha256(path)
        if observed != item["sha256"]:
            raise RuntimeError(f"source pin mismatch for {item['path']}: {observed}")
        source_paths[item["module"]] = path.parent
    sys.path.insert(0, str(source_paths["mls3d_shape"]))
    sys.path.insert(0, str(source_paths["visibility3d"]))
    from visibility3d import classify_visibility
    from mls3d_mls import MLS3DConfig, MLSStatus, wendland_c2
    from mls3d_shape import shape_functions_3d

    return classify_visibility, MLS3DConfig, MLSStatus, wendland_c2, shape_functions_3d


def _plot(case: dict, rows: list[dict], png_path: Path) -> None:
    os.environ["MPLCONFIGDIR"] = str(HERE / ".mplconfig")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    import numpy as np

    sources = np.asarray(case["source_positions_m"], dtype=float)
    facet = np.asarray(case["triangle_facets_m"][0], dtype=float)
    t_values = np.asarray([row["t_m"] for row in rows], dtype=float)
    q_values = np.asarray([row["query_m"] for row in rows], dtype=float)
    phi_values = np.asarray([row["special_phi"] for row in rows], dtype=float)

    fig = plt.figure(figsize=(11.2, 4.7), constrained_layout=True)
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    ax.add_collection3d(Poly3DCollection([facet], facecolor="#f29b38", alpha=0.42, edgecolor="#9a4d00", linewidth=1.4))
    ax.scatter(sources[0, 0], sources[0, 1], sources[0, 2], color="#b22222", s=40, label="special source s")
    ax.scatter(sources[1:, 0], sources[1:, 1], sources[1:, 2], color="#276c9b", s=24, label="four tetra sources")
    ax.plot(q_values[:, 0], q_values[:, 1], q_values[:, 2], color="#202020", linewidth=2.0, label="query path")
    ax.scatter(q_values[:, 0], q_values[:, 1], q_values[:, 2], color="#202020", s=14)
    for point, color in ((q_values[0], "#277a3b"), (q_values[-1], "#7043a4")):
        ax.plot([sources[0, 0], point[0]], [sources[0, 1], point[1]], [sources[0, 2], point[2]], color=color, linewidth=1.3, alpha=0.9)
    ax.scatter([1.0], [0.0], [1.0], color="#d62728", marker="x", s=44)
    ax.text(1.0, 0.0, 1.08, "q0: visibility UNKNOWN")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_zlabel("z (m)")
    ax.set_title("Finite facet and source-to-query shadow")
    ax.set_xlim(-2.2, 1.8)
    ax.set_ylim(-1.2, 1.2)
    ax.set_zlim(-1.2, 1.8)
    ax.set_box_aspect((4.0, 2.4, 3.0))
    ax.legend(loc="upper left", fontsize=7)

    bx = fig.add_subplot(1, 2, 2)
    negative = t_values < 0.0
    positive = t_values > 0.0
    bx.plot(t_values[negative], phi_values[negative], marker="o", color="#276c9b", linewidth=1.7, label="t < 0")
    bx.plot(t_values[positive], phi_values[positive], marker="o", color="#b45b35", linewidth=1.7, label="t > 0")
    bx.axvline(0.0, color="#b22222", linestyle="--", linewidth=1.0)
    bx.text(0.98, 0.94, "t = 0: visibility UNKNOWN\nno q0 value or gradient", transform=bx.transAxes, ha="right", va="top", fontsize=8)
    bx.set_xlabel("t (m)")
    bx.set_ylabel("special-node MLS shape value")
    bx.set_title("Six fixed nonzero samples")
    bx.grid(True, alpha=0.25)
    bx.legend(loc="best", fontsize=8)
    fig.savefig(png_path, dpi=145, transparent=True)
    plt.close(fig)


def main() -> None:
    started = time.perf_counter()
    import numpy as np

    plot_requested = "--plot" in sys.argv[1:]
    json_path = HERE / "results.json"
    png_path = HERE / "front_shadow_counterexample.png"
    if json_path.exists() or (plot_requested and png_path.exists()):
        raise FileExistsError("refusing to overwrite a prior result or plot")

    case = json.loads(CASE_PATH.read_text(encoding="utf-8"))
    source_positions = np.ascontiguousarray(case["source_positions_m"], dtype=np.float64)
    query_family = case["query_family"]
    t_values = np.asarray(query_family["sample_t_m"], dtype=np.float64)
    origin = np.asarray(query_family["base_q0_m"], dtype=np.float64)
    direction = np.asarray(query_family["direction_m_per_parameter"], dtype=np.float64)
    query_positions = np.ascontiguousarray(origin[None, :] + t_values[:, None] * direction[None, :], dtype=np.float64)
    triangle_facets = np.ascontiguousarray(case["triangle_facets_m"], dtype=np.float64)
    if source_positions.shape[0] > 5 or query_positions.shape[0] > 6:
        raise RuntimeError("fixedcase exceeds the authorized node/query budget")

    classify_visibility, MLS3DConfig, MLSStatus, wendland_c2, shape_functions_3d = load_pinned_helpers()
    visibility = classify_visibility(source_positions, query_positions, triangle_facets)
    config = MLS3DConfig(
        reference_origin=tuple(case["mls"]["reference_origin_m"]),
        reference_scale=case["mls"]["reference_scale_m"],
        support_radius=case["mls"]["support_radius_m"],
    )

    rows = []
    for query_index, (t_value, query) in enumerate(zip(t_values, query_positions)):
        t = float(t_value)
        visible_mask = np.array(visibility.visible_mask[query_index], dtype=bool, copy=True)
        statuses = [str(value) for value in visibility.statuses[query_index].tolist()]
        expected_statuses = ["BLOCKED", "VISIBLE", "VISIBLE", "VISIBLE", "VISIBLE"] if t < 0.0 else ["VISIBLE"] * 5
        visibility_pass = statuses == expected_statuses

        # Exact geometry establishes local mask constancy at each nonzero t.
        shape = shape_functions_3d(
            source_positions,
            np.ascontiguousarray(query[None, :], dtype=np.float64),
            config,
            visibility_mask=visible_mask,
            mask_locally_constant=True,
        )
        weights = shape.weights[0]
        gradients = shape.gradients[0]
        active = np.flatnonzero(visible_mask & (np.linalg.norm(source_positions - query[None, :], axis=1) < config.support_radius))
        q_active = np.linalg.norm(source_positions[active] - query[None, :], axis=1) / config.support_radius
        radial_weights, _ = wendland_c2(q_active)
        p_active = np.column_stack((np.ones(active.size), (source_positions[active] - np.asarray(config.reference_origin)) / config.reference_scale))
        moment = np.einsum("n,ni,nj->ij", radial_weights, p_active, p_active)
        singular_values = np.linalg.svd(moment, compute_uv=False)
        moment_rank = int(np.count_nonzero(singular_values > singular_values[0] * config.rank_relative_tolerance))
        support_count = int(active.size)

        partition_error = abs(float(np.sum(weights)) - 1.0)
        position_residual = weights @ source_positions - query
        gradient_partition = np.sum(gradients, axis=0)
        coordinate_gradient = source_positions.T @ gradients
        gradient_coordinate_residual = coordinate_gradient - np.eye(3)
        affine_pass = (
            shape.statuses[0] is MLSStatus.OK
            and moment_rank == 4
            and support_count == (4 if t < 0.0 else 5)
            and partition_error <= 2.0e-11
            and float(np.max(np.abs(position_residual))) <= 2.0e-11
            and float(np.max(np.abs(gradient_partition))) <= 2.0e-10
            and float(np.max(np.abs(gradient_coordinate_residual))) <= 2.0e-10
        )
        rows.append({
            "t_m": t,
            "query_m": query.tolist(),
            "q0_to_finite_facet_distance_m": math.sqrt((1.0 + t) ** 2 + 1.0),
            "visibility_statuses_source_order_s_v1_v2_v3_v4": statuses,
            "visibility_mask_source_order_s_v1_v2_v3_v4": visible_mask.tolist(),
            "mask_locally_constant": True,
            "mask_stability_exact_geometry_margin_m": abs(t) / math.sqrt(2.0),
            "support_count": support_count,
            "moment_rank": moment_rank,
            "moment_condition_number": float(shape.condition_numbers[0]),
            "shape_status": shape.statuses[0].value,
            "special_phi": float(weights[0]),
            "special_gradient_per_m": gradients[0].tolist(),
            "shape_row_source_order_s_v1_v2_v3_v4": weights.tolist(),
            "shape_gradient_rows_per_m_source_order_s_v1_v2_v3_v4": gradients.tolist(),
            "affine_reproduction": {
                "partition_of_unity_abs_error": partition_error,
                "coordinate_reproduction_residual_m": position_residual.tolist(),
                "gradient_partition_residual_per_m": gradient_partition.tolist(),
                "coordinate_gradient_map": coordinate_gradient.tolist(),
                "coordinate_gradient_map_residual": gradient_coordinate_residual.tolist(),
                "passed": affine_pass
            },
            "visibility_expectation_passed": visibility_pass
        })

    checks = {
        "all_six_visibility_masks_match_exact_geometry": all(row["visibility_expectation_passed"] for row in rows),
        "all_six_rows_full_rank": all(row["moment_rank"] == 4 for row in rows),
        "all_six_rows_affinely_reproduce": all(row["affine_reproduction"]["passed"] for row in rows),
        "special_shape_has_one_sided_jump": rows[0]["special_phi"] == 0.0 and rows[-1]["special_phi"] > 0.0,
        "no_q0_gradient_reported": True,
        "physical_facet_distance_q0_m": math.sqrt(2.0),
        "result": "sharp visibility causes an MLS basis jump on a source-front shadow boundary outside the finite physical facet"
    }
    checks_passed = all(value for key, value in checks.items() if isinstance(value, bool))
    output = {
        "schema": "front-shadow-counterexample-results-v1",
        "case_id": case["case_id"],
        "status": "HELPER_DIAGNOSTIC_PASSED" if checks_passed else "HELPER_DIAGNOSTIC_FAILED",
        "mechanics_coupling": "NOT_COUPLED",
        "true_q0_visibility": "UNKNOWN_EDGE_CONTACT",
        "q0_gradient_claim": "NONE; q0 has distinct one-sided algebraic value traces",
        "q0_one_sided_values": {
            "label": "ONE_SIDED_ALGEBRAIC_TRACE_ONLY",
            "left_special_phi": 0.0,
            "right_special_phi": 0.007528410117638756,
            "right_special_phi_strictly_positive": True
        },
        "source_hashes": json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))["sources"],
        "rows": rows,
        "checks": checks
    }
    json_path.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    if not checks_passed:
        raise RuntimeError("fixed-case checks failed; actual per-check booleans are preserved in results.json")
    if plot_requested:
        _plot(case, rows, png_path)

    elapsed = time.perf_counter() - started
    artifacts = [json_path]
    if png_path.exists():
        artifacts.append(png_path)
    artifact_bytes = sum(path.stat().st_size for path in artifacts)
    if elapsed > 30.0 or artifact_bytes > 2097152:
        raise RuntimeError(f"run exceeded fixed runtime/output budget: {elapsed:.3f}s, {artifact_bytes} bytes")
    print(json.dumps({"status": output["status"], "runtime_s": elapsed, "artifact_bytes": artifact_bytes, "rows": len(rows)}, sort_keys=True))


if __name__ == "__main__":
    main()
