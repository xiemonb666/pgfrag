"""Render the pinned twelve-step opening-control result from its saved JSON values."""
import argparse
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def resolve(path):
    value = Path(path)
    return value if value.is_absolute() else ROOT / value


def number(value, label):
    result = float(value)
    if not math.isfinite(result):
        raise RuntimeError(f"BLOCKED: nonfinite saved {label}")
    return result


def load_saved(case_path, pins_path, token):
    case = json.loads(case_path.read_text(encoding="utf-8"))
    pins = json.loads(pins_path.read_text(encoding="utf-8"))
    auth, go = case["execution_authorization"], pins["root_go"]
    if auth.get("authorized") is not True or auth.get("authorized_by") != "ROOT": raise RuntimeError("BLOCKED: root render authorization is absent")
    if (pins.get("execution_authorized") is not True or go.get("authorized") is not True
            or not token or token != go.get("token")):
        raise RuntimeError("BLOCKED: root render token does not match SOURCE_PINS")
    refs = pins["references"]
    for name, path in (("renderer", Path(__file__).resolve()), ("case", case_path)):
        if refs[name].get("sha256") != sha(path): raise RuntimeError(f"BLOCKED: pinned {name} source changed")
    if any(not resolve(rec["path"]).is_file() or sha(resolve(rec["path"])) != rec["sha256"] for rec in refs.values()):
        raise RuntimeError("BLOCKED: one or more source references changed")
    result_ref = case["input_result"]
    if (not result_ref.get("sha256") or result_ref["path"] != pins["input_result"].get("path") or result_ref["sha256"] != pins["input_result"].get("sha256")):
        raise RuntimeError("BLOCKED: actual result SHA-256 is not consistently pinned")
    result_path = resolve(result_ref["path"])
    if not result_path.is_file() or sha(result_path) != result_ref["sha256"]: raise RuntimeError("BLOCKED: actual result is absent or changed")
    result = json.loads(result_path.read_text(encoding="utf-8"))
    if (result.get("schema") != result_ref["expected_schema"] or result.get("status") != result_ref["expected_status"] or result.get("classification") != result_ref["expected_classification"]):
        raise RuntimeError("BLOCKED: result schema, diagnostic status, or classification changed")
    if result.get("source_hashes") != case["producer_sources"]["result_source_hashes"]: raise RuntimeError("BLOCKED: result producer source bindings differ from the fixed case")
    trajectory = result.get("trajectory", [])
    contract = case["trajectory_contract"]
    fractions = contract["expected_Q_fraction_Qstar"]
    if len(trajectory) != contract["step_count"] or result.get("trajectory_fraction_Qstar") != fractions: raise RuntimeError("BLOCKED: saved trajectory length or Qstar schedule changed")
    interfaces = result.get("saved_input_binding", {}).get("interfaces_reconstructed_from_saved_Hplus_Hminus", [])
    base = trajectory[0].get("points", []) if trajectory else []
    if len(interfaces) != 16 or len(base) != 16: raise RuntimeError("BLOCKED: expected sixteen saved interface QPs")
    point_ids = [tuple(point.get("point_id", [])) for point in base]
    if any(len(point_id) != 2 for point_id in point_ids) or len(set(point_ids)) != 16: raise RuntimeError("BLOCKED: saved QP identities are incomplete or duplicated")
    if any(point_ids[i][1] != face.get("quadrature_id") for i, face in enumerate(interfaces)): raise RuntimeError("BLOCKED: interface/QP order mismatch")
    alphas = [number(p["analytic_alpha"], "alpha") for p in base]
    for i, row in enumerate(trajectory):
        if row.get("step_index") != i or not math.isclose(number(row["Q_fraction_of_Qstar"], "Q fraction"), fractions[i], rel_tol=0.0, abs_tol=1e-14):
            raise RuntimeError(f"BLOCKED: saved step order changed at row {i}")
        number(row["Q_m"], "Q")
        points = row.get("points", [])
        if len(points) != 16 or [tuple(p.get("point_id", [])) for p in points] != point_ids: raise RuntimeError(f"BLOCKED: saved QP identity/order changed at step {i}")
        for j, p in enumerate(points):
            xyz = [number(v, "point coordinate") for v in p["point_xyz_m"]]
            if len(xyz) != 3 or abs(xyz[1]) > 1e-14 or p["point_xyz_m"] != base[j]["point_xyz_m"]: raise RuntimeError(f"BLOCKED: saved point coordinates changed at step {i}")
            number(p["opening_m"], "opening"); number(p["traction_pa"], "traction")
            damage = number(p["state"]["damage"], "damage")
            if not 0.0 <= damage <= 1.0: raise RuntimeError(f"BLOCKED: saved damage outside [0,1] at step {i}")
        gamma = row.get("cohesive_generalized_internal_force_N", [])
        if len(gamma) != 93: raise RuntimeError(f"BLOCKED: cohesive force vector length changed at step {i}")
        for key in ("control_applied_force_positive_opening_N", "cohesive_stored_energy_j",
                    "cohesive_dissipated_energy_j", "cohesive_work_increment_independent_law_j"):
            number(row[key], key)
        number(gamma[82], "cohesive gamma[82]")
    partial = trajectory[3]
    if partial.get("partial_damage_point_count", 0) <= 0: raise RuntimeError("BLOCKED: saved step 3 has no partial damage")
    if any(trajectory[i].get("all_points_fully_damaged") is not True for i in contract["full_damage_step_indices_to_annotate"]):
        raise RuntimeError("BLOCKED: a declared full-damage step is not fully damaged")
    order = sorted(range(16), key=lambda i: (alphas[i], point_ids[i][1]))
    selected = [order[0], order[len(order) // 2], order[-1]]
    return case, result, result_path, result_ref["sha256"], selected


def render(case, result, output_path, selected):
    # Plotting imports stay below every authorization, source, result, and history gate.
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    trajectory = result["trajectory"]
    steps = list(range(12))
    cmap = plt.get_cmap("viridis", 12)
    step_norm = matplotlib.colors.Normalize(vmin=-0.5, vmax=11.5)
    fig = plt.figure(figsize=(13, 9), constrained_layout=True)
    grid = fig.add_gridspec(5, 3, height_ratios=(0.16, 0.22, 1, 1, 0.14), width_ratios=(1, 1, 0.065))
    header, legend_ax = fig.add_subplot(grid[0, :]), fig.add_subplot(grid[1, :]); a0, a1 = fig.add_subplot(grid[2, 0]), fig.add_subplot(grid[2, 1])
    step_cax = fig.add_subplot(grid[2, 2]); a2, a3 = fig.add_subplot(grid[3, 0]), fig.add_subplot(grid[3, 1])
    damage_cax, footer = fig.add_subplot(grid[3, 2]), fig.add_subplot(grid[4, :])
    for ax in (header, legend_ax, footer):
        ax.axis("off")
    header.text(0.5, 0.72, "FIXED_FACE | TOY OPENING HISTORY", ha="center", va="center", weight="bold")
    header.text(0.5, 0.18, "NOT A FRACTURE OR FRAGMENTATION MODEL | NO_GROWTH | NO_MATERIAL_QUALIFICATION", ha="center", va="center")
    qp_colors, qp_markers, qp_styles = ("#0072B2", "#D55E00", "#009E73"), ("o", "s", "^"), ("-", "--", ":")
    qp_names = ("Minimum α QP", "Upper-median α QP", "Maximum α QP")
    handles = []
    for rank, index in enumerate(selected):
        p = trajectory[0]["points"][index]
        label = f"{qp_names[rank]} ({p['point_id'][1]}, α={float(p['analytic_alpha']):.3g})"
        handles.append(Line2D([], [], color=qp_colors[rank], linestyle=qp_styles[rank], marker=qp_markers[rank], markerfacecolor="#999999", label=label))
    handles.extend([Line2D([], [], color="#222222", marker="o", label="Cohesive γ[82] internal force"),
                    Line2D([], [], color="#CC79A7", marker="s", label="Total positive-opening control force"),
                    Line2D([], [], color="#0072B2", marker="o", label="Cohesive stored energy"),
                    Line2D([], [], color="#D55E00", marker="s", label="Cohesive dissipated energy"),
                    Line2D([], [], color="#009E73", marker="^", label="Cumulative independent cohesive law work")])
    legend_ax.legend(handles=handles, loc="center", ncol=4, frameon=False, fontsize=7.5)
    for rank, index in enumerate(selected):
        opening = [float(row["points"][index]["opening_m"]) for row in trajectory]; traction = [float(row["points"][index]["traction_pa"]) for row in trajectory]
        a0.plot(opening, traction, color=qp_colors[rank], linestyle=qp_styles[rank], linewidth=1.1)
        a0.scatter(opening, traction, c=steps, cmap=cmap, norm=step_norm, marker=qp_markers[rank], s=38,
                   edgecolors="black", linewidths=0.3)
    a0.set_title("(a) Three real QP traction-opening histories"); a0.set_xlabel("Normal opening δ (m)"); a0.set_ylabel("Traction t (Pa)"); a0.grid(alpha=0.2)
    step_mappable = plt.cm.ScalarMappable(norm=step_norm, cmap=cmap)
    fig.colorbar(step_mappable, cax=step_cax, ticks=steps, label="Saved step index (zero-based)"); Q = [float(row["Q_m"]) for row in trajectory]
    gamma = [float(row["cohesive_generalized_internal_force_N"][82]) for row in trajectory]; total = [float(row["control_applied_force_positive_opening_N"]) for row in trajectory]
    a1.plot(Q, gamma, color="#222222", linestyle="-", linewidth=1.0); a1.scatter(Q, gamma, c=steps, cmap=cmap, norm=step_norm, marker="o", edgecolors="#222222", linewidths=0.3)
    a1.plot(Q, total, color="#CC79A7", linestyle="--", linewidth=1.0); a1.scatter(Q, total, c=steps, cmap=cmap, norm=step_norm, marker="s", edgecolors="#222222", linewidths=0.3)
    a1.set_title("(b) Cohesive γ[82] vs total control"); a1.set_xlabel("Saved prescribed control Q (m)"); a1.set_ylabel("Signed force (N)"); a1.grid(alpha=0.2)
    stored = [float(row["cohesive_stored_energy_j"]) * 1e6 for row in trajectory]
    dissipated = [float(row["cohesive_dissipated_energy_j"]) * 1e6 for row in trajectory]
    law_work, cumulative = [], 0.0
    for row in trajectory: cumulative += float(row["cohesive_work_increment_independent_law_j"]); law_work.append(cumulative * 1e6)
    a2.plot(steps, stored, "o-", color="#0072B2"); a2.plot(steps, dissipated, "s-", color="#D55E00"); a2.plot(steps, law_work, "^-", color="#009E73")
    a2.set_title("(c) Cohesive-only energy and law work"); a2.set_xlabel("Saved step index (zero-based)"); a2.set_ylabel("Cohesive energy / path work (µJ)"); a2.set_xticks(steps); a2.grid(alpha=0.2)
    map_points = trajectory[3]["points"]
    x = [float(p["point_xyz_m"][0]) for p in map_points]; z = [float(p["point_xyz_m"][2]) for p in map_points]; damage = [float(p["state"]["damage"]) for p in map_points]
    sc = a3.scatter(x, z, c=damage, cmap="viridis", vmin=0.0, vmax=1.0, s=54, edgecolors="black", linewidths=0.35)
    a3.set_title("(d) Saved step 3 partial damage (16 QPs)"); a3.set_xlabel("x (m)"); a3.set_ylabel("z (m)")
    a3.set_xlim(-case["geometry"]["radius_m"], 0.0); a3.set_ylim(-case["geometry"]["length_m"] / 2.0, case["geometry"]["length_m"] / 2.0)
    a3.set_aspect("equal", adjustable="box"); a3.grid(alpha=0.2)
    fig.colorbar(sc, cax=damage_cax, label="Saved damage d (—)", ticks=(0, 0.25, 0.5, 0.75, 1))
    footer.text(0.5, 0.5, "All QPs fully damaged at steps 8, 10 and 11 (step 10 is Q=0); saved cohesive-only values, fixed face, no growth.", ha="center", va="center", fontsize=8)
    if output_path.exists(): raise RuntimeError("BLOCKED: refusing to overwrite saved visual output")
    fig.savefig(output_path, dpi=140, format="png", facecolor="white"); plt.close(fig)
    if output_path.stat().st_size > case["limits"]["max_combined_output_bytes"]: raise RuntimeError("BLOCKED: PNG alone exceeds output byte limit")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--pins", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--root-go", required=True)
    args = parser.parse_args()
    case, result, result_path, expected_sha, selected = load_saved(args.case, args.pins, args.root_go)
    render(case, result, args.output, selected)
    if sha(result_path) != expected_sha:
        raise RuntimeError("BLOCKED: input result changed during render")
    print(json.dumps({"status": "RENDERED_SAVED_OPENING_HISTORY_ONLY", "result_sha256": expected_sha,
                      "output_bytes": args.output.stat().st_size}, sort_keys=True))


if __name__ == "__main__":
    raise SystemExit(main())
