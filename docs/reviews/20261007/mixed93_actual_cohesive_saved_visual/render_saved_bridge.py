"""Render only the pinned mixed93 cohesive result after root authorization."""
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


def load_saved(case_path, pins_path, token):
    case = json.loads(case_path.read_text(encoding="utf-8"))
    pins = json.loads(pins_path.read_text(encoding="utf-8"))
    auth, go = case["execution_authorization"], pins["root_go"]
    if auth.get("authorized") is not True or auth.get("authorized_by") != "ROOT":
        raise RuntimeError("BLOCKED: root render authorization is absent")
    if (pins.get("execution_authorized") is not True or go.get("authorized") is not True
            or not token or token != go.get("token")):
        raise RuntimeError("BLOCKED: root render token does not match SOURCE_PINS")
    refs = pins.get("references", {})
    if refs.get("renderer", {}).get("sha256") != sha(Path(__file__).resolve()) or refs.get("case", {}).get("sha256") != sha(case_path):
        raise RuntimeError("BLOCKED: renderer or case differs from its reviewed source pin")
    worker = case["producer_worker_source"]
    if pins.get("producer_worker_source") != worker or sha(resolve(worker["path"])) != worker["sha256"]:
        raise RuntimeError("BLOCKED: numeric worker source differs from its fixed pin")
    ref = case["input_result"]
    if not ref.get("sha256") or ref["path"] != pins["input_result"].get("path") or ref["sha256"] != pins["input_result"].get("sha256"):
        raise RuntimeError("BLOCKED: actual-result SHA-256 is not consistently pinned")
    result_path = resolve(ref["path"])
    if not result_path.is_file() or sha(result_path) != ref["sha256"]:
        raise RuntimeError("BLOCKED: accepted numeric result is absent or changed")
    result = json.loads(result_path.read_text(encoding="utf-8"))
    if (result.get("schema") != ref["expected_schema"] or result.get("status") != ref["expected_status"]
            or result.get("classification") != ref["expected_classification"]):
        raise RuntimeError("BLOCKED: result schema, diagnostic status, or classification changed")
    gates = result.get("gates", {})
    if not gates or any(gate.get("passed") is not True for gate in gates.values()):
        raise RuntimeError("BLOCKED: numeric worker did not pass every saved gate")
    evidence = result.get("evidence", {})
    interface = evidence.get("interface", [])
    baseline = evidence.get("baseline", {}).get("points", [])
    coupled = evidence.get("coupled", {}).get("points", [])
    count = case["geometry"]["expected_saved_qps"]
    if not (len(interface) == len(baseline) == len(coupled) == count):
        raise RuntimeError("BLOCKED: saved interface/baseline/coupled point counts differ")
    delta0 = float(case["toy_law"]["delta0_m"])
    for i, (face, before, after) in enumerate(zip(interface, baseline, coupled)):
        qid = face.get("quadrature_id")
        if before.get("point_id", [None, None])[-1] != qid or after.get("point_id", [None, None])[-1] != qid:
            raise RuntimeError(f"BLOCKED: saved point order mismatch at row {i}")
        xyz = face.get("point_xyz_m", [])
        if len(xyz) != 3 or not all(math.isfinite(float(v)) for v in xyz) or abs(float(xyz[1])) > 1.0e-14:
            raise RuntimeError(f"BLOCKED: saved QP coordinate invalid at row {i}")
        area = float(face.get("one_side_area_weight_m2", 0.0))
        if face.get("normal_unit_lower_to_upper") != case["geometry"]["normal_unit_lower_to_upper"] or not math.isfinite(area) or area <= 0.0:
            raise RuntimeError(f"BLOCKED: saved QP normal or one-side area is invalid at row {i}")
        for point in (before, after):
            opening, traction = float(point["opening_m"]), float(point["traction_pa"])
            if not math.isfinite(opening + traction) or not (0.0 < opening < delta0) or traction < 0.0:
                raise RuntimeError(f"BLOCKED: saved point is outside the virgin elastic opening range at row {i}")
    solution = result.get("solution", {})
    if not all(math.isfinite(float(solution[name])) for name in ("q_change_norm_m", "internal_force_norm_n")):
        raise RuntimeError("BLOCKED: saved q-change/force summary is absent or non-finite")
    return case, result, result_path, ref["sha256"]


def render(case, result, output_path):
    # Matplotlib is deliberately imported only after pins, root GO, and saved data pass.
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    interface = result["evidence"]["interface"]
    baseline = result["evidence"]["baseline"]["points"]
    coupled = result["evidence"]["coupled"]["points"]
    x = [float(row["point_xyz_m"][0]) for row in interface]
    z = [float(row["point_xyz_m"][2]) for row in interface]
    bopen = [float(row["opening_m"]) for row in baseline]
    copen = [float(row["opening_m"]) for row in coupled]
    law, geo = case["toy_law"], case["geometry"]
    radius, length = geo["radius_m"], geo["length_m"]
    delta0, k0 = law["delta0_m"], law["k0_pa_per_m"]
    observed_max = max(bopen + copen)
    opening_max = 1.1 * observed_max
    fig = plt.figure(figsize=(12, 9), constrained_layout=True)
    grid = fig.add_gridspec(5, 3, height_ratios=(0.14, 0.12, 1.0, 1.0, 0.12), width_ratios=(1.0, 1.0, 0.18))
    header, legend_ax = fig.add_subplot(grid[0, :]), fig.add_subplot(grid[1, :]); a0, a1 = fig.add_subplot(grid[2, 0]), fig.add_subplot(grid[2, 1])
    cbar_ax = fig.add_subplot(grid[2, 2]); a2, a3 = fig.add_subplot(grid[3, 0]), fig.add_subplot(grid[3, 1])
    note_ax, footer = fig.add_subplot(grid[3, 2]), fig.add_subplot(grid[4, :])
    header.axis("off"); legend_ax.axis("off"); note_ax.axis("off"); footer.axis("off")
    header.text(0.5, 0.70, "FIXED_CRACK_TOY_ELASTIC_DIAGNOSTIC", ha="center", va="center", weight="bold")
    header.text(0.5, 0.18, "NO_GROWTH  |  NO_MATERIAL_CALIBRATION", ha="center", va="center")
    handles = [Line2D([], [], marker="o", linestyle="none", color="#444444", label="Baseline QPs"),
               Line2D([], [], marker="s", linestyle="none", color="#444444", label="Coupled-final QPs"),
               Line2D([], [], color="#555555", label="Virgin law t=K0δ"),
               Line2D([], [], marker="o", color="#1B7837", label="Internal share +t n"),
               Line2D([], [], marker="x", linestyle="--", color="#762A83", label="Physical upper face −t n"),
               Line2D([], [], marker="s", linestyle=":", color="#2166AC", label="Physical lower face +t n")]
    legend_ax.legend(handles=handles, loc="center", ncol=3, frameon=False, fontsize=8)
    s0 = a0.scatter(x, z, c=bopen, cmap="viridis", vmin=0.0, vmax=opening_max, marker="o", s=64, edgecolors="black", linewidths=0.35)
    a1.scatter(x, z, c=copen, cmap="viridis", vmin=0.0, vmax=opening_max, marker="s", s=64, edgecolors="black", linewidths=0.35)
    for ax, title in ((a0, "(a) Saved baseline opening"), (a1, "(b) Saved coupled-final opening")):
        ax.set_title(title); ax.set_xlabel("x (m)"); ax.set_ylabel("z (m)")
        ax.set_xlim(-radius, 0.0); ax.set_ylim(-length / 2.0, length / 2.0); ax.set_aspect("equal", adjustable="box")
        ax.grid(alpha=0.2)
    fig.colorbar(s0, cax=cbar_ax, label="Normal opening δ (m)"); a2.set_title("(c) Saved traction-opening points")
    a2.plot([0.0, observed_max], [0.0, k0 * observed_max], color="#555555")
    a2.scatter(bopen, [float(p["traction_pa"]) for p in baseline], marker="o", color="#444444")
    a2.scatter(copen, [float(p["traction_pa"]) for p in coupled], marker="s", color="#444444")
    data_tmax = max(float(p["traction_pa"]) for p in baseline + coupled)
    a2.set_xlim(0.0, opening_max); a2.set_ylim(0.0, 1.1 * max(k0 * observed_max, data_tmax))
    a2.set_xlabel("Normal opening δ (m)"); a2.set_ylabel("Cohesive traction t (Pa)"); a2.grid(alpha=0.2)
    sol = result["solution"]
    a3.set_title("(d) Coupled-final signed normal tractions")
    normal = interface[0]["normal_unit_lower_to_upper"]
    qp = list(range(1, len(interface) + 1)); internal, upper, lower = [], [], []
    project = lambda vector: sum(float(f) * float(n) for f, n in zip(vector, normal))
    for face, point in zip(interface, coupled):
        area = float(face["one_side_area_weight_m2"])
        internal.append(project(point["primitive_positive_face_internal_force_n"]) / area)
        upper.append(project(point["physical_upper_face_traction_n_per_m2"])); lower.append(project(point["physical_lower_face_traction_n_per_m2"]))
    a3.plot(qp, internal, "o-", color="#1B7837")
    a3.plot(qp, upper, "x--", color="#762A83")
    a3.plot(qp, lower, "s:", color="#2166AC")
    a3.axhline(0.0, color="black", linewidth=0.7); a3.set_xlabel("Saved QP order"); a3.set_ylabel("Normal projection (Pa)"); a3.grid(alpha=0.2)
    note_ax.text(0.02, 0.92, f"Observed δ max={observed_max:.3g} m\nδ0={delta0:.3g} m\nActual points below δ0; virgin-only.\nNo softening is shown.\n\nSaved ‖q_c−q_0‖₂={sol['q_change_norm_m']:.3g} m\nSaved ‖fΓ‖₂={sol['internal_force_norm_n']:.3g} N",
                 va="top", fontsize=8, wrap=True)
    footer.text(0.5, 0.5, "Saved 16-QP values only; fixed crack, no growth, and no material calibration.", ha="center", va="center", fontsize=9)
    if output_path.exists(): raise RuntimeError("BLOCKED: refusing to overwrite the saved visual output")
    fig.savefig(output_path, dpi=140, format="png", facecolor="white"); plt.close(fig)
    if output_path.stat().st_size > case["limits"]["max_combined_output_bytes"]: raise RuntimeError("BLOCKED: PNG alone exceeds the fixed output byte budget")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--pins", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--root-go", required=True)
    args = parser.parse_args()
    case, result, result_path, expected_sha = load_saved(args.case, args.pins, args.root_go)
    render(case, result, args.output)
    if sha(result_path) != expected_sha:
        raise RuntimeError("BLOCKED: actual result changed during render")
    print(json.dumps({"status": "RENDERED_SAVED_DIAGNOSTIC_ONLY", "result_sha256": expected_sha,
                      "output_bytes": args.output.stat().st_size}, sort_keys=True))


if __name__ == "__main__":
    raise SystemExit(main())
