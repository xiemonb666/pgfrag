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
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    interface = result["evidence"]["interface"]
    baseline = result["evidence"]["baseline"]["points"]
    coupled = result["evidence"]["coupled"]["points"]
    x = [float(row["point_xyz_m"][0]) for row in interface]
    z = [float(row["point_xyz_m"][2]) for row in interface]
    bopen = [float(row["opening_m"]) for row in baseline]
    copen = [float(row["opening_m"]) for row in coupled]
    law, radius, length = case["toy_law"], case["geometry"]["radius_m"], case["geometry"]["length_m"]
    delta0, k0, tpeak = law["delta0_m"], law["k0_pa_per_m"], law["tpeak_pa"]
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    a0, a1, a2, a3 = axes.ravel()
    s0 = a0.scatter(x, z, c=bopen, cmap="viridis", vmin=0.0, vmax=delta0, marker="o", s=64, edgecolors="black", linewidths=0.35)
    a1.scatter(x, z, c=copen, cmap="viridis", vmin=0.0, vmax=delta0, marker="s", s=64, edgecolors="black", linewidths=0.35)
    for ax, title in ((a0, "(a) Saved baseline opening"), (a1, "(b) Saved coupled-final opening")):
        ax.set_title(title); ax.set_xlabel("x (m)"); ax.set_ylabel("z (m)")
        ax.set_xlim(-radius, 0.0); ax.set_ylim(-length / 2.0, length / 2.0); ax.set_aspect("equal", adjustable="box")
        ax.grid(alpha=0.2)
    fig.colorbar(s0, ax=[a0, a1], label="Normal opening δ (m)", fraction=0.045, pad=0.04)
    a2.set_title("(c) Saved traction-opening points")
    a2.plot([0.0, delta0], [0.0, k0 * delta0], color="#555555", label="Virgin law t=K0δ; shown only to δ0")
    a2.axvline(delta0, color="#333333", linestyle="--", linewidth=1.0, label=f"δ0={delta0:.3g} m")
    a2.scatter(bopen, [float(p["traction_pa"]) for p in baseline], marker="o", color="#0072B2", label="Baseline QPs")
    a2.scatter(copen, [float(p["traction_pa"]) for p in coupled], marker="s", color="#D55E00", label="Coupled-final QPs")
    data_tmax = max(float(p["traction_pa"]) for p in baseline + coupled)
    a2.set_xlim(0.0, delta0 * 1.08); a2.set_ylim(0.0, max(tpeak, data_tmax) * 1.12)
    a2.set_xlabel("Normal opening δ (m)"); a2.set_ylabel("Cohesive traction t (Pa)"); a2.grid(alpha=0.2); a2.legend(fontsize=8)
    sol = result["solution"]
    a2.text(0.03, 0.96, f"Saved ‖q_c−q_0‖₂={sol['q_change_norm_m']:.3g} m\nSaved ‖fΓ‖₂={sol['internal_force_norm_n']:.3g} N",
            transform=a2.transAxes, va="top", fontsize=8)
    a3.set_title("(d) Coupled-final signed normal tractions")
    normal = interface[0]["normal_unit_lower_to_upper"]
    qp = list(range(1, len(interface) + 1)); internal, upper, lower = [], [], []
    for face, point in zip(interface, coupled):
        area = float(face["one_side_area_weight_m2"])
        internal.append(sum(float(f) * float(n) for f, n in zip(point["primitive_positive_face_internal_force_n"], normal)) / area)
        upper.append(sum(float(f) * float(n) for f, n in zip(point["physical_upper_face_traction_n_per_m2"], normal)))
        lower.append(sum(float(f) * float(n) for f, n in zip(point["physical_lower_face_traction_n_per_m2"], normal)))
    a3.plot(qp, internal, "o-", label="Internal share +t n")
    a3.plot(qp, upper, "x--", label="Physical upper face −t n")
    a3.plot(qp, lower, "s:", label="Physical lower face +t n")
    a3.axhline(0.0, color="black", linewidth=0.7); a3.set_xlabel("Saved QP order"); a3.set_ylabel("Normal projection (Pa)")
    a3.grid(alpha=0.2); a3.legend(fontsize=8)
    fig.suptitle("FIXED_CRACK_TOY_ELASTIC_DIAGNOSTIC\nNO_GROWTH  |  NO_MATERIAL_CALIBRATION", y=0.98)
    fig.subplots_adjust(top=0.87, bottom=0.14, hspace=0.36, wspace=0.3)
    fig.text(0.5, 0.035, "Saved 16-QP result only; no crack extension or post-δ0 softening is plotted.", ha="center", fontsize=9)
    if output_path.exists():
        raise RuntimeError("BLOCKED: refusing to overwrite the saved visual output")
    fig.savefig(output_path, dpi=140, format="png", facecolor="white")
    plt.close(fig)
    if output_path.stat().st_size > case["limits"]["max_combined_output_bytes"]:
        raise RuntimeError("BLOCKED: PNG alone exceeds the fixed output byte budget")


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
