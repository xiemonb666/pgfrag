from collections import defaultdict
import json
import os
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(HERE / ".mplconfig"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle

ROOT = HERE.parents[4]
EXP = ROOT / "pgfrag/experiments/moving_front_provider_20261006"
result = json.loads((EXP / "moving_front_check_result.run02.json").read_text(encoding="utf-8"))
case_doc = json.loads((EXP / "moving_front_cases.json").read_text(encoding="utf-8"))
assert result["status"] == "PASS" and len(result["states"]) == 2
cases = {float(c["front_x_m"]): c for c in case_doc["moving_front_cases"]}
plots = []
for state in result["states"]:
    case = cases[float(state["front_x_m"])]
    audit = state["query_audit"]
    assert len(audit) == 25
    rows = [q for q in audit if q["point_m"][2] == 0.0]
    plots.append((state, case, rows))
assert sum(len(rows) for _, _, rows in plots) <= 50
xs, ys = [], []
for state, case, rows in plots:
    xs += [q["point_m"][0] for q in rows]
    ys += [q["point_m"][1] for q in rows]
    xlo, xhi = case["finite_crack"]["x_interval_m"]
    xs += [xlo, xhi]
    ys += [case["finite_crack"]["y_m"]]
    r = state["measured_core_radius_m"]
    fx, fy = case["front_segment"]["x_m"], case["front_segment"]["y_m"]
    xs += [fx-r, fx+r]
    ys += [fy-r, fy+r]
pad_x = max((max(xs)-min(xs))*0.025, 0.00025)
pad_y = max((max(ys)-min(ys))*0.10, 0.00025)
fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharex=True, sharey=True)
for ax, (state, case, rows) in zip(axes, plots):
    scale = 1000.0
    xlo, xhi = case["finite_crack"]["x_interval_m"]
    cy = case["finite_crack"]["y_m"]
    ax.plot([xlo*scale, xhi*scale], [cy*scale, cy*scale], color="#343a40", lw=2.2, label="crack projection Γ")
    fx, fy = case["front_segment"]["x_m"], case["front_segment"]["y_m"]
    ax.scatter([fx*scale], [fy*scale], marker="D", s=52, color="#d62728", zorder=5, label="front projection")
    for radius, color, style, label in [(state["measured_core_radius_m"], "#9467bd", ":", "measured core"), *[(r, "#1f77b4", "--", "probe radius") for r in state["core_probe_radii_m"]]]:
        ax.add_patch(Circle((fx*scale, fy*scale), radius*scale, fill=False, ec=color, ls=style, lw=1.15, alpha=0.9, label=label))
    for status, marker, color in [("OK", "o", "#2474a6"), ("UNKNOWN", "x", "#d1495b")]:
        group = [q for q in rows if q["status"] == status]
        if group:
            ax.scatter([q["point_m"][0]*scale for q in group], [q["point_m"][1]*scale for q in group], marker=marker, s=23, color=color, zorder=4, label=status)
    coloc = defaultdict(list)
    for q in rows:
        p = q["point_m"]
        if q.get("trace_side"):
            coloc[(p[0], p[1])].append(q)
    for idx, ((x, y), traces) in enumerate(coloc.items()):
        if len(traces) < 2:
            continue
        ax.scatter([x*scale], [y*scale], s=88, facecolors="none", edgecolors="#6a3d9a", lw=1.5, zorder=6)
        stem = re.sub(r"-(?:plus|minus)$", "", sorted(q["label"] for q in traces)[0])
        sides = "/".join(sorted({q["trace_side"] for q in traces}))
        ax.annotate(f"{stem}\n{sides} traces share XY", (x*scale, y*scale), xytext=(8, 16 if idx % 2 == 0 else -26), textcoords="offset points", fontsize=7, color="#542788", arrowprops={"arrowstyle":"-", "color":"#6a3d9a", "lw":0.7})
    ok = sum(q["status"] == "OK" for q in rows)
    unknown = sum(q["status"] == "UNKNOWN" for q in rows)
    ax.set_title(f"front x={fx*scale:g} mm | z=0 rows={len(rows)} (OK {ok}, UNKNOWN {unknown})")
    ax.set_xlabel("x [mm]")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim((min(xs)-pad_x)*scale, (max(xs)+pad_x)*scale)
    ax.set_ylim((min(ys)-pad_y)*scale, (max(ys)+pad_y)*scale)
    ax.grid(True, lw=0.45, alpha=0.35)
axes[0].set_ylabel("y [mm]")
handles = [Line2D([], [], marker="o", ls="", color="#2474a6", label="OK query"), Line2D([], [], marker="x", ls="", color="#d1495b", label="UNKNOWN query"), Line2D([], [], color="#343a40", lw=2.2, label="Γ XY projection"), Line2D([], [], marker="D", ls="", color="#d62728", label="front projection"), Line2D([], [], color="#9467bd", ls=":", label="measured core"), Line2D([], [], color="#1f77b4", ls="--", label="saved probe radii")]
fig.legend(handles=handles, loc="upper center", ncol=6, frameon=False, bbox_to_anchor=(0.5, 0.96), fontsize=8)
fig.suptitle("Moving-front saved query audit: actual XY locations only", y=1.02, fontsize=12)
fig.text(0.5, 0.01, "Only saved z=0 query rows are plotted. Purple rings identify co-located PLUS/MINUS trace queries; crack and front are XY projections. No field interpolation or new solve.", ha="center", fontsize=8)
fig.tight_layout(rect=[0, 0.04, 1, 0.90])
fig.savefig(HERE / "saved_query_audit_xy.png", dpi=180, bbox_inches="tight")
plt.close(fig)