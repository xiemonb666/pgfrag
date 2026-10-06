"""Render-only view of focused_run01/numeric_result.json; never imports or calls the law/assembler."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RUN = ROOT / "pgfrag/tmp/review_20261006/cohesive_interface_assembly_20261006/focused_run01"
BRANCHES = ("elastic", "soft", "unload")
LABELS = {"elastic": "Elastic load", "soft": "Softening load", "unload": "Unload"}
COLORS = {"elastic": "#2878b5", "soft": "#c44536", "unload": "#4b8b3b"}


def main() -> None:
    source = Path(sys.argv[1]) if len(sys.argv) > 1 else RUN / "numeric_result.json"
    output = Path(sys.argv[2]) if len(sys.argv) > 2 else RUN / "toy_numeric_diagnostic.png"
    result = json.loads(source.read_text(encoding="utf-8"))
    if result.get("status") != "PASS" or "not vessel/physical validation" not in result.get("kind", ""):
        raise ValueError("expected the passed toy-checker result JSON")
    fig, (traction_ax, fd_ax) = plt.subplots(1, 2, figsize=(13, 5.5), constrained_layout=True)
    for key in BRANCHES:
        row = result["baseline"][key]
        x, y = row["opening_m"], row["traction_pa"]
        traction_ax.scatter(x, y, s=72, color=COLORS[key], label=LABELS[key], zorder=3)
        traction_ax.annotate(f"{LABELS[key]}\nk = {row['k_alg_pa_per_m']:.3g} Pa/m", (x, y), xytext=(7, 7), textcoords="offset points", fontsize=8)
    traction_ax.set(xlabel="Normal opening δ [m]", ylabel="Normal traction t [Pa]", title="Trial baseline points (history-dependent)")
    traction_ax.grid(True, alpha=0.25)
    traction_ax.legend(frameon=False, loc="best")

    y = list(range(len(BRANCHES)))
    error = [max(result["fd_branches"][key]["max_abs_error_n_per_m"], 1.0e-30) for key in BRANCHES]
    limit = [result["fd_branches"][key]["predeclared_limit_n_per_m"] for key in BRANCHES]
    fd_ax.barh([i - 0.18 for i in y], error, height=0.34, color="#2878b5", label="FD max |error|")
    fd_ax.barh([i + 0.18 for i in y], limit, height=0.34, color="#b8c6d1", label="Predeclared limit")
    fd_ax.set_xscale("log")
    fd_ax.set_yticks(y)
    fd_ax.set_yticklabels([LABELS[key] for key in BRANCHES])
    fd_ax.set(xlabel="Directional tangent error / limit [N/m] (log scale)", title="Finite-difference tangent checks")
    fd_ax.grid(True, axis="x", which="both", alpha=0.25)
    fd_ax.legend(frameon=False, loc="lower right")
    broken = result["fully_broken"]
    fig.suptitle("Mode-I cohesive interface — toy mathematical diagnostic only")
    fig.text(0.5, -0.015, f"Fully broken: dissipated {broken['dissipated_energy_j']:.6g} J / expected {broken['expected_gc_energy_j']:.6g} J; ledger residual {broken['energy_ledger_abs_residual_j']:.3g} J", ha="center", fontsize=9)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(output)


if __name__ == "__main__":
    main()
