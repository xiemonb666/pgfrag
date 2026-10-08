#!/usr/bin/env python3
"""Plot existing B raw CSV only; never launches or changes a native solver run."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "pgfrag"
SNAPSHOT = BASE / "tmp" / "review_20261004" / "pd_loaded_diagnostics" / "snapshots" / "affine_threshold_b_run01"
DIRECTED = SNAPSHOT / "affine_native_bond_states.csv"
CANONICAL = SNAPSHOT / "affine_canonical_pair_states.csv"
OUT = BASE / "tmp" / "review_20261005" / "native_resume" / "native_b_affine_threshold_audit.png"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    canonical_rows = read_rows(CANONICAL)
    directed_rows = read_rows(DIRECTED)
    state_order = ["rest", "below", "above"]
    colors = {"rest": "#54718c", "below": "#31817b", "above": "#c97a2b"}
    values_by_state: dict[str, list[tuple[float, bool]]] = defaultdict(list)
    for row in canonical_rows:
        state = row["state_id"]
        normalized = float(row["actual_stretch"]) / float(row["native_epsilon_c"])
        broken = row["active_after"].strip().lower() == "false"
        values_by_state[state].append((normalized, broken))

    directed_by_state_pair: dict[str, dict[tuple[int, int], dict[tuple[int, int], float]]] = defaultdict(lambda: defaultdict(dict))
    for row in directed_rows:
        state = row["state_id"]
        source = int(row["source_point_id"])
        neighbor = int(row["neighbor_point_id"])
        pair = (int(row["canonical_i_id"]), int(row["canonical_j_id"]))
        normalized = float(row["actual_stretch"]) / float(row["native_epsilon_c"])
        directed_by_state_pair[state][pair][(source, neighbor)] = normalized

    pair_counts = {state: len(values_by_state[state]) for state in state_order}
    broken_counts = {state: sum(broken for _, broken in values_by_state[state]) for state in state_order}
    bins = np.linspace(-0.02, 1.62, 67)
    plt.rcParams.update({"font.size": 10, "axes.titleweight": "semibold", "figure.facecolor": "white"})
    fig, (ax_hist, ax_recip) = plt.subplots(1, 2, figsize=(13.2, 5.3), constrained_layout=True)

    for state in ("rest", "below"):
        values = [value for value, _ in values_by_state[state]]
        ax_hist.hist(values, bins=bins, color=colors[state], alpha=0.58,
                     label=f"{state}: {len(values):,} pairs")
    above_active = [value for value, broken in values_by_state["above"] if not broken]
    above_broken = [value for value, broken in values_by_state["above"] if broken]
    ax_hist.hist(above_active, bins=bins, color="#da9c53", alpha=0.72,
                 label=f"above, still active: {len(above_active):,}")
    ax_hist.hist(above_broken, bins=bins, color="#b64d45", alpha=0.88,
                 label=f"above, native transition: {len(above_broken):,}")
    ax_hist.axvline(1.0, color="#25364a", linestyle="--", linewidth=1.5,
                    label=r"native threshold $epsilon_c$")
    ax_hist.set_title("Canonical pair stretch under prescribed affine map")
    ax_hist.set_xlabel(r"actual stretch / native $epsilon_c$")
    ax_hist.set_ylabel("canonical pair count")
    ax_hist.set_xlim(-0.02, 1.62)
    ax_hist.grid(axis="y", color="#dce2e9", linewidth=0.7)
    ax_hist.legend(frameon=False, loc="upper left", fontsize=8.5)
    ax_hist.text(0.98, 0.95,
                 "native rule: break when stretch > εc and fail_permit",
                 transform=ax_hist.transAxes, ha="right", va="top", fontsize=8.5,
                 bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "#d6dde5"})

    max_abs_reciprocal_delta = 0.0
    for state in state_order:
        xs = []
        ys = []
        for pair, directions in directed_by_state_pair[state].items():
            low, high = pair
            forward = directions[(low, high)]
            reverse = directions[(high, low)]
            xs.append(forward)
            ys.append(reverse)
            max_abs_reciprocal_delta = max(max_abs_reciprocal_delta, abs(forward - reverse))
        ax_recip.scatter(xs, ys, s=8, alpha=0.16, color=colors[state],
                         label=f"{state}: {len(xs):,} reciprocal pairs")
    ax_recip.plot([-0.03, 1.58], [-0.03, 1.58], color="#25364a", linestyle="--", linewidth=1.25)
    ax_recip.set_title("Directed pair measures: low-id → high-id vs reverse")
    ax_recip.set_xlabel(r"stretch / $epsilon_c$ (low-id → high-id)")
    ax_recip.set_ylabel(r"stretch / $epsilon_c$ (high-id → low-id)")
    ax_recip.set_xlim(-0.03, 1.58)
    ax_recip.set_ylim(-0.03, 1.58)
    ax_recip.set_aspect("equal", adjustable="box")
    ax_recip.grid(color="#dce2e9", linewidth=0.7)
    ax_recip.legend(frameon=False, loc="upper left", fontsize=8.5)
    ax_recip.text(0.97, 0.05,
                  f"max |forward − reverse| = {max_abs_reciprocal_delta:.1e}\n"
                  f"break counts (rest / below / above): 0 / 0 / {broken_counts['above']:,} pairs",
                  transform=ax_recip.transAxes, ha="right", va="bottom", fontsize=8.5,
                  bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "#d6dde5"})

    fig.suptitle("PG-Frag B: native CriticalStretch threshold diagnostic", fontsize=15, fontweight="bold")
    fig.text(0.5, -0.015,
             "F = diag(1, 1, 1 + η), applied to every point from pristine coordinates; no time integration, no free fracture. "
             "OSB / CriticalStretch prototype values are unvalidated. Source: existing B raw CSVs.",
             ha="center", va="top", fontsize=8.4, color="#4c5866")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"figure={OUT}")
    print(f"canonical_pair_counts={pair_counts}")
    print(f"broken_pair_counts={broken_counts}")
    print(f"max_abs_reciprocal_delta={max_abs_reciprocal_delta:.17g}")


if __name__ == "__main__":
    main()
