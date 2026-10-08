"""Run a bounded illustrative history and save review CSV/JSON/PNGs."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
PGFRAG_ROOT = HERE.parents[1]
OUT = PGFRAG_ROOT / "tmp" / "review_20261006" / "cohesive_scalar_law_20261006"
sys.path.insert(0, str(HERE))

from cohesive_scalar_law import (  # noqa: E402
    BilinearCohesiveMaterial,
    initial_state,
    traction_pa,
    update_state,
)


def _path(material: BilinearCohesiveMaterial) -> list[float]:
    df = material.deltaf_m
    # 66 load states, 32 unload states, 19 reload/post-peak states, 32 final
    # unload states: 149 samples total including the initial state.
    load = [0.90 * df * i / 65 for i in range(1, 66)]
    unload = [0.90 * df * (1.0 - i / 32) for i in range(1, 33)]
    reload = [0.96 * df * i / 19 for i in range(1, 20)]
    final_unload = [0.96 * df * (1.0 - i / 32) for i in range(1, 33)]
    return [0.0, *load, *unload, *reload, *final_unload]


def main() -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    OUT.mkdir(parents=True, exist_ok=True)
    material = BilinearCohesiveMaterial(
        k0_pa_per_m=1.0e12,
        tpeak_pa=1.0e7,
        gc_j_per_m2=100.0,
    )
    state = initial_state(material)
    rows: list[dict[str, float | int]] = []
    peak_traction = 0.0
    min_dissipation_increment = 0.0
    previous_dissipation = state.dissipated_energy_j_per_m2
    for step, opening in enumerate(_path(material)):
        if step:
            state = update_state(material, state, opening)
        value = traction_pa(material, state)
        peak_traction = max(peak_traction, value)
        min_dissipation_increment = min(
            min_dissipation_increment,
            state.dissipated_energy_j_per_m2 - previous_dissipation,
        )
        previous_dissipation = state.dissipated_energy_j_per_m2
        rows.append(
            {
                "step": step,
                "opening_m": state.opening_m,
                "kappa_m": state.kappa_m,
                "damage": state.damage,
                "traction_pa": value,
                "stored_energy_j_per_m2": state.stored_energy_j_per_m2,
                "dissipated_energy_j_per_m2": state.dissipated_energy_j_per_m2,
                "external_work_j_per_m2": state.external_work_j_per_m2,
            }
        )

    csv_path = OUT / "cohesive_history.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    ledger_error = max(
        abs(row["external_work_j_per_m2"] - row["stored_energy_j_per_m2"] - row["dissipated_energy_j_per_m2"])
        for row in rows
    )
    receipt = {
        "artifact_kind": "scalar_mode_i_bilinear_cohesive_review_fixture",
        "implementation": "cohesive_scalar_law.py",
        "status": "PASS" if ledger_error < 1.0e-8 and min_dissipation_increment >= -1.0e-10 else "FAIL",
        "fixture_warning": "Illustrative software fixture only; not a ceramic material calibration.",
        "units": {
            "opening": "m",
            "k0": "Pa/m",
            "tpeak": "Pa",
            "gc_and_energy_per_area": "J/m^2",
        },
        "parameters": {
            "k0_pa_per_m": material.k0_pa_per_m,
            "tpeak_pa": material.tpeak_pa,
            "gc_j_per_m2": material.gc_j_per_m2,
            "delta0_m": material.delta0_m,
            "deltaf_m": material.deltaf_m,
        },
        "limits": {
            "compression": "UNKNOWN / unsupported",
            "mixed_mode": "UNKNOWN / unsupported",
            "rate_dependence": "not represented",
            "front_coupling_or_solver": "not implemented",
        },
        "reference_basis": [
            {
                "citation": "Zhang, Z. and Paulino, G.H. (2005), Cohesive zone modeling of dynamic failure in homogeneous and functionally graded materials, International Journal of Plasticity 21, 1195-1254, doi:10.1016/j.ijplas.2004.06.009.",
                "url": "https://paulino.princeton.edu/journal_papers/2005/IJP_05_CohesiveZoneModeling.pdf",
                "use": "Primary source for a traction-based bilinear cohesive law and use of maximum separation history for irreversibility; the present API narrows it to scalar opening-only.",
            },
            {
                "citation": "Alfano, M. et al. (2009), Mode I fracture of adhesive joints using tailored cohesive zone models, International Journal of Fracture 157, 193-204, doi:10.1007/s10704-009-9344-9.",
                "url": "https://paulino.princeton.edu/journal_papers/2009/IJF_09_ModeIFracture.pdf",
                "use": "Equations (5)-(6) state the initial linear branch, linear softening branch, and Gc = sigma_cr delta_f / 2 area relation; used as a directly inspectable formula reference.",
                "note": "This is a peer-reviewed implementation/formulation reference, not the original first proposal.",
            },
            {
                "citation": "Hillerborg, A., Modeer, M. and Petersson, P.-E. (1976), Analysis of crack formation and crack growth in concrete by means of fracture mechanics and finite elements, Cement and Concrete Research 6, 773-782.",
                "url": "https://www.codil.or.kr/filebank/original/RK/OTMCRK500148/OTMCRK500148.pdf",
                "use": "Foundational fictitious-crack/fracture-energy concept; cited for context, not as the precise modern bilinear parameterization.",
            },
        ],
        "run_receipt": {
            "sample_count": len(rows),
            "sample_limit": 256,
            "samples_bounded": len(rows) <= 256,
            "peak_traction_pa": peak_traction,
            "final_damage": rows[-1]["damage"],
            "final_kappa_m": rows[-1]["kappa_m"],
            "final_dissipated_energy_j_per_m2": rows[-1]["dissipated_energy_j_per_m2"],
            "maximum_abs_energy_ledger_error_j_per_m2": ledger_error,
            "minimum_dissipation_increment_j_per_m2": min_dissipation_increment,
        },
    }
    (OUT / "run_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")

    steps = [row["step"] for row in rows]
    openings_um = [row["opening_m"] * 1.0e6 for row in rows]
    kappas_um = [row["kappa_m"] * 1.0e6 for row in rows]
    tractions_mpa = [row["traction_pa"] / 1.0e6 for row in rows]

    fig, ax = plt.subplots(figsize=(7.0, 4.5), constrained_layout=True)
    ax.plot(openings_um, tractions_mpa, color="#146c94", linewidth=1.8, marker="o", markersize=2.1)
    ax.set(xlabel="Opening δ (µm)", ylabel="Traction (MPa)", title="Bilinear cohesive traction history")
    ax.grid(True, alpha=0.25)
    fig.savefig(OUT / "traction_history.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.0, 4.5), constrained_layout=True)
    ax.plot(steps, openings_um, label="current δ", color="#e07a1f", linewidth=1.8)
    ax.plot(steps, kappas_um, label="history κ = max δ", color="#5a6fbb", linewidth=1.8)
    ax.set(xlabel="State update", ylabel="Opening (µm)", title="Opening history and irreversible maximum")
    ax.legend(frameon=False)
    ax.grid(True, alpha=0.25)
    fig.savefig(OUT / "opening_history.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.0, 4.5), constrained_layout=True)
    ax.plot(steps, [r["stored_energy_j_per_m2"] for r in rows], label="stored ψ", linewidth=1.8)
    ax.plot(steps, [r["dissipated_energy_j_per_m2"] for r in rows], label="dissipated D", linewidth=1.8)
    ax.plot(steps, [r["external_work_j_per_m2"] for r in rows], label="external work W", linewidth=1.8)
    ax.set(xlabel="State update", ylabel="Energy per area (J/m²)", title="Cohesive energy ledger")
    ax.legend(frameon=False)
    ax.grid(True, alpha=0.25)
    fig.savefig(OUT / "energy_ledger.png", dpi=160)
    plt.close(fig)

    print(json.dumps({"status": receipt["status"], "csv": str(csv_path), "json": str(OUT / "run_receipt.json"), "samples": len(rows), "ledger_error": ledger_error}, indent=2))


if __name__ == "__main__":
    main()
