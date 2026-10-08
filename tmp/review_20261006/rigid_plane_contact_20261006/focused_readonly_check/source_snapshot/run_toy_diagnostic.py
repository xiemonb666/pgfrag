"""One bounded validation and 1-D Velocity-Verlet bounce diagnostic."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
from PIL import Image, ImageDraw

from rigid_plane_contact import RigidPlanePointContact


HERE = Path(__file__).resolve().parent
DEFAULT_OUT = HERE.parents[1] / "tmp" / "review_20261006" / "rigid_plane_contact_20261006"


def _expect_rejected(label: str, call: Any) -> str:
    try:
        call()
    except (TypeError, ValueError, FloatingPointError):
        return label
    raise AssertionError(f"expected validation failure: {label}")


FOCUSED_CHECK_DETAILS: dict[str, Any] = {}


def run_self_checks() -> list[str]:
    plane = RigidPlanePointContact([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 100.0)
    checked: list[str] = []

    assert plane.signed_gap_m([0.0, 0.25, 0.0]) == 0.25
    assert np.array_equal(plane.force_N([0.0, 0.25, 0.0]), np.zeros(3))
    assert plane.potential_J([0.0, 0.25, 0.0]) == 0.0
    checked.append("positive_gap_is_force_free_and_has_zero_potential")

    penetration = plane.force_N([0.0, -0.02, 0.0])
    assert np.allclose(penetration, [0.0, 2.0, 0.0])
    assert math.isclose(plane.potential_J([0.0, -0.02, 0.0]), 0.02, rel_tol=0.0, abs_tol=1e-15)
    checked.append("penetration_force_and_energy_match_linear_penalty")

    angle = 0.37
    c, s = math.cos(angle), math.sin(angle)
    # Rotate about z so the plane normal itself changes orientation.
    rotation = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    original_pos = np.array([0.1, -0.02, 0.3])
    rotated = RigidPlanePointContact(
        rotation @ plane.plane_point_m,
        rotation @ plane.unit_normal,
        plane.stiffness_N_per_m,
    )
    assert not np.allclose(rotated.unit_normal, plane.unit_normal)
    original_gap = plane.signed_gap_m(original_pos)
    rotated_gap = rotated.signed_gap_m(rotation @ original_pos)
    original_force = plane.force_N(original_pos)
    rotated_force = rotated.force_N(rotation @ original_pos)
    rotated_expected_force = rotation @ original_force
    original_potential = plane.potential_J(original_pos)
    rotated_potential = rotated.potential_J(rotation @ original_pos)
    assert math.isclose(rotated_gap, original_gap, abs_tol=1e-14)
    assert np.allclose(rotated_force, rotated_expected_force, atol=1e-13)
    assert math.isclose(rotated_potential, original_potential, abs_tol=1e-14)
    FOCUSED_CHECK_DETAILS["rotation_covariance"] = {
        "rotation_axis": "z",
        "angle_rad": angle,
        "original_unit_normal": plane.unit_normal.tolist(),
        "rotated_unit_normal": rotated.unit_normal.tolist(),
        "original_gap_m": original_gap,
        "rotated_gap_m": rotated_gap,
        "gap_abs_error_m": abs(rotated_gap - original_gap),
        "original_force_N": original_force.tolist(),
        "rotated_force_N": rotated_force.tolist(),
        "rotated_expected_force_N": rotated_expected_force.tolist(),
        "force_covariance_error_N": float(np.linalg.norm(rotated_force - rotated_expected_force)),
        "original_potential_J": original_potential,
        "rotated_potential_J": rotated_potential,
        "potential_abs_error_J": abs(rotated_potential - original_potential),
        "orientation_actually_changed": not bool(np.allclose(rotated.unit_normal, plane.unit_normal)),
    }
    checked.append("rigid_rotation_covariance_of_gap_force_and_energy")

    unlock_attempts = []
    for field_name in ("plane_point_m", "unit_normal"):
        array = getattr(plane, field_name)
        try:
            array.setflags(write=True)
        except ValueError as exc:
            unlock_attempts.append({
                "field": field_name,
                "setflags_write_true": "rejected",
                "exception_type": type(exc).__name__,
                "still_writeable": bool(array.flags.writeable),
            })
        else:
            array.setflags(write=False)
            raise AssertionError(f"{field_name} unexpectedly allowed setflags(write=True)")
    assert len(unlock_attempts) == 2 and all(not item["still_writeable"] for item in unlock_attempts)
    FOCUSED_CHECK_DETAILS["immutable_plane_arrays"] = unlock_attempts
    checked.append("plane_point_and_normal_cannot_be_unlocked")

    _expect_rejected("zero_normal_rejected", lambda: RigidPlanePointContact([0, 0, 0], [0, 0, 0], 100))
    _expect_rejected("nonunit_normal_rejected", lambda: RigidPlanePointContact([0, 0, 0], [0, 2, 0], 100))
    _expect_rejected("wrong_shape_normal_rejected", lambda: RigidPlanePointContact([0, 0, 0], [0, 1], 100))
    _expect_rejected("nonfinite_normal_rejected", lambda: RigidPlanePointContact([0, 0, 0], [0, float("nan"), 0], 100))
    checked.append("invalid_normals_rejected")

    _expect_rejected("zero_stiffness_rejected", lambda: RigidPlanePointContact([0, 0, 0], [0, 1, 0], 0))
    _expect_rejected("unit_string_rejected", lambda: RigidPlanePointContact([0, 0, 0], [0, 1, 0], "100 N/m"))
    _expect_rejected("nonfinite_plane_point_rejected", lambda: RigidPlanePointContact([0, float("inf"), 0], [0, 1, 0], 100))
    _expect_rejected("nonfinite_query_position_rejected", lambda: plane.force_N([0, float("nan"), 0]))
    checked.append("invalid_si_scalar_and_nonfinite_inputs_rejected")

    FOCUSED_CHECK_DETAILS["passed"] = checked.copy()
    return checked


def _energy_rows(mass_kg: float, velocity_m_per_s: float, y_m: float, force_N: float) -> tuple[float, float, float]:
    kinetic = 0.5 * mass_kg * velocity_m_per_s * velocity_m_per_s
    penetration = min(y_m, 0.0)
    potential = 0.5 * K_N_per_m * penetration * penetration
    return kinetic, potential, kinetic + potential


MASS_KG = 1.0
K_N_per_m = 100.0
INITIAL_GAP_M = 0.2
INITIAL_VELOCITY_M_PER_S = -2.0
END_TIME_S = 0.48


def simulate(case_id: str, dt_s: float) -> tuple[list[dict[str, float | int | str]], dict[str, float | str | int]]:
    plane = RigidPlanePointContact([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], K_N_per_m)
    nsteps_float = END_TIME_S / dt_s
    nsteps = int(round(nsteps_float))
    if not math.isclose(nsteps * dt_s, END_TIME_S, rel_tol=0.0, abs_tol=1e-13):
        raise ValueError("END_TIME_S must be an integer multiple of dt_s")

    y = INITIAL_GAP_M
    v = INITIAL_VELOCITY_M_PER_S
    time_s = 0.0
    force = float(plane.force_N([0.0, y, 0.0])[1])
    impulse_Ns = 0.0
    rows: list[dict[str, float | int | str]] = []

    def append_row(step: int) -> None:
        gap = plane.signed_gap_m([0.0, y, 0.0])
        kinetic, potential, total = _energy_rows(MASS_KG, v, y, force)
        rows.append({
            "case_id": case_id,
            "step": step,
            "time_s": time_s,
            "y_m": y,
            "signed_gap_m": gap,
            "vy_m_per_s": v,
            "contact_force_y_N": force,
            "contact_potential_J": potential,
            "kinetic_J": kinetic,
            "total_energy_J": total,
            "cumulative_contact_impulse_Ns": impulse_Ns,
        })

    append_row(0)
    for step in range(1, nsteps + 1):
        y_new = y + v * dt_s + 0.5 * (force / MASS_KG) * dt_s * dt_s
        time_new = step * dt_s
        force_new = float(plane.force_N([0.0, y_new, 0.0])[1])
        v_new = v + 0.5 * (force + force_new) / MASS_KG * dt_s
        impulse_Ns += 0.5 * (force + force_new) * dt_s
        y, v, force, time_s = y_new, v_new, force_new, time_new
        append_row(step)

    initial_energy = 0.5 * MASS_KG * INITIAL_VELOCITY_M_PER_S**2
    max_energy_error_J = max(abs(float(r["total_energy_J"]) - initial_energy) for r in rows)
    final_energy_error_J = float(rows[-1]["total_energy_J"]) - initial_energy
    point_delta_momentum_Ns = MASS_KG * (v - INITIAL_VELOCITY_M_PER_S)
    impulse_balance_error_Ns = point_delta_momentum_Ns - impulse_Ns

    hit_t = INITIAL_GAP_M / abs(INITIAL_VELOCITY_M_PER_S)
    contact_duration = math.pi * math.sqrt(MASS_KG / K_N_per_m)
    exact_exit_t = hit_t + contact_duration
    exact_exit_v = abs(INITIAL_VELOCITY_M_PER_S)
    exact_impulse_Ns = 2.0 * MASS_KG * abs(INITIAL_VELOCITY_M_PER_S)
    exit_time_s, exit_velocity_m_per_s = _detect_exit(rows)
    metrics: dict[str, float | str | int] = {
        "case_id": case_id,
        "dt_s": dt_s,
        "step_count": nsteps,
        "exact_first_contact_time_s": hit_t,
        "exact_contact_duration_s": contact_duration,
        "exact_rebound_time_s": exact_exit_t,
        "detected_rebound_time_s": exit_time_s,
        "rebound_time_abs_error_s": abs(exit_time_s - exact_exit_t),
        "exact_rebound_velocity_m_per_s": exact_exit_v,
        "detected_rebound_velocity_m_per_s": exit_velocity_m_per_s,
        "rebound_velocity_abs_error_m_per_s": abs(exit_velocity_m_per_s - exact_exit_v),
        "exact_elastic_contact_impulse_Ns": exact_impulse_Ns,
        "numerical_trapezoid_contact_impulse_Ns": impulse_Ns,
        "contact_impulse_abs_error_vs_analytic_Ns": abs(impulse_Ns - exact_impulse_Ns),
        "point_delta_momentum_Ns": point_delta_momentum_Ns,
        "impulse_balance_error_Ns": impulse_balance_error_Ns,
        "initial_total_energy_J": initial_energy,
        "final_total_energy_J": float(rows[-1]["total_energy_J"]),
        "final_total_energy_signed_error_J": final_energy_error_J,
        "max_abs_total_energy_error_J": max_energy_error_J,
        "max_abs_total_energy_relative_error": max_energy_error_J / initial_energy,
        "final_gap_m": float(rows[-1]["signed_gap_m"]),
        "final_velocity_m_per_s": v,
        "diagnostic_note": "energy errors are measured Verlet residuals, not a conservation pass/fail claim",
    }
    return rows, metrics


def _detect_exit(rows: list[dict[str, float | int | str]]) -> tuple[float, float]:
    previous = rows[0]
    for current in rows[1:]:
        g0, g1 = float(previous["signed_gap_m"]), float(current["signed_gap_m"])
        if g0 < 0.0 <= g1:
            fraction = -g0 / (g1 - g0)
            t0, t1 = float(previous["time_s"]), float(current["time_s"])
            v0, v1 = float(previous["vy_m_per_s"]), float(current["vy_m_per_s"])
            return t0 + fraction * (t1 - t0), v0 + fraction * (v1 - v0)
        previous = current
    raise AssertionError("no rebound exit found in the bounded trajectory")


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_plot(path: Path, cases: list[tuple[list[dict[str, Any]], dict[str, Any]]]) -> None:
    width, height = 1280, 760
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    draw.text((36, 18), "Fixed rigid plane contact: actual Velocity-Verlet samples", fill="black")

    def panel(rect: tuple[int, int, int, int], title: str, key: str, ylabel: str, ylow: float, yhigh: float) -> None:
        left, top, right, bottom = rect
        draw.rectangle(rect, outline=(80, 80, 80), width=1)
        draw.text((left + 8, top + 7), title, fill="black")
        plot_l, plot_t, plot_r, plot_b = left + 70, top + 42, right - 20, bottom - 36
        draw.line((plot_l, plot_b, plot_r, plot_b), fill="black", width=1)
        draw.line((plot_l, plot_t, plot_l, plot_b), fill="black", width=1)
        for fraction in (0.0, 0.5, 1.0):
            yy = plot_b - fraction * (plot_b - plot_t)
            draw.line((plot_l, yy, plot_r, yy), fill=(225, 225, 225), width=1)
            label_value = ylow + fraction * (yhigh - ylow)
            draw.text((left + 4, int(yy - 7)), f"{label_value:.4g}", fill=(60, 60, 60))
        draw.text((plot_l, bottom - 23), "time (s)", fill=(30, 30, 30))
        draw.text((left + 5, top + 23), ylabel, fill=(30, 30, 30))
        colors = [(25, 105, 190), (225, 95, 30)]
        for index, (rows, metrics) in enumerate(cases):
            pts = []
            for row in rows:
                x = plot_l + float(row["time_s"]) / END_TIME_S * (plot_r - plot_l)
                value = float(row[key])
                y = plot_b - (value - ylow) / (yhigh - ylow) * (plot_b - plot_t)
                y = min(max(y, plot_t), plot_b)
                pts.append((int(x), int(y)))
            draw.line(pts, fill=colors[index], width=2)
            label = f"dt={float(metrics['dt_s']):g}s"
            draw.line((plot_r - 185, top + 14 + 19 * index, plot_r - 160, top + 14 + 19 * index), fill=colors[index], width=3)
            draw.text((plot_r - 153, top + 7 + 19 * index), label, fill=colors[index])

    min_y = min(float(row["y_m"]) for rows, _ in cases for row in rows)
    max_y = max(float(row["y_m"]) for rows, _ in cases for row in rows)
    panel((24, 58, 1256, 390), "Point trajectory normal coordinate", "y_m", "y (m)", min_y - 0.02, max_y + 0.02)
    initial_energy = 0.5 * MASS_KG * INITIAL_VELOCITY_M_PER_S**2
    errors = [abs(float(row["total_energy_J"]) - initial_energy) for rows, _ in cases for row in rows]
    ycap = max(max(errors) * 1.12, 1e-12)
    panel((24, 410, 1256, 742), "Total energy error (measured residual; zero is reference)", "total_energy_error_J", "|E-E0| (J)", 0.0, ycap)
    # The second panel is the absolute energy error; keep the source trajectory rows unchanged.
    image.save(path, format="PNG")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    checks = run_self_checks()
    # Deliberately only two levels and 720 integration steps in total.
    cases = [("dt_0p002", 0.002), ("dt_0p001", 0.001)]
    result_rows: list[tuple[list[dict[str, Any]], dict[str, Any]]] = []
    all_rows: list[dict[str, Any]] = []
    metrics: list[dict[str, Any]] = []
    for case_id, dt_s in cases:
        rows, summary = simulate(case_id, dt_s)
        result_rows.append((rows, summary))
        all_rows.extend(rows)
        metrics.append(summary)

    coarse, fine = metrics
    convergence = {
        "energy_error_metric": "max absolute |(Kinetic + penalty Potential) - initial total energy| over saved Verlet samples",
        "coarse_dt_s": coarse["dt_s"],
        "fine_dt_s": fine["dt_s"],
        "coarse_max_abs_energy_error_J": coarse["max_abs_total_energy_error_J"],
        "fine_max_abs_energy_error_J": fine["max_abs_total_energy_error_J"],
        "fine_over_coarse_energy_error_ratio": float(fine["max_abs_total_energy_error_J"]) / float(coarse["max_abs_total_energy_error_J"]),
        "momentum_impulse_balance_abs_errors_Ns": [abs(float(item["impulse_balance_error_Ns"])) for item in metrics],
        "step_total": sum(int(item["step_count"]) for item in metrics),
        "interpretation": "reporting numerical residual and refinement ratio only; no energy-conservation pass/fail declaration",
    }
    exact = {
        "force_and_energy_basis": "U(g)=0.5*K*min(g,0)^2 and F=-dU/dx=-K*min(g,0)*n, for a fixed plane and unit normal n",
        "bounce_basis": "free flight to g=0 followed by m*y_ddot=-K*y while y<0; contact duration pi*sqrt(m/K), ideal exit speed equals entry speed for the elastic model",
        "parameters_si": {"mass_kg": MASS_KG, "stiffness_N_per_m": K_N_per_m, "initial_gap_m": INITIAL_GAP_M, "initial_velocity_m_per_s": INITIAL_VELOCITY_M_PER_S, "end_time_s": END_TIME_S},
        "dt_levels_s": [0.002, 0.001],
        "dt_ratio": 0.5,
        "total_steps_cap": 2000,
        "no_damping_or_friction": True,
    }
    payload = {
        "schema": "rigid_plane_point_contact_toy.v1",
        "status": "completed_bounded_diagnostic",
        "checks": {"count": len(checks), "passed": checks},
        "analytical_reference": exact,
        "cases": metrics,
        "refinement_and_balance": convergence,
        "limitations": [
            "This is a fixed infinite plane against one point mass with a scalar spring penalty.",
            "It is an independent force and potential API, not a vessel/fragment collision engine or a production integration.",
            "It does not implement damping, friction, moving planes, body rotation, self-contact, broad phase, or multiple contacts.",
            "Runtime scalar/vector validation enforces bare finite SI-number inputs and a unit normal; it cannot infer whether a caller mislabeled an ordinary bare number's physical unit.",
            "Energy errors are measured discrete Velocity-Verlet residuals; the error trend is reported and not promoted to a conservation pass/fail claim.",
        ],
    }
    write_csv(outdir / "actual_trajectory_force_energy_impulse.csv", all_rows)
    error_rows = []
    for item in metrics:
        error_rows.append({"case_id": item["case_id"], "dt_s": item["dt_s"], "metric": "detected_rebound_time_s", "measured": item["detected_rebound_time_s"], "reference": item["exact_rebound_time_s"], "absolute_error": item["rebound_time_abs_error_s"]})
        error_rows.append({"case_id": item["case_id"], "dt_s": item["dt_s"], "metric": "detected_rebound_velocity_m_per_s", "measured": item["detected_rebound_velocity_m_per_s"], "reference": item["exact_rebound_velocity_m_per_s"], "absolute_error": item["rebound_velocity_abs_error_m_per_s"]})
        error_rows.append({"case_id": item["case_id"], "dt_s": item["dt_s"], "metric": "max_abs_total_energy_error_J", "measured": item["max_abs_total_energy_error_J"], "reference": 0.0, "absolute_error": item["max_abs_total_energy_error_J"]})
        error_rows.append({"case_id": item["case_id"], "dt_s": item["dt_s"], "metric": "impulse_balance_error_Ns", "measured": item["impulse_balance_error_Ns"], "reference": 0.0, "absolute_error": abs(float(item["impulse_balance_error_Ns"]))})
    with (outdir / "comparison_errors.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(error_rows[0]))
        writer.writeheader()
        writer.writerows(error_rows)

    # Derive plotting field from the saved, actual total energy samples.
    for rows, _ in result_rows:
        for row in rows:
            row["total_energy_error_J"] = abs(float(row["total_energy_J"]) - 0.5 * MASS_KG * INITIAL_VELOCITY_M_PER_S**2)
    write_plot(outdir / "contact_trajectory_energy.png", result_rows)
    (outdir / "diagnostic_summary.json").write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    (outdir / "check_receipt.json").write_text(json.dumps({"status": "PASS", "checks": checks, "note": "Only the listed focused API checks and one bounded toy run were executed."}, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"outdir": str(outdir), "step_total": convergence["step_total"], "checks": len(checks), "refinement_and_balance": convergence, "cases": metrics}, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
