"""One bounded source-only geometry check for the safe-v2 classifier."""
from __future__ import annotations

import copy
import json
import math
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from support_topology import (  # noqa: E402
    FiniteCrackGeometry,
    GeometryInputError,
    classify_support,
)

ROOT = HERE.parents[2]
CASE_PATH = ROOT / "pgfrag/experiments/local_front_support_contract_20261006/fixed_case.json"


def main() -> None:
    case = json.loads(CASE_PATH.read_text(encoding="utf-8"))
    tolerance = case["topology"]["tolerance_m"]
    h = case["support"]["radius_m"]
    box_z_max = case["geometry"]["domain_box_m"][2][1]
    calls = 0

    def classify(name: str, point: tuple[float, float, float], geometry: FiniteCrackGeometry):
        nonlocal calls
        calls += 1
        return classify_support(name, point, geometry)

    # Original H counterexample: a 0.5-tolerance protrusion must fail closed.
    for name, overrun in (("below_tol_overhang", 0.5 * tolerance),
                          ("above_tol_overhang", 1.5 * tolerance)):
        case_with_end = copy.deepcopy(case)
        case_with_end["geometry"]["crack"]["z_bounds_m"][1] = box_z_max - overrun
        cert = classify(name, (-0.008, 0.004, 0.004), FiniteCrackGeometry.from_case(case_with_end))
        assert cert["status"] == "UNKNOWN", (name, cert["reason_code"])
        assert cert["geometric_role"] == "UNKNOWN", (name, cert["geometric_role"])

    # A nonzero interior margin smaller than tolerance is also uncertain.
    geometry = FiniteCrackGeometry.from_case(case)
    rho = math.sqrt(h * h - 0.004 * 0.004)
    target_z_hi = box_z_max - 0.5 * tolerance
    near_inside_source = (-0.008, 0.004, target_z_hi - rho)
    near_inside = classify("below_tol_inside_margin", near_inside_source, geometry)
    assert near_inside["status"] == "UNKNOWN", near_inside["reason_code"]

    # The identical section with a margin comfortably above tolerance is usable.
    target_z_hi = box_z_max - 2.0 * tolerance
    safe_inside_source = (-0.008, 0.004, target_z_hi - rho)
    safe_inside = classify("above_tol_inside_margin", safe_inside_source, geometry)
    assert safe_inside["status"] == "CERTIFIED", safe_inside["reason_code"]
    assert safe_inside["geometric_role"] == "HEAVISIDE_ELIGIBLE"
    assert safe_inside["basis_role"] == "HEAVISIDE"

    # Exact box-clipped z endpoint remains legal and certifies H.
    exact_endpoint = classify("exact_box_crack_endpoint", (-0.008, 0.004, 0.008), geometry)
    assert exact_endpoint["status"] == "CERTIFIED", exact_endpoint["reason_code"]
    assert exact_endpoint["geometric_role"] == "HEAVISIDE_ELIGIBLE"
    assert exact_endpoint["plane_section_bounds_m"]["z_max_m"] == 0.008
    assert exact_endpoint["plane_section_inside_finite_crack"] is True

    # Original front false-positive counterexample: reject a finite patch wholly
    # outside Omega even when it is less than one tolerance beyond the box.
    outside_front_case = copy.deepcopy(case)
    outside_front_case["geometry"]["crack"]["z_bounds_m"] = [
        box_z_max + 0.1 * tolerance,
        box_z_max + 0.9 * tolerance,
    ]
    try:
        FiniteCrackGeometry.from_case(outside_front_case)
    except GeometryInputError:
        outside_front_rejected = True
    else:
        outside_front_rejected = False
    assert outside_front_rejected

    # A source center beyond the closed physical box is UNKNOWN at any distance.
    out_of_domain = classify("source_outside_box", (0.004, 0.004, box_z_max + 0.5 * tolerance), geometry)
    assert out_of_domain["status"] == "UNKNOWN"
    assert out_of_domain["reason_code"] == "source_center_outside_domain"

    # Recompute the frozen 5 x 4 x 5 lattice without provider or manifest code.
    axes = case["sources"]["axes_m"]
    role_counts = {"FRONT": 0, "HEAVISIDE_ELIGIBLE": 0, "BASE": 0, "UNKNOWN": 0}
    dofs = 0
    source_y_nonzero = True
    for ix, x in enumerate(axes["x"]):
        for iy, y in enumerate(axes["y"]):
            for iz, z in enumerate(axes["z"]):
                source_id = f"frozen-{ix}-{iy}-{iz}"
                cert = classify(source_id, (x, y, z), geometry)
                role_counts[cert["geometric_role"]] += 1
                dofs += len(cert["basis_roles"])
                source_y_nonzero = source_y_nonzero and y != 0.0
                assert cert["status"] == "CERTIFIED", (source_id, cert["reason_code"])

    expected_counts = {"FRONT": 30, "HEAVISIDE_ELIGIBLE": 10, "BASE": 60, "UNKNOWN": 0}
    assert role_counts == expected_counts, role_counts
    assert dofs == 140, dofs
    assert source_y_nonzero
    assert calls == 106, calls

    print(json.dumps({
        "result": "PASS",
        "classification_queries": calls,
        "elapsed_limit_seconds": 10,
        "overhang_below_and_above_tolerance": "UNKNOWN",
        "nonzero_interior_gap_below_tolerance": "UNKNOWN",
        "interior_gap_above_tolerance": "HEAVISIDE",
        "exact_box_clipped_endpoint": "HEAVISIDE",
        "finite_front_outside_box_within_tolerance": "REJECTED",
        "source_outside_closed_box": "UNKNOWN",
        "frozen_role_counts": role_counts,
        "frozen_scalar_dofs": dofs,
        "all_source_y_nonzero": source_y_nonzero,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
