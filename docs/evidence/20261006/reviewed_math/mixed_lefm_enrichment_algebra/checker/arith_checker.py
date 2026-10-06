"""One-shot, bounded arithmetic checker for mixed_enrichment_algebra.py."""
import importlib.util
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[4]
SOURCE = ROOT / "pgfrag" / "experiments" / "mixed_lefm_enrichment_algebra_20261006" / "mixed_enrichment_algebra.py"
SPEC = importlib.util.spec_from_file_location("mixed_enrichment_under_review", SOURCE)
if SPEC is None or SPEC.loader is None: raise RuntimeError("cannot load pinned arithmetic helper")
api = importlib.util.module_from_spec(SPEC); sys.modules[SPEC.name] = api
SPEC.loader.exec_module(api)

I = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
R, H = 0.65, 1.0e-5
P, N0, G = (0.37, 0.29, -0.12), 0.63, (0.40, -0.25, 0.30)
COORDS, CALLS, CHECKS = set(), 0, []

def check(ok, message):
    if not ok:
        raise AssertionError(message)
def near(a, b, scale, rel, message):
    check(abs(a - b) <= rel * max(scale, abs(a), abs(b)),
          f"{message}: {a!r} vs {b!r}")
def reserve(point):
    global CALLS
    if CALLS >= 24:
        raise AssertionError("helper build-call budget exceeded")
    COORDS.add(tuple(float(v) for v in point))
    if len(COORDS) > 12:
        raise AssertionError("distinct-coordinate budget exceeded")
    CALLS += 1
def world(point, *, ids=("probe",), values=(N0,), grads=(G,), q=I, tip=(0.0, 0.0, 0.0)):
    reserve(point)
    return api.build_mixed_enriched_row(
        point_world_m=point, tip_world_m=tip, source_ids=ids,
        base_values=values, base_gradients_world_per_m=grads,
        axes_world_from_local=q, reference_length_m=R)
def local(point, *, trace=None, ids=("probe",), values=(N0,), grads=(G,), q=I):
    reserve(point)
    return api.build_mixed_enriched_row(
        local_point_from_tip_m=point, source_ids=ids, base_values=values,
        base_gradients_world_per_m=grads, axes_world_from_local=q,
        reference_length_m=R, branch_trace=trace)
def item(row, source, component, kind, index):
    hits = [x for x in row.contributions if (x.source_id, x.component, x.basis_kind, x.basis_index)
            == (source, component, kind, index)]
    check(len(hits) == 1, f"expected one contribution for {source}/{component}/{kind}/{index}")
    return hits[0]
def mv(a, x):
    return tuple(sum(a[i][j] * x[j] for j in range(3)) for i in range(3))
def direct_b(component, gradient):
    d = [[0.0, 0.0, 0.0] for _ in range(3)]; d[component][:] = gradient
    e = [[0.5 * (d[i][j] + d[j][i]) for j in range(3)] for i in range(3)]
    return (e[0][0], e[1][1], e[2][2], 2.0 * e[0][1], 2.0 * e[1][2], 2.0 * e[0][2])

center = world(P)
for axis in range(3):
    for sign in (-1.0, 1.0):
        point = tuple(P[j] + (sign * H if j == axis else 0.0) for j in range(3))
        varied_n = N0 + sign * H * G[axis]
        row = world(point, values=(varied_n,))
        if sign < 0.0:
            minus = row
        else:
            plus = row
    for branch in range(4):
        fd = (plus.branch_values[branch] - minus.branch_values[branch]) / (2.0 * H)
        near(fd, center.branch_gradients_world_per_m[branch][axis], 1.0 / R, 2.0e-8,
             f"branch {branch + 1} central derivative axis {axis}")
        fp = item(plus, "probe", "ux", "front_branch", branch + 1)
        fm = item(minus, "probe", "ux", "front_branch", branch + 1)
        fc = (fp.shape_value - fm.shape_value) / (2.0 * H)
        gc = item(center, "probe", "ux", "front_branch", branch + 1)
        near(fc, gc.shape_gradient_world_per_m[axis], 1.0 / R, 2.0e-8,
             f"product branch {branch + 1} central derivative axis {axis}")
CHECKS.append("off_cut_central_differences")

a, b = 0.37, -0.51
ca, sa, cb, sb = math.cos(a), math.sin(a), math.cos(b), math.sin(b)
Q = ((ca * cb, -sa, ca * sb), (sa * cb, ca, sa * sb), (-sb, 0.0, cb))
tip = (0.13, -0.27, 0.41)
rot_point = tuple(tip[i] + mv(Q, P)[i] for i in range(3))
rot = world(rot_point, q=Q, tip=tip, grads=(mv(Q, G),))
for j in range(4):
    near(rot.branch_values[j], center.branch_values[j], 1.0, 2.0e-11, f"rotated branch value {j + 1}")
    expected = mv(Q, center.branch_gradients_world_per_m[j])
    for k in range(3):
        near(rot.branch_gradients_world_per_m[j][k], expected[k], 1.0 / R, 2.0e-11,
             f"rotated branch gradient {j + 1}/{k}")
    for comp in ("ux", "uy", "uz"):
        g0 = item(center, "probe", comp, "front_branch", j + 1).shape_gradient_world_per_m
        gr = item(rot, "probe", comp, "front_branch", j + 1).shape_gradient_world_per_m
        er = mv(Q, g0)
        for k in range(3):
            near(gr[k], er[k], 1.0 / R, 2.0e-11, f"rotated product gradient {j + 1}/{comp}/{k}")
CHECKS.append("proper_frame_rotation_covariance")
four_r = world((4.0 * P[0], 4.0 * P[1], P[2]))
for j in range(4):
    near(four_r.branch_values[j], 2.0 * center.branch_values[j], 1.0, 2.0e-12,
         f"radial value homogeneity {j + 1}")
    for k in range(3):
        near(four_r.branch_gradients_world_per_m[j][k],
             0.5 * center.branch_gradients_world_per_m[j][k], 1.0 / R, 2.0e-12,
             f"radial gradient homogeneity {j + 1}/{k}")
CHECKS.append("radial_homogeneity")

for contribution in center.contributions:
    expected = direct_b(("ux", "uy", "uz").index(contribution.component),
                        contribution.shape_gradient_world_per_m)
    for got, want in zip(contribution.B_column_per_m, expected):
        near(got, want, 1.0 / R, 2.0e-12, "direct symmetric-gradient B column")
CHECKS.append("direct_symmetric_gradient_B_columns")
ids, vals, grads = ("zeta", "alpha"), (0.8, -0.15), ((0.2, 0.3, -0.1), (-0.4, 0.1, 0.25))
pair1 = world(P, ids=ids, values=vals, grads=grads)
pair2 = world(P, ids=ids[::-1], values=vals[::-1], grads=grads[::-1])
check(pair1.coefficient_ids == pair2.coefficient_ids and pair1.contributions == pair2.contributions,
      "paired source permutation changed sorted semantic rows")
many_ids = tuple(f"synthetic-{i:02d}" for i in range(27))
many_values = tuple(N0 + i * 0.01 for i in range(27))
many_grads = tuple((G[0] + i * 0.001, G[1], G[2]) for i in range(27))
many = world(P, ids=many_ids, values=many_values, grads=many_grads)
expected_ids = {json.dumps([sid, comp, kind, index], ensure_ascii=True, separators=(",", ":"))
                for sid in many_ids for kind, index in (("base", 0), *(('front_branch', k) for k in range(1, 5)))
                for comp in ("ux", "uy", "uz")}
check(len(many.coefficient_ids) == 405 and len(set(many.coefficient_ids)) == 405
      and set(many.coefficient_ids) == expected_ids, "27-source semantic inventory is not exact/unique")
check(api.mixed_dof_count(27) == 405, "27-source DOF count changed")
CHECKS.append("paired_permutation_and_27_source_semantic_inventory")

front = local((0.0, 0.0, 0.0))
check(front.status is api.RowStatus.UNKNOWN_FRONT and not front.coefficient_ids and not front.contributions,
      "front must be UNKNOWN_FRONT with no partial contributions")
cut_point = (-0.31, 0.0, 0.17)
ambiguous = local(cut_point)
upper, lower = local(cut_point, trace="upper"), local(cut_point, trace="lower")
check(ambiguous.status is api.RowStatus.UNKNOWN_BRANCH_FACE and not ambiguous.coefficient_ids
      and not ambiguous.contributions, "unselected cut trace must fail closed")
check(upper.status is api.RowStatus.OK and lower.status is api.RowStatus.OK
      and upper.theta_rad == math.pi and lower.theta_rad == -math.pi
      and len(upper.coefficient_ids) == 15 and len(lower.coefficient_ids) == 15,
      "explicit upper/lower face traces must be complete and distinct")
near(upper.branch_values[0], -lower.branch_values[0], 1.0, 2.0e-12, "upper/lower first-branch trace")
CHECKS.append("front_cut_and_explicit_traces")

def must_reject(fn, label):
    try:
        fn()
    except (api.EnrichmentContractError, FloatingPointError):
        return
    raise AssertionError(f"{label} was not rejected")
must_reject(lambda: world(P, q=((-1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))),
            "improper crack frame")
must_reject(lambda: world(P, values=(float("nan"),)), "nonfinite base value")
CHECKS.append("improper_frame_and_nonfinite_rejection")
check(CALLS == 18 and len(COORDS) == 11, f"actual budget use unexpected: {CALLS} calls/{len(COORDS)} coordinates")
print(json.dumps({"status": "PASS", "checks": CHECKS, "distinct_query_coordinates": len(COORDS),
                  "helper_build_calls": CALLS, "scope": "arithmetic only; no MLS/provider/solver/mechanics"},
                 ensure_ascii=True, sort_keys=True, separators=(",", ":")))
