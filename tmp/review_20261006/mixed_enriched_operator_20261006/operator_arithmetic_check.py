"""Prepared tiny independent checks; do not run before a separate root GO."""
import importlib.util
import json
import math
from pathlib import Path
import struct
import sys
from dataclasses import replace

ROOT = Path(__file__).resolve().parents[4]
def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None: raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec); sys.modules[name] = module; spec.loader.exec_module(module)
    return module

HELP = load("mixed_lefm_helper_review", ROOT / "pgfrag/experiments/mixed_lefm_enrichment_algebra_20261006/mixed_enrichment_algebra.py")
OP = load("mixed_enriched_operator_review", ROOT / "pgfrag/experiments/mixed_enriched_operator_20261006/mixed_enriched_operator.py")
IDS = ("a", "b", "c", "d")
POS = {"a": (0.0, 0.0, 0.0), "b": (1.0, 0.0, 0.0), "c": (0.0, 1.0, 0.0), "d": (0.0, 0.0, 1.0)}
GRAD = {"a": (-1.0, -1.0, -1.0), "b": (1.0, 0.0, 0.0), "c": (0.0, 1.0, 0.0), "d": (0.0, 0.0, 1.0)}
I = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
R = 0.7
POINTS = ((0.21, 0.24, 0.17), (0.13, 0.32, 0.19))
WEIGHTS = (0.2, 0.3)
CHECKS = []

def check(ok, message):
    if not ok: raise AssertionError(message)

def close(a, b, scale=1.0, tol=3.0e-11):
    check(abs(a - b) <= tol * max(scale, abs(a), abs(b)), f"{a!r} != {b!r}")

def row_at(point, source_ids=IDS):
    x, y, z = point
    bary = {"a": 1.0-x-y-z, "b": x, "c": y, "d": z}
    raw = HELP.build_mixed_enriched_row(point_world_m=point, tip_world_m=(0.0, 0.0, 0.0),
        axes_world_from_local=I, reference_length_m=R, source_ids=source_ids,
        base_values=tuple(bary[s] for s in source_ids),
        base_gradients_world_per_m=tuple(GRAD[s] for s in source_ids))
    return OP.adapt_mixed_shape_row(raw, source_ids)

def branch_values_gradients(x, y):
    r = math.hypot(x, y); theta = math.atan2(y, x)
    sh, ch, st, ct = math.sin(theta/2), math.cos(theta/2), math.sin(theta), math.cos(theta)
    h = (sh, ch, sh*st, ch*st)
    hp = (0.5*ch, -0.5*sh, 0.5*ch*st+sh*ct, -0.5*sh*st+ch*ct)
    q, s = math.sqrt(r/R), math.sqrt(r/R)/r
    values = tuple(q*v for v in h)
    gradients = tuple((s*(0.5*ct*h[i]-st*hp[i]), s*(0.5*st*h[i]+ct*hp[i]), 0.0) for i in range(4))
    return values, gradients

def direct_b(component, gradient):
    d = [[0.0]*3 for _ in range(3)]; d[component][:] = gradient
    e = [[0.5*(d[i][j]+d[j][i]) for j in range(3)] for i in range(3)]
    return (e[0][0],e[1][1],e[2][2],2*e[0][1],2*e[1][2],2*e[0][2])

def independent_c(young, poisson):
    mu = young/(2*(1+poisson)); lam = young*poisson/((1+poisson)*(1-2*poisson))
    c = [[0.0]*6 for _ in range(6)]
    for i in range(3):
        for j in range(3): c[i][j] = lam + (2*mu if i == j else 0.0)
    for i in range(3,6): c[i][i] = mu
    return tuple(tuple(row) for row in c)

def independent_stiffness(rows, weights, c):
    n = len(rows[0].contributions); out = [0.0]*(n*n)
    for row, weight in zip(rows, weights):
        bcols = [direct_b(("ux","uy","uz").index(item.component), item.shape_gradient_world_per_m) for item in row.contributions]
        for i, bi in enumerate(bcols):
            for j, bj in enumerate(bcols):
                out[i*n+j] += weight*sum(bi[a]*sum(c[a][b]*bj[b] for b in range(6)) for a in range(6))
    return tuple(out)

def check_nulls(matrix_bytes, modes):
    n, count = modes.dof_count, modes.mode_count
    k = struct.unpack(f"<{n*n}d", matrix_bytes); v = struct.unpack(f"<{count*n}d", modes.values_le_f64)
    ratios = []
    for m in range(count):
        x = v[m*n:(m+1)*n]
        bound = max(sum(abs(k[i*n+j]*x[j]) for j in range(n)) for i in range(n))
        residual = max(abs(sum(k[i*n+j]*x[j] for j in range(n))) for i in range(n))
        ratio = residual/max(bound, 1.0e-300)
        check(ratio <= 2.0e-11, f"mode {m} residual ratio {ratio:.4e} is not null")
        ratios.append(ratio)
    return tuple(ratios)

rows = tuple(row_at(p) for p in POINTS)
for row in rows:
    for item in row.contributions:
        expected = direct_b(("ux", "uy", "uz").index(item.component), item.shape_gradient_world_per_m)
        for a, b in zip(item.B_column_per_m, expected): close(a, b, 1.0/R)
CHECKS.append("independent_B_reconstruction")

point, row = POINTS[0], rows[0]
offset = (0.3, -0.2, 0.4)
affine = ((0.2, -0.1, 0.3), (-0.3, 0.4, 0.1), (0.15, 0.2, -0.25))
front = ((0.12, -0.08, 0.05), (0.07, 0.11, -0.09), (-0.04, 0.06, 0.13), (0.09, -0.03, 0.08))
base_coeff = {(s,c): offset[c]+sum(affine[c][j]*POS[s][j] for j in range(3)) for s in IDS for c in range(3)}
u, grad = [0.0]*3, [[0.0]*3 for _ in range(3)]
for item in row.contributions:
    c = ("ux", "uy", "uz").index(item.component)
    coef = base_coeff[(item.owner_id,c)] if item.scope == "SOURCE_BASE" else front[item.basis_index-1][c]
    u[c] += item.shape_value*coef
    for j in range(3): grad[c][j] += item.shape_gradient_world_per_m[j]*coef
fv, fg = branch_values_gradients(point[0], point[1])
for c in range(3):
    eu = offset[c]+sum(affine[c][j]*point[j] for j in range(3))+sum(front[b][c]*fv[b] for b in range(4))
    for j in range(3): close(grad[c][j], affine[c][j]+sum(front[b][c]*fg[b][j] for b in range(4)), 1.0/R)
    close(u[c], eu)
CHECKS.append("mixed_affine_plus_four_global_branch_representation")

ids27 = tuple(f"s{i:02d}" for i in range(27))
check(OP.scalar_basis_count(27) == 31 and OP.vector_dof_count(27) == 93, "93-DOF global layout count")
check(OP.scalar_basis_count(27, OP.PER_SOURCE_DIAGNOSTIC_LAYOUT) == 135
      and OP.vector_dof_count(27, OP.PER_SOURCE_DIAGNOSTIC_LAYOUT) == 405, "405 diagnostic layout count")
check(len(set(OP.semantic_coefficient_ids(ids27))) == 93
      and len(set(OP.semantic_coefficient_ids(ids27, OP.PER_SOURCE_DIAGNOSTIC_LAYOUT))) == 405,
      "semantic layouts must expose every unique ID")
CHECKS.append("semantic_layout_counts_and_unique_IDs")

C = OP.isotropic_constitutive_matrix_pa(8.0e10, 0.23)
C_REF = independent_c(8.0e10, 0.23)
for i in range(36): close(C[i//6][i%6], C_REF[i//6][i%6], 8.0e10)
vpoints = tuple(OP.VolumePoint(w, r) for w, r in zip(WEIGHTS, rows))
K = OP.assemble_raw_stiffness(IDS, (vpoints,), C)
perm = IDS[::-1]
prows = tuple(row_at(p, perm) for p in POINTS)
Kperm = OP.assemble_raw_stiffness(perm, (tuple(OP.VolumePoint(w,r) for w,r in zip(WEIGHTS,prows)),), C)
check(type(K.values_le_f64) is bytes and K.layout == OP.GLOBAL_FRONT_LAYOUT and K.dof_count == 24,
      "global raw K must be immutable and semantically sized")
check(K.coefficient_ids == Kperm.coefficient_ids and K.values_le_f64 == Kperm.values_le_f64,
      "paired source permutation changed global raw K")
K_RAW = struct.unpack(f"<{K.dof_count*K.dof_count}d", K.values_le_f64); K_REF = independent_stiffness(rows, WEIGHTS, C_REF)
K_SCALE = max(abs(value) for value in K_REF); K_ERROR = max(abs(a-b) for a,b in zip(K_RAW, K_REF)); K_REL_ERROR = K_ERROR/K_SCALE
K_ASYMMETRY = max(abs(K_RAW[i*24+j]-K_RAW[j*24+i]) for i in range(24) for j in range(24)); K_REL_ASYMMETRY = K_ASYMMETRY/K_SCALE
check(K_REL_ERROR <= 3.0e-11, f"independent K max error {K_ERROR:.4e}, relative {K_REL_ERROR:.4e}")
check(K_REL_ASYMMETRY <= 3.0e-11, f"raw K relative asymmetry {K_REL_ASYMMETRY:.4e}")
CHECKS.append("independent_24DOF_raw_K_and_unsymmetrized_asymmetry")
RIGID_RATIOS = check_nulls(K.values_le_f64, OP.rigid_body_mode_coefficients(IDS, tuple(POS[s] for s in IDS)))
CHECKS.append("raw_K_and_six_physical_rigid_nulls")

diag_rows = tuple(OP.adapt_mixed_shape_row(HELP.build_mixed_enriched_row(
    point_world_m=p, tip_world_m=(0.0,0.0,0.0), axes_world_from_local=I, reference_length_m=R,
    source_ids=IDS, base_values=tuple({"a":1-p[0]-p[1]-p[2],"b":p[0],"c":p[1],"d":p[2]}[s] for s in IDS),
    base_gradients_world_per_m=tuple(GRAD[s] for s in IDS)), IDS, OP.PER_SOURCE_DIAGNOSTIC_LAYOUT) for p in POINTS)
Kd = OP.assemble_raw_stiffness(IDS, (tuple(OP.VolumePoint(w,r) for w,r in zip(WEIGHTS,diag_rows)),), C,
                               OP.PER_SOURCE_DIAGNOSTIC_LAYOUT)
DIAG_RIGID_RATIOS = check_nulls(Kd.values_le_f64, OP.rigid_body_mode_coefficients(
    IDS, tuple(POS[s] for s in IDS), OP.PER_SOURCE_DIAGNOSTIC_LAYOUT))
BRANCH_RATIOS = check_nulls(Kd.values_le_f64, OP.branch_dependency_mode_coefficients(IDS, tuple(POS[s] for s in IDS)))
check(K.dof_count == 24 and Kd.dof_count == 60, "4-source matrix dimensions must be 24 and 60")
CHECKS.append("405_diagnostic_six_additional_branch_null_candidates")

force = OP.assemble_raw_force(IDS, ((OP.SurfaceTractionPoint(0.4, (3.0, -2.0, 5.0), row),),))
fvect = struct.unpack(f"<{force.dof_count}d", force.values_le_f64)
expected_force = tuple(0.4*item.shape_value*(3.0,-2.0,5.0)[("ux","uy","uz").index(item.component)]
                       for item in row.contributions)
check(force.units == "N" and fvect == expected_force, "raw traction force or units mismatch")
CHECKS.append("raw_surface_traction_force")

def reject(fn, label):
    try: fn()
    except (OP.MixedOperatorInputError, FloatingPointError): return
    raise AssertionError(f"{label} did not fail closed")

bad_ids = replace(row, coefficient_ids=("wrong",)+row.coefficient_ids[1:])
reject(lambda: OP.assemble_raw_stiffness(IDS, ((OP.VolumePoint(0.2,bad_ids),),), C), "ID mismatch")
reject(lambda: OP.assemble_raw_stiffness(IDS, ((OP.VolumePoint(float("nan"),row),),), C), "nonfinite weight")
ambiguous = HELP.build_mixed_enriched_row(local_point_from_tip_m=(-0.3,0.0,0.1), source_ids=IDS,
    base_values=(0.25,0.25,0.25,0.25), base_gradients_world_per_m=tuple(GRAD[s] for s in IDS),
    axes_world_from_local=I, reference_length_m=R)
reject(lambda: OP.adapt_mixed_shape_row(ambiguous, IDS), "ambiguous cut")
CHECKS.append("ambiguous_nonfinite_and_ID_mismatch_fail_closed")
print(json.dumps({"status":"PASS_ARITHMETIC_ONLY","checks":CHECKS,
 "metrics":{"K_max_abs_error":K_ERROR,"K_relative_max_error":K_REL_ERROR,"raw_K_max_abs_asymmetry":K_ASYMMETRY,
 "raw_K_relative_asymmetry":K_REL_ASYMMETRY,"K_scale_max_abs":K_SCALE,"physical_null_ratios":RIGID_RATIOS,
 "diagnostic_rigid_null_ratios":DIAG_RIGID_RATIOS,"diagnostic_branch_null_ratios":BRANCH_RATIOS,
 "matrix_dof_counts":[K.dof_count,Kd.dof_count]},
 "27_source_inventory":"unique IDs only; no 93-DOF matrix or spectrum","physical_eligibility":"NOT_CERTIFIED; 93-DOF spectrum not counted."},sort_keys=True,separators=(",",":")))
