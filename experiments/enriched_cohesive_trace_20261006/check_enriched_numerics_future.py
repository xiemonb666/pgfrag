"""Future micro-checker source only; do not run until root approves its source."""

from __future__ import annotations

import importlib.util
import json
import math
import sys
import time
from pathlib import Path
from unittest.mock import patch


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("_enriched_trace_future_check", HERE / "enriched_cohesive_trace.py")
if SPEC is None or SPEC.loader is None: raise ImportError(HERE)
trace = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = trace
SPEC.loader.exec_module(trace)
IDS = ("base.tx", "base.ty", "base.tz", "base.rx", "base.ry", "base.rz", *(f"globalfront.f{i:02}" for i in range(18)))
MATERIAL = trace.BilinearCohesiveMaterial(1.0e12, 1.0e7, 100.0)
SCALAR = trace._DERIVATIVE._SCALAR
VIRGIN = SCALAR.initial_state(MATERIAL)
AREA, NORMAL, POINT_M = 2.0, (1./3., 2./3., 2./3.), (0.2, -0.1, 0.3)
ROTATION_SCALE_M, A = 1.0, (0.5, -0.25, 1.0)
FD_STEP_M, FD_ATOL_N_PER_M, FD_RTOL = 1.0e-10, 50.0, 2.0e-5
MAX_SECONDS, CALLS = 15.0, 0


def dot(a, b): return sum(x*y for x,y in zip(a,b))


def point(qid="q0", *, bad_rigid=False, slip=False):
    plus, minus = [[0.0]*24 for _ in range(3)], [[0.0]*24 for _ in range(3)]
    for c in range(3): plus[c][c] = minus[c][c] = 1.0
    rigid = ((0.,-POINT_M[2],POINT_M[1]), (POINT_M[2],0.,-POINT_M[0]), (-POINT_M[1],POINT_M[0],0.))
    for mode in range(3):
        for c in range(3): plus[c][3+mode] = minus[c][3+mode] = rigid[mode][c]/ROTATION_SCALE_M
    for j, aj in enumerate(A):
        for c in range(3):
            plus[c][6+j], minus[c][6+j] = .5*NORMAL[c]*aj, -.5*NORMAL[c]*aj
    if slip:
        slip_direction = (0., 2./3., -2./3.)
        for c in range(3): plus[c][9], minus[c][9] = .5*slip_direction[c], -.5*slip_direction[c]
    expected = ((1.,0.,0.), (0.,1.,0.), (0.,0.,1.), *rigid)
    if bad_rigid: expected = (*expected[:3], (expected[3][0]+.01, expected[3][1], expected[3][2]), *expected[4:])
    return trace.EnrichedFacePoint("toy-front", qid, trace.TraceFamily.GLOBALFRONT, IDS, plus, minus, NORMAL, POINT_M, POINT_M, expected, AREA)


P = point()
CANDIDATES = trace.RigidTraceCandidates(IDS, tuple(tuple(float(i == j)*ROTATION_SCALE_M for i in range(24)) for j in range(6)))


def coefficients(opening):
    q = [0.0]*24
    q[6] = opening/A[0]
    return tuple(q)


def state_map(points, state): return {p.point_id: state for p in points}


def invoke(q, *, points=(P,), states=None, candidates=CANDIDATES):
    global CALLS
    CALLS += 1
    history = state_map(points,VIRGIN) if states is None else states
    return trace.assemble_enriched_trial(MATERIAL,points,q,history,candidates)


def require(ok, message):
    if not ok: raise AssertionError(message)


def rejects(error, operation):
    try: operation()
    except error: return
    raise AssertionError(f"expected {error.__name__}")


def close(a,b,atol,rtol): return abs(a-b) <= atol+rtol*max(abs(a),abs(b))


def directional_fd(opening,committed,label):
    base = invoke(coefficients(opening),states=state_map((P,),committed))
    fd = {}
    for j, aj in zip(range(6,9), A):
        qplus = list(coefficients(opening))
        qplus[j] += FD_STEP_M
        plus = invoke(tuple(qplus),states=state_map((P,),committed))
        finite = tuple((x-y)/FD_STEP_M for x,y in zip(plus.generalized_internal_force_n,base.generalized_internal_force_n))
        analytic = tuple(base.tangent_n_per_m[i][j] for i in range(24))
        expected = tuple(AREA*base.point_results[0].k_alg_pa_per_m*A[i-6]*aj if 6 <= i <= 8 else 0.0 for i in range(24))
        error = max(abs(x-y) for x,y in zip(finite,analytic))
        reference = max(abs(x) for x in analytic)
        limit = FD_ATOL_N_PER_M+FD_RTOL*reference
        require(error <= limit and all(close(x,y,FD_ATOL_N_PER_M,FD_RTOL) for x,y in zip(analytic,expected)), f"{label}, column {j}: tangent mismatch")
        fd[str(j)] = {"finite_difference_force_column_n_per_m":finite,"assembled_tangent_column_n_per_m":analytic,"independent_area_k_a_outer_a_column_n_per_m":expected,"max_abs_error_n_per_m":error,"reference_max_abs_n_per_m":reference,"predeclared_limit_n_per_m":limit}
    require(base.point_results[0].tangent_status is trace.TangentStatus.DIFFERENTIABLE, f"{label}: kink")
    return base, fd


def baseline_record(result):
    p = result.point_results[0]
    aidx = range(6,9)
    action_reaction = max(abs(x+y) for x,y in zip(p.positive_face_force_at_reference_n,p.negative_face_force_at_reference_n))
    return {"opening_m":p.opening_m,"traction_pa":p.traction_pa,"k_alg_pa_per_m":p.k_alg_pa_per_m,"generalized_force_n":result.generalized_internal_force_n,"expected_area_t_a_n":tuple(AREA*p.traction_pa*x for x in A),"enriched_tangent_3x3_n_per_m":tuple(tuple(result.tangent_n_per_m[i][j] for j in aidx) for i in aidx),"expected_area_k_a_outer_a_n_per_m":tuple(tuple(AREA*p.k_alg_pa_per_m*x*y for y in A) for x in A),"positive_face_force_at_reference_n":p.positive_face_force_at_reference_n,"negative_face_force_at_reference_n":p.negative_face_force_at_reference_n,"action_reaction_max_residual_n":action_reaction,"stored_energy_j":result.stored_energy_j,"dissipated_energy_j":result.dissipated_energy_j,"external_work_j":result.external_work_j}


def main():
    started = time.perf_counter()
    d0, df = MATERIAL.delta0_m, MATERIAL.deltaf_m
    elastic, fd_elastic = directional_fd(.40*d0,VIRGIN,"elastic")
    soft, fd_soft = directional_fd(.75*df,VIRGIN,"softening")
    peak = SCALAR.update_state(MATERIAL,VIRGIN,.75*df)
    unload, fd_unload = directional_fd(.35*df,peak,"unload")
    for result in (elastic,soft,unload):
        p=result.point_results[0]
        expected=tuple(AREA*p.traction_pa*x for x in A)
        require(all(close(result.generalized_internal_force_n[6+i],expected[i],1.e-10,1.e-12) for i in range(3)), "independent A*t*a force mismatch")
        require(all(abs(sum(f*q for f,q in zip(result.generalized_internal_force_n,c))) < 1.e-10 for c in CANDIDATES.candidate_vectors_m), "rigid virtual work must vanish")
        face=p.positive_face_force_at_reference_n
        require(all(x == -y for x,y in zip(face,p.negative_face_force_at_reference_n)), "paired reference tractions must be equal/opposite")
    require(soft.point_results[0].k_alg_pa_per_m < 0., "post-peak tangent must stay negative")
    repeated=invoke(coefficients(.75*df))
    require(repeated == soft and VIRGIN.update_count == 0, "repeat trial changed result or committed history")
    broken=invoke(coefficients(1.10*df))
    require(close(broken.dissipated_energy_j,AREA*MATERIAL.gc_j_per_m2,1.e-10,1.e-12), "single-side area must produce A*Gc")
    require(close(broken.external_work_j,broken.stored_energy_j+broken.dissipated_energy_j,1.e-10,1.e-12), "full-break energy ledger mismatch")
    wrong=trace.RigidTraceCandidates(IDS[::-1],CANDIDATES.candidate_vectors_m)
    rejects(ValueError,lambda:invoke(coefficients(.4*d0),candidates=wrong))
    bad=point("bad-rigid",bad_rigid=True)
    rejects(ValueError,lambda:invoke(coefficients(0.),points=(bad,)))
    rejects(trace.UnsupportedEnrichedModeError,lambda:invoke(coefficients(-1.e-6)))
    slip_point=point("slip",slip=True)
    qslip=[0.0]*24
    qslip[9]=1.e-6
    rejects(trace.UnsupportedEnrichedModeError,lambda:invoke(tuple(qslip),points=(slip_point,)))
    rejects(ValueError,lambda:invoke(coefficients(0.),points=(P,P)))
    rejects(ValueError,lambda:invoke(coefficients(0.),states={}))
    with patch.object(trace,"_allocate_tangent",side_effect=AssertionError("tangent allocated before 16-QP cap")):
        many=tuple(point(f"q{i}") for i in range(17))
        rejects(ValueError,lambda:invoke(coefficients(0.),points=many,states=state_map(many,VIRGIN)))
    require(CALLS == 21, f"expected exactly 21 assembly calls, got {CALLS}")
    elapsed=time.perf_counter()-started
    require(elapsed <= MAX_SECONDS,"future micro-check exceeded its declared 15-second budget")
    result={"status":"PASS","kind":"toy mathematical micro-check only; not physical validation","assembly_calls":CALLS,"fd_scheme":"forward finite differences at differentiable points","fd_step_m":FD_STEP_M,"fd_abs_tol_n_per_m":FD_ATOL_N_PER_M,"fd_rel_tol":FD_RTOL,"branches":{"elastic":{"baseline":baseline_record(elastic),"fd_columns":fd_elastic},"softening":{"baseline":baseline_record(soft),"fd_columns":fd_soft},"unload":{"baseline":baseline_record(unload),"fd_columns":fd_unload}},"fully_broken":{"stored_energy_j":broken.stored_energy_j,"dissipated_energy_j":broken.dissipated_energy_j,"expected_area_gc_j":AREA*MATERIAL.gc_j_per_m2,"external_work_j":broken.external_work_j},"rigid_virtual_work_checks":6,"semantic_order_rejection":"PASS","bad_rigid_reproduction_rejection":"PASS","mode_rejections":["compression","tangential slip"],"duplicate_point_id_rejection":"PASS","history_id_mismatch_rejection":"PASS","cap_rejection_qps":17,"budget_s":MAX_SECONDS,"elapsed_s":elapsed}
    out=HERE.parents[1]/"tmp/review_20261006/enriched_cohesive_trace_20261006/future_run01/numeric_result.json"
    if out.exists(): raise RuntimeError("refusing to overwrite future_run01 numeric result")
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2))


if __name__ == "__main__": main()
