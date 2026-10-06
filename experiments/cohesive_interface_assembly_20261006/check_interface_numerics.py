"""Toy math checker for root authorization; not ceramic/physical validation; do not run yet."""

from __future__ import annotations

import hashlib, importlib.util, json, math, time
from pathlib import Path
import sys
from unittest.mock import patch


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "cohesive_interface_assembly.py"
SCALAR_PATH = HERE.parent / "cohesive_scalar_law_20261006" / "cohesive_scalar_law.py"
SCALAR_PIN = "78433638608AD4216979607412FE0D4966F1423706C47F734F50DD2A8121FA17"
ASSEMBLY_PIN = "d2bf0dd79507e5ab26efe69b5489c3284b7c9118656ab49fffe942303ba11854"
if hashlib.sha256(SCALAR_PATH.read_bytes()).hexdigest().upper() != SCALAR_PIN: raise RuntimeError("frozen scalar-law pin changed")
if hashlib.sha256(SOURCE.read_bytes()).hexdigest() != ASSEMBLY_PIN: raise RuntimeError("reviewed interface source pin changed")
spec = importlib.util.spec_from_file_location("_cohesive_interface_trial_check", SOURCE)
if spec is None or spec.loader is None:
    raise ImportError(SOURCE)
assembly = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = assembly
spec.loader.exec_module(assembly)
scalar = assembly._SCALAR
material = scalar.BilinearCohesiveMaterial(1.0e12, 1.0e7, 100.0)
area = 2.0
NODE_POSITIONS = ((0.,0.,0.), (1.,0.,0.), (0.,0.,0.), (1.,0.,0.))
FORCE_ATOL_N, FORCE_RTOL = 1.0e-12, 1.0e-10
MOMENT_ATOL_NM, MOMENT_RTOL = 1.0e-12, 1.0e-10
FD_STEP_M, FD_ATOL_N_PER_M, FD_RTOL = 1.0e-10, 50.0, 2.0e-5
ENERGY_ATOL_J, ENERGY_RTOL = 1.0e-10, 1.0e-12
MAX_SECONDS = 10.0
CALLS = 0


def point(qid="q0", node_count=4):
    ndof = 3*node_count
    plus, minus = [[0.]*ndof for _ in range(3)], [[0.]*ndof for _ in range(3)]
    if node_count == 4:
        plus_nodes, minus_nodes, q = ((2,.5),(3,.5)), ((0,.5),(1,.5)), (.5,0.,0.)
    else:
        plus_nodes, minus_nodes, q = ((node_count-1,1.),), ((0,1.),), (0.,0.,0.)
    for c in range(3):
        for node, weight in plus_nodes: plus[c][3*node+c] = weight
        for node, weight in minus_nodes: minus[c][3*node+c] = weight
    return assembly.CohesiveQuadraturePoint("toy", qid, plus, minus, (0.,0.,1.), q, q, area)


P = point()
VIRGIN = scalar.initial_state(material)


def displacement(opening, slip=0.):
    d = [0.]*12
    for node in (2,3): d[3*node] = slip; d[3*node+2] = opening
    return tuple(d)


def states_for(points, state): return {p.point_id: state for p in points}


def invoke(d, *, points=(P,), positions=NODE_POSITIONS, states=None):
    global CALLS
    CALLS += 1
    history = states_for(points, VIRGIN) if states is None else states
    return assembly.assemble_trial(material, points, positions, d, history)


def close(a, b, atol, rtol): return abs(a-b) <= atol+rtol*max(abs(a),abs(b))


def require(ok, message):
    if not ok: raise AssertionError(message)


def rejects(error, operation):
    try: operation()
    except error: return
    raise AssertionError(f"expected {error.__name__}")


def vector_close(a, b, atol, rtol):
    require(len(a)==len(b) and all(close(x,y,atol,rtol) for x,y in zip(a,b)), "vector mismatch")


def matvec(matrix, vector): return tuple(sum(a*b for a,b in zip(row,vector)) for row in matrix)


def cross(a,b): return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def resultant_and_torque(force, positions):
    nodal = tuple(tuple(force[3*i+c] for c in range(3)) for i in range(len(positions)))
    resultant = tuple(sum(f[c] for f in nodal) for c in range(3))
    torque = tuple(sum(cross(x,f)[c] for x,f in zip(positions,nodal)) for c in range(3))
    return resultant, torque


def max_abs(values): return max((abs(x) for x in values), default=0.0)


def directional_fd(opening, committed, label):
    direction = tuple(1. if i in (8,11) else 0. for i in range(12))
    base = invoke(displacement(opening), states=states_for((P,), committed))
    plus = invoke(displacement(opening+FD_STEP_M), states=states_for((P,), committed))
    minus = invoke(displacement(opening-FD_STEP_M), states=states_for((P,), committed))
    finite = tuple((a-b)/(2*FD_STEP_M) for a,b in zip(plus.internal_force_n,minus.internal_force_n))
    analytic = matvec(base.tangent_n_per_m, direction)
    error, reference = max_abs(a-b for a,b in zip(finite,analytic)), max_abs(analytic)
    limit = FD_ATOL_N_PER_M+FD_RTOL*reference
    require(error <= limit, f"{label}: FD error {error} exceeds {limit}")
    require(base.point_results[0].tangent_status is assembly.TangentStatus.DIFFERENTIABLE, f"{label} must be away from kinks")
    return base, {"max_abs_error_n_per_m":error,"reference_max_abs_n_per_m":reference,"predeclared_limit_n_per_m":limit}


def branch_record(result):
    point_result=result.point_results[0]
    return {"opening_m":point_result.opening_m,"traction_pa":point_result.traction_pa,"k_alg_pa_per_m":point_result.k_alg_pa_per_m,"stored_energy_j":result.stored_energy_j,"dissipated_energy_j":result.dissipated_energy_j,"external_work_j":result.external_work_j}


def main():
    global CALLS
    started = time.perf_counter()
    counts = {"scalar_update":0, "scalar_traction":0}
    update, traction = assembly._scalar_update, assembly._scalar_traction
    assembly._scalar_update = lambda *a,**k: (counts.__setitem__("scalar_update",counts["scalar_update"]+1) or update(*a,**k))
    assembly._scalar_traction = lambda *a,**k: (counts.__setitem__("scalar_traction",counts["scalar_traction"]+1) or traction(*a,**k))
    try:
        d0, df = material.delta0_m, material.deltaf_m
        fd_reports = {}
        pre, fd_reports["elastic"] = directional_fd(.40*d0, VIRGIN, "elastic loading")
        soft, fd_reports["soft"] = directional_fd(.75*df, VIRGIN, "softening loading")
        require(soft.point_results[0].k_alg_pa_per_m < 0.0, "post-peak tangent must be explicitly negative")
        committed_peak = update(material, VIRGIN, .75*df)
        unload, fd_reports["unload"] = directional_fd(.35*df, committed_peak, "unloading")
        repeated = invoke(displacement(.75*df), states=states_for((P,), VIRGIN))
        require(repeated == soft and VIRGIN.update_count == 0, "repeat trial changed output or committed history")

        plus_resultant, plus_torque = resultant_and_torque(soft.point_results[0].positive_face_force_n, NODE_POSITIONS)
        minus_resultant, minus_torque = resultant_and_torque(soft.point_results[0].negative_face_force_n, NODE_POSITIONS)
        total_resultant, total_torque = resultant_and_torque(soft.internal_force_n, NODE_POSITIONS)
        action_reaction_residual = max_abs(a+b for a,b in zip(plus_resultant,minus_resultant))
        pair_torque_residual = max_abs(a+b for a,b in zip(plus_torque,minus_torque))
        vector_close(plus_resultant, soft.point_results[0].positive_face_resultant_n, FORCE_ATOL_N, FORCE_RTOL)
        vector_close(minus_resultant, soft.point_results[0].negative_face_resultant_n, FORCE_ATOL_N, FORCE_RTOL)
        vector_close(plus_resultant, tuple(-x for x in minus_resultant), FORCE_ATOL_N, FORCE_RTOL)
        vector_close(tuple(plus_torque[i]+minus_torque[i] for i in range(3)), (0.,0.,0.), MOMENT_ATOL_NM, MOMENT_RTOL)
        vector_close(total_resultant, (0.,0.,0.), FORCE_ATOL_N, FORCE_RTOL)
        vector_close(total_torque, (0.,0.,0.), MOMENT_ATOL_NM, MOMENT_RTOL)

        failed = invoke(displacement(1.10*df))
        vector_close((failed.dissipated_energy_j,), (area*material.gc_j_per_m2,), ENERGY_ATOL_J, ENERGY_RTOL)
        vector_close((failed.external_work_j,), (failed.stored_energy_j+failed.dissipated_energy_j,), ENERGY_ATOL_J, ENERGY_RTOL)
        failed_point=failed.point_results[0]
        failed_plus, failed_minus = resultant_and_torque(failed_point.positive_face_force_n,NODE_POSITIONS), resultant_and_torque(failed_point.negative_face_force_n,NODE_POSITIONS)
        failed_total = resultant_and_torque(failed.internal_force_n,NODE_POSITIONS)
        failed_action_reaction_residual=max_abs(a+b for a,b in zip(failed_plus[0],failed_minus[0]))
        failed_torque_residual=max(max_abs(a+b for a,b in zip(failed_plus[1],failed_minus[1])),max_abs(failed_total[1]))
        rejects(assembly.UnsupportedInterfaceModeError, lambda: invoke(displacement(-1.e-6)))
        rejects(assembly.UnsupportedInterfaceModeError, lambda: invoke(displacement(1.e-6, slip=1.e-6)))
        rejects(ValueError, lambda: invoke(displacement(0.), points=(P,P)))
        rejects(ValueError, lambda: invoke(displacement(0.), states={}))

        cap_allocator = patch.object(assembly, "_allocate_tangent", side_effect=AssertionError("allocated before cap rejection"))
        with cap_allocator:
            many = tuple(point(f"q{i}") for i in range(17))
            rejects(ValueError, lambda: invoke(displacement(0.), points=many, states=states_for(many,VIRGIN)))
            big_positions = tuple((0.,0.,0.) for _ in range(13))
            big_point = point("qbig", node_count=13)
            rejects(ValueError, lambda: invoke((0.,)*39, points=(big_point,), positions=big_positions, states=states_for((big_point,),VIRGIN)))
        require(CALLS == 17, f"expected 17 assembly calls, got {CALLS}")
        require(counts == {"scalar_update":11,"scalar_traction":11}, f"unexpected scalar calls: {counts}")
        elapsed = time.perf_counter()-started
        require(elapsed <= MAX_SECONDS, f"numeric checker exceeded {MAX_SECONDS}s budget")
        print(json.dumps({"status":"PASS","kind":"toy mathematical checks only; not vessel/physical validation","assembly_calls":CALLS,"scalar_updates_in_assembly":counts["scalar_update"],"scalar_traction_calls":counts["scalar_traction"],"direct_scalar_history_preparation_calls":1,"main_qp_count":1,"oversize_rejection_qps":17,"oversize_rejection_nodes":13,"fd_branches":fd_reports,"baseline":{"elastic":branch_record(pre),"soft":branch_record(soft),"unload":branch_record(unload)},"soft_action_reaction_max_residual_n":action_reaction_residual,"soft_paired_torque_max_residual_nm":pair_torque_residual,"fully_broken":{"opening_m":failed_point.opening_m,"stored_energy_j":failed.stored_energy_j,"dissipated_energy_j":failed.dissipated_energy_j,"external_work_j":failed.external_work_j,"expected_gc_energy_j":area*material.gc_j_per_m2,"energy_ledger_abs_residual_j":abs(failed.external_work_j-failed.stored_energy_j-failed.dissipated_energy_j),"action_reaction_max_residual_n":failed_action_reaction_residual,"torque_max_residual_nm":failed_torque_residual},"fd_step_m":FD_STEP_M,"fd_atol_n_per_m":FD_ATOL_N_PER_M,"fd_rtol":FD_RTOL,"force_abs_tol_n":FORCE_ATOL_N,"force_rel_tol":FORCE_RTOL,"moment_abs_tol_nm":MOMENT_ATOL_NM,"moment_rel_tol":MOMENT_RTOL,"energy_atol_j":ENERGY_ATOL_J,"energy_rtol":ENERGY_RTOL,"budget_s":MAX_SECONDS,"elapsed_s":elapsed},indent=2))
    finally:
        assembly._scalar_update, assembly._scalar_traction = update, traction


if __name__ == "__main__": main()
