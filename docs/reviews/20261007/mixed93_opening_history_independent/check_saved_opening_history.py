#!/usr/bin/env python3
"""Independent stdlib audit of one saved 93-DOF Mode-I opening history."""
import argparse, base64, hashlib, json, math, struct, sys; from pathlib import Path; HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[3]
CASE_PATH, PINS_PATH = HERE / "case.json", HERE / "SOURCE_PINS.json"

class GateFailure(Exception):
    def __init__(self, gate, detail):
        self.gate, self.detail = gate, detail; super().__init__(gate + ": " + str(detail))
def need(ok, gate, detail):
    if not ok:
        raise GateFailure(gate, detail)

def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def bad_constant(value):
    raise ValueError("non-finite JSON constant: " + value)

def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"), parse_constant=bad_constant)

def num(x, label):
    need(isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(float(x)), "finite", label); return float(x)

def near(a, b, atol, rtol, gate):
    a, b = num(a, gate), num(b, gate)
    need(abs(a-b) <= atol + rtol*max(abs(a), abs(b)), gate,
         {"actual": a, "expected": b, "abs_error": abs(a-b)})

def close(a, b, atol, rtol): return abs(a-b) <= atol + rtol*max(abs(a), abs(b))

def enear(a, b, case, gate):
    near(a, b, case["energy_contract"]["absolute_tolerance_j"],
         case["energy_contract"]["relative_tolerance"], gate)

def vnear(a, b, atol, rtol, gate):
    need(isinstance(a, list) and isinstance(b, list) and len(a) == len(b), gate, "vector shape")
    for i, (x, y) in enumerate(zip(a, b)):
        near(x, y, atol, rtol, gate + "[" + str(i) + "]")

def mnear(a, b, rows, cols, gate):
    need(isinstance(a, list) and len(a) == rows, gate, "matrix row count"); need(isinstance(b, list) and len(b) == rows, gate, "reference matrix row count")
    for i in range(rows):
        need(isinstance(a[i], list) and len(a[i]) == cols, gate, {"row": i}); vnear(a[i], b[i], 0.0, 0.0, gate + "[" + str(i) + "]")

def dot(a, b):
    return sum(x*y for x, y in zip(a, b))

def mv(a, x):
    return [dot(row, x) for row in a]

def tmv(a, x):
    return [sum(a[i][j]*x[i] for i in range(len(a))) for j in range(len(a[0]))]

def inf(x):
    return max((abs(v) for v in x), default=0.0)

def matrix(a, rows, cols, label):
    need(isinstance(a, list) and len(a) == rows, label, "row count")
    for i, row in enumerate(a):
        need(isinstance(row, list) and len(row) == cols, label, {"row": i})
        for j, x in enumerate(row):
            num(x, label + "[" + str(i) + "," + str(j) + "]")

def envelope(kappa, law):
    k, d0, df, tp = law["k0_pa_per_m"], law["delta0_m"], law["deltaf_m"], law["tpeak_pa"]
    if kappa <= d0:
        return 0.5*k*kappa*kappa
    if kappa < df:
        s = kappa-d0; return 0.5*k*d0*d0 + tp*(s-0.5*s*s/(df-d0))
    return law["gc_j_per_m2"]

def damage(kappa, law):
    d0, df = law["delta0_m"], law["deltaf_m"]
    if kappa <= d0: return 0.0
    if kappa >= df: return 1.0
    return min(1.0, max(0.0, df*(kappa-d0)/(kappa*(df-d0))))

def path_work(start, end, old_k, old_d, law):
    if start == end: return 0.0
    sign, lo, hi = (1.0, start, end) if end > start else (-1.0, end, start); cut = min(hi, max(lo, old_k)); fixed = 0.5*law["k0_pa_per_m"]*(1-old_d)*(cut*cut-lo*lo)
    soft = envelope(hi, law)-envelope(max(lo, old_k), law) if hi > old_k else 0.0; return sign*(fixed+soft)

def split_work(a, b, old_k, old_d, law):
    lo, hi = min(a,b), max(a,b); cuts = [a,b] + [x for x in (0.5*(a+b),old_k,law["delta0_m"],law["deltaf_m"]) if lo < x < hi]; cuts.sort(reverse=b < a)
    return sum(path_work(x,y,old_k,old_d,law) for x,y in zip(cuts,cuts[1:]))

def decode(blob, count, label):
    raw = base64.b64decode(blob, validate=True); need(len(raw) == count*8, "binary_length", label); values = [v[0] for v in struct.iter_unpack("<d", raw)]
    need(all(math.isfinite(x) for x in values), "binary_finite", label); return raw, values

def source_contract(case, pins):
    for name, path in (("checker", HERE/"check_saved_opening_history.py"),
                       ("case", CASE_PATH), ("readme", HERE/"README.md"), ("design", HERE/"DESIGN.md")):
        need(sha(path) == pins["source_files"][name]["sha256"], "own_source_pin", name)
    for group in ("source_pins", "producer_support_files"):
        for name, ref in case[group].items():
            p = ROOT/ref["path"]; need(p.is_file() and sha(p) == ref["sha256"], "input_source_pin", name)
    need(str(ROOT) == case["workspace_root_literal"], "workspace_root", str(ROOT))

def receipt_gate(case, receipt_path, result_path, result_sha):
    rec = read_json(receipt_path); need(rec.get("schema") == case["producer"]["capture_receipt_schema"], "receipt_schema", rec.get("schema"))
    need(rec.get("status") == "CAPTURED" and rec.get("actual_returncode") == 0 and rec.get("timed_out") is False
         and rec.get("child_started") is True, "receipt_status", rec.get("status"))
    need(rec.get("child_root_go_token_forwarded_exactly") is True and rec.get("pre_post_hashes_match") is True
         and rec.get("checker_run") is False, "receipt_integrity", "GO/source/checker flags")
    ri = rec.get("result", {})
    need(ri.get("exists") is True and ri.get("sha256") == result_sha
         and ri.get("bytes") == Path(result_path).stat().st_size, "receipt_result_binding", ri)
    pre, post = rec.get("source_hashes_pre", {}), rec.get("source_hashes_post", {}); need(bool(pre) and pre == post, "receipt_pre_post", "changed or empty source hash map")
    known = {"capture_case":case["producer"]["capture_case_sha256"],
             "capture_driver":case["producer"]["capture_wrapper_sha256"],
             "opening_case":case["producer"]["case_sha256"], "opening_driver":case["producer"]["driver_sha256"],
             "disabled_snapshot":case["producer"]["source_pins_disabled_snapshot_sha256"]}
    ref_alias = {"producer_case":"case", "producer_driver":"driver"}
    for name, ref in case["source_pins"].items():
        known["reference:"+ref_alias.get(name,name)] = ref["sha256"]
    for ref in case["producer_support_files"].values():
        known[ref["receipt_key"]] = ref["sha256"]
    for name, digest in known.items():
        e = pre.get(name,{}); need(e.get("exists") is True and e.get("sha256") == digest, "receipt_source_hash", name)
    for name in ("trace_primitive","cohesive_scalar_law","cohesive_derivative"):
        e = pre.get("trace_dependency:"+name,{})
        need(e.get("exists") is True and e.get("sha256") == case["source_pins"][name]["sha256"],
             "receipt_trace_hash", name)
    run = Path(result_path).parent; sizes = {p.relative_to(run).as_posix():p.stat().st_size for p in run.rglob("*") if p.is_file()}; total = sum(sizes.values())
    need(sizes == rec.get("file_tally_bytes") and total == rec.get("combined_bytes_including_receipt")
         and total <= rec.get("combined_file_limit_bytes",-1) and rec.get("combined_file_limit_pass") is True,
         "receipt_output_budget", total)
    return rec

def audit(case, result, result_path, receipt_path):
    b = case["pending_root_binding"]
    for key in ("opening_result_path","opening_result_sha256","opening_receipt_path","opening_receipt_sha256","independent_output_path"):
        need(b.get(key) is not None, "root_binding", key+" is pending")
    need(result_path.resolve() == (ROOT/b["opening_result_path"]).resolve()
         and receipt_path.resolve() == (ROOT/b["opening_receipt_path"]).resolve(),
         "root_binding_paths", "input paths differ")
    rs = sha(result_path)
    need(rs == b["opening_result_sha256"] and sha(receipt_path) == b["opening_receipt_sha256"],
         "opening_input_hash", "result or receipt SHA mismatch")
    receipt_gate(case, receipt_path, result_path, rs)
    need(result.get("schema") == case["producer"]["result_schema"] and result.get("status") == "DIAGNOSTIC_ONLY"
         and result.get("classification") == case["classification"], "result_identity", "schema/status/class")
    need(result.get("execution_authorized_by_root_go") is True, "producer_go", "missing execution marker"); counts = result["counts"]
    for k,v in case["expected_counts"].items():
        need(counts.get(k) == v, "producer_count", {k:counts.get(k),"expected":v})
    src = dict(case["source_pins"])
    for name in ("producer_case","producer_driver"):
        src[{"producer_case":"case","producer_driver":"driver"}[name]] = src.pop(name)
    for k,ref in src.items():
        need(result.get("source_hashes",{}).get(k) == ref["sha256"], "result_source_hash", k)
    ep = ROOT/case["elastic_input"]["path"]; need(sha(ep) == case["elastic_input"]["sha256"], "elastic_input_sha", str(ep)); elastic = read_json(ep)
    need(elastic.get("schema") == case["elastic_input"]["schema"] and elastic.get("status") == case["elastic_input"]["status"],
         "elastic_input_identity", str(ep))
    n, ng = case["operator"]["coefficient_count"], case["operator"]["physical_gauge_count"]
    cidx, geom, law = case["operator"]["control_index_zero_based"], case["geometry"], case["material_fixture"]
    tol, ftol, fr = case["tolerances"], case["tolerances"]["force_abs_n"], case["tolerances"]["force_rel"]; need(close_material(law), "material_identity", law)
    cache, bind, resp = elastic["evidence"]["cache"], result["saved_input_binding"], result["unit_control_response"]; K, C = cache["K_raw_N_per_m"], cache["C_phys"]
    matrix(K,n,n,"K"); matrix(C,ng,n,"C")
    mnear(bind["K_raw_N_per_m"],K,n,n,"bound_K"); mnear(bind["C_phys"],C,ng,n,"bound_C")
    need(bind["coefficient_ids"] == cache["semantic_ids"] and bind["control_semantic_id"] == cache["semantic_ids"][cidx],
         "semantic_ids","coefficient order")
    need(bind["control_index_zero_based"] == cidx, "control_index", bind["control_index_zero_based"])
    vnear(bind["physical_external_force_N"],[0.0]*n,0.0,0.0,"physical_external_force"); vnear(bind["physical_gauge_target_m"],[0.0]*ng,0.0,0.0,"physical_gauge_target")
    mnear(bind["source_positions_xyz_m"],cache["source_positions_xyz_m"],27,3,"source_positions")
    mnear(bind["rigid_coefficient_modes"],cache["rigid_coefficient_modes"],6,n,"rigid_modes")
    mnear(resp["K_raw_N_per_m"],K,n,n,"response_K"); mnear(resp["C_phys"],C,ng,n,"response_C")
    qu, lu = resp["q_per_Q"], resp["physical_gauge_lambda_per_Q"]; need(len(qu)==n and len(lu)==ng, "unit_response_shape", "q/lambda")
    near(qu[cidx],1.0,1e-12,1e-10,"unit_control"); unit_res = [a+b for a,b in zip(mv(K,qu),tmv(C,lu))]; need(inf(mv(C,qu)) <= tol["gauge_inf_m"], "unit_gauge", inf(mv(C,qu)))
    need(inf(unit_res[:cidx]+unit_res[cidx+1:]) <= tol["unit_force_abs_n_per_m"]+fr*max(inf(unit_res),1.0),
         "unit_free_equilibrium",inf(unit_res[:cidx]+unit_res[cidx+1:]))
    near(resp["bulk_applied_control_force_per_Q_N_per_m"],unit_res[cidx],tol["unit_force_abs_n_per_m"],fr,"unit_reaction"); interfaces = elastic["evidence"]["interface"]
    need(len(interfaces)==geom["quadrature_points"]==16,"interface_count",len(interfaces)); bound = bind["interfaces_reconstructed_from_saved_Hplus_Hminus"]
    need(len(bound)==16,"bound_interface_count",len(bound)); nvec, areas, alphas, trace = geom["normal_lower_to_upper"], [], [], []
    for i,(src,saved) in enumerate(zip(interfaces,bound)):
        hp,hm=src["Hplus"],src["Hminus"]; matrix(hp,3,n,"Hplus"); matrix(hm,3,n,"Hminus")
        J=[[hp[r][j]-hm[r][j] for j in range(n)] for r in range(3)]
        need(saved["quadrature_id"]==src["quadrature_id"],"saved_qp_id",i); vnear(saved["point_xyz_m"],src["point_xyz_m"],0.0,0.0,"point_coordinate")
        vnear(saved["normal_lower_to_upper"],nvec,0.0,0.0,"normal"); mnear(saved["Hplus"],hp,3,n,"saved_Hplus"); mnear(saved["Hminus"],hm,3,n,"saved_Hminus")
        mnear(saved["J_reconstructed"],J,3,n,"J"); nr=[sum(nvec[r]*J[r][j] for r in range(3)) for j in range(n)]
        vnear(saved["normal_row"],nr,tol["analytic_normal_abs"],0.0,"normal_row"); x=num(src["point_xyz_m"][0],"x"); R=geom["radius_m"]
        need(-R<=x<=0.0,"interface_x",x); alpha=2.0*math.sqrt(-x/R)
        near(nr[cidx],alpha,tol["analytic_normal_abs"],1e-14,"analytic_normal_control")
        for j,value in enumerate(nr):
            if j != cidx: near(value,0.0,tol["analytic_normal_abs"],0.0,"analytic_normal_other["+str(j)+"]")
        near(saved["analytic_alpha"],alpha,tol["analytic_normal_abs"],1e-14,"alpha_bound"); x,y,z=src["point_xyz_m"]
        expected=[(1.,0.,0.),(0.,1.,0.),(0.,0.,1.),(0.,-z,y),(z,0.,-x),(-y,x,0.)]
        source_expected=src["expected_rigid_displacements_m"]
        mnear(saved["expected_rigid_displacements_m"],source_expected,6,3,"rigid_expected_source")
        for j,(saved_disp,physical_disp) in enumerate(zip(saved["expected_rigid_displacements_m"],expected)):
            vnear(saved_disp,list(physical_disp),tol["trace_abs_m"],tol["trace_rel"],"rigid_expected_physical["+str(j)+"]")
        for mode,disp in zip(cache["rigid_coefficient_modes"],expected):
            vnear(mv(hp,mode),list(disp),tol["trace_abs_m"],tol["trace_rel"],"Hplus_rigid"); vnear(mv(hm,mode),list(disp),tol["trace_abs_m"],tol["trace_rel"],"Hminus_rigid")
            need(inf(mv(J,mode))<=tol["trace_abs_m"],"jump_rigid_annihilation",i)
        area=num(src["one_side_area_weight_m2"],"area"); need(area>0.0,"positive_area",i); near(saved["one_side_area_weight_m2"],area,0.0,0.0,"bound_area")
        areas.append(area); alphas.append(alpha); trace.append((J,nr,area,alpha,src["quadrature_id"]))
    ids=[x["quadrature_id"] for x in interfaces]; need(len({json.dumps(x,sort_keys=True) for x in ids})==16,"unique_qp_ids",ids)
    near(sum(areas),geom["one_side_area_m2"],1e-15,1e-12,"one_side_area"); near(geom["one_side_area_m2"],geom["radius_m"]*geom["length_m"],1e-15,1e-12,"R_times_L")
    need(len(resp["alpha_i"])==16,"response_alpha_count",len(resp["alpha_i"]))
    for a,bv in zip(resp["alpha_i"],alphas): near(a,bv,tol["analytic_normal_abs"],1e-14,"response_alpha")
    qstar=law["deltaf_m"]/min(alphas); near(result["Qstar_m"],qstar,1e-18,1e-12,"Qstar"); frac=case["trajectory_fraction_Qstar"]
    need(result["trajectory_fraction_Qstar"]==frac and len(result["trajectory"])==12,"trajectory","frozen path"); binary=result["trajectory_binary_le_f64"]
    need(binary.get("format")=="IEEE-754 little-endian float64, base64" and binary.get("q_physics_shape")==[12,n]
         and binary.get("damage_shape")==[12,16],"binary_schema","format/shape")
    qb,qflat=decode(binary["q_physics"],12*n,"q"); db,df=decode(binary["damage"],192,"damage")
    kb,kflat=decode(binary["kappa_m"],192,"kappa"); ob,oflat=decode(binary["opening_m"],192,"opening")
    qbulk=dot(qu,mv(K,qu)); need(qbulk>0.0,"bulk_stiffness",qbulk)
    previous=[[0.0,0.0,0.0,0.0,0.0] for _ in range(16)]; previous_q=previous_bulk=cumulative=0.0; max_res=max_suberr=0.0; checkpoints={}; partial_ref=None; full_ref=None
    for s,row in enumerate(result["trajectory"]):
        need(row["step_index"]==s,"step_order",row.get("step_index")); Q=frac[s]*qstar; near(row["Q_fraction_of_Qstar"],frac[s],0.0,0.0,"fraction")
        near(row["Q_m"],Q,1e-18,1e-12,"Q"); q=[Q*x for x in qu]; vnear(row["q_physics_m"],q,tol["trace_abs_m"],tol["trace_rel"],"q")
        for j,x in enumerate(row["q_physics_m"]):
            need(struct.pack("<d",num(x,"qjson"))==qb[(s*n+j)*8:(s*n+j+1)*8],"binary_q",(s,j))
        bulk=mv(K,q); gauge=[Q*x for x in tmv(C,lu)]
        vnear(row["bulk_internal_generalized_force_N"],bulk,ftol,fr,"bulk_force"); vnear(row["physical_gauge_generalized_force_N"],gauge,ftol,fr,"gauge_force")
        vnear(row["bulk_plus_gauge_residual_N"],[a+b for a,b in zip(bulk,gauge)],ftol,fr,"bulk_gauge"); gamma=[0.0]*n; st=ds=work=prim=pcum=0.0; states=[]; deltas=[]; dam=[]
        pts=row["points"]; need(len(pts)==16,"point_count",s)
        for i,(pt,item) in enumerate(zip(pts,trace)):
            J,nr,area,alpha,qid=item; need(isinstance(pt["point_id"],list) and pt["point_id"][-1]==qid,"point_id",(s,i))
            jump=[sum(J[r][j]*q[j] for j in range(n)) for r in range(3)]; delta=dot(nvec,jump); need(delta>=-tol["trace_abs_m"],"negative_opening",(s,i,delta))
            tang=[jump[r]-delta*nvec[r] for r in range(3)]
            need(math.sqrt(dot(tang,tang))<=tol["trace_abs_m"]+tol["trace_rel"]*max(abs(delta),1e-30),
                 "tangential_slip",(s,i))
            near(delta,alpha*Q,tol["trace_abs_m"],tol["trace_rel"],"trace"); near(pt["opening_m"],delta,tol["trace_abs_m"],tol["trace_rel"],"opening")
            near(pt["normal_opening_recomputed_m"],delta,tol["trace_abs_m"],tol["trace_rel"],"normal_opening")
            vnear(pt["jump_vector_m"],jump,tol["trace_abs_m"],tol["trace_rel"],"jump"); vnear(pt["tangential_jump_vector_m"],tang,tol["trace_abs_m"],tol["trace_rel"],"tangent")
            near(pt["analytic_alpha"],alpha,tol["analytic_normal_abs"],1e-14,"point_alpha"); od,ok,odmg,os,oD=previous[i]; k=max(ok,delta); d=damage(k,law)
            t=(1-d)*law["k0_pa_per_m"]*delta; S=0.5*(1-d)*law["k0_pa_per_m"]*delta*delta
            D=max(0.0,envelope(k,law)-0.5*(1-d)*law["k0_pa_per_m"]*k*k); w=path_work(od,delta,ok,odmg,law); ws=split_work(od,delta,ok,odmg,law)
            enear(area*w,area*ws,case,"piecewise_split_work"); max_suberr=max(max_suberr,abs(area*(w-ws)))
            state=pt["state"]; near(state["opening_m"],delta,tol["trace_abs_m"],tol["trace_rel"],"state_opening")
            near(state["kappa_m"],k,tol["trace_abs_m"],tol["trace_rel"],"state_kappa"); near(state["damage"],d,tol["damage_abs"],1e-12,"state_damage")
            enear(area*state["stored_energy_j_per_m2"],area*S,case,"stored_per_area"); enear(area*state["dissipated_energy_j_per_m2"],area*D,case,"diss_per_area")
            need(state["update_count"]==s+1,"update_count",(s,i)); near(pt["traction_pa"],t,1e-9,2e-10,"traction")
            need(k+1e-18>=ok and d+tol["damage_abs"]>=odmg and D+1e-16>=oD,"irreversible",(s,i))
            for val,raw,label in ((state["damage"],db,"damage"),(state["kappa_m"],kb,"kappa"),
                                  (state["opening_m"],ob,"opening")):
                need(struct.pack("<d",num(val,label))==raw[(s*16+i)*8:(s*16+i+1)*8],"binary_json_match",(s,i,label))
            pf=[area*nr[j]*t for j in range(n)]; vnear(pt["generalized_internal_force_N"],pf,ftol,fr,"qp_generalized_force"); gamma=[a+b for a,b in zip(gamma,pf)]
            pos=[area*t*x for x in nvec]; vnear(pt["positive_face_internal_force_N"],pos,ftol,fr,"positive_face")
            vnear(pt["negative_face_internal_force_N"],[-x for x in pos],ftol,fr,"negative_face")
            vnear(pt["upper_physical_traction_N_per_m2"],[-t*x for x in nvec],1e-8,1e-10,"upper_traction")
            vnear(pt["lower_physical_traction_N_per_m2"],[t*x for x in nvec],1e-8,1e-10,"lower_traction")
            enear(pt["cohesive_work_increment_primitive_j"],area*w,case,"point_primitive_work"); enear(pt["cohesive_work_increment_independent_law_j"],area*w,case,"point_law_work")
            enear(area*state["cumulative_external_work_j_per_m2"],area*(S+D),case,"cumulative_point_work")
            st+=area*S; ds+=area*D; work+=area*w; prim+=pt["cohesive_work_increment_primitive_j"]
            pcum+=area*state["cumulative_external_work_j_per_m2"]; states.append([delta,k,d,S,D]); deltas.append(delta); dam.append(d)
        force=[a+b+c for a,b,c in zip(bulk,gauge,gamma)]; vnear(row["cohesive_generalized_internal_force_N"],gamma,ftol,fr,"cohesive_force"); reaction=force[cidx]
        near(row["control_applied_force_positive_opening_N"],reaction,ftol,fr,"control_reaction"); need(reaction>=-(ftol+fr*max(abs(reaction),1.0)),"control_force_sign",reaction)
        residual=force[:cidx]+force[cidx+1:]; max_res=max(max_res,inf(residual)); need(inf(residual)<=ftol+fr*max(abs(reaction),inf(force),1.0),"free_force_balance",inf(residual))
        near(row["noncontrol_force_inf_N"],inf(residual),ftol,fr,"saved_force_residual")
        for field,value in (("cohesive_stored_energy_j",st),("cohesive_dissipated_energy_j",ds),
                            ("cohesive_work_increment_independent_law_j",work),
                            ("cohesive_work_increment_primitive_j",prim),
                            ("primitive_aggregate_cumulative_external_work_j",pcum)):
            enear(row[field],value,case,field)
        U=0.5*dot(q,bulk); Uw=0.5*qbulk*(Q*Q-previous_q*previous_q)
        enear(row["bulk_stored_energy_j"],U,case,"bulk_storage"); enear(row["bulk_work_increment_exact_j"],Uw,case,"exact_bulk_work"); mid=0.5*(Q+previous_q)
        enear(Uw,0.5*qbulk*(mid*mid-previous_q*previous_q)+0.5*qbulk*(Q*Q-mid*mid),case,"bulk_subdivision")
        W=Uw+work; enear(row["control_work_increment_exact_j"],W,case,"total_work")
        oldS=sum(areas[i]*previous[i][3] for i in range(16)); oldD=sum(areas[i]*previous[i][4] for i in range(16))
        balance=(U-previous_bulk)+(st-oldS)+(ds-oldD); enear(W,balance,case,"global_energy_balance")
        cumulative+=W; enear(row["control_work_cumulative_exact_j"],cumulative,case,"cumulative_control_work")
        near(row["energy_balance_relative_error"],abs(W-balance)/max(abs(W),abs(balance),1e-300),2e-12,2e-9,"recorded_energy_error")
        if s==3:
            need(any(0<d<1 for d in dam),"partial_damage_step3",dam); partial_ref=[x[:] for x in states]
        if s in (4,5):
            need(all(close(states[i][1],partial_ref[i][1],1e-18,0.0) and close(states[i][2],partial_ref[i][2],2e-12,0.0)
                     and close(states[i][4],partial_ref[i][4],1e-16,0.0) for i in range(16)),"partial_history_preserved",s)
            if s==4: need(all(deltas[i]<partial_ref[i][0] for i in range(16)),"unload_path",s)
            else: need(all(close(deltas[i],partial_ref[i][0],tol["trace_abs_m"],tol["trace_rel"]) for i in range(16)),
                       "reload_path",s)
            need(any(0<d<1 for d in dam),"partial_unload_reload",s)
        if s==8:
            need(all(d==1.0 for d in dam),"all_full_step8",dam); full_ref=[x[:] for x in states]
            checkpoints["full_failure_dissipation_j"]=ds
        if s==10:
            need(all(d==1.0 for d in dam) and all(x==0.0 for x in deltas),"zero_opening_no_heal",s)
            need(all(close(states[i][1],full_ref[i][1],1e-18,0.0)
                     and close(states[i][4],full_ref[i][4],1e-16,0.0) for i in range(16)),"full_history_step10",s)
        if s==11:
            need(all(d==1.0 and pt["traction_pa"]==0.0 for d,pt in zip(dam,pts)),"full_reload_zero_traction",s)
            need(all(close(states[i][1],full_ref[i][1],1e-18,0.0)
                     and close(states[i][4],full_ref[i][4],1e-16,0.0) for i in range(16)),"full_history_step11",s)
        if s in (3,4,5,8,10,11):
            checkpoints["step_"+str(s)]={"Q_fraction":frac[s],"partial_count":sum(0<d<1 for d in dam),
                "min_damage":min(dam),"max_damage":max(dam),"stored_energy_j":st,"dissipated_energy_j":ds,
                "cohesive_work_increment_j":work,"control_reaction_N":reaction}
        previous,previous_q,previous_bulk=states,Q,U
    finalD=checkpoints["full_failure_dissipation_j"]; areaGc=geom["one_side_area_m2"]*law["gc_j_per_m2"]
    enear(finalD,areaGc,case,"one_side_area_Gc"); near(finalD,1e-6,1e-22,2e-9,"total_Gc")
    need(checkpoints["step_10"]["dissipated_energy_j"]==checkpoints["step_8"]["dissipated_energy_j"]
         and checkpoints["step_11"]["dissipated_energy_j"]==checkpoints["step_8"]["dissipated_energy_j"],
         "full_nonhealing","dissipation changed after full failure")
    return {"schema":"mixed93_opening_history_independent_result.v1","status":"ACCEPTED_DIAGNOSTIC_ONLY",
        "classification":case["classification"],"execution_scope":"saved JSON only; no solver, NPZ, provider, or author imports",
        "limitations":["toy material fixture, not ceramic calibration","no crack growth or front advance",
                       "path subdivision is analytic bookkeeping, not solver convergence"],
        "bindings":{"opening_result_sha256":rs,"opening_receipt_sha256":sha(receipt_path),
                    "elastic_result_sha256":case["elastic_input"]["sha256"]},
        "gates":{"capture_receipt":"PASS","binary_json_consistency":"PASS","interfaces_and_opening":"PASS",
                 "unit_control_equilibrium":"PASS","pointwise_law_and_history":"PASS","force_balance_and_signs":"PASS",
                 "exact_path_work_and_energy":"PASS","required_checkpoints":"PASS","one_side_area_Gc":"PASS"},
        "counts":{"steps":12,"quadrature_points":16,"point_state_checks":192,"solver_or_provider_calls":0,"npz_content_reads":0},
        "metrics":{"one_side_area_m2":sum(areas),"full_failure_dissipation_j":finalD,
                   "maximum_noncontrol_force_residual_n":max_res,"maximum_branch_subdivision_energy_error_j":max_suberr,
                   "energy_tolerance":case["energy_contract"]},"checkpoints":checkpoints}

def close_material(law):
    return (close(law["delta0_m"],law["tpeak_pa"]/law["k0_pa_per_m"],1e-18,1e-14)
            and close(law["deltaf_m"],2*law["gc_j_per_m2"]/law["tpeak_pa"],1e-18,1e-14)
            and law["deltaf_m"]>law["delta0_m"])

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__); p.add_argument("--execute",action="store_true"); p.add_argument("--root-go")
    p.add_argument("--result",type=Path); p.add_argument("--receipt",type=Path)
    a=p.parse_args(argv)
    if not a.execute:
        print(json.dumps({"status":"EXECUTION_DISABLED","npz_opened":False,"author_imports":False})); return 2
    try:
        case,pins=read_json(CASE_PATH),read_json(PINS_PATH); go=pins.get("root_go",{})
        need(pins.get("execution_authorized") is True and go.get("authorized") is True
             and isinstance(go.get("token"),str) and a.root_go==go["token"],"root_go","exact independent-checker GO required")
        source_contract(case,pins); need(a.result is not None and a.receipt is not None,"cli_inputs","bound result and receipt required")
        root=ROOT.resolve(); rp=a.result.resolve(); cp=a.receipt.resolve()
        rp.relative_to(root); cp.relative_to(root)
        out=(ROOT/case["pending_root_binding"]["independent_output_path"]).resolve(); out.relative_to(root)
        need(not out.exists(),"exclusive_output",str(out)); result=read_json(rp); audited=audit(case,result,rp,cp)
        raw=(json.dumps(audited,indent=2,sort_keys=True,allow_nan=False)+"\n").encode("utf-8")
        with out.open("xb") as f: f.write(raw)
        print(json.dumps({"status":audited["status"],"output":str(out),"bytes":len(raw),
                          "opening_result_sha256":audited["bindings"]["opening_result_sha256"]},allow_nan=False))
        return 0
    except Exception as exc:
        failure={"schema":"mixed93_opening_history_independent_result.v1","status":"NOT_ACCEPTED",
                 "classification":"DIAGNOSTIC_ONLY","failure":{"gate":getattr(exc,"gate","audit"),
                 "type":type(exc).__name__,"detail":getattr(exc,"detail",str(exc))},
                 "limitations":["failed saved JSON audit","no material or growth qualification"]}
        print(json.dumps(failure,allow_nan=False)); return 1
if __name__=="__main__":
    sys.exit(main())
