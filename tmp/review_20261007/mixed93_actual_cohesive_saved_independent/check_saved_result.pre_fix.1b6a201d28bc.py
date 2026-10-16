"""Root-GO-gated independent audit of one frozen JSON/NPZ result; never a runner."""
from __future__ import annotations
import argparse, hashlib, json, math
from pathlib import Path
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
CASE = json.loads((HERE / "case.json").read_text(encoding="utf-8"))
PINS = json.loads((HERE / "pins.json").read_text(encoding="utf-8"))
def sha(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda:stream.read(1<<20),b""): h.update(block)
    return h.hexdigest()
def require(ok: bool, message: str) -> None:
    if not ok: raise ValueError(message)
def close(np, label, a, b, rtol, atol=0.0) -> None:
    require(bool(np.allclose(a,b,rtol=rtol,atol=atol)),f"{label}: mismatch")
def hash_pins() -> dict:
    found = {}
    for name, ref in CASE["source_pins"].items():
        path=ROOT/ref["path"]; require(path.is_file(),f"missing pinned source {name}")
        actual=sha(path); require(actual==ref["sha256"],f"pinned source changed: {name}"); found[name]=actual
    cache = CASE["expected_hashes"]
    npz_path=ROOT/CASE["producer_npz"]; require(npz_path.is_file(),"accepted cached NPZ missing")
    npz_hash=sha(npz_path); require(npz_hash==cache["producer_npz"],"accepted cached NPZ SHA-256 mismatch")
    found["producer_npz"] = npz_hash
    return found
def point_check(np, point, q, ids, weights, case, label):
    g = case["geometry"]; mat = case["material_fixture"]
    x,y,z=point["point_xyz_m"]; area=float(point["one_side_area_weight_m2"])
    n=np.asarray(point["normal_unit_lower_to_upper"],dtype=float)
    H=np.asarray(point["Hplus"],dtype=float); Hm=np.asarray(point["Hminus"],dtype=float)
    J=np.asarray(point["J"],dtype=float); modes=np.asarray(point["expected_rigid_displacements_m"],dtype=float)
    require(H.shape == (3,93) and Hm.shape == (3,93) and J.shape == (3,93), f"{label} trace shape")
    close(np, f"{label} normal", n, [0.,1.,0.], 0., 1e-15)
    close(np, f"{label} jump identity", J, H-Hm, 0., 2e-14)
    a = math.sqrt(-x / g["R_m"])
    center=np.asarray(CASE["source_center_xyz_m"],dtype=float)
    expected_modes=np.asarray([np.eye(3)[i] for i in range(3)]+
        [np.cross(np.eye(3)[i],np.asarray(point["point_xyz_m"],dtype=float)-center) for i in range(3)])
    close(np,f"{label} geometric rigid displacement candidates",modes,expected_modes,0.,2e-13)
    expectedJ = np.zeros((3,93)); expectedJ[0,81] = expectedJ[1,82] = expectedJ[2,83] = 2*a
    close(np, f"{label} exact F1 J", J, expectedJ, 0., 2e-13)
    require(abs(y)<=1e-15 and -g["R_m"]<x<0. and abs(z)<g["L_m"]/2,f"{label} point outside open quadrature sheet")
    expected_weight=np.asarray(weights,dtype=float); eH=np.zeros((3,93)); eHm=np.zeros((3,93))
    base_cols = set()
    for si, sid in enumerate(ids):
        for ci, comp in enumerate(("ux","uy","uz")):
            key = json.dumps(["SOURCE_BASE",sid,comp,0], separators=(",",":"))
            col = CASE["semantic_index"][key]
            base_cols.add(col)
            eH[ci,col]=eHm[ci,col]=expected_weight[si]
    for ci, comp in enumerate(("ux","uy","uz")):
        col = CASE["semantic_index"][json.dumps(["GLOBAL_FRONT","global",comp,1],separators=(",",":"))]
        eH[ci,col]=a; eHm[ci,col]=-a
    require(len(base_cols) == 81 and base_cols == set(range(81)), f"{label} semantic base registry")
    close(np,f"{label} upper trace map",H,eH,0.,2e-13)
    close(np,f"{label} lower trace map",Hm,eHm,0.,2e-13)
    return J, n, area, modes
def audit_points(np, rows, q, rigid, weights, case, label, tol):
    require(len(rows) == 16, f"{label}: expected 16 quadrature rows")
    ids = CASE["source_ids"]
    total_area = 0.; total_energy = 0.; total_diss = 0.; opens=[]; max_rigid=0.; max_local=0.
    for qi, row in enumerate(rows):
        J,n,area,modes=point_check(np,row,q,ids,weights[qi],case,label)
        expected = modes.T
        max_rigid=max(max_rigid,float(np.max(np.abs(np.asarray(row["Hplus"])@rigid.T-expected))),
                      float(np.max(np.abs(np.asarray(row["Hminus"])@rigid.T-expected))),
                      float(np.max(np.abs(J@rigid.T))))
        opening=float(n@J@q); traction=float(row["traction_pa"]); k=float(row["k_alg_pa_per_m"])
        close(np,f"{label} ledger opening",row["opening_m"],opening,2e-9,1e-15)
        close(np,f"{label} recomputed normal opening",row["normal_opening_recomputed_m"],opening,2e-9,1e-15)
        close(np,f"{label} jump vector",row["jump_vector_m"],J@q,2e-9,1e-15)
        close(np,f"{label} tangential jump",row["tangential_jump_vector_m"],J@q-opening*n,0.,tol["rigid_jump_m"])
        state=row["state"]; expected_t=case["material_fixture"]["K0_Pa_per_m"]*opening
        require(0. < opening < case["material_fixture"]["delta0_m"], f"{label} opening outside virgin elastic range")
        require(traction <= case["material_fixture"]["tpeak_Pa"] and traction >= 0., f"{label} traction range")
        close(np,f"{label} traction",traction,expected_t,2e-9,1e-8)
        close(np,f"{label} tangent",k,case["material_fixture"]["K0_Pa_per_m"],2e-12)
        require(row["tangent_status"]=="differentiable" and not row["nondifferentiable_reasons"],f"{label} unexpected law kink")
        close(np,f"{label} state opening",state["opening_m"],opening,2e-9,1e-15)
        close(np,f"{label} state history",state["kappa_m"],opening,2e-9,1e-15)
        require(float(state["damage"]) == 0. and float(state["dissipated_energy_j_per_m2"]) == 0., f"{label} damage/dissipation")
        energy=0.5*case["material_fixture"]["K0_Pa_per_m"]*opening*opening
        close(np,f"{label} point stored energy",state["stored_energy_j_per_m2"],energy,2e-9,tol["absolute_energy_density_J_per_m2"])
        close(np,f"{label} point external work",state["external_work_j_per_m2"],energy,2e-9,tol["absolute_energy_density_J_per_m2"])
        local=np.asarray(row["generalized_internal_force_n"])
        expected_local=area*(J.T@n)*traction
        close(np,f"{label} local generalized force",local,expected_local,2e-9,1e-9)
        positive=np.asarray(row["primitive_positive_face_internal_force_n"])
        upper=np.asarray(row["physical_upper_face_traction_n_per_m2"])
        lower=np.asarray(row["physical_lower_face_traction_n_per_m2"])
        close(np,f"{label} upper physical traction sign",upper,-traction*n,2e-12,1e-10)
        close(np,f"{label} lower physical traction sign",lower,traction*n,2e-12,1e-10)
        close(np,f"{label} internal face force",positive,area*traction*n,2e-12,1e-10)
        close(np,f"{label} face action/reaction",area*upper+positive,[0.,0.,0.],0.,1e-10)
        max_local=max(max_local,float(np.max(np.abs(local-expected_local))))
        opens.append(opening); total_area+=area; total_energy+=area*energy
        total_diss+=area*float(state["dissipated_energy_j_per_m2"])
    require(max_rigid <= tol["rigid_displacement_m"],f"{label} rigid trace/jump residual")
    return {"area_m2":total_area,"open_min_m":min(opens),"open_max_m":max(opens),
            "stored_energy_J":total_energy,"dissipated_energy_J":total_diss,
            "max_rigid_residual_m":max_rigid,"max_local_force_residual_N":max_local}
GL_NODES=CASE["gauss_legendre_4"]["nodes"]
GL_WEIGHTS=CASE["gauss_legendre_4"]["weights"]
def run(result_path: Path, receipt_path: Path, token: str) -> dict:
    require(sha(HERE/"check_saved_result.py")==PINS["saved_check_sha256"] and sha(HERE/"check_source_contract.py")==PINS["source_check_sha256"] and sha(HERE/"case.json")==PINS["case_sha256"],"independent source package changed after review")
    require(PINS["execution_authorized"] is True and PINS["root_go"].get("authorized") is True and isinstance(PINS["root_go"].get("token"),str) and token==PINS["root_go"]["token"],"independent saved audit is disabled until exact root GO")
    require(PINS["result"]["sha256"] and PINS["receipt"]["sha256"], "root must freeze result and receipt hashes first")
    require(result_path==ROOT/PINS["result"]["path"] and receipt_path==ROOT/PINS["receipt"]["path"],"input paths differ from frozen pins")
    source_hashes=hash_pins()
    current_worker_pins=json.loads((ROOT/PINS["author_numeric_source_pins"]["path"]).read_text(encoding="utf-8"))
    require(current_worker_pins.get("execution_authorized") is True and current_worker_pins.get("root_go",{}).get("authorized") is True and current_worker_pins["root_go"].get("token")==PINS["author_numeric_source_pins"]["producer_root_go_token"],"producer numerical run authorization changed")
    for name,ref in CASE["source_pins"].items():
        require(current_worker_pins.get("references",{}).get(name)==ref,f"worker pin reference changed: {name}")
    binding=current_worker_pins.get("accepted_cache_binding",{})
    require(binding.get("arrays_sha256")==CASE["expected_hashes"]["producer_npz"] and binding.get("audit_sha256")==CASE["expected_hashes"]["recovery_audit"] and binding.get("receipt_sha256")==CASE["expected_hashes"]["recovery_receipt"] and binding.get("actual_array_keys")==CASE["cache_array_keys"],"worker accepted-cache binding changed")
    result_sha,receipt_sha=sha(result_path),sha(receipt_path)
    require(result_sha==PINS["result"]["sha256"] and receipt_sha==PINS["receipt"]["sha256"], "saved worker artifacts changed")
    r=json.loads(result_path.read_text(encoding="utf-8")); receipt=json.loads(receipt_path.read_text(encoding="utf-8"))
    require(r.get("schema")=="mixed93_actual_cohesive_bridge_numeric_result.v1" and r.get("status")=="ACCEPTED_DIAGNOSTIC_ONLY" and r.get("classification")=="FIXED_MODE_I_ELASTIC_COHESIVE_BRIDGE_ONLY_NOT_MATERIAL_OR_GROWTH_QUALIFICATION" and r.get("npz_opened") is True and r.get("numeric_modules_imported") is True,"worker result status/classification/flags")
    require(receipt.get("outer_return_code")==0 and receipt.get("timed_out") is False and receipt.get("launch_error") is None, "outer receipt indicates failed/timeout launch")
    require(receipt.get("source_pins_equal") is True and receipt.get("numeric_result_status")==r["status"], "outer receipt source/result binding")
    require(receipt.get("files",{}).get("result.json",{}).get("sha256")==result_sha and receipt.get("files",{}).get("result.json",{}).get("bytes")==result_path.stat().st_size,"receipt does not bind exact result bytes")
    require(receipt.get("source_pre_pins",{}).get("references")==receipt.get("source_post_pins",{}).get("references"),"receipt source pre/post pins differ")
    require(r.get("inputs",{}).get("npz_sha256")==CASE["expected_hashes"]["producer_npz"],"result input cache binding mismatch")
    for name in CASE["required_gates"]:
        gate=r.get("gates",{}).get(name)
        require(isinstance(gate,dict) and gate.get("passed") is True and gate.get("status")=="PASS",f"producer gate not PASS: {name}")
    require(r["counts"].get("new_mls_queries")==16 and r["counts"].get("mls_queries_requested")==16 and r["counts"].get("mls_query_attempts")==1 and r["counts"].get("interface_qps")==16 and r["counts"].get("primitive_assemblies")==3 and r["counts"].get("primitive_successful_returns")==3 and r["counts"].get("coupled_kkt_solves")==1 and r["counts"].get("coupled_kkt_successful_returns")==1,"worker call budget/count mismatch")
    import numpy as np
    npz_path=ROOT/CASE["producer_npz"]
    with np.load(npz_path,allow_pickle=False) as archive:
        arrays={name:np.asarray(archive[key]) for name,key in CASE["cache_array_keys"].items()}
    c=r["evidence"]["cache"]; tol=CASE["tolerances"]; geo=CASE["geometry"]
    require(c["arrays_sha256"]==CASE["expected_hashes"]["producer_npz"],"cache evidence SHA mismatch")
    require(c["actual_array_keys"]==CASE["cache_array_keys"] and c["independent_audit_sha256"]==CASE["expected_hashes"]["recovery_audit"] and c["recovery_receipt_sha256"]==CASE["expected_hashes"]["recovery_receipt"],"cache audit/receipt/key binding")
    def strings(a): return [v.decode("utf-8") if isinstance(v,bytes) else str(v) for v in a.tolist()]
    for key,field in (("K","K_raw_N_per_m"),("physical_constraint_matrix","C_phys"),("rigid_vectors","rigid_coefficient_modes"),("q_solution","q_solution_KI100_m"),("q_exact","q_exact_KI100_m_diagnostic_only"),("external_force","f_ext_KI100_N"),("physical_gauge_target","g_phys_KI100_m"),("physical_lambda","lambda_phys_KI100")):
        close(np,f"cached {key} vs evidence",arrays[key],c[field],0.,0.)
    source_ids=strings(arrays["source_ids"]); semantic_ids=strings(arrays["semantic_ids"])
    CASE["source_ids"]=source_ids
    require(source_ids==c["source_ids"]==r["inputs"]["source_ids"] and semantic_ids==c["semantic_ids"] and r["inputs"]["q_solution_scale"]==CASE["linear_scale"] and r["inputs"]["q_exact_kept_diagnostic_only"] is True,"cached ID/scale binding mismatch")
    close(np,"cached source positions",arrays["source_positions"],c["source_positions_xyz_m"],0.,0.)
    require(arrays["K"].shape==(93,93) and arrays["physical_constraint_matrix"].shape==(6,93) and arrays["rigid_vectors"].shape==(6,93),"cached matrix dimensions changed")
    producer_case=json.loads((ROOT/CASE["source_pins"]["case"]["path"]).read_text(encoding="utf-8"))
    source_binding=producer_case["source_binding"]
    require(source_ids==source_binding["source_ids"],"cached source order differs from frozen producer case")
    close(np,"cached source coordinates vs frozen case",arrays["source_positions"],source_binding["source_positions_xyz_m"],0.,0.)
    require(semantic_ids==producer_case["registry"]["semantic_coefficient_ids_in_order"],"cached semantic registry differs from frozen producer case")
    scale=float(CASE["linear_scale"]); q0=arrays["q_solution"]*scale; f=arrays["external_force"]*scale; g=arrays["physical_gauge_target"]*scale; lam0=arrays["physical_lambda"]*scale; K=arrays["K"]
    close(np,"saved q0",c["q0_KI1_m"],q0,0.,1e-18)
    close(np,"saved scaled force",c["f_ext_KI1_N"],f,0.,1e-18)
    close(np,"saved scaled gauge",c["g_phys_KI1_m"],g,0.,1e-18)
    close(np,"saved scaled lambda",c["lambda_phys_KI1"],lam0,0.,1e-18)
    C=arrays["physical_constraint_matrix"]
    cached_weak=K@q0+C.T@lam0-f; cached_gauge=C@q0-g
    cached_rel=max(float(np.linalg.norm(cached_weak,np.inf))/max(float(np.linalg.norm(K,np.inf)*np.linalg.norm(q0,np.inf)+np.linalg.norm(C.T@lam0,np.inf)+np.linalg.norm(f,np.inf)),1e-300),
                   float(np.linalg.norm(cached_gauge,np.inf))/max(float(np.linalg.norm(C,np.inf)*np.linalg.norm(q0,np.inf)+np.linalg.norm(g,np.inf)),1e-300))
    require(cached_rel<=tol["kkt_backward_error"],"cached scaled q/lambda do not satisfy saved physical KKT")
    require(len(semantic_ids)==93 and len(source_ids)==27,"cached registry dimensions")
    semantic_index={s:i for i,s in enumerate(semantic_ids)}
    require(semantic_index[json.dumps(["GLOBAL_FRONT","global","uy",1],separators=(",",":"))]==82,"opening semantic index changed")
    require(all(semantic_index[json.dumps(["SOURCE_BASE",sid,comp,0],separators=(",",":"))]==3*ni+ci for ni,sid in enumerate(source_ids) for ci,comp in enumerate(("ux","uy","uz"))),"base semantic order changed")
    require(all(semantic_index[json.dumps(["GLOBAL_FRONT","global",comp,b],separators=(",",":"))]==81+3*(b-1)+ci for b in range(1,5) for ci,comp in enumerate(("ux","uy","uz"))),"front semantic order changed")
    CASE["source_ids"]=source_ids; CASE["semantic_index"]=semantic_index
    rigid=arrays["rigid_vectors"]
    require(rigid.shape==(6,93),"six rigid candidates required")
    center=np.asarray(producer_case["rigid_modes"]["geometric_center_xyz_m"],dtype=float)
    CASE["source_center_xyz_m"]=center.tolist()
    close(np,"cached source geometric center",np.mean(arrays["source_positions"],axis=0),center,0.,1e-14)
    expected_rigid=np.zeros((6,93)); axes=np.eye(3)
    for ni,sid in enumerate(source_ids):
        X=arrays["source_positions"][ni]
        for ci,comp in enumerate(("ux","uy","uz")):
            col=semantic_index[json.dumps(["SOURCE_BASE",sid,comp,0],separators=(",",":"))]
            expected_rigid[ci,col]=1.
            for ri,omega in enumerate(axes):
                expected_rigid[3+ri,col]=np.cross(omega,X-center)[ci]
    close(np,"independently constructed six rigid coefficient vectors",rigid,expected_rigid,0.,1e-15)
    points=r["evidence"]["interface"]; require(len(points)==16 and len({p["quadrature_id"] for p in points})==16,"saved interface/history IDs")
    area_sum=0.; baseline_rows=r["evidence"]["baseline"]["points"]; coupled=r["evidence"]["coupled"]
    final_rows=coupled["points"]
    for qi,p in enumerate(points):
        ix,iz=qi//4,qi%4; x=geo["R_m"]/2*(GL_NODES[ix]-1); z=geo["L_m"]/2*GL_NODES[iz]
        expected_w=geo["R_m"]*geo["L_m"]/4*GL_WEIGHTS[ix]*GL_WEIGHTS[iz]
        close(np,"interface quadrature point",p["point_xyz_m"],[x,0.,z],0.,2e-14)
        close(np,"interface one-side area",p["one_side_area_weight_m2"],expected_w,0.,2e-16)
        area_sum+=p["one_side_area_weight_m2"]
    require(abs(area_sum-geo["one_side_area_m2"])<=2e-16,"one-side area must equal R*L exactly once")
    G=np.zeros((93,93)); G[82,82]=2*CASE["material_fixture"]["K0_Pa_per_m"]*geo["R_m"]*geo["L_m"]
    close(np,"analytic cohesive KG",r["evidence"]["analytic_KG_N_per_m"],G,0.,1e-12)
    require(G[82,82]==20000. and np.count_nonzero(G)==1,"analytic KG closed form")
    mls=r["evidence"]["mls"]; mls_weights=np.asarray(mls["weights_query_source"],dtype=float)
    require(mls_weights.shape==(16,27) and bool(np.isfinite(mls_weights).all()),"saved MLS query weights shape/finiteness")
    close(np,"saved MLS query coordinates",mls["query_xyz_m"],[p["point_xyz_m"] for p in points],0.,2e-14)
    close(np,"saved MLS partition of unity",np.sum(mls_weights,axis=1),np.ones(16),0.,1e-10)
    trace_rows=r["evidence"]["interface"]
    baseline_rows=[dict(trace_rows[i],**baseline_rows[i]) for i in range(16)]
    require(all(baseline_rows[i]["point_id"][-1]==trace_rows[i]["quadrature_id"] for i in range(16)),"baseline history/point IDs do not bind unique QPs")
    base_check=audit_points(np,baseline_rows,q0,rigid,mls_weights,CASE,"baseline",tol)
    qc=np.asarray(r["solution"]["q_c_m"]); lam=np.asarray(r["solution"]["lambda_c"])
    final_rows=[dict(trace_rows[i],**final_rows[i]) for i in range(16)]
    require(all(final_rows[i]["point_id"][-1]==trace_rows[i]["quadrature_id"] for i in range(16)),"final history/point IDs do not bind unique QPs")
    final_check=audit_points(np,final_rows,qc,rigid,mls_weights,CASE,"coupled-final",tol)
    require(abs(base_check["area_m2"]-area_sum)<=2e-16 and abs(final_check["area_m2"]-area_sum)<=2e-16,"baseline/final quadrature double-count or omission")
    close(np,"saved coupled q",coupled["q_c_m"],qc,0.,0.)
    base=r["evidence"]["baseline"]
    close(np,"baseline force",base["generalized_internal_force_N"],G@q0,2e-9,1e-8)
    close(np,"baseline tangent",base["tangent_N_per_m"],G,2e-12,1e-10)
    U0=.5*float(q0@G@q0)
    close(np,"baseline stored ledger",base["stored_energy_j"],U0,2e-9,tol["absolute_energy_J"])
    close(np,"baseline path integral ledger",base["primitive_external_work_path_integral_j"],U0,2e-9,tol["absolute_energy_J"])
    require(base["work_interpretation"]=="not_incremental_work","baseline path integral mislabeled")
    require(base["dissipated_energy_j"]==0. and abs(base_check["stored_energy_J"]-U0)<=tol["absolute_energy_J"],"baseline energy/dissipation")
    close(np,"baseline work pairing",float(q0@np.asarray(base["generalized_internal_force_N"])),
          sum(float(p["one_side_area_weight_m2"])*float(row["opening_m"])*float(row["traction_pa"]) for p,row in zip(points,baseline_rows)),2e-9,tol["absolute_energy_J"])
    close(np,"coupled K+KG",coupled["K_plus_KG_N_per_m"],K+G,2e-12,1e-8)
    close(np,"coupled external force",coupled["f_ext_N"],f,0.,1e-16)
    close(np,"coupled gauge",coupled["g_phys_m"],g,0.,1e-16)
    close(np,"reported baseline q",r["solution"]["q0_m"],q0,0.,1e-18)
    close(np,"coupled multiplier",coupled["lambda_c"],lam,0.,0.)
    A=K+G; weak=A@qc+C.T@lam-f; gauge=C@qc-g
    weak_rel=float(np.linalg.norm(weak,np.inf))/max(float(np.linalg.norm(A,np.inf)*np.linalg.norm(qc,np.inf)+np.linalg.norm(C.T@lam,np.inf)+np.linalg.norm(f,np.inf)),1e-300)
    gauge_rel=float(np.linalg.norm(gauge,np.inf))/max(float(np.linalg.norm(C,np.inf)*np.linalg.norm(qc,np.inf)+np.linalg.norm(g,np.inf)),1e-300)
    ubulk=.5*float(qc@K@qc); ugamma=.5*float(qc@G@qc); pair=float(qc@f); gp=float(g@lam)
    energy_rel=abs(2*(ubulk+ugamma)-pair+gp)/max(abs(2*(ubulk+ugamma)),abs(pair),abs(gp),1e-300)
    require(max(weak_rel,gauge_rel,energy_rel)<=tol["kkt_backward_error"],"independent coupled KKT/energy identity")
    close(np,"reported bulk energy",coupled["U_bulk_j"],ubulk,2e-9,tol["absolute_energy_J"])
    close(np,"reported cohesive energy",coupled["U_gamma_j"],ugamma,2e-9,tol["absolute_energy_J"])
    close(np,"reported force pairing",coupled["external_force_displacement_pairing_j"],pair,2e-9,tol["absolute_energy_J"])
    close(np,"reported gauge pairing",coupled["gauge_multiplier_pairing_j"],gp,2e-9,tol["absolute_energy_J"])
    close(np,"coupled internal force",coupled["generalized_internal_force_N"],G@qc,2e-9,1e-8)
    close(np,"coupled tangent",coupled["tangent_N_per_m"],G,2e-12,1e-10)
    close(np,"coupled primitive path integral",coupled["primitive_external_work_path_integral_j"],ugamma,2e-9,tol["absolute_energy_J"])
    require("not incremental" in coupled["energy_identity_interpretation"] and "not global" in coupled["primitive_external_work_interpretation"],"coupled energy work interpretation changed")
    close(np,"top-level bulk energy",r["energy"]["U_bulk_j"],ubulk,2e-9,tol["absolute_energy_J"])
    close(np,"top-level cohesive energy",r["energy"]["U_gamma_j"],ugamma,2e-9,tol["absolute_energy_J"])
    close(np,"top-level force pairing",r["energy"]["external_force_displacement_pairing_j"],pair,2e-9,tol["absolute_energy_J"])
    close(np,"top-level gauge pairing",r["energy"]["gauge_multiplier_pairing_j"],gp,2e-9,tol["absolute_energy_J"])
    require(r["energy"]["identity_is_not_incremental_work_or_loading_path_assessment"] is True,"energy interpretation must remain algebraic")
    require(abs(final_check["stored_energy_J"]-ugamma)<=tol["absolute_energy_J"] and final_check["dissipated_energy_J"]==0.,"final virgin cohesive energy ledger")
    dq=float(np.linalg.norm(qc-q0)); fn=float(np.linalg.norm(G@qc))
    require(dq>128*np.finfo(float).eps*float(np.linalg.norm(q0)) and fn>128*np.finfo(float).eps*float(np.linalg.norm(f)),"coupling did not change q/force")
    close(np,"reported q change norm",r["solution"]["q_change_norm_m"],dq,2e-9,1e-18)
    close(np,"reported internal force norm",r["solution"]["internal_force_norm_n"],fn,2e-9,1e-10)
    zgate=r["gates"]["zero_tangent_control"]["evidence"]
    require(zgate["synthetic_fully_damaged_history"] is True and zgate["state_discarded"] is True and zgate["physical_solve_input"] is False
            and zgate["force_norm_N"]==0. and zgate["tangent_norm_N_per_m"]==0.,"zero-tangent control provenance")
    return {"schema":"mixed93_actual_cohesive_saved_independent_result.v1","status":"PASS_DIAGNOSTIC_ONLY",
            "classification":r["classification"],"source_hashes":source_hashes,"result_sha256":result_sha,
            "receipt_sha256":receipt_sha,"npz_sha256":CASE["expected_hashes"]["producer_npz"],
            "gate_status":{"producer":{n:r["gates"][n]["status"] for n in CASE["required_gates"]},"independent":{"artifact_binding":"PASS","cached_kkt":"PASS","interface_area_and_traces":"PASS","six_rigid_modes":"PASS","baseline_traction_force_energy":"PASS","synthetic_zero_control":"PASS","coupled_kkt_and_energy":"PASS","mechanical_change":"PASS"}},
            "independent_checks":{"all_producer_gates_PASS":True,"cache_bytes_match_npz":True,
            "cache_q_lambda_kkt_relative_error":cached_rel,"quadrature_points":16,"one_side_area_m2":area_sum,
            "rigid_trace_checks_baseline_and_final":True,"analytic_KG_82_82_N_per_m":float(G[82,82]),
            "baseline":base_check,"final":final_check,"coupled_weak_relative_error":weak_rel,
            "coupled_gauge_relative_error":gauge_rel,"coupled_energy_identity_relative_error":energy_rel,
            "q_change_norm_m":float(np.linalg.norm(qc-q0)),"zero_tangent_control_is_synthetic":True,
            "spectrum_recomputed":False,"producer_spectrum_gates_only":True},
            "limitations":["fixed crack; no growth","no incremental external-work claim","toy material fixture only; no calibration or qualification","zero-tangent control is synthetic","no solver/eigensolve/provider/author module call"]}
def main() -> int:
    p=argparse.ArgumentParser(description="one root-authorized saved-result audit; never runs the producer")
    p.add_argument("--execute",action="store_true"); p.add_argument("--root-go"); p.add_argument("--output",type=Path)
    a=p.parse_args()
    if not a.execute or a.output is None: print(json.dumps({"status":"EXECUTION_DISABLED","npz_opened":False,"result_read":False})); return 2
    try:
        out=run(ROOT/PINS["result"]["path"],ROOT/PINS["receipt"]["path"],a.root_go or "")
        a.output.parent.mkdir(parents=True,exist_ok=True)
        with a.output.open("x",encoding="utf-8",newline="\n") as f: f.write(json.dumps(out,indent=2,sort_keys=True,allow_nan=False)+"\n")
        print(json.dumps({"status":out["status"],"output":str(a.output)})); return 0
    except Exception as exc:
        print(json.dumps({"status":"FAIL_CLOSED","error":f"{type(exc).__name__}: {exc}"})); return 3
if __name__=="__main__": raise SystemExit(main())
