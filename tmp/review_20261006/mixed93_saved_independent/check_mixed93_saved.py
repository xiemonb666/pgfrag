"""Independent mixed93 cache reader: --preflight is stdlib metadata only;
--audit is the separately authorized saved-array readback. Never imports the
producer, queries a field/provider, creates samples, or solves a system."""
from __future__ import annotations
import argparse, ast, hashlib, json, math, time, zipfile
from pathlib import Path, PurePosixPath

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
RUN = ROOT / "pgfrag/tmp/review_20261006/mixed93_cached_lefm_solve_20261006"
FROZEN = RUN / "failure_source_freeze"
CAPTURE = RUN / "execution_capture"
ARRAYS = RUN / "run01/arrays.npz"
FAILED_RESULT = RUN / "run01/result.json"
RESULT = ROOT / "pgfrag/tmp/review_20261006/tip_aligned_quadrature_study_20261006/run01/result.json"
UPSTREAM = RESULT.with_name("arrays.npz")
OUT = HERE / "audit.json"
IDS = ("pure_I", "pure_II", "pure_III", "mixed_70_40_100")
COMP = ("ux", "uy", "uz")
EPS = 2.220446049250313e-16
PINS = {
    "frozen_entry":"61e88d32a6d5543c891e92c8996b66f1562db6a20781db061c6f266d0dffbfa8", "case":"67ff4df51a0269e88b8d107a2a95929c3990a7cd9605585a31a2df583434cc1c",
    "source_pins":"c168fc96691b95ef6a4cfa54c1d17a6898abce08d5100d72cdf4e8ec0ddb0df4", "execution_receipt":"32bfd3dcedcf2275768c482cfd936f12482756ef3e97ebe73176167d25e4d2a1",
    "stdout":"e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855", "stderr":"d9b1de8626c37b363db57db5ef6ab2f5e585efbf122be36be2a46b2c75a39ecf",
    "failed_arrays":"3e1e6fd17751e4e06113010f20ea33c7b69914ca90431e26bbf76124951b5d41", "upstream_result":"2510eabd3e3406f243d9fe900ac5af3cfdf17071cdbe50383093de0a13adff4b",
    "upstream_arrays":"a9de2e3eaf88ba0eb130fca6f160d5d327a963bb00203266b7390cadbe6693ea", "upstream_driver":"a719854abbbbb06d881fb4fe7cc80e75e878109a1db001302d96e02045358ee9",
    "field_source":"ed90063029d69649fbb8b4a0785946cfdd7e7b6a22e483513ff27f447745dacd", "enrichment_source":"466d0c01445f497e46a49314a03e89e92c62a204a9665cf01c3d2005f2dc99e3",
    "semantic_source":"f28f3b55e1ee15212c72ffba75bfe5a4b4f324a7fd4921d63dd47cbd840d7e6f",
}
SOURCE_PATHS = {
    "field_source": ROOT / "pgfrag/experiments/mixed_lefm_benchmark_20261006/mixed_lefm_field.py",
    "enrichment_source": ROOT / "pgfrag/experiments/mixed_lefm_enrichment_algebra_20261006/mixed_enrichment_algebra.py", "semantic_source": ROOT / "pgfrag/experiments/mixed_enriched_operator_20261006/mixed_enriched_operator.py",
}

def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""): h.update(chunk)
    return h.hexdigest()

def read_npy_headers(path):
    """Return member name/shape/dtype from headers; do not read NPY payloads."""
    out = {}
    with zipfile.ZipFile(path) as zf:
        for info in zf.infolist():
            if not info.filename.endswith(".npy"): raise ValueError("non-NPY archive member")
            with zf.open(info) as f:
                if f.read(6) != b"\x93NUMPY": raise ValueError("bad NPY magic")
                version = tuple(f.read(2))
                if version not in ((1, 0), (2, 0), (3, 0)): raise ValueError("unsupported NPY header")
                n = int.from_bytes(f.read(2 if version == (1, 0) else 4), "little")
                enc = "utf-8" if version == (3, 0) else "latin1"
                header = ast.literal_eval(f.read(n).decode(enc).strip())
            key = PurePosixPath(info.filename).stem
            if key in out: raise ValueError("duplicate NPZ key: " + key)
            if header.get("fortran_order"): raise ValueError("unexpected Fortran-order array: " + key)
            out[key] = {"shape": tuple(header["shape"]), "dtype": header["descr"], "nbytes": info.file_size}
    return out

def verify_inputs():
    paths = {
        "frozen_entry": FROZEN / "mixed93_cached_lefm_solve.py",
        "case": FROZEN / "case.authorized.json",
        "source_pins": FROZEN / "SOURCE_PINS.authorized.json",
        "execution_receipt": CAPTURE / "execution_receipt.json",
        "stdout": CAPTURE / "stdout.bin", "stderr": CAPTURE / "stderr.bin",
        "failed_arrays": ARRAYS, "upstream_result": RESULT, "upstream_arrays": UPSTREAM,
        "upstream_driver": ROOT / "pgfrag/experiments/tip_aligned_quadrature_study_20261006/tip_aligned_quadrature_study.py",
        **SOURCE_PATHS,
    }
    got = {k: sha(v) for k, v in paths.items()}
    for key, expected in PINS.items():
        if got.get(key) != expected: raise ValueError(f"SHA pin mismatch {key}: {got.get(key)}")
    case = json.loads(paths["case"].read_text(encoding="utf-8"))
    pins = json.loads(paths["source_pins"].read_text(encoding="utf-8"))
    receipt = json.loads(paths["execution_receipt"].read_text(encoding="utf-8"))
    upstream = json.loads(RESULT.read_text(encoding="utf-8"))
    if case.get("case_id") != "mixed93_cached_lefm_solve_20261006": raise ValueError("wrong frozen case")
    if case.get("execution_authorization", {}).get("authorized") is not True: raise ValueError("missing frozen authorization")
    if case.get("load_cases", []) and [x["case_id"] for x in case["load_cases"]] != list(IDS): raise ValueError("four load-case order changed")
    if pins.get("entry_sha256") != got["frozen_entry"] or pins.get("current_case_sha256") != got["case"]: raise ValueError("frozen manifest/source mismatch")
    if pins.get("references", {}).get("cached_result", {}).get("sha256") != got["upstream_result"]: raise ValueError("manifest upstream JSON mismatch")
    if pins.get("references", {}).get("cached_arrays", {}).get("sha256") != got["upstream_arrays"]: raise ValueError("manifest upstream NPZ mismatch")
    if upstream.get("schema") != "tip_aligned_full3d_quadrature_study_result.v1" or upstream.get("status") != "OK_DIAGNOSTIC_ONLY_UNQUALIFIED": raise ValueError("unexpected upstream cache status")
    if receipt.get("schema") != "mixed93_one_shot_execution_capture.v1" or receipt.get("actual_process_return_code") != 1 or receipt.get("timed_out") is not False: raise ValueError("execution receipt is not the failed one-shot record")
    if receipt.get("elapsed_seconds") != 0.6208375999704003 or receipt.get("input_hashes_unchanged_during_run") is not True: raise ValueError("failed-run timing/pin evidence changed")
    if receipt.get("pre_run_eight_pin_hashes") != receipt.get("post_run_eight_pin_hashes"): raise ValueError("producer input pins changed during failed run")
    if receipt.get("output_files", {}).get("result_json", {}).get("exists") is not False: raise ValueError("unexpected producer JSON result exists")
    if FAILED_RESULT.exists(): raise ValueError("failed producer result JSON unexpectedly exists")
    produced = receipt.get("output_files", {}).get("arrays_npz", {})
    if (produced.get("sha256"), produced.get("bytes")) != (got["failed_arrays"], 1569232): raise ValueError("failed-run array output pin mismatch")
    linked={"entry":"frozen_entry","case":"case","source_pins_manifest":"source_pins","cached_arrays":"upstream_arrays","cached_result":"upstream_result","global_semantic_operator":"semantic_source","mixed_enrichment_algebra":"enrichment_source","mixed_lefm_field":"field_source"}
    if any(receipt["pre_run_eight_pin_hashes"].get(k)!=got[v] for k,v in linked.items()): raise ValueError("receipt does not bind frozen source/case/dependencies")
    if upstream.get("source_sha256", {}).get("driver") != got["upstream_driver"]: raise ValueError("upstream JSON does not bind pinned study driver")
    source_ids = upstream.get("registry", {}).get("source_ids", [])
    if len(source_ids) != 27 or source_ids != sorted(source_ids) or len(set(source_ids)) != 27: raise ValueError("upstream source registry is not sorted 27-source registry")
    if upstream.get("registry", {}).get("component_order") != list(COMP): raise ValueError("upstream component order changed")
    if upstream.get("registry", {}).get("basis_order") != ["base", "front_f1_unshifted"]: raise ValueError("upstream cached base basis order changed")
    pre = read_npy_headers(ARRAYS); up = read_npy_headers(UPSTREAM)
    require_shapes(pre, up)
    if len(pre) != 146 or len(up) != 179: raise ValueError("NPZ member count differs from frozen schema")
    return paths, got, case, pins, receipt, upstream, pre, up

def require_shapes(arr, upstream):
    shapes = {"K_N_per_m":(93,93), "K_diagonally_scaled":(93,93), "constitutive_matrix_Pa":(6,6), "source_positions_xyz_m":(27,3), "source_ids":(27,), "semantic_ids":(93,),
        "rigid_coefficient_modes":(6,93), "physical_volume_moment_matrix":(6,93), "physical_gauge_rigid_coupling":(6,6), "volume_xyz_m":(512,3), "volume_unique_row_ids":(512,),
        "volume_weights_m3":(512,), "volume_basis_values":(512,31), "volume_basis_gradients_per_m":(512,31,3), "boundary_xyz_m":(768,3), "boundary_unique_row_ids":(768,),
        "boundary_weights_m2":(768,), "boundary_normals":(768,3), "boundary_basis_values":(768,31), "boundary_basis_gradients_per_m":(768,31,3),
        "offgrid_xyz_m":(128,3), "offgrid_unique_row_ids":(128,), "offgrid_basis_values":(128,31), "offgrid_basis_gradients_per_m":(128,31,3),
        "Cbar_scaled":(6,93), "Chat_scaled":(6,93), "scaled_KKT_matrix":(99,99)}
    for key, shape in shapes.items():
        if key not in arr or arr[key]["shape"] != shape: raise ValueError(f"cached NPZ header mismatch {key}")
    cached_shapes={"source_positions_xyz_m":(27,3), "unique_query_positions_xyz_m":(4224,3), "base_mls_values":(4224,27), "base_mls_gradients_per_m":(4224,27,3),
        "cloud_a0_volume32_xyz_m":(512,3), "cloud_a0_volume32_weights":(512,), "cloud_a0_volume32_unique_row_ids":(512,), "cloud_side_xyz_m":(256,3), "cloud_side_weights":(256,),
        "cloud_side_normals":(256,3), "cloud_side_unique_row_ids":(256,), "cloud_a0_cap16_xyz_m":(512,3), "cloud_a0_cap16_weights":(512,),
        "cloud_a0_cap16_normals":(512,3), "cloud_a0_cap16_unique_row_ids":(512,), "cloud_offgrid_xyz_m":(128,3), "cloud_offgrid_unique_row_ids":(128,)}
    for key, shape in cached_shapes.items():
        if key not in upstream or upstream[key]["shape"] != shape: raise ValueError(f"upstream NPZ header mismatch {key}")
    case_shapes=(("q_exact",(93,)),("q_solution",(93,)),("lambda_physical",(6,)),("lambda_scaled",(6,)),("kkt_rhs_scaled",(99,)),
        ("kkt_solution_scaled",(99,)),("gauge_rhs_scaled",(6,)),("physical_gauge_target_m",(6,)),("analytic_volume_moments_m",(6,)),("boundary_force_N",(93,)),("boundary_traction_Pa",(768,3)))
    for case_id in IDS:
        for key, shape in case_shapes:
            if arr.get(f"{case_id}_{key}",{}).get("shape") != shape: raise ValueError(f"missing/bad case array {case_id}_{key}")
        for group,n in (("volume",512),("boundary",768),("offgrid",128)):
            for kind,shape in (("exact_u",(n,3)),("exact_strain",(n,6)),("exact_stress",(n,6)),("solution_u_m",(n,3)),("solution_strain_m_per_m",(n,6)),("solution_stress_Pa",(n,6))):
                key=f"{case_id}_{group}_{kind}"
                if arr.get(key,{}).get("shape") != shape: raise ValueError(f"missing/bad case field {key}")

def preflight():
    _, got, _, _, receipt, _, arrays, cached = verify_inputs()
    return {"status":"SOURCE_AND_CACHED_SCHEMA_READY_NO_NUMERICAL_ARRAYS_READ","root":str(ROOT),"pins":{k:got[k] for k in ("frozen_entry","case","execution_receipt","failed_arrays","upstream_result","upstream_arrays")},
        "npz":{"failed_members":len(arrays),"upstream_members":len(cached)},"producer":{"return_code":receipt["actual_process_return_code"],"elapsed_seconds":receipt["elapsed_seconds"]},"case_ids":list(IDS),"dof_count":93,"kkt_order":99,"source_count":27,"global_branch_count":4,"phase":"no NumPy import; NPY payloads untouched"}

def slot(s, c): return 3*s+c if s < 27 else 81+3*(s-27)+c

def voigt_B(np, gradients):
    n=gradients.shape[0]; b=np.zeros((n,6,93),dtype=np.float64)
    for s in range(31):
        gx,gy,gz=(gradients[:,s,k] for k in range(3))
        for c in range(3):
            j=slot(s,c)
            if c==0: b[:,0,j],b[:,3,j],b[:,5,j]=gx,gy,gz
            elif c==1: b[:,1,j],b[:,3,j],b[:,4,j]=gy,gx,gz
            else: b[:,2,j],b[:,4,j],b[:,5,j]=gz,gy,gx
    return b

def displacement(np, values, q):
    u=np.zeros((len(values),3),dtype=np.float64)
    for s in range(31):
        for c in range(3): u[:,c] += values[:,s]*q[slot(s,c)]
    return u

def stiffness(np, gradients, weights, D):
    K=np.zeros((93,93),dtype=np.float64)
    for start in range(0,len(weights),256):
        stop=min(start+256,len(weights)); B=voigt_B(np,gradients[start:stop])
        K += np.einsum("qai,ab,qbj,q->ij",B,D,B,weights[start:stop],optimize=True)
    return K

def tensor6(np, v):
    t=np.zeros((len(v),3,3),dtype=np.float64)
    t[:,0,0],t[:,1,1],t[:,2,2]=v[:,0],v[:,1],v[:,2]
    t[:,0,1]=t[:,1,0]=v[:,3]; t[:,1,2]=t[:,2,1]=v[:,4]; t[:,0,2]=t[:,2,0]=v[:,5]; return t

def force_from_saved(np, values, weights, normals, stress):
    tr=np.einsum("qij,qj->qi",tensor6(np,stress),normals)
    f=np.zeros(93,dtype=np.float64)
    for s in range(31):
        for c in range(3): f[slot(s,c)]=np.dot(weights,values[:,s]*tr[:,c])
    return f,tr

def physical_C(np, values, xyz, weights, center, radius):
    volume=float(np.sum(weights)); C=np.zeros((6,93),dtype=np.float64); rel=xyz-np.asarray(center)
    for s in range(31):
        pw=weights*values[:,s]/volume
        for c in range(3):
            j=slot(s,c); C[c,j]=float(np.sum(pw))
            if c==0: C[4,j]=float(np.dot(pw,rel[:,2]/radius)); C[5,j]=float(np.dot(pw,-rel[:,1]/radius))
            elif c==1: C[3,j]=float(np.dot(pw,-rel[:,2]/radius)); C[5,j]=float(np.dot(pw,rel[:,0]/radius))
            else: C[3,j]=float(np.dot(pw,rel[:,1]/radius)); C[4,j]=float(np.dot(pw,-rel[:,0]/radius))
    return C,volume

def exact_q(np, ids, material, load, r):
    e=float(material["E_Pa"]); nu=float(material["nu"]); r=float(r)
    mu=e/(2*(1+nu)); kap=3-4*nu; root=math.sqrt(r/(2*math.pi))/(2*mu)
    ai,aii=float(load[0])*root,float(load[1])*root; aiii=float(load[2])/mu*math.sqrt(2*r/math.pi)
    branches={1:(aii*(kap+1),ai*(kap+1),aiii),2:(ai*(kap-1),aii*(1-kap),0.),3:(ai,aii,0.),4:(aii,-ai,0.)}
    q=np.zeros(93,dtype=np.float64)
    for j,token in enumerate(ids):
        scope,owner,component,index=json.loads(str(token))
        if scope=="GLOBAL_FRONT": q[j]=branches[index][{"ux":0,"uy":1,"uz":2}[component]]
        elif scope!="SOURCE_BASE" or index!=0: raise ValueError("unknown semantic token")
    return q

def relative(np, a, b): return float(np.linalg.norm(a-b)/max(float(np.linalg.norm(b)),np.finfo(float).tiny))

def weighted_relative(np, a, b, weights, cw=None):
    d=(a-b)**2; r=b**2
    if cw is not None: d=d*cw; r=r*cw
    return float(math.sqrt(max(0.,float(np.sum(weights[:,None]*d))/max(float(np.sum(weights[:,None]*r)),np.finfo(float).tiny))))

def expected_tokens(source_ids): return [json.dumps(["SOURCE_BASE",s,c,0],separators=(",",":")) for s in source_ids for c in COMP]+[json.dumps(["GLOBAL_FRONT","global",c,b],separators=(",",":")) for b in range(1,5) for c in COMP]

def run_audit():
    t0=time.monotonic(); own_pre=sha(Path(__file__)); paths, pre, case, manifest, receipt, upstream, _, _=verify_inputs()
    if OUT.exists(): raise FileExistsError(f"refusing to overwrite {OUT}")
    import numpy as np
    material=dict(case["material"]); geom=case["geometry"]
    radius=float(geom["radius_m"]); length=float(geom["length_m"]); center=np.asarray(case["rigid_gauge"]["center_xyz_m"],dtype=np.float64)
    e=float(material["E_Pa"]); nu=float(material["nu"]); mu=e/(2*(1+nu)); lam=e*nu/((1+nu)*(1-2*nu))
    D=np.zeros((6,6),dtype=np.float64)
    for i in range(3):
        for j in range(3): D[i,j]=lam+(2*mu if i==j else 0.)
    D[3,3]=D[4,4]=D[5,5]=mu
    out={}; cache={}
    with np.load(ARRAYS,allow_pickle=False) as z, np.load(UPSTREAM,allow_pickle=False) as c:
        source_ids=list(upstream["registry"]["source_ids"])
        if list(map(str,z["source_ids"].tolist()))!=source_ids: raise ValueError("saved source ID order differs from frozen upstream registry")
        ids=list(map(str,z["semantic_ids"].tolist()))
        if ids!=expected_tokens(source_ids): raise ValueError("93-DOF semantic tokens/layout differ from base-then-global-branch contract")
        if not np.array_equal(z["source_positions_xyz_m"],c["source_positions_xyz_m"]): raise ValueError("source positions differ from pinned MLS cache")
        region_data={
            "volume":(z["volume_xyz_m"],z["volume_unique_row_ids"],z["volume_weights_m3"],None,z["volume_basis_values"],z["volume_basis_gradients_per_m"],c["cloud_a0_volume32_xyz_m"],c["cloud_a0_volume32_unique_row_ids"],c["cloud_a0_volume32_weights"],None),
            "boundary":(z["boundary_xyz_m"],z["boundary_unique_row_ids"],z["boundary_weights_m2"],z["boundary_normals"],z["boundary_basis_values"],z["boundary_basis_gradients_per_m"],np.concatenate((c["cloud_side_xyz_m"],c["cloud_a0_cap16_xyz_m"])),np.concatenate((c["cloud_side_unique_row_ids"],c["cloud_a0_cap16_unique_row_ids"])),np.concatenate((c["cloud_side_weights"],c["cloud_a0_cap16_weights"])),np.concatenate((c["cloud_side_normals"],c["cloud_a0_cap16_normals"]))),
            "offgrid":(z["offgrid_xyz_m"],z["offgrid_unique_row_ids"],np.ones(128),None,z["offgrid_basis_values"],z["offgrid_basis_gradients_per_m"],c["cloud_offgrid_xyz_m"],c["cloud_offgrid_unique_row_ids"],np.ones(128),None),
        }
        for name,(xyz,rid,w,norm,val,grad,refxyz,refid,refw,refnorm) in region_data.items():
            if not np.array_equal(rid,refid) or not np.array_equal(xyz,refxyz): raise ValueError(f"cached {name} cloud/row token mismatch")
            if not np.array_equal(w,refw) or (norm is not None and not np.array_equal(norm,refnorm)): raise ValueError(f"cached {name} weights/normals mismatch")
            if not np.array_equal(val[:,:27],c["base_mls_values"][rid]) or not np.array_equal(grad[:,:27,:],c["base_mls_gradients_per_m"][rid]): raise ValueError(f"cached {name} base MLS rows/order mismatch")
            if not np.isfinite(xyz).all() or not np.isfinite(w).all() or np.any(w<=0): raise ValueError(f"invalid {name} coordinates/weights")
            if norm is not None and float(np.max(np.abs(np.linalg.norm(norm,axis=1)-1)))>1e-10: raise ValueError(f"non-unit {name} normals")
            cache[name]={"xyz":xyz,"weights":w,"normals":norm,"values":val,"gradients":grad,"B":voigt_B(np,grad)}
        K=stiffness(np,cache["volume"]["gradients"],cache["volume"]["weights"],D)
        Ksaved=z["K_N_per_m"]; s=1/np.sqrt(np.diag(K)); Ks=s[:,None]*K*s[None,:]
        evals=np.linalg.eigvalsh(Ks); eigscale=float(np.max(np.abs(evals))); eigzero=64*EPS*93*eigscale
        eigcounts={"zero":int(np.count_nonzero(np.abs(evals)<=eigzero)),"positive":int(np.count_nonzero(evals>eigzero)),"negative":int(np.count_nonzero(evals < -eigzero)),"tol":float(eigzero)}
        C,V=physical_C(np,cache["volume"]["values"],cache["volume"]["xyz"],cache["volume"]["weights"],center,radius)
        Csaved=z["physical_volume_moment_matrix"]
        Z=np.zeros((93,6),dtype=np.float64); src=z["source_positions_xyz_m"]-np.mean(z["source_positions_xyz_m"],axis=0)
        for i,p in enumerate(src):
            for cidx in range(3): Z[3*i+cidx,cidx]=1.
            for axis in range(3):
                omega=np.eye(3)[axis]; cr=np.cross(omega,p)
                for cidx in range(3): Z[3*i+cidx,3+axis]=cr[cidx]
        rigid_null=float(np.linalg.norm(Ks@(Z/s[:,None]),ord=np.inf)/max(np.linalg.norm(Ks,ord=np.inf)*np.linalg.norm(Z/s[:,None],ord=np.inf),np.finfo(float).tiny))
        coupling=C@Z; sv=np.linalg.svd(coupling,compute_uv=False); rank=int(np.count_nonzero(sv>128*EPS*93*float(sv[0])))
        kdiag={"D_saved_rel":relative(np,D,z["constitutive_matrix_Pa"]),"relative_rebuild":relative(np,K,Ksaved),"raw_asymmetry":float(np.linalg.norm(K-K.T)/max(np.linalg.norm(K),np.finfo(float).tiny)),"scaled_asymmetry":float(np.linalg.norm(Ks-Ks.T)/max(np.linalg.norm(Ks),np.finfo(float).tiny)),"scaled_saved_rel":relative(np,Ks,z["K_diagonally_scaled"]),"spectrum_counts":eigcounts,"saved_eigenvalues_rel":relative(np,evals,z["eigenvalues_scaled"]),"volume_m3":V,"analytic_cylinder_volume_m3":math.pi*radius**2*length,"rigid_saved_rel":relative(np,Z.T,z["rigid_coefficient_modes"]),"rigid_scaled_null_rel":rigid_null,"physical_C_saved_rel":relative(np,C,Csaved),"coupling_saved_rel":relative(np,coupling,z["physical_gauge_rigid_coupling"]),"physical_rigid_coupling_rank":rank}
        b=cache["boundary"]; gmat=C; row_norm=np.linalg.norm(gmat*s[None,:],axis=1); row_scale=1/row_norm; Cbar=gmat*s[None,:]; Chat=row_scale[:,None]*Cbar
        KKT=np.zeros((99,99)); KKT[:93,:93]=Ks; KKT[:93,93:]=Chat.T; KKT[93:,:93]=Chat
        case_metrics={}; gates={"matrix_and_layout":True,"spectrum_and_rigid_modes":True,"four_case_algebra_and_fields":True}
        case_loads={x["case_id"]:x["K_Pa_sqrt_m"] for x in case["load_cases"]}
        for cid in IDS:
            qe=exact_q(np,ids,material,case_loads[cid],radius); qes=z[cid+"_q_exact"]
            if relative(np,qe,qes)>1e-11: gates["matrix_and_layout"]=False
            g=C@qe; gu=z[cid+"_physical_gauge_target_m"]
            analytic=z[cid+"_analytic_volume_moments_m"]
            u_exact_vol=z[cid+"_volume_exact_u"]
            xyz=cache["volume"]["xyz"]-center; wv=cache["volume"]["weights"]; direct=np.concatenate((np.einsum("q,qi->i",wv,u_exact_vol)/V,np.einsum("q,qi->i",wv,np.cross(xyz,u_exact_vol))/(V*radius)))
            qg_rel=relative(np,g,gu); moment_rel=relative(np,direct,analytic); target_rel=relative(np,g,direct)
            f,tr=force_from_saved(np,b["values"],b["weights"],b["normals"],z[cid+"_boundary_exact_stress"])
            fsave=z[cid+"_boundary_force_N"]; trsave=z[cid+"_boundary_traction_Pa"]
            f_rel=relative(np,f,fsave); tr_rel=relative(np,tr,trsave)
            field_metrics={}; exact_gate=True; sol_gate=True; eng={}
            for grp,data in cache.items():
                val=data["values"]; B=data["B"]; weights=data["weights"]
                strain_e=np.einsum("qai,i->qa",B,qe); stress_e=strain_e@D.T; u_e=displacement(np,val,qe)
                refs=(z[f"{cid}_{grp}_exact_u"],z[f"{cid}_{grp}_exact_strain"],z[f"{cid}_{grp}_exact_stress"])
                rec=(u_e,strain_e,stress_e); exrel=[relative(np,a,r) for a,r in zip(rec,refs)]
                abs_floor=[]
                for a,r in zip(rec,refs):
                    tol=1e-11*float(np.linalg.norm(r))+64*EPS*float(np.max(np.abs(r)))*math.sqrt(max(1,np.size(r)))
                    abs_floor.append(float(np.linalg.norm(a-r))<=tol)
                exact_gate=exact_gate and all(abs_floor)
                q=z[cid+"_q_solution"]
                strain_s=np.einsum("qai,i->qa",B,q); stress_s=strain_s@D.T; u_s=displacement(np,val,q)
                refs_s=(z[f"{cid}_{grp}_solution_u_m"],z[f"{cid}_{grp}_solution_strain_m_per_m"],z[f"{cid}_{grp}_solution_stress_Pa"])
                rec_s=(u_s,strain_s,stress_s); saved_rel=[relative(np,a,r) for a,r in zip(rec_s,refs_s)]
                if max(saved_rel)>1e-10: sol_gate=False
                wfield=weights if grp!="offgrid" else np.ones(len(val))
                wc=(None,np.asarray((1.,1.,1.,.5,.5,.5)),np.asarray((1.,1.,1.,2.,2.,2.)))
                err=[weighted_relative(np,a,r,wfield,cw) for a,r,cw in zip(rec_s,refs,wc)]
                eng[grp]={"displacement_rel_L2":err[0],"strain_tensor_rel_L2":err[1],"stress_tensor_rel_L2":err[2]}
                field_metrics[grp]={"q_exact_reconstruction_rel_u_strain_stress":exrel,"q_solution_vs_saved_rel_u_strain_stress":saved_rel}
            Uhalf=float(.5*q@(K@q)); qf=float(q@f); glam=float(g@z[cid+"_lambda_physical"]); identity=2*Uhalf-qf+glam
            e_rel=abs(identity)/max(abs(2*Uhalf),abs(qf),abs(glam),np.finfo(float).tiny)
            weak=K@q+C.T@z[cid+"_lambda_physical"]-f
            weakrel=float(np.linalg.norm(weak,ord=np.inf)/max(np.linalg.norm(K,ord=np.inf)*np.linalg.norm(q,ord=np.inf)+np.linalg.norm(C.T@z[cid+"_lambda_physical"],ord=np.inf)+np.linalg.norm(f,ord=np.inf),np.finfo(float).tiny))
            gres=C@q-g; g_rel=float(np.linalg.norm(gres,ord=np.inf)/max(np.linalg.norm(C,ord=np.inf)*np.linalg.norm(q,ord=np.inf)+np.linalg.norm(g,ord=np.inf),np.finfo(float).tiny))
            rhs=np.concatenate((s*f,row_scale*g)); sol=z[cid+"_kkt_solution_scaled"]; kkt_res=KKT@sol-rhs
            kktrel=float(np.linalg.norm(kkt_res,ord=np.inf)/max(np.linalg.norm(KKT,ord=np.inf)*np.linalg.norm(sol,ord=np.inf)+np.linalg.norm(rhs,ord=np.inf),np.finfo(float).tiny))
            restore_q=s*sol[:93]; restore_l=row_scale*sol[93:]
            kkt={"backward_rel":kktrel,"saved_rhs_rel":relative(np,rhs,z[cid+"_kkt_rhs_scaled"]),"saved_gauge_rhs_rel":relative(np,row_scale*g,z[cid+"_gauge_rhs_scaled"]),"restored_q_rel":relative(np,restore_q,q),"restored_lambda_rel":relative(np,restore_l,z[cid+"_lambda_physical"]),"saved_KKT_matrix_rel":relative(np,KKT,z["scaled_KKT_matrix"]),"saved_Cbar_rel":relative(np,Cbar,z["Cbar_scaled"]),"saved_Chat_rel":relative(np,Chat,z["Chat_scaled"]),"saved_row_scale_rel":relative(np,row_scale,z["physical_gauge_row_scale"])}
            kkt_saved_ok=max(v for k,v in kkt.items() if k!="backward_rel")<=1e-10
            declared=[eng["volume"]["displacement_rel_L2"],eng["offgrid"]["displacement_rel_L2"],eng["volume"]["strain_tensor_rel_L2"],eng["volume"]["stress_tensor_rel_L2"],eng["offgrid"]["stress_tensor_rel_L2"]]
            field_ok=exact_gate and sol_gate and max(declared)<=.02
            case_ok=(qg_rel<=1e-11 and moment_rel<=1e-11 and target_rel<=1e-11 and f_rel<=1e-11 and tr_rel<=1e-11 and kktrel<=1e-9 and kkt_saved_ok and weakrel<=1e-9 and g_rel<=1e-9 and e_rel<=1e-9 and field_ok)
            gates["four_case_algebra_and_fields"] &= bool(case_ok)
            case_metrics[cid]={"q_exact_formula_rel":relative(np,qe,qes),"Cq_exact_vs_saved_g_rel":qg_rel,"analytic_saved_u_moment_vs_saved_g_rel":moment_rel,"Cq_vs_analytic_moment_rel":target_rel,"g_norm_m":float(np.linalg.norm(g)),"g_nonzero":bool(np.linalg.norm(g)>0),"f_rebuild_rel":f_rel,"traction_rebuild_rel":tr_rel,"KKT":kkt,"KKT_saved_algebra_gate":bool(kkt_saved_ok),"physical_weak_rel":weakrel,"physical_gauge_rel":g_rel,"U_J":Uhalf,"q_dot_f_J":qf,"g_dot_lambda_J":glam,"energy_identity_defect_J":identity,"energy_identity_rel":e_rel,"engineering_field_rel_L2":eng,"field_reconstruction_checks":field_metrics,"exact_representation_gate":bool(exact_gate),"saved_solution_field_reconstruction_gate":bool(sol_gate),"declared_engineering_field_gate":bool(max(declared)<=.02),"declared_engineering_max_rel":float(max(declared))}
        gates["matrix_and_layout"] &= bool(kdiag["D_saved_rel"]<=1e-12 and kdiag["relative_rebuild"]<=1e-10 and kdiag["scaled_saved_rel"]<=1e-10 and kdiag["physical_C_saved_rel"]<=1e-10 and kdiag["rigid_saved_rel"]<=1e-12 and kdiag["coupling_saved_rel"]<=1e-10 and kdiag["saved_eigenvalues_rel"]<=1e-10 and max(kdiag["raw_asymmetry"],kdiag["scaled_asymmetry"])<=1e-12)
        gates["spectrum_and_rigid_modes"] &= bool(kdiag["spectrum_counts"]["zero"]==6 and kdiag["spectrum_counts"]["positive"]==87 and kdiag["spectrum_counts"]["negative"]==0 and kdiag["rigid_scaled_null_rel"]<=128*EPS*93 and kdiag["physical_rigid_coupling_rank"]==6)
        if any(not v for v in gates.values()): status="BLOCKED_CACHED_NUMERIC_OR_STRUCTURAL_MISMATCH"
        else: status="RECOVERED_CACHE_INDEPENDENTLY_CONSISTENT_PHYSICS_UNQUALIFIED"
        post={k:sha(v) for k,v in paths.items()}
        if pre!=post or own_pre!=sha(Path(__file__)): status="BLOCKED_INPUT_CHANGED_DURING_AUDIT"
        report={"schema":"mixed93_saved_independent_audit.v1","status":status,"dof_count":93,"kkt_order":99,"source_count":27,"global_branch_count":4,
            "input_paths":{"failed_arrays_npz":str(ARRAYS),"producer_result_json":str(FAILED_RESULT),"upstream_result_json":str(RESULT),"upstream_arrays_npz":str(UPSTREAM)},"input_bytes":{"failed_arrays_npz":ARRAYS.stat().st_size},
            "saved_cache_state":"RECOVERED_CACHE_FROM_FAILED_SERIALIZATION_DIAGNOSTIC_ONLY","producer_run":{"return_code":1,"elapsed_seconds":receipt["elapsed_seconds"],"result_json_path":str(FAILED_RESULT),"result_json_present":False,"failure":"metrics JSON serialization failed after NPZ save; result JSON absent"},
            "physical_fracture_status":"UNKNOWN_NOT_ASSESSED","interpretation":case["classification"],
            "pins_pre":pre,"pins_post":post,"reviewer_sha256":own_pre,"source_pins_manifest_sha256":pre["source_pins"],
            "matrix_layout_spectrum_rigid_audit":kdiag,"consistency_gates":gates,"case_audits":case_metrics,
            "checks_scope":"Saved arrays only: independent raw K/f/C/g, physical moments, scaled 99x99 KKT residual/restoration, reconstructed fields; no provider, MLS, quadrature regeneration, author import, or solve.","elapsed_seconds":time.monotonic()-t0}
    payload=json.dumps(report,indent=2,allow_nan=False)+"\n"
    if len(payload.encode("utf-8"))>4*1024*1024: raise ValueError("audit JSON exceeds 4 MiB")
    OUT.write_text(payload,encoding="utf-8")
    print(json.dumps({"status":status,"audit":str(OUT),"elapsed_seconds":report["elapsed_seconds"]},allow_nan=False))
    return 0 if status=="RECOVERED_CACHE_INDEPENDENTLY_CONSISTENT_PHYSICS_UNQUALIFIED" else 2

def main():
    p=argparse.ArgumentParser(description=__doc__); mode=p.add_mutually_exclusive_group(required=True); mode.add_argument("--preflight",action="store_true"); mode.add_argument("--audit",action="store_true"); args=p.parse_args()
    if args.preflight: print(json.dumps(preflight(),indent=2,allow_nan=False)); return 0
    return run_audit()

if __name__=="__main__": raise SystemExit(main())
