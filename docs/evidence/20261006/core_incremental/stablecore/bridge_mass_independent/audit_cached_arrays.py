"""Independent cached-array mass, generalized-spectrum, and bridge readback."""
from __future__ import annotations
import hashlib, json, math, time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[4]
OUT = Path(__file__).resolve().parent
MASS_PATH = ROOT / "pgfrag/tmp/review_20261006/full3d_consistent_mass_replay_20261006/run01/mass_result.json"
OP_PATH = ROOT / "tmp/review_20261006/full3d_cracked_elastic_operator_20261006/run01/result.json"
BRIDGE_PATH = ROOT / "pgfrag/tmp/review_20261006/full3d_static_solve_bridge_20261006/run02/static_solve_bridge_result.json"
BRIDGE_RUN = ROOT / "pgfrag/tmp/review_20261006/full3d_static_solve_bridge_20261006/run02/outer_execution_receipt.json"
BRIDGE_SCRIPT = ROOT / "pgfrag/experiments/full3d_static_solve_bridge_20261006/static_solve_bridge.py"
BRIDGE_CASE = ROOT / "pgfrag/experiments/full3d_static_solve_bridge_20261006/static_solve_case.json"
PINS_PATH = ROOT / "tmp/review_20261006/full3d_cracked_elastic_operator_20261006/SOURCE_PINS.json"
CASE_PATH = ROOT / "pgfrag/experiments/full3d_consistent_mass_replay_20261006/fixed_case.json"
PIN = {"mass": "326df86f867eec47141fce58a748eb68403853adce3a2d8b3b59d8f992ff15e1",
       "operator": "3b31e7d40c1d8e4b5d829f0b093fee4edadc966bde4086a4b3f4ba48fdf89cd2",
       "bridge": "ef66783e2961dd77f34785a99e1b397210ccc0beb6e8fa2979678aa4992a0ce3"}
ND, NB, NV, NQ, NS = 162, 54, 256, 384, 27
COMP, KINDS = ("ux", "uy", "uz"), ("base", "front_f1_unshifted")

def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()

def read(path: Path, expected: str | None = None) -> tuple[dict, bytes]:
    raw = path.read_bytes()
    if len(raw) > 8 * 1024 * 1024 or (expected and sha(raw) != expected):
        raise RuntimeError(f"size/hash pin failure: {path}")
    obj = json.loads(raw.decode("utf-8"), parse_constant=lambda x: (_ for _ in ()).throw(ValueError(x)))
    if not isinstance(obj, dict): raise RuntimeError(f"JSON root is not an object: {path}")
    return obj, raw

def check(ok: bool, message: str) -> None:
    if not ok: raise RuntimeError(message)

def ar(obj: dict, key: str, shape: tuple[int, ...]) -> np.ndarray:
    a = np.asarray(obj[key], dtype=np.float64)
    if a.shape != shape or not np.isfinite(a).all(): raise RuntimeError(f"invalid {key}: {a.shape}")
    return np.ascontiguousarray(a)

def relmax(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.max(np.abs(a-b)) / max(float(np.max(np.abs(b))), np.finfo(float).tiny))

def strain(grad: np.ndarray) -> np.ndarray:
    e = np.empty((len(grad), 6)); e[:, :3] = grad[:, (0, 1, 2), (0, 1, 2)]
    e[:, 3] = grad[:,0,1] + grad[:,1,0]; e[:,4] = grad[:,1,2] + grad[:,2,1]
    e[:, 5] = grad[:,0,2] + grad[:,2,0]
    return e

def tensor(v: np.ndarray) -> np.ndarray:
    s = np.zeros((len(v), 3, 3)); s[:,0,0],s[:,1,1],s[:,2,2] = v[:,0],v[:,1],v[:,2]
    s[:,0,1]=s[:,1,0]=v[:,3]; s[:,1,2]=s[:,2,1]=v[:,4]; s[:,0,2]=s[:,2,0]=v[:,5]
    return s

def reconstruct(coef: np.ndarray, keys: list[dict], rows: dict, val: np.ndarray,
                grad: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    u=np.zeros((NV,3)); g=np.zeros((NV,3,3)); ci={c:i for i,c in enumerate(COMP)}
    for slot,key in enumerate(keys):
        r=rows[(key["source_id"],key["basis_kind"])]; c=ci[key["component"]]
        u[:,c] += coef[slot]*val[:NV,r]; g[:,c,:] += coef[slot]*grad[:NV,r,:]
    return u,g

def plots(me, rho_v, direct, coeff, expected, ee, es, se, ss) -> list[str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dr=np.diag(direct)/np.diag(expected); cr=np.diag(coeff)/np.diag(expected)
    fig,ax=plt.subplots(1,2,figsize=(9,3.5),constrained_layout=True)
    ax[0].semilogy(np.arange(1,len(me)+1),me,color="#245b83",lw=1.4)
    ax[0].set(xlabel="Sorted eigenvalue index",ylabel="Eigenvalue of M/(rho V)",title="Saved mass spectrum")
    ax[0].grid(True,which="both",alpha=.25); x=np.arange(6); w=.36
    ax[1].bar(x-w/2,dr,w,label="Volume quadrature"); ax[1].bar(x+w/2,cr,w,label="Cached coefficients",alpha=.8)
    ax[1].axhline(1,color="black",ls="--",lw=1); ax[1].set_xticks(x,("TX","TY","TZ","RX","RY","RZ"))
    ax[1].set(ylabel="Diagonal / analytic reference",title="Six rigid kinetic diagonals")
    ax[1].grid(True,axis="y",alpha=.25); ax[1].legend(fontsize=7)
    fig.suptitle("Mass audit - rho=1 kg/m3 toy scale; not calibrated")
    p1=OUT/"mass_spectrum_rigid_kinetic.png"; fig.savefig(p1,dpi=90); plt.close(fig)
    labels=("exx","eyy","ezz","gxy","gyz","gxz"); colors=plt.cm.tab10(np.arange(6))
    fig,ax=plt.subplots(1,2,figsize=(8,4),constrained_layout=True)
    for j,label in enumerate(labels):
        ax[0].scatter(ee[:,j],es[:,j],s=8,alpha=.4,color=colors[j],label=label)
        ax[1].scatter(se[:,j],ss[:,j],s=8,alpha=.4,color=colors[j],label=label.replace("e","s").replace("g","t"))
    for a,xdata,ydata,title,unit in ((ax[0],ee,es,"Engineering strain","1"),(ax[1],se,ss,"Full Voigt stress","Pa")):
        lo=min(float(xdata.min()),float(ydata.min())); hi=max(float(xdata.max()),float(ydata.max()))
        pad=max((hi-lo)*.04,1e-15); a.plot([lo-pad,hi+pad],[lo-pad,hi+pad],"k--",lw=1)
        a.set(xlabel=f"Exact ({unit})",ylabel=f"Solved ({unit})",title=title); a.grid(True,alpha=.2); a.legend(fontsize=6,ncol=2)
    fig.suptitle("Saved 256 volume rows - solved diagnostic vs manufactured reference")
    p2=OUT/"bridge_strain_stress_scatter.png"; fig.savefig(p2,dpi=90); plt.close(fig)
    return [str(p1),str(p2)]

def main() -> int:
    start=time.monotonic()
    for name in ("audit_result.json","mass_spectrum_rigid_kinetic.png","bridge_strain_stress_scatter.png"):
        if (OUT/name).exists(): raise RuntimeError(f"refusing overwrite/retry: {name}")
    dense=ND*ND*8; vb=NQ*NB*8; gb=NQ*NB*3*8; fb=NV*3*6*8
    raster=max(810*315*4,720*360*4)
    array_est=16*dense+8*vb+4*gb+8*fb+4*NB*NB*8+raster
    check(array_est<=8*1024*1024,"conservative array-plus-plot estimate exceeds 8 MiB")
    mass,mraw=read(MASS_PATH,PIN["mass"]); op,oraw=read(OP_PATH,PIN["operator"])
    bridge,braw=read(BRIDGE_PATH,PIN["bridge"]); run,rraw=read(BRIDGE_RUN); pins,praw=read(PINS_PATH); case,craw=read(CASE_PATH)
    check(mass.get("status")=="MASS_MATRIX_PASS_UNQUALIFIED" and mass.get("solve_allowed") is False,"mass status/gate mismatch")
    check(op.get("case_id")=="modeiii_pu_full3d_elastic_operator_20261006","operator identity mismatch")
    check(bridge.get("status")=="SOLVED_DIAGNOSTIC_ONLY_UNQUALIFIED","bridge status mismatch")
    check(bridge["source_result"]["sha256"]==PIN["operator"]==mass["input"]["sha256"],"upstream pin mismatch")
    check(sha(praw)==mass["input"]["operator_source_pins_sha256"],"source-pins hash mismatch")
    check(sha(craw)==mass["input"]["mass_case_sha256"],"mass case hash mismatch")
    check(run["actual_return_code"]==0 and run["primary_result"]["sha256"]==PIN["bridge"],"bridge outer receipt mismatch")
    check(sha(BRIDGE_SCRIPT.read_bytes())==run["pre_run_pins"]["source_sha256"]==run["post_run_pins"]["source_sha256"],"bridge source pin changed")
    check(sha(BRIDGE_CASE.read_bytes())==run["pre_run_pins"]["case_sha256"]==run["post_run_pins"]["case_sha256"],"bridge case pin changed")
    check(case.get("case_id")=="modeiii_pu_full3d_consistent_mass_replay_20261006"
          and case.get("input",{}).get("sha256")==PIN["operator"],"mass case identity/input mismatch")
    for item in [pins["entrypoint"],*pins["dependencies"]]:
        check(sha((ROOT/item["path"]).read_bytes())==item["sha256"],f"source pin changed: {item['path']}")
    replay=op["static_solve_replay"]; reg=replay["registry"]; sources=[f"node-{i:03d}" for i in range(NS)]
    check(reg["source_ids"]==sources,"source ordering mismatch")
    keys=reg["ordered_keys"]; exp=[(s,c,b) for s in sources for c in COMP for b in KINDS]
    check(len(keys)==ND and len(reg["ordered_canonical_ids"])==ND,"162 semantic keys missing")
    for i,(k,(s,c,b)) in enumerate(zip(keys,exp)):
        canon=json.dumps([s,c,b,0 if b=="base" else 1],separators=(",",":"))
        check((k["slot"],k["source_id"],k["component"],k["basis_kind"],k["basis_index"],k["canonical_id"])==
              (i,s,c,b,0 if b=="base" else 1,canon) and reg["ordered_canonical_ids"][i]==canon,f"semantic slot {i} mismatch")
    rl=replay["basis_row_ids"]; rows={(r["source_id"],r["basis_kind"]):i for i,r in enumerate(rl)}
    check(len(rl)==NB and len(rows)==NB and [(r["source_id"],r["basis_kind"]) for r in rl]==
          [(s,b) for s in sources for b in KINDS],"54 scalar basis rows mismatch")
    q=replay["quadrature"]; basis=replay["basis"]; opdata=replay["operator"]
    val=ar(basis,"values",(NQ,NB)); grad=ar(basis,"gradients_per_m",(NQ,NB,3))
    w=ar(q,"volume_weights_m3",(NV,)); check(bool(np.all(w>0)),"volume weights are not positive")
    pos=ar(q,"positions_xyz_m",(NQ,3))[:NV]
    rho=float(mass["density_kg_per_m3"]); check(rho==1.0,"toy density changed")
    volume=float(w.sum()); rho_v=rho*volume
    slots=np.asarray([rows[(k["source_id"],k["basis_kind"])] for k in keys])
    comps=np.asarray([COMP.index(k["component"]) for k in keys])
    scalar=(val[:NV].T*(rho*w))@val[:NV]; same=comps[:,None]==comps[None,:]
    M=scalar[np.ix_(slots,slots)]*same
    Ms=ar({"M":mass["mass_matrix_kg"]},"M",(ND,ND))
    merr=float(np.linalg.norm(M-Ms,"fro")/np.linalg.norm(Ms,"fro"))
    masym=float(np.linalg.norm(M-M.T,"fro")/np.linalg.norm(M,"fro"))
    me=np.linalg.eigvalsh(M/rho_v); mes=ar(mass["spectrum"],"eigenvalues_M_over_rhoV",(ND,))
    ematch=float(np.max(np.abs(me-mes))/max(float(np.max(np.abs(mes))),1e-300)); L=np.linalg.cholesky(M)
    center=ar(replay["sources"],"rigid_mode_center_xyz_m",(3,)); src=ar(replay["sources"],"positions_xyz_m",(NS,3))
    centererr=float(np.max(np.abs(src.mean(axis=0)-center))); phys=np.zeros((NV,3,6)); eye=np.eye(3)
    for a in range(3):
        phys[:,a,a]=1.; phys[:,:,3+a]=np.cross(eye[a],pos-center)
    modes=ar(opdata,"rigid_mode_coefficients",(ND,6)); cr=np.zeros((NB,3,6))
    for i,k in enumerate(keys): cr[rows[(k["source_id"],k["basis_kind"])],COMP.index(k["component"])]=modes[i]
    prebuilt=np.einsum("qr,rcm->qcm",val[:NV],cr,optimize=True)
    ms=np.array([1.,1.,1.,100.,100.,100.]); recr=float(np.max(np.abs((prebuilt-phys)*ms[None,None,:])))
    direct=np.einsum("q,qcm,qcn->mn",rho*w,phys,phys,optimize=True); coeff=modes.T@M@modes
    kinerr=float(np.max(np.abs((direct-coeff)*ms[:,None]*ms[None,:]/rho_v)))
    geo=op["fixed_case"]["geometry"]; radius=float(geo["radius_m"]); length=float(geo["length_m"])
    inertia=rho_v*np.array([radius**2/4+length**2/12]*2+[radius**2/2])
    expect=np.diag([rho_v,rho_v,rho_v,*inertia])
    anerr=float(np.max(np.abs((direct-expect)*ms[:,None]*ms[None,:]/rho_v)))
    volerr=abs(volume-math.pi*radius**2*length)/(math.pi*radius**2*length)

    K=ar(opdata,"K_N_per_m",(ND,ND)); C=ar(bridge["gauge"],"C_dimensionless",(6,ND))
    f=ar({"f":bridge["solve"]["raw_applied_load_N"]},"f",(ND,))
    d=ar({"d":bridge["solve"]["coefficients_m"]},"d",(ND,))
    lam=ar({"l":bridge["solve"]["lagrange_multipliers_N"]},"l",(6,))
    g=ar({"g":bridge["gauge"]["g_target_m"]},"g",(6,))
    khash=sha(np.asarray(K,dtype="<f8",order="C").tobytes(order="C"))
    eq=K@d+C.T@lam-f; con=C@d-g
    weak=float(np.linalg.norm(eq)/max(np.linalg.norm(K@d)+np.linalg.norm(C.T@lam)+np.linalg.norm(f),1e-300))
    gauge=float(np.linalg.norm(con)/max(np.linalg.norm(g),np.linalg.norm(C@d),1e-300))
    U=.5*float(d@K@d); W=float(d@f); gl=float(g@lam); ident=2*U-W+gl
    left=np.linalg.solve(L,K); A=np.linalg.solve(L,left.T).T
    aasy=float(np.linalg.norm(A-A.T,"fro")/np.linalg.norm(A,"fro"))
    ek=np.linalg.eigvalsh(A,UPLO="L"); escale=max(float(np.max(np.abs(ek))),np.finfo(float).tiny)
    ktol=max(64*np.finfo(float).eps*ND,1e-12)*escale
    nzero=int(np.count_nonzero(np.abs(ek)<=ktol)); nneg=int(np.count_nonzero(ek < -ktol))
    omega=math.sqrt(float(ek[-1])); young=float(op["fixed_case"]["material_and_manufactured_field"]["young_modulus_Pa"])
    dt=2/omega; nondim=omega*radius/math.sqrt(young/rho)
    rr=L.T@modes; rres=np.linalg.norm(A@rr,axis=0)/np.maximum(escale*np.linalg.norm(rr,axis=0),1e-300)

    ef=bridge["solved_field_at_saved_volume_rows_only"]
    us,gs=reconstruct(d,keys,rows,val,grad)
    ec=ar(bridge["manufactured_reference_separate_from_solution"],"exact_coefficients_m",(ND,))
    ue,ge=reconstruct(ec,keys,rows,val,grad)
    su=ar(ef,"solved_displacement_m",(NV,3)); se=ar(ef,"solved_engineering_strain_dimensionless",(NV,6))
    eu=ar(ef,"exact_displacement_m",(NV,3)); ee=ar(ef,"exact_engineering_strain_dimensionless",(NV,6))
    es=strain(gs); ecalc=strain(ge); esave=strain(ar(replay["exact_field"],"displacement_gradient_dimensionless",(NQ,3,3))[:NV])
    mat=op["fixed_case"]["material_and_manufactured_field"]; E=float(mat["young_modulus_Pa"]); nu=float(mat["poisson_ratio"])
    mu=E/(2*(1+nu)); lame=E*nu/((1+nu)*(1-2*nu)); D=np.zeros((6,6))
    D[:3,:3]=lame; D[0,0]+=2*mu; D[1,1]+=2*mu; D[2,2]+=2*mu; D[3,3]=D[4,4]=D[5,5]=mu
    sv=es@D.T; sve=ecalc@D.T; sv0=esave@D.T
    st=ar(replay["exact_field"],"stress_Pa",(NQ,3,3))[:NV]
    stress_exact=tensor(sv0); stress_solved=tensor(sv)
    sterr=float(np.linalg.norm(stress_exact-st)/np.linalg.norm(st))
    svsaved=np.column_stack([st[:,0,0],st[:,1,1],st[:,2,2],st[:,0,1],st[:,1,2],st[:,0,2]])
    checks={
        "mass_rebuild_matches_saved":merr<=1e-10,"mass_symmetry":masym<=1e-12,"mass_spectrum_matches_saved":ematch<=1e-8,
        "unmodified_mass_cholesky":True,"rigid_field_reconstruction":recr<=1e-8,"all_six_kinetic_forms":kinerr<=1e-8,
        "analytic_cylinder_six_form":anerr<=1e-10,"volume_matches_cylinder":volerr<=1e-12,
        "bridge_raw_K_pin":khash==bridge["solve"]["raw_K_sha256_float64_little_endian_C_order"],
        "bridge_weak_equilibrium":weak<1e-10,"bridge_physical_gauge":gauge<1e-10,"bridge_energy_identity":abs(ident)<1e-15,
        "bridge_field_reconstruction":relmax(us,su)<1e-10 and relmax(es,se)<1e-10,
        "exact_field_reconstruction":relmax(ue,eu)<1e-10 and relmax(ecalc,ee)<1e-10,"exact_stress_full_voigt":sterr<1e-10,
        "generalized_six_rigid_modes":nzero==6 and nneg==0 and aasy<1e-8}
    images=plots(me,rho_v,direct,coeff,expect,ee,es,svsaved,sv)
    report={
      "schema":"independent_cached_mass_bridge_audit.v1",
      "status":"CACHED_ARRAYS_AUDIT_PASS_UNQUALIFIED" if all(checks.values()) else "CACHED_ARRAYS_AUDIT_REVIEW_REQUIRED",
      "classification":"NUMERICAL_DIAGNOSTIC_ONLY_NOT_PHYSICAL_FRACTURE_OR_MATERIAL_PASS","checks":checks,
      "inputs":{"mass_result":{"path":str(MASS_PATH),"sha256":sha(mraw)},"operator_result":{"path":str(OP_PATH),"sha256":sha(oraw)},
        "bridge_result":{"path":str(BRIDGE_PATH),"sha256":sha(braw)},"bridge_run_receipt_sha256":sha(rraw),
        "bridge_script_sha256":sha(BRIDGE_SCRIPT.read_bytes()),"bridge_case_sha256":sha(BRIDGE_CASE.read_bytes()),
        "operator_source_pins_sha256":sha(praw),"mass_case_sha256":sha(craw),
        "bridge_outer_execution":{"return_code":run["actual_return_code"],"elapsed_seconds":run["elapsed_seconds"],
          "stdout_sha256":run["stdout"]["sha256"],"stderr_sha256":run["stderr"]["sha256"],
          "input_sha256":run["pre_run_pins"]["input_sha256"],"source_sha256":run["pre_run_pins"]["source_sha256"],
          "case_sha256":run["pre_run_pins"]["case_sha256"]},
        "script":{"path":str(Path(__file__).resolve()),"sha256":sha(Path(__file__).read_bytes())}},
      "mass":{"density_kg_per_m3":rho,"density_role":mass["density_role"],"volume_m3":volume,"rhoV_kg":rho_v,
        "independent_matrix_sha256_float64_le_C":sha(np.asarray(M,dtype="<f8",order="C").tobytes(order="C")),
        "relative_frobenius_error_vs_saved":merr,"relative_asymmetry":masym,"eigen_min_M_over_rhoV":float(me[0]),
        "eigen_max_M_over_rhoV":float(me[-1]),"condition_estimate":float(me[-1]/me[0]),"spectrum_saved_max_relative_error":ematch,
        "translation_mass_diagonal_kg":np.diag(direct)[:3].tolist(),"rigid_kinetic_direct_6x6":direct.tolist(),
        "rigid_kinetic_from_coefficients_6x6":coeff.tolist(),"analytic_cylinder_inertia_kg_m2":inertia.tolist(),
        "analytic_six_form_scaled_max_error":anerr,"rigid_reconstruction_scaled_max_error":recr,
        "six_kinetic_forms_scaled_max_error":kinerr,"source_center_error_m":centererr,"volume_relative_error":volerr},
      "cached_K_over_M_spectrum":{"units":"s^-2","method":"A=L^-1 K L^-T from raw K and unmodified Cholesky(M); no truncation/regularization",
        "relative_asymmetry_A":aasy,"eigen_min":float(ek[0]),"eigen_max":float(ek[-1]),"near_zero_tolerance":ktol,
        "near_zero_count":nzero,"significant_negative_count":nneg,"first_12_signed_eigenvalues_s^-2":ek[:12].tolist(),
        "rigid_mode_residuals":rres.tolist(),"omega_max_rad_per_s":omega,"central_difference_dt_upper_bound_s":dt,
        "dimensionless_omega_R_over_sqrt_E_over_rho":nondim,
        "dt_interpretation":"toy rho and synthetic E only; not a physical time-step recommendation"},
      "bridge_cached_readback":{"status":bridge["status"],"weak_equilibrium_relative":weak,"gauge_relative":gauge,"energy_U_J":U,
        "external_work_J":W,"g_dot_lambda_J":gl,"energy_identity_defect_J":ident,
        "recorded_energy_identity_defect_J":bridge["solve"]["energy_identity_defect_2U_minus_W_plus_g_dot_lambda_J"],
        "displacement_reconstruction_relative":relmax(us,su),"strain_reconstruction_relative":relmax(es,se),
        "exact_displacement_reconstruction_relative":relmax(ue,eu),"exact_strain_reconstruction_relative":relmax(ecalc,ee),
        "exact_full_tensor_stress_relative_error":sterr,"solved_vs_exact_stress_relative_error":float(np.linalg.norm(tensor(sv)-st)/np.linalg.norm(st)),
        "images":images},
      "array_budget":{"conservative_explicit_peak_including_one_plot_raster_bytes":array_est,"cap_bytes":8388608,
        "JSON_objects_and_BLAS_workspace_excluded":True},
      "elapsed_seconds":time.monotonic()-start,
      "limitations":["one fixed quadrature rule is not a convergence study","rho=1 kg/m3 and synthetic E are not calibrated ceramic properties",
        "generalized spectrum and time-step bound are cached linear diagnostics, not a dynamics simulation",
        "bridge output remains an unqualified manufactured-field diagnostic"]}
    (OUT/"audit_result.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({"status":report["status"],"elapsed_seconds":report["elapsed_seconds"],"checks":checks,
      "omega_max_rad_per_s":omega,"dt_upper_bound_s":dt,"images":images},separators=(",",":")))
    return 0 if all(checks.values()) else 3

if __name__=="__main__":
    try: raise SystemExit(main())
    except Exception as exc:
        import sys
        print(f"AUDIT_BLOCKED: {type(exc).__name__}: {exc}",file=sys.stderr); raise SystemExit(2)
