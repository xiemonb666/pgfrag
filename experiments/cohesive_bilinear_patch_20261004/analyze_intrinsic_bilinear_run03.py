from __future__ import annotations
import json, re, hashlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

RAW = Path(r"C:\Users\xiemo\Desktop\ShardForge\relicweave-synlab\tmp\review_20261004\implementation\cohesive_bilinear_patch_20261004\controlled_intrinsic_bilinear_modeI_201_run03")
EXP = Path(r"C:\Users\xiemo\Desktop\ShardForge\pgfrag\experiments\cohesive_bilinear_patch_20261004")
MANIFEST = EXP / "manifest_bilinear_modeI_201_run03.json"
RECEIPT = RAW.parent / "controlled_intrinsic_bilinear_modeI_201_run03_launch_receipts" / "run03_execution_receipt.json"
DERIVED = RAW.parent / "controlled_intrinsic_bilinear_modeI_201_run03_analysis"
ENERGY_KEYS = ("kinetic", "potential", "dissipated", "reversible", "cohesive contact")
N = 201
def sha256(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def write_json(p, x): p.write_text(json.dumps(x, indent=2, allow_nan=False) + "\n", encoding="utf-8")
def analytic(s, sigma, gc, d0):
    dc = 2.0 * gc / sigma
    return sigma*s/d0 if s <= d0 else (sigma*(dc-s)/(dc-d0) if s <= dc else 0.0)
def ids(paths):
    found=[]
    for p in paths:
        m=re.search(r"_(\d{4})\.(?:pvtu|vtu)$", p.name)
        if m: found.append(int(m.group(1)))
    return sorted(found)
def main():
    if DERIVED.exists(): raise SystemExit("refusing existing derived output: " + str(DERIVED))
    manifest=json.loads(MANIFEST.read_text(encoding="utf-8"))
    status=json.loads((RAW/"run_status.json").read_text(encoding="utf-8"))
    started=json.loads((RAW/"run_started.json").read_text(encoding="utf-8"))
    receipt=json.loads(RECEIPT.read_text(encoding="utf-8")) if RECEIPT.is_file() else {}
    trace=[json.loads(x) for x in (RAW/"trace.jsonl").read_text(encoding="utf-8").splitlines() if x.strip()]
    with np.load(RAW/"native_history_partial_or_complete.npz", allow_pickle=False) as z:
        a={k:z[k].copy() for k in z.files}
    with np.load(RAW/"initial_state.npz", allow_pickle=False) as z:
        initial={k:z[k].copy() for k in z.files}
    sigma=float(manifest["material"]["cohesive"]["sigma_c_Pa"])
    gc=float(manifest["material"]["cohesive"]["G_c_J_m2"])
    d0=float(manifest["material"]["cohesive"]["delta_0_m"])
    dc=float(manifest["material"]["cohesive"]["delta_c_theory_m"])
    pfiles=sorted((RAW/"paraview").rglob("*.pvtu"))
    vfiles=sorted((RAW/"paraview").rglob("*.vtu"))
    pv_ids, v_ids=ids(pfiles),ids(vfiles)
    sfiles=sorted((RAW/"samples").glob("point_*.npz"))
    smap={int(m.group(1)):p for p in sfiles if (m:=re.fullmatch(r"point_(\d{3})\.npz",p.name))}
    trace_ids=[int(r["point_index"]) for r in trace]
    expected=list(range(N))
    # Cross-check all 201 closed point files against the aggregate arrays.
    sample_mismatch=[]
    maps={"opening":"opening","tractions":"tractions","damage":"damage","delta_max":"delta_max",
          "displacement":"displacement","internal_force":"internal_force","mass":"mass",
          "rest_nodes":"rest_nodes","bulk_tets":"bulk_tets","cohesive_connectivity":"cohesive_connectivity",
          "energies_J":"energies_J","mass_component_sums_kg":"mass_component_sums_kg",
          "point_index":"point_index","solver_step":"solver_step","commanded_opening_m":"commanded_opening_m"}
    for i in sorted(smap):
        with np.load(smap[i], allow_pickle=False) as z: row={k:z[k].copy() for k in z.files}
        for sk,ak in maps.items():
            sv=row[sk]
            av=a[ak][i]
            if sk in ("damage","delta_max"): sv=sv[:,0]
            if not np.array_equal(sv,av): sample_mismatch.append({"point":i,"field":sk})
    # Arrays and indexing, retaining native signs and values.
    idx=a["point_index"]; steps=a["solver_step"]; cmd=a["commanded_opening_m"]
    opening=a["opening"]; traction=a["tractions"]; damage=a["damage"]; delta=a["delta_max"]
    energy=a["energies_J"]; mass=a["mass"]; rest=a["rest_nodes"]; tets=a["bulk_tets"]; coh=a["cohesive_connectivity"]
    all_finite=all(np.isfinite(a[k]).all() for k in ("opening","tractions","damage","delta_max","energies_J","mass","rest_nodes","bulk_tets","cohesive_connectivity","displacement","internal_force"))
    trace_counts=np.asarray([[r["counts"][k] for k in ("nodes","bulk_tets","cohesive_3d6")] for r in trace],dtype=np.int64)
    pre=np.asarray([float(r["prescribed_readback_error_m"]) for r in trace])
    post=np.asarray([float(r["post_solve_boundary_readback_error_m"]) for r in trace])
    blocked=np.asarray([r.get("native_blocked_dofs_all_true_pre_solve") is True for r in trace])
    counts_ok=bool(np.array_equal(trace_counts,np.tile([54,96,8],(len(trace),1))) and np.array_equal(a["counts"],trace_counts))
    rest_ok=bool(np.array_equal(rest,np.repeat(initial["rest_nodes"][None,:,:],N,axis=0)))
    tet_ok=bool(np.array_equal(tets,np.repeat(initial["bulk_tets"][None,:,:],N,axis=0)))
    coh_ok=bool(np.array_equal(coh,np.repeat(initial["cohesive_connectivity"][None,:,:],N,axis=0)))
    coh0=initial["cohesive_connectivity"]
    face_pair_ok=bool(np.array_equal(initial["rest_nodes"][coh0[:,:3]],initial["rest_nodes"][coh0[:,3:]]))
    mass_sums=a["mass_component_sums_kg"]
    mass_ok=bool(np.isfinite(mass).all() and (mass>0).all() and np.all(np.abs(mass_sums-4e-5)<=1e-12))
    bc_ok=bool(blocked.all() and pre.max(initial=0)<=1e-12 and post.max(initial=0)<=1e-12)
    target=np.zeros_like(opening); target[:,:,2]=cmd[:,None]
    opening_err=np.abs(opening-target)
    opening_ok=bool(opening_err.max()<=1e-12)
    expected_t=np.asarray([analytic(float(s),sigma,gc,d0) for s in cmd])
    traction_error=traction.copy(); traction_error[:,:,0:2]=traction_error[:,:,0:2]
    traction_error[:,:,2]-=expected_t[:,None]
    traction_abs=float(np.linalg.norm(traction_error,axis=2).max())
    traction_ok=traction_abs<=0.01*sigma
    d=damage.copy(); dm=delta.copy()
    damage_ok=bool(np.isfinite(d).all() and (d>=0).all() and (d<=1).all() and np.all(np.diff(d,axis=0)>=-1e-12))
    delta_ok=bool(np.isfinite(dm).all() and np.all(np.diff(dm,axis=0)>=-1e-12))
    energies_ok=bool(np.isfinite(energy).all())
    diss=energy[:,2]; rev=energy[:,3]
    area_cell=np.asarray([0.5*np.linalg.norm(np.cross(initial["rest_nodes"][f[1]]-initial["rest_nodes"][f[0]],initial["rest_nodes"][f[2]]-initial["rest_nodes"][f[0]])) for f in coh0])
    area=float(area_cell.sum()); gca=gc*area
    finalD=float(diss[-1]); d_rel=abs(finalD-gca)/gca if gca else None
    endpoint_t=float(max(np.max(np.abs(traction[0,:,2])),np.max(np.abs(traction[-1,:,2]))))
    d_dot=np.einsum("tci,tci->tc",opening,traction)
    work=0.0
    for i in range(1,N):
        work+=float(np.sum(area_cell*np.einsum("ci,ci->c",0.5*(traction[i]+traction[i-1]),opening[i]-opening[i-1])))
    hist_complete=(status.get("status")=="RUN_201_POINTS_COMPLETE" and len(trace)==N and len(smap)==N and len(pv_ids)==N and len(v_ids)==N and trace_ids==expected and sorted(smap)==expected and pv_ids==expected and v_ids==expected and np.array_equal(idx,expected) and np.array_equal(steps,np.arange(1,N+1)))
    gates={
      "closed_201_trace_sample_and_native_frames":hist_complete,
      "aggregate_exactly_matches_201_closed_samples":len(sample_mismatch)==0,
      "finite_native_arrays":bool(all_finite),
      "counts_54_96_8_all_steps":counts_ok,
      "rest_tets_cohesive_connectivity_exact_all_steps":bool(rest_ok and tet_ok and coh_ok),
      "paired_rest_faces_exact":face_pair_ok,
      "positive_mass_and_4e-5kg_per_component":mass_ok,
      "fresh_blocked_and_prescribed_readback_le_1e-12m":bc_ok,
      "signed_opening_vector_matches_command_le_1e-12m":opening_ok,
      "signed_traction_vector_matches_analytic_le_1pct_sigma":bool(traction_ok),
      "native_damage_finite_bounded_and_monotone":damage_ok,
      "native_delta_max_finite_and_monotone":delta_ok,
      "raw_dissipated_energy_min_ge_minus_1e-12J":bool(energies_ok and diss.min()>=-1e-12),
      "final_dissipation_vs_GcA_within_1pct":bool(energies_ok and d_rel is not None and d_rel<=0.01),
      "final_reversible_abs_le_1e-12J":bool(energies_ok and abs(float(rev[-1]))<=1e-12),
      "final_damage_all_cells_one_within_1e-12":bool(np.max(np.abs(d[-1]-1.0))<=1e-12),
      "final_delta_max_matches_command_le_1e-12m":bool(np.max(np.abs(dm[-1]-cmd[-1]))<=1e-12),
      "initial_and_final_signed_traction_near_zero_le_1pct_sigma":endpoint_t<=0.01*sigma}
    raw_times=[{"before_s":r.get("time_before_step_s"),"after_s":r.get("time_after_step_s")} for r in trace]
    DERIVED.mkdir(parents=False)
    # 4 panels: signed response, native signed energies, and actual-scale rest faces.
    fig,axs=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
    ax=axs[0,0]
    for c in range(8): ax.plot(opening[:,c,2]*1e9,traction[:,c,2]/1e6,color="steelblue",alpha=.35,lw=.8)
    x=np.linspace(0,1.2*dc,301)
    ax.plot(x*1e9,np.asarray([analytic(float(s),sigma,gc,d0) for s in x])/1e6,"--",color="darkorange",lw=1.6,label="analytic bilinear reference")
    ax.plot(opening[:,:,2].mean(axis=1)*1e9,traction[:,:,2].mean(axis=1)/1e6,color="navy",lw=1.6,label="native signed z, 8-cell mean")
    ax.set(xlabel="native opening z (nm)",ylabel="native traction z (MPa)",title="Signed Mode-I response"); ax.grid(alpha=.25); ax.legend(fontsize=8)
    ax=axs[0,1]
    xnm=cmd*1e9
    for j,label in ((2,"dissipated raw signed"),(3,"reversible raw signed"),(4,"cohesive contact raw"),(1,"bulk potential raw"),(0,"kinetic raw")):
        ax.plot(xnm,energy[:,j]*1e6,label=label,lw=1)
    ax.axhline(gca*1e6,color="black",ls=":",label="Gc × rest area")
    ax.set(xlabel="commanded opening (nm)",ylabel="native energy (µJ)",title="Native energy history; raw signs"); ax.grid(alpha=.25); ax.legend(fontsize=7)
    for ax,side,color,faces in ((axs[1,0],"A","seagreen",coh0[:,:3]),(axs[1,1],"B","darkorange",coh0[:,3:])):
        polys=[(initial["rest_nodes"][f,:2]*1000).tolist() for f in faces]
        ax.add_collection(PolyCollection(polys,facecolors=color,edgecolors="black",alpha=.5,linewidths=.6))
        xy=np.concatenate([initial["rest_nodes"][f,:2] for f in faces])*1000
        ax.scatter(xy[:,0],xy[:,1],s=4,c="black"); ax.set_xlim(xy[:,0].min()-.1,xy[:,0].max()+.1); ax.set_ylim(xy[:,1].min()-.1,xy[:,1].max()+.1)
        ax.set_aspect("equal"); ax.set(xlabel="x (mm)",ylabel="y (mm)",title=f"Rest interface side {side}; actual scale"); ax.grid(alpha=.2)
    fig.suptitle("UNVALIDATED controlled intrinsic bilinear toy; rest geometry only")
    png=DERIVED/"run03_bilinear_modeI_review.png"; fig.savefig(png,dpi=160); plt.close(fig)
    summary={
      "analysis_status":"RUN03_DERIVED_NUMERICAL_EVIDENCE",
      "run_id":manifest["run_id"],"native_status":status.get("status"),"solver_steps":status.get("solver_steps_completed"),
      "execution_receipt":{"python_exit_code":receipt.get("python_exit_code_file"),"wsl_host_exit_code":receipt.get("wsl_host_exit_code"),"powershell_exit_code":receipt.get("powershell_orchestration_script_exit_code"),"timeout_wrapper_exit_code":receipt.get("timeout_wrapper_exit_code"),"elapsed_wall_s":receipt.get("elapsed_s"),"clock_note":"Trace clock fields preserved raw; null values remain unknown."},
      "indices":{"trace":trace_ids,"samples":sorted(smap),"pvtu":pv_ids,"vtu":v_ids,"aggregate":idx.tolist(),"solver_steps":steps.tolist(),"counts":[len(trace),len(smap),len(pv_ids),len(v_ids)]},
      "gates":{k:("PASS" if v else "FAIL") for k,v in gates.items()},
      "measurements":{
        "thresholds":{"opening_abs_m":1e-12,"traction_vector_abs_Pa":0.01*sigma,"dissipated_min_J":-1e-12,"final_D_vs_GcA_relative":0.01,"final_reversible_abs_J":1e-12},
        "max_signed_opening_vector_error_m":float(opening_err.max()),
        "max_signed_traction_vector_error_Pa":traction_abs,"max_native_xy_traction_Pa":float(np.linalg.norm(traction[:,:,:2],axis=2).max()),
        "opening_z_first_last_m":opening[[0,-1],:,2].tolist(),"traction_z_first_last_Pa":traction[[0,-1],:,2].tolist(),
        "T_dot_opening_first_last_Pa_m_raw_signed":d_dot[[0,-1]].tolist(),
        "energy_keys":list(ENERGY_KEYS),"initial_energy_J_raw":energy[0].tolist(),"final_energy_J_raw":energy[-1].tolist(),
        "dissipated_min_raw_J":float(diss.min()),"dissipated_min_point_index":int(idx[np.argmin(diss)]),
        "GcA_reference_J":gca,"final_D_relative_error_vs_GcA":d_rel,"final_reversible_J_raw":float(rev[-1]),
        "independent_signed_area_traction_opening_work_J":work,"work_is_diagnostic_not_native_energy":True,
        "area_m2":area,"cell_area_m2":area_cell.tolist(),
        "native_damage_representation":"8 cohesive-cell averages, scalar/cell; not quadrature-point data",
        "damage_initial_peak_final_per_cell":[d[0].tolist(),d.max(axis=0).tolist(),d[-1].tolist()],
        "damage_final_cell_spread":float(np.ptp(d[-1])),"delta_max_initial_final_m_per_cell":[dm[0].tolist(),dm[-1].tolist()],
        "delta_max_final_matches_command_max_abs_m":float(np.max(np.abs(dm[-1]-cmd[-1]))),
        "all_native_arrays_finite":bool(all_finite),"sample_aggregate_mismatch_count":len(sample_mismatch),
        "counts_unique":np.unique(trace_counts,axis=0).tolist(),"all_frame_rest_tets_coh_exact":[rest_ok,tet_ok,coh_ok],
        "mass_component_sums_first_last_kg":[mass_sums[0].tolist(),mass_sums[-1].tolist()],
        "all_blocked_pre":bool(blocked.all()),"max_BC_readback_errors_m":[float(pre.max()),float(post.max())],
        "time_fields_first_last_raw":[raw_times[0],raw_times[-1]],
        "native_energy_key_note":"Uses 'cohesive contact'; no fallback to a generic 'contact' key.",
        "field_units_note":"VTK fields do not carry physical units; m/Pa interpretation is semantic inference from manifest and runner contract."},
      "energy_balance":"NOT_APPLICABLE: all XYZ DOFs on both half blocks are prescribed.",
      "scope":"UNVALIDATED toy constitutive patch only; controlled all-DOF Mode-I loading. Not physical calibration, free fracture, mesh-bias evidence, R2, or PRETRAIN_ELIGIBLE.",
      "outputs":{"png":str(png),"summary_json":str(DERIVED/"run03_analysis.json"),"analysis_source":str(Path(__file__)),"analysis_source_sha256":sha256(Path(__file__))}}
    write_json(DERIVED/"run03_analysis.json",summary)
    print(json.dumps({"analysis_status":summary["analysis_status"],"gates":summary["gates"],"final_D_J":finalD,"GcA_J":gca,"PNG":str(png)},allow_nan=False,sort_keys=True))
if __name__=="__main__": main()
