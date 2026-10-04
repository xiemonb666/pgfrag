from pathlib import Path
import hashlib, json, os, sys, traceback
import akantu as aka
import numpy as np

EXP=Path(__file__).resolve().parent
INPUT=EXP/"input"
OUT=Path("/mnt/c/Users/xiemo/Desktop/ShardForge/relicweave-synlab/tmp/review_20261004/implementation/cohesive_bilinear_patch_20261004/preflight_dump_v3")
TOL=1e-12
OBS={}
STATE={}
def cp(x): return np.asarray(x).copy()
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def tri_area(x): return 0.5*float(np.linalg.norm(np.cross(x[1]-x[0],x[2]-x[0])))
def input_faces(p):
    q=Path(p).read_text().splitlines(); a=q.index("$Elements")+2; b=q.index("$EndElements"); out=[]
    for line in q[a:b]:
        f=list(map(int,line.split()))
        if f[1]==2 and f[2]>=1 and f[3]==2: out.append([x-1 for x in f[5:8]])
    return out
def save_native_npz():
    if STATE: np.savez_compressed(OUT/"preflight_native_state.npz",**STATE)
def init():
    aka.parseInput(str(INPUT/"material.dat"))
    mesh=aka.Mesh(3); mesh.read(str(INPUT/"block_2x2x4mm.msh"))
    groups=[{"name":g.getName(),"element_count":len(g)} for g in mesh.elementGroups()]
    OBS["physical_groups"]=groups
    OBS["interface_group_element_count"]=next((g["element_count"] for g in groups if g["name"]=="interface"),0)
    OBS["solver_steps"]=0
    if OBS["interface_group_element_count"]!=8: raise RuntimeError("mesh-read interface group count is not 8")
    n0=cp(mesh.getNodes()); t0=cp(mesh.getConnectivity(aka._tetrahedron_4))
    STATE.update(input_nodes=n0,input_bulk_tets=t0)
    OBS["input_counts"]={"nodes":len(n0),"bulk_tets":len(t0)}
    if n0.shape!=(45,3) or t0.shape!=(96,4): raise RuntimeError(f"input counts unexpected: {n0.shape}, {t0.shape}")
    model=aka.SolidMechanicsModelCohesive(mesh)
    model.initFull(_analysis_method=aka._explicit_lumped_mass,_is_extrinsic=False)
    model.assembleMassLumped()
    n=cp(mesh.getNodes()); t=cp(mesh.getConnectivity(aka._tetrahedron_4)); c=cp(mesh.getConnectivity(aka._cohesive_3d_6)); mass=cp(model.getMass())
    STATE.update(native_rest_nodes=n,native_bulk_tets=t,cohesive_connectivity=c,mass=mass)
    OBS["native_counts"]={"nodes":len(n),"bulk_tets":len(t),"cohesive_3d6":len(c)}
    OBS["native_count_gate"]=(len(n),len(t),len(c))==(54,96,8)
    mass_comp=np.asarray(mass.sum(axis=0) if mass.ndim==2 else [mass.sum()],dtype=float)
    energies={k:float(model.getEnergy(k)) for k in ("kinetic","potential","dissipated","reversible","cohesive contact")}
    OBS["mass"]={"shape":list(mass.shape),"all_positive":bool(np.all(mass>0)),"dof_sums_kg":mass_comp.tolist()}
    OBS["energies_J"]=energies
    save_native_npz()
    if not OBS["native_count_gate"]: raise RuntimeError(f"native counts unexpected: {OBS['native_counts']}")
    if len(t)!=len(t0) or t.shape[1]!=4: raise RuntimeError("native bulk cell/local-ordinal lineage unavailable")
    lineage={}
    for old,new in zip(t0,t):
        for src,dst in zip(old,new):
            src,dst=int(src),int(dst)
            if not np.array_equal(n[dst],n0[src]): raise RuntimeError("bulk cell/local-ordinal rest coordinate mismatch")
            if dst in lineage and lineage[dst]!=src: raise RuntimeError("conflicting bulk cell/local-ordinal ancestry")
            lineage[dst]=src
    OBS["exact_bulk_cell_local_ordinal_rest_map"]=True
    if len(c)!=8 or c.shape[1]!=6: raise RuntimeError("expected eight native cohesive 3d6 elements")
    pairset=set(); ordinal_equal=True
    for row in c:
        a,b=[int(x) for x in row[:3]],[int(x) for x in row[3:]]
        la,lb=[lineage[x] for x in a],[lineage[x] for x in b]
        if sorted(la)!=sorted(lb): raise RuntimeError("cohesive side ancestry sets do not match")
        if la!=lb: ordinal_equal=False; raise RuntimeError("cohesive side ancestors differ by local ordinal")
        for src in la:
            ia=a[la.index(src)]; ib=b[lb.index(src)]
            if not np.array_equal(n[ia],n[ib]): raise RuntimeError("duplicated rest coordinates not bit-equal")
            pairset.add(tuple(sorted((ia,ib))))
    OBS["cohesive_side_ordinal_ancestry_equal"]=ordinal_equal
    OBS["exact_duplicate_pairs"]={"count":len(pairset),"coordinates_bit_equal":True}
    low=set(map(int,t[n[t].mean(axis=1)[:,2]<0.002].ravel()))
    high=set(map(int,t[n[t].mean(axis=1)[:,2]>0.002].ravel()))
    OBS["centroid_half_partition"]={"lower_nodes":len(low),"upper_nodes":len(high),"disjoint":not bool(low&high),"covers_54":(low|high)==set(range(54))}
    f0=input_faces(INPUT/"block_2x2x4mm.msh")
    ar0=sum(tri_area(n0[x]) for x in f0); arA=sum(tri_area(n[x[:3]]) for x in c); arB=sum(tri_area(n[x[3:]]) for x in c)
    OBS["areas_m2"]={"input_interface_face_count":len(f0),"input_interface":ar0,"native_side_A":arA,"native_side_B":arB}
    gates={"physical_interface_group_8":OBS["interface_group_element_count"]==8,
      "native_54_96_8":OBS["native_count_gate"],"input_interface_faces_8_area_4e-6":len(f0)==8 and abs(ar0-4e-6)<=TOL,
      "exact_bulk_cell_local_ordinal_coordinates":OBS["exact_bulk_cell_local_ordinal_rest_map"],
      "cohesive_sides_ordinal_ancestry_equal":ordinal_equal,
      "nine_duplicate_pairs_exact":len(pairset)==9,
      "rest_side_areas_4e-6":abs(arA-4e-6)<=TOL and abs(arB-4e-6)<=TOL,
      "centroid_halves_disjoint_cover_all":not bool(low&high) and (low|high)==set(range(54)),
      "mass_positive_and_4e-5_each_dof":bool(np.all(mass>0)) and len(mass_comp)==3 and bool(np.all(np.abs(mass_comp-4e-5)<=TOL)),
      "initial_energies_zero":all(v==0.0 for v in energies.values())}
    OBS["qa_gates"]=gates
    OBS["all_qa_gates_pass"]=all(gates.values())
    OBS["preflight_measurements_saved_before_dump"]=True
    (OUT/"partial_measurements_before_dump.json").write_text(json.dumps(OBS,indent=2)+"\n")
    save_native_npz()
    if not OBS["all_qa_gates_pass"]: raise RuntimeError("one or more initialization QA gates failed")
    fields=["opening","tractions","damage","delta max"]
    OBS["native_dumper_fields_requested"]=fields
    for name in fields: model.addDumpFieldToDumper("cohesive elements",name)
    os.chdir(OUT)
    before={p.resolve() for p in OUT.iterdir() if p.is_file()}
    model.dump("cohesive elements",0)
    after=[p for p in OUT.iterdir() if p.is_file() and p.resolve() not in before]
    dump_files=[{"name":p.name,"bytes":p.stat().st_size} for p in sorted(after)]
    OBS["native_dump"]={"status":"called","dumper":"cohesive elements","step":0,"files":dump_files}
    if not dump_files or any(x["bytes"]<=0 for x in dump_files): raise RuntimeError("native dumper created no nonempty file")
    OBS["status"]="DUMP_PREFLIGHT_PASS_NO_SOLVER_STEPS"
    OBS["provenance"]={"mesh_sha256":sha(INPUT/"block_2x2x4mm.msh"),"material_sha256":sha(INPUT/"material.dat"),
      "manifest_sha256":sha(EXP/"manifest_dump_v3.json"),"driver_sha256":sha(__file__),"python":sys.version.split()[0],
      "wheel_to_tag_correspondence":"UNKNOWN"}
    return OBS
def main():
    OUT.mkdir(parents=True,exist_ok=True)
    if (OUT/"preflight_result.json").exists(): raise SystemExit("v3 output already has a result; refusing overwrite")
    os.chdir(OUT)
    result={}
    try:
        result=init()
    except Exception as e:
        result={"status":"FAIL_NO_RETRY","error":repr(e),"traceback":traceback.format_exc(),**OBS}
    finally:
        result["argv"] = sys.argv
        save_native_npz()
        (OUT/"preflight_result.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,sort_keys=True))
    return 0 if result.get("status")=="DUMP_PREFLIGHT_PASS_NO_SOLVER_STEPS" else 2
if __name__=="__main__": raise SystemExit(main())