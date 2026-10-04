from pathlib import Path
import hashlib, json, sys, traceback
import akantu as aka
import numpy as np

EXP = Path(__file__).resolve().parent
INPUT = EXP / "input"
OUT = Path("/mnt/c/Users/xiemo/Desktop/ShardForge/relicweave-synlab/tmp/review_20261004/implementation/cohesive_bilinear_patch_20261004")
TOL = 1e-12
def cp(a): return np.asarray(a).copy()
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def area(x):
    return 0.5 * float(np.linalg.norm(np.cross(x[1]-x[0], x[2]-x[0])))
def faces_from_msh(p):
    lines=Path(p).read_text().splitlines(); a=lines.index("$Elements")+2; b=lines.index("$EndElements")
    out=[]
    for line in lines[a:b]:
        f=list(map(int,line.split()))
        if f[1]==2 and f[2]>=1 and f[3]==2: out.append([n-1 for n in f[5:8]])
    return out
def init():
    aka.parseInput(str(INPUT/"material.dat"))
    mesh=aka.Mesh(3); mesh.read(str(INPUT/"block_2x2x4mm.msh"))
    mesh.createGroupsFromMeshData("physical_names")
    n0=cp(mesh.getNodes()); t0=cp(mesh.getConnectivity(aka._tetrahedron_4))
    if n0.shape!=(45,3) or t0.shape!=(96,4): raise RuntimeError(f"input counts {n0.shape} {t0.shape}")
    model=aka.SolidMechanicsModelCohesive(mesh)
    model.initFull(_analysis_method=aka._explicit_lumped_mass,_is_extrinsic=False)
    model.assembleMassLumped()
    n=cp(mesh.getNodes()); t=cp(mesh.getConnectivity(aka._tetrahedron_4)); c=cp(mesh.getConnectivity(aka._cohesive_3d_6))
    if (len(n),len(t),len(c))!=(54,96,8): raise RuntimeError(f"native counts nodes/tet/coh={(len(n),len(t),len(c))}")
    lineage={}
    if len(t)!=len(t0): raise RuntimeError("bulk element count/order changed")
    for old,new in zip(t0,t):
        for src,dst in zip(old,new):
            src,dst=int(src),int(dst)
            if not np.array_equal(n[dst],n0[src]): raise RuntimeError("cell/local-ordinal ancestry coordinate mismatch")
            if dst in lineage and lineage[dst]!=src: raise RuntimeError("conflicting native cell/local-ordinal lineage")
            lineage[dst]=src
    pairs=set(); ordinal_equal=True
    for row in c:
        a,b=[int(x) for x in row[:3]],[int(x) for x in row[3:]]
        la,lb=[lineage[x] for x in a],[lineage[x] for x in b]
        if sorted(la)!=sorted(lb): raise RuntimeError("cohesive sides do not pair by exact input node identity")
        if la!=lb: ordinal_equal=False; raise RuntimeError("cohesive side ancestry differs by local ordinal")
        for src in la:
            ia=a[la.index(src)]; ib=b[lb.index(src)]
            if not np.array_equal(n[ia],n[ib]): raise RuntimeError("duplicated rest coordinates differ")
            pairs.add(tuple(sorted((ia,ib))))
    low=set(map(int,t[n[t].mean(axis=1)[:,2]<0.002].ravel()))
    high=set(map(int,t[n[t].mean(axis=1)[:,2]>0.002].ravel()))
    if low & high or low|high!=set(range(54)): raise RuntimeError("centroid side node sets overlap or do not cover all native nodes")
    f0=faces_from_msh(INPUT/"block_2x2x4mm.msh"); ar0=sum(area(n0[x]) for x in f0)
    arA=sum(area(n[x[:3]]) for x in c); arB=sum(area(n[x[3:]]) for x in c)
    mat=model.getMaterial(1); op=cp(mat.getOpening(aka._cohesive_3d_6,aka._not_ghost)); tr=cp(mat.getTraction(aka._cohesive_3d_6,aka._not_ghost))
    mass=cp(model.getMass())
    mass_components=[float(mass[:,i].sum()) for i in range(mass.shape[1])] if mass.ndim==2 else [float(mass.sum())]
    energies={k:float(model.getEnergy(k)) for k in ("kinetic","potential","dissipated","reversible","cohesive contact")}
    gates={"native_96tet_8coh_54nodes":True,"input_interface_8faces_area":len(f0)==8 and abs(ar0-4e-6)<=TOL,
           "centroid_halves_disjoint_cover54":True,"duplicate_pairs_9_exact_rest":len(pairs)==9,
           "side_areas_4e-6m2":abs(arA-4e-6)<=TOL and abs(arB-4e-6)<=TOL,
           "cohesive_side_ancestry_equal_by_ordinal":ordinal_equal,"opening_traction_native_shapes":op.shape==(24,3) and tr.shape==(24,3)}
    np.savez_compressed(OUT/"preflight_native_state.npz",rest_nodes=n,bulk_tets=t,cohesive_connectivity=c,mass=mass,opening=op,traction=tr)
    if not all(gates.values()): raise RuntimeError(f"init gates failed: {gates}")
    return {"status":"PREFLIGHT_PASS_NO_SOLVER_STEPS","solver_steps":0,"solver_attempts":0,"gates":gates,"counts":{"input_nodes":len(n0),"native_nodes":len(n),"bulk_tets":len(t),"cohesive_3d6":len(c),"lower_nodes":len(low),"upper_nodes":len(high),"duplicate_pairs":len(pairs),"cohesive_side_ancestry_equal_by_ordinal":ordinal_equal,"input_interface_faces":len(f0)},
      "areas_m2":{"input_interface":ar0,"cohesive_side_A":arA,"cohesive_side_B":arB},
      "native_initial":{"opening_shape":list(op.shape),"traction_shape":list(tr.shape),"opening_max_abs_m":float(np.max(np.abs(op))),"traction_max_abs_Pa":float(np.max(np.abs(tr))),"mass_shape":list(mass.shape),"mass_kg_by_dof":mass_components,"energies_J":energies},
      "provenance":{"mesh_sha256":sha(INPUT/"block_2x2x4mm.msh"),"material_sha256":sha(INPUT/"material.dat"),"driver_sha256":sha(__file__),"python":sys.version.split()[0],"wheel_source_correspondence":"UNKNOWN"}}
def main():
    OUT.mkdir(parents=True,exist_ok=True)
    try:
        result=init(); result["argv"]=sys.argv
    except Exception as e:
        result={"status":"FAIL_NO_RETRY","error":repr(e),"traceback":traceback.format_exc(),"argv":sys.argv}
    (OUT/"preflight_result.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,sort_keys=True))
    return 0 if result["status"]=="PREFLIGHT_PASS_NO_SOLVER_STEPS" else 2
if __name__=="__main__": raise SystemExit(main())