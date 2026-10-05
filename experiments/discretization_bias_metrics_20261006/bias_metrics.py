#!/usr/bin/env python3
"""Compare paired fracture descriptors; never treat fixtures as physical evidence."""
from __future__ import annotations
import argparse, csv, json, math
from pathlib import Path

LEN = {"mm": 1.0, "m": 1000.0, "um": 0.001}
TIME = {"s": 1.0, "ms": 0.001, "us": 1e-6}
MASS = {"kg": 1.0, "g": 0.001}
ENERGY = {"J": 1.0, "mJ": 0.001, "kJ": 1000.0}
QUAL = ("geometry", "native_run", "material", "load", "energy_ledger", "crack_surface")
CAP = 512


class InputError(ValueError):
    pass


def get(obj, path):
    cur = obj
    for part in path.split("."):
        if not isinstance(cur, dict) or part not in cur:
            raise InputError("missing " + path)
        cur = cur[part]
    return cur


def number(value, label):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(float(value)):
        raise InputError(label + " must be finite numeric")
    return float(value)


def digest(value, label):
    if not isinstance(value, str) or len(value.strip()) < 8:
        raise InputError(label + " must be an 8+ character hash/reference")
    return value.strip()


def text(value, label):
    if not isinstance(value, str) or not value.strip():
        raise InputError(label + " must be nonempty text")
    return value.strip()


def scale(units, table, label):
    if units not in table:
        raise InputError("unsupported " + label + " units: " + repr(units))
    return table[units]


def frame(d):
    f = get(d, "physical_frame")
    frame_id = text(get(f, "common_frame_id"), "common_frame_id")
    fac = scale(get(f, "coordinate_units"), LEN, "coordinate")
    m = get(f, "local_to_common")
    if not isinstance(m, list) or len(m) != 4 or any(not isinstance(r, list) or len(r) != 4 for r in m):
        raise InputError("local_to_common must be a 4x4 matrix")
    m = [[number(v, "local_to_common") for v in row] for row in m]
    if any(abs(m[3][j] - v) > 1e-10 for j, v in enumerate((0., 0., 0., 1.))):
        raise InputError("local_to_common must be homogeneous")
    r = [row[:3] for row in m[:3]]
    for i in range(3):
        for j in range(3):
            dot = sum(r[i][k] * r[j][k] for k in range(3))
            if abs(dot - (i == j)) > 1e-8:
                raise InputError("frame rotation must be orthonormal")
    det = r[0][0]*(r[1][1]*r[2][2]-r[1][2]*r[2][1])-r[0][1]*(r[1][0]*r[2][2]-r[1][2]*r[2][0])+r[0][2]*(r[1][0]*r[2][1]-r[1][1]*r[2][0])
    if abs(det - 1.) > 1e-8:
        raise InputError("frame transform must have rotation determinant +1")
    return frame_id, m, fac


def validate(d):
    try:
        if get(d, "schema") != "pgfrag.discretization_bias.descriptor.v1":
            raise InputError("unsupported descriptor schema")
        text(get(d, "run_id"), "run_id")
        if get(d, "evidence.kind") not in ("FIXTURE", "NATIVE_PHYSICAL_RESULT"):
            raise InputError("evidence.kind must be FIXTURE or NATIVE_PHYSICAL_RESULT")
        get(d, "evidence.scope")
        frame(d)
        text(get(d, "geometry.reference"), "geometry.reference")
        digest(get(d, "geometry.sha256"), "geometry.sha256")
        digest(get(d, "load.sha256"), "load.sha256")
        load_units = get(d, "load.units")
        if not isinstance(load_units, str) or not load_units.strip():
            raise InputError("load.units must be a nonempty unit description")
        digest(get(d, "material.parameter_sha256"), "material.parameter_sha256")
        digest(get(d, "material.field_sha256"), "material.field_sha256")
        if isinstance(get(d, "material.random_field_seed"), bool) or not isinstance(get(d, "material.random_field_seed"), int):
            raise InputError("material.random_field_seed must be an integer (0 for uniform)")
        get(d, "material.units")
        digest(get(d, "solver.build_sha256"), "solver.build_sha256")
        digest(get(d, "solver.failure_law_sha256"), "solver.failure_law_sha256")
        text(get(d, "solver.surface_correction"), "solver.surface_correction")
        get(d, "solver.run_status")
        digest(get(d, "discretization.mesh_sha256"), "mesh_sha256")
        for key, ukey, table in (("h","h_units",LEN),("delta","delta_units",LEN),("dt","dt_units",TIME)):
            if number(get(d, "discretization."+key), key) <= 0:
                raise InputError("discretization."+key+" must be positive")
            scale(get(d, "discretization."+ukey), table, key)
        text(get(d, "discretization.orientation_id"), "orientation_id")
        text(get(d, "discretization.translation_id"), "translation_id")
        if not isinstance(get(d, "fracture.fracture_occurred"), bool) or not isinstance(get(d, "fracture.history_complete"), bool):
            raise InputError("fracture flags must be boolean")
        digest(get(d, "fracture.surface_extractor_sha256"), "surface_extractor_sha256")
        digest(get(d, "fracture.path_sampling_policy_sha256"), "path_sampling_policy_sha256")
        get(d, "energy_terms")
        return []
    except InputError as exc:
        return [str(exc)]


def dims(d):
    q = get(d, "discretization")
    return {k: number(q[k],k)*scale(q[uk], tab,k) for k,uk,tab in (
        ("h","h_units",LEN),("delta","delta_units",LEN),("dt","dt_units",TIME))}


def compatible(a,b,factor):
    fixed = ("physical_frame.common_frame_id","geometry.reference","geometry.sha256",
             "load.sha256","load.units","material.parameter_sha256","material.field_sha256",
             "material.random_field_seed","material.units","solver.build_sha256",
             "solver.failure_law_sha256","solver.surface_correction",
             "fracture.surface_extractor_sha256","fracture.path_sampling_policy_sha256")
    for p in fixed:
        try:
            x,y=get(a,p),get(b,p)
        except InputError as exc:
            return "UNKNOWN",[str(exc)]
        if x != y:
            return "INCOMPARABLE",[f"fixed invariant differs: {p}"]
    am,bm=get(a,"discretization"),get(b,"discretization")
    if factor not in ("orientation","translation","h","dt","delta"):
        return "UNKNOWN",[f"unsupported factor: {factor}"]
    av,bv=dims(a),dims(b)
    ids_a=(am["orientation_id"],am["translation_id"])
    ids_b=(bm["orientation_id"],bm["translation_id"])
    reasons=[]
    if factor in ("orientation","translation"):
        idx=0 if factor=="orientation" else 1
        if ids_a[idx]==ids_b[idx]: reasons.append(f"{factor}_id did not change")
        if ids_a[1-idx]!=ids_b[1-idx]: reasons.append("other placement factor changed")
        if am["mesh_sha256"]==bm["mesh_sha256"]: reasons.append("mesh_sha256 identical; remesh not evidenced")
        if any(not math.isclose(av[k],bv[k],rel_tol=1e-10,abs_tol=1e-12) for k in av): reasons.append("h, delta, and dt must stay fixed")
    else:
        if ids_a!=ids_b: reasons.append("orientation/translation changed")
        if factor=="h":
            if math.isclose(av["h"],bv["h"],rel_tol=1e-10,abs_tol=1e-12): reasons.append("h did not change")
            if not math.isclose(av["delta"],bv["delta"],rel_tol=1e-10,abs_tol=1e-12) or not math.isclose(av["dt"],bv["dt"],rel_tol=1e-10,abs_tol=1e-12): reasons.append("delta and dt must stay fixed for h pairing")
            if am["mesh_sha256"]==bm["mesh_sha256"]: reasons.append("mesh hash identical during h refinement")
        elif factor=="dt":
            if math.isclose(av["dt"],bv["dt"],rel_tol=1e-10,abs_tol=1e-12): reasons.append("dt did not change")
            if any(not math.isclose(av[k],bv[k],rel_tol=1e-10,abs_tol=1e-12) for k in ("h","delta")): reasons.append("h and delta must stay fixed for dt pairing")
            if am["mesh_sha256"]!=bm["mesh_sha256"]: reasons.append("mesh hash changed during dt pairing")
        else:
            if math.isclose(av["delta"],bv["delta"],rel_tol=1e-10,abs_tol=1e-12): reasons.append("delta did not change")
            if any(not math.isclose(av[k],bv[k],rel_tol=1e-10,abs_tol=1e-12) for k in ("h","dt")): reasons.append("h and dt must stay fixed for delta pairing")
            if am["mesh_sha256"]!=bm["mesh_sha256"]: reasons.append("mesh hash changed during delta pairing")
        if a["physical_frame"]["local_to_common"] != b["physical_frame"]["local_to_common"]:
            reasons.append("frame transform changed outside orientation/translation pairing")
    if reasons: return "INCOMPARABLE",reasons
    try: frame(a); frame(b)
    except InputError as exc: return "UNKNOWN",[str(exc)]
    return "COMPARABLE",[]


def transform(points,d,units,label):
    if not isinstance(points,list) or not points: raise InputError(label+" must be a nonempty point array")
    _,m,frame_scale=frame(d); s=scale(units,LEN,"point")
    out=[]
    for p in points:
        if not isinstance(p,list) or len(p)!=3: raise InputError(label+" points must be 3-vectors")
        v=[number(x,label)*s for x in p]
        out.append(tuple(sum(m[i][j]*v[j] for j in range(3))+m[i][3]*frame_scale for i in range(3)))
    return out


def capped(points):
    pts=sorted(set(points))
    if len(pts)<=CAP: return pts
    return [pts[round(i*(len(pts)-1)/(CAP-1))] for i in range(CAP)]


def nearest(a,b):
    return [math.sqrt(min(sum((x[k]-y[k])**2 for k in range(3)) for y in b)) for x in a]


def quantile(v,p):
    s=sorted(v); x=p*(len(s)-1); lo=int(x); hi=math.ceil(x)
    return s[lo]+(s[hi]-s[lo])*(x-lo)


def path_metrics(a,b):
    pa=get(a,"fracture.path_samples_local"); pb=get(b,"fracture.path_samples_local")
    aa=capped(transform(pa,a,get(a,"fracture.path_units"),"path"))
    bb=capped(transform(pb,b,get(b,"fracture.path_units"),"path"))
    ab,ba=nearest(aa,bb),nearest(bb,aa)
    return {"path_symmetric_p95_mm":quantile(ab+ba,.95),"path_symmetric_hausdorff_mm":max(max(ab),max(ba))}


def area(d):
    s=get(d,"fracture.surface"); factor=scale(get(s,"coordinate_units"),LEN,"surface")
    _,m,frame_scale=frame(d); verts={}
    for v in get(s,"vertices"):
        key=str(get(v,"id"))
        if key in verts: raise InputError("duplicate surface vertex id")
        p=[number(x,"surface vertex")*factor for x in get(v,"xyz")]
        verts[key]=tuple(sum(m[i][j]*p[j] for j in range(3))+m[i][3]*frame_scale for i in range(3))
    total=0.
    for face in get(s,"triangles"):
        if not isinstance(face,list) or len(face)!=3: raise InputError("triangles must contain three vertex IDs")
        try: a,b,c=(verts[str(k)] for k in face)
        except KeyError as exc: raise InputError("triangle references missing vertex ID") from exc
        u=[b[i]-a[i] for i in range(3)]; v=[c[i]-a[i] for i in range(3)]
        cross=(u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0])
        total+=.5*math.sqrt(sum(x*x for x in cross))
    if not math.isfinite(total) or total<=0: raise InputError("surface area must be positive")
    return total


def energies(d):
    raw=get(d,"energy_terms")
    if not isinstance(raw,dict) or not raw: raise InputError("energy_terms must be nonempty")
    return {k:number(get(v,"value"),"energy."+k)*scale(get(v,"units"),ENERGY,"energy") for k,v in raw.items()}


def mass_shares(d):
    f=get(d,"fragment_evidence")
    kind=get(f,"kind")
    if kind!="closed_surface_geometric_components":
        raise InputError(f"{kind!r} is not geometric fragment evidence")
    val=get(f,"validation")
    if get(val,"closed_surfaces") is not True or get(val,"pairwise_nonoverlap") is not True:
        raise InputError("geometric fragment surfaces are not closed and nonoverlapping")
    coverage=number(get(val,"parent_coverage_error"),"parent_coverage_error")
    parent=number(get(f,"parent_mass"),"parent_mass")*scale(get(f,"parent_mass_units"),MASS,"mass")
    if parent<=0 or coverage<0: raise InputError("invalid parent mass/coverage")
    ids=set(); shares=[]
    for c in get(f,"components"):
        cid=str(get(c,"id"))
        if cid in ids: raise InputError("duplicate fragment id")
        ids.add(cid); digest(get(c,"surface_sha256"),"fragment surface_sha256")
        mass=number(get(c,"mass"),"fragment mass")*scale(get(c,"mass_units"),MASS,"mass")
        if mass<=0: raise InputError("fragment mass must be positive")
        shares.append(mass/parent)
    if not shares: raise InputError("no geometric components")
    return shares,abs(sum(shares)-1.),coverage


def w1(a,b):
    x,y=sorted(a),sorted(b)
    if not x or not y: raise InputError("empty mass distribution")
    i=j=0; rx,ry=1/len(x),1/len(y); total=0.
    while i<len(x) and j<len(y):
        w=min(rx,ry); total+=abs(x[i]-y[j])*w; rx-=w; ry-=w
        if rx<=1e-14:
            i+=1
            if i<len(x): rx=1/len(x)
        if ry<=1e-14:
            j+=1
            if j<len(y): ry=1/len(y)
    return total


def metrics(a,b):
    out={}; missing=[]
    for name,fn in (
        ("path",lambda:path_metrics(a,b)),
        ("surface",lambda:{"surface_area_abs_delta_mm2":abs(area(a)-area(b))}),
        ("initiation",lambda:{"initiation_time_abs_delta_s":abs(
            number(get(a,"fracture.initiation_time"),"initiation_time")*scale(get(a,"fracture.initiation_time_units"),TIME,"time")-
            number(get(b,"fracture.initiation_time"),"initiation_time")*scale(get(b,"fracture.initiation_time_units"),TIME,"time"))})):
        try: out.update(fn())
        except (InputError,KeyError,TypeError) as exc: missing.append(f"{name}: {exc}")
    try:
        ea,eb=energies(a),energies(b)
        if ea.keys()!=eb.keys(): raise InputError("energy term names differ")
        out.update({f"energy.{k}_abs_delta_J":abs(ea[k]-eb[k]) for k in ea})
    except (InputError,KeyError,TypeError) as exc: missing.append("energy: "+str(exc))
    try:
        sa,ca,va=mass_shares(a); sb,cb,vb=mass_shares(b)
        out.update({"fragment_mass_w1_fraction":w1(sa,sb),
                    "fragment_mass_closure_abs_delta":abs(ca-cb),
                    "reference_fragment_mass_closure_error":ca,
                    "candidate_fragment_mass_closure_error":cb,
                    "reference_parent_coverage_error":va,
                    "candidate_parent_coverage_error":vb})
    except (InputError,KeyError,TypeError) as exc: missing.append("fragment_mass_distribution: "+str(exc))
    return out,missing


def tolerance_status(values,tolerances):
    if not isinstance(tolerances,dict) or not tolerances: return "UNKNOWN",["no explicit caller tolerances"]
    missing=[]; exceeded=[]
    for k,v in values.items():
        t=tolerances.get(k)
        if not isinstance(t,dict): missing.append("no tolerance for "+k); continue
        try: maximum=number(get(t,"max_abs"),"tolerance.max_abs"); basis=text(get(t,"basis"),"tolerance basis")
        except InputError as exc: missing.append(str(exc)); continue
        if maximum<0: missing.append("negative tolerance for "+k)
        elif v>maximum: exceeded.append(f"{k}={v:.9g} exceeds {maximum:.9g} ({basis})")
    if exceeded: return "OUT_OF_TOLERANCE",exceeded+missing
    if missing: return "UNKNOWN",missing
    return "WITHIN_TOLERANCES",[]


def physical_errors(d):
    errors=[]
    if get(d,"evidence.kind")!="NATIVE_PHYSICAL_RESULT": errors.append("not NATIVE_PHYSICAL_RESULT")
    if get(d,"solver.run_status")!="COMPLETE": errors.append("native run not COMPLETE")
    q=d.get("qualification")
    if not isinstance(q,dict): return errors+["qualification attestations missing"]
    errors.extend("qualification."+k+" is not PASS" for k in QUAL if q.get(k)!="PASS")
    return errors


def compare(desc,pair):
    rid,cid=pair.get("reference_run_id"),pair.get("candidate_run_id")
    out={"pair_id":pair.get("pair_id","unnamed"),"factor":pair.get("factor"),
         "reference_run_id":rid,"candidate_run_id":cid,"pairing_status":"UNKNOWN",
         "assessment_status":"UNKNOWN","metrics":{},"unavailable_metrics":[],"reasons":[]}
    if rid not in desc or cid not in desc:
        out["reasons"]=["unknown descriptor run_id"]; return out
    a,b=desc[rid],desc[cid]
    errs=validate(a)+validate(b)
    if errs: out["reasons"]=sorted(set(errs)); return out
    status,why=compatible(a,b,str(pair.get("factor","")))
    out["pairing_status"]=status
    if status!="COMPARABLE": out["assessment_status"]=status; out["reasons"]=why; return out
    for d,label in ((a,"reference"),(b,"candidate")):
        if get(d,"solver.run_status")!="COMPLETE": out["reasons"].append(label+" run incomplete")
        if get(d,"fracture.fracture_occurred") is not True: out["reasons"].append(label+" contains no fracture event")
        if get(d,"fracture.history_complete") is not True: out["reasons"].append(label+" history incomplete")
    if out["reasons"]: return out
    try: vals,missing=metrics(a,b)
    except Exception as exc:
        out["reasons"]=[f"metric error {type(exc).__name__}: {exc}"]; return out
    out["metrics"],out["unavailable_metrics"]=vals,missing
    fixtures=any(get(d,"evidence.kind")=="FIXTURE" or get(d,"evidence.scope")=="FIXTURE" for d in (a,b))
    tstatus,treasons=tolerance_status(vals,pair.get("tolerances"))
    if missing: out["reasons"].extend(missing)
    if fixtures:
        out["assessment_status"]="FIXTURE_ONLY"
        out["reasons"].append("fixture metrics test evaluator only; no bias-removal claim")
    elif physical_errors(a) or physical_errors(b):
        out["assessment_status"]="UNKNOWN"; out["reasons"].extend(physical_errors(a)+physical_errors(b))
    elif missing:
        out["assessment_status"]="UNKNOWN"
    elif tstatus=="OUT_OF_TOLERANCE":
        out["assessment_status"]=tstatus; out["reasons"].extend(treasons)
    elif tstatus=="UNKNOWN":
        out["assessment_status"]="UNKNOWN"; out["reasons"].extend(treasons)
    else:
        out["assessment_status"]="WITHIN_TOLERANCES_NOT_BIAS_ELIMINATED"
        out["reasons"].append("within caller tolerances is not proof of zero bias or A12 eligibility")
    return out


def run_study(path):
    study=json.loads(path.read_text(encoding="utf-8"))
    if study.get("schema")!="pgfrag.discretization_bias.study.v1": raise InputError("unsupported study schema")
    if not isinstance(study.get("descriptors"),list) or not study["descriptors"]: raise InputError("descriptors required")
    desc={}; errors=[]
    for rel in study["descriptors"]:
        try:
            d=json.loads((path.parent/rel).read_text(encoding="utf-8")); rid=str(get(d,"run_id"))
            if rid in desc: errors.append("duplicate run_id "+rid)
            else: desc[rid]=d; errors.extend(rid+": "+e for e in validate(d))
        except (OSError,json.JSONDecodeError,InputError) as exc: errors.append(str(rel)+": "+str(exc))
    if not isinstance(study.get("pairs"),list) or not study["pairs"]: raise InputError("pairs required")
    rows=[compare(desc,p) for p in study["pairs"]]
    fixture=any(d.get("evidence",{}).get("kind")=="FIXTURE" for d in desc.values())
    overall="FIXTURE_ONLY" if fixture else overall_status(rows)
    return {"schema":"pgfrag.discretization_bias.report.v1","study_id":study.get("study_id",path.stem),
            "evidence_scope":"FIXTURE_ONLY" if fixture else "NATIVE_CANDIDATE",
            "overall_status":overall,"admission":"NOT_ADMITTED","bias_eliminated":False,
            "input_errors":errors,"pair_results":rows,
            "statement":"Fixture metrics validate evaluator plumbing only; they do not show mesh bias is absent." if fixture
            else "This report does not by itself qualify physical truth, R2, or training eligibility."}


def overall_status(rows):
    states=[r["assessment_status"] for r in rows]
    if "INCOMPARABLE" in states: return "INCOMPARABLE"
    if "OUT_OF_TOLERANCE" in states: return "OUT_OF_TOLERANCE"
    if any(s!="WITHIN_TOLERANCES_NOT_BIAS_ELIMINATED" for s in states): return "UNKNOWN"
    return "WITHIN_TOLERANCES_NOT_BIAS_ELIMINATED"


CSV_FIELDS=("pair_id","factor","reference_run_id","candidate_run_id","pairing_status","assessment_status",
 "path_symmetric_p95_mm","path_symmetric_hausdorff_mm","surface_area_abs_delta_mm2","initiation_time_abs_delta_s",
 "fragment_mass_w1_fraction","fragment_mass_closure_abs_delta","energy_abs_deltas_J","reasons")


def write_outputs(report,outdir):
    outdir.mkdir(parents=True,exist_ok=True)
    (outdir/"metrics_report.json").write_text(json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    with (outdir/"pair_metrics.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=CSV_FIELDS); w.writeheader()
        for r in report["pair_results"]:
            m=r["metrics"]; e={k:v for k,v in m.items() if k.startswith("energy.")}
            w.writerow({"pair_id":r["pair_id"],"factor":r.get("factor"),"reference_run_id":r.get("reference_run_id"),
             "candidate_run_id":r.get("candidate_run_id"),"pairing_status":r["pairing_status"],"assessment_status":r["assessment_status"],
             "path_symmetric_p95_mm":m.get("path_symmetric_p95_mm",""),
             "path_symmetric_hausdorff_mm":m.get("path_symmetric_hausdorff_mm",""),
             "surface_area_abs_delta_mm2":m.get("surface_area_abs_delta_mm2",""),
             "initiation_time_abs_delta_s":m.get("initiation_time_abs_delta_s",""),
             "fragment_mass_w1_fraction":m.get("fragment_mass_w1_fraction",""),
             "fragment_mass_closure_abs_delta":m.get("fragment_mass_closure_abs_delta",""),
             "energy_abs_deltas_J":json.dumps(e,sort_keys=True),"reasons":" | ".join(r["reasons"])})
    png_name = "fixture_metrics.png" if report["evidence_scope"] == "FIXTURE_ONLY" else "metric_comparison.png"
    plot(report,outdir/png_name)


def plot(report,path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows=[r for r in report["pair_results"] if r["metrics"]]
    panels=(("path_symmetric_p95_mm","Symmetric path p95 (mm)"),
            ("surface_area_abs_delta_mm2","Crack area difference (mm^2)"),
            ("initiation_time_abs_delta_s","Initiation time difference (s)"),
            ("fragment_mass_w1_fraction","Fragment mass W1 (fraction)"))
    fig,axes=plt.subplots(2,2,figsize=(10,7),constrained_layout=True)
    for ax,(key,title) in zip(axes.flat,panels):
        vals=[r["metrics"].get(key) for r in rows]; idx=[i for i,x in enumerate(vals) if isinstance(x,(int,float))]
        if idx:
            ax.bar(idx,[vals[i] for i in idx],color="#4477aa")
            ax.set_xticks(idx,[rows[i]["pair_id"] for i in idx],rotation=30,ha="right")
        else: ax.text(.5,.5,"No available metric",ha="center",va="center",transform=ax.transAxes)
        ax.set_title(title); ax.grid(axis="y",alpha=.25)
    title = ("FIXTURE ONLY - evaluator plumbing; no mesh-bias qualification"
             if report["evidence_scope"] == "FIXTURE_ONLY"
             else "Native candidate metrics - review required; no automatic physical admission")
    fig.suptitle(title)
    fig.savefig(path,dpi=150); plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--study",required=True,type=Path)
    p.add_argument("--output-dir",required=True,type=Path)
    a=p.parse_args()
    try:
        report=run_study(a.study.resolve()); write_outputs(report,a.output_dir.resolve())
    except (OSError,json.JSONDecodeError,InputError) as exc: p.error(str(exc))
    print(json.dumps({"study_id":report["study_id"],"overall_status":report["overall_status"],
                      "pair_count":len(report["pair_results"]),"output_dir":str(a.output_dir.resolve())}))
    return 0 if not report["input_errors"] else 2


if __name__=="__main__": raise SystemExit(main())

