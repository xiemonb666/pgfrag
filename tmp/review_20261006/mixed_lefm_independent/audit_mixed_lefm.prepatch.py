#!/usr/bin/env python3
"""Fixed, source-pinned independent checks for the mixed LEFM field.

Uses only the standard library. It imports the target only after its SHA256
matches --expected-source-sha. Stress equations and plane-strain displacement
equations are transcribed independently from Bower, Applied Mechanics of
Solids, section 9.3.1. This is an elastic field audit, not calibration or a
crack-growth check.
"""
import argparse, csv, hashlib, importlib.util, json, math, struct, sys, time, zlib
from pathlib import Path

TOL = 1.0e-10
RADII = (0.002, 0.008)
ANGLES = (-2.4, -1.1, 0.37, 1.6, 2.7)
MODES = (("pure_I", 100.0, 0.0, 0.0), ("pure_II", 0.0, 80.0, 0.0),
         ("pure_III", 0.0, 0.0, 60.0), ("mixed", 100.0, 80.0, 60.0))

def flat(v):
    if isinstance(v, (tuple, list)):
        return [x for y in v for x in flat(y)]
    return [float(v)]

def mag(v):
    return max((abs(x) for x in flat(v)), default=0.0)

def equal(got, want, label, scale=None):
    a, b = flat(got), flat(want)
    if len(a) != len(b):
        raise AssertionError(label + ": shape mismatch")
    s = max(mag(want), 1.0e-300) if scale is None else max(scale, mag(want), 1.0e-300)
    for i, (x, y) in enumerate(zip(a, b)):
        if not math.isfinite(x) or not math.isfinite(y) or abs(x-y) > TOL*s:
            raise AssertionError(f"{label}[{i}]: {x!r} != {y!r} (scale={s:g})")

def reference(x, y, ki, kii, kiii, young=1.0e6, nu=0.3):
    """Bower Cartesian stresses and plane-strain displacement, local xyz."""
    r, t = math.hypot(x, y), math.atan2(y, x)
    sh, ch = math.sin(t/2), math.cos(t/2)
    s3, c3 = math.sin(1.5*t), math.cos(1.5*t)
    f = 1.0 / math.sqrt(2.0*math.pi*r)
    xx = f*(ki*ch*(1-sh*s3) - kii*sh*(2+ch*c3))
    yy = f*(ki*ch*(1+sh*s3) + kii*ch*sh*c3)
    xy = f*(ki*ch*sh*c3 + kii*ch*(1-sh*s3))
    xz, yz = -f*kiii*sh, f*kiii*ch
    stress = ((xx, xy, xz), (xy, yy, yz), (xz, yz, nu*(xx+yy)))
    mu = young/(2*(1+nu)); root = math.sqrt(r/(2*math.pi))
    a, ap = 1-2*nu+sh*sh, 0.5*math.sin(t)
    b, bp = 3-2*nu-sh*sh, -0.5*math.sin(t)
    c, cp = -1+2*nu+sh*sh, 0.5*math.sin(t)
    hs = ((0, a*ch, ap*ch-0.5*a*sh),
          (0, a*sh, ap*sh+0.5*a*ch),
          (1, b*sh, bp*sh+0.5*b*ch),
          (1, c*ch, cp*ch-0.5*c*sh),
          (2, sh, 0.5*ch))
    coeffs = (ki/(mu*math.sqrt(2*math.pi)), ki/(mu*math.sqrt(2*math.pi)),
              kii/(mu*math.sqrt(2*math.pi)), kii/(mu*math.sqrt(2*math.pi)),
              kiii/mu*math.sqrt(2/math.pi))
    u, grad = [0.0]*3, [[0.0]*3 for _ in range(3)]
    for (component, h, hp), q in zip(hs, coeffs):
        u[component] += q*math.sqrt(r)*h
        factor = q/math.sqrt(r)
        grad[component][0] += factor*(0.5*math.cos(t)*h-math.sin(t)*hp)
        grad[component][1] += factor*(0.5*math.sin(t)*h+math.cos(t)*hp)
    return tuple(u), tuple(tuple(row) for row in grad), stress

def branch_reference(x, y, ref):
    r, t = math.hypot(x, y), math.atan2(y, x)
    sh, ch, st, ct = math.sin(t/2), math.cos(t/2), math.sin(t), math.cos(t)
    h = (sh, ch, sh*st, ch*st)
    hp = (0.5*ch, -0.5*sh, 0.5*ch*st+sh*ct, -0.5*sh*st+ch*ct)
    q = math.sqrt(r/ref); g = []
    for v, d in zip(h, hp):
        g.append((q/r*(0.5*ct*v-st*d), q/r*(0.5*st*v+ct*d), 0.0))
    return tuple(q*v for v in h), tuple(g)

def mv(a, v):
    return tuple(sum(a[i][j]*v[j] for j in range(3)) for i in range(3))

def mm(a, b):
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(3)) for j in range(3)) for i in range(3))

def tr(a):
    return tuple(tuple(a[j][i] for j in range(3)) for i in range(3))

def rodrigues():
    n = math.sqrt(14); x, y, z = 1/n, 2/n, 3/n; c, s = math.cos(0.67), math.sin(0.67)
    return ((c+x*x*(1-c), x*y*(1-c)-z*s, x*z*(1-c)+y*s),
            (y*x*(1-c)+z*s, c+y*y*(1-c), y*z*(1-c)-x*s),
            (z*x*(1-c)-y*s, z*y*(1-c)+x*s, c+z*z*(1-c)))

def png_chunk(tag, data):
    return struct.pack(">I", len(data))+tag+data+struct.pack(">I", zlib.crc32(tag+data)&0xffffffff)

def write_quiver(path, rows, radius):
    w = h = 480; cx = cy = w//2; rp = 190
    pix = bytearray([255])*(w*h*3)
    def put(x, y, color):
        if 0 <= x < w and 0 <= y < h:
            i = (y*w+x)*3; pix[i:i+3] = bytes(color)
    def line(x0, y0, x1, y1, color):
        dx, dy = abs(x1-x0), -abs(y1-y0); sx = 1 if x0<x1 else -1; sy = 1 if y0<y1 else -1; err = dx+dy
        while True:
            put(x0,y0,color)
            if x0==x1 and y0==y1: break
            e = 2*err
            if e>=dy: err+=dy; x0+=sx
            if e<=dx: err+=dx; y0+=sy
    def xy(x,y): return cx+round(x/radius*rp), cy-round(y/radius*rp)
    for i in range(361):
        a = 2*math.pi*i/360; put(*xy(radius*math.cos(a), radius*math.sin(a)), (30,30,30))
    line(cx-rp,cy,cx+rp,cy,(180,180,180)); line(cx,cy-rp,cx,cy+rp,(180,180,180))
    for x in range(-rp,0,12): line(cx+x,cy,cx+min(x+6,0),cy,(210,210,210))
    umax = max(math.hypot(r["u"][0],r["u"][1]) for r in rows)
    for row in rows:
        x,y = row["point"][:2]; x0,y0 = xy(x,y); put(x0,y0,(20,70,180))
        if umax:
            dx = row["u"][0]/umax*0.20*radius; dy = row["u"][1]/umax*0.20*radius
            x1,y1 = xy(x+dx,y+dy); line(x0,y0,x1,y1,(190,35,35))
            ang = math.atan2(y1-y0,x1-x0)
            for da in (-0.65,0.65):
                line(x1,y1,round(x1-5*math.cos(ang+da)),round(y1-5*math.sin(ang+da)),(190,35,35))
    raw = b"".join(b"\0"+bytes(pix[y*w*3:(y+1)*w*3]) for y in range(h))
    data = (b"\x89PNG\r\n\x1a\n"+png_chunk(b"IHDR",struct.pack(">2I5B",w,h,8,2,0,0,0))+
            png_chunk(b"IDAT",zlib.compress(raw,9))+png_chunk(b"IEND",b""))
    path.write_bytes(data)

def run(args):
    start = time.monotonic(); source = args.source.resolve(); digest = hashlib.sha256(source.read_bytes()).hexdigest()
    if digest.lower() != args.expected_source_sha.lower():
        raise SystemExit(f"source SHA mismatch: observed {digest}; no import performed")
    spec = importlib.util.spec_from_file_location("pinned_mixed_lefm", source)
    mod = importlib.util.module_from_spec(spec); sys.modules[spec.name] = mod; spec.loader.exec_module(mod)
    calls = [0]; rows = []; cases = []
    for name, ki, kii, kiii in MODES:
        cases.append((name, mod.MixedLEFMCase(0.01,0.02,ki,kii,kiii,1.0e6,0.3)))
    def evaluate(local, case, trace=None):
        calls[0] += 1
        return mod.evaluate_local_point(local, case=case, branch_trace=trace)
    for name, case in cases:
        for r in RADII:
            for angle in ANGLES:
                p = (r*math.cos(angle),r*math.sin(angle),0.0); f=evaluate(p,case)
                if f.status != mod.FieldStatus.OK: raise AssertionError(f"unexpected status {name}: {f.status}")
                u,g,s = reference(p[0],p[1],case.K_I_Pa_sqrt_m,case.K_II_Pa_sqrt_m,case.K_III_Pa_sqrt_m)
                equal(f.displacement_world_m,u,name+" displacement")
                equal(f.displacement_gradient_world,g,name+" analytic gradient")
                equal(f.stress_tensor_world_Pa,s,name+" Bower Cartesian stress",mag(s))
                rows.append({"case":name,"r":r,"theta":angle,"point":p,"u":f.displacement_world_m,
                             "gradient":f.displacement_gradient_world,"stress":f.stress_tensor_world_Pa})
    # Homogeneity: quadrupling r doubles u and halves each gradient component.
    for name in (m[0] for m in MODES):
        for j, angle in enumerate(ANGLES):
            a=next(q for q in rows if q["case"]==name and q["r"]==RADII[0] and q["theta"]==angle)
            b=next(q for q in rows if q["case"]==name and q["r"]==RADII[1] and q["theta"]==angle)
            equal(b["u"],tuple(2*x for x in a["u"]),name+" sqrt(r) homogeneity",mag(b["u"]))
            equal(a["gradient"],tuple(tuple(2*x for x in row) for row in b["gradient"]),name+" gradient homogeneity",mag(a["gradient"]))
    # All four normalized branch values and their local analytic gradients.
    for r in RADII:
        for angle in ANGLES:
            x,y=r*math.cos(angle),r*math.sin(angle); wantv,wantg=branch_reference(x,y,0.01)
            got=mod.branch_functions_local(x,y,reference_length_m=0.01)
            equal(got.values,wantv,"branch values"); equal(got.gradients_per_m,wantg,"branch gradients")
    # Crack-face traction and the signed jump direction for every mode mixture.
    for name,case in cases:
        up=evaluate((-0.004,0.0,0.0),case,"upper"); lo=evaluate((-0.004,0.0,0.0),case,"lower")
        if up.status!=mod.FieldStatus.OK or lo.status!=mod.FieldStatus.OK: raise AssertionError("face trace status")
        scale=mag(up.stress_tensor_world_Pa)
        for f in (up,lo): equal((f.stress_tensor_world_Pa[0][1],f.stress_tensor_world_Pa[1][1],f.stress_tensor_world_Pa[2][1]),(0,0,0),name+" zero face traction",scale)
        mu=case.shear_modulus_Pa; kappa=3-4*case.poisson_ratio; root=math.sqrt(0.004/(2*math.pi))
        jump=(case.K_II_Pa_sqrt_m/mu*root*(kappa+1),case.K_I_Pa_sqrt_m/mu*root*(kappa+1),
              2*case.K_III_Pa_sqrt_m/mu*math.sqrt(2*0.004/math.pi))
        equal(tuple(up.displacement_world_m[i]-lo.displacement_world_m[i] for i in range(3)),jump,name+" actual crack-face jump",mag(jump))
    # Exact front and ambiguous face stay unknown; malformed data fail closed.
    if evaluate((0.0,0.0,0.0),cases[3][1]).status!=mod.FieldStatus.UNKNOWN_FRONT: raise AssertionError("r=0 must be UNKNOWN_FRONT")
    if evaluate((-0.004,0.0,0.0),cases[3][1]).status!=mod.FieldStatus.UNKNOWN_BRANCH_FACE: raise AssertionError("unqualified face must be UNKNOWN_BRANCH_FACE")
    try: evaluate((math.nan,0.0,0.0),cases[3][1])
    except ValueError: pass
    else: raise AssertionError("non-finite point accepted")
    for label,fn in (("non-finite material",lambda: mod.MixedLEFMCase(.01,.02,0,0,0,1e6,math.nan)),
                     ("improper frame",lambda: mod.CrackFrame(axes_world_from_local=((1,0,0),(0,1,0),(0,0,-1))))):
        try: fn()
        except ValueError: pass
        else: raise AssertionError(label+" accepted")
    # Pure Mode III identities from the existing field_reference.py METHOD, evaluated directly.
    for row in (q for q in rows if q["case"]=="pure_III"):
        x,y=row["point"][:2]; r=math.hypot(x,y); t=math.atan2(y,x); sh,ch=math.sin(t/2),math.cos(t/2)
        case=cases[2][1]; mu=case.shear_modulus_Pa; amp=60/mu*math.sqrt(2/math.pi)
        gradf=(-sh/(2*math.sqrt(r)),ch/(2*math.sqrt(r)),0.0)
        equal(row["u"],(0,0,amp*math.sqrt(r)*sh),"Mode III reference displacement")
        equal(row["stress"][0][2:3]+row["stress"][1][2:3],(mu*amp*gradf[0],mu*amp*gradf[1]),"Mode III reference shear")
    # Rodrigues rotation (axis 1,2,3; angle .67) and local API covariance.
    q=rodrigues(); frame=mod.CrackFrame(axes_world_from_local=q); old=cases[3][1]
    rotated=mod.MixedLEFMCase(old.radius_m,old.length_m,old.K_I_Pa_sqrt_m,old.K_II_Pa_sqrt_m,
                              old.K_III_Pa_sqrt_m,old.young_modulus_Pa,old.poisson_ratio,frame)
    baseline=next(v for v in rows if v["case"]=="mixed" and v["r"]==RADII[0] and v["theta"]==ANGLES[0])
    rf=evaluate(baseline["point"],rotated)
    equal(rf.displacement_world_m,mv(q,baseline["u"]),"rotated displacement",mag(rf.displacement_world_m))
    for key in ("gradient","stress"):
        local=baseline[key]; transformed=mm(mm(q,local),tr(q))
        world=rf.displacement_gradient_world if key=="gradient" else rf.stress_tensor_world_Pa
        equal(world,transformed,"rotated "+key,mag(world))
    localg=baseline["gradient"]
    localstrain=tuple(tuple(0.5*(localg[i][j]+localg[j][i]) for j in range(3)) for i in range(3))
    equal(rf.strain_tensor_world,mm(mm(q,localstrain),tr(q)),"rotated strain",mag(rf.strain_tensor_world))
    if calls[0]>160: raise AssertionError("field-call budget exceeded")
    # Small review artifacts are derived only from the sampled actual outputs.
    out=args.out_dir.resolve(); out.mkdir(parents=True,exist_ok=True)
    csvpath=out/"audit_points.csv"
    with csvpath.open("w",newline="",encoding="utf-8") as f:
        w=csv.writer(f); w.writerow(("case","r_m","theta_rad","x_m","y_m","ux_m","uy_m","uz_m","sxx_Pa","syy_Pa","szz_Pa","sxy_Pa","syz_Pa","sxz_Pa"))
        for a in rows:
            s=a["stress"]; w.writerow((a["case"],a["r"],a["theta"],*a["point"][:2],*a["u"],s[0][0],s[1][1],s[2][2],s[0][1],s[1][2],s[0][2]))
    mixed=[a for a in rows if a["case"]=="mixed"]; png=out/"mixed_displacement_quiver.png"; write_quiver(png,mixed,0.01)
    elapsed=time.monotonic()-start
    if elapsed>10: raise AssertionError(f"runtime budget exceeded: {elapsed:.3f}s")
    report={"classification":"INDEPENDENT_ANALYTIC_FIELD_AUDIT","source_sha256":digest,"stress_reference":"Bower, Applied Mechanics of Solids, §9.3.1 Cartesian stresses; plane-strain displacement","mode_iii_reference":"field_reference.py METHOD formulas evaluated independently","fixed_cases":[{"name":n,"KI":a,"KII":b,"KIII":c,"E_Pa":1e6,"nu":0.3} for n,a,b,c in MODES],"radii_m":RADII,"angles_rad":ANGLES,"checks":"all passed","field_evaluations":calls[0],"elapsed_seconds":elapsed,"tolerance_scaled":TOL,"artifacts":{"csv":csvpath.name,"quiver_png":png.name,"quiver":"equal aspect; actual mixed displacement directions; arrows scaled to 0.20R at the largest in-plane displacement"},"scope":"elastic analytic verification only; no calibration, solver, provider, or crack-growth claim"}
    jpath=out/"audit_report.json"; jpath.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    total=sum(p.stat().st_size for p in (csvpath,png,jpath))
    if total>8*1024*1024: raise AssertionError(f"artifact byte budget exceeded: {total}")
    print(json.dumps({"result":"PASS","field_evaluations":calls[0],"elapsed_seconds":elapsed,"artifact_bytes":total,"source_sha256":digest},sort_keys=True))

def main():
    here=Path(__file__).resolve(); default=here.parents[3]/"experiments"/"mixed_lefm_benchmark_20261006"/"mixed_lefm_field.py"
    p=argparse.ArgumentParser(description=__doc__); p.add_argument("--source",type=Path,default=default)
    p.add_argument("--expected-source-sha",required=True); p.add_argument("--out-dir",type=Path,default=here.parent)
    run(p.parse_args())

if __name__=="__main__": main()
