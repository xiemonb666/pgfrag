"""Render the audit's saved spectrum/ledger arrays only; performs no reassembly."""
import json, math
from pathlib import Path
from PIL import Image, ImageDraw

HERE = Path(__file__).resolve().parent
r = json.loads((HERE/"audit.json").read_text(encoding="utf-8"))
if r["status"] != "PASS_ALGEBRAIC_REANALYSIS_ONLY": raise ValueError("audit did not pass algebra checks")
s, e = r["spectrum"], r["energy_ledger_and_G_comparison"]
old, new = [float(x) for x in s["baseline_eigenvalues"]], [float(x) for x in s["extended_eigenvalues"]]
w,h=1040,690; im=Image.new("RGB",(w,h),"white"); d=ImageDraw.Draw(im)
left,right,top,bottom=92,30,72,224; xmax,ymax=w-right,h-bottom
d.text((left,16),"Baseline vs fixed-load shifted-front stiffness spectrum",fill=(20,20,20))
d.text((left,40),"diagonally scaled eigenvalues; log10(|lambda| / max |lambda|), orange/purple = near-zero modes",fill=(55,55,55))
maxabs=max(max(abs(x) for x in old),max(abs(x) for x in new),1e-300)
def points(e):
    z=[max(-16.0,min(0.0,math.log10(max(abs(x)/maxabs,1e-16)))) for x in e]
    return [(round(left+i/(len(e)-1)*(xmax-left)),round(top+(0-v)/16*(ymax-top))) for i,v in enumerate(z)]
p0,p1=points(old),points(new)
for tick in (0,-4,-8,-12,-16):
    y=round(top+(0-tick)/16*(ymax-top)); d.line((left,y,xmax,y),fill=(225,225,225)); d.text((49,y-7),str(tick),fill=(60,60,60))
d.line((left,top,left,ymax),fill=(25,25,25),width=2); d.line((left,ymax,xmax,ymax),fill=(25,25,25),width=2)
d.line(p0,fill=(40,100,170),width=2); d.line(p1,fill=(195,70,60),width=2)
for eig, pts, tol, color in ((old,p0,s["baseline_zero_tolerance"],(230,140,20)),(new,p1,s["extended_zero_tolerance"],(120,45,170))):
    for i,v in enumerate(eig):
        if abs(v)<=s["extended_zero_tolerance"]:
            x,y=pts[i]; d.ellipse((x-3,y-3,x+3,y+3),fill=color)
d.text((left,ymax+13),"eigenvalue rank (sorted)",fill=(30,30,30))
d.text((left,ymax+37),"Blue: baseline operator K     Red: independently refronted K     Six rigid modes each",fill=(40,80,140))
d.text((left,ymax+63),f"U: {e['baseline']['U_J']:.6e} -> {e['extended']['U_J']:.6e} J     W: {e['baseline']['W_J']:.6e} -> {e['extended']['W_J']:.6e} J",fill=(25,25,25))
d.text((left,ymax+86),f"Pi: {e['baseline']['Pi_U_minus_W_J']:.6e} -> {e['extended']['Pi_U_minus_W_J']:.6e} J     deltaPi={e['delta_Pi_new_minus_old_J']:.6e} J",fill=(25,25,25))
d.text((left,ymax+109),f"one-face L*da={e['one_sided_area_L_times_delta_a_m2']:.3e} m^2     Gfinite={e['finite_change_G_J_per_m2']:.6e}     Gan={e['analytic_local_modeIII_G_J_per_m2']:.6e} J/m^2",fill=(25,25,25))
d.text((left,ymax+132),f"finite/local ratio={e['finite_over_analytic']:.3f}; unqualified finite-change diagnostic",fill=(110,45,35))
im.save(HERE/"spectrum_energy_ledger.png",format="PNG",optimize=True)
print(json.dumps({"status":"RENDERED_AUDIT_SPECTRUM_LEDGER","output":"spectrum_energy_ledger.png"}))
