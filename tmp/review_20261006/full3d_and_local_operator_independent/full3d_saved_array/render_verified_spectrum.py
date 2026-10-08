"""Render only the saved eigenspectrum already checked by the independent audit."""
from __future__ import annotations

import json
from pathlib import Path

from PIL import Image, ImageDraw

HERE = Path(__file__).resolve().parent
RUN = Path(r"C:\Users\xiemo\Desktop\ShardForge\tmp\review_20261006\full3d_cracked_elastic_operator_20261006\run01\result.json")
saved = json.loads(RUN.read_text(encoding="utf-8"))
audit = json.loads((HERE / "audit.json").read_text(encoding="utf-8"))
eig = [float(x) for x in saved["spectrum"]["symmetric_eigenvalues"]]
if len(eig) != 162 or not all(__import__("math").isfinite(v) for v in eig):
    raise ValueError("expected 162 finite cached eigenvalues")
if audit["status"] != "PASS_SAVED_ARRAY_ALGEBRA_ONLY" or audit["eigenvalue_relative_l2_delta"] > 5e-11:
    raise ValueError("saved eigenspectrum is not independently verified")
ztol = float(audit["zero_tolerance_N_per_m"])
w, h = 960, 560
im = Image.new("RGB", (w, h), "white")
d = ImageDraw.Draw(im)
left, right, top, bottom = 94, 28, 62, 88
xmax, ymax = w - right, h - bottom
d.text((left, 16), "Full3D eigenvalues of independently rebuilt K", fill=(20, 20, 20))
d.text((left, 38), "absolute magnitude / max |lambda| (log10); orange = 6 near-zero modes", fill=(55, 55, 55))
maxabs = max(abs(v) for v in eig)
logs = [max(-16.0, min(0.0, __import__("math").log10(max(abs(v) / maxabs, 1e-16)))) for v in eig]
points = [(round(left + i / (len(eig)-1) * (xmax-left)), round(top + (0-logv) / 16 * (ymax-top))) for i, logv in enumerate(logs)]
for tick in (0, -4, -8, -12, -16):
    y = round(top + (0-tick) / 16 * (ymax-top))
    d.line((left, y, xmax, y), fill=(222, 222, 222), width=1)
    d.text((48, y-7), str(tick), fill=(65, 65, 65))
d.line((left, top, left, ymax), fill=(20, 20, 20), width=2)
d.line((left, ymax, xmax, ymax), fill=(20, 20, 20), width=2)
d.line(points, fill=(35, 89, 145), width=2)
for i, value in enumerate(eig):
    if abs(value) <= ztol:
        x, y = points[i]
        d.ellipse((x-4, y-4, x+4, y+4), fill=(230, 120, 20))
d.text((left, ymax+20), "eigenvalue rank (sorted)", fill=(30, 30, 30))
d.text((left, ymax+43), f"min={min(eig):.3e} N/m    max={max(eig):.6g} N/m    zero tolerance={ztol:.3e} N/m", fill=(50, 50, 50))
im.save(HERE / "spectrum.png", format="PNG", optimize=True)
print(json.dumps({"status": "RENDERED_VERIFIED_CACHED_SPECTRUM", "eigenvalue_count": len(eig), "zero_modes_marked": sum(abs(v) <= ztol for v in eig)}))
