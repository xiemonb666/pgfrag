import datetime as dt, hashlib, json, math, pathlib, sys, time
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = pathlib.Path(r"C:\Users\xiemo\Desktop\ShardForge")
OUT = ROOT / "pgfrag/tmp/review_20261006/solved_j_fsum_independent"
RCP = ROOT / "pgfrag/tmp/review_20261006/solved_field_j_coupling_20261006/runs/run01/run_receipt.json"
PRI = ROOT / "pgfrag/tmp/review_20261006/modeiii_pu_static_solve_v3_20261006/runs/run01/static_solve_diagnostic.json"
STAGE_SHA = "1377e428a12f53bbc9dc9537c4dbf5d0c3f67a937bcdeceee6ddc6c98dfdfac9"
PRI_SHA = "39cb2dabe767998baf5e84b038f39510d6ef7db10bcd6c77298b480045cb24d3"
MU, JEXACT = 384615.3846153846, 0.013

def need(ok, msg):
    if not ok: raise ValueError(msg)
def num(x, where):
    need(type(x) in (int, float) and math.isfinite(float(x)), "bad number: " + where)
    return float(x)
def vec3(x, where):
    need(type(x) is list and len(x) == 3, "bad 3-list: " + where)
    return [num(y, where) for y in x]
def sha(raw): return hashlib.sha256(raw).hexdigest()
def rec(path, raw): return {"path": str(path), "bytes": len(raw), "sha256": sha(raw)}
def save(path, obj): path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")

started, t0 = dt.datetime.now(dt.timezone.utc), time.perf_counter()
rc_raw = RCP.read_bytes(); rc = json.loads(rc_raw)
need(type(rc) is dict, "receipt root must be dict")
stage_path = pathlib.Path(rc["worker_staging_path"])
st_raw, pri_raw = stage_path.read_bytes(), PRI.read_bytes()
inputs = {"source_run_receipt": rec(RCP, rc_raw), "stage": rec(stage_path, st_raw), "primary": rec(PRI, pri_raw)}
need(inputs["stage"]["sha256"] == STAGE_SHA and rc.get("worker_staging_sha256") == STAGE_SHA, "stage hash mismatch")
need(inputs["primary"]["sha256"] == PRI_SHA, "primary hash mismatch")
need(rc.get("worker_staging_bytes") == len(st_raw), "stage byte count mismatch")
need(sum(x["bytes"] for x in inputs.values()) <= 2 * 1024 * 1024, "input byte budget exceeded")
st, pri = json.loads(st_raw), json.loads(pri_raw)
need(type(st) is dict and type(pri) is dict, "JSON roots must be dicts")
need(rc.get("result_status") == "UNKNOWN_WORKER_RESULT_CONTRACT", "parent gate status changed")
ids = pri["registry"]["ordered_canonical_coefficient_ids"]
need(type(ids) is list and len(ids) == 54 and all(type(x) is str for x in ids), "primary IDs invalid")
plevs = {x["name"]: x["solution"]["coefficients_m"] for x in pri["volume_levels"]}
need(len(plevs) == 2, "expected two primary levels")
for name, coeff in plevs.items():
    need(type(coeff) is list and len(coeff) == 54, "primary coefficients invalid: " + name)
    for i, x in enumerate(coeff): num(x, f"primary {name}[{i}]")
bindings = []
for lev in st["source_solution"]["levels"]:
    name, bind = lev["name"], lev["semantic_coefficient_binding"]
    need(type(bind) is dict and bind["ordered_canonical_coefficient_ids"] == ids, "binding IDs differ: " + name)
    coeff = bind["coefficients_m_in_provider_slot_order"]
    need(type(coeff) is list and len(coeff) == 54, "binding coeffs invalid: " + name)
    for i, x in enumerate(coeff): num(x, f"binding {name}[{i}]")
    need(name in plevs and coeff == plevs[name], "binding coefficients differ: " + name)
    bindings.append({"volume_level": name, "coefficient_count": 54, "ids_equal": True, "coefficients_exactly_equal": True})
need(len(bindings) == 2, "expected two coefficient bindings")

contours = st["contours"]
need(type(contours) is list and len(contours) == 8, "expected eight contours")
paths, contour_by_id, total, max_delta = [], {}, 0, 0.0
for c in contours:
    need(type(c) is dict and type(c.get("contour_id")) is str, "contour invalid")
    cid, semi = c["contour_id"], c["ellipse_semiaxes_m"]
    a, b = num(semi["a"], cid + ".a"), num(semi["b"], cid + ".b")
    coords, normals, weights = c["coordinates_m"], c["normals_outward"], c["weights_ds_m"]
    need(a > 0 and b > 0 and type(coords) is list and type(normals) is list and type(weights) is list, "contour arrays invalid: " + cid)
    n = len(coords)
    need(n > 0 and len(normals) == n and len(weights) == n and c["sample_count"] == n, "contour lengths differ: " + cid)
    for i in range(n): vec3(coords[i], f"{cid}.coord[{i}]"); vec3(normals[i], f"{cid}.normal[{i}]"); num(weights[i], f"{cid}.weight[{i}]")
    by = c["by_solved_volume_level"]
    need(type(by) is dict and set(by) == set(plevs), "contour levels invalid: " + cid)
    row = {"contour_id": cid, "ellipse_semiaxes_m": {"a": a, "b": b}, "rotation_rad": num(c["rotation_rad"], cid + ".rotation"), "quadrature_order_per_smooth_arc": num(c["quadrature_order_per_smooth_arc"], cid + ".order"), "sample_count": n, "levels": {}}
    for name, lev in by.items():
        grads, uz = lev["solved_grad_u_z_m_per_m_full_xyz"], lev["solved_uz_m"]
        need(type(grads) is list and len(grads) == n and type(uz) is list and len(uz) == n, "solved lengths differ: " + cid + "/" + name)
        for i in range(n): vec3(grads[i], f"{cid}/{name}.grad[{i}]"); num(uz[i], f"{cid}/{name}.uz[{i}]")
        integ = lev["integrator_result"]
        need(type(integ) is dict and integ.get("status") == "OK", "integrator status invalid: " + cid + "/" + name)
        samples = integ["samples"]
        need(type(samples) is list and len(samples) == n, "sample count differs: " + cid + "/" + name)
        reported = num(integ["J_J_per_m2"], "stored J " + cid + "/" + name)
        terms, stored = [], []
        for i, s in enumerate(samples):
            need(type(s) is dict, "sample row invalid")
            g, nv = vec3(s["grad_u_z_m_per_m"], "sample grad"), vec3(s["normal_outward"], "sample normal")
            vec3(s["coordinate_m"], "sample coordinate")
            ds = num(s["weight_ds_m"], "sample weight")
            need(s["coordinate_m"] == coords[i] and s["normal_outward"] == normals[i] and s["weight_ds_m"] == weights[i] and s["grad_u_z_m_per_m"] == grads[i], "saved arrays disagree: " + cid + "/" + name)
            gx, gy, nx, ny = g[0], g[1], nv[0], nv[1]
            terms.append((MU * 0.5 * (gx * gx + gy * gy) * nx - MU * (gx * nx + gy * ny) * gx) * ds)
            stored.append(num(s["weighted_contribution_J_per_m2"], "stored sample contribution"))
        calc, stored_sum = math.fsum(terms), math.fsum(stored)
        delta = calc - reported; max_delta = max(max_delta, abs(delta))
        need(math.isclose(calc, reported, rel_tol=1e-10, abs_tol=1e-12) and math.isclose(calc, stored_sum, rel_tol=1e-10, abs_tol=1e-12), "independent J mismatch: " + cid + "/" + name)
        row["levels"][name] = {"sample_count": n, "recomputed_J_J_per_m2": calc, "stored_J_J_per_m2": reported, "difference_from_stored_J_J_per_m2": delta, "sum_stored_sample_contributions_J_per_m2": stored_sum, "relative_error_to_analytic_pct": (calc / JEXACT - 1) * 100, "analytic_reference_J_per_m2": JEXACT}
    paths.append(row); contour_by_id[cid] = c; total += n
need(total == 224, "expected 224 stored positions per level")

zrows = []
for z in st["z_reflection_controls"]["points_by_solved_volume_level"]:
    need(type(z["points"]) is list and len(z["points"]) == 12, "z control pair count invalid")
    g = num(z["summary"]["maximum_abs_grad_z_m_per_m_at_controls"], "z max grad")
    zrows.append({"volume_level": z["volume_level"], "pair_count": 12, "maximum_abs_grad_z_m_per_m_at_controls": g})
need(len(zrows) == 2 and all(z["maximum_abs_grad_z_m_per_m_at_controls"] > 0 for z in zrows), "expected nonzero z controls")
q12 = [abs(p["levels"]["volume_gl4x16x4"]["relative_error_to_analytic_pct"]) for p in paths if int(p["quadrature_order_per_smooth_arc"]) == 12]
q16 = [abs(p["levels"]["volume_gl4x16x4"]["relative_error_to_analytic_pct"]) for p in paths if int(p["quadrature_order_per_smooth_arc"]) == 16]
qd = [abs(p["levels"]["volume_gl4x16x4"]["recomputed_J_J_per_m2"] - p["levels"]["volume_gl4x8x4"]["recomputed_J_J_per_m2"]) for p in paths]
summary = {"path_count": 8, "stored_positions_per_level": total, "saved_gradient_samples_across_two_levels": total * 2, "maximum_abs_difference_from_stored_J_J_per_m2": max_delta, "q12_abs_error_to_analytic_pct_range": [min(q12), max(q12)], "q16_abs_error_to_analytic_pct_range": [min(q16), max(q16)], "maximum_abs_q12_q16_difference_J_per_m2": max(qd)}
results = {"schema": "solved_j_fsum_independent_review.v1", "formula": "fsum((mu*0.5*(gx^2+gy^2)*nx-mu*(gx*nx+gy*ny)*gx)*ds)", "shear_modulus_Pa": MU, "analytic_reference_J_per_m2": JEXACT, "inputs": inputs, "parent_gate_status_preserved": rc["result_status"], "semantic_bindings": bindings, "summary": summary, "paths": paths, "z_reflection_controls": zrows, "interpretation": "saved midsection contour functional; nonzero grad_z controls do not establish a full-domain 2D reduction or physical G"}
result_path = OUT / "numeric_results.json"; save(result_path, results)

fig, (axg, axj) = plt.subplots(1, 2, figsize=(16, 7), gridspec_kw={"width_ratios": [1, 1.25]})
for p in paths:
    c = contour_by_id[p["contour_id"]]; pts = sorted(c["coordinates_m"], key=lambda x: math.atan2(x[1], x[0]))
    xs, ys = [x[0] * 1000 for x in pts], [x[1] * 1000 for x in pts]
    color = "#2878B5" if p["ellipse_semiaxes_m"]["a"] < .004 else "#D55E00"
    line = "-" if p["rotation_rad"] == 0 else "--"; marker = "o" if int(p["quadrature_order_per_smooth_arc"]) == 12 else "^"
    axg.plot(xs + [xs[0]], ys + [ys[0]], color=color, ls=line, lw=1, alpha=.65)
    axg.scatter(xs, ys, color=color, marker=marker, s=13, alpha=.65)
handles = [Line2D([0],[0],color="#2878B5",lw=2,label="small ellipse"), Line2D([0],[0],color="#D55E00",lw=2,label="large ellipse"), Line2D([0],[0],color="#444",ls="-",label="rotation 0.00"), Line2D([0],[0],color="#444",ls="--",label="rotation 0.37"), Line2D([0],[0],color="#444",marker="o",ls="None",label="12 / arc"), Line2D([0],[0],color="#444",marker="^",ls="None",label="16 / arc")]
axg.legend(handles=handles, fontsize=8, loc="upper right"); axg.set_aspect("equal", adjustable="datalim")
axg.set_xlabel("x (mm)"); axg.set_ylabel("y (mm)"); axg.set_title("Stored contour coordinates (224 points)"); axg.grid(alpha=.2)
x = list(range(8)); j12, j16, labels = [], [], []
for p in paths:
    j12.append(p["levels"]["volume_gl4x8x4"]["recomputed_J_J_per_m2"] * 1000)
    j16.append(p["levels"]["volume_gl4x16x4"]["recomputed_J_J_per_m2"] * 1000)
    size = "S" if p["ellipse_semiaxes_m"]["a"] < .004 else "L"
    labels.append(f"{size} r{p['rotation_rad']:.2f} q{int(p['quadrature_order_per_smooth_arc'])}")
axj.vlines(x, j12, j16, color="#777", lw=1, alpha=.65)
axj.scatter(x, j12, marker="o", color="#4477AA", label="volume_gl4x8x4")
axj.scatter(x, j16, marker="s", color="#CC6677", label="volume_gl4x16x4")
axj.axhline(JEXACT * 1000, color="#222", ls="--", label="analytic reference 13.0")
axj.set_xticks(x); axj.set_xticklabels(labels, rotation=35, ha="right", fontsize=8)
axj.set_ylabel("J (mJ/m²)"); axj.set_title("Saved-array J by the 8 contour paths"); axj.grid(axis="y", alpha=.25); axj.legend(fontsize=8)
fig.suptitle("Independent math.fsum review of saved solved-field rows", fontsize=14); fig.tight_layout()
png_path = OUT / "solved_j_review.png"; fig.savefig(png_path, dpi=170, bbox_inches="tight"); plt.close(fig)
(OUT / "source_run_receipt.json").write_bytes(rc_raw)
review = ("# Saved-array J 独立复核\n\n"
    f"按stage现成8条contour的坐标、法向、权重和梯度，每档224点依给定式以math.fsum重算；两档共448样点。样点字段与保存数组逐项一致；54项IDs和系数binding均与primary逐项相同。\n\n"
    f"对照解析值0.013 J/m²：每弧12点绝对误差{min(q12):.3f}–{max(q12):.3f}%，每弧16点{min(q16):.3f}–{max(q16):.3f}%；逐路径重算与保存J最大差{max_delta:.3g} J/m²。数值和8路径图见附件。\n\n"
    f"原coupling receipt保持UNKNOWN_WORKER_RESULT_CONTRACT，原件未改并逐字节副本归档。此前审阅器两次JSON解析失败，2.4737s未产数值；这是解析失败，不是物理失败。细档控制点|grad_z|max={zrows[1]['maximum_abs_grad_z_m_per_m_at_controls']:.12g}非零，因此这里只报告midsection functional，不据此称全域二维或实际G。\n")
need(len(review) <= 900, "REVIEW.md exceeds 900 characters")
review_path = OUT / "REVIEW.md"; review_path.write_text(review, encoding="utf-8")

finished = dt.datetime.now(dt.timezone.utc); elapsed = time.perf_counter() - t0
paths_in = {"source_run_receipt": RCP, "stage": stage_path, "primary": PRI}
after = {k: rec(p, p.read_bytes()) for k, p in paths_in.items()}
unchanged = {k: after[k]["sha256"] == inputs[k]["sha256"] for k in inputs}
need(all(unchanged.values()), "an input hash changed during review")
outputs = {}
for name, path in {"numeric_results": result_path, "review": review_path, "plot": png_path, "source_receipt_copy": OUT / "source_run_receipt.json", "review_script": pathlib.Path(__file__)}.items():
    outputs[name] = rec(path, path.read_bytes())
run = {"schema": "solved_j_independent_review_receipt.v1", "started_at_utc": started.isoformat(), "finished_at_utc": finished.isoformat(), "elapsed_seconds_excluding_python_imports": elapsed, "within_10s_processing_budget": elapsed <= 10, "input_bytes_total": sum(x["bytes"] for x in inputs.values()), "interpreter": str(pathlib.Path(sys.executable).resolve()), "input_hashes_before": inputs, "input_hashes_after": after, "inputs_unchanged": unchanged, "parent_gate_preserved": rc["result_status"], "review_character_count": len(review), "outputs": outputs}
save(OUT / "review_receipt.json", run)
print(json.dumps({"elapsed_seconds": elapsed, "points_per_level": total, "max_stored_J_delta": max_delta, "q12_error_pct": summary["q12_abs_error_to_analytic_pct_range"], "q16_error_pct": summary["q16_abs_error_to_analytic_pct_range"], "review_chars": len(review), "inputs_unchanged": unchanged}, allow_nan=False))