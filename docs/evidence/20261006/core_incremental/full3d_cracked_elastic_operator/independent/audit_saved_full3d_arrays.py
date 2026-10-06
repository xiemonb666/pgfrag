"""Independent algebraic readback from saved Full3D JSON only; no solver/provider imports."""
import hashlib, json, math
from collections import Counter; from pathlib import Path
import numpy as np; from PIL import Image, ImageDraw
ROOT = Path(r"C:\Users\xiemo\Desktop\ShardForge")
INPUT = ROOT / "tmp/review_20261006/full3d_cracked_elastic_operator_20261006/run01/result.json"
PINS = ROOT / "tmp/review_20261006/full3d_cracked_elastic_operator_20261006/SOURCE_PINS.json"
OUT = ROOT / "pgfrag/tmp/review_20261006/full3d_and_local_operator_independent/full3d_saved_array"
INPUT_SHA = "3b31e7d40c1d8e4b5d829f0b093fee4edadc966bde4086a4b3f4ba48fdf89cd2"
ENTRY_SHA = "e0f0cbbf0e60f65ce573ae8eed163e56c9b40a58be85b251445fd988e251936a"
REL_TOL = 5.0e-11
RIGID_TOL = 1.0e-10
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
relerr = lambda a, b: float(np.linalg.norm(a-b)/max(float(np.linalg.norm(b)), np.finfo(float).tiny))
finite = lambda a: bool(np.isfinite(a).all())
def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    input_sha = sha(INPUT)
    if input_sha != INPUT_SHA:
        raise ValueError(f"saved input SHA changed: {input_sha}")
    d = json.loads(INPUT.read_text(encoding="utf-8"))
    pins = json.loads(PINS.read_text(encoding="utf-8"))
    entry_path = ROOT / pins["entrypoint"]["path"]
    if pins["entrypoint"]["sha256"].lower() != ENTRY_SHA or sha(entry_path) != ENTRY_SHA:
        raise ValueError("Full3D frozen entrypoint hash mismatch")
    expected_deps = {x["name"]: x["sha256"] for x in pins["dependencies"]}
    if d["source_sha256"] != expected_deps:
        raise ValueError("result dependency pins differ from SOURCE_PINS")
    if (d["status"], d["assembly_status"], d["mechanical_eligibility"], d["solve_allowed"]) != ("OK_DIAGNOSTIC_ONLY_UNQUALIFIED", "OK", "DIAGNOSTIC_PASS_UNQUALIFIED", False):
        raise ValueError("saved result status is not the pinned successful diagnostic status")
    replay = d["static_solve_replay"]
    op, reg, src, quad = replay["operator"], replay["registry"], replay["sources"], replay["quadrature"]
    basis, boundary, exact = replay["basis"], replay["boundary"], replay["exact_field"]
    K0 = np.asarray(op["K_N_per_m"], dtype=np.float64)
    f0 = np.asarray(op["boundary_load_N"], dtype=np.float64)
    grad = np.asarray(basis["gradients_per_m"], dtype=np.float64)
    vals = np.asarray(basis["values"], dtype=np.float64)
    xyz = np.asarray(quad["positions_xyz_m"], dtype=np.float64)
    vw = np.asarray(quad["volume_weights_m3"], dtype=np.float64)
    bn = np.asarray(quad["boundary_normals"], dtype=np.float64)
    aw = np.asarray(quad["boundary_area_weights_m2"], dtype=np.float64)
    tr = np.asarray(boundary["traction_Pa"], dtype=np.float64)
    src_xyz = np.asarray(src["positions_xyz_m"], dtype=np.float64)
    keys, ids = reg["ordered_keys"], reg["ordered_canonical_ids"]
    basis_ids = replay["basis_row_ids"]
    nsrc, nq, nb, ndof, nvol, nbound = 27, 384, 54, 162, 256, 128
    shapes = (K0.shape == (ndof,ndof) and f0.shape == (ndof,) and grad.shape == (nq,nb,3) and vals.shape == (nq,nb)
              and xyz.shape == (nq,3) and vw.shape == (nvol,) and bn.shape == (nbound,3) and aw.shape == (nbound,)
              and tr.shape == (nbound,3) and src_xyz.shape == (nsrc,3) and len(keys) == len(ids) == ndof and len(basis_ids) == nb)
    if not shapes or not all(finite(a) for a in (K0, f0, grad, vals, xyz, vw, bn, aw, tr, src_xyz)):
        raise ValueError("saved replay has a malformed or non-finite numeric array")
    if np.any(vw <= 0.0) or np.any(aw <= 0.0):
        raise ValueError("saved quadrature contains a non-positive weight")
    if [k["slot"] for k in keys] != list(range(ndof)) or [k["canonical_id"] for k in keys] != ids:
        raise ValueError("canonical registry slots and IDs are not aligned")
    counts = Counter((k["source_id"], k["basis_kind"]) for k in basis_ids)
    if len(counts) != nb or Counter(k["basis_kind"] for k in basis_ids) != {"base": 27, "front_f1_unshifted": 27}:
        raise ValueError("cached scalar basis row IDs are not the 27 BASE + 27 global F1 rows")
    basis_row = {(k["source_id"], k["basis_kind"]): i for i, k in enumerate(basis_ids)}; comp = {"ux":0,"uy":1,"uz":2}
    material = d["fixed_case"]["material_and_manufactured_field"]
    E, nu = float(material["young_modulus_Pa"]), float(material["poisson_ratio"])
    mu = E / (2.0 * (1.0 + nu))
    lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    D = np.zeros((6, 6), dtype=np.float64)
    D[:3, :3] = lam
    D[np.arange(3), np.arange(3)] += 2.0 * mu
    D[3, 3] = D[4, 4] = D[5, 5] = mu
    K = np.zeros((ndof, ndof), dtype=np.float64)
    for q in range(nvol):
        B = np.zeros((6, ndof), dtype=np.float64)
        for slot, key in enumerate(keys):
            g = grad[q, basis_row[(key["source_id"], key["basis_kind"])]]
            c = comp[key["component"]]
            if c == 0:
                B[0, slot], B[3, slot], B[5, slot] = g[0], g[1], g[2]
            elif c == 1:
                B[1, slot], B[3, slot], B[4, slot] = g[1], g[0], g[2]
            else:
                B[2, slot], B[4, slot], B[5, slot] = g[2], g[1], g[0]
        K += vw[q] * (B.T @ D @ B)
    f = np.zeros(ndof, dtype=np.float64)
    for slot, key in enumerate(keys):
        brow, c = basis_row[(key["source_id"], key["basis_kind"])], comp[key["component"]]
        f[slot] = np.sum(vals[nvol:, brow] * tr[:, c] * aw, dtype=np.float64)
    k_rel, f_rel = relerr(K, K0), relerr(f, f0)
    eig = np.linalg.eigvalsh(0.5 * (K + K.T))
    eig0 = np.asarray(d["spectrum"]["symmetric_eigenvalues"], dtype=np.float64)
    if eig0.shape != (ndof,) or not finite(eig0):
        raise ValueError("saved eigenvalue vector is malformed or non-finite")
    eig_rel = relerr(eig, eig0)
    ztol = float(d["spectrum"]["zero_tolerance"])
    spectrum_counts = {"near_zero": int(np.count_nonzero(np.abs(eig)<=ztol)), "positive": int(np.count_nonzero(eig>ztol)), "significant_negative": int(np.count_nonzero(eig < -ztol))}
    center = src_xyz.mean(axis=0)
    R = np.zeros((ndof, 6), dtype=np.float64)
    axes = np.eye(3)
    source_row = {sid: i for i, sid in enumerate(reg["source_ids"])}
    for slot, key in enumerate(keys):
        if key["basis_kind"] != "base":
            continue
        c = comp[key["component"]]
        R[slot, c] = 1.0
        x = src_xyz[source_row[key["source_id"]]] - center
        for a in range(3):
            R[slot, 3 + a] = np.cross(axes[a], x)[c]
    R_saved = np.asarray(op["rigid_mode_coefficients"], dtype=np.float64)
    if R_saved.shape != (ndof, 6) or not finite(R_saved): raise ValueError("saved rigid coefficient matrix is malformed or non-finite")
    rigid_rel = relerr(R, R_saved)
    knorm = float(np.linalg.norm(K, ord="fro"))
    rigid_resid = [float(np.linalg.norm(K @ R[:, i]) / max(knorm * np.linalg.norm(R[:, i]), np.finfo(float).tiny)) for i in range(6)]
    coeff = np.asarray(op["manufactured_coefficients_m"], dtype=np.float64)
    if coeff.shape != (ndof,) or not finite(coeff): raise ValueError("preset manufactured coefficient vector is malformed or non-finite")
    preset_energy = float(0.5 * coeff @ K @ coeff)
    saved_energy = float(d["field_and_energy"]["coefficient_quadratic_energy_J"])
    energy_rel = abs(preset_energy - saved_energy) / max(abs(saved_energy), np.finfo(float).tiny)
    bxyz = xyz[nvol:]
    def wrench(a: int, b: int):
        df, dm = tr[a:b]*aw[a:b,None], np.cross(bxyz[a:b],tr[a:b])*aw[a:b,None]
        return np.sum(df,axis=0), np.sum(dm,axis=0), float(np.linalg.norm(df,axis=1).sum()), float(np.linalg.norm(dm,axis=1).sum())
    regions = {"side": wrench(0, 64), "top": wrench(64, 96), "bottom": wrench(96, 128), "total": wrench(0, 128)}
    net_center_m = regions["total"][1] - np.cross(center, regions["total"][0])
    work = R.T @ f
    net_center_wrench = np.concatenate((regions["total"][0], net_center_m))
    saved_wrench = d["integrated_boundary_wrench_and_rigid_virtual_work"]
    stored_work = np.asarray(saved_wrench["rigid_vectors_transpose_load_vector"], dtype=np.float64)
    if stored_work.shape != (6,) or not finite(stored_work): raise ValueError("saved rigid virtual work malformed")
    force_scale, moment_scale = regions["total"][2], regions["total"][3]
    work_rel = max(float(np.linalg.norm(work[:3]-stored_work[:3])/force_scale), float(np.linalg.norm(work[3:]-stored_work[3:])/moment_scale))
    work_balance_rel = max(float(np.linalg.norm(work[:3]-net_center_wrench[:3])/force_scale), float(np.linalg.norm(work[3:]-net_center_wrench[3:])/moment_scale))
    region_report = {}
    for name, (ff, mm, fscale, mscale) in regions.items():
        if name == "total":
            rf = np.asarray(saved_wrench["integrated_force_N"], dtype=np.float64)
            rm = np.asarray(saved_wrench["integrated_moment_about_global_origin_Nm"], dtype=np.float64)
        else:
            rw = saved_wrench[f"{name}_integrated_wrench"]
            rf, rm = np.asarray(rw["force_N"], dtype=np.float64), np.asarray(rw["moment_about_global_origin_Nm"], dtype=np.float64)
        region_report[name] = {"force_N": ff.tolist(), "moment_Nm": mm.tolist(), "saved_scaled_delta": float(max(np.linalg.norm(ff-rf)/fscale,np.linalg.norm(mm-rm)/mscale))}
    checks = {"K_reassembly_rel_le_5e-11": k_rel<=REL_TOL, "f_reassembly_rel_le_5e-11": f_rel<=REL_TOL,
        "eigenvalues_rel_le_5e-11": eig_rel<=REL_TOL, "rigid_vectors_match_le_5e-11": rigid_rel<=REL_TOL,
        "six_rigid_residuals_le_1e-10": max(rigid_resid)<=RIGID_TOL, "preset_energy_rel_le_5e-11": energy_rel<=REL_TOL,
        "rigid_work_and_wrenches_match_le_5e-11": max(work_rel,work_balance_rel,max(x["saved_scaled_delta"] for x in region_report.values()))<=REL_TOL,
        "spectrum_counts_6_156_0": spectrum_counts=={"near_zero":6,"positive":156,"significant_negative":0}}
    result = {
        "schema": "full3d_saved_array_independent_audit.v1",
        "status": "PASS_SAVED_ARRAY_ALGEBRA_ONLY" if all(checks.values()) else "FAIL_SAVED_ARRAY_ALGEBRA",
        "input_sha256": input_sha,
        "source_entrypoint_sha256": sha(entry_path),
        "source_pin_manifest_sha256": sha(PINS),
        "independent_method": "NumPy reconstruction from saved JSON arrays; no provider/assembler/reference/solver imports",
        "dimensions": {"source_count":nsrc,"scalar_basis_rows":nb,"semantic_dofs":ndof,"volume_points":nvol,"boundary_points":nbound},
        "array_finiteness": True,
        "positive_volume_and_area_weights": True,
        "semantic_basis_counts": {str(k): int(v) for k, v in sorted(Counter(x["basis_kind"] for x in basis_ids).items())},
        "K_relative_frobenius_delta": k_rel,
        "f_relative_l2_delta": f_rel,
        "K_relative_asymmetry": float(np.linalg.norm(K-K.T, ord="fro") / knorm),
        "eigenvalue_relative_l2_delta": eig_rel,
        "eigenvalue_min_N_per_m": float(eig[0]),
        "eigenvalue_max_N_per_m": float(eig[-1]),
        "recomputed_spectrum_counts": spectrum_counts,
        "zero_tolerance_N_per_m": ztol,
        "rigid_vector_relative_delta": rigid_rel,
        "rigid_residuals_recomputed": rigid_resid,
        "preset_manufactured_coefficient_quadratic_energy_J": preset_energy,
        "saved_preset_coefficient_quadratic_energy_J": saved_energy,
        "preset_energy_relative_delta": energy_rel,
        "coefficient_interpretation": "prescribed manufactured coefficients; not a solution obtained from K u=f",
        "recomputed_rigid_virtual_work": work.tolist(),
        "rigid_work_saved_relative_delta": work_rel,
        "rigid_work_vs_integrated_wrench_relative_delta": work_balance_rel,
        "regional_and_total_wrenches": region_report,
        "net_force_N": regions["total"][0].tolist(),
        "net_moment_about_global_origin_Nm": regions["total"][1].tolist(),
        "checks": checks,
        "execution_scope": "cached arrays only; no provider, original assembler, solver, field reference, CAD, or Git",
    }
    (OUT / "audit.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    width, height = 960, 560
    im = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(im)
    x0, x1, y0, y1 = 88, 30, 48, height-70
    xmax, ymax = width-x1, height-y1
    draw.text((x0, 14), "Full3D spectrum from independently rebuilt K", fill=(20, 20, 20))
    maxabs = max(float(np.max(np.abs(eig))), np.finfo(float).tiny)
    yv = np.clip(np.log10(np.maximum(np.abs(eig)/maxabs, 1e-16)), -16.0, 0.0)
    pts = [(int(x0+i/(len(eig)-1)*(xmax-x0)), int(y0+(0.0-float(v))/16.0*(ymax-y0))) for i,v in enumerate(yv)]
    for ytick in (0,-4,-8,-12,-16):
        yy=int(y0+(0-ytick)/16.0*(ymax-y0)); draw.line((x0,yy,xmax,yy),fill=(220,220,220),width=1); draw.text((42,yy-6),str(ytick),fill=(70,70,70))
    draw.line((x0,y0,x0,ymax),fill=(20,20,20),width=2); draw.line((x0,ymax,xmax,ymax),fill=(20,20,20),width=2)
    draw.line(pts,fill=(35,89,145),width=2)
    for i,v in enumerate(eig):
        if abs(float(v)) <= ztol:
            x,y=pts[i]; draw.ellipse((x-3,y-3,x+3,y+3),fill=(230,120,20))
    draw.text((x0,ymax+14),"eigenvalue rank (sorted)",fill=(30,30,30))
    im.save(OUT/"spectrum.png",format="PNG",optimize=True)
    print(json.dumps({"status":result["status"],"checks":checks,"Krel":k_rel,"frel":f_rel,"eigrel":eig_rel,"rigid_max":max(rigid_resid),"energy_rel":energy_rel},allow_nan=False))
    return 0 if all(checks.values()) else 1
if __name__ == "__main__":
    raise SystemExit(main())
